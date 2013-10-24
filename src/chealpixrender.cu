#include <cstdio>
//#include <cmath>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include <healpix_base.h>
#include <healpix_map.h>

#include "cuda.h"
#include "math_functions.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "driver_functions.h"
#include "chealpixrender.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define PIOVER2  (0.5*M_PI)
#define PI M_PI
#define TWOPI (2.0*M_PI)
#define Z0 (2.0/3.0)
#define TWOTHIRD Z0
#define NS_MAX 8192


float * d_map;
renderpart * d_parts;
int nside_;

cudaError_t initializeCUDA(int nside, int numofparts){

    int count = 0;
	int i = 0;
    cudaError_t cudaStatus;

	cudaGetDeviceCount(&count);
	if(count == 0) {
		fprintf(stderr, "There is no device.\n");
		return cudaErrorNotReady;
	}

	for(i = 0; i < count; i++) {
		cudaDeviceProp prop;
		if(cudaGetDeviceProperties(&prop, i) == cudaSuccess) {
			if(prop.major >= 1) {
				break;
			}
		}
	}
	if(i == count) {
		fprintf(stderr, "There is no device supporting CUDA.\n");
		return cudaErrorNotReady;
	}
	cudaSetDevice(i);

	printf("CUDA initialized.\n");

    nside_ = nside;
    int npix = 12 * nside * nside;
    cudaStatus = cudaMalloc((void**)&d_map, npix * sizeof(float));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating HEALPix map memory!\n");
        return cudaStatus;
    }
    
    //clear the memory
    cudaMemset(d_map, 0, npix * sizeof(float));
    
    cudaStatus = cudaMalloc((void**)&d_parts, numofparts * sizeof(renderpart));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating Render Particles memory!\n");
        return cudaStatus;
    }
    return cudaStatus;
}


cudaError_t calculateMapByGPU(renderpart * parts, int num_of_parts){
    int blocks = num_of_parts;
    if(num_of_parts == 0){
        return cudaSuccess;
    }
    //cuda mem copy
    //copy particle data to GPU
    cudaError_t cudaStatus = cudaMemcpy(d_parts, parts, num_of_parts * sizeof(renderpart), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying particle data to device!\n");
        return cudaStatus;
    }
    
    calcfluxGPU<<<blocks, NUM_THREADS_PER_BLOCK>>>(nside_,
                                                   d_map,
                                                   num_of_parts,
                                                   d_parts);
    
    cudaStatus = cudaThreadSynchronize();
    if (cudaStatus != cudaSuccess) {
        printf("cudaThreadSynchronize error: %s\n", cudaGetErrorString(cudaStatus));
    }
    return cudaStatus;
}

cudaError_t getCUDAMap(float * map){
    int npix = 12 * nside_ * nside_;
    cudaError_t cudaStatus = cudaMemcpy(map, d_map, npix * sizeof(float), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying map data to host!\n");

    }
    return cudaStatus;
}

void cudaCleaingUp(){
    cudaFree(d_map);
	cudaFree(d_parts);
}

__device__ float SPHKenerl(float d2){
	return exp(-0.5 * d2 / 0.333);
}

//(x1, y1, z1), (x2, y2, z2) must be normalized
__device__ float flux(healpix_par &par, float x1, float y1, float z1, float x2, float y2, float z2, float r_angle){
    float prod = x1 * x2 + y1 * y2 + z1 * z2;
    if(prod > 1) prod = 1;
    float atheta = acos(prod);
    if(atheta > (par.theta_per_pix / 2.0 + r_angle)){
        return 0;
    }
    float d2 =  atheta / r_angle;
    d2 *= d2;
    return SPHKenerl(d2);
}


// 27 flops
__device__ void angle2pix(healpix_par &par,
                          float z,
                          float phi,
                          int & iring,
                          int & icol,
                          int & ipix){
    
    int jp, jm;
    float  tt, tp, tmp;
    int ir, ip, kshift;
    
    float za = fabs(z);
    if( phi >= TWOPI)  phi = phi - TWOPI;
    if (phi < 0.)     phi = phi + TWOPI;
    tt = phi / PIOVER2;//  ! in [0,4)
    
    if( za <= Z0 ) {
        jp = (int)floor(par.nside*(0.5 + tt - z*0.75)); /*index of ascending edge line*/
        jm = (int)floor(par.nside*(0.5 + tt + z*0.75)); /*index of descending edge line*/
        ir = par.nside + 1 + jp - jm;// ! in {1,2n+1} (ring number counted from z=2/3)
        iring = ir + par.nside; // ! in {1,3n+1} (ring number counted from z=1)
        
        kshift = 0;
        if (ir % 2==0.) kshift = 1;// ! kshift=1 if ir even, 0 otherwise
        ip = ( ( jp+jm - par.nside + kshift + 1 ) / 2 );// ! in {0,4n)
        if( ip>par.nl4 ) ip = ip - par.nl4;
        icol = ip;
        ipix = par.ncap + par.nl4*(ir-1) + ip ;
    }else {
        tp = tt - floor(tt);//      !MOD(tt,1.d0)
        tmp = sqrt( 3.*(1. - za) );
        
        jp = (int)floor( par.nside * tp * tmp );// ! increasing edge line index
        jm = (int)floor( par.nside * (1. - tp) * tmp );// ! decreasing edge line index
        
        ir = jp + jm + 1;//        ! ring number counted from the closest pole
        ip = (int)floor( tt * ir );// ! in {0,4*ir)
        
        if( ip>4*ir ) ip = ip - 4*ir;
        
        ipix = 2*ir*(ir-1) + ip;
        iring = ir;
        icol = ip;
        if( z<=0. ) {
            ipix = par.npix - 2*ir*(ir+1) + ip;
            iring = par.nl4 - ir;
        }
    }
}

//6 flops
__device__ void pix2vec(healpix_par &par, int r, int c, float &x, float &y, float &z, float &ct, float &phi){
    float sintheta;
    if(r <= par.nside){
        ct  = 1.0 - 4.0 * r * r / (float) par.npix;
        phi = (c + 0.5) / (2.0 * r) * M_PI;
    }else if (r < par.nl3){
        ct  = 2.0 / 3.0 * (2.0 * par.nside - r) / (float) par.nside;
        phi = (c + 0.5 * (1 - (r + par.nside) % 2)) / (2.0 * par.nside) * M_PI;
    }else{
        float cr = par.nl4 - r;
        ct = -1.0 + 4.0 * cr * cr / (float)par.npix;
        phi = (c + 0.5) / (2 * cr) * M_PI;
    }
    sintheta = sqrt(1 - ct * ct);
    x = sintheta * cos(phi);
    y = sintheta * sin(phi);
    z = ct;
}



//4 flops
__device__ int cr2pix(healpix_par &par, int col, int ring){
    if(ring <= par.nside){
        return 2 * ring * (ring - 1) + col;
    }else if(ring  < par.nl3){
        return par.nl4 * ring - par.nl2 * par.nside - par.nl2 + col;
    }else{
        int r = par.nl4 - ring;
        return par.npix - (2 * r * (r + 1)) + col;
    }
}

//4 flops
__device__ int pix2icol(healpix_par &par, int ring, int pix){
    return pix - cr2pix(par, 0, ring);
}

//3 flops
__device__ int pix2ring(healpix_par &par, int ipix){
    if(ipix < par.ncap + par.nl4){
        return ((sqrt(2.0 * ipix + 1) + 1) / 2);
    }else if (ipix < 10 * par.nsidesq - par.nl4){
        return ((ipix - 2 * par.nsidesq - par.nl2)/par.nl4) + par.nside + 1;
    }else{
        //floor((sqrt(2*(12 * ns * ns - pix - 1) + 1) + 1) / 2)
        return par.nl4 - floor((sqrt(2.0*(par.npix - 1- ipix) + 1) + 1) / 2);
    }
}

//4 flops
//get the ring num of certain z
__device__ int ring_above (long nside_, float z){
    float az=abs(z);
    if (az>TWOTHIRD) // polar caps
    {
        int iring = int(nside_*sqrt(3*(1-az)));
        return (z>0) ? iring : 4*nside_-iring-1;
    }
    else // ----- equatorial region ---------
        return int(nside_*(2-1.5*z));
}

//count flops
__global__ void calcfluxGPU(int nside,
                       float * map,
                       int numOfParts,
                       renderpart * parts){
    healpix_par pars;
    __shared__ float pixval[NUM_THREADS_PER_BLOCK];
    
    int icol, iring, dc, rmin, rmax, dr, ipix;
    bool isNorthPolarIn = false;
    bool isSouthPolarIn = false;
    int numPix, npixNorthPole = 0, npixSouthPole = 0, totalPix;
    int numPixPerThread;
    int startPix = 0;
    float z;
    float x, y;
    float  phi;
    float norm = 0;
    
    if(blockIdx.x >= numOfParts){
        return;
    }
    
    renderpart particle = parts[blockIdx.x];

    
    //6 flops
    pars.nl2 = 2*nside;
    pars.nl4 = 4*nside;
    pars.nl3 = 3*nside;
    pars.nside = nside;
    pars.nsidesq = nside * nside;
    pars.ncap  = pars.nl2*(nside-1);
            // ! number of pixels in the north polar cap
    pars.npix  = 12*nside*nside;
    pars.theta_per_pix = 2 * M_PI / pars.nl4;
    x = particle.x;
    y = particle.y;
    z = particle.z;
    //theta = particle.theta;
    phi = particle.phi;
    float za = fabs(z);
    
    //10 flops
    float rlat1 = particle.theta - particle.angular_radius;
    float zmax = cos(rlat1);
    rmin = ring_above (nside, zmax) + 1;
    float rlat2 = particle.theta + particle.angular_radius;
    float zmin = cos(rlat2);
    rmax = ring_above (nside, zmin) + 1;
    
    //27 flops
    //get the phi and col of theta and phi
    angle2pix(pars, z, particle.phi, iring, icol, ipix);
    
    dc = particle.angular_radius / pars.theta_per_pix + 1;
    dr = rmax - rmin;
    
    numPix = (2 * dc + 1) * dr;
    isNorthPolarIn = (rlat1 <= 0.0);
    isSouthPolarIn = (rlat2 >= M_PI);
    
    if(isNorthPolarIn)
            npixNorthPole = 2 * rmin * (rmin - 1);
    
    if(isSouthPolarIn)
            npixSouthPole = pars.npix - 2 * (pars.nl4 - rmax) * (pars.nl4 - rmax - 1);

    totalPix = numPix + npixNorthPole + npixSouthPole;
    
    numPixPerThread = totalPix / NUM_THREADS_PER_BLOCK
        + (totalPix % NUM_THREADS_PER_BLOCK == 0) ? 0 : 1;
    startPix = 0;
    
    bool isIgnored = false;

    int c0;
    int p = 0;
    int pr = 0;
    int pc = 0;
    float weight = 0;
    int k = 0;
    for(int i = 0; i < numPixPerThread; i++){
        pixval[threadIdx.x] = 0;
        __syncthreads();
        
        k = threadIdx.x + startPix;
        if( k > totalPix){
            break;
        }
        
        //calculate the pixel id
        if (k < npixNorthPole){
            p=k;
            pr = pix2ring(pars, p);
            pc = pix2icol(pars, pr, p);
        }
        else if (k < npixNorthPole + npixSouthPole){
            p = pars.npix - (k - npixNorthPole) - 1;
            pr = pix2ring(pars, p);
            pc = pix2icol(pars, pr, p);
        }else{
            int np = k - npixNorthPole - npixSouthPole;
            pr = np / (2 * dc +1)+rmin;
            if(pr < 1 || pr > pars.nl4){
                isIgnored = true;
            }else{
                if(pr <= pars.nside){
                    c0 = (int)(2 * (phi - particle.angular_radius) * pr / M_PI) - 1;
                }
                else if(pr < pars.nl3){
                    c0 = int(phi - particle.angular_radius)
                        / pars.theta_per_pix - 1;
                }else{
                    c0 = (int)(2 * (phi - particle.angular_radius) *
                        (pars.nl4 - pr) / M_PI) - 1;
                }
                pc = np % (2 * dc +1)+c0;
                if((pc < 0) ||
                   (pr <= pars.nside && pc >= 4 * pr) ||
                   (pr < pars.nl3 && pc >= pars.nl4) ||
                   (pr >= pars.nl3 && pc > 4 * (pars.nl4 - pr))){
                    isIgnored = true;
                }else{
                    p = cr2pix(pars, pc, pr);
                }
            }

        }
 
        //calculate the value
        if(!isIgnored){
            float x1, y1, z1, ct, phi;
            pix2vec(pars, pr, pc, x1, y1, z1, ct, phi);
            weight = pixval[i] = flux(pars, x1, y1, z1,
                            x, y, z,
                            particle.angular_radius);
        }
        
        
        //calculate the norm (reduce-sweeping algorithm)
        //up-sweeping
        int halfThreadNum = NUM_THREADS_PER_BLOCK / 2;
        while(halfThreadNum > 0){
            if(threadIdx.x < halfThreadNum){
                pixval[threadIdx.x] += pixval[halfThreadNum + threadIdx.x];
            }
            __syncthreads();
            halfThreadNum /= 2;
        }
        //down-sweeping
        halfThreadNum = 1;
        while(halfThreadNum < NUM_THREADS_PER_BLOCK){
            if(threadIdx.x < halfThreadNum){
                pixval[halfThreadNum + threadIdx.x] = pixval[threadIdx.x];
            }
            __syncthreads();
            halfThreadNum *= 2;
        }
        norm = pixval[threadIdx.x];
        
        //calculated the result and record them to the global memory
        if(!isIgnored){
            if(numPixPerThread == 1){
                //map[p]  = weight / norm;
                atomicAdd(map + p, weight * particle.flux / norm);
            }
        }
        
        //continue to next block
        startPix += NUM_THREADS_PER_BLOCK;
    }
    
    if(numPixPerThread == 1){
        return;
    }
    
    startPix = 0;
    for(int i = 0; i < numPixPerThread; i++){
        //calculate the pixel id
        if (k < npixNorthPole){
            p=k;
            pr = pix2ring(pars, p);
        }
        else if (k < npixNorthPole + npixSouthPole){
            p = pars.npix - (k - npixNorthPole) - 1;
            pr = pix2ring(pars, p);
            
        }else{
            int np = k - npixNorthPole - npixSouthPole;
            pr = np / (2 * dc +1)+rmin;
            if(pr < 1 || pr > pars.nl4){
                isIgnored = true;
            }else{
                if(pr <= pars.nside){
                    c0 = (int)(2 * (phi - particle.angular_radius) * pr / M_PI) - 1;
                }
                else if(pr < pars.nl3){
                    c0 = int(phi - particle.angular_radius)
                    / pars.theta_per_pix - 1;
                }else{
                    c0 = (int)(2 * (phi - particle.angular_radius) *
                               (pars.nl4 - pr) / M_PI) - 1;
                }
                pc = np % (2 * dc +1)+c0;
                if((pc < 0) ||
                   (pr <= pars.nside && pc >= 4 * pr) ||
                   (pr < pars.nl3 && pc >= pars.nl4) ||
                   (pr >= pars.nl3 && pc > 4 * (pars.nl4 - pr))){
                    isIgnored = true;
                }else{
                    p = cr2pix(pars, pc, pr);
                }
            }
            
        }
        
        //calculate the value
        if(!isIgnored){
            float x1, y1, z1, ct, phi;
            pix2vec(pars, pr, pc, x1, y1, z1, ct, phi);
            weight = flux(pars, x1, y1, z1,
                          x, y, z,
                          particle.angular_radius);
            atomicAdd(map + p, weight * particle.flux / norm);
        }

    }
}
