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


MAPTYPE * d_map;
renderpart * d_parts;
healpix_par params;
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
    params.setup(nside);
    
    int npix = 12 * nside * nside;
    cudaStatus = cudaMalloc((void**)&d_map, npix * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating HEALPix map memory!\n");
        return cudaStatus;
    }
    
    //clear the memory
    cudaMemset(d_map, 0, npix * sizeof(MAPTYPE));
    
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
    
    calcfluxGPU<<<blocks, NUM_THREADS_PER_BLOCK>>>(params,
                                                   d_map,
                                                   num_of_parts,
                                                   d_parts);
    

    return cudaStatus;
}

cudaError_t getCUDAMap(MAPTYPE * map){

	cudaError_t cudaStatus = cudaThreadSynchronize();
    if (cudaStatus != cudaSuccess) {
        printf("cudaThreadSynchronize error: %s\n", cudaGetErrorString(cudaStatus));
		return cudaStatus;
    }

    int npix = 12 * nside_ * nside_;
    cudaStatus = cudaMemcpy(map, d_map, npix * sizeof(MAPTYPE), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying map data to host!\n");

    }
    return cudaStatus;
}

void cudaCleaingUp(){
    cudaFree(d_map);
	cudaFree(d_parts);
}

void renderpart::setup(int nside){
    healpix_par params_;
    params_.setup(nside);
    
    //the ring number
    int iring;
    int icol;
   
    if( phi >= TWOPI)  phi = phi - TWOPI;
    if (phi < 0.)     phi = phi + TWOPI;    
    


    int ipix, dr;
    MAPTYPE rlat1 = theta - 2.0 * angular_radius;
    MAPTYPE zmax = cos(rlat1);
    rmin = ring_above (nside, zmax) - 1;
    MAPTYPE rlat2 = theta + 2.0 * angular_radius;
    MAPTYPE zmin = cos(rlat2);
    rmax = ring_above (nside, zmin) + 1;
    angle2pix(params_, z, phi, iring, icol, ipix);
    
    dc = (int)(4.0 * angular_radius / params.theta_per_pix +0.5) + 2;
    dr = rmax - rmin;
    
    numPix = (2 * dc + 1) * dr;
    
    bool isNorthPolarIn = (rlat1 <= 0.0);
    bool isSouthPolarIn = (rlat2 >= M_PI);
    
    npixNorthPole = 0;
    if(isNorthPolarIn)
        npixNorthPole = 2 *
            rmin * (rmin - 1);
    
    npixSouthPole = 0;
    if(isSouthPolarIn)
        npixSouthPole = params_.npix -
            2 * (params_.nl4 - rmax) *
            (params_.nl4 - rmax - 1);
    
    
}


__device__ MAPTYPE SPHKenerl(MAPTYPE d2){
    //test
	return exp(-0.5 * d2 / 0.333);
    //return exp(- 10 * d2 );

}

//(x1, y1, z1), (x2, y2, z2) must be normalized
__device__ MAPTYPE flux(healpix_par &par, MAPTYPE x1, MAPTYPE y1, MAPTYPE z1, MAPTYPE x2, MAPTYPE y2, MAPTYPE z2, MAPTYPE r_angle){
    MAPTYPE prod = x1 * x2 + y1 * y2 + z1 * z2;
    //test
    //return 1.0;

    if(prod > 1) prod = 1;
    if(prod < -1) prod = -1;
    MAPTYPE atheta = acos(prod);
    if(atheta > (par.theta_per_pix / 2.0 + r_angle * 2)){
        return 0;
    }

    //test
    //return 1.0;

    MAPTYPE d2 =  atheta / r_angle;
    d2 *= d2;
    return SPHKenerl(d2);
}


// 27 flops
__host__ __device__ void angle2pix(healpix_par &par,
                          MAPTYPE z,
                          MAPTYPE phi,
                          int & iring,
                          int & icol,
                          int & ipix){
    
    int jp, jm;
    MAPTYPE  tt, tp, tmp;
    int ir, ip, kshift;
    
    MAPTYPE za = fabs(z);
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
__host__ __device__ void pix2vec(healpix_par &par, int r, int c, MAPTYPE &x, MAPTYPE &y, MAPTYPE &z, MAPTYPE &ct, MAPTYPE &phi){
    MAPTYPE sintheta;
    if(r <= par.nside){
        ct  = 1.0 - 4.0 * r * r / (MAPTYPE) par.npix;
        phi = (c + 0.5) / (2.0 * r) * M_PI;
    }else if (r < par.nl3){
        ct  = 2.0 / 3.0 * (2.0 * par.nside - r) / (MAPTYPE) par.nside;
        phi = (c + 0.5 * (1 - (r + par.nside) % 2)) / (2.0 * par.nside) * M_PI;
    }else{
        MAPTYPE cr = par.nl4 - r;
        ct = -1.0 + 4.0 * cr * cr / (MAPTYPE)par.npix;
        phi = (c + 0.5) / (2 * cr) * M_PI;
    }
    sintheta = sqrt(1 - ct * ct);
    x = sintheta * cos(phi);
    y = sintheta * sin(phi);
    z = ct;
}



//4 flops
__host__ __device__ int cr2pix(healpix_par &par, int col, int ring){
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
__host__ __device__ int pix2icol(healpix_par &par, int ring, int pix){
    return pix - cr2pix(par, 0, ring);
}

//3 flops
__host__ __device__ int pix2ring(healpix_par &par, int ipix){
    if(ipix < par.ncap + par.nl4){
        return (int)((sqrt(2.0 * ipix + 1.0) + 1.0) / 2.0);
    }else if (ipix < 10 * par.nsidesq - par.nl4){
        return ((ipix - 2 * par.nsidesq - par.nl2)/par.nl4) + par.nside + 1;
    }else{
        //floor((sqrt(2*(12 * ns * ns - pix - 1) + 1) + 1) / 2)
        return par.nl4 - (int) floor((sqrt(2.0*(par.npix - 1.0 - ipix) + 1.0) + 1.0) / 2.0);
    }
}

//4 flops
//get the ring num of certain z
__host__ __device__ int ring_above (long nside_, MAPTYPE z){
    MAPTYPE az=abs(z);
    if (az>TWOTHIRD) // polar caps
    {
        int iring = (int)(nside_*sqrt(3*(1-az)));
        return (z>0) ? iring : (4*nside_-iring-1);
    }
    else // ----- equatorial region ---------
        return (int)(nside_*(2-1.5f*z));
}

__device__ void getPixIdThread(
                healpix_par &params, renderpart &particle,
                int k, int dc, int rmin, MAPTYPE phi,
                int  npixNorthPole, int npixSouthPole,
                int &pr, int &pc, int &p){
        int c0 = 0;
        p = -1;
        if (k < (npixNorthPole)){
                
                p=k;
                pr = pix2ring(params, p);
                pc = pix2icol(params, pr, p);
                
        }
        else if (k < (npixNorthPole + npixSouthPole)){
                
                p = params.npix - (k - npixNorthPole) - 1;
                pr = pix2ring(params, p);
                pc = pix2icol(params, pr, p);
                
        }else{
                int np = k - npixNorthPole - npixSouthPole;
                pr = np / (2 * dc +1)+rmin;
                if(pr < 1 || pr > params.nl4){
                        p = -1;
                }else{
                    int  npixatthisring = params.nl4;
                        //c0 = (int)((phi - 2.0 * particle.angular_radius)
                        //             / params.theta_per_pix - 1);
                    if(pr <= params.nside){
                        //c0 = (int)floor(2 * (phi - (2.0 * particle.angular_radius)) * pr / M_PI
                        //                - 0.5)
                        //                + 1;
                        npixatthisring = 4 * pr;
                    }
                    else if(pr < params.nl3){
                        //c0 = (int)floor((phi - 2.0 * particle.angular_radius)
                        //           / params.theta_per_pix - 0.5) + 1;
                        npixatthisring = params.nl4;
                    }else{
                        //c0 = (int)floor(2 * (phi - 2.0 * particle.angular_radius) *
                        //           (params.nl4 - pr) / M_PI 
                        //           - 0.5) + 1;
                        npixatthisring = 4 * (params.nl4 - pr);
                    }
				
					c0 = floor(npixatthisring/TWOPI*(phi) - 0.5) - dc;
                    //int ip_hi = floor(npixatthisring/TWOPI*(phi0+dphi) - 0.5);


                    pc = np % (2 * dc +1)+c0;
                    if(pc < 0) pc += npixatthisring;
                    if(pc > npixatthisring) pc = pc % npixatthisring;
                    p = cr2pix(params, pc, pr);
                    
              }
       }
       if((p < 0) || (p > params.npix)){
              p = -1;
       }
                
}


__device__ double atomicAdd(double* address, double val)
{
	unsigned long long int* address_as_ull =
                                         (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
       assumed = old;
       old = atomicCAS(address_as_ull, assumed,
       __double_as_longlong(val +__longlong_as_double(assumed)));
       } while (assumed != old);
    return __longlong_as_double(old);
}



//count flops
__global__ void calcfluxGPU(
                            healpix_par params,
                            MAPTYPE * map,
                            int numOfParts,
                            renderpart * parts
                            ){

    __shared__ MAPTYPE pixval[NUM_THREADS_PER_BLOCK];
    __shared__ renderpart listOfParticles[NUM_THREADS_PER_BLOCK];
    renderpart particle;
    
    if(blockIdx.x >= numOfParts){
        return;
    }
    
    ////////////////////////Read Particles/////////////////////////
    //listOfParticles[blockIdx.x] = particle;
    if(threadIdx.x == 0){
        particle = parts[blockIdx.x];
        listOfParticles[0] = particle;
    }
    __syncthreads();
    
    //down-sweeping
    int halfThreadNum = 1;
    while(halfThreadNum < NUM_THREADS_PER_BLOCK){
        if(threadIdx.x < halfThreadNum){
            listOfParticles[halfThreadNum + threadIdx.x]
                = listOfParticles[threadIdx.x];
        }
        __syncthreads();
        halfThreadNum *= 2;
    }
    particle = listOfParticles[threadIdx.x];
    ////////////////////////////////////////////////////////////////
    
    int rmin;
    int dc;
    int numPix, npixNorthPole, npixSouthPole, totalPix;
    
    int numPixPerThread;
    int startPix = 0;
    
    MAPTYPE x, y, z;
    
    MAPTYPE  phi;
    MAPTYPE norm = 0;

    
    dc = particle.dc;
    x = particle.x;
    y = particle.y;
    z = particle.z;
    phi = particle.phi;
    numPix = particle.numPix;
    npixNorthPole = particle.npixNorthPole;
    npixSouthPole = particle.npixSouthPole;
    rmin = particle.rmin;
    
    MAPTYPE za = fabs(z);
    
    totalPix = numPix + npixNorthPole + npixSouthPole;
    
    
    numPixPerThread = totalPix / NUM_THREADS_PER_BLOCK
        + ((totalPix % NUM_THREADS_PER_BLOCK == 0) ? 0 : 1);
    
    startPix = 0;
    

    //int c0;
    int p = 0;
    int pr = 0;
    int pc = 0;
    MAPTYPE weight = 0;
    int k = 0; 
    norm = 0.0;
    
    for(int i = 0; i < numPixPerThread; i++){
        weight = 0;

        k = threadIdx.x + startPix;
        if( k < totalPix){
        
            //calculate the pixel id
                                
            getPixIdThread(params, particle,
                    k, dc, rmin, phi,
                    npixNorthPole, npixSouthPole,
                    pr, pc, p);
                
            if(p != -1){ 
            
                //calculate the value
                MAPTYPE x1, y1, z1, ct, phi1;
                pix2vec(params, pr, pc, x1, y1, z1, ct, phi1);
                //test
                weight = flux(params, x1, y1, z1,
                          x, y, z,
                          particle.angular_radius);
                //pixval[threadIdx.x] = weight;
            }
        }
        pixval[threadIdx.x] = weight;
        
        /////////////////////Calculating Norm////////////////////////
        //calculate the norm (reduce-sweeping algorithm)
        //up-sweeping
        __syncthreads();
        halfThreadNum = NUM_THREADS_PER_BLOCK / 2;
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
        norm += pixval[threadIdx.x];
        ////////////////////////////////////////////////////////////
        
        //calculated the result and record them to the global memory
        
        if(numPixPerThread == 1){
			if((norm > 0) && (p != -1)){
				//test
                MAPTYPE fofp = weight * particle.flux / norm;
                atomicAdd(map + p, fofp);
				//atomicAdd(map + p, weight / norm);
				//map[p] = particle.flux;
			}
        }

        
        //continue to next block
        startPix += NUM_THREADS_PER_BLOCK;
    }
    
    if(numPixPerThread <= 1){
        return;
    }
    
    startPix = 0;
    for(int i = 0; i < numPixPerThread; i++){
        k = threadIdx.x + startPix;
        if( k >= totalPix){
        }else{
            
            //calculate the pixel id
            getPixIdThread(params, particle,
                    k, dc, rmin, phi,
                    npixNorthPole, npixSouthPole,
                    pr, pc, p);
            
            if(p != -1 ){
                //calculate the value
                MAPTYPE x1, y1, z1, ct, phi1;
                pix2vec(params, pr, pc, x1, y1, z1, ct, phi1);
            
			    if(norm > 0){
                    //test 
                    weight = flux(params, x1, y1, z1,
                          x, y, z,
                          particle.angular_radius);

				    //test
            	    atomicAdd(map + p, weight * particle.flux / norm);
				    //atomicAdd(map + p, weight / norm);
			    }

            }
        }
        
        startPix += NUM_THREADS_PER_BLOCK;

    }
}
