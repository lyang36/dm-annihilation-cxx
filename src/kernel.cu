#include <cstdio>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "cuda.h"
#include "math_functions.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "driver_functions.h"
#include "kernel.h"

namespace kernel_space{
	MAPTYPE * healpixX_GPU;
	MAPTYPE * healpixY_GPU;
	MAPTYPE * healpixZ_GPU;
	
	MAPTYPE * map_GPU;
	
	DMParticle * parts_GPU;

	//the normlized vector
	MAPTYPE * norm_GPU;
	//this is a zero vector, use to zerolize the GPU vec
	MAPTYPE * norm_CPU;
	
	int Npix_;
	int memParts_;
	
}

using namespace kernel_space;

//atomicAdd for double
__device__ double atomicAdd(double* address, double val)
{
	unsigned long long int* address_as_ull =
                                         (unsigned long long int*)address;
        unsigned long long int old = *address_as_ull, assumed;
            do {
                            assumed = old;
                            old = atomicCAS(address_as_ull, assumed,
                                                                    __double_as_longlong(val +
                                                                                                           __longlong_as_double(assumed)));
                                } while (assumed != old);
                return __longlong_as_double(old);
}

__device__ double SPHKenerl(MAPTYPE d2){
	return exp(-0.5 * d2 / 0.333);
}

__device__ MAPTYPE dotProd(MAPTYPE ax, MAPTYPE ay, MAPTYPE az, MAPTYPE bx, MAPTYPE by, MAPTYPE bz){
	return ax * ax + ay * ay + az * az;
}

//use atomicAdd, may affect perfermance
__global__ void calculateNorm(int Npix, MAPTYPE * healpixX, MAPTYPE * healpixY, MAPTYPE * healpixZ, 
		int numParts, DMParticle * parts, MAPTYPE * norm){
	int pix = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(pix >= Npix){
		return;
	}
	
	int i = 0;
	for(i = 0; i < numParts; i ++){
		if(parts[i].eps < 0){
			continue;
		}
		
		MAPTYPE prod = dotProd(healpixX[pix], healpixY[pix], healpixZ[pix]
		                      , parts[i].velx, parts[i].vely, parts[i].velz);
		
		//could add more constraints here
		if(prod < 0){
			continue;
		}
		
		MAPTYPE d2 = acos(prod) / parts[i].posy;
		if(d2 > 2){
			continue;
		}
		
		d2 = d2 * d2;
		MAPTYPE weight = SPHKenerl(d2);
		atomicAdd(&(norm[i]), weight);
	}
}

//no synchronized, very fast
__global__ void calculateMap(int Npix, MAPTYPE * healpixX, MAPTYPE * healpixY, MAPTYPE * healpixZ, 
		int numParts, DMParticle * parts, MAPTYPE * norm, MAPTYPE * map){
	int pix = blockIdx.x * blockDim.x + threadIdx.x;
	
	if(pix >= Npix){
		return;
	}
	
	int i = 0;
	for(i = 0; i < numParts; i ++){
		if(parts[i].eps < 0){
			continue;
		}
		
		MAPTYPE prod = dotProd(healpixX[pix], healpixY[pix], healpixZ[pix]
		                      , parts[i].velx, parts[i].vely, parts[i].velz);
		
		//could add more constraints here
		if(prod < 0){
			continue;
		}
		
		MAPTYPE d2 = acos(prod) / parts[i].posy;
		if(d2 > 2){
			continue;
		}
		
		d2 = d2 * d2;
		MAPTYPE weight = SPHKenerl(d2);
		map[pix] += weight / norm[i] * parts[i].mass;
	}
}

cudaError_t zeroLizeNorm(){
    //zerolize the norm
    cudaError_t cudaStatus = cudaMemcpy(norm_GPU, norm_CPU, memParts_ * sizeof(MAPTYPE), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX X!\n");
    }
    return cudaStatus;
}

cudaError_t initializeCUDA(MAPTYPE * healpixX, MAPTYPE * healpixY, MAPTYPE * healpixZ, int Npix, int memParts){
    cudaError_t cudaStatus;
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        return cudaStatus;
    }
    Npix_ = Npix;
    memParts_ = memParts;
    

    // Allocate GPU map.
    cudaStatus = cudaMalloc((void**)&map_GPU, Npix_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating map memory!");
        return cudaStatus;
    }
    
    // Allocate HEALPIX X.
    cudaStatus = cudaMalloc((void**)&healpixX_GPU, Npix_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating healpix x memory!");
        return cudaStatus;
    }
    // Allocate HEALPIX Y.
    cudaStatus = cudaMalloc((void**)&healpixY_GPU, Npix_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating healpix y memory!");
        return cudaStatus;
    }
    // Allocate HEALPIX Z.
    cudaStatus = cudaMalloc((void**)&healpixZ_GPU, Npix_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating healpix z memory!");
        return cudaStatus;
    }
    
    
    // Allocate GPU Particles.
    cudaStatus = cudaMalloc((void**)&parts_GPU, memParts_ * sizeof(DMParticle));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating Particles memory!");
        return cudaStatus;
    }
    
    // Allocate GPU Particle Norm.
    cudaStatus = cudaMalloc((void**)&norm_GPU, memParts_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating Particles norm memory!");
        return cudaStatus;
    }
    
    norm_CPU = (MAPTYPE *) calloc(memParts_, sizeof(MAPTYPE));
    zeroLizeNorm();
       
    //copy the HEALPIX X data to GPU
    cudaStatus = cudaMemcpy(healpixX_GPU, healpixX, Npix_ * sizeof(MAPTYPE), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX X!\n");
        return cudaStatus;
    }
    
    //copy the HEALPIX Y data to GPU
    cudaStatus = cudaMemcpy(healpixY_GPU, healpixY, Npix_ * sizeof(MAPTYPE), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX Y!\n");
        return cudaStatus;
    }
	
    //copy the HEALPIX Z data to GPU
    cudaStatus = cudaMemcpy(healpixZ_GPU, healpixZ, Npix_ * sizeof(MAPTYPE), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX Z!\n");
        return cudaStatus;
    }
    return cudaStatus;	
}


void cudaCleaingUp(){
	cudaFree(healpixX_GPU);
	cudaFree(healpixY_GPU);
	cudaFree(healpixZ_GPU);
	cudaFree(norm_GPU);
	cudaFree(map_GPU);
	cudaFree(parts_GPU);
	
	free(norm_CPU);
}

cudaError_t calulateMapWithCUDA(MAPTYPE * map, DMParticle * parts, int numParts){
	int blocksize = 512;
	int gridsize = Npix_ / blocksize + 1;
	
	zeroLizeNorm();
    //copy the Map data to GPU
    cudaError_t cudaStatus = cudaMemcpy(map_GPU, map, Npix_ * sizeof(MAPTYPE), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX X!\n");
        return cudaStatus;
    }
    
    //copy the particle data to GPU
    cudaStatus = cudaMemcpy(parts_GPU, parts, memParts_ * sizeof(DMParticle), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX X!\n");
        return cudaStatus;
    }
    
    calculateNorm<<<gridsize, blocksize>>>(Npix_, healpixX_GPU, healpixY_GPU, healpixZ_GPU, 
    		numParts, parts_GPU, norm_GPU);
    cudaStatus = cudaThreadSynchronize();
    if( cudaStatus != cudaSuccess){
        fprintf(stderr,"cudaThreadSynchronize error -- sync Norm: %s\n", cudaGetErrorString(cudaStatus));
        return cudaStatus;
    }
    
    calculateMap<<<gridsize, blocksize>>>(Npix_, healpixX_GPU, healpixY_GPU, healpixZ_GPU, 
    		numParts, parts_GPU, norm_GPU, map_GPU);
    cudaStatus = cudaThreadSynchronize();
    if( cudaStatus != cudaSuccess){
        fprintf(stderr,"cudaThreadSynchronize error -- sync map: %s\n", cudaGetErrorString(cudaStatus));
        return cudaStatus;
    }
    
    //copy map back
     cudaStatus = cudaMemcpy(map, map_GPU, Npix_ * sizeof(MAPTYPE), cudaMemcpyDeviceToHost);
     if (cudaStatus != cudaSuccess) {
         fprintf(stderr, "cudaMemcpy failed -- copying Map Back!\n");
         return cudaStatus;
     }

     return cudaStatus;
}
