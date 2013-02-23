#include <cstdio>
#include <cmath>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include <healpix_base.h>
#include <healpix_map.h>

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
    //healpix coordinates for super pixels in GPU
    MAPTYPE * heal_superx_GPU;
    MAPTYPE * heal_supery_GPU;
    MAPTYPE * heal_superz_GPU;
	
	MAPTYPE * map_GPU;
	
	DMParticle * parts_GPU;

	//the normlized vector
	MAPTYPE * norm_GPU;
	//this is a zero vector, use to zerolize the GPU vec
	MAPTYPE * norm_CPU;
	
	int Npix_;
    int Nside_;
	int memParts_;
    int NpixCoase_;
    MAPTYPE dOmegaSuper_;        //the dOmega of the super pixels
    
    
    //healpix pix coordinates
    MAPTYPE * healx;
    MAPTYPE * healy;
    MAPTYPE * healz;
    
    //healpix coordinates for super pixels
    MAPTYPE * heal_superx;
    MAPTYPE * heal_supery;
    MAPTYPE * heal_superz;
    
    //int * host_part_list;       //particle list for each super pixel
    int * dev_part_list;        //particle list for device for each super pixel
    int * host_part_num_list;   //particle number list for super pixels
    int * dev_part_num_list;    //particle number list for super pixels
	
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
       __double_as_longlong(val +__longlong_as_double(assumed)));
       } while (assumed != old);
    return __longlong_as_double(old);
}

__device__ double SPHKenerl(MAPTYPE d2){
	return exp(-0.5 * d2 / 0.333);
}

__device__ MAPTYPE dotProd(MAPTYPE ax, MAPTYPE ay, MAPTYPE az, MAPTYPE bx, MAPTYPE by, MAPTYPE bz){
    MAPTYPE  pd = ax * bx + ay * by + az * bz;
    if(pd > 1.0) pd = 1.0;
	return pd;
}

//use atomicAdd, may affect perfermance
__global__ void calculateNorm(int Npix, MAPTYPE * healpixX, MAPTYPE * healpixY, MAPTYPE * healpixZ,
                int numParts, DMParticle * parts, MAPTYPE * norm, int * numlist, int * partlist){
	//int pix = blockIdx.x * blockDim.x + threadIdx.x;
    int superpix = blockIdx.x;
    int pix = superpix * blockDim.x + threadIdx.x;
	
	if(pix >= Npix){
		return;
	}
    
	int start_pt = 0;
    if(superpix == 0){
        start_pt = 0;
    }else{
        start_pt = numlist[superpix - 1];
    }
    int end_pt = numlist[superpix];
    
	for(int i = start_pt; i < end_pt; i ++){
        int part_pt = partlist[i];
        
		if(parts[part_pt].eps >= 0){
            MAPTYPE prod = dotProd(healpixX[pix], healpixY[pix], healpixZ[pix]
		                      , parts[part_pt].velx, parts[part_pt].vely, parts[part_pt].velz);
		    if(prod >= parts[part_pt].posz){
		
		        MAPTYPE d2 = acos(prod) / parts[part_pt].posy;
		        //continue;
	            d2 = d2 * d2;
		        MAPTYPE weight = SPHKenerl(d2);
		        //testing...
                //norm[i] += weight;
                atomicAdd(&(norm[part_pt]), weight);
            }
        }
	}
}

//no synchronized, very fast
__global__ void calculateMap(int Npix, MAPTYPE * healpixX, MAPTYPE * healpixY, MAPTYPE * healpixZ,
                             int numParts, DMParticle * parts, MAPTYPE * norm, MAPTYPE * map,
                             int * numlist, int * partlist){
	//int pix = blockIdx.x * blockDim.x + threadIdx.x;
    int superpix = blockIdx.x;
    int pix = superpix * blockDim.x + threadIdx.x;
	
	if(pix >= Npix){
		return;
	}
    
	int start_pt = 0;
    if(superpix == 0){
        start_pt = 0;
    }else{
        start_pt = numlist[superpix - 1];
    }
    int end_pt = numlist[superpix];
    
	for(int i = start_pt; i < end_pt; i ++){
        int part_pt = partlist[i];
        
		if(parts[part_pt].eps >= 0){
            MAPTYPE prod = dotProd(healpixX[pix], healpixY[pix], healpixZ[pix]
                                   , parts[part_pt].velx, parts[part_pt].vely, parts[part_pt].velz);
		    if(prod >= parts[part_pt].posz){
                
		        MAPTYPE d2 = acos(prod) / parts[part_pt].posy;
		        //continue;
	            d2 = d2 * d2;
		        MAPTYPE weight = SPHKenerl(d2);
		        //testing...
                 map[pix] += weight * parts[part_pt].mass / norm[part_pt];
            }
        }
	}
    
    
    
}

//calculate how many particles are there in each coase pixel
//dOmega2--2*dOmega, cos2DOmega--cos(2*dOmega)
__global__ void calcuatePartsNum(int NpixCoase,  MAPTYPE * healpixSuperX, MAPTYPE * healpixSuperY, MAPTYPE * healpixSuperZ,
        int numParts, DMParticle * parts, int * numlist,
        MAPTYPE dOmega2, MAPTYPE cos2DOmega){
    int pix = blockIdx.x * blockDim.x + threadIdx.x;
    
    //clear them to be zero
    numlist[pix] = 0;
    __syncthreads();
    
    //search for particles 
    for(int i = 0; i < numParts; i++){
        if(parts[i].eps >= 0.0){
            MAPTYPE prod = dotProd(healpixSuperX[pix], healpixSuperY[pix], healpixSuperZ[pix]
                               , parts[i].velx, parts[i].vely, parts[i].velz);
            if(prod >= cos2DOmega){
                numlist[pix] ++;
            }else{
                MAPTYPE theta = acos(prod);
                if(theta <= dOmega2 + parts[i].posy){
                    numlist[pix] ++;
                }
            }
        }
    }
    //done
}


//calculate the particle list in each coase pixel
//dOmega2--2*dOmega, cos2DOmega--cos(2*dOmega)
__global__ void calcuatePartsList(int NpixCoase,  MAPTYPE * healpixSuperX, MAPTYPE * healpixSuperY, MAPTYPE * healpixSuperZ,
                                 int numParts, DMParticle * parts, int * numlist, int * partList,
                                 MAPTYPE dOmega2, MAPTYPE cos2DOmega){
    int pix = blockIdx.x * blockDim.x + threadIdx.x;
    
    int plist_ct = 0;
    int start_pt;
    if(pix == 0){
        start_pt = 0;
    }else{
        start_pt = numlist[pix - 1]
    }
    //search for particles
    for(int i = 0; i < numParts; i++){
        if(parts[i].eps >= 0.0){
            MAPTYPE prod = dotProd(healpixSuperX[pix], healpixSuperY[pix], healpixSuperZ[pix]
                                   , parts[i].velx, parts[i].vely, parts[i].velz);
            if(prod >= cos2DOmega){
                partList[start_pt + plist_ct] = i;
                plist_ct ++;
            }else{
                MAPTYPE theta = acos(prod);
                if(theta <= dOmega2 + parts[i].posy){
                    partList[start_pt + plist_ct] = i;
                    plist_ct ++;
                }
            }
        }
    }
    //done
}


cudaError_t zeroLizeNorm(){
    //zerolize the norm
    cudaError_t cudaStatus = cudaMemcpy(norm_GPU, norm_CPU, memParts_ * sizeof(MAPTYPE), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX X!\n");
    }
    return cudaStatus;
}

cudaError_t initializeHealpix(){

    
    Npix_ = 12 * Nside_ * Nside_;
    int NsideSuper = Nside_ / BLOCK_N_DIVIDER;
    
    //the actuall base
    Healpix_Base base(Nside_, NEST, SET_NSIDE);
    
    //the coased healpix base
    Healpix_Base super_base(NsideSuper, NEST, SET_NSIDE);
    NpixCoase_ = 12 * NsideSuper * NsideSuper;
    dOmegaSuper_ = 4 * PI / NpixCoase_;
    
    //healpix pix coordinates
    healx = new MAPTYPE[Npix_];
    healy = new MAPTYPE[Npix_];
    healz = new MAPTYPE[Npix_];
    
    //healpix coordinates for super pixels
    heal_superx = new MAPTYPE[NpixCoase];
    heal_supery = new MAPTYPE[NpixCoase];
    heal_superz = new MAPTYPE[NpixCoase];
    
    
    
    
    for(int i = 0; i < Npix_; i++){
    	vec3 this_vec = base.pix2vec(i);
    	healx[i] = this_vec.x;
    	healy[i] = this_vec.y;
    	healz[i] = this_vec.z;
    }
    
    for(int i = 0; i < NpixCoase_; i++){
    	vec3 this_vec = super_base.pix2vec(i);
    	heal_superx[i] = this_vec.x;
    	heal_supery[i] = this_vec.y;
    	heal_superz[i] = this_vec.z;
    }
    
    
    cudaError_t cudaStatus;
    
    host_part_num_list = new int[NpixCoase];
    cudaStatus = cudaMalloc((void**)&dev_part_num_list, NpixCoase * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating dev_part_num_list memory!\n");
        return cudaStatus;
    }
    
    
    // Allocate GPU map.
    cudaStatus = cudaMalloc((void**)&map_GPU, Npix_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating map memory!\n");
        return cudaStatus;
    }
    
    // Allocate HEALPIX X.
    cudaStatus = cudaMalloc((void**)&healpixX_GPU, Npix_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating healpix x memory!\n");
        return cudaStatus;
    }
    // Allocate HEALPIX Y.
    cudaStatus = cudaMalloc((void**)&healpixY_GPU, Npix_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating healpix y memory!\n");
        return cudaStatus;
    }
    // Allocate HEALPIX Z.
    cudaStatus = cudaMalloc((void**)&healpixZ_GPU, Npix_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating healpix z memory!\n");
        return cudaStatus;
    }
    
    // Allocate HEALPIX SUPER X.
    cudaStatus = cudaMalloc((void**)&heal_superx_GPU, NpixCoase_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating healpix SUPER x memory!\n");
        return cudaStatus;
    }
    // Allocate HEALPIX SUPER Y.
    cudaStatus = cudaMalloc((void**)&heal_supery_GPU, NpixCoase_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating healpix SUPER y memory!\n");
        return cudaStatus;
    }
    // Allocate HEALPIX SUPER Z.
    cudaStatus = cudaMalloc((void**)&heal_superz_GPU, NpixCoase_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating healpix SUPER z memory!\n");
        return cudaStatus;
    }
    
    //copy the HEALPIX X data to GPU
    cudaStatus = cudaMemcpy(healpixX_GPU, healx, Npix_ * sizeof(MAPTYPE), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX X!\n");
        return cudaStatus;
    }
    
    //copy the HEALPIX Y data to GPU
    cudaStatus = cudaMemcpy(healpixY_GPU, healy, Npix_ * sizeof(MAPTYPE), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX Y!\n");
        return cudaStatus;
    }
	
    //copy the HEALPIX Z data to GPU
    cudaStatus = cudaMemcpy(healpixZ_GPU, healz, Npix_ * sizeof(MAPTYPE), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX Z!\n");
        return cudaStatus;
    }

    
    //copy the HEALPIX Super X data to GPU
    cudaStatus = cudaMemcpy(heal_superx_GPU, heal_superx, NpixCoase_ * sizeof(MAPTYPE), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX Super X!\n");
        return cudaStatus;
    }
    
    //copy the HEALPIX Super Y data to GPU
    cudaStatus = cudaMemcpy(heal_supery_GPU, heal_supery, NpixCoase_ * sizeof(MAPTYPE), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX Super Y!\n");
        return cudaStatus;
    }
	
    //copy the HEALPIX Super Z data to GPU
    cudaStatus = cudaMemcpy(heal_superz_GPU, heal_superz, NpixCoase_ * sizeof(MAPTYPE), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX Super Z!\n");
        return cudaStatus;
    }
}

cudaError_t initializeCUDA(int Nside, int memParts){
    cudaError_t cudaStatus;
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?\n");
        return cudaStatus;
    }
    memParts_ = memParts;
    Nside_ = Nside;
    
    
    //initialze HEALPIX
    cudaStatus = initializeHealpix();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "Initialize HEALPIX FAILED!\n");
        return cudaStatus;
    }

    
    // Allocate GPU Particles.
    cudaStatus = cudaMalloc((void**)&parts_GPU, memParts_ * sizeof(DMParticle));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating Particles memory!\n");
        return cudaStatus;
    }
    
    // Allocate GPU Particle Norm.
    cudaStatus = cudaMalloc((void**)&norm_GPU, memParts_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating Particles norm memory!\n");
        return cudaStatus;
    }
    
    norm_CPU = (MAPTYPE *) calloc(memParts_, sizeof(MAPTYPE));
    zeroLizeNorm();
       
    return cudaStatus;	
}

void cudaCleaingUp(){
	cudaFree(healpixX_GPU);
	cudaFree(healpixY_GPU);
	cudaFree(healpixZ_GPU);
	cudaFree(norm_GPU);
	cudaFree(map_GPU);
	cudaFree(parts_GPU);
    
    cudaFree(heal_superx_GPU);
	cudaFree(heal_supery_GPU);
	cudaFree(heal_superz_GPU);
    
    cudaFree(dev_part_list);
    cudaFree(dev_part_num_list);
	
	free(norm_CPU);
}

//step 1: calculate how many particles are there in each super pixel, and allocating memory
cudaError_t calulatePartsNumListH(int numParts){
    cudaError_t cudaStatus;
    
    int blocksize = 512;
	int gridsize = NpixCoase_ / blocksize + 1;
    calcuatePartsNum<<<gridsize, blocksize>>>(NpixCoase_, heal_superx_GPU, heal_supery_GPU, heal_superz_GPU,
                                              numParts, parts_GPU, dev_part_num_list, dOmegaSuper_ * 2, cos(dOmegaSuper_ * 2));
    
    cudaStatus = cudaThreadSynchronize();
    if( cudaStatus != cudaSuccess){
        fprintf(stderr,"Sync particle num list calculating error: %s\n", cudaGetErrorString(cudaStatus));
        return cudaStatus;
    }
    
    //copy the num_list back
    cudaStatus = cudaMemcpy(host_part_num_list, dev_part_num_list, NpixCoase_ * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying num_list Back!\n");
        return cudaStatus;
    }
    
    //calculating the total memory needed for particle list
    for(int i = 1; i < NpixCoase_; i ++){
        host_part_num_list[i] += host_part_num_list[i - 1];
    }
    
    int memory_need_for_part_list = host_part_num_list[NpixCoase_ - 1];
    cudaFree(dev_part_list);
    
    //allocating memory for particle list
    cudaStatus = cudaMalloc((void**)&dev_part_list, memory_need_for_part_list * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating dev_part_list memory!\n");
        return cudaStatus;
    }
    
    //copy particle num list to GPU
    cudaStatus = cudaMemcpy(dev_part_num_list, host_part_num_list, NpixCoase_ * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying num_list Back!\n");
        return cudaStatus;
    }
}

//step 2: calculate the particle list in each
cudaError_t calulatePartsListH(int numParts){
    int blocksize = 512;
	int gridsize = NpixCoase_ / blocksize + 1;
    calcuatePartsList<<<gridsize, blocksize>>>(NpixCoase_, heal_superx_GPU, heal_supery_GPU, heal_superz_GPU,
                                              numParts, parts_GPU, dev_part_num_list, dev_part_list,
                                               dOmegaSuper_ * 2, cos(dOmegaSuper_ * 2));
    cudaStatus = cudaThreadSynchronize();
    if( cudaStatus != cudaSuccess){
        fprintf(stderr,"Sync particle list calculating error: %s\n", cudaGetErrorString(cudaStatus));
        return cudaStatus;
    }
    
}

cudaError_t calulateMapH(int numParts){
    int blocksize = BLOCK_N_DIVIDER * BLOCK_N_DIVIDER;
	int gridsize = Npix_ / blocksize + 1;
    

    //printf("Start Norm kernel...\n");
    calculateNorm<<<gridsize, blocksize>>>(Npix_, healpixX_GPU, healpixY_GPU, healpixZ_GPU,
                                           numParts, parts_GPU, norm_GPU, dev_part_num_list, dev_part_list);
    cudaStatus = cudaThreadSynchronize();
    if( cudaStatus != cudaSuccess){
        fprintf(stderr,"cudaThreadSynchronize error -- sync Norm: %s\n", cudaGetErrorString(cudaStatus));
        return cudaStatus;
    }
    
    //printf("Start Map kernel...\n");
    calculateMap<<<gridsize, blocksize>>>(Npix_, healpixX_GPU, healpixY_GPU, healpixZ_GPU,
                                          numParts, parts_GPU, norm_GPU, map_GPU,
                                          dev_part_num_list, dev_part_list);
    cudaStatus = cudaThreadSynchronize();
    if( cudaStatus != cudaSuccess){
        fprintf(stderr,"cudaThreadSynchronize error -- sync map: %s\n", cudaGetErrorString(cudaStatus));
        return cudaStatus;
    }
    
}

cudaError_t calulateMapWithCUDA(MAPTYPE * map, DMParticle * parts, int numParts){
    //copy the particle data to GPU
    cudaStatus = cudaMemcpy(parts_GPU, parts, memParts_ * sizeof(DMParticle), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying HEALPIX X!\n");
        return cudaStatus;
    }
    
    zeroLizeNorm();
    //copy the Map data to GPU
    cudaError_t cudaStatus = cudaMemcpy(map_GPU, map, Npix_ * sizeof(MAPTYPE), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying MAP!\n");
        return cudaStatus;
    }
    //step1:
    calulatePartsNumListH(numParts);
    //step2:
    calulatePartsListH(numParts);
    //step3:
    calulateMapH(numParts);
    
    //copy map back
     cudaStatus = cudaMemcpy(map, map_GPU, Npix_ * sizeof(MAPTYPE), cudaMemcpyDeviceToHost);
     if (cudaStatus != cudaSuccess) {
         fprintf(stderr, "cudaMemcpy failed -- copying Map Back!\n");
         return cudaStatus;
     }

     return cudaStatus;
}
