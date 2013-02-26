#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "cuda.h"
#include "math_functions.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "driver_functions.h"
#include "../kernel.h"

#include "halocore.h"
#include "cudahalocore.h"
#define RADIUS_RATIO (18.2 / 402.0)
#define SAT_RADIUS (1.0 / 400000.0 / RADIUS_RATIO)

namespace cuda_halo_core{
    static const string halo_core_file_name_ = "./halocore/VL2_Halos.data";
	HaloData * halos_;
    int numHalos_;
	int memNumParts_;;

	HaloData * dev_halos;
	DMParticle * dev_dmparts;
	MAPTYPE * dev_result;
};

using namespace cuda_halo_core;
using namespace std;

__device__ MAPTYPE coreFunc(float x, float y, float z, 
                float hxc, float hyc, float hzc, float radius){
    MAPTYPE xx = (x - hxc);
    MAPTYPE yy = (y - hyc);
    MAPTYPE zz = (z - hzc);
    MAPTYPE r = sqrt(xx*xx + yy*yy + zz*zz);
    MAPTYPE ratio = 1.0;
    if( r < radius){
        if( r > SAT_RADIUS * radius)
            ratio = pow(r / radius, 0.6);
            //ratio = 0.0;    //remove the center
        else
            ratio = pow(r / SAT_RADIUS / radius, 0.6);
    }
    return ratio;
}



__global__ void calculateCoreCorrectionGPU(int numParts, int numHalos,
	HaloData * halos, DMParticle * parts, MAPTYPE * results){
	int id = blockIdx.x * blockDim.x + threadIdx.x;
	if(id >= numParts){
		return;	
	}
	MAPTYPE correction = 1.0;
	MAPTYPE xp = parts[id].posx;
	MAPTYPE yp = parts[id].posy;
	MAPTYPE zp = parts[id].posz;
	for(int i = 0; i < numHalos; i++){
		correction *= coreFunc(xp, yp, zp, halos[i].xc, 
                        halos[i].yc,
                        halos[i].zc, 
                        halos[i].radius * RADIUS_RATIO / 40000.0);	
	}
	results[id] = correction;
}


cudaError_t applyingCore(DMParticle * parts, int numParts, MAPTYPE * result){
	int blocksize = 512;
	int gridsize = numParts / blocksize + 1;

	cudaError_t cudaStatus = cudaMemcpy(dev_dmparts, parts, numParts * sizeof(DMParticle), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying particles!\n");
        return cudaStatus;
    }

	calculateCoreCorrectionGPU<<<gridsize, blocksize>>>(numParts, numHalos_,
		dev_halos, dev_dmparts, dev_result);

	cudaStatus = cudaThreadSynchronize();
    if( cudaStatus != cudaSuccess){
        fprintf(stderr,"Sync core correction calculating error: %s\n", cudaGetErrorString(cudaStatus));
        return cudaStatus;
    }

	 //copy result back
    cudaStatus = cudaMemcpy(result, dev_result, numParts * sizeof(MAPTYPE), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
         fprintf(stderr, "cudaMemcpy failed -- copying core_result Back!\n");
        
    }	
	return cudaStatus;
}

cudaError_t initializeCore(DMParticle * dmparts, int numParts){
	memNumParts_ = numParts;
	
	int number_of_lines = 0;
    std::string line;
    std::ifstream counterfile(halo_core_file_name_.c_str());

    while (std::getline(counterfile, line))
        ++number_of_lines;
    counterfile.close();
    if(number_of_lines == 0){
        printf("Are you input a empty halo file?\n");
        exit(0);
    }
    numHalos_ = number_of_lines;
    halos_ = new HaloData[number_of_lines];
    
    std::ifstream haloFile(halo_core_file_name_.c_str());
    int i = 0;
    while(haloFile.good()){
        haloFile >> halos_[i].xc;
        haloFile >> halos_[i].yc;
        haloFile >> halos_[i].zc;
        haloFile >> halos_[i].radius;
        i++;
    }
    haloFile.close();

	dev_dmparts = dmparts;
	//allocating GPU memory and copy data to it
	  // Allocate GPU map.
    cudaError_t cudaStatus = cudaMalloc((void**)&dev_halos, numHalos_ * sizeof(HaloData));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating halocore memory!\n");
        return cudaStatus;
    }

	cudaStatus = cudaMemcpy(dev_halos, halos_, numHalos_ * sizeof(HaloData), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed -- copying Halocore data!\n");
        return cudaStatus;
    }

	delete halos_;

	//allocating result memory
	cudaStatus = cudaMalloc((void**)&dev_result, memNumParts_ * sizeof(MAPTYPE));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed -- allocating halocore result memory!\n");
        return cudaStatus;
    }

	return cudaStatus;
}

void clearUpHaloCore(){
	cudaFree(dev_halos);
	cudaFree(dev_result);
}

