
#ifndef __LY__MAP__KERNEL__
#define __LY__MAP__KERNEL__
#include "datatypes.h"

#define BLOCK_N_DIVIDER 16			//each block has 16*16=256 threads 



cudaError_t initializeCUDA(int Nside, int memParts);
cudaError_t calulateMapWithCUDA(MAPTYPE * map, DMParticle * parts, int numParts);
void cudaCleaingUp();
void setCudaPixelList(vector<int> * parPixList);

//get the pointer of the particles in GPU -- basically used for add core
DMParticle * getDevParticle();

#endif
