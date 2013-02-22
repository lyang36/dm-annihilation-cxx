#ifndef __LY__MAP__KERNEL__
#define __LY__MAP__KERNEL__
#include "datatypes.h"

cudaError_t initializeCUDA(MAPTYPE * healpixX, MAPTYPE * healpixY, MAPTYPE * healpixZ, int Npix, int memParts);
cudaError_t calulateMapWithCUDA(MAPTYPE * map, DMParticle * parts, int numParts);
void cudaCleaingUp();

#endif