#ifndef __CUDA__HALO__CORE__
#define __CUDA__HALO__CORE__
#include "../datatypes.h"
#include "../kernel.h"
#include "halocore.h"

//cuda must be initialized before call any of these functions
cudaError_t initializeCore(DMParticle * dev_parts, int numParts);
//when applied the result, return back the results
cudaError_t applyingCore(DMParticle * parts, int numParts, MAPTYPE * result);

void clearUpHaloCore();
#endif
