#ifndef __DATATYPE__LY
#define __DATATYPE__LY
//the particle for visualization
#include "tipsydefs.h"

class DMParticle {
public:
    float mass;
    float dens;
    float hsmooth;
    float posx;
    float posy;
    float posz;
    float velx;
    float vely;
    float velz;
    float eps;
    float phi;
};

#endif

