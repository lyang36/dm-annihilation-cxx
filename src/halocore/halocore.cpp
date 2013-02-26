#include <cmath>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "../datatypes.h"
#include "halocore.h"

using namespace std;

#define SAT_RADIUS (1.0/40000000.0)

HaloCore::HaloCore(){
    
}

MAPTYPE HaloCore::getCorrection(float x, float y, float z){
    const MAPTYPE xc = 3.612030e-01;
    const MAPTYPE yc = 2.106580e-01;
    const MAPTYPE zc = 6.877260e-03;
    MAPTYPE r = coreFunc(x, y, z, xc, yc, zc, 20.0/40000.0);
    return r;
}


MAPTYPE HaloCore::coreFunc(float x, float y, float z, 
                float hxc, float hyc, float hzc, float radius){
    MAPTYPE xx = (x - hxc);
    MAPTYPE yy = (y - hyc);
    MAPTYPE zz = (z - hzc);
    MAPTYPE r = sqrt(xx*xx + yy*yy + zz*zz);
    MAPTYPE ratio = 1.0;
    if( r < radius){
        if( r > SAT_RADIUS)
            ratio = pow(r / radius, 0.6);
        else
            ratio = pow(r / SAT_RADIUS, 0.6);
    }
    return ratio;
}

