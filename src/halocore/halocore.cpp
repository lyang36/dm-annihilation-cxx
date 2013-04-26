#include <cmath>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <fstream>
#include "../datatypes.h"
#include "halocore.h"

using namespace std;

#define RADIUS_RATIO (4.0 / 402.0)
#define SAT_RADIUS (1.0 / 40000.0 / RADIUS_RATIO)

static const string halo_core_file_name_ = "./halocore/VL2_Halos.data";

HaloCore::HaloCore(){
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
    printf("There are %d halos in consideration!\n", numHalos_);

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
}

HaloCore::~HaloCore(){
    delete halos_;
}

MAPTYPE HaloCore::getCorrection(float x, float y, float z){
    MAPTYPE xp = x;
    MAPTYPE yp = y;
    MAPTYPE zp = z;

    MAPTYPE r = 1.0;
    for(int i = 0; i < numHalos_; i++){
        r *= coreFunc(xp, yp, zp, halos_[i].xc, 
                        halos_[i].yc,
                        halos_[i].zc, 
                        halos_[i].radius * RADIUS_RATIO / 40000.0);
        /*printf("%f %f %f %f %f %f %f\n", xp, yp, zp,
                        halos_[i].xc,
                        halos_[i].yc,
                        halos_[i].zc,
                        halos_[i].radius * RADIUS_RATIO);*/
    }
    //printf("%f\n", r);   
    return r;
}


MAPTYPE HaloCore::coreFunc(float x, float y, float z, 
                float hxc, float hyc, float hzc, float radius){
    MAPTYPE xx = (x - hxc);
    MAPTYPE yy = (y - hyc);
    MAPTYPE zz = (z - hzc);
    MAPTYPE r = sqrt(xx*xx + yy*yy + zz*zz);
    MAPTYPE ratio = 1.0;
    //if(r < 0.002)
        //printf("%f %f\n", r, radius);
    if( r < radius){
        if( r > SAT_RADIUS * radius)
            ratio = pow(r / radius, 0.6);
        else
            ratio = pow(r / SAT_RADIUS / radius, 0.6);
    }
    return ratio;
}

