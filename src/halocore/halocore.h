#ifndef __HALO__CORE__
#define __HALO__CORE__
//input a list of halo's x, y, z, radius and core function, calculate flux correction of each particle
#include "../datatypes.h"

class HaloData{
public:
    float xc;
    float yc;
    float zc;
    float mass;
    float radius;
};

class HaloCore{
public:
    HaloCore(std::string halofile = "./halocore/VL2_Halos.data");
    ~HaloCore();
    MAPTYPE getCorrection(float x, float y, float z);
private:
    std::string halo_core_file_name_;
    std::string haloDataFiles_; 
    MAPTYPE coreFunc(float x, float y, float z, 
                    float hxc, float hyc, float hzc, float radius);
    HaloData * halos_;
    int numHalos_;
};

#endif
