#include <string>
#include <cstring>
#include <cmath>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <vector>
#include <string>
#include "../tipsydefs.h"
#include "../datareader.h"
#include "vl2read.h"

using namespace std;

const int VL2Reader::INITMEMSIZE = 128 * 1024;

VL2Reader::VL2Reader(std::string filename, double xc,
                     double yc, double zc, double boxsize,
                     std::string maskfile){
    
    loadParts("", filename,  xc,  yc,  zc,  boxsize, maskfile, true, false);
}

VL2Reader::VL2Reader(std::string basename, std::string filename,
                     double xc, double yc, double zc, double boxsize,
                     std::string maskfile, bool isSingleFile){
    loadParts(basename, filename,  xc,  yc,  zc,
              boxsize, maskfile, isSingleFile, true);
}

void VL2Reader::loadParts(std::string basename, std::string filename, double xc, double yc, double zc, double boxsize, std::string maskfile, bool isSingleFile, bool isTipsy){
    
    density_.reserve(INITMEMSIZE);
    velX_.reserve(INITMEMSIZE);
    velY_.reserve(INITMEMSIZE);
    velZ_.reserve(INITMEMSIZE);
    //return;

    bool isMask = false;
    if(maskfile != ""){
        isMask = true;
    }
    //return;    
    DataReader *preader;
    if(!isTipsy){
        preader = new DataReader(filename, isMask, maskfile);
    }else{
        preader = new DataReader(basename, filename, isMask, maskfile, isSingleFile);
    }
    DataReader &reader = *preader;
    reader.open();
    
    //printf("%s %d %d\n", filename.c_str(), isMask, reader.hasNext());
    
    while(reader.hasNext()){
        //printf("ok\n");
        DMParticle * ps = reader.getBuf();
        for(int i = 0; i < reader.getMemparts(); i++){
            DMParticle & part = ps[i];
            if((abs(part.posx - xc) < boxsize / 2) &&
               (abs(part.posy - yc) < boxsize / 2) &&
               (abs(part.posz - zc) < boxsize / 2)){
                
                density_.push_back(part.dens);
                velX_.push_back(part.velx);
                velY_.push_back(part.vely);
                velZ_.push_back(part.velz);
            }
        }
        reader.loadBuffer();
    }
    reader.close();

    
}



std::vector<float> & VL2Reader::getDensity(){
    return density_;
}
std::vector<float> & VL2Reader::getVelocityX(){
    return velX_;
}
std::vector<float> & VL2Reader::getVelocityY(){
    return velY_;
}
std::vector<float> & VL2Reader::getVelocityZ(){
    return velZ_;
}


//test
/*
int main(){
    VL2Reader reader("/data1/data/data_float_all.bin", 3.612030e-01, 2.106580e-01, 6.877260e-03, 0.00001);
    for(int i = 0; i < reader.getDensity().size(); i++){
        printf("%f\n", reader.getDensity()[i]);
    }
}*/
