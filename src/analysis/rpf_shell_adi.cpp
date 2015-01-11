/***********************************************
 *  * Measure the radial profile of the simulation
 *   * from a center, using shell bins
 *    * Outputs are in code units.
 *     *
 *      * Author: Lin Yang
 *       * Date: Feb/2014
 ********************************************/

//fits and HEALPIX
//measure the density field using shell method 

#include <sstream>
#include <string>
#include <cmath>
#include <cstring>


#include <rpc/types.h>
#include <rpc/xdr.h>

#include "datareader.h"
#include "datatypes.h"
#include "../halocore/halocore.h"


#ifndef PI
#define PI 3.141592653589793238462643383279
#endif


using namespace std;

int main(int argc, const char **argv){
    //printf("%d\n", argc);
    if(!(argc == 8 || argc == 7 || argc == 9)){
        printf("Usage: %s <fitsfilename>, <num of bins>, <radius>, <x>, <y>, <z> [maskfile] [-core]\n", argv[0]);
        exit(1);
    }

    int numbins;
    double radius;
    double theta, phi;
    double * databins;
    //double * angbins;
    double * countbins;
    double dr, x, y, z;
    string filename;
    string maskfile;
    bool isMask;
    bool isCore;
    
    isCore = false;
    if(argc == 7){
        isMask = false;
        isCore = false;
    }else if(argc == 8){
        //printf("%s\n", argv[7]);
        if(strstr (argv[7], "-core") != NULL){
            isMask = false;
            isCore = true;
        }else{
            isMask = true;
            isCore = false;
        }
    }else if(argc == 9){
        isMask = true;
        isCore = true;
    }
    
    stringstream ss, ss1, ss2, ss3, ss4, ss5, ss6;
    ss << argv[1];
    ss >> filename;
    ss1 << argv[2];
    ss1 >> numbins;
    ss2 << argv[3];
    ss2 >> radius;
    ss3 << argv[4];
    ss3 >> x;
    ss4 << argv[5];
    ss4 >> y;
    ss5 << argv[6];
    ss5 >> z;
    if(isMask){
        ss6 << argv[7];
        ss6 >> maskfile;
    }
    
    
    
    databins = new double[numbins];
    countbins = new double[numbins];
    
    for(int i = 0; i < numbins; i++){
        databins[i] = 0;
        countbins[i] = 0;
    }
    
    fprintf(stderr, "Filename: %s\n", filename.c_str());
    fprintf(stderr, "Bins: %d\n", numbins);
    fprintf(stderr, "Radius: %f (code units)\n", radius);
    fprintf(stderr, "x y z: %f %f %f (code units)\n", x, y, z);
    if(isMask){
        fprintf(stderr, "Mask: %s\n", maskfile.c_str());
    }
    //radius /= 40.0;
    
    dr = radius / numbins;
    
    DataReader reader(filename, isMask, maskfile);
    reader.open();

    static HaloCore halocore;

    double minsigmaav = 1.0e99;
    double maxsigmaav = 0.0;
    double totalmass = 0.0; 

    int cts = 0;
    while(reader.hasNext()){
        DMParticle * ps = reader.getBuf();
        
        for(int i = 0; i < reader.getMemparts(); i++){
            DMParticle & part = ps[i];
            if(part.sigmav > maxsigmaav){
                maxsigmaav = part.sigmav;
            }
            if(part.sigmav < minsigmaav){
                minsigmaav = part.sigmav;
            }
            
            double r = sqrt((part.posx - x) *(part.posx - x) +
                            (part.posy - y) *(part.posy - y) +
                            (part.posz - z) *(part.posz - z));
            //r *= (40000.0);
            if((r < radius) && (part.dens > 0)){
                double corr = 1.0;
                if(isCore){
                    corr = halocore.getCorrection(part.posx, part.posy, part.posz);
                }
                
                corr = (14.6534 * pow(1 + 0.0355872 * r, 1.76))/
                       (pow(1 + 0.10888 * pow(r, 0.76), 3.34211) * pow(r, 0.13));



                int ind = r / dr;
                databins[ind] += corr * part.mass / (4 * PI / 3
                                              * (pow((ind+1.0) * dr, 3) - pow((ind) * dr, 3)));//part.dens * corr;
                countbins[ind] ++;
                totalmass += part.mass * corr;
            }
            cts ++;
        }
        reader.loadBuffer();
        //printf("%d\n", cts);
    }
    //printf("\n");
    
    fprintf(stderr, "Sigma V range: [%f %e], total mass: %10e\n", minsigmaav, maxsigmaav, totalmass);
    
    reader.close();
    
    for(int i = 0; i < numbins; i++){
        if(countbins[i] != 0){
            //databins[i] /= countbins[i];
        }else{
            databins[i] = 0.0;
        }
        
        printf("%f %e\n", i * dr, databins[i]);
    }
    
    delete[] databins;
    delete[] countbins;
    
    return 0;
}
