/***********************************************
 * Measure the radial profile of the simulation
 * from a center, using shell bins
 * Outputs are in code units.
 *
 * Author: Lin Yang
 * Date: Feb/2014
 * ********************************************/



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

void printUsage(const char * s){
    printf("Usage:\n"
           "%s\n"
           "-b <baseName> for tipsyfile only\n"
           "-f <fitsfilename>\n"
           "-n <num of bins>\n"
           "-r <radius (code units)>\n"
           "-c <x>, <y>, <z> (code units)\n"
           "-m <maskfile>\n"
           "-core\n", s);
}


int main(int argc, const char **argv){
    int numbins = 100;
    double radius;
    double theta, phi;
    double * databins;
    //double * angbins;
    double * countbins;
    double dr, x, y, z;
    string baseName = "";
    string filename = "";
    string maskfile = "";
    bool isMask = false;
    bool isCore = false;
    bool isTipsy = false;
    
    int m=1;
    //printf("ok!\n");
    while (m<argc)
    {
        string arg = argv[m];
        //printf("ok!\n");
        stringstream ss, ss1, ss2;
        if (arg == "-b") {
            isTipsy = true;
            baseName = argv[m+1];
            m+=1;
        }else if (arg == "-f") {
            filename = argv[m+1];
            m+=1;
        }else if (arg == "-n") {
            ss << argv[m+1];
            ss >> numbins;
            m+=1;
        }else if (arg == "-r") {
            ss << argv[m+1];
            ss >> radius;
            m+=1;
        }else if (arg == "-c") {
            ss << argv[m+1];
            ss >> x;
            m++;
            
            ss1 << argv[m+1];
            ss1 >> y;
            m++;
            
            ss2 << argv[m+1];
            ss2 >> z;
            m++;
        }else if (arg == "-m") {
            isMask = true;
            maskfile = argv[m+1];
            m+=1;
        }else if (arg == "-core") {
            isCore = true;
        }else{
            printUsage(argv[0]);
            exit(1);
        }
        m++;
    }
   
   
    if(filename == ""){
        printUsage(argv[0]);
        exit(1);
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
    
    DataReader *preader;
    if(!isTipsy){
        preader = new DataReader(filename, isMask, maskfile);
    }else{
        preader = new DataReader(baseName, filename, isMask, maskfile, true);
    }
    
    DataReader &reader = *preader;
    reader.open();

    HaloCore * phalocore;
    if(isCore){
        phalocore = new HaloCore();
    }

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
            if((r < radius)){// && (part.dens > 0)){
                double corr = 1.0;
                // Correlation
                if(isCore){
                    corr = phalocore->getCorrection(part.posx, part.posy, part.posz);
                }
                int ind = r / dr;
                databins[ind] += part.mass / (4 * PI / 3
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
    delete phalocore;
    
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
    delete preader;

    return 0;
}
