// Measure the radial bulk velocity profile
// Definition:
// For each shell
// v = \sum(m_i*(v_i \cdot \hat{r})) / \sum{m_i}

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
           "-r <radius>\n"
           "-c <x>, <y>, <z>\n"
           "-m <maskfile>\n"
           "-core\n", s);
}


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
    double * massbins;
    double * densbins;
    double dr, x, y, z;
    string filename = "";
    string maskfile = "";
    string baseName = "";
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
    massbins = new double[numbins];
    countbins = new double[numbins];
    densbins = new double[numbins];
    
    for(int i = 0; i < numbins; i++){
        databins[i] = 0;
        countbins[i] = 0;
        massbins[i] = 0;
        densbins[i] = 0;
    }
    
    fprintf(stderr, "Filename: %s\n", filename.c_str());
    fprintf(stderr, "Bins: %d\n", numbins);
    fprintf(stderr, "Radius: %f (kpc/h)\n", radius);
    fprintf(stderr, "x y z: %f %f %f (40 Mpc/h)\n", x, y, z);
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
            
            double rx = part.posx - x;
            double ry = part.posy - y;
            double rz = part.posz - z;
            double r = sqrt(rx * rx + ry * ry + rz * rz);
            rx /= r;
            ry /= r;
            rz /= r;
            
            r *= (40000.0);
            
            if((r < radius) && (part.dens > 0)){
                double corr = 1.0;
                if(isCore){
                    corr = halocore.getCorrection(part.posx, part.posy, part.posz);
                }
                
                int ind = r / dr;
                databins[ind] += part.mass *
                        (part.velx * rx + part.vely * ry + part.velz * rz);
                
                countbins[ind] ++;
                massbins[ind] += part.mass;
                densbins[ind] += part.mass / (4 * PI / 3
                                              * (pow((ind+1.0) * dr, 3) - pow((ind) * dr, 3)));
                
                totalmass += part.mass * corr;
            }
            cts ++;
            //if(cts > 1000000)
            //      break;
        }

        //test
        //break;

        reader.loadBuffer();
        //printf("%d\n", cts);
    }
    //printf("\n");
    
    fprintf(stderr, "Sigma V range: [%f %f], total mass: %10e\n", minsigmaav, maxsigmaav, totalmass);
    
    reader.close();
    
    for(int i = 0; i < numbins; i++){
        if(massbins[i] != 0){
            databins[i] /= massbins[i];
            //densbins[i] = massbins[i] / (4 * PI / 3) / (pow(i * dr + dr, 3) - pow(i * dr, 3));
        }else{
            databins[i] = 0.0;
            //densbins[i] = 0.0;
        }
        
        printf("%f %e %e %e %d\n", i * dr, databins[i], densbins[i], massbins[i], (int)countbins[i]);
    }
    
    delete[] databins;
    delete[] countbins;
    delete[] massbins;
    delete[] densbins;
    
    return 0;
}
