// Sample a number of particles from a given file with smooth data
// (the data is the same with the one given to the DM gammaray producer)

#include <sstream>
#include <string>
#include <cmath>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include <ctime>       /* time */


#include <rpc/types.h>
#include <rpc/xdr.h>

#include "datareader.h"
#include "datatypes.h"


#ifndef PI
#define PI 3.141592653589793238462643383279
#endif


using namespace std;

void printUsage(const char * s){
    fprintf(stderr, "Usage:\n"
           "%s\n"
           "-b <baseName> for tipsyfile only\n"
           "-f <fitsFilename>\n"
           "-density <min (0.0)> <max (100.0)> <binsize (1.0)>\n"
           "-velocity <min (-100.0)> <max (100.0)> <binsize 1.0>\n"
           "-box <xc> <yc> <zc> <size>\n"
           "-m <maskFile>\n"
           "-o <outputFilename>\n"
           "Output the distribution of velocity and density. \n"
           , s);
}


int main(int argc, const char **argv){
    int numParts = 0;           //number of particles to be sampled
    string filename = "";
    string maskfile = "";
    string baseName = "";
    string outputfilename = "";
    
    double mindens = 0.0;
    double maxdens = 100.0;
    double densBin = 1.0;
    
    double minvel = -100.0;
    double maxvel = 100.0;
    double velBin = 1.0;
    
    double xc, yc, zc, boxsize;
    
    bool isMask = false;
    bool isTipsy = false;

    
    int m=1;
    //printf("ok!\n");
    while (m<argc)
    {
        string arg = argv[m];
        //printf("ok!\n");
        stringstream ss, ss1, ss2, ss3;
        if (arg == "-b") {
            isTipsy = true;
            baseName = argv[m+1];
            m+=1;
        }else if (arg == "-f") {
            filename = argv[m+1];
            m+=1;
        }else if (arg == "-n") {
            ss << argv[m+1];
            ss >> numParts;
            m+=1;
        }else if (arg == "-velocity") {
            ss << argv[m+1];
            ss >> minvel;
            m++;
            ss1 << argv[m+1];
            ss1 >> maxvel;
            m++;
            ss2 << argv[m+1];
            ss2 >> velBin;
            m++;
        }else if (arg == "-density") {
            ss << argv[m+1];
            ss >> mindens;
            m++;
            ss1 << argv[m+1];
            ss1 >> maxdens;
            m++;
            ss2 << argv[m+1];
            ss2 >> densBin;
            m++;
        }else if (arg == "-box") {
            ss << argv[m+1];
            ss >> xc;
            m++;
            ss1 << argv[m+1];
            ss1 >> yc;
            m++;
            ss2 << argv[m+1];
            ss2 >> zc;
            m++;
            ss3 << argv[m+1];
            ss3 >> boxsize;
            m++;
        }else if (arg == "-m") {
            isMask = true;
            maskfile = argv[m+1];
            m+=1;
        }else if (arg == "-o") {
            outputfilename = argv[m+1];
            m+=1;
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

    
    
    int numDensBin = (maxdens - mindens) / densBin;
    int numVelXYZBin = (maxvel - minvel) / velBin;
    double maxVelMag = sqrt(pow(max(abs(maxvel), abs(minvel)), 2) * 3);
    int numVelMBin = maxVelMag / velBin;
    
    int * densBinCounts = new int[numDensBin];
    int * velXBinCounts = new int[numVelXYZBin];
    int * velYBinCounts = new int[numVelXYZBin];
    int * velZBinCounts = new int[numVelXYZBin];
    int * velMBinCounts = new int[numVelMBin];
    

    
    DataReader *preader;
    if(!isTipsy){
        preader = new DataReader(filename, isMask, maskfile);
    }else{
        preader = new DataReader(baseName, filename, isMask, maskfile, true);
    }
    DataReader &reader = *preader;
    
    reader.open();

    int totalNumParts = reader.getPartNum();
    double rejectingprob = ((double) numParts) / ((double) totalNumParts);
    
    
    fprintf(stderr, "Filename: %s\n", filename.c_str());
    fprintf(stderr, "Velocity: %f %f %f\n", minvel, maxvel, velBin);
    fprintf(stderr, "Density: %f %f %f\n", mindens, maxdens, densBin);
    fprintf(stderr, "Box: %f %f %f %f\n", xc, yc, zc, boxsize);
    fprintf(stderr, "Num Particles: %d\n", totalNumParts);
    
    if(isMask){
        fprintf(stderr, "Mask: %s\n", maskfile.c_str());
    }
    
    
    double measureMaxDens;
    double measureMaxVel, measureMinVel;
    
    //clear all the histogram to be zero
    for(int i = 0; i < numDensBin; i++){
        densBinCounts[i] = 0;
    }
    for(int i = 0; i < numVelXYZBin; i++){
        velXBinCounts[i] = velYBinCounts[i] = velZBinCounts[i] = 0;
    }
    for(int i = 0; i < numVelMBin; i++){
        velMBinCounts[i] = 0;
    }
    
    int cts = 0;
    while(reader.hasNext()){
        DMParticle * ps = reader.getBuf();
        for(int i = 0; i < reader.getMemparts(); i++){
            DMParticle & part = ps[i];
            if((abs(part.posx - xc) < boxsize / 2) &&
               (abs(part.posy - yc) < boxsize / 2) &&
               (abs(part.posz - zc) < boxsize / 2)){
                
                double dens = part.dens;
                if((dens < maxdens) && (dens >= mindens)){
                    int densInd = (dens - mindens) / densBin;
                    densBinCounts[densInd] ++;
                }
                
                if((part.velx < maxvel) && (part.velx >= minvel)){
                    int velInd = (part.velx - minvel) / velBin;
                    velXBinCounts[velInd] ++;
                }
                
                if((part.vely < maxvel) && (part.vely >= minvel)){
                    int velInd = (part.vely - minvel) / velBin;
                    velYBinCounts[velInd] ++;
                }
                
                if((part.velz < maxvel) && (part.velz >= minvel)){
                    int velInd = (part.velz - minvel) / velBin;
                    velZBinCounts[velInd] ++;
                }
                
                double vel = sqrt(part.velx * part.velx + part.vely * part.vely + part.velz);
                
                if((vel < maxVelMag) && (vel >= 0)){
                    int velInd = (vel) / velBin;
                    velMBinCounts[velInd] ++;
                }
                cts ++;
            }
        }
        reader.loadBuffer();
        //printf("%d\n", cts);
    }
    reader.close();
    
    //output the results:
    fprintf(stdout, "Density, Counts\n");
    for(int i = 0; i < numDensBin; i++){
        fprintf(stdout, "%f %d\n", i * densBin + mindens, densBinCounts[i]);
    }
    
    fprintf(stdout, "\nVelocity, VelX Counts, VelY Counts, VelZ Counts\n");
    for(int i = 0; i < numVelXYZBin; i++){
        fprintf(stdout, "%f %d %d %d\n", i * velBin + minvel,
                velXBinCounts[i],
                velYBinCounts[i],
                velZBinCounts[i]);
    }
    
    fprintf(stdout, "\nVelocity Magnitude, Counts\n");
    for(int i = 0; i < numVelMBin; i++){
        fprintf(stdout, "%f %d\n", i * velBin, velMBinCounts[i]);
    }
    
    delete[] densBinCounts;
    delete[] velXBinCounts;
    delete[] velYBinCounts;
    delete[] velZBinCounts;
    delete[] velMBinCounts;
    return 0;
}
