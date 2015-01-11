//SHOW the data in ascii format, output at stdout

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

#include "../datareader.h"
#include "../datatypes.h"


#ifndef PI
#define PI 3.141592653589793238462643383279
#endif


using namespace std;

void printUsage(const char * s){
    fprintf(stderr, "Usage:\n"
           "%s\n"
           "-b <baseName> for tipsyfile only\n"
           "-f <fitsFilename>\n"
           "-box <xc> <yc> <zc> <size>\n"
           "-m <maskFile>\n"
           "Output the data in ascii format. \n"
           , s);
}


int main(int argc, const char **argv){
    int numParts = 0;           //number of particles to be sampled
    string filename = "";
    string maskfile = "";
    string baseName = "";
    string outputfilename = "";
    
    
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

    
    DataReader *preader;
    if(!isTipsy){
        preader = new DataReader(filename, isMask, maskfile);
    }else{
        preader = new DataReader(baseName, filename, isMask, maskfile, true);
    }
    DataReader &reader = *preader;
    
    reader.open();

    int totalNumParts = reader.getPartNum();
    while(reader.hasNext()){
        DMParticle * ps = reader.getBuf();
        for(int i = 0; i < reader.getMemparts(); i++){
            DMParticle & part = ps[i];
            if((abs(part.posx - xc) < boxsize / 2) &&
               (abs(part.posy - yc) < boxsize / 2) &&
               (abs(part.posz - zc) < boxsize / 2)){
                
                double dens = part.dens;
                fprintf(stdout, "%0.9e %0.9e %0.9e %0.9e\n", part.dens, part.velx, part.vely, part.velz);
            }
        }
        reader.loadBuffer();
        //printf("%d\n", cts);
    }
    reader.close();
    return 0;
}
