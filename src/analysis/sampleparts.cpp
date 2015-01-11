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
    printf("Usage:\n"
           "%s\n"
           "-b <baseName> for tipsyfile only\n"
           "-f <fitsFilename>\n"
           "-n <numberOfParts>\n"
           "-m <maskFile>\n"
           "-o <outputFilename>\n"
           "Output the sampleling particles, with dens=-1 if masked out\n"
           , s);
}


int main(int argc, const char **argv){
    int numParts = 0;           //number of particles to be sampled
    string filename = "";
    string maskfile = "";
    string baseName = "";
    string outputfilename = "";
    bool isMask = false;
    bool isTipsy = false;
    
    int m=1;
    //printf("ok!\n");
    while (m<argc)
    {
        string arg = argv[m];
        //printf("ok!\n");
        stringstream ss;
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

    
    fprintf(stderr, "Filename: %s\n", filename.c_str());
    fprintf(stderr, "Sampling Parts: %d\n", numParts);
    fprintf(stderr, "Output file: %s\n", outputfilename.c_str());

    if(isMask){
        fprintf(stderr, "Mask: %s\n", maskfile.c_str());
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
    double rejectingprob = ((double) numParts) / ((double) totalNumParts);
    
    reader.close();
    delete preader;
    
    int partCounts = 0;
    
    //open the output file for writting
    fstream outputStream(outputfilename.c_str(), ios::binary | ios::out);
    if(!outputStream.good()){
        printf("Cannot write output file!\n");
        exit(1);
    }
    
    //write the number of particles in the file
    outputStream.write((char *) &numParts, sizeof(int));
    
    //initialize random generater
    srand (time(NULL));
    
    //starting sampling the particles
    while(partCounts < numParts){
        
        if(!isTipsy){
            preader = new DataReader(filename, isMask, maskfile);
        }else{
            preader = new DataReader(baseName, filename, isMask, maskfile, true);
        }
        DataReader &reader = *preader;
        reader.open();

        
        while(reader.hasNext()){
            DMParticle * ps = reader.getBuf();
            for(int i = 0; i < reader.getMemparts(); i++){
                double prob = ((double) rand()) / ((double) RAND_MAX);
                if(prob < rejectingprob){
                    outputStream.write((char *) &(ps[i]), sizeof(DMParticle));
                    //printf("Parts: %f, %f, %d\n", prob, rejectingprob, i);
                    partCounts ++;
                    if(partCounts >= numParts){
                        goto SAMPLINGOVER;
                    }
                }
            }
            reader.loadBuffer();
        }
    SAMPLINGOVER:
        reader.close();
        delete preader;
    }
    
    printf("Finished!\n");
    outputStream.close();
    return 0;
}
