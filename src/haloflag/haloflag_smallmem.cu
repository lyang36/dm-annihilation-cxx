#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/binary_search.h>
#include <thrust/pair.h>

#define IGNORE_FIRST_N 1        //ignore the first n halo? host?
#define IGNORE_LAST_N 100000000        //ignore the halo after this number?
#define GPU_MEM 10000000

using namespace std;

string index_file =  "/home/lyang/data/vl2b.00400.r200.index";
string ahf_part_file = "/home/lyang/halodata/vl_400_rhovesc.z0.000.AHF_particles";
string output_file = "vl2b.00400.r200.ahf.haloflags";



//get flags
void getFlag(int * particles_, char * flags_, int numParts_){
    int numHalos = 0;
    for(int i = 0; i < numParts_; i++){
        flags_[i] = 0;
    }
    
    thrust::device_vector<int> dev_val(GPU_MEM);
    int * haloParticles_ = new int[GPU_MEM];
    
    ifstream haloInputFile_(ahf_part_file.c_str());
    haloInputFile_ >> numHalos;
    int numPartsRead_ = 0;
    for(int i = 0; (i < numHalos) && (i <= IGNORE_LAST_N); i ++){
        int numHaloParts;
        haloInputFile_ >> numHaloParts;
        for(int j = 0; j < numHaloParts; j++){
            int partindex;
            haloInputFile_ >> partindex;
            if(i >= IGNORE_FIRST_N){
                haloParticles_[numPartsRead_] = partindex;
                numPartsRead_ ++;
            }
            
            if(numPartsRead_ >= GPU_MEM){
                printf("Start testing %n halo particles...\n", numPartsRead_);
                //start filling the tags
                //step 1: sorting
                thrust::copy(haloParticles_, haloParticles_ + numPartsRead_, dev_val.begin());
                thrust::sort(dev_val.begin(), dev_val.end());
                
                //step 2: testing
                //test every particle whether it's in the array
                for(int k = 0; k < numParts_; k ++){
                    if(flags_[k] == 0){
                        if(thrust::binary_search(dev_val.begin(), dev_val.end(), particles_[k])){
                            flags_[k] = 1;
                        }
                    }
                }
                numPartsRead_ = 0;
            }
        }
    }
    
    if(numPartsRead_ > 0){
        printf("Start testing %d halo particles...\n", numPartsRead_);
        //start filling the tags
        //step 1: sorting
        thrust::copy(haloParticles_, haloParticles_ + numPartsRead_, dev_val.begin());
        thrust::sort(dev_val.begin(), dev_val.end());
        
        //step 2: testing
        //test every particle whether it's in the array
        for(int k = 0; k < numParts_; k ++){
            if(flags_[k] == 0){
                if(thrust::binary_search(dev_val.begin(), dev_val.end(), particles_[k])){
                    flags_[k] = 1;
                }
            }
        }
        numPartsRead_ = 0;
    }
    printf("\n");
    haloInputFile_.close();
    
    delete haloParticles_;
}


int main(int argc, const char **argv){
    int m=1;
    //bool verbose = false;
    int * particles_;
    char * flags_;
    int numParts_ = 0;
    
    while (m<argc)
    {
        string arg = argv[m];
        if (arg == "-index") { index_file = argv[m+1]; m+=1;}
        if (arg == "-ahf") { ahf_part_file = argv[m+1]; m+=1;}
        if (arg == "-output") { output_file = argv[m+1]; m+=1;}
        //else if (arg == "-verbose") { verbose = true;}
        else{
            cout << "Usage:" << endl;
            exit(0);
        }
        m++;
    }
    
    ifstream dataInputFile_;
    dataInputFile_.open(index_file.c_str(), ios::binary);
    if(!dataInputFile_.good()){
        printf("Datafile error: %s !\n", index_file.c_str());
        exit(1);
    }
    dataInputFile_.read((char*)&numParts_, sizeof(int));
    cout << "Particles: " << numParts_ << endl;
    
    particles_ = new int[numParts_];
    //printf("ok\n");
    flags_ = new char[numParts_];
    //printf("ok1\n");
    
    dataInputFile_.read((char *) particles_, sizeof(int) * numParts_);
    dataInputFile_.close();
    
    getFlag(particles_, flags_, numParts_);
    
    //output
    printf("Output the result...\n");
    ofstream dataOutputStream_(output_file.c_str(), ios::binary);
    dataOutputStream_.write((char *) &numParts_, sizeof(int));
    dataOutputStream_.write((char *) flags_, sizeof(char) * numParts_);
    dataOutputStream_.close();
    printf("Finished...\n");
    
    delete particles_;
    delete flags_;
}
