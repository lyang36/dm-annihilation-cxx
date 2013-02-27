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

using namespace std;

string index_file =  "/home/lyang/data/vl2b.00400.r200.index";
string ahf_part_file = "/home/lyang/halodata/vl_400_rhovesc.z0.000.AHF_particles";
string output_file = "vl2b.00400.r200.ahf.haloflags";

int * sorted_key_;

//get flags
void getFlag(int * particles_, char * flags_, int numParts_){
    int numHalos = 0;    
    thrust::device_vector<int> dev_key(particles_, particles_ + numParts_);    

    //sorted_key_ = new int[numParts_];
    sorted_key_ = particles_;

    for(int i = 0; i < numParts_; i++){
        flags_[i] = 0;
        sorted_key_[i] = i;
    }

    //thrust::device_vector<int> dev_key(particles_,  particles_+ numParts_);
    printf("ok2.5\n");
    thrust::device_vector<int> dev_val(sorted_key_, sorted_key_ + numParts_);
    printf("ok3\n");    

    thrust::sort_by_key(dev_key.begin(), dev_key.end(), dev_val.begin());
    printf("ok4\n");
    thrust::copy(dev_val.begin(), dev_val.end(), sorted_key_);

    ifstream haloInputFile_(ahf_part_file.c_str());
    haloInputFile_ >> numHalos;
    for(int i = 0; i < numHalos; i ++){
        int numHaloParts;
        haloInputFile_ >> numHaloParts;
        for(int j = 0; j < numHaloParts; j++){
            int partindex;
            haloInputFile_ >> partindex;
            if(i >= IGNORE_FIRST_N){
                thrust::pair<thrust::device_vector<int>::iterator, thrust::device_vector<int>::iterator> ret
                     = thrust::equal_range(dev_key.begin(), dev_key.end(), partindex);
                if(ret.first != ret.second){
                    int ind = (ret.first - dev_key.begin());
                    flags_[sorted_key_[ind]] = 1;
                }
            }
        }
        if(i % 100 == 0){
            printf(".");
        }
    }
    printf("\n");
    //delete sorted_key_;
    haloInputFile_.close();
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
    printf("ok\n");
    flags_ = new char[numParts_];
    printf("ok1\n");    

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
