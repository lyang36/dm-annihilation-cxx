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

#define GPU_MEM 100000000

using namespace std;

//the particle index in the r200 file
bool isIndexBin = true;
string index_file_txt =  "";//"/home/lyang/data/vl2b.00400.r200.index";
string index_file_bin =  "";// "/home/lyang/data/vl2b.00400.r200.index";

//the particle file in the AHF particles
bool isAHFPartFileBin = false;
string ahf_part_file_txt = "";//"/home/lyang/halodata/vl_400_rhovesc.z0.000.AHF_particles";
string ahf_part_file_bin = "";//"/home/lyang/halodata/vl_400_rhovesc.z0.000.AHF_particles";

//the output halo flags (a binary file)
string output_file = "";//"vl2b.00400.r200.ahf.haloflags";

//halo ids to be selected
//the first int is a number of total number to be selected
//followed by the list of halo ids
bool isHaloIdsBin = false;
string haloids_to_be_selected_bin = "to_be_seleted.ids";
string haloids_to_be_selected_txt = "";

int * haloParticles_;
int * haloIds_;
int * searchParts_;
int * searchIndex_;
bool * searchResult_;
char * flags_;
int * particles_;

//bool verbose = false;
int numParts_ = 0;
int numOfHalos_ = 0;

void getSearchRes(int requiredSearchPartNum, int numPartsRead_,
                  thrust::device_vector<int> &dev_searchParts_,
                  thrust::device_vector<int> &dev_searchResult_,
                  thrust::device_vector<int> &dev_val){
    //do the search
    thrust::copy(searchParts_, searchParts_ + requiredSearchPartNum, dev_searchParts_.begin());
    thrust::binary_search(dev_val.begin(), dev_val.begin() + numPartsRead_,
                          dev_searchParts_.begin(),
                          dev_searchParts_.begin() + requiredSearchPartNum,
                          dev_searchResult_.begin());
    thrust::copy(dev_searchResult_.begin(), dev_searchResult_.begin() + requiredSearchPartNum, searchResult_);
    for(int l = 0; l < requiredSearchPartNum; l++){
        if(searchResult_[l]){
            flags_[searchIndex_[l]] = 1;
        }
    }
}

void doSearch(int numPartsRead_, thrust::device_vector<int> &dev_searchParts_,
              thrust::device_vector<int> &dev_searchResult_,
              thrust::device_vector<int> &dev_val){
    printf("Start testing %d halo particles...\n", numPartsRead_);
    //start filling the tags
    //step 1: sorting
    printf("Sorting ...\n");
    thrust::copy(haloParticles_, haloParticles_ + numPartsRead_, dev_val.begin());
    thrust::sort(dev_val.begin(), dev_val.begin() + numPartsRead_);
    
    //step 2: testing
    printf("Searching ...\n");
    //test every particle whether it's in the array
    int requiredSearchPartNum = 0;
    for(int k = 0; k < numParts_; k ++){
        if(flags_[k] == 0){
            searchParts_[requiredSearchPartNum] = particles_[k];
            searchIndex_[requiredSearchPartNum] = k;
            requiredSearchPartNum ++;
        }
        if(requiredSearchPartNum >= GPU_MEM){
            getSearchRes(requiredSearchPartNum, numPartsRead_,
                         dev_searchParts_, dev_searchResult_, dev_val);
            requiredSearchPartNum = 0;
        }
    }
    if(requiredSearchPartNum > 0){
        getSearchRes(requiredSearchPartNum, numPartsRead_,
                     dev_searchParts_, dev_searchResult_, dev_val);
        requiredSearchPartNum = 0;
    }

}


//get flags
void getFlag(){
    
    thrust::device_vector<int> dev_searchParts_(GPU_MEM);
    thrust::device_vector<int> dev_searchResult_(GPU_MEM);
    thrust::device_vector<int> dev_val(GPU_MEM);
    thrust::device_vector<int> dev_searchHaloIds_(numOfHalos_);
    
    //thrust::binary_search(input.begin(), input.end(), 0, thrust::less<int>()); // returns true
    thrust::copy(haloIds_, haloIds_ + numOfHalos_, dev_searchHaloIds_.begin());
    thrust::sort(dev_searchHaloIds_.begin(), dev_searchHaloIds_.begin() + numOfHalos_);
    
    
    //int currentHalo = 0;
    int totalNumHalos = 0;
    for(int i = 0; i < numParts_; i++){
        flags_[i] = 0;
    }
    
    haloParticles_ = new int[GPU_MEM];
    searchParts_ = new int[GPU_MEM];
    searchIndex_ = new int[GPU_MEM];
    searchResult_ = new bool[GPU_MEM];
    

    ifstream haloInputFile_;
    if(isAHFPartFileBin){
        haloInputFile_.open(ahf_part_file_bin.c_str(), ios::binary | ios::in);
        haloInputFile_.read((char *)&totalNumHalos, sizeof(int));
    }else{
        haloInputFile_.open(ahf_part_file_txt.c_str());
        haloInputFile_ >> totalNumHalos;
    }
    
    if(!haloInputFile_.good()){
        printf("AHF Particle File Error!\n");
        exit(1);
    }
    
    int numPartsRead_ = 0;
    printf("Start reading halo particles...\n", numOfHalos_);
    

    for(int i = 0; i < totalNumHalos; i ++){
        int numHaloParts;
        if(isAHFPartFileBin){
            haloInputFile_.read((char *) &numHaloParts, sizeof(int));
        }else{
            haloInputFile_ >> numHaloParts;
        }

        
        if(thrust::binary_search(dev_searchHaloIds_.begin(),
                                 dev_searchHaloIds_.end(),
                                 i,
                                 thrust::less<int>())){
            printf("Halo: %d, Particles: %d.\n", i, numHaloParts);
            
            for(int j = 0; j < numHaloParts; j++){
                int partindex;
                int ch;
                
                if(isAHFPartFileBin){
                    haloInputFile_.read((char *) &partindex, sizeof(int));
                    haloInputFile_.read((char *) &ch, sizeof(int));
                }else{
                    haloInputFile_ >> partindex;
                    haloInputFile_ >> ch;
                }
                haloParticles_[numPartsRead_] = partindex;
                numPartsRead_ ++;
                
                if(numPartsRead_ >= GPU_MEM){
                    doSearch(numPartsRead_, dev_searchParts_, dev_searchResult_, dev_val);
                    numPartsRead_ = 0;
                }
            }
        }else{
            string line;
            if(isAHFPartFileBin){
                haloInputFile_.seekg(sizeof(int) * numHaloParts, ios_base::cur);
            }else{
                for(int j = 0; j < numHaloParts; j++){
                    //haloInputFile_ >> partindex;
                    //haloInputFile_ >> ch;
                    getline(haloInputFile_, line);
                }
            }
            
        }

    }
    
    if(numPartsRead_ > 0){
        doSearch(numPartsRead_, dev_searchParts_, dev_searchResult_, dev_val);
        numPartsRead_ = 0;
    }
    printf("\n");
    haloInputFile_.close();
    
    delete haloParticles_;
    delete searchParts_;
    delete searchIndex_;
    delete searchResult_;
}

void printUsage(const char * name){
   printf("%s \n%s \n%s \n%s \n%s\n", name,
	"-[bin/txt]_index <index file>",
	"-[bin/txt]_ahf <AHF particle output file>",
	"-[bin/txt]_haloid <ids of halo to be selected>",
	"-output <outputfile>");
}

int main(int argc, const char **argv){
    int m=1;

    if(argc < 9){
    	printUsage(argv[0]);
	exit(1);
    } 
    while (m<argc)
    {
        string arg = argv[m];
        if (arg == "-bin_index") {
            isIndexBin = true;
            index_file_bin = argv[m+1];
            m+=1;
        }else if (arg == "-txt_index") {
            isIndexBin = false;
            index_file_txt = argv[m+1];
            m+=1;
        }else if (arg == "-bin_ahf") {
            isAHFPartFileBin = true;
            ahf_part_file_bin = argv[m+1];
            m+=1;
        }else if (arg == "-txt_ahf") {
            isAHFPartFileBin = false;
            ahf_part_file_txt = argv[m+1];
            m+=1;
        }else if (arg == "-output") {
            output_file = argv[m+1];
            m+=1;
        }else if (arg == "-bin_haloid") {
            isHaloIdsBin = true;
            haloids_to_be_selected_bin = argv[m+1];
            m+=1;
        }else if (arg == "-txt_haloid") {
            isHaloIdsBin = false;
            haloids_to_be_selected_txt = argv[m+1];
            m+=1;
        }
        //else if (arg == "-verbose") { verbose = true;}
        else{
            printUsage(argv[0]);
            exit(0);
        }
        m++;
    }
    
    ifstream dataInputFile_;
    if(isIndexBin){
        dataInputFile_.open(index_file_bin.c_str(), ios::binary | ios::in);
        if(!dataInputFile_.good()){
            printf("Datafile error: %s !\n", index_file_bin.c_str());
            exit(1);
        }
        
        dataInputFile_.read((char*)&numParts_, sizeof(int));
    }else{
        dataInputFile_.open(index_file_txt.c_str(), ios::in);
        if(!dataInputFile_.good()){
            printf("Datafile error: %s !\n", index_file_txt.c_str());
            exit(1);
        }
        
        dataInputFile_ >>numParts_;
    }
    cout << "Particles: " << numParts_ << endl;
    
    particles_ = new int[numParts_];
    //printf("ok\n");
    flags_ = new char[numParts_];
    //printf("ok1\n");
    if(isIndexBin){
        dataInputFile_.read((char *) particles_, sizeof(int) * numParts_);
        dataInputFile_.close();
    }else{
        for(int i = 0; i < numParts_; i++){
            dataInputFile_ >> particles_[i];
        }
        dataInputFile_.close();
    }
  
    ifstream haloidsStream_;
    //printf("%d %s\n", isHaloIdsBin, )
    if(isHaloIdsBin){
        haloidsStream_.open(haloids_to_be_selected_bin.c_str(), ios::binary);
        if(!haloidsStream_.good()){
            printf("Halo Id error: %s!\n", haloids_to_be_selected_bin.c_str());
            exit(1);
        }
        haloidsStream_.read((char *) &numOfHalos_, sizeof(int));
    }else{
        haloidsStream_.open(haloids_to_be_selected_txt.c_str());
        if(!haloidsStream_.good()){
            printf("Halo Id error: %s!\n", haloids_to_be_selected_txt.c_str());
            exit(1);
        }
        haloidsStream_ >> numOfHalos_;
    }
    
    printf("Number of Halos: %d\n", numOfHalos_);
    haloIds_ = new int[numOfHalos_];
    
    if(isHaloIdsBin){
        haloidsStream_.read((char *) haloIds_, sizeof(int) * numOfHalos_);
        haloidsStream_.close();
    }else{
        for(int i = 0; i < numOfHalos_; i++){
            haloidsStream_ >> haloIds_[i];
        }
        haloidsStream_.close();
    }
    
    getFlag();
    
    //output
    printf("Output the result...\n");
    ofstream dataOutputStream_(output_file.c_str(), ios::binary);
    dataOutputStream_.write((char *) &numParts_, sizeof(int));
    dataOutputStream_.write((char *) flags_, sizeof(char) * numParts_);
    dataOutputStream_.close();
    printf("Finished...\n");
    
    delete particles_;
    delete flags_;
    delete haloIds_;
}
