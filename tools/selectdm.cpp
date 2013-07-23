#include <string>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <algorithm>    // std::sort

#include <rpc/types.h>
#include <rpc/xdr.h>

#include "../src/tipsydefs.h"

using namespace std;

int main(int argn, char ** argv){
    bool istext = false;
    if(argn != 4 && argn !=5){
        printf("Usage: selectpart partilce_file selected_partilce_id_file outfile [-txt]\n");
        printf("To use a text selected particle file, use \"-txt\"");
        exit(1);
    }
    
    if(argn == 5){
        if(strcmp(argv[4], "-txt") != 0){
            printf("Usage: selectpart partilce_file selected_partilce_id_file outfile [-txt]\n");
            printf("To use a text selected particle file, use \"-txt\"");
            exit(1);
        }else{
            istext = true;
        }
    }
    
    string partfile = argv[1];
    string seletfile = argv[2];
    string outfile = argv[3];
    
    
    Pdm dp;
    tipsy_header header;
    XDR xdr;
    FILE *fp;
    FILE *ofp;
    
    int * partlist;
    Pdm * all_dm_particles;
    int partnum = 0;
    int partind = 0;
    
    if(!istext){
        fstream lstr(seletfile.c_str(), ios::binary | ios::in);
        if(lstr.good()){
            lstr.read((char *) &partnum, sizeof(int));
            printf("Will select %d partilces from simulation results.\n", partnum);
            partlist = new int[partnum];
            all_dm_particles = new Pdm[partnum];
            lstr.read((char * ) partlist, sizeof(int) * partnum);
        }else{
            printf("Input list file not good!\n");
            exit(1);
        }
        lstr.close();
    }else{
        string line;
        ifstream myfile0(seletfile.c_str());
        while (std::getline(myfile0, line)){
            ++partnum;
        }
        myfile0.close();
        
        //test
        printf("Particle numbers: %d\n", partnum);
        
        ifstream myfile1(seletfile.c_str());
        partlist = new int[partnum];
        for(int i = 0; i < partnum; i++){
            myfile1 >> partlist[i];
        }
        myfile1.close();

    }
    
    //sort the particle list
    printf("Sorting IDs ...\n");
    std::sort (partlist, partlist + partnum);
    
    read_tipsy_header(partfile.c_str(), &header);
    printf("Number of particles = %d\n",header.nbodies);
    printf("Number of DM particles = %d\n",header.ndark);
    
    fp = fopen(partfile.c_str(),"r");
    xdrstdio_create(&xdr,fp,XDR_DECODE);
    int status = xdr_header(&xdr,&header);
    
    
    ofp = fopen(outfile.c_str(), "w");
    
    printf("Reading particles...\n");
    for(int i=0;i<header.ndark && partind < partnum; i++) {
        status = xdr_dark(&xdr,&dp);
        
        if (status != 1) {
            fprintf(stderr,"Error reading dark particle from input file.\n");
            exit(1);
        }
        
        
        if(i == partlist[partind]){
            all_dm_particles[partind] = dp;
            partind ++;
        }
    }
    
    header.nbodies = partnum;
    header.ndark = partnum;
    header.nsph = 0;
    header.nstar = 0;
    printf("Write to output...\n");
    write_tipsy_file(outfile.c_str(), header, all_dm_particles, partnum, 0.0);
    printf("Finished.\n");
    delete partlist;
    delete all_dm_particles;
}
