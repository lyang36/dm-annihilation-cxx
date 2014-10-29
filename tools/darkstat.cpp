/*****************************************************
 * Show the statistics of the dark matter part in the
 * simulation of a tipsyfile.
 *
 * Author: Lin Yang
 * Date: Oct 2014
 *****************************************************/



#include <string>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <algorithm>    // std::sort
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "../src/tipsydefs.h"


#define MAX_NUM_MASS_TYPES 10

using namespace std;

double masses[MAX_NUM_MASS_TYPES];             // this vector records the mass of particles
int num_of_particles[MAX_NUM_MASS_TYPES];      // this vector records the number of each mass
int num_of_mass_types;

void initMassRecorder(){
    for(int i = 0; i < MAX_NUM_MASS_TYPES; i++){
        masses[i] = 0.0;
        num_of_particles[i] = 0;
    }
    num_of_mass_types = 0;
}

void addParticle(double mass){
    bool find_mass = false;
    for(int i = 0; i < num_of_mass_types; i++){
        if(masses[i] == mass){
            num_of_particles[i] ++;
            find_mass = true;
        }
    }
    if(!find_mass){
        num_of_mass_types ++;
        if(num_of_mass_types > MAX_NUM_MASS_TYPES){
            num_of_mass_types = MAX_NUM_MASS_TYPES;
        }
        masses[num_of_mass_types - 1] = mass;
        num_of_particles[num_of_mass_types - 1] ++;
    }
}


int main(int argn, char ** argv){
    if(argn != 2){
        printf("Usage: %s tipsyfile\n", argv[0]);
        exit(1);
    }
    

    string partfile = argv[1]; 
    
    Pdm dp;
    tipsy_header header;
    XDR xdr;
    
    FILE *fp;
    read_tipsy_header(partfile.c_str(), &header);
    printf("Number of particles =         %d\n",header.nbodies);
    printf("Number of DM particles =      %d\n",header.ndark);
    
    fp = fopen(partfile.c_str(),"r");
    
    if( NULL == fp ){
        fprintf(stderr, "DATAFILE ERROR!\n");
        exit(1);
    }
    
    xdrstdio_create(&xdr,fp,XDR_DECODE);
    int status = xdr_header(&xdr,&header);
    
    if(header.nsph != 0){
        fseek(fp, sizeof(double) + 6*sizeof(float) 
                + header.nsph * 12 * sizeof(float),
                0);
    }

    for(int i=0; (i<header.ndark); i++) {
        status = xdr_dark(&xdr,&dp);
        
        if (status != 1) {
            fprintf(stderr,"Error reading dark particle from input file.\n");
            exit(1);
        }
        addParticle(dp.mass);
    }
    
    printf("# of masses: %d", num_of_mass_types);
    for(int i = 0; i < num_of_mass_types; i++){
        printf("Mass_%d = %10.10e, #parts = %d\n", masses[i], num_of_particles[i]);
    }
    printf("\n");
    
    fclose(fp);
    
    return(0);
}
