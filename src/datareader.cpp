#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "tipsydefs.h"
#include "readfiles.h"
#include "datareader.h"

#define OUTPUT_RESULT_

DataReader::DataReader(string basedir, string basename){
    particle_file_name = basedir + "/" + basename;
    density_file_name = particle_file_name + ".den32";
    hsmooth_file_name = particle_file_name + ".hsm32";
    particles = NULL;
    numParts_ = 0;
}

DataReader::~DataReader(){
    if(particles != NULL){
        delete particles;
    }
}

bool DataReader::readParticle(){
    int loopi;
   
    tipsy_header tips;
    XDRFileReader::read_tipsyheader(particle_file_name, &tips);
    numParts_ = tips.ndark;
    Pdm * dmpart;
    printf("Reading %d particles...\n", numParts_);
    
    XDRFileReader::read_particles(particle_file_name, dmpart);
    //copy the data to the DMparticles array
    printf("Particles read.\n");
    particles = new DMParticle[numParts_];
    for(loopi = 0; loopi < numParts_; loopi++){
        particles[loopi].mass = dmpart[loopi].mass;
        particles[loopi].posx = dmpart[loopi].pos[0];
        particles[loopi].posy = dmpart[loopi].pos[1];
        particles[loopi].posz = dmpart[loopi].pos[2];
        //printf("%d\n", loopi);
        particles[loopi].velx = dmpart[loopi].vel[0];
        particles[loopi].vely = dmpart[loopi].vel[1];
        particles[loopi].velz = dmpart[loopi].vel[2];
        particles[loopi].eps = dmpart[loopi].eps;
        particles[loopi].phi = dmpart[loopi].phi;
    }
    //printf("fdsfjals\n");
    delete dmpart;
    
    printf("Reading density...\n");
    //reader scalars
    float * density;
    int testn = 0;
    XDRFileReader::read_scalar(density_file_name, density, testn);
    if(testn != numParts_){
        printf("Number of particles not consistent!\n");
        exit(1);
    }
    for(loopi = 0; loopi < numParts_; loopi++){
        particles[loopi].dens = density[loopi];
    }
    delete density;
    
    //reader scalars
    printf("Reading hsmooth...");
    float * hsmooth;
    //int testn = 0;
    XDRFileReader::read_scalar(density_file_name, hsmooth, testn);
    if(testn != numParts_){
        printf("Number of particles not consistent!\n");
        exit(1);
    }
    for(loopi = 0; loopi < numParts_; loopi++){
        particles[loopi].hsmooth = hsmooth[loopi];
    }
    delete hsmooth;
    
#ifdef OUTPUT_RESULT_
    for(loopi = 0; loopi < numParts_; loopi++){
        printf("Particle %d:{%e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e}\n",
               loopi,
               particles[loopi].mass,
               particles[loopi].dens,
               particles[loopi].hsmooth,
               particles[loopi].posx,
               particles[loopi].posy,
               particles[loopi].posz,
               particles[loopi].velx,
               particles[loopi].vely,
               particles[loopi].velz,
               particles[loopi].eps,
               particles[loopi].phi);
    }
#endif
    return true;
}

DMParticle * DataReader::getParticles(){
    return particles;
}

int DataReader::getNumParts(){
    return numParts_;
}
