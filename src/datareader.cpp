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
   
    Tipsyheader tips;
    XDRFileReader::read_tipsyheader(particle_file_name, &tips);
    numParts_ = tips.ndark;
    Pdm * dmpart;
    XDRFileReader::read_particles(particle_file_name, dmpart);
    
    //copy the data to the DMparticles array
    particles = new DMParticle[numParts_];
    for(loopi = 0; loopi < numParts_; loopi++){
        particles[i].mass = dmpart[i].mass;
        particles[i].posx = dmpart[i].pos[0];
        particles[i].posy = dmpart[i].pos[1];
        particles[i].posz = dmpart[i].pos[2];
        particles[i].velx = dmpart[i].vel[0];
        particles[i].vely = dmpart[i].vel[1];
        particles[i].velz = dmpart[i].vel[2];
        particles[i].eps = dmpart[i].eps;
        particles[i].phi = dmpart[i].phi;
    }
    delete dmpart;
    
    //reader scalars
    float * density;
    int testn = 0;
    XDRFileReader::read_scalar(density_file_name, density, testn);
    if(testn != numParts_){
        printf("Number of particles not consistent!\n");
        exit(1);
    }
    for(loopi = 0; loopi < numParts_; loopi++){
        particles[i].dens = density[i];
    }
    delete density;
    
    //reader scalars
    float * hsmooth;
    int testn = 0;
    XDRFileReader::read_scalar(density_file_name, hsmooth, testn);
    if(testn != numParts_){
        printf("Number of particles not consistent!\n");
        exit(1);
    }
    for(loopi = 0; loopi < numParts_; loopi++){
        particles[i].hsmooth = hsmooth[i];
    }
    delete hsmooth;
    
#ifdef OUTPUT_RESULT_
    for(loopi = 0; loopi < numParts_; loopi++){
        printf("Particle %d:{%e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e}\n",
               loopi,
               particles[i].mass,
               particles[i].dens,
               particles[i].hsmooth,
               particles[i].posx,
               particles[i].posy,
               particles[i].posz,
               particles[i].velx,
               particles[i].vely,
               particles[i].velz,
               particles[i].eps,
               particles[i].phi)
    }
#endif
    return true;
}

DMParticle * DataReader::getParticles(){
    return particles;
}

int DMParticle::getNumParts(){
    return numParts_;
}
