#ifndef __DATATYPE__LY
#define __DATATYPE__LY
//the particle for visualization
#include "tipsydefs.h"

class DMParticle {
public:
    float mass;
    float dens;
    float hsmooth;
    float posx;
    float posy;
    float posz;
    float velx;
    float vely;
    float velz;
    float eps;
    float phi;
    
    void setPdm(Pdm dp){
        mass = dp.mass;
        posx = dp.pos[0];
        posy = dp.pos[1];
        posz = dp.pos[2];
        velx = dp.vel[0];
        vely = dp.vel[1];
        velz = dp.vel[2];
        eps = dp.eps;
        phi = dp.phi;
    };
    
    DMParticle & operator=(const DMParticle &rhs){
        mass = rhs.mass;
        dens = rhs.dens;
        hsmooth = rhs.hsmooth;
        posx = rhs.posx;
        posy = rhs.posy;
        posz = rhs.posz;
        velx = rhs.velx;
        vely = rhs.vely;
        velz = rhs.velz;
        eps = rhs.eps;
        phi = rhs.phi;
    };
    DMParticle(const DMParticle &rhs){
        mass = rhs.mass;
        dens = rhs.dens;
        hsmooth = rhs.hsmooth;
        posx = rhs.posx;
        posy = rhs.posy;
        posz = rhs.posz;
        velx = rhs.velx;
        vely = rhs.vely;
        velz = rhs.velz;
        eps = rhs.eps;
        phi = rhs.phi;
    };
    
    DMParticle(){};
};

#endif

