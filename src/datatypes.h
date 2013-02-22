#ifndef __DATATYPE__LY
#define __DATATYPE__LY
//the particle for visualization
#include <string>
#include "tipsydefs.h"
using namespace std;

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

const double PI = M_PI;

//define the map's type
typedef double MAPTYPE;

//float can lead to incorrect result
//typedef float MAPTYPE;


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
        return (*this);
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
        //return *this;
    };
    
    DMParticle(){};
};



class Info{
public:
	double Omega_M;
	double h100;
	double n_s;
	double sigma_8;
	double Lbox_in_Mpc;
	double particle_masses[10];
	long particle_numbers[10];
	double centerpos[3];
	double centervel[3];
	double Mvir_in_Msun;
	double M200_in_Msun;
	double M200crit_in_Msun;
	double Rvir_in_Mpc;
	double R200_in_Mpc;
	double R200crit_in_Mpc;
	double Vmax_in_kms;
	double RVmax_in_kpc;
	double rconverg_in_kpc;
	double params_NFW[3];
	double params_GNFW[3];
	double params_Einasto[3];
};


class Codeunits{
public:
	double mass_to_cgs;
	double mass_to_Msun;
	double length_to_Mpc;
	double length_to_cgs;
	double time_to_cgs;
	double velocity_to_cgs;
	double density_to_cgs;
	double annihilation_flux_to_cgs;
};


class Units{
public:
	static const double Msun_in_cgs = 1.9889212e+33;
	static const double Mpc_in_cgs = 3.08568025e+24;
	static const double eV_in_cgs = 1.602e-12;
};

class Natconst{
public:    
	double h100;
	static const double G_in_cgs = 6.67259e-8;
	double rho_crit_in_cgs;
	double Omega_m;
	double Delta_vir;
    double Rvir_MW_in_Mpc;// = 0.275;
	static const double c_in_cgs = 2.99792458e10;
};


class Halo{
public:
	double Mvir_in_Msun;
	double M200_in_Msun;
	double M200crit_in_Msun;
	double Rvir_in_Mpc;
	double R200_in_Mpc;
	double R200crit_in_Mpc;
	double Vmax_in_kms;
	double RVmax_in_kpc;
	double rconverg_in_kpc;
	double params_NFW[3];
	double params_GNFW[3];
	double params_Einasto[3];
	double shape_r;
	double shape;
};

class Dm{
public:
	double M_in_GeV;// = 46.0;
	double M_in_cgs;// = 0.0;
	double sigma_v_in_cgs;// = 5.0e-26;
};

class Params{
public:
	double z;
	double cpos[3];
	double cvel[3];
	double opos[3];
	double otheta;
	double ophi;
	double Lbox_in_Mpc;
	int    particle_numbers[10];
	double particle_masses[10];
	double particle_masses_in_Msun[10];
};


class Map{
public:
	string projection;
	int Nside;
	int Npix;
	double dOmega;
	double theta0;
    void setNside(int nside){
        Nside = nside;
        Npix = 12 * Nside * Nside;
        dOmega = 4 * PI / Npix;
        theta0 = 0.0;
    };
};

#endif

