/***********************************************
 * This file defines the datatypes used in this 
 * code. They include DM particles, the information
 * of halos etc.
 * Some of them are criticle to produce maps, while
 * others are not.
 *
 * Author: Lin F. Yang
 * Date: Feb. 2014
 ***********************************************/

#ifndef __DATATYPE__LY
#define __DATATYPE__LY
#include <string>
#include "tipsydefs.h"
using namespace std;

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#ifndef MAX_NUM_MASS_RES 
#define MAX_NUM_MASS_RES 10            // max number of mass resolutions
#endif 

const double PI = M_PI;

/* Define the map's type               */
typedef double MAPTYPE;

/* Float can lead to incorrect result  */
/* because of the precision            */
//typedef float MAPTYPE;



/************DMParticle******************
 * The required information of a 
 * DM particle to use to produce
 * Maps. Mass, position, velocity
 * are information directly from
 * the simulation. Density, hsmooth,
 * sigmav(velocity dispersion) are
 * calculated using the smooth code
 ****************************************/
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
    float sigmav;                       //the sigma of the velocity dispersion
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
        sigmav = rhs.sigmav;
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
        sigmav = rhs.sigmav;
        eps = rhs.eps;
        phi = rhs.phi;
    };
    
    DMParticle(){};
};


/***********The infomation of a halo *********/
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



/**************************************************/
/* The information of the halo to be used to produ-
 * ce maps.
 **************************************************/
class Info{
public:
    /******Required information******/
	double Omega_M;                                // Omega_m
	double h100;                                   // H or unit 100 m/s/Mpc
	double n_s;                                    // spectura index of the initial power spectrum
	double sigma_8;                                // sigma_8
	double Lbox_in_Mpc;                            // Simulation boxsize in Mpc
    double mass_to_Msun;                            // Simulation mass unit in msun
	double particle_masses[MAX_NUM_MASS_RES];      // Particle masses
    long   particle_numbers[MAX_NUM_MASS_RES];
	double centerpos[3];                           // center postion to be observed

    
    /**********Optional information***************/
    /* Asign arbitration infor would not
     * affect the final result.
     *********************************************/
    double centervel[3];
    Halo haloInfo;
	/*double Mvir_in_Msun;
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
	double params_Einasto[3];*/
};


/********The simulation units*******/
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


/********Some constant units*******/
class Units{
public:
	static const double Msun_in_cgs;    // = 1.9889212e+33;
	static const double Mpc_in_cgs;     // = 3.08568025e+24;
	static const double eV_in_cgs;      // = 1.602e-12;
};



/********Some Nature units*******/
class Natconst{
public:    
	double h100;
	static const double G_in_cgs;       // = 6.67259e-8;
	double rho_crit_in_cgs;
	double Omega_m;
    
    
	double Delta_vir;                   // the viral overdensity at some z for
                                        // a spherical collapse model
    
    double Rvir_MW_in_Mpc;              // = 0.275;
	static const double c_in_cgs;       // = 2.99792458e10;
};




class Dm{
public:
	double M_in_GeV;                    // = 46.0;
	double M_in_cgs;                    // = 0.0;
	double sigma_v_in_cgs;              // = 5.0e-26;
};



class Params{
public:
	double z;                           // redshift
	double cpos[3];                     // center pos
	double cvel[3];                     // center velocity, NOT used
	double opos[3];
	double otheta;                      // theta and phi of the observer
	double ophi;
	double Lbox_in_Mpc;                 // Box Size
	int    particle_numbers[MAX_NUM_MASS_RES];        // #particles of particle masses
	double particle_masses[MAX_NUM_MASS_RES];         // Particle masses
	double particle_masses_in_Msun[MAX_NUM_MASS_RES];
};


/***********Parameters of a map*********/
class Map{
public:
    
	string projection;                  // what is the projection to be used
                                        // Currently, only mollweide is implemented
    
	int Nside;                          // Nside of the HEALPIX
	int Npix;                           // #of pixels in the map
    
	double dOmega;                      // anuglar resolution of the map
	double theta0;                      // the smallest theta in the map, = 0.0
    
    void setNside(int nside){
        Nside = nside;
        Npix = 12 * Nside * Nside;
        dOmega = 4 * PI / Npix;
        theta0 = 0.0;
    };
};

#endif

