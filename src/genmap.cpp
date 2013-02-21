#include <string>
#include <cstring>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <cmath>

//fits and HEALPIX
#include <healpix_base.h>
#include <healpix_map.h>
#include <arr.h>
#include <fitsio.h>
#include <fitshandle.h>
#include <healpix_map_fitsio.h>

#include "datatypes.h"
#include "mapgenerator.h"

//this file implements the mapgenerator's gen function


//input x, y, z
//output costheta, phi
void calc_angles( double xpos, double ypos, double zpos, double &distances,
                  Parameters * params, double & costheta, double &phi){
	double vec[] = { xpos, ypos, zpos};
	double temp[] = {0, 0, 0};

    for(int i = 0; i < 3; i ++){
        temp[i] = 0;
        for(int j = 0; j < 3; j ++){
            temp[i] += params->rotmatrix[i][j] * vec[j];
        }
    }
	
	double xcoord = vec[0];
	double ycoord = vec[1];
	double zcoord = vec[2];
    
	costheta = zcoord / distances;
    
	phi = atan2( ycoord, xcoord );
	
	if( phi < 0 ){
		phi += 2.0 * PI;
	}
    
	//a few adjustments...
    
	//one more rotation by 180 degrees...
	phi -= PI;
    
	//force phi to lie between 0 and 2*pi
	if( phi < 0 ){
		phi = 2.0 * PI + phi;
	}
    
	if( phi > 2 * PI){
		phi = phi - 2.0 * PI;
	}
}


void ang2vec(double theta, double phi, vec3 *vec) {
    
    double sz;
    double PI=M_PI;
    
    if( theta<0. || theta>PI) {
        fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
        exit(1);
    }
    
    sz = sin(theta);
    
    vec->x = sz * cos(phi);
    vec->y = sz * sin(phi);
    vec->z = cos(theta);
    
}


void MapGenerator::start(){
    double distances;
    double fluxes;
    double angular_radius;
    double costheta;
	double theta;
	double phi;
    double * weight;
    
	int Nside = par_->map.Nside;
	int Npix = par_->map.Npix;
    double dOmega = par_->map.dOmega;
    double theta0 = acos( 1.0 - dOmega/(2.0*3.141592) );
    
    Healpix_Base base(Nside, RING, SET_NSIDE);
    
    int count=0;
    map_ = (double *) calloc(Npix, sizeof(double));
    weight = (double *) calloc(Npix, sizeof(double));
    int Np = reader_->getPartNum();
    
    int rec = Np / 50;
    int ll = 0;
    cout << "Creating map!!!" << endl;
	cout << "---10---20---30---40---50---60---70---80---90--100%\n";
    

    //remove low resolution particles
    float hires_particle_mass = 1.0e9f;
    for(int i = 0; i < 10; i ++){
        if(par_->params.particle_masses[i] > 0 &&
           par_->params.particle_masses[i] < hires_particle_mass){
            hires_particle_mass = par_->params.particle_masses[i];
        }
    }
    //= min(master.params.particle_masses)
    
    DMParticle current_part;
    double oposx = par_->params.opos[0];
    double oposy = par_->params.opos[1];
    double oposz = par_->params.opos[2];
    
    for( int i = 0; i< Np; i++){
        if(i % rec == 0){
            cout << "#";
            cout.flush();
		}
        reader_->readParticle(& current_part);
        
        //printf("%e %e %e %e %e %e\n", current_part.mass, current_part.dens, 
        //                current_part.hsmooth, current_part.posx,
        //                current_part.posy, current_part.posz);
        //ignore the low resolution mass
        if(current_part.mass >= hires_particle_mass * 1.1){
            continue;
        }
        
        distances = sqrt((current_part.posx-oposx) * (current_part.posx-oposx)
                        +(current_part.posy-oposy) * (current_part.posy-oposy)
                        +(current_part.posz-oposz) * (current_part.posz-oposz));
        
        fluxes = par_->codeunits.annihilation_flux_to_cgs * current_part.dens * current_part.mass / (4.0 * PI * distances * distances);
        
        calc_angles(current_part.posx-oposx, current_part.posy-oposy,
                    current_part.posz-oposz, distances, (par_),
                    costheta, phi);
        
        theta = acos(costheta);
		angular_radius = current_part.hsmooth / distances;
        
        
        pointing p(theta,phi);
        vec3 vec;
        ang2vec(theta,phi,&vec);
        
        int pix = base.ang2pix(p);
        if( 2.0*angular_radius < theta0 ) {
            count++;
            map_[pix] += fluxes;
            continue;
        }
        
        vector<int> pix_list;
        base.query_disc(p, 2.0*angular_radius, pix_list);
        
        int npix_disc = pix_list.size();
        
        if(npix_disc < 2) {
            map_[pix] += fluxes;
            continue;
        }
        
        //printf("angular_radius: %e, disc size: %d, dist: %f\n", 2.0*angular_radius, 
        //                npix_disc, distances);
        // get here only if the particle covers more than one pixel
        double weight_norm = 0.0;
        for(int j=0; j<npix_disc; j++) {
            int this_pix = pix_list[j];
            vec3 this_vec = base.pix2vec(this_pix);
            
            double d2 = acos( dotprod(this_vec,vec) ) / angular_radius;
            d2 = d2*d2;
            weight[j] = exp(-0.5 * d2 / 0.333);
            weight_norm += weight[j];		    
        }
        
        // apply weighted flux to map
        for(int j=0;j<weight_norm;j++) {
            // first normalize weight array
            weight[j] /= weight_norm;
            map_[pix_list[j]] += weight[j] * fluxes;
            
        }  // loop over pixels covered by this particle
        
    }
    isFinished_ = true;
    cout << "finished!." << endl;
    
    double unit_factor = pow(pow((par_ -> natconst.c_in_cgs), 2) /
                             (par_->units.eV_in_cgs * 1.0e9), 2) / (par_->units.Mpc_in_cgs * 1.0e-3);
    for(int i = 0; i < Npix; i++){
        map_[i] *= unit_factor  / par_->map.dOmega;
        printf("%e\n", map_[i]);
    }
    
    free(weight);
}
