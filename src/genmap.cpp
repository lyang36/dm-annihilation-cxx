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
#include "anglefuncs.h"

//this file implements the mapgenerator's gen function
void MapGenerator::start(){
    MAPTYPE distances;
    MAPTYPE fluxes;
    MAPTYPE angular_radius;
    MAPTYPE costheta;
	MAPTYPE theta;
	MAPTYPE phi;
    MAPTYPE * weight;
    
    MAPTYPE unit_factor = pow(pow((par_ -> natconst.c_in_cgs), 2) /
                              (par_->units.eV_in_cgs * 1.0e9), 2) / (par_->units.Mpc_in_cgs * 1.0e-3)
                            * par_->codeunits.annihilation_flux_to_cgs;
    
	int Nside = par_->map.Nside;
	int Npix = par_->map.Npix;
    MAPTYPE dOmega = par_->map.dOmega;
    MAPTYPE theta0 = acos( 1.0 - dOmega/(2.0*3.141592) );
    
    Healpix_Base base(Nside, RING, SET_NSIDE);

    map_ = (MAPTYPE *) calloc(Npix, sizeof(MAPTYPE));
    weight = (MAPTYPE *) calloc(Npix, sizeof(MAPTYPE));
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
    
    DMParticle * parts;//current_part;
    MAPTYPE oposx = par_->params.opos[0];
    MAPTYPE oposy = par_->params.opos[1];
    MAPTYPE oposz = par_->params.opos[2];
    
    int count = 0;
    while(reader_->hasNext()){
    	parts = reader_->getBuf();
    	
		for( int i = 0; i< reader_->getMemparts(); i++){
			if(count % rec == 0){
				cout << "#";
				cout.flush();
			}
			count ++;
			//reader_->readParticle(& current_part);
			DMParticle &current_part = parts[i];
			
			//ignore the low resolution mass
			if(current_part.mass >= hires_particle_mass * 1.1){
				continue;
			}
			
			distances = sqrt((current_part.posx-oposx) * (current_part.posx-oposx)
							+(current_part.posy-oposy) * (current_part.posy-oposy)
							+(current_part.posz-oposz) * (current_part.posz-oposz));
			
			fluxes = unit_factor * current_part.dens * current_part.mass / (4.0 * PI * distances * distances);
			
			calc_angles(current_part.posx-oposx, current_part.posy-oposy,
						current_part.posz-oposz, distances, par_,
						costheta, phi);
			
			theta = acos(costheta);
			angular_radius = current_part.hsmooth / distances;
			
			
			pointing p(theta,phi);
			vec3 vec;
			ang2vec(theta,phi,&vec);
			
			int pix = base.ang2pix(p);
			if( 2.0*angular_radius < theta0 ) {
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
			
			/*static int mm = 0;
	
			if(mm < 100) {printf("%e %e %e %e %e %e\n", current_part.mass, current_part.dens,
						   current_part.hsmooth, current_part.posx,
						   current_part.posy, current_part.posz);
						   mm ++;}
	
			if(mm < 100) printf("angular_radius: %e, disc size: %d, dist: %f\n", 2.0*angular_radius, 
							npix_disc, distances);*/
			// get here only if the particle covers more than one pixel
			
			MAPTYPE weight_norm = 0.0;
			for(int j=0; j<npix_disc; j++) {
				int this_pix = pix_list[j];
				vec3 this_vec = base.pix2vec(this_pix);
				
				MAPTYPE d2 = acos( dotprod(this_vec,vec) ) / angular_radius;
				d2 = d2*d2;
				weight[j] = exp(-0.5 * d2 / 0.333);
				weight_norm += weight[j];		    
			}
			
			// apply weighted flux to map
			for(int j=0;j<npix_disc;j++) {
				// first normalize weight array
				weight[j] = weight[j] / weight_norm;
				map_[pix_list[j]] += weight[j] * fluxes;
				
			}  // loop over pixels covered by this particle
			/*if(mm < 100) printf("pix %d flux %e Map %e weightnorm %e\n", pix, fluxes, map_[pix], weight_norm);*/
		}
		reader_->loadBuffer();
    }
    isFinished_ = true;
    cout << "\nFinished!." << endl;
    
    for(int i = 0; i < Npix; i++){
        map_[i] /= par_->map.dOmega;
    }
    
    //printf("%e\n%e\n%e", map_[0], map_[100], map_[10000]);


    free(weight);
}
