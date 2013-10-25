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

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "datatypes.h"
#include "mapgenerator.h"
#include "anglefuncs.h"
#include "chealpixrender.h"
#include "flux.h"


void MapGenerator::start(){
    double distances;
    double fluxes;
    double angular_radius;
    double costheta;
	double theta;
	double phi;
    
    float unit_factor = pow(pow((par_ -> natconst.c_in_cgs), 2) /
                              (par_->units.eV_in_cgs * 1.0e9), 2) / (par_->units.Mpc_in_cgs * 1.0e-3)
                            * par_->codeunits.annihilation_flux_to_cgs;
    
	int Nside = par_->map.Nside;
	int Npix = par_->map.Npix;
    float dOmega = par_->map.dOmega;
    float theta0 = acos( 1.0 - dOmega/(2.0*3.141592) );
   
    
    //int NpixCoase_;
    //NpixCoase_ = Npix / BLOCK_N_DIVIDER / BLOCK_N_DIVIDER;
    //the actuall base
    Healpix_Base base(Nside, NEST, SET_NSIDE);
    
    map_ = (double *) calloc(Npix, sizeof(double));
    //use map1_ to get the value from device
    float * map1_ = (float *) calloc(Npix, sizeof(float));
    int Np = reader_->getPartNum();
    int rec = Np / par_->memParts / 50;
    if(rec < 1) rec = 1;

    //remove low resolution particles
    float hires_particle_mass = 1.0e9f;
    for(int i = 0; i < 10; i ++){
        if(par_->params.particle_masses[i] > 0 &&
           par_->params.particle_masses[i] < hires_particle_mass){
            hires_particle_mass = par_->params.particle_masses[i];
        }
    }
    
    DMParticle * parts;//current_part;
    renderpart * renderparts;
    renderparts = new renderpart[par_->memParts];
    
    float oposx = par_->params.opos[0];
    float oposy = par_->params.opos[1];
    float oposz = par_->params.opos[2];
    

    cudaError_t status;
    
    //printf("ok3\n");
    //initialize
    status = initializeCUDA(Nside, par_->memParts);
    //printf("ok4\n");

    if(status != cudaSuccess){
    	printf("CUDA initialize error!\n");
    	exit(1);
    }
    
    int count = 0;
    int rendercount = 0;
   

    cout << "Creating map!!!" << endl;
	cout << "---10---20---30---40---50---60---70---80---90--100%\n";
    
    
    int t_count = 0;
    while(reader_->hasNext()){
    	parts = reader_->getBuf();
    	count += reader_->getMemparts();
        if((count / par_->memParts) % rec == 0){
                cout << "#";
                cout.flush();
        }

		for( int i = 0; i< reader_->getMemparts(); i++){
			
            /////Particle calculating
			DMParticle &current_part = parts[i];
			
			current_part.eps = 0;
			
			//ignore the low resolution mass
			if(current_part.mass >= hires_particle_mass * 1.1 ||
               current_part.dens < 0){
				current_part.eps = -1;
				continue;
			}
			
			distances = sqrt((current_part.posx-oposx) * (current_part.posx-oposx)
							+(current_part.posy-oposy) * (current_part.posy-oposy)
							+(current_part.posz-oposz) * (current_part.posz-oposz));
			
			//fluxes = unit_factor * current_part.dens * current_part.mass / (4.0 * PI * distances * distances);
			fluxes = getflux(par_, current_part, distances);
            

			calc_angles(current_part.posx-oposx, current_part.posy-oposy,
						current_part.posz-oposz, distances, par_,
						costheta, phi);
            
			theta = acos(costheta);
			angular_radius = current_part.hsmooth / distances;
			
			
			pointing p(theta,phi);
			vec3 vec;
			ang2vec(theta,phi,&vec);
			
            
            t_count ++;
            
			if( 2.0*angular_radius < theta0 ) {
				int pix = base.ang2pix(p);

                //test
				map_[pix] += fluxes;

				current_part.eps = -1;
			}else{
                renderparts[rendercount].theta = theta;
                renderparts[rendercount].phi = phi;
                renderparts[rendercount].x = vec.x;
                renderparts[rendercount].y = vec.y;
                renderparts[rendercount].z = vec.z;
                renderparts[rendercount].angular_radius = angular_radius;

                //test
                renderparts[rendercount].flux = fluxes * 1.0e10;


                renderparts[rendercount].setup(Nside);
                rendercount ++;
            }

		}
		
		//pass particles to CUDA, calculating
		//printf("Ok0\n");
        status = calculateMapByGPU(renderparts, rendercount);
	    if(status != cudaSuccess){
	    	printf("CUDA initialize error!\n");
	    	exit(1);
	    }
		rendercount = 0;
        //printf("Ok2\n");
		reader_->loadBuffer();
    }
    
    getCUDAMap(map1_);
    
    isFinished_ = true;

    cudaCleaingUp();
    
    cout << "\nFinished, " << t_count <<" particles rendered!." << endl;
    
    
    //printf("dOmega->%f\n", par_->map.dOmega);
    for(int i = 0; i < Npix; i++){
        map_[i] += map1_[i] / 1.0e10;
		
		//test
        map_[i] /= par_->map.dOmega;
    }

    delete map1_;
    delete renderparts;
}
