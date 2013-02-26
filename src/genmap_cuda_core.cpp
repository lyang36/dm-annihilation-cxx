//this file copies all the code from genmap_cuda but add applying core function
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
#include "kernel.h"
#include "flux.h"
#include "halocore/cudahalocore.h"

//this file implements the mapgenerator's gen function
//Step 1: Divide the high resolution Healpix grid to lower resolution one, calculate how many particles are in each pixel
//Step 2: Calculate the tetra list for each pixel
//Step 3: Calculate the normalization for each particle
//Step 4: Calculate the final value for each high resolution pixel
void MapGenerator::start(){
    MAPTYPE distances;
    MAPTYPE fluxes;
    MAPTYPE angular_radius;
    MAPTYPE costheta;
	MAPTYPE theta;
	MAPTYPE phi;
    
    MAPTYPE unit_factor = pow(pow((par_ -> natconst.c_in_cgs), 2) /
                              (par_->units.eV_in_cgs * 1.0e9), 2) / (par_->units.Mpc_in_cgs * 1.0e-3)
                            * par_->codeunits.annihilation_flux_to_cgs;
    
	int Nside = par_->map.Nside;
	int Npix = par_->map.Npix;
    MAPTYPE dOmega = par_->map.dOmega;
    MAPTYPE theta0 = acos( 1.0 - dOmega/(2.0*3.141592) );
   
    
    //int NpixCoase_;
    //NpixCoase_ = Npix / BLOCK_N_DIVIDER / BLOCK_N_DIVIDER;
    //the actuall base
    Healpix_Base base(Nside, NEST, SET_NSIDE);
    Healpix_Base superBase(Nside / BLOCK_N_DIVIDER, RING, SET_NSIDE);
    
    map_ = (MAPTYPE *) calloc(Npix, sizeof(MAPTYPE));
    int Np = reader_->getPartNum();
    
    //int rec = Np / 50;
    //int ll = 0;
    cout << "Creating map!!!" << endl;
	//cout << "---10---20---30---40---50---60---70---80---90--100%\n";
    

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
    

    cudaError_t status;
    //initialize
    status = initializeCUDA(Nside, par_->memParts);
    vector<int> * pixellist = new vector<int>[par_->memParts];
    setCudaPixelList(pixellist);
    
    if(status != cudaSuccess){
    	printf("CUDA initialize error!\n");
    	exit(1);
    }
    
    status = initializeCore(getDevParticle(), par_->memParts);
	if(status != cudaSuccess){
    	printf("HaloCore initialize error!\n");
    	exit(1);
    }
	MAPTYPE * coreCorr = new MAPTYPE[par_ -> memParts];
    int count = 0;
	
    
    while(reader_->hasNext()){
    	parts = reader_->getBuf();
    	count += reader_->getMemparts();
        printf("Particles: %d\n", count);
    	//first step -- filter the particles
		applyingCore(parts, reader_->getMemparts(), coreCorr);
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
			fluxes = getflux(par_, current_part, distances) * coreCorr[i];
			//printf("%f\n", coreCorr[i]);

			calc_angles(current_part.posx-oposx, current_part.posy-oposy,
						current_part.posz-oposz, distances, par_,
						costheta, phi);
			
			theta = acos(costheta);
			angular_radius = current_part.hsmooth / distances;
			
			
			pointing p(theta,phi);
			vec3 vec;
			ang2vec(theta,phi,&vec);
			
			//put the calculated value into particle data
			current_part.mass = fluxes;
			current_part.dens = distances;
			current_part.hsmooth = theta;
			current_part.posx = phi;
			current_part.posy = angular_radius;
			current_part.posz = cos(2.0 * angular_radius);
            current_part.velx = vec.x;
			current_part.vely = vec.y;
			current_part.velz = vec.z;
			
			if( 2.0*angular_radius < theta0 ) {
				int pix = base.ang2pix(p);
				map_[pix] += fluxes;
				current_part.eps = -1;
			}
			
			//query the pixles in the superBase
			pixellist[i].clear();
			superBase.query_disc_inclusive(p, 2.0 * angular_radius, pixellist[i]);
		}
		
		//pass particles, map to CUDA, calculating
		status = calulateMapWithCUDA(map_, parts, reader_->getMemparts());
	    if(status != cudaSuccess){
	    	printf("CUDA initialize error!\n");
	    	exit(1);
	    }
		
		reader_->loadBuffer();
    }
    isFinished_ = true;
    
    //cleaing up CUDA
    cudaCleaingUp();
	clearUpHaloCore();
    
    cout << "\nFinished!." << endl;
    
    MAPTYPE * ringmap = new MAPTYPE[Npix];

    for(int i = 0; i < Npix; i++){
        map_[i] /= par_->map.dOmega;
        ringmap[base.nest2ring(i)] = map_[i];
    }
    for(int i = 0; i < Npix; i++){
        map_[i] = ringmap[i];
    }

    delete ringmap;
    printf("%e\n%e\n%e", map_[0], map_[100], map_[10000]);
}
