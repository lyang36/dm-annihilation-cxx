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
#include "flux.h"

#include "opengl/render.h"

//this file implements the mapgenerator's gen function
void MapGenerator::start(){
    MAPTYPE distances;
    MAPTYPE fluxes;
    MAPTYPE angular_radius;
    MAPTYPE costheta;
	MAPTYPE theta;
	MAPTYPE phi;
    RenderParticle * renderParts = NULL;
    
    //MAPTYPE unit_factor = pow(pow((par_ -> natconst.c_in_cgs), 2) /
    //                          (par_->units.eV_in_cgs * 1.0e9), 2) / (par_->units.Mpc_in_cgs * 1.0e-3)
    //                        * par_->codeunits.annihilation_flux_to_cgs;
    
	int Nside = par_->map.Nside;
	int Npix = par_->map.Npix;
    MAPTYPE dOmega = par_->map.dOmega;
    MAPTYPE theta0 = acos( 1.0 - dOmega/(2.0*3.141592) );
    map_ = new double[Npix];

    Healpix_Base base(Nside, RING, SET_NSIDE);

    int Np = reader_->getPartNum();
    
    int rec = Np / 50;
    int ll = 0;

    

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
    
    
    render oglRender(*par_, Nside, 256, Nside);
    
    int count = 0;
    
    cout << "Creating map!!!" << endl;
	cout << "---10---20---30---40---50---60---70---80---90--100%\n";
    while(reader_->hasNext()){
    	parts = reader_->getBuf();
    	
        renderParts = (RenderParticle *) realloc(renderParts, sizeof(RenderParticle) * reader_->getMemparts());
		for( int i = 0; i< reader_->getMemparts(); i++){
			if(count % rec == 0){
				cout << "#";
				cout.flush();
			}
			count ++;
			DMParticle &current_part = parts[i];
			
			//ignore the low resolution mass
			if(current_part.mass >= hires_particle_mass * 1.1 
             ||current_part.dens < 0){
				continue;
			}
			
			distances = sqrt((current_part.posx-oposx) * (current_part.posx-oposx)
							+(current_part.posy-oposy) * (current_part.posy-oposy)
							+(current_part.posz-oposz) * (current_part.posz-oposz));
			
			//fluxes = unit_factor * current_part.dens * current_part.mass / (4.0 * PI * distances * distances);
			fluxes = getflux(par_, current_part, distances);
            
            //test
            fluxes /= (distances * distances);
            //printf("%e\n", fluxes);
            
            renderParts[i].x = current_part.posx;
            renderParts[i].y = current_part.posy;
            renderParts[i].z = current_part.posz;
            
            //printf("First: %f %f %f\n", renderParts[i].x, renderParts[i].y, renderParts[i].z);
            
            renderParts[i].densityfac1 = fluxes;
            renderParts[i].densityfac2 = 1.0;
            renderParts[i].hsmooth = current_part.hsmooth;
			
		}
        //printf("ok\n");
        oglRender.rend(renderParts, reader_->getMemparts());
		reader_->loadBuffer();
    }
    isFinished_ = true;
    cout << "\nFinished!." << endl;
    double * hmp = oglRender.getHealPixMap();
    
    //test
    double totalC = 0;
    
    for(int i = 0; i < Npix; i++){
        totalC += hmp[i];
        
        map_[i] = hmp[i] / par_->map.dOmega;;
    }
    
    //test
    printf("Total Conts: %f\n", totalC);
    
    
    free(renderParts);
}
