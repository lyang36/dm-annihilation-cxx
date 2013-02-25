#include <string>
#include <cstring>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <cmath>

#include "flux.h"

#define SAT_V 0.000991235   // 1km/hr

MAPTYPE getflux(Parameters * par_, DMParticle & current_part, MAPTYPE distances){
       MAPTYPE unit_factor = pow(pow((par_ -> natconst.c_in_cgs), 2) /
                              (par_->units.eV_in_cgs * 1.0e9), 2) / (par_->units.Mpc_in_cgs * 1.0e-3)
                            * par_->codeunits.annihilation_flux_to_cgs;
       MAPTYPE fluxes;
       MAPTYPE v = current_part.sigmav;//sqrt(current_part.velx * current_part.velx +
                   //     current_part.vely * current_part.vely +
                   //     current_part.velz * current_part.velz);
       if(v < SAT_V) v = SAT_V; 

       MAPTYPE sophomer_v =  par_->natconst.c_in_cgs / (v * par_->codeunits.velocity_to_cgs);

       fluxes = unit_factor * sophomer_v * current_part.dens 
               * current_part.mass / (4.0 * PI * distances * distances);
       return fluxes;
}
