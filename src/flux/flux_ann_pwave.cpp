#include <string>
#include <cstring>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <cmath>

#include "../flux.h"

/**
 *The p-wave annihilatoin, <sigma v> \propto v^2
 */

MAPTYPE getflux(Parameters * par_, DMParticle & current_part, MAPTYPE distances){
       MAPTYPE unit_factor = pow(pow((par_ -> natconst.c_in_cgs), 2) /
                              (par_->units.eV_in_cgs * 1.0e9), 2) / (par_->units.Mpc_in_cgs * 1.0e-3)
                            * par_->codeunits.annihilation_flux_to_cgs;
       MAPTYPE fluxes;
       MAPTYPE v = current_part.sigmav;//sqrt(current_part.velx * current_part.velx +
                   //     current_part.vely * current_part.vely +
                   //     current_part.velz * current_part.velz);
       //if(v < SAT_V) v = SAT_V;

       MAPTYPE somer_v = (v * par_->codeunits.velocity_to_cgs) / par_->natconst.c_in_cgs;
    
        somer_v = somer_v * somer_v;

       fluxes = unit_factor * somer_v * current_part.dens
               * current_part.mass / (4.0 * PI * distances * distances);
       return fluxes;
}
