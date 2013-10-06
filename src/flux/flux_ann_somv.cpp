#include <string>
#include <cstring>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <cmath>

#include "../flux.h"

/**
 *The first somerfield enhancement, <sigma v> \propto 1 / v
 */


#define SAT_V 1  // 1km/s

MAPTYPE getflux(Parameters * par_, DMParticle & current_part, MAPTYPE distances){
       MAPTYPE unit_factor = pow(pow((par_ -> natconst.c_in_cgs), 2) /
                              (par_->units.eV_in_cgs * 1.0e9), 2) / (par_->units.Mpc_in_cgs * 1.0e-3)
                            * par_->codeunits.annihilation_flux_to_cgs;
       MAPTYPE fluxes;
       MAPTYPE v = current_part.sigmav * par_->codeunits.velocity_to_cgs;
    
       if(v < SAT_V / 100000.0)
           v = SAT_V / 100000.0;

       MAPTYPE somer_v =  par_->natconst.c_in_cgs / v;

       fluxes = unit_factor * somer_v * current_part.dens 
               * current_part.mass / (4.0 * PI * distances * distances);
       return fluxes;
}
