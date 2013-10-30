#include <string>
#include <cstring>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <cmath>

#include "../flux.h"


MAPTYPE getflux(Parameters * par_, DMParticle & current_part, MAPTYPE distances){
       MAPTYPE unit_factor = pow(pow((par_ -> natconst.c_in_cgs), 2) /
                              (par_->units.eV_in_cgs * 1.0e9), 2) / (par_->units.Mpc_in_cgs * 1.0e-3)
                            * par_->codeunits.annihilation_flux_to_cgs;
       
       //distance to the halo center
       double dx = current_part.posx - par_->cpos[0];
       double dy = current_part.posy - par_->cpos[1];
       double dz = current_part.posz - par_->cpos[2];    
       //convert distance to kpc
       double r = sqrt(dx * dx + dy * dy + dz * dz) * (par_->codeunits.length_to_Mpc) / 1000;
        
       double factor = (14.6534 * pow(1 + 0.0355872 * r, 1.76))/
               (pow(1 + 0.10888 * pow(r, 0.76), 3.34211) * pow(r, 0.13));
       
       
       MAPTYPE fluxes = unit_factor * current_part.dens * current_part.mass / (4.0 * PI * distances * distances);
       
       factor *= (factor * factor); 

       return fluxes;

}
