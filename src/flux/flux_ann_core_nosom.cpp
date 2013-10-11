#include <string>
#include <cstring>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <cmath>

#include "../flux.h"
#include "../halocore/halocore.h"

MAPTYPE getflux(Parameters * par_, DMParticle & current_part, MAPTYPE distances){
       static HaloCore halocore;
        
       MAPTYPE unit_factor = pow(pow((par_ -> natconst.c_in_cgs), 2) /
                              (par_->units.eV_in_cgs * 1.0e9), 2) / (par_->units.Mpc_in_cgs * 1.0e-3)
                            * par_->codeunits.annihilation_flux_to_cgs;
       MAPTYPE r = halocore.getCorrection(current_part.posx, current_part.posy, current_part.posz);      
       MAPTYPE fluxes = r * r * unit_factor * current_part.dens * current_part.mass / (4.0 * PI * distances * distances);
       return fluxes;
}
