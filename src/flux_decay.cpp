#include <string>
#include <cstring>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <cmath>

#include "flux.h"
//decay
//With units: GeV.kpc/cm^6
MAPTYPE getflux(Parameters * par_, DMParticle & current_part, MAPTYPE distances){
    MAPTYPE decay_flux_to_cgs = (par_->codeunits.mass_to_cgs / pow(par_->codeunits.length_to_cgs, 2));
    MAPTYPE unit_factor = pow((par_ -> natconst.c_in_cgs), 2) /
                              (par_->units.eV_in_cgs * 1.0e9) / (par_->units.Mpc_in_cgs * 1.0e-3)
                            * decay_flux_to_cgs;
    MAPTYPE fluxes = unit_factor * current_part.mass / (4.0 * PI * distances * distances);
    return fluxes;
}
