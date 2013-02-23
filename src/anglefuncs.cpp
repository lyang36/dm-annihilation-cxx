#include <cmath>
#include <rpc/types.h>
#include <rpc/xdr.h>

//fits and HEALPIX
#include <healpix_base.h>
#include <healpix_map.h>

#include "mapgenerator.h"
#include "anglefuncs.h"

using namespace std;
//input x, y, z
//output costheta, phi
void calc_angles( MAPTYPE xpos, MAPTYPE ypos, MAPTYPE zpos, MAPTYPE &distances,
                 Parameters * params, MAPTYPE & costheta, MAPTYPE &phi){
	MAPTYPE vec[] = { xpos, ypos, zpos};
	MAPTYPE temp[] = {0, 0, 0};
    
    for(int i = 0; i < 3; i ++){
        temp[i] = 0;
        for(int j = 0; j < 3; j ++){
            temp[i] += params->rotmatrix[j][i] * vec[j];
        }
    }
	
	MAPTYPE xcoord = vec[0];
	MAPTYPE ycoord = vec[1];
	MAPTYPE zcoord = vec[2];
    
	costheta = zcoord / distances;
    
	phi = atan2( ycoord, xcoord );
	
	if( phi < 0 ){
		phi += 2.0 * PI;
	}
    
	//a few adjustments...
    
	//one more rotation by 180 degrees...
	phi -= PI;
    
	//force phi to lie between 0 and 2*pi
	if( phi < 0 ){
		phi = 2.0 * PI + phi;
	}
    
	if( phi > 2 * PI){
		phi = phi - 2.0 * PI;
	}
}


void ang2vec(MAPTYPE theta, MAPTYPE phi, vec3 *vec) {
    
    MAPTYPE sz;
    MAPTYPE PI=M_PI;
    
    if( theta<0. || theta>PI) {
        fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
        exit(1);
    }
    
    sz = sin(theta);
    
    vec->x = sz * cos(phi);
    vec->y = sz * sin(phi);
    vec->z = cos(theta);
    
}
