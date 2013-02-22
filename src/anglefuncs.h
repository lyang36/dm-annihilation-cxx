#ifndef __ANGLEFUNC__
#define __ANGLEFUNC__
//input x, y, z
//output costheta, phi
#include "mapgenerator.h"
void calc_angles( MAPTYPE xpos, MAPTYPE ypos, MAPTYPE zpos, MAPTYPE &distances,
                 Parameters * params, MAPTYPE & costheta, MAPTYPE &phi);

void ang2vec(MAPTYPE theta, MAPTYPE phi, vec3 *vec);

#endif