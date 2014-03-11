/*****************************************************************************/
//particle flux shader
//convert the particles to point sprite
/*input postion and parameter of exactly a particle*/
/*outout: xc, yc, r, dtheta, */
//Author Lin F. Yang
//01/18/2014
/*****************************************************************************/
#version 120

//Constants
#define PI 3.1415926535897932
#define e32 4.4816890703380648226           //e^(3/2)

#define USEACTURALNORM

/******************************Parameters*************************************/
// rotation matrix
uniform mat3 rotmatrix;

// the observer position
uniform vec3 opos;

// the geometric factor defined as:
// geometry factor {size(projected size of the equater), viewportsize, maxpointsize }
// controlling the point sprite stuctures, determined by the
// pointsprite features, supported by hardware
uniform vec3 geofac;

// whether use the norm map? true: 1 else:0
// uniform int usenormmap;
/*****************************************************************************/

// the radius of the particle circle and the coordianate
// this variable will be passed to the fragment shader.
// vec4(newsize, npvec.x, npvec.y, npvec.z);
varying vec4 particle;


/*****************************************************************************/
// calculate the profile of a particle
// using the position vector r1, and angular size dtheta
// the position vector is a point on the unit sphere
float profile(vec3 r1,float dtheta){
    // get the particle position on the sphere
    vec3 r0 = vec3(particle.gba);
    
    // calculate the costheta between the queried position and particle position
    float costheta = dot(r0, r1)/(length(r0)*length(r1));
    
    // clamp the costherta to be between -1 and 1
    costheta = clamp(costheta, -1.0, 1.0);
    
    // use tayler seriers to calculate the angle^2
    // acos may too much error
    float t2 = 2.0 * ( 1.0 - costheta) + 1.0/3.0*(costheta - 1.0)*(costheta - 1.0) - 4.0/45.0 * (costheta - 1.0) *(costheta - 1.0)*(costheta - 1.0);
    
    // alternative method
    //float t2 = acos(costheta);
    //t2 = t2*t2;
   
    // calculate the final value
    float d2 = (t2 / dtheta / dtheta);
    
    //if(d2 > 1.0){
    //    return 0.0;
    //}
    return exp(- 1.5 * d2);
}

// reverse stereoprojection, given a point on the tangential plane, output the point on the sphere
vec3 prev(vec2 xy){
    float r2 = xy.x*xy.x + xy.y*xy.y;
    return vec3(2.0 * xy.x/(1.0 + r2), 2.0 * xy.y/(1.0 + r2), (r2 - 1.0)/(r2 + 1.0));
}


//projected profile
float profPRJ(vec3 r1, float dtheta){
    //test
    //return 1.0;
    
    return (1.0 - r1.z) * (1.0 - r1.z) * profile(r1, dtheta);
}


// calculate the normalization of a particle by counting all the pixels
// this is very accurate, but very slow to get the result
float calc_norm(vec2 svec, float newsize, float dtheta){
    float norm = 0.0;
    
    vec2 coor = svec * geofac.y / 2.0;
    
    
    float x=0.0;
    float y=0.0;
    for(x = 0.0; x < newsize; x++){
        for(y = 0.0; y < newsize; y++){
            // the pixel coordinates is at the center of the pixel
            float px = (x+0.5)/newsize;
            float py = (y+0.5)/newsize;
            px = 2.0*(px-0.5); // -1...1
            py = 2.0*(py-0.5);
            vec2 xy = vec2(px, py);
            float u = dot(xy, xy);
            
            // reject the points that are outside the circle
            if (u > 1.0)
                continue;
            
            // calculate the actual coordinates on the projected plane
            vec2 xyp = xy * (newsize / 2.0) + coor;
            vec2 xyr = xyp / (geofac.y / 2.0);
            //float pr2 = dot(xyr, xyr);
            
            // calculate the flux
            norm += profPRJ(prev(xyr), dtheta);
            //4.0/(1.0+pr2)/(1.0+pr2) * profile(prev(xyr), dtheta);
        }
    }
    return 1.0/norm;
}

// use an analytical method to calculate the normalization, fast, but not very accurate
float calc_norm1(float theta0){
    return 2.0*(-1.0 + e32) * PI * theta0 * theta0 / (3.0 * e32) -
        (PI * (-5.0 + 2.0 * e32) * theta0 * theta0 * theta0 * theta0) / (27.0 * e32)
        +(PI * (-29.0 + 8.0 * e32) * PI * theta0 * theta0 * theta0 * theta0 * theta0 * theta0) / ( 1620.0 * e32);
}

void main(){
    
    // the new position to be calculated
    vec4 newpos;
   
    // the position vector of the target particle in the frame of the observer
    // input x, y z of the particle
    // transform it to the stereoprojection plane
    vec3 pvec = vec3(gl_Vertex.x, gl_Vertex.y, gl_Vertex.z) - opos;

    
    // the parameters of a partilce is stored in its color
    // mass, density and hsmooth
    vec4 parameter = vec4(gl_Color);
    
    // anular radias of a particle
    float dtheta;
    
    // calculated point size of the final object
    float dsize;
 
    // the distance from the observer to the particle
    float distance = length(pvec);
    
    // parameter.b is hsmooth, so this is the half angular size
    // used in the gaussian parameters
    dtheta = parameter.b / distance;    //2.186
    float angdsize = 2.0 * dtheta;
    
    // rotation and normalize the position vector
    vec3 npvec = normalize(rotmatrix * pvec);
    
    // rotate 180 degree to match the c++ rotation
    // npvec.x = - npvec.x;
    
    // caltucle the angle of the observed particle
    float costheta = clamp(npvec.z, -1.0, 1.0);//dot(npvec, nzaxis);
    
    // costheta -3.0 / 5.0;
    // now caluce theta from costheta
    float theta = acos(costheta);      //0.955
    
    
    
    if((theta > (PI / 2.0) || (theta + angdsize) >= PI / 2.0)
       && (angdsize < PI / 2.0))
    {   // if the particle is in the lower sphere, do the following
        // also ignore those particles that are too large
        
        // calculate sintheta
        float sintheta = sin(theta);
        
        float sinphi;
        float cosphi;
        
        
        if(sintheta < 1.0e-8 ){
            sinphi = 0.0;
            cosphi = 1.0;
        }else{
            sinphi = npvec.y/sintheta;//newpvec.y / sintheta;
            cosphi = npvec.x/sintheta;//newpvec.x / sintheta;
        }
        
        // calculate the flux, encoded by parameter.g*parameter.r
        float flux = parameter.r;// / (4.0 * PI * distance * distance);
        
        // the output parameters of a point sprite
        float xc, yc, r;
    
        // transform the vertex:
        // stereoproject a circle from a sphere to the plane
        float sintpr = sin(theta + angdsize);
        float costpr = cos(theta + angdsize);
        float sintmr = sin(theta - angdsize);
        float costmr = cos(theta - angdsize);
        float a = sintpr/(1.0-costpr);
        float b = sintmr/(1.0-costmr);
        r = -(a - b)/2.0;
        float prho = (a + b)/2.0;
        xc = prho * cosphi;
        yc = prho * sinphi;
        
        // calculate the point size
        //float newsize = floor(r *geofac.y); ///!!!!!!!!!!!!!!!!
        float newsize = ceil(r *geofac.y);

        // calculate the actuall point position on the screen
        newpos = vec4(xc * geofac.x, yc * geofac.x, 0.0, 1.0);
        
        // if the point size been calculated is too large than the largest allowed,
        // clamp it to the limited value
        if(newsize > geofac.z){
            dsize = geofac.z / newsize * r;
            newsize = geofac.z;
        }else{
            dsize = r;
        }
        
        // if the point size is less than 1, asign it to be 1
        if(newsize < 1.0){
            newsize = 1.0;
        }
        
        // output the pointsize to the fragment shader
        gl_PointSize = newsize;  //point size
    


        // particle must be written before fhe nomal fac
        // particle = vec4(dsize, npvec.x, npvec.y, npvec.z);
        particle = vec4(newsize, npvec.x, npvec.y, npvec.z);
        
        // the normalization of the particle
        float normfac;
        
        // the angularsize^2
        float d2 = dtheta * dtheta;
        
        // calculate normfac
        {
            //if(usenormmap == 0 && newsize != 1.0){
            if(newsize != 1.0){
#ifdef USEACTURALNORM
                //use actual norm
                normfac = calc_norm(vec2(xc, yc), newsize, dtheta);
#else
                //use analytical norm
				normfac = 1.0 / (geofac.y * geofac.y) / calc_norm1(dtheta) * 4.0;
#endif
            }else{
                normfac = 1.0;
            }
        }
        
        //Must add another vector (xc, yc)
        //flux = 1.0;
        
        // the color is used to transfor the value of the center and norm
        gl_FrontColor = vec4(xc, yc, flux * normfac , dtheta);

        // textura
        // gl_TexCoord[0] = gl_MultiTexCoord0;
    }else{
        // if the particles are in the upper sphere, remove this particle
        
        gl_PointSize = 1.0;  //point size
        newpos = vec4(0.0, 0.0, 0.0, 1.0);
        gl_FrontColor = vec4(0, 0, 0, 0);
        // gl_TexCoord[0] = gl_MultiTexCoord0;
    }
    
    // calculate the final position of the particle on the screen
    gl_Position = gl_ModelViewProjectionMatrix * newpos;
    
    
    gl_TexCoord[0] = gl_MultiTexCoord0;
    
}
