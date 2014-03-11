/*****************************************************************************/
// The fragment shader for a point sprite to draw projected gaussian profile
// Author: Lin F. Yang
// Date: 01/18/2014
/*****************************************************************************/
#version 120

//Constants
#define PI 3.1415926535897932

/******************************Parameters*************************************/
// orthsize, windowsize, pointsize
uniform vec3 geofac;

// The normalization map
uniform sampler2D normmap;

// whether use the norm map? true: 1 else:0
// uniform int usenormmap;

// The particle parameter passing from flux.vert
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
    
    /*if(d2 > 1.0){
        return 0.0;
    }*/
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


void main(){
    
    // the new size of the pointsprite
    float newsize = particle.r;
    
    // the center coordinates of the point sprite
    vec2 xyc = gl_Color.rg;
    
    // convert the coordinate to actuall screen coordinates
    vec2 coor = xyc *geofac.y / 2.0;
    
    // initialize
    float flux = 0.0;
    
    // flux fac is the normalization calculated in vert
    float fluxfac = gl_Color.b;
    
    // rasterizing the point sprite
	if(newsize > 1.0){
		if(fluxfac > 0.0){
            // the angular radius
			float dtheta = gl_Color.a;
            
            // get current relative coordinates of the pixel inside the pointsprite
			vec2 p = floor(newsize * vec2(gl_PointCoord.s, gl_PointCoord.t));
			
            // because of each coordinates of the pixels
            // is the lower left corner, must be converted to be the
            // center of the pixel, then convert to the coordinates
            // to be in a circle of radias 1
			p = (p+0.5) / newsize;
			p = 2.0*(p-0.5);
            
            // if the pixel is outside the unit circle, discard it
            // this is because of each point sprite is actually
            // a ractangle.
			float u = dot(p, p);
			if (u > 1.0)
                discard;
			
            // calculate the actuall coordinates of a pixel on the screen
			vec2 xyp = p * (newsize / 2.0) + coor;
			vec2 xyr = xyp / (geofac.y / 2.0);
            
            // calculate the r-coordinates
			float pr2 = dot(xyr, xyr);
            
            // use the actual norm
			flux = fluxfac  * profPRJ(prev(xyr), dtheta);
            
            //profile(prev(xyr), dtheta) * 4.0/(1.0+pr2)/(1.0+pr2);
            
            //use analytical norm
            //flux = fluxfac  * profile(prev(xyr), dtheta) * 4.0/(1.0+pr2)/(1.0+pr2);
			//flux = fluxfac;
			
			
            // output the flux into the red component of color
			gl_FragColor = vec4(flux, 0, 0, 1.0);
		}else{
			discard;
		}
	}else{
		gl_FragColor = vec4(fluxfac, 0, 0, 1.0);
	}
}
