/*****************************************************************************/
// The fragment shader for a point sprite to draw projected gaussian profile
// Author: Lin F. Yang
// Date: 01/18/2014
/*****************************************************************************/
//Constants
#define PI 3.1415926535897932

/******************************Parameters*************************************/
// orthsize, windowsize, pointsize
uniform vec3 geofac;

// The normalization map
uniform sampler2D normmap;

// whether use the norm map? true: 1 else:0
uniform int usenormmap;

// The particle parameter passing from flux.vert
// vec4(newsize, npvec.x, npvec.y, npvec.z);
varying vec4 particle;


// The output profile of a particle
// This is very important, must be checked
float profile(vec3 r1,float dtheta){
    // get the position of the partilce on the sky
    vec3 r0 = vec3(particle.gba);
    
    // calcualte the angle between the pixel and the particle
    float costheta = dot(r0, r1)/(length(r0)*length(r1));
    

    costheta = clamp(costheta, -1.0, 1.0);
    
    // use Taylor seriers
    // acos may has too much error
    float t2 = 2.0 * ( 1.0 - costheta) + 1.0/3.0*(costheta - 1.0)*(costheta - 1.0) - 4.0/45.0 * (costheta - 1.0) *(costheta - 1.0)*(costheta - 1.0);
    //costheta = clamp(costheta, -1.0, 1.0);
    //float t2 = acos(costheta);
    //t2 = t2*t2;
    
    
    float d2 = clamp(t2 / dtheta / dtheta, 0.0, 1.0);
    //float d2 = (t2 / dtheta / dtheta, 0.0, 1.0);
    
    if(t2 > 1.0){
        return 0.0;
    }
    if(t2 < 0.0){
        t2 = 0.0;
    }
    return exp(- 1.5 * d2);         //here comes the problems
}

// reverse stereoprojection
vec3 prev(vec2 xy){
    float r2 = xy.x*xy.x + xy.y*xy.y;
    return vec3(2.0 * xy.x/(1.0 + r2), 2.0 * xy.y/(1.0 + r2), (r2 - 1.0)/(r2 + 1.0));
}

// stereoproject a profile to the plane
float projprofile(vec2 xy, float fc, float dtheta){
    return fc * profile(prev(xy), dtheta);
}

void main(){
    
    // the new size of the halo
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
			vec2 p = floor(newsize * vec2(gl_TexCoord[0].s,gl_TexCoord[0].t));
			
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
			if (u > 1.0) discard;
			
            // calculate the actuall coordinates of a pixel on the screen
			vec2 xyp = p * (newsize / 2.0) + coor;
			vec2 xyr = xyp / (geofac.y / 2.0);
            
            // calculate the r-coordinates
			float pr2 = dot(xyr, xyr);
            
            // use the actual norm
			flux = fluxfac  * profile(prev(xyr), dtheta) * 4.0/(1.0+pr2)/(1.0+pr2);
            
            //use analytical norm
            //flux = fluxfac  * profile(prev(xyr), dtheta) * 4.0/(1.0+pr2)/(1.0+pr2);
			//flux = fluxfac;
			
            // use a normalization map, stored in the texture
            if(usenormmap == 1){
				float r0 = sqrt(pr2);
				float r = newsize / geofac.z;
				float norm = (texture2D(normmap, vec2(r0, r))).r;
				flux = flux / norm;
			}
			
            // output the flux into the red component of color
			gl_FragColor = vec4(flux, 0, 0, 1.0);
		}else{
			discard;
		}
	}else{
		gl_FragColor = vec4(fluxfac, 0, 0, 1.0);
	}
}
