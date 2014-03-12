#include <iostream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <fstream>

//fits and HEALPIX
#include <healpix_base.h>
#include <healpix_map.h>
#include <arr.h>
#include <fitsio.h>
#include <fitshandle.h>
#include <healpix_map_fitsio.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "buffers.h"
#include "shaders.h"
#include "render.h"

#include "types.h"

#ifdef __APPLE__
#include <glew.h>
#include <GLUT/glut.h> // darwin uses glut.h rather than GL/glut.h
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif


// the margin to handle the boundary problem
#define CANVASMARGIN 0.0

// set to the edge of the OGL convas, to avoid edge problem

using namespace std;

//colorBuffer * CB, *CBL, *CBU;       //for final rendering
static unsigned int WSIZE, POINTSIZE;
string picfile;
bool isonscreen = false;

void render::openGLInit(){
    int argc = 1;
    char *argv[1] = {(char*)"Particle Render"};
    
    // initialize glut
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
    glutInitWindowSize(WSIZE, WSIZE);
    glutCreateWindow("Dark Matter GammaRay rendering!");
    
    // initialize glut
    glewExperimental = GL_TRUE;
    glewInit();
    
    // hide the window
    glutHideWindow();
    
    
    // check whether GLSL is supported
    if (GLEW_ARB_vertex_shader && GLEW_ARB_fragment_shader)
		printf("Ready for GLSL\n");
    else {
		printf("No GLSL support\n");
		exit(1);
    }
    
    //set up shaders and buffers
    if(!initialized){
        
        fshaderL = new fluxShader();
        fbufferL = new buffer(WSIZE, WSIZE);
        fbufferL->setBuffer();
        
        fshaderU = new fluxShader();
        fbufferU = new buffer(WSIZE, WSIZE);
        fbufferU->setBuffer();
        
        initialized = true;
    }
    
    //initialize enviroment
    init();
    
    //setup blending
    glEnable (GL_BLEND);
    glBlendFunc (GL_ONE,GL_ONE);    //blending
    
}




render::render(Parameters &par, int imsize, int pointSize, int nside){
    par_ = & par;
    
    windowSize = imsize;
    pointSize = pointSize;
    
    WSIZE = windowSize;
    POINTSIZE = pointSize;
    
    initialized = false;
    orthsize = 1.0;
    totalFlux = 0.0;
    
    nside_ = nside;
    
    fluxmapL = new float[WSIZE*WSIZE];
    fluxmapU = new float[WSIZE*WSIZE];
    healpixmap_ = new double[12 * nside * nside];
    
    openGLInit();
    
    init();
    
    
    //setup shaders L and U
    {
        //begin shader
        fshaderL->begin();
        //fshaderL->setIsUseRotm(IS_USE_NORM);
        
        //setup shader parameters
        fshaderL->setgeofac3f(orthsize, windowSize, pointSize);
        
        //now use the rotatation matrix
        fshaderL->setopos(*par_);
        fshaderL->setrotmatrix(*par_, false);
        fshaderL->end();
        
        //begin shader
        fshaderU->begin();
        //fshaderU->setIsUseRotm(IS_USE_NORM);
        
        //setup shader parameters
        fshaderU->setgeofac3f(orthsize, windowSize, pointSize);
        
        //now use the rotatation matrix
        fshaderL->setopos(*par_);
        fshaderL->setrotmatrix(*par_, true);
        fshaderU->end();
    }
    

    
    
    fbufferL->bindBuf();
    {
        //clear color
        glClearColor (0.0, 0.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }
    fbufferL->unbindBuf();
    
    
    fbufferU->bindBuf();
    {
        //clear color
        glClearColor (0.0, 0.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }
    fbufferU->unbindBuf();
};

render::~render(){
    delete[] fluxmapL;
    delete[] fluxmapU;
    delete[] healpixmap_;
    
    delete fshaderL;
    delete fshaderU;
    delete fbufferU;
    delete fbufferL;
}

void render::init(){
    
    // enable texture
    glEnable(GL_TEXTURE_2D);
    glDisable(GL_DEPTH_TEST);

    // setup some generic opengl options
    // disable clamp
    glClampColorARB(GL_CLAMP_VERTEX_COLOR_ARB, GL_FALSE);
    glClampColorARB(GL_CLAMP_FRAGMENT_COLOR_ARB, GL_FALSE);
    glClampColorARB(GL_CLAMP_READ_COLOR_ARB, GL_FALSE);

    // enable point sprite
    glEnable(GL_POINT_SPRITE);
    glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
    // enable gl_PointSize in vertex shader
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
   
    
    // set Point Size
    GLfloat sizes[2];
    glGetFloatv(GL_ALIASED_POINT_SIZE_RANGE, sizes);
	glPointParameterfARB( GL_POINT_FADE_THRESHOLD_SIZE_ARB, 1000.0 );
	glPointParameterfARB( GL_POINT_SIZE_MIN_ARB, sizes[0]);
	glPointParameterfARB( GL_POINT_SIZE_MAX_ARB, sizes[1]);

    // set the coordinates to be lower left to avoid distortion
    glPointParameteri(GL_POINT_SPRITE_COORD_ORIGIN, GL_LOWER_LEFT);

    glPointSize(pointSize);
    
    
    //set viewport
    glViewport(0,0,WSIZE, WSIZE);
    
    //setup matrix
    // setup orthgonal matrix to projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-orthsize, orthsize, -orthsize, orthsize, -100, 100);
    
    // setup identity to the modelview
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void render::drawFlux(RenderParticle * fluxdata, int numParts){

        GLfloat * vetexarray = (GLfloat *) fluxdata;
        glEnableClientState (GL_VERTEX_ARRAY);
        glEnableClientState (GL_COLOR_ARRAY);
        glColorPointer (3, GL_FLOAT, 6*sizeof(GLfloat), &(vetexarray[3]));
        glVertexPointer (3, GL_FLOAT, 6*sizeof(GLfloat), &(vetexarray[0]));

        //test
        //for(int i = 0; i < numParts)
        //printf("%e %e %e %e %e %e\n", vetexarray[0], vetexarray[1], vetexarray[2],
        //                vetexarray[3], vetexarray[4], vetexarray[5]);
    
        //lower sphere
        fbufferL->bindBuf();
        fshaderL->begin();
        {
            glDrawArrays(GL_POINTS, 0, numParts);
        }
        fshaderL->end();
        fbufferL->unbindBuf();
        
        //upper sphere
        fbufferU->bindBuf();
        fshaderU->begin();
        {
            glDrawArrays(GL_POINTS, 0, numParts);
        }
        fshaderU->end();
        fbufferU->unbindBuf();
        
        glDisableClientState (GL_VERTEX_ARRAY);
        glDisableClientState (GL_COLOR_ARRAY);
    


    glBindTexture(GL_TEXTURE_2D, 0);

    glPopAttrib();
}
 
void render::readFluxMap(){
    glPixelStorei(GL_PACK_ALIGNMENT, 4);  


    fbufferL -> bindTex();
    glGetTexImage(GL_TEXTURE_2D,0,GL_RED,GL_FLOAT,fluxmapL);
    fbufferL->unbindTex();
    
    fbufferU -> bindTex();
    glGetTexImage(GL_TEXTURE_2D,0,GL_RED,GL_FLOAT,fluxmapU);
    fbufferU->unbindTex();
}

void render::rend(RenderParticle * fluxdata, int numParts){
    drawFlux(fluxdata, numParts);
}


//get the pixel value
double render::_getPixel(double x, double y, bool isUp){
    double d = (double) WSIZE / 2 - CANVASMARGIN;
    isUp = _getTexCoord(x, y, d, isUp);
    int ind = (int)floor((WSIZE / 2 - y - 1) * WSIZE + x + WSIZE / 2);
    
    if(isUp){
        return fluxmapU[ind];
    }
    else{
        return fluxmapL[ind];
    }
}


// return is up down
bool render::_getTexCoord(double &x, double &y, double d, double isupdown){
    double cx = (x + 0.5) / d;
    double cy = (y + 0.5) / d;
    if(cx*cx + cy*cy <= 1.0){
        return isupdown;
    }else{
        double dnom = 1.0 + cx * cx + cy * cy;
        double vx = (2 * cx) / dnom;
        double vy = (2 * cy) / dnom;
        double vz = (dnom - 2) / dnom;
        vz = - vz;
        cx = vx / (1 - vz);
        cy = vy / (1 - vz);
        
        x = floor(cx * d);
        y = floor(cy * d);
        return !isupdown;
    }
};

double render::_getpixflux(double x, double y, double z){
    // calculate is upper or lower sphere
    bool isupshere = false;
    if(z > 0.0){
        isupshere = true;
        
        // reverse the z-direction
        z = - z;
    }
    
    double d = WSIZE / 2.0 - CANVASMARGIN;
    
    //bilinear interpolation
    double pxc = x/(1-z);
    double pyc = y/(1-z);
    
    double xc = (pxc) * d;
    double yc = (pyc) * d;
    
    double x1 = floor(xc);
    
    double x2 = x1 + 1;
    
    double y1 = floor(yc);
    
    double y2 = y1 + 1;
    
    //return _getPixel(x1, y1, isupshere);
    
    double f11, f12, f21, f22;
    f11 = _getPixel(x1, y1, isupshere);
    f12 = _getPixel(x1, y2, isupshere);
    f21 = _getPixel(x2, y1, isupshere);
    f22 = _getPixel(x2, y2, isupshere);
    
    double flux = 0;
    double fr1 = (x2 - xc) / (x2 - x1) * f11 + (xc - x1) / (x2 - x1) * f21;
    double fr2 = (x2 - xc) / (x2 - x1) * f12 + (xc - x1) / (x2 - x1) * f22;
    flux = (y2 - yc) / (y2 - y1) * fr1 + (yc - y1) / (y2 - y1) * fr2;
    return flux;
}

void render::convertToHealpixMap(){
    double * healmap;
    int nside = nside_;
    //int pixss = WSIZE * WSIZE;
    int npix = nside * nside * 12;
    healmap = healpixmap_;
    
    double _rffmin = 1.0e36;
    double _rffmax = 0.0;
    double total_f = 0.0;
    
    Healpix_Base base(nside, RING, SET_NSIDE);

    
    for(int i = 0; i < npix; i++){
        //double x, y,r, factor;
        vec3 this_vec = base.pix2vec(i);

        //make a rotation on y to match the data stored in memory
        this_vec.x = - this_vec.x;
        double flux = _getpixflux(this_vec.x, this_vec.y, this_vec.z);
        
        if(this_vec.z > 0){
            this_vec.z = - this_vec.z;
        }
        double pr = sqrt(1 - this_vec.z * this_vec.z)/(1-this_vec.z);
        healmap[i] = flux / (4.0 / (1 + pr*pr)/(1 + pr*pr))
            * WSIZE * WSIZE / 4.0 *  par_->map.dOmega;;
        
    }
}

void render::clear(){
    memset(fluxmapL, 0, WSIZE * WSIZE * sizeof(float));
    memset(fluxmapU, 0, WSIZE * WSIZE * sizeof(float));
    
    fbufferL->bindBuf();
    {
        //clear color
        glClearColor (0.0, 0.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }
    fbufferL->unbindBuf();
    
    fbufferU->bindBuf();
    {
        //clear color
        glClearColor (0.0, 0.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }
}


double* render::getHealPixMap(){
    readFluxMap();
    convertToHealpixMap();
    return healpixmap_;
}

float* render::getUpSphere(){
    readFluxMap();
    return fluxmapU;
}
float* render::getLowerSphere(){
    return fluxmapL;
}
