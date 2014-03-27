#ifndef __LY__RENDER__
#define __LY__RENDER__
#include "buffers.h"
#include "shaders.h"
#include "../parameters.h"

//the class for rendering a particle
struct RenderParticle{
    // w is used for other flags
    // e.g. the flag for subhalo/halo
    float x, y, z, w;
    
    //final flux is calculated  as fluxfac1 * fluxfac2 * fluxfac3/ 4*pi*r^2
    float fluxfac1;
    float fluxfac2;
    float fluxfac3;
    
     //smooth radius
    float hsmooth;
};


class render{
public:
    //DataReader * reader;
    //Parameter * params;
    
    // physfac, scalefac, usemask: please read the shader
    render(Parameters &par, int imsize, int pointSize = 256, int nside = 512,
           int physFac = 0, int isUseMask = 0, float scalefac =  1.0);
    ~render();
    
    double * getHealPixMap();
    float * getUpSphere();
    float * getLowerSphere();

    
    void rend(RenderParticle * fluxdata, int numParts);
    void clear();       //clear the current map
    
private:
    
    Parameters * par_;
    
    //L for the lower sphere
    //U for the upper sphere
    fluxShader * fshaderL;
//    colorShader * cshaderL;
    buffer * fbufferL;
//    colorBuffer * cbufferL;
    
    fluxShader * fshaderU;
//    colorShader * cshaderU;
    buffer * fbufferU;
//    colorBuffer * cbufferU;
    
    int nside_;
    
    
    //GLuint textureIni;  //initial empty texture

    unsigned int windowSize;      //2^m
    unsigned int pointSize;
    float orthsize;                //orthogonal view port size
    
    void init();
    void drawFlux(RenderParticle * fluxdata, int numParts);
                                  //    void drawImage();

    bool initialized;             //is initialized?
    
    void openGLInit();
    
    float * fluxmapL, *fluxmapU;   //L and U map
    double * healpixmap_;
    
    void readFluxMap();           //read map from the GPU
    
    //get the flux at (x, y, z) on the sphere using linear interpolation
    double _getpixflux(double x, double y, double z);
    double _getPixel(double x, double y, bool isUp);
    bool _getTexCoord(double &x, double &y, double d, double isupdown);
    
    void convertToHealpixMap();
    
    double totalFlux;

};

#endif
