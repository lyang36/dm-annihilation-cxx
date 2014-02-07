#ifndef __LY__RENDER__
#define __LY__RENDER__
#include "buffers.h"
#include "shaders.h"
#include "../parameters.h"

//the class for rendering a particle
struct RenderParticle{
    float x, y, z;
    
    //final flux is calculated  as densityfac1 * densityfac2 / 4*pi*r^2
    float densityfac1;
    float densityfac2;      //angular radius
    float hsmooth;
};


class render{
public:
    //DataReader * reader;
    //Parameter * params;
    
    render(Parameters &par, int imsize, int pointSize = 256, int nside = 512);
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
    fluxBuffer * fbufferL;
//    colorBuffer * cbufferL;
    
    fluxShader * fshaderU;
//    colorShader * cshaderU;
    fluxBuffer * fbufferU;
//    colorBuffer * cbufferU;
    
    int nside_;
    
    
    GLuint textureIni;  //initial empty texture

    unsigned int windowSize;      //2^m
    unsigned int pointSize;
    REAL orthsize;                //orthogonal view port size
    
    void init();
    void drawFlux(RenderParticle * fluxdata, int numParts);
                                  //    void drawImage();

    bool initialized;             //is initialized?
    
    void openGLInit();
    
    REAL * fluxmapL, *fluxmapU;   //L and U map
    double * healpixmap_;
    
    void readFluxMap();           //read map from the GPU
    
    double _getpixflux(int x1, int y1, bool isupshere);
    void convertToHealpixMap();
    
    double totalFlux;

};

#endif