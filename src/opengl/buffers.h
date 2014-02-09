/*****************************************************************************/
// buffers.h
// defines the pixel buffer for
// buffer rendering or texture
//  Author: Lin F. Yang
//  01/18/2014
/*****************************************************************************/

#ifndef __LY_BUFFERS__
#define __LY_BUFFERS__
#include <string>
#ifdef __APPLE__
#include <glew.h>
#include <GLUT/glut.h> // darwin uses glut.h rather than GL/glut.h
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif
#include "shaders.h"

using namespace std;

class buffer{               //buffers
protected:
    GLuint textureId;       //texture binding
    bool fboUsed;           //fbo ok?
    GLuint fboId, rboId;    //framebuffer and renderbuffer ID
    unsigned int twidth, theight; //the width and height of the buffer
    void checkbuffer();     //check whether buffer is ok
    bool  issetup;
    
public:     
    buffer(unsigned int w, unsigned int h);
    ~buffer();
    GLuint getTex();         //get texture ID
    void bindTex();          //bind the texture
    void unbindTex();        //unbind texture
    void setTex(GLuint tex);       //set a texture, bind it to the buffer
    void genTex();           //generate a new texture
    void attachTex();        //attach the texture
    void genBuffer();        //generate the frame buffer and render buffer
                             //and attach them together and bind buffer
    void setBuffer();        //automatically set up all the things
    
    void bindBuf();          //bind buffer to render
    void unbindBuf();        //unbind buffer from render
    
    void readTex(void * tex);//read the data from the GPU
    void saveImg(string filename);  //save the picture to file
};

class fluxBuffer:public buffer{       //buffer for additing flux
private:
    GLuint normtex;          //map used to calculate map flux
    int normMapRes;          //map resolution, in pixels.
    int normPointSize;       //max point size of the point sprite, in pixels
    void loadnorm();         //calculate normalization text or from file
    string normfile;         // = "norm.dat";
    float * normtextbuf;
                             //bool isUseMap;
public:
    
    void setNormTex();       //load norm and bind texture
                             //(using the normalization map)
    
                             //set the normalization map resolution
    void setMapRes(int m, int n);
    
                             //initialization
    fluxBuffer(unsigned int w, unsigned int h);
    
    ~fluxBuffer();
    
};                 


#endif
