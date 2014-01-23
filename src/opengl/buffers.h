#ifndef __LY_BUFFERS__
#define __LY_BUFFERS__
#include <string>
//#include <glew.h>
//#include <GLUT/glut.h>
#ifdef __APPLE__
#include <glew.h>
#include <GLUT/glut.h> // darwin uses glut.h rather than GL/glut.h
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif
#include "shaders.h"

using namespace std;

class buffer{              //buffers
protected:
    GLuint textureId;       //texture binding
    bool fboUsed;          //fbo ok?
    GLuint fboId, rboId;   //framebuffer and renderbuffer ID
    unsigned int twidth, theight; //the width and height of the buffer
    void checkbuffer();     //check whether buffer is ok
    bool  issetup;
    
public:     
    buffer(unsigned int w, unsigned int h);
    ~buffer();
    GLuint getTex();     //get texture ID
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
    string normfile;// = "norm.dat";
    float * normtextbuf;
    //bool isUseMap;
public:
    //load norm and bind texture
    void setNormTex(){
        normtextbuf = new float[normMapRes * normMapRes];
        loadnorm();
        glEnable(GL_TEXTURE_2D);
        glGenTextures(1, &normtex);
        glBindTexture(GL_TEXTURE_2D, normtex);
        
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, normMapRes, normMapRes, 0, GL_RED, GL_FLOAT, normtextbuf);
        //for(int i = 0; i < normMapRes*normMapRes; i++){
        //    printf("%f  ", normtextbuf[i]);
        //}
        // set its parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
        
        // setup some generic opengl options
        glClampColorARB(GL_CLAMP_VERTEX_COLOR_ARB, GL_FALSE);
        glClampColorARB(GL_CLAMP_FRAGMENT_COLOR_ARB, GL_FALSE);
        glClampColorARB(GL_CLAMP_READ_COLOR_ARB, GL_FALSE);
        
    };
    
    void setMapRes(int m, int n){
        normMapRes = m;
        normPointSize = n;
    };
    
    fluxBuffer(unsigned int w, unsigned int h):buffer(w,h){
        normtex = 100;
        normMapRes = 1024;
        normPointSize = 256;
        normfile = "norm.dat";
        normtextbuf = NULL;
        //isUseMap = false;
    };//:buffer(w,h);
    
    //bool isBufferLoaded(){
    //    return isUseMap;
    //}
    
    ~fluxBuffer(){
        if(normtextbuf != NULL){
            delete normtextbuf;
            normtextbuf = NULL;
        }
    }
    
};                 


#endif
