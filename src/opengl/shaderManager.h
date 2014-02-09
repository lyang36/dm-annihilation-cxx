/*****************************************************************************/
// shaderManager.h
// defines the shader manager
// manages a shader
//  Author: Lin F. Yang
//  01/18/2014
/*****************************************************************************/

#ifndef __LY_SHADERMANAGER__
#define __LY_SHADERMANAGER__
#include <cstdlib>
#include <cstdio>
#ifdef __APPLE__
#include <glew.h>
#include <GLUT/glut.h> // darwin uses glut.h rather than GL/glut.h
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif

class shaderObj{
public: 
    GLhandleARB vertexShader;
    GLhandleARB fragmentShader;
    GLhandleARB progmObj;
    void begin();
    void end();
    
};

class shaderManager{
private:
    void printLog(GLhandleARB obj);
public:
    shaderObj * loadShaderFile(const char * vfile, const char * ffile);
    shaderObj * loadShader(const char *vf, const char * ff);
    static char * textFileRead(const char * fn, int * length = NULL);
    
};


#endif
