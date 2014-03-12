/*****************************************************************************/
// shaders.h
// defines the actual shader
// manages a shader
//  Author: Lin F. Yang
//  01/18/2014
/*****************************************************************************/

#ifndef __LY_PARTICLESHADER__
#define __LY_PARTICLESHADER__
#include <string>
#include <cstdio>
#include <stdlib.h>
#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#endif
#include "../types.h"
#include "shaderManager.h"
#include "../parameters.h"

using namespace std;

class Shader{
protected:
    shaderObj * shader;
	string vfile, ffile;
    bool good;
    
public:
    bool is_good(){
        return good;
    };
    
    //start shader
    void begin(){
        shader->begin();
    };
    
    //end shader
    void end(){
        shader->end();
    };
    
    void setfile(string v, string f){
        vfile =v;
        ffile =f;
    };
    
};



// the flux shader, used to render the DM gamma flux
class fluxShader:public Shader{
public:
    static const string VERTEX_SHADER_SRC;
    static const string FRAGMENT_SHADER_SRC;
    
    //set opos
    void setopos3f(REAL x, REAL y, REAL z);
    void setopos3fv(REAL * opos);
    void setopos(Parameters &pars);
    
    //set geometric factor
    void setgeofac3f(REAL x, REAL y, REAL z);
    void setgeofac3fv(REAL * geo);
    
    void setrotmatrix(Parameters &pars, bool updown);
    
    void setIsUseRotm(bool isrotm){
        isRotm = isrotm;
    };
    
    //load shader from file
    fluxShader();
    
private:
    //GLint xaxisloc, yaxisloc, zaxisloc;
    GLint rotmloc;
    GLint oposloc;
    GLint geofacloc;
    
    //setup rotation matrix, 9 variables
    //updown: true/up; false/down
    void setrotm(Parameters &pars, bool updown);
    
    void loadUniform();

    REAL * opos;

    REAL rotmatrix[9];
    
    //set whether use the rotation matrix
    bool isRotm;
};



#endif
