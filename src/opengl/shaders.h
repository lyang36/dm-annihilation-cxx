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
    
    void begin(){
        shader->begin();
    };           //start shader
    
    void end(){
        shader->end();
    };//end shader
    
    void setfile(string v, string f){
        vfile =v;
        ffile =f;
    };
    
    
};


class fluxShader:public Shader{

public:
    static const string VERTEX_SHADER_SRC;
    static const string FRAGMENT_SHADER_SRC;
    
    
    void setopos3f(REAL x, REAL y, REAL z); //set opos
    void setopos3fv(REAL * opos);
    void setopos(Parameters &pars);
    
    void setgeofac3f(REAL x, REAL y, REAL z); //set opos
    void setgeofac3fv(REAL * geo);
    //void setrotmatrix(REAL * alignvec, REAL * obsvec, REAL * centervec, bool updown);     //setup rotation matrix
    
    void setrotmatrix(Parameters &pars, bool updown);
                                            //alignvec has 3 vars
    

    
    void setusenormmap(bool isnm);          //whether use the norm map
    void setIsUseRotm(bool isrotm){
        isRotm = isrotm;
    };
    
    
    //load shader from file
    fluxShader(){
        setIsUseRotm(false);
        shaderManager SM;
        shader = SM.loadShader(VERTEX_SHADER_SRC.c_str(), FRAGMENT_SHADER_SRC.c_str());
        if(shader == 0){
            printf("Shader error!\n");
            exit(1);
        }
        if(shader->progmObj == 0){
            printf("Shader error!\n");
            exit(1);
        }else{
            good = true;
        }
        
        if(good){
            loadUniform();
        }
        
        isRotm = true;
    };
    
private:
    //GLint xaxisloc, yaxisloc, zaxisloc;
    GLint rotmloc;
    GLint oposloc;
    GLint geofacloc;
    GLint isnormmaploc;
    void setrotm(Parameters &pars, bool updown);  //setup rotation matrix, 9 variables
    //updown: true/up; false/down
    void loadUniform();
    //REAL * align_vec;
    REAL * opos;
    //REAL * cpos;
    REAL rotmatrix[9];
    
    bool isRotm;                //set whether use the rotation matrix
};



#endif
