#include <cmath>
#include <rpc/types.h>
#include <rpc/xdr.h>
//#include "types.h"
#include "shaders.h"

using namespace std;
void fluxShader::loadUniform(){

    rotmloc = glGetUniformLocationARB(shader->progmObj,"rotmatrix");
    //printf("rotmloc, %d\n", rotmloc);
    
    oposloc = glGetUniformLocationARB(shader->progmObj,"opos");
    //printf("oloc, %d\n", oposloc);

    geofacloc = glGetUniformLocationARB(shader->progmObj,"geofac");
    //printf("gloc, %d\n", geofacloc);
    
    physfacLoc = glGetUniformLocationARB(shader->progmObj,"physfac");
    //printf("gloc, %d\n", geofacloc);
    
    scaleFacLoc = glGetUniformLocationARB(shader->progmObj,"scaleFac");
    //printf("gloc, %d\n", geofacloc);
    
    isUseMaskLoc = glGetUniformLocationARB(shader->progmObj,"isUseMask");
    //printf("gloc, %d\n", geofacloc);

    
}

fluxShader::fluxShader(){
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

void fluxShader::setopos(Parameters &pars){
    //printf("%f %f %f\n", pars.params.opos[0], pars.params.opos[1], pars.params.opos[2]);
    setopos3f(pars.params.opos[0], pars.params.opos[1], pars.params.opos[2]);
}


void fluxShader::setrotm(Parameters &pars, bool updown){
    for (int i = 0; i<3; i++){
        for (int j =0; j<3; j++){
            rotmatrix[i * 3 + j] = pars.rotmatrix[i][j];
        }
    }
    
    if(updown){
        rotmatrix[2] *= -1;
        rotmatrix[5] *= -1;
        rotmatrix[8] *= -1;
    }
}

//void fluxShader::setrotmatrix(float * alignvec, float * obsvec, float * centvec, bool updown){
void fluxShader::setrotmatrix(Parameters &pars, bool updown){
    setrotm(pars, updown);
    glUniformMatrix3fvARB(rotmloc, 1, GL_FALSE, this->rotmatrix);
    
}

void fluxShader::setopos3fv(float * ya){
    glUniform3fARB(oposloc, ya[0], ya[1], ya[2]);
}


void fluxShader::setPhysicsFac(int physfac){
    glUniform1fARB(physfacLoc, physfac);
}

// wether use mask? 0 no, 1 yes
// if yes, use position.w as mask
// if w > 0.5, the point keep
// if w < 0.5, the point dispose
void fluxShader::setUseMask(int isUseMask){
    glUniform1iARB(isUseMaskLoc, isUseMask);
}

// set a scale factor appling to flux
void fluxShader::setScaleFac(float scaleFac){
    glUniform1fARB(scaleFacLoc, scaleFac);
}


void fluxShader::setopos3f(float x, float y, float z){
    //printf("Opos: %f %f %f\n", x, y, z);
    glUniform3fARB(oposloc, x, y, z);
}

void fluxShader::setgeofac3fv(float * ya){
    glUniform3fARB(geofacloc, ya[0], ya[1], ya[2]);
}

void fluxShader::setgeofac3f(float x, float y, float z){
    glUniform3fARB(geofacloc, x, y, z);
}



