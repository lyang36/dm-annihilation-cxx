#include <stdlib.h>
#include <cmath>

#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>


#ifdef __APPLE__
#include <glew.h>
#include <GLUT/glut.h> // darwin uses glut.h rather than GL/glut.h
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "buffers.h"

using namespace std;

void buffer::checkbuffer(){
    // check FBO status
    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    
    if(status != GL_FRAMEBUFFER_COMPLETE){
        fboUsed = false;
        printf("No frame buffer!\n");
        exit(-1);
    }
}

buffer::buffer(unsigned int w, unsigned int h){
    twidth = w;
    theight = h;
}


buffer::~buffer(){
    glDeleteFramebuffers(1, &fboId);
    glDeleteRenderbuffers(1, &rboId);
    glDeleteTextures(1, &textureId);
}

GLuint buffer::getTex(){
    return textureId;
}

void buffer::bindTex(){
    glBindTexture(GL_TEXTURE_2D, textureId);
}

void buffer::unbindTex(){
    glBindTexture(GL_TEXTURE_2D, 0);
}

void buffer::setTex(GLuint tex){
    textureId = tex;
}

void buffer::genTex(){
    glGenTextures(1, &textureId);
    bindTex();
    
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, twidth, theight, 0,
                 GL_RGBA, GL_UNSIGNED_BYTE, 0);
    
    // set its parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    // setup some generic opengl options
    glClampColorARB(GL_CLAMP_VERTEX_COLOR_ARB, GL_FALSE);
    glClampColorARB(GL_CLAMP_FRAGMENT_COLOR_ARB, GL_FALSE);
    glClampColorARB(GL_CLAMP_READ_COLOR_ARB, GL_FALSE);
    
    unbindTex();

}


void buffer::attachTex(){
    // attach the texture to FBO color attachment point
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D, textureId, 0);

}


void buffer::genBuffer(){
    // create a renderbuffer object to store depth info
    //GLuint rboId;
    glGenRenderbuffers(1, &rboId);
    glBindRenderbuffer(GL_RENDERBUFFER, rboId);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT,
                          twidth, theight);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    
    // create a framebuffer object
    //GLuint fboId;
    glGenFramebuffers(1, &fboId);
    //glBindFramebuffer(GL_FRAMEBUFFER, fboId);
    bindBuf();
    // attach the renderbuffer to depth attachment point
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, rboId);
    //unbindBuf();
    printf("Buffer generated!\n");
    //printf("frame: %d, render: %d \n", fboId, rboId);
}

void buffer::setBuffer(){
    genTex();
    genBuffer();
    attachTex();
    checkbuffer();
    unbindBuf();
}

void buffer::bindBuf(){

    glBindFramebuffer(GL_FRAMEBUFFER, fboId);
}

void buffer::unbindBuf(){
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void buffer::readTex(void* tex){
    
}

void buffer::saveImg(string filename){
    
}



float profile(float * r, float * r0, float dtheta){ 
    
    float costheta = r[0] * r0[0] + r[1] * r0[1] + r[2] * r0[2];
    if(costheta > 1.0){
        costheta = 1.0;
        //printf("bbb\n");
    }
    if(costheta < -1.0){
        costheta = -1.0;
    }

    float t2 = acos(costheta);
    t2 = t2*t2;
    float d2 = t2 / dtheta / dtheta;
    /*if(exp(- 1.5 * d2) == 0.0){
        printf("d2 %f, dtheta %f\n", d2, dtheta);
    }*/
    return exp(- 1.5 * d2);
}

//reverse stereoprojection
void prev(float * xy, float * r){
    float r2 = xy[0]*xy[0] + xy[1]*xy[1];
    r[0] = 2.0 * xy[0]/(1.0 + r2);
    r[1] = 2.0 * xy[1]/(1.0 + r2);
    r[2] = (r2 - 1.0)/(r2 + 1.0);
}

float projprofile(float * xy, float * r0, float fc, float dtheta){
    float r[3];
    prev(xy, r);
    /*if(profile(r, r0, dtheta) == 0.0){
        printf("r %f %f %f, r0 %f %f %f\n", r[0], r[1], r[2], r0[0], r0[1], r0[2]);
    }*/
    return fc * profile(r, r0, dtheta);
}


void fluxBuffer::loadnorm(){
    //check whether the file is exsisted
    //if so, read from file
    
    ifstream infile;
    infile.open(normfile.c_str(), ios::in | ios::binary);
    if(infile.good())
    {
        int tempsize;
        infile.read((char *) &tempsize, sizeof(int));
        if(tempsize != normMapRes){
            //the file do not match
            printf("Current norm file %s do not match required, please delete or move it.\n", normfile.c_str());
            infile.close();
            exit(1);
        }
        printf("Loading norm map from file %s ...\n", normfile.c_str());
        infile.read((char *) normtextbuf, sizeof(float) * normMapRes * normMapRes);
        printf("%s Loaded ...\n", normfile.c_str());
        infile.close();
        return;
    }else{
        float * tempnorm = new float[normPointSize * normMapRes/2];
        infile.close();
        printf("Normalization file not found, generating ...\n");
        //otherwise, create the file
        float r0,r1;
        float vr0[3];
        //float xy[2];
        float dtheta, theta;
        float norm = 0;

        int _ck = 0;

        for(int i = 0; i < normPointSize; i ++){
            for(int j = 0; j < normMapRes / 2.0; j ++){

                if(_ck % (normPointSize*normMapRes / 2 / 20) == 0 ){
                    printf("%d%% ", _ck * 100/(normPointSize*normMapRes / 2)+1);
                    cout.flush();
                }
                _ck++;
                r0 = ((float) j / ((float) (normMapRes/2.0)));
                r1 = ((float) i / ((float) (normMapRes)));
                if(i == 0){
                    r1 = (1.0 / ((float) (normMapRes)));
                }

                float theta1 = acos((-1.0 + (r1+r0)*(r1+r0))/(1.0 + (r1+r0)*(r1+r0)));
                float theta2 = acos((-1.0 + (r0-r1)*(r0-r1))/(1.0 + (r0-r1)*(r0-r1)));
            
                dtheta = abs(theta1 - theta2) / 2.0;
                theta = abs(theta1 + theta2) / 2.0;
                if(((r1+r0)) * ((r0-r1)) < 0.0){
                    dtheta =abs(2*PI - theta1 - theta2) / 2.0;
                    theta = abs(2*PI - theta2 + theta1) / 2.0;
                    }
            
                vr0[0] = sin(theta);
                vr0[1] = 0.0;
                vr0[2] = cos(theta);


                float coorx = j;
                float coory = 0;
                float x = 0;
                float y = 0;
                float newsize = i;
                if(newsize < 1){
                    newsize = 1;
                }
                norm = 0;
                int n = 0;
                for(x = 0.0; x < newsize; x++){
                    for(y = 0.0; y < newsize; y++){
                        float px = (x+0.5)/newsize;
                        float py = (y+0.5)/newsize;
                        px = 2.0*(px-0.5); // -1...1
                        py = 2.0*(py-0.5);
                        float u = px * px + py * py;
                        if (u > 1.0)
                            continue;
                        float xyr[2];
                        xyr[0] = (px * (newsize / 2.0) + coorx) / (normMapRes / 2.0);
                        xyr[1] = (py * (newsize / 2.0) + coory) / (normMapRes / 2.0);
                        //xyr[0] = (coorx) / (normMapRes / 2.0);
                        //xyr[1] = (coory) / (normMapRes / 2.0);

                        float pr2 = xyr[0] * xyr[0] + xyr[1] * xyr[1];
                        norm += 4.0/(1.0+pr2)/(1.0+pr2) * projprofile(xyr, vr0, 1.0, dtheta);
                        n++;
                        //float t =  4.0/(1.0+pr2)/(1.0+pr2) * projprofile(xyr, vr0, 1.0, dtheta);
                        /*if(t == 0.0) {
                            printf("%f %f %f %f %f %f\n", xyr[0], xyr[1], vr0[0], vr0[1], vr0[2], dtheta);
                        }*/
                    }
                }

                tempnorm[i * normMapRes/2 + j]=norm;

            }

        }
        
        
        for(int i = 0; i < normMapRes; i ++){
            for(int j = 0; j < normMapRes; j ++){
                int row = i * normPointSize / normMapRes;
                int col = j / 2;
                normtextbuf[i * normMapRes + j] = tempnorm[row * normMapRes/2 + col];
                //if(i == 0){
                //    normtextbuf[i * normMapRes + j] = 4.0;
                //}
            }
        }
                
        printf("Saving to file %s ... \n", normfile.c_str());
        ofstream oufile;
        oufile.open(normfile.c_str(), ios::out | ios::binary);
        oufile.write((char * ) &normMapRes, sizeof(int));
        oufile.write((char * ) normtextbuf, normMapRes * normMapRes * sizeof(float));
        oufile.close();
        delete tempnorm;
    }
     
}

//18885004802
