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

// TODO
void buffer::readTex(void* tex){
    
}

// TODO
void buffer::saveImg(string filename){
    
}
