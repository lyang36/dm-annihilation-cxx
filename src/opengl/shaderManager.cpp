#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <sys/types.h>

#include <sys/types.h>
#include <stdlib.h>
#ifdef __APPLE__
#include <glew.h>
#include <GLUT/glut.h> // darwin uses glut.h rather than GL/glut.h
#include <unistd.h>
#endif

//#ifdef _WIN32 || _WIN64
#if defined(_WIN32) || defined(_WIN64)
#include <GL/glew.h>
#include <GL/glut.h>
#include <io.h>
#endif

#ifdef __linux__
#include <GL/glew.h>
#include <GL/glut.h>
#include <unistd.h>
#endif
#include "shaderManager.h"

void shaderObj::begin(){
    
    if (progmObj == 0) {
        printf("Shader incorrect!!\n");
        exit(1);
    };
    
    glUseProgram(progmObj);
    
}


void shaderObj::end(){

    glUseProgram(0);

}

char * shaderManager::textFileRead(const char * fn){
    
    
    FILE *fp;
    char * content = NULL;
    
    int f, count;
    
    f = open(fn, O_RDONLY);
    
    count = lseek(f, 0, SEEK_END);
    
    close(f);
    
    if (fn != NULL) {
        fp = fopen(fn,"rt");
        
        if (fp != NULL) {
            
            
            if (count > 0) {
                content = (char *)malloc(sizeof(char) * (count+1));
                count = fread(content,sizeof(char),count,fp);
                content[count] = '\0';
            }
            fclose(fp);
        }else{
            printf("no shader file: %s\n", fn);
            exit(1);
        }
    }
    return content;
}

void shaderManager::printLog(GLhandleARB obj){
    int infologLength = 0;
    int charsWritten  = 0;
    char *infoLog;
    
    glGetObjectParameterivARB(obj, GL_OBJECT_INFO_LOG_LENGTH_ARB,
                              &infologLength);
    if (infologLength > 1)
    {
        infoLog = (char *)malloc(infologLength);
        glGetInfoLogARB(obj, infologLength, &charsWritten, infoLog);
        printf("%s\n",infoLog);
        free(infoLog);
        exit(1);
    }
}

shaderObj * shaderManager::loadShaderFile(const char * vfile, const char * ffile){
    char *vs,*fs;
    vs = textFileRead(vfile);
	fs = textFileRead(ffile);
    //const char* vvs = vs;
    //const char* ffs = fs;
    shaderObj * temp = loadShader(vs, fs);
    free(vs);free(fs);
    return temp;
}

shaderObj * shaderManager::loadShader(const char *vf, const char * ff){
    
    GLhandleARB v,f,p; 
    shaderObj * sobj = new shaderObj();
    v = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
	f = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);
    
    glShaderSourceARB(v, 1, &vf, NULL);
	glShaderSourceARB(f, 1, &ff, NULL);
    

    
	glCompileShaderARB(v);
	glCompileShaderARB(f);
    
	printLog(v);
	printLog(f);
    
    p = glCreateProgramObjectARB();
    
	glAttachObjectARB(p,f);
	glAttachObjectARB(p,v);
    
	glLinkProgramARB(p);
    
    sobj->vertexShader = v;
    sobj->fragmentShader = f;
    sobj->progmObj = p;
    return sobj; 
    
}
