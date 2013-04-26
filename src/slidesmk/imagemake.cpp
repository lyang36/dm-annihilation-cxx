/**
 * Make a slide of any halo
 */
#include <sstream>
#include <string>
#include <cmath>
#include <cstring>


#include <rpc/types.h>
#include <rpc/xdr.h>

#include "datareader.h"
#include "datatypes.h"
#include "../point.h"
#include "../halocore/halocore.h"


#ifndef PI
#define PI 3.141592653589793238462643383279
#endif


using namespace std;

void getJetColor(float value, float &r, float &g, float &b) {
    float fourValue = 4 * value;
    r = min(fourValue - 1.5, -fourValue + 4.5);
    g = min(fourValue - 0.5, -fourValue + 3.5);
    b = min(fourValue + 0.5, -fourValue + 2.5);
    if(r > 1.0) r = 1.0;
    if(r < 0.0) r = 0.0;
    if(g > 1.0) g = 1.0;
    if(g < 0.0) g = 0.0;
    if(b > 1.0) b = 1.0;
    if(b < 0.0) b = 0.0;
}

void getcolorImge(float *value, float * colorimg, int wsize, int hsize,
                  float immin = 0.0, float immax = -1.0){
    float r, g, b;
    
    if(immax < immin){
        float max = 0.0, min = 1.0e20;
        for(int i = 0; i < wsize * hsize; i++){
            if(max < value[i]){
                max = value[i];
            }
            if(min > value[i] && value[i] > 0.0){
                min = value[i];
            }
            //test
            //if(value[i] == 0.0){
            //    printf("%d %d %d", i, i%windowSize, i/windowSize);
            //}
        }
        immin = min;
        immax = max;
    }

	if(immin == immax){
		immin = immax / 1.0e5;
	}
    
    
	printf("Min = %f  Max = %f\n", immin, immax);
    
    float x = log(immax) - log(immin);
    for(int i = 0; i < wsize * hsize; i++){
        float v = (log(value[i]) - log(immin)) / x;
        getJetColor(v, r, g, b);
        colorimg[3 * i] = r;
        colorimg[3 * i + 1] = g;
        colorimg[3 * i + 2] = b;
    }
}


void savePic(string picfile, float * img, int wsize, int hsize){
    printf("Save pic to file: %s ...\n", picfile.c_str());
    //save to file
    ofstream fScreenshot(picfile.c_str(), ios::out | ios::binary);
    if (!fScreenshot.good()) {
        printf("Open pic file error.\n");
    }else{
        unsigned char TGAheader[12]={0,0,2,0,0,0,0,0,0,0,0,0};
        unsigned char header[6] = {wsize%256,wsize/256,
            hsize%256,hsize/256,24,0};
        //convert to BGR format
        unsigned char temp;
        unsigned int i = 0;
        unsigned char * pixels = new unsigned char[wsize * hsize * 3];
        while (i < wsize*hsize*3)
        {
            //temp = pixels[i];       //grab blue
            pixels[i] = (int)(img[i+2] * 256);//assign red to blue
            pixels[i+2] = (int)(img[i] * 256);//temp;     //assign blue to red
            pixels[i+1] = (int)(img[i+1] * 256);
            i += 3;     //skip to next blue byte
        }
        fScreenshot.write((char *)TGAheader, 12);
        fScreenshot.write((char *)header, 6);
        fScreenshot.write((char *)pixels, wsize*hsize*3);
        fScreenshot.close();
        printf("Saving ok.\n");
        delete pixels;
    }
}


int main(int argc, const char **argv){
    //printf("%d\n", argc);
    if(!(argc == 2 || argc == 5 || argc == 8)){
        printf("Usage: %s <basename> [min = min] [max = max]\n", argv[0]);
        exit(1);
    }
    
    string basename;
    float imgmin = 0.0;
    float imgmax = -1.0;
    int imagesize;
    
    
    stringstream ss;
    ss << argv[1];
    ss >> basename;
    if(argc >= 5){
        if((strstr (argv[2], "min") != NULL) &&
           (strstr (argv[3], "=") != NULL)){
            stringstream ss1;
            ss1 << argv[4];
            ss1 >> imgmin;
        }else if((strstr (argv[2], "max") != NULL) &&
                 (strstr (argv[3], "=") != NULL)){
            stringstream ss1;
            ss1 << argv[4];
            ss1 >> imgmax;
        }
        else{
            printf("Usage: %s <basename> [min = min] [max = max]\n", argv[0]);
            exit(1);
        }
    }
    
    if(argc >= 8){
        if((strstr (argv[5], "min") != NULL) &&
           (strstr (argv[6], "=") != NULL)){
            stringstream ss1;
            ss1 << argv[7];
            ss1 >> imgmin;
        }else if((strstr (argv[5], "max") != NULL) &&
                 (strstr (argv[6], "=") != NULL)){
            stringstream ss1;
            ss1 << argv[7];
            ss1 >> imgmax;
        }
        else{
            printf("Usage: %s <basename> [min = min] [max = max]\n", argv[0]);
            exit(1);
        }
    }
    
    
    
    fprintf(stderr, "Base Name: %s\n", basename.c_str());

    string ifiledens = basename + ".dens";
    string ifilevdis = basename + ".vdis";
    
    fstream infile(ifiledens.c_str(), ios::in | ios::binary);
    if(infile.good()){
        infile.read((char *) & imagesize, sizeof(int));
    }else{
        printf("File Error!!\n");
        exit(1);
    }
    
    int numbins = imagesize * imagesize;
    double * databins = new double[numbins];
    double * vdisp = new double[numbins];
    
    if(infile.good()){
        infile.read((char *) databins, sizeof(double) * numbins);
    }else{
        printf("File Error!!\n");
        exit(1);
    }
    
    infile.close();
    
    fstream infile1(ifilevdis.c_str(), ios::in | ios::binary);
    if(infile1.good()){
        infile1.read((char *) & imagesize, sizeof(int));
    }else{
        printf("File Error!!\n");
        exit(1);
    }
    
    if(imagesize * imagesize != numbins){
        printf("Error: file not match!\n");
        exit(1);
    }
    if(infile1.good()){
        infile1.read((char *) vdisp, sizeof(double) * numbins);
    }else{
        printf("File Error!!\n");
        exit(1);
    }
    infile1.close();
    
    fprintf(stderr, "Image Size: %d\n", imagesize);
    fprintf(stderr, "Min = %f Max = %f\n", imgmin, imgmax);
    
    //produce the image
     float *rhosq = new float[imagesize * imagesize];
     float *rhosqv = new float[imagesize * imagesize];
     float *rhosqv2 = new float[imagesize * imagesize];
     float *imag = new float[imagesize * imagesize * 3];
     
     //produce rho-square image
    float normfacsq = 0.0;
    float normfacsqv = 0.0;
    float normfacsqv2 = 0.0;
    for(int i = 0; i < numbins; i++){
        rhosq[i] = databins[i] * databins[i];
        if(vdisp[i] != 0){
            rhosqv[i] = rhosq[i] / vdisp[i];
            rhosqv2[i] = rhosqv[i] / vdisp[i];
        }
        
        if(rhosq[i] > normfacsq){
            normfacsq = rhosq[i];
        }
        
        if(rhosqv[i] > normfacsqv){
            normfacsqv = rhosqv[i];
        }
        
        if(rhosqv2[i] > normfacsqv2){
            normfacsqv2 = rhosqv2[i];
        }
    }
    for(int i = 0; i < numbins; i++){
        rhosq[i] /= normfacsq;
        rhosqv[i] /= normfacsqv;
        rhosqv2[i] /= normfacsqv2;
    }
    
    //printf("%f, %f, %f\n", normfacsq, normfacsqv, normfacsqv2);
    
    printf("Get Rho-Sqr image\n");
    getcolorImge(rhosq, imag, imagesize, imagesize, imgmin, imgmax);
    savePic(basename + "_rsqr.tga", imag, imagesize, imagesize);
    
    delete rhosq, rhosqv, rhosqv2, imag;
    delete databins;
    delete vdisp;
    
    return 0;
}
