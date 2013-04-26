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

int main(int argc, const char **argv){
    //printf("%d\n", argc);
    if(!(argc == 11 || argc == 12 || argc == 13)){
        printf("Usage: %s <fitsfilename>, <outbase>, <radius>, <imagesize>, <x>, <y>, <z> <theta, phi, psi> [maskfile] [-core]\n", argv[0]);
        exit(1);
    }
    
    //radius is for the size of the half plane
    double radius;
    double dz;
    
    //the image center is the origin
    int imagesize;
    
    //the new coordinates
    double theta, phi, psi;
    
    //data and countbins
    double * databins;
    double * vdisp;
    double * countbins;

    //the halo center point
    double x, y, z;
    
    string filename;
    string outbase;
    string maskfile;
    
    bool isMask;
    bool isCore;
    
    isCore = false;
    if(argc == 11){
        isMask = false;
        isCore = false;
    }else if(argc == 12){
        //printf("%s\n", argv[7]);
        if(strstr (argv[11], "-core") != NULL){
            isMask = false;
            isCore = true;
        }else{
            isMask = true;
            isCore = false;
        }
    }else if(argc == 13){
        isMask = true;
        isCore = true;
    }
    
    stringstream ss, ss1, ss2,
                ss3, ss4, ss5,
                ss6, ss7, ss8,
                ss9, ss10, ss11;
    ss << argv[1];
    ss >> filename;
    ss11 << argv[2];
    ss11 >> outbase;
    ss1 << argv[3];
    ss1 >> radius;
    ss2 << argv[4];
    ss2 >> imagesize;
    ss3 << argv[5];
    ss3 >> x;
    ss4 << argv[6];
    ss4 >> y;
    ss5 << argv[7];
    ss5 >> z;
    ss6 << argv[8];
    ss6 >> theta;
    ss7 << argv[9];
    ss7 >> phi;
    ss8 << argv[10];
    ss8 >> psi;
    if(isMask){
        ss10 << argv[11];
        ss10 >> maskfile;
    }
    
    int numbins = imagesize * imagesize;
    databins = new double[numbins];
    countbins = new double[numbins];
    vdisp = new double[numbins];

    for(int i = 0; i < numbins; i++){
        databins[i] = 0;
        countbins[i] = 0;
    }
    
    fprintf(stderr, "Filename: %s\n", filename.c_str());
    fprintf(stderr, "Output Base: %s\n", outbase.c_str());
    fprintf(stderr, "Imagesize: %d\n", imagesize);
    fprintf(stderr, "Radius: %f (Mpc)\n", radius);
    fprintf(stderr, "x y z: %f %f %f (40 Mpc)\n", x, y, z);
    fprintf(stderr, "theta phi psi: %f %f %f\n", theta, phi, psi);
    if(isMask){
        fprintf(stderr, "Mask: %s\n", maskfile.c_str());
    }
    radius /= 40.0;
    dz = radius / 10.0;
    
    DataReader reader(filename, isMask, maskfile);
    reader.open();
    
    //static HaloCore halocore;
    
    double minsigmaav = 1.0e99;
    double maxsigmaav = 0.0;
    
    int cts = 0;
    
    Point e1, e2, e3, e1p, e2p;
    
    e1.x = cos(theta) * cos(phi);
    e1.y = cos(theta) * sin(phi);
    e1.z = -sin(theta);
    
    e2.x = -sin(phi);
    e2.y = cos(phi);
    e2.z = 0;
    
    e3.x = sin(theta) * cos(phi);
    e3.y = sin(theta) * sin(phi);
    e3.z = cos(theta);
    
    e1p = e1 * cos(psi) + e2 * sin(psi);
    e2p = e1 * (- sin(psi)) + e2 * cos(psi);
    
    while(reader.hasNext()){
        DMParticle * ps = reader.getBuf();
        
        for(int i = 0; i < reader.getMemparts(); i++){
            DMParticle & part = ps[i];
            Point vec;
            vec.x = part.posx - x;
            vec.y = part.posy - y;
            vec.z = part.posz - z;
            
            double ix = vec.dot(e1p);
            double iy = vec.dot(e2p);
            double iz = vec.dot(e3);
            

            part.hsmooth *= 3;
            //printf("%f %f %f %f %f\n", part.posx, part.posy, part.posz, iz, dz);
            if((abs(iz) <= part.hsmooth)){
                //printf("%f %f\n", iz, dz);
                
                //calculate the configuration on the plane
                double imr = sqrt(part.hsmooth * part.hsmooth - iz * iz);
                double xmin, ymin, xmax, ymax;
                int ixmin, iymin, ixmax, iymax;
                xmin = (ix - imr);
                ymin = (iy - imr);
                xmax = (ix + imr);
                ymax = (iy + imr);
                ixmin = floor((xmin / radius + 1) * imagesize / 2);
                iymin = floor((ymin / radius + 1) * imagesize / 2);
                ixmax = ceil((xmax / radius + 1) * imagesize / 2);
                iymax = ceil((ymax / radius + 1) * imagesize / 2);
                
                //printf("%d %d %d %d\n", ixmin, ixmax, iymin, iymax);
                
                if((((xmin > - radius) && (xmin < radius))
                   || ((xmax > - radius) && (xmax < radius)))
                   && (((ymax > - radius) && (ymax < radius))
                   || ((ymin > - radius) && (ymin < radius)))){
                    for(int imi = ixmin ;
                        imi < ixmax;
                        imi++){
                        for(int imj = iymin;
                            imj < iymax;
                            imj ++){
                            if((imi < imagesize && imi >= 0)
                               && (imj >=0 && imj < imagesize)){
                                //calculate the radius
                                double cx, cy;
                                cx = (imi - imagesize / 2) * radius / (imagesize / 2);
                                cy = (imj - imagesize / 2) * radius / (imagesize / 2);
                                //printf("%f %f %f %f %f %f\n", cx, cy, xmin, xmax, ymin, ymax);
                                if((cx - ix) * (cx - ix) +
                                   (cy - iy) * (cy - iy)
                                    < imr * imr){
                                    databins[imi + imj * imagesize] += part.dens;
                                    vdisp[imi + imj * imagesize] += part.sigmav;
                                    countbins[imi + imj * imagesize] ++;
                                }
                            }
                        }
                        
                    }
                }
            }
            cts ++;
        }
        reader.loadBuffer();
        //printf("%d\n", cts);
    }
    //printf("\n");
    
   //fprintf(stderr, "Sigma V range: [%f %f]\n", minsigmaav, maxsigmaav);
    
    reader.close();
    
    for(int i = 0; i < numbins; i++){
        if(countbins[i] != 0){
            databins[i] /= countbins[i];
            vdisp[i] /= countbins[i];
        }else{
            databins[i] = 0.0;
            vdisp[i] = 0.0;
        }
        
        //printf("%f %f\n", i * dr * 40, databins[i]);
    }
    
    string ofiledens = outbase + ".dens";
    string ofilevdis = outbase + ".vdis";
    
    fstream outfile(ofiledens.c_str(), ios::out | ios::binary);
    outfile.write((char *) & imagesize, sizeof(int));
    outfile.write((char *) databins, sizeof(double) * numbins);
    outfile.close();
    
    fstream outfile1(ofilevdis.c_str(), ios::out | ios::binary);
    outfile1.write((char *) & imagesize, sizeof(int));
    outfile1.write((char *) vdisp, sizeof(double) * numbins);
    outfile1.close();
    
    /* //produce the image
    float *values = new float[imagesize * imagesize];
    float *imag = new float[imagesize * imagesize * 3];
    float imgmin = 1.0e-8;
    float imgmax = 1.0;
    
    //produce rho-square image
    for(int i = 0; i < numbins; i++){
        if(countbins[i] != 0){
            databins[i] /= countbins[i];
            vdisp[i] /= countbins[i];
        }else{
            databins[i] = 0.0;
            vdisp[i] = 0.0;
        }
    }*/
    delete databins;
    delete countbins;
    delete vdisp;
     
    return 0;
}
