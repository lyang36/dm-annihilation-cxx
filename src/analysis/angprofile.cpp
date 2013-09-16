//fits and HEALPIX
#include <sstream>
#include <string>
#include <cmath>

#include <healpix_base.h>
#include <healpix_map.h>
#include <arr.h>
#include <fitsio.h>
#include <fitshandle.h>
#include <healpix_map_fitsio.h>

#ifndef PI
#define PI 3.141592653589793238462643383279
#endif

using namespace std;

int main(int argc, const char **argv){
    if(argc != 6){
        printf("Usage: %s <fitsfilename>, <num of bins>, <angular radius (degree)>, <theta>, <phi>\n", argv[0]);
        exit(1);
    }
    
    int numbins;
    double anglerad;
    double theta, phi;
    double * databins;
    //double * angbins;
    double * countbins;
    double dtheta;
    string filename;
    
    stringstream ss, ss1, ss2, ss3, ss4;
    ss << argv[1];
    ss >> filename;
    ss1 << argv[2];
    ss1 >> numbins;
    ss2 << argv[3];
    ss2 >> anglerad;
    ss3 << argv[4];
    ss3 >> theta;
    ss4 << argv[5];
    ss4 >> phi;
    
    databins = new double[numbins];
    //angbins = new double[numbins + 1];
    countbins = new double[numbins];
    
    fprintf(stderr, "Filename: %s\n", filename.c_str());
    fprintf(stderr, "Bins: %d\n", numbins);
    fprintf(stderr, "Angular Radius: %f\n", anglerad);
    fprintf(stderr, "Theta: %f\n", theta);
    fprintf(stderr, "Phi: %f\n", phi);
    
    anglerad = PI / 180.0 * anglerad;
    theta = PI / 180.0 * theta;
    phi = PI / 180.0 * phi;
    
    for(int i = 0; i < numbins; i++){
        databins[i] = 0;
        //angbins[i] = (double) i * anglerad / (double) numbins;
        countbins[i] = 0;
    }
    
    dtheta = anglerad / (double) numbins;
    
    Healpix_Map<float> inputmap;
    //fitshandle fitsreader;
    //fitsreader.open(filename);
    read_Healpix_map_from_fits<float>(filename, inputmap, 1,2);
    //fitsreader.close();
    
    vec3 centvec;
    centvec.x = sin(theta)*cos(phi);
    centvec.y = sin(theta)*sin(phi);
    centvec.z = cos(theta);
    
    for(int i = 0; i < inputmap.Npix(); i++){
        vec3 cvec = inputmap.pix2vec(i);
        double angle = acos(cvec.x * centvec.x + cvec.y * centvec.y + cvec.z * centvec.z);
        if(angle < anglerad){
            int ind = angle / dtheta;
            databins[ind] += inputmap[i];
            countbins[ind] += 1;
        }
    }
    
    for(int i = 0; i < numbins; i++){
        if(countbins[i] != 0){
            databins[i] /= countbins[i];
        }else{
            databins[i] = 0.0;
        }
        
        printf("%f %f\n", i * dtheta, databins[i]);
    }
    
    delete databins;
    delete countbins;
}
