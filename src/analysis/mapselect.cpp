//fits and HEALPIX
#include <sstream>
#include <fstream>
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
    if(argc != 7){
        printf("Usage: %s <fitsfilename>, <resolution>, <angular radius>, <theta>, <phi>, <outputfile>\n", argv[0]);
        exit(1);
    }
    
    int res;
    double anglerad;
    double theta, phi;
    double * databins;
    //double * angbins;
    double * countbins;
    //double dtheta;
    string filename;
    string outputfile;
    
    stringstream ss, ss1, ss2, ss3, ss4, ss5;
    ss << argv[1];
    ss >> filename;
    ss1 << argv[2];
    ss1 >> res;
    ss2 << argv[3];
    ss2 >> anglerad;
    ss3 << argv[4];
    ss3 >> theta;
    ss4 << argv[5];
    ss4 >> phi;
    ss5 << argv[6];
    ss5 >> outputfile;
    
    databins = new double[res * res];
    //angbins = new double[numbins + 1];
    countbins = new double[res * res];
    
    fprintf(stderr, "Filename: %s\n", filename.c_str());
    fprintf(stderr, "Resolution: %d\n", res);
    fprintf(stderr, "Angular Radius: %f\n", anglerad);
    fprintf(stderr, "Theta: %f\n", theta);
    fprintf(stderr, "Phi: %f\n", phi);
    
    anglerad = PI / 180.0 * anglerad;
    theta = PI / 180.0 * theta;
    phi = PI / 180.0 * phi;
    
    for(int i = 0; i < res * res; i++){
        databins[i] = 0;
        //angbins[i] = (double) i * anglerad / (double) numbins;
        countbins[i] = 0;
    }
    
    //dtheta = anglerad / (double) numbins;
    
    Healpix_Map<float> inputmap;
    //fitshandle fitsreader;
    //fitsreader.open(filename);
    read_Healpix_map_from_fits<float>(filename, inputmap, 1,2);
    //fitsreader.close();
    
    vec3 centvec;
    centvec.x = sin(theta)*cos(phi);
    centvec.y = sin(theta)*sin(phi);
    centvec.z = cos(theta);
    
    vec3 xvec;//e_theta
    vec3 yvec;//e_phi
    
    xvec.x = cos(theta)*cos(phi);
    xvec.y = cos(theta)*sin(phi);
    xvec.z = -sin(theta);
    
    yvec.x = -sin(phi);
    yvec.y = cos(phi);
    yvec.z = 0;
    
    /*for(int i = 0; i < inputmap.Npix(); i++){
        vec3 cvec = inputmap.pix2vec(i);
        
        double angle = acos(cvec.x * centvec.x + cvec.y * centvec.y + cvec.z * centvec.z);
        
        if(angle < anglerad){
            
            //int ind = angle / dtheta;
            double gr = sin(angle) / sin(anglerad);
            double stcp = cvec.x * xvec.x + cvec.y * xvec.y + cvec.z * xvec.z;
            double stsp = cvec.x * yvec.x + cvec.y * yvec.y + cvec.z * yvec.z;
            double gtheta;
            
            if(stcp == 0.0){
                gtheta = PI / 2.0;
            }else{
                gtheta = atan2(stsp, stcp) + PI;
            }
            double indr = gr * (res / 2);
            int indx = indr * cos(gtheta) + res / 2;
            int indy = indr * sin(gtheta) + res / 2;
            int ind = indx + indy * res;
            
            databins[ind] += inputmap[i];
            countbins[ind] += 1;
        }
    }*/
    
    
    
    for(int i = 0; i < res * res; i++){
        if(countbins[i] != 0){
            databins[i] /= countbins[i];
        }else{
            databins[i] = 0.0;
        }
    }
    

    
    fstream outfile(outputfile.c_str(), ios::out | ios::binary);
    outfile.write((char *) & res, sizeof(int));
    outfile.write((char *) databins, sizeof(double) * res * res);
    outfile.close();
    
    delete databins;
    delete countbins;
    
}