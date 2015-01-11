// measure the clumpniness in a shell
// get the result of <rho^2>/<rho>^2 of each shell

#include <sstream>
#include <string>
#include <cmath>
#include <cstring>


#include <rpc/types.h>
#include <rpc/xdr.h>

#include "datareader.h"
#include "datatypes.h"
#include "../halocore/halocore.h"


#ifndef PI
#define PI 3.141592653589793238462643383279
#endif


using namespace std;

int main(int argc, const char **argv){

    int numbins = 100;
    double * rhobins;
    double * rhosqbins;
    double * countbins;
    
    double dr, x, y, z;

    dr = 4.0; //kpc
    x = 3.612030e-01;
    y = 2.106580e-01;
    z = 6.877260e-03;
    
    double simsize = 40000.0; //simulation size
    
    string filename = "/data1/data/data_float_all.bin";
    //string filename = "/Users/lyang/data/data_samp_1e6.bin";

    
    rhobins = new double[numbins];
    rhosqbins = new double[numbins];
    countbins = new double[numbins];
    
    for(int i = 0; i < numbins; i++){
        rhobins[i] = 0;
        rhosqbins[i] = 0;
        countbins[i] = 0;
    }
    
    //fprintf(stderr, "Filename: %s\n", filename.c_str());
    //fprintf(stderr, "Bins: %d\n", numbins);
    //fprintf(stderr, "x y z: %f %f %f (40 Mpc/h)\n", x, y, z);
    
    DataReader reader(filename);
    reader.open();

    int cts = 0;
    while(reader.hasNext()){
        DMParticle * ps = reader.getBuf();
        
        for(int i = 0; i < reader.getMemparts(); i++){
            DMParticle & part = ps[i];
            
            double r = sqrt((part.posx - x) *(part.posx - x) +
                            (part.posy - y) *(part.posy - y) +
                            (part.posz - z) *(part.posz - z));
            r *= simsize;
            int inds = floor(r / dr);
            if((inds < numbins) && (part.dens > 0)){
                rhobins[inds] += part.dens;
                rhosqbins[inds] += part.dens * part.dens;
                countbins[inds] ++;
            }
            cts ++;
        }
        reader.loadBuffer();
    }
    reader.close();
    
    for(int i = 0; i < numbins; i++){
        if(countbins[i] != 0){
            rhobins[i] /= countbins[i];
            rhosqbins[i] /= countbins[i];
        }else{
            rhobins[i] = 0.0;
            rhosqbins[i] = 0.0;
        }
        
        printf("%f %f %f\n", i * dr, rhobins[i], rhosqbins[i]);
    }
    
    delete[] rhobins;
    delete[] rhosqbins;
    delete[] countbins;
    
    return 0;
}
