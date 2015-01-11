// measure the clumpniness distribution in a shell
// given a shell, with the healpix setting, each pixel is a box
// the clumpiness is measured inside the box
// get the result of int <rho^2> int <rho>^2 of each shell
// int <rho^2>  \approx rho m
// int <rho> \approx mass

#include <sstream>
#include <string>
#include <cmath>
#include <cstring>


#include <rpc/types.h>
#include <rpc/xdr.h>


//fits and HEALPIX
#include <healpix_base.h>
#include <healpix_map.h>
#include <arr.h>
#include <fitsio.h>
#include <fitshandle.h>
#include <healpix_map_fitsio.h>


#include "datareader.h"
#include "datatypes.h"
#include "../halocore/halocore.h"


#ifndef PI
#define PI 3.141592653589793238462643383279
#endif



using namespace std;

void printUsage(const char * s){ 
    printf("Usage:\n"
           "%s\n"
           "-d <datafile>\n"
           "-r <radius (kpc)>\n"
           "-c <x>, <y>, <z> (code units)\n"
           "-m <maskfile>\n"
           "-core\n" 
           "-log <rmin> (using log bins)\n", s); 
}




int main(int argc, const char **argv){

    int numbins = 100;
    double * rhobins;
    double * rhosqbins;
    double * countbins;
    
    
    double shellr;       //shell radius
    double boxsize, angularsize; //boxsize
    int nside, npix;      //nside, determined by the boxsize and r
    double x, y, z; //center coordinates

    double * rhosq_map;
    double * rho_map;
    double * count_map;

    shellr = 8.0; //kpc
    boxsize = 0.5; //kpc
    angularsize = boxsize / shellr;
    x = 3.612030e-01;
    y = 2.106580e-01;
    z = 6.877260e-03;
    npix = PI * 4.0 / angularsize / angularsize;
    nside = sqrt(npix / 12);
    if(nside == 0){
        nside = 1;
    }
    npix = 12 * nside * nside;
    printf("Nside = %d, Npix = %d.\n", nside, npix);

    double simsize = 40000.0; //simulation size, kpc
    
    string filename = "/data1/data/data_float_all.bin";
    //string filename = "/Users/lyang/data/data_samp_1e6.bin";
    string outfile_rhomap = "./rhomap.fits";
    string outfile_rhosqmap = "./rhosqmap.fits";
    
       
    Healpix_Base base(nside, RING, SET_NSIDE);
    rhosq_map = (double *) malloc(sizeof(double) * npix);
    rho_map = (double *) malloc(sizeof(double) * npix);
    count_map = (double *) malloc(sizeof(double) * npix);

    for(int i = 0; i < npix; i++){
        rhosq_map[i] = 0;
        rho_map[i] = 0;
        count_map[i] = 0;
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
            
            double px, py, pz;
            px = part.posx - x;
            py = part.posy - y;
            pz = part.posz - z;
            double r = sqrt( px * px + py * py + pz * pz);
            r *= simsize;
            
            if((r <= (boxsize / 2.0 + shellr)) && 
                    (r > (shellr - boxsize / 2.0)) ){
                // calculate the healpix id
           
    
    
                vec3 vec(px, py, pz);
                vec.Normalize();
                pointing p(vec);
    
                int pix = base.ang2pix(p);

                rho_map[pix] += part.mass;
                rhosq_map[pix] += part.dens * part.mass;
                count_map[pix] += 1.0;
            }

            cts ++;
        }
        reader.loadBuffer();
    }
    reader.close();
    
    //save to file
    arr<double> map_rho (rho_map, npix);
    arr<double> map_rhosq (rhosq_map, npix);
    Healpix_Map<double> outputmap_rho (map_rho, RING);
    Healpix_Map<double> outputmap_rhosq (map_rhosq, RING);
  
    {
        cout << "Creating fits, rho_map!..." << endl;
        fitshandle fitswriter;
        fitswriter.create(outfile_rhomap);
        write_Healpix_map_to_fits<double>(fitswriter, outputmap_rho, PLANCK_FLOAT64 );
        fitswriter.close();
    }
  
    {
        cout << "Creating fits, rhosq_map!..." << endl;
        fitshandle fitswriter;
        fitswriter.create(outfile_rhosqmap);
        write_Healpix_map_to_fits<double>(fitswriter, outputmap_rhosq, PLANCK_FLOAT64 );
        fitswriter.close();
    }

    free(rho_map);
    free(rhosq_map);
    free(count_map);
    
    return 0;
}
