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


#include "VL2_CONSTANT.h"
#include "datareader.h"
#include "datatypes.h"
#include "../halocore/halocore.h"


#ifndef PI
#define PI 3.141592653589793238462643383279
#endif



using namespace std;

void savetoFits(double * map, int npix, string fn);


void printUsage(const char * s){ 
    printf("Usage:\n"
           "%s\n"
           "-d <datafile>\n"
           "-r <radius (kpc)>\n"
           "-b <boxsize (kpc)>\n"
           "-o <outputfile prefix>\n"
           "-c <x>, <y>, <z> (code units, center coordinates)\n"
           "-m <maskfile>\n"
           , s); 
}




int main(int argc, const char **argv){

    
    
    double shellr;               //shell radius
    double boxsize, angularsize; //boxsize
    int nside, npix;             //nside, determined by the boxsize and r
    double x, y, z;              //center coordinates

    double * rhosq_map;
    double * rho_map;
    double * rhosq_v_map;        //somfeld 1/v map
    double * rhosq_v2_map;       //somfeld 1/v^2 map
    double * rho_v_map;          //rho1/vdisp map
    double * rho_v2_map;         //rho1/vdisp^2 map
    double * count_map;

    bool isMask = false;
    string maskfile = "";

    shellr = 8.0; //kpc
    boxsize = 0.5; //kpc
    x = 3.612030e-01;
    y = 2.106580e-01;
    z = 6.877260e-03;

    string filename = "/data1/data/data_float_all.bin";
    string outfile_rhomap = "./rhomap.fits";
    string outfile_rhosqmap = "./rhosqmap.fits";
    string outfile_rhosq_v_map = "./rhosqmap_v.fits";
    string outfile_rhosq_v2_map = "./rhosqmap_v2.fits";
    string outfile_rho_v2_map = "./rhomap_v2.fits";
    string outfile_rho_v_map = "./rhomap_v.fits";





    // read from the command line
    int m=1;
    while (m<argc)
    {
        string arg = argv[m];
        //printf("ok!\n");
        stringstream ss, ss1, ss2;
        if (arg == "-d") {
            filename = argv[m + 1]; 
            m+=1;
        }else if (arg == "-r") {
            ss << argv[m+1];
            ss >> shellr;
            m+=1;
        }else if (arg == "-b") {
            ss << argv[m+1];
            ss >> boxsize;
            m+=1;
        }else if (arg == "-c") {
            ss << argv[m+1];
            ss >> x;
            m++;

            ss1 << argv[m+1];
            ss1 >> y;
            m++;

            ss2 << argv[m+1];
            ss2 >> z;
            m++;
        }else if (arg == "-m") {
            isMask = true;
            maskfile = argv[m+1];
            m+=1;
        }else if (arg == "-o"){
            string tmpst = "";
            outfile_rhomap = tmpst + argv[m+1] + "_rho.fits";

            outfile_rhosqmap = tmpst + argv[m+1] + "_rhosq.fits";

            outfile_rho_v_map = tmpst + argv[m+1] + "_rho_v.fits";
            outfile_rho_v2_map = tmpst + argv[m+1] + "_rho_v2.fits";
            
            outfile_rhosq_v_map = tmpst + argv[m+1] + "_rhosq_v.fits";
            outfile_rhosq_v2_map = tmpst + argv[m+1] + "_rhosq_v2.fits";
            m++;
        }
        else{
            printUsage(argv[0]);
            exit(1);
        }
        m++;
    }
 
    printf("Measure Clumpiness. Output file int<rho>,\n"
            "and int<rho^2>, in units of M_sun/Mpc**3\n");
    printf("Datafile: %s\n", filename.c_str());
    printf("Shell Radius: %f (Kpc)\n", shellr);
    printf("Boxsize: %f (Kpc)\n", boxsize);
    printf("Center: %f %f %f (code)\n", x, y, z);
    printf("Rho Map: %s\n", outfile_rhomap.c_str());
    printf("Rho Sq Map: %s\n", outfile_rhosqmap.c_str());
    printf("Rho <v> Map: %s\n", outfile_rho_v_map.c_str());
    printf("Rho <v2> Map: %s\n", outfile_rho_v2_map.c_str());
    printf("Rho sq <v> Map: %s\n", outfile_rhosq_v_map.c_str());
    printf("Rhosq <v2> Map: %s\n", outfile_rhosq_v2_map.c_str());
    // start computing
    angularsize = boxsize / shellr;
    npix = PI * 4.0 / angularsize / angularsize;
    nside = pow(2, ceil(log(sqrt(npix / 12.0)) / log(2.0)));
    if(nside == 0){
        nside = 1;
    }
    npix = 12 * nside * nside;
    printf("Nside = %d, Npix = %d.\n", nside, npix);

    double simsize = BOXSIZE_KPC; //simulation size, kpc
    
    //string filename = "/Users/lyang/data/data_samp_1e6.bin";
    
       
    Healpix_Base base(nside, RING, SET_NSIDE);
    rhosq_map = (double *) malloc(sizeof(double) * npix);
    rhosq_v_map = (double *) malloc(sizeof(double) * npix);
    rhosq_v2_map = (double *) malloc(sizeof(double) * npix);

    rho_map = (double *) malloc(sizeof(double) * npix);
    rho_v_map = (double *) malloc(sizeof(double) * npix);
    rho_v2_map = (double *) malloc(sizeof(double) * npix);

    count_map = (double *) malloc(sizeof(double) * npix);

    for(int i = 0; i < npix; i++){
        rhosq_map[i] = 0;
        rhosq_v_map[i] = 0;
        rhosq_v2_map[i] = 0;

        rho_map[i] = 0;
        rho_v_map[i] = 0;
        rho_v2_map[i] = 0;
        
        count_map[i] = 0;
    }
    
    //fprintf(stderr, "Filename: %s\n", filename.c_str());
    //fprintf(stderr, "Bins: %d\n", numbins);
    //fprintf(stderr, "x y z: %f %f %f (40 Mpc/h)\n", x, y, z);
    
    DataReader reader(filename);
    reader.open();

    double sd, su, rd, ru;         // the area of the bottom and top
    rd = shellr + boxsize / 2.0;
    ru = shellr - boxsize / 2.0;
    sd = 4.0 * PI / ((double) npix) * pow(rd, 2);
    su = 4.0 * PI / ((double) npix) * pow(ru, 2);
    double volumesize = 1.0 / 3.0 * (rd * sd  - ru * su);
    printf("Volume size: %f (kpc^3)\n", volumesize);

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
            r *= simsize;                           // in Kpc
            
            if((r <= (boxsize / 2.0 + shellr)) && 
                    (r > (shellr - boxsize / 2.0)) ){
                // calculate the healpix id 
    
                vec3 vec(px, py, pz);
                vec.Normalize();
                pointing p(vec);
    
                int pix = base.ang2pix(p);

                double intrho =  part.mass * MASS_UNIT / volumesize * 1.0e9;       // M_sun / Mpc^3 
                double intrhosq =  part.dens * part.mass * 
                     RHO_C * H * H * MASS_UNIT/ volumesize * 1.0e9; // (M_sun / Mpc^3)^2

                rho_map[pix] += intrho;                          
                
                rho_v_map[pix] += intrho / part.sigmav;
                rho_v2_map[pix] += intrho / part.sigmav / part.sigmav;
                
                rhosq_map[pix] += intrhosq;
                rhosq_v_map[pix] += intrhosq / part.sigmav;
                rhosq_v2_map[pix] += intrhosq / part.sigmav / part.sigmav;
                


                count_map[pix] += 1.0;
            }

            cts ++;
        }
        reader.loadBuffer();
    }
    reader.close();
    
    //save to file


    savetoFits(rho_map, npix, outfile_rhomap);
    savetoFits(rhosq_map, npix, outfile_rhosqmap);

    savetoFits(rho_v_map, npix, outfile_rho_v_map);
    savetoFits(rho_v2_map, npix, outfile_rho_v2_map);

    savetoFits(rhosq_v_map, npix, outfile_rhosq_v_map);
    savetoFits(rhosq_v2_map, npix, outfile_rhosq_v2_map);
    
    free(rho_map);
    free(rho_v_map);
    free(rho_v2_map);

    free(rhosq_map);
    free(rhosq_v_map);
    free(rhosq_v2_map);

    free(count_map);
    
    return 0;
}

void savetoFits(double * map, int npix, string fn){
   
    ifstream ifile(fn.c_str());
    if (ifile) {
    // The file exists, and is open for input
        cout << "File Exists! Owerite!" << endl;
        remove(fn.c_str());
    }
    ifile.close();


    arr<double> map_rho (map, npix);
    Healpix_Map<double> outputmap_rho (map_rho, RING);
    cout << "Creating fits" << fn <<"!..." << endl;
    fitshandle fitswriter;
    fitswriter.create(fn);
    write_Healpix_map_to_fits<double>(fitswriter, outputmap_rho, PLANCK_FLOAT64 );
    fitswriter.close();
}
