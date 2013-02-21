#include <string>
#include <rpc/types.h>
#include <rpc/xdr.h>

//fits and HEALPIX
#include <healpix_base.h>
#include <healpix_map.h>
#include <arr.h>
#include <fitsio.h>
#include <fitshandle.h>
#include <healpix_map_fitsio.h>


#include "mapgenerator.h"

using namespace std;

MapGenerator::MapGenerator(Parameters * par, DataReader * reader){
    par_ = par;
    reader_ = reader;
    map_ = NULL;
}

MapGenerator::~MapGenerator(){
    if(map_!=NULL){
        free(map_);
    }
}

void MapGenerator::saveToFits(string filename){
    if(!isFinished()){
        printf("Map not generated yet");
        return;
    }
    
    //write a fits file
    int npix = 12 * (par_->map.Nside) * (par_->map.Nside);
    arr<double> maparr (map_, npix);
    Healpix_Map<double> outputmap (maparr, RING);
    string fits_filename = par_-> outputFileName;
    ifstream ifile(fits_filename.c_str());
    if (ifile) {
        // The file exists, and is open for input
        cout << "File exists! Owerite?" << endl;
        char t;
        cin >> t;
        if( t == 'y'){
            remove( fits_filename.c_str() );
            fitshandle fitswriter;
            fitswriter.create(fits_filename);
            write_Healpix_map_to_fits<double>(fitswriter, outputmap, PLANCK_FLOAT64 );
            fitswriter.close();
        }
    }else{
        cout << "Creating fits!..." << endl;
        fitshandle fitswriter;
        fitswriter.create(fits_filename);
        write_Healpix_map_to_fits<double>(fitswriter, outputmap, PLANCK_FLOAT64 );
        fitswriter.close();
    }

}

bool MapGenerator::isFinished(){
    return isFinished_;
}

double * MapGenerator::getMap(){
    return map_;
}

int MapGenerator::getNside(){
    return par_->map.Nside;
}