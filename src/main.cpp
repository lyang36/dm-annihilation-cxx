#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <rpc/types.h>
#include <rpc/xdr.h>

#include "parameters.h"
#include "datareader.h"
#include "mapgenerator.h"
#include "datatypes.h"

using namespace std;

int main(int argc, const char **argv){
    int m=1;
    string basedir = "/home/lyang/files/r200";
    string basename =  "vl2b.00400.r200";
    string info_file =  "VL2_info.txt";
    string fits_file = "skymap.fits";
    string par_file = "configure.ini";
    while (m<argc)
    {
        string arg = argv[m];
        if (arg == "-base") { basedir = argv[m+1]; m+=2;}
        if (arg == "-bname") { basename = argv[m+1]; m+=2;}
        if (arg == "-info") { info_file = argv[m+1]; m+=2;}
        if (arg == "-fits") { fits_file = argv[m+1]; m+=2;}
        if (arg == "-conf") { par_file = argv[m+1]; m+=2;}
        if (arg == "-help") {
            cout << "Usage:" << endl;
            cout << "-base name: setup basedir of the data. Default: /home/lyang/files/r200" << endl;
            cout << "-bname name: setup basename of the data. Default: vl2b.00400.r200" << endl;
            cout << "-info name: setup info file path. Default: VL2_info.txt" << endl;
            cout << "-fits name: setup output fitsname. Default: skymap.fits" << endl;
            cout << "-conf name: setup configure filename. Default: configure.ini" << endl;
            cout << "-help: watching help" << endl;
            exit(0);
        }
    }
    cout << "BASE DIR: " << basedir << endl;
    cout << "BASE NAME: " << basename << endl;
    cout << "INFO FILE: " << info_file << endl;
    cout << "OUTPUT FITS: " << fits_file << endl;
    cout << "CONFIGUR FILE: " << par_file << endl;
    
    Parameters par(info_file, par_file);
    DataReader reader(basedir, basename);
    MapGenerator generator(&par, &reader);
    printf("Start generating dark matter map:...\n");
    generator.start();
    generator.saveToFits(fits_file);
    printf("Finished. Map saved as fits file in %s.\n", fits_file.c_str());
}
