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
    bool verbose = false;
    string info_file =  "VL2_info.ini";
    while (m<argc)
    {
        string arg = argv[m];
        if (arg == "-info") { info_file = argv[m+1]; m+=1;}
        else if (arg == "-verbose") { verbose = true;}
        else{
            cout << "Usage:" << endl;
            cout << "-info name: setup info file path. Default: VL2_info.ini" << endl;
            cout << "-help: watching help" << endl;
            exit(0);
        }
        m++;
    }
    if(verbose)
        cout << "INFO FILE: " << info_file << endl;
    
    Parameters par(info_file);

    if(verbose)
        par.printParameters();

    DataReader * reader;

    if(par.isNative)
        reader = new DataReader(par.nativeDatafile);
    else
        reader = new DataReader(par.baseDir, par.baseName);
    reader->setTest(par.testNum);

    MapGenerator generator(&par, reader);

    printf("Start generating dark matter map:...\n");
    generator.start();
    generator.saveToFits(par.outputFileName);
    printf("Finished....\nMap saved as fits file in %s.\n", par.outputFileName.c_str());

    delete reader;
}
