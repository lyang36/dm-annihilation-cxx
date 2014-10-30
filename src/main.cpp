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

    if(par.isNative){
        reader = new DataReader(par.nativeDatafile, par.isMask, par.maskFileName);
    }else{
        //reader = new DataReader(par.baseDir, par.baseName, par.isMask, par.maskFileName);
        reader = new DataReader(
                par.partFile, 
                par.densFile,
                par.hsmFile,
                par.dspFile,
                par.isMask,
                par.maskFileName);
    }
    reader->setTest(par.testNum);
    reader->setBuf(par.memParts);
    reader->open();

    MapGenerator generator(&par, reader);

    printf("Start generating dark matter map:...\n");
    generator.start();
    reader->close();

    generator.saveToFits(par.outputFileName, verbose);
    printf("Finished....\nMap saved as fits file in %s.\n", par.outputFileName.c_str());
   
    string command = "map2tga " + par.outputFileName + " " + par.outputFileName + ".tga -log -bar -lon 180";
    if(0 == system(command.c_str()))
        //Do nothing;

    delete reader;
}
