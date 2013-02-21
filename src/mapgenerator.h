#ifndef __MAP_GENERATOR__LY
#define __MAP_GENERATOR__LY
#include "parameters.h"
#include "datareader.h"
class MapGenerator{
public:
    MapGenerator(Parameters * par, DataReader * reader);
    void start();
    void saveToFits(std::string fits_file);
    
private:
    Parameters * par_;
    DataReader * reader_;
};

#endif

