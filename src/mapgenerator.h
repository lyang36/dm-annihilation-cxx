#ifndef __MAP_GENERATOR__LY
#define __MAP_GENERATOR__LY
#include "parameters.h"
#include "datareader.h"



class MapGenerator{
public:
    MapGenerator(Parameters * par, DataReader * reader);
    ~MapGenerator();
    void start();
    void saveToFits(std::string fits_file, bool isVerbose = false);
    bool isFinished();
    MAPTYPE * getMap();
    int getNside();
    
    
    
private:
    Parameters * par_;
    DataReader * reader_;
    bool isFinished_;
    MAPTYPE * map_;
    
};

#endif

