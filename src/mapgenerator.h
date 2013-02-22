#ifndef __MAP_GENERATOR__LY
#define __MAP_GENERATOR__LY
#include "parameters.h"
#include "datareader.h"

//define the map's type
typedef double MAPTYPE;

//float can lead to incorrect result
//typedef float MAPTYPE;

class MapGenerator{
public:
    MapGenerator(Parameters * par, DataReader * reader);
    ~MapGenerator();
    void start();
    void saveToFits(std::string fits_file);
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

