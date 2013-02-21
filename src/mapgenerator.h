#ifndef __MAP_GENERATOR__LY
#define __MAP_GENERATOR__LY
class MapGenerator{
public:
    MapGenerator(Parameters * par, DataReader * reader);
    void start();
    void saveToFits(std::string fits_file);
};

#endif

