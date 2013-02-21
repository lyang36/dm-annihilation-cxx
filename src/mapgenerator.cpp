#include <string>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "mapgenerator.h"

using namespace std;

MapGenerator::MapGenerator(Parameters * par, DataReader * reader){
    par_ = par;
    reader_ = reader;
}

void MapGenerator::start(){
    
}

void MapGenerator::saveToFits(string filename){
    
}