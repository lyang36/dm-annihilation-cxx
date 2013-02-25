#ifndef __PARAMETERS__LY__
#define __PARAMETERS__LY__
#include "datatypes.h"
class Parameters{
public:
    Parameters(std::string info_file);
	
    //the simulation parameters
    Info info;
	Units units;
	Natconst natconst;
	Halo halo;
	Codeunits codeunits;
	Params params;
    Dm dm;
    Map map;
	double rotmatrix[3][3];
    
    //particles in the memory
    int memParts;
    
    //use how many particles to test?
    int testNum;
    
    //output all the parameters
    void printParameters();
    
    //for filename
    string baseDir, baseName;
    string nativeDatafile;
    string maskFileName;
    
    //for output
    string outputFileName;
    
    //whether use native data?
    bool isNative;
    bool isMask;    //whether a mask file is set up
    
private:
    double align_vector[3];
    double obs_pos[3];

    void setupUnits();
    void setupParams();
    void setupHalo();
    void setupRotation();
    void setupOthers();
    
    void initialize();
    
    double redshift;
    

};

#endif
