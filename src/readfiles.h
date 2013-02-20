#ifndef __READFILES__
#define __READFILES__
#include <string>
#include "tipsydefs.h"
#include "datatypes.h"

/**
 * This class read particles in to the array, all the arrays should alread be allocated memory
 */

using namespace std;

class XDRFileReader{
	
public:
    //function reader the tipsheader
    static void read_tipsyheader(string filename, tipsy_header * tipsyheader);
    
    //Create the vector and read darkmatter particles from XDR files
    static void read_particles(string filename, Pdm * particles);
    
    //Reader scalers from the the file, return back the number of particles in np
    static void read_scalar(string filename, float * data, int &np);
};

#endif
