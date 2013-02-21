#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "../tipsydefs.h"
#include "../parameters.h"
using namespace std;
int main(){
        string infofile = "VL2_info.ini";
        Parameters params(infofile);
        params.printParameters();
        return 0;
}

