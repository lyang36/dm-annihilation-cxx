#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "../tipsydefs.h"
#include "../datareader.h"

int main(int argc, const char **argv){
    int m=1;
    string basedir = "/home/lyang/files/r200";
    string basename =  "vl2b.00400.r200";
    string output = "data.bin";
    while (m<argc)
    {
        string arg = argv[m];
        if (arg == "-base") { basedir = argv[m+1]; m++;}
        else if (arg == "-bname") { basename = argv[m+1]; m++;}
        else if (arg == "-output") {output = argv[m+1]; m++; }
        else {
            cout << "Usage:" << endl;
            cout << "-base name: setup basedir of the data. Default: /home/lyang/files/r200" << endl;
            cout << "-bname name: setup basename of the data. Default: vl2b.00400.r200" << endl;
            cout << "-output name: setup output file. Default: data.bin" << endl;
            cout << "-help: watching help" << endl;
            exit(0);
        }
        m ++;
    }
    
    DataReader reader(basedir, basename);
    reader.open();
    
    cout << "writing data to files ... " << endl;
    int nums = reader.getPartNum();
    cout << "There are " << nums << " particles" << endl;
    
    ofstream myFile (output.c_str(), ios::out | ios::binary);;
    myFile.write ((char*)&nums, sizeof (nums));
    
    if(myFile.good()){
        cout << "starting ... "<< endl;
        int nn = 0;
        while(reader.hasNext()){
            nn += reader.getMemparts();
            printf("Particles: %d\n", nn);
            DMParticle *dp = reader.getBuf();
            //printf("%e %e %e %e %e %e %e %e %e %e %e\n", dp[0].mass, dp[0].dens,
            //                dp[0].hsmooth, dp[0].posx, dp[0].posy, dp[0].posz,
            //                dp[0].velx, dp[0].vely, dp[0].velz, dp[0].eps, dp[0].phi);
            cout.flush();
            myFile.write((char *) reader.getBuf(), sizeof(DMParticle) * reader.getMemparts()) ;
            reader.loadBuffer();
        }
        printf("\n");
        cout << "ending... Output file: " << output << endl;
        myFile.close();
    }
    reader.close();
    return 0;
}
