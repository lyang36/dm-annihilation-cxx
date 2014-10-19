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
    //string basedir = "/home/lyang/files/r200";
    //string basename =  "vl2b.00400.r200";
    string partfile = "";
    string densfile = "";
    string hsmfile = "";
    string dspfile = "";
    string output = "data.bin";
    while (m<argc)
    {
        string arg = argv[m];
        if (arg == "-partfile") {partfile = argv[m+1]; m++;}
        else if (arg == "-densfile") {densfile = argv[m+1]; m++;}
        else if (arg == "-hsmfile") {hsmfile = argv[m+1]; m++;}
        else if (arg == "-dspfile") {dspfile = argv[m+1]; m++;}
        else if (arg == "-output") {output = argv[m+1]; m++; }
        else {
            cout << "Usage:" << endl;
            cout << "-partfile: XDR datafile" << endl;
            cout << "-densfile: density file, output by smooth" << endl;
            cout << "-hsmfile: hsmooth radius file, output by smooth" << endl;
            cout << "-dspfile: velocity disp file, output by smooth" << endl;
            cout << "-output name: setup output file. Default: data.bin" << endl;
            cout << "-help: watching help" << endl;
            exit(0);
        }
        m ++;
    }
    
    DataReader reader(partfile, densfile, hsmfile, dspfile);
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
