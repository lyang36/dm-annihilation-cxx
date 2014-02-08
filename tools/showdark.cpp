#include <string>
#include <cstring>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <algorithm>    // std::sort
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "../src/tipsydefs.h"

using namespace std;

int main(int argn, char ** argv){
    if(argn != 4){
        printf("Usage: showdark file startind endind\n");
        exit(1);
    }
    
    string partfile = argv[1];

    int startind = 0;
    int endind = 0;
    
    stringstream ss0, ss1;
    ss0 << argv[2];
    ss0 >> startind;
    ss1 << argv[3];
    ss1 >> endind;
    
    Pdm dp;
    tipsy_header header;
    XDR xdr;
    
    FILE *fp;
    read_tipsy_header(partfile.c_str(), &header);
    printf("Number of particles = %d\n",header.nbodies);
    printf("Number of DM particles = %d\n",header.ndark);
    
    fp = fopen(partfile.c_str(),"r");
    xdrstdio_create(&xdr,fp,XDR_DECODE);
    int status = xdr_header(&xdr,&header);
    
    for(int i=0; (i<header.ndark) && (i < endind); i++) {
        status = xdr_dark(&xdr,&dp);
        
        if(i < startind){
                continue;
        }
        if (status != 1) {
            fprintf(stderr,"Error reading dark particle from input file.\n");
            exit(1);
        }
        printf("%d > \n", i);
        printf("mass = %e\n", dp.mass);
        printf("pos = [%f %f %f]\n", dp.pos[0], dp.pos[1], dp.pos[2]);
        printf("vel = [%f %f %f]\n", dp.vel[0], dp.vel[1], dp.vel[2]);
        printf("eps = %f\n", dp.eps);
        printf("phi = %f\n", dp.phi);
        printf("\n");
    }
    printf("\n");
    
    fclose(fp);
    
    return(0);
}
