/*************************************************************/
//					Class: ReadFiles
//			Read particle data from files
//
/**************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <rpc/xdr.h>
#include <rpc/types.h>
#include "tipsydefs.h"
#include "readfiles.h"


using namespace std;


void XDRFileReader::read_tipsyheader(string filename, tipsy_header * tipsyheader){
	char fi[255];
	strcpy(fi, filename.c_str());
	read_tipsy_header(fi, tipsyheader);
}


void XDRFileReader::read_particles(string filename, Pdm * particles){
	long i,j;
	FILE * fp;
	XDR xdrs;
	char fi[255];
    
    tipsy_header tips;
    
	//keywords
	//positions keyword overrides the {xyz}pos keywords
    
	strcpy(fi, filename.c_str());
	fp=fopen(fi,"r");
	if(!fp) {
		fprintf(stderr,"Problem opening file: %s\n",filename.c_str());
		exit(1);
		return ;
	}
	
	xdrstdio_create(&xdrs,fp,XDR_DECODE);
    
	//read header
	int status = xdr_header(&xdrs, &tips);
    int Nparticles = tips.ndark;
    
    //create the vector
    particles = new Pdm[Nparticles];
    
	//buid all particles, this is a slow way
	for(i=0; i < Nparticles; i++) {
		status = xdr_dark(&xdrs,&(particles[i]));
		if (status != 1) {
			fprintf(stderr,"Error reading dark particle from input file.\n");
			exit(1);
		}
	}
    
	xdr_destroy(&xdrs);
	fclose(fp);
}


void XDRFileReader::read_scalar(string filename, float * scalar, int &np){
    ifstream myFile ( filename.c_str(), ios::in | ios::binary);
	myFile.read( reinterpret_cast<char*>( &np ), sizeof np );
    
    scalar = new float[np];
    
    if (!myFile.read (reinterpret_cast<char*>( scalar ), (sizeof(float)) * np )) {
		cout << "READ ERROR!"  << endl;
		exit(1);
    }
	myFile.close();

}

