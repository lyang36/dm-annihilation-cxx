#include <fstream>
#include <iostream>
#include <string>
#include <cstring>

#include <rpc/types.h>
#include <rpc/xdr.h>

#include "../tipsydefs.h"
#include "../datareader.h"
int main(){
        string basedir = "/home/lyang/files/r200";
        string basename =  "vl2b.00400.r200";
        string nativefile = "/home/lyang/data/data_float_all.bin";
        DataReader readerxdr(basedir, basename);
        DataReader readerna(nativefile);
        readerxdr.open();
        readerna.open();

        while(readerxdr.hasNext()){
            DMParticle * dm1 = readerxdr.getBuf();
            DMParticle * dm2 = readerna.getBuf();
            printf("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",
                            dm1->mass - dm2->mass,
                            dm1->dens - dm2->dens,
                            dm1->hsmooth - dm2->hsmooth,
                            dm1->posx - dm2->posx,
                            dm1->posy - dm2->posy,
                            dm1->posz - dm2->posz,
                            dm1->velx - dm2->velx,
                            dm1->vely - dm2->vely,
                            dm1->velz - dm2->velz,
                            dm1->eps - dm2->eps,
                            dm1->phi - dm2->phi);
            readerxdr.loadBuffer();
            readerna.loadBuffer();
        }
        readerxdr.close();
        readerna.close();
}

