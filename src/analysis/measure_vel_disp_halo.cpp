/**
 * Measure the velocity dispersion of a list of halos
 * Author: Lin F. Yang
 * */
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <limits>
#include <ios>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>

#include <rpc/types.h>
#include <rpc/xdr.h>

#include "../tipsydefs.h"

using namespace std;


/**Given a list of IDs, return the halo catalog.
 * fn is the AHF halo catalog file
 * the ids should be sorted
 * */
vector<set<int> > loadHaloCataLog(vector<int> ids, string fn);

/**
 * Load halo particles from the tipsyfile
 */
vector<vector<Pdm> > loadHaloParticles(vector<set<int> > haloPartIds, string fn);

/**
 * Load the selected halo ids from the file
 * the file is a txt file contains a list 
 * of IDs. The first integer includes the
 * number of ids in it.
 * The returned id will be sorted
 * */
vector<int> loadSelectedHaloIds(string fn);


int main(int argc, const char **argv){
   string datafilename = 
       "/grouphome/gwln2/data2/dss001/vialactea/snapshots/vl2b.00400";
   string halocatalog = 
       "/grouphome/gwln2/home/lyang/data/VL2Halos/vl_400_rhovesc.z0.000.AHF_particles";
   string haloselected = 
       "/home/lyang/iPythonNoteBooks/darkmatter_clumpiness/subhalos_6_10.txt";

   /*
    * The format in ths file is: 
    * #num of halos
    * #ID Velocitydisp
    * */
   string outputfile_name = 
       "/home/lyang/iPythonNoteBooks/darkmatter_clumpiness/haloveldisp_6_10.txt";
    


   vector<int> selectedIds = loadSelectedHaloIds(haloselected);

   vector<set<int> > haloPartIds = loadHaloCataLog(selectedIds, halocatalog);

   vector<vector<Pdm> > haloParts = loadHaloParticles(haloPartIds, datafilename);


   vector<double> average_velmag(haloPartIds.size(), 0.0);
   vector<double> average_velsq(haloPartIds.size(), 0.0);

   vector<double> average_velx(haloPartIds.size(), 0.0);
   vector<double> average_vely(haloPartIds.size(), 0.0);
   vector<double> average_velz(haloPartIds.size(), 0.0);
   
   vector<double> average_velsqx(haloPartIds.size(), 0.0);
   vector<double> average_velsqy(haloPartIds.size(), 0.0);
   vector<double> average_velsqz(haloPartIds.size(), 0.0);


   cout << "Measuring velocity dispersion...\n";

   for(unsigned int i = 0; i < haloPartIds.size(); i++){
       for(unsigned int j = 0; j < haloParts[i].size(); j++){
           double velx, vely, velz;
           
           velx = haloParts[i][j].vel[0];
           vely = haloParts[i][j].vel[1];
           velz = haloParts[i][j].vel[2];

           double velsq = velx * velx + vely * vely + velz * velz;
           
           average_velmag[i] += sqrt(velsq);
           average_velsq[i] += velsq;

           average_velx[i] += velx;
           average_vely[i] += vely;
           average_velz[i] += velz;

           average_velsqx[i] += velx * velx;
           average_velsqy[i] += vely * vely;
           average_velsqz[i] += velz * velz;
       }
   }

   ofstream oputfile(outputfile_name.c_str());
   if(oputfile.good()){
       oputfile << "#num of halos\n"; 
       oputfile << "#ID <vx> <vy> <vz> <vx^2> <vy^2> <vz^2> <v> <v^2>\n";
       for(unsigned int i = 0; i < haloPartIds.size(); i++){
            oputfile << average_velx[i] << " ";
            oputfile << average_vely[i] << " ";
            oputfile << average_velz[i] << " ";

            oputfile << average_velsqx[i] << " ";
            oputfile << average_velsqy[i] << " ";
            oputfile << average_velsqz[i] << " ";

            oputfile << average_velmag[i] << " ";
            oputfile << average_velsq[i] << endl;

       }
   }
   oputfile.close();

}


int getLineInt(ifstream & fs){
    int tmp;
    string tmpline;
    getline(fs, tmpline);
    stringstream tmplinestm(tmpline);
    tmplinestm >> tmp;
    return tmp;
}

// read from halocatalog
vector<set<int> > loadHaloCataLog(vector<int> ids, string fn){
    cout << "Load the selected halo catalogs.\n";

    ifstream halocatalog(fn.c_str());
    vector<set<int> > haloPartIds;
    if(halocatalog.good()){
        int numHalos;
        numHalos = getLineInt(halocatalog);
        unsigned int halopt = 0;
        for(int i = 0; i < numHalos; i++){
            int numparts;
            numparts = getLineInt(halocatalog);
            
            {
                printf("\r %2.2f%% %d %d            ", 
                        (double)(i) / (double)(numHalos) * 100.0, 
                        i, numparts);
                cout.flush();
            }
            
            
            set<int> halopids;
            for(int j = 0; j < numparts; j++){
                string line;
                //getline(halocatalog, line);
                //stringstream linestream(line);
                if(i != ids[halopt]){
                    halocatalog.ignore(numeric_limits<streamsize>::max(), 
                           '\n');
                }else{
                    getline(halocatalog, line);
                    stringstream linestream(line);
                    int id, ptype;
                    linestream >> id;
                    linestream >> ptype;
                    halopids.insert(id);
                }
            }

            if(i == ids[halopt]){
                haloPartIds.push_back(halopids);
                halopt ++;
            }

            if(halopt >= ids.size()){
                break;
            }
        }
        printf("\n");
        halocatalog.close();
    }else{
        printf("Catalog File Error!\n");
        exit(1);
    }


    return haloPartIds;
}


vector<vector<Pdm> > loadHaloParticles(vector<set<int> > haloPartIds, string fn){
    cout << "Load the halo particles.\n";

    vector<vector<Pdm> > particles(haloPartIds.size());
    
    //int id = 0;
    set<int> totalSet;

    for(unsigned int i = 0; i < haloPartIds.size(); i++){
        totalSet.insert(haloPartIds[i].begin(), haloPartIds[i].end());
    }
    
    tipsy_header header;
    XDR xdr;
    Pdm dp;
    FILE *fp;

    read_tipsy_header(fn.c_str(), &header);
    printf("Number of particles = %d\n",header.nbodies);
    printf("Number of DM particles = %d\n",header.ndark);

    fp = fopen(fn.c_str(),"r");
    xdrstdio_create(&xdr,fp,XDR_DECODE);
    int status = xdr_header(&xdr,&header);
    for(int i=0; i<header.ndark; i++) {
        if (status != 1) {
            fprintf(stderr,"Error reading dark particle from input file.\n");
            exit(1);

        }
        status = xdr_dark(&xdr,&dp);
        if(i % 10000 == 0){
            printf("\r %2.2f%%", 100.0 * (double)(i) / (double)(header.ndark));
            cout.flush();
        }
        if(totalSet.end() != totalSet.find(i)){
            //insert into the list
            for(unsigned int j = 0; j < haloPartIds.size(); j++){
                if(haloPartIds[j].end() != haloPartIds[j].find(i)){
                    particles[j].push_back(dp);
                }
            }
        }
    }
    printf("\n"); 
    fclose(fp);

    return particles;

}

vector<int> loadSelectedHaloIds(string fn){
    cout << "Load the selected IDs.\n";
    vector<int> selecteHaloIds;
    ifstream ist(fn.c_str());
    if(!ist.good()){
        printf("The selected halo ID file corrupted?!\n");
        exit(1);
    }


    int numhalos;
    int id;

    ist >> numhalos;
    selecteHaloIds.reserve(numhalos);

    for(int i = 0; i < numhalos; i ++){
        ist >> id;
        selecteHaloIds.push_back(id);
    }
    sort(selecteHaloIds.begin(), selecteHaloIds.end());

    return selecteHaloIds;
}
