#include <fstream>
#include <string>
#include <sstream>

using namespace std;

int main(int argn, char ** args){
    if(argn != 4){
        printf("Usage: id2mark totaln idfile markfile\n");
        printf("totaln: the total number of particles in the file\n");
        printf("idfile: the particle ids of the halo\n");
        printf("markfile: the output ascii tipsy array\n");
    }
    
    int totaln;
    string idfile;
    string markfile;
    stringstream ss;
    int * partlist;
    int partnum = 0;
    int partind = 0;
    
    ss << args[1];
    ss >> totaln;
    
    idfile = args[2];
    markfile = args[3];
    
    {
        string line;
        ifstream myfile0(idfile.c_str());
        while (std::getline(myfile0, line)){
            ++partnum;
        }
        myfile0.close();
        
        //test
        //printf("Particle numbers: %d\n", partnum);
        printf("Will select %d partilces from simulation results.\n", partnum);
        
        ifstream myfile1(idfile.c_str());
        partlist = new int[partnum];
        for(int i = 0; i < partnum; i++){
            myfile1 >> partlist[i];
        }
        myfile1.close();
    }
    
    //sort the particle list
    printf("Sorting IDs ...\n");
    std::sort (partlist, partlist + partnum);
    
    ofstream markstream(markfile.c_str());
    markstream << totaln;
    for(int i = 0; i < totaln; i++){
        if((partind < partnum) && (i == partlist[partind])){
            markstream << "\n1";
            partind ++;
        }else{
            markstream << "\n0";
        }
    }
    markstream.close();
    delete partlist;
    return 0;
}