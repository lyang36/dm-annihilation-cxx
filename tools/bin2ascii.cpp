/**
 * convert a native binary tipsy array to a ascii file
 * the format is:
 * [32 bit integer]  //for number of partilces
 * [value[i]]        //the value, in certain format
 */

#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>

#define CHAR 0
#define DOUBLE 1
#define FLOAT 2
#define INTEGER 3

using namespace std;

void printusage(){
    printf("Usage: bin2ascii binfile outfile <[-c]/[-d]/[-f]/[-i]>\n");
    printf("ascfile: the input binary file\n");
    printf("outfile: the output ascii file\n");
    printf("-c: char, -d: double, -f: float, -i: integer\n");
}

int main(int argn, char ** args){
    int format = CHAR;
    if(argn != 4){
        printusage();
        exit(1);
    }

    string infile = args[1];
    string outfile = args[2];
    if(strcmp(args[3], "-c") == 0){
        format = CHAR;
    }else if(strcmp(args[3], "-d") == 0){
        format = DOUBLE;
    }else if(strcmp(args[3], "-f") == 0){
        format = FLOAT;
    }else if(strcmp(args[3], "-i") == 0){
        format = INTEGER;
    }
    
    ifstream inputstr(infile.c_str(), ios::binary | ios::in);
    if(!inputstr.good()){
        printf("Input file corrupted!\n");
        exit(1);
    }

    int partnums;
    char valc;
    double vald;
    float valf;
    int vali;

    inputstr.read((char *) &partnums, sizeof(int));
    printf("Number of Particles: %d\n", partnums);

    ofstream outputstr(outfile.c_str());

    if(!outputstr.good()){
        printf("Cannot open/create output file!\n");
        exit(1);
    }

    outputstr << partnums;
    for(int i = 0; i < partnums; i++){
        int tempu;
        switch (format) {
            case CHAR:
                inputstr.read((char *) &valc, sizeof(char));
                tempu = valc;
                outputstr << "\n";
                outputstr << tempu;
                break;
                
            case DOUBLE:
                inputstr.read((char *) &vald, sizeof(double));
                outputstr << "\n";
                outputstr << vald;
                break;
                
            case FLOAT:
                inputstr.read((char *) &valf, sizeof(float));
                outputstr << "\n";
                outputstr << valf;
                break;
                
            case INTEGER:
                inputstr.read((char *) &vali, sizeof(int));
                outputstr << "\n";
                outputstr << vali;
                break;
                
            default:
                break;
        }
    }
    
    outputstr.close();
    inputstr.close();
    
    return 0;
}
