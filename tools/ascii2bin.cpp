/**
* convert a ascii tipsy array to a native binary file
* the format is:
* [32 bit integer]  //for number of partilces
* [value[i]]        //the value, in certain format
*/


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
    printf("Usage: ascii2bin ascfile outfile <[-c]/[-d]/[-f]/[-i]>\n");
    printf("ascfile: the input ascii file\n");
    printf("outfile: the output file\n");
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
    
    ifstream inputstr(infile.c_str());
    if(!inputstr.good()){
        printf("Input file corrupted!\n");
        exit(1);
    }
    
    int partnums;
    char valc;
    double vald;
    float valf;
    int vali;
    
    inputstr >> partnums;
    printf("Number of Particles: %d\n", partnums);
    
    ofstream outputstr(outfile.c_str(), ios::binary | ios::out);
    if(!outputstr.good()){
        printf("Cannot open/create output file!\n");
        exit(1);
    }
    
    outputstr.write((char *) &partnums, sizeof(int));
    
    for(int i = 0; i < partnums; i++){
        switch (format) {
            case CHAR:
                int intec;
                inputstr >> intec;
                valc = intec;
                outputstr.write((char *) &valc, sizeof(char));
                //printf("%d\n", valc);
                break;
            case DOUBLE:
                inputstr >> vald;
                outputstr.write((char *) &vald, sizeof(double));
                //printf("%f\n", vald);
                break;
            case FLOAT:
                inputstr >> valf;
                outputstr.write((char *) &valf, sizeof(float));
                //printf("%f\n", valf);
                break;
            case INTEGER:
                inputstr >> vali;
                outputstr.write((char *) &vali, sizeof(int));
                //printf("%d\n", vali);
                break;
            default:
                break;
        }
    }
    
    outputstr.close();
    inputstr.close();
    
    return 0;
    
}
