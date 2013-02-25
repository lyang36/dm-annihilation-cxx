/* DataReader
 * Read data from the data file
 * set up a data file path
 * open it
 * use it untill done
 * Lin Yang 06/05/2012
 */

#ifndef __LY__DATAREADER__
#define __LY__DATAREADER__
#include "datatypes.h"
#include <string>
#include <fstream>

using namespace std;
class DataReader{
public:
    int getPartNum(){
        if(testnum_ != -1){
            return testnum_;
        }else{
            return partNumbers_;
        }  
    };
    
    //to use this reader to read the native file
    DataReader(string path, bool isMask = false, string maskfile = ""); //set path
    
    //to use this reader to read the XDR file
    //it should hase a basedir, and a basename
    DataReader(string basedir, string basename, bool isMask = false, string maskfile = "");
    
    void setTest(int tn);                           //set test: not read all the particles
    
    void loadBuffer();                              //load buffer
    DMParticle * getBuf();                         //get buffer
    int getMemparts(){
        return memParts_;
    };
    void move2bufEnd();                             //move the memCursor to then end of the memory
    unsigned int getBufSize();                      //get mparts
    
    bool open();                                    //open the file and load the buffer
    
    bool isOpen();                                  //check is open?
    void setBuf(unsigned int mem);                  //set up the buffer
    bool hasNext();                                 //whether the data is over?
    bool readParticle(DMParticle * part);          //read a particle from the file
    
    void close();
    
private:
    bool isNative_;
    bool isMask_;
    
    //for XDR file
    string basedir_, basename_;
    
    string particle_file_name;
    string density_file_name;
    string hsmooth_file_name;
    string sigmav_file_name;
    string mask_file_name;   

    tipsy_header tips_;
    FILE * fp_;
	XDR xdrs_;
    ifstream densStream_;
    ifstream hsmoothStream_;
    ifstream sigmavStream_;
    ifstream maskStream_;

    //for native file
    string filePath_;//file path
    
    unsigned int partNumbers_;//particle numbers in the file
    unsigned int readCursor_; //the read Cursor of the file
    bool isopen_d;  //is open?
    unsigned int memBuffer_;  //how many particles should be load in the buffer
    //default 10^6
    
    unsigned int memCursor_;  //reading from the memory cursor
    unsigned int memParts_;   //particles in membuffer
    DMParticle * buffer_; //particle buffer
    ifstream dataInputFile_;//input stream
    int testnum_;   //for testing
};
#endif
