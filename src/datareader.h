#ifndef __DATAREADER__LY__
#define __DATAREADER__LY__
using namespace std;
class DataReader{
public:
    DataReader(std::string basedir, string basename);
    ~DataReader();
    
    //read all the particles into memory
    bool readParticle();
    
    //get the all the particles
    DMParticle * getParticles();
    
    //get the number of the parts
    int getNumParts();
private:
    DMParticle * particles;
    string particle_file_name;
    string density_file_name;
    string hsmooth_file_name;
    int numParts_;
    
};

#endif
