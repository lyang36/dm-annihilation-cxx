#ifndef __LY__VL2READ__
#define __LY__VL2READ__



class VL2Reader{
public:
    // Read the native file
    VL2Reader(std::string filename, double xc, double yc, double zc, double boxsize, std::string maskfile="");
    
    // If isSingleFile = true, then read the only the tipsy file, and density are all set to be zero
    VL2Reader(std::string basename, std::string filename, double xc, double yc, double zc, double boxsize, std::string maskfile="", bool isSingleFile = false);
    
    std::vector<float> & getDensity();
    std::vector<float> & getVelocityX();
    std::vector<float> & getVelocityY();
    std::vector<float> & getVelocityZ();
    
    static const int INITMEMSIZE;
    
private:
    std::vector<float> density_;
    std::vector<float> velX_;
    std::vector<float> velY_;
    std::vector<float> velZ_;
    
    void loadParts(std::string basename, std::string filename, double xc, double yc, double zc, double boxsize, std::string maskfile, bool isSingleFile, bool isTipsy);
};
#endif