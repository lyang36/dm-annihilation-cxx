#ifndef __KERNEL__LY
#define NUM_THREADS_PER_BLOCK 512

class renderpart{
public:
    //these information must be setup
    float theta;
    float phi;
    float x, y, z;
    float angular_radius;
    float flux;
    
    
    //the ring number
    int iring;
    int icol;
    
    //the difference in col
    int rmin;
    int rmax;
    int numPix;
    int npixNorthPole;
    int npixSouthPole;
    
    //giving theta, phi, x, y, z, angular_radius
    void setup(int nside);
    
};



class healpix_par{
public:
    int nside;
    int nl2;
    int nl3;
    int nl4;
    int nsidesq;
    int ncap;
    int npix;
    float theta_per_pix;
    void setup(int nside){
        this->nl2 = 2*nside;
        this->nl4 = 4*nside;
        this->nl3 = 3*nside;
        this->nside = nside;
        this->nsidesq = nside * nside;
        
        // ! number of pixels in the north polar cap
        this->ncap  = this->nl2*(nside-1);
        this->npix  = 12*nside*nside;
        this->theta_per_pix = 2 * M_PI / this->nl4;
    };
};

__global__ void calcfluxGPU(healpix_par pars,
                            float * map,
                            int numOfParts,
                            renderpart * parts);

cudaError_t initializeCUDA(int nside, int numofparts);
cudaError_t calculateMapByGPU(renderpart * parts, int num_of_parts);
//get the map from the device
cudaError_t getCUDAMap(float * map);
void cudaCleaingUp();
int ring_above (long nside_, float z);
int pix2ring(healpix_par &par, int ipix);
int pix2icol(healpix_par &par, int ring, int pix);
int cr2pix(healpix_par &par, int col, int ring);
void pix2vec(healpix_par &par, int r, int c, float &x, float &y, float &z, float &ct, float &phi);
void angle2pix(healpix_par &par,
               float z,
               float phi,
               int & iring,
               int & icol,
               int & ipix);
#endif
