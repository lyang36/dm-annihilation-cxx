#ifndef __KERNEL__LY
#define NUM_THREADS_PER_BLOCK 512

class renderpart{
public:
    float theta;
    float phi;
    float x, y, z;
    float angular_radius;
    float flux;
};



struct healpix_par{
    int nside;
    int nl2;
    int nl3;
    int nl4;
    int nsidesq;
    int ncap;
    int npix;
    float theta_per_pix;
};

__global__ void calcfluxGPU(int nside,
                       float * map,
                       int numOfParts,
                       renderpart * parts);

cudaError_t initializeCUDA(int nside, int numofparts);
cudaError_t calculateMapByGPU(renderpart * parts, int num_of_parts);
//get the map from the device
cudaError_t getMap(float * map);
void cudaCleaingUp();
#endif