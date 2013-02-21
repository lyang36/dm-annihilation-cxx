#ifndef __TIPSDEF__
#define __TIPSDEF__

#define POS_FLOAT
#ifdef POS_FLOAT
  #define Real float
#else
  #define Real double
#endif

#define FAILURE 1
#define SUCCESS 0
//#include <rpc/types.h>
//#include <rpc/rpc.h>

struct gas_particle {
  float mass;
  Real pos[3];
  float vel[3];
  float rho;
  float temp;
  float hsmooth;
  float metals ;
  float phi ;
} ;

struct dark_particle {
    float mass;
    Real pos[3];
    float vel[3];
    float eps;
    float phi ;
} ;

struct star_particle {
    float mass;
    Real pos[3];
    float vel[3];
    float metals ;
    float tform ;
    float eps;
    float phi ;
} ;

struct dump {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
} ;

typedef struct gas_particle Pgas;
typedef struct dark_particle Pdm;
typedef struct star_particle Pstar;
typedef struct dump tipsy_header;




int xdr_header(XDR *xdrs, tipsy_header *header);
int xdr_dark(XDR *xdrs,struct dark_particle *p);

int read_tipsy_header(char *filename, tipsy_header *header);
int read_tipsy_dm_particles(char *filename, tipsy_header header, Pdm *all_dm_particles);
int write_tipsy_file(char *filename, tipsy_header header, Pdm *all_dm_particles, int Nhires, float hires_mass);

#endif

