#include<stdlib.h>
#include<stdio.h>

#include <rpc/types.h>
#include <rpc/xdr.h>

#include "tipsydefs.h"

int read_tipsy_header(char *filename, tipsy_header *header) {
  int status;
  FILE *fp;
  XDR xdr;

  fp=fopen(filename,"r");
  if(!fp) {
    fprintf(stderr,"Problem opening file: %s\n",filename);
    return FAILURE;
  }

  /* Use the following if you ever come across a TIPSY file in native
     binary format (as opposed to XDR). */

/*   if(fread(header,sizeof(tipsy_header),1,fp) == 0) { */
/*     fprintf(stderr,"Problem reading tipsy header.\n"); */
/*     return FAILURE; */
/*   } */

  xdrstdio_create(&xdr,fp,XDR_DECODE);
  status = xdr_header(&xdr,header);
  if(status != 1) {
    fprintf(stderr,"Problem reading tipsy header.\n");
    return FAILURE;
  }
  xdr_destroy(&xdr);
  fclose(fp);

  return SUCCESS;
}


int read_tipsy_dm_particles(char *filename, tipsy_header header, Pdm *all_dm_particles) {
  FILE *fp;
  int i,j,status;
  Pdm dp;
  XDR xdr;

  fp=fopen(filename,"r");
  if(!fp) {
    fprintf(stderr,"Problem opening file: %s\n",filename);
    return FAILURE;
  }

  /* Use the following if you ever come across a TIPSY file in native
     binary format (as opposed to XDR). */

/*   status=fseek(fp,sizeof(header)+header.nsph*sizeof(Psph),SEEK_SET); */
/*   fread(all_dm_particles,sizeof(Pdm),header.ndark,fp); */


  xdrstdio_create(&xdr,fp,XDR_DECODE);
  status = xdr_header(&xdr,&header);  

  for(i=0;i<header.ndark;i++) {
    status = xdr_dark(&xdr,&dp);
    if (status != 1) {
      fprintf(stderr,"Error reading dark particle from input file.\n");
      exit(1);
    }
    for(j=0;j<3;j++) {
      all_dm_particles[i].pos[j] = dp.pos[j];
      all_dm_particles[i].vel[j] = dp.vel[j];
    }
    all_dm_particles[i].mass = dp.mass;
    all_dm_particles[i].eps = dp.eps;
    all_dm_particles[i].phi = dp.phi;
  }
  xdr_destroy(&xdr);

  fclose(fp);

  return SUCCESS;
}



int write_tipsy_file(char *filename, tipsy_header header, Pdm *all_dm_particles, int Nhires, float hires_mass) {
  FILE *fp;
  int i,status;
  int Ntot;
  XDR xdr;

  Ntot = header.ndark;

  header.nbodies = Nhires;
  header.ndark = Nhires;
  header.nsph = 0;
  header.nstar = 0;

  fp=fopen(filename,"w");
  if(!fp) {
    fprintf(stderr,"Problem creating file: %s\n",filename);
    return FAILURE;
  }

  xdrstdio_create(&xdr,fp,XDR_ENCODE);
  status = xdr_header(&xdr,&header);  

  for(i=0;i<Ntot;i++) {

    if(all_dm_particles[i].mass > hires_mass*1.001) continue;

    status = xdr_dark(&xdr,&(all_dm_particles[i]));
    if (status != 1) {
      fprintf(stderr,"Error writing dark particle to output file.\n");
      exit(1);
    }
  }
  xdr_destroy(&xdr);
  
  fclose(fp);

  return SUCCESS;
}



int xdr_header(XDR *xdrs, tipsy_header *header) {
  int pad=0;
  
  if (xdr_double(xdrs,&header->time) != TRUE) return 0;
  if (xdr_int(xdrs,&header->nbodies) != TRUE) return 0;
  if (xdr_int(xdrs,&header->ndim) != TRUE) return 0;
  if (xdr_int(xdrs,&header->nsph) != TRUE) return 0;
  if (xdr_int(xdrs,&header->ndark) != TRUE) return 0;
  if (xdr_int(xdrs,&header->nstar) != TRUE) return 0;
  if (xdr_int(xdrs,&pad) != TRUE) return 0;
  return 1;
}

int xdr_dark(XDR *xdrs,struct dark_particle *p) {

  if (xdr_float(xdrs,&p->mass) != TRUE) return 0;
#ifdef POS_FLOAT
  if (xdr_float(xdrs,&p->pos[0]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->pos[1]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->pos[2]) != TRUE) return 0;
#else
  if (xdr_double(xdrs,&p->pos[0]) != TRUE) return 0;
  if (xdr_double(xdrs,&p->pos[1]) != TRUE) return 0;
  if (xdr_double(xdrs,&p->pos[2]) != TRUE) return 0;
#endif
  if (xdr_float(xdrs,&p->vel[0]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->vel[1]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->vel[2]) != TRUE) return 0;
  if (xdr_float(xdrs,&p->eps) != TRUE) return 0;
  if (xdr_float(xdrs,&p->phi) != TRUE) return 0;


  return 1;
}  
