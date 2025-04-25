/* driver program for van Hove gdiff program */
#include <stdio.h>
#include <stdlib.h>
void gdiff_();
main(argc,argv)
int argc;
char* argv[];
{
  char atnam1[9],atnam2[9],fname[41],title[80];
  int dl_fmt,nxtim,ncrs,nxstr,nxrad,nxatms,natms,numfil;
  int nconf,isampl,iocrs,status;
  double tstep,rcut,cell[9],avcell[9],bcell[10];

  char *name,*newnam;
  int *imd,*msm;
  double *xyz,*vel,*frc,*rdf,*weight,*chge;
  double *gcross,*gzero,*xy0,*xy1,*xy2;
  
  status=0;
  
  /* input parameters from command line */

  dl_fmt=atoi(argv[1]);
  strcpy(atnam1,argv[2]);
  strcpy(atnam2,argv[3]);
  strcpy(fname,argv[4]);
  natms=atoi(argv[5]);
  nconf=atoi(argv[6]);
  ncrs=atoi(argv[7]);
  isampl=atoi(argv[8]);
  iocrs=atoi(argv[9]);
  numfil=atoi(argv[10]);
  nxstr=atoi(argv[11]);
  rcut=atof(argv[12]);
  nxrad=ncrs;
  nxtim=ncrs;
  nxatms=natms;
  name=(char *)malloc(natms*8);
  imd=(int *)malloc(nxtim*sizeof(int));
  msm=(int *)malloc(nxtim*sizeof(int));
  newnam=(char *)malloc(nxstr*8);
  gcross=(double *)malloc(nxtim*nxrad*sizeof(double));
  gzero=(double *)malloc(3*nxtim*nxstr*sizeof(double));
  xy0=(double *)malloc(3*nxstr*sizeof(double));
  xy1=(double *)malloc(3*nxstr*sizeof(double));
  xy2=(double *)malloc(3*nxstr*sizeof(double));
  xyz=(double *)malloc(3*natms*sizeof(double));
  vel=(double *)malloc(3*natms*sizeof(double));
  frc=(double *)malloc(3*natms*sizeof(double));
  weight=(double *)malloc(natms*sizeof(double));
  chge=(double *)malloc(natms*sizeof(double));

  /* now call fortran gdiff routine */  

  gdiff_(&dl_fmt,atnam1,atnam2,fname,title,name,newnam,&nxtim,&ncrs,
	&nxstr,&nxrad,&nxatms,&natms,&nconf,&isampl,&iocrs,&numfil,
	&tstep,&rcut,imd,msm,avcell,gcross,gzero,xyz,vel,frc,weight,
	 chge,cell,bcell,xy0,xy1,xy2,&status);

  exit(0);

}

