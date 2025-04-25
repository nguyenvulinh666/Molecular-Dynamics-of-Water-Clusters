/* driver program for density correlation program */
#include <stdio.h>
#include <stdlib.h>
void dencor_();
main(argc,argv)
int argc;
char* argv[];
{
  char fname[41],title[80];
  int lden,lcor,ltim,lchg,mcore,msp,natm,ncon,kmax,ntime,ngap;
  int status,numfil,klim,m1,m2,m3;
  double *xxx,*yyy,*zzz,*chge,*space;
  
  /* input parameters from command line */

  strcpy(fname,argv[1]);
  lden=atoi(argv[2]);
  lcor=atoi(argv[3]);
  ltim=atoi(argv[4]);
  lchg=atoi(argv[5]);
  natm=atoi(argv[6]);
  ncon=atoi(argv[7]);
  kmax=atoi(argv[8]);
  ntime=atoi(argv[9]);
  ngap=atoi(argv[10]);
  numfil=atoi(argv[11]);
  msp=natm;
  klim=((2*kmax+1)*(2*kmax+1)*(2*kmax+1)-1)/2;
  m1=6*msp*(kmax+2)+6*klim;
  m2=4*klim*ntime+2*klim+2*ntime;
  m3=ntime*(2*klim+19);
  mcore=m1;
  if(m2>mcore)mcore=m2;
  if(m3>mcore)mcore=m3;
  xxx=(double *)malloc(msp*sizeof(double));
  yyy=(double *)malloc(msp*sizeof(double));
  zzz=(double *)malloc(msp*sizeof(double));
  chge=(double *)malloc(msp*sizeof(double));
  space=(double *)malloc(mcore*sizeof(double));

  /* now call fortran dencor routine */  
  dencor_(&lden,&lcor,&ltim,&lchg,fname,&mcore,&msp,&natm,&ncon,
          &kmax,&ntime,&ngap,&numfil,xxx,yyy,zzz,chge,space,&status);

  exit(0);
  
}

