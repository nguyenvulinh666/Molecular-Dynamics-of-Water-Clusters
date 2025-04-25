/*
*******************************************************************************
*                                                                             *
*                                PLUMED                                       *
*   A Portable Plugin for Free Energy Calculations with Molecular Dynamics    *
*                              VERSION 1.2                                    *
*                                                                             *
*******************************************************************************
*
*  
*  Copyright (c) 2010 The PLUMED team.
*  See http://merlino.mi.infn.it/plumed for more information. 
*
*  This file is part of PLUMED.
*
*  PLUMED is free software: you can redistribute it and/or modify
*  it under the terms of the GNU Lesser General Public License as 
*  published by the Free Software Foundation, either version 3 of 
*  the License, or (at your option) any later version.
*
*  PLUMED is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
*  You should have received a copy of the GNU Lesser General
*  Public License along with PLUMED.  
*  If not, see <http://www.gnu.org/licenses/>.
*
*  For more info, see:  http://merlino.mi.infn.it/plumed
*  or subscribe to plumed-users@googlegroups.com
*
*/
#include "metadyn.h"

void PREFIX restraint(struct mtd_data_s *mtd_data)
{
  int i_c, ncv, nth, ntrh, ntp, rpxm, ntwg;                    // indexes for cycles and time

  hills.Vhills = cvw.Vwall = Vext = Vrecon = 0.;               // Hills,  Walls and external potential energy initialization

  colvar.it=mtd_data->istep;                                   // update microdynamics step
  if(colvar.it==mtd_data->istep_old){
    logical.not_same_step=0;
  } else {
    logical.not_same_step=1;
    mtd_data->istep_old=colvar.it;
  }

  ncv  = colvar.nconst;                                                                           	                  // number of CVs
  ntp  = (logical.not_same_step)&&(logical.print)&&(!(colvar.it%colvar.nt_print));			                  // have I got to print COLVAR?
  nth  = ( (logical.not_same_step)&&(logical.do_hills)&&(!(colvar.it%hills.nt_hills))&&(!firstTime) )
         || (hills.max_height>0);              									          // have I got to add HILL?
  ntrh = (logical.not_same_step)&&(logical.do_hills)&&(!(colvar.it%hills.nr_hills))&&(logical.do_walkers)&&(!firstTime) ; // period in steps to read HILLS
  ntwg = (logical.write_grid)&&(logical.not_same_step)&&(!(colvar.it%grid.w_stride))&&(!firstTime);                       // write GRID on file

  #ifdef GROMACS3
  set_pbc(&mtd_data->metapbc, mtd_data->cell);                                              // set pbc
  #endif
  #ifdef GROMACS4
  set_pbc(&mtd_data->metapbc, mtd_data->ePBC, mtd_data->cell);                                              // set pbc
  #endif

  if(logical.rpxm) rpxm = !(colvar.it%mtd_data->repl_ex_nst) && (logical.not_same_step);                    // have I got a replica exchange trial
  if(logical.debug){						// debug is a GROMACS global which identifies
    test_derivatives(mtd_data);		                       // the use of -debug options in mdrun
    #ifdef RECONMETAD
       if( reconOn==1 ){test_recon_derivatives(mtd_data);}    // Reconnaissance metadynamics version of test derivatives
    #endif  
    EXIT();			        // allow for only one step of dynamics
  }

// eventually align atoms
  if(colvar.align_atoms){
     int i,iatom1,iatom2;
     real distance[3],dummy;
     for(i=0;i<colvar.align_atoms-1;i++){
       iatom1=colvar.align_list[i];
       iatom2=colvar.align_list[i+1];
       minimal_image(mtd_data->pos[iatom2],mtd_data->pos[iatom1],&dummy,distance);
       mtd_data->pos[iatom2][0]=mtd_data->pos[iatom1][0]+distance[0];
       mtd_data->pos[iatom2][1]=mtd_data->pos[iatom1][1]+distance[1];
       mtd_data->pos[iatom2][2]=mtd_data->pos[iatom1][2]+distance[2];
     }
   };

  // this cycle is intended to calculate CVs values and derivatives

  for(i_c=0;i_c<ncv;i_c++){
    colvar.ff_hills[i_c] = 0.;                 	// initialization hills forces
    cvw.fwall[i_c]       = 0.;  		// initialization walls forces
    fext[i_c]            = 0.;                  // initialization external forces

//    if((!logical.always[i_c])&&(logical.rpxm)&&(!rpxm)&&(!ntp)) continue;
    if((!logical.always[i_c])&&(!ntp)) {   //if a variable is used only for printing purposes and it isn't the time to print then
      if(logical.rpxm) {
        if(!rpxm) continue;  //if we are doing bias-exchange and it isn't the time for an exchange we could avoid to calculate the variable
      } else continue; //if we are not doing bias-exchange we could avoid to calculate the variable
      // but if we are doing bias-exchange and it is time to try an axchange we must calculate all the variables
    }

    switch(colvar.type_s[i_c]){
      // geometric CVs
      case 1: dist_restraint(i_c, mtd_data); break;			// DISTANCE
      case 2: mindist_restraint(i_c, mtd_data); break;               	// MINDIST
      case 3: coord_restraint(i_c, mtd_data); break;	              	// COORD
      case 4: angle_restraint(i_c, mtd_data); break;	                // ANGLE
      case 5: torsion_restraint(i_c, mtd_data); break;                  // TORSION
      case 6: alfabeta_restraint(i_c, mtd_data); break;                	// ALPHA-BETA
      // interaction CVs
      case 7: hbonds_restraint(i_c, mtd_data); break;                   // HBONDS
      case 8: dipole_restraint(i_c, mtd_data); break;       		// DIPOLE
      // conformations CVs
      case 11: radgyr_restraint(i_c, mtd_data); break;   	       	// RGYR
      case 14: rmsdtor_restraint(i_c, mtd_data); break;	               	// RMSDTOR
      case 16: dihcor_restraint(i_c, mtd_data); break;                 	// DIHEDRAL-COR
      // water CVs
      case 20: waterbridge_restraint(i_c, mtd_data); break;            	// WATERBRIDGE
      // trajectory CVs
      case 30: spath_restraint(i_c, mtd_data); break;                   // S_MAPPATH
      case 31: zpath_restraint(i_c, mtd_data); break;                   // Z_MAPPATH
      case 32: position_restraint(i_c, mtd_data); break;                // ATOM POSITION
      case 33: elstpot_restraint(i_c, mtd_data); break;                 // ELSTPOT POSITION
      case 34: puckering_restraint(i_c, mtd_data); break;               // PUCKERING
      case 35: energy_restraint(i_c, mtd_data); break;                  // ENERGY
      case 36: helix_restraint(i_c, mtd_data); break;                   // HELIX
      case 37: alpharmsd_restraint(i_c, mtd_data); break;               // ALPHARMSD
      case 38: antibetarmsd_restraint(i_c, mtd_data); break;            // ANTIBETARMSD
      case 39: parabetarmsd_restraint(i_c, mtd_data); break;            // PARABETARMSD
      case 40: energy_restraint(i_c, mtd_data); break;                  // ENERGY P + S
    }

#ifdef PATHREF_FINDIFF
    fprintf(mtd_data->fplog,"|---CALLING THE TEST \n");
    //finite difference tests over reference frame derivatives
        if(colvar.type_s[i_c]==30){pathref_findiff(i_c,mtd_data);}
        if(colvar.type_s[i_c]==31){pathref_findiff(i_c,mtd_data);}
    fprintf(mtd_data->fplog,"|---END OF CALL \n");
    EXIT();
#endif

  }

  mtd_data->time=colvar.it*(mtd_data->dt)+mtd_data->time_offset;

  if(logical.commit) commit_analysis();     // committors analysis

  // this is the really dynamics code in which we calculate hills forces and then forces on CVs.
  if(logical.do_hills){
    hills_force();						// compute hills force and energy
    if(logical.widthadapt) hills_adapt();			// to adapt gaussian width
    if(nth) {
      hills_add(mtd_data);                                      // add HILL
    } 
    if(ntrh) {
       read_hills(mtd_data,0,hills.first_read);              	// is it time to read_hills?
       hills.first_read = 0;
    }
    if(ntwg)  grid_write_tofile(&grid);                         // write GRID on file
  }

  cvw.Vwall=soft_walls_engine(colvar.ss0,cvw.fwall);                // Wall potential

  Vext=ext_forces_engine(colvar.ss0,&extpot,fext);              // External potential

  cvw.Vwall+=steer_engine(colvar.ss0,cvw.fwall);                // Wall potential
  
  constraint_engine(mtd_data->dt);            // shake on the cvs 
  
  if(logical.do_steerplan)steerplan_engine();               // steerplan to plan adaptive us and more! 

  apply_forces(mtd_data);
//  for(i_c=0;i_c<ncv;i_c++) apply_force(i_c, mtd_data);          // add force coming from hills, restraint, walls... 

  if(logical.debug_derivatives) debug_derivatives(mtd_data,ntp);

#ifndef RECON_DRIVER                                            // Never do anything to colvar file if we are doing RECON_DRIVER
  if(ntp) print_colvar_enercv(mtd_data->time);	        	// dump COLVAR
#endif

  if(firstTime)firstTime = 0;			                // the first PluMed step 

  if(colvar.pg.nlist!=0)calc_projections( &(colvar.pg)); 
}

//------------------------------------------------------------------------------------------

real PREFIX soft_walls_engine(real* ss0,real* force){
  // Wall potential: V_wall(s) = sigma*((s-s0+offset)/redux)**n
  // WARNING: n must be even, sigma>0, redux>0; offset>0; s0 can be positive or negative
  real uscale,lscale,uexp,lexp;
  real V;
  int i;
  V=0.;
  for(i=0;i<colvar.nconst;i++){
    if(force) force[i]=0.0;
    if(logical.upper[i]) {                                                      // if there is a soft wall on this cv
      uscale = (ss0[i]-cvw.upper[i]+cvw.uoff[i])/cvw.ueps[i];                   // calculates the position on the wall 
      if(uscale>0.) {
        uexp = (real) cvw.uexp[i];
        V+=cvw.sigma[i]*pow(uscale, cvw.uexp[i]);
        if(force) force[i]+=-(cvw.sigma[i]/cvw.ueps[i])*uexp*pow(uscale, cvw.uexp[i]-1);
      }
    }
    if(logical.lower[i]) {                                                      // if there is a soft wall on this cv
      lscale = (ss0[i]-cvw.lower[i]-cvw.loff[i])/cvw.leps[i];                   // calculates the position on the wall
      if(lscale<0.) {
        lexp = (real) cvw.lexp[i];
        V+=cvw.lsigma[i]*pow(lscale, cvw.lexp[i]);
        if(force) force[i]+=-(cvw.lsigma[i]/cvw.leps[i])*lexp*pow(lscale, cvw.lexp[i]-1);
      }
    }
  };
  return V;
};
  
//-------------------------------------------------------------------------------------------

real PREFIX steer_engine(real* ss0,real* force)
{ 
  real V;
  int i_c;
  V=0.0;
  for(i_c=0;i_c<colvar.nconst;i_c++)if(logical.steer[i_c]){

  real tmp;
  if(firstTime){
   if ( cvsteer.impose_start[i_c] == 0 ) { 
          cvsteer.pos[i_c] = cvsteer.start[i_c] = ss0[i_c]; } 
   else {    cvsteer.pos[i_c] = cvsteer.start[i_c] ; }   

   fprintf(mtd_data.fplog,"|- STEERING %d CV : STARTVALUE %f\n",i_c+1,cvsteer.start[i_c]);
   fflush(mtd_data.fplog);
 
   if(cvsteer.max[i_c] < cvsteer.start[i_c]){
    cvsteer.sign[i_c] = -1; 
   } else {
    cvsteer.sign[i_c] = +1;
   } 
   // check if you're there since the beginning
   if(cvsteer.pos[i_c]==cvsteer.max[i_c]) {
          fprintf(mtd_data.fplog,"|- STEERING %d CV ARRIVED TO TARGET POINT %f in %d STEPS\n",i_c+1,cvsteer.max[i_c],colvar.it);
          fflush(mtd_data.fplog);
   } 
#ifdef STANDALONE
   FILE *file;
   char filename[100];
   sprintf(filename, "STEER.%d.rst",i_c);
   file = fopen(filename,"w");
   fprintf(file,"%lf  %lf",cvsteer.start[i_c],cvsteer.pos[i_c]);
   fclose(file);
#endif
  } 
  else{ // increase when you're not at the first step 
#ifdef STANDALONE 
        FILE *file;
        char *str, stringa[800],filename[100];
        // open the file 
        sprintf(filename, "STEER.%d.rst",i_c); 
        file = fopen(filename,"r");
        if(file==NULL){
          char buf[1024];
          sprintf(buf,"Cannot read %s  : EXITING\n",filename);
          plumed_error(buf);
        }else{
          str = fgets(stringa, 800, file); 
          sscanf(str, "%lf %lf",&cvsteer.start[i_c],&cvsteer.pos[i_c]);  
          fprintf(mtd_data.fplog,"|- STEERING %d CV RESTARTED FROM POINT: %f  STARTED: %f\n",i_c+1,cvsteer.pos[i_c],cvsteer.start[i_c]);  
          fclose(file);
          if(cvsteer.max[i_c] < cvsteer.start[i_c]){
            cvsteer.sign[i_c] = -1; 
          } else {
            cvsteer.sign[i_c] = +1;
          } 
        }
        if((logical.not_same_step) && fabs(cvsteer.pos[i_c]-cvsteer.start[i_c])<fabs(cvsteer.max[i_c]-cvsteer.start[i_c])) {
                cvsteer.pos[i_c] += cvsteer.sign[i_c] * cvsteer.delta[i_c] / 1000.0;
        }
        sprintf(filename, "STEER.%d.rst",i_c); 
        file = fopen(filename,"w");
        fprintf(file,"%lf  %lf",cvsteer.start[i_c],cvsteer.pos[i_c]);
        fclose(file); 
#else
      if((logical.not_same_step) && fabs(cvsteer.pos[i_c]-cvsteer.start[i_c])<fabs(cvsteer.max[i_c]-cvsteer.start[i_c])) {
              cvsteer.pos[i_c] += cvsteer.sign[i_c] * cvsteer.delta[i_c] / 1000.0;
      }
#endif
  } 
  if(fabs(cvsteer.pos[i_c]-cvsteer.start[i_c])>fabs(cvsteer.max[i_c]-cvsteer.start[i_c])) {
    cvsteer.pos[i_c] = cvsteer.max[i_c];
    fprintf(mtd_data.fplog,"|- STEERING %d CV ARRIVED TO TARGET POINT %f in %d STEPS\n",i_c+1,cvsteer.max[i_c],colvar.it);  
    fflush(mtd_data.fplog); 
  }
      /* HERE PUT THE PERIODICITY YOU NEED!!!!!!!! */ 
  tmp=ss0[i_c]-cvsteer.pos[i_c];
  if(colvar.type_s[i_c]==5 || ( colvar.type_s[i_c]==34 && colvar.type[i_c]==2 )) {
                   if(tmp > M_PI)
                     tmp -= 2.*M_PI;
                   if(tmp < -M_PI)
                    tmp += 2.*M_PI;
  } 
  real ff=-cvsteer.spring[i_c]*tmp-cvsteer.slope[i_c];
  V+=0.5*cvsteer.spring[i_c]*tmp*tmp+cvsteer.slope[i_c]*tmp;
  if(force) force[i_c] += ff;
  }
  return V;
}

//---------------------------------------------------------------------------------------------

void PREFIX abmd_cv(int i_c)
{
  real force;

  abmd.now[i_c] = (colvar.ss0[i_c]-abmd.exp[i_c])*(colvar.ss0[i_c]-abmd.exp[i_c]);
  if(abmd.now[i_c]<abmd.min[i_c]) abmd.min[i_c] = abmd.now[i_c];
  else {
    force = -2.*abmd.spring[i_c]*(abmd.now[i_c]-abmd.min[i_c])*(colvar.ss0[i_c]-abmd.exp[i_c]);
    cvw.fwall[i_c] += force;
    cvw.Vwall += 0.5*abmd.spring[i_c]*(abmd.now[i_c]-abmd.min[i_c])*(abmd.now[i_c]-abmd.min[i_c]);
  }
}

//---------------------------------------------------------------------------------------------

real PREFIX ext_forces_engine(real* ss0, struct grid_s *grid, real* force)
{
  real Vex = 0.0;

  if(logical.do_external) Vex=grid_getstuff(grid,ss0,force);

  return Vex;
}

//------------------------------------------------------------------------------

void PREFIX init_print_colvar_enercv()
{
/*
  In this routine we print the headers for the COLVAR file.
*/
  int i_c;
  FILE* cv_file;
  cv_file = fopen((mtd_data.ionode?mtd_data.colfilen:"/dev/null"), "a");
  fprintf(cv_file, "%s", "#! FIELDS");
  fprintf(cv_file, "%s", " time");
  for(i_c=0;i_c<colvar.nconst;i_c++) fprintf(cv_file, " cv%i", i_c+1);
  fprintf(cv_file, " vbias vwall vext");
  for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.steer[i_c]){fprintf(cv_file," XX XX RST%d",i_c+1);} }
  for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.abmd[i_c]){fprintf(cv_file," XX XX ABMD%d ",i_c+1);} }
  fprintf(cv_file,"\n");
  fclose(cv_file);
}

void PREFIX print_colvar_enercv(real time_s)
{
    int i_c; real consRecon=0.;
    static FILE *cv_file=NULL;

// This will allow us to build a conserved quantity for reconnaissance metadynamics
#ifdef RECONMETAD
    if( reconOn==1 ) consRecon=getCons_recon( myreconObj );
#endif

    if(!cv_file) cv_file = fopen((mtd_data.ionode?mtd_data.colfilen:"/dev/null"), "a");
/*
ATTENTION: all the quantities written here should be consistent with the header written in init_print_colvar_enercv()
In this way the new colvar format will work properly
To keep the new fmt equivalent to the old one (at least for a transition time)
we leave here also the "comments" such as "RST 3". These numbers are labeled with XX
in the header, and at some point will be removed.
*/
    fprintf(cv_file, "%10.3f", time_s);
    for(i_c=0;i_c<colvar.nconst;i_c++) fprintf(cv_file, "   %14.9f", colvar.ss0[i_c]);
    fprintf(cv_file, "   %14.9f   %14.9f   %14.9f   %14.9f   %14.9f", hills.Vhills/mtd_data.eunit, cvw.Vwall/mtd_data.eunit, Vext/mtd_data.eunit, Vrecon/mtd_data.eunit, consRecon/mtd_data.eunit );
    for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.steer[i_c]){fprintf(cv_file," RST %d %14.9f ",i_c+1,cvsteer.pos[i_c] );} }
    for(i_c=0;i_c<colvar.nconst;i_c++){  if(logical.abmd[i_c]){fprintf(cv_file," ABMD %d %14.9f ",i_c+1, abmd.min[i_c] );} }
    if(logical.do_steerplan)fprintf(cv_file,"%s",steerplan.log)  ;
    fprintf(cv_file, "\n");
    fflush(cv_file);
#ifdef STANDALONE 
    fclose(cv_file);
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------

void PREFIX commit_analysis()
{
  int i, a, b, ix;
  FILE *commit_file;

  if(firstTime){
    commit_file = fopen((mtd_data.ionode?"COMMIT":"/dev/null"), "a");
    for(i=0;i<colvar.nconst;i++) fprintf(commit_file, "   %14.7f", colvar.ss0[i]);
    fclose(commit_file);
  }

  a = 0;
  b = 0;
  for(i=0;i<commit.ncv;i++){
    ix=commit.index[i];
    if(commit.Amin[ix]<colvar.ss0[ix] && colvar.ss0[ix]<commit.Amax[ix]) a++;
    if(commit.Bmin[ix]<colvar.ss0[ix] && colvar.ss0[ix]<commit.Bmax[ix]) b++;
  }

  if(a==commit.ncv || b==commit.ncv) {
    fprintf(mtd_data.fplog, "|- SYSTEM HAS REACHED AN ENDING REGION\n");
    logical.commit = 0;
    commit_file = fopen((mtd_data.ionode?"COMMIT":"/dev/null"), "a");
    if(a==commit.ncv) fprintf(commit_file, " A \n");
    if(b==commit.ncv) fprintf(commit_file, " B \n");
    fclose(commit_file);
    EXIT();
  }
}

//----------------------------------------------------------------------------------------------------------------------
void PREFIX apply_forces(struct mtd_data_s *mtd_data )
{
  real ddr, uscale, lscale, dsdt, nor, fact;
  int ix, i, iat, wall,i_c;
  #ifdef NAMD
  Vector f;   
  #endif 

#ifdef RECONMETAD
   real recon_der[colvar.nconst];
   if( reconOn==1 ){
      double welikedupes_s[colvar.nconst], welikedupes_d[colvar.nconst];
      for(i=0;i<colvar.nconst;++i) welikedupes_s[i]=colvar.ss0[i];
      real ene; ene=dometa_recon(colvar.it, colvar.nconst, welikedupes_s, welikedupes_d, myreconObj);
      Vrecon+=ene;  //!TODO make sure this is giving expected results
      for(i=0;i<colvar.nconst;++i) recon_der[i]=welikedupes_d[i]; 
   } 
#endif

// set to zero all forces
  for(i=0;i<mtd_data->natoms;i++){
    mtd_data->force[i][0] = 0.0; 
    mtd_data->force[i][1] = 0.0;
    mtd_data->force[i][2] = 0.0;
  }  

  for(i_c=0;i_c<colvar.nconst;i_c++){

    if(logical.abmd[i_c]) abmd_cv(i_c);				
    
//    printf("In Restraint FORCE i_c colvar.ff_hills %d %f \n",i_c,colvar.ff_hills[i_c]);
    ddr = colvar.ff_hills[i_c] + cvw.fwall[i_c] + fext[i_c];	// hills, soft wall and external contribution
//    printf("In Restraint FORCE ddr %f \n",ddr);
#ifdef RECONMETAD
    if( reconOn==1 ) ddr+= recon_der[i_c];                      // Reconnaissance metadynamics forces
#endif

    if(logical.cnstr[i_c]) ddr-=cvcnstr.lambdadt2[i_c];         // constraint

    for(i=0;i<colvar.natoms[i_c];i++) {
      iat = colvar.cvatoms[i_c][i];
#ifdef NAMD 
       f.x=ddr*colvar.myder[i_c][i][0];
       f.y=ddr*colvar.myder[i_c][i][1];
       f.z=ddr*colvar.myder[i_c][i][2];
       addForce(iat,f);
#else
      mtd_data->force[iat][0] += ddr*colvar.myder[i_c][i][0];     // PluMeD forces
      mtd_data->force[iat][1] += ddr*colvar.myder[i_c][i][1];
      mtd_data->force[iat][2] += ddr*colvar.myder[i_c][i][2];
#endif    
    }
  }
}
//
// this engine SHAKEs  the needed d.o.f. (one by one, not interwined)  
// based on Frenkel's book
//
real PREFIX constraint_engine(real tstep){
  real V;
  int go,i,i_c,iat,niter,j;
  real ***posc,***newposc,***velc,***oldposc,***startder;
  FILE *fp;
  char *str,string[100]; 
  real **hack_pos1;
  real lambdadt2,tmp,totlambdadt2;
  char buf[200]; 
  // check whether constraining else skip
  go=0;
 
  for(i_c=0;i_c<colvar.nconst;i_c++){
     if(logical.cnstr[i_c]){go=1;}
  }
  if(go==0){return 0.;}

      //sprintf(buf,"|- ENTERING CONSTRAINING NST %d DT %f \n",mtd_data.istep,tstep );
      //plumed_warn(buf);
 
  // first step: do the mallocs 
#ifdef STANDALONE
    cvcnstr.oldforce=(real *)malloc(mtd_data.natoms*3*sizeof(real));
#endif

#ifndef STANDALONE
  if(mtd_data.istep==0){
#endif
    cvcnstr.newposc  =(real ***)malloc(colvar.nconst*sizeof(real **));
    cvcnstr.oldposc  =(real ***)malloc(colvar.nconst*sizeof(real **));
    cvcnstr.velc     =(real ***)malloc(colvar.nconst*sizeof(real **));
    cvcnstr.posc     =(real ***)malloc(colvar.nconst*sizeof(real **));
    cvcnstr.startder =(real ***)malloc(colvar.nconst*sizeof(real **));
    cvcnstr.go       =(int *)malloc(colvar.nconst*sizeof(int));
    cvcnstr.oldcv    =(real *)malloc(colvar.nconst*sizeof(real));
    for(i_c=0;i_c<colvar.nconst;i_c++){
       if(logical.cnstr[i_c]==1){
         cvcnstr.posc[i_c]       =float_2d_array_alloc(mtd_data.natoms,3); 
         cvcnstr.newposc[i_c]    =float_2d_array_alloc(mtd_data.natoms,3); 
         cvcnstr.oldposc[i_c]    =float_2d_array_alloc(mtd_data.natoms,3); 
         cvcnstr.velc[i_c]       =float_2d_array_alloc(mtd_data.natoms,3); 
         cvcnstr.startder[i_c]   =float_2d_array_alloc(mtd_data.natoms,3); 
       }   
    } 
#ifndef STANDALONE
  }
#endif
  // alias those names 
  posc   =cvcnstr.posc;
  newposc=cvcnstr.newposc;
  oldposc=cvcnstr.oldposc;
  velc=cvcnstr.velc;
  startder=cvcnstr.startder;
  // assign positions
  for(i_c=0;i_c<colvar.nconst;i_c++){
       if(logical.cnstr[i_c]==1){
          cvcnstr.go[i_c]=1;
          cvcnstr.oldcv[i_c]=colvar.ss0[i_c];
          for(i=0;i<colvar.natoms[i_c];i++){
               iat = colvar.cvatoms[i_c][i];
               posc[i_c][iat][0]=mtd_data.pos[iat][0]; 
               posc[i_c][iat][1]=mtd_data.pos[iat][1]; 
               posc[i_c][iat][2]=mtd_data.pos[iat][2]; 
               startder[i_c][i][0]=colvar.myder[i_c][i][0]; 
               startder[i_c][i][1]=colvar.myder[i_c][i][1]; 
               startder[i_c][i][2]=colvar.myder[i_c][i][2]; 
           //    fprintf(mtd_data.fplog,"|- %d D %f %f %f\n",i,colvar.myder[i_c][i][0],colvar.myder[i_c][i][1],colvar.myder[i_c][i][2]);
           //    fprintf(mtd_data.fplog,"|- %d D1 %f %f %f\n",i,startder[i_c][i][0],startder[i_c][i][1],startder[i_c][i][2]);
          }
       }
  }
  if(mtd_data.istep==0){
  // first step: put velocity to zero, put old position to the actual one 
    for(i_c=0;i_c<colvar.nconst;i_c++){
       if(logical.cnstr[i_c]==1){
              for(i=0;i<colvar.natoms[i_c];i++){
                   iat = colvar.cvatoms[i_c][i];
                   oldposc[i_c][iat][0]=posc[i_c][iat][0]; 
                   oldposc[i_c][iat][1]=posc[i_c][iat][1]; 
                   oldposc[i_c][iat][2]=posc[i_c][iat][2]; 
              }
       }
     }
  }
   
#ifdef STANDALONE
  // if standalone check if velocity is available from yout software
  fp = fopen("f.xyz","r");
//  sprintf(buf,"|- TAKING FORCE FILE \n");
//  plumed_warn(buf);
 
  if( fp ) {
    // exists
    i=0;
    while(1){
       str=fgets(string,100,fp);
       if(str==NULL)break;
       if(feof(fp))break;
   //    plumed_warn(str);
       str = strtok( string, " \t" );
       cvcnstr.oldforce[i*3+0]=atof(str);
       str=strtok(NULL," \t");
       cvcnstr.oldforce[i*3+1]=atof(str);
       str=strtok(NULL," \t");
       cvcnstr.oldforce[i*3+2]=atof(str);
       str=strtok(NULL," \t");
       i++; 
    } 
    if(i!=mtd_data.natoms){ plumed_error("|- FORCE FILE DOES NOT MATCH THE NUMBER OF ATOMS ");}
    fclose(fp);
  } else {
    plumed_error("|- NEED FOR FORCE FILE  \n");
  } 
//  sprintf(buf," TAKING VEL FILE \n");
//  plumed_warn(buf);
 
  fp = fopen("v.xyz","r");
  if( fp ) {
    // exists
    for(i_c=0;i_c<colvar.nconst;i_c++){
       if(logical.cnstr[i_c]==1){
          i=0;
          while(1){
             str=fgets(string,100,fp);
             if(str==NULL)break;
             if(feof(fp))break;
          //   plumed_warn(str);
             str = strtok( string, " \t" );
             cvcnstr.velc[i_c][i][0]=atof(str);
             str=strtok(NULL," \t");
             cvcnstr.velc[i_c][i][0]=atof(str);
             str=strtok(NULL," \t");
             cvcnstr.velc[i_c][i][0]=atof(str);
             str=strtok(NULL," \t");
             i++;
          } 
          if(i!=mtd_data.natoms){ plumed_error("|- VELOCITY FILE DOES NOT MATCH THE NUMBER OF ATOMS \n");}
          rewind(fp);
       } 
    }
    fclose(fp);
  } else {
    // doesnt exist ! old positions?
//    sprintf(buf," TAKING OLDPOS FILE \n");
//    plumed_warn(buf);
 
    fp = fopen("oldx.xyz","r");
    if( fp ) {
      for(i_c=0;i_c<colvar.nconst;i_c++){
          if(logical.cnstr[i_c]==1){
             i=0;
             while(1){
                str=fgets(string,100,fp);
                if(str==NULL)break;
                if(feof(fp))break;
             //   plumed_warn(str);
                str = strtok( string, " \t" );
                oldposc[i_c][i][0]=atof(str);
                str=strtok(NULL," \t");
                oldposc[i_c][i][0]=atof(str);
                str=strtok(NULL," \t");
                oldposc[i_c][i][0]=atof(str);
                str=strtok(NULL," \t");
                i++;
             } 
             if(i!=mtd_data.natoms){ plumed_error("|- VELOCITY FILE DOES NOT MATCH THE NUMBER OF ATOMS \n");}
          }
       }
       fclose(fp);
       for(i_c=0;i_c<colvar.nconst;i_c++){
           if(logical.cnstr[i_c]==1){
                  for(i=0;i<colvar.natoms[i_c];i++){
                       iat = colvar.cvatoms[i_c][i];
                       velc[i_c][iat][0]=(posc[i_c][iat][0]-oldposc[i_c][iat][0])/tstep; 
                       velc[i_c][iat][1]=(posc[i_c][iat][1]-oldposc[i_c][iat][1])/tstep; 
                       velc[i_c][iat][2]=(posc[i_c][iat][2]-oldposc[i_c][iat][2])/tstep; 
                       //fprintf(mtd_data.fplog,"|- %d V %f %f %f\n",i,velc[i_c][iat][0],velc[i_c][iat][1],velc[i_c][iat][2]);
                       //fprintf(mtd_data.fplog,"|- %d P %f %f %f\n",i,posc[i_c][iat][0],posc[i_c][iat][1],posc[i_c][iat][2]);
                       // save the old constraint derivative: pay attention on the derivative ordering
                  }
           }
       } 
    } else {
       plumed_error("|- CONSTRAINT NOT POSSIBLE WITHOUT EITHER OLD POSITION OR ACTUAL VELOCITIES \n");
    }
  }
#else
  // calculate velocity through finite difference 
  for(i_c=0;i_c<colvar.nconst;i_c++){
      if(logical.cnstr[i_c]==1){
             for(i=0;i<colvar.natoms[i_c];i++){
                  iat = colvar.cvatoms[i_c][i];
                  velc[i_c][iat][0]=(posc[i_c][iat][0]-oldposc[i_c][iat][0])/tstep; 
                  velc[i_c][iat][1]=(posc[i_c][iat][1]-oldposc[i_c][iat][1])/tstep; 
                  velc[i_c][iat][2]=(posc[i_c][iat][2]-oldposc[i_c][iat][2])/tstep; 
                //  fprintf(mtd_data.fplog,"|- %d V %f %f %f\n",i,velc[i_c][iat][0],velc[i_c][iat][1],velc[i_c][iat][2]);
                //  fprintf(mtd_data.fplog,"|- %d P %f %f %f\n",i,posc[i_c][iat][0],posc[i_c][iat][1],posc[i_c][iat][2]);
                  // save the old constraint derivative: pay attention on the derivative ordering
             }
      }
  }  
#endif
  for(i_c=0;i_c<colvar.nconst;i_c++){
     if(logical.cnstr[i_c]==1){

           if(cvcnstr.verbose[i_c]==1){sprintf(buf," DOING CONSTRAINTS ON CV %d ",i_c+1);
           plumed_warn(buf);}
           niter=0;
           totlambdadt2=0.;
           while(cvcnstr.go[i_c]){// check on threshold value
              tmp=2*cvcnstr.spring[i_c]*(cvcnstr.oldcv[i_c]-cvcnstr.pos[i_c]);
              for(i=0;i<colvar.natoms[i_c];i++){
                  iat = colvar.cvatoms[i_c][i];
                  // make a velocity verlet free evolution (assume the units are congruent)
                  //  r_n= r_o +v dt - 0.5 * dt**2 * F
                  newposc[i_c][iat][0]=posc[i_c][iat][0]+velc[i_c][iat][0]*tstep+ 0.5 * tstep *tstep * 
#if defined (AMBER)  || defined (STANDALONE)
                  (cvcnstr.oldforce[iat*3+0])    /   mtd_data.mass[iat]       
#else
           sprintf(buf," CONSTRAINTS NOT IMPLEMENTED IN THIS PROGRAM ");
           plumed_error(buf);
#endif            
                  - totlambdadt2*tmp*startder[i_c][i][0]/mtd_data.mass[iat];
                  
                  newposc[i_c][iat][1]=posc[i_c][iat][1]+velc[i_c][iat][1]*tstep+ 0.5 * tstep *tstep * 
#if defined (AMBER)  || defined (STANDALONE)
                  (cvcnstr.oldforce[iat*3+1])    /   mtd_data.mass[iat]       
#else
           sprintf(buf," CONSTRAINTS NOT IMPLEMENTED IN THIS PROGRAM ");
           plumed_error(buf);
#endif 
                  - totlambdadt2*tmp*startder[i_c][i][1]/mtd_data.mass[iat];
                  
                  newposc[i_c][iat][2]=posc[i_c][iat][2]+velc[i_c][iat][2]*tstep+ 0.5 * tstep *tstep * 
#if defined (AMBER)  || defined (STANDALONE)
                  (cvcnstr.oldforce[iat*3+2])   /    mtd_data.mass[iat]         
#else
           sprintf(buf," CONSTRAINTS NOT IMPLEMENTED IN THIS PROGRAM ");
           plumed_error(buf);
#endif 
                  - totlambdadt2*tmp*startder[i_c][i][2]/mtd_data.mass[iat];
              }
              //calculate the constraint value
              
              // do an hack : unlock the address from mtd_data.pos and redirict to an identical vector newposc[i_c]    
              hack_pos1=mtd_data.pos; // remember this address
              mtd_data.pos=newposc[i_c]; //calculate now with this coordinates 
              //fprintf(mtd_data.fplog,"|- PREVIOUS CV %f\n",colvar.ss0[i_c] ); 
              switch(colvar.type_s[i_c]){
                  // geometric CVs
                  case 1: dist_restraint(i_c, &mtd_data); break;			// DISTANCE
                  case 2: mindist_restraint(i_c, &mtd_data); break;               	// MINDIST
                  case 3: coord_restraint(i_c, &mtd_data); break;	              	// COORD
                  case 4: angle_restraint(i_c, &mtd_data); break;	                // ANGLE
                  case 5: torsion_restraint(i_c, &mtd_data); break;                  // TORSION
                  case 6: alfabeta_restraint(i_c, &mtd_data); break;                	// ALPHA-BETA
                  // interaction CVs
                  case 7: hbonds_restraint(i_c, &mtd_data); break;                   // HBONDS
                  case 8: dipole_restraint(i_c, &mtd_data); break;       		// DIPOLE
                  // conformations CVs
                  case 11: radgyr_restraint(i_c, &mtd_data); break;   	       	// RGYR
                  case 14: rmsdtor_restraint(i_c, &mtd_data); break;	               	// RMSDTOR
                  case 16: dihcor_restraint(i_c, &mtd_data); break;                 	// DIHEDRAL-COR
                  // water CVs
                  case 20: waterbridge_restraint(i_c, &mtd_data); break;            	// WATERBRIDGE
                  // trajectory CVs
                  case 30: spath_restraint(i_c, &mtd_data); break;                   // S_MAPPATH
                  case 31: zpath_restraint(i_c, &mtd_data); break;                   // Z_MAPPATH
                  case 32: position_restraint(i_c, &mtd_data); break;                // ATOM POSITION
                  case 33: elstpot_restraint(i_c, &mtd_data); break;                 // ELSTPOT POSITION
                  case 34: puckering_restraint(i_c, &mtd_data); break;               // PUCKERING
                  case 35: energy_restraint(i_c, &mtd_data); break;                  // ENERGY
                  case 36: helix_restraint(i_c, &mtd_data); break;                   // HELIX
                  case 37: alpharmsd_restraint(i_c, &mtd_data); break;               // ALPHARMSD
                  case 38: antibetarmsd_restraint(i_c, &mtd_data); break;            // ANTIBETARMSD
                  case 39: parabetarmsd_restraint(i_c, &mtd_data); break;            // PARABETARMSD
                  case 40: energy_restraint(i_c, &mtd_data); break;                //   ENERGYPS
              }
              mtd_data.pos=hack_pos1;  
                
              //calculate lambda*deltat*deltat 
              tmp=0.;
              for(i=0;i<colvar.natoms[i_c];i++){
                 tmp+=startder[i_c][i][0]*colvar.myder[i_c][i][0];
                 tmp+=startder[i_c][i][1]*colvar.myder[i_c][i][1];
                 tmp+=startder[i_c][i][2]*colvar.myder[i_c][i][2];
              }   
              // calculate new lambda
              lambdadt2=(0.25/cvcnstr.spring[i_c])*(colvar.ss0[i_c]-cvcnstr.pos[i_c])/((cvcnstr.oldcv[i_c]-cvcnstr.pos[i_c])*tmp);
              // apply it  
              totlambdadt2+=lambdadt2;
              // die?
              tmp=sqrt((colvar.ss0[i_c]-cvcnstr.pos[i_c])*(colvar.ss0[i_c]-cvcnstr.pos[i_c]));  
              //fprintf(mtd_data.fplog,"|- totlambdadt2 %f lambdadt2 %f tmp %f POS %f CV %f NITER %d\n",totlambdadt2,lambdadt2,tmp,cvcnstr.pos[i_c],colvar.ss0[i_c],niter);
              if(tmp<cvcnstr.delta[i_c])cvcnstr.go[i_c]=0;
              niter++;
              if(niter>cvcnstr.maxiter[i_c])cvcnstr.go[i_c]=0;  
           }    
           cvcnstr.lambdadt2[i_c]=totlambdadt2*cvcnstr.spring[i_c]*(cvcnstr.oldcv[i_c]-cvcnstr.pos[i_c])/(tstep*tstep); 
           if(cvcnstr.verbose[i_c]==1){sprintf(buf," DONE CONSTRAINT ON CV %d LAMBDA %f ITER %d",i_c+1,cvcnstr.lambdadt2[i_c],niter);
           plumed_warn(buf);}
     }
  }
  // put the actual position as the old position
  for(i_c=0;i_c<colvar.nconst;i_c++){
     if(logical.cnstr[i_c]==1){
        colvar.ss0[i_c]=cvcnstr.oldcv[i_c];
        for(i=0;i<colvar.natoms[i_c];i++){
             iat = colvar.cvatoms[i_c][i];
             oldposc[i_c][iat][0]=posc[i_c][iat][0]; 
             oldposc[i_c][iat][1]=posc[i_c][iat][1]; 
             oldposc[i_c][iat][2]=posc[i_c][iat][2]; 
             colvar.myder[i_c][i][0]=startder[i_c][i][0];
             colvar.myder[i_c][i][1]=startder[i_c][i][1];
             colvar.myder[i_c][i][2]=startder[i_c][i][2];
 
        }
     }
  }
 
  
      //sprintf(buf,"|- EXITING CONSTRAINING \n");
      //plumed_warn(buf);
 
  V=0.;
//  EXIT();
  return V;
};
 
void PREFIX steerplan_engine()
{ 
   int i,j,k,l,m,thisstage,nextstage,mycv; 
   real time,x1,x2,k1,k2,deltat,tmp,force;
   time=mtd_data.time;
   //fprintf(mtd_data.fplog,"|- STEERPLAN_ON TIME %f FIRSTTIME %d\n",time,firstTime);

   //char string2[200],string[200];
   //for(i=0;i<steerplan.totstages;i++){
   //    fprintf(mtd_data.fplog,"|- TIME: %12.6f  ",steerplan.actions[i].t);  
   //    for(j=0;j<steerplan.ncvs;j++){
   //      if(steerplan.actions[i].activecv[j].wildcardpos){
   //       sprintf(string,"*");
   //      }else {
   //       sprintf(string,"%12.6f",steerplan.actions[i].activecv[j].pos);
   //      }
   //      if(steerplan.actions[i].activecv[j].wildcardk){
   //       sprintf(string2,"*");
   //      }else {
   //       sprintf(string2,"%12.6f",steerplan.actions[i].activecv[j].k);
   //      }
   //      if(steerplan.actions[i].activecv[j].k<0){
   //           fprintf(mtd_data.fplog," CV %3d SKIPPING THIS STAGE... ",(steerplan.actions[i].activecv[j].ncv+1));
   //      }else{
   //           fprintf(mtd_data.fplog," CV %3d KAPPA %12s POS %12s ",steerplan.actions[i].activecv[j].ncv+1,string2,string);
   //      }
   //    }
   //    fprintf(mtd_data.fplog,"\n");  
   //}
 
   // firsttime: define the initial stage and the time for the nextone 
   // which stage I am? New one? calculate new params
#ifdef STANDALONE 
        FILE *file;
        char filename[100],sstr[100] ;
        char *str, stringa[800];
        char buf[1024];
        // open the file 
        if(!firstTime){
          sprintf(filename, "STEERPLAN.rst"); 
          file = fopen(filename,"r");
          if(file==NULL){
            sprintf(buf,"Cannot read %s  : EXITING\n",filename);
            plumed_error(buf);
          }else{
            str = fgets(sstr, 100, file); 
            sscanf(str,"%lf",&steerplan.nextstage_time);
            for(i=0;i<steerplan.ncvs;i++){
                 str = fgets(sstr, 100, file); 
                 sscanf(sstr,"%lf %lf %lf %lf %d",&steerplan.actualcv[i].k,&steerplan.actualcv[i].kv,&steerplan.actualcv[i].x0,&steerplan.actualcv[i].v,&steerplan.actualcv[i].type);
            }
            fclose(file);
          }
        }
#endif
   if(time>=steerplan.nextstage_time||firstTime){
       for(i=0;i<steerplan.totstages;i++){
           if(steerplan.actions[i].t>=time){steerplan.current_stage=i;
              if((i+1)!=steerplan.totstages){steerplan.nextstage_time=steerplan.actions[i+1].t;}
              else{steerplan.nextstage_time=1.e9;} 
           break;} 
       }
       fprintf(mtd_data.fplog,"|- CURRENT STAGE %d  \n",steerplan.current_stage);
       // setup the run
       thisstage=steerplan.current_stage;
       for(i=0;i<steerplan.ncvs;i++){
          mycv=steerplan.actions[thisstage].activecv[i].ncv;
          // is this a skipped stage for this cv?
          if (steerplan.actions[thisstage].activecv[i].k<0){// bypass and keep the values 
                 fprintf(mtd_data.fplog,"|- STEERPLAN CV  %d :NOTHING TO UPDATE TO AT THIS STAGE. I'LL GO ON AS BEFORE\n",steerplan.actions[thisstage].activecv[i].ncv+1);
               // if firstime you don't specify any action this is taken as k=0 pos=* type=1 v=0.  action
               // case1: coming  from nowhere: default is put k=0 (no effective potential)
              if(firstTime){
                       steerplan.actualcv[i].k=0.;
                       steerplan.actualcv[i].v=0.;
                       steerplan.actualcv[i].kv=0.;
                       steerplan.actualcv[i].x0=colvar.ss0[mycv];
                       steerplan.actualcv[i].type=1;
                       // keeep nowhere
              }
          }else{ // scan to the next stages to find the following pinpoints if any
                 fprintf(mtd_data.fplog,"|- STEERPLAN CV  %d SCANNING FOR THE NEXT ACTION\n",steerplan.actions[thisstage].activecv[i].ncv+1);
                 // case2: going to  nowhere: default is keep the old restraint (static: v=0.)
                 // this happens when the list is over or when all the k are <0  
                 if(thisstage==steerplan.totstages-1){
             //    fprintf(mtd_data.fplog,"1 THISSTAGE %d TOTSTAGES %d\n",thisstage,steerplan.totstages);
                       if(!steerplan.actions[thisstage].activecv[i].wildcardk){steerplan.actualcv[i].k=steerplan.actions[thisstage].activecv[i].k;}
                       steerplan.actualcv[i].v=0.;
                       steerplan.actualcv[i].kv=0.;
                       if(steerplan.actions[thisstage].activecv[i].wildcardpos){
                                steerplan.actualcv[i].x0=colvar.ss0[mycv]; 
                       }
                       else{ steerplan.actualcv[i].x0=steerplan.actions[thisstage].activecv[i].pos;}
                       steerplan.actualcv[i].type=steerplan.actions[thisstage].activecv[i].type;
                 }else{ 
                       for(j=thisstage+1;j<steerplan.totstages;j++){
                    //       fprintf(mtd_data.fplog,"XXCV %d K %f WC %d\n",mycv+1,steerplan.actions[j].activecv[i].k,steerplan.actions[j].activecv[i].wildcardk);
                           if((steerplan.actions[j].activecv[i].k>=0.) || (steerplan.actions[j].activecv[i].wildcardk==1)){
                           break;} // got it!!!!
                       } 
                       if(j==steerplan.totstages){//no more actions to take! 
                           if(!steerplan.actions[thisstage].activecv[i].wildcardk){steerplan.actualcv[i].k=steerplan.actions[thisstage].activecv[i].k;}
                           steerplan.actualcv[i].v=0.;
                           steerplan.actualcv[i].kv=0.;
                           if(steerplan.actions[thisstage].activecv[i].wildcardpos){
                               steerplan.actualcv[i].x0=colvar.ss0[mycv];
                           }
                           else{ steerplan.actualcv[i].x0=steerplan.actions[thisstage].activecv[i].pos;}
                           steerplan.actualcv[i].type=steerplan.actions[thisstage].activecv[i].type;
                    //       fprintf(mtd_data.fplog,"2 THISSTAGE %d TOTSTAGES %d\n",thisstage,steerplan.totstages);
                       }else{ 
                    //      fprintf(mtd_data.fplog,"3 THISSTAGE %d TOTSTAGES %d\n",thisstage,steerplan.totstages);
                          deltat=(steerplan.actions[j].t-steerplan.actions[thisstage].t)/mtd_data.dt;              
                          // copy the actual center (wildcard?) : start from the actual real cv 
                          if(steerplan.actions[thisstage].activecv[i].wildcardpos){steerplan.actualcv[i].x0=colvar.ss0[mycv];}
                          else{ steerplan.actualcv[i].x0=steerplan.actions[thisstage].activecv[i].pos;}
                          // copy the next center (wildcard?) : start from the previous center!!!
                          if(steerplan.actions[j].activecv[i].wildcardpos){
                              steerplan.actualcv[i].v=0.;   
                              steerplan.actualcv[i].nowhere=1;
                          }
                          else { 
                             x1=steerplan.actualcv[i].x0; 
                             x2=steerplan.actions[j].activecv[i].pos  ;
                             steerplan.actualcv[i].v=(x2-x1)/deltat   ;
                             steerplan.actualcv[i].nowhere=0;
                          };
 
                           // copy the actual k (wildcard? dont do anything ) 
                          if(!steerplan.actions[thisstage].activecv[i].wildcardk){steerplan.actualcv[i].k=steerplan.actions[thisstage].activecv[i].k;}
                          // copy the next k (wildcard?) : keep this one 
                          if(steerplan.actions[j].activecv[i].wildcardk){k2=steerplan.actions[thisstage].activecv[i].k;}
                          else { k2=steerplan.actions[j].activecv[i].k  ;};
                          k1=steerplan.actualcv[i].k;
                          steerplan.actualcv[i].kv=(k2-k1)/deltat   ;
                          steerplan.actualcv[i].type=steerplan.actions[thisstage].activecv[i].type;
                     //     fprintf(mtd_data.fplog,"XXX CV %d X1 %f X2 %f K1 %f K2 %f \n",steerplan.actions[thisstage].activecv[i].ncv+1,x1,x2,k1,k2); 
                       } 
                 } 
                  
          }
       }
   }else{// evolove normally
      for(i=0;i<steerplan.ncvs;i++){
            steerplan.actualcv[i].k=steerplan.actualcv[i].k+steerplan.actualcv[i].kv;
            steerplan.actualcv[i].x0=steerplan.actualcv[i].x0+steerplan.actualcv[i].v;
            if(steerplan.actualcv[i].nowhere){
                     mycv=steerplan.actions[steerplan.current_stage].activecv[i].ncv;
                     steerplan.actualcv[i].x0=colvar.ss0[mycv]; 
            }

      }      
   }
#ifdef STANDALONE
       sprintf(filename,"STEERPLAN.rst");
       file = fopen(filename,"w");
       fprintf(file,"%lf\n",steerplan.nextstage_time);
       for(i=0;i<steerplan.ncvs;i++){
            fprintf(file,"%lf %lf %lf %lf %d\n",steerplan.actualcv[i].k,steerplan.actualcv[i].kv,steerplan.actualcv[i].x0,steerplan.actualcv[i].v,steerplan.actualcv[i].type);
       }
       fclose(file);
#endif

   // create a log string
   strcpy(steerplan.log," STP ");
   char cvlog[100]; 
   for(i=0;i<steerplan.ncvs;i++){
          mycv=steerplan.actions[steerplan.current_stage].activecv[i].ncv;
          //fprintf(mtd_data.fplog,"|- STEERPLAN CV  %d : \n",steerplan.actions[steerplan.current_stage].activecv[i].ncv+1);
          //fprintf(mtd_data.fplog,"|- X  %12.6f XV %12.6f  K %12.6f KV %12.6f  \n",steerplan.actualcv[i].x0,steerplan.actualcv[i].v,steerplan.actualcv[i].k,steerplan.actualcv[i].kv);
          sprintf(cvlog," CV %2d X %12.6f K %12.6f T %1d ",mycv+1,steerplan.actualcv[i].x0,steerplan.actualcv[i].k,steerplan.actualcv[i].type);
          strcat(steerplan.log,cvlog);
   }  
   //fprintf(mtd_data.fplog,"%s\n",steerplan.log);

    // now the plan is fine: do what you need!!

   for(i=0;i<steerplan.ncvs;i++){
      mycv=steerplan.actions[steerplan.current_stage].activecv[i].ncv; 
      tmp=colvar.ss0[mycv]-steerplan.actualcv[i].x0;
      if(steerplan.actualcv[i].type==2){
         // positive: wall only if positive
         if(tmp<0.)tmp=0.;
      } else if (steerplan.actualcv[i].type==3){
         // negative: wall only if negative 
         if(tmp>0.)tmp=0.;
      } 
      if(colvar.type_s[mycv]==5 || ( colvar.type_s[mycv]==34 && colvar.type[mycv]==2 )){
                       if(tmp > M_PI)
                         tmp -= 2.*M_PI;
                       if(tmp < -M_PI)
                        tmp += 2.*M_PI;
      } 
      force = -steerplan.actualcv[i].k*tmp;
      cvw.fwall[mycv] += force;
      cvw.Vwall += -0.5*force*tmp; 
   } 

//   fprintf(mtd_data.fplog,"|- STEERPLAN_ON END\n");
}

//---------------------------------------------------------------------------------------------


