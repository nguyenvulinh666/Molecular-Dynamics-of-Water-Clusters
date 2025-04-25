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

  void PREFIX test_derivatives(struct mtd_data_s *mtd_data)
{
  rvec testforce;
  real teststep = 0.0001, invstep = 10000., val1, val2, diff, analder;
  int i_c, ix, i, j, iat;

  for(i_c=0;i_c<colvar.nconst;i_c++) {
    for(i=0;i<colvar.natoms[i_c];i++) {
      iat = colvar.cvatoms[i_c][i];
      for(ix=0;ix<3;ix++) {
        mtd_data->pos[iat][ix] += teststep;
        switch(colvar.type_s[i_c]){
          case 1: dist_restraint(i_c, mtd_data); break;		
          case 2: mindist_restraint(i_c, mtd_data); break;                   
          case 3: coord_restraint(i_c, mtd_data); break;                     
          case 4: angle_restraint(i_c, mtd_data); break;                     
          case 5: torsion_restraint(i_c, mtd_data); break;           
          case 6: alfabeta_restraint(i_c, mtd_data); break;                  
          case 7: hbonds_restraint(i_c, mtd_data); break;                    
          case 8: dipole_restraint(i_c, mtd_data); break;            
          case 11: radgyr_restraint(i_c, mtd_data); break;           
          case 14: rmsdtor_restraint(i_c, mtd_data); break;                  
          case 16: dihcor_restraint(i_c, mtd_data); break;
          case 20: waterbridge_restraint(i_c, mtd_data); break;              
          case 30: spath_restraint(i_c, mtd_data); break;
          case 31: zpath_restraint(i_c, mtd_data); break;
          case 32: position_restraint(i_c, mtd_data); break;
          case 33: elstpot_restraint(i_c, mtd_data); break;
          case 34: puckering_restraint(i_c, mtd_data); break;
          case 36: helix_restraint(i_c, mtd_data); break;
          case 37: alpharmsd_restraint(i_c, mtd_data); break;
          case 38: antibetarmsd_restraint(i_c, mtd_data); break;
          case 39: parabetarmsd_restraint(i_c, mtd_data); break;
          //case 40: camshift_restraint(i_c, mtd_data); break;
        }
        val1 = colvar.ss0[i_c];
        mtd_data->pos[iat][ix] += -2.*teststep;
        switch(colvar.type_s[i_c]){
          case 1: dist_restraint(i_c, mtd_data); break;
          case 2: mindist_restraint(i_c, mtd_data); break;           
          case 3: coord_restraint(i_c, mtd_data); break;                
          case 4: angle_restraint(i_c, mtd_data); break;                 
          case 5: torsion_restraint(i_c, mtd_data); break;              
          case 6: alfabeta_restraint(i_c, mtd_data); break;             
          case 7: hbonds_restraint(i_c, mtd_data); break;               
          case 8: dipole_restraint(i_c, mtd_data); break;               
          case 11: radgyr_restraint(i_c, mtd_data); break;              
          case 14: rmsdtor_restraint(i_c, mtd_data); break;          
          case 16: dihcor_restraint(i_c, mtd_data); break;              
          case 20: waterbridge_restraint(i_c, mtd_data); break;      
          case 30: spath_restraint(i_c, mtd_data); break;
          case 31: zpath_restraint(i_c, mtd_data); break;
          case 32: position_restraint(i_c, mtd_data); break;
          case 33: elstpot_restraint(i_c, mtd_data); break;
          case 34: puckering_restraint(i_c, mtd_data); break;
          case 36: helix_restraint(i_c, mtd_data); break;
          case 37: alpharmsd_restraint(i_c, mtd_data); break;
          case 38: antibetarmsd_restraint(i_c, mtd_data); break;
          case 39: parabetarmsd_restraint(i_c, mtd_data); break;
          //case 40: camshift_restraint(i_c, mtd_data); break;
        }
        val2 = colvar.ss0[i_c];
        testforce[ix] = 0.5*((val1*invstep)-(val2*invstep));
        if(colvar.type_s[i]==5 || ( colvar.type_s[i]==34 && colvar.type[i]==2 )) {
           diff = val1-val2;
           if(diff>M_PI) diff-=2.*M_PI;
           else if(diff<-M_PI) diff+=2.*M_PI;
           testforce[ix] = 0.5*(diff*invstep);
        }
        mtd_data->pos[iat][ix] += teststep;
        switch(colvar.type_s[i_c]){
          case 1: dist_restraint(i_c, mtd_data); break;
          case 2: mindist_restraint(i_c, mtd_data); break;                               
          case 3: coord_restraint(i_c, mtd_data); break;                
          case 4: angle_restraint(i_c, mtd_data); break;                
          case 5: torsion_restraint(i_c, mtd_data); break;              
          case 6: alfabeta_restraint(i_c, mtd_data); break;             
          case 7: hbonds_restraint(i_c, mtd_data); break;               
          case 8: dipole_restraint(i_c, mtd_data); break;               
          case 11: radgyr_restraint(i_c, mtd_data); break;              
          case 14: rmsdtor_restraint(i_c, mtd_data); break;                              
          case 16: dihcor_restraint(i_c, mtd_data); break;              
          case 20: waterbridge_restraint(i_c, mtd_data); break;                          
          case 30: spath_restraint(i_c, mtd_data); break;
          case 31: zpath_restraint(i_c, mtd_data); break;
          case 32: position_restraint(i_c, mtd_data); break;
          case 33: elstpot_restraint(i_c, mtd_data); break;
          case 34: puckering_restraint(i_c, mtd_data); break;
          case 36: helix_restraint(i_c, mtd_data); break;
          case 37: alpharmsd_restraint(i_c, mtd_data); break;
          case 38: antibetarmsd_restraint(i_c, mtd_data); break;
          case 39: parabetarmsd_restraint(i_c, mtd_data); break;
          //case 40: camshift_restraint(i_c, mtd_data); break;
        }
// analytic derivative: the same atom may be given in input more than once
        analder=0.;
        for(j=0;j<colvar.natoms[i_c];j++) {
          if(colvar.cvatoms[i_c][j]==iat) {
            analder+=colvar.myder[i_c][j][ix];
          }
        }
        printf("atom %5i[%i] ** analytic %15.10f ** numeric %15.10f *** DELTA %15.10f\n",
                iat, ix, analder, testforce[ix], analder-testforce[ix]);
      }  
    }
  }
}

#ifdef RECONMETAD  
void PREFIX test_recon_derivatives(struct mtd_data_s *mtd_data)
{
  rvec testforce;
  real teststep = 0.0001, invstep = 10000., diff, analder; // Qanalder;
  double val1, val2 ;
  int i_c,i_oc, ix, i, j, iat, it;
  double welikedupes_s[colvar.nconst], recon_der[colvar.nconst];

  // Initialize recon metadynamics plugin for testing
  for(i_c=0;i_c<colvar.nconst;i_c++) { // Loop over all cvs
      switch(colvar.type_s[i_c]){
        case 1: dist_restraint(i_c, mtd_data); break;
        case 2: mindist_restraint(i_c, mtd_data); break;
        case 3: coord_restraint(i_c, mtd_data); break;
        case 4: angle_restraint(i_c, mtd_data); break;
        case 5: torsion_restraint(i_c, mtd_data); break;
        case 6: alfabeta_restraint(i_c, mtd_data); break;
        case 7: hbonds_restraint(i_c, mtd_data); break;
        case 8: dipole_restraint(i_c, mtd_data); break;
        case 11: radgyr_restraint(i_c, mtd_data); break;
        case 14: rmsdtor_restraint(i_c, mtd_data); break;
        case 16: dihcor_restraint(i_c, mtd_data); break;
        case 20: waterbridge_restraint(i_c, mtd_data); break;
        case 30: spath_restraint(i_c, mtd_data); break;
        case 31: zpath_restraint(i_c, mtd_data); break;
        case 32: position_restraint(i_c, mtd_data); break;
        case 33: elstpot_restraint(i_c, mtd_data); break;
        case 34: puckering_restraint(i_c, mtd_data); break;
        case 36: helix_restraint(i_c, mtd_data); break;
        case 37: alpharmsd_restraint(i_c, mtd_data); break;
        case 38: antibetarmsd_restraint(i_c, mtd_data); break;
        case 39: parabetarmsd_restraint(i_c, mtd_data); break;
        //case 40: camshift_restraint(i_c, mtd_data); break;
      }
  }

  for(it=0;it<colvar.nconst;++it) welikedupes_s[it]=colvar.ss0[it];
  setup_test_recon( colvar.nconst, welikedupes_s, myreconObj );

  for(i_oc=0;i_oc<colvar.nconst;i_oc++) {    // Loop over all atoms in all cvs
      for(i=0;i<colvar.natoms[i_oc];i++) {   //  "    "    "   "    "   "   "
          iat = colvar.cvatoms[i_oc][i];
          for(ix=0;ix<3;ix++) {              // Loop over all coordinates 
              mtd_data->pos[iat][ix] += teststep;
              for(i_c=0;i_c<colvar.nconst;i_c++) { // Loop over all cvs
                  switch(colvar.type_s[i_c]){
                    case 1: dist_restraint(i_c, mtd_data); break;
                    case 2: mindist_restraint(i_c, mtd_data); break;
                    case 3: coord_restraint(i_c, mtd_data); break;
                    case 4: angle_restraint(i_c, mtd_data); break;
                    case 5: torsion_restraint(i_c, mtd_data); break;
                    case 6: alfabeta_restraint(i_c, mtd_data); break;
                    case 7: hbonds_restraint(i_c, mtd_data); break;
                    case 8: dipole_restraint(i_c, mtd_data); break;
                    case 11: radgyr_restraint(i_c, mtd_data); break;
                    case 14: rmsdtor_restraint(i_c, mtd_data); break;
                    case 16: dihcor_restraint(i_c, mtd_data); break;
                    case 20: waterbridge_restraint(i_c, mtd_data); break;
                    case 30: spath_restraint(i_c, mtd_data); break;
                    case 31: zpath_restraint(i_c, mtd_data); break;
                    case 32: position_restraint(i_c, mtd_data); break;
                    case 33: elstpot_restraint(i_c, mtd_data); break;
                    case 34: puckering_restraint(i_c, mtd_data); break;
                    case 36: helix_restraint(i_c, mtd_data); break;
                    case 37: alpharmsd_restraint(i_c, mtd_data); break;
                    case 38: antibetarmsd_restraint(i_c, mtd_data); break;
                    case 39: parabetarmsd_restraint(i_c, mtd_data); break;
                    //case 40: camshift_restraint(i_c, mtd_data); break;
                  }
              }
              // Now P and Q using recon plugin
              for(it=0;it<colvar.nconst;++it) welikedupes_s[it]=colvar.ss0[it];
              recon_calc_cv( colvar.nconst, welikedupes_s, &val1, myreconObj);

              mtd_data->pos[iat][ix] += -2.*teststep;
              for(i_c=0;i_c<colvar.nconst;i_c++) {
                  switch(colvar.type_s[i_c]){
                   case 1: dist_restraint(i_c, mtd_data); break;
                   case 2: mindist_restraint(i_c, mtd_data); break;
                   case 3: coord_restraint(i_c, mtd_data); break;
                   case 4: angle_restraint(i_c, mtd_data); break;
                   case 5: torsion_restraint(i_c, mtd_data); break;
                   case 6: alfabeta_restraint(i_c, mtd_data); break;
                   case 7: hbonds_restraint(i_c, mtd_data); break;
                   case 8: dipole_restraint(i_c, mtd_data); break;
                   case 11: radgyr_restraint(i_c, mtd_data); break;
                   case 14: rmsdtor_restraint(i_c, mtd_data); break;
                   case 16: dihcor_restraint(i_c, mtd_data); break;
                   case 20: waterbridge_restraint(i_c, mtd_data); break;
                   case 30: spath_restraint(i_c, mtd_data); break;
                   case 31: zpath_restraint(i_c, mtd_data); break;
                   case 32: position_restraint(i_c, mtd_data); break;
                   case 33: elstpot_restraint(i_c, mtd_data); break;
                   case 34: puckering_restraint(i_c, mtd_data); break;
                   case 36: helix_restraint(i_c, mtd_data); break;
                   case 37: alpharmsd_restraint(i_c, mtd_data); break;
                   case 38: antibetarmsd_restraint(i_c, mtd_data); break;
                   case 39: parabetarmsd_restraint(i_c, mtd_data); break;
                   //case 40: camshift_restraint(i_c, mtd_data); break;
                  }
              }
              // Now P and Q using recon plugin
              for(it=0;it<colvar.nconst;++it) welikedupes_s[it]=colvar.ss0[it];
              recon_calc_cv( colvar.nconst, welikedupes_s, &val2, myreconObj); 

              testforce[ix] = 0.5*((val1*invstep)-(val2*invstep));
            
              mtd_data->pos[iat][ix] += teststep;
              for(i_c=0;i_c<colvar.nconst;i_c++) {
                  switch(colvar.type_s[i_c]){
                    case 1: dist_restraint(i_c, mtd_data); break;
                    case 2: mindist_restraint(i_c, mtd_data); break;
                    case 3: coord_restraint(i_c, mtd_data); break;
                    case 4: angle_restraint(i_c, mtd_data); break;
                    case 5: torsion_restraint(i_c, mtd_data); break;
                    case 6: alfabeta_restraint(i_c, mtd_data); break;
                    case 7: hbonds_restraint(i_c, mtd_data); break;
                    case 8: dipole_restraint(i_c, mtd_data); break;
                    case 11: radgyr_restraint(i_c, mtd_data); break;
                    case 14: rmsdtor_restraint(i_c, mtd_data); break;
                    case 16: dihcor_restraint(i_c, mtd_data); break;
                    case 20: waterbridge_restraint(i_c, mtd_data); break;
                    case 30: spath_restraint(i_c, mtd_data); break;
                    case 31: zpath_restraint(i_c, mtd_data); break;
                    case 32: position_restraint(i_c, mtd_data); break;
                    case 33: elstpot_restraint(i_c, mtd_data); break;
                    case 34: puckering_restraint(i_c, mtd_data); break;
                    case 36: helix_restraint(i_c, mtd_data); break;
                    case 37: alpharmsd_restraint(i_c, mtd_data); break;
                    case 38: antibetarmsd_restraint(i_c, mtd_data); break;
                    case 39: parabetarmsd_restraint(i_c, mtd_data); break;
                    //case 40: camshift_restraint(i_c, mtd_data); break;
                  }
              }
              for(it=0;it<colvar.nconst;++it) welikedupes_s[it]=colvar.ss0[it];
              recon_derivatives( colvar.nconst,  welikedupes_s, recon_der , myreconObj );

              analder=0;                  
              for(i_c=0;i_c<colvar.nconst;i_c++) { 
                  for(j=0;j<colvar.natoms[i_c];j++) {
                      if(colvar.cvatoms[i_c][j]==iat) {
                        analder+=recon_der[i_c]*colvar.myder[i_c][j][ix];
                      }
                  }
              }
              printf("Force atom %5i[%i] ** analytic %15.10f ** numeric %15.10f *** DELTA %15.10f\n",
                iat, ix, analder, testforce[ix], analder-testforce[ix]);
        }
     }
  }

}
#endif

void PREFIX debug_derivatives(struct mtd_data_s *mtd_data,int print){
  static int first=1;
  static rvec** pos=NULL;
  static rvec** der=NULL;
  static rvec** previous_pos=NULL;
  static rvec** previous_der=NULL;
  static real* integral;
  int icv,i,j,ncv,iat;
  static FILE* debug_file=NULL;

  ncv = colvar.nconst;
  if(first){
    snew(integral,ncv);
    snew(pos,ncv);
    snew(der,ncv);
    snew(previous_pos,ncv);
    snew(previous_der,ncv);
    for(icv=0;icv<ncv;icv++){
      integral[icv]=colvar.ss0[icv];
      snew(pos[icv],colvar.natoms[icv]);
      snew(der[icv],colvar.natoms[icv]);
      snew(previous_pos[icv],colvar.natoms[icv]);
      snew(previous_der[icv],colvar.natoms[icv]);
    }
  }

  for(icv=0;icv<ncv;icv++){
    for(i=0;i<colvar.natoms[icv];i++){
      iat = colvar.cvatoms[icv][i];
      for(j=0;j<3;j++){
        der[icv][i][j]=colvar.myder[icv][i][j];
        pos[icv][i][j]=mtd_data->pos[iat][j];
      }
    }
  }

  if(!first){
    for(icv=0;icv<ncv;icv++){
      for(i=0;i<colvar.natoms[icv];i++){
        for(j=0;j<3;j++){
          integral[icv]+=
          0.5*(der[icv][i][j]+previous_der[icv][i][j])*(pos[icv][i][j]-previous_pos[icv][i][j]);
        }
      }
    }
  }

  if(print){
    if(!debug_file) debug_file = fopen((mtd_data->ionode?"DEBUG":"/dev/null"), "w");
    fprintf(debug_file, "%10.3f", mtd_data->time);
    for(icv=0;icv<colvar.nconst;icv++) fprintf(debug_file, "   %10.8f", integral[icv]-colvar.ss0[icv]);
    fprintf(debug_file, "\n");
    fflush(debug_file);
  }

  for(icv=0;icv<ncv;icv++){
    for(i=0;i<colvar.natoms[icv];i++){
      for(j=0;j<3;j++){
        previous_der[icv][i][j]=der[icv][i][j];
        previous_pos[icv][i][j]=pos[icv][i][j];
      }
    }
  }

  first=0;
};
