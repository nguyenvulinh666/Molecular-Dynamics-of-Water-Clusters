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

void PREFIX coord_restraint(int i_c, struct mtd_data_s *mtd_data)
{
  int firstAtom, secondAtom, i, j, ix;
  int nn = colvar.nn[i_c], mm = colvar.mm[i_c];
  rvec rij;
  real num, iden, mod_rij, rdist, func, dfunc, rNdist, rMdist;
  real ncoord, r_0 = colvar.r_0[i_c], d_0 = colvar.d_0[i_c];
  real threshold;

  threshold=pow(0.00001,1./(nn-mm));

  ncoord = 0.;                                                                    
  for(i=0;i<colvar.natoms[i_c];i++) for(ix=0;ix<3;ix++) colvar.myder[i_c][i][ix] = 0.;

  for(i=0;i<colvar.list[i_c][0];i++) {                                           // sum over CoordNumber(i)
    firstAtom = colvar.cvatoms[i_c][i];
    for(j=colvar.list[i_c][0];j<colvar.natoms[i_c];j++) {
      if(colvar.logic[i_c] == 1 && (j-colvar.list[i_c][0]) != i) continue;
      secondAtom = colvar.cvatoms[i_c][j];
      if(colvar.cell_pbc[i_c]){
        minimal_image(mtd_data->pos[firstAtom], mtd_data->pos[secondAtom], &mod_rij, rij);
      } else {
        rij[0] = mtd_data->pos[firstAtom][0]-mtd_data->pos[secondAtom][0];
        rij[1] = mtd_data->pos[firstAtom][1]-mtd_data->pos[secondAtom][1];
        rij[2] = mtd_data->pos[firstAtom][2]-mtd_data->pos[secondAtom][2];
        mod_rij  = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
      };
      rdist = (mod_rij-d_0)/r_0;
      /* analitic limit of the switching function */
      if(rdist<=0.){
       ncoord+=1.;
       dfunc=0.;
      }else if(rdist>0.999999 && rdist<1.000001){
       ncoord+=nn/mm;
       dfunc=0.5*nn*(nn-mm)/mm;
      }else if(rdist>threshold){
       dfunc=0.;
      }else{
        rNdist = pow(rdist, nn-1);
        rMdist = pow(rdist, mm-1);
        num = 1.-rNdist*rdist;
        iden = 1./(1.-rMdist*rdist);
        func = num*iden;
        ncoord += func;
        dfunc = ((-nn*rNdist*iden)+(func*(iden*mm)*rMdist))/(mod_rij*colvar.r_0[i_c]);
      }
      for(ix=0;ix<3;ix++) {
        colvar.myder[i_c][i][ix] += +dfunc*rij[ix];
        colvar.myder[i_c][j][ix] += -dfunc*rij[ix];
      }
    }
  }
  colvar.ss0[i_c] = ncoord;

}

//--------------------------------------------------------------------------------------------------------------

int PREFIX read_coord(char **word, int count, t_plumed_input *input, FILE *fplog)
{
  int i, iw, iat, j, help;
  double r_0, d_0;
  double delta = 0.0;
  char string[400];
  real threshold, value;

  help=0;
  d_0=0.;

  colvar.cell_pbc[count]=1; // default is PBC

  iw = seek_word(word,"LIST");
  if(iw>=0){   
             j=plumed_get_group(word[iw+1],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][0]=j;
             j=plumed_get_group(word[iw+2],&colvar.cvatoms[count],colvar.natoms[count],input,fplog);
             colvar.natoms[count]+=j;
             colvar.list[count][1]=j;
  } else{ fprintf(fplog,"|- NEEDED LIST KEYWORD FOR COORD\n"); help=1;}
  iw=seek_word(word,"PAIR");
  if(iw>=0) { colvar.logic[count] = 1;} 
  iw=seek_word(word,"SIGMA");
  if(iw>=0){ sscanf(word[iw+1],"%lf", &delta);
             colvar.delta_r[count]  = (real) delta; }
  iw=seek_word(word,"NN");
  if(iw>=0) { sscanf(word[iw+1],"%i", &colvar.nn[count]); } else { fprintf(fplog,"|- NEEDED NN KEYWORD FOR COORD\n"); help=1;}
  iw=seek_word(word,"MM");
  if(iw>=0) { sscanf(word[iw+1],"%i", &colvar.mm[count]);} else { fprintf(fplog,"|- NEEDED MM KEYWORD FOR COORD\n"); help=1;}
  iw=seek_word(word,"R_0");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &r_0); } else { fprintf(fplog,"|- NEEDED R_0 KEYWORD FOR COORD\n"); help=1;}
  iw=seek_word(word,"D_0");
  if(iw>=0) { sscanf(word[iw+1],"%lf", &d_0); }
  iw=seek_word(word,"NOPBC");
  if(iw>=0) {colvar.cell_pbc[count] = 0;}
  iw=seek_word(word,"PBC");
  if(iw>=0) {colvar.cell_pbc[count] = 1;}

  if(help){
          fprintf(fplog, "\n-COORD CV: WRONG SYNTAX\n");
          fprintf(fplog, "e.g.:     \n");
          fprintf(fplog, "COORD LIST <g1> <g2>  NN 6 MM 12 R_0 0.75 D_0 3. SIGMA 1.0 \n");
          fprintf(fplog, "         g1->    \n");
          fprintf(fplog, "         5 1 6    \n");
          fprintf(fplog, "         g1<-    \n");
          fprintf(fplog, "       \n");
          fprintf(fplog, "         g2->    \n");
          fprintf(fplog, "         LOOP 6921 21786 5 \n");
          fprintf(fplog, "         g2<-    \n");
          fprintf(fplog, "       \n");
          fprintf(stderr, "PluMed dead with errors: check log file  \n");
          EXIT(); 
  }

  colvar.r_0[count]      = (real) r_0;
  colvar.d_0[count]      = (real) d_0;
  colvar.type_s[count]   = 3;

  fprintf(fplog, "%1i-COORDINATION NUMBER OF (1st SET: %i ATOMS) WRT (2nd SET: %i ATOMS); ", count+1, colvar.list[count][0], colvar.list[count][1]);
  if(colvar.cell_pbc[count]) fprintf(fplog, " PBC ON");
  else                       fprintf(fplog, " PBC OFF");
  if (logical.do_hills) fprintf(fplog," SIGMA %f\n",colvar.delta_r[count]);
  else fprintf(fplog,"\n"); 
  fprintf(fplog, "|--FUNCTIONAL FORM: (1-((d_ij-d_0)/r_0)^n) / (1-((d_ij-d_0)/r_0)^m) \n");
  fprintf(fplog, "|--PARAMETERS: n= %i m= %i r_0= %f d_0= %f\n", colvar.nn[count], colvar.mm[count], colvar.r_0[count], colvar.d_0[count]);
  threshold=pow(0.00001,1./(colvar.nn[count]-colvar.mm[count]));
  value=(1.-pow(threshold,colvar.nn[count]))/(1.-pow(threshold,colvar.mm[count]));
  fprintf(fplog, "|--CUTOFF VALUE: %f\n",value);
  fprintf(fplog, "|--CUTOFF DISTANCE: %f\n",threshold*r_0+d_0);
  if(colvar.logic[count]) fprintf(fplog,"|--PAIR COORDINATION ACTIVE \n"); 
  if(colvar.logic[count] == 1 && colvar.list[count][0] != colvar.list[count][1]) {
    fprintf(fplog, "Error!! If PAIR is active <g1> anf <g2> must have the same number of atoms \n");
    EXIT();    
  }

  iat=0;
  fprintf(fplog,"|- 1st SET MEMBERS: ");
  for(i=0;i<colvar.list[count][0];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n");
  fprintf(fplog,"|- 2nd SET MEMBERS: ");
  for(i=0;i<colvar.list[count][1];i++){fprintf(fplog," %d ",colvar.cvatoms[count][iat++]+1);if((i+1)%20==0)fprintf(fplog,"\n                    ");}fprintf(fplog,"\n\n");

  snew(colvar.myder[count], colvar.natoms[count]);

  return colvar.natoms[count]; 
}
