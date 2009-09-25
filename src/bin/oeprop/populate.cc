/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "prototypes.h"
#include "globals.h"

namespace psi { namespace oeprop {

void populate()
{
  int i,j,k,l,count,jmax;
  int ao_ind,ao_i_1st,ao_i_last,ao_j_1st,ao_j_last;
  int atom_i,atom_j;
  double **aopop;		/* AO population matrix */
  double **aospop;              /* AO spin population matrix */
  double *gop;			/* Gross orbital population */
  double *gsop;                 /* Gross spin orbital population */
  double **batm, *qatm;		/* Mulliken's atomic bond and gross atomic population */
  double tmp;
  
  aopop = init_matrix(nbfao,nbfao);
  for(i=0;i<nbfao;i++)
    for(j=0;j<=i;j++)
      aopop[j][i] = (aopop[i][j] = Ptot[ioff[i]+j]*S[ioff[i]+j]);

  if (spin_prop) {
    aospop = init_matrix(nbfao,nbfao);
    for(i=0;i<nbfao;i++)
      for(j=0;j<=i;j++)
        aospop[j][i] = (aospop[i][j] = Pspin[ioff[i]+j]*S[ioff[i]+j]);
  }

  
  if (print_lvl > PRINTAOPOPLEVEL) {
    fprintf(outfile,"  AO population matrix :\n");
    print_mat(aopop,nbfao,nbfao,outfile);
    fprintf(outfile,"\n\n");
  }

  gop = init_array(nbfao);
  for(i=0;i<nbfao;i++) {
    tmp = 0.0;
    for(j=0;j<nbfao;j++)
      tmp += aopop[i][j];
    gop[i] = tmp;
  }

  if (spin_prop) {
    gsop = init_array(nbfao);
    for(i=0;i<nbfao;i++) {
      tmp = 0.0;
      for(j=0;j<nbfao;j++)
        tmp += aospop[i][j];
      gsop[i] = tmp;
    }
  }
  
  if (!spin_prop) {
    fprintf(outfile," -Gross orbital populations :\n\n");
    fprintf(outfile,"  Center   AO          L          Q(AO)\n");
    fprintf(outfile,"  ------  ----        ---      -----------\n");
  }
  else {
    fprintf(outfile," -Gross orbital and spin-orbital populations :\n\n");

    fprintf(outfile,"  Center   AO          L          Q(AO)       AlphQ(AO)     BetaQ(AO)\n");
    fprintf(outfile,"  ------  ----        ---      -----------   -----------   -----------\n");
  }
  count = 0;
  for(i=0;i<nshell;i++) {
    jmax = stype[i]*(stype[i]+1)/2;
    for(j=0;j<stype[i]*(stype[i]+1)/2;j++) {
      ao_ind = sloc[i]+j;
      if (!spin_prop)
        fprintf(outfile,"%6d  %5d         %2d       %11.8lf\n",
		snuc[i],ao_ind,stype[i]-1,gop[ao_ind-1]);
      else
	fprintf(outfile,"%6d  %5d         %2d       %11.8lf   %11.8lf   %11.8lf\n",
		snuc[i],ao_ind,stype[i]-1,gop[ao_ind-1],
		(gop[ao_ind-1]+gsop[ao_ind-1])/2,
		(gop[ao_ind-1]-gsop[ao_ind-1])/2);
    }
  }
  fprintf(outfile,"\n\n");
  
  qatm = init_array(natom);
  qnet = init_array(natom);
  batm = init_matrix(natom,natom);
  for(i=0;i<nshell;i++) {
    atom_i = snuc[i]-1;
    ao_i_1st = sloc[i]-1;
    ao_i_last = ao_i_1st + stype[i]*(stype[i]+1)/2 - 1;
    for(k=ao_i_1st;k<=ao_i_last;k++)
      qatm[atom_i] += gop[k];
    for(j=0;j<=i;j++) {
      atom_j = snuc[j]-1;
      ao_j_1st = sloc[j]-1;
      ao_j_last = ao_j_1st + stype[j]*(stype[j]+1)/2 - 1;
      for(k=ao_i_1st;k<=ao_i_last;k++)
        for(l=ao_j_1st;l<=ao_j_last;l++)
          batm[atom_i][atom_j] += aopop[k][l];
      batm[atom_j][atom_i] = batm[atom_i][atom_j];
    }
    qnet[atom_i] = zvals[atom_i] - qatm[atom_i];
  }

  fprintf(outfile," -Atomic bond populations :\n");
  print_mat(batm,natom,natom,outfile);
  fprintf(outfile,"\n\n");

  fprintf(outfile," -Gross atomic populations and net charges :\n\n");
  fprintf(outfile,"    Center    Atomic Population    Net Charge\n");
  fprintf(outfile,"    ------    -----------------    ----------\n");
  for(i=0;i<natom;i++)
    fprintf(outfile,"%8d          %10.6lf       %+10.6lf\n",i+1,qatm[i],qnet[i]);
  fprintf(outfile,"\n\n");


		/* Cleaning up */
  
  free_matrix(aopop,nbfao);
  free(gop);
  free(qatm);
  free_matrix(batm,natom);
  if (spin_prop) {
    free_matrix(aospop,nbfao);
    free(gsop);
  }
}

}} // namespace psi::oeprop
