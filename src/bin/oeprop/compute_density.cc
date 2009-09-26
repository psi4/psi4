/*! 
** \file
** \ingroup OEPROP
** \brief Compute the density for Hartree-Fock wavefunctions
*/

#define EXTERN
#include "includes.h"
#include "globals.h"
#include "prototypes.h"

namespace psi { namespace oeprop {

void compute_density()
{
  int i,j,k,l,count,open_1st;
  int nalpha, nbeta;
  double *occ;		/* occupation vector - nmo long */
  int *occ_a,*occ_b,*dummy; /* alpha and beta occupation strings */
  double *alpha;
  double *beta;
  double **ccvecs;	/* vectors of coupling coeffs alpha and beta, packed */
  double occ_tcscf[2];	/* occupations for first and second 
  			   open-shell orbitals in a closed-shell TCSCF case */
  double tmp_a, tmp_b;
  
  
  			/* Obtaining occupation vector */
  occ = init_array(nmo);
  openirrs = 0;
  openmos = 0;
  for (i=0;i<nirreps;i++) {
    if (openpi[i] != 0)
      openirrs++;
    openmos += openpi[i];
  }

	/* Reading in coupling coefficients for TCSCF wavefunction */

  if (iopen < 0) {
    ccvecs = chkpt_rd_ccvecs();
    alpha = ccvecs[0];
    beta = ccvecs[1];
    if (print_lvl >= PRINTCCOEFFLEVEL) {
      fprintf(outfile,"  Coupling coefficient vectors :\n");
      fprintf(outfile,"  i  j\t  alpha \t   beta \n");
      for (i=0;i<openirrs;i++)
        for (j=0;j<=i;j++)
          fprintf(outfile,"%3d%3d\t%8.2f\t%8.2f\n",
                          i+1,j+1,alpha[ioff[i]+j],beta[ioff[i]+j]);
      fprintf(outfile,"\n\n");
    }
  }


		/* Computing occupations */

  if (iopen >= 0) {	/* single-reference case */
    count = 0;
    for (i=0;i<nirreps;i++) {
      for (j=0;j<clsdpi[i];j++)
        occ[count++] = 2.0;
      for (j=0;j<openpi[i];j++)
        occ[count++] = 1.0;
      count += orbspi[i] - openpi[i] - clsdpi[i];
    }
  }
  else {		/* TCSCF for closed shells */
    occ_tcscf[0] = 2/(1-alpha[0]);	/* occupation for the first conf */
    occ_tcscf[1] = 2/(1-alpha[2]);	/* occupation for the second one */
    count = 0;
    k = 0;
    for (i=0;i<nirreps;i++) {
      for (j=0;j<clsdpi[i];j++)
        occ[count++] = 2.0;
      for (j=0;j<openpi[i];j++)
        occ[count++] = occ_tcscf[k++];
      count += orbspi[i] - openpi[i] - clsdpi[i];
    }
  }

  Ptot = init_array(natri);


		/* Computing total density and spin density if required */

  if ((iopen >= 0) && spin_prop) {
    occ_a = init_int_array(nirreps);
    occ_b = init_int_array(nirreps);
    for (i=0;i<nirreps;i++)
      occ_a[i] = occ_b[i] = clsdpi[i];
    if (iopen) {
      for (i=0;i<nirreps;i++)
        if (openpi[i] != 0) {
          occ_a[i] += openpi[i];
          open_1st = i;
          break;
        }
      count = 1;
      for (i=open_1st+1;i<nirreps;i++)
        if (openpi[i] != 0) {
          if ((int)beta[ioff[count]] == 3)
            occ_b[i] += openpi[i];
          else if ((int)beta[ioff[count]] == -1)
	    occ_a[i] += openpi[i];
	  else
            throw PsiException("Can't assign alpha and beta occupation vectors!\nPossibly, you provided erroneous BETA string", __FILE__, __LINE__);
//	    punt("Can't assign alpha and beta occupation vectors!\nPossibly, you provided erroneous BETA string");
          count++;
        }
    }
      
      nalpha = nbeta = 0;
      for(i=0;i<nirreps;i++) {
        nalpha += occ_a[i];
        nbeta += occ_b[i];
      }
      if (nalpha < nbeta) {	/* Interchanging alpha and beta */
        count = nalpha;
        nalpha = nbeta;
        nbeta = count;
        dummy = occ_a;
        occ_a = occ_b;
        occ_b = dummy;
      }
      if (print_lvl >= PRINTOCCUPLEVEL) {
        fprintf(outfile,"  Nalpha = %d, Nbeta = %d\n",nalpha,nbeta);
        fprintf(outfile,"\n  Irrep #   alpha   beta\n");
        for (i=0;i<nirreps;i++)
          fprintf(outfile,"  %5d   %5d   %5d\n",i+1,occ_a[i],occ_b[i]);
        fprintf(outfile,"\n\n");
      }
      
      		/* Computing spin and total densities */
      Pspin = init_array(natri);
      for(i=0;i<nbfao;i++)
        for(j=0;j<=i;j++) {
          count = 0;
          tmp_a = tmp_b = 0.0;
          for(k=0;k<nirreps;k++) {
				/****** Note : For UHF case use corresponding 
				eigenvectors for alpha and beta summation ****/
            for(l=0; l < occ_a[k]; l++)
              tmp_a += scf_evec_ao[i][count+l]*scf_evec_ao[j][count+l];
            for(l=0; l < occ_b[k]; l++)
              tmp_b += scf_evec_ao[i][count+l]*scf_evec_ao[j][count+l];
            count += orbspi[k];
          }
          Ptot[ioff[i]+j] = tmp_a + tmp_b;
          Pspin[ioff[i]+j] = tmp_a - tmp_b;
        }
      if (print_lvl >= PRINTOPDMLEVEL) {
        fprintf(outfile,"\tTotal density matrix in AO basis :\n");
        print_array(Ptot,nbfao,outfile);
        fprintf(outfile,"\n\n");
        fprintf(outfile,"\tSpin density matrix in AO basis :\n");
        print_array(Pspin,nbfao,outfile);
        fprintf(outfile,"\n\n");
      }
  }
  
  else {	/* Computing only total density */
    for(i=0;i<nbfao;i++)
      for(j=0;j<=i;j++) {
        count = 0;
        for(k=0;k<nirreps;k++) {
          for(l=0; l < (clsdpi[k]+openpi[k]); l++) {
            Ptot[ioff[i]+j] += occ[count]*scf_evec_ao[i][count]*scf_evec_ao[j][count];
            count++;
          }
          count += orbspi[k] - openpi[k] - clsdpi[k];
        }
      }
  
    if (print_lvl >= PRINTOPDMLEVEL) {
      fprintf(outfile,"\tTotal density matrix in AO basis :\n");
      print_array(Ptot,nbfao,outfile);
      fprintf(outfile,"\n\n");
    }
  }
  

		/* Cleaning up */

  free(occ);
  if ((iopen >= 0) && spin_prop) {
    free(occ_a);
    free(occ_b);
  }
  free(clsdpi);
  free(openpi);
/*  free(alpha);
  free(beta); */

}

}} // namespace psi::oeprop
