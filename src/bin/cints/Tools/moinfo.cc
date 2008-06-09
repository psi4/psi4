/*! \file moinfo.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<cmath>
#include<libciomr/libciomr.h>
#include<libchkpt/chkpt.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"small_fns.h"

namespace psi { namespace CINTS {

/*-------------------------------
  Explicit function declarations
 -------------------------------*/
static void print_ccoefs(void);

void init_moinfo()
{
  int i, j;
  int irrep, iopen;
  int openirrs;
  double **ccvecs, *alpha, *beta;
  double **scf_evec_so;

  MOInfo.Escf = chkpt_rd_escf();
  MOInfo.Eref = 0.0;
  MOInfo.Ecorr = 0.0;

  /* CDS: I revised this stuff about correlation and SCF energies */

  /* If SCF, can say Eref = Escf (I guess...) */
  if (strcmp(UserOptions.wfn,"SCF")==0) {
    MOInfo.Eref = MOInfo.Escf;
  }
  else {
    MOInfo.Eref = chkpt_rd_eref();
  }

  /* Note: this init_moinfo() routine is not always called!  We'll need
     to re-grab some of the above energies in other subroutines to be
     positive we have them on-hand */

  MOInfo.num_mo = chkpt_rd_nmo();
  MOInfo.orbspi = chkpt_rd_orbspi();
  MOInfo.clsdpi = chkpt_rd_clsdpi();
  MOInfo.openpi = chkpt_rd_openpi();
  MOInfo.virtpi = init_int_array(Symmetry.nirreps);
  for(irrep=0;irrep<Symmetry.nirreps;irrep++)
    MOInfo.virtpi[irrep] = MOInfo.orbspi[irrep] - MOInfo.clsdpi[irrep] - MOInfo.openpi[irrep];
  MOInfo.scf_evals[0] = chkpt_rd_evals();
  scf_evec_so = chkpt_rd_scf();
  MOInfo.scf_evec[0] = block_matrix(MOInfo.num_mo,BasisSet.num_ao);
  mmult(Symmetry.usotao,1,scf_evec_so,0,MOInfo.scf_evec[0],1,BasisSet.num_ao,Symmetry.num_so,MOInfo.num_mo,0);
  free_block(scf_evec_so);
  if (UserOptions.reftype == uhf) {
    MOInfo.scf_evals[1] = chkpt_rd_beta_evals();
    scf_evec_so = chkpt_rd_beta_scf();
    MOInfo.scf_evec[1] = block_matrix(MOInfo.num_mo,BasisSet.num_ao);
    mmult(Symmetry.usotao,1,scf_evec_so,0,MOInfo.scf_evec[1],1,BasisSet.num_ao,Symmetry.num_so,MOInfo.num_mo,0);
    free_block(scf_evec_so);
  }

  /*--- Check the validity of the checkpoint file ---*/
  iopen = chkpt_rd_iopen();
  switch (UserOptions.reftype) {
  case rhf:     if (iopen != 0) throw std::domain_error("Content of checkpoint file inconsistent with REFERENCE\n"); break;
  case uhf:     if (iopen != 0) throw std::domain_error("Content of checkpoint file inconsistent with REFERENCE\n"); break;
  case rohf:    if (iopen <= 0) throw std::domain_error("Content of checkpoint file inconsistent with REFERENCE\n"); break;
  case twocon:  if (iopen >= 0) throw std::domain_error("Content of checkpoint file inconsistent with REFERENCE\n"); break;
  }

  /*--- Number of d.-o. MOs and s.-o. MOs ---*/
  MOInfo.ndocc = 0;
  openirrs = 0;
  MOInfo.nsocc = 0;
  for (i=0;i<Symmetry.nirreps;i++) {
    MOInfo.ndocc += MOInfo.clsdpi[i];
    if (MOInfo.openpi[i] != 0)
      openirrs++;
    MOInfo.nsocc += MOInfo.openpi[i];
  }
  MOInfo.nuocc = MOInfo.num_mo - MOInfo.ndocc - MOInfo.nsocc;
  
  /*--- Number of closed and open shells (no virtuals) ---*/
  MOInfo.num_moshells = MOInfo.num_openmoshells = ( -1 + (int)sqrt(1.0+8.0*abs(iopen)) )/2;    /* number of open shells */
  if (MOInfo.ndocc > 0)  /* add closed shells as well */
    MOInfo.num_moshells++;

  
  /*--- Read in open-shell coupling coeffcients ---*/
  ccvecs = chkpt_rd_ccvecs();
  if (iopen != 0) {    /*--- NOTE! These are Pitzer's coupling constants (a and b).
			 To get Yamaguchi's constants (alpha and beta) use this:
			 alpha = (1-a)/2  beta = (b-1)/4
			---*/
    alpha = ccvecs[0];
    beta = ccvecs[1];
  }

  if (UserOptions.reftype == twocon) {
    MOInfo.tcscf_occ[0] = 2.0/(1.0-alpha[0]);
    MOInfo.tcscf_occ[1] = 2.0/(1.0-alpha[2]);
  }

  if (UserOptions.reftype == rohf || UserOptions.reftype == twocon) {
    /*--- Form square matrices of coupling coeffcients ---*/
    MOInfo.Alpha = block_matrix(MOInfo.num_moshells,MOInfo.num_moshells);
    MOInfo.Beta  = block_matrix(MOInfo.num_moshells,MOInfo.num_moshells);
    /* Put alpha's and beta's in Alpha and Beta */
    if (MOInfo.ndocc > 0) { /* There are closed shells */
      /* Closed-Closed CCs */
      MOInfo.Alpha[0][0] = 2.0;
      MOInfo.Beta[0][0] = -1.0;
      if (UserOptions.reftype == rohf) {           /*--- Highspin and opeh-shell singlet cases ---*/
	/* Closed-Open CCs */
	for(i=1;i<MOInfo.num_moshells;i++) {
	  MOInfo.Alpha[0][i] = MOInfo.Alpha[i][0] = 1.0;
	  MOInfo.Beta[0][i] = MOInfo.Beta[i][0] = -0.5;
	}
	/* Open-Open blocks */
	for(i=0;i<MOInfo.num_openmoshells;i++)
	  for(j=0;j<=i;j++) {
	    MOInfo.Alpha[i+1][j+1] = MOInfo.Alpha[j+1][i+1] = (1.0 - alpha[ioff[i]+j]) * 0.5;
	    MOInfo.Beta[i+1][j+1] = MOInfo.Beta[j+1][i+1] = (beta[ioff[i]+j] + 1.0) * -0.25;
	  }
      }
      else if (UserOptions.reftype == twocon) {      /*--- TCSCF ---*/
	MOInfo.Alpha[0][1] = MOInfo.Alpha[1][0] = MOInfo.tcscf_occ[0];
	MOInfo.Beta[0][1] = MOInfo.Beta[1][0] = (-0.5) * MOInfo.tcscf_occ[0];
	MOInfo.Alpha[0][2] = MOInfo.Alpha[2][0] = MOInfo.tcscf_occ[1];
	MOInfo.Beta[0][2] = MOInfo.Beta[2][0] = (-0.5) * MOInfo.tcscf_occ[1];
	MOInfo.Alpha[1][1] = MOInfo.tcscf_occ[0] * 0.5;
	MOInfo.Alpha[2][2] = MOInfo.tcscf_occ[1] * 0.5;
	MOInfo.Beta[1][2] = MOInfo.Beta[2][1] = sqrt(MOInfo.tcscf_occ[0]*MOInfo.tcscf_occ[1]) * (-0.5);
      }
    }
    else { /* There are no closed shells */
      if (UserOptions.reftype == rohf) {           /*--- Highspin and opeh-shell singlet cases ---*/
	/* Open-Open blocks */
	for(i=0;i<MOInfo.num_openmoshells;i++)
	  for(j=0;j<=i;j++) {
	    MOInfo.Alpha[i][j] = MOInfo.Alpha[j][i] = (1.0 - alpha[ioff[i]+j]) * 0.5;
	    MOInfo.Beta[i][j] = MOInfo.Beta[j][i] = (beta[ioff[i]+j] - 1.0) * 0.25;
	  }
      }
      else if (UserOptions.reftype == twocon) {      /*--- TCSCF ---*/
	MOInfo.Alpha[0][0] = MOInfo.tcscf_occ[0] * 0.5;
	MOInfo.Alpha[1][1] = MOInfo.tcscf_occ[1] * 0.5;
	MOInfo.Beta[0][1] = MOInfo.Beta[1][0] = sqrt(MOInfo.tcscf_occ[0]*MOInfo.tcscf_occ[1]) * (-0.5);
      }
    }
    print_ccoefs();
  }

  free_block(ccvecs);

  return;
}


void cleanup_moinfo()
{
  free(MOInfo.scf_evals[0]);
  free_block(MOInfo.scf_evec[0]);
  if (UserOptions.reftype == uhf) {
    free(MOInfo.scf_evals[1]);
    free_block(MOInfo.scf_evec[1]);
  }
  free(MOInfo.orbspi);
  free(MOInfo.clsdpi);
  free(MOInfo.openpi);
  free(MOInfo.virtpi);
  if (UserOptions.reftype == rohf || UserOptions.reftype == twocon) {
    free_block(MOInfo.Alpha);
    free_block(MOInfo.Beta);
  }

  return;
}


/*!----------------------
  Print coupling coeffs
 ----------------------*/
void print_ccoefs()
{
  int i,j;

  if (UserOptions.print_lvl >= PRINT_CCOEFF) {
    fprintf(outfile,"  -Yamaguchi's coupling coefficients :\n\n");
    fprintf(outfile,"    i  j    alpha     beta \n");
    fprintf(outfile,"   -------------------------\n");
    for (i=0;i<MOInfo.num_moshells;i++)
      for (j=0;j<=i;j++)
	fprintf(outfile,"   %2d %2d    %5.3f    %5.3f\n",
		i+1,j+1,MOInfo.Alpha[i][j],MOInfo.Beta[i][j]);
    fprintf(outfile,"\n\n");
  }

  return;
}
};};
