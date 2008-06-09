/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <psifiles.h>

#include "input.h"
#include "defines.h"
#define EXTERN
#include "global.h"

namespace {
  void get_oeints();
  double **canon_orthog(double **);
  double **core_vector(double **);
  double **projector(double **, double **);
  double **project(double **, double **, double **);
  void finish_projection();
  
  double **S11; /* Overlap matrix 11 - for the new basis */
  double **H11; /* Core Hamiltonian matrix */
}

namespace psi { namespace input {

/*-------------------------------------------------------
  This function projects the old MOs  from the old basis
  2 onto the new basis 1 :
  
  1) compute overlap S11 between for the new basis
     and X1 - the transformation matrix to the orthogonal
     basis 1' (using canonical orthogonalization with
     elimination of linear dependencies)
  2) compute overlap S21 between the old and new bases
  3) project the old vector C1 = X1 X1(transp) S21 C2
  4) add dummy virtual MOs - those do not
     contribute to the first density anyway
 -------------------------------------------------------*/
void oldcalc_projection()
{
  int i, so, mo;
  double **S12;                 /* Overlap matrix 12 - between the bases */
  double **X1;                  /* transforms from the new basis to the orthogonal new basis
				   (canonical orthogonalization, see Szabo and Ostlund) */
  double **P12;                 /* Projection operator from the old to the new basis */
  double **CV1;                 /* Core Hamiltonian eigenvector -- to be used to provide
				   virtual MOs for the projected vector */
  double **local;               /* localized SCF MOs */

  if (max_angmom > Oldcalc.max_angmom)
      init_gto(max_angmom);
  else
      init_gto(Oldcalc.max_angmom);

  /*--- Call CINTS and get one-electron integrals for new basis ---*/
  get_oeints();
  X1 = canon_orthog(S11);
  free_block(S11);

  CV1 = core_vector(X1);
  free_block(H11);
  
  S12 = overlap_new_old();

  P12 = projector(S12,X1);

  if (dont_project_mos) {
    /* Check that the new calculation can use the old eigenvector */
    for(i=0;i<nirreps;i++)
      if (num_so_per_irrep[i] != Oldcalc.sopi[i] ||
	  orbspi[i] != Oldcalc.orbspi[i])
	punt("Symmetry block structure of old MOs is not suitable for the new calculation. Cannot use old MOs.");

    if (Oldcalc.spinrestr_ref) {
      scf_evect_so = block_matrix(num_so,num_mo);
      for(so=0;so<num_so;so++)
	for(mo=0;mo<num_mo;mo++)
	  scf_evect_so[so][mo] = Oldcalc.scf_evect_so[so][mo];
    }
    else {
      scf_evect_so_alpha = block_matrix(num_so,num_mo);
      for(so=0;so<num_so;so++)
	for(mo=0;mo<num_mo;mo++)
	  scf_evect_so_alpha[so][mo] = Oldcalc.scf_evect_so_alpha[so][mo];
      scf_evect_so_beta = block_matrix(num_so,num_mo);
      for(so=0;so<num_so;so++)
	for(mo=0;mo<num_mo;mo++)
	  scf_evect_so_beta[so][mo] = Oldcalc.scf_evect_so_beta[so][mo];
    }

    /* copy localized MOs, too, if available */
    if(Oldcalc.local != NULL) {
      scf_evect_local = block_matrix(num_so,num_mo);
      for(so=0;so<num_so;so++)
	for(mo=0;mo<num_mo;mo++)
	  scf_evect_local[so][mo] = Oldcalc.local[so][mo];
    }
    else scf_evect_local = NULL;
  }
  else {
    if (Oldcalc.spinrestr_ref)
      scf_evect_so = project(P12,Oldcalc.scf_evect_so,CV1);
    else {
      scf_evect_so_alpha = project(P12,Oldcalc.scf_evect_so_alpha,CV1);
      scf_evect_so_beta = project(P12,Oldcalc.scf_evect_so_beta,CV1);
    }
  }
  
  free_block(S12);
  free_block(P12);
  free_block(X1);

  /*--- misc. stuff like copying escf, ref, compute mxcoef etc. ---*/
  finish_projection();

  if (print_lvl >= DEBUGPRINT) {
      if (spinrestr_ref) {
	  fprintf(outfile,"  -Projected eigenvector (in SO basis):\n");
	  print_mat(scf_evect_so,num_so,num_mo,outfile);
      }
      else {
	  fprintf(outfile,"  -Projected alpha eigenvector (in SO basis):\n");
	  print_mat(scf_evect_so_alpha,num_so,num_mo,outfile);
	  fprintf(outfile,"  -Projected beta eigenvector (in SO basis):\n");
	  print_mat(scf_evect_so_beta,num_so,num_mo,outfile);
      }
  }
      
  return;
}

}} // namespace psi::input

namespace {

using namespace psi;
using namespace psi::input;

/*--------------------------------------------------------
  This function calls CINTS to get one-electron integrals
 --------------------------------------------------------*/
void get_oeints()
{
  int stat;
  int ntri = num_so*(num_so+1)/2;
  int i, j, count;
  double *ints;
  double e_fzc;

  /* clear out the old PSIF_OEI */
  psio_open(PSIF_OEI, PSIO_OPEN_OLD);
  psio_close(PSIF_OEI, 0);

  stat = system("cints --oeints");
  switch (stat) {
  case 0:
      /* CINTS ran successfully - continue */
      break;

  default:
      /* Something went wrong */
      punt("System call to CINTS failed. Check to see if it's in your PATH");
  }

  ints = init_array(ntri);

  /* S integrals */
  stat = iwl_rdone(PSIF_OEI,PSIF_SO_S,ints,ntri,0,0,outfile);
  S11 = block_matrix(num_so,num_so);
  count = 0;
  for(i=0;i<num_so;i++)
      for(j=0;j<=i;j++)
	  S11[i][j] = S11[j][i] = ints[count++];

  
  /* T integrals */
  stat = iwl_rdone(PSIF_OEI,PSIF_SO_T,ints,ntri,0,0,outfile);
  H11 = block_matrix(num_so,num_so);
  count = 0;
  for(i=0;i<num_so;i++)
      for(j=0;j<=i;j++)
	  H11[i][j] = H11[j][i] = ints[count++];

  /* V integrals */
  stat = iwl_rdone(PSIF_OEI,PSIF_SO_V,ints,ntri,1,0,outfile);
  count = 0;
  for(i=0;i<num_so;i++)
      for(j=0;j<=i;j++) {
	  H11[i][j] += ints[count++];
	  H11[j][i] = H11[i][j];
      }

  if (print_lvl >= DEBUGPRINT) {
      fprintf(outfile,"  -Overlap matrix in the new basis (in SO basis):\n");
      print_mat(S11,num_so,num_so,outfile);
      fprintf(outfile,"  -Core Hamiltonian matrix in the new basis (in SO basis):\n");
      print_mat(H11,num_so,num_so,outfile);
  }
  
  free(ints);
}


double **canon_orthog(double **S11)
{
  int so, so1, so2, oo, mo;
  int irrep, so_offset, oo_offset, blksz;
  int mostart;
  double **X;
  double **symblk, **transmat;
  double *evals, **evecs;
  double feval, min_eval, sahalf;

  /*---------------------------------------------
    Diagonalize sym. blocks of S and find the
    number of lin.-independent functions in each
   ---------------------------------------------*/
  symblk = block_matrix(num_so,num_so);
  transmat = block_matrix(num_so,num_so);
  evals = init_array(num_so);
  evecs = block_matrix(num_so,num_so);
  orbspi = init_int_array(nirreps);
  min_eval = 100000.0;
  so_offset = 0;
  oo_offset = 0;
  for(irrep=0;irrep<nirreps;irrep++) {
      blksz = num_so_per_irrep[irrep];
      if (blksz) {
	  for(so1=0;so1<blksz;so1++)
	      for(so2=0;so2<blksz;so2++)
		  symblk[so1][so2] = S11[so_offset+so1][so_offset+so2];

	  sq_rsp(blksz,blksz,symblk,evals,1,evecs,1.0E-14);
	  if (min_eval > evals[0]) min_eval = evals[0];
	  /* count the number of linearly-independent orthogonal orbitals */
	  for(oo=0;oo<blksz;oo++) {
	      feval = fabs(evals[oo]);
	      if (feval > LINDEP_CUTOFF)
		  orbspi[irrep]++;
	  }
	  
	  /* construct the transformation matrix to the orthogonal basis */
	  mostart = blksz-orbspi[irrep];
	  for(oo=0,mo=mostart;mo<blksz;oo++,mo++) {
	      sahalf = 1.0/sqrt(evals[mo]);
	      for(so=0;so<blksz;so++) {
		  transmat[so+so_offset][oo+oo_offset] = evecs[so][mo]*sahalf;
	      }
	  }
	  so_offset += blksz;
	  oo_offset += orbspi[irrep];
	  
      }
  }
  free_block(evecs);
  free(evals);
  free_block(symblk);

  num_mo = 0;
  for(irrep=0;irrep<nirreps;irrep++)
      num_mo += orbspi[irrep];
  
  X = block_matrix(num_so,num_mo);
  for(so=0;so<num_so;so++)
      for(oo=0;oo<num_mo;oo++)
	  X[so][oo] = transmat[so][oo];
  free_block(transmat);

  return X;
}


double **core_vector(double **X1)
{
  int mo1, mo2, mo_offset, blksz, irrep;
  double **temp1, **temp2, *evals, **evect, **evect_so;

  temp1 = block_matrix(num_mo,num_so);
  temp2 = block_matrix(num_mo,num_mo);
  
  mmult(X1,1,H11,0,temp1,0,num_mo,num_so,num_so,0);
  mmult(temp1,0,X1,0,temp2,0,num_mo,num_so,num_mo,0);

  evals = init_array(num_mo);
  evect = block_matrix(num_mo,num_mo);
  mo_offset = 0;
  for(irrep=0;irrep<nirreps;irrep++) {
      blksz = orbspi[irrep];
      if (blksz) {
	  for(mo1=0;mo1<blksz;mo1++)
	      for(mo2=0;mo2<blksz;mo2++)
		  temp1[mo1][mo2] = temp2[mo_offset+mo1][mo_offset+mo2];
	  sq_rsp(blksz,blksz,temp1,evals,1,temp2,1.0E-14);
	  for(mo1=0;mo1<blksz;mo1++)
	      for(mo2=0;mo2<blksz;mo2++)
		  evect[mo_offset+mo1][mo_offset+mo2] = temp2[mo1][mo2];
      }
      mo_offset += blksz;
  }
  free_block(temp1);
  free_block(temp2);
  free(evals);

  evect_so = block_matrix(num_so,num_mo);
  mmult(X1,0,evect,0,evect_so,0,num_so,num_mo,num_mo,0);
  free_block(evect);

  return evect_so;
}

double **projector(double **S12, double **X1)
{
  double **tmpmat, **P12;

  tmpmat = block_matrix(num_mo,Oldcalc.num_so);
  mmult(X1,1,S12,0,tmpmat,0,
	num_mo,num_so,Oldcalc.num_so,0);
  P12 = block_matrix(num_so,Oldcalc.num_so);
  mmult(X1,0,tmpmat,0,P12,0,
	num_so,num_mo,Oldcalc.num_so,0);
  free_block(tmpmat);

  return P12;
}


/*-----------------------------------------------------------
  This function projects the old eigenvector onto new basis.
  It uses virtual orbitals from core Hamiltonian eigenvector
  (to avoid problems in CSCF)
 -----------------------------------------------------------*/

double **project(double **P12, double **v2, double **core_evect)
{
  int irrep, mo, mo1, mo2, so;
  int nocc, nso, nmo, soffset, moffset1, moffset2;
  double **evect, **v1;

  evect = block_matrix(num_so,Oldcalc.num_mo);
  mmult(P12,0,v2,0,evect,0,
	num_so,Oldcalc.num_so,Oldcalc.num_mo,0);
 
  v1 = block_matrix(num_so,num_mo);
  soffset = moffset1 = moffset2 = 0;
  for(irrep=0;irrep<nirreps;irrep++) {
      nocc = Oldcalc.clsdpi[irrep] + Oldcalc.openpi[irrep];
      nso = num_so_per_irrep[irrep];
      nmo = orbspi[irrep];
      for(mo=0,mo1=moffset1,mo2=moffset2;mo<nmo;mo++,mo1++,mo2++) {
	  if (mo < nocc)
	      for(so=soffset;so<nso+soffset;so++)
		  v1[so][mo1] = evect[so][mo2];
	  else
	      for(so=soffset;so<nso+soffset;so++)
		  v1[so][mo1] = core_evect[so][mo1];
      }
      moffset1 += orbspi[irrep];
      moffset2 += Oldcalc.orbspi[irrep];
      soffset += num_so_per_irrep[irrep];
  }
  free_block(evect);

  return v1;
}


void finish_projection()
{
  int irrep;
  
  ref = Oldcalc.ref;
  spinrestr_ref = Oldcalc.spinrestr_ref;
  escf = Oldcalc.escf;
  iopen = Oldcalc.iopen;
  clsdpi = init_int_array(nirreps);
  openpi = init_int_array(nirreps);
  mxcoef = 0;
  num_so_typs = 0;

  for(irrep=0;irrep<nirreps;irrep++) {
      mxcoef += num_so_per_irrep[irrep]*orbspi[irrep];
      clsdpi[irrep] = Oldcalc.clsdpi[irrep];
      openpi[irrep] = Oldcalc.openpi[irrep];
      if (orbspi[irrep]) num_so_typs++;
  }

  if (print_lvl > 0 && !dont_project_mos) {
      fprintf(outfile,"  -MO projection\n");
      fprintf(outfile,"    MO projection is complete.\n\n");
  }

  return;
}

} // namespace
