/*! \file deriv2.cc
    \ingroup CINTS
    \brief Driver for second derivative integrals.
*/
#include<cstdio>
#include<cstdlib>

#include<libipv1/ip_lib.h>
#include<cmath>
#include<libciomr/libciomr.h>
#include<libpsio/psio.h>
#include<libint/libint.h>
#include<libiwl/iwl.h>
#include<libqt/qt.h>
#include<psifiles.h>
#include <Tools/prints.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"moinfo.h"
#include"compute_scf_opdm.h"
#include"read_gen_opdm.h"
#include"enuc_deriv2.h"
#include"oe_deriv2.h"
#include"te_deriv2_scf.h"
#include"file11.h"

#define DO_NUC 1
#define DO_OE 1
#define DO_TE 1
namespace psi {
  namespace cints {
  //! Driver for integral second derivatives.
  void deriv2(void)
  {
  int i, j, ij, coord, ntri, natom, size;
  char *label;
  double **tmpmat, **my_grad, *outbuf;
  double **C;

  /*--- Hessian in the canonical frame ---*/
  Hess = block_matrix(Molecule.num_atoms*3,Molecule.num_atoms*3);

  /* Derivative Fock and overlap matrices */
  F = (double ***) malloc(Molecule.num_atoms * 3 * sizeof(double **));
  S = (double ***) malloc(Molecule.num_atoms * 3 * sizeof(double **));
  HDS = (double ***) malloc(Molecule.num_atoms * 3 * sizeof(double **));
  for(i=0; i < Molecule.num_atoms*3; i++) {
    F[i] = block_matrix(BasisSet.num_ao, BasisSet.num_ao);
    S[i] = block_matrix(BasisSet.num_ao, BasisSet.num_ao);
    HDS[i] = block_matrix(BasisSet.num_ao, BasisSet.num_ao);
  }

  init_moinfo();
  compute_scf_opdm();

#if DO_NUC
  enuc_deriv2();
#endif
#if DO_OE
  oe_deriv2();
#endif



  /* Multiply the F's by 2 -- accounts for orbital population factor */
  for(coord=0; coord < Molecule.num_atoms*3; coord++) {
    for(i=0; i < BasisSet.num_ao; i++)
      for(j=0; j < BasisSet.num_ao; j++)
	F[coord][i][j] *= 2.0;
  }

#if DO_TE
  te_deriv2_scf();
  /* te_deriv2_scf_symm(); */
#endif

  /* divide the whole thing by 2 */
  for(coord=0; coord < Molecule.num_atoms*3; coord++) {
    for(i=0; i < BasisSet.num_ao; i++)
      for(j=0; j < BasisSet.num_ao; j++)
	F[coord][i][j] *= 0.5;
  }

  /* symmetrize the F's and S's */
  for(coord=0; coord < Molecule.num_atoms*3; coord++) {
    if (UserOptions.print_lvl >= PRINT_OEDERIV) {
      fprintf(outfile, "AO-basis Overlap Derivs (Pre-Symm) (%d)", coord);
      print_mat(S[coord],BasisSet.num_ao,BasisSet.num_ao,outfile);
    }

    for(i=0; i < BasisSet.num_ao; i++) {
      for(j=0; j <= i; j++) {
	if(i!=j) {
	  F[coord][i][j] = F[coord][j][i] = 0.5 * (F[coord][i][j] + F[coord][j][i]);
	  S[coord][i][j] = S[coord][j][i] = 0.5 * (S[coord][i][j] + S[coord][j][i]);
	}
      }
    }

    if (UserOptions.print_lvl >= PRINT_OEDERIV) {
      fprintf(outfile, "AO-basis Overlap Derivs (%d)", coord);
      print_mat(S[coord],BasisSet.num_ao,BasisSet.num_ao,outfile);
    }
  }

  /* Transform Fock and Overlap derivatives to the MO basis */
  tmpmat = block_matrix(MOInfo.num_mo, BasisSet.num_ao);
  for(i=0; i < Molecule.num_atoms*3; i++) {
    C_DGEMM('n','n',MOInfo.num_mo,BasisSet.num_ao,BasisSet.num_ao,1.0,
	    &(MOInfo.scf_evec[0][0][0]),BasisSet.num_ao,&(F[i][0][0]),BasisSet.num_ao,
	    0.0,&(tmpmat[0][0]),BasisSet.num_ao);
    C_DGEMM('n','t',MOInfo.num_mo,MOInfo.num_mo,BasisSet.num_ao,1.0,
	    &(tmpmat[0][0]),BasisSet.num_ao,&(MOInfo.scf_evec[0][0][0]),BasisSet.num_ao,
	    0.0,&(F[i][0][0]),BasisSet.num_ao);

    C_DGEMM('n','n',MOInfo.num_mo,BasisSet.num_ao,BasisSet.num_ao,1.0,
	    &(MOInfo.scf_evec[0][0][0]),BasisSet.num_ao,&(S[i][0][0]),BasisSet.num_ao,
	    0.0,&(tmpmat[0][0]),BasisSet.num_ao);
    C_DGEMM('n','t',MOInfo.num_mo,MOInfo.num_mo,BasisSet.num_ao,1.0,
	    &(tmpmat[0][0]),BasisSet.num_ao,&(MOInfo.scf_evec[0][0][0]),BasisSet.num_ao,
	    0.0,&(S[i][0][0]),BasisSet.num_ao);
  }
  free_block(tmpmat);

  /* write derivative Fock and overlap derivatives to disk */
  ntri = MOInfo.num_mo * (MOInfo.num_mo + 1)/2;
  outbuf = init_array(ntri);
  label = (char *) malloc(PSIO_KEYLEN * sizeof(char));
  for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';
  for(coord=0; coord < Molecule.num_atoms*3; coord++) {

    for(i=0, ij=0; i < MOInfo.num_mo; i++) {
      for(j=0; j <= i; j++, ij++) { /* Lower triangles only */
	outbuf[ij] = F[coord][i][j];
      }
    }

    sprintf(label, "MO-basis Fock Derivs (%d)", coord);
    iwl_wrtone(PSIF_OEI, label, ntri, outbuf);
    if (UserOptions.print_lvl >= PRINT_OEDERIV) {
      fprintf(outfile,"  -%s\n",label);
      print_mat(F[coord],MOInfo.num_mo,MOInfo.num_mo,outfile);
    }
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';
  }
  for(coord=0; coord < Molecule.num_atoms*3; coord++) {

    for(i=0, ij=0; i < MOInfo.num_mo; i++) {
      for(j=0; j <= i; j++, ij++) { /* Lower triangles only */
	outbuf[ij] = S[coord][i][j];
      }
    }

    sprintf(label, "MO-basis Overlap Derivs (%d)", coord);
    iwl_wrtone(PSIF_OEI, label, ntri, outbuf);
    if (UserOptions.print_lvl >= PRINT_OEDERIV) {
      fprintf(outfile,"  -%s\n",label);
      print_mat(S[coord],MOInfo.num_mo,MOInfo.num_mo,outfile);
    }
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';
  }
  /* Write half-differentiated overlap integrals to disk */
  size = BasisSet.num_ao * BasisSet.num_ao;
  for(coord=0; coord < Molecule.num_atoms*3; coord++) {

    sprintf(label, "AO-basis Half-Diff Overlap (%d)", coord);
    psio_open(PSIF_OEI, PSIO_OPEN_OLD);
    psio_write_entry(PSIF_OEI, label, (char *) HDS[coord][0], size*sizeof(double));
    psio_close(PSIF_OEI, 1);
    if (UserOptions.print_lvl >= PRINT_OEDERIV) {
      fprintf(outfile,"  -%s\n",label);
      print_mat(HDS[coord],BasisSet.num_ao,BasisSet.num_ao,outfile);
    }
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';
  }
  free(outbuf);
  free(label);

  /* print_atommat("Skeleton contribution to the molecular Hessian (a.u.)",Hess); */
  natom = Molecule.num_atoms;
  psio_open(PSIF_DERINFO, PSIO_OPEN_NEW);
  psio_write_entry(PSIF_DERINFO, "Skeleton Hessian", (char *) Hess[0], natom*3*natom*3*sizeof(double));
  psio_close(PSIF_DERINFO, 1);

  for(i=0; i < Molecule.num_atoms*3; i++) {
    free_block(F[i]);
    free_block(S[i]);
    free_block(HDS[i]);
  }
  free(F);  free(S);  free(HDS);
  free_block(Hess);
  cleanup_moinfo();

  return;
}

}}
