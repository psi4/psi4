/*! \file
    \ingroup CIS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libdpd/dpd.h>
#include <libpsio/psio.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cis {

void mp2(void);
void Fab_build(void);
void Fij_build(void);
void Fkc_build(int, int, enum Spin);
void v_build(int, int, enum Spin);
void Z_build(int, int, enum Spin);
int U_build(int, int, double, enum Spin);
void denom(int irrep, double root);

void d_corr(void)
{
  int h, root, errcod, ij, ab;
  char lbl[32];
  double e1, e2_AA, e2_BB, e2_AB, e2;
  double **singlet_d, **singlet_weakp, **triplet_d, **uhf_d, *pair_energy;
  dpdfile2 F, B, B_A, B_B, V, V_A, V_B;
  dpdbuf4 D, U, Z;

  /* compute the first-order ground-state wfn and MP2 energy */
  mp2();

  /* Compute intermediates Fab and Fij */
  Fab_build();
  Fij_build();

  if(params.ref == 0) { /** RHF **/
    singlet_d = (double **) malloc(moinfo.nirreps * sizeof(double *));
    for(h=0; h < moinfo.nirreps; h++)
      singlet_d[h] = init_array(params.rpi[h]);
    moinfo.singlet_d = singlet_d;

    if(params.local) {
      singlet_weakp = (double **) malloc(moinfo.nirreps * sizeof(double *));
      for(h=0; h < moinfo.nirreps; h++)
	singlet_weakp[h] = init_array(params.rpi[h]);
      moinfo.singlet_weakp = singlet_weakp;
      pair_energy = init_array(local.nocc*local.nocc);
    }

    triplet_d = (double **) malloc(moinfo.nirreps * sizeof(double *));
    for(h=0; h < moinfo.nirreps; h++)
      triplet_d[h] = init_array(params.rpi[h]);
    moinfo.triplet_d = triplet_d;
  }
  else if(params.ref == 2) { /** UHF **/
    uhf_d = (double **) malloc(moinfo.nirreps * sizeof(double *));
    for(h=0; h < moinfo.nirreps; h++)
      uhf_d[h] = init_array(params.rpi[h]);
    moinfo.uhf_d = uhf_d;
  }

  /** loop over the desired roots **/
  for(h=0; h < moinfo.nirreps; h++) {

    if(params.ref == 0) { /** RHF **/

      for(root=0; root < params.rpi[h]; root++) {

	/* build the Fkc intermediate for this root */
	Fkc_build(h, root, singlet);

	/* build the v1 intermediate */
	v_build(h, root, singlet);

	/* e1 = b_i^a v_i^a */
	sprintf(lbl, "BIA(%d)[%d] singlet", root, h);
	dpd_file2_init(&B, CC_OEI, h, 0, 1, lbl);
	sprintf(lbl, "VIA[%d]", h);
	dpd_file2_init(&V, CC_MISC, h, 0, 1, lbl);
	e1 = dpd_file2_dot(&V, &B);
	dpd_file2_close(&V);
	dpd_file2_close(&B);

	/* compute the first-order excited-state doubles:
	   u_ij^ab = {<ab||cj> b_i^c - <ab||ci> b_j^c + <ka||ij> b_k^b - <kb||ij> b_k^a}/(D_ij^ab + w) */
	Z_build(h, root, singlet);

	/* generate the necessary denominators */
	denom(h, moinfo.singlet_evals[h][root]);

	errcod = U_build(h, root, moinfo.singlet_evals[h][root], singlet);

	if(!errcod) {

	  /* e2 = 1/4 [ u_ij^ab * ( <ab||cj> b_i^c - <ab||ci> b_j^c + <ka||ij> b_k^b - <kb||ij> b_k^a ) */
	  sprintf(lbl, "UIjAb[%d]", h);
	  dpd_buf4_init(&U, CC_MISC, h, 0, 5, 0, 5, 0, lbl);
	  sprintf(lbl, "(ZIjAb - 1/2 ZIjbA)[%d]", h);
	  dpd_buf4_init(&Z, CC_MISC, h, 0, 5, 0, 5, 0, lbl);
	  e2 = dpd_buf4_dot(&Z, &U);
	  dpd_buf4_close(&U);
	  dpd_buf4_close(&Z);

	  /*	fprintf(outfile, "Singlet: irrep = %d; root = %d; e1 = %20.14f\n", h, root, e1); */
	  /*	  fprintf(outfile, "Singlet: irrep = %d; root = %d; e2 = %20.14f\n", h, root, e2); */
	  singlet_d[h][root] = e1 + e2;

	  if(params.local) { /* compute the weak-pair E2 correction */
	    sprintf(lbl, "UIjAb[%d]", h);
	    dpd_buf4_init(&U, CC_MISC, h, 0, 5, 0, 5, 0, lbl);
	    sprintf(lbl, "(ZIjAb - 1/2 ZIjbA)[%d]", h);
	    dpd_buf4_init(&Z, CC_MISC, h, 0, 5, 0, 5, 0, lbl);

	    dpd_buf4_mat_irrep_init(&U, 0);
	    dpd_buf4_mat_irrep_rd(&U, 0);
	    dpd_buf4_mat_irrep_init(&Z, 0);
	    dpd_buf4_mat_irrep_rd(&Z, 0);

	    singlet_weakp[h][root] = 0.0;
	    for(ij=0; ij < local.nocc*local.nocc; ij++) {
	      if(local.weak_pairs[ij]) {
		for(ab=0; ab < local.nvir*local.nvir; ab++) {
		  singlet_weakp[h][root] += Z.matrix[0][ij][ab] * U.matrix[0][ij][ab];
		}
	      }
	      pair_energy[ij] = 0.0;
	      for(ab=0; ab < local.nvir*local.nvir; ab++)
		pair_energy[ij] += Z.matrix[0][ij][ab] * U.matrix[0][ij][ab];
	    }

	    dpd_buf4_mat_irrep_close(&Z, 0);
	    dpd_buf4_mat_irrep_close(&U, 0);
	    dpd_buf4_close(&Z);
	    dpd_buf4_close(&U);

	    sprintf(lbl, "CIS(D) Pair Energies (%d)", root);
	    psio_write_entry(CC_INFO, lbl, (char *) pair_energy, local.nocc*local.nocc*sizeof(double));
	  }
	}
	else singlet_d[h][root] = 0.0;
      }
      
      for(root=0; root < params.rpi[h]; root++) {

	/*
	  Fkc_build(h, root, triplet);

	  v_build(h, root, triplet);

	  sprintf(lbl, "BIA(%d)[%d] triplet", root, h);
	  dpd_file2_init(&B, CC_OEI, h, 0, 1, lbl);
	  sprintf(lbl, "VIA[%d]", h);
	  dpd_file2_init(&V, CC_MISC, h, 0, 1, lbl);
	  e1 = dpd_file2_dot(&V, &B);
	  dpd_file2_close(&V);
	  dpd_file2_close(&B);

	  Z_build(h, root, triplet);

	  denom(h, moinfo.triplet_evals[h][root]);

	  sprintf(lbl, "ZIjAb[%d]", h);
	  dpd_buf4_init(&Z, CC_MISC, h, 0, 5, 0, 5, 0, lbl);
	  sprintf(lbl, "UIjAb[%d]", h);
	  dpd_buf4_copy(&Z, CC_MISC, lbl);
	  dpd_buf4_close(&Z);

	  sprintf(lbl, "UIjAb[%d]", h);
	  dpd_buf4_init(&U, CC_MISC, h, 0, 5, 0, 5, 0, lbl);
	  sprintf(lbl, "dIjAb[%d]", h);
	  dpd_buf4_init(&D, CC_MISC, h, 0, 5, 0, 5, 0, lbl);
	  dpd_buf4_dirprd(&D, &U);
	  dpd_buf4_close(&D);
	  dpd_buf4_close(&U);

	  sprintf(lbl, "UIjAb[%d]", h);
	  dpd_buf4_init(&U, CC_MISC, h, 0, 5, 0, 5, 0, lbl);
	  sprintf(lbl, "ZIjAb[%d]", h);
	  dpd_buf4_init(&Z, CC_MISC, h, 0, 5, 0, 5, 0, lbl);
	  e2 = dpd_buf4_dot(&Z, &U);
	  dpd_buf4_close(&U);
	  dpd_buf4_close(&Z);
	*/

	/*	fprintf(outfile, "Triplet: irrep = %d; root = %d; e1 = %20.14f\n", h, root, e1); */
	/* fprintf(outfile, "Triplet: irrep = %d; root = %d; e1 = %20.14f\n", h, root, e2); */

	/*	triplet_d[h][root] = e1 + e2; */
      }

    }
    else if(params.ref == 2) { /** UHF **/

      for(root=0; root < params.rpi[h]; root++) {

	/* build the Fkc intermediate for this root */
	Fkc_build(h, root, uhf);

	/* build the v1 intermediate */
	v_build(h, root, uhf);

	/* e1 = b_i^a v_i^a */
	sprintf(lbl, "BIA(%d)[%d]", root, h);
	dpd_file2_init(&B_A, CC_OEI, h, 0, 1, lbl);
	sprintf(lbl, "Bia(%d)[%d]", root, h);
	dpd_file2_init(&B_B, CC_OEI, h, 2, 3, lbl);
	sprintf(lbl, "VIA[%d]", h);
	dpd_file2_init(&V_A, CC_MISC, h, 0, 1, lbl);
	sprintf(lbl, "Via[%d]", h);
	dpd_file2_init(&V_B, CC_MISC, h, 2, 3, lbl);
	e1 = dpd_file2_dot(&V_A, &B_A);
	e1 += dpd_file2_dot(&V_B, &B_B);
	dpd_file2_close(&V_A);
	dpd_file2_close(&V_B);
	dpd_file2_close(&B_A);
	dpd_file2_close(&B_B);

	/* compute the first-order excited-state doubles:
	   z_ij^ab = {<ab||cj> b_i^c - <ab||ci> b_j^c + <ka||ij> b_k^b - <kb||ij> b_k^a}/(D_ij^ab + w) */
	Z_build(h, root, uhf);

	/* generate the necessary denominators */
	denom(h, moinfo.uhf_evals[h][root]);

	U_build(h, root, moinfo.uhf_evals[h][root], uhf);

	/* e2 = 1/4 [ u_ij^ab * ( <ab||cj> b_i^c - <ab||ci> b_j^c + <ka||ij> b_k^b - <kb||ij> b_k^a ) */
	sprintf(lbl, "UIJAB[%d]", h);
	dpd_buf4_init(&U, CC_MISC, h, 2, 7, 2, 7, 0, lbl);
	sprintf(lbl, "ZIJAB[%d]", h);
	dpd_buf4_init(&Z, CC_MISC, h, 2, 7, 2, 7, 0, lbl);
	e2_AA = dpd_buf4_dot(&Z, &U);
	dpd_buf4_close(&U);
	dpd_buf4_close(&Z);

	sprintf(lbl, "Uijab[%d]", h);
	dpd_buf4_init(&U, CC_MISC, h, 12, 17, 12, 17, 0, lbl);
	sprintf(lbl, "Zijab[%d]", h);
	dpd_buf4_init(&Z, CC_MISC, h, 12, 17, 12, 17, 0, lbl);
	e2_BB = dpd_buf4_dot(&Z, &U);
	dpd_buf4_close(&U);
	dpd_buf4_close(&Z);

	sprintf(lbl, "UIjAb[%d]", h);
	dpd_buf4_init(&U, CC_MISC, h, 22, 28, 22, 28, 0, lbl);
	sprintf(lbl, "ZIjAb[%d]", h);
	dpd_buf4_init(&Z, CC_MISC, h, 22, 28, 22, 28, 0, lbl);
	e2_AB = dpd_buf4_dot(&Z, &U);
	dpd_buf4_close(&U);
	dpd_buf4_close(&Z);

	uhf_d[h][root] = e1 + e2_AA + e2_BB + e2_AB;
      }
    }
  }

  if(params.local) free(pair_energy);
}

}} // namespace psi::cis
