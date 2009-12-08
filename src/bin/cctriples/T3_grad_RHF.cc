/*! \file
    \ingroup CCTRIPLES
    \brief Computes T3-dependent terms needed in cclambda and
    ccdensity for (T) contributions to the CCSD(T) energy gradient.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cctriples {

    void T3_grad_RHF(void)
    {
      int h, nirreps;
      int I, J, K, A, B, C, D, L;
      int i, j, k, a, b, c, d, l;
      int ij, ji, ik, ki, jk, kj;
      int ab, ba, ac, ca, bc, cb;
      int il, jl, kl;
      int ad, bd, cd;
      int la, lb, lc;
      int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
      int Gij, Gji, Gik, Gki, Gjk, Gkj, Gijk;
      int Gid, Gjd, Gkd, Gil, Gjl, Gkl;
      int Gab, Gba, Gac, Gca, Gbc, Gcb;
      int nrows, ncols, nlinks;
      int *occpi, *virtpi, *occ_off, *vir_off;
      double t_ia, t_jb, t_kc, D_jkbc, D_ikac, D_ijab;
      double f_ia, f_jb, f_kc, t_jkbc, t_ikac, t_ijab;
      double dijk, denom;
      double ***W0, ***W1, ***M, ***V;
      dpdbuf4 T2, Fints, Eints, Dints, S2, F2ints;
      dpdfile2 fIJ, fAB, fIA, T1, S1;
      int nijk, mijk;
      FILE *ijkfile;
      double value, value1, value2, ET=0.0;
      int id;
      double **Z;

      nirreps = moinfo.nirreps;
      occpi = moinfo.occpi; virtpi = moinfo.virtpi;
      occ_off = moinfo.occ_off;
      vir_off = moinfo.vir_off;

      dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
      dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
      dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
      dpd_file2_mat_init(&fIJ);
      dpd_file2_mat_init(&fAB);
      dpd_file2_mat_init(&fIA);
      dpd_file2_mat_rd(&fIJ);
      dpd_file2_mat_rd(&fAB);
      dpd_file2_mat_rd(&fIA);

      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_mat_init(&T1);
      dpd_file2_mat_rd(&T1);

      dpd_file2_init(&S1, CC_OEI, 0, 0, 1, "SIA");
      dpd_file2_mat_init(&S1);

      dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_buf4_init(&S2, CC_MISC, 0, 0, 5, 0, 5, 0, "SIjAb");
      dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
      dpd_buf4_init(&F2ints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
      dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
      dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_init(&T2, h);
	dpd_buf4_mat_irrep_rd(&T2, h);

	dpd_buf4_mat_irrep_init(&S2, h);

	dpd_buf4_mat_irrep_init(&Eints, h);
	dpd_buf4_mat_irrep_rd(&Eints, h);

	dpd_buf4_mat_irrep_init(&F2ints, h);
	dpd_buf4_mat_irrep_rd(&F2ints, h);

	dpd_buf4_mat_irrep_init(&Dints, h);
	dpd_buf4_mat_irrep_rd(&Dints, h);
      }

      /* Compute the number of IJK combinations */
      /* For now, we need all combinations for gradients */
      nijk = 0;
      for(Gi=0; Gi < nirreps; Gi++)
	for(Gj=0; Gj < nirreps; Gj++)
	  for(Gk=0; Gk < nirreps; Gk++)
	    for(i=0; i < occpi[Gi]; i++) {
	      I = occ_off[Gi] + i;
	      for(j=0; j < occpi[Gj]; j++) {
		J = occ_off[Gj] + j;
		for(k=0; k < occpi[Gk]; k++) {
		  K = occ_off[Gk] + k;

		  nijk++;
		}
	      }
	    }

      ffile(&ijkfile,"ijk.dat", 0);
      fprintf(ijkfile, "Number of IJK combintions: %d\n", nijk);
      fprintf(ijkfile, "\nCurrent IJK Combination: ");

      W0 = (double ***) malloc(nirreps * sizeof(double **));
      W1 = (double ***) malloc(nirreps * sizeof(double **));
      M = (double ***) malloc(nirreps * sizeof(double **));
      V = (double ***) malloc(nirreps * sizeof(double **));

      mijk = 0;
      for(Gi=0; Gi < nirreps; Gi++) {
	for(Gj=0; Gj < nirreps; Gj++) {
	  for(Gk=0; Gk < nirreps; Gk++) {

	    Gkj = Gjk = Gk ^ Gj;
	    Gji = Gij = Gi ^ Gj;
	    Gik = Gki = Gi ^ Gk;

	    Gijk = Gi ^ Gj ^ Gk;

	    for(i=0; i < occpi[Gi]; i++) {
	      I = occ_off[Gi] + i;
	      for(j=0; j < occpi[Gj]; j++) {
		J = occ_off[Gj] + j;
		for(k=0; k < occpi[Gk]; k++) {
		  K = occ_off[Gk] + k;

		  mijk++;
		  fprintf(ijkfile, "%d\n", mijk);
		  fflush(ijkfile);

		  ij = T2.params->rowidx[I][J];
		  ji = T2.params->rowidx[J][I];
		  ik = T2.params->rowidx[I][K];
		  ki = T2.params->rowidx[K][I];
		  jk = T2.params->rowidx[J][K];
		  kj = T2.params->rowidx[K][J];

		  dijk = 0.0;
		  if(fIJ.params->rowtot[Gi])
		    dijk += fIJ.matrix[Gi][i][i];
		  if(fIJ.params->rowtot[Gj])
		    dijk += fIJ.matrix[Gj][j][j];
		  if(fIJ.params->rowtot[Gk])
		    dijk += fIJ.matrix[Gk][k][k];

		  /* Malloc space for the W intermediate */
		  timer_on("malloc");
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;

		    W0[Gab] = dpd_block_matrix(Fints.params->coltot[Gab],virtpi[Gc]);
		    W1[Gab] = dpd_block_matrix(Fints.params->coltot[Gab],virtpi[Gc]);
		    V[Gab] = dpd_block_matrix(Fints.params->coltot[Gab],virtpi[Gc]);
		    M[Gab] = dpd_block_matrix(Fints.params->coltot[Gab], virtpi[Gc]);
		  }
		  timer_off("malloc");

		  /**** Build T3(c) for curent i,j,k;  Result stored in W0 ****/

		  timer_on("N7 Terms");

		  /* +F_idab * t_kjcd */
		  for(Gd=0; Gd < nirreps; Gd++) {

		    Gab = Gid = Gi ^ Gd;
		    Gc = Gkj ^ Gd;

		    /* Set up F integrals */
		    Fints.matrix[Gid] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gid]);
		    dpd_buf4_mat_irrep_rd_block(&Fints, Gid, Fints.row_offset[Gid][I], virtpi[Gd]);

		    /* Set up T2 amplitudes */
		    cd = T2.col_offset[Gkj][Gc];

		    /* Set up multiplication parameters */
		    nrows = Fints.params->coltot[Gid];
		    ncols = virtpi[Gc];
		    nlinks = virtpi[Gd];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			      &(Fints.matrix[Gid][0][0]), nrows, 
			      &(T2.matrix[Gkj][kj][cd]), nlinks, 0.0,
			      &(W0[Gab][0][0]), ncols);

		    dpd_free_block(Fints.matrix[Gid], virtpi[Gd], Fints.params->coltot[Gid]);
		  }

		  /* -E_jklc * t_ilab */
		  for(Gl=0; Gl < nirreps; Gl++) {

		    Gab = Gil = Gi ^ Gl;
		    Gc = Gjk ^ Gl;

		    /* Set up E integrals */
		    lc = Eints.col_offset[Gjk][Gl];

		    /* Set up T2 amplitudes */
		    il = T2.row_offset[Gil][I];

		    /* Set up multiplication parameters */
		    nrows = T2.params->coltot[Gil];
		    ncols = virtpi[Gc];
		    nlinks = occpi[Gl];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			      &(T2.matrix[Gil][il][0]), nrows,
			      &(Eints.matrix[Gjk][jk][lc]), ncols, 1.0,
			      &(W0[Gab][0][0]), ncols);
		  }

		  /* Sort W[ab][c] --> W[ac][b] */
		  dpd_3d_sort(W0, W1, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx, 
			      Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
			      vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, acb, 0);

		  /* +F_idac * t_jkbd */
		  for(Gd=0; Gd < nirreps; Gd++) {

		    Gac = Gid = Gi ^ Gd;
		    Gb = Gjk ^ Gd;

		    Fints.matrix[Gid] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gid]);
		    dpd_buf4_mat_irrep_rd_block(&Fints, Gid, Fints.row_offset[Gid][I], virtpi[Gd]);

		    bd = T2.col_offset[Gjk][Gb];

		    nrows = Fints.params->coltot[Gid];
		    ncols = virtpi[Gb];
		    nlinks = virtpi[Gd];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			      &(Fints.matrix[Gid][0][0]), nrows, 
			      &(T2.matrix[Gjk][jk][bd]), nlinks, 1.0,
			      &(W1[Gac][0][0]), ncols);

		    dpd_free_block(Fints.matrix[Gid], virtpi[Gd], Fints.params->coltot[Gid]);
		  }

		  /* -E_kjlb * t_ilac */
		  for(Gl=0; Gl < nirreps; Gl++) {

		    Gac = Gil = Gi ^ Gl;
		    Gb = Gkj ^ Gl;

		    lb = Eints.col_offset[Gkj][Gl];

		    il = T2.row_offset[Gil][I];

		    nrows = T2.params->coltot[Gil];
		    ncols = virtpi[Gb];
		    nlinks = occpi[Gl];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			      &(T2.matrix[Gil][il][0]), nrows,
			      &(Eints.matrix[Gkj][kj][lb]), ncols, 1.0,
			      &(W1[Gac][0][0]), ncols);
		  }

		  /* Sort W[ac][b] --> W[ca][b] */
		  dpd_3d_sort(W1, W0, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx, 
			      Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
			      vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, bac, 0);

		  /* +F_kdca * t_jibd */
		  for(Gd=0; Gd < nirreps; Gd++) {

		    Gca = Gkd = Gk ^ Gd;
		    Gb = Gji ^ Gd;

		    Fints.matrix[Gkd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gkd]);
		    dpd_buf4_mat_irrep_rd_block(&Fints, Gkd, Fints.row_offset[Gkd][K], virtpi[Gd]);

		    bd = T2.col_offset[Gji][Gb];

		    nrows = Fints.params->coltot[Gkd];
		    ncols = virtpi[Gb];
		    nlinks = virtpi[Gd];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			      &(Fints.matrix[Gkd][0][0]), nrows, 
			      &(T2.matrix[Gji][ji][bd]), nlinks, 1.0,
			      &(W0[Gca][0][0]), ncols);

		    dpd_free_block(Fints.matrix[Gkd], virtpi[Gd], Fints.params->coltot[Gkd]);
		  }

		  /* -E_ijlb * t_klca */
		  for(Gl=0; Gl < nirreps; Gl++) {

		    Gca = Gkl = Gk ^ Gl;
		    Gb = Gij ^ Gl;

		    lb = Eints.col_offset[Gij][Gl];

		    kl = T2.row_offset[Gkl][K];

		    nrows = T2.params->coltot[Gkl];
		    ncols = virtpi[Gb];
		    nlinks = occpi[Gl];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			      &(T2.matrix[Gkl][kl][0]), nrows,
			      &(Eints.matrix[Gij][ij][lb]), ncols, 1.0,
			      &(W0[Gca][0][0]), ncols);
		  }

		  /* Sort W[ca][b] --> W[cb][a] */
		  dpd_3d_sort(W0, W1, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx, 
			      Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
			      vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, acb, 0);

		  /* +F_kdcb * t_ijad */
		  for(Gd=0; Gd < nirreps; Gd++) {

		    Gcb = Gkd = Gk ^ Gd;
		    Ga = Gij ^ Gd;

		    Fints.matrix[Gkd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gkd]);
		    dpd_buf4_mat_irrep_rd_block(&Fints, Gkd, Fints.row_offset[Gkd][K], virtpi[Gd]);

		    ad = T2.col_offset[Gij][Ga];

		    nrows = Fints.params->coltot[Gkd];
		    ncols = virtpi[Ga];
		    nlinks = virtpi[Gd];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			      &(Fints.matrix[Gkd][0][0]), nrows, 
			      &(T2.matrix[Gij][ij][ad]), nlinks, 1.0,
			      &(W1[Gcb][0][0]), ncols);

		    dpd_free_block(Fints.matrix[Gkd], virtpi[Gd], Fints.params->coltot[Gkd]);
		  }

		  /* -E_jila * t_klcb */
		  for(Gl=0; Gl < nirreps; Gl++) {

		    Gcb = Gkl = Gk ^ Gl;
		    Ga = Gji ^ Gl;

		    la = Eints.col_offset[Gji][Gl];

		    kl = T2.row_offset[Gkl][K];

		    nrows = T2.params->coltot[Gkl];
		    ncols = virtpi[Ga];
		    nlinks = occpi[Gl];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			      &(T2.matrix[Gkl][kl][0]), nrows,
			      &(Eints.matrix[Gji][ji][la]), ncols, 1.0,
			      &(W1[Gcb][0][0]), ncols);
		  }

		  /* Sort W[cb][a] --> W[bc][a] */
		  dpd_3d_sort(W1, W0, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx, 
			      Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
			      vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, bac, 0);

		  /* +F_jdbc * t_ikad */
		  for(Gd=0; Gd < nirreps; Gd++) {

		    Gbc = Gjd = Gj ^ Gd;
		    Ga = Gik ^ Gd;

		    Fints.matrix[Gjd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gjd]);
		    dpd_buf4_mat_irrep_rd_block(&Fints, Gjd, Fints.row_offset[Gjd][J], virtpi[Gd]);

		    ad = T2.col_offset[Gik][Ga];

		    nrows = Fints.params->coltot[Gjd];
		    ncols = virtpi[Ga];
		    nlinks = virtpi[Gd];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			      &(Fints.matrix[Gjd][0][0]), nrows, 
			      &(T2.matrix[Gik][ik][ad]), nlinks, 1.0,
			      &(W0[Gbc][0][0]), ncols);

		    dpd_free_block(Fints.matrix[Gjd], virtpi[Gd], Fints.params->coltot[Gjd]);
		  }

		  /* -E_kila * t_jlbc */
		  for(Gl=0; Gl < nirreps; Gl++) {

		    Gbc = Gjl = Gj ^ Gl;
		    Ga = Gki ^ Gl;

		    la = Eints.col_offset[Gki][Gl];

		    jl = T2.row_offset[Gjl][J];

		    nrows = T2.params->coltot[Gjl];
		    ncols = virtpi[Ga];
		    nlinks = occpi[Gl];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			      &(T2.matrix[Gjl][jl][0]), nrows,
			      &(Eints.matrix[Gki][ki][la]), ncols, 1.0,
			      &(W0[Gbc][0][0]), ncols);
		  }

		  /* Sort W[bc][a] --> W[ba][c] */
		  dpd_3d_sort(W0, W1, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx, 
			      Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
			      vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, acb, 0);

		  /* +F_jdba * t_kicd */
		  for(Gd=0; Gd < nirreps; Gd++) {

		    Gba = Gjd = Gj ^ Gd;
		    Gc = Gki ^ Gd;

		    Fints.matrix[Gjd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gjd]);
		    dpd_buf4_mat_irrep_rd_block(&Fints, Gjd, Fints.row_offset[Gjd][J], virtpi[Gd]);

		    cd = T2.col_offset[Gki][Gc];

		    nrows = Fints.params->coltot[Gjd];
		    ncols = virtpi[Gc];
		    nlinks = virtpi[Gd];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			      &(Fints.matrix[Gjd][0][0]), nrows, 
			      &(T2.matrix[Gki][ki][cd]), nlinks, 1.0,
			      &(W1[Gba][0][0]), ncols);

		    dpd_free_block(Fints.matrix[Gjd], virtpi[Gd], Fints.params->coltot[Gjd]);
		  }

		  /* -E_iklc * t_jlba */
		  for(Gl=0; Gl < nirreps; Gl++) {

		    Gba = Gjl = Gj ^ Gl;
		    Gc = Gik ^ Gl;

		    lc = Eints.col_offset[Gik][Gl];

		    jl = T2.row_offset[Gjl][J];

		    nrows = T2.params->coltot[Gjl];
		    ncols = virtpi[Gc];
		    nlinks = occpi[Gl];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			      &(T2.matrix[Gjl][jl][0]), nrows,
			      &(Eints.matrix[Gik][ik][lc]), ncols, 1.0,
			      &(W1[Gba][0][0]), ncols);
		  }

		  /* Sort W[ba][c] --> W[ab][c] */
		  dpd_3d_sort(W1, W0, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx, 
			      Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
			      vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, bac, 0);

		  timer_off("N7 Terms");

		  /**** T3(c) complete ****/


		  /**** Build T3(d) for current i,j,k; Result stored in V ****/

		  timer_on("T3(d) Terms");

		  for(Gab=0; Gab < nirreps; Gab++) {

		    Gc = Gab ^ Gijk;

		    for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {

		      A = Fints.params->colorb[Gab][ab][0];
		      Ga = Fints.params->rsym[A];
		      a = A - vir_off[Ga];
		      B = Fints.params->colorb[Gab][ab][1];
		      Gb = Fints.params->ssym[B];
		      b = B - vir_off[Gb];

		      Gbc = Gb ^ Gc;
		      Gac = Ga ^ Gc;

		      for(c=0; c < virtpi[Gc]; c++) {
			C = vir_off[Gc] + c;

			bc = Dints.params->colidx[B][C];
			ac = Dints.params->colidx[A][C];

			/* +t_ia * D_jkbc + f_ia * t_jkbc */
			if(Gi == Ga && Gjk == Gbc) {
			  t_ia = D_jkbc = 0.0;

			  if(T1.params->rowtot[Gi] && T1.params->coltot[Gi]) {
			    t_ia = T1.matrix[Gi][i][a];
			    f_ia = fIA.matrix[Gi][i][a];
			  }

			  if(Dints.params->rowtot[Gjk] && Dints.params->coltot[Gjk]) {
			    D_jkbc = Dints.matrix[Gjk][jk][bc];
			    t_jkbc = T2.matrix[Gjk][jk][bc];
			  }

			  V[Gab][ab][c] += t_ia * D_jkbc + f_ia * t_jkbc;

			}

			/* +t_jb * D_ikac */
			if(Gj == Gb && Gik == Gac) {
			  t_jb = D_ikac = 0.0;

			  if(T1.params->rowtot[Gj] && T1.params->coltot[Gj]) {
			    t_jb = T1.matrix[Gj][j][b];
			    f_jb = fIA.matrix[Gj][j][b];
			  }

			  if(Dints.params->rowtot[Gik] && Dints.params->coltot[Gik]) {
			    D_ikac = Dints.matrix[Gik][ik][ac];
			    t_ikac = T2.matrix[Gik][ik][ac];
			  }

			  V[Gab][ab][c] += t_jb * D_ikac + f_jb * t_ikac;
			}

			/* +t_kc * D_ijab */
			if(Gk == Gc && Gij == Gab) {
			  t_kc = D_ijab = 0.0;

			  if(T1.params->rowtot[Gk] && T1.params->coltot[Gk]) {
			    t_kc = T1.matrix[Gk][k][c];
			    f_kc = fIA.matrix[Gk][k][c];
			  }

			  if(Dints.params->rowtot[Gij] && Dints.params->coltot[Gij]) {
			    D_ijab = Dints.matrix[Gij][ij][ab];
			    t_ijab = T2.matrix[Gij][ij][ab];
			  }

			  V[Gab][ab][c] += t_kc * D_ijab + f_kc * t_ijab;
			}

		      } /* c */
		    } /* ab */
		  } /* Gab */

		  timer_off("T3(d) Terms");

		  /**** T3(d) complete ****/

		  /**** Compute (T) Energy as a Test ****/
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk; Gba = Gab;
		    for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {
		      A = Fints.params->colorb[Gab][ab][0];
		      Ga = Fints.params->rsym[A];
		      a = A - vir_off[Ga];
		      B = Fints.params->colorb[Gab][ab][1];
		      Gb = Fints.params->ssym[B];
		      b = B - vir_off[Gb];

		      Gac = Gca = Ga ^ Gc;  Gbc = Gcb = Gb ^ Gc;

		      ba = Dints.params->colidx[B][A];

		      for(c=0; c < virtpi[Gc]; c++) {
			C = vir_off[Gc] + c;

			ac = Dints.params->colidx[A][C];
			ca = Dints.params->colidx[C][A];
			bc = Dints.params->colidx[B][C];
			cb = Dints.params->colidx[C][B];

			denom = dijk;
			if(fAB.params->rowtot[Ga])
			  denom -= fAB.matrix[Ga][a][a];
			if(fAB.params->rowtot[Gb])
			  denom -= fAB.matrix[Gb][b][b];
			if(fAB.params->rowtot[Gc])
			  denom -= fAB.matrix[Gc][c][c];

			value1 = W0[Gab][ab][c] + V[Gab][ab][c] - W0[Gcb][cb][a] - V[Gcb][cb][a];
			value2 = 4 * W0[Gab][ab][c] + W0[Gbc][bc][a] + W0[Gca][ca][b];

			ET += value1 * value2 / (3.0 * denom);

		      } /* c */
		    } /* ab */
		  } /* Gab */

		  /**** (T) Energy Contribution Complete ****/

		  /**** Compute S1 needed for lambda equations ****/
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk; Gba = Gab;
		    for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {
		      A = Fints.params->colorb[Gab][ab][0];
		      Ga = Fints.params->rsym[A];
		      a = A - vir_off[Ga];
		      B = Fints.params->colorb[Gab][ab][1];
		      Gb = Fints.params->ssym[B];
		      b = B - vir_off[Gb];

		      Gac = Gca = Ga ^ Gc;  Gbc = Gcb = Gb ^ Gc;

		      ba = Dints.params->colidx[B][A];

		      if(Gi == Ga && S1.params->rowtot[Gi] && S1.params->coltot[Gi]) {
			for(c=0; c < virtpi[Gc]; c++) {
			  C = vir_off[Gc] + c;

			  ac = Dints.params->colidx[A][C];
			  ca = Dints.params->colidx[C][A];
			  bc = Dints.params->colidx[B][C];
			  cb = Dints.params->colidx[C][B];

			  denom = dijk;
			  if(fAB.params->rowtot[Ga])
			    denom -= fAB.matrix[Ga][a][a];
			  if(fAB.params->rowtot[Gb])
			    denom -= fAB.matrix[Gb][b][b];
			  if(fAB.params->rowtot[Gc])
			    denom -= fAB.matrix[Gc][c][c];

			  /* TJ Lee's expression */
			  value = (4 * W0[Gab][ab][c] + W0[Gbc][bc][a] + W0[Gca][ca][b]
				   -3 * W0[Gcb][cb][a] - 2 * W0[Gac][ac][b] - W0[Gba][ba][c])/denom;

			  S1.matrix[Gi][i][a] += 0.5 * Dints.matrix[Gjk][jk][bc] * value;

			} /* c */
		      } /* Gi == Ga && S1 rows and S1 cols */
		    } /* ab */
		  } /* Gab */

		  /**** S1 contributions complete ****/

		  /* Build M3 array */

		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk; Gba = Gab;
		    for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {
		      A = Fints.params->colorb[Gab][ab][0];
		      Ga = Fints.params->rsym[A];
		      a = A - vir_off[Ga];
		      B = Fints.params->colorb[Gab][ab][1];
		      Gb = Fints.params->ssym[B];
		      b = B - vir_off[Gb];
		      Gac = Gca = Ga ^ Gc;
		      Gbc = Gcb = Gb ^ Gc;
		      ba = Dints.params->colidx[B][A];
		      for(c=0; c < virtpi[Gc]; c++) {
			C = vir_off[Gc] + c;
			ac = Dints.params->colidx[A][C];
			ca = Dints.params->colidx[C][A];
			bc = Dints.params->colidx[B][C];
			cb = Dints.params->colidx[C][B];
			denom = dijk;
			if(fAB.params->rowtot[Ga]) denom -= fAB.matrix[Ga][a][a];
			if(fAB.params->rowtot[Gb]) denom -= fAB.matrix[Gb][b][b];
			if(fAB.params->rowtot[Gc]) denom -= fAB.matrix[Gc][c][c];

			/* TJ Lee Expression */
// 			M[Gab][ab][c] = 6 * (8 * W0[Gab][ab][c] + W0[Gbc][bc][a] + W0[Gca][ca][b]
// 					     - 4 * W0[Gcb][cb][a] - 4 * W0[Gac][ac][b] - 4 * W0[Gba][ba][c]
// 					     + 4 * V[Gab][ab][c] + V[Gbc][bc][a] + V[Gca][ca][b]
// 					     - 2 * V[Gcb][cb][a] - 2 * V[Gac][ac][b] - 2 * V[Gba][ba][c])/denom;

                        /* GE Scuseria Expression */
			M[Gab][ab][c] = (4 * (2 * W0[Gab][ab][c] + V[Gab][ab][c]) + (2 * W0[Gbc][bc][a] + V[Gbc][bc][a])
					 + (2 * W0[Gca][ca][b] + V[Gca][ca][b]) - 2 * (2 * W0[Gcb][cb][a] + V[Gcb][cb][a])
					 - (2 * W0[Gac][ac][b] + V[Gac][ac][b]) - 2 * (2 * W0[Gba][ba][c] + V[Gba][ba][c]))/denom;
		      }
		    }
		  }

		  /**** Compute S2 needed for lambda equations ****/
		  /* S_kjcd --> t_ijkabc <id|ab> */
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Gid = Gab = Gi ^ Gd;
		    Gc = Gkj ^ Gd;

		    nrows = virtpi[Gc];
		    ncols = virtpi[Gd];
		    nlinks = Fints.params->coltot[Gid];

		    if(nrows && ncols && nlinks) {
		      id = Fints.row_offset[Gid][I];
		      Fints.matrix[Gid] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gid]);
		      dpd_buf4_mat_irrep_rd_block(&Fints, Gid, id, virtpi[Gd]);
		      cd = S2.col_offset[Gkj][Gc];

		      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, M[Gab][0], nrows, 
			      Fints.matrix[Gid][0], nlinks, 1.0, &(S2.matrix[Gkj][kj][cd]), ncols);

		      dpd_free_block(Fints.matrix[Gid], virtpi[Gd], Fints.params->coltot[Gid]);
		    }
		  } /* Gd */

		  timer_on("malloc");
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    dpd_free_block(W0[Gab],Fints.params->coltot[Gab],virtpi[Gc]);
		    dpd_free_block(W1[Gab],Fints.params->coltot[Gab],virtpi[Gc]);
		    dpd_free_block(V[Gab],Fints.params->coltot[Gab],virtpi[Gc]);
		    dpd_free_block(M[Gab],Fints.params->coltot[Gab],virtpi[Gc]);
		  }
		  timer_off("malloc");

		} /* k */
	      } /* j */
	    } /* i */

	  } /* Gk */
	} /* Gj */
      } /* Gi */

      fprintf(outfile, "\tE(T) = %20.14f\n", ET);

      free(W0); free(W1); free(V); free(M);

      fclose(ijkfile);

      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_wrt(&S2, h);
	dpd_buf4_mat_irrep_close(&S2, h);
	dpd_buf4_mat_irrep_close(&T2, h);
	dpd_buf4_mat_irrep_close(&Eints, h);
	dpd_buf4_mat_irrep_close(&F2ints, h);
	dpd_buf4_mat_irrep_close(&Dints, h);
      }
      dpd_buf4_print(&S2, outfile, 1);
      dpd_buf4_close(&S2);
      dpd_buf4_close(&T2);
      dpd_buf4_close(&Eints);
      dpd_buf4_close(&Dints);
      dpd_buf4_close(&F2ints);
      dpd_buf4_close(&Fints);

      dpd_file2_mat_wrt(&S1);
      dpd_file2_mat_close(&S1);
      dpd_file2_print(&S1, outfile);
      dpd_file2_close(&S1);

      dpd_file2_mat_close(&T1);
      dpd_file2_close(&T1);

      dpd_file2_mat_close(&fIJ);
      dpd_file2_mat_close(&fAB);
      dpd_file2_mat_close(&fIA);
      dpd_file2_close(&fIJ);
      dpd_file2_close(&fAB);
      dpd_file2_close(&fIA);
    }


  }}
