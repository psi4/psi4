/*
 * @BEGIN LICENSE
 *  
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCTRIPLES
    \brief Computes T3-dependent terms needed in cclambda and
    ccdensity for (T) contributions to the CCSD(T) energy gradient.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi { namespace cctriples {

    void T3_grad_RHF(void)
    {
      int h, nirreps;
      int I, J, K, A, B, C, D, L;
      int i, j, k, a, b, c, d, l;
      int ij, ji, ik, ki, jk, kj;
      int ab, ba, ac, ca, bc, cb, dc;
      int la, al, lb, lc;
      int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
      int Gij, Gji, Gik, Gki, Gjk, Gkj, Gijk, Gabc;
      int Gid, Gib, Gjd, Gkd, Gil, Gjl, Gkl;
      int Gab, Gaj, Gba, Gac, Gca, Gbc, Gcb, Gcd;
      int Gad, Gal, Gcl, Gbd, Gbl;
      int ad, kd, lk, il, jl, kl, cd, bd, di, id, ib, dj, jd, dk, aj;
      int nrows, ncols, nlinks;
      int *occpi, *virtpi, *occ_off, *vir_off;
      double t_ia, t_jb, t_kc, D_jkbc, D_ikac, D_ijab;
      double f_ia, f_jb, f_kc, t_jkbc, t_ikac, t_ijab;
      double dijk, dabc, denom;
      double ***W0, ***W1, ***M, ***E, ***V, ***W2, ***Mi, ***RW;
      dpdbuf4 T2, Fints, Eints, Dints, S2, S2i, F2ints, T2i, F3ints, E2ints, si, kappa, chi;
      dpdfile2 fIJ, fAB, fIA, T1, S1, DAB, DIJ;
      int nijk, mijk, nabc, mabc;
      double value, value1, value2, ET1=0.0, ET2=0.0;
      double **Z;

      nirreps = moinfo.nirreps;
      occpi = moinfo.occpi; virtpi = moinfo.virtpi;
      occ_off = moinfo.occ_off;
      vir_off = moinfo.vir_off;

      global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
      global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
      global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
      global_dpd_->file2_mat_init(&fIJ);
      global_dpd_->file2_mat_init(&fAB);
      global_dpd_->file2_mat_init(&fIA);
      global_dpd_->file2_mat_rd(&fIJ);
      global_dpd_->file2_mat_rd(&fAB);
      global_dpd_->file2_mat_rd(&fIA);

      global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
      global_dpd_->file2_mat_init(&T1);
      global_dpd_->file2_mat_rd(&T1);

      global_dpd_->file2_init(&S1, PSIF_CC_OEI, 0, 0, 1, "SIA(T)");
      global_dpd_->file2_mat_init(&S1);

      global_dpd_->file2_init(&DAB, PSIF_CC_OEI, 0, 1, 1, "DAB(T)");
      global_dpd_->file2_mat_init(&DAB);

      global_dpd_->file2_init(&DIJ, PSIF_CC_OEI, 0, 0, 0, "DIJ(T)");
      global_dpd_->file2_mat_init(&DIJ);

      global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      global_dpd_->buf4_sort(&T2, PSIF_CC_TAMPS, rspq, 5, 0, "tAbIj");
      global_dpd_->buf4_close(&T2);

      global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      global_dpd_->buf4_init(&T2i, PSIF_CC_TAMPS, 0, 5, 0, 5, 0, 0, "tAbIj");


      global_dpd_->buf4_init(&S2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "SIjAb(T)");
      global_dpd_->buf4_init(&S2i, PSIF_CC_TMP, 0, 5, 0, 5, 0, 0, "SAbIj(T)");


      global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
      global_dpd_->buf4_sort(&Fints, PSIF_CC_FINTS, rspq, 5, 10, "F <bc|ia>");
      global_dpd_->buf4_close(&Fints);

      global_dpd_->buf4_init(&F3ints, PSIF_CC_FINTS, 0, 5, 10, 5, 10, 0, "F <bc|ia>");
      global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");



      global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
      global_dpd_->buf4_init(&E2ints, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");


      global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

      global_dpd_->buf4_init(&si, PSIF_CC_TMP, 0, 10, 5, 10, 5, 0, "SIbAd(T)");
      global_dpd_->buf4_init(&kappa, PSIF_CC_TMP, 0, 0, 10, 0, 10, 0, "SIjLa(T)");
      global_dpd_->buf4_init(&chi, PSIF_CC_TMP, 0, 0, 5, 0, 5, 0, "SJkBc(T)");


      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&T2, h);
	global_dpd_->buf4_mat_irrep_rd(&T2, h);

	global_dpd_->buf4_mat_irrep_init(&T2i, h);
	global_dpd_->buf4_mat_irrep_rd(&T2i, h);

	global_dpd_->buf4_mat_irrep_init(&S2, h);
	global_dpd_->buf4_mat_irrep_init(&S2i, h);

	global_dpd_->buf4_mat_irrep_init(&Eints, h);
	global_dpd_->buf4_mat_irrep_rd(&Eints, h);

	global_dpd_->buf4_mat_irrep_init(&E2ints, h);
	global_dpd_->buf4_mat_irrep_rd(&E2ints, h);
	
	global_dpd_->buf4_mat_irrep_init(&F3ints, h);
	global_dpd_->buf4_mat_irrep_rd(&F3ints, h);

	global_dpd_->buf4_mat_irrep_init(&Fints, h);
	global_dpd_->buf4_mat_irrep_rd(&Fints, h);

	global_dpd_->buf4_mat_irrep_init(&Dints, h);
	global_dpd_->buf4_mat_irrep_rd(&Dints, h);

	global_dpd_->buf4_mat_irrep_init(&si, h);
	
	global_dpd_->buf4_mat_irrep_init(&kappa, h);
	
	global_dpd_->buf4_mat_irrep_init(&chi, h);
      }

      /* Compute the number of IJK combinations */
      /* For now, we need all combinations for gradients */
      /*nijk = 0;
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
      //boost::shared_ptr<OutFile> printer(new OutFile("ijk.dat",TRUNCATE));
      //ffile(&ijkfile,"ijk.dat", 0);
      //printer->Printf( "Number of IJK combinations: %d\n", nijk);
      //printer->Printf( "\nCurrent IJK Combination: ");*/

      W0 = (double ***) malloc(nirreps * sizeof(double **));
      W1 = (double ***) malloc(nirreps * sizeof(double **));
      W2 = (double ***) malloc(nirreps * sizeof(double **));
      M = (double ***) malloc(nirreps * sizeof(double **));
      E = (double ***) malloc(nirreps * sizeof(double **));
      V = (double ***) malloc(nirreps * sizeof(double **));
      Mi = (double ***) malloc(nirreps * sizeof(double **));
      RW = (double ***) malloc(nirreps * sizeof(double **));

      //mijk = 0;
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

		  //mijk++;
		  //printer->Printf( "%d\n", mijk);


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
		  timer_on("malloc)");
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;

		    W0[Gab] = global_dpd_->dpd_block_matrix(Fints.params->coltot[Gab],virtpi[Gc]);
		    W1[Gab] = global_dpd_->dpd_block_matrix(Fints.params->coltot[Gab],virtpi[Gc]);
		    V[Gab] = global_dpd_->dpd_block_matrix(Fints.params->coltot[Gab],virtpi[Gc]);
		    M[Gab] = global_dpd_->dpd_block_matrix(Fints.params->coltot[Gab], virtpi[Gc]);
		    E[Gab] = global_dpd_->dpd_block_matrix(Fints.params->coltot[Gab], virtpi[Gc]);
		  }
		  for(Ga=0; Ga < nirreps; Ga++) {
              	   Gbc = Ga ^ Gijk;
                   Mi[Ga] = global_dpd_->dpd_block_matrix(virtpi[Ga],Fints.params->coltot[Gbc]);
                   RW[Ga] = global_dpd_->dpd_block_matrix(virtpi[Ga],Fints.params->coltot[Gbc]);
                  }

		  timer_off("malloc)");

		  /**** Build T3(c) for curent i,j,k;  Result stored in W0 ****/

		  timer_on("N7 Terms(ijk)");

		  /* +F_idab * t_kjcd */
		  for(Gd=0; Gd < nirreps; Gd++) {

		    Gab = Gid = Gi ^ Gd;
		    Gc = Gkj ^ Gd;

		    /* Set up F integrals */
		    Fints.matrix[Gid] = global_dpd_->dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gid]);
		    global_dpd_->buf4_mat_irrep_rd_block(&Fints, Gid, Fints.row_offset[Gid][I], virtpi[Gd]);

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

		    global_dpd_->free_dpd_block(Fints.matrix[Gid], virtpi[Gd], Fints.params->coltot[Gid]);
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
          global_dpd_->sort_3d(W0, W1, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx,
			      Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
			      vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, acb, 0);

		  /* +F_idac * t_jkbd */
		  for(Gd=0; Gd < nirreps; Gd++) {

		    Gac = Gid = Gi ^ Gd;
		    Gb = Gjk ^ Gd;

		    Fints.matrix[Gid] = global_dpd_->dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gid]);
		    global_dpd_->buf4_mat_irrep_rd_block(&Fints, Gid, Fints.row_offset[Gid][I], virtpi[Gd]);

		    bd = T2.col_offset[Gjk][Gb];

		    nrows = Fints.params->coltot[Gid];
		    ncols = virtpi[Gb];
		    nlinks = virtpi[Gd];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			      &(Fints.matrix[Gid][0][0]), nrows, 
			      &(T2.matrix[Gjk][jk][bd]), nlinks, 1.0,
			      &(W1[Gac][0][0]), ncols);

		    global_dpd_->free_dpd_block(Fints.matrix[Gid], virtpi[Gd], Fints.params->coltot[Gid]);
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
          global_dpd_->sort_3d(W1, W0, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx,
			      Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
			      vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, bac, 0);

		  /* +F_kdca * t_jibd */
		  for(Gd=0; Gd < nirreps; Gd++) {

		    Gca = Gkd = Gk ^ Gd;
		    Gb = Gji ^ Gd;

		    Fints.matrix[Gkd] = global_dpd_->dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gkd]);
		    global_dpd_->buf4_mat_irrep_rd_block(&Fints, Gkd, Fints.row_offset[Gkd][K], virtpi[Gd]);

		    bd = T2.col_offset[Gji][Gb];

		    nrows = Fints.params->coltot[Gkd];
		    ncols = virtpi[Gb];
		    nlinks = virtpi[Gd];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			      &(Fints.matrix[Gkd][0][0]), nrows, 
			      &(T2.matrix[Gji][ji][bd]), nlinks, 1.0,
			      &(W0[Gca][0][0]), ncols);

		    global_dpd_->free_dpd_block(Fints.matrix[Gkd], virtpi[Gd], Fints.params->coltot[Gkd]);
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
          global_dpd_->sort_3d(W0, W1, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx,
			      Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
			      vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, acb, 0);

		  /* +F_kdcb * t_ijad */
		  for(Gd=0; Gd < nirreps; Gd++) {

		    Gcb = Gkd = Gk ^ Gd;
		    Ga = Gij ^ Gd;

		    Fints.matrix[Gkd] = global_dpd_->dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gkd]);
		    global_dpd_->buf4_mat_irrep_rd_block(&Fints, Gkd, Fints.row_offset[Gkd][K], virtpi[Gd]);

		    ad = T2.col_offset[Gij][Ga];

		    nrows = Fints.params->coltot[Gkd];
		    ncols = virtpi[Ga];
		    nlinks = virtpi[Gd];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			      &(Fints.matrix[Gkd][0][0]), nrows, 
			      &(T2.matrix[Gij][ij][ad]), nlinks, 1.0,
			      &(W1[Gcb][0][0]), ncols);

		    global_dpd_->free_dpd_block(Fints.matrix[Gkd], virtpi[Gd], Fints.params->coltot[Gkd]);
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
          global_dpd_->sort_3d(W1, W0, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx,
			      Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
			      vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, bac, 0);

		  /* +F_jdbc * t_ikad */
		  for(Gd=0; Gd < nirreps; Gd++) {

		    Gbc = Gjd = Gj ^ Gd;
		    Ga = Gik ^ Gd;

		    Fints.matrix[Gjd] = global_dpd_->dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gjd]);
		    global_dpd_->buf4_mat_irrep_rd_block(&Fints, Gjd, Fints.row_offset[Gjd][J], virtpi[Gd]);

		    ad = T2.col_offset[Gik][Ga];

		    nrows = Fints.params->coltot[Gjd];
		    ncols = virtpi[Ga];
		    nlinks = virtpi[Gd];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			      &(Fints.matrix[Gjd][0][0]), nrows, 
			      &(T2.matrix[Gik][ik][ad]), nlinks, 1.0,
			      &(W0[Gbc][0][0]), ncols);

		    global_dpd_->free_dpd_block(Fints.matrix[Gjd], virtpi[Gd], Fints.params->coltot[Gjd]);
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
          global_dpd_->sort_3d(W0, W1, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx,
			      Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
			      vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, acb, 0);

		  /* +F_jdba * t_kicd */
		  for(Gd=0; Gd < nirreps; Gd++) {

		    Gba = Gjd = Gj ^ Gd;
		    Gc = Gki ^ Gd;

		    Fints.matrix[Gjd] = global_dpd_->dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gjd]);
		    global_dpd_->buf4_mat_irrep_rd_block(&Fints, Gjd, Fints.row_offset[Gjd][J], virtpi[Gd]);

		    cd = T2.col_offset[Gki][Gc];

		    nrows = Fints.params->coltot[Gjd];
		    ncols = virtpi[Gc];
		    nlinks = virtpi[Gd];

		    if(nrows && ncols && nlinks)
		      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			      &(Fints.matrix[Gjd][0][0]), nrows, 
			      &(T2.matrix[Gki][ki][cd]), nlinks, 1.0,
			      &(W1[Gba][0][0]), ncols);

		    global_dpd_->free_dpd_block(Fints.matrix[Gjd], virtpi[Gd], Fints.params->coltot[Gjd]);
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
          global_dpd_->sort_3d(W1, W0, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx,
			      Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
			      vir_off, vir_off, virtpi, vir_off, Fints.params->colidx, bac, 0);

      	          timer_off("N7 Terms(ijk)");

		  /**** T3(c) complete ****/


		  /**** Build T3(d) for current i,j,k; Result stored in V ****/

		  timer_on("T3(d) Terms(ijk)");

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

		  timer_off("T3(d) Terms(ijk)");

		  /**** T3(d) complete ****/

		  /**** Compute (T) Energy as a Test ****/
		  /*for(Gab=0; Gab < nirreps; Gab++) {
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
			ET1 += value1 * value2 / (3.0 * denom);

		      } // c 
		    } // ab 
		  } // Gab 
	          */	
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

			  value = (4 * W0[Gab][ab][c] + W0[Gbc][bc][a] + W0[Gca][ca][b]
			  	   -3 * W0[Gcb][cb][a] - 2 * W0[Gac][ac][b] - W0[Gba][ba][c])/denom;
			 
			   S1.matrix[Gi][i][a] += 0.5 * Dints.matrix[Gjk][jk][bc] * value;

			} /* c */
		      } /* Gi == Ga && S1 rows and S1 cols */
		    } /* ab */
		  } /* Gab */

		  /**** S1 contributions complete ****/

		  /* Build M3, E3 and chi arrays */

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
 			M[Gab][ab][c] = 2 * (8 * W0[Gab][ab][c] + 2 * W0[Gbc][bc][a] + 2 * W0[Gca][ca][b]
 					     - 4 * W0[Gcb][cb][a] - 4 * W0[Gac][ac][b] - 4 * W0[Gba][ba][c]
 					     + 4 * V[Gab][ab][c] + V[Gbc][bc][a] + V[Gca][ca][b]
 					     - 2 * V[Gcb][cb][a] - 2 * V[Gac][ac][b] - 2 * V[Gba][ba][c])/denom;


                        E[Gab][ab][c] = (8 * W0[Gab][ab][c] + 2 * W0[Gbc][bc][a] + 2 * W0[Gca][ca][b]
                                             - 4 * W0[Gcb][cb][a] - 4 * W0[Gac][ac][b] - 4 * W0[Gba][ba][c]
                                             + 8 * V[Gab][ab][c] + 2 * V[Gbc][bc][a] + 2 * V[Gca][ca][b]
                                             - 4 * V[Gcb][cb][a] - 4 * V[Gac][ac][b] - 4 * V[Gba][ba][c])/(denom * denom);

			Mi[Ga][a][cb] = M[Gab][ab][c];

			RW[Ga][a][bc] = (4 * W0[Gab][ab][c] +  W0[Gbc][bc][a] + W0[Gca][ca][b] - 3 * W0[Gcb][cb][a] - 2 * W0[Gac][ac][b] - W0[Gba][ba][c])/denom ;


		      }
		    }
		  }


                for(Gbc=0; Gbc < nirreps; Gbc++) {
                    Ga = Gbc ^ Gijk;
                    if(Gi == Ga) {
                      for(bc=0; bc < Fints.params->coltot[Gbc]; bc++) {
                        for(a=0; a < virtpi[Ga]; a++) {
                          A = vir_off[Ga] + a;
                          if(T1.params->rowtot[Gi] && T1.params->coltot[Gi])
                            chi.matrix[Gjk][jk][bc] += 0.5 * RW[Ga][a][bc] * T1.matrix[Gi][i][a];
                        }
                      }
                    }

                  } /* Gab */

		  /**** Compute S2 needed for lambda equations ****/
		  /* S_kjcd --> M_ijkabc <id|ab>  */

		
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Gid = Gab = Gi ^ Gd;
		    Gc = Gkj ^ Gd;

		    nrows = virtpi[Gc];
		    ncols = virtpi[Gd];
		    nlinks = Fints.params->coltot[Gid];

		    if(nrows && ncols && nlinks) {
		      id = Fints.row_offset[Gid][I];
		      Fints.matrix[Gid] = global_dpd_->dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gid]);
		      global_dpd_->buf4_mat_irrep_rd_block(&Fints, Gid, id, virtpi[Gd]);
		      cd = S2.col_offset[Gkj][Gc];

		      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, M[Gab][0], nrows, 
			      Fints.matrix[Gid][0], nlinks, 1.0, &(S2.matrix[Gkj][kj][cd]), ncols);

		      global_dpd_->free_dpd_block(Fints.matrix[Gid], virtpi[Gd], Fints.params->coltot[Gid]);
		    }
		  } /* Gd */

		/* E3 --> D_ab += e_a\delta_ab */
		
		 for (Gab =0; Gab < nirreps; Gab++){
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < Fints.params->coltot[Gab]; ab++){
		        A = Fints.params->colorb[Gab][ab][0];
                        Ga = Fints.params->rsym[A];
                        a = A - vir_off[Ga];
		        for(c=0; c < virtpi[Gc]; c++){
                          DAB.matrix[Ga][a][a] += 0.5 * W0[Gab][ab][c] * E[Gab][ab][c]	;	
		    }
		}
	    }		


		/* M3 --> si(id,ab) = \sum_jkc M[ab][c](ijk) * t_kjcd*/
		for(Gc=0; Gc < nirreps; Gc++){
		    Gd = Gkj ^ Gc;
		    Gab = Gc ^ Gijk;
 
		    id = si.row_offset[Gab][I]; 	
		    cd = T2.col_offset[Gkj][Gc];
		    nrows = virtpi[Gd];
		    ncols = si.params->coltot[Gab];
		    nlinks = virtpi[Gc];	
                  
		    if (nrows && ncols && nlinks)	 
		     C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, &(T2.matrix[Gkj][kj][cd]), nrows,
                              M[Gab][0], nlinks, 1.0, &(si.matrix[Gab][id][0]), ncols);

			}

                  /* M3 --> kappa(ji,la) = \sum_kbc M[a][bc](ijk) * t_klcb*/
		  /* t_klcb  * Mi[a][cb]  */
		  for(Gcb=0; Gcb < nirreps; Gcb++){
		     Ga = Gcb ^ Gijk;
		     Gl = Gk ^ Gcb;
		     Gkl = Gcb;	
		     	
		     kl = T2.row_offset[Gkl][K];
                     la = kappa.col_offset[Gji][Gl];
		     nrows = occpi[Gl]; 	
	             ncols = virtpi[Ga];
		     nlinks = T2.params->coltot[Gcb];

		    if (nrows && ncols && nlinks)	 
		     C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, &(T2.matrix[Gkl][kl][0]), nlinks, 
                              Mi[Ga][0], nlinks, 1.0, &(kappa.matrix[Gji][ji][la]), ncols);

		}

			
		  timer_on("malloc");
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    global_dpd_->free_dpd_block(W0[Gab],Fints.params->coltot[Gab],virtpi[Gc]);
		    global_dpd_->free_dpd_block(W1[Gab],Fints.params->coltot[Gab],virtpi[Gc]);
		    global_dpd_->free_dpd_block(V[Gab],Fints.params->coltot[Gab],virtpi[Gc]);
		    global_dpd_->free_dpd_block(M[Gab],Fints.params->coltot[Gab],virtpi[Gc]);
		    global_dpd_->free_dpd_block(E[Gab],Fints.params->coltot[Gab],virtpi[Gc]);
		  }
                  for(Ga=0; Ga < nirreps; Ga++) {
                     Gbc = Ga ^ Gijk;
                     global_dpd_->free_dpd_block(Mi[Ga], virtpi[Ga],Fints.params->coltot[Gbc]);
                     global_dpd_->free_dpd_block(RW[Ga], virtpi[Ga],Fints.params->coltot[Gbc]);
                     }
		  timer_off("malloc");

		} /* k */
	      } /* j */
	    } /* i */

	  } /* Gk */
	} /* Gj */
      } /* Gi */

      //outfile->Printf( "\tE(T)(ijk) = %20.14f\n", ET1);

      for(h=0; h < nirreps; h++) {
   	     global_dpd_->buf4_mat_irrep_wrt(&S2, h);
       	     global_dpd_->buf4_mat_irrep_close(&S2, h);
	}
       	     global_dpd_->buf4_close(&S2);

      for(h=0; h < nirreps; h++) {
       	     global_dpd_->buf4_mat_irrep_wrt(&si, h);
       	     global_dpd_->buf4_mat_irrep_close(&si, h);
	}
             global_dpd_->buf4_close(&si);

      for(h=0; h < nirreps; h++) {
      	     global_dpd_->buf4_mat_irrep_wrt(&kappa, h);
       	     global_dpd_->buf4_mat_irrep_close(&kappa, h);
	}
             global_dpd_->buf4_close(&kappa);

      for(h=0; h < nirreps; h++) {
       	     global_dpd_->buf4_mat_irrep_wrt(&chi, h);
       	     global_dpd_->buf4_mat_irrep_close(&chi, h);
	}
             global_dpd_->buf4_close(&chi);

      global_dpd_->file2_mat_wrt(&S1);
      global_dpd_->file2_mat_close(&S1);
      global_dpd_->file2_close(&S1);

      global_dpd_->file2_mat_wrt(&DAB);
      global_dpd_->file2_mat_close(&DAB);
      global_dpd_->file2_close(&DAB);

	/* Compute the number of ABC combinations */
      	/* For now, we need all combinations for gradients */
      /*nabc = 0;
      for(Ga=0; Ga < nirreps; Ga++)
        for(Gb=0; Gb < nirreps; Gb++)
          for(Gc=0; Gc < nirreps; Gc++)
            for(a=0; a < virtpi[Ga]; a++) {
              A = vir_off[Ga] + a;
              for(b=0; b < virtpi[Gb]; b++) {
                B = vir_off[Gb] + b;
                for(c=0; c < virtpi[Gc]; c++) {
                  C = vir_off[Gc] + c;

                  nabc++;
                }
              }
            }

      boost::shared_ptr<OutFile> printer1(new OutFile("abc.dat",TRUNCATE));
      //ffile(&abcfile,"abc.dat", 0);
      printer1->Printf( "Number of ABC combinations: %d\n", nabc);
      printer1->Printf( "\nCurrent ABC Combination: ");*/
    

      value1 = 0.0;
      value2 = 0.0;

      //mabc = 0;
      for(Ga=0; Ga < nirreps; Ga++) {
        for(Gb=0; Gb < nirreps; Gb++) {
          for(Gc=0; Gc < nirreps; Gc++) {

            Gcb = Gbc = Gc ^ Gb;
            Gba = Gab = Ga ^ Gb;
            Gac = Gca = Ga ^ Gc;

            Gabc = Ga ^ Gb ^ Gc;

            for(a=0; a < virtpi[Ga]; a++) {
              A = vir_off[Ga] + a;
              for(b=0; b < virtpi[Gb]; b++) {
                B = vir_off[Gb] + b;
                for(c=0; c < virtpi[Gc]; c++) {
                  C = vir_off[Gc] + c;

                  //mabc++;
                  //printer1->Printf( "%d\n", mabc);


                  ab = F3ints.params->rowidx[A][B];
                  ba = F3ints.params->rowidx[B][A];
                  ac = F3ints.params->rowidx[A][C];
                  ca = F3ints.params->rowidx[C][A];
                  bc = F3ints.params->rowidx[B][C];
                  cb = F3ints.params->rowidx[C][B];

                  dabc = 0.0;
                  if(fAB.params->rowtot[Ga])
                    dabc -= fAB.matrix[Ga][a][a];
                  if(fAB.params->rowtot[Gb])
                    dabc -= fAB.matrix[Gb][b][b];
                  if(fAB.params->rowtot[Gc])
                    dabc -= fAB.matrix[Gc][c][c];

		

                  /* Malloc space for the W intermediate */
                  timer_on("malloc");
                  for(Gij=0; Gij < nirreps; Gij++) {
                    Gk = Gij ^ Gabc;
                    W0[Gij] = global_dpd_->dpd_block_matrix(Eints.params->rowtot[Gij], occpi[Gk]);
                    W1[Gij] = global_dpd_->dpd_block_matrix(Eints.params->rowtot[Gij], occpi[Gk]);
                    V[Gij] = global_dpd_->dpd_block_matrix(Eints.params->rowtot[Gij], occpi[Gk]);
                    M[Gij] = global_dpd_->dpd_block_matrix(Eints.params->rowtot[Gij], occpi[Gk]);
                    E[Gij] = global_dpd_->dpd_block_matrix(Eints.params->rowtot[Gij], occpi[Gk]);
                  }
                  timer_off("malloc");

              /**** Build T3(c) for curent a,b,c  Result stored in W0 ****/

                  timer_on("N7 Terms(abc)");


		/* t_adij F_cbkd */
                  for(Gd=0; Gd < nirreps; Gd++) {

                    Gij = Gad = Ga ^ Gd;
                    Gk = Gd ^ Gcb;

                 /* Set up  T2i amplitudes */ 
                  ad = T2i.row_offset[Gad][A];        

                 /* Set up  F3ints integrals*/ 
                  kd = F3ints.col_offset[Gcb][Gk];

                    /* Set up multiplication parameters */
                    nrows = T2i.params->coltot[Gad];
                    ncols = occpi[Gk];
                    nlinks = virtpi[Gd];

                    if(nrows && ncols && nlinks && W0[Gij] != NULL)
                      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                              &(T2i.matrix[Gad][ad][0]), nrows,
                              &(F3ints.matrix[Gcb][cb][kd]), nlinks, 0.0,
                              &(W0[Gij][0][0]), ncols);
              }

              /* -E_alij t_cbkl */

                for(Gl=0; Gl < nirreps; Gl++){
                                      
                  Gij = Gal = Ga ^ Gl;
                  Gk = Gl ^ Gcb;
              
              /* Set up  T2i amplitudes     */        
                   kl = T2i.col_offset[Gcb][Gk];
      
                /* Set up multiplication parameters */
                    nrows = E2ints.params->coltot[Gal];
                    ncols = occpi[Gk];
                    nlinks = occpi[Gl];

                    if(nrows && ncols && nlinks && W0[Gij] != NULL){

            	  E2ints.matrix[Gal] = global_dpd_->dpd_block_matrix(occpi[Gl], E2ints.params->coltot[Gal]);
                  global_dpd_->buf4_mat_irrep_rd_block(&E2ints, Gal, E2ints.row_offset[Gal][A], occpi[Gl]);   

                      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                              &(E2ints.matrix[Gal][0][0]), nrows,
                              &(T2i.matrix[Gcb][cb][kl]), nlinks, 1.0,
                              &(W0[Gij][0][0]), ncols);

                    global_dpd_->free_dpd_block(E2ints.matrix[Gal], occpi[Gl], E2ints.params->coltot[Gal]);  
		}
    }
              /* Sort W[ij][k] --> W[ik][j] */
         global_dpd_->sort_3d(W0, W1, nirreps, Gabc, Eints.params->rowtot, Eints.params->rowidx,
                              Eints.params->roworb, Eints.params->psym, Eints.params->qsym,
                              occ_off, occ_off, occpi, occ_off, Eints.params->rowidx, acb, 0);

                  /* t_adik F_bcjd */
                  for(Gd=0; Gd < nirreps; Gd++) {

                    Gik = Gad = Ga ^ Gd;
                    Gj = Gd ^ Gbc;

                   /* Set up  T2i amplitudes */
                  ad = T2i.row_offset[Gad][A];

                   /* Set up  F3ints integrals*/
                    jd = F3ints.col_offset[Gbc][Gj];

                    /* Set up multiplication parameters */
                    nrows = T2i.params->coltot[Gad];
                    ncols = occpi[Gj];
                    nlinks = virtpi[Gd];

                    if(nrows && ncols && nlinks && W1[Gik] != NULL)
                      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                              &(T2i.matrix[Gad][ad][0]), nrows,
                              &(F3ints.matrix[Gbc][bc][jd]), nlinks, 1.0,
                              &(W1[Gik][0][0]), ncols);

                }

                /* -E_alik t_bcjl */

                  for(Gl=0; Gl < nirreps; Gl++){

                    Gik = Gal = Ga ^ Gl;
                    Gj = Gl ^ Gbc;

                /* Set up  T2i amplitudes     */
                    jl = T2i.col_offset[Gbc][Gj];

                /* Set up multiplication parameters */
                    nrows = E2ints.params->coltot[Gal];
                    ncols = occpi[Gj];
                    nlinks = occpi[Gl];

 			if(nrows && ncols && nlinks && W1[Gik] != NULL){

                    E2ints.matrix[Gal] = global_dpd_->dpd_block_matrix(occpi[Gl], E2ints.params->coltot[Gal]);
                    global_dpd_->buf4_mat_irrep_rd_block(&E2ints, Gal, E2ints.row_offset[Gal][A], occpi[Gl]);

                      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                              &(E2ints.matrix[Gal][0][0]), nrows,
                              &(T2i.matrix[Gbc][bc][jl]), nlinks, 1.0,
                              &(W1[Gik][0][0]), ncols);

                    global_dpd_->free_dpd_block(E2ints.matrix[Gal], occpi[Gl], E2ints.params->coltot[Gal]);
			}
                   }
                /* Sort W[ik][j] --> W[ki][j] */
         global_dpd_->sort_3d(W1, W0, nirreps, Gabc, Eints.params->rowtot, Eints.params->rowidx,
                              Eints.params->roworb, Eints.params->psym, Eints.params->qsym,
                              occ_off, occ_off, occpi, occ_off, Eints.params->rowidx, bac, 0);

              /* t_cdki Fbajd */

                  for(Gd=0; Gd < nirreps; Gd++) {

                    Gki = Gcd = Gc ^ Gd;
                    Gj = Gd ^ Gba;

                   /* Set up T2i amplitudes*/
                    cd = T2i.row_offset[Gcd][C];
                   /* Set up  F3ints integrals*/
                    jd = F3ints.col_offset[Gba][Gj];

                    /* Set up multiplication parameters */
                    nrows = T2i.params->coltot[Gcd];
                    ncols = occpi[Gj];
                    nlinks = virtpi[Gd];

                    if(nrows && ncols && nlinks && W0[Gki] != NULL)
                      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                              &(T2i.matrix[Gcd][cd][0]), nrows,
                              &(F3ints.matrix[Gba][ba][jd]), nlinks, 1.0,
                              &(W0[Gki][0][0]), ncols);

                }

                /* -E_clki t_bajl */

                  for(Gl=0; Gl < nirreps; Gl++){

                    Gki = Gcl = Gc ^ Gl;
                    Gj = Gl ^ Gba;

		/*	* Set up  T2i amplitudes     */
                    jl = T2i.col_offset[Gba][Gj];

                /* Set up multiplication parameters */
                    nrows = E2ints.params->coltot[Gcl];
                    ncols = occpi[Gj];
                    nlinks = occpi[Gl];

                    if(nrows && ncols && nlinks && W0[Gki] != NULL){

                    E2ints.matrix[Gcl] = global_dpd_->dpd_block_matrix(occpi[Gl], E2ints.params->coltot[Gcl]);
                    global_dpd_->buf4_mat_irrep_rd_block(&E2ints, Gcl, E2ints.row_offset[Gcl][C], occpi[Gl]);

                      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                              &(E2ints.matrix[Gcl][0][0]), nrows,
                              &(T2i.matrix[Gba][ba][jl]), nlinks, 1.0,
                              &(W0[Gki][0][0]), ncols);

                    global_dpd_->free_dpd_block(E2ints.matrix[Gcl], occpi[Gl], E2ints.params->coltot[Gcl]);
		}
                        }
                /* Sort W[ki][j] --> W[kj][i] */
         global_dpd_->sort_3d(W0, W1, nirreps, Gabc, Eints.params->rowtot, Eints.params->rowidx,
                              Eints.params->roworb, Eints.params->psym, Eints.params->qsym,
                              occ_off, occ_off, occpi, occ_off, Eints.params->rowidx, acb, 0);


                /* t_cdkj * F_abid */

		for(Gd=0; Gd < nirreps; Gd++) {

         	     Gkj = Gcd = Gc ^ Gd;
                     Gi = Gd ^ Gab;
 
                    /* Set up  T2i amplitudes*/
                    cd = T2i.row_offset[Gcd][C];
 
                    /* Set up F3ints integrals*/
                     id = F3ints.col_offset[Gab][Gi];
 
                     /* Set up multiplication parameters */
                     nrows = T2i.params->coltot[Gcd];
                     ncols = occpi[Gi];
                     nlinks = virtpi[Gd];
 
 
                     if(nrows && ncols && nlinks && W1[Gkj] != NULL)
                       C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                               &(T2i.matrix[Gcd][cd][0]), nrows,
                               &(F3ints.matrix[Gab][ab][id]), nlinks, 1.0, 
                               &(W1[Gkj][0][0]), ncols);
                 }
 
               /* -E_clkj * t_abil */
 
                   for(Gl=0; Gl < nirreps; Gl++){
 
                     Gkj = Gcl = Gc ^ Gl;
                     Gi = Gl ^ Gab;
 
                 /* Set up  T2i amplitudes     */
                     il = T2i.col_offset[Gab][Gi];
 
                 /* Set up multiplication parameters */
                     nrows = E2ints.params->coltot[Gcl];
                     ncols = occpi[Gi];
                     nlinks = occpi[Gl];
 
                     if(nrows && ncols && nlinks && W1[Gkj] != NULL){

                     E2ints.matrix[Gcl] = global_dpd_->dpd_block_matrix(occpi[Gl], E2ints.params->coltot[Gcl]);
                     global_dpd_->buf4_mat_irrep_rd_block(&E2ints, Gcl, E2ints.row_offset[Gcl][C], occpi[Gl]);
 
                       C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                               &(E2ints.matrix[Gcl][0][0]), nrows,
                               &(T2i.matrix[Gab][ab][il]), nlinks, 1.0,
                               &(W1[Gkj][0][0]), ncols);
 
                     global_dpd_->free_dpd_block(E2ints.matrix[Gcl], occpi[Gl], E2ints.params->coltot[Gcl]);
                           }
                 }


         /* Sort W[kj][i] --> W[jk][i] */
         global_dpd_->sort_3d(W1, W0, nirreps, Gabc, Eints.params->rowtot, Eints.params->rowidx,
                              Eints.params->roworb, Eints.params->psym, Eints.params->qsym,
                              occ_off, occ_off, occpi, occ_off, Eints.params->rowidx, bac, 0);

               /* t_bdjk * F_acid */

                  for(Gd=0; Gd < nirreps; Gd++) {

                    Gjk = Gbd = Gb ^ Gd;
                    Gi = Gd ^ Gac;

                   /* Set up  T2i amplitudes*/
                   bd = T2i.row_offset[Gbd][B];       

                   /* Set up  F3ints integrals*/
                     id = F3ints.col_offset[Gac][Gi];

                    /* Set up multiplication parameters */
                    nrows = T2i.params->coltot[Gbd];
                    ncols = occpi[Gi];
                    nlinks = virtpi[Gd];

                    if(nrows && ncols && nlinks && W0[Gjk] != NULL)
                      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                              &(T2i.matrix[Gbd][bd][0]), nrows,
                              &(F3ints.matrix[Gac][ac][id]), nlinks, 1.0,
                              &(W0[Gjk][0][0]), ncols);


                }

                /* -E_bljk * t_acil */

                  for(Gl=0; Gl < nirreps; Gl++){

                    Gjk = Gbl = Gb ^ Gl;
                    Gi = Gl ^ Gac;

                /* Set up  T2i amplitudes     */
                    il = T2i.col_offset[Gac][Gi];

                /* Set up multiplication parameters */
                    nrows = E2ints.params->coltot[Gbl];
                    ncols = occpi[Gi];
                    nlinks = occpi[Gl];

                    if(nrows && ncols && nlinks && W0[Gjk] != NULL){

                    E2ints.matrix[Gbl] = global_dpd_->dpd_block_matrix(occpi[Gl], E2ints.params->coltot[Gbl]);
                    global_dpd_->buf4_mat_irrep_rd_block(&E2ints, Gbl, E2ints.row_offset[Gbl][B], occpi[Gl]);

                      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                              &(E2ints.matrix[Gbl][0][0]), nrows,
                              &(T2i.matrix[Gac][ac][il]), nlinks, 1.0,
                              &(W0[Gjk][0][0]), ncols);

                    global_dpd_->free_dpd_block(E2ints.matrix[Gbl], occpi[Gl], E2ints.params->coltot[Gbl]);
                        }
		}
              /* Sort W[jk][i] --> W[ji][k] */
         global_dpd_->sort_3d(W0, W1, nirreps, Gabc, Eints.params->rowtot, Eints.params->rowidx,
                              Eints.params->roworb, Eints.params->psym, Eints.params->qsym,
                              occ_off, occ_off, occpi, occ_off, Eints.params->rowidx, acb, 0);

               /* t_bdji * F_cakd */

                  for(Gd=0; Gd < nirreps; Gd++) {

                    Gji = Gbd = Gb ^ Gd;
                    Gk = Gd ^ Gca;

                   /* Set up T2i amplitudes*/
                  bd = T2i.row_offset[Gbd][B];        


                   /* Set up F3ints integrals*/
                    kd = F3ints.col_offset[Gca][Gk];

                    /* Set up multiplication parameters */
                    nrows = T2i.params->coltot[Gbd];
                    ncols = occpi[Gk];
                    nlinks = virtpi[Gd];

                    if(nrows && ncols && nlinks && W1[Gji] != NULL)
                      C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                              &(T2i.matrix[Gbd][bd][0]), nrows,
                              &(F3ints.matrix[Gca][ca][kd]), nlinks, 1.0,
                              &(W1[Gji][0][0]), ncols);
                  
                }
                    
                /* -E_blji * t_cakl */

                  for(Gl=0; Gl < nirreps; Gl++){
                    
                    Gji = Gbl = Gb ^ Gl;
                    Gk = Gl ^ Gca;

                /* Set up  T2i amplitudes    */
                    kl = T2i.col_offset[Gca][Gk];

                /* Set up multiplication parameters */
                    nrows = E2ints.params->coltot[Gbl];
                    ncols = occpi[Gk];
                    nlinks = occpi[Gl];

                    if(nrows && ncols && nlinks && W1[Gji] != NULL){

                    E2ints.matrix[Gbl] = global_dpd_->dpd_block_matrix(occpi[Gl], E2ints.params->coltot[Gbl]);
                    global_dpd_->buf4_mat_irrep_rd_block(&E2ints, Gbl, E2ints.row_offset[Gbl][B], occpi[Gl]);

                      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
                              &(E2ints.matrix[Gbl][0][0]), nrows,
                            &(T2i.matrix[Gca][ca][kl]), nlinks, 1.0,
                              &(W1[Gji][0][0]), ncols);

                    global_dpd_->free_dpd_block(E2ints.matrix[Gbl], occpi[Gl], E2ints.params->coltot[Gbl]);
                        }
		}

               /* Sort W[ji][k] --> W[ij][k] */
          global_dpd_->sort_3d(W1, W0, nirreps, Gabc, Eints.params->rowtot, Eints.params->rowidx,
                              Eints.params->roworb, Eints.params->psym, Eints.params->qsym,
                              occ_off, occ_off, occpi, occ_off, Eints.params->rowidx, bac, 0);

                  timer_off("N7 Terms(abc)");

                  /**** T3(c) complete ****/

                  /**** Build T3(d) for current a,b,c; Result stored in V ****/

                  timer_on("T3(d) Terms(abc)");

              for(Gij=0; Gij < nirreps; Gij++) {

                    Gk = Gij ^ Gabc;

                    for(ij=0; ij < Eints.params->rowtot[Gij]; ij++) {

                      I = Eints.params->roworb[Gij][ij][0];
                      Gi = Eints.params->psym[I];
                      i = I - occ_off[Gi];
                      J = Eints.params->roworb[Gij][ij][1];
                      Gj = Eints.params->qsym[J];
                      j = J - occ_off[Gj];

                      Gjk = Gj ^ Gk;
                      Gik = Gi ^ Gk;

                      for(k=0; k < occpi[Gk]; k++) {
                        K = occ_off[Gk] + k;

                        jk = Dints.params->rowidx[J][K];
                        ik = Dints.params->rowidx[I][K];

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

                          V[Gij][ij][k] += t_ia * D_jkbc + f_ia * t_jkbc;

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

			 V[Gij][ij][k] += t_jb * D_ikac + f_jb * t_ikac;

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

                          V[Gij][ij][k] += t_kc * D_ijab + f_kc * t_ijab;
                        }

                      } /* k */
                    } /* ij */
                  } /* Gij */
                  
                 timer_off("T3(d) Terms(abc)");



                 /**** Compute (T) Energy as a Test ****/
                  /*for(Gij=0; Gij < nirreps; Gij++) {
                    Gk = Gij ^ Gabc; Gji = Gij;
                    for(ij=0; ij < Eints.params->rowtot[Gij]; ij++) {
                      I = Eints.params->roworb[Gij][ij][0];
                      Gi = Eints.params->psym[I];
                      i = I - occ_off[Gi];
                      J = Eints.params->roworb[Gij][ij][1];
                      Gj = Eints.params->qsym[J];
                      j = J - occ_off[Gj];

                      Gik = Gki = Gi ^ Gk;  Gjk = Gkj = Gj ^ Gk;

                      ji = Dints.params->rowidx[J][I];

                      for(k=0; k < occpi[Gk]; k++) {
                        K = occ_off[Gk] + k;

                        ik = Eints.params->rowidx[I][K];
                        ki = Eints.params->rowidx[K][I];
                        jk = Eints.params->rowidx[J][K];
                        kj = Eints.params->rowidx[K][J];

                        denom = dabc;
                        if(fIJ.params->rowtot[Gi])
                          denom += fIJ.matrix[Gi][i][i];
                        if(fIJ.params->rowtot[Gj])
                          denom += fIJ.matrix[Gj][j][j];
                        if(fIJ.params->rowtot[Gk])
                          denom += fIJ.matrix[Gk][k][k];
			

                          value1 = W0[Gij][ij][k] + V[Gij][ij][k] - W0[Gkj][kj][i] - V[Gkj][kj][i];
                          value2 = 4 * W0[Gij][ij][k] + W0[Gjk][jk][i] + W0[Gki][ki][j];
                          ET2 += value1 * value2 / (3.0 * denom);

                      } // k 
                    } // ij 
                  } // Gij
		*/

                  /**** (T) Energy Contribution Complete ****/

 		/* Build E3 array */

		for(Gij=0; Gij < nirreps; Gij++) {
                    Gk = Gij ^ Gabc; Gji = Gij;
                    for(ij=0; ij < Eints.params->rowtot[Gij]; ij++) {
                      I = Eints.params->roworb[Gij][ij][0];
                      Gi = Eints.params->psym[I];
                      i = I - occ_off[Gi];
                      J = Eints.params->roworb[Gij][ij][1];
                      Gj = Eints.params->qsym[J];
                      j = J - occ_off[Gj];

                      Gik = Gki = Gi ^ Gk;  Gjk = Gkj = Gj ^ Gk;

                      ji = Dints.params->rowidx[J][I];

                      for(k=0; k < occpi[Gk]; k++) {
                        K = occ_off[Gk] + k;

                        ik = Eints.params->rowidx[I][K];
                        ki = Eints.params->rowidx[K][I];
                        jk = Eints.params->rowidx[J][K];
                        kj = Eints.params->rowidx[K][J];

                        denom = dabc;
                        if(fIJ.params->rowtot[Gi])
                          denom += fIJ.matrix[Gi][i][i];
                        if(fIJ.params->rowtot[Gj])
                          denom += fIJ.matrix[Gj][j][j];
                        if(fIJ.params->rowtot[Gk])
                          denom += fIJ.matrix[Gk][k][k];

                        /* TJ Lee Expression */
                        M[Gij][ij][k] = 2 * (8 * W0[Gij][ij][k] + 2 * W0[Gjk][jk][i] + 2 * W0[Gki][ki][j]
                                             - 4 * W0[Gkj][kj][i] - 4 * W0[Gik][ik][j] - 4 * W0[Gji][ji][k]
                                             + 4 * V[Gij][ij][k] + V[Gjk][jk][i] + V[Gki][ki][j]
                                             - 2 * V[Gkj][kj][i] - 2 * V[Gik][ik][j] - 2 * V[Gji][ji][k])/denom;


                        E[Gij][ij][k] =  (8 * W0[Gij][ij][k] + 2 * W0[Gjk][jk][i] + 2 * W0[Gki][ki][j]
                                             - 4 * W0[Gkj][kj][i] - 4 * W0[Gik][ik][j] - 4 * W0[Gji][ji][k]
                                             + 2 * (4 * V[Gij][ij][k] + V[Gjk][jk][i] + V[Gki][ki][j]
                                             - 2 * V[Gkj][kj][i] - 2 * V[Gik][ik][j] - 2 * V[Gji][ji][k]))/(denom * denom);

                      }
                    }
                  }



                 for(Gij=0;Gij < nirreps; Gij++){
		    Gk = Gij ^ Gabc;
		   for(ij=0; ij< Eints.params->rowtot[Gij]; ij++){
		    I = Eints.params->roworb[Gij][ij][0];
		    Gi = Eints.params->psym[I];	
		    i = I - occ_off[Gi];
		     for (k=0; k < occpi[Gk]; k++){
		     	DIJ.matrix[Gi][i][i] -= 0.5 * W0[Gij][ij][k] * E[Gij][ij][k];
		    } /* Gij */
		} /* if Gi == Gj .....*/
	     } /* ij */ 


                  /**** Compute S2 needed for lambda equations ****/

		   /*  S_cbkj --> sum_ial -M_[il][k](abc) <aj|il>*/

                  for(Gil=0; Gil < nirreps; Gil++) {
                       Gk = Gil ^ Gabc;
		       Gj = Ga ^ Gil ;
		       Gaj = Gil;
		       Gkj = Gcb;	
				
                      nrows = occpi[Gk];
                      ncols = occpi[Gj];
                      nlinks = E2ints.params->coltot[Gil];
                      aj = E2ints.row_offset[Gaj][A];
		      kj = S2i.col_offset[Gcb][Gk];

                    if(nrows && ncols && nlinks){
		         E2ints.matrix[Gil] = global_dpd_->dpd_block_matrix(occpi[Gj], E2ints.params->coltot[Gil]);
                         global_dpd_->buf4_mat_irrep_rd_block(&E2ints, Gil, aj, occpi[Gj]);

                      C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, M[Gil][0], nrows,
                              E2ints.matrix[Gil][0], nlinks, 1.0, &(S2i.matrix[Gkj][cb][kj]), ncols);

                       global_dpd_->free_dpd_block(E2ints.matrix[Gil], occpi[Gj], E2ints.params->coltot[Gil]);

                    }
                  } /* Gil */
	           
                  timer_on("malloc");
                  for(Gij=0; Gij < nirreps; Gij++) {
                    Gk = Gij ^ Gabc;
			global_dpd_->free_dpd_block(W0[Gij],Eints.params->rowtot[Gij],occpi[Gk]);
			global_dpd_->free_dpd_block(W1[Gij],Eints.params->rowtot[Gij],occpi[Gk]);
			global_dpd_->free_dpd_block(V[Gij],Eints.params->rowtot[Gij],occpi[Gk]);
			global_dpd_->free_dpd_block(M[Gij],Eints.params->rowtot[Gij],occpi[Gk]);
          	       }

                timer_off("malloc");

                } /* c */
              } /* b */
            } /* a */

          } /* Gc */
        } /* Gb */
      } /* Ga */


      //outfile->Printf( "\tE(T)(ijk) = %20.14f\n", ET1);
      //outfile->Printf( "\tE(T)(abc) = %20.14f\n", ET2);

      free(W0); free(W1); free(V); free(M); free(E);
      free(W2); free(Mi); free(RW);


      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&S2i, h);
	global_dpd_->buf4_mat_irrep_close(&S2i, h);
	global_dpd_->buf4_mat_irrep_close(&T2, h);
	global_dpd_->buf4_mat_irrep_close(&Eints, h);
	global_dpd_->buf4_mat_irrep_close(&E2ints, h);
	global_dpd_->buf4_mat_irrep_close(&F3ints, h);
	global_dpd_->buf4_mat_irrep_close(&Dints, h);
      }
      global_dpd_->buf4_close(&S2i);
      global_dpd_->buf4_close(&T2);
      global_dpd_->buf4_close(&Eints);
      global_dpd_->buf4_close(&E2ints);
      global_dpd_->buf4_close(&Dints);
      global_dpd_->buf4_close(&F3ints);

      global_dpd_->file2_mat_wrt(&DIJ);
      global_dpd_->file2_mat_close(&DIJ);
      global_dpd_->file2_close(&DIJ);

      global_dpd_->file2_mat_close(&T1);
      global_dpd_->file2_close(&T1);

      global_dpd_->file2_mat_close(&fIJ);
      global_dpd_->file2_mat_close(&fAB);
      global_dpd_->file2_mat_close(&fIA);
      global_dpd_->file2_close(&fIJ);
      global_dpd_->file2_close(&fAB);
      global_dpd_->file2_close(&fIA);


        // symmetrizing wrt interchange of i,j and a,b    
        

      global_dpd_->buf4_init(&S2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "SIjAb(T)");
      global_dpd_->buf4_init(&S2i, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "SIjAb(T)");
      global_dpd_->buf4_sort_axpy(&S2i, PSIF_CC_MISC, qpsr, 0, 5, "SIjAb(T)", 1);
      global_dpd_->buf4_close(&S2);
      global_dpd_->buf4_close(&S2i);

      global_dpd_->buf4_init(&S2, PSIF_CC_TMP, 0, 5, 0, 5, 0, 0, "SAbIj(T)");
      global_dpd_->buf4_init(&S2i, PSIF_CC_TMP, 0, 5, 0, 5, 0, 0, "SAbIj(T)");
      global_dpd_->buf4_sort_axpy(&S2i, PSIF_CC_TMP, qpsr, 5, 0, "SAbIj(T)", 1);
      global_dpd_->buf4_close(&S2);
      global_dpd_->buf4_close(&S2i);
      
      global_dpd_->buf4_init(&S2, PSIF_CC_TMP, 0, 0, 5, 0, 5, 0, "SJkBc(T)");
      global_dpd_->buf4_init(&S2i, PSIF_CC_TMP, 0, 0, 5, 0, 5, 0, "SJkBc(T)");
      global_dpd_->buf4_sort_axpy(&S2i, PSIF_CC_TMP, qpsr, 0, 5, "SJkBc(T)", 1);
      global_dpd_->buf4_close(&S2);
      global_dpd_->buf4_close(&S2i); 
   
      // Adjusting to PSI4's definiton of density matrices

      global_dpd_->buf4_init(&S2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "SIjAb(T)");
      global_dpd_->buf4_init(&S2i, PSIF_CC_TMP, 0, 5, 0, 5, 0, 0, "SAbIj(T)");
      global_dpd_->buf4_sort_axpy(&S2i, PSIF_CC_MISC, rspq, 0, 5, "SIjAb(T)", 1);
      global_dpd_->buf4_close(&S2i);
      global_dpd_->buf4_close(&S2);

       
      global_dpd_->buf4_init(&S2, PSIF_CC_MISC, 0, 0, 5, 0, 5, 0, "SIjAb(T)");
      global_dpd_->buf4_scmcopy(&S2, PSIF_CC_TMP, "1/3 SIjAb + 1/6 SIjBa", 1.0/3.0);	
      global_dpd_->buf4_sort_axpy(&S2, PSIF_CC_TMP, pqsr, 0, 5, "1/3 SIjAb + 1/6 SIjBa", 1.0/6.0);
      global_dpd_->buf4_close(&S2);


      global_dpd_->buf4_init(&S2, PSIF_CC_TMP, 0, 0, 5, 0, 5, 0, "1/3 SIjAb + 1/6 SIjBa");
      global_dpd_->buf4_copy(&S2, PSIF_CC_MISC, "SIjAb(T)");
      global_dpd_->buf4_close(&S2);
     

      global_dpd_->file2_init(&DIJ, PSIF_CC_OEI, 0, 0, 0, "DIJ(T)");
      global_dpd_->file2_scm(&DIJ, 0.5);
      global_dpd_->file2_close(&DIJ);

      global_dpd_->file2_init(&DAB, PSIF_CC_OEI, 0, 1, 1, "DAB(T)");
      global_dpd_->file2_scm(&DAB, 0.5);
      global_dpd_->file2_close(&DAB);

      global_dpd_->buf4_init(&si, PSIF_CC_TMP, 0, 10, 5, 10, 5, 0, "SIbAd(T)");
      global_dpd_->buf4_scmcopy(&si, PSIF_CC_TMP, "GIcAb(T)", 1.0/3.0);	
      global_dpd_->buf4_sort_axpy(&si, PSIF_CC_TMP, pqsr, 10, 5, "GIcAb(T)", 1.0/6.0);
      global_dpd_->buf4_close(&si);

      global_dpd_->buf4_init(&si, PSIF_CC_TMP, 0, 10, 5, 10, 5, 0, "GIcAb(T)");
      global_dpd_->buf4_sort(&si, PSIF_CC_FINTS, qpsr, 11, 5, "GCiAb(T)");
      global_dpd_->buf4_close(&si);

      
      global_dpd_->buf4_init(&kappa, PSIF_CC_TMP, 0, 0, 10, 0, 10, 0, "SIjLa(T)");
      global_dpd_->buf4_scmcopy(&kappa, PSIF_CC_EINTS, "GIjKa(T)", 1.0/3.0);
      global_dpd_->buf4_sort_axpy(&kappa, PSIF_CC_EINTS, qprs, 0, 10, "GIjKa(T)", 1.0/6.0);
      global_dpd_->buf4_close(&kappa);

      global_dpd_->buf4_init(&S2, PSIF_CC_TMP, 0, 0, 5, 0, 5, 0, "SJkBc(T)");
      global_dpd_->buf4_scmcopy(&S2, PSIF_CC_TMP, "2/3 SJkBc + 1/3 SJkCb", 2.0/3.0);
      global_dpd_->buf4_sort_axpy(&S2, PSIF_CC_TMP, pqsr, 0, 5, "2/3 SJkBc + 1/3 SJkCb", 1.0/3.0);
      global_dpd_->buf4_close(&S2);


      global_dpd_->buf4_init(&S2, PSIF_CC_TMP, 0, 0, 5, 0, 5, 0, "2/3 SJkBc + 1/3 SJkCb");
      global_dpd_->buf4_copy(&S2, PSIF_CC_DINTS, "GIjAb(T)");
      global_dpd_->buf4_close(&S2);
      
    }


  }}
