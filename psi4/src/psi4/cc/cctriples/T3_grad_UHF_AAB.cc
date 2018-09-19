/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCTRIPLES
    \brief Enter brief description of file here
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cctriples {

extern void T3_UHF_AAB(double ***W, double ***V, int disc, int nirreps,
    int I, int Gi, int J, int Gj, int K, int Gk, dpdbuf4 *T2AA, dpdbuf4 *T2AB,
    dpdbuf4 *T2BA, dpdbuf4 *FAA, dpdbuf4 *FAB, dpdbuf4 *FBA, dpdbuf4 *EAA,
    dpdbuf4 *EAB, dpdbuf4 *EBA, dpdfile2 *T1A, dpdfile2 *T1B, dpdbuf4 *DAA,
    dpdbuf4 *DAB, dpdfile2 *fIA, dpdfile2 *fia, dpdfile2 *fIJ, dpdfile2 *fij,
    dpdfile2 *fAB, dpdfile2 *fab, int *aoccpi, int *aocc_off, int *boccpi,
    int *bocc_off, int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off,
    double omega);

extern void T3_UHF_AAB_abc(double ***W, double ***V, int disc, int nirreps,
    int I, int Gi, int J, int Gj, int K, int Gk, dpdbuf4 *T2AA, dpdbuf4 *T2AB,
    dpdbuf4 *T2BA, dpdbuf4 *FAA, dpdbuf4 *FAB, dpdbuf4 *FBA, dpdbuf4 *EAA,
    dpdbuf4 *EAB, dpdbuf4 *EBA, dpdfile2 *T1A, dpdfile2 *T1B, dpdbuf4 *DAA,
    dpdbuf4 *DAB, dpdfile2 *fIA, dpdfile2 *fia, dpdfile2 *fIJ, dpdfile2 *fij,
    dpdfile2 *fAB, dpdfile2 *fab, int *aoccpi, int *aocc_off, int *boccpi,
    int *bocc_off, int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off,
    double omega);

    double T3_grad_UHF_AAB(void)
    {
      int cnt;
      int h, nirreps;
      int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
      int Gji, Gij, Gjk, Gkj, Gik, Gki, Gijk;
      int Gab, Gbc, Gac, Gcb, Gcd;
      int Gid, Gjd, Gkd;
      int Gil, Gjl, Gkl, Gli, Glk;
      int I, J, K, L, A, B, C, D;
      int i, j, k, l, a, b, c, d;
      int ij, ji, ik, ki, jk, kj;
      int ab, ba, ac, ca, bc, cb;
      int dc, ad, bd, da;
      int lc, la, lb;
      int id, jd, kd;
      int il, jl, kl, li, lk;
      int *aoccpi, *avirtpi, *aocc_off, *avir_off;
      int *boccpi, *bvirtpi, *bocc_off, *bvir_off;
      double value_c, value_d, dijk, denom, ET;
      int nrows, ncols, nlinks;
      dpdbuf4 T2AB, T2AA, T2BA;
      dpdbuf4 FAAints, FABints, FBAints;
      dpdbuf4 EAAints, EABints, EBAints;
      dpdbuf4 DAAints, DABints;
      dpdfile2 T1A, T1B, fIJ, fij, fAB, fab, fIA, fia;
      dpdfile2 S1A, S1B, DAB, Dab, DIJ, Dij;
      dpdbuf4 S2AA, S2AB, GIJAB, GIjAb, GIJKA, GIjKa, GiJkA, GIDAB, GIdAb, GiDaB;
      double ***WABc, ***VABc;
      double ***XABc, ***Y1, ***Y2;
      double **Z;

      nirreps = moinfo.nirreps;
      aoccpi = moinfo.aoccpi;
      avirtpi = moinfo.avirtpi;
      aocc_off = moinfo.aocc_off;
      avir_off = moinfo.avir_off;
      boccpi = moinfo.boccpi;
      bvirtpi = moinfo.bvirtpi;
      bocc_off = moinfo.bocc_off;
      bvir_off = moinfo.bvir_off;

      double ***WIJk = (double ***) malloc(nirreps * sizeof(double **));
      double ***VIJk = (double ***) malloc(nirreps * sizeof(double **));

      global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
      global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
      global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
      global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
      global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
      global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");

      global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
      global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");

      global_dpd_->buf4_init(&T2AA, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
      global_dpd_->buf4_init(&T2AB, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
      global_dpd_->buf4_init(&T2BA, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

      global_dpd_->buf4_init(&FAAints, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
      global_dpd_->buf4_init(&FABints, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
      global_dpd_->buf4_init(&FBAints, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");

      global_dpd_->buf4_init(&EAAints, PSIF_CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
      global_dpd_->buf4_init(&EABints, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
      global_dpd_->buf4_init(&EBAints, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");

      global_dpd_->buf4_init(&DAAints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
      global_dpd_->buf4_init(&DABints, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");

      global_dpd_->file2_init(&S1A, PSIF_CC_OEI, 0, 0, 1, "SIA");
      global_dpd_->file2_mat_init(&S1A);
      global_dpd_->file2_mat_rd(&S1A);
      global_dpd_->file2_init(&S1B, PSIF_CC_OEI, 0, 2, 3, "Sia");
      global_dpd_->file2_mat_init(&S1B);
      global_dpd_->file2_mat_rd(&S1B);

      global_dpd_->file2_init(&DAB, PSIF_CC_OEI, 0, 1, 1, "DAB");
      global_dpd_->file2_mat_init(&DAB);
      global_dpd_->file2_mat_rd(&DAB);
      global_dpd_->file2_init(&Dab, PSIF_CC_OEI, 0, 3, 3, "Dab");
      global_dpd_->file2_mat_init(&Dab);
      global_dpd_->file2_mat_rd(&Dab);

      global_dpd_->buf4_init(&S2AA, PSIF_CC_MISC, 0, 0, 5, 2, 7, 0, "SIJAB");
      global_dpd_->buf4_init(&S2AB, PSIF_CC_MISC, 0, 22, 28, 22, 28, 0, "SIjAb");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&S2AA, h);
	global_dpd_->buf4_mat_irrep_rd(&S2AA, h);
	global_dpd_->buf4_mat_irrep_init(&S2AB, h);
      }

      global_dpd_->buf4_init(&GIJAB, PSIF_CC_GAMMA, 0, 0, 5, 2, 7, 0, "GIJAB");
      global_dpd_->buf4_init(&GIjAb, PSIF_CC_GAMMA, 0, 22, 28, 22, 28, 0, "GIjAb");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&GIJAB, h);
	global_dpd_->buf4_mat_irrep_rd(&GIJAB, h);
	global_dpd_->buf4_mat_irrep_init(&GIjAb, h);
      }

      global_dpd_->buf4_init(&GIJKA, PSIF_CC_GAMMA, 0, 0, 20, 2, 20, 0, "GIJKA");
      global_dpd_->buf4_init(&GIjKa, PSIF_CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
      global_dpd_->buf4_init(&GiJkA, PSIF_CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&GIJKA, h);
	global_dpd_->buf4_mat_irrep_rd(&GIJKA, h);
	global_dpd_->buf4_mat_irrep_init(&GIjKa, h);
	global_dpd_->buf4_mat_irrep_init(&GiJkA, h);
      }

      global_dpd_->buf4_init(&GIDAB, PSIF_CC_GAMMA, 0, 20, 5, 20, 7, 0, "GIDAB");
      global_dpd_->buf4_init(&GIdAb, PSIF_CC_GAMMA, 0, 24, 28, 24, 28, 0, "GIdAb");
      global_dpd_->buf4_init(&GiDaB, PSIF_CC_GAMMA, 0, 27, 29, 27, 29, 0, "GiDaB");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&GIDAB, h);
	global_dpd_->buf4_mat_irrep_rd(&GIDAB, h);
	global_dpd_->buf4_mat_irrep_init(&GIdAb, h);
	global_dpd_->buf4_mat_irrep_init(&GiDaB, h);
      }

      ET = 0.0;

      WABc = (double ***) malloc(nirreps * sizeof(double **));
      VABc = (double ***) malloc(nirreps * sizeof(double **));
      XABc = (double ***) malloc(nirreps * sizeof(double **));
      Y1 = (double ***) malloc(nirreps * sizeof(double **));
      Y2 = (double ***) malloc(nirreps * sizeof(double **));

      for(Gi=0; Gi < nirreps; Gi++) {
	for(Gj=0; Gj < nirreps; Gj++) {
	  for(Gk=0; Gk < nirreps; Gk++) {

	    Gij = Gji = Gi ^ Gj;
	    Gjk = Gkj = Gj ^ Gk;
	    Gik = Gki = Gi ^ Gk;

	    Gijk = Gi ^ Gj ^ Gk;

	    for(Gab=0; Gab < nirreps; Gab++) {
	      Gc = Gab ^ Gijk;

	      WABc[Gab] = global_dpd_->dpd_block_matrix(FAAints.params->coltot[Gab], bvirtpi[Gc]);
	      VABc[Gab] = global_dpd_->dpd_block_matrix(FAAints.params->coltot[Gab], bvirtpi[Gc]);
	      XABc[Gab] = global_dpd_->dpd_block_matrix(FAAints.params->coltot[Gab], bvirtpi[Gc]);
	    }

	    for(Ga=0; Ga < nirreps; Ga++) {
	      Gbc = Ga ^ Gijk;
	      Y1[Ga] = global_dpd_->dpd_block_matrix(avirtpi[Ga], FBAints.params->coltot[Gbc]); /* alpha-beta-alpha */
	      Y2[Ga] = global_dpd_->dpd_block_matrix(avirtpi[Ga], FABints.params->coltot[Gbc]); /* alpha-alpha-beta */
	    }

	    for(i=0; i < aoccpi[Gi]; i++) {
	      I = aocc_off[Gi] + i;
	      for(j=0; j < aoccpi[Gj]; j++) {
		J = aocc_off[Gj] + j;
		for(k=0; k < boccpi[Gk]; k++) {
		  K = bocc_off[Gk] + k;

		  T3_UHF_AAB(WABc, VABc, 1, nirreps, I, Gi, J, Gj, K, Gk, &T2AA, &T2AB, &T2BA,
			     &FAAints, &FABints, &FBAints, &EAAints, &EABints, &EBAints,
			     &T1A, &T1B, &DAAints, &DABints, &fIA, &fia, &fIJ, &fij, &fAB, &fab,
			     aoccpi, aocc_off, boccpi, bocc_off, avirtpi, avir_off, bvirtpi, bvir_off, 0.0);

		  global_dpd_->file2_mat_init(&fIJ);
		  global_dpd_->file2_mat_init(&fij);
		  global_dpd_->file2_mat_init(&fAB);
		  global_dpd_->file2_mat_init(&fab);
		  global_dpd_->file2_mat_init(&fIA);
		  global_dpd_->file2_mat_init(&fia);
		  global_dpd_->file2_mat_rd(&fIJ);
		  global_dpd_->file2_mat_rd(&fij);
		  global_dpd_->file2_mat_rd(&fAB);
		  global_dpd_->file2_mat_rd(&fab);
		  global_dpd_->file2_mat_rd(&fIA);
		  global_dpd_->file2_mat_rd(&fia);
		  global_dpd_->file2_mat_init(&T1A);
		  global_dpd_->file2_mat_rd(&T1A);
		  global_dpd_->file2_mat_init(&T1B);
		  global_dpd_->file2_mat_rd(&T1B);
		  for(h=0; h < nirreps; h++) {
		    global_dpd_->buf4_mat_irrep_init(&T2AA, h);
		    global_dpd_->buf4_mat_irrep_rd(&T2AA, h);

		    global_dpd_->buf4_mat_irrep_init(&T2AB, h);
		    global_dpd_->buf4_mat_irrep_rd(&T2AB, h);

		    global_dpd_->buf4_mat_irrep_init(&T2BA, h);
		    global_dpd_->buf4_mat_irrep_rd(&T2BA, h);

		    global_dpd_->buf4_mat_irrep_init(&EAAints, h);
		    global_dpd_->buf4_mat_irrep_rd(&EAAints, h);

		    global_dpd_->buf4_mat_irrep_init(&EABints, h);
		    global_dpd_->buf4_mat_irrep_rd(&EABints, h);

		    global_dpd_->buf4_mat_irrep_init(&EBAints, h);
		    global_dpd_->buf4_mat_irrep_rd(&EBAints, h);

		    global_dpd_->buf4_mat_irrep_init(&DAAints, h);
		    global_dpd_->buf4_mat_irrep_rd(&DAAints, h);

		    global_dpd_->buf4_mat_irrep_init(&DABints, h);
		    global_dpd_->buf4_mat_irrep_rd(&DABints, h);
		  }

		  ij = EAAints.params->rowidx[I][J];
		  ji = EAAints.params->rowidx[J][I];
		  jk = EABints.params->rowidx[J][K];
		  kj = EBAints.params->rowidx[K][J];
		  ik = EABints.params->rowidx[I][K];
		  ki = EBAints.params->rowidx[K][I];

		  dijk = 0.0;
		  if(fIJ.params->rowtot[Gi]) dijk += fIJ.matrix[Gi][i][i];
		  if(fIJ.params->rowtot[Gj]) dijk += fIJ.matrix[Gj][j][j];
		  if(fij.params->rowtot[Gk]) dijk += fij.matrix[Gk][k][k];


		  /**** Apply denominators and compute AAB part of (T) as a test ****/
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;

		    for(ab=0; ab < FAAints.params->coltot[Gab]; ab++) {
		      A = FAAints.params->colorb[Gab][ab][0];
		      Ga = FAAints.params->rsym[A];
		      a = A - avir_off[Ga];
		      B = FAAints.params->colorb[Gab][ab][1];
		      Gb = FAAints.params->ssym[B];
		      b = B - avir_off[Gb];

		      for(c=0; c < bvirtpi[Gc]; c++) {
			C = bvir_off[Gc] + c;

			denom = dijk;
			if(fAB.params->rowtot[Ga]) denom -= fAB.matrix[Ga][a][a];
			if(fAB.params->rowtot[Gb]) denom -= fAB.matrix[Gb][b][b];
			if(fab.params->rowtot[Gc]) denom -= fab.matrix[Gc][c][c];

			ET += WABc[Gab][ab][c] * (WABc[Gab][ab][c] + VABc[Gab][ab][c]) * denom;

		      } /* c */
		    } /* ab */
		  } /* Gab */

		  /**** T3 --> S1 ****/

		  /* S_IA = <Jk|Bc> t(c)_IJkABc */
		  /* S_kc = 1/4 <IJ||AB> t(c)_IJkABc */
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FAAints.params->coltot[Gab]; ab++) {
		      A = FAAints.params->colorb[Gab][ab][0];
		      Ga = FAAints.params->rsym[A];
		      a = A - avir_off[Ga];
		      B = FAAints.params->colorb[Gab][ab][1];
		      Gb = FAAints.params->ssym[B];
		      b = B - avir_off[Gb];
		      Gbc = Gb ^ Gc;
		      Gac = Ga ^ Gc;
		      for(c=0; c < bvirtpi[Gc]; c++) {
			C = bvir_off[Gc] + c;
			bc = DABints.params->colidx[B][C];

			if(Gi==Ga && S1A.params->rowtot[Gi] && S1A.params->coltot[Gi])
			  S1A.matrix[Gi][i][a] += WABc[Gab][ab][c] * DABints.matrix[Gjk][jk][bc];

			if(Gk==Gc && S1B.params->rowtot[Gk] && S1B.params->coltot[Gk])
			  S1B.matrix[Gk][k][c] += 0.25 * WABc[Gab][ab][c] * DAAints.matrix[Gij][ij][ab];

		      } /* c */
		    } /* ab */
		  } /* Gab */

		  /**** T3 --> S1 Complete ****/

		  /**** T3 --> S2 ****/

		  /*** Build X_IJkABc = 2 W_IJkABc + V_IJkABc ***/
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FAAints.params->coltot[Gab]; ab++) {
		      for(c=0; c < bvirtpi[Gc]; c++) {
			XABc[Gab][ab][c] = 2 * WABc[Gab][ab][c] + VABc[Gab][ab][c];
		      }
		    }
		  }
		  /*** X_IJkABc Complete ***/

		  /*** Sort X(AB,c) to Y(A,cB) ***/
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FAAints.params->coltot[Gab]; ab++) {
		      A = FAAints.params->colorb[Gab][ab][0];
		      B = FAAints.params->colorb[Gab][ab][1];
		      Ga = FAAints.params->rsym[A];
		      a = A - avir_off[Ga];
		      for(c=0; c < bvirtpi[Gc]; c++) {
			C = bvir_off[Gc] + c;
			cb = FBAints.params->colidx[C][B];
			Y1[Ga][a][cb] = XABc[Gab][ab][c];
		      }
		    }
		  }
		  /*** S_JIDA <-- +t_IJkABc W_kDcB ***/
		  /*** S_JIAD <-- -t_IJkABc W_kDcB ***/
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Ga = Gd ^ Gji;
		    Gkd = Gcb = Gk ^ Gd;
		    kd = FBAints.row_offset[Gkd][K];
		    nrows = avirtpi[Gd];
		    ncols = avirtpi[Ga];
		    nlinks = FBAints.params->coltot[Gkd];
		    if(nrows && ncols && nlinks) {
		      FBAints.matrix[Gkd] = global_dpd_->dpd_block_matrix(nrows, nlinks);
		      global_dpd_->buf4_mat_irrep_rd_block(&FBAints, Gkd, kd, nrows);
		      Z = block_matrix(nrows, ncols);

		      C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, FBAints.matrix[Gkd][0], nlinks,
			      Y1[Ga][0], nlinks, 0.0, Z[0], ncols);

		      for(d=0; d < avirtpi[Gd]; d++) {
			D = avir_off[Gd] + d;
			for(a=0; a < avirtpi[Ga]; a++) {
			  A = avir_off[Ga] + a;
			  ad = S2AA.params->colidx[A][D];
			  da = S2AA.params->colidx[D][A];
			  S2AA.matrix[Gji][ji][da] += Z[d][a];
			  S2AA.matrix[Gji][ji][ad] -= Z[d][a];
			}
		      }

		      global_dpd_->free_dpd_block(FBAints.matrix[Gkd], nrows, nlinks);
		      free_block(Z);
		    } /* nrows && ncols && nlinks */
		  } /* Gd */

		    /*** S_LIAB <-- +t_IJkABc <Jk|Lc> ***/
		    /*** S_ILAB <-- -t_IJkABc <Jk|Lc> ***/
		  for(Gl=0; Gl < nirreps; Gl++) {
		    Gli = Gab = Gl ^ Gi;
		    Gc = Gab ^ Gijk;

		    nrows = aoccpi[Gl];
		    ncols = FAAints.params->coltot[Gab];
		    nlinks = bvirtpi[Gc];

		    if(nrows && ncols && nlinks) {
		      lc = EABints.col_offset[Gjk][Gl];
		      Z = block_matrix(nrows, ncols);
		      C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, &(EABints.matrix[Gjk][jk][lc]), nlinks,
			      XABc[Gab][0], nlinks, 0.0, Z[0], ncols);
		      for(l=0; l < nrows; l++) {
			L = aocc_off[Gl] + l;
			li = S2AA.params->rowidx[L][I];
			il = S2AA.params->rowidx[I][L];
			for(ab=0; ab < ncols; ab++) {
			  S2AA.matrix[Gli][li][ab] += Z[l][ab];
			  S2AA.matrix[Gli][il][ab] -= Z[l][ab];
			}
		      }
		      free_block(Z);
		    } /* nrows && ncols && nlinks */
		  } /* Gl */

		    /* S_JkDc <-- 1/2 <ID||AB> X_IJkABc */
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Gid = Gab = Gi ^ Gd;
		    Gc = Gab ^ Gijk;

		    nrows = avirtpi[Gd];
		    ncols = bvirtpi[Gc];
		    nlinks = FAAints.params->coltot[Gid];
		    if(nrows && ncols && nlinks) {
		      id = FAAints.row_offset[Gid][I];
		      FAAints.matrix[Gid] = global_dpd_->dpd_block_matrix(nrows, nlinks);
		      global_dpd_->buf4_mat_irrep_rd_block(&FAAints, Gid, id, nrows);
		      Z = block_matrix(nrows, ncols);
		      C_DGEMM('n', 'n', nrows, ncols, nlinks, 0.5, FAAints.matrix[Gid][0], nlinks,
			      XABc[Gab][0], ncols, 0.0, Z[0], ncols);

		      for(d=0; d < nrows; d++) {
			D = avir_off[Gd] + d;
			for(c=0; c < ncols; c++) {
			  C = bvir_off[Gc] + c;
			  dc = S2AB.params->colidx[D][C];
			  S2AB.matrix[Gjk][jk][dc] += Z[d][c];
			}
		      }

		      global_dpd_->free_dpd_block(FAAints.matrix[Gid], nrows, nlinks);
		      free_block(Z);
		    } /* nrows && ncols && nlinks */
		  } /* Gd */

		    /* S_JkBd <-- X_IJkABc <Id|Ac> */
		    /* sort X(AB,c) to Y2(B,Ac) */
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FAAints.params->coltot[Gab]; ab++) {
		      A = FAAints.params->colorb[Gab][ab][0];
		      B = FAAints.params->colorb[Gab][ab][1];
		      Gb = FAAints.params->ssym[B];
		      b = B - avir_off[Gb];
		      for(c=0; c < bvirtpi[Gc]; c++) {
			C = bvir_off[Gc] + c;
			ac = FABints.params->colidx[A][C];
			Y2[Gb][b][ac] = XABc[Gab][ab][c];
		      }
		    }
		  }

		  for(Gd=0; Gd < nirreps; Gd++) {
		    Gid = Gac = Gi ^ Gd;
		    Gb = Gac ^ Gijk;

		    nrows = avirtpi[Gb];
		    ncols = bvirtpi[Gd];
		    nlinks = FABints.params->coltot[Gid];

		    if(nrows && ncols && nlinks) {
		      id = FABints.row_offset[Gid][I];
		      FABints.matrix[Gid] = global_dpd_->dpd_block_matrix(ncols, nlinks);
		      global_dpd_->buf4_mat_irrep_rd_block(&FABints, Gid, id, ncols);
		      Z = block_matrix(nrows, ncols);
		      C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, Y2[Gb][0], nlinks,
			      FABints.matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

		      for(b=0; b < nrows; b++) {
			B = avir_off[Gb] + b;
			for(d=0; d < ncols; d++) {
			  D = bvir_off[Gd] + d;
			  bd = S2AB.params->colidx[B][D];
			  S2AB.matrix[Gjk][jk][bd] += Z[b][d];
			}
		      }

		      global_dpd_->free_dpd_block(FABints.matrix[Gid], ncols, nlinks);
		      free_block(Z);

		    } /* nrows && ncols && nlinks */
		  } /* Gd */

		  /* S_LkBc <-- 1/2 <IJ||LA> X_IJkABc */
		  /* sort X(AB,c) to Y2(A,Bc) */
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FAAints.params->coltot[Gab]; ab++) {
		      A = FAAints.params->colorb[Gab][ab][0];
		      B = FAAints.params->colorb[Gab][ab][1];
		      Ga = FAAints.params->rsym[A];
		      a = A - avir_off[Ga];
		      for(c=0; c < bvirtpi[Gc]; c++) {
			C = bvir_off[Gc] + c;
			bc = S2AB.params->colidx[B][C];
			Y2[Ga][a][bc] = XABc[Gab][ab][c];
		      } /* c */
		    } /* ab */
		  } /* Gab */

		  for(Gl=0; Gl < nirreps; Gl++) {
		    Glk = Gbc = Gl ^ Gk;
		    Ga = Gbc ^ Gijk;

		    nrows = aoccpi[Gl];
		    ncols = S2AB.params->coltot[Glk];
		    nlinks = avirtpi[Ga];
		    if(nrows && ncols && nlinks) {
		      la = EAAints.col_offset[Gij][Gl];
		      Z = global_dpd_->dpd_block_matrix(nrows, ncols);
		      C_DGEMM('n', 'n', nrows, ncols, nlinks, 0.5, &(EAAints.matrix[Gij][ij][la]), nlinks,
			      Y2[Ga][0], ncols, 0.0, Z[0], ncols);
		      for(l=0; l < nrows; l++) {
			L = aocc_off[Gl] + l;
			lk = S2AB.params->rowidx[L][K];
			for(bc=0; bc < ncols; bc++) {
			  S2AB.matrix[Glk][lk][bc] += Z[l][bc];
			}
		      }

		      global_dpd_->free_dpd_block(Z, nrows, ncols);
		    } /* nrows && ncols && nlinks */
		  } /* Gl */

		  /* S_IlBc <-- <kJ|lA> X_IJkABc */
		  for(Gl=0; Gl < nirreps; Gl++) {
		    Gil = Gbc = Gi ^ Gl;
		    Ga = Gbc ^ Gijk;

		    nrows = boccpi[Gl];
		    ncols = S2AB.params->coltot[Gil];
		    nlinks = avirtpi[Ga];
		    if(nrows && ncols && nlinks) {
		      la = EBAints.col_offset[Gjk][Gl];
		      Z = global_dpd_->dpd_block_matrix(nrows, ncols);
		      C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, &(EBAints.matrix[Gjk][kj][la]), nlinks,
			      Y2[Ga][0], ncols, 0.0, Z[0], ncols);
		      for(l=0; l < nrows; l++) {
			L = bocc_off[Gl] + l;
			il = S2AB.params->rowidx[I][L];
			for(bc=0; bc < ncols; bc++) {
			  S2AB.matrix[Gil][il][bc] += Z[l][bc];
			}
		      }
		      global_dpd_->free_dpd_block(Z, nrows, ncols);
		    } /* nrows && ncols && nlinks */
		  } /* Gl */

		  /**** T3 --> S2 Complete ****/

		  /**** T3 --> DAB ****/
		  for(Ga=0; Ga < nirreps; Ga++) {
		    Gb = Ga;
		    Gcd = Ga ^ Gijk;
		    for(Gc=0; Gc < nirreps; Gc++) {
		      Gd = Gc ^ Gcd;
		      Gac = Gbc = Ga ^ Gc;
		      for(a=0; a < avirtpi[Ga]; a++) {
			A = avir_off[Ga] + a;
			for(b=0; b < avirtpi[Gb]; b++) {
			  B = avir_off[Gb] + b;
			  for(c=0; c < avirtpi[Gc]; c++) {
			    C = avir_off[Gc] + c;
			    ac = FAAints.params->colidx[A][C];
			    bc = FAAints.params->colidx[B][C];
			    for(d=0; d < bvirtpi[Gd]; d++) {
			      DAB.matrix[Ga][b][a] += 0.5 * WABc[Gac][ac][d] * (WABc[Gbc][bc][d] + VABc[Gbc][bc][d]);
			    } /* d */
			  } /* c */
			} /* b */
		      } /* a */
		    } /* Gc */
		  } /* Ga */

		  /**** T3 --> DAB complete ****/

		  /**** T3 --> Dab ****/

		  for(Gc=0; Gc < nirreps; Gc++) {
		    Gd = Gc;
		    Gab = Gc ^ Gijk;
		    for(ab=0; ab < FAAints.params->coltot[Gab]; ab++) {
		      for(c=0; c < bvirtpi[Gc]; c++) {
			for(d=0; d < bvirtpi[Gd]; d++) {
			  Dab.matrix[Gc][d][c] += 0.25 * WABc[Gab][ab][c] * (WABc[Gab][ab][d] + VABc[Gab][ab][d]);
			}
		      }
		    } /* ab */
		  } /* Gc */

		  /**** T3 --> Dab complete ****/

		  /* T3 --> GIJAB ****/

		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    if(Gk == Gc) {
		      for(ab=0; ab < FAAints.params->coltot[Gab]; ab++) {
			for(c=0; c < bvirtpi[Gc]; c++) {
			  C = bvir_off[Gc] + c;
			  if(T1B.params->rowtot[Gk] && T1B.params->coltot[Gk])
			    GIJAB.matrix[Gij][ij][ab] += WABc[Gab][ab][c] * T1B.matrix[Gk][k][c];
			}
		      }
		    }
		  } /* Gab */

		  /**** T3 --> GIJAB complete ****/

		  /**** T3 --> GIjAb ****/
		  /* Sort W(AB,c) --> Y2(A,Bc) */
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FAAints.params->coltot[Gab]; ab++) {
		      A = FAAints.params->colorb[Gab][ab][0];
		      B = FAAints.params->colorb[Gab][ab][1];
		      Ga = FAAints.params->rsym[A];
		      a = A - avir_off[Ga];
		      for(c=0; c < bvirtpi[Gc]; c++) {
			C = bvir_off[Gc] + c;
			bc = S2AB.params->colidx[B][C];
			Y2[Ga][a][bc] = WABc[Gab][ab][c];
		      } /* c */
		    } /* ab */
		  } /* Gab */

		  Ga = Gi; Gbc = Ga ^ Gijk;
		  if(T1A.params->rowtot[Gi] && T1A.params->coltot[Gi]) {
		    for(a=0; a < avirtpi[Ga]; a++) {
		      for(bc=0; bc < GIjAb.params->coltot[Gbc]; bc++) {
			GIjAb.matrix[Gjk][jk][bc] += Y2[Ga][a][bc] * T1A.matrix[Gi][i][a];
		      }
		    }
		  }

		  /**** T3 --> GiJaB complete ****/

		  /**** T3 --> GIJKA ****/
		  /* Sort W(AB,c) --> Y1(A,cB) */
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FAAints.params->coltot[Gab]; ab++) {
		      A = FAAints.params->colorb[Gab][ab][0];
		      B = FAAints.params->colorb[Gab][ab][1];
		      Ga = FAAints.params->rsym[A];
		      a = A - avir_off[Ga];
		      for(c=0; c < bvirtpi[Gc]; c++) {
			C = bvir_off[Gc] + c;
			cb = T2BA.params->colidx[C][B];
			Y1[Ga][a][cb] = 2 * WABc[Gab][ab][c] + VABc[Gab][ab][c];
		      } /* c */
		    } /* ab */
		  } /* Gab */

		  /* G_IJLA <-- t_kLcB Y_IJkABc */
		  for(Gl=0; Gl < nirreps; Gl++) {
		    Ga = Gl ^ Gij;
		    Gkl = Gcb = Gk ^ Gl;

		    nrows = aoccpi[Gl];
		    ncols = avirtpi[Ga];
		    nlinks = T2BA.params->coltot[Gcb];
		    if(nrows && ncols && nlinks) {
		      kl = T2BA.row_offset[Gkl][K];
		      la = GIJKA.col_offset[Gij][Gl];
		      C_DGEMM('n','t', nrows, ncols, nlinks, 1.0, T2BA.matrix[Gkl][kl], nlinks,
			      Y1[Ga][0], nlinks, 1.0, &(GIJKA.matrix[Gij][ij][la]), ncols);
		    }
		  } /* Gl */

		  /**** T3 --> GIJKA complete ****/

		  /**** T3 --> GIjKa ****/
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FAAints.params->coltot[Gab]; ab++) {
		      for(c=0; c < bvirtpi[Gc]; c++) {
			XABc[Gab][ab][c] = 2 * WABc[Gab][ab][c] + VABc[Gab][ab][c];
		      } /* c */
		    } /* ab */
		  } /* Gab */

		  /* GIkLc <-- 1/2 t_JLAB X_IJkABc */
		  for(Gl=0; Gl < nirreps; Gl++) {
		    Gc = Gl ^ Gik;
		    Gab = Gjl = Gj ^ Gl;
		    nrows = aoccpi[Gl];
		    ncols = bvirtpi[Gc];
		    nlinks = T2AA.params->coltot[Gjl];
		    if(nrows && ncols && nlinks) {
		      jl = T2AA.row_offset[Gjl][J];
		      lc = GIjKa.col_offset[Gik][Gl];
		      C_DGEMM('n','n', nrows, ncols, nlinks, 0.5, T2AA.matrix[Gjl][jl], nlinks,
			      XABc[Gab][0], ncols, 1.0, &(GIjKa.matrix[Gik][ik][lc]), ncols);
		    }
		  } /* Gl */

		  /**** T3 --> GIjKa complete ****/

		  /**** T3 --> GiJkA ****/
		  /* Sort W(AB,c) --> Y2(A,Bc) */
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FAAints.params->coltot[Gab]; ab++) {
		      A = FAAints.params->colorb[Gab][ab][0];
		      B = FAAints.params->colorb[Gab][ab][1];
		      Ga = FAAints.params->rsym[A];
		      a = A - avir_off[Ga];
		      for(c=0; c < bvirtpi[Gc]; c++) {
			C = bvir_off[Gc] + c;
			bc = S2AB.params->colidx[B][C];
			Y2[Ga][a][bc] = 2 * WABc[Gab][ab][c] + VABc[Gab][ab][c];
		      } /* c */
		    } /* ab */
		  } /* Gab */

		  /* G_kIlA <-- -t_JlBc X_IJkABc **/
		  for(Gl=0; Gl < nirreps; Gl++) {
		    Ga = Gki ^ Gl;
		    Gjl = Gbc = Gj ^ Gl;
		    nrows = boccpi[Gl];
		    ncols = avirtpi[Ga];
		    nlinks = T2AB.params->coltot[Gbc];
		    if(nrows && ncols && nlinks) {
		      jl = T2AB.row_offset[Gjl][J];
		      la = GiJkA.col_offset[Gki][Gl];
		      C_DGEMM('n','t', nrows, ncols, nlinks, -1.0, T2AB.matrix[Gjl][jl], nlinks,
			      Y2[Ga][0], nlinks, 1.0, &(GiJkA.matrix[Gki][ki][la]), ncols);
		    }
		  } /* Gl */

		  /**** T3 --> GiJkA complete ****/

		  /* GIDAB <-- -t_JkDc X_IJkABc */
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Gab = Gid = Gi ^ Gd;
		    Gc = Gjk ^ Gd;

		    nrows = avirtpi[Gd];
		    ncols = GIDAB.params->coltot[Gid];
		    nlinks = bvirtpi[Gc];
		    if(nrows && ncols && nlinks) {
		      id = GIDAB.row_offset[Gid][I];
		      dc = T2AB.col_offset[Gjk][Gd];
		      C_DGEMM('n','t',nrows, ncols, nlinks, -1.0, &(T2AB.matrix[Gjk][jk][dc]), nlinks,
			      XABc[Gab][0], nlinks, 1.0, GIDAB.matrix[Gid][id], ncols);
		    }
		  }
		  /*** T3 --> GIDAB complete ***/

		  /* GIdBc <-- t_JkAd t_IJkABc */
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Ga = Gd ^ Gjk;
		    Gid = Gi ^ Gd;

		    nrows = bvirtpi[Gd];
		    ncols = GIdAb.params->coltot[Gid];
		    nlinks = avirtpi[Ga];
		    if(nrows && ncols && nlinks) {
		      ad = T2AB.col_offset[Gjk][Ga];
		      id = GIdAb.row_offset[Gid][I];
		      C_DGEMM('t','n',nrows, ncols, nlinks, -1.0, &(T2AB.matrix[Gjk][jk][ad]), nrows,
			      Y2[Ga][0], ncols, 1.0, GIdAb.matrix[Gid][id], ncols);
		    }
		  }
		  /*** T3 --> GIdAb complete ***/

		  /* GkDcA <-- -1/2 t_IJAD t_IJkABc */
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Ga = Gd ^ Gij;
		    Gkd = Gk ^ Gd;

		    nrows = avirtpi[Gd];
		    ncols = GiDaB.params->coltot[Gkd];
		    nlinks = avirtpi[Ga];
		    if(nrows && ncols && nlinks) {
		      ad = T2AA.col_offset[Gij][Ga];
		      kd = GiDaB.row_offset[Gkd][K];
		      C_DGEMM('t','n', nrows, ncols, nlinks, 0.5, &(T2AA.matrix[Gij][ij][ad]), nrows,
			      Y1[Ga][0], ncols, 1.0, GiDaB.matrix[Gkd][kd], ncols);
		    }
		  }
		  /*** T3 --> GiDaB complete ***/

		  for(h=0; h < nirreps; h++) {
		    global_dpd_->buf4_mat_irrep_close(&T2AA, h);
		    global_dpd_->buf4_mat_irrep_close(&T2AB, h);
		    global_dpd_->buf4_mat_irrep_close(&T2BA, h);
		    global_dpd_->buf4_mat_irrep_close(&EAAints, h);
		    global_dpd_->buf4_mat_irrep_close(&EABints, h);
		    global_dpd_->buf4_mat_irrep_close(&EBAints, h);
		    global_dpd_->buf4_mat_irrep_close(&DAAints, h);
		    global_dpd_->buf4_mat_irrep_close(&DABints, h);
		  }
		  global_dpd_->file2_mat_close(&T1A);
		  global_dpd_->file2_mat_close(&T1B);
		  global_dpd_->file2_mat_close(&fIJ);
		  global_dpd_->file2_mat_close(&fij);
		  global_dpd_->file2_mat_close(&fAB);
		  global_dpd_->file2_mat_close(&fab);
		  global_dpd_->file2_mat_close(&fIA);
		  global_dpd_->file2_mat_close(&fia);

		} /* k */
	      } /* j */
	    } /* i */

	    for(Gab=0; Gab < nirreps; Gab++) {
	      Gc = Gab ^ Gijk;
	      global_dpd_->free_dpd_block(WABc[Gab], FAAints.params->coltot[Gab], bvirtpi[Gc]);
	      global_dpd_->free_dpd_block(VABc[Gab], FAAints.params->coltot[Gab], bvirtpi[Gc]);
	      global_dpd_->free_dpd_block(XABc[Gab], FAAints.params->coltot[Gab], bvirtpi[Gc]);
	    }
	    for(Ga=0; Ga < nirreps; Ga++) {
	      Gbc = Ga ^ Gijk;
	      global_dpd_->free_dpd_block(Y1[Ga], avirtpi[Ga], FBAints.params->coltot[Gbc]);
	      global_dpd_->free_dpd_block(Y2[Ga], avirtpi[Ga], FABints.params->coltot[Gbc]);
	    }

	  } /* Gk */
	} /* Gj */
      } /* Gi */

      ET *= 0.25;

      free(WABc);
      free(VABc);
      free(XABc);
      free(Y1); free(Y2);

      global_dpd_->file2_mat_wrt(&DAB);
      global_dpd_->file2_mat_close(&DAB);
      global_dpd_->file2_close(&DAB);
      global_dpd_->file2_mat_wrt(&Dab);
      global_dpd_->file2_mat_close(&Dab);
      global_dpd_->file2_close(&Dab);

      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&S2AA, h);
	global_dpd_->buf4_mat_irrep_close(&S2AA, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&S2AB, h);
	global_dpd_->buf4_mat_irrep_close(&S2AB, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&GIJAB, h);
	global_dpd_->buf4_mat_irrep_close(&GIJAB, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&GIjAb, h);
	global_dpd_->buf4_mat_irrep_close(&GIjAb, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&GIJKA, h);
	global_dpd_->buf4_mat_irrep_close(&GIJKA, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&GIjKa, h);
	global_dpd_->buf4_mat_irrep_close(&GIjKa, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&GiJkA, h);
	global_dpd_->buf4_mat_irrep_close(&GiJkA, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&GIDAB, h);
	global_dpd_->buf4_mat_irrep_close(&GIDAB, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&GIdAb, h);
	global_dpd_->buf4_mat_irrep_close(&GIdAb, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&GiDaB, h);
	global_dpd_->buf4_mat_irrep_close(&GiDaB, h);
      }
      global_dpd_->buf4_close(&S2AA);
      global_dpd_->buf4_close(&S2AB);
      global_dpd_->buf4_close(&GIJAB);
      global_dpd_->buf4_close(&GIjAb);
      global_dpd_->buf4_close(&GIJKA);
      global_dpd_->buf4_close(&GIjKa);
      global_dpd_->buf4_close(&GiJkA);
      global_dpd_->buf4_close(&GIDAB);
      global_dpd_->buf4_close(&GIdAb);
      global_dpd_->buf4_close(&GiDaB);

      global_dpd_->file2_mat_wrt(&S1A);
      global_dpd_->file2_mat_close(&S1A);
      global_dpd_->file2_close(&S1A);
      global_dpd_->file2_mat_wrt(&S1B);
      global_dpd_->file2_mat_close(&S1B);
      global_dpd_->file2_close(&S1B);

      global_dpd_->buf4_close(&T2AA);
      global_dpd_->buf4_close(&T2AB);
      global_dpd_->buf4_close(&T2BA);
      global_dpd_->buf4_close(&FAAints);
      global_dpd_->buf4_close(&FABints);
      global_dpd_->buf4_close(&FBAints);
      global_dpd_->buf4_close(&EAAints);
      global_dpd_->buf4_close(&EABints);
      global_dpd_->buf4_close(&EBAints);
      global_dpd_->buf4_close(&DAAints);
      global_dpd_->buf4_close(&DABints);

      global_dpd_->file2_close(&T1A);
      global_dpd_->file2_close(&T1B);
      global_dpd_->file2_close(&fIJ);
      global_dpd_->file2_close(&fij);
      global_dpd_->file2_close(&fAB);
      global_dpd_->file2_close(&fab);
      global_dpd_->file2_close(&fIA);
      global_dpd_->file2_close(&fia);

      /*** T3 --> DIJ and Dij ***/

      global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
      global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
      global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
      global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
      global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
      global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
      global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
      global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");

      global_dpd_->buf4_init(&T2AA, PSIF_CC_TAMPS, 0, 5, 0, 7, 2, 0, "tABIJ");
      global_dpd_->buf4_init(&T2AB, PSIF_CC_TAMPS, 0, 28, 22, 28, 22, 0, "tAbIj");
      global_dpd_->buf4_init(&T2BA, PSIF_CC_TAMPS, 0, 29, 23, 29, 23, 0, "taBiJ");
      global_dpd_->buf4_init(&FAAints, PSIF_CC_FINTS, 0, 5, 20, 7, 20, 0, "F <BC||IA>");
      global_dpd_->buf4_init(&FBAints, PSIF_CC_FINTS, 0, 29, 27, 29, 27, 0, "F <bC|iA>");
      global_dpd_->buf4_init(&FABints, PSIF_CC_FINTS, 0, 28, 24, 28, 24, 0, "F <Bc|Ia>");
      global_dpd_->buf4_init(&EAAints, PSIF_CC_EINTS, 0, 21, 0, 21, 2, 0, "E <AK||IJ> (AK, I>J)");
      global_dpd_->buf4_init(&EABints, PSIF_CC_EINTS, 0, 25, 22, 25, 22, 0, "E <aK|Ij>");
      global_dpd_->buf4_init(&EBAints, PSIF_CC_EINTS, 0, 26, 23, 26, 23, 0, "E <Ak|iJ>");
      global_dpd_->buf4_init(&DAAints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
      global_dpd_->buf4_init(&DABints, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");

      global_dpd_->file2_init(&DIJ, PSIF_CC_OEI, 0, 0, 0, "DIJ");
      global_dpd_->file2_mat_init(&DIJ);
      global_dpd_->file2_mat_rd(&DIJ);
      global_dpd_->file2_init(&Dij, PSIF_CC_OEI, 0, 2, 2, "Dij");
      global_dpd_->file2_mat_init(&Dij);
      global_dpd_->file2_mat_rd(&Dij);

      int Gabc;
      for (Ga=0; Ga < nirreps; ++Ga) {
        for (a=0; a<avirtpi[Ga]; ++a) {
          A = avir_off[Ga] + a;
          for (Gb=0; Gb < nirreps; ++Gb) {
            for (b=0; b<avirtpi[Gb]; ++b) {
              B = avir_off[Gb] + b;
              for (Gc=0; Gc < nirreps; ++Gc) {
                for (c=0; c < bvirtpi[Gc]; ++c) {
                  C = bvir_off[Gc] + c;
                  Gabc = Ga ^ Gb ^ Gc;
                  //Allocate the memory for connected and disconnected triples
                  for (Gij=0; Gij < nirreps; ++Gij) {
                    Gk = Gij ^ Gabc;
                    WIJk[Gij] = global_dpd_->dpd_block_matrix(T2AA.params->coltot[Gij], boccpi[Gk]);
                    VIJk[Gij] = global_dpd_->dpd_block_matrix(T2AA.params->coltot[Gij], boccpi[Gk]);
                  }
                  T3_UHF_AAB_abc(WIJk, VIJk, 1, nirreps, A, Ga, B, Gb, C, Gc,
                      &T2AA, &T2AB, &T2BA, &FAAints, &FABints, &FBAints, &EAAints,
                      &EABints, &EBAints, &T1A, &T1B, &DAAints, &DABints, &fIA,
                      &fia, &fIJ, &fij, &fAB, &fab, aoccpi, aocc_off, boccpi,
                      bocc_off, avirtpi, avir_off, bvirtpi, bvir_off, 0.0);

                  for(Gi=0; Gi < nirreps; Gi++) {
                    Gj = Gi;
                    Gkl = Gi ^ Gabc;
                    for(Gk=0; Gk < nirreps; Gk++) {
                      Gl = Gk ^ Gkl;
                      Gik = Gjk = Gi ^ Gk;
                      for(i=0; i < aoccpi[Gi]; i++) {
                        I = aocc_off[Gi] + i;
                        for(j=0; j < aoccpi[Gj]; j++) {
                          J = aocc_off[Gj] + j;
                          for(k=0; k < aoccpi[Gk]; k++) {
                            K = aocc_off[Gk] + k;
                            ik = T2AA.params->colidx[I][K];
                            jk = T2AA.params->colidx[J][K];
                            for(l=0; l < boccpi[Gl]; l++) {
                              DIJ.matrix[Gi][i][j] -= 0.5 * WIJk[Gik][ik][l] * (WIJk[Gjk][jk][l] + VIJk[Gjk][jk][l]);
                            } /* l */
                          } /* k */
                        } /* j */
                      } /* i */
                    } /* Gk */
                  } /* Gi */

                  for(Gi=0; Gi < nirreps; Gi++) {
                    Gj = Gi;
                    Gkl = Gi ^ Gabc;
                    for(kl=0; kl < T2AA.params->coltot[Gkl]; kl++) {
                      for(i=0; i < boccpi[Gi]; i++) {
                        for(j=0; j < boccpi[Gj]; j++) {
                          Dij.matrix[Gi][i][j] -= 0.25 * WIJk[Gkl][kl][i] * (WIJk[Gkl][kl][j] + VIJk[Gkl][kl][j]);
                        } /* j */
                      } /* i */
                    } /* kl */
                  } /* Gi */

                  //Deallocate the memory for connected and disconnected triples
                  for (Gij=0; Gij < nirreps; ++Gij) {
                    Gk = Gij ^ Gabc;
                    global_dpd_->free_dpd_block(WIJk[Gij], T2AA.params->coltot[Gij], boccpi[Gk]);
                    global_dpd_->free_dpd_block(VIJk[Gij], T2AA.params->coltot[Gij], boccpi[Gk]);
                  }
                }
              }
            }
          }
        }
      }

      global_dpd_->file2_mat_wrt(&DIJ);
      global_dpd_->file2_mat_close(&DIJ);
      global_dpd_->file2_close(&DIJ);
      global_dpd_->file2_mat_wrt(&Dij);
      global_dpd_->file2_mat_close(&Dij);
      global_dpd_->file2_close(&Dij);

      global_dpd_->file2_close(&fij);
      global_dpd_->file2_close(&fIJ);
      global_dpd_->file2_close(&fab);
      global_dpd_->file2_close(&fAB);

      global_dpd_->buf4_close(&T2AA);
      global_dpd_->buf4_close(&T2AB);
      global_dpd_->buf4_close(&T2BA);
      global_dpd_->buf4_close(&EAAints);
      global_dpd_->buf4_close(&EABints);
      global_dpd_->buf4_close(&EBAints);
      global_dpd_->buf4_close(&FAAints);
      global_dpd_->buf4_close(&FABints);
      global_dpd_->buf4_close(&FBAints);
      global_dpd_->buf4_close(&DAAints);
      global_dpd_->buf4_close(&DABints);

      return ET;
    }


  }} // namespace psi::CCTRIPLES
