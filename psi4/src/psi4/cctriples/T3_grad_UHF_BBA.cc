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
    \brief Enter brief description of file here
*/
 #include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
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

    double T3_grad_UHF_BBA(void)
    {
      int cnt;
      int h, nirreps;
      int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
      int Gji, Gij, Gjk, Gkj, Gik, Gki, Gijk;
      int Gab, Gbc, Gac, Gca, Gba, Gcb, Gcd;
      int Gid, Gjd, Gkd;
      int Gil, Gjl, Gkl, Gli, Glk;
      int I, J, K, L, A, B, C, D;
      int i, j, k, l, a, b, c, d;
      int ij, ji, ik, ki, jk, kj;
      int ab, ba, ac, ca, bc, cb;
      int cd, bd, ad, db, dc, da;
      int lc, lb, la;
      int id, jd, kd;
      int il, jl, kl, li, lk;
      int *aoccpi, *avirtpi, *aocc_off, *avir_off;
      int *boccpi, *bvirtpi, *bocc_off, *bvir_off;
      double value_c, value_d, dijk, denom, ET;
      int nrows, ncols, nlinks;
      dpdbuf4 T2AB, T2BB, T2BA;
      dpdbuf4 FBBints, FABints, FBAints;
      dpdbuf4 EBBints, EABints, EBAints;
      dpdbuf4 DBBints, DBAints;
      dpdfile2 T1A, T1B, fIJ, fij, fAB, fab, fIA, fia;
      dpdfile2 S1A, S1B, DAB, Dab, DIJ, Dij;
      dpdbuf4 S2BB, S2BA, Gijab, GiJaB, Gijka, GIjKa, GiJkA, Gidab, GiDaB, GIdAb;
      double ***WabC, ***VabC;
      double ***XabC, ***Y1, ***Y2;
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

      double ***WijK = (double ***) malloc(nirreps * sizeof(double **));
      double ***VijK = (double ***) malloc(nirreps * sizeof(double **));

      global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
      global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
      global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
      global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
      global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
      global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");

      global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
      global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");

      global_dpd_->buf4_init(&T2BB, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
      global_dpd_->buf4_init(&T2AB, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
      global_dpd_->buf4_init(&T2BA, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

      global_dpd_->buf4_init(&FBBints, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
      global_dpd_->buf4_init(&FABints, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
      global_dpd_->buf4_init(&FBAints, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");

      global_dpd_->buf4_init(&EBBints, PSIF_CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
      global_dpd_->buf4_init(&EABints, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
      global_dpd_->buf4_init(&EBAints, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");

      global_dpd_->buf4_init(&DBBints, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
      global_dpd_->buf4_init(&DBAints, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");

      global_dpd_->file2_init(&S1A, PSIF_CC_OEI, 0, 0, 1, "SIA");
      global_dpd_->file2_mat_init(&S1A);
      global_dpd_->file2_mat_rd(&S1A);
      global_dpd_->file2_init(&S1B, PSIF_CC_OEI, 0, 2, 3, "Sia");
      global_dpd_->file2_mat_init(&S1B);
      global_dpd_->file2_mat_rd(&S1B);

      global_dpd_->buf4_init(&S2BB, PSIF_CC_MISC, 0, 10, 15, 12, 17, 0, "Sijab");
      global_dpd_->buf4_init(&S2BA, PSIF_CC_MISC, 0, 23, 29, 23, 29, 0, "SiJaB");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&S2BB, h);
  	global_dpd_->buf4_mat_irrep_rd(&S2BB, h);
	global_dpd_->buf4_mat_irrep_init(&S2BA, h);
      }

      global_dpd_->file2_init(&DAB, PSIF_CC_OEI, 0, 1, 1, "DAB");
      global_dpd_->file2_mat_init(&DAB);
      global_dpd_->file2_mat_rd(&DAB);
      global_dpd_->file2_init(&Dab, PSIF_CC_OEI, 0, 3, 3, "Dab");
      global_dpd_->file2_mat_init(&Dab);
      global_dpd_->file2_mat_rd(&Dab);

      global_dpd_->buf4_init(&Gijab, PSIF_CC_GAMMA, 0, 10, 15, 12, 17, 0, "Gijab");
      global_dpd_->buf4_init(&GiJaB, PSIF_CC_GAMMA, 0, 23, 29, 23, 29, 0, "GiJaB");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&Gijab, h);
	global_dpd_->buf4_mat_irrep_rd(&Gijab, h);
	global_dpd_->buf4_mat_irrep_init(&GiJaB, h);
      }

      global_dpd_->buf4_init(&Gijka, PSIF_CC_GAMMA, 0, 10, 30, 12, 30, 0, "Gijka");
      global_dpd_->buf4_init(&GIjKa, PSIF_CC_GAMMA, 0, 22, 24, 22, 24, 0, "GIjKa");
      global_dpd_->buf4_init(&GiJkA, PSIF_CC_GAMMA, 0, 23, 27, 23, 27, 0, "GiJkA");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&Gijka, h);
	global_dpd_->buf4_mat_irrep_rd(&Gijka, h);
	global_dpd_->buf4_mat_irrep_init(&GIjKa, h);
	global_dpd_->buf4_mat_irrep_rd(&GIjKa, h);
	global_dpd_->buf4_mat_irrep_init(&GiJkA, h);
	global_dpd_->buf4_mat_irrep_rd(&GiJkA, h);
      }

      global_dpd_->buf4_init(&Gidab, PSIF_CC_GAMMA, 0, 30, 15, 30, 17, 0, "Gidab");
      global_dpd_->buf4_init(&GiDaB, PSIF_CC_GAMMA, 0, 27, 29, 27, 29, 0, "GiDaB");
      global_dpd_->buf4_init(&GIdAb, PSIF_CC_GAMMA, 0, 24, 28, 24, 28, 0, "GIdAb");
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_init(&Gidab, h);
	global_dpd_->buf4_mat_irrep_rd(&Gidab, h);
	global_dpd_->buf4_mat_irrep_init(&GiDaB, h);
	global_dpd_->buf4_mat_irrep_rd(&GiDaB, h);
	global_dpd_->buf4_mat_irrep_init(&GIdAb, h);
	global_dpd_->buf4_mat_irrep_rd(&GIdAb, h);
      }

      ET = 0.0;

      WabC = (double ***) malloc(nirreps * sizeof(double **));
      VabC = (double ***) malloc(nirreps * sizeof(double **));
      XabC = (double ***) malloc(nirreps * sizeof(double **));
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

	      WabC[Gab] = global_dpd_->dpd_block_matrix(FBBints.params->coltot[Gab], avirtpi[Gc]);
	      VabC[Gab] = global_dpd_->dpd_block_matrix(FBBints.params->coltot[Gab], avirtpi[Gc]);
	      XabC[Gab] = global_dpd_->dpd_block_matrix(FBBints.params->coltot[Gab], avirtpi[Gc]);
	    }

	    for(Ga=0; Ga < nirreps; Ga++) {
	      Gbc = Ga ^ Gijk;
	      Y1[Ga] = global_dpd_->dpd_block_matrix(bvirtpi[Ga], FABints.params->coltot[Gbc]); /* beta-alpha-beta */
	      Y2[Ga] = global_dpd_->dpd_block_matrix(bvirtpi[Ga], FBAints.params->coltot[Gbc]); /* beta-beta-alpha */
	    }

	    for(i=0; i < boccpi[Gi]; i++) {
	      I = bocc_off[Gi] + i;
	      for(j=0; j < boccpi[Gj]; j++) {
		J = bocc_off[Gj] + j;
		for(k=0; k < aoccpi[Gk]; k++) {
		  K = aocc_off[Gk] + k;

		  T3_UHF_AAB(WabC, VabC, 1, nirreps, I, Gi, J, Gj, K, Gk, &T2BB, &T2BA, &T2AB,
			     &FBBints, &FBAints, &FABints, &EBBints, &EBAints, &EABints,
			     &T1B, &T1A, &DBBints, &DBAints, &fia, &fIA, &fij, &fIJ, &fab, &fAB,
			     boccpi, bocc_off, aoccpi, aocc_off, bvirtpi, bvir_off, avirtpi, avir_off, 0.0);

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
		    global_dpd_->buf4_mat_irrep_init(&T2BB, h);
		    global_dpd_->buf4_mat_irrep_rd(&T2BB, h);

		    global_dpd_->buf4_mat_irrep_init(&T2AB, h);
		    global_dpd_->buf4_mat_irrep_rd(&T2AB, h);

		    global_dpd_->buf4_mat_irrep_init(&T2BA, h);
		    global_dpd_->buf4_mat_irrep_rd(&T2BA, h);

		    global_dpd_->buf4_mat_irrep_init(&EBBints, h);
		    global_dpd_->buf4_mat_irrep_rd(&EBBints, h);

		    global_dpd_->buf4_mat_irrep_init(&EABints, h);
		    global_dpd_->buf4_mat_irrep_rd(&EABints, h);

		    global_dpd_->buf4_mat_irrep_init(&EBAints, h);
		    global_dpd_->buf4_mat_irrep_rd(&EBAints, h);

		    global_dpd_->buf4_mat_irrep_init(&DBBints, h);
		    global_dpd_->buf4_mat_irrep_rd(&DBBints, h);

		    global_dpd_->buf4_mat_irrep_init(&DBAints, h);
		    global_dpd_->buf4_mat_irrep_rd(&DBAints, h);
		  }

		  ij = EBBints.params->rowidx[I][J];
		  ji = EBBints.params->rowidx[J][I];
		  jk = EBAints.params->rowidx[J][K];
		  kj = EABints.params->rowidx[K][J];
		  ik = EBAints.params->rowidx[I][K];
		  ki = EABints.params->rowidx[K][I];

		  dijk = 0.0;
		  if(fij.params->rowtot[Gi]) dijk += fij.matrix[Gi][i][i];
		  if(fij.params->rowtot[Gj]) dijk += fij.matrix[Gj][j][j];
		  if(fIJ.params->rowtot[Gk]) dijk += fIJ.matrix[Gk][k][k];

		  /*** Compute BBA contribution to (T) as a test ***/
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;

		    for(ab=0; ab < FBBints.params->coltot[Gab]; ab++) {
		      A = FBBints.params->colorb[Gab][ab][0];
		      Ga = FBBints.params->rsym[A];
		      a = A - bvir_off[Ga];
		      B = FBBints.params->colorb[Gab][ab][1];
		      Gb = FBBints.params->ssym[B];
		      b = B - bvir_off[Gb];

		      for(c=0; c < avirtpi[Gc]; c++) {
			C = avir_off[Gc] + c;

			denom = dijk;
			if(fab.params->rowtot[Ga]) denom -= fab.matrix[Ga][a][a];
			if(fab.params->rowtot[Gb]) denom -= fab.matrix[Gb][b][b];
			if(fAB.params->rowtot[Gc]) denom -= fAB.matrix[Gc][c][c];

			ET += WabC[Gab][ab][c] * (WabC[Gab][ab][c] + VabC[Gab][ab][c]) * denom;

		      } /* c */
		    } /* ab */
		  } /* Gab */

		  /**** T3 --> S1 ****/

		  /* S_ia = <jK|bC> t(c)_ijKabC */
		  /* S_KC = 1/4 <ij||ab> t(c)_ijKabC */
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FBBints.params->coltot[Gab]; ab++) {
		      A = FBBints.params->colorb[Gab][ab][0];
		      Ga = FBBints.params->rsym[A];
		      a = A - bvir_off[Ga];
		      B = FBBints.params->colorb[Gab][ab][1];
		      Gb = FBBints.params->ssym[B];
		      b = B - bvir_off[Gb];
		      Gbc = Gb ^ Gc;
		      Gac = Ga ^ Gc;
		      for(c=0; c < avirtpi[Gc]; c++) {
			C = avir_off[Gc] + c;
			bc = DBAints.params->colidx[B][C];

			if(Gi==Ga && S1B.params->rowtot[Gi] && S1B.params->coltot[Gi])
			  S1B.matrix[Gi][i][a] += WabC[Gab][ab][c] * DBAints.matrix[Gjk][jk][bc];

			if(Gk==Gc && S1A.params->rowtot[Gk] && S1A.params->coltot[Gk])
			  S1A.matrix[Gk][k][c] += 0.25 * WabC[Gab][ab][c] * DBBints.matrix[Gij][ij][ab];

		      } /* c */
		    } /* ab */
		  } /* Gab */

		  /**** T3 --> S1 Complete ****/

		  /**** T3 --> S2 ****/

		  /*** Build X_ijKabC = 2 W_ijKabC + V_ijKabC ***/
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FBBints.params->coltot[Gab]; ab++) {
		      for(c=0; c < avirtpi[Gc]; c++) {
			XabC[Gab][ab][c] = 2 * WabC[Gab][ab][c] + VabC[Gab][ab][c];
		      }
		    }
		  }
		  /*** X_ijkabC Complete ***/

		  /*** Sort X(ab,C) to Y(a,Cb) ***/
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FBBints.params->coltot[Gab]; ab++) {
		      A = FBBints.params->colorb[Gab][ab][0];
		      B = FBBints.params->colorb[Gab][ab][1];
		      Ga = FBBints.params->rsym[A];
		      a = A - bvir_off[Ga];
		      for(c=0; c < avirtpi[Gc]; c++) {
			C = avir_off[Gc] + c;
			cb = FABints.params->colidx[C][B];
			Y1[Ga][a][cb] = XabC[Gab][ab][c];
		      }
		    }
		  }
		  /*** S_jida <-- +t_ijKabD W_KdCb ***/
		  /*** S_jiad <-- -t_ijKabC W_KdCb ***/
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Ga = Gd ^ Gji;
		    Gkd = Gcb = Gk ^ Gd;
		    kd = FABints.row_offset[Gkd][K];
		    nrows = bvirtpi[Gd];
		    ncols = bvirtpi[Ga];
		    nlinks = FABints.params->coltot[Gkd];
		    if(nrows && ncols && nlinks) {
		      FABints.matrix[Gkd] = global_dpd_->dpd_block_matrix(nrows, nlinks);
		      global_dpd_->buf4_mat_irrep_rd_block(&FABints, Gkd, kd, nrows);
		      Z = block_matrix(nrows, ncols);

		      C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, FABints.matrix[Gkd][0], nlinks,
			      Y1[Ga][0], nlinks, 0.0, Z[0], ncols);

		      for(d=0; d < bvirtpi[Gd]; d++) {
			D = bvir_off[Gd] + d;
			for(a=0; a < bvirtpi[Ga]; a++) {
			  A = bvir_off[Ga] + a;
			  ad = S2BB.params->colidx[A][D];
			  da = S2BB.params->colidx[D][A];
			  S2BB.matrix[Gji][ji][da] += Z[d][a];
			  S2BB.matrix[Gji][ji][ad] -= Z[d][a];
			}
		      }

		      global_dpd_->free_dpd_block(FABints.matrix[Gkd], nrows, nlinks);
		      free_block(Z);
		    } /* nrows && ncols && nlinks */
		  } /* Gd */

		    /*** S_liab <-- +t_ijKabC <jK|lC> ***/
		    /*** S_ilab <-- -t_ijKabC <jK|lC> ***/
		  for(Gl=0; Gl < nirreps; Gl++) {
		    Gli = Gab = Gl ^ Gi;
		    Gc = Gab ^ Gijk;

		    nrows = boccpi[Gl];
		    ncols = FBBints.params->coltot[Gab];
		    nlinks = avirtpi[Gc];

		    if(nrows && ncols && nlinks) {
		      lc = EBAints.col_offset[Gjk][Gl];
		      Z = block_matrix(nrows, ncols);
		      C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, &(EBAints.matrix[Gjk][jk][lc]), nlinks,
			      XabC[Gab][0], nlinks, 0.0, Z[0], ncols);
		      for(l=0; l < nrows; l++) {
			L = bocc_off[Gl] + l;
			li = S2BB.params->rowidx[L][I];
			il = S2BB.params->rowidx[I][L];
			for(ab=0; ab < ncols; ab++) {
			  S2BB.matrix[Gli][li][ab] += Z[l][ab];
			  S2BB.matrix[Gli][il][ab] -= Z[l][ab];
			}
		      }
		      free_block(Z);
		    } /* nrows && ncols && nlinks */
		  } /* Gl */

		    /* S_jKdC <-- 1/2 <id||ab> X_ijKabC */
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Gid = Gab = Gi ^ Gd;
		    Gc = Gab ^ Gijk;

		    nrows = bvirtpi[Gd];
		    ncols = avirtpi[Gc];
		    nlinks = FBBints.params->coltot[Gid];
		    if(nrows && ncols && nlinks) {
		      id = FBBints.row_offset[Gid][I];
		      FBBints.matrix[Gid] = global_dpd_->dpd_block_matrix(nrows, nlinks);
		      global_dpd_->buf4_mat_irrep_rd_block(&FBBints, Gid, id, nrows);
		      Z = block_matrix(nrows, ncols);
		      C_DGEMM('n', 'n', nrows, ncols, nlinks, 0.5, FBBints.matrix[Gid][0], nlinks,
			      XabC[Gab][0], ncols, 0.0, Z[0], ncols);

		      for(d=0; d < nrows; d++) {
			D = bvir_off[Gd] + d;
			for(c=0; c < ncols; c++) {
			  C = avir_off[Gc] + c;
			  dc = S2BA.params->colidx[D][C];
			  S2BA.matrix[Gjk][jk][dc] += Z[d][c];
			}
		      }

		      global_dpd_->free_dpd_block(FBBints.matrix[Gid], nrows, nlinks);
		      free_block(Z);
		    } /* nrows && ncols && nlinks */
		  } /* Gd */

		    /* S_jKbD <-- X_ijKabC <iD|aC> */
		    /* sort X(ab,C) to Y2(b,aC) */
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FBBints.params->coltot[Gab]; ab++) {
		      A = FBBints.params->colorb[Gab][ab][0];
		      B = FBBints.params->colorb[Gab][ab][1];
		      Gb = FBBints.params->ssym[B];
		      b = B - bvir_off[Gb];
		      for(c=0; c < avirtpi[Gc]; c++) {
			C = avir_off[Gc] + c;
			ac = FBAints.params->colidx[A][C];
			Y2[Gb][b][ac] = XabC[Gab][ab][c];
		      }
		    }
		  }

		  for(Gd=0; Gd < nirreps; Gd++) {
		    Gid = Gac = Gi ^ Gd;
		    Gb = Gac ^ Gijk;

		    nrows = bvirtpi[Gb];
		    ncols = avirtpi[Gd];
		    nlinks = FBAints.params->coltot[Gid];

		    if(nrows && ncols && nlinks) {
		      id = FBAints.row_offset[Gid][I];
		      FBAints.matrix[Gid] = global_dpd_->dpd_block_matrix(ncols, nlinks);
		      global_dpd_->buf4_mat_irrep_rd_block(&FBAints, Gid, id, ncols);
		      Z = block_matrix(nrows, ncols);
		      C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, Y2[Gb][0], nlinks,
			      FBAints.matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

		      for(b=0; b < nrows; b++) {
			B = bvir_off[Gb] + b;
			for(d=0; d < ncols; d++) {
			  D = avir_off[Gd] + d;
			  bd = S2BA.params->colidx[B][D];
			  S2BA.matrix[Gjk][jk][bd] += Z[b][d];
			}
		      }

		      global_dpd_->free_dpd_block(FBAints.matrix[Gid], ncols, nlinks);
		      free_block(Z);

		    } /* nrows && ncols && nlinks */
		  } /* Gd */

		  /* S_lKbC <-- 1/2 <ij||la> X_ijKabC */
		  /* sort X(ab,C) to Y2(a,bC) */
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FBBints.params->coltot[Gab]; ab++) {
		      A = FBBints.params->colorb[Gab][ab][0];
		      B = FBBints.params->colorb[Gab][ab][1];
		      Ga = FBBints.params->rsym[A];
		      a = A - bvir_off[Ga];
		      for(c=0; c < avirtpi[Gc]; c++) {
			C = avir_off[Gc] + c;
			bc = S2BA.params->colidx[B][C];
			Y2[Ga][a][bc] = XabC[Gab][ab][c];
		      } /* c */
		    } /* ab */
		  } /* Gab */

		  for(Gl=0; Gl < nirreps; Gl++) {
		    Glk = Gbc = Gl ^ Gk;
		    Ga = Gbc ^ Gijk;

		    nrows = boccpi[Gl];
		    ncols = S2BA.params->coltot[Glk];
		    nlinks = bvirtpi[Ga];
		    if(nrows && ncols && nlinks) {
		      la = EBBints.col_offset[Gij][Gl];
		      Z = global_dpd_->dpd_block_matrix(nrows, ncols);
		      C_DGEMM('n', 'n', nrows, ncols, nlinks, 0.5, &(EBBints.matrix[Gij][ij][la]), nlinks,
			      Y2[Ga][0], ncols, 0.0, Z[0], ncols);
		      for(l=0; l < nrows; l++) {
			L = bocc_off[Gl] + l;
			lk = S2BA.params->rowidx[L][K];
			for(bc=0; bc < ncols; bc++) {
			  S2BA.matrix[Glk][lk][bc] += Z[l][bc];
			}
		      }

		      global_dpd_->free_dpd_block(Z, nrows, ncols);
		    } /* nrows && ncols && nlinks */
		  } /* Gl */

		  /* S_iLbC <-- <Kj|La> X_ijKabC */
		  for(Gl=0; Gl < nirreps; Gl++) {
		    Gil = Gbc = Gi ^ Gl;
		    Ga = Gbc ^ Gijk;

		    nrows = aoccpi[Gl];
		    ncols = S2BA.params->coltot[Gil];
		    nlinks = bvirtpi[Ga];
		    if(nrows && ncols && nlinks) {
		      la = EABints.col_offset[Gjk][Gl];
		      Z = global_dpd_->dpd_block_matrix(nrows, ncols);
		      C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, &(EABints.matrix[Gjk][kj][la]), nlinks,
			      Y2[Ga][0], ncols, 0.0, Z[0], ncols);
		      for(l=0; l <nrows; l++) {
			L = aocc_off[Gl] + l;
			il = S2BA.params->rowidx[I][L];
			for(bc=0; bc < ncols; bc++) {
			  S2BA.matrix[Gil][il][bc] += Z[l][bc];
			}
		      }
		      global_dpd_->free_dpd_block(Z, nrows, ncols);
		    } /* nrows && ncols && nlinks */
		  } /* Gl */

		  /**** T3 --> S2 Complete ****/

		  /**** T3 --> DAB ****/

		  for(Gc=0; Gc < nirreps; Gc++) {
		    Gd = Gc;
		    Gab = Gc ^ Gijk;
		    for(ab=0; ab < FBBints.params->coltot[Gab]; ab++) {
		      for(c=0; c < avirtpi[Gc]; c++) {
			for(d=0; d < avirtpi[Gd]; d++) {
			  DAB.matrix[Gc][d][c] += 0.25 * WabC[Gab][ab][c] * (WabC[Gab][ab][d] + VabC[Gab][ab][d]);
			}
		      }
		    } /* ab */
		  } /* Gc */

		  /**** T3 --> DAB complete ****/

		  /**** T3 --> Dab ****/
		  for(Ga=0; Ga < nirreps; Ga++) {
		    Gb = Ga;
		    Gcd = Ga ^ Gijk;
		    for(Gc=0; Gc < nirreps; Gc++) {
		      Gd = Gc ^ Gcd;
		      Gac = Gbc = Ga ^ Gc;
		      for(a=0; a < bvirtpi[Ga]; a++) {
			A = bvir_off[Ga] + a;
			for(b=0; b < bvirtpi[Gb]; b++) {
			  B = bvir_off[Gb] + b;
			  for(c=0; c < bvirtpi[Gc]; c++) {
			    C = bvir_off[Gc] + c;
			    ac = FBBints.params->colidx[A][C];
			    bc = FBBints.params->colidx[B][C];
			    for(d=0; d < avirtpi[Gd]; d++) {
			      Dab.matrix[Ga][b][a] += 0.5 * WabC[Gac][ac][d] * (WabC[Gbc][bc][d] + VabC[Gbc][bc][d]);
			    } /* d */
			  } /* c */
			} /* b */
		      } /* a */
		    } /* Gc */
		  } /* Ga */

		  /**** T3 --> Dab complete ****/

		  /**** T3 --> Gijab ****/

		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    if(Gk == Gc) {
		      for(ab=0; ab < FBBints.params->coltot[Gab]; ab++) {
			for(c=0; c < avirtpi[Gc]; c++) {
			  C = avir_off[Gc] + c;
			  if(T1A.params->rowtot[Gk] && T1A.params->coltot[Gk])
			    Gijab.matrix[Gij][ij][ab] += WabC[Gab][ab][c] * T1A.matrix[Gk][k][c];
			}
		      }
		    }
		  } /* Gab */

		  /**** T3 --> Gijab complete ****/

		  /**** T3 --> GiJaB ****/
		  /* Sort W(ab,C) --> Y2(a,bC) */
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FBBints.params->coltot[Gab]; ab++) {
		      A = FBBints.params->colorb[Gab][ab][0];
		      B = FBBints.params->colorb[Gab][ab][1];
		      Ga = FBBints.params->rsym[A];
		      a = A - bvir_off[Ga];
		      for(c=0; c < avirtpi[Gc]; c++) {
			C = avir_off[Gc] + c;
			bc = S2BA.params->colidx[B][C];
			Y2[Ga][a][bc] = WabC[Gab][ab][c];
		      } /* c */
		    } /* ab */
		  } /* Gab */

		  Ga = Gi; Gbc = Ga ^ Gijk;
		  if(T1B.params->rowtot[Gi] && T1B.params->coltot[Gi]) {
		    for(a=0; a < bvirtpi[Ga]; a++) {
		      for(bc=0; bc < GiJaB.params->coltot[Gbc]; bc++) {
			GiJaB.matrix[Gjk][jk][bc] += Y2[Ga][a][bc] * T1B.matrix[Gi][i][a];
		      }
		    }
		  }

		  /**** T3 --> GiJaB complete ****/

		  /**** T3 --> Gijka ****/
		  /* Sort W(AB,c) --> Y1(A,cB) */
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FBBints.params->coltot[Gab]; ab++) {
		      A = FBBints.params->colorb[Gab][ab][0];
		      B = FBBints.params->colorb[Gab][ab][1];
		      Ga = FBBints.params->rsym[A];
		      a = A - bvir_off[Ga];
		      for(c=0; c < avirtpi[Gc]; c++) {
			C = avir_off[Gc] + c;
			cb = T2AB.params->colidx[C][B];
			Y1[Ga][a][cb] = 2 * WabC[Gab][ab][c] + VabC[Gab][ab][c];
		      } /* c */
		    } /* ab */
		  } /* Gab */

		  /* G_ijla <-- t_KlCb Y_ijKabC */
		  for(Gl=0; Gl < nirreps; Gl++) {
		    Ga = Gl ^ Gij;
		    Gkl = Gcb = Gk ^ Gl;

		    nrows = boccpi[Gl];
		    ncols = bvirtpi[Ga];
		    nlinks = T2AB.params->coltot[Gcb];
		    if(nrows && ncols && nlinks) {
		      kl = T2AB.row_offset[Gkl][K];
		      la = Gijka.col_offset[Gij][Gl];
		      C_DGEMM('n','t', nrows, ncols, nlinks, 1.0, T2AB.matrix[Gkl][kl], nlinks,
			      Y1[Ga][0], nlinks, 1.0, &(Gijka.matrix[Gij][ij][la]), ncols);
		    }
		  } /* Gl */

		  /**** T3 --> Gijka complete ****/

		  /**** T3 --> GiJkA ****/
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FBBints.params->coltot[Gab]; ab++) {
		      for(c=0; c < avirtpi[Gc]; c++) {
			XabC[Gab][ab][c] = 2 * WabC[Gab][ab][c] + VabC[Gab][ab][c];
		      } /* c */
		    } /* ab */
		  } /* Gab */

		  /* GiKlC <-- 1/2 t_jlab X_ijKabC */
		  for(Gl=0; Gl < nirreps; Gl++) {
		    Gc = Gl ^ Gik;
		    Gab = Gjl = Gj ^ Gl;
		    nrows = boccpi[Gl];
		    ncols = avirtpi[Gc];
		    nlinks = T2BB.params->coltot[Gjl];
		    if(nrows && ncols && nlinks) {
		      jl = T2BB.row_offset[Gjl][J];
		      lc = GiJkA.col_offset[Gik][Gl];
		      C_DGEMM('n','n', nrows, ncols, nlinks, 0.5, T2BB.matrix[Gjl][jl], nlinks,
			      XabC[Gab][0], ncols, 1.0, &(GiJkA.matrix[Gik][ik][lc]), ncols);
		    }
		  } /* Gl */

		  /**** T3 --> GiJkA complete ****/

		  /**** T3 --> GIjKa ****/
		  /* Sort W(ab,C) --> Y2(a,bC) */
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < FBBints.params->coltot[Gab]; ab++) {
		      A = FBBints.params->colorb[Gab][ab][0];
		      B = FBBints.params->colorb[Gab][ab][1];
		      Ga = FBBints.params->rsym[A];
		      a = A - bvir_off[Ga];
		      for(c=0; c < avirtpi[Gc]; c++) {
			C = avir_off[Gc] + c;
			bc = S2BA.params->colidx[B][C];
			Y2[Ga][a][bc] = 2 * WabC[Gab][ab][c] + VabC[Gab][ab][c];
		      } /* c */
		    } /* ab */
		  } /* Gab */

		  /* G_KiLa <-- -t_jLbC X_ijKabC **/
		  for(Gl=0; Gl < nirreps; Gl++) {
		    Ga = Gki ^ Gl;
		    Gjl = Gbc = Gj ^ Gl;
		    nrows = aoccpi[Gl];
		    ncols = bvirtpi[Ga];
		    nlinks = T2BA.params->coltot[Gbc];
		    if(nrows && ncols && nlinks) {
		      jl = T2BA.row_offset[Gjl][J];
		      la = GIjKa.col_offset[Gki][Gl];
		      C_DGEMM('n','t', nrows, ncols, nlinks, -1.0, T2BA.matrix[Gjl][jl], nlinks,
			      Y2[Ga][0], nlinks, 1.0, &(GIjKa.matrix[Gki][ki][la]), ncols);
		    }
		  } /* Gl */

		  /**** T3 --> GIjKa complete ****/

		  /* Gidab <-- -t_jKdC X_ijKabC */
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Gab = Gid = Gi ^ Gd;
		    Gc = Gjk ^ Gd;

		    nrows = bvirtpi[Gd];
		    ncols = Gidab.params->coltot[Gid];
		    nlinks = avirtpi[Gc];
		    if(nrows && ncols && nlinks) {
		      id = Gidab.row_offset[Gid][I];
		      dc = T2BA.col_offset[Gjk][Gd];
		      C_DGEMM('n','t',nrows, ncols, nlinks, -1.0, &(T2BA.matrix[Gjk][jk][dc]), nlinks,
			      XabC[Gab][0], nlinks, 1.0, Gidab.matrix[Gid][id], ncols);
		    }
		  }
		  /*** T3 --> Gidab complete ****/

		  /* GiDbC <-- t_jKaD t_ijKabC */
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Ga = Gd ^ Gjk;
		    Gid = Gi ^ Gd;

		    nrows = avirtpi[Gd];
		    ncols = GiDaB.params->coltot[Gid];
		    nlinks = bvirtpi[Ga];
		    if(nrows && ncols && nlinks) {
		      ad = T2BA.col_offset[Gjk][Ga];
		      id = GiDaB.row_offset[Gid][I];
		      C_DGEMM('t','n',nrows, ncols, nlinks, -1.0, &(T2BA.matrix[Gjk][jk][ad]), nrows,
			      Y2[Ga][0], ncols, 1.0, GiDaB.matrix[Gid][id], ncols);
		    }
		  }
		  /*** T3 --> GiDaB complete ***/

		  /* GKdCa <-- -1/2 t_ijad t_ijKabC */
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Ga = Gd ^ Gij;
		    Gkd = Gk ^ Gd;

		    nrows = bvirtpi[Gd];
		    ncols = GIdAb.params->coltot[Gkd];
		    nlinks = bvirtpi[Ga];
		    if(nrows && ncols && nlinks) {
		      ad = T2BB.col_offset[Gij][Ga];
		      kd = GIdAb.row_offset[Gkd][K];
		      C_DGEMM('t','n', nrows, ncols, nlinks, 0.5, &(T2BB.matrix[Gij][ij][ad]), nrows,
			      Y1[Ga][0], ncols, 1.0, GIdAb.matrix[Gkd][kd], ncols);
		    }
		  }
		  /*** T3 --> GIdAb complete ***/

		  for(h=0; h < nirreps; h++) {
		    global_dpd_->buf4_mat_irrep_close(&T2BB, h);
		    global_dpd_->buf4_mat_irrep_close(&T2AB, h);
		    global_dpd_->buf4_mat_irrep_close(&T2BA, h);
		    global_dpd_->buf4_mat_irrep_close(&EBBints, h);
		    global_dpd_->buf4_mat_irrep_close(&EABints, h);
		    global_dpd_->buf4_mat_irrep_close(&EBAints, h);
		    global_dpd_->buf4_mat_irrep_close(&DBBints, h);
		    global_dpd_->buf4_mat_irrep_close(&DBAints, h);
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
	      global_dpd_->free_dpd_block(WabC[Gab], FBBints.params->coltot[Gab], avirtpi[Gc]);
	      global_dpd_->free_dpd_block(VabC[Gab], FBBints.params->coltot[Gab], avirtpi[Gc]);
	      global_dpd_->free_dpd_block(XabC[Gab], FBBints.params->coltot[Gab], avirtpi[Gc]);
	    }
	    for(Ga=0; Ga < nirreps; Ga++) {
	      Gbc = Ga ^ Gijk;
	      global_dpd_->free_dpd_block(Y1[Ga], bvirtpi[Ga], FABints.params->coltot[Gbc]);
	      global_dpd_->free_dpd_block(Y2[Ga], bvirtpi[Ga], FBAints.params->coltot[Gbc]);
	    }

	  } /* Gk */
	} /* Gj */
      } /* Gi */

      ET *= 0.25;

      free(WabC);
      free(VabC);
      free(XabC);
      free(Y1); free(Y2);

      global_dpd_->file2_mat_wrt(&DAB);
      global_dpd_->file2_mat_close(&DAB);
      global_dpd_->file2_close(&DAB);
      global_dpd_->file2_mat_wrt(&Dab);
      global_dpd_->file2_mat_close(&Dab);
      global_dpd_->file2_close(&Dab);

      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&S2BB, h);
	global_dpd_->buf4_mat_irrep_close(&S2BB, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&S2BA, h);
	global_dpd_->buf4_mat_irrep_close(&S2BA, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&Gijab, h);
	global_dpd_->buf4_mat_irrep_close(&Gijab, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&GiJaB, h);
	global_dpd_->buf4_mat_irrep_close(&GiJaB, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&Gijka, h);
	global_dpd_->buf4_mat_irrep_close(&Gijka, h);
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
	global_dpd_->buf4_mat_irrep_wrt(&Gidab, h);
	global_dpd_->buf4_mat_irrep_close(&Gidab, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&GiDaB, h);
	global_dpd_->buf4_mat_irrep_close(&GiDaB, h);
      }
      for(h=0; h < nirreps; h++) {
	global_dpd_->buf4_mat_irrep_wrt(&GIdAb, h);
	global_dpd_->buf4_mat_irrep_close(&GIdAb, h);
      }
      global_dpd_->buf4_close(&S2BB);
      /* Combine SIjAb and SiJaB */
      global_dpd_->buf4_sort_axpy(&S2BA, PSIF_CC_MISC, qpsr, 22, 28, "SIjAb", 1);
      global_dpd_->buf4_close(&S2BA);
      global_dpd_->buf4_close(&Gijab);
      /* Combine GIjAb and GiJaB */
      global_dpd_->buf4_sort_axpy(&GiJaB, PSIF_CC_GAMMA, qpsr, 22, 28, "GIjAb", 1);
      global_dpd_->buf4_close(&GiJaB);
      global_dpd_->buf4_close(&Gijka);
      global_dpd_->buf4_close(&GIjKa);
      global_dpd_->buf4_close(&GiJkA);
      global_dpd_->buf4_close(&Gidab);
      global_dpd_->buf4_close(&GiDaB);
      global_dpd_->buf4_close(&GIdAb);

      global_dpd_->file2_mat_wrt(&S1A);
      global_dpd_->file2_mat_close(&S1A);
      global_dpd_->file2_close(&S1A);
      global_dpd_->file2_mat_wrt(&S1B);
      global_dpd_->file2_mat_close(&S1B);
      global_dpd_->file2_close(&S1B);

      global_dpd_->buf4_close(&T2BB);
      global_dpd_->buf4_close(&T2AB);
      global_dpd_->buf4_close(&T2BA);
      global_dpd_->buf4_close(&FBBints);
      global_dpd_->buf4_close(&FABints);
      global_dpd_->buf4_close(&FBAints);
      global_dpd_->buf4_close(&EBBints);
      global_dpd_->buf4_close(&EABints);
      global_dpd_->buf4_close(&EBAints);
      global_dpd_->buf4_close(&DBBints);
      global_dpd_->buf4_close(&DBAints);

      global_dpd_->file2_close(&T1A);
      global_dpd_->file2_close(&T1B);
      global_dpd_->file2_close(&fIJ);
      global_dpd_->file2_close(&fij);
      global_dpd_->file2_close(&fAB);
      global_dpd_->file2_close(&fab);
      global_dpd_->file2_close(&fIA);
      global_dpd_->file2_close(&fia);


      /*** T3 -> DIJ and Dij ***/

      global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
      global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
      global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
      global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
      global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
      global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
      global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
      global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");

      global_dpd_->buf4_init(&T2BB, PSIF_CC_TAMPS, 0, 15, 10, 17, 12, 0, "tabij");
      global_dpd_->buf4_init(&T2AB, PSIF_CC_TAMPS, 0, 28, 22, 28, 22, 0, "tAbIj");
      global_dpd_->buf4_init(&T2BA, PSIF_CC_TAMPS, 0, 29, 23, 29, 23, 0, "taBiJ");
      global_dpd_->buf4_init(&FBBints, PSIF_CC_FINTS, 0, 15, 30, 17, 30, 0, "F <bc||ia>");
      global_dpd_->buf4_init(&FBAints, PSIF_CC_FINTS, 0, 29, 27, 29, 27, 0, "F <bC|iA>");
      global_dpd_->buf4_init(&FABints, PSIF_CC_FINTS, 0, 28, 24, 28, 24, 0, "F <Bc|Ia>");
      global_dpd_->buf4_init(&EBBints, PSIF_CC_EINTS, 0, 31, 10, 31, 12, 0, "E <ak||ij> (ak, i>j)");
      global_dpd_->buf4_init(&EABints, PSIF_CC_EINTS, 0, 25, 22, 25, 22, 0, "E <aK|Ij>");
      global_dpd_->buf4_init(&EBAints, PSIF_CC_EINTS, 0, 26, 23, 26, 23, 0, "E <Ak|iJ>");
      global_dpd_->buf4_init(&DBBints, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
      global_dpd_->buf4_init(&DBAints, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");

      global_dpd_->file2_init(&DIJ, PSIF_CC_OEI, 0, 0, 0, "DIJ");
      global_dpd_->file2_mat_init(&DIJ);
      global_dpd_->file2_mat_rd(&DIJ);
      global_dpd_->file2_init(&Dij, PSIF_CC_OEI, 0, 2, 2, "Dij");
      global_dpd_->file2_mat_init(&Dij);
      global_dpd_->file2_mat_rd(&Dij);

      int Gabc;
      for (Ga=0; Ga < nirreps; ++Ga) {
        for (a=0; a<bvirtpi[Ga]; ++a) {
          A = bvir_off[Ga] + a;
          for (Gb=0; Gb < nirreps; ++Gb) {
            for (b=0; b<bvirtpi[Gb]; ++b) {
              B = bvir_off[Gb] + b;
              for (Gc=0; Gc < nirreps; ++Gc) {
                for (c=0; c < avirtpi[Gc]; ++c) {
                  C = avir_off[Gc] + c;
                  Gabc = Ga ^ Gb ^ Gc;
                  //Allocate the memory for connected and disconnected triples
                  for (Gij=0; Gij < nirreps; ++Gij) {
                    Gk = Gij ^ Gabc;
                    WijK[Gij] = global_dpd_->dpd_block_matrix(T2BB.params->coltot[Gij], aoccpi[Gk]);
                    VijK[Gij] = global_dpd_->dpd_block_matrix(T2BB.params->coltot[Gij], aoccpi[Gk]);
                  }

                  T3_UHF_AAB_abc(WijK, VijK, 1, nirreps, A, Ga, B, Gb, C, Gc,
                      &T2BB, &T2BA, &T2AB, &FBBints, &FBAints, &FABints, &EBBints,
                      &EBAints, &EABints, &T1B, &T1A, &DBBints, &DBAints, &fia,
                      &fIA, &fij, &fIJ, &fab, &fAB, boccpi, bocc_off, aoccpi,
                      aocc_off, bvirtpi, bvir_off, avirtpi, avir_off, 0.0);

                  for(Gi=0; Gi < nirreps; Gi++) {
                    Gj = Gi;
                    Gkl = Gi ^ Gabc;
                    for(Gk=0; Gk < nirreps; Gk++) {
                      Gl = Gk ^ Gkl;
                      Gik = Gjk = Gi ^ Gk;
                      for(i=0; i < boccpi[Gi]; i++) {
                        I = bocc_off[Gi] + i;
                        for(j=0; j < boccpi[Gj]; j++) {
                          J = bocc_off[Gj] + j;
                          for(k=0; k < boccpi[Gk]; k++) {
                            K = bocc_off[Gk] + k;
                            ik = T2BB.params->colidx[I][K];
                            jk = T2BB.params->colidx[J][K];
                            for(l=0; l < aoccpi[Gl]; l++) {
                              Dij.matrix[Gi][i][j] -= 0.5 * WijK[Gik][ik][l] * (WijK[Gjk][jk][l] + VijK[Gjk][jk][l]);
                            } /* l */
                          } /* k */
                        } /* j */
                      } /* i */
                    } /* Gk */
                  } /* Gi */

                  for(Gi=0; Gi < nirreps; Gi++) {
                    Gj = Gi;
                    Gkl = Gi ^ Gabc;
                    for(kl=0; kl < T2BB.params->coltot[Gkl]; kl++) {
                      for(i=0; i < aoccpi[Gi]; i++) {
                        for(j=0; j < aoccpi[Gj]; j++) {
                          DIJ.matrix[Gi][i][j] -= 0.25 * WijK[Gkl][kl][i] * (WijK[Gkl][kl][j] + VijK[Gkl][kl][j]);
                        } /* j */
                      } /* i */
                    } /* kl */
                  } /* Gi */

                  //Deallocate the memory for connected and disconnected triples
                  for (Gij=0; Gij < nirreps; ++Gij) {
                    Gk = Gij ^ Gabc;
                    global_dpd_->free_dpd_block(WijK[Gij], T2BB.params->coltot[Gij], aoccpi[Gk]);
                    global_dpd_->free_dpd_block(VijK[Gij], T2BB.params->coltot[Gij], aoccpi[Gk]);
                  }

                } // c
              } // Gc
            } // b
          } // Gb
        } // a
      } // Ga

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

      global_dpd_->buf4_close(&T2BB);
      global_dpd_->buf4_close(&T2AB);
      global_dpd_->buf4_close(&T2BA);
      global_dpd_->buf4_close(&EBBints);
      global_dpd_->buf4_close(&EABints);
      global_dpd_->buf4_close(&EBAints);
      global_dpd_->buf4_close(&FBBints);
      global_dpd_->buf4_close(&FABints);
      global_dpd_->buf4_close(&FBAints);
      global_dpd_->buf4_close(&DBBints);
      global_dpd_->buf4_close(&DBAints);

      return ET;

    }

  }} // namespace psi::CCTRIPLES
