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

/*! \defgroup CCTRIPLES CCTRIPLES: Evaluate triple excitations */

#include <cstdio>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cctriples {

double ET_AAA(void)
{
  int cnt;
  int h, nirreps;
  int Gi, Gj, Gk, Ga, Gb, Gc, Ge, Gm;
  int Gji, Gij, Gjk, Gik, Gbc, Gac, Gba;
  int I, J, K, A, B, C, E, M;
  int i, j, k, a, b, c, e, m;
  int ij, ji, ik, jk, bc, ac, ba;
  int im, jm, km, ma, mb, mc;
  int ae, be, ce, ke, ie, je;
  int *occpi, *virtpi, *occ_off, *vir_off;
  double value_c, value_d, denom, ET_AAA;
  double t_ijae, t_ijbe, t_ijce, t_jkae, t_jkbe, t_jkce, t_ikae, t_ikbe, t_ikce;
  double F_kebc, F_keac, F_keba, F_iebc, F_ieac, F_ieba, F_jebc, F_jeac, F_jeba;
  double t_imbc, t_imac, t_imba, t_jmbc, t_jmac, t_jmba, t_kmbc, t_kmac, t_kmba;
  double E_jkma, E_jkmb, E_jkmc, E_ikma, E_ikmb, E_ikmc, E_jima, E_jimb, E_jimc;
  double t_ia, t_ib, t_ic, t_ja, t_jb, t_jc, t_ka, t_kb, t_kc;
  double D_jkbc, D_jkac, D_jkba, D_ikbc, D_ikac, D_ikba, D_jibc, D_jiac, D_jiba;
  dpdbuf4 T2, Fints, Eints, Dints;
  dpdfile2 fIJ, fAB, T1;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off;
  vir_off = moinfo.vir_off;

  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_mat_init(&fIJ);
  global_dpd_->file2_mat_init(&fAB);
  global_dpd_->file2_mat_rd(&fIJ);
  global_dpd_->file2_mat_rd(&fAB);

  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&T2, h);
    global_dpd_->buf4_mat_irrep_rd(&T2, h);

    global_dpd_->buf4_mat_irrep_init(&Fints, h);
    global_dpd_->buf4_mat_irrep_rd(&Fints, h);

    global_dpd_->buf4_mat_irrep_init(&Eints, h);
    global_dpd_->buf4_mat_irrep_rd(&Eints, h);

    global_dpd_->buf4_mat_irrep_init(&Dints, h);
    global_dpd_->buf4_mat_irrep_rd(&Dints, h);
  }

  cnt = 0;
  ET_AAA = 0.0;

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gji = Gi ^ Gj;
      for(Gk=0; Gk < nirreps; Gk++) {
	Gjk = Gj ^ Gk;
	Gik = Gi ^ Gk;

	for(Ga=0; Ga < nirreps; Ga++) {
	  for(Gb=0; Gb < nirreps; Gb++) {
	    Gc = Gi ^ Gj ^ Gk ^ Ga ^ Gb;

	    Gbc = Gb^Gc;
	    Gac = Ga^Gc;
	    Gba = Gb^Ga;

	    for(i=0; i < occpi[Gi]; i++) {
	      I = occ_off[Gi] + i;
	      for(j=0; j < occpi[Gj]; j++) {
		J = occ_off[Gj] + j;
		for(k=0; k < occpi[Gk]; k++) {
		  K = occ_off[Gk] + k;

		  ij = T2.params->rowidx[I][J];
		  ji = T2.params->rowidx[J][I];
		  jk = T2.params->rowidx[J][K];
		  ik = T2.params->rowidx[I][K];

		  for(a=0; a < virtpi[Ga]; a++) {
		    A = vir_off[Ga] + a;
		    for(b=0; b < virtpi[Gb]; b++) {
		      B = vir_off[Gb] + b;
		      for(c=0; c < virtpi[Gc]; c++) {
			C = vir_off[Gc] + c;

			bc = Fints.params->colidx[B][C];
			ac = Fints.params->colidx[A][C];
			ba = Fints.params->colidx[B][A];

			value_c = 0.0;

			/** <ov||vv> --> connected triples **/

                        /* -t_jkae * F_iebc */
			Ge = Gj ^ Gk ^ Ga;
			for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

                          ae = T2.params->colidx[A][E];
                          ie = Fints.params->rowidx[I][E];

                          t_jkae = F_iebc = 0.0;

                          if(T2.params->rowtot[Gjk] && T2.params->coltot[Gjk])
			    t_jkae = T2.matrix[Gjk][jk][ae];

                          if(Fints.params->rowtot[Gbc] && Fints.params->coltot[Gbc])
			    F_iebc = Fints.matrix[Gbc][ie][bc];

                          value_c -= t_jkae * F_iebc;
			}

                        /* +t_jkbe * F_ieac */
			Ge = Gj ^ Gk ^ Gb;
			for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

                          be = T2.params->colidx[B][E];
                          ie = Fints.params->rowidx[I][E];

                          t_jkbe = F_ieac = 0.0;

                          if(T2.params->rowtot[Gjk] && T2.params->coltot[Gjk])
			    t_jkbe = T2.matrix[Gjk][jk][be];

                          if(Fints.params->rowtot[Gac] && Fints.params->coltot[Gac])
			    F_ieac = Fints.matrix[Gac][ie][ac];

                          value_c += t_jkbe * F_ieac;
			}

                        /* +t_jkce * F_ieba */
			Ge = Gj ^ Gk ^ Gc;
			for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

                          ce = T2.params->colidx[C][E];
                          ie = Fints.params->rowidx[I][E];

                          t_jkce = F_ieba = 0.0;

                          if(T2.params->rowtot[Gjk] && T2.params->coltot[Gjk])
			    t_jkce = T2.matrix[Gjk][jk][ce];

                          if(Fints.params->rowtot[Gba] && Fints.params->coltot[Gba])
			    F_ieba = Fints.matrix[Gba][ie][ba];

                          value_c += t_jkce * F_ieba;
			}

                        /* +t_ikae * F_jebc */
			Ge = Gi ^ Gk ^ Ga;
			for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

                          ae = T2.params->colidx[A][E];
                          je = Fints.params->rowidx[J][E];

                          t_ikae = F_jebc = 0.0;

                          if(T2.params->rowtot[Gik] && T2.params->coltot[Gik])
			    t_ikae = T2.matrix[Gik][ik][ae];

                          if(Fints.params->rowtot[Gbc] && Fints.params->coltot[Gbc])
			    F_jebc = Fints.matrix[Gbc][je][bc];

                          value_c += t_ikae * F_jebc;
			}

                        /* -t_ikbe * F_jeac */
			Ge = Gi ^ Gk ^ Gb;
			for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

                          be = T2.params->colidx[B][E];
                          je = Fints.params->rowidx[J][E];

                          t_ikbe = F_jeac = 0.0;

                          if(T2.params->rowtot[Gik] && T2.params->coltot[Gik])
			    t_ikbe = T2.matrix[Gik][ik][be];

                          if(Fints.params->rowtot[Gac] && Fints.params->coltot[Gac])
			    F_jeac = Fints.matrix[Gac][je][ac];

                          value_c -= t_ikbe * F_jeac;
			}

                        /* -t_ikce * F_jeba */
			Ge = Gi ^ Gk ^ Gc;
			for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

                          ce = T2.params->colidx[C][E];
                          je = Fints.params->rowidx[J][E];

                          t_ikce = F_jeba = 0.0;

                          if(T2.params->rowtot[Gik] && T2.params->coltot[Gik])
			    t_ikce = T2.matrix[Gik][ik][ce];

                          if(Fints.params->rowtot[Gba] && Fints.params->coltot[Gba])
			    F_jeba = Fints.matrix[Gba][je][ba];

                          value_c -= t_ikce * F_jeba;
			}

                        /* -t_ijae * F_kebc */
			Ge = Gi ^ Gj ^ Ga;
			for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  ae = T2.params->colidx[A][E];
			  ke = Fints.params->rowidx[K][E];

                          t_ijae = F_kebc = 0.0;

                          if(T2.params->rowtot[Gij] && T2.params->coltot[Gij])
			    t_ijae = T2.matrix[Gij][ij][ae];

                          if(Fints.params->rowtot[Gbc] && Fints.params->coltot[Gbc])
			    F_kebc = Fints.matrix[Gbc][ke][bc];

                          value_c -= t_ijae * F_kebc;
			}

                        /* +t_ijbe * F_keac */
			Ge = Gi ^ Gj ^ Gb;
			for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

                          be = T2.params->colidx[B][E];
                          ke = Fints.params->rowidx[K][E];

                          t_ijbe = F_keac = 0.0;

                          if(T2.params->rowtot[Gij] && T2.params->coltot[Gij])
			    t_ijbe = T2.matrix[Gij][ij][be];

                          if(Fints.params->rowtot[Gac] && Fints.params->coltot[Gac])
			    F_keac = Fints.matrix[Gac][ke][ac];

                          value_c += t_ijbe * F_keac;
			}

                        /* +t_ijce * F_keba */
			Ge = Gi ^ Gj ^ Gc;
			for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

                          ce = T2.params->colidx[C][E];
                          ke = Fints.params->rowidx[K][E];

                          t_ijce = F_keba = 0.0;

                          if(T2.params->rowtot[Gij] && T2.params->coltot[Gij])
			    t_ijce = T2.matrix[Gij][ij][ce];

                          if(Fints.params->rowtot[Gba] && Fints.params->coltot[Gba])
			    F_keba = Fints.matrix[Gba][ke][ba];

                          value_c += t_ijce * F_keba;
			}

			/** <oo||ov> --> connected triples **/

                        /* -t_imbc * E_jkma */
			Gm = Gi ^ Gb ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          im = T2.params->rowidx[I][M];
                          ma = Eints.params->colidx[M][A];

                          t_imbc = E_jkma = 0.0;

                          if(T2.params->rowtot[Gbc] && T2.params->coltot[Gbc])
			    t_imbc = T2.matrix[Gbc][im][bc];

                          if(Eints.params->rowtot[Gjk] && Eints.params->coltot[Gjk])
			    E_jkma = Eints.matrix[Gjk][jk][ma];

                          value_c -= t_imbc * E_jkma;
			}

                        /* +t_imac * E_jkmb */
			Gm = Gi ^ Ga ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          im = T2.params->rowidx[I][M];
                          mb = Eints.params->colidx[M][B];

                          t_imac = E_jkmb = 0.0;

                          if(T2.params->rowtot[Gac] && T2.params->coltot[Gac])
			    t_imac = T2.matrix[Gac][im][ac];

                          if(Eints.params->rowtot[Gjk] && Eints.params->coltot[Gjk])
			    E_jkmb = Eints.matrix[Gjk][jk][mb];

                          value_c += t_imac * E_jkmb;
			}

                        /* +t_imba * E_jkmc */
			Gm = Gi ^ Gb ^ Ga;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          im = T2.params->rowidx[I][M];
                          mc = Eints.params->colidx[M][C];

                          t_imba = E_jkmc = 0.0;

                          if(T2.params->rowtot[Gba] && T2.params->coltot[Gba])
			    t_imba = T2.matrix[Gba][im][ba];

                          if(Eints.params->rowtot[Gjk] && Eints.params->coltot[Gjk])
			    E_jkmc = Eints.matrix[Gjk][jk][mc];

                          value_c += t_imba * E_jkmc;
			}

                        /* +t_jmbc * E_ikma */
			Gm = Gj ^ Gb ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          jm = T2.params->rowidx[J][M];
                          ma = Eints.params->colidx[M][A];

                          t_jmbc = E_ikma = 0.0;

                          if(T2.params->rowtot[Gbc] && T2.params->coltot[Gbc])
			    t_jmbc = T2.matrix[Gbc][jm][bc];

                          if(Eints.params->rowtot[Gik] && Eints.params->coltot[Gik])
			    E_ikma = Eints.matrix[Gik][ik][ma];

                          value_c += t_jmbc * E_ikma;
			}

                        /* -t_jmac * E_ikmb */
			Gm = Gj ^ Ga ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          jm = T2.params->rowidx[J][M];
                          mb = Eints.params->colidx[M][B];

                          t_jmac = E_ikmb = 0.0;

                          if(T2.params->rowtot[Gac] && T2.params->coltot[Gac])
			    t_jmac = T2.matrix[Gac][jm][ac];

                          if(Eints.params->rowtot[Gik] && Eints.params->coltot[Gik])
			    E_ikmb = Eints.matrix[Gik][ik][mb];

                          value_c -= t_jmac * E_ikmb;
			}

                        /* -t_jmba * E_ikmc */
			Gm = Gj ^ Gb ^ Ga;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          jm = T2.params->rowidx[J][M];
                          mc = Eints.params->colidx[M][C];

                          t_jmba = E_ikmc = 0.0;

                          if(T2.params->rowtot[Gba] && T2.params->coltot[Gba])
			    t_jmba = T2.matrix[Gba][jm][ba];

                          if(Eints.params->rowtot[Gik] && Eints.params->coltot[Gik])
			    E_ikmc = Eints.matrix[Gik][ik][mc];

                          value_c -= t_jmba * E_ikmc;
			}

                        /* +t_kmbc * E_jima */
			Gm = Gk ^ Gb ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          km = T2.params->rowidx[K][M];
                          ma = Eints.params->colidx[M][A];

                          t_kmbc = E_jima = 0.0;

                          if(T2.params->rowtot[Gbc] && T2.params->coltot[Gbc])
			    t_kmbc = T2.matrix[Gbc][km][bc];

                          if(Eints.params->rowtot[Gji] && Eints.params->coltot[Gji])
			    E_jima = Eints.matrix[Gji][ji][ma];

                          value_c += t_kmbc * E_jima;
			}

                        /* -t_kmac * E_jimb */
			Gm = Gk ^ Ga ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          km = T2.params->rowidx[K][M];
                          mb = Eints.params->colidx[M][B];

                          t_kmac = E_jimb = 0.0;

                          if(T2.params->rowtot[Gac] && T2.params->coltot[Gac])
			    t_kmac = T2.matrix[Gac][km][ac];

                          if(Eints.params->rowtot[Gji] && Eints.params->coltot[Gji])
			    E_jimb = Eints.matrix[Gji][ji][mb];

                          value_c -= t_kmac * E_jimb;
			}

                        /* -t_kmba * E_jimc */
			Gm = Gk ^ Gb ^ Ga;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          km = T2.params->rowidx[K][M];
                          mc = Eints.params->colidx[M][C];

                          t_kmba = E_jimc = 0.0;

                          if(T2.params->rowtot[Gba] && T2.params->coltot[Gba])
			    t_kmba = T2.matrix[Gba][km][ba];

                          if(Eints.params->rowtot[Gji] && Eints.params->coltot[Gji])
			    E_jimc = Eints.matrix[Gji][ji][mc];

                          value_c -= t_kmba * E_jimc;
			}

			/** disconnected triples **/

			value_d = 0.0;

                        /* +t_ia * D_jkbc */
			if(Gi == Ga && Gjk == Gbc) {
			  t_ia = D_jkbc = 0.0;

                          if(T1.params->rowtot[Gi] && T1.params->coltot[Gi])
			    t_ia = T1.matrix[Gi][i][a];

                          if(Dints.params->rowtot[Gjk] && Dints.params->coltot[Gjk])
			    D_jkbc = Dints.matrix[Gjk][jk][bc];

                          value_d += t_ia * D_jkbc;
			}

			/* -t_ib * D_jkac */
			if(Gi == Gb && Gjk == Gac) {
			  t_ib = D_jkac = 0.0;

                          if(T1.params->rowtot[Gi] && T1.params->coltot[Gi])
			    t_ib = T1.matrix[Gi][i][b];

                          if(Dints.params->rowtot[Gjk] && Dints.params->coltot[Gjk])
			    D_jkac = Dints.matrix[Gjk][jk][ac];

                          value_d -= t_ib * D_jkac;
			}

			/* -t_ic * D_jkba */
			if(Gi == Gc && Gjk == Gba) {
			  t_ic = D_jkba = 0.0;

                          if(T1.params->rowtot[Gi] && T1.params->coltot[Gi])
			    t_ic = T1.matrix[Gi][i][c];

                          if(Dints.params->rowtot[Gjk] && Dints.params->coltot[Gjk])
			    D_jkba = Dints.matrix[Gjk][jk][ba];

                          value_d -= t_ic * D_jkba;
			}

                        /* -t_ja * D_ikbc */
			if(Gj == Ga && Gik == Gbc) {
			  t_ja = D_ikbc = 0.0;

                          if(T1.params->rowtot[Gj] && T1.params->coltot[Gj])
			    t_ja = T1.matrix[Gj][j][a];

                          if(Dints.params->rowtot[Gik] && Dints.params->coltot[Gik])
			    D_ikbc = Dints.matrix[Gik][ik][bc];

                          value_d -= t_ja * D_ikbc;
			}

			/* +t_jb * D_ikac */
			if(Gj == Gb && Gik == Gac) {
			  t_jb = D_ikac = 0.0;

                          if(T1.params->rowtot[Gj] && T1.params->coltot[Gj])
			    t_jb = T1.matrix[Gj][j][b];

                          if(Dints.params->rowtot[Gik] && Dints.params->coltot[Gik])
			    D_ikac = Dints.matrix[Gik][ik][ac];

                          value_d += t_jb * D_ikac;
			}

			/* +t_jc * D_ikba */
			if(Gj == Gc && Gik == Gba) {
			  t_jc = D_ikba = 0.0;

                          if(T1.params->rowtot[Gj] && T1.params->coltot[Gj])
			    t_jc = T1.matrix[Gj][j][c];

                          if(Dints.params->rowtot[Gik] && Dints.params->coltot[Gik])
			    D_ikba = Dints.matrix[Gik][ik][ba];

                          value_d += t_jc * D_ikba;
			}

                        /* -t_ka * D_jibc */
			if(Gk == Ga && Gji == Gbc) {
			  t_ka = D_jibc = 0.0;

                          if(T1.params->rowtot[Gk] && T1.params->coltot[Gk])
			    t_ka = T1.matrix[Gk][k][a];

                          if(Dints.params->rowtot[Gji] && Dints.params->coltot[Gji])
			    D_jibc = Dints.matrix[Gji][ji][bc];

                          value_d -= t_ka * D_jibc;
			}

			/* +t_kb * D_jiac */
			if(Gk == Gb && Gji == Gac) {
			  t_kb = D_jiac = 0.0;

                          if(T1.params->rowtot[Gk] && T1.params->coltot[Gk])
			    t_kb = T1.matrix[Gk][k][b];

                          if(Dints.params->rowtot[Gji] && Dints.params->coltot[Gji])
			    D_jiac = Dints.matrix[Gji][ji][ac];

                          value_d += t_kb * D_jiac;
			}

			/* +t_kc * D_jiba */
			if(Gk == Gc && Gji == Gba) {
			  t_kc = D_jiba = 0.0;

                          if(T1.params->rowtot[Gk] && T1.params->coltot[Gk])
			    t_kc = T1.matrix[Gk][k][c];

                          if(Dints.params->rowtot[Gji] && Dints.params->coltot[Gji])
			    D_jiba = Dints.matrix[Gji][ji][ba];

                          value_d += t_kc * D_jiba;
			}

			/*
			if(fabs(value_c) > 1e-7) {
			  cnt++;
			  outfile->Printf( "%d %d %d %d %d %d %20.14f\n", I, J, K, A, B, C, value_c);
			}
			*/

			/* Compute the Fock denominator */
			denom = 0.0;
			if(fIJ.params->rowtot[Gi])
			  denom += fIJ.matrix[Gi][i][i];
			if(fIJ.params->rowtot[Gj])
			  denom += fIJ.matrix[Gj][j][j];
			if(fIJ.params->rowtot[Gk])
			  denom += fIJ.matrix[Gk][k][k];
			if(fAB.params->rowtot[Ga])
			  denom -= fAB.matrix[Ga][a][a];
			if(fAB.params->rowtot[Gb])
			  denom -= fAB.matrix[Gb][b][b];
			if(fAB.params->rowtot[Gc])
			  denom -= fAB.matrix[Gc][c][c];

			ET_AAA += (value_d + value_c) * value_c / denom;


		      } /* c */
		    } /* b */
		  } /* a */

		} /* k */
	      } /* j */
	    } /* i */

	  } /* Gb */
	} /* Ga */

      } /* Gk */
    } /* Gj */
  } /* Gi */

  /*  outfile->Printf( "cnt = %d\n", cnt); */
  ET_AAA /= 36.0;
  /*  outfile->Printf( "ET_AAA = %20.14f\n", ET_AAA); */

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_close(&T2, h);
    global_dpd_->buf4_mat_irrep_close(&Fints, h);
    global_dpd_->buf4_mat_irrep_close(&Eints, h);
    global_dpd_->buf4_mat_irrep_close(&Dints, h);
  }

  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&Fints);
  global_dpd_->buf4_close(&Eints);
  global_dpd_->buf4_close(&Dints);

  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);

  global_dpd_->file2_mat_close(&fIJ);
  global_dpd_->file2_mat_close(&fAB);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fAB);

  return ET_AAA;
}

}} // namespace psi::CCTRIPLES
