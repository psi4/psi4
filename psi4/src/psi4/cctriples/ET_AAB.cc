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
#include <cstdio>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cctriples {

double ET_AAB(void)
{
  int cnt;
  int h, nirreps;
  int Gi, Gj, Gk, Ga, Gb, Gc, Ge, Gm;
  int Gji, Gij, Gjk, Gik, Gbc, Gac, Gba;
  int I, J, K, A, B, C, E, M;
  int i, j, k, a, b, c, e, m;
  int ij, ji, ik, ki, jk, kj;
  int ab, ba, ac, ca, bc, cb;
  int ae, be, ec, ke, ie, je;
  int im, jm, km, mk, ma, mb, mc;
  int *occpi, *virtpi, *occ_off, *vir_off;
  double value_c, value_d, denom, ET_AAB;
  double t_ijae, t_ijbe, t_jkae, t_jkbe, t_jkec, t_ikae, t_ikbe, t_ikec;
  double F_kecb, F_keca, F_iebc, F_ieac, F_ieab, F_jebc, F_jeac, F_jeab;
  double t_imbc, t_imac, t_imba, t_jmbc, t_jmac, t_jmba, t_mkbc, t_mkac;
  double E_kjma, E_kjmb, E_jkmc, E_kima, E_kimb, E_ikmc, E_jima, E_jimb;
  double t_ia, t_ib, t_ja, t_jb, t_kc;
  double D_jkbc, D_jkac, D_ikbc, D_ikac, D_jiba;
  dpdbuf4 T2AB, T2AA, Faints, Fints, Eaints, Eints, Daints, Dints;
  dpdfile2 T1A, T1B, fIJ, fij, fAB, fab;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off;
  vir_off = moinfo.vir_off;

  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");
  global_dpd_->file2_mat_init(&fIJ);
  global_dpd_->file2_mat_init(&fij);
  global_dpd_->file2_mat_init(&fAB);
  global_dpd_->file2_mat_init(&fab);
  global_dpd_->file2_mat_rd(&fIJ);
  global_dpd_->file2_mat_rd(&fij);
  global_dpd_->file2_mat_rd(&fAB);
  global_dpd_->file2_mat_rd(&fab);

  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1A);
  global_dpd_->file2_mat_rd(&T1A);
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->file2_mat_init(&T1B);
  global_dpd_->file2_mat_rd(&T1B);

  global_dpd_->buf4_init(&T2AA, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  global_dpd_->buf4_init(&T2AB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->buf4_init(&Faints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->buf4_init(&Eaints, PSIF_CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  global_dpd_->buf4_init(&Eints, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  global_dpd_->buf4_init(&Daints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&T2AA, h);
    global_dpd_->buf4_mat_irrep_rd(&T2AA, h);

    global_dpd_->buf4_mat_irrep_init(&T2AB, h);
    global_dpd_->buf4_mat_irrep_rd(&T2AB, h);

    global_dpd_->buf4_mat_irrep_init(&Faints, h);
    global_dpd_->buf4_mat_irrep_rd(&Faints, h);

    global_dpd_->buf4_mat_irrep_init(&Fints, h);
    global_dpd_->buf4_mat_irrep_rd(&Fints, h);

    global_dpd_->buf4_mat_irrep_init(&Eaints, h);
    global_dpd_->buf4_mat_irrep_rd(&Eaints, h);

    global_dpd_->buf4_mat_irrep_init(&Eints, h);
    global_dpd_->buf4_mat_irrep_rd(&Eints, h);

    global_dpd_->buf4_mat_irrep_init(&Daints, h);
    global_dpd_->buf4_mat_irrep_rd(&Daints, h);

    global_dpd_->buf4_mat_irrep_init(&Dints, h);
    global_dpd_->buf4_mat_irrep_rd(&Dints, h);
  }

  cnt = 0;
  ET_AAB = 0.0;

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {

	Gij = Gji = Gi ^ Gj;
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

		  ij = T2AA.params->rowidx[I][J];
		  ji = T2AA.params->rowidx[J][I];
		  jk = T2AA.params->rowidx[J][K];
		  kj = T2AA.params->rowidx[K][J];
		  ik = T2AA.params->rowidx[I][K];
		  ki = T2AA.params->rowidx[K][I];

		  for(a=0; a < virtpi[Ga]; a++) {
		    A = vir_off[Ga] + a;
		    for(b=0; b < virtpi[Gb]; b++) {
		      B = vir_off[Gb] + b;
		      for(c=0; c < virtpi[Gc]; c++) {
			C = vir_off[Gc] + c;

			ab = Fints.params->colidx[A][B];
			ba = Fints.params->colidx[B][A];
			bc = Fints.params->colidx[B][C];
			cb = Fints.params->colidx[C][B];
			ac = Fints.params->colidx[A][C];
			ca = Fints.params->colidx[C][A];

			value_c = 0.0;

			/** <ov||vv> --> connected triples **/

                        /* -t_JkAe * F_IeBc */
                        Ge = Gj ^ Gk ^ Ga;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  ae = T2AB.params->colidx[A][E];
			  ie = Fints.params->rowidx[I][E];

                          t_jkae = F_iebc = 0.0;

                          if(T2AB.params->rowtot[Gjk] && T2AB.params->coltot[Gjk])
                            t_jkae = T2AB.matrix[Gjk][jk][ae];

                          if(Fints.params->rowtot[Gbc] && Fints.params->coltot[Gbc])
			    F_iebc = Fints.matrix[Gbc][ie][bc];

                          value_c -= t_jkae * F_iebc;
                        }

                        /* +t_JkBe * F_IeAc */
                        Ge = Gj ^ Gk ^ Gb;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  be = T2AB.params->colidx[B][E];
			  ie = Fints.params->rowidx[I][E];

                          t_jkbe = F_ieac = 0.0;

                          if(T2AB.params->rowtot[Gjk] && T2AB.params->coltot[Gjk])
                            t_jkbe = T2AB.matrix[Gjk][jk][be];

                          if(Fints.params->rowtot[Gac] && Fints.params->coltot[Gac])
			    F_ieac = Fints.matrix[Gac][ie][ac];

                          value_c += t_jkbe * F_ieac;
                        }

                        /* +t_JkEc * F_IEAB */
                        Ge = Gj ^ Gk ^ Gc;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  ec = T2AB.params->colidx[E][C];
			  ie = Faints.params->rowidx[I][E];

                          t_jkec = F_ieab = 0.0;

                          if(T2AB.params->rowtot[Gjk] && T2AB.params->coltot[Gjk])
                            t_jkec = T2AB.matrix[Gjk][jk][ec];

                          if(Faints.params->rowtot[Gba] && Faints.params->coltot[Gba])
			    F_ieab = Faints.matrix[Gba][ie][ab];

                          value_c += t_jkec * F_ieab;
                        }

                        /* +t_IkAe * F_JeBc */
                        Ge = Gi ^ Gk ^ Ga;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  ae = T2AB.params->colidx[A][E];
			  je = Fints.params->rowidx[J][E];

                          t_ikae = F_jebc = 0.0;

                          if(T2AB.params->rowtot[Gik] && T2AB.params->coltot[Gik])
                            t_ikae = T2AB.matrix[Gik][ik][ae];

                          if(Fints.params->rowtot[Gbc] && Fints.params->coltot[Gbc])
			    F_jebc = Fints.matrix[Gbc][je][bc];

                          value_c += t_ikae * F_jebc;
                        }

                        /* -t_IkBe * F_JeAc */
                        Ge = Gi ^ Gk ^ Gb;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  be = T2AB.params->colidx[B][E];
			  je = Fints.params->rowidx[J][E];

                          t_ikbe = F_jeac = 0.0;

                          if(T2AB.params->rowtot[Gik] && T2AB.params->coltot[Gik])
                            t_ikbe = T2AB.matrix[Gik][ik][be];

                          if(Fints.params->rowtot[Gac] && Fints.params->coltot[Gac])
			    F_jeac = Fints.matrix[Gac][je][ac];

                          value_c -= t_ikbe * F_jeac;
                        }

                        /* -t_IkEc * F_JEAB */
                        Ge = Gi ^ Gk ^ Gc;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  ec = T2AB.params->colidx[E][C];
			  je = Faints.params->rowidx[J][E];

                          t_ikec = F_jeab = 0.0;

                          if(T2AB.params->rowtot[Gik] && T2AB.params->coltot[Gik])
                            t_ikec = T2AB.matrix[Gik][ik][ec];

                          if(Faints.params->rowtot[Gba] && Faints.params->coltot[Gba])
			    F_jeab = Faints.matrix[Gba][je][ab];

                          value_c -= t_ikec * F_jeab;
                        }

                        /* +t_IJAE * F_kEcB */
                        Ge = Gi ^ Gj ^ Ga;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  ae = T2AA.params->colidx[A][E];
			  ke = Fints.params->rowidx[K][E];

                          t_ijae = F_kecb = 0.0;

                          if(T2AA.params->rowtot[Gij] && T2AA.params->coltot[Gij])
                            t_ijae = T2AA.matrix[Gij][ij][ae];

                          if(Fints.params->rowtot[Gbc] && Fints.params->coltot[Gbc])
			    F_kecb = Fints.matrix[Gbc][ke][cb];

                          value_c += t_ijae * F_kecb;
                        }

                        /* -t_IJBE * F_kEcA */
                        Ge = Gi ^ Gj ^ Gb;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  be = T2AA.params->colidx[B][E];
			  ke = Fints.params->rowidx[K][E];

                          t_ijbe = F_keca = 0.0;

                          if(T2AA.params->rowtot[Gij] && T2AA.params->coltot[Gij])
                            t_ijbe = T2AA.matrix[Gij][ij][be];

                          if(Fints.params->rowtot[Gac] && Fints.params->coltot[Gac])
			    F_keca = Fints.matrix[Gac][ke][ca];

                          value_c -= t_ijbe * F_keca;
                        }

			/** <oo||ov> --> connected triples **/

                        /* +t_ImBc * E_kJmA */
			Gm = Gi ^ Gb ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          im = T2AB.params->rowidx[I][M];
                          ma = Eints.params->colidx[M][A];

                          t_imbc = E_kjma = 0.0;

                          if(T2AB.params->rowtot[Gbc] && T2AB.params->coltot[Gbc])
			    t_imbc = T2AB.matrix[Gbc][im][bc];

                          if(Eints.params->rowtot[Gjk] && Eints.params->coltot[Gjk])
			    E_kjma = Eints.matrix[Gjk][kj][ma];

                          value_c += t_imbc * E_kjma;
			}

                        /* -t_ImAc * E_kJmB */
			Gm = Gi ^ Ga ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          im = T2AB.params->rowidx[I][M];
                          mb = Eints.params->colidx[M][B];

                          t_imac = E_kjmb = 0.0;

                          if(T2AB.params->rowtot[Gac] && T2AB.params->coltot[Gac])
			    t_imac = T2AB.matrix[Gac][im][ac];

                          if(Eints.params->rowtot[Gjk] && Eints.params->coltot[Gjk])
			    E_kjmb = Eints.matrix[Gjk][kj][mb];

                          value_c -= t_imac * E_kjmb;
			}

                        /* +t_IMBA * E_JkMc */
			Gm = Gi ^ Gb ^ Ga;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          im = T2AA.params->rowidx[I][M];
                          mc = Eints.params->colidx[M][C];

                          t_imba = E_jkmc = 0.0;

                          if(T2AA.params->rowtot[Gba] && T2AA.params->coltot[Gba])
			    t_imba = T2AA.matrix[Gba][im][ba];

                          if(Eints.params->rowtot[Gjk] && Eints.params->coltot[Gjk])
			    E_jkmc = Eints.matrix[Gjk][jk][mc];

                          value_c += t_imba * E_jkmc;
			}

                        /* -t_JmBc * E_kImA */
			Gm = Gj ^ Gb ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          jm = T2AB.params->rowidx[J][M];
                          ma = Eints.params->colidx[M][A];

                          t_jmbc = E_kima = 0.0;

                          if(T2AB.params->rowtot[Gbc] && T2AB.params->coltot[Gbc])
			    t_jmbc = T2AB.matrix[Gbc][jm][bc];

                          if(Eints.params->rowtot[Gik] && Eints.params->coltot[Gik])
			    E_kima = Eints.matrix[Gik][ki][ma];

                          value_c -= t_jmbc * E_kima;
			}

                        /* +t_JmAc * E_kImB */
			Gm = Gj ^ Ga ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          jm = T2AB.params->rowidx[J][M];
                          mb = Eints.params->colidx[M][B];

                          t_jmac = E_kimb = 0.0;

                          if(T2AB.params->rowtot[Gac] && T2AB.params->coltot[Gac])
			    t_jmac = T2AB.matrix[Gac][jm][ac];

                          if(Eints.params->rowtot[Gik] && Eints.params->coltot[Gik])
			    E_kimb = Eints.matrix[Gik][ki][mb];

                          value_c += t_jmac * E_kimb;
			}

                        /* -t_JMBA * E_IkMc */
			Gm = Gj ^ Gb ^ Ga;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          jm = T2AA.params->rowidx[J][M];
                          mc = Eints.params->colidx[M][C];

                          t_jmba = E_ikmc = 0.0;

                          if(T2AA.params->rowtot[Gba] && T2AA.params->coltot[Gba])
			    t_jmba = T2AA.matrix[Gba][jm][ba];

                          if(Eints.params->rowtot[Gik] && Eints.params->coltot[Gik])
			    E_ikmc = Eints.matrix[Gik][ik][mc];

                          value_c -= t_jmba * E_ikmc;
			}

                        /* -t_MkBc * E_JIMA */
			Gm = Gk ^ Gb ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          mk = T2AB.params->rowidx[M][K];
                          ma = Eints.params->colidx[M][A];

                          t_mkbc = E_jima = 0.0;

                          if(T2AB.params->rowtot[Gbc] && T2AB.params->coltot[Gbc])
			    t_mkbc = T2AB.matrix[Gbc][mk][bc];

                          if(Eaints.params->rowtot[Gji] && Eaints.params->coltot[Gji])
			    E_jima = Eaints.matrix[Gji][ji][ma];

                          value_c -= t_mkbc * E_jima;
			}

                        /* +t_MkAc * E_JIMB */
			Gm = Gk ^ Ga ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          mk = T2AB.params->rowidx[M][K];
                          mb = Eints.params->colidx[M][B];

                          t_mkac = E_jimb = 0.0;

                          if(T2AB.params->rowtot[Gac] && T2AB.params->coltot[Gac])
			    t_mkac = T2AB.matrix[Gac][mk][ac];

                          if(Eaints.params->rowtot[Gji] && Eaints.params->coltot[Gji])
			    E_jimb = Eaints.matrix[Gji][ji][mb];

                          value_c += t_mkac * E_jimb;
			}

			/** disconnected triples **/

			value_d = 0.0;

                        /* +t_IA * D_JkBc */
			if(Gi == Ga && Gjk == Gbc) {
			  t_ia = D_jkbc = 0.0;

                          if(T1A.params->rowtot[Gi] && T1A.params->coltot[Gi])
			    t_ia = T1A.matrix[Gi][i][a];

                          if(Dints.params->rowtot[Gjk] && Dints.params->coltot[Gjk])
			    D_jkbc = Dints.matrix[Gjk][jk][bc];

                          value_d += t_ia * D_jkbc;
			}

			/* -t_IB * D_JkAc */
			if(Gi == Gb && Gjk == Gac) {
			  t_ib = D_jkac = 0.0;

                          if(T1A.params->rowtot[Gi] && T1A.params->coltot[Gi])
			    t_ib = T1A.matrix[Gi][i][b];

                          if(Dints.params->rowtot[Gjk] && Dints.params->coltot[Gjk])
			    D_jkac = Dints.matrix[Gjk][jk][ac];

                          value_d -= t_ib * D_jkac;
			}

                        /* -t_JA * D_IkBc */
			if(Gj == Ga && Gik == Gbc) {
			  t_ja = D_ikbc = 0.0;

                          if(T1A.params->rowtot[Gj] && T1A.params->coltot[Gj])
			    t_ja = T1A.matrix[Gj][j][a];

                          if(Dints.params->rowtot[Gik] && Dints.params->coltot[Gik])
			    D_ikbc = Dints.matrix[Gik][ik][bc];

                          value_d -= t_ja * D_ikbc;
			}

			/* +t_JB * D_IkAc */
			if(Gj == Gb && Gik == Gac) {
			  t_jb = D_ikac = 0.0;

                          if(T1A.params->rowtot[Gj] && T1A.params->coltot[Gj])
			    t_jb = T1A.matrix[Gj][j][b];

                          if(Dints.params->rowtot[Gik] && Dints.params->coltot[Gik])
			    D_ikac = Dints.matrix[Gik][ik][ac];

                          value_d += t_jb * D_ikac;
			}

			/* +t_kc * D_JIBA */
			if(Gk == Gc && Gji == Gba) {
			  t_kc = D_jiba = 0.0;

                          if(T1B.params->rowtot[Gk] && T1B.params->coltot[Gk])
			    t_kc = T1B.matrix[Gk][k][c];

                          if(Daints.params->rowtot[Gji] && Daints.params->coltot[Gji])
			    D_jiba = Daints.matrix[Gji][ji][ba];

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
			if(fij.params->rowtot[Gk])
			  denom += fij.matrix[Gk][k][k];
			if(fAB.params->rowtot[Ga])
			  denom -= fAB.matrix[Ga][a][a];
			if(fAB.params->rowtot[Gb])
			  denom -= fAB.matrix[Gb][b][b];
			if(fab.params->rowtot[Gc])
			  denom -= fab.matrix[Gc][c][c];

			ET_AAB += (value_d + value_c) * value_c / denom;

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
  ET_AAB /= 4.0;
  /*  outfile->Printf( "ET_AAB = %20.14f\n", ET_AAB); */

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_close(&T2AA, h);
    global_dpd_->buf4_mat_irrep_close(&T2AB, h);
    global_dpd_->buf4_mat_irrep_close(&Faints, h);
    global_dpd_->buf4_mat_irrep_close(&Fints, h);
    global_dpd_->buf4_mat_irrep_close(&Eaints, h);
    global_dpd_->buf4_mat_irrep_close(&Eints, h);
    global_dpd_->buf4_mat_irrep_close(&Daints, h);
    global_dpd_->buf4_mat_irrep_close(&Dints, h);
  }

  global_dpd_->buf4_close(&T2AA);
  global_dpd_->buf4_close(&T2AB);
  global_dpd_->buf4_close(&Faints);
  global_dpd_->buf4_close(&Fints);
  global_dpd_->buf4_close(&Eaints);
  global_dpd_->buf4_close(&Eints);
  global_dpd_->buf4_close(&Daints);
  global_dpd_->buf4_close(&Dints);

  global_dpd_->file2_mat_close(&T1A);
  global_dpd_->file2_close(&T1A);
  global_dpd_->file2_mat_close(&T1B);
  global_dpd_->file2_close(&T1B);

  global_dpd_->file2_mat_close(&fIJ);
  global_dpd_->file2_mat_close(&fij);
  global_dpd_->file2_mat_close(&fAB);
  global_dpd_->file2_mat_close(&fab);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fij);
  global_dpd_->file2_close(&fAB);
  global_dpd_->file2_close(&fab);

  return ET_AAB;
}

}} // namespace psi::CCTRIPLES
