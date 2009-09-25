/*! \file
    \ingroup CCTRIPLES
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cctriples {

double ET_ABB(void)
{
  int cnt;
  int h, nirreps;
  int Gi, Gj, Gk, Ga, Gb, Gc, Ge, Gm;
  int Gji, Gij, Gjk, Gik, Gbc, Gac, Gba, Gab;
  int I, J, K, A, B, C, E, M;
  int i, j, k, a, b, c, e, m;
  int ij, ji, ik, ki, jk;
  int ab, ba, ac, ca, bc, cb;
  int ae, be, eb, ce, ec, ie, je, ke;
  int im, jm, mj, km, mk, ma, mb, mc;
  int *occpi, *virtpi, *occ_off, *vir_off;
  double value_c, value_d, denom, ET_ABB;
  double t_ijae, t_ijeb, t_ijec, t_jkbe, t_jkce, t_ikae, t_ikeb, t_ikec;
  double F_kebc, F_keca, F_keba, F_ieac, F_ieab, F_jebc, F_jeca, F_jeba;
  double t_imac, t_imab, t_jmbc, t_mjac, t_mjab, t_kmbc, t_mkac, t_mkab;
  double E_jkmb, E_jkmc, E_kima, E_ikmb, E_ikmc, E_jima, E_ijmb, E_ijmc;
  double t_ia, t_jb, t_jc, t_kb, t_kc;
  double D_jkbc, D_ikac, D_ikab, D_ijac, D_ijab;
  dpdbuf4 T2AB, T2BB, Faints, Fints, Eaints, Eints, Daints, Dints;
  dpdfile2 T1A, T1B, fIJ, fij, fAB, fab;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off;
  vir_off = moinfo.vir_off;

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fij");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_init(&fij);
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_init(&fab);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_mat_rd(&fij);
  dpd_file2_mat_rd(&fAB);
  dpd_file2_mat_rd(&fab);

  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1A);
  dpd_file2_mat_rd(&T1A);
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_file2_mat_init(&T1B);
  dpd_file2_mat_rd(&T1B);

  dpd_buf4_init(&T2BB, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&Faints, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&Eaints, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&Daints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2BB, h);
    dpd_buf4_mat_irrep_rd(&T2BB, h);

    dpd_buf4_mat_irrep_init(&T2AB, h);
    dpd_buf4_mat_irrep_rd(&T2AB, h);

    dpd_buf4_mat_irrep_init(&Faints, h);
    dpd_buf4_mat_irrep_rd(&Faints, h);

    dpd_buf4_mat_irrep_init(&Fints, h);
    dpd_buf4_mat_irrep_rd(&Fints, h);

    dpd_buf4_mat_irrep_init(&Eaints, h);
    dpd_buf4_mat_irrep_rd(&Eaints, h);

    dpd_buf4_mat_irrep_init(&Eints, h);
    dpd_buf4_mat_irrep_rd(&Eints, h);

    dpd_buf4_mat_irrep_init(&Daints, h);
    dpd_buf4_mat_irrep_rd(&Daints, h);

    dpd_buf4_mat_irrep_init(&Dints, h);
    dpd_buf4_mat_irrep_rd(&Dints, h);
  }

  cnt = 0;
  ET_ABB = 0.0;

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

		  ij = T2BB.params->rowidx[I][J];
		  ji = T2BB.params->rowidx[J][I];
		  jk = T2BB.params->rowidx[J][K];
		  ik = T2BB.params->rowidx[I][K];
		  ki = T2BB.params->rowidx[K][I];

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

                        /* +t_jkbe * F_IeAc */
                        Ge = Gj ^ Gk ^ Gb;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  be = T2BB.params->colidx[B][E];
			  ie = Fints.params->rowidx[I][E];

                          t_jkbe = F_ieac = 0.0;

                          if(T2BB.params->rowtot[Gjk] && T2BB.params->coltot[Gjk])
                            t_jkbe = T2BB.matrix[Gjk][jk][be];

                          if(Fints.params->rowtot[Gac] && Fints.params->coltot[Gac])
			    F_ieac = Fints.matrix[Gac][ie][ac];
 
                          value_c += t_jkbe * F_ieac;
                        }

                        /* -t_jkce * F_IeAb */
                        Ge = Gj ^ Gk ^ Gc;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  ce = T2BB.params->colidx[C][E];
			  ie = Faints.params->rowidx[I][E];

                          t_jkce = F_ieab = 0.0;

                          if(T2BB.params->rowtot[Gjk] && T2BB.params->coltot[Gjk])
                            t_jkce = T2BB.matrix[Gjk][jk][ce];

                          if(Fints.params->rowtot[Gba] && Fints.params->coltot[Gba])
			    F_ieab = Fints.matrix[Gba][ie][ab];
 
                          value_c -= t_jkce * F_ieab;
                        }

                        /* +t_IkAe * F_jebc */
                        Ge = Gi ^ Gk ^ Ga;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  ae = T2AB.params->colidx[A][E];
			  je = Faints.params->rowidx[J][E];

                          t_ikae = F_jebc = 0.0;

                          if(T2AB.params->rowtot[Gik] && T2AB.params->coltot[Gik])
                            t_ikae = T2AB.matrix[Gik][ik][ae];

                          if(Faints.params->rowtot[Gbc] && Faints.params->coltot[Gbc])
			    F_jebc = Faints.matrix[Gbc][je][bc];
 
                          value_c += t_ikae * F_jebc;
                        }

                        /* -t_IkEb * F_jEcA */
                        Ge = Gi ^ Gk ^ Gb;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  eb = T2AB.params->colidx[E][B];
			  je = Fints.params->rowidx[J][E];

                          t_ikeb = F_jeca = 0.0;

                          if(T2AB.params->rowtot[Gik] && T2AB.params->coltot[Gik])
                            t_ikeb = T2AB.matrix[Gik][ik][eb];

                          if(Fints.params->rowtot[Gac] && Fints.params->coltot[Gac])
			    F_jeca = Fints.matrix[Gac][je][ca];
 
                          value_c -= t_ikeb * F_jeca;
                        }

                        /* +t_IkEc * F_jEbA */
                        Ge = Gi ^ Gk ^ Gc;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  ec = T2AB.params->colidx[E][C];
			  je = Fints.params->rowidx[J][E];

                          t_ikec = F_jeba = 0.0;

                          if(T2AB.params->rowtot[Gik] && T2AB.params->coltot[Gik])
                            t_ikec = T2AB.matrix[Gik][ik][ec];

                          if(Fints.params->rowtot[Gba] && Fints.params->coltot[Gba])
			    F_jeba = Fints.matrix[Gba][je][ba];
 
                          value_c += t_ikec * F_jeba;
                        }

                        /* -t_IjAe * F_kebc */
                        Ge = Gi ^ Gj ^ Ga;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  ae = T2AB.params->colidx[A][E];
			  ke = Faints.params->rowidx[K][E];

                          t_ijae = F_kebc = 0.0;

                          if(T2AB.params->rowtot[Gij] && T2AB.params->coltot[Gij])
                            t_ijae = T2AB.matrix[Gij][ij][ae];

                          if(Faints.params->rowtot[Gbc] && Faints.params->coltot[Gbc])
			    F_kebc = Faints.matrix[Gbc][ke][bc];
 
                          value_c -= t_ijae * F_kebc;
                        }

                        /* +t_IjEb * F_kEcA */
                        Ge = Gi ^ Gj ^ Gb;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  eb = T2AB.params->colidx[E][B];
			  ke = Fints.params->rowidx[K][E];

                          t_ijeb = F_keca = 0.0;

                          if(T2AB.params->rowtot[Gij] && T2AB.params->coltot[Gij])
                            t_ijeb = T2AB.matrix[Gij][ij][eb];

                          if(Fints.params->rowtot[Gac] && Fints.params->coltot[Gac])
			    F_keca = Fints.matrix[Gac][ke][ca];
 
                          value_c += t_ijeb * F_keca;
                        }

                        /* -t_IjEc * F_kEbA */
                        Ge = Gi ^ Gj ^ Gc;
                        for(e=0; e < virtpi[Ge]; e++) {
                          E = vir_off[Ge] + e;

			  ec = T2AB.params->colidx[E][C];
			  ke = Fints.params->rowidx[K][E];

                          t_ijec = F_keba = 0.0;

                          if(T2AB.params->rowtot[Gij] && T2AB.params->coltot[Gij])
                            t_ijec = T2AB.matrix[Gij][ij][ec];

                          if(Fints.params->rowtot[Gba] && Fints.params->coltot[Gba])
			    F_keba = Fints.matrix[Gba][ke][ba];
 
                          value_c -= t_ijec * F_keba;
                        }

			/** <oo||ov> --> connected triples **/

                        /* +t_ImAc * E_jkmb */
			Gm = Gi ^ Ga ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          im = T2AB.params->rowidx[I][M];
                          mb = Eaints.params->colidx[M][B];

                          t_imac = E_jkmb = 0.0;

                          if(T2AB.params->rowtot[Gac] && T2AB.params->coltot[Gac])
			    t_imac = T2AB.matrix[Gac][im][ac];

                          if(Eaints.params->rowtot[Gjk] && Eaints.params->coltot[Gjk])
			    E_jkmb = Eaints.matrix[Gjk][jk][mb];

                          value_c += t_imac * E_jkmb;
			}

                        /* -t_ImAb * E_jkmc */
			Gm = Gi ^ Gb ^ Ga;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          im = T2AB.params->rowidx[I][M];
                          mc = Eaints.params->colidx[M][C];

                          t_imab = E_jkmc = 0.0;

                          if(T2AB.params->rowtot[Gba] && T2AB.params->coltot[Gba])
			    t_imab = T2AB.matrix[Gba][im][ab];

                          if(Eaints.params->rowtot[Gjk] && Eaints.params->coltot[Gjk])
			    E_jkmc = Eaints.matrix[Gjk][jk][mc];

                          value_c -= t_imab * E_jkmc;
			}

                        /* -t_jmbc * E_kImA */
			Gm = Gj ^ Gb ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          jm = T2BB.params->rowidx[J][M];
                          ma = Eints.params->colidx[M][A];

                          t_jmbc = E_kima = 0.0;

                          if(T2BB.params->rowtot[Gbc] && T2BB.params->coltot[Gbc])
			    t_jmbc = T2BB.matrix[Gbc][jm][bc];

                          if(Eints.params->rowtot[Gik] && Eints.params->coltot[Gik])
			    E_kima = Eints.matrix[Gik][ki][ma];

                          value_c -= t_jmbc * E_kima;
			}

                        /* +t_MjAc * E_IkMb */
			Gm = Gj ^ Ga ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          mj = T2AB.params->rowidx[M][J];
                          mb = Eints.params->colidx[M][B];

                          t_mjac = E_ikmb = 0.0;

                          if(T2AB.params->rowtot[Gac] && T2AB.params->coltot[Gac])
			    t_mjac = T2AB.matrix[Gac][mj][ac];

                          if(Eints.params->rowtot[Gik] && Eints.params->coltot[Gik])
			    E_ikmb = Eints.matrix[Gik][ik][mb];

                          value_c += t_mjac * E_ikmb;
			}

                        /* -t_MjAb * E_IkMc */
			Gm = Gj ^ Gb ^ Ga;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          mj = T2AB.params->rowidx[M][J];
                          mc = Eints.params->colidx[M][C];

                          t_mjab = E_ikmc = 0.0;

                          if(T2AB.params->rowtot[Gba] && T2AB.params->coltot[Gba])
			    t_mjab = T2AB.matrix[Gba][mj][ab];

                          if(Eints.params->rowtot[Gik] && Eints.params->coltot[Gik])
			    E_ikmc = Eints.matrix[Gik][ik][mc];

                          value_c -= t_mjab * E_ikmc;
			}

                        /* +t_kmbc * E_jImA */
			Gm = Gk ^ Gb ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          km = T2BB.params->rowidx[K][M];
                          ma = Eints.params->colidx[M][A];

                          t_kmbc = E_jima = 0.0;

                          if(T2BB.params->rowtot[Gbc] && T2BB.params->coltot[Gbc])
			    t_kmbc = T2BB.matrix[Gbc][km][bc];

                          if(Eints.params->rowtot[Gji] && Eints.params->coltot[Gji])
			    E_jima = Eints.matrix[Gji][ji][ma];

                          value_c += t_kmbc * E_jima;
			}

                        /* -t_MkAc * E_IjMb */
			Gm = Gk ^ Ga ^ Gc;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          mk = T2AB.params->rowidx[M][K];
                          mb = Eints.params->colidx[M][B];

                          t_mkac = E_ijmb = 0.0;

                          if(T2AB.params->rowtot[Gac] && T2AB.params->coltot[Gac])
			    t_mkac = T2AB.matrix[Gac][mk][ac];

                          if(Eints.params->rowtot[Gji] && Eints.params->coltot[Gji])
			    E_ijmb = Eints.matrix[Gji][ij][mb];

                          value_c -= t_mkac * E_ijmb;
			}

                        /* +t_MkAb * E_IjMc */
			Gm = Gk ^ Gb ^ Ga;
			for(m=0; m < occpi[Gm]; m++) {
                          M = occ_off[Gm] + m;

                          mk = T2AB.params->rowidx[M][K];
                          mc = Eints.params->colidx[M][C];

                          t_mkab = E_ijmc = 0.0;

                          if(T2AB.params->rowtot[Gba] && T2AB.params->coltot[Gba])
			    t_mkab = T2AB.matrix[Gba][mk][ab];

                          if(Eints.params->rowtot[Gji] && Eints.params->coltot[Gji])
			    E_ijmc = Eints.matrix[Gji][ij][mc];

                          value_c += t_mkab * E_ijmc;
			}

			/** disconnected triples **/

			value_d = 0.0;

                        /* +t_IA * D_jkbc */
			if(Gi == Ga && Gjk == Gbc) {
			  t_ia = D_jkbc = 0.0;

                          if(T1A.params->rowtot[Gi] && T1A.params->coltot[Gi])
			    t_ia = T1A.matrix[Gi][i][a];

                          if(Daints.params->rowtot[Gjk] && Daints.params->coltot[Gjk])
			    D_jkbc = Daints.matrix[Gjk][jk][bc];

                          value_d += t_ia * D_jkbc;
			}

			/* +t_jb * D_IkAc */
			if(Gj == Gb && Gik == Gac) {
			  t_jb = D_ikac = 0.0;

                          if(T1B.params->rowtot[Gj] && T1B.params->coltot[Gj])
			    t_jb = T1B.matrix[Gj][j][b];

                          if(Dints.params->rowtot[Gik] && Dints.params->coltot[Gik])
			    D_ikac = Dints.matrix[Gik][ik][ac];

                          value_d += t_jb * D_ikac;
			}

			/* -t_jc * D_IkAb */
			if(Gj == Gc && Gik == Gba) {
			  t_jc = D_ikab = 0.0;

                          if(T1B.params->rowtot[Gj] && T1B.params->coltot[Gj])
			    t_jc = T1B.matrix[Gj][j][c];

                          if(Dints.params->rowtot[Gik] && Dints.params->coltot[Gik])
			    D_ikab = Dints.matrix[Gik][ik][ab];

                          value_d -= t_jc * D_ikab;
			}

			/* -t_kb * D_IjAc */
			if(Gk == Gb && Gji == Gac) {
			  t_kb = D_ijac = 0.0;

                          if(T1B.params->rowtot[Gk] && T1B.params->coltot[Gk])
			    t_kb = T1B.matrix[Gk][k][b];

                          if(Dints.params->rowtot[Gji] && Dints.params->coltot[Gji])
			    D_ijac = Dints.matrix[Gji][ij][ac];

                          value_d -= t_kb * D_ijac;
			}

			/* +t_kc * D_IjAb */
			if(Gk == Gc && Gji == Gba) {
			  t_kc = D_ijab = 0.0;

                          if(T1B.params->rowtot[Gk] && T1B.params->coltot[Gk])
			    t_kc = T1B.matrix[Gk][k][c];

                          if(Dints.params->rowtot[Gji] && Dints.params->coltot[Gji])
			    D_ijab = Dints.matrix[Gji][ij][ab];

                          value_d += t_kc * D_ijab;
			}

			/*
			if(fabs(value_c) > 1e-7) {
			  cnt++;
			  fprintf(outfile, "%d %d %d %d %d %d %20.14f\n", I, J, K, A, B, C, value_c);
			}
			*/

			/* Compute the Fock denominator */
			denom = 0.0;
			if(fIJ.params->rowtot[Gi])
			  denom += fIJ.matrix[Gi][i][i];
			if(fij.params->rowtot[Gj])
			  denom += fij.matrix[Gj][j][j];
			if(fij.params->rowtot[Gk])
			  denom += fij.matrix[Gk][k][k];
			if(fAB.params->rowtot[Ga])
			  denom -= fAB.matrix[Ga][a][a];
			if(fab.params->rowtot[Gb])
			  denom -= fab.matrix[Gb][b][b];
			if(fab.params->rowtot[Gc])
			  denom -= fab.matrix[Gc][c][c];

			ET_ABB += (value_d + value_c) * value_c / denom;


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

  /*  fprintf(outfile, "cnt = %d\n", cnt); */
  ET_ABB /= 4.0;
  /*  fprintf(outfile, "ET_ABB = %20.14f\n", ET_ABB); */


  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&T2BB, h);
    dpd_buf4_mat_irrep_close(&T2AB, h);
    dpd_buf4_mat_irrep_close(&Faints, h);
    dpd_buf4_mat_irrep_close(&Fints, h);
    dpd_buf4_mat_irrep_close(&Eaints, h);
    dpd_buf4_mat_irrep_close(&Eints, h);
    dpd_buf4_mat_irrep_close(&Daints, h);
    dpd_buf4_mat_irrep_close(&Dints, h);
  }

  dpd_buf4_close(&T2BB);
  dpd_buf4_close(&T2AB);
  dpd_buf4_close(&Faints);
  dpd_buf4_close(&Fints);
  dpd_buf4_close(&Eaints);
  dpd_buf4_close(&Eints);
  dpd_buf4_close(&Daints);
  dpd_buf4_close(&Dints);

  dpd_file2_mat_close(&T1A);
  dpd_file2_close(&T1A);
  dpd_file2_mat_close(&T1B);
  dpd_file2_close(&T1B);

  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fij);
  dpd_file2_mat_close(&fAB);
  dpd_file2_mat_close(&fab);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fij);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fab);

  return ET_ABB;
}

}} // namespace psi::cctriples
