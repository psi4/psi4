/*! \file
    \ingroup CCTRIPLES
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"



namespace psi { namespace cctriples {

extern void T3_UHF_AAA(double ***W, double ***V, int disc, int nirreps, 
    int I, int Gi, int J, int Gj, int K, int Gk, dpdbuf4 *C2, dpdbuf4 *F, 
    dpdbuf4 *E, dpdfile2 *C1, dpdbuf4 *D, dpdfile2 *fIA, dpdfile2 *fIJ, 
    dpdfile2 *fAB, int *occpi, int *occ_off, int *virtpi, int *vir_off, 
    double omega);

extern void T3_UHF_AAA_abc(double ***W, double ***V, int disc, int nirreps,
    int A, int Ga, int B, int Gb, int C, int Gc, dpdbuf4 *C2, dpdbuf4 *F,
    dpdbuf4 *E, dpdfile2 *C1, dpdbuf4 *D, dpdfile2 *fIA, dpdfile2 *fIJ,
    dpdfile2 *fAB, int *occpi, int *occ_off, int *virtpi, int *vir_off,
    double omega);

    double T3_grad_UHF_AAA(void)
    {
      int h, nirreps;
      int *occpi, *virtpi, *occ_off, *vir_off;
      int i, j, k, a, b, c, d, l;
      int I, J, K, A, B, C, D, L;
      int ij, ji, ik, ki, jk, kj;
      int ab, ba, ac, ca, bc, cb;
      int il, jl, kl, li;
      int id, jd, kd;
      int ad, bd, cd, dc;
      int la, lb, lc;
      int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
      int Gij, Gji, Gik, Gki, Gjk, Gkj, Gijk;
      int Gid, Gjd, Gkd, Gil, Gjl, Gkl, Gli;
      int Gab, Gba, Gac, Gca, Gbc, Gcb, Gcd;
      int ncols, nrows, nlinks;
      double value_c, value_d, dijk, denom, ET;
      double t_ia, t_ib, t_ic, t_ja, t_jb, t_jc, t_ka, t_kb, t_kc;
      double f_ia, f_ib, f_ic, f_ja, f_jb, f_jc, f_ka, f_kb, f_kc;
      double D_jkbc, D_jkac, D_jkba, D_ikbc, D_ikac, D_ikba, D_jibc, D_jiac, D_jiba;
      double t_jkbc, t_jkac, t_jkba, t_ikbc, t_ikac, t_ikba, t_jibc, t_jiac, t_jiba;
      dpdbuf4 T2, Fints, Eints, Dints, S2, GIJAB, GIJKA, GIDAB;
      dpdfile2 fIJ, fAB, fIA, T1, S1, DAB, DIJ;
      double ***WABC, ***VABC, ***XABC, ***Y;
      double **Z;
      double ***WIJK = (double ***) malloc(nirreps * sizeof(double **));
      double ***VIJK = (double ***) malloc(nirreps * sizeof(double **));

      nirreps = moinfo.nirreps;
      occpi = moinfo.aoccpi; 
      virtpi = moinfo.avirtpi;
      occ_off = moinfo.aocc_off;
      vir_off = moinfo.avir_off;

      dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
      dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
      dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");

      dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
      dpd_buf4_init(&Fints, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
      dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
      dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");

      dpd_file2_init(&S1, CC_OEI, 0, 0, 1, "SIA");
      dpd_file2_mat_init(&S1); 
      dpd_buf4_init(&S2, CC_MISC, 0, 0, 5, 2, 7, 0, "SIJAB");
      for(h=0; h < nirreps; h++)
	dpd_buf4_mat_irrep_init(&S2, h);

      dpd_file2_init(&DAB, CC_OEI, 0, 1, 1, "DAB");
      dpd_file2_mat_init(&DAB);

      dpd_buf4_init(&GIJAB, CC_GAMMA, 0, 0, 5, 2, 7, 0, "GIJAB");
      for(h=0; h < nirreps; h++)
	dpd_buf4_mat_irrep_init(&GIJAB, h);

      dpd_buf4_init(&GIJKA, CC_GAMMA, 0, 0, 20, 2, 20, 0, "GIJKA");
      for(h=0; h < nirreps; h++)
	dpd_buf4_mat_irrep_init(&GIJKA, h);

      dpd_buf4_init(&GIDAB, CC_GAMMA, 0, 20, 5, 20, 7, 0, "GIDAB");
      for(h=0; h < nirreps; h++)
	dpd_buf4_mat_irrep_init(&GIDAB, h);

      WABC = (double ***) malloc(nirreps * sizeof(double **));
      VABC = (double ***) malloc(nirreps * sizeof(double **));
      XABC = (double ***) malloc(nirreps * sizeof(double **));
      Y = (double ***) malloc(nirreps * sizeof(double **));

      ET = 0.0;

      for(Gi=0; Gi < nirreps; Gi++) {
	for(Gj=0; Gj < nirreps; Gj++) {
	  for(Gk=0; Gk < nirreps; Gk++) {

	    Gij = Gji = Gi ^ Gj;
	    Gjk = Gkj = Gj ^ Gk;
	    Gik = Gki = Gi ^ Gk;

	    Gijk = Gi ^ Gj ^ Gk;

	    for(Gab=0; Gab < nirreps; Gab++) {
	      Gc = Gab ^ Gijk;
	      WABC[Gab] = dpd_block_matrix(Fints.params->coltot[Gab], virtpi[Gc]);
	      VABC[Gab] = dpd_block_matrix(Fints.params->coltot[Gab], virtpi[Gc]);
	      XABC[Gab] = dpd_block_matrix(Fints.params->coltot[Gab], virtpi[Gc]);
	    }
	    for(Ga=0; Ga < nirreps; Ga++) {
	      Gbc = Ga ^ Gijk;
	      Y[Ga] = dpd_block_matrix(virtpi[Ga],Fints.params->coltot[Gbc]);
	    }

	    for(i=0; i < occpi[Gi]; i++) {
	      I = occ_off[Gi] + i;
	      for(j=0; j < occpi[Gj]; j++) {
		J = occ_off[Gj] + j;
		for(k=0; k < occpi[Gk]; k++) {
		  K = occ_off[Gk] + k;

		  T3_UHF_AAA(WABC, VABC, 1, nirreps, I, Gi, J, Gj, K, Gk, &T2, &Fints, &Eints,
			     &T1, &Dints, &fIA, &fIJ, &fAB, occpi, occ_off, virtpi, vir_off, 0.0);

		  dpd_file2_mat_init(&T1);
		  dpd_file2_mat_rd(&T1);
		  dpd_file2_mat_init(&fIJ);
		  dpd_file2_mat_rd(&fIJ);
		  dpd_file2_mat_init(&fAB);
		  dpd_file2_mat_rd(&fAB);
		  dpd_file2_mat_init(&fIA);
		  dpd_file2_mat_rd(&fIA);
		  for(h=0; h < nirreps; h++) {
		    dpd_buf4_mat_irrep_init(&T2, h);
		    dpd_buf4_mat_irrep_rd(&T2, h);
		    dpd_buf4_mat_irrep_init(&Dints, h);
		    dpd_buf4_mat_irrep_rd(&Dints, h);
		    dpd_buf4_mat_irrep_init(&Eints, h);
		    dpd_buf4_mat_irrep_rd(&Eints, h);
		  }

		  ij = Eints.params->rowidx[I][J];
		  ji = Eints.params->rowidx[J][I];
		  jk = Eints.params->rowidx[J][K];
		  kj = Eints.params->rowidx[K][J];
		  ik = Eints.params->rowidx[I][K];
		  ki = Eints.params->rowidx[K][I];

		  dijk = 0.0;
		  if(fIJ.params->rowtot[Gi]) dijk += fIJ.matrix[Gi][i][i];
		  if(fIJ.params->rowtot[Gj]) dijk += fIJ.matrix[Gj][j][j];
		  if(fIJ.params->rowtot[Gk]) dijk += fIJ.matrix[Gk][k][k];

		  /**** Compute AAA part of (T) as a test ****/

		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {
		      A = Fints.params->colorb[Gab][ab][0];
		      Ga = Fints.params->rsym[A];
		      a = A - vir_off[Ga];
		      B = Fints.params->colorb[Gab][ab][1];
		      Gb = Fints.params->ssym[B];
		      b = B - vir_off[Gb];

		      for(c=0; c < virtpi[Gc]; c++) {
			C = vir_off[Gc] + c;
			denom = dijk;
			if(fAB.params->rowtot[Ga]) denom -= fAB.matrix[Ga][a][a];
			if(fAB.params->rowtot[Gb]) denom -= fAB.matrix[Gb][b][b];
			if(fAB.params->rowtot[Gc]) denom -= fAB.matrix[Gc][c][c];

			ET += WABC[Gab][ab][c] * (WABC[Gab][ab][c]+VABC[Gab][ab][c]) * denom;

		      } /* c */
		    } /* ab */
		  } /* Gab */

		  /**** Denominators and energy test complete ****/

		  /**** T3 --> S1 ****/

		  /* S_ia = 1/4 <jk||bc> t(c)_ijkabc */
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

			if(Gi==Ga && S1.params->rowtot[Gi] && S1.params->coltot[Gi])
			  S1.matrix[Gi][i][a] += 0.25 * WABC[Gab][ab][c] * Dints.matrix[Gjk][jk][bc];

		      } /* c */
		    } /* ab */
		  } /* Gab */

		  /**** T3 --> S1 Complete ****/

		  /**** Build Xijkabc = 2 Wijkabc + Vijkabc ****/

		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {
		      for(c=0; c < virtpi[Gc]; c++) {
			XABC[Gab][ab][c] = 2 * WABC[Gab][ab][c] + VABC[Gab][ab][c];
		      }
		    } /* ab */
		  } /* Gab */

		  /**** Xijkabc complete ****/

		  /**** T3 --> S2 ****/
		  /* S_JKDC <-- +1/2 <ID||AB> [2 W_IJKABC + V_IJKABC] */
		  /* S_JKCD <-- -1/2 <ID||AB> [2 W_IJKABC + V_IJKABC] */
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Gc = Gd ^ Gjk;
		    Gid = Gab = Gi ^ Gd;
		    nrows = virtpi[Gd];
		    ncols = virtpi[Gc];
		    nlinks = Fints.params->coltot[Gid];
		    if(nrows && ncols && nlinks) {
		      id = Fints.row_offset[Gid][I];
		      Fints.matrix[Gid] = dpd_block_matrix(nrows,Fints.params->coltot[Gid]);
		      dpd_buf4_mat_irrep_rd_block(&Fints, Gid, id, nrows);
		      Z = block_matrix(nrows, ncols);

		      C_DGEMM('n','n', nrows, ncols, nlinks, 0.5, Fints.matrix[Gid][0], nlinks,
			      XABC[Gab][0], ncols, 0.0, Z[0], ncols);

		      for(d=0; d < virtpi[Gd]; d++) {
			D = vir_off[Gd] + d;
			for(c=0; c < virtpi[Gc]; c++) {
			  C = vir_off[Gc] + c;
			  cd = S2.params->colidx[C][D];
			  dc = S2.params->colidx[D][C];
			  S2.matrix[Gjk][jk][dc] += Z[d][c];
			  S2.matrix[Gjk][jk][cd] -= Z[d][c];
			}
		      }
		      dpd_free_block(Fints.matrix[Gid], nrows, Fints.params->coltot[Gid]);
		      free_block(Z);
		    } /* if nrows && ncols && nlinks */
		  } /* Gd */

		  /* S_LIAB <-- +1/2 <JK||LC> [2 W_IJKABC + V_IJKABC] */
		  /* S_ILAB <-- -1/2 <JK||LC> [2 W_IJKABC + V_IJKABC] */
		  for(Gl=0; Gl < nirreps; Gl++) {
		    Gli = Gab = Gl ^ Gi;
		    Gc = Gab ^ Gijk;
		    lc = Eints.col_offset[Gjk][Gl];
		    nrows = occpi[Gl];
		    ncols = Fints.params->coltot[Gab];
		    nlinks = virtpi[Gc];
		    if(nrows && ncols && nlinks) {
		      Z = block_matrix(nrows, ncols);
		      C_DGEMM('n','t',nrows, ncols, nlinks, 0.5, &(Eints.matrix[Gjk][jk][lc]), nlinks,
			      XABC[Gab][0], nlinks, 0.0, Z[0], ncols);
		      for(l=0; l < occpi[Gl]; l++) {
			L = occ_off[Gl] + l;
			li = S2.params->rowidx[L][I];
			il = S2.params->rowidx[I][L];
			for(ab=0; ab < ncols; ab++) {
			  S2.matrix[Gli][li][ab] += Z[l][ab];
			  S2.matrix[Gli][il][ab] -= Z[l][ab];
			}
		      }
		      free_block(Z);
		    } /* nrows && ncols && nlinks */
		  } /* Gm */

		  /**** T3 --> S2 complete ****/

		  /**** T3 --> DAB ****/
		  for(Ga=0; Ga < nirreps; Ga++) {
		    Gb = Ga;
		    Gcd = Ga ^ Gijk;
		    for(Gc=0; Gc < nirreps; Gc++) {
		      Gd = Gc ^ Gcd;
		      Gac = Gbc = Ga ^ Gc;
		      for(a=0; a < virtpi[Ga]; a++) {
			A = vir_off[Ga] + a;
			for(b=0; b < virtpi[Gb]; b++) {
			  B = vir_off[Gb] + b;
			  for(c=0; c < virtpi[Gc]; c++) {
			    C = vir_off[Gc] + c;
			    ac = Fints.params->colidx[A][C];
			    bc = Fints.params->colidx[B][C];
			    for(d=0; d < virtpi[Gd]; d++) {
			      DAB.matrix[Ga][a][b] += (1.0/12.0) * WABC[Gac][ac][d] * (WABC[Gbc][bc][d] + VABC[Gbc][bc][d]);
			    } /* d */
			  } /* c */
			} /* b */
		      } /* a */
		    } /* Gc */
		  } /* Ga */

		  /**** T3 --> DAB complete ****/

		  /* T3 --> GIJAB ****/

		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    if(Gk == Gc) {
		      for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {
			for(c=0; c < virtpi[Gc]; c++) {
			  C = vir_off[Gc] + c;
			  if(T1.params->rowtot[Gk] && T1.params->coltot[Gk])
			    GIJAB.matrix[Gij][ij][ab] += WABC[Gab][ab][c] * T1.matrix[Gk][k][c];
			}
		      }
		    }

		  } /* Gab */

		  /**** T3 --> GIJAB complete ****/

		  /**** T3 --> GIJKA ****/
		  /**** Build Xijkabc = 2 * Wijkabc + Vijkabc ****/

		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {
		      A = Fints.params->colorb[Gab][ab][0];
		      Ga = Fints.params->rsym[A];
		      a = A - vir_off[Ga];
		      B = Fints.params->colorb[Gab][ab][1];
		      Gb = Fints.params->ssym[B];
		      b = B - vir_off[Gb];

		      for(c=0; c < virtpi[Gc]; c++) {
			C = vir_off[Gc] + c;
			bc = Fints.params->colidx[B][C];
			Y[Ga][a][bc] = 2 * WABC[Gab][ab][c] + VABC[Gab][ab][c];
		      }
		    } /* ab */
		  } /* Gab */

		  /**** Xijkabc complete ****/

		  /* G_IJLA = -1/2 t_KLBC Y_IJKABC */
		  for(Gl=0; Gl < nirreps; Gl++) {
		    Ga = Gl ^ Gij;
		    Gkl = Gbc = Gl ^ Gk;

		    nrows = occpi[Gl];
		    ncols = virtpi[Ga];
		    nlinks = T2.params->coltot[Gkl];
		    if(nrows && ncols && nlinks) {
		      kl = T2.row_offset[Gkl][K];
		      la = GIJKA.col_offset[Gij][Gl];
		      C_DGEMM('n','t', nrows, ncols, nlinks, -0.5, T2.matrix[Gkl][kl], nlinks,
			      Y[Ga][0], nlinks, 1.0, &(GIJKA.matrix[Gij][ij][la]), ncols);
		    }
		  } /* Gl */

		  /**** T3 --> GIJKA complete ****/

		  /* GIDAB = 1/2 t_JKCD X_IJKABC */
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Gab = Gid = Gi ^ Gd;
		    Gc = Gjk ^ Gd;

		    nrows = virtpi[Gd];
		    ncols = GIDAB.params->coltot[Gid];
		    nlinks = virtpi[Gc];
		    if(nrows && ncols && nlinks) {
		      id = GIDAB.row_offset[Gid][I];
		      cd = T2.col_offset[Gjk][Gc];
		      C_DGEMM('t','t', nrows, ncols, nlinks, 0.5, &(T2.matrix[Gjk][jk][cd]), nrows,
			      XABC[Gab][0], nlinks, 1.0, GIDAB.matrix[Gid][id], ncols);
		    }

		  }
		  /**** T3 --> GCIAB complete ****/

		  dpd_file2_mat_close(&T1);
		  dpd_file2_mat_close(&fIJ);
		  dpd_file2_mat_close(&fAB);
		  dpd_file2_mat_close(&fIA);
		  for(h=0; h < nirreps; h++) {
		    dpd_buf4_mat_irrep_close(&T2, h);
		    dpd_buf4_mat_irrep_close(&Dints, h);
		    dpd_buf4_mat_irrep_close(&Eints, h);
		  }

		} /* K */
	      } /* J */
	    } /* I */

	    for(Gab=0; Gab < nirreps; Gab++) {
	      Gc = Gab ^ Gijk;
	      dpd_free_block(WABC[Gab], Fints.params->coltot[Gab], virtpi[Gc]);
	      dpd_free_block(VABC[Gab], Fints.params->coltot[Gab], virtpi[Gc]);
	      dpd_free_block(XABC[Gab], Fints.params->coltot[Gab], virtpi[Gc]);
	    }
	    for(Ga=0; Ga < nirreps; Ga++) {
	      Gbc = Ga ^ Gijk;
	      dpd_free_block(Y[Ga], virtpi[Ga],Fints.params->coltot[Gbc]);
	    }

	  } /* Gk */
	} /* Gj */
      } /* Gi */

      ET *= (1.0/36.0);

      free(WABC);
      free(VABC);
      free(XABC);
      free(Y);

      dpd_file2_mat_wrt(&DAB);
      dpd_file2_mat_close(&DAB);
      dpd_file2_close(&DAB);

      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_wrt(&S2, h);
	dpd_buf4_mat_irrep_close(&S2, h);
      }
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_wrt(&GIJAB, h);
	dpd_buf4_mat_irrep_close(&GIJAB, h);
      }
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_wrt(&GIJKA, h);
	dpd_buf4_mat_irrep_close(&GIJKA, h);
      }
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_wrt(&GIDAB, h);
	dpd_buf4_mat_irrep_close(&GIDAB, h);
      }
      dpd_buf4_close(&S2);
      dpd_buf4_close(&GIJAB);
      dpd_buf4_close(&GIJKA);
      dpd_buf4_close(&GIDAB);

      dpd_file2_mat_wrt(&S1);
      dpd_file2_mat_close(&S1);
      dpd_file2_close(&S1);

      dpd_buf4_close(&T2);
      dpd_buf4_close(&Fints);
      dpd_buf4_close(&Eints);
      dpd_buf4_close(&Dints);

      dpd_file2_close(&T1);
      dpd_file2_close(&fIJ);
      dpd_file2_close(&fAB);
      dpd_file2_close(&fIA);


      /** T3 --> DIJ **/
      dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
      dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
      dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
      dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");

      dpdbuf4 EAAints;
      dpd_buf4_init(&EAAints, CC_EINTS, 0, 21, 0, 21, 2, 0, "E <AK||IJ> (AK, I>J)");
      dpdbuf4 FAAints;
      dpd_buf4_init(&FAAints, CC_FINTS, 0, 5, 20, 7, 20, 0, "F <BC||IA>");
      dpdbuf4 T2AA;
      dpd_buf4_init(&T2AA, CC_TAMPS, 0, 5, 0, 7, 2, 0, "tABIJ");
      dpdbuf4 DAAints;
      dpd_buf4_init(&DAAints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");

      dpd_file2_init(&DIJ, CC_OEI, 0, 0, 0, "DIJ");
      dpd_file2_mat_init(&DIJ);

      int Gabc;
      for (Ga=0; Ga < nirreps; ++Ga) {
        for (a=0; a<virtpi[Ga]; ++a) {
          A = vir_off[Ga] + a;
          for (Gb=0; Gb < nirreps; ++Gb) {
            for (b=0; b<virtpi[Gb]; ++b) {
              B = vir_off[Gb] + b;
              for (Gc=0; Gc < nirreps; ++Gc) {
                for (c=0; c < virtpi[Gc]; ++c) {
                  C = vir_off[Gc] + c;
                  Gabc = Ga ^ Gb ^ Gc;
                  //Allocate the memory for connected and disconnected triples
                  for (Gij=0; Gij < nirreps; ++Gij) {
                    Gk = Gij ^ Gabc;
                    WIJK[Gij] = dpd_block_matrix(T2AA.params->coltot[Gij], occpi[Gk]);
                    VIJK[Gij] = dpd_block_matrix(T2AA.params->coltot[Gij], occpi[Gk]);
                  }

                  T3_UHF_AAA_abc(WIJK, VIJK, 1, nirreps, A, Ga, B, Gb, C, Gc,
                      &T2AA, &FAAints, &EAAints, &T1, &DAAints, &fIA, &fIJ, &fAB,
                      occpi, occ_off, virtpi, vir_off, 0.0);

                  /**** T3 --> DIJ ****/
                  for(Gi=0; Gi < nirreps; Gi++) {
                    Gj = Gi;
                    Gkl = Gi ^ Gabc;
                    for(Gk=0; Gk < nirreps; Gk++) {
                      Gl = Gk ^ Gkl;
                      Gik = Gjk = Gi ^ Gk;
                      for(i=0; i < occpi[Gi]; i++) {
                        I = occ_off[Gi] + i;
                        for(j=0; j < occpi[Gj]; j++) {
                          J = occ_off[Gj] + j;
                          for(k=0; k < occpi[Gk]; k++) {
                            K = occ_off[Gk] + k;
                            ik = T2AA.params->colidx[I][K];
                            jk = T2AA.params->colidx[J][K];
                            for(l=0; l < occpi[Gl]; l++) {
                              DIJ.matrix[Gi][i][j] -= (1.0/12.0) * WIJK[Gik][ik][l] * (WIJK[Gjk][jk][l] + VIJK[Gjk][jk][l]);
                            } /* l */
                          } /* k */
                        } /* j */
                      } /* i */
                    } /* Gk */
                  } /* Gi */

                  for (Gij=0; Gij < nirreps; ++Gij) {
                    Gk = Gij ^ Gabc;
                    dpd_free_block(WIJK[Gij], T2AA.params->coltot[Gij], occpi[Gk]);
                    dpd_free_block(VIJK[Gij], T2AA.params->coltot[Gij], occpi[Gk]);
                  }

                }
              }
            }
          }
        }
      }

      dpd_file2_mat_wrt(&DIJ);
      dpd_file2_mat_close(&DIJ);
      dpd_file2_close(&DIJ);
      dpd_buf4_close(&EAAints);
      dpd_buf4_close(&FAAints);
      dpd_buf4_close(&T2AA);
      dpd_file2_close(&T1);
      dpd_file2_close(&fIJ);
      dpd_file2_close(&fIA);
      dpd_file2_close(&fAB);

      /** T3 --> DIJ complete **/

      return ET;

    } /* void T3_grad_UHF_AAA() */

  } /* namespace */
}
