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

    double T3_grad_UHF_BBB(void)
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
      dpdbuf4 T2, Fints, Eints, Dints, S2, Gijab, Gijka, Gidab;
      dpdfile2 fIJ, fAB, fIA, T1, S1, Dab, Dij;
      dpdfile2 fij, fab, fia;
      double ***WABC, ***VABC, ***XABC, ***Y;
      double **Z;
      double ***WIJK = (double ***) malloc(nirreps * sizeof(double **));
      double ***VIJK = (double ***) malloc(nirreps * sizeof(double **));

      nirreps = moinfo.nirreps;
      occpi = moinfo.boccpi; 
      virtpi = moinfo.bvirtpi;
      occ_off = moinfo.bocc_off;
      vir_off = moinfo.bvir_off;

      dpd_file2_init(&fIJ, CC_OEI, 0, 2, 2, "fij");
      dpd_file2_init(&fAB, CC_OEI, 0, 3, 3, "fab");
      dpd_file2_init(&fIA, CC_OEI, 0, 2, 3, "fia");
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");

      dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
      dpd_buf4_init(&Fints, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
      dpd_buf4_init(&Eints, CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
      dpd_buf4_init(&Dints, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");

      dpd_file2_init(&S1, CC_OEI, 0, 2, 3, "Sia");
      dpd_file2_mat_init(&S1); 
      dpd_buf4_init(&S2, CC_MISC, 0, 10, 15, 12, 17, 0, "Sijab");
      for(h=0; h < nirreps; h++)
	dpd_buf4_mat_irrep_init(&S2, h);

      dpd_file2_init(&Dab, CC_OEI, 0, 3, 3, "Dab");
      dpd_file2_mat_init(&Dab); 

      dpd_buf4_init(&Gijab, CC_GAMMA, 0, 10, 15, 12, 17, 0, "Gijab");
      for(h=0; h < nirreps; h++)
	dpd_buf4_mat_irrep_init(&Gijab, h);

      dpd_buf4_init(&Gijka, CC_GAMMA, 0, 10, 30, 12, 30, 0, "Gijka");
      for(h=0; h < nirreps; h++)
	dpd_buf4_mat_irrep_init(&Gijka, h);

      dpd_buf4_init(&Gidab, CC_GAMMA, 0, 30, 15, 30, 17, 0, "Gidab");
      for(h=0; h < nirreps; h++)
	dpd_buf4_mat_irrep_init(&Gidab, h);

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

		  /**** Compute BBB part of (T) as a test ****/

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

		  /**** T3 contributions to S1 ****/

		  /* S_ijab = 1/4 <jk||bc> t(c)_ijkabc */
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
		    XABC[Gab] = dpd_block_matrix(Fints.params->coltot[Gab], virtpi[Gc]);
		    for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {
		      for(c=0; c < virtpi[Gc]; c++) {
			XABC[Gab][ab][c] = 2 * WABC[Gab][ab][c] + VABC[Gab][ab][c];
		      }
		    } /* ab */
		  } /* Gab */

		  /**** Xijkabc complete ****/

		  /**** T3 --> S2 ****/
		  /* S_jkdc <-- +1/2 <id||ab> [2 W_ijkabc + V_ijkabc] */
		  /* S_jkcd <-- -1/2 <id||ab> [2 W_ijkabc + V_ijkabc] */
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

		  /* S_liab <-- +1/2 <jk||lc> [2 W_ijkabc + V_ijkabc] */
		  /* S_ilab <-- -1/2 <jk||lc> [2 W_ijkabc + V_ijkabc] */
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
			      Dab.matrix[Ga][a][b] += (1.0/12.0) * WABC[Gac][ac][d] * (WABC[Gbc][bc][d] + VABC[Gbc][bc][d]);
			    } /* d */
			  } /* c */
			} /* b */
		      } /* a */
		    } /* Gc */
		  } /* Ga */

		  /**** T3 --> DAB complete ****/

		  /**** T3 --> Gijab ****/
		  /** This can be simplified */
		  for(Gab=0; Gab < nirreps; Gab++) {
		    Gc = Gab ^ Gijk;
		    if(Gk == Gc) {
		      for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {
			for(c=0; c < virtpi[Gc]; c++) {
			  C = vir_off[Gc] + c;
			  if(T1.params->rowtot[Gk] && T1.params->coltot[Gk])
			    Gijab.matrix[Gij][ij][ab] += WABC[Gab][ab][c] * T1.matrix[Gk][k][c];
			}
		      }
		    }

		  } /* Gab */
		  /**** T3 --> Gijab complete ****/

		  /**** T3 --> Gijka ****/
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

		  /* G_ijla = -1/2 t_klbc Y_ijkabc */
		  for(Gl=0; Gl < nirreps; Gl++) {
		    Ga = Gl ^ Gij;
		    Gkl = Gbc = Gl ^ Gk;

		    nrows = occpi[Gl];
		    ncols = virtpi[Ga];
		    nlinks = T2.params->coltot[Gkl];
		    if(nrows && ncols && nlinks) {
		      kl = T2.row_offset[Gkl][K];
		      la = Gijka.col_offset[Gij][Gl];
		      C_DGEMM('n','t', nrows, ncols, nlinks, -0.5, T2.matrix[Gkl][kl], nlinks,
			      Y[Ga][0], nlinks, 1.0, &(Gijka.matrix[Gij][ij][la]), ncols);
		    }
		  } /* Gl */


		  /**** T3 --> Gijka complete ****/

		  /* Gidab = 1/2 t_jkcd X_ijkabc */
		  for(Gd=0; Gd < nirreps; Gd++) {
		    Gab = Gid = Gi ^ Gd;
		    Gc = Gjk ^ Gd;

		    nrows = virtpi[Gd];
		    ncols = Gidab.params->coltot[Gid];
		    nlinks = virtpi[Gc];
		    if(nrows && ncols && nlinks) {
		      id = Gidab.row_offset[Gid][I];
		      cd = T2.col_offset[Gjk][Gc];
		      C_DGEMM('t','t', nrows, ncols, nlinks, 0.5, &(T2.matrix[Gjk][jk][cd]), nrows,
			      XABC[Gab][0], nlinks, 1.0, Gidab.matrix[Gid][id], ncols);
		    }

		  }
		  /**** T3 --> Gciab complete ****/

		  dpd_file2_mat_close(&T1);
		  dpd_file2_mat_close(&fIJ);
		  dpd_file2_mat_close(&fAB);
		  dpd_file2_mat_close(&fIA);
		  for(h=0; h < nirreps; h++) {
		    dpd_buf4_mat_irrep_close(&T2, h);
		    dpd_buf4_mat_irrep_close(&Eints, h);
		    dpd_buf4_mat_irrep_close(&Dints, h);
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

      dpd_file2_mat_wrt(&Dab);
      dpd_file2_mat_close(&Dab);
      dpd_file2_close(&Dab);

      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_wrt(&S2, h);
	dpd_buf4_mat_irrep_close(&S2, h);
      }
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_wrt(&Gijab, h);
	dpd_buf4_mat_irrep_close(&Gijab, h);
      }
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_wrt(&Gijka, h);
	dpd_buf4_mat_irrep_close(&Gijka, h);
      }
      for(h=0; h < nirreps; h++) {
	dpd_buf4_mat_irrep_wrt(&Gidab, h);
	dpd_buf4_mat_irrep_close(&Gidab, h);
      }
      dpd_buf4_close(&Gijab);
      dpd_buf4_close(&Gijka);
      dpd_buf4_close(&Gidab);
      dpd_buf4_close(&S2);

      dpd_buf4_close(&T2);
      dpd_buf4_close(&Fints);
      dpd_buf4_close(&Eints);
      dpd_buf4_close(&Dints);

      dpd_file2_mat_wrt(&S1);
      dpd_file2_mat_close(&S1);
      dpd_file2_close(&S1);

      dpd_file2_close(&T1);
      dpd_file2_close(&fIJ);
      dpd_file2_close(&fAB);
      dpd_file2_close(&fIA);

      /*** T3 --> Dij ***/

      dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
      dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");
      dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
      dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");

      dpdbuf4 EBBints;
      dpd_buf4_init(&EBBints, CC_EINTS, 0, 31, 10, 31, 12, 0, "E <ak||ij> (ak, i>j)");
      dpdbuf4 FBBints;
      dpd_buf4_init(&FBBints, CC_FINTS, 0, 15, 30, 17, 30, 0, "F <bc||ia>");
      dpdbuf4 T2BB;
      dpd_buf4_init(&T2BB, CC_TAMPS, 0, 15, 10, 17, 12, 0, "tabij");
      dpdbuf4 DBBints;
      dpd_buf4_init(&DBBints, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");

      dpd_file2_init(&Dij, CC_OEI, 0, 2, 2, "Dij");
      dpd_file2_mat_init(&Dij);

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
                    WIJK[Gij] = dpd_block_matrix(T2BB.params->coltot[Gij], occpi[Gk]);
                    VIJK[Gij] = dpd_block_matrix(T2BB.params->coltot[Gij], occpi[Gk]);
                  }

                  T3_UHF_AAA_abc(WIJK, VIJK, 1, nirreps, A, Ga, B, Gb, C, Gc,
                      &T2BB, &FBBints, &EBBints, &T1, &DBBints, &fia, &fij, &fab,
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
                            ik = T2BB.params->colidx[I][K];
                            jk = T2BB.params->colidx[J][K];
                            for(l=0; l < occpi[Gl]; l++) {
                              Dij.matrix[Gi][i][j] -= (1.0/12.0) * WIJK[Gik][ik][l] * (WIJK[Gjk][jk][l] + VIJK[Gjk][jk][l]);
                            } /* l */
                          } /* k */
                        } /* j */
                      } /* i */
                    } /* Gk */
                  } /* Gi */

                  //Deallocate the memory for connected and disconnected triples
                  for (Gij=0; Gij < nirreps; ++Gij) {
                    Gk = Gij ^ Gabc;
                    dpd_free_block(WIJK[Gij], T2BB.params->coltot[Gij], occpi[Gk]);
                    dpd_free_block(VIJK[Gij], T2BB.params->coltot[Gij], occpi[Gk]);
                  }

                } // c
              } // Gc
            } // b
          } // Gb
        } // a
      } // Ga

      dpd_file2_mat_wrt(&Dij);
      dpd_file2_mat_close(&Dij);
      dpd_file2_close(&Dij);

      dpd_buf4_close(&EBBints);
      dpd_buf4_close(&FBBints);
      dpd_buf4_close(&DBBints);
      dpd_buf4_close(&T2BB);
      dpd_file2_close(&T1);
      dpd_file2_close(&fij);
      dpd_file2_close(&fia);
      dpd_file2_close(&fab);

      return ET;

    } /* void T3_grad_UHF_BBB() */

  } /* namespace */
}
