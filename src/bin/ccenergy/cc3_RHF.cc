/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void cc3_RHF_obsolete(void)
{
  int h, nirreps;
  int *occpi, *virtpi, *occ_off, *vir_off;
  int Gi, Gj, Gk, Gl, Ga, Gb, Gc, Gd;
  int i, j, k, l, a, b, c, d;
  int I, J, K, L, A, B, C, D;
  int kj, jk, ji, ij, ik, ki;
  int Gkj, Gjk, Gji, Gij, Gik, Gki;
  int Gijk;
  int ab, ba, ac, ca, bc, cb;
  int Gab, Gba, Gac, Gca, Gbc, Gcb;
  int id, jd, kd, ad, bd, cd;
  int il, jl, kl, la, lb, lc, li, lk;
  int da, di, dj, dk;
  int Gad, Gdi, Gdj, Gdk, Glc, Gli, Glk;
  double value, F_val, t_val, E_val;
  double dijk, denom;
  double value_ia, value_ka, denom_ia, denom_ka;
  dpdfile2 fIJ, fIJ2, fAB, fAB2, t1, t1new, Fme;
  dpdbuf4 T2, T2new, F, E, Dints, Wamef, Wmnie;
  double t2norm, t3norm;
  double ***T3, ***W3;
  int nv;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off;
  vir_off = moinfo.vir_off;

  /* these are sent to T3 function */
  dpd_file2_init(&fIJ2, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB2, CC_OEI, 0, 1, 1, "fAB");

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_mat_rd(&fAB);

  dpd_file2_init(&t1new, CC_OEI, 0, 0, 1, "CC3 tIA");  /* T3->T1 increment */
  dpd_file2_mat_init(&t1new);
  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_mat_init(&Fme);
  dpd_file2_mat_rd(&Fme);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

  dpd_buf4_init(&T2new, CC_MISC, 0, 0, 5, 0, 5, 0, "CC3 tIjAb");
  dpd_buf4_scm(&T2new, 0);
  dpd_buf4_init(&F, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
  dpd_buf4_init(&E, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
  dpd_buf4_init(&Wamef, CC3_HET1, 0, 11, 5, 11, 5, 0, "CC3 WAmEf (Am,Ef)");
  dpd_buf4_init(&Wmnie, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");
  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  for(h=0; h < nirreps; h++) {

    dpd_buf4_mat_irrep_init(&Dints, h);
    dpd_buf4_mat_irrep_rd(&Dints, h);

    dpd_buf4_mat_irrep_init(&Wamef, h);
    dpd_buf4_mat_irrep_rd(&Wamef, h);

    dpd_buf4_mat_irrep_init(&Wmnie, h);
    dpd_buf4_mat_irrep_rd(&Wmnie, h);

    dpd_buf4_mat_irrep_init(&T2new, h);
  }

  for(h=0,nv=0; h < nirreps; h++) nv += virtpi[h];

  W3 = (double ***) malloc(nirreps * sizeof(double **));
  T3 = init_3d_array(nv, nv, nv);

  t3norm = 0.0;
  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gi ^ Gj;
      for(Gk=0; Gk < nirreps; Gk++) {
        Gijk = Gi ^ Gj ^ Gk;

        /* allocate memory for all irrep blocks of (ab,c) */
        for(Gab=0; Gab < nirreps; Gab++) {
          Gc = Gab ^ Gijk;
          W3[Gab] = dpd_block_matrix(F.params->coltot[Gab], virtpi[Gc]);
        }

        for(i=0; i < occpi[Gi]; i++) {
          I = occ_off[Gi] + i;
          for(j=0; j < occpi[Gj]; j++) {
            J = occ_off[Gj] + j;
            for(k=0; k < occpi[Gk]; k++) {
              K = occ_off[Gk] + k;

	      Gkj = Gjk = Gk ^ Gj;
	      Gji = Gij = Gi ^ Gj;
	      Gik = Gki = Gi ^ Gk;
	      Gijk = Gi ^ Gj ^ Gk;

	      ij = T2.params->rowidx[I][J];
	      ji = T2.params->rowidx[J][I];
	      ik = T2.params->rowidx[I][K];
	      ki = T2.params->rowidx[K][I];
	      jk = T2.params->rowidx[J][K];
	      kj = T2.params->rowidx[K][J];

             /* T3_RHF(W3, nirreps, I, Gi, J, Gj, K, Gk, &T2, &F, &E,
                     &fIJ2, &fAB2, occpi, occ_off, virtpi, vir_off, 0.0); */


              /* sort (ab,c) into T3[A][B][C] */
              for(Gab=0; Gab < nirreps; Gab++) {
                Gc = Gab ^ Gijk;
                for (ab=0; ab < F.params->coltot[Gab]; ++ab) {
                  A = F.params->colorb[Gab][ab][0];
                  B = F.params->colorb[Gab][ab][1];
                  Ga = F.params->rsym[A];
                  Gb = Ga ^ Gab;
                  /* a = A - vir_off[Ga];
                  b = B - vir_off[Gb]; */
                  for (c=0; c<virtpi[Gc]; ++c) {
                    C = c + vir_off[Gc];
                    T3[A][B][C] = W3[Gab][ab][c];
                  }
                }
              }

              /* if (!I && !J && !K) {
                fprintf(outfile,"T3[0][0][0][a][b][c] fast code\n");
                for (a=0;a<nv;++a)
                  for (b=0;b<nv;++b)
                    for (c=0;c<nv;++c)
                      if ( fabs(T3[a][b][c]) > 1.0E-10 )
                      fprintf(outfile,"A: %d B: %d C: %d -> %15.10lf\n", a,b,c,T3[a][b][c]);
              } */


              /*** T3 --> T1 Contributions ***/

              for(Ga=0; Ga < nirreps; Ga++) {
                for(a=0; a < virtpi[Ga]; a++) {
                  A = vir_off[Ga] + a;

                  value_ia = 0.0;
                  value_ka = 0.0;
                  for(Gb=0; Gb < nirreps; Gb++) {
                    for(b=0; b < virtpi[Gb]; b++) {
                      B = vir_off[Gb] + b;

                      Gc = Gijk ^ Ga ^ Gb;
                      Gbc = Gb ^ Gc;

                      for(c=0; c < virtpi[Gc]; c++) {
                        C = vir_off[Gc] + c;

                        bc = Dints.params->colidx[B][C];

                        if(Gi == Ga && Gjk == Gbc) {

                          if(Dints.params->rowtot[Gjk] && Dints.params->coltot[Gjk])
                            value_ia += T3[A][B][C] * Dints.matrix[Gjk][jk][bc];
                        }

                        if(Gk == Ga && Gji == Gbc) {

                          if(Dints.params->rowtot[Gji] && Dints.params->coltot[Gji])
                            value_ka -= T3[A][B][C] * Dints.matrix[Gji][ji][bc];
                        }

                      } /* c */
                    } /* b */
                  } /* Gb */

                  denom_ia = denom_ka = 0.0;
                  if(fIJ.params->rowtot[Gi])
                    denom_ia += fIJ.matrix[Gi][i][i];

                  if(fIJ.params->rowtot[Gk])
                    denom_ka += fIJ.matrix[Gk][k][k];

                  if(fAB.params->rowtot[Ga]) {
                    denom_ia -= fAB.matrix[Ga][a][a];
                    denom_ka -= fAB.matrix[Ga][a][a];
                  }
                  value_ia /= denom_ia;
                  value_ka /= denom_ka;

                  if(t1new.params->rowtot[Gi] && t1new.params->coltot[Gi])
                    t1new.matrix[Gi][i][a] += value_ia;
                  if(t1new.params->rowtot[Gk] && t1new.params->coltot[Gk])
                    t1new.matrix[Gk][k][a] += value_ka;

                } /* a */
              } /* Ga */

	      /*** T3 --> T2 Contributions ***/

	      for(Ga=0; Ga < nirreps; Ga++) {
	       	for(Gb=0; Gb < nirreps; Gb++) {
		  Gc = Gijk ^ Ga ^ Gb;
		  Gab = Ga ^ Gb;

		  for(a=0; a < virtpi[Ga]; a++) {
		    A = vir_off[Ga] + a;
		    for(b=0; b < virtpi[Gb]; b++) {
		      B = vir_off[Gb] + b;

                      ab = T2new.params->colidx[A][B];
                      ba = T2new.params->colidx[B][A];

		      if(Gij == Gab && Gk == Gc) { 

			denom = 0.0;
		       	if(fIJ.params->rowtot[Gi])
			  denom += fIJ.matrix[Gi][i][i];
		       	if(fIJ.params->rowtot[Gj])
			  denom += fIJ.matrix[Gj][j][j];
		       	if(fAB.params->rowtot[Ga])
			  denom -= fAB.matrix[Ga][a][a];
		       	if(fAB.params->rowtot[Gb])
			  denom -= fAB.matrix[Gb][b][b];

			value = 0.0;
		       	for(c=0; c < virtpi[Gc]; c++) {
			  C = vir_off[Gc] + c;
			  if(Fme.params->rowtot[Gk] && Fme.params->coltot[Gc])
			    value += T3[A][B][C] * Fme.matrix[Gk][k][c];
			}
			value /= denom;

			if(T2new.params->rowtot[Gij] && T2new.params->coltot[Gab]) {
			  T2new.matrix[Gij][ij][ab] += value;
			  T2new.matrix[Gij][ji][ba] += value;
			}
		      }

                      if(Gjk == Gab && Gi == Gc) {

                        denom = 0.0;
                        if(fIJ.params->rowtot[Gk])
                          denom += fIJ.matrix[Gk][k][k];
                        if(fIJ.params->rowtot[Gj])
                          denom += fIJ.matrix[Gj][j][j];
                        if(fAB.params->rowtot[Ga])
                          denom -= fAB.matrix[Ga][a][a];
                        if(fAB.params->rowtot[Gb])
                          denom -= fAB.matrix[Gb][b][b];

                        value = 0.0;
                        for(c=0; c < virtpi[Gc]; c++) {
			  C = vir_off[Gc] + c;
                          if(Fme.params->rowtot[Gi] && Fme.params->coltot[Gc])
                            value -= T3[A][B][C] * Fme.matrix[Gi][i][c];
			}
                        value /= denom;

                        if(T2new.params->rowtot[Gjk] && T2new.params->coltot[Gab]) {
                          T2new.matrix[Gjk][kj][ab] += value;
                          T2new.matrix[Gjk][jk][ba] += value;
                        }
                      }

		    } /* b */
		  } /* a */
	 	} /* Gb */
	      } /* Ga */

	      for(Gd=0; Gd < nirreps; Gd++) {
		Gdi = Gd ^ Gi;
		Gdj = Gd ^ Gj;
		Gdk = Gd ^ Gk;
	       	for(Ga=0; Ga < nirreps; Ga++) {
		  Gad = Ga ^ Gd;

		  for(d=0; d < virtpi[Gd]; d++) {
		    D = vir_off[Gd] + d;

		    di = Wamef.params->rowidx[D][I];
		    dj = Wamef.params->rowidx[D][J];
		    dk = Wamef.params->rowidx[D][K];

		    for(a=0; a < virtpi[Ga]; a++) {
		      A = vir_off[Ga] + a;

		      ad = T2new.params->colidx[A][D];
		      da = T2new.params->colidx[D][A];

		      if(Gij == Gad) {
			value = 0.0;
			for(Gb=0; Gb < nirreps; Gb++) {
			  Gc = Gijk ^ Ga ^ Gb;
			  Gbc = Gb ^ Gc;
			  if(Gdk == Gbc) {
			    for(b=0; b < virtpi[Gb]; b++) {
			      B = vir_off[Gb] + b;
			      for(c=0; c < virtpi[Gc]; c++) {
			       	C = vir_off[Gc] + c;

				bc = Wamef.params->colidx[B][C];

				if(Wamef.params->rowtot[Gdk] && Wamef.params->coltot[Gdk])
				  value += 2.0 * T3[A][B][C] * Wamef.matrix[Gdk][dk][bc];
			      }
			    }
			  }
			}

			denom = 0.0;
			if(fIJ.params->rowtot[Gi])
			  denom += fIJ.matrix[Gi][i][i];
			if(fIJ.params->rowtot[Gj])
			  denom += fIJ.matrix[Gj][j][j];
			if(fAB.params->rowtot[Ga])
			  denom -= fAB.matrix[Ga][a][a];
			if(fAB.params->rowtot[Gd])
			  denom -= fAB.matrix[Gd][d][d];
			value /= denom;

			if(T2new.params->rowtot[Gij] && T2new.params->coltot[Gij]) {
			  T2new.matrix[Gij][ij][ad] += value;
			  T2new.matrix[Gij][ji][da] += value;
			}
		      }

		      if(Gjk == Gad) {
                        value = 0.0; 
                        for(Gb=0; Gb < nirreps; Gb++) {
                          Gc = Gijk ^ Ga ^ Gb;
                          Gbc = Gb ^ Gc; 
                          if(Gdi == Gbc) {
                            for(b=0; b < virtpi[Gb]; b++) {
                              B = vir_off[Gb] + b;
                              for(c=0; c < virtpi[Gc]; c++) {
                                C = vir_off[Gc] + c;

				bc = Wamef.params->colidx[B][C];

                                if(Wamef.params->rowtot[Gdi] && Wamef.params->coltot[Gdi])
                                  value -= T3[A][B][C] * Wamef.matrix[Gdi][di][bc];
                              }   
                            } 
                          } 
                        } 
                        
                        denom = 0.0;
                        if(fIJ.params->rowtot[Gk])
                          denom += fIJ.matrix[Gk][k][k];
                        if(fIJ.params->rowtot[Gj])
                          denom += fIJ.matrix[Gj][j][j];
                        if(fAB.params->rowtot[Ga])
                          denom -= fAB.matrix[Ga][a][a];
                        if(fAB.params->rowtot[Gd])
                          denom -= fAB.matrix[Gd][d][d];
                        value /= denom;
                        
                        if(T2new.params->rowtot[Gjk] && T2new.params->coltot[Gjk]) {
                          T2new.matrix[Gjk][kj][ad] += value;
                          T2new.matrix[Gjk][jk][da] += value;
                        } 

		      }

		      if(Gik == Gad) {
                        value = 0.0;
                        for(Gb=0; Gb < nirreps; Gb++) {
                          Gc = Gijk ^ Ga ^ Gb;
                          Gbc = Gb ^ Gc;
                          if(Gdj == Gbc) {
                            for(b=0; b < virtpi[Gb]; b++) {
                              B = vir_off[Gb] + b;
                              for(c=0; c < virtpi[Gc]; c++) {
                                C = vir_off[Gc] + c;

				bc = Wamef.params->colidx[B][C];

                                if(Wamef.params->rowtot[Gdj] && Wamef.params->coltot[Gdj])
                                  value -= T3[A][B][C] * Wamef.matrix[Gdj][dj][bc];
                              }
                            }
                          }
                        }

                        denom = 0.0;
                        if(fIJ.params->rowtot[Gk])
                          denom += fIJ.matrix[Gk][k][k];
                        if(fIJ.params->rowtot[Gi])
                          denom += fIJ.matrix[Gi][i][i];
                        if(fAB.params->rowtot[Ga])
                          denom -= fAB.matrix[Ga][a][a];
                        if(fAB.params->rowtot[Gd])
                          denom -= fAB.matrix[Gd][d][d];
                        value /= denom;

                        if(T2new.params->rowtot[Gik] && T2new.params->coltot[Gik]) {
                          T2new.matrix[Gki][ki][da] += value;
                          T2new.matrix[Gki][ik][ad] += value;
                        }
		      }

		    } /* a */
		  } /* d */
		} /* Ga */
	      } /* Gd */

	      for(Gl=0; Gl < nirreps; Gl++) {
		Gli = Gl ^ Gi;
		Glk = Gl ^ Gk;

		for(Ga=0; Ga < nirreps; Ga++) {
		  for(Gb=0; Gb < nirreps; Gb++) {
		    Gab = Ga ^ Gb;

		    if(Gli == Gab) {

		      for(l=0; l < occpi[Gl]; l++) {
			L = occ_off[Gl] + l;

			li = T2new.params->rowidx[L][I];
			il = T2new.params->rowidx[I][L];

			for(a=0; a < virtpi[Ga]; a++) {
			  A = vir_off[Ga] + a;
			  for(b=0; b < virtpi[Gb]; b++) {
			    B = vir_off[Gb] + b;

			    ab = T2new.params->colidx[A][B];
			    ba = T2new.params->colidx[B][A];

			    denom = 0.0;
			    if(fIJ.params->rowtot[Gl])
			      denom += fIJ.matrix[Gl][l][l];
			    if(fIJ.params->rowtot[Gi])
			      denom += fIJ.matrix[Gi][i][i];
			    if(fAB.params->rowtot[Ga])
			      denom -= fAB.matrix[Ga][a][a];
			    if(fAB.params->rowtot[Gb])
			      denom -= fAB.matrix[Gb][b][b];

			    value = 0.0;
			    for(Gc=0; Gc < nirreps; Gc++) {
			      Glc = Gl ^ Gc;

			      if(Gjk == Glc) {

				for(c=0; c < virtpi[Gc]; c++) {
				  C = vir_off[Gc] + c;
				  lc = Wmnie.params->colidx[L][C];

				  if(Wmnie.params->rowtot[Gjk] && Wmnie.params->coltot[Gjk])
				    value -= T3[A][B][C] * Wmnie.matrix[Gjk][jk][lc];
				}
			      }

			    } /* Gc */

			    value *= 2.0;
			    value /= denom;

			    if(T2new.params->rowtot[Gli] && T2new.params->coltot[Gli]) {
			      T2new.matrix[Gli][li][ba] += value;
			      T2new.matrix[Gli][il][ab] += value;
			    }

			  } /* b */
			} /* a */
		      } /* l */

		    } /* if(Gli == Gab) */

		    if(Glk == Gab) {

		      for(l=0; l < occpi[Gl]; l++) {
			L = occ_off[Gl] + l;

			lk = T2new.params->rowidx[L][K];
			kl = T2new.params->rowidx[K][L];

			for(a=0; a < virtpi[Ga]; a++) {
			  A = vir_off[Ga] + a;
			  for(b=0; b < virtpi[Gb]; b++) {
			    B = vir_off[Gb] + b;

			    ab = T2new.params->colidx[A][B];
			    ba = T2new.params->colidx[B][A];

			    denom = 0.0;
			    if(fIJ.params->rowtot[Gl])
			      denom += fIJ.matrix[Gl][l][l];
			    if(fIJ.params->rowtot[Gk])
			      denom += fIJ.matrix[Gk][k][k];
			    if(fAB.params->rowtot[Ga])
			      denom -= fAB.matrix[Ga][a][a];
			    if(fAB.params->rowtot[Gb])
			      denom -= fAB.matrix[Gb][b][b];

			    value = 0.0;
			    for(Gc=0; Gc < nirreps; Gc++) {
			      Glc = Gl ^ Gc;

			      if(Gji == Glc) {

				for(c=0; c < virtpi[Gc]; c++) {
				  C = vir_off[Gc] + c;
				  lc = Wmnie.params->colidx[L][C];

				  if(Wmnie.params->rowtot[Gji] && Wmnie.params->coltot[Gji])
				    value += T3[A][B][C] * Wmnie.matrix[Gji][ji][lc];
				}
			      }

			    } /* Gc */

			    value /= denom;

			    if(T2new.params->rowtot[Glk] && T2new.params->coltot[Glk]) {
			      T2new.matrix[Glk][lk][ba] += value;
			      T2new.matrix[Glk][kl][ab] += value;
			    }

			  } /* b */
			} /* a */
		      } /* l */

		    } /* if(Glk==Gab) */

		    if(Gli == Gab) {

		      for(l=0; l < occpi[Gl]; l++) {
			L = occ_off[Gl] + l;

			li = T2new.params->rowidx[L][I];
			il = T2new.params->rowidx[I][L];

			for(a=0; a < virtpi[Ga]; a++) {
			  A = vir_off[Ga] + a;
			  for(b=0; b < virtpi[Gb]; b++) {
			    B = vir_off[Gb] + b;

			    ab = T2new.params->colidx[A][B];
			    ba = T2new.params->colidx[B][A];

			    denom = 0.0;
			    if(fIJ.params->rowtot[Gl])
			      denom += fIJ.matrix[Gl][l][l];
			    if(fIJ.params->rowtot[Gi])
			      denom += fIJ.matrix[Gi][i][i];
			    if(fAB.params->rowtot[Ga])
			      denom -= fAB.matrix[Ga][a][a];
			    if(fAB.params->rowtot[Gb])
			      denom -= fAB.matrix[Gb][b][b];

			    value = 0.0;
			    for(Gc=0; Gc < nirreps; Gc++) {
			      Glc = Gl ^ Gc;

			      if(Gjk == Glc) {

				for(c=0; c < virtpi[Gc]; c++) {
				  C = vir_off[Gc] + c;
				  lc = Wmnie.params->colidx[L][C];

				  if(Wmnie.params->rowtot[Gjk] && Wmnie.params->coltot[Gjk])
				    value += T3[A][B][C] * Wmnie.matrix[Gjk][kj][lc];
				}
			      }

			    } /* Gc */

			    value /= denom;

			    if(T2new.params->rowtot[Gli] && T2new.params->coltot[Gli]) {
			      T2new.matrix[Gli][li][ba] += value;
			      T2new.matrix[Gli][il][ab] += value;
			    }

			  } /* b */
			} /* a */
		      } /* l */

		    } /* if(Gli==Gab) */

		  } /* Gb */
		} /* Ga */
	      } /* Gl */

	    } /* k */
	  } /* j */
	} /* i */

        for(Gab=0; Gab < nirreps; Gab++) {
          /* This will need to change for non-totally-symmetric cases */
          Gc = Gab ^ Gijk;
          dpd_free_block(W3[Gab], F.params->coltot[Gab], virtpi[Gc]);
        }
      } /* Gk */
    } /* Gj */
  } /* Gi */

  /*
  fprintf(outfile, "t3norm = %20.10f\n", sqrt(t3norm));
  */
  free_3d_array(T3, nv, nv);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&Dints, h);
    dpd_buf4_mat_irrep_close(&Wamef, h);
    dpd_buf4_mat_irrep_close(&Wmnie, h);

    dpd_buf4_mat_irrep_wrt(&T2new, h);
    dpd_buf4_mat_irrep_close(&T2new, h);
  }

  dpd_buf4_close(&T2);
  dpd_buf4_close(&F);
  dpd_buf4_close(&E);
  dpd_buf4_close(&Dints);
  dpd_buf4_close(&Wamef);
  dpd_buf4_close(&Wmnie);

  dpd_file2_mat_wrt(&t1new);
  dpd_file2_mat_close(&t1new);
  /*
  fprintf(outfile, "triples t1norm = %20.10f\n", sqrt(dpd_file2_dot_self(&t1new)));
  fprintf(outfile, "triples t2norm = %20.10f\n", sqrt(dpd_buf4_dot_self(&T2new)));
  */

  /* Update the amplitudes */
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
  dpd_buf4_axpy(&T2new, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2new);

  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "New tIA");
  dpd_file2_axpy(&t1new, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1new);

  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fAB);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);

  dpd_file2_close(&fIJ2);
  dpd_file2_close(&fAB2);
}
}} // namespace psi::ccenergy
