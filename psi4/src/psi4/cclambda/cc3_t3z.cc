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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
/* cc3_t3z(): Builds the following intermediates, which are required
** for the T3-->L1 contributions in the CC3 model:
**
** Zifln = +1/2 t3(lmn,def) * <im||de>
**
** Zdfan = -1/2 t3(lmn,def) * <lm||ae>
**
** Translating indices in each expression to use t3(ijk,abc):
**
** Zmcik = +1/2 t3(ijk,abc) * <mj||ab>
**
** Zacek = -1/2 t3(ijk,abc) * <ij||eb>
**
** The current version is coded ROHF-like, but works only for RHF.
**
** -TDC, 7/04
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void cc3_t3z_RHF_AAA(void);
void cc3_t3z_RHF_AAB(void);

void cc3_t3z(void)
{
  if(params.ref == 0) { /** RHF **/
    cc3_t3z_RHF_AAA();
    cc3_t3z_RHF_AAB();
  }
  else if(params.ref == 2) { /** UHF **/
    /* TBD */
  }
}

void cc3_t3z_RHF_AAA(void)
{
  int h, nirreps;
  int *occ_off, *occpi;
  int *vir_off, *virtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc, Gab;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int ab;
  double ***W1, ***W2;
  dpdbuf4 T2, E, F;
  dpdfile2 fIJ, fAB;
  dpdbuf4 ZIFLN, ZDFAN;
  dpdbuf4 Dints;
  int Gm, Gmj, Gmc, mj, mc;
  int m, M, ik;
  int nrows, ncols, nlinks;
  double *Z;
  int Gij, ij, Gca, ca, Ge, Gke, ke;
  int EE, e, eb;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  occ_off = moinfo.occ_off;
  virtpi = moinfo.virtpi;
  vir_off = moinfo.vir_off;

  global_dpd_->buf4_init(&ZIFLN, PSIF_CC3_MISC, 0, 10, 0, 10, 0, 0, "CC3 ZIFLN");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&ZIFLN, h);
  }

  global_dpd_->buf4_init(&ZDFAN, PSIF_CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZDFAN (NA,FD)");
  global_dpd_->buf4_scm(&ZDFAN, 0.0);

  global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Dints, h);
    global_dpd_->buf4_mat_irrep_rd(&Dints, h);
  }

  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  global_dpd_->buf4_init(&F, PSIF_CC3_HET1, 0, 10, 5, 10, 7, 0, "CC3 WABEI (IE,B>A)");
  global_dpd_->buf4_init(&E, PSIF_CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 WMBIJ (I>J,MB)");

  /* target T3 amplitudes go in here */
  W1 = (double ***) malloc(nirreps * sizeof(double **));
  W2 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gi ^ Gj;
      for(Gk=0; Gk < nirreps; Gk++) {
	Gijk = Gi ^ Gj ^ Gk;

	for(Gab=0; Gab < nirreps; Gab++) {
	  Gc = Gab ^ Gijk; /* totally symmetric */
	  W1[Gab] = global_dpd_->dpd_block_matrix(F.params->coltot[Gab], virtpi[Gc]);
	}
	for(Gb=0; Gb < nirreps; Gb++) {
	  Gca = Gb ^ Gijk; /* totally symmtric */
	  W2[Gb] = global_dpd_->dpd_block_matrix(virtpi[Gb], F.params->coltot[Gca]);
	}

	for(i=0; i < occpi[Gi]; i++) {
	  I = occ_off[Gi] + i;
	  for(j=0; j < occpi[Gj]; j++) {
	    J = occ_off[Gj] + j;
	    for(k=0; k < occpi[Gk]; k++) {
	      K = occ_off[Gk] + k;

	      global_dpd_->T3_AAA(W1, nirreps, I, Gi, J, Gj, K, Gk, &T2, &F, &E, &fIJ, &fAB,
		     occpi, occ_off, virtpi, vir_off, 0.0);

	      /* Z_MCIK <-- 1/2 t_IJKABC <MJ||AB> */

	      ik = ZIFLN.params->colidx[I][K];

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gab = Gmj = Gm ^ Gj; /* totally symmetric? */
		Gc = Gab ^ Gijk;     /* totally symmetric? */
		Gmc = Gm ^ Gc;

		nrows = Dints.params->coltot[Gmj];
		ncols = virtpi[Gc];
		Z = init_array(ncols);

		if(nrows && ncols) {
		  for(m=0; m < occpi[Gm]; m++) {
		    M = occ_off[Gm] + m;
		    mj = Dints.params->rowidx[M][J];
		    C_DGEMV('t', nrows, ncols, 0.5, W1[Gab][0], ncols, Dints.matrix[Gmj][mj], 1, 0.0, Z, 1);
		    for(c=0; c < virtpi[Gc]; c++) {
		      C = vir_off[Gc] + c;
		      mc = ZIFLN.params->rowidx[M][C];
		      ZIFLN.matrix[Gmc][mc][ik] += Z[c];
		    }
		  }
		}
		free(Z);
	      }

	      /* Z_ACEK (KE,CA) <-- -1/2 t_IJKABC <IJ||EB> */

	      ij = Dints.params->rowidx[I][J];
	      /* sort W(AB,C) to W(B,CA) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk; /* totally symmetric */
		for(ab=0; ab < F.params->coltot[Gab]; ab++) {
		  A = F.params->colorb[Gab][ab][0];
		  B = F.params->colorb[Gab][ab][1];
		  Gb = F.params->ssym[B];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    ca = F.params->colidx[C][A];
		    W2[Gb][b][ca] = W1[Gab][ab][c];
		  }
		}
	      }

	      for(Ge=0; Ge < nirreps; Ge++) {
		Gb = Ge ^ Gij; /* totally symmetric */
		Gca = Gke = Gk ^ Ge; /* totally symmetric */

		nrows = virtpi[Ge];
		ncols = ZDFAN.params->coltot[Gke];
		nlinks = virtpi[Gb];

		eb = Dints.col_offset[Gij][Ge];
		ke = ZDFAN.row_offset[Gke][K];
		ZDFAN.matrix[Gke] = global_dpd_->dpd_block_matrix(nrows, ncols);
		global_dpd_->buf4_mat_irrep_rd_block(&ZDFAN, Gke, ke, nrows);

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 'n', nrows, ncols, nlinks, -0.5, &(Dints.matrix[Gij][ij][eb]), nlinks,
			  W2[Gb][0], ncols, 1.0, ZDFAN.matrix[Gke][0], ncols);

		global_dpd_->buf4_mat_irrep_wrt_block(&ZDFAN, Gke, ke, nrows);
		global_dpd_->free_dpd_block(ZDFAN.matrix[Gke], nrows, ncols);
	      }

	    } /* k */
	  } /* j */
	} /* i */

	for(Gab=0; Gab < nirreps; Gab++) {
	  Gc = Gab ^ Gijk; /* totally symmetric */
	  global_dpd_->free_dpd_block(W1[Gab], F.params->coltot[Gab], virtpi[Gc]);
	}
	for(Gb=0; Gb < nirreps; Gb++) {
	  Gca = Gb ^ Gijk; /* totally symmtric */
	  global_dpd_->free_dpd_block(W2[Gb], virtpi[Gb], F.params->coltot[Gca]);
	}

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(W1);
  free(W2);

  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&T2);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fAB);

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_wrt(&ZIFLN, h);
    global_dpd_->buf4_mat_irrep_close(&ZIFLN, h);
  }
  global_dpd_->buf4_close(&ZIFLN);

  global_dpd_->buf4_close(&ZDFAN);

  for(h=0; h < nirreps; h++) global_dpd_->buf4_mat_irrep_close(&Dints, h);
  global_dpd_->buf4_close(&Dints);
}

void cc3_t3z_RHF_AAB(void)
{
  int h, nirreps;
  int *occ_off, *occpi;
  int *vir_off, *virtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc, Gab;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int ab;
  double ***W1, ***W2;
  dpdbuf4 T2AA, T2AB, T2BA, EAA, EAB, EBA, FAA, FAB, FBA;
  dpdfile2 fIJ, fAB, fij, fab;
  dpdbuf4 Dints, DAAints;
  dpdbuf4 ZIFLN , ZIfLn, ZDFAN, ZDfAn;
  int ji, kj;
  int Gbc, Gca, bc, ca;
  int Gmk, Gmi, mk, mi, Gma, ma, Gmb, mb;
  int Gm, m, M;
  double *Z;
  int nrows, ncols, nlinks;
  int Gik, ik, Gki, ki;
  int ba, Ge, Gje, ec, je, ea;
  int Gmj, mj, Gmc, mc;
  int Gij, ij, Gke, ke, eb;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  occ_off = moinfo.occ_off;
  virtpi = moinfo.virtpi;
  vir_off = moinfo.vir_off;

  global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_init(&DAAints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Dints, h);
    global_dpd_->buf4_mat_irrep_rd(&Dints, h);

    global_dpd_->buf4_mat_irrep_init(&DAAints, h);
    global_dpd_->buf4_mat_irrep_rd(&DAAints, h);
  }

  global_dpd_->buf4_init(&ZIFLN, PSIF_CC3_MISC, 0, 10, 0, 10, 0, 0, "CC3 ZIFLN");
  global_dpd_->buf4_init(&ZIfLn, PSIF_CC3_MISC, 0, 10, 0, 10, 0, 0, "CC3 ZIfLn");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&ZIFLN, h);
    global_dpd_->buf4_mat_irrep_rd(&ZIFLN, h);

    global_dpd_->buf4_mat_irrep_init(&ZIfLn, h);
  }

  global_dpd_->buf4_init(&ZDFAN, PSIF_CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZDFAN (NA,FD)");
  global_dpd_->buf4_init(&ZDfAn, PSIF_CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZDfAn (nA,fD)");
  global_dpd_->buf4_scm(&ZDfAn, 0.0);

  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");

  global_dpd_->buf4_init(&T2AA, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  global_dpd_->buf4_init(&T2AB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->buf4_init(&T2BA, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  global_dpd_->buf4_init(&FAA, PSIF_CC3_HET1, 0, 10, 5, 10, 7, 0, "CC3 WABEI (IE,B>A)");
  global_dpd_->buf4_init(&FAB, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WaBeI (Ie,Ba)");
  global_dpd_->buf4_init(&FBA, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
  global_dpd_->buf4_init(&EAA, PSIF_CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 WMBIJ (I>J,MB)");
  global_dpd_->buf4_init(&EAB, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
  global_dpd_->buf4_init(&EBA, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WmBiJ (iJ,mB)");

  /* target T3 amplitudes go in here */
  W1 = (double ***) malloc(nirreps * sizeof(double **));
  W2 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gi ^ Gj;
      for(Gk=0; Gk < nirreps; Gk++) {
	Gijk = Gi ^ Gj ^ Gk;
	Gik = Gki = Gi ^ Gk;

	for(Gab=0; Gab < nirreps; Gab++) {
	  Gc = Gab ^ Gijk; /* totally symmtric */
	  W1[Gab] = global_dpd_->dpd_block_matrix(FAA.params->coltot[Gab], virtpi[Gc]);
	}
	for(Ga=0; Ga < nirreps; Ga++) {
	  Gbc = Ga ^ Gijk; /* totally symmtric */
	  W2[Ga] = global_dpd_->dpd_block_matrix(virtpi[Ga], FAA.params->coltot[Gbc]);
	}

	for(i=0; i < occpi[Gi]; i++) {
	  I = occ_off[Gi] + i;
	  for(j=0; j < occpi[Gj]; j++) {
	    J = occ_off[Gj] + j;
	    for(k=0; k < occpi[Gk]; k++) {
	      K = occ_off[Gk] + k;

	      global_dpd_->T3_AAB(W1, nirreps, I, Gi, J, Gj, K, Gk, &T2AA, &T2AB, &T2BA,
		     &FAA, &FAB, &FBA, &EAA, &EAB, &EBA, &fIJ, &fij, &fAB, &fab,
		     occpi, occ_off, occpi, occ_off, virtpi, vir_off, virtpi, vir_off, 0.0);

	      /* Z_MAJI <-- t_ABcIJk <Mk|Bc> */
	      ji = ZIFLN.params->colidx[J][I];
	      /* sort W(AB,c) --> W(A,Bc) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk; /* totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Ga = FAA.params->rsym[A];
		  Gb = FAA.params->ssym[B];
		  a = A - vir_off[Ga];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    bc = FAA.params->colidx[B][C];
		    W2[Ga][a][bc] = W1[Gab][ab][c];
		  }
		}
	      }

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gbc = Gmk = Gm ^ Gk;  /* totally symmetric */
		Ga = Gbc ^ Gijk;      /* totally symmetric */
		Gma = Gm ^ Ga;

		nrows = virtpi[Ga];
		ncols = Dints.params->coltot[Gmk];

		Z = init_array(nrows);

		if(nrows && ncols) {
		  for(m=0; m < occpi[Gm]; m++) {
		    M = occ_off[Gm] + m;
		    mk = Dints.params->rowidx[M][K];
		    C_DGEMV('n', nrows, ncols, 1.0, W2[Ga][0], ncols, Dints.matrix[Gmk][mk], 1, 0.0, Z, 1);
		    for(a=0; a < virtpi[Ga]; a++) {
		      A = vir_off[Ga] + a;
		      ma = ZIFLN.params->rowidx[M][A];
		      ZIFLN.matrix[Gma][ma][ji] += Z[a];
		    }
		  }
		}

		free(Z);
	      }

	      /* ZMbKj <-- t_ijKabC <Mi|Ca> */
	      kj = ZIfLn.params->colidx[K][J];
	      /* sort W(ab,C) to W(b,Ca) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk; /* totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Ga = FAA.params->rsym[A];
		  Gb = FAA.params->ssym[B];
		  a = A - vir_off[Ga];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    ca = FAA.params->colidx[C][A];
		    W2[Gb][b][ca] = W1[Gab][ab][c];
		  }
		}
	      }

	      for(Gm=0; Gm < nirreps; Gm++) {
		Gca = Gmi = Gm ^ Gi;  /* totally symmetric */
		Gb = Gca ^ Gijk;      /* totally symmetric */
		Gmb = Gm ^ Gb;

		nrows = virtpi[Gb];
		ncols = Dints.params->coltot[Gmi];

		Z = init_array(nrows);

		if(nrows && ncols) {
		  for(m=0; m < occpi[Gm]; m++) {
		    M = occ_off[Gm] + m;
		    mi = Dints.params->rowidx[M][I];
		    C_DGEMV('n', nrows, ncols, 1.0, W2[Gb][0], ncols, Dints.matrix[Gmi][mi], 1, 0.0, Z, 1);
		    for(b=0; b < virtpi[Gb]; b++) {
		      B = vir_off[Gb] + b;
		      mb = ZIfLn.params->rowidx[M][B];
		      ZIfLn.matrix[Gmb][mb][kj] += Z[b];
		    }
		  }
		}

		free(Z);
	      }

	      /* Z_McIk <-- 1/2 t_IJkABc <MJ||AB> */
	      ik = ZIfLn.params->colidx[I][K];

	      for(Gm=0; Gm < nirreps; Gm++) {

		Gab = Gmj = Gm ^ Gj; /* totally symmetric */
		Gc = Gab ^ Gijk; /* totally symmetric */
		Gmc = Gm ^ Gc;

		nrows = DAAints.params->coltot[Gmj];
		ncols = virtpi[Gc];
		Z = init_array(ncols);

		if(nrows && ncols) {
		  for(m=0; m < occpi[Gm]; m++) {
		    M = occ_off[Gm] + m;
		    mj = DAAints.params->rowidx[M][J];
		    C_DGEMV('t', nrows, ncols, 0.5, W1[Gab][0], ncols, DAAints.matrix[Gmj][mj], 1, 0.0, Z, 1);
		    for(c=0; c < virtpi[Gc]; c++) {
		      C = vir_off[Gc] + c;
		      mc = ZIfLn.params->rowidx[M][C];
		      ZIfLn.matrix[Gmc][mc][ik] += Z[c];
		    }
		  }
		}
		free(Z);
	      }

	      /* Z_ABEJ (JE,BA) <-- - t_IJkABc <Ik|Ec> */
	      ik = Dints.params->rowidx[I][K];
	      /* sort W(AB,c) to W(c,BA) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk; /* totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  ba = FAA.params->colidx[B][A];
		  for(c=0; c < virtpi[Gc]; c++) {
		    W2[Gc][c][ba] = W1[Gab][ab][c];
		  }
		}
	      }

	      for(Ge=0; Ge < nirreps; Ge++) {
		Gc = Gik ^ Ge; /* totally symmetric */
		Gje = Gj ^ Ge;

		nrows = virtpi[Ge];
		ncols = ZDFAN.params->coltot[Gje];
		nlinks = virtpi[Gc];

		ec = Dints.col_offset[Gik][Ge];
		je = ZDFAN.row_offset[Gje][J];
		ZDFAN.matrix[Gje] = global_dpd_->dpd_block_matrix(nrows, ncols);
		global_dpd_->buf4_mat_irrep_rd_block(&ZDFAN, Gje, je, nrows);

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 'n', nrows, ncols, nlinks, -1.0, &(Dints.matrix[Gik][ik][ec]), nlinks,
			  W2[Gc][0], ncols, 1.0, ZDFAN.matrix[Gje][0], ncols);

		global_dpd_->buf4_mat_irrep_wrt_block(&ZDFAN, Gje, je, nrows);
		global_dpd_->free_dpd_block(ZDFAN.matrix[Gje], nrows, ncols);
	      }

	      /* Z_CbEj <-- - t_ijKabC <Ki|Ea> */
	      ki = Dints.params->rowidx[K][I];
	      /* sort W(ab,C) to W(a,bC) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk; /* totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Ga = FAA.params->rsym[A];
		  a = A - vir_off[Ga];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    bc = FAA.params->colidx[B][C];
		    W2[Ga][a][bc] = W1[Gab][ab][c];
		  }
		}
	      }

	      for(Ge=0; Ge < nirreps; Ge++) {
		Ga = Ge ^ Gki; /* totally symmetric */
		Gje = Gj ^ Ge;

		nrows = virtpi[Ge];
		ncols = ZDfAn.params->coltot[Gje];
		nlinks = virtpi[Ga];

		ea = Dints.col_offset[Gki][Ge];
		je = ZDfAn.row_offset[Gje][J];
		ZDfAn.matrix[Gje] = global_dpd_->dpd_block_matrix(nrows, ncols);
		global_dpd_->buf4_mat_irrep_rd_block(&ZDfAn, Gje, je, nrows);

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 'n', nrows, ncols, nlinks, -1.0, &(Dints.matrix[Gki][ki][ea]), nlinks,
			  W2[Ga][0], ncols, 1.0, ZDfAn.matrix[Gje][0], ncols);

		global_dpd_->buf4_mat_irrep_wrt_block(&ZDfAn, Gje, je, nrows);
		global_dpd_->free_dpd_block(ZDfAn.matrix[Gje], nrows, ncols);
	      }

	      /* Z_AcEk <-- -1/2 t_IJkABc <IJ||EB> */
	      ij = DAAints.params->rowidx[I][J];
	      /* sort W(AB,C) to W(B,CA) */
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk; /* totally symmetric */
		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Gb = FAA.params->ssym[B];
		  b = B - vir_off[Gb];
		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    ca = FAA.params->colidx[C][A];
		    W2[Gb][b][ca] = W1[Gab][ab][c];
		  }
		}
	      }

	      for(Ge=0; Ge < nirreps; Ge++) {
		Gb = Ge ^ Gij; /* totally symmetric */
		Gca = Gke = Gk ^ Ge; /* totally symmetric */

		nrows = virtpi[Ge];
		ncols = ZDfAn.params->coltot[Gke];
		nlinks = virtpi[Gb];

		eb = DAAints.col_offset[Gij][Ge];
		ke = ZDfAn.row_offset[Gke][K];
		ZDfAn.matrix[Gke] = global_dpd_->dpd_block_matrix(nrows, ncols);
		global_dpd_->buf4_mat_irrep_rd_block(&ZDfAn, Gke, ke, nrows);

		if(nrows && ncols && nlinks)
		  C_DGEMM('n', 'n', nrows, ncols, nlinks, -0.5, &(DAAints.matrix[Gij][ij][eb]), nlinks,
			  W2[Gb][0], ncols, 1.0, ZDfAn.matrix[Gke][0], ncols);

		global_dpd_->buf4_mat_irrep_wrt_block(&ZDfAn, Gke, ke, nrows);
		global_dpd_->free_dpd_block(ZDfAn.matrix[Gke], nrows, ncols);
	      }

	    } /* k */
	  } /* j */
	} /* i */

	for(Gab=0; Gab < nirreps; Gab++) {
	  Gc = Gab ^ Gijk;  /* totally symmetric */
	  global_dpd_->free_dpd_block(W1[Gab], FAA.params->coltot[Gab], virtpi[Gc]);
	}
	for(Ga=0; Ga < nirreps; Ga++) {
	  Gbc = Ga ^ Gijk; /* totally symmtric */
	  global_dpd_->free_dpd_block(W2[Ga], virtpi[Ga], FAA.params->coltot[Gbc]);
	}

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(W1);
  free(W2);

  global_dpd_->buf4_close(&EAA);
  global_dpd_->buf4_close(&EAB);
  global_dpd_->buf4_close(&EBA);
  global_dpd_->buf4_close(&FAA);
  global_dpd_->buf4_close(&FAB);
  global_dpd_->buf4_close(&FBA);
  global_dpd_->buf4_close(&T2AA);
  global_dpd_->buf4_close(&T2AB);
  global_dpd_->buf4_close(&T2BA);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fAB);
  global_dpd_->file2_close(&fij);
  global_dpd_->file2_close(&fab);

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_wrt(&ZIFLN, h);
    global_dpd_->buf4_mat_irrep_close(&ZIFLN, h);

    global_dpd_->buf4_mat_irrep_wrt(&ZIfLn, h);
    global_dpd_->buf4_mat_irrep_close(&ZIfLn, h);
  }
  global_dpd_->buf4_close(&ZIFLN);
  global_dpd_->buf4_close(&ZIfLn);

  global_dpd_->buf4_sort(&ZDFAN, PSIF_CC3_MISC, qpsr, 11, 5, "CC3 ZDFAN (AN,DF)");
  global_dpd_->buf4_close(&ZDFAN);
  global_dpd_->buf4_sort(&ZDfAn, PSIF_CC3_MISC, qpsr, 11, 5, "CC3 ZDfAn (An,Df)");
  global_dpd_->buf4_close(&ZDfAn);

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_close(&Dints, h);
    global_dpd_->buf4_mat_irrep_close(&DAAints, h);
  }
  global_dpd_->buf4_close(&Dints);
  global_dpd_->buf4_close(&DAAints);
}


}} // namespace psi::cclambda
