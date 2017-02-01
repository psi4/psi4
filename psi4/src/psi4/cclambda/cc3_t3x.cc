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
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void cc3_t3x(void)
{
  if(params.ref == 0) {
    int h, nirreps;
    int *occ_off, *occpi;
    int *vir_off, *virtpi;
    int Gi, Gj, Gk, Gijk;
    int Ga, Gb, Gc, Gab;
    int i, j, k, I, J, K;
    int a, b, c, A, B, C;
    int ab;
    double ***W1;
    dpdbuf4 T2, E, F, T2AA, T2AB, T2BA, EAA, EAB, EBA, FAA, FAB, FBA;
    dpdfile2 fIJ, fAB, fij, fab;
    dpdfile2 XLD;
    dpdbuf4 L2, L2AB;
    int Gij, ij, Gbc, bc, Gjk, jk;
    int nrows, ncols;
    int **W_offset, offset;

    nirreps = moinfo.nirreps;
    occpi = moinfo.occpi;
    occ_off = moinfo.occ_off;
    virtpi = moinfo.virtpi;
    vir_off = moinfo.vir_off;

    W_offset = init_int_matrix(nirreps, nirreps);
    for(Gab=0; Gab < nirreps; Gab++) {
      for(Ga=0,offset=0; Ga < nirreps; Ga++) {
	Gb = Ga ^ Gab;
	W_offset[Gab][Ga] = offset;
	offset += virtpi[Ga] * virtpi[Gb];
      }
    }

    global_dpd_->file2_init(&XLD, PSIF_CC3_MISC, 0, 0, 1, "CC3 XLD");
    global_dpd_->file2_mat_init(&XLD);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, 0, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_init(&L2AB, PSIF_CC_LAMBDA, 0, 0, 5, 0, 5, 0, "LIjAb");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&L2, h);
      global_dpd_->buf4_mat_irrep_rd(&L2, h);

      global_dpd_->buf4_mat_irrep_init(&L2AB, h);
      global_dpd_->buf4_mat_irrep_rd(&L2AB, h);
    }

    global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
    global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&F, PSIF_CC3_HET1, 0, 10, 5, 10, 7, 0, "CC3 WABEI (IE,B>A)");
    global_dpd_->buf4_init(&E, PSIF_CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 WMBIJ (I>J,MB)");

    T2AA = T2;
    global_dpd_->buf4_init(&T2AB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&T2BA, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    FAA = F;
    global_dpd_->buf4_init(&FAB, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WaBeI (Ie,Ba)");
    global_dpd_->buf4_init(&FBA, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
    EAA = E;
    global_dpd_->buf4_init(&EAB, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
    global_dpd_->buf4_init(&EBA, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WmBiJ (iJ,mB)");

    /* target T3 amplitudes go in here */
    W1 = (double ***) malloc(nirreps * sizeof(double **));

    for(Gi=0; Gi < nirreps; Gi++) {
      for(Gj=0; Gj < nirreps; Gj++) {
	Gij = Gi ^ Gj;
	for(Gk=0; Gk < nirreps; Gk++) {
	  Gijk = Gi ^ Gj ^ Gk;
	  Gjk = Gj ^ Gk;

	  for(Gab=0; Gab < nirreps; Gab++) {
	    Gc = Gab ^ Gijk; /* totally symmetric */
	    W1[Gab] = global_dpd_->dpd_block_matrix(F.params->coltot[Gab], virtpi[Gc]);
	  }

	  for(i=0; i < occpi[Gi]; i++) {
	    I = occ_off[Gi] + i;
	    for(j=0; j < occpi[Gj]; j++) {
	      J = occ_off[Gj] + j;
	      for(k=0; k < occpi[Gk]; k++) {
		K = occ_off[Gk] + k;

		global_dpd_->T3_AAA(W1, nirreps, I, Gi, J, Gj, K, Gk, &T2, &F, &E, &fIJ, &fAB,
		       occpi, occ_off, virtpi, vir_off, 0.0);

		/* X_KC <-- 1/4 t_IJKABC <IJ||AB> */

		Gc = Gk;    /* assumes T1 is totally symmetric */
		Gab = Gij;  /* assumes <ij||ab> is totally symmetric */

		ij = L2.params->rowidx[I][J];

		nrows = L2.params->coltot[Gij];
		ncols = virtpi[Gc];

		if(nrows && ncols)
		  C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, L2.matrix[Gij][ij], 1,
			  1.0, XLD.matrix[Gk][k], 1);

		global_dpd_->T3_AAB(W1, nirreps, I, Gi, J, Gj, K, Gk, &T2AA, &T2AB, &T2BA,
		       &FAA, &FAB, &FBA, &EAA, &EAB, &EBA, &fIJ, &fij, &fAB, &fab,
		       occpi, occ_off, occpi, occ_off, virtpi, vir_off, virtpi, vir_off, 0.0);

		/* t_IA <-- t_IJkABc <Jk|Bc> */

		Ga = Gi;   /* assumes T1 is totally symmetric */
		Gbc = Gjk; /* assumes <jk|bc> is totally symmetric */

		jk = L2AB.params->rowidx[J][K];

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gb = Ga ^ Gab;
		  Gc = Gb ^ Gbc;

		  ab = W_offset[Gab][Ga];
		  bc = L2AB.col_offset[Gjk][Gb];

		  nrows = virtpi[Ga];
		  ncols = virtpi[Gb] * virtpi[Gc];

		  if(nrows && ncols)
		    C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][ab], ncols, &(L2AB.matrix[Gjk][jk][bc]), 1,
			    1.0, XLD.matrix[Gi][i], 1);

		}
		/* t_KC <-- 1/4 t_ijKabC <ij||ab> */

		Gc = Gk;  /* assumes T1 is totally symmetric */
		Gab = Gij; /* assumes <ij||ab> is totally symmetric */

		ij = L2.params->rowidx[I][J];

		nrows = L2.params->coltot[Gij];
		ncols = virtpi[Gc];

		if(nrows && ncols)
		  C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, L2.matrix[Gij][ij], 1,
			  1.0, XLD.matrix[Gk][k], 1);


	      } /* k */
	    } /* j */
	  } /* i */

	  for(Gab=0; Gab < nirreps; Gab++) {
	    Gc = Gab ^ Gijk; /* totally symmetric */
	    global_dpd_->free_dpd_block(W1[Gab], F.params->coltot[Gab], virtpi[Gc]);
	  }
	} /* Gk */
      } /* Gj */
    } /* Gi */

    free(W1);

    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&fIJ);
    global_dpd_->file2_close(&fAB);

    global_dpd_->buf4_close(&EAB);
    global_dpd_->buf4_close(&EBA);
    global_dpd_->buf4_close(&FAB);
    global_dpd_->buf4_close(&FBA);
    global_dpd_->buf4_close(&T2AB);
    global_dpd_->buf4_close(&T2BA);
    global_dpd_->file2_close(&fij);
    global_dpd_->file2_close(&fab);

    free_int_matrix(W_offset);

    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_close(&L2, h);
      global_dpd_->buf4_mat_irrep_close(&L2AB, h);
    }
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&L2AB);

    global_dpd_->file2_mat_wrt(&XLD);
    global_dpd_->file2_close(&XLD);
  }
}

}} // namespace psi::cclambda
