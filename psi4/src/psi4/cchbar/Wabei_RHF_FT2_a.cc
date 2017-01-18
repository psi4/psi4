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
    \ingroup CCHBAR
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

namespace psi { namespace cchbar {

/* Wabei_RHF_FT2_a(): Computes the following contribution to the Wabei HBAR
** matrix elements for RHF references:
**
** Z_aeib = <am|ef> [ 2 t_mi^fb - t_mi^bf ] - <mf|ae> t_mi^fb
**
** This is the final term of the Wabef HBAR element's contribution to Wabei.
** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).]
**
** TDC, March 2004
*/

void Wabei_RHF_FT2_a(void)
{
  int h, nirreps;
  int Gam, Ga, Gm;
  int Gef, Gmf, Ge, Gf, Gae;
  int a, A, e, E, m, M, f, FF;
  int mf, am , ef;
  int *virtpi, *vir_off, *occpi, *occ_off;
  double ***W;
  double value;
  dpdbuf4 F, T2, Z;

  nirreps = moinfo.nirreps;
  virtpi = moinfo.virtpi;
  vir_off = moinfo.vir_off;
  occpi = moinfo.occpi;
  occ_off = moinfo.occ_off;

  /* Zaeib <-- - <mf|ae> t_mi^fb */
  /**** this contraction still requires storage of a full <ia|bc> set *****/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 5, 10, 5, 10, 0, "Z(AE,ib)");
  global_dpd_->contract444(&F, &T2, &Z, 1, 1, -1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /* Zaeib <-- <am|ef> [ 2 t_mi^fb - t_mi^bf ] */

   global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, 0, 5, 10, 5, 10, 0, "Z(AE,ib)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&T2, h);
    global_dpd_->buf4_mat_irrep_rd(&T2, h);
  }

  W = (double ***) malloc(nirreps * sizeof(double **));

  for(Ga=0; Ga < nirreps; Ga++) {

    /* allocate space for the input integrals */
    for(Gm=0; Gm < nirreps; Gm++) {
      Gam = Ga^Gm;
      F.matrix[Gam] = global_dpd_->dpd_block_matrix(occpi[Gm], F.params->coltot[Gam]);
    }

    /* allocate space for the reordered integrals, W[e][mf] and target, Z[e][ib] */
    for(Ge=0; Ge < nirreps; Ge++) {
      Gae = Ga^Ge; Gmf = Gae;
      W[Ge] = global_dpd_->dpd_block_matrix(virtpi[Ge], Z.params->coltot[Gmf]);
      Z.matrix[Gae] = global_dpd_->dpd_block_matrix(virtpi[Ge], Z.params->coltot[Gmf]);
    }

    for(a=0; a < virtpi[Ga]; a++) {
      A = vir_off[Ga] + a;

      for(Gm=0; Gm < nirreps; Gm++) {
	Gam = Ga^Gm;
	Gef = Gam;

 	/* read all <am|ef> integrals for given orbital index a --> F[m][ef] */
	global_dpd_->buf4_mat_irrep_rd_block(&F, Gam, F.params->start13[Gam][A], occpi[Gm]);
      }

      for(Ge=0; Ge < nirreps; Ge++) {
	Gae = Ga^Ge; Gmf = Gae;

	/* sort F[m][ef] --> W[e][mf] */
	for(e=0; e < virtpi[Ge]; e++) {
	  E = vir_off[Ge] + e;

	  for(Gm=0; Gm < nirreps; Gm++) {
	    Gf = Gmf^Gm;
	    Gam = Ga^Gm;

	    for(m=0; m < occpi[Gm]; m++) {
	      M = occ_off[Gm] + m;

	      for(f=0; f < virtpi[Gf]; f++) {
		FF = vir_off[Gf] + f;

		mf = Z.params->colidx[M][FF];
		am = F.params->rowidx[A][M];
		ef = F.params->colidx[E][FF];

		W[Ge][e][mf] = F.matrix[Gam][am-F.params->start13[Gam][A]][ef];
	      } /* f */
	    } /* m */
	  } /* Gm */
	} /* e */

	/* read all existing Z(ae,ib) for a given orbital index a -->  Z[e][ib] */
	global_dpd_->buf4_mat_irrep_rd_block(&Z, Gae, Z.params->start13[Gae][A], virtpi[Ge]);

	/* contract W[e][mf] * T[mf][ib] -> Z[e][ib] */
	if(virtpi[Ge] && Z.params->coltot[Gmf])
	  C_DGEMM('n','n', virtpi[Ge], Z.params->coltot[Gmf], Z.params->coltot[Gmf], 1.0,
		  W[Ge][0], Z.params->coltot[Gmf], T2.matrix[Gmf][0], Z.params->coltot[Gmf],
		  1.0, Z.matrix[Gae][0], Z.params->coltot[Gmf]);

	/* write out the new Z(ae,ib) */
	global_dpd_->buf4_mat_irrep_wrt_block(&Z, Gae, Z.params->start13[Gae][A], virtpi[Ge]);

      } /* Ge */

    } /* a */

    for(Gm=0; Gm < nirreps; Gm++) {
      Gam = Ga^Gm;
      global_dpd_->free_dpd_block(F.matrix[Gam], occpi[Gm], F.params->coltot[Gam]);
    }

    for(Ge=0; Ge < nirreps; Ge++) {
      Gae = Ga^Ge; Gmf = Gae;
      global_dpd_->free_dpd_block(W[Ge], virtpi[Ge], Z.params->coltot[Gmf]);
      global_dpd_->free_dpd_block(Z.matrix[Gae], virtpi[Ge], Z.params->coltot[Gmf]);
    }

  } /* Ga */

  for(h=0; h < nirreps; h++) global_dpd_->buf4_mat_irrep_close(&T2, h);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_HBAR, prqs, 11, 5, "WAbEi (Ei,Ab)", 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&F);


}

}} // namespace psi::cchbar
