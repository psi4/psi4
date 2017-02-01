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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <stdio.h>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

    /* DUMP_RHF(): Mulliken-order the RHF-CC two-electron density and
    ** dump it to a file for subsequent backtransformation.  Basically
    ** all we have to do is swap indices two and three, e.g.
    **
    ** G'(pr,qs) = G(pq,rs)
    **
    ** In order for the Mulliken-ordered density to be valid for the
    ** backtransformation algorithm used in TRANSQT, the final density
    ** must have eight-fold permutational symmetry like the original
    ** integrals, so symmetrization of some components is necessary.
    **
    ** Note that the CC-ordered indices must be translated inside
    ** dpd_buf4_dump() to QT ordering, as required by TRANSQT, hence
    ** the qt_occ and qt_vir arrays.
    **
    ** The bk_pack option to dpd_buf4_dump() indicates whether the
    ** Mulliken-ordered twopdm component is bra-ket symmetric
    ** already.  If so, then only the lower-triangle needs to be sent
    ** to TRANSQT2.
    **
    ** I don't remember why I have the swap23 option in
    ** dpd_buf4_dump(), particularly since I do an explicit sort of
    ** the density components using DPD functions.  It would actually
    ** be better to let the dump function handle this part to reduce
    ** disk usage.
    **
    ** Finally, I note that Gibja is multiplied by 0.5 before it is
    ** dumped to disk.  I don't remember why this is necessary.
    **
    ** TDC, last updated 2/08
    */

    void dump_RHF(struct iwlbuf *OutBuf, struct RHO_Params rho_params)
    {
      int nirreps, nmo, nfzv;
      int *qt_occ, *qt_vir;
      int h, row, col, p, q, r, s;
      dpdbuf4 G;

      qt_occ = moinfo.qt_occ;  qt_vir = moinfo.qt_vir;
      nirreps = moinfo.nirreps;
      nmo = moinfo.nmo;
      nfzv = moinfo.nfzv;

      psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
      psio_write_entry(PSIF_MO_OPDM, rho_params.opdm_lbl, (char *) moinfo.opdm[0],
		       sizeof(double)*(nmo-nfzv)*(nmo-nfzv));
      psio_close(PSIF_MO_OPDM, 1);

      if (!params.onepdm) {
	psio_open(PSIF_MO_LAG, PSIO_OPEN_OLD);
	psio_write_entry(PSIF_MO_LAG, "MO-basis Lagrangian", (char *) moinfo.I[0],
			 sizeof(double)*nmo*nmo);
	psio_close(PSIF_MO_LAG, 1);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prqs, 0, 0, "G(IK,JL)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "G(IK,JL)");
	global_dpd_->buf4_dump(&G, OutBuf, qt_occ, qt_occ, qt_occ, qt_occ, 1, 0);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prqs, 0, 10, "G(IK,JA)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "G(IK,JA)");
	global_dpd_->buf4_dump(&G, OutBuf, qt_occ, qt_occ, qt_occ, qt_vir, 0, 0);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP9, prqs, 10, 10, "G(IA,JB)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP9, 0, 10, 10, 10, 10, 0, "G(IA,JB)");
	global_dpd_->buf4_symm(&G);
	global_dpd_->buf4_dump(&G, OutBuf, qt_occ, qt_vir, qt_occ, qt_vir, 1, 0);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIbJa");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prqs, 0, 5, "G(IJ,AB)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "G(IJ,AB)");
	//        dpd_buf4_scm(&G, 0.5);  /* why is this necessary? */
	global_dpd_->buf4_dump(&G, OutBuf, qt_occ, qt_occ, qt_vir, qt_vir, 0, 0);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prqs, 5, 10, "G(ca,IB)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 5, 10, 5, 10, 0, "G(ca,IB)");
	global_dpd_->buf4_dump(&G, OutBuf, qt_vir, qt_vir, qt_occ, qt_vir, 0, 0);
	global_dpd_->buf4_close(&G);

	global_dpd_->buf4_init(&G, PSIF_CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
	global_dpd_->buf4_sort(&G, PSIF_CC_TMP0, prqs, 5, 5, "G(AC,BD)");
	global_dpd_->buf4_close(&G);
	global_dpd_->buf4_init(&G, PSIF_CC_TMP0, 0, 5, 5, 5, 5, 0, "G(AC,BD)");
	global_dpd_->buf4_dump(&G, OutBuf, qt_vir, qt_vir, qt_vir, qt_vir, 1, 0);
	global_dpd_->buf4_close(&G);

      }
    }

  }} // namespace psi::ccdensity
