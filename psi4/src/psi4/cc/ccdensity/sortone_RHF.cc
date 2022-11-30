/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {
#include "psi4/psifiles.h"

/*
** sortone_rhf(): Place all the components of the 1pdm into a large
** matrix, O (moinfo.opdm), which we also symmetrize by computing Opq
** = 1/2 (Opq + Oqp).  This matrix is later written to disk in dump()
** for subsequent backtransformation.  Note that the components of the
** 1pdm computed into the DIJ, Dij, DAB, Dab, DAI, Dai, DIA, and Dia
** matrices remain non-symmetric (e.g., DIJ neq DJI).
**
** This version doesn't work with frozen orbitals yet.
**
** TDC, 1/03
**
** NB: For now, I just multiply the components by two to account for
** spin adaptation.  The factors of two should later move into the I
** build itself.
**
** TDC, 2/2008
*/

void sortone_RHF(const struct RHO_Params& rho_params) {
    dpdfile2 D;

    Slice occ_slice(moinfo.frdocc, moinfo.frdocc + moinfo.occpi);
    Slice vir_slice(moinfo.frdocc + moinfo.occpi, moinfo.orbspi - moinfo.fruocc);

    /* O = block_matrix(nmo-nfzc,nmo-nfzc); */
    Matrix O(moinfo.orbspi, moinfo.orbspi);

    /* Sort A components first */
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
    Matrix temp(&D);
    O.set_block(occ_slice, temp);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
    temp = Matrix(&D);
    O.set_block(vir_slice, temp);
    global_dpd_->file2_close(&D);

    /* Note that this component of the density is stored occ-vir */
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
    temp = Matrix(&D);
    temp = *temp.transpose();
    O.set_block(vir_slice, occ_slice, temp);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    temp = Matrix(&D);
    O.set_block(occ_slice, vir_slice, temp);
    global_dpd_->file2_close(&D);

    O.hermitivitize();
    O.scale(2);
    moinfo.opdm = O;
}
}
}  // namespace psi
