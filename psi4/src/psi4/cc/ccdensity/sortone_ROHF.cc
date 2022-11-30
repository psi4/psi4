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
** sortone_rohf(): Place all the components of the 1pdm into a large
** matrix, O (moinfo.opdm), which we also symmetrize by computing Opq
** = 1/2 (Opq + Oqp).  This matrix is later written to disk in dump()
** for subsequent backtransformation.  Note that the components of the
** 1pdm computed into the DIJ, Dij, DAB, Dab, DAI, Dai, DIA, and Dia
** matrices remain non-symmetric (e.g., DIJ neq DJI).
**
** This version doesn't work with frozen orbitals yet.
**
** TDC, 1/03
*/

void sortone_ROHF(const struct RHO_Params& rho_params) {
    dpdfile2 D;

    /* Sort A components first */
    Slice aocc_in_full(moinfo.frdocc, moinfo.frdocc + moinfo.occpi);
    Slice avir_in_vir(Dimension(moinfo.nirreps), moinfo.virtpi - moinfo.openpi);
    Slice avir_in_full(moinfo.frdocc + moinfo.occpi, moinfo.orbspi - moinfo.fruocc);
    Slice aocc_in_occ(Dimension(moinfo.nirreps), moinfo.occpi);

    Matrix O_a(moinfo.orbspi, moinfo.orbspi);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
    Matrix temp(&D);
    O_a.set_block(aocc_in_full, temp);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
    temp = Matrix(&D);
    auto temp2 = temp.get_block(avir_in_vir);
    O_a.set_block(avir_in_full, *temp2);
    global_dpd_->file2_close(&D);

    /* Note that this component of the density is stored occ-vir */
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
    temp = Matrix(&D);
    temp2 = temp.get_block(aocc_in_occ, avir_in_vir);
    temp2 = temp2->transpose();
    O_a.set_block(avir_in_full, aocc_in_full, temp2);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    temp = Matrix(&D);
    temp2 = temp.get_block(aocc_in_occ, avir_in_vir);
    O_a.set_block(aocc_in_full, avir_in_full, temp2);
    global_dpd_->file2_close(&D);

    /* Sort B components */
    // vir is stored in DPD as UOCC then SOCC.
    // The standard order used by Psi is SOCC then UOCC.
    // This inconsistency is why we need to treat uocc and socc separately.
    Slice docc_in_full(moinfo.frdocc, moinfo.frdocc + moinfo.occpi - moinfo.openpi);
    Slice docc_in_occ(Dimension(moinfo.nirreps), moinfo.occpi - moinfo.openpi);
    Slice socc_in_full(moinfo.frdocc + moinfo.clsdpi, moinfo.frdocc + moinfo.occpi);
    Slice socc_in_vir(moinfo.virtpi - moinfo.openpi, moinfo.virtpi);
    Slice uocc_in_vir(Dimension(moinfo.nirreps), moinfo.virtpi - moinfo.openpi);
    Slice uocc_in_full(moinfo.frdocc + moinfo.occpi, moinfo.orbspi - moinfo.fruocc);

    Matrix O_b(moinfo.orbspi, moinfo.orbspi);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.Dij_lbl);
    temp = Matrix(&D);
    temp2 = temp.get_block(docc_in_occ);
    O_b.set_block(docc_in_full, *temp2);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.Dab_lbl);
    temp = Matrix(&D);
    temp2 = temp.get_block(uocc_in_vir);
    O_b.set_block(uocc_in_full, *temp2);
    temp2 = temp.get_block(socc_in_vir);
    O_b.set_block(socc_in_full, *temp2);
    temp2 = temp.get_block(uocc_in_vir, socc_in_vir);
    O_b.set_block(uocc_in_full, socc_in_full, *temp2);
    temp2 = temp.get_block(socc_in_vir, uocc_in_vir);
    O_b.set_block(socc_in_full, uocc_in_full, *temp2);
    global_dpd_->file2_close(&D);

    /* Note that this component of the density is stored occ-vir */
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dai_lbl);
    temp = Matrix(&D);
    temp2 = temp.transpose();
    auto temp3 = temp2->get_block(socc_in_vir, docc_in_occ);
    O_b.set_block(socc_in_full, docc_in_full, *temp3);
    temp3 = temp2->get_block(uocc_in_vir, docc_in_occ);
    O_b.set_block(uocc_in_full, docc_in_full, *temp3);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.Dia_lbl);
    temp = Matrix(&D);
    temp2 = temp.get_block(docc_in_occ, socc_in_vir);
    O_b.set_block(docc_in_full, socc_in_full, *temp2);
    temp2 = temp.get_block(docc_in_occ, uocc_in_vir);
    O_b.set_block(docc_in_full, uocc_in_full, *temp2);
    global_dpd_->file2_close(&D);

    /* Symmetrize the onepdm */

    O_a.hermitivitize();
    O_b.hermitivitize();
    moinfo.opdm_a = O_a;
    moinfo.opdm_b = O_b;
    moinfo.opdm = Matrix(O_a);
    moinfo.opdm.add(O_b);
}
}
}  // namespace psi
