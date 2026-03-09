/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
#include "psi4/libiwl/iwl.h"
#include "psi4/libmints/wavefunction.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {

/* SORTI_ROHF(): Place all the components of the ROHF Lagrangian into
** a large matrix, I (moinfo.I), which we also symmetrize by computing
** Ipq = 1/2 (Ipq + Iqp).  This matrix is later written to disk in
** dump() for subsequent backtransformation.  Note that some of the
** components of the Lagrangian computed into the IIJ, Iij, IIA, and
** Iia matrices remain non-symmetric (e.g., IIJ neq IJI).  I re-used
** my sortone.c code here, so don't let some of the variable names
** confuse you. */

void sortI_ROHF(Wavefunction& wfn) {
    dpdfile2 D;

    Slice aocc_in_full(moinfo.frdocc, moinfo.frdocc + moinfo.occpi);
    Slice avir_in_vir(Dimension(moinfo.nirreps), moinfo.virtpi - moinfo.openpi);
    Slice avir_in_full(moinfo.frdocc + moinfo.occpi, moinfo.orbspi - moinfo.fruocc);
    Slice aocc_in_occ(Dimension(moinfo.nirreps), moinfo.occpi);

    auto O = std::make_shared<Matrix>("Lagrangian matrix", moinfo.orbspi, moinfo.orbspi);

    /* Sort alpha components first */
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "I(I,J)");
    Matrix temp(&D);
    O->set_block(aocc_in_full, temp);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "I'AB");
    temp = Matrix(&D);
    auto temp2 = temp.get_block(avir_in_vir, avir_in_vir);
    O->set_block(avir_in_full, *temp2);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "I(I,A)");
    temp = Matrix(&D);
    temp2 = temp.get_block(aocc_in_occ, avir_in_vir);
    O->set_block(aocc_in_full, avir_in_full, temp2);
    temp2 = temp2->transpose();
    O->set_block(avir_in_full, aocc_in_full, temp2);
    global_dpd_->file2_close(&D);

    /* Sort beta components */
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

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "I(i,j)");
    temp = Matrix(&D);
    temp2 = temp.get_block(docc_in_occ);
    O_b.set_block(docc_in_full, *temp2);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "I'ab");
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

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "I(i,a)");
    temp = Matrix(&D);
    temp2 = temp.get_block(docc_in_occ, socc_in_vir);
    O_b.set_block(docc_in_full, socc_in_full, *temp2);
    temp2 = temp2->transpose();
    O_b.set_block(socc_in_full, docc_in_full, *temp2);
    temp2 = temp.get_block(docc_in_occ, uocc_in_vir);
    O_b.set_block(docc_in_full, uocc_in_full, *temp2);
    temp2 = temp2->transpose();
    O_b.set_block(uocc_in_full, docc_in_full, *temp2);
    global_dpd_->file2_close(&D);

    O->add(O_b);
    O->hermitivitize();
    O->scale(-1.0);

    wfn.set_lagrangian(linalg::triplet(wfn.Ca(), O, wfn.Ca(), false, false, true));
}

}  // namespace ccdensity
}  // namespace psi
