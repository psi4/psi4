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
** sortone_uhf(): Place all the components of the 1pdm into two
** spin-factored matrices, O_a (moinfo.opdm_a) and O_b
** (moinfo.opdm_b), which we also symmetrize by computing Opq = 1/2
** (Opq + Oqp).  These matrices are later written to disk in dump()
** for subsequent backtransformation.  Note that the components of the
** 1pdm computed into the DIJ, Dij, DAB, Dab, DAI, Dai, DIA, and Dia
** matrices remain non-symmetric (e.g., DIJ neq DJI).
**
** This will not work at present for frozen orbitals!
**
** TDC, 1/03
*/

void sortone_UHF(const struct RHO_Params& rho_params) {
    dpdfile2 D;

    Slice aocc_slice(moinfo.frdocc, moinfo.frdocc + moinfo.aoccpi);
    Slice avir_slice(moinfo.frdocc + moinfo.aoccpi, moinfo.orbspi - moinfo.fruocc);
    Slice bocc_slice(moinfo.frdocc, moinfo.frdocc + moinfo.boccpi);
    Slice bvir_slice(moinfo.frdocc + moinfo.boccpi, moinfo.orbspi - moinfo.fruocc);

    Matrix O_a(moinfo.orbspi, moinfo.orbspi);
    Matrix O_b(moinfo.orbspi, moinfo.orbspi);

    /* Sort A components first */
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, rho_params.DIJ_lbl);
    Matrix temp(&D);
    O_a.set_block(aocc_slice, temp);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, rho_params.DAB_lbl);
    temp = Matrix(&D);
    O_a.set_block(avir_slice, temp);
    global_dpd_->file2_close(&D);

    /* Note that this component of the density is stored occ-vir */
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DAI_lbl);
    temp = Matrix(&D);
    temp = *temp.transpose();
    O_a.set_block(avir_slice, aocc_slice, temp);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, rho_params.DIA_lbl);
    temp = Matrix(&D);
    O_a.set_block(aocc_slice, avir_slice, temp);
    global_dpd_->file2_close(&D);

    /* Sort B components */
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, rho_params.Dij_lbl);
    temp = Matrix(&D);
    O_b.set_block(bocc_slice, temp);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, rho_params.Dab_lbl);
    temp = Matrix(&D);
    O_b.set_block(bvir_slice, temp);
    global_dpd_->file2_close(&D);

    /* Note that this component of the density is stored occ-vir */
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dai_lbl);
    temp = Matrix(&D);
    temp = *temp.transpose();
    O_b.set_block(bvir_slice, bocc_slice, temp);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, rho_params.Dia_lbl);
    temp = Matrix(&D);
    O_b.set_block(bocc_slice, bvir_slice, temp);
    global_dpd_->file2_close(&D);

    O_a.hermitivitize();
    O_b.hermitivitize();
    moinfo.opdm_a = O_a;
    moinfo.opdm_b = O_b;
}
}
}  // namespace psi
