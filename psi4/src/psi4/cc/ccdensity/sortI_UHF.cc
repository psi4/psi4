/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

/* SORTI_UHF(): Place all the components of the UHF Lagrangian into a
** large matrix, I (moinfo.I), which we also symmetrize by computing
** Ipq = 1/2 (Ipq + Iqp).  This matrix is later written to disk in
** dump() for subsequent backtransformation.  Note that some of the
** components of the Lagrangian computed into the IIJ, Iij, IIA, and
** Iia matrices remain non-symmetric (e.g., IIJ neq IJI).  I re-used
** my sortone.c code here, so don't let some of the variable names
** confuse you. */

void sortI_UHF(Wavefunction& wfn) {
    dpdfile2 D;

    Slice aocc_slice(moinfo.frdocc, moinfo.frdocc + moinfo.aoccpi);
    Slice avir_slice(moinfo.frdocc + moinfo.aoccpi, moinfo.orbspi - moinfo.fruocc);
    Slice bocc_slice(moinfo.frdocc, moinfo.frdocc + moinfo.boccpi);
    Slice bvir_slice(moinfo.frdocc + moinfo.boccpi, moinfo.orbspi - moinfo.fruocc);

    auto O = std::make_shared<Matrix>("Lagrangian matrix", moinfo.orbspi, moinfo.orbspi);
    Matrix O_b(moinfo.orbspi, moinfo.orbspi);

    /* Sort alpha components first */
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "I(I,J)");
    Matrix temp(&D);
    O->set_block(aocc_slice, temp);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "I'AB");
    temp = Matrix(&D);
    O->set_block(avir_slice, temp);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "I(I,A)");
    temp = Matrix(&D);
    O->set_block(aocc_slice, avir_slice, temp);
    temp = *temp.transpose();
    O->set_block(avir_slice, aocc_slice, temp);
    global_dpd_->file2_close(&D);

    /* Sort beta components */
    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 2, "I(i,j)");
    temp = Matrix(&D);
    O_b.set_block(bocc_slice, temp);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 3, 3, "I'ab");
    temp = Matrix(&D);
    O_b.set_block(bvir_slice, temp);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 2, 3, "I(i,a)");
    temp = Matrix(&D);
    O_b.set_block(bocc_slice, bvir_slice, temp);
    temp = *temp.transpose();
    O_b.set_block(bvir_slice, bocc_slice, temp);
    global_dpd_->file2_close(&D);

    O->hermitivitize();
    O_b.hermitivitize();

    O->back_transform(wfn.Ca());
    O_b.back_transform(wfn.Cb());
    O->add(O_b);
    O->scale(-1.0);
    wfn.set_lagrangian(O);
}

}  // namespace ccdensity
}  // namespace psi
