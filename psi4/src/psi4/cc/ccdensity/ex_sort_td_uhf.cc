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
#include "psi4/liboptions/liboptions.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {

void ex_sort_td_uhf(char hand, int Tirrep) {
    dpdfile2 D;

    Matrix a_mat(moinfo.orbspi, moinfo.orbspi, Tirrep);
    Matrix b_mat(moinfo.orbspi, moinfo.orbspi, Tirrep);

    Slice aocc_slice(moinfo.frdocc, moinfo.frdocc + moinfo.aoccpi);
    Slice avir_slice(moinfo.frdocc + moinfo.aoccpi, moinfo.orbspi - moinfo.fruocc);
    Slice bocc_slice(moinfo.frdocc, moinfo.frdocc + moinfo.boccpi);
    Slice bvir_slice(moinfo.frdocc + moinfo.boccpi, moinfo.orbspi - moinfo.fruocc);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, Tirrep, 0, 0, "LTDIJ");
    Matrix temp_mat(&D);
    a_mat.set_block(aocc_slice, temp_mat);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, Tirrep, 1, 1, "LTDAB");
    temp_mat = Matrix(&D);
    a_mat.set_block(avir_slice, temp_mat);
    global_dpd_->file2_close(&D);

    /* Note that this component of the density is stored occ-vir */
    global_dpd_->file2_init(&D, PSIF_CC_TMP, Tirrep, 0, 1, "LTDAI");
    temp_mat = Matrix(&D);
    temp_mat = *temp_mat.transpose();
    a_mat.set_block(avir_slice, aocc_slice, temp_mat);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, Tirrep, 0, 1, "LTDIA");
    temp_mat = Matrix(&D);
    a_mat.set_block(aocc_slice, avir_slice, temp_mat);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, Tirrep, 2, 2, "LTDij");
    temp_mat = Matrix(&D);
    b_mat.set_block(bocc_slice, temp_mat);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, Tirrep, 3, 3, "LTDab");
    temp_mat = Matrix(&D);
    b_mat.set_block(bvir_slice, temp_mat);
    global_dpd_->file2_close(&D);

    /* Note that this component of the density is stored occ-vir */
    global_dpd_->file2_init(&D, PSIF_CC_TMP, Tirrep, 2, 3, "LTDai");
    temp_mat = Matrix(&D);
    temp_mat = *temp_mat.transpose();
    b_mat.set_block(bvir_slice, bocc_slice, temp_mat);
    global_dpd_->file2_close(&D);

    global_dpd_->file2_init(&D, PSIF_CC_TMP, Tirrep, 2, 3, "LTDia");
    temp_mat = Matrix(&D);
    b_mat.set_block(bocc_slice, bvir_slice, temp_mat);
    global_dpd_->file2_close(&D);

    if (hand == 'l') {
        moinfo.ltd_a_mat = std::move(a_mat);
        moinfo.ltd_b_mat = std::move(b_mat);
    } else if (hand == 'r') {
        moinfo.rtd_a_mat = std::move(a_mat);
        moinfo.rtd_b_mat = std::move(b_mat);
    } else
        throw PsiException("ccdensity: error", __FILE__, __LINE__);

    return;
}

}  // namespace ccdensity
}  // namespace psi
