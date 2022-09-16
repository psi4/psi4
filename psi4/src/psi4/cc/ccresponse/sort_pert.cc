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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstring>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

/* sort_pert(): Sorts the specified MO-basis one-electron property
** integrals into CC ordering for use in building the
** similarity-transformed integrals and certain components of the
** total linear response function.
**
** NB: Some integrals are antisymmetric (e.g. L or P integrals), and
** others are symmetric (e.g. Mu integrals), so we must be careful in
** this and subsequent routines.
**
** TDC, 10/05
*/

void write_blocks(const Matrix& mat) {
    Slice occ_slice(Dimension(moinfo.nirreps), moinfo.act_occpi);
    Slice vir_slice(moinfo.act_occpi, moinfo.act_pi);

    dpdfile2 f;

    auto block = mat.get_block(occ_slice);
    auto lbl = mat.name() + "_IJ";
    global_dpd_->file2_init(&f, PSIF_CC_OEI, mat.symmetry(), 0, 0, lbl);
    block->write_to_dpdfile2(&f);
    global_dpd_->file2_close(&f);

    block = mat.get_block(vir_slice);
    lbl = mat.name() + "_AB";
    global_dpd_->file2_init(&f, PSIF_CC_OEI, mat.symmetry(), 1, 1, lbl);
    block->write_to_dpdfile2(&f);
    global_dpd_->file2_close(&f);

    block = mat.get_block(occ_slice, vir_slice);
    lbl = mat.name() + "_IA";
    global_dpd_->file2_init(&f, PSIF_CC_OEI, mat.symmetry(), 0, 1, lbl);
    block->write_to_dpdfile2(&f);
    global_dpd_->file2_close(&f);
}

}  // namespace ccresponse
}  // namespace psi
