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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "Params.h"
#include "psi4/cc/ccwave.h"
#include "psi4/libqt/qt.h"

namespace psi {
namespace ccenergy {

void CCEnergyWavefunction::cc2_t2_build() {
    DT2();

    if ((params_.ref == 0) || params_.t2_coupled) { /** RHF or ROHF with coupled T2's **/
        timer_on("fT2");
        cc2_faeT2();
        cc2_fmiT2();
        if (params_.print & 2) status("f -> T2", "outfile");
        timer_off("fT2");
    }

    timer_on("WmbijT2");
    cc2_WmbijT2();
    if (params_.print & 2) status("Wmbij -> T2", "outfile");
    timer_off("WmbijT2");

    timer_on("WabeiT2");
    cc2_WabeiT2();
    if (params_.print & 2) status("Wabei -> T2", "outfile");
    timer_off("WabeiT2");
}
}  // namespace ccenergy
}  // namespace psi
