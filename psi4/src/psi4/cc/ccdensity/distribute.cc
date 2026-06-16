/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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
#include "psi4/libiwl/iwl_reader.h"
#include "psi4/libiwl/iwl_writer.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#include "globals.h"

namespace psi {
namespace ccdensity {

void classify(int p, int q, int r, int s, double value, IWLWriter &ABuf, IWLWriter &BBuf, IWLWriter &CBuf,
              IWLWriter &DBuf, IWLWriter &EBuf, IWLWriter &FBuf);

void distribute() {
    double tolerance = params.tolerance;

    IWLWriter ABuf(_default_psio_lib_, 90, tolerance);
    IWLWriter BBuf(_default_psio_lib_, 91, tolerance);
    IWLWriter CBuf(_default_psio_lib_, 92, tolerance);
    IWLWriter DBuf(_default_psio_lib_, 93, tolerance);
    IWLWriter EBuf(_default_psio_lib_, 94, tolerance);
    IWLWriter FBuf(_default_psio_lib_, 95, tolerance);

    IWLReader eri(_default_psio_lib_, PSIF_MO_TEI);
    for (const auto &integral : eri) {
        /* Check integral into each class */
        classify(integral.p, integral.q, integral.r, integral.s, integral.value, ABuf, BBuf, CBuf, DBuf, EBuf,
                 FBuf);
    }
    // IWLReader/IWLWriter close (and the writers flush the last buffer) on scope exit.
}

}  // namespace ccdensity
}  // namespace psi
