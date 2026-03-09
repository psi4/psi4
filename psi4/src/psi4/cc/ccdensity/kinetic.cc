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
    \brief Computes the kinetic energy and the virial ratio for CC wave functions.
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {

void kinetic(std::shared_ptr<Wavefunction> wfn) {
    /* RHF/ROHF only for now */
    if (params.ref == 2) return;

    /*** Transform the kinetic energy integrals to the MO basis ***/
    auto T = wfn->mintshelper()->so_kinetic();
    T->transform(moinfo.Ca);

    /*** Contract the correlated kinetic energy ***/
    auto tcorr = T->vector_dot(moinfo.opdm);

    /*** Recall the SCF kinetic energy ***/
    auto tref = wfn->scalar_variable("HF KINETIC ENERGY");

    /*** Compute the virial ratios ***/
    auto ttot = tcorr + tref;
    auto vtot = moinfo.eref + moinfo.ecc - ttot;
    auto vref = moinfo.eref - tref;
    auto vcorr = moinfo.ecc - tcorr;

    outfile->Printf("\n\tVirial Theorem Data:\n");
    outfile->Printf("\t--------------------\n");
    outfile->Printf("\tKinetic energy (ref)   = %20.15f\n", tref);
    outfile->Printf("\tKinetic energy (corr)  = %20.15f\n", tcorr);
    outfile->Printf("\tKinetic energy (total) = %20.15f\n", ttot);

    outfile->Printf("\t-V/T (ref)             = %20.15f\n", -vref / tref);
    outfile->Printf("\t-V/T (corr)            = %20.15f\n", -vcorr / tcorr);
    outfile->Printf("\t-V/T (total)           = %20.15f\n", -vtot / ttot);

    wfn->set_scalar_variable("CC CORRELATION KINETIC ENERGY", tcorr);
    wfn->set_scalar_variable("CC CORRELATION POTENTIAL ENERGY", vcorr);
    wfn->set_scalar_variable("CC CORRELATION VIRIAL RATIO", -vcorr / tcorr);
    wfn->set_scalar_variable("CC VIRIAL RATIO", -vtot / ttot);
}
}  // namespace ccdensity
}  // namespace psi
