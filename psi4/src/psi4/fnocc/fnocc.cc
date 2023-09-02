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

#include "ccsd.h"
#include "frozen_natural_orbitals.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libqt/qt.h"

namespace psi {
namespace fnocc {

SharedWavefunction fnocc(SharedWavefunction ref_wfn, Options &options) {
    std::shared_ptr<Wavefunction> wfn;

    if (!options.get_bool("DFCC")) {
        // frozen natural orbital ccsd(t)
        if (options.get_bool("NAT_ORBS")) {
            auto fno = std::make_shared<FrozenNO>(ref_wfn, options);
            fno->ComputeNaturalOrbitals();
            wfn = (std::shared_ptr<Wavefunction>)fno;

        } else {
            wfn = ref_wfn;
        }

        // transform integrals
        tstart();
        outfile->Printf("        ==> Transform all two-electron integrals <==\n");
        outfile->Printf("\n");

        std::vector<std::shared_ptr<MOSpace> > spaces;
        spaces.push_back(MOSpace::all);
        std::shared_ptr<IntegralTransform> ints = std::make_shared<IntegralTransform>(
            wfn, spaces, IntegralTransform::TransformationType::Restricted, IntegralTransform::OutputType::IWLOnly,
            IntegralTransform::MOOrdering::QTOrder, IntegralTransform::FrozenOrbitals::OccAndVir, false);
        ints->set_dpd_id(0);
        ints->set_keep_iwl_so_ints(true);
        ints->set_keep_dpd_so_ints(true);
        timer_on("FNOCC: Initializing Integrals");
        ints->initialize();
        timer_off("FNOCC: Initializing Integrals");
        timer_on("FNOCC: Two Elec Int Trans.");
        ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
        timer_off("FNOCC: Two Elec Int Trans.");
        tstop();

        if (!options.get_bool("RUN_CEPA")) {
            auto ccsd = std::make_shared<CoupledCluster>(wfn, options);
            timer_on("FNOCC: CC energy");
            ccsd->compute_energy();
            timer_off("FNOCC: CC energy");
            return ccsd;
        } else {
            auto cepa = std::make_shared<CoupledPair>(wfn, options);
            timer_on("FNOCC: CEPA energy");
            cepa->compute_energy();
            timer_off("FNOCC: CEPA energy");
            return cepa;
        }

    } else {
        tstart();

        outfile->Printf("\n\n");
        outfile->Printf("        *******************************************************\n");
        outfile->Printf("        *                                                     *\n");
        outfile->Printf("        *                       DF-CCSD                       *\n");
        outfile->Printf("        *                 Density-fitted CCSD                 *\n");
        outfile->Printf("        *                                                     *\n");
        outfile->Printf("        *                   Eugene DePrince                   *\n");
        outfile->Printf("        *                                                     *\n");
        outfile->Printf("        *******************************************************\n");
        outfile->Printf("\n\n");

        // three-index integrals are generated/read by fno class
        auto fno = std::make_shared<DFFrozenNO>(ref_wfn, options);
        fno->ThreeIndexIntegrals();
        if (options.get_bool("NAT_ORBS")) {
            fno->ComputeNaturalOrbitals();
            wfn = (std::shared_ptr<Wavefunction>)fno;
        } else {
            wfn = ref_wfn;
        }
// ccsd(t)!

#ifdef GPUCC
        auto ccsd = std::make_shared<GPUDFCoupledCluster>(wfn, options);
        timer_on("FNOCC: CC energy");
        ccsd->compute_energy();
        timer_off("FNOCC: CC energy");
#else
        auto ccsd = std::make_shared<DFCoupledCluster>(wfn, options);
        timer_on("FNOCC: CC energy");
        ccsd->compute_energy();
        timer_off("FNOCC: CC energy");
#endif

        tstop();

        return ccsd;
    }

    // return wfn;
}  // end fnocc
}
}  // end namespaces
