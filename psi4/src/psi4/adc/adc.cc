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

#include "psi4/psi4-dec.h"

#include "psi4/liboptions/liboptions.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/adc/adc.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace adc {

ADCWfn::ADCWfn(SharedWavefunction ref_wfn, Options& options) : Wavefunction(options) {
    shallow_copy(ref_wfn);
    reference_wavefunction_ = ref_wfn;
    module_ = "adc";

    std::vector<std::string> irreps_ = molecule_->irrep_labels();
    aoccpi_ = nalphapi_ - frzcpi_;
    boccpi_ = nbetapi_ - frzcpi_;
    avirpi_ = nmopi_ - frzvpi_ - nalphapi_;
    bvirpi_ = nmopi_ - frzvpi_ - nbetapi_;

    int naocc = aoccpi_.sum();
    int nbocc = boccpi_.sum();
    int navir = avirpi_.sum();
    int nbvir = bvirpi_.sum();

    int aoccount = 0, boccount = 0, avircount = 0, bvircount = 0;

    aocce_ = new double[naocc];
    bocce_ = new double[nbocc];
    avire_ = new double[navir];
    bvire_ = new double[nbvir];

    nopen_ = soccpi().sum();
    if (!psio_) {
        throw PSIEXCEPTION("The wavefunction passed in lacks a PSIO object, crashing ADC. See GitHub issue #1851.");
    }
    if (nopen_) throw PSIEXCEPTION("Open-shell calculation has not been implemented yet!");

    aoccount = 0, boccount = 0, avircount = 0, bvircount = 0;
    for (int h = 0; h < nirrep_; h++) {
        for (int a = frzcpi_[h]; a < nalphapi_[h]; a++) aocce_[aoccount++] = epsilon_a_->get(h, a);

        for (int b = frzcpi_[h]; b < nbetapi_[h]; b++) bocce_[boccount++] = epsilon_b_->get(h, b);

        for (int a = nalphapi_[h]; a < nmopi_[h] - frzvpi_[h]; a++)
            avire_[avircount++] = epsilon_a_->get(h, a);

        for (int b = nbetapi_[h]; b < nmopi_[h] - frzvpi_[h]; b++) bvire_[bvircount++] = epsilon_b_->get(h, b);
    }

    const auto docc = doccpi();
    const auto socc = soccpi();
    outfile->Printf("\n\n\tIrrep  Core  Docc  Socc  aOcc  aVir  bOcc  bVir  FVir\n");
    outfile->Printf("\t*****************************************************\n");
    for (int h = 0; h < nirrep_; h++) {
        outfile->Printf("\t %3s   %3d   %3d   %3d   %3d   %3d   %3d   %3d   %3d\n", irreps_[h].c_str(), frzcpi_[h],
                        docc[h], socc[h], aoccpi_[h], avirpi_[h], boccpi_[h], bvirpi_[h], frzvpi_[h]);
    }
    outfile->Printf("\t*****************************************************\n\n");

    conv_ = options_.get_double("NEWTON_CONVERGENCE");
    norm_tol_ = options_.get_double("NORM_TOLERANCE");
    pole_max_ = options_.get_int("POLE_MAXITER");
    sem_max_ = options_.get_int("SEM_MAXITER");
    num_amps_ = options_.get_int("NUM_AMPS_PRINT");

    if (options_["ROOTS_PER_IRREP"].size() > 0) {
        int i = options_["ROOTS_PER_IRREP"].size();

        if (i != nirrep_) {
            outfile->Printf("dim of states_per_irrep vector must be %d\n", nirrep_);
            throw PsiException("adc input comparison error ROOTS_PER_IRREP and nirrep_", __FILE__, __LINE__);
        }
        rpi_ = Dimension(options_.get_int_vector("ROOTS_PER_IRREP"));
    } else {
        rpi_ = Dimension(std::vector<int>(nirrep_, 1));
    }

    // Setting up dimensions for each irrep block and totoal dimension of S manifold.
    nxs_ = 0;
    nxspi_ = new int[nirrep_];
    poles_ = (struct pole**)malloc(nirrep_ * sizeof(struct pole*));
    for (int h = 0; h < nirrep_; h++) {
        nxspi_[h] = 0;
        poles_[h] = (struct pole*)malloc(rpi_[h] * sizeof(struct pole));
        for (int Go = 0; Go < nirrep_; Go++) {
            nxspi_[h] += aoccpi_[Go] * avirpi_[Go ^ h];
        }
        nxs_ += nxspi_[h];
    }

    outfile->Printf("\t==> Input Parameters <==\n");
    outfile->Printf("\tNEWTON_CONV = %3g, NORM_TOL = %3g\n", conv_, norm_tol_);
    outfile->Printf("\tPOLE_MAX    = %3d, SEM_MAX  = %3d\n\n", pole_max_, sem_max_);

    outfile->Printf("\tNXS           = %d\n", nxs_);
    //    outfile->Printf( "\tIRREP_XYZ     = [");
    //    for(int i = 0;i < 3;i++)
    //        outfile->Printf( " %3s ", irreps_[irrep_axis_[i]]);
    //    outfile->Printf( "]\n");
    outfile->Printf("\tNXS_PER_IRREP = [");
    for (int i = 0; i < nirrep_; i++) {
        outfile->Printf(" %d ", nxspi_[i]);
    }
    outfile->Printf("]\n");

    if (DEBUG_) {
        outfile->Printf("Debagging mode...\n");
        outfile->Printf("\tNMO   = %3d, NXS = %3d\n", nmo_, nxs_);
        outfile->Printf("\tNOPEN = %3d\n", nopen_);
        outfile->Printf("\tROOTS_PER_IRREP = [");
        for (int i = 0; i < nirrep_; i++) outfile->Printf("%3d", rpi_[i]);
        outfile->Printf(" ]\n");
    }
}

ADCWfn::~ADCWfn() {}

void ADCWfn::release_mem() {
    free(poles_);
    delete _ints;
    delete[] aocce_;
    delete[] avire_;
    delete[] bocce_;
    delete[] bvire_;

    // omega_guess_.reset();
}
}
}  // End Namespaces
