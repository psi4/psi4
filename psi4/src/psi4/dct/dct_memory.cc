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

#include "dct.h"
#include "psi4/psifiles.h"
#include <vector>
#include <cmath>
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/libdpd/dpd.h"

namespace psi {
namespace dct {

/**
 * Reads the orbital information that can be determined before the SCF procedure
 * and initializes SO matrices.
 */
void DCTSolver::init() {
    nso_ = reference_wavefunction_->nso();
    nirrep_ = reference_wavefunction_->nirrep();
    nmo_ = reference_wavefunction_->nmo();
    enuc_ = reference_wavefunction_->molecule()->nuclear_repulsion_energy(
        reference_wavefunction_->get_dipole_field_strength());
    scf_energy_ = reference_wavefunction_->energy();
    ntriso_ = nso_ * (nso_ + 1) / 2;
    frzcpi_ = reference_wavefunction_->frzcpi();
    frzvpi_ = reference_wavefunction_->frzvpi();
    nmopi_ = reference_wavefunction_->nmopi();
    nsopi_ = reference_wavefunction_->nsopi();
    nalphapi_ = naoccpi_ = reference_wavefunction_->nalphapi();
    nbetapi_ = nboccpi_ = reference_wavefunction_->nbetapi();
    navirpi_ = (nmopi_ - naoccpi_ - frzvpi_);
    nbvirpi_ = (nmopi_ - nboccpi_ - frzvpi_);
    nalpha_ = naoccpi_.sum();
    nbeta_ = nboccpi_.sum();
    navir_ = navirpi_.sum();
    nbvir_ = nbvirpi_.sum();

    auto zero = Dimension(nirrep_);
    slices_ = {
        {"SO",  Slice(zero, nsopi_)},
        {"MO",  Slice(zero, nmopi_)},
        {"ACTIVE_OCC_A",  Slice(frzcpi_, nalphapi_)},
        {"ACTIVE_OCC_B",  Slice(frzcpi_, nbetapi_)},
        {"ACTIVE_VIR_A",  Slice(nalphapi_, nalphapi_ + navirpi_)},
        {"ACTIVE_VIR_B",  Slice(nbetapi_, nbetapi_ + nbvirpi_)}
    };

    aocc_c_ = std::make_shared<Matrix>("Alpha Occupied MO Coefficients", nirrep_, nsopi_, naoccpi_);
    bocc_c_ = std::make_shared<Matrix>("Beta Occupied MO Coefficients", nirrep_, nsopi_, nboccpi_);
    avir_c_ = std::make_shared<Matrix>("Alpha Virtual MO Coefficients", nirrep_, nsopi_, navirpi_);
    bvir_c_ = std::make_shared<Matrix>("Beta Virtual MO Coefficients", nirrep_, nsopi_, nbvirpi_);
    scf_error_a_ = std::make_shared<Matrix>("Alpha SCF Error Vector", nirrep_, nsopi_, nsopi_);
    scf_error_b_ = std::make_shared<Matrix>("Beta SCF Error Vector", nirrep_, nsopi_, nsopi_);
    Fa_ = reference_wavefunction_->Fa()->clone();
    Fb_ = reference_wavefunction_->Fb()->clone();
    Ftilde_a_ = std::make_shared<Matrix>("Alpha MO Ftilde Matrix", nirrep_, nmopi_, nmopi_);
    Ftilde_b_ = std::make_shared<Matrix>("Beta MO Ftilde Matrix", nirrep_, nmopi_, nmopi_);
    Ca_ = std::make_shared<Matrix>("Alpha MO Coefficients", nirrep_, nsopi_, nsopi_);
    Cb_ = std::make_shared<Matrix>("Beta MO Coefficients", nirrep_, nsopi_, nsopi_);
    moFa_ = std::make_shared<Matrix>("Alpha MO Fock Matrix", nirrep_, nmopi_, nmopi_);
    moFb_ = std::make_shared<Matrix>("Beta MO Fock Matrix", nirrep_, nmopi_, nmopi_);
    old_ca_ = std::make_shared<Matrix>("Old Alpha MO Coefficients", nirrep_, nsopi_, nsopi_);
    old_cb_ = std::make_shared<Matrix>("Old Beta MO Coefficients", nirrep_, nsopi_, nsopi_);
    kappa_so_a_ = std::make_shared<Matrix>("Alpha Kappa Matrix", nirrep_, nsopi_, nsopi_);
    kappa_so_b_ = std::make_shared<Matrix>("Beta Kappa Matrix", nirrep_, nsopi_, nsopi_);
    ao_s_ = std::make_shared<Matrix>("SO Basis Overlap Integrals", nirrep_, nsopi_, nsopi_);
    so_h_ = Matrix("SO basis one-electron integrals", nirrep_, nsopi_, nsopi_);
    s_half_inv_ = std::make_shared<Matrix>("SO Basis Inverse Square Root Overlap Matrix", nirrep_, nsopi_, nsopi_);
    epsilon_a_ = std::make_shared<Vector>(nirrep_, nsopi_);
    epsilon_b_ = std::make_shared<Vector>(nirrep_, nsopi_);
    kappa_mo_a_ = std::make_shared<Matrix>("MO basis Kappa (Alpha)", nirrep_, nmopi_, nmopi_);
    kappa_mo_b_ = std::make_shared<Matrix>("MO basis Kappa (Beta)", nirrep_, nmopi_, nmopi_);
    tau_so_a_ = std::make_shared<Matrix>("Alpha Tau Matrix", nirrep_, nsopi_, nsopi_);
    tau_so_b_ = std::make_shared<Matrix>("Beta Tau Matrix", nirrep_, nsopi_, nsopi_);
    aocc_tau_ = Matrix("MO basis Tau (Alpha Occupied)", nirrep_, naoccpi_, naoccpi_);
    bocc_tau_ = Matrix("MO basis Tau (Beta Occupied)", nirrep_, nboccpi_, nboccpi_);
    avir_tau_ = Matrix("MO basis Tau (Alpha Virtual)", nirrep_, navirpi_, navirpi_);
    bvir_tau_ = Matrix("MO basis Tau (Beta Virtual)", nirrep_, nbvirpi_, nbvirpi_);

    // Compute MO offsets
    aocc_off_ = new int[nirrep_];
    avir_off_ = new int[nirrep_];
    double ocount = naoccpi_[0];
    double vcount = navirpi_[0];
    aocc_off_[0] = 0;
    avir_off_[0] = 0;
    for (int h = 1; h < nirrep_; h++) {
        aocc_off_[h] = ocount;
        ocount += naoccpi_[h];
        avir_off_[h] = vcount;
        vcount += navirpi_[h];
    }

    bocc_off_ = new int[nirrep_];
    bvir_off_ = new int[nirrep_];
    ocount = nboccpi_[0];
    vcount = nbvirpi_[0];
    bocc_off_[0] = 0;
    bvir_off_[0] = 0;
    for (int h = 1; h < nirrep_; h++) {
        bocc_off_[h] = ocount;
        ocount += nboccpi_[h];
        bvir_off_[h] = vcount;
        vcount += nbvirpi_[h];
    }

    // Quadratically-convergent algorithm or orbital-optimized methods
    if (options_.get_str("ALGORITHM") == "QC" || orbital_optimized_) {
        orbital_gradient_a_ = std::make_shared<Matrix>("MO basis Orbital Gradient (Alpha)", nirrep_, nmopi_, nmopi_);
        orbital_gradient_b_ = std::make_shared<Matrix>("MO basis Orbital Gradient (Beta)", nirrep_, nmopi_, nmopi_);
        Xtotal_a_ = std::make_shared<Matrix>("Generator of the orbital rotations w.r.t. reference orbitals (Alpha)",
                                             nirrep_, nmopi_, nmopi_);
        Xtotal_b_ = std::make_shared<Matrix>("Generator of the orbital rotations w.r.t. reference orbitals (Beta)",
                                             nirrep_, nmopi_, nmopi_);
    }

    if (options_.get_str("ALGORITHM") == "QC") {
        // The number of IDPs is set to zero in the beginning
        nidp_ = 0;

        dim_orbitals_ = nalpha_ * navir_ + nbeta_ * nbvir_;
        dim_cumulant_ = (nalpha_ * (nalpha_ - 1) / 2) * (navir_ * (navir_ - 1) / 2);
        dim_cumulant_ += (nalpha_ * nbeta_) * (navir_ * nbvir_);
        dim_cumulant_ += (nbeta_ * (nbeta_ - 1) / 2) * (nbvir_ * (nbvir_ - 1) / 2);
        dim_ = dim_orbitals_ + dim_cumulant_;

        lookup_orbitals_ = std::vector<bool>(dim_orbitals_, false);
        lookup_cumulant_ = std::vector<bool>(dim_cumulant_, false);
    }

    // Fill up Kappa array
    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < naoccpi_[h]; ++i) {
            for (int j = 0; j < naoccpi_[h]; ++j) {
                kappa_mo_a_->set(h, i, j, (i == j ? 1.0 : 0.0));
            }
        }
        for (int i = 0; i < nboccpi_[h]; ++i) {
            for (int j = 0; j < nboccpi_[h]; ++j) {
                kappa_mo_b_->set(h, i, j, (i == j ? 1.0 : 0.0));
            }
        }
    }

    // Construct the Core Hamiltonian
    so_h_.add(*mintshelper()->so_kinetic()->clone());
    so_h_.add(*mintshelper()->so_potential()->clone());

    // Grab the SO overlap matrix from MintsHelper
    ao_s_->copy(mintshelper_->so_overlap());

    // Form S^(-1/2) matrix
    Matrix eigvec(nirrep_, nsopi_, nsopi_);
    Matrix eigtemp(nirrep_, nsopi_, nsopi_);
    Matrix eigtemp2(nirrep_, nsopi_, nsopi_);
    Vector eigval(nirrep_, nsopi_);
    ao_s_->diagonalize(eigvec, eigval);
    // Convert the eigenvales to 1/sqrt(eigenvalues)
    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < nsopi_[h]; ++i) {
            double scale = 1.0 / sqrt(eigval.get(h, i));
            eigval.set(h, i, scale);
        }
    }
    eigtemp2.set_diagonal(eigval);
    eigtemp.gemm(false, true, 1.0, eigtemp2, eigvec, 0.0);
    s_half_inv_->gemm(false, false, 1.0, eigvec, eigtemp, 0.0);
}

/**
 * Frees up the memory sequestered by the init_moinfo() and read_checkpoint() routines.
 */
void DCTSolver::finalize() {
    psio_->close(PSIF_DCT_DPD, 1);
    // We need to manually reset _ints to trigger freeing PSIF_DCT_DPD.
    // Having two at once is forbidden.
    _ints.reset();

    aocc_c_.reset();
    bocc_c_.reset();
    avir_c_.reset();
    bvir_c_.reset();
    scf_error_a_.reset();
    scf_error_b_.reset();
    Fa_.reset();
    Fb_.reset();
    old_ca_.reset();
    old_cb_.reset();
    kappa_so_a_.reset();
    kappa_so_b_.reset();
    ao_s_.reset();
    s_half_inv_.reset();
}

}  // namespace dct
}  // namespace psi
