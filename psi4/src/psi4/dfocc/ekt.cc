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

// Latest revision on April 38, 2013.
#include <cstdio>
#include <fstream>
#include <cmath>
#include "ekt.h"

namespace psi {
namespace dfoccwave {

/********************************************************************************************/
/************************** 1d array ********************************************************/
/********************************************************************************************/
Ektip::Ektip(std::string name, int nocc, int norb, const SharedTensor2d& GFock, const SharedTensor2d& Gamma,
             double scale_gf, double scale_ps) {
    name_ = name;
    nocc_ = nocc;
    norb_ = norb;
    scale_ = scale_gf;
    scale2_ = scale_ps;
    nvir_ = norb_ - nocc_;
    cutoff_ = 1.0E-10;

    // malloc
    GF_ = std::make_shared<Tensor2d>("MO-basis GFM", norb_, norb_);
    GF_->copy(GFock);
    G1_ = std::make_shared<Tensor2d>("MO-basis OPDM", norb_, norb_);
    G1_->copy(Gamma);
    Uvec_ = std::make_shared<Tensor2d>("Uvec", norb_, norb_);
    Uvecp_ = std::make_shared<Tensor2d>("Uvec Prime", norb_, norb_);
    G1half_ = std::make_shared<Tensor2d>("G1^1/2", norb_, norb_);
    temp_ = std::make_shared<Tensor2d>("Temp", norb_, norb_);
    GFp_ = std::make_shared<Tensor2d>("GFM Prime", norb_, norb_);
    PS_ = std::make_shared<Tensor2d>("Pole Strength", norb_, norb_);
    GCt_ = std::make_shared<Tensor2d>("Alpha C'*gamma ", norb_, norb_);

    // 1D
    eocc_ = std::make_shared<Tensor1d>("epsilon <I|J>", nocc_);
    eorb_ = std::make_shared<Tensor1d>("epsilon <P|Q>", norb_);
    diagG1_ = std::make_shared<Tensor1d>("Diag G1", norb_);
    ps_vec_ = std::make_shared<Tensor1d>("pole strength vector", norb_);
    ps_occ_ = std::make_shared<Tensor1d>("occupied pole strength vector", nocc_);

    // compute ekt
    compute_ektip();

}  //

Ektip::~Ektip() {
    GF_.reset();
    G1_.reset();
    Uvec_.reset();
    Uvecp_.reset();
    G1half_.reset();
    temp_.reset();
    GFp_.reset();
    PS_.reset();
    GCt_.reset();

    eocc_.reset();
    eorb_.reset();
    diagG1_.reset();
    ps_vec_.reset();
    ps_occ_.reset();
}  //

void Ektip::compute_ektip() {
    // scale
    GF_->scale(scale_);

    // Make sure GFM is symmetric
    GF_copy_ = std::make_shared<Tensor2d>("MO-basis GFM", norb_, norb_);
    GF_copy_->copy(GF_);
    GFt_ = std::make_shared<Tensor2d>("MO-basis GFM", norb_, norb_);
    GFt_->trans(GF_);
    GF_copy_->add(GFt_);
    GF_copy_->scale(0.5);
    GF_->copy(GF_copy_);
    GFt_.reset();
    GF_copy_.reset();

    // Symm OPDM
    G1_copy_ = std::make_shared<Tensor2d>("MO-basis OPDM", norb_, norb_);
    G1_copy_->copy(G1_);
    G1t_ = std::make_shared<Tensor2d>("MO-basis OPDM", norb_, norb_);
    G1t_->trans(G1_);
    G1_copy_->add(G1t_);
    G1_copy_->scale(0.5);
    G1_->copy(G1_copy_);
    G1t_.reset();
    G1_copy_.reset();

    // Diagonalize OPDM
    G1_->diagonalize(Uvec_, diagG1_, cutoff_);

    // Make sure all eigenvalues are positive
    for (int i = 0; i < norb_; ++i) {
        if (diagG1_->get(i) < 0.0) diagG1_->set(i, -1.0 * diagG1_->get(i));
    }

    // Form g^(-1/2)
    for (int i = 0; i < norb_; ++i) {
        diagG1_->set(i, 1 / std::sqrt(diagG1_->get(i)));
    }

    for (int i = 0; i < norb_; ++i) {
        G1half_->set(i, i, diagG1_->get(i));
    }

    temp_->gemm(false, true, G1half_, Uvec_, 1.0, 0.0);
    G1half_->gemm(false, false, Uvec_, temp_, 1.0, 0.0);

    // Build GFock prime matrix
    temp_->gemm(true, false, G1half_, GF_, 1.0, 0.0);
    GFp_->gemm(false, false, temp_, G1half_, 1.0, 0.0);

    // Diagonalize GFock to get orbital energies
    GFp_->diagonalize(Uvecp_, eorb_, cutoff_);
    Uvec_->gemm(false, false, G1half_, Uvecp_, 1.0, 0.0);

    // Pole strength
    temp_->gemm(false, false, G1_, Uvec_, 1.0, 0.0);
    GCt_->trans(temp_);
    PS_->gemm(false, false, GCt_, temp_, 1.0, 0.0);
    // scale for RHF ref where G_pq = G_PQ + G_pq
    PS_->scale(scale2_);
    for (int i = 0; i < norb_; ++i) {
        ps_vec_->set(i, PS_->get(i, i));
    }

    // Sort pole strength
    // Sort to descending order
    for (int i = 0; i < norb_; ++i) {
        for (int j = norb_ - 1; j > i; --j) {
            if (ps_vec_->get(j - 1) < ps_vec_->get(j)) {
                double dum = eorb_->get(j - 1);
                eorb_->set(j - 1, eorb_->get(j));
                eorb_->set(j, dum);

                double dum2 = ps_vec_->get(j - 1);
                ps_vec_->set(j - 1, ps_vec_->get(j));
                ps_vec_->set(j, dum2);
            }
        }
    }

    // Re-Sort occupied orbitals to energy order
    // Copy
    for (int i = 0; i < nocc_; ++i) {
        eocc_->set(i, eorb_->get(i));
        ps_occ_->set(i, ps_vec_->get(i));
    }

    // Sort to ascending order
    for (int i = 0; i < nocc_; ++i) {
        for (int j = nocc_ - 1; j > i; --j) {
            if (eocc_->get(j - 1) > eocc_->get(j)) {
                double dum = eocc_->get(j - 1);
                eocc_->set(j - 1, eocc_->get(j));
                eocc_->set(j, dum);

                double dum2 = ps_occ_->get(j - 1);
                ps_occ_->set(j - 1, ps_occ_->get(j));
                ps_occ_->set(j, dum2);
            }
        }
    }

}  //

}  // namespace dfoccwave
}  // namespace psi
