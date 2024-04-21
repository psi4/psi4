/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "thc_eri.h"
#include "basisset.h"
#include "mintshelper.h"
#include "integral.h"

#include "psi4/libqt/qt.h"
#include "psi4/libfock/cubature.h"
#include "psi4/lib3index/dftensor.h"

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

THC_Computer::THC_Computer(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> primary, Options& options) :
    molecule_(molecule), primary_(primary), options_(options) {}

THC_Computer::~THC_Computer() {}

LS_THC_Computer::LS_THC_Computer(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> primary, 
                                    std::shared_ptr<BasisSet> auxiliary, Options& options) : THC_Computer(molecule, primary, options) {
    auxiliary_ = auxiliary;
}

LS_THC_Computer::~LS_THC_Computer() {}

SharedMatrix LS_THC_Computer::build_E_exact() {
    /**
     * Parrish et al. 2012 Procedure 2
    */

    throw PSIEXCEPTION("Error: Build E exact has not been implemented for LS-THC! Buy poor Andy a coffee!");
}

SharedMatrix LS_THC_Computer::build_E_df() {
    /**
     * Parrish et al. 2012 Procedure 3
    */
    
    size_t nbf = primary_->nbf();
    size_t naux = auxiliary_->nbf();
    size_t rank = x1_->nrow();

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    auto zero = BasisSet::zero_ao_basis_set();
    IntegralFactory factory(auxiliary_, zero, primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri_computers(nthreads);

    eri_computers[0] = std::shared_ptr<TwoBodyAOInt>(factory.eri());
#pragma omp parallel for
    for (int thread = 1; thread < nthreads; ++thread) {
        eri_computers[thread] = std::shared_ptr<TwoBodyAOInt>(eri_computers[0]->clone());
    }

    size_t nshell_aux = auxiliary_->nshell();
    size_t nshellpair = eri_computers[0]->shell_pairs().size();
    size_t nshelltriplet = nshell_aux * nshellpair;

    SharedMatrix E_temp = std::make_shared<Matrix>(rank, naux);

#pragma omp parallel for
    for (size_t MNP = 0; MNP < nshelltriplet; ++MNP) {

        size_t MN = MNP % nshellpair;
        size_t P = MNP / nshellpair;
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        auto bra = eri_computers[thread]->shell_pairs()[MN];
        size_t M = bra.first;
        size_t N = bra.second;

        int np = auxiliary_->shell(P).nfunction();
        int pstart = auxiliary_->shell(P).function_index();
        int nm = primary_->shell(M).nfunction();
        int mstart = primary_->shell(M).function_index();
        int nn = primary_->shell(N).nfunction();
        int nstart = primary_->shell(N).function_index();

        // TODO: Implement screening to make this more efficient
        eri_computers[thread]->compute_shell(P, 0, M, N);
        const auto &buffer = eri_computers[thread]->buffers()[0];

        double prefactor = (M == N) ? 1.0 : 2.0;

        for (size_t r = 0; r < rank; ++r) {
            for (size_t p = pstart, index = 0; p < pstart + np; ++p) {
                for (size_t m = mstart; m < mstart + nm; ++m) {
                    for (size_t n = nstart; n < nstart + nn; ++n, ++index) {
#pragma omp atomic
                        (*E_temp)(r, p) += prefactor * buffer[index] * (*x1_)(r, m) * (*x1_)(r, n);
                    } // end n
                } // end m
            } // end p
        } // end r
    } // end MNP

    FittingMetric J_metric_obj(auxiliary_, true);
    J_metric_obj.form_fitting_metric();
    auto J_metric = J_metric_obj.get_metric();
    
    int nremoved = 0;
    auto Jinv = J_metric->pseudoinverse(1.0e-14, nremoved);

    return linalg::triplet(E_temp, Jinv, E_temp, false, false, true);
}

void LS_THC_Computer::compute_thc_factorization() {
    /**
     * Parrish et al. 2012 Procedure 1
    */

    size_t nbf = primary_->nbf();

    timer_on("LS-THC: Build Grid");

    DFTGrid grid(molecule_, primary_, options_);
    size_t npoints = grid.npoints();

    x1_ = std::make_shared<Matrix>(npoints, nbf);

    size_t point_idx = 0;
    for (auto& block : grid.blocks()) {
        auto w = block->w();
        auto x = block->x();
        auto y = block->y();
        auto z = block->z();

#pragma omp parallel for
        for (size_t p = 0; p < block->npoints(); ++p) {
            primary_->compute_phi(&(*x1_)(point_idx + p, 0), x[p], y[p], z[p]);
            x1_->scale_row(0, point_idx + p, std::pow(w[p], 0.25));
        }

        point_idx += block->npoints();
    }

    x2_ = x1_;
    x3_ = x1_;
    x4_ = x1_;

    timer_off("LS-THC: Build Grid");
    
    timer_on("LS-THC: Form E");

    SharedMatrix E_PQ = options_.get_bool("LS_THC_DF") ? build_E_df() : build_E_exact();

    timer_off("LS-THC: Form E");

    timer_on("LS-THC: Form S");

    SharedMatrix S_Qq = std::make_shared<Matrix>(npoints, npoints);
    SharedMatrix S_temp = linalg::doublet(x1_, x1_, false, true);

#pragma omp parallel for
    for (size_t p = 0; p < npoints; ++p) {
        for (size_t q = 0; q < npoints; ++q) {
            (*S_Qq)(p, q) = (*S_temp)(p, q) * (*S_temp)(p, q);
        }
    }

    int nremoved = 0;
    SharedMatrix S_Qq_inv = S_Qq->pseudoinverse(1.0e-10, nremoved);

    timer_off("LS-THC: Form S");

    timer_on("LS-THC: Form Z");

    Z_PQ_ = linalg::triplet(S_Qq_inv, E_PQ, S_Qq_inv);

    timer_off("LS-THC: Form Z");
}

}