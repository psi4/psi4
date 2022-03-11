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

#ifndef _psi_src_lib_libmints_potentialint_h_
#define _psi_src_lib_libmints_potentialint_h_

#include "psi4/libmints/gshell.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/potential.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include "libint2/engine.h"

namespace psi {

class GaussianShell;
class SphericalTransform;
class BasisSet;

/**
 * This is a modifed version of the potential integral code that can parallelize loops over external points
 * and can inline code to process integrals into the inner loops, via functors.
 *
 * By defining the compute function of integral to be a template class, we can write classes (functors)
 * that will be inlined into the innermost loops, allowing us to do different tasks without re-writing
 * the code or having to make function calls. (AS)
 *
 * NB: This code must be specified in the .h file in order for the compiler to properly in-line the functors. (TDC)
 */
class PCMPotentialInt : public PotentialInt {

    std::vector<std::unique_ptr<libint2::Engine>> engines_;
   public:
    PCMPotentialInt(std::vector<SphericalTransform> &, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                    int deriv = 0);
    /// Drives the loops over all shell pairs, to compute integrals
    template <typename PCMPotentialIntFunctor>
    void compute(PCMPotentialIntFunctor &functor);
};

template <typename PCMPotentialIntFunctor>
void PCMPotentialInt::compute(PCMPotentialIntFunctor &functor) {
    bool bs1_equiv_bs2 = bs1_ == bs2_;
    size_t npoints = Zxyz_.size();
    size_t nthreads = engines_.size();
    functor.initialize(nthreads);
#pragma omp parallel for schedule(static) num_threads(nthreads)
    for (int point = 0; point < npoints; ++point) {
        size_t rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        engines_[rank]->set_params(std::vector<std::pair<double, std::array<double, 3>>>{Zxyz_[point]});
        for (const auto &pair : shellpairs_) {
            int p1 = pair.first;
            int p2 = pair.second;

            int ni = bs1_->shell(p1).nfunction();
            int nj = bs2_->shell(p2).nfunction();
            int i_offset = bs1_->shell_to_basis_function(p1);
            int j_offset = bs2_->shell_to_basis_function(p2);

            // Compute the shell
            engines_[rank]->compute(bs1_->l2_shell(p1), bs2_->l2_shell(p2));
            // For each integral that we got put in its contribution
            const double *location = engines_[rank]->results()[0];
            for (int p = 0; p < ni; ++p) {
                for (int q = 0; q < nj; ++q) {
                    functor(p + i_offset, q + j_offset, point, *location, rank);
                    if (bs1_equiv_bs2 && p1 != p2) {
                        functor(q + j_offset, p + i_offset, point, *location, rank);
                    }
                    ++location;
                }
            }
        }
    }
    functor.finalize(nthreads);
}

class PrintIntegralsFunctor {
   public:
    void initialize(int num_threads) {}
    void finalize(int num_threads) {}
    /**
     * A functor, to be used with PCMPotentialInt, that just prints the integrals out for debugging
     */
    void operator()(int bf1, int bf2, int center, double integral, int thread) {
        outfile->Printf("bf1: %3d bf2 %3d center (%5d) integral %16.10f\n", bf1, bf2, center, integral);
    }
};

class ContractOverDensityFunctor {
    /**
     * A functor, to be used with PCMPotentialInt, that just contracts potential integrals and the
     * density matrix, over the basis function indices, giving the charge expectation value.
     */
   protected:
    /// Pointer to the density matrix.
    double **pD_;
    /// The array of charges
    double *charges_;

   public:
    void initialize(int num_threads) {}
    void finalize(int num_threads) {}
    ContractOverDensityFunctor(size_t /*ncenters*/, double *charges, SharedMatrix D)
        : pD_(D->pointer()), charges_(charges) {}
    void operator()(int bf1, int bf2, int center, double integral, int thread) { charges_[center] += pD_[bf1][bf2] * integral; }
};

class ContractOverChargesFunctor {
    /**
     * A functor, to be used with PCMPotentialInt, that just contracts potential integrals over charges,
     * leaving a contribution to the Fock matrix
     */
   protected:
    /// Smart pointer to the Fock matrix contribution
    SharedMatrix F_;
    /// The array of charges
    const double *charges_;
    /// A copy of the Fock matrix for each thread to accumulate into
    std::vector<SharedMatrix> tempFock;

   public:
    ContractOverChargesFunctor(const double *charges, SharedMatrix F) : F_(F), charges_(charges) {
        if (F->rowdim() != F->coldim()) throw PSIEXCEPTION("Invalid Fock matrix in ContractOverCharges");
    }
    void initialize(int num_threads) {
        tempFock.clear();
        F_->zero();
        for (int thread = 0; thread < num_threads; ++thread) {
            tempFock.push_back(F_->clone());
        }
    }
    void finalize(int num_threads) {
        for (int thread = 0; thread < num_threads; ++thread) {
            F_->add(tempFock[thread]);
        }
    }
    void operator()(int bf1, int bf2, int center, double integral, int thread) {
        double **pF = tempFock[thread]->pointer();
        pF[bf1][bf2] += integral * charges_[center];
    }
};

}  // namespace psi
#endif
