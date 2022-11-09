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

#include "psi4/libmints/potentialint.h"
#include "psi4/libpsi4util/process.h"
#include "libint2/engine.h"

namespace psi {

PCMPotentialInt::PCMPotentialInt(std::vector<SphericalTransform>& trans, std::shared_ptr<BasisSet> bs1,
                                 std::shared_ptr<BasisSet> /* bs2 */, int /* deriv */)
    : PotentialInt(trans, bs1, bs1) {

    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());

    int nthreads = Process::environment.get_n_threads();
    for(int thread=0; thread < nthreads; ++thread) {
        engines_.push_back(std::make_unique<libint2::Engine>(libint2::Operator::nuclear, max_nprim, max_am, 0));
    }
}

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

template void PCMPotentialInt::compute(ContractOverChargesFunctor&);
template void PCMPotentialInt::compute(PrintIntegralsFunctor&);
template void PCMPotentialInt::compute(ContractOverDensityFunctor&);

}  // namespace psi
