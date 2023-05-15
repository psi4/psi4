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

#include "psi4/libqt/qt.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/integral.h"
#include "psi4/lib3index/dftensor.h"
#include "jk.h"
#include "SplitJK.h"

#include <unordered_set>
#include <vector>
#include <map>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {

DirectDFJ::DirectDFJ(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options) : SplitJK(primary, options), auxiliary_(auxiliary) {
    timer_on("DirectDFJ: Setup");

    // => General Setup <= //

    // thread count
    nthreads_ = 1;
#ifdef _OPENMP
    nthreads_ = Process::environment.get_n_threads();
#endif

    // pre-compute coulomb fitting metric
    timer_on("DirectDFJ: DIRECTDFJ Coulomb Metric");

    FittingMetric J_metric_obj(auxiliary_, true);
    J_metric_obj.form_fitting_metric();
    J_metric_ = J_metric_obj.get_metric();

    timer_off("DirectDFJ: DIRECTDFJ Coulomb Metric");

    timer_off("DirectDFJ: Setup");
}

DirectDFJ::~DirectDFJ() {}

size_t DirectDFJ::num_computed_shells() {
    return num_computed_shells_;
}

void DirectDFJ::print_header() const {
    if (print_) {
        outfile->Printf("\n");
        outfile->Printf("  ==> DF-DirJ: Integral-Direct Density-Fitted J <==\n\n");

        outfile->Printf("    J Screening Cutoff:%11.0E\n", cutoff_);
    }
}

// build the J matrix using Weigend's integral-direct density fitting algorithm
// algorithm is in Figure 1 of https://doi.org/10.1039/B204199P
void DirectDFJ::build_G_component(std::vector<std::shared_ptr<Matrix>>& D, std::vector<std::shared_ptr<Matrix>>& J,
    std::vector<std::shared_ptr<TwoBodyAOInt> >& eri_computers) {

    timer_on("Setup");

    // => Sizing <= //
    int njk = D.size();
    int nbf = primary_->nbf();
    int nshell = primary_->nshell();
    int nbf_aux = auxiliary_->nbf();
    int nshell_aux = auxiliary_->nshell();

    // benchmarking
    size_t nshellpair = eri_computers[0]->shell_pairs().size();
    size_t nshelltriplet = nshell_aux * nshellpair;
    size_t computed_triplets1 = 0, computed_triplets2 = 0;

    // screening threshold
    double thresh2 = options_.get_double("INTS_TOLERANCE") * options_.get_double("INTS_TOLERANCE");

    // per-thread G Vector buffers (for accumulating thread contributions to G)
    // G is the contraction of the density matrix with the 3-index ERIs
    std::vector<std::vector<SharedVector>> GT(njk, std::vector<SharedVector>(nthreads_));

    // H is the contraction of G with the inverse coulomb metric
    std::vector<SharedVector> H(njk);

    // per-thread J Matrix buffers (for accumulating thread contributions to J)
    std::vector<std::vector<SharedMatrix>> JT(njk, std::vector<SharedMatrix>(nthreads_));

    // initialize per-thread objects
    for(size_t jki = 0; jki < njk; jki++) {
        for(size_t thread = 0; thread < nthreads_; thread++) {
            JT[jki][thread] = std::make_shared<Matrix>(nbf, nbf);
            GT[jki][thread] = std::make_shared<Vector>(nbf_aux);
        }
        H[jki] = std::make_shared<Vector>(nbf_aux);
    }

    // diagonal shell maxima of J_metric_ for screening
    std::vector<double> J_metric_shell_diag(nshell_aux, 0.0);
    for (size_t s = 0; s < nshell_aux; s++) {
        int bf_start = auxiliary_->shell(s).function_index();
        int bf_end = bf_start + auxiliary_->shell(s).nfunction();
        for (size_t bf = bf_start; bf < bf_end; bf++) {
            J_metric_shell_diag[s] = std::max(J_metric_shell_diag[s], J_metric_->get(bf, bf));
        }
    }

    // shell maxima of D for screening
    Matrix Dshell(nshell, nshell);
    auto Dshellp = Dshell.pointer();

    for(size_t M = 0; M < nshell; M++) {
        int nm = primary_->shell(M).nfunction();
        int mstart = primary_->shell(M).function_index();
        for(size_t N = 0; N < nshell; N++) {
            int nn = primary_->shell(N).nfunction();
            int nstart = primary_->shell(N).function_index();
            for(size_t jki = 0; jki < njk; jki++) {
                auto Dp = D[jki]->pointer();
                for(size_t m = mstart; m < mstart + nm; m++) {
                    for(size_t n = nstart; n < nstart + nn; n++) {
                        Dshellp[M][N] = std::max(Dshellp[M][N], std::abs(Dp[m][n]));
                    }
                }
            }
        }
    }

    timer_off("Setup");

    //  => First Contraction <= //

    // contract D with three-index DF ERIs to get G:
    // G_{p} = D_{mn} * (mn|p)
    // G_{p} correlates to gamma_P in Figure 1 of Weigend's paper

    timer_on("ERI1");
#pragma omp parallel for schedule(guided) num_threads(nthreads_) reduction(+ : computed_triplets1)
    for (size_t MNP = 0; MNP < nshelltriplet; MNP++) {

        size_t MN = MNP % nshellpair;
        size_t P = MNP / nshellpair;
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        auto bra = eri_computers[rank]->shell_pairs()[MN];
        size_t M = bra.first;
        size_t N = bra.second;
        if(Dshellp[M][N] * Dshellp[M][N] * J_metric_shell_diag[P] * eri_computers[rank]->shell_pair_value(M,N) < thresh2) {
            continue;
        }
        computed_triplets1++;
        int np = auxiliary_->shell(P).nfunction();
        int pstart = auxiliary_->shell(P).function_index();
        int nm = primary_->shell(M).nfunction();
        int mstart = primary_->shell(M).function_index();
        int nn = primary_->shell(N).nfunction();
        int nstart = primary_->shell(N).function_index();
        eri_computers[rank]->compute_shell(P, 0, M, N);
        const auto & buffer = eri_computers[rank]->buffers()[0];

        for(size_t jki = 0; jki < njk; jki++) {

            auto GTp = GT[jki][rank]->pointer();
            auto Dp = D[jki]->pointer();

            for (int p = pstart, index = 0; p < pstart + np; p++) {
                for (int m = mstart; m < mstart + nm; m++) {
                    for (int n = nstart; n < nstart + nn; n++, index++) {
                        GTp[p] += buffer[index] * Dp[m][n];
                        if (N != M) GTp[p] += buffer[index] * Dp[n][m];
                    }
                }
            }

        }

    }

    timer_off("ERI1");

    //  => Second Contraction <= //

    //  linear solve for H:
    //  G_{p} = H_{q} (q|p)
    //  H_{p} correlates to gamma_Q in Figure 1 of Weigend's paper

    timer_on("Metric");

    std::vector<int> ipiv(nbf_aux);

    for(size_t jki = 0; jki < njk; jki++) {
        for(size_t thread = 0; thread < nthreads_; thread++) {
            H[jki]->add(*GT[jki][thread]);
        }
        C_DGESV(nbf_aux, 1, J_metric_->clone()->pointer()[0], nbf_aux, ipiv.data(), H[jki]->pointer(), nbf_aux);
    }


    // I believe C_DSYSV should be faster than C_GESV, but I've found the opposite to be true.
    // This performance issue should be investigated, but is not consequential here.
    // The cost of either linear solve is dwarfed by the actual integral computation.
    //
    //std::vector<double> work(3 * nbf_aux);
    //int errcode = C_DSYSV('U', nbf_aux, 1, J_metric_->clone()->pointer()[0], nbf_aux, ipiv.data(), H[jki]->pointer(), nbf_aux, work.data(), 3 * nbf_aux);

    // shell maxima of H for screening
    Vector H_shell_max(nshell_aux);
    auto H_shell_maxp = H_shell_max.pointer();

    for(size_t jki = 0; jki < njk; jki++) {

        auto Hp = H[jki]->pointer();

        for (int P = 0; P < nshell_aux; P++) {
            int np = auxiliary_->shell(P).nfunction();
            int pstart = auxiliary_->shell(P).function_index();
            for (int p = pstart; p < pstart + np; p++) {
                H_shell_maxp[P] = std::max(H_shell_maxp[P], std::abs(Hp[p]));
            }
        }

    }

    timer_off("Metric");

    //  => Third Contraction <= //

    // contract H with three-index DF ERIs to get J
    // J_{mn} = H_{p} (mn|p)
    // J_{mn} correlates to J_[uv] in Figure 1 of Weigend's paper

    timer_on("ERI2");

#pragma omp parallel for schedule(guided) num_threads(nthreads_) reduction(+ : computed_triplets2)
    for (size_t MNP = 0; MNP < nshelltriplet; MNP++) {

        size_t MN = MNP % nshellpair;
        size_t P = MNP / nshellpair;
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        auto bra = eri_computers[rank]->shell_pairs()[MN];
        size_t M = bra.first;
        size_t N = bra.second;
        if(H_shell_maxp[P] * H_shell_maxp[P] * J_metric_shell_diag[P] * eri_computers[rank]->shell_pair_value(M,N) < thresh2) {
            continue;
        }
        computed_triplets2++;
        int np = auxiliary_->shell(P).nfunction();
        int pstart = auxiliary_->shell(P).function_index();
        int nm = primary_->shell(M).nfunction();
        int mstart = primary_->shell(M).function_index();
        int nn = primary_->shell(N).nfunction();
        int nstart = primary_->shell(N).function_index();

        eri_computers[rank]->compute_shell(P, 0, M, N);
        const auto & buffer = eri_computers[rank]->buffers()[0];

        for(size_t jki = 0; jki < njk; jki++) {

            auto JTp = JT[jki][rank]->pointer();
            auto Hp = H[jki]->pointer();

            for (int p = pstart, index = 0; p < pstart + np; p++) {
                for (int m = mstart; m < mstart + nm; m++) {
                    for (int n = nstart; n < nstart + nn; n++, index++) {
                        JTp[m][n] += buffer[index] * Hp[p];
                        if (N != M) JTp[n][m] += buffer[index] * Hp[p];
                    }
                }
            }

        }

    }

    timer_off("ERI2");

    num_computed_shells_ = computed_triplets1 + computed_triplets2;
 
    for(size_t jki = 0; jki < njk; jki++) {
        for (size_t thread = 0; thread < nthreads_; thread++) {
            J[jki]->add(JT[jki][thread]);
        }
        J[jki]->hermitivitize();
    }

}

}  // namespace psi
