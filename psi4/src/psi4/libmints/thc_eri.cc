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
                                    Options& options) : THC_Computer(molecule, primary, options), auxiliary_(nullptr) {}

LS_THC_Computer::LS_THC_Computer(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, 
                                    Options& options) : THC_Computer(molecule, primary, options), auxiliary_(auxiliary) {}

LS_THC_Computer::~LS_THC_Computer() {}

void LS_THC_Computer::print_header() {
    outfile->Printf("   --------------------------------------------\n");
    outfile->Printf("       Least-Squares Tensor Hypercontraction   \n");
    outfile->Printf("                   by Andy Jiang               \n");
    outfile->Printf("              doi: 10.1063/1.4768233           \n");
    outfile->Printf("   --------------------------------------------\n\n");
    outfile->Printf("  Detailed LS-THC thresholds and cutoffs:\n");
    outfile->Printf("    LS_THC_SPHERICAL_POINTS   =   %3d \n", options_.get_int("LS_THC_SPHERICAL_POINTS"));
    outfile->Printf("    LS_THC_RADIAL_POINTS      =   %3d \n", options_.get_int("LS_THC_RADIAL_POINTS"));
    outfile->Printf("    LS_THC_BASIS_TOLERANCE    = %6.3e \n", options_.get_double("LS_THC_BASIS_TOLERANCE"));
    outfile->Printf("    LS_THC_WEIGHTS_TOLERANCE  = %6.3e \n", options_.get_double("LS_THC_WEIGHTS_TOLERANCE"));
    outfile->Printf("    LS_THC_S_EPSILON          = %6.3e \n", options_.get_double("LS_THC_S_EPSILON"));
    outfile->Printf("    Using DF?                     %3s \n", use_df_ ? "Yes" : "No");
    outfile->Printf("\n\n");
}

SharedMatrix LS_THC_Computer::build_E_exact() {
    /**
     * Parrish et al. 2012 Procedure 2
    */

    size_t nbf = primary_->nbf();
    size_t rank = x1_->nrow();

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    IntegralFactory factory(primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri_computers(nthreads);

    eri_computers[0] = std::shared_ptr<TwoBodyAOInt>(factory.eri());
#pragma omp parallel for
    for (int thread = 1; thread < nthreads; ++thread) {
        eri_computers[thread] = std::shared_ptr<TwoBodyAOInt>(eri_computers[0]->clone());
    }

    size_t nshellpair = eri_computers[0]->shell_pairs().size();

    SharedMatrix E_PQ = std::make_shared<Matrix>(rank, rank);

    // E^{PQ} = x_{m}^{P}x_{n}^{P}(mn|rs)x_{r}^{Q}x_{s}^{Q} (Parrish 2012 eq. 33)
#pragma omp parallel for
    for (size_t MN = 0; MN < nshellpair; ++MN) {

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        auto bra = eri_computers[thread]->shell_pairs()[MN];
        size_t M = bra.first;
        size_t N = bra.second;

        int nm = primary_->shell(M).nfunction();
        int mstart = primary_->shell(M).function_index();
        int nn = primary_->shell(N).nfunction();
        int nstart = primary_->shell(N).function_index();

        SharedMatrix E_temp = std::make_shared<Matrix>(rank, nm * nn);

        double prefactor1 = (M == N) ? 1.0 : 2.0;

        for (size_t RS = 0; RS < nshellpair; ++RS) {

            auto ket = eri_computers[thread]->shell_pairs()[RS];
            size_t R = ket.first;
            size_t S = ket.second;

            int nr = primary_->shell(R).nfunction();
            int rstart = primary_->shell(R).function_index();
            int ns = primary_->shell(S).nfunction();
            int sstart = primary_->shell(S).function_index();

            // TODO: Implement screening to make this more efficient
            eri_computers[thread]->compute_shell(M, N, R, S);
            const auto &buffer = eri_computers[thread]->buffers()[0];

            double prefactor2 = (R == S) ? 1.0 : 2.0;

            // This loop is forming the intermediate (E')^{Q}_{mn} = (mn|rs)x_{r}^{Q}x_{s}^{Q}
            for (size_t p = 0; p < rank; ++p) {
                for (size_t dm = 0, index = 0; dm < nm; ++dm) {
                    for (size_t dn = 0; dn < nn; ++dn) {
                        for (size_t r = rstart; r < rstart + nr; ++r) {
                            for (size_t s = sstart; s < sstart + ns; ++s, ++index) {
                                (*E_temp)(p, dm * nn + dn) += prefactor2 * buffer[index] * (*x1_)(p, r) * (*x1_)(p, s);
                            } // end s
                        } // end r
                    } // end dn
                } // end dm
            } // end p
        } // end RS

        // Final contractions E^{PQ} = x_{m}^{P}x_{n}^{P}(E')^{Q}_{mn}
        for (size_t p = 0; p < rank; ++p) {
            for (size_t q = 0; q < rank; ++q) {
                double e_pq_cont = 0.0;
                for (size_t m = mstart; m < mstart + nm; ++m) {
                    size_t dm = m - mstart;
                    for (size_t n = nstart; n < nstart + nn; ++n) {
                        size_t dn = n - nstart;
                        e_pq_cont += prefactor1 * (*E_temp)(q, dm * nn + dn) * (*x1_)(p, m) * (*x1_)(p, n);
                    } // end n
                } // end m
#pragma omp atomic
                (*E_PQ)(p, q) += e_pq_cont;
            } // end q
        } // end p
    } // end MN

    return E_PQ;
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

    // E^{IJ} = x_{m}^{I}x_{n}^{I}(mn|P)(P|Q)^{-1}(Q|rs)x_{r}^{J}x_{s}^{J} (Parrish 2012 eq. 34-35)
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

        // This loop is forming the intermediate (E')^{IP} = x_{m}^{I}x_{n}^{I}(P|mn)
        for (size_t r = 0; r < rank; ++r) {
            for (size_t p = pstart, index = 0; p < pstart + np; ++p) {
                double e_temp_cont = 0.0;
                for (size_t m = mstart; m < mstart + nm; ++m) {
                    for (size_t n = nstart; n < nstart + nn; ++n, ++index) {
                        e_temp_cont += prefactor * buffer[index] * (*x1_)(r, m) * (*x1_)(r, n);
                    } // end n
                } // end m
#pragma omp atomic
                (*E_temp)(r, p) += e_temp_cont;
            } // end p
        } // end r
    } // end MNP

    // Final contractions: (E)^{IJ} = (E')^{IP}(P|Q)^{-1}(E')^{JQ}
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

    use_df_ = options_.get_bool("LS_THC_DF");
    
    if (!options_["LS_THC_DF"].has_changed() && !auxiliary_) {
        outfile->Printf("    Warning: No auxiliary basis function specified for least-squares tensor hypercontraction... decomposing exact integrals!\n\n");
        use_df_ = false;
    } else if (options_.get_bool("LS_THC_DF") && !auxiliary_) {
        throw PSIEXCEPTION("    Attempting to do LS-THC with DF integrals w/o specifed auxiliary basis set :(\n\n");
    }

    print_header();

    size_t nbf = primary_->nbf();

    timer_on("LS-THC: Build Grid");

    // Create a grid for LS_THC
    std::map<std::string, std::string> grid_init_str_options = {
        {"DFT_PRUNING_SCHEME", options_.get_str("LS_THC_PRUNING_SCHEME")},
        {"DFT_RADIAL_SCHEME",  "TREUTLER"},
        {"DFT_NUCLEAR_SCHEME", "TREUTLER"},
        {"DFT_GRID_NAME",      ""},
        {"DFT_BLOCK_SCHEME",   "OCTREE"},
    };
    std::map<std::string, int> grid_init_int_options = {
        {"DFT_SPHERICAL_POINTS", options_.get_int("LS_THC_SPHERICAL_POINTS")}, 
        {"DFT_RADIAL_POINTS",    options_.get_int("LS_THC_RADIAL_POINTS")},
        {"DFT_BLOCK_MIN_POINTS", 100},
        {"DFT_BLOCK_MAX_POINTS", 256},
    };
    std::map<std::string, double> grid_init_float_options = {
        {"DFT_BASIS_TOLERANCE",   options_.get_double("LS_THC_BASIS_TOLERANCE")}, 
        {"DFT_BS_RADIUS_ALPHA",   1.0},
        {"DFT_PRUNING_ALPHA",     1.0},
        {"DFT_BLOCK_MAX_RADIUS",  3.0},
        {"DFT_WEIGHTS_TOLERANCE", options_.get_double("LS_THC_WEIGHTS_TOLERANCE")},
    };
    auto grid = DFTGrid(molecule_, primary_, grid_init_int_options, grid_init_str_options, grid_init_float_options, options_);
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
            x1_->scale_row(0, point_idx + p, std::pow(std::abs(w[p]), 0.25));
        }

        point_idx += block->npoints();
    }

    x2_ = x1_;
    x3_ = x1_;
    x4_ = x1_;

    timer_off("LS-THC: Build Grid");
    
    timer_on("LS-THC: Form E");

    SharedMatrix E_PQ = use_df_ ? build_E_df() : build_E_exact();

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
    SharedMatrix S_Qq_inv = S_Qq->pseudoinverse(options_.get_double("LS_THC_S_EPSILON"), nremoved);

    timer_off("LS-THC: Form S");

    timer_on("LS-THC: Form Z");

    Z_PQ_ = linalg::triplet(S_Qq_inv, E_PQ, S_Qq_inv);

    timer_off("LS-THC: Form Z");

    outfile->Printf("    Tensor Hypercontraction Complete! \n\n");
    outfile->Printf("    Number of Grid Points (Rank) : %6d (%3d per atom)\n\n", npoints, npoints / molecule_->natom());
    outfile->Printf("    Memory Required to Store Exact Integrals : %6.3f [GiB]\n", (nbf * (nbf + 1) / 2) * (nbf * (nbf + 1) / 2 + 1) / 2 * pow(2.0, -30) * sizeof(double));
    if (auxiliary_) outfile->Printf("    Memory Required to Store DF Integrals : %6.3f [GiB]\n", auxiliary_->nbf() * nbf * (nbf + 1) / 2 * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    Memory Required in LS-THC factored form  : %6.3f [GiB]\n\n", (npoints * (npoints + nbf)) * pow(2.0, -30) * sizeof(double));
}

}