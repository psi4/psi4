
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

#include "jk_grad.h"

#include "psi4/lib3index/3index.h"
#include "psi4/libqt/qt.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/electrostatic.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"

#include <vector>
#include <map>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif


using namespace psi;

namespace psi {
namespace scfgrad {

Matrix compute_numeric_overlap(const DFTGrid &grid, const std::shared_ptr<BasisSet> &primary) {

    // DOI 10.1063/1.3646921, EQ. 9

    int nbf = primary->nbf();
    BasisFunctions bf_computer(primary, grid.max_points(), grid.max_functions());
    Matrix S_num("Numerical Overlap", nbf, nbf);
    auto S_nump = S_num.pointer();

    // This loop could be parallelized over blocks of grid points. However, the cost of the loop is
    // so small (< 10 seconds for a 200 heavy atom system), parallelism isn't necessary

    for (const auto &block : grid.blocks()) {

        // grid points in this block
        int npoints_block = block->npoints();
        int nbf_block = block->local_nbf();
        auto w = block->w();

        // compute basis functions at these grid points
        bf_computer.compute_functions(block);
        auto point_values = bf_computer.basis_values()["PHI"];

        // resize the buffer of basis function values
        Matrix X_block("phi_g,u", npoints_block, nbf_block);  // points x nbf_block
        auto X_blockp = X_block.pointer();
        for (size_t p = 0; p < npoints_block; p++) {
            for (size_t k = 0; k < nbf_block; k++) {
                X_blockp[p][k] = point_values->get(p, k) * std::sqrt(w[p]);
            }
        }

        // significant basis functions at these grid points
        const auto &bf_map = block->functions_local_to_global();

        auto S_num_block = linalg::doublet(X_block, X_block, true, false);
        auto S_num_blockp = S_num_block.pointer();

        for (size_t mu_local = 0; mu_local < nbf_block; mu_local++) {
            size_t mu = bf_map[mu_local];
            for (size_t nu_local = 0; nu_local < nbf_block; nu_local++) {
                size_t nu = bf_map[nu_local];
                S_nump[mu][nu] += S_num_blockp[mu_local][nu_local];
            }
        }

    }

    S_num.hermitivitize();

    return S_num;

}

Matrix compute_esp_bound(const BasisSet &primary) {

    // DOI 10.1016/j.chemphys.2008.10.036, EQ. 20
    // This is a pretty loose ESP bound, which should eventually be swapped out for something tighter
    // The bound is also only based on the overlap between the basis functions, not the distance 
    // between the basis functions and the grid point.

    int nshell = primary.nshell();

    Matrix esp_bound("Shell Integral Bound", nshell, nshell);
    auto esp_boundp = esp_bound.pointer();

    auto dist = primary.molecule()->distance_matrix();
    auto distp = dist.pointer();

    for (size_t s1=0; s1 < nshell; s1++) {
        int c1 = primary.shell_to_center(s1);
        int np1 = primary.shell(s1).nprimitive();
        for (size_t s2=0; s2 < nshell; s2++) {
            int c2 = primary.shell_to_center(s2);
            int np2 = primary.shell(s2).nprimitive();

            double r2 = distp[c1][c2] * distp[c1][c2] ;
            for(size_t pi1 = 0; pi1 < np1; pi1++) {
                for(size_t pi2 = 0; pi2 < np2; pi2++) {
                    double exp1 = primary.shell(s1).exp(pi1);
                    double exp2 = primary.shell(s2).exp(pi2);
                    double coef1 = primary.shell(s1).coef(pi1);
                    double coef2 = primary.shell(s2).coef(pi2);
                    esp_boundp[s1][s2] += coef1 * coef2 * std::exp(-1 * r2 * exp1 * exp2 / (exp1 + exp2)) * 2 * M_PI / (exp1 + exp2);
                }
            }

            esp_boundp[s1][s2] = std::abs(esp_boundp[s1][s2]);
        }
    }

    return esp_bound;

}

DFJCOSKGrad::DFJCOSKGrad(int deriv, std::shared_ptr<MintsHelper> mints, Options& options) :
    JKGrad(deriv,mints->get_basisset("ORBITAL")), auxiliary_(mints->get_basisset("DF_BASIS_SCF")), mints_(mints), options_(options)
{
    common_init(); 
}
DFJCOSKGrad::~DFJCOSKGrad() {}
void DFJCOSKGrad::common_init() {
    ints_num_threads_ = 1;
#ifdef _OPENMP
    ints_num_threads_ = Process::environment.get_n_threads();
#endif

    // => Direct Density-Fitted Coulomb Setup <= //

    // pre-compute coulomb fitting metric
    /*
    timer_on("Coulomb Metric");
    FittingMetric J_metric_obj(auxiliary_, true);
    J_metric_obj.form_fitting_metric();
    J_metric_ = J_metric_obj.get_metric();
    timer_off("Coulomb Metric");

    // pre-construct per-thread TwoBodyAOInt objects for computing 3-index ERIs
    timer_on("ERI Computers");
    eri_computers_.resize(nthreads_);
    auto zero = BasisSet::zero_ao_basis_set();
    IntegralFactory rifactory(auxiliary_, zero, primary_, primary_);
    eri_computers_[0] = std::shared_ptr<TwoBodyAOInt>(rifactory.eri());
    for(int rank = 1; rank < nthreads_; rank++) {
        eri_computers_[rank] = std::shared_ptr<TwoBodyAOInt>(eri_computers_.front()->clone());
    }
    timer_off("ERI Computers");
    */

    // => Chain of Spheres Exchange Setup <= //

    timer_on("Grid Construction");

    // TODO: specify bool "DFT_REMOVE_DISTANT_POINTS" in the DFTGrid constructors

    // Create a large DFTGrid for the gradient
    std::map<std::string, std::string> grid_final_str_options = {
        {"DFT_PRUNING_SCHEME", options_.get_str("COSX_PRUNING_SCHEME")},
        {"DFT_RADIAL_SCHEME",  "TREUTLER"},
        {"DFT_NUCLEAR_SCHEME", "TREUTLER"},
        {"DFT_GRID_NAME",      ""},
        {"DFT_BLOCK_SCHEME",   "OCTREE"},
    };
    std::map<std::string, int> grid_final_int_options = {
        {"DFT_SPHERICAL_POINTS", options_.get_int("COSX_SPHERICAL_POINTS_FINAL")}, 
        {"DFT_RADIAL_POINTS",    options_.get_int("COSX_RADIAL_POINTS_FINAL")},
        {"DFT_BLOCK_MIN_POINTS", 100},
        {"DFT_BLOCK_MAX_POINTS", 256},
    };
    std::map<std::string, double> grid_final_float_options = {
        {"DFT_BASIS_TOLERANCE",   options_.get_double("COSX_BASIS_TOLERANCE")}, 
        {"DFT_BS_RADIUS_ALPHA",   1.0},
        {"DFT_PRUNING_ALPHA",     1.0},
        {"DFT_BLOCK_MAX_RADIUS",  3.0},
        {"DFT_WEIGHTS_TOLERANCE", 1e-15},
    };
    grid_final_ = std::make_shared<DFTGrid>(primary_->molecule(), primary_, grid_final_int_options, grid_final_str_options, grid_final_float_options, options_);

    timer_off("Grid Construction");

    // => Overlap Fitting Metric <= //

    // Fit an overlap metric (Q) for both grids to reduce numerical error

    // DOI 10.1063/1.3646921, EQ. 18
    // Note: the above reference defines Q as S_an @ S_num^{-1} @ X
    // Here, Q refers to just S_ @ S_num^{-1} (no X)
    // This Q is contracted with X later to agree with the literature definition

    timer_on("Numeric Overlap");

    // compute the numeric overlap matrix for each grid
    auto S_num = compute_numeric_overlap(*grid_final_, primary_ );

    timer_off("Numeric Overlap");

    timer_on("Analytic Overlap");

    // compute the analytic overlap matrix
    MintsHelper helper(primary_, options_);
    auto S_an = helper.ao_overlap();

    timer_off("Analytic Overlap");

    // form the overlap metric (Q) for each grid

    timer_on("Overlap Metric Solve");

    int nbf = primary_->nbf();
    std::vector<int> ipiv(nbf);

    // solve: Q_final_ = S_an @ S_num_final_^{-1}
    Q_final_ = S_an->clone();
    C_DGESV(nbf, nbf, S_num.pointer()[0], nbf, ipiv.data(), Q_final_->pointer()[0], nbf);

    timer_off("Overlap Metric Solve");
    outfile->Printf("COSX gradients successfully initialized.\n");
}
void DFJCOSKGrad::print_header() const {
    std::string screen_type = options_.get_str("SCREENING");
    if (print_) {
        outfile->Printf("  ==> DFJCOSKGrad: Chain-of-Spheres SCF Gradients <==\n\n");

        outfile->Printf("    Gradient:          %11d\n", deriv_);
        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:             %11.3E\n", omega_);
        outfile->Printf("    Integrals threads: %11d\n", ints_num_threads_);
        outfile->Printf("    Memory [MiB]:       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf("    J Screening Type:   %11s\n", screen_type.c_str());
        outfile->Printf("    J Screening Cutoff: %11.0E\n", cutoff_);
        outfile->Printf("    K Screening Cutoff: %11.0E\n", options_.get_double("COSX_INTS_TOLERANCE"));
        outfile->Printf("    K Density Cutoff:   %11.0E\n", options_.get_double("COSX_DENSITY_TOLERANCE"));
        outfile->Printf("    K Basis Cutoff:     %11.0E\n", options_.get_double("COSX_BASIS_TOLERANCE"));
        outfile->Printf("    K Overlap Fitting:  %11s\n", (options_.get_bool("COSX_OVERLAP_FITTING") ? "Yes" : "No"));
        outfile->Printf("\n");
    }
}
void DFJCOSKGrad::compute_gradient() {
    if (!do_J_ && !do_K_ && !do_wK_) return;

    if (!(Ca_ && Cb_ && Da_ && Db_ && Dt_))
        throw PSIEXCEPTION("Occupation/Density not set");
    
    // => Set up gradients <= //
    int natom = primary_->molecule()->natom();
    gradients_.clear();
    if (do_J_) {
        gradients_["Coulomb"] = std::make_shared<Matrix>("Coulomb Gradient", natom, 3);
    }
    if (do_K_) {
        gradients_["Exchange"] = std::make_shared<Matrix>("Exchange Gradient", natom, 3);
    }
    if (do_wK_) {
        gradients_["Exchange,LR"] = std::make_shared<Matrix>("Exchange,LR Gradient", natom, 3);
    }

    if (do_J_) {
        build_JGrad();
        outfile->Printf("coulomb gradients successfully calculated.\n");
    }
    if (do_K_) {
        build_KGrad();
        outfile->Printf("exchange gradients successfully calculated.\n");
    }
    if (do_wK_) {
        throw PSIEXCEPTION("No wK for DFJCOSK.");
    }

    // Printing
    gradients_["Coulomb"]->print();
    if (do_K_) {
        gradients_["Exchange"]->print();
     }
}



void DFJCOSKGrad::build_KGrad() {
    if (!do_K_) return;
    timer_on("K Grad");

    auto nthreads_ = ints_num_threads_;
    std::vector<SharedMatrix> D;
    D.push_back(Dt_);

    // => Sizing <= //
    int njk = D.size();
    int nbf = primary_->nbf();
    int nshell = primary_->nshell();
    int natom = primary_->molecule()->natom();

    // => Knobs <= //
    //double gscreen = options_.get_double("COSX_INTS_TOLERANCE");
    double kscreen = options_.get_double("COSX_INTS_TOLERANCE");
    double dscreen = options_.get_double("COSX_DENSITY_TOLERANCE");
    bool overlap_fitted = options_.get_bool("COSX_OVERLAP_FITTING");

    // use a small DFTGrid grid (and overlap metric) for early SCF iterations
    // otherwise use a large DFTGrid
    auto grid = grid_final_;
    auto Q = Q_final_;

    // => Initialization <= //

    // per-thread ElectrostaticInt object (for computing one-electron "pseudospectral" integrals)
    std::vector<std::shared_ptr<ElectrostaticInt>> int_computers(nthreads_);

    // per-thread BasisFunctions object (for computing basis function values at grid points)
    std::vector<std::shared_ptr<BasisFunctions>> bf_computers(nthreads_);

    // per-thread K Matrix buffers (for accumulating thread contributions to K)
    std::vector<std::vector<SharedMatrix>> KT(njk, std::vector<SharedMatrix>(nthreads_));
    std::vector<std::vector<SharedMatrix>> KTgrad(njk, std::vector<SharedMatrix>(nthreads_));

    // initialize per-thread objects
    IntegralFactory factory(primary_);
    for(size_t thread = 0; thread < nthreads_; thread++) {
        int_computers[thread] = std::shared_ptr<ElectrostaticInt>(static_cast<ElectrostaticInt *>(factory.electrostatic().release()));
        bf_computers[thread] = std::make_shared<BasisFunctions>(primary_, grid->max_points(), grid->max_functions());
        bf_computers[thread]->set_deriv(1);
        for(size_t jki = 0; jki < njk; jki++) {
            KT[jki][thread] = std::make_shared<Matrix>(nbf, nbf);
            KTgrad[jki][thread] = std::make_shared<Matrix>(natom, 3);
        }
    }

    // precompute bounds for the one-electron integrals
    auto esp_bound = compute_esp_bound(*primary_);
    auto esp_boundp = esp_bound.pointer();

    // inter-atom and inter-shell distances [Bohr]
    auto dist = primary_->molecule()->distance_matrix();
    auto shell_dist = std::make_shared<Matrix>(nshell, nshell);
    for(size_t s1 = 0; s1 < nshell; s1++) {
        size_t c1 = primary_->shell_to_center(s1);
        for(size_t s2 = 0; s2 < nshell; s2++) {
            size_t c2 = primary_->shell_to_center(s2);
            shell_dist->set(s1, s2, dist.get(c1, c2));
        }
    }

    // extent of each basis shell [Bohr]
    auto shell_extents = grid->extents()->shell_extents();

    // map of shell pairs with overlapping extents
    std::vector<std::vector<int>> shell_extent_map(nshell);
    for(size_t s1 = 0; s1 < nshell; s1++) {
        for(size_t s2 = 0; s2 < nshell; s2++) {
            if (shell_dist->get(s1, s2) <= shell_extents->get(s2) + shell_extents->get(s1)) {
                shell_extent_map[s1].push_back(s2);
            }
        }
    }

    // => Integral Computation <= //

    // benchmarking statistics
    size_t int_shells_total = 0;
    size_t int_shells_computed = 0;

    timer_on("Grid Loop");

    // The primary COSK loop over blocks of grid points
#pragma omp parallel for schedule(dynamic) num_threads(nthreads_) reduction(+ : int_shells_total, int_shells_computed)
    for (size_t bi = 0; bi < grid->blocks().size(); bi++) {

        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        // grid points in this block
        auto block = grid->blocks()[bi];
        int npoints_block = block->npoints();
        auto x = block->x();
        auto y = block->y();
        auto z = block->z();
        auto w = block->w();

        // significant basis functions and shells at these grid points
        // significance determined via basis extent
        const auto &bf_map = block->functions_local_to_global();
        const auto &shell_map = block->shells_local_to_global();
        int nbf_block = bf_map.size();
        int ns_block = shell_map.size();

        // lists of all basis functions and shells
        //
        // The use of these "all" lists adds O(N^2) cost to the COSK grid loop (w/ small prefactor)
        // This cost is negligible relative to the esp integral computation, which is O(N) (w/ a much larger prefactor),
        // but future COSK work could remove this potential bottleneck
        std::vector<int> bf_map_all;
        std::vector<int> shell_map_all;
        for (size_t bf = 0; bf < nbf; bf++) bf_map_all.push_back(bf);
        for (size_t s = 0; s < nshell; s++) shell_map_all.push_back(s);
        int nbf_block_all = bf_map_all.size();
        int ns_block_all = shell_map_all.size();

        // => Bookkeeping <= //

        // map index in shell_map_all to first index in bf_map_all
        std::vector<int> shell_map_all_to_bf_map_all;
        if (shell_map_all.size() > 0) {
            shell_map_all_to_bf_map_all.push_back(0);
            for(size_t shell_map_ind = 0; (shell_map_ind + 1) < shell_map_all.size(); shell_map_ind++) {
                size_t MU = shell_map_all[shell_map_ind];
                shell_map_all_to_bf_map_all.push_back(primary_->shell(MU).nfunction() + shell_map_all_to_bf_map_all.back());
            }
        }
        
        // map index in shell_map to first index in bf_map
        std::vector<int> shell_map_to_bf_map;
        if (shell_map.size() > 0) {
            shell_map_to_bf_map.push_back(0);
            for(size_t shell_map_ind = 0; (shell_map_ind + 1) < shell_map.size(); shell_map_ind++) {
                size_t MU = shell_map[shell_map_ind];
                shell_map_to_bf_map.push_back(primary_->shell(MU).nfunction() + shell_map_to_bf_map.back());
            }
        }

        // map back from global shell index to local shell index
        std::map<size_t, size_t> shell_map_inv;
        for (size_t shell_ind = 0; shell_ind < shell_map.size(); shell_ind++) {
            shell_map_inv[shell_map[shell_ind]] = shell_ind;
        }

        // => Process Density Matrix <= //

        // significant cols of D for this grid block
        std::vector<SharedMatrix> D_block(njk);
        for(size_t jki = 0; jki < njk; jki++) {
            D_block[jki] = std::make_shared<Matrix>(nbf_block_all, nbf_block);
        }

        for(size_t jki = 0; jki < njk; jki++) {
            auto Dp = D[jki]->pointer();
            auto D_blockp = D_block[jki]->pointer();
            for (size_t tau_ind = 0; tau_ind < nbf_block_all; tau_ind++) {
                size_t tau = bf_map_all[tau_ind];
                for (size_t kappa_ind = 0; kappa_ind < nbf_block; kappa_ind++) {
                    size_t kappa = bf_map[kappa_ind];
                    D_blockp[tau_ind][kappa_ind] = Dp[tau][kappa];
                }
            }
        }

        // shell-pair maxima of D_block
        auto D_block_shell = std::make_shared<Matrix>(ns_block_all, ns_block);
        auto D_block_shellp = D_block_shell->pointer();

        for (size_t TAU_ind = 0; TAU_ind < ns_block_all; TAU_ind++) {
            size_t TAU = shell_map_all[TAU_ind];
            size_t tau_start = shell_map_all_to_bf_map_all[TAU_ind];
            size_t num_tau = primary_->shell(TAU).nfunction();
            for (size_t KAPPA_ind = 0; KAPPA_ind < ns_block; KAPPA_ind++) {
                size_t KAPPA = shell_map[KAPPA_ind];
                size_t kappa_start = shell_map_to_bf_map[KAPPA_ind];
                size_t num_kappa = primary_->shell(KAPPA).nfunction();
                for(size_t jki = 0; jki < njk; jki++) {
                    auto D_blockp = D_block[jki]->pointer();
                    for (size_t bf1 = tau_start; bf1 < tau_start + num_tau; bf1++) {
                        for (size_t bf2 = kappa_start; bf2 < kappa_start + num_kappa; bf2++) {
                            D_block_shellp[TAU_ind][KAPPA_ind] = std::max(D_block_shellp[TAU_ind][KAPPA_ind], std::abs(D_blockp[bf1][bf2]));
                        }
                    }
                }
            }
        }

        // significant TAU shells determined from sparsity of the density matrix 
        // i.e. KAPPA -> TAU sparsity. Refered to by Neese as a "p-junction"
        std::vector<int> shell_map_tau;

        for(size_t TAU = 0; TAU < ns_block_all; TAU++) {
            for(size_t KAPPA_ind = 0; KAPPA_ind < ns_block; KAPPA_ind++) {
                size_t KAPPA = shell_map[KAPPA_ind];
                if (D_block_shellp[TAU][KAPPA_ind] > dscreen) {
                    shell_map_tau.push_back(TAU);
                    break;
                }
            }
        }

        // => X Matrix <= //

        // DOI 10.1016/j.chemphys.2008.10.036, EQ. 4

        // compute basis functions at these grid points
        if(rank == 0) timer_on("Basis Derivs");
        bf_computers[rank]->compute_functions(block);
        auto point_values = bf_computers[rank]->basis_values()["PHI"];
        auto point_values_dx = bf_computers[rank]->basis_values()["PHI_X"];
        auto point_values_dy = bf_computers[rank]->basis_values()["PHI_Y"];
        auto point_values_dz = bf_computers[rank]->basis_values()["PHI_Z"];
        if(rank ==0) timer_off("Basis Derivs");

        // resize the buffer of basis function values
        auto X_block = std::make_shared<Matrix>(npoints_block, nbf_block);  // points x nbf_block
        auto X_blockp = X_block->pointer();
        for (size_t p = 0; p < npoints_block; p++) {
            for (size_t k = 0; k < nbf_block; k++) {
                X_blockp[p][k] = point_values->get(p, k) * std::sqrt(w[p]);
            }
        }

        // absmax of X matrix over basis functions (row maximum) needed for screening
        Vector X_block_bfmax(npoints_block);
        auto X_block_bfmaxp = X_block_bfmax.pointer();
        for (size_t p = 0; p < npoints_block; p++) {
            for (size_t k = 0; k < nbf_block; k++) {
                X_block_bfmaxp[p] = std::max(X_block_bfmaxp[p], std::abs(X_blockp[p][k]));
            }
        }

        double X_block_max = X_block->absmax();

        // => F Matrix <= //

        // DOI 10.1016/j.chemphys.2008.10.036, EQ. 6

        // contract density with basis functions values at these grid points
        std::vector<SharedMatrix> F_block(njk);
        for(size_t jki = 0; jki < njk; jki++) {
            F_block[jki] = linalg::doublet(X_block, D_block[jki], false, true);
        }
        
        // shell maxima of F_block
        auto F_block_shell = std::make_shared<Matrix>(npoints_block, nshell);
        auto F_block_shellp = F_block_shell->pointer();

        // grid point maxima of F_block_gmax
        auto F_block_gmax = std::make_shared<Vector>(nshell);
        auto F_block_gmaxp = F_block_gmax->pointer();

        for (size_t p = 0; p < npoints_block; p++) {
            for (size_t TAU_local = 0; TAU_local < shell_map_all.size(); TAU_local++) {
                size_t TAU = shell_map_all[TAU_local];
                size_t num_tau = primary_->shell(TAU).nfunction();
                size_t tau_start = shell_map_all_to_bf_map_all[TAU_local];
                for(size_t jki = 0; jki < njk; jki++) {
                    auto F_blockp = F_block[jki]->pointer();
                    for (size_t tau = tau_start; tau < tau_start + num_tau; tau++) {
                        F_block_shellp[p][TAU_local] = std::max(F_block_shellp[p][TAU_local], std::abs(F_blockp[p][tau]));
                        F_block_gmaxp[TAU_local] = std::max(F_block_gmaxp[TAU_local], std::abs(F_blockp[p][tau]));
                    }
                }
            }
        }

        // => Q Matrix <= //

        // DOI 10.1063/1.3646921, EQ. 18

        // slice of overlap metric (Q) made up of significant basis functions at this grid point
        auto Q_block = std::make_shared<Matrix>(nbf_block, nbf_block);
        for(size_t mu_local = 0; mu_local < nbf_block; mu_local++) {
            size_t mu = bf_map[mu_local];
            for(size_t nu_local = 0; nu_local < nbf_block; nu_local++) {
                size_t nu = bf_map[nu_local];
                Q_block->set(mu_local, nu_local, Q->get(mu, nu));
            }
        }

        // now Q_block agrees with EQ. 18 (see note about Q_init_ and Q_final_ in common_init())
        Q_block = linalg::doublet(X_block, Q_block, false, true);

        // => G Matrix <= //

        // DOI 10.1016/j.chemphys.2008.10.036, EQ. 7

        std::vector<SharedMatrix> G_block(njk);
        for(size_t jki = 0; jki < njk; jki++) {
            G_block[jki] = std::make_shared<Matrix>(nbf_block_all, npoints_block);
        }

        if(rank == 0) timer_on("ESP Integrals");

        const auto & int_buff = int_computers[rank]->buffers()[0];

        // calculate A_NU_TAU at all grid points in this block
        // contract A_NU_TAU with F_TAU to get G_NU
        for (size_t TAU : shell_map_tau) {
            const size_t num_tau = primary_->shell(TAU).nfunction();
            const size_t tau_start = shell_map_all_to_bf_map_all[TAU];
            const size_t center_TAU = primary_->shell_to_center(TAU);
            const double x_TAU = primary_->molecule()->x(center_TAU);
            const double y_TAU = primary_->molecule()->y(center_TAU);
            const double z_TAU = primary_->molecule()->z(center_TAU);

            // TAU -> NU sparity determined by shell extents
            for (size_t NU : shell_extent_map[TAU]) {
                const size_t num_nu = primary_->shell(NU).nfunction();
                const size_t nu_start = shell_map_all_to_bf_map_all[NU];
                const size_t center_NU = primary_->shell_to_center(NU);
                const double x_NU = primary_->molecule()->x(center_NU);
                const double y_NU = primary_->molecule()->y(center_NU);
                const double z_NU = primary_->molecule()->z(center_NU);

                // is this value of NU also a possible value of TAU for this grid block?
                // i.e. can we use permutational symmetry of this (NU|TAU) integral shell pair?
                bool symm = (NU != TAU) && std::binary_search(shell_map_tau.begin(), shell_map_tau.end(), NU);

                // we've already done these integrals
                if (symm && TAU > NU) continue;

                // benchmarking
                int_shells_total += npoints_block;

                // can we screen the whole block over K_uv = (X_ug (A_vtg (F_tg)) upper bound?
                double k_bound = X_block_max * esp_boundp[NU][TAU] * F_block_gmaxp[TAU];
                if (symm) k_bound = std::max(k_bound, X_block_max * esp_boundp[TAU][NU] * F_block_gmaxp[NU]);
                if (k_bound < kscreen) continue;

                for (size_t g = 0; g < npoints_block; g++) {

                    // grid-point specific screening
                    // account for the distance between the grid point and the shell pair
                    double dist_TAU_g = std::sqrt((x_TAU - x[g]) * (x_TAU - x[g]) + (y_TAU - y[g]) * (y_TAU - y[g]) + (z_TAU - z[g]) * (z_TAU - z[g]));
                    double dist_NU_g = std::sqrt((x_NU - x[g]) * (x_NU - x[g]) + (y_NU - y[g]) * (y_NU - y[g]) + (z_NU - z[g]) * (z_NU - z[g]));
                    double dist_NUTAU_g = std::min(dist_TAU_g - shell_extents->get(TAU), dist_NU_g - shell_extents->get(NU));
                    double dist_decay = 1.0 / std::max(1.0, dist_NUTAU_g);

                    // can we screen this single point over K_uv = (X_ug (A_vtg (F_tg))) upper bound?
                    k_bound = X_block_bfmaxp[g] * esp_boundp[NU][TAU] * dist_decay * F_block_shellp[g][TAU];
                    if (symm) k_bound = std::max(k_bound, X_block_bfmaxp[g] * esp_boundp[TAU][NU] * dist_decay * F_block_shellp[g][NU]);
                    if (k_bound < kscreen) continue;

                    // calculate pseudospectral integral shell pair (A_NU_TAU) at gridpoint g
                    int_computers[rank]->set_origin({x[g], y[g], z[g]});
                    int_computers[rank]->compute_shell(NU, TAU);

                    // benchmarking
                    int_shells_computed++;

                    // contract A_nu_tau with F_tau to get contribution to G_nu 
                    // symmetry permitting, also contract A_nu_tau with F_nu to get contribution to G_tau
                    for(size_t jki = 0; jki < njk; jki++) {
                        auto F_blockp = F_block[jki]->pointer();
                        auto G_blockp = G_block[jki]->pointer();
                        for (size_t nu = nu_start, index = 0; nu < (nu_start + num_nu); ++nu) {
                            for (size_t tau = tau_start; tau < (tau_start + num_tau); ++tau, index++) {
                                G_blockp[nu][g] += int_buff[index] * F_blockp[g][tau];
                                if (symm) G_blockp[tau][g] += int_buff[index] * F_blockp[g][nu];
                            }
                        }
                    }

                }
            }

        }

        if(rank == 0) timer_off("ESP Integrals");

        // Contract X (or Q if overlap fitting) with G to get contribution to K
        for(size_t jki = 0; jki < njk; jki++) {
            SharedMatrix KT_block;
            if (overlap_fitted) {
                KT_block = linalg::doublet(Q_block, G_block[jki], true, true);
            } else {
                KT_block = linalg::doublet(X_block, G_block[jki], true, true);
            }
            auto KT_blockp = KT_block->pointer();
            auto KTp = KT[jki][rank]->pointer();
            for(size_t mu_ind = 0; mu_ind < bf_map.size(); mu_ind++) {
                size_t mu = bf_map[mu_ind];
                for(size_t nu_ind = 0; nu_ind < bf_map_all.size(); nu_ind++) {
                    size_t nu = bf_map_all[nu_ind];
                    KTp[mu][nu] += KT_blockp[mu_ind][nu_ind];
                }
            }
        }

        if (rank == 0) timer_on("Contract K Gradient");

        for (size_t jki = 0; jki < njk; jki++) {

        // which atom centers have basis functions in this block?
        std::vector<int> centers_in_block;

        for (size_t k = 0; k < nbf_block; k++) {
            auto center = primary_->function_to_center(bf_map[k]);
            if (std::find(centers_in_block.begin(), centers_in_block.end(), center) == centers_in_block.end()){
                centers_in_block.push_back(center);
            }
        }

        //for (auto center : centers_in_block) {
            //  mask X_dx matrices:
        auto Xdx_block = std::make_shared<Matrix>(npoints_block, nbf_block);  // points x nbf_block
        auto Xdy_block = std::make_shared<Matrix>(npoints_block, nbf_block);  // points x nbf_block
        auto Xdz_block = std::make_shared<Matrix>(npoints_block, nbf_block);  // points x nbf_block
        auto Xdx_blockp = Xdx_block->pointer();
        auto Xdy_blockp = Xdy_block->pointer();
        auto Xdz_blockp = Xdz_block->pointer();
        if (rank == 0) timer_on("Form Xdx");
        for (size_t p = 0; p < npoints_block; p++) {
            for (size_t k = 0; k < nbf_block; k++) {
                Xdx_blockp[p][k] = point_values_dx->get(p, k) * std::sqrt(w[p]);
                Xdy_blockp[p][k] = point_values_dy->get(p, k) * std::sqrt(w[p]);
                Xdz_blockp[p][k] = point_values_dz->get(p, k) * std::sqrt(w[p]);
            }
        }
        if (rank == 0) timer_off("Form Xdx");
        if (rank == 0) timer_on("Form KT_dx");
        //  contract with G_block:
        SharedMatrix KT_dx_block;
        SharedMatrix KT_dy_block;
        SharedMatrix KT_dz_block;
        if (overlap_fitted) {
            //Qdx_block = linalg::doublet(Xdx_block, Q_block, false, true);
            //Qdy_block = linalg::doublet(Xdy_block, Q_block, false, true);
            //Qdz_block = linalg::doublet(Xdz_block, Q_block, false, true);
            KT_dx_block = linalg::doublet(Xdx_block, G_block[jki], true, true);
            KT_dy_block = linalg::doublet(Xdy_block, G_block[jki], true, true);
            KT_dz_block = linalg::doublet(Xdz_block, G_block[jki], true, true);
        } else {
            KT_dx_block = linalg::doublet(Xdx_block, G_block[jki], true, true);
            KT_dy_block = linalg::doublet(Xdy_block, G_block[jki], true, true);
            KT_dz_block = linalg::doublet(Xdz_block, G_block[jki], true, true);
        }
        if (rank == 0) timer_off("Form KT_dx");
        if (rank == 0) timer_on("Contract Density");
        //  contract with density, load into derivative matrix
        auto D_blockp = D_block[jki]->pointer();
        auto KT_dx_blockp = KT_dx_block->pointer();
        auto KT_dy_blockp = KT_dy_block->pointer();
        auto KT_dz_blockp = KT_dz_block->pointer();
        auto KTgradp = KTgrad[jki][rank]->pointer();
        for (size_t mu = 0; mu < nbf_block_all; mu++) {
            for (size_t nu = 0; nu < nbf_block; nu++) {
                auto center = primary_->function_to_center(bf_map[nu]);
                KTgradp[center][0] -= D_blockp[mu][nu] * KT_dx_blockp[nu][mu];
                KTgradp[center][1] -= D_blockp[mu][nu] * KT_dy_blockp[nu][mu];
                KTgradp[center][2] -= D_blockp[mu][nu] * KT_dz_blockp[nu][mu];
            }
        }
        if (rank == 0) timer_off("Contract Density");
    }

    if (rank == 0) timer_off("Contract K Gradient");

    }

    timer_off("Grid Loop");

    // Reduce per-thread contributions
    for(size_t jki = 0; jki < njk; jki++) {
        for (size_t thread = 0; thread < nthreads_; thread++) {
            //K[jki]->add(KT[jki][thread]);
            gradients_["Exchange"]->add(KTgrad[jki][thread]);
        }
        //if (lr_symmetric_) {
        //    K[jki]->hermitivitize();
        //}
    }

    if (bench_) {
        auto mode = std::ostream::app;
        PsiOutStream printer("bench.dat", mode);
        size_t ints_per_atom = int_shells_computed  / (size_t) natom;
        printer.Printf("COSK ESP Shells: %zu,%zu,%zu\n", ints_per_atom, int_shells_computed, int_shells_total);
    }

    timer_off("K Grad");
}

void DFJCOSKGrad::build_JGrad() {
    if (!do_J_) return;
    timer_on("J Grad");
    int nthread_df = ints_num_threads_;

    timer_on("contract c_A");
    // contract c_A
    // c_A = (B|pq) Dt_pq

    // => Sizing <= //

    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    // => Integrals <= //

    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(ints_num_threads_);
    eri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)); //eri(1) ??
    for (int t = 1; t < ints_num_threads_; t++) {
        eri[t] = std::shared_ptr<TwoBodyAOInt>(eri.front()->clone());
    }

    const std::vector<std::pair<int, int>>& shell_pairs = eri[0]->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //

    int max_rows;
    int maxP = auxiliary_->max_function_per_shell();
    size_t row_cost = 0L;
    row_cost += nso * (size_t)nso;
    size_t rows = memory_ / row_cost;
    rows = (rows > naux ? naux : rows);
    rows = (rows < maxP ? maxP : rows);
    max_rows = (int)rows;

    // => Block Sizing <= //

    std::vector<int> Pstarts;
    int counter = 0;
    Pstarts.push_back(0);
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        int nP = auxiliary_->shell(P).nfunction();
        if (counter + nP > max_rows) {
            counter = 0;
            Pstarts.push_back(P);
        }
        counter += nP;
    }
    Pstarts.push_back(auxiliary_->nshell());

    // => Temporary Buffers <= //

    SharedVector c;
    double* cp;

    c = std::make_shared<Vector>("c", naux);
    cp = c->pointer();

    SharedMatrix Amn;
    double** Amnp;
    Amn = std::make_shared<Matrix>("Amn", max_rows, nso * (size_t)nso);
    Amnp = Amn->pointer();

    double** Dtp = Dt_->pointer();

    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {
        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop = Pstarts[block + 1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop = (Pstop == auxiliary_->nshell() ? naux : auxiliary_->shell(Pstop).function_index());
        int np = pstop - pstart;

        // > Clear Integrals Register < //
        ::memset((void*)Amnp[0], '\0', sizeof(double) * np * nso * nso);

        // > Integrals < //
#pragma omp parallel for schedule(dynamic) num_threads(nthread_df)
        for (long int PMN = 0L; PMN < static_cast<long>(NP) * npairs; PMN++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int P = PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            eri[thread]->compute_shell(P, 0, M, N);

            const double* buffer = eri[thread]->buffer();

            int nP = auxiliary_->shell(P).nfunction();
            int oP = auxiliary_->shell(P).function_index() - pstart;

            int nM = primary_->shell(M).nfunction();
            int oM = primary_->shell(M).function_index();

            int nN = primary_->shell(N).nfunction();
            int oN = primary_->shell(N).function_index();

            for (int p = oP; p < nP + oP; p++) {
                for (int m = oM; m < nM + oM; m++) {
                    for (int n = oN; n < nN + oN; n++) {
                        Amnp[p][(m) * nso + (n)] = Amnp[p][(n) * nso + (m)] = *buffer++;
                    }
                }
            }
        }

        // > (A|mn) D_mn -> c_A < //
        if (do_J_) {
            C_DGEMV('N', np, nso * (size_t)nso, 1.0, Amnp[0], nso * (size_t)nso, Dtp[0], 1, 0.0, &cp[pstart], 1);
        }
    }

    timer_off("contract c_A");
    timer_on("contract d_A");
    // contract d_A
    // d_A = (A|B)^-1 c_A

    auto metric = std::make_shared<FittingMetric>(auxiliary_, true);
    metric->form_full_eig_inverse(condition_);
    SharedMatrix J = metric->get_metric();
    double** Jp = J->pointer();

    auto d = std::make_shared<Vector>("d", naux);
    double* dp = d->pointer();

    C_DGEMV('N', naux, naux, 1.0, Jp[0], naux, cp, 1, 0.0, dp, 1);


    timer_off("contract d_A");
    timer_on("metric grad");
    // metrix grad term
    // J = 0.5 * (A|B)^x d_B d_A

    std::map<std::string, SharedMatrix> densities;

    auto D = std::make_shared<Matrix>("D", naux, naux);
    auto Dp = D->pointer();
    C_DGER(naux, naux, 1, dp, 1, dp, 1, Dp[0], naux);
    densities["Coulomb"] = D;

    auto results = mints_->metric_grad(densities, "DF_BASIS_SCF");

    for (const auto& kv: results) {
        gradients_[kv.first] = kv.second;
    }
 

    timer_off("metric grad");
    timer_on("eri grad");
    // eri grad term
    // J += (A|pq)^x d_A Dt_pq

    int natom = primary_->molecule()->natom();
    std::vector<SharedMatrix> Jtemps;
    for (int t = 0; t < ints_num_threads_; t++) {
        if (do_J_) {
            Jtemps.push_back(std::make_shared<Matrix>("Jtemp", natom, 3));
        }
    }

    for (int block = 0; block < Pstarts.size() - 1; block++) {
        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop = Pstarts[block + 1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop = (Pstop == auxiliary_->nshell() ? naux : auxiliary_->shell(Pstop).function_index());
        int np = pstop - pstart;

        // > Integrals < //
#pragma omp parallel for schedule(dynamic) num_threads(nthread_df)
        for (long int PMN = 0L; PMN < static_cast<long>(NP) * npairs; PMN++) {
            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int P = PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            eri[thread]->compute_shell_deriv1(P, 0, M, N);

            const double* buffer = eri[thread]->buffer();

            int nP = auxiliary_->shell(P).nfunction();
            int cP = auxiliary_->shell(P).ncartesian();
            int aP = auxiliary_->shell(P).ncenter();
            int oP = auxiliary_->shell(P).function_index() - pstart;

            int nM = primary_->shell(M).nfunction();
            int cM = primary_->shell(M).ncartesian();
            int aM = primary_->shell(M).ncenter();
            int oM = primary_->shell(M).function_index();

            int nN = primary_->shell(N).nfunction();
            int cN = primary_->shell(N).ncartesian();
            int aN = primary_->shell(N).ncenter();
            int oN = primary_->shell(N).function_index();

            const auto& buffers = eri[thread]->buffers();
            const double* Px = buffers[0];
            const double* Py = buffers[1];
            const double* Pz = buffers[2];
            const double* Mx = buffers[3];
            const double* My = buffers[4];
            const double* Mz = buffers[5];
            const double* Nx = buffers[6];
            const double* Ny = buffers[7];
            const double* Nz = buffers[8];

            double perm = (M == N ? 1.0 : 2.0);

            double** grad_Jp;

            if (do_J_) {
                grad_Jp = Jtemps[thread]->pointer();
            }

            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        //  J^x = (A|pq)^x d_A Dt_pq
                        if (do_J_) {
                            double Ival = 1.0 * perm * dp[p + oP + pstart] * Dtp[m + oM][n + oN];
                            grad_Jp[aP][0] += Ival * (*Px);
                            grad_Jp[aP][1] += Ival * (*Py);
                            grad_Jp[aP][2] += Ival * (*Pz);
                            grad_Jp[aM][0] += Ival * (*Mx);
                            grad_Jp[aM][1] += Ival * (*My);
                            grad_Jp[aM][2] += Ival * (*Mz);
                            grad_Jp[aN][0] += Ival * (*Nx);
                            grad_Jp[aN][1] += Ival * (*Ny);
                            grad_Jp[aN][2] += Ival * (*Nz);
                        }

                        Px++;
                        Py++;
                        Pz++;
                        Mx++;
                        My++;
                        Mz++;
                        Nx++;
                        Ny++;
                        Nz++;
                    }
                }
            }
        }
    }

    timer_off("eri grad");
    // => Temporary Gradient Reduction <= //

    if (do_J_) {
        for (int t = 0; t < ints_num_threads_; t++) {
            gradients_["Coulomb"]->add(Jtemps[t]);
        }
    }
    timer_off("J Grad");
}

void DFJCOSKGrad::compute_hessian(){
    throw PSIEXCEPTION("DFJCOSK analytic Hessians not implemented.");
}

}
}