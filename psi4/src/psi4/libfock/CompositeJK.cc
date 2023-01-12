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

#include "jk.h"
#include "psi4/libqt/qt.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/electrostatic.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/integral.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/lib3index/dftensor.h"

#include <unordered_set>
#include <vector>
#include <map>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {

Matrix compute_numeric_overlap(const DFTGrid &grid, const std::shared_ptr<BasisSet> &primary) {

    // DOI 10.1063/1.3646921, EQ. 9

    // note that the S_num matrix is defined slightly differently in our code
    // to account for the possibility of negative grid weights
    // here, we use S_num = X_sign*(X_nosign)^T
    // where:
    //   1. X_nosign uses sqrt(abs(w)) instead of sqrt(w) for the X matrix
    //   2. X_sign uses sign(w) * sqrt(abs(w)), where sign returns the sign of w, instead of sqrt(w) for the X matrix

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

        // lambda for returning sign of double
        auto sign = [ ](double val) {
            return (val >= 0.0) ? 1.0 : -1.0;
        };

        // resize the buffer of basis function values
        Matrix X_block_nosign("phi_g,u", npoints_block, nbf_block);  // points x nbf_block
        Matrix X_block_sign("phi_g,u", npoints_block, nbf_block);  // points x nbf_block

        auto X_block_nosignp = X_block_nosign.pointer();
        auto X_block_signp = X_block_sign.pointer();
        for (size_t p = 0; p < npoints_block; p++) {
            for (size_t k = 0; k < nbf_block; k++) {
                X_block_nosignp[p][k] = point_values->get(p, k) * std::sqrt(std::abs(w[p]));
                X_block_signp[p][k] = sign(w[p])*X_block_nosignp[p][k];
            }
        }

        // significant basis functions at these grid points
        const auto &bf_map = block->functions_local_to_global();

        auto S_num_block = linalg::doublet(X_block_sign, X_block_nosign, true, false);
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
CompositeJK::CompositeJK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options) : JK(primary), auxiliary_(auxiliary), options_(options) {
    timer_on("CompositeJK: Setup");
    common_init(); 
    timer_off("CompositeJK: Setup");
}

CompositeJK::~CompositeJK() {}

void CompositeJK::common_init() {

    // => General Setup <= //

    // thread count
    nthreads_ = 1;
#ifdef _OPENMP
    nthreads_ = Process::environment.get_n_threads();
#endif

    // incremental Fock build
    incfock_ = options_.get_bool("INCFOCK");
    incfock_count_ = 0;
    do_incfock_iter_ = false;
    if (options_.get_int("INCFOCK_FULL_FOCK_EVERY") <= 0) {
        throw PSIEXCEPTION("Invalid input for option INCFOCK_FULL_FOCK_EVERY (<= 0)");
    }

    computed_shells_per_iter_["Quartets"] = {};
    
    // derive separate J+K algorithms from scf_type
    auto jk_type = options_.get_str("SCF_TYPE");
    auto j_type = jk_type.substr(0, jk_type.find("+"));
    k_type_ = jk_type.substr(jk_type.find("+") + 1, jk_type.length());

    // occurs if no composite K algorithm was specified; useful for LDA/GGA DFT runs
    if (k_type_ == j_type_) {
      k_type_ = "NONE";
    }

    // other options
    density_screening_ = options_.get_str("SCREENING") == "DENSITY";
    set_cutoff(options_.get_double("INTS_TOLERANCE"));
    early_screening_ = k_type_ == "COSX" ? true : false;

    // pre-construct per-thread TwoBodyAOInt objects for computing 3- and 4-index ERIs
    timer_on("CompositeJK: ERI Computers");
    
    auto zero = BasisSet::zero_ao_basis_set();

    // initialize 4-Center ERIs
    eri_computers_["4-Center"].emplace({});
    eri_computers_["4-Center"].resize(nthreads_);

    IntegralFactory factory(primary_, primary_, primary_, primary_);
    eri_computers_["4-Center"][0] = std::shared_ptr<TwoBodyAOInt>(factory.eri());

    // create each threads' ERI computers
    for(int rank = 1; rank < nthreads_; rank++) {
        eri_computers_["4-Center"][rank] = std::shared_ptr<TwoBodyAOInt>(eri_computers_["4-Center"].front()->clone());
    }

    timer_off("CompositeJK: ERI Computers");

    // => Set up separate J algorithm <= //

    // Direct DF-J
    if (j_type_ == "DFDIRJ") {
	    // initialize SplitJK algo
	    j_algo_ = std::make_shared<DirectDFJ>();

    	// create 3-center ERIs
        eri_computers_["3-Center"].emplace({});
        eri_computers_["3-Center"].resize(nthreads_);

        computed_shells_per_iter_["Triplets"] = {};
        
        IntegralFactory rifactory(auxiliary_, zero, primary_, primary_);
        eri_computers_["3-Center"][0] = std::shared_ptr<TwoBodyAOInt>(rifactory.eri());

        for(int rank = 1; rank < nthreads_; rank++) {
            eri_computers_["3-Center"][rank] = std::shared_ptr<TwoBodyAOInt>(eri_computers_["3-Center"].front()->clone());
        }
    } else {
        throw PSIEXCEPTION("Invalid Composite J algorithm selected!");
    }

    // => Set up separate K algorithm <= //

    // Linear Exchange (LinK)
    if (k_type_ == "LINK") {
	    k_algo_ = std::make_shared<LinK>(primary_, options_);

    // Chain-of-Spheres Exchange (COSX)
    } else if (k_type_ == "COSX") {
	    // null k_algo_ for now
	    k_algo_ = nullptr; 

        timer_on("CompositeJK: COSX Grid Construction");

        // TODO: specify bool "DFT_REMOVE_DISTANT_POINTS" in the DFTGrid constructors

        // Create a small DFTGrid for the initial SCF iterations
        std::map<std::string, std::string> grid_init_str_options = {
            {"DFT_PRUNING_SCHEME", options_.get_str("COSX_PRUNING_SCHEME")},
            {"DFT_RADIAL_SCHEME",  "TREUTLER"},
            {"DFT_NUCLEAR_SCHEME", "TREUTLER"},
            {"DFT_GRID_NAME",      ""},
            {"DFT_BLOCK_SCHEME",   "OCTREE"},
        };
        std::map<std::string, int> grid_init_int_options = {
            {"DFT_SPHERICAL_POINTS", options_.get_int("COSX_SPHERICAL_POINTS_INITIAL")},
            {"DFT_RADIAL_POINTS",    options_.get_int("COSX_RADIAL_POINTS_INITIAL")},
            {"DFT_BLOCK_MIN_POINTS", 100},
            {"DFT_BLOCK_MAX_POINTS", 256},
        };
        std::map<std::string, double> grid_init_float_options = {
            {"DFT_BASIS_TOLERANCE",   options_.get_double("COSX_BASIS_TOLERANCE")},
            {"DFT_BS_RADIUS_ALPHA",   1.0},
            {"DFT_PRUNING_ALPHA",     1.0},
            {"DFT_BLOCK_MAX_RADIUS",  3.0},
            {"DFT_WEIGHTS_TOLERANCE", 1e-15},
        };
        grid_init_ = std::make_shared<DFTGrid>(primary_->molecule(), primary_, grid_init_int_options, grid_init_str_options, grid_init_float_options, options_);

        // Create a large DFTGrid for the final SCF iteration
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

        // Print out warning if grid with negative grid weights is used 
        // Original Nesse COSX formulation does not support negative grid weights
        // which can happen with certain grid configurations
        // the Psi4 COSX implementation is slightly modified to work with negative grid weights
        // See https://github.com/psi4/psi4/issues/2890 for discussion
        auto warning_printed_init = false;
        for (const auto &init_block : grid_init_->blocks()) {
            const auto w = init_block->w();
            for (int ipoint = 0; ipoint < init_block->npoints(); ++ipoint) {
                if (w[ipoint] < 0.0) {
	                outfile->Printf("  INFO: The definition of the current initial grid includes negative weights, which the original COSX formulation does not support!\n    If this is of concern, please choose another initial grid through adjusting either COSX_PRUNING_SCHEME or COSX_SPHERICAL_POINTS_INITIAL.\n\n");
	                warning_printed_init = true;
		            break;
	            }
            }
	        if (warning_printed_init) break;
        }

        auto warning_printed_final = false;
        for (const auto &final_block : grid_final_->blocks()) {
            const auto w = final_block->w();
            for (int ipoint = 0; ipoint < final_block->npoints(); ++ipoint) {
                if (w[ipoint] < 0.0) {
	                outfile->Printf("  INFO: The definition of the current final grid includes negative weights, which the original COSX formulation does not support!\n    If this is of concern, please choose another final grid through adjusting either COSX_PRUNING_SCHEME or COSX_SPHERICAL_POINTS_FINAL.\n\n");
                    warning_printed_final = true;
	                break;
	            }
            }
	        if (warning_printed_final) break;
        }

        timer_off("CompositeJK: COSX Grid Construction");

        // => Overlap Fitting Metric <= //

        // Fit an overlap metric (Q) for both grids to reduce numerical error

        // DOI 10.1063/1.3646921, EQ. 18
        // Note: the above reference defines Q as S_an @ S_num^{-1} @ X
        // Here, Q refers to just S_ @ S_num^{-1} (no X)
        // This Q is contracted with X later to agree with the literature definition

        timer_on("CompositeJK: COSX Numeric Overlap");

        // compute the numeric overlap matrix for each grid
        auto S_num_init = compute_numeric_overlap(*grid_init_, primary_);
        auto S_num_final = compute_numeric_overlap(*grid_final_, primary_ );

        timer_off("CompositeJK: COSX Numeric Overlap");

        timer_on("CompositeJK: COSX Analytic Overlap");

        // compute the analytic overlap matrix
        MintsHelper helper(primary_, options_);
        auto S_an = helper.ao_overlap();

        timer_off("CompositeJK: COSX Analytic Overlap");

        // form the overlap metric (Q) for each grid

        timer_on("CompositeJK: COSX Overlap Metric Solve");

        int nbf = primary_->nbf();
        std::vector<int> ipiv(nbf);

        // solve: Q_init_ = S_an @ S_num_init_^{-1}
        Q_init_ = S_an->clone();
        C_DGESV(nbf, nbf, S_num_init.pointer()[0], nbf, ipiv.data(), Q_init_->pointer()[0], nbf);

        // solve: Q_final_ = S_an @ S_num_final_^{-1}
        Q_final_ = S_an->clone();
        C_DGESV(nbf, nbf, S_num_final.pointer()[0], nbf, ipiv.data(), Q_final_->pointer()[0], nbf);

        timer_off("CompositeJK: COSX Overlap Metric Solve");

    // Do nothing special if no composite K algorithm
    } else if (k_type_ == "NONE") {
        ;
    } else {
        throw PSIEXCEPTION("Invalid Composite K algorithm selected!");
    }
}

void CompositeJK::set_do_K(bool do_K) {
    // if doing K, we need an associated composite K build algorithm
    if (do_K && k_type_ == "NONE") {
        std::string error_message = "No composite K build algorithm was specified, but K matrix is required for current method! Please specify a composite K build algorithm by setting SCF_TYPE to ";
        error_message += j_type_;
        error_message += "+{K_ALGO}.";

        throw PSIEXCEPTION(error_message);
    } else if (!do_K && k_type_ != "NONE") {
        std::string info_message = "  INFO: A K algorithm (";
        info_message += k_type_;
        info_message += ") was specified in SCF_TYPE, but the current method does not use a K matrix!\n";
        info_message += "  Thus, the specified K algorithm will be unused.\n\n";

        outfile->Printf(info_message);
    }

    do_K_ = do_K;
}

size_t CompositeJK::num_computed_shells() {
    return num_computed_shells_;
}

size_t CompositeJK::memory_estimate() {
    return 0;  // Memory is O(N^2), which psi4 counts as effectively 0
}

void CompositeJK::print_header() const {
    std::string screen_type = options_.get_str("SCREENING");
    if (print_) {
        outfile->Printf("  ==> CompositeJK: Mix-and-Match J+K Algorithm Combos <==\n\n");

        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        if (do_J_) outfile->Printf("    J algorithm:       %11s\n", j_algo_->name().c_str());
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        if (do_K_) outfile->Printf("    K algorithm:       %11s\n", k_type_.c_str());
        outfile->Printf("    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:             %11.3E\n", omega_);
        outfile->Printf("    Integrals threads: %11d\n", nthreads_);
        outfile->Printf("    Memory [MiB]:      %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf("    Incremental Fock:  %11s\n", (incfock_ ? "Yes" : "No"));
        outfile->Printf("    Screening Type:    %11s\n", screen_type.c_str());

        if (do_J_) {
            j_algo_->print_header();
	    }
        if (do_K_) {
            if (k_type_ == "LINK") { k_algo_->print_header(); }
            else if (k_type_ == "COSX") { print_COSX_header(); }
        }
        outfile->Printf("\n");
    }
}

void CompositeJK::print_COSX_header() const {
    if (print_) {
        outfile->Printf("\n");
        outfile->Printf("  ==> COSX: Chain-of-Spheres Semi-Numerical K <==\n\n");

        outfile->Printf("    K Screening Cutoff: %11.0E\n", options_.get_double("COSX_INTS_TOLERANCE"));
        outfile->Printf("    K Density Cutoff:   %11.0E\n", options_.get_double("COSX_DENSITY_TOLERANCE"));
        outfile->Printf("    K Basis Cutoff:     %11.0E\n", options_.get_double("COSX_BASIS_TOLERANCE"));
        outfile->Printf("    K Overlap Fitting:  %11s\n", (options_.get_bool("COSX_OVERLAP_FITTING") ? "Yes" : "No"));
    }
}

void CompositeJK::preiterations() {}

void CompositeJK::incfock_setup() {
    if (do_incfock_iter_) {
        auto njk = D_ao_.size();

        // If there is no previous pseudo-density, this iteration is normal
        if (initial_iteration_ || D_prev_.size() != njk) {
            initial_iteration_ = true;

            D_ref_ = D_ao_;
            zero();
        } else { // Otherwise, the iteration is incremental
            for (size_t jki = 0; jki < njk; jki++) {
                D_ref_[jki] = D_ao_[jki]->clone();
                D_ref_[jki]->subtract(D_prev_[jki]);
            }
        }
    } else {
        D_ref_ = D_ao_;
        zero();
    }
}

void CompositeJK::incfock_postiter() {
    // Save a copy of the density for the next iteration
    D_prev_.clear();
    for(auto const &Di : D_ao_) {
        D_prev_.push_back(Di->clone());
    }
}

void CompositeJK::compute_JK() {
    // wK not supported in CompositeJK yet
    // range-separated semi-numerical exchange needs https://github.com/psi4/psi4/pull/2473
    if (do_wK_) throw PSIEXCEPTION("CompositeJK algorithms do not support wK integrals yet!");

    // explicit setup of Incfock for this SCF iteration
    if (incfock_) {
        timer_on("CompositeJK: INCFOCK Preprocessing");

        auto reset = options_.get_int("INCFOCK_FULL_FOCK_EVERY");
        auto incfock_conv = options_.get_double("INCFOCK_CONVERGENCE");
        auto Dnorm = Process::environment.globals["SCF D NORM"];
        // Do IFB on this iteration?
        do_incfock_iter_ = (Dnorm >= incfock_conv) && !initial_iteration_ && (incfock_count_ % reset != reset - 1);

        if (!initial_iteration_ && (Dnorm >= incfock_conv)) incfock_count_ += 1;

        incfock_setup();
        
        timer_off("CompositeJK: INCFOCK Preprocessing");
    } else {
        D_ref_ = D_ao_;
        zero();
    }

    // update ERI engine density matrices for density screening
    if (density_screening_) {
        for (auto eri_computer : eri_computers_["4-Center"]) {
            eri_computer->update_density(D_ref_);
        }
    }

    // => Perform matrix calculations <= //

    // Coulomb Matrix
    if (do_J_) {
        timer_on("CompositeJK: J");

	    j_algo_->build_G_component(D_ref_, J_ao_, eri_computers_["3-Center"]);

        timer_off("CompositeJK: J");
    }

    // Exchange Matrix
    if (do_K_) {
        timer_on("CompositeJK: K");

        // LinK
        if (k_type_ == "LINK") {
	        k_algo_->build_G_component(D_ref_, K_ao_, eri_computers_["4-Center"]);
        // COSX
        } else if (k_type_ == "COSX") {
            build_COSK(D_ref_, K_ao_);
        }

        timer_off("CompositeJK: K");
    }

    // => Finalize Incremental Fock if required <= //

    if (incfock_) {
        timer_on("CompositeJK: INCFOCK Postprocessing");
        incfock_postiter();
        timer_off("CompositeJK: INCFOCK Postprocessing");
    }

    if (initial_iteration_) initial_iteration_ = false;
}

void CompositeJK::postiterations() {}

// build the K matrix using Neeses's Chain-of-Spheres Exchange algorithm
// algorithm is originally proposed in https://doi.org/10.1016/j.chemphys.2008.10.036 
// overlap fitting is discussed in https://doi.org/10.1063/1.3646921
void CompositeJK::build_COSK(std::vector<std::shared_ptr<Matrix>>& D, std::vector<std::shared_ptr<Matrix>>& K) {

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
    auto grid = early_screening_ ? grid_init_ : grid_final_;
    auto Q = early_screening_ ? Q_init_ : Q_final_;

    // => Initialization <= //

    // per-thread ElectrostaticInt object (for computing one-electron "pseudospectral" integrals)
    std::vector<std::shared_ptr<ElectrostaticInt>> int_computers(nthreads_);

    // per-thread BasisFunctions object (for computing basis function values at grid points)
    std::vector<std::shared_ptr<BasisFunctions>> bf_computers(nthreads_);

    // per-thread K Matrix buffers (for accumulating thread contributions to K)
    std::vector<std::vector<SharedMatrix>> KT(njk, std::vector<SharedMatrix>(nthreads_));

    // initialize per-thread objects
    IntegralFactory factory(primary_);
    for(size_t thread = 0; thread < nthreads_; thread++) {
        int_computers[thread] = std::shared_ptr<ElectrostaticInt>(static_cast<ElectrostaticInt *>(factory.electrostatic().release()));
        bf_computers[thread] = std::make_shared<BasisFunctions>(primary_, grid->max_points(), grid->max_functions());
        for(size_t jki = 0; jki < njk; jki++) {
            KT[jki][thread] = std::make_shared<Matrix>(nbf, nbf);
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
    num_computed_shells_ = 0L;
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
        // as discussed in section 3.1 of DOI 10.1016/j.chemphys.2008.10.036
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

        // note that the X matrix is defined slightly differently in our code
        // to account for the possibility of negative grid weights
        // here, we define X using sqrt(abs(w)) instead of sqrt(w)

	// compute basis functions at these grid points
        bf_computers[rank]->compute_functions(block);
        auto point_values = bf_computers[rank]->basis_values()["PHI"];

        // resize the buffer of basis function values
        auto X_block = std::make_shared<Matrix>(npoints_block, nbf_block);  // points x nbf_block

        auto X_blockp = X_block->pointer();
        for (size_t p = 0; p < npoints_block; p++) {
            for (size_t k = 0; k < nbf_block; k++) {
                X_blockp[p][k] = point_values->get(p, k) * std::sqrt(std::abs(w[p]));
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
        // algorithm can be found in Scheme 1 of DOI 10.1016/j.chemphys.2008.10.036
        std::vector<SharedMatrix> G_block(njk);
        for(size_t jki = 0; jki < njk; jki++) {
            G_block[jki] = std::make_shared<Matrix>(nbf_block_all, npoints_block);
        }

        if(rank == 0) timer_on("ESP Integrals");

        const auto & int_buff = int_computers[rank]->buffers()[0];

        // lambda for returning sign of double
        // needed for formation of G
        auto sign = [ ](double val) {
            return (val >= 0.0) ? 1.0 : -1.0;
        };

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
                    // we fold sign(w) into the formation of G to correct for the modified definition of X
                    for(size_t jki = 0; jki < njk; jki++) {
                        auto F_blockp = F_block[jki]->pointer();
                        auto G_blockp = G_block[jki]->pointer();
                        for (size_t nu = nu_start, index = 0; nu < (nu_start + num_nu); ++nu) {
                            for (size_t tau = tau_start; tau < (tau_start + num_tau); ++tau, index++) {
                                G_blockp[nu][g] += sign(w[g]) * int_buff[index] * F_blockp[g][tau];
                                if (symm) G_blockp[tau][g] += sign(w[g]) * int_buff[index] * F_blockp[g][nu];
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

    }

    timer_off("Grid Loop");

    // Reduce per-thread contributions
    for(size_t jki = 0; jki < njk; jki++) {
        for (size_t thread = 0; thread < nthreads_; thread++) {
            K[jki]->add(KT[jki][thread]);
        }
        if (lr_symmetric_) {
            K[jki]->hermitivitize();
        }
    }

    num_computed_shells_ = int_shells_computed;
    if (get_bench()) {
        computed_shells_per_iter_["Quartets"].push_back(num_computed_shells());
    }
}

}  // namespace psi
