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

#include <memory>
#include <unordered_set>
#include <vector>
#include <map>
#include <algorithm>
#include <cctype>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {

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
    auto k_type = jk_type.substr(jk_type.find("+") + 1, jk_type.length());

    // occurs if no composite K algorithm was specified; useful for LDA/GGA DFT runs
    if (k_type == j_type) {
      k_type = "NONE";
    }

    // other options
    density_screening_ = options_.get_str("SCREENING") == "DENSITY";
    set_cutoff(options_.get_double("INTS_TOLERANCE"));

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

    // DF-DirJ
    if (j_type == "DFDIRJ") {
        // initialize SplitJK algo
        j_algo_ = std::make_shared<DirectDFJ>(primary_, auxiliary_, options_);

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

    // LinK
    if (k_type == "LINK") {
        k_algo_ = std::make_shared<LinK>(primary_, options_);

    // COSX
    } else if (k_type == "COSX") {
        k_algo_ = std::make_shared<COSK>(primary_, options_);

    // sn-LinK (via GauXC) 
    } else if (k_type == "SNLINK") {
        k_algo_ = std::make_shared<snLinK>(primary_, options_);
 
    // No K algorithm specified in SCF_TYPE
    } else if (k_type == "NONE") {
        k_algo_ = nullptr;
    
    // Invalid K algorithm selected
    } else {
        throw PSIEXCEPTION("Invalid Composite K algorithm selected!");
    }
}

void CompositeJK::set_do_K(bool do_K) {
    // if doing K, we need an associated composite K build algorithm
    if (do_K && k_algo_ == nullptr) {
        std::string error_message = "No composite K build algorithm was specified, but K matrix is required for current method! Please specify a composite K build algorithm by setting SCF_TYPE to ";
        error_message += j_algo_->name();
        error_message += "+{K_ALGO}.";

        throw PSIEXCEPTION(error_message);
    } else if (!do_K && k_algo_ != nullptr) {
        std::string info_message = "  INFO: A K algorithm (";
        info_message += k_algo_->name();
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
        if (do_K_) outfile->Printf("    K algorithm:       %11s\n", k_algo_->name().c_str());
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
            k_algo_->print_header();
        }
        outfile->Printf("\n");
    }
}

bool CompositeJK::shell_significant(int M, int N, int R, int S, 
    const std::shared_ptr<TwoBodyAOInt> ints, 
    const std::vector<SharedMatrix>& D) 
{
    throw PSIEXCEPTION("CompositeJK::shell_significant() must be called per-algorithm!");
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

    // set compute()-specific parameters
    j_algo_->set_lr_symmetric(lr_symmetric_);
    
    if (do_K_) {
        k_algo_->set_lr_symmetric(lr_symmetric_);
    }

    // explicit setup of Incfock for this SCF iteration
    if (incfock_) {
        timer_on("CompositeJK: INCFOCK Preprocessing");

        auto reset = options_.get_int("INCFOCK_FULL_FOCK_EVERY");
        auto incfock_conv = options_.get_double("INCFOCK_CONVERGENCE");
        auto Dnorm = Process::environment.globals["SCF D NORM"];
        // Do IFB on this iteration?
        do_incfock_iter_ = (Dnorm >= incfock_conv) && !initial_iteration_ && (incfock_count_ % reset != reset - 1);

        if (k_algo_->name() == "sn-LinK") {
            auto k_algo_derived = std::dynamic_pointer_cast<snLinK>(k_algo_); 
            k_algo_derived->set_incfock_iter(do_incfock_iter_);
        }

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
        timer_on("CompositeJK: " + j_algo_->name());

        j_algo_->build_G_component(D_ref_, J_ao_, eri_computers_["3-Center"]);

        if (get_bench()) {
            computed_shells_per_iter_["Triplets"].push_back(j_algo_->num_computed_shells());
        }
 
        timer_off("CompositeJK: " + j_algo_->name());
    }

    // Exchange Matrix
    if (do_K_) {
        timer_on("CompositeJK: " + k_algo_->name());

        if (k_algo_->name() == "COSX") {
            std::string gridname = get_COSX_grid();
            timer_on("COSX " + gridname + " Grid");
        }

        k_algo_->build_G_component(D_ref_, K_ao_, eri_computers_["4-Center"]);

        if (get_bench()) {
            computed_shells_per_iter_["Quartets"].push_back(k_algo_->num_computed_shells());
        }

        if (k_algo_->name() == "COSX") {
            std::string gridname = get_COSX_grid();
            timer_off("COSX " + gridname + " Grid");
        }

        timer_off("CompositeJK: " + k_algo_->name());
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

// => Method-specific knobs go here <= //

void CompositeJK::set_COSX_grid(std::string current_grid) { 
    if (k_algo_->name() == "COSX") {
        auto k_algo_derived = std::dynamic_pointer_cast<COSK>(k_algo_); 
        k_algo_derived->set_grid(current_grid); 
    } else {
        throw PSIEXCEPTION("CompositeJK::set_COSX_grid() was called, but COSX is not selected in SCF_TYPE!");
    }
}

std::string CompositeJK::get_COSX_grid() { 
    if (k_algo_->name() == "COSX") {
        auto k_algo_derived = std::dynamic_pointer_cast<COSK>(k_algo_); 
        return k_algo_derived->get_grid(); 
    } else {
        throw PSIEXCEPTION("CompositeJK::get_COSX_grid() was called, but COSX is not selected in SCF_TYPE!");
    }
}

int CompositeJK::get_snLinK_max_am() { 
    if (k_algo_->name() == "sn-LinK") {
        auto k_algo_derived = std::dynamic_pointer_cast<snLinK>(k_algo_); 
        return k_algo_derived->get_max_am(); 
    } else {
        throw PSIEXCEPTION("CompositeJK::get_snLinK_max_am() was called, but snLinK is not selected in SCF_TYPE!");
    }
}

}  // namespace psi
