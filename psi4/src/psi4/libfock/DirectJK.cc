/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "psi4/lib3index/3index.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libiwl/iwl.hpp"
#include "jk.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/lib3index/cholesky.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"

#include <algorithm>
#include <span>
#include <utility>
#include <vector>
#include "psi4/libpsi4util/PsiOutStream.h"
#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

#ifdef USING_BrianQC

#include <use_brian_wrapper.h>
#include <brian_macros.h>
#include <brian_common.h>
#include <brian_scf.h>
#include <brian_cphf.h>

extern void checkBrian();
extern BrianCookie brianCookie;
extern bool brianEnable;
extern bool brianEnableDFT;
extern bool brianCPHFFlag;
extern bool brianCPHFLeftSideFlag;
extern brianInt brianRestrictionType;

#endif

using namespace psi;

namespace psi {

DirectJK::DirectJK(std::shared_ptr<BasisSet> primary, Options& options) : JK(primary), options_(options) { common_init(); }
DirectJK::~DirectJK() {}
void DirectJK::common_init() {
    df_ints_num_threads_ = 1;

#ifdef _OPENMP
    df_ints_num_threads_ = Process::environment.get_n_threads();
#endif

    incfock_ = options_.get_bool("INCFOCK");
    incfock_count_ = 0;
    do_incfock_iter_ = false;
    if (options_.get_int("INCFOCK_FULL_FOCK_EVERY") <= 0) {
        throw PSIEXCEPTION("Invalid input for option INCFOCK_FULL_FOCK_EVERY (<= 0)");
    }

    // other options
    auto screening_type = options_.get_str("SCREENING");
    density_screening_ = screening_type == "DENSITY";
    computed_shells_per_iter_["Quartets"] = {};
}

size_t DirectJK::num_computed_shells() { 
    return num_computed_shells_; 
}

size_t DirectJK::memory_estimate() {
    // TODO: return an accurate value.
    return 0;  // Effectively
}

void DirectJK::print_header() const {
    std::string screen_type = options_.get_str("SCREENING");
    if (print_) {
        outfile->Printf("  ==> DirectJK: Integral-Direct J/K Matrices <==\n\n");

        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:             %11.3E\n", omega_);
        outfile->Printf("    Integrals threads: %11d\n", df_ints_num_threads_);
        // outfile->Printf( "    Memory [MiB]:      %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf("    Screening Type:    %11s\n", screen_type.c_str());
        outfile->Printf("    Screening Cutoff:  %11.0E\n", cutoff_);
        outfile->Printf("    Incremental Fock:  %11s\n", incfock_ ? "Yes" : "No");
        outfile->Printf("\n");
    }
}
void DirectJK::preiterations() {

#ifdef USING_BrianQC
    if (brianEnable) {
        double threshold = cutoff_ * (brianCPHFFlag ? 1e-3 : 1e-0); // CPHF needs higher precision
        brianCOMSetPrecisionThresholds(&brianCookie, &threshold);
        checkBrian();
    }
#endif
}

void DirectJK::incfock_setup() {
    if (do_incfock_iter_) {
        size_t njk = D_ao_.size();

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

void DirectJK::incfock_postiter() {
    // Save a copy of the density for the next iteration
    D_prev_.clear();
    for(auto const &Di : D_ao_) {
        D_prev_.push_back(Di->clone());
    }
}

void DirectJK::compute_JK() {
   
#ifdef USING_BrianQC
    if (brianEnable) {
        // zero out J, K, and wK matrices
        zero();
        
        brianBool computeCoulomb = (do_J_ ? BRIAN_TRUE : BRIAN_FALSE);
        brianBool computeExchange = ((do_K_ || do_wK_) ? BRIAN_TRUE : BRIAN_FALSE);

        if (do_wK_ and not brianEnableDFT) {
            throw PSIEXCEPTION("Currently, BrianQC cannot compute range-separated exact exchange when Psi4 is handling the DFT terms");
        }

        if (not brianCPHFFlag) {
            if (!lr_symmetric_) {
                throw PSIEXCEPTION("Currently, BrianQC's non-CPHF Fock building only works with symmetric densities");
            }

            // BrianQC only computes the sum of all Coulomb contributions.
            // For ROHF, the matrices are not the alpha and beta densities, but
            // the doubly and singly occupied densities, and the weight of
            // the first Coulomb contribution must be two. Currently, we
            // achieve this by scaling the doubly occupied density
            // before building, and doing the reverse for the results.
            // We also restore the original density in case it is still needed.
            if (brianRestrictionType == BRIAN_RESTRICTION_TYPE_ROHF) {
                D_ao_[0]->scale(2.0);
            }

            double* exchangeAlpha = nullptr;
            double* exchangeBeta = nullptr;
            if (do_K_) {
                exchangeAlpha = K_ao_[0]->get_pointer();
                exchangeBeta = (D_ao_.size() > 1) ? K_ao_[1]->get_pointer() : nullptr;
            } else if (do_wK_) {
                exchangeAlpha = wK_ao_[0]->get_pointer();
                exchangeBeta = (D_ao_.size() > 1) ? wK_ao_[1]->get_pointer() : nullptr;
            }

            brianSCFBuildFockRepulsion(&brianCookie,
                &computeCoulomb,
                &computeExchange,
                D_ao_[0]->get_pointer(0),
                (D_ao_.size() > 1 ? D_ao_[1]->get_pointer() : nullptr),
                (do_J_ ? J_ao_[0]->get_pointer() : nullptr),
                exchangeAlpha,
                exchangeBeta
            );
            checkBrian();

            // BrianQC computes the sum of all Coulomb contributions into
            // J_ao_[0], so all other contributions must be zeroed out for
            // the sum to be correct. For RHF/RKS, Psi4 expects J_ao_[0]
            // to contain the alpha contribution instead of the total, so
            // we halve it.
            if (do_J_) {
                if (brianRestrictionType == BRIAN_RESTRICTION_TYPE_RHF) {
                    J_ao_[0]->scale(0.5);
                }

                for (size_t ind = 1; ind < J_ao_.size(); ind++) {
                    J_ao_[ind]->zero();
                }
            }

            if (brianRestrictionType == BRIAN_RESTRICTION_TYPE_ROHF) {
                D_ao_[0]->scale(0.5);

                if (do_J_) {
                    J_ao_[0]->scale(0.5);
                }

                if (do_K_) {
                    K_ao_[0]->scale(0.5);
                }

                if (do_wK_) {
                    wK_ao_[0]->scale(0.5);
                }
            }
        } else {
            brianInt maxSegmentSize;
            brianCPHFMaxSegmentSize(&brianCookie, &maxSegmentSize);

            brianInt densityCount = (brianRestrictionType == BRIAN_RESTRICTION_TYPE_RHF) ? 1 : 2;
            if (D_ao_.size() % densityCount != 0) {
                throw PSIEXCEPTION("Invalid number of density matrices for CPHF");
            }

            brianInt derivativeCount = D_ao_.size() / densityCount;

            for (brianInt segmentStartIndex = 0; segmentStartIndex < derivativeCount; segmentStartIndex += maxSegmentSize) {
                brianInt segmentSize = std::min(maxSegmentSize, derivativeCount - segmentStartIndex);

                std::vector<std::vector<SharedMatrix>> pseudoDensitySymmetrized(densityCount);
                std::vector<std::vector<const double*>> pseudoDensityPointers(densityCount);
                std::vector<std::vector<double*>> pseudoExchangePointers(densityCount);
                for (brianInt densityIndex = 0; densityIndex < densityCount; densityIndex++) {
                    pseudoDensitySymmetrized[densityIndex].resize(segmentSize);
                    pseudoDensityPointers[densityIndex].resize(segmentSize, nullptr);
                    pseudoExchangePointers[densityIndex].resize(segmentSize, nullptr);
                    for (brianInt i = 0; i < segmentSize; i++) {
                        // Psi4's code computing the left- and right-hand side CPHF terms use different indexing conventions
                        brianInt psi4Index = brianCPHFLeftSideFlag ? (densityIndex * derivativeCount + segmentStartIndex + i) : ((segmentStartIndex + i) * densityCount + densityIndex);

                        pseudoDensitySymmetrized[densityIndex][i] = D_ao_[psi4Index]->clone();
                        pseudoDensitySymmetrized[densityIndex][i]->add(D_ao_[psi4Index]->transpose());
                        pseudoDensitySymmetrized[densityIndex][i]->scale(0.5);
                        pseudoDensityPointers[densityIndex][i] = pseudoDensitySymmetrized[densityIndex][i]->get_pointer();

                        if (do_K_) {
                            pseudoExchangePointers[densityIndex][i] = K_ao_[psi4Index]->get_pointer();
                        } else if (do_wK_) {
                            pseudoExchangePointers[densityIndex][i] = wK_ao_[psi4Index]->get_pointer();
                        }
                    }
                }

                std::vector<double*> pseudoCoulombPointers(segmentSize, nullptr);
                for (brianInt i = 0; i < segmentSize; i++) {
                    if (do_J_) {
                        // we always write the total coulomb into the densityIndex == 0 matrix, and later divide it if necessary
                        brianInt psi4Index = brianCPHFLeftSideFlag ? (0 * derivativeCount + segmentStartIndex + i) : ((segmentStartIndex + i) * densityCount + 0);
                        pseudoCoulombPointers[i] = J_ao_[psi4Index]->get_pointer();
                    }
                }

                brianCPHFBuildRepulsion(&brianCookie,
                    &computeCoulomb,
                    &computeExchange,
                    &segmentSize,
                    pseudoDensityPointers[0].data(),
                    (densityCount > 1) ? pseudoDensityPointers[1].data() : nullptr,
                    pseudoCoulombPointers.data(),
                    pseudoExchangePointers[0].data(),
                    (densityCount > 1) ? pseudoExchangePointers[1].data() : nullptr
                );
                checkBrian();

                // BrianQC computes the sum of all Coulomb contributions into
                // J_ao_[0], so all other contributions must be zeroed out for
                // the sum to be correct. For RHF/RKS, Psi4 expects J_ao_[0]
                // to contain the alpha contribution instead of the total, so
                // we halve it.
                if (do_J_) {
                    for (brianInt i = 0; i < segmentSize; i++) {
                        if (brianRestrictionType == BRIAN_RESTRICTION_TYPE_RHF) {
                            brianInt psi4Index = brianCPHFLeftSideFlag ? (0 * derivativeCount + segmentStartIndex + i) : ((segmentStartIndex + i) * densityCount + 0);
                            J_ao_[psi4Index]->scale(0.5);
                        }

                        for (brianInt densityIndex = 1; densityIndex < densityCount; densityIndex++) {
                            brianInt psi4Index = brianCPHFLeftSideFlag ? (densityIndex * derivativeCount + segmentStartIndex + i) : ((segmentStartIndex + i) * densityCount + densityIndex);
                            J_ao_[psi4Index]->zero();
                        }
                    }
                }
            }
        }

        return;
    }
#endif

    if (incfock_) {
        timer_on("DirectJK: INCFOCK Preprocessing");
        int reset = options_.get_int("INCFOCK_FULL_FOCK_EVERY");
        double incfock_conv = options_.get_double("INCFOCK_CONVERGENCE");
        double Dnorm = Process::environment.globals["SCF D NORM"];
        // Do IFB on this iteration?
        do_incfock_iter_ = (Dnorm >= incfock_conv) && !initial_iteration_ && (incfock_count_ % reset != reset - 1);
        
        if (!initial_iteration_ && (Dnorm >= incfock_conv)) incfock_count_ += 1;
        
        incfock_setup();
	
        timer_off("DirectJK: INCFOCK Preprocessing");
    } else {
        D_ref_ = D_ao_;
        zero();
    }

    auto factory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);
    
    // Passed in as a dummy when J (and/or K) is not built
    std::vector<SharedMatrix> temp;

    if (do_wK_) {
        std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega_)));
        if (density_screening_) ints[0]->update_density(D_ref_);
        for (int thread = 1; thread < df_ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(ints[0]->clone()));
        }
        if (do_J_) {
            build_JK_matrices(ints, D_ref_, J_ao_, wK_ao_);
        } else {
            build_JK_matrices(ints, D_ref_, temp, wK_ao_);
        }
    }

    if (do_J_ || do_K_) {
        std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri()));
        if (density_screening_) ints[0]->update_density(D_ref_);
        for (int thread = 1; thread < df_ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(ints[0]->clone()));
        }
        if (do_J_ && do_K_) {
            build_JK_matrices(ints, D_ref_, J_ao_, K_ao_);
        } else if (do_J_) {
            build_JK_matrices(ints, D_ref_, J_ao_, temp);
        } else {
            build_JK_matrices(ints, D_ref_, temp, K_ao_);
        }
    }

    if (incfock_) {
        timer_on("DirectJK: INCFOCK Postprocessing");
        incfock_postiter();
        timer_off("DirectJK: INCFOCK Postprocessing");
    }

    if (initial_iteration_) initial_iteration_ = false;
}
void DirectJK::postiterations() {}

namespace {

// Minimal integer range [first, last) used in range-for loops. This replaces
// std::views::iota, which fails to compile under clang with libstdc++ < 13.
// Whenever Psi4's minimum GCC version is bumped to at least 14, we should
// revert back to using std::views::iota.
struct integer_range {
    std::size_t first;
    std::size_t last;

    struct iterator {
        std::size_t value;
        std::size_t operator*() const { return value; }
        iterator& operator++() { ++value; return *this; }
        bool operator==(const iterator& other) const { return value == other.value; }
        bool operator!=(const iterator& other) const { return value != other.value; }
    };

    iterator begin() const { return {first}; }
    iterator end() const { return {last}; }
};

// Helper method to remove the indexing madness. Produces iterable list of iterators.
// Converts Something like [0,1,4,5] into [[0, 1), [1, 4), [4, 5)].
auto partition(std::span<const std::size_t> atom_to_shell) {
    std::vector<integer_range> result;
    result.reserve(atom_to_shell.size());
    for (std::size_t atom = 0; atom + 1 < atom_to_shell.size(); ++atom)
        result.push_back(integer_range{atom_to_shell[atom], atom_to_shell[atom + 1]});
    return result;
}

// Like above but enumerated. Converts Something like [0,1,4,5] into
//   [{0, [0, 1)}, {1, [1, 4)}, {2, [4, 5)}].
auto partition_with_idx(std::span<const std::size_t> atom_to_shell) {
    std::vector<std::pair<std::size_t, integer_range>> result;
    result.reserve(atom_to_shell.size());
    for (std::size_t atom = 0; atom + 1 < atom_to_shell.size(); ++atom)
        result.emplace_back(atom, integer_range{atom_to_shell[atom], atom_to_shell[atom + 1]});
    return result;
}

}  // namespace

void DirectJK::build_JK_matrices(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints, const std::vector<SharedMatrix>& D,
                        std::vector<SharedMatrix>& J, std::vector<SharedMatrix>& K) {

    bool build_J = (!J.empty());
    bool build_K = (!K.empty());

    if (!build_J && !build_K) return;
    
    timer_on("build_JK_matrices()");

    // => Zeroing... <= //
    
    // Ideally, this wouldnt be here at all
    // It would be better covered in incfock_setup()
    // But removing this causes a couple of tests to fail for some reason
    
    if (!do_incfock_iter_) {
        for (auto& Jmat : J) {
            Jmat->zero();
        }
    
        for (auto& Kmat : K) {
            Kmat->zero();
        }
    }

    // => Sizing <= //

    int nshell = primary_->nshell();
    int nthread = df_ints_num_threads_;

    // => Atomic Task Blocking <= //
    // One task = all shells on one center. Shells stay in basis order (identity map).

    // shell index that starts each center
    std::vector<size_t> atom_to_shell;
    // basis function index that starts each shell
    std::vector<size_t> shell_to_basis;

    /*******************************************************
     * Consider H2O with STO-3G...                         *
     *                                                     *
     *   3 atoms:       H  ,         O           ,   H     *
     *   5 shells:    [ S ], [ S , S ,    P     ], [ S ]   *
     *   7 functions: [[s]], [[s],[s],[px,py,pz]], [[s]]   *
     *                                                     *
     * atom_to_shell  [ 0,           1,              4, 5] *
     * shell_to_basis [ 0,     1,  2,     3,         6, 7] *
     *                                                     *
     *******************************************************/

    // => Helper functions that need internal variables <= //

    auto nfunctions_in_shell = [this](const size_t& shell_idx) {
        return this->primary_->shell(shell_idx).nfunction();
    };

    auto function_index_of_shell = [this](const size_t& shell_idx) {
        return this->primary_->shell(shell_idx).function_index();
    };

    // Returns the basis function index that begins shell ``shell_idx`` where 0
    // indicates the first basis function of the first shell of atom ``task``.
    auto basis_index_of_shell_from_atom = [&atom_to_shell, &shell_to_basis](const size_t& shell_idx, const size_t& task) {
        return shell_to_basis[shell_idx] - shell_to_basis[atom_to_shell[task]];
    };

    // Number of basis functions on center ``task``.
    auto nfunctions_on_center = [&atom_to_shell, &shell_to_basis](const size_t& task) {
        return shell_to_basis[atom_to_shell[task + 1]] - shell_to_basis[atom_to_shell[task]];
    };

    // Returns a view of shell indices for given center idx.
    auto shells_on_center = [&atom_to_shell](size_t task) {
        return integer_range{atom_to_shell[task], atom_to_shell[task + 1]};
    };

    // => Welcome to the jungle <= //

    {
        shell_to_basis.push_back(0);

        size_t total_nfuncs = 0;
        int atomic_ind = -1;
        for (int P = 0; P < nshell; P++) {
            const auto& shell = primary_->shell(P);

            total_nfuncs += shell.nfunction();
            shell_to_basis.push_back(total_nfuncs);

            if (shell.ncenter() > atomic_ind) {
                atom_to_shell.push_back(P);
                atomic_ind++;
            }
        }
        atom_to_shell.push_back(nshell);
    }

    // Largest number of functions for any center. Used to allocate the
    // temporary J and K contraction buffers.
    size_t max_nfuncs_per_center = 0;
    for (auto centers : partition(atom_to_shell)) {
        size_t size = 0;
        for (size_t shell : centers) {
            size += primary_->shell(shell).nfunction();
        }
        max_nfuncs_per_center = std::max(max_nfuncs_per_center, size);
    }

    if (debug_) {
        outfile->Printf("  ==> DirectJK: Task Blocking <==\n\n");
        for (auto [atom, shells] : partition_with_idx(atom_to_shell)) {
            outfile->Printf("  Task: %3zu, Task Start: %4zu, Task End: %4zu\n", atom, atom_to_shell[atom],
                            atom_to_shell[atom + 1]);
            for (size_t shell : shells) {
                int size = primary_->shell(shell).nfunction();
                int off = primary_->shell(shell).function_index();
                int off2 = shell_to_basis[shell];
                outfile->Printf("    Index %4zu, Shell: %4zu, Size: %4d, Offset: %4d, Offset2: %4d\n", shell, shell,
                                size, off, off2);
            }
        }
        outfile->Printf("\n");
    }

    // => Significant Task Pairs (PQ|-style <= //

    std::vector<std::pair<size_t, size_t>> significant_pairs;
    for (auto [Patom, Pshells] : partition_with_idx(atom_to_shell)) {
        for (auto [Qatom, Qshells] : partition_with_idx(atom_to_shell)) {
            if (Qatom > Patom) continue;
            bool found_significant_pair = false;
            for (size_t Ps : Pshells) {
                for (size_t Qs : Qshells) {
                    if (ints[0]->shell_pair_significant(Ps, Qs)) {
                        found_significant_pair = true;
                        significant_pairs.emplace_back(Patom, Qatom);
                        break;
                    }
                }
                if (found_significant_pair) break;
            }
        }
    }
    size_t n_pair = significant_pairs.size();

    // => Intermediate Buffers <= //

    // Intermediate J buffer per thread
    std::vector<std::vector<SharedMatrix>> JT;
    if (build_J) {
        for (int thread = 0; thread < nthread; thread++) {
            std::vector<SharedMatrix> J2;
            for (size_t ind = 0; ind < D.size(); ind++) {
                // The factor of 2 comes from exploiting ERI permutational symmetry
                J2.push_back(std::make_shared<Matrix>("JT", 2 * max_nfuncs_per_center, max_nfuncs_per_center));
            }
            JT.push_back(J2);
        }
    }

    // Intermediate K buffer per thread
    std::vector<std::vector<SharedMatrix>> KT;
    if (build_K) {
        for (int thread = 0; thread < nthread; thread++) {
            std::vector<SharedMatrix > K2;
            for (size_t ind = 0; ind < D.size(); ind++) {
                // The factor of 4 or 8 comes from exploiting ERI permutational symmetry
                K2.push_back(std::make_shared<Matrix>("KT", (lr_symmetric_ ? 4 : 8) * max_nfuncs_per_center,
                                                      max_nfuncs_per_center));
            }
            KT.push_back(K2);
        }
    }
    
    // => Benchmarks <= //

    num_computed_shells_ = 0L;
    size_t computed_shells = 0L;

// ==> Master Task Loop <== //

#pragma omp parallel for num_threads(nthread) schedule(dynamic) collapse(2) reduction(+ : computed_shells)
    for (size_t pq_pair = 0; pq_pair < n_pair; ++pq_pair) {
    for (size_t rs_pair = 0; rs_pair < n_pair; ++rs_pair) {
        auto [Patom, Qatom] = significant_pairs[pq_pair];
        auto [Ratom, Satom] = significant_pairs[rs_pair];

        // GOTCHA! Thought this should be RStask > PQtask, but
        // H2/3-21G: Task (10|11) gives valid quartets (30|22) and (31|22)
        // This is an artifact that multiple shells on each task allow
        // for for the Ptask's index to possibly trump any RStask pair,
        // regardless of Qtask's index
        if (Ratom > Patom) continue;

        // printf("Task: %2zu %2zu %2zu %2zu\n", Patom, Qatom, Ratom, Satom);

        size_t Patom_size = nfunctions_on_center(Patom);
        size_t Qatom_size = nfunctions_on_center(Qatom);
        size_t Ratom_size = nfunctions_on_center(Ratom);
        size_t Satom_size = nfunctions_on_center(Satom);

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // => Master shell quartet loops <= //

        bool touched = false;
        for (auto Ps : shells_on_center(Patom)) {
            for (auto Qs : shells_on_center(Qatom)) {
                if (Qs > Ps) continue;
                if (!ints[0]->shell_pair_significant(Ps, Qs)) continue;
                for (auto Rs : shells_on_center(Ratom)) {
                    for (auto Ss : shells_on_center(Satom)) {
                        if (Ss > Rs) continue;
                        if (Rs * nshell + Ss > Ps * nshell + Qs) continue;
                        if (!ints[0]->shell_pair_significant(Rs, Ss)) continue;
                        if (!ints[0]->shell_significant(Ps, Qs, Rs, Ss)) continue;

                        // printf("Quartet: %2zu %2zu %2zu %2zu\n", Ps, Qs, Rs, Ss);
                        // if (thread == 0) timer_on("JK: Ints");
                        if (ints[thread]->compute_shell(Ps, Qs, Rs, Ss) == 0)
                            continue;  // No integrals in this shell quartet
                        computed_shells++;
                        // if (thread == 0) timer_off("JK: Ints");

                        const double* buffer = ints[thread]->buffer();

                        size_t Psize = nfunctions_in_shell(Ps);
                        size_t Qsize = nfunctions_in_shell(Qs);
                        size_t Rsize = nfunctions_in_shell(Rs);
                        size_t Ssize = nfunctions_in_shell(Ss);

                        size_t Pao = function_index_of_shell(Ps);
                        size_t Qao = function_index_of_shell(Qs);
                        size_t Rao = function_index_of_shell(Rs);
                        size_t Sao = function_index_of_shell(Ss);

                        size_t Plo = basis_index_of_shell_from_atom(Ps, Patom);
                        size_t Qlo = basis_index_of_shell_from_atom(Qs, Qatom);
                        size_t Rlo = basis_index_of_shell_from_atom(Rs, Ratom);
                        size_t Slo = basis_index_of_shell_from_atom(Ss, Satom);

                        // if (thread == 0) timer_on("JK: GEMV");
                        for (size_t ind = 0; ind < D.size(); ind++) {
                            double** Dp = D[ind]->pointer();
                            double** JTp; 
                            if (build_J) JTp = JT[thread][ind]->pointer();
                            double** KTp;
                            if (build_K) KTp = KT[thread][ind]->pointer();
                            const double* buffer2 = buffer;

                            if (!touched) {
                                if (build_J) {
                                    ::memset((void*)JTp[0L * max_nfuncs_per_center], '\0', Patom_size * Qatom_size * sizeof(double));
                                    ::memset((void*)JTp[1L * max_nfuncs_per_center], '\0', Ratom_size * Satom_size * sizeof(double));
                                }

                                if (build_K) {
                                    ::memset((void*)KTp[0L * max_nfuncs_per_center], '\0', Patom_size * Ratom_size * sizeof(double));
                                    ::memset((void*)KTp[1L * max_nfuncs_per_center], '\0', Patom_size * Satom_size * sizeof(double));
                                    ::memset((void*)KTp[2L * max_nfuncs_per_center], '\0', Qatom_size * Ratom_size * sizeof(double));
                                    ::memset((void*)KTp[3L * max_nfuncs_per_center], '\0', Qatom_size * Satom_size * sizeof(double));
                                    if (!lr_symmetric_) {
                                        ::memset((void*)KTp[4L * max_nfuncs_per_center], '\0', Ratom_size * Patom_size * sizeof(double));
                                        ::memset((void*)KTp[5L * max_nfuncs_per_center], '\0', Satom_size * Patom_size * sizeof(double));
                                        ::memset((void*)KTp[6L * max_nfuncs_per_center], '\0', Ratom_size * Qatom_size * sizeof(double));
                                        ::memset((void*)KTp[7L * max_nfuncs_per_center], '\0', Satom_size * Qatom_size * sizeof(double));
                                    }
                                }
                            }

                            // Intermediate Contraction Pointers
                            double* J1p;
                            double* J2p;
                            double* K1p;
                            double* K2p;
                            double* K3p;
                            double* K4p;
                            double* K5p;
                            double* K6p;
                            double* K7p;
                            double* K8p;

                            if (build_J) {
                                J1p = JTp[0L * max_nfuncs_per_center];
                                J2p = JTp[1L * max_nfuncs_per_center];
                            }

                            if (build_K) {
                                K1p = KTp[0L * max_nfuncs_per_center];
                                K2p = KTp[1L * max_nfuncs_per_center];
                                K3p = KTp[2L * max_nfuncs_per_center];
                                K4p = KTp[3L * max_nfuncs_per_center];
                                if (!lr_symmetric_) {
                                    K5p = KTp[4L * max_nfuncs_per_center];
                                    K6p = KTp[5L * max_nfuncs_per_center];
                                    K7p = KTp[6L * max_nfuncs_per_center];
                                    K8p = KTp[7L * max_nfuncs_per_center];
                                }
                            }

                            double prefactor = 1.0;
                            if (Ps == Qs) prefactor *= 0.5;
                            if (Rs == Ss) prefactor *= 0.5;
                            if (Ps == Rs && Qs == Ss) prefactor *= 0.5;

                            for (size_t p = 0; p < Psize; p++) {
                                for (size_t q = 0; q < Qsize; q++) {
                                    for (size_t r = 0; r < Rsize; r++) {
                                        for (size_t s = 0; s < Ssize; s++) {
                                            if (build_J) {
                                                J1p[(p + Plo) * Qatom_size + q + Qlo] +=
                                                    prefactor * (Dp[r + Rao][s + Sao] + Dp[s + Sao][r + Rao]) *
                                                    (*buffer2);
                                                J2p[(r + Rlo) * Satom_size + s + Slo] +=
                                                    prefactor * (Dp[p + Pao][q + Qao] + Dp[q + Qao][p + Pao]) *
                                                    (*buffer2);
                                            }

                                            if (build_K) {
                                                K1p[(p + Plo) * Ratom_size + r + Rlo] +=
                                                    prefactor * (Dp[q + Qao][s + Sao]) * (*buffer2);
                                                K2p[(p + Plo) * Satom_size + s + Slo] +=
                                                    prefactor * (Dp[q + Qao][r + Rao]) * (*buffer2);
                                                K3p[(q + Qlo) * Ratom_size + r + Rlo] +=
                                                    prefactor * (Dp[p + Pao][s + Sao]) * (*buffer2);
                                                K4p[(q + Qlo) * Satom_size + s + Slo] +=
                                                    prefactor * (Dp[p + Pao][r + Rao]) * (*buffer2);
                                                if (!lr_symmetric_) {
                                                    K5p[(r + Rlo) * Patom_size + p + Plo] +=
                                                        prefactor * (Dp[s + Sao][q + Qao]) * (*buffer2);
                                                    K6p[(s + Slo) * Patom_size + p + Plo] +=
                                                        prefactor * (Dp[r + Rao][q + Qao]) * (*buffer2);
                                                    K7p[(r + Rlo) * Qatom_size + q + Qlo] +=
                                                        prefactor * (Dp[s + Sao][p + Pao]) * (*buffer2);
                                                    K8p[(s + Slo) * Qatom_size + q + Qlo] +=
                                                        prefactor * (Dp[r + Rao][p + Pao]) * (*buffer2);
                                                }
                                            }
                                            
                                            buffer2++;
                                        }
                                    }
                                }
                            }
                        }
                        touched = true;
                        // if (thread == 0) timer_off("JK: GEMV");
                    }
                }
            }
        }  // End Shell Quartets

        if (!touched) continue;

        // => Stripe out <= //
        if (build_J) {
	    for (auto& JTmat : JT[thread]) {
                JTmat->scale(2.0);
	    }
        }
        
        if (build_K && lr_symmetric_) {
	    for (auto& KTmat : KT[thread]) {
                KTmat->scale(2.0);
	    }
        }

        // if (thread == 0) timer_on("JK: Atomic");
        for (size_t ind = 0; ind < D.size(); ind++) {
            double** JTp;
            double** KTp;
            double** Jp;
            double** Kp;

            if (build_J) {
                JTp = JT[thread][ind]->pointer();
                Jp = J[ind]->pointer();
            }
            
            if (build_K) {
                KTp = KT[thread][ind]->pointer();
                Kp = K[ind]->pointer();
            }

            double* J1p;
            double* J2p;
            double* K1p;
            double* K2p;
            double* K3p;
            double* K4p;
            double* K5p;
            double* K6p;
            double* K7p;
            double* K8p;

            if (build_J) {
                J1p = JTp[0L * max_nfuncs_per_center];
                J2p = JTp[1L * max_nfuncs_per_center];
            }

            if (build_K) {
                K1p = KTp[0L * max_nfuncs_per_center];
                K2p = KTp[1L * max_nfuncs_per_center];
                K3p = KTp[2L * max_nfuncs_per_center];
                K4p = KTp[3L * max_nfuncs_per_center];
                if (!lr_symmetric_) {
                    K5p = KTp[4L * max_nfuncs_per_center];
                    K6p = KTp[5L * max_nfuncs_per_center];
                    K7p = KTp[6L * max_nfuncs_per_center];
                    K8p = KTp[7L * max_nfuncs_per_center];
                }
            }

            if (build_J) {

                // > J_PQ < //

                for (auto Ps : shells_on_center(Patom)) {
                    for (auto Qs : shells_on_center(Qatom)) {
                        size_t Psize = nfunctions_in_shell(Ps);
                        size_t Qsize = nfunctions_in_shell(Qs);
                        size_t Pao = function_index_of_shell(Ps);
                        size_t Qao = function_index_of_shell(Qs);
                        size_t Plo = basis_index_of_shell_from_atom(Ps, Patom);
                        size_t Qlo = basis_index_of_shell_from_atom(Qs, Qatom);
                        for (int p = 0; p < Psize; p++) {
                            for (int q = 0; q < Qsize; q++) {
#pragma omp atomic
                                Jp[p + Pao][q + Qao] += J1p[(p + Plo) * Qatom_size + q + Qlo];
                            }
                        }
                    }
                }

                // > J_RS < //

                for (auto Rs : shells_on_center(Ratom)) {
                    for (auto Ss : shells_on_center(Satom)) {
                        size_t Rsize = nfunctions_in_shell(Rs);
                        size_t Ssize = nfunctions_in_shell(Ss);
                        size_t Rao = function_index_of_shell(Rs);
                        size_t Sao = function_index_of_shell(Ss);
                        size_t Rlo = basis_index_of_shell_from_atom(Rs, Ratom);
                        size_t Slo = basis_index_of_shell_from_atom(Ss, Satom);
                        for (int r = 0; r < Rsize; r++) {
                            for (int s = 0; s < Ssize; s++) {
#pragma omp atomic
                                Jp[r + Rao][s + Sao] += J2p[(r + Rlo) * Satom_size + s + Slo];
                            }
                        }
                    }
                }
            }

            if (build_K) {

                // > K_PR < //

                for (auto Ps : shells_on_center(Patom)) {
                    for (auto Rs : shells_on_center(Ratom)) {
                        size_t Psize = nfunctions_in_shell(Ps);
                        size_t Rsize = nfunctions_in_shell(Rs);
                        size_t Pao = function_index_of_shell(Ps);
                        size_t Rao = function_index_of_shell(Rs);
                        size_t Plo = basis_index_of_shell_from_atom(Ps, Patom);
                        size_t Rlo = basis_index_of_shell_from_atom(Rs, Ratom);
                        for (int p = 0; p < Psize; p++) {
                            for (int r = 0; r < Rsize; r++) {
#pragma omp atomic
                                Kp[p + Pao][r + Rao] += K1p[(p + Plo) * Ratom_size + r + Rlo];
                                if (!lr_symmetric_) {
#pragma omp atomic
                                    Kp[r + Rao][p + Pao] += K5p[(r + Rlo) * Patom_size + p + Plo];
                                }
                            }
                        }
                    }
                }

                // > K_PS < //

                for (auto Ps : shells_on_center(Patom)) {
                    for (auto Ss : shells_on_center(Satom)) {
                        size_t Psize = nfunctions_in_shell(Ps);
                        size_t Ssize = nfunctions_in_shell(Ss);
                        size_t Pao = function_index_of_shell(Ps);
                        size_t Sao = function_index_of_shell(Ss);
                        size_t Plo = basis_index_of_shell_from_atom(Ps, Patom);
                        size_t Slo = basis_index_of_shell_from_atom(Ss, Satom);
                        for (int p = 0; p < Psize; p++) {
                            for (int s = 0; s < Ssize; s++) {
#pragma omp atomic
                                Kp[p + Pao][s + Sao] += K2p[(p + Plo) * Satom_size + s + Slo];
                                if (!lr_symmetric_) {
#pragma omp atomic
                                    Kp[s + Sao][p + Pao] += K6p[(s + Slo) * Patom_size + p + Plo];
                                }
                            }
                        }
                    }
                }

                // > K_QR < //

                for (auto Qs : shells_on_center(Qatom)) {
                    for (auto Rs : shells_on_center(Ratom)) {
                        size_t Qsize = nfunctions_in_shell(Qs);
                        size_t Rsize = nfunctions_in_shell(Rs);
                        size_t Qao = function_index_of_shell(Qs);
                        size_t Rao = function_index_of_shell(Rs);
                        size_t Qlo = basis_index_of_shell_from_atom(Qs, Qatom);
                        size_t Rlo = basis_index_of_shell_from_atom(Rs, Ratom);
                        for (int q = 0; q < Qsize; q++) {
                            for (int r = 0; r < Rsize; r++) {
#pragma omp atomic
                                Kp[q + Qao][r + Rao] += K3p[(q + Qlo) * Ratom_size + r + Rlo];
                                if (!lr_symmetric_) {
#pragma omp atomic
                                    Kp[r + Rao][q + Qao] += K7p[(r + Rlo) * Qatom_size + q + Qlo];
                                }
                            }
                        }
                    }
                }

                // > K_QS < //

                for (auto Qs : shells_on_center(Qatom)) {
                    for (auto Ss : shells_on_center(Satom)) {
                        size_t Qsize = nfunctions_in_shell(Qs);
                        size_t Ssize = nfunctions_in_shell(Ss);
                        size_t Qao = function_index_of_shell(Qs);
                        size_t Sao = function_index_of_shell(Ss);
                        size_t Qlo = basis_index_of_shell_from_atom(Qs, Qatom);
                        size_t Slo = basis_index_of_shell_from_atom(Ss, Satom);
                        for (int q = 0; q < Qsize; q++) {
                            for (int s = 0; s < Ssize; s++) {
#pragma omp atomic
                                Kp[q + Qao][s + Sao] += K4p[(q + Qlo) * Satom_size + s + Slo];
                                if (!lr_symmetric_) {
#pragma omp atomic
                                    Kp[s + Sao][q + Qao] += K8p[(s + Slo) * Qatom_size + q + Qlo];
                                }
                            }
                        }
                    }
                }
            }

        }  // End stripe out
        // if (thread == 0) timer_off("JK: Atomic");

    }
    }  // End master task list

    for (auto& Jmat : J) {
        Jmat->hermitivitize();
    }

    if (lr_symmetric_) {
        for (auto& Kmat : K) {
            Kmat->hermitivitize();
        }
    }

    num_computed_shells_ = computed_shells;
    if (get_bench()) {
        computed_shells_per_iter_["Quartets"].push_back(num_computed_shells());
    }

    timer_off("build_JK_matrices()");
}

}  // namespace psi
