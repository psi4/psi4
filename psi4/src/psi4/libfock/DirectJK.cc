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
#include <limits>
#include <sstream>
#include <unordered_set>
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
    density_screening_ = options_.get_str("SCREENING") == "DENSITY";

    computed_shells_per_iter_["Quartets"] = {};
    
    set_cutoff(options_.get_double("INTS_TOLERANCE"));
}
size_t DirectJK::num_computed_shells() { 
    return num_computed_shells_; 
}
size_t DirectJK::memory_estimate() {
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

bool DirectJK::shell_significant(int M, int N, int R, int S, 
    const std::shared_ptr<TwoBodyAOInt> ints, 
    const std::vector<SharedMatrix>& D) 
{
    if (density_screening_) {
        if ((ints != nullptr) && !D.empty()) {
            // Maximum density matrix equation
            double max_density = 0.0;

            std::array<std::string, 2> one_density_jk_references = { "RHF", "RKS" };
            bool is_one_density_ref = std::any_of(
                one_density_jk_references.cbegin(),
                one_density_jk_references.cend(),
                [&](std::string one_density_jk_reference) { return jk_reference_.find(one_density_jk_reference) != std::string::npos; }
            );

            std::array<std::string, 3> two_density_jk_references = { "UHF", "UKS", "ROHF" }; 
            bool is_two_density_ref = std::any_of(
                two_density_jk_references.cbegin(),
                two_density_jk_references.cend(),
                [&](std::string two_density_jk_reference) { return jk_reference_.find(two_density_jk_reference) != std::string::npos; }
            );

            // Equation 6 (RHF/RKS Case)
            if (is_one_density_ref) {
                max_density = std::max({4.0 * ints->shell_pair_max_density(0, M, N), 4.0 * ints->shell_pair_max_density(0, R, S),
                    ints->shell_pair_max_density(0, M, R), ints->shell_pair_max_density(0, M, S),
                    ints->shell_pair_max_density(0, N, R), ints->shell_pair_max_density(0, N, S)});

            // UHF/UKS/ROHF Case
            } else if (is_two_density_ref) { 
                // J-like terms
                double D_MN = ints->shell_pair_max_density(0, M, N) + ints->shell_pair_max_density(1, M, N);
                double D_RS = ints->shell_pair_max_density(0, R, S) + ints->shell_pair_max_density(1, R, S);

                // K-like terms
                double D_MR = std::max(ints->shell_pair_max_density(0, M, R), ints->shell_pair_max_density(1, M, R));
                double D_MS = std::max(ints->shell_pair_max_density(0, M, S), ints->shell_pair_max_density(1, M, S));
                double D_NR = std::max(ints->shell_pair_max_density(0, N, R), ints->shell_pair_max_density(1, N, R));
                double D_NS = std::max(ints->shell_pair_max_density(0, N, S), ints->shell_pair_max_density(1, N, S));

                max_density = std::max({2.0 * D_MN, 2.0 * D_RS, D_MR, D_MS, D_NR, D_NS});
            // throw, because we are using incompatible reference wave function
            } else {
                std::string error_message = "DirectJK shell significance tests being performed with invalid reference wave function type " + jk_reference_ + "! DirectJK::shell_significant() is only allowed with RHF/RKS, UHF/UKS, and ROHF reference wave functions.";
                throw PSIEXCEPTION(error_message);
            }

            return (ints->shell_ceiling2(M, N, R, S) * max_density * max_density >= cutoff_*cutoff_);
        } else {
            throw PSIEXCEPTION("Tests for significant shell quartets using density screening are being conducted, but the ints and/or D variables are undefined. Check your DirectJK::shell_significant function arguments");
        }
    } else {
        if (ints != nullptr) {
            return ints->shell_significant(M, N, R, S);
        } else {
            throw PSIEXCEPTION("Tests for significant shell quartets are being conducted, but the ints variable is undefined. Check your DirectJK::shell_significant function arguments");
        }
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
        for (int thread = 0; thread < df_ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega_)));
            if (density_screening_) ints[thread]->update_density(D_ref_);
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

    // => Task Blocking <= //

    std::vector<int> task_shells;
    std::vector<int> task_starts;

    // > Atomic Blocking < //

    int atomic_ind = -1;
    for (int P = 0; P < nshell; P++) {
        if (primary_->shell(P).ncenter() > atomic_ind) {
            task_starts.push_back(P);
            atomic_ind++;
        }
        task_shells.push_back(P);
    }
    task_starts.push_back(nshell);

    // < End Atomic Blocking > //

    size_t ntask = task_starts.size() - 1;

    std::vector<int> task_offsets;
    task_offsets.push_back(0);
    for (int P2 = 0; P2 < primary_->nshell(); P2++) {
        task_offsets.push_back(task_offsets[P2] + primary_->shell(task_shells[P2]).nfunction());
    }

    size_t max_task = 0L;
    for (size_t task = 0; task < ntask; task++) {
        size_t size = 0L;
        for (int P2 = task_starts[task]; P2 < task_starts[task + 1]; P2++) {
            size += primary_->shell(task_shells[P2]).nfunction();
        }
        max_task = (max_task >= size ? max_task : size);
    }

    if (debug_) {
        outfile->Printf("  ==> DirectJK: Task Blocking <==\n\n");
        for (size_t task = 0; task < ntask; task++) {
            outfile->Printf("  Task: %3zu, Task Start: %4d, Task End: %4d\n", task, task_starts[task],
                            task_starts[task + 1]);
            for (int P2 = task_starts[task]; P2 < task_starts[task + 1]; P2++) {
                int P = task_shells[P2];
                int size = primary_->shell(P).nfunction();
                int off = primary_->shell(P).function_index();
                int off2 = task_offsets[P2];
                outfile->Printf("    Index %4d, Shell: %4d, Size: %4d, Offset: %4d, Offset2: %4d\n", P2, P, size, off,
                                off2);
            }
        }
        outfile->Printf("\n");
    }

    // => Significant Task Pairs (PQ|-style <= //

    std::vector<std::pair<int, int> > task_pairs;
    for (size_t Ptask = 0; Ptask < ntask; Ptask++) {
        for (size_t Qtask = 0; Qtask < ntask; Qtask++) {
            if (Qtask > Ptask) continue;
            bool found = false;
            for (int P2 = task_starts[Ptask]; P2 < task_starts[Ptask + 1]; P2++) {
                for (int Q2 = task_starts[Qtask]; Q2 < task_starts[Qtask + 1]; Q2++) {
                    int P = task_shells[P2];
                    int Q = task_shells[Q2];
                    if (ints[0]->shell_pair_significant(P, Q)) {
                        found = true;
                        task_pairs.push_back(std::pair<int, int>(Ptask, Qtask));
                        break;
                    }
                }
                if (found) break;
            }
        }
    }
    size_t ntask_pair = task_pairs.size();
    size_t ntask_pair2 = ntask_pair * ntask_pair;

    // => Intermediate Buffers <= //

    // Intermediate J buffer per thread
    std::vector<std::vector<SharedMatrix>> JT;
    if (build_J) {
        for (int thread = 0; thread < nthread; thread++) {
            std::vector<SharedMatrix> J2;
            for (size_t ind = 0; ind < D.size(); ind++) {
                // The factor of 2 comes from exploiting ERI permutational symmetry
                J2.push_back(std::make_shared<Matrix>("JT", 2 * max_task, max_task));
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
                K2.push_back(std::make_shared<Matrix>("KT", (lr_symmetric_ ? 4 : 8) * max_task, max_task));
            }
            KT.push_back(K2);
        }
    }
    
    // => Benchmarks <= //

    num_computed_shells_ = 0L;
    size_t computed_shells = 0L;

// ==> Master Task Loop <== //

#pragma omp parallel for num_threads(nthread) schedule(dynamic) reduction(+ : computed_shells)
    for (size_t task = 0L; task < ntask_pair2; task++) {
        size_t task1 = task / ntask_pair;
        size_t task2 = task % ntask_pair;

        int Ptask = task_pairs[task1].first;
        int Qtask = task_pairs[task1].second;
        int Rtask = task_pairs[task2].first;
        int Stask = task_pairs[task2].second;

        // GOTCHA! Thought this should be RStask > PQtask, but
        // H2/3-21G: Task (10|11) gives valid quartets (30|22) and (31|22)
        // This is an artifact that multiple shells on each task allow
        // for for the Ptask's index to possibly trump any RStask pair,
        // regardless of Qtask's index
        if (Rtask > Ptask) continue;

        // printf("Task: %2d %2d %2d %2d\n", Ptask, Qtask, Rtask, Stask);

        int nPtask = task_starts[Ptask + 1] - task_starts[Ptask];
        int nQtask = task_starts[Qtask + 1] - task_starts[Qtask];
        int nRtask = task_starts[Rtask + 1] - task_starts[Rtask];
        int nStask = task_starts[Stask + 1] - task_starts[Stask];

        int P2start = task_starts[Ptask];
        int Q2start = task_starts[Qtask];
        int R2start = task_starts[Rtask];
        int S2start = task_starts[Stask];

        int dPsize = task_offsets[P2start + nPtask] - task_offsets[P2start];
        int dQsize = task_offsets[Q2start + nQtask] - task_offsets[Q2start];
        int dRsize = task_offsets[R2start + nRtask] - task_offsets[R2start];
        int dSsize = task_offsets[S2start + nStask] - task_offsets[S2start];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // => Master shell quartet loops <= //

        bool touched = false;
        for (int P2 = P2start; P2 < P2start + nPtask; P2++) {
            for (int Q2 = Q2start; Q2 < Q2start + nQtask; Q2++) {
                if (Q2 > P2) continue;
                int P = task_shells[P2];
                int Q = task_shells[Q2];
                if (!ints[0]->shell_pair_significant(P, Q)) continue;
                for (int R2 = R2start; R2 < R2start + nRtask; R2++) {
                    for (int S2 = S2start; S2 < S2start + nStask; S2++) {
                        if (S2 > R2) continue;
                        int R = task_shells[R2];
                        int S = task_shells[S2];
                        if (R2 * nshell + S2 > P2 * nshell + Q2) continue;
                        if (!ints[0]->shell_pair_significant(R, S)) continue;
                        if (!shell_significant(P, Q, R, S, ints[thread], D)) continue;

                        // printf("Quartet: %2d %2d %2d %2d\n", P, Q, R, S);
                        // if (thread == 0) timer_on("JK: Ints");
                        if (ints[thread]->compute_shell(P, Q, R, S) == 0)
                            continue;  // No integrals in this shell quartet
                        computed_shells++;
                        // if (thread == 0) timer_off("JK: Ints");

                        const double* buffer = ints[thread]->buffer();

                        int Psize = primary_->shell(P).nfunction();
                        int Qsize = primary_->shell(Q).nfunction();
                        int Rsize = primary_->shell(R).nfunction();
                        int Ssize = primary_->shell(S).nfunction();

                        int Poff = primary_->shell(P).function_index();
                        int Qoff = primary_->shell(Q).function_index();
                        int Roff = primary_->shell(R).function_index();
                        int Soff = primary_->shell(S).function_index();

                        int Poff2 = task_offsets[P2] - task_offsets[P2start];
                        int Qoff2 = task_offsets[Q2] - task_offsets[Q2start];
                        int Roff2 = task_offsets[R2] - task_offsets[R2start];
                        int Soff2 = task_offsets[S2] - task_offsets[S2start];

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
                                    ::memset((void*)JTp[0L * max_task], '\0', dPsize * dQsize * sizeof(double));
                                    ::memset((void*)JTp[1L * max_task], '\0', dRsize * dSsize * sizeof(double));
                                }

                                if (build_K) {
                                    ::memset((void*)KTp[0L * max_task], '\0', dPsize * dRsize * sizeof(double));
                                    ::memset((void*)KTp[1L * max_task], '\0', dPsize * dSsize * sizeof(double));
                                    ::memset((void*)KTp[2L * max_task], '\0', dQsize * dRsize * sizeof(double));
                                    ::memset((void*)KTp[3L * max_task], '\0', dQsize * dSsize * sizeof(double));
                                    if (!lr_symmetric_) {
                                        ::memset((void*)KTp[4L * max_task], '\0', dRsize * dPsize * sizeof(double));
                                        ::memset((void*)KTp[5L * max_task], '\0', dSsize * dPsize * sizeof(double));
                                        ::memset((void*)KTp[6L * max_task], '\0', dRsize * dQsize * sizeof(double));
                                        ::memset((void*)KTp[7L * max_task], '\0', dSsize * dQsize * sizeof(double));
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
                                J1p = JTp[0L * max_task];
                                J2p = JTp[1L * max_task];
                            }

                            if (build_K) {
                                K1p = KTp[0L * max_task];
                                K2p = KTp[1L * max_task];
                                K3p = KTp[2L * max_task];
                                K4p = KTp[3L * max_task];
                                if (!lr_symmetric_) {
                                    K5p = KTp[4L * max_task];
                                    K6p = KTp[5L * max_task];
                                    K7p = KTp[6L * max_task];
                                    K8p = KTp[7L * max_task];
                                }
                            }

                            double prefactor = 1.0;
                            if (P == Q) prefactor *= 0.5;
                            if (R == S) prefactor *= 0.5;
                            if (P == R && Q == S) prefactor *= 0.5;

                            for (int p = 0; p < Psize; p++) {
                                for (int q = 0; q < Qsize; q++) {
                                    for (int r = 0; r < Rsize; r++) {
                                        for (int s = 0; s < Ssize; s++) {
                                            if (build_J) {
                                                J1p[(p + Poff2) * dQsize + q + Qoff2] +=
                                                    prefactor * (Dp[r + Roff][s + Soff] + Dp[s + Soff][r + Roff]) *
                                                    (*buffer2);
                                                J2p[(r + Roff2) * dSsize + s + Soff2] +=
                                                    prefactor * (Dp[p + Poff][q + Qoff] + Dp[q + Qoff][p + Poff]) *
                                                    (*buffer2);
                                            }
                                            
                                            if (build_K) {
                                                K1p[(p + Poff2) * dRsize + r + Roff2] +=
                                                    prefactor * (Dp[q + Qoff][s + Soff]) * (*buffer2);
                                                K2p[(p + Poff2) * dSsize + s + Soff2] +=
                                                    prefactor * (Dp[q + Qoff][r + Roff]) * (*buffer2);
                                                K3p[(q + Qoff2) * dRsize + r + Roff2] +=
                                                    prefactor * (Dp[p + Poff][s + Soff]) * (*buffer2);
                                                K4p[(q + Qoff2) * dSsize + s + Soff2] +=
                                                    prefactor * (Dp[p + Poff][r + Roff]) * (*buffer2);
                                                if (!lr_symmetric_) {
                                                    K5p[(r + Roff2) * dPsize + p + Poff2] +=
                                                        prefactor * (Dp[s + Soff][q + Qoff]) * (*buffer2);
                                                    K6p[(s + Soff2) * dPsize + p + Poff2] +=
                                                        prefactor * (Dp[r + Roff][q + Qoff]) * (*buffer2);
                                                    K7p[(r + Roff2) * dQsize + q + Qoff2] +=
                                                        prefactor * (Dp[s + Soff][p + Poff]) * (*buffer2);
                                                    K8p[(s + Soff2) * dQsize + q + Qoff2] +=
                                                        prefactor * (Dp[r + Roff][p + Poff]) * (*buffer2);
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
                J1p = JTp[0L * max_task];
                J2p = JTp[1L * max_task];
            }

            if (build_K) {
                K1p = KTp[0L * max_task];
                K2p = KTp[1L * max_task];
                K3p = KTp[2L * max_task];
                K4p = KTp[3L * max_task];
                if (!lr_symmetric_) {
                    K5p = KTp[4L * max_task];
                    K6p = KTp[5L * max_task];
                    K7p = KTp[6L * max_task];
                    K8p = KTp[7L * max_task];
                }
            }

            if (build_J) {

                // > J_PQ < //

                for (int P2 = 0; P2 < nPtask; P2++) {
                    for (int Q2 = 0; Q2 < nQtask; Q2++) {
                        int P = task_shells[P2start + P2];
                        int Q = task_shells[Q2start + Q2];
                        int Psize = primary_->shell(P).nfunction();
                        int Qsize = primary_->shell(Q).nfunction();
                        int Poff = primary_->shell(P).function_index();
                        int Qoff = primary_->shell(Q).function_index();
                        int Poff2 = task_offsets[P2 + P2start] - task_offsets[P2start];
                        int Qoff2 = task_offsets[Q2 + Q2start] - task_offsets[Q2start];
                        for (int p = 0; p < Psize; p++) {
                            for (int q = 0; q < Qsize; q++) {
#pragma omp atomic
                                Jp[p + Poff][q + Qoff] += J1p[(p + Poff2) * dQsize + q + Qoff2];
                            }
                        }
                    }
                }

                // > J_RS < //

                for (int R2 = 0; R2 < nRtask; R2++) {
                    for (int S2 = 0; S2 < nStask; S2++) {
                        int R = task_shells[R2start + R2];
                        int S = task_shells[S2start + S2];
                        int Rsize = primary_->shell(R).nfunction();
                        int Ssize = primary_->shell(S).nfunction();
                        int Roff = primary_->shell(R).function_index();
                        int Soff = primary_->shell(S).function_index();
                        int Roff2 = task_offsets[R2 + R2start] - task_offsets[R2start];
                        int Soff2 = task_offsets[S2 + S2start] - task_offsets[S2start];
                        for (int r = 0; r < Rsize; r++) {
                            for (int s = 0; s < Ssize; s++) {
#pragma omp atomic
                                Jp[r + Roff][s + Soff] += J2p[(r + Roff2) * dSsize + s + Soff2];
                            }
                        }
                    }
                }
            }

            if (build_K) {

                // > K_PR < //

                for (int P2 = 0; P2 < nPtask; P2++) {
                    for (int R2 = 0; R2 < nRtask; R2++) {
                        int P = task_shells[P2start + P2];
                        int R = task_shells[R2start + R2];
                        int Psize = primary_->shell(P).nfunction();
                        int Rsize = primary_->shell(R).nfunction();
                        int Poff = primary_->shell(P).function_index();
                        int Roff = primary_->shell(R).function_index();
                        int Poff2 = task_offsets[P2 + P2start] - task_offsets[P2start];
                        int Roff2 = task_offsets[R2 + R2start] - task_offsets[R2start];
                        for (int p = 0; p < Psize; p++) {
                            for (int r = 0; r < Rsize; r++) {
#pragma omp atomic
                                Kp[p + Poff][r + Roff] += K1p[(p + Poff2) * dRsize + r + Roff2];
                                if (!lr_symmetric_) {
#pragma omp atomic
                                    Kp[r + Roff][p + Poff] += K5p[(r + Roff2) * dPsize + p + Poff2];
                                }
                            }
                        }
                    }
                }

                // > K_PS < //

                for (int P2 = 0; P2 < nPtask; P2++) {
                    for (int S2 = 0; S2 < nStask; S2++) {
                        int P = task_shells[P2start + P2];
                        int S = task_shells[S2start + S2];
                        int Psize = primary_->shell(P).nfunction();
                        int Ssize = primary_->shell(S).nfunction();
                        int Poff = primary_->shell(P).function_index();
                        int Soff = primary_->shell(S).function_index();
                        int Poff2 = task_offsets[P2 + P2start] - task_offsets[P2start];
                        int Soff2 = task_offsets[S2 + S2start] - task_offsets[S2start];
                        for (int p = 0; p < Psize; p++) {
                            for (int s = 0; s < Ssize; s++) {
#pragma omp atomic
                                Kp[p + Poff][s + Soff] += K2p[(p + Poff2) * dSsize + s + Soff2];
                                if (!lr_symmetric_) {
#pragma omp atomic
                                    Kp[s + Soff][p + Poff] += K6p[(s + Soff2) * dPsize + p + Poff2];
                                }
                            }
                        }
                    }
                }

                // > K_QR < //

                for (int Q2 = 0; Q2 < nQtask; Q2++) {
                    for (int R2 = 0; R2 < nRtask; R2++) {
                        int Q = task_shells[Q2start + Q2];
                        int R = task_shells[R2start + R2];
                        int Qsize = primary_->shell(Q).nfunction();
                        int Rsize = primary_->shell(R).nfunction();
                        int Qoff = primary_->shell(Q).function_index();
                        int Roff = primary_->shell(R).function_index();
                        int Qoff2 = task_offsets[Q2 + Q2start] - task_offsets[Q2start];
                        int Roff2 = task_offsets[R2 + R2start] - task_offsets[R2start];
                        for (int q = 0; q < Qsize; q++) {
                            for (int r = 0; r < Rsize; r++) {
#pragma omp atomic
                                Kp[q + Qoff][r + Roff] += K3p[(q + Qoff2) * dRsize + r + Roff2];
                                if (!lr_symmetric_) {
#pragma omp atomic
                                    Kp[r + Roff][q + Qoff] += K7p[(r + Roff2) * dQsize + q + Qoff2];
                                }
                            }
                        }
                    }
                }

                // > K_QS < //

                for (int Q2 = 0; Q2 < nQtask; Q2++) {
                    for (int S2 = 0; S2 < nStask; S2++) {
                        int Q = task_shells[Q2start + Q2];
                        int S = task_shells[S2start + S2];
                        int Qsize = primary_->shell(Q).nfunction();
                        int Ssize = primary_->shell(S).nfunction();
                        int Qoff = primary_->shell(Q).function_index();
                        int Soff = primary_->shell(S).function_index();
                        int Qoff2 = task_offsets[Q2 + Q2start] - task_offsets[Q2start];
                        int Soff2 = task_offsets[S2 + S2start] - task_offsets[S2start];
                        for (int q = 0; q < Qsize; q++) {
                            for (int s = 0; s < Ssize; s++) {
#pragma omp atomic
                                Kp[q + Qoff][s + Soff] += K4p[(q + Qoff2) * dSsize + s + Soff2];
                                if (!lr_symmetric_) {
#pragma omp atomic
                                    Kp[s + Soff][q + Qoff] += K8p[(s + Soff2) * dQsize + q + Qoff2];
                                }
                            }
                        }
                    }
                }
            }

        }  // End stripe out
        // if (thread == 0) timer_off("JK: Atomic");

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
