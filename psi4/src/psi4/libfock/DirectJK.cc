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

    density_screening_ = options_.get_str("SCREENING") == "DENSITY";
    linK_ = options_.get_bool("DO_LINK");

    if (options_["LINK_INTS_TOLERANCE"].has_changed()) {
        linK_ints_cutoff_ = options_.get_double("LINK_INTS_TOLERANCE");
    } else {
        linK_ints_cutoff_ = options_.get_double("INTS_TOLERANCE");
    }
    
    set_cutoff(options_.get_double("INTS_TOLERANCE"));
}
size_t DirectJK::num_computed_shells() { 
    if (linK_) {
	//no bench data returned if LinK is enabled - to come in a future update
	return JK::num_computed_shells(); 
    } else {
	return num_computed_shells_; 
    }
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
        outfile->Printf("    LinK:              %11s\n", linK_ ? "Yes" : "No");
        outfile->Printf("\n");
    }
    if (linK_) outfile->Printf("    WARNING: LinK is still under development and should not be used!\n\n");
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

void DirectJK::compute_JK() {

    // zero out J, K, and wK matrices
    zero();

#ifdef USING_BrianQC
    if (brianEnable) {
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

    if (incfock_) incfock_setup();

    auto factory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);
    
    std::vector<SharedMatrix>& D_ref = (perform_incfock_ ? delta_D_ao_ : D_ao_);
    std::vector<SharedMatrix>& J_ref = (perform_incfock_ ? delta_J_ao_ : J_ao_);
    std::vector<SharedMatrix>& K_ref = (perform_incfock_ ? delta_K_ao_ : K_ao_);
    std::vector<SharedMatrix>& wK_ref = (perform_incfock_ ? delta_wK_ao_ : wK_ao_);

    // Passed in as a dummy when J (and/or K) is not built
    std::vector<SharedMatrix> temp;

    if (do_wK_) {
        std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
        for (int thread = 0; thread < df_ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega_)));
            if (density_screening_) ints[thread]->update_density(D_ref);
        }
        if (do_J_) {
            build_JK_matrices(ints, D_ref, J_ref, wK_ref);
        } else {
            build_JK_matrices(ints, D_ref, temp, wK_ref);
        }
    }

    if (do_J_ || do_K_) {
        std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri()));
        if (density_screening_) ints[0]->update_density(D_ref);
        for (int thread = 1; thread < df_ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(ints[0]->clone()));
        }
        if (do_J_ && do_K_) {
            if (linK_) {
                // NOTE: For the time being, there is no expected performance gain from LinK
                // due to a lack of a fast J algorithm to complement LinK
                build_linK(ints, D_ref, K_ref);
                build_JK_matrices(ints, D_ref, J_ref, temp);
            } else {
                build_JK_matrices(ints, D_ref, J_ref, K_ref);
            }
            
        } else if (do_J_) {
            build_JK_matrices(ints, D_ref, J_ref, temp);
        } else {
            if (linK_) {
                build_linK(ints, D_ref, K_ref);
            } else {
                build_JK_matrices(ints, D_ref, temp, K_ref);
            }
        }
    }

    if (incfock_) incfock_postiter();
    
}
void DirectJK::postiterations() {}

void DirectJK::build_JK_matrices(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints, const std::vector<SharedMatrix>& D,
                        std::vector<SharedMatrix>& J, std::vector<SharedMatrix>& K) {

    bool build_J = (!J.empty());
    bool build_K = (!K.empty());

    if (!build_J && !build_K) return;
    
    timer_on("build_JK_matrices()");

    // => Zeroing <= //
    for (auto& Jmat : J) {
        Jmat->zero();
    }
    
    for (auto& Kmat : K) {
        Kmat->zero();
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
            outfile->Printf("  Task: %3d, Task Start: %4d, Task End: %4d\n", task, task_starts[task],
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
                        if (!ints[0]->shell_significant(P, Q, R, S)) continue;

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
        Jmat->scale(2.0);
        Jmat->hermitivitize();
    }

    if (lr_symmetric_) {
        for (auto& Kmat : K) {
            Kmat->scale(2.0);
            Kmat->hermitivitize();
        }
    }

    num_computed_shells_ = computed_shells;
    if (get_bench()) {
        computed_shells_per_iter_.push_back(num_computed_shells());
    }

    timer_off("build_JK_matrices()");
}

// To follow this code, compare with figure 1 of DOI: 10.1063/1.476741
void DirectJK::build_linK(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints, const std::vector<SharedMatrix>& D,
                  std::vector<SharedMatrix>& K) {

    if (!lr_symmetric_) {
        throw PSIEXCEPTION("Non-symmetric K matrix builds are currently not supported in the LinK algorithm.");
    }

    timer_on("build_linK()");

    // ==> Prep Auxiliary Quantities <== //

    // => Zeroing <= //
    for (auto& Kmat : K) {
        Kmat->zero();
    }

    // => Sizing <= //
    int nshell = primary_->nshell();
    int nbf = primary_->nbf();
    int nthread = df_ints_num_threads_;

    // => Atom Blocking <= //
    std::vector<int> shell_endpoints_for_atom;
    std::vector<int> basis_endpoints_for_shell;

    int atomic_ind = -1;
    for (int P = 0; P < nshell; P++) {
        if (primary_->shell(P).ncenter() > atomic_ind) {
            shell_endpoints_for_atom.push_back(P);
            atomic_ind++;
        }
        basis_endpoints_for_shell.push_back(primary_->shell_to_basis_function(P));
    }
    shell_endpoints_for_atom.push_back(nshell);
    basis_endpoints_for_shell.push_back(nbf);

    size_t natom = shell_endpoints_for_atom.size() - 1;

    size_t max_functions_per_atom = 0L;
    for (size_t atom = 0; atom < natom; atom++) {
        size_t size = 0L;
        for (int P = shell_endpoints_for_atom[atom]; P < shell_endpoints_for_atom[atom + 1]; P++) {
            size += primary_->shell(P).nfunction();
        }
        max_functions_per_atom = std::max(max_functions_per_atom, size);
    }

    if (debug_) {
        outfile->Printf("  ==> LinK: Atom Blocking <==\n\n");
        for (size_t atom = 0; atom < natom; atom++) {
            outfile->Printf("  Atom: %3d, Atom Start: %4d, Atom End: %4d\n", atom, shell_endpoints_for_atom[atom],
                            shell_endpoints_for_atom[atom + 1]);
            for (int P = shell_endpoints_for_atom[atom]; P < shell_endpoints_for_atom[atom + 1]; P++) {
                int size = primary_->shell(P).nfunction();
                int off = primary_->shell(P).function_index();
                int off2 = basis_endpoints_for_shell[P];
                outfile->Printf("    Shell: %4d, Size: %4d, Offset: %4d, Offset2: %4d\n", P, size, off,
                                off2);
            }
        }
        outfile->Printf("\n");
    }

    // ==> Prep Atom Pairs <== //
    // Atom-pair blocking inherited from DirectJK code
    // TODO: Test shell-pair blocking

    std::vector<std::pair<int, int>> atom_pairs;
    for (size_t Patom = 0; Patom < natom; Patom++) {
        for (size_t Qatom = 0; Qatom <= Patom; Qatom++) {
            bool found = false;
            for (int P = shell_endpoints_for_atom[Patom]; P < shell_endpoints_for_atom[Patom + 1]; P++) {
                for (int Q = shell_endpoints_for_atom[Qatom]; Q < shell_endpoints_for_atom[Qatom + 1]; Q++) {
                    if (ints[0]->shell_pair_significant(P, Q)) {
                        found = true;
                        atom_pairs.emplace_back(Patom, Qatom);
                        break;
                    }
                }
                if (found) break;
            }
        }
    }

    // ==> Prep Bra-Bra Shell Pairs <== //

    // A comparator used for sorting integral screening values
    auto screen_compare = [](const std::pair<int, double> &a, 
                                    const std::pair<int, double> &b) { return a.second > b.second; };

    std::vector<std::vector<int>> significant_bras(nshell);
    double max_integral = ints[0]->max_integral();

#pragma omp parallel for
    for (size_t P = 0; P < nshell; P++) {
        std::vector<std::pair<int, double>> PQ_shell_values;
        for (size_t Q = 0; Q < nshell; Q++) {
            double pq_pq = std::sqrt(ints[0]->shell_ceiling2(P, Q, P, Q));
            double schwarz_value = std::sqrt(pq_pq * max_integral);
            if (schwarz_value >= cutoff_) {
                PQ_shell_values.emplace_back(Q, schwarz_value);
            }
        }
        std::sort(PQ_shell_values.begin(), PQ_shell_values.end(), screen_compare);

        for (const auto& value : PQ_shell_values) {
            significant_bras[P].push_back(value.first);
        }
    }

    // ==> Prep Bra-Ket Shell Pairs <== //

    // => Calculate Shell Ceilings <= //
    std::vector<double> shell_ceilings(nshell, 0.0);

    // sqrt(Umax|Umax) in Ochsenfeld Eq. 3
#pragma omp parallel for
    for (int P = 0; P < nshell; P++) {
        for (int Q = 0; Q <= P; Q++) {
            double val = std::sqrt(ints[0]->shell_ceiling2(P, Q, P, Q));
            shell_ceilings[P] = std::max(shell_ceilings[P], val);
#pragma omp critical
            shell_ceilings[Q] = std::max(shell_ceilings[Q], val);
        }
    }

    std::vector<std::vector<int>> significant_kets(nshell);

    // => Use shell ceilings to compute significant ket-shells for each bra-shell <= //
#pragma omp parallel for
    for (size_t P = 0; P < nshell; P++) {
        std::vector<std::pair<int, double>> PR_shell_values;
        for (size_t R = 0; R < nshell; R++) {
            double screen_val = shell_ceilings[P] * shell_ceilings[R] * ints[0]->shell_pair_max_density(P, R);
            if (screen_val >= linK_ints_cutoff_) {
                PR_shell_values.emplace_back(R, screen_val);
            }
        }
        std::sort(PR_shell_values.begin(), PR_shell_values.end(), screen_compare);

        for (const auto& value : PR_shell_values) {
            significant_kets[P].push_back(value.first);
        }
    }

    size_t natom_pair = atom_pairs.size();

    // ==> Intermediate Buffers <== //

    // Temporary buffers used during the K contraction process to
    // Take full advantage of permutational symmetry of ERIs
    std::vector<std::vector<SharedMatrix>> KT;

    // To prevent race conditions, give every thread a buffer
    for (int thread = 0; thread < nthread; thread++) {
        std::vector<SharedMatrix> K2;
        for (size_t ind = 0; ind < D.size(); ind++) {
            // (pq|rs) can be contracted into Kpr, Kps, Kqr, Kqs (hence the 4)
            K2.push_back(std::make_shared<Matrix>("KT (linK)", 4 * max_functions_per_atom, nbf));
        }
        KT.push_back(K2);
    }

    // Number of computed shell quartets is tracked for benchmarking purposes
    size_t computed_shells = 0L;

    // ==> Integral Formation Loop <== //

#pragma omp parallel for num_threads(nthread) schedule(dynamic) reduction(+ : computed_shells)
    for (size_t ipair = 0L; ipair < natom_pair; ipair++) { // O(N) shell-pairs in asymptotic limit

        int Patom = atom_pairs[ipair].first;
        int Qatom = atom_pairs[ipair].second;
        
        // Number of shells per atom
        int nPshell = shell_endpoints_for_atom[Patom + 1] - shell_endpoints_for_atom[Patom];
        int nQshell = shell_endpoints_for_atom[Qatom + 1] - shell_endpoints_for_atom[Qatom];

        // First shell per atom
        int Pstart = shell_endpoints_for_atom[Patom];
        int Qstart = shell_endpoints_for_atom[Qatom];

        // Number of basis functions per atom
        int nPbasis = basis_endpoints_for_shell[Pstart + nPshell] - basis_endpoints_for_shell[Pstart];
        int nQbasis = basis_endpoints_for_shell[Qstart + nQshell] - basis_endpoints_for_shell[Qstart];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // Keep track of contraction indices for stripeout (Towards end of this function)
        std::vector<std::unordered_set<int>> P_stripeout_list(nPshell);
        std::vector<std::unordered_set<int>> Q_stripeout_list(nQshell);

        bool touched = false;
        for (int P = Pstart; P < Pstart + nPshell; P++) {
            for (int Q = Qstart; Q < Qstart + nQshell; Q++) {

                if (Q > P) continue;
                if (!ints[0]->shell_pair_significant(P, Q)) continue;

                int dP = P - Pstart;
                int dQ = Q - Qstart;

                // => Formation of Significant Shell Pair List ML <= //

                // Significant ket shell pairs RS for bra shell pair PQ
                // represents the merge of ML_P and ML_Q (mini-lists) as defined in Oschenfeld
                // Unordered set structure allows for automatic merging as new elements are added
                std::unordered_set<int> ML_PQ;

                // Form ML_P as part of ML_PQ
                for (const int R : significant_kets[P]) {
                    bool is_significant = false;
                    for (const int S : significant_bras[R]) {
                        double screen_val = ints[0]->shell_pair_max_density(P, R) * std::sqrt(ints[0]->shell_ceiling2(P, Q, R, S));

                        if (screen_val >= linK_ints_cutoff_) {
                            if (!is_significant) is_significant = true;
                            int RS = (R >= S) ? (R * nshell + S) : (S * nshell + R);
                            if (RS > P * nshell + Q) continue;
                            ML_PQ.emplace(RS);
                            Q_stripeout_list[dQ].emplace(S);
                        }
                        else break;
                    }
                    if (!is_significant) break;
                }

                // Form ML_Q as part of ML_PQ
                for (const int R : significant_kets[Q]) {
                    bool is_significant = false;
                    for (const int S : significant_bras[R]) {
                        double screen_val = ints[0]->shell_pair_max_density(Q, R) * std::sqrt(ints[0]->shell_ceiling2(P, Q, R, S));

                        if (screen_val >= linK_ints_cutoff_) {
                            if (!is_significant) is_significant = true;
                            int RS = (R >= S) ? (R * nshell + S) : (S * nshell + R);
                            if (RS > P * nshell + Q) continue;
                            ML_PQ.emplace(RS);
                            P_stripeout_list[dP].emplace(S);
                        }
                        else break;
                    }
                    if (!is_significant) break;
                }

                // Loop over significant RS pairs
                for (const int RS : ML_PQ) {

                    int R = RS / nshell;
                    int S = RS % nshell;

                    if (!ints[0]->shell_pair_significant(R, S)) continue;
                    if (!ints[0]->shell_significant(P, Q, R, S)) continue;

                    if (ints[thread]->compute_shell(P, Q, R, S) == 0)
                        continue;
                    computed_shells++;

                    const double* buffer = ints[thread]->buffer();

                    // Number of basis functions in shells P, Q, R, S
                    int shell_P_nfunc = primary_->shell(P).nfunction();
                    int shell_Q_nfunc = primary_->shell(Q).nfunction();
                    int shell_R_nfunc = primary_->shell(R).nfunction();
                    int shell_S_nfunc = primary_->shell(S).nfunction();

                    // Basis Function Starting index for shell
                    int shell_P_start = primary_->shell(P).function_index();
                    int shell_Q_start = primary_->shell(Q).function_index();
                    int shell_R_start = primary_->shell(R).function_index();
                    int shell_S_start = primary_->shell(S).function_index();

                    // Basis Function offset from first basis function in the atom
                    int shell_P_offset = basis_endpoints_for_shell[P] - basis_endpoints_for_shell[Pstart];
                    int shell_Q_offset = basis_endpoints_for_shell[Q] - basis_endpoints_for_shell[Qstart];

                    for (size_t ind = 0; ind < D.size(); ind++) {
                        double** Kp = K[ind]->pointer();
                        double** Dp = D[ind]->pointer();
                        double** KTp = KT[thread][ind]->pointer();
                        const double* buffer2 = buffer;

                        if (!touched) {
                            ::memset((void*)KTp[0L * max_functions_per_atom], '\0', nPbasis * nbf * sizeof(double));
                            ::memset((void*)KTp[1L * max_functions_per_atom], '\0', nPbasis * nbf * sizeof(double));
                            ::memset((void*)KTp[2L * max_functions_per_atom], '\0', nQbasis * nbf * sizeof(double));
                            ::memset((void*)KTp[3L * max_functions_per_atom], '\0', nQbasis * nbf * sizeof(double));
                        }

                        // Four pointers needed for PR, PS, QR, QS
                        double* K1p = KTp[0L * max_functions_per_atom];
                        double* K2p = KTp[1L * max_functions_per_atom];
                        double* K3p = KTp[2L * max_functions_per_atom];
                        double* K4p = KTp[3L * max_functions_per_atom];

                        double prefactor = 1.0;
                        if (P == Q) prefactor *= 0.5;
                        if (R == S) prefactor *= 0.5;
                        if (P == R && Q == S) prefactor *= 0.5;

                        // => Computing integral contractions to K buffers <= //
                        for (int p = 0; p < shell_P_nfunc; p++) {
                            for (int q = 0; q < shell_Q_nfunc; q++) {
                                for (int r = 0; r < shell_R_nfunc; r++) {
                                    for (int s = 0; s < shell_S_nfunc; s++) {

                                        K1p[(p + shell_P_offset) * nbf + r + shell_R_start] +=
                                            prefactor * (Dp[q + shell_Q_start][s + shell_S_start]) * (*buffer2);
                                        K2p[(p + shell_P_offset) * nbf + s + shell_S_start] +=
                                            prefactor * (Dp[q + shell_Q_start][r + shell_R_start]) * (*buffer2);
                                        K3p[(q + shell_Q_offset) * nbf + r + shell_R_start] +=
                                            prefactor * (Dp[p + shell_P_start][s + shell_S_start]) * (*buffer2);
                                        K4p[(q + shell_Q_offset) * nbf + s + shell_S_start] +=
                                            prefactor * (Dp[p + shell_P_start][r + shell_R_start]) * (*buffer2);

                                        buffer2++;
                                    }
                                }
                            }
                        }
                    }
                    touched = true;
                }
            }
        }

        // => Master shell quartet loops <= //

        if (!touched) continue;

        // => Stripe out (Writing to K matrix) <= //

        for (size_t ind = 0; ind < D.size(); ind++) {
            double** KTp = KT[thread][ind]->pointer();
            double** Kp = K[ind]->pointer();

            double* K1p = KTp[0L * max_functions_per_atom];
            double* K2p = KTp[1L * max_functions_per_atom];
            double* K3p = KTp[2L * max_functions_per_atom];
            double* K4p = KTp[3L * max_functions_per_atom];

            // K_PR and K_PS
            for (int P = Pstart; P < Pstart + nPshell; P++) {
                int dP = P - Pstart;
                int shell_P_start = primary_->shell(P).function_index();
                int shell_P_nfunc = primary_->shell(P).nfunction();
                int shell_P_offset = basis_endpoints_for_shell[P] - basis_endpoints_for_shell[Pstart];
                for (const int S : P_stripeout_list[dP]) {
                    int shell_S_start = primary_->shell(S).function_index();
                    int shell_S_nfunc = primary_->shell(S).nfunction();

                    for (int p = 0; p < shell_P_nfunc; p++) {
                        for (int s = 0; s < shell_S_nfunc; s++) {
#pragma omp atomic
                            Kp[shell_P_start + p][shell_S_start + s] += K1p[(p + shell_P_offset) * nbf + s + shell_S_start];
#pragma omp atomic
                            Kp[shell_P_start + p][shell_S_start + s] += K2p[(p + shell_P_offset) * nbf + s + shell_S_start];
                        }
                    }

                }
            }

            // K_QR and K_QS
            for (int Q = Qstart; Q < Qstart + nQshell; Q++) {
                int dQ = Q - Qstart;
                int shell_Q_start = primary_->shell(Q).function_index();
                int shell_Q_nfunc = primary_->shell(Q).nfunction();
                int shell_Q_offset = basis_endpoints_for_shell[Q] - basis_endpoints_for_shell[Qstart];
                for (const int S : Q_stripeout_list[dQ]) {
                    int shell_S_start = primary_->shell(S).function_index();
                    int shell_S_nfunc = primary_->shell(S).nfunction();

                    for (int q = 0; q < shell_Q_nfunc; q++) {
                        for (int s = 0; s < shell_S_nfunc; s++) {
#pragma omp atomic
                            Kp[shell_Q_start + q][shell_S_start + s] += K3p[(q + shell_Q_offset) * nbf + s + shell_S_start];
#pragma omp atomic
                            Kp[shell_Q_start + q][shell_S_start + s] += K4p[(q + shell_Q_offset) * nbf + s + shell_S_start];
                        }
                    }

                }
            }

        }  // End stripe out

    }  // End master task list

    for (auto& Kmat : K) {
        Kmat->scale(2.0);
        Kmat->hermitivitize();
    }

    if (bench_) {
        auto mode = std::ostream::app;
        auto printer = PsiOutStream("bench.dat", mode);
        size_t ntri = nshell * (nshell + 1L) / 2L;
        size_t possible_shells = ntri * (ntri + 1L) / 2L;
        printer.Printf("(LinK) Computed %20zu Shell Quartets out of %20zu, (%11.3E ratio)\n", computed_shells,
                        possible_shells, computed_shells / (double)possible_shells);
    }
    timer_off("build_linK()");
}

}  // namespace psi
