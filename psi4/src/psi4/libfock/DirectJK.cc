/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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
#include "psi4/libmints/sieve.h"
#include "psi4/libiwl/iwl.hpp"
#include "jk.h"
//#include "jk_independent.h"
//#include "link.h"
//#include "direct_screening.h"
//#include "cubature.h"
//#include "points.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/lib3index/cholesky.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"

#include <sstream>
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

DirectJK::DirectJK(std::shared_ptr<BasisSet> primary) : JK(primary) { common_init(); }
DirectJK::~DirectJK() {}
void DirectJK::common_init() {
    df_ints_num_threads_ = 1;
#ifdef _OPENMP
    df_ints_num_threads_ = Process::environment.get_n_threads();
#endif
}
size_t DirectJK::memory_estimate() {
    return 0;  // Effectively
}
void DirectJK::print_header() const {
    if (print_) {
        outfile->Printf("  ==> DirectJK: Integral-Direct J/K Matrices <==\n\n");

        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:             %11.3E\n", omega_);
        outfile->Printf("    Integrals threads: %11d\n", df_ints_num_threads_);
        // outfile->Printf( "    Memory [MiB]:      %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf("    Schwarz Cutoff:    %11.0E\n\n", cutoff_);
    }
}
void DirectJK::preiterations() {
    sieve_ = std::make_shared<ERISieve>(primary_, cutoff_, do_csam_);
    
#ifdef USING_BrianQC
    if (brianEnable) {
        double threshold = cutoff_ * (brianCPHFFlag ? 1e-3 : 1e-0); // CPHF needs higher precision
        brianCOMSetPrecisionThresholds(&brianCookie, &threshold);
        checkBrian();
    }
#endif
}
void DirectJK::compute_JK() {
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
                
                std::vector<std::vector<std::shared_ptr<Matrix>>> pseudoDensitySymmetrized(densityCount);
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

    auto factory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);

    if (do_wK_) {
        std::vector<std::shared_ptr<TwoBodyAOInt> > ints;
        for (int thread = 0; thread < df_ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega_)));
        }
        // TODO: Fast K algorithm
        if (do_J_) {
            build_JK(ints, D_ao_, J_ao_, wK_ao_);
        } else {
            std::vector<std::shared_ptr<Matrix> > temp;
            for (size_t i = 0; i < D_ao_.size(); i++) {
                temp.push_back(std::make_shared<Matrix>("temp", primary_->nbf(), primary_->nbf()));
            }
            build_JK(ints, D_ao_, temp, wK_ao_);
        }
    }

    if (do_J_ || do_K_) {
        std::vector<std::shared_ptr<TwoBodyAOInt> > ints;
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri()));
        for (int thread = 1; thread < df_ints_num_threads_; thread++) {
            if (ints[0]->cloneable())
                ints.push_back(std::shared_ptr<TwoBodyAOInt>(ints[0]->clone()));
            else
                ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri()));
        }
        if (do_J_ && do_K_) {
            build_JK(ints, D_ao_, J_ao_, K_ao_);
        } else if (do_J_) {
            std::vector<std::shared_ptr<Matrix> > temp;
            for (size_t i = 0; i < D_ao_.size(); i++) {
                temp.push_back(std::make_shared<Matrix>("temp", primary_->nbf(), primary_->nbf()));
            }
            build_JK(ints, D_ao_, J_ao_, temp);
        } else {
            std::vector<std::shared_ptr<Matrix> > temp;
            for (size_t i = 0; i < D_ao_.size(); i++) {
                temp.push_back(std::make_shared<Matrix>("temp", primary_->nbf(), primary_->nbf()));
            }
            build_JK(ints, D_ao_, temp, K_ao_);
        }
    }
}
void DirectJK::postiterations() { sieve_.reset(); }
void DirectJK::build_JK(std::vector<std::shared_ptr<TwoBodyAOInt> >& ints, std::vector<std::shared_ptr<Matrix> >& D,
                        std::vector<std::shared_ptr<Matrix> >& J, std::vector<std::shared_ptr<Matrix> >& K) {

    // => Zeroing <= //
    for (size_t ind = 0; ind < J.size(); ind++) {
        J[ind]->zero();
    }
    for (size_t ind = 0; ind < K.size(); ind++) {
        K[ind]->zero();
    }

    int nthread = df_ints_num_threads_;

    /*
     * This version of the code processes a single shell block at a time.  The previous version
     * computed a batch comprising all shells on each atom - we may want to consider using larger
     * blocks in this version after profiling.  For now, this is the simplest implementation.
     */
    size_t max_task = primary_->max_function_per_shell();
    // => Intermediate Buffers <= //

    std::vector<std::vector<std::shared_ptr<Matrix> > > JKT;
    for (int thread = 0; thread < nthread; thread++) {
        std::vector<std::shared_ptr<Matrix> > JK2;
        for (size_t ind = 0; ind < D.size(); ind++) {
            JK2.push_back(std::make_shared<Matrix>("JKT", (lr_symmetric_ ? 6 : 10) * max_task, max_task));
        }
        JKT.push_back(JK2);
    }

    size_t computed_shells = 0L;
    // shell pair blocks
    auto blocksPQ = ints[0]->get_blocks12();
    auto blocksRS = ints[0]->get_blocks34();
    bool use_batching = blocksPQ != blocksRS;

#pragma omp parallel for schedule(dynamic) num_threads(nthread)
    // loop over all the blocks of (P>=Q|
    for (size_t blockPQ_idx = 0; blockPQ_idx < blocksPQ.size(); blockPQ_idx++) {
        const auto& blockPQ = blocksPQ[blockPQ_idx];
#ifdef _OPENMP
        const int rank = omp_get_thread_num();
#else
        const int rank = 0;
#endif
        // loop over all the blocks of |R>=S)
        int loop_start = use_batching ? 0 : blockPQ_idx;
        for (int blockRS_idx = loop_start; blockRS_idx < blocksRS.size(); ++blockRS_idx) {
            const auto& blockRS = blocksRS[blockRS_idx];

            // This is where we want to screen with density and schwarz-like screening

            // compute the integrals and continue if none were computed
            ints[rank]->compute_shell_blocks(blockPQ_idx, blockRS_idx);
            const double* block_start = ints[rank]->buffer();

            // Loop over all of the P,Q,R,S shells within the blocks.  We have P>=Q, R>=S and PQ<=RS.
            for (const auto& pairPQ : blockPQ) {
                const auto &P = pairPQ.first;
                const auto &Q = pairPQ.second;
                const auto& Pshell = primary_->shell(P);
                const auto& Qshell = primary_->shell(Q);
                const auto& Psize = Pshell.nfunction();
                const auto& Qsize = Qshell.nfunction();
                const auto& Poff = Pshell.function_index();
                const auto& Qoff = Qshell.function_index();

                for (const auto& pairRS : blockRS) {
                    const auto &R = pairRS.first;
                    const auto &S = pairRS.second;
                    const auto & Rshell = primary_->shell(R);
                    const auto & Sshell = primary_->shell(S);
                    int Rsize = Rshell.nfunction();
                    int Ssize = Sshell.nfunction();
                    int Roff = Rshell.function_index();
                    int Soff = Sshell.function_index();

                    size_t block_size = Psize * Qsize * Rsize * Ssize;
                    // When there are chunks of shellpairs in RS, we need to make sure
                    // we filter out redundant combinations.  This should probably be done
                    // by having a block of RS generated for each PQ at list build time.
                    if (use_batching && ((P > R) || (P == R && Q > S))) {
                        block_start += block_size;
                        continue;
                    }

                    // if (thread == 0) timer_on("JK: GEMV");
                    for (size_t ind = 0; ind < D.size(); ind++) {
                        double** Dp = D[ind]->pointer();
                        double** JKTp = JKT[rank][ind]->pointer();

                        double* J1p = JKTp[0L * max_task];
                        double* J2p = JKTp[1L * max_task];
                        double* K1p = JKTp[2L * max_task];
                        double* K2p = JKTp[3L * max_task];
                        double* K3p = JKTp[4L * max_task];
                        double* K4p = JKTp[5L * max_task];
                        double* K5p;
                        double* K6p;
                        double* K7p;
                        double* K8p;
                        if (!lr_symmetric_) {
                            K5p = JKTp[6L * max_task];
                            K6p = JKTp[7L * max_task];
                            K7p = JKTp[8L * max_task];
                            K8p = JKTp[9L * max_task];
                        }

                        double prefactor = 1.0;
                        if (P == Q) prefactor *= 0.5;
                        if (R == S) prefactor *= 0.5;
                        if (P == R && Q == S) prefactor *= 0.5;

                        std::fill_n(JKTp[0L * max_task], Psize * Qsize, 0.0);
                        std::fill_n(JKTp[1L * max_task], Rsize * Ssize, 0.0);
                        std::fill_n(JKTp[2L * max_task], Psize * Rsize, 0.0);
                        std::fill_n(JKTp[3L * max_task], Psize * Ssize, 0.0);
                        std::fill_n(JKTp[4L * max_task], Qsize * Rsize, 0.0);
                        std::fill_n(JKTp[5L * max_task], Qsize * Ssize, 0.0);
                        if (!lr_symmetric_) {
                            std::fill_n(JKTp[6L * max_task], Rsize * Psize, 0.0);
                            std::fill_n(JKTp[7L * max_task], Ssize * Psize, 0.0);
                            std::fill_n(JKTp[8L * max_task], Rsize * Qsize, 0.0);
                            std::fill_n(JKTp[9L * max_task], Ssize * Qsize, 0.0);
                        }
                        const double* int_ptr = block_start;
                        for (int p = 0; p < Psize; p++) {
                            for (int q = 0; q < Qsize; q++) {
                                for (int r = 0; r < Rsize; r++) {
                                    for (int s = 0; s < Ssize; s++) {
                                        double val = prefactor * (*int_ptr++);
                                        J1p[p * Qsize + q] += val * (Dp[r + Roff][s + Soff] + Dp[s + Soff][r + Roff]);
                                        J2p[r * Ssize + s] += val * (Dp[p + Poff][q + Qoff] + Dp[q + Qoff][p + Poff]);
                                        K1p[p * Rsize + r] += val * Dp[q + Qoff][s + Soff];
                                        K2p[p * Ssize + s] += val * Dp[q + Qoff][r + Roff];
                                        K3p[q * Rsize + r] += val * Dp[p + Poff][s + Soff];
                                        K4p[q * Ssize + s] += val * Dp[p + Poff][r + Roff];
                                        if (!lr_symmetric_) {
                                            K5p[r * Psize + p] += val * Dp[s + Soff][q + Qoff];
                                            K6p[s * Psize + p] += val * Dp[r + Roff][q + Qoff];
                                            K7p[r * Qsize + q] += val * Dp[s + Soff][p + Poff];
                                            K8p[s * Qsize + q] += val * Dp[r + Roff][p + Poff];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    block_start += block_size;

                    // => Stripe out <= //

                    // if (thread == 0) timer_on("JK: Atomic");
                    for (size_t ind = 0; ind < D.size(); ind++) {
                        double** JKTp = JKT[rank][ind]->pointer();
                        double** Jp = J[ind]->pointer();
                        double** Kp = K[ind]->pointer();

                        double* J1p = JKTp[0L * max_task];
                        double* J2p = JKTp[1L * max_task];
                        double* K1p = JKTp[2L * max_task];
                        double* K2p = JKTp[3L * max_task];
                        double* K3p = JKTp[4L * max_task];
                        double* K4p = JKTp[5L * max_task];
                        double* K5p;
                        double* K6p;
                        double* K7p;
                        double* K8p;
                        if (!lr_symmetric_) {
                            K5p = JKTp[6L * max_task];
                            K6p = JKTp[7L * max_task];
                            K7p = JKTp[8L * max_task];
                            K8p = JKTp[9L * max_task];
                        }

                        // > J_PQ < //
                        for (int p = 0; p < Psize; p++) {
                            for (int q = 0; q < Qsize; q++) {
#pragma omp atomic
                                Jp[p + Poff][q + Qoff] += J1p[p * Qsize + q];
                            }
                        }

                        // > J_RS < //

                        for (int r = 0; r < Rsize; r++) {
                            for (int s = 0; s < Ssize; s++) {
#pragma omp atomic
                                Jp[r + Roff][s + Soff] += J2p[r * Ssize + s];
                            }
                        }

                        // > K_PR < //

                        for (int p = 0; p < Psize; p++) {
                            for (int r = 0; r < Rsize; r++) {
#pragma omp atomic
                                Kp[p + Poff][r + Roff] += K1p[p * Rsize + r];
                                if (!lr_symmetric_) {
#pragma omp atomic
                                    Kp[r + Roff][p + Poff] += K5p[r * Psize + p];
                                }
                            }
                        }

                        // > K_PS < //

                        for (int p = 0; p < Psize; p++) {
                            for (int s = 0; s < Ssize; s++) {
#pragma omp atomic
                                Kp[p + Poff][s + Soff] += K2p[p * Ssize + s];
                                if (!lr_symmetric_) {
#pragma omp atomic
                                    Kp[s + Soff][p + Poff] += K6p[s * Psize + p];
                                }
                            }
                        }

                        // > K_QR < //

                        for (int q = 0; q < Qsize; q++) {
                            for (int r = 0; r < Rsize; r++) {
#pragma omp atomic
                                Kp[q + Qoff][r + Roff] += K3p[q * Rsize + r];
                                if (!lr_symmetric_) {
#pragma omp atomic
                                    Kp[r + Roff][q + Qoff] += K7p[r * Qsize + q];
                                }
                            }
                        }
                        // > K_QS < //

                        for (int q = 0; q < Qsize; q++) {
                            for (int s = 0; s < Ssize; s++) {
#pragma omp atomic
                                Kp[q + Qoff][s + Soff] += K4p[q * Ssize + s];
                                if (!lr_symmetric_) {
#pragma omp atomic
                                    Kp[s + Soff][q + Qoff] += K8p[s * Qsize + q];
                                }
                            }
                        }
                    }
                }  // pairRS
            }      // pairPQ
        }          // blockRS
    }              // blockPQ

    for (size_t ind = 0; ind < D.size(); ind++) {
        J[ind]->scale(2.0);
        J[ind]->hermitivitize();
        if (lr_symmetric_) {
            K[ind]->scale(2.0);
            K[ind]->hermitivitize();
        }
    }

    if (bench_) {
        auto mode = std::ostream::app;
        auto printer = std::make_shared<PsiOutStream>("bench.dat", mode);
        size_t ntri = nshell * (nshell + 1L) / 2L;
        size_t possible_shells = ntri * (ntri + 1L) / 2L;
        printer->Printf("Computed %20zu Shell Quartets out of %20zu, (%11.3E ratio)\n", computed_shells,
                        possible_shells, computed_shells / (double)possible_shells);
    }
    
//     for (size_t ind = 0; ind < J.size(); ind++) {
//         outfile->Printf("J[%d] psi:\n", ind);
//         J[ind]->print();
//         outfile->Printf("K[%d] psi:\n", ind);
//         K[ind]->print();
//     }
}

#if 0
DirectJK::DirectJK(std::shared_ptr<BasisSet> primary) :
   JK(primary)
{
    common_init();
}
DirectJK::~DirectJK()
{
}
void DirectJK::common_init()
{
}
void DirectJK::print_header() const
{
    if (print_) {
        outfile->Printf( "  ==> DirectJK: Integral-Direct J/K Matrices <==\n\n");

        outfile->Printf( "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf( "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf( "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_)
            outfile->Printf( "    Omega:             %11.3E\n", omega_);
        outfile->Printf( "    OpenMP threads:    %11d\n", omp_nthread_);
        outfile->Printf( "    Memory [MiB]:      %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf( "    Schwarz Cutoff:    %11.0E\n\n", cutoff_);
    }
}
void DirectJK::preiterations()
{
    sieve_ = std::make_shared<ERISieve>(primary_, cutoff_);
    factory_= std::make_shared<IntegralFactory>(primary_,primary_,primary_,primary_);
    eri_.clear();
    for (int thread = 0; thread < omp_nthread_; thread++) {
        eri_.push_back(std::shared_ptr<TwoBodyAOInt>(factory_->erd_eri()));
    }
}
void DirectJK::compute_JK()
{
    // Correctness always counts
    const double* buffer = eri_[0]->buffer();
    for (int M = 0; M < primary_->nshell(); ++M) {
    for (int N = 0; N < primary_->nshell(); ++N) {
    for (int R = 0; R < primary_->nshell(); ++R) {
    for (int S = 0; S < primary_->nshell(); ++S) {

        if(eri_[0]->compute_shell(M,N,R,S) == 0)
            continue; // No integrals were computed here

        int nM = primary_->shell(M).nfunction();
        int nN = primary_->shell(N).nfunction();
        int nR = primary_->shell(R).nfunction();
        int nS = primary_->shell(S).nfunction();

        int sM = primary_->shell(M).function_index();
        int sN = primary_->shell(N).function_index();
        int sR = primary_->shell(R).function_index();
        int sS = primary_->shell(S).function_index();

        for (int oM = 0, index = 0; oM < nM; oM++) {
        for (int oN = 0; oN < nN; oN++) {
        for (int oR = 0; oR < nR; oR++) {
        for (int oS = 0; oS < nS; oS++, index++) {

            double val = buffer[index];

            int m = oM + sM;
            int n = oN + sN;
            int r = oR + sR;
            int s = oS + sS;

            if (do_J_) {
                for (int N = 0; N < J_ao_.size(); N++) {
                    J_ao_[N]->add(0,m,n, D_ao_[N]->get(0,r,s)*val);
                }
            }

            if (do_K_) {
                for (int N = 0; N < K_ao_.size(); N++) {
                    K_ao_[N]->add(0,m,s, D_ao_[N]->get(0,n,r)*val);
                }
            }

        }}}}

    }}}}

    // Faster version, not finished
    /**
    sieve_->set_sieve(cutoff_);
    const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();
    size_t nMN = shell_pairs.size();
    size_t nMNRS = nMN * nMN;
    int nthread = eri_.size();

    #pragma omp parallel for schedule(dynamic,30) num_threads(nthread)
    for (size_t index = 0L; index < nMNRS; ++index) {

        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif

        const double* buffer = eri_[thread]->buffer();

        size_t MN = index / nMN;
        size_t RS = index % nMN;
        if (MN < RS) continue;

        int M = shell_pairs[MN].first;
        int N = shell_pairs[MN].second;
        int R = shell_pairs[RS].first;
        int S = shell_pairs[RS].second;

        eri_[thread]->compute_shell(M,N,R,S);

        int nM = primary_->shell(M)->nfunction();
        int nN = primary_->shell(N)->nfunction();
        int nR = primary_->shell(R)->nfunction();
        int nS = primary_->shell(S)->nfunction();

        int sM = primary_->shell(M)->function_index();
        int sN = primary_->shell(N)->function_index();
        int sR = primary_->shell(R)->function_index();
        int sS = primary_->shell(S)->function_index();

        for (int oM = 0, index = 0; oM < nM; oM++) {
        for (int oN = 0; oN < nN; oN++) {
        for (int oR = 0; oR < nR; oR++) {
        for (int oS = 0; oS < nS; oS++, index++) {

            int m = oM + sM;
            int n = oN + sN;
            int r = oR + sR;
            int s = oS + sS;

            if ((n > m) || (s > r) || ((r*(r+1) >> 1) + s > (m*(m+1) >> 1) + n)) continue;

            double val = buffer[index];

            if (do_J_) {
                for (int N = 0; N < J_ao_.size(); N++) {
                    double** Dp = D_ao_[N]->pointer();
                    double** Jp = J_ao_[N]->pointer();

                    // I've given you all the unique ones
                    // Make sure to use #pragma omp atomic
                    // TODO
                }
            }

            if (do_K_) {
                for (int N = 0; N < K_ao_.size(); N++) {
                    double** Dp = D_ao_[N]->pointer();
                    double** Kp = J_ao_[N]->pointer();

                    // I've given you all the unique ones
                    // Make sure to use #pragma omp atomic
                    // TODO
                }
            }

        }}}}

    }
    **/
}

#endif
}  // namespace psi
