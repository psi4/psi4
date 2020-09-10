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
#include <set>
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
        std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
        for (int thread = 0; thread < df_ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega_)));
        }
        // TODO: Fast K algorithm
        if (do_J_) {
            build_JK(ints, D_ao_, J_ao_, wK_ao_);
        } else {
            std::vector<std::shared_ptr<Matrix>> temp;
            for (size_t i = 0; i < D_ao_.size(); i++) {
                temp.push_back(std::make_shared<Matrix>("temp", primary_->nbf(), primary_->nbf()));
            }
            build_JK(ints, D_ao_, temp, wK_ao_);
        }
    }

    if (do_J_ || do_K_) {
        std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri()));
        for (int thread = 1; thread < df_ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(ints[0]->clone()));
        }
        if (do_J_ && do_K_) {
            build_JK(ints, D_ao_, J_ao_, K_ao_);
        } else if (do_J_) {
            std::vector<std::shared_ptr<Matrix>> temp;
            for (size_t i = 0; i < D_ao_.size(); i++) {
                temp.push_back(std::make_shared<Matrix>("temp", primary_->nbf(), primary_->nbf()));
            }
            build_JK(ints, D_ao_, J_ao_, temp);
        } else {
            std::vector<std::shared_ptr<Matrix>> temp;
            for (size_t i = 0; i < D_ao_.size(); i++) {
                temp.push_back(std::make_shared<Matrix>("temp", primary_->nbf(), primary_->nbf()));
            }
            build_JK(ints, D_ao_, temp, K_ao_);
        }
    }
}
void DirectJK::postiterations() {}

void DirectJK::build_JK(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints, std::vector<std::shared_ptr<Matrix>>& D,
                        std::vector<std::shared_ptr<Matrix>>& J, std::vector<std::shared_ptr<Matrix>>& K) {
    // => Zeroing <= //
    for (size_t ind = 0; ind < J.size(); ind++) {
        J[ind]->zero();
    }
    for (size_t ind = 0; ind < K.size(); ind++) {
        K[ind]->zero();
    }

    int nthread = df_ints_num_threads_;

    // Get function pair blocking info from the integral engine. We don't use batching in the ket here
    // because the experimental code that implemented it showed a speedup in the generation of the ints
    // but a fairly dramatic slowdown in the speed of their consumption in building the Fock matrix.
    auto blocksPQ = ints[0]->get_blocks12();
    auto blocksRS = ints[0]->get_blocks12();

    // => Find the groups of PQ/RS blocks that have the same atoms involved <= //
    std::vector<std::array<int, 10>> PQtask_offsets;
    int current_atom_p = -1;
    int current_atom_q = -1;
    int dim_p = 0;
    int dim_q = 0;
    int max_p = 0;
    int max_q = 0;
    int min_p = std::numeric_limits<int>::max();
    int min_q = std::numeric_limits<int>::max();
    int current_offset = 0;
    int atom_p;
    int atom_q;
    for (size_t blockPQ_idx = 0; blockPQ_idx < blocksPQ.size(); blockPQ_idx++) {
        int shell_p = blocksPQ[blockPQ_idx][0].first;
        int shell_q = blocksPQ[blockPQ_idx][0].second;
        atom_p = primary_->shell(shell_p).ncenter();
        atom_q = primary_->shell(shell_q).ncenter();
        if (current_atom_p != atom_p || current_atom_q != atom_q) {
            if (blockPQ_idx != 0) {
                dim_p = 0;
                dim_q = 0;
                for (int pshell = min_p; pshell <= max_p; ++pshell) dim_p += primary_->shell(pshell).nfunction();
                for (int qshell = min_q; qshell <= max_q; ++qshell) dim_q += primary_->shell(qshell).nfunction();
                PQtask_offsets.emplace_back(std::array<int, 10>{
                    {atom_p, atom_q, current_offset, static_cast<int>(blockPQ_idx), min_p, min_q, max_p, max_q, dim_p, dim_q}});
            }
            current_offset = blockPQ_idx;
            current_atom_p = atom_p;
            current_atom_q = atom_q;
            max_p = 0;
            max_q = 0;
            min_p = std::numeric_limits<int>::max();
            min_q = std::numeric_limits<int>::max();
        }
        min_p = std::min(min_p, shell_p);
        min_q = std::min(min_q, shell_q);
        max_p = std::max(max_p, shell_p);
        max_q = std::max(max_q, shell_q);
    }
    dim_p = 0;
    dim_q = 0;
    for (int pshell = min_p; pshell <= max_p; ++pshell) dim_p += primary_->shell(pshell).nfunction();
    for (int qshell = min_q; qshell <= max_q; ++qshell) dim_q += primary_->shell(qshell).nfunction();
    PQtask_offsets.emplace_back(std::array<int, 10>{
        {atom_p, atom_q, current_offset, static_cast<int>(blocksPQ.size()), min_p, min_q, max_p, max_q, dim_p, dim_q}});
    std::vector<std::array<int, 10>> RStask_offsets = PQtask_offsets;

    // => The maximum amount of storage associated with any atom =< //
    std::vector<int> nfunctions(primary_->molecule()->natom(), 0);
    for (int shellnum = 0; shellnum < primary_->nshell(); ++shellnum) {
        const auto& shell = primary_->shell(shellnum);
        nfunctions[shell.ncenter()] += shell.nfunction();
    }
    size_t max_task = *std::max_element(nfunctions.begin(), nfunctions.end());

    // => Intermediate Buffers <= //

    std::vector<std::vector<std::shared_ptr<Matrix>>> JKT;
    for (int thread = 0; thread < nthread; thread++) {
        std::vector<std::shared_ptr<Matrix>> JK2;
        for (size_t ind = 0; ind < D.size(); ind++) {
            JK2.push_back(std::make_shared<Matrix>("JKT", (lr_symmetric_ ? 6 : 10) * max_task, max_task));
        }
        JKT.push_back(JK2);
    }

    size_t computed_shells = 0L;
    size_t numPQtasks = PQtask_offsets.size();
    size_t numRStasks = RStask_offsets.size();
    size_t numPQRStasks = numPQtasks * numRStasks;
#pragma omp parallel for schedule(dynamic) num_threads(nthread)  reduction(+ : computed_shells)
    for (size_t PQRStask = 0; PQRStask < numPQRStasks; ++PQRStask) {
#ifdef _OPENMP
        const int rank = omp_get_thread_num();
#else
        const int rank = 0;
#endif
        size_t PQtask = PQRStask / numRStasks;
        size_t RStask = PQRStask % numRStasks;
        const auto& PQBlockInfo = PQtask_offsets[PQtask];
        const auto& RSBlockInfo = RStask_offsets[RStask];
        const auto& atomP = PQBlockInfo[0];
        const auto& atomQ = PQBlockInfo[1];
        const auto& atomR = RSBlockInfo[0];
        const auto& atomS = RSBlockInfo[1];
        if (atomR > atomP) continue;
        const auto& buffers = ints[rank]->buffers();
        bool touched = false;
        const auto& firstPQ = PQBlockInfo[2];
        const auto& firstRS = RSBlockInfo[2];
        const auto& lastPQ = PQBlockInfo[3];
        const auto& lastRS = RSBlockInfo[3];
        const auto& firstP = PQBlockInfo[4];
        const auto& firstQ = PQBlockInfo[5];
        const auto& firstR = RSBlockInfo[4];
        const auto& firstS = RSBlockInfo[5];
        const auto& lastP = PQBlockInfo[6];
        const auto& lastQ = PQBlockInfo[7];
        const auto& lastR = RSBlockInfo[6];
        const auto& lastS = RSBlockInfo[7];
        const auto& totalPsize = PQBlockInfo[8];
        const auto& totalQsize = PQBlockInfo[9];
        const auto& totalRsize = RSBlockInfo[8];
        const auto& totalSsize = RSBlockInfo[9];
        const auto& firstPshell = primary_->shell(firstP);
        const auto& firstQshell = primary_->shell(firstQ);
        const auto& firstRshell = primary_->shell(firstR);
        const auto& firstSshell = primary_->shell(firstS);
        const auto firstPoff = firstPshell.function_index();
        const auto firstQoff = firstQshell.function_index();
        const auto firstRoff = firstRshell.function_index();
        const auto firstSoff = firstSshell.function_index();

        // loop over all the blocks of (P>=Q| belonging to the bra atom pair
        for (size_t blockPQ_idx = firstPQ; blockPQ_idx < lastPQ; blockPQ_idx++) {
            const auto& blockPQ = blocksPQ[blockPQ_idx];
            // Loop over all of the P,Q,R,S shells within the blocks.  We have P>=Q, R>=S and PQ<=RS.
            const auto& pairPQ = blockPQ[0];
            int P = pairPQ.first;
            int Q = pairPQ.second;
            const auto& Pshell = primary_->shell(P);
            const auto& Qshell = primary_->shell(Q);
            const auto& Pam = Pshell.am();
            const auto& Qam = Qshell.am();
            const auto& Psize = Pshell.nfunction();
            const auto& Qsize = Qshell.nfunction();
            const auto& Poff = Pshell.function_index();
            const auto& Qoff = Qshell.function_index();
            const int relPoff = Poff - firstPoff;
            const int relQoff = Qoff - firstQoff;
            // loop over all the blocks of |R>=S) belonging to the ket atom pair
            for (int blockRS_idx = firstRS; blockRS_idx < lastRS; ++blockRS_idx) {
                if (blockRS_idx > blockPQ_idx) continue;
                const auto& blockRS = blocksRS[blockRS_idx];

                const auto& pairRS = blockRS[0];
                int R = pairRS.first;
                int S = pairRS.second;
                const auto& Rshell = primary_->shell(R);
                const auto& Sshell = primary_->shell(S);
                const auto& Ram = Rshell.am();
                const auto& Sam = Sshell.am();
                const auto& Rsize = Rshell.nfunction();
                const auto& Ssize = Sshell.nfunction();
                const auto& Roff = Rshell.function_index();
                const auto& Soff = Sshell.function_index();
                const int relRoff = Roff - firstRoff;
                const int relSoff = Soff - firstSoff;

                size_t block_size = Psize * Qsize * Rsize * Ssize;
                // compute the integrals and continue if none were computed
                // if (rank == 0) timer_off("JK: INTS");
                if (!ints[rank]->shell_significant(P, Q, R, S)) continue;
                if (ints[rank]->compute_shell(P, Q, R, S) == 0) continue;
                computed_shells++;

                double prefactor = 1.0;
                if (P == Q) prefactor *= 0.5;
                if (R == S) prefactor *= 0.5;
                if (P == R && Q == S) prefactor *= 0.5;

                // if (rank == 0) timer_on("JK: GEMV");
                for (size_t ind = 0; ind < D.size(); ind++) {
                    const double* int_ptr = ints[rank]->buffers()[0];
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

                    if (!touched) {
                        std::fill_n(JKTp[0L * max_task], totalPsize * totalQsize, 0.0);
                        std::fill_n(JKTp[1L * max_task], totalRsize * totalSsize, 0.0);
                        std::fill_n(JKTp[2L * max_task], totalPsize * totalRsize, 0.0);
                        std::fill_n(JKTp[3L * max_task], totalPsize * totalSsize, 0.0);
                        std::fill_n(JKTp[4L * max_task], totalQsize * totalRsize, 0.0);
                        std::fill_n(JKTp[5L * max_task], totalQsize * totalSsize, 0.0);
                        if (!lr_symmetric_) {
                            std::fill_n(JKTp[6L * max_task], totalRsize * totalPsize, 0.0);
                            std::fill_n(JKTp[7L * max_task], totalSsize * totalPsize, 0.0);
                            std::fill_n(JKTp[8L * max_task], totalRsize * totalQsize, 0.0);
                            std::fill_n(JKTp[9L * max_task], totalSsize * totalQsize, 0.0);
                        }
                    }
                    for (int p = 0; p < Psize; p++) {
                        const int Prel = p + relPoff;
                        const int Pabs = p + Poff;
                        for (int q = 0; q < Qsize; q++) {
                            const int Qrel = q + relQoff;
                            const int Qabs = q + Qoff;
                            for (int r = 0; r < Rsize; r++) {
                                const int Rrel = r + relRoff;
                                const int Rabs = r + Roff;
                                for (int s = 0; s < Ssize; s++) {
                                    const int Srel = s + relSoff;
                                    const int Sabs = s + Soff;
                                    double val = prefactor * (*int_ptr++);
                                    J1p[Prel * totalQsize + Qrel] += val * (Dp[Rabs][Sabs] + Dp[Sabs][Rabs]);
                                    J2p[Rrel * totalSsize + Srel] += val * (Dp[Pabs][Qabs] + Dp[Qabs][Pabs]);
                                    K1p[Prel * totalRsize + Rrel] += val * Dp[Qabs][Sabs];
                                    K2p[Prel * totalSsize + Srel] += val * Dp[Qabs][Rabs];
                                    K3p[Qrel * totalRsize + Rrel] += val * Dp[Pabs][Sabs];
                                    K4p[Qrel * totalSsize + Srel] += val * Dp[Pabs][Rabs];
                                    if (!lr_symmetric_) {
                                        K5p[Rrel * totalPsize + Prel] += val * Dp[Sabs][Qabs];
                                        K6p[Srel * totalPsize + Prel] += val * Dp[Rabs][Qabs];
                                        K7p[Rrel * totalQsize + Qrel] += val * Dp[Sabs][Pabs];
                                        K8p[Srel * totalQsize + Qrel] += val * Dp[Rabs][Pabs];
                                    }
                                }
                            }
                        }
                    }
                }
                touched = true;
                // if (rank == 0) timer_off("JK: GEMV");
            }  // pairRS
        }      // pairPQ

        // => Accumulate atom blocks into full J and K matrices <= //
        if (!touched) continue;
        // if (rank == 0) timer_on("JK: Atomic");
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
            for (int P = firstP; P <= lastP; ++P) {
                for (int Q = firstQ; Q <= lastQ; ++Q) {
                    const auto& Pshell = primary_->shell(P);
                    const int Psize = Pshell.nfunction();
                    const int Poff = Pshell.function_index();
                    const int relPoff = Poff - firstPoff;
                    const auto& Qshell = primary_->shell(Q);
                    const int Qsize = Qshell.nfunction();
                    const int Qoff = Qshell.function_index();
                    const int relQoff = Qoff - firstQoff;
                    for (int p = 0; p < Psize; p++) {
                        const int pRel = p + relPoff;
                        const int pAbs = p + Poff;
                        for (int q = 0; q < Qsize; q++) {
                            const int qRel = q + relQoff;
                            const int qAbs = q + Qoff;
#pragma omp atomic
                            Jp[pAbs][qAbs] += J1p[pRel * totalQsize + qRel];
                        }
                    }
                }
            }

            // > J_RS < //

            for (int R = firstR; R <= lastR; ++R) {
                for (int S = firstS; S <= lastS; ++S) {
                    const auto& Rshell = primary_->shell(R);
                    const int Rsize = Rshell.nfunction();
                    const int Roff = Rshell.function_index();
                    const int relRoff = Roff - firstRoff;
                    const auto& Sshell = primary_->shell(S);
                    const int Ssize = Sshell.nfunction();
                    const int Soff = Sshell.function_index();
                    const int relSoff = Soff - firstSoff;
                    for (int r = 0; r < Rsize; r++) {
                        const int rRel = r + relRoff;
                        const int rAbs = r + Roff;
                        for (int s = 0; s < Ssize; s++) {
                            const int sRel = s + relSoff;
                            const int sAbs = s + Soff;
#pragma omp atomic
                            Jp[rAbs][sAbs] += J2p[rRel * totalSsize + sRel];
                        }
                    }
                }
            }

            // > K_PR < //

            for (int P = firstP; P <= lastP; ++P) {
                for (int R = firstR; R <= lastR; ++R) {
                    const auto& Pshell = primary_->shell(P);
                    const int Psize = Pshell.nfunction();
                    const int Poff = Pshell.function_index();
                    const int relPoff = Poff - firstPoff;
                    const auto& Rshell = primary_->shell(R);
                    const int Rsize = Rshell.nfunction();
                    const int Roff = Rshell.function_index();
                    const int relRoff = Roff - firstRoff;
                    for (int p = 0; p < Psize; p++) {
                        const int pRel = p + relPoff;
                        const int pAbs = p + Poff;
                        for (int r = 0; r < Rsize; r++) {
                            const int rRel = r + relRoff;
                            const int rAbs = r + Roff;
#pragma omp atomic
                            Kp[pAbs][rAbs] += K1p[pRel * totalRsize + rRel];
                            if (!lr_symmetric_) {
#pragma omp atomic
                                Kp[rAbs][pAbs] += K5p[rRel * totalPsize + pRel];
                            }
                        }
                    }
                }
            }

            // > K_PS < //

            for (int P = firstP; P <= lastP; ++P) {
                for (int S = firstS; S <= lastS; ++S) {
                    const auto& Pshell = primary_->shell(P);
                    const int Psize = Pshell.nfunction();
                    const int Poff = Pshell.function_index();
                    const int relPoff = Poff - firstPoff;
                    const auto& Sshell = primary_->shell(S);
                    const int Ssize = Sshell.nfunction();
                    const int Soff = Sshell.function_index();
                    const int relSoff = Soff - firstSoff;
                    for (int p = 0; p < Psize; p++) {
                        const int pRel = p + relPoff;
                        const int pAbs = p + Poff;
                        for (int s = 0; s < Ssize; s++) {
                            const int sRel = s + relSoff;
                            const int sAbs = s + Soff;
#pragma omp atomic
                            Kp[pAbs][sAbs] += K2p[pRel * totalSsize + sRel];
                            if (!lr_symmetric_) {
#pragma omp atomic
                                Kp[sAbs][pAbs] += K6p[sRel * totalPsize + pRel];
                            }
                        }
                    }
                }
            }

            // > K_QR < //

            for (int Q = firstQ; Q <= lastQ; ++Q) {
                for (int R = firstR; R <= lastR; ++R) {
                    const auto& Qshell = primary_->shell(Q);
                    const int Qsize = Qshell.nfunction();
                    const int Qoff = Qshell.function_index();
                    const int relQoff = Qoff - firstQoff;
                    const auto& Rshell = primary_->shell(R);
                    const int Rsize = Rshell.nfunction();
                    const int Roff = Rshell.function_index();
                    const int relRoff = Roff - firstRoff;
                    for (int q = 0; q < Qsize; q++) {
                        const int qRel = q + relQoff;
                        const int qAbs = q + Qoff;
                        for (int r = 0; r < Rsize; r++) {
                            const int rRel = r + relRoff;
                            const int rAbs = r + Roff;
#pragma omp atomic
                            Kp[qAbs][rAbs] += K3p[qRel * totalRsize + rRel];
                            if (!lr_symmetric_) {
#pragma omp atomic
                                Kp[rAbs][qAbs] += K7p[rRel * totalQsize + qRel];
                            }
                        }
                    }
                }
            }
            // > K_QS < //

            for (int Q = firstQ; Q <= lastQ; ++Q) {
                for (int S = firstS; S <= lastS; ++S) {
                    const auto& Qshell = primary_->shell(Q);
                    const int Qsize = Qshell.nfunction();
                    const int Qoff = Qshell.function_index();
                    const int relQoff = Qoff - firstQoff;
                    const auto& Sshell = primary_->shell(S);
                    const int Ssize = Sshell.nfunction();
                    const int Soff = Sshell.function_index();
                    const int relSoff = Soff - firstSoff;
                    for (int q = 0; q < Qsize; q++) {
                        const int qRel = q + relQoff;
                        const int qAbs = q + Qoff;
                        for (int s = 0; s < Ssize; s++) {
                            const int sRel = s + relSoff;
                            const int sAbs = s + Soff;
#pragma omp atomic
                            Kp[qAbs][sAbs] += K4p[qRel * totalSsize + sRel];
                            if (!lr_symmetric_) {
#pragma omp atomic
                                Kp[sAbs][qAbs] += K8p[sRel * totalQsize + qRel];
                            }
                        }
                    }
                }
            }
        }
        // if (rank == 0) timer_off("JK: Atomic");
    }  // PQRS task
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
        size_t ntri = primary_->nshell() * (primary_->nshell() + 1L) / 2L;
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

}  // namespace psi
