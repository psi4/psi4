
/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#include "psi4/libqt/qt.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/process.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USING_BrianQC

#include <use_brian_wrapper.h>
#include <brian_macros.h>
#include <brian_common.h>
#include <brian_geom_opt.h>

extern void checkBrian();
extern BrianCookie brianCookie;
extern bool brianEnable;
extern brianInt brianRestrictionType;

#endif

using namespace psi;

namespace psi {
namespace scfgrad {

DirectJKGrad::DirectJKGrad(int deriv, std::shared_ptr<BasisSet> primary) : JKGrad(deriv, primary) { common_init(); }
DirectJKGrad::~DirectJKGrad() {}
void DirectJKGrad::common_init() {
    ints_num_threads_ = 1;
#ifdef _OPENMP
    ints_num_threads_ = Process::environment.get_n_threads();
#endif
}
void DirectJKGrad::print_header() const {
    if (print_) {
        outfile->Printf("  ==> DirectJKGrad: Integral-Direct SCF Gradients <==\n\n");

        outfile->Printf("    Gradient:          %11d\n", deriv_);
        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:             %11.3E\n", omega_);
        outfile->Printf("    Integrals threads: %11d\n", ints_num_threads_);
        outfile->Printf("    Schwarz Cutoff:    %11.0E\n", cutoff_);
        outfile->Printf("\n");
    }
}
void DirectJKGrad::compute_gradient() {
    if (!do_J_ && !do_K_ && !do_wK_) return;

    if (!(Ca_ && Cb_ && Da_ && Db_ && Dt_))
        throw PSIEXCEPTION("Occupation/Density not set");
    
#ifdef USING_BrianQC
    if (brianEnable) {
        brianBool computeCoulomb = (do_J_ ? BRIAN_TRUE : BRIAN_FALSE);
        brianBool computeExchange = ((do_K_ || do_wK_) ? BRIAN_TRUE : BRIAN_FALSE);
        bool betaFlag = (brianRestrictionType != BRIAN_RESTRICTION_TYPE_RHF);
        
        std::shared_ptr<Matrix> Jgrad, Kgrada, Kgradb;
        if (computeCoulomb) {
            Jgrad = std::make_shared<Matrix>("Coulomb Gradient", primary_->molecule()->natom(), 3);
        }
        if (computeExchange) {
            Kgrada = std::make_shared<Matrix>("Exchange Gradient", primary_->molecule()->natom(), 3);
            if (betaFlag) {
                Kgradb = std::make_shared<Matrix>("Exchange Gradient beta", primary_->molecule()->natom(), 3);
            }
        }
        
        brianOPTBuildGradientRepulsionDeriv(&brianCookie,
            &computeCoulomb,
            &computeExchange,
            Da_->get_pointer(),
            (betaFlag ? Db_->get_pointer() : nullptr),
            (computeCoulomb ? Jgrad->get_pointer() : nullptr),
            (computeExchange ? Kgrada->get_pointer() : nullptr),
            ((computeExchange && betaFlag) ? Kgradb->get_pointer() : nullptr)
        );
        
        if (computeExchange) {
            if (betaFlag) {
                Kgrada->add(Kgradb);
            } else {
                Kgrada->scale(2.0);
            }
        }
        
        gradients_.clear();
        
        if (do_J_) {
            gradients_["Coulomb"] = Jgrad;
        }
        
        if (do_K_) {
            gradients_["Exchange"] = Kgrada;
            
            if (do_wK_) {
                gradients_["Exchange,LR"] = std::make_shared<Matrix>("Exchange,LR Gradient", primary_->molecule()->natom(), 3);
            }
        } else if (do_wK_) {
            gradients_["Exchange,LR"] = Kgrada;
        }
        
        return;
    }
#endif

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

    auto factory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);

    if (do_J_ || do_K_) {
        std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
        for (int thread = 0; thread < ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri(1)));
        }
        std::map<std::string, std::shared_ptr<Matrix>> vals = compute1(ints);
        if (do_J_) {
            gradients_["Coulomb"]->copy(vals["J"]);
            // gradients_["Coulomb"]->print();
        }
        if (do_K_) {
            gradients_["Exchange"]->copy(vals["K"]);
            // gradients_["Exchange"]->print();
        }
    }
    if (do_wK_) {
        std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
        for (int thread = 0; thread < ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega_, 1)));
        }
        std::map<std::string, std::shared_ptr<Matrix>> vals = compute1(ints);
        gradients_["Exchange,LR"]->copy(vals["K"]);
        // gradients_["Exchange,LR"]->print();
    }
}
std::map<std::string, std::shared_ptr<Matrix>> DirectJKGrad::compute1(
    std::vector<std::shared_ptr<TwoBodyAOInt>>& ints) {
    int nthreads = ints.size();

    int natom = primary_->molecule()->natom();

    std::vector<std::shared_ptr<Matrix>> Jgrad;
    std::vector<std::shared_ptr<Matrix>> Kgrad;
    for (int thread = 0; thread < nthreads; thread++) {
        Jgrad.push_back(std::make_shared<Matrix>("JGrad", natom, 3));
        Kgrad.push_back(std::make_shared<Matrix>("KGrad", natom, 3));
    }

    double** Dtp = Dt_->pointer();
    double** Dap = Da_->pointer();
    double** Dbp = Db_->pointer();

    size_t computed_shells = 0L;
    // shell pair blocks
    auto blocksPQ = ints[0]->get_blocks12();
    auto blocksRS = ints[0]->get_blocks34();
    bool use_batching = ints[0]->maximum_block_size() > 1;

#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    // loop over all the blocks of (P>=Q|
    for (size_t blockPQ_idx = 0; blockPQ_idx < blocksPQ.size(); blockPQ_idx++) {
        const auto& blockPQ = blocksPQ[blockPQ_idx];
#ifdef _OPENMP
        const int rank = omp_get_thread_num();
#else
        const int rank = 0;
#endif
        double** Jp = Jgrad[rank]->pointer();
        double** Kp = Kgrad[rank]->pointer();
        // loop over all the blocks of |R>=S)
        size_t start = ints[rank]->first_RS_shell_block(blockPQ_idx);
        for (int blockRS_idx = start; blockRS_idx < blocksRS.size(); ++blockRS_idx) {
            const auto& blockRS = blocksRS[blockRS_idx];

            if (!ints[rank]->shell_block_significant(blockPQ_idx, blockRS_idx)) continue;

            // compute the integrals and continue if none were computed
            ints[rank]->compute_shell_blocks_deriv1(blockPQ_idx, blockRS_idx);
            const auto& buffers = ints[rank]->buffers();

            const double* pAx = buffers[0];
            const double* pAy = buffers[1];
            const double* pAz = buffers[2];
            const double* pBx = buffers[3];
            const double* pBy = buffers[4];
            const double* pBz = buffers[5];
            const double* pCx = buffers[6];
            const double* pCy = buffers[7];
            const double* pCz = buffers[8];
            const double* pDx = buffers[9];
            const double* pDy = buffers[10];
            const double* pDz = buffers[11];

            // Loop over all of the P,Q,R,S shells within the blocks.  We have P>=Q, R>=S and PQ<=RS.
            for (const auto& pairPQ : blockPQ) {
                const auto& P = pairPQ.first;
                const auto& Q = pairPQ.second;
                const auto& Pshell = primary_->shell(P);
                const auto& Qshell = primary_->shell(Q);
                const auto Pam = Pshell.am();
                const auto Qam = Qshell.am();
                const auto& Psize = Pshell.nfunction();
                const auto& Qsize = Qshell.nfunction();
                const auto& Poff = Pshell.function_index();
                const auto& Qoff = Qshell.function_index();
                const auto& Pcenter = Pshell.ncenter();
                const auto& Qcenter = Qshell.ncenter();

                for (const auto& pairRS : blockRS) {
                    const auto& R = pairRS.first;
                    const auto& S = pairRS.second;
                    const auto& Rshell = primary_->shell(R);
                    const auto& Sshell = primary_->shell(S);
                    const auto Ram = Rshell.am();
                    const auto Sam = Sshell.am();
                    const auto& Rsize = Rshell.nfunction();
                    const auto& Ssize = Sshell.nfunction();
                    const auto& Roff = Rshell.function_index();
                    const auto& Soff = Sshell.function_index();
                    const auto& Rcenter = Rshell.ncenter();
                    const auto& Scenter = Sshell.ncenter();

                    double prefactor = 1.0;
                    if (P != Q) prefactor *= 2.0;
                    if (R != S) prefactor *= 2.0;
                    if (P != R || Q != S) prefactor *= 2.0;

                    size_t block_size = (size_t)Psize * Qsize * Rsize * Ssize;

                    // When there are chunks of shellpairs in RS, we need to make sure
                    // we filter out redundant combinations.
                    if (use_batching && Pam == Ram && Qam == Sam && ((P > R) || (P == R && Q > S))) {
                        pAx += block_size;
                        pAy += block_size;
                        pAz += block_size;
                        pBx += block_size;
                        pBy += block_size;
                        pBz += block_size;
                        pCx += block_size;
                        pCy += block_size;
                        pCz += block_size;
                        pDx += block_size;
                        pDy += block_size;
                        pDz += block_size;
                        continue;
                    }
                    double val;
                    double Dpq, Drs;
                    size_t delta;
                    double Ax, Ay, Az;
                    double Bx, By, Bz;
                    double Cx, Cy, Cz;
                    double Dx, Dy, Dz;

                    // => Coulomb Term <= //

                    Ax = 0.0;
                    Ay = 0.0;
                    Az = 0.0;
                    Bx = 0.0;
                    By = 0.0;
                    Bz = 0.0;
                    Cx = 0.0;
                    Cy = 0.0;
                    Cz = 0.0;
                    Dx = 0.0;
                    Dy = 0.0;
                    Dz = 0.0;
                    delta = 0L;
                    for (int p = 0; p < Psize; p++) {
                        for (int q = 0; q < Qsize; q++) {
                            for (int r = 0; r < Rsize; r++) {
                                for (int s = 0; s < Ssize; s++) {
                                    Dpq = Dtp[p + Poff][q + Qoff];
                                    Drs = Dtp[r + Roff][s + Soff];
                                    val = prefactor * Dpq * Drs;
                                    Ax += val * pAx[delta];
                                    Ay += val * pAy[delta];
                                    Az += val * pAz[delta];
                                    Bx += val * pBx[delta];
                                    By += val * pBy[delta];
                                    Bz += val * pBz[delta];
                                    Cx += val * pCx[delta];
                                    Cy += val * pCy[delta];
                                    Cz += val * pCz[delta];
                                    Dx += val * pDx[delta];
                                    Dy += val * pDy[delta];
                                    Dz += val * pDz[delta];
                                    delta++;
                                }
                            }
                        }
                    }

                    Jp[Pcenter][0] += Ax;
                    Jp[Pcenter][1] += Ay;
                    Jp[Pcenter][2] += Az;
                    Jp[Qcenter][0] += Bx;
                    Jp[Qcenter][1] += By;
                    Jp[Qcenter][2] += Bz;
                    Jp[Rcenter][0] += Cx;
                    Jp[Rcenter][1] += Cy;
                    Jp[Rcenter][2] += Cz;
                    Jp[Scenter][0] += Dx;
                    Jp[Scenter][1] += Dy;
                    Jp[Scenter][2] += Dz;

                    // => Exchange Term <= //

                    Ax = 0.0;
                    Ay = 0.0;
                    Az = 0.0;
                    Bx = 0.0;
                    By = 0.0;
                    Bz = 0.0;
                    Cx = 0.0;
                    Cy = 0.0;
                    Cz = 0.0;
                    Dx = 0.0;
                    Dy = 0.0;
                    Dz = 0.0;
                    delta = 0L;
                    for (int p = 0; p < Psize; p++) {
                        for (int q = 0; q < Qsize; q++) {
                            for (int r = 0; r < Rsize; r++) {
                                for (int s = 0; s < Ssize; s++) {
                                    val = 0.0;
                                    Dpq = Dap[p + Poff][r + Roff];
                                    Drs = Dap[q + Qoff][s + Soff];
                                    val += prefactor * Dpq * Drs;
                                    Dpq = Dap[p + Poff][s + Soff];
                                    Drs = Dap[q + Qoff][r + Roff];
                                    val += prefactor * Dpq * Drs;
                                    Dpq = Dbp[p + Poff][r + Roff];
                                    Drs = Dbp[q + Qoff][s + Soff];
                                    val += prefactor * Dpq * Drs;
                                    Dpq = Dbp[p + Poff][s + Soff];
                                    Drs = Dbp[q + Qoff][r + Roff];
                                    val += prefactor * Dpq * Drs;
                                    val *= 0.5;
                                    Ax += val * pAx[delta];
                                    Ay += val * pAy[delta];
                                    Az += val * pAz[delta];
                                    Bx += val * pBx[delta];
                                    By += val * pBy[delta];
                                    Bz += val * pBz[delta];
                                    Cx += val * pCx[delta];
                                    Cy += val * pCy[delta];
                                    Cz += val * pCz[delta];
                                    Dx += val * pDx[delta];
                                    Dy += val * pDy[delta];
                                    Dz += val * pDz[delta];
                                    delta++;
                                }
                            }
                        }
                    }

                    Kp[Pcenter][0] += Ax;
                    Kp[Pcenter][1] += Ay;
                    Kp[Pcenter][2] += Az;
                    Kp[Qcenter][0] += Bx;
                    Kp[Qcenter][1] += By;
                    Kp[Qcenter][2] += Bz;
                    Kp[Rcenter][0] += Cx;
                    Kp[Rcenter][1] += Cy;
                    Kp[Rcenter][2] += Cz;
                    Kp[Scenter][0] += Dx;
                    Kp[Scenter][1] += Dy;
                    Kp[Scenter][2] += Dz;

                    pAx += block_size;
                    pAy += block_size;
                    pAz += block_size;
                    pBx += block_size;
                    pBy += block_size;
                    pBz += block_size;
                    pCx += block_size;
                    pCy += block_size;
                    pCz += block_size;
                    pDx += block_size;
                    pDy += block_size;
                    pDz += block_size;
                }  // pairRS
            }      // pairPQ
        }          // blockRS
    }              // blockPQ

    for (int thread = 1; thread < nthreads; thread++) {
        Jgrad[0]->add(Jgrad[thread]);
        Kgrad[0]->add(Kgrad[thread]);
    }

    Jgrad[0]->scale(0.5);
    Kgrad[0]->scale(0.5);

    std::map<std::string, std::shared_ptr<Matrix>> val;
    val["J"] = Jgrad[0];
    val["K"] = Kgrad[0];
    return val;
}
void DirectJKGrad::compute_hessian() {
    if (!do_J_ && !do_K_ && !do_wK_) return;

    if (!(Ca_ && Cb_ && Da_ && Db_ && Dt_)) throw PSIEXCEPTION("Occupation/Density not set");

    // => Set up hessians <= //
    int natom = primary_->molecule()->natom();
    hessians_.clear();
    if (do_J_) {
        hessians_["Coulomb"] = std::make_shared<Matrix>("Coulomb Hessian", 3 * natom, 3 * natom);
    }
    if (do_K_) {
        hessians_["Exchange"] = std::make_shared<Matrix>("Exchange Hessian", 3 * natom, 3 * natom);
    }
    if (do_wK_) {
        hessians_["Exchange,LR"] = std::make_shared<Matrix>("Exchange,LR Hessian", 3 * natom, 3 * natom);
    }

    auto factory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);

    if (do_J_ || do_K_) {
        std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
        for (int thread = 0; thread < ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri(2)));
        }
        std::map<std::string, std::shared_ptr<Matrix>> vals = compute2(ints);
        if (do_J_) {
            hessians_["Coulomb"]->copy(vals["J"]);
        }
        if (do_K_) {
            hessians_["Exchange"]->copy(vals["K"]);
        }
    }
    if (do_wK_) {
        std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
        for (int thread = 0; thread < ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega_, 2)));
        }
        std::map<std::string, std::shared_ptr<Matrix>> vals = compute2(ints);
        hessians_["Exchange,LR"]->copy(vals["K"]);
    }
}
std::map<std::string, std::shared_ptr<Matrix>> DirectJKGrad::compute2(
    std::vector<std::shared_ptr<TwoBodyAOInt>>& ints) {
    int nthreads = ints.size();

    int natom = primary_->molecule()->natom();

    std::vector<std::shared_ptr<Matrix>> Jhess;
    std::vector<std::shared_ptr<Matrix>> Khess;
    for (int thread = 0; thread < nthreads; thread++) {
        Jhess.push_back(std::make_shared<Matrix>("JHess", 3 * natom, 3 * natom));
        Khess.push_back(std::make_shared<Matrix>("KHess", 3 * natom, 3 * natom));
    }

    double** Dtp = Dt_->pointer();
    double** Dap = Da_->pointer();
    double** Dbp = Db_->pointer();

    size_t computed_shells = 0L;
    // shell pair blocks
    auto blocksPQ = ints[0]->get_blocks12();
    auto blocksRS = ints[0]->get_blocks34();
    bool use_batching = blocksPQ != blocksRS;

#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
    // loop over all the blocks of (P>=Q|
    for (size_t blockPQ_idx = 0; blockPQ_idx < blocksPQ.size(); blockPQ_idx++) {
        const auto& blockPQ = blocksPQ[blockPQ_idx];
#ifdef _OPENMP
        const int rank = omp_get_thread_num();
#else
        const int rank = 0;
#endif
        double** Jp = Jhess[rank]->pointer();
        double** Kp = Khess[rank]->pointer();
        const auto& buffers = ints[rank]->buffers();
        // loop over all the blocks of |R>=S)
        int loop_start = use_batching ? 0 : blockPQ_idx;
        for (int blockRS_idx = loop_start; blockRS_idx < blocksRS.size(); ++blockRS_idx) {
            const auto& blockRS = blocksRS[blockRS_idx];

            // This is where we want to screen with density and schwarz-like screening

            // compute the integrals and continue if none were computed
            ints[rank]->compute_shell_blocks_deriv2(blockPQ_idx, blockRS_idx);

            std::array<const double*, 78> bufptrs;
            for (int buf = 0; buf < 78; ++buf) bufptrs[buf] = buffers[buf];

            // Loop over all of the P,Q,R,S shells within the blocks.  We have P>=Q, R>=S and PQ<=RS.
            for (const auto& pairPQ : blockPQ) {
                const auto& P = pairPQ.first;
                const auto& Q = pairPQ.second;
                const auto& Pshell = primary_->shell(P);
                const auto& Qshell = primary_->shell(Q);
                const auto& Psize = Pshell.nfunction();
                const auto& Qsize = Qshell.nfunction();
                const auto& Poff = Pshell.function_index();
                const auto& Qoff = Qshell.function_index();
                const auto& Pcenter = Pshell.ncenter();
                const auto& Qcenter = Qshell.ncenter();

                for (const auto& pairRS : blockRS) {
                    const auto& R = pairRS.first;
                    const auto& S = pairRS.second;
                    const auto& Rshell = primary_->shell(R);
                    const auto& Sshell = primary_->shell(S);
                    const auto& Rsize = Rshell.nfunction();
                    const auto& Ssize = Sshell.nfunction();
                    const auto& Roff = Rshell.function_index();
                    const auto& Soff = Sshell.function_index();
                    const auto& Rcenter = Rshell.ncenter();
                    const auto& Scenter = Sshell.ncenter();

                    size_t block_size = (size_t)Psize * Qsize * Rsize * Ssize;

                    // When there are chunks of shellpairs in RS, we need to make sure
                    // we filter out redundant combinations.  This should probably be done
                    // by having a block of RS generated for each PQ at list build time.
                    if (use_batching && ((P > R) || (P == R && Q > S))) {
                        for (int buf = 0; buf < 78; ++buf) bufptrs[buf] += block_size;
                        continue;
                    }
                    double PQscale = Pcenter == Qcenter ? 2.0 : 1.0;
                    double PRscale = Pcenter == Rcenter ? 2.0 : 1.0;
                    double PSscale = Pcenter == Scenter ? 2.0 : 1.0;
                    double QRscale = Qcenter == Rcenter ? 2.0 : 1.0;
                    double QSscale = Qcenter == Scenter ? 2.0 : 1.0;
                    double RSscale = Rcenter == Scenter ? 2.0 : 1.0;

                    int Px = 3 * Pcenter + 0;
                    int Py = 3 * Pcenter + 1;
                    int Pz = 3 * Pcenter + 2;

                    int Qx = 3 * Qcenter + 0;
                    int Qy = 3 * Qcenter + 1;
                    int Qz = 3 * Qcenter + 2;

                    int Rx = 3 * Rcenter + 0;
                    int Ry = 3 * Rcenter + 1;
                    int Rz = 3 * Rcenter + 2;

                    int Sx = 3 * Scenter + 0;
                    int Sy = 3 * Scenter + 1;
                    int Sz = 3 * Scenter + 2;

                    double prefactor = 0.5;
                    if (P != Q) prefactor *= 2.0;
                    if (R != S) prefactor *= 2.0;
                    if (P != R || Q != S) prefactor *= 2.0;

                    double val;
                    std::array<double, 78> contributions;
                    double Dpq, Drs;
                    size_t delta;

                    // => Coulomb Term <= //
                    delta = 0L;
                    contributions.fill(0.0);
                    for (int p = 0; p < Psize; p++) {
                        for (int q = 0; q < Qsize; q++) {
                            for (int r = 0; r < Rsize; r++) {
                                for (int s = 0; s < Ssize; s++) {
                                    Dpq = Dtp[p + Poff][q + Qoff];
                                    Drs = Dtp[r + Roff][s + Soff];
                                    val = prefactor * Dpq * Drs;
                                    for (int buf = 0; buf < 78; ++buf) {
                                        contributions[buf] += val * bufptrs[buf][delta];
                                    }
                                    delta++;
                                }
                            }
                        }
                    }

                    Jp[Px][Px] += contributions[0];
                    Jp[Px][Py] += contributions[1];
                    Jp[Px][Pz] += contributions[2];
                    Jp[Px][Qx] += contributions[3] * PQscale;
                    Jp[Px][Qy] += contributions[4];
                    Jp[Px][Qz] += contributions[5];
                    Jp[Px][Rx] += contributions[6] * PRscale;
                    Jp[Px][Ry] += contributions[7];
                    Jp[Px][Rz] += contributions[8];
                    Jp[Px][Sx] += contributions[9] * PSscale;
                    Jp[Px][Sy] += contributions[10];
                    Jp[Px][Sz] += contributions[11];
                    Jp[Py][Py] += contributions[12];
                    Jp[Py][Pz] += contributions[13];
                    Jp[Py][Qx] += contributions[14];
                    Jp[Py][Qy] += contributions[15] * PQscale;
                    Jp[Py][Qz] += contributions[16];
                    Jp[Py][Rx] += contributions[17];
                    Jp[Py][Ry] += contributions[18] * PRscale;
                    Jp[Py][Rz] += contributions[19];
                    Jp[Py][Sx] += contributions[20];
                    Jp[Py][Sy] += contributions[21] * PSscale;
                    Jp[Py][Sz] += contributions[22];
                    Jp[Pz][Pz] += contributions[23];
                    Jp[Pz][Qx] += contributions[24];
                    Jp[Pz][Qy] += contributions[25];
                    Jp[Pz][Qz] += contributions[26] * PQscale;
                    Jp[Pz][Rx] += contributions[27];
                    Jp[Pz][Ry] += contributions[28];
                    Jp[Pz][Rz] += contributions[29] * PRscale;
                    Jp[Pz][Sx] += contributions[30];
                    Jp[Pz][Sy] += contributions[31];
                    Jp[Pz][Sz] += contributions[32] * PSscale;
                    Jp[Qx][Qx] += contributions[33];
                    Jp[Qx][Qy] += contributions[34];
                    Jp[Qx][Qz] += contributions[35];
                    Jp[Qx][Rx] += contributions[36] * QRscale;
                    Jp[Qx][Ry] += contributions[37];
                    Jp[Qx][Rz] += contributions[38];
                    Jp[Qx][Sx] += contributions[39] * QSscale;
                    Jp[Qx][Sy] += contributions[40];
                    Jp[Qx][Sz] += contributions[41];
                    Jp[Qy][Qy] += contributions[42];
                    Jp[Qy][Qz] += contributions[43];
                    Jp[Qy][Rx] += contributions[44];
                    Jp[Qy][Ry] += contributions[45] * QRscale;
                    Jp[Qy][Rz] += contributions[46];
                    Jp[Qy][Sx] += contributions[47];
                    Jp[Qy][Sy] += contributions[48] * QSscale;
                    Jp[Qy][Sz] += contributions[49];
                    Jp[Qz][Qz] += contributions[50];
                    Jp[Qz][Rx] += contributions[51];
                    Jp[Qz][Ry] += contributions[52];
                    Jp[Qz][Rz] += contributions[53] * QRscale;
                    Jp[Qz][Sx] += contributions[54];
                    Jp[Qz][Sy] += contributions[55];
                    Jp[Qz][Sz] += contributions[56] * QSscale;
                    Jp[Rx][Rx] += contributions[57];
                    Jp[Rx][Ry] += contributions[58];
                    Jp[Rx][Rz] += contributions[59];
                    Jp[Rx][Sx] += contributions[60] * RSscale;
                    Jp[Rx][Sy] += contributions[61];
                    Jp[Rx][Sz] += contributions[62];
                    Jp[Ry][Ry] += contributions[63];
                    Jp[Ry][Rz] += contributions[64];
                    Jp[Ry][Sx] += contributions[65];
                    Jp[Ry][Sy] += contributions[66] * RSscale;
                    Jp[Ry][Sz] += contributions[67];
                    Jp[Rz][Rz] += contributions[68];
                    Jp[Rz][Sx] += contributions[69];
                    Jp[Rz][Sy] += contributions[70];
                    Jp[Rz][Sz] += contributions[71] * RSscale;
                    Jp[Sx][Sx] += contributions[72];
                    Jp[Sx][Sy] += contributions[73];
                    Jp[Sx][Sz] += contributions[74];
                    Jp[Sy][Sy] += contributions[75];
                    Jp[Sy][Sz] += contributions[76];
                    Jp[Sz][Sz] += contributions[77];

                    // => Exchange Term <= //
                    delta = 0L;
                    contributions.fill(0);
                    for (int p = 0; p < Psize; p++) {
                        for (int q = 0; q < Qsize; q++) {
                            for (int r = 0; r < Rsize; r++) {
                                for (int s = 0; s < Ssize; s++) {
                                    val = 0.0;
                                    Dpq = Dap[p + Poff][r + Roff];
                                    Drs = Dap[q + Qoff][s + Soff];
                                    val += prefactor * Dpq * Drs;
                                    Dpq = Dap[p + Poff][s + Soff];
                                    Drs = Dap[q + Qoff][r + Roff];
                                    val += prefactor * Dpq * Drs;
                                    Dpq = Dbp[p + Poff][r + Roff];
                                    Drs = Dbp[q + Qoff][s + Soff];
                                    val += prefactor * Dpq * Drs;
                                    Dpq = Dbp[p + Poff][s + Soff];
                                    Drs = Dbp[q + Qoff][r + Roff];
                                    val += prefactor * Dpq * Drs;
                                    val *= 0.5;
                                    for (int buf = 0; buf < 78; ++buf) {
                                        contributions[buf] += val * bufptrs[buf][delta];
                                    }
                                    delta++;
                                }
                            }
                        }
                    }

                    Kp[Px][Px] += contributions[0];
                    Kp[Px][Py] += contributions[1];
                    Kp[Px][Pz] += contributions[2];
                    Kp[Px][Qx] += contributions[3] * PQscale;
                    Kp[Px][Qy] += contributions[4];
                    Kp[Px][Qz] += contributions[5];
                    Kp[Px][Rx] += contributions[6] * PRscale;
                    Kp[Px][Ry] += contributions[7];
                    Kp[Px][Rz] += contributions[8];
                    Kp[Px][Sx] += contributions[9] * PSscale;
                    Kp[Px][Sy] += contributions[10];
                    Kp[Px][Sz] += contributions[11];
                    Kp[Py][Py] += contributions[12];
                    Kp[Py][Pz] += contributions[13];
                    Kp[Py][Qx] += contributions[14];
                    Kp[Py][Qy] += contributions[15] * PQscale;
                    Kp[Py][Qz] += contributions[16];
                    Kp[Py][Rx] += contributions[17];
                    Kp[Py][Ry] += contributions[18] * PRscale;
                    Kp[Py][Rz] += contributions[19];
                    Kp[Py][Sx] += contributions[20];
                    Kp[Py][Sy] += contributions[21] * PSscale;
                    Kp[Py][Sz] += contributions[22];
                    Kp[Pz][Pz] += contributions[23];
                    Kp[Pz][Qx] += contributions[24];
                    Kp[Pz][Qy] += contributions[25];
                    Kp[Pz][Qz] += contributions[26] * PQscale;
                    Kp[Pz][Rx] += contributions[27];
                    Kp[Pz][Ry] += contributions[28];
                    Kp[Pz][Rz] += contributions[29] * PRscale;
                    Kp[Pz][Sx] += contributions[30];
                    Kp[Pz][Sy] += contributions[31];
                    Kp[Pz][Sz] += contributions[32] * PSscale;
                    Kp[Qx][Qx] += contributions[33];
                    Kp[Qx][Qy] += contributions[34];
                    Kp[Qx][Qz] += contributions[35];
                    Kp[Qx][Rx] += contributions[36] * QRscale;
                    Kp[Qx][Ry] += contributions[37];
                    Kp[Qx][Rz] += contributions[38];
                    Kp[Qx][Sx] += contributions[39] * QSscale;
                    Kp[Qx][Sy] += contributions[40];
                    Kp[Qx][Sz] += contributions[41];
                    Kp[Qy][Qy] += contributions[42];
                    Kp[Qy][Qz] += contributions[43];
                    Kp[Qy][Rx] += contributions[44];
                    Kp[Qy][Ry] += contributions[45] * QRscale;
                    Kp[Qy][Rz] += contributions[46];
                    Kp[Qy][Sx] += contributions[47];
                    Kp[Qy][Sy] += contributions[48] * QSscale;
                    Kp[Qy][Sz] += contributions[49];
                    Kp[Qz][Qz] += contributions[50];
                    Kp[Qz][Rx] += contributions[51];
                    Kp[Qz][Ry] += contributions[52];
                    Kp[Qz][Rz] += contributions[53] * QRscale;
                    Kp[Qz][Sx] += contributions[54];
                    Kp[Qz][Sy] += contributions[55];
                    Kp[Qz][Sz] += contributions[56] * QSscale;
                    Kp[Rx][Rx] += contributions[57];
                    Kp[Rx][Ry] += contributions[58];
                    Kp[Rx][Rz] += contributions[59];
                    Kp[Rx][Sx] += contributions[60] * RSscale;
                    Kp[Rx][Sy] += contributions[61];
                    Kp[Rx][Sz] += contributions[62];
                    Kp[Ry][Ry] += contributions[63];
                    Kp[Ry][Rz] += contributions[64];
                    Kp[Ry][Sx] += contributions[65];
                    Kp[Ry][Sy] += contributions[66] * RSscale;
                    Kp[Ry][Sz] += contributions[67];
                    Kp[Rz][Rz] += contributions[68];
                    Kp[Rz][Sx] += contributions[69];
                    Kp[Rz][Sy] += contributions[70];
                    Kp[Rz][Sz] += contributions[71] * RSscale;
                    Kp[Sx][Sx] += contributions[72];
                    Kp[Sx][Sy] += contributions[73];
                    Kp[Sx][Sz] += contributions[74];
                    Kp[Sy][Sy] += contributions[75];
                    Kp[Sy][Sz] += contributions[76];
                    Kp[Sz][Sz] += contributions[77];

                    for (auto& buf : bufptrs) buf += block_size;
                }  // pairRS
            }      // pairPQ
        }          // blockRS
    }              // blockPQ

    for (int thread = 1; thread < nthreads; thread++) {
        Jhess[0]->add(Jhess[thread]);
        Khess[0]->add(Khess[thread]);
    }
    int dim = Jhess[0]->rowdim();
    double** Jp = Jhess[0]->pointer();
    double** Kp = Khess[0]->pointer();
    for (int row = 0; row < dim; ++row) {
        for (int col = 0; col < row; ++col) {
            Jp[row][col] = Jp[col][row] = (Jp[row][col] + Jp[col][row]);
            Kp[row][col] = Kp[col][row] = (Kp[row][col] + Kp[col][row]);
        }
    }
    Jhess[0]->print();
    Khess[0]->print();

    std::map<std::string, std::shared_ptr<Matrix>> val;
    val["J"] = Jhess[0];
    val["K"] = Khess[0];
    return val;
}
}
}
