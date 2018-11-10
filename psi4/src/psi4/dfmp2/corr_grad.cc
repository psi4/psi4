/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#include "corr_grad.h"

#include "psi4/libqt/qt.h"
#include "psi4/lib3index/3index.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/psifiles.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"

#include "psi4/libmints/molecule.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/vector.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {
namespace dfmp2 {

CorrGrad::CorrGrad(std::shared_ptr<BasisSet> primary) : primary_(primary) { common_init(); }
CorrGrad::~CorrGrad() {}
std::shared_ptr<CorrGrad> CorrGrad::build_CorrGrad(std::shared_ptr<BasisSet> primary,
                                                   std::shared_ptr<BasisSet> auxiliary) {
    Options& options = Process::environment.options;

    if (options.get_str("SCF_TYPE").find("DF") != std::string::npos) {
        DFCorrGrad* jk = new DFCorrGrad(primary, auxiliary);

        if (options["INTS_TOLERANCE"].has_changed()) jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["PRINT"].has_changed()) jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed()) jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed()) jk->set_bench(options.get_int("BENCH"));
        if (options["DF_FITTING_CONDITION"].has_changed())
            jk->set_condition(options.get_double("DF_FITTING_CONDITION"));
        if (options["DF_INTS_NUM_THREADS"].has_changed())
            jk->set_df_ints_num_threads(options.get_int("DF_INTS_NUM_THREADS"));

        return std::shared_ptr<CorrGrad>(jk);

    } else {
        throw PSIEXCEPTION("CorrGrad::build_CorrGrad: Unknown SCF Type");
    }
}
void CorrGrad::common_init() {
    print_ = 1;
    debug_ = 0;
    bench_ = 0;

    memory_ = 32000000L;
    nthreads_ = 1;
#ifdef _OPENMP
    nthreads_ = Process::environment.get_n_threads();
#endif

    cutoff_ = 0.0;
}
DFCorrGrad::DFCorrGrad(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary)
    : CorrGrad(primary), auxiliary_(auxiliary) {
    common_init();
}
DFCorrGrad::~DFCorrGrad() {}
void DFCorrGrad::common_init() {
    df_ints_num_threads_ = 1;
#ifdef _OPENMP
    df_ints_num_threads_ = Process::environment.get_n_threads();
#endif
    condition_ = 1.0E-12;
    unit_a_ = 105;
    unit_b_ = 106;
    unit_c_ = 107;
    psio_ = PSIO::shared_object();
}
void DFCorrGrad::print_header() const {
    if (print_) {
        outfile->Printf("  ==> DFCorrGrad: Density-Fitted Correlated Gradients <==\n\n");

        outfile->Printf("    OpenMP threads:    %11d\n", nthreads_);
        outfile->Printf("    Integrals threads: %11d\n", df_ints_num_threads_);
        outfile->Printf("    Memory [GiB]:      %11.3f\n", ((double) memory_ * 8.0) / (1024.0 * 1024.0 * 1024.0));
        outfile->Printf("    Schwarz Cutoff:    %11.0E\n", cutoff_);
        outfile->Printf("    Fitting Condition: %11.0E\n\n", condition_);

        outfile->Printf("   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level("outfile", print_);
    }
}
void DFCorrGrad::compute_gradient() {
    if (!(Ca_ && Cb_ && Da_ && Db_ && Dt_)) throw PSIEXCEPTION("Occupation/Density not set");

    // => Set up gradients <= //
    int natom = primary_->molecule()->natom();
    gradients_.clear();
    gradients_["Coulomb"] = std::make_shared<Matrix>("Coulomb Gradient", natom, 3);
    gradients_["Exchange"] = std::make_shared<Matrix>("Exchange Gradient", natom, 3);

    // => Build ERI Sieve <= //
    sieve_ = std::make_shared<ERISieve>(primary_, cutoff_);

    // => Open temp files <= //
    psio_->open(unit_a_, PSIO_OPEN_NEW);
    psio_->open(unit_b_, PSIO_OPEN_NEW);
    psio_->open(unit_c_, PSIO_OPEN_NEW);

    // => Gradient Construction: Get in there and kill 'em all! <= //

    timer_on("CorrGrad: Amn");
    build_Amn_terms();
    timer_off("CorrGrad: Amn");

    timer_on("CorrGrad: AB");
    build_AB_inv_terms();
    timer_off("CorrGrad: AB");

    timer_on("CorrGrad: UV");
    build_UV_terms();
    timer_off("CorrGrad: UV");

    timer_on("CorrGrad: ABx");
    build_AB_x_terms();
    timer_off("CorrGrad: ABx");

    timer_on("CorrGrad: Amnx");
    build_Amn_x_terms();
    timer_off("CorrGrad: Amnx");

    // => Close temp files <= //
    psio_->close(unit_a_, 0);
    psio_->close(unit_b_, 0);
    psio_->close(unit_c_, 0);
}
void DFCorrGrad::build_Amn_terms() {
    // => Sizing <= //

    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];
    int la = La_->colspi()[0];
    int lb = Lb_->colspi()[0];
    int ra = Ra_->colspi()[0];
    int rb = Rb_->colspi()[0];

    int nlr = 0;
    nlr = (la > nlr ? la : nlr);
    nlr = (lb > nlr ? lb : nlr);
    nlr = (ra > nlr ? ra : nlr);
    nlr = (rb > nlr ? rb : nlr);

    bool restricted = (Ca_ == Cb_);

    const std::vector<std::pair<int, int> >& shell_pairs = sieve_->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //

    int max_rows;
    int maxP = auxiliary_->max_function_per_shell();
    size_t row_cost = 0L;
    row_cost += nso * (size_t)nso;
    row_cost += nso * (size_t)na;
    row_cost += na * (size_t)nlr;
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

    auto c = std::make_shared<Vector>("c", naux);
    double* cp = c->pointer();
    auto d = std::make_shared<Vector>("d", naux);
    double* dp = d->pointer();

    SharedMatrix Amn;
    SharedMatrix Ami;
    SharedMatrix Aij;

    double** Amnp;
    double** Amip;
    double** Aijp;

    Amn = std::make_shared<Matrix>("Amn", max_rows, nso * (size_t)nso);
    Amnp = Amn->pointer();

    Ami = std::make_shared<Matrix>("Ami", max_rows, nso * (size_t)na);
    Aij = std::make_shared<Matrix>("Aij", max_rows, na * (size_t)nlr);
    Amip = Ami->pointer();
    Aijp = Aij->pointer();

    double** Dtp = Dt_->pointer();
    double** Ptp = Pt_->pointer();
    double** Cap = Ca_->pointer();
    double** Cbp = Cb_->pointer();
    double** Lap = La_->pointer();
    double** Lbp = Lb_->pointer();
    double** Rap = Ra_->pointer();
    double** Rbp = Rb_->pointer();

    // => Prestripe <= //

    psio_address next_Aila = PSIO_ZERO;
    psio_address next_Ailb = PSIO_ZERO;
    psio_address next_Aira = PSIO_ZERO;
    psio_address next_Airb = PSIO_ZERO;

    if (true) {
        if (la) {
            for (int p = 0; p < naux; p++) {
                psio_->write(unit_a_, "(A|il)", (char*)Aijp[0], sizeof(double) * na * la, next_Aila, &next_Aila);
            }
        }
        if (ra) {
            for (int p = 0; p < naux; p++) {
                psio_->write(unit_a_, "(A|ir)", (char*)Aijp[0], sizeof(double) * na * ra, next_Aira, &next_Aira);
            }
        }
    }
    if (!restricted) {
        if (lb) {
            for (int p = 0; p < naux; p++) {
                psio_->write(unit_b_, "(A|il)", (char*)Aijp[0], sizeof(double) * nb * lb, next_Ailb, &next_Ailb);
            }
        }
        if (rb) {
            for (int p = 0; p < naux; p++) {
                psio_->write(unit_b_, "(A|ir)", (char*)Aijp[0], sizeof(double) * nb * rb, next_Airb, &next_Airb);
            }
        }
    }

    next_Aila = PSIO_ZERO;
    next_Ailb = PSIO_ZERO;
    next_Aira = PSIO_ZERO;
    next_Airb = PSIO_ZERO;

    // => Integrals <= //

    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->eri()));
    }

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
#pragma omp parallel for schedule(dynamic) num_threads(df_ints_num_threads_)
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {
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

            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        Amnp[p + oP][(m + oM) * nso + (n + oN)] = Amnp[p + oP][(n + oN) * nso + (m + oM)] = *buffer++;
                    }
                }
            }
        }

        // > (A|mn) D_mn -> c_A < //
        C_DGEMV('N', np, nso * (size_t)nso, 1.0, Amnp[0], nso * (size_t)nso, Dtp[0], 1, 0.0, &cp[pstart], 1);
        C_DGEMV('N', np, nso * (size_t)nso, 1.0, Amnp[0], nso * (size_t)nso, Ptp[0], 1, 0.0, &dp[pstart], 1);

        // > Alpha < //
        if (true) {
            // > (A|mn) C_ni -> (A|mi) < //
            C_DGEMM('N', 'N', np * (size_t)nso, na, nso, 1.0, Amnp[0], nso, Cap[0], na, 0.0, Amip[0], na);

            if (la) {
// > (A|mi) C_mj -> (A|ij) < //
#pragma omp parallel for num_threads(nthreads_)
                for (int p = 0; p < np; p++) {
                    C_DGEMM('T', 'N', na, la, nso, 1.0, Amip[p], na, Lap[0], la, 0.0, &Aijp[0][p * (size_t)na * la],
                            la);
                }

                // > Stripe < //
                psio_->write(unit_a_, "(A|il)", (char*)Aijp[0], sizeof(double) * np * na * la, next_Aila, &next_Aila);
            }

            if (ra) {
// > (A|mi) C_mj -> (A|ij) < //
#pragma omp parallel for num_threads(nthreads_)
                for (int p = 0; p < np; p++) {
                    C_DGEMM('T', 'N', na, ra, nso, 1.0, Amip[p], na, Rap[0], ra, 0.0, &Aijp[0][p * (size_t)na * ra],
                            ra);
                }

                // > Stripe < //
                psio_->write(unit_a_, "(A|ir)", (char*)Aijp[0], sizeof(double) * np * na * ra, next_Aira, &next_Aira);
            }
        }

        // > Beta < //
        if (!restricted) {
            // > (A|mn) C_ni -> (A|mi) < //
            C_DGEMM('N', 'N', np * (size_t)nso, nb, nso, 1.0, Amnp[0], nso, Cbp[0], nb, 0.0, Amip[0], na);

            if (lb) {
// > (A|mi) C_mj -> (A|ij) < //
#pragma omp parallel for num_threads(nthreads_)
                for (int p = 0; p < np; p++) {
                    C_DGEMM('T', 'N', nb, lb, nso, 1.0, Amip[p], na, Lbp[0], lb, 0.0, &Aijp[0][p * (size_t)nb * lb],
                            lb);
                }

                // > Stripe < //
                psio_->write(unit_b_, "(A|il)", (char*)Aijp[0], sizeof(double) * np * nb * lb, next_Ailb, &next_Ailb);
            }

            if (rb) {
// > (A|mi) C_mj -> (A|ij) < //
#pragma omp parallel for num_threads(nthreads_)
                for (int p = 0; p < np; p++) {
                    C_DGEMM('T', 'N', nb, rb, nso, 1.0, Amip[p], na, Rbp[0], rb, 0.0, &Aijp[0][p * (size_t)nb * rb],
                            rb);
                }

                // > Stripe < //
                psio_->write(unit_b_, "(A|ir)", (char*)Aijp[0], sizeof(double) * np * nb * la, next_Airb, &next_Airb);
            }
        }
    }

    psio_->write_entry(unit_c_, "c", (char*)cp, sizeof(double) * naux);
    psio_->write_entry(unit_c_, "d", (char*)dp, sizeof(double) * naux);
}
void DFCorrGrad::build_AB_inv_terms() {
    // => Sizing <= //

    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];
    int la = La_->colspi()[0];
    int lb = Lb_->colspi()[0];
    int ra = Ra_->colspi()[0];
    int rb = Rb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    // => Fitting Metric Full Inverse <= //

    auto metric = std::make_shared<FittingMetric>(auxiliary_, true);
    metric->form_full_eig_inverse();
    SharedMatrix J = metric->get_metric();
    double** Jp = J->pointer();

    // => d_A = (A|B)^{-1} c_B <= //
    auto c = std::make_shared<Vector>("c", naux);
    auto d = std::make_shared<Vector>("d", naux);
    double* cp = c->pointer();
    double* dp = d->pointer();

    psio_->read_entry(unit_c_, "c", (char*)cp, sizeof(double) * naux);

    C_DGEMV('N', naux, naux, 1.0, Jp[0], naux, cp, 1, 0.0, dp, 1);

    psio_->write_entry(unit_c_, "c", (char*)dp, sizeof(double) * naux);

    psio_->read_entry(unit_c_, "d", (char*)cp, sizeof(double) * naux);

    C_DGEMV('N', naux, naux, 1.0, Jp[0], naux, cp, 1, 0.0, dp, 1);

    psio_->write_entry(unit_c_, "d", (char*)dp, sizeof(double) * naux);

    if (true) {
        if (la) {
            fitting_helper(J, unit_a_, "(A|il)", naux, na * (size_t)la, memory_);
        }
        if (ra) {
            fitting_helper(J, unit_a_, "(A|ir)", naux, na * (size_t)ra, memory_);
        }
    }

    if (!restricted) {
        if (lb) {
            fitting_helper(J, unit_b_, "(A|il)", naux, nb * (size_t)lb, memory_);
        }
        if (rb) {
            fitting_helper(J, unit_b_, "(A|ir)", naux, nb * (size_t)rb, memory_);
        }
    }
}
void DFCorrGrad::fitting_helper(SharedMatrix J, size_t file, const std::string& label, size_t naux, size_t nij,
                                size_t memory) {
    int max_cols;
    size_t effective_memory = memory - 1L * naux * naux;
    size_t col_cost = 2L * naux;
    size_t cols = effective_memory / col_cost;
    cols = (cols > nij ? nij : cols);
    cols = (cols < 1L ? 1L : cols);
    max_cols = (int)cols;

    auto Aij = std::make_shared<Matrix>("Aij", naux, max_cols);
    auto Bij = std::make_shared<Matrix>("Bij", naux, max_cols);
    double** Aijp = Aij->pointer();
    double** Bijp = Bij->pointer();
    double** Jp = J->pointer();

    psio_address next_Aijb = PSIO_ZERO;

    for (long int ij = 0L; ij < nij; ij += max_cols) {
        int ncols = (ij + max_cols >= nij ? nij - ij : max_cols);

        // > Read < //
        for (int Q = 0; Q < naux; Q++) {
            next_Aijb = psio_get_address(PSIO_ZERO, sizeof(double) * (Q * (size_t)nij + ij));
            psio_->read(file, label.c_str(), (char*)Aijp[Q], sizeof(double) * ncols, next_Aijb, &next_Aijb);
        }

        // > GEMM <//
        C_DGEMM('N', 'N', naux, ncols, naux, 1.0, Jp[0], naux, Aijp[0], max_cols, 0.0, Bijp[0], max_cols);

        // > Stripe < //
        for (int Q = 0; Q < naux; Q++) {
            next_Aijb = psio_get_address(PSIO_ZERO, sizeof(double) * (Q * (size_t)nij + ij));
            psio_->write(file, label.c_str(), (char*)Bijp[Q], sizeof(double) * ncols, next_Aijb, &next_Aijb);
        }
    }
}
void DFCorrGrad::build_UV_terms() {
    // => Sizing <= //

    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];
    int la = La_->colspi()[0];
    int lb = Lb_->colspi()[0];
    int ra = Ra_->colspi()[0];
    int rb = Rb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    auto V = std::make_shared<Matrix>("W", naux, naux);
    double** Vp = V->pointer();

    // => V < = //

    // > Alpha < //
    if (true) {
        if (la) {
            UV_helper(V, 1.0, unit_a_, "(A|il)", naux, na * (size_t)la, memory_);
        }
        if (ra) {
            UV_helper(V, -1.0, unit_a_, "(A|ir)", naux, na * (size_t)ra, memory_);
        }
    }

    if (!restricted) {
        if (lb) {
            UV_helper(V, 1.0, unit_b_, "(A|il)", naux, nb * (size_t)lb, memory_);
        }
        if (rb) {
            UV_helper(V, -1.0, unit_b_, "(A|ir)", naux, nb * (size_t)rb, memory_);
        }
    } else {
        V->scale(2.0);
    }
    psio_->write_entry(unit_c_, "V", (char*)Vp[0], sizeof(double) * naux * naux);
}
void DFCorrGrad::UV_helper(SharedMatrix V, double c, size_t file, const std::string& label, size_t naux, size_t nij,
                           size_t memory) {
    int max_rows;
    size_t effective_memory = memory - 1L * naux * naux;
    size_t row_cost = 2L * nij;
    size_t rows = memory_ / row_cost;
    rows = (rows > naux ? naux : rows);
    rows = (rows < 1L ? 1L : rows);
    max_rows = (int)rows;

    // => Temporary Buffers <= //

    auto Aij = std::make_shared<Matrix>("Aij", max_rows, nij);
    auto Bij = std::make_shared<Matrix>("Bij", max_rows, nij);
    double** Aijp = Aij->pointer();
    double** Bijp = Bij->pointer();
    double** Vp = V->pointer();

    psio_address next_Aij = PSIO_ZERO;
    for (int P = 0; P < naux; P += max_rows) {
        psio_address next_Bij = PSIO_ZERO;
        int nP = (P + max_rows >= naux ? naux - P : max_rows);
        psio_->read(file, label.c_str(), (char*)Aijp[0], sizeof(double) * nP * nij, next_Aij, &next_Aij);
        for (int Q = 0; Q < naux; Q += max_rows) {
            int nQ = (Q + max_rows >= naux ? naux - Q : max_rows);
            psio_->read(file, label.c_str(), (char*)Bijp[0], sizeof(double) * nQ * nij, next_Bij, &next_Bij);

            C_DGEMM('N', 'T', nP, nQ, nij, c, Aijp[0], nij, Bijp[0], nij, 1.0, &Vp[P][Q], naux);
        }
    }
}
void DFCorrGrad::build_AB_x_terms() {
    // => Sizing <= //

    int natom = primary_->molecule()->natom();
    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();

    // => Forcing Terms/Gradients <= //
    SharedMatrix V;
    SharedMatrix W;
    SharedVector c;
    SharedVector d;

    double** Vp;
    double** Wp;
    double* cp;
    double* dp;

    c = std::make_shared<Vector>("c", naux);
    cp = c->pointer();
    psio_->read_entry(unit_c_, "c", (char*)cp, sizeof(double) * naux);
    // c->print();

    d = std::make_shared<Vector>("d", naux);
    dp = d->pointer();
    psio_->read_entry(unit_c_, "d", (char*)dp, sizeof(double) * naux);
    // d->print();

    V = std::make_shared<Matrix>("V", naux, naux);
    Vp = V->pointer();
    psio_->read_entry(unit_c_, "V", (char*)Vp[0], sizeof(double) * naux * naux);

    // => Integrals <= //

    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), auxiliary_,
                                                       BasisSet::zero_ao_basis_set());
    std::vector<std::shared_ptr<TwoBodyAOInt> > Jint;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        Jint.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)));
    }

    // => Temporary Gradients <= //

    std::vector<SharedMatrix> Jtemps;
    std::vector<SharedMatrix> Ktemps;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        Jtemps.push_back(std::make_shared<Matrix>("Jtemp", natom, 3));
        Ktemps.push_back(std::make_shared<Matrix>("Ktemp", natom, 3));
    }

    std::vector<std::pair<int, int> > PQ_pairs;
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        for (int Q = 0; Q <= P; Q++) {
            PQ_pairs.push_back(std::pair<int, int>(P, Q));
        }
    }

#pragma omp parallel for schedule(dynamic) num_threads(df_ints_num_threads_)
    for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {
        int P = PQ_pairs[PQ].first;
        int Q = PQ_pairs[PQ].second;

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        Jint[thread]->compute_shell_deriv1(P, 0, Q, 0);
        const double* buffer = Jint[thread]->buffer();

        int nP = auxiliary_->shell(P).nfunction();
        int cP = auxiliary_->shell(P).ncartesian();
        int aP = auxiliary_->shell(P).ncenter();
        int oP = auxiliary_->shell(P).function_index();

        int nQ = auxiliary_->shell(Q).nfunction();
        int cQ = auxiliary_->shell(Q).ncartesian();
        int aQ = auxiliary_->shell(Q).ncenter();
        int oQ = auxiliary_->shell(Q).function_index();

        int ncart = cP * cQ;
        const double* Px = buffer + 0 * ncart;
        const double* Py = buffer + 1 * ncart;
        const double* Pz = buffer + 2 * ncart;
        const double* Qx = buffer + 3 * ncart;
        const double* Qy = buffer + 4 * ncart;
        const double* Qz = buffer + 5 * ncart;

        double perm = (P == Q ? 1.0 : 2.0);

        double** grad_Jp;
        double** grad_Kp;

        grad_Jp = Jtemps[thread]->pointer();
        grad_Kp = Ktemps[thread]->pointer();

        for (int p = 0; p < nP; p++) {
            for (int q = 0; q < nQ; q++) {
                double Uval = 0.5 * perm * (0.5 * (cp[p + oP] * dp[q + oQ] + cp[q + oQ] * dp[p + oP]));
                grad_Jp[aP][0] -= Uval * (*Px);
                grad_Jp[aP][1] -= Uval * (*Py);
                grad_Jp[aP][2] -= Uval * (*Pz);
                grad_Jp[aQ][0] -= Uval * (*Qx);
                grad_Jp[aQ][1] -= Uval * (*Qy);
                grad_Jp[aQ][2] -= Uval * (*Qz);

                double Vval = 0.5 * perm * Vp[p + oP][q + oQ];
                grad_Kp[aP][0] -= Vval * (*Px);
                grad_Kp[aP][1] -= Vval * (*Py);
                grad_Kp[aP][2] -= Vval * (*Pz);
                grad_Kp[aQ][0] -= Vval * (*Qx);
                grad_Kp[aQ][1] -= Vval * (*Qy);
                grad_Kp[aQ][2] -= Vval * (*Qz);

                Px++;
                Py++;
                Pz++;
                Qx++;
                Qy++;
                Qz++;
            }
        }
    }

    // => Temporary Gradient Reduction <= //

    // gradients_["Coulomb"]->zero();
    // gradients_["Exchange"]->zero();

    for (int t = 0; t < df_ints_num_threads_; t++) {
        gradients_["Coulomb"]->add(Jtemps[t]);
        gradients_["Exchange"]->add(Ktemps[t]);
    }

    // gradients_["Coulomb"]->print();
    // gradients_["Exchange"]->print();
}
void DFCorrGrad::build_Amn_x_terms() {
    // => Sizing <= //

    int natom = primary_->molecule()->natom();
    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];
    int la = La_->colspi()[0];
    int lb = Lb_->colspi()[0];
    int ra = Ra_->colspi()[0];
    int rb = Rb_->colspi()[0];

    int nlr = 0;
    nlr = (la > nlr ? la : nlr);
    nlr = (lb > nlr ? lb : nlr);
    nlr = (ra > nlr ? ra : nlr);
    nlr = (rb > nlr ? rb : nlr);

    bool restricted = (Ca_ == Cb_);

    const std::vector<std::pair<int, int> >& shell_pairs = sieve_->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //

    int max_rows;
    int maxP = auxiliary_->max_function_per_shell();
    size_t row_cost = 0L;
    row_cost += nso * (size_t)nso;
    row_cost += nso * (size_t)na;
    row_cost += na * (size_t)nlr;
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
    SharedVector d;
    double* dp;

    c = std::make_shared<Vector>("c", naux);
    cp = c->pointer();
    psio_->read_entry(unit_c_, "c", (char*)cp, sizeof(double) * naux);

    d = std::make_shared<Vector>("d", naux);
    dp = d->pointer();
    psio_->read_entry(unit_c_, "d", (char*)dp, sizeof(double) * naux);

    SharedMatrix Jmn;
    SharedMatrix Ami;
    SharedMatrix Aij;

    double** Jmnp;
    double** Amip;
    double** Aijp;

    Jmn = std::make_shared<Matrix>("Jmn", max_rows, nso * (size_t)nso);
    Ami = std::make_shared<Matrix>("Ami", max_rows, nso * (size_t)na);
    Aij = std::make_shared<Matrix>("Aij", max_rows, na * (size_t)nlr);
    Jmnp = Jmn->pointer();
    Amip = Ami->pointer();
    Aijp = Aij->pointer();

    double** Dtp = Dt_->pointer();
    double** Ptp = Pt_->pointer();
    double** Cap = Ca_->pointer();
    double** Cbp = Cb_->pointer();
    double** Lap = La_->pointer();
    double** Lbp = Lb_->pointer();
    double** Rap = Ra_->pointer();
    double** Rbp = Rb_->pointer();

    psio_address next_Aila = PSIO_ZERO;
    psio_address next_Ailb = PSIO_ZERO;
    psio_address next_Aira = PSIO_ZERO;
    psio_address next_Airb = PSIO_ZERO;

    // => Integrals <= //

    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)));
    }

    // => Temporary Gradients <= //

    std::vector<SharedMatrix> Jtemps;
    std::vector<SharedMatrix> Ktemps;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        Jtemps.push_back(std::make_shared<Matrix>("Jtemp", natom, 3));
        Ktemps.push_back(std::make_shared<Matrix>("Ktemp", natom, 3));
    }

    // => R/U doubling factor <= //

    double factor = (restricted ? 2.0 : 1.0);

    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {
        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop = Pstarts[block + 1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop = (Pstop == auxiliary_->nshell() ? naux : auxiliary_->shell(Pstop).function_index());
        int np = pstop - pstart;

        // => J_mn^A <= //

        // > Alpha < //
        if (true) {
            Ami->zero();

            if (la) {
                // > Stripe < //
                psio_->read(unit_a_, "(A|il)", (char*)Aijp[0], sizeof(double) * np * na * la, next_Aila, &next_Aila);

// > (A|ij) C_mi -> (A|mj) < //
#pragma omp parallel for num_threads(nthreads_)
                for (int P = 0; P < np; P++) {
                    C_DGEMM('N', 'T', nso, na, la, 1.0, Lap[0], la, &Aijp[0][P * (size_t)na * la], la, 1.0, Amip[P],
                            na);
                }
            }

            if (ra) {
                // > Stripe < //
                psio_->read(unit_a_, "(A|ir)", (char*)Aijp[0], sizeof(double) * np * na * ra, next_Aira, &next_Aira);

// > (A|ij) C_mi -> (A|mj) < //
#pragma omp parallel for num_threads(nthreads_)
                for (int P = 0; P < np; P++) {
                    C_DGEMM('N', 'T', nso, na, ra, -1.0, Rap[0], ra, &Aijp[0][P * (size_t)na * ra], ra, 1.0, Amip[P],
                            na);
                }
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N', 'T', np * (size_t)nso, nso, na, factor, Amip[0], na, Cap[0], na, 0.0, Jmnp[0], nso);
        }

        // > Beta < //
        if (!restricted) {
            Ami->zero();

            if (lb) {
                // > Stripe < //
                psio_->read(unit_b_, "(A|il)", (char*)Aijp[0], sizeof(double) * np * nb * lb, next_Ailb, &next_Ailb);

// > (A|ij) C_mi -> (A|mj) < //
#pragma omp parallel for num_threads(nthreads_)
                for (int P = 0; P < np; P++) {
                    C_DGEMM('N', 'T', nso, nb, lb, 1.0, Lbp[0], lb, &Aijp[0][P * (size_t)nb * lb], lb, 1.0, Amip[P],
                            na);
                }
            }

            if (rb) {
                // > Stripe < //
                psio_->read(unit_b_, "(A|ir)", (char*)Aijp[0], sizeof(double) * np * nb * rb, next_Airb, &next_Airb);

// > (A|ij) C_mi -> (A|mj) < //
#pragma omp parallel for num_threads(nthreads_)
                for (int P = 0; P < np; P++) {
                    C_DGEMM('N', 'T', nso, nb, rb, -1.0, Rbp[0], rb, &Aijp[0][P * (size_t)nb * rb], rb, 1.0, Amip[P],
                            na);
                }
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N', 'T', np * (size_t)nso, nso, nb, factor, Amip[0], na, Cbp[0], nb, 0.0, Jmnp[0], nso);
        }

// > Integrals < //
#pragma omp parallel for schedule(dynamic) num_threads(df_ints_num_threads_)
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {
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

            int ncart = cP * cM * cN;
            const double* Px = buffer + 0 * ncart;
            const double* Py = buffer + 1 * ncart;
            const double* Pz = buffer + 2 * ncart;
            const double* Mx = buffer + 3 * ncart;
            const double* My = buffer + 4 * ncart;
            const double* Mz = buffer + 5 * ncart;
            const double* Nx = buffer + 6 * ncart;
            const double* Ny = buffer + 7 * ncart;
            const double* Nz = buffer + 8 * ncart;

            double perm = (M == N ? 1.0 : 2.0);

            double** grad_Jp;
            double** grad_Kp;

            grad_Jp = Jtemps[thread]->pointer();
            grad_Kp = Ktemps[thread]->pointer();

            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        double Ival = 1.0 * perm *
                                      (0.5 * (dp[p + oP + pstart] * Dtp[m + oM][n + oN] +
                                              cp[p + oP + pstart] * Ptp[m + oM][n + oN]));
                        grad_Jp[aP][0] += Ival * (*Px);
                        grad_Jp[aP][1] += Ival * (*Py);
                        grad_Jp[aP][2] += Ival * (*Pz);
                        grad_Jp[aM][0] += Ival * (*Mx);
                        grad_Jp[aM][1] += Ival * (*My);
                        grad_Jp[aM][2] += Ival * (*Mz);
                        grad_Jp[aN][0] += Ival * (*Nx);
                        grad_Jp[aN][1] += Ival * (*Ny);
                        grad_Jp[aN][2] += Ival * (*Nz);

                        double Jval =
                            1.0 * perm *
                            (0.5 * (Jmnp[p + oP][(m + oM) * nso + (n + oN)] + Jmnp[p + oP][(n + oN) * nso + (m + oM)]));
                        grad_Kp[aP][0] += Jval * (*Px);
                        grad_Kp[aP][1] += Jval * (*Py);
                        grad_Kp[aP][2] += Jval * (*Pz);
                        grad_Kp[aM][0] += Jval * (*Mx);
                        grad_Kp[aM][1] += Jval * (*My);
                        grad_Kp[aM][2] += Jval * (*Mz);
                        grad_Kp[aN][0] += Jval * (*Nx);
                        grad_Kp[aN][1] += Jval * (*Ny);
                        grad_Kp[aN][2] += Jval * (*Nz);

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

    // => Temporary Gradient Reduction <= //

    // gradients_["Coulomb"]->zero();
    // gradients_["Exchange"]->zero();

    for (int t = 0; t < df_ints_num_threads_; t++) {
        gradients_["Coulomb"]->add(Jtemps[t]);
        gradients_["Exchange"]->add(Ktemps[t]);
    }

    // gradients_["Coulomb"]->print();
    // gradients_["Exchange"]->print();
}

}  // namespace dfmp2
}  // namespace psi
