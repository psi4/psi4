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
#include "psi4/lib3index/3index.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/vector.h"
#include "psi4/liboptions/liboptions.h"
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

JKGrad::JKGrad(int deriv, std::shared_ptr<BasisSet> primary) :
    deriv_(deriv), primary_(primary)
{
    common_init();
}
JKGrad::~JKGrad()
{
}
std::shared_ptr<JKGrad> JKGrad::build_JKGrad(int deriv, std::shared_ptr<MintsHelper> mints)
{
    Options& options = Process::environment.options;

    if (options.get_str("SCF_TYPE") == "DFDIRJ+COSX") {
        DFJCOSKGrad* jk = new DFJCOSKGrad(deriv, mints, options);

        if (options["INTS_TOLERANCE"].has_changed())
            jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["PRINT"].has_changed())
            jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed())
            jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed())
            jk->set_bench(options.get_int("BENCH"));
        jk->set_condition(options.get_double("DF_FITTING_CONDITION"));
        if (options["DF_INTS_NUM_THREADS"].has_changed())
            jk->set_ints_num_threads(options.get_int("DF_INTS_NUM_THREADS"));

        return std::shared_ptr<JKGrad>(jk);
    } else if (options.get_str("SCF_TYPE").find("DF") != std::string::npos) {
        DFJKGrad* jk = new DFJKGrad(deriv, mints);

        if (options["INTS_TOLERANCE"].has_changed())
            jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["PRINT"].has_changed())
            jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed())
            jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed())
            jk->set_bench(options.get_int("BENCH"));
        jk->set_condition(options.get_double("DF_FITTING_CONDITION"));
        if (options["DF_INTS_NUM_THREADS"].has_changed())
            jk->set_df_ints_num_threads(options.get_int("DF_INTS_NUM_THREADS"));

        return std::shared_ptr<JKGrad>(jk);
    } else if (options.get_str("SCF_TYPE") == "DIRECT" || options.get_str("SCF_TYPE") == "PK" || options.get_str("SCF_TYPE") == "OUT_OF_CORE") {

        DirectJKGrad* jk = new DirectJKGrad(deriv, mints->get_basisset("ORBITAL"));

        if (options["INTS_TOLERANCE"].has_changed())
            jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["PRINT"].has_changed())
            jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed())
            jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed())
            jk->set_bench(options.get_int("BENCH"));
        // TODO: rename every DF case
        if (options["DF_INTS_NUM_THREADS"].has_changed())
            jk->set_ints_num_threads(options.get_int("DF_INTS_NUM_THREADS"));

        return std::shared_ptr<JKGrad>(jk);

    } else {
        throw PSIEXCEPTION("JKGrad::build_JKGrad: Unknown SCF Type");
    }
}
void JKGrad::common_init() {
    print_ = 1;
    debug_ = 0;
    bench_ = 0;

    memory_ = 32000000L;
    omp_num_threads_ = 1;
#ifdef _OPENMP
    omp_num_threads_ = Process::environment.get_n_threads();
#endif

    cutoff_ = Process::environment.options.get_double("INTS_TOLERANCE");

    do_J_ = true;
    do_K_ = true;
    do_wK_ = false;
    omega_ = 0.0;
}
DFJKGrad::DFJKGrad(int deriv, std::shared_ptr<MintsHelper> mints) :
    JKGrad(deriv,mints->get_basisset("ORBITAL")), auxiliary_(mints->get_basisset("DF_BASIS_SCF")), mints_(mints)
{
    common_init();
}
DFJKGrad::~DFJKGrad() {}
void DFJKGrad::common_init() {
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
void DFJKGrad::print_header() const {
    if (print_) {
        outfile->Printf("  ==> DFJKGrad: Density-Fitted SCF Gradients <==\n\n");

        outfile->Printf("    Gradient:          %11d\n", deriv_);
        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:             %11.3E\n", omega_);
        outfile->Printf("    OpenMP threads:    %11d\n", omp_num_threads_);
        outfile->Printf("    Integrals threads: %11d\n", df_ints_num_threads_);
        outfile->Printf("    Memory [MiB]:      %11ld\n", (memory_ * 8L) / (1024L * 1024L));
        outfile->Printf("    Schwarz Cutoff:    %11.0E\n", cutoff_);
        outfile->Printf("    Fitting Condition: %11.0E\n\n", condition_);

        outfile->Printf("   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level("outfile", print_);
    }
}
void DFJKGrad::compute_gradient() {
    if (!do_J_ && !do_K_ && !do_wK_) return;

    if (!(Ca_ && Cb_ && Da_ && Db_ && Dt_)) throw PSIEXCEPTION("Occupation/Density not set");

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
        // throw PSIEXCEPTION("Exchange,LR gradients are not currently available with DF.");
        gradients_["Exchange,LR"] = std::make_shared<Matrix>("Exchange,LR Gradient", natom, 3);
    }

#ifdef USING_BrianQC
    if (brianEnable) {
        double threshold = cutoff_ * 1e-2;
        brianCOMSetPrecisionThresholds(&brianCookie, &threshold);
        checkBrian();
    }
#endif

    // => Open temp files <= //
    psio_->open(unit_a_, PSIO_OPEN_NEW);
    psio_->open(unit_b_, PSIO_OPEN_NEW);
    psio_->open(unit_c_, PSIO_OPEN_NEW);

    // => Gradient Construction: Get in there and kill 'em all! <= //

    // Using

    // d Minv          d M
    // ------ = -Minv ----- Minv
    //   dx             dx

    // we get

    // d (mn|w|A) Minv[A][B] (B|rs)
    // ---------------------------- =
    //             dx
    //                              x
    //                      (mn|w|A)  Minv[A][B] (B|rs)

    //                                            x
    //                   -  (mn|w|A) Minv[A][B] M[B][C] Minv[C][D] (D|rs)

    //                                                 x
    //                   +  (mn|w|A)  Minv[A][B] (B|rs)

    // c_B = (B|pq) D_pq
    // (B|ij) = (B|pq) C_pi C_qi
    timer_on("JKGrad: Amn");
    build_Amn_terms();
    timer_off("JKGrad: Amn");

    // (B|w|ij) = (B|w|pq) C_pi C_qi
    timer_on("JKGrad: Awmn");
    build_Amn_lr_terms();
    timer_off("JKGrad: Awmn");

    // c_A = (A|B)^-1 c_B
    // (A|ij) = (A|B)^-1 (B|ij)
    // (A|w|ij) = (A|B)^-1 (B|w|ij)
    timer_on("JKGrad: AB");
    build_AB_inv_terms();
    timer_off("JKGrad: AB");

    // W_AB = (A|ij) (B|ij)
    // V_AB = 0.5 * [(A|w|ij) (B|ij) + (A|ij) (B|w|ij)]
    timer_on("JKGrad: UV");
    build_UV_terms();
    timer_off("JKGrad: UV");

    //  J^x = 0.5 * (A|B)^x c_B c_A
    //  K^x = 0.5 * (A|B)^x W_AB
    // wK^x = 0.5 * (A|B)^x V_AB
    timer_on("JKGrad: ABx");
    build_AB_x_terms();
    timer_off("JKGrad: ABx");

    //  (A|pq)   = (A|ij) C_ip C_iq
    //  (A|w|pq) = (A|w|ij) C_ip C_iq
    //  J^x = (A|pq)^x d_A Dt_pq
    //  K^x = (A|pq)^x (A|pq)
    // wK^x = 0.5 * (A|pq)^x (A|w|pq)
    // wK^x = 0.5 * (A|w|pq)^x (A|pq)
    timer_on("JKGrad: Amnx");
    build_Amn_x_terms();
    timer_off("JKGrad: Amnx");

    //  J^x  = 0.5 * (C|D)^x [(C|B)^-1 (B|pq) D_tpq] [(D|A)^-1 (A|rs) Dt_rs]
    //  J^x += (A|pq)^x D_pq  [(A|B)^-1 (B|rs) Dt_rs]

    //  K^x  = 0.5 * (A|B)^x (ij|A)(A|C)^-1(C|B)^-1(B|ij)
    //  K^x += (A|pq)^x (A|pq)

    // wK^x  = 0.5 * (A|B)^x [(ij|A)(A|C)^-1] [(C|B)^-1(B|w|ij)]
    // wK^x += 0.5 * (A|pq)^x (A|w|pq)
    // wK^x += 0.5 * (A|w|pq)^x (A|pq)

    // Printing
    // gradients_["Coulomb"]->print();
    // if (do_K_) {
    //     gradients_["Exchange"]->print();
    // }
    // if (do_wK_) {
    //     gradients_["Exchange,LR"]->print();
    // }

    // => Close temp files <= //
    psio_->close(unit_a_, 0);
    psio_->close(unit_b_, 0);
    psio_->close(unit_c_, 0);
}
void DFJKGrad::build_Amn_terms() {
    // => Sizing <= //

    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    // => Integrals <= //

    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(df_ints_num_threads_);
    eri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    for (int t = 1; t < df_ints_num_threads_; t++) {
        eri[t] = std::shared_ptr<TwoBodyAOInt>(eri.front()->clone());
    }

    const std::vector<std::pair<int, int>>& shell_pairs = eri[0]->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //

    int max_rows;
    int maxP = auxiliary_->max_function_per_shell();
    size_t row_cost = 0L;
    row_cost += nso * (size_t)nso;
    if (do_K_ || do_wK_) {
        row_cost += nso * (size_t)na;
        row_cost += na * (size_t)na;
    }
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

    if (do_J_) {
        c = std::make_shared<Vector>("c", naux);
        cp = c->pointer();
    }

    SharedMatrix Amn;
    SharedMatrix Ami;
    SharedMatrix Aij;

    double** Amnp;
    double** Amip;
    double** Aijp;

    if (true) {
        Amn = std::make_shared<Matrix>("Amn", max_rows, nso * (size_t)nso);
        Amnp = Amn->pointer();
    }
    if (do_K_ || do_wK_) {
        Ami = std::make_shared<Matrix>("Ami", max_rows, nso * (size_t)na);
        Aij = std::make_shared<Matrix>("Aij", max_rows, na * (size_t)na);
        Amip = Ami->pointer();
        Aijp = Aij->pointer();
    }

    double** Dtp = Dt_->pointer();
    double** Cap = Ca_->pointer();
    double** Cbp = Cb_->pointer();

    psio_address next_Aija = PSIO_ZERO;
    psio_address next_Aijb = PSIO_ZERO;

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
        int nthread_df = df_ints_num_threads_;
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

            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        Amnp[p + oP][(m + oM) * nso + (n + oN)] = Amnp[p + oP][(n + oN) * nso + (m + oM)] = *buffer++;
                    }
                }
            }
        }

        // > (A|mn) D_mn -> c_A < //
        if (do_J_) {
            C_DGEMV('N', np, nso * (size_t)nso, 1.0, Amnp[0], nso * (size_t)nso, Dtp[0], 1, 0.0, &cp[pstart], 1);
        }

        // > Alpha < //
        if (do_K_ || do_wK_) {
            // > (A|mn) C_ni -> (A|mi) < //
            C_DGEMM('N', 'N', np * (size_t)nso, na, nso, 1.0, Amnp[0], nso, Cap[0], na, 0.0, Amip[0], na);

            // > (A|mi) C_mj -> (A|ij) < //
#pragma omp parallel for
            for (int p = 0; p < np; p++) {
                C_DGEMM('T', 'N', na, na, nso, 1.0, Amip[p], na, Cap[0], na, 0.0, &Aijp[0][p * (size_t)na * na], na);
            }

            // > Stripe < //
            psio_->write(unit_a_, "(A|ij)", (char*)Aijp[0], sizeof(double) * np * na * na, next_Aija, &next_Aija);
        }

        // > Beta < //
        if (!restricted && (do_K_ || do_wK_)) {
            // skip if there are no beta electrons
            if (nb > 0){
                // > (A|mn) C_ni -> (A|mi) < //
                    C_DGEMM('N', 'N', np * (size_t)nso, nb, nso, 1.0, Amnp[0], nso, Cbp[0], nb, 0.0, Amip[0], na);

                // > (A|mi) C_mj -> (A|ij) < //
#pragma omp parallel for
                for (int p = 0; p < np; p++) {
                    C_DGEMM('T', 'N', nb, nb, nso, 1.0, Amip[p], na, Cbp[0], nb, 0.0, &Aijp[0][p * (size_t)nb * nb], nb);
                }
            }
            // > Stripe < //
            psio_->write(unit_b_, "(A|ij)", (char*)Aijp[0], sizeof(double) * np * nb * nb, next_Aijb, &next_Aijb);
        }
    }

    if (do_J_) {
        psio_->write_entry(unit_c_, "c", (char*)cp, sizeof(double) * naux);
    }
}
void DFJKGrad::build_Amn_lr_terms() {
    if (!do_wK_) return;

    // => Sizing <= //

    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    // => Integrals <= //

    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_)));
    }

    const std::vector<std::pair<int, int>>& shell_pairs = eri[0]->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //

    int max_rows;
    int maxP = auxiliary_->max_function_per_shell();
    size_t row_cost = 0L;
    row_cost += nso * (size_t)nso;
    row_cost += nso * (size_t)na;
    row_cost += na * (size_t)na;
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

    SharedMatrix Amn;
    SharedMatrix Ami;
    SharedMatrix Aij;

    double** Amnp;
    double** Amip;
    double** Aijp;

    Amn = std::make_shared<Matrix>("Amn", max_rows, nso * (size_t)nso);
    Ami = std::make_shared<Matrix>("Ami", max_rows, nso * (size_t)na);
    Aij = std::make_shared<Matrix>("Aij", max_rows, na * (size_t)na);

    Amnp = Amn->pointer();
    Amip = Ami->pointer();
    Aijp = Aij->pointer();

    double** Cap = Ca_->pointer();
    double** Cbp = Cb_->pointer();

    psio_address next_Aija = PSIO_ZERO;
    psio_address next_Aijb = PSIO_ZERO;

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
        int nthread_df = df_ints_num_threads_;
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

            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        Amnp[p + oP][(m + oM) * nso + (n + oN)] = Amnp[p + oP][(n + oN) * nso + (m + oM)] = *buffer++;
                    }
                }
            }
        }

        // > Alpha < //
        if (true) {
            // > (A|mn) C_ni -> (A|mi) < //
            C_DGEMM('N', 'N', np * (size_t)nso, na, nso, 1.0, Amnp[0], nso, Cap[0], na, 0.0, Amip[0], na);

            // > (A|mi) C_mj -> (A|ij) < //
#pragma omp parallel for
            for (int p = 0; p < np; p++) {
                C_DGEMM('T', 'N', na, na, nso, 1.0, Amip[p], na, Cap[0], na, 0.0, &Aijp[0][p * (size_t)na * na], na);
            }

            // > Stripe < //
            psio_->write(unit_a_, "(A|w|ij)", (char*)Aijp[0], sizeof(double) * np * na * na, next_Aija, &next_Aija);
        }

        // > Beta < //
        if (!restricted) {
            // > (A|mn) C_ni -> (A|mi) < //
            C_DGEMM('N', 'N', np * (size_t)nso, nb, nso, 1.0, Amnp[0], nso, Cbp[0], nb, 0.0, Amip[0], na);

            // > (A|mi) C_mj -> (A|ij) < //
#pragma omp parallel for
            for (int p = 0; p < np; p++) {
                C_DGEMM('T', 'N', nb, nb, nso, 1.0, Amip[p], na, Cbp[0], nb, 0.0, &Aijp[0][p * (size_t)nb * nb], nb);
            }

            // > Stripe < //
            psio_->write(unit_b_, "(A|w|ij)", (char*)Aijp[0], sizeof(double) * np * nb * nb, next_Aijb, &next_Aijb);
        }
    }
}
void DFJKGrad::build_AB_inv_terms() {
    // => Sizing <= //

    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    // => Fitting Metric Full Inverse <= //

    auto metric = std::make_shared<FittingMetric>(auxiliary_, true);
    metric->form_full_eig_inverse(condition_);
    SharedMatrix J = metric->get_metric();
    double** Jp = J->pointer();

    // => d_A = (A|B)^{-1} c_B <= //
    if (do_J_) {
        auto c = std::make_shared<Vector>("c", naux);
        auto d = std::make_shared<Vector>("d", naux);
        double* cp = c->pointer();
        double* dp = d->pointer();

        psio_->read_entry(unit_c_, "c", (char*)cp, sizeof(double) * naux);

        C_DGEMV('N', naux, naux, 1.0, Jp[0], naux, cp, 1, 0.0, dp, 1);

        psio_->write_entry(unit_c_, "c", (char*)dp, sizeof(double) * naux);
    }

    if (!(do_K_ || do_wK_)) return;

    int max_cols;
    size_t effective_memory = memory_ - 1L * naux * naux;
    size_t col_cost = 2L * naux;
    size_t cols = effective_memory / col_cost;
    cols = (cols > na * (size_t)na ? na * (size_t)na : cols);
    cols = (cols < na ? na : cols);
    max_cols = (int)cols;

    auto Aij = std::make_shared<Matrix>("Aij", naux, max_cols);
    auto Bij = std::make_shared<Matrix>("Bij", naux, max_cols);
    double** Aijp = Aij->pointer();
    double** Bijp = Bij->pointer();

    // K and wK 3 index temps
    std::vector<std::string> buffers;
    buffers.push_back("(A|ij)");
    if (do_wK_) buffers.push_back("(A|w|ij)");

    // Units and sizing for alpha/beta
    std::vector<std::pair<size_t, size_t>> us_vec;
    us_vec.push_back(std::make_pair(unit_a_, na));
    if (!restricted) us_vec.push_back(std::make_pair(unit_b_, nb));

    // Transform all three index buffers (A|B)(B|ij) -> (A|ij)
    for (const auto& buff_name : buffers) {
        for (const auto& us : us_vec) {
            size_t unit_name = us.first;
            size_t nmo_size = us.second;
            size_t nmo_size2 = nmo_size * nmo_size;
            // printf("%s | %zu %zu\n", buff_name.c_str(), unit_name, nmo_size);

            psio_address next_Aija = PSIO_ZERO;

            for (long int ij = 0L; ij < nmo_size2; ij += max_cols) {
                int ncols = (ij + max_cols >= nmo_size2 ? nmo_size2 - ij : max_cols);

                // > Read < //
                for (int Q = 0; Q < naux; Q++) {
                    next_Aija = psio_get_address(PSIO_ZERO, sizeof(double) * (Q * (size_t)nmo_size2 + ij));
                    psio_->read(unit_name, buff_name.c_str(), (char*)Aijp[Q], sizeof(double) * ncols, next_Aija,
                                &next_Aija);
                }

                // > GEMM <//
                C_DGEMM('N', 'N', naux, ncols, naux, 1.0, Jp[0], naux, Aijp[0], max_cols, 0.0, Bijp[0], max_cols);

                // > Stripe < //
                for (int Q = 0; Q < naux; Q++) {
                    next_Aija = psio_get_address(PSIO_ZERO, sizeof(double) * (Q * (size_t)nmo_size2 + ij));
                    psio_->write(unit_name, buff_name.c_str(), (char*)Bijp[Q], sizeof(double) * ncols, next_Aija,
                                 &next_Aija);
                }
            }
        }
    }
}
void DFJKGrad::build_UV_terms() {
    if (!(do_K_ || do_wK_)) return;

    // => Sizing <= //

    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    auto V = std::make_shared<Matrix>("W", naux, naux);
    double** Vp = V->pointer();

    // => Memory Constraints <= //

    int max_rows;
    size_t effective_memory = memory_ - 1L * naux * naux;
    size_t row_cost = 2L * na * (size_t)na;
    size_t rows = memory_ / row_cost;
    rows = (rows > naux ? naux : rows);
    rows = (rows < 1L ? 1L : rows);
    max_rows = (int)rows;

    // => Temporary Buffers <= //

    auto Aij = std::make_shared<Matrix>("Aij", max_rows, na * (size_t)na);
    auto Bij = std::make_shared<Matrix>("Bij", max_rows, na * (size_t)na);
    double** Aijp = Aij->pointer();
    double** Bijp = Bij->pointer();

    // => V < = //

    // > Alpha < //
    if (true) {
        psio_address next_Aij = PSIO_ZERO;
        for (int P = 0; P < naux; P += max_rows) {
            psio_address next_Bij = PSIO_ZERO;
            int nP = (P + max_rows >= naux ? naux - P : max_rows);
            psio_->read(unit_a_, "(A|ij)", (char*)Aijp[0], sizeof(double) * nP * na * na, next_Aij, &next_Aij);
            for (int Q = 0; Q < naux; Q += max_rows) {
                int nQ = (Q + max_rows >= naux ? naux - Q : max_rows);
                psio_->read(unit_a_, "(A|ij)", (char*)Bijp[0], sizeof(double) * nQ * na * na, next_Bij, &next_Bij);

                C_DGEMM('N', 'T', nP, nQ, na * (size_t)na, 1.0, Aijp[0], na * (size_t)na, Bijp[0], na * (size_t)na, 0.0,
                        &Vp[P][Q], naux);
            }
        }
    }
    // > Beta < //
    if (!restricted) {
        psio_address next_Aij = PSIO_ZERO;
        for (int P = 0; P < naux; P += max_rows) {
            psio_address next_Bij = PSIO_ZERO;
            int nP = (P + max_rows >= naux ? naux - P : max_rows);
            psio_->read(unit_b_, "(A|ij)", (char*)Aijp[0], sizeof(double) * nP * nb * nb, next_Aij, &next_Aij);
            for (int Q = 0; Q < naux; Q += max_rows) {
                int nQ = (Q + max_rows >= naux ? naux - Q : max_rows);
                psio_->read(unit_b_, "(A|ij)", (char*)Bijp[0], sizeof(double) * nQ * nb * nb, next_Bij, &next_Bij);

                C_DGEMM('N', 'T', nP, nQ, nb * (size_t)nb, 1.0, Aijp[0], nb * (size_t)nb, Bijp[0], nb * (size_t)nb, 1.0,
                        &Vp[P][Q], naux);
            }
        }
    } else {
        V->scale(2.0);
    }
    psio_->write_entry(unit_c_, "V", (char*)Vp[0], sizeof(double) * naux * naux);

    if (!do_wK_) return;

    // => W < = //
    V->zero();

    // > Alpha < //
    if (true) {
        psio_address next_Aij = PSIO_ZERO;
        for (int P = 0; P < naux; P += max_rows) {
            psio_address next_Bij = PSIO_ZERO;
            int nP = (P + max_rows >= naux ? naux - P : max_rows);
            psio_->read(unit_a_, "(A|ij)", (char*)Aijp[0], sizeof(double) * nP * na * na, next_Aij, &next_Aij);
            for (int Q = 0; Q < naux; Q += max_rows) {
                int nQ = (Q + max_rows >= naux ? naux - Q : max_rows);
                psio_->read(unit_a_, "(A|w|ij)", (char*)Bijp[0], sizeof(double) * nQ * na * na, next_Bij, &next_Bij);

                C_DGEMM('N', 'T', nP, nQ, na * (size_t)na, 1.0, Aijp[0], na * (size_t)na, Bijp[0], na * (size_t)na, 0.0,
                        &Vp[P][Q], naux);
            }
        }
    }
    // > Beta < //
    if (!restricted) {
        psio_address next_Aij = PSIO_ZERO;
        for (int P = 0; P < naux; P += max_rows) {
            psio_address next_Bij = PSIO_ZERO;
            int nP = (P + max_rows >= naux ? naux - P : max_rows);
            psio_->read(unit_b_, "(A|ij)", (char*)Aijp[0], sizeof(double) * nP * nb * nb, next_Aij, &next_Aij);
            for (int Q = 0; Q < naux; Q += max_rows) {
                int nQ = (Q + max_rows >= naux ? naux - Q : max_rows);
                psio_->read(unit_b_, "(A|w|ij)", (char*)Bijp[0], sizeof(double) * nQ * nb * nb, next_Bij, &next_Bij);

                C_DGEMM('N', 'T', nP, nQ, nb * (size_t)nb, 1.0, Aijp[0], nb * (size_t)nb, Bijp[0], nb * (size_t)nb, 1.0,
                        &Vp[P][Q], naux);
            }
        }
    } else {
        V->scale(2.0);
    }
    V->hermitivitize();
    psio_->write_entry(unit_c_, "W", (char*)Vp[0], sizeof(double) * naux * naux);
}

void DFJKGrad::build_AB_x_terms()
{
    auto naux = auxiliary_->nbf();

    std::map<std::string, SharedMatrix> densities;

    if (do_J_) {
        auto d = std::make_shared<Vector>("d", naux);
        auto dp = d->pointer();
        psio_->read_entry(unit_c_, "c", (char*) dp, sizeof(double) * naux);
        auto D = std::make_shared<Matrix>("D", naux, naux);
        auto Dp = D->pointer();
        C_DGER(naux, naux, 1, dp, 1, dp, 1, Dp[0], naux);
        densities["Coulomb"] = D;
    }
    if (do_K_) {
        auto V = std::make_shared<Matrix>("V", naux, naux);
        auto Vp = V->pointer();
        psio_->read_entry(unit_c_, "V", (char*) Vp[0], sizeof(double) * naux * naux);
        densities["Exchange"] = V;
    }
    if (do_wK_) {
        auto W = std::make_shared<Matrix>("W", naux, naux);
        auto Wp = W->pointer();
        psio_->read_entry(unit_c_, "W", (char*) Wp[0], sizeof(double) * naux * naux);
        densities["Exchange,LR"] = W;
    }

    auto results = mints_->metric_grad(densities, "DF_BASIS_SCF");

    for (const auto& kv: results) {
        gradients_[kv.first] = kv.second;
    }
}
void DFJKGrad::build_Amn_x_terms() {
    // => Sizing <= //

    int natom = primary_->molecule()->natom();
    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    // => Integrals <= //

    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(df_ints_num_threads_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> omega_eri(df_ints_num_threads_);
    eri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri(1));
    if (do_wK_) {
        omega_eri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_, 1));
    }
    for (int t = 1; t < df_ints_num_threads_; t++) {
        eri[t] = std::shared_ptr<TwoBodyAOInt>(eri.front()->clone());
        if (do_wK_) {
            omega_eri[t] = std::shared_ptr<TwoBodyAOInt>(omega_eri.front()->clone());
        }
    }

    const std::vector<std::pair<int, int>>& shell_pairs = eri[0]->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //

    int max_rows;
    if (do_K_ || do_wK_) {
        int maxP = auxiliary_->max_function_per_shell();
        size_t row_cost = 0L;
        row_cost += nso * (size_t)nso;
        if (do_wK_) {
            row_cost += nso * (size_t)nso;
        }
        row_cost += nso * (size_t)na;
        row_cost += na * (size_t)na;
        size_t rows = memory_ / row_cost;
        rows = (rows > naux ? naux : rows);
        rows = (rows < maxP ? maxP : rows);
        max_rows = (int)rows;
    } else {
        max_rows = auxiliary_->nshell();
    }

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

    SharedVector d;
    double* dp;

    if (do_J_) {
        d = std::make_shared<Vector>("d", naux);
        dp = d->pointer();
        psio_->read_entry(unit_c_, "c", (char*)dp, sizeof(double) * naux);
    }

    SharedMatrix Kmn;
    SharedMatrix wKmn;
    SharedMatrix Ami;
    SharedMatrix Aij;

    double** Kmnp;
    double** wKmnp;
    double** Amip;
    double** Aijp;

    if (do_K_ || do_wK_) {
        Kmn = std::make_shared<Matrix>("Kmn", max_rows, nso * (size_t)nso);
        Ami = std::make_shared<Matrix>("Ami", max_rows, nso * (size_t)na);
        Aij = std::make_shared<Matrix>("Aij", max_rows, na * (size_t)na);
        Kmnp = Kmn->pointer();
        Amip = Ami->pointer();
        Aijp = Aij->pointer();
    }
    if (do_wK_) {
        wKmn = std::make_shared<Matrix>("wKmn", max_rows, nso * (size_t)nso);
        wKmnp = wKmn->pointer();
    }

    double** Dtp = Dt_->pointer();
    double** Cap = Ca_->pointer();
    double** Cbp = Cb_->pointer();

    psio_address next_Aija = PSIO_ZERO;
    psio_address next_Aijb = PSIO_ZERO;
    psio_address next_Awija = PSIO_ZERO;
    psio_address next_Awijb = PSIO_ZERO;

    // => Temporary Gradients <= //

    std::vector<SharedMatrix> Jtemps;
    std::vector<SharedMatrix> Ktemps;
    std::vector<SharedMatrix> wKtemps;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        if (do_J_) {
            Jtemps.push_back(std::make_shared<Matrix>("Jtemp", natom, 3));
        }
        if (do_K_) {
            Ktemps.push_back(std::make_shared<Matrix>("Ktemp", natom, 3));
        }
        if (do_wK_) {
            wKtemps.push_back(std::make_shared<Matrix>("wKtemp", natom, 3));
        }
    }

    // => R/U doubling factor <= //

    double factor = (restricted ? 2.0 : 1.0);

    // => Figure out required transforms <= //

    // unit, disk buffer name, nmo_size, psio_address, output_buffer
    std::vector<std::tuple<size_t, std::string, double**, size_t, psio_address*, double**>> transforms;
    if (do_K_ || do_wK_) {
        transforms.push_back(std::make_tuple(unit_a_, "(A|ij)", Cap, na, &next_Aija, Kmnp));
        // skip if there are no beta electrons
        if (!restricted && nb > 0) {
            transforms.push_back(std::make_tuple(unit_b_, "(A|ij)", Cbp, nb, &next_Aijb, Kmnp));
        }
    }
    if (do_wK_) {
        transforms.push_back(std::make_tuple(unit_a_, "(A|w|ij)", Cap, na, &next_Awija, wKmnp));
        // skip if there are no beta electrons        
        if (!restricted && nb > 0) {
            transforms.push_back(std::make_tuple(unit_b_, "(A|w|ij)", Cbp, nb, &next_Awijb, wKmnp));
        }
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

        // => J_mn^A <= //

        if (do_K_ || do_wK_) Kmn->zero();
        if (do_wK_) wKmn->zero();
        for (const auto& trans : transforms) {
            // > Unpack transform < //
            size_t unit = std::get<0>(trans);
            std::string buffer = std::get<1>(trans);
            double** Cp = std::get<2>(trans);
            size_t nmo = std::get<3>(trans);
            psio_address* address = std::get<4>(trans);
            double** retp = std::get<5>(trans);

            size_t nmo2 = nmo * nmo;

            // > Stripe < //
            psio_->read(unit, buffer.c_str(), (char*)Aijp[0], sizeof(double) * np * nmo2, *address, address);

            // > (A|ij) C_mi -> (A|mj) < //
#pragma omp parallel for
            for (int P = 0; P < np; P++) {
                C_DGEMM('N', 'N', nso, nmo, nmo, 1.0, Cp[0], nmo, &Aijp[0][P * nmo2], nmo, 0.0, Amip[P], na);
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N', 'T', np * (size_t)nso, nso, nmo, factor, Amip[0], na, Cp[0], nmo, 1.0, retp[0], nso);
        }

        // > Integrals < //
        int nthread_df = df_ints_num_threads_;
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
            double** grad_Kp;
            double** grad_wKp;

            if (do_J_) {
                grad_Jp = Jtemps[thread]->pointer();
            }
            if (do_K_) {
                grad_Kp = Ktemps[thread]->pointer();
            }
            if (do_wK_) {
                grad_wKp = wKtemps[thread]->pointer();
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

                        //  K^x = (A|pq)^x (A|pq)
                        if (do_K_) {
                            double Kval = 1.0 * perm * Kmnp[p + oP][(m + oM) * nso + (n + oN)];
                            grad_Kp[aP][0] += Kval * (*Px);
                            grad_Kp[aP][1] += Kval * (*Py);
                            grad_Kp[aP][2] += Kval * (*Pz);
                            grad_Kp[aM][0] += Kval * (*Mx);
                            grad_Kp[aM][1] += Kval * (*My);
                            grad_Kp[aM][2] += Kval * (*Mz);
                            grad_Kp[aN][0] += Kval * (*Nx);
                            grad_Kp[aN][1] += Kval * (*Ny);
                            grad_Kp[aN][2] += Kval * (*Nz);
                        }

                        // wK^x = 0.5 * (A|pq)^x (A|w|pq)
                        if (do_wK_) {
                            double wKval = 0.5 * perm * wKmnp[p + oP][(m + oM) * nso + (n + oN)];
                            grad_wKp[aP][0] += wKval * (*Px);
                            grad_wKp[aP][1] += wKval * (*Py);
                            grad_wKp[aP][2] += wKval * (*Pz);
                            grad_wKp[aM][0] += wKval * (*Mx);
                            grad_wKp[aM][1] += wKval * (*My);
                            grad_wKp[aM][2] += wKval * (*Mz);
                            grad_wKp[aN][0] += wKval * (*Nx);
                            grad_wKp[aN][1] += wKval * (*Ny);
                            grad_wKp[aN][2] += wKval * (*Nz);
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

            //  wK^x = 0.5 * (A|w|pq)^x (A|pq)
            if (do_wK_) {
                omega_eri[thread]->compute_shell_deriv1(P, 0, M, N);
                const double* buffer = omega_eri[thread]->buffer();

                const auto buffers = omega_eri[thread]->buffers();
                const double* Px = buffers[0];
                const double* Py = buffers[1];
                const double* Pz = buffers[2];
                const double* Mx = buffers[3];
                const double* My = buffers[4];
                const double* Mz = buffers[5];
                const double* Nx = buffers[6];
                const double* Ny = buffers[7];
                const double* Nz = buffers[8];

                for (int p = 0; p < nP; p++) {
                    for (int m = 0; m < nM; m++) {
                        for (int n = 0; n < nN; n++) {
                            double wKval = 0.5 * perm * Kmnp[p + oP][(m + oM) * nso + (n + oN)];
                            grad_wKp[aP][0] += wKval * (*Px);
                            grad_wKp[aP][1] += wKval * (*Py);
                            grad_wKp[aP][2] += wKval * (*Pz);
                            grad_wKp[aM][0] += wKval * (*Mx);
                            grad_wKp[aM][1] += wKval * (*My);
                            grad_wKp[aM][2] += wKval * (*Mz);
                            grad_wKp[aN][0] += wKval * (*Nx);
                            grad_wKp[aN][1] += wKval * (*Ny);
                            grad_wKp[aN][2] += wKval * (*Nz);
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
    }

    // => Temporary Gradient Reduction <= //

    // gradients_["Coulomb"]->zero();
    // gradients_["Exchange"]->zero();
    // gradients_["Exchange,LR"]->zero();

    if (do_J_) {
        for (int t = 0; t < df_ints_num_threads_; t++) {
            gradients_["Coulomb"]->add(Jtemps[t]);
        }
    }
    if (do_K_) {
        for (int t = 0; t < df_ints_num_threads_; t++) {
            gradients_["Exchange"]->add(Ktemps[t]);
        }
    }
    if (do_wK_) {
        for (int t = 0; t < df_ints_num_threads_; t++) {
            gradients_["Exchange,LR"]->add(wKtemps[t]);
        }
    }

    // gradients_["Coulomb"]->print();
    // gradients_["Exchange"]->print();
    // gradients_["Exchange,LR"]->print();
}

void DFJKGrad::compute_hessian() {
    // clang-format off
    /*
     * If we define Minv as the inverse metric matrix, and use the identity
     *
     *    d Minv          d M
     *    ------ = -Minv ----- Minv
     *      dx             dx
     *
     * then we get
     *
     *    d (mn|A) Minv[A][B] (B|rs)           x                                          x
     *    -------------------------- = 2 (mn|A)  Minv[A][B] (B|rs)  -  (mn|A) Minv[A][B] M[B][C] Minv[C][D] (D|rs)
     *                dx
     *
     * whose 2 terms we call term1 and term2, respectively.  Indices {m,n,r,s} refer to AOs,
     * while {A,B,C,D...} are aux basis indices.  The superscripts are just a shorthand for derivatives.  Expanding to get second derivatives, we have
     *
     *    d term1           xy                               x             y                                    x                 y
     *    ------- = 2 (mn|A)   Minv[A][B] (B|rs)  -  2 (mn|A)  Minv[A][B] M[B][C] Minv[C][D] (D|rs)  +  2 (mn|A) Minv[A][B] (B|rs)
     *       dy
     *
     *    d term2            y             x                                                 y                  x                                                 xy
     *    ------- = -2 (mn|A)  Minv[A][B] M[B][C] Minv[C][D] (D|rs)  +  2 (mn|A) Minv[A][B] M[B][C] Minv[C][D] M[D][E] Minv[E][F] (F|rs)  -  2 (mn|A) Minv[A][B] M[B][C] Minv[C][D] (D|rs)
     *       dy
     *
     * Note that the second term from term1 and the first from term2 are the same, leaving us with 5 terms to implement.  The code below is a first attempt at
     * this, and uses intermediates that were the first thing that came to mind.  The code needs to be adapted to run out of core (c.f. DFJKGrad::build_gradient())
     * and should call routines like build_Amn_terms() to generate intermediates, rather than re-coding that stuff here.  We also need to add UHF capabilities.
     *
     * Andy Simmonett (07/16)
     *
     */
     // clang-format on

    // => Set up hessians <= //
    int natom = primary_->molecule()->natom();
    hessians_.clear();
    double** JHessp = nullptr;
    double** KHessp = nullptr;
    if (do_J_) {
        hessians_["Coulomb"] = std::make_shared<Matrix>("Coulomb Hessian", 3 * natom, 3 * natom);
        JHessp = hessians_["Coulomb"]->pointer();
    }
    if (do_K_) {
        hessians_["Exchange"] = std::make_shared<Matrix>("Exchange Hessian", 3 * natom, 3 * natom);
        KHessp = hessians_["Exchange"]->pointer();
    }
    if (do_wK_) {
        hessians_["Exchange,LR"] = std::make_shared<Matrix>("Exchange,LR Hessian", 3 * natom, 3 * natom);
    }

    bool same_ab = (Ca_ == Cb_) ? true : false;

    std::shared_ptr<Molecule> mol = primary_->molecule();

    int np = auxiliary_->nbf();
    int nso = primary_->nbf();
    int nauxshell = auxiliary_->nshell();
    int nshell = primary_->nshell();
    int natoms = mol->natom();

    double** Dtp = Dt_->pointer();
    double** Cap = Ca_->pointer();
    double** Cbp = Cb_->pointer();

    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];
    auto metric = std::make_shared<FittingMetric>(auxiliary_, true);
    metric->form_full_eig_inverse(condition_);
    SharedMatrix PQ = metric->get_metric();
    double** PQp = PQ->pointer();

    auto c = std::make_shared<Vector>("c[A] = (mn|A) D[m][n]", np);
    double *cp = c->pointer();
    auto dc = std::make_shared<Matrix>("dc[x][A] = (mn|A)^x D[m][n]",  3*natoms, np);
    double **dcp = dc->pointer();

    auto dAa_ij = std::make_shared<Matrix>("dAij[x][A,i,j] = (mn|A)^x C[m][i] C[n][j]",  3*natoms, np*na*na);
    double **dAa_ijp = dAa_ij->pointer();
    auto dAb_ij = std::make_shared<Matrix>("dAij[x][A,i,j] = (mn|A)^x C[m][i] C[n][j]",  3*natoms, np*nb*nb);
    double **dAb_ijp = dAb_ij->pointer();

    auto d = std::make_shared<Vector>("d[A] = Minv[A][B] C[B]", np);
    double *dp = d->pointer();
    auto dd = std::make_shared<Matrix>("dd[x][B] = dc[x][A] Minv[A][B]", 3*natoms, np);
    double **ddp = dd->pointer();
    auto de = std::make_shared<Matrix>("de[x][A] = (A|B)^x d[B] ", 3*natoms, np);
    double **dep = de->pointer();
    auto dea_ij = std::make_shared<Matrix>("deij[x][A,i,j] = (A|B)^x Bij[B,i,j]", 3*natoms, np*na*na);
    double **dea_ijp = dea_ij->pointer();
    auto deb_ij = std::make_shared<Matrix>("deij[x][A,i,j] = (A|B)^x Bij[B,i,j]", 3*natoms, np*nb*nb);
    double **deb_ijp = deb_ij->pointer();

    // Build some integral factories
    auto Pmnfactory = std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_);
    auto PQfactory = std::make_shared<IntegralFactory>(auxiliary_, BasisSet::zero_ao_basis_set(), auxiliary_,
                                                       BasisSet::zero_ao_basis_set());
    std::shared_ptr<TwoBodyAOInt> Pmnint(Pmnfactory->eri(2));
    std::shared_ptr<TwoBodyAOInt> PQint(PQfactory->eri(2));
    auto Amn = std::make_shared<Matrix>("(A|mn)", np, nso*nso);
    auto Aa_mi = std::make_shared<Matrix>("(A|mi)", np, nso*na);
    auto Aa_ij = std::make_shared<Matrix>("(A|ij)", np, na*na);
    auto Ba_ij = std::make_shared<Matrix>("Minv[B][A] (A|ij)", np, na*na);
    auto Ba_im = std::make_shared<Matrix>("Minv[B][A] (A|im)", np, nso*na);
    auto Ba_mn = std::make_shared<Matrix>("Minv[B][A] (A|mn)", np, nso*nso);
    auto Da_PQ = std::make_shared<Matrix>("B(P|ij) B(Q|ij)", np, np);
    double **Amnp = Amn->pointer();
    double **Aa_mip = Aa_mi->pointer();
    double **Aa_ijp = Aa_ij->pointer();
    double **Ba_ijp = Ba_ij->pointer();
    double **Ba_imp = Ba_im->pointer();
    double **Ba_mnp = Ba_mn->pointer();
    double **Da_PQp = Da_PQ->pointer();

    auto Ab_mi = std::make_shared<Matrix>("(A|mi)", np, nso*nb);
    auto Ab_ij = std::make_shared<Matrix>("(A|ij)", np, nb*nb);
    auto Bb_ij = std::make_shared<Matrix>("Minv[B][A] (A|ij)", np, nb*nb);
    auto Bb_im = std::make_shared<Matrix>("Minv[B][A] (A|im)", np, nso*nb);
    auto Bb_mn = std::make_shared<Matrix>("Minv[B][A] (A|mn)", np, nso*nso);
    auto Db_PQ = std::make_shared<Matrix>("B(P|ij) B(Q|ij)", np, np);
    double **Ab_mip = Ab_mi->pointer();
    double **Ab_ijp = Ab_ij->pointer();
    double **Bb_ijp = Bb_ij->pointer();
    double **Bb_imp = Bb_im->pointer();
    double **Bb_mnp = Bb_mn->pointer();
    double **Db_PQp = Db_PQ->pointer();

    for (int P = 0; P < nauxshell; ++P) {
        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();
        for (int M = 0; M < nshell; ++M) {
            int nM = primary_->shell(M).nfunction();
            int oM = primary_->shell(M).function_index();
            for (int N = 0; N < nshell; ++N) {
                int nN = primary_->shell(N).nfunction();
                int oN = primary_->shell(N).function_index();

                Pmnint->compute_shell(P, 0, M, N);
                const double* buffer = Pmnint->buffer();

                for (int p = oP; p < oP+nP; p++) {
                    for (int m = oM; m < oM+nM; m++) {
                        for (int n = oN; n < oN+nN; n++) {
                            Amnp[p][m*nso+n] = (*buffer++);
                        }
                    }
                }
            }
        }
    }
    // First alpha
    // (A|mj) = (A|mn) C[n][j]
    C_DGEMM('N','N',np*(size_t)nso,na,nso,1.0,Amnp[0],nso,Cap[0],na,0.0,Aa_mip[0],na);
    // (A|ij) = (A|mj) C[m][i]
    #pragma omp parallel for
    for (int p = 0; p < np; p++) {
        C_DGEMM('T','N',na,na,nso,1.0,Aa_mip[p],na,Cap[0],na,0.0,&Aa_ijp[0][p * (size_t) na * na],na);
    }
    // Beta
    if (!same_ab){
        // (A|mj) = (A|mn) C[n][j]
        C_DGEMM('N','N',np*(size_t)nso,nb,nso,1.0,Amnp[0],nso,Cbp[0],nb,0.0,Ab_mip[0],nb);
        // (A|ij) = (A|mj) C[m][i]
        #pragma omp parallel for
        for (int p = 0; p < np; p++) {
            C_DGEMM('T','N',nb,nb,nso,1.0,Ab_mip[p],nb,Cbp[0],nb,0.0,&Ab_ijp[0][p * (size_t) nb * nb],nb);
        }
    }
    // c[A] = (A|mn) D[m][n]
    C_DGEMV('N', np, nso*(size_t)nso, 1.0, Amnp[0], nso*(size_t)nso, Dtp[0], 1, 0.0, cp, 1);
    // (A|mj) = (A|mn) C[n][j]
    C_DGEMM('N','N',np*(size_t)nso,na,nso,1.0,Amnp[0],nso,Cap[0],na,0.0,Aa_mip[0],na);
    // (A|ij) = (A|mj) C[m][i]
    #pragma omp parallel for
    for (int p = 0; p < np; p++) {
        C_DGEMM('T','N',na,na,nso,1.0,Aa_mip[p],na,Cap[0],na,0.0,&Aa_ijp[0][p * (size_t) na * na],na);
    }

    // d[A] = Minv[A][B] c[B]
    C_DGEMV('n', np, np, 1.0, PQp[0], np, cp, 1, 0.0, dp, 1);

    // Alpha
    // B[B][i,j] = Minv[A][B] (A|ij)
    C_DGEMM('n','n', np, na*na, np, 1.0, PQp[0], np, Aa_ijp[0], na*na, 0.0, Ba_ijp[0], na*na);
    // B[B][i,n] = B[B][i,j] C[n][j]
    C_DGEMM('N', 'T', np*(size_t)na, nso, na, 1.0, Ba_ijp[0], na, Cap[0], na, 0.0, Ba_imp[0], nso);
    // B[B][m,n] = C[m][i] B[B][i,n]
    #pragma omp parallel for
    for (int p = 0; p < np; p++) {
        C_DGEMM('n', 'n', nso, nso, na, 1.0, Cap[0], na, Ba_imp[p], nso, 0.0, Ba_mnp[p], nso);
    }
    // D[A][B] = B[A][ij] B[B][ij]
    C_DGEMM('n','t', np, np, na*na, 1.0, Ba_ijp[0], na*na, Ba_ijp[0], na*na, 0.0, Da_PQp[0], np);

    // Beta
    if(!same_ab){
        // B[B][i,j] = Minv[A][B] (A|ij)
        C_DGEMM('n','n', np, nb*nb, np, 1.0, PQp[0], np, Ab_ijp[0], nb*nb, 0.0, Bb_ijp[0], nb*nb);
        // B[B][i,n] = B[B][i,j] C[n][j]
        C_DGEMM('N', 'T', np*(size_t)nb, nso, nb, 1.0, Bb_ijp[0], nb, Cbp[0], nb, 0.0, Bb_imp[0], nso);
        // B[B][m,n] = C[m][i] B[B][i,n]
        #pragma omp parallel for
        for (int p = 0; p < np; p++) {
            C_DGEMM('n', 'n', nso, nso, nb, 1.0, Cbp[0], nb, Bb_imp[p], nso, 0.0, Bb_mnp[p], nso);
        }
        // D[A][B] = B[A][ij] B[B][ij]
        C_DGEMM('n','t', np, np, nb*nb, 1.0, Bb_ijp[0], nb*nb, Bb_ijp[0], nb*nb, 0.0, Db_PQp[0], np);
    }

    int maxp = auxiliary_->max_function_per_shell();
    int maxm = primary_->max_function_per_shell();
    auto Ta = std::make_shared<Matrix>("Ta", maxp, maxm*na);
    double **Tap = Ta->pointer();
    auto Tb = std::make_shared<Matrix>("Tb", maxp, maxm*nb);
    double **Tbp = Tb->pointer();

    for (int P = 0; P < nauxshell; ++P) {
        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();
        int Pcenter = auxiliary_->shell(P).ncenter();
        int Pncart = auxiliary_->shell(P).ncartesian();
        int Px = 3 * Pcenter + 0;
        int Py = 3 * Pcenter + 1;
        int Pz = 3 * Pcenter + 2;
        for (int M = 0; M < nshell; ++M) {
            int nM = primary_->shell(M).nfunction();
            int oM = primary_->shell(M).function_index();
            int Mcenter = primary_->shell(M).ncenter();
            int Mncart = primary_->shell(M).ncartesian();
            int mx = 3 * Mcenter + 0;
            int my = 3 * Mcenter + 1;
            int mz = 3 * Mcenter + 2;
            for (int N = 0; N < nshell; ++N) {
                int nN = primary_->shell(N).nfunction();
                int oN = primary_->shell(N).function_index();
                int Ncenter = primary_->shell(N).ncenter();
                int Nncart = primary_->shell(N).ncartesian();
                int nx = 3 * Ncenter + 0;
                int ny = 3 * Ncenter + 1;
                int nz = 3 * Ncenter + 2;

                Pmnint->compute_shell_deriv1(P, 0, M, N);
                const double* buffer = Pmnint->buffer();
                const auto& buffers = Pmnint->buffers();
                const double* PxBuf = buffers[0];
                const double* PyBuf = buffers[1];
                const double* PzBuf = buffers[2];
                const double* mxBuf = buffers[3];
                const double* myBuf = buffers[4];
                const double* mzBuf = buffers[5];
                const double* nxBuf = buffers[6];
                const double* nyBuf = buffers[7];
                const double* nzBuf = buffers[8];

                size_t delta = 0L;
                // Terms for J intermediates
                // dc[x][A] = D[m][n] (A|mn)^x
                for (int p = oP; p < oP + nP; p++) {
                    for (int m = oM; m < oM + nM; m++) {
                        for (int n = oN; n < oN + nN; n++) {
                            double Cpmn = Dtp[m][n];
                            dcp[Px][p] += Cpmn * PxBuf[delta];
                            dcp[Py][p] += Cpmn * PyBuf[delta];
                            dcp[Pz][p] += Cpmn * PzBuf[delta];
                            dcp[mx][p] += Cpmn * mxBuf[delta];
                            dcp[my][p] += Cpmn * myBuf[delta];
                            dcp[mz][p] += Cpmn * mzBuf[delta];
                            dcp[nx][p] += Cpmn * nxBuf[delta];
                            dcp[ny][p] += Cpmn * nyBuf[delta];
                            dcp[nz][p] += Cpmn * nzBuf[delta];
                            ++delta;
                        }
                    }
                }
                // Terms for K intermediates
                // dAij[x][p,i,j] <- (p|mn)^x C[m][i] C[n][j]
                //
                // implemented as
                //
                // T[p][m,j] <- (p|mn) C[n][j]
                // dAij[x][p,i,j] <- C[m][i] T[p][m,j]
                if(do_K_) {
                    // Alpha
                    C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, const_cast<double*>(PxBuf), nN, Cap[oN], na, 0.0, Tap[0], na);
#pragma omp parallel for
                    for(int p = 0; p < nP; ++p)
                        C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tap[0]+p*(nM*na), na, 1.0, &dAa_ijp[Px][(p+oP)*na*na], na);
                    C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, const_cast<double*>(PyBuf), nN, Cap[oN], na, 0.0, Tap[0], na);
#pragma omp parallel for
                    for(int p = 0; p < nP; ++p)
                        C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tap[0]+p*(nM*na), na, 1.0, &dAa_ijp[Py][(p+oP)*na*na], na);
                    C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, const_cast<double*>(PzBuf), nN, Cap[oN], na, 0.0, Tap[0], na);
#pragma omp parallel for
                    for(int p = 0; p < nP; ++p)
                        C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tap[0]+p*(nM*na), na, 1.0, &dAa_ijp[Pz][(p+oP)*na*na], na);
                    C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, const_cast<double*>(mxBuf), nN, Cap[oN], na, 0.0, Tap[0], na);
#pragma omp parallel for
                    for(int p = 0; p < nP; ++p)
                        C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tap[0]+p*(nM*na), na, 1.0, &dAa_ijp[mx][(p+oP)*na*na], na);
                    C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, const_cast<double*>(myBuf), nN, Cap[oN], na, 0.0, Tap[0], na);
#pragma omp parallel for
                    for(int p = 0; p < nP; ++p)
                        C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tap[0]+p*(nM*na), na, 1.0, &dAa_ijp[my][(p+oP)*na*na], na);
                    C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, const_cast<double*>(mzBuf), nN, Cap[oN], na, 0.0, Tap[0], na);
#pragma omp parallel for
                    for(int p = 0; p < nP; ++p)
                        C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tap[0]+p*(nM*na), na, 1.0, &dAa_ijp[mz][(p+oP)*na*na], na);
                    C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, const_cast<double*>(nxBuf), nN, Cap[oN], na, 0.0, Tap[0], na);
#pragma omp parallel for
                    for(int p = 0; p < nP; ++p)
                        C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tap[0]+p*(nM*na), na, 1.0, &dAa_ijp[nx][(p+oP)*na*na], na);
                    C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, const_cast<double*>(nyBuf), nN, Cap[oN], na, 0.0, Tap[0], na);
#pragma omp parallel for
                    for(int p = 0; p < nP; ++p)
                        C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tap[0]+p*(nM*na), na, 1.0, &dAa_ijp[ny][(p+oP)*na*na], na);
                    C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, const_cast<double*>(nzBuf), nN, Cap[oN], na, 0.0, Tap[0], na);
#pragma omp parallel for
                    for(int p = 0; p < nP; ++p)
                        C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tap[0]+p*(nM*na), na, 1.0, &dAa_ijp[nz][(p+oP)*na*na], na);

                    // Beta
                    if (!same_ab){
                        C_DGEMM('n', 'n', nP*nM, nb, nN, 1.0, const_cast<double*>(PxBuf), nN, Cbp[oN], nb, 0.0, Tbp[0], nb);
#pragma omp parallel for
                        for(int p = 0; p < nP; ++p)
                            C_DGEMM('t', 'n', nb, nb, nM, 1.0, Cbp[oM], nb, Tbp[0]+p*(nM*nb), nb, 1.0, &dAb_ijp[Px][(p+oP)*nb*nb], nb);
                        C_DGEMM('n', 'n', nP*nM, nb, nN, 1.0, const_cast<double*>(PyBuf), nN, Cbp[oN], nb, 0.0, Tbp[0], nb);
#pragma omp parallel for
                        for(int p = 0; p < nP; ++p)
                            C_DGEMM('t', 'n', nb, nb, nM, 1.0, Cbp[oM], nb, Tbp[0]+p*(nM*nb), nb, 1.0, &dAb_ijp[Py][(p+oP)*nb*nb], nb);
                        C_DGEMM('n', 'n', nP*nM, nb, nN, 1.0, const_cast<double*>(PzBuf), nN, Cbp[oN], nb, 0.0, Tbp[0], nb);
#pragma omp parallel for
                        for(int p = 0; p < nP; ++p)
                            C_DGEMM('t', 'n', nb, nb, nM, 1.0, Cbp[oM], nb, Tbp[0]+p*(nM*nb), nb, 1.0, &dAb_ijp[Pz][(p+oP)*nb*nb], nb);
                        C_DGEMM('n', 'n', nP*nM, nb, nN, 1.0, const_cast<double*>(mxBuf), nN, Cbp[oN], nb, 0.0, Tbp[0], nb);
#pragma omp parallel for
                        for(int p = 0; p < nP; ++p)
                            C_DGEMM('t', 'n', nb, nb, nM, 1.0, Cbp[oM], nb, Tbp[0]+p*(nM*nb), nb, 1.0, &dAb_ijp[mx][(p+oP)*nb*nb], nb);
                        C_DGEMM('n', 'n', nP*nM, nb, nN, 1.0, const_cast<double*>(myBuf), nN, Cbp[oN], nb, 0.0, Tbp[0], nb);
#pragma omp parallel for
                        for(int p = 0; p < nP; ++p)
                            C_DGEMM('t', 'n', nb, nb, nM, 1.0, Cbp[oM], nb, Tbp[0]+p*(nM*nb), nb, 1.0, &dAb_ijp[my][(p+oP)*nb*nb], nb);
                        C_DGEMM('n', 'n', nP*nM, nb, nN, 1.0, const_cast<double*>(mzBuf), nN, Cbp[oN], nb, 0.0, Tbp[0], nb);
#pragma omp parallel for
                        for(int p = 0; p < nP; ++p)
                            C_DGEMM('t', 'n', nb, nb, nM, 1.0, Cbp[oM], nb, Tbp[0]+p*(nM*nb), nb, 1.0, &dAb_ijp[mz][(p+oP)*nb*nb], nb);
                        C_DGEMM('n', 'n', nP*nM, nb, nN, 1.0, const_cast<double*>(nxBuf), nN, Cbp[oN], nb, 0.0, Tbp[0], nb);
#pragma omp parallel for
                        for(int p = 0; p < nP; ++p)
                            C_DGEMM('t', 'n', nb, nb, nM, 1.0, Cbp[oM], nb, Tbp[0]+p*(nM*nb), nb, 1.0, &dAb_ijp[nx][(p+oP)*nb*nb], nb);
                        C_DGEMM('n', 'n', nP*nM, nb, nN, 1.0, const_cast<double*>(nyBuf), nN, Cbp[oN], nb, 0.0, Tbp[0], nb);
#pragma omp parallel for
                        for(int p = 0; p < nP; ++p)
                            C_DGEMM('t', 'n', nb, nb, nM, 1.0, Cbp[oM], nb, Tbp[0]+p*(nM*nb), nb, 1.0, &dAb_ijp[ny][(p+oP)*nb*nb], nb);
                        C_DGEMM('n', 'n', nP*nM, nb, nN, 1.0, const_cast<double*>(nzBuf), nN, Cbp[oN], nb, 0.0, Tbp[0], nb);
#pragma omp parallel for
                        for(int p = 0; p < nP; ++p)
                            C_DGEMM('t', 'n', nb, nb, nM, 1.0, Cbp[oM], nb, Tbp[0]+p*(nM*nb), nb, 1.0, &dAb_ijp[nz][(p+oP)*nb*nb], nb);
                    }
                }
            }
        }
    }

    // dd[x][A] = dc[x][B] Minv[B][A]
    C_DGEMM('N', 'N', 3 * natoms, np, np, 1.0, dcp[0], np, PQp[0], np, 0.0, ddp[0], np);

    for (int P = 0; P < nauxshell; ++P) {
        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();
        int Pcenter = auxiliary_->shell(P).ncenter();
        int Pncart = auxiliary_->shell(P).ncartesian();
        int Px = 3 * Pcenter + 0;
        int Py = 3 * Pcenter + 1;
        int Pz = 3 * Pcenter + 2;
        for (int Q = 0; Q < nauxshell; ++Q) {
            int nQ = auxiliary_->shell(Q).nfunction();
            int oQ = auxiliary_->shell(Q).function_index();
            int Qcenter = auxiliary_->shell(Q).ncenter();
            int Qncart = auxiliary_->shell(Q).ncartesian();
            int Qx = 3 * Qcenter + 0;
            int Qy = 3 * Qcenter + 1;
            int Qz = 3 * Qcenter + 2;

            //size_t stride = static_cast<size_t>(Pncart) * Qncart;

            PQint->compute_shell_deriv1(P, 0, Q, 0);
            const auto& buffers = PQint->buffers();
            const double* Pxbuf = buffers[0];
            const double* Pybuf = buffers[1];
            const double* Pzbuf = buffers[2];
            const double* Qxbuf = buffers[3];
            const double* Qybuf = buffers[4];
            const double* Qzbuf = buffers[5];

            size_t delta = 0L;
            // J term intermediates
            // de[x][A] = (A|B)^x d[B]
            for (int p = oP; p < oP + nP; p++) {
                for (int q = oQ; q < oQ + nQ; q++) {
                    double dq = dp[q];
                    dep[Px][p] += dq * Pxbuf[delta];
                    dep[Py][p] += dq * Pybuf[delta];
                    dep[Pz][p] += dq * Pzbuf[delta];
                    dep[Qx][p] += dq * Qxbuf[delta];
                    dep[Qy][p] += dq * Qybuf[delta];
                    dep[Qz][p] += dq * Qzbuf[delta];
                    ++delta;
                }
            }
            // K term intermediates
            // deij[x][A,i,j] <- (A|B)^x Bij[B,i,j]
            if(do_K_){
                C_DGEMM('n', 'n', nP, na*na, nQ, 1.0, const_cast<double*>(Pxbuf), nQ, Ba_ijp[oQ], na*na, 1.0, &dea_ijp[Px][oP*na*na], na*na);
                C_DGEMM('n', 'n', nP, na*na, nQ, 1.0, const_cast<double*>(Pybuf), nQ, Ba_ijp[oQ], na*na, 1.0, &dea_ijp[Py][oP*na*na], na*na);
                C_DGEMM('n', 'n', nP, na*na, nQ, 1.0, const_cast<double*>(Pzbuf), nQ, Ba_ijp[oQ], na*na, 1.0, &dea_ijp[Pz][oP*na*na], na*na);
                C_DGEMM('n', 'n', nP, na*na, nQ, 1.0, const_cast<double*>(Qxbuf), nQ, Ba_ijp[oQ], na*na, 1.0, &dea_ijp[Qx][oP*na*na], na*na);
                C_DGEMM('n', 'n', nP, na*na, nQ, 1.0, const_cast<double*>(Qybuf), nQ, Ba_ijp[oQ], na*na, 1.0, &dea_ijp[Qy][oP*na*na], na*na);
                C_DGEMM('n', 'n', nP, na*na, nQ, 1.0, const_cast<double*>(Qzbuf), nQ, Ba_ijp[oQ], na*na, 1.0, &dea_ijp[Qz][oP*na*na], na*na);

                if (!same_ab){
                    C_DGEMM('n', 'n', nP, nb*nb, nQ, 1.0, const_cast<double*>(Pxbuf), nQ, Bb_ijp[oQ], nb*nb, 1.0, &deb_ijp[Px][oP*nb*nb], nb*nb);
                    C_DGEMM('n', 'n', nP, nb*nb, nQ, 1.0, const_cast<double*>(Pybuf), nQ, Bb_ijp[oQ], nb*nb, 1.0, &deb_ijp[Py][oP*nb*nb], nb*nb);
                    C_DGEMM('n', 'n', nP, nb*nb, nQ, 1.0, const_cast<double*>(Pzbuf), nQ, Bb_ijp[oQ], nb*nb, 1.0, &deb_ijp[Pz][oP*nb*nb], nb*nb);
                    C_DGEMM('n', 'n', nP, nb*nb, nQ, 1.0, const_cast<double*>(Qxbuf), nQ, Bb_ijp[oQ], nb*nb, 1.0, &deb_ijp[Qx][oP*nb*nb], nb*nb);
                    C_DGEMM('n', 'n', nP, nb*nb, nQ, 1.0, const_cast<double*>(Qybuf), nQ, Bb_ijp[oQ], nb*nb, 1.0, &deb_ijp[Qy][oP*nb*nb], nb*nb);
                    C_DGEMM('n', 'n', nP, nb*nb, nQ, 1.0, const_cast<double*>(Qzbuf), nQ, Bb_ijp[oQ], nb*nb, 1.0, &deb_ijp[Qz][oP*nb*nb], nb*nb);
                }
            }
        }
    }

    for (int P = 0; P < nauxshell; ++P) {
        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();
        int Pcenter = auxiliary_->shell(P).ncenter();
        int Pncart = auxiliary_->shell(P).ncartesian();
        int Px = 3 * Pcenter + 0;
        int Py = 3 * Pcenter + 1;
        int Pz = 3 * Pcenter + 2;
        for (int M = 0; M < nshell; ++M) {
            int nM = primary_->shell(M).nfunction();
            int oM = primary_->shell(M).function_index();
            int Mcenter = primary_->shell(M).ncenter();
            int Mncart = primary_->shell(M).ncartesian();
            int mx = 3 * Mcenter + 0;
            int my = 3 * Mcenter + 1;
            int mz = 3 * Mcenter + 2;
            for (int N = 0; N < nshell; ++N) {
                int nN = primary_->shell(N).nfunction();
                int oN = primary_->shell(N).function_index();
                int Ncenter = primary_->shell(N).ncenter();
                int Nncart = primary_->shell(N).ncartesian();
                int nx = 3 * Ncenter + 0;
                int ny = 3 * Ncenter + 1;
                int nz = 3 * Ncenter + 2;

                Pmnint->compute_shell_deriv2(P, 0, M, N);
                const auto& buffers = Pmnint->buffers();
                const double* PxPxBuf = buffers[0];
                const double* PxPyBuf = buffers[1];
                const double* PxPzBuf = buffers[2];
                const double* PxmxBuf = buffers[3];
                const double* PxmyBuf = buffers[4];
                const double* PxmzBuf = buffers[5];
                const double* PxnxBuf = buffers[6];
                const double* PxnyBuf = buffers[7];
                const double* PxnzBuf = buffers[8];
                const double* PyPyBuf = buffers[9];
                const double* PyPzBuf = buffers[10];
                const double* PymxBuf = buffers[11];
                const double* PymyBuf = buffers[12];
                const double* PymzBuf = buffers[13];
                const double* PynxBuf = buffers[14];
                const double* PynyBuf = buffers[15];
                const double* PynzBuf = buffers[16];
                const double* PzPzBuf = buffers[17];
                const double* PzmxBuf = buffers[18];
                const double* PzmyBuf = buffers[19];
                const double* PzmzBuf = buffers[20];
                const double* PznxBuf = buffers[21];
                const double* PznyBuf = buffers[22];
                const double* PznzBuf = buffers[23];
                const double* mxmxBuf = buffers[24];
                const double* mxmyBuf = buffers[25];
                const double* mxmzBuf = buffers[26];
                const double* mxnxBuf = buffers[27];
                const double* mxnyBuf = buffers[28];
                const double* mxnzBuf = buffers[29];
                const double* mymyBuf = buffers[30];
                const double* mymzBuf = buffers[31];
                const double* mynxBuf = buffers[32];
                const double* mynyBuf = buffers[33];
                const double* mynzBuf = buffers[34];
                const double* mzmzBuf = buffers[35];
                const double* mznxBuf = buffers[36];
                const double* mznyBuf = buffers[37];
                const double* mznzBuf = buffers[38];
                const double* nxnxBuf = buffers[39];
                const double* nxnyBuf = buffers[40];
                const double* nxnzBuf = buffers[41];
                const double* nynyBuf = buffers[42];
                const double* nynzBuf = buffers[43];
                const double* nznzBuf = buffers[44];

                double Pmscale = Pcenter == Mcenter ? 2.0 : 1.0;
                double Pnscale = Pcenter == Ncenter ? 2.0 : 1.0;
                double mnscale = Mcenter == Ncenter ? 2.0 : 1.0;

                double PxPx = 0.0, PxPy = 0.0, PxPz = 0.0, PyPy = 0.0, PyPz = 0.0, PzPz = 0.0;
                double mxmx = 0.0, mxmy = 0.0, mxmz = 0.0, mymy = 0.0, mymz = 0.0, mzmz = 0.0;
                double nxnx = 0.0, nxny = 0.0, nxnz = 0.0, nyny = 0.0, nynz = 0.0, nznz = 0.0;
                double Pxmx = 0.0, Pxmy = 0.0, Pxmz = 0.0, Pymx = 0.0, Pymy = 0.0, Pymz = 0.0, Pzmx = 0.0, Pzmy = 0.0,
                       Pzmz = 0.0;
                double Pxnx = 0.0, Pxny = 0.0, Pxnz = 0.0, Pynx = 0.0, Pyny = 0.0, Pynz = 0.0, Pznx = 0.0, Pzny = 0.0,
                       Pznz = 0.0;
                double mxnx = 0.0, mxny = 0.0, mxnz = 0.0, mynx = 0.0, myny = 0.0, mynz = 0.0, mznx = 0.0, mzny = 0.0,
                       mznz = 0.0;
                size_t delta = 0L;
                for (int p = oP; p < oP + nP; p++) {
                    for (int m = oM; m < oM + nM; m++) {
                        for (int n = oN; n < oN + nN; n++) {
                            double Cpmn = 2.0 * dp[p] * Dtp[m][n];
                            PxPx += Cpmn * PxPxBuf[delta];
                            PxPy += Cpmn * PxPyBuf[delta];
                            PxPz += Cpmn * PxPzBuf[delta];
                            Pxmx += Cpmn * PxmxBuf[delta];
                            Pxmy += Cpmn * PxmyBuf[delta];
                            Pxmz += Cpmn * PxmzBuf[delta];
                            Pxnx += Cpmn * PxnxBuf[delta];
                            Pxny += Cpmn * PxnyBuf[delta];
                            Pxnz += Cpmn * PxnzBuf[delta];
                            PyPy += Cpmn * PyPyBuf[delta];
                            PyPz += Cpmn * PyPzBuf[delta];
                            Pymx += Cpmn * PymxBuf[delta];
                            Pymy += Cpmn * PymyBuf[delta];
                            Pymz += Cpmn * PymzBuf[delta];
                            Pynx += Cpmn * PynxBuf[delta];
                            Pyny += Cpmn * PynyBuf[delta];
                            Pynz += Cpmn * PynzBuf[delta];
                            PzPz += Cpmn * PzPzBuf[delta];
                            Pzmx += Cpmn * PzmxBuf[delta];
                            Pzmy += Cpmn * PzmyBuf[delta];
                            Pzmz += Cpmn * PzmzBuf[delta];
                            Pznx += Cpmn * PznxBuf[delta];
                            Pzny += Cpmn * PznyBuf[delta];
                            Pznz += Cpmn * PznzBuf[delta];
                            mxmx += Cpmn * mxmxBuf[delta];
                            mxmy += Cpmn * mxmyBuf[delta];
                            mxmz += Cpmn * mxmzBuf[delta];
                            mxnx += Cpmn * mxnxBuf[delta];
                            mxny += Cpmn * mxnyBuf[delta];
                            mxnz += Cpmn * mxnzBuf[delta];
                            mymy += Cpmn * mymyBuf[delta];
                            mymz += Cpmn * mymzBuf[delta];
                            mynx += Cpmn * mynxBuf[delta];
                            myny += Cpmn * mynyBuf[delta];
                            mynz += Cpmn * mynzBuf[delta];
                            mzmz += Cpmn * mzmzBuf[delta];
                            mznx += Cpmn * mznxBuf[delta];
                            mzny += Cpmn * mznyBuf[delta];
                            mznz += Cpmn * mznzBuf[delta];
                            nxnx += Cpmn * nxnxBuf[delta];
                            nxny += Cpmn * nxnyBuf[delta];
                            nxnz += Cpmn * nxnzBuf[delta];
                            nyny += Cpmn * nynyBuf[delta];
                            nynz += Cpmn * nynzBuf[delta];
                            nznz += Cpmn * nznzBuf[delta];
                            ++delta;
                        }
                    }
                }
                JHessp[Px][Px] += PxPx;
                JHessp[Px][Py] += PxPy;
                JHessp[Px][Pz] += PxPz;
                JHessp[Px][mx] += Pmscale * Pxmx;
                JHessp[Px][my] += Pxmy;
                JHessp[Px][mz] += Pxmz;
                JHessp[Px][nx] += Pnscale * Pxnx;
                JHessp[Px][ny] += Pxny;
                JHessp[Px][nz] += Pxnz;
                JHessp[Py][Py] += PyPy;
                JHessp[Py][Pz] += PyPz;
                JHessp[Py][mx] += Pymx;
                JHessp[Py][my] += Pmscale * Pymy;
                JHessp[Py][mz] += Pymz;
                JHessp[Py][nx] += Pynx;
                JHessp[Py][ny] += Pnscale * Pyny;
                JHessp[Py][nz] += Pynz;
                JHessp[Pz][Pz] += PzPz;
                JHessp[Pz][mx] += Pzmx;
                JHessp[Pz][my] += Pzmy;
                JHessp[Pz][mz] += Pmscale * Pzmz;
                JHessp[Pz][nx] += Pznx;
                JHessp[Pz][ny] += Pzny;
                JHessp[Pz][nz] += Pnscale * Pznz;
                JHessp[mx][mx] += mxmx;
                JHessp[mx][my] += mxmy;
                JHessp[mx][mz] += mxmz;
                JHessp[mx][nx] += mnscale * mxnx;
                JHessp[mx][ny] += mxny;
                JHessp[mx][nz] += mxnz;
                JHessp[my][my] += mymy;
                JHessp[my][mz] += mymz;
                JHessp[my][nx] += mynx;
                JHessp[my][ny] += mnscale * myny;
                JHessp[my][nz] += mynz;
                JHessp[mz][mz] += mzmz;
                JHessp[mz][nx] += mznx;
                JHessp[mz][ny] += mzny;
                JHessp[mz][nz] += mnscale * mznz;
                JHessp[nx][nx] += nxnx;
                JHessp[nx][ny] += nxny;
                JHessp[nx][nz] += nxnz;
                JHessp[ny][ny] += nyny;
                JHessp[ny][nz] += nynz;
                JHessp[nz][nz] += nznz;

                if (do_K_) {
                    // K terms
                    // Loop through alpha/beta terms
                    std::vector<double**> Blist = {Ba_mnp};
                    if (!same_ab) Blist.push_back(Bb_mnp);
                    for (auto& B : Blist){
                        PxPx=0.0; PxPy=0.0; PxPz=0.0; PyPy=0.0; PyPz=0.0; PzPz=0.0;
                        mxmx=0.0; mxmy=0.0; mxmz=0.0; mymy=0.0; mymz=0.0; mzmz=0.0;
                        nxnx=0.0; nxny=0.0; nxnz=0.0; nyny=0.0; nynz=0.0; nznz=0.0;
                        Pxmx=0.0; Pxmy=0.0; Pxmz=0.0; Pymx=0.0; Pymy=0.0; Pymz=0.0; Pzmx=0.0; Pzmy=0.0; Pzmz=0.0;
                        Pxnx=0.0; Pxny=0.0; Pxnz=0.0; Pynx=0.0; Pyny=0.0; Pynz=0.0; Pznx=0.0; Pzny=0.0; Pznz=0.0;
                        mxnx=0.0; mxny=0.0; mxnz=0.0; mynx=0.0; myny=0.0; mynz=0.0; mznx=0.0; mzny=0.0; mznz=0.0;
                        delta = 0L;
                        for (int p = oP; p < oP+nP; p++) {
                            for (int m = oM; m < oM+nM; m++) {
                                for (int n = oN; n < oN+nN; n++) {
                                    double Cpmn = 2.0 * B[p][m*nso+n];
                                    PxPx += Cpmn * PxPxBuf[delta];
                                    PxPy += Cpmn * PxPyBuf[delta];
                                    PxPz += Cpmn * PxPzBuf[delta];
                                    Pxmx += Cpmn * PxmxBuf[delta];
                                    Pxmy += Cpmn * PxmyBuf[delta];
                                    Pxmz += Cpmn * PxmzBuf[delta];
                                    Pxnx += Cpmn * PxnxBuf[delta];
                                    Pxny += Cpmn * PxnyBuf[delta];
                                    Pxnz += Cpmn * PxnzBuf[delta];
                                    PyPy += Cpmn * PyPyBuf[delta];
                                    PyPz += Cpmn * PyPzBuf[delta];
                                    Pymx += Cpmn * PymxBuf[delta];
                                    Pymy += Cpmn * PymyBuf[delta];
                                    Pymz += Cpmn * PymzBuf[delta];
                                    Pynx += Cpmn * PynxBuf[delta];
                                    Pyny += Cpmn * PynyBuf[delta];
                                    Pynz += Cpmn * PynzBuf[delta];
                                    PzPz += Cpmn * PzPzBuf[delta];
                                    Pzmx += Cpmn * PzmxBuf[delta];
                                    Pzmy += Cpmn * PzmyBuf[delta];
                                    Pzmz += Cpmn * PzmzBuf[delta];
                                    Pznx += Cpmn * PznxBuf[delta];
                                    Pzny += Cpmn * PznyBuf[delta];
                                    Pznz += Cpmn * PznzBuf[delta];
                                    mxmx += Cpmn * mxmxBuf[delta];
                                    mxmy += Cpmn * mxmyBuf[delta];
                                    mxmz += Cpmn * mxmzBuf[delta];
                                    mxnx += Cpmn * mxnxBuf[delta];
                                    mxny += Cpmn * mxnyBuf[delta];
                                    mxnz += Cpmn * mxnzBuf[delta];
                                    mymy += Cpmn * mymyBuf[delta];
                                    mymz += Cpmn * mymzBuf[delta];
                                    mynx += Cpmn * mynxBuf[delta];
                                    myny += Cpmn * mynyBuf[delta];
                                    mynz += Cpmn * mynzBuf[delta];
                                    mzmz += Cpmn * mzmzBuf[delta];
                                    mznx += Cpmn * mznxBuf[delta];
                                    mzny += Cpmn * mznyBuf[delta];
                                    mznz += Cpmn * mznzBuf[delta];
                                    nxnx += Cpmn * nxnxBuf[delta];
                                    nxny += Cpmn * nxnyBuf[delta];
                                    nxnz += Cpmn * nxnzBuf[delta];
                                    nyny += Cpmn * nynyBuf[delta];
                                    nynz += Cpmn * nynzBuf[delta];
                                    nznz += Cpmn * nznzBuf[delta];
                                    ++delta;
                                }
                            }
                        }
                        KHessp[Px][Px] += PxPx;
                        KHessp[Px][Py] += PxPy;
                        KHessp[Px][Pz] += PxPz;
                        KHessp[Px][mx] += Pmscale*Pxmx;
                        KHessp[Px][my] += Pxmy;
                        KHessp[Px][mz] += Pxmz;
                        KHessp[Px][nx] += Pnscale*Pxnx;
                        KHessp[Px][ny] += Pxny;
                        KHessp[Px][nz] += Pxnz;
                        KHessp[Py][Py] += PyPy;
                        KHessp[Py][Pz] += PyPz;
                        KHessp[Py][mx] += Pymx;
                        KHessp[Py][my] += Pmscale*Pymy;
                        KHessp[Py][mz] += Pymz;
                        KHessp[Py][nx] += Pynx;
                        KHessp[Py][ny] += Pnscale*Pyny;
                        KHessp[Py][nz] += Pynz;
                        KHessp[Pz][Pz] += PzPz;
                        KHessp[Pz][mx] += Pzmx;
                        KHessp[Pz][my] += Pzmy;
                        KHessp[Pz][mz] += Pmscale*Pzmz;
                        KHessp[Pz][nx] += Pznx;
                        KHessp[Pz][ny] += Pzny;
                        KHessp[Pz][nz] += Pnscale*Pznz;
                        KHessp[mx][mx] += mxmx;
                        KHessp[mx][my] += mxmy;
                        KHessp[mx][mz] += mxmz;
                        KHessp[mx][nx] += mnscale*mxnx;
                        KHessp[mx][ny] += mxny;
                        KHessp[mx][nz] += mxnz;
                        KHessp[my][my] += mymy;
                        KHessp[my][mz] += mymz;
                        KHessp[my][nx] += mynx;
                        KHessp[my][ny] += mnscale*myny;
                        KHessp[my][nz] += mynz;
                        KHessp[mz][mz] += mzmz;
                        KHessp[mz][nx] += mznx;
                        KHessp[mz][ny] += mzny;
                        KHessp[mz][nz] += mnscale*mznz;
                        KHessp[nx][nx] += nxnx;
                        KHessp[nx][ny] += nxny;
                        KHessp[nx][nz] += nxnz;
                        KHessp[ny][ny] += nyny;
                        KHessp[ny][nz] += nynz;
                        KHessp[nz][nz] += nznz;
                    }
                }
            }
        }
    }

    for (int P = 0; P < nauxshell; ++P) {
        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();
        int Pcenter = auxiliary_->shell(P).ncenter();
        int Pncart = auxiliary_->shell(P).ncartesian();
        int Px = 3 * Pcenter + 0;
        int Py = 3 * Pcenter + 1;
        int Pz = 3 * Pcenter + 2;
        for (int Q = 0; Q < nauxshell; ++Q) {
            int nQ = auxiliary_->shell(Q).nfunction();
            int oQ = auxiliary_->shell(Q).function_index();
            int Qcenter = auxiliary_->shell(Q).ncenter();
            int Qncart = auxiliary_->shell(Q).ncartesian();
            int Qx = 3 * Qcenter + 0;
            int Qy = 3 * Qcenter + 1;
            int Qz = 3 * Qcenter + 2;

            PQint->compute_shell_deriv2(P, 0, Q, 0);
            const auto& buffers = PQint->buffers();
            const double* PxPxBuf = buffers[0];
            const double* PxPyBuf = buffers[1];
            const double* PxPzBuf = buffers[2];
            const double* PxQxBuf = buffers[3];
            const double* PxQyBuf = buffers[4];
            const double* PxQzBuf = buffers[5];
            const double* PyPyBuf = buffers[6];
            const double* PyPzBuf = buffers[7];
            const double* PyQxBuf = buffers[8];
            const double* PyQyBuf = buffers[9];
            const double* PyQzBuf = buffers[10];
            const double* PzPzBuf = buffers[11];
            const double* PzQxBuf = buffers[12];
            const double* PzQyBuf = buffers[13];
            const double* PzQzBuf = buffers[14];
            const double* QxQxBuf = buffers[15];
            const double* QxQyBuf = buffers[16];
            const double* QxQzBuf = buffers[17];
            const double* QyQyBuf = buffers[18];
            const double* QyQzBuf = buffers[19];
            const double* QzQzBuf = buffers[20];

            double PQscale = Pcenter == Qcenter ? 2.0 : 1.0;

            double PxPx = 0.0, PxPy = 0.0, PxPz = 0.0, PyPy = 0.0, PyPz = 0.0, PzPz = 0.0;
            double QxQx = 0.0, QxQy = 0.0, QxQz = 0.0, QyQy = 0.0, QyQz = 0.0, QzQz = 0.0;
            double PxQx = 0.0, PxQy = 0.0, PxQz = 0.0, PyQx = 0.0, PyQy = 0.0, PyQz = 0.0, PzQx = 0.0, PzQy = 0.0,
                   PzQz = 0.0;
            size_t delta = 0L;
            for (int p = oP; p < oP + nP; p++) {
                for (int q = oQ; q < oQ + nQ; q++) {
                    double dAdB = -dp[p] * dp[q];
                    PxPx += dAdB * PxPxBuf[delta];
                    PxPy += dAdB * PxPyBuf[delta];
                    PxPz += dAdB * PxPzBuf[delta];
                    PxQx += dAdB * PxQxBuf[delta];
                    PxQy += dAdB * PxQyBuf[delta];
                    PxQz += dAdB * PxQzBuf[delta];
                    PyPy += dAdB * PyPyBuf[delta];
                    PyPz += dAdB * PyPzBuf[delta];
                    PyQx += dAdB * PyQxBuf[delta];
                    PyQy += dAdB * PyQyBuf[delta];
                    PyQz += dAdB * PyQzBuf[delta];
                    PzPz += dAdB * PzPzBuf[delta];
                    PzQx += dAdB * PzQxBuf[delta];
                    PzQy += dAdB * PzQyBuf[delta];
                    PzQz += dAdB * PzQzBuf[delta];
                    QxQx += dAdB * QxQxBuf[delta];
                    QxQy += dAdB * QxQyBuf[delta];
                    QxQz += dAdB * QxQzBuf[delta];
                    QyQy += dAdB * QyQyBuf[delta];
                    QyQz += dAdB * QyQzBuf[delta];
                    QzQz += dAdB * QzQzBuf[delta];
                    ++delta;
                }
            }
            JHessp[Px][Px] += PxPx;
            JHessp[Px][Py] += PxPy;
            JHessp[Px][Pz] += PxPz;
            JHessp[Px][Qx] += PQscale * PxQx;
            JHessp[Px][Qy] += PxQy;
            JHessp[Px][Qz] += PxQz;
            JHessp[Py][Py] += PyPy;
            JHessp[Py][Pz] += PyPz;
            JHessp[Py][Qx] += PyQx;
            JHessp[Py][Qy] += PQscale * PyQy;
            JHessp[Py][Qz] += PyQz;
            JHessp[Pz][Pz] += PzPz;
            JHessp[Pz][Qx] += PzQx;
            JHessp[Pz][Qy] += PzQy;
            JHessp[Pz][Qz] += PQscale * PzQz;
            JHessp[Qx][Qx] += QxQx;
            JHessp[Qx][Qy] += QxQy;
            JHessp[Qx][Qz] += QxQz;
            JHessp[Qy][Qy] += QyQy;
            JHessp[Qy][Qz] += QyQz;
            JHessp[Qz][Qz] += QzQz;

            if (do_K_) {
                // K terms
                std::vector<double**> Dlist = {Da_PQp};
                if (!same_ab) Dlist.push_back(Db_PQp);
                for (auto& DPQp : Dlist){
                    PxPx=0.0; PxPy=0.0; PxPz=0.0; PyPy=0.0; PyPz=0.0; PzPz=0.0;
                    QxQx=0.0; QxQy=0.0; QxQz=0.0; QyQy=0.0; QyQz=0.0; QzQz=0.0;
                    PxQx=0.0; PxQy=0.0; PxQz=0.0; PyQx=0.0; PyQy=0.0; PyQz=0.0; PzQx=0.0; PzQy=0.0; PzQz=0.0;
                    delta = 0L;
                    for (int p = oP; p < oP+nP; p++) {
                        for (int q = oQ; q < oQ+nQ; q++) {
                            double dAdB = -DPQp[p][q];
                            PxPx += dAdB * PxPxBuf[delta];
                            PxPy += dAdB * PxPyBuf[delta];
                            PxPz += dAdB * PxPzBuf[delta];
                            PxQx += dAdB * PxQxBuf[delta];
                            PxQy += dAdB * PxQyBuf[delta];
                            PxQz += dAdB * PxQzBuf[delta];
                            PyPy += dAdB * PyPyBuf[delta];
                            PyPz += dAdB * PyPzBuf[delta];
                            PyQx += dAdB * PyQxBuf[delta];
                            PyQy += dAdB * PyQyBuf[delta];
                            PyQz += dAdB * PyQzBuf[delta];
                            PzPz += dAdB * PzPzBuf[delta];
                            PzQx += dAdB * PzQxBuf[delta];
                            PzQy += dAdB * PzQyBuf[delta];
                            PzQz += dAdB * PzQzBuf[delta];
                            QxQx += dAdB * QxQxBuf[delta];
                            QxQy += dAdB * QxQyBuf[delta];
                            QxQz += dAdB * QxQzBuf[delta];
                            QyQy += dAdB * QyQyBuf[delta];
                            QyQz += dAdB * QyQzBuf[delta];
                            QzQz += dAdB * QzQzBuf[delta];
                            ++delta;
                        }

                    }
                    KHessp[Px][Px] += PxPx;
                    KHessp[Px][Py] += PxPy;
                    KHessp[Px][Pz] += PxPz;
                    KHessp[Px][Qx] += PQscale*PxQx;
                    KHessp[Px][Qy] += PxQy;
                    KHessp[Px][Qz] += PxQz;
                    KHessp[Py][Py] += PyPy;
                    KHessp[Py][Pz] += PyPz;
                    KHessp[Py][Qx] += PyQx;
                    KHessp[Py][Qy] += PQscale*PyQy;
                    KHessp[Py][Qz] += PyQz;
                    KHessp[Pz][Pz] += PzPz;
                    KHessp[Pz][Qx] += PzQx;
                    KHessp[Pz][Qy] += PzQy;
                    KHessp[Pz][Qz] += PQscale*PzQz;
                    KHessp[Qx][Qx] += QxQx;
                    KHessp[Qx][Qy] += QxQy;
                    KHessp[Qx][Qz] += QxQz;
                    KHessp[Qy][Qy] += QyQy;
                    KHessp[Qy][Qz] += QyQz;
                    KHessp[Qz][Qz] += QzQz;
                }
            }
        }
    }

    // Add permutational symmetry components missing from the above
    for (int i = 0; i < 3 * natoms; ++i) {
        for (int j = 0; j < i; ++j) {
            JHessp[i][j] = JHessp[j][i] = (JHessp[i][j] + JHessp[j][i]);
            if (do_K_) {
                KHessp[i][j] = KHessp[j][i] = (KHessp[i][j] + KHessp[j][i]);
            }
        }
    }

    // Stitch all the intermediates together to form the actual Hessian contributions

    auto tmp1 = std::make_shared<Matrix>("Tmp1", np, np*np);
    double **ptmp1 = tmp1->pointer();

    auto tmp_a = std::make_shared<Matrix>("Tmp [P][i,j]", np, na*na);
    double **ptmp_a = tmp_a->pointer();
    auto tmp_b = std::make_shared<Matrix>("Tmp [P][i,j]", np, nb*nb);
    double **ptmp_b = tmp_b->pointer();

    for (int x = 0; x < 3 * natoms; ++x) {
        for (int y = 0; y < 3 * natoms; ++y) {
            // J terms
            JHessp[x][y] += 2.0*C_DDOT(np, ddp[x], 1, dcp[y], 1);
            JHessp[x][y] -= 4.0*C_DDOT(np, ddp[x], 1, dep[y], 1);
            C_DGEMV('n', np, np, 1.0, PQp[0], np, dep[y], 1, 0.0, ptmp1[0], 1);
            JHessp[x][y] += 2.0*C_DDOT(np, dep[x], 1, ptmp1[0], 1);

            if (do_K_) {
                // K terms
                C_DGEMM('n', 'n', np, na*na, np,  1.0, PQp[0], np, dAa_ijp[y], na*na, 0.0, ptmp_a[0], na*na);
                KHessp[x][y] += 2.0*C_DDOT(static_cast<size_t> (np)*na*na, dAa_ijp[x], 1, ptmp_a[0], 1);
                C_DGEMM('n', 'n', np, na*na, np,  1.0, PQp[0], np, dea_ijp[y], na*na, 0.0, ptmp_a[0], na*na);
                KHessp[x][y] -= 4.0*C_DDOT(static_cast<size_t> (np)*na*na, dAa_ijp[x], 1, ptmp_a[0], 1);
                KHessp[x][y] += 2.0*C_DDOT(static_cast<size_t> (np)*na*na, dea_ijp[x], 1, ptmp_a[0], 1);
                if (!same_ab){
                    C_DGEMM('n', 'n', np, nb*nb, np,  1.0, PQp[0], np, dAb_ijp[y], nb*nb, 0.0, ptmp_b[0], nb*nb);
                    KHessp[x][y] += 2.0*C_DDOT(static_cast<size_t> (np)*nb*nb, dAb_ijp[x], 1, ptmp_b[0], 1);
                    C_DGEMM('n', 'n', np, nb*nb, np,  1.0, PQp[0], np, deb_ijp[y], nb*nb, 0.0, ptmp_b[0], nb*nb);
                    KHessp[x][y] -= 4.0*C_DDOT(static_cast<size_t> (np)*nb*nb, dAb_ijp[x], 1, ptmp_b[0], 1);
                    KHessp[x][y] += 2.0*C_DDOT(static_cast<size_t> (np)*nb*nb, deb_ijp[x], 1, ptmp_b[0], 1);
                }
            }
        }
    }

    // Make sure the newly added components are symmetric
    for (int i = 0; i < 3 * natoms; ++i) {
        for (int j = 0; j < i; ++j) {
            JHessp[i][j] = JHessp[j][i] = 0.5 * (JHessp[i][j] + JHessp[j][i]);
            if (do_K_) {
                KHessp[i][j] = KHessp[j][i] = 0.5 * (KHessp[i][j] + KHessp[j][i]);
            }
        }
    }

    hessians_["Coulomb"]->scale(0.5);
    if (do_K_ && !same_ab) {
        hessians_["Exchange"]->scale(0.5);
    }
}

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

}  // namespace scfgrad
}  // namespace psi
