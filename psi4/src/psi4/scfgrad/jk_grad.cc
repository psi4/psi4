/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */


#include "psi4/libmints/sieve.h"
#include "psi4/libqt/qt.h"
#include "psi4/lib3index/3index.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/psifiles.h"
#include "jk_grad.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/vector.h"
#ifdef _OPENMP
#include <omp.h>
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
std::shared_ptr<JKGrad> JKGrad::build_JKGrad(int deriv, std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary)
{
    Options& options = Process::environment.options;

    if (options.get_str("SCF_TYPE") == "DF") {

        DFJKGrad* jk = new DFJKGrad(deriv,primary,auxiliary);

        if (options["INTS_TOLERANCE"].has_changed())
            jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["PRINT"].has_changed())
            jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed())
            jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed())
            jk->set_bench(options.get_int("BENCH"));
        if (options["DF_FITTING_CONDITION"].has_changed())
            jk->set_condition(options.get_double("DF_FITTING_CONDITION"));
        if (options["DF_INTS_NUM_THREADS"].has_changed())
            jk->set_df_ints_num_threads(options.get_int("DF_INTS_NUM_THREADS"));

        return std::shared_ptr<JKGrad>(jk);
    } else if (options.get_str("SCF_TYPE") == "DIRECT" || options.get_str("SCF_TYPE") == "PK" || options.get_str("SCF_TYPE") == "OUT_OF_CORE") {

        DirectJKGrad* jk = new DirectJKGrad(deriv,primary);

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
void JKGrad::common_init()
{
    print_ = 1;
    debug_ = 0;
    bench_ = 0;

    memory_ = 32000000L;
    omp_num_threads_ = 1;
#ifdef _OPENMP
    omp_num_threads_ = Process::environment.get_n_threads();
#endif

    cutoff_ = 0.0;

    do_J_ = true;
    do_K_ = true;
    do_wK_ = false;
    omega_ = 0.0;
}
DFJKGrad::DFJKGrad(int deriv, std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary) :
    JKGrad(deriv,primary), auxiliary_(auxiliary)
{
    common_init();
}
DFJKGrad::~DFJKGrad()
{
}
void DFJKGrad::common_init()
{
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
void DFJKGrad::print_header() const
{
    if (print_) {
        outfile->Printf( "  ==> DFJKGrad: Density-Fitted SCF Gradients <==\n\n");

        outfile->Printf( "    Gradient:          %11d\n", deriv_);
        outfile->Printf( "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf( "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf( "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_)
            outfile->Printf( "    Omega:             %11.3E\n", omega_);
        outfile->Printf( "    OpenMP threads:    %11d\n", omp_num_threads_);
        outfile->Printf( "    Integrals threads: %11d\n", df_ints_num_threads_);
        outfile->Printf( "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        outfile->Printf( "    Schwarz Cutoff:    %11.0E\n", cutoff_);
        outfile->Printf( "    Fitting Condition: %11.0E\n\n", condition_);

        outfile->Printf( "   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level("outfile", print_);
    }
}
void DFJKGrad::compute_gradient()
{
    if (!do_J_ && !do_K_ && !do_wK_)
        return;

    if (!(Ca_ && Cb_ && Da_ && Db_ && Dt_))
        throw PSIEXCEPTION("Occupation/Density not set");

    // => Set up gradients <= //
    int natom = primary_->molecule()->natom();
    gradients_.clear();
    if (do_J_) {
        gradients_["Coulomb"] = SharedMatrix(new Matrix("Coulomb Gradient",natom,3));
    }
    if (do_K_) {
        gradients_["Exchange"] = SharedMatrix(new Matrix("Exchange Gradient",natom,3));
    }
    if (do_wK_) {
        throw PSIEXCEPTION("Exchange,LR gradients are not currently available with DF.");
        gradients_["Exchange,LR"] = SharedMatrix(new Matrix("Exchange,LR Gradient",natom,3));
    }

    // => Build ERI Sieve <= //
    sieve_ = std::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff_));

    // => Open temp files <= //
    psio_->open(unit_a_, PSIO_OPEN_NEW);
    psio_->open(unit_b_, PSIO_OPEN_NEW);
    psio_->open(unit_c_, PSIO_OPEN_NEW);

    // => Gradient Construction: Get in there and kill 'em all! <= //

    timer_on("JKGrad: Amn");
    build_Amn_terms();
    timer_off("JKGrad: Amn");

    timer_on("JKGrad: Awmn");
    build_Amn_lr_terms();
    timer_off("JKGrad: Awmn");

    timer_on("JKGrad: AB");
    build_AB_inv_terms();
    timer_off("JKGrad: AB");

    timer_on("JKGrad: UV");
    build_UV_terms();
    timer_off("JKGrad: UV");

    timer_on("JKGrad: ABx");
    build_AB_x_terms();
    timer_off("JKGrad: ABx");

    timer_on("JKGrad: Amnx");
    build_Amn_x_terms();
    timer_off("JKGrad: Amnx");

    timer_on("JKGrad: Awmnx");
    build_Amn_x_lr_terms();
    timer_off("JKGrad: Awmnx");

    // => Close temp files <= //
    psio_->close(unit_a_, 0);
    psio_->close(unit_b_, 0);
    psio_->close(unit_c_, 0);
}
void DFJKGrad::build_Amn_terms()
{
    // => Sizing <= //

    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //

    int max_rows;
    int maxP = auxiliary_->max_function_per_shell();
    ULI row_cost = 0L;
    row_cost += nso * (ULI) nso;
    if (do_K_ || do_wK_) {
        row_cost += nso * (ULI) na;
        row_cost += na * (ULI) na;
    }
    ULI rows = memory_ / row_cost;
    rows = (rows > naux ? naux : rows);
    rows = (rows < maxP ? maxP : rows);
    max_rows = (int) rows;

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
        c = SharedVector(new Vector("c", naux));
        cp = c->pointer();
    }

    SharedMatrix Amn;
    SharedMatrix Ami;
    SharedMatrix Aij;

    double** Amnp;
    double** Amip;
    double** Aijp;

    if (true) {
        Amn = SharedMatrix(new Matrix("Amn", max_rows, nso * (ULI) nso));
        Amnp = Amn->pointer();
    }
    if (do_K_ || do_wK_) {
        Ami = SharedMatrix(new Matrix("Ami", max_rows, nso * (ULI) na));
        Aij = SharedMatrix(new Matrix("Aij", max_rows, na * (ULI) na));
        Amip = Ami->pointer();
        Aijp = Aij->pointer();
    }

    double** Dtp = Dt_->pointer();
    double** Cap = Ca_->pointer();
    double** Cbp = Cb_->pointer();

    psio_address next_Aija = PSIO_ZERO;
    psio_address next_Aijb = PSIO_ZERO;

    // => Integrals <= //

    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->eri()));
    }

    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {

        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop  = Pstarts[block+1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop  = (Pstop == auxiliary_->nshell() ? naux : auxiliary_->shell(Pstop ).function_index());
        int np = pstop - pstart;

        // > Clear Integrals Register < //
        ::memset((void*) Amnp[0], '\0', sizeof(double) * np * nso * nso);

        // > Integrals < //
        int nthread_df = df_ints_num_threads_;
#pragma omp parallel for schedule(dynamic) num_threads(nthread_df)
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {

            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int P =  PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            eri[thread]->compute_shell(P,0,M,N);

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
                        Amnp[p + oP][(m + oM) * nso + (n + oN)] =
                                Amnp[p + oP][(n + oN) * nso + (m + oM)] =
                                *buffer++;
                    }
                }
            }

        }

        // > (A|mn) D_mn -> c_A < //
        if (do_J_) {
            C_DGEMV('N',np,nso*(ULI)nso,1.0,Amnp[0],nso*(ULI)nso,Dtp[0],1,0.0,&cp[pstart],1);
        }

        // > Alpha < //
        if (do_K_ || do_wK_) {
            // > (A|mn) C_ni -> (A|mi) < //
            C_DGEMM('N','N',np*(ULI)nso,na,nso,1.0,Amnp[0],nso,Cap[0],na,0.0,Amip[0],na);

            // > (A|mi) C_mj -> (A|ij) < //
#pragma omp parallel for
            for (int p = 0; p < np; p++) {
                C_DGEMM('T','N',na,na,nso,1.0,Amip[p],na,Cap[0],na,0.0,&Aijp[0][p * (ULI) na * na],na);
            }

            // > Stripe < //
            psio_->write(unit_a_, "(A|ij)", (char*) Aijp[0], sizeof(double) * np * na * na, next_Aija, &next_Aija);
        }

        // > Beta < //
        if (!restricted && (do_K_ || do_wK_)) {
            // > (A|mn) C_ni -> (A|mi) < //
            C_DGEMM('N','N',np*(ULI)nso,nb,nso,1.0,Amnp[0],nso,Cbp[0],nb,0.0,Amip[0],na);

            // > (A|mi) C_mj -> (A|ij) < //
#pragma omp parallel for
            for (int p = 0; p < np; p++) {
                C_DGEMM('T','N',nb,nb,nso,1.0,Amip[p],na,Cbp[0],nb,0.0,&Aijp[0][p * (ULI) nb * nb],nb);
            }

            // > Stripe < //
            psio_->write(unit_b_, "(A|ij)", (char*) Aijp[0], sizeof(double) * np * nb * nb, next_Aijb, &next_Aijb);
        }
    }

    if (do_J_) {
        psio_->write_entry(unit_c_, "c", (char*) cp, sizeof(double) * naux);
    }
}
void DFJKGrad::build_Amn_lr_terms()
{
    if (!do_wK_) return;

    // => Sizing <= //

    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //

    int max_rows;
    int maxP = auxiliary_->max_function_per_shell();
    ULI row_cost = 0L;
    row_cost += nso * (ULI) nso;
    row_cost += nso * (ULI) na;
    row_cost += na * (ULI) na;
    ULI rows = memory_ / row_cost;
    rows = (rows > naux ? naux : rows);
    rows = (rows < maxP ? maxP : rows);
    max_rows = (int) rows;

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

    Amn = SharedMatrix(new Matrix("Amn", max_rows, nso * (ULI) nso));
    Ami = SharedMatrix(new Matrix("Ami", max_rows, nso * (ULI) na));
    Aij = SharedMatrix(new Matrix("Aij", max_rows, na * (ULI) na));

    Amnp = Amn->pointer();
    Amip = Ami->pointer();
    Aijp = Aij->pointer();

    double** Cap = Ca_->pointer();
    double** Cbp = Cb_->pointer();

    psio_address next_Aija = PSIO_ZERO;
    psio_address next_Aijb = PSIO_ZERO;

    // => Integrals <= //

    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_)));
    }

    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {

        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop  = Pstarts[block+1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop  = (Pstop == auxiliary_->nshell() ? naux : auxiliary_->shell(Pstop ).function_index());
        int np = pstop - pstart;

        // > Clear Integrals Register < //
        ::memset((void*) Amnp[0], '\0', sizeof(double) * np * nso * nso);

        // > Integrals < //
        int nthread_df = df_ints_num_threads_;
#pragma omp parallel for schedule(dynamic) num_threads(nthread_df)
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {

            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int P =  PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            eri[thread]->compute_shell(P,0,M,N);

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
                        Amnp[p + oP][(m + oM) * nso + (n + oN)] =
                                Amnp[p + oP][(n + oN) * nso + (m + oM)] =
                                *buffer++;
                    }
                }
            }

        }

        // > Alpha < //
        if (true) {
            // > (A|mn) C_ni -> (A|mi) < //
            C_DGEMM('N','N',np*(ULI)nso,na,nso,1.0,Amnp[0],nso,Cap[0],na,0.0,Amip[0],na);

            // > (A|mi) C_mj -> (A|ij) < //
#pragma omp parallel for
            for (int p = 0; p < np; p++) {
                C_DGEMM('T','N',na,na,nso,1.0,Amip[p],na,Cap[0],na,0.0,&Aijp[0][p * (ULI) na * na],na);
            }

            // > Stripe < //
            psio_->write(unit_a_, "(A|w|ij)", (char*) Aijp[0], sizeof(double) * np * na * na, next_Aija, &next_Aija);
        }

        // > Beta < //
        if (!restricted) {
            // > (A|mn) C_ni -> (A|mi) < //
            C_DGEMM('N','N',np*(ULI)nso,nb,nso,1.0,Amnp[0],nso,Cbp[0],nb,0.0,Amip[0],na);

            // > (A|mi) C_mj -> (A|ij) < //
#pragma omp parallel for
            for (int p = 0; p < np; p++) {
                C_DGEMM('T','N',nb,nb,nso,1.0,Amip[p],na,Cbp[0],nb,0.0,&Aijp[0][p * (ULI) nb * nb],nb);
            }

            // > Stripe < //
            psio_->write(unit_b_, "(A|w|ij)", (char*) Aijp[0], sizeof(double) * np * nb * nb, next_Aijb, &next_Aijb);
        }
    }
}
void DFJKGrad::build_AB_inv_terms()
{

    // => Sizing <= //

    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    // => Fitting Metric Full Inverse <= //

    std::shared_ptr<FittingMetric> metric(new FittingMetric(auxiliary_, true));
    metric->form_full_eig_inverse();
    SharedMatrix J = metric->get_metric();
    double** Jp = J->pointer();

    // => d_A = (A|B)^{-1} c_B <= //
    if (do_J_) {
        SharedVector c(new Vector("c", naux));
        SharedVector d(new Vector("d", naux));
        double* cp = c->pointer();
        double* dp = d->pointer();

        psio_->read_entry(unit_c_, "c", (char*) cp, sizeof(double) * naux);

        C_DGEMV('N',naux,naux,1.0,Jp[0],naux,cp,1,0.0,dp,1);

        psio_->write_entry(unit_c_, "c", (char*) dp, sizeof(double) * naux);
    }

    if (!(do_K_ || do_wK_))
        return;

    int max_cols;
    ULI effective_memory = memory_ - 1L * naux * naux;
    ULI col_cost = 2L * naux;
    ULI cols = effective_memory / col_cost;
    cols = (cols > na * (ULI) na ? na * (ULI) na : cols);
    cols = (cols < na ? na : cols);
    max_cols = (int) cols;

    SharedMatrix Aij(new Matrix("Aij", naux, max_cols));
    SharedMatrix Bij(new Matrix("Bij", naux, max_cols));
    double** Aijp = Aij->pointer();
    double** Bijp = Bij->pointer();

    // > Alpha < //
    if (true) {
        psio_address next_Aija = PSIO_ZERO;

        for (long int ij = 0L; ij < na *(ULI) na; ij += max_cols) {
            int ncols = (ij + max_cols >= na * (ULI) na ? na * (ULI) na - ij : max_cols);

            // > Read < //
            for (int Q = 0; Q < naux; Q++) {
                next_Aija = psio_get_address(PSIO_ZERO,sizeof(double) * (Q * (ULI) na * na + ij));
                psio_->read(unit_a_,"(A|ij)",(char*) Aijp[Q], sizeof(double) * ncols, next_Aija, &next_Aija);
            }

            // > GEMM <//
            C_DGEMM('N','N',naux,ncols,naux,1.0,Jp[0],naux,Aijp[0],max_cols,0.0,Bijp[0],max_cols);

            // > Stripe < //
            for (int Q = 0; Q < naux; Q++) {
                next_Aija = psio_get_address(PSIO_ZERO,sizeof(double) * (Q * (ULI) na * na + ij));
                psio_->write(unit_a_,"(A|ij)",(char*) Bijp[Q], sizeof(double) * ncols, next_Aija, &next_Aija);
            }

        }
    }

    // > Beta < //
    if (!restricted) {
        psio_address next_Aijb = PSIO_ZERO;

        for (long int ij = 0L; ij < nb *(ULI) nb; ij += max_cols) {
            int ncols = (ij + max_cols >= nb * (ULI) nb ? nb * (ULI) nb - ij : max_cols);

            // > Read < //
            for (int Q = 0; Q < naux; Q++) {
                next_Aijb = psio_get_address(PSIO_ZERO,sizeof(double) * (Q * (ULI) nb * nb + ij));
                psio_->read(unit_b_,"(A|ij)",(char*) Aijp[Q], sizeof(double) * ncols, next_Aijb, &next_Aijb);
            }

            // > GEMM <//
            C_DGEMM('N','N',naux,ncols,naux,1.0,Jp[0],naux,Aijp[0],max_cols,0.0,Bijp[0],max_cols);

            // > Stripe < //
            for (int Q = 0; Q < naux; Q++) {
                next_Aijb = psio_get_address(PSIO_ZERO,sizeof(double) * (Q * (ULI) nb * nb + ij));
                psio_->write(unit_b_,"(A|ij)",(char*) Bijp[Q], sizeof(double) * ncols, next_Aijb, &next_Aijb);
            }

        }
    }

    if (!do_wK_)
        return;

    // > Alpha < //
    if (true) {
        psio_address next_Aija = PSIO_ZERO;

        for (long int ij = 0L; ij < na *(ULI) na; ij += max_cols) {
            int ncols = (ij + max_cols >= na * (ULI) na ? na * (ULI) na - ij : max_cols);

            // > Read < //
            for (int Q = 0; Q < naux; Q++) {
                next_Aija = psio_get_address(PSIO_ZERO,sizeof(double) * (Q * (ULI) na * na + ij));
                psio_->read(unit_a_,"(A|w|ij)",(char*) Aijp[Q], sizeof(double) * ncols, next_Aija, &next_Aija);
            }

            // > GEMM <//
            C_DGEMM('N','N',naux,ncols,naux,1.0,Jp[0],naux,Aijp[0],max_cols,0.0,Bijp[0],max_cols);

            // > Stripe < //
            for (int Q = 0; Q < naux; Q++) {
                next_Aija = psio_get_address(PSIO_ZERO,sizeof(double) * (Q * (ULI) na * na + ij));
                psio_->write(unit_a_,"(A|w|ij)",(char*) Bijp[Q], sizeof(double) * ncols, next_Aija, &next_Aija);
            }

        }
    }

    // > Beta < //
    if (!restricted) {
        psio_address next_Aijb = PSIO_ZERO;

        for (long int ij = 0L; ij < nb *(ULI) nb; ij += max_cols) {
            int ncols = (ij + max_cols >= nb * (ULI) nb ? nb * (ULI) nb - ij : max_cols);

            // > Read < //
            for (int Q = 0; Q < naux; Q++) {
                next_Aijb = psio_get_address(PSIO_ZERO,sizeof(double) * (Q * (ULI) nb * nb + ij));
                psio_->read(unit_b_,"(A|w|ij)",(char*) Aijp[Q], sizeof(double) * ncols, next_Aijb, &next_Aijb);
            }

            // > GEMM <//
            C_DGEMM('N','N',naux,ncols,naux,1.0,Jp[0],naux,Aijp[0],max_cols,0.0,Bijp[0],max_cols);

            // > Stripe < //
            for (int Q = 0; Q < naux; Q++) {
                next_Aijb = psio_get_address(PSIO_ZERO,sizeof(double) * (Q * (ULI) nb * nb + ij));
                psio_->write(unit_b_,"(A|w|ij)",(char*) Bijp[Q], sizeof(double) * ncols, next_Aijb, &next_Aijb);
            }

        }
    }
}
void DFJKGrad::build_UV_terms()
{
    if (!(do_K_ || do_wK_))
        return;

    // => Sizing <= //

    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    SharedMatrix V = SharedMatrix(new Matrix("W", naux, naux));
    double** Vp = V->pointer();

    // => Memory Constraints <= //

    int max_rows;
    ULI effective_memory = memory_ - 1L * naux * naux;
    ULI row_cost = 2L * na * (ULI) na;
    ULI rows = memory_ / row_cost;
    rows = (rows > naux ? naux : rows);
    rows = (rows < 1L ? 1L : rows);
    max_rows = (int) rows;

    // => Temporary Buffers <= //

    SharedMatrix Aij(new Matrix("Aij", max_rows, na*(ULI)na));
    SharedMatrix Bij(new Matrix("Bij", max_rows, na*(ULI)na));
    double** Aijp = Aij->pointer();
    double** Bijp = Bij->pointer();

    // => V < = //

    // > Alpha < //
    if (true) {
        psio_address next_Aij = PSIO_ZERO;
        for (int P = 0; P < naux; P += max_rows) {
            psio_address next_Bij = PSIO_ZERO;
            int nP = (P + max_rows >= naux ? naux - P : max_rows);
            psio_->read(unit_a_,"(A|ij)",(char*) Aijp[0], sizeof(double)*nP*na*na, next_Aij, &next_Aij);
            for (int Q = 0; Q < naux; Q += max_rows) {
                int nQ = (Q + max_rows >= naux ? naux - Q : max_rows);
                psio_->read(unit_a_,"(A|ij)",(char*) Bijp[0], sizeof(double)*nQ*na*na, next_Bij, &next_Bij);

                C_DGEMM('N','T',nP,nQ,na*(ULI)na,1.0,Aijp[0],na*(ULI)na,Bijp[0],na*(ULI)na,0.0,&Vp[P][Q],naux);
            }
        }
    }
    // > Beta < //
    if (!restricted) {
        psio_address next_Aij = PSIO_ZERO;
        for (int P = 0; P < naux; P += max_rows) {
            psio_address next_Bij = PSIO_ZERO;
            int nP = (P + max_rows >= naux ? naux - P : max_rows);
            psio_->read(unit_b_,"(A|ij)",(char*) Aijp[0], sizeof(double)*nP*nb*nb, next_Aij, &next_Aij);
            for (int Q = 0; Q < naux; Q += max_rows) {
                int nQ = (Q + max_rows >= naux ? naux - Q : max_rows);
                psio_->read(unit_b_,"(A|ij)",(char*) Bijp[0], sizeof(double)*nQ*nb*nb, next_Bij, &next_Bij);

                C_DGEMM('N','T',nP,nQ,nb*(ULI)nb,1.0,Aijp[0],nb*(ULI)nb,Bijp[0],nb*(ULI)nb,1.0,&Vp[P][Q],naux);
            }
        }
    } else {
        V->scale(2.0);
    }
    psio_->write_entry(unit_c_,"V",(char*) Vp[0], sizeof(double) * naux * naux);

    if (!do_wK_)
        return;

    // => W < = //

    // > Alpha < //
    if (true) {
        psio_address next_Aij = PSIO_ZERO;
        for (int P = 0; P < naux; P += max_rows) {
            psio_address next_Bij = PSIO_ZERO;
            int nP = (P + max_rows >= naux ? naux - P : max_rows);
            psio_->read(unit_a_,"(A|ij)",(char*) Aijp[0], sizeof(double)*nP*na*na, next_Aij, &next_Aij);
            for (int Q = 0; Q < naux; Q += max_rows) {
                int nQ = (Q + max_rows >= naux ? naux - Q : max_rows);
                psio_->read(unit_a_,"(A|w|ij)",(char*) Bijp[0], sizeof(double)*nQ*na*na, next_Bij, &next_Bij);

                C_DGEMM('N','T',nP,nQ,na*(ULI)na,1.0,Aijp[0],na*(ULI)na,Bijp[0],na*(ULI)na,0.0,&Vp[P][Q],naux);
            }
        }
    }
    // > Beta < //
    if (!restricted) {
        psio_address next_Aij = PSIO_ZERO;
        for (int P = 0; P < naux; P += max_rows) {
            psio_address next_Bij = PSIO_ZERO;
            int nP = (P + max_rows >= naux ? naux - P : max_rows);
            psio_->read(unit_b_,"(A|ij)",(char*) Aijp[0], sizeof(double)*nP*nb*nb, next_Aij, &next_Aij);
            for (int Q = 0; Q < naux; Q += max_rows) {
                int nQ = (Q + max_rows >= naux ? naux - Q : max_rows);
                psio_->read(unit_b_,"(A|w|ij)",(char*) Bijp[0], sizeof(double)*nQ*nb*nb, next_Bij, &next_Bij);

                C_DGEMM('N','T',nP,nQ,nb*(ULI)nb,1.0,Aijp[0],nb*(ULI)nb,Bijp[0],nb*(ULI)nb,1.0,&Vp[P][Q],naux);
            }
        }
    } else {
        V->scale(2.0);
    }
    psio_->write_entry(unit_c_,"W",(char*) Vp[0], sizeof(double) * naux * naux);
}
void DFJKGrad::build_AB_x_terms()
{

    // => Sizing <= //

    int natom = primary_->molecule()->natom();
    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();

    // => Forcing Terms/Gradients <= //
    SharedMatrix V;
    SharedMatrix W;
    SharedVector d;

    double** Vp;
    double** Wp;
    double*  dp;

    if (do_J_) {
        d = SharedVector(new Vector("d", naux));
        dp = d->pointer();
        psio_->read_entry(unit_c_, "c", (char*) dp, sizeof(double) * naux);
    }
    if (do_K_) {
        V = SharedMatrix(new Matrix("V", naux, naux));
        Vp = V->pointer();
        psio_->read_entry(unit_c_, "V", (char*) Vp[0], sizeof(double) * naux * naux);
    }
    if (do_wK_) {
        W = SharedMatrix(new Matrix("W", naux, naux));
        Wp = W->pointer();
        psio_->read_entry(unit_c_, "W", (char*) Wp[0], sizeof(double) * naux * naux);
    }

    // => Integrals <= //

    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_,BasisSet::zero_ao_basis_set(),auxiliary_,BasisSet::zero_ao_basis_set()));
    std::vector<std::shared_ptr<TwoBodyAOInt> > Jint;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        Jint.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)));
    }

    // => Temporary Gradients <= //

    std::vector<SharedMatrix> Jtemps;
    std::vector<SharedMatrix> Ktemps;
    std::vector<SharedMatrix> wKtemps;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        if (do_J_) {
            Jtemps.push_back(SharedMatrix(new Matrix("Jtemp", natom, 3)));
        }
        if (do_K_) {
            Ktemps.push_back(SharedMatrix(new Matrix("Ktemp", natom, 3)));
        }
        if (do_wK_) {
            wKtemps.push_back(SharedMatrix(new Matrix("wKtemp", natom, 3)));
        }
    }

    std::vector<std::pair<int,int> > PQ_pairs;
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        for (int Q = 0; Q <= P; Q++) {
            PQ_pairs.push_back(std::pair<int,int>(P,Q));
        }
    }

    int nthread_df = df_ints_num_threads_;
#pragma omp parallel for schedule(dynamic) num_threads(nthread_df)
    for (long int PQ = 0L; PQ < PQ_pairs.size(); PQ++) {

        int P = PQ_pairs[PQ].first;
        int Q = PQ_pairs[PQ].second;

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        Jint[thread]->compute_shell_deriv1(P,0,Q,0);
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
        const double *Px = buffer + 0*ncart;
        const double *Py = buffer + 1*ncart;
        const double *Pz = buffer + 2*ncart;
        const double *Qx = buffer + 3*ncart;
        const double *Qy = buffer + 4*ncart;
        const double *Qz = buffer + 5*ncart;

        double perm = (P == Q ? 1.0 : 2.0);

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
            for (int q = 0; q < nQ; q++) {

                if (do_J_) {
                    double Uval = 0.5 * perm * dp[p + oP] * dp[q + oQ];
                    grad_Jp[aP][0] -= Uval * (*Px);
                    grad_Jp[aP][1] -= Uval * (*Py);
                    grad_Jp[aP][2] -= Uval * (*Pz);
                    grad_Jp[aQ][0] -= Uval * (*Qx);
                    grad_Jp[aQ][1] -= Uval * (*Qy);
                    grad_Jp[aQ][2] -= Uval * (*Qz);
                }

                if (do_K_) {
                    double Vval = 0.5 * perm * Vp[p + oP][q + oQ];
                    grad_Kp[aP][0] -= Vval * (*Px);
                    grad_Kp[aP][1] -= Vval * (*Py);
                    grad_Kp[aP][2] -= Vval * (*Pz);
                    grad_Kp[aQ][0] -= Vval * (*Qx);
                    grad_Kp[aQ][1] -= Vval * (*Qy);
                    grad_Kp[aQ][2] -= Vval * (*Qz);
                }

                if (do_wK_) {
                    double Wval = 0.5 * perm * (0.5 * (Wp[p + oP][q + oQ] + Wp[q + oQ][p + oP]));
                    grad_wKp[aP][0] -= Wval * (*Px);
                    grad_wKp[aP][1] -= Wval * (*Py);
                    grad_wKp[aP][2] -= Wval * (*Pz);
                    grad_wKp[aQ][0] -= Wval * (*Qx);
                    grad_wKp[aQ][1] -= Wval * (*Qy);
                    grad_wKp[aQ][2] -= Wval * (*Qz);
                }

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

    //gradients_["Coulomb"]->zero();
    //gradients_["Exchange"]->zero();
    //gradients_["Exchange,LR"]->zero();

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

    //gradients_["Coulomb"]->print();
    //gradients_["Exchange"]->print();
    //gradients_["Exchange,LR"]->print();
}
void DFJKGrad::build_Amn_x_terms()
{

    // => Sizing <= //

    int natom = primary_->molecule()->natom();
    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //

    int max_rows;
    if (do_K_ || do_wK_) {
        int maxP = auxiliary_->max_function_per_shell();
        ULI row_cost = 0L;
        row_cost += nso * (ULI) nso;
        if (do_wK_) {
            row_cost += nso * (ULI) nso;
        }
        row_cost += nso * (ULI) na;
        row_cost += na * (ULI) na;
        ULI rows = memory_ / row_cost;
        rows = (rows > naux ? naux : rows);
        rows = (rows < maxP ? maxP : rows);
        max_rows = (int) rows;
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
        d = SharedVector(new Vector("d", naux));
        dp = d->pointer();
        psio_->read_entry(unit_c_, "c", (char*) dp, sizeof(double) * naux);
    }

    SharedMatrix Jmn;
    SharedMatrix Kmn;
    SharedMatrix Ami;
    SharedMatrix Aij;

    double** Jmnp;
    double** Kmnp;
    double** Amip;
    double** Aijp;

    if (do_K_ || do_wK_) {
        Jmn = SharedMatrix(new Matrix("Jmn", max_rows, nso * (ULI) nso));
        Ami = SharedMatrix(new Matrix("Ami", max_rows, nso * (ULI) na));
        Aij = SharedMatrix(new Matrix("Aij", max_rows, na * (ULI) na));
        Jmnp = Jmn->pointer();
        Amip = Ami->pointer();
        Aijp = Aij->pointer();
    }
    if (do_wK_) {
        Kmn = SharedMatrix(new Matrix("Kmn", max_rows, nso * (ULI) nso));
        Kmnp = Kmn->pointer();
    }

    double** Dtp = Dt_->pointer();
    double** Cap = Ca_->pointer();
    double** Cbp = Cb_->pointer();

    psio_address next_Aija = PSIO_ZERO;
    psio_address next_Aijb = PSIO_ZERO;
    psio_address next_Awija = PSIO_ZERO;
    psio_address next_Awijb = PSIO_ZERO;

    // => Integrals <= //

    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)));
    }

    // => Temporary Gradients <= //

    std::vector<SharedMatrix> Jtemps;
    std::vector<SharedMatrix> Ktemps;
    std::vector<SharedMatrix> wKtemps;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        if (do_J_) {
            Jtemps.push_back(SharedMatrix(new Matrix("Jtemp", natom, 3)));
        }
        if (do_K_) {
            Ktemps.push_back(SharedMatrix(new Matrix("Ktemp", natom, 3)));
        }
        if (do_wK_) {
            wKtemps.push_back(SharedMatrix(new Matrix("wKtemp", natom, 3)));
        }
    }

    // => R/U doubling factor <= //

    double factor = (restricted ? 2.0 : 1.0);

    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {

        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop  = Pstarts[block+1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop  = (Pstop == auxiliary_->nshell() ? naux : auxiliary_->shell(Pstop ).function_index());
        int np = pstop - pstart;

        // => J_mn^A <= //

        // > Alpha < //
        if (do_K_ || do_wK_) {

            // > Stripe < //
            psio_->read(unit_a_, "(A|ij)", (char*) Aijp[0], sizeof(double) * np * na * na, next_Aija, &next_Aija);

            // > (A|ij) C_mi -> (A|mj) < //
#pragma omp parallel for
            for (int P = 0; P < np; P++) {
                C_DGEMM('N','N',nso,na,na,1.0,Cap[0],na,&Aijp[0][P * (ULI) na * na],na,0.0,Amip[P],na);
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N','T',np * (ULI) nso, nso, na, factor, Amip[0], na, Cap[0], na, 0.0, Jmnp[0], nso);
        }

        // > Beta < //
        if (!restricted && (do_K_ || do_wK_)) {

            // > Stripe < //
            psio_->read(unit_b_, "(A|ij)", (char*) Aijp[0], sizeof(double) * np * nb * nb, next_Aijb, &next_Aijb);

            // > (A|ij) C_mi -> (A|mj) < //
#pragma omp parallel for
            for (int P = 0; P < np; P++) {
                C_DGEMM('N','N',nso,nb,nb,1.0,Cbp[0],nb,&Aijp[0][P* (ULI) nb * nb],nb,0.0,Amip[P],na);
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N','T',np * (ULI) nso, nso, nb, 1.0, Amip[0], na, Cbp[0], nb, 1.0, Jmnp[0], nso);
        }

        // => K_mn^A <= //

        // > Alpha < //
        if (do_wK_) {

            // > Stripe < //
            psio_->read(unit_a_, "(A|w|ij)", (char*) Aijp[0], sizeof(double) * np * na * na, next_Awija, &next_Awija);

            // > (A|ij) C_mi -> (A|mj) < //
#pragma omp parallel for
            for (int P = 0; P < np; P++) {
                C_DGEMM('N','N',nso,na,na,1.0,Cap[0],na,&Aijp[0][P * (ULI) na * na],na,0.0,Amip[P],na);
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N','T',np * (ULI) nso, nso, na, factor, Amip[0], na, Cap[0], na, 0.0, Kmnp[0], nso);
        }

        // > Beta < //
        if (!restricted && do_wK_) {

            // > Stripe < //
            psio_->read(unit_b_, "(A|w|ij)", (char*) Aijp[0], sizeof(double) * np * nb * nb, next_Awijb, &next_Awijb);

            // > (A|ij) C_mi -> (A|mj) < //
#pragma omp parallel for
            for (int P = 0; P < np; P++) {
                C_DGEMM('N','N',nso,nb,nb,1.0,Cbp[0],nb,&Aijp[0][P* (ULI) nb * nb],nb,0.0,Amip[P],na);
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N','T',np * (ULI) nso, nso, nb, 1.0, Amip[0], na, Cbp[0], nb, 1.0, Kmnp[0], nso);
        }

        // > Integrals < //
        int nthread_df = df_ints_num_threads_;
#pragma omp parallel for schedule(dynamic) num_threads(nthread_df)
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {

            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int P =  PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            eri[thread]->compute_shell_deriv1(P,0,M,N);

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
            const double *Px = buffer + 0*ncart;
            const double *Py = buffer + 1*ncart;
            const double *Pz = buffer + 2*ncart;
            const double *Mx = buffer + 3*ncart;
            const double *My = buffer + 4*ncart;
            const double *Mz = buffer + 5*ncart;
            const double *Nx = buffer + 6*ncart;
            const double *Ny = buffer + 7*ncart;
            const double *Nz = buffer + 8*ncart;

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

                        if (do_K_) {
                            double Jval = 1.0 * perm * Jmnp[p + oP][(m + oM) * nso + (n + oN)];
                            grad_Kp[aP][0] += Jval * (*Px);
                            grad_Kp[aP][1] += Jval * (*Py);
                            grad_Kp[aP][2] += Jval * (*Pz);
                            grad_Kp[aM][0] += Jval * (*Mx);
                            grad_Kp[aM][1] += Jval * (*My);
                            grad_Kp[aM][2] += Jval * (*Mz);
                            grad_Kp[aN][0] += Jval * (*Nx);
                            grad_Kp[aN][1] += Jval * (*Ny);
                            grad_Kp[aN][2] += Jval * (*Nz);
                        }

                        if (do_wK_) {
                            double Kval = 0.5 * perm * Kmnp[p + oP][(m + oM) * nso + (n + oN)];
                            grad_wKp[aP][0] += Kval * (*Px);
                            grad_wKp[aP][1] += Kval * (*Py);
                            grad_wKp[aP][2] += Kval * (*Pz);
                            grad_wKp[aM][0] += Kval * (*Mx);
                            grad_wKp[aM][1] += Kval * (*My);
                            grad_wKp[aM][2] += Kval * (*Mz);
                            grad_wKp[aN][0] += Kval * (*Nx);
                            grad_wKp[aN][1] += Kval * (*Ny);
                            grad_wKp[aN][2] += Kval * (*Nz);
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
        }
    }

    // => Temporary Gradient Reduction <= //

    //gradients_["Coulomb"]->zero();
    //gradients_["Exchange"]->zero();
    //gradients_["Exchange,LR"]->zero();

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

    //gradients_["Coulomb"]->print();
    //gradients_["Exchange"]->print();
    //gradients_["Exchange,LR"]->print();
}
void DFJKGrad::build_Amn_x_lr_terms()
{
    if (!do_wK_) return;

    // => Sizing <= //

    int natom = primary_->molecule()->natom();
    int nso = primary_->nbf();
    int naux = auxiliary_->nbf();
    int na = Ca_->colspi()[0];
    int nb = Cb_->colspi()[0];

    bool restricted = (Ca_ == Cb_);

    const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();
    int npairs = shell_pairs.size();

    // => Memory Constraints <= //

    int max_rows;
    int maxP = auxiliary_->max_function_per_shell();
    ULI row_cost = 0L;
    row_cost += nso * (ULI) nso;
    row_cost += nso * (ULI) na;
    row_cost += na * (ULI) na;
    ULI rows = memory_ / row_cost;
    rows = (rows > naux ? naux : rows);
    rows = (rows < maxP ? maxP : rows);
    max_rows = (int) rows;

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

    SharedMatrix Jmn;
    SharedMatrix Ami;
    SharedMatrix Aij;

    double** Jmnp;
    double** Amip;
    double** Aijp;

    Jmn = SharedMatrix(new Matrix("Jmn", max_rows, nso * (ULI) nso));
    Ami = SharedMatrix(new Matrix("Ami", max_rows, nso * (ULI) na));
    Aij = SharedMatrix(new Matrix("Aij", max_rows, na * (ULI) na));
    Jmnp = Jmn->pointer();
    Amip = Ami->pointer();
    Aijp = Aij->pointer();

    double** Cap = Ca_->pointer();
    double** Cbp = Cb_->pointer();

    psio_address next_Aija = PSIO_ZERO;
    psio_address next_Aijb = PSIO_ZERO;

    // => Integrals <= //

    std::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<std::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(std::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_,1)));
    }

    // => Temporary Gradients <= //

    std::vector<SharedMatrix> wKtemps;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        wKtemps.push_back(SharedMatrix(new Matrix("wKtemp", natom, 3)));
    }

    // => R/U doubling factor <= //

    double factor = (restricted ? 2.0 : 1.0);

    // => Master Loop <= //

    for (int block = 0; block < Pstarts.size() - 1; block++) {

        // > Sizing < //

        int Pstart = Pstarts[block];
        int Pstop  = Pstarts[block+1];
        int NP = Pstop - Pstart;

        int pstart = auxiliary_->shell(Pstart).function_index();
        int pstop  = (Pstop == auxiliary_->nshell() ? naux : auxiliary_->shell(Pstop ).function_index());
        int np = pstop - pstart;

        // => J_mn^A <= //

        // > Alpha < //
        if (true) {

            // > Stripe < //
            psio_->read(unit_a_, "(A|ij)", (char*) Aijp[0], sizeof(double) * np * na * na, next_Aija, &next_Aija);

            // > (A|ij) C_mi -> (A|mj) < //
#pragma omp parallel for
            for (int P = 0; P < np; P++) {
                C_DGEMM('N','N',nso,na,na,1.0,Cap[0],na,&Aijp[0][P * (ULI) na * na],na,0.0,Amip[P],na);
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N','T',np * (ULI) nso, nso, na, factor, Amip[0], na, Cap[0], na, 0.0, Jmnp[0], nso);
        }

        // > Beta < //
        if (!restricted) {

            // > Stripe < //
            psio_->read(unit_b_, "(A|ij)", (char*) Aijp[0], sizeof(double) * np * nb * nb, next_Aijb, &next_Aijb);

            // > (A|ij) C_mi -> (A|mj) < //
#pragma omp parallel for
            for (int P = 0; P < np; P++) {
                C_DGEMM('N','N',nso,nb,nb,1.0,Cbp[0],nb,&Aijp[0][P* (ULI) nb * nb],nb,0.0,Amip[P],na);
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N','T',np * (ULI) nso, nso, nb, 1.0, Amip[0], na, Cbp[0], nb, 1.0, Jmnp[0], nso);
        }

        // > Integrals < //
        int nthread_df = df_ints_num_threads_;
#pragma omp parallel for schedule(dynamic) num_threads(nthread_df)
        for (long int PMN = 0L; PMN < NP * npairs; PMN++) {

            int thread = 0;
#ifdef _OPENMP
            thread = omp_get_thread_num();
#endif

            int P =  PMN / npairs + Pstart;
            int MN = PMN % npairs;
            int M = shell_pairs[MN].first;
            int N = shell_pairs[MN].second;

            eri[thread]->compute_shell_deriv1(P,0,M,N);

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
            const double *Px = buffer + 0*ncart;
            const double *Py = buffer + 1*ncart;
            const double *Pz = buffer + 2*ncart;
            const double *Mx = buffer + 3*ncart;
            const double *My = buffer + 4*ncart;
            const double *Mz = buffer + 5*ncart;
            const double *Nx = buffer + 6*ncart;
            const double *Ny = buffer + 7*ncart;
            const double *Nz = buffer + 8*ncart;

            double perm = (M == N ? 1.0 : 2.0);

            double** grad_wKp = wKtemps[thread]->pointer();

            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {

                        double Kval = 0.5 * perm * Jmnp[p + oP][(m + oM) * nso + (n + oN)];
                        grad_wKp[aP][0] += Kval * (*Px);
                        grad_wKp[aP][1] += Kval * (*Py);
                        grad_wKp[aP][2] += Kval * (*Pz);
                        grad_wKp[aM][0] += Kval * (*Mx);
                        grad_wKp[aM][1] += Kval * (*My);
                        grad_wKp[aM][2] += Kval * (*Mz);
                        grad_wKp[aN][0] += Kval * (*Nx);
                        grad_wKp[aN][1] += Kval * (*Ny);
                        grad_wKp[aN][2] += Kval * (*Nz);

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

    //gradients_["Exchange,LR"]->zero();

    for (int t = 0; t < df_ints_num_threads_; t++) {
        gradients_["Exchange,LR"]->add(wKtemps[t]);
    }

    //gradients_["Exchange,LR"]->print();
}
void DFJKGrad::compute_hessian()
{

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

    // => Set up hessians <= //
    int natom = primary_->molecule()->natom();
    hessians_.clear();
    if (do_J_) {
        hessians_["Coulomb"] = SharedMatrix(new Matrix("Coulomb Hessian",3*natom,3*natom));
    }
    if (do_K_) {
        hessians_["Exchange"] = SharedMatrix(new Matrix("Exchange Hessian",3*natom,3*natom));
    }
    if (do_wK_) {
        hessians_["Exchange,LR"] = SharedMatrix(new Matrix("Exchange,LR Hessian",3*natom,3*natom));
    }


    std::shared_ptr<Molecule> mol = primary_->molecule();

    int np = auxiliary_->nbf();
    int nso = primary_->nbf();
    int nauxshell = auxiliary_->nshell();
    int nshell = primary_->nshell();
    int natoms = mol->natom();
    double **JHessp = hessians_["Coulomb"]->pointer();
    double **KHessp = hessians_["Exchange"]->pointer();

    double **Dtp = Dt_->pointer();
    double **Cap = Ca_->pointer();
    double **Cbp = Cb_->pointer();

    int na = Ca_->colspi()[0];
    std::shared_ptr<FittingMetric> metric(new FittingMetric(auxiliary_, true));
    metric->form_full_eig_inverse();
    SharedMatrix PQ = metric->get_metric();
    double** PQp = PQ->pointer();

    SharedVector c(new Vector("c[A] = (mn|A) D[m][n]", np));
    double *cp = c->pointer();
    SharedMatrix dc(new Matrix("dc[x][A] = (mn|A)^x D[m][n]",  3*natoms, np));
    double **dcp = dc->pointer();
    SharedMatrix dAij(new Matrix("dAij[x][A,i,j] = (mn|A)^x C[m][i] C[n][j]",  3*natoms, np*na*na));
    double **dAijp = dAij->pointer();
    SharedVector d(new Vector("d[A] = Minv[A][B] C[B]", np));
    double *dp = d->pointer();
    SharedMatrix dd(new Matrix("dd[x][B] = dc[x][A] Minv[A][B]", 3*natoms, np));
    double **ddp = dd->pointer();
    SharedMatrix de(new Matrix("de[x][A] = (A|B)^x d[B] ", 3*natoms, np));
    double **dep = de->pointer();
    SharedMatrix deij(new Matrix("deij[x][A,i,j] = (A|B)^x Bij[B,i,j]", 3*natoms, np*na*na));
    double **deijp = deij->pointer();

    // Build some integral factories
    std::shared_ptr<IntegralFactory> Pmnfactory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::shared_ptr<IntegralFactory> PQfactory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), auxiliary_, BasisSet::zero_ao_basis_set()));
    std::shared_ptr<TwoBodyAOInt> Pmnint(Pmnfactory->eri(2));
    std::shared_ptr<TwoBodyAOInt> PQint(PQfactory->eri(2));
    SharedMatrix Amn(new Matrix("(A|mn)", np, nso*nso));
    SharedMatrix Ami(new Matrix("(A|mi)", np, nso*na));
    SharedMatrix Aij(new Matrix("(A|ij)", np, na*na));
    SharedMatrix Bij(new Matrix("Minv[B][A] (A|ij)", np, na*na));
    SharedMatrix Bim(new Matrix("Minv[B][A] (A|im)", np, nso*na));
    SharedMatrix Bmn(new Matrix("Minv[B][A] (A|mn)", np, nso*nso));
    SharedMatrix DPQ(new Matrix("B(P|ij) B(Q|ij)", np, np));
    double **Amnp = Amn->pointer();
    double **Amip = Ami->pointer();
    double **Aijp = Aij->pointer();
    double **Bijp = Bij->pointer();
    double **Bimp = Bim->pointer();
    double **Bmnp = Bmn->pointer();
    double **DPQp = DPQ->pointer();


    for (int P = 0; P < nauxshell; ++P){
        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();
        for(int M = 0; M < nshell; ++M){
            int nM = primary_->shell(M).nfunction();
            int oM = primary_->shell(M).function_index();
            for(int N = 0; N < nshell; ++N){
                int nN = primary_->shell(N).nfunction();
                int oN = primary_->shell(N).function_index();

                Pmnint->compute_shell(P,0,M,N);
                const double* buffer = Pmnint->buffer();

                for (int p = oP; p < oP+nP; p++) {
                    for (int m = oM; m < oM+nM; m++) {
                        for (int n = oN; n < oN+nN; n++) {
                            Amnp[p][m*nso+n] += (*buffer++);
                        }
                    }
                }
                // c[A] = (A|mn) D[m][n]
                C_DGEMV('N', np, nso*(ULI)nso, 1.0, Amnp[0], nso*(ULI)nso, Dtp[0], 1, 0.0, cp, 1);
                // (A|mj) = (A|mn) C[n][j]
                C_DGEMM('N','N',np*(ULI)nso,na,nso,1.0,Amnp[0],nso,Cap[0],na,0.0,Amip[0],na);
                // (A|ij) = (A|mj) C[m][i]
                #pragma omp parallel for
                for (int p = 0; p < np; p++) {
                    C_DGEMM('T','N',na,na,nso,1.0,Amip[p],na,Cap[0],na,0.0,&Aijp[0][p * (ULI) na * na],na);
                }

            }
        }
    }

    // d[A] = Minv[A][B] c[B]
    C_DGEMV('n', np, np, 1.0, PQp[0], np, cp, 1, 0.0, dp, 1);
    // B[B][i,j] = Minv[A][B] (A|ij)
    C_DGEMM('n','n', np, na*na, np, 1.0, PQp[0], np, Aijp[0], na*na, 0.0, Bijp[0], na*na);
    // B[B][i,n] = B[B][i,j] C[n][j]
    C_DGEMM('N', 'T', np*(ULI)na, nso, na, 1.0, Bijp[0], na, Cap[0], na, 0.0, Bimp[0], nso);
    // B[B][m,n] = C[m][i] B[B][i,n]
    #pragma omp parallel for
    for (int p = 0; p < np; p++) {
        C_DGEMM('n', 'n', nso, nso, na, 1.0, Cap[0], na, Bimp[p], nso, 0.0, Bmnp[p], nso);
    }
    // D[A][B] = B[A][ij] B[B][ij]
    C_DGEMM('n','t', np, np, na*na, 1.0, Bijp[0], na*na, Bijp[0], na*na, 0.0, DPQp[0], np);

    int maxp = auxiliary_->max_function_per_shell();
    int maxm = primary_->max_function_per_shell();
    SharedMatrix T(new Matrix("T", maxp, maxm*na));
    double **Tp = T->pointer();

    for (int P = 0; P < nauxshell; ++P){
        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();
        int Pcenter = auxiliary_->shell(P).ncenter();
        int Pncart = auxiliary_->shell(P).ncartesian();
        int Px = 3 * Pcenter + 0;
        int Py = 3 * Pcenter + 1;
        int Pz = 3 * Pcenter + 2;
        for(int M = 0; M < nshell; ++M){
            int nM = primary_->shell(M).nfunction();
            int oM = primary_->shell(M).function_index();
            int Mcenter = primary_->shell(M).ncenter();
            int Mncart = primary_->shell(M).ncartesian();
            int mx = 3 * Mcenter + 0;
            int my = 3 * Mcenter + 1;
            int mz = 3 * Mcenter + 2;
            for(int N = 0; N < nshell; ++N){
                int nN = primary_->shell(N).nfunction();
                int oN = primary_->shell(N).function_index();
                int Ncenter = primary_->shell(N).ncenter();
                int Nncart = primary_->shell(N).ncartesian();
                int nx = 3 * Ncenter + 0;
                int ny = 3 * Ncenter + 1;
                int nz = 3 * Ncenter + 2;

                size_t stride = Pncart * Mncart * Nncart;

                Pmnint->compute_shell_deriv1(P,0,M,N);
                const double* buffer = Pmnint->buffer();

                size_t delta = 0L;
                // Terms for J intermediates
                // dc[x][A] = D[m][n] (A|mn)^x
                for (int p = oP; p < oP+nP; p++) {
                    for (int m = oM; m < oM+nM; m++) {
                        for (int n = oN; n < oN+nN; n++) {
                            double Cpmn = Dtp[m][n];
                            dcp[Px][p] += Cpmn * buffer[0 * stride + delta];
                            dcp[Py][p] += Cpmn * buffer[1 * stride + delta];
                            dcp[Pz][p] += Cpmn * buffer[2 * stride + delta];
                            dcp[mx][p] += Cpmn * buffer[3 * stride + delta];
                            dcp[my][p] += Cpmn * buffer[4 * stride + delta];
                            dcp[mz][p] += Cpmn * buffer[5 * stride + delta];
                            dcp[nx][p] += Cpmn * buffer[6 * stride + delta];
                            dcp[ny][p] += Cpmn * buffer[7 * stride + delta];
                            dcp[nz][p] += Cpmn * buffer[8 * stride + delta];
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
                double *ptr = const_cast<double*>(buffer);
                C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, ptr+0*stride, nN, Cap[oN], na, 0.0, Tp[0], na);
                #pragma omp parallel for
                for(int p = 0; p < nP; ++p)
                    C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tp[0]+p*(nM*na), na, 1.0, &dAijp[Px][(p+oP)*na*na], na);
                C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, ptr+1*stride, nN, Cap[oN], na, 0.0, Tp[0], na);
                #pragma omp parallel for
                for(int p = 0; p < nP; ++p)
                    C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tp[0]+p*(nM*na), na, 1.0, &dAijp[Py][(p+oP)*na*na], na);
                C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, ptr+2*stride, nN, Cap[oN], na, 0.0, Tp[0], na);
                #pragma omp parallel for
                for(int p = 0; p < nP; ++p)
                    C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tp[0]+p*(nM*na), na, 1.0, &dAijp[Pz][(p+oP)*na*na], na);
                C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, ptr+3*stride, nN, Cap[oN], na, 0.0, Tp[0], na);
                #pragma omp parallel for
                for(int p = 0; p < nP; ++p)
                    C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tp[0]+p*(nM*na), na, 1.0, &dAijp[mx][(p+oP)*na*na], na);
                C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, ptr+4*stride, nN, Cap[oN], na, 0.0, Tp[0], na);
                #pragma omp parallel for
                for(int p = 0; p < nP; ++p)
                    C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tp[0]+p*(nM*na), na, 1.0, &dAijp[my][(p+oP)*na*na], na);
                C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, ptr+5*stride, nN, Cap[oN], na, 0.0, Tp[0], na);
                #pragma omp parallel for
                for(int p = 0; p < nP; ++p)
                    C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tp[0]+p*(nM*na), na, 1.0, &dAijp[mz][(p+oP)*na*na], na);
                C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, ptr+6*stride, nN, Cap[oN], na, 0.0, Tp[0], na);
                #pragma omp parallel for
                for(int p = 0; p < nP; ++p)
                    C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tp[0]+p*(nM*na), na, 1.0, &dAijp[nx][(p+oP)*na*na], na);
                C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, ptr+7*stride, nN, Cap[oN], na, 0.0, Tp[0], na);
                #pragma omp parallel for
                for(int p = 0; p < nP; ++p)
                    C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tp[0]+p*(nM*na), na, 1.0, &dAijp[ny][(p+oP)*na*na], na);
                C_DGEMM('n', 'n', nP*nM, na, nN, 1.0, ptr+8*stride, nN, Cap[oN], na, 0.0, Tp[0], na);
                #pragma omp parallel for
                for(int p = 0; p < nP; ++p)
                    C_DGEMM('t', 'n', na, na, nM, 1.0, Cap[oM], na, Tp[0]+p*(nM*na), na, 1.0, &dAijp[nz][(p+oP)*na*na], na);

            }
        }
    }

    // dd[x][A] = dc[x][B] Minv[B][A]
    C_DGEMM('N', 'N', 3*natoms, np, np, 1.0, dcp[0], np, PQp[0], np, 0.0, ddp[0], np);

    for (int P = 0; P < nauxshell; ++P){
        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();
        int Pcenter = auxiliary_->shell(P).ncenter();
        int Pncart = auxiliary_->shell(P).ncartesian();
        int Px = 3 * Pcenter + 0;
        int Py = 3 * Pcenter + 1;
        int Pz = 3 * Pcenter + 2;
        for(int Q = 0; Q < nauxshell; ++Q){
            int nQ = auxiliary_->shell(Q).nfunction();
            int oQ = auxiliary_->shell(Q).function_index();
            int Qcenter = auxiliary_->shell(Q).ncenter();
            int Qncart = auxiliary_->shell(Q).ncartesian();
            int Qx = 3 * Qcenter + 0;
            int Qy = 3 * Qcenter + 1;
            int Qz = 3 * Qcenter + 2;

            size_t stride = Pncart * Qncart;

            PQint->compute_shell_deriv1(P,0,Q,0);
            const double* buffer = PQint->buffer();

            size_t delta = 0L;
            // J term intermediates
            // de[x][A] = (A|B)^x d[B]
            for (int p = oP; p < oP+nP; p++) {
                for (int q = oQ; q < oQ+nQ; q++) {
                    double dq = dp[q];
                    dep[Px][p] += dq * buffer[0 * stride + delta];
                    dep[Py][p] += dq * buffer[1 * stride + delta];
                    dep[Pz][p] += dq * buffer[2 * stride + delta];
                    dep[Qx][p] += dq * buffer[3 * stride + delta];
                    dep[Qy][p] += dq * buffer[4 * stride + delta];
                    dep[Qz][p] += dq * buffer[5 * stride + delta];
                    ++delta;
                }
            }
            // K term intermediates
            // deij[x][A,i,j] <- (A|B)^x Bij[B,i,j]
            double *ptr = const_cast<double*>(buffer);
            C_DGEMM('n', 'n', nP, na*na, nQ, 1.0, ptr+0*stride, nQ, Bijp[oQ], na*na, 1.0, &deijp[Px][oP*na*na], na*na);
            C_DGEMM('n', 'n', nP, na*na, nQ, 1.0, ptr+1*stride, nQ, Bijp[oQ], na*na, 1.0, &deijp[Py][oP*na*na], na*na);
            C_DGEMM('n', 'n', nP, na*na, nQ, 1.0, ptr+2*stride, nQ, Bijp[oQ], na*na, 1.0, &deijp[Pz][oP*na*na], na*na);
            C_DGEMM('n', 'n', nP, na*na, nQ, 1.0, ptr+3*stride, nQ, Bijp[oQ], na*na, 1.0, &deijp[Qx][oP*na*na], na*na);
            C_DGEMM('n', 'n', nP, na*na, nQ, 1.0, ptr+4*stride, nQ, Bijp[oQ], na*na, 1.0, &deijp[Qy][oP*na*na], na*na);
            C_DGEMM('n', 'n', nP, na*na, nQ, 1.0, ptr+5*stride, nQ, Bijp[oQ], na*na, 1.0, &deijp[Qz][oP*na*na], na*na);

        }
    }


    for (int P = 0; P < nauxshell; ++P){
        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();
        int Pcenter = auxiliary_->shell(P).ncenter();
        int Pncart = auxiliary_->shell(P).ncartesian();
        int Px = 3 * Pcenter + 0;
        int Py = 3 * Pcenter + 1;
        int Pz = 3 * Pcenter + 2;
        for(int M = 0; M < nshell; ++M){
            int nM = primary_->shell(M).nfunction();
            int oM = primary_->shell(M).function_index();
            int Mcenter = primary_->shell(M).ncenter();
            int Mncart = primary_->shell(M).ncartesian();
            int mx = 3 * Mcenter + 0;
            int my = 3 * Mcenter + 1;
            int mz = 3 * Mcenter + 2;
            for(int N = 0; N < nshell; ++N){
                int nN = primary_->shell(N).nfunction();
                int oN = primary_->shell(N).function_index();
                int Ncenter = primary_->shell(N).ncenter();
                int Nncart = primary_->shell(N).ncartesian();
                int nx = 3 * Ncenter + 0;
                int ny = 3 * Ncenter + 1;
                int nz = 3 * Ncenter + 2;

                size_t stride = Pncart * Mncart * Nncart;

                Pmnint->compute_shell_deriv2(P,0,M,N);
                const double* buffer = Pmnint->buffer();

                double Pmscale = Pcenter == Mcenter ? 2.0 : 1.0;
                double Pnscale = Pcenter == Ncenter ? 2.0 : 1.0;
                double mnscale = Mcenter == Ncenter ? 2.0 : 1.0;

                double PxPx=0.0, PxPy=0.0, PxPz=0.0, PyPy=0.0, PyPz=0.0, PzPz=0.0;
                double mxmx=0.0, mxmy=0.0, mxmz=0.0, mymy=0.0, mymz=0.0, mzmz=0.0;
                double nxnx=0.0, nxny=0.0, nxnz=0.0, nyny=0.0, nynz=0.0, nznz=0.0;
                double Pxmx=0.0, Pxmy=0.0, Pxmz=0.0, Pymx=0.0, Pymy=0.0, Pymz=0.0, Pzmx=0.0, Pzmy=0.0, Pzmz=0.0;
                double Pxnx=0.0, Pxny=0.0, Pxnz=0.0, Pynx=0.0, Pyny=0.0, Pynz=0.0, Pznx=0.0, Pzny=0.0, Pznz=0.0;
                double mxnx=0.0, mxny=0.0, mxnz=0.0, mynx=0.0, myny=0.0, mynz=0.0, mznx=0.0, mzny=0.0, mznz=0.0;
                size_t delta = 0L;
                for (int p = oP; p < oP+nP; p++) {
                    for (int m = oM; m < oM+nM; m++) {
                        for (int n = oN; n < oN+nN; n++) {
                            double Cpmn = 2.0 * dp[p] * Dtp[m][n];
                            PxPx += Cpmn * buffer[9  * stride + delta];
                            PxPy += Cpmn * buffer[10 * stride + delta];
                            PxPz += Cpmn * buffer[11 * stride + delta];
                            Pxmx += Cpmn * buffer[12 * stride + delta];
                            Pxmy += Cpmn * buffer[13 * stride + delta];
                            Pxmz += Cpmn * buffer[14 * stride + delta];
                            Pxnx += Cpmn * buffer[15 * stride + delta];
                            Pxny += Cpmn * buffer[16 * stride + delta];
                            Pxnz += Cpmn * buffer[17 * stride + delta];
                            PyPy += Cpmn * buffer[18 * stride + delta];
                            PyPz += Cpmn * buffer[19 * stride + delta];
                            Pymx += Cpmn * buffer[20 * stride + delta];
                            Pymy += Cpmn * buffer[21 * stride + delta];
                            Pymz += Cpmn * buffer[22 * stride + delta];
                            Pynx += Cpmn * buffer[23 * stride + delta];
                            Pyny += Cpmn * buffer[24 * stride + delta];
                            Pynz += Cpmn * buffer[25 * stride + delta];
                            PzPz += Cpmn * buffer[26 * stride + delta];
                            Pzmx += Cpmn * buffer[27 * stride + delta];
                            Pzmy += Cpmn * buffer[28 * stride + delta];
                            Pzmz += Cpmn * buffer[29 * stride + delta];
                            Pznx += Cpmn * buffer[30 * stride + delta];
                            Pzny += Cpmn * buffer[31 * stride + delta];
                            Pznz += Cpmn * buffer[32 * stride + delta];
                            mxmx += Cpmn * buffer[33 * stride + delta];
                            mxmy += Cpmn * buffer[34 * stride + delta];
                            mxmz += Cpmn * buffer[35 * stride + delta];
                            mxnx += Cpmn * buffer[36 * stride + delta];
                            mxny += Cpmn * buffer[37 * stride + delta];
                            mxnz += Cpmn * buffer[38 * stride + delta];
                            mymy += Cpmn * buffer[39 * stride + delta];
                            mymz += Cpmn * buffer[40 * stride + delta];
                            mynx += Cpmn * buffer[41 * stride + delta];
                            myny += Cpmn * buffer[42 * stride + delta];
                            mynz += Cpmn * buffer[43 * stride + delta];
                            mzmz += Cpmn * buffer[44 * stride + delta];
                            mznx += Cpmn * buffer[45 * stride + delta];
                            mzny += Cpmn * buffer[46 * stride + delta];
                            mznz += Cpmn * buffer[47 * stride + delta];
                            nxnx += Cpmn * buffer[48 * stride + delta];
                            nxny += Cpmn * buffer[49 * stride + delta];
                            nxnz += Cpmn * buffer[50 * stride + delta];
                            nyny += Cpmn * buffer[51 * stride + delta];
                            nynz += Cpmn * buffer[52 * stride + delta];
                            nznz += Cpmn * buffer[53 * stride + delta];
                            ++delta;
                        }
                    }
                }
                JHessp[Px][Px] += PxPx;
                JHessp[Px][Py] += PxPy;
                JHessp[Px][Pz] += PxPz;
                JHessp[Px][mx] += Pmscale*Pxmx;
                JHessp[Px][my] += Pxmy;
                JHessp[Px][mz] += Pxmz;
                JHessp[Px][nx] += Pnscale*Pxnx;
                JHessp[Px][ny] += Pxny;
                JHessp[Px][nz] += Pxnz;
                JHessp[Py][Py] += PyPy;
                JHessp[Py][Pz] += PyPz;
                JHessp[Py][mx] += Pymx;
                JHessp[Py][my] += Pmscale*Pymy;
                JHessp[Py][mz] += Pymz;
                JHessp[Py][nx] += Pynx;
                JHessp[Py][ny] += Pnscale*Pyny;
                JHessp[Py][nz] += Pynz;
                JHessp[Pz][Pz] += PzPz;
                JHessp[Pz][mx] += Pzmx;
                JHessp[Pz][my] += Pzmy;
                JHessp[Pz][mz] += Pmscale*Pzmz;
                JHessp[Pz][nx] += Pznx;
                JHessp[Pz][ny] += Pzny;
                JHessp[Pz][nz] += Pnscale*Pznz;
                JHessp[mx][mx] += mxmx;
                JHessp[mx][my] += mxmy;
                JHessp[mx][mz] += mxmz;
                JHessp[mx][nx] += mnscale*mxnx;
                JHessp[mx][ny] += mxny;
                JHessp[mx][nz] += mxnz;
                JHessp[my][my] += mymy;
                JHessp[my][mz] += mymz;
                JHessp[my][nx] += mynx;
                JHessp[my][ny] += mnscale*myny;
                JHessp[my][nz] += mynz;
                JHessp[mz][mz] += mzmz;
                JHessp[mz][nx] += mznx;
                JHessp[mz][ny] += mzny;
                JHessp[mz][nz] += mnscale*mznz;
                JHessp[nx][nx] += nxnx;
                JHessp[nx][ny] += nxny;
                JHessp[nx][nz] += nxnz;
                JHessp[ny][ny] += nyny;
                JHessp[ny][nz] += nynz;
                JHessp[nz][nz] += nznz;


                // K terms
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
                            double Cpmn = 2.0 * Bmnp[p][m*nso+n];
                            PxPx += Cpmn * buffer[9  * stride + delta];
                            PxPy += Cpmn * buffer[10 * stride + delta];
                            PxPz += Cpmn * buffer[11 * stride + delta];
                            Pxmx += Cpmn * buffer[12 * stride + delta];
                            Pxmy += Cpmn * buffer[13 * stride + delta];
                            Pxmz += Cpmn * buffer[14 * stride + delta];
                            Pxnx += Cpmn * buffer[15 * stride + delta];
                            Pxny += Cpmn * buffer[16 * stride + delta];
                            Pxnz += Cpmn * buffer[17 * stride + delta];
                            PyPy += Cpmn * buffer[18 * stride + delta];
                            PyPz += Cpmn * buffer[19 * stride + delta];
                            Pymx += Cpmn * buffer[20 * stride + delta];
                            Pymy += Cpmn * buffer[21 * stride + delta];
                            Pymz += Cpmn * buffer[22 * stride + delta];
                            Pynx += Cpmn * buffer[23 * stride + delta];
                            Pyny += Cpmn * buffer[24 * stride + delta];
                            Pynz += Cpmn * buffer[25 * stride + delta];
                            PzPz += Cpmn * buffer[26 * stride + delta];
                            Pzmx += Cpmn * buffer[27 * stride + delta];
                            Pzmy += Cpmn * buffer[28 * stride + delta];
                            Pzmz += Cpmn * buffer[29 * stride + delta];
                            Pznx += Cpmn * buffer[30 * stride + delta];
                            Pzny += Cpmn * buffer[31 * stride + delta];
                            Pznz += Cpmn * buffer[32 * stride + delta];
                            mxmx += Cpmn * buffer[33 * stride + delta];
                            mxmy += Cpmn * buffer[34 * stride + delta];
                            mxmz += Cpmn * buffer[35 * stride + delta];
                            mxnx += Cpmn * buffer[36 * stride + delta];
                            mxny += Cpmn * buffer[37 * stride + delta];
                            mxnz += Cpmn * buffer[38 * stride + delta];
                            mymy += Cpmn * buffer[39 * stride + delta];
                            mymz += Cpmn * buffer[40 * stride + delta];
                            mynx += Cpmn * buffer[41 * stride + delta];
                            myny += Cpmn * buffer[42 * stride + delta];
                            mynz += Cpmn * buffer[43 * stride + delta];
                            mzmz += Cpmn * buffer[44 * stride + delta];
                            mznx += Cpmn * buffer[45 * stride + delta];
                            mzny += Cpmn * buffer[46 * stride + delta];
                            mznz += Cpmn * buffer[47 * stride + delta];
                            nxnx += Cpmn * buffer[48 * stride + delta];
                            nxny += Cpmn * buffer[49 * stride + delta];
                            nxnz += Cpmn * buffer[50 * stride + delta];
                            nyny += Cpmn * buffer[51 * stride + delta];
                            nynz += Cpmn * buffer[52 * stride + delta];
                            nznz += Cpmn * buffer[53 * stride + delta];
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

    for (int P = 0; P < nauxshell; ++P){
        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();
        int Pcenter = auxiliary_->shell(P).ncenter();
        int Pncart = auxiliary_->shell(P).ncartesian();
        int Px = 3 * Pcenter + 0;
        int Py = 3 * Pcenter + 1;
        int Pz = 3 * Pcenter + 2;
        for(int Q = 0; Q < nauxshell; ++Q){
            int nQ = auxiliary_->shell(Q).nfunction();
            int oQ = auxiliary_->shell(Q).function_index();
            int Qcenter = auxiliary_->shell(Q).ncenter();
            int Qncart = auxiliary_->shell(Q).ncartesian();
            int Qx = 3 * Qcenter + 0;
            int Qy = 3 * Qcenter + 1;
            int Qz = 3 * Qcenter + 2;

            size_t stride = Pncart * Qncart;

            PQint->compute_shell_deriv2(P,0,Q,0);
            const double* buffer = PQint->buffer();

            double PQscale = Pcenter == Qcenter ? 2.0 : 1.0;

            double PxPx=0.0, PxPy=0.0, PxPz=0.0, PyPy=0.0, PyPz=0.0, PzPz=0.0;
            double QxQx=0.0, QxQy=0.0, QxQz=0.0, QyQy=0.0, QyQz=0.0, QzQz=0.0;
            double PxQx=0.0, PxQy=0.0, PxQz=0.0, PyQx=0.0, PyQy=0.0, PyQz=0.0, PzQx=0.0, PzQy=0.0, PzQz=0.0;
            size_t delta = 0L;
            for (int p = oP; p < oP+nP; p++) {
                for (int q = oQ; q < oQ+nQ; q++) {
                    double dAdB = -dp[p]*dp[q];
                    PxPx += dAdB * buffer[9  * stride + delta];
                    PxPy += dAdB * buffer[10 * stride + delta];
                    PxPz += dAdB * buffer[11 * stride + delta];
                    PxQx += dAdB * buffer[12 * stride + delta];
                    PxQy += dAdB * buffer[13 * stride + delta];
                    PxQz += dAdB * buffer[14 * stride + delta];
                    PyPy += dAdB * buffer[18 * stride + delta];
                    PyPz += dAdB * buffer[19 * stride + delta];
                    PyQx += dAdB * buffer[20 * stride + delta];
                    PyQy += dAdB * buffer[21 * stride + delta];
                    PyQz += dAdB * buffer[22 * stride + delta];
                    PzPz += dAdB * buffer[26 * stride + delta];
                    PzQx += dAdB * buffer[27 * stride + delta];
                    PzQy += dAdB * buffer[28 * stride + delta];
                    PzQz += dAdB * buffer[29 * stride + delta];
                    QxQx += dAdB * buffer[33 * stride + delta];
                    QxQy += dAdB * buffer[34 * stride + delta];
                    QxQz += dAdB * buffer[35 * stride + delta];
                    QyQy += dAdB * buffer[39 * stride + delta];
                    QyQz += dAdB * buffer[40 * stride + delta];
                    QzQz += dAdB * buffer[44 * stride + delta];
                    ++delta;
                }

            }
            JHessp[Px][Px] += PxPx;
            JHessp[Px][Py] += PxPy;
            JHessp[Px][Pz] += PxPz;
            JHessp[Px][Qx] += PQscale*PxQx;
            JHessp[Px][Qy] += PxQy;
            JHessp[Px][Qz] += PxQz;
            JHessp[Py][Py] += PyPy;
            JHessp[Py][Pz] += PyPz;
            JHessp[Py][Qx] += PyQx;
            JHessp[Py][Qy] += PQscale*PyQy;
            JHessp[Py][Qz] += PyQz;
            JHessp[Pz][Pz] += PzPz;
            JHessp[Pz][Qx] += PzQx;
            JHessp[Pz][Qy] += PzQy;
            JHessp[Pz][Qz] += PQscale*PzQz;
            JHessp[Qx][Qx] += QxQx;
            JHessp[Qx][Qy] += QxQy;
            JHessp[Qx][Qz] += QxQz;
            JHessp[Qy][Qy] += QyQy;
            JHessp[Qy][Qz] += QyQz;
            JHessp[Qz][Qz] += QzQz;

            // K terms
            PxPx=0.0; PxPy=0.0; PxPz=0.0; PyPy=0.0; PyPz=0.0; PzPz=0.0;
            QxQx=0.0; QxQy=0.0; QxQz=0.0; QyQy=0.0; QyQz=0.0; QzQz=0.0;
            PxQx=0.0; PxQy=0.0; PxQz=0.0; PyQx=0.0; PyQy=0.0; PyQz=0.0; PzQx=0.0; PzQy=0.0; PzQz=0.0;
            delta = 0L;
            for (int p = oP; p < oP+nP; p++) {
                for (int q = oQ; q < oQ+nQ; q++) {
                    double dAdB = -DPQp[p][q];
                    PxPx += dAdB * buffer[9  * stride + delta];
                    PxPy += dAdB * buffer[10 * stride + delta];
                    PxPz += dAdB * buffer[11 * stride + delta];
                    PxQx += dAdB * buffer[12 * stride + delta];
                    PxQy += dAdB * buffer[13 * stride + delta];
                    PxQz += dAdB * buffer[14 * stride + delta];
                    PyPy += dAdB * buffer[18 * stride + delta];
                    PyPz += dAdB * buffer[19 * stride + delta];
                    PyQx += dAdB * buffer[20 * stride + delta];
                    PyQy += dAdB * buffer[21 * stride + delta];
                    PyQz += dAdB * buffer[22 * stride + delta];
                    PzPz += dAdB * buffer[26 * stride + delta];
                    PzQx += dAdB * buffer[27 * stride + delta];
                    PzQy += dAdB * buffer[28 * stride + delta];
                    PzQz += dAdB * buffer[29 * stride + delta];
                    QxQx += dAdB * buffer[33 * stride + delta];
                    QxQy += dAdB * buffer[34 * stride + delta];
                    QxQz += dAdB * buffer[35 * stride + delta];
                    QyQy += dAdB * buffer[39 * stride + delta];
                    QyQz += dAdB * buffer[40 * stride + delta];
                    QzQz += dAdB * buffer[44 * stride + delta];
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


    // Add permutational symmetry components missing from the above
    for(int i = 0; i < 3*natoms; ++i){
        for(int j = 0; j < i; ++j){
            JHessp[i][j] = JHessp[j][i] = (JHessp[i][j] + JHessp[j][i]);
            KHessp[i][j] = KHessp[j][i] = (KHessp[i][j] + KHessp[j][i]);
        }
    }


    // Stitch all the intermediates together to form the actual Hessian contributions
    SharedMatrix tmp(new Matrix("Tmp [P][i,j]", np, na*na));
    double **ptmp = tmp->pointer();

    for(int x = 0; x < 3*natoms; ++x){
        for(int y = 0; y < 3*natoms; ++y){
            // J terms
            JHessp[x][y] += 2.0*C_DDOT(np, ddp[x], 1, dcp[y], 1);
            JHessp[x][y] -= 4.0*C_DDOT(np, ddp[x], 1, dep[y], 1);
            C_DGEMV('n', np, np, 1.0, PQp[0], np, dep[y], 1, 0.0, ptmp[0], 1);
            JHessp[x][y] += 2.0*C_DDOT(np, dep[x], 1, ptmp[0], 1);

            // K terms
            C_DGEMM('n', 'n', np, na*na, np,  1.0, PQp[0], np, dAijp[y], na*na, 0.0, ptmp[0], na*na);
            KHessp[x][y] += 2.0*C_DDOT(np*na*na, dAijp[x], 1, ptmp[0], 1);
            C_DGEMM('n', 'n', np, na*na, np,  1.0, PQp[0], np, deijp[y], na*na, 0.0, ptmp[0], na*na);
            KHessp[x][y] -= 4.0*C_DDOT(np*na*na, dAijp[x], 1, ptmp[0], 1);
            KHessp[x][y] += 2.0*C_DDOT(np*na*na, deijp[x], 1, ptmp[0], 1);
        }
    }

    // Make sure the newly added components are symmetric
    for(int i = 0; i < 3*natoms; ++i){
        for(int j = 0; j < i; ++j){
            JHessp[i][j] = JHessp[j][i] = 0.5*(JHessp[i][j] + JHessp[j][i]);
            KHessp[i][j] = KHessp[j][i] = 0.5*(KHessp[i][j] + KHessp[j][i]);
        }
    }

    hessians_["Coulomb"]->scale(0.5);
}

DirectJKGrad::DirectJKGrad(int deriv, std::shared_ptr<BasisSet> primary) :
    JKGrad(deriv,primary)
{
    common_init();
}
DirectJKGrad::~DirectJKGrad()
{
}
void DirectJKGrad::common_init()
{
    ints_num_threads_ = 1;
#ifdef _OPENMP
    ints_num_threads_ = Process::environment.get_n_threads();
#endif
}
void DirectJKGrad::print_header() const
{
    if (print_) {
        outfile->Printf( "  ==> DirectJKGrad: Integral-Direct SCF Gradients <==\n\n");

        outfile->Printf( "    Gradient:          %11d\n", deriv_);
        outfile->Printf( "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf( "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf( "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_)
            outfile->Printf( "    Omega:             %11.3E\n", omega_);
        outfile->Printf( "    Integrals threads: %11d\n", ints_num_threads_);
        outfile->Printf( "    Schwarz Cutoff:    %11.0E\n", cutoff_);
        outfile->Printf( "\n");
    }
}
void DirectJKGrad::compute_gradient()
{
    if (!do_J_ && !do_K_ && !do_wK_)
        return;

    if (!(Ca_ && Cb_ && Da_ && Db_ && Dt_))
        throw PSIEXCEPTION("Occupation/Density not set");

    // => Set up gradients <= //
    int natom = primary_->molecule()->natom();
    gradients_.clear();
    if (do_J_) {
        gradients_["Coulomb"] = SharedMatrix(new Matrix("Coulomb Gradient",natom,3));
    }
    if (do_K_) {
        gradients_["Exchange"] = SharedMatrix(new Matrix("Exchange Gradient",natom,3));
    }
    if (do_wK_) {
        gradients_["Exchange,LR"] = SharedMatrix(new Matrix("Exchange,LR Gradient",natom,3));
    }

    // => Build ERI Sieve <= //
    sieve_ = std::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff_));

    std::shared_ptr<IntegralFactory> factory(new IntegralFactory(primary_,primary_,primary_,primary_));

    if (do_J_ || do_K_) {
        std::vector<std::shared_ptr<TwoBodyAOInt> > ints;
        for (int thread = 0; thread < ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri(1)));
        }
        std::map<std::string, std::shared_ptr<Matrix> > vals = compute1(ints);
        if (do_J_) {
            gradients_["Coulomb"]->copy(vals["J"]);
        }
        if (do_K_) {
            gradients_["Exchange"]->copy(vals["K"]);
        }
    }
    if (do_wK_) {
        std::vector<std::shared_ptr<TwoBodyAOInt> > ints;
        for (int thread = 0; thread < ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega_,1)));
        }
        std::map<std::string, std::shared_ptr<Matrix> > vals = compute1(ints);
        gradients_["Exchange,LR"]->copy(vals["K"]);
    }
}
std::map<std::string, std::shared_ptr<Matrix> > DirectJKGrad::compute1(std::vector<std::shared_ptr<TwoBodyAOInt> >& ints)
{
    int nthreads = ints.size();

    int natom = primary_->molecule()->natom();

    std::vector<std::shared_ptr<Matrix> > Jgrad;
    std::vector<std::shared_ptr<Matrix> > Kgrad;
    for (int thread = 0; thread < nthreads; thread++) {
        Jgrad.push_back(SharedMatrix(new Matrix("JGrad",natom,3)));
        Kgrad.push_back(SharedMatrix(new Matrix("KGrad",natom,3)));
    }

    const std::vector<std::pair<int, int> >& shell_pairs = sieve_->shell_pairs();
    size_t npairs = shell_pairs.size();
    size_t npairs2 = npairs * npairs;

    double** Dtp = Dt_->pointer();
    double** Dap = Da_->pointer();
    double** Dbp = Db_->pointer();

#pragma omp parallel for num_threads(nthreads) schedule(dynamic)
    for (size_t index = 0L; index < npairs2; index++) {

        size_t PQ = index / npairs;
        size_t RS = index % npairs;

        if (RS > PQ) continue;

        int P = shell_pairs[PQ].first;
        int Q = shell_pairs[PQ].second;
        int R = shell_pairs[RS].first;
        int S = shell_pairs[RS].second;

        if (!sieve_->shell_significant(P,Q,R,S)) continue;

        //outfile->Printf("(%d,%d,%d,%d)\n", P,Q,R,S);

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        ints[thread]->compute_shell_deriv1(P,Q,R,S);

        const double* buffer = ints[thread]->buffer();

        double** Jp = Jgrad[thread]->pointer();
        double** Kp = Kgrad[thread]->pointer();

        int Psize = primary_->shell(P).nfunction();
        int Qsize = primary_->shell(Q).nfunction();
        int Rsize = primary_->shell(R).nfunction();
        int Ssize = primary_->shell(S).nfunction();

        int Pncart = primary_->shell(P).ncartesian();
        int Qncart = primary_->shell(Q).ncartesian();
        int Rncart = primary_->shell(R).ncartesian();
        int Sncart = primary_->shell(S).ncartesian();

        int Poff = primary_->shell(P).function_index();
        int Qoff = primary_->shell(Q).function_index();
        int Roff = primary_->shell(R).function_index();
        int Soff = primary_->shell(S).function_index();

        int Pcenter = primary_->shell(P).ncenter();
        int Qcenter = primary_->shell(Q).ncenter();
        int Rcenter = primary_->shell(R).ncenter();
        int Scenter = primary_->shell(S).ncenter();

        double prefactor = 1.0;
        if (P != Q)   prefactor *= 2.0;
        if (R != S)   prefactor *= 2.0;
        if (PQ != RS) prefactor *= 2.0;

        size_t stride = Pncart * Qncart * Rncart * Sncart;

        double val;
        double Dpq, Drs;
        size_t delta;
        double Ax, Ay, Az;
        double Bx, By, Bz;
        double Cx, Cy, Cz;
        double Dx, Dy, Dz;

        // => Coulomb Term <= //

        Ax = 0.0; Ay = 0.0; Az = 0.0;
        Bx = 0.0; By = 0.0; Bz = 0.0;
        Cx = 0.0; Cy = 0.0; Cz = 0.0;
        Dx = 0.0; Dy = 0.0; Dz = 0.0;
        delta = 0L;
        for (int p = 0; p < Psize; p++) {
            for (int q = 0; q < Qsize; q++) {
                for (int r = 0; r < Rsize; r++) {
                    for (int s = 0; s < Ssize; s++) {
                        Dpq = Dtp[p + Poff][q + Qoff];
                        Drs = Dtp[r + Roff][s + Soff];
                        val = prefactor * Dpq * Drs;
                        Ax += val * buffer[0 * stride + delta];
                        Ay += val * buffer[1 * stride + delta];
                        Az += val * buffer[2 * stride + delta];
                        Cx += val * buffer[3 * stride + delta];
                        Cy += val * buffer[4 * stride + delta];
                        Cz += val * buffer[5 * stride + delta];
                        Dx += val * buffer[6 * stride + delta];
                        Dy += val * buffer[7 * stride + delta];
                        Dz += val * buffer[8 * stride + delta];
                        delta++;
                    }
                }
            }
        }
        Bx = -(Ax + Cx + Dx);
        By = -(Ay + Cy + Dy);
        Bz = -(Az + Cz + Dz);

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

        Ax = 0.0; Ay = 0.0; Az = 0.0;
        Bx = 0.0; By = 0.0; Bz = 0.0;
        Cx = 0.0; Cy = 0.0; Cz = 0.0;
        Dx = 0.0; Dy = 0.0; Dz = 0.0;
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
                        Ax += val * buffer[0 * stride + delta];
                        Ay += val * buffer[1 * stride + delta];
                        Az += val * buffer[2 * stride + delta];
                        Cx += val * buffer[3 * stride + delta];
                        Cy += val * buffer[4 * stride + delta];
                        Cz += val * buffer[5 * stride + delta];
                        Dx += val * buffer[6 * stride + delta];
                        Dy += val * buffer[7 * stride + delta];
                        Dz += val * buffer[8 * stride + delta];
                        delta++;
                    }
                }
            }
        }
        Bx = -(Ax + Cx + Dx);
        By = -(Ay + Cy + Dy);
        Bz = -(Az + Cz + Dz);

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

    }

    for (int thread = 1; thread < nthreads; thread++) {
        Jgrad[0]->add(Jgrad[thread]);
        Kgrad[0]->add(Kgrad[thread]);
    }

    Jgrad[0]->scale(0.5);
    Kgrad[0]->scale(0.5);

    std::map<std::string, std::shared_ptr<Matrix> > val;
    val["J"] = Jgrad[0];
    val["K"] = Kgrad[0];
    return val;
}
void DirectJKGrad::compute_hessian()
{
    if (!do_J_ && !do_K_ && !do_wK_)
        return;

    if (!(Ca_ && Cb_ && Da_ && Db_ && Dt_))
        throw PSIEXCEPTION("Occupation/Density not set");

    // => Set up hessians <= //
    int natom = primary_->molecule()->natom();
    hessians_.clear();
    if (do_J_) {
        hessians_["Coulomb"] = SharedMatrix(new Matrix("Coulomb Hessian",3*natom,3*natom));
    }
    if (do_K_) {
        hessians_["Exchange"] = SharedMatrix(new Matrix("Exchange Hessian",3*natom,3*natom));
    }
    if (do_wK_) {
        hessians_["Exchange,LR"] = SharedMatrix(new Matrix("Exchange,LR Hessian",3*natom,3*natom));
    }

    // => Build ERI Sieve <= //
    sieve_ = std::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff_));

    std::shared_ptr<IntegralFactory> factory(new IntegralFactory(primary_,primary_,primary_,primary_));

    if (do_J_ || do_K_) {
        std::vector<std::shared_ptr<TwoBodyAOInt> > ints;
        for (int thread = 0; thread < ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri(2)));
        }
        std::map<std::string, std::shared_ptr<Matrix> > vals = compute2(ints);
        if (do_J_) {
            hessians_["Coulomb"]->copy(vals["J"]);
        }
        if (do_K_) {
            hessians_["Exchange"]->copy(vals["K"]);
        }
    }
    if (do_wK_) {
        std::vector<std::shared_ptr<TwoBodyAOInt> > ints;
        for (int thread = 0; thread < ints_num_threads_; thread++) {
            ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega_,2)));
        }
        std::map<std::string, std::shared_ptr<Matrix> > vals = compute2(ints);
        hessians_["Exchange,LR"]->copy(vals["K"]);
    }
}
std::map<std::string, std::shared_ptr<Matrix> > DirectJKGrad::compute2(std::vector<std::shared_ptr<TwoBodyAOInt> >& ints)
{
    int nthreads = ints.size();

    int natom = primary_->molecule()->natom();

    std::vector<std::shared_ptr<Matrix> > Jhess;
    std::vector<std::shared_ptr<Matrix> > Khess;
    for (int thread = 0; thread < nthreads; thread++) {
        Jhess.push_back(SharedMatrix(new Matrix("JHess",3*natom,3*natom)));
        Khess.push_back(SharedMatrix(new Matrix("KHess",3*natom,3*natom)));
    }

    const std::vector<std::pair<int, int> >& shell_pairs = sieve_->shell_pairs();
    size_t npairs = shell_pairs.size();
    size_t npairs2 = npairs * npairs;

    double** Dtp = Dt_->pointer();
    double** Dap = Da_->pointer();
    double** Dbp = Db_->pointer();

#pragma omp parallel for num_threads(nthreads) schedule(dynamic)
    for (size_t index = 0L; index < npairs2; index++) {

        size_t PQ = index / npairs;
        size_t RS = index % npairs;

        if (RS > PQ) continue;

        int P = shell_pairs[PQ].first;
        int Q = shell_pairs[PQ].second;
        int R = shell_pairs[RS].first;
        int S = shell_pairs[RS].second;

        if (!sieve_->shell_significant(P,Q,R,S)) continue;

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        ints[thread]->compute_shell_deriv2(P,Q,R,S);

        const double* buffer = ints[thread]->buffer();
        double** Jp = Jhess[thread]->pointer();
        double** Kp = Khess[thread]->pointer();


        int Psize = primary_->shell(P).nfunction();
        int Qsize = primary_->shell(Q).nfunction();
        int Rsize = primary_->shell(R).nfunction();
        int Ssize = primary_->shell(S).nfunction();

        int Pncart = primary_->shell(P).ncartesian();
        int Qncart = primary_->shell(Q).ncartesian();
        int Rncart = primary_->shell(R).ncartesian();
        int Sncart = primary_->shell(S).ncartesian();

        int Poff = primary_->shell(P).function_index();
        int Qoff = primary_->shell(Q).function_index();
        int Roff = primary_->shell(R).function_index();
        int Soff = primary_->shell(S).function_index();

        int Pcenter = primary_->shell(P).ncenter();
        int Qcenter = primary_->shell(Q).ncenter();
        int Rcenter = primary_->shell(R).ncenter();
        int Scenter = primary_->shell(S).ncenter();

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

        double prefactor = 4.0;
        if (P == Q)   prefactor *= 0.5;
        if (R == S)   prefactor *= 0.5;
        if (PQ == RS) prefactor *= 0.5;

        size_t stride = Pncart * Qncart * Rncart * Sncart;

        double val;
        double Dpq, Drs;
        size_t delta;

        // => Coulomb Term <= //

        double AxAx=0.0, AxAy=0.0, AxAz=0.0, AyAy=0.0, AyAz=0.0, AzAz=0.0;
        double BxBx=0.0, BxBy=0.0, BxBz=0.0, ByBy=0.0, ByBz=0.0, BzBz=0.0;
        double CxCx=0.0, CxCy=0.0, CxCz=0.0, CyCy=0.0, CyCz=0.0, CzCz=0.0;
        double DxDx=0.0, DxDy=0.0, DxDz=0.0, DyDy=0.0, DyDz=0.0, DzDz=0.0;
        double AxBx=0.0, AxBy=0.0, AxBz=0.0, AyBx=0.0, AyBy=0.0, AyBz=0.0, AzBx=0.0, AzBy=0.0, AzBz=0.0;
        double AxCx=0.0, AxCy=0.0, AxCz=0.0, AyCx=0.0, AyCy=0.0, AyCz=0.0, AzCx=0.0, AzCy=0.0, AzCz=0.0;
        double AxDx=0.0, AxDy=0.0, AxDz=0.0, AyDx=0.0, AyDy=0.0, AyDz=0.0, AzDx=0.0, AzDy=0.0, AzDz=0.0;
        double BxCx=0.0, BxCy=0.0, BxCz=0.0, ByCx=0.0, ByCy=0.0, ByCz=0.0, BzCx=0.0, BzCy=0.0, BzCz=0.0;
        double BxDx=0.0, BxDy=0.0, BxDz=0.0, ByDx=0.0, ByDy=0.0, ByDz=0.0, BzDx=0.0, BzDy=0.0, BzDz=0.0;
        double CxDx=0.0, CxDy=0.0, CxDz=0.0, CyDx=0.0, CyDy=0.0, CyDz=0.0, CzDx=0.0, CzDy=0.0, CzDz=0.0;

        delta = 0L;
        for (int p = 0; p < Psize; p++) {
            for (int q = 0; q < Qsize; q++) {
                for (int r = 0; r < Rsize; r++) {
                    for (int s = 0; s < Ssize; s++) {
                        Dpq = Dtp[p + Poff][q + Qoff];
                        Drs = Dtp[r + Roff][s + Soff];
                        val = prefactor * Dpq * Drs;
                        AxAx += val * buffer[9  * stride + delta];
                        AxAy += val * buffer[10 * stride + delta];
                        AxAz += val * buffer[11 * stride + delta];
                        AxCx += val * buffer[12 * stride + delta];
                        AxCy += val * buffer[13 * stride + delta];
                        AxCz += val * buffer[14 * stride + delta];
                        AxDx += val * buffer[15 * stride + delta];
                        AxDy += val * buffer[16 * stride + delta];
                        AxDz += val * buffer[17 * stride + delta];
                        AyAy += val * buffer[18 * stride + delta];
                        AyAz += val * buffer[19 * stride + delta];
                        AyCx += val * buffer[20 * stride + delta];
                        AyCy += val * buffer[21 * stride + delta];
                        AyCz += val * buffer[22 * stride + delta];
                        AyDx += val * buffer[23 * stride + delta];
                        AyDy += val * buffer[24 * stride + delta];
                        AyDz += val * buffer[25 * stride + delta];
                        AzAz += val * buffer[26 * stride + delta];
                        AzCx += val * buffer[27 * stride + delta];
                        AzCy += val * buffer[28 * stride + delta];
                        AzCz += val * buffer[29 * stride + delta];
                        AzDx += val * buffer[30 * stride + delta];
                        AzDy += val * buffer[31 * stride + delta];
                        AzDz += val * buffer[32 * stride + delta];
                        CxCx += val * buffer[33 * stride + delta];
                        CxCy += val * buffer[34 * stride + delta];
                        CxCz += val * buffer[35 * stride + delta];
                        CxDx += val * buffer[36 * stride + delta];
                        CxDy += val * buffer[37 * stride + delta];
                        CxDz += val * buffer[38 * stride + delta];
                        CyCy += val * buffer[39 * stride + delta];
                        CyCz += val * buffer[40 * stride + delta];
                        CyDx += val * buffer[41 * stride + delta];
                        CyDy += val * buffer[42 * stride + delta];
                        CyDz += val * buffer[43 * stride + delta];
                        CzCz += val * buffer[44 * stride + delta];
                        CzDx += val * buffer[45 * stride + delta];
                        CzDy += val * buffer[46 * stride + delta];
                        CzDz += val * buffer[47 * stride + delta];
                        DxDx += val * buffer[48 * stride + delta];
                        DxDy += val * buffer[49 * stride + delta];
                        DxDz += val * buffer[50 * stride + delta];
                        DyDy += val * buffer[51 * stride + delta];
                        DyDz += val * buffer[52 * stride + delta];
                        DzDz += val * buffer[53 * stride + delta];
                        delta++;
                    }
                }
            }
        }

        // Translational invariance relationships
        AxBx = -(AxAx + AxCx + AxDx);
        AxBy = -(AxAy + AxCy + AxDy);
        AxBz = -(AxAz + AxCz + AxDz);
        AyBx = -(AxAy + AyCx + AyDx);
        AyBy = -(AyAy + AyCy + AyDy);
        AyBz = -(AyAz + AyCz + AyDz);
        AzBx = -(AxAz + AzCx + AzDx);
        AzBy = -(AyAz + AzCy + AzDy);
        AzBz = -(AzAz + AzCz + AzDz);
        BxCx = -(AxCx + CxCx + CxDx);
        BxCy = -(AxCy + CxCy + CyDx);
        BxCz = -(AxCz + CxCz + CzDx);
        ByCx = -(AyCx + CxCy + CxDy);
        ByCy = -(AyCy + CyCy + CyDy);
        ByCz = -(AyCz + CyCz + CzDy);
        BzCx = -(AzCx + CxCz + CxDz);
        BzCy = -(AzCy + CyCz + CyDz);
        BzCz = -(AzCz + CzCz + CzDz);
        BxDx = -(AxDx + CxDx + DxDx);
        BxDy = -(AxDy + CxDy + DxDy);
        BxDz = -(AxDz + CxDz + DxDz);
        ByDx = -(AyDx + CyDx + DxDy);
        ByDy = -(AyDy + CyDy + DyDy);
        ByDz = -(AyDz + CyDz + DyDz);
        BzDx = -(AzDx + CzDx + DxDz);
        BzDy = -(AzDy + CzDy + DyDz);
        BzDz = -(AzDz + CzDz + DzDz);

        BxBx = AxAx + AxCx + AxDx
                + AxCx + CxCx + CxDx
                + AxDx + CxDx + DxDx;
        ByBy = AyAy + AyCy + AyDy
                + AyCy + CyCy + CyDy
                + AyDy + CyDy + DyDy;
        BzBz = AzAz + AzCz + AzDz
                + AzCz + CzCz + CzDz
                + AzDz + CzDz + DzDz;
        BxBy = AxAy + AxCy + AxDy
                + AyCx + CxCy + CxDy
                + AyDx + CyDx + DxDy;
        BxBz = AxAz + AxCz + AxDz
                + AzCx + CxCz + CxDz
                + AzDx + CzDx + DxDz;
        ByBz = AyAz + AyCz + AyDz
                + AzCy + CyCz + CyDz
                + AzDy + CzDy + DyDz;

        Jp[Px][Px] += AxAx;
        Jp[Px][Py] += AxAy;
        Jp[Px][Pz] += AxAz;
        Jp[Px][Qx] += PQscale*AxBx;
        Jp[Px][Qy] += AxBy;
        Jp[Px][Qz] += AxBz;
        Jp[Px][Rx] += PRscale*AxCx;
        Jp[Px][Ry] += AxCy;
        Jp[Px][Rz] += AxCz;
        Jp[Px][Sx] += PSscale*AxDx;
        Jp[Px][Sy] += AxDy;
        Jp[Px][Sz] += AxDz;
        Jp[Py][Py] += AyAy;
        Jp[Py][Pz] += AyAz;
        Jp[Py][Qx] += AyBx;
        Jp[Py][Qy] += PQscale*AyBy;
        Jp[Py][Qz] += AyBz;
        Jp[Py][Rx] += AyCx;
        Jp[Py][Ry] += PRscale*AyCy;
        Jp[Py][Rz] += AyCz;
        Jp[Py][Sx] += AyDx;
        Jp[Py][Sy] += PSscale*AyDy;
        Jp[Py][Sz] += AyDz;
        Jp[Pz][Pz] += AzAz;
        Jp[Pz][Qx] += AzBx;
        Jp[Pz][Qy] += AzBy;
        Jp[Pz][Qz] += PQscale*AzBz;
        Jp[Pz][Rx] += AzCx;
        Jp[Pz][Ry] += AzCy;
        Jp[Pz][Rz] += PRscale*AzCz;
        Jp[Pz][Sx] += AzDx;
        Jp[Pz][Sy] += AzDy;
        Jp[Pz][Sz] += PSscale*AzDz;
        Jp[Qx][Qx] += BxBx;
        Jp[Qx][Qy] += BxBy;
        Jp[Qx][Qz] += BxBz;
        Jp[Qx][Rx] += QRscale*BxCx;
        Jp[Qx][Ry] += BxCy;
        Jp[Qx][Rz] += BxCz;
        Jp[Qx][Sx] += QSscale*BxDx;
        Jp[Qx][Sy] += BxDy;
        Jp[Qx][Sz] += BxDz;
        Jp[Qy][Qy] += ByBy;
        Jp[Qy][Qz] += ByBz;
        Jp[Qy][Rx] += ByCx;
        Jp[Qy][Ry] += QRscale*ByCy;
        Jp[Qy][Rz] += ByCz;
        Jp[Qy][Sx] += ByDx;
        Jp[Qy][Sy] += QSscale*ByDy;
        Jp[Qy][Sz] += ByDz;
        Jp[Qz][Qz] += BzBz;
        Jp[Qz][Rx] += BzCx;
        Jp[Qz][Ry] += BzCy;
        Jp[Qz][Rz] += QRscale*BzCz;
        Jp[Qz][Sx] += BzDx;
        Jp[Qz][Sy] += BzDy;
        Jp[Qz][Sz] += QSscale*BzDz;
        Jp[Rx][Rx] += CxCx;
        Jp[Rx][Ry] += CxCy;
        Jp[Rx][Rz] += CxCz;
        Jp[Rx][Sx] += RSscale*CxDx;
        Jp[Rx][Sy] += CxDy;
        Jp[Rx][Sz] += CxDz;
        Jp[Ry][Ry] += CyCy;
        Jp[Ry][Rz] += CyCz;
        Jp[Ry][Sx] += CyDx;
        Jp[Ry][Sy] += RSscale*CyDy;
        Jp[Ry][Sz] += CyDz;
        Jp[Rz][Rz] += CzCz;
        Jp[Rz][Sx] += CzDx;
        Jp[Rz][Sy] += CzDy;
        Jp[Rz][Sz] += RSscale*CzDz;
        Jp[Sx][Sx] += DxDx;
        Jp[Sx][Sy] += DxDy;
        Jp[Sx][Sz] += DxDz;
        Jp[Sy][Sy] += DyDy;
        Jp[Sy][Sz] += DyDz;
        Jp[Sz][Sz] += DzDz;

        // => Exchange Term <= //

        AxAx=0.0; AxAy=0.0; AxAz=0.0; AyAy=0.0; AyAz=0.0; AzAz=0.0;
        BxBx=0.0; BxBy=0.0; BxBz=0.0; ByBy=0.0; ByBz=0.0; BzBz=0.0;
        CxCx=0.0; CxCy=0.0; CxCz=0.0; CyCy=0.0; CyCz=0.0; CzCz=0.0;
        DxDx=0.0; DxDy=0.0; DxDz=0.0; DyDy=0.0; DyDz=0.0; DzDz=0.0;
        AxBx=0.0; AxBy=0.0; AxBz=0.0; AyBx=0.0; AyBy=0.0; AyBz=0.0; AzBx=0.0; AzBy=0.0; AzBz=0.0;
        AxCx=0.0; AxCy=0.0; AxCz=0.0; AyCx=0.0; AyCy=0.0; AyCz=0.0; AzCx=0.0; AzCy=0.0; AzCz=0.0;
        AxDx=0.0; AxDy=0.0; AxDz=0.0; AyDx=0.0; AyDy=0.0; AyDz=0.0; AzDx=0.0; AzDy=0.0; AzDz=0.0;
        BxCx=0.0; BxCy=0.0; BxCz=0.0; ByCx=0.0; ByCy=0.0; ByCz=0.0; BzCx=0.0; BzCy=0.0; BzCz=0.0;
        BxDx=0.0; BxDy=0.0; BxDz=0.0; ByDx=0.0; ByDy=0.0; ByDz=0.0; BzDx=0.0; BzDy=0.0; BzDz=0.0;
        CxDx=0.0; CxDy=0.0; CxDz=0.0; CyDx=0.0; CyDy=0.0; CyDz=0.0; CzDx=0.0; CzDy=0.0; CzDz=0.0;


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
                        AxAx += val * buffer[9 * stride + delta];
                        AxAy += val * buffer[10 * stride + delta];
                        AxAz += val * buffer[11 * stride + delta];
                        AxCx += val * buffer[12 * stride + delta];
                        AxCy += val * buffer[13 * stride + delta];
                        AxCz += val * buffer[14 * stride + delta];
                        AxDx += val * buffer[15 * stride + delta];
                        AxDy += val * buffer[16 * stride + delta];
                        AxDz += val * buffer[17 * stride + delta];
                        AyAy += val * buffer[18 * stride + delta];
                        AyAz += val * buffer[19 * stride + delta];
                        AyCx += val * buffer[20 * stride + delta];
                        AyCy += val * buffer[21 * stride + delta];
                        AyCz += val * buffer[22 * stride + delta];
                        AyDx += val * buffer[23 * stride + delta];
                        AyDy += val * buffer[24 * stride + delta];
                        AyDz += val * buffer[25 * stride + delta];
                        AzAz += val * buffer[26 * stride + delta];
                        AzCx += val * buffer[27 * stride + delta];
                        AzCy += val * buffer[28 * stride + delta];
                        AzCz += val * buffer[29 * stride + delta];
                        AzDx += val * buffer[30 * stride + delta];
                        AzDy += val * buffer[31 * stride + delta];
                        AzDz += val * buffer[32 * stride + delta];
                        CxCx += val * buffer[33 * stride + delta];
                        CxCy += val * buffer[34 * stride + delta];
                        CxCz += val * buffer[35 * stride + delta];
                        CxDx += val * buffer[36 * stride + delta];
                        CxDy += val * buffer[37 * stride + delta];
                        CxDz += val * buffer[38 * stride + delta];
                        CyCy += val * buffer[39 * stride + delta];
                        CyCz += val * buffer[40 * stride + delta];
                        CyDx += val * buffer[41 * stride + delta];
                        CyDy += val * buffer[42 * stride + delta];
                        CyDz += val * buffer[43 * stride + delta];
                        CzCz += val * buffer[44 * stride + delta];
                        CzDx += val * buffer[45 * stride + delta];
                        CzDy += val * buffer[46 * stride + delta];
                        CzDz += val * buffer[47 * stride + delta];
                        DxDx += val * buffer[48 * stride + delta];
                        DxDy += val * buffer[49 * stride + delta];
                        DxDz += val * buffer[50 * stride + delta];
                        DyDy += val * buffer[51 * stride + delta];
                        DyDz += val * buffer[52 * stride + delta];
                        DzDz += val * buffer[53 * stride + delta];
                        delta++;
                    }
                }
            }
        }

        // Translational invariance relationships
        AxBx = -(AxAx + AxCx + AxDx);
        AxBy = -(AxAy + AxCy + AxDy);
        AxBz = -(AxAz + AxCz + AxDz);
        AyBx = -(AxAy + AyCx + AyDx);
        AyBy = -(AyAy + AyCy + AyDy);
        AyBz = -(AyAz + AyCz + AyDz);
        AzBx = -(AxAz + AzCx + AzDx);
        AzBy = -(AyAz + AzCy + AzDy);
        AzBz = -(AzAz + AzCz + AzDz);
        BxCx = -(AxCx + CxCx + CxDx);
        BxCy = -(AxCy + CxCy + CyDx);
        BxCz = -(AxCz + CxCz + CzDx);
        ByCx = -(AyCx + CxCy + CxDy);
        ByCy = -(AyCy + CyCy + CyDy);
        ByCz = -(AyCz + CyCz + CzDy);
        BzCx = -(AzCx + CxCz + CxDz);
        BzCy = -(AzCy + CyCz + CyDz);
        BzCz = -(AzCz + CzCz + CzDz);
        BxDx = -(AxDx + CxDx + DxDx);
        BxDy = -(AxDy + CxDy + DxDy);
        BxDz = -(AxDz + CxDz + DxDz);
        ByDx = -(AyDx + CyDx + DxDy);
        ByDy = -(AyDy + CyDy + DyDy);
        ByDz = -(AyDz + CyDz + DyDz);
        BzDx = -(AzDx + CzDx + DxDz);
        BzDy = -(AzDy + CzDy + DyDz);
        BzDz = -(AzDz + CzDz + DzDz);

        BxBx = AxAx + AxCx + AxDx
                + AxCx + CxCx + CxDx
                + AxDx + CxDx + DxDx;
        ByBy = AyAy + AyCy + AyDy
                + AyCy + CyCy + CyDy
                + AyDy + CyDy + DyDy;
        BzBz = AzAz + AzCz + AzDz
                + AzCz + CzCz + CzDz
                + AzDz + CzDz + DzDz;
        BxBy = AxAy + AxCy + AxDy
                + AyCx + CxCy + CxDy
                + AyDx + CyDx + DxDy;
        BxBz = AxAz + AxCz + AxDz
                + AzCx + CxCz + CxDz
                + AzDx + CzDx + DxDz;
        ByBz = AyAz + AyCz + AyDz
                + AzCy + CyCz + CyDz
                + AzDy + CzDy + DyDz;

        Kp[Px][Px] += AxAx;
        Kp[Px][Py] += AxAy;
        Kp[Px][Pz] += AxAz;
        Kp[Px][Qx] += PQscale*AxBx;
        Kp[Px][Qy] += AxBy;
        Kp[Px][Qz] += AxBz;
        Kp[Px][Rx] += PRscale*AxCx;
        Kp[Px][Ry] += AxCy;
        Kp[Px][Rz] += AxCz;
        Kp[Px][Sx] += PSscale*AxDx;
        Kp[Px][Sy] += AxDy;
        Kp[Px][Sz] += AxDz;
        Kp[Py][Py] += AyAy;
        Kp[Py][Pz] += AyAz;
        Kp[Py][Qx] += AyBx;
        Kp[Py][Qy] += PQscale*AyBy;
        Kp[Py][Qz] += AyBz;
        Kp[Py][Rx] += AyCx;
        Kp[Py][Ry] += PRscale*AyCy;
        Kp[Py][Rz] += AyCz;
        Kp[Py][Sx] += AyDx;
        Kp[Py][Sy] += PSscale*AyDy;
        Kp[Py][Sz] += AyDz;
        Kp[Pz][Pz] += AzAz;
        Kp[Pz][Qx] += AzBx;
        Kp[Pz][Qy] += AzBy;
        Kp[Pz][Qz] += PQscale*AzBz;
        Kp[Pz][Rx] += AzCx;
        Kp[Pz][Ry] += AzCy;
        Kp[Pz][Rz] += PRscale*AzCz;
        Kp[Pz][Sx] += AzDx;
        Kp[Pz][Sy] += AzDy;
        Kp[Pz][Sz] += PSscale*AzDz;
        Kp[Qx][Qx] += BxBx;
        Kp[Qx][Qy] += BxBy;
        Kp[Qx][Qz] += BxBz;
        Kp[Qx][Rx] += QRscale*BxCx;
        Kp[Qx][Ry] += BxCy;
        Kp[Qx][Rz] += BxCz;
        Kp[Qx][Sx] += QSscale*BxDx;
        Kp[Qx][Sy] += BxDy;
        Kp[Qx][Sz] += BxDz;
        Kp[Qy][Qy] += ByBy;
        Kp[Qy][Qz] += ByBz;
        Kp[Qy][Rx] += ByCx;
        Kp[Qy][Ry] += QRscale*ByCy;
        Kp[Qy][Rz] += ByCz;
        Kp[Qy][Sx] += ByDx;
        Kp[Qy][Sy] += QSscale*ByDy;
        Kp[Qy][Sz] += ByDz;
        Kp[Qz][Qz] += BzBz;
        Kp[Qz][Rx] += BzCx;
        Kp[Qz][Ry] += BzCy;
        Kp[Qz][Rz] += QRscale*BzCz;
        Kp[Qz][Sx] += BzDx;
        Kp[Qz][Sy] += BzDy;
        Kp[Qz][Sz] += QSscale*BzDz;
        Kp[Rx][Rx] += CxCx;
        Kp[Rx][Ry] += CxCy;
        Kp[Rx][Rz] += CxCz;
        Kp[Rx][Sx] += RSscale*CxDx;
        Kp[Rx][Sy] += CxDy;
        Kp[Rx][Sz] += CxDz;
        Kp[Ry][Ry] += CyCy;
        Kp[Ry][Rz] += CyCz;
        Kp[Ry][Sx] += CyDx;
        Kp[Ry][Sy] += RSscale*CyDy;
        Kp[Ry][Sz] += CyDz;
        Kp[Rz][Rz] += CzCz;
        Kp[Rz][Sx] += CzDx;
        Kp[Rz][Sy] += CzDy;
        Kp[Rz][Sz] += RSscale*CzDz;
        Kp[Sx][Sx] += DxDx;
        Kp[Sx][Sy] += DxDy;
        Kp[Sx][Sz] += DxDz;
        Kp[Sy][Sy] += DyDy;
        Kp[Sy][Sz] += DyDz;
        Kp[Sz][Sz] += DzDz;
    }

    for (int thread = 1; thread < nthreads; thread++) {
        Jhess[0]->add(Jhess[thread]);
        Khess[0]->add(Khess[thread]);
    }
    int dim = Jhess[0]->rowdim();
    double **Jp = Jhess[0]->pointer();
    double **Kp = Khess[0]->pointer();
    for (int row = 0; row < dim; ++row){
        for (int col = 0; col < row; ++col){
            Jp[row][col] = Jp[col][row] = (Jp[row][col] + Jp[col][row]);
            Kp[row][col] = Kp[col][row] = (Kp[row][col] + Kp[col][row]);
        }
    }

    std::map<std::string, std::shared_ptr<Matrix> > val;
    val["J"] = Jhess[0];
    val["K"] = Khess[0];
    return val;
}

}} // Namespaces
