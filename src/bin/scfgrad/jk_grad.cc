#include <libmints/mints.h>
#include <libmints/sieve.h>
#include <libqt/qt.h>
#include <lib3index/3index.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include <psifiles.h>
#include "jk_grad.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;

namespace psi {
namespace scfgrad {

JKGrad::JKGrad(int deriv, boost::shared_ptr<BasisSet> primary) :
    deriv_(deriv), primary_(primary)
{
    common_init();
}
JKGrad::~JKGrad()
{
}
boost::shared_ptr<JKGrad> JKGrad::build_JKGrad(int deriv)
{
    Options& options = Process::environment.options;
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    boost::shared_ptr<BasisSet> primary = BasisSet::construct(parser, Process::environment.molecule(), "BASIS");

    if (options.get_str("SCF_TYPE") == "DF") {

        boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, primary->molecule(), "DF_BASIS_SCF");

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

        return boost::shared_ptr<JKGrad>(jk);
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

        return boost::shared_ptr<JKGrad>(jk);

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
        omp_num_threads_ = omp_get_max_threads();
    #endif

    cutoff_ = 0.0;

    do_J_ = true;
    do_K_ = true;
    do_wK_ = false;
    omega_ = 0.0;
}
DFJKGrad::DFJKGrad(int deriv, boost::shared_ptr<BasisSet> primary, boost::shared_ptr<BasisSet> auxiliary) :
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
        df_ints_num_threads_ = omp_get_max_threads();
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
        fprintf(outfile, "  ==> DFJKGrad: Density-Fitted SCF Gradients <==\n\n");

        fprintf(outfile, "    Gradient:          %11d\n", deriv_);
        fprintf(outfile, "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        fprintf(outfile, "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        fprintf(outfile, "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_)
            fprintf(outfile, "    Omega:             %11.3E\n", omega_);
        fprintf(outfile, "    OpenMP threads:    %11d\n", omp_num_threads_);
        fprintf(outfile, "    Integrals threads: %11d\n", df_ints_num_threads_);
        fprintf(outfile, "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        fprintf(outfile, "    Schwarz Cutoff:    %11.0E\n", cutoff_);
        fprintf(outfile, "    Fitting Condition: %11.0E\n\n", condition_);

        fprintf(outfile, "   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level(outfile, print_);
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
        gradients_["Exchange,LR"] = SharedMatrix(new Matrix("Exchange,LR Gradient",natom,3));
    }

    // => Build ERI Sieve <= //
    sieve_ = boost::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff_));

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

    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->eri()));
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

    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_)));
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

    boost::shared_ptr<FittingMetric> metric(new FittingMetric(auxiliary_, true));
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

    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_,BasisSet::zero_ao_basis_set(),auxiliary_,BasisSet::zero_ao_basis_set()));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > Jint;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        Jint.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)));
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

    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->eri(1)));
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

    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > eri;
    for (int t = 0; t < df_ints_num_threads_; t++) {
        eri.push_back(boost::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_,1)));
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
    throw PSIEXCEPTION("Not implemented");
}

DirectJKGrad::DirectJKGrad(int deriv, boost::shared_ptr<BasisSet> primary) :
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
        ints_num_threads_ = omp_get_max_threads();
    #endif
}
void DirectJKGrad::print_header() const
{
    if (print_) {
        fprintf(outfile, "  ==> DirectJKGrad: Integral-Direct SCF Gradients <==\n\n");

        fprintf(outfile, "    Gradient:          %11d\n", deriv_);
        fprintf(outfile, "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        fprintf(outfile, "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        fprintf(outfile, "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_)
            fprintf(outfile, "    Omega:             %11.3E\n", omega_);
        fprintf(outfile, "    Integrals threads: %11d\n", ints_num_threads_);
        fprintf(outfile, "    Schwarz Cutoff:    %11.0E\n", cutoff_);
        fprintf(outfile, "\n");
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
    sieve_ = boost::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff_));

    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(primary_,primary_,primary_,primary_));

    if (do_J_ || do_K_) {
        std::vector<boost::shared_ptr<TwoBodyAOInt> > ints;
        for (int thread = 0; thread < ints_num_threads_; thread++) {
            ints.push_back(boost::shared_ptr<TwoBodyAOInt>(factory->eri(1)));
        }
        std::map<std::string, boost::shared_ptr<Matrix> > vals = compute1(ints);
        if (do_J_) {
            gradients_["Coulomb"]->copy(vals["J"]);
        }
        if (do_K_) {
            gradients_["Exchange"]->copy(vals["K"]);
        }
    }
    if (do_wK_) {
        std::vector<boost::shared_ptr<TwoBodyAOInt> > ints;
        for (int thread = 0; thread < ints_num_threads_; thread++) {
            ints.push_back(boost::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega_,1)));
        }
        std::map<std::string, boost::shared_ptr<Matrix> > vals = compute1(ints);
        gradients_["Exchange,LR"]->copy(vals["K"]);
    }
}
std::map<std::string, boost::shared_ptr<Matrix> > DirectJKGrad::compute1(std::vector<boost::shared_ptr<TwoBodyAOInt> >& ints)
{
    int nthreads = ints.size();

    int natom = primary_->molecule()->natom();

    std::vector<boost::shared_ptr<Matrix> > Jgrad;
    std::vector<boost::shared_ptr<Matrix> > Kgrad;
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

        //fprintf(outfile,"(%d,%d,%d,%d)\n", P,Q,R,S);

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
        }}}}
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
        }}}}
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

    std::map<std::string, boost::shared_ptr<Matrix> > val;
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
    sieve_ = boost::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff_));

    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(primary_,primary_,primary_,primary_));

    if (do_J_ || do_K_) {
        std::vector<boost::shared_ptr<TwoBodyAOInt> > ints;
        for (int thread = 0; thread < ints_num_threads_; thread++) {
            ints.push_back(boost::shared_ptr<TwoBodyAOInt>(factory->eri(2)));
        }
        std::map<std::string, boost::shared_ptr<Matrix> > vals = compute2(ints);
        if (do_J_) {
            hessians_["Coulomb"]->copy(vals["J"]);
        }
        if (do_K_) {
            hessians_["Exchange"]->copy(vals["K"]);
        }
    }
    if (do_wK_) {
        std::vector<boost::shared_ptr<TwoBodyAOInt> > ints;
        for (int thread = 0; thread < ints_num_threads_; thread++) {
            ints.push_back(boost::shared_ptr<TwoBodyAOInt>(factory->erf_eri(omega_,2)));
        }
        std::map<std::string, boost::shared_ptr<Matrix> > vals = compute2(ints);
        hessians_["Exchange,LR"]->copy(vals["K"]);
    }
}
std::map<std::string, boost::shared_ptr<Matrix> > DirectJKGrad::compute2(std::vector<boost::shared_ptr<TwoBodyAOInt> >& ints)
{
    int nthreads = ints.size();

    int natom = primary_->molecule()->natom();

    std::vector<boost::shared_ptr<Matrix> > Jhess;
    std::vector<boost::shared_ptr<Matrix> > Khess;
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

        //fprintf(outfile,"(%d,%d,%d,%d)\n", P,Q,R,S);

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
               
        double prefactor = 1.0;
        if (P != Q)   prefactor *= 2.0;
        if (R != S)   prefactor *= 2.0;
        if (PQ != RS) prefactor *= 2.0;

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

        //double BxAx=0.0, BxAy=0.0, BxAz=0.0, ByAx=0.0, ByAy=0.0, ByAz=0.0, BzAx=0.0, BzAy=0.0, BzAz=0.0;
        double BxCx=0.0, BxCy=0.0, BxCz=0.0, ByCx=0.0, ByCy=0.0, ByCz=0.0, BzCx=0.0, BzCy=0.0, BzCz=0.0;
        double BxDx=0.0, BxDy=0.0, BxDz=0.0, ByDx=0.0, ByDy=0.0, ByDz=0.0, BzDx=0.0, BzDy=0.0, BzDz=0.0;

        //double CxAx=0.0, CxAy=0.0, CxAz=0.0, CyAx=0.0, CyAy=0.0, CyAz=0.0, CzAx=0.0, CzAy=0.0, CzAz=0.0;
        //double CxBx=0.0, CxBy=0.0, CxBz=0.0, CyBx=0.0, CyBy=0.0, CyBz=0.0, CzBx=0.0, CzBy=0.0, CzBz=0.0;
        double CxDx=0.0, CxDy=0.0, CxDz=0.0, CyDx=0.0, CyDy=0.0, CyDz=0.0, CzDx=0.0, CzDy=0.0, CzDz=0.0;

        //double DxAx=0.0, DxAy=0.0, DxAz=0.0, DyAx=0.0, DyAy=0.0, DyAz=0.0, DzAx=0.0, DzAy=0.0, DzAz=0.0;
        //double DxBx=0.0, DxBy=0.0, DxBz=0.0, DyBx=0.0, DyBy=0.0, DyBz=0.0, DzBx=0.0, DzBy=0.0, DzBz=0.0;
        //double DxCx=0.0, DxCy=0.0, DxCz=0.0, DyCx=0.0, DyCy=0.0, DyCz=0.0, DzCx=0.0, DzCy=0.0, DzCz=0.0;

        delta = 0L;
        for (int p = 0; p < Psize; p++) {
        for (int q = 0; q < Qsize; q++) {
        for (int r = 0; r < Rsize; r++) {
        for (int s = 0; s < Ssize; s++) {
            Dpq = Dtp[p + Poff][q + Qoff];
            Drs = Dtp[r + Roff][s + Soff];
            val = prefactor * Dpq * Drs;
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
        }}}}
        AxBx = -(AxAx + AxCx + AxDx);
        AxBy = -(AxAy + AxCy + AxDy);
        AxBz = -(AxAz + AxCz + AxDz);
        AyBy = -(AyAy + AyCy + AyDy);
        AyBz = -(AyAz + AyCz + AyDz);
        AzBz = -(AzAz + AzCz + AzDz);
        BxCx = -(AxCx + CxCx + CxDx);        
        BxCy = -(AxCy + CxCy + CxDy);        
        BxCz = -(AxCz + CxCz + CxDz);        
        ByCy = -(AyCy + CyCy + CyDy);        
        ByCz = -(AyCz + CyCz + CyDz);        
        BzCz = -(AzCz + CzCz + CzDz);        
        BxDx = -(AxDx + CxDx + DxDx);        
        BxDy = -(AxDy + CxDy + DxDy);        
        BxDz = -(AxDz + CxDz + DxDz);        
        ByDy = -(AyDy + CyDy + DyDy);        
        ByDz = -(AyDz + CyDz + DyDz);        
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
        Jp[Py][Px] += AxAy;
        Jp[Py][Py] += AyAy;
        Jp[Py][Pz] += AyAz;
        Jp[Pz][Px] += AxAz;
        Jp[Pz][Py] += AyAz;
        Jp[Pz][Pz] += AzAz;

        Jp[Qx][Qx] += BxBx;
        Jp[Qx][Qy] += BxBy;
        Jp[Qx][Qz] += BxBz;
        Jp[Qy][Qx] += BxBy;
        Jp[Qy][Qy] += ByBy;
        Jp[Qy][Qz] += ByBz;
        Jp[Qz][Qx] += BxBz;
        Jp[Qz][Qy] += ByBz;
        Jp[Qz][Qz] += BzBz;

        Jp[Rx][Rx] += CxCx;
        Jp[Rx][Ry] += CxCy;
        Jp[Rx][Rz] += CxCz;
        Jp[Ry][Rx] += CxCy;
        Jp[Ry][Ry] += CyCy;
        Jp[Ry][Rz] += CyCz;
        Jp[Rz][Rx] += CxCz;
        Jp[Rz][Ry] += CyCz;
        Jp[Rz][Rz] += CzCz;

        Jp[Sx][Sx] += DxDx;
        Jp[Sx][Sy] += DxDy;
        Jp[Sx][Sz] += DxDz;
        Jp[Sy][Sx] += DxDy;
        Jp[Sy][Sy] += DyDy;
        Jp[Sy][Sz] += DyDz;
        Jp[Sz][Sx] += DxDz;
        Jp[Sz][Sy] += DyDz;
        Jp[Sz][Sz] += DzDz;

        Jp[Px][Qx] += AxBx;
        Jp[Px][Qy] += AxBy;
        Jp[Px][Qz] += AxBz;
        Jp[Py][Qx] += AyBx;
        Jp[Py][Qy] += AyBy;
        Jp[Py][Qz] += AyBz;
        Jp[Pz][Qx] += AzBx;
        Jp[Pz][Qy] += AzBy;
        Jp[Pz][Qz] += AzBz;

        Jp[Qx][Px] += AxBx;
        Jp[Qy][Px] += AxBy;
        Jp[Qz][Px] += AxBz;
        Jp[Qx][Py] += AyBx;
        Jp[Qy][Py] += AyBy;
        Jp[Qz][Py] += AyBz;
        Jp[Qx][Pz] += AzBx;
        Jp[Qy][Pz] += AzBy;
        Jp[Qz][Pz] += AzBz;

        Jp[Px][Rx] += AxCx;
        Jp[Px][Ry] += AxCy;
        Jp[Px][Rz] += AxCz;
        Jp[Py][Rx] += AyCx;
        Jp[Py][Ry] += AyCy;
        Jp[Py][Rz] += AyCz;
        Jp[Pz][Rx] += AzCx;
        Jp[Pz][Ry] += AzCy;
        Jp[Pz][Rz] += AzCz;

        Jp[Rx][Px] += AxCx;
        Jp[Ry][Px] += AxCy;
        Jp[Rz][Px] += AxCz;
        Jp[Rx][Py] += AyCx;
        Jp[Ry][Py] += AyCy;
        Jp[Rz][Py] += AyCz;
        Jp[Rx][Pz] += AzCx;
        Jp[Ry][Pz] += AzCy;
        Jp[Rz][Pz] += AzCz;

        Jp[Px][Sx] += AxDx;
        Jp[Px][Sy] += AxDy;
        Jp[Px][Sz] += AxDz;
        Jp[Py][Sx] += AyDx;
        Jp[Py][Sy] += AyDy;
        Jp[Py][Sz] += AyDz;
        Jp[Pz][Sx] += AzDx;
        Jp[Pz][Sy] += AzDy;
        Jp[Pz][Sz] += AzDz;

        Jp[Sx][Px] += AxDx;
        Jp[Sy][Px] += AxDy;
        Jp[Sz][Px] += AxDz;
        Jp[Sx][Py] += AyDx;
        Jp[Sy][Py] += AyDy;
        Jp[Sz][Py] += AyDz;
        Jp[Sx][Pz] += AzDx;
        Jp[Sy][Pz] += AzDy;
        Jp[Sz][Pz] += AzDz;

        Jp[Qx][Rx] += BxCx;
        Jp[Qx][Ry] += BxCy;
        Jp[Qx][Rz] += BxCz;
        Jp[Qy][Rx] += ByCx;
        Jp[Qy][Ry] += ByCy;
        Jp[Qy][Rz] += ByCz;
        Jp[Qz][Rx] += BzCx;
        Jp[Qz][Ry] += BzCy;
        Jp[Qz][Rz] += BzCz;

        Jp[Rx][Qx] += BxCx;
        Jp[Ry][Qx] += BxCy;
        Jp[Rz][Qx] += BxCz;
        Jp[Rx][Qy] += ByCx;
        Jp[Ry][Qy] += ByCy;
        Jp[Rz][Qy] += ByCz;
        Jp[Rx][Qz] += BzCx;
        Jp[Ry][Qz] += BzCy;
        Jp[Rz][Qz] += BzCz;

        Jp[Qx][Sx] += BxDx;
        Jp[Qx][Sy] += BxDy;
        Jp[Qx][Sz] += BxDz;
        Jp[Qy][Sx] += ByDx;
        Jp[Qy][Sy] += ByDy;
        Jp[Qy][Sz] += ByDz;
        Jp[Qz][Sx] += BzDx;
        Jp[Qz][Sy] += BzDy;
        Jp[Qz][Sz] += BzDz;

        Jp[Sx][Qx] += BxDx;
        Jp[Sy][Qx] += BxDy;
        Jp[Sz][Qx] += BxDz;
        Jp[Sx][Qy] += ByDx;
        Jp[Sy][Qy] += ByDy;
        Jp[Sz][Qy] += ByDz;
        Jp[Sx][Qz] += BzDx;
        Jp[Sy][Qz] += BzDy;
        Jp[Sz][Qz] += BzDz;

        Jp[Rx][Sx] += CxDx;
        Jp[Rx][Sy] += CxDy;
        Jp[Rx][Sz] += CxDz;
        Jp[Ry][Sx] += CyDx;
        Jp[Ry][Sy] += CyDy;
        Jp[Ry][Sz] += CyDz;
        Jp[Rz][Sx] += CzDx;
        Jp[Rz][Sy] += CzDy;
        Jp[Rz][Sz] += CzDz;

        Jp[Sx][Rx] += CxDx;
        Jp[Sy][Rx] += CxDy;
        Jp[Sz][Rx] += CxDz;
        Jp[Sx][Ry] += CyDx;
        Jp[Sy][Ry] += CyDy;
        Jp[Sz][Ry] += CyDz;
        Jp[Sx][Rz] += CzDx;
        Jp[Sy][Rz] += CzDy;
        Jp[Sz][Rz] += CzDz;

        // => Exchange Term <= //

        AxAx=0.0; AxAy=0.0; AxAz=0.0; AyAy=0.0; AyAz=0.0; AzAz=0.0;
        BxBx=0.0; BxBy=0.0; BxBz=0.0; ByBy=0.0; ByBz=0.0; BzBz=0.0;
        CxCx=0.0; CxCy=0.0; CxCz=0.0; CyCy=0.0; CyCz=0.0; CzCz=0.0;
        DxDx=0.0; DxDy=0.0; DxDz=0.0; DyDy=0.0; DyDz=0.0; DzDz=0.0;

        AxBx=0.0; AxBy=0.0; AxBz=0.0; AyBx=0.0; AyBy=0.0; AyBz=0.0; AzBx=0.0; AzBy=0.0; AzBz=0.0;
        AxCx=0.0; AxCy=0.0; AxCz=0.0; AyCx=0.0; AyCy=0.0; AyCz=0.0; AzCx=0.0; AzCy=0.0; AzCz=0.0;
        AxDx=0.0; AxDy=0.0; AxDz=0.0; AyDx=0.0; AyDy=0.0; AyDz=0.0; AzDx=0.0; AzDy=0.0; AzDz=0.0;

        //BxAx=0.0; BxAy=0.0; BxAz=0.0; ByAx=0.0; ByAy=0.0; ByAz=0.0; BzAx=0.0; BzAy=0.0; BzAz=0.0;
        BxCx=0.0; BxCy=0.0; BxCz=0.0; ByCx=0.0; ByCy=0.0; ByCz=0.0; BzCx=0.0; BzCy=0.0; BzCz=0.0;
        BxDx=0.0; BxDy=0.0; BxDz=0.0; ByDx=0.0; ByDy=0.0; ByDz=0.0; BzDx=0.0; BzDy=0.0; BzDz=0.0;

        //CxAx=0.0; CxAy=0.0; CxAz=0.0; CyAx=0.0; CyAy=0.0; CyAz=0.0; CzAx=0.0; CzAy=0.0; CzAz=0.0;
        //CxBx=0.0; CxBy=0.0; CxBz=0.0; CyBx=0.0; CyBy=0.0; CyBz=0.0; CzBx=0.0; CzBy=0.0; CzBz=0.0;
        CxDx=0.0; CxDy=0.0; CxDz=0.0; CyDx=0.0; CyDy=0.0; CyDz=0.0; CzDx=0.0; CzDy=0.0; CzDz=0.0;

        //DxAx=0.0; DxAy=0.0; DxAz=0.0; DyAx=0.0; DyAy=0.0; DyAz=0.0; DzAx=0.0; DzAy=0.0; DzAz=0.0;
        //DxBx=0.0; DxBy=0.0; DxBz=0.0; DyBx=0.0; DyBy=0.0; DyBz=0.0; DzBx=0.0; DzBy=0.0; DzBz=0.0;
        //DxCx=0.0; DxCy=0.0; DxCz=0.0; DyCx=0.0; DyCy=0.0; DyCz=0.0; DzCx=0.0; DzCy=0.0; DzCz=0.0;

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
        }}}}
        AxBx = -(AxAx + AxCx + AxDx);
        AxBy = -(AxAy + AxCy + AxDy);
        AxBz = -(AxAz + AxCz + AxDz);
        AyBy = -(AyAy + AyCy + AyDy);
        AyBz = -(AyAz + AyCz + AyDz);
        AzBz = -(AzAz + AzCz + AzDz);
        BxCx = -(AxCx + CxCx + CxDx);        
        BxCy = -(AxCy + CxCy + CxDy);        
        BxCz = -(AxCz + CxCz + CxDz);        
        ByCy = -(AyCy + CyCy + CyDy);        
        ByCz = -(AyCz + CyCz + CyDz);        
        BzCz = -(AzCz + CzCz + CzDz);        
        BxDx = -(AxDx + CxDx + DxDx);        
        BxDy = -(AxDy + CxDy + DxDy);        
        BxDz = -(AxDz + CxDz + DxDz);        
        ByDy = -(AyDy + CyDy + DyDy);        
        ByDz = -(AyDz + CyDz + DyDz);        
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
        Kp[Py][Px] += AxAy;
        Kp[Py][Py] += AyAy;
        Kp[Py][Pz] += AyAz;
        Kp[Pz][Px] += AxAz;
        Kp[Pz][Py] += AyAz;
        Kp[Pz][Pz] += AzAz;

        Kp[Qx][Qx] += BxBx;
        Kp[Qx][Qy] += BxBy;
        Kp[Qx][Qz] += BxBz;
        Kp[Qy][Qx] += BxBy;
        Kp[Qy][Qy] += ByBy;
        Kp[Qy][Qz] += ByBz;
        Kp[Qz][Qx] += BxBz;
        Kp[Qz][Qy] += ByBz;
        Kp[Qz][Qz] += BzBz;

        Kp[Rx][Rx] += CxCx;
        Kp[Rx][Ry] += CxCy;
        Kp[Rx][Rz] += CxCz;
        Kp[Ry][Rx] += CxCy;
        Kp[Ry][Ry] += CyCy;
        Kp[Ry][Rz] += CyCz;
        Kp[Rz][Rx] += CxCz;
        Kp[Rz][Ry] += CyCz;
        Kp[Rz][Rz] += CzCz;

        Kp[Sx][Sx] += DxDx;
        Kp[Sx][Sy] += DxDy;
        Kp[Sx][Sz] += DxDz;
        Kp[Sy][Sx] += DxDy;
        Kp[Sy][Sy] += DyDy;
        Kp[Sy][Sz] += DyDz;
        Kp[Sz][Sx] += DxDz;
        Kp[Sz][Sy] += DyDz;
        Kp[Sz][Sz] += DzDz;

        Kp[Px][Qx] += AxBx;
        Kp[Px][Qy] += AxBy;
        Kp[Px][Qz] += AxBz;
        Kp[Py][Qx] += AyBx;
        Kp[Py][Qy] += AyBy;
        Kp[Py][Qz] += AyBz;
        Kp[Pz][Qx] += AzBx;
        Kp[Pz][Qy] += AzBy;
        Kp[Pz][Qz] += AzBz;

        Kp[Qx][Px] += AxBx;
        Kp[Qy][Px] += AxBy;
        Kp[Qz][Px] += AxBz;
        Kp[Qx][Py] += AyBx;
        Kp[Qy][Py] += AyBy;
        Kp[Qz][Py] += AyBz;
        Kp[Qx][Pz] += AzBx;
        Kp[Qy][Pz] += AzBy;
        Kp[Qz][Pz] += AzBz;

        Kp[Px][Rx] += AxCx;
        Kp[Px][Ry] += AxCy;
        Kp[Px][Rz] += AxCz;
        Kp[Py][Rx] += AyCx;
        Kp[Py][Ry] += AyCy;
        Kp[Py][Rz] += AyCz;
        Kp[Pz][Rx] += AzCx;
        Kp[Pz][Ry] += AzCy;
        Kp[Pz][Rz] += AzCz;

        Kp[Rx][Px] += AxCx;
        Kp[Ry][Px] += AxCy;
        Kp[Rz][Px] += AxCz;
        Kp[Rx][Py] += AyCx;
        Kp[Ry][Py] += AyCy;
        Kp[Rz][Py] += AyCz;
        Kp[Rx][Pz] += AzCx;
        Kp[Ry][Pz] += AzCy;
        Kp[Rz][Pz] += AzCz;

        Kp[Px][Sx] += AxDx;
        Kp[Px][Sy] += AxDy;
        Kp[Px][Sz] += AxDz;
        Kp[Py][Sx] += AyDx;
        Kp[Py][Sy] += AyDy;
        Kp[Py][Sz] += AyDz;
        Kp[Pz][Sx] += AzDx;
        Kp[Pz][Sy] += AzDy;
        Kp[Pz][Sz] += AzDz;

        Kp[Sx][Px] += AxDx;
        Kp[Sy][Px] += AxDy;
        Kp[Sz][Px] += AxDz;
        Kp[Sx][Py] += AyDx;
        Kp[Sy][Py] += AyDy;
        Kp[Sz][Py] += AyDz;
        Kp[Sx][Pz] += AzDx;
        Kp[Sy][Pz] += AzDy;
        Kp[Sz][Pz] += AzDz;

        Kp[Qx][Rx] += BxCx;
        Kp[Qx][Ry] += BxCy;
        Kp[Qx][Rz] += BxCz;
        Kp[Qy][Rx] += ByCx;
        Kp[Qy][Ry] += ByCy;
        Kp[Qy][Rz] += ByCz;
        Kp[Qz][Rx] += BzCx;
        Kp[Qz][Ry] += BzCy;
        Kp[Qz][Rz] += BzCz;

        Kp[Rx][Qx] += BxCx;
        Kp[Ry][Qx] += BxCy;
        Kp[Rz][Qx] += BxCz;
        Kp[Rx][Qy] += ByCx;
        Kp[Ry][Qy] += ByCy;
        Kp[Rz][Qy] += ByCz;
        Kp[Rx][Qz] += BzCx;
        Kp[Ry][Qz] += BzCy;
        Kp[Rz][Qz] += BzCz;

        Kp[Qx][Sx] += BxDx;
        Kp[Qx][Sy] += BxDy;
        Kp[Qx][Sz] += BxDz;
        Kp[Qy][Sx] += ByDx;
        Kp[Qy][Sy] += ByDy;
        Kp[Qy][Sz] += ByDz;
        Kp[Qz][Sx] += BzDx;
        Kp[Qz][Sy] += BzDy;
        Kp[Qz][Sz] += BzDz;

        Kp[Sx][Qx] += BxDx;
        Kp[Sy][Qx] += BxDy;
        Kp[Sz][Qx] += BxDz;
        Kp[Sx][Qy] += ByDx;
        Kp[Sy][Qy] += ByDy;
        Kp[Sz][Qy] += ByDz;
        Kp[Sx][Qz] += BzDx;
        Kp[Sy][Qz] += BzDy;
        Kp[Sz][Qz] += BzDz;

        Kp[Rx][Sx] += CxDx;
        Kp[Rx][Sy] += CxDy;
        Kp[Rx][Sz] += CxDz;
        Kp[Ry][Sx] += CyDx;
        Kp[Ry][Sy] += CyDy;
        Kp[Ry][Sz] += CyDz;
        Kp[Rz][Sx] += CzDx;
        Kp[Rz][Sy] += CzDy;
        Kp[Rz][Sz] += CzDz;

        Kp[Sx][Rx] += CxDx;
        Kp[Sy][Rx] += CxDy;
        Kp[Sz][Rx] += CxDz;
        Kp[Sx][Ry] += CyDx;
        Kp[Sy][Ry] += CyDy;
        Kp[Sz][Ry] += CyDz;
        Kp[Sx][Rz] += CzDx;
        Kp[Sy][Rz] += CzDy;
        Kp[Sz][Rz] += CzDz;
    }

    for (int thread = 1; thread < nthreads; thread++) {
        Jhess[0]->add(Jhess[thread]);
        Khess[0]->add(Khess[thread]);
    }

    Jhess[0]->scale(0.5);
    Khess[0]->scale(0.5);

    std::map<std::string, boost::shared_ptr<Matrix> > val;
    val["J"] = Jhess[0];
    val["K"] = Khess[0];
    return val;
}

}} // Namespaces
