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

JKGrad::JKGrad(boost::shared_ptr<BasisSet> primary) :
    primary_(primary)
{
    common_init();
}
JKGrad::~JKGrad()
{
}
boost::shared_ptr<JKGrad> JKGrad::build_JKGrad()
{
    Options& options = Process::environment.options;
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    boost::shared_ptr<BasisSet> primary = BasisSet::construct(parser, Process::environment.molecule(), "BASIS");

    if (options.get_str("SCF_TYPE") == "DF") {

        if (options.get_str("DF_BASIS_SCF") == "") {
            primary->molecule()->set_basis_all_atoms(options.get_str("BASIS") + "-JKFIT", "DF_BASIS_SCF");
        }
        boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, primary->molecule(), "DF_BASIS_SCF");

        DFJKGrad* jk = new DFJKGrad(primary,auxiliary);

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
DFJKGrad::DFJKGrad(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<BasisSet> auxiliary) :
    JKGrad(primary), auxiliary_(auxiliary)
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
        #pragma omp parallel for schedule(dynamic) num_threads(df_ints_num_threads_)
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
            #pragma omp parallel for num_threads(omp_num_threads_)
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
            #pragma omp parallel for num_threads(omp_num_threads_)
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
        #pragma omp parallel for schedule(dynamic) num_threads(df_ints_num_threads_)
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
            #pragma omp parallel for num_threads(omp_num_threads_)
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
            #pragma omp parallel for num_threads(omp_num_threads_)
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
        psio_address next_Bij = PSIO_ZERO;

        for (int P = 0; P < naux; P += max_rows) {
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
        psio_address next_Bij = PSIO_ZERO;

        for (int P = 0; P < naux; P += max_rows) {
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
        psio_address next_Bij = PSIO_ZERO;

        for (int P = 0; P < naux; P += max_rows) {
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
        psio_address next_Bij = PSIO_ZERO;

        for (int P = 0; P < naux; P += max_rows) {
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

    #pragma omp parallel for schedule(dynamic) num_threads(df_ints_num_threads_)
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
                    double Wval = 0.5 * perm * Wp[p + oP][q + oQ];
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
            #pragma omp parallel for num_threads(omp_num_threads_)
            for (int P = 0; P < np; P++) {
                C_DGEMM('N','N',nso,na,na,1.0,Cap[0],na,&Aijp[0][P * (ULI) na * na],na,0.0,Amip[P],na);
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N','T',np * (ULI) nso, nso, na, 1.0, Amip[0], na, Cap[0], na, 0.0, Jmnp[0], nso);
        }

        // > Beta < //
        if (!restricted && (do_K_ || do_wK_)) {

            // > Stripe < //
            psio_->read(unit_b_, "(A|ij)", (char*) Aijp[0], sizeof(double) * np * nb * nb, next_Aijb, &next_Aijb);

            // > (A|ij) C_mi -> (A|mj) < //
            #pragma omp parallel for num_threads(omp_num_threads_)
            for (int P = 0; P < np; P++) {
                C_DGEMM('N','N',nso,nb,nb,1.0,Cbp[0],nb,&Aijp[0][P* (ULI) nb * nb],nb,0.0,Amip[P],na);
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N','T',np * (ULI) nso, nso, nb, 1.0, Amip[0], na, Cbp[0], nb, 1.0, Jmnp[0], nso);
        } else if (do_K_ || do_wK_) {
            Jmn->scale(2.0);
        }

        // => K_mn^A <= // 

        // > Alpha < //
        if (do_wK_) {

            // > Stripe < //
            psio_->read(unit_a_, "(A|w|ij)", (char*) Aijp[0], sizeof(double) * np * na * na, next_Awija, &next_Awija);

            // > (A|ij) C_mi -> (A|mj) < //
            #pragma omp parallel for num_threads(omp_num_threads_)
            for (int P = 0; P < np; P++) {
                C_DGEMM('N','N',nso,na,na,1.0,Cap[0],na,&Aijp[0][P * (ULI) na * na],na,0.0,Amip[P],na);
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N','T',np * (ULI) nso, nso, na, 1.0, Amip[0], na, Cap[0], na, 0.0, Kmnp[0], nso);
        }

        // > Beta < //
        if (!restricted && do_wK_) {

            // > Stripe < //
            psio_->read(unit_b_, "(A|w|ij)", (char*) Aijp[0], sizeof(double) * np * nb * nb, next_Awijb, &next_Awijb);

            // > (A|ij) C_mi -> (A|mj) < //
            #pragma omp parallel for num_threads(omp_num_threads_)
            for (int P = 0; P < np; P++) {
                C_DGEMM('N','N',nso,nb,nb,1.0,Cbp[0],nb,&Aijp[0][P* (ULI) nb * nb],nb,0.0,Amip[P],na);
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N','T',np * (ULI) nso, nso, nb, 1.0, Amip[0], na, Cbp[0], nb, 1.0, Kmnp[0], nso);
        } else if (do_wK_) {
            Kmn->scale(2.0);
        }

        // > Integrals < //
        #pragma omp parallel for schedule(dynamic) num_threads(df_ints_num_threads_)
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
            #pragma omp parallel for num_threads(omp_num_threads_)
            for (int P = 0; P < np; P++) {
                C_DGEMM('N','N',nso,na,na,1.0,Cap[0],na,&Aijp[0][P * (ULI) na * na],na,0.0,Amip[P],na);
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N','T',np * (ULI) nso, nso, na, 1.0, Amip[0], na, Cap[0], na, 0.0, Jmnp[0], nso);
        }

        // > Beta < //
        if (!restricted) {

            // > Stripe < //
            psio_->read(unit_b_, "(A|ij)", (char*) Aijp[0], sizeof(double) * np * nb * nb, next_Aijb, &next_Aijb);

            // > (A|ij) C_mi -> (A|mj) < //
            #pragma omp parallel for num_threads(omp_num_threads_)
            for (int P = 0; P < np; P++) {
                C_DGEMM('N','N',nso,nb,nb,1.0,Cbp[0],nb,&Aijp[0][P* (ULI) nb * nb],nb,0.0,Amip[P],na);
            }

            // > (A|mj) C_nj -> (A|mn) < //
            C_DGEMM('N','T',np * (ULI) nso, nso, nb, 1.0, Amip[0], na, Cbp[0], nb, 1.0, Jmnp[0], nso);
        } else { 
            Jmn->scale(2.0);
        }

        // > Integrals < //
        #pragma omp parallel for schedule(dynamic) num_threads(df_ints_num_threads_)
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

    for (int t = 0; t < df_ints_num_threads_; t++) {
        gradients_["Exchange,LR"]->add(wKtemps[t]);
    }
}

}} // Namespaces
