#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>
#include <psiconfig.h>

#include "hf.h"
#include "rhf.h"
#include "uhf.h"
#include "rohf.h"
#include "df.h"
#include <lib3index/3index.h>

//MKL Header
#ifdef HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif

#include <libmints/mints.h>

using namespace std;
using namespace psi;

namespace psi { namespace scf {

DFHF::DFHF(shared_ptr<BasisSet> basis, shared_ptr<PSIO> psio, Options& opt) :
    primary_(basis), psio_(psio), options_(opt)
{
    common_init();
}
DFHF::~DFHF()
{
}
void DFHF::common_init()
{
    is_initialized_ = false;
    is_jk_ = false;
    restricted_ = false;
    is_disk_ = false;

    memory_ = Process::environment.get_memory() / 8L;
    memory_ = (unsigned long int) (0.7 * memory_);

    // Build auxiliary basis from options
    zero_ = BasisSet::zero_ao_basis_set();
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    auxiliary_ = BasisSet::construct(parser, primary_->molecule(), "RI_BASIS_SCF");
    parser.reset();
   
    schwarz_ = shared_ptr<SchwarzSieve>(new SchwarzSieve(primary_, options_.get_double("SCHWARZ_CUTOFF")));
    Jinv_ = shared_ptr<FittingMetric>(new FittingMetric(auxiliary_));
}
void DFHF::initialize()
{
    if (is_initialized_) return;
    is_initialized_ = true;

    // Make a memory decision here 
    int ntri = schwarz_->get_nfun_pairs();
    ULI three_memory = (ULI)primary_->nbf()*ntri;
    ULI two_memory = (ULI)auxiliary_->nbf()*auxiliary_->nbf(); 

    is_disk_ = (three_memory + two_memory < memory_ ? false : true);

    if (is_jk_) {
        if (is_disk_) 
            initialize_JK_disk();
        else
            initialize_JK_core();    
    } else {
        if (is_disk_) 
            initialize_J_disk();
        else
            initialize_J_core();    
    }
}
void DFHF::initialize_JK_core()
{
    int ntri = schwarz_->get_nfun_pairs();
    ULI three_memory = (ULI)primary_->nbf()*ntri;
    ULI two_memory = (ULI)auxiliary_->nbf()*auxiliary_->nbf(); 

    int nthread = 1;
    #ifdef _OPENMP
        if (options_.get_int("RI_INTS_NUM_THREADS") == 0) {
            nthread = omp_get_max_threads();
        } else {
            nthread = options_.get_int("RI_INTS_NUM_THREADS");
        }
    #endif
    int rank = 0;
    
    Qmn_ = shared_ptr<Matrix>(new Matrix("Qmn (Fitted Integrals)", 
        auxiliary_->nbf(), ntri)); 
    double** Qmnp = Qmn_->pointer();

    //Get a TEI for each thread
    shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero_, primary_, primary_));
    const double **buffer = new const double*[nthread];
    shared_ptr<TwoBodyAOInt> *eri = new shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri[Q] = shared_ptr<TwoBodyAOInt>(rifactory->eri());
        buffer[Q] = eri[Q]->buffer();
    }

    long int* schwarz_shell_pairs = schwarz_->get_schwarz_shells_reverse();
    long int* schwarz_fun_pairs = schwarz_->get_schwarz_funs_reverse();
    int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
    int index;
    //The integrals (A|mn)
    timer_on("(A|mn)");
    #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
    for (MU=0; MU < primary_->nshell(); ++MU) {
        #ifdef _OPENMP
            rank = omp_get_thread_num();
            //fprintf(outfile,"  Thread %d doing MU = %d",rank,MU); fflush(outfile);
        #endif
        nummu = primary_->shell(MU)->nfunction();
        for (NU=0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU)->nfunction();
            if (schwarz_shell_pairs[MU*(MU+1)/2+NU] > -1) {
                for (Pshell=0; Pshell < auxiliary_->nshell(); ++Pshell) {
                    numP = auxiliary_->shell(Pshell)->nfunction();
                    eri[rank]->compute_shell(Pshell, 0, MU, NU);
                    for (mu=0 ; mu < nummu; ++mu) {
                        omu = primary_->shell(MU)->function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = primary_->shell(NU)->function_index() + nu;
                            if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] > -1) {
                                for (P=0; P < numP; ++P) {
                                    PHI = auxiliary_->shell(Pshell)->function_index() + P;
                                    Qmnp[PHI][schwarz_fun_pairs[omu*(omu+1)/2+onu]] = buffer[rank][P*nummu*numnu + mu*numnu + nu];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    timer_off("(A|mn)");

    delete []buffer;
    delete []eri;
   
    if (options_.get_str("FITTING_TYPE") == "EIG") {
        Jinv_->form_eig_inverse();
    } else {
        throw PSIEXCEPTION("Fitting Metric type is not implemented.");
    } 
    
    double** Jinvp = Jinv_->get_metric()->pointer();
    ULI max_cols = (memory_-three_memory-two_memory) / auxiliary_->nbf();
    if (max_cols < 1)
        max_cols = 1;
    if (max_cols > ntri)
        max_cols = ntri;
    shared_ptr<Matrix> temp(new Matrix("Qmn buffer", auxiliary_->nbf(), max_cols));
    double** tempp = temp->pointer();   

    int nblocks = ntri / max_cols;
    if ((ULI)nblocks*max_cols != ntri) nblocks++;

    int ncol = 0;
    int col = 0;
    for (int block = 0; block < nblocks; block++) {
        ncol = max_cols;
        if (col + ncol > ntri)
            ncol = ntri - col; 
        
        C_DGEMM('N','N',auxiliary_->nbf(), ncol, auxiliary_->nbf(), 1.0,
            Jinvp[0], auxiliary_->nbf(), &Qmnp[0][col], ntri, 0.0, 
            tempp[0], max_cols);

        for (int Q = 0; Q < auxiliary_->nbf(); Q++) {
            C_DCOPY(ncol, tempp[Q], 1, &Qmnp[Q][col], 1);
        }

        col += ncol;
    }    
    
    Jinv_.reset();
}
void DFHF::initialize_JK_disk()
{
}
void DFHF::initialize_J_core()
{
    int ntri = schwarz_->get_nfun_pairs();
    ULI three_memory = (ULI)primary_->nbf()*ntri;
    ULI two_memory = (ULI)auxiliary_->nbf()*auxiliary_->nbf(); 

    int nthread = 1;
    #ifdef _OPENMP
        if (options_.get_int("RI_INTS_NUM_THREADS") == 0) {
            nthread = omp_get_max_threads();
        } else {
            nthread = options_.get_int("RI_INTS_NUM_THREADS");
        }
    #endif
    int rank = 0;
    
    Qmn_ = shared_ptr<Matrix>(new Matrix("Qmn (Fitted Integrals)", 
        auxiliary_->nbf(), ntri)); 
    double** Qmnp = Qmn_->pointer();

    //Get a TEI for each thread
    shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero_, primary_, primary_));
    const double **buffer = new const double*[nthread];
    shared_ptr<TwoBodyAOInt> *eri = new shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri[Q] = shared_ptr<TwoBodyAOInt>(rifactory->eri());
        buffer[Q] = eri[Q]->buffer();
    }

    long int* schwarz_shell_pairs = schwarz_->get_schwarz_shells_reverse();
    long int* schwarz_fun_pairs = schwarz_->get_schwarz_funs_reverse();
    int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
    int index;
    //The integrals (A|mn)
    timer_on("(A|mn)");
    #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
    for (MU=0; MU < primary_->nshell(); ++MU) {
        #ifdef _OPENMP
            rank = omp_get_thread_num();
            //fprintf(outfile,"  Thread %d doing MU = %d",rank,MU); fflush(outfile);
        #endif
        nummu = primary_->shell(MU)->nfunction();
        for (NU=0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU)->nfunction();
            if (schwarz_shell_pairs[MU*(MU+1)/2+NU] > -1) {
                for (Pshell=0; Pshell < auxiliary_->nshell(); ++Pshell) {
                    numP = auxiliary_->shell(Pshell)->nfunction();
                    eri[rank]->compute_shell(Pshell, 0, MU, NU);
                    for (mu=0 ; mu < nummu; ++mu) {
                        omu = primary_->shell(MU)->function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = primary_->shell(NU)->function_index() + nu;
                            if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] > -1) {
                                for (P=0; P < numP; ++P) {
                                    PHI = auxiliary_->shell(Pshell)->function_index() + P;
                                    Qmnp[PHI][schwarz_fun_pairs[omu*(omu+1)/2+onu]] = buffer[rank][P*nummu*numnu + mu*numnu + nu];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    timer_off("(A|mn)");

    delete []buffer;
    delete []eri;
  
    // TODO respect pivoting 
    Jinv_->form_cholesky_factor();
}
void DFHF::initialize_J_disk()
{
}
void DFHF::compute_JK_block(shared_ptr<Matrix> Qmn, int nrows)
{
    int nbf = primary_->nbf();
    int naux = nrows;
    int nalpha = nalpha_[0];
    int nbeta = 0;
    if (!restricted_)
        nbeta = nbeta_[0]; 
    int ntri = schwarz_->get_nfun_pairs();

    double** Qmnp = Qmn->pointer();

    shared_ptr<Matrix> Dt(new Matrix("D Total",nbf,nbf)); 
    double** Dtp = Dt->pointer(); 
    if (restricted_) {
        Dt->copy(Da_);
    } else {
        Dt->copy(Da_);
        Dt->add(Db_);
    }  
    double* Dtri = new double[ntri];
    double* Jtri = new double[ntri];
    int* schwarz_funs = schwarz_->get_schwarz_funs();
    for (int munu = 0; munu < ntri; munu++) {
        int mu = schwarz_funs[2*munu];
        int nu = schwarz_funs[2*munu + 1];
        double perm = (mu == nu ? 1.0 : 2.0);
        Dtri[munu] = perm * Dtp[mu][nu];
    }     
    Dt.reset();
        
    double *dQ = new double[naux];
    C_DGEMV('N', naux, ntri, 1.0, Qmnp[0], ntri, Dtri, 1, 0.0, dQ, 1);
    C_DGEMV('T', naux, ntri, 1.0, Qmnp[0], ntri, dQ, 1, 0.0, Jtri, 1);  

    double** Jp = Ja_->pointer();
    for (int munu = 0; munu < ntri; munu++) {
        int mu = schwarz_funs[2*munu];
        int nu = schwarz_funs[2*munu + 1];
        Jp[mu][nu] += Jtri[munu];
        if (mu != nu)
            Jp[nu][mu] += Jtri[munu];
    }     
    
    delete[] dQ;
    delete[] Jtri;
    delete[] Dtri;

    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif
    int rank = 0;

    shared_ptr<Matrix> E(new Matrix("E_im^Q", nbf, nalpha*naux));
    double** Ep = E->pointer();
    double** Cap = Ca_->pointer();
    double** Kap = Ka_->pointer(); 
    double** Cbp;
    double** Kbp;
    if (!restricted_) {
        Cbp = Cb_->pointer();
        Kbp = Kb_->pointer(); 
    }

    //QS temp matrix for DGEMM
    double*** QS = new double**[nthread];
    for (int T = 0; T < nthread; T++)
        QS[T] = block_matrix(naux,nbf);
    // Temp matrix for sparse DGEMM if sieve exists
    double*** Ctemp = new double**[nthread];
    for (int T = 0; T < nthread; T++)
        Ctemp[T] = block_matrix(nalpha,nbf);
    // Index array for non-canonical ordering of mn
    int** m_ij_indices = init_int_matrix(nbf,nbf);
    // Index array of n for given m (in order of above)
    int** n_indices = init_int_matrix(nbf,nbf);
    // sizes of above for schwarz sieve
    int* index_sizes = init_int_array(nbf);

    for (int ij = 0; ij<ntri; ij++) {
        int m = schwarz_funs[2*ij];
        int n = schwarz_funs[2*ij+1];
        
        m_ij_indices[m][index_sizes[m]] = ij;
        n_indices[m][index_sizes[m]] = n;
        index_sizes[m]++;
        if (m != n){
            m_ij_indices[n][index_sizes[n]] = ij;
            n_indices[n][index_sizes[n]] = m;
            index_sizes[n]++;
        }
    }

    #ifdef HAVE_MKL
        int mkl_nthread = mkl_get_max_threads();
        mkl_set_num_threads(1);
    #endif

    int m, n , ij, index;
    #pragma omp parallel for private (m, n , ij, index, rank) schedule (dynamic)
    for (m = 0; m<nbf; m++) {

        rank = 0;
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif

        int n, ij;
        for (index = 0; index<index_sizes[m]; index++) {
            ij = m_ij_indices[m][index];
            n = n_indices[m][index];

            C_DCOPY(naux,&Qmnp[0][ij],ntri,&QS[rank][0][index],nbf);
            C_DCOPY(nalpha,Cap[n],1,&Ctemp[rank][0][index],nbf);
        }

        C_DGEMM('N','T',nalpha,naux,index_sizes[m],1.0,Ctemp[rank][0],nbf,QS[rank][0],nbf, 0.0, Ep[m], naux);
    }

    #ifdef HAVE_MKL
        mkl_set_num_threads(mkl_nthread);
    #endif

    C_DGEMM('N','T',nbf,nbf,naux*nalpha,1.0,Ep[0],naux*nalpha,Ep[0],naux*nalpha,1.0,Kap[0], nbf);

    if (!restricted_) {
        #ifdef HAVE_MKL
            int mkl_nthread = mkl_get_max_threads();
            mkl_set_num_threads(1);
        #endif

        #pragma omp parallel for private (m, n , ij, index, rank) schedule (dynamic)
        for (m = 0; m<nbf; m++) {

            rank = 0;
            #ifdef _OPENMP
                rank = omp_get_thread_num();
            #endif

            int n, ij;
            for (index = 0; index<index_sizes[m]; index++) {
                ij = m_ij_indices[m][index];
                n = n_indices[m][index];
                C_DCOPY(naux,&Qmnp[0][ij],ntri,&QS[rank][0][index],nbf);
                C_DCOPY(nbeta,Cbp[n],1,&Ctemp[rank][0][index],nbf);
            }

            C_DGEMM('N','T',nbeta,naux,index_sizes[m],1.0,Ctemp[rank][0],nbf,QS[rank][0],nbf, 0.0, Ep[m], naux);
        }

        #ifdef HAVE_MKL
            mkl_set_num_threads(mkl_nthread);
        #endif

        C_DGEMM('N','T',nbf,nbf,naux*nbeta,1.0,Ep[0],naux*nbeta,Ep[0],naux*nalpha,1.0,Kbp[0], nbf);
    }
    for (int thread = 0; thread < nthread; thread++) {
        free_block(QS[thread]);
        free_block(Ctemp[thread]);
    }
    delete[] QS;
    delete[] Ctemp; 
    free(m_ij_indices[0]);
    free(m_ij_indices);
    free(n_indices[0]);
    free(n_indices);
    free(index_sizes);
}
void DFHF::compute_J_core()
{
    int nbf = primary_->nbf();
    int naux = auxiliary_->nbf();
    int ntri = schwarz_->get_nfun_pairs();

    double** Qmnp = Qmn_->pointer();

    shared_ptr<Matrix> Dt(new Matrix("D Total",nbf,nbf)); 
    double** Dtp = Dt->pointer(); 
    if (restricted_) {
        Dt->copy(Da_);
    } else {
        Dt->copy(Da_);
        Dt->add(Db_);
    }  
    double* Dtri = new double[ntri];
    double* Jtri = new double[ntri];
    int* schwarz_funs = schwarz_->get_schwarz_funs();
    for (int munu = 0; munu < ntri; munu++) {
        int mu = schwarz_funs[2*munu];
        int nu = schwarz_funs[2*munu + 1];
        double perm = (mu == nu ? 1.0 : 2.0);
        Dtri[munu] = perm * Dtp[mu][nu];
    }     
    Dt.reset();
        
    double *dQ = new double[naux];
    C_DGEMV('N', naux, ntri, 1.0, Qmnp[0], ntri, Dtri, 1, 0.0, dQ, 1);
    // TODO pivoting
    C_DPOTRS('L', naux, 1, Jinv_->get_metric()->pointer()[0], naux, dQ, naux);  
    C_DGEMV('T', naux, ntri, 1.0, Qmnp[0], ntri, dQ, 1, 0.0, Jtri, 1);  

    double** Jp = Ja_->pointer();
    for (int munu = 0; munu < ntri; munu++) {
        int mu = schwarz_funs[2*munu];
        int nu = schwarz_funs[2*munu + 1];
        Jp[mu][nu] += Jtri[munu];
        if (mu != nu)
            Jp[nu][mu] += Jtri[munu];
    }     
    
    delete[] dQ;
    delete[] Jtri;
    delete[] Dtri;
}
void DFHF::form_J_DF_RHF()
{
    initialize(); 

    if (is_disk_) {
    } else {
        compute_J_core();
    }
}
void DFHF::form_JK_DF_RHF()
{
    initialize(); 
    
    if (is_disk_) {
    } else {
        // TODO compute blocking
        compute_JK_block(Qmn_, auxiliary_->nbf());             
    }
}
void DFHF::form_J_DF_UHF()
{
    initialize(); 

    if (is_disk_) {
    } else {
        compute_J_core();
    }
}
void DFHF::form_JK_DF_UHF()
{
    initialize(); 

    if (is_disk_) {
    } else {
        // TODO compute blocking
        compute_JK_block(Qmn_, auxiliary_->nbf());             
    }
}

}}
