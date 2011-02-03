/***************************************************************************
*
*       df.cc in psi4/src/lib/libscf_solver
*       By Rob Parrish, CCMST Georgia Tech
*       robparrish@gmail.com
*       14 June 2010
*
*       Canonical (delocalized) density fitting routines for SCF
*
*       This code is heavily based on BLAS calls, and may be easily
*       threaded simply by using threaded BLAS/LAPACK. MKL in particular
*       is extremely efficient, as loops with internal DGEMMS may be
*       parallelized with OpenMP, and the DGEMMS may be performed with one
*       thread each to avoid thrashing
*
*
*****************************************************************************/


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

#include "hf.h"
#include "rhf.h"
#include "uhf.h"
#include "rohf.h"

//MKL Header
#ifdef _MKL
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

void HF::form_B()
{
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                        FORM B
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    //Welcome to form_B, responsible for creation of the three-index tensor with embedded fitting metric
    //on core or disk.

    //Make sure we're in the right spot
    if (print_)
        fprintf(outfile, "\n  Computing Integrals using Density Fitting\n");
    //TODO: Add support for molecular symmetry
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Must run in C1 for now.\n"); fflush(outfile);
        abort();
    }

    //Grab norbs and ndocc and get the ri basis up to the class scope
    int norbs = basisset_->nbf();
    int ndocc = doccpi_[0];
    //Amount of memory available for DF Algorithm, in doubles
    df_memory_ = (long)(memory_/sizeof(double)*(1.0-MEMORY_SAFETY_FACTOR));

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                        RESTART?
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (options_.get_bool("RI_SCF_RESTART"))
    {
        //restart!! Use existing 3-index tensor on disk
        if (print_)
            fprintf(outfile,"\n  Attempting to restart existing DF-SCF computation\n"); fflush(outfile);

        //First read in tensor sizes to set the bookkeeping up
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->read(PSIF_DFSCF_BJ,"N_TRI",(char *) &(ntri_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->read(PSIF_DFSCF_BJ,"N_TRI_NAIVE",(char *) &(ntri_naive_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->read(PSIF_DFSCF_BJ,"N_AUX_RAW",(char *) &(naux_raw_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->read(PSIF_DFSCF_BJ,"N_AUX_FIN",(char *) &(naux_fin_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);

        //Use ri_pair_mu_ and ri_pair_nu_ and ri_back_map_ to keep track of things
        //Across schwarz sieve and unfortunate shell indexing
        ri_pair_nu_ = init_int_array(ntri_naive_);
        ri_pair_mu_ = init_int_array(ntri_naive_);
        ri_back_map_ = init_int_array(norbs*(norbs+1)/2);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->read(PSIF_DFSCF_BJ,"RI_PAIR_BACK",(char *) &(ri_back_map_[0]),sizeof(int)*(norbs*(norbs+1)/2),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->read(PSIF_DFSCF_BJ,"RI_PAIR_MU",(char *) &(ri_pair_mu_[0]),sizeof(int)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->read(PSIF_DFSCF_BJ,"RI_PAIR_NU",(char *) &(ri_pair_nu_[0]),sizeof(int)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);

        //Now determine the storage type. It might change if you switch machines or change memory available
        string storage_type;
        storage_type = options_.get_str("RI_SCF_STORAGE");

        //Detrmine storage algorithm based on user input
        //Size of the three-index tensor
        unsigned long memA = ntri_naive_*(long)naux_fin_;
        //Size of the fitting metric
        unsigned long memJ = naux_fin_*(long)naux_fin_;
        if (storage_type == "CORE")
            df_storage_ = core;
        else if (storage_type == "DISK")
            df_storage_ = disk;
        else if (storage_type == "DEFAULT")
        {
            //set df_storage_ based on available memory
            //memJ here is padding for other uses, we may need more in the end
            if (memA+memJ+memJ<df_memory_)
                df_storage_ = core; //Full in-core
            else
                df_storage_ = disk; //Disk
        }

        if (df_storage_ == core)
            fprintf(outfile,"\n  Density Fitting Algorithm proceeding on Core.\n");
        else if (df_storage_ == disk)
            fprintf(outfile,"\n  Density Fitting Algorithm proceeding on Disk\n");
        fflush(outfile);

        unsigned long int max_rows;
        if (df_storage_ == core)
            max_rows = ((df_memory_-naux_raw_*(long)ntri_naive_)/((long)(norbs*ndocc)));
        else
            max_rows = ((df_memory_)/((long)(norbs*ndocc+ntri_naive_)));

        if (max_rows > naux_fin_)
            max_rows = naux_fin_;
        if (max_rows < 1)
            max_rows = 1;

        if (print_>1) {
            fprintf(outfile,"\n  DF MEMORY USAGE:\n");
            fprintf(outfile,"  Total memory:           %15ld doubles, %15ld bytes.\n",memory_/sizeof(double),memory_);
            fprintf(outfile,"  DF memory:              %15ld doubles, %15ld bytes.\n",df_memory_,df_memory_*sizeof(double));
            fprintf(outfile,"  3-index storage:        %15ld doubles, %15ld bytes.\n",memA,memA*sizeof(double));
            fprintf(outfile,"  Exchange size:          %15ld doubles, %15ld bytes.\n",ndocc*naux_fin_*(long)norbs,ndocc*naux_fin_*(long)norbs*sizeof(double));
            fprintf(outfile,"  J memory:               %15ld doubles, %15ld bytes.\n",memJ,memJ*sizeof(double));
            fprintf(outfile,"  Max rows (iterations):  %15ld rows.\n",max_rows);
            fprintf(outfile,"  Nso:                    %15d functions.\n",norbs);
            fprintf(outfile,"  Npairs:                 %15d pairs.\n",ntri_naive_);
            fprintf(outfile,"  Ndocc:                  %15d orbitals.\n",ndocc);
            fprintf(outfile,"  Naux raw:               %15d functions.\n",naux_raw_);
            fprintf(outfile,"  Naux finished:          %15d functions.\n",naux_fin_);
            fprintf(outfile,"\n");
            fflush(outfile);
        }
        if (df_storage_ == core) {
            //We need the three-index tensor in the core
            //fprintf(outfile,"  n_tri_ %d, n_tri_naive_ %d, naux_fin_ %d\n",ntri_,ntri_naive_,naux_fin_); fflush(outfile);
            B_ia_P_ = block_matrix(naux_fin_,ntri_naive_);
            next_PSIF_DFSCF_BJ = PSIO_ZERO;
            psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(B_ia_P_[0][0]),sizeof(double)*ntri_naive_*naux_fin_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        }
        psio_->close(PSIF_DFSCF_BJ,1); //we'll need to reuse this guy (probably)
        return;
    }

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                    SCHWARZ SIEVE
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //Form the schwarz sieve
    timer_on("Schwarz Sieve");

    int sig_fun_pairs = 0;
    int sig_shell_pairs = 0;

    int *schwarz_shell_pairs;
    int *schwarz_fun_pairs;
    if (schwarz_ > 0.0) {

        schwarz_shell_pairs = init_int_array(basisset_->nshell()*(basisset_->nshell()+1)/2);
        schwarz_fun_pairs = init_int_array(norbs*(norbs+1)/2);
        double* max_shell_val = init_array(basisset_->nshell()*(basisset_->nshell()+1)/2);;
        double* max_fun_val = init_array(norbs*(norbs+1)/2);
        double max_global_val = 0.0;

        IntegralFactory schwarzfactory(basisset_,basisset_,basisset_,basisset_);
        shared_ptr<TwoBodyAOInt> eri = shared_ptr<TwoBodyAOInt>(schwarzfactory.eri());
        const double *buffer = eri->buffer();

        int MU, NU, mu, nu,omu,onu, nummu, numnu, index;
        int MUNU = 0;
        int munu = 0;
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU, ++MUNU) {
                numnu = basisset_->shell(NU)->nfunction();
                eri->compute_shell(MU,NU,MU,NU);
                for (mu=0; mu < nummu; ++mu) {
                    omu = basisset_->shell(MU)->function_index() + mu;
                    for (nu=0; nu < numnu; ++nu) {
                        onu = basisset_->shell(NU)->function_index() + nu;

                        if (omu>=onu) {
                            index = mu*(numnu*nummu*numnu+numnu)+nu*(nummu*numnu+1);
                            //int check = mu*numnu*nummu*numnu+nu*nummu*numnu+mu*numnu+nu;
                            //fprintf(outfile,"   Index = %d, (%d %d| %d %d) = %20.15f\n",index,omu,onu,omu,onu, buffer[index] );
                            if (max_global_val<abs(buffer[index]))
                                max_global_val = abs(buffer[index]);
                            if (max_shell_val[MUNU]<abs(buffer[index]))
                                max_shell_val[MUNU] = abs(buffer[index]);
                            if (max_fun_val[omu*(omu+1)/2+onu]<abs(buffer[index]))
                                max_fun_val[omu*(omu+1)/2+onu] = abs(buffer[index]);
                        }
                    }
                }
            }
        }
        for (int ij = 0; ij < norbs*(norbs+1)/2; ij ++)
            if (max_fun_val[ij]*max_global_val>=schwarz_*schwarz_){
                schwarz_fun_pairs[ij] = 1;
                sig_fun_pairs++;
            }
        for (int ij = 0; ij < basisset_->nshell()*(basisset_->nshell()+1)/2; ij ++)
            if (max_shell_val[ij]*max_global_val>=schwarz_*schwarz_){
                schwarz_shell_pairs[ij] = 1;
                sig_shell_pairs++;
            }

        //for (int i = 0, ij = 0; i<norbs; i++)
            //for (int j = 0; j<=i; j++, ij++)
                //fprintf(outfile,"   Function pair %d = (%d,%d), Max val %14.10f, Max Integral %14.10f, Significant %s\n",ij,i,j,max_fun_val[ij],max_fun_val[ij]*max_global_val,(schwarz_fun_pairs[ij])?"YES":"NO");
        //fprintf(outfile,"\n  Shell Pair Schwarz Sieve, schwarz_ = %14.10f:\n",schwarz_);
        //for (int i = 0, ij = 0; i<basisset_->nshell(); i++)
            //for (int j = 0; j<=i; j++, ij++)
                //fprintf(outfile,"   Shell pair %d = (%d,%d), Max val %14.10f, Max Integral %14.10f, Significant %s\n",ij,i,j,max_shell_val[ij],max_shell_val[ij]*max_global_val,(schwarz_shell_pairs[ij])?"YES":"NO");
        //fprintf(outfile, "\n");

        free(max_fun_val);
        free(max_shell_val);

        ntri_naive_ = sig_fun_pairs; //Matrix size for most of the algorithm
        ntri_ = ntri_naive_; //For now!

    } else {
        ntri_ = norbs*(norbs+1)/2; //Yeah, eat it
        ntri_naive_ = norbs*(norbs+1)/2;
        schwarz_shell_pairs = init_int_array(basisset_->nshell()*(basisset_->nshell()+1)/2);
        schwarz_fun_pairs = init_int_array(norbs*(norbs+1)/2);
        for (int ij = 0; ij < basisset_->nshell()*(basisset_->nshell()+1)/2; ij++)
            schwarz_shell_pairs[ij] = 1;
        for (int ij = 0; ij < ntri_; ij++)
            schwarz_fun_pairs[ij] = 1;
    }

    //Use ri_pair_mu_ and ri_pair_nu_ (and the back map) to keep track of things
    //Across schwarz sieve and unfortunate shell indexing
    ri_pair_nu_ = init_int_array(ntri_naive_);
    ri_pair_mu_ = init_int_array(ntri_naive_);
    ri_back_map_ = init_int_array(norbs*(norbs+1)/2);

    {//<<< Drop out of scope

        //Set up schwarz bookkeeping
        int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
        int index;
        for (MU=0, index = 0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU) {
                numnu = basisset_->shell(NU)->nfunction();
                if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                    for (mu=0 ; mu < nummu; ++mu) {
                        omu = basisset_->shell(MU)->function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = basisset_->shell(NU)->function_index() + nu;
                            if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                ri_pair_mu_[index] = omu;
                                ri_pair_nu_[index] = onu;
                                ri_back_map_[omu*(omu+1)/2+onu] = index;
                                index++;
                            }
                        }
                    }
                }
            }
        }
    } //back into scope

    timer_off("Schwarz Sieve");

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                  FORM FITTING METRIC
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (options_.get_str("RI_FITTING_TYPE") == "RAW")
        form_Wm12_raw();
    else if (options_.get_str("RI_FITTING_TYPE") == "FINISHED")
        form_Wm12_fin();
    else if (options_.get_str("RI_FITTING_TYPE") == "CHOLESKY")
        form_Wp12_chol();

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                    DETERMINE STORAGE
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //Detrmine storage algorithm based on user input
    //This far down as values may change due to schwarz sieve and
    //fitting metric considerations

    //Size of the three-index tensor
    unsigned long memA = ntri_naive_*(long)naux_raw_;
    //Size of the fitting metric
    unsigned long memJ = naux_raw_*(long)naux_fin_;

    string storage_type;
    storage_type = options_.get_str("RI_SCF_STORAGE");

    if (storage_type == "CORE")
        df_storage_ = core;
    else if (storage_type == "DISK")
        df_storage_ = disk;
    else if (storage_type == "DEFAULT")
    {
        //set df_storage_ semi-heuristically based on available memory
        //The second memJ is padding for transformations, we may need more
        if (memA+memJ+memJ<df_memory_)
            df_storage_ = core; //Full in-core
        else
            df_storage_ = disk; //Disk
    }

    if (df_storage_ == core)
        fprintf(outfile,"\n  Density Fitting Algorithm proceeding on Core.\n");
    else if (df_storage_ == disk)
        fprintf(outfile,"\n  Density Fitting Algorithm proceeding on Disk\n");
    fflush(outfile);

    //Zero basis
    shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                        FORM B (FINALLY)
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //Form the AO tensor (A|mn), transform to (B|mn) by embedding fitting metric W
    timer_on("Overall (Q|mn)");

    double three_index_cutoff = options_.get_double("THREE_INDEX_CUTOFF");

    //Threading values (defaults to single thread unless OpenMP exists)
    //The user may specify fewer threads for this due to memory requirements
    int nthread = 1;
    #ifdef _OPENMP
        if (options_.get_int("RI_INTS_NUM_THREADS") == 0) {
            nthread = omp_get_max_threads();
           // printf("OMP says %d threads\n",nthread);
        } else {
            nthread = options_.get_int("RI_INTS_NUM_THREADS");
        }
    #endif
    int rank = 0;

    if (df_storage_ == core)
    {
        //Build (A|mn) on core, and then transform in place using as large of a buffer as possible

        //Overall containers
        IntegralFactory rifactory(basisset_, basisset_, ribasis_, zero);
        B_ia_P_ = block_matrix(naux_raw_,ntri_naive_);

        //Get a TEI for each thread
        const double **buffer = new const double*[nthread];
        shared_ptr<TwoBodyAOInt> *eri = new shared_ptr<TwoBodyAOInt>[nthread];
        for (int Q = 0; Q<nthread; Q++) {
            eri[Q] = shared_ptr<TwoBodyAOInt>(rifactory.eri());
            buffer[Q] = eri[Q]->buffer();
        }

        int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
        int index;
        //The integrals (A|mn)
        timer_on("(A|mn)");
        #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            #ifdef _OPENMP
                rank = omp_get_thread_num();
                //fprintf(outfile,"  Thread %d doing MU = %d",rank,MU); fflush(outfile);
            #endif
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU) {
                numnu = basisset_->shell(NU)->nfunction();
                if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                    for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
                        numP = ribasis_->shell(Pshell)->nfunction();
                        eri[rank]->compute_shell(MU, NU, Pshell, 0);
                        for (mu=0 ; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                            for (nu=0; nu < numnu; ++nu) {
                                onu = basisset_->shell(NU)->function_index() + nu;
                                if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                    for (P=0; P < numP; ++P) {
                                        PHI = ribasis_->shell(Pshell)->function_index() + P;
                                        B_ia_P_[PHI][ri_back_map_[omu*(omu+1)/2+onu]] = buffer[rank][mu*numnu*numP+nu*numP+P];
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
        //print_mat(B_ia_P_, naux_raw_,ntri_ ,outfile);

        //Columns per batch for integral fitting
        unsigned long int max_cols = ((df_memory_-memA-memJ)/((long)naux_raw_));
        if (max_cols > ntri_naive_)
            max_cols = ntri_naive_;
        if (max_cols < 1)
            max_cols = 1; //You need to be able to spare that much at least!

        //Rows per batch for iterations (for diagnostics)
        unsigned long int max_rows = ((df_memory_-naux_raw_*(long)ntri_naive_)/((long)(norbs*ndocc)));
        if (max_rows > naux_fin_)
            max_rows = naux_fin_;
        if (max_rows < 1)
            max_rows = 1;

        if (print_>1) {
            fprintf(outfile,"\n  DF MEMORY USAGE:\n");
            fprintf(outfile,"  Total memory:           %15ld doubles, %15ld bytes.\n",memory_/sizeof(double),memory_);
            fprintf(outfile,"  DF memory:              %15ld doubles, %15ld bytes.\n",df_memory_,df_memory_*sizeof(double));
            fprintf(outfile,"  3-index storage:        %15ld doubles, %15ld bytes.\n",memA,memA*sizeof(double));
            fprintf(outfile,"  Exchange size:          %15ld doubles, %15ld bytes.\n",ndocc*naux_fin_*(long)norbs,ndocc*naux_fin_*(long)norbs*sizeof(double));
            fprintf(outfile,"  J memory:               %15ld doubles, %15ld bytes.\n",memJ,memJ*sizeof(double));
            fprintf(outfile,"  Max cols (integrals):   %15ld columns.\n",max_cols);
            fprintf(outfile,"  Max rows (iterations):  %15ld rows.\n",max_rows);
            fprintf(outfile,"  Nso:                    %15d functions.\n",norbs);
            fprintf(outfile,"  Npairs:                 %15d pairs.\n",ntri_naive_);
            fprintf(outfile,"  Ndocc:                  %15d orbitals.\n",ndocc);
            fprintf(outfile,"  Naux raw:               %15d functions.\n",naux_raw_);
            fprintf(outfile,"  Naux finished:          %15d functions.\n",naux_fin_);
            fprintf(outfile,"\n");
            fflush(outfile);
        }


        //fprintf(outfile,"  Max cols %d\n",max_cols);
        //fprintf(outfile,"  (A|mn):");
        //for (int ij = 0; ij < ntri_naive_; ij++) {
        //    fprintf(outfile,"  Column %d: mu = %d, nu = %d\n",ij+1,ri_pair_mu_[ij],ri_pair_nu_[ij]);
        //}

        //print_mat(B_ia_P_,naux_raw_,ntri_naive_,outfile);

    // Transformation buffer
        double** Temp1;
        if (options_.get_str("RI_FITTING_TYPE") == "RAW" || options_.get_str("RI_FITTING_TYPE") == "FINISHED")
            Temp1 = block_matrix(naux_raw_,max_cols);
        else
            Temp1 = block_matrix(max_cols,naux_raw_);

        timer_on("(Q|mn)");
        for (int index = 0; index<ntri_naive_; index+=max_cols)
    {
            int cols = max_cols;
            if (index+cols>=ntri_naive_)
                cols = ntri_naive_-index;

            if (options_.get_str("RI_FITTING_TYPE") == "RAW" || options_.get_str("RI_FITTING_TYPE") == "FINISHED") {

                #pragma omp parallel for schedule (static)
                for (int r = 0; r<naux_raw_; r++) {
                    C_DCOPY(cols,&(B_ia_P_[r][index]),1,&(Temp1[r][0]),1);
                }

                C_DGEMM('T','N',naux_fin_,cols,naux_raw_,1.0, Winv_[0], naux_fin_, Temp1[0], max_cols,0.0, &B_ia_P_[0][index],ntri_naive_);
            } else if (options_.get_str("RI_FITTING_TYPE") == "CHOLESKY") {

                #pragma omp parallel for schedule(static)
                for (int A = 0; A< naux_raw_; A++) {
                    C_DCOPY(cols, &B_ia_P_[A][index], 1, &Temp1[0][A], naux_raw_);
                }

                int info = C_DPOTRS('U',naux_raw_,cols,&Winv_[0][0],naux_raw_,&Temp1[0][0],naux_raw_);

                #pragma omp parallel for schedule(static)
                for (int A = 0; A< naux_raw_; A++) {
                    C_DCOPY(cols, &Temp1[0][A], naux_raw_, &B_ia_P_[A][index], 1);
                }

            }
    }
        timer_off("(Q|mn)");

        free_block(Temp1);

        //fprintf(outfile,"  (B|mn)");
        //print_mat(B_ia_P_,naux_fin_,ntri_naive_,outfile);
        timer_on("3-Sieve");

        if (three_index_cutoff>0.0) {
            int left =  0;
            int right = 0;
            bool negligible_col;
            for (right = 0; right<ntri_naive_; right++) {
                negligible_col = true;
                for (int Q = 0; Q<naux_fin_; Q++) {
                    if (fabs(B_ia_P_[Q][right])>three_index_cutoff) {
                        negligible_col = false;
                        break;
                    }
                }
                if (!negligible_col) {
                    C_DCOPY(naux_fin_,&B_ia_P_[0][right],ntri_naive_,&B_ia_P_[0][left],ntri_naive_);
                    ri_pair_mu_[left] = ri_pair_mu_[right];
                    ri_pair_nu_[left] = ri_pair_nu_[right];
                    left++;
                } else {
                    ntri_--;
                }
            }
        }
        timer_off("3-Sieve");
        //fprintf(outfile,"  B is here:\n");
        //print_mat(B_ia_P_, naux_fin_,ntri_naive_ ,outfile);
        //for (int i = 0; i<ntri_naive_; i++)
        //    fprintf(outfile,"  i = %d, mu = %d, nu = %d\n",i,ri_pair_mu_[i],ri_pair_nu_[i]);

        if (options_.get_bool("RI_SCF_SAVE"))
        {
            write_B();
        }
    }
    else if (df_storage_ == disk)
    {

        //Open the BJ file and prestripe it to avoid wrong block errors
        timer_on("(B|mn) Prestripe");
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        double *Prestripe = init_array(ntri_naive_);
        for (int Q = 0; Q < naux_fin_; Q++) {
                psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(Prestripe[0]),sizeof(double)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        }
        free(Prestripe);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        timer_off("(B|mn) Prestripe");

        int pass = 0;

        //fprintf(outfile, "  Striped"); fflush(outfile);

        //Get an ERI object for the AO three-index integrals
        IntegralFactory rifactory(basisset_, basisset_, ribasis_,zero);
        //Get a TEI for each thread
        const double **buffer = new const double*[nthread];
        shared_ptr<TwoBodyAOInt> *eri = new shared_ptr<TwoBodyAOInt>[nthread];
        for (int Q = 0; Q<nthread; Q++) {
            eri[Q] = shared_ptr<TwoBodyAOInt>(rifactory.eri());
            buffer[Q] = eri[Q]->buffer();
        }

        //Determine the maximum nubmer of functions and pairs in the AO basis
        int maxfun = 0;
        for (int m = 0; m<basisset_->nshell(); m++)
            if (maxfun<basisset_->shell(m)->nfunction())
                maxfun=basisset_->shell(m)->nfunction();
        int maxpairs = maxfun*maxfun;

        //Find maximum allowed memory block size
        unsigned long int max_cols = ((df_memory_-memJ)/((long)(naux_raw_+naux_fin_)));
        if (max_cols > ntri_naive_ + maxpairs - 1)
            max_cols = ntri_naive_ + maxpairs - 1;
        if (max_cols < maxpairs)
            max_cols = maxpairs; //Gotta give me something to work with

        //Max rows per block for iterations (for diagnostics)
        unsigned long int max_rows = ((df_memory_)/((long)(norbs*ndocc+ntri_naive_)));
        if (max_rows > naux_fin_)
            max_rows = naux_fin_;
        if (max_rows < 1)
            max_rows = 1; //Without a row, I can't work

       if (print_>1) {
            fprintf(outfile,"\n  DF MEMORY USAGE:\n");
            fprintf(outfile,"  Total memory:           %15ld doubles, %15ld bytes.\n",memory_/sizeof(double),memory_);
            fprintf(outfile,"  DF memory:              %15ld doubles, %15ld bytes.\n",df_memory_,df_memory_*sizeof(double));
            fprintf(outfile,"  3-index storage:        %15ld doubles, %15ld bytes.\n",memA,memA*sizeof(double));
            fprintf(outfile,"  Exchange size:          %15ld doubles, %15ld bytes.\n",ndocc*naux_fin_*(long)norbs,ndocc*naux_fin_*(long)norbs*sizeof(double));
            fprintf(outfile,"  J memory:               %15ld doubles, %15ld bytes.\n",memJ,memJ*sizeof(double));
            fprintf(outfile,"  Max cols (integrals):   %15ld columns.\n",max_cols);
            fprintf(outfile,"  Max rows (iterations):  %15ld rows.\n",max_rows);
            fprintf(outfile,"  Nso:                    %15d functions.\n",norbs);
            fprintf(outfile,"  Npairs:                 %15d pairs.\n",ntri_naive_);
            fprintf(outfile,"  Ndocc:                  %15d orbitals.\n",ndocc);
            fprintf(outfile,"  Naux raw:               %15d functions.\n",naux_raw_);
            fprintf(outfile,"  Naux finished:          %15d functions.\n",naux_fin_);
            fprintf(outfile,"\n");
            fflush(outfile);
        }

        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        //
        //                                   NEW DISK ALGORITHM (THREADED)
        //
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        //Indexing
        int l_index = 0; //Runing schwarz map index
        int disk_index = 0; //Runing disk offset (dense start)
        //How many physical blocks could there be (I'll be safe, just for malloc)?
        int nblocks = ntri_naive_/(max_cols-maxpairs)*2;

        //What MUNU does the block start on?
        int* block_starts = init_int_array(nblocks);
        //How many shell pairs is the block
        int* block_sizes = init_int_array(nblocks);
        //What schwarzed omuonu does the block start on?
        int* porous_starts = init_int_array(nblocks);
        //How many schwarzed function pairs is the block?
        int* porous_sizes = init_int_array(nblocks);

        //Lets be safe and make a MU backmap
        int* MU_map = init_int_array(basisset_->nshell()*(basisset_->nshell()+1)/2);

        {//<<< Drop out of scope

            //Set up block bookkeeping
            nblocks = 0;
            int MU,NU,mu,nu,nummu,numnu,omu,onu, index;
            int shell_pairs = 0;
            int fun_pairs = 0;
            int shell_tot = 0;
            int fun_tot = 0;
            int block = 0;
            for (MU=0, index = 0; MU < basisset_->nshell(); ++MU) {
                nummu = basisset_->shell(MU)->nfunction();
                for (NU=0; NU <= MU; ++NU) {
                    MU_map[shell_tot] = MU;
                    shell_pairs++;
                    shell_tot++;
                    numnu = basisset_->shell(NU)->nfunction();
                    if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                        for (mu=0 ; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                            for (nu=0; nu < numnu; ++nu) {
                                onu = basisset_->shell(NU)->function_index() + nu;
                                if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                    fun_pairs++;
                                    fun_tot++;
                                }
                            }
                        }
                    }
                    if (MU == basisset_->nshell() - 1 && NU == MU) {
                        //Last block end
                        nblocks++;
                        block_sizes[block] = shell_pairs;
                        porous_sizes[block] = fun_pairs;
                        block++;
                    } else if (fun_pairs+maxpairs > max_cols) {
                        //Out of space
                        nblocks++;
                        block_sizes[block] = shell_pairs;
                        porous_sizes[block] = fun_pairs;
                        block_starts[block+1] = shell_tot;
                        porous_starts[block+1] = fun_tot;
                        shell_pairs = 0;
                        fun_pairs = 0;
                        block++;
                    }
                }
            }
        } //back into scope

        //fprintf(outfile,"  Number of blocks %d\n",nblocks);
        //for (int block = 0; block<nblocks; block++) {
        //    fprintf(outfile,"  Block %d: block_starts = %d, block sizes = %d, porous starts = %d, porous sizes = %d\n",block,block_starts[block],block_sizes[block],porous_starts[block],porous_sizes[block]);
        //}
        //fflush(outfile);

        //Allocate the fitted and unfitted blocks for three-index integrals
        double ** Amn;
        if (options_.get_str("RI_FITTING_TYPE") == "CHOLESKY")
            Amn = block_matrix(max_cols,naux_raw_); //Raw integrals
        else
            Amn = block_matrix(naux_raw_,max_cols); //Raw integrals
        double **Bmn = block_matrix(naux_fin_,max_cols); //Fitted integrals

        //Loop over blocks of index MUNU
        for (int block = 0; block<nblocks; block++) {

            //Need all the indexes
            int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
            int index;
            //Compute all the integrals for this block and place in buffer Amn
            timer_on("(A|mn)");
            if (options_.get_str("RI_FITTING_TYPE") == "CHOLESKY") {
                #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
                for (int block_address = block_starts[block]; block_address<block_starts[block]+block_sizes[block]; block_address++) {
                    #ifdef _OPENMP
                        rank = omp_get_thread_num();
                    #endif
                    //Where are we in MU/NU? Using the canonical backmap definition.
                    MU = MU_map[block_address];
                    NU = block_address-MU*(MU+1)/2;

                    //fprintf(outfile,"  MU = %d, NU = %d, address = %d",MU,NU,block_address);

                    //How big is the shell pair?
                    nummu = basisset_->shell(MU)->nfunction();
                    numnu = basisset_->shell(NU)->nfunction();
                    //Gogo integrals! (A|mn)
                    if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                        for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
                            numP = ribasis_->shell(Pshell)->nfunction();
                            eri[rank]->compute_shell(MU, NU, Pshell, 0);
                            for (mu=0 ; mu < nummu; ++mu) {
                                omu = basisset_->shell(MU)->function_index() + mu;
                                for (nu=0; nu < numnu; ++nu) {
                                    onu = basisset_->shell(NU)->function_index() + nu;
                                    if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                        for (P=0; P < numP; ++P) {
                                            PHI = ribasis_->shell(Pshell)->function_index() + P;
                                            Amn[ri_back_map_[omu*(omu+1)/2+onu]-porous_starts[block]][PHI] = buffer[rank][mu*numnu*numP+nu*numP+P];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
                for (int block_address = block_starts[block]; block_address<block_starts[block]+block_sizes[block]; block_address++) {
                    #ifdef _OPENMP
                        rank = omp_get_thread_num();
                    #endif
                    //Where are we in MU/NU? Using the canonical backmap definition.
                    MU = MU_map[block_address];
                    NU = block_address-MU*(MU+1)/2;

                    //fprintf(outfile,"  MU = %d, NU = %d, address = %d",MU,NU,block_address);

                    //How big is the shell pair?
                    nummu = basisset_->shell(MU)->nfunction();
                    numnu = basisset_->shell(NU)->nfunction();
                    //Gogo integrals! (A|mn)
                    if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                        for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
                            numP = ribasis_->shell(Pshell)->nfunction();
                            eri[rank]->compute_shell(MU, NU, Pshell, 0);
                            for (mu=0 ; mu < nummu; ++mu) {
                                omu = basisset_->shell(MU)->function_index() + mu;
                                for (nu=0; nu < numnu; ++nu) {
                                    onu = basisset_->shell(NU)->function_index() + nu;
                                    if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                        for (P=0; P < numP; ++P) {
                                            PHI = ribasis_->shell(Pshell)->function_index() + P;
                                            Amn[PHI][ri_back_map_[omu*(omu+1)/2+onu]-porous_starts[block]] = buffer[rank][mu*numnu*numP+nu*numP+P];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            timer_off("(A|mn)");

            //Embed fitting metric
            //(B|mn) = (A|mn)Winv_BA
            timer_on("(Q|mn)");
            if (options_.get_str("RI_FITTING_TYPE") == "CHOLESKY") {

                int info = C_DPOTRS('U',naux_raw_,porous_sizes[block],&Winv_[0][0],naux_raw_,&Amn[0][0],naux_raw_);

                #pragma omp parallel for schedule(static)
                for (int A = 0; A< naux_raw_; A++) {
                    C_DCOPY(porous_sizes[block], &Amn[0][A], naux_raw_, &Bmn[A][0], 1);
                }
            }
            else {
                C_DGEMM('T','N',naux_fin_,porous_sizes[block],naux_raw_,1.0, Winv_[0], naux_fin_, Amn[0], max_cols,0.0, Bmn[0],max_cols);
            }
            timer_off("(Q|mn)");

            //Use the three index sieve to compact the tensor
            int dense_size = 0;
            timer_on("3-Sieve");
            if (three_index_cutoff > 0.0) {
                for (int pair = 0; pair < porous_sizes[block]; pair++) {
                    bool sig = false;
                    for (int Q = 0; Q<naux_fin_;Q++) {
                        if (fabs(Bmn[Q][pair])>=three_index_cutoff) {
                             sig = true;
                             break;
                        }
                    }
                    if (sig) {
                        C_DCOPY(naux_fin_,&Bmn[0][pair],max_cols,&Amn[0][0]+dense_size,max_cols);
                        dense_size++;
                        ri_pair_mu_[l_index] = ri_pair_mu_[porous_starts[block]+pair];
                        ri_pair_nu_[l_index] = ri_pair_nu_[porous_starts[block]+pair];
                        l_index++;
                    } else {
                        ntri_--;
                    }
                }
            } else {
                dense_size = porous_sizes[block];
            }
            timer_off("3-Sieve");

            //Flush the bastards
            //Write the tensor out with the correct striping (mn is the fast index, but we only have blocks)
            //NOTE: If three_index_cutoff > 0.0, the tensor is in Amn, otherwise it is in Bmn
            //This allows for threading of the three index sieve
            timer_on("(Q|mn) Write");
            for (int Q = 0; Q < naux_fin_; Q++) {
                next_PSIF_DFSCF_BJ = psio_get_address(PSIO_ZERO,(ULI)(Q*(ULI)ntri_naive_*sizeof(double)+disk_index*sizeof(double)));
                psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *)((three_index_cutoff > 0.0)?&Amn[0][0]+Q*(ULI)ntri_naive_:&Bmn[Q][0]),sizeof(double)*dense_size,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            }
            timer_off("(Q|mn) Write");

            //Update indexing
            disk_index+=dense_size;
            //fprintf(outfile,"  disk index = %d, dense size = %d, porous sizes = %d, block = %d\n",disk_index,dense_size,porous_sizes[block],block);
        }

        //frees
        delete []buffer;
        delete []eri;
        free_block(Amn);
        free_block(Bmn);
        free(porous_sizes);
        free(porous_starts);
        free(block_sizes);
        free(block_starts);
        free(MU_map);

        if (print_>1)
            fprintf(outfile,"\n  Through DF integrals on disk.\n");
        fflush(outfile);

        //fprintf(outfile,"  B is here:\n");
        //next_PSIF_DFSCF_BJ = PSIO_ZERO;
        //double** Bhack = block_matrix(naux_fin_,ntri_naive_);
        //psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *)&Bhack[0][0],sizeof(double)*ntri_naive_*naux_fin_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        //print_mat(Bhack,naux_fin_,ntri_naive_,outfile);
        //free(Bhack);

        //fprintf(outfile,"  Indexing-fun pairs:\n");
        //for (int i = 0; i<ntri_naive_; i++)
        //    fprintf(outfile,"  i = %d, mu = %d, nu = %d, munu = %d\n",i,ri_pair_mu_[i],ri_pair_nu_[i], ri_pair_mu_[i]*(ri_pair_mu_[i]+1)/2+ri_pair_nu_[i]);
        //fprintf(outfile,"  Indexing-back map:\n");
        //for (int i = 0; i<norbs*(norbs+1)/2; i++)
        //    fprintf(outfile,"  munu = %d, index = %d\n",i,ri_back_map_[i]);


        //Write the restart data, it's cheap
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->write(PSIF_DFSCF_BJ,"RI_PAIR_BACK",(char *) &(ri_back_map_[0]),sizeof(int)*(norbs*(norbs+1)/2),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->write(PSIF_DFSCF_BJ,"RI_PAIR_MU",(char *) &(ri_pair_mu_[0]),sizeof(int)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->write(PSIF_DFSCF_BJ,"RI_PAIR_NU",(char *) &(ri_pair_nu_[0]),sizeof(int)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->write(PSIF_DFSCF_BJ,"N_TRI",(char *) &(ntri_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->write(PSIF_DFSCF_BJ,"N_TRI_NAIVE",(char *) &(ntri_naive_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->write(PSIF_DFSCF_BJ,"N_AUX_RAW",(char *) &(naux_raw_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->write(PSIF_DFSCF_BJ,"N_AUX_FIN",(char *) &(naux_fin_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        psio_->close(PSIF_DFSCF_BJ,1); //We'll reuse this methinks
    }
    timer_off("Overall (Q|mn)");

    free(schwarz_fun_pairs);
    free(schwarz_shell_pairs);
    free_block(Winv_);

    if (print_>1) {
        if (schwarz_) {
            fprintf(outfile,"\n  Function Pair Schwarz Sieve, Cutoff = %14.10E:\n",schwarz_);
            fprintf(outfile,"  %d out of %d basis function pairs removed, %8.5f%% attenuation.\n",norbs*(norbs+1)/2-sig_fun_pairs,norbs*(norbs+1)/2,100.0*(norbs*(norbs+1)/2-sig_fun_pairs)/(1.0*norbs*(norbs+1)/2));
            int pairs = basisset_->nshell()*(basisset_->nshell()+1)/2;
            fprintf(outfile,"  %d out of %d basis shell pairs removed, %8.5f%% attenuation.\n",pairs-sig_shell_pairs,pairs,100.0*(pairs-sig_shell_pairs)/(1.0*pairs));
        }
        if (three_index_cutoff) {
            int attenuation = ntri_naive_-ntri_;
            fprintf(outfile,"  Direct Three-Index Tensor Sieve, Cutoff = %14.10E:\n",three_index_cutoff);
            fprintf(outfile,"  %d of %d (remaining) basis function pairs removed, %8.5f%% attenuation.\n",attenuation,ntri_naive_, 100.0*attenuation/(1.0*ntri_naive_));
        }
        if (schwarz_>0.0 || three_index_cutoff>0.0)
            fprintf(outfile,"  After sieving, %d out of %d basis function pairs remain, %8.5f%% attenuation.\n\n",ntri_,norbs*(norbs+1)/2,100.0*(1.0-ntri_/(1.0*norbs*(norbs+1)/2)));
    }
    fflush(outfile);
}
void HF::write_B()
{
    if (print_>1)
        fprintf(outfile,"  Saving three-index tensor to file 97 for restart purposes.\n");

    psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
    psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;

    psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(B_ia_P_[0][0]),sizeof(double)*ntri_naive_*naux_fin_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);

    int norbs = basisset_->nbf();

    next_PSIF_DFSCF_BJ = PSIO_ZERO;
    psio_->write(PSIF_DFSCF_BJ,"RI_PAIR_BACK",(char *) &(ri_back_map_[0]),sizeof(int)*(norbs*(norbs+1)/2),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    next_PSIF_DFSCF_BJ = PSIO_ZERO;
    psio_->write(PSIF_DFSCF_BJ,"RI_PAIR_MU",(char *) &(ri_pair_mu_[0]),sizeof(int)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    next_PSIF_DFSCF_BJ = PSIO_ZERO;
    psio_->write(PSIF_DFSCF_BJ,"RI_PAIR_NU",(char *) &(ri_pair_nu_[0]),sizeof(int)*ntri_naive_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    next_PSIF_DFSCF_BJ = PSIO_ZERO;
    psio_->write(PSIF_DFSCF_BJ,"N_TRI",(char *) &(ntri_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    next_PSIF_DFSCF_BJ = PSIO_ZERO;
    psio_->write(PSIF_DFSCF_BJ,"N_TRI_NAIVE",(char *) &(ntri_naive_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    next_PSIF_DFSCF_BJ = PSIO_ZERO;
    psio_->write(PSIF_DFSCF_BJ,"N_AUX_RAW",(char *) &(naux_raw_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    next_PSIF_DFSCF_BJ = PSIO_ZERO;
    psio_->write(PSIF_DFSCF_BJ,"N_AUX_FIN",(char *) &(naux_fin_),sizeof(int),next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    psio_->close(PSIF_DFSCF_BJ,1);
}
void HF::free_B()
{
    if (df_storage_ == core)
        free_block(B_ia_P_);
    free(ri_pair_mu_);
    free(ri_pair_nu_);
    free(ri_back_map_);
}
void HF::form_W()
{
    //It takes a lot of work to get a null basis with Psi4!
    shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    //if (options_.get_bool("NO_INPUT") == false) {
    //    ribasis_ =shared_ptr<BasisSet>(new BasisSet(chkpt_, "DF_BASIS_SCF"));
    //} else {
        shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
        ribasis_ = BasisSet::construct(parser, molecule_, options_.get_str("RI_BASIS_SCF"));
    //}
    naux_raw_ = ribasis_->nbf();
    naux_fin_ = naux_raw_;

    if (print_>5) {
        basisset_->print(outfile);
        ribasis_->print(outfile);
        fflush(outfile);
    }

    timer_on("Form J");
    // J Matrix Time!
    W_ = block_matrix(naux_raw_, naux_raw_);

    //Next form J in the raw basis,
    IntegralFactory rifactory_J(ribasis_, zero, ribasis_, zero);
    shared_ptr<TwoBodyAOInt> Jint(rifactory_J.eri());

    // Integral buffer
    const double *Jbuffer = Jint->buffer();

    // J_{MN} = (0M|N0)
    int index = 0;

    for (int MU=0; MU < ribasis_->nshell(); ++MU) {
        int nummu = ribasis_->shell(MU)->nfunction();

        for (int NU=0; NU <= MU; ++NU) {
            int numnu = ribasis_->shell(NU)->nfunction();

            Jint->compute_shell(MU, 0, NU, 0);

            index = 0;
            for (int mu=0; mu < nummu; ++mu) {
                int omu = ribasis_->shell(MU)->function_index() + mu;

                for (int nu=0; nu < numnu; ++nu, ++index) {
                    int onu = ribasis_->shell(NU)->function_index() + nu;

                    W_[omu][onu] = Jbuffer[index];
                    W_[onu][omu] = Jbuffer[index];
                }
            }
        }
    }
    if (print_>5) {
        fprintf(outfile,"\nJ (raw):\n");
        print_mat(W_,naux_raw_,naux_raw_,outfile);
        fflush(outfile);
    }

    timer_off("Form J");
}
void HF::form_W_Poisson()
{
    //It takes a lot of work to get a null basis with Psi4!
    shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
//    if (options_.get_bool("NO_INPUT") == false) {
//        ribasis_ =shared_ptr<BasisSet>(new BasisSet(chkpt_, "DF_BASIS_SCF"));
//        poissonbasis_ =shared_ptr<BasisSet>(new BasisSet(chkpt_, "POISSON_BASIS_SCF"));
//    } else {
        shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
        ribasis_ = BasisSet::construct(parser, molecule_, options_.get_str("RI_BASIS_SCF"));
        poissonbasis_ = BasisSet::construct(parser, molecule_, options_.get_str("POISSON_BASIS_SCF"));
//    }
    naux_raw_ = ribasis_->nbf() + poissonbasis_->nbf();
    naux_fin_ = naux_raw_;

    int nri = ribasis_->nbf();
    int npois = poissonbasis_->nbf();

    if (print_>5) {
        basisset_->print(outfile);
        ribasis_->print(outfile);
        poissonbasis_->print(outfile);
        fflush(outfile);
    }

    timer_on("Form J");
    // J Matrix Time!
    W_ = block_matrix(naux_raw_, naux_raw_);

    // == (A|B) Block == //
    IntegralFactory rifactory_J(ribasis_, zero, ribasis_, zero);
    shared_ptr<TwoBodyAOInt> Jint(rifactory_J.eri());

    // Integral buffer
    const double *Jbuffer = Jint->buffer();

    // J_{MN} = (0M|N0)
    int index = 0;

    for (int MU=0; MU < ribasis_->nshell(); ++MU) {
        int nummu = ribasis_->shell(MU)->nfunction();

        for (int NU=0; NU <= MU; ++NU) {
            int numnu = ribasis_->shell(NU)->nfunction();

            Jint->compute_shell(MU, 0, NU, 0);

            index = 0;
            for (int mu=0; mu < nummu; ++mu) {
                int omu = ribasis_->shell(MU)->function_index() + mu;

                for (int nu=0; nu < numnu; ++nu, ++index) {
                    int onu = ribasis_->shell(NU)->function_index() + nu;

                    W_[omu][onu] = Jbuffer[index];
                    W_[onu][omu] = Jbuffer[index];
                }
            }
        }
    }

    // == (AB) Block == //
    IntegralFactory rifactory_RP(ribasis_, poissonbasis_,  zero, zero);
    shared_ptr<OneBodySOInt> O(rifactory_RP.so_overlap());
    SharedMatrix Omat(new Matrix(1, &nri, &npois));
    O->compute(Omat);

    for (int omu = 0; omu < ribasis_->nbf(); omu++) {
        for (int onu = ribasis_->nbf(); onu < ribasis_->nbf() + npois; onu++) {
            W_[omu][onu] = Omat->get(0,omu, onu - ribasis_->nbf());
            W_[onu][omu] = W_[omu][onu];
        }
    }

    // == (A|T|B) Block == //
    IntegralFactory rifactory_P(poissonbasis_, poissonbasis_,  zero, zero);
    shared_ptr<OneBodySOInt> T(rifactory_P.so_kinetic());
    SharedMatrix Tmat(new Matrix(1, &npois, &npois));
    T->compute(Tmat);

    for (int omu = ribasis_->nbf(); omu < ribasis_->nbf() + npois; omu++) {
        for (int onu = ribasis_->nbf(); onu < ribasis_->nbf() + npois; onu++) {
            W_[omu][onu] = 1.0/(2.0*M_PI)*Tmat->get(0,omu - ribasis_->nbf(), onu - ribasis_->nbf());
            W_[onu][omu] = W_[omu][onu];
        }
    }


    if (print_>5) {
        fprintf(outfile,"\nJ (raw):\n");
        print_mat(W_,naux_raw_,naux_raw_,outfile);
        fflush(outfile);
    }

    timer_off("Form J");
}
void HF::form_Wp12_chol()
{
    if (print_>1)
        fprintf(outfile,"  Fitting metric conditioned by Cholesky decomposition.\n");
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                   FORM CHOL(J^+1/2)
    //            NEW ALGORITHM: CHOLESKY SOLVER
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // Build RI basis, indexing, and W_
    if (options_.get_str("SCF_TYPE") == "POISSON")
        form_W_Poisson();
    else
        form_W();

    timer_on("Form J^+1/2");

    // Form J^-1/2
    // First, diagonalize J
    // the C_DSYEV call replaces the original matrix J with its eigenvectors
    double* eigval = init_array(naux_raw_);
    int lwork = naux_raw_ * 3;
    double* work = init_array(lwork);
    int stat = C_DSYEV('v','u',naux_raw_,W_[0],naux_raw_,eigval, work,lwork);
    if (stat != 0) {
        fprintf(outfile, "C_DSYEV failed\n");
        exit(PSI_RETURN_FAILURE);
    }
    free(work);

    // Now Jp contains the eigenvectors of the original Jp
    // Copy Jp to Jp_copy
    double **J_copy = block_matrix(naux_raw_, naux_raw_);
    C_DCOPY(naux_raw_*naux_raw_,W_[0],1,J_copy[0],1);

    // Now form Jp^{-1/2} = U(T)*j'^{-1/2}*U,
    // where j'^{-1/2} is the diagonal matrix of the inverse square roots
    // of the eigenvalues, and U is the matrix of eigenvectors of J'
    double min_J = eigval[0];
    double max_J = eigval[naux_raw_-1];

    #pragma omp parallel for schedule (static)
    for (int ind=0; ind<naux_raw_; ind++) {
        //All eignvalues should be approved for use!
        eigval[ind] = sqrt(eigval[ind]);
        // scale one set of eigenvectors by the diagonal elements jp^{-1/2}
        C_DSCAL(naux_raw_, eigval[ind], W_[ind], 1);
    }
    free(eigval);

    if (print_>1) {
        fprintf(outfile,"  Smallest eigenvalue in the raw fitting metric is %14.10E.\n",min_J);
        fprintf(outfile,"  Largest eigenvalue in the raw fitting metric is %14.10E.\n",max_J);
        fprintf(outfile,"  Condition number of the raw fitting metric is %14.10E.\n",max_J/min_J);
        fflush(outfile);
    }

    // J'^-1/2 = Jp_copy(T) * Jp
    Winv_ = block_matrix(naux_raw_,naux_raw_);
    C_DGEMM('t','n',naux_raw_,naux_raw_,naux_raw_,1.0,
            J_copy[0],naux_raw_,W_[0],naux_raw_,0.0,Winv_[0],naux_raw_);

    //fprintf(outfile,"W^+1/2:\n");
    //print_mat(Winv_,naux_raw_,naux_raw_,outfile);
    //fflush(outfile);

    free_block(J_copy);
    free_block(W_);

    timer_off("Form J^+1/2");

    timer_on("J Cholesky");

    //Choleskify Winv_ (results returned in place, should not be touched)
    stat = C_DPOTRF('U',naux_raw_,Winv_[0],naux_raw_);
    if (stat < 0) {
        fprintf(outfile, "Cholesky argument %d is bad.\n",stat);
        exit(PSI_RETURN_FAILURE);
    }
    else if (stat > 0) {
        fprintf(outfile, "W^+1/2 is not positive definite.\n");
        exit(PSI_RETURN_FAILURE);
    }

    if (debug_) {
        fprintf(outfile,"  L:\n");
        print_mat(Winv_,naux_raw_,naux_fin_,outfile);
    }

    timer_off("J Cholesky");

}
void HF::form_Wm12_raw()
{
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                      FORM J^-1/2
    //               OLD ALGORITHM: RAW FITTING
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // Build RI basis, indexing, and W_
    if (options_.get_str("SCF_TYPE") == "POISSON")
        form_W_Poisson();
    else
        form_W();

    timer_on("Form J^+1/2");

    // Form J^-1/2
    // First, diagonalize J
    // the C_DSYEV call replaces the original matrix J with its eigenvectors
    double* eigval = init_array(naux_raw_);
    int lwork = naux_raw_ * 3;
    double* work = init_array(lwork);
    int stat = C_DSYEV('v','u',naux_raw_,W_[0],naux_raw_,eigval, work,lwork);
    if (stat != 0) {
        fprintf(outfile, "C_DSYEV failed\n");
        exit(PSI_RETURN_FAILURE);
    }
    free(work);

    // Now Jp contains the eigenvectors of the original Jp
    // Copy Jp to Jp_copy
    double **J_copy = block_matrix(naux_raw_, naux_raw_);
    C_DCOPY(naux_raw_*naux_raw_,W_[0],1,J_copy[0],1);

    // Now form Jp^{-1/2} = U(T)*j'^{-1/2}*U,
    // where j'^{-1/2} is the diagonal matrix of the inverse square roots
    // of the eigenvalues, and U is the matrix of eigenvectors of J'
    double min_J = eigval[0];
    double max_J = eigval[naux_raw_-1];

    #pragma omp parallel for schedule (static)
    for (int ind=0; ind<naux_raw_; ind++) {
        if (eigval[ind] < 1.0E-10)
            eigval[ind] = 0.0;
        else {
            eigval[ind] = 1.0 / sqrt(eigval[ind]);
        }
        // scale one set of eigenvectors by the diagonal elements j^{-1/2}
        C_DSCAL(naux_raw_, eigval[ind], W_[ind], 1);
    }
    free(eigval);

    if (print_>1) {
        fprintf(outfile,"  Smallest eigenvalue in the raw fitting metric is %14.10E.\n",min_J);
        fprintf(outfile,"  Largest eigenvalue in the raw fitting metric is %14.10E.\n",max_J);
        fprintf(outfile,"  Condition number of the raw fitting metric is %14.10E.\n",max_J/min_J);
        fflush(outfile);
    }

    // J'^-1/2 = Jp_copy(T) * Jp
    Winv_ = block_matrix(naux_raw_,naux_raw_);
    C_DGEMM('t','n',naux_raw_,naux_raw_,naux_raw_,1.0,
            J_copy[0],naux_raw_,W_[0],naux_raw_,0.0,Winv_[0],naux_raw_);

    //fprintf(outfile,"W^+1/2:\n");
    //print_mat(Winv_,naux_raw_,naux_raw_,outfile);
    //fflush(outfile);

    free_block(J_copy);
    free_block(W_);

    timer_off("Form J^+1/2");
}
void HF::form_Wm12_fin()
{
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                      FORM J^-1/2
    //        NEW ALGORITHM: CANONICAL ORTHOGONALIZATION
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // Build RI basis, indexing, and W_
    if (options_.get_str("SCF_TYPE") == "POISSON") {
        //form_W_Poisson();
        throw std::domain_error("Poisson fitting with canonical orthogonalization not implemented yet");
    } else
        form_W();

    timer_on("Form X");
    //OK, integrals time. start with the fitting matrix J

    // Create integral factory for J_S (Overlap Matrix in form_B)
    shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    IntegralFactory rifactory_JS(ribasis_, ribasis_, zero, zero);
    OneBodySOInt *S = rifactory_JS.so_overlap();
    MatrixFactory matJS;
    matJS.init_with(1,&naux_raw_,&naux_raw_);

    //Put the integrals in a good old double**
    double** SJ;

    {
        SharedMatrix S_J = shared_ptr<Matrix>(matJS.create_matrix("S_J"));
        //Compute those integrals
        S->compute(S_J);
        delete S;
        //S_J->print(outfile);
        SJ = S_J->to_block_matrix();
    }


    // Form X
    // First, diagonalize S_J
    // the C_DSYEV call replaces the original matrix SJ with its eigenvectors
    double* s_eigval = init_array(naux_raw_);
    int lswork = naux_raw_ * 3;
    double* swork = init_array(lswork);
    int stat = C_DSYEV('v','u',naux_raw_,SJ[0],naux_raw_,s_eigval, swork,lswork);
    if (stat != 0) {
        fprintf(outfile, "C_DSYEV failed\n");
        exit(PSI_RETURN_FAILURE);
    }
    free(swork);

    //Find the condition number and largest/smallest eigenvalues of SJ
    //The largest one is what is important
    double min_SJ = s_eigval[0];
    double max_SJ = s_eigval[naux_raw_-1];

    //Print some stuff
    double cond_SJ = max_SJ/min_SJ;
    if (print_>1) {
        fprintf(outfile,"  Min J overlap eigenvalue is %14.10E.\n",min_SJ);
        fprintf(outfile,"  Max J overlap eigenvalue is %14.10E.\n",max_SJ);
        fprintf(outfile,"  J overlap condition number is %14.10E.\n",cond_SJ);
        fflush(outfile);
    }

    //Find the user's allowed condition number
    double cond_SJ_tol = options_.get_double("RI_MAX_COND");
    //and compute the smallest allowed eigenvlue in the overlap matrix
    double cutoff_SJ = max_SJ/cond_SJ_tol;

    //Go through and scale columns of eigenvectors by allowed s^-1/2
    int linear_dependencies = 0;
    double new_min_SJ = s_eigval[naux_raw_-1]; //should be the biggest
    for (int i = 0; i<naux_raw_; i++) {
        if (s_eigval[i] <= cutoff_SJ) {
            //OK this eigenvalue does not contribute to the subspace completeness
            //and moreover it kills my condition number
            s_eigval[i] = 0.0; //This isn't officially a big deal, but it's proper etiquette
            naux_fin_--;
            linear_dependencies++;
        } else {
            if (new_min_SJ > s_eigval[i])
                new_min_SJ = s_eigval[i];
            //OK, this eigenvalue is approved for launch
            s_eigval[i] = 1.0 / sqrt(s_eigval[i]);
        }
    }

    int ind;
    #pragma omp parallel for private (ind) schedule (dynamic)
    for (int ind = 0; ind<naux_raw_; ind++) {

        //Scale the column out to produce the porous X matrix in place
        //(Columns are rows in this messed up LAPACK-ness
        if (s_eigval[ind] > 0.0)
            C_DSCAL(naux_raw_, s_eigval[ind], SJ[ind], 1);
    }

    //fprintf(outfile,"X (in S ordering):\n");
    //print_mat(SJ,naux_raw_,naux_raw_,outfile);
    //fflush(outfile);


    //So how'd that go?
    if (print_>1) {
        fprintf(outfile, "  %d of %d finished auxiliary functions eliminated by canonical orthogonalization.\n",linear_dependencies, naux_raw_);
        fprintf(outfile, "  Resultant condition number in J overlap matrix is %14.10E.\n",max_SJ/new_min_SJ);
        fflush(outfile);
    }

    //The famed transformation matrix
    double **X = block_matrix(naux_raw_,naux_fin_);

    //We'll put it in canonical (descending eigenvalue) order
    //This is after all, canonical. Well...as canonical as
    //representing one subspace with another subspace cast
    //onto a third subspace gets
    #pragma omp parallel for
    for (int m = 0; m<naux_raw_; m++) {
        for (int n = 0; n<naux_fin_; n++) {
            X[m][n] = SJ[naux_raw_-n-1][m];
        }
    }

    //fprintf(outfile,"X:\n");
    //print_mat(X,naux_raw_,naux_raw_,outfile);
    //fflush(outfile);

    free(s_eigval);
    free_block(SJ);

    if (print_>5) {
        fprintf(outfile,"  X (The J -> J' change of basis matrix):\n");
        print_mat(X,naux_raw_,naux_fin_,outfile);
        fflush(outfile);
    }

    timer_off("Form X");

    //fprintf(outfile,"J:\n");
    //print_mat(J,naux_raw_,naux_raw_,outfile);
    //fflush(outfile);

    //Tansform the J matrix to J'
    double **Jtemp1 = block_matrix(naux_fin_,naux_raw_);
    double **Jp = block_matrix(naux_fin_,naux_fin_);

    //Temp = X'*J
    C_DGEMM('T','N',naux_fin_,naux_raw_,naux_raw_,1.0,X[0],naux_fin_,W_[0],naux_raw_,0.0,Jtemp1[0],naux_raw_);

    //J' = Temp*X
    C_DGEMM('N','N',naux_fin_,naux_fin_,naux_raw_,1.0,Jtemp1[0],naux_raw_,X[0],naux_fin_,0.0,Jp[0],naux_fin_);

    //fprintf(outfile,"J':\n");
    //print_mat(Jp,naux_fin_,naux_fin_,outfile);
    //fflush(outfile);

    free_block(Jtemp1);

    //For diagnostics purposes, we may want to know
    //how bad the raw J was
    if (options_.get_bool("FIND_RAW_J_COND")) {
        double** J2 = block_matrix(naux_raw_,naux_raw_);
        C_DCOPY(naux_raw_*naux_raw_,W_[0],1,J2[0],1);

        //Find eigenvalues only
        double* eigval = init_array(naux_raw_);
        int lwork = naux_raw_ * 3;
        double* work = init_array(lwork);
        stat = C_DSYEV('n','u',naux_raw_,J2[0],naux_raw_,eigval, work,lwork);
        if (stat != 0) {
            fprintf(outfile, "C_DSYEV failed\n");
            exit(PSI_RETURN_FAILURE);
        }

        free(work);

        double min_J = eigval[0];
        double max_J = eigval[naux_raw_-1];
        free(eigval);

        fprintf(outfile,"  Smallest eigenvalue in the raw fitting metric is %14.10E.\n",min_J);
        fprintf(outfile,"  Largest eigenvalue in the raw fitting metric is %14.10E.\n",max_J);
        fprintf(outfile,"  Condition number of the raw fitting metric is %14.10E.\n",max_J/min_J);
        fflush(outfile);

        free_block(J2);
    }
    free_block(W_);

    //Find J'^-1/2
    timer_on("Form J^-1/2");

    //double **J_copy2 = block_matrix(naux_fin_, naux_fin_);
    //C_DCOPY(naux_fin_*naux_fin_,Jp[0],1,J_copy2[0],1);

    // Form Jp^-1/2
    // First, diagonalize Jp
    // the C_DSYEV call replaces the original matrix Jp with its eigenvectors
    double* eigval = init_array(naux_fin_);
    int lwork = naux_fin_ * 3;
    double* work = init_array(lwork);
    stat = C_DSYEV('v','u',naux_fin_,Jp[0],naux_fin_,eigval, work,lwork);
    if (stat != 0) {
        fprintf(outfile, "C_DSYEV failed\n");
        exit(PSI_RETURN_FAILURE);
    }
    free(work);

    // Now Jp contains the eigenvectors of the original Jp
    // Copy Jp to Jp_copy
    double **Jp_copy = block_matrix(naux_fin_, naux_fin_);
    C_DCOPY(naux_fin_*naux_fin_,Jp[0],1,Jp_copy[0],1);

    // Now form Jp^{-1/2} = U(T)*j'^{-1/2}*U,
    // where j'^{-1/2} is the diagonal matrix of the inverse square roots
    // of the eigenvalues, and U is the matrix of eigenvectors of J'
    double min_J = eigval[0];
    double max_J = eigval[naux_fin_-1];

    #pragma omp parallel for private (ind) schedule (static)
    for (ind=0; ind<naux_fin_; ind++) {
        //All eignvalues should be approved for use!
        eigval[ind] = 1.0 / sqrt(eigval[ind]);
        // scale one set of eigenvectors by the diagonal elements jp^{-1/2}
        C_DSCAL(naux_fin_, eigval[ind], Jp[ind], 1);
    }
    free(eigval);

    if (print_>1) {
        fprintf(outfile,"  Smallest eigenvalue in the finished fitting metric is %14.10E.\n",min_J);
        fprintf(outfile,"  Largest eigenvalue in the finished fitting metric is %14.10E.\n",max_J);
        fprintf(outfile,"  Condition number of the finished fitting metric is %14.10E.\n",max_J/min_J);
        fflush(outfile);
    }

    // J'^-1/2 = Jp_copy(T) * Jp
    double** Jp_mhalf = block_matrix(naux_fin_,naux_fin_);
    C_DGEMM('t','n',naux_fin_,naux_fin_,naux_fin_,1.0,
            Jp_copy[0],naux_fin_,Jp[0],naux_fin_,0.0,Jp_mhalf[0],naux_fin_);

    //fprintf(outfile,"J'^-1/2:\n");
    //print_mat(Jp_mhalf,naux_fin_,naux_fin_,outfile);
    //fflush(outfile);

    free_block(Jp_copy);
    free_block(Jp);

    timer_off("Form J^-1/2");
    timer_on("Form W");

    //Winv_AB' = X_AC'*J_CB'^-1/2
    Winv_ = block_matrix(naux_raw_,naux_fin_);

    C_DGEMM('N','N',naux_raw_,naux_fin_,naux_fin_,1.0,X[0],naux_fin_,Jp_mhalf[0],naux_fin_,0.0,Winv_[0],naux_fin_);

    //fprintf(outfile,"W:\n");
    //print_mat(W,naux_raw_,naux_fin_,outfile);
    //fflush(outfile);
    free_block(Jp_mhalf);
    free_block(X);
    timer_off("Form W");
}
void RHF::form_G_from_RI()
{
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                        FORM G
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    timer_on("Overall G");
    //Get norbs
    int norbs = basisset_->nbf();
    //Zero the J matrix
    if (J_is_required_)
        J_->zero();
    //Zero the K matrix
    if (K_is_required_)
        K_->zero();

    //D_->print(outfile);
    //C_->print(outfile);

    //Rearrange the D matrix as a vector in terms of ri_pair indices
    //Off diagonal elements get 2x weight due to permutational symmetry
    double* DD = init_array(ntri_);

    for (int ij = 0; ij<ntri_; ij++) {
        DD[ij] = D_->get(0,ri_pair_mu_[ij],ri_pair_nu_[ij]);
        if (ri_pair_mu_[ij] != ri_pair_nu_[ij])
            DD[ij] *= 2.0;
            //only A irrep at the moment!!
    }
    //Get the C matrix (exchange messes things up)
    int ndocc = doccpi_[0];
    if (iteration_ == 0 && options_.get_str("GUESS") == "SAD" && options_.get_str("SAD_C") == "CHOLESKY")
        ndocc = sad_nocc_[0];
    double** Cocc = block_matrix(ndocc,norbs);
    for (int i=0; i<norbs; i++) {
        for (int j=0; j<ndocc; j++)
            Cocc[j][i] = C_->get(0,i,j);
        //only A irrep at the moment!!
    }

    //fprintf(outfile,"  ndocc = %d\n",ndocc);
    //C_->print(outfile);

    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                        CORE ALGORITHM
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (df_storage_ == core) {
        //B is in core, DGEMM everything
        if (J_is_required_) {
            /* COULOMB PART */
            //Coulomb convolution vector
            double *L = init_array(naux_fin_);
            //Temporary J matrix
            double *J = init_array(ntri_);
            //DGEMV -> L:
            //L_Q = (Q|ls)*D_{ls}
            timer_on("J DDOT");
            C_DGEMV('N',naux_fin_,ntri_,1.0,B_ia_P_[0],ntri_naive_,DD,1,0.0,L,1);
            timer_off("J DDOT");
            //DGEMV -> J:
            //J_{mn} = L_Q(Q|mn)
            timer_on("J DAXPY");
            C_DGEMV('T',naux_fin_,ntri_,1.0,B_ia_P_[0],ntri_naive_,L,1,0.0,J,1);
            timer_off("J DAXPY");
            //Put everything in J_
            for (int ij = 0; ij < ntri_; ij++) {
                J_->set(0,ri_pair_mu_[ij],ri_pair_nu_[ij],J[ij]);
                J_->set(0,ri_pair_nu_[ij],ri_pair_mu_[ij],J[ij]);
            }

            free(L);
            free(J);
        }
        if (K_is_required_) {
            timer_on("Form E");
            /* EXCHANGE PART */
            unsigned long int max_rows = ((df_memory_-naux_raw_*(long)ntri_naive_)/((long)(norbs*ndocc)));
            if (max_rows > naux_fin_)
                max_rows = naux_fin_;
            if (max_rows < 1)
                max_rows = 1;

            //fprintf(outfile,"  max_rows = %d\n",max_rows);
            int nthread = 1;
            #ifdef _OPENMP
                nthread = omp_get_max_threads();
                //printf("Now OMP says %d threads\n",nthread);
            #endif

            //E exchange matrix
            double** E = block_matrix(norbs, ndocc*max_rows);
            //Temporary K matrix
            double** K = block_matrix(norbs,norbs);
            //QS temp matrix for DGEMM
            double*** QS = new double**[nthread];
            for (int T = 0; T < nthread; T++)
                QS[T] = block_matrix(max_rows,norbs);
            // Temp matrix for sparse DGEMM if sieve exists
            double*** Ctemp = new double**[nthread];
            for (int T = 0; T < nthread; T++)
                Ctemp[T] = block_matrix(ndocc,norbs);
            // Index array for non-canonical ordering of mn
            int** m_ij_indices = init_int_matrix(norbs,norbs);
            // Index array of n for given m (in order of above)
            int** n_indices = init_int_matrix(norbs,norbs);
            // sizes of above for schwarz sieve
            int* index_sizes = init_int_array(norbs);

            timer_on("Initial E Indexing");

            for (int ij = 0; ij<ntri_; ij++) {
                int m = ri_pair_mu_[ij];
                int n = ri_pair_nu_[ij];
                m_ij_indices[m][index_sizes[m]] = ij;
                n_indices[m][index_sizes[m]] = n;
                index_sizes[m]++;
                if (m != n){
                    m_ij_indices[n][index_sizes[n]] = ij;
                    n_indices[n][index_sizes[n]] = m;
                    index_sizes[n]++;
                }
            }

            timer_off("Initial E Indexing");

            //CUTOFF STUFF (KINDA LOOSE)
            //What is the smallest overlap element worth computing exchange for
            double cutoff = options_.get_double("OVERLAP_CUTOFF");
            //Contribution to the exchange matrix
            double contribution;
            //Partial contribution K matrix
            double** Ktemp;
            int att_elements = norbs*(norbs+1)/2;

            //MASTER LOOP
            int current_rows;
            for (int row = 0; row<naux_fin_; row += max_rows) {
                current_rows = max_rows;
                if (row+current_rows>=naux_fin_)
                    current_rows = naux_fin_-row;

                timer_on("E SORT");

                #ifdef _MKL
                    int mkl_nthreads = mkl_get_max_threads();
                    //printf("MKL says %d threads\n", mkl_nthreads);
                    mkl_set_num_threads(1);
                #endif

                int m, n , ij, index, rank;

                #pragma omp parallel for private (m, n , ij, index, rank) schedule (dynamic)
                for (m = 0; m<norbs; m++) {

                    rank = 0;
                    #ifdef _OPENMP
                        rank = omp_get_thread_num();
                    #endif

                    //Find out where the m's are!
                    //timer_on("Intermediate Indexing");
                    int n, ij;
                    for (index = 0; index<index_sizes[m]; index++) {
                        ij = m_ij_indices[m][index];
                        n = n_indices[m][index];

                        //fprintf(outfile,"  index = %d ij = %d \n",index,ij); fflush(outfile);//fprintf(outfile,"(ij, mu) = (%d, %d)\n",ij,mu);
                        C_DCOPY(current_rows,&B_ia_P_[row][ij],ntri_naive_,&QS[rank][0][index],norbs);
                        //for (int Q = 0; Q<naux_fin_; Q++) {
                        //     QS[Q][index] = B_ia_P_[Q][ij];
                             //fprintf(outfile," ij = %d, mu = %d, Q = %d, val = %14.10f\n",ij,mu,Q,QS[Q][mu]);
                        //}
                        C_DCOPY(ndocc,&Cocc[0][n],norbs,&Ctemp[rank][0][index],norbs);
                        //for (int o = 0; o<ndocc; o++) {
                        //    Ctemp[o][index] = Cocc[o][n];
                        //}
                    }
                    //timer_off("Intermediate Indexing");
                    //timer_on("DGEMM 2");

                    //E_im^Q = C_in (Q|mn). Watch the lda on the E, or you'll cut if off in the last block!
                    C_DGEMM('N','T',ndocc,current_rows,index_sizes[m],1.0,Ctemp[rank][0],norbs,QS[rank][0],norbs, 0.0, E[m], current_rows);
                    //timer_off("DGEMM 2");
                }

                #ifdef _MKL
                    mkl_set_num_threads(mkl_nthreads);
                #endif

                timer_off("E SORT");
                timer_on("E DGEMM");

                //K_{mn} = E_{im}^QE_{in}^Q

                if (cutoff > 0.0) {
                    for (int i = 0; i< norbs; i++)
                        for (int j = 0; j <= i; j++) {
                            //Is S big enough?
                            if (abs(S_->get(0,i,j))>=cutoff) {
                                contribution = C_DDOT(current_rows*ndocc,&E[i][0],1,&E[j][0],1);
                                K[i][j] += contribution;
                                if (i != j)
                                    K[j][i] += contribution;
                                if (row == 0)
                                    att_elements--;
                            }
                        }
                } else {
                    //DGEMM, usually faster regardless
                    //There it is
                    if (row == 0)
                        Ktemp = block_matrix(norbs,norbs);
                    //If the K matrix is sparse due to low overlap or density
                    //Form the matrix elements individually
                    C_DGEMM('N','T',norbs,norbs,current_rows*ndocc,1.0,E[0],max_rows*ndocc,E[0],max_rows*ndocc, 0.0, Ktemp[0], norbs);

                    C_DAXPY(norbs*norbs,1.0,&Ktemp[0][0],1,&K[0][0],1);
                }
                timer_off("E DGEMM");
            }

            //print_mat(K,norbs,norbs,outfile); fflush(outfile);

            for (int i = 0; i < norbs; i++)
                for (int j = 0; j<=i; j++) {
                    K_->set(0,j,i,K[i][j]);
                    K_->set(0,i,j,K[i][j]);
                }

            if (cutoff > 0.0 && print_>2) {
                fprintf(outfile,"  K matrix was built with with an overlap cutoff of %14.10E\n",cutoff);
                fprintf(outfile,"  %d out of %d elements were attenuated due to overlap, %8.5f%% attenuation.\n",att_elements,norbs*(norbs+1)/2,100.0*att_elements/(1.0*norbs*(norbs+1)/2));
            } else if (cutoff == 0.0) {
                free_block(Ktemp);
            }

            //Frees
            free_block(E);
            free_block(K);
            free(m_ij_indices[0]);
            free(m_ij_indices);
            free(n_indices[0]);
            free(n_indices);
            free(index_sizes);

            for (int T = 0; T<nthread; T++) {
                free_block(Ctemp[T]);
                free_block(QS[T]);
            }
            delete []Ctemp;
            delete []QS;

            timer_off("Form E");

            //fprintf(outfile,"\n E: \n");
            //print_mat(E,norbs,ndocc*naux_fin_,outfile);

            //K_->print(outfile);
            fflush(outfile);
        }
    } else {
        //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        //
        //                      DISK ALGORITHM
        //
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        //B is on disk, E will be done in blocks on core
        //B_ia_P_ stores multiple aux basis function rows of
        //the three index tensor

        //How many rows to read?
        unsigned long int max_rows = ((df_memory_)/((long)(norbs*ndocc+ntri_naive_)));
    if (max_rows > naux_fin_)
            max_rows = naux_fin_;
        if (max_rows < 1)
            max_rows = 1; //Without a row, I can't work

        //fprintf(outfile,"  max_rows = %d\n",max_rows); fflush(outfile);

        //FILE STUFF
        //Transformed integrals are stored in PSIF_DFSCF_BJ
        //Which had better exist at this point
        //timer_on("Open B");
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        //timer_off("Open B");
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        //Row height per read, so that we can tune this value
        int rows_per_read = options_.get_int("ROWS_PER_READ");
        if (rows_per_read == 0)
            rows_per_read = max_rows;

        //Chunk of three-index tensor
        B_ia_P_ = block_matrix(max_rows,ntri_naive_);

        //COULOMB STUFF
        //Coulomb convolution vector, done element by element
        double *L; //Density coefficients
        double *J; //J register
        double *Jtemp; //J contribution
        if (J_is_required_){
            L = init_array(max_rows);
            Jtemp = init_array(ntri_);
            J = init_array(ntri_);
        }

        int nthread = 1;
        #ifdef _OPENMP
            nthread = omp_get_max_threads();
        #endif
        //EXCHANGE STUFF
        double **E;
        double **K;
        double ***QS;
        double ***Ctemp;
        int** m_ij_indices;
        int** n_indices;
        int* index_sizes;
        if (K_is_required_) {
            //E exchange matrix
            E = block_matrix(norbs, ndocc*max_rows);
            //Temporary K matrix
            K = block_matrix(norbs,norbs);
            //QS temp matrix for DGEMM
            QS = new double**[nthread];
            for (int T = 0; T < nthread; T++)
                QS[T] = block_matrix(max_rows,norbs);
            // Temp matrix for sparse DGEMM if sieve exists
            Ctemp = new double**[nthread];
            for (int T = 0; T < nthread; T++)
                Ctemp[T] = block_matrix(ndocc,norbs);
            // Index array for non-canonical ordering of mn
            m_ij_indices = init_int_matrix(norbs,norbs);
            // Index array of n for given m (in order of above)
            n_indices = init_int_matrix(norbs,norbs);
            // sizes of above for schwarz sieve
            index_sizes = init_int_array(norbs);

            timer_on("Initial E Indexing");

            for (int ij = 0; ij<ntri_; ij++) {
                int m = ri_pair_mu_[ij];
                int n = ri_pair_nu_[ij];
                m_ij_indices[m][index_sizes[m]] = ij;
                n_indices[m][index_sizes[m]] = n;
                index_sizes[m]++;
                if (m != n){
                    m_ij_indices[n][index_sizes[n]] = ij;
                    n_indices[n][index_sizes[n]] = m;
                    index_sizes[n]++;
                }
            }

            timer_off("Initial E Indexing");
        }

        //CUTOFF STUFF (KINDA LOOSE)
        //What is the smallest overlap element worth computing exchange for
        double cutoff = options_.get_double("OVERLAP_CUTOFF");
        //Contribution to the exchange matrix
        double contribution;
        //Partial contribution K matrix
        double** Ktemp;
        int att_elements = norbs*(norbs+1)/2;

        //MASTER LOOP
        int current_rows,offset;
        //int pass = 0;
        for (int row = 0; row <naux_fin_; row+=max_rows)
        {

            //Setup block size
            //Careful Starfox
        current_rows = max_rows;
        if (row+max_rows>=naux_fin_)
        current_rows = naux_fin_-row;
            //Read max_rows of the (B|mn) tensor in, place in B_ia_P
            //fprintf(outfile,"\n  Pass %d, current_rows = %d\n",pass,current_rows);
            //pass++;
            timer_on("Read B");

            //New method read in a few rows at a time
            int block_height = rows_per_read;
            for (int Q = 0; Q<current_rows; Q+=rows_per_read) {
                if (Q+rows_per_read>=current_rows)
                    block_height = current_rows-Q;
                psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(B_ia_P_[Q][0]),sizeof(double)*ntri_naive_*(ULI)block_height,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);

            }

            //fprintf(outfile,"  Read pass = %d. Current rows = %d. Max rows = %d, Ntri_naive = %d, Ntri = %d.\n",pass++,current_rows,max_rows,ntri_naive_,ntri_);

            //double **Bhack = block_matrix(current_rows,ntri_);
            //C_DCOPY(current_rows*ntri_,in_buffer,1,&Bhack[0][0],1);
            //fprintf(outfile,"  B:\n");
            //print_mat(Bhack,current_rows,ntri_,outfile);

            timer_off("Read B");

            //Do J here with blocked DGEMV
            if (J_is_required_) {
                //DGEMV -> L:
                //L_Q = (Q|ls)*D_{ls}
                timer_on("J DDOT");
                C_DGEMV('N',current_rows,ntri_,1.0,B_ia_P_[0],ntri_naive_,DD,1,0.0,L,1);
                timer_off("J DDOT");
                //DGEMV -> J:
                //J_{mn} = L_Q(Q|mn)
                timer_on("J DAXPY");
                C_DGEMV('T',current_rows,ntri_,1.0,B_ia_P_[0],ntri_naive_,L,1,0.0,Jtemp,1);
                timer_off("J DAXPY");

                C_DAXPY(ntri_,1.0,Jtemp,1,J,1);

            }

            //Do K here just like the core algorithm
            if (K_is_required_) {
                timer_on("E SORT");

                #ifdef _MKL
                    int mkl_nthreads = mkl_get_max_threads();
                    mkl_set_num_threads(1);
                #endif

                int m, n , ij3, index, rank;

                #pragma omp parallel for private (m, n , ij3, index, rank) schedule (dynamic)
                for (m = 0; m<norbs; m++) {

                    rank = 0;
                    #ifdef _OPENMP
                        rank = omp_get_thread_num();
                    #endif
                    //Find out where the m's are!
                    //timer_on("Intermediate Indexing");
                    for (index = 0; index<index_sizes[m]; index++) {
                        ij3 = m_ij_indices[m][index];
                        n = n_indices[m][index];


                        //fprintf(outfile,"  index = %d ij = %d \n",index,ij); fflush(outfile);//fprintf(outfile,"(ij, mu) = (%d, %d)\n",ij,mu);
                        C_DCOPY(current_rows,&B_ia_P_[0][ij3],ntri_naive_,&QS[rank][0][index],norbs);
                        //for (int Q = 0; Q<naux_fin_; Q++) {
                        //     QS[Q][index] = B_ia_P_[Q][ij];
                             //fprintf(outfile," ij = %d, mu = %d, Q = %d, val = %14.10f\n",ij,mu,Q,QS[Q][mu]);
                        //}
                        C_DCOPY(ndocc,&Cocc[0][n],norbs,&Ctemp[rank][0][index],norbs);
                        //for (int o = 0; o<ndocc; o++) {
                        //    Ctemp[o][index] = Cocc[o][n];
                        //}
                    }
                    //timer_off("Intermediate Indexing");
                    //timer_on("DGEMM 2");

                    //E_im^Q = C_in (Q|mn). Watch the lda on the E, or you'll cut if off in the last block!
                    C_DGEMM('N','T',ndocc,current_rows,index_sizes[m],1.0,Ctemp[rank][0],norbs,QS[rank][0],norbs, 0.0, E[m], current_rows);
                    //timer_off("DGEMM 2");
                    //timer_on("Final Indexing");
                }
                #ifdef _MKL
                    mkl_set_num_threads(mkl_nthreads);
                #endif

                timer_off("E SORT");
                timer_on("E DGEMM");

                //K_{mn} = E_{im}^QE_{in}^Q

                if (cutoff > 0.0) {
                    for (int i = 0; i< norbs; i++)
                        for (int j = 0; j <= i; j++) {
                            //Is S big enough?
                            if (abs(S_->get(0,i,j))>=cutoff) {
                                contribution = C_DDOT(current_rows*ndocc,&E[i][0],1,&E[j][0],1);
                                K[i][j] += contribution;
                                if (i != j)
                                    K[j][i] += contribution;
                                if (row == 0)
                                    att_elements--;
                            }
                        }
                } else {
                    //DGEMM, usually faster regardless
                    //There it is
                    if (row == 0)
                        Ktemp = block_matrix(norbs,norbs);
                    //If the K matrix is sparse due to low overlap or density
                    //Form the matrix elements individually
                    C_DGEMM('N','T',norbs,norbs,current_rows*ndocc,1.0,E[0],max_rows*ndocc,E[0],max_rows*ndocc, 0.0, Ktemp[0], norbs);

                    C_DAXPY(norbs*norbs,1.0,&Ktemp[0][0],1,&K[0][0],1);
                }
                timer_off("E DGEMM");
            }
        }
        psio_->close(PSIF_DFSCF_BJ,1);
        free_block(B_ia_P_);
        /* Form J */
        if (J_is_required_) {
            for (int ij2 = 0; ij2 < ntri_; ij2++) {
                J_->set(0,ri_pair_mu_[ij2],ri_pair_nu_[ij2],J[ij2]);
                J_->set(0,ri_pair_nu_[ij2],ri_pair_mu_[ij2],J[ij2]);
            }
            free(J);
            free(Jtemp);
            free(L);
        }
        /* Form K */
        if (K_is_required_) {
            for (int i = 0; i < norbs; i++)
                for (int j = 0; j<=i; j++) {
                    K_->set(0,j,i,K[i][j]);
                    K_->set(0,i,j,K[i][j]);
            }
            if (cutoff > 0.0 && print_>2) {
                fprintf(outfile,"  K matrix was built with with an overlap cutoff of %14.10E\n",cutoff);
                fprintf(outfile,"  %d out of %d elements were attenuated due to overlap, %8.5f%% attenuation.\n",att_elements,norbs*(norbs+1)/2,100.0*att_elements/(1.0*norbs*(norbs+1)/2));
            } else if (cutoff == 0.0) {
                free_block(Ktemp);
            }

            //Frees
            free_block(E);
            free_block(K);
            free(m_ij_indices[0]);
            free(m_ij_indices);
            free(n_indices[0]);
            free(n_indices);
            free(index_sizes);
            for (int T = 0; T<nthread; T++) {
                free_block(Ctemp[T]);
                free_block(QS[T]);
            }
            delete []Ctemp;
            delete []QS;
        }
    }

    free(DD);
    free_block(Cocc);

    //J_->print(outfile);
    //K_->print(outfile);

    /* FORM G_ */
    //This method takes one extra O(N^2) scale,
    //but preserves J and K in place
    G_->zero();
    G_->add(J_);
    G_->scale(-2.0);
    if (K_is_required_)
        G_->add(K_);
    G_->scale(-1.0);
    //G_->print(outfile);
    timer_off("Overall G");

}
void RHF::form_J_from_RI()
{
    fprintf(stderr, "RKS RI  Not implemented yet!\n");
    abort();
    /**
    J_->zero();
    int norbs = basisset_->nbf();
    double** D = D_->to_block_matrix();
    //print_mat(D, norbs, norbs,outfile);

    double **J = block_matrix(norbs, norbs);

    double* D2 = init_array(norbs*(norbs+1)/2);
    for (int i = 0, ij = 0; i<norbs; i++) {
        for (int j = 0; j<=i; ij++, j++)
        {
            D2[ij] = (i==j?1.0:2.0)*D[i][j];
        }
    }

    double *L = init_array(naux_fin_);
    double *Gtemp = init_array(norbs*(norbs+1)/2);

    //B_ia_P_ in core
    if (df_storage_ == core||df_storage_ == double_core)
    {
        for (int i=0; i<naux_fin_; i++) {
            L[i]=C_DDOT(norbs*(norbs+1)/2,D2,1,B_ia_P_[i],1);
        }

        C_DGEMM('T','N',1,norbs*(norbs+1)/2,naux_fin_,1.0,L,1,B_ia_P_[0],norbs*(norbs+1)/2, 0.0, Gtemp, norbs*(norbs+1)/2);
        free(D2);
    }
    //B_ia_P_ on disk
    else
    {
        double *DD = init_array(norbs*(norbs+1)/2);
        for (int ij = 0; ij<norbs*(norbs+1)/2; ij++)
            DD[ij] = D2[ioff[ri_pair_mu_[ij]]+ri_pair_nu_[ij]];
        free(D2);

        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        double *in_buffer = init_array(norbs*(norbs+1)/2);
        for (int i=0; i<naux_fin_; i++) {
            psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            L[i]=C_DDOT(norbs*(norbs+1)/2,DD,1,in_buffer,1);
        }

        free(DD);
        psio_->close(PSIF_DFSCF_BJ,1);

        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        double *G2 = init_array(norbs*(norbs+1)/2);
        register double LL;
        for (int Q = 0; Q<naux_fin_; Q++)
        {
            psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            LL = L[Q];
            for (int ij = 0; ij<norbs*(norbs+1)/2; ij++)
                G2[ij]+=LL*in_buffer[ij];
        }
        free(in_buffer);
        psio_->close(PSIF_DFSCF_BJ,1);

        for (int ij = 0; ij<norbs*(norbs+1)/2; ij++)
            Gtemp[ioff[ri_pair_mu_[ij]]+ri_pair_nu_[ij]] = G2[ij];
        free(G2);
    }

    free(L);

    for (int i = 0, ij=0; i<norbs; i++) {
        for (int j = 0; j<=i; ij++,j++)
        {
            J[i][j] = Gtemp[ij];
            J[j][i] = Gtemp[ij];
        }
    }
    //fprintf(outfile, "\nJ:\n");
    //print_mat(J,norbs,norbs,outfile); fflush(outfile);
    free(Gtemp);
    free_block(D);

    for (int i=0; i<norbs; i++) {
        for (int j=0; j<=i; j++) {
            J_->set(0,i,j,J[i][j]);
            if (i!= j)
                J_->set(0,j,i,J[i][j]);
        }
    }
    //fprintf(outfile,"\n");
    //J_->print();
    free_block(J);
    **/
}
void RHF::form_K_from_RI()
{
    fprintf(stderr, "RKS RI  Not implemented yet!\n");
    abort();
    /**
    K_->zero();
    int norbs = basisset_->nbf();
    int ndocc = doccpi_[0];
    double **K = block_matrix(norbs, norbs);
    double** Cocc = block_matrix(ndocc,norbs);
    for (int i=0; i<norbs; i++) {
        for (int j=0; j<ndocc; j++)
            Cocc[j][i] = C_->get(0,i,j);
    }
    //fprintf(outfile,"\nC:\n");
    //print_mat(Cocc,ndocc,norbs,outfile);
    double** B_im_Q;
    //B_ia_P in core, B_im_Q in core
    if (df_storage_ == core||df_storage_ == double_core)
    {
        B_im_Q = block_matrix(norbs, ndocc*naux_fin_);
        double** QS = block_matrix(naux_fin_,norbs);
        double** Temp = block_matrix(ndocc,naux_fin_);

        //print_mat(B_ia_P_,naux_fin_,norbs*(norbs+1)/2,outfile);
        //fprintf(outfile,"\nYo\n");
        //print_mat(Cocc,ndocc,norbs,outfile);
        for (int m = 0; m<norbs; m++) {
            for (int Q = 0; Q<naux_fin_; Q++) {
                for (int s = 0; s<norbs; s++) {
                    QS[Q][s] = B_ia_P_[Q][((s>=m)?ioff[s]+m:ioff[m]+s)];
                }
            }
            C_DGEMM('N','T',ndocc,naux_fin_,norbs,1.0,Cocc[0],norbs,QS[0],norbs, 0.0, Temp[0], naux_fin_);
            //print_mat(Temp,ndocc,naux_fin_,outfile);
            for (int Q = 0; Q<naux_fin_; Q++) {
                for (int i = 0; i<ndocc; i++) {
                    B_im_Q[m][i+Q*ndocc] = Temp[i][Q];
                }
            }
        }
        //print_mat(B_im_Q,norbs,ndocc*naux_fin_,outfile);
        free_block(QS);
        free_block(Temp);
    }
    //B_ia_P in disk, B_im_Q in core
    else if (df_storage_ == flip_B_disk || df_storage_ == k_incore)
    {
        double *in_buffer = init_array(norbs*(norbs+1)/2);
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        B_im_Q = block_matrix(norbs, ndocc*naux_fin_);

        int mu, nu;

        for (int Q = 0; Q<naux_fin_; Q++)
        {
            psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
            {
                mu = ri_pair_mu_[ij];
                nu = ri_pair_nu_[ij];
                for (int i = 0; i<ndocc; i++)
                {
                    B_im_Q[mu][i+Q*ndocc]+=Cocc[i][nu]*in_buffer[ij];
                    if (mu != nu)
                        B_im_Q[nu][i+Q*ndocc]+=Cocc[i][mu]*in_buffer[ij];
                }
            }
        }

        psio_->close(PSIF_DFSCF_BJ,1);
        free(in_buffer);
        //B_ia_P in disk, B_im_Q in disk
    } else {
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_->open(PSIF_DFSCF_K,PSIO_OPEN_NEW);
        psio_address next_PSIF_DFSCF_K = PSIO_ZERO;

        double *in_buffer = init_array(norbs*(norbs+1)/2);
        double *out_buffer = init_array(norbs*ndocc);

        int mu, nu;
        for (int Q = 0; Q<naux_fin_; Q++)
        {
            for (int im = 0; im<ndocc*norbs; im++)
                out_buffer[im] = 0.0;

            psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
            {
                mu = ri_pair_mu_[ij];
                nu = ri_pair_nu_[ij];
                for (int i = 0; i<ndocc; i++)
                {
                    out_buffer[mu*ndocc+i]+=Cocc[i][nu]*in_buffer[ij];
                    if (mu != nu)
                        out_buffer[nu*ndocc+i]+=Cocc[i][mu]*in_buffer[ij];
                }
            }
            psio_->write(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(out_buffer[0]),sizeof(double)*norbs*ndocc,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
        }

        free(in_buffer);
        free(out_buffer);

        psio_->close(PSIF_DFSCF_BJ,1);
        psio_->close(PSIF_DFSCF_K,1);
    }

    free_block(Cocc);

    if (df_storage_ == core ||df_storage_ == double_core|| df_storage_ == flip_B_disk || df_storage_ == k_incore)
    {
        C_DGEMM('N','T',norbs,norbs,naux_fin_*ndocc,1.0,B_im_Q[0],naux_fin_*ndocc,B_im_Q[0],naux_fin_*ndocc, 0.0, K[0], norbs);
        free_block(B_im_Q);
    } else {
        psio_->open(PSIF_DFSCF_K,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_K = PSIO_ZERO;

        double *in_buffer = init_array(norbs*ndocc);

        for (int Q = 0; Q<naux_fin_; Q++)
        {
            psio_->read(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(in_buffer[0]),sizeof(double)*norbs*ndocc,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);

            for (int m = 0; m<norbs; m++)
                for (int n = 0; n<=m; n++)
                    for (int i = 0; i<ndocc; i++)
                    {
                K[m][n]+=in_buffer[m*ndocc+i]*in_buffer[n*ndocc+i];
                K[n][m] = K[m][n];
            }
        }

        free(in_buffer);
        psio_->close(PSIF_DFSCF_K,1);
    }

    //fprintf(outfile, "\nK:\n");
    //print_mat(K,norbs,norbs,outfile);

    for (int i=0; i<norbs; i++) {
        for (int j=0; j<=i; j++) {
            K_->set(0,i,j,K[i][j]);
            if (i!= j)
                K_->set(0,j,i,K[i][j]);
        }
    }
    //fprintf(outfile,"\n");
    //K_->print();
    free_block(K);**/
}
void UHF::form_G_from_RI()
{
    fprintf(stderr, "UHF RI Not implemented yet!\n");
    abort();
    /**
    * This guy needs an update

    int norbs = basisset_->nbf();
    int nalpha = nalphapi_[0];
    //int nbeta = nbetapi_[0];

    double** Da = Da_->to_block_matrix();
    double** Db = Db_->to_block_matrix();

    Ga_->zero();
    Gb_->zero();

    //fprintf(outfile,"  Arrival in form G from RI \n");fflush(outfile);

    double **J = block_matrix(norbs, norbs);

    double* D2 = init_array(norbs*(norbs+1)/2);
    for (int i = 0, ij = 0; i<norbs; i++) {
        for (int j = 0; j<=i; ij++, j++)
        {
            D2[ij] = (i==j?1.0:2.0)*Da[i][j];
            D2[ij] += (i==j?1.0:2.0)*Db[i][j];
        }
    }

    double *L = init_array(naux_fin_);
    double *Gtemp = init_array(norbs*(norbs+1)/2);

    //B_ia_P_ in core
    if (df_storage_ == core||df_storage_ == double_core)
    {

        for (int i=0; i<naux_fin_; i++) {
            L[i]=C_DDOT(norbs*(norbs+1)/2,D2,1,B_ia_P_[i],1);
        }


        C_DGEMM('T','N',1,norbs*(norbs+1)/2,naux_fin_,1.0,L,1,B_ia_P_[0],norbs*(norbs+1)/2, 0.0, Gtemp, norbs*(norbs+1)/2);
        free(D2);
    }
    //B_ia_P_ on disk
    else
    {
        double *DD = init_array(norbs*(norbs+1)/2);
        for (int ij = 0; ij<norbs*(norbs+1)/2; ij++)
        {
            DD[ij] = D2[ioff[ri_pair_mu_[ij]]+ri_pair_nu_[ij]];
            //DD[ij] += D2[ioff[ri_pair_mu_[ij]]+ri_pair_nu_[ij]];
        }
        free(D2);

        psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        double *in_buffer = init_array(norbs*(norbs+1)/2);
        for (int i=0; i<naux_fin_; i++) {
            psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            L[i]=C_DDOT(norbs*(norbs+1)/2,DD,1,in_buffer,1);
        }

        free(DD);
        psio_close(PSIF_DFSCF_BJ,1);

        psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        double *G2 = init_array(norbs*(norbs+1)/2);
        register double LL;
        for (int Q = 0; Q<naux_fin_; Q++)
        {
            psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            LL = L[Q];
            for (int ij = 0; ij<norbs*(norbs+1)/2; ij++)
            {
                G2[ij]+=LL*in_buffer[ij];
            }
        }
        free(in_buffer);
        psio_close(PSIF_DFSCF_BJ,1);

        for (int ij = 0; ij<norbs*(norbs+1)/2; ij++)
        {
            Gtemp[ioff[ri_pair_mu_[ij]]+ri_pair_nu_[ij]] = G2[ij];
        }
        free(G2);
    }

    free(L);

    for (int i = 0, ij=0; i<norbs; i++) {
        for (int j = 0; j<=i; ij++,j++)
        {
            J[i][j] = Gtemp[ij];
            J[j][i] = Gtemp[ij];
        }
    }
    //fprintf(outfile, "\nJ:\n");
    //print_mat(J,norbs,norbs,outfile); fflush(outfile);
    free(Gtemp);
    free_block(Da);
    free_block(Db);

    double** B_im_Q;
    double** Cocc = block_matrix(nalpha,norbs);

    //First Pass: Ca_
    for (int i=0; i<norbs; i++) {
        for (int j=0; j<nalpha; j++)
            Cocc[j][i] = Ca_->get(0,i,j);
    }
    //B_ia_P in core, B_im_Q in core
    if (df_storage_ == core||df_storage_ == double_core)
    {
        B_im_Q = block_matrix(norbs, nalpha*naux_fin_);
        double** QS = block_matrix(naux_fin_,norbs);
        double** Temp = block_matrix(nalpha,naux_fin_);

        //print_mat(B_ia_P_,naux_fin_,norbs*(norbs+1)/2,outfile);
        //fprintf(outfile,"\nYo\n");
        //print_mat(Cocc,nalpha,norbs,outfile);
        for (int m = 0; m<norbs; m++) {
            for (int Q = 0; Q<naux_fin_; Q++) {
                for (int s = 0; s<norbs; s++) {
                    QS[Q][s] = B_ia_P_[Q][((s>=m)?ioff[s]+m:ioff[m]+s)];
                }
            }
            //C_DGEMM('N','T',ndocc,naux_fin_,norbs,1.0,Cocc[0],norbs,QS[0],norbs, 0.0, Temp[0], naux_fin_);
            C_DGEMM('N','T',nalpha,naux_fin_,norbs,1.0,Cocc[0],norbs,QS[0],norbs, 0.0, Temp[0], naux_fin_);
            for (int Q = 0; Q<naux_fin_; Q++) {
                for (int i = 0; i<nalpha; i++) {
                    B_im_Q[m][i+Q*nalpha] = Temp[i][Q];
                }
            }
        }
        free_block(QS);
        free_block(Temp);

    }

    //B_ia_P in disk, B_im_Q in core
    else if (df_storage_ == flip_B_disk || df_storage_ == k_incore)
    {
        double *in_buffer = init_array(norbs*(norbs+1)/2);
        psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        B_im_Q = block_matrix(norbs, nalpha*naux_fin_);

        int mu, nu;

        for (int Q = 0; Q<naux_fin_; Q++)
        {
            psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
            {
                mu = ri_pair_mu_[ij];
                nu = ri_pair_nu_[ij];
                for (int i = 0; i<nalpha; i++)
                {
                    B_im_Q[mu][i+Q*nalpha]+=Cocc[i][nu]*in_buffer[ij];
                    if (mu != nu)
                        B_im_Q[nu][i+Q*nalpha]+=Cocc[i][mu]*in_buffer[ij];
                }
            }
        }

        psio_close(PSIF_DFSCF_BJ,1);
        free(in_buffer);
        //B_ia_P in disk, B_im_Q in disk
    } else {
        psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        psio_open(PSIF_DFSCF_K,PSIO_OPEN_NEW);
        psio_address next_PSIF_DFSCF_K = PSIO_ZERO;

        double *in_buffer = init_array(norbs*(norbs+1)/2);
        double *out_buffer = init_array(norbs*nalpha);

        int mu, nu;
        for (int Q = 0; Q<naux_fin_; Q++)
        {
            for (int im = 0; im<nalpha*norbs; im++)
                out_buffer[im] = 0.0;

            psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
            {
                mu = ri_pair_mu_[ij];
                nu = ri_pair_nu_[ij];
                for (int i = 0; i<nalpha; i++)
                {
                    out_buffer[mu*nalpha+i]+=Cocc[i][nu]*in_buffer[ij];
                    if (mu != nu)
                        out_buffer[nu*nalpha+i]+=Cocc[i][mu]*in_buffer[ij];
                }
            }
            psio_write(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(out_buffer[0]),sizeof(double)*norbs*nalpha,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
        }

        free(in_buffer);
        free(out_buffer);

        psio_close(PSIF_DFSCF_BJ,1);
        psio_close(PSIF_DFSCF_K,1);
    }
    free_block(Cocc);
    //fprintf(outfile,"  3-Index A Formed \n"); fflush(outfile);
    double **Ka = block_matrix(norbs, norbs);

    if (df_storage_ == core ||df_storage_ == double_core|| df_storage_ == flip_B_disk || df_storage_ == k_incore)
    {
        C_DGEMM('N','T',norbs,norbs,naux_fin_*nalpha,1.0,B_im_Q[0],naux_fin_*nalpha,B_im_Q[0],naux_fin_*nalpha, 0.0, Ka[0], norbs);
        free_block(B_im_Q);
    } else {
        psio_open(PSIF_DFSCF_K,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_K = PSIO_ZERO;

        double *in_buffer = init_array(norbs*nalpha);

        for (int Q = 0; Q<naux_fin_; Q++)
        {
            psio_read(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(in_buffer[0]),sizeof(double)*norbs*nalpha,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);

            for (int m = 0; m<norbs; m++)
                for (int n = 0; n<=m; n++)
                    for (int i = 0; i<nalpha; i++)
                    {
                Ka[m][n]+=in_buffer[m*nalpha+i]*in_buffer[n*nalpha+i];
                Ka[n][m] = Ka[m][n];
            }
        }

        free(in_buffer);
        psio_close(PSIF_DFSCF_K,0);
    }
    //fprintf(outfile,"  3-Index A Used \n"); fflush(outfile);

    //Second Pass: Cb_
    double **Kb; //Has to be defined outside the no work group
    if (nbeta_ > 0)
    {//Do work!


    Cocc = block_matrix(nbeta_,norbs);
    for (int i=0; i<norbs; i++) {
            for (int j=0; j<nbeta_; j++)
                Cocc[j][i] = Cb_->get(0,i,j);
        }

        //B_ia_P in core, B_im_Q in core
        if (df_storage_ == core||df_storage_ == double_core)
        {
            B_im_Q = block_matrix(norbs, nbeta_*naux_fin_);
            double** QS = block_matrix(naux_fin_,norbs);
            double** Temp = block_matrix(nbeta_,naux_fin_);

            //print_mat(B_ia_P_,naux_fin_,norbs*(norbs+1)/2,outfile);
            //fprintf(outfile,"\nYo\n");
            //print_mat(Cocc,nbeta_,norbs,outfile);
            for (int m = 0; m<norbs; m++) {
                for (int Q = 0; Q<naux_fin_; Q++) {
                    for (int s = 0; s<norbs; s++) {
                        QS[Q][s] = B_ia_P_[Q][((s>=m)?ioff[s]+m:ioff[m]+s)];
                    }
                }
                C_DGEMM('N','T',nbeta_,naux_fin_,norbs,1.0,Cocc[0],norbs,QS[0],norbs, 0.0, Temp[0], naux_fin_);
                //print_mat(Temp,nbeta_,naux_fin_,outfile);
                for (int Q = 0; Q<naux_fin_; Q++) {
                    for (int i = 0; i<nbeta_; i++) {
                        B_im_Q[m][i+Q*nbeta_] = Temp[i][Q];
                    }
                }
            }
            free_block(QS);
            free_block(Temp);
        }

        //B_ia_P in disk, B_im_Q in core
        else if (df_storage_ == flip_B_disk || df_storage_ == k_incore)
        {
            double *in_buffer = init_array(norbs*(norbs+1)/2);
            psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
            psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
            B_im_Q = block_matrix(norbs, nbeta_*naux_fin_);

            int mu, nu;

            for (int Q = 0; Q<naux_fin_; Q++)
            {
            psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
            {
                    mu = ri_pair_mu_[ij];
                    nu = ri_pair_nu_[ij];
                    for (int i = 0; i<nbeta_; i++)
                    {
                        B_im_Q[mu][i+Q*nbeta_]+=Cocc[i][nu]*in_buffer[ij];
                        if (mu != nu)
                            B_im_Q[nu][i+Q*nbeta_]+=Cocc[i][mu]*in_buffer[ij];
                    }
            }
            }

            psio_close(PSIF_DFSCF_BJ,1);
            free(in_buffer);
            //B_ia_P in disk, B_im_Q in disk
        } else {
            psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
            psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
            psio_open(PSIF_DFSCF_K,PSIO_OPEN_NEW);
            psio_address next_PSIF_DFSCF_K = PSIO_ZERO;

            double *in_buffer = init_array(norbs*(norbs+1)/2);
            double *out_buffer = init_array(norbs*nbeta_);

            int mu, nu;
            for (int Q = 0; Q<naux_fin_; Q++)
            {
            for (int im = 0; im<nbeta_*norbs; im++)
                    out_buffer[im] = 0.0;

            psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
            {
                    mu = ri_pair_mu_[ij];
                    nu = ri_pair_nu_[ij];
                    for (int i = 0; i<nbeta_; i++)
                    {
                        out_buffer[mu*nbeta_+i]+=Cocc[i][nu]*in_buffer[ij];
                        if (mu != nu)
                            out_buffer[nu*nbeta_+i]+=Cocc[i][mu]*in_buffer[ij];
                    }
            }
            psio_write(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(out_buffer[0]),sizeof(double)*norbs*nbeta_,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
            }

            free(in_buffer);
            free(out_buffer);

            psio_close(PSIF_DFSCF_BJ,1);
            psio_close(PSIF_DFSCF_K,1);
        }

        free_block(Cocc);
        //fprintf(outfile,"  3-Index B Formed \n"); fflush(outfile);
    Kb = block_matrix(norbs, norbs);

    if (df_storage_ == core ||df_storage_ == double_core|| df_storage_ == flip_B_disk || df_storage_ == k_incore)
        {
            C_DGEMM('N','T',norbs,norbs,naux_fin_*nbeta_,1.0,B_im_Q[0],naux_fin_*nbeta_,B_im_Q[0],naux_fin_*nbeta_, 0.0, Kb[0], norbs);
            free_block(B_im_Q);
        } else {
            psio_open(PSIF_DFSCF_K,PSIO_OPEN_OLD);
            psio_address next_PSIF_DFSCF_K = PSIO_ZERO;

            double *in_buffer = init_array(norbs*nbeta_);

            for (int Q = 0; Q<naux_fin_; Q++)
            {
            psio_read(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(in_buffer[0]),sizeof(double)*norbs*nbeta_,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);

            for (int m = 0; m<norbs; m++)
                    for (int n = 0; n<=m; n++)
                        for (int i = 0; i<nbeta_; i++)
                        {
                    Kb[m][n]+=in_buffer[m*nbeta_+i]*in_buffer[n*nbeta_+i];
                    Kb[n][m] = Kb[m][n];
                }
            }

            free(in_buffer);
            psio_close(PSIF_DFSCF_K,0);
        }
    } //Stop Work!


    for (int i=0; i<norbs; i++) {
        for (int j=0; j<=i; j++) {
            Ga_->set(0,i,j,J[i][j]-Ka[i][j]);
            if (nbeta_ > 0)
                Gb_->set(0,i,j,J[i][j]-Kb[i][j]);
            if (i!= j)
            {
                Ga_->set(0,j,i,J[i][j]-Ka[i][j]);
                if (nbeta_ > 0)
                    Gb_->set(0,j,i,J[i][j]-Kb[i][j]);
            }
        }
    }
    //Ga_->print();
    //Gb_->print();
    //fprintf(outfile,"\n");
    //G_.print();
    free_block(J);
    free_block(Ka);
    if (nbeta_ > 0)
        free_block(Kb);

    **/

}
void ROHF::form_G_from_RI()
{
    fprintf(stderr, "ROHF RI Not implemented yet!\n");
    abort();
}

}}
