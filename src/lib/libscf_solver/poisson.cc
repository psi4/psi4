/***************************************************************************
*
*       poisson.cc in psi4/src/lib/libscf_solver
*       By Rob Parrish, CCMST Georgia Tech
*       robparrish@gmail.com
*       22 November 2010
*
*       Poisson density fitting routines for SCF
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

void HF::form_B_Poisson()
{
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //                    FORM B POISSON
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    //Welcome to form_B, responsible for creation of the three-index tensor with embedded fitting metric
    //on core or disk.

    //Make sure we're in the right spot
    if (print_)
        fprintf(outfile, "\n  Computing Integrals using Poisson-Based Density Fitting\n");
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
        IntegralFactory poissonfactory(basisset_, basisset_, poissonbasis_, zero);
        B_ia_P_ = block_matrix(naux_raw_,ntri_naive_);

        //Get a TEI for each thread
        const double **buffer = new const double*[nthread];
        shared_ptr<TwoBodyAOInt> *eri = new shared_ptr<TwoBodyAOInt>[nthread];
        for (int Q = 0; Q<nthread; Q++) {
            eri[Q] = shared_ptr<TwoBodyAOInt>(rifactory.eri());
            buffer[Q] = eri[Q]->buffer();
        }
        //Get a TEI for each thread
        const double **pbuffer = new const double*[nthread];
        shared_ptr<ThreeCenterOverlapInt> *o3 = new shared_ptr<ThreeCenterOverlapInt>[nthread];
        for (int Q = 0; Q<nthread; Q++) {
            o3[Q] = shared_ptr<ThreeCenterOverlapInt>(poissonfactory.overlap_3c());
            pbuffer[Q] = o3[Q]->buffer();
        }

        int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
        int index;
        //The integrals (A|mn)
        timer_on("(A|mn)");
        #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            #ifdef _OPENMP
                rank = omp_get_thread_num();
                //fprintf(outfile,"  Thread %d doing MU = %d\n",rank,MU); fflush(outfile);
            #endif
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU) {
                numnu = basisset_->shell(NU)->nfunction();
                // == GAUSSIAN PART (FOR MULTIPOLES)
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
                // == POISSON PART (FOR SPEED)
                if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                    for (Pshell=0; Pshell < poissonbasis_->nshell(); ++Pshell) {
                        numP = poissonbasis_->shell(Pshell)->nfunction();
                        o3[rank]->compute_shell(MU, NU, Pshell);
                        for (mu=0 ; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                            for (nu=0; nu < numnu; ++nu) {
                                onu = basisset_->shell(NU)->function_index() + nu;
                                if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                    for (P=0; P < numP; ++P) {
                                        PHI = poissonbasis_->shell(Pshell)->function_index() + P;
                                        B_ia_P_[PHI + ribasis_->nbf()][ri_back_map_[omu*(omu+1)/2+onu]] = pbuffer[rank][mu*numnu*numP+nu*numP+P];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        timer_off("(A|mn)");

        //fprintf(outfile,"(A|mn)");
        //print_mat(B_ia_P_, naux_fin_,ntri_naive_ ,outfile);

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

        //fprintf(outfile,"(Q|mn)");
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

        //Get an O3 object for the AO three-index Poisson integrals
        IntegralFactory poissonfactory(basisset_, basisset_, poissonbasis_,zero);
        //Get a TEI for each thread
        const double **pbuffer = new const double*[nthread];
        shared_ptr<ThreeCenterOverlapInt> *o3 = new shared_ptr<ThreeCenterOverlapInt>[nthread];
        for (int Q = 0; Q<nthread; Q++) {
            o3[Q] = shared_ptr<ThreeCenterOverlapInt>(poissonfactory.overlap_3c());
            pbuffer[Q] = o3[Q]->buffer();
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
        //    fprintf(outfile,"  Block %d: block_starts = %d, block sizes = %d, porous starts = %d, porous sizes = %d\n",\
                block,block_starts[block],block_sizes[block],porous_starts[block],porous_sizes[block]);
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
                    //Gogo integrals! (A|mn) Gaussian
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
                    //Gogo integrals! (Amn) Poisson
                    if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                        for (Pshell=0; Pshell < poissonbasis_->nshell(); ++Pshell) {
                            numP = poissonbasis_->shell(Pshell)->nfunction();
                            o3[rank]->compute_shell(MU, NU, Pshell);
                            for (mu=0 ; mu < nummu; ++mu) {
                                omu = basisset_->shell(MU)->function_index() + mu;
                                for (nu=0; nu < numnu; ++nu) {
                                    onu = basisset_->shell(NU)->function_index() + nu;
                                    if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                        for (P=0; P < numP; ++P) {
                                            PHI = poissonbasis_->shell(Pshell)->function_index() + P;
                                            Amn[ri_back_map_[omu*(omu+1)/2+onu]-porous_starts[block]][PHI + ribasis_->nbf()] = \
                                                pbuffer[rank][mu*numnu*numP+nu*numP+P];
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
                    //Gogo integrals! (A|mn) Gaussian
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
                    //Gogo integrals! (Amn) Poisson
                    if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                        for (Pshell=0; Pshell < poissonbasis_->nshell(); ++Pshell) {
                            numP = poissonbasis_->shell(Pshell)->nfunction();
                            o3[rank]->compute_shell(MU, NU, Pshell);
                            for (mu=0 ; mu < nummu; ++mu) {
                                omu = basisset_->shell(MU)->function_index() + mu;
                                for (nu=0; nu < numnu; ++nu) {
                                    onu = basisset_->shell(NU)->function_index() + nu;
                                    if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                        for (P=0; P < numP; ++P) {
                                            PHI = poissonbasis_->shell(Pshell)->function_index() + P;
                                            Amn[PHI + ribasis_->nbf()][ri_back_map_[omu*(omu+1)/2+onu]-porous_starts[block]] =
                                                pbuffer[rank][mu*numnu*numP+nu*numP+P];
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

}}
