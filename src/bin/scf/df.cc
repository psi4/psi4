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
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include "hf.h"
#include "rhf.h"

#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/molecule.h>



using namespace std;
using namespace psi;

namespace psi { namespace scf {

void HF::form_B()
{
    fprintf(outfile, "\n  Computing Integrals using RI Basis\n");
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Must run in C1 for now.\n");
        abort();
    } 
    int norbs = basisset_->nbf(); 
    shared_ptr<BasisSet> ribasis_ =shared_ptr<BasisSet>(new BasisSet(chkpt_, "DF_BASIS"));
    ri_nbf_ = ribasis_->nbf();
    //ribasis_->print();

    fprintf(outfile, "\n  Memory Requirements:    (ab|P)    (ab|P)(PQ)^(-1/2)    Exchange Tensor    Max in Form B    Max in Form G");
    fprintf(outfile, "\n  --------------------------------------------------------------------------------------------------------");
    long memA = norbs*(norbs+1)/2*(long)ri_nbf_;
    long memB = memA;
    int ndocc = doccpi_[0];
    long memC = norbs*ndocc*(long)ri_nbf_;
    fprintf(outfile, "\n  Doubles:          %14ld %14ld      %14ld    %14ld %14ld",memA,memB,memC,memA+memB,memA+memC);
    fprintf(outfile, "\n  MiB:               %14ld %14ld      %14ld    %14ld %14ld",memA*8/1000000,memB*8/1000000,memC*8/1000000,(memA+memB)*8/1000000,(memA+memC)*8/1000000);
    fflush(outfile);

    //set df_storage_ based on available memory
    if (((long)((memA+memB)*(1.0+MEMORY_SAFETY_FACTOR)))<(memory_/sizeof(double)))
        df_storage_ = full; //Full in-core, including both (ab|P) tensors
    else if (((long)((memC)*(1.0+MEMORY_SAFETY_FACTOR)))<(memory_/sizeof(double)))
        df_storage_ = k_incore; //K only in-core
    else
        df_storage_ = disk;

    //fprintf(outfile,"\n Memory required in bytes: %f\n",(memA+memB)*(1.0+MEMORY_SAFETY_FACTOR));
    //fprintf(outfile,"\n Memory required in bytes: %d\n",(int)((memA+memB)*(1.0+MEMORY_SAFETY_FACTOR)));
    //fprintf(outfile,"\n Memory available in doubles: %d\n",memory_/sizeof(double));
    //fprintf(outfile,"\n Memory available in bytes: %d\n",memory_);
    //fprintf(outfile,"\n Memory safety factor: %f\n",MEMORY_SAFETY_FACTOR);
    df_storage_ = flip_B_disk;

    if (df_storage_ == full)
        fprintf(outfile,"\n\n  Density Fitting Algorithm proceeding In Core\n"); 
    else if (df_storage_ == flip_B_core)
        fprintf(outfile,"\n\n  Density Fitting Algorithm proceeding with K In Core, B on core, Transpose on core/disk.\n");
    else if (df_storage_ == flip_B_disk)
        fprintf(outfile,"\n\n  Density Fitting Algorithm proceeding with K In Core, B on disk, Transpose on core/disk.\n");
    else if (df_storage_ == k_incore)
        fprintf(outfile,"\n\n  Density Fitting Algorithm proceeding with K In Core, B on disk, Transpose on disk.\n");
    else if (df_storage_ == disk)
        fprintf(outfile,"\n\n  Density Fitting Algorithm proceeding on Disk\n"); 
    fflush(outfile);

    //TODO: Add cases for [(ab|P) on memory, (ab|Q) on disk->(ab|Q) on memory (ij|K) on disk] and 
    //[(ab|P) on memory, (ab|Q) on disk->(ab|Q) on memory (ij|K) on memory]

    shared_ptr<BasisSet> zero = BasisSet::zero_basis_set();

    // Create integral factory
    IntegralFactory rifactory_J(ribasis_, zero, ribasis_, zero);
    TwoBodyInt* Jint = rifactory_J.eri();
    double **J = block_matrix(ri_nbf_, ri_nbf_);
    double **J_mhalf = block_matrix(ri_nbf_, ri_nbf_);
    const double *Jbuffer = Jint->buffer();

#ifdef TIME_SCF
    timer_init();
    timer_on("Form J");
#endif

    int index = 0;

#ifdef OMP
#pragma omp parallel for 
#endif
    for (int MU=0; MU < ribasis_->nshell(); ++MU) {
        int nummu = ribasis_->shell(MU)->nfunction();

        for (int NU=0; NU < ribasis_->nshell(); ++NU) {
            int numnu = ribasis_->shell(NU)->nfunction();

            Jint->compute_shell(MU, 0, NU, 0);

            index = 0;
            for (int mu=0; mu < nummu; ++mu) {
                int omu = ribasis_->shell(MU)->function_index() + mu;

                for (int nu=0; nu < numnu; ++nu, ++index) {
                    int onu = ribasis_->shell(NU)->function_index() + nu;

                    J[omu][onu] = Jbuffer[index];
                }
            }
        }
    }
    //fprintf(outfile,"\nJ:\n");
    //print_mat(J,ri_nbf_,ri_nbf_,outfile);

    // Form J^-1/2
    // First, diagonalize J
    // the C_DSYEV call replaces the original matrix J with its eigenvectors
    double* eigval = init_array(ri_nbf_);
    int lwork = ri_nbf_ * 3;
    double* work = init_array(lwork);
    int stat = C_DSYEV('v','u',ri_nbf_,J[0],ri_nbf_,eigval, work,lwork);
    if (stat != 0) {
        fprintf(outfile, "C_DSYEV failed\n");
        exit(PSI_RETURN_FAILURE);
    }
    free(work);

    // Now J contains the eigenvectors of the original J
    // Copy J to J_copy
    double **J_copy = block_matrix(ri_nbf_, ri_nbf_);
    C_DCOPY(ri_nbf_*ri_nbf_,J[0],1,J_copy[0],1); 

    // Now form J^{-1/2} = U(T)*j^{-1/2}*U,
    // where j^{-1/2} is the diagonal matrix of the inverse square roots
    // of the eigenvalues, and U is the matrix of eigenvectors of J
    for (int i=0; i<ri_nbf_; i++) {
        if (eigval[i] < 1.0E-10)
            eigval[i] = 0.0;
        else 
            eigval[i] = 1.0 / sqrt(eigval[i]);

        // scale one set of eigenvectors by the diagonal elements j^{-1/2}
        C_DSCAL(ri_nbf_, eigval[i], J[i], 1);
    }
    free(eigval);

    // J_mhalf = J_copy(T) * J
    C_DGEMM('t','n',ri_nbf_,ri_nbf_,ri_nbf_,1.0,
        J_copy[0],ri_nbf_,J[0],ri_nbf_,0.0,J_mhalf[0],ri_nbf_);

    free_block(J);
    free_block(J_copy);

    //fprintf(outfile,"\nJmhalf:\n");
    //print_mat(J_mhalf,ri_nbf_,ri_nbf_,outfile);

#ifdef TIME_SCF
    timer_off("Form J");
    timer_on("Form ao_p_ia");
#endif
    double **ao_p_ia;
    if (df_storage_ == full||df_storage_ == flip_B_core||df_storage_ == flip_B_disk)
    {
        IntegralFactory rifactory(ribasis_, zero,basisset_, basisset_);
        TwoBodyInt* eri = rifactory.eri();
        const double *buffer = eri->buffer();
        ao_p_ia = block_matrix(ri_nbf_,basisset_->nbf()*(basisset_->nbf()+1)/2); 

        int numPshell,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
        for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
            numPshell = ribasis_->shell(Pshell)->nfunction();
            for (MU=0; MU < basisset_->nshell(); ++MU) {
                nummu = basisset_->shell(MU)->nfunction();
                for (NU=0; NU <= MU; ++NU) {
                    numnu = basisset_->shell(NU)->nfunction();
                    eri->compute_shell(Pshell, 0, MU, NU);
                    for (P=0, index=0; P < numPshell; ++P) {
                        PHI = ribasis_->shell(Pshell)->function_index() + P;
                        for (mu=0; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                            for (nu=0; nu < numnu; ++nu, ++index) {
                                onu = basisset_->shell(NU)->function_index() + nu;
                                if(omu>=onu)
                                    ao_p_ia[PHI][ioff[omu]+onu]= buffer[index];
                            } 
                        }
                    } // end loop over P in Pshell
                } // end loop over NU shell
            } // end loop over MU shell
            // now we've gone through all P, mu, nu for a given Pshell
        } // end loop over P shells; done with forming MO basis (P|ia)'s*/
        fprintf(outfile,"\n  Through ao_p_ia in core\n"); fflush(outfile);
        //fprintf(outfile,"\nao_p_ia:\n");
        //print_mat(ao_p_ia, ri_nbf_,norbs*(norbs+1)/2 ,outfile);
    }
    else
    {
        psio_->open(PSIF_DFSCF_B,PSIO_OPEN_NEW);

        IntegralFactory rifactory(ribasis_, zero, basisset_, basisset_);
        TwoBodyInt* eri = rifactory.eri();
        const double *buffer = eri->buffer();
        ao_p_ia = block_matrix(ri_nbf_,1); 
        double **storage;
        ri_pair_nu_ = init_int_array(basisset_->nbf()*(basisset_->nbf()+1)/2);
        ri_pair_mu_ = init_int_array(basisset_->nbf()*(basisset_->nbf()+1)/2);
        int row, npairs;

        psio_address next_PSIF_DFSCF_B = PSIO_ZERO;

        int pair_index = 0;

        int numPshell,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU) {
                numnu = basisset_->shell(NU)->nfunction();

                if (MU != NU)
                    npairs = nummu*numnu;
                else
                    npairs = nummu*(nummu+1)/2;
                //fprintf(outfile,"\n  Starting Computing Quartet (%d %d| P)",MU,NU); fflush(outfile);
                if (MU != 0 || NU !=0)
                    free_block(storage);
                storage = block_matrix(ri_nbf_,npairs);
                //fprintf(outfile,"\n  Memory Allocated for Quartet (%d %d| P)",MU,NU); fflush(outfile);

                for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
                    numPshell = ribasis_->shell(Pshell)->nfunction();

                    eri->compute_shell(Pshell, 0, MU, NU);

                    for (P=0, index=0; P < numPshell; ++P) {
                        PHI = ribasis_->shell(Pshell)->function_index() + P;
                        row = 0;
                        for (mu=0; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                            for (nu=0; nu < numnu; ++nu, ++index) {
                                onu = basisset_->shell(NU)->function_index() + nu;
                                if(omu>=onu)
                                {
                                    storage[PHI][row++]= buffer[index];
                                    if (Pshell == 0)
                                    {
                                        ri_pair_nu_[pair_index] = onu;
                                        ri_pair_mu_[pair_index++] = omu;
                                    }
                                } 
                            } 
                        }
                    } // end loop over P in Pshell

                }
                //fprintf(outfile,"\n  Finished Computing Quartet (%d %d| P)",MU,NU); fflush(outfile);
                //print_mat(storage,ri_nbf_,npairs,outfile);

                row = 0;
                for (mu=0; mu < nummu; ++mu) {
                    omu = basisset_->shell(MU)->function_index() + mu;
                    for (nu=0; nu < numnu; ++nu) {
                        onu = basisset_->shell(NU)->function_index() + nu;
                        if(omu>=onu)
                        {
                            //fprintf(outfile,"\n  WE here!"); fflush(outfile);
                            for (P = 0; P<ri_nbf_; P++)
                                ao_p_ia[P][0] = storage[P][row];
                            //fprintf(outfile,"\n  Finished Transposing Quartet (%d %d| P)\n",MU,NU); fflush(outfile);
                            psio_->write(PSIF_DFSCF_B,"B Three-Index Integrals",(char *) &(ao_p_ia[0][0]),sizeof(double)*ri_nbf_,next_PSIF_DFSCF_B,&next_PSIF_DFSCF_B);
                            row++;
                        }
                    }
                }
                //fprintf(outfile,"\n  Finished Writing Quartet (%d %d| P)\n",MU,NU); fflush(outfile);
            }
        }

        free(storage);
        free(ao_p_ia);
        fprintf(outfile,"\n  Through B on disk."); fflush(outfile);
        psio_->close(PSIF_DFSCF_B,1);
    }

#ifdef TIME_SCF
    timer_off("Form ao_p_ia");
    timer_on("Form B_ia^P");
#endif
    //ao_p_ia in core, B_ia_P in core
    if (df_storage_ == full)
    {
        // ao_p_ia has integrals
        // B_ia^P = Sum_Q (i a | Q) (J^-1/2)_QP
        B_ia_P_ = block_matrix(ri_nbf_,norbs*(norbs+1)/2);

        C_DGEMM('N','N',ri_nbf_,norbs*(norbs+1)/2,ri_nbf_,
            1.0, J_mhalf[0], ri_nbf_, ao_p_ia[0], norbs*(norbs+1)/2,
            0.0, B_ia_P_[0], norbs*(norbs+1)/2);
        //fprintf(outfile,"\nB_p_ia:\n");
        //print_mat(B_ia_P_, ri_nbf_,norbs*(norbs+1)/2 ,outfile);
        free_block(ao_p_ia);
        free_block(J_mhalf);
    }
    //ao_p_ia in core, B_ia_P on disk => B_ia_P on core
    else if (df_storage_ == flip_B_core )
    {
        //flip B by blocks to disk and send to B_ia_P
        double **buffer = block_matrix(ri_nbf_,1);
        double **temp = block_matrix(ri_nbf_,1);
        psio_->open(PSIF_DFSCF_BJI,PSIO_OPEN_NEW);
        psio_address next_PSIF_DFSCF_BJI = PSIO_ZERO;
        for (int ij = 0 ; ij < norbs*(norbs+1)/2; ij++)
        {
            for (int Q = 0; Q<ri_nbf_; Q++)
                temp[Q][0] = ao_p_ia[Q][ij];

            C_DGEMM('N','N',ri_nbf_,1,ri_nbf_,1.0, J_mhalf[0], ri_nbf_, &(temp[0][0]), 1,0.0, &(buffer[0][0]), 1);

            psio_->write(PSIF_DFSCF_BJI,"BJ Three-Index Integrals",(char *) &(buffer[0][0]),sizeof(double)*ri_nbf_,next_PSIF_DFSCF_BJI,&next_PSIF_DFSCF_BJI);
        }
        psio_->close(PSIF_DFSCF_BJI,1);
        free(ao_p_ia);
        free(temp);

        B_ia_P_ = block_matrix(ri_nbf_,norbs*(norbs+1)/2);
        psio_->open(PSIF_DFSCF_BJI,PSIO_OPEN_OLD);
        next_PSIF_DFSCF_BJI = PSIO_ZERO;
        for (int ij = 0 ; ij < norbs*(norbs+1)/2; ij++)
        {
            psio_->read(PSIF_DFSCF_BJI,"BJ Three-Index Integrals",(char *) &(buffer[0][0]),sizeof(double)*ri_nbf_,next_PSIF_DFSCF_BJI,&next_PSIF_DFSCF_BJI);

            for (int Q = 0; Q<ri_nbf_; Q++)
                B_ia_P_[Q][ij] = buffer[Q][0];
        }
        psio_->close(PSIF_DFSCF_BJI,0);
        free(buffer);
    }
    //ao_p_ia in core, B_ia_P on disk => B_ia_P on disk
    else if (df_storage_ == flip_B_disk )
    {
    //flip B by blocks to disk and leave it there
        double **buffer = block_matrix(ri_nbf_,1);
        double **temp = block_matrix(ri_nbf_,1);
        psio_->open(PSIF_DFSCF_BJI,PSIO_OPEN_NEW);
        psio_address next_PSIF_DFSCF_BJI = PSIO_ZERO;
        for (int ij = 0 ; ij < norbs*(norbs+1)/2; ij++)
        {
            for (int Q = 0; Q<ri_nbf_; Q++)
                temp[Q][0] = ao_p_ia[Q][ij];

            C_DGEMM('N','N',ri_nbf_,1,ri_nbf_,1.0, J_mhalf[0], ri_nbf_, &(temp[0][0]), 1,0.0, &(buffer[0][0]), 1);

            psio_->write(PSIF_DFSCF_BJI,"BJ Three-Index Integrals",(char *) &(buffer[0][0]),sizeof(double)*ri_nbf_,next_PSIF_DFSCF_BJI,&next_PSIF_DFSCF_BJI);
        }
        psio_->close(PSIF_DFSCF_BJI,1);
        free(ao_p_ia);
        free(temp);

        B_ia_P_ = block_matrix(ri_nbf_,norbs*(norbs+1)/2);
        psio_->open(PSIF_DFSCF_BJI,PSIO_OPEN_OLD);
        next_PSIF_DFSCF_BJI = PSIO_ZERO;
        for (int ij = 0 ; ij < norbs*(norbs+1)/2; ij++)
        {
            psio_->read(PSIF_DFSCF_BJI,"BJ Three-Index Integrals",(char *) &(buffer[0][0]),sizeof(double)*ri_nbf_,next_PSIF_DFSCF_BJI,&next_PSIF_DFSCF_BJI);

            for (int Q = 0; Q<ri_nbf_; Q++)
                B_ia_P_[Q][ij] = buffer[Q][0];
        }
        psio_->close(PSIF_DFSCF_BJI,0);
        free(buffer);

    //RESTRIPE
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        for (int Q = 0; Q<ri_nbf_; Q++)
        {
            psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(B_ia_P_[Q][0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
        }
        psio_->close(PSIF_DFSCF_BJ,1);
        free(B_ia_P_);

        ri_pair_mu_ = init_int_array(norbs*(norbs+1)/2);
        ri_pair_nu_ = init_int_array(norbs*(norbs+1)/2);
        for (int i = 0, ij = 0; i< norbs; i++) {
            for (int j = 0; j<=i; j++, ij++)
            {
                ri_pair_mu_[ij] = i;
                ri_pair_nu_[ij] = j;
            }
        }
    }
    //ao_p_ia in disk, B_ia_P on disk => B_ia_P on disk
    else 
    {
        psio_->open(PSIF_DFSCF_B,PSIO_OPEN_OLD);
        psio_->open(PSIF_DFSCF_BJI,PSIO_OPEN_NEW);
        psio_address next_PSIF_DFSCF_B = PSIO_ZERO;
        psio_address next_PSIF_DFSCF_BJI = PSIO_ZERO;

        double **in_buffer = block_matrix(ri_nbf_,1);
        double **out_buffer = block_matrix(ri_nbf_,1);

        for (int ij = 0; ij<norbs*(norbs+1)/2; ij++)
        {
            psio_->read(PSIF_DFSCF_B,"B Three-Index Integrals",(char *) &(in_buffer[0][0]),sizeof(double)*ri_nbf_,next_PSIF_DFSCF_B,&next_PSIF_DFSCF_B);

            C_DGEMM('N','N',ri_nbf_,1,ri_nbf_,1.0, J_mhalf[0], ri_nbf_, in_buffer[0], 1,0.0, &(out_buffer[0][0]), 1);

            psio_->write(PSIF_DFSCF_BJI,"BJ Three-Index Integrals",(char *) &(out_buffer[0][0]),sizeof(double)*ri_nbf_,next_PSIF_DFSCF_BJI,&next_PSIF_DFSCF_BJI);
        }
        free(out_buffer);

        psio_->close(PSIF_DFSCF_BJI,1);
        psio_->close(PSIF_DFSCF_B,0);
        fprintf(outfile,"\n  Through BJ on disk."); fflush(outfile);

    //RESTRIPE
        psio_->open(PSIF_DFSCF_BJI,PSIO_OPEN_OLD);
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
        next_PSIF_DFSCF_BJI = PSIO_ZERO;
        psio_address next_PSIF_DFSCF_BJ;

        int max_cols = (memory_/sizeof(double))/((1.0+MEMORY_SAFETY_FACTOR)*ri_nbf_);
        if (max_cols > norbs*(norbs+1)/2)
            max_cols = norbs*(norbs+1)/2;

        //max_cols = 100;
        double **buffer = block_matrix(ri_nbf_,max_cols);

        int buf_ind = 0;
        ULI global_offset = 0;
        for (int ij = 0; ij < norbs*(norbs+1)/2; ij++)
        {
            psio_->read(PSIF_DFSCF_BJI,"BJ Three-Index Integrals",(char *) &(in_buffer[0][0]),sizeof(double)*ri_nbf_,next_PSIF_DFSCF_BJI,&next_PSIF_DFSCF_BJI);
            //fprintf(outfile,"\n  Read in pair %d",ij); fflush(outfile);
            for (int Q = 0; Q<ri_nbf_; Q++)
            {
                buffer[Q][buf_ind] = in_buffer[Q][0];
            }
            buf_ind++;
            //fprintf(outfile,"\n  Moved Pair to position %d in buffer",buf_ind); fflush(outfile);
            if (buf_ind == max_cols || ij == norbs*(norbs+1)/2-1)
            {
                fprintf(outfile,"\n  Writing %d pairs to disk",buf_ind); fflush(outfile);
                for (int Q = 0; Q<ri_nbf_; Q++)
                {
                    //fprintf(outfile,"\n  Working on Q %d",Q); fflush(outfile);
                    next_PSIF_DFSCF_BJ = psio_get_address(PSIO_ZERO,(ULI)(Q*norbs*(ULI)(norbs+1)/2*sizeof(double)+global_offset*sizeof(double)));
                    //fprintf(outfile,"\n  Address Acquired"); fflush(outfile);
                    //for (int K = 0; K<buf_ind; K++)
                        //out_buffer[K][0] = buffer[Q][K];
                    //fprintf(outfile,"\n  out_buffer transposed"); fflush(outfile);
                    //errcode = psio_write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(out_buffer[0][0]),sizeof(double)*buf_ind,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
                    psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(buffer[Q][0]),sizeof(double)*buf_ind,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
                    //fprintf(outfile,"\n  Entry Written"); fflush(outfile);
                }
                global_offset+=buf_ind;
                buf_ind=0;
            } 
        }

        free(in_buffer);
        free(buffer);
        psio_->close(PSIF_DFSCF_BJI,0);
        psio_->close(PSIF_DFSCF_BJ,1);

        fprintf(outfile,"\n  BJ Restriped on disk.\n"); fflush(outfile);
    }
#ifdef TIME_SCF
    timer_off("Form B_ia^P");
#endif 
}

void RHF::form_G_from_RI()
{
    int norbs = basisset_->nbf(); 
    double** D = D_->to_block_matrix();

    G_->zero();

    //print_mat(D, norbs, norbs,outfile);

#ifdef TIME_SCF
    timer_on("Form Coulomb");
#endif

    double **J = block_matrix(norbs, norbs);

    double* D2 = init_array(norbs*(norbs+1)/2);
    for (int i = 0, ij = 0; i<norbs; i++) {
        for (int j = 0; j<=i; ij++, j++)
        {
            D2[ij] = (i==j?1.0:2.0)*D[i][j];
        }
    }

    double *L = init_array(ri_nbf_);
    double *Gtemp = init_array(norbs*(norbs+1)/2);

    //B_ia_P_ in core
    if (df_storage_ == full||df_storage_ == flip_B_core)
    {
        for (int i=0; i<ri_nbf_; i++) {
            L[i]=C_DDOT(norbs*(norbs+1)/2,D2,1,B_ia_P_[i],1);
        }

        C_DGEMM('T','N',1,norbs*(norbs+1)/2,ri_nbf_,1.0,L,1,B_ia_P_[0],norbs*(norbs+1)/2, 0.0, Gtemp, norbs*(norbs+1)/2);    
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
        for (int i=0; i<ri_nbf_; i++) {
            psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            L[i]=C_DDOT(norbs*(norbs+1)/2,DD,1,in_buffer,1);
        }

        free(DD);
        psio_->close(PSIF_DFSCF_BJ,1);

        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        next_PSIF_DFSCF_BJ = PSIO_ZERO;
        double *G2 = init_array(norbs*(norbs+1)/2);
        register double LL;
        for (int Q = 0; Q<ri_nbf_; Q++)
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
            J[i][j] = 2.0*Gtemp[ij];
            J[j][i] = 2.0*Gtemp[ij];
        }
    }
    //fprintf(outfile, "\nJ:\n");
    //print_mat(J,norbs,norbs,outfile); fflush(outfile);
    free(Gtemp);
    free_block(D);

#ifdef TIME_SCF
    timer_off("Form Coulomb");
    timer_on("Form Exchange");
#endif
    int ndocc = doccpi_[0];
    double** Cocc = block_matrix(ndocc,norbs);
    for (int i=0; i<norbs; i++) {
        for (int j=0; j<ndocc; j++)
            Cocc[j][i] = C_->get(0,i,j);
    }
    double** B_im_Q;
    //B_ia_P in core, B_im_Q in core
    if (df_storage_ == full||df_storage_ == flip_B_core)
    {
        B_im_Q = block_matrix(norbs, ndocc*ri_nbf_);
        double** QS = block_matrix(ri_nbf_,norbs);
        double** Temp = block_matrix(ndocc,ri_nbf_);

        //print_mat(B_ia_P_,ri_nbf_,norbs*(norbs+1)/2,outfile);
        //fprintf(outfile,"\nYo\n");
        //print_mat(Cocc,ndocc,norbs,outfile);
        for (int m = 0; m<norbs; m++) {
            for (int Q = 0; Q<ri_nbf_; Q++) {
                for (int s = 0; s<norbs; s++) {
                    QS[Q][s] = B_ia_P_[Q][((s>=m)?ioff[s]+m:ioff[m]+s)];
                }
            }
            C_DGEMM('N','T',ndocc,ri_nbf_,norbs,1.0,Cocc[0],norbs,QS[0],norbs, 0.0, Temp[0], ri_nbf_);
            //print_mat(Temp,ndocc,ri_nbf_,outfile);
            for (int Q = 0; Q<ri_nbf_; Q++) {
                for (int i = 0; i<ndocc; i++) {
                    B_im_Q[m][i+Q*ndocc] = Temp[i][Q];
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
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        B_im_Q = block_matrix(norbs, ndocc*ri_nbf_);

        int mu, nu;

        for (int Q = 0; Q<ri_nbf_; Q++)
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
        for (int Q = 0; Q<ri_nbf_; Q++)
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

    double **K = block_matrix(norbs, norbs);

    if (df_storage_ == full || df_storage_ == flip_B_core || df_storage_ == flip_B_disk || df_storage_ == k_incore)
    {
        C_DGEMM('N','T',norbs,norbs,ri_nbf_*ndocc,1.0,B_im_Q[0],ri_nbf_*ndocc,B_im_Q[0],ri_nbf_*ndocc, 0.0, K[0], norbs);
        free_block(B_im_Q);
    } else {
        psio_->open(PSIF_DFSCF_K,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_K = PSIO_ZERO;

        double *in_buffer = init_array(norbs*ndocc);

        for (int Q = 0; Q<ri_nbf_; Q++)
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
#ifdef TIME_SCF
    timer_off("Form Exchange");
#endif

    for (int i=0; i<norbs; i++) {
        for (int j=0; j<=i; j++) {
            G_->add(0,i,j,J[i][j]-K[i][j]);
            if (i!= j)
                G_->add(0,j,i,J[i][j]-K[i][j]);
        }
    }
    //fprintf(outfile,"\n");
    //G_.print();
    free_block(J);
    free_block(K);
}
}}
