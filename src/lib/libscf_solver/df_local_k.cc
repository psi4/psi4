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
#include "uhf.h"
#include "rohf.h"

#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/molecule.h>

using namespace std;
using namespace psi;

namespace psi { namespace scf {

void HF::form_A()
{
    fprintf(outfile, "\n  Computing Integrals using Density Fitting\n");
    //TODO: Add support for molecular symmetry
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Must run in C1 for now.\n"); fflush(outfile);
        abort();
    } 
    int norbs = basisset_->nbf(); 
    shared_ptr<BasisSet> ribasis_ =shared_ptr<BasisSet>(new BasisSet(chkpt_, "DF_BASIS_SCF"));
    ri_nbf_ = ribasis_->nbf();
    
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
        double max_global_val;

        IntegralFactory schwarzfactory(basisset_,basisset_,basisset_,basisset_);
        shared_ptr<TwoBodyInt> eri = shared_ptr<TwoBodyInt>(schwarzfactory.eri());
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

    timer_off("Schwarz Sieve");
    //Size of the three-index tensor
    unsigned long memA = ntri_*(long)ri_nbf_;

    string storage_type;
    storage_type = options_.get_str("RI_STORAGE");

    if (storage_type == "IN_CORE")
        df_storage_ = full;
    else if (storage_type == "DISK")
        df_storage_ = disk;
    else if (storage_type == "DEFAULT")
    {
    	//set df_storage_ semi-heuristically based on available memory
    	if (((long)((memA)*(1.0+MEMORY_SAFETY_FACTOR)))<(memory_/sizeof(double)))
            df_storage_ = full; //K only in-core
    	else
            df_storage_ = disk; //Disk
    } else {
        throw std::domain_error("For Form A, IN_CORE and DISK are the only valid storage options");
    }	

    if (df_storage_ == full)
        fprintf(outfile,"\n  Density Fitting Algorithm proceeding In Core.\n"); 
    else if (df_storage_ == disk)
        fprintf(outfile,"\n  Density Fitting Algorithm proceeding on Disk.\n"); 
    fflush(outfile);

    //It takes a lot of work to get a null basis with Psi4 
    shared_ptr<BasisSet> zero = BasisSet::zero_basis_set();
    
    // Create integral factory for J (Fitting Matrix in form_B)
    IntegralFactory rifactory_J(ribasis_, zero, ribasis_, zero);
    shared_ptr<TwoBodyInt> Jint = shared_ptr<TwoBodyInt>(rifactory_J.eri());

    // Integral buffer
    const double *Jbuffer = Jint->buffer();

    // J Matrix (will later hold Jinv, 
    //LAPACK is in place
    Jinv_ = block_matrix(ri_nbf_, ri_nbf_);
    
    timer_on("Form J Matrix;");
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

                    Jinv_[omu][onu] = Jbuffer[index];
                    Jinv_[onu][omu] = Jbuffer[index];
                }
            }
        }
    }
    //fprintf(outfile,"\nJ:\n"); fflush(outfile);
    //print_mat(Jinv_,ri_nbf_,ri_nbf_,outfile);
    timer_off("Form J Matrix;");
    timer_on("Form J^-1");
    
    //Pivot indices
    int* ipiv = init_int_array(ri_nbf_);

    //LU Factorization (in place)
    int LUerror = C_DGETRF(ri_nbf_,ri_nbf_,Jinv_[0],ri_nbf_,ipiv);
    if (LUerror != 0)
        throw std::domain_error("J Matrix LU Decomposition Failed!");

    //Work array
    double *work = init_array(3*ri_nbf_);   
   
    //Perform the inversion 
    int Ierror = C_DGETRI(ri_nbf_,Jinv_[0],ri_nbf_,ipiv,work,3*ri_nbf_);
    if (Ierror != 0)
        throw std::domain_error("J Matrix Inversion Failed!");

    free(work);
    free(ipiv); 

    timer_off("Form J^-1");

    //fprintf(outfile,"\nJinv_:\n"); fflush(outfile);
    //print_mat(Jinv_,ri_nbf_,ri_nbf_,outfile);
    
    timer_on("Overall (A|mn)");
    
    //Use ri_pair_mu_ and ri_pair_nu_ to keep track of things
    //Across schwarz sieve and unfortunate shell indexing
    ri_pair_nu_ = init_int_array(ntri_naive_);
    ri_pair_mu_ = init_int_array(ntri_naive_);
  
    //double three_index_cutoff = options_.get_double("THREE_INDEX_CUTOFF");
 
    if (df_storage_ == full)
    {
    	IntegralFactory rifactory(basisset_, basisset_, ribasis_, zero);
        shared_ptr<TwoBodyInt> eri = shared_ptr<TwoBodyInt>(rifactory.eri());
        const double *buffer = eri->buffer();
        A_ia_P_ = block_matrix(ri_nbf_,ntri_naive_); 

        int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
        int start_index, delta_index, l_index;
        start_index = 0;
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU) {
                numnu = basisset_->shell(NU)->nfunction();
                if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                    delta_index = 0;
                    for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
                        numP = ribasis_->shell(Pshell)->nfunction();
                        timer_on("(A|mn) Integrals");
                        eri->compute_shell(MU, NU, Pshell, 0);
                        timer_off("(A|mn) Integrals");
                        l_index = start_index;
                        for (mu=0 ; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                            for (nu=0; nu < numnu; ++nu) {
                                onu = basisset_->shell(NU)->function_index() + nu;
                                if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                    for (P=0; P < numP; ++P) {
                                        PHI = ribasis_->shell(Pshell)->function_index() + P;
                                        A_ia_P_[PHI][l_index]= buffer[mu*numnu*numP+nu*numP+P];
                                    }
                                    if (Pshell == 0) {
                                        delta_index++;
                                        ri_pair_mu_[l_index] = omu;
                                        ri_pair_nu_[l_index] = onu;
                                    }
                                    l_index++;
                                } 
                            }
                        }
                    }
                    start_index+=delta_index;
                }
            } 
        }
        //print_mat(A_ia_P, ri_nbf_,ntri_naive_ ,outfile);
        //for (int left = 0; left<ntri_; left++)
        //    fprintf(outfile,"  %d pair: (%d, %d)\n",left,ri_pair_mu_[left],ri_pair_nu_[left]);

        //fflush(outfile);
    } 
    else if (df_storage_ == disk)
    {
        psio_->open(PSIF_DFSCF_BJI,PSIO_OPEN_NEW);
        psio_address next_PSIF_DFSCF_BJI = PSIO_ZERO;
        
        IntegralFactory rifactory(basisset_, basisset_, ribasis_,zero);
        shared_ptr<TwoBodyInt> eri = shared_ptr<TwoBodyInt>(rifactory.eri());
        const double *buffer = eri->buffer();
        int maxPfun = 0;
        for (int m = 0; m<ribasis_->nshell(); m++)
            if (maxPfun<ribasis_->shell(m)->nfunction())
                maxPfun=ribasis_->shell(m)->nfunction();
        double *Temp1 = init_array(ntri_*maxPfun);

        int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
        int start_index, delta_index, l_index,r_index, s_index;
        start_index = 0; l_index = 0, r_index = 0;
        for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
            numP = ribasis_->shell(Pshell)->nfunction();
            l_index = 0;
            for (MU=0; MU < basisset_->nshell(); ++MU) {
                nummu = basisset_->shell(MU)->nfunction();
                for (NU=0; NU <= MU; ++NU) {
                    numnu = basisset_->shell(NU)->nfunction();
                    //fprintf(outfile, "  MU = %d, NU = %d, Sig = %d\n",MU,NU,schwarz_shell_pairs[MU*(MU+1)/2+NU]); fflush(outfile);
                    if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                        timer_on("(A|mn) Integrals");
                        eri->compute_shell(MU, NU, Pshell, 0);
                        timer_off("(A|mn) Integrals");
                        for (P=0; P < numP; ++P) {
                            PHI = ribasis_->shell(Pshell)->function_index() + P;
                            for (mu=0 ; mu < nummu; ++mu) {
                                omu = basisset_->shell(MU)->function_index() + mu;
                                for (nu=0; nu < numnu; ++nu) {
                                    onu = basisset_->shell(NU)->function_index() + nu;
                                    if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                        Temp1[l_index]= buffer[mu*numnu*numP+nu*numP+P];
                                    }
                                    if (PHI == 0) {
                                        ri_pair_mu_[l_index] = omu;
                                        ri_pair_nu_[l_index] = onu;
                                    }
                                    l_index++;
                                } 
                            }
                        }
                    }
                }
            }
            timer_on("(A|mn) disk");
            psio_->write(PSIF_DFSCF_BJI,"(A|mn) Three-Index Integrals",(char *) Temp1,sizeof(double)*ntri_*numP,next_PSIF_DFSCF_BJI,&next_PSIF_DFSCF_BJI);
            timer_off("(A|mn) disk");
        } 
        //for (int left = 0; left<ntri_; left++)
         //   fprintf(outfile,"  %d pair: (%d, %d)\n",left,ri_pair_mu_[left],ri_pair_nu_[left]);

        free(Temp1);
        //fprintf(outfile,"\n  Through B on disk."); fflush(outfile);
        psio_->close(PSIF_DFSCF_BJI,1);
        /**timer_on("(B|mn) restripe");
        
    	//RESTRIPE
        psio_->open(PSIF_DFSCF_BJI,PSIO_OPEN_OLD);
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
	next_PSIF_DFSCF_BJI = PSIO_ZERO;
	psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
	
	double *Temp = init_array(ntri_);
	for (int Q = 0; Q < ri_nbf_; Q++)
            psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(Temp[0]),sizeof(double)*ntri_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
	free(Temp);	

        int max_cols = (memory_/sizeof(double))/((1.0+MEMORY_SAFETY_FACTOR)*ri_nbf_);
        if (max_cols > ntri_)
            max_cols = ntri_;
	double *in_buffer = init_array(ri_nbf_);
        //max_cols = 100;
        double **buffer2 = block_matrix(ri_nbf_,max_cols);

        int buf_ind = 0;
        ULI global_offset = 0;
        for (int ij = 0; ij < ntri_; ij++)
        {
            psio_->read(PSIF_DFSCF_BJI,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*ri_nbf_,next_PSIF_DFSCF_BJI,&next_PSIF_DFSCF_BJI);
            //fprintf(outfile,"\n  Read in pair %d",ij); fflush(outfile);
            for (int Q = 0; Q<ri_nbf_; Q++)
            {
                buffer2[Q][buf_ind] = in_buffer[Q];
            }
            buf_ind++;
            //fprintf(outfile,"\n  Moved Pair to position %d in buffer",buf_ind); fflush(outfile);
            if (buf_ind == max_cols || ij == ntri_-1)
            {
                //fprintf(outfile,"\n  Writing %d pairs to disk",buf_ind); fflush(outfile);
                for (int Q = 0; Q<ri_nbf_; Q++)
                {
                    //fprintf(outfile,"\n  Working on Q %d",Q); fflush(outfile);
                    next_PSIF_DFSCF_BJ = psio_get_address(PSIO_ZERO,(ULI)(Q*(ULI)ntri_*sizeof(double)+global_offset*sizeof(double)));
                    //fprintf(outfile,"\n  Address Acquired"); fflush(outfile);
                    psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(buffer2[Q][0]),sizeof(double)*buf_ind,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
                    //fprintf(outfile,"\n  Entry Written"); fflush(outfile);
                }
                global_offset+=buf_ind;
                buf_ind=0;
            } 
        }

        free(in_buffer);
        free_block(buffer2);
        psio_->close(PSIF_DFSCF_BJI,0);
        psio_->close(PSIF_DFSCF_BJ,1);
        timer_off("(B|mn) restripe");**/

        //fprintf(outfile,"\n  BJ Restriped on disk.\n"); fflush(outfile);
    }
    timer_off("Overall (A|mn)");

    if (schwarz_) {
        fprintf(outfile,"\n  Function Pair Schwarz Sieve, Cutoff = %14.10E:\n",schwarz_);
        fprintf(outfile,"  %d out of %d basis function pairs removed, %8.5f%% attenuation.\n",norbs*(norbs+1)/2-sig_fun_pairs,norbs*(norbs+1)/2,100.0*(norbs*(norbs+1)/2-sig_fun_pairs)/(1.0*norbs*(norbs+1)/2));
        int pairs = basisset_->nshell()*(basisset_->nshell()+1)/2;
        fprintf(outfile,"  %d out of %d basis shell pairs removed, %8.5f%% attenuation.\n",pairs-sig_shell_pairs,pairs,100.0*(pairs-sig_shell_pairs)/(1.0*pairs));
    }
    /**
    if (three_index_cutoff) {
        int attenuation = ntri_naive_-ntri_;
        fprintf(outfile,"  Direct Three-Index Tensor Sieve, Cutoff = %14.10E:\n",three_index_cutoff);
        fprintf(outfile,"  %d of %d (remaining) basis function pairs removed, %8.5f%% attenuation.\n",attenuation,ntri_naive_, 100.0*attenuation/(1.0*ntri_naive_));
    }
    if (schwarz_>0.0 || three_index_cutoff>0.0)
        fprintf(outfile,"  After sieving, %d out of %d basis function pairs remain, %8.5f%% attenuation.\n\n",ntri_,norbs*(norbs+1)/2,100.0*(1.0-ntri_/(1.0*norbs*(norbs+1)/2)));**/
}
/**
void HF::write_B()
{
    psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
    psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
            
    for (int Q = 0; Q<ri_nbf_; Q++) {
        psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(B_ia_P_[Q][0]),sizeof(double)*ntri_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
                
     }
     psio_->close(PSIF_DFSCF_BJ,1);
}**/
void HF::free_A()
{
    if (df_storage_ == full)
        free(A_ia_P_);
    free(ri_pair_mu_);
    free(ri_pair_nu_);
    free_block(Jinv_);
}
void RHF::form_G_from_RI_local_K()
{
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
    double* DD;
    
    if (J_is_required_) {
        DD = init_array(ntri_); 
        for (int ij = 0; ij<ntri_; ij++) {
            DD[ij] = D_->get(0,ri_pair_mu_[ij],ri_pair_nu_[ij]); 
            if (ri_pair_mu_[ij] != ri_pair_nu_[ij])
                DD[ij] *= 2.0;
                //only A irrep at the moment!!
        }
    }
    //Get the C matrix (exchange messes things up)
    int ndocc = doccpi_[0];
    double** Cocc;
    if (K_is_required_) {
        Cocc = block_matrix(ndocc,norbs);
        for (int i=0; i<norbs; i++) {
            for (int j=0; j<ndocc; j++)
                Cocc[j][i] = C_->get(0,i,j);
            //only A irrep at the moment!!
        }
    }
    if (df_storage_ == full) {
        if (J_is_required_) {
            //As the great Ryann Kopacka would say, Play Your Game!
            //Oh, and see R. Polly et. al., J. Chem. Phys. 102(21-22), pp. 2311-2321, DOI 10.1080/0026897042000274801
            double* c = init_array(ri_nbf_);
            double* d = init_array(ri_nbf_);
            double *J = init_array(ntri_);
            //DGEMV -> L:
            //L_Q = (Q|ls)*D_{ls}
            timer_on("J DDOT");
            C_DGEMV('N',ri_nbf_,ntri_,1.0,A_ia_P_[0],ntri_naive_,DD,1,0.0,c,1);
            timer_off("J DDOT");

            //d_A = J^{-1}*c = J^{-1}_{AB}c_B
            timer_on("J Fitting");
            C_DGEMV('N',ri_nbf_,ri_nbf_,1.0,Jinv_[0],ri_nbf_,c,1,0.0,d,1);
            timer_off("J Fitting");

            //DGEMV -> J:
            //J_{mn} = L_Q(Q|mn)
            timer_on("J DAXPY");
            C_DGEMV('T',ri_nbf_,ntri_,1.0,A_ia_P_[0],ntri_naive_,d,1,0.0,J,1);
            timer_off("J DAXPY");
            //Put everything in J_
            for (int ij = 0; ij < ntri_; ij++) {
                J_->set(0,ri_pair_mu_[ij],ri_pair_nu_[ij],J[ij]);
                J_->set(0,ri_pair_nu_[ij],ri_pair_mu_[ij],J[ij]);
            }
            free(J);
            free(d);
            free(c);
        }
        if (K_is_required_) {
        }
    } else if (df_storage_ == disk) {
        //Not implemented yet
    } else {
        throw std::domain_error("Storage algorithm is not correct for L-DF-SCF");
    }

    if (J_is_required_) 
        free(DD);
    if (K_is_required_) 
        free_block(Cocc);
    
    //J_->print(outfile);
    //K_->print(outfile);

    /* FORM G_ */
    //This method takes one extra O(N^2) scale,
    //but preserves J and K in place
    G_->zero();
    G_->add(J_);
    G_->scale(-2.0);
    G_->add(K_);
    G_->scale(-1.0);
    //G_->print(outfile);
    timer_off("Overall G");
}
void RHF::propagate_local_mos()
{
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Must run in C1 for now.\n"); fflush(outfile);
        abort();
    } 
    int norbs = basisset_->nbf(); 
    //Get the C matrix (occupied only) 
    int ndocc = doccpi_[0];
    double **Cocc = block_matrix(norbs,ndocc);
    for (int i=0; i<norbs; i++) 
        for (int j=0; j<ndocc; j++)
            Cocc[i][j] = C_->get(0,i,j);
    
    double **T = block_matrix(norbs,norbs);
    {
        SharedMatrix temp(factory_.create_matrix("Temp"));
        //T = S_{mn}L_{ni}
        temp->gemm(false,false,1.0,S_,Lref_,0.0);

        for (int i=0; i<norbs; i++) 
            for (int j=0; j<norbs; j++)
                T[i][j] = temp->get(0,i,j);
    }

    double **W = block_matrix(ndocc,norbs); 
    //W = W_{oi} = C'T = C_{om}T_{ni}
    C_DGEMM('T','N',ndocc,norbs,norbs,1.0,Cocc[0],ndocc,T[0],norbs,0.0,W[0],norbs);
    free_block(T);
    
    double **Q = block_matrix(norbs,norbs);
    //Q = Q_{ii} = W'W = W_{io}W_{oi}
    C_DGEMM('T','N',norbs,norbs,ndocc,1.0,W[0],norbs,W[0],norbs,0.0,Q[0],norbs);

    double **Q12 = block_matrix(norbs,norbs);
    //Q^{-1/2} = Q_{ii}^{-1/2} via eigendecomposition

    // Form Q^-1/2
    // First, diagonalize Q
    // the C_DSYEV call replaces the original matrix Q with its eigenvectors
    double* eigval = init_array(norbs);
    int lwork = norbs * 3;
    double* work = init_array(lwork);
    int stat = C_DSYEV('v','u',norbs,Q[0],norbs,eigval, work,lwork);
    if (stat != 0) {
        fprintf(outfile, "C_DSYEV failed\n");
        exit(PSI_RETURN_FAILURE);
    }
    free(work);

    // Now Q contains the eigenvectors of the original Q
    // Copy Q to Q_copy
    double **Q_copy = block_matrix(norbs, norbs);
    C_DCOPY(norbs*norbs,Q[0],1,Q_copy[0],1); 

    // Now form Q^{-1/2} = U(T)*q^{-1/2}*U,
    // where q^{-1/2} is the diagonal matrix of the inverse square roots
    // of the eigenvalues, and U is the matrix of eigenvectors of Q
    for (int i=0; i<norbs; i++) {
        if (fabs(eigval[i]) < 1.0E-10)
            eigval[i] = 0.0;
        else 
            eigval[i] = 1.0 / sqrt(eigval[i]);

        // scale one set of eigenvectors by the diagonal elements q^{-1/2}
        C_DSCAL(norbs, eigval[i], Q[i], 1);
    }
    free(eigval);

    // Q_mhalf = Q_copy(T) * J
    C_DGEMM('t','n',norbs,norbs,norbs,1.0,
            Q_copy[0],norbs,Q[0],norbs,0.0,Q12[0],norbs);

    free_block(Q);
    free_block(Q_copy);
    //Or you could Q12 = Q^{1/2} in MATLAB!

    double **U = block_matrix(ndocc,ndocc);
    //U = U_{oo} = WQ^{-1/2} = W_{oi}Q_{ii}^{-1/2}
    C_DGEMM('N','N',ndocc,ndocc,norbs,1.0,W[0],norbs,Q12[0],norbs,0.0,U[0],ndocc);
    free_block(W);
    free_block(Q12);

    double **L = block_matrix(norbs,ndocc);

    //L = L_{mi} = CU = C_{mo}U_{oo}
    C_DGEMM('N','N',norbs,ndocc,ndocc,1.0,Cocc[0],ndocc,U[0],ndocc,0.0,L[0],ndocc);

    free_block(Cocc);
    free_block(U);
    
    for (int i=0; i<norbs; i++) 
        for (int j=0; j<ndocc; j++)
            L_->set(0,i,j,L[i][j]);
    free_block(L);

    L_->print(outfile);
}
void RHF::fully_localize_mos()
{
  //Pipek-Mizey procedure ripped from DF-MP2, which in turn was ripped from 
  //Localize in Psi3
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Must run in C1 for now.\n"); fflush(outfile);
        abort();
    } 

  //This one runs in O(N^3) after work in Summer 2009
  //Old localize had an O(N^5) step in the Givens rotations
  int iter, s, t, A, k, l, iold, max;
  int i, j, ij, am, atom, shell_length, offset;
  int ntri, *soccpi, *stype, *snuc;
  //double  **LCtmp, **F_occ;
  //double *scratch, *evals;
  //int *orb_order, *orb_boolean;
  double P, PiiA, Pst, Pss, Ptt, Ast, Bst, AB;
  double Uss, Utt, Ust, Uts, LCks, LCkt, **V; //**V, **VV;
  double cos4a, alpha, alphamax, alphalast, conv;

  int norbs = basisset_->nbf();
  int nshell = basisset_->nshell();
  int natom = basisset_->molecule()->natom();
  int nocc = doccpi_[0];
  double **C = block_matrix(norbs,norbs);
  for (int i = 0; i<norbs; i++)
    for (int j = 0; j<norbs; j++)
      C[i][j] = C_->get(0,i,j);
  
  double **S = block_matrix(norbs,norbs);
  for (int i = 0; i<norbs; i++)
    for (int j = 0 ; j<norbs; j++)
      S[i][j] = S_->get(0,i,j);

  fprintf(outfile, "\nC Matrix in the AO basis:\n");
  print_mat(C, norbs, norbs, outfile);

  //evals = get_evals();
  bool puream = basisset_->has_puream();
  snuc = init_int_array(nshell);
  for (i = 0; i<nshell; i++)
    snuc[i] = basisset_->shell(i)->ncenter();

  fprintf(outfile, "Overlap Matrix");
  print_mat(S, norbs, norbs, outfile);

  // Compute the length of each AM block
  int *l_length = init_int_array(basisset_->max_am());
  l_length[0] = 1;
  for(l=1; l < (basisset_->max_am()); l++) {
    if(puream) l_length[l] = 2 * l + 1;
    else l_length[l] = l_length[l-1] + l + 1;
  }

  // Set up the atom->AO and AO->atom lookup arrays
  int* aostart = init_int_array(natom);
  int* aostop = init_int_array(natom);
  int* aostart_shell = init_int_array(natom);
  int* aostop_shell = init_int_array(natom);
  stype = init_int_array(nshell);
  for (i=0; i<nshell; i++)
    stype[i] = basisset_->shell(i)->am();  

  for(i=0,atom=-1,offset=0; i < nshell; i++) {
    am = stype[i] - 1;                  // am is the angular momentum of the orbital
    shell_length = l_length[am];        // shell_length is the number of obritals in each shell

    if(atom != snuc[i]-1) {             // snuc is the nucleus that the shell belongs to
      if(atom != -1) {
        aostop[atom] = offset-1;
        aostop_shell[atom] = i-1;
      }
      atom = snuc[i]-1;
      aostart[atom] = offset;
      aostart_shell[atom] = i;
    }

    offset += shell_length;
  }
  aostop[atom] = offset-1;
  aostop_shell[atom] = i-1;

  int* ao2atom = init_int_array(norbs);
  for(i=0; i < natom; i++)
    for(j=aostart[i]; j <= aostop[i]; j++) {
      ao2atom[j] = i;                   // ao2atom is the atom number that the AO is located on
    }

  fprintf(outfile, "\tNumber of doubly occupied orbitals: %d\n\n", nocc);

  fprintf(outfile, "\tIter     Pop. Localization   Max. Rotation Angle       Conv\n");
  fprintf(outfile, "\t------------------------------------------------------------\n");

  V = block_matrix(nocc, nocc);
  double* tvec = init_array(nocc);
  double* svec = init_array(nocc);
  
  for(iter=0; iter < 100; iter++) {

  P = 0.0;
  for(i=0; i < nocc; i++) {
    for(A=0; A < natom; A++) {
      PiiA = 0.0;

      for(l=aostart[A]; l <= aostop[A]; l++)
        for(k=0; k < norbs; k++)
          PiiA += C[k][i] * C[l][i] * S[k][l];

      P += PiiA * PiiA;
    }
  }

  // Compute 2x2 rotations for Pipek-Mezey lo.zation
  alphamax = 0.0;

  for(s=0; s < nocc; s++) {
    for(t=0; t < s; t++) {

      Ast = Bst = 0.0;

      for(A=0; A < natom; A++) {

        Pst = Pss = Ptt = 0.0;

        for(l=aostart[A]; l <= aostop[A]; l++) {
          for(k=0; k < norbs; k++) {
            Pst += 0.5 * (C[k][s] * C[l][t] +
                          C[l][s] * C[k][t]) * S[k][l];            // Eqn 31 (JCP 90, 4916)

            Pss += C[k][s] * C[l][s] * S[k][l];                    // Eqn 31 (JCP 90, 4916)

            Ptt += C[k][t] * C[l][t] * S[k][l];                    // Eqn 31 (JCP 90, 4916)
          }
        }

        Ast += Pst * Pst - 0.25 * (Pss - Ptt) * (Pss - Ptt);               // Eqn 29A (JCP 90, 4916)
        Bst += Pst * (Pss - Ptt);                                  // Eqn 29B (JCP 90, 4916)

      } // A-loop

      // Compute the rotation angle
      AB = Ast * Ast + Bst * Bst;

        if(fabs(AB) > 0.0) {
          cos4a = -Ast/sqrt(AB);                                     // Eqn 13b (JCP 90, 4916)
          alpha = 0.25 * acos(cos4a) * (Bst > 0 ? 1 : -1);
        }
        else alpha = 0.0;

        // Keep up with the maximum 2x2 rotation angle
        alphamax = (fabs(alpha) > alphamax ? alpha : alphamax);


        Uss = cos(alpha);                                            // Eqn 10a/b (JCP 90, 4916)
        Utt = cos(alpha);                                               // Eqn 10a/b (JCP 90, 4916)
        Ust = sin(alpha);                                               // Eqn 10a/b (JCP 90, 4916)
        Uts = -Ust;                                                  // Eqn 10a/b (JCP 90, 4916)

        // Now do the rotation
        for(k=0; k < norbs; k++) {
          LCks = C[k][s];
          LCkt = C[k][t];
          C[k][s] = Uss * LCks + Ust * LCkt;
          C[k][t] = Uts * LCks + Utt * LCkt;
        }
        for (i = 0; i< nocc; i++)
        {
          svec[i] = Uss*V[i][s] + Ust*V[i][t];
          tvec[i] = Uts*V[i][s] + Utt*V[i][t];
        }
        for (i = 0; i< nocc; i++)
        {
          V[i][s] = svec[i]; 
          V[i][t] = tvec[i];
        }
      } // t-loop
    } // s-loop

    conv = fabs(alphamax) - fabs(alphalast);
    fprintf(outfile, "\t%4d  %20.10f  %20.10f  %4.3e\n", iter, P, alphamax, conv);
    if((iter > 2) && ((fabs(conv) < 1E-12) || alphamax == 0.0)) break;
    alphalast = alphamax;

    fflush(outfile);

  } // iter-loop
  free(tvec);
  free(svec);
  free_block(V);

  fprintf(outfile, "\nC Matrix in the LO basis:\n");
  print_mat(C, norbs, norbs, outfile);

  //free(evals);
  free(aostart);
  free(aostop);
  free(aostart_shell);
  free(aostop_shell);
  free(ao2atom);
  free(l_length);

    for (int i = 0; i<norbs; i++)
        for (int j = 0; j<norbs; j++) {
            L_->set(0,i,j,C[i][j]);
        }
    Lref_->copy(L_);
    free_block(C);
}
void RHF::localized_Lodwin_charges()
{
    L_->print(outfile);

    //Compute Lodwin atomic charges of localized orbitals (for L-DF-SCF)
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Must run in C1 for now.\n"); fflush(outfile);
        abort();
    } 
    int ndocc = doccpi_[0];
    int norbs = basisset_->nbf();
    int natom = basisset_->molecule()->natom();
   
    Shalf_->print(outfile); 
    Sphalf_->print(outfile);

    double **Q = block_matrix(norbs,ndocc);
    {
        SharedMatrix temp(factory_.create_matrix("Temp"));
        //Q = Q_{mo} = S^{+1/2}C = S_{mn}^{+1/2}C_{mo}
        temp->gemm(false,false,1.0,Sphalf_,L_,0.0);
        for (int m = 0; m<norbs; m++)
            for (int o = 0; o<ndocc; o++)
                Q[m][o] = temp->get(0,m,o);
        temp->print(outfile);
    }    
    //Q^2 = 2*[Q_{mo}]^2 //two for spin degeneracy
    for (int m = 0; m<norbs; m++)
        for (int o = 0; o<ndocc; o++)
            Q[m][o] = 2.0*Q[m][o]*Q[m][o];
    
    //Assume I_ is already initialized at natom x ndocc    
    zero_mat(I_,natom,ndocc);
    
    int atom_index, nummu, mu1,m;
    for (int MU = 0; MU<basisset_->nshell(); MU++) {
        atom_index = basisset_->shell(MU)->ncenter();
        nummu = basisset_->shell(MU)->nfunction();
        mu1 = basisset_->shell(MU)->function_index();
        for (int mu = 0; mu<nummu; mu++)
            m = mu1+mu;
            for (int o = 0; o<ndocc; o++)
                I_[atom_index][o] += Q[m][o];
    }
    
    free_block(Q);
    
    fprintf(outfile,"  Lodwin Atomic Charges");
    print_mat(I_,natom,ndocc,outfile);

}
}}
