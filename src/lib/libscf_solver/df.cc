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

void HF::form_B()
{
    fprintf(outfile, "\n  Computing Integrals using Density Fitting\n");
    //TODO: Add support for molecular symmetry
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Must run in C1 for now.\n"); fflush(outfile);
        abort();
    } 
    int norbs = basisset_->nbf(); 
    ribasis_ =shared_ptr<BasisSet>(new BasisSet(chkpt_, "DF_BASIS_SCF"));
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
    int ndocc = doccpi_[0];
    //Size of the exchange tensor
    unsigned long memC = norbs*ndocc*(long)ri_nbf_;

    string storage_type;
    storage_type = options_.get_str("RI_STORAGE");

    if (storage_type == "DOUBLE_IN_CORE")
        df_storage_ = double_full;
    else if (storage_type == "IN_CORE")
        df_storage_ = full;
    else if (storage_type == "FLIP_B_DISK")
        df_storage_ = flip_B_disk;
    else if (storage_type == "K_IN_CORE")
        df_storage_ = k_incore;
    else if (storage_type == "DISK")
        df_storage_ = disk;
    else if (storage_type == "DEFAULT")
    {
    	//set df_storage_ semi-heuristically based on available memory
    	if (((long)((memA+memA)*(1.0+MEMORY_SAFETY_FACTOR)))<(memory_/sizeof(double)))
            df_storage_ = double_full; //Double in-core, including both (ab|P) tensors
    	else if (((long)((memA+memC)*(1.0+MEMORY_SAFETY_FACTOR)))<(memory_/sizeof(double)))
            df_storage_ = full; //Full in-core, including both (ab|P) tensors
        else if (((long)((memA)*(1.0+MEMORY_SAFETY_FACTOR)))<(memory_/sizeof(double)))
            if (((long)((memC)*(1.0+MEMORY_SAFETY_FACTOR)))<(memory_/sizeof(double)))
                df_storage_ = flip_B_disk; //Transpose B using disk scratch and core, leave it on disk
    	    else
                df_storage_ = disk; //TODO If Sieve makes B small, might be smaller than E
        else if (((long)((memC)*(1.0+MEMORY_SAFETY_FACTOR)))<(memory_/sizeof(double)))
            df_storage_ = k_incore; //K only in-core
    	else
            df_storage_ = disk; //Disk
    }	

    if (df_storage_ == double_full)
        fprintf(outfile,"\n  Density Fitting Algorithm proceeding In Core, Fast transform.\n"); 
    if (df_storage_ == full)
        fprintf(outfile,"\n  Density Fitting Algorithm proceeding In Core, Slow Transform.\n"); 
    else if (df_storage_ == flip_B_disk)
        fprintf(outfile,"\n  Density Fitting Algorithm proceeding with K In Core, B On Disk, Transpose in core.\n");
    else if (df_storage_ == k_incore)
        fprintf(outfile,"\n  Density Fitting Algorithm proceeding with K In Core, B on Disk.\n");
    else if (df_storage_ == disk)
        fprintf(outfile,"\n  Density Fitting Algorithm proceeding on Disk\n"); 
    fflush(outfile);

    //It takes a lot of work to get a null basis with Psi4 
    shared_ptr<BasisSet> zero = BasisSet::zero_basis_set();
    
    // Create integral factory for J (Fitting Matrix in form_B)
    IntegralFactory rifactory_J(ribasis_, zero, ribasis_, zero);
    shared_ptr<TwoBodyInt> Jint = shared_ptr<TwoBodyInt>(rifactory_J.eri());

    // Integral buffer
    const double *Jbuffer = Jint->buffer();

    // J Matrix
    double **J = block_matrix(ri_nbf_, ri_nbf_);
    // J^{-1/2}
    double **J_mhalf = block_matrix(ri_nbf_, ri_nbf_);
    
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

                    J[omu][onu] = Jbuffer[index];
                    J[onu][omu] = Jbuffer[index];
                }
            }
        }
    }
    //fprintf(outfile,"\nJ:\n"); fflush(outfile);
    //print_mat(J,ri_nbf_,ri_nbf_,outfile);
    timer_off("Form J Matrix;");
    timer_on("Form J^-1/2;");

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
    timer_off("Form J^-1/2;");


    //fprintf(outfile,"\nJmhalf:\n"); fflush(outfile);
    //print_mat(J_mhalf,ri_nbf_,ri_nbf_,outfile);
    timer_on("Overall (B|mn)");
    
    //Use ri_pair_mu_ and ri_pair_nu_ to keep track of things
    //Across schwarz sieve and unfortunate shell indexing
    ri_pair_nu_ = init_int_array(ntri_naive_);
    ri_pair_mu_ = init_int_array(ntri_naive_);
  
    double three_index_cutoff = options_.get_double("THREE_INDEX_CUTOFF");
 
    if (df_storage_ == double_full)
    {
    	IntegralFactory rifactory(basisset_, basisset_, ribasis_, zero);
        shared_ptr<TwoBodyInt> eri = shared_ptr<TwoBodyInt>(rifactory.eri());
        const double *buffer = eri->buffer();
        double** A_ia_P = block_matrix(ri_nbf_,ntri_naive_); 

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
                        timer_on("(B|mn) Integrals");
                        eri->compute_shell(MU, NU, Pshell, 0);
                        timer_off("(B|mn) Integrals");
                        l_index = start_index;
                        for (mu=0 ; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                            for (nu=0; nu < numnu; ++nu) {
                                onu = basisset_->shell(NU)->function_index() + nu;
                                if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                    for (P=0; P < numP; ++P) {
                                        PHI = ribasis_->shell(Pshell)->function_index() + P;
                                        A_ia_P[PHI][l_index]= buffer[mu*numnu*numP+nu*numP+P];
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
        
 
        B_ia_P_ = block_matrix(ri_nbf_,ntri_naive_); 
        timer_on("(B|mn) Transform");
        C_DGEMM('N','N',ri_nbf_,ntri_naive_,ri_nbf_,1.0, J_mhalf[0], ri_nbf_, A_ia_P[0], ntri_naive_,
                0.0, B_ia_P_[0], ntri_naive_);
        timer_off("(B|mn) Transform");
        
	free_block(A_ia_P);
        //print_mat(B_ia_P_, ri_nbf_,ntri_naive_ ,outfile);

        if (three_index_cutoff>0.0) {
            int left =  0;
            int right = 0;
            bool negligible_col;
            for (right = 0; right<ntri_naive_; right++) {
                negligible_col = true;
                for (int Q = 0; Q<ri_nbf_; Q++) {
                    if (fabs(B_ia_P_[Q][right])>three_index_cutoff) {
                        negligible_col = false;
                        break;
                    }
                }
                if (!negligible_col) {
                    for (int Q = 0; Q<ri_nbf_; Q++)
                        B_ia_P_[Q][left] = B_ia_P_[Q][right];
                    ri_pair_mu_[left] = ri_pair_mu_[right];
                    ri_pair_nu_[left] = ri_pair_nu_[right]; 
                    left++; 
                } else {
                    ntri_--;
                }
            }
        }
        //print_mat(B_ia_P_, ri_nbf_,ntri_ ,outfile);
        //for (int left = 0; left<ntri_; left++)
        //    fprintf(outfile,"  %d pair: (%d, %d)\n",left,ri_pair_mu_[left],ri_pair_nu_[left]);

        //fflush(outfile);
    } 
    if (df_storage_ == full||df_storage_ == flip_B_disk)
    {	
    	IntegralFactory rifactory(basisset_, basisset_, ribasis_, zero);
        shared_ptr<TwoBodyInt> eri = shared_ptr<TwoBodyInt>(rifactory.eri());
        const double *buffer = eri->buffer();
        B_ia_P_ = block_matrix(ri_nbf_,ntri_naive_); 
        
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
                        timer_on("(B|mn) Integrals");
                        eri->compute_shell(MU, NU, Pshell, 0);
                        timer_off("(B|mn) Integrals");
                        l_index = start_index;
                        for (mu=0 ; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                            for (nu=0; nu < numnu; ++nu) {
                                onu = basisset_->shell(NU)->function_index() + nu;
                                if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                        
                                    for (P=0; P < numP; ++P) {
                                        PHI = ribasis_->shell(Pshell)->function_index() + P;
                                        B_ia_P_[PHI][l_index]= buffer[mu*numnu*numP+nu*numP+P];
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

        double **Temp1;
	double **Temp2;
	bool allocated = false;
	int unit = ri_nbf_; //May want to make this smaller for memory!!
	for (index = 0; index<ntri_naive_; index+=ri_nbf_)
	{
            int cols = unit;
            if (index+unit>=ntri_naive_) {
                cols = ntri_naive_-index;
                if (allocated) {
                    free_block(Temp1);
                    free_block(Temp2);
                    allocated = false;
                }
            }

            if (!allocated) {
                Temp1 = block_matrix(ri_nbf_,cols);
                Temp2 = block_matrix(ri_nbf_,cols);
            }
            for (int r = 0; r<ri_nbf_; r++)
                for (int c = index; c<index+cols; c++)
                    Temp1[r][c-index] = B_ia_P_[r][c];

            timer_on("(B|mn) Transform");
            C_DGEMM('N','N',ri_nbf_,cols,ri_nbf_,1.0, J_mhalf[0], ri_nbf_, Temp1[0], cols,0.0, Temp2[0],cols);
            timer_off("(B|mn) Transform");

            for (int r = 0; r<ri_nbf_; r++)
                for (int c = index; c<index+cols; c++)
                    B_ia_P_[r][c] = Temp2[r][c-index];


	}
	free_block(Temp1);
	free_block(Temp2);

        //print_mat(B_ia_P_,ri_nbf_,ntri_naive_,outfile);

        if (three_index_cutoff>0.0) {
            int left =  0;
            int right = 0;
            bool negligible_col;
            for (right = 0; right<ntri_naive_; right++) {
                negligible_col = true;
                for (int Q = 0; Q<ri_nbf_; Q++) {
                    if (fabs(B_ia_P_[Q][right])>three_index_cutoff) {
                        negligible_col = false;
                        break;
                    }
                }
                if (!negligible_col) {
                    for (int Q = 0; Q<ri_nbf_; Q++)
                        B_ia_P_[Q][left] = B_ia_P_[Q][right];
                    ri_pair_mu_[left] = ri_pair_mu_[right];
                    ri_pair_nu_[left] = ri_pair_nu_[right]; 
                    left++; 
                } else {
                    ntri_--;
                }
            }
        }
        //print_mat(B_ia_P_, ri_nbf_,ntri_ ,outfile);
        //for (int left = 0; left<ntri_; left++)
            //fprintf(outfile,"  %d pair: (%d, %d)\n",left,ri_pair_mu_[left],ri_pair_nu_[left]);

        fflush(outfile);

        if (df_storage_ == flip_B_disk)
        {
            write_B();
            free_block(B_ia_P_);
        }
    }
    else if (df_storage_ == k_incore||df_storage_ == disk)
    {
        psio_->open(PSIF_DFSCF_BJI,PSIO_OPEN_NEW);
        psio_address next_PSIF_DFSCF_BJI = PSIO_ZERO;
        
        IntegralFactory rifactory(basisset_, basisset_, ribasis_,zero);
        shared_ptr<TwoBodyInt> eri = shared_ptr<TwoBodyInt>(rifactory.eri());
        const double *buffer = eri->buffer();
        int maxfun = 0;
        for (int m = 0; m<basisset_->nshell(); m++)
            if (maxfun<basisset_->shell(m)->nfunction())
                maxfun=basisset_->shell(m)->nfunction();
        int maxpairs = maxfun*maxfun;
        double **Temp1 = block_matrix(ri_nbf_,maxpairs);
        double **Temp2 = block_matrix(ri_nbf_,maxpairs);
	double *Temp3 = init_array(ri_nbf_);

        int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
        int start_index, delta_index, l_index,r_index, s_index;
        start_index = 0; l_index = 0, r_index = 0;
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU) {
                numnu = basisset_->shell(NU)->nfunction();
                //fprintf(outfile, "  MU = %d, NU = %d, Sig = %d\n",MU,NU,schwarz_shell_pairs[MU*(MU+1)/2+NU]); fflush(outfile);
                if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                    delta_index = 0;
                    for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
                        numP = ribasis_->shell(Pshell)->nfunction();
                        timer_on("(B|mn) Integrals");
                        eri->compute_shell(MU, NU, Pshell, 0);
                        timer_off("(B|mn) Integrals");
                        s_index = 0;
                        for (mu=0 ; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                            for (nu=0; nu < numnu; ++nu) {
                                onu = basisset_->shell(NU)->function_index() + nu;
                                if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] == 1) {
                                    for (P=0; P < numP; ++P) {
                                        PHI = ribasis_->shell(Pshell)->function_index() + P;
                                        Temp1[PHI][s_index]= buffer[mu*numnu*numP+nu*numP+P];
                                    }
                                    //if (Pshell == 0)
                                    //    fprintf(outfile,"  %MU = %d, NU = %d, mu = %d, nu = %d, omu = %d, onu = %d, l_index = %d, r_index = %d, s_index = %d\n",MU,NU,mu,nu,omu,onu, l_index,r_index,s_index);  
                                    if (Pshell == 0) {
                                        delta_index++;
                                        ri_pair_mu_[r_index] = omu;
                                        ri_pair_nu_[r_index] = onu;
                                        r_index++;
                                    }
                                    s_index++;
                                } 
                            }
                        }
                    }
                    //Delta index is the filled length of the Temp1 matrix 
                    //Transformation by fitting metric!
                    C_DGEMM('N','N',ri_nbf_,delta_index,ri_nbf_,1.0, J_mhalf[0], ri_nbf_, Temp1[0], maxpairs,0.0, Temp2[0],maxpairs);
        
                    for (int pair = 0; pair < delta_index; pair++) {
                        bool sig = false;
                        for (int Q = 0; Q<ri_nbf_;Q++) {
                            if (fabs(Temp2[Q][pair])>=three_index_cutoff)
                                sig = true;
                            Temp3[Q] = Temp2[Q][pair];
                        }
                        //fprintf(outfile,"  Pair = %d, Sig = %s\n",pair,(sig?"Yes":"No"));
                        if (sig) { 
                            timer_on("(B|mn) disk");
                            psio_->write(PSIF_DFSCF_BJI,"BJ Three-Index Integrals",(char *) Temp3,sizeof(double)*ri_nbf_,next_PSIF_DFSCF_BJI,&next_PSIF_DFSCF_BJI);
                            timer_off("(B|mn) disk");
                            ri_pair_mu_[l_index] = ri_pair_mu_[start_index+pair];
                            ri_pair_nu_[l_index] = ri_pair_nu_[start_index+pair];
                            l_index++; 
                        } else {
                            ntri_--;
                        }
                    }
                    start_index+=delta_index;
                }
            }
        } 
        //for (int left = 0; left<ntri_; left++)
         //   fprintf(outfile,"  %d pair: (%d, %d)\n",left,ri_pair_mu_[left],ri_pair_nu_[left]);

        free_block(Temp1);
        free_block(Temp2);
	free(Temp3);
        //fprintf(outfile,"\n  Through B on disk."); fflush(outfile);
        psio_->close(PSIF_DFSCF_BJI,1);
        timer_on("(B|mn) restripe");
        
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
        timer_off("(B|mn) restripe");

        //fprintf(outfile,"\n  BJ Restriped on disk.\n"); fflush(outfile);
    }
    timer_off("Overall (B|mn)");
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
void HF::write_B()
{
    psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
    psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
            
    for (int Q = 0; Q<ri_nbf_; Q++) {
        psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(B_ia_P_[Q][0]),sizeof(double)*ntri_,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
                
     }
     psio_->close(PSIF_DFSCF_BJ,1);
}
void HF::free_B()
{
    if (df_storage_ == full||df_storage_ == double_full)
        free(B_ia_P_);
    free(ri_pair_mu_);
    free(ri_pair_nu_);
}
void RHF::form_G_from_RI()
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
    double* DD = init_array(ntri_);
    
    for (int ij = 0; ij<ntri_; ij++) {
        DD[ij] = D_->get(0,ri_pair_mu_[ij],ri_pair_nu_[ij]); 
        if (ri_pair_mu_[ij] != ri_pair_nu_[ij])
            DD[ij] *= 2.0;
            //only A irrep at the moment!!
    }
    //Get the C matrix (exchange messes things up)
    int ndocc = doccpi_[0];
    double** Cocc = block_matrix(ndocc,norbs);
    for (int i=0; i<norbs; i++) {
        for (int j=0; j<ndocc; j++)
            Cocc[j][i] = C_->get(0,i,j);
        //only A irrep at the moment!!
    }

    if (df_storage_ == full || df_storage_ == double_full) {
        //B is in core, E will be in core, DGEMM everything
        if (J_is_required_) {
            /* COULOMB PART */
            //Coulomb convolution vector
            double *L = init_array(ri_nbf_);
            //Temporary J matrix
            double *J = init_array(ntri_);
            //DGEMV -> L:
            //L_Q = (Q|ls)*D_{ls}
            timer_on("J DDOT");
            C_DGEMV('N',ri_nbf_,ntri_,1.0,B_ia_P_[0],ntri_naive_,DD,1,0.0,L,1);
            timer_off("J DDOT");
            //DGEMV -> J:
            //J_{mn} = L_Q(Q|mn)
            timer_on("J DAXPY");
            C_DGEMV('T',ri_nbf_,ntri_,1.0,B_ia_P_[0],ntri_naive_,L,1,0.0,J,1);
            timer_off("J DAXPY");
            //Put everything in J_
            for (int ij = 0; ij < ntri_; ij++) {
                J_->set(0,ri_pair_mu_[ij],ri_pair_nu_[ij],J[ij]);
                J_->set(0,ri_pair_nu_[ij],ri_pair_mu_[ij],J[ij]);
            }

            //J_->print(outfile);

            free(L);
            free(J);
            //J_->print(outfile);
        }
        if (K_is_required_) {
            timer_on("Form E");
            /* EXCHANGE PART */
            //E exchange matrix
            double** E = block_matrix(norbs, ndocc*ri_nbf_);
            //QS temp matrix for DGEMM
            double** QS = block_matrix(ri_nbf_,norbs);
            //Temp matrix for DGEMM
            double** Temp = block_matrix(ndocc,ri_nbf_);
            // Temp matrix for sparse DGEMM if sieve exists
            double** Ctemp = block_matrix(ndocc,norbs);
            // Index array for non-canonical ordering of mn
            int* m_indices = init_int_array(norbs);

            //Exchange tensor E
            for (int m = 0; m<norbs; m++) {
                //Find out where the m's are!
                int index = 0;
                for (int ij =0; ij<ntri_; ij++) {
                    if (ri_pair_mu_[ij] == m || ri_pair_nu_[ij] == m) 
                        m_indices[index++] = ij; //ns where the current m hits
                }

                int n_m = index;
                /**
                fprintf(outfile,"\n  nu = %d, ij indices are: ",nu);
                for (int k = 0; k<norbs; k++)
                    fprintf(outfile,"%d ",nu_indices[k]);
                fprintf(outfile,"\n");
                **/
                int n, ij;

                for (index = 0; index<n_m; index++) {
                    ij = m_indices[index];
                    if (ri_pair_nu_[ij] == m)
                        n = ri_pair_mu_[ij];
                    else
                        n = ri_pair_nu_[ij];
                
                    //fprintf(outfile,"  index = %d ij = %d \n",index,ij); fflush(outfile);//fprintf(outfile,"(ij, mu) = (%d, %d)\n",ij,mu);
                    for (int Q = 0; Q<ri_nbf_; Q++) {
                         QS[Q][index] = B_ia_P_[Q][ij];
                         //fprintf(outfile," ij = %d, mu = %d, Q = %d, val = %14.10f\n",ij,mu,Q,QS[Q][mu]);
                    }
                    for (int o = 0; o<ndocc; o++) {
                        Ctemp[o][index] = Cocc[o][n];
                    }
                }
                
                C_DGEMM('N','T',ndocc,ri_nbf_,n_m,1.0,Ctemp[0],norbs,QS[0],norbs, 0.0, Temp[0], ri_nbf_);
                
                int offset;
                for (int Q = 0; Q<ri_nbf_; Q++) {
                    offset = Q*ndocc;
                    for (int i = 0; i<ndocc; i++) {
                        E[m][i+offset] = Temp[i][Q];
                    }
                }
            }    
            free_block(Ctemp);
            free(m_indices);
            timer_off("Form E");
            
            //fprintf(outfile,"\n E: \n");
            //print_mat(E,norbs,ndocc*ri_nbf_,outfile);
            
            free_block(Temp);
            free_block(QS);
            timer_on("E DGEMM");

            //K_{mn} = E_{im}^QE_{in}^Q

            //What is the smallest overlap element worth computing exchange for
            double cutoff = options_.get_double("OVERLAP_CUTOFF");
            double contribution;
    
            int att_elements = norbs*(norbs+1)/2;

            if (cutoff > 0.0) {
                //If the K matrix is sparse due to low overlap or density
                //Form the matrix elements individually
                for (int i = 0, ij = 0; i< norbs; i++)
                    for (int j = 0; j <= i; j++, ij++) {
                        //Is S big enough?
                        if (abs(S_->get(0,i,j))>=cutoff) {
                            contribution = C_DDOT(ri_nbf_*ndocc,&E[i][0],1,&E[j][0],1);
                            K_->set(0,j,i,contribution);
                            K_->set(0,i,j,contribution);
                            att_elements--;
                        }
                    }
                fprintf(outfile,"  K matrix was built with with an overlap cutoff of %14.10E\n",cutoff);  
                fprintf(outfile,"  %d out of %d elements were attenuated due to overlap, %8.5f%% attenuation.\n",att_elements,norbs*(norbs+1)/2,100.0*att_elements/(1.0*norbs*(norbs+1)/2));  
            } else {
                //DGEMM, usually faster regardless, maybe if the filling of the overlap matrix is <1/6 this will be better
                //Temporary K matrix    
                double** K = block_matrix(norbs,norbs);
                //There it is
                C_DGEMM('N','T',norbs,norbs,ri_nbf_*ndocc,1.0,E[0],ri_nbf_*ndocc,E[0],ri_nbf_*ndocc, 0.0, K[0], norbs);
                for (int i = 0; i < norbs; i++)
                    for (int j = 0; j<=i; j++) {
                        K_->set(0,j,i,K[i][j]);
                        K_->set(0,i,j,K[i][j]);
                }
                free_block(K);
            }
            free_block(E);
            timer_off("E DGEMM");
            //K_->print(outfile);
        }
    }
    else if (df_storage_ == flip_B_disk || df_storage_ == k_incore) {
        //B is on disk, E will be in core, Single disk pass
        //in_buffer stores multiple aux basis function rows of 
        //the three index tensor
        int max_rows = floor(((memory_/sizeof(double))-norbs*ndocc*ri_nbf_)/((1.0+MEMORY_SAFETY_FACTOR)*ntri_));
	if (max_rows>ri_nbf_)
            max_rows = ri_nbf_;
        
        //Row height per read, so that we can tune this value
        int rows_per_read = options_.get_int("ROWS_PER_READ"); 
        
        double *in_buffer = init_array(max_rows*ntri_);
        //Transformed integrals are stored in PSIF_DFSCF_BJ
        //Which better exist at this point
        timer_on("Open B");
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        timer_off("Open B");
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;

        //Coulomb convolution vector, done element by element
        double L;
        double *J;
        if (J_is_required_){
            J = init_array(ntri_);
        }

        //Three index tensor, in core
        double **E;
        double **QS;
        int *m_indices;
        double **Temp; 
        double **Ctemp;
        double **K;
        if (K_is_required_){
            Ctemp = block_matrix(ndocc,norbs);
            E = block_matrix(norbs, ndocc*ri_nbf_);
            QS = block_matrix(max_rows,norbs);
            m_indices = init_int_array(norbs);
            Temp = block_matrix(ndocc,max_rows);
        }

        int mu, nu;
        int current_rows,offset;
        for (int row = 0; row <ri_nbf_; row+=max_rows)
        {
	    current_rows = max_rows;
	    if (row+max_rows>ri_nbf_)
		current_rows = ri_nbf_-row;
            //Read max_rows of the (B|mn) tensor in, place in in_buffer
            timer_on("Read B");
            
            //Old method, one go, fast reads but sometime Linux makes a 
            //huge buffer -> swaps and lots of hard faults
            //psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2*current_rows,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            //New method read in a few rows at a time
            int block_height = rows_per_read;
            for (int Q = 0; Q<current_rows; Q+=rows_per_read) {
                if (Q+rows_per_read>current_rows)
                    block_height = current_rows-Q;
                psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[Q*ntri_]),sizeof(double)*ntri_*block_height,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
                
            }
            
            timer_off("Read B");
            /**
            double** tempB = block_matrix(current_rows,norbs*(norbs+1)>>1);
            C_DCOPY(current_rows*norbs*(norbs+1)>>1,&in_buffer[0],1,&tempB[0][0],1);
            fprintf(outfile,"\n  B Temp: \n");
            print_mat(tempB,current_rows,norbs*(norbs+1)>>1,outfile);
            free_block(tempB);
            fprintf(outfile,"\n  Block indices are: ",nu);
            for (int k = 0; k<norbs*(norbs+1)>>1; k++)
                fprintf(outfile,"(%d, %d) ",ri_pair_mu_[k],ri_pair_nu_[k]);
            fprintf(outfile,"\n");
            **/

            if (J_is_required_) {
                for (int Q = row; Q< row+current_rows; Q++) {
		    offset = (Q-row)*ntri_;
                    /* COULOMB PART */
                    //L_Q = (Q|ls)D_{ls}
                    timer_on("J DDOT");
                    L = C_DDOT(ntri_,&in_buffer[offset],1,DD,1);
                    timer_off("J DDOT");
                    //J_{mn} += (Q|mn)L_Q
                    timer_on("J DAXPY");
                    C_DAXPY(ntri_,L,&in_buffer[offset],1,J,1);
                    timer_off("J DAXPY");
                }
            }
            if (K_is_required_) {
                timer_on("Form E");
                //Exchange tensor E
                for (int m = 0; m<norbs; m++) {
                    //Find out where the m's are!
                    int index = 0;
                    for (int ij =0; ij<ntri_; ij++) {
                        if (ri_pair_mu_[ij] == m || ri_pair_nu_[ij] == m) 
                            m_indices[index++] = ij; //ns where the current m hits
                    }

                    int n_m = index;
                    /**
                    fprintf(outfile,"\n  nu = %d, ij indices are: ",nu);
                    for (int k = 0; k<norbs; k++)
                        fprintf(outfile,"%d ",nu_indices[k]);
                    fprintf(outfile,"\n");
                    **/
                    int n, ij;

                    for (index = 0; index<n_m; index++) {
                        ij = m_indices[index];
                        if (ri_pair_nu_[ij] == m)
                            n = ri_pair_mu_[ij];
                        else
                            n = ri_pair_nu_[ij];
                
                         //fprintf(outfile,"  index = %d ij = %d \n",index,ij); fflush(outfile);//fprintf(outfile,"(ij, mu) = (%d, %d)\n",ij,mu);
                        for (int Q = 0; Q<current_rows; Q++) {
                             QS[Q][index] = in_buffer[Q*ntri_+ij];
                            //fprintf(outfile," ij = %d, mu = %d, Q = %d, val = %14.10f\n",ij,mu,Q,QS[Q][mu]);
                        }
                        for (int o = 0; o<ndocc; o++) {
                            Ctemp[o][index] = Cocc[o][n];
                        }
                    }
                
                    C_DGEMM('N','T',ndocc,current_rows,n_m,1.0,Ctemp[0],norbs,QS[0],norbs, 0.0, Temp[0], max_rows);
                
                    int offset;
                    for (int Q = 0; Q<current_rows; Q++) {
                        offset = (Q+row)*ndocc;
                        for (int i = 0; i<ndocc; i++) {
                            E[m][i+offset] = Temp[i][Q];
                        }
                    }
                }    
                timer_off("Form E");
            }
               
            //fprintf(outfile,"\n E: \n");
            //print_mat(E,norbs,ndocc*max_rows,outfile);
            
            
        }
        free_block(Ctemp);
        free_block(QS);
        free(m_indices);
        free_block(Temp);
        free(in_buffer);
        psio_->close(PSIF_DFSCF_BJ,1);

        /* Form J and */
        if (J_is_required_) {
            for (int ij = 0; ij < ntri_; ij++) {
                J_->set(0,ri_pair_mu_[ij],ri_pair_nu_[ij],J[ij]);
                J_->set(0,ri_pair_nu_[ij],ri_pair_mu_[ij],J[ij]);
            }
        }
        if (J_is_required_) {
            free(J);
        }
        
        timer_on("E DGEMM");
        /* Exchange Tensor DGEMM */
        if (K_is_required_) {
            //What is the smallest overlap element worth computing exchange for
            double cutoff = options_.get_double("OVERLAP_CUTOFF");
            double contribution;
    
            int att_elements = norbs*(norbs+1)/2;

            if (cutoff > 0.0) {
                //If the K matrix is sparse due to low overlap or density
                //Form the matrix elements individually
                for (int i = 0, ij = 0; i< norbs; i++)
                    for (int j = 0; j <= i; j++, ij++) {
                        //Is S big enough?
                        if (abs(S_->get(0,i,j))>=cutoff) {
                            contribution = C_DDOT(ri_nbf_*ndocc,&E[i][0],1,&E[j][0],1);
                            K_->set(0,j,i,contribution);
                            K_->set(0,i,j,contribution);
                            att_elements--;
                        }
                    }
                fprintf(outfile,"  K matrix was built with with an overlap cutoff of %14.10E\n",cutoff);  
                fprintf(outfile,"  %d out of %d elements were attenuated due to overlap, %8.5f%% attenuation.\n",att_elements,norbs*(norbs+1)/2,100.0*att_elements/(1.0*norbs*(norbs+1)/2));  
            } else {
                //DGEMM, usually faster regardless, maybe if the filling of the overlap matrix is <1/6 this will be better
                //Temporary K matrix    
                double** K = block_matrix(norbs,norbs);
                //There it is
                C_DGEMM('N','T',norbs,norbs,ri_nbf_*ndocc,1.0,E[0],ri_nbf_*ndocc,E[0],ri_nbf_*ndocc, 0.0, K[0], norbs);
                for (int i = 0; i < norbs; i++)
                    for (int j = 0; j<=i; j++) {
                        K_->set(0,j,i,K[i][j]);
                        K_->set(0,i,j,K[i][j]);
                }
                free_block(K);
            }
            free_block(E);
        }
        timer_off("E DGEMM");

    }
    else {
        //B is on disk, K will be in disk, Single disk pass
        //B is on disk, E will be in disk, Single disk pass
        //in_buffer stores multiple aux basis function rows of 
        //the three index tensor
        int max_rows = floor(((memory_/sizeof(double)))/((1.0+MEMORY_SAFETY_FACTOR)*(ntri_+ndocc*norbs)));
	if (max_rows>ri_nbf_)
            max_rows = ri_nbf_;

        //Row height per read, so that we can tune this value
        int rows_per_read = options_.get_int("ROWS_PER_READ"); 
        
        double *in_buffer = init_array(max_rows*ntri_);
        //Transformed integrals are stored in PSIF_DFSCF_BJ
        //Which better exist at this point
        timer_on("Open B");
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        timer_off("Open B");
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
        
        psio_address next_PSIF_DFSCF_K;
        if (K_is_required_) {
            //Exchange matrix is stored in PSIF_DFSCF_K
            //Which is created and destroyed in each iteration
            psio_->open(PSIF_DFSCF_K,PSIO_OPEN_NEW);
            next_PSIF_DFSCF_K = PSIO_ZERO;
        }
        //Coulomb convolution vector, done element by element
        double L;
        double *J;
        if (J_is_required_){
            J = init_array(ntri_);
        }

        //Three index tensor, buffered to disk
        double *out_buffer;
        double **QS;
        double **Ctemp;
        int *m_indices;
        double **Temp;
        double **K;
        if (K_is_required_){
            out_buffer = init_array(norbs*ndocc*max_rows);
            QS = block_matrix(max_rows,norbs);
            m_indices = init_int_array(norbs);
            Temp = block_matrix(ndocc,max_rows);
            Ctemp = block_matrix(ndocc,norbs);
        }

        int mu, nu;
        int current_rows,offset;
        for (int row = 0; row <ri_nbf_; row+=max_rows)
        {
	    current_rows = max_rows;
	    if (row+max_rows>ri_nbf_)
		current_rows = ri_nbf_-row;
            //Read a max_rows of the (B|mn) tensor in, place in in_buffer
            timer_on("Read B");
            
            //Old method, one go, fast reads but sometime Linux makes a 
            //huge buffer -> swaps and lots of hard faults
            //psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2*current_rows,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            //New method read in a few rows at a time
            int block_height = rows_per_read;
            for (int Q = 0; Q<current_rows; Q+=rows_per_read) {
                if (Q+rows_per_read>current_rows)
                    block_height = current_rows-Q;
                psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[Q*ntri_]),sizeof(double)*ntri_*block_height,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
                
            }

            //double **Bhack = block_matrix(ri_nbf_,ntri_);
            //C_DCOPY(ri_nbf_*ntri_,in_buffer,1,&Bhack[0][0],1);
            //fprintf(outfile,"  B:\n");
            //print_mat(Bhack,ri_nbf_,ntri_,outfile);
            
            timer_off("Read B");
            for (int Q = row; Q< row+current_rows; Q++) {
		//offset in three-index tensor
		int offset_B = (Q-row)*ntri_;
		//offset in E tensor
		int offset_E = (Q-row)*norbs*ndocc;

                if (J_is_required_) {
                    /* COULOMB PART */
                    //L_Q = (Q|ls)D_{ls}
                    timer_on("J DDOT");
                    L = C_DDOT(ntri_,&in_buffer[offset_B],1,DD,1);
                    timer_off("J DDOT");
                    //J_{mn} += (Q|mn)L_Q
                    timer_on("J DAPXY");
                    C_DAXPY(ntri_,L,&in_buffer[offset_B],1,J,1);
                    timer_off("J DAPXY");
                } 
                /*
                timer_on("Form E");
                if (K_is_required_) {
                    // EXCHANGE TENSOR 
                    memset(&out_buffer[offset_E],0,norbs*ndocc*sizeof(double));
                    for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
                    {
                        mu = ri_pair_mu_[ij];
                        nu = ri_pair_nu_[ij];
                        for (int i = 0; i<ndocc; i++)
                        {
                            out_buffer[mu*ndocc+i+offset_E]+=Cocc[i][nu]*in_buffer[ij+offset_B];
                            if (mu != nu)
                                out_buffer[nu*ndocc+i+offset_E]+=Cocc[i][mu]*in_buffer[ij+offset_B];
                        }
                    }
                }
                timer_off("Form E");**/
            }
            
            if (K_is_required_) {
                timer_on("Form E");
                //Exchange tensor E
                for (int m = 0; m<norbs; m++) {
                    //Find out where the m's are!
                    int index = 0;
                    for (int ij =0; ij<ntri_; ij++) {
                        if (ri_pair_mu_[ij] == m || ri_pair_nu_[ij] == m) 
                            m_indices[index++] = ij; //ns where the current m hits
                    }

                    int n_m = index;
                    /**
                    fprintf(outfile,"\n  nu = %d, ij indices are: ",nu);
                    for (int k = 0; k<norbs; k++)
                        fprintf(outfile,"%d ",nu_indices[k]);
                    fprintf(outfile,"\n");
                    **/
                    int n, ij;

                    for (index = 0; index<n_m; index++) {
                        ij = m_indices[index];
                        if (ri_pair_nu_[ij] == m)
                            n = ri_pair_mu_[ij];
                        else
                            n = ri_pair_nu_[ij];
                
                         //fprintf(outfile,"  index = %d ij = %d \n",index,ij); fflush(outfile);//fprintf(outfile,"(ij, mu) = (%d, %d)\n",ij,mu);
                        for (int Q = 0; Q<current_rows; Q++) {
                             QS[Q][index] = in_buffer[Q*ntri_+ij];
                            //fprintf(outfile," ij = %d, mu = %d, Q = %d, val = %14.10f\n",ij,mu,Q,QS[Q][mu]);
                        }
                        for (int o = 0; o<ndocc; o++) {
                            Ctemp[o][index] = Cocc[o][n];
                        }
                    }
                
                    C_DGEMM('N','T',ndocc,current_rows,n_m,1.0,Ctemp[0],norbs,QS[0],norbs, 0.0, Temp[0], max_rows);
                
                    int delta;
                    for (int Q = row; Q<row+current_rows; Q++) {
                        delta = (Q-row)*ndocc*norbs;
                        C_DCOPY(ndocc,&Temp[0][Q-row],max_rows,&out_buffer[delta+m*ndocc],1);
                    }
                }    
                timer_off("Form E");
                

                timer_on("Write E");
            
                //Old method, one go, fast reads but sometime Linux makes a 
                //huge buffer -> swaps and lots of hard faults
                //psio_->write(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(out_buffer[0]),sizeof(double)*norbs*ndocc*current_rows,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
                //New method read in a few rows at a time
                int block_height = rows_per_read;
                for (int Q = 0; Q<current_rows; Q+=rows_per_read) {
                    if (Q+rows_per_read>current_rows)
                        block_height = current_rows-Q;
                    psio_->write(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(out_buffer[Q*ndocc*norbs]),sizeof(double)*norbs*ndocc*block_height,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
                
                }
            
                timer_off("Write E");
            }
        }
        free(in_buffer);
        /* Form J */
        if (J_is_required_) {
            for (int ij2 = 0; ij2 < ntri_; ij2++) {
                J_->set(0,ri_pair_mu_[ij2],ri_pair_nu_[ij2],J[ij2]);
                J_->set(0,ri_pair_nu_[ij2],ri_pair_mu_[ij2],J[ij2]);
            }
        }
        if (J_is_required_) {
            free(J);
        }
        psio_->close(PSIF_DFSCF_BJ,1);
        if (K_is_required_) {
            free(out_buffer);
            free_block(QS);
            free_block(Ctemp);
            free(m_indices);
            free_block(Temp);
            psio_->close(PSIF_DFSCF_K,1);
        }

        /* Exchange Tensor DGEMM */
        if (K_is_required_) {
            psio_->open(PSIF_DFSCF_K,PSIO_OPEN_OLD);
            next_PSIF_DFSCF_K = PSIO_ZERO;

            K = block_matrix(norbs,norbs);
            //Large blocks implemented
            max_rows = floor(((memory_/sizeof(double)))/((1.0+MEMORY_SAFETY_FACTOR)*(2.0*ndocc*norbs)));
	    if (max_rows>ri_nbf_)
                max_rows = ri_nbf_;

            //What is the smallest overlap element worth computing exchange for
            double cutoff = options_.get_double("OVERLAP_CUTOFF");
            double contribution;
    
            in_buffer = init_array(norbs*(norbs+1)/2*max_rows);
            double **E = block_matrix(norbs,max_rows*ndocc);
            double **Ktemp = block_matrix(norbs,norbs);
            for (int row = 0; row <ri_nbf_; row+=max_rows)
            {
	        current_rows = max_rows;
	        if (row+max_rows>ri_nbf_)
		    current_rows = ri_nbf_-row;
                //Read max_rows of the exchange tensor in, place in in_buffer
                timer_on("E Read");
            
                //Old method, one go, fast reads but sometime Linux makes a 
                //huge buffer -> swaps and lots of hard faults
                //psio_->read(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(in_buffer[0]),sizeof(double)*norbs*ndocc*current_rows,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
                //New method read in a few rows at a time
                int block_height = rows_per_read;
                for (int Q = 0; Q<current_rows; Q+=rows_per_read) {
                    if (Q+rows_per_read>current_rows)
                        block_height = current_rows-Q;
                    psio_->read(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(in_buffer[Q*ndocc*norbs]),sizeof(double)*norbs*ndocc*block_height,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
                }
                timer_off("E Read");
                
                timer_on("E DGEMM");

                /**
                double** tempB = block_matrix(current_rows,norbs*ndocc);
                C_DCOPY(current_rows*norbs*ndocc,&in_buffer[0],1,&tempB[0][0],1);
                fprintf(outfile,"\n  E Temp: \n");
                print_mat(tempB,current_rows,norbs*ndocc,outfile);
                free_block(tempB);
                **/

                //Setup nice, get a semi-E matrix (blocked by aux basis)
                for (int m = 0; m<norbs; m++) 
                    for (int Q = 0; Q<current_rows; Q++)
                        C_DCOPY(ndocc,&(in_buffer[Q*ndocc*norbs+m*ndocc]),1,&(E[m][Q*ndocc]),1);

                //fprintf(outfile, "  E:\n");
                //print_mat(E,norbs,ndocc*ri_nbf_,outfile);                
    
                //Here's the actual DGEMM
                if (cutoff > 0.0) {
                    int att_elements = norbs*(norbs+1)/2;
                    
                    //If the K matrix is sparse due to low overlap or density
                    //Form the matrix elements individually
                    for (int i = 0, ij = 0; i< norbs; i++)
                        for (int j = 0; j <= i; j++, ij++) {
                            //Is S big enough?
                            if (abs(S_->get(0,i,j))>=cutoff) {
                                contribution = C_DDOT(current_rows*ndocc,&E[i][0],1,&E[j][0],1);
                                K[i][j] += contribution;
                                if (i != j)
                                    K[j][i] += contribution;
                                att_elements--;
                            }
                        }
                    fprintf(outfile,"  K matrix was built with with an overlap cutoff of %14.10E\n",cutoff);  
                    fprintf(outfile,"  %d out of %d elements were attenuated due to overlap, %8.5f%% attenuation.\n",att_elements,norbs*(norbs+1)/2,100.0*att_elements/(1.0*norbs*(norbs+1)/2));  
                } else {
                    //DGEMM, usually faster regardless, maybe if the filling of the overlap matrix is <1/6 this will be better
                    //Temporary K matrix    
                    //There it is
                    C_DGEMM('N','T',norbs,norbs,current_rows*ndocc,1.0,E[0],ri_nbf_*ndocc,E[0],ri_nbf_*ndocc, 0.0, Ktemp[0], norbs);

                    C_DAXPY(norbs*norbs,1.0,&Ktemp[0][0],1,&K[0][0],1);
                }
                                
                timer_off("E DGEMM"); 

                /**
                for (int i = 0; i<norbs; i++)
                    for (int j = 0; j<=i; j++) {
                        K[i][j] += C_DDOT(current_rows*ndocc,&(in_buffer[i]),norbs,&(in_buffer[j]),norbs);
                        if (i!=j)
                            K[j][i] += K[i][j];
                    }**/
                /**
                for (int Q = row; Q< row+current_rows; Q++) {
		    offset = (Q-row)*norbs*ndocc;

                    for (int m = 0; m<norbs; m++)
                        for (int n = 0; n<=m; n++)
                            for (int i = 0; i<ndocc; i++) {
                        K[m][n]+=in_buffer[m*ndocc+i+offset]*in_buffer[n*ndocc+i+offset];
                        K[n][m] = K[m][n];
                    }
                }**/
            }
            free(in_buffer);
            free_block(E);
            free_block(Ktemp);
            psio_->close(PSIF_DFSCF_K,0);
            
            for (int i = 0; i < norbs; i++)
                for (int j = 0; j<=i; j++) {
                    K_->set(0,j,i,K[i][j]);
                    K_->set(0,i,j,K[i][j]);
            }
            free_block(K);
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
    G_->add(K_);
    G_->scale(-1.0);
    //G_->print(outfile);
    timer_off("Overall G");

}
void RHF::form_J_from_RI()
{
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

    double *L = init_array(ri_nbf_);
    double *Gtemp = init_array(norbs*(norbs+1)/2);

    //B_ia_P_ in core
    if (df_storage_ == full||df_storage_ == double_full)
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
}
void RHF::form_K_from_RI()
{
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
    if (df_storage_ == full||df_storage_ == double_full)
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
        //print_mat(B_im_Q,norbs,ndocc*ri_nbf_,outfile);
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

    if (df_storage_ == full ||df_storage_ == double_full|| df_storage_ == flip_B_disk || df_storage_ == k_incore)
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

    for (int i=0; i<norbs; i++) {
        for (int j=0; j<=i; j++) {
            K_->set(0,i,j,K[i][j]);
            if (i!= j)
                K_->set(0,j,i,K[i][j]);
        }
    }
    //fprintf(outfile,"\n");
    //K_->print();
    free_block(K);
}
void UHF::form_G_from_RI()
{
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

    double *L = init_array(ri_nbf_);
    double *Gtemp = init_array(norbs*(norbs+1)/2);
    
    //B_ia_P_ in core
    if (df_storage_ == full||df_storage_ == double_full)
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
        {
            DD[ij] = D2[ioff[ri_pair_mu_[ij]]+ri_pair_nu_[ij]];
            //DD[ij] += D2[ioff[ri_pair_mu_[ij]]+ri_pair_nu_[ij]];
        }
    	free(D2);
    	
    	psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
    	psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
    	double *in_buffer = init_array(norbs*(norbs+1)/2);
    	for (int i=0; i<ri_nbf_; i++) {
            psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            L[i]=C_DDOT(norbs*(norbs+1)/2,DD,1,in_buffer,1);
      	}

    	free(DD);
    	psio_close(PSIF_DFSCF_BJ,1);
    	
    	psio_open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
    	next_PSIF_DFSCF_BJ = PSIO_ZERO;
    	double *G2 = init_array(norbs*(norbs+1)/2);
    	register double LL;
    	for (int Q = 0; Q<ri_nbf_; Q++)
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
    if (df_storage_ == full||df_storage_ == double_full)
    {
    	B_im_Q = block_matrix(norbs, nalpha*ri_nbf_);
    	double** QS = block_matrix(ri_nbf_,norbs);
    	double** Temp = block_matrix(nalpha,ri_nbf_);

        //print_mat(B_ia_P_,ri_nbf_,norbs*(norbs+1)/2,outfile);
        //fprintf(outfile,"\nYo\n");
        //print_mat(Cocc,nalpha,norbs,outfile);
        for (int m = 0; m<norbs; m++) {
            for (int Q = 0; Q<ri_nbf_; Q++) {
                for (int s = 0; s<norbs; s++) {
                    QS[Q][s] = B_ia_P_[Q][((s>=m)?ioff[s]+m:ioff[m]+s)];
                }
            }
            //C_DGEMM('N','T',ndocc,ri_nbf_,norbs,1.0,Cocc[0],norbs,QS[0],norbs, 0.0, Temp[0], ri_nbf_);
            C_DGEMM('N','T',nalpha,ri_nbf_,norbs,1.0,Cocc[0],norbs,QS[0],norbs, 0.0, Temp[0], ri_nbf_);
            for (int Q = 0; Q<ri_nbf_; Q++) {
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
    	B_im_Q = block_matrix(norbs, nalpha*ri_nbf_);
    	
    	int mu, nu;
    	
    	for (int Q = 0; Q<ri_nbf_; Q++)
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
    	for (int Q = 0; Q<ri_nbf_; Q++)
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

    if (df_storage_ == full ||df_storage_ == double_full|| df_storage_ == flip_B_disk || df_storage_ == k_incore)
    {
    	C_DGEMM('N','T',norbs,norbs,ri_nbf_*nalpha,1.0,B_im_Q[0],ri_nbf_*nalpha,B_im_Q[0],ri_nbf_*nalpha, 0.0, Ka[0], norbs);
    	free_block(B_im_Q);
    } else {
    	psio_open(PSIF_DFSCF_K,PSIO_OPEN_OLD);
    	psio_address next_PSIF_DFSCF_K = PSIO_ZERO;
    	
    	double *in_buffer = init_array(norbs*nalpha);
    	
    	for (int Q = 0; Q<ri_nbf_; Q++)
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
        if (df_storage_ == full||df_storage_ == double_full)
        {
            B_im_Q = block_matrix(norbs, nbeta_*ri_nbf_);
            double** QS = block_matrix(ri_nbf_,norbs);
            double** Temp = block_matrix(nbeta_,ri_nbf_);

            //print_mat(B_ia_P_,ri_nbf_,norbs*(norbs+1)/2,outfile);
            //fprintf(outfile,"\nYo\n");
            //print_mat(Cocc,nbeta_,norbs,outfile);
            for (int m = 0; m<norbs; m++) {
                for (int Q = 0; Q<ri_nbf_; Q++) {
                    for (int s = 0; s<norbs; s++) {
                        QS[Q][s] = B_ia_P_[Q][((s>=m)?ioff[s]+m:ioff[m]+s)];
                    }
                }
                C_DGEMM('N','T',nbeta_,ri_nbf_,norbs,1.0,Cocc[0],norbs,QS[0],norbs, 0.0, Temp[0], ri_nbf_);
                //print_mat(Temp,nbeta_,ri_nbf_,outfile);
                for (int Q = 0; Q<ri_nbf_; Q++) {
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
            B_im_Q = block_matrix(norbs, nbeta_*ri_nbf_);

            int mu, nu;

            for (int Q = 0; Q<ri_nbf_; Q++)
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
            for (int Q = 0; Q<ri_nbf_; Q++)
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

 	if (df_storage_ == full ||df_storage_ == double_full|| df_storage_ == flip_B_disk || df_storage_ == k_incore)
        {
            C_DGEMM('N','T',norbs,norbs,ri_nbf_*nbeta_,1.0,B_im_Q[0],ri_nbf_*nbeta_,B_im_Q[0],ri_nbf_*nbeta_, 0.0, Kb[0], norbs);
            free_block(B_im_Q);
        } else {
            psio_open(PSIF_DFSCF_K,PSIO_OPEN_OLD);
            psio_address next_PSIF_DFSCF_K = PSIO_ZERO;

            double *in_buffer = init_array(norbs*nbeta_);

            for (int Q = 0; Q<ri_nbf_; Q++)
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

}
void ROHF::form_G_from_RI()
{
    fprintf(stderr, "ROHF RI Not implemented yet!\n");
    abort();
}

}}
