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

#include <libmints/mints.h>

#include "hf.h"
#include "rhf.h"
#include "uhf.h"
#include "rohf.h"

using namespace std;
using namespace psi;

namespace psi { namespace scf {

void HF::form_A()
{
    fprintf(outfile, "\n  Computing Integrals using Density Fitting\n");
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Local SCF must run in C1.\n"); fflush(outfile);
        abort();
    } 
    int norbs = basisset_->nbf(); 
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
    ribasis_ = BasisSet::construct(parser, molecule_, options_.get_str("RI_BASIS_SCF"));
    naux_fin_ = ribasis_->nbf();
    
    //Form the schwarz sieve
    timer_on("Schwarz Sieve");

    sig_fun_pairs_ = 0;
    sig_shell_pairs_ = 0;

    schwarz_shell_pairs_ = init_int_array(basisset_->nshell()*(basisset_->nshell()+1)/2);
    schwarz_fun_pairs_ = init_int_array(norbs*(norbs+1)/2);
    if (schwarz_ > 0.0) {
        
        double* max_shell_val = init_array(basisset_->nshell()*(basisset_->nshell()+1)/2);
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
                schwarz_fun_pairs_[ij] = 1;
                sig_fun_pairs_++;
            }
        for (int ij = 0; ij < basisset_->nshell()*(basisset_->nshell()+1)/2; ij ++)
            if (max_shell_val[ij]*max_global_val>=schwarz_*schwarz_){
                schwarz_shell_pairs_[ij] = 1;
                sig_shell_pairs_++;
            }

        if (print_>4) {        
            for (int i = 0, ij = 0; i<norbs; i++)
                for (int j = 0; j<=i; j++, ij++)
                    fprintf(outfile,"   Function pair %d = (%d,%d), Max val %14.10f, Max Integral %14.10f, Significant %s\n",ij,i,j,max_fun_val[ij],max_fun_val[ij]*max_global_val,(schwarz_fun_pairs_[ij])?"YES":"NO");
            fprintf(outfile,"\n  Shell Pair Schwarz Sieve, schwarz_ = %14.10f:\n",schwarz_);
            for (int i = 0, ij = 0; i<basisset_->nshell(); i++)
                for (int j = 0; j<=i; j++, ij++)
                    fprintf(outfile,"   Shell pair %d = (%d,%d), Max val %14.10f, Max Integral %14.10f, Significant %s\n",ij,i,j,max_shell_val[ij],max_shell_val[ij]*max_global_val,(schwarz_shell_pairs_[ij])?"YES":"NO");
        fprintf(outfile, "\n");
        }

        free(max_fun_val);
        free(max_shell_val);
    
        ntri_naive_ = sig_fun_pairs_; //Matrix size for most of the algorithm
        ntri_ = ntri_naive_; //For now!

    } else {
        ntri_ = norbs*(norbs+1)/2; //Yeah, eat it 
        ntri_naive_ = norbs*(norbs+1)/2; 
        sig_fun_pairs_ = ntri_;
        sig_shell_pairs_ = basisset_->nshell()*(basisset_->nshell()+1)/2;
        for (int ij = 0; ij < basisset_->nshell()*(basisset_->nshell()+1)/2; ij++)
            schwarz_shell_pairs_[ij] = 1;
        for (int ij = 0; ij < ntri_; ij++)
            schwarz_fun_pairs_[ij] = 1;
    }

    timer_off("Schwarz Sieve");
    //Size of the three-index tensor
    unsigned long memA = ntri_*(long)naux_fin_;

    string storage_type;
    storage_type = options_.get_str("RI_SCF_STORAGE");

    if (storage_type == "CORE")
        df_storage_ = core;
    else if (storage_type == "DISK")
        df_storage_ = disk;
    else if (storage_type == "DEFAULT")
    {
    	//set df_storage_ semi-heuristically based on available memory
    	if (((long)((memA)*(1.0+MEMORY_SAFETY_FACTOR)))<(memory_/sizeof(double)))
            df_storage_ = core; //K only in-core
    	else
            df_storage_ = disk; //Disk
    } else {
        throw std::domain_error("For Form A, IN_CORE and DISK are the only valid storage options");
    }	

    if (df_storage_ == core)
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
    Winv_ = block_matrix(naux_fin_, naux_fin_);
    
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

                    Winv_[omu][onu] = Jbuffer[index];
                    Winv_[onu][omu] = Jbuffer[index];
                }
            }
        }
    }

    W_ = block_matrix(naux_fin_,naux_fin_);
    
    //Copy J to W_ for later use (elements of K)
    C_DCOPY(naux_fin_*naux_fin_,Winv_[0],1,W_[0],1);

    if (print_>4) {
        fprintf(outfile,"\nJ:\n"); fflush(outfile);
        print_mat(Winv_,naux_fin_,naux_fin_,outfile);
    }

    timer_off("Form J Matrix;");
    timer_on("Form J^-1");
    
    //fprintf(outfile,"  J\n");
    //print_mat(Winv_,naux_fin_,naux_fin_,outfile);

    //Cholesky facotrization (in place)
    int CholError = C_DPOTRF('L',naux_fin_,Winv_[0],naux_fin_);
    if (CholError !=0 )
        throw std::domain_error("J Matrix Cholesky Decomposition Failed!");
    
    //fprintf(outfile,"  Jchol\n");
    //print_mat(Winv_,naux_fin_,naux_fin_,outfile);
    
    //Inversion (in place)
    int IError = C_DPOTRI('L',naux_fin_,Winv_[0],naux_fin_);
    if (IError !=0 )
        throw std::domain_error("J Matrix Inversion Failed!");

    //LAPACK is smart and all, only uses half of the thing
    for (int m = 0; m<naux_fin_; m++)
        for (int n = 0; n<m; n++)
            Winv_[m][n] = Winv_[n][m]; 

    if (print_>4) {
        fprintf(outfile,"  J^-1\n");
        print_mat(Winv_,naux_fin_,naux_fin_,outfile);
    }

    timer_off("Form J^-1");

    timer_on("Overall (A|mn)");
    
    //Use ri_pair_mu_ and ri_pair_nu_ to keep track of things
    //Across schwarz sieve and unfortunate shell indexing
    ri_pair_nu_ = init_int_array(ntri_naive_);
    ri_pair_mu_ = init_int_array(ntri_naive_);
  
    //double three_index_cutoff = options_.get_double("THREE_INDEX_CUTOFF");
 
    if (df_storage_ == core)
    {
    	IntegralFactory rifactory(basisset_, basisset_, ribasis_, zero);
        shared_ptr<TwoBodyInt> eri = shared_ptr<TwoBodyInt>(rifactory.eri());
        const double *buffer = eri->buffer();
        A_ia_P_ = block_matrix(naux_fin_,ntri_naive_); 

        int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
        int start_index, delta_index, l_index;
        start_index = 0;
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU) {
                numnu = basisset_->shell(NU)->nfunction();
                if (schwarz_shell_pairs_[MU*(MU+1)/2+NU] == 1) {
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
                                if(omu>=onu && schwarz_fun_pairs_[omu*(omu+1)/2+onu] == 1) {
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
        if (print_>4) {
            fprintf(outfile,"  3-Index Tensor:\n");
            print_mat(A_ia_P_, naux_fin_,ntri_naive_ ,outfile);
            fprintf(outfile,"\n");
        }
        if (print_>5) {
            fprintf(outfile,"  Pair Indices:\n");
            for (int left = 0; left<ntri_; left++)
                fprintf(outfile,"  %d pair: (%d, %d)\n",left,ri_pair_mu_[left],ri_pair_nu_[left]);
            fprintf(outfile,"\n");
        }
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
                    //fprintf(outfile, "  MU = %d, NU = %d, Sig = %d\n",MU,NU,schwarz_shell_pairs_[MU*(MU+1)/2+NU]); fflush(outfile);
                    if (schwarz_shell_pairs_[MU*(MU+1)/2+NU] == 1) {
                        timer_on("(A|mn) Integrals");
                        eri->compute_shell(MU, NU, Pshell, 0);
                        timer_off("(A|mn) Integrals");
                        for (P=0; P < numP; ++P) {
                            PHI = ribasis_->shell(Pshell)->function_index() + P;
                            for (mu=0 ; mu < nummu; ++mu) {
                                omu = basisset_->shell(MU)->function_index() + mu;
                                for (nu=0; nu < numnu; ++nu) {
                                    onu = basisset_->shell(NU)->function_index() + nu;
                                    if(omu>=onu && schwarz_fun_pairs_[omu*(omu+1)/2+onu] == 1) {
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
        if (print_>5) {
            fprintf(outfile,"  Pair Indices:\n");
            for (int left = 0; left<ntri_; left++)
                fprintf(outfile,"  %d pair: (%d, %d)\n",left,ri_pair_mu_[left],ri_pair_nu_[left]);
            fprintf(outfile,"\n");
        }

        free(Temp1);
        //fprintf(outfile,"\n  Through B on disk."); fflush(outfile);
        psio_->close(PSIF_DFSCF_BJI,1);
    }
    timer_off("Overall (A|mn)");

    ri_back_map_ = init_int_array(norbs*(norbs+1)/2);
    for (int ij =0; ij<norbs*(norbs+1)/2; ij++)
        ri_back_map_[ij] = -1;
    for (int ij = 0 ; ij<ntri_; ij++)
        ri_back_map_[ri_pair_mu_[ij]*(ri_pair_mu_[ij]+1)/2+ri_pair_nu_[ij]] = ij;

    if (schwarz_) {
        fprintf(outfile,"\n  Function Pair Schwarz Sieve, Cutoff = %14.10E:\n",schwarz_);
        fprintf(outfile,"  %d out of %d basis function pairs removed, %8.5f%% attenuation.\n",norbs*(norbs+1)/2-sig_fun_pairs_,norbs*(norbs+1)/2,100.0*(norbs*(norbs+1)/2-sig_fun_pairs_)/(1.0*norbs*(norbs+1)/2));
        int pairs = basisset_->nshell()*(basisset_->nshell()+1)/2;
        fprintf(outfile,"  %d out of %d basis shell pairs removed, %8.5f%% attenuation.\n",pairs-sig_shell_pairs_,pairs,100.0*(pairs-sig_shell_pairs_)/(1.0*pairs));
    }
}
void HF::free_A()
{
    if (df_storage_ == core)
        free_block(A_ia_P_);
    free(ri_pair_mu_);
    free(ri_pair_nu_);
    free(ri_back_map_);
    free(schwarz_fun_pairs_);
    free(schwarz_shell_pairs_);
    free_block(W_);
    free_block(Winv_);
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
    
    //Rearrange the D matrix as a vector in terms of ri_pair indices
    //Off diagonal elements get 2x weight due to permutational symmetry
    double* DD;
    //Number of double occupied orbitals, LOCAL occupation matrix (ndocc x norbs, weird, but it won't matter)
    int ndocc = doccpi_[0];
    double** Cocc;
    
    if (J_is_required_) {
        DD = init_array(ntri_); 
        for (int ij = 0; ij<ntri_; ij++) {
            DD[ij] = D_->get(0,ri_pair_mu_[ij],ri_pair_nu_[ij]); 
            if (ri_pair_mu_[ij] != ri_pair_nu_[ij])
                DD[ij] *= 2.0;
                //only A irrep at the moment!!
        }
    }
    //The localization and domain ID code is the same for all storage types
    if (K_is_required_) {
        //Localize the canonical orbitals (or propagate guess)
        if ((iteration_-1)%options_.get_int("STEPS_PER_LOCALIZE") == 0) {
            timer_on("Pipek-Mekey");
            fully_localize_mos();      
            timer_off("Pipek-Mekey");
        } else {
            timer_on("Localization Update");
            propagate_local_mos();
            timer_off("Localization Update");
        }
        //Form the domain bookmarks
        timer_on("Lodwin Analysis");
        localized_Lodwin_charges(); 
        timer_off("Lodwin Analysis");
        //Form the domain bookmarks
        timer_on("Form Domains");
        form_domains(); 
        timer_off("Form Domains");
        
        //Get the L matrix (The local one!!!)
        Cocc = block_matrix(ndocc,norbs);
        for (int i=0; i<norbs; i++) {
            for (int j=0; j<ndocc; j++)
                Cocc[j][i] = L_->get(0,i,j);
            //only A irrep at the moment!!
        }
    }
    if (print_>2)
        L_->print(outfile);
    if (df_storage_ == core) {
        if (J_is_required_) {
            //Oh, and see R. Polly et. al., J. Chem. Phys. 102(21-22), pp. 2311-2321, DOI 10.1080/0026897042000274801
            double* c = init_array(naux_fin_);
            double* d = init_array(naux_fin_);
            double *J = init_array(ntri_);
            //DGEMV -> L:
            //c_A = (A|ls)*D_{ls}
            timer_on("J DDOT");
            C_DGEMV('N',naux_fin_,ntri_,1.0,A_ia_P_[0],ntri_naive_,DD,1,0.0,c,1);
            timer_off("J DDOT");

            //d_A = J^{-1}*c = J^{-1}_{AB}c_B
            timer_on("J Fitting");
            C_DGEMV('N',naux_fin_,naux_fin_,1.0,Winv_[0],naux_fin_,c,1,0.0,d,1);
            timer_off("J Fitting");

            //DGEMV -> J:
            //J_{mn} = d_B(B|mn)
            timer_on("J DAXPY");
            C_DGEMV('T',naux_fin_,ntri_,1.0,A_ia_P_[0],ntri_naive_,d,1,0.0,J,1);
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
            //Hold K as we add across i
            double **Ktemp = block_matrix(norbs,norbs);
            //Contribution to K from the current i
            double **K = block_matrix(max_domain_size_,max_domain_size_);
            //Half transformed 3-index tensor (by i)
            double **E = block_matrix(max_fit_size_,max_domain_size_);
            //Fitting coefficients for density
            double **D = block_matrix(max_fit_size_,max_domain_size_);
            //Local J and J^-1 matrix
            double **J = block_matrix(max_fit_size_,max_fit_size_);
            //QS is a temp matrix to allow for DGEMV into E by C_iv
            double **QS = block_matrix(max_fit_size_,max_domain_size_);
            //C_local is a local vector of significant n corresponding to a given m and i
            double *C_local = init_array(max_domain_size_); 

            //Main loop over i. Working in the MO basis is awesome
            for (int i = 0; i<ndocc; i++) {

                //fprintf(outfile,"  MO Orbital %d:\n",i);

                //Form J for this domain
                timer_on("Form Local J");
               for (int A = 0, Pl = 0; A<domain_atoms_[i]; A++)
                    for (int P = fit_fun_start_[i][A]; P<fit_fun_start_[i][A]+fit_fun_length_[i][A]; P++, Pl++)    
                        for (int B = 0, Ql = 0; B<domain_atoms_[i]; B++)
                            for (int Q = fit_fun_start_[i][B]; Q<fit_fun_start_[i][B]+fit_fun_length_[i][B]; Q++, Ql++)
                                 J[Pl][Ql] = W_[P][Q]; 
                timer_off("Form Local J");

                //fprintf(outfile,"  Local J:\n");
                //print_mat(J,fit_size_[i],fit_size_[i],outfile);
                //fprintf(outfile,"\n");
  
                //Invert J for this domain
                timer_on("Local J^-1");
                //Cholesky facotrization (in place)
                int CholError = C_DPOTRF('L',fit_size_[i],J[0],max_fit_size_);
                if (CholError !=0 )
                    throw std::domain_error("J Matrix Cholesky Decomposition Failed!");
    
                //Inversion (in place)
                int IError = C_DPOTRI('L',fit_size_[i],J[0],max_fit_size_);
                if (IError !=0 )
                    throw std::domain_error("J Matrix Inversion Failed!");

                //LAPACK is smart and all, only uses half of the thing
                for (int m = 0; m<fit_size_[i]; m++)
                    for (int n = 0; n<m; n++)
                        J[m][n] = J[n][m]; 
                //J now contains inverse local fitting metric 
                timer_off("Local J^-1");
                
                //fprintf(outfile,"  Local J^-1:\n");
                //print_mat(J,fit_size_[i],fit_size_[i],outfile);
                //fprintf(outfile,"\n");

                //Form E (Yeah, life's rough)
                timer_on("Local E");
                int counter, ij;
                for (int A = 0, ml = 0; A<domain_atoms_[i]; A++)
                    for (int m = domain_fun_start_[i][A]; m<domain_fun_start_[i][A]+domain_fun_length_[i][A]; m++, ml++) {
                        counter = 0;
                        for (int B =0, nl = 0; B<domain_atoms_[i]; B++)
                            for (int n = domain_fun_start_[i][B]; n<domain_fun_start_[i][B]+domain_fun_length_[i][B]; n++, nl++) {
                                if (n>m)
                                    ij = n*(n+1)/2+m;
                                else
                                    ij = m*(m+1)/2+n;
                                if (ri_back_map_[ij] >= 0) {

                                    C_local[counter] = Cocc[i][n];                                    
                                    for (int C = 0, Pl = 0; C<domain_atoms_[i]; C++)
                                        for (int P = fit_fun_start_[i][C]; P<fit_fun_start_[i][C]+fit_fun_length_[i][C]; P++, Pl++)
                                            QS[Pl][counter] = A_ia_P_[P][ri_back_map_[ij]];
                                    counter++;
                                }
                            }
                        //fprintf(outfile,"  m = %d:\n",m); 
                        //fprintf(outfile,"  C_local:\n");
                        //for (int dumm = 0; dumm<counter; dumm++)
                            //fprintf(outfile,"    C[%d] = %14.10f\n",dumm+1,C_local[dumm]);
                        //fprintf(outfile, "\n");

                        //fprintf(outfile,"  QS local:\n");
                        //print_mat(QS,fit_size_[i],domain_size_[i],outfile);
                        //fprintf(outfile,"\n");                       

                        //The number of significant n corresponding to this m is in counter
                        C_DGEMV('N',fit_size_[i],counter,1.0,QS[0],max_domain_size_,C_local,1,0.0,&E[0][ml],max_domain_size_);   
                    }
                timer_off("Local E");

                //fprintf(outfile,"  Local E:\n");
                //print_mat(E,fit_size_[i],domain_size_[i],outfile);
                //fprintf(outfile,"\n");

                //Form D
                timer_on("Form Local D");
                C_DGEMM('N','N',fit_size_[i],domain_size_[i],fit_size_[i],1.0,J[0],max_fit_size_,E[0],max_domain_size_,0.0,D[0],max_domain_size_);
                timer_off("Form Local D");

                //fprintf(outfile,"  Local D:\n");
                //print_mat(D,fit_size_[i],domain_size_[i],outfile);
                //fprintf(outfile,"\n");
                
                //Form K contributions
                timer_on("Local E DGEMM");
                C_DGEMM('T','N',domain_size_[i],domain_size_[i],fit_size_[i],1.0,E[0],max_domain_size_,D[0],max_domain_size_,0.0,K[0],max_domain_size_);
                timer_off("Local E DGEMM");

                //fprintf(outfile,"  Local K:\n");
                //print_mat(K,domain_size_[i],domain_size_[i],outfile);
                //fprintf(outfile,"\n");
                
                //Move the local K matrix back up to global scope
                for (int A=0, ml=0; A<domain_atoms_[i];A++)
                    for (int m=domain_fun_start_[i][A]; m<domain_fun_start_[i][A]+domain_fun_length_[i][A]; m++, ml++)  
                        for (int B=0, nl=0; B<domain_atoms_[i];B++)
                            for (int n=domain_fun_start_[i][B]; n<domain_fun_start_[i][B]+domain_fun_length_[i][B]; n++, nl++)
                                Ktemp[m][n] += K[ml][nl];  
            }
            //Load K_ and you're done
            for (int m = 0; m<norbs; m++)
                for (int n = 0; n<norbs; n++)
                    K_->set(0,m,n,Ktemp[m][n]);
        
            free(C_local);
            free_block(J);
            free_block(QS);
            free_block(E);
            free_block(D);
            free_block(K);
            free_block(Ktemp);
        }
    } else if (df_storage_ == disk) {
        //Not implemented yet
    } else {
        throw std::domain_error("Sepecified storage algorithm is not correct for L-DF-SCF");
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
    //At this point, I need a pint.
}
void RHF::propagate_local_mos()
{
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Must run in C1 for now.\n"); fflush(outfile);
        abort();
    } 
    if (print_ > 1)
        fprintf(outfile,"\n  Computing extrapolation of local orbitals.\n");
    
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

    //L_->print(outfile);
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

  //fprintf(outfile, "\C Matrix in the AO basis:\n");
  //print_mat(C, norbs, ndocc, outfile);

  //evals = get_evals();
  bool puream = basisset_->has_puream();
  snuc = init_int_array(nshell);
  for (i = 0; i<nshell; i++)
    snuc[i] = basisset_->shell(i)->ncenter();

  //for (int q = 0; q< natom; q++)
  //  fprintf(outfile,"  Snuc %d is %d",q,snuc[q]);
  //fflush(outfile);
  
  //fprintf(outfile, "Overlap Matrix");
  //print_mat(S, norbs, norbs, outfile);

  // Compute the length of each AM block
  int *l_length = init_int_array(basisset_->max_am()+1);
  l_length[0] = 1;
  for(l=1; l < (basisset_->max_am())+1; l++) {
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

  int nfun;
  int m = 0;
  for (int MU = 0; MU < nshell; MU++) {
   atom = snuc[MU] ; //Use c++ indexing !
   am = stype[MU];
   nfun = l_length[am];

   //fprintf(outfile, "  MU = %d, atom = %d, am = %d, nfun = %d, m = %d\n", MU, atom, am, nfun, m);
   
   if (aostart[atom] == 0 && atom != 0) {
    aostart_shell[atom] = MU;
    aostart[atom] = m;
   }
   m += nfun;
   aostop_shell[atom] = MU;
   aostop[atom] = m-1;
  }

  int* ao2atom = init_int_array(norbs);
  for(i=0; i < natom; i++)
    for(j=aostart[i]; j <= aostop[i]; j++) {
      ao2atom[j] = i;                   // ao2atom is the atom number that the AO is located on
    }

  //fprintf(outfile, "\tNumber of doubly occupied orbitals: %d\n\n", nocc);
  if (print_>2) {
    fprintf(outfile, "\n  Pipek-Mezey Localization Procedure:\n\n");

    fprintf(outfile, "\tIter     Pop. Localization   Max. Rotation Angle       Conv   Rotations\n");
    fprintf(outfile, "\t-----------------------------------------------------------------------\n");
  }
  V = block_matrix(nocc, nocc);
  double* tvec = init_array(nocc);
  double* svec = init_array(nocc);

  int max_iter = 200;
  double conv_tol = 1E-12;
  double angle_tol = 1E-9;
  int tot_rotations, rotations; 

  alphalast = 0.0;
  tot_rotations = 0;
  
  for(iter=0; iter < max_iter; iter++) {

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
  rotations = 0;

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

        if (fabs(alpha) > angle_tol) {
         rotations++;
         Uss = cos(alpha);                                            // Eqn 10a/b (JCP 90, 4916)
         Utt = Uss;                                               // Eqn 10a/b (JCP 90, 4916)
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
        }
      } // t-loop
    } // s-loop

    conv = fabs(alphamax) - fabs(alphalast);
    if (print_>2)
        fprintf(outfile, "\t%4d  %20.10f  %20.10f  %4.3e  %8d\n", iter, P, alphamax, conv, rotations);

    tot_rotations+=rotations;

    if((iter > 2) && ((fabs(conv) < conv_tol) || alphamax == 0.0)) break;
    alphalast = alphamax;

    fflush(outfile);

  } // iter-loop
  free(tvec);
  free(svec);
  free_block(V);

  if (print_>2) {
    fprintf(outfile,"\n  %d total rotations performed.\n",tot_rotations);
  }
  //fprintf(outfile, "\nC Matrix in the LO basis:\n");
  //print_mat(C, norbs, norbs, outfile);

  //free(evals);
  free(snuc);
  free(stype);
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
    free_block(S);
}
void RHF::localized_Lodwin_charges()
{
    //L_->print(outfile);

    //Compute Lodwin atomic charges of localized orbitals (for L-DF-SCF)
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Must run in C1 for now.\n"); fflush(outfile);
        abort();
    } 
    int ndocc = doccpi_[0];
    int norbs = basisset_->nbf();
    int natom = basisset_->molecule()->natom();
   
    //Shalf_->print(outfile); 
    //Sphalf_->print(outfile);

    double **Q = block_matrix(norbs,ndocc);
    {
        SharedMatrix temp(factory_.create_matrix("Temp"));
        //Q = Q_{mo} = S^{+1/2}C = S_{mn}^{+1/2}C_{mo}
        temp->gemm(false,false,1.0,Sphalf_,L_,0.0);
        for (int m = 0; m<norbs; m++)
            for (int o = 0; o<ndocc; o++)
                Q[m][o] = temp->get(0,m,o);
        //temp->print(outfile);
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
    
    if (print_>2) {
        fprintf(outfile,"  Lodwin Atomic Charges (atoms x occ):\n");
        print_mat(I_,natom,ndocc,outfile);
        fprintf(outfile,"\n");
    }

}
void RHF::form_domain_bookkeeping()
{
    //Sets up all the crazy tables for local domains
    //Must be called after form_A, as ribasis_ must be 
    //Initialized

    //Some sizes for clarity
    int ndocc = doccpi_[0];
    int norbs = basisset_->nbf();
    int natom = basisset_->molecule()->natom();
   
    //Number of contributing atoms per i (hopefully<<natom)
    domain_atoms_ = init_int_array(ndocc);

    //Array of primary domain shell starts per i
    domain_shell_start_ = (int**)malloc(ndocc*sizeof(int*));
    int* dummy_dss = init_int_array(ndocc*natom);
    for (int i = 0; i<ndocc; i++)
        domain_shell_start_[i] = &dummy_dss[i*natom];
 
    //Array of primary domain function starts per i
    domain_fun_start_ = (int**)malloc(ndocc*sizeof(int*));
    int* dummy_dfs = init_int_array(ndocc*natom);
    for (int i = 0; i<ndocc; i++)
        domain_fun_start_[i] = &dummy_dfs[i*natom];
 
    //Array of primary domain shell lengths per i
    domain_shell_length_ = (int**)malloc(ndocc*sizeof(int*));
    int* dummy_dsl = init_int_array(ndocc*natom);
    for (int i = 0; i<ndocc; i++)
        domain_shell_length_[i] = &dummy_dsl[i*natom];
 
    //Array of primary domain function lengths per i
    domain_fun_length_ = (int**)malloc(ndocc*sizeof(int*));
    int* dummy_dfl = init_int_array(ndocc*natom);
    for (int i = 0; i<ndocc; i++)
        domain_fun_length_[i] = &dummy_dfl[i*natom];
    
    //Array of primary fit shell starts per i
    fit_shell_start_ = (int**)malloc(ndocc*sizeof(int*));
    int* dummy_ass = init_int_array(ndocc*natom);
    for (int i = 0; i<ndocc; i++)
        fit_shell_start_[i] = &dummy_ass[i*natom]; //redundant no?
 
    //Array of primary fit function starts per i
    fit_fun_start_ = (int**)malloc(ndocc*sizeof(int*));
    int* dummy_afs = init_int_array(ndocc*natom);
    for (int i = 0; i<ndocc; i++)
        fit_fun_start_[i] = &dummy_afs[i*natom];
 
    //Array of primary fit shell lengths per i
    fit_shell_length_ = (int**)malloc(ndocc*sizeof(int*));
    int* dummy_asl = init_int_array(ndocc*natom);
    for (int i = 0; i<ndocc; i++)
        fit_shell_length_[i] = &dummy_asl[i*natom];
 
    //Array of primary fit function lengths per i
    fit_fun_length_ = (int**)malloc(ndocc*sizeof(int*));
    int* dummy_afl = init_int_array(ndocc*natom);
    for (int i = 0; i<ndocc; i++)
        fit_fun_length_[i] = &dummy_afl[i*natom];
 
    //Domain pair sizes (hopefully << norbs*(norbs+1)/2)
    domain_pairs_ = init_int_array(ndocc);
    //Biggest domain pair size, for memory allocation
    max_domain_pairs_ = 0;
    
    //Domain sizes (hopefully << norbs)
    domain_size_ = init_int_array(ndocc);
    //Biggest domain size, for memory allocation
    max_domain_size_ = 0;

    //Fitting domain sizes (hopefully << naux_fin_)
    fit_size_ = init_int_array(ndocc);
    //Biggest fitting domain size, for memory allocation
    max_fit_size_ = 0;

    //Has the atomic domain changed for this i?
    //If not, might be able to save a lot of work
    domain_changed_ = (bool*)malloc(ndocc*sizeof(bool));
    
    //Flag-style array for primary/extended domains;
    atom_domains_ = (int**)malloc(natom*sizeof(int*));
    int* dummy3 = init_int_array(ndocc*natom);
    for (int A = 0; A<natom; A++)
        atom_domains_[A] = &dummy3[A*ndocc];
    
    //Flag-style array for OLD primary/extended domains;
    old_atom_domains_ = (int**)malloc(natom*sizeof(int*));
    int* dummy4 = init_int_array(ndocc*natom);
    for (int A = 0; A<natom; A++)
        old_atom_domains_[A] = &dummy4[A*ndocc];
   
    //atom to primary basis shell start index;
    primary_shell_start_ = init_int_array(natom);
    //atom to primary function start index;
    primary_fun_start_ = init_int_array(natom);
    //atom to primary basis shell length;
    primary_shell_length_ = init_int_array(natom);
    //atom to primary function length;
    primary_fun_length_ = init_int_array(natom);
    int this_atom = 0;
    int this_shell_length = 0; 
    int this_fun_length = 0; 
    for (int MU = 0; MU< basisset_->nshell(); MU++) {
        if (basisset_->shell(MU)->ncenter() == this_atom+1) {
            primary_fun_length_[this_atom] = this_fun_length;
            primary_fun_start_[this_atom+1] = basisset_->shell(MU)->function_index();
            primary_shell_length_[this_atom] = this_shell_length;
            primary_shell_start_[this_atom+1] = MU;
            this_fun_length = basisset_->shell(MU)->nfunction();
            this_shell_length = 1;
            this_atom++;
        } else {
            this_shell_length++;
            this_fun_length += basisset_->shell(MU)->nfunction();
        }
    }
    primary_fun_length_[this_atom] = this_fun_length;
    primary_shell_length_[this_atom] = this_shell_length;

    //atom to aux basis shell start index;
    aux_shell_start_ = init_int_array(natom);
    //atom to aux function start index;
    aux_fun_start_ = init_int_array(natom);
    //atom to aux basis shell length;
    aux_shell_length_ = init_int_array(natom);
    //atom to aux function length;
    aux_fun_length_ = init_int_array(natom);
    this_atom = 0; 
    this_shell_length = 0; 
    this_fun_length = 0; 
    for (int MU = 0; MU< ribasis_->nshell(); MU++) {
        if (ribasis_->shell(MU)->ncenter() == this_atom+1) {
            aux_fun_length_[this_atom] = this_fun_length;
            aux_fun_start_[this_atom+1] = ribasis_->shell(MU)->function_index();
            aux_shell_length_[this_atom] = this_shell_length;
            aux_shell_start_[this_atom+1] = MU;
            this_fun_length = ribasis_->shell(MU)->nfunction();
            this_shell_length = 1;
            this_atom++;
        } else {
            this_shell_length++;
            this_fun_length += ribasis_->shell(MU)->nfunction();
        }
    }
    aux_fun_length_[this_atom] = this_fun_length;
    aux_shell_length_[this_atom] = this_shell_length;
} 
void RHF::form_domains()
{
    //Lodwin Atomic charge cutoff for primary domains
    double Q_cutoff = options_.get_double("CHARGE_CUTOFF");
    
    //Extended pair domain radius
    double R_ext_ang = options_.get_double("R_EXT");
    double R_ext = R_ext_ang*1.889725989; //1.889752989 [bohr/ang]

    //Some sizes for clarity
    int ndocc = doccpi_[0];
    int norbs = basisset_->nbf();
    int natom = basisset_->molecule()->natom();

    //ATOM DOMAIN IDENTIFICATION

    //PRIMARY DOMAIN 
    //For each i, set all elements of atom_domains_ to 1
    //if the Lodwin charge is big enough
    memset((void*)atom_domains_[0],'\0',ndocc*natom*sizeof(int));
    for (int i = 0; i<ndocc; i++)
        for (int N = 0 ; N<natom; N++)
            if (I_[N][i] >= Q_cutoff)
                atom_domains_[N][i] = 1; //Primary domain

    //EXTENDED DOMAIN 
    //For each i, set all elements of atom_domains_ to 2
    //if they are not already 1 AND are within R_ext bohr
    //of a primary domain atom
    for (int i = 0; i<ndocc; i++)
        for (int N = 0 ; N<natom; N++)
            if (I_[N][i] < Q_cutoff)
                for (int M = 0; M<natom; M++)
                    if (I_[M][i] >= Q_cutoff)
                        if ((basisset_->molecule()->xyz(M)).distance(basisset_->molecule()->xyz(N)) <= R_ext) {
                            atom_domains_[N][i] = 2; //Extended domain
                            continue;
                        }

    //Check for changes
    memset((void*)domain_changed_,'\0',ndocc*sizeof(bool));
    for (int i = 0; i<ndocc; i++)
        for (int N = 0; N<natom; N++) {
            if ((atom_domains_[N][i] > 0) && (old_atom_domains_[N][i] == 0)) {
                domain_changed_[i] = true; continue;
            } else if ((atom_domains_[N][i] == 0) && (old_atom_domains_[N][i] > 0)) {
                domain_changed_[i] = true; continue;
            }
        } 

    //Set old_atom_domains_ for next round
    memcpy((void*)old_atom_domains_[0],(void*)atom_domains_[0],ndocc*natom*sizeof(int));
    
    //BASIS DOMAIN IDENTIFICATION
    memset((void*)domain_atoms_,'\0',ndocc*sizeof(int));
    int counter;
    for (int i = 0; i<ndocc; i++) {
        counter = 0;
        for (int N = 0; N<natom; N++)
            if (atom_domains_[N][i] > 0) {
                domain_atoms_[i]++;
                domain_shell_start_[i][counter] = primary_shell_start_[N];
                domain_fun_start_[i][counter] = primary_fun_start_[N];
                domain_shell_length_[i][counter] = primary_shell_length_[N];
                domain_fun_length_[i][counter] = primary_fun_length_[N];
                fit_shell_start_[i][counter] = aux_shell_start_[N];
                fit_fun_start_[i][counter] = aux_fun_start_[N];
                fit_shell_length_[i][counter] = aux_shell_length_[N];
                fit_fun_length_[i][counter] = aux_fun_length_[N];
                counter++; 
            }
    }
    //Put -1 in non-used indices (will help with debugging)
    for (int i = 0; i<ndocc; i++) {
        for (int N=domain_atoms_[i]; N<natom; N++) {
            domain_shell_start_[i][N] = -1;
            domain_shell_length_[i][N] = -1;
            domain_fun_start_[i][N] = -1;
            domain_fun_length_[i][N] = -1;
            fit_shell_start_[i][N] = -1;
            fit_shell_length_[i][N] = -1;
            fit_fun_start_[i][N] = -1;
            fit_fun_length_[i][N] = -1;
        }
    }
            
    //PRIMARY BASIS DOMAIN SIZE IDENTIFICATION 
    memset((void*)domain_pairs_,'\0',ndocc*sizeof(int));
    memset((void*)domain_size_,'\0',ndocc*sizeof(int));
    memset((void*)fit_size_,'\0',ndocc*sizeof(int));
    max_domain_pairs_ = 0;
    max_domain_size_ = 0;
    max_fit_size_ = 0;
    for (int i = 0; i<ndocc; i++) {
        //Find domain pairs w/ schwarz sieve
        for (int A = 0; A<domain_atoms_[i]; A++)
            for (int m = domain_fun_start_[i][A]; m<domain_fun_start_[i][A]+domain_fun_length_[i][A]; m++)
                for (int B=0; B<=A; B++)
                    for (int n = domain_fun_start_[i][B]; n<domain_fun_start_[i][B]+domain_fun_length_[i][B] && n<=m; n++)
                        if (ri_back_map_[m*(m+1)/2+n] >= 0)
                            domain_pairs_[i]++; 

        //Domain size is just total number of primary basis functions
        //Same for fit_size_
        for (int A = 0; A<domain_atoms_[i]; A++) {
            domain_size_[i]+=domain_fun_length_[i][A];
            fit_size_[i]+=fit_fun_length_[i][A];
        }

        if (max_domain_pairs_ < domain_pairs_[i])
            max_domain_pairs_ = domain_pairs_[i];
        if (max_domain_size_ < domain_size_[i])
            max_domain_size_ = domain_size_[i];
        if (max_fit_size_ < fit_size_[i])
            max_fit_size_ = fit_size_[i];
    }

    if (print_>2) {
        fprintf(outfile,"  Atomic Domain Selection (atoms x occ):\n");
        print_int_mat(atom_domains_,natom,ndocc,outfile);
        fprintf(outfile,"\n");
        
        fprintf(outfile,"  Number of atoms in domain:\n");
        for (int i = 0; i<ndocc; i++)
            fprintf(outfile,"    Orbital %d, %d atoms\n",i+1,domain_atoms_[i]);
        fprintf(outfile,"\n");

    }

    if (print_>3) {
        fprintf(outfile,"  Domain changes?:\n");
        for (int i = 0; i<ndocc; i++)
            fprintf(outfile,"    Orbital %d, %s\n",i+1,(domain_changed_[i]?"Yes":"No"));
        fprintf(outfile,"\n");
    
        fprintf(outfile,"  Primary Basis Domain Shell Starts (occ x involved atoms):\n");
        print_int_mat(domain_shell_start_,ndocc,natom,outfile);
        fprintf(outfile,"\n");

        fprintf(outfile,"  Primary Basis Domain Shell Lengths (occ x involved atoms):\n");
        print_int_mat(domain_shell_length_,ndocc,natom,outfile);
        fprintf(outfile,"\n");

        fprintf(outfile,"  Primary Basis Domain Function Starts (occ x involved atoms):\n");
        print_int_mat(domain_fun_start_,ndocc,natom,outfile);
        fprintf(outfile,"\n");

        fprintf(outfile,"  Primary Basis Domain Function Lengths (occ x involved atoms):\n");
        print_int_mat(domain_fun_length_,ndocc,natom,outfile);
        fprintf(outfile,"\n");

        fprintf(outfile,"  Auxiliary Basis Domain Shell Starts (occ x involved atoms):\n");
        print_int_mat(fit_shell_start_,ndocc,natom,outfile);
        fprintf(outfile,"\n");

        fprintf(outfile,"  Auxiliary Basis Domain Shell Lengths (occ x involved atoms):\n");
        print_int_mat(fit_shell_length_,ndocc,natom,outfile);
        fprintf(outfile,"\n");

        fprintf(outfile,"  Auxiliary Basis Domain Function Starts (occ x involved atoms):\n");
        print_int_mat(fit_fun_start_,ndocc,natom,outfile);
        fprintf(outfile,"\n");

        fprintf(outfile,"  Auxiliary Basis Domain Function Lengths (occ x involved atoms):\n");
        print_int_mat(fit_fun_length_,ndocc,natom,outfile);
        fprintf(outfile,"\n");

        fprintf(outfile,"  Number of primary basis function pairs in domain:\n");
        for (int i = 0; i<ndocc; i++)
            fprintf(outfile,"    Orbital %d, %d function pairs\n",i+1,domain_pairs_[i]);
        fprintf(outfile,"\n");

        fprintf(outfile,"  Number of primary basis functions in domain:\n");
        for (int i = 0; i<ndocc; i++)
        fprintf(outfile,"    Orbital %d, %d functions\n",i+1,domain_size_[i]);
        fprintf(outfile,"\n");

        fprintf(outfile,"  Number of auxiliary basis functions in domain:\n");
        for (int i = 0; i<ndocc; i++)
        fprintf(outfile,"    Orbital %d, %d functions\n",i+1,fit_size_[i]);
        fprintf(outfile,"\n");
    }
    if (print_>2) {
        fprintf(outfile,"  Maximum number of primary basis function pairs for a domain is %d\n",max_domain_pairs_);
        fprintf(outfile,"  Maximum number of primary basis functions for a domain is %d\n",max_domain_size_);
        fprintf(outfile,"  Maximum number of auxiliary basis functions for a domain is %d\n",max_fit_size_);
    }
}
void RHF::free_domain_bookkeeping()
{
    //frees (int** pointers are tricky)
    free(domain_shell_start_[0]);
    free(domain_shell_start_);
    free(domain_shell_length_[0]);
    free(domain_shell_length_);
    free(domain_fun_start_[0]);
    free(domain_fun_start_);
    free(domain_fun_length_[0]);
    free(domain_fun_length_);
    free(fit_shell_start_[0]);
    free(fit_shell_start_);
    free(fit_shell_length_[0]);
    free(fit_shell_length_);
    free(fit_fun_start_[0]);
    free(fit_fun_start_);
    free(fit_fun_length_[0]);
    free(fit_fun_length_);
    
    free(primary_shell_start_);
    free(primary_shell_length_);
    free(primary_fun_start_);
    free(primary_fun_length_);
    free(aux_shell_start_);
    free(aux_shell_length_);
    free(aux_fun_start_);
    free(aux_fun_length_);

    free(old_atom_domains_[0]);
    free(old_atom_domains_);
    free(atom_domains_[0]);
    free(atom_domains_);
    
    free(domain_changed_);   
    free(domain_atoms_);
    free(domain_pairs_);
    free(domain_size_);
    free(fit_size_);
}
}}
