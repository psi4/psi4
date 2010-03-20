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
    fprintf(outfile, "\n  Computing Integrals using Density Fitting");
    if (factory_.nirreps() != 1)
    {
        fprintf(outfile,"Must run in C1 for now.\n"); fflush(outfile);
        abort();
    } 
    int norbs = basisset_->nbf(); 
    shared_ptr<BasisSet> ribasis_ =shared_ptr<BasisSet>(new BasisSet(chkpt_, "DF_BASIS"));
    ri_nbf_ = ribasis_->nbf();
    //ribasis_->print();
	/*
    fprintf(outfile, "\n  Memory Requirements:    (ab|P)    (ab|P)(PQ)^(-1/2)    Exchange Tensor    Max in Form B    Max in Form G");
    fprintf(outfile, "\n  --------------------------------------------------------------------------------------------------------"); */
    unsigned long memA = norbs*(norbs+1)/2*(long)ri_nbf_;
    int ndocc = doccpi_[0];
    unsigned long memC = norbs*ndocc*(long)ri_nbf_;
    //fprintf(outfile, "\n  Doubles:          %14ld %14ld      %14ld    %14ld %14ld",memA,memB,memC,memA+memB,memA+memC);
    //fprintf(outfile, "\n  MiB:               %14ld %14ld      %14ld    %14ld %14ld",memA*8/1000000,memB*8/1000000,memC*8/1000000,(memA+memB)*8/1000000,(memA+memC)*8/1000000);
    //fflush(outfile);

	string storage_type;
	storage_type = options_.get_str("RI_STORAGE");
	
	if (storage_type == "IN_CORE_DOUBLE")
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
        	df_storage_ = flip_B_disk;	//Transpose B using disk scratch and core, leave it on disk
    	else if (((long)((memC)*(1.0+MEMORY_SAFETY_FACTOR)))<(memory_/sizeof(double)))
       		df_storage_ = k_incore; //K only in-core
    	else
        	df_storage_ = disk;
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
    
    shared_ptr<BasisSet> zero = BasisSet::zero_basis_set();
    
    // Create integral factory
    IntegralFactory rifactory_J(ribasis_, zero, ribasis_, zero);
    TwoBodyInt* Jint = rifactory_J.eri();
    double **J = block_matrix(ri_nbf_, ri_nbf_);
    double **J_mhalf = block_matrix(ri_nbf_, ri_nbf_);
    const double *Jbuffer = Jint->buffer();
    
timer_on("Form J Matrix;");

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
    if (df_storage_ == double_full)
    {
    	IntegralFactory rifactory(basisset_, basisset_, ribasis_, zero);
        TwoBodyInt* eri = rifactory.eri();
        const double *buffer = eri->buffer();
        double** A_ia_P = block_matrix(ri_nbf_,basisset_->nbf()*(basisset_->nbf()+1)/2); 
        int numPshell,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu, npairs, start;
        
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU) {
                numnu = basisset_->shell(NU)->nfunction();
                for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
            		numPshell = ribasis_->shell(Pshell)->nfunction();
timer_on("(B|mn) Integrals");
                    eri->compute_shell(MU, NU, Pshell, 0);
timer_off("(B|mn) Integrals");
                    for (mu=0, index=0; mu < nummu; ++mu) {
                        omu = basisset_->shell(MU)->function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = basisset_->shell(NU)->function_index() + nu;
                            for (P=0; P < numPshell; ++P, ++index) {
                            	if(omu>=onu)
                            	{	
                        			PHI = ribasis_->shell(Pshell)->function_index() + P;
                                    A_ia_P[PHI][ioff[omu]+onu]= buffer[index];
                                }
                            } 
                        }
                    } // end loop over P in Pshell
                } // end loop over NU shell
            } // end loop over MU shell
            // now we've gone through all P, mu, nu for a given Pshell
        } // end loop over P shells; done with forming MO basis (P|ia)'s
        B_ia_P_ = block_matrix(ri_nbf_,basisset_->nbf()*(basisset_->nbf()+1)/2); 
timer_on("(B|mn) Transform");
        C_DGEMM('N','N',ri_nbf_,basisset_->nbf()*(basisset_->nbf()+1)/2,ri_nbf_,1.0, J_mhalf[0], ri_nbf_, A_ia_P[0], basisset_->nbf()*(basisset_->nbf()+1)/2,
             0.0, B_ia_P_[0], basisset_->nbf()*(basisset_->nbf()+1)/2);
timer_off("(B|mn) Transform");
        
	free(A_ia_P);
        //print_mat(B_ia_P_, ri_nbf_,norbs*(norbs+1)/2 ,outfile);
    } 
    if (df_storage_ == full||df_storage_ == flip_B_disk)
    {	
    	IntegralFactory rifactory(basisset_, basisset_, ribasis_, zero);
        TwoBodyInt* eri = rifactory.eri();
        const double *buffer = eri->buffer();
        B_ia_P_ = block_matrix(ri_nbf_,basisset_->nbf()*(basisset_->nbf()+1)/2); 
        int numPshell,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu, npairs, start;
        
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU) {
                numnu = basisset_->shell(NU)->nfunction();
                for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
            		numPshell = ribasis_->shell(Pshell)->nfunction();
timer_on("(B|mn) Integrals");
                    eri->compute_shell(MU, NU, Pshell, 0);
timer_off("(B|mn) Integrals");
                    for (mu=0, index=0; mu < nummu; ++mu) {
                        omu = basisset_->shell(MU)->function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = basisset_->shell(NU)->function_index() + nu;
                            for (P=0; P < numPshell; ++P, ++index) {
                            	if(omu>=onu)
                            	{	
                        			PHI = ribasis_->shell(Pshell)->function_index() + P;
                                    B_ia_P_[PHI][ioff[omu]+onu]= buffer[index];
                                }
                            } 
                        }
                    } // end loop over P in Pshell
                } // end loop over NU shell
            } // end loop over MU shell
            // now we've gone through all P, mu, nu for a given Pshell
        } // end loop over P shells; done with forming MO basis (P|ia)'s
	int start_index = 0;
        double **Temp1;
	double **Temp2;
	bool allocated = false;
	int unit = ri_nbf_; //May want to make this smaller for memory!!
	for (int index = 0; index<norbs*(norbs+1)/2; index+=ri_nbf_)
	{
		int cols = unit;
		if (index+unit>=norbs*(norbs+1)/2) {
			cols = norbs*(norbs+1)/2-index;
			if (allocated) {
				free(Temp1);
				free(Temp2);
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

        if (df_storage_ == flip_B_disk)
        {
        	write_B();
        	free(B_ia_P_);
        	ri_pair_nu_ = init_int_array(norbs*(norbs+1)/2);
        	ri_pair_mu_ = init_int_array(norbs*(norbs+1)/2);
        	for (int i=0, ij=0; i<norbs; i++)
        		for (int j=0; j<=i; j++, ij++)
        		{
        			ri_pair_mu_[ij] = i;
        			ri_pair_nu_[ij] = j;
        		}
        }
    	
        //print_mat(B_ia_P_, ri_nbf_,norbs*(norbs+1)/2 ,outfile);
    }
    else if (df_storage_ == k_incore||df_storage_ == disk)
    {
        psio_->open(PSIF_DFSCF_BJI,PSIO_OPEN_NEW);
        IntegralFactory rifactory(ribasis_, zero, basisset_, basisset_);
        TwoBodyInt* eri = rifactory.eri();
        const double *buffer = eri->buffer();
        double **Temp1;
	double *Temp2 = init_array(ri_nbf_);
	double *Temp3 = init_array(ri_nbf_);
	ri_pair_nu_ = init_int_array(basisset_->nbf()*(basisset_->nbf()+1)/2);
        ri_pair_mu_ = init_int_array(basisset_->nbf()*(basisset_->nbf()+1)/2);
        int row, npairs;

        psio_address next_PSIF_DFSCF_BJI = PSIO_ZERO;

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
                if (MU != 0 || NU !=0) {
		    free_block(Temp1);
		}
                Temp1 = block_matrix(ri_nbf_,npairs);
                //fprintf(outfile,"\n  Memory Allocated for Quartet (%d %d| P)",MU,NU); fflush(outfile);

                for (Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
                    numPshell = ribasis_->shell(Pshell)->nfunction();
timer_on("(B|mn) integrals");

                    eri->compute_shell(Pshell, 0, MU, NU);
timer_off("(B|mn) integrals");

                    for (P=0, index=0; P < numPshell; ++P) {
                        PHI = ribasis_->shell(Pshell)->function_index() + P;
                        row = 0;
                        for (mu=0; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                            for (nu=0; nu < numnu; ++nu, ++index) {
                                onu = basisset_->shell(NU)->function_index() + nu;
                                if(omu>=onu)
                                {
                                    Temp1[PHI][row++]= buffer[index];
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
	        for (int c = 0; c<npairs; c++) {	
timer_on("(B|mn) Transform");
			for (int r=0; r<ri_nbf_; r++)
				Temp2[r] = Temp1[r][c];	
	        	C_DGEMV('N',ri_nbf_,ri_nbf_,1.0, J_mhalf[0], ri_nbf_, Temp2,1,0.0, Temp3,1);
timer_off("(B|mn) Transform");

timer_on("(B|mn) disk");
               		psio_->write(PSIF_DFSCF_BJI,"BJ Three-Index Integrals",(char *) Temp3,sizeof(double)*ri_nbf_,next_PSIF_DFSCF_BJI,&next_PSIF_DFSCF_BJI);
timer_off("(B|mn) disk");
			
		}
		
                //fprintf(outfile,"\n  Finished Writing Quartet (%d %d| P)\n",MU,NU); fflush(outfile);
            }
        }

        free_block(Temp1);
        free(Temp2);
	free(Temp3);
        //fprintf(outfile,"\n  Through B on disk."); fflush(outfile);
        psio_->close(PSIF_DFSCF_BJI,1);
timer_on("(B|mn) restripe");
        
    	//RESTRIPE
        psio_->open(PSIF_DFSCF_BJI,PSIO_OPEN_OLD);
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
	next_PSIF_DFSCF_BJI = PSIO_ZERO;
	psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
	
	double *Temp = init_array(norbs*(norbs+1)/2);
	for (int Q = 0; Q < ri_nbf_; Q++)
		psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(Temp[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
	free(Temp);	

        int max_cols = (memory_/sizeof(double))/((1.0+MEMORY_SAFETY_FACTOR)*ri_nbf_);
        if (max_cols > norbs*(norbs+1)/2)
            max_cols = norbs*(norbs+1)/2;
	double **in_buffer = block_matrix(ri_nbf_,1);
        //max_cols = 100;
        double **buffer2 = block_matrix(ri_nbf_,max_cols);

        int buf_ind = 0;
        ULI global_offset = 0;
        for (int ij = 0; ij < norbs*(norbs+1)/2; ij++)
        {
            psio_->read(PSIF_DFSCF_BJI,"BJ Three-Index Integrals",(char *) &(in_buffer[0][0]),sizeof(double)*ri_nbf_,next_PSIF_DFSCF_BJI,&next_PSIF_DFSCF_BJI);
            //fprintf(outfile,"\n  Read in pair %d",ij); fflush(outfile);
            for (int Q = 0; Q<ri_nbf_; Q++)
            {
                buffer2[Q][buf_ind] = in_buffer[Q][0];
            }
            buf_ind++;
            //fprintf(outfile,"\n  Moved Pair to position %d in buffer",buf_ind); fflush(outfile);
            if (buf_ind == max_cols || ij == norbs*(norbs+1)/2-1)
            {
                //fprintf(outfile,"\n  Writing %d pairs to disk",buf_ind); fflush(outfile);
                for (int Q = 0; Q<ri_nbf_; Q++)
                {
                    //fprintf(outfile,"\n  Working on Q %d",Q); fflush(outfile);
                    next_PSIF_DFSCF_BJ = psio_get_address(PSIO_ZERO,(ULI)(Q*norbs*(ULI)(norbs+1)/2*sizeof(double)+global_offset*sizeof(double)));
                    //fprintf(outfile,"\n  Address Acquired"); fflush(outfile);
                    psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(buffer2[Q][0]),sizeof(double)*buf_ind,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
                    //fprintf(outfile,"\n  Entry Written"); fflush(outfile);
                }
                global_offset+=buf_ind;
                buf_ind=0;
            } 
        }

        free(in_buffer);
        free(buffer2);
        psio_->close(PSIF_DFSCF_BJI,0);
        psio_->close(PSIF_DFSCF_BJ,1);
timer_off("(B|mn) restripe");

        //fprintf(outfile,"\n  BJ Restriped on disk.\n"); fflush(outfile);
    }
timer_off("Overall (B|mn)");

}
void HF::write_B()
{
    int norbs = basisset_->nbf();
    psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_NEW);
    psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;
    for (int Q = 0; Q<ri_nbf_; Q++)
    {
        psio_->write(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(B_ia_P_[Q][0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
    }
    psio_->close(PSIF_DFSCF_BJ,1);
}
void HF::free_B()
{
	if (df_storage_ == full||df_storage_ == double_full)
    		free(B_ia_P_);
    else 
    {
    	free(ri_pair_mu_);
    	free(ri_pair_nu_);
    }
}
void RHF::form_G_from_RI()
{
timer_on("Overall G");
    //Get norbs
    int norbs = basisset_->nbf();     
    //Zero the J matrix
    J_->zero();
    //Zero the K matrix
    K_->zero();

    //D_->print(outfile);
    //C_->print(outfile);
    
    //Rearrange the D matix as a vector in terms of ri_pair indices
    //Off diagonal elements get 2x weight due to permutational symmetry
    double* DD = init_array(norbs*(norbs+1)/2);
    
    if (df_storage_ == full || df_storage_ == double_full) {
        for (int i = 0, ij = 0; i < norbs; i++)
            for (int j = 0; j<=i; j++, ij++) {
            DD[ij] = (i!=j?2.0:1.0)*D_->get(0,i,j); 
            //only A irrep at the moment!!
        }
    } else {
        for (int ij = 0; ij<norbs*(norbs+1)/2; ij++) {
            DD[ij] = D_->get(0,ri_pair_mu_[ij],ri_pair_nu_[ij]); 
            if (ri_pair_mu_[ij] != ri_pair_nu_[ij])
                DD[ij] *= 2.0;
                //only A irrep at the moment!!
        }
    }
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
            double *J = init_array(norbs*(norbs+1)/2);
            //DGEMV -> L:
            //L_Q = (Q|ls)*D_{ls}
            C_DGEMV('N',ri_nbf_,norbs*(norbs+1)/2,1.0,B_ia_P_[0],norbs*(norbs+1)/2,DD,1,0.0,L,1);
            //DGEMV -> J:
            //J_{mn} = L_Q(Q|mn)
            C_DGEMV('T',ri_nbf_,norbs*(norbs+1)/2,1.0,B_ia_P_[0],norbs*(norbs+1)/2,L,1,0.0,J,1);
            //Put everything in J_
            for (int i = 0, ij = 0; i < norbs; i++)
                for (int j = 0; j<=i; j++, ij++) {
                    J_->set(0,i,j,J[ij]);
                    J_->set(0,j,i,J[ij]);
                }

            free(L);
            free(J);
            //J_->print(outfile);
        }
        if (K_is_required_) {
            /* EXCHANGE PART */
            //E exchange matrix
            double** E = block_matrix(norbs, ndocc*ri_nbf_);
            //QS temp matrix for DGEMM
            double** QS = block_matrix(ri_nbf_,norbs);
            //Temp matrix for DGEMM
            double** Temp = block_matrix(ndocc,ri_nbf_);
            //Temporary K matrix
            double** K = block_matrix(norbs,norbs);

            //Exchange tensor E
            for (int m = 0; m<norbs; m++) {
                for (int Q = 0; Q<ri_nbf_; Q++) {
                    for (int s = 0; s<norbs; s++) {
                        QS[Q][s] = B_ia_P_[Q][((s>=m)?ioff[s]+m:ioff[m]+s)];
                    }
                }
                //E_{il}^Q = (Q|ln) C_{in}
                C_DGEMM('N','T',ndocc,ri_nbf_,norbs,1.0,Cocc[0],norbs,QS[0],norbs, 0.0, Temp[0], ri_nbf_);
                //print_mat(Temp,ndocc,ri_nbf_,outfile);
                int offset;
                for (int Q = 0; Q<ri_nbf_; Q++) {
                    offset = Q*ndocc;
                    for (int i = 0; i<ndocc; i++) {
                        E[m][i+offset] = Temp[i][Q];
                    }
                }
            }
            free_block(Temp);
            free_block(QS);

            //K_{mn} = E_{im}^QE_{in}^Q
            C_DGEMM('N','T',norbs,norbs,ri_nbf_*ndocc,1.0,E[0],ri_nbf_*ndocc,E[0],ri_nbf_*ndocc, 0.0, K[0], norbs);
            free_block(E);
            for (int i = 0; i < norbs; i++)
                for (int j = 0; j<=i; j++) {
                    K_->set(0,j,i,K[i][j]);
                    K_->set(0,i,j,K[i][j]);
            }
            free_block(K);
            //K_->print(outfile);
        }
    }
    else if (df_storage_ == flip_B_disk || df_storage_ == k_incore) {
        //B is on disk, E will be in core, Single disk pass
        //in_buffer stores single aux basis function rows of 
        //the three index tensor
        //TODO use larger blocks
        double *in_buffer = init_array(norbs*(norbs+1)/2);
        //Transformed integrals are stored in PSIF_DFSCF_BJ
        //Which better exist at this point
        psio_->open(PSIF_DFSCF_BJ,PSIO_OPEN_OLD);
        psio_address next_PSIF_DFSCF_BJ = PSIO_ZERO;

        //Coulomb convolution vector, done element by element
        double L;
        double *J;
        if (J_is_required_){
            J = init_array(norbs*(norbs+1)/2);
        }

        //Three index tensor, in core
        double **E;
        double **K;
        if (K_is_required_){
            E = block_matrix(norbs, ndocc*ri_nbf_);
        }

        int mu, nu;
        for (int Q = 0; Q<ri_nbf_; Q++)
        {
            //Read a single row of the (B|mn) tensor in, place in in_buffer
            psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
            
            if (J_is_required_) {
                /* COULOMB PART */
                //L_Q = (Q|ls)D_{ls}
                L = C_DDOT(norbs*(norbs+1)/2,in_buffer,1,DD,1);
                //J_{mn} += (Q|mn)L_Q
                C_DAXPY(norbs*(norbs+1)/2,L,in_buffer,1,J,1);
            } 
            if (K_is_required_) {
                /* EXCHANGE TENSOR */
                for (int ij = 0 ; ij<norbs*(norbs+1)/2; ij++)
                {
                    mu = ri_pair_mu_[ij];
                    nu = ri_pair_nu_[ij];
                    for (int i = 0; i<ndocc; i++)
                    {
                        E[mu][i+Q*ndocc]+=Cocc[i][nu]*in_buffer[ij];
                        if (mu != nu)
                            E[nu][i+Q*ndocc]+=Cocc[i][mu]*in_buffer[ij];
                    }
                }
            }
        }
        free(in_buffer);
        psio_->close(PSIF_DFSCF_BJ,1);

        /* Exchange Tensor DGEMM */
        if (K_is_required_) {
            K = block_matrix(norbs,norbs);
            C_DGEMM('N','T',norbs,norbs,ri_nbf_*ndocc,1.0,E[0],ri_nbf_*ndocc,E[0],ri_nbf_*ndocc, 0.0, K[0], norbs);
            free_block(E);
        }

        /* Form J and K */
        if (J_is_required_) {
            for (int i = 0, ij = 0; i < norbs; i++)
                for (int j = 0; j<=i; j++, ij++) {
                    J_->set(0,ri_pair_mu_[ij],ri_pair_nu_[ij],J[ij]);
                    J_->set(0,ri_pair_nu_[ij],ri_pair_mu_[ij],J[ij]);
                }
        }
        if (K_is_required_) {
            for (int i = 0; i < norbs; i++)
                for (int j = 0; j <=i; j++) {
                    K_->set(0,i,j,K[i][j]);
                    K_->set(0,j,i,K[i][j]);
                }
        }

        /* Frees */
        if (J_is_required_) {
            free(J);
        }
        if (K_is_required_) {
            free_block(K);
        }
    }
    else {
        //B is on disk, K will be in disk, Single disk pass
        //B is on disk, E will be in core, Single disk pass
        //in_buffer stores single aux basis function rows of 
        //the three ndex tensor
        //TODO use larger blocks
        double *in_buffer = init_array(norbs*(norbs+1)/2);
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
            J = init_array(norbs*(norbs+1)/2);
        }

        //Three index tensor, buffered to disk
        double *out_buffer;
        double **K;
        if (K_is_required_){
            out_buffer = init_array(norbs*ndocc);
        }

        int mu, nu;
        for (int Q = 0; Q<ri_nbf_; Q++)
        {
            //Read a single row of the (B|mn) tensor in, place in in_buffer
timer_on("Read B");
            psio_->read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
timer_off("Read B");
            
            if (J_is_required_) {
                /* COULOMB PART */
                //L_Q = (Q|ls)D_{ls}
timer_on("J DDOT");
                L = C_DDOT(norbs*(norbs+1)/2,in_buffer,1,DD,1);
timer_off("J DDOT");
                //J_{mn} += (Q|mn)L_Q
timer_on("J DAPXY");
                C_DAXPY(norbs*(norbs+1)/2,L,in_buffer,1,J,1);
timer_off("J DAPXY");
            } 
timer_on("Form E");
            if (K_is_required_) {
                /* EXCHANGE TENSOR */
                memset(out_buffer,0,norbs*ndocc*sizeof(double));
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
timer_on("Write E");
                psio_->write(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(out_buffer[0]),sizeof(double)*norbs*ndocc,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
timer_off("Write E");
            }
timer_off("Form E");
        }
        if (K_is_required_) {
            free(out_buffer);
            psio_->close(PSIF_DFSCF_K,1);
        }
        free(in_buffer);
        psio_->close(PSIF_DFSCF_BJ,1);
timer_on("Exchange DGEMM");

        /* Exchange Tensor DGEMM */
        if (K_is_required_) {
            psio_->open(PSIF_DFSCF_K,PSIO_OPEN_OLD);
            next_PSIF_DFSCF_K = PSIO_ZERO;

            K = block_matrix(norbs,norbs);
            in_buffer = init_array(norbs*ndocc);

            for (int Q = 0; Q<ri_nbf_; Q++) {
timer_on("E Read");
                psio_->read(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(in_buffer[0]),sizeof(double)*norbs*ndocc,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
timer_off("E Read");

                for (int m = 0; m<norbs; m++)
                    for (int n = 0; n<=m; n++)
                        for (int i = 0; i<ndocc; i++) {
                            K[m][n]+=in_buffer[m*ndocc+i]*in_buffer[n*ndocc+i];
                            K[n][m] = K[m][n];
                        }
            }
            free(in_buffer);
            psio_->close(PSIF_DFSCF_K,0);
        } 
timer_off("Exchange DGEMM");
        /* Form J and K */
        if (J_is_required_) {
            for (int i = 0, ij = 0; i < norbs; i++)
                for (int j = 0; j<=i; j++, ij++) {
                    J_->set(0,ri_pair_mu_[ij],ri_pair_nu_[ij],J[ij]);
                    J_->set(0,ri_pair_nu_[ij],ri_pair_mu_[ij],J[ij]);
                }
        }
        if (K_is_required_) {
            for (int i = 0; i < norbs; i++)
                for (int j = 0; j <=i; j++) {
                    K_->set(0,i,j,K[i][j]);
                    K_->set(0,j,i,K[i][j]);
                }
        }

        /* Frees */
        if (J_is_required_) {
            free(J);
        }
        if (K_is_required_) {
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
	int nbeta = nbetapi_[0];
	
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
    			int errcode = psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
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
    		int errcode = psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
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
    		int errcode = psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
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
    			
    		int errcode = psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
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
    		errcode = psio_write(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(out_buffer[0]),sizeof(double)*norbs*nalpha,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
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
    		int errcode = psio_read(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(in_buffer[0]),sizeof(double)*norbs*nalpha,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
    		
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
    		int errcode = psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
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
    			
    		int errcode = psio_read(PSIF_DFSCF_BJ,"BJ Three-Index Integrals",(char *) &(in_buffer[0]),sizeof(double)*norbs*(norbs+1)/2,next_PSIF_DFSCF_BJ,&next_PSIF_DFSCF_BJ);
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
    		errcode = psio_write(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(out_buffer[0]),sizeof(double)*norbs*nbeta_,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
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
    		int errcode = psio_read(PSIF_DFSCF_K,"Exchange Tensor",(char *) &(in_buffer[0]),sizeof(double)*norbs*nbeta_,next_PSIF_DFSCF_K,&next_PSIF_DFSCF_K);
    		
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
