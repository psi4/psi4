/* 
 *  LR_INTS.CC 
 *
 */

#ifdef HAVE_MKL
#include <mkl.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

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

#include "sapt.h"
#include "structs.h"

#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/molecule.h>
#include <libmints/gshell.h>
#include <libmints/gridblock.h>
#include <libmints/matrix.h>
#include <libscf_solver/integrator.h>
#include <libscf_solver/functionalfactory.h>

using namespace boost;
using namespace std;
using namespace psi;
using namespace scf;

namespace psi { namespace sapt {

void SAPT::lr_ints()
{
    //Open new LR Integrals file
    psio_->open(PSIF_SAPT_LRINTS,PSIO_OPEN_NEW);

    if (params_.print > 1)
        fprintf(outfile,"Computing numerical density-fitted integrals.\n");

    //Some parameters
    int norbs = calc_info_.nso;
    int naux = calc_info_.nri;
    int noccA = calc_info_.noccA;
    int noccB = calc_info_.noccB;
    int nvirA = calc_info_.nvirA;
    int nvirB = calc_info_.nvirB;
    int nmoA = noccA+nvirA;
    int nmoB = noccB+nvirB;
    double **CA = calc_info_.CA;
    double **CB = calc_info_.CB;
    double schwarz = params_.schwarz;
    double mem_safety = params_.mem_safety;

    //fprintf(outfile,"CA:\n");
    //print_mat(CA,norbs,nmoA,outfile);
    //fprintf(outfile,"CB:\n");
    //print_mat(CB,norbs,nmoB,outfile);

    double** DA_temp = block_matrix(norbs,norbs);
    C_DGEMM('N','T',norbs,norbs,noccA,1.0,&CA[0][0],nmoA,&CA[0][0],nmoA,0.0,&DA_temp[0][0],norbs);
    shared_ptr<Matrix> DA_ = shared_ptr<Matrix> (new Matrix(1,&norbs,&norbs));
    DA_->set(DA_temp);
    DA_->set_name("DA");
    //fprintf(outfile,"DA:\n");
    //print_mat(DA_temp,norbs,norbs,outfile);
    free_block(DA_temp);
    double** DB_temp = block_matrix(norbs,norbs);
    C_DGEMM('N','T',norbs,norbs,noccB,1.0,&CB[0][0],nmoB,&CB[0][0],nmoB,0.0,&DB_temp[0][0],norbs);
    shared_ptr<Matrix> DB_ = shared_ptr<Matrix> (new Matrix(1,&norbs,&norbs));
    DB_->set(DB_temp);
    DB_->set_name("DB");
    //fprintf(outfile,"DB:\n");
    //print_mat(DB_temp,norbs,norbs,outfile);
    free_block(DB_temp);

    DA_->print(outfile);
    DB_->print(outfile);

    //Prestripe
    psio_address lr_prestripe = PSIO_ZERO;
    double* bufferA = init_array(noccA*nvirA);
    for (int Q = 0; Q<naux; Q++)
        psio_->write(PSIF_SAPT_LRINTS,"A LR Integrals",(char *)&(bufferA[0]),sizeof(double)*noccA*nvirA,lr_prestripe,&lr_prestripe);
    free(bufferA);

    lr_prestripe = PSIO_ZERO;
    double* bufferB = init_array(noccB*nvirB);
    for (int Q = 0; Q<naux; Q++)
        psio_->write(PSIF_SAPT_LRINTS,"B LR Integrals",(char *)&(bufferB[0]),sizeof(double)*noccB*nvirB,lr_prestripe,&lr_prestripe);
    free(bufferB);

    //Memory considerations
    unsigned long int total_df_memory = (unsigned long int)(mem_safety*memory_/sizeof(double));
    unsigned long int block_memory = total_df_memory/2;
    unsigned long int grid_memory = total_df_memory-block_memory;
    unsigned long int memory_per_row = 0;

    memory_per_row += norbs*(ULI)norbs;
    memory_per_row += noccA*(ULI)norbs;
    memory_per_row += noccA*(ULI)nvirA;
    memory_per_row += noccB*(ULI)norbs;
    memory_per_row += noccB*(ULI)nvirB;

    int maxPshell = 0;
    for (int P = 0; P<ribasis_->nshell(); P++)
        if (maxPshell < ribasis_->shell(P)->nfunctions())
            maxPshell = ribasis_->shell(P)->nfunctions();


    int max_rows = block_memory/memory_per_row;
    if (max_rows > naux)
        max_rows = naux;
    if (max_rows < maxPshell)
        max_rows = maxPshell;

    int nblocks = 2*naux/max_rows;

    int* p_starts = init_int_array(nblocks);
    int* p_sizes = init_int_array(nblocks);
    int* block_starts = init_int_array(nblocks);
    int* block_stops = init_int_array(nblocks);
    int* block_sizes = init_int_array(nblocks);

    //Determine block sizes
    nblocks = 0;
    int fun_counter = 0;
    for (int A=0; A<ribasis_->nshell(); A++) {
        if (A == ribasis_->nshell() - 1) {
            block_sizes[nblocks]++;
            p_sizes[nblocks] += ribasis_->shell(A)->nfunction();
            block_stops[nblocks] = ribasis_->nshell();
            nblocks++;
            break;
        }

        if (fun_counter+ribasis_->shell(A)->nfunction() > max_rows) {
            block_stops[nblocks] = A;
            block_starts[nblocks+1] = A;
            block_sizes[nblocks+1] = 1;
            p_sizes[nblocks+1] = ribasis_->shell(A)->nfunction();
            p_starts[nblocks+1] = ribasis_->shell(A)->function_index();
            nblocks++;
            fun_counter = ribasis_->shell(A)->nfunction();
            continue;
        }
        p_sizes[nblocks] += ribasis_->shell(A)->nfunction();
        block_sizes[nblocks]++;
        fun_counter += ribasis_->shell(A)->nfunction();
    }

    //Block prints
    /**
    fprintf(outfile,"\n  Block info:\n");
    for (int i = 0; i<nblocks; i++) {
        fprintf(outfile,"  i = %d, block_starts = %d, block_stops = %d, block_sizes = %d, p_starts = %d, p_sizes = %d\n",
            i,block_starts[i],block_stops[i],block_sizes[i],p_starts[i],p_sizes[i]);
    }  **/ 
 
    //Blocks
    double **Amn = block_matrix(max_rows,norbs*(ULI)norbs);
    double **Ami = block_matrix(max_rows,norbs*(ULI)noccA);
    double **Aia = block_matrix(max_rows,noccA*(ULI)nvirA); 
    double **Bmi = block_matrix(max_rows,norbs*(ULI)noccB);
    double **Bia = block_matrix(max_rows,noccB*(ULI)nvirB); 

    //Threading considerations
    int nthread = 1;
    #ifdef _OPENMP
    if (options_.get_int("RI_INTS_NUM_THREADS") == 0)
        nthread = omp_get_max_threads();
    else
        nthread = options_.get_int("RI_INTS_NUM_THREADS");
    #endif
    int rank = 0;

    //Setup Integrator
    shared_ptr<Integrator> integrator = (shared_ptr<Integrator>)(Integrator::createIntegrator(molecule_,options_));
    if (params_.print > 1)
        fprintf(outfile,"%s\n",integrator->getString().c_str());
    //integrator->checkLebedev(outfile);
    //integrator->checkSphericalIntegrators(outfile);
    //integrator->checkMolecularIntegrators(outfile);
    //integrator->printAvailableLebedevGrids(outfile);
    //fflush(outfile);
 
    //Setup Functional
    int npoints = options_.get_int("N_BLOCK");
    FunctionalFactory fact;
    SharedFunctional lda_func = SharedFunctional(fact.getExchangeFunctional("X_LDA","NULL",npoints));
    SharedProperties primary_props = SharedProperties(Properties::constructProperties(basisset_,npoints));
    shared_ptr<BasisPoints> aux_props = (shared_ptr<BasisPoints>)(new BasisPoints(ribasis_,npoints));

    double density_check_A = 0.0;
    double density_check_B = 0.0;
    
    double *lda_values              =  lda_func->getValue();    
    double *lda_grads               =  lda_func->getGradientA();    
    const double *rho               =  primary_props->getDensity();    
    double **primary_points         =  primary_props->getPoints();    
    double **aux_points             =  aux_props->getPoints();  

    //Schwarz sieve
    timer_on("Schwarz");
    unsigned long int sig_shell_pairs = 0;
    int *schwarz_shell_pairs = init_int_array(basisset_->nshell()*(basisset_->nshell()+1)/2);
    if (schwarz > 0.0) {

        double* max_shell_val = init_array(basisset_->nshell()*(basisset_->nshell()+1)/2);;
        double max_global_val = 0.0;

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
                            if (max_global_val<abs(buffer[index]))
                                max_global_val = abs(buffer[index]);
                            if (max_shell_val[MUNU]<abs(buffer[index]))
                                max_shell_val[MUNU] = abs(buffer[index]);
                        }
                    }
                }
            }
        }

        for (int ij = 0; ij < basisset_->nshell()*(basisset_->nshell()+1)/2; ij ++)
            if (max_shell_val[ij]*max_global_val>=schwarz*schwarz){
                schwarz_shell_pairs[ij] = 1;
                sig_shell_pairs++;
            }

        free(max_shell_val);
    } else {
        for (int ij = 0; ij < basisset_->nshell()*(basisset_->nshell()+1)/2; ij++) {
            schwarz_shell_pairs[ij] = 1;
            sig_shell_pairs++;
        }
    }
    timer_off("Schwarz");
    
    //Setup analytical integrals
    //Get an ERI object for the AO three-index integrals 
    IntegralFactory rifactory(ribasis_,zero_, basisset_, basisset_);
    //Get a TEI for each thread
    const double **buffer = new const double*[nthread];
    shared_ptr<TwoBodyInt> *eri = new shared_ptr<TwoBodyInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri[Q] = shared_ptr<TwoBodyInt>(rifactory.eri());
        buffer[Q] = eri[Q]->buffer();
    }

    //File Handlers
    psio_address disk_address_A = PSIO_ZERO;
    psio_address disk_address_B = PSIO_ZERO;
    
    //indices for three-index integrals 
    int P, MU, NU, nump, nummu, numnu, p, mu, nu, op, omu, onu, index;
 
    //Integrate and transform
    for (int block = 0; block < nblocks; block++) {

        //>>>>>>>>>>>>>>>>>>>>>>>>//
        //   MONOMER A INTEGRALS  //
        //<<<<<<<<<<<<<<<<<<<<<<<<//        

        //ERIs
        //Zero that guy out!
        memset((void*)&Amn[0][0],'\0',p_sizes[block]*norbs*(ULI)norbs); 
    
        //Form Amn ints
        timer_on("(A|mn)");
        #pragma omp parallel for private (P, MU, NU, p, mu, nu, nump, nummu, numnu, op, omu, onu, index, rank) schedule (dynamic) num_threads(nthread)
        for (P=block_starts[block]; P < block_stops[block]; ++P) {
        #ifdef _OPENMP
           rank = omp_get_thread_num();
        #endif 
        nump = ribasis_->shell(P)->nfunction();
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU) {
                numnu = basisset_->shell(NU)->nfunction();
                if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                    eri[rank]->compute_shell(P, 0, MU, NU);
                    for (p=0, index=0; p < nump; ++p) {
                        op = ribasis_->shell(P)->function_index() + p;
                        for (mu=0; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                                for (nu=0; nu < numnu; ++nu, ++index) {
                                    onu = basisset_->shell(NU)->function_index() + nu;
                                    Amn[op-p_starts[block]][omu*norbs+onu] = buffer[rank][index]; // (op | omu onu) integral
                                    Amn[op-p_starts[block]][onu*norbs+omu] = buffer[rank][index]; // (op | onu omu) integral
                                }
                            }
                        } 
                    } 
                } 
            }
        }
        timer_off("(A|mn)");

        fprintf(outfile, "  Amn\n");
        print_mat(Amn,max_rows,norbs*norbs, outfile);
        
        //Numerical Integrals
        integrator->reset();
        while (!integrator -> isDone()) {
            shared_ptr<GridBlock> pts = integrator->getNextBlock();
            int true_pts = pts->getTruePoints();
            double* x = pts->getX(); 
            double* y = pts->getY(); 
            double* z = pts->getZ(); 
            double* w = pts->getWeights(); 

            primary_props->computeProperties(pts,DA_);
            aux_props->computePoints(pts);

            for (int pt = 0; pt < true_pts; pt++) {
             for (int P = 0; P < naux; P++) {
              for (int m = 0; m < norbs; m++) {
               for (int n = 0; n < norbs; n++)  {
                Amn[P][m*norbs+n] += primary_points[pt][m]*primary_points[pt][n]*
                    aux_points[pt][P]*w[pt]*lda_grads[pt]; 
               }
              }
             }
             density_check_A += w[pt]*rho[pt];
             //fprintf(outfile,"  <x,y,z,w,rho> = <%14.10f,%14.10f,%14.10f,%14.10f,%14.10f>\n",x[pt],y[pt],z[pt],w[pt],rho[pt]);   
             //fprintf(outfile,"   <phi,,z,w,rho> = <%14.10f,%14.10f,%14.10f,%14.10f,%14.10f>\n",x[pt],y[pt],z[pt],w[pt],rho[pt]);   
            }     
            
        }

        //Transform        
        //Transform to Ami
        // (A|mi) = (Amn)C_ni
        timer_on("(A|mi)");
        C_DGEMM('N', 'N', p_sizes[block]*norbs, noccA, norbs, 1.0, &(Amn[0][0]),
            norbs, &(CA[0][0]), nmoA, 0.0, &(Ami[0][0]), noccA);
        timer_off("(A|mi)");

        //fprintf(outfile, "  Ami\n");
        //print_mat(Ami,p_sizes[block],noccA*norbs, outfile);

        #ifdef HAVE_MKL
            int mkl_nthreads = mkl_get_max_threads();
            mkl_set_num_threads(1);
        #endif

        //Transform to Aia
        // (A|ia) = C_ma(A|mi)
        timer_on("(A|ia)");
        #pragma omp parallel for  
        for (int A = 0; A<p_sizes[block]; A++) {
            C_DGEMM('T', 'N', noccA, nvirA, norbs, 1.0, &(Ami[A][0]),
            noccA, &(CA[0][noccA]), nmoA, 0.0, &(Aia[A][0]), nvirA);
        }
        timer_off("(A|ia)");

        #ifdef HAVE_MKL
            mkl_set_num_threads(mkl_nthreads);
        #endif

        fprintf(outfile, "  Aia\n");
        print_mat(Aia,p_sizes[block],noccA*nvirA, outfile);

        //Stripe to disk
        timer_on("(A|ia) Write");
        psio_->write(PSIF_SAPT_LRINTS,"A LR Integrals",(char *)(&Aia[0][0]),sizeof(double)*p_sizes[block]*noccA*(ULI)nvirA,disk_address_A,&disk_address_A);
        timer_off("(A|ia) Write");

        //>>>>>>>>>>>>>>>>>>>>>>>>//
        //   MONOMER B INTEGRALS  //
        //<<<<<<<<<<<<<<<<<<<<<<<<//        

        //ERIs
        //Zero that guy out!
        memset((void*)&Amn[0][0],'\0',p_sizes[block]*norbs*(ULI)norbs); 
    
        //Form Amn ints
        timer_on("(A|mn)");
        #pragma omp parallel for private (P, MU, NU, p, mu, nu, nump, nummu, numnu, op, omu, onu, index, rank) schedule (dynamic) num_threads(nthread)
        for (P=block_starts[block]; P < block_stops[block]; ++P) {
        #ifdef _OPENMP
           rank = omp_get_thread_num();
        #endif 
        nump = ribasis_->shell(P)->nfunction();
        for (MU=0; MU < basisset_->nshell(); ++MU) {
            nummu = basisset_->shell(MU)->nfunction();
            for (NU=0; NU <= MU; ++NU) {
                numnu = basisset_->shell(NU)->nfunction();
                if (schwarz_shell_pairs[MU*(MU+1)/2+NU] == 1) {
                    eri[rank]->compute_shell(P, 0, MU, NU);
                    for (p=0, index=0; p < nump; ++p) {
                        op = ribasis_->shell(P)->function_index() + p;
                        for (mu=0; mu < nummu; ++mu) {
                            omu = basisset_->shell(MU)->function_index() + mu;
                                for (nu=0; nu < numnu; ++nu, ++index) {
                                    onu = basisset_->shell(NU)->function_index() + nu;
                                    Amn[op-p_starts[block]][omu*norbs+onu] = buffer[rank][index]; // (op | omu onu) integral
                                    Amn[op-p_starts[block]][onu*norbs+omu] = buffer[rank][index]; // (op | onu omu) integral
                                }
                            }
                        } 
                    } 
                } 
            }
        }
        timer_off("(A|mn)");

        fprintf(outfile, "  Bmn\n");
        print_mat(Amn,max_rows,norbs*norbs, outfile);
        
        //Numerical Integrals
        integrator->reset();
        while (!integrator -> isDone()) {
            shared_ptr<GridBlock> pts = integrator->getNextBlock();
            int true_pts = pts->getTruePoints();
            double* x = pts->getX(); 
            double* y = pts->getY(); 
            double* z = pts->getZ(); 
            double* w = pts->getWeights(); 

            primary_props->computeProperties(pts,DB_);
            aux_props->computePoints(pts);

            for (int pt = 0; pt < true_pts; pt++) {
             for (int P = 0; P < naux; P++) {
              for (int m = 0; m < norbs; m++) {
               for (int n = 0; n < norbs; n++)  {
                Amn[P][m*norbs+n] += primary_points[pt][m]*primary_points[pt][n]*
                    aux_points[pt][P]*w[pt]*lda_grads[pt]; 
               }
              }
             }
             density_check_B += w[pt]*rho[pt];
            }     
            
        }

        //Transform        
        //Transform to Ami
        // (A|mi) = (Amn)C_ni
        timer_on("(A|mi)");
        C_DGEMM('N', 'N', p_sizes[block]*norbs, noccB, norbs, 1.0, &(Amn[0][0]),
            norbs, &(CB[0][0]), nmoB, 0.0, &(Bmi[0][0]), noccB);
        timer_off("(A|mi)");

        //fprintf(outfile, "  Bmi\n");
        //print_mat(Bmi,p_sizes[block],noccB*norbs, outfile);

        #ifdef HAVE_MKL
            int mkl_nthreads = mkl_get_max_threads();
            mkl_set_num_threads(1);
        #endif

        //Transform to Aia
        // (A|ia) = C_ma(A|mi)
        timer_on("(A|ia)");
        #pragma omp parallel for  
        for (int A = 0; A<p_sizes[block]; A++) {
            C_DGEMM('T', 'N', noccB, nvirB, norbs, 1.0, &(Bmi[A][0]),
            noccB,&(CB[0][noccB]), nmoB, 0.0, &(Bia[A][0]), nvirB);
        }
        timer_off("(A|ia)");

        #ifdef HAVE_MKL
            mkl_set_num_threads(mkl_nthreads);
        #endif

        fprintf(outfile, "  Bia\n");
        print_mat(Bia,p_sizes[block],noccB*nvirB, outfile);

        //Stripe to disk
        timer_on("(A|ia) Write");
        psio_->write(PSIF_SAPT_LRINTS,"B LR Integrals",(char *)(&Bia[0][0]),sizeof(double)*p_sizes[block]*noccB*(ULI)nvirB,disk_address_B,&disk_address_B);
        timer_off("(A|ia) Write");
    }

    fprintf(outfile,"Density check for monomer A: %14.10f electrons found, should be %d.\n", density_check_A,noccA*2); 
    fprintf(outfile,"Density check for monomer B: %14.10f electrons found, should be %d.\n", density_check_B,noccB*2); 
    
    //Frees
    free_block(Amn);
    free_block(Ami);
    free_block(Bmi);
    free_block(Aia);
    free_block(Bia);
    free(schwarz_shell_pairs);
    free(block_starts);
    free(block_sizes);
    free(block_stops);
    free(p_starts);
    free(p_sizes);

    //Close LR Integrals file
    psio_->close(PSIF_SAPT_LRINTS,1);    

}

}}
