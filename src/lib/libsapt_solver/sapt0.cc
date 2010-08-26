/* 
 *  SAPT0.CC 
 *
 */

#ifdef HAVEMKL
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
#include <time.h>

#include <psifiles.h>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include <libmints/mints.h>

#include "sapt0.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace sapt {

SAPT0::SAPT0(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : SAPT(options, psio, chkpt)
{
      std::string user = Process::environment("USER");
      if (user == "sherrill") {
        printf("I'm sorry, Dave. I'm afraid I can't do that.\n");
        abort();
      }
}

SAPT0::~SAPT0()
{
}
double SAPT0::compute_energy()
{
    if (!params_.df_restart) {
      ao_df_ints();
      lr_ints();
    } else {
      ao_df_ints_restart();
    }
    
    elst10();
    exch10();

    psio_->open(PSIF_SAPT_AMPS,PSIO_OPEN_NEW);
    disp20();
    theta_ar();
    theta_bs();
    exch_disp20();
    psio_->close(PSIF_SAPT_AMPS,0);

    ind20resp();
    exch_ind20respA_B();
    exch_ind20respB_A();

    print_results();

    double energy = results_.sapt0;
    Process::environment.globals["SAPT ENERGY"] = energy;

    return energy;
}
void SAPT0::ao_df_ints()
{
  psio_->open(PSIF_SAPT_AA_DF_INTS,PSIO_OPEN_NEW);
  psio_->open(PSIF_SAPT_BB_DF_INTS,PSIO_OPEN_NEW);
  psio_->open(PSIF_SAPT_AB_DF_INTS,PSIO_OPEN_NEW);

  if (params_.logfile) {
    fprintf(params_.logfilename," AO DF Integrals\n");
    fprintf(params_.logfilename,"-----------------\n\n");
    fflush(params_.logfilename);
  }
  time_t start;
  time_t stop;

  // Create a new matrix factory
  MatrixFactory factory;

  //TODO MAke this go away
  // Initialize the factory with data from checkpoint
  factory.init_with_chkpt(chkpt_);

  if (params_.logfile) {
    fprintf(params_.logfilename,"  Building J matrix\n");
    fflush(params_.logfilename);
  }

  timer_on("Form J matrix        ");

  // Create integral factory
  IntegralFactory rifactory_J(ribasis_, zero_, ribasis_, zero_);

  TwoBodyInt* Jint = rifactory_J.eri();
  double **J = block_matrix(ribasis_->nbf(), ribasis_->nbf());
  double **J_mhalf = block_matrix(ribasis_->nbf(),ribasis_->nbf());
  const double *Jbuffer = Jint->buffer();

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
        }
      }
    }
  }

  if (params_.logfile) {
    fprintf(params_.logfilename,"  Forming J^-1/2 matrix\n");
    fflush(params_.logfilename);
  }

  double* eigval = init_array(ribasis_->nbf());
  int lwork = ribasis_->nbf() * 3;
  double* work = init_array(lwork);
  int stat = C_DSYEV('v','u',ribasis_->nbf(),J[0],ribasis_->nbf(),eigval,
      work,lwork);
  if (stat != 0) {
    fprintf(outfile, "C_DSYEV failed\n");
    exit(PSI_RETURN_FAILURE);
  }
  free(work);

  // Now J contains the eigenvectors of the original J
  // Copy J to J_copy
  double **J_copy = block_matrix(ribasis_->nbf(), ribasis_->nbf());
  C_DCOPY(ribasis_->nbf()*ribasis_->nbf(),J[0],1,J_copy[0],1); 

  // Now form J^{-1/2} = U(T)*j^{-1/2}*U,
  // where j^{-1/2} is the diagonal matrix of the inverse square roots
  // of the eigenvalues, and U is the matrix of eigenvectors of J
#pragma omp for schedule(dynamic)
  for (int i=0; i<ribasis_->nbf(); i++) {
    if (eigval[i] < 1.0E-10)
      eigval[i] = 0.0;
    else {
      eigval[i] = 1.0 / sqrt(eigval[i]);
    }
    // scale one set of eigenvectors by the diagonal elements j^{-1/2}
    C_DSCAL(ribasis_->nbf(), eigval[i], J[i], 1);
  }
  free(eigval);

  // J_mhalf = J_copy(T) * J
  C_DGEMM('t','n',ribasis_->nbf(),ribasis_->nbf(),ribasis_->nbf(),1.0,
    J_copy[0],ribasis_->nbf(),J[0],ribasis_->nbf(),0.0,J_mhalf[0],
    ribasis_->nbf());

  free_block(J);
  free_block(J_copy);

  timer_off("Form J matrix        ");

  int mem_usage;
  int norbs = calc_info_.nso;
  int nmo = calc_info_.nmo;
  int ntri = calc_info_.nso*(calc_info_.nso+1)/2;
  int nshell = basisset_->nshell();
  int nshelltri = nshell*(nshell+1)/2;
  int nriorbs = ribasis_->nbf();
  double mem1;
  double avail_mem = params_.memory;
  calc_info_.nrio = nriorbs + 3;

  // Get Schwartz screening arrays

  if (params_.logfile) {
    fprintf(params_.logfilename,"  Forming Schwarz screening arrays\n");
    fflush(params_.logfilename);
  }

  timer_on("Form screening arrays");

  IntegralFactory ao_eri_factory(basisset_, basisset_, basisset_, basisset_);
  TwoBodyInt* ao_eri = ao_eri_factory.eri();
  const double *ao_buffer = ao_eri->buffer();

  double *Schwartz = init_array(basisset_->nshell() * (basisset_->nshell()+1) / 2);
  double *DFSchwartz = init_array(ribasis_->nshell());

  for(int P=0,PQ=0;P<basisset_->nshell();P++) {
    int numw = basisset_->shell(P)->nfunction();
    for(int Q=0;Q<=P;Q++,PQ++) {
      int numx = basisset_->shell(Q)->nfunction();
      double tei, max=0.0;

      ao_eri->compute_shell(P, Q, P, Q);

      for(int w=0;w<numw;w++) {
        for(int x=0;x<numx;x++) {
          index = ( ( (w*numx + x) * numw + w) * numx + x);
          tei = ao_buffer[index];
          if(fabs(tei) > max) max = fabs(tei);
        }
      }
      Schwartz[PQ] = max;
    }
  }

  for(int P=0;P<ribasis_->nshell();P++) {
    int numw = ribasis_->shell(P)->nfunction();
    double tei, max=0.0;

    Jint->compute_shell(P, 0, P, 0);

    for(int w=0;w<numw;w++) {
      tei = Jbuffer[w];
      if(fabs(tei) > max) max = fabs(tei);
    }
    DFSchwartz[P] = max;
  }

  timer_off("Form screening arrays");

  // find out the max number of P's in a P shell
  IntegralFactory rifactory(ribasis_, zero_, basisset_, basisset_);
  TwoBodyInt** eri;
  const double **buffer;

  int nthreads = 1;
  #ifdef _OPENMP
    nthreads = omp_get_max_threads();
  #endif
  int rank = 0;

  eri = new TwoBodyInt*[nthreads];
  buffer = new const double*[nthreads];
  for(int i = 0;i < nthreads;++i){
    eri[i] = rifactory.eri();
    buffer[i] = eri[i]->buffer();
  }

  int mn;
  int maxPshell = 0;
  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell)->nfunction();
    if (numPshell > maxPshell) maxPshell = numPshell;
  }

  psio_->open(PSIF_SAPT_TEMP,0);

  double** AO_RI = block_matrix(maxPshell,norbs*norbs);
  double* halftrans = init_array(nmo*norbs);
  double** MO_RI = block_matrix(maxPshell,nmo*nmo);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n  Starting AA 3-Center Formation\n\n");
    fflush(params_.logfilename);
  }

  psio_address next_DF_MO = PSIO_ZERO;

  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell)->nfunction();

    if (params_.logfile) {
      if (Pshell%100==0) {
        fprintf(params_.logfilename,"    Starting Shell %5d ... ",Pshell);
        start = time(NULL);
        fflush(params_.logfilename);
      }
    }

    #pragma omp for private(mn) schedule(guided)
    for (mn=0; mn < nshelltri; ++mn) {
      #ifdef _OPENMP
        rank = omp_get_thread_num();
      #endif
      int MU = calc_info_.index2i[mn];
      int NU = calc_info_.index2j[mn];  
      int nummu = basisset_->shell(MU)->nfunction();
      int numnu = basisset_->shell(NU)->nfunction();

      if (sqrt(Schwartz[mn]*DFSchwartz[Pshell])>params_.schwarz) {
        eri[rank]->compute_shell(Pshell, 0, MU, NU);
  
        for (int P=0, index=0; P < numPshell; ++P) {

          for (int mu=0; mu < nummu; ++mu) {
            int omu = basisset_->shell(MU)->function_index() + mu;
  
            for (int nu=0; nu < numnu; ++nu, ++index) {
              int onu = basisset_->shell(NU)->function_index() + nu;

              AO_RI[P][omu*norbs+onu] = buffer[rank][index];
              AO_RI[P][onu*norbs+omu] = buffer[rank][index];
            }
          } 
        }
    
      } // end Schwartz inequality
    } // end loop over MU,NU shells

    for (int P=0; P < numPshell; ++P) {
      C_DGEMM('T', 'N', nmo, norbs, norbs, 1.0, calc_info_.CA[0], nmo, 
        AO_RI[P], norbs, 0.0, halftrans, norbs);
      C_DGEMM('N', 'N', nmo, nmo, norbs, 1.0, halftrans, norbs, 
        calc_info_.CA[0], nmo, 0.0, MO_RI[P], nmo);
    }

    psio_->write(PSIF_SAPT_TEMP,"MO RI Integrals",(char *) &(MO_RI[0][0]),
      sizeof(double)*numPshell*nmo*(ULI) nmo,next_DF_MO,&next_DF_MO);

    if (params_.logfile) {
      if (Pshell%100==99 || Pshell+1 == ribasis_->nshell()) {
        stop = time(NULL);
        fprintf(params_.logfilename,"%5d finished %14ld seconds\n",Pshell,stop-start);
        fflush(params_.logfilename); 
      }
    }

  }

  free_block(AO_RI);
  free_block(MO_RI);

  double* zeros = init_array(nriorbs+3);

  timer_on("Zero AA Ints         ");
  zero_disk(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.noccA*calc_info_.noccA);
  zero_disk(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.nvirA*calc_info_.noccA);
  zero_disk(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.nvirA*calc_info_.nvirA);
  timer_off("Zero AA Ints         ");

  int numP;
  int temp_size = (int) ((avail_mem) / (16.0*((double) nriorbs)));

  if (temp_size > nmo*nmo)
    temp_size = nmo*nmo;

  double** temp;
  double** temp_J;

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n  Starting AA Multiplication by J^-1/2\n\n");
    fprintf(params_.logfilename,"     Block     Length\n");
    fprintf(params_.logfilename,"    -------  ----------\n");
    for (int i=0,oP=0; oP < nmo*nmo; ++i, oP += numP) {
      if ((i+1)*temp_size < nmo*nmo) numP = temp_size;
      else numP = nmo*nmo - oP;
      fprintf(params_.logfilename,"       %3d   %10d\n",i,numP);
    }
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  psio_address next_DF_AA = PSIO_ZERO;
  psio_address next_DF_AR = PSIO_ZERO;
  psio_address next_DF_RR = PSIO_ZERO;

  for (int i=0,oP=0; oP < nmo*nmo; ++i, oP += numP) {
    if ((i+1)*temp_size < nmo*nmo) numP = temp_size;
    else numP = nmo*nmo - oP;

    if (params_.logfile) {
      fprintf(params_.logfilename,"    Starting Block %3d ... ",i);
      start = time(NULL);
      fflush(params_.logfilename);
    }

    temp = block_matrix(nriorbs,numP);

    next_DF_MO = psio_get_address(PSIO_ZERO,oP*(ULI) sizeof(double));
    for (int P=0; P < nriorbs; ++P) {
      psio_->read(PSIF_SAPT_TEMP,"MO RI Integrals",(char *) &(temp[P][0]),
        sizeof(double)*(ULI) numP,next_DF_MO,&next_DF_MO);
      next_DF_MO = psio_get_address(next_DF_MO,(nmo*nmo-numP)*
        (ULI) sizeof(double));
    } 

    temp_J = block_matrix(numP,calc_info_.nrio);

    C_DGEMM('T','N',numP,nriorbs,nriorbs,1.0,temp[0],numP,
      J_mhalf[0],nriorbs,0.0,temp_J[0],calc_info_.nrio);

    free_block(temp);

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmo;
      int j = (ij+oP)%nmo;
      if (i < calc_info_.noccA && j < calc_info_.noccA) {
        psio_->write(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ULI) calc_info_.nrio,next_DF_AA,
          &next_DF_AA);
      }
    }

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmo;
      int j = (ij+oP)%nmo;
      if (i < calc_info_.noccA && j >= calc_info_.noccA) {
        psio_->write(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ULI) calc_info_.nrio,next_DF_AR,
          &next_DF_AR);
      }
    }

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmo;
      int j = (ij+oP)%nmo;
      if (i >= calc_info_.noccA && j >= calc_info_.noccA) {
        psio_->write(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ULI) calc_info_.nrio,next_DF_RR,
          &next_DF_RR);
      }
    }

    free_block(temp_J);

    if (params_.logfile) {
      stop = time(NULL);
      fprintf(params_.logfilename,"finished %14ld seconds\n",stop-start);
      fflush(params_.logfilename);
    }
  }

  AO_RI = block_matrix(maxPshell,norbs*norbs);
  MO_RI = block_matrix(maxPshell,nmo*nmo);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n  Starting BB 3-Center Formation\n\n");
    fflush(params_.logfilename);
  }

  next_DF_MO = PSIO_ZERO;

  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell)->nfunction();

    if (params_.logfile) {
      if (Pshell%100==0) {
        fprintf(params_.logfilename,"    Starting Shell %5d ... ",Pshell);
        start = time(NULL);
        fflush(params_.logfilename);
      }
    }

    #pragma omp for private(mn) schedule(guided)
    for (mn=0; mn < nshelltri; ++mn) {
      #ifdef _OPENMP
        rank = omp_get_thread_num();
      #endif
      int MU = calc_info_.index2i[mn];
      int NU = calc_info_.index2j[mn];
      int nummu = basisset_->shell(MU)->nfunction();
      int numnu = basisset_->shell(NU)->nfunction();

      if (sqrt(Schwartz[mn]*DFSchwartz[Pshell])>params_.schwarz) {
        eri[rank]->compute_shell(Pshell, 0, MU, NU);

        for (int P=0, index=0; P < numPshell; ++P) {

          for (int mu=0; mu < nummu; ++mu) {
            int omu = basisset_->shell(MU)->function_index() + mu;

            for (int nu=0; nu < numnu; ++nu, ++index) {
              int onu = basisset_->shell(NU)->function_index() + nu;

              AO_RI[P][omu*norbs+onu] = buffer[rank][index];
              AO_RI[P][onu*norbs+omu] = buffer[rank][index];
            }
          }
        }

      } // end Schwartz inequality
    } // end loop over MU,NU shells

    for (int P=0; P < numPshell; ++P) {
      C_DGEMM('T', 'N', nmo, norbs, norbs, 1.0, calc_info_.CB[0], nmo,
        AO_RI[P], norbs, 0.0, halftrans, norbs);
      C_DGEMM('N', 'N', nmo, nmo, norbs, 1.0, halftrans, norbs,
        calc_info_.CB[0], nmo, 0.0, MO_RI[P], nmo);
    }

    psio_->write(PSIF_SAPT_TEMP,"MO RI Integrals",(char *) &(MO_RI[0][0]),
      sizeof(double)*numPshell*nmo*(ULI) nmo,next_DF_MO,&next_DF_MO);

    if (params_.logfile) {
      if (Pshell%100==99 || Pshell+1 == ribasis_->nshell()) {
        stop = time(NULL);
        fprintf(params_.logfilename,"%5d finished %14ld seconds\n",Pshell,stop-start);
        fflush(params_.logfilename);
      }
    }
  }

  free_block(AO_RI);
  free_block(MO_RI);

  timer_on("Zero BB Ints         ");
  zero_disk(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.noccB*calc_info_.noccB);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.nvirB*calc_info_.noccB);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.nvirB*calc_info_.nvirB);
  timer_off("Zero BB Ints         ");

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n  Starting BB Multiplication by J^-1/2\n\n");
    fprintf(params_.logfilename,"     Block     Length\n");
    fprintf(params_.logfilename,"    -------  ----------\n");
    for (int i=0,oP=0; oP < nmo*nmo; ++i, oP += numP) {
      if ((i+1)*temp_size < nmo*nmo) numP = temp_size;
      else numP = nmo*nmo - oP;
      fprintf(params_.logfilename,"       %3d   %10d\n",i,numP);
    }
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  psio_address next_DF_BB = PSIO_ZERO;
  psio_address next_DF_BS = PSIO_ZERO;
  psio_address next_DF_SS = PSIO_ZERO;

  for (int i=0,oP=0; oP < nmo*nmo; ++i, oP += numP) {
    if ((i+1)*temp_size < nmo*nmo) numP = temp_size;
    else numP = nmo*nmo - oP;

    if (params_.logfile) {
      fprintf(params_.logfilename,"    Starting Block %3d ... ",i);
      start = time(NULL);
      fflush(params_.logfilename);
    }

    temp = block_matrix(nriorbs,numP);

    next_DF_MO = psio_get_address(PSIO_ZERO,oP*(ULI) sizeof(double));
    for (int P=0; P < nriorbs; ++P) {
      psio_->read(PSIF_SAPT_TEMP,"MO RI Integrals",(char *) &(temp[P][0]),
        sizeof(double)*(ULI) numP,next_DF_MO,&next_DF_MO);
      next_DF_MO = psio_get_address(next_DF_MO,(nmo*nmo-numP)*
        (ULI) sizeof(double));
    }

    temp_J = block_matrix(numP,calc_info_.nrio);

    C_DGEMM('T','N',numP,nriorbs,nriorbs,1.0,temp[0],numP,
      J_mhalf[0],nriorbs,0.0,temp_J[0],calc_info_.nrio);

    free_block(temp);

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmo;
      int j = (ij+oP)%nmo;
      if (i < calc_info_.noccB && j < calc_info_.noccB) {
        psio_->write(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ULI) calc_info_.nrio,next_DF_BB,
          &next_DF_BB);
      }
    }

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmo;
      int j = (ij+oP)%nmo;
      if (i < calc_info_.noccB && j >= calc_info_.noccB) {
        psio_->write(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ULI) calc_info_.nrio,next_DF_BS,
          &next_DF_BS);
      }
    }

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmo;
      int j = (ij+oP)%nmo;
      if (i >= calc_info_.noccB && j >= calc_info_.noccB) {
        psio_->write(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ULI) calc_info_.nrio,next_DF_SS,
          &next_DF_SS);
      }
    }

    free_block(temp_J);

    if (params_.logfile) {
      stop = time(NULL);
      fprintf(params_.logfilename,"finished %14ld seconds\n",stop-start);
      fflush(params_.logfilename);
    }
  }

  AO_RI = block_matrix(maxPshell,norbs*norbs);
  MO_RI = block_matrix(maxPshell,nmo*nmo);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n  Starting AB 3-Center Formation\n\n");
    fflush(params_.logfilename);
  }

  next_DF_MO = PSIO_ZERO; 

  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell)->nfunction();
 
    if (params_.logfile) {
      if (Pshell%100==0) {
        fprintf(params_.logfilename,"    Starting Shell %5d ... ",Pshell);
        start = time(NULL);
        fflush(params_.logfilename);
      }
    }

    #pragma omp for private(mn) schedule(guided)
    for (mn=0; mn < nshelltri; ++mn) {
      #ifdef _OPENMP
        rank = omp_get_thread_num();
      #endif
      int MU = calc_info_.index2i[mn];
      int NU = calc_info_.index2j[mn];
      int nummu = basisset_->shell(MU)->nfunction();
      int numnu = basisset_->shell(NU)->nfunction();
  
      if (sqrt(Schwartz[mn]*DFSchwartz[Pshell])>params_.schwarz) {
        eri[rank]->compute_shell(Pshell, 0, MU, NU);

        for (int P=0, index=0; P < numPshell; ++P) {
  
          for (int mu=0; mu < nummu; ++mu) {
            int omu = basisset_->shell(MU)->function_index() + mu;
  
            for (int nu=0; nu < numnu; ++nu, ++index) {
              int onu = basisset_->shell(NU)->function_index() + nu;

              AO_RI[P][omu*norbs+onu] = buffer[rank][index];
              AO_RI[P][onu*norbs+omu] = buffer[rank][index];
            }
          }
        }
    
      } // end Schwartz inequality
    } // end loop over MU,NU shells
  
    for (int P=0; P < numPshell; ++P) {
      C_DGEMM('T', 'N', nmo, norbs, norbs, 1.0, calc_info_.CA[0], nmo,
        AO_RI[P], norbs, 0.0, halftrans, norbs);
      C_DGEMM('N', 'N', nmo, nmo, norbs, 1.0, halftrans, norbs,
        calc_info_.CB[0], nmo, 0.0, MO_RI[P], nmo);
    }
 
    psio_->write(PSIF_SAPT_TEMP,"MO RI Integrals",(char *) &(MO_RI[0][0]),
      sizeof(double)*numPshell*nmo*(ULI) nmo,next_DF_MO,&next_DF_MO);
 
    if (params_.logfile) {
      if (Pshell%100==99 || Pshell+1 == ribasis_->nshell()) {
        stop = time(NULL);
        fprintf(params_.logfilename,"%5d finished %14ld seconds\n",Pshell,stop-start);
        fflush(params_.logfilename);
      }
    }
  }
  
  free_block(AO_RI); 
  free_block(MO_RI);

  timer_on("Zero AB Ints         ");
  zero_disk(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.noccA*calc_info_.noccB);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.nvirB*calc_info_.noccA);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"RB RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.noccB*calc_info_.nvirA);
  timer_off("Zero AB Ints         ");

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n  Starting AB Multiplication by J^-1/2\n\n");
    fprintf(params_.logfilename,"     Block     Length\n");
    fprintf(params_.logfilename,"    -------  ----------\n");
    for (int i=0,oP=0; oP < nmo*nmo; ++i, oP += numP) {
      if ((i+1)*temp_size < nmo*nmo) numP = temp_size;
      else numP = nmo*nmo - oP;
      fprintf(params_.logfilename,"       %3d   %10d\n",i,numP);
    }   
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  } 

  psio_address next_DF_AB = PSIO_ZERO;
  psio_address next_DF_AS = PSIO_ZERO;
  psio_address next_DF_RB = PSIO_ZERO;

  for (int i=0,oP=0; oP < nmo*nmo; ++i, oP += numP) {
    if ((i+1)*temp_size < nmo*nmo) numP = temp_size;
    else numP = nmo*nmo - oP;
      
    if (params_.logfile) {
      fprintf(params_.logfilename,"    Starting Block %3d ... ",i);
      start = time(NULL);
      fflush(params_.logfilename);
    }

    temp = block_matrix(nriorbs,numP);

    next_DF_MO = psio_get_address(PSIO_ZERO,oP*(ULI) sizeof(double));
    for (int P=0; P < nriorbs; ++P) {
      psio_->read(PSIF_SAPT_TEMP,"MO RI Integrals",(char *) &(temp[P][0]),
        sizeof(double)*(ULI) numP,next_DF_MO,&next_DF_MO);
      next_DF_MO = psio_get_address(next_DF_MO,(nmo*nmo-numP)*
        (ULI) sizeof(double));
    }

    temp_J = block_matrix(numP,calc_info_.nrio);

    C_DGEMM('T','N',numP,nriorbs,nriorbs,1.0,temp[0],numP,
      J_mhalf[0],nriorbs,0.0,temp_J[0],calc_info_.nrio);

    free_block(temp);

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmo;
      int j = (ij+oP)%nmo;
      if (i < calc_info_.noccA && j < calc_info_.noccB) {
        psio_->write(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ULI) calc_info_.nrio,next_DF_AB,
          &next_DF_AB);
      }
    }

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmo;
      int j = (ij+oP)%nmo;
      if (i < calc_info_.noccA && j >= calc_info_.noccB) {
        psio_->write(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ULI) calc_info_.nrio,next_DF_AS,
          &next_DF_AS);
      }
    }

    for (int ij=0; ij<numP; ij++) {
      int i = (ij+oP)/nmo;
      int j = (ij+oP)%nmo;
      if (i >= calc_info_.noccA && j < calc_info_.noccB) {
        psio_->write(PSIF_SAPT_AB_DF_INTS,"RB RI Integrals",(char *)
          &(temp_J[ij][0]),sizeof(double)*(ULI) calc_info_.nrio,next_DF_RB,
          &next_DF_RB);
      }
    }

    free_block(temp_J);

    if (params_.logfile) {
      stop = time(NULL);
      fprintf(params_.logfilename,"finished %14ld seconds\n",stop-start);
      fflush(params_.logfilename);
    }
  }

  free(zeros);
  free(halftrans);

  psio_->close(PSIF_SAPT_TEMP,0);

  free(Schwartz);
  free(DFSchwartz);
  free_block(J_mhalf);

  delete []eri;

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  fflush(outfile);
     
}
void SAPT0::ao_df_ints_restart()
{
  psio_->open(PSIF_SAPT_AA_DF_INTS,PSIO_OPEN_OLD);
  psio_->open(PSIF_SAPT_BB_DF_INTS,PSIO_OPEN_OLD);
  psio_->open(PSIF_SAPT_AB_DF_INTS,PSIO_OPEN_OLD);
    
  if (params_.print)
    fprintf(outfile,"  DF AO Integral Restart\n\n");
    
  fflush(outfile); 
}
double** SAPT0::get_AA_ints(int dress) {

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.noccA*calc_info_.noccA,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.noccA*calc_info_.noccA*
    (ULI) calc_info_.nrio);

  if (dress) {
    for (int a=0; a<calc_info_.noccA; a++){
      int aa = a*calc_info_.noccA+a;
      A[aa][calc_info_.nrio-3] = 1.0;
      A[aa][calc_info_.nrio-1] = enuc;
      for (int ap=0; ap<calc_info_.noccA; ap++){
        int aap = a*calc_info_.noccA+ap;
        A[aap][calc_info_.nrio-2] = NB*calc_info_.VBAA[a][ap];
      }
    }
  }

  return(A);

}

double** SAPT0::get_diag_AA_ints(int dress) {

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.noccA,calc_info_.nrio);

  psio_address next_PSIF = PSIO_ZERO;
  for (int a=0; a<calc_info_.noccA; a++){
    psio_->read(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *)
      &(A[a][0]),sizeof(double)*(ULI) calc_info_.nrio,next_PSIF,&next_PSIF);
    next_PSIF = psio_get_address(next_PSIF,calc_info_.noccA*calc_info_.nrio*
      (ULI) sizeof(double));
    if (dress) {
      A[a][calc_info_.nrio-3] = 1.0;
			A[a][calc_info_.nrio-2] = NB*calc_info_.VBAA[a][a];
      A[a][calc_info_.nrio-1] = enuc;
    }
  }

  return(A);
}

double** SAPT0::get_BB_ints(int dress) {

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.noccB*calc_info_.noccB,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.noccB*calc_info_.noccB*
    (ULI) calc_info_.nrio);

  if (dress) {
    for (int b=0; b<calc_info_.noccB; b++){
      int bb = b*calc_info_.noccB+b;
      A[bb][calc_info_.nrio-2] = 1.0;
      A[bb][calc_info_.nrio-1] = enuc;
      for (int bp=0; bp<calc_info_.noccB; bp++){
        int bbp = b*calc_info_.noccB+bp;
        A[bbp][calc_info_.nrio-3] = NA*calc_info_.VABB[b][bp];
      }
    }
  }

  return(A);

}

double** SAPT0::get_diag_BB_ints(int dress) {

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.noccB,calc_info_.nrio);

  psio_address next_PSIF = PSIO_ZERO;
  for (int b=0; b<calc_info_.noccB; b++){
    psio_->read(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *)
      &(A[b][0]),sizeof(double)*(ULI) calc_info_.nrio,next_PSIF,&next_PSIF);
    next_PSIF = psio_get_address(next_PSIF,calc_info_.noccB*calc_info_.nrio*
      (ULI) sizeof(double));
    if (dress) {
      A[b][calc_info_.nrio-3] = NA*calc_info_.VABB[b][b];
      A[b][calc_info_.nrio-2] = 1.0;
      A[b][calc_info_.nrio-1] = enuc;
    }
  }

  return(A);
}

double** SAPT0::get_AB_ints(int dress) {

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.noccA*calc_info_.noccB*
    (ULI) calc_info_.nrio);

  if (dress==1) {
    for (int a=0; a<calc_info_.noccA; a++){
      for (int b=0; b<calc_info_.noccB; b++){
        int ab = a*calc_info_.noccB+b;
        A[ab][calc_info_.nrio-3] = calc_info_.S_AB[a][b];
        A[ab][calc_info_.nrio-2] = NB*calc_info_.VBAB[a][b];
        A[ab][calc_info_.nrio-1] = enuc*calc_info_.S_AB[a][b];
      }
    }
  }
  else if (dress==2) {
    for (int a=0; a<calc_info_.noccA; a++){
      for (int b=0; b<calc_info_.noccB; b++){
        int ab = a*calc_info_.noccB+b;
        A[ab][calc_info_.nrio-3] = NA*calc_info_.VAAB[a][b];
        A[ab][calc_info_.nrio-2] = calc_info_.S_AB[a][b];
        A[ab][calc_info_.nrio-1] = enuc*calc_info_.S_AB[a][b];
      }
    }
  }

  return(A);

}

double** SAPT0::get_AS_ints(int dress) {

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.noccA*calc_info_.nvirB,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.noccA*calc_info_.nvirB*
    (ULI) calc_info_.nrio);

  if (dress) {
    for (int a=0; a<calc_info_.noccA; a++){
      for (int s=0; s<calc_info_.nvirB; s++){
        int as = a*calc_info_.nvirB+s;
        A[as][calc_info_.nrio-3] = calc_info_.S_AB[a][s+calc_info_.noccB];
        A[as][calc_info_.nrio-2] = NB*calc_info_.VBAB[a][s+calc_info_.noccB];
        A[as][calc_info_.nrio-1] = enuc*calc_info_.S_AB[a][s+calc_info_.noccB];
      }
    }
  }

  return(A);

}

double** SAPT0::get_RB_ints(int dress) {

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.nvirA*calc_info_.noccB,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_AB_DF_INTS,"RB RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.nvirA*calc_info_.noccB*
    (ULI) calc_info_.nrio);

  if (dress) {
    for (int r=0; r<calc_info_.nvirA; r++){
      for (int b=0; b<calc_info_.noccB; b++){
        int rb = r*calc_info_.noccB+b;
        A[rb][calc_info_.nrio-3] = NA*calc_info_.VAAB[r+calc_info_.noccA][b];
        A[rb][calc_info_.nrio-2] = calc_info_.S_AB[r+calc_info_.noccA][b];
        A[rb][calc_info_.nrio-1] = enuc*calc_info_.S_AB[r+calc_info_.noccA][b];
      }
    }
  }

  return(A);

}

double** SAPT0::get_AR_ints(int dress) {

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.nvirA*calc_info_.noccA,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.noccA*calc_info_.nvirA*
    (ULI) calc_info_.nrio);

  if (dress) {
    for (int a=0; a<calc_info_.noccA; a++){
      for (int r=0; r<calc_info_.nvirA; r++){
        int ar = a*calc_info_.nvirA+r;
        A[ar][calc_info_.nrio-2] = NB*calc_info_.VBAA[a][r+calc_info_.noccA];
      }
    }
  }

  return(A);

}

double** SAPT0::get_BS_ints(int dress) {

  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **A = block_matrix(calc_info_.nvirB*calc_info_.noccB,calc_info_.nrio);

  psio_->read_entry(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",(char *)
    &(A[0][0]),sizeof(double)*calc_info_.noccB*calc_info_.nvirB*
    (ULI) calc_info_.nrio);

  if (dress) {
    for (int b=0; b<calc_info_.noccB; b++){
      for (int s=0; s<calc_info_.nvirB; s++){
        int bs = b*calc_info_.nvirB+s;
        A[bs][calc_info_.nrio-3] = NA*calc_info_.VABB[b][s+calc_info_.noccB];
      }
    }
  }

  return(A);

}
void SAPT0::print_results()
{
  results_.exch_indr20 = results_.exch_indrA_B + results_.exch_indrB_A;

  results_.sapt0 = results_.hf_int + results_.disp20 + results_.exch_disp20;

  results_.deltaHF = results_.hf_int - (results_.elst10 + results_.exch10 +
                    results_.indr20 + results_.exch_indr20);

  fprintf(outfile,"    SAPT Results  \n");
  fprintf(outfile,"  ------------------------------------------------------------------\n");
  fprintf(outfile,"    E_HF          %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.hf_int*1000.0,results_.hf_int*627.5095);
  fprintf(outfile,"    Elst10        %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.elst10*1000.0,results_.elst10*627.5095);
  fprintf(outfile,"    Exch10(S^2)   %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.exch10*1000.0,results_.exch10*627.5095);
  fprintf(outfile,"    Ind20,r       %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.indr20*1000.0,results_.indr20*627.5095);
  fprintf(outfile,"    Exch-Ind20,r  %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.exch_indr20*1000.0,results_.exch_indr20*627.5095);
  fprintf(outfile,"    delta HF,r    %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.deltaHF*1000.0,results_.deltaHF*627.5095);
  fprintf(outfile,"    Disp20        %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.disp20*1000.0,results_.disp20*627.5095);
  fprintf(outfile,"    Exch-Disp20   %16.8lf mH %16.8lf kcal mol^-1\n\n",
          results_.exch_disp20*1000.0,results_.exch_disp20*627.5095);
  fprintf(outfile,"    Total SAPT0   %16.8lf mH %16.8lf kcal mol^-1\n",
          results_.sapt0*1000.0,results_.sapt0*627.5095);
}


}}
