/* This handles the 1e, 2e and DF integrals for all SAPT jobs */

#ifdef _MKL
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
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include <libmints/mints.h>

#include "sapt2b.h"
#include "sapt.h"

namespace psi { namespace sapt {

void SAPT::compute_integrals()
{
  oetrans();
  df_ints();
  w_ints();
}

void SAPT2B::df_ints()
{
  psio_->open(PSIF_SAPT_AA_DF_INTS,PSIO_OPEN_NEW);
  psio_->open(PSIF_SAPT_BB_DF_INTS,PSIO_OPEN_NEW);
  psio_->open(PSIF_SAPT_AB_DF_INTS,PSIO_OPEN_NEW);

  // Create integral factory
  IntegralFactory rifactory_J(ribasis_, zero_, ribasis_, zero_);

  TwoBodyAOInt* Jint = rifactory_J.eri();
//double **calc_info_.J = block_matrix(ribasis_->nbf(), ribasis_->nbf());
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

//C_DCOPY(ribasis_->nbf()*ribasis_->nbf(),&(J[0][0]),1,
//  &(calc_info_.J[0][0]),1);

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

  int mem_usage;
  int norbs = calc_info_.nso;
  int nmo = calc_info_.nmo;
  int ntri = calc_info_.nso*(calc_info_.nso+1)/2;
  int nshell = basisset_->nshell();
  int nshelltri = nshell*(nshell+1)/2;
  int nriorbs = ribasis_->nbf();
  calc_info_.nrio = nriorbs + 3;

  // Get Schwartz screening arrays

  IntegralFactory ao_eri_factory(basisset_, basisset_, basisset_, basisset_);
  TwoBodyAOInt* ao_eri = ao_eri_factory.eri();
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

  // find out the max number of P's in a P shell
  IntegralFactory rifactory(ribasis_, zero_, basisset_, basisset_);
  TwoBodyAOInt** eri;
  const double **buffer;

  int nthreads = 1;
  #ifdef _OPENMP
    nthreads = omp_get_max_threads();
  #endif
  int rank = 0;

  eri = new TwoBodyAOInt*[nthreads];
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

  psio_address next_DF_MO = PSIO_ZERO;
  psio_address next_bare_AR = PSIO_ZERO;

  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell)->nfunction();

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

  }

  free_block(AO_RI);
  free_block(MO_RI);

  double* zeros = init_array(nriorbs+3);

  zero_disk(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.noccA*calc_info_.noccA);
  zero_disk(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.nvirA*calc_info_.noccA);
  zero_disk(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.nvirA*calc_info_.nvirA);

  int numP;
  long int temp_size = params_.memory / (2*sizeof(double)*(long int) nriorbs);

  if (temp_size > nmo*nmo)
    temp_size = nmo*nmo;

  double** temp;
  double** temp_J;

  psio_address next_DF_AA = PSIO_ZERO;
  psio_address next_DF_AR = PSIO_ZERO;
  psio_address next_DF_RR = PSIO_ZERO;

  for (int i=0,oP=0; oP < nmo*nmo; ++i, oP += numP) {
    if ((i+1)*temp_size < nmo*nmo) numP = temp_size;
    else numP = nmo*nmo - oP;

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
  }

  AO_RI = block_matrix(maxPshell,norbs*norbs);
  MO_RI = block_matrix(maxPshell,nmo*nmo);

  next_DF_MO = PSIO_ZERO;
  psio_address next_bare_BS = PSIO_ZERO;

  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell)->nfunction();

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

  }

  free_block(AO_RI);
  free_block(MO_RI);

  zero_disk(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.noccB*calc_info_.noccB);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.nvirB*calc_info_.noccB);
  zero_disk(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.nvirB*calc_info_.nvirB);

  psio_address next_DF_BB = PSIO_ZERO;
  psio_address next_DF_BS = PSIO_ZERO;
  psio_address next_DF_SS = PSIO_ZERO;

  for (int i=0,oP=0; oP < nmo*nmo; ++i, oP += numP) {
    if ((i+1)*temp_size < nmo*nmo) numP = temp_size;
    else numP = nmo*nmo - oP;

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
  }

  AO_RI = block_matrix(maxPshell,norbs*norbs);
  MO_RI = block_matrix(maxPshell,nmo*nmo);

  next_DF_MO = PSIO_ZERO; 

  for (int Pshell=0; Pshell < ribasis_->nshell(); ++Pshell) {
    int numPshell = ribasis_->shell(Pshell)->nfunction();
 
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
 
  }
  
  free_block(AO_RI); 
  free_block(MO_RI);

  zero_disk(PSIF_SAPT_AB_DF_INTS,"AB RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.noccA*calc_info_.noccB);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.nvirB*calc_info_.noccA);
  zero_disk(PSIF_SAPT_AB_DF_INTS,"RB RI Integrals",(char *) &(zeros[0]),
            calc_info_.nrio,calc_info_.noccB*calc_info_.nvirA);

  psio_address next_DF_AB = PSIO_ZERO;
  psio_address next_DF_AS = PSIO_ZERO;
  psio_address next_DF_RB = PSIO_ZERO;

  for (int i=0,oP=0; oP < nmo*nmo; ++i, oP += numP) {
    if ((i+1)*temp_size < nmo*nmo) numP = temp_size;
    else numP = nmo*nmo - oP;
      
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
  }

  free(zeros);
  free(halftrans);

  free_block(J_mhalf);

  psio_->close(PSIF_SAPT_TEMP,0);

  free(Schwartz);
  free(DFSchwartz);

  delete []eri;

  fflush(outfile);
}

void SAPT2B::oetrans()
{
  double **Sij, **SiB;

  Sij = block_matrix(calc_info_.nso,calc_info_.nso);
  SiB = block_matrix(calc_info_.nso,calc_info_.nmo);

  for (int i=0,ij=0; i<calc_info_.nso; i++) {
    for (int j=0; j<=i; j++,ij++) {
      Sij[i][j] = calc_info_.S[ij];
      Sij[j][i] = calc_info_.S[ij];
      }}

  free(calc_info_.S);

  C_DGEMM('N','N',calc_info_.nso,calc_info_.nmo,calc_info_.nso,1.0,
          &(Sij[0][0]),calc_info_.nso,&(calc_info_.CB[0][0]),calc_info_.nmo,
          0.0,&(SiB[0][0]),calc_info_.nmo);

  free_block(Sij);
  calc_info_.S_AB = block_matrix(calc_info_.nmo,calc_info_.nmo);

  C_DGEMM('T','N',calc_info_.nmo,calc_info_.nmo,calc_info_.nso,1.0,
          &(calc_info_.CA[0][0]),calc_info_.nmo,&(SiB[0][0]),calc_info_.nmo,
          0.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo);

  free_block(SiB);

  double **VAij, **VBij;
  double **VAiB, **VBAj;

  VAij = block_matrix(calc_info_.nso,calc_info_.nso);
  VAiB = block_matrix(calc_info_.nso,calc_info_.nmo);

  for (int i=0,ij=0; i<calc_info_.nso; i++) {
    for (int j=0; j<=i; j++,ij++) {
      VAij[i][j] = calc_info_.VA[ij];
      VAij[j][i] = calc_info_.VA[ij];
      }}

  free(calc_info_.VA);

  C_DGEMM('N','N',calc_info_.nso,calc_info_.nmo,calc_info_.nso,1.0,
          &(VAij[0][0]),calc_info_.nso,&(calc_info_.CB[0][0]),calc_info_.nmo,
          0.0,&(VAiB[0][0]),calc_info_.nmo);

  free_block(VAij);

  calc_info_.VABB = block_matrix(calc_info_.nmo,calc_info_.nmo);

  C_DGEMM('T','N',calc_info_.nmo,calc_info_.nmo,calc_info_.nso,1.0,
          &(calc_info_.CB[0][0]),calc_info_.nmo,&(VAiB[0][0]),calc_info_.nmo,
          0.0,&(calc_info_.VABB[0][0]),calc_info_.nmo);

  calc_info_.VAAB = block_matrix(calc_info_.nmo,calc_info_.nmo);

  C_DGEMM('T','N',calc_info_.nmo,calc_info_.nmo,calc_info_.nso,1.0,
          &(calc_info_.CA[0][0]),calc_info_.nmo,&(VAiB[0][0]),calc_info_.nmo,
          0.0,&(calc_info_.VAAB[0][0]),calc_info_.nmo);

  free_block(VAiB);

  VBij = block_matrix(calc_info_.nso,calc_info_.nso);
  VBAj = block_matrix(calc_info_.nmo,calc_info_.nso);

  for (int i=0,ij=0; i<calc_info_.nso; i++) {
    for (int j=0; j<=i; j++,ij++) {
      VBij[i][j] = calc_info_.VB[ij];
      VBij[j][i] = calc_info_.VB[ij];
      }}

  free(calc_info_.VB);

  C_DGEMM('T','N',calc_info_.nmo,calc_info_.nso,calc_info_.nso,1.0,
          &(calc_info_.CA[0][0]),calc_info_.nmo,&(VBij[0][0]),calc_info_.nso,
          0.0,&(VBAj[0][0]),calc_info_.nso);

  free_block(VBij);

  calc_info_.VBAA = block_matrix(calc_info_.nmo,calc_info_.nmo);

  C_DGEMM('N','N',calc_info_.nmo,calc_info_.nmo,calc_info_.nso,1.0,
          &(VBAj[0][0]),calc_info_.nso,&(calc_info_.CA[0][0]),calc_info_.nmo,
          0.0,&(calc_info_.VBAA[0][0]),calc_info_.nmo);

  calc_info_.VBAB = block_matrix(calc_info_.nmo,calc_info_.nmo);

  C_DGEMM('N','N',calc_info_.nmo,calc_info_.nmo,calc_info_.nso,1.0,
          &(VBAj[0][0]),calc_info_.nso,&(calc_info_.CB[0][0]),calc_info_.nmo,
          0.0,&(calc_info_.VBAB[0][0]),calc_info_.nmo);

  free_block(VBAj);
}

void SAPT2B::w_ints()
{
  double **B_p_A = get_diag_AA_ints(1);
  double **B_p_B = get_diag_BB_ints(1);

  calc_info_.diagAA = init_array(calc_info_.nrio);
  calc_info_.diagBB = init_array(calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA; a++){
    C_DAXPY(calc_info_.nrio,1.0,&(B_p_A[a][0]),1,calc_info_.diagAA,1);
  }

  for(int b=0; b<calc_info_.noccB; b++){
    C_DAXPY(calc_info_.nrio,1.0,&(B_p_B[b][0]),1,calc_info_.diagBB,1);
  }

  free_block(B_p_A);
  free_block(B_p_B);

  if (workflow_.W_ov) {

    calc_info_.WBAR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

    for(int a=0; a<calc_info_.noccA; a++){
      C_DAXPY(calc_info_.nvirA,1.0,&(calc_info_.VBAA[a][calc_info_.noccA]),1,
        &(calc_info_.WBAR[a][0]),1);
    }

    double **B_p_AR = get_AR_ints(0);

    C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nri,2.0,
      &(B_p_AR[0][0]),calc_info_.nrio,calc_info_.diagBB,1,1.0,
      &(calc_info_.WBAR[0][0]),1);

    free_block(B_p_AR);

    calc_info_.WABS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

    for(int b=0; b<calc_info_.noccB; b++){
      C_DAXPY(calc_info_.nvirB,1.0,&(calc_info_.VABB[b][calc_info_.noccB]),1,
        &(calc_info_.WABS[b][0]),1);
    }

    double **B_p_BS = get_BS_ints(0);

    C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nri,2.0,
      &(B_p_BS[0][0]),calc_info_.nrio,calc_info_.diagAA,1,1.0,
      &(calc_info_.WABS[0][0]),1);

    free_block(B_p_BS);

  }

  if (workflow_.W_oo) {

    calc_info_.WBAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

    for(int a=0; a<calc_info_.noccA; a++){
      C_DAXPY(calc_info_.noccA,1.0,&(calc_info_.VBAA[a][0]),1,
        &(calc_info_.WBAA[a][0]),1);
    }

    double **B_p_AA = get_AA_ints(0);

    C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nri,2.0,
      &(B_p_AA[0][0]),calc_info_.nrio,calc_info_.diagBB,1,1.0,
      &(calc_info_.WBAA[0][0]),1);

    free_block(B_p_AA);

    calc_info_.WABB = block_matrix(calc_info_.noccB,calc_info_.noccB);

    for(int b=0; b<calc_info_.noccB; b++){
      C_DAXPY(calc_info_.noccB,1.0,&(calc_info_.VABB[b][0]),1,
        &(calc_info_.WABB[b][0]),1);
    }

    double **B_p_BB = get_BB_ints(0);

    C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nri,2.0,
      &(B_p_BB[0][0]),calc_info_.nrio,calc_info_.diagAA,1,1.0,
      &(calc_info_.WABB[0][0]),1);

    free_block(B_p_BB);

  }

  if (workflow_.W_vv) {

    calc_info_.WBRR = block_matrix(calc_info_.nvirA,calc_info_.nvirA);

    for(int r=0; r<calc_info_.nvirA; r++){
      C_DAXPY(calc_info_.nvirA,1.0,
        &(calc_info_.VBAA[r+calc_info_.noccA][calc_info_.noccA]),1,
        &(calc_info_.WBRR[r][0]),1);
    }

    double **B_p_RR = get_RR_ints(0);

    C_DGEMV('n',calc_info_.nvirA*calc_info_.nvirA,calc_info_.nri,2.0,
      &(B_p_RR[0][0]),calc_info_.nrio,calc_info_.diagBB,1,1.0,
      &(calc_info_.WBRR[0][0]),1);

    free_block(B_p_RR);

    calc_info_.WASS = block_matrix(calc_info_.nvirB,calc_info_.nvirB);

    for(int s=0; s<calc_info_.nvirB; s++){
      C_DAXPY(calc_info_.nvirB,1.0,
        &(calc_info_.VABB[s+calc_info_.noccB][calc_info_.noccB]),1,
        &(calc_info_.WASS[s][0]),1);
    }

    double **B_p_SS = get_SS_ints(0);

    C_DGEMV('n',calc_info_.nvirB*calc_info_.nvirB,calc_info_.nri,2.0,
      &(B_p_SS[0][0]),calc_info_.nrio,calc_info_.diagAA,1,1.0,
      &(calc_info_.WASS[0][0]),1);

    free_block(B_p_SS);

  }
  
}

}}
