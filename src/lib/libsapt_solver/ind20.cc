#include "sapt0.h"
#include "sapt2.h"

#include <boost/shared_ptr.hpp>
#include <libpsio/aiohandler.h>

using namespace boost;

namespace psi { namespace sapt {

void SAPT0::ind20()
{
  double **tAR = block_matrix(noccA_,nvirA_);

  for (int a=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++) {
      tAR[a][r] = wBAR_[a][r]/(evalsA_[a] - evalsA_[r+noccA_]);
  }}

  double indA_B = 2.0*C_DDOT(noccA_*nvirA_,tAR[0],1,wBAR_[0],1);

  if (no_response_) {
    CHFA_ = block_matrix(noccA_,nvirA_);
    C_DCOPY(noccA_*nvirA_,tAR[0],1,CHFA_[0],1);
  }

  free_block(tAR);

  double **tBS = block_matrix(noccB_,nvirB_);

  for (int b=0; b<noccB_; b++) {
    for (int s=0; s<nvirB_; s++) {
      tBS[b][s] = wABS_[b][s]/(evalsB_[b] - evalsB_[s+noccB_]);
  }}

  double indB_A = 2.0*C_DDOT(noccB_*nvirB_,tBS[0],1,wABS_[0],1);

  if (no_response_) {
    CHFB_ = block_matrix(noccB_,nvirB_);
    C_DCOPY(noccB_*nvirB_,tBS[0],1,CHFB_[0],1);
  }

  free_block(tBS);

  e_ind20_ = indA_B + indB_A;

  if (print_) {
    fprintf(outfile,"    Ind20 (A<-B)        = %18.12lf H\n",indA_B);
    fprintf(outfile,"    Ind20 (B<-A)        = %18.12lf H\n",indB_A);
    fprintf(outfile,"    Ind20               = %18.12lf H\n",e_ind20_);
    fflush(outfile);
  }
}

void SAPT0::ind20r()
{
  if (aio_cphf_) {
    ind20rA_B_aio();
    ind20rB_A_aio();
  }
  else {
    ind20rA_B();
    ind20rB_A();
  }

  double indA_B, indB_A;

  indA_B = 2.0*C_DDOT(noccA_*nvirA_,CHFA_[0],1,wBAR_[0],1);
  indB_A = 2.0*C_DDOT(noccB_*nvirB_,CHFB_[0],1,wABS_[0],1);

  e_ind20_ = indA_B + indB_A;

  if (print_) {
    fprintf(outfile,"    Ind20,r (A<-B)      = %18.12lf H\n",indA_B);
    fprintf(outfile,"    Ind20,r (B<-A)      = %18.12lf H\n",indB_A);
    fprintf(outfile,"    Ind20,r             = %18.12lf H\n",e_ind20_);
    fflush(outfile);
  }
}

void SAPT0::ind20rA_B()
{
  time_t start = time(NULL);
  time_t stop;
  int iter=0;
  double E_old, E;
  double conv,dE;
  double *tAR_old = init_array(noccA_*nvirA_);
  double *tAR_new = init_array(noccA_*nvirA_);

  double alpha, beta;
  double *R_old = init_array(noccA_*nvirA_);
  double *R_new = init_array(noccA_*nvirA_);
  double *z_old = init_array(noccA_*nvirA_);
  double *z_new = init_array(noccA_*nvirA_);
  double *Ax = init_array(noccA_*nvirA_);

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  for (int a=0, ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++, ar++) {
      tAR_old[ar] = wBAR_[a][r]/(evalsA_[a] - evalsA_[r+noccA_]);
  }}

  E_old = 0.0;
  conv = 1.0;

  if (debug_) {
    fprintf(outfile,"\n    Maxiter = %d\n",maxiter_);
    fprintf(outfile,"    D converge = %lE\n",d_conv_);
    fprintf(outfile,"    E converge = %lE\n",e_conv_);
  }

  if (print_)
    fprintf(outfile,"\n    Iter     Energy (mH)           dE (mH)          Residual      Time (s)\n");

  SAPTDFInts C_p_AA = set_C_AA();
  SAPTDFInts C_p_AR = set_C_AR();
  SAPTDFInts C_p_RR = set_C_RR();

  double *X = init_array(ndf_);
  double **xAA = block_matrix(nthreads,noccA_*noccA_);
  double **xAR = block_matrix(nthreads,noccA_*nvirA_);
  double **xRR = block_matrix(nthreads,nvirA_*nvirA_);
  double **tAR_dump = block_matrix(nthreads,noccA_*nvirA_);

  do {

    memset(&(Ax[0]),'\0',sizeof(double)*noccA_*nvirA_);
    memset(&(tAR_dump[0][0]),'\0',sizeof(double)*nthreads*noccA_*nvirA_);

    Iterator AR_iter = get_iterator(mem_,&C_p_AR);

    for (int i=0, off=0; i<AR_iter.num_blocks; i++) {
      read_block(&AR_iter,&C_p_AR);

      C_DGEMV('n',AR_iter.curr_size,noccA_*nvirA_,1.0,&(C_p_AR.B_p_[0][0]),
        noccA_*nvirA_,tAR_old,1,0.0,&(X[0]),1);
      C_DGEMV('t',AR_iter.curr_size,noccA_*nvirA_,-4.0,&(C_p_AR.B_p_[0][0]),
        noccA_*nvirA_,&(X[0]),1,1.0,Ax,1);

#pragma omp parallel
{
#pragma omp for private(rank)
      for (int j=0; j<AR_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

        C_DGEMM('N','T',noccA_,noccA_,nvirA_,1.0,&(C_p_AR.B_p_[j][0]),nvirA_,
          tAR_old,nvirA_,0.0,xAA[rank],noccA_);
        C_DGEMM('N','N',noccA_,nvirA_,noccA_,1.0,xAA[rank],noccA_,
          &(C_p_AR.B_p_[j][0]),nvirA_,1.0,tAR_dump[rank],nvirA_);
      }
}
      off += AR_iter.curr_size;
    }

    C_p_AR.clear();

    Iterator RR_iter = get_iterator(mem_,&C_p_AA,&C_p_RR);

    for (int i=0, off=0; i<RR_iter.num_blocks; i++) {
      read_block(&RR_iter,&C_p_AA,&C_p_RR);

#pragma omp parallel
{
#pragma omp for private(rank)
      for (int j=0; j<RR_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

        for (int a=0, ab=0; a<nvirA_; a++) {
          for (int b=0; b<=a; b++) {
            xRR[rank][a*nvirA_+b] = C_p_RR.B_p_[j][ab];
            xRR[rank][b*nvirA_+a] = C_p_RR.B_p_[j][ab++];
        }}

        C_DGEMM('N','N',noccA_,nvirA_,nvirA_,1.0,tAR_old,nvirA_,
          &(xRR[rank][0]),nvirA_,0.0,xAR[rank],nvirA_);
        C_DGEMM('N','N',noccA_,nvirA_,noccA_,1.0,&(C_p_AA.B_p_[j][0]),noccA_,
          xAR[rank],nvirA_,1.0,tAR_dump[rank],nvirA_);
      }
}
      off += RR_iter.curr_size;
    }

    C_p_AA.clear();
    C_p_RR.clear();

    for (int n=0; n<nthreads; n++)
      C_DAXPY(noccA_*nvirA_,1.0,tAR_dump[n],1,Ax,1);

    for (int a=0, ar=0; a<noccA_; a++) {
      for (int r=0; r<nvirA_; r++, ar++) {
        Ax[ar] += (evalsA_[a] - evalsA_[r+noccA_])*tAR_old[ar];
    }}

    if (!iter) {
      C_DCOPY(noccA_*nvirA_,wBAR_[0],1,R_old,1);
      C_DAXPY(noccA_*nvirA_,-1.0,Ax,1,R_old,1);
      for (int a=0, ar=0; a<noccA_; a++) {
        for (int r=0; r<nvirA_; r++, ar++) {
          z_old[ar] = R_old[ar]/(evalsA_[a] - evalsA_[r+noccA_]);
      }}
      C_DCOPY(noccA_*nvirA_,tAR_old,1,tAR_new,1);
      C_DCOPY(noccA_*nvirA_,z_old,1,tAR_old,1);
    }
    else {
      alpha = C_DDOT(noccA_*nvirA_,R_old,1,z_old,1);
      alpha /= C_DDOT(noccA_*nvirA_,tAR_old,1,Ax,1);
      C_DAXPY(noccA_*nvirA_,alpha,tAR_old,1,tAR_new,1);
      C_DCOPY(noccA_*nvirA_,R_old,1,R_new,1);
      C_DAXPY(noccA_*nvirA_,-alpha,Ax,1,R_new,1);
      for (int a=0, ar=0; a<noccA_; a++) {
        for (int r=0; r<nvirA_; r++, ar++) {
          z_new[ar] = R_new[ar]/(evalsA_[a] - evalsA_[r+noccA_]);
      }}
      beta = C_DDOT(noccA_*nvirA_,R_new,1,z_new,1); 
      beta /= C_DDOT(noccA_*nvirA_,R_old,1,z_old,1); 
      C_DSCAL(noccA_*nvirA_,beta,tAR_old,1);
      C_DAXPY(noccA_*nvirA_,1.0,z_new,1,tAR_old,1);
      C_DCOPY(noccA_*nvirA_,z_new,1,z_old,1);
      C_DCOPY(noccA_*nvirA_,R_new,1,R_old,1);
    }

    E = 2.0*C_DDOT(noccA_*nvirA_,tAR_new,1,&(wBAR_[0][0]),1);

    conv = C_DDOT(noccA_*nvirA_,R_old,1,R_old,1);
    conv = sqrt(conv);
    dE = E_old-E;

    iter++;
    stop = time(NULL);
    if (print_) {
      fprintf(outfile,"    %4d %16.8lf %17.9lf %17.9lf    %10ld\n",
        iter,E*1000.0,dE*1000.0,conv*1000.0,stop-start);
      fflush(outfile);
    }
    E_old = E;
  }
  while((conv > d_conv_ || fabs(dE) > e_conv_) && iter < maxiter_);

  if ((conv <= d_conv_) && (fabs(dE) <= e_conv_)) {
    if (print_)
      fprintf(outfile,"\n    CHF Iterations converged\n\n");
    }
  else {
    fprintf(outfile,"\n    CHF Iterations did not converge\n\n");
    }

  CHFA_ = block_matrix(noccA_,nvirA_);
  C_DCOPY(noccA_*nvirA_,tAR_new,1,CHFA_[0],1);

  free(tAR_new);
  free(tAR_old);
  free_block(tAR_dump);
  free(X);
  free(R_old);
  free(R_new);
  free(z_old);
  free(z_new);
  free(Ax);
  free_block(xAA);
  free_block(xAR);
  free_block(xRR);
}

void SAPT0::ind20rB_A()
{
  time_t start = time(NULL);
  time_t stop;
  int iter=0;
  double E_old, E;
  double conv,dE;
  double *tBS_old = init_array(noccB_*nvirB_);
  double *tBS_new = init_array(noccB_*nvirB_);

  double alpha, beta;
  double *R_old = init_array(noccB_*nvirB_);
  double *R_new = init_array(noccB_*nvirB_);
  double *z_old = init_array(noccB_*nvirB_);
  double *z_new = init_array(noccB_*nvirB_);
  double *Ax = init_array(noccB_*nvirB_);

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  for (int b=0, bs=0; b<noccB_; b++) {
    for (int s=0; s<nvirB_; s++, bs++) {
      tBS_old[bs] = wABS_[b][s]/(evalsB_[b] - evalsB_[s+noccB_]);
  }}

  E_old = 0.0;
  conv = 1.0;

  if (debug_) {
    fprintf(outfile,"    Maxiter = %d\n",maxiter_);
    fprintf(outfile,"    D converge = %lE\n",d_conv_);
    fprintf(outfile,"    E converge = %lE\n",e_conv_);
  }

  if (print_)
    fprintf(outfile,"\n    Iter     Energy (mH)           dE (mH)          Residual      Time (s)\n");

  SAPTDFInts C_p_BB = set_C_BB();
  SAPTDFInts C_p_BS = set_C_BS();
  SAPTDFInts C_p_SS = set_C_SS();

  double *X = init_array(ndf_);
  double **xBB = block_matrix(nthreads,noccB_*noccB_);
  double **xBS = block_matrix(nthreads,noccB_*nvirB_);
  double **xSS = block_matrix(nthreads,nvirB_*nvirB_);
  double **tBS_dump = block_matrix(nthreads,noccB_*nvirB_);

  do {

    memset(&(Ax[0]),'\0',sizeof(double)*noccB_*nvirB_);
    memset(&(tBS_dump[0][0]),'\0',sizeof(double)*nthreads*noccB_*nvirB_);

    Iterator BS_iter = get_iterator(mem_,&C_p_BS);

    for (int i=0, off=0; i<BS_iter.num_blocks; i++) {
      read_block(&BS_iter,&C_p_BS);

      C_DGEMV('n',BS_iter.curr_size,noccB_*nvirB_,1.0,&(C_p_BS.B_p_[0][0]),
        noccB_*nvirB_,tBS_old,1,0.0,&(X[0]),1);
      C_DGEMV('t',BS_iter.curr_size,noccB_*nvirB_,-4.0,&(C_p_BS.B_p_[0][0]),
        noccB_*nvirB_,&(X[0]),1,1.0,Ax,1);

#pragma omp parallel
{
#pragma omp for private(rank)
      for (int j=0; j<BS_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

        C_DGEMM('N','T',noccB_,noccB_,nvirB_,1.0,&(C_p_BS.B_p_[j][0]),nvirB_,
          tBS_old,nvirB_,0.0,xBB[rank],noccB_);
        C_DGEMM('N','N',noccB_,nvirB_,noccB_,1.0,xBB[rank],noccB_,
          &(C_p_BS.B_p_[j][0]),nvirB_,1.0,tBS_dump[rank],nvirB_);
      }
}
      off += BS_iter.curr_size;
    }

    C_p_BS.clear();

    Iterator SS_iter = get_iterator(mem_,&C_p_BB,&C_p_SS);

    for (int i=0, off=0; i<SS_iter.num_blocks; i++) {
      read_block(&SS_iter,&C_p_BB,&C_p_SS);

#pragma omp parallel
{
#pragma omp for private(rank)
      for (int j=0; j<SS_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

        for (int a=0, ab=0; a<nvirB_; a++) {
          for (int b=0; b<=a; b++) {
            xSS[rank][a*nvirB_+b] = C_p_SS.B_p_[j][ab];
            xSS[rank][b*nvirB_+a] = C_p_SS.B_p_[j][ab++];
        }}

        C_DGEMM('N','N',noccB_,nvirB_,nvirB_,1.0,tBS_old,nvirB_,
          &(xSS[rank][0]),nvirB_,0.0,xBS[rank],nvirB_);
        C_DGEMM('N','N',noccB_,nvirB_,noccB_,1.0,&(C_p_BB.B_p_[j][0]),noccB_,
          xBS[rank],nvirB_,1.0,tBS_dump[rank],nvirB_);
      }
}
      off += SS_iter.curr_size;
    }

    C_p_BB.clear();
    C_p_SS.clear();

    for (int n=0; n<nthreads; n++)
      C_DAXPY(noccB_*nvirB_,1.0,tBS_dump[n],1,Ax,1);

    for (int b=0, bs=0; b<noccB_; b++) {
      for (int s=0; s<nvirB_; s++, bs++) {
        Ax[bs] += (evalsB_[b] - evalsB_[s+noccB_])*tBS_old[bs];
    }}

    if (!iter) {
      C_DCOPY(noccB_*nvirB_,wABS_[0],1,R_old,1);
      C_DAXPY(noccB_*nvirB_,-1.0,Ax,1,R_old,1);
      for (int b=0, bs=0; b<noccB_; b++) {
        for (int s=0; s<nvirB_; s++, bs++) {
          z_old[bs] = R_old[bs]/(evalsB_[b] - evalsB_[s+noccB_]);
      }}
      C_DCOPY(noccB_*nvirB_,tBS_old,1,tBS_new,1);
      C_DCOPY(noccB_*nvirB_,z_old,1,tBS_old,1);
    }
    else {
      alpha = C_DDOT(noccB_*nvirB_,R_old,1,z_old,1);
      alpha /= C_DDOT(noccB_*nvirB_,tBS_old,1,Ax,1);
      C_DAXPY(noccB_*nvirB_,alpha,tBS_old,1,tBS_new,1);
      C_DCOPY(noccB_*nvirB_,R_old,1,R_new,1);
      C_DAXPY(noccB_*nvirB_,-alpha,Ax,1,R_new,1);
      for (int b=0, bs=0; b<noccB_; b++) {
        for (int s=0; s<nvirB_; s++, bs++) {
          z_new[bs] = R_new[bs]/(evalsB_[b] - evalsB_[s+noccB_]);
      }}
      beta = C_DDOT(noccB_*nvirB_,R_new,1,z_new,1); 
      beta /= C_DDOT(noccB_*nvirB_,R_old,1,z_old,1); 
      C_DSCAL(noccB_*nvirB_,beta,tBS_old,1);
      C_DAXPY(noccB_*nvirB_,1.0,z_new,1,tBS_old,1);
      C_DCOPY(noccB_*nvirB_,z_new,1,z_old,1);
      C_DCOPY(noccB_*nvirB_,R_new,1,R_old,1);
    }

    E = 2.0*C_DDOT(noccB_*nvirB_,tBS_new,1,&(wABS_[0][0]),1);

    conv = C_DDOT(noccB_*nvirB_,R_old,1,R_old,1);
    conv = sqrt(conv);
    dE = E_old-E;

    iter++;
    stop = time(NULL);
    if (print_) {
      fprintf(outfile,"    %4d %16.8lf %17.9lf %17.9lf    %10ld\n",
        iter,E*1000.0,dE*1000.0,conv*1000.0,stop-start);
      fflush(outfile);
    }
    E_old = E;
  }
  while((conv > d_conv_ || fabs(dE) > e_conv_) && iter < maxiter_);

  if ((conv <= d_conv_) && (fabs(dE) <= e_conv_)) {
    if (print_)
      fprintf(outfile,"\n    CHF Iterations converged\n\n");
    }
  else {
    fprintf(outfile,"\n    CHF Iterations did not converge\n\n");
    }

  CHFB_ = block_matrix(noccB_,nvirB_);
  C_DCOPY(noccB_*nvirB_,tBS_new,1,CHFB_[0],1);

  free(tBS_new);
  free(tBS_old);
  free_block(tBS_dump);
  free(X);
  free(R_old);
  free(R_new);
  free(z_old);
  free(z_new);
  free(Ax);
  free_block(xBB);
  free_block(xBS);
  free_block(xSS);
}

void SAPT0::ind20rA_B_aio()
{
  time_t start = time(NULL);
  time_t stop;
  int iter=0;
  double E_old, E;
  double conv,dE;
  double *tAR_old = init_array(noccA_*nvirA_);
  double *tAR_new = init_array(noccA_*nvirA_);

  double alpha, beta;
  double *R_old = init_array(noccA_*nvirA_);
  double *R_new = init_array(noccA_*nvirA_);
  double *z_old = init_array(noccA_*nvirA_);
  double *z_new = init_array(noccA_*nvirA_);
  double *Ax = init_array(noccA_*nvirA_);

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  for (int a=0, ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++, ar++) {
      tAR_old[ar] = wBAR_[a][r]/(evalsA_[a] - evalsA_[r+noccA_]);
  }}

  E_old = 0.0;
  conv = 1.0;

  if (debug_) {
    fprintf(outfile,"\n    Maxiter = %d\n",maxiter_);
    fprintf(outfile,"    D converge = %lE\n",d_conv_);
    fprintf(outfile,"    E converge = %lE\n",e_conv_);
  }

  if (print_)
    fprintf(outfile,"\n    Iter     Energy (mH)           dE (mH)          Residual      Time (s)\n");

  SAPTDFInts C_p_AR = set_C_AR();

  double *X = init_array(ndf_);
  double **xAA = block_matrix(nthreads,noccA_*noccA_);
  double **xAR = block_matrix(nthreads,noccA_*nvirA_);
  double **xRR = block_matrix(nthreads,nvirA_*nvirA_);
  double **tAR_dump = block_matrix(nthreads,noccA_*nvirA_);

  long int block_length = mem_ / (2L*noccA_*noccA_ + 2L*nvirA_*(nvirA_+1)/2);
  if (block_length > ndf_) block_length = ndf_;

  int num_blocks = ndf_ / block_length;
  if (ndf_ % block_length) num_blocks++;

  boost::shared_ptr<AIOHandler> aio(new AIOHandler(psio_));

  double **C_p_AA[2];
  double **C_p_RR[2];

  psio_address next_C_p_AA;
  psio_address next_C_p_RR;

  do {

    memset(&(Ax[0]),'\0',sizeof(double)*noccA_*nvirA_);
    memset(&(tAR_dump[0][0]),'\0',sizeof(double)*nthreads*noccA_*nvirA_);

    Iterator AR_iter = get_iterator(mem_,&C_p_AR);

    for (int i=0, off=0; i<AR_iter.num_blocks; i++) {
      read_block(&AR_iter,&C_p_AR);

      C_DGEMV('n',AR_iter.curr_size,noccA_*nvirA_,1.0,&(C_p_AR.B_p_[0][0]),
        noccA_*nvirA_,tAR_old,1,0.0,&(X[0]),1);
      C_DGEMV('t',AR_iter.curr_size,noccA_*nvirA_,-4.0,&(C_p_AR.B_p_[0][0]),
        noccA_*nvirA_,&(X[0]),1,1.0,Ax,1);

#pragma omp parallel
{
#pragma omp for private(rank)
      for (int j=0; j<AR_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

        C_DGEMM('N','T',noccA_,noccA_,nvirA_,1.0,&(C_p_AR.B_p_[j][0]),nvirA_,
          tAR_old,nvirA_,0.0,xAA[rank],noccA_);
        C_DGEMM('N','N',noccA_,nvirA_,noccA_,1.0,xAA[rank],noccA_,
          &(C_p_AR.B_p_[j][0]),nvirA_,1.0,tAR_dump[rank],nvirA_);
      }
}
      off += AR_iter.curr_size;
    }

    C_p_AR.clear();

    C_p_AA[0] = block_matrix(block_length,noccA_*noccA_);
    C_p_AA[1] = block_matrix(block_length,noccA_*noccA_);
    C_p_RR[0] = block_matrix(block_length,nvirA_*(nvirA_+1)/2);
    C_p_RR[1] = block_matrix(block_length,nvirA_*(nvirA_+1)/2);

    next_C_p_AA = PSIO_ZERO;
    next_C_p_RR = PSIO_ZERO;

    psio_->read(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *)
      &(C_p_AA[0][0][0]),sizeof(double)*block_length*noccA_*noccA_,
      next_C_p_AA,&next_C_p_AA);
    psio_->read(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",(char *)
      &(C_p_RR[0][0][0]),sizeof(double)*block_length*nvirA_*(nvirA_+1)/2,
      next_C_p_RR,&next_C_p_RR);

    for (int i=0; i<num_blocks; i++) {

      if (i < num_blocks-1) {
        int read_length = block_length;
        if (i == num_blocks-2 && ndf_ % block_length)
          read_length = ndf_ % block_length;
          aio->read(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",(char *)
            &(C_p_AA[(i+1)%2][0][0]),sizeof(double)*read_length*noccA_*noccA_,
            next_C_p_AA,&next_C_p_AA);
          aio->read(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",(char *)
            &(C_p_RR[(i+1)%2][0][0]),sizeof(double)*read_length*nvirA_*
            (nvirA_+1)/2,next_C_p_RR,&next_C_p_RR);
      }

      int loopsize = block_length;
      if (i == num_blocks-1 && ndf_ % block_length)
        loopsize = ndf_ % block_length;

#pragma omp parallel
{
#pragma omp for private(rank)
      for (int j=0; j<loopsize; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

        for (int a=0, ab=0; a<nvirA_; a++) {
          for (int b=0; b<=a; b++) {
            xRR[rank][a*nvirA_+b] = C_p_RR[i%2][j][ab];
            xRR[rank][b*nvirA_+a] = C_p_RR[i%2][j][ab++];
        }}

        C_DGEMM('N','N',noccA_,nvirA_,nvirA_,1.0,tAR_old,nvirA_,
          &(xRR[rank][0]),nvirA_,0.0,xAR[rank],nvirA_);
        C_DGEMM('N','N',noccA_,nvirA_,noccA_,1.0,&(C_p_AA[i%2][j][0]),noccA_,
          xAR[rank],nvirA_,1.0,tAR_dump[rank],nvirA_);
      }
}

      if (i < num_blocks-1)
        aio->synchronize();
    }

    free_block(C_p_AA[0]);
    free_block(C_p_AA[1]);
    free_block(C_p_RR[0]);
    free_block(C_p_RR[1]);

    for (int n=0; n<nthreads; n++)
      C_DAXPY(noccA_*nvirA_,1.0,tAR_dump[n],1,Ax,1);

    for (int a=0, ar=0; a<noccA_; a++) {
      for (int r=0; r<nvirA_; r++, ar++) {
        Ax[ar] += (evalsA_[a] - evalsA_[r+noccA_])*tAR_old[ar];
    }}

    if (!iter) {
      C_DCOPY(noccA_*nvirA_,wBAR_[0],1,R_old,1);
      C_DAXPY(noccA_*nvirA_,-1.0,Ax,1,R_old,1);
      for (int a=0, ar=0; a<noccA_; a++) {
        for (int r=0; r<nvirA_; r++, ar++) {
          z_old[ar] = R_old[ar]/(evalsA_[a] - evalsA_[r+noccA_]);
      }}
      C_DCOPY(noccA_*nvirA_,tAR_old,1,tAR_new,1);
      C_DCOPY(noccA_*nvirA_,z_old,1,tAR_old,1);
    }
    else {
      alpha = C_DDOT(noccA_*nvirA_,R_old,1,z_old,1);
      alpha /= C_DDOT(noccA_*nvirA_,tAR_old,1,Ax,1);
      C_DAXPY(noccA_*nvirA_,alpha,tAR_old,1,tAR_new,1);
      C_DCOPY(noccA_*nvirA_,R_old,1,R_new,1);
      C_DAXPY(noccA_*nvirA_,-alpha,Ax,1,R_new,1);
      for (int a=0, ar=0; a<noccA_; a++) {
        for (int r=0; r<nvirA_; r++, ar++) {
          z_new[ar] = R_new[ar]/(evalsA_[a] - evalsA_[r+noccA_]);
      }}
      beta = C_DDOT(noccA_*nvirA_,R_new,1,z_new,1); 
      beta /= C_DDOT(noccA_*nvirA_,R_old,1,z_old,1); 
      C_DSCAL(noccA_*nvirA_,beta,tAR_old,1);
      C_DAXPY(noccA_*nvirA_,1.0,z_new,1,tAR_old,1);
      C_DCOPY(noccA_*nvirA_,z_new,1,z_old,1);
      C_DCOPY(noccA_*nvirA_,R_new,1,R_old,1);
    }

    E = 2.0*C_DDOT(noccA_*nvirA_,tAR_new,1,&(wBAR_[0][0]),1);

    conv = C_DDOT(noccA_*nvirA_,R_old,1,R_old,1);
    conv = sqrt(conv);
    dE = E_old-E;

    iter++;
    stop = time(NULL);
    if (print_) {
      fprintf(outfile,"    %4d %16.8lf %17.9lf %17.9lf    %10ld\n",
        iter,E*1000.0,dE*1000.0,conv*1000.0,stop-start);
      fflush(outfile);
    }
    E_old = E;
  }
  while((conv > d_conv_ || fabs(dE) > e_conv_) && iter < maxiter_);

  if ((conv <= d_conv_) && (fabs(dE) <= e_conv_)) {
    if (print_)
      fprintf(outfile,"\n    CHF Iterations converged\n\n");
    }
  else {
    fprintf(outfile,"\n    CHF Iterations did not converge\n\n");
    }

  CHFA_ = block_matrix(noccA_,nvirA_);
  C_DCOPY(noccA_*nvirA_,tAR_new,1,CHFA_[0],1);

  free(tAR_new);
  free(tAR_old);
  free_block(tAR_dump);
  free(X);
  free(R_old);
  free(R_new);
  free(z_old);
  free(z_new);
  free(Ax);
  free_block(xAA);
  free_block(xAR);
  free_block(xRR);
}

void SAPT0::ind20rB_A_aio()
{
  time_t start = time(NULL);
  time_t stop;
  int iter=0;
  double E_old, E;
  double conv,dE;
  double *tBS_old = init_array(noccB_*nvirB_);
  double *tBS_new = init_array(noccB_*nvirB_);

  double alpha, beta;
  double *R_old = init_array(noccB_*nvirB_);
  double *R_new = init_array(noccB_*nvirB_);
  double *z_old = init_array(noccB_*nvirB_);
  double *z_new = init_array(noccB_*nvirB_);
  double *Ax = init_array(noccB_*nvirB_);

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  for (int b=0, bs=0; b<noccB_; b++) {
    for (int s=0; s<nvirB_; s++, bs++) {
      tBS_old[bs] = wABS_[b][s]/(evalsB_[b] - evalsB_[s+noccB_]);
  }}

  E_old = 0.0;
  conv = 1.0;

  if (debug_) {
    fprintf(outfile,"    Maxiter = %d\n",maxiter_);
    fprintf(outfile,"    D converge = %lE\n",d_conv_);
    fprintf(outfile,"    E converge = %lE\n",e_conv_);
  }

  if (print_)
    fprintf(outfile,"\n    Iter     Energy (mH)           dE (mH)          Residual      Time (s)\n");

  SAPTDFInts C_p_BS = set_C_BS();

  double *X = init_array(ndf_);
  double **xBB = block_matrix(nthreads,noccB_*noccB_);
  double **xBS = block_matrix(nthreads,noccB_*nvirB_);
  double **xSS = block_matrix(nthreads,nvirB_*nvirB_);
  double **tBS_dump = block_matrix(nthreads,noccB_*nvirB_);

  long int block_length = mem_ / (2L*noccB_*noccB_ + 2L*nvirB_*(nvirB_+1)/2);
  if (block_length > ndf_) block_length = ndf_;

  int num_blocks = ndf_ / block_length;
  if (ndf_ % block_length) num_blocks++;

  boost::shared_ptr<AIOHandler> aio(new AIOHandler(psio_));

  double **C_p_BB[2];
  double **C_p_SS[2];

  psio_address next_C_p_BB;
  psio_address next_C_p_SS;

  do {

    memset(&(Ax[0]),'\0',sizeof(double)*noccB_*nvirB_);
    memset(&(tBS_dump[0][0]),'\0',sizeof(double)*nthreads*noccB_*nvirB_);

    Iterator BS_iter = get_iterator(mem_,&C_p_BS);

    for (int i=0, off=0; i<BS_iter.num_blocks; i++) {
      read_block(&BS_iter,&C_p_BS);

      C_DGEMV('n',BS_iter.curr_size,noccB_*nvirB_,1.0,&(C_p_BS.B_p_[0][0]),
        noccB_*nvirB_,tBS_old,1,0.0,&(X[0]),1);
      C_DGEMV('t',BS_iter.curr_size,noccB_*nvirB_,-4.0,&(C_p_BS.B_p_[0][0]),
        noccB_*nvirB_,&(X[0]),1,1.0,Ax,1);

#pragma omp parallel
{
#pragma omp for private(rank)
      for (int j=0; j<BS_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

        C_DGEMM('N','T',noccB_,noccB_,nvirB_,1.0,&(C_p_BS.B_p_[j][0]),nvirB_,
          tBS_old,nvirB_,0.0,xBB[rank],noccB_);
        C_DGEMM('N','N',noccB_,nvirB_,noccB_,1.0,xBB[rank],noccB_,
          &(C_p_BS.B_p_[j][0]),nvirB_,1.0,tBS_dump[rank],nvirB_);
      }
}
      off += BS_iter.curr_size;
    }

    C_p_BS.clear();

    C_p_BB[0] = block_matrix(block_length,noccB_*noccB_);
    C_p_BB[1] = block_matrix(block_length,noccB_*noccB_);
    C_p_SS[0] = block_matrix(block_length,nvirB_*(nvirB_+1)/2);
    C_p_SS[1] = block_matrix(block_length,nvirB_*(nvirB_+1)/2);

    next_C_p_BB = PSIO_ZERO;
    next_C_p_SS = PSIO_ZERO;

    psio_->read(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *)
      &(C_p_BB[0][0][0]),sizeof(double)*block_length*noccB_*noccB_,
      next_C_p_BB,&next_C_p_BB);
    psio_->read(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",(char *)
      &(C_p_SS[0][0][0]),sizeof(double)*block_length*nvirB_*(nvirB_+1)/2,
      next_C_p_SS,&next_C_p_SS);

    for (int i=0; i<num_blocks; i++) {

      if (i < num_blocks-1) {
        int read_length = block_length;
        if (i == num_blocks-2 && ndf_ % block_length)
          read_length = ndf_ % block_length;
          aio->read(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",(char *)
            &(C_p_BB[(i+1)%2][0][0]),sizeof(double)*read_length*noccB_*noccB_,
            next_C_p_BB,&next_C_p_BB);
          aio->read(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",(char *)
            &(C_p_SS[(i+1)%2][0][0]),sizeof(double)*read_length*nvirB_*
            (nvirB_+1)/2,next_C_p_SS,&next_C_p_SS);
      }

      int loopsize = block_length;
      if (i == num_blocks-1 && ndf_ % block_length)
        loopsize = ndf_ % block_length;

#pragma omp parallel
{
#pragma omp for private(rank)
      for (int j=0; j<loopsize; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

        for (int a=0, ab=0; a<nvirB_; a++) {
          for (int b=0; b<=a; b++) {
            xSS[rank][a*nvirB_+b] = C_p_SS[i%2][j][ab];
            xSS[rank][b*nvirB_+a] = C_p_SS[i%2][j][ab++];
        }}

        C_DGEMM('N','N',noccB_,nvirB_,nvirB_,1.0,tBS_old,nvirB_,
          &(xSS[rank][0]),nvirB_,0.0,xBS[rank],nvirB_);
        C_DGEMM('N','N',noccB_,nvirB_,noccB_,1.0,&(C_p_BB[i%2][j][0]),noccB_,
          xBS[rank],nvirB_,1.0,tBS_dump[rank],nvirB_);
      }
}

      if (i < num_blocks-1)
        aio->synchronize();
    }

    free_block(C_p_BB[0]);
    free_block(C_p_BB[1]);
    free_block(C_p_SS[0]);
    free_block(C_p_SS[1]);

    for (int n=0; n<nthreads; n++)
      C_DAXPY(noccB_*nvirB_,1.0,tBS_dump[n],1,Ax,1);

    for (int b=0, bs=0; b<noccB_; b++) {
      for (int s=0; s<nvirB_; s++, bs++) {
        Ax[bs] += (evalsB_[b] - evalsB_[s+noccB_])*tBS_old[bs];
    }}

    if (!iter) {
      C_DCOPY(noccB_*nvirB_,wABS_[0],1,R_old,1);
      C_DAXPY(noccB_*nvirB_,-1.0,Ax,1,R_old,1);
      for (int b=0, bs=0; b<noccB_; b++) {
        for (int s=0; s<nvirB_; s++, bs++) {
          z_old[bs] = R_old[bs]/(evalsB_[b] - evalsB_[s+noccB_]);
      }}
      C_DCOPY(noccB_*nvirB_,tBS_old,1,tBS_new,1);
      C_DCOPY(noccB_*nvirB_,z_old,1,tBS_old,1);
    }
    else {
      alpha = C_DDOT(noccB_*nvirB_,R_old,1,z_old,1);
      alpha /= C_DDOT(noccB_*nvirB_,tBS_old,1,Ax,1);
      C_DAXPY(noccB_*nvirB_,alpha,tBS_old,1,tBS_new,1);
      C_DCOPY(noccB_*nvirB_,R_old,1,R_new,1);
      C_DAXPY(noccB_*nvirB_,-alpha,Ax,1,R_new,1);
      for (int b=0, bs=0; b<noccB_; b++) {
        for (int s=0; s<nvirB_; s++, bs++) {
          z_new[bs] = R_new[bs]/(evalsB_[b] - evalsB_[s+noccB_]);
      }}
      beta = C_DDOT(noccB_*nvirB_,R_new,1,z_new,1); 
      beta /= C_DDOT(noccB_*nvirB_,R_old,1,z_old,1); 
      C_DSCAL(noccB_*nvirB_,beta,tBS_old,1);
      C_DAXPY(noccB_*nvirB_,1.0,z_new,1,tBS_old,1);
      C_DCOPY(noccB_*nvirB_,z_new,1,z_old,1);
      C_DCOPY(noccB_*nvirB_,R_new,1,R_old,1);
    }

    E = 2.0*C_DDOT(noccB_*nvirB_,tBS_new,1,&(wABS_[0][0]),1);

    conv = C_DDOT(noccB_*nvirB_,R_old,1,R_old,1);
    conv = sqrt(conv);
    dE = E_old-E;

    iter++;
    stop = time(NULL);
    if (print_) {
      fprintf(outfile,"    %4d %16.8lf %17.9lf %17.9lf    %10ld\n",
        iter,E*1000.0,dE*1000.0,conv*1000.0,stop-start);
      fflush(outfile);
    }
    E_old = E;
  }
  while((conv > d_conv_ || fabs(dE) > e_conv_) && iter < maxiter_);

  if ((conv <= d_conv_) && (fabs(dE) <= e_conv_)) {
    if (print_)
      fprintf(outfile,"\n    CHF Iterations converged\n\n");
    }
  else {
    fprintf(outfile,"\n    CHF Iterations did not converge\n\n");
    }

  CHFB_ = block_matrix(noccB_,nvirB_);
  C_DCOPY(noccB_*nvirB_,tBS_new,1,CHFB_[0],1);

  free(tBS_new);
  free(tBS_old);
  free_block(tBS_dump);
  free(X);
  free(R_old);
  free(R_new);
  free(z_old);
  free(z_new);
  free(Ax);
  free_block(xBB);
  free_block(xBS);
  free_block(xSS);
}

void SAPT2::ind20r()
{
  CHFA_ = block_matrix(noccA_,nvirA_);

  cphf_solver(CHFA_,wBAR_,evalsA_,PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
    "AR RI Integrals","RR RI Integrals",noccA_,nvirA_);

  CHFB_ = block_matrix(noccB_,nvirB_);

  cphf_solver(CHFB_,wABS_,evalsB_,PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
    "BS RI Integrals","SS RI Integrals",noccB_,nvirB_);

  double indA_B, indB_A;

  indA_B = 2.0*C_DDOT(noccA_*nvirA_,CHFA_[0],1,wBAR_[0],1);
  indB_A = 2.0*C_DDOT(noccB_*nvirB_,CHFB_[0],1,wABS_[0],1);

  e_ind20_ = indA_B + indB_A;

  if (print_) {
    fprintf(outfile,"    Ind20,r (A<-B)      = %18.12lf H\n",indA_B);
    fprintf(outfile,"    Ind20,r (B<-A)      = %18.12lf H\n",indB_A);
    fprintf(outfile,"    Ind20,r             = %18.12lf H\n",e_ind20_);
    fflush(outfile);
  }
}

void SAPT2::cphf_solver(double **tAR, double **wBAR, double *evals,
  int intfile, const char *AAints, const char *ARints, const char *RRints,
  int nocc, int nvir)
{
  double **B_p_AR = block_matrix(nocc*nvir,ndf_+3);

  psio_->read_entry(intfile,ARints,(char *) &(B_p_AR[0][0]),
    sizeof(double)*nocc*nvir*(ndf_+3));

  double **Amat = block_matrix(nocc*nvir,nocc*nvir);

  C_DGEMM('N','T',nocc*nvir,nocc*nvir,ndf_,-4.0,B_p_AR[0],ndf_+3,B_p_AR[0],
    ndf_+3,0.0,Amat[0],nocc*nvir);

  for (int a=0, ar=0; a<nocc; a++) {
    for (int r=0; r<nvir; r++, ar++) {
      C_DGEMM('N','T',nocc,nvir,ndf_,1.0,B_p_AR[r],nvir*(ndf_+3),B_p_AR[a*nvir],
        ndf_+3,1.0,Amat[ar],nvir);
  }}

  free_block(B_p_AR);

  double **B_p_AA = block_matrix(nocc*nocc,ndf_+3);
  double **B_p_R = block_matrix(nvir,ndf_+3);

  psio_->read_entry(intfile,AAints,(char *) &(B_p_AA[0][0]),
    sizeof(double)*nocc*nocc*(ndf_+3));

  psio_address next_PSIF = PSIO_ZERO;

  for (int r=0; r<nvir; r++) {
    psio_->read(intfile,RRints,(char *) &(B_p_R[0][0]),
      sizeof(double)*nvir*(ndf_+3),next_PSIF,&next_PSIF);
    for (int a=0; a<nocc; a++) {
      int ar = a*nvir+r;
      C_DGEMM('N','T',nocc,nvir,ndf_,1.0,B_p_AA[a*nocc],ndf_+3,B_p_R[0],
        ndf_+3,1.0,Amat[ar],nvir);
  }}

  free_block(B_p_AA);
  free_block(B_p_R);

  for (int a=0, ar=0; a<nocc; a++) {
    for (int r=0; r<nvir; r++, ar++) {
      Amat[ar][ar] += (evals[a] - evals[r+nocc]);
  }}

  int* ipiv = init_int_array(nocc*nvir);

  C_DCOPY(nocc*nvir,wBAR[0],1,tAR[0],1);
  C_DGESV(nocc*nvir,1,Amat[0],nocc*nvir,ipiv,tAR[0],nocc*nvir);

  free(ipiv);
  free_block(Amat);
}

}}

