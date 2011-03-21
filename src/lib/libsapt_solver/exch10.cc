#include "sapt0.h"
    
namespace psi { namespace sapt {
  
void SAPT0::exch10_s2()
{ 
  double ex1 = 0.0, ex2 = 0.0, ex3 = 0.0, ex4 = 0.0, ex5 = 0.0, ex6 = 0.0;

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  SAPTDFInts A_p_AA = set_A_AA();
  SAPTDFInts B_p_BB = set_B_BB();
  SAPTDFInts A_p_AB = set_A_AB();
  SAPTDFInts B_p_AB = set_B_AB();

  Iterator E1_iter = get_iterator(mem_,&A_p_AB,&B_p_AB);

  for (int i=0; i<E1_iter.num_blocks; i++) {
    read_block(&E1_iter,&A_p_AB,&B_p_AB);
    ex1 += C_DDOT((long int) E1_iter.curr_size*noccA_*noccB_,A_p_AB.B_p_[0],1,
      B_p_AB.B_p_[0],1); 
  }

  A_p_AB.clear();
  B_p_AB.clear();

  double *Ap_diag = init_array(ndf_+3);
  double **X_AA = block_matrix(nthreads,noccA_*noccA_);

  Iterator E2_iter = get_iterator(mem_,&A_p_AA,&B_p_AB);

  for (int i=0, off=0; i<E2_iter.num_blocks; i++) {
    read_block(&E2_iter,&A_p_AA,&B_p_AB);

#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex2)
    for (int j=0; j<E2_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif 

      C_DGEMM('N','T',noccA_,noccA_,noccB_,1.0,B_p_AB.B_p_[j],noccB_,
        sAB_[0],nmo_,0.0,X_AA[rank],noccA_);
      ex2 -= C_DDOT(noccA_*noccA_,X_AA[rank],1,A_p_AA.B_p_[j],1);
      for (int a=0,aa=0; a<noccA_; a++, aa += noccA_+1) {
        Ap_diag[j+off] += X_AA[rank][aa];
    }}
}

    off += E2_iter.curr_size;
  }

  ex2 += 2.0*C_DDOT(ndf_+3,Ap_diag,1,diagAA_,1);

  free(Ap_diag);
  free_block(X_AA);

  A_p_AA.clear();
  B_p_AB.done();

  double *Bp_diag = init_array(ndf_+3);
  double **X_BB = block_matrix(nthreads,noccB_*noccB_);

  Iterator E3_iter = get_iterator(mem_,&A_p_AB,&B_p_BB);

  for (int i=0, off=0; i<E3_iter.num_blocks; i++) {
    read_block(&E3_iter,&A_p_AB,&B_p_BB);

#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex3)
    for (int j=0; j<E3_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('T','N',noccB_,noccB_,noccA_,1.0,A_p_AB.B_p_[j],noccB_,
        sAB_[0],nmo_,0.0,X_BB[rank],noccB_);
      ex3 -= C_DDOT(noccB_*noccB_,X_BB[rank],1,B_p_BB.B_p_[j],1);
      for (int b=0,bb=0; b<noccB_; b++, bb += noccB_+1) {
        Bp_diag[j+off] += X_BB[rank][bb];
    }}
}

    off += E3_iter.curr_size;
  }

  ex3 += 2.0*C_DDOT(ndf_+3,Bp_diag,1,diagBB_,1);

  free(Bp_diag);
  free_block(X_BB);

  A_p_AB.done();
  B_p_BB.clear();

  double **S_AA = block_matrix(noccA_,noccA_);

  C_DGEMM('N','T',noccA_,noccA_,noccB_,1.0,&(sAB_[0][0]),nmo_,&(sAB_[0][0]),
    nmo_,0.0,&(S_AA[0][0]),noccA_);

  double **S_BB = block_matrix(noccB_,noccB_);

  C_DGEMM('T','N',noccB_,noccB_,noccA_,1.0,&(sAB_[0][0]),nmo_,&(sAB_[0][0]),
    nmo_,0.0,&(S_BB[0][0]),noccB_);

  double **A_AB = block_matrix(nthreads,noccA_*noccB_);
  double **B_AB = block_matrix(nthreads,noccA_*noccB_);
  double *AA_ints = init_array(ndf_+3);
  double *BB_ints = init_array(ndf_+3);
  Iterator E4_iter = get_iterator(mem_,&A_p_AA,&B_p_BB);

  for (int i=0, off=0; i<E4_iter.num_blocks; i++) {
    read_block(&E4_iter,&A_p_AA,&B_p_BB);

    C_DGEMV('n',E4_iter.curr_size,noccA_*noccA_,1.0,&(A_p_AA.B_p_[0][0]),
      noccA_*noccA_,S_AA[0],1,0.0,&(AA_ints[off]),1);
    C_DGEMV('n',E4_iter.curr_size,noccB_*noccB_,1.0,&(B_p_BB.B_p_[0][0]),
      noccB_*noccB_,S_BB[0],1,0.0,&(BB_ints[off]),1);

#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex6)
    for (int j=0; j<E4_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,A_p_AA.B_p_[j],noccA_,
        sAB_[0],nmo_,0.0,A_AB[rank],noccB_);
      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,sAB_[0],nmo_,
        B_p_BB.B_p_[j],noccB_,0.0,B_AB[rank],noccB_);
      ex6 += C_DDOT(noccA_*noccB_,A_AB[rank],1,B_AB[rank],1);
    }
}
    off += E4_iter.curr_size;
  }

  ex4 = 2.0*C_DDOT(ndf_+3,BB_ints,1,diagAA_,1);
  ex5 = 2.0*C_DDOT(ndf_+3,AA_ints,1,diagBB_,1);

  A_p_AA.done();
  B_p_BB.done();

  free_block(S_AA);
  free_block(S_BB);

  free(AA_ints);
  free(BB_ints);
  free_block(A_AB);
  free_block(B_AB);

  e_exch10_s2_ = -2.0*(ex1+ex2+ex3-ex4-ex5+ex6);

  if (debug_) {
    fprintf(outfile,"\n    Ex1                 = %18.12lf H\n",ex1);
    fprintf(outfile,"    Ex2                 = %18.12lf H\n",ex2);
    fprintf(outfile,"    Ex3                 = %18.12lf H\n",ex3);
    fprintf(outfile,"    Ex4                 = %18.12lf H\n",ex4);
    fprintf(outfile,"    Ex5                 = %18.12lf H\n",ex5);
    fprintf(outfile,"    Ex6                 = %18.12lf H\n\n",ex6);
  }

  if (print_) {
    fprintf(outfile,"    Exch10 (S^2)        = %18.12lf H\n",e_exch10_s2_);
    fflush(outfile);
  }
}

void SAPT0::exch10()
{
  double ex1=0, ex2=0, ex3=0, ex4=0, ex5=0, ex6=0, ex7=0, ex8=0, ex9=0;

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  double **P = block_matrix(noccA_+noccB_,noccA_+noccB_);

  for (int i=0; i<noccA_+noccB_; i++)
    P[i][i] = 1.0;

  for (int a=0; a<noccA_; a++) {
    for (int b=0; b<noccB_; b++) {
      P[a][b+noccA_] = sAB_[a][b];
      P[b+noccA_][a] = sAB_[a][b];
  }}

  C_DPOTRF('L',noccA_+noccB_,P[0],noccA_+noccB_);
  C_DPOTRI('L',noccA_+noccB_,P[0],noccA_+noccB_);

  for (int i=0; i<noccA_+noccB_; i++)
    P[i][i] -= 1.0;

  double **pAA = block_matrix(noccA_,noccA_);
  double **pBB = block_matrix(noccB_,noccB_);
  double **pAB = block_matrix(noccA_,noccB_);

  for (int a1=0; a1<noccA_; a1++) {
    for (int a2=0; a2<noccA_; a2++) {
      if (a2 > a1) { pAA[a1][a2] = P[a1][a2]; }
      else { pAA[a1][a2] = P[a2][a1]; }
  }}

  for (int b1=0; b1<noccB_; b1++) {
    for (int b2=0; b2<noccB_; b2++) {
      if (b2 > b1) { pBB[b1][b2] = P[b1+noccA_][b2+noccA_]; }
      else { pBB[b1][b2] = P[b2+noccA_][b1+noccA_]; }
  }}

  for (int a=0; a<noccA_; a++) {
    for (int b=0; b<noccB_; b++) {
      pAB[a][b] = P[a][b+noccA_];
  }}

  free_block(P);

  double *W = init_array(ndf_+3);
  double *X = init_array(ndf_+3);
  double *Y = init_array(ndf_+3);
  double *Z = init_array(ndf_+3);

  double **xAB = block_matrix(nthreads,noccA_*noccB_);
  double **yAB = block_matrix(nthreads,noccA_*noccB_);

  SAPTDFInts A_p_AA = set_A_AA();
  SAPTDFInts B_p_BB = set_B_BB();
  SAPTDFInts A_p_AB = set_A_AB();
  SAPTDFInts B_p_AB = set_B_AB();

  Iterator E1_iter = get_iterator(mem_,&A_p_AB,&B_p_AB);

  for (int i=0, off=0; i<E1_iter.num_blocks; i++) {
    read_block(&E1_iter,&A_p_AB,&B_p_AB);
    ex1 += C_DDOT((long int) E1_iter.curr_size*noccA_*noccB_,A_p_AB.B_p_[0],1,
      B_p_AB.B_p_[0],1);

    C_DGEMV('n',E1_iter.curr_size,noccA_*noccB_,1.0,&(B_p_AB.B_p_[0][0]),
      noccA_*noccB_,pAB[0],1,0.0,&(W[off]),1);

    C_DGEMV('n',E1_iter.curr_size,noccA_*noccB_,1.0,&(A_p_AB.B_p_[0][0]),
      noccA_*noccB_,pAB[0],1,0.0,&(X[off]),1);

#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex2,ex3,ex8)
    for (int j=0; j<E1_iter.curr_size; j++) {

#ifdef _OPENMP 
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,pAA[0],noccA_,
        A_p_AB.B_p_[j],noccB_,0.0,xAB[rank],noccB_);
      ex2 += C_DDOT(noccA_*noccB_,xAB[rank],1,B_p_AB.B_p_[j],1);

      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,B_p_AB.B_p_[j],noccB_,
        pBB[0],noccB_,0.0,yAB[rank],noccB_);
      ex3 += C_DDOT(noccA_*noccB_,yAB[rank],1,A_p_AB.B_p_[j],1);

      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,xAB[rank],noccB_,
        pBB[0],noccB_,0.0,yAB[rank],noccB_);
      ex8 += C_DDOT(noccA_*noccB_,yAB[rank],1,B_p_AB.B_p_[j],1);
    }
}
    off += E1_iter.curr_size;
  }

  ex4 -= 2.0*C_DDOT(ndf_+3,W,1,diagAA_,1);
  ex5 -= 2.0*C_DDOT(ndf_+3,X,1,diagBB_,1);
  ex9 -= 2.0*C_DDOT(ndf_+3,W,1,X,1);

  A_p_AB.clear();
  B_p_AB.clear();

  Iterator E2_iter = get_iterator(mem_,&A_p_AA,&B_p_BB);
  
  for (int i=0, off=0; i<E2_iter.num_blocks; i++) {
    read_block(&E2_iter,&A_p_AA,&B_p_BB);

    C_DGEMV('n',E2_iter.curr_size,noccA_*noccA_,1.0,&(A_p_AA.B_p_[0][0]),
      noccA_*noccA_,pAA[0],1,0.0,&(Y[off]),1);

    C_DGEMV('n',E2_iter.curr_size,noccB_*noccB_,1.0,&(B_p_BB.B_p_[0][0]),
      noccB_*noccB_,pBB[0],1,0.0,&(Z[off]),1);

#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex9)
    for (int j=0; j<E2_iter.curr_size; j++) {

#ifdef _OPENMP 
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,A_p_AA.B_p_[j],noccA_,
        pAB[0],noccB_,0.0,xAB[rank],noccB_);
      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,pAB[0],noccB_,
        B_p_BB.B_p_[j],noccB_,0.0,yAB[rank],noccB_);
      ex9 += C_DDOT(noccA_*noccB_,xAB[rank],1,yAB[rank],1);
    }
}
    off += E2_iter.curr_size;
  }
  
  ex2 += -2.0*C_DDOT(ndf_+3,Y,1,diagBB_,1);
  ex3 += -2.0*C_DDOT(ndf_+3,Z,1,diagAA_,1);
  ex6 += -2.0*C_DDOT(ndf_+3,X,1,Z,1);
  ex7 += -2.0*C_DDOT(ndf_+3,W,1,Y,1);
  ex8 += -2.0*C_DDOT(ndf_+3,Y,1,Z,1);

  A_p_AA.clear();
  B_p_BB.clear();

  Iterator E3_iter = get_iterator(mem_,&A_p_AA,&B_p_AB);
  
  for (int i=0; i<E3_iter.num_blocks; i++) {
    read_block(&E3_iter,&A_p_AA,&B_p_AB);

#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex4,ex7)
    for (int j=0; j<E3_iter.curr_size; j++) {

#ifdef _OPENMP 
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,A_p_AA.B_p_[j],noccA_,
        pAB[0],noccB_,0.0,xAB[rank],noccB_);
      ex4 += C_DDOT(noccA_*noccB_,xAB[rank],1,B_p_AB.B_p_[j],1);

      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,pAA[0],noccA_,
        xAB[rank],noccB_,0.0,yAB[rank],noccB_);
      ex7 += C_DDOT(noccA_*noccB_,yAB[rank],1,B_p_AB.B_p_[j],1);
    }
}
  }

  A_p_AA.done();
  B_p_AB.done();

  Iterator E4_iter = get_iterator(mem_,&A_p_AB,&B_p_BB);
  
  for (int i=0; i<E4_iter.num_blocks; i++) {
    read_block(&E4_iter,&A_p_AB,&B_p_BB);
    
#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex5,ex6)
    for (int j=0; j<E4_iter.curr_size; j++) {

#ifdef _OPENMP 
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,pAB[0],noccB_,
        B_p_BB.B_p_[j],noccB_,0.0,xAB[rank],noccB_);
      ex5 += C_DDOT(noccA_*noccB_,xAB[rank],1,A_p_AB.B_p_[j],1);

      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,xAB[rank],noccB_,
        pBB[0],noccB_,0.0,yAB[rank],noccB_);
      ex6 += C_DDOT(noccA_*noccB_,yAB[rank],1,A_p_AB.B_p_[j],1);
    }
}
  } 

  A_p_AB.done();
  B_p_BB.done();

  free(W);
  free(X);
  free(Y);
  free(Z);
  free_block(xAB);
  free_block(yAB);

  free_block(pAA);
  free_block(pBB);
  free_block(pAB);

  e_exch10_ = -2.0*(ex1+ex2+ex3+ex4+ex5+ex6+ex7+ex8+ex9);

  if (debug_) {
    fprintf(outfile,"\n    Ex1                 = %18.12lf H\n",ex1);
    fprintf(outfile,"    Ex2                 = %18.12lf H\n",ex2);
    fprintf(outfile,"    Ex3                 = %18.12lf H\n",ex3);
    fprintf(outfile,"    Ex4                 = %18.12lf H\n",ex4);
    fprintf(outfile,"    Ex5                 = %18.12lf H\n",ex5);
    fprintf(outfile,"    Ex6                 = %18.12lf H\n",ex6);
    fprintf(outfile,"    Ex7                 = %18.12lf H\n",ex7);
    fprintf(outfile,"    Ex8                 = %18.12lf H\n",ex8);
    fprintf(outfile,"    Ex9                 = %18.12lf H\n\n",ex9);
  }

  if (print_) {
    fprintf(outfile,"    Exch10              = %18.12lf H\n",e_exch10_);
    fflush(outfile);
  }
}

}}
