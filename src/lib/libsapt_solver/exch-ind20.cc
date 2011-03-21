#include "sapt0.h"

namespace psi { namespace sapt {

void SAPT0::exch_ind20A_B()
{
  double ex1 = 0.0, ex2 = 0.0, ex3 = 0.0, ex4 = 0.0, ex5 = 0.0, ex6 = 0.0;
  double ex7 = 0.0, ex8 = 0.0, ex9 = 0.0, ex10 = 0.0, ex11 = 0.0, ex12 = 0.0;
  double ex13 = 0.0, exind = 0.0;

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  double **S_AB = block_matrix(noccA_,noccB_);
  double **S_RB = block_matrix(nvirA_,noccB_);

  for (int a=0; a<noccA_; a++) {
    C_DCOPY(noccB_,sAB_[a],1,S_AB[a],1);
  }

  for (int r=0; r<nvirA_; r++) {
    C_DCOPY(noccB_,sAB_[r+noccA_],1,S_RB[r],1);
  }

  SAPTDFInts A_p_AA = set_A_AA();
  SAPTDFInts B_p_BB = set_B_BB();
  SAPTDFInts A_p_AB = set_A_AB();
  SAPTDFInts B_p_AB = set_B_AB();
  SAPTDFInts A_p_AR = set_A_AR();
  SAPTDFInts B_p_RB = set_B_RB();

  double **xAR = block_matrix(nthreads,noccA_*nvirA_);

  Iterator E1_iter = get_iterator(mem_,&A_p_AB,&B_p_RB);

  for (int i=0; i<E1_iter.num_blocks; i++) {
    read_block(&E1_iter,&A_p_AB,&B_p_RB);

#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex1)
    for (int j=0; j<E1_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','T',noccA_,nvirA_,noccB_,1.0,A_p_AB.B_p_[j],noccB_,
        B_p_RB.B_p_[j],noccB_,0.0,xAR[rank],nvirA_);
      ex1 += C_DDOT(noccA_*nvirA_,xAR[rank],1,CHFA_[0],1);
    }
}
  }


  free_block(xAR);

  A_p_AB.clear();
  B_p_RB.clear();

  double **xRB = block_matrix(nvirA_,noccB_);
  double **yRB = block_matrix(nvirA_,noccB_);

  double **xAB = block_matrix(nthreads,noccA_*noccB_);
  double **yAB = block_matrix(nthreads,noccA_*noccB_);

  Iterator E2_iter = get_iterator(mem_,&A_p_AA,&B_p_RB);

  for (int i=0,off=0; i<E2_iter.num_blocks; i++) {
    read_block(&E2_iter,&A_p_AA,&B_p_RB);
    C_DGEMV('t',E2_iter.curr_size,nvirA_*noccB_,1.0,&(B_p_RB.B_p_[0][0]),
      nvirA_*noccB_,&(diagAA_[off]),1,1.0,xRB[0],1);

#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex5)
    for (int j=0; j<E2_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',noccA_,noccB_,nvirA_,1.0,CHFA_[0],nvirA_,
        &(B_p_RB.B_p_[j][0]),noccB_,0.0,xAB[rank],noccB_);
      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,&(A_p_AA.B_p_[j][0]),noccA_,
        S_AB[0],noccB_,0.0,yAB[rank],noccB_);
      ex5 -= C_DDOT(noccA_*noccB_,xAB[rank],1,yAB[rank],1);
    }
}
    off += E2_iter.curr_size;
  }

  C_DGEMM('T','N',nvirA_,noccB_,noccA_,1.0,CHFA_[0],nvirA_,S_AB[0],noccB_,
    0.0,&(yRB[0][0]),noccB_);

  ex2 = 2.0*C_DDOT(nvirA_*noccB_,xRB[0],1,yRB[0],1);

  free_block(xAB);
  free_block(yAB);
  free_block(xRB);

  A_p_AA.clear();
  B_p_RB.clear();

  double *X = init_array(ndf_+3);
  double *Y = init_array(ndf_+3);
  xRB = block_matrix(nthreads,nvirA_*noccB_);

  Iterator E3_iter = get_iterator(mem_,&A_p_AR,&B_p_AB);

  for (int i=0,off=0; i<E3_iter.num_blocks; i++) {
    read_block(&E3_iter,&A_p_AR,&B_p_AB);
    C_DGEMV('n',E3_iter.curr_size,noccA_*nvirA_,1.0,&(A_p_AR.B_p_[0][0]),
      noccA_*nvirA_,CHFA_[0],1,0.0,&(X[off]),1);
    C_DGEMV('n',E3_iter.curr_size,noccA_*noccB_,1.0,&(B_p_AB.B_p_[0][0]),
      noccA_*noccB_,S_AB[0],1,0.0,&(Y[off]),1);
#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex4)
    for (int j=0; j<E3_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('T','N',nvirA_,noccB_,noccA_,1.0,A_p_AR.B_p_[j],nvirA_,
        &(B_p_AB.B_p_[j][0]),noccB_,0.0,xRB[rank],noccB_);
      ex4 -= C_DDOT(nvirA_*noccB_,xRB[rank],1,yRB[0],1);
    }
}
    off += E3_iter.curr_size;
  }

  ex3 = 2.0*C_DDOT(ndf_+3,X,1,Y,1);

  free(X);
  free(Y);
  free_block(xRB);
  free_block(yRB);

  A_p_AR.clear();
  B_p_AB.clear();

  xAB = block_matrix(noccA_,noccB_);
  yAB = block_matrix(noccA_,noccB_);
  double **zAB = block_matrix(nthreads,noccA_*noccB_);

  C_DGEMM('N','N',noccA_,noccB_,nvirA_,1.0,CHFA_[0],nvirA_,S_RB[0],noccB_,
    0.0,&(xAB[0][0]),noccB_);

  Iterator E4_iter = get_iterator(mem_,&A_p_AB,&B_p_BB);

  for (int i=0,off=0; i<E4_iter.num_blocks; i++) {
    read_block(&E4_iter,&A_p_AB,&B_p_BB);
    C_DGEMV('t',E4_iter.curr_size,noccA_*noccB_,1.0,&(A_p_AB.B_p_[0][0]),
      noccA_*noccB_,&(diagBB_[off]),1,1.0,yAB[0],1);
#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex7)
    for (int j=0; j<E4_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,A_p_AB.B_p_[j],noccB_,
        &(B_p_BB.B_p_[j][0]),noccB_,0.0,zAB[rank],noccB_);
      ex7 += -1.0*C_DDOT(noccA_*noccB_,zAB[rank],1,xAB[0],1);
    }
}
    off += E4_iter.curr_size;
  }

  ex6 = 2.0*C_DDOT(noccA_*noccB_,yAB[0],1,xAB[0],1);

  A_p_AB.clear();
  B_p_BB.clear();

  free_block(yAB);

  double **xAA = block_matrix(noccA_,noccA_);
  double **xBB = block_matrix(noccB_,noccB_);
  double **yAA = block_matrix(noccA_,noccA_);
  double **yBB = block_matrix(noccB_,noccB_);
  yAB = block_matrix(nthreads,noccA_*noccB_);

  C_DGEMM('N','T',noccA_,noccA_,noccB_,1.0,xAB[0],noccB_,S_AB[0],noccB_,
    0.0,&(xAA[0][0]),noccA_);

  C_DGEMM('T','N',noccB_,noccB_,noccA_,1.0,xAB[0],noccB_,S_AB[0],noccB_,
    0.0,&(xBB[0][0]),noccB_);

  Iterator E5_iter = get_iterator(mem_,&A_p_AA,&B_p_BB);

  for (int i=0,off=0; i<E5_iter.num_blocks; i++) {
    read_block(&E5_iter,&A_p_AA,&B_p_BB);
    C_DGEMV('t',E5_iter.curr_size,noccA_*noccA_,1.0,&(A_p_AA.B_p_[0][0]),
      noccA_*noccA_,&(diagBB_[off]),1,1.0,yAA[0],1);
    C_DGEMV('t',E5_iter.curr_size,noccB_*noccB_,1.0,&(B_p_BB.B_p_[0][0]),
      noccB_*noccB_,&(diagAA_[off]),1,1.0,yBB[0],1);
#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex12)
    for (int j=0; j<E5_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,xAB[0],noccB_,
        B_p_BB.B_p_[j],noccB_,0.0,yAB[rank],noccB_);
      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,A_p_AA.B_p_[j],noccA_,
        S_AB[0],noccB_,0.0,zAB[rank],noccB_);
      ex12 += C_DDOT(noccA_*noccB_,yAB[rank],1,zAB[rank],1);
    }
}
    off += E5_iter.curr_size;
  }

  ex8 = -2.0*C_DDOT(noccB_*noccB_,xBB[0],1,yBB[0],1);
  ex10 = -2.0*C_DDOT(noccA_*noccA_,xAA[0],1,yAA[0],1);

  free_block(xAA);
  free_block(xBB);
  free_block(yAA);
  free_block(yBB);
  free_block(xAB);

  A_p_AA.clear();
  B_p_BB.clear();

  double **sAA = block_matrix(noccA_,noccA_);
  double **sBB = block_matrix(noccB_,noccB_);
  xAR = block_matrix(noccA_,nvirA_);
  double **yAR = block_matrix(noccA_,nvirA_);
  double **xBR = block_matrix(noccB_,nvirA_);
  X = init_array(ndf_+3);
  Y = init_array(ndf_+3);

  C_DGEMM('N','T',noccA_,noccA_,noccB_,1.0,S_AB[0],noccB_,S_AB[0],noccB_,
    0.0,&(sAA[0][0]),noccA_);

  C_DGEMM('T','N',noccB_,noccB_,noccA_,1.0,S_AB[0],noccB_,S_AB[0],noccB_,
    0.0,&(sBB[0][0]),noccB_);

  C_DGEMM('N','N',noccA_,nvirA_,noccA_,1.0,sAA[0],noccA_,CHFA_[0],nvirA_,
    0.0,&(xAR[0][0]),nvirA_);

  C_DGEMM('T','N',noccB_,nvirA_,noccA_,1.0,S_AB[0],noccB_,CHFA_[0],nvirA_,
    0.0,&(xBR[0][0]),nvirA_);

  Iterator E6_iter = get_iterator(mem_,&A_p_AR,&B_p_BB);

  for (int i=0,off=0; i<E6_iter.num_blocks; i++) {
    read_block(&E6_iter,&A_p_AR,&B_p_BB);
    C_DGEMV('n',E6_iter.curr_size,noccA_*nvirA_,1.0,&(A_p_AR.B_p_[0][0]),
      noccA_*nvirA_,CHFA_[0],1,0.0,&(X[off]),1);
    C_DGEMV('n',E6_iter.curr_size,noccB_*noccB_,1.0,&(B_p_BB.B_p_[0][0]),
      noccB_*noccB_,sBB[0],1,0.0,&(Y[off]),1);
    C_DGEMV('t',E6_iter.curr_size,noccA_*nvirA_,1.0,&(A_p_AR.B_p_[0][0]),
      noccA_*nvirA_,&(diagBB_[off]),1,1.0,yAR[0],1);
#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex13)
    for (int j=0; j<E6_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,S_AB[0],noccB_,
        B_p_BB.B_p_[j],noccB_,0.0,yAB[rank],noccB_);
      C_DGEMM('N','T',noccA_,noccB_,nvirA_,1.0,A_p_AR.B_p_[j],nvirA_,
        xBR[0],nvirA_,0.0,zAB[rank],noccB_);
      ex13 += C_DDOT(noccA_*noccB_,yAB[rank],1,zAB[rank],1);
    }
}
    off += E6_iter.curr_size;
  }

  ex9 = -2.0*C_DDOT(ndf_+3,X,1,Y,1);
  ex11 = -2.0*C_DDOT(noccA_*nvirA_,xAR[0],1,yAR[0],1);

  A_p_AR.clear();
  B_p_BB.clear();

  free_block(S_AB);
  free_block(S_RB);
  free_block(yAB);
  free_block(zAB);
  free_block(sAA);
  free_block(sBB);
  free_block(xAR);
  free_block(yAR);
  free_block(xBR);
  free(X);
  free(Y);

  exind = -2.0*(ex1+ex2+ex3+ex4+ex5+ex6+ex7+ex8+ex9+ex10+ex11+ex12+ex13);
  e_exch_ind20_ = exind;

  if (debug_) {
    fprintf(outfile,"\n    Ex1                 = %18.12lf H\n",ex1);
    fprintf(outfile,"    Ex2                 = %18.12lf H\n",ex2);
    fprintf(outfile,"    Ex3                 = %18.12lf H\n",ex3);
    fprintf(outfile,"    Ex4                 = %18.12lf H\n",ex4);
    fprintf(outfile,"    Ex5                 = %18.12lf H\n",ex5);
    fprintf(outfile,"    Ex6                 = %18.12lf H\n",ex6);
    fprintf(outfile,"    Ex7                 = %18.12lf H\n",ex7);
    fprintf(outfile,"    Ex8                 = %18.12lf H\n",ex8);
    fprintf(outfile,"    Ex9                 = %18.12lf H\n",ex9);
    fprintf(outfile,"    Ex10                = %18.12lf H\n",ex10);
    fprintf(outfile,"    Ex11                = %18.12lf H\n",ex11);
    fprintf(outfile,"    Ex12                = %18.12lf H\n",ex12);
    fprintf(outfile,"    Ex13                = %18.12lf H\n\n",ex13);
  }

  if (print_) {
    fprintf(outfile,"    Exch-Ind20,r (A<-B) = %18.12lf H\n",exind);
    fflush(outfile);
  }
}

void SAPT0::exch_ind20B_A()
{
  double ex1 = 0.0, ex2 = 0.0, ex3 = 0.0, ex4 = 0.0, ex5 = 0.0, ex6 = 0.0;
  double ex7 = 0.0, ex8 = 0.0, ex9 = 0.0, ex10 = 0.0, ex11 = 0.0, ex12 = 0.0;
  double ex13 = 0.0, exind = 0.0;

  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  double **S_AB = block_matrix(noccA_,noccB_);
  double **S_AS = block_matrix(noccA_,nvirB_);

  for (int a=0; a<noccA_; a++) {
    C_DCOPY(noccB_,sAB_[a],1,S_AB[a],1);
  }

  for (int a=0; a<noccA_; a++) {
    C_DCOPY(nvirB_,&(sAB_[a][noccB_]),1,S_AS[a],1);
  }

  SAPTDFInts A_p_AA = set_A_AA();
  SAPTDFInts B_p_BB = set_B_BB();
  SAPTDFInts A_p_AB = set_A_AB();
  SAPTDFInts B_p_AB = set_B_AB();
  SAPTDFInts B_p_BS = set_B_BS();
  SAPTDFInts A_p_AS = set_A_AS();

  double **xBS = block_matrix(nthreads,noccB_*nvirB_);

  Iterator E1_iter = get_iterator(mem_,&B_p_AB,&A_p_AS);

  for (int i=0; i<E1_iter.num_blocks; i++) {
    read_block(&E1_iter,&B_p_AB,&A_p_AS);
#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex1)
    for (int j=0; j<E1_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('T','N',noccB_,nvirB_,noccA_,1.0,B_p_AB.B_p_[j],noccB_,
        A_p_AS.B_p_[j],nvirB_,0.0,xBS[rank],nvirB_);
      ex1 += C_DDOT(noccB_*nvirB_,xBS[rank],1,CHFB_[0],1);
    }
}
  }

  free_block(xBS);

  B_p_AB.clear();
  A_p_AS.clear();

  double **xAS = block_matrix(noccA_,nvirB_);
  double **yAS = block_matrix(noccA_,nvirB_);

  double **xAB = block_matrix(nthreads,noccA_*noccB_);
  double **yAB = block_matrix(nthreads,noccA_*noccB_);

  Iterator E2_iter = get_iterator(mem_,&B_p_BB,&A_p_AS);

  for (int i=0,off=0; i<E2_iter.num_blocks; i++) {
    read_block(&E2_iter,&B_p_BB,&A_p_AS);
    C_DGEMV('t',E2_iter.curr_size,noccA_*nvirB_,1.0,&(A_p_AS.B_p_[0][0]),
      noccA_*nvirB_,&(diagBB_[off]),1,1.0,xAS[0],1);
#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex5)
    for (int j=0; j<E2_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','T',noccA_,noccB_,nvirB_,1.0,&(A_p_AS.B_p_[j][0]),nvirB_,
        CHFB_[0],nvirB_,0.0,xAB[rank],noccB_);
      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,S_AB[0],noccB_,
        &(B_p_BB.B_p_[j][0]),noccB_,0.0,yAB[rank],noccB_);
      ex5 -= C_DDOT(noccA_*noccB_,xAB[rank],1,yAB[rank],1);
    }
}
    off += E2_iter.curr_size;
  }

  C_DGEMM('N','N',noccA_,nvirB_,noccB_,1.0,S_AB[0],noccB_,CHFB_[0],nvirB_,
    0.0,&(yAS[0][0]),nvirB_);

  ex2 = 2.0*C_DDOT(noccA_*nvirB_,xAS[0],1,yAS[0],1);

  free_block(xAB);
  free_block(yAB);
  free_block(xAS);

  B_p_BB.clear();
  A_p_AS.clear();

  double *X = init_array(ndf_+3);
  double *Y = init_array(ndf_+3);
  xAS = block_matrix(nthreads,noccA_*nvirB_);

  Iterator E3_iter = get_iterator(mem_,&B_p_BS,&A_p_AB);

  for (int i=0,off=0; i<E3_iter.num_blocks; i++) {
    read_block(&E3_iter,&B_p_BS,&A_p_AB);
    C_DGEMV('n',E3_iter.curr_size,noccB_*nvirB_,1.0,&(B_p_BS.B_p_[0][0]),
      noccB_*nvirB_,CHFB_[0],1,0.0,&(X[off]),1);
    C_DGEMV('n',E3_iter.curr_size,noccA_*noccB_,1.0,&(A_p_AB.B_p_[0][0]),
      noccA_*noccB_,S_AB[0],1,0.0,&(Y[off]),1);
#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex4)
    for (int j=0; j<E3_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',noccA_,nvirB_,noccB_,1.0,&(A_p_AB.B_p_[j][0]),noccB_,
        B_p_BS.B_p_[j],nvirB_,0.0,xAS[rank],nvirB_);
      ex4 -= C_DDOT(noccA_*nvirB_,xAS[rank],1,yAS[0],1);
    }
}
    off += E3_iter.curr_size;
  }

  ex3 = 2.0*C_DDOT(ndf_+3,X,1,Y,1);

  free(X);
  free(Y);
  free_block(xAS);
  free_block(yAS);

  B_p_BS.clear();
  A_p_AB.clear();

  xAB = block_matrix(noccA_,noccB_);
  yAB = block_matrix(noccA_,noccB_);
  double **zAB = block_matrix(nthreads,noccA_*noccB_);

  C_DGEMM('N','T',noccA_,noccB_,nvirB_,1.0,S_AS[0],nvirB_,CHFB_[0],nvirB_,
    0.0,&(xAB[0][0]),noccB_);

  Iterator E4_iter = get_iterator(mem_,&A_p_AA,&B_p_AB);

  for (int i=0,off=0; i<E4_iter.num_blocks; i++) {
    read_block(&E4_iter,&A_p_AA,&B_p_AB);
    C_DGEMV('t',E4_iter.curr_size,noccA_*noccB_,1.0,&(B_p_AB.B_p_[0][0]),
      noccA_*noccB_,&(diagAA_[off]),1,1.0,yAB[0],1);
#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex7)
    for (int j=0; j<E4_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,A_p_AA.B_p_[j],noccA_,
        &(B_p_AB.B_p_[j][0]),noccB_,0.0,zAB[rank],noccB_);
      ex7 += -1.0*C_DDOT(noccA_*noccB_,zAB[rank],1,xAB[0],1);
    }
}
    off += E4_iter.curr_size;
  }

  ex6 = 2.0*C_DDOT(noccA_*noccB_,yAB[0],1,xAB[0],1);

  free_block(yAB);
  yAB = block_matrix(nthreads,noccA_*noccB_);

  A_p_AA.clear();
  B_p_AB.clear();

  double **xAA = block_matrix(noccA_,noccA_);
  double **xBB = block_matrix(noccB_,noccB_);
  double **yAA = block_matrix(noccA_,noccA_);
  double **yBB = block_matrix(noccB_,noccB_);

  C_DGEMM('N','T',noccA_,noccA_,noccB_,1.0,xAB[0],noccB_,S_AB[0],noccB_,
    0.0,&(xAA[0][0]),noccA_);

  C_DGEMM('T','N',noccB_,noccB_,noccA_,1.0,xAB[0],noccB_,S_AB[0],noccB_,
    0.0,&(xBB[0][0]),noccB_);

  Iterator E5_iter = get_iterator(mem_,&A_p_AA,&B_p_BB);

  for (int i=0,off=0; i<E5_iter.num_blocks; i++) {
    read_block(&E5_iter,&A_p_AA,&B_p_BB);
    C_DGEMV('t',E5_iter.curr_size,noccA_*noccA_,1.0,&(A_p_AA.B_p_[0][0]),
      noccA_*noccA_,&(diagBB_[off]),1,1.0,yAA[0],1);
    C_DGEMV('t',E5_iter.curr_size,noccB_*noccB_,1.0,&(B_p_BB.B_p_[0][0]),
      noccB_*noccB_,&(diagAA_[off]),1,1.0,yBB[0],1);
#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex12)
    for (int j=0; j<E5_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,S_AB[0],noccB_,
        B_p_BB.B_p_[j],noccB_,0.0,yAB[rank],noccB_);
      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,A_p_AA.B_p_[j],noccA_,
        xAB[0],noccB_,0.0,zAB[rank],noccB_);
      ex12 += C_DDOT(noccA_*noccB_,yAB[rank],1,zAB[rank],1);
    }
}
    off += E5_iter.curr_size;
  }

  ex8 = -2.0*C_DDOT(noccB_*noccB_,xBB[0],1,yBB[0],1);
  ex10 = -2.0*C_DDOT(noccA_*noccA_,xAA[0],1,yAA[0],1);

  free_block(xAA);
  free_block(xBB);
  free_block(yAA);
  free_block(yBB);
  free_block(xAB);

  A_p_AA.clear();
  B_p_BB.clear();

  double **sAA = block_matrix(noccA_,noccA_);
  double **sBB = block_matrix(noccB_,noccB_);
  xBS = block_matrix(noccB_,nvirB_);
  double **yBS = block_matrix(noccB_,nvirB_);
  xAS = block_matrix(noccA_,nvirB_);
  X = init_array(ndf_+3);
  Y = init_array(ndf_+3);

  C_DGEMM('N','T',noccA_,noccA_,noccB_,1.0,S_AB[0],noccB_,S_AB[0],noccB_,
    0.0,&(sAA[0][0]),noccA_);

  C_DGEMM('T','N',noccB_,noccB_,noccA_,1.0,S_AB[0],noccB_,S_AB[0],noccB_,
    0.0,&(sBB[0][0]),noccB_);

  C_DGEMM('N','N',noccB_,nvirB_,noccB_,1.0,sBB[0],noccB_,CHFB_[0],nvirB_,
    0.0,&(xBS[0][0]),nvirB_);

  C_DGEMM('N','N',noccA_,nvirB_,noccB_,1.0,S_AB[0],noccB_,CHFB_[0],nvirB_,
    0.0,&(xAS[0][0]),nvirB_);

  Iterator E6_iter = get_iterator(mem_,&A_p_AA,&B_p_BS);

  for (int i=0,off=0; i<E6_iter.num_blocks; i++) {
    read_block(&E6_iter,&A_p_AA,&B_p_BS);
    C_DGEMV('n',E6_iter.curr_size,noccB_*nvirB_,1.0,&(B_p_BS.B_p_[0][0]),
      noccB_*nvirB_,CHFB_[0],1,0.0,&(X[off]),1);
    C_DGEMV('n',E6_iter.curr_size,noccA_*noccA_,1.0,&(A_p_AA.B_p_[0][0]),
      noccA_*noccA_,sAA[0],1,0.0,&(Y[off]),1);
    C_DGEMV('t',E6_iter.curr_size,noccB_*nvirB_,1.0,&(B_p_BS.B_p_[0][0]),
      noccB_*nvirB_,&(diagAA_[off]),1,1.0,yBS[0],1);
#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:ex13)
    for (int j=0; j<E6_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,A_p_AA.B_p_[j],noccA_,
        S_AB[0],noccB_,0.0,yAB[rank],noccB_);
      C_DGEMM('N','T',noccA_,noccB_,nvirB_,1.0,xAS[0],nvirB_,
        B_p_BS.B_p_[j],nvirB_,0.0,zAB[rank],noccB_);
      ex13 += C_DDOT(noccA_*noccB_,yAB[rank],1,zAB[rank],1);
    }
}
    off += E6_iter.curr_size;
  }

  ex9 = -2.0*C_DDOT(ndf_+3,X,1,Y,1);
  ex11 = -2.0*C_DDOT(noccB_*nvirB_,xBS[0],1,yBS[0],1);

  A_p_AA.clear();
  B_p_BS.clear();

  free_block(S_AB);
  free_block(S_AS);
  free_block(yAB);
  free_block(zAB);
  free_block(sAA);
  free_block(sBB);
  free_block(xBS);
  free_block(yBS);
  free_block(xAS);
  free(X);
  free(Y);

  exind = -2.0*(ex1+ex2+ex3+ex4+ex5+ex6+ex7+ex8+ex9+ex10+ex11+ex12+ex13);
  e_exch_ind20_ += exind;

  if (debug_) {
    fprintf(outfile,"\n    Ex1                 = %18.12lf H\n",ex1);
    fprintf(outfile,"    Ex2                 = %18.12lf H\n",ex2);
    fprintf(outfile,"    Ex3                 = %18.12lf H\n",ex3);
    fprintf(outfile,"    Ex4                 = %18.12lf H\n",ex4);
    fprintf(outfile,"    Ex5                 = %18.12lf H\n",ex5);
    fprintf(outfile,"    Ex6                 = %18.12lf H\n",ex6);
    fprintf(outfile,"    Ex7                 = %18.12lf H\n",ex7);
    fprintf(outfile,"    Ex8                 = %18.12lf H\n",ex8);
    fprintf(outfile,"    Ex9                 = %18.12lf H\n",ex9);
    fprintf(outfile,"    Ex10                = %18.12lf H\n",ex10);
    fprintf(outfile,"    Ex11                = %18.12lf H\n",ex11);
    fprintf(outfile,"    Ex12                = %18.12lf H\n",ex12);
    fprintf(outfile,"    Ex13                = %18.12lf H\n\n",ex13);
  }

  if (print_) {
    fprintf(outfile,"    Exch-Ind20,r (B<-A) = %18.12lf H\n",exind);
    fprintf(outfile,"    Exch-Ind20,r        = %18.12lf H\n",e_exch_ind20_);
    fflush(outfile);
  }
} 

}}

