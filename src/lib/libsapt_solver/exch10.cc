#include "sapt0.h"
    
namespace psi { namespace sapt {
  
void SAPT0::exch10_s2()
{ 
  double ex1 = 0.0, ex2 = 0.0, ex3 = 0.0, ex4 = 0.0, ex5 = 0.0, ex6 = 0.0;

  SAPTDFInts A_p_AA = set_A_AA();
  SAPTDFInts B_p_BB = set_B_BB();
  SAPTDFInts A_p_AB = set_A_AB();
  SAPTDFInts B_p_AB = set_B_AB();

  Iterator E1_iter = get_iterator(mem_,&A_p_AB,&B_p_AB);

  for (int i=0; i<E1_iter.num_blocks; i++) {
    read_block(&E1_iter,&A_p_AB,&B_p_AB);
    ex1 += C_DDOT(E1_iter.curr_size*noccA_*noccB_,A_p_AB.B_p_[0],1,
      B_p_AB.B_p_[0],1); 
  }

  A_p_AB.clear();
  B_p_AB.clear();

  double *Ap_diag = init_array(ndf_+3);
  double **X_AA = block_matrix(noccA_,noccA_);

  Iterator E2_iter = get_iterator(mem_,&A_p_AA,&B_p_AB);

  for (int i=0, off=0; i<E2_iter.num_blocks; i++) {
    read_block(&E2_iter,&A_p_AA,&B_p_AB);
    for (int j=0; j<E2_iter.curr_size; j++) {
      C_DGEMM('N','T',noccA_,noccA_,noccB_,1.0,B_p_AB.B_p_[j],noccB_,
        sAB_[0],nmo_,0.0,X_AA[0],noccA_);
      ex2 -= C_DDOT(noccA_*noccA_,X_AA[0],1,A_p_AA.B_p_[j],1);
      for (int a=0; a<noccA_; a++) {
        Ap_diag[j+off] += X_AA[a][a];
    }}
    off += E2_iter.curr_size;
  }

  ex2 += 2.0*C_DDOT(ndf_+3,Ap_diag,1,diagAA_,1);

  free(Ap_diag);
  free_block(X_AA);

  A_p_AA.clear();
  B_p_AB.done();

  double *Bp_diag = init_array(ndf_+3);
  double **X_BB = block_matrix(noccB_,noccB_);

  Iterator E3_iter = get_iterator(mem_,&A_p_AB,&B_p_BB);

  for (int i=0, off=0; i<E3_iter.num_blocks; i++) {
    read_block(&E3_iter,&A_p_AB,&B_p_BB);
    for (int j=0; j<E3_iter.curr_size; j++) {
      C_DGEMM('T','N',noccB_,noccB_,noccA_,1.0,A_p_AB.B_p_[j],noccB_,
        sAB_[0],nmo_,0.0,X_BB[0],noccB_);
      ex3 -= C_DDOT(noccB_*noccB_,X_BB[0],1,B_p_BB.B_p_[j],1);
      for (int b=0; b<noccB_; b++) {
        Bp_diag[j+off] += X_BB[b][b];
    }}
    off += E3_iter.curr_size;
  }

  ex3 += 2.0*C_DDOT(ndf_+3,Bp_diag,1,diagBB_,1);

  free(Bp_diag);
  free_block(X_BB);

  A_p_AB.done();
  B_p_BB.clear();

  double **S_BB = block_matrix(noccB_,noccB_);

  C_DGEMM('T','N',noccB_,noccB_,noccA_,1.0,&(sAB_[0][0]),nmo_,&(sAB_[0][0]),
    nmo_,0.0,&(S_BB[0][0]),noccB_);

  double *BB_ints = init_array(ndf_+3);

  Iterator E4_iter = get_iterator(mem_,&B_p_BB);

  for (int i=0,off=0; i<E4_iter.num_blocks; i++) {
    read_block(&E4_iter,&B_p_BB);

    C_DGEMV('n',E4_iter.curr_size,noccB_*noccB_,1.0,&(B_p_BB.B_p_[0][0]),
      noccB_*noccB_,S_BB[0],1,0.0,&(BB_ints[off]),1);

    off += E4_iter.curr_size;
  }

  ex4 = 2.0*C_DDOT(ndf_+3,BB_ints,1,diagAA_,1);

  B_p_BB.clear();
  free_block(S_BB);
  free(BB_ints);

  double **S_AA = block_matrix(noccA_,noccA_);

  C_DGEMM('N','T',noccA_,noccA_,noccB_,1.0,&(sAB_[0][0]),nmo_,&(sAB_[0][0]),
    nmo_,0.0,&(S_AA[0][0]),noccA_);

  double *AA_ints = init_array(ndf_+3);

  Iterator E5_iter = get_iterator(mem_,&A_p_AA);

  for (int i=0,off=0; i<E5_iter.num_blocks; i++) {
    read_block(&E5_iter,&A_p_AA);

    C_DGEMV('n',E5_iter.curr_size,noccA_*noccA_,1.0,&(A_p_AA.B_p_[0][0]),
      noccA_*noccA_,S_AA[0],1,0.0,&(AA_ints[off]),1);
    
    off += E5_iter.curr_size;
  }

  ex5 = 2.0*C_DDOT(ndf_+3,AA_ints,1,diagBB_,1);

  A_p_AA.clear();
  free_block(S_AA);
  free(AA_ints);

  double **A_AB = block_matrix(noccA_,noccB_);
  double **B_AB = block_matrix(noccA_,noccB_);

  Iterator E6_iter = get_iterator(mem_,&A_p_AA,&B_p_BB);

  for (int i=0, off=0; i<E6_iter.num_blocks; i++) {
    read_block(&E6_iter,&A_p_AA,&B_p_BB);
    for (int j=0; j<E6_iter.curr_size; j++) {
      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,A_p_AA.B_p_[j],noccA_,
        sAB_[0],nmo_,0.0,A_AB[0],noccB_);
      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,sAB_[0],nmo_,
        B_p_BB.B_p_[j],noccB_,0.0,B_AB[0],noccB_);
      ex6 += C_DDOT(noccA_*noccB_,A_AB[0],1,B_AB[0],1);
    }
    off += E6_iter.curr_size;
  }

  A_p_AA.done();
  B_p_BB.done();
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

  double *X = init_array(ndf_+3);
  double *Y = init_array(ndf_+3);

  double **xAB = block_matrix(noccA_,noccB_);
  double **yAB = block_matrix(noccA_,noccB_);

  SAPTDFInts A_p_AA = set_A_AA();
  SAPTDFInts B_p_BB = set_B_BB();
  SAPTDFInts A_p_AB = set_A_AB();
  SAPTDFInts B_p_AB = set_B_AB();

  Iterator E1_iter = get_iterator(mem_,&A_p_AB,&B_p_AB);

  for (int i=0; i<E1_iter.num_blocks; i++) {
    read_block(&E1_iter,&A_p_AB,&B_p_AB);
    ex1 += C_DDOT(E1_iter.curr_size*noccA_*noccB_,A_p_AB.B_p_[0],1,
      B_p_AB.B_p_[0],1);
  }

  A_p_AB.clear();
  B_p_AB.clear();

  Iterator E2_iter = get_iterator(mem_,&A_p_AA);
  
  for (int i=0, off=0; i<E2_iter.num_blocks; i++) {
    read_block(&E2_iter,&A_p_AA);

    C_DGEMV('n',E2_iter.curr_size,noccA_*noccA_,1.0,&(A_p_AA.B_p_[0][0]),
      noccA_*noccA_,pAA[0],1,0.0,&(X[off]),1);

    off += E2_iter.curr_size;
  }
  
  ex2 = -2.0*C_DDOT(ndf_+3,X,1,diagBB_,1);

  A_p_AA.clear();

  Iterator E3_iter = get_iterator(mem_,&A_p_AB,&B_p_AB);

  for (int i=0, off=0; i<E3_iter.num_blocks; i++) {
    read_block(&E3_iter,&A_p_AB,&B_p_AB);
    for (int j=0; j<E3_iter.curr_size; j++) {
      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,pAA[0],noccA_,
        A_p_AB.B_p_[j],noccB_,0.0,xAB[0],noccB_);
      ex2 += C_DDOT(noccA_*noccB_,xAB[0],1,B_p_AB.B_p_[j],1);
  }}

  A_p_AB.clear();
  B_p_AB.clear();

  Iterator E4_iter = get_iterator(mem_,&B_p_BB);

  for (int i=0, off=0; i<E4_iter.num_blocks; i++) {
    read_block(&E4_iter,&B_p_BB);

    C_DGEMV('n',E4_iter.curr_size,noccB_*noccB_,1.0,&(B_p_BB.B_p_[0][0]),
      noccB_*noccB_,pBB[0],1,0.0,&(X[off]),1);

    off += E4_iter.curr_size;
  }

  ex3 = -2.0*C_DDOT(ndf_+3,X,1,diagAA_,1);

  B_p_BB.clear();

  Iterator E5_iter = get_iterator(mem_,&A_p_AB,&B_p_AB);

  for (int i=0, off=0; i<E5_iter.num_blocks; i++) {
    read_block(&E5_iter,&A_p_AB,&B_p_AB);
    for (int j=0; j<E5_iter.curr_size; j++) {
      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,B_p_AB.B_p_[j],noccB_,
        pBB[0],noccB_,0.0,xAB[0],noccB_);
      ex3 += C_DDOT(noccA_*noccB_,xAB[0],1,A_p_AB.B_p_[j],1);
  }}

  A_p_AB.clear();
  B_p_AB.clear();

  Iterator E6_iter = get_iterator(mem_,&A_p_AA,&B_p_AB);

  for (int i=0, off=0; i<E6_iter.num_blocks; i++) {
    read_block(&E6_iter,&A_p_AA,&B_p_AB);

    C_DGEMV('n',E6_iter.curr_size,noccA_*noccB_,1.0,&(B_p_AB.B_p_[0][0]),
      noccA_*noccB_,pAB[0],1,0.0,&(X[off]),1);

    for (int j=0; j<E6_iter.curr_size; j++) {
      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,A_p_AA.B_p_[j],noccA_,
        pAB[0],noccB_,0.0,xAB[0],noccB_);
      ex4 += C_DDOT(noccA_*noccB_,xAB[0],1,B_p_AB.B_p_[j],1);
    }
    off += E6_iter.curr_size;
  }

  ex4 -= 2.0*C_DDOT(ndf_+3,X,1,diagAA_,1);

  A_p_AA.clear();
  B_p_AB.clear();

  Iterator E7_iter = get_iterator(mem_,&A_p_AB,&B_p_BB);

  for (int i=0, off=0; i<E7_iter.num_blocks; i++) {
    read_block(&E7_iter,&A_p_AB,&B_p_BB);
    
    C_DGEMV('n',E7_iter.curr_size,noccA_*noccB_,1.0,&(A_p_AB.B_p_[0][0]),
      noccA_*noccB_,pAB[0],1,0.0,&(X[off]),1);
    
    for (int j=0; j<E7_iter.curr_size; j++) {
      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,pAB[0],noccB_,
        B_p_BB.B_p_[j],noccB_,0.0,xAB[0],noccB_);
      ex5 += C_DDOT(noccA_*noccB_,xAB[0],1,A_p_AB.B_p_[j],1);
    }
    off += E7_iter.curr_size;
  }

  ex5 -= 2.0*C_DDOT(ndf_+3,X,1,diagBB_,1);

  A_p_AB.clear();
  B_p_BB.clear();

  Iterator E8_iter = get_iterator(mem_,&A_p_AB);

  for (int i=0, off=0; i<E8_iter.num_blocks; i++) {
    read_block(&E8_iter,&A_p_AB);

    C_DGEMV('n',E8_iter.curr_size,noccA_*noccB_,1.0,&(A_p_AB.B_p_[0][0]),
      noccA_*noccB_,pAB[0],1,0.0,&(X[off]),1);

    off += E8_iter.curr_size;
  }

  A_p_AB.clear();

  Iterator E9_iter = get_iterator(mem_,&B_p_BB);

  for (int i=0, off=0; i<E9_iter.num_blocks; i++) {
    read_block(&E9_iter,&B_p_BB);

    C_DGEMV('n',E9_iter.curr_size,noccB_*noccB_,1.0,&(B_p_BB.B_p_[0][0]),
      noccB_*noccB_,pBB[0],1,0.0,&(Y[off]),1);

    off += E9_iter.curr_size;
  }

  ex6 = -2.0*C_DDOT(ndf_+3,X,1,Y,1);

  B_p_BB.clear();

  Iterator E10_iter = get_iterator(mem_,&A_p_AB,&B_p_BB);

  for (int i=0, off=0; i<E10_iter.num_blocks; i++) {
    read_block(&E10_iter,&A_p_AB,&B_p_BB);
    for (int j=0; j<E10_iter.curr_size; j++) {
      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,pAB[0],noccB_,
        B_p_BB.B_p_[j],noccB_,0.0,xAB[0],noccB_);
      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,xAB[0],noccB_,
        pBB[0],noccB_,0.0,yAB[0],noccB_);
      ex6 += C_DDOT(noccA_*noccB_,yAB[0],1,A_p_AB.B_p_[j],1);
  }}

  A_p_AB.clear();
  B_p_BB.clear();

  Iterator E11_iter = get_iterator(mem_,&A_p_AA);

  for (int i=0, off=0; i<E11_iter.num_blocks; i++) {
    read_block(&E11_iter,&A_p_AA);

    C_DGEMV('n',E11_iter.curr_size,noccA_*noccA_,1.0,&(A_p_AA.B_p_[0][0]),
      noccA_*noccA_,pAA[0],1,0.0,&(X[off]),1);

    off += E11_iter.curr_size;
  }

  A_p_AA.clear();

  Iterator E12_iter = get_iterator(mem_,&B_p_AB);

  for (int i=0, off=0; i<E12_iter.num_blocks; i++) {
    read_block(&E12_iter,&B_p_AB);

    C_DGEMV('n',E12_iter.curr_size,noccA_*noccB_,1.0,&(B_p_AB.B_p_[0][0]),
      noccA_*noccB_,pAB[0],1,0.0,&(Y[off]),1);

    off += E12_iter.curr_size;
  }

  ex7 = -2.0*C_DDOT(ndf_+3,X,1,Y,1);

  B_p_AB.clear();

  Iterator E13_iter = get_iterator(mem_,&A_p_AA,&B_p_AB);

  for (int i=0, off=0; i<E13_iter.num_blocks; i++) {
    read_block(&E13_iter,&A_p_AA,&B_p_AB);
    for (int j=0; j<E13_iter.curr_size; j++) {
      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,A_p_AA.B_p_[j],noccA_,
        pAB[0],noccB_,0.0,xAB[0],noccB_);
      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,pAA[0],noccA_,
        xAB[0],noccB_,0.0,yAB[0],noccB_);
      ex7 += C_DDOT(noccA_*noccB_,yAB[0],1,B_p_AB.B_p_[j],1);
  }}
  
  A_p_AA.clear();
  B_p_AB.clear();

  Iterator E14_iter = get_iterator(mem_,&A_p_AA);

  for (int i=0, off=0; i<E14_iter.num_blocks; i++) {
    read_block(&E14_iter,&A_p_AA);

    C_DGEMV('n',E14_iter.curr_size,noccA_*noccA_,1.0,&(A_p_AA.B_p_[0][0]),
      noccA_*noccA_,pAA[0],1,0.0,&(X[off]),1);

    off += E14_iter.curr_size;
  }

  A_p_AA.clear();

  Iterator E15_iter = get_iterator(mem_,&B_p_BB);

  for (int i=0, off=0; i<E15_iter.num_blocks; i++) {
    read_block(&E15_iter,&B_p_BB);

    C_DGEMV('n',E15_iter.curr_size,noccB_*noccB_,1.0,&(B_p_BB.B_p_[0][0]),
      noccB_*noccB_,pBB[0],1,0.0,&(Y[off]),1);

    off += E15_iter.curr_size;
  }

  ex8 = -2.0*C_DDOT(ndf_+3,X,1,Y,1);

  B_p_BB.clear();

  Iterator E16_iter = get_iterator(mem_,&A_p_AB,&B_p_AB);

  for (int i=0, off=0; i<E16_iter.num_blocks; i++) {
    read_block(&E16_iter,&A_p_AB,&B_p_AB);
    for (int j=0; j<E16_iter.curr_size; j++) {
      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,pAA[0],noccA_,
        A_p_AB.B_p_[j],noccB_,0.0,xAB[0],noccB_);
      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,xAB[0],noccB_,
        pBB[0],noccB_,0.0,yAB[0],noccB_);
      ex8 += C_DDOT(noccA_*noccB_,yAB[0],1,B_p_AB.B_p_[j],1);
  }}

  A_p_AB.clear();
  B_p_AB.clear();

  Iterator E17_iter = get_iterator(mem_,&A_p_AB);

  for (int i=0, off=0; i<E17_iter.num_blocks; i++) {
    read_block(&E17_iter,&A_p_AB);

    C_DGEMV('n',E17_iter.curr_size,noccA_*noccB_,1.0,&(A_p_AB.B_p_[0][0]),
      noccA_*noccB_,pAB[0],1,0.0,&(X[off]),1);

    off += E17_iter.curr_size;
  }

  A_p_AB.done();

  Iterator E18_iter = get_iterator(mem_,&B_p_AB);

  for (int i=0, off=0; i<E18_iter.num_blocks; i++) {
    read_block(&E18_iter,&B_p_AB);

    C_DGEMV('n',E18_iter.curr_size,noccA_*noccB_,1.0,&(B_p_AB.B_p_[0][0]),
      noccA_*noccB_,pAB[0],1,0.0,&(Y[off]),1);

    off += E18_iter.curr_size;
  }

  ex9 = -2.0*C_DDOT(ndf_+3,X,1,Y,1);

  B_p_AB.done();

  Iterator E19_iter = get_iterator(mem_,&A_p_AA,&B_p_BB);

  for (int i=0, off=0; i<E19_iter.num_blocks; i++) {
    read_block(&E19_iter,&A_p_AA,&B_p_BB);
    for (int j=0; j<E19_iter.curr_size; j++) {
      C_DGEMM('N','N',noccA_,noccB_,noccA_,1.0,A_p_AA.B_p_[j],noccA_,
        pAB[0],noccB_,0.0,xAB[0],noccB_);
      C_DGEMM('N','N',noccA_,noccB_,noccB_,1.0,pAB[0],noccB_,
        B_p_BB.B_p_[j],noccB_,0.0,yAB[0],noccB_);
      ex9 += C_DDOT(noccA_*noccB_,xAB[0],1,yAB[0],1);
  }}

  A_p_AA.done();
  B_p_BB.done();

  free(X);
  free(Y);
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
