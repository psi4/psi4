#include "sapt0.h"
#include "sapt2.h"

namespace psi { namespace sapt {

void SAPT0::exch_disp20_n5()
{
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  arbs();
  v1();
  q12();

  double e_disp20 = 0.0;
  double v_1 = 0.0, q_12 = 0.0;

  bool A_in_core = false;
  bool B_in_core = false;

  long int mem = mem_;

  if (mem > 2L*aoccB_*nvirA_*(ndf_+3) + 2L*aoccA_*nvirB_*(ndf_+3) 
    + 1L*aoccA_*nvirA_*(ndf_) + 1L*aoccB_*nvirB_*(ndf_)) {
    A_in_core = true;
    B_in_core = true;
  }
  else if (mem > 2L*aoccA_*nvirB_*(ndf_+3) + 1L*aoccA_*nvirA_*(ndf_)
    + 2L*nvirA_*(ndf_+3) + 1L*nvirB_*(ndf_)) {
    A_in_core = true;
  }

  int A_chunk;
  int B_chunk;
  int A_blocks;
  int B_blocks;

  if (A_in_core && B_in_core) {
    A_chunk = aoccA_;
    B_chunk = aoccB_;
    A_blocks = 1;
    B_blocks = 1;
  }
  else if (A_in_core) {
    long int avail_mem = mem - (2L*aoccA_*nvirB_*(ndf_+3) 
      + 1L*aoccA_*nvirA_*(ndf_));
    A_chunk = aoccA_;
    A_blocks = 1;
    B_chunk = avail_mem / (2L*nvirA_*(ndf_+3) + 1L*nvirB_*ndf_);
    if (B_chunk > aoccB_) B_chunk = aoccB_;
    B_blocks = aoccB_ / B_chunk;
    if (aoccB_ % B_chunk) B_blocks++;
  }
  else if (mem > 2L*nvirB_*(ndf_+3) + 1L*nvirA_*(ndf_) + 2L*nvirA_*(ndf_+3) 
    + 1L*nvirB_*(ndf_)) {
    long int avail_mem = mem - (2L*nvirA_*(ndf_+3) + 1L*nvirB_*(ndf_));
    B_chunk = 1;
    B_blocks = aoccB_;
    A_chunk = mem / (2L*nvirB_*(ndf_+3) + 1L*nvirA_*ndf_);
    if (A_chunk > aoccA_) A_chunk = aoccA_;
    A_blocks = aoccA_ / A_chunk;
    if (aoccA_ % A_chunk) A_blocks++;
  }
  else {
    throw PsiException("Not enough memory", __FILE__,__LINE__);
  }

  if (debug_) {
    fprintf(outfile,"\n    A in core = %d\n",A_in_core);
    fprintf(outfile,"    B in core = %d\n",B_in_core);
    fprintf(outfile,"    A read in %3d blocks with length %d\n",
      A_blocks,A_chunk);
    fprintf(outfile,"    B read in %3d blocks with length %d\n\n",
      B_blocks,B_chunk);
  }

  double **tabRS = block_matrix(nthreads,nvirA_*nvirB_);
  double **vabRS = block_matrix(nthreads,nvirA_*nvirB_);

  double **T_AR = block_matrix(A_chunk*nvirA_,ndf_);
  double **T_BS = block_matrix(B_chunk*nvirB_,ndf_);
  double **V_BR = block_matrix(B_chunk*nvirA_,ndf_+3);
  double **V_AS = block_matrix(A_chunk*nvirB_,ndf_+3);
  double **Q_BR = block_matrix(B_chunk*nvirA_,ndf_+3);
  double **Q_AS = block_matrix(A_chunk*nvirB_,ndf_+3);

  if (A_in_core) {

    psio_->read_entry(PSIF_SAPT_TEMP,"AR RI Integrals",(char *)
      &(T_AR[0][0]),sizeof(double)*A_chunk*nvirA_*ndf_);
    psio_->read_entry(PSIF_SAPT_TEMP,"V1 AS RI Integrals",(char *)
      &(V_AS[0][0]),sizeof(double)*A_chunk*nvirB_*(ndf_+3));
    psio_->read_entry(PSIF_SAPT_TEMP,"Q12 AS RI Integrals",(char *)
      &(Q_AS[0][0]),sizeof(double)*A_chunk*nvirB_*(ndf_+3));

  }

  if (B_in_core) {

    psio_->read_entry(PSIF_SAPT_TEMP,"BS RI Integrals",(char *)
      &(T_BS[0][0]),sizeof(double)*B_chunk*nvirB_*ndf_);
    psio_->read_entry(PSIF_SAPT_TEMP,"V1 BR RI Integrals",(char *)
      &(V_BR[0][0]),sizeof(double)*B_chunk*nvirA_*(ndf_+3));
    psio_->read_entry(PSIF_SAPT_TEMP,"Q12 BR RI Integrals",(char *)
      &(Q_BR[0][0]),sizeof(double)*B_chunk*nvirA_*(ndf_+3));

  }

  psio_address next_T_AR = PSIO_ZERO;
  psio_address next_T_BS = PSIO_ZERO;
  psio_address next_V_BR = PSIO_ZERO;
  psio_address next_V_AS = PSIO_ZERO;
  psio_address next_Q_BR = PSIO_ZERO;
  psio_address next_Q_AS = PSIO_ZERO;

  for (int a=0,amax=0; a<A_blocks; a++) {

    int A_length = -amax;
    amax += A_chunk;
    if (amax > aoccA_) amax = aoccA_;
    A_length += amax;

    if (!A_in_core) {
      psio_->read(PSIF_SAPT_TEMP,"AR RI Integrals",(char *)
        &(T_AR[0][0]),sizeof(double)*A_length*nvirA_*ndf_,next_T_AR,
        &next_T_AR);
      psio_->read(PSIF_SAPT_TEMP,"V1 AS RI Integrals",(char *)
        &(V_AS[0][0]),sizeof(double)*A_length*nvirB_*(ndf_+3),next_V_AS,
        &next_V_AS);
      psio_->read(PSIF_SAPT_TEMP,"Q12 AS RI Integrals",(char *)
        &(Q_AS[0][0]),sizeof(double)*A_length*nvirB_*(ndf_+3),next_Q_AS,
        &next_Q_AS);
    }

    for (int b=0,bmax=0; b<B_blocks; b++) {

      int B_length = -bmax;
      bmax += B_chunk;
      if (bmax > aoccB_) bmax = aoccB_;
      B_length += bmax;

      if (!B_in_core) {
        psio_->read(PSIF_SAPT_TEMP,"BS RI Integrals",(char *)
          &(T_BS[0][0]),sizeof(double)*B_length*nvirB_*ndf_,next_T_BS,
          &next_T_BS);
        psio_->read(PSIF_SAPT_TEMP,"V1 BR RI Integrals",(char *)
          &(V_BR[0][0]),sizeof(double)*B_length*nvirA_*(ndf_+3),next_V_BR,
          &next_V_BR);
        psio_->read(PSIF_SAPT_TEMP,"Q12 BR RI Integrals",(char *)
          &(Q_BR[0][0]),sizeof(double)*B_length*nvirA_*(ndf_+3),next_Q_BR,
          &next_Q_BR);
      }

#pragma omp parallel
{
#pragma omp for private(rank) reduction(+:e_disp20,v_1,q_12)
      for (int ab=0; ab<A_length*B_length; ab++) {

#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        int arel = ab / B_length;
        int brel = ab % B_length;
        int aabs = amax - A_length + arel;
        int babs = bmax - B_length + brel;

        C_DGEMM('N','T',nvirA_,nvirB_,ndf_,1.0,T_AR[arel*nvirA_],ndf_,
          T_BS[brel*nvirB_],ndf_,0.0,tabRS[rank],nvirB_);

        for (int r=0,rs=0; r<nvirA_; r++) {
          for (int s=0; s<nvirB_; s++,rs++) {
            double denom =  evalsA_[aabs+foccA_] + evalsB_[babs+foccB_]
              - evalsA_[r+noccA_] - evalsB_[s+noccB_];
            double tval = tabRS[rank][rs];
            tabRS[rank][rs] /= denom;
            e_disp20 += 4.0*tval*tval/denom;
        }}

        C_DGEMM('N','T',nvirA_,nvirB_,ndf_+3,1.0,V_BR[brel*nvirA_],ndf_+3,
          V_AS[arel*nvirB_],ndf_+3,0.0,vabRS[rank],nvirB_);

        v_1 += C_DDOT(nvirA_*nvirB_,tabRS[rank],1,vabRS[rank],1);

        C_DGEMM('N','T',nvirA_,nvirB_,ndf_+3,1.0,Q_BR[brel*nvirA_],ndf_+3,
          Q_AS[arel*nvirB_],ndf_+3,0.0,vabRS[rank],nvirB_);
  
        q_12 += C_DDOT(nvirA_*nvirB_,tabRS[rank],1,vabRS[rank],1);

      }
}
    }
  
    next_T_BS = PSIO_ZERO;
    next_V_BR = PSIO_ZERO;
    next_Q_BR = PSIO_ZERO;

  }

  free_block(tabRS);
  free_block(vabRS);

  free_block(T_AR);
  free_block(T_BS);
  free_block(V_BR);
  free_block(V_AS);
  free_block(Q_BR);
  free_block(Q_AS);

  e_disp20_ = e_disp20;
  e_exch_disp20_ = -2.0*(v_1+q_12);

  if (print_) {
    fprintf(outfile,"    Disp20              = %18.12lf H\n",e_disp20_);
    fflush(outfile);
  }

  if (debug_) {
    fprintf(outfile,"\n    V1 + H2 + H4 + Q9   = %18.12lf H\n",v_1);
    fprintf(outfile,"    Q12                 = %18.12lf H\n",q_12);
    fflush(outfile);
  }
}

void SAPT0::exch_disp20_n4()
{
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  double h_1 = 0.0, h_2 = 0.0, h_3 = 0.0, h_4 = 0.0;
  double q_1 = 0.0, q_2 = 0.0, q_3 = 0.0, q_4 = 0.0; 
  double q_5 = 0.0, q_6 = 0.0, q_7 = 0.0, q_8 = 0.0;
  double q_10 = 0.0, q_11 = 0.0, q_13 = 0.0, q_14 = 0.0;

  theta_ar();
  theta_bs();
  if (debug_) test_theta();

  h_2 = h2();
	h_4 = h4();
  q_2 = q2();
  q_6 = q6();
  q_13 = q13();
  q_14 = q14();

  if (debug_) {
    fprintf(outfile,"    H2                  = %18.12lf H\n",h_2);
    fprintf(outfile,"    H4                  = %18.12lf H\n",h_4);
    fprintf(outfile,"    Q2                  = %18.12lf H\n",q_2);
    fprintf(outfile,"    Q6                  = %18.12lf H\n",q_6);
    fprintf(outfile,"    Q13                 = %18.12lf H\n",q_13);
    fprintf(outfile,"    Q14                 = %18.12lf H\n",q_14);
    fflush(outfile);
  }

  h1();
  h3();
  q1();
  q3();
  q5();
  q7();
  q10();
  q11();

  double **sAS = block_matrix(aoccA_,nvirB_);
  double **sRB = block_matrix(nvirA_,aoccB_);
  double **sAR = block_matrix(aoccA_,nvirA_);
  double **sBS = block_matrix(aoccB_,nvirB_);

  for (int a=0; a<aoccA_; a++) {
    C_DCOPY(nvirB_,&(sAB_[a+foccA_][noccB_]),1,sAS[a],1);
  }

  for (int r=0; r<nvirA_; r++) {
    C_DCOPY(aoccB_,&(sAB_[r+noccA_][foccB_]),1,sRB[r],1);
  } 

  C_DGEMM('N','T',aoccA_,nvirA_,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(sAB_[noccA_][0]),nmoB_,0.0,sAR[0],nvirA_);

  C_DGEMM('T','N',aoccB_,nvirB_,noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(sAB_[0][noccB_]),nmoB_,0.0,sBS[0],nvirB_);

  double **H1_RB = block_matrix(nvirA_,aoccB_);
  double **H3_AS = block_matrix(aoccA_,nvirB_);
  double **Q1_AS = block_matrix(aoccA_,nvirB_);
  double **Q3_AS = block_matrix(aoccA_,nvirB_);
  double **Q4_BS = block_matrix(aoccB_,nvirB_);
  double **Q5_RB = block_matrix(nvirA_,aoccB_);
  double **Q7_RB = block_matrix(nvirA_,aoccB_);
  double **Q8_AR = block_matrix(aoccA_,nvirA_);
  double **Q10_AS = block_matrix(aoccA_,nvirB_);
  double **Q11_RB = block_matrix(nvirA_,aoccB_);

  psio_->read_entry(PSIF_SAPT_TEMP,"H1 RB Array",(char *) &(H1_RB[0][0]),
    sizeof(double)*nvirA_*aoccB_);

  psio_->read_entry(PSIF_SAPT_TEMP,"H3 AS Array",(char *) &(H3_AS[0][0]),
    sizeof(double)*aoccA_*nvirB_);

  psio_->read_entry(PSIF_SAPT_TEMP,"Q1 AS Array",(char *) &(Q1_AS[0][0]),
    sizeof(double)*aoccA_*nvirB_);

  psio_->read_entry(PSIF_SAPT_TEMP,"Q3 AS Array",(char *) &(Q3_AS[0][0]),
    sizeof(double)*aoccA_*nvirB_);

  psio_->read_entry(PSIF_SAPT_TEMP,"Q4 BS Array",(char *) &(Q4_BS[0][0]),
    sizeof(double)*aoccB_*nvirB_);

  psio_->read_entry(PSIF_SAPT_TEMP,"Q5 RB Array",(char *) &(Q5_RB[0][0]),
    sizeof(double)*nvirA_*aoccB_);

  psio_->read_entry(PSIF_SAPT_TEMP,"Q7 RB Array",(char *) &(Q7_RB[0][0]),
    sizeof(double)*nvirA_*aoccB_);

  psio_->read_entry(PSIF_SAPT_TEMP,"Q8 AR Array",(char *) &(Q8_AR[0][0]),
    sizeof(double)*aoccA_*nvirA_);

  psio_->read_entry(PSIF_SAPT_TEMP,"Q10 AS Array",(char *) &(Q10_AS[0][0]),
    sizeof(double)*aoccA_*nvirB_);

  psio_->read_entry(PSIF_SAPT_TEMP,"Q11 RB Array",(char *) &(Q11_RB[0][0]),
    sizeof(double)*nvirA_*aoccB_);

  SAPTDFInts B_p_AR = set_act_C_AR();
  SAPTDFInts B_p_BS = set_act_C_BS();
  Iterator B_ARBS_iter = get_iterator(mem_/2,&B_p_AR,&B_p_BS);

  double **xAB = block_matrix(nthreads,aoccA_*aoccB_);
  double **yAB = block_matrix(nthreads,aoccA_*aoccB_);

  double **T_p_AR = block_matrix(B_ARBS_iter.block_size[0],aoccA_*nvirA_);
  double **T_p_BS = block_matrix(B_ARBS_iter.block_size[0],aoccB_*nvirB_);

  for (int j=0; j<B_ARBS_iter.num_blocks; j++) {
    read_block(&B_ARBS_iter,&B_p_AR,&B_p_BS);

    for (int i=0; i<nvec_; i++) {

      C_DCOPY((long int) B_ARBS_iter.block_size[j]*aoccA_*nvirA_,
        B_p_AR.B_p_[0],1,T_p_AR[0],1);
      C_DCOPY((long int) B_ARBS_iter.block_size[j]*aoccB_*nvirB_,
        B_p_BS.B_p_[0],1,T_p_BS[0],1);

#pragma omp parallel
{
#pragma omp for 
        for (int ar=0; ar<aoccA_*nvirA_; ar++) {
          double scale = dAR_[i][ar];
          C_DSCAL(B_ARBS_iter.curr_size,scale,&(T_p_AR[0][ar]),aoccA_*nvirA_);
        }

#pragma omp for 
        for (int bs=0; bs<aoccB_*nvirB_; bs++) {
          double scale = dBS_[i][bs];
          C_DSCAL(B_ARBS_iter.curr_size,scale,&(T_p_BS[0][bs]),aoccB_*nvirB_);
        }

#pragma omp for private(rank) reduction(+:h_1)
      for (int p=0; p<B_ARBS_iter.curr_size; p++) {

#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR[p],nvirA_,
          H1_RB[0],aoccB_,0.0,xAB[rank],aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,sAS[0],nvirB_,
          T_p_BS[p],nvirB_,0.0,yAB[rank],aoccB_);
        h_1 += C_DDOT(aoccA_*aoccB_,xAB[rank],1,yAB[rank],1);
      }

#pragma omp for private(rank) reduction(+:h_3)
      for (int p=0; p<B_ARBS_iter.curr_size; p++) {
        
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR[p],nvirA_,
          sRB[0],aoccB_,0.0,xAB[rank],aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,H3_AS[0],nvirB_,
          T_p_BS[p],nvirB_,0.0,yAB[rank],aoccB_);
        h_3 += C_DDOT(aoccA_*aoccB_,xAB[rank],1,yAB[rank],1);
      }

#pragma omp for private(rank) reduction(+:q_1)
      for (int p=0; p<B_ARBS_iter.curr_size; p++) {

#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR[p],nvirA_,
          sRB[0],aoccB_,0.0,xAB[rank],aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,Q1_AS[0],nvirB_,
          T_p_BS[p],nvirB_,0.0,yAB[rank],aoccB_);
        q_1 += C_DDOT(aoccA_*aoccB_,xAB[rank],1,yAB[rank],1);
      } 

#pragma omp for private(rank) reduction(+:q_3)
      for (int p=0; p<B_ARBS_iter.curr_size; p++) {

#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR[p],nvirA_,
          sRB[0],aoccB_,0.0,xAB[rank],aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,Q3_AS[0],nvirB_,
          T_p_BS[p],nvirB_,0.0,yAB[rank],aoccB_);
        q_3 -= 2.0*C_DDOT(aoccA_*aoccB_,xAB[rank],1,yAB[rank],1);
      }

#pragma omp for reduction(+:q_4)
      for (int p=0; p<B_ARBS_iter.curr_size; p++) {
        q_4 += 4.0*C_DDOT(aoccA_*nvirA_,T_p_AR[p],1,sAR[0],1)*
          C_DDOT(aoccB_*nvirB_,T_p_BS[p],1,Q4_BS[0],1);
      }

#pragma omp for private(rank) reduction(+:q_5)
      for (int p=0; p<B_ARBS_iter.curr_size; p++) {

#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR[p],nvirA_,
          Q5_RB[0],aoccB_,0.0,xAB[rank],aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,sAS[0],nvirB_,
          T_p_BS[p],nvirB_,0.0,yAB[rank],aoccB_);
        q_5 += C_DDOT(aoccA_*aoccB_,xAB[rank],1,yAB[rank],1);
      }

#pragma omp for private(rank) reduction(+:q_7)
      for (int p=0; p<B_ARBS_iter.curr_size; p++) { 
        
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR[p],nvirA_,
          Q7_RB[0],aoccB_,0.0,xAB[rank],aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,sAS[0],nvirB_,
          T_p_BS[p],nvirB_,0.0,yAB[rank],aoccB_);
        q_7 -= 2.0*C_DDOT(aoccA_*aoccB_,xAB[rank],1,yAB[rank],1);
      }

#pragma omp for reduction(+:q_8)
      for (int p=0; p<B_ARBS_iter.curr_size; p++) {
        q_8 += 4.0*C_DDOT(aoccA_*nvirA_,T_p_AR[p],1,Q8_AR[0],1)*
          C_DDOT(aoccB_*nvirB_,T_p_BS[p],1,sBS[0],1);
      }

#pragma omp for private(rank) reduction(+:q_10)
      for (int p=0; p<B_ARBS_iter.curr_size; p++) {

#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR[p],nvirA_,
          sRB[0],aoccB_,0.0,xAB[rank],aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,Q10_AS[0],nvirB_,
          T_p_BS[p],nvirB_,0.0,yAB[rank],aoccB_);
        q_10 -= 2.0*C_DDOT(aoccA_*aoccB_,xAB[rank],1,yAB[rank],1);
      }

#pragma omp for private(rank) reduction(+:q_11)
      for (int p=0; p<B_ARBS_iter.curr_size; p++) {

#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR[p],nvirA_,
          Q11_RB[0],aoccB_,0.0,xAB[rank],aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,sAS[0],nvirB_,
          T_p_BS[p],nvirB_,0.0,yAB[rank],aoccB_);
        q_11 -= 2.0*C_DDOT(aoccA_*aoccB_,xAB[rank],1,yAB[rank],1);
      }
}

    }
  }

  if (debug_) {
    fprintf(outfile,"    H1                  = %18.12lf H\n",h_1);
    fprintf(outfile,"    H3                  = %18.12lf H\n",h_3);
    fprintf(outfile,"    Q1                  = %18.12lf H\n",q_1);
    fprintf(outfile,"    Q3                  = %18.12lf H\n",q_3);
    fprintf(outfile,"    Q4                  = %18.12lf H\n",q_4);
    fprintf(outfile,"    Q5                  = %18.12lf H\n",q_5);
    fprintf(outfile,"    Q7                  = %18.12lf H\n",q_7);
    fprintf(outfile,"    Q8                  = %18.12lf H\n",q_8);
    fprintf(outfile,"    Q10                 = %18.12lf H\n",q_10);
    fprintf(outfile,"    Q11                 = %18.12lf H\n\n",q_11);
  }

  free_block(sAS);
  free_block(sRB);
  free_block(sAR);
  free_block(sBS);
  free_block(H1_RB);
  free_block(H3_AS);
  free_block(Q1_AS);
  free_block(Q3_AS);
  free_block(Q4_BS);
  free_block(Q5_RB);
  free_block(Q7_RB);
  free_block(Q8_AR);
  free_block(Q10_AS);
  free_block(Q11_RB);

  free_block(xAB);
  free_block(yAB);

  B_p_AR.done();
  B_p_BS.done();

  free_block(T_p_AR);
  free_block(T_p_BS);

  e_exch_disp20_ += 2.0*(h_1+h_2+h_3+h_4+q_1+q_2+q_3+q_4+q_5+q_6+q_7+q_8+
    q_10+q_11+q_13+q_14);

  if (print_) {
    fprintf(outfile,"    Exch-Disp20         = %18.12lf H\n",e_exch_disp20_);
    fflush(outfile);
  }
}

void SAPT0::theta_ar() 
{
  long int avail_mem = mem_ - (long int) nvec_*ndf_*(ndf_+3);

  if ((long int) 3*aoccB_*nvirB_ > avail_mem) 
    throw PsiException("Not enough memory", __FILE__,__LINE__);

  SAPTDFInts B_p_BS = set_act_B_BS();
  Iterator B_BS_iter = get_iterator(avail_mem/3,&B_p_BS);

  SAPTDFInts C_p_BS = set_act_C_BS();
  Iterator C_BS_iter = get_iterator(avail_mem/3,&C_p_BS);

  double **yPQ = block_matrix(nvec_,(long int) ndf_*(ndf_+3));
  double **T_p_BS = block_matrix(C_BS_iter.block_size[0],aoccB_*nvirB_);

  for (int j=0,offB=0; j<B_BS_iter.num_blocks; j++) {
    read_block(&B_BS_iter,&B_p_BS);

    for (int k=0,offC=0; k<C_BS_iter.num_blocks; k++) {
      read_block(&C_BS_iter,&C_p_BS);

      for (int i=0; i<nvec_; i++) {

        C_DCOPY((long int) C_BS_iter.block_size[k]*aoccB_*nvirB_,
          C_p_BS.B_p_[0],1,T_p_BS[0],1);

#pragma omp parallel
{
#pragma omp for 
        for (int bs=0; bs<aoccB_*nvirB_; bs++) {
          double scale = dBS_[i][bs];
          C_DSCAL(C_BS_iter.curr_size,scale,&(T_p_BS[0][bs]),
            aoccB_*nvirB_);
        }
}
        C_DGEMM('N','T',C_BS_iter.curr_size,B_BS_iter.curr_size,
          aoccB_*nvirB_,1.0,T_p_BS[0],aoccB_*nvirB_,B_p_BS.B_p_[0],
          aoccB_*nvirB_,0.0,&(yPQ[i][(long int) offC*(ndf_+3)+offB]),
          ndf_+3);

      }
      offC += C_BS_iter.curr_size;
    }
    C_p_BS.rewind();
    C_BS_iter.rewind();
    offB += B_BS_iter.curr_size;
  }

  B_p_BS.done();
  C_p_BS.done();

  free_block(T_p_BS);

  bool in_core = false;

  if (avail_mem > 1L*aoccA_*nvirA_*(ndf_+3) + 2L*aoccA_*nvirA_*ndf_) {
    in_core = true;
  }
  else if (avail_mem > 2L*nvirA_*ndf_ + 1L*nvirA_*(ndf_+3)) {
    in_core = false;
  }
  else {
    throw PsiException("Not enough memory", __FILE__,__LINE__);
  }

  int chunk;
  int blocks;

  if (in_core) {
    chunk = aoccA_;
    blocks = 1;
  }
  else {
    chunk = avail_mem / (2L*nvirA_*ndf_ + 1L*nvirA_*(ndf_+3));
    if (chunk > aoccA_) chunk = aoccA_;
    blocks = aoccA_ / chunk;
    if (aoccA_ % chunk) blocks++;
  }

  double **B_AR = block_matrix(chunk*nvirA_,ndf_);
  double **L_AR = block_matrix(chunk*nvirA_,ndf_);
  double **T_AR = block_matrix(chunk*nvirA_,ndf_+3);
  double *temp = init_array(chunk*nvirA_);

  if (in_core) {
    psio_->read_entry(PSIF_SAPT_TEMP,"AR RI Integrals",(char *)
      &(B_AR[0][0]),sizeof(double)*chunk*nvirA_*ndf_);
  }

  psio_address next_B_AR = PSIO_ZERO;
  psio_address next_T_AR = PSIO_ZERO;

  zero_disk(PSIF_SAPT_TEMP,"Theta AR Intermediate",ndf_+3,aoccA_*nvirA_);

  for (int a=0,amax=0; a<blocks; a++) {

    int length = -amax;
    int amin = amax;
    amax += chunk;
    if (amax > aoccA_) amax = aoccA_;
    length += amax;

    if (!in_core) {
      psio_->read(PSIF_SAPT_TEMP,"AR RI Integrals",(char *)
        &(B_AR[0][0]),sizeof(double)*length*nvirA_*ndf_,next_B_AR,
        &next_B_AR);
    }

    memset(&(T_AR[0][0]),'\0',sizeof(double)*length*nvirA_*(ndf_+3));

    for (int i=0; i<nvec_; i++) {

      C_DCOPY((long int) ndf_*length*nvirA_,B_AR[0],1,L_AR[0],1);

#pragma omp parallel
{
#pragma omp for 
      for (int ar=amin*nvirA_; ar<amax*nvirA_; ar++) {
          double scale = dAR_[i][ar];
          C_DSCAL(ndf_,scale,L_AR[ar-amin*nvirA_],1);
        }
}

      C_DGEMM('N','N',length*nvirA_,ndf_+3,ndf_,1.0,L_AR[0],ndf_,yPQ[i],ndf_+3,
        1.0,T_AR[0],ndf_+3);
 
    }
    for (int P=0; P<ndf_+3; P++) {
      next_T_AR = psio_get_address(PSIO_ZERO,sizeof(double)*P*aoccA_*nvirA_
        +sizeof(double)*amin*nvirA_);
      C_DCOPY((long int) length*nvirA_,&(T_AR[0][P]),ndf_+3,temp,1);
      psio_->write(PSIF_SAPT_TEMP,"Theta AR Intermediate",(char *)
        &(temp[0]),sizeof(double)*length*nvirA_,next_T_AR,&next_T_AR);
    }
  }

  free_block(B_AR);
  free_block(L_AR);
  free_block(T_AR);
  free(temp);

  if (debug_)
    psio_->write_entry(PSIF_SAPT_TEMP,"Y PQ Intermediate",(char *)
      &(yPQ[0][0]),sizeof(double)*nvec_*ndf_*(ndf_+3));

  free_block(yPQ);
}

void SAPT0::theta_bs()
{
  long int avail_mem = mem_ - (long int) nvec_*ndf_*(ndf_+3);

  if ((long int) 3*aoccA_*nvirA_ > avail_mem) 
    throw PsiException("Not enough memory", __FILE__,__LINE__);

  SAPTDFInts B_p_AR = set_act_A_AR();
  Iterator B_AR_iter = get_iterator(avail_mem/3,&B_p_AR);

  SAPTDFInts C_p_AR = set_act_C_AR();
  Iterator C_AR_iter = get_iterator(avail_mem/3,&C_p_AR);

  double **xPQ = block_matrix(nvec_,(long int) ndf_*(ndf_+3));
  double **T_p_AR = block_matrix(C_AR_iter.block_size[0],aoccA_*nvirA_);

  for (int j=0,offB=0; j<B_AR_iter.num_blocks; j++) {
    read_block(&B_AR_iter,&B_p_AR);

    for (int k=0,offC=0; k<C_AR_iter.num_blocks; k++) {
      read_block(&C_AR_iter,&C_p_AR);

      for (int i=0; i<nvec_; i++) {

        C_DCOPY((long int) C_AR_iter.block_size[k]*aoccA_*nvirA_,
          C_p_AR.B_p_[0],1,T_p_AR[0],1);

#pragma omp parallel
{
#pragma omp for 
        for (int ar=0; ar<aoccA_*nvirA_; ar++) {
          double scale = dAR_[i][ar];
          C_DSCAL(C_AR_iter.curr_size,scale,&(T_p_AR[0][ar]),
            aoccA_*nvirA_);
        }
}
        C_DGEMM('N','T',C_AR_iter.curr_size,B_AR_iter.curr_size,
          aoccA_*nvirA_,1.0,T_p_AR[0],aoccA_*nvirA_,B_p_AR.B_p_[0],
          aoccA_*nvirA_,0.0,&(xPQ[i][(long int) offC*(ndf_+3)+offB]),
          ndf_+3);

      }
      offC += C_AR_iter.curr_size;
    }
    C_p_AR.rewind();
    C_AR_iter.rewind();
    offB += B_AR_iter.curr_size;
  }

  B_p_AR.done();
  C_p_AR.done();

  free_block(T_p_AR);

  bool in_core = false;

  if (avail_mem > 1L*aoccB_*nvirB_*(ndf_+3) + 2L*aoccB_*nvirB_*ndf_) {
    in_core = true;
  }
  else if (avail_mem > 2L*nvirB_*ndf_ + 1L*nvirB_*(ndf_+3)) {
    in_core = false;
  }
  else {
    throw PsiException("Not enough memory", __FILE__,__LINE__);
  }

  int chunk;
  int blocks;

  if (in_core) {
    chunk = aoccB_;
    blocks = 1;
  }
  else {
    chunk = avail_mem / (2L*nvirB_*ndf_ + 1L*nvirB_*(ndf_+3));
    if (chunk > aoccB_) chunk = aoccB_;
    blocks = aoccB_ / chunk;
    if (aoccB_ % chunk) blocks++;
  }

  double **B_BS = block_matrix(chunk*nvirB_,ndf_);
  double **L_BS = block_matrix(chunk*nvirB_,ndf_);
  double **T_BS = block_matrix(chunk*nvirB_,ndf_+3);
  double *temp = init_array(chunk*nvirB_);

  if (in_core) {
    psio_->read_entry(PSIF_SAPT_TEMP,"BS RI Integrals",(char *)
      &(B_BS[0][0]),sizeof(double)*chunk*nvirB_*ndf_);
  }

  psio_address next_B_BS = PSIO_ZERO;
  psio_address next_T_BS = PSIO_ZERO;

  zero_disk(PSIF_SAPT_TEMP,"Theta BS Intermediate",ndf_+3,aoccB_*nvirB_);

  for (int b=0,bmax=0; b<blocks; b++) {

    int length = -bmax;
    int bmin = bmax;
    bmax += chunk;
    if (bmax > aoccB_) bmax = aoccB_;
    length += bmax;

    if (!in_core) {
      psio_->read(PSIF_SAPT_TEMP,"BS RI Integrals",(char *)
        &(B_BS[0][0]),sizeof(double)*length*nvirB_*ndf_,next_B_BS,
        &next_B_BS);
    }

    memset(&(T_BS[0][0]),'\0',sizeof(double)*length*nvirB_*(ndf_+3));

    for (int i=0; i<nvec_; i++) {

      C_DCOPY((long int) ndf_*length*nvirB_,B_BS[0],1,L_BS[0],1);

#pragma omp parallel
{
#pragma omp for 
      for (int bs=bmin*nvirB_; bs<bmax*nvirB_; bs++) {
          double scale = dBS_[i][bs];
          C_DSCAL(ndf_,scale,L_BS[bs-bmin*nvirB_],1);
        }
}

      C_DGEMM('N','N',length*nvirB_,ndf_+3,ndf_,1.0,L_BS[0],ndf_,xPQ[i],ndf_+3,
        1.0,T_BS[0],ndf_+3);
 
    }
    for (int P=0; P<ndf_+3; P++) {
      next_T_BS = psio_get_address(PSIO_ZERO,sizeof(double)*P*aoccB_*nvirB_
        +sizeof(double)*bmin*nvirB_);
      C_DCOPY((long int) length*nvirB_,&(T_BS[0][P]),ndf_+3,temp,1);
      psio_->write(PSIF_SAPT_TEMP,"Theta BS Intermediate",(char *)
        &(temp[0]),sizeof(double)*length*nvirB_,next_T_BS,&next_T_BS);
    }
  }

  free_block(B_BS);
  free_block(L_BS);
  free_block(T_BS);
  free(temp);

  if (debug_)
    psio_->write_entry(PSIF_SAPT_TEMP,"X PQ Intermediate",(char *)
      &(xPQ[0][0]),sizeof(double)*nvec_*ndf_*(ndf_+3));

  free_block(xPQ);
}

void SAPT0::test_theta()
{
/*
  double **B_AR = block_matrix(aoccA_*nvirA_,ndf_);
  double **T_AR = block_matrix(aoccA_*nvirA_,ndf_+3);

  psio_->read_entry(PSIF_SAPT_TEMP,"AR RI Integrals",(char *)
    &(B_AR[0][0]),sizeof(double)*aoccA_*nvirA_*ndf_);

  psio_->read_entry(PSIF_SAPT_TEMP,"Theta AR Intermediate",(char *)
    &(T_AR[0][0]),sizeof(double)*aoccA_*nvirA_*(ndf_+3));

  double disp_ar = 0.0;	

  for (int ar=0; ar<aoccA_*nvirA_; ar++) {
    disp_ar += -4.0*C_DDOT(ndf_,B_AR[ar],1,T_AR[ar],1);
  }

  free_block(B_AR);
  free_block(T_AR);

  double **B_BS = block_matrix(aoccB_*nvirB_,ndf_);
  double **T_BS = block_matrix(aoccB_*nvirB_,ndf_+3);

  psio_->read_entry(PSIF_SAPT_TEMP,"BS RI Integrals",(char *)
    &(B_BS[0][0]),sizeof(double)*aoccB_*nvirB_*ndf_);

  psio_->read_entry(PSIF_SAPT_TEMP,"Theta BS Intermediate",(char *)
    &(T_BS[0][0]),sizeof(double)*aoccB_*nvirB_*(ndf_+3));

  double disp_bs = 0.0;

  for (int bs=0; bs<aoccB_*nvirB_; bs++) {
    disp_bs += -4.0*C_DDOT(ndf_,B_BS[bs],1,T_BS[bs],1);
  }

  free_block(B_BS);
  free_block(T_BS);
*/
  double **xPQ = block_matrix(nvec_,ndf_*(ndf_+3));
  double **yPQ = block_matrix(nvec_,ndf_*(ndf_+3));

  psio_->read_entry(PSIF_SAPT_TEMP,"X PQ Intermediate",(char *)
    &(xPQ[0][0]),sizeof(double)*nvec_*ndf_*(ndf_+3));
  psio_->read_entry(PSIF_SAPT_TEMP,"Y PQ Intermediate",(char *)
    &(yPQ[0][0]),sizeof(double)*nvec_*ndf_*(ndf_+3));

  double disp_xy = -4.0*C_DDOT((long int) nvec_*ndf_*(ndf_+3),xPQ[0],1,
    yPQ[0],1);

  free_block(xPQ); 
  free_block(yPQ);

  fprintf(outfile,"    Disp20 (xPQ yPQ)    = %18.12lf H\n",disp_xy);
//fprintf(outfile,"    Disp20 (Theta AR)   = %18.12lf H\n",disp_ar);
//fprintf(outfile,"    Disp20 (Theta BS)   = %18.12lf H\n",disp_bs);
}

void SAPT0::arbs()
{
  zero_disk(PSIF_SAPT_TEMP,"AR RI Integrals",aoccA_*nvirA_,ndf_);
  zero_disk(PSIF_SAPT_TEMP,"BS RI Integrals",aoccB_*nvirB_,ndf_);

  SAPTDFInts T_p_AR = set_act_C_AR();
  Iterator T_AR_iter = get_iterator(mem_/2,&T_p_AR);

  double **X_AR_p = block_matrix(aoccA_*nvirA_,T_AR_iter.block_size[0]);

  psio_address next_T_AR = PSIO_ZERO;

  for (int i=0,off=0; i<T_AR_iter.num_blocks; i++) {
    read_block(&T_AR_iter,&T_p_AR);
    for (int p=0; p<T_AR_iter.curr_size; p++) {
      C_DCOPY((long int) aoccA_*nvirA_,&(T_p_AR.B_p_[p][0]),1,&(X_AR_p[0][p]),
        T_AR_iter.block_size[0]);
    }

    int skip = ndf_ - T_AR_iter.curr_size;

    next_T_AR = psio_get_address(PSIO_ZERO,sizeof(double)*off);

    for (int ar=0; ar<aoccA_*nvirA_; ar++) {
      psio_->write(PSIF_SAPT_TEMP,"AR RI Integrals",(char *)
        &(X_AR_p[ar][0]),sizeof(double)*T_AR_iter.curr_size,next_T_AR,
        &next_T_AR);
      next_T_AR = psio_get_address(next_T_AR,sizeof(double)*skip);
    }

    off += T_AR_iter.curr_size;
  }

  free_block(X_AR_p);
  T_p_AR.done();

  SAPTDFInts T_p_BS = set_act_C_BS();
  Iterator T_BS_iter = get_iterator(mem_/2,&T_p_BS);

  double **X_BS_p = block_matrix(aoccB_*nvirB_,T_BS_iter.block_size[0]);

  psio_address next_T_BS = PSIO_ZERO;

  for (int i=0,off=0; i<T_BS_iter.num_blocks; i++) {
    read_block(&T_BS_iter,&T_p_BS);
    for (int p=0; p<T_BS_iter.curr_size; p++) {
      C_DCOPY((long int) aoccB_*nvirB_,&(T_p_BS.B_p_[p][0]),1,&(X_BS_p[0][p]),
        T_BS_iter.block_size[0]);
    }

    int skip = ndf_ - T_BS_iter.curr_size; 

    next_T_BS = psio_get_address(PSIO_ZERO,sizeof(double)*off);

    for (int bs=0; bs<aoccB_*nvirB_; bs++) {
      psio_->write(PSIF_SAPT_TEMP,"BS RI Integrals",(char *)
        &(X_BS_p[bs][0]),sizeof(double)*T_BS_iter.curr_size,next_T_BS,
        &next_T_BS);
      next_T_BS = psio_get_address(next_T_BS,sizeof(double)*skip);
    }
  
    off += T_BS_iter.curr_size;
  }
  
  free_block(X_BS_p);
  T_p_BS.done();
}

void SAPT0::v1()
{
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  zero_disk(PSIF_SAPT_TEMP,"V1 AS RI Integrals",aoccA_*nvirB_,ndf_+3);
  zero_disk(PSIF_SAPT_TEMP,"V1 BR RI Integrals",aoccB_*nvirA_,ndf_+3);

  SAPTDFInts A_p_AA = set_A_AA();
  SAPTDFInts A_p_AS = set_A_AS();
  Iterator AS_iter = get_iterator(mem_/2,&A_p_AA,&A_p_AS);
  double **xAS = block_matrix(nthreads,aoccA_*nvirB_);
  double **X_AS_p = block_matrix(aoccA_*nvirB_,AS_iter.block_size[0]);

  psio_address next_A_AS = PSIO_ZERO;

  for (int i=0,off=0; i<AS_iter.num_blocks; i++) {
    read_block(&AS_iter,&A_p_AA,&A_p_AS);
#pragma omp parallel
{
#pragma omp for private(rank)
    for (int p=0; p<AS_iter.curr_size; p++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',aoccA_,nvirB_,noccA_,-1.0,
        &(A_p_AA.B_p_[p][foccA_*noccA_]),noccA_,&(sAB_[0][noccB_]),nmoB_,
        0.0,xAS[rank],nvirB_);
      C_DCOPY((long int) aoccA_*nvirB_,&(A_p_AS.B_p_[p][foccA_*nvirB_]),1,
        &(X_AS_p[0][p]),AS_iter.block_size[0]);
      C_DAXPY((long int) aoccA_*nvirB_,1.0,xAS[rank],1,&(X_AS_p[0][p]),
        AS_iter.block_size[0]); 
    }
}
    int skip = ndf_ + 3 - AS_iter.curr_size;

    next_A_AS = psio_get_address(PSIO_ZERO,sizeof(double)*off);

    for (int as=0; as<aoccA_*nvirB_; as++) {
      psio_->write(PSIF_SAPT_TEMP,"V1 AS RI Integrals",(char *)
        &(X_AS_p[as][0]),sizeof(double)*AS_iter.curr_size,next_A_AS,
        &next_A_AS);
      next_A_AS = psio_get_address(next_A_AS,sizeof(double)*skip);
    }

    off += AS_iter.curr_size;
  }

  free_block(xAS);
  free_block(X_AS_p);
  A_p_AA.done();
  A_p_AS.done();

  SAPTDFInts B_p_BB = set_B_BB();
  SAPTDFInts B_p_RB = set_B_RB();
  Iterator BR_iter = get_iterator(mem_/2,&B_p_BB,&B_p_RB);
  double **xBR = block_matrix(nthreads,aoccB_*nvirA_);
  double **X_BR_p = block_matrix(aoccB_*nvirA_,BR_iter.block_size[0]);

  psio_address next_B_BR = PSIO_ZERO;

  for (int i=0,off=0; i<BR_iter.num_blocks; i++) {
    read_block(&BR_iter,&B_p_BB,&B_p_RB);
#pragma omp parallel
{
#pragma omp for private(rank)
    for (int p=0; p<BR_iter.curr_size; p++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','T',aoccB_,nvirA_,noccB_,-1.0,
        &(B_p_BB.B_p_[p][foccB_*noccB_]),noccB_,&(sAB_[noccA_][0]),nmoB_,
        0.0,xBR[rank],nvirA_);
      for (int b=0; b<aoccB_; b++) {
        C_DCOPY(nvirA_,&(B_p_RB.B_p_[p][b+foccB_]),noccB_,
          &(X_BR_p[b*nvirA_][p]),BR_iter.block_size[0]);
      }
      C_DAXPY(aoccB_*nvirA_,1.0,xBR[rank],1,&(X_BR_p[0][p]),
        BR_iter.block_size[0]);
    }
}
    int skip = ndf_ + 3 - BR_iter.curr_size;

    next_B_BR = psio_get_address(PSIO_ZERO,sizeof(double)*off);

    for (int br=0; br<aoccB_*nvirA_; br++) {
      psio_->write(PSIF_SAPT_TEMP,"V1 BR RI Integrals",(char *)
        &(X_BR_p[br][0]),sizeof(double)*BR_iter.curr_size,next_B_BR,
        &next_B_BR);
      next_B_BR = psio_get_address(next_B_BR,sizeof(double)*skip);
    }

    off += BR_iter.curr_size;
  }

  free_block(xBR);
  free_block(X_BR_p);
  B_p_BB.done();
  B_p_RB.done();
}

void SAPT0::h1()
{
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  SAPTDFInts B_p_RB = set_B_RB();
  Iterator RB_iter = get_iterator(mem_,&B_p_RB);

  double *tRB = init_array(nvirA_*noccB_);

  for (int i=0,off=0; i<RB_iter.num_blocks; i++) {
    read_block(&RB_iter,&B_p_RB);
    C_DGEMV('t',RB_iter.curr_size,nvirA_*noccB_,2.0,&(B_p_RB.B_p_[0][0]),
      nvirA_*noccB_,&(diagAA_[off]),1,1.0,tRB,1);
    off += RB_iter.curr_size;
  }

  B_p_RB.done();

  double **xRB = block_matrix(nthreads,nvirA_*aoccB_);

  for (int r=0; r<nvirA_; r++) {
    C_DCOPY(aoccB_,&(tRB[r*noccB_ + foccB_]),1,&(xRB[0][r*aoccB_]),1);
  }

  free(tRB);

  SAPTDFInts A_p_AR = set_A_AR();
  SAPTDFInts B_p_AB = set_B_AB();
  Iterator ARAB_iter = get_iterator(mem_,&A_p_AR,&B_p_AB);

  for (int i=0; i<ARAB_iter.num_blocks; i++) {
    read_block(&ARAB_iter,&A_p_AR,&B_p_AB);
#pragma omp parallel
{
#pragma omp for private(rank)
    for (int j=0; j<ARAB_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('T','N',nvirA_,aoccB_,noccA_,-1.0,A_p_AR.B_p_[j],nvirA_,
        &(B_p_AB.B_p_[j][foccB_]),noccB_,1.0,xRB[rank],aoccB_);
    }
}
  }

  for (int n=1; n<nthreads; n++)
    C_DAXPY(nvirA_*aoccB_,1.0,xRB[n],1,xRB[0],1);

  A_p_AR.done();
  B_p_AB.done();

  psio_->write_entry(PSIF_SAPT_TEMP,"H1 RB Array",(char *) &(xRB[0][0]),
    sizeof(double)*nvirA_*aoccB_);

  free_block(xRB);
}

double SAPT0::h2()
{
  double h_2 = 0.0;

  SAPTDFInts B_p_AB = set_B_AB();
  Iterator AB_iter = get_iterator(mem_,&B_p_AB);

  double **xBS = block_matrix(aoccB_,nvirB_);
  double **yBS = block_matrix(aoccB_,nvirB_);
  psio_address next_BS = PSIO_ZERO;
  psio_address next_T_BS = PSIO_ZERO;

  for (int i=0; i<AB_iter.num_blocks; i++) {
    read_block(&AB_iter,&B_p_AB);
    for (int p=0; p<AB_iter.curr_size; p++) {
      C_DGEMM('T','N',aoccB_,nvirB_,noccA_,1.0,&(B_p_AB.B_p_[p][foccB_]),
        noccB_,&(sAB_[0][noccB_]),nmoB_,0.0,xBS[0],nvirB_);
      psio_->write(PSIF_SAPT_TEMP,"H2 BS RI Integrals",(char *)
        &(xBS[0][0]),sizeof(double)*aoccB_*nvirB_,next_BS,&next_BS);
      psio_->read(PSIF_SAPT_TEMP,"Theta BS Intermediate",(char *)
        &(yBS[0][0]),sizeof(double)*aoccB_*nvirB_,next_T_BS,&next_T_BS);
      h_2 += 2.0*C_DDOT(aoccB_*nvirB_,xBS[0],1,yBS[0],1);
    }
  }

  free_block(xBS);
  free_block(yBS);
  B_p_AB.done();

  return(h_2);
}

void SAPT0::h3()
{
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  SAPTDFInts A_p_AS = set_A_AS();
  Iterator AS_iter = get_iterator(mem_,&A_p_AS);

  double *tAS = init_array(noccA_*nvirB_);

  for (int i=0,off=0; i<AS_iter.num_blocks; i++) {
    read_block(&AS_iter,&A_p_AS);
    C_DGEMV('t',AS_iter.curr_size,noccA_*nvirB_,2.0,&(A_p_AS.B_p_[0][0]),
      noccA_*nvirB_,&(diagBB_[off]),1,1.0,tAS,1);
    off += AS_iter.curr_size;
  }

  A_p_AS.done();

  double **xAS = block_matrix(nthreads,aoccA_*nvirB_);

  C_DCOPY(aoccA_*nvirB_,&(tAS[foccA_*nvirB_]),1,xAS[0],1);

  free(tAS);

  SAPTDFInts A_p_AB = set_A_AB();
  SAPTDFInts B_p_BS = set_B_BS();
  Iterator ABBS_iter = get_iterator(mem_,&A_p_AB,&B_p_BS);

  for (int i=0; i<ABBS_iter.num_blocks; i++) {
    read_block(&ABBS_iter,&A_p_AB,&B_p_BS);
#pragma omp parallel
{
#pragma omp for private(rank)
    for (int j=0; j<ABBS_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',aoccA_,nvirB_,noccB_,-1.0,
        &(A_p_AB.B_p_[j][foccA_*noccB_]),noccB_,B_p_BS.B_p_[j],nvirB_,
        1.0,xAS[rank],nvirB_);
    }
}
  }

  for (int n=1; n<nthreads; n++)
    C_DAXPY(aoccA_*nvirB_,1.0,xAS[n],1,xAS[0],1);

  psio_->write_entry(PSIF_SAPT_TEMP,"H3 AS Array",(char *) &(xAS[0][0]),
    sizeof(double)*aoccA_*nvirB_);

  free_block(xAS);
}

double SAPT0::h4()
{
  double h_4 = 0.0;

  SAPTDFInts A_p_AB = set_A_AB();
  Iterator AB_iter = get_iterator(mem_,&A_p_AB);

  double **xAR = block_matrix(aoccA_,nvirA_);
  double **yAR = block_matrix(aoccA_,nvirA_);
  psio_address next_AR = PSIO_ZERO;
  psio_address next_T_AR = PSIO_ZERO;

  for (int i=0; i<AB_iter.num_blocks; i++) {
    read_block(&AB_iter,&A_p_AB);
    for (int p=0; p<AB_iter.curr_size; p++) {
      C_DGEMM('N','T',aoccA_,nvirA_,noccB_,1.0,
        &(A_p_AB.B_p_[p][foccA_*noccB_]),noccB_,&(sAB_[noccA_][0]),nmoB_,
        0.0,xAR[0],nvirA_);
      psio_->write(PSIF_SAPT_TEMP,"H4 AR RI Integrals",(char *)
        &(xAR[0][0]),sizeof(double)*aoccA_*nvirA_,next_AR,&next_AR);
      psio_->read(PSIF_SAPT_TEMP,"Theta AR Intermediate",(char *)
        &(yAR[0][0]),sizeof(double)*aoccA_*nvirA_,next_T_AR,&next_T_AR);
      h_4 += 2.0*C_DDOT(aoccA_*nvirA_,xAR[0],1,yAR[0],1);
    }
  }

  free_block(xAR);
  free_block(yAR);
  A_p_AB.done();

  return(h_4);
}

void SAPT0::q1()
{
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  SAPTDFInts A_p_AA = set_A_AA();
  SAPTDFInts B_p_BS = set_B_BS();
  Iterator AABS_iter = get_iterator(mem_,&A_p_AA,&B_p_BS);
 
  double **xAB = block_matrix(nthreads,aoccA_*noccB_);
  double **xAS = block_matrix(nthreads,aoccA_*nvirB_);
 
  for (int i=0; i<AABS_iter.num_blocks; i++) {
    read_block(&AABS_iter,&A_p_AA,&B_p_BS);
#pragma omp parallel
{
#pragma omp for private(rank)
    for (int j=0; j<AABS_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',aoccA_,noccB_,noccA_,1.0,
        &(A_p_AA.B_p_[j][foccA_*noccA_]),noccA_,&(sAB_[0][0]),nmoB_,
        0.0,xAB[rank],noccB_);
      C_DGEMM('N','N',aoccA_,nvirB_,noccB_,1.0,xAB[rank],noccB_,
        B_p_BS.B_p_[j],nvirB_,1.0,xAS[rank],nvirB_);
    }
}
  }

  for (int n=1; n<nthreads; n++)
    C_DAXPY(aoccA_*nvirB_,1.0,xAS[n],1,xAS[0],1);

  psio_->write_entry(PSIF_SAPT_TEMP,"Q1 AS Array",(char *) &(xAS[0][0]),
    sizeof(double)*aoccA_*nvirB_);

  free_block(xAB);
  free_block(xAS);
}

double SAPT0::q2()
{
  double q_2 = 0.0;

  double *sAR = init_array(noccA_*nvirA_);

  C_DGEMM('N','T',noccA_,nvirA_,noccB_,1.0,&(sAB_[0][0]),nmoB_,
    &(sAB_[noccA_][0]),nmoB_,1.0,sAR,nvirA_);

  SAPTDFInts A_p_AA = set_A_AA();
  Iterator AA_iter = get_iterator(mem_,&A_p_AA);

  double *xAR = init_array(aoccA_*nvirA_);
  double *yAR = init_array(aoccA_*nvirA_);
  psio_address next_AR = PSIO_ZERO;
  psio_address next_T_AR = PSIO_ZERO;

  for (int i=0; i<AA_iter.num_blocks; i++) {
    read_block(&AA_iter,&A_p_AA);
    for (int j=0; j<AA_iter.curr_size; j++) {
      C_DGEMM('N','N',aoccA_,nvirA_,noccA_,1.0,
        &(A_p_AA.B_p_[j][foccA_*noccA_]),noccA_,sAR,nvirA_,
        0.0,xAR,nvirA_);
      psio_->write(PSIF_SAPT_TEMP,"Q2 AR RI Integrals",(char *)
        &(xAR[0]),sizeof(double)*aoccA_*nvirA_,next_AR,&next_AR);
      psio_->read(PSIF_SAPT_TEMP,"Theta AR Intermediate",(char *)
        &(yAR[0]),sizeof(double)*aoccA_*nvirA_,next_T_AR,&next_T_AR);
      q_2 -= 2.0*C_DDOT(aoccA_*nvirA_,xAR,1,yAR,1);
    }
  }

  free(sAR);
  free(xAR);
  free(yAR);

  return(q_2);
}

void SAPT0::q3()
{
  SAPTDFInts B_p_BS = set_B_BS();
  Iterator BS_iter = get_iterator(mem_,&B_p_BS);

  double *xBS = init_array(noccB_*nvirB_);

  for (int i=0,off=0; i<BS_iter.num_blocks; i++) {
    read_block(&BS_iter,&B_p_BS);
    C_DGEMV('t',BS_iter.curr_size,noccB_*nvirB_,1.0,&(B_p_BS.B_p_[0][0]),
      noccB_*nvirB_,&(diagAA_[off]),1,1.0,xBS,1);
    off += BS_iter.curr_size;
  }

  double *xAS = init_array(aoccA_*nvirB_);

  C_DGEMM('N','N',aoccA_,nvirB_,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
    xBS,nvirB_,0.0,xAS,nvirB_);

  psio_->write_entry(PSIF_SAPT_TEMP,"Q3 AS Array",(char *) &(xAS[0]),
    sizeof(double)*aoccA_*nvirB_);

  psio_->write_entry(PSIF_SAPT_TEMP,"Q4 BS Array",(char *) 
    &(xBS[foccB_*nvirB_]),sizeof(double)*aoccB_*nvirB_);

  free(xBS);
  free(xAS);
}

void SAPT0::q5()
{
  int nthreads = 1; 
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  SAPTDFInts A_p_AR = set_A_AR();
  SAPTDFInts B_p_BB = set_B_BB();
  Iterator ARBB_iter = get_iterator(mem_,&A_p_AR,&B_p_BB);
 
  double **xAB = block_matrix(nthreads,noccA_*aoccB_);
  double **xRB = block_matrix(nthreads,nvirA_*aoccB_);
 
  for (int i=0; i<ARBB_iter.num_blocks; i++) {
    read_block(&ARBB_iter,&A_p_AR,&B_p_BB);
#pragma omp parallel
{
#pragma omp for private(rank)
    for (int j=0; j<ARBB_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','T',noccA_,aoccB_,noccB_,1.0,&(sAB_[0][0]),nmoB_,
        &(B_p_BB.B_p_[j][foccB_*noccB_]),noccB_,0.0,xAB[rank],aoccB_);
      C_DGEMM('T','N',nvirA_,aoccB_,noccA_,1.0,
        A_p_AR.B_p_[j],nvirA_,xAB[rank],aoccB_,1.0,xRB[rank],aoccB_);
    }
}
  }

  for (int n=1; n<nthreads; n++)
    C_DAXPY(nvirA_*aoccB_,1.0,xRB[n],1,xRB[0],1);

  psio_->write_entry(PSIF_SAPT_TEMP,"Q5 RB Array",(char *) &(xRB[0][0]),
    sizeof(double)*nvirA_*aoccB_);

  free_block(xAB);
  free_block(xRB);
}

double SAPT0::q6()
{
  double q_6 = 0.0;
  double *sBS = init_array(noccB_*nvirB_);

  C_DGEMM('T','N',noccB_,nvirB_,noccA_,1.0,&(sAB_[0][0]),nmoB_,
    &(sAB_[0][noccB_]),nmoB_,1.0,sBS,nvirB_);

  SAPTDFInts B_p_BB = set_B_BB();
  Iterator BB_iter = get_iterator(mem_,&B_p_BB);

  double *xBS = init_array(aoccB_*nvirB_);
  double *yBS = init_array(aoccB_*nvirB_);
  psio_address next_BS = PSIO_ZERO;
  psio_address next_T_BS = PSIO_ZERO;

  for (int i=0; i<BB_iter.num_blocks; i++) {
    read_block(&BB_iter,&B_p_BB);
    for (int j=0; j<BB_iter.curr_size; j++) {
      C_DGEMM('N','N',aoccB_,nvirB_,noccB_,1.0,
        &(B_p_BB.B_p_[j][foccB_*noccB_]),noccB_,sBS,nvirB_,
        0.0,xBS,nvirB_);
      psio_->write(PSIF_SAPT_TEMP,"Q6 BS RI Integrals",(char *)
        &(xBS[0]),sizeof(double)*aoccB_*nvirB_,next_BS,&next_BS);
      psio_->read(PSIF_SAPT_TEMP,"Theta BS Intermediate",(char *)
        &(yBS[0]),sizeof(double)*aoccB_*nvirB_,next_T_BS,&next_T_BS);
      q_6 -= 2.0*C_DDOT(aoccB_*nvirB_,xBS,1,yBS,1);
    }
  }

  free(sBS);
  free(xBS);
  free(yBS);

  return(q_6);
}

void SAPT0::q7()
{
  SAPTDFInts A_p_AR = set_A_AR();
  Iterator AR_iter = get_iterator(mem_,&A_p_AR);

  double *xAR = init_array(noccA_*nvirA_);

  for (int i=0,off=0; i<AR_iter.num_blocks; i++) {
    read_block(&AR_iter,&A_p_AR);
    C_DGEMV('t',AR_iter.curr_size,noccA_*nvirA_,1.0,&(A_p_AR.B_p_[0][0]),
      noccA_*nvirA_,&(diagBB_[off]),1,1.0,xAR,1);
    off += AR_iter.curr_size;
  }

  double *xRB = init_array(nvirA_*aoccB_);

  C_DGEMM('T','N',nvirA_,aoccB_,noccA_,1.0,xAR,nvirA_,&(sAB_[0][foccB_]),nmoB_,
    0.0,xRB,aoccB_);

  psio_->write_entry(PSIF_SAPT_TEMP,"Q7 RB Array",(char *) &(xRB[0]),
    sizeof(double)*nvirA_*aoccB_);

  psio_->write_entry(PSIF_SAPT_TEMP,"Q8 AR Array",(char *) 
    &(xAR[foccA_*nvirA_]),sizeof(double)*aoccA_*nvirA_);

  free(xAR);
  free(xRB);
}   

void SAPT0::q10()
{
  SAPTDFInts A_p_AA = set_A_AA();
  Iterator AA_iter = get_iterator(mem_,&A_p_AA);

  double *xAA = init_array(noccA_*noccA_);

  for (int i=0,off=0; i<AA_iter.num_blocks; i++) {
    read_block(&AA_iter,&A_p_AA);
    C_DGEMV('t',AA_iter.curr_size,noccA_*noccA_,1.0,&(A_p_AA.B_p_[0][0]),
      noccA_*noccA_,&(diagBB_[off]),1,1.0,xAA,1);
    off += AA_iter.curr_size;
  }

  double *xAS = init_array(aoccA_*nvirB_);

  C_DGEMM('N','N',aoccA_,nvirB_,noccA_,1.0,&(xAA[foccA_*noccA_]),noccA_,
    &(sAB_[0][noccB_]),nmoB_,0.0,xAS,nvirB_);

  psio_->write_entry(PSIF_SAPT_TEMP,"Q10 AS Array",(char *) &(xAS[0]),
    sizeof(double)*aoccA_*nvirB_);

  free(xAA);
  free(xAS);
}

void SAPT0::q11()
{
  SAPTDFInts B_p_BB = set_B_BB();
  Iterator BB_iter = get_iterator(mem_,&B_p_BB);

  double *xBB = init_array(noccB_*noccB_);

  for (int i=0,off=0; i<BB_iter.num_blocks; i++) {
    read_block(&BB_iter,&B_p_BB);
    C_DGEMV('t',BB_iter.curr_size,noccB_*noccB_,1.0,&(B_p_BB.B_p_[0][0]),
      noccB_*noccB_,&(diagAA_[off]),1,1.0,xBB,1);
    off += BB_iter.curr_size;
  }

  double *xRB = init_array(nvirA_*aoccB_);

  C_DGEMM('N','T',nvirA_,aoccB_,noccB_,1.0,&(sAB_[noccA_][0]),nmoB_,
    &(xBB[foccB_*noccB_]),noccB_,0.0,xRB,aoccB_);

  psio_->write_entry(PSIF_SAPT_TEMP,"Q11 RB Array",(char *) &(xRB[0]),
    sizeof(double)*nvirA_*aoccB_);

  free(xBB);
  free(xRB);
}

void SAPT0::q12()
{
  int nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  int rank = 0;

  zero_disk(PSIF_SAPT_TEMP,"Q12 AS RI Integrals",aoccA_*nvirB_,ndf_+3);
  zero_disk(PSIF_SAPT_TEMP,"Q12 BR RI Integrals",aoccB_*nvirA_,ndf_+3);

  SAPTDFInts A_p_AR = set_A_AR();
  Iterator AR_iter = get_iterator(mem_/2,&A_p_AR);
 
  double **xBR = block_matrix(nthreads,aoccB_*nvirA_); 
  double **X_BR_p = block_matrix(aoccB_*nvirA_,AR_iter.block_size[0]);
  psio_address next_B_BR = PSIO_ZERO;

  for (int i=0,off=0; i<AR_iter.num_blocks; i++) {
    read_block(&AR_iter,&A_p_AR);
#pragma omp parallel
{
#pragma omp for private(rank)
    for (int j=0; j<AR_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('T','N',aoccB_,nvirA_,noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
        &(A_p_AR.B_p_[j][0]),nvirA_,0.0,xBR[rank],nvirA_);
      for (int b=0; b<aoccB_; b++) {
        C_DCOPY(nvirA_,&(xBR[rank][b*nvirA_]),1,&(X_BR_p[b*nvirA_][j]),
          AR_iter.block_size[0]);
      }
    }
}

    int skip = ndf_ + 3 - AR_iter.curr_size;

    next_B_BR = psio_get_address(PSIO_ZERO,sizeof(double)*off);

    for (int br=0; br<aoccB_*nvirA_; br++) {
      psio_->write(PSIF_SAPT_TEMP,"Q12 BR RI Integrals",(char *)
        &(X_BR_p[br][0]),sizeof(double)*AR_iter.curr_size,next_B_BR,
        &next_B_BR);
      next_B_BR = psio_get_address(next_B_BR,sizeof(double)*skip);
    }

    off += AR_iter.curr_size;
  }

  free_block(xBR);
  free_block(X_BR_p);
  A_p_AR.done();

  SAPTDFInts B_p_BS = set_B_BS();
  Iterator BS_iter = get_iterator(mem_/2,&B_p_BS);

  double **xAS = block_matrix(nthreads,aoccA_*nvirB_); 
  double **X_AS_p = block_matrix(aoccA_*nvirB_,BS_iter.block_size[0]);
  psio_address next_AS = PSIO_ZERO;

  for (int i=0,off=0; i<BS_iter.num_blocks; i++) {
    read_block(&BS_iter,&B_p_BS);
#pragma omp parallel
{
#pragma omp for private(rank)
    for (int j=0; j<BS_iter.curr_size; j++) {

#ifdef _OPENMP
      rank = omp_get_thread_num();
#endif

      C_DGEMM('N','N',aoccA_,nvirB_,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
        &(B_p_BS.B_p_[j][0]),nvirB_,0.0,xAS[rank],nvirB_);
      for (int a=0; a<aoccA_; a++) {
        C_DCOPY(nvirB_,&(xAS[rank][a*nvirB_]),1,&(X_AS_p[a*nvirB_][j]),
          BS_iter.block_size[0]);
      }
    }
}
    int skip = ndf_ + 3 - BS_iter.curr_size;

    next_AS = psio_get_address(PSIO_ZERO,sizeof(double)*off);

    for (int as=0; as<aoccA_*nvirB_; as++) {
      psio_->write(PSIF_SAPT_TEMP,"Q12 AS RI Integrals",(char *)
        &(X_AS_p[as][0]),sizeof(double)*BS_iter.curr_size,next_AS,
        &next_AS);
      next_AS = psio_get_address(next_AS,sizeof(double)*skip);
    }
    
    off += BS_iter.curr_size;
  }

  free_block(xAS);
  free_block(X_AS_p);
}

double SAPT0::q13()
{
  double q_13 = 0.0;
  double **sBB = block_matrix(aoccB_,noccB_);

  C_DGEMM('T','N',aoccB_,noccB_,noccA_,1.0,&(sAB_[0][foccB_]),nmoB_,
    &(sAB_[0][0]),nmoB_,0.0,sBB[0],noccB_);

  SAPTDFInts B_p_BS = set_B_BS();
  Iterator BS_iter = get_iterator(mem_,&B_p_BS);

  double *xBS = init_array(aoccB_*nvirB_);
  double *yBS = init_array(aoccB_*nvirB_);
  psio_address next_BS = PSIO_ZERO;
  psio_address next_T_BS = PSIO_ZERO;

  for (int i=0; i<BS_iter.num_blocks; i++) {
    read_block(&BS_iter,&B_p_BS);
    for (int j=0; j<BS_iter.curr_size; j++) {
      C_DGEMM('N','N',aoccB_,nvirB_,noccB_,1.0,sBB[0],noccB_,
        &(B_p_BS.B_p_[j][0]),nvirB_,0.0,xBS,nvirB_);
      psio_->write(PSIF_SAPT_TEMP,"Q13 BS RI Integrals",(char *)
        &(xBS[0]),sizeof(double)*aoccB_*nvirB_,next_BS,&next_BS);
      psio_->read(PSIF_SAPT_TEMP,"Theta BS Intermediate",(char *)
        &(yBS[0]),sizeof(double)*aoccB_*nvirB_,next_T_BS,&next_T_BS);
      q_13 -= 2.0*C_DDOT(aoccB_*nvirB_,xBS,1,yBS,1);
    }
  }

  free(xBS);
  free(yBS);
  free_block(sBB);

  return(q_13);
}

double SAPT0::q14()
{
  double q_14 = 0.0;
  double **sAA = block_matrix(aoccA_,noccA_);

  C_DGEMM('N','T',aoccA_,noccA_,noccB_,1.0,&(sAB_[foccA_][0]),nmoB_,
    &(sAB_[0][0]),nmoB_,0.0,sAA[0],noccA_);

  SAPTDFInts A_p_AR = set_A_AR();
  Iterator AR_iter = get_iterator(mem_,&A_p_AR);

  double *xAR = init_array(aoccA_*nvirA_);
  double *yAR = init_array(aoccA_*nvirA_);
  psio_address next_AR = PSIO_ZERO;
  psio_address next_T_AR = PSIO_ZERO;

  for (int i=0; i<AR_iter.num_blocks; i++) {
    read_block(&AR_iter,&A_p_AR);
    for (int j=0; j<AR_iter.curr_size; j++) {
      C_DGEMM('N','N',aoccA_,nvirA_,noccA_,1.0,sAA[0],noccA_,
        &(A_p_AR.B_p_[j][0]),nvirA_,0.0,xAR,nvirA_);
      psio_->write(PSIF_SAPT_TEMP,"Q14 AR RI Integrals",(char *)
        &(xAR[0]),sizeof(double)*aoccA_*nvirA_,next_AR,&next_AR);
      psio_->read(PSIF_SAPT_TEMP,"Theta AR Intermediate",(char *)
        &(yAR[0]),sizeof(double)*aoccA_*nvirA_,next_T_AR,&next_T_AR);
      q_14 -= 2.0*C_DDOT(aoccA_*nvirA_,xAR,1,yAR,1);
    }
  }

  free(xAR);
  free(yAR);
  free_block(sAA);

  return(q_14);
}

void SAPT2::exch_disp20()
{
  double **yARBS = block_matrix(noccA_*nvirA_,noccB_*nvirB_);

  double **B_p_AA = get_AA_ints(1);
  double **B_p_BB = get_BB_ints(1);
  double **B_p_AB = get_AB_ints(1);
  double **C_p_AB = get_AB_ints(2);

  double **B_p_AS = get_AS_ints(1);
  double **B_p_RB = get_RB_ints(1);
  double **B_p_AR = get_AR_ints(1);
  double **B_p_BS = get_BS_ints(1);

  double **X_AA = block_matrix(noccA_,noccA_);
  double **X_BB = block_matrix(noccB_,noccB_);
  double **X_AR = block_matrix(noccA_,nvirA_);
  double **X_BS = block_matrix(noccB_,nvirB_);
  double **X_AS = block_matrix(noccA_,nvirB_);
  double **X_RB = block_matrix(nvirA_,noccB_);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','T',nvirA_*noccB_,nvirB_,
      (ndf_+3),1.0,&(B_p_RB[0][0]),(ndf_+3),
      &(B_p_AS[a*nvirB_][0]),(ndf_+3),0.0,
      &(yARBS[a*nvirA_][0]),nvirB_);
  }

  double **D_p_AR = block_matrix(noccA_*nvirA_,
    (ndf_+3));

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','N',nvirA_,(ndf_+3),noccB_,1.0,
      &(sAB_[noccA_][0]),nmoB_,
      &(B_p_AB[a*noccB_][0]),(ndf_+3),0.0,
      &(D_p_AR[a*nvirA_][0]),(ndf_+3));
  }

  C_DGEMM('N','T',noccA_*nvirA_,noccB_*
    nvirB_,(ndf_+3),2.0,&(D_p_AR[0][0]),(ndf_+3),
    &(B_p_BS[0][0]),(ndf_+3),1.0,&(yARBS[0][0]),noccB_*
    nvirB_);

  free_block(D_p_AR);

  memset(&(X_AS[0][0]),'\0',sizeof(double)*noccA_*nvirB_);

  for (int b=0; b<noccB_; b++) {
    C_DGEMM('N','T',noccA_,nvirB_,(ndf_+3),1.0,
      &(B_p_AB[b][0]),noccB_*(ndf_+3),
      &(B_p_BS[b*nvirB_][0]),(ndf_+3),1.0,&(X_AS[0][0]),
      nvirB_);
  }

  for (int a=0,ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0; b<noccB_; b++) {
        C_DAXPY(nvirB_,-sAB_[r+noccA_][b],
          &(X_AS[a][0]),1,&(yARBS[ar][b*nvirB_]),1);
  }}}

  double **D_p_RB = block_matrix(nvirA_*noccB_,
    (ndf_+3));

  C_DGEMM('N','N',nvirA_,noccB_*(ndf_+3),
    noccB_,1.0,&(sAB_[noccA_][0]),
    nmoB_,&(B_p_BB[0][0]),noccB_*(ndf_+3),0.0,
    &(D_p_RB[0][0]),noccB_*(ndf_+3));

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','T',nvirA_*noccB_,nvirB_,
      (ndf_+3),-1.0,&(D_p_RB[0][0]),(ndf_+3),
      &(B_p_AS[a*nvirB_][0]),(ndf_+3),1.0,
      &(yARBS[a*nvirA_][0]),nvirB_);
  }

  free_block(D_p_RB);

  C_DGEMV('n',noccA_*nvirB_,(ndf_+3),1.0,
    &(B_p_AS[0][0]),(ndf_+3),diagBB_,1,0.0,&(X_AS[0][0]),1);

  for (int a=0,ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0; b<noccB_; b++) {
        C_DAXPY(nvirB_,2.0*sAB_[r+noccA_][b],
          &(X_AS[a][0]),1,&(yARBS[ar][b*nvirB_]),1);
  }}}

  double **D_p_BS = block_matrix(noccB_*nvirB_,
    (ndf_+3));

  for (int b=0; b<noccB_; b++) {
    C_DGEMM('T','N',nvirB_,(ndf_+3),noccA_,1.0,
      &(sAB_[0][noccB_]),nmoB_,
      &(C_p_AB[b][0]),noccB_*(ndf_+3),0.0,
      &(D_p_BS[b*nvirB_][0]),(ndf_+3));
  }

  C_DGEMM('N','T',noccA_*nvirA_,noccB_*
    nvirB_,(ndf_+3),2.0,&(B_p_AR[0][0]),(ndf_+3),
    &(D_p_BS[0][0]),(ndf_+3),1.0,&(yARBS[0][0]),noccB_*
    nvirB_);

  free_block(D_p_BS);

  memset(&(X_RB[0][0]),'\0',sizeof(double)*nvirA_*noccB_);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','T',nvirA_,noccB_,(ndf_+3),1.0,
      &(B_p_AR[a*nvirA_][0]),(ndf_+3),
      &(C_p_AB[a*noccB_][0]),(ndf_+3),1.0,&(X_RB[0][0]),
      noccB_);
  }

  for (int a=0,ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0; b<noccB_; b++) {
        C_DAXPY(nvirB_,-X_RB[r][b],
          &(sAB_[a][noccB_]),1,
          &(yARBS[ar][b*nvirB_]),1);
  }}}

  double **D_p_AS = block_matrix(noccA_*nvirB_,
    (ndf_+3));

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('T','N',nvirB_,(ndf_+3),noccA_,1.0,
      &(sAB_[0][noccB_]),nmoB_,
      &(B_p_AA[a*noccA_][0]),(ndf_+3),0.0,
      &(D_p_AS[a*nvirB_][0]),(ndf_+3));
  }

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','T',nvirA_*noccB_,nvirB_,
      (ndf_+3),-1.0,&(B_p_RB[0][0]),(ndf_+3),
      &(D_p_AS[a*nvirB_][0]),(ndf_+3),1.0,
      &(yARBS[a*nvirA_][0]),nvirB_);
  }

  free_block(D_p_AS);

  C_DGEMV('n',nvirA_*noccB_,(ndf_+3),1.0,
    &(B_p_RB[0][0]),(ndf_+3),diagAA_,1,0.0,&(X_RB[0][0]),1);

  for (int a=0,ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0; b<noccB_; b++) {
        C_DAXPY(nvirB_,2.0*X_RB[r][b],
          &(sAB_[a][noccB_]),1,
          &(yARBS[ar][b*nvirB_]),1);
  }}}

  C_DGEMM('T','N',noccB_,noccB_,noccA_,
    1.0,&(sAB_[0][0]),nmoB_,&(sAB_[0][0]),
    nmoB_,0.0,&(X_BB[0][0]),noccB_);

  D_p_BS = block_matrix(noccB_*nvirB_,(ndf_+3));

  C_DGEMM('N','N',noccB_,nvirB_*(ndf_+3),
    noccB_,1.0,&(X_BB[0][0]),noccB_,&(B_p_BS[0][0]),
    nvirB_*(ndf_+3),0.0,&(D_p_BS[0][0]),nvirB_*
    (ndf_+3));

  C_DGEMM('N','T',noccA_*nvirA_,noccB_*
    nvirB_,(ndf_+3),-2.0,&(B_p_AR[0][0]),(ndf_+3),
    &(D_p_BS[0][0]),(ndf_+3),1.0,&(yARBS[0][0]),noccB_*
    nvirB_);

  free_block(D_p_BS);

  C_DGEMM('T','N',noccB_,nvirB_,noccA_,
    1.0,&(sAB_[0][0]),nmoB_,
    &(sAB_[0][noccB_]),nmoB_,0.0,&(X_BS[0][0]),
    nvirB_);

  D_p_BS = block_matrix(noccB_*nvirB_,(ndf_+3));

  for (int b=0; b<noccB_; b++) {
    C_DGEMM('T','N',nvirB_,(ndf_+3),noccB_,1.0,
      &(X_BS[0][0]),nvirB_,&(B_p_BB[b*noccB_][0]),
      (ndf_+3),0.0,&(D_p_BS[b*nvirB_][0]),(ndf_+3));
  }

  C_DGEMM('N','T',noccA_*nvirA_,noccB_*
    nvirB_,(ndf_+3),-2.0,&(B_p_AR[0][0]),(ndf_+3),
    &(D_p_BS[0][0]),(ndf_+3),1.0,&(yARBS[0][0]),noccB_*
    nvirB_);

  free_block(D_p_BS);

  C_DGEMV('n',noccA_*nvirA_,(ndf_+3),1.0,
    &(B_p_AR[0][0]),(ndf_+3),diagBB_,1,0.0,&(X_AR[0][0]),1);

  for (int a=0,ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      C_DAXPY(noccB_*nvirB_,4.0*X_AR[a][r],&(X_BS[0][0]),1,
        &(yARBS[ar][0]),1);
  }}

  C_DGEMM('N','T',noccA_,noccA_,noccB_,
    1.0,&(sAB_[0][0]),nmoB_,&(sAB_[0][0]),
    nmoB_,0.0,&(X_AA[0][0]),noccA_);

  D_p_AR = block_matrix(noccA_*nvirA_,(ndf_+3));

  C_DGEMM('N','N',noccA_,nvirA_*(ndf_+3),
    noccA_,1.0,&(X_AA[0][0]),noccA_,&(B_p_AR[0][0]),
    nvirA_*(ndf_+3),0.0,&(D_p_AR[0][0]),nvirA_*
    (ndf_+3));

  C_DGEMM('N','T',noccA_*nvirA_,noccB_*
    nvirB_,(ndf_+3),-2.0,&(D_p_AR[0][0]),(ndf_+3),
    &(B_p_BS[0][0]),(ndf_+3),1.0,&(yARBS[0][0]),noccB_*
    nvirB_);

  free_block(D_p_AR);

  C_DGEMM('N','T',noccA_,nvirA_,noccB_,
    1.0,&(sAB_[0][0]),nmoB_,
    &(sAB_[noccA_][0]),nmoB_,0.0,&(X_AR[0][0]),
    nvirA_);

  D_p_AR = block_matrix(noccA_*nvirA_,(ndf_+3));

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('T','N',nvirA_,(ndf_+3),noccA_,1.0,
      &(X_AR[0][0]),nvirA_,&(B_p_AA[a*noccA_][0]),
      (ndf_+3),0.0,&(D_p_AR[a*nvirA_][0]),(ndf_+3));
  }

  C_DGEMM('N','T',noccA_*nvirA_,noccB_*
    nvirB_,(ndf_+3),-2.0,&(D_p_AR[0][0]),(ndf_+3),
    &(B_p_BS[0][0]),(ndf_+3),1.0,&(yARBS[0][0]),noccB_*
    nvirB_);

  free_block(D_p_AR);

  C_DGEMV('n',noccB_*nvirB_,(ndf_+3),1.0,
    &(B_p_BS[0][0]),(ndf_+3),diagAA_,1,0.0,&(X_BS[0][0]),1);

  for (int a=0,ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      C_DAXPY(noccB_*nvirB_,4.0*X_AR[a][r],&(X_BS[0][0]),1,
        &(yARBS[ar][0]),1);
  }}

  D_p_AS = block_matrix(noccA_*nvirB_,(ndf_+3));
  D_p_RB = block_matrix(nvirA_*noccB_,(ndf_+3));

  C_DGEMM('N','N',noccA_,nvirB_*(ndf_+3),
    noccB_,1.0,&(sAB_[0][0]),nmoB_,
    &(B_p_BS[0][0]),nvirB_*(ndf_+3),0.0,&(D_p_AS[0][0]),
    nvirB_*(ndf_+3));

  for (int r=0; r<nvirA_; r++) {
    C_DGEMM('T','N',noccB_,(ndf_+3),noccA_,1.0,
      &(sAB_[0][0]),nmoB_,&(B_p_AR[r][0]),
      nvirA_*(ndf_+3),0.0,&(D_p_RB[r*noccB_][0]),
      (ndf_+3));
  }

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','T',nvirA_*noccB_,nvirB_,
      (ndf_+3),1.0,&(D_p_RB[0][0]),(ndf_+3),
      &(D_p_AS[a*nvirB_][0]),(ndf_+3),1.0,
      &(yARBS[a*nvirA_][0]),nvirB_);
  }

  C_DGEMM('N','T',nvirA_,noccB_,noccB_*
    (ndf_+3),1.0,&(D_p_RB[0][0]),noccB_*(ndf_+3),
    &(B_p_BB[0][0]),noccB_*(ndf_+3),0.0,&(X_RB[0][0]),
    noccB_);

  memset(&(X_AS[0][0]),'\0',sizeof(double)*noccA_*nvirB_);

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','T',noccA_,nvirB_,(ndf_+3),1.0,
      &(B_p_AA[a*noccA_][0]),(ndf_+3),
      &(D_p_AS[a*nvirB_][0]),(ndf_+3),1.0,&(X_AS[0][0]),
      nvirB_);
  }

  for (int a=0,ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0; b<noccB_; b++) {
        C_DAXPY(nvirB_,X_RB[r][b],
          &(sAB_[a][noccB_]),1,
          &(yARBS[ar][b*nvirB_]),1);
  }}}

  for (int a=0,ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0; b<noccB_; b++) {
        C_DAXPY(nvirB_,sAB_[r+noccA_][b],
          &(X_AS[a][0]),1,&(yARBS[ar][b*nvirB_]),1);
  }}}

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('T','N',nvirB_,(ndf_+3),noccA_,1.0,
      &(sAB_[0][noccB_]),nmoB_,
      &(B_p_AA[a*noccA_][0]),(ndf_+3),0.0,
      &(D_p_AS[a*nvirB_][0]),(ndf_+3));
  }

  C_DGEMM('N','N',nvirA_,noccB_*(ndf_+3),
    noccB_,1.0,&(sAB_[noccA_][0]),
    nmoB_,&(B_p_BB[0][0]),noccB_*(ndf_+3),0.0,
    &(D_p_RB[0][0]),noccB_*(ndf_+3));

  for (int a=0; a<noccA_; a++) {
    C_DGEMM('N','T',nvirA_*noccB_,nvirB_,
      (ndf_+3),1.0,&(D_p_RB[0][0]),(ndf_+3),
      &(D_p_AS[a*nvirB_][0]),(ndf_+3),1.0,
      &(yARBS[a*nvirA_][0]),nvirB_);
  }

  free_block(D_p_AS);
  free_block(D_p_RB);

  C_DGEMV('n',noccA_*nvirA_,(ndf_+3),1.0,
    &(B_p_AR[0][0]),(ndf_+3),diagBB_,1,0.0,&(X_AR[0][0]),1);

  C_DGEMV('n',noccB_*nvirB_,(ndf_+3),1.0,
    &(B_p_BS[0][0]),(ndf_+3),diagAA_,1,0.0,&(X_BS[0][0]),1);

  C_DGEMM('T','N',nvirA_,noccB_,noccA_,1.0,
    &(X_AR[0][0]),nvirA_,&(sAB_[0][0]),nmoB_,
    0.0,&(X_RB[0][0]),noccB_);

  C_DGEMM('N','N',noccA_,nvirB_,noccB_,1.0,
    &(sAB_[0][0]),nmoB_,&(X_BS[0][0]),nvirB_,
    0.0,&(X_AS[0][0]),nvirB_);

  for (int a=0,ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0; b<noccB_; b++) {
        C_DAXPY(nvirB_,-2.0*X_RB[r][b],
          &(sAB_[a][noccB_]),1,
          &(yARBS[ar][b*nvirB_]),1);
  }}}

  for (int a=0,ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0; b<noccB_; b++) {
        C_DAXPY(nvirB_,-2.0*sAB_[r+noccA_][b],
          &(X_AS[a][0]),1,&(yARBS[ar][b*nvirB_]),1);
  }}}

  C_DGEMV('n',noccA_*noccA_,(ndf_+3),1.0,
    &(B_p_AA[0][0]),(ndf_+3),diagBB_,1,0.0,&(X_AA[0][0]),1);

  C_DGEMV('n',noccB_*noccB_,(ndf_+3),1.0,
    &(B_p_BB[0][0]),(ndf_+3),diagAA_,1,0.0,&(X_BB[0][0]),1);

  C_DGEMM('N','N',noccA_,nvirB_,noccA_,1.0,
    &(X_AA[0][0]),noccA_,&(sAB_[0][noccB_]),
    nmoB_,0.0,&(X_AS[0][0]),nvirB_);

  C_DGEMM('N','N',nvirA_,noccB_,noccB_,1.0,
    &(sAB_[noccA_][0]),nmoB_,&(X_BB[0][0]),
    noccB_,0.0,&(X_RB[0][0]),noccB_);

  for (int a=0,ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0; b<noccB_; b++) {
        C_DAXPY(nvirB_,-2.0*X_RB[r][b],
          &(sAB_[a][noccB_]),1,
          &(yARBS[ar][b*nvirB_]),1);
  }}}

  for (int a=0,ar=0; a<noccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      for (int b=0; b<noccB_; b++) {
        C_DAXPY(nvirB_,-2.0*sAB_[r+noccA_][b],
          &(X_AS[a][0]),1,&(yARBS[ar][b*nvirB_]),1);
  }}}

  double **tARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
  
  psio_->read_entry(PSIF_SAPT_AMPS,"tARBS Amplitudes",(char *) tARBS[0],
    sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_); 

  e_exch_disp20_ = 0.0;

  for (int a=0,ar=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++,ar++) {
      e_exch_disp20_ -= 2.0*C_DDOT((long int) aoccB_*nvirB_,
        tARBS[ar],1,&(yARBS[(a+foccA_)*nvirA_+r][foccB_*nvirB_]),1);
  }}

  free_block(tARBS);
  free_block(B_p_AA);
  free_block(B_p_BB);
  free_block(B_p_AB);
  free_block(C_p_AB);
  free_block(B_p_AS);
  free_block(B_p_RB);
  free_block(B_p_AR);
  free_block(B_p_BS);
  free_block(X_AA);
  free_block(X_BB);
  free_block(X_AR);
  free_block(X_BS);
  free_block(X_AS);
  free_block(X_RB);

  if (print_) {
    fprintf(outfile,"    Exch-Disp20         = %18.12lf H\n",e_exch_disp20_);
    fflush(outfile);
  }

  psio_->write_entry(PSIF_SAPT_AMPS,"Exch-Disp V_ARBS",(char *) yARBS[0],
    sizeof(double)*noccA_*nvirA_*noccB_*nvirB_);

  free_block(yARBS);
}

}}
