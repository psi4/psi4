#include "sapt0.h"

namespace psi { namespace sapt {

void SAPT0::exch_disp20()
{
  psio_->open(PSIF_SAPT_TEMP,PSIO_OPEN_NEW);

  v1();
  h1();
  h2();
  h3();
  h4();
  q1();
  q2();
  q3();
  q5();
  q6();
  q7();
  q10();
  q11();
  q12();
  q13();
  q14();

	double v_1 = 0.0;
  double h_1 = 0.0, h_2 = 0.0, h_3 = 0.0, h_4 = 0.0;
  double q_1 = 0.0, q_2 = 0.0, q_3 = 0.0, q_4 = 0.0; 
  double q_5 = 0.0, q_6 = 0.0, q_7 = 0.0, q_8 = 0.0;
  double q_9 = 0.0, q_10 = 0.0, q_11 = 0.0, q_12 = 0.0; 
  double q_13 = 0.0, q_14 = 0.0;

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

  C_DGEMM('N','T',aoccA_,nvirA_,noccB_,1.0,&(sAB_[foccA_][0]),nmo_,
    &(sAB_[noccA_][0]),nmo_,0.0,sAR[0],nvirA_);

  C_DGEMM('T','N',aoccB_,nvirB_,noccA_,1.0,&(sAB_[0][foccB_]),nmo_,
    &(sAB_[0][noccB_]),nmo_,0.0,sBS[0],nvirB_);

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

  SAPTDFInts T_p_AR = set_act_C_AR();
  SAPTDFInts T_p_BS = set_act_C_BS();
  Iterator T_ARBS_iter = get_iterator(mem_/3,&T_p_AR,&T_p_BS);

  SAPTDFInts A_p_AS = set_act_A_AS();
  SAPTDFInts B_p_RB = set_act_B_RB();
  Iterator ASRB_iter = get_iterator(mem_/3,&A_p_AS,&B_p_RB);

  SAPTDFInts A_p_AR = set_act_A_AR();
  SAPTDFInts B_p_BS = set_act_B_BS();
  Iterator ARBS_iter = get_iterator(mem_/3,&A_p_AR,&B_p_BS);

  SAPTDFInts H2_p_BS = set_H2_BS();
  SAPTDFInts H2_p_AS = set_H2_AS();
  SAPTDFInts H4_p_AR = set_H4_AR();
  SAPTDFInts H4_p_RB = set_H4_RB();
  SAPTDFInts Q2_p_AR = set_Q2_AR();
  SAPTDFInts Q6_p_BS = set_Q6_BS();
  SAPTDFInts Q12_p_AS = set_Q12_AS();
  SAPTDFInts Q12_p_RB = set_Q12_RB();
  SAPTDFInts Q13_p_BS = set_Q13_BS();
  SAPTDFInts Q14_p_AR = set_Q14_AR();

  H2_p_BS.B_p_ = B_p_BS.B_p_;
  H2_p_AS.B_p_ = A_p_AS.B_p_;
  H4_p_AR.B_p_ = A_p_AR.B_p_;
  H4_p_RB.B_p_ = B_p_RB.B_p_;
  Q2_p_AR.B_p_ = A_p_AR.B_p_;
  Q6_p_BS.B_p_ = B_p_BS.B_p_;
  Q12_p_AS.B_p_ = A_p_AS.B_p_;
  Q12_p_RB.B_p_ = B_p_RB.B_p_;
  Q13_p_BS.B_p_ = B_p_BS.B_p_;
  Q14_p_AR.B_p_ = A_p_AR.B_p_;

  Iterator H2_ASRB_iter = get_iterator(mem_/3,&H2_p_AS,&B_p_RB,false);
  Iterator H2_ARBS_iter = get_iterator(mem_/3,&A_p_AR,&H2_p_BS,false);
  Iterator H4_ASRB_iter = get_iterator(mem_/3,&A_p_AS,&H4_p_RB,false);
  Iterator H4_ARBS_iter = get_iterator(mem_/3,&H4_p_AR,&B_p_BS,false);
  Iterator Q2_ARBS_iter = get_iterator(mem_/3,&Q2_p_AR,&B_p_BS,false);
  Iterator Q6_ARBS_iter = get_iterator(mem_/3,&A_p_AR,&Q6_p_BS,false);
  Iterator Q9_ASRB_iter = get_iterator(mem_/3,&H2_p_AS,&H4_p_RB,false);
  Iterator Q12_ASRB_iter = get_iterator(mem_/3,&Q12_p_AS,&Q12_p_RB,false);
  Iterator Q13_ARBS_iter = get_iterator(mem_/3,&A_p_AR,&Q13_p_BS,false);
  Iterator Q14_ARBS_iter = get_iterator(mem_/3,&Q14_p_AR,&B_p_BS,false);

  double *xAB = init_array(aoccA_*aoccB_);
  double *yAB = init_array(aoccA_*aoccB_);

  double *xPQ = init_array(T_ARBS_iter.block_size[0]*
    H2_ARBS_iter.block_size[0]);
  double *yPQ = init_array(T_ARBS_iter.block_size[0]*
    H2_ARBS_iter.block_size[0]);

  for (int j=0; j<T_ARBS_iter.num_blocks; j++) {
    read_block(&T_ARBS_iter,&T_p_AR,&T_p_BS);

    for (int i=0; i<nvec_; i++) {

      for (int a=0,ar=0; a<aoccA_; a++) { 
        for (int r=0; r<nvirA_; r++,ar++) {
          double scale = dAR_[i][ar];
          if (i) scale /= dAR_[i-1][ar];
          C_DSCAL(T_ARBS_iter.curr_size,scale,&(T_p_AR.B_p_[0][ar]),
            aoccA_*nvirA_);
         }
      }

      for (int b=0,bs=0; b<aoccB_; b++) { 
        for (int s=0; s<nvirB_; s++,bs++) {
          double scale = dBS_[i][bs];
          if (i) scale /= dBS_[i-1][bs];
          C_DSCAL(T_ARBS_iter.curr_size,scale,&(T_p_BS.B_p_[0][bs]),
            aoccB_*nvirB_);
        }
      }

      for (int k=0; k<ASRB_iter.num_blocks; k++) {
        read_block(&ASRB_iter,&A_p_AS,&B_p_RB);

        for (int p=0; p<T_ARBS_iter.curr_size; p++) {
          for (int q=0; q<ASRB_iter.curr_size; q++) {
            C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR.B_p_[p],nvirA_,
              B_p_RB.B_p_[q],aoccB_,0.0,xAB,aoccB_);
            C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,A_p_AS.B_p_[q],nvirB_,
              T_p_BS.B_p_[p],nvirB_,0.0,yAB,aoccB_);
            v_1 += C_DDOT(aoccA_*aoccB_,xAB,1,yAB,1);
          }
        }
      }
      A_p_AS.rewind();
      B_p_RB.rewind();
      ASRB_iter.rewind();

      for (int p=0; p<T_ARBS_iter.curr_size; p++) {
        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR.B_p_[p],nvirA_,
          H1_RB[0],aoccB_,0.0,xAB,aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,sAS[0],nvirB_,
          T_p_BS.B_p_[p],nvirB_,0.0,yAB,aoccB_);
        h_1 += C_DDOT(aoccA_*aoccB_,xAB,1,yAB,1);
      }

      for (int k=0; k<H2_ARBS_iter.num_blocks; k++) {
        read_block(&H2_ARBS_iter,&A_p_AR,&H2_p_BS);

        C_DGEMM('N','T',T_ARBS_iter.curr_size,H2_ARBS_iter.curr_size,
          aoccA_*nvirA_,1.0,T_p_AR.B_p_[0],aoccA_*nvirA_,A_p_AR.B_p_[0],
          aoccA_*nvirA_,0.0,xPQ,H2_ARBS_iter.curr_size);

        C_DGEMM('N','T',T_ARBS_iter.curr_size,H2_ARBS_iter.curr_size,
          aoccB_*nvirB_,1.0,T_p_BS.B_p_[0],aoccB_*nvirB_,H2_p_BS.B_p_[0],
          aoccB_*nvirB_,0.0,yPQ,H2_ARBS_iter.curr_size);

        h_2 += 2.0*C_DDOT(T_ARBS_iter.curr_size*H2_ARBS_iter.curr_size,
          xPQ,1,yPQ,1);
      }
      A_p_AR.rewind();
      H2_p_BS.rewind();
      H2_ARBS_iter.rewind();

      for (int k=0; k<H2_ASRB_iter.num_blocks; k++) {
        read_block(&H2_ASRB_iter,&H2_p_AS,&B_p_RB);

        for (int p=0; p<T_ARBS_iter.curr_size; p++) {
          for (int q=0; q<H2_ASRB_iter.curr_size; q++) {
            C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR.B_p_[p],nvirA_,
              B_p_RB.B_p_[q],aoccB_,0.0,xAB,aoccB_);
            C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,H2_p_AS.B_p_[q],nvirB_,
              T_p_BS.B_p_[p],nvirB_,0.0,yAB,aoccB_);
            h_2 -= C_DDOT(aoccA_*aoccB_,xAB,1,yAB,1);
          }
        }
      }
      H2_p_AS.rewind();
      B_p_RB.rewind();
      H2_ASRB_iter.rewind();

      for (int p=0; p<T_ARBS_iter.curr_size; p++) {
        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR.B_p_[p],nvirA_,
          sRB[0],aoccB_,0.0,xAB,aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,H3_AS[0],nvirB_,
          T_p_BS.B_p_[p],nvirB_,0.0,yAB,aoccB_);
        h_3 += C_DDOT(aoccA_*aoccB_,xAB,1,yAB,1);
      }

      for (int k=0; k<H4_ARBS_iter.num_blocks; k++) {
        read_block(&H4_ARBS_iter,&H4_p_AR,&B_p_BS);

        C_DGEMM('N','T',T_ARBS_iter.curr_size,H4_ARBS_iter.curr_size,
          aoccA_*nvirA_,1.0,T_p_AR.B_p_[0],aoccA_*nvirA_,H4_p_AR.B_p_[0],
          aoccA_*nvirA_,0.0,xPQ,H4_ARBS_iter.curr_size);

        C_DGEMM('N','T',T_ARBS_iter.curr_size,H4_ARBS_iter.curr_size,
          aoccB_*nvirB_,1.0,T_p_BS.B_p_[0],aoccB_*nvirB_,B_p_BS.B_p_[0],
          aoccB_*nvirB_,0.0,yPQ,H4_ARBS_iter.curr_size);

        h_4 += 2.0*C_DDOT(T_ARBS_iter.curr_size*H4_ARBS_iter.curr_size,
          xPQ,1,yPQ,1);
      }
      H4_p_AR.rewind();
      B_p_BS.rewind();
      H4_ARBS_iter.rewind();

      for (int k=0; k<H4_ASRB_iter.num_blocks; k++) {
        read_block(&H4_ASRB_iter,&A_p_AS,&H4_p_RB);

        for (int p=0; p<T_ARBS_iter.curr_size; p++) {
          for (int q=0; q<H4_ASRB_iter.curr_size; q++) {
            C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR.B_p_[p],nvirA_,
              H4_p_RB.B_p_[q],aoccB_,0.0,xAB,aoccB_);
            C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,A_p_AS.B_p_[q],nvirB_,
              T_p_BS.B_p_[p],nvirB_,0.0,yAB,aoccB_);
            h_4 -= C_DDOT(aoccA_*aoccB_,xAB,1,yAB,1);
          }
        }
      }
      A_p_AS.rewind();
      H4_p_RB.rewind();
      H4_ASRB_iter.rewind();

      for (int p=0; p<T_ARBS_iter.curr_size; p++) {
        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR.B_p_[p],nvirA_,
          sRB[0],aoccB_,0.0,xAB,aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,Q1_AS[0],nvirB_,
          T_p_BS.B_p_[p],nvirB_,0.0,yAB,aoccB_);
        q_1 += C_DDOT(aoccA_*aoccB_,xAB,1,yAB,1);
      } 

      for (int k=0; k<Q2_ARBS_iter.num_blocks; k++) {
        read_block(&Q2_ARBS_iter,&Q2_p_AR,&B_p_BS);

        C_DGEMM('N','T',T_ARBS_iter.curr_size,Q2_ARBS_iter.curr_size,
          aoccA_*nvirA_,1.0,T_p_AR.B_p_[0],aoccA_*nvirA_,Q2_p_AR.B_p_[0],
          aoccA_*nvirA_,0.0,xPQ,Q2_ARBS_iter.curr_size);

        C_DGEMM('N','T',T_ARBS_iter.curr_size,Q2_ARBS_iter.curr_size,
          aoccB_*nvirB_,1.0,T_p_BS.B_p_[0],aoccB_*nvirB_,B_p_BS.B_p_[0],
          aoccB_*nvirB_,0.0,yPQ,Q2_ARBS_iter.curr_size);

        q_2 -= 2.0*C_DDOT(T_ARBS_iter.curr_size*Q2_ARBS_iter.curr_size,
          xPQ,1,yPQ,1);
      }
      Q2_p_AR.rewind();
      B_p_BS.rewind();
      Q2_ARBS_iter.rewind();

      for (int p=0; p<T_ARBS_iter.curr_size; p++) {
        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR.B_p_[p],nvirA_,
          sRB[0],aoccB_,0.0,xAB,aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,Q3_AS[0],nvirB_,
          T_p_BS.B_p_[p],nvirB_,0.0,yAB,aoccB_);
        q_3 -= 2.0*C_DDOT(aoccA_*aoccB_,xAB,1,yAB,1);
      }

      for (int p=0; p<T_ARBS_iter.curr_size; p++) {
        q_4 += 4.0*C_DDOT(aoccA_*nvirA_,T_p_AR.B_p_[p],1,sAR[0],1)*
          C_DDOT(aoccB_*nvirB_,T_p_BS.B_p_[p],1,Q4_BS[0],1);
      }

      for (int p=0; p<T_ARBS_iter.curr_size; p++) {
        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR.B_p_[p],nvirA_,
          Q5_RB[0],aoccB_,0.0,xAB,aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,sAS[0],nvirB_,
          T_p_BS.B_p_[p],nvirB_,0.0,yAB,aoccB_);
        q_5 += C_DDOT(aoccA_*aoccB_,xAB,1,yAB,1);
      }

      for (int k=0; k<Q6_ARBS_iter.num_blocks; k++) {
        read_block(&Q6_ARBS_iter,&A_p_AR,&Q6_p_BS);

        C_DGEMM('N','T',T_ARBS_iter.curr_size,Q6_ARBS_iter.curr_size,
          aoccA_*nvirA_,1.0,T_p_AR.B_p_[0],aoccA_*nvirA_,A_p_AR.B_p_[0],
          aoccA_*nvirA_,0.0,xPQ,Q6_ARBS_iter.curr_size);

        C_DGEMM('N','T',T_ARBS_iter.curr_size,Q6_ARBS_iter.curr_size,
          aoccB_*nvirB_,1.0,T_p_BS.B_p_[0],aoccB_*nvirB_,Q6_p_BS.B_p_[0],
          aoccB_*nvirB_,0.0,yPQ,Q6_ARBS_iter.curr_size);

        q_6 -= 2.0*C_DDOT(T_ARBS_iter.curr_size*Q6_ARBS_iter.curr_size,
          xPQ,1,yPQ,1);
      }
      A_p_AR.rewind();
      Q6_p_BS.rewind();
      Q6_ARBS_iter.rewind();

      for (int p=0; p<T_ARBS_iter.curr_size; p++) { 
        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR.B_p_[p],nvirA_,
          Q7_RB[0],aoccB_,0.0,xAB,aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,sAS[0],nvirB_,
          T_p_BS.B_p_[p],nvirB_,0.0,yAB,aoccB_);
        q_7 -= 2.0*C_DDOT(aoccA_*aoccB_,xAB,1,yAB,1);
      }

      for (int p=0; p<T_ARBS_iter.curr_size; p++) {
        q_8 += 4.0*C_DDOT(aoccA_*nvirA_,T_p_AR.B_p_[p],1,Q8_AR[0],1)*
          C_DDOT(aoccB_*nvirB_,T_p_BS.B_p_[p],1,sBS[0],1);
      }

      for (int k=0; k<Q9_ASRB_iter.num_blocks; k++) {
        read_block(&Q9_ASRB_iter,&H2_p_AS,&H4_p_RB);

        for (int p=0; p<T_ARBS_iter.curr_size; p++) {
          for (int q=0; q<Q9_ASRB_iter.curr_size; q++) {
            C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR.B_p_[p],nvirA_,
              H4_p_RB.B_p_[q],aoccB_,0.0,xAB,aoccB_);
            C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,H2_p_AS.B_p_[q],nvirB_,
              T_p_BS.B_p_[p],nvirB_,0.0,yAB,aoccB_);
            q_9 += C_DDOT(aoccA_*aoccB_,xAB,1,yAB,1);
          }
        }
      }
      H2_p_AS.rewind();
      H4_p_RB.rewind();
      Q9_ASRB_iter.rewind();

      for (int p=0; p<T_ARBS_iter.curr_size; p++) {
        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR.B_p_[p],nvirA_,
          sRB[0],aoccB_,0.0,xAB,aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,Q10_AS[0],nvirB_,
          T_p_BS.B_p_[p],nvirB_,0.0,yAB,aoccB_);
        q_10 -= 2.0*C_DDOT(aoccA_*aoccB_,xAB,1,yAB,1);
      }

      for (int p=0; p<T_ARBS_iter.curr_size; p++) {
        C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR.B_p_[p],nvirA_,
          Q11_RB[0],aoccB_,0.0,xAB,aoccB_);
        C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,sAS[0],nvirB_,
          T_p_BS.B_p_[p],nvirB_,0.0,yAB,aoccB_);
        q_11 -= 2.0*C_DDOT(aoccA_*aoccB_,xAB,1,yAB,1);
      }

      for (int k=0; k<Q12_ASRB_iter.num_blocks; k++) {
        read_block(&Q12_ASRB_iter,&Q12_p_AS,&Q12_p_RB);

        for (int p=0; p<T_ARBS_iter.curr_size; p++) {
          for (int q=0; q<Q12_ASRB_iter.curr_size; q++) {
            C_DGEMM('N','N',aoccA_,aoccB_,nvirA_,1.0,T_p_AR.B_p_[p],nvirA_,
              Q12_p_RB.B_p_[q],aoccB_,0.0,xAB,aoccB_);
            C_DGEMM('N','T',aoccA_,aoccB_,nvirB_,1.0,Q12_p_AS.B_p_[q],nvirB_,
              T_p_BS.B_p_[p],nvirB_,0.0,yAB,aoccB_);
            q_12 += C_DDOT(aoccA_*aoccB_,xAB,1,yAB,1);
          }
        }
      }
      Q12_p_AS.rewind();
      Q12_p_RB.rewind();
      Q12_ASRB_iter.rewind();

      for (int k=0; k<Q13_ARBS_iter.num_blocks; k++) {
        read_block(&Q13_ARBS_iter,&A_p_AR,&Q13_p_BS);

        C_DGEMM('N','T',T_ARBS_iter.curr_size,Q13_ARBS_iter.curr_size,
          aoccA_*nvirA_,1.0,T_p_AR.B_p_[0],aoccA_*nvirA_,A_p_AR.B_p_[0],
          aoccA_*nvirA_,0.0,xPQ,Q13_ARBS_iter.curr_size);

        C_DGEMM('N','T',T_ARBS_iter.curr_size,Q13_ARBS_iter.curr_size,
          aoccB_*nvirB_,1.0,T_p_BS.B_p_[0],aoccB_*nvirB_,Q13_p_BS.B_p_[0],
          aoccB_*nvirB_,0.0,yPQ,Q13_ARBS_iter.curr_size);

        q_13 -= 2.0*C_DDOT(T_ARBS_iter.curr_size*Q13_ARBS_iter.curr_size,
          xPQ,1,yPQ,1);
      }
      A_p_AR.rewind();
      Q13_p_BS.rewind();
      Q13_ARBS_iter.rewind();

      for (int k=0; k<Q14_ARBS_iter.num_blocks; k++) {
        read_block(&Q14_ARBS_iter,&Q14_p_AR,&B_p_BS);

        C_DGEMM('N','T',T_ARBS_iter.curr_size,Q14_ARBS_iter.curr_size,
          aoccA_*nvirA_,1.0,T_p_AR.B_p_[0],aoccA_*nvirA_,Q14_p_AR.B_p_[0],
          aoccA_*nvirA_,0.0,xPQ,Q14_ARBS_iter.curr_size);

        C_DGEMM('N','T',T_ARBS_iter.curr_size,Q14_ARBS_iter.curr_size,
          aoccB_*nvirB_,1.0,T_p_BS.B_p_[0],aoccB_*nvirB_,B_p_BS.B_p_[0],
          aoccB_*nvirB_,0.0,yPQ,Q14_ARBS_iter.curr_size);

        q_14 -= 2.0*C_DDOT(T_ARBS_iter.curr_size*Q14_ARBS_iter.curr_size,
          xPQ,1,yPQ,1);
      }
      Q14_p_AR.rewind();
      B_p_BS.rewind();
      Q14_ARBS_iter.rewind();

    }
  }

  if (debug_) {
    fprintf(outfile,"\n    V1                  = %18.12lf H\n",v_1);
    fprintf(outfile,"    H1                  = %18.12lf H\n",h_1);
    fprintf(outfile,"    H2                  = %18.12lf H\n",h_2);
    fprintf(outfile,"    H3                  = %18.12lf H\n",h_3);
    fprintf(outfile,"    H4                  = %18.12lf H\n",h_4);
    fprintf(outfile,"    Q1                  = %18.12lf H\n",q_1);
    fprintf(outfile,"    Q2                  = %18.12lf H\n",q_2);
    fprintf(outfile,"    Q3                  = %18.12lf H\n",q_3);
    fprintf(outfile,"    Q4                  = %18.12lf H\n",q_4);
    fprintf(outfile,"    Q5                  = %18.12lf H\n",q_5);
    fprintf(outfile,"    Q6                  = %18.12lf H\n",q_6);
    fprintf(outfile,"    Q7                  = %18.12lf H\n",q_7);
    fprintf(outfile,"    Q8                  = %18.12lf H\n",q_8);
    fprintf(outfile,"    Q9                  = %18.12lf H\n",q_9);
    fprintf(outfile,"    Q10                 = %18.12lf H\n",q_10);
    fprintf(outfile,"    Q11                 = %18.12lf H\n",q_11);
    fprintf(outfile,"    Q12                 = %18.12lf H\n",q_12);
    fprintf(outfile,"    Q13                 = %18.12lf H\n",q_13);
    fprintf(outfile,"    Q14                 = %18.12lf H\n\n",q_14);
  }

  H2_p_BS.B_p_ = NULL;
  H2_p_AS.B_p_ = NULL;
  H4_p_AR.B_p_ = NULL;
  H4_p_RB.B_p_ = NULL;
  Q2_p_AR.B_p_ = NULL;
  Q6_p_BS.B_p_ = NULL;
  Q12_p_AS.B_p_ = NULL;
  Q12_p_RB.B_p_ = NULL;
  Q13_p_BS.B_p_ = NULL;
  Q14_p_AR.B_p_ = NULL;

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

  free(xAB);
  free(yAB);
  free(xPQ);
  free(yPQ);

  e_exch_disp20_ = 2.0*(v_1+h_1+h_2+h_3+h_4+q_1+q_2+q_3+q_4+q_5+q_6+q_7
    +q_8+q_9+q_10+q_11+q_12+q_13+q_14);

  if (print_) {
    fprintf(outfile,"    Exch-Disp20         = %18.12lf H\n",e_exch_disp20_);
    fflush(outfile);
  }

  psio_->close(PSIF_SAPT_TEMP,0);
}

void SAPT0::v1()
{
  SAPTDFInts A_p_AS = set_A_AS();
  Iterator AS_iter = get_iterator(mem_,&A_p_AS);
  double **xAS = block_matrix(aoccA_,nvirB_);

  psio_address next_A_AS = PSIO_ZERO;

  for (int i=0; i<AS_iter.num_blocks; i++) {
    read_block(&AS_iter,&A_p_AS);
    for (int p=0; p<AS_iter.curr_size; p++) {
      C_DCOPY(aoccA_*nvirB_,&(A_p_AS.B_p_[p][foccA_*nvirB_]),1,xAS[0],1);
      psio_->write(PSIF_SAPT_TEMP,"Active AS RI Integrals",(char *) 
        &(xAS[0][0]),sizeof(double)*aoccA_*nvirB_,next_A_AS,&next_A_AS);
    }
  }

  free_block(xAS);
  A_p_AS.done();

  SAPTDFInts B_p_RB = set_B_RB();
  Iterator RB_iter = get_iterator(mem_,&B_p_RB);
  double **xRB = block_matrix(nvirA_,aoccB_);

  psio_address next_B_RB = PSIO_ZERO;

  for (int i=0; i<RB_iter.num_blocks; i++) {
    read_block(&RB_iter,&B_p_RB);
    for (int p=0; p<RB_iter.curr_size; p++) {
      for (int r=0; r<nvirA_; r++) {
        C_DCOPY(aoccB_,&(B_p_RB.B_p_[p][r*noccB_+foccB_]),1,xRB[r],1);
      }
      psio_->write(PSIF_SAPT_TEMP,"Active RB RI Integrals",(char *) 
        &(xRB[0][0]),sizeof(double)*nvirA_*aoccB_,next_B_RB,&next_B_RB);
    }
  }

  free_block(xRB);
  B_p_RB.done();
}

void SAPT0::h1()
{
  SAPTDFInts B_p_RB = set_act_B_RB();
  Iterator RB_iter = get_iterator(mem_,&B_p_RB);

  double *xRB = init_array(nvirA_*aoccB_);

  for (int i=0,off=0; i<RB_iter.num_blocks; i++) {
    read_block(&RB_iter,&B_p_RB);
    C_DGEMV('t',RB_iter.curr_size,nvirA_*aoccB_,2.0,&(B_p_RB.B_p_[0][0]),
      nvirA_*aoccB_,&(diagAA_[off]),1,1.0,xRB,1);
    off += RB_iter.curr_size;
  }

  B_p_RB.done();

  SAPTDFInts A_p_AR = set_A_AR();
  SAPTDFInts B_p_AB = set_B_AB();
  Iterator ARAB_iter = get_iterator(mem_,&A_p_AR,&B_p_AB);

  for (int i=0; i<ARAB_iter.num_blocks; i++) {
    read_block(&ARAB_iter,&A_p_AR,&B_p_AB);
    for (int j=0; j<ARAB_iter.curr_size; j++) {
      C_DGEMM('T','N',nvirA_,aoccB_,noccA_,-1.0,A_p_AR.B_p_[j],nvirA_,
        &(B_p_AB.B_p_[j][foccB_]),noccB_,1.0,xRB,aoccB_);
    }
  }

  A_p_AR.done();
  B_p_AB.done();

  psio_->write_entry(PSIF_SAPT_TEMP,"H1 RB Array",(char *) &(xRB[0]),
    sizeof(double)*nvirA_*aoccB_);

  free(xRB);
}

void SAPT0::h2()
{
  SAPTDFInts B_p_AB = set_B_AB();
  Iterator AB_iter = get_iterator(mem_,&B_p_AB);

  double **xBS = block_matrix(aoccB_,nvirB_);
  psio_address next_BS = PSIO_ZERO;

  for (int i=0; i<AB_iter.num_blocks; i++) {
    read_block(&AB_iter,&B_p_AB);
    for (int p=0; p<AB_iter.curr_size; p++) {
      C_DGEMM('T','N',aoccB_,nvirB_,noccA_,1.0,&(B_p_AB.B_p_[p][foccB_]),
        noccB_,&(sAB_[0][noccB_]),nmo_,0.0,xBS[0],nvirB_);
      psio_->write(PSIF_SAPT_TEMP,"H2 BS RI Integrals",(char *)
        &(xBS[0][0]),sizeof(double)*aoccB_*nvirB_,next_BS,&next_BS);
    }
  }

  free_block(xBS);
  B_p_AB.done();

  SAPTDFInts A_p_AA = set_A_AA();
  Iterator AA_iter = get_iterator(mem_,&A_p_AA);

  double **xAS = block_matrix(aoccA_,nvirB_);
  psio_address next_AS = PSIO_ZERO;

  for (int i=0; i<AA_iter.num_blocks; i++) {
    read_block(&AA_iter,&A_p_AA);
    for (int p=0; p<AA_iter.curr_size; p++) {
      C_DGEMM('N','N',aoccA_,nvirB_,noccA_,1.0,
        &(A_p_AA.B_p_[p][foccA_*noccA_]),noccA_,&(sAB_[0][noccB_]),nmo_,
        0.0,xAS[0],nvirB_);
      psio_->write(PSIF_SAPT_TEMP,"H2 AS RI Integrals",(char *)
        &(xAS[0][0]),sizeof(double)*aoccA_*nvirB_,next_AS,&next_AS);
    }
  }

  free_block(xAS);
  A_p_AA.done();
}

void SAPT0::h3()
{
  SAPTDFInts A_p_AS = set_act_A_AS();
  Iterator AS_iter = get_iterator(mem_,&A_p_AS);

  double *xAS = init_array(aoccA_*nvirB_);

  for (int i=0,off=0; i<AS_iter.num_blocks; i++) {
    read_block(&AS_iter,&A_p_AS);
    C_DGEMV('t',AS_iter.curr_size,aoccA_*nvirB_,2.0,&(A_p_AS.B_p_[0][0]),
      aoccA_*nvirB_,&(diagBB_[off]),1,1.0,xAS,1);
    off += AS_iter.curr_size;
  }

  A_p_AS.done();

  SAPTDFInts A_p_AB = set_A_AB();
  SAPTDFInts B_p_BS = set_B_BS();
  Iterator ABBS_iter = get_iterator(mem_,&A_p_AB,&B_p_BS);

  for (int i=0; i<ABBS_iter.num_blocks; i++) {
    read_block(&ABBS_iter,&A_p_AB,&B_p_BS);
    for (int j=0; j<ABBS_iter.curr_size; j++) {
      C_DGEMM('N','N',aoccA_,nvirB_,noccB_,-1.0,
        &(A_p_AB.B_p_[j][foccA_*noccB_]),noccB_,B_p_BS.B_p_[j],nvirB_,
        1.0,xAS,nvirB_);
    }
  }

  psio_->write_entry(PSIF_SAPT_TEMP,"H3 AS Array",(char *) &(xAS[0]),
    sizeof(double)*aoccA_*nvirB_);

  free(xAS);
}

void SAPT0::h4()
{
  SAPTDFInts A_p_AB = set_A_AB();
  Iterator AB_iter = get_iterator(mem_,&A_p_AB);

  double **xAR = block_matrix(aoccA_,nvirA_);
  psio_address next_AR = PSIO_ZERO;

  for (int i=0; i<AB_iter.num_blocks; i++) {
    read_block(&AB_iter,&A_p_AB);
    for (int p=0; p<AB_iter.curr_size; p++) {
      C_DGEMM('N','T',aoccA_,nvirA_,noccB_,1.0,
        &(A_p_AB.B_p_[p][foccA_*noccB_]),noccB_,&(sAB_[noccA_][0]),nmo_,
        0.0,xAR[0],nvirA_);
      psio_->write(PSIF_SAPT_TEMP,"H4 AR RI Integrals",(char *)
        &(xAR[0][0]),sizeof(double)*aoccA_*nvirA_,next_AR,&next_AR);
    }
  }

  free_block(xAR);
  A_p_AB.done();

  SAPTDFInts B_p_BB = set_B_BB();
  Iterator BB_iter = get_iterator(mem_,&B_p_BB);

  double **xRB = block_matrix(nvirA_,aoccB_);
  psio_address next_RB = PSIO_ZERO;

  for (int i=0; i<BB_iter.num_blocks; i++) {
    read_block(&BB_iter,&B_p_BB);
    for (int p=0; p<BB_iter.curr_size; p++) {
      C_DGEMM('N','N',nvirA_,aoccB_,noccB_,1.0,&(sAB_[noccA_][0]),nmo_,
        &(B_p_BB.B_p_[p][foccB_]),noccB_,0.0,xRB[0],aoccB_);
      psio_->write(PSIF_SAPT_TEMP,"H4 RB RI Integrals",(char *)
        &(xRB[0][0]),sizeof(double)*nvirA_*aoccB_,next_RB,&next_RB);
    }
  }

  free_block(xRB);
  B_p_BB.done();
}

void SAPT0::q1()
{
  SAPTDFInts A_p_AA = set_A_AA();
  SAPTDFInts B_p_BS = set_B_BS();
  Iterator AABS_iter = get_iterator(mem_,&A_p_AA,&B_p_BS);
 
  double *xAB = init_array(aoccA_*noccB_);
  double *xAS = init_array(aoccA_*nvirB_);
 
  for (int i=0; i<AABS_iter.num_blocks; i++) {
    read_block(&AABS_iter,&A_p_AA,&B_p_BS);
    for (int j=0; j<AABS_iter.curr_size; j++) {
      C_DGEMM('N','N',aoccA_,noccB_,noccA_,1.0,
        &(A_p_AA.B_p_[j][foccA_*noccA_]),noccA_,&(sAB_[0][0]),nmo_,
        0.0,xAB,noccB_);
      C_DGEMM('N','N',aoccA_,nvirB_,noccB_,1.0,xAB,noccB_,
        B_p_BS.B_p_[j],nvirB_,1.0,xAS,nvirB_);
    }
  }

  psio_->write_entry(PSIF_SAPT_TEMP,"Q1 AS Array",(char *) &(xAS[0]),
    sizeof(double)*aoccA_*nvirB_);

  free(xAB);
  free(xAS);
}

void SAPT0::q2()
{
  double *sAR = init_array(noccA_*nvirA_);

  C_DGEMM('N','T',noccA_,nvirA_,noccB_,1.0,&(sAB_[0][0]),nmo_,
    &(sAB_[noccA_][0]),nmo_,1.0,sAR,nvirA_);

  SAPTDFInts A_p_AA = set_A_AA();
  Iterator AA_iter = get_iterator(mem_,&A_p_AA);

  double *xAR = init_array(aoccA_*nvirA_);
  psio_address next_AR = PSIO_ZERO;

  for (int i=0; i<AA_iter.num_blocks; i++) {
    read_block(&AA_iter,&A_p_AA);
    for (int j=0; j<AA_iter.curr_size; j++) {
      C_DGEMM('N','N',aoccA_,nvirA_,noccA_,1.0,
        &(A_p_AA.B_p_[j][foccA_*noccA_]),noccA_,sAR,nvirA_,
        0.0,xAR,nvirA_);
      psio_->write(PSIF_SAPT_TEMP,"Q2 AR RI Integrals",(char *)
        &(xAR[0]),sizeof(double)*aoccA_*nvirA_,next_AR,&next_AR);
    }
  }

  free(sAR);
  free(xAR);
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

  C_DGEMM('N','N',aoccA_,nvirB_,noccB_,1.0,&(sAB_[foccA_][0]),nmo_,
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
  SAPTDFInts A_p_AR = set_A_AR();
  SAPTDFInts B_p_BB = set_B_BB();
  Iterator ARBB_iter = get_iterator(mem_,&A_p_AR,&B_p_BB);
 
  double *xAB = init_array(noccA_*aoccB_);
  double *xRB = init_array(nvirA_*aoccB_);
 
  for (int i=0; i<ARBB_iter.num_blocks; i++) {
    read_block(&ARBB_iter,&A_p_AR,&B_p_BB);
    for (int j=0; j<ARBB_iter.curr_size; j++) {
      C_DGEMM('N','T',noccA_,aoccB_,noccB_,1.0,&(sAB_[0][0]),nmo_,
        &(B_p_BB.B_p_[j][foccB_*noccB_]),noccB_,0.0,xAB,aoccB_);
      C_DGEMM('T','N',nvirA_,aoccB_,noccA_,1.0,
        A_p_AR.B_p_[j],nvirA_,xAB,aoccB_,1.0,xRB,aoccB_);
    }
  }

  psio_->write_entry(PSIF_SAPT_TEMP,"Q5 RB Array",(char *) &(xRB[0]),
    sizeof(double)*nvirA_*aoccB_);

  free(xAB);
  free(xRB);
}

void SAPT0::q6()
{
  double *sBS = init_array(noccB_*nvirB_);

  C_DGEMM('T','N',noccB_,nvirB_,noccA_,1.0,&(sAB_[0][0]),nmo_,
    &(sAB_[0][noccB_]),nmo_,1.0,sBS,nvirB_);

  SAPTDFInts B_p_BB = set_B_BB();
  Iterator BB_iter = get_iterator(mem_,&B_p_BB);

  double *xBS = init_array(aoccB_*nvirB_);
  psio_address next_BS = PSIO_ZERO;

  for (int i=0; i<BB_iter.num_blocks; i++) {
    read_block(&BB_iter,&B_p_BB);
    for (int j=0; j<BB_iter.curr_size; j++) {
      C_DGEMM('N','N',aoccB_,nvirB_,noccB_,1.0,
        &(B_p_BB.B_p_[j][foccB_*noccB_]),noccB_,sBS,nvirB_,
        0.0,xBS,nvirB_);
      psio_->write(PSIF_SAPT_TEMP,"Q6 BS RI Integrals",(char *)
        &(xBS[0]),sizeof(double)*aoccB_*nvirB_,next_BS,&next_BS);
    }
  }

  free(sBS);
  free(xBS);
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

  C_DGEMM('T','N',nvirA_,aoccB_,noccA_,1.0,xAR,nvirA_,&(sAB_[0][foccB_]),nmo_,
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
    &(sAB_[0][noccB_]),nmo_,0.0,xAS,nvirB_);

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

  C_DGEMM('N','T',nvirA_,aoccB_,noccB_,1.0,&(sAB_[noccA_][0]),nmo_,
    &(xBB[foccB_*noccB_]),noccB_,0.0,xRB,aoccB_);

  psio_->write_entry(PSIF_SAPT_TEMP,"Q11 RB Array",(char *) &(xRB[0]),
    sizeof(double)*nvirA_*aoccB_);

  free(xBB);
  free(xRB);
}

void SAPT0::q12()
{
  SAPTDFInts A_p_AR = set_A_AR();
  Iterator AR_iter = get_iterator(mem_,&A_p_AR);
  
  double *xRB = init_array(nvirA_*aoccB_);
  psio_address next_RB = PSIO_ZERO;

  for (int i=0; i<AR_iter.num_blocks; i++) {
    read_block(&AR_iter,&A_p_AR);
    for (int j=0; j<AR_iter.curr_size; j++) {
      C_DGEMM('T','N',nvirA_,aoccB_,noccA_,1.0,
        &(A_p_AR.B_p_[j][0]),nvirA_,&(sAB_[0][foccB_]),nmo_,0.0,xRB,aoccB_);
      psio_->write(PSIF_SAPT_TEMP,"Q12 RB RI Integrals",(char *)
        &(xRB[0]),sizeof(double)*nvirA_*aoccB_,next_RB,&next_RB);
    }
  }

  free(xRB);
  A_p_AR.done();

  SAPTDFInts B_p_BS = set_B_BS();
  Iterator BS_iter = get_iterator(mem_,&B_p_BS);

  double *xAS = init_array(aoccA_*nvirB_);
  psio_address next_AS = PSIO_ZERO;

  for (int i=0; i<BS_iter.num_blocks; i++) {
    read_block(&BS_iter,&B_p_BS);
    for (int j=0; j<BS_iter.curr_size; j++) {
      C_DGEMM('N','N',aoccA_,nvirB_,noccB_,1.0,&(sAB_[foccA_][0]),nmo_,
        &(B_p_BS.B_p_[j][0]),nvirB_,0.0,xAS,nvirB_);
      psio_->write(PSIF_SAPT_TEMP,"Q12 AS RI Integrals",(char *)
        &(xAS[0]),sizeof(double)*aoccA_*nvirB_,next_AS,&next_AS);
    }
  }

  free(xAS);
}

void SAPT0::q13()
{
  double **sBB = block_matrix(aoccB_,noccB_);

  C_DGEMM('T','N',aoccB_,noccB_,noccA_,1.0,&(sAB_[0][foccB_]),nmo_,
    &(sAB_[0][0]),nmo_,0.0,sBB[0],noccB_);

  SAPTDFInts B_p_BS = set_B_BS();
  Iterator BS_iter = get_iterator(mem_,&B_p_BS);

  double *xBS = init_array(aoccB_*nvirB_);
  psio_address next_BS = PSIO_ZERO;

  for (int i=0; i<BS_iter.num_blocks; i++) {
    read_block(&BS_iter,&B_p_BS);
    for (int j=0; j<BS_iter.curr_size; j++) {
      C_DGEMM('N','N',aoccB_,nvirB_,noccB_,1.0,sBB[0],noccB_,
        &(B_p_BS.B_p_[j][0]),nvirB_,0.0,xBS,nvirB_);
      psio_->write(PSIF_SAPT_TEMP,"Q13 BS RI Integrals",(char *)
        &(xBS[0]),sizeof(double)*aoccB_*nvirB_,next_BS,&next_BS);
    }
  }

  free(xBS);
  free_block(sBB);
}

void SAPT0::q14()
{
  double **sAA = block_matrix(aoccA_,noccA_);

  C_DGEMM('N','T',aoccA_,noccA_,noccB_,1.0,&(sAB_[foccA_][0]),nmo_,
    &(sAB_[0][0]),nmo_,0.0,sAA[0],noccA_);

  SAPTDFInts A_p_AR = set_A_AR();
  Iterator AR_iter = get_iterator(mem_,&A_p_AR);

  double *xAR = init_array(aoccA_*nvirA_);
  psio_address next_AR = PSIO_ZERO;

  for (int i=0; i<AR_iter.num_blocks; i++) {
    read_block(&AR_iter,&A_p_AR);
    for (int j=0; j<AR_iter.curr_size; j++) {
      C_DGEMM('N','N',aoccA_,nvirA_,noccA_,1.0,sAA[0],noccA_,
        &(A_p_AR.B_p_[j][0]),nvirA_,0.0,xAR,nvirA_);
      psio_->write(PSIF_SAPT_TEMP,"Q14 AR RI Integrals",(char *)
        &(xAR[0]),sizeof(double)*aoccA_*nvirA_,next_AR,&next_AR);
    }
  }

  free(xAR);
  free_block(sAA);
}

}}
