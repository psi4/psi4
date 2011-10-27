#include "sapt0.h"
#include "sapt2.h"

namespace psi { namespace sapt {
/*
 *
 * This dispersion evaluation will work with an arbitrary amount of memory,
 * however, it is not optimal if everything would fit in core.
 *
 * This exists only for debugging purposes to test the accuracy of Laplace
 * denominators.
 *
 */
void SAPT0::disp20()
{
  long int avail_mem = mem_ - (long int) ndf_*ndf_;

  SAPTDFInts B_p_AR = set_act_C_AR();
  SAPTDFInts B_p_BS = set_act_C_BS();
  Iterator B_ARBS_iter = get_iterator(avail_mem/3,&B_p_AR,&B_p_BS);

  SAPTDFInts C_p_AR = set_act_C_AR();
  SAPTDFInts C_p_BS = set_act_C_BS();
  Iterator C_ARBS_iter = get_iterator(avail_mem/3,&C_p_AR,&C_p_BS);

  double *xPQ = init_array((long int) B_ARBS_iter.block_size[0]*
    C_ARBS_iter.block_size[0]);
  double *yPQ = init_array((long int) B_ARBS_iter.block_size[0]*
    C_ARBS_iter.block_size[0]);

  double **T_p_AR = block_matrix(C_ARBS_iter.block_size[0],aoccA_*nvirA_);
  double **T_p_BS = block_matrix(C_ARBS_iter.block_size[0],aoccB_*nvirB_);

  e_disp20_ = 0.0;
  for (int j=0; j<B_ARBS_iter.num_blocks; j++) {
    read_block(&B_ARBS_iter,&B_p_AR,&B_p_BS);

    for (int k=0; k<C_ARBS_iter.num_blocks; k++) {
      read_block(&C_ARBS_iter,&C_p_AR,&C_p_BS);

      for (int i=0; i<nvec_; i++) {

        C_DCOPY(C_ARBS_iter.block_size[k]*aoccA_*nvirA_,C_p_AR.B_p_[0],1,
          T_p_AR[0],1);
        C_DCOPY(C_ARBS_iter.block_size[k]*aoccB_*nvirB_,C_p_BS.B_p_[0],1,
          T_p_BS[0],1);

#pragma omp parallel
{
#pragma omp for
        for (int ar=0; ar<aoccA_*nvirA_; ar++) {
          double scale = dAR_[i][ar];
          C_DSCAL(C_ARBS_iter.curr_size,scale,&(T_p_AR[0][ar]),aoccA_*nvirA_);
        }
      
#pragma omp for 
        for (int bs=0; bs<aoccB_*nvirB_; bs++) {
          double scale = dBS_[i][bs];
          C_DSCAL(C_ARBS_iter.curr_size,scale,&(T_p_BS[0][bs]),aoccB_*nvirB_);
        }
}

        C_DGEMM('N','T',B_ARBS_iter.curr_size,C_ARBS_iter.curr_size,
          aoccA_*nvirA_,2.0,B_p_AR.B_p_[0],aoccA_*nvirA_,T_p_AR[0],
          aoccA_*nvirA_,0.0,xPQ,C_ARBS_iter.curr_size);

        C_DGEMM('N','T',B_ARBS_iter.curr_size,C_ARBS_iter.curr_size,
          aoccB_*nvirB_,2.0,B_p_BS.B_p_[0],aoccB_*nvirB_,T_p_BS[0],
          aoccB_*nvirB_,0.0,yPQ,C_ARBS_iter.curr_size);

        e_disp20_ -= C_DDOT(B_ARBS_iter.curr_size*C_ARBS_iter.curr_size,
          xPQ,1,yPQ,1);
      }
    }
    C_p_AR.rewind();
    C_p_BS.rewind();
    C_ARBS_iter.rewind();
  }

  B_p_AR.done();
  C_p_AR.done();
  B_p_BS.done();
  C_p_BS.done();

  free(xPQ);
  free(yPQ);

  free_block(T_p_AR);
  free_block(T_p_BS);

  if (print_) {
    fprintf(outfile,"    Disp20              = %18.12lf H\n",e_disp20_);
    fflush(outfile);
  }
}

void SAPT2::disp20()
{
  double **B_p_AR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    foccA_,noccA_,0,nvirA_);
  double **B_p_BS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
    foccB_,noccB_,0,nvirB_);
  double **vARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);

  C_DGEMM('N','T',aoccA_*nvirA_,aoccB_*nvirB_,ndf_,1.0,B_p_AR[0],ndf_+3,
    B_p_BS[0],ndf_+3,0.0,vARBS[0],aoccB_*nvirB_);

  free_block(B_p_AR);
  free_block(B_p_BS);

  double **tARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);

  psio_->read_entry(PSIF_SAPT_AMPS,"tARBS Amplitudes",(char *) tARBS[0],
    sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_); 

  e_disp20_ = 4.0*C_DDOT((long int) aoccA_*nvirA_*aoccB_*nvirB_,vARBS[0],1,
    tARBS[0],1);

  if (print_) {
    fprintf(outfile,"    Disp20              = %18.12lf H\n",e_disp20_);
    fflush(outfile);
  }

  free_block(tARBS);
  free_block(vARBS);

  if (nat_orbs_) {
    double **C_p_AR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AR NO RI Integrals",
      foccA_,noccA_,0,no_nvirA_);
    double **C_p_BS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BS NO RI Integrals",
      foccB_,noccB_,0,no_nvirB_);
    vARBS = block_matrix(aoccA_*no_nvirA_,aoccB_*no_nvirB_);

    C_DGEMM('N','T',aoccA_*no_nvirA_,aoccB_*no_nvirB_,ndf_,1.0,C_p_AR[0],
      ndf_+3,C_p_BS[0],ndf_+3,0.0,vARBS[0],aoccB_*no_nvirB_);
  
    free_block(C_p_AR);
    free_block(C_p_BS);

    e_no_disp20_ = 0.0;

    for (int a=0, ar=0; a<aoccA_; a++){
      for (int r=0; r<no_nvirA_; r++, ar++){
        for (int b=0, bs=0; b<aoccB_; b++){
          for (int s=0; s<no_nvirB_; s++, bs++){
            double tval = vARBS[ar][bs];
            e_no_disp20_ += 4.0*tval*tval / (no_evalsA_[a+foccA_]
              +no_evalsB_[b+foccB_]-no_evalsA_[r+noccA_]-no_evalsB_[s+noccB_]);
    }}}}

    free_block(vARBS);

    if (print_) {
      fprintf(outfile,"    Disp20 (NO)         = %18.12lf H\n",e_no_disp20_);
      fflush(outfile);
    }
  }
}
/*
 * This version is probably worse...unless there isn't much memory available
 *
void SAPT0::disp20()
{
  boost::shared_ptr<Vector> evals_aoccA(new Vector(aoccA_));
  boost::shared_ptr<Vector> evals_virA(new Vector(nvirA_));
  boost::shared_ptr<Vector> evals_aoccB(new Vector(aoccB_));
  boost::shared_ptr<Vector> evals_virB(new Vector(nvirB_));

  for (int a=0; a<aoccA_; a++)
    evals_aoccA->set(0,a,evalsA_[a+foccA_]);
  for (int r=0; r<nvirA_; r++)
    evals_virA->set(0,r,evalsA_[r+noccA_]);
  for (int b=0; b<aoccB_; b++)
    evals_aoccB->set(0,b,evalsB_[b+foccB_]);
  for (int s=0; s<nvirB_; s++)
    evals_virB->set(0,s,evalsB_[s+noccB_]);

  denom_ = boost::shared_ptr<SAPTLaplaceDenominator>(new 
    SAPTLaplaceDenominator(evals_aoccA,evals_virA,evals_aoccB,evals_virB,
    options_.get_double("DENOMINATOR_DELTA"),debug_));

  SharedMatrix tauAR = denom_->denominatorA();
  SharedMatrix tauBS = denom_->denominatorB();

  dAR_ = tauAR->pointer();
  dBS_ = tauBS->pointer();

  nvec_ = denom_->nvector();

  psio_->open(PSIF_SAPT_TEMP,PSIO_OPEN_NEW);

  long int avail_mem = mem_ - (long int) ndf_*ndf_;

  psio_address next_xPQ = PSIO_ZERO;
  psio_address next_yPQ = PSIO_ZERO;

  SAPTDFInts B_p_AR = set_act_C_AR();
  SAPTDFInts C_p_AR = set_act_C_AR();
  Iterator B_AR_iter = get_iterator(avail_mem/2,&B_p_AR);
  Iterator C_AR_iter = get_iterator(avail_mem/2,&C_p_AR);

  double **xPQ = block_matrix(ndf_,ndf_);

  for (int i=0; i<nvec_; i++) {
    for (int j=0, Boff=0; j<B_AR_iter.num_blocks; j++) {
      read_block(&B_AR_iter,&B_p_AR);

      for (int a=0,ar=0; a<aoccA_; a++) {
        for (int r=0; r<nvirA_; r++,ar++) {
          C_DSCAL(B_AR_iter.curr_size,dAR_[i][ar],&(B_p_AR.B_p_[0][ar]),
            aoccA_*nvirA_);
        }
      }

      for (int k=0, Coff=0; k<C_AR_iter.num_blocks; k++) {
        read_block(&C_AR_iter,&C_p_AR);

        C_DGEMM('N','T',B_AR_iter.curr_size,C_AR_iter.curr_size,aoccA_*nvirA_,
          2.0,B_p_AR.B_p_[0],aoccA_*nvirA_,C_p_AR.B_p_[0],aoccA_*nvirA_,
          0.0,&(xPQ[Boff][Coff]),ndf_);

        Coff += C_AR_iter.curr_size;
      }
      C_p_AR.rewind();
      C_AR_iter.rewind();
      Boff += B_AR_iter.curr_size;
    }
    B_p_AR.rewind();
    B_AR_iter.rewind();

    psio_->write(PSIF_SAPT_TEMP,"X PQ Matrices",(char *) &(xPQ[0][0]),
      sizeof(double)*ndf_*ndf_,next_xPQ,&next_xPQ);
  }

  B_p_AR.done();
  C_p_AR.done();

  SAPTDFInts B_p_BS = set_act_C_BS();
  SAPTDFInts C_p_BS = set_act_C_BS();
  Iterator B_BS_iter = get_iterator(avail_mem/2,&B_p_BS);
  Iterator C_BS_iter = get_iterator(avail_mem/2,&C_p_BS);

  for (int i=0; i<nvec_; i++) {
     for (int j=0, Boff=0; j<B_BS_iter.num_blocks; j++) {
      read_block(&B_BS_iter,&B_p_BS);

      for (int b=0,bs=0; b<aoccB_; b++) {
        for (int s=0; s<nvirB_; s++,bs++) {
          C_DSCAL(B_BS_iter.curr_size,dBS_[i][bs],&(B_p_BS.B_p_[0][bs]),
            aoccB_*nvirB_);
        }
      }

      for (int k=0, Coff=0; k<C_BS_iter.num_blocks; k++) {
        read_block(&C_BS_iter,&C_p_BS);

        C_DGEMM('N','T',B_BS_iter.curr_size,C_BS_iter.curr_size,aoccB_*nvirB_,
          2.0,B_p_BS.B_p_[0],aoccB_*nvirB_,C_p_BS.B_p_[0],aoccB_*nvirB_,
          0.0,&(xPQ[Boff][Coff]),ndf_);

        Coff += C_BS_iter.curr_size;
      }
      C_p_BS.rewind();
      C_BS_iter.rewind();
      Boff += B_BS_iter.curr_size;
    }
    B_p_BS.rewind();
    B_BS_iter.rewind();

    psio_->write(PSIF_SAPT_TEMP,"Y PQ Matrices",(char *) &(xPQ[0][0]),
      sizeof(double)*ndf_*ndf_,next_yPQ,&next_yPQ);
  }

  B_p_BS.done();
  C_p_BS.done();

  next_xPQ = PSIO_ZERO;
  next_yPQ = PSIO_ZERO;

  double **yPQ = block_matrix(ndf_,ndf_);

  if (debug_)
    fprintf(outfile,"\n");

  e_disp20_ = 0.0;
  for (int i=0; i<nvec_; i++) {
    psio_->read(PSIF_SAPT_TEMP,"X PQ Matrices",(char *) &(xPQ[0][0]),
      sizeof(double)*ndf_*ndf_,next_xPQ,&next_xPQ);
    psio_->read(PSIF_SAPT_TEMP,"Y PQ Matrices",(char *) &(yPQ[0][0]),
      sizeof(double)*ndf_*ndf_,next_yPQ,&next_yPQ);
    double tval = C_DDOT(ndf_*ndf_,xPQ[0],1,yPQ[0],1);
    e_disp20_ -= tval;
    if (debug_)
      fprintf(outfile,"    Disp %2d             = %18.12lf H\n",i+1,-tval);
  }

  if (debug_)
    fprintf(outfile,"\n");

  free_block(xPQ);
  free_block(yPQ);

  if (print_) {
    fprintf(outfile,"    Disp20              = %18.12lf H\n",e_disp20_);
    fflush(outfile);
  }

  psio_->close(PSIF_SAPT_TEMP,0);
}
*/
}}
