#include "sapt2p3.h"

using namespace boost;

namespace psi { namespace sapt {

void SAPT2p3::ind30()
{
  double **tAR = block_matrix(noccA_,nvirA_);

  psio_->read_entry(PSIF_SAPT_AMPS,"Ind30 uAR Amplitudes", (char *) tAR[0],
    sizeof(double)*noccA_*nvirA_);

  double indA_B = 2.0*C_DDOT(noccA_*nvirA_,tAR[0],1,wBAR_[0],1);

  free_block(tAR);

  double **tBS = block_matrix(noccB_,nvirB_);

  psio_->read_entry(PSIF_SAPT_AMPS,"Ind30 uBS Amplitudes", (char *) tBS[0],
    sizeof(double)*noccB_*nvirB_);

  double indB_A = 2.0*C_DDOT(noccB_*nvirB_,tBS[0],1,wABS_[0],1);

  free_block(tBS);

  e_ind30_ = indA_B + indB_A;

  if (debug_) {
    fprintf(outfile,"\n    Ind30_1             = %18.12lf H\n",indA_B);
    fprintf(outfile,"    Ind30_2             = %18.12lf H\n",indB_A);
  }
  if (print_) {
    fprintf(outfile,"    Ind30               = %18.12lf H\n",e_ind30_);
    fflush(outfile);
  }
}

}}

