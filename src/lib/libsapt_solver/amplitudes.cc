#include "sapt2.h"

namespace psi { namespace sapt {

void SAPT2::amplitudes()
{
}

void SAPT2::tOVOV(int intfileA, const char *ARlabel, int foccA, int noccA,
  int nvirA, double *evalsA, int intfileB, const char *BSlabel, int foccB,
  int noccB, int nvirB, double *evalsB, int ampout, const char *amplabel)
{
  int aoccA = noccA - foccA;
  int aoccB = noccB - foccB;

  double **B_p_AR = get_DF_ints(intfileA,ARlabel,foccA,noccA,0,nvirA);
  double **B_p_BS = get_DF_ints(intfileB,BSlabel,foccB,noccB,0,nvirB);
  double **tARBS = block_matrix(aoccA*nvirA,aoccB*nvirB);

  C_DGEMM('N','T',aoccA*nvirA,aoccB*nvirB,ndf_,1.0,B_p_AR[0],ndf_+3,
    B_p_BS[0],ndf_+3,0.0,tARBS[0],aoccB*nvirB);

  for (int a=0, ar=0; a<aoccA; a++){
    for (int r=0; r<nvirA; r++, ar++){
      for (int b=0, bs=0; b<aoccB; b++){
        for (int s=0; s<nvirB; s++, bs++){
          tARBS[ar][bs] /= evalsA[a+foccA]+evalsB[b+foccB]-evalsA[r+noccA]-
            evalsB[s+noccB];
  }}}}

  psio_->write_entry(ampout,amplabel,(char *) tARBS[0],sizeof(double)*aoccA*
    nvirA*aoccB*nvirB); 
}

}}
