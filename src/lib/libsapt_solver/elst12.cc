#include "sapt2.h"

namespace psi { namespace sapt {

void SAPT2::elst12()
{
  double e_elst120 = elst120(wBAA_,wBRR_,CHFA_,PSIF_SAPT_AMPS,
    "pAA Density Matrix","pRR Density Matrix","Y2 AR Amplitudes",
    foccA_,noccA_,nvirA_);
  double e_elst102 = elst120(wABB_,wASS_,CHFB_,PSIF_SAPT_AMPS,
    "pBB Density Matrix","pSS Density Matrix","Y2 BS Amplitudes",
    foccB_,noccB_,nvirB_);

  e_elst12_ = e_elst120 + e_elst102;

  if (debug_) {
    fprintf(outfile,"    Elst120,r           = %18.12lf H\n",e_elst120);
    fprintf(outfile,"    Elst102,r           = %18.12lf H\n",e_elst102);
  }

  if (print_) {
    fprintf(outfile,"    Elst12,r            = %18.12lf H\n",e_elst12_);
    fflush(outfile);
  }
}

double SAPT2::elst120(double **wBAA, double **wBRR, double **CHFA, int ampfile, 
  const char *pAAlabel, const char *pRRlabel, const char *Ylabel, int foccA, 
  int noccA, int nvirA)
{
  int aoccA = noccA - foccA;
  
  double **pAA = block_matrix(aoccA,aoccA);
  psio_->read_entry(ampfile,pAAlabel,(char *) pAA[0],
    sizeof(double)*aoccA*aoccA);

  double **pRR = block_matrix(nvirA,nvirA);
  psio_->read_entry(ampfile,pRRlabel,(char *) pRR[0],
    sizeof(double)*nvirA*nvirA);

  double **yAR = block_matrix(aoccA,nvirA); 
  psio_->read_entry(ampfile,Ylabel,(char *) yAR[0],
    sizeof(double)*aoccA*nvirA);

  double elst = 0.0;

  for (int a=0; a<aoccA; a++) {
    elst -= 2.0*C_DDOT(aoccA,pAA[a],1,&(wBAA[a+foccA][foccA]),1);
  }

  elst += 2.0*C_DDOT(nvirA*nvirA,pRR[0],1,wBRR[0],1);
  elst += 4.0*C_DDOT(aoccA*nvirA,yAR[0],1,CHFA[foccA],1);

  free_block(pAA);
  free_block(pRR);
  free_block(yAR);

  return(elst);
}

}}
