#include "sapt2.h"
#include "sapt2p.h"
#include "sapt2p3.h"

namespace psi { namespace sapt {

void SAPT2::amplitudes()
{
  tOVOV(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",foccA_,noccA_,nvirA_,evalsA_,
    PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",foccA_,noccA_,nvirA_,evalsA_,
    PSIF_SAPT_AMPS,"tARAR Amplitudes");
  tOVOV(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",foccB_,noccB_,nvirB_,evalsB_,
    PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",foccB_,noccB_,nvirB_,evalsB_,
    PSIF_SAPT_AMPS,"tBSBS Amplitudes");

  tOVOV(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",foccA_,noccA_,nvirA_,evalsA_,
    PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",foccB_,noccB_,nvirB_,evalsB_,
    PSIF_SAPT_AMPS,"tARBS Amplitudes");

  pOOpVV(PSIF_SAPT_AMPS,"tARAR Amplitudes","tARAR Amplitudes",aoccA_,nvirA_,
    PSIF_SAPT_AMPS,"pAA Density Matrix","pRR Density Matrix");
  pOOpVV(PSIF_SAPT_AMPS,"tBSBS Amplitudes","tBSBS Amplitudes",aoccB_,nvirB_,
    PSIF_SAPT_AMPS,"pBB Density Matrix","pSS Density Matrix");

  if (nat_orbs_) {
    natural_orbitalify(PSIF_SAPT_AMPS,"pRR Density Matrix",evalsA_,foccA_,
      noccA_,nvirA_,'A');
    natural_orbitalify(PSIF_SAPT_AMPS,"pSS Density Matrix",evalsB_,foccB_,
      noccB_,nvirB_,'B');
    natural_orbitalify_df_ints();
    tOVOV(PSIF_SAPT_AA_DF_INTS,"AR NO RI Integrals",foccA_,noccA_,no_nvirA_,
      no_evalsA_,PSIF_SAPT_AA_DF_INTS,"AR NO RI Integrals",foccA_,noccA_,
      no_nvirA_,no_evalsA_,PSIF_SAPT_AMPS,"tARAR NO Amplitudes");
    tOVOV(PSIF_SAPT_BB_DF_INTS,"BS NO RI Integrals",foccB_,noccB_,no_nvirB_,
      no_evalsB_,PSIF_SAPT_BB_DF_INTS,"BS NO RI Integrals",foccB_,noccB_,
      no_nvirB_,no_evalsB_,PSIF_SAPT_AMPS,"tBSBS NO Amplitudes");
    if (print_) fprintf(outfile,"\n");
  }

  theta(PSIF_SAPT_AMPS,"tARAR Amplitudes",'N',true,aoccA_,nvirA_,aoccA_,nvirA_,
    "AR RI Integrals",PSIF_SAPT_AMPS,"Theta AR Intermediates");
  theta(PSIF_SAPT_AMPS,"tBSBS Amplitudes",'N',true,aoccB_,nvirB_,aoccB_,nvirB_,
    "BS RI Integrals",PSIF_SAPT_AMPS,"Theta BS Intermediates");

  theta(PSIF_SAPT_AMPS,"tARBS Amplitudes",'N',false,aoccA_,nvirA_,
    aoccB_,nvirB_,"BS RI Integrals",PSIF_SAPT_AMPS,"T AR Intermediates");
  theta(PSIF_SAPT_AMPS,"tARBS Amplitudes",'T',false,aoccA_,nvirA_,
    aoccB_,nvirB_,"AR RI Integrals",PSIF_SAPT_AMPS,"T BS Intermediates");

  Y2(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
    "RR RI Integrals",PSIF_SAPT_AMPS,"pAA Density Matrix","pRR Density Matrix",
    "Theta AR Intermediates",foccA_,noccA_,nvirA_,evalsA_,PSIF_SAPT_AMPS,
    "Y2 AR Amplitudes","T2 AR Amplitudes");
  Y2(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
    "SS RI Integrals",PSIF_SAPT_AMPS,"pBB Density Matrix","pSS Density Matrix",
    "Theta BS Intermediates",foccB_,noccB_,nvirB_,evalsB_,PSIF_SAPT_AMPS,
    "Y2 BS Amplitudes","T2 BS Amplitudes");

  if (nat_orbs_t2_) {
    t2OVOV(PSIF_SAPT_AMPS,"tARAR Amplitudes","tARAR NO Amplitudes",
      "Theta AR Intermediates",PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
      "AR RI Integrals","RR RI Integrals","RR NO RI Integrals",
      foccA_,noccA_,nvirA_,no_nvirA_,evalsA_,no_CA_,PSIF_SAPT_AMPS,
      "t2ARAR Amplitudes");
    t2OVOV(PSIF_SAPT_AMPS,"tBSBS Amplitudes","tBSBS NO Amplitudes",
      "Theta BS Intermediates",PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
      "BS RI Integrals","SS RI Integrals","SS NO RI Integrals",
      foccB_,noccB_,nvirB_,no_nvirB_,evalsB_,no_CB_,PSIF_SAPT_AMPS,
      "t2BSBS Amplitudes");
  }
  else {
    t2OVOV(PSIF_SAPT_AMPS,"tARAR Amplitudes","Theta AR Intermediates",
      PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
      "RR RI Integrals",foccA_,noccA_,nvirA_,evalsA_,PSIF_SAPT_AMPS,
      "t2ARAR Amplitudes");
    t2OVOV(PSIF_SAPT_AMPS,"tBSBS Amplitudes","Theta BS Intermediates",
      PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
      "SS RI Integrals",foccB_,noccB_,nvirB_,evalsB_,PSIF_SAPT_AMPS,
      "t2BSBS Amplitudes");
  }

  theta(PSIF_SAPT_AMPS,"t2ARAR Amplitudes",'N',true,aoccA_,nvirA_,aoccA_,
    nvirA_,"AR RI Integrals",PSIF_SAPT_AMPS,"Theta 2 AR Intermediates");
  theta(PSIF_SAPT_AMPS,"t2BSBS Amplitudes",'N',true,aoccB_,nvirB_,aoccB_,
    nvirB_,"BS RI Integrals",PSIF_SAPT_AMPS,"Theta 2 BS Intermediates");
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

  free_block(B_p_AR);
  free_block(B_p_BS);
  free_block(tARBS);
}

void SAPT2::pOOpVV(int ampfile, const char *amplabel, const char *thetalabel,
  int aoccA, int nvirA, int ampout, const char *OOlabel, const char *VVlabel)
{
  double **tARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  psio_->read_entry(ampfile,amplabel,(char *) tARAR[0],sizeof(double)*aoccA*
    nvirA*aoccA*nvirA);

  double **thetaARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  psio_->read_entry(ampfile,thetalabel,(char *) thetaARAR[0],
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);
  antisym(thetaARAR,aoccA,nvirA);

  double **pOO = block_matrix(aoccA,aoccA);
  double **pVV = block_matrix(nvirA,nvirA);

  C_DGEMM('N','T',aoccA,aoccA,nvirA*aoccA*nvirA,1.0,tARAR[0],nvirA*aoccA*nvirA,
    thetaARAR[0],nvirA*aoccA*nvirA,0.0,pOO[0],aoccA);

  C_DGEMM('T','N',nvirA,nvirA,aoccA*nvirA*aoccA,1.0,tARAR[0],nvirA,
    thetaARAR[0],nvirA,0.0,pVV[0],nvirA);

  psio_->write_entry(ampout,OOlabel,(char *) pOO[0],
    sizeof(double)*aoccA*aoccA);
  psio_->write_entry(ampout,VVlabel,(char *) pVV[0],
    sizeof(double)*nvirA*nvirA);

  free_block(tARAR);
  free_block(thetaARAR);
  free_block(pOO);
  free_block(pVV);
}

void SAPT2::theta(int ampfile, const char *amplabel, const char trans,
  bool antisymmetrized, int aoccA, int nvirA, int aoccB, int nvirB,
  const char *OVlabel, int ampout, const char *thetalabel)
{
  double **tARBS = block_matrix(aoccA*nvirA,aoccB*nvirB);
  psio_->read_entry(ampfile,amplabel,(char *) tARBS[0],sizeof(double)*aoccA*
    nvirA*aoccB*nvirB);

  if (antisymmetrized)
    antisym(tARBS,aoccA,nvirA);

  double **B_p_OV;
  if (strcmp(OVlabel, "AR RI Integrals"))
    B_p_OV = get_AR_ints(1,foccA_);
  else if (strcmp(OVlabel, "BS RI Integrals"))
    B_p_OV = get_BS_ints(1,foccB_);
  else
    throw PsiException("Those integrals do not exist",__FILE__,__LINE__);

  if (trans == 'N' || trans == 'n') {
    double **T_p_OV = block_matrix(aoccA*nvirA,ndf_+3);

    C_DGEMM('N','N',aoccA*nvirA,ndf_+3,aoccB*nvirB,1.0,tARBS[0],aoccB*nvirB,
      B_p_OV[0],ndf_+3,0.0,T_p_OV[0],ndf_+3);

    psio_->write_entry(ampout,thetalabel,(char *) T_p_OV[0],
      sizeof(double)*aoccA*nvirA*(ndf_+3));

    free_block(T_p_OV);
  }
  else if (trans == 'T' || trans == 't') {
    double **T_p_OV = block_matrix(aoccB*nvirB,ndf_+3);

    C_DGEMM('T','N',aoccB*nvirB,ndf_+3,aoccA*nvirA,1.0,tARBS[0],aoccB*nvirB,
      B_p_OV[0],ndf_+3,0.0,T_p_OV[0],ndf_+3);

    psio_->write_entry(ampout,thetalabel,(char *) T_p_OV[0],
      sizeof(double)*aoccB*nvirB*(ndf_+3));

    free_block(T_p_OV);
  }
  else
    throw PsiException("You want me to do what to that matrix?",
       __FILE__,__LINE__);

  free_block(tARBS);
  free_block(B_p_OV);
}

void SAPT2::Y2(int intfile, const char *AAlabel, const char *ARlabel,
  const char *RRlabel, int ampfile, const char *pAAlabel, const char *pRRlabel,
  const char *thetalabel, int foccA, int noccA, int nvirA, double *evals,
  int ampout, const char *Ylabel, const char *tlabel)
{
  int aoccA = noccA - foccA;

  double **yAR = block_matrix(aoccA,nvirA);
  double **tAR = block_matrix(aoccA,nvirA);

  Y2_3(yAR,intfile,AAlabel,RRlabel,ampfile,thetalabel,foccA,noccA,nvirA);

  C_DCOPY(aoccA*nvirA,yAR[0],1,tAR[0],1);

  for (int a=0; a<aoccA; a++) {
    for (int r=0; r<nvirA; r++) {
      tAR[a][r] /= evals[a+foccA] - evals[r+noccA];
  }}

  psio_->write_entry(ampfile,tlabel,(char *) tAR[0],
    sizeof(double)*aoccA*nvirA);

  free_block(tAR);

  Y2_1(yAR,intfile,ARlabel,RRlabel,ampfile,pRRlabel,foccA,noccA,nvirA);
  Y2_2(yAR,intfile,AAlabel,ARlabel,ampfile,pAAlabel,foccA,noccA,nvirA);

  psio_->write_entry(ampout,Ylabel,(char *) yAR[0],
    sizeof(double)*aoccA*nvirA);

  free_block(yAR);
}

void SAPT2::Y2_1(double **yAR, int intfile, const char *ARlabel,
  const char *RRlabel, int ampfile, const char *pRRlabel, int foccA,
  int noccA, int nvirA)
{
  int aoccA = noccA - foccA;

  double **pRR = block_matrix(nvirA,nvirA);
  psio_->read_entry(ampfile,pRRlabel,(char *) pRR[0],
    sizeof(double)*nvirA*nvirA);

  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);
  double **B_p_RR = get_DF_ints(intfile,RRlabel,0,nvirA,0,nvirA);

  double *X = init_array(ndf_+3);
  double **C_p_AR = block_matrix(aoccA*nvirA,ndf_+3);

  C_DGEMV('t',nvirA*nvirA,ndf_+3,1.0,B_p_RR[0],ndf_+3,pRR[0],1,0.0,X,1);

  for (int a=0; a<aoccA; a++) {
    C_DGEMM('T','N',nvirA,ndf_+3,nvirA,1.0,pRR[0],nvirA,B_p_AR[a*nvirA],ndf_+3,
      0.0,C_p_AR[a*nvirA],ndf_+3);
  }

  C_DGEMV('n',aoccA*nvirA,ndf_+3,2.0,B_p_AR[0],ndf_+3,X,1,1.0,yAR[0],1);

  C_DGEMM('N','T',aoccA,nvirA,nvirA*(ndf_+3),-1.0,C_p_AR[0],nvirA*(ndf_+3),
    B_p_RR[0],nvirA*(ndf_+3),1.0,yAR[0],nvirA);

  free(X);
  free_block(pRR);
  free_block(B_p_AR);
  free_block(C_p_AR);
  free_block(B_p_RR);
}

void SAPT2::Y2_2(double **yAR, int intfile, const char *AAlabel,
  const char *ARlabel, int ampfile, const char *pAAlabel, int foccA,
  int noccA, int nvirA)
{
  int aoccA = noccA - foccA;

  double **pAA = block_matrix(aoccA,aoccA);
  psio_->read_entry(ampfile,pAAlabel,(char *) pAA[0],
    sizeof(double)*aoccA*aoccA);

  double **B_p_AA = get_DF_ints(intfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);

  double *X = init_array(ndf_+3);
  double **C_p_AA = block_matrix(aoccA*aoccA,ndf_+3);

  C_DGEMV('t',aoccA*aoccA,ndf_+3,1.0,B_p_AA[0],ndf_+3,pAA[0],1,0.0,X,1);

  C_DGEMM('N','N',aoccA,aoccA*(ndf_+3),aoccA,1.0,pAA[0],aoccA,B_p_AA[0],
    aoccA*(ndf_+3),0.0,C_p_AA[0],aoccA*(ndf_+3));

  C_DGEMV('n',aoccA*nvirA,ndf_+3,-2.0,B_p_AR[0],ndf_+3,X,1,1.0,yAR[0],1);

  for (int a=0; a<aoccA; a++) {
    C_DGEMM('N','T',aoccA,nvirA,ndf_+3,1.0,C_p_AA[a*aoccA],ndf_+3,
      B_p_AR[a*nvirA],ndf_+3,1.0,yAR[0],nvirA);
  }

  free(X);
  free_block(pAA);
  free_block(B_p_AA);
  free_block(C_p_AA);
  free_block(B_p_AR);
}

void SAPT2::Y2_3(double **yAR, int intfile, const char *AAlabel,
  const char *RRlabel, int ampfile, const char *thetalabel, int foccA,
  int noccA, int nvirA)
{
  int aoccA = noccA - foccA;

  double **T_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  psio_->read_entry(ampfile,thetalabel,(char *) T_p_AR[0],
    sizeof(double)*aoccA*nvirA*(ndf_+3));

  double **B_p_AA = get_DF_ints(intfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_RR = get_DF_ints(intfile,RRlabel,0,nvirA,0,nvirA);

  for (int a=0; a<aoccA; a++) {
    C_DGEMM('N','T',aoccA,nvirA,ndf_+3,-1.0,B_p_AA[a*aoccA],ndf_+3,
      T_p_AR[a*nvirA],ndf_+3,1.0,yAR[0],nvirA);
  }

  C_DGEMM('N','T',aoccA,nvirA,nvirA*(ndf_+3),1.0,T_p_AR[0],nvirA*(ndf_+3),
    B_p_RR[0],nvirA*(ndf_+3),1.0,yAR[0],nvirA);

  free_block(B_p_AA);
  free_block(T_p_AR);
  free_block(B_p_RR);
}

void SAPT2::t2OVOV(int ampfile, const char *tlabel, const char *thetalabel,
  int intfile, const char *AAlabel, const char *ARlabel, const char *RRlabel,
  int foccA, int noccA, int nvirA, double *evalsA, int ampout,
  const char *t2label)
{
  int aoccA = noccA - foccA;

  double *t2ARAR = init_array((long int) aoccA*nvirA*aoccA*nvirA);

  double **vAARR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  double **B_p_AA = get_DF_ints(intfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_RR = get_DF_ints(intfile,RRlabel,0,nvirA,0,nvirA);

  for (int a=0, ar=0; a<aoccA; a++) {
    for (int r=0; r<nvirA; r++, ar++) {
      C_DGEMM('N','T',aoccA,nvirA,ndf_+3,1.0,B_p_AA[a*aoccA],ndf_+3,
        B_p_RR[r*nvirA],ndf_+3,0.0,vAARR[ar],nvirA);
  }}

  free_block(B_p_AA);
  free_block(B_p_RR);

  double *tARAR = init_array((long int) aoccA*nvirA*aoccA*nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) tARAR,
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);

  OVOpVp_to_OVpOpV(tARAR,aoccA,nvirA);

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,aoccA*nvirA,-1.0,tARAR,aoccA*nvirA,
    vAARR[0],aoccA*nvirA,0.0,t2ARAR,aoccA*nvirA);

  OVOpVp_to_OVpOpV(tARAR,aoccA,nvirA);
  OVOpVp_to_OVpOpV(t2ARAR,aoccA,nvirA);

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,aoccA*nvirA,-1.0,tARAR,aoccA*nvirA,
    vAARR[0],aoccA*nvirA,1.0,t2ARAR,aoccA*nvirA);

  free_block(vAARR);

  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);
  double **T_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  psio_->read_entry(ampfile,thetalabel,(char *) T_p_AR[0],
    sizeof(double)*aoccA*nvirA*(ndf_+3));

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,ndf_+3,1.0,B_p_AR[0],ndf_+3,
    T_p_AR[0],ndf_+3,1.0,t2ARAR,aoccA*nvirA);

  free_block(B_p_AR);
  free_block(T_p_AR);

  ijkl_to_ikjl(tARAR,aoccA,nvirA,aoccA,nvirA);
  ijkl_to_ikjl(t2ARAR,aoccA,nvirA,aoccA,nvirA);

  double **vAAAA = block_matrix(aoccA*aoccA,aoccA*aoccA);
  B_p_AA = get_DF_ints(intfile,AAlabel,foccA,noccA,foccA,noccA);

  for (int a=0, aa1=0; a<aoccA; a++) {
    for (int a1=0; a1<aoccA; a1++, aa1++) {
      C_DGEMM('N','T',aoccA,aoccA,ndf_+3,1.0,B_p_AA[a*aoccA],ndf_+3,
        B_p_AA[a1*aoccA],ndf_+3,0.0,vAAAA[aa1],aoccA);
  }}

  free_block(B_p_AA);

  C_DGEMM('N','N',aoccA*aoccA,nvirA*nvirA,aoccA*aoccA,0.5,vAAAA[0],aoccA*aoccA,
    tARAR,nvirA*nvirA,1.0,t2ARAR,nvirA*nvirA);

  free_block(vAAAA);

  B_p_RR = get_DF_ints(intfile,RRlabel,0,nvirA,0,nvirA);
  double **xRRR = block_matrix(nvirA*nvirA,nvirA);

  for (int r=0; r < nvirA; r++) {
    C_DGEMM('N','T',nvirA,nvirA*nvirA,ndf_+3,1.0,&(B_p_RR[r*nvirA][0]),
      ndf_+3,&(B_p_RR[0][0]),ndf_+3,0.0,&(xRRR[0][0]),nvirA*nvirA);
    C_DGEMM('N','T',aoccA*aoccA,nvirA*nvirA,nvirA,0.5,&(tARAR[r*nvirA]),
      nvirA*nvirA,&(xRRR[0][0]),nvirA,1.0,&(t2ARAR[0]),nvirA*nvirA);
  }

  free(tARAR);
  free_block(B_p_RR);
  free_block(xRRR);

  ijkl_to_ikjl(t2ARAR,aoccA,aoccA,nvirA,nvirA);
  symmetrize(t2ARAR,aoccA,nvirA);

  for(int a=0; a<aoccA; a++) {
    for(int r=0; r<nvirA; r++) {
      for(int a1=0; a1<aoccA; a1++) {
        for(int r1=0; r1<nvirA; r1++) {
          long int ar = a*nvirA + r;
          long int a1r1 = a1*nvirA + r1;
          long int ara1r1 = ar*aoccA*nvirA + a1r1;
          t2ARAR[ara1r1] /= evalsA[a+foccA] + evalsA[a1+foccA]
            - evalsA[r+noccA] - evalsA[r1+noccA];
  }}}}

  psio_->write_entry(ampout,t2label,(char *) t2ARAR,sizeof(double)*aoccA*
    nvirA*aoccA*nvirA);

  free(t2ARAR);
}

void SAPT2::t2OVOV(int ampfile, const char *tlabel, const char *no_tlabel,
  const char *thetalabel, int intfile, const char *AAlabel,
  const char *ARlabel, const char *RRlabel, const char *no_RRlabel, int foccA,
  int noccA, int nvirA, int no_nvirA, double *evalsA, double **CA, int ampout,
  const char *t2label)
{
  int aoccA = noccA - foccA;

  double *t2ARAR = init_array((long int) aoccA*nvirA*aoccA*nvirA);

  double **vAARR = block_matrix(aoccA*nvirA,aoccA*nvirA);
  double **B_p_AA = get_DF_ints(intfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_RR = get_DF_ints(intfile,RRlabel,0,nvirA,0,nvirA);

  for (int a=0, ar=0; a<aoccA; a++) {
    for (int r=0; r<nvirA; r++, ar++) {
      C_DGEMM('N','T',aoccA,nvirA,ndf_+3,1.0,B_p_AA[a*aoccA],ndf_+3,
        B_p_RR[r*nvirA],ndf_+3,0.0,vAARR[ar],nvirA);
  }}

  free_block(B_p_AA);
  free_block(B_p_RR);

  double *tARAR = init_array((long int) aoccA*nvirA*aoccA*nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) tARAR,
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);

  OVOpVp_to_OVpOpV(tARAR,aoccA,nvirA);

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,aoccA*nvirA,-1.0,tARAR,aoccA*nvirA,
    vAARR[0],aoccA*nvirA,0.0,t2ARAR,aoccA*nvirA);

  OVOpVp_to_OVpOpV(tARAR,aoccA,nvirA);
  OVOpVp_to_OVpOpV(t2ARAR,aoccA,nvirA);

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,aoccA*nvirA,-1.0,tARAR,aoccA*nvirA,
    vAARR[0],aoccA*nvirA,1.0,t2ARAR,aoccA*nvirA);

  free_block(vAARR);

  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);
  double **T_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  psio_->read_entry(ampfile,thetalabel,(char *) T_p_AR[0],
    sizeof(double)*aoccA*nvirA*(ndf_+3));

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,ndf_+3,1.0,B_p_AR[0],ndf_+3,
    T_p_AR[0],ndf_+3,1.0,t2ARAR,aoccA*nvirA);

  free_block(B_p_AR);
  free_block(T_p_AR);

  ijkl_to_ikjl(tARAR,aoccA,nvirA,aoccA,nvirA);
  ijkl_to_ikjl(t2ARAR,aoccA,nvirA,aoccA,nvirA);

  double **vAAAA = block_matrix(aoccA*aoccA,aoccA*aoccA);
  B_p_AA = get_DF_ints(intfile,AAlabel,foccA,noccA,foccA,noccA);

  for (int a=0, aa1=0; a<aoccA; a++) {
    for (int a1=0; a1<aoccA; a1++, aa1++) {
      C_DGEMM('N','T',aoccA,aoccA,ndf_+3,1.0,B_p_AA[a*aoccA],ndf_+3,
        B_p_AA[a1*aoccA],ndf_+3,0.0,vAAAA[aa1],aoccA);
  }}

  free_block(B_p_AA);

  C_DGEMM('N','N',aoccA*aoccA,nvirA*nvirA,aoccA*aoccA,0.5,vAAAA[0],aoccA*aoccA,
    tARAR,nvirA*nvirA,1.0,t2ARAR,nvirA*nvirA);

  free(tARAR);
  free_block(vAAAA);

  double **t2AArr = block_matrix(aoccA*aoccA,no_nvirA*no_nvirA);
  double *xRR = init_array(nvirA*no_nvirA);

  for (int a=0,aaa=0; a<aoccA; a++) {
    for (int aa=0; aa<aoccA; aa++,aaa++) {
      long int aarr = (long int) aaa*nvirA*nvirA;
      C_DGEMM('N','N',nvirA,no_nvirA,nvirA,1.0,&(t2ARAR[aarr]),nvirA,CA[0],
        no_nvirA,0.0,xRR,no_nvirA);
      C_DGEMM('T','N',no_nvirA,no_nvirA,nvirA,1.0,CA[0],no_nvirA,
        xRR,no_nvirA,0.0,t2AArr[aaa],no_nvirA);
  }}

  free(t2ARAR);

  double *tArAr = init_array((long int) aoccA*aoccA*no_nvirA*no_nvirA);
  psio_->read_entry(ampfile,no_tlabel,(char *) tArAr,
    sizeof(double)*aoccA*no_nvirA*aoccA*no_nvirA);
  ijkl_to_ikjl(tArAr,aoccA,no_nvirA,aoccA,no_nvirA);

  B_p_RR = get_DF_ints(intfile,no_RRlabel,0,no_nvirA,0,no_nvirA);
  double **xRRR = block_matrix(no_nvirA*no_nvirA,no_nvirA);

  for (int r=0; r < no_nvirA; r++) {
    C_DGEMM('N','T',no_nvirA,no_nvirA*no_nvirA,ndf_+3,1.0,
      &(B_p_RR[r*no_nvirA][0]),ndf_+3,&(B_p_RR[0][0]),ndf_+3,0.0,
      &(xRRR[0][0]),no_nvirA*no_nvirA);
    C_DGEMM('N','T',aoccA*aoccA,no_nvirA*no_nvirA,no_nvirA,0.5,
      &(tArAr[r*no_nvirA]),no_nvirA*no_nvirA,&(xRRR[0][0]),no_nvirA,
      1.0,&(t2AArr[0][0]),no_nvirA*no_nvirA);
  }

  free(tArAr);
  free_block(B_p_RR);
  free_block(xRRR);

  t2ARAR = init_array((long int) aoccA*aoccA*nvirA*nvirA);

  for (int a=0,aaa=0; a<aoccA; a++) {
    for (int aa=0; aa<aoccA; aa++,aaa++) {
      long int aarr = (long int) aaa*nvirA*nvirA;
      C_DGEMM('N','N',nvirA,no_nvirA,no_nvirA,1.0,CA[0],no_nvirA,
        t2AArr[aaa],no_nvirA,0.0,xRR,no_nvirA);
      C_DGEMM('N','T',nvirA,nvirA,no_nvirA,1.0,xRR,no_nvirA,
        CA[0],no_nvirA,0.0,&(t2ARAR[aarr]),nvirA);
  }}

  free(xRR);
  free_block(t2AArr);

  ijkl_to_ikjl(t2ARAR,aoccA,aoccA,nvirA,nvirA);
  symmetrize(t2ARAR,aoccA,nvirA);

  for(int a=0; a<aoccA; a++) {
    for(int r=0; r<nvirA; r++) {
      for(int a1=0; a1<aoccA; a1++) {
        for(int r1=0; r1<nvirA; r1++) {
          long int ar = a*nvirA + r;
          long int a1r1 = a1*nvirA + r1;
          long int ara1r1 = ar*aoccA*nvirA + a1r1;
          t2ARAR[ara1r1] /= evalsA[a+foccA] + evalsA[a1+foccA]
            - evalsA[r+noccA] - evalsA[r1+noccA];
  }}}}

  psio_->write_entry(ampout,t2label,(char *) t2ARAR,sizeof(double)*aoccA*
    nvirA*aoccA*nvirA);

  free(t2ARAR);
}

void SAPT2::OVOpVp_to_OVpOpV(double *tARAR, int nocc, int nvir)
{
  for(int a=0; a<nocc; a++) {
    for(int r=0; r<nvir; r++) {
      for(int a1=0; a1<a; a1++) {
        for(int r1=0; r1<nvir; r1++) {
          long int ar = a*nvir + r;
          long int a1r1 = a1*nvir + r1;
          long int a1r = a1*nvir + r;
          long int ar1 = a*nvir + r1;
          long int ara1r1 = ar*nocc*nvir + a1r1;
          long int a1rar1 = a1r*nocc*nvir + ar1;
          double tval = tARAR[ara1r1];
          tARAR[ara1r1] = tARAR[a1rar1];
          tARAR[a1rar1] = tval;
  }}}}
}

void SAPT2::ijkl_to_ikjl(double *tARAR, int ilength, int jlength, int klength,
  int llength)
{
  double *X = init_array(jlength*klength);

  for(int i=0; i<ilength; i++) {
    for(int l=0; l<llength; l++) {
      long int ijkl = i*jlength*klength*llength + l;
      C_DCOPY(jlength*klength,&(tARAR[ijkl]),llength,X,1);
      for(int j=0; j<jlength; j++) {
        for(int k=0; k<klength; k++) {
          long int ik = i*klength + k;
          long int jl = j*llength + l;
          long int jk = j*klength + k;
          long int ikjl = ik*jlength*llength + jl;
          tARAR[ikjl] = X[jk];
      }}
  }}

  free(X);
}

void SAPT2::symmetrize(double *tARAR, int nocc, int nvir)
{
  for(int ar=0; ar<nocc*nvir; ar++) {
    for (int a1r1=0; a1r1<=ar; a1r1++) {
      long int ara1r1 = ar*nocc*nvir + a1r1;
      long int a1r1ar = a1r1*nocc*nvir + ar;
      tARAR[ara1r1] = tARAR[a1r1ar] = tARAR[ara1r1] + tARAR[a1r1ar];
  }}
}

void SAPT2p::amplitudes()
{
  tOVOV(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",foccA_,noccA_,nvirA_,evalsA_,
    PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",foccA_,noccA_,nvirA_,evalsA_,
    PSIF_SAPT_AMPS,"tARAR Amplitudes");
  tOVOV(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",foccB_,noccB_,nvirB_,evalsB_,
    PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",foccB_,noccB_,nvirB_,evalsB_,
    PSIF_SAPT_AMPS,"tBSBS Amplitudes");

  tOVOV(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",foccA_,noccA_,nvirA_,evalsA_,
    PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",foccB_,noccB_,nvirB_,evalsB_,
    PSIF_SAPT_AMPS,"tARBS Amplitudes");

  pOOpVV(PSIF_SAPT_AMPS,"tARAR Amplitudes","tARAR Amplitudes",aoccA_,nvirA_,
    PSIF_SAPT_AMPS,"pAA Density Matrix","pRR Density Matrix");
  pOOpVV(PSIF_SAPT_AMPS,"tBSBS Amplitudes","tBSBS Amplitudes",aoccB_,nvirB_,
    PSIF_SAPT_AMPS,"pBB Density Matrix","pSS Density Matrix");

  if (nat_orbs_) {
    natural_orbitalify(PSIF_SAPT_AMPS,"pRR Density Matrix",evalsA_,foccA_,
      noccA_,nvirA_,'A');
    natural_orbitalify(PSIF_SAPT_AMPS,"pSS Density Matrix",evalsB_,foccB_,
      noccB_,nvirB_,'B');
    natural_orbitalify_df_ints();
    tOVOV(PSIF_SAPT_AA_DF_INTS,"AR NO RI Integrals",foccA_,noccA_,no_nvirA_,
      no_evalsA_,PSIF_SAPT_AA_DF_INTS,"AR NO RI Integrals",foccA_,noccA_,
      no_nvirA_,no_evalsA_,PSIF_SAPT_AMPS,"tARAR NO Amplitudes");
    tOVOV(PSIF_SAPT_BB_DF_INTS,"BS NO RI Integrals",foccB_,noccB_,no_nvirB_,
      no_evalsB_,PSIF_SAPT_BB_DF_INTS,"BS NO RI Integrals",foccB_,noccB_,
      no_nvirB_,no_evalsB_,PSIF_SAPT_AMPS,"tBSBS NO Amplitudes");
    if (print_) fprintf(outfile,"\n");
  }

  theta(PSIF_SAPT_AMPS,"tARAR Amplitudes",'N',true,aoccA_,nvirA_,aoccA_,nvirA_,
    "AR RI Integrals",PSIF_SAPT_AMPS,"Theta AR Intermediates");
  theta(PSIF_SAPT_AMPS,"tBSBS Amplitudes",'N',true,aoccB_,nvirB_,aoccB_,nvirB_,
    "BS RI Integrals",PSIF_SAPT_AMPS,"Theta BS Intermediates");

  theta(PSIF_SAPT_AMPS,"tARBS Amplitudes",'N',false,aoccA_,nvirA_,
    aoccB_,nvirB_,"BS RI Integrals",PSIF_SAPT_AMPS,"T AR Intermediates");
  theta(PSIF_SAPT_AMPS,"tARBS Amplitudes",'T',false,aoccA_,nvirA_,
    aoccB_,nvirB_,"AR RI Integrals",PSIF_SAPT_AMPS,"T BS Intermediates");

  Y2(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
    "RR RI Integrals",PSIF_SAPT_AMPS,"pAA Density Matrix","pRR Density Matrix",
    "Theta AR Intermediates",foccA_,noccA_,nvirA_,evalsA_,PSIF_SAPT_AMPS,
    "Y2 AR Amplitudes","T2 AR Amplitudes");
  Y2(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
    "SS RI Integrals",PSIF_SAPT_AMPS,"pBB Density Matrix","pSS Density Matrix",
    "Theta BS Intermediates",foccB_,noccB_,nvirB_,evalsB_,PSIF_SAPT_AMPS,
    "Y2 BS Amplitudes","T2 BS Amplitudes");

  if (nat_orbs_t2_) {
    t2OVOV(PSIF_SAPT_AMPS,"tARAR Amplitudes","tARAR NO Amplitudes",
      "Theta AR Intermediates",PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
      "AR RI Integrals","RR RI Integrals","RR NO RI Integrals",
      foccA_,noccA_,nvirA_,no_nvirA_,evalsA_,no_CA_,PSIF_SAPT_AMPS,
      "t2ARAR Amplitudes");
    t2OVOV(PSIF_SAPT_AMPS,"tBSBS Amplitudes","tBSBS NO Amplitudes",
      "Theta BS Intermediates",PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
      "BS RI Integrals","SS RI Integrals","SS NO RI Integrals",
      foccB_,noccB_,nvirB_,no_nvirB_,evalsB_,no_CB_,PSIF_SAPT_AMPS,
      "t2BSBS Amplitudes");
  }
  else {
    t2OVOV(PSIF_SAPT_AMPS,"tARAR Amplitudes","Theta AR Intermediates",
      PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
      "RR RI Integrals",foccA_,noccA_,nvirA_,evalsA_,PSIF_SAPT_AMPS,
      "t2ARAR Amplitudes");
    t2OVOV(PSIF_SAPT_AMPS,"tBSBS Amplitudes","Theta BS Intermediates",
      PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
      "SS RI Integrals",foccB_,noccB_,nvirB_,evalsB_,PSIF_SAPT_AMPS,
      "t2BSBS Amplitudes");
  }

  theta(PSIF_SAPT_AMPS,"t2ARAR Amplitudes",'N',true,aoccA_,nvirA_,aoccA_,
    nvirA_,"AR RI Integrals",PSIF_SAPT_AMPS,"Theta 2 AR Intermediates");
  theta(PSIF_SAPT_AMPS,"t2BSBS Amplitudes",'N',true,aoccB_,nvirB_,aoccB_,
    nvirB_,"BS RI Integrals",PSIF_SAPT_AMPS,"Theta 2 BS Intermediates");

  gARARxtARBS(PSIF_SAPT_AMPS,"tARBS Amplitudes",'N',PSIF_SAPT_AA_DF_INTS,
    "AA RI Integrals","AR RI Integrals","RR RI Integrals",foccA_,noccA_,
    nvirA_,foccB_,noccB_,nvirB_,PSIF_SAPT_AMPS,"gARAR x tARBS");
  gARARxtARBS(PSIF_SAPT_AMPS,"tARBS Amplitudes",'T',PSIF_SAPT_BB_DF_INTS,
    "BB RI Integrals","BS RI Integrals","SS RI Integrals",foccB_,noccB_,
    nvirB_,foccA_,noccA_,nvirA_,PSIF_SAPT_AMPS,"gBSBS x tARBS");
}

void SAPT2p::gARARxtARBS(int ampfile, const char *tlabel, const char trans,
  int intfile, const char *AAlabel, const char *ARlabel, const char *RRlabel,
  int foccA, int noccA, int nvirA, int foccB, int noccB, int nvirB, int ampout,
  const char *labelout)
{
  int aoccA = noccA - foccA;
  int aoccB = noccB - foccB;

  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);
  double **gARAR = block_matrix(aoccA*nvirA,aoccA*nvirA);

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,ndf_+3,2.0,B_p_AR[0],ndf_+3,
    B_p_AR[0],ndf_+3,0.0,gARAR[0],aoccA*nvirA);

  free_block(B_p_AR);

  double **B_p_AA = get_DF_ints(intfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_RR = get_DF_ints(intfile,RRlabel,0,nvirA,0,nvirA);

  for (int a=0, ar=0; a<aoccA; a++){
    for (int r=0; r<nvirA; r++, ar++){
      C_DGEMM('N','T',aoccA,nvirA,ndf_+3,-1.0,B_p_AA[a*aoccA],ndf_+3,
        B_p_RR[r*nvirA],ndf_+3,1.0,gARAR[ar],nvirA);
  }}

  free_block(B_p_AA);
  free_block(B_p_RR);

  double **tARBS;
  double **xARBS;

  if (trans == 'n' || trans == 'N') {
    tARBS = block_matrix(aoccA*nvirA,aoccB*nvirB);
    xARBS = block_matrix(aoccA*nvirA,aoccB*nvirB);
    psio_->read_entry(ampfile,tlabel,(char *) tARBS[0],
      sizeof(double)*aoccA*nvirA*aoccB*nvirB);

    C_DGEMM('N','N',aoccA*nvirA,aoccB*nvirB,aoccA*nvirA,1.0,gARAR[0],
      aoccA*nvirA,tARBS[0],aoccB*nvirB,0.0,xARBS[0],aoccB*nvirB);
  }
  else if (trans == 't' || trans == 'T') {
    tARBS = block_matrix(aoccB*nvirB,aoccA*nvirA);
    xARBS = block_matrix(aoccB*nvirB,aoccA*nvirA);
    psio_->read_entry(ampfile,tlabel,(char *) tARBS[0],
      sizeof(double)*aoccA*nvirA*aoccB*nvirB);

    C_DGEMM('N','N',aoccB*nvirB,aoccA*nvirA,aoccA*nvirA,1.0,tARBS[0],
      aoccA*nvirA,gARAR[0],aoccA*nvirA,0.0,xARBS[0],aoccA*nvirA);
  }
  else
    throw PsiException("You want me to do what to that matrix?",
       __FILE__,__LINE__);

  psio_->write_entry(ampout,labelout,(char *) xARBS[0],
    sizeof(double)*aoccA*nvirA*aoccB*nvirB);

  free_block(gARAR);
  free_block(tARBS);
  free_block(xARBS);
}

void SAPT2p3::amplitudes()
{
  tOVOV(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",foccA_,noccA_,nvirA_,evalsA_,
    PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",foccA_,noccA_,nvirA_,evalsA_,
    PSIF_SAPT_AMPS,"tARAR Amplitudes");
  tOVOV(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",foccB_,noccB_,nvirB_,evalsB_,
    PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",foccB_,noccB_,nvirB_,evalsB_,
    PSIF_SAPT_AMPS,"tBSBS Amplitudes");

  tOVOV(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",foccA_,noccA_,nvirA_,evalsA_,
    PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",foccB_,noccB_,nvirB_,evalsB_,
    PSIF_SAPT_AMPS,"tARBS Amplitudes");

  pOOpVV(PSIF_SAPT_AMPS,"tARAR Amplitudes","tARAR Amplitudes",aoccA_,nvirA_,
    PSIF_SAPT_AMPS,"pAA Density Matrix","pRR Density Matrix");
  pOOpVV(PSIF_SAPT_AMPS,"tBSBS Amplitudes","tBSBS Amplitudes",aoccB_,nvirB_,
    PSIF_SAPT_AMPS,"pBB Density Matrix","pSS Density Matrix");

  if (nat_orbs_) {
    natural_orbitalify(PSIF_SAPT_AMPS,"pRR Density Matrix",evalsA_,foccA_,
      noccA_,nvirA_,'A');
    natural_orbitalify(PSIF_SAPT_AMPS,"pSS Density Matrix",evalsB_,foccB_,
      noccB_,nvirB_,'B');
    natural_orbitalify_df_ints();
    tOVOV(PSIF_SAPT_AA_DF_INTS,"AR NO RI Integrals",foccA_,noccA_,no_nvirA_,
      no_evalsA_,PSIF_SAPT_AA_DF_INTS,"AR NO RI Integrals",foccA_,noccA_,
      no_nvirA_,no_evalsA_,PSIF_SAPT_AMPS,"tARAR NO Amplitudes");
    tOVOV(PSIF_SAPT_BB_DF_INTS,"BS NO RI Integrals",foccB_,noccB_,no_nvirB_,
      no_evalsB_,PSIF_SAPT_BB_DF_INTS,"BS NO RI Integrals",foccB_,noccB_,
      no_nvirB_,no_evalsB_,PSIF_SAPT_AMPS,"tBSBS NO Amplitudes");
    if (print_) fprintf(outfile,"\n");
  }

  theta(PSIF_SAPT_AMPS,"tARAR Amplitudes",'N',true,aoccA_,nvirA_,aoccA_,nvirA_,
    "AR RI Integrals",PSIF_SAPT_AMPS,"Theta AR Intermediates");
  theta(PSIF_SAPT_AMPS,"tBSBS Amplitudes",'N',true,aoccB_,nvirB_,aoccB_,nvirB_,
    "BS RI Integrals",PSIF_SAPT_AMPS,"Theta BS Intermediates");

  theta(PSIF_SAPT_AMPS,"tARBS Amplitudes",'N',false,aoccA_,nvirA_,
    aoccB_,nvirB_,"BS RI Integrals",PSIF_SAPT_AMPS,"T AR Intermediates");
  theta(PSIF_SAPT_AMPS,"tARBS Amplitudes",'T',false,aoccA_,nvirA_,
    aoccB_,nvirB_,"AR RI Integrals",PSIF_SAPT_AMPS,"T BS Intermediates");

  Y2(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
    "RR RI Integrals",PSIF_SAPT_AMPS,"pAA Density Matrix","pRR Density Matrix",
    "Theta AR Intermediates",foccA_,noccA_,nvirA_,evalsA_,PSIF_SAPT_AMPS,
    "Y2 AR Amplitudes","T2 AR Amplitudes");
  Y2(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
    "SS RI Integrals",PSIF_SAPT_AMPS,"pBB Density Matrix","pSS Density Matrix",
    "Theta BS Intermediates",foccB_,noccB_,nvirB_,evalsB_,PSIF_SAPT_AMPS,
    "Y2 BS Amplitudes","T2 BS Amplitudes");

  if (nat_orbs_t2_) {
    t2OVOV(PSIF_SAPT_AMPS,"tARAR Amplitudes","tARAR NO Amplitudes",
      "Theta AR Intermediates",PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
      "AR RI Integrals","RR RI Integrals","RR NO RI Integrals",
      foccA_,noccA_,nvirA_,no_nvirA_,evalsA_,no_CA_,PSIF_SAPT_AMPS,
      "t2ARAR Amplitudes");
    t2OVOV(PSIF_SAPT_AMPS,"tBSBS Amplitudes","tBSBS NO Amplitudes",
      "Theta BS Intermediates",PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
      "BS RI Integrals","SS RI Integrals","SS NO RI Integrals",
      foccB_,noccB_,nvirB_,no_nvirB_,evalsB_,no_CB_,PSIF_SAPT_AMPS,
      "t2BSBS Amplitudes");
  }
  else {
    t2OVOV(PSIF_SAPT_AMPS,"tARAR Amplitudes","Theta AR Intermediates",
      PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
      "RR RI Integrals",foccA_,noccA_,nvirA_,evalsA_,PSIF_SAPT_AMPS,
      "t2ARAR Amplitudes");
    t2OVOV(PSIF_SAPT_AMPS,"tBSBS Amplitudes","Theta BS Intermediates",
      PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
      "SS RI Integrals",foccB_,noccB_,nvirB_,evalsB_,PSIF_SAPT_AMPS,
      "t2BSBS Amplitudes");
  }

  theta(PSIF_SAPT_AMPS,"t2ARAR Amplitudes",'N',true,aoccA_,nvirA_,aoccA_,
    nvirA_,"AR RI Integrals",PSIF_SAPT_AMPS,"Theta 2 AR Intermediates");
  theta(PSIF_SAPT_AMPS,"t2BSBS Amplitudes",'N',true,aoccB_,nvirB_,aoccB_,
    nvirB_,"BS RI Integrals",PSIF_SAPT_AMPS,"Theta 2 BS Intermediates");

  gARARxtARBS(PSIF_SAPT_AMPS,"tARBS Amplitudes",'N',PSIF_SAPT_AA_DF_INTS,
    "AA RI Integrals","AR RI Integrals","RR RI Integrals",foccA_,noccA_,
    nvirA_,foccB_,noccB_,nvirB_,PSIF_SAPT_AMPS,"gARAR x tARBS");
  gARARxtARBS(PSIF_SAPT_AMPS,"tARBS Amplitudes",'T',PSIF_SAPT_BB_DF_INTS,
    "BB RI Integrals","BS RI Integrals","SS RI Integrals",foccB_,noccB_,
    nvirB_,foccA_,noccA_,nvirA_,PSIF_SAPT_AMPS,"gBSBS x tARBS");

  pOOpVV(PSIF_SAPT_AMPS,"tARAR Amplitudes","t2ARAR Amplitudes",aoccA_,nvirA_,
    PSIF_SAPT_AMPS,"qAA Density Matrix","qRR Density Matrix");
  pOOpVV(PSIF_SAPT_AMPS,"tBSBS Amplitudes","t2BSBS Amplitudes",aoccB_,nvirB_,
    PSIF_SAPT_AMPS,"qBB Density Matrix","qSS Density Matrix");

  pOOpVV(PSIF_SAPT_AMPS,"t2ARAR Amplitudes","tARAR Amplitudes",aoccA_,nvirA_,
    PSIF_SAPT_AMPS,"qbarAA Density Matrix","qbarRR Density Matrix");
  pOOpVV(PSIF_SAPT_AMPS,"t2BSBS Amplitudes","tBSBS Amplitudes",aoccB_,nvirB_,
    PSIF_SAPT_AMPS,"qbarBB Density Matrix","qbarSS Density Matrix");

  Y3(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","AR RI Integrals",
    "RR RI Integrals",PSIF_SAPT_AMPS,"qAA Density Matrix",
    "qRR Density Matrix","qbarAA Density Matrix","qbarRR Density Matrix",
    "Theta 2 AR Intermediates","tARAR Amplitudes",foccA_,noccA_,nvirA_,
    evalsA_,PSIF_SAPT_AMPS,"Y3 AR Amplitudes");
  Y3(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","BS RI Integrals",
    "SS RI Integrals",PSIF_SAPT_AMPS,"qBB Density Matrix",
    "qSS Density Matrix","qbarBB Density Matrix","qbarSS Density Matrix",
    "Theta 2 BS Intermediates","tBSBS Amplitudes",foccB_,noccB_,nvirB_,
    evalsB_,PSIF_SAPT_AMPS,"Y3 BS Amplitudes");

  if (third_order_) {

    ind30_amps(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",PSIF_SAPT_BB_DF_INTS,
      "BS RI Integrals",wBAA_,wBAR_,wBRR_,wABS_,noccA_,nvirA_,evalsA_,
      noccB_,nvirB_,evalsB_,PSIF_SAPT_AMPS,"Ind30 uAR Amplitudes");
    ind30_amps(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",PSIF_SAPT_AA_DF_INTS,
      "AR RI Integrals",wABB_,wABS_,wASS_,wBAR_,noccB_,nvirB_,evalsB_,
      noccA_,nvirA_,evalsA_,PSIF_SAPT_AMPS,"Ind30 uBS Amplitudes");

    inddisp30_amps();

    disp30_amps(PSIF_SAPT_AMPS,"tARBS Amplitudes",PSIF_SAPT_AA_DF_INTS,
      "AA RI Integrals","RR RI Integrals",PSIF_SAPT_BB_DF_INTS,
      "BB RI Integrals","SS RI Integrals",foccA_,noccA_,nvirA_,evalsA_,
      foccB_,noccB_,nvirB_,evalsB_,PSIF_SAPT_AMPS,"Disp30 uARBS Amplitudes");

  }
}

void SAPT2p3::Y3(int intfile, const char *AAlabel, const char *ARlabel,
  const char *RRlabel, int ampfile, const char *qAAlabel, const char *qRRlabel,
  const char *qbarAAlabel, const char *qbarRRlabel, const char *thetalabel,
  const char *tlabel, int foccA, int noccA, int nvirA, double *evals,
  int ampout, const char *Ylabel)
{
  int aoccA = noccA - foccA;

  double **yAR = block_matrix(aoccA,nvirA);

  Y2_1(yAR,intfile,ARlabel,RRlabel,ampfile,qRRlabel,foccA,noccA,nvirA);
  Y2_2(yAR,intfile,AAlabel,ARlabel,ampfile,qAAlabel,foccA,noccA,nvirA);
  Y2_1(yAR,intfile,ARlabel,RRlabel,ampfile,qbarRRlabel,foccA,noccA,nvirA);
  Y2_2(yAR,intfile,AAlabel,ARlabel,ampfile,qbarAAlabel,foccA,noccA,nvirA);
  Y2_3(yAR,intfile,AAlabel,RRlabel,ampfile,thetalabel,foccA,noccA,nvirA);
  Y3_1(yAR,intfile,AAlabel,ARlabel,ampfile,tlabel,foccA,noccA,nvirA);
  Y3_2(yAR,intfile,ARlabel,RRlabel,ampfile,tlabel,foccA,noccA,nvirA);
  Y3_3(yAR,intfile,AAlabel,ARlabel,RRlabel,ampfile,tlabel,foccA,noccA,nvirA);
  Y3_4(yAR,intfile,AAlabel,ARlabel,RRlabel,ampfile,tlabel,foccA,noccA,nvirA);

  psio_->write_entry(ampout,Ylabel,(char *) yAR[0],
    sizeof(double)*aoccA*nvirA);

  free_block(yAR);
}

void SAPT2p3::Y3_1(double **yAR, int intfile, const char *AAlabel,
  const char *ARlabel, int ampfile, const char *tlabel, int foccA,
  int noccA, int nvirA)
{
  int aoccA = noccA - foccA;

  double *tARAR = init_array((long int) aoccA*nvirA*aoccA*nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) tARAR,
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);
  ijkl_to_ikjl(tARAR,aoccA,nvirA,aoccA,nvirA);

  double *yAAAA = init_array((long int) aoccA*aoccA*aoccA*aoccA);

  C_DGEMM('N','T',aoccA*aoccA,aoccA*aoccA,nvirA*nvirA,1.0,tARAR,nvirA*nvirA,
    tARAR,nvirA*nvirA,0.0,yAAAA,aoccA*aoccA);

  free(tARAR);

  double **B_p_AA = get_DF_ints(intfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);

  double **AAAR = block_matrix(aoccA*aoccA,aoccA*nvirA);

  C_DGEMM('N','T',aoccA*aoccA,aoccA*nvirA,ndf_+3,1.0,B_p_AA[0],ndf_+3,
    B_p_AR[0],ndf_+3,0.0,AAAR[0],aoccA*nvirA);

  free_block(B_p_AA);
  free_block(B_p_AR);

  double **gAAAR = block_matrix(aoccA*aoccA,aoccA*nvirA);

  for(int a=0; a<aoccA; a++) {
    for(int a1=0; a1<aoccA; a1++) {
      for(int a2=0; a2<aoccA; a2++) {
        for(int r=0; r<nvirA; r++) {
          int ar = a*nvirA+r;
          int a1r = a1*nvirA+r;
          int a2r = a2*nvirA+r;
          int aa1 = a*aoccA+a1;
          int a1a2 = a1*aoccA+a2;
          int a2a1 = a2*aoccA+a1;
          gAAAR[a1a2][ar] = 2.0*AAAR[aa1][a2r] - AAAR[a2a1][ar];
  }}}}

  C_DGEMM('N','N',aoccA,nvirA,aoccA*aoccA*aoccA,1.0,yAAAA,aoccA*aoccA*aoccA,
    gAAAR[0],nvirA,1.0,yAR[0],nvirA);

  free(yAAAA);
  free_block(AAAR);
  free_block(gAAAR);
}

void SAPT2p3::Y3_2(double **yAR, int intfile, const char *ARlabel,
  const char *RRlabel, int ampfile, const char *tlabel, int foccA,
  int noccA, int nvirA)
{
  int aoccA = noccA - foccA;
  int virtri = nvirA*(nvirA+1)/2;

  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);
  double **AAAR = block_matrix(aoccA,aoccA*aoccA*nvirA);
  double **RRR = block_matrix(virtri,nvirA);
  double **xRRR = block_matrix(nvirA,nvirA*nvirA);
  double **X_RR = block_matrix(nvirA,nvirA);

  double *tARAR = init_array((long int) aoccA*nvirA*aoccA*nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) tARAR,
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);
  ijkl_to_ikjl(tARAR,aoccA,nvirA,aoccA,nvirA);

  double **B_p_RR = block_matrix(virtri,ndf_+3);

  psio_address next_DF_RR = PSIO_ZERO;

  for (int r1=0,r1r2=0; r1 < nvirA; r1++) {
  for (int r2=0; r2 <= r1; r2++,r1r2++) {
    next_DF_RR = psio_get_address(PSIO_ZERO,
      sizeof(double)*(r1*nvirA+r2)*(ndf_+3));
    psio_->read(intfile,RRlabel,(char *) &(B_p_RR[r1r2][0]),
      sizeof(double)*(ndf_+3),next_DF_RR,&next_DF_RR);
  }}

  for(int a=0; a<aoccA; a++) {
    C_DGEMM('N','T',virtri,nvirA,ndf_+3,1.0,B_p_RR[0],ndf_+3,B_p_AR[a*nvirA],
      ndf_+3,0.0,&(RRR[0][0]),nvirA);

    for(int r=0; r<nvirA; r++) {
      for(int r1=0,r1r2=0; r1<nvirA; r1++) {
        for(int r2=0; r2<nvirA; r2++,r1r2++) {
          xRRR[r][r1r2] = RRR[INDEX(r,r1)][r2] - 2.0*RRR[INDEX(r,r2)][r1];
    }}}

    C_DGEMM('N','T',aoccA*aoccA,nvirA,nvirA*nvirA,1.0,&(tARAR[0]),nvirA*nvirA,
      &(xRRR[0][0]),nvirA*nvirA,0.0,&(AAAR[a][0]),nvirA);
  }

  for(int a1=0,a1a2=0; a1<aoccA; a1++) {
    for(int a2=0; a2<aoccA; a2++,a1a2++) {
      C_DCOPY(nvirA*nvirA,&(tARAR[(long int) a1a2*nvirA*nvirA]),1,
        &(X_RR[0][0]),1);
      for(int r=0; r<nvirA; r++) {
        C_DCOPY(nvirA,&(X_RR[0][r]),nvirA,
          &(tARAR[(long int) a1a2*nvirA*nvirA+r*nvirA]),1);
      }
  }}

  C_DGEMM('N','N',aoccA,nvirA,aoccA*aoccA*nvirA,1.0,&(AAAR[0][0]),
    aoccA*aoccA*nvirA,&(tARAR[0]),nvirA,1.0,&(yAR[0][0]),nvirA);

  free_block(B_p_AR);
  free_block(B_p_RR);
  free_block(AAAR);
  free_block(RRR);
  free_block(xRRR);
  free_block(X_RR);
  free(tARAR);
}

void SAPT2p3::Y3_3(double **yAR, int intfile, const char *AAlabel,
  const char *ARlabel, const char *RRlabel, int ampfile, const char *tlabel,
  int foccA, int noccA, int nvirA)
{
  int aoccA = noccA - foccA;

  double *tARAR = init_array((long int) aoccA*nvirA*aoccA*nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) tARAR,
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);

  double *thetaARAR = init_array((long int) aoccA*nvirA*aoccA*nvirA);
  C_DCOPY((long int) aoccA*nvirA*aoccA*nvirA,tARAR,1,thetaARAR,1);
  antisym(thetaARAR,aoccA,nvirA);

  double *yARAR = init_array((long int) aoccA*nvirA*aoccA*nvirA);

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,aoccA*nvirA,1.0,thetaARAR,
    aoccA*nvirA,tARAR,aoccA*nvirA,0.0,yARAR,aoccA*nvirA);

  C_DCOPY((long int) aoccA*nvirA*aoccA*nvirA,tARAR,1,thetaARAR,1);
  OVOpVp_to_OVpOpV(tARAR,aoccA,nvirA);

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,aoccA*nvirA,-1.0,thetaARAR,
    aoccA*nvirA,tARAR,aoccA*nvirA,1.0,yARAR,aoccA*nvirA);

  free(tARAR);
  free(thetaARAR);

  double **B_p_AA = get_DF_ints(intfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);
  double **B_p_RR = get_DF_ints(intfile,RRlabel,0,nvirA,0,nvirA);

  double **Y_p_AR = block_matrix(aoccA*nvirA,ndf_+3);

  C_DGEMM('N','N',aoccA*nvirA,ndf_+3,aoccA*nvirA,1.0,yARAR,aoccA*nvirA,
    B_p_AR[0],ndf_+3,0.0,Y_p_AR[0],ndf_+3);

  C_DGEMM('N','T',aoccA,nvirA,nvirA*(ndf_+3),2.0,Y_p_AR[0],nvirA*(ndf_+3),
    B_p_RR[0],nvirA*(ndf_+3),1.0,yAR[0],nvirA);

  for (int a=0; a<aoccA; a++) {
    C_DGEMM('N','T',aoccA,nvirA,ndf_+3,-2.0,B_p_AA[a*aoccA],ndf_+3,
      Y_p_AR[a*nvirA],ndf_+3,1.0,yAR[0],nvirA);
  }

  free_block(Y_p_AR);

  double **Y_p_AA = block_matrix(aoccA*aoccA,ndf_+3);

  ijkl_to_ikjl(yARAR,aoccA,nvirA,aoccA,nvirA);

  C_DGEMM('N','N',aoccA*aoccA,ndf_+3,nvirA*nvirA,1.0,yARAR,nvirA*nvirA,
    B_p_RR[0],ndf_+3,0.0,Y_p_AA[0],ndf_+3);

  for (int a=0; a<aoccA; a++) {
    C_DGEMM('N','T',aoccA,nvirA,ndf_+3,-1.0,Y_p_AA[a*aoccA],ndf_+3,
      B_p_AR[a*nvirA],ndf_+3,1.0,yAR[0],nvirA);
  }

  free_block(Y_p_AA);

  C_DGEMM('T','N',nvirA*nvirA,ndf_+3,aoccA*aoccA,1.0,yARAR,nvirA*nvirA,
    B_p_AA[0],ndf_+3,0.0,B_p_RR[0],ndf_+3);

  C_DGEMM('N','T',aoccA,nvirA,nvirA*(ndf_+3),1.0,B_p_AR[0],nvirA*(ndf_+3),
    B_p_RR[0],nvirA*(ndf_+3),1.0,yAR[0],nvirA);

  free(yARAR);
  free_block(B_p_AA);
  free_block(B_p_AR);
  free_block(B_p_RR);
}

void SAPT2p3::Y3_4(double **yAR, int intfile, const char *AAlabel,
  const char *ARlabel, const char *RRlabel, int ampfile, const char *tlabel,
  int foccA, int noccA, int nvirA)
{
  int aoccA = noccA - foccA;

  double *tARAR = init_array((long int) aoccA*nvirA*aoccA*nvirA);
  psio_->read_entry(ampfile,tlabel,(char *) tARAR,
    sizeof(double)*aoccA*nvirA*aoccA*nvirA);
  OVOpVp_to_OVpOpV(tARAR,aoccA,nvirA);

  double *yARAR = init_array((long int) aoccA*nvirA*aoccA*nvirA);

  C_DGEMM('N','T',aoccA*nvirA,aoccA*nvirA,aoccA*nvirA,1.0,tARAR,
    aoccA*nvirA,tARAR,aoccA*nvirA,0.0,yARAR,aoccA*nvirA);

  free(tARAR);

  double **B_p_AA = get_DF_ints(intfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_AR = get_DF_ints(intfile,ARlabel,foccA,noccA,0,nvirA);
  double **B_p_RR = get_DF_ints(intfile,RRlabel,0,nvirA,0,nvirA);

  double **Y_p_AR = block_matrix(aoccA*nvirA,ndf_+3);

  C_DGEMM('N','N',aoccA*nvirA,ndf_+3,aoccA*nvirA,1.0,yARAR,aoccA*nvirA,
    B_p_AR[0],ndf_+3,0.0,Y_p_AR[0],ndf_+3);

  C_DGEMM('N','T',aoccA,nvirA,nvirA*(ndf_+3),1.0,Y_p_AR[0],nvirA*(ndf_+3),
    B_p_RR[0],nvirA*(ndf_+3),1.0,yAR[0],nvirA);

  for (int a=0; a<aoccA; a++) {
    C_DGEMM('N','T',aoccA,nvirA,ndf_+3,-1.0,B_p_AA[a*aoccA],ndf_+3,
      Y_p_AR[a*nvirA],ndf_+3,1.0,yAR[0],nvirA);
  }

  free_block(Y_p_AR);

  double **Y_p_AA = block_matrix(aoccA*aoccA,ndf_+3);

  ijkl_to_ikjl(yARAR,aoccA,nvirA,aoccA,nvirA);

  C_DGEMM('N','N',aoccA*aoccA,ndf_+3,nvirA*nvirA,1.0,yARAR,nvirA*nvirA,
    B_p_RR[0],ndf_+3,0.0,Y_p_AA[0],ndf_+3);

  for (int a=0; a<aoccA; a++) {
    C_DGEMM('N','T',aoccA,nvirA,ndf_+3,-2.0,Y_p_AA[a*aoccA],ndf_+3,
      B_p_AR[a*nvirA],ndf_+3,1.0,yAR[0],nvirA);
  }

  free_block(Y_p_AA);

  C_DGEMM('T','N',nvirA*nvirA,ndf_+3,aoccA*aoccA,1.0,yARAR,nvirA*nvirA,
    B_p_AA[0],ndf_+3,0.0,B_p_RR[0],ndf_+3);

  C_DGEMM('N','T',aoccA,nvirA,nvirA*(ndf_+3),2.0,B_p_AR[0],nvirA*(ndf_+3),
    B_p_RR[0],nvirA*(ndf_+3),1.0,yAR[0],nvirA);

  free(yARAR);
  free_block(B_p_AA);
  free_block(B_p_AR);
  free_block(B_p_RR);
}

void SAPT2p3::ind30_amps(int AAfile, const char *ARlabel, int BBfile,
  const char *BSlabel, double **wBAA, double **wBAR, double **wBRR,
  double **wABS, int noccA, int nvirA, double *evalsA, int noccB, int nvirB,
  double *evalsB, int ampout, const char *amplabel)
{
  double **sAR = block_matrix(noccA,nvirA);
  double **sBS = block_matrix(noccB,nvirB);

  for (int a=0; a<noccA; a++) {
    for (int r=0; r<nvirA; r++) {
      sAR[a][r] = wBAR[a][r] / (evalsA[a] - evalsA[r+noccA]);
  }}

  for (int b=0; b<noccB; b++) {
    for (int s=0; s<nvirB; s++) {
      sBS[b][s] = wABS[b][s] / (evalsB[b] - evalsB[s+noccB]);
  }}

  double **uAR = block_matrix(noccA,nvirA);

  C_DGEMM('N','T',noccA,nvirA,nvirA,1.0,sAR[0],nvirA,wBRR[0],nvirA,0.0,uAR[0],
    nvirA);

  C_DGEMM('N','N',noccA,nvirA,noccA,-1.0,wBAA[0],noccA,sAR[0],nvirA,1.0,uAR[0],
    nvirA);

  double **B_p_AR = get_DF_ints(AAfile,ARlabel,0,noccA,0,nvirA);
  double **B_p_BS = get_DF_ints(BBfile,BSlabel,0,noccB,0,nvirB);

  double *X = init_array(ndf_+3);

  C_DGEMV('t',noccB*nvirB,ndf_+3,1.0,B_p_BS[0],ndf_+3,sBS[0],1,0.0,X,1);
  C_DGEMV('n',noccA*nvirA,ndf_+3,2.0,B_p_AR[0],ndf_+3,X,1,1.0,uAR[0],1);

  free(X);

  double **tARBS = block_matrix(noccA*nvirA,noccB*nvirB);

  C_DGEMM('N','T',noccA*nvirA,noccB*nvirB,ndf_+3,1.0,B_p_AR[0],ndf_+3,
    B_p_BS[0],ndf_+3,0.0,tARBS[0],noccB*nvirB);

  free_block(B_p_AR);
  free_block(B_p_BS);

  for (int a=0, ar=0; a<noccA; a++) {
    for (int r=0; r<nvirA; r++, ar++) {
      for (int b=0, bs=0; b<noccB; b++) {
        for (int s=0; s<nvirB; s++, bs++) {
          tARBS[ar][bs] /= evalsA[a] + evalsB[b] - evalsA[r+noccA]
            - evalsB[s+noccB];
  }}}}

 C_DGEMV('n',noccA*nvirA,noccB*nvirB,2.0,tARBS[0],noccB*nvirB,wABS[0],1,1.0,
   uAR[0],1);

  free_block(tARBS);
  free_block(sAR);
  free_block(sBS);

  for (int a=0; a<noccA; a++) {
    for (int r=0; r<nvirA; r++) {
      uAR[a][r] /= evalsA[a] - evalsA[r+noccA];
  }}

  psio_->write_entry(ampout,amplabel,(char *) uAR[0],
    sizeof(double)*noccA*nvirA);

  free_block(uAR);
}

void SAPT2p3::inddisp30_amps()
{
  inddisp30_ov(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","RR RI Integrals",
    PSIF_SAPT_AMPS,"T AR Intermediates",foccA_,noccA_,nvirA_,evalsA_,
    PSIF_SAPT_AMPS,"IndDisp30 uAR Amplitudes");
  inddisp30_ov(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","SS RI Integrals",
    PSIF_SAPT_AMPS,"T BS Intermediates",foccB_,noccB_,nvirB_,evalsB_,
    PSIF_SAPT_AMPS,"IndDisp30 uBS Amplitudes");
  inddisp30_ovov();
}

void SAPT2p3::inddisp30_ov(int AAfile, const char *AAlabel,
  const char *RRlabel, int ampfile, const char *Tlabel, int foccA, int noccA,
  int nvirA, double *evalsA, int ampout, const char *amplabel)
{
  int aoccA = noccA - foccA;

  double **B_p_AA = get_DF_ints(AAfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_RR = get_DF_ints(AAfile,RRlabel,0,nvirA,0,nvirA);

  double **T_p_AR = block_matrix(aoccA*nvirA,ndf_+3);
  psio_->read_entry(ampfile,Tlabel,(char *) T_p_AR[0],
    sizeof(double)*aoccA*nvirA*(ndf_+3));

  double **uAR = block_matrix(aoccA,nvirA);

  C_DGEMM('N','T',aoccA,nvirA,nvirA*(ndf_+3),2.0,T_p_AR[0],nvirA*(ndf_+3),
    B_p_RR[0],nvirA*(ndf_+3),0.0,uAR[0],nvirA);

  for (int a=0; a<aoccA; a++) {
    C_DGEMM('N','T',aoccA,nvirA,ndf_+3,-2.0,B_p_AA[a*aoccA],ndf_+3,
      T_p_AR[a*nvirA],ndf_+3,1.0,uAR[0],nvirA);
  }

  free_block(B_p_AA);
  free_block(B_p_RR);
  free_block(T_p_AR);

  for (int a=0; a<aoccA; a++) {
    for (int r=0; r<nvirA; r++) {
      uAR[a][r] /= evalsA[a+foccA] - evalsA[r+noccA];
  }}

  psio_->write_entry(ampout,amplabel,(char *) uAR[0],
    sizeof(double)*aoccA*nvirA);

  free_block(uAR);
}

void SAPT2p3::inddisp30_ovov()
{
  double **sAR = block_matrix(aoccA_,nvirA_);
  double **sBS = block_matrix(aoccB_,nvirB_);

  for (int a=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++) {
      sAR[a][r] = wBAR_[a+foccA_][r] / (evalsA_[a+foccA_] - evalsA_[r+noccA_]);
  }}

  for (int b=0; b<aoccB_; b++) {
    for (int s=0; s<nvirB_; s++) {
      sBS[b][s] = wABS_[b+foccB_][s] / (evalsB_[b+foccB_] - evalsB_[s+noccB_]);
  }}

  double **uARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);

  double **B_p_AA = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
    foccA_,noccA_,foccA_,noccA_);
  double **B_p_RR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"RR RI Integrals",
    0,nvirA_,0,nvirA_);

  double **X_p_AR = block_matrix(aoccA_*nvirA_,ndf_+3);
  double **B_p_BS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
    foccB_,noccB_,0,nvirB_);

  C_DGEMM('N','N',aoccA_,nvirA_*(ndf_+3),nvirA_,1.0,sAR[0],nvirA_,
    B_p_RR[0],nvirA_*(ndf_+3),0.0,X_p_AR[0],nvirA_*(ndf_+3));

  for(int a=0; a<aoccA_; a++) {
    C_DGEMM('T','N',nvirA_,ndf_+3,aoccA_,-1.0,sAR[0],nvirA_,
      B_p_AA[a*aoccA_],ndf_+3,1.0,X_p_AR[a*nvirA_],ndf_+3);
  }

  C_DGEMM('N','T',aoccA_*nvirA_,aoccB_*nvirB_,ndf_+3,1.0,X_p_AR[0],ndf_+3,
    B_p_BS[0],ndf_+3,0.0,uARBS[0],aoccB_*nvirB_);

  free_block(B_p_AA);
  free_block(B_p_RR);
  free_block(X_p_AR);
  free_block(B_p_BS);

  double **B_p_BB = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
    foccB_,noccB_,foccB_,noccB_);
  double **B_p_SS = get_DF_ints(PSIF_SAPT_BB_DF_INTS,"SS RI Integrals",
    0,nvirB_,0,nvirB_);

  double **X_p_BS = block_matrix(aoccB_*nvirB_,ndf_+3);
  double **B_p_AR = get_DF_ints(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    foccA_,noccA_,0,nvirA_);

  C_DGEMM('N','N',aoccB_,nvirB_*(ndf_+3),nvirB_,1.0,sBS[0],nvirB_,
    B_p_SS[0],nvirB_*(ndf_+3),0.0,X_p_BS[0],nvirB_*(ndf_+3));

  for(int b=0; b<aoccB_; b++) {
    C_DGEMM('T','N',nvirB_,ndf_+3,aoccB_,-1.0,sBS[0],nvirB_,
      B_p_BB[b*aoccB_],ndf_+3,1.0,X_p_BS[b*nvirB_],ndf_+3);
  }

  C_DGEMM('N','T',aoccA_*nvirA_,aoccB_*nvirB_,ndf_+3,1.0,B_p_AR[0],ndf_+3,
    X_p_BS[0],ndf_+3,1.0,uARBS[0],aoccB_*nvirB_);

  free_block(B_p_BB);
  free_block(B_p_SS);
  free_block(B_p_AR);
  free_block(X_p_BS);

  double **tARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
  psio_->read_entry(PSIF_SAPT_AMPS,"tARBS Amplitudes",(char *) tARBS[0],
    sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);

  C_DGEMM('N','N',aoccA_,nvirA_*aoccB_*nvirB_,aoccA_,-1.0,
    &(wBAA_[foccA_][foccA_]),noccA_,tARBS[0],nvirA_*aoccB_*nvirB_,1.0,
    uARBS[0],nvirA_*aoccB_*nvirB_);

  for (int a=0; a<aoccA_; a++) {
    C_DGEMM('N','N',nvirA_,aoccB_*nvirB_,nvirA_,1.0,wBRR_[0],nvirA_,
      tARBS[a*nvirA_],aoccB_*nvirB_,1.0,uARBS[a*nvirA_],aoccB_*nvirB_);
  }

  for (int a=0, ar=0; a<aoccA_; a++) {
    for (int r=0; r<nvirA_; r++, ar++) {
      C_DGEMM('N','N',aoccB_,nvirB_,aoccB_,-1.0,&(wABB_[foccB_][foccB_]),
        noccB_,tARBS[ar],nvirB_,1.0,uARBS[ar],nvirB_);
  }}

  C_DGEMM('N','N',aoccA_*nvirA_*aoccB_,nvirB_,nvirB_,1.0,tARBS[0],nvirB_,
    wASS_[0],nvirB_,1.0,uARBS[0],nvirB_);

  free_block(tARBS);
  free_block(sAR);
  free_block(sBS);

  for (int a=0, ar=0; a<aoccA_; a++){
    for (int r=0; r<nvirA_; r++, ar++){
      for (int b=0, bs=0; b<aoccB_; b++){
        for (int s=0; s<nvirB_; s++, bs++){
          uARBS[ar][bs] /= evalsA_[a+foccA_]+evalsB_[b+foccB_]
            -evalsA_[r+noccA_]-evalsB_[s+noccB_];
  }}}}

  psio_->write_entry(PSIF_SAPT_AMPS,"IndDisp30 uARBS Amplitudes",
    (char *) uARBS[0],sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);

  free_block(uARBS);
}

void SAPT2p3::disp30_amps(int ampfile, const char *amplabel, int AAintfile,
  const char *AAlabel, const char *RRlabel, int BBintfile,
  const char *BBlabel, const char *SSlabel, int foccA, int noccA, int nvirA,
  double *evalsA, int foccB, int noccB, int nvirB, double *evalsB,
  int ampout, const char *tlabel)
{
  int aoccA = noccA - foccA;
  int aoccB = noccB - foccB;

  double **tARBS = block_matrix(aoccA_*nvirA_,aoccB_*nvirB_);
  double **tABRS = block_matrix(aoccA*aoccB,nvirA*nvirB);

  psio_->read_entry(ampfile,amplabel,(char *) tARBS[0],
    sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);

    for (int a=0, ar=0; a < aoccA; a++) {
    for (int r=0; r < nvirA; r++, ar++) {
      for (int b=0, bs=0; b < aoccB; b++) {
      for (int s=0; s < nvirB; s++, bs++) {
        int ab = a*aoccB + b;
        int rs = r*nvirB + s;
        tABRS[ab][rs] = tARBS[ar][bs];
      }}
    }}

  free_block(tARBS);

  double **t2ABRS = block_matrix(aoccA*aoccB,nvirA*nvirB);

  double **B_p_RR = get_DF_ints(AAintfile,RRlabel,0,nvirA,0,nvirA);
  double **B_p_SS = get_DF_ints(BBintfile,SSlabel,0,nvirB,0,nvirB);

  double **X_RS = block_matrix(nvirA,nvirB*nvirB);

  for (int r=0; r < nvirA; r++) {
    C_DGEMM('N','T',nvirA,nvirB*nvirB,ndf_+3,1.0,&(B_p_RR[r*nvirA][0]),
      ndf_+3,&(B_p_SS[0][0]),ndf_+3,0.0,&(X_RS[0][0]),nvirB*nvirB);
    C_DGEMM('N','T',aoccA*aoccB,nvirA*nvirB,nvirB,1.0,&(tABRS[0][r*nvirB]),
      nvirA*nvirB,&(X_RS[0][0]),nvirB,1.0,&(t2ABRS[0][0]),nvirA*nvirB);
  }

  free_block(B_p_RR);
  free_block(B_p_SS);
  free_block(X_RS);

  double **B_p_AA = get_DF_ints(AAintfile,AAlabel,foccA,noccA,foccA,noccA);
  double **B_p_BB = get_DF_ints(BBintfile,BBlabel,foccB,noccB,foccB,noccB);

  double **ABAB = block_matrix(aoccA*aoccB,aoccA*aoccB);

  for (int a=0, ab=0; a < aoccA; a++) {
    for (int b=0; b < aoccB; b++,ab++) {
      C_DGEMM('N','T',aoccA,aoccB,ndf_+3,1.0,&(B_p_AA[a*aoccA][0]),
        ndf_+3,&(B_p_BB[b*aoccB][0]),ndf_+3,0.0,&(ABAB[ab][0]),aoccB);
  }}

  free_block(B_p_AA);
  free_block(B_p_BB);

  C_DGEMM('N','N',aoccA*aoccB,nvirA*nvirB,aoccA*aoccB,1.0,&(ABAB[0][0]),
    aoccA*aoccB,&(tABRS[0][0]),nvirA*nvirB,1.0,&(t2ABRS[0][0]),nvirA*nvirB);

  free_block(ABAB);

  double **tBRAS = block_matrix(aoccB*nvirA,aoccA*nvirB);

    for (int a=0, ab=0; a < aoccA; a++) {
    for (int b=0; b < aoccB; b++, ab++) {
      for (int r=0, rs=0; r < nvirA; r++) {
      for (int s=0; s < nvirB; s++, rs++) {
        int br = b*nvirA + r;
        int as = a*nvirB + s;
        tBRAS[br][as] = tABRS[ab][rs];
      }}
    }}

  free_block(tABRS);

  double **t2BRAS = block_matrix(aoccB*nvirA,aoccA*nvirB);

    for (int a=0, ab=0; a < aoccA; a++) {
    for (int b=0; b < aoccB; b++, ab++) {
      for (int r=0, rs=0; r < nvirA; r++) {
      for (int s=0; s < nvirB; s++, rs++) {
        int br = b*nvirA + r;
        int as = a*nvirB + s;
        t2BRAS[br][as] = t2ABRS[ab][rs];
      }}
    }}

  free_block(t2ABRS);

  B_p_BB = get_DF_ints(BBintfile,BBlabel,foccB,noccB,foccB,noccB);
  B_p_RR = get_DF_ints(AAintfile,RRlabel,0,nvirA,0,nvirA);

  double **BRBR = block_matrix(aoccB*nvirA,aoccB*nvirA);

  for (int b=0, br=0; b < aoccB; b++) {
    for (int r=0; r < nvirA; r++, br++) {
      C_DGEMM('N','T',aoccB,nvirA,ndf_+3,1.0,&(B_p_BB[b*aoccB][0]),
        ndf_+3,&(B_p_RR[r*nvirA][0]),ndf_+3,0.0,&(BRBR[br][0]),nvirA);
  }}

  free_block(B_p_BB);
  free_block(B_p_RR);

  C_DGEMM('N','N',aoccB*nvirA,aoccA*nvirB,aoccB*nvirA,-1.0,&(BRBR[0][0]),
    aoccB*nvirA,&(tBRAS[0][0]),aoccA*nvirB,1.0,&(t2BRAS[0][0]),aoccA*nvirB);

  free_block(BRBR);

  B_p_AA = get_DF_ints(AAintfile,AAlabel,foccA,noccA,foccA,noccA);
  B_p_SS = get_DF_ints(BBintfile,SSlabel,0,nvirB,0,nvirB);

  double **ASAS = block_matrix(aoccA*nvirB,aoccA*nvirB);

  for (int a=0, as=0; a < aoccA; a++) {
    for (int s=0; s < nvirB; s++, as++) {
      C_DGEMM('N','T',aoccA,nvirB,ndf_+3,1.0,&(B_p_AA[a*aoccA][0]),
        ndf_+3,&(B_p_SS[s*nvirB][0]),ndf_+3,0.0,&(ASAS[as][0]),nvirB);
  }}

  free_block(B_p_AA);
  free_block(B_p_SS);

  C_DGEMM('N','N',aoccB*nvirA,aoccA*nvirB,aoccA*nvirB,-1.0,&(tBRAS[0][0]),
    aoccA*nvirB,&(ASAS[0][0]),aoccA*nvirB,1.0,&(t2BRAS[0][0]),aoccA*nvirB);

  free_block(ASAS);
  free_block(tBRAS);

  tARBS = block_matrix(aoccA*nvirA,aoccB*nvirB);

  for (int a=0,ar=0; a < aoccA; a++) {
    for (int r=0; r < nvirA; r++,ar++) {
      for (int b=0,bs=0; b < aoccB; b++) {
        for (int s=0; s < nvirB; s++,bs++) {
          int br = b*nvirA + r;
          int as = a*nvirB + s;
          double denom = evalsA[a+foccA]+evalsB[b+foccB]-
                  evalsA[r+noccA]-evalsB[s+noccB];
          tARBS[ar][bs] = t2BRAS[br][as]/denom;
  }}}}

  free_block(t2BRAS);

  psio_->write_entry(ampout,tlabel,(char *) tARBS[0],
    sizeof(double)*aoccA_*nvirA_*aoccB_*nvirB_);

  free_block(tARBS);
}

}}
