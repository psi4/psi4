/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "sapt2p.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"


namespace psi { namespace sapt {

void SAPT2p::disp2ccd() //!
{
  int occA = noccA_ - foccA_;
  int occB = noccB_ - foccB_;

  if (print_) {
    outfile->Printf("\n==> CCD Dispersion <==\n\n");
  }

  std::shared_ptr<Matrix> mo2noA;
  std::shared_ptr<Matrix> mo2noB;
  if (nat_orbs_v4_) {
    double cutoff = options_.get_double("OCC_TOLERANCE");
    mo2noA = mo2no(PSIF_SAPT_AMPS,"pRR Density Matrix", nvirA_,cutoff);
    mo2noB = mo2no(PSIF_SAPT_AMPS,"pSS Density Matrix", nvirB_,cutoff);
    //if (print_) {
    //  outfile->Printf("Natural Orbital (vv|vv): Cutoff %11.3E\n", cutoff);
    //  outfile->Printf("Monomer A: %4d of %4d nvir used.\n", mo2noA->colspi()[0], mo2noA->rowspi()[0]);
    //  outfile->Printf("Monomer B: %4d of %4d nvir used.\n", mo2noB->colspi()[0], mo2noB->rowspi()[0]);
    //  outfile->Printf("\n");
    //}
  }

  if (print_) {
    outfile->Printf("Beginning Monomer A CCD\n\n");

  }

  // Calculate monomer A CCD amplitudes
timer_on("CCD Prep           ");
  ccd_prep("T ARAR Amplitudes","Theta ARAR Amplitudes","G ARAR Integrals",
    "G ARRA Integrals","AAAA Integrals","ARAR Integrals","AARR Integrals",
    "RRRR+ Integrals","RRRR- Integrals",PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
    "AR RI Integrals","RR RI Integrals",evalsA_,noccA_,
    nvirA_,foccA_,mo2noA,
    "T2 ARAR Amplitudes");
timer_off("CCD Prep           ");
  ccd_iterate("T ARAR Amplitudes","T ARAR Error","Theta ARAR Amplitudes",
    "G ARAR Integrals","G ARRA Integrals","AAAA Integrals","ARAR Integrals",
    "AARR Integrals","RRRR+ Integrals","RRRR- Integrals",evalsA_,
    noccA_,nvirA_,foccA_,mo2noA);

  if (print_) {
    outfile->Printf("Beginning Monomer B CCD\n\n");

  }
  // Calculate monomer B CCD amplitudes
timer_on("CCD Prep           ");
  ccd_prep("T BSBS Amplitudes","Theta BSBS Amplitudes","G BSBS Integrals",
    "G BSSB Integrals","BBBB Integrals","BSBS Integrals","BBSS Integrals",
    "SSSS+ Integrals","SSSS- Integrals",PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
    "BS RI Integrals","SS RI Integrals",evalsB_,noccB_,
    nvirB_,foccB_,mo2noB,
    "T2 BSBS Amplitudes");
timer_off("CCD Prep           ");
  ccd_iterate("T BSBS Amplitudes","T BSBS Error","Theta BSBS Amplitudes",
    "G BSBS Integrals","G BSSB Integrals","BBBB Integrals","BSBS Integrals",
    "BBSS Integrals","SSSS+ Integrals","SSSS- Integrals",evalsB_,
    noccB_,nvirB_,foccB_,mo2noB);

  if (print_) {
    outfile->Printf("Beginning CCD Dispersion Amplitude Formation\n\n");

  }
  // Calculate dispersion CCD amplitudes
timer_on("CCD Disp Prep      ");
  r_ccd_prep("T ARBS Amplitudes","ARBS Integrals","T ARBS (ARBS)",
    "T ARAR Amplitudes","Theta ARAR Amplitudes","T BSBS Amplitudes",
    "Theta BSBS Amplitudes","G ARAR Integrals","G BSBS Integrals",
    "Theta x G ARAR","T x G ARAR","T x G AA","T x G RR","Theta x G BSBS",
    "T x G BSBS","T x G BB","T x G SS",PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",evalsA_,evalsB_,
    noccA_,nvirA_,foccA_,noccB_,
    nvirB_,foccB_);
timer_off("CCD Disp Prep      ");
  double disp2_ccd = r_ccd_iterate("T ARBS Amplitudes","T ARBS Error",
    "T ARBS (ARBS)","G ARRA Integrals","G BSSB Integrals","Theta x G ARAR",
    "T x G AA","T x G RR","Theta x G BSBS","T x G BB","T x G SS",
    "ARBS Integrals",evalsA_,evalsB_,noccA_,
    nvirA_,foccA_,noccB_,nvirB_,
    foccB_);

  double **ARBS = block_matrix(occA*nvirA_,occB*nvirB_); //!
  double **BSAR = block_matrix(occB*nvirB_,occA*nvirA_); //!

  psio_->read_entry(PSIF_SAPT_CCD,"ARBS Integrals",(char *) &(ARBS[0][0]),
    occA*nvirA_*occB*nvirB_*(ULI) sizeof(double));

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<nvirA_; r1++,a1r1++) {
    for(int b1=0,b1s1=0; b1<occB; b1++) {
    for(int s1=0; s1<nvirB_; s1++,b1s1++) {
      BSAR[b1s1][a1r1] = ARBS[a1r1][b1s1];
  }}}}

  psio_->write_entry(PSIF_SAPT_CCD,"BSAR Integrals",(char *) &(BSAR[0][0]),
    occA*nvirA_*occB*nvirB_*(ULI) sizeof(double));

  psio_->read_entry(PSIF_SAPT_CCD,"T ARBS Amplitudes",(char *) &(ARBS[0][0]),
    occA*nvirA_*occB*nvirB_*(ULI) sizeof(double));

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<nvirA_; r1++,a1r1++) {
    for(int b1=0,b1s1=0; b1<occB; b1++) {
    for(int s1=0; s1<nvirB_; s1++,b1s1++) {
      BSAR[b1s1][a1r1] = ARBS[a1r1][b1s1];
  }}}}

  psio_->write_entry(PSIF_SAPT_CCD,"T BSAR Amplitudes",(char *) &(BSAR[0][0]),
    occA*nvirA_*occB*nvirB_*(ULI) sizeof(double));

  free_block(ARBS); //!
  free_block(BSAR); //!

  if (print_) {
    outfile->Printf("Beginning Intramonomer A CCD Dispersion\n\n");

  }
  // Calculate intramonomer dispersion CCD amplitudes for monomer A
timer_on("CCD Intra Prep     ");
  s_ccd_prep("S ARAR Amplitudes","S ARAR (tARBS)","T ARAR Amplitudes",
    "Theta ARAR Amplitudes","T ARBS Amplitudes","G BSBS Integrals",
    "ARBS Integrals",evalsA_,noccA_,nvirA_,
    foccA_,noccB_,nvirB_,foccB_);
timer_off("CCD Intra Prep     ");
  disp2_ccd += s_ccd_iterate("S ARAR Amplitudes","S ARAR Error",
    "S ARAR (tARBS)","T ARAR Amplitudes","G ARAR Integrals","G ARRA Integrals",
    "AAAA Integrals","ARAR Integrals","AARR Integrals","RRRR+ Integrals",
    "RRRR- Integrals","Theta x G ARAR","T x G ARAR","T x G AA","T x G RR",
    evalsA_,noccA_,nvirA_,foccA_,
    mo2noA);

  if (print_) {
    outfile->Printf("Beginning Intramonomer B CCD Dispersion\n\n");

  }
  // Calculate intramonomer dispersion CCD amplitudes for monomer B
timer_on("CCD Intra Prep     ");
  s_ccd_prep("S BSBS Amplitudes","S BSBS (tBSAR)","T BSBS Amplitudes",
    "Theta BSBS Amplitudes","T BSAR Amplitudes","G ARAR Integrals",
    "BSAR Integrals",evalsB_,noccB_,nvirB_,
    foccB_,noccA_,nvirA_,foccA_);
timer_off("CCD Intra Prep     ");
  disp2_ccd += s_ccd_iterate("S BSBS Amplitudes","S BSBS Error",
    "S BSBS (tBSAR)","T BSBS Amplitudes","G BSBS Integrals","G BSSB Integrals",
    "BBBB Integrals","BSBS Integrals","BBSS Integrals","SSSS+ Integrals",
    "SSSS- Integrals","Theta x G BSBS","T x G BSBS","T x G BB","T x G SS",
    evalsB_,noccB_,nvirB_,foccB_,
    mo2noB);

  e_disp2d_ccd_ = disp2_ccd;
  outfile->Printf("    Disp2 (CCD)         = %18.12lf [Eh]\n",e_disp2d_ccd_);

  // => (S) <= //

  disp_s_prep("T AR Amplitudes","T(BS) AR","Theta ARAR Amplitudes",
    "T ARBS Amplitudes",PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
    "AR RI Integrals","RR RI Integrals",PSIF_SAPT_BB_DF_INTS,
    "BS RI Integrals",evalsA_,noccA_,nvirA_,
    foccA_,noccB_,nvirB_,foccB_);
  double d220s = disp220s(PSIF_SAPT_CCD,"T AR Amplitudes","T(BS) AR",
    PSIF_SAPT_AA_DF_INTS,"AA RI Integrals","RR RI Integrals",
    foccA_,noccA_,nvirA_);

  disp_s_prep("T BS Amplitudes","T(AR) BS","Theta BSBS Amplitudes",
    "T BSAR Amplitudes",PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
    "BS RI Integrals","SS RI Integrals",PSIF_SAPT_AA_DF_INTS,
    "AR RI Integrals",evalsB_,noccB_,nvirB_,
    foccB_,noccA_,nvirA_,foccA_);
  double d202s = disp220s(PSIF_SAPT_CCD,"T BS Amplitudes","T(AR) BS",
    PSIF_SAPT_BB_DF_INTS,"BB RI Integrals","SS RI Integrals",
    foccB_,noccB_,nvirB_);

  e_disp22s_ccd_ = d220s + d202s;

  if (print_) {
    outfile->Printf("\n    Disp220 (S)         = %18.12lf mEh\n",d220s*1000.0);
    outfile->Printf("    Disp202 (S)         = %18.12lf mEh\n",d202s*1000.0);
    outfile->Printf("    Disp22 (S)          = %18.12lf mEh\n\n",e_disp22s_ccd_*1000.0);

  }
}
void SAPT2p::r_ccd_prep(const char *TARBS, const char *ARBS, const char *CA_RBS, const char *TARAR, //!
  const char *ThetaARAR, const char *TBSBS, const char *ThetaBSBS, const char *GARAR, const char *GBSBS,
  const char *XARAR, const char *YARAR, const char *XAA, const char *XRR, const char *XBSBS, const char *YBSBS,
  const char *XBB, const char *XSS, int AAfile, const char *ARints, int BBfile, const char *BSints,
  double *evalsA_, double *evalsB_, int noccA_, int virA, int foccA_, int noccB_,
  int virB, int foccB_)
{
  int occA = noccA_ - foccA_;
  int occB = noccB_ - foccB_;

  double **B_p_AR = get_DF_ints_nongimp(AAfile,ARints,foccA_,noccA_,0,virA);
  double **B_p_BS = get_DF_ints_nongimp(BBfile,BSints,foccB_,noccB_,0,virB);
  double **vARBS = block_matrix(occA*virA,occB*virB); //!

  C_DGEMM('N','T',occA*virA,occB*virB,ndf_,1.0,&(B_p_AR[0][0]),
    ndf_,&(B_p_BS[0][0]),ndf_,0.0,&(vARBS[0][0]),
    occB*virB);

  psio_->write_entry(PSIF_SAPT_CCD,ARBS,(char *) &(vARBS[0][0]),
    occA*virA*occB*virB*(ULI) sizeof(double));

  free_block(B_p_AR);
  free_block(B_p_BS);

  double **tARBS = block_matrix(occA*virA,occB*virB); //!

  C_DCOPY((size_t)occA*virA*occB*virB,&(vARBS[0][0]),1,&(tARBS[0][0]),1);

  double **thARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,ThetaARAR,(char *) &(thARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occB*virB,occA*virA,1.0,thARAR[0],occA*virA,
    vARBS[0],occB*virB,1.0,tARBS[0],occB*virB);

  double **thBSBS = block_matrix(occB*virB,occB*virB); //!

  psio_->read_entry(PSIF_SAPT_CCD,ThetaBSBS,(char *) &(thBSBS[0][0]),
    occB*virB*occB*virB*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occB*virB,occB*virB,1.0,vARBS[0],occB*virB,
    thBSBS[0],occB*virB,1.0,tARBS[0],occB*virB);

  double **xARBS = block_matrix(occA*virA,occB*virB); //!

  C_DGEMM('N','N',occA*virA,occB*virB,occA*virA,1.0,thARAR[0],occA*virA,
    vARBS[0],occB*virB,0.0,xARBS[0],occB*virB);

  C_DGEMM('N','N',occA*virA,occB*virB,occB*virB,1.0,xARBS[0],occB*virB,
    thBSBS[0],occB*virB,1.0,tARBS[0],occB*virB);

  free_block(xARBS); //!
  free_block(vARBS); //!

  psio_->write_entry(PSIF_SAPT_CCD,CA_RBS,(char *) &(tARBS[0][0]),
    occA*virA*occB*virB*(ULI) sizeof(double));

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int b1=0,b1s1=0; b1<occB; b1++) {
    for(int s1=0; s1<virB; s1++,b1s1++) {
      double denom = evalsA_[a1+foccA_]+evalsB_[b1+foccB_]-evalsA_[r1+noccA_]-
        evalsB_[s1+noccB_];
      tARBS[a1r1][b1s1] /= denom;
  }}}}

  psio_->write_entry(PSIF_SAPT_CCD,TARBS,(char *) &(tARBS[0][0]),
    occA*virA*occB*virB*(ULI) sizeof(double));

  free_block(tARBS); //!

  double **xARAR = block_matrix(occA*virA,occA*virA); //!
  double **gARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,GARAR,(char *) &(gARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,1.0,thARAR[0],occA*virA,
    gARAR[0],occA*virA,0.0,xARAR[0],occA*virA);

  psio_->write_entry(PSIF_SAPT_CCD,XARAR,(char *) &(xARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  free_block(xARAR);  //!
  free_block(thARAR); //!

  double **xBSBS = block_matrix(occB*virB,occB*virB); //!
  double **gBSBS = block_matrix(occB*virB,occB*virB); //!

  psio_->read_entry(PSIF_SAPT_CCD,GBSBS,(char *) &(gBSBS[0][0]),
    occB*virB*occB*virB*(ULI) sizeof(double));

  C_DGEMM('N','N',occB*virB,occB*virB,occB*virB,1.0,thBSBS[0],occB*virB,
    gBSBS[0],occB*virB,0.0,xBSBS[0],occB*virB);

  psio_->write_entry(PSIF_SAPT_CCD,XBSBS,(char *) &(xBSBS[0][0]),
    occB*virB*occB*virB*(ULI) sizeof(double));

  free_block(xBSBS);  //!
  free_block(thBSBS); //!

  double **yARAR = block_matrix(occA*virA,occA*virA); //!
  double **tARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,TARAR,(char *) &(tARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,1.0,tARAR[0],occA*virA,
    gARAR[0],occA*virA,0.0,yARAR[0],occA*virA);

  psio_->write_entry(PSIF_SAPT_CCD,YARAR,(char *) &(yARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  free_block(yARAR); //!

  double **xAA = block_matrix(occA,occA);
  double **xRR = block_matrix(virA,virA);

  C_DGEMM('N','T',occA,occA,virA*occA*virA,1.0,tARAR[0],virA*occA*virA,
    gARAR[0],virA*occA*virA,0.0,xAA[0],occA);

  C_DGEMM('T','N',virA,virA,occA*virA*occA,1.0,gARAR[0],virA,
    tARAR[0],virA,0.0,xRR[0],virA);

  psio_->write_entry(PSIF_SAPT_CCD,XAA,(char *) &(xAA[0][0]),
    occA*occA*(ULI) sizeof(double));
  psio_->write_entry(PSIF_SAPT_CCD,XRR,(char *) &(xRR[0][0]),
    virA*virA*(ULI) sizeof(double));

  free_block(tARAR); //!
  free_block(gARAR); //!
  free_block(xAA);
  free_block(xRR);

  double **yBSBS = block_matrix(occB*virB,occB*virB); //!
  double **tBSBS = block_matrix(occB*virB,occB*virB); //!

  psio_->read_entry(PSIF_SAPT_CCD,TBSBS,(char *) &(tBSBS[0][0]),
    occB*virB*occB*virB*(ULI) sizeof(double));

  C_DGEMM('N','N',occB*virB,occB*virB,occB*virB,1.0,tBSBS[0],occB*virB,
    gBSBS[0],occB*virB,0.0,yBSBS[0],occB*virB);

  psio_->write_entry(PSIF_SAPT_CCD,YBSBS,(char *) &(yBSBS[0][0]),
    occB*virB*occB*virB*(ULI) sizeof(double));

  free_block(yBSBS); //!

  double **xBB = block_matrix(occB,occB);
  double **xSS = block_matrix(virB,virB);

  C_DGEMM('N','T',occB,occB,virB*occB*virB,1.0,tBSBS[0],virB*occB*virB,
    gBSBS[0],virB*occB*virB,0.0,xBB[0],occB);

  C_DGEMM('T','N',virB,virB,occB*virB*occB,1.0,gBSBS[0],virB,
    tBSBS[0],virB,0.0,xSS[0],virB);

  psio_->write_entry(PSIF_SAPT_CCD,XBB,(char *) &(xBB[0][0]),
    occB*occB*(ULI) sizeof(double));
  psio_->write_entry(PSIF_SAPT_CCD,XSS,(char *) &(xSS[0][0]),
    virB*virB*(ULI) sizeof(double));

  free_block(tBSBS); //!
  free_block(gBSBS); //!
  free_block(xBB);
  free_block(xSS);
}

double SAPT2p::r_ccd_energy(const char *TARBS, const char *ARBS, int occA, int virA,  //!
  int occB, int virB)
{
  double **vARBS = block_matrix(occA*virA,occB*virB); //!
  double **tARBS = block_matrix(occA*virA,occB*virB); //!

  psio_->read_entry(PSIF_SAPT_CCD,ARBS,(char *) &(vARBS[0][0]),
    occA*virA*occB*virB*(ULI) sizeof(double));
  psio_->read_entry(PSIF_SAPT_CCD,TARBS,(char *) &(tARBS[0][0]),
    occA*virA*occB*virB*(ULI) sizeof(double));

  double energy = C_DDOT((size_t)occA*virA*occB*virB,&(vARBS[0][0]),1,
    &(tARBS[0][0]),1);

  free_block(vARBS); //!
  free_block(tARBS); //!

  return(energy);
}

double SAPT2p::r_ccd_iterate(const char *TARBS, const char *TARBSerr, const char *CA_RBS,  //!
  const char *GARRA, const char *GBSSB, const char *XARAR, const char *XAA, const char *XRR, const char *XBSBS,
  const char *XBB, const char *XSS, const char *ARBS, double *evalsA_, double *evalsB_, int noccA_,
  int virA, int foccA_, int noccB_, int virB, int foccB_)
{
  int occA = noccA_ - foccA_;
  int occB = noccB_ - foccB_;

  if (print_) {
    outfile->Printf("Iter      Energy [mEh]        dE [mEh]           RMS [mEh]\n");

  }

  int iter = 1;
  double E_old=0.0, E_new=0.0, RMS=0.0;

  SAPTDIIS diis(PSIF_SAPT_CCD,TARBS,TARBSerr,(size_t)occA*virA*occB*virB,
    max_ccd_vecs_,psio_);

  do {
    E_new = r_ccd_energy(TARBS,ARBS,occA,virA,occB,virB);
    outfile->Printf("%4d %16.8lf %17.9lf %17.9lf",iter,E_new*4000.0,
      (E_old-E_new)*4000.0,RMS*4000.0);


    if (iter > 1 && (4000.0*fabs(E_old-E_new) < ccd_e_conv_ &&
      4000.0*RMS < ccd_t_conv_)) {
      if (iter > min_ccd_vecs_) {
        outfile->Printf("  DIIS\n");
      }
      break;
    }
    E_old = E_new;

timer_on("CCD Disp Amps      ");
    RMS = r_ccd_amplitudes(TARBS,TARBSerr,CA_RBS,GARRA,GBSSB,XARAR,XAA,XRR,
      XBSBS,XBB,XSS,evalsA_,evalsB_,noccA_,virA,foccA_,noccB_,virB,foccB_);
timer_off("CCD Disp Amps      ");

    diis.store_vectors();
    if (iter > min_ccd_vecs_) {
      diis.get_new_vector();
      outfile->Printf("  DIIS\n");
    }
    else {
      outfile->Printf("\n");
    }

    iter++;
  }
  while(iter<ccd_maxiter_+1);

  outfile->Printf("\n");

  return(4.0*E_new);
}

double SAPT2p::r_ccd_amplitudes(const char *TARBS, const char *TARBSerr, const char *CA_RBS,  //!
  const char *GARRA, const char *GBSSB, const char *XARAR, const char *XAA, const char *XRR, const char *XBSBS,
  const char *XBB, const char *XSS, double *evalsA_, double *evalsB_, int noccA_, int virA,
  int foccA_, int noccB_, int virB, int foccB_)
{
  int occA = noccA_ - foccA_;
  int occB = noccB_ - foccB_;

  double **t2ARBS = block_matrix(occA*virA,occB*virB); //!

  psio_->read_entry(PSIF_SAPT_CCD,CA_RBS,(char *) &(t2ARBS[0][0]),
    occA*virA*occB*virB*(ULI) sizeof(double));

  double **tARBS = block_matrix(occA*virA,occB*virB); //!

  psio_->read_entry(PSIF_SAPT_CCD,TARBS,(char *) &(tARBS[0][0]),
    occA*virA*occB*virB*(ULI) sizeof(double));

  double **gARRA = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,GARRA,(char *) &(gARRA[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occB*virB,occA*virA,1.0,gARRA[0],occA*virA,
    tARBS[0],occB*virB,1.0,t2ARBS[0],occB*virB);

  free_block(gARRA); //!

  double **gBSSB = block_matrix(occB*virB,occB*virB); //!

  psio_->read_entry(PSIF_SAPT_CCD,GBSSB,(char *) &(gBSSB[0][0]),
    occB*virB*occB*virB*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occB*virB,occB*virB,1.0,tARBS[0],occB*virB,
    gBSSB[0],occB*virB,1.0,t2ARBS[0],occB*virB);

  free_block(gBSSB); //!

  double **xARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,XARAR,(char *) &(xARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  double **xAA = block_matrix(occA,occA);
  double **xRR = block_matrix(virA,virA);

  psio_->read_entry(PSIF_SAPT_CCD,XAA,(char *) &(xAA[0][0]),
    occA*occA*(ULI) sizeof(double));
  psio_->read_entry(PSIF_SAPT_CCD,XRR,(char *) &(xRR[0][0]),
    virA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occB*virB,occA*virA,1.0,xARAR[0],occA*virA,
    tARBS[0],occB*virB,1.0,t2ARBS[0],occB*virB);

  free_block(xARAR);  //!

  C_DGEMM('N','N',occA,virA*occB*virB,occA,-1.0,xAA[0],occA,
    tARBS[0],virA*occB*virB,1.0,t2ARBS[0],virA*occB*virB);

  for(int a1=0; a1<occA; a1++) {
    C_DGEMM('T','N',virA,occB*virB,virA,-1.0,xRR[0],virA,
      tARBS[a1*virA],occB*virB,1.0,t2ARBS[a1*virA],occB*virB);
  }

  free_block(xAA);
  free_block(xRR);

  double **xBSBS = block_matrix(occB*virB,occB*virB); //!

  psio_->read_entry(PSIF_SAPT_CCD,XBSBS,(char *) &(xBSBS[0][0]),
    occB*virB*occB*virB*(ULI) sizeof(double));

  double **xBB = block_matrix(occB,occB);
  double **xSS = block_matrix(virB,virB);

  psio_->read_entry(PSIF_SAPT_CCD,XBB,(char *) &(xBB[0][0]),
    occB*occB*(ULI) sizeof(double));
  psio_->read_entry(PSIF_SAPT_CCD,XSS,(char *) &(xSS[0][0]),
    virB*virB*(ULI) sizeof(double));

  C_DGEMM('N','T',occA*virA,occB*virB,occB*virB,1.0,tARBS[0],occB*virB,
    xBSBS[0],occB*virB,1.0,t2ARBS[0],occB*virB);

  free_block(xBSBS); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    C_DGEMM('N','N',occB,virB,occB,-1.0,xBB[0],occB,
      tARBS[a1r1],virB,1.0,t2ARBS[a1r1],virB);
  }}

  C_DGEMM('N','N',occA*virA*occB,virB,virB,-1.0,tARBS[0],virB,
    xSS[0],virB,1.0,t2ARBS[0],virB);

  free_block(xBB);
  free_block(xSS);

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int b1=0,b1s1=0; b1<occB; b1++) {
    for(int s1=0; s1<virB; s1++,b1s1++) {
      double denom = evalsA_[a1+foccA_]+evalsB_[b1+foccB_]-evalsA_[r1+noccA_]-
        evalsB_[s1+noccB_];
      t2ARBS[a1r1][b1s1] /= denom;
  }}}}

  psio_->write_entry(PSIF_SAPT_CCD,TARBS,(char *) &(t2ARBS[0][0]),
    occA*virA*occB*virB*(ULI) sizeof(double));

  C_DAXPY((size_t)occA*virA*occB*virB,-1.0,tARBS[0],1,t2ARBS[0],1);

  double RMS = C_DDOT((size_t)occA*virA*occB*virB,t2ARBS[0],1,t2ARBS[0],1);
  RMS /= (double) ((size_t)occA*virA*occB*virB);

  psio_->write_entry(PSIF_SAPT_CCD,TARBSerr,(char *) &(t2ARBS[0][0]),
    occA*virA*occB*virB*(ULI) sizeof(double));

  free_block(tARBS);  //!
  free_block(t2ARBS); //!

  return(sqrt(RMS));
}

void SAPT2p::s_ccd_prep(const char *SARAR, const char *CA_RAR, const char *TARAR, const char *ThetaARAR, //!
  const char *TARBS, const char *GBSBS, const char *ARBS, double *evalsA_, int noccA_, int virA,
  int foccA_, int noccB_, int virB, int foccB_)
{
  int occA = noccA_ - foccA_;
  int occB = noccB_ - foccB_;

  double **sARAR = block_matrix(occA*virA,occA*virA); //!
  double **tARBS = block_matrix(occA*virA,occB*virB); //!
  double **vARBS = block_matrix(occA*virA,occB*virB); //!

  psio_->read_entry(PSIF_SAPT_CCD,TARBS,(char *) &(tARBS[0][0]),
    occA*virA*occB*virB*(ULI) sizeof(double));

  psio_->read_entry(PSIF_SAPT_CCD,ARBS,(char *) &(vARBS[0][0]),
    occA*virA*occB*virB*(ULI) sizeof(double));

  C_DGEMM('N','T',occA*virA,occA*virA,occB*virB,2.0,tARBS[0],occB*virB,
    vARBS[0],occB*virB,0.0,sARAR[0],occA*virA);

  double **xARAR = block_matrix(occA*virA,occA*virA); //!

  C_DGEMM('N','T',occA*virA,occA*virA,occB*virB,1.0,tARBS[0],occB*virB,
    vARBS[0],occB*virB,0.0,xARAR[0],occA*virA);

  double **thARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,ThetaARAR,(char *) &(thARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,2.0,xARAR[0],occA*virA,
    thARAR[0],occA*virA,1.0,sARAR[0],occA*virA);

  free_block(xARAR); //!
  free_block(thARAR); //!

  double **xAA = block_matrix(occA,occA);
  double **xRR = block_matrix(virA,virA);

  for(int a1=0; a1<occA; a1++) {
    for(int a2=0; a2<occA; a2++) {
      xAA[a1][a2] = C_DDOT(virA*occB*virB,tARBS[a1*virA],1,vARBS[a2*virA],1);
  }}

  for(int a1=0; a1<occA; a1++) {
    for(int r1=0; r1<virA; r1++) {
      for(int r2=0; r2<virA; r2++) {
        int a1r1 = a1*virA + r1;
        int a1r2 = a1*virA + r2;
        xRR[r1][r2] += C_DDOT(occB*virB,tARBS[a1r1],1,vARBS[a1r2],1);
  }}}

  free_block(vARBS); //!

  double **tARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,TARAR,(char *) &(tARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA,virA*occA*virA,occA,-2.0,xAA[0],occA,
    tARAR[0],virA*occA*virA,1.0,sARAR[0],virA*occA*virA);

  C_DGEMM('N','T',occA*virA*occA,virA,virA,-2.0,tARAR[0],virA,
    xRR[0],virA,1.0,sARAR[0],virA);

  free_block(xAA);
  free_block(xRR);
  free_block(tARAR); //!

  double **gBSBS = block_matrix(occB*virB,occB*virB); //!
  double **xARBS = block_matrix(occA*virA,occB*virB); //!

  psio_->read_entry(PSIF_SAPT_CCD,GBSBS,(char *) &(gBSBS[0][0]),
    occB*virB*occB*virB*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occB*virB,occB*virB,1.0,tARBS[0],occB*virB,
    gBSBS[0],occB*virB,0.0,xARBS[0],occB*virB);

  C_DGEMM('N','T',occA*virA,occA*virA,occB*virB,1.0,xARBS[0],occB*virB,
    tARBS[0],occB*virB,1.0,sARAR[0],occA*virA);

  free_block(gBSBS); //!
  free_block(xARBS); //!
  free_block(tARBS); //!

  psio_->write_entry(PSIF_SAPT_CCD,CA_RAR,(char *) &(sARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  for(int a1r1=0; a1r1<occA*virA; a1r1++) {
    for(int a2r2=0; a2r2<a1r1; a2r2++) {
      double tval = sARAR[a1r1][a2r2] + sARAR[a2r2][a1r1];
      sARAR[a1r1][a2r2] = tval;
      sARAR[a2r2][a1r1] = tval;
  }}

  for(int a1r1=0; a1r1<occA*virA; a1r1++) {
    sARAR[a1r1][a1r1] *= 2.0;
  }

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      double denom = evalsA_[a1+foccA_]+evalsA_[a2+foccA_]-evalsA_[r1+noccA_]-
        evalsA_[r2+noccA_];
      sARAR[a1r1][a2r2] /= denom;
  }}}}

  psio_->write_entry(PSIF_SAPT_CCD,SARAR,(char *) &(sARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  free_block(sARAR); //!
}

double SAPT2p::s_ccd_iterate(const char *SARAR, const char *SARARerr, const char *CA_RAR, const char *TARAR, //!
  const char *GARAR, const char *GARRA, const char *AAAA, const char *ARAR, const char *AARR, const char *RRRRp,
  const char *RRRRm, const char *XARAR, const char *YARAR, const char *XAA, const char *XRR, double *evalsA_,
  int noccA_, int virA, int foccA_, std::shared_ptr<Matrix> mo2no)
{
  int occA = noccA_ - foccA_;

  if (print_) {
    outfile->Printf("Iter      Energy [mEh]        dE [mEh]           RMS [mEh]\n");

  }

  int iter = 1;
  double E_old=0.0, E_new=0.0, RMS=0.0;

  SAPTDIIS diis(PSIF_SAPT_CCD,SARAR,SARARerr,(size_t)occA*virA*occA*virA,
    max_ccd_vecs_,psio_);

  do {
    E_new = ccd_energy(SARAR,GARAR,occA,virA);
    outfile->Printf("%4d %16.8lf %17.9lf %17.9lf",iter,E_new*1000.0,
      (E_old-E_new)*1000.0,RMS*1000.0);


    if (iter > 1 && (1000.0*fabs(E_old-E_new) < ccd_e_conv_ &&
      1000.0*RMS < ccd_t_conv_)) {
      if (iter > min_ccd_vecs_) {
        outfile->Printf("  DIIS\n");
      }
      break;
    }
    E_old = E_new;

timer_on("CCD Intra Amps     ");
    RMS = s_ccd_amplitudes(SARAR,SARARerr,CA_RAR,TARAR,GARAR,GARRA,AAAA,ARAR,
      AARR,RRRRp,RRRRm,XARAR,YARAR,XAA,XRR,evalsA_,noccA_,virA,foccA_,mo2no);
timer_off("CCD Intra Amps     ");

    diis.store_vectors();
    if (iter > min_ccd_vecs_) {
      diis.get_new_vector();
      outfile->Printf("  DIIS\n");
    }
    else {
      outfile->Printf("\n");
    }

    iter++;
  }
  while(iter<ccd_maxiter_+1);

  outfile->Printf("\n");

  return(E_new);
}

double SAPT2p::s_ccd_amplitudes(const char *SARAR, const char *SARARerr, const char *CA_RAR, //!
  const char *TARAR, const char *GARAR, const char *GARRA, const char *AAAA, const char *ARAR, const char *AARR,
  const char *RRRRp, const char *RRRRm, const char *XARAR, const char *YARAR, const char *XAA, const char *XRR,
  double *evalsA_, int noccA_, int virA, int foccA_, std::shared_ptr<Matrix> mo2no)
{
  int occA = noccA_ - foccA_;

  double **s2ARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,CA_RAR,(char *) &(s2ARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  double **gARRA = block_matrix(occA*virA,occA*virA); //!
  double **sARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,GARRA,(char *) &(gARRA[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));
  psio_->read_entry(PSIF_SAPT_CCD,SARAR,(char *) &(sARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,1.0,gARRA[0],occA*virA,
    sARAR[0],occA*virA,1.0,s2ARAR[0],occA*virA);

  free_block(gARRA); //!

  double **s2AARR = vvvv_ccd(SARAR,RRRRp,RRRRm,occA,virA,mo2no); //!!
  double **sAARR = block_matrix(occA*occA,virA*virA); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1a2 = a1*occA + a2;
      int r1r2 = r1*virA + r2;
      sAARR[a1a2][r1r2] = sARAR[a1r1][a2r2];
  }}}}

  double **vAAAA = block_matrix(occA*occA,occA*occA);

  psio_->read_entry(PSIF_SAPT_CCD,AAAA,(char *) &(vAAAA[0][0]),
    occA*occA*occA*occA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*occA,virA*virA,occA*occA,0.5,vAAAA[0],occA*occA,
    sAARR[0],virA*virA,0.5,s2AARR[0],virA*virA);

  free_block(vAAAA);

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1a2 = a1*occA + a2;
      int r1r2 = r1*virA + r2;
      s2ARAR[a1r1][a2r2] += s2AARR[a1a2][r1r2];
  }}}}

  free_block(sAARR);  //!
  free_block(s2AARR); //!

  double **sARRA = block_matrix(occA*virA,occA*virA); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      sARRA[a1r2][a2r1] = sARAR[a1r1][a2r2];
  }}}}

  double **vARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,ARAR,(char *) &(vARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,-1.0,sARRA[0],occA*virA,
    vARAR[0],occA*virA,1.0,s2ARAR[0],occA*virA);

  double **vAARR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,AARR,(char *) &(vAARR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,-1.0,sARRA[0],occA*virA,
    vAARR[0],occA*virA,0.0,vARAR[0],occA*virA);

  free_block(vAARR); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      s2ARAR[a1r1][a2r2] += vARAR[a1r2][a2r1];
  }}}}

  free_block(vARAR); //!

  double **xARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,XARAR,(char *) &(xARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,1.0,xARAR[0],occA*virA,
    sARAR[0],occA*virA,1.0,s2ARAR[0],occA*virA);

  psio_->read_entry(PSIF_SAPT_CCD,YARAR,(char *) &(xARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,-1.0,xARAR[0],occA*virA,
    sARRA[0],occA*virA,1.0,s2ARAR[0],occA*virA);

  free_block(sARRA); //!
  free_block(xARAR); //!

  double **gARAR = block_matrix(occA*virA,occA*virA); //!
  double **tARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,GARAR,(char *) &(gARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  psio_->read_entry(PSIF_SAPT_CCD,TARAR,(char *) &(tARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  double **xAA = block_matrix(occA,occA);
  double **xRR = block_matrix(virA,virA);

  psio_->read_entry(PSIF_SAPT_CCD,XAA,(char *) &(xAA[0][0]),
    occA*occA*(ULI) sizeof(double));
  psio_->read_entry(PSIF_SAPT_CCD,XRR,(char *) &(xRR[0][0]),
    virA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA,virA*occA*virA,occA,-1.0,xAA[0],occA,
    sARAR[0],virA*occA*virA,1.0,s2ARAR[0],virA*occA*virA);

  C_DGEMM('N','N',occA*virA*occA,virA,virA,-1.0,sARAR[0],virA,
    xRR[0],virA,1.0,s2ARAR[0],virA);

  C_DGEMM('N','T',occA,occA,virA*occA*virA,1.0,sARAR[0],virA*occA*virA,
    gARAR[0],virA*occA*virA,0.0,xAA[0],occA);

  C_DGEMM('T','N',virA,virA,occA*virA*occA,1.0,gARAR[0],virA,
    sARAR[0],virA,0.0,xRR[0],virA);

  free_block(gARAR); //!

  C_DGEMM('N','N',occA,virA*occA*virA,occA,-1.0,xAA[0],occA,
    tARAR[0],virA*occA*virA,1.0,s2ARAR[0],virA*occA*virA);

  C_DGEMM('N','N',occA*virA*occA,virA,virA,-1.0,tARAR[0],virA,
    xRR[0],virA,1.0,s2ARAR[0],virA);

  free_block(xAA);
  free_block(xRR);

  sARRA = block_matrix(occA*virA,occA*virA); //!
  double **tARRA = block_matrix(occA*virA,occA*virA); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      sARRA[a1r2][a2r1] = sARAR[a1r1][a2r2];
      tARRA[a1r2][a2r1] = tARAR[a1r1][a2r2];
  }}}}

  free_block(tARAR); //!
  free_block(sARAR); //!

  xARAR = block_matrix(occA*virA,occA*virA); //!
  vARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,ARAR,(char *) &(vARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,1.0,tARRA[0],occA*virA,
    vARAR[0],occA*virA,0.0,xARAR[0],occA*virA);

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,1.0,xARAR[0],occA*virA,
    sARRA[0],occA*virA,1.0,s2ARAR[0],occA*virA);

  free_block(xARAR); //!

  double **vARRA = block_matrix(occA*virA,occA*virA); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      vARRA[a1r2][a2r1] = vARAR[a1r1][a2r2];
  }}}}

  free_block(vARAR); //!
  xARAR = block_matrix(occA*virA,occA*virA); //!

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,1.0,tARRA[0],occA*virA,
    vARRA[0],occA*virA,0.0,xARAR[0],occA*virA);

  free_block(vARRA); //!
  vARAR = block_matrix(occA*virA,occA*virA); //!

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,1.0,xARAR[0],occA*virA,
    sARRA[0],occA*virA,0.0,vARAR[0],occA*virA);

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      s2ARAR[a1r1][a2r2] += vARAR[a1r2][a2r1];
  }}}}

  free_block(xARAR); //!
  free_block(vARAR); //!

  vARRA = block_matrix(occA*virA,occA*virA); //!
  vAARR = block_matrix(occA*occA,virA*virA); //!
  psio_->read_entry(PSIF_SAPT_CCD,ARAR,(char *) &(vARRA[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));
  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1a2 = a1*occA + a2;
      int r1r2 = r1*virA + r2;
      vAARR[a1a2][r1r2] = vARRA[a1r1][a2r2];
  }}}}

  free_block(vARRA); //!

  sAARR = block_matrix(occA*occA,virA*virA); //!
  psio_->read_entry(PSIF_SAPT_CCD,SARAR,(char *) &(sARRA[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));
  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1a2 = a1*occA + a2;
      int r1r2 = r1*virA + r2;
      sAARR[a1a2][r1r2] = sARRA[a1r1][a2r2];
  }}}}

  free_block(sARRA); //!

  double **tAARR = block_matrix(occA*occA,virA*virA); //!
  psio_->read_entry(PSIF_SAPT_CCD,TARAR,(char *) &(tARRA[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1a2 = a1*occA + a2;
      int r1r2 = r1*virA + r2;
      tAARR[a1a2][r1r2] = tARRA[a1r1][a2r2];
  }}}}

  free_block(tARRA); //!

  double **xAAAA = block_matrix(occA*occA,occA*occA);
  s2AARR = block_matrix(occA*occA,virA*virA); //!

  C_DGEMM('N','T',occA*occA,occA*occA,virA*virA,1.0,tAARR[0],virA*virA,
    vAARR[0],virA*virA,0.0,xAAAA[0],occA*occA);

  C_DGEMM('N','N',occA*occA,virA*virA,occA*occA,0.5,xAAAA[0],occA*occA,
    sAARR[0],virA*virA,0.0,s2AARR[0],virA*virA);

  C_DGEMM('N','T',occA*occA,occA*occA,virA*virA,1.0,sAARR[0],virA*virA,
    vAARR[0],virA*virA,0.0,xAAAA[0],occA*occA);

  C_DGEMM('N','N',occA*occA,virA*virA,occA*occA,0.5,xAAAA[0],occA*occA,
    tAARR[0],virA*virA,1.0,s2AARR[0],virA*virA);

  free_block(vAARR); //!
  free_block(sAARR); //!
  free_block(tAARR); //!
  free_block(xAAAA);

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1a2 = a1*occA + a2;
      int r1r2 = r1*virA + r2;
      s2ARAR[a1r1][a2r2] += s2AARR[a1a2][r1r2];
  }}}}

  free_block(s2AARR); //!

  for(int a1r1=0; a1r1<occA*virA; a1r1++) {
    for(int a2r2=0; a2r2<a1r1; a2r2++) {
      double tval = s2ARAR[a1r1][a2r2] + s2ARAR[a2r2][a1r1];
      s2ARAR[a1r1][a2r2] = tval;
      s2ARAR[a2r2][a1r1] = tval;
  }}

  for(int a1r1=0; a1r1<occA*virA; a1r1++) {
    s2ARAR[a1r1][a1r1] *= 2.0;
  }

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      double denom = evalsA_[a1+foccA_]+evalsA_[a2+foccA_]-evalsA_[r1+noccA_]-
        evalsA_[r2+noccA_];
      s2ARAR[a1r1][a2r2] /= denom;
  }}}}

  sARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,SARAR,(char *) &(sARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  psio_->write_entry(PSIF_SAPT_CCD,SARAR,(char *) &(s2ARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DAXPY((size_t)occA*virA*occA*virA,-1.0,sARAR[0],1,s2ARAR[0],1);
  double RMS = C_DDOT((size_t)occA*virA*occA*virA,s2ARAR[0],1,s2ARAR[0],1);
  RMS /= (double) ((size_t)occA*virA*occA*virA);

  psio_->write_entry(PSIF_SAPT_CCD,SARARerr,(char *) &(s2ARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  free_block(sARAR);  //!
  free_block(s2ARAR); //!

  return (sqrt(RMS));
}

void SAPT2p::disp_s_prep(const char *TAR, const char *TpAR, const char *ThetaARAR, //!
  const char *TARBS, int AAfile, const char *AAints, const char *ARints, const char *RRints,
  int BBfile, const char *BSints, double *evalsA_, int noccA_, int virA, int foccA_,
  int noccB_, int virB, int foccB_)
{
  int occA = noccA_ - foccA_;
  int occB = noccB_ - foccB_;

  double **thARAR = block_matrix(occA*virA,occA*virA); //!
  double **B_p_AR = get_DF_ints_nongimp(AAfile,ARints,foccA_,noccA_,0,virA);
  double **T_p_AR = block_matrix(occA*virA,ndf_);

  psio_->read_entry(PSIF_SAPT_CCD,ThetaARAR,(char *) &(thARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*virA,ndf_,occA*virA,1.0,thARAR[0],occA*virA,
    B_p_AR[0],ndf_,0.0,T_p_AR[0],ndf_);

  free_block(thARAR); //!
  free_block(B_p_AR);

  double **T = block_matrix(occA,virA);

  double **B_p_AA = get_DF_ints_nongimp(AAfile,AAints,foccA_,noccA_,foccA_,noccA_);
  double **B_p_RR = get_DF_ints_nongimp(AAfile,RRints,0,virA,0,virA);

  C_DGEMM('N','T',occA,virA,virA*ndf_,1.0,&(T_p_AR[0][0]),
    virA*ndf_,&(B_p_RR[0][0]),virA*ndf_,0.0,&(T[0][0]),
    virA);

  for (int a=0; a<occA; a++) {
    C_DGEMM('N','T',occA,virA,ndf_,-1.0,&(B_p_AA[a*occA][0]),
      ndf_,&(T_p_AR[a*virA][0]),ndf_,1.0,&(T[0][0]),
      virA);
  }

  for (int a=0; a<occA; a++) {
    for (int r=0; r<virA; r++) {
      double denom = evalsA_[a+foccA_] - evalsA_[r+noccA_];
      T[a][r] /= denom;
  }}

  write_IJKL(T,PSIF_SAPT_CCD,TAR,occA,virA);

  free_block(B_p_AA);
  free_block(B_p_RR);

  double **tARBS = block_matrix(occA*virA,occB*virB); //!
  double **B_p_BS = get_DF_ints_nongimp(BBfile,BSints,foccB_,noccB_,0,virB);

  psio_->read_entry(PSIF_SAPT_CCD,TARBS,(char *) &(tARBS[0][0]),
    occA*virA*occB*virB*(ULI) sizeof(double));

  // Mickey mouse with ndf + 3
  double** T_p_AR2 = block_matrix(occA*virA,ndf_+3);

  C_DGEMM('N','N',occA*virA,ndf_,occB*virB,1.0,tARBS[0],occB*virB,
    B_p_BS[0],ndf_,0.0,T_p_AR2[0],ndf_+3);

  write_IJKL(T_p_AR2,PSIF_SAPT_CCD,TpAR,occA*virA,ndf_+3);

  free_block(T_p_AR);
  free_block(tARBS);  //!
  free_block(B_p_BS);
}

void SAPT2p::vvvv_prep(const char *RRRRp, const char *RRRRm, //!
  double** B_p_RR, int virA, int ndf,
  std::shared_ptr<Matrix> mo2no)
{
  timer_on("v^4 Prep           ");

  std::shared_ptr<Matrix> B2;
  std::shared_ptr<Matrix> B3;
  double** B_p_RR2;
  int virA2;
  if (mo2no) {
    virA2 = mo2no->colspi()[0];
    B3 = std::shared_ptr<Matrix>(new Matrix("B3",virA2 * virA,  ndf));
    B2 = std::shared_ptr<Matrix>(new Matrix("B2",virA2 * virA2, ndf));
    double** B_p_RR3 = B3->pointer();
    B_p_RR2 = B2->pointer();
    double** Vp = mo2no->pointer();
    C_DGEMM('T','N',virA2,virA*ndf,virA,1.0,Vp[0],virA2,B_p_RR[0],virA*ndf,0.0,B_p_RR3[0],virA*ndf);
    for (int a = 0; a < virA2; a++) {
        C_DGEMM('T','N',virA2,ndf,virA,1.0,Vp[0],virA2,B_p_RR3[a*virA],ndf,0.0,B_p_RR2[a*virA2],ndf);
    }
    B3.reset();
  } else {
    B_p_RR2 = B_p_RR;
    virA2 = virA;
  }

  int virtri = virA2*(virA2+1)/2;
  int svirtri = virA2*(virA2-1)/2;
  double **RR = block_matrix(virA2,virA2);
  double *xRR = init_array(virtri);
  double *yRR = init_array(svirtri);

  zero_disk(PSIF_SAPT_CCD,RRRRp,virtri,virtri);
  zero_disk(PSIF_SAPT_CCD,RRRRm,svirtri,svirtri);

  psio_address next_RRRRp = PSIO_ZERO;
  psio_address next_RRRRm = PSIO_ZERO;

  for (int r1=0; r1 < virA2; r1++) {
  for (int r2=0; r2 <= r1; r2++) {

    C_DGEMM('N','T',virA2,virA2,ndf,1.0,&(B_p_RR2[r1*virA2][0]),
      ndf,&(B_p_RR2[r2*virA2][0]),ndf,0.0,
      &(RR[0][0]),virA2);

    for (int r3=0; r3 < virA2; r3++) {
    for (int r4=0; r4 <= r3; r4++) {
      int r3r4 = INDEX(r3,r4);
      xRR[r3r4] = RR[r3][r4] + RR[r4][r3];
    }}
    psio_->write(PSIF_SAPT_CCD,RRRRp,(char *) &(xRR[0]),
      virtri*(ULI) sizeof(double),next_RRRRp,&next_RRRRp);

    if (r1 != r2) {
      for (int r3=0; r3 < virA2; r3++) {
      for (int r4=0; r4 < r3; r4++) {
        int r3r4 = INDEX(r3-1,r4);
        yRR[r3r4] = RR[r3][r4] - RR[r4][r3];
      }}
      psio_->write(PSIF_SAPT_CCD,RRRRm,(char *) &(yRR[0]),
        svirtri*(ULI) sizeof(double),next_RRRRm,&next_RRRRm);
    }

  }}

  free(xRR);
  free(yRR);
  free_block(RR);

  timer_off("v^4 Prep           ");
}
double **SAPT2p::vvvv_ccd(const char *TARAR, const char *RRRRp, const char *RRRRm, //!
  int occA, int virA, std::shared_ptr<Matrix> mo2no)
{
  timer_on("v^4 Term           ");

  double **tARAR;
  int virA2;
  if (mo2no) {
    double** Vp = mo2no->pointer();
    virA2 = mo2no->colspi()[0];
    double** T1 = block_matrix(occA*virA,occA*virA);
    psio_->read_entry(PSIF_SAPT_CCD,TARAR,(char *) &(T1[0][0]),
      occA*virA*occA*virA*(ULI) sizeof(double));
    double** T2 = block_matrix(occA*virA,occA*virA2);
    C_DGEMM('N','N',occA*(size_t)virA*occA,virA2,virA,1.0,T1[0],virA,Vp[0],virA2,0.0,T2[0],virA2);
    free_block(T1);
    tARAR = block_matrix(occA*virA2,occA*virA2);
    #pragma omp parallel for
    for (int i = 0; i < occA; i++) {
        C_DGEMM('T','N',virA2,occA*virA2,virA,1.0,Vp[0],virA2,T2[i*virA],occA*virA2,0.0,tARAR[i*virA2],occA*virA2);
    }
    free_block(T2);
  } else {
    virA2 = virA;
    tARAR = block_matrix(occA*virA,occA*virA);
    psio_->read_entry(PSIF_SAPT_CCD,TARAR,(char *) &(tARAR[0][0]),
      occA*virA*occA*virA*(ULI) sizeof(double));
  }

  int occtri = occA*(occA+1)/2;
  int virtri = virA2*(virA2+1)/2;
  int svirtri = virA2*(virA2-1)/2;

  double **tpAARR = block_matrix(occtri,virtri);
  double **tmAARR = block_matrix(occtri,svirtri);

  for(int a1=0; a1<occA; a1++) {
  for(int a2=0; a2<=a1; a2++) {
    for(int r1=0; r1<virA2; r1++) {
    for(int r2=0; r2<=r1; r2++) {
      int a1a2 = INDEX(a1,a2);
      int r1r2 = INDEX(r1-1,r2);
      int s1s2 = INDEX(r1,r2);
      int a1r1 = a1*virA2 + r1;
      int a2r2 = a2*virA2 + r2;
      int a2r1 = a2*virA2 + r1;
      int a1r2 = a1*virA2 + r2;
      if (r1 != r2) {
        tpAARR[a1a2][s1s2] = 0.5*tARAR[a1r1][a2r2];
        tpAARR[a1a2][s1s2] += 0.5*tARAR[a2r1][a1r2];
        tmAARR[a1a2][r1r2] = 0.5*tARAR[a1r1][a2r2];
        tmAARR[a1a2][r1r2] -= 0.5*tARAR[a2r1][a1r2];
      }
      else {
        tpAARR[a1a2][s1s2] = 0.25*tARAR[a1r1][a2r2];
        tpAARR[a1a2][s1s2] += 0.25*tARAR[a2r1][a1r2];
      }
  }}}}

  free_block(tARAR);

  int blocksize;
  int loopsize;

  if (virA2 % 2 == 0) {
    blocksize = virA2+1;
    loopsize = virtri/blocksize;
  }
  else {
    blocksize = virA2;
    loopsize = virtri/blocksize;
  }

  double **sAARR = block_matrix(occtri,virtri);
  double **vRRRRp[2];
  vRRRRp[0] = block_matrix(blocksize,virtri);
  vRRRRp[1] = block_matrix(blocksize,virtri);

  psio_address next_RRRRp = PSIO_ZERO;

    std::shared_ptr<AIOHandler> aio(new AIOHandler(psio_));

    psio_->read(PSIF_SAPT_CCD,RRRRp,(char *) &(vRRRRp[0][0][0]),
        blocksize*virtri*(ULI) sizeof(double),next_RRRRp,&next_RRRRp);

  for(int r_read=0; r_read<loopsize; r_read++) {
    if (r_read < loopsize-1)
      aio->read(PSIF_SAPT_CCD,RRRRp,(char *) &(vRRRRp[(r_read+1)%2][0][0]),
        blocksize*virtri*(ULI) sizeof(double),next_RRRRp,&next_RRRRp);

    C_DGEMM('N','T',occtri,blocksize,virtri,1.0,tpAARR[0],virtri,
      vRRRRp[r_read%2][0],virtri,1.0,&(sAARR[0][r_read*blocksize]),virtri);

        if (r_read < loopsize-1)
            aio->synchronize();
    }

  free_block(vRRRRp[0]);
  free_block(vRRRRp[1]);
  free_block(tpAARR);

  if (virA2 % 2 == 0) {
    blocksize = virA2-1;
    loopsize = svirtri/blocksize;
  }
  else {
    blocksize = virA2;
    loopsize = svirtri/blocksize;
  }

  double **aAARR = block_matrix(occtri,svirtri);
  double **vRRRRm[2];
  vRRRRm[0] = block_matrix(blocksize,svirtri);
  vRRRRm[1] = block_matrix(blocksize,svirtri);

  psio_address next_RRRRm = PSIO_ZERO;

  psio_->read(PSIF_SAPT_CCD,RRRRm,(char *) &(vRRRRm[0][0][0]),
        blocksize*svirtri*(ULI) sizeof(double),next_RRRRm,&next_RRRRm);

  for(int r_read=0; r_read<loopsize; r_read++) {
    if (r_read < loopsize-1)
      aio->read(PSIF_SAPT_CCD,RRRRm,(char *) &(vRRRRm[(r_read+1)%2][0][0]),
        blocksize*svirtri*(ULI) sizeof(double),next_RRRRm,&next_RRRRm);

    C_DGEMM('N','T',occtri,blocksize,svirtri,1.0,tmAARR[0],svirtri,
      vRRRRm[r_read%2][0],svirtri,1.0,&(aAARR[0][r_read*blocksize]),svirtri);

    if (r_read < loopsize-1)
      aio->synchronize();
  }

  free_block(vRRRRm[0]);
  free_block(vRRRRm[1]);
  free_block(tmAARR);

  double **t2AARR = block_matrix(occA*occA,virA2*virA2);

  for(int a1=0,a1a2=0; a1<occA; a1++) {
  for(int a2=0; a2<occA; a2++,a1a2++) {
    int b1b2 = INDEX(a1,a2);
    for(int r3=0,r3r4=0; r3<virA2; r3++) {
    for(int r4=0; r4<virA2; r4++,r3r4++) {
      int s3s4 = INDEX(r3,r4);
      t2AARR[a1a2][r3r4] = sAARR[b1b2][s3s4];
  }}}}

  for(int a1=0; a1<occA; a1++) {
  for(int a2=0; a2<a1; a2++) {
    int a1a2 = a1*occA + a2;
    int a2a1 = a2*occA + a1;
    int b1b2 = INDEX(a1,a2);
    for(int r3=0; r3<virA2; r3++) {
    for(int r4=0; r4<r3; r4++) {
      int r3r4 = r3*virA2 + r4;
      int r4r3 = r4*virA2 + r3;
      int s3s4 = INDEX(r3-1,r4);
      t2AARR[a1a2][r3r4] += aAARR[b1b2][s3s4];
      t2AARR[a2a1][r4r3] += aAARR[b1b2][s3s4];
      t2AARR[a1a2][r4r3] -= aAARR[b1b2][s3s4];
      t2AARR[a2a1][r3r4] -= aAARR[b1b2][s3s4];
  }}}}

  for(int a1=0; a1<occA; a1++) {
    int a1a1 = a1*occA + a1;
    int b1b1 = INDEX(a1,a1);
    for(int r3=0; r3<virA2; r3++) {
    for(int r4=0; r4<r3; r4++) {
      int r3r4 = r3*virA2 + r4;
      int r4r3 = r4*virA2 + r3;
      int s3s4 = INDEX(r3-1,r4);
      t2AARR[a1a1][r3r4] += aAARR[b1b1][s3s4];
      t2AARR[a1a1][r4r3] -= aAARR[b1b1][s3s4];
  }}}

  free_block(sAARR);
  free_block(aAARR);

  double** ret;
  if (mo2no) {
    double** Vp = mo2no->pointer();
    double** T1 = block_matrix(occA*occA,virA*virA2);
    #pragma omp parallel for
    for (size_t ij = 0L; ij < occA*occA; ij++) {
        C_DGEMM('N','N',virA,virA2,virA2,1.0,Vp[0],virA2,t2AARR[ij],virA2,0.0,T1[ij],virA2);
    }
    free_block(t2AARR);
    ret = block_matrix(occA*occA,virA*virA);
    C_DGEMM('N','T',occA*(size_t)occA*virA,virA,virA2,1.0,T1[0],virA2,Vp[0],virA2,0.0,ret[0],virA);
    free_block(T1);
  } else {
    ret = t2AARR;
  }

  timer_off("v^4 Term           ");
  return(ret);
}

void SAPT2p::ccd_prep(const char *TARAR, const char *ThetaARAR, const char *GARAR, const char *GARRA, //!
  const char *AAAA, const char *ARAR, const char *AARR, const char *RRRRp, const char *RRRRm,
  int DFfile, const char *AAints, const char *ARints, const char *RRints, double *evals,
  int noccA, int virA, int foccA, std::shared_ptr<Matrix> mo2no, const char *T2ARAR)
{
  int occA = noccA - foccA;

  double **vAAAA = block_matrix(occA*occA,occA*occA);
  double **B_p_AA = get_DF_ints_nongimp(DFfile,AAints,foccA,noccA,foccA,noccA);

  for(int a1=0,a1a2=0; a1<occA; a1++) {
  for(int a2=0; a2<occA; a2++,a1a2++) {
    C_DGEMM('N','T',occA,occA,ndf_,1.0,&(B_p_AA[a1*occA][0]),
      ndf_,&(B_p_AA[a2*occA][0]),ndf_,0.0,
      &(vAAAA[a1a2][0]),occA);
  }}

  psio_->write_entry(PSIF_SAPT_CCD,AAAA,(char *) &(vAAAA[0][0]),
    occA*occA*occA*occA*(ULI) sizeof(double));
  free_block(vAAAA);


  double **B_p_RR = get_DF_ints_nongimp(DFfile,RRints,0,virA,0,virA);

  vvvv_prep(RRRRp,RRRRm,B_p_RR,virA,ndf_,mo2no); //!!

  double **vAARR = block_matrix(occA*virA,occA*virA); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    C_DGEMM('N','T',occA,virA,ndf_,1.0,&(B_p_AA[a1*occA][0]),
      ndf_,&(B_p_RR[r1*virA][0]),ndf_,0.0,
      &(vAARR[a1r1][0]),virA);
  }}

  psio_->write_entry(PSIF_SAPT_CCD,AARR,(char *) &(vAARR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  free_block(B_p_AA);
  free_block(B_p_RR);

  double **B_p_AR = get_DF_ints_nongimp(DFfile,ARints,foccA,noccA,0,virA);
  double **vARAR = block_matrix(occA*virA,occA*virA); //!

  C_DGEMM('N','T',occA*virA,occA*virA,ndf_,1.0,&(B_p_AR[0][0]),
    ndf_,&(B_p_AR[0][0]),ndf_,0.0,&(vARAR[0][0]),
    occA*virA);

  psio_->write_entry(PSIF_SAPT_CCD,ARAR,(char *) &(vARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  free_block(B_p_AR);

  double **gARAR = block_matrix(occA*virA,occA*virA); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      gARAR[a1r1][a2r2] = 2.0*vARAR[a1r1][a2r2] - vAARR[a1r1][a2r2];
  }}}}

  psio_->write_entry(PSIF_SAPT_CCD,GARRA,(char *) &(gARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  free_block(vAARR); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      gARAR[a1r1][a2r2] = 2.0*vARAR[a1r1][a2r2] - vARAR[a1r2][a2r1];
  }}}}

  psio_->write_entry(PSIF_SAPT_CCD,GARAR,(char *) &(gARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      double denom = evals[a1+foccA]+evals[a2+foccA]-evals[r1+noccA]-
        evals[r2+noccA];
      vARAR[a1r1][a2r2] /= denom;
  }}}}

  psio_->write_entry(PSIF_SAPT_CCD,TARAR,(char *) &(vARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      gARAR[a1r1][a2r2] = 2.0*vARAR[a1r1][a2r2] - vARAR[a1r2][a2r1];
  }}}}

  psio_->write_entry(PSIF_SAPT_CCD,ThetaARAR,(char *) &(gARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  free_block(vARAR); //!
  free_block(gARAR); //!

  double **t2AARR = vvvv_ccd(TARAR,RRRRp,RRRRm,occA,virA,mo2no); //!!
  double **tARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,TARAR,(char *) &(tARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  double **tAARR = block_matrix(occA*occA,virA*virA); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1a2 = a1*occA + a2;
      int r1r2 = r1*virA + r2;
      tAARR[a1a2][r1r2] = tARAR[a1r1][a2r2];
  }}}}

  vAAAA = block_matrix(occA*occA,occA*occA);

  psio_->read_entry(PSIF_SAPT_CCD,AAAA,(char *) &(vAAAA[0][0]),
    occA*occA*occA*occA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*occA,virA*virA,occA*occA,0.5,vAAAA[0],occA*occA,
    tAARR[0],virA*virA,0.5,t2AARR[0],virA*virA);

  free_block(tAARR); //!
  free_block(vAAAA);

  double **t2ARAR = block_matrix(occA*virA,occA*virA); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1a2 = a1*occA + a2;
      int r1r2 = r1*virA + r2;
      t2ARAR[a1r1][a2r2] = t2AARR[a1a2][r1r2];
  }}}}

  free_block(t2AARR); //!

  vAARR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,AARR,(char *) &(vAARR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','T',occA*virA,occA*virA,occA*virA,-1.0,tARAR[0],occA*virA,
    vAARR[0],occA*virA,1.0,t2ARAR[0],occA*virA);

  double **tARRA = block_matrix(occA*virA,occA*virA); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      tARRA[a1r2][a2r1] = tARAR[a1r1][a2r2];
  }}}}

  C_DGEMM('N','T',occA*virA,occA*virA,occA*virA,-1.0,tARRA[0],occA*virA,
    vAARR[0],occA*virA,0.0,tARAR[0],occA*virA);

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      t2ARAR[a1r1][a2r2] += tARAR[a1r2][a2r1];
  }}}}

  free_block(tARRA); //!

  vARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,ARAR,(char *) &(vAARR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  psio_->read_entry(PSIF_SAPT_CCD,ThetaARAR,(char *) &(tARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','T',occA*virA,occA*virA,occA*virA,1.0,tARAR[0],occA*virA,
    vARAR[0],occA*virA,1.0,t2ARAR[0],occA*virA);

  free_block(vAARR); //!
  free_block(vARAR); //!
  free_block(tARAR); //!

  for(int a1r1=0; a1r1<occA*virA; a1r1++) {
    for(int a2r2=0; a2r2<a1r1; a2r2++) {
      double tval = t2ARAR[a1r1][a2r2] + t2ARAR[a2r2][a1r1];
      t2ARAR[a1r1][a2r2] = tval;
      t2ARAR[a2r2][a1r1] = tval;
  }}

  for(int a1r1=0; a1r1<occA*virA; a1r1++) {
    t2ARAR[a1r1][a1r1] *= 2.0;
  }

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      double denom = evals[a1+foccA]+evals[a2+foccA]-evals[r1+noccA]-
        evals[r2+noccA];
      t2ARAR[a1r1][a2r2] /= denom;
  }}}}

  psio_->write_entry(PSIF_SAPT_CCD,T2ARAR,(char *) &(t2ARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  free_block(t2ARAR); //!
}

void SAPT2p::ccd_iterate(const char *TARAR, const char *TARARerr, const char *ThetaARAR, //!
  const char *GARAR, const char *GARRA, const char *AAAA, const char *ARAR, const char *AARR, const char *RRRRp,
  const char *RRRRm, double *evals, int noccA, int virA, int foccA, std::shared_ptr<Matrix> mo2no)
{
  int occA = noccA - foccA;

  if (print_) {
    outfile->Printf("Iter       Energy (H)          dE (H)             RMS (H)\n");

  }

  int iter = 1;
  double E_old=0.0, E_new=0.0, RMS=0.0;

  SAPTDIIS diis(PSIF_SAPT_CCD,TARAR,TARARerr,(size_t)occA*virA*occA*virA,
    max_ccd_vecs_,psio_);

  do {
    E_new = ccd_energy(TARAR,GARAR,occA,virA);
    outfile->Printf("%4d %16.8lf %17.9lf %17.9lf",iter,E_new,E_old-E_new,RMS);


    if (iter > 1 && (fabs(E_old-E_new) < ccd_e_conv_ &&
      RMS < ccd_t_conv_)) {
      if (iter > min_ccd_vecs_) {
        outfile->Printf("  DIIS\n");
      }
      break;
    }

    E_old = E_new;

timer_on("CCD Amps           ");
    RMS = ccd_amplitudes(TARAR,TARARerr,ThetaARAR,GARAR,GARRA,AAAA,ARAR,AARR,
      RRRRp,RRRRm,evals,noccA,virA,foccA,mo2no);
timer_off("CCD Amps           ");

    diis.store_vectors();
    if (iter > min_ccd_vecs_) {
      diis.get_new_vector();
      outfile->Printf("  DIIS\n");
    }
    else {
      outfile->Printf("\n");
    }

    iter++;
  }
  while(iter<ccd_maxiter_+1);

  outfile->Printf("\n");
}

double SAPT2p::ccd_energy(const char *TARAR, const char *GARAR, int occA, int virA) //!
{
  double **gARAR = block_matrix(occA*virA,occA*virA); //!
  double **tARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,GARAR,(char *) &(gARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));
  psio_->read_entry(PSIF_SAPT_CCD,TARAR,(char *) &(tARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  double energy = C_DDOT((size_t)occA*virA*occA*virA,&(gARAR[0][0]),1,
    &(tARAR[0][0]),1);

  free_block(gARAR); //!
  free_block(tARAR); //!

  return(energy);
}

double SAPT2p::ccd_amplitudes(const char *TARAR, const char *TARARerr, const char *ThetaARAR, //!
  const char *GARAR, const char *GARRA, const char *AAAA, const char *ARAR, const char *AARR, const char *RRRRp,
  const char *RRRRm, double *evals, int noccA, int virA, int foccA, std::shared_ptr<Matrix> mo2no)
{
  int occA = noccA - foccA;

  double **t2ARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,ARAR,(char *) &(t2ARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DSCAL((size_t)occA*virA*occA*virA,0.5,&(t2ARAR[0][0]),1);

  double **gARRA = block_matrix(occA*virA,occA*virA); //!
  double **tARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,GARRA,(char *) &(gARRA[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));
  psio_->read_entry(PSIF_SAPT_CCD,TARAR,(char *) &(tARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

    C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,1.0,gARRA[0],occA*virA,
      tARAR[0],occA*virA,1.0,t2ARAR[0],occA*virA);

  free_block(gARRA); //!

  double **t2AARR = vvvv_ccd(TARAR,RRRRp,RRRRm,occA,virA,mo2no); //!!
  double **tAARR = block_matrix(occA*occA,virA*virA); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1a2 = a1*occA + a2;
      int r1r2 = r1*virA + r2;
      tAARR[a1a2][r1r2] = tARAR[a1r1][a2r2];
  }}}}

  double **vAAAA = block_matrix(occA*occA,occA*occA);

  psio_->read_entry(PSIF_SAPT_CCD,AAAA,(char *) &(vAAAA[0][0]),
    occA*occA*occA*occA*(ULI) sizeof(double));

  C_DGEMM('N','N',occA*occA,virA*virA,occA*occA,0.5,vAAAA[0],occA*occA,
    tAARR[0],virA*virA,0.5,t2AARR[0],virA*virA);

  free_block(vAAAA);

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1a2 = a1*occA + a2;
      int r1r2 = r1*virA + r2;
      t2ARAR[a1r1][a2r2] += t2AARR[a1a2][r1r2];
  }}}}

  free_block(tAARR); //!
  free_block(t2AARR); //!

  double **tARRA = block_matrix(occA*virA,occA*virA); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      tARRA[a1r2][a2r1] = tARAR[a1r1][a2r2];
  }}}}

  double **vARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,ARAR,(char *) &(vARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

    C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,-1.0,tARRA[0],occA*virA,
      vARAR[0],occA*virA,1.0,t2ARAR[0],occA*virA);

  double **vAARR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,AARR,(char *) &(vAARR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

    C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,-1.0,tARRA[0],occA*virA,
      vAARR[0],occA*virA,0.0,vARAR[0],occA*virA);

  free_block(vAARR); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      t2ARAR[a1r1][a2r2] += vARAR[a1r2][a2r1];
  }}}}

  free_block(vARAR); //!

  double **gARAR = block_matrix(occA*virA,occA*virA); //!
  double **xARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,GARAR,(char *) &(gARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','T',occA*virA,occA*virA,occA*virA,1.0,tARAR[0],occA*virA,
    gARAR[0],occA*virA,0.0,xARAR[0],occA*virA);

  // => Moved up from below free xARAR to reduce memory <= //

  double **xAA = block_matrix(occA,occA);
  double **xRR = block_matrix(virA,virA);

  C_DGEMM('N','T',occA,occA,virA*occA*virA,1.0,tARAR[0],virA*occA*virA,
    gARAR[0],virA*occA*virA,0.0,xAA[0],occA);

  C_DGEMM('T','N',virA,virA,occA*virA*occA,1.0,gARAR[0],virA,
    tARAR[0],virA,0.0,xRR[0],virA);

  free_block(gARAR); //!

  C_DGEMM('N','N',occA,virA*occA*virA,occA,-1.0,xAA[0],occA,
    tARAR[0],virA*occA*virA,1.0,t2ARAR[0],virA*occA*virA);

  C_DGEMM('N','N',occA*virA*occA,virA,virA,-1.0,tARAR[0],virA,
    xRR[0],virA,1.0,t2ARAR[0],virA);

  free_block(xAA);
  free_block(xRR);

  // => End Move <= //

  double **thARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,ThetaARAR,(char *) &(thARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  C_DGEMM('N','T',occA*virA,occA*virA,occA*virA,0.5,thARAR[0],occA*virA,
    xARAR[0],occA*virA,1.0,t2ARAR[0],occA*virA);

  free_block(thARAR); //!

  C_DGEMM('N','T',occA*virA,occA*virA,occA*virA,-0.5,xARAR[0],occA*virA,
    tARRA[0],occA*virA,1.0,t2ARAR[0],occA*virA);

  vARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,ARAR,(char *) &(vARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

    C_DGEMM('N','T',occA*virA,occA*virA,occA*virA,1.0,tARRA[0],occA*virA,
      vARAR[0],occA*virA,0.0,xARAR[0],occA*virA);

    C_DGEMM('N','T',occA*virA,occA*virA,occA*virA,0.5,xARAR[0],occA*virA,
      tARRA[0],occA*virA,1.0,t2ARAR[0],occA*virA);

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      xARAR[a1r2][a2r1] = vARAR[a1r1][a2r2];
  }}}}

  free_block(vARAR); //! // Reduces memory

  double **yARAR = block_matrix(occA*virA,occA*virA); //!

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,1.0,tARRA[0],occA*virA,
    xARAR[0],occA*virA,0.0,yARAR[0],occA*virA);

  C_DGEMM('N','N',occA*virA,occA*virA,occA*virA,0.5,yARAR[0],occA*virA,
    tARRA[0],occA*virA,0.0,xARAR[0],occA*virA);

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      t2ARAR[a1r1][a2r2] += xARAR[a1r2][a2r1];
  }}}}

  free_block(yARAR); //!
  free_block(xARAR); //!
  free_block(tARRA); //!

  vARAR = block_matrix(occA*virA,occA*virA); //! // Reduces memory

  psio_->read_entry(PSIF_SAPT_CCD,ARAR,(char *) &(vARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double)); // Reduces memory

  tAARR = block_matrix(occA*occA,virA*virA); //!
  vAARR = block_matrix(occA*occA,virA*virA); //!
  double **tAAAA = block_matrix(occA*occA,occA*occA);

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1a2 = a1*occA + a2;
      int r1r2 = r1*virA + r2;
      vAARR[a1a2][r1r2] = vARAR[a1r1][a2r2];
      tAARR[a1a2][r1r2] = tARAR[a1r1][a2r2];
  }}}}

  free_block(vARAR); //!
  free_block(tARAR); //!

  C_DGEMM('N','T',occA*occA,occA*occA,virA*virA,1.0,vAARR[0],virA*virA,
    tAARR[0],virA*virA,0.0,tAAAA[0],occA*occA);

  C_DGEMM('N','N',occA*occA,virA*virA,occA*occA,0.5,tAAAA[0],occA*occA,
    tAARR[0],virA*virA,0.0,vAARR[0],virA*virA);

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1a2 = a1*occA + a2;
      int r1r2 = r1*virA + r2;
      t2ARAR[a1r1][a2r2] += vAARR[a1a2][r1r2];
  }}}}

  free_block(vAARR); //!
  free_block(tAARR); //!
  free_block(tAAAA);

  for(int a1r1=0; a1r1<occA*virA; a1r1++) {
    for(int a2r2=0; a2r2<a1r1; a2r2++) {
      double tval = t2ARAR[a1r1][a2r2] + t2ARAR[a2r2][a1r1];
      t2ARAR[a1r1][a2r2] = tval;
      t2ARAR[a2r2][a1r1] = tval;
  }}

  for(int a1r1=0; a1r1<occA*virA; a1r1++) {
    t2ARAR[a1r1][a1r1] *= 2.0;
  }

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      double denom = evals[a1+foccA]+evals[a2+foccA]-evals[r1+noccA]-
        evals[r2+noccA];
      t2ARAR[a1r1][a2r2] /= denom;
  }}}}

  double **th2ARAR = block_matrix(occA*virA,occA*virA); //!

  for(int a1=0,a1r1=0; a1<occA; a1++) {
  for(int r1=0; r1<virA; r1++,a1r1++) {
    for(int a2=0,a2r2=0; a2<occA; a2++) {
    for(int r2=0; r2<virA; r2++,a2r2++) {
      int a1r2 = a1*virA + r2;
      int a2r1 = a2*virA + r1;
      th2ARAR[a1r1][a2r2] = 2.0*t2ARAR[a1r1][a2r2] - t2ARAR[a1r2][a2r1];
  }}}}

  tARAR = block_matrix(occA*virA,occA*virA); //!

  psio_->read_entry(PSIF_SAPT_CCD,TARAR,(char *) &(tARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  psio_->write_entry(PSIF_SAPT_CCD,TARAR,(char *) &(t2ARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  psio_->write_entry(PSIF_SAPT_CCD,ThetaARAR,(char *) &(th2ARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  free_block(th2ARAR); //!

  C_DAXPY((size_t)occA*virA*occA*virA,-1.0,t2ARAR[0],1,tARAR[0],1);

  double RMS = C_DDOT((size_t)occA*virA*occA*virA,tARAR[0],1,tARAR[0],1);
  RMS /= (double) ((size_t)occA*virA*occA*virA);

  psio_->write_entry(PSIF_SAPT_CCD,TARARerr,(char *) &(tARAR[0][0]),
    occA*virA*occA*virA*(ULI) sizeof(double));

  free_block(tARAR); //!
  free_block(t2ARAR); //!

  return(sqrt(RMS));
}

// => Utility Crap <= //

double **SAPT2p::read_IJKL(int filenum, char *label, int length_IJ,
  int length_KL)
{
  double **A = block_matrix(length_IJ,length_KL);

  psio_->read_entry(filenum,label,(char *) &(A[0][0]),
                  sizeof(double)*length_IJ*(ULI) length_KL);

  return(A);
}

void SAPT2p::write_IJKL(double **A, int filenum, const char *label, int length_IJ,
  int length_KL)
{
  psio_->write_entry(filenum,label,(char *) &(A[0][0]),
                  sizeof(double)*length_IJ*(ULI) length_KL);

  free_block(A);
}

SAPTDIIS::SAPTDIIS(int ampfile, const char *amplabel, const char *errlabel, size_t length,
                   int maxvec, std::shared_ptr<PSIO> psio) : psio_(psio), vec_label_(amplabel), err_label_(errlabel)
{
    diis_file_ = 56;
    psio_->open(diis_file_,0);

    max_diis_vecs_ = maxvec;

    filenum_ = ampfile;

    vec_length_ = length;

    curr_vec_ = 0;
    num_vecs_ = 0;
}

SAPTDIIS::~SAPTDIIS()
{
    psio_->close(diis_file_,0);
}

void SAPTDIIS::store_vectors()
{
    char *diis_vec_label = get_vec_label(curr_vec_);
    char *diis_err_label = get_err_label(curr_vec_);
    curr_vec_ = (curr_vec_+1)%max_diis_vecs_;
    num_vecs_++;
    if (num_vecs_ > max_diis_vecs_) num_vecs_ = max_diis_vecs_;

    double *vec = init_array(vec_length_);

    psio_->read_entry(filenum_,vec_label_,(char *) &(vec[0]),
      vec_length_*(ULI) sizeof(double));
    psio_->write_entry(diis_file_,diis_vec_label,(char *) &(vec[0]),
      vec_length_*(ULI) sizeof(double));

    psio_->read_entry(filenum_,err_label_,(char *) &(vec[0]),
      vec_length_*(ULI) sizeof(double));
    psio_->write_entry(diis_file_,diis_err_label,(char *) &(vec[0]),
      vec_length_*(ULI) sizeof(double));

    free(vec);
    free(diis_vec_label);
    free(diis_err_label);
}

void SAPTDIIS::get_new_vector()
{
    int *ipiv;
    double *Cvec;
    double **Bmat;

    ipiv = init_int_array(num_vecs_+1);
    Bmat = block_matrix(num_vecs_+1,num_vecs_+1);
    Cvec = (double *) malloc((num_vecs_+1)*sizeof(double));

    double *vec_i = init_array(vec_length_);
    double *vec_j = init_array(vec_length_);

    for (int i=0; i<num_vecs_; i++) {
      char *err_label_i = get_err_label(i);
      psio_->read_entry(diis_file_,err_label_i,(char *) &(vec_i[0]),
        vec_length_*(ULI) sizeof(double));
      for (int j=0; j<=i; j++) {
        char *err_label_j = get_err_label(j);
        psio_->read_entry(diis_file_,err_label_j,(char *) &(vec_j[0]),
          vec_length_*(ULI) sizeof(double));
        Bmat[i][j] = Bmat[j][i] = C_DDOT(vec_length_,vec_i,1,vec_j,1);
        free(err_label_j);
      }
      free(err_label_i);
    }

    for (int i=0; i<num_vecs_; i++) {
      Bmat[num_vecs_][i] = -1.0;
      Bmat[i][num_vecs_] = -1.0;
      Cvec[i] = 0.0;
    }

    Bmat[num_vecs_][num_vecs_] = 0.0;
    Cvec[num_vecs_] = -1.0;

    C_DGESV(num_vecs_+1,1,&(Bmat[0][0]),num_vecs_+1,&(ipiv[0]),&(Cvec[0]),
      num_vecs_+1);

    memset(vec_j,'\0',sizeof(double)*vec_length_);

    for (int i=0; i<num_vecs_; i++) {
      char *vec_label_i = get_vec_label(i);
      psio_->read_entry(diis_file_,vec_label_i,(char *) &(vec_i[0]),
        vec_length_*(ULI) sizeof(double));
      C_DAXPY(vec_length_,Cvec[i],vec_i,1,vec_j,1);
      free(vec_label_i);
    }

    psio_->write_entry(filenum_,vec_label_,(char *) &(vec_j[0]),
      vec_length_*(ULI) sizeof(double));

    free(vec_i);
    free(vec_j);

    free(ipiv);
    free(Cvec);
    free_block(Bmat);
}

char *SAPTDIIS::get_err_label(int num)
{
    char *label = (char *) malloc(16*sizeof(char));
    sprintf(label,"Error vector %2d",num);
    return(label);
}

char *SAPTDIIS::get_vec_label(int num)
{
    char *label = (char *) malloc(10*sizeof(char));
    sprintf(label,"Vector %2d",num);
    return(label);
}

std::shared_ptr<Matrix> SAPT2p::mo2no(int ampfile, const char* VV_opdm, int nvir, double cutoff)
{
    std::shared_ptr<Matrix> D(new Matrix("D", nvir, nvir));
    double** Dp = D->pointer();
    psio_->read_entry(ampfile,VV_opdm,(char *) Dp[0],
      sizeof(double)*nvir*nvir);
    D->scale(2.0);

    std::shared_ptr<Matrix> V(new Matrix("V", nvir, nvir));
    std::shared_ptr<Vector> d(new Vector("d", nvir));
    D->diagonalize(V,d,descending);
    D.reset();

    int nno = 0;
    double* dp = d->pointer();
    for (int a=0; a<nvir; a++) {
        if (dp[a] > cutoff) {
            nno++;
        }
    }

    if (options_.get_bool("BENCH")) {
        FILE* fh = fopen(VV_opdm, "w");
        fwrite((void*) dp, sizeof(double), nvir, fh);
        fclose(fh);
    }

    std::shared_ptr<Matrix> U(new Matrix("U",nvir,nno));
    double** Up = U->pointer();
    double** Vp = V->pointer();
    int offset = 0;
    for (int a=0; a<nvir; a++) {
        if (dp[a] > cutoff) {
            C_DCOPY(nvir,&Vp[0][a],nvir,&Up[0][offset],nno);
            offset++;
        }
    }

    return U;
}

}}
