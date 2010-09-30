/* This function calculates the Disp20 energy */

#define EXTERN

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <cstring>
#include <iostream>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <libipv1/ip_data.gbl>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "sapt2p3.h"

namespace psi { namespace sapt {

void SAPT2p3::elst13()
{
  double elst130, elst130_1, elst130_2, elst130_3;
  double elst103, elst103_1, elst103_2, elst103_3;

  if (params_.print) {
    fprintf(outfile,"Begining Elst13 Calculation\n\n");
    fflush(outfile);
  }

  double **xRR = block_matrix(calc_info_.nvirA,calc_info_.nvirA);
  double **xAA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  double **thetaARAR = read_IJKL(PSIF_SAPT_AMPS,"T2 ARAR Antisym Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*calc_info_.nvirA);
  double **tARAR = read_IJKL(PSIF_SAPT_AMPS,"T ARAR Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccA*calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccA*
    calc_info_.nvirA*calc_info_.nvirA,1.0,&(thetaARAR[0][0]),
    calc_info_.noccA*calc_info_.nvirA*calc_info_.nvirA,&(tARAR[0][0]),
    calc_info_.noccA*calc_info_.nvirA*calc_info_.nvirA,0.0,&(xAA[0][0]),
    calc_info_.noccA);

  C_DGEMM('T','N',calc_info_.nvirA,calc_info_.nvirA,calc_info_.noccA*
    calc_info_.noccA*calc_info_.nvirA,1.0,&(thetaARAR[0][0]),
    calc_info_.nvirA,&(tARAR[0][0]),calc_info_.nvirA,0.0,&(xRR[0][0]),
    calc_info_.nvirA);

  free_block(thetaARAR);
  free_block(tARAR);

  double **yAR = read_IJKL(PSIF_SAPT_AMPS,"Y3 AR Amplitudes",calc_info_.noccA,
    calc_info_.nvirA);

  elst130_1 = 4.0*C_DDOT(calc_info_.nvirA*calc_info_.nvirA,xRR[0],1,
    calc_info_.WBRR[0],1);
  elst130_2 = -4.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,xAA[0],1,
    calc_info_.WBAA[0],1);
  elst130_3 = 4.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA,yAR[0],1,
    calc_info_.CHFA[0],1);

  elst130 = elst130_1 + elst130_2 + elst130_3;

  if (params_.print) {
    fprintf(outfile,"elstr130           = %18.12lf  H\n",elst130);
    fflush(outfile);
  }

  free_block(xRR);
  free_block(xAA);
  free_block(yAR);

  double **xSS = block_matrix(calc_info_.nvirB,calc_info_.nvirB);
  double **xBB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  double **thetaBSBS = read_IJKL(PSIF_SAPT_AMPS,"T2 BSBS Antisym Amplitudes",
    calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB);
  double **tBSBS = read_IJKL(PSIF_SAPT_AMPS,"T BSBS Amplitudes",
    calc_info_.noccB*calc_info_.nvirB,calc_info_.noccB*calc_info_.nvirB);

  C_DGEMM('N','T',calc_info_.noccB,calc_info_.noccB,calc_info_.noccB*
    calc_info_.nvirB*calc_info_.nvirB,1.0,&(thetaBSBS[0][0]),
    calc_info_.noccB*calc_info_.nvirB*calc_info_.nvirB,&(tBSBS[0][0]),
    calc_info_.noccB*calc_info_.nvirB*calc_info_.nvirB,0.0,&(xBB[0][0]),
    calc_info_.noccB);

  C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nvirB,calc_info_.noccB*
    calc_info_.noccB*calc_info_.nvirB,1.0,&(thetaBSBS[0][0]),
    calc_info_.nvirB,&(tBSBS[0][0]),calc_info_.nvirB,0.0,&(xSS[0][0]),
    calc_info_.nvirB);

  free_block(thetaBSBS);
  free_block(tBSBS);

  double **yBS = read_IJKL(PSIF_SAPT_AMPS,"Y3 BS Amplitudes",calc_info_.noccB,
    calc_info_.nvirB);

  elst103_1 = 4.0*C_DDOT(calc_info_.nvirB*calc_info_.nvirB,xSS[0],1,
    calc_info_.WASS[0],1);
  elst103_2 = -4.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,xBB[0],1,
    calc_info_.WABB[0],1);
  elst103_3 = 4.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB,yBS[0],1,
    calc_info_.CHFB[0],1);

  elst103 = elst103_1 + elst103_2 + elst103_3;

  if (params_.print) {
    fprintf(outfile,"elstr103           = %18.12lf  H\n\n",elst103);
    fflush(outfile);
  }

  free_block(xSS);
  free_block(xBB);
  free_block(yBS);

  results_.elst13 = elst130 + elst103;
}

}}
