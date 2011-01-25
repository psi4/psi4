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
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "sapt2.h"

namespace psi { namespace sapt {

void SAPT2::elst12()
{
  double elst120, elst120_1, elst120_2, elst120_3;
  double elst102, elst102_1, elst102_2, elst102_3;

  if (params_.print) {
    fprintf(outfile,"Begining Elst12 Calculation\n\n");
    fflush(outfile);
  }

  double **xRR = read_IJKL(PSIF_SAPT_AMPS,"RR MP2 OPDM",calc_info_.nvirA,
    calc_info_.nvirA);
  double **xAA = read_IJKL(PSIF_SAPT_AMPS,"AA MP2 OPDM",calc_info_.noccA,
    calc_info_.noccA);
  double **yAR = read_IJKL(PSIF_SAPT_AMPS,"Y2 AR Amplitudes",calc_info_.noccA,
    calc_info_.nvirA);

  elst120_1 = 2.0*C_DDOT(calc_info_.nvirA*calc_info_.nvirA,xRR[0],1,
    calc_info_.WBRR[0],1);
  elst120_2 = -2.0*C_DDOT(calc_info_.noccA*calc_info_.noccA,xAA[0],1,
    calc_info_.WBAA[0],1);
  elst120_3 = 4.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA,yAR[0],1,
    calc_info_.CHFA[0],1);

  elst120 = elst120_1 + elst120_2 + elst120_3;

  if (params_.print) {
    fprintf(outfile,"elstr120           = %18.12lf  H\n",elst120);
    fflush(outfile);
  }

  free_block(xRR);
  free_block(xAA);
  free_block(yAR);

  double **xSS = read_IJKL(PSIF_SAPT_AMPS,"SS MP2 OPDM",calc_info_.nvirB,
    calc_info_.nvirB);
  double **xBB = read_IJKL(PSIF_SAPT_AMPS,"BB MP2 OPDM",calc_info_.noccB,
    calc_info_.noccB);
  double **yBS = read_IJKL(PSIF_SAPT_AMPS,"Y2 BS Amplitudes",calc_info_.noccB,
    calc_info_.nvirB);

  elst102_1 = 2.0*C_DDOT(calc_info_.nvirB*calc_info_.nvirB,xSS[0],1,
    calc_info_.WASS[0],1);
  elst102_2 = -2.0*C_DDOT(calc_info_.noccB*calc_info_.noccB,xBB[0],1,
    calc_info_.WABB[0],1);
  elst102_3 = 4.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB,yBS[0],1,
    calc_info_.CHFB[0],1);

  elst102 = elst102_1 + elst102_2 + elst102_3;

  if (params_.print) {
    fprintf(outfile,"elstr102           = %18.12lf  H\n\n",elst102);
    fflush(outfile);
  }

  free_block(xSS);
  free_block(xBB);
  free_block(yBS);

  results_.elst12 = elst120 + elst102;
}

}}
