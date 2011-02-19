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
#include "sapt2p.h"

namespace psi { namespace sapt {

void SAPT2p::disp21()
{
  double d210, d201;
  if (params_.print)
    fprintf(outfile,"Begining Disp21 Calculation\n\n");

  d210 = disp210();

  if (params_.print) {
    fprintf(outfile,"disp210            = %18.12lf  H\n",d210);
    fflush(outfile);
  }

  d201 = disp201();

  if (params_.print) {
    fprintf(outfile,"disp201            = %18.12lf  H\n\n",d201);
    fflush(outfile);
  }

  results_.disp21 = d210 + d201;
}

double SAPT2p::disp210()
{
  double energy = 0.0;

  double **gARBS = read_IJKL(PSIF_SAPT_AMPS,"GxT(AR) ARBS",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);
  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"Theta(AR) AR",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);
  double **B_p_BS = get_BS_ints(0);

  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nri,2.0,&(T_p_AR[0][0]),calc_info_.nrio,
    &(B_p_BS[0][0]),calc_info_.nrio,1.0,&(gARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  free_block(T_p_AR);
  free_block(B_p_BS);

  double **tARBS = read_IJKL(PSIF_SAPT_AMPS,"T ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  energy += 4.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.noccB*
    calc_info_.nvirB,tARBS[0],1,gARBS[0],1);

  free_block(gARBS);
  free_block(tARBS);

  return(energy);
}

double SAPT2p::disp201()
{
  double energy = 0.0;

  double **gARBS = read_IJKL(PSIF_SAPT_AMPS,"GxT(BS) ARBS",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);
  double **B_p_AR = get_AR_ints(0);
  double **T_p_BS = read_IJKL(PSIF_SAPT_AMPS,"Theta(BS) BS",calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nrio);
  
  C_DGEMM('N','T',calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*
    calc_info_.nvirB,calc_info_.nri,2.0,&(B_p_AR[0][0]),calc_info_.nrio,
    &(T_p_BS[0][0]),calc_info_.nrio,1.0,&(gARBS[0][0]),calc_info_.noccB*
    calc_info_.nvirB);

  free_block(B_p_AR);
  free_block(T_p_BS);

  double **tARBS = read_IJKL(PSIF_SAPT_AMPS,"T ARBS Amplitudes",
    calc_info_.noccA*calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  energy += 4.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.noccB*
    calc_info_.nvirB,tARBS[0],1,gARBS[0],1);
  
  free_block(gARBS);
  free_block(tARBS);

  return(energy);
}

}}
