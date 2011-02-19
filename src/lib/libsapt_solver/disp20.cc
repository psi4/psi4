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
#include "sapt0.h"
#include "sapt2p.h"

namespace psi { namespace sapt {

void SAPT0::disp20()
{
  double energy=0.0;
  long int avail_mem = params_.memory - sizeof(double)*calc_info_.noccB*
    calc_info_.nvirB*(long int) calc_info_.nrio;

  if (params_.print)
    fprintf(outfile,"Begining Disp20 Calculation\n");

  long int temp_size = avail_mem / (sizeof(double) * (2*calc_info_.noccB*
    calc_info_.nvirB + (long int) calc_info_.nrio));

  if (temp_size < 1) {
    fprintf(outfile,"Not enough memory in Disp20\n\n");
    exit(0);
  }

  if (temp_size > calc_info_.noccA*calc_info_.nvirA) 
    temp_size = calc_info_.noccA*calc_info_.nvirA;

  double **B_p_AR = block_matrix(temp_size,calc_info_.nrio);
  double **B_p_BS = get_BS_ints(0);
  double **ARBS = block_matrix(temp_size,calc_info_.noccB*calc_info_.nvirB);
  double **tARBS = block_matrix(temp_size,calc_info_.noccB*calc_info_.nvirB);

  int blocks = (calc_info_.noccA*calc_info_.nvirA)/temp_size;
  if ((calc_info_.noccA*calc_info_.nvirA)%temp_size) blocks++;

  if (params_.print) {
    fprintf(outfile,"T2 ARBS Amplitudes formed in %d chunks\n\n",blocks);
    fflush(outfile);
  }

  psio_address next_PSIF_DF_AR = PSIO_ZERO;
  psio_address next_PSIF_tARBS = PSIO_ZERO;

  for (int t_ar=0; t_ar<blocks; t_ar++) {
    int ar_start = temp_size*t_ar;
    int ar_stop = temp_size*(t_ar+1);
    if (ar_stop > calc_info_.noccA*calc_info_.nvirA)
      ar_stop = calc_info_.noccA*calc_info_.nvirA;

    psio_->read(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *) &(B_p_AR[0][0]),
      sizeof(double)*(ar_stop-ar_start)*(ULI) calc_info_.nrio,next_PSIF_DF_AR,
      &next_PSIF_DF_AR);

    C_DGEMM('N','T',ar_stop-ar_start,calc_info_.noccB*calc_info_.nvirB,
      calc_info_.nrio,1.0,&(B_p_AR[0][0]),calc_info_.nrio,&(B_p_BS[0][0]),
      calc_info_.nrio,0.0,&(ARBS[0][0]),calc_info_.noccB*calc_info_.nvirB);
    
    #pragma omp for schedule(static)
      for (int ar=ar_start; ar<ar_stop; ar++) {
        int a = ar/calc_info_.nvirA;
        int r = ar%calc_info_.nvirA + calc_info_.noccA;
        for (int b=0,bs=0; b<calc_info_.noccB; b++) {
          for (int s=calc_info_.noccB; s<calc_info_.nmo; s++,bs++) {
            double denom = calc_info_.evalsA[a]+calc_info_.evalsB[b]-
                           calc_info_.evalsA[r]-calc_info_.evalsB[s];
            tARBS[ar-ar_start][bs] = ARBS[ar-ar_start][bs]/denom;
        }}
      }

    energy += C_DDOT((ar_stop-ar_start)*calc_info_.noccB*calc_info_.nvirB,
              &(ARBS[0][0]),1,&(tARBS[0][0]),1);

    psio_->write(PSIF_SAPT_AMPS,"T ARBS Amplitudes",(char *) &(tARBS[0][0]),
      sizeof(double)*(ar_stop-ar_start)*calc_info_.noccB*
      (ULI) calc_info_.nvirB,next_PSIF_tARBS,&next_PSIF_tARBS);
  }

  free_block(ARBS);
  free_block(tARBS);
  free_block(B_p_AR);
  free_block(B_p_BS);

  results_.disp20 = energy*4.0;
  
  if (params_.print) {
    fprintf(outfile,"Dispersion Energy     = %18.12lf mH\n\n",energy*4000.0);
    fflush(outfile);
  }
}

void SAPT2p::disp20()
{
  double energy=0.0;

  if (params_.print)
    fprintf(outfile,"Begining Disp20 Calculation\n\n");

  double **T_p_AR = read_IJKL(PSIF_SAPT_AMPS,"T(BS) AR",calc_info_.noccA*
    calc_info_.nvirA,calc_info_.nrio);
  double **B_p_AR = get_AR_ints(1);

  energy = C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nrio,
    B_p_AR[0],1,T_p_AR[0],1);

  free_block(B_p_AR);
  free_block(T_p_AR);

  results_.disp20 = energy*4.0;
  
  if (params_.print) {
    fprintf(outfile,"Dispersion Energy     = %18.12lf mH\n\n",energy*4000.0);
    fflush(outfile);
  }
}

}}
