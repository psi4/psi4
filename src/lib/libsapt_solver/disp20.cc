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
#include "sapt0.h"


namespace psi { namespace sapt {

void SAPT0::disp20()
{
  if (params_.logfile) {
    fprintf(params_.logfilename," Disp20\n");
    fprintf(params_.logfilename,"--------\n\n");
    fflush(params_.logfilename);
  }
  time_t start;
  time_t stop;

  double energy=0.0;
  double avail_mem = params_.memory - params_.dfbs_mem; 

  fprintf(outfile,"Begining Disp20 Calculation\n");

  int temp_size = (int) ((avail_mem) / ((double) sizeof(double) * ((double) 
                  (2*calc_info_.noccB*calc_info_.nvirB + calc_info_.nrio))));

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

  fprintf(outfile,"T2 ARBS Amplitudes formed in %d chunks\n\n",blocks);
  fflush(outfile);

  if (params_.logfile) {
    fprintf(params_.logfilename,"     Block      Start       Stop\n");
    fprintf(params_.logfilename,"    -------  ----------  ----------\n");
    for (int t_ar=0; t_ar<blocks; t_ar++) {
      int ar_start = temp_size*t_ar;
      int ar_stop = temp_size*(t_ar+1);
      if (ar_stop > calc_info_.noccA*calc_info_.nvirA)
        ar_stop = calc_info_.noccA*calc_info_.nvirA;
      fprintf(params_.logfilename,"       %3d   %10d  %10d\n",t_ar,ar_start,ar_stop);
    }
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  psio_address next_PSIF_DF_AR = PSIO_ZERO;
  psio_address next_PSIF_tARBS = PSIO_ZERO;

  for (int t_ar=0; t_ar<blocks; t_ar++) {
    int ar_start = temp_size*t_ar;
    int ar_stop = temp_size*(t_ar+1);
    if (ar_stop > calc_info_.noccA*calc_info_.nvirA)
      ar_stop = calc_info_.noccA*calc_info_.nvirA;

    if (params_.logfile) {
      fprintf(params_.logfilename,"    Starting Block %3d ... ",t_ar);
      start = time(NULL);
      fflush(params_.logfilename);
    }

    psio_read(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *) &(B_p_AR[0][0]),
      sizeof(double)*(ar_stop-ar_start)*calc_info_.nrio,next_PSIF_DF_AR,
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

    psio_write(PSIF_SAPT_AMPS,"T2 ARBS Amplitudes",(char *) &(tARBS[0][0]),
      sizeof(double)*(ar_stop-ar_start)*calc_info_.noccB*calc_info_.nvirB,
      next_PSIF_tARBS,&next_PSIF_tARBS);

    if (params_.logfile) {
      stop = time(NULL);
      fprintf(params_.logfilename,"finished %14ld seconds\n",stop-start);
      fflush(params_.logfilename);
    }
  }

  free_block(ARBS);
  free_block(tARBS);
  free_block(B_p_AR);
  free_block(B_p_BS);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  results_.disp20 = energy*4.0;
  
  fprintf(outfile,"Dispersion Energy     = %18.12lf mH\n\n",energy*4000.0);
  fflush(outfile);
}
void SAPT0::theta_ar()
{
  if (params_.logfile) {
    fprintf(params_.logfilename," Theta (BS) AR\n");
    fprintf(params_.logfilename,"---------------\n\n");
    fflush(params_.logfilename);
  }
  time_t start;
  time_t stop;

  double NA = 1.0 / ((double) calc_info_.NA);
  double avail_mem = params_.memory - params_.dfbs_mem; 

  fprintf(outfile,"Forming Theta (BS) AR Intermediates\n");

  int temp_size = (int) ((avail_mem) / (8.0*((double) (calc_info_.noccB*
                  calc_info_.nvirB + calc_info_.nrio))));

  if (temp_size < 1) {
    fprintf(outfile,"Not enough memory in Theta AR\n\n");
    exit(0);
  }

  if (temp_size > calc_info_.noccA*calc_info_.nvirA) 
    temp_size = calc_info_.noccA*calc_info_.nvirA;

  double **B_p_BS = get_BS_ints(1);
  double **theta_AR = block_matrix(temp_size,calc_info_.nrio);
  double **tARBS = block_matrix(temp_size,calc_info_.noccB*calc_info_.nvirB);

  int blocks = (calc_info_.noccA*calc_info_.nvirA)/temp_size;
  if ((calc_info_.noccA*calc_info_.nvirA)%temp_size) blocks++;

  fprintf(outfile,"T2 ARBS Amplitudes read in %d chunks\n\n",blocks);
  fflush(outfile);

  if (params_.logfile) {
    fprintf(params_.logfilename,"     Block      Start       Stop\n");
    fprintf(params_.logfilename,"    -------  ----------  ----------\n");
    for (int t_ar=0; t_ar<blocks; t_ar++) {
      int ar_start = temp_size*t_ar;
      int ar_stop = temp_size*(t_ar+1);
      if (ar_stop > calc_info_.noccA*calc_info_.nvirA)
        ar_stop = calc_info_.noccA*calc_info_.nvirA;
      fprintf(params_.logfilename,"       %3d   %10d  %10d\n",t_ar,ar_start,ar_stop);
    }
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  psio_address next_PSIF_tARBS = PSIO_ZERO;
  psio_address next_PSIF_theta = PSIO_ZERO;

  for (int t_ar=0; t_ar<blocks; t_ar++) {
    int ar_start = temp_size*t_ar;
    int ar_stop = temp_size*(t_ar+1);
    if (ar_stop > calc_info_.noccA*calc_info_.nvirA)
      ar_stop = calc_info_.noccA*calc_info_.nvirA;

    if (params_.logfile) {
      fprintf(params_.logfilename,"    Starting Block %3d ... ",t_ar);
      start = time(NULL);
      fflush(params_.logfilename);
    }

    psio_read(PSIF_SAPT_AMPS,"T2 ARBS Amplitudes",(char *) &(tARBS[0][0]),
      sizeof(double)*(ar_stop-ar_start)*calc_info_.noccB*calc_info_.nvirB,
      next_PSIF_tARBS,&next_PSIF_tARBS);

    C_DGEMM('N','N',(ar_stop-ar_start),calc_info_.nrio,calc_info_.noccB*
      calc_info_.nvirB,1.0,&(tARBS[0][0]),calc_info_.noccB*calc_info_.nvirB,
      &(B_p_BS[0][0]),calc_info_.nrio,0.0,&(theta_AR[0][0]),calc_info_.nrio);

    psio_write(PSIF_SAPT_AMPS,"Theta (BS) AR",(char *) &(theta_AR[0][0]),
      sizeof(double)*(ar_stop-ar_start)*calc_info_.nrio,next_PSIF_theta,
      &next_PSIF_theta);

    if (params_.logfile) {
      stop = time(NULL);
      fprintf(params_.logfilename,"finished %14ld seconds\n",stop-start);
      fflush(params_.logfilename);
    }
  }

  free_block(theta_AR);
  free_block(tARBS);
  free_block(B_p_BS);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }
}

void SAPT0::theta_bs()
{
  if (params_.logfile) {
    fprintf(params_.logfilename," Theta (AR) BS\n");
    fprintf(params_.logfilename,"---------------\n\n");
    fflush(params_.logfilename);
  }
  time_t start;
  time_t stop;

  double NB = 1.0 / ((double) calc_info_.NB);
  double avail_mem = params_.memory - params_.dfbs_mem; 

  fprintf(outfile,"Forming Theta (AR) BS Intermediates\n");

  int temp_size = (int) ((avail_mem) / (8.0*((double) (calc_info_.noccB*
                  calc_info_.nvirB + calc_info_.nrio))));

  if (temp_size < 1) {
    fprintf(outfile,"Not enough memory in Theta BS\n\n");
    exit(0);
  }

  if (temp_size > calc_info_.noccA*calc_info_.nvirA) 
    temp_size = calc_info_.noccA*calc_info_.nvirA;

  double **B_p_AR = block_matrix(temp_size,calc_info_.nrio);
  double **theta_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio);
  double **tARBS = block_matrix(temp_size,calc_info_.noccB*calc_info_.nvirB);

  int blocks = (calc_info_.noccA*calc_info_.nvirA)/temp_size;
  if ((calc_info_.noccA*calc_info_.nvirA)%temp_size) blocks++;

  fprintf(outfile,"T2 ARBS Amplitudes read in %d chunks\n\n",blocks);
  fflush(outfile);

  if (params_.logfile) {
    fprintf(params_.logfilename,"     Block      Start       Stop\n");
    fprintf(params_.logfilename,"    -------  ----------  ----------\n");
    for (int t_ar=0; t_ar<blocks; t_ar++) {
      int ar_start = temp_size*t_ar;
      int ar_stop = temp_size*(t_ar+1);
      if (ar_stop > calc_info_.noccA*calc_info_.nvirA)
        ar_stop = calc_info_.noccA*calc_info_.nvirA;
      fprintf(params_.logfilename,"       %3d   %10d  %10d\n",t_ar,ar_start,ar_stop);
    }
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  psio_address next_PSIF_DF_AR = PSIO_ZERO;
  psio_address next_PSIF_tARBS = PSIO_ZERO;

  for (int t_ar=0; t_ar<blocks; t_ar++) {
    int ar_start = temp_size*t_ar;
    int ar_stop = temp_size*(t_ar+1);
    if (ar_stop > calc_info_.noccA*calc_info_.nvirA)
      ar_stop = calc_info_.noccA*calc_info_.nvirA;

    if (params_.logfile) {
      fprintf(params_.logfilename,"    Starting Block %3d ... ",t_ar);
      start = time(NULL);
      fflush(params_.logfilename);
    }

    psio_read(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",(char *) &(B_p_AR[0][0]),
      sizeof(double)*(ar_stop-ar_start)*calc_info_.nrio,next_PSIF_DF_AR,
      &next_PSIF_DF_AR);

    for (int ar=ar_start; ar<ar_stop; ar++){
      int a = ar/calc_info_.nvirA;
      int r = ar%calc_info_.nvirA+calc_info_.noccA;
      B_p_AR[ar-ar_start][calc_info_.nrio-2] = NB*calc_info_.VBAA[a][r];
    }

    psio_read(PSIF_SAPT_AMPS,"T2 ARBS Amplitudes",(char *) &(tARBS[0][0]),
      sizeof(double)*(ar_stop-ar_start)*calc_info_.noccB*calc_info_.nvirB,
      next_PSIF_tARBS,&next_PSIF_tARBS);

    C_DGEMM('T','N',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,
      (ar_stop-ar_start),1.0,&(tARBS[0][0]),calc_info_.noccB*calc_info_.nvirB,
      &(B_p_AR[0][0]),calc_info_.nrio,1.0,&(theta_BS[0][0]),calc_info_.nrio);

    if (params_.logfile) {
      stop = time(NULL);
      fprintf(params_.logfilename,"finished %14ld seconds\n",stop-start);
      fflush(params_.logfilename);
    }
  }

  psio_write_entry(PSIF_SAPT_AMPS,"Theta (AR) BS",(char *) &(theta_BS[0][0]),
    sizeof(double)*calc_info_.noccB*calc_info_.nvirB*calc_info_.nrio);

  free_block(theta_BS);
  free_block(tARBS);
  free_block(B_p_AR);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }
}

}}
