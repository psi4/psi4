/* This function calculates the Exch-Disp20 energy */

#include <cstdio>
#include <cstdlib>
#include <cmath>
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

void SAPT0::exch_disp20()
{
  if (params_.logfile) {
    fprintf(params_.logfilename," Exch-Disp20\n");
    fprintf(params_.logfilename,"-------------\n\n");
    fflush(params_.logfilename);
  }

  fprintf(outfile,"Begining Exch-Disp20 Calculation\n\n");
  fflush(outfile);

  if (params_.logfile) { 
    fprintf(params_.logfilename,"  Starting Exch-Disp 1\n\n"); fflush(params_.logfilename); 
  }
  results_.exch_disp20 = exch_disp_1();
  if (params_.logfile) { 
    fprintf(params_.logfilename,"  Starting Exch-Disp 2\n\n"); fflush(params_.logfilename); 
  }
  results_.exch_disp20 += exch_disp_2();
  if (params_.logfile) { 
    fprintf(params_.logfilename,"  Starting Exch-Disp 3\n\n"); fflush(params_.logfilename); 
  }
  results_.exch_disp20 += exch_disp_3();
  if (params_.logfile) { 
    fprintf(params_.logfilename,"  Starting Exch-Disp 4\n\n"); fflush(params_.logfilename); 
  }
  results_.exch_disp20 += exch_disp_4();
  if (params_.logfile) { 
    fprintf(params_.logfilename,"  Starting Exch-Disp 5\n\n"); fflush(params_.logfilename); 
  }
  results_.exch_disp20 += exch_disp_5();
  if (params_.logfile) { 
    fprintf(params_.logfilename,"  Starting Exch-Disp 6\n\n"); fflush(params_.logfilename); 
  }
  results_.exch_disp20 += exch_disp_6();
  if (params_.logfile) { 
    fprintf(params_.logfilename,"  Starting Exch-Disp 7\n\n"); fflush(params_.logfilename); 
  }
  results_.exch_disp20 += exch_disp_7();
  results_.exch_disp20 *= -2.0;

  fprintf(outfile,"\nExch_Disp Energy = %18.12lf  H\n\n",results_.exch_disp20);
  fflush(outfile);
}

double SAPT0::exch_disp_1()
{
  time_t start;
  time_t stop;
  double vtot = 0.0;
  double h1 = 0.0;
  double h3 = 0.0;
  double enuc, NA, NB;
  
  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **Y_RB = block_matrix(calc_info_.nvirA,calc_info_.noccB);

  if (params_.logfile) { 
    fprintf(params_.logfilename,"    Starting H1\n"); fflush(params_.logfilename); 
  }
  H1(Y_RB);
  if (params_.logfile) {
    fprintf(params_.logfilename,"    Starting Q5\n"); fflush(params_.logfilename);         
  }
  Q5(Y_RB);
  if (params_.logfile) {
    fprintf(params_.logfilename,"    Starting Q7\n"); fflush(params_.logfilename);         
  }
  Q7(Y_RB);
  if (params_.logfile) {
    fprintf(params_.logfilename,"    Starting Q11\n"); fflush(params_.logfilename);         
  }
  Q11(Y_RB);

  double **Y_AS = block_matrix(calc_info_.noccA,calc_info_.nvirB);

  if (params_.logfile) {
    fprintf(params_.logfilename,"    Starting H3\n"); fflush(params_.logfilename);         
  }
  H3(Y_AS);
  if (params_.logfile) {
    fprintf(params_.logfilename,"    Starting Q1\n"); fflush(params_.logfilename);         
  }
  Q1(Y_AS);
  if (params_.logfile) {
    fprintf(params_.logfilename,"    Starting Q3\n"); fflush(params_.logfilename);         
  }
  Q3(Y_AS);
  if (params_.logfile) {
    fprintf(params_.logfilename,"    Starting Q10\n"); fflush(params_.logfilename);         
  }
  Q10(Y_AS);

  if (params_.logfile) { fprintf(params_.logfilename,"\n"); fflush(params_.logfilename); }

  double **B_p_RB = get_RB_ints(1);
  double **B_p_aS = block_matrix(calc_info_.nvirB,calc_info_.nrio);
  double **vRBaS = block_matrix(calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB);
  double **taRBS = block_matrix(calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  psio_address next_PSIF_tARBS = PSIO_ZERO;
  psio_address next_PSIF_DF_AS = PSIO_ZERO;

  for (int a=0; a<calc_info_.noccA; a++) {

    if (params_.logfile) {
      fprintf(params_.logfilename,"    Starting Occupied A %3d ... ",a);
      start = time(NULL);
      fflush(params_.logfilename);
    }

    psio_read(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",(char *) &(B_p_aS[0][0]),
      sizeof(double)*calc_info_.nvirB*calc_info_.nrio,next_PSIF_DF_AS,
      &next_PSIF_DF_AS);

    for (int s=0; s<calc_info_.nvirB; s++){
      B_p_aS[s][calc_info_.nrio-3] = calc_info_.S_AB[a][s+calc_info_.noccB];
      B_p_aS[s][calc_info_.nrio-2] = NB*calc_info_.VBAB[a][s+calc_info_.noccB];
      B_p_aS[s][calc_info_.nrio-1] = enuc*calc_info_.S_AB[a][s+calc_info_.noccB];
    }

    psio_read(PSIF_SAPT_AMPS,"T2 ARBS Amplitudes",(char *) &(taRBS[0][0]),
      sizeof(double)*calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      next_PSIF_tARBS,&next_PSIF_tARBS);

    C_DGEMM('N','T',calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB,
      calc_info_.nrio,1.0,&(B_p_RB[0][0]),calc_info_.nrio,&(B_p_aS[0][0]),
      calc_info_.nrio,0.0,&(vRBaS[0][0]),calc_info_.nvirB);

    vtot += C_DDOT(calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      &(vRBaS[0][0]),1,&(taRBS[0][0]),1);

    for (int r=0; r<calc_info_.nvirA; r++) {
      for (int b=0; b<calc_info_.noccB; b++) {
        h1 += Y_RB[r][b]*C_DDOT(calc_info_.nvirB,&(taRBS[r][b*calc_info_.nvirB]),1,
              &(calc_info_.S_AB[a][calc_info_.noccB]),1);
        h3 += calc_info_.S_AB[r+calc_info_.noccA][b]*C_DDOT(calc_info_.nvirB,
              &(taRBS[r][b*calc_info_.nvirB]),1,&(Y_AS[a][0]),1);
    }}      

    if (params_.logfile) {
      stop = time(NULL);
      fprintf(params_.logfilename,"finished %14ld seconds\n",stop-start);
      fflush(params_.logfilename);
    }

  }


  free_block(Y_RB);
  free_block(Y_AS);
  free_block(B_p_RB);
  free_block(B_p_aS);
  free_block(vRBaS);
  free_block(taRBS);

  fprintf(outfile,"VTOT        Energy = %18.12lf  H\n",vtot);
  fprintf(outfile,"H1          Energy = %18.12lf  H\n",h1);
  fprintf(outfile,"H3          Energy = %18.12lf  H\n",h3);
  fflush(outfile);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  return(vtot + h1 + h3);
}

double SAPT0::exch_disp_2()
{
  double h2 = 0.0;

  double **C_p_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio);

  if (params_.logfile) {
    fprintf(params_.logfilename,"    Starting H2\n"); fflush(params_.logfilename);
  }
  H2(C_p_BS);
  if (params_.logfile) {
    fprintf(params_.logfilename,"    Starting Q6\n"); fflush(params_.logfilename);
  }
  Q6(C_p_BS);
  if (params_.logfile) {
    fprintf(params_.logfilename,"    Starting Q13\n"); fflush(params_.logfilename);
  }
  Q13(C_p_BS);

  if (params_.logfile) {
    fprintf(params_.logfilename,"    Dot with Theta\n"); fflush(params_.logfilename);
  }

  double **theta_BS = block_matrix(calc_info_.noccB*calc_info_.nvirB,
    calc_info_.nrio);

  psio_read_entry(PSIF_SAPT_AMPS,"Theta (AR) BS",(char *) &(theta_BS[0][0]),
    sizeof(double)*calc_info_.noccB*calc_info_.nvirB*calc_info_.nrio);

  double **X_BS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][calc_info_.noccB]),
    calc_info_.nmo,0.0,&(X_BS[0][0]),calc_info_.nvirB);

  double **Z_BS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,4.0,
    &(theta_BS[0][0]),calc_info_.nrio,calc_info_.diagBB,1,0.0,&(Z_BS[0][0]),1);

  h2 = C_DDOT(calc_info_.noccB*calc_info_.nvirB*calc_info_.nrio,&(C_p_BS[0][0]),1,
    &(theta_BS[0][0]),1);
 
  h2 += C_DDOT(calc_info_.noccB*calc_info_.nvirB,&(Z_BS[0][0]),1,&(X_BS[0][0]),1);
 
  free_block(theta_BS);
  free_block(C_p_BS);
  free_block(X_BS);
  free_block(Z_BS);

  fprintf(outfile,"H2          Energy = %18.12lf  H\n",h2);
  fflush(outfile);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  return(h2);
}

double SAPT0::exch_disp_3()
{
  double h4 = 0.0;

  double **C_p_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio);

  if (params_.logfile) {
    fprintf(params_.logfilename,"    Starting H4\n"); fflush(params_.logfilename);
  }
  H4(C_p_AR);
  if (params_.logfile) {
    fprintf(params_.logfilename,"    Starting Q2\n"); fflush(params_.logfilename);
  }
  Q2(C_p_AR);
  if (params_.logfile) {
    fprintf(params_.logfilename,"    Starting Q14\n"); fflush(params_.logfilename);
  }
  Q14(C_p_AR);

  if (params_.logfile) {
    fprintf(params_.logfilename,"    Dot with Theta\n"); fflush(params_.logfilename);
  }

  double **theta_AR = block_matrix(calc_info_.noccA*calc_info_.nvirA,
    calc_info_.nrio);

  psio_read_entry(PSIF_SAPT_AMPS,"Theta (BS) AR",(char *) &(theta_AR[0][0]),
    sizeof(double)*calc_info_.noccA*calc_info_.nvirA*calc_info_.nrio);

  double **X_AR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[calc_info_.noccA][0]),
    calc_info_.nmo,0.0,&(X_AR[0][0]),calc_info_.nvirA);

  double **Z_AR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,4.0,
    &(theta_AR[0][0]),calc_info_.nrio,calc_info_.diagAA,1,0.0,&(Z_AR[0][0]),1);

  h4 = C_DDOT(calc_info_.noccA*calc_info_.nvirA*calc_info_.nrio,&(C_p_AR[0][0]),1,
    &(theta_AR[0][0]),1);

  h4 += C_DDOT(calc_info_.noccA*calc_info_.nvirA,&(Z_AR[0][0]),1,&(X_AR[0][0]),1);

  free_block(theta_AR);
  free_block(C_p_AR);
  free_block(X_AR);
  free_block(Z_AR);

  fprintf(outfile,"H4          Energy = %18.12lf  H\n",h4);
  fflush(outfile);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  return(h4);
}

double SAPT0::exch_disp_4()
{
  time_t start;
  time_t stop;
  double h2 = 0.0;

  double **B_p_AA = get_AA_ints(1);
  double **C_p_aS = block_matrix(calc_info_.nvirB,calc_info_.nrio);
  double **B_p_RB = get_RB_ints(1);
  double **vRBaS = block_matrix(calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB);
  double **taRBS = block_matrix(calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  psio_address next_PSIF_tARBS = PSIO_ZERO;

  for (int a=0; a<calc_info_.noccA; a++) {
  
    if (params_.logfile) {
      fprintf(params_.logfilename,"    Starting Occupied A %3d ... ",a);
      start = time(NULL);
      fflush(params_.logfilename);
    }

    psio_read(PSIF_SAPT_AMPS,"T2 ARBS Amplitudes",(char *) &(taRBS[0][0]),
      sizeof(double)*calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      next_PSIF_tARBS,&next_PSIF_tARBS);

    C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.noccA,
      1.0,&(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
      &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,0.0,&(C_p_aS[0][0]),
      calc_info_.nrio);

    C_DGEMM('N','T',calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB,
      calc_info_.nrio,1.0,&(B_p_RB[0][0]),calc_info_.nrio,&(C_p_aS[0][0]),
      calc_info_.nrio,0.0,&(vRBaS[0][0]),calc_info_.nvirB);

    h2 -= C_DDOT(calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      &(vRBaS[0][0]),1,&(taRBS[0][0]),1);

    if (params_.logfile) {
      stop = time(NULL);
      fprintf(params_.logfilename,"finished %14ld seconds\n",stop-start);
      fflush(params_.logfilename);
    }
  }

  free_block(B_p_AA);
  free_block(C_p_aS);
  free_block(B_p_RB);
  free_block(vRBaS);
  free_block(taRBS);

  fprintf(outfile,"H2*         Energy = %18.12lf  H\n",h2);
  fflush(outfile);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  return(h2);
}

double SAPT0::exch_disp_5()
{
  time_t start;
  time_t stop;
  double h4 = 0.0;
  double enuc, NA, NB;

  NA = 1.0 / ((double) calc_info_.NA);
  NB = 1.0 / ((double) calc_info_.NB);
  enuc = sqrt((calc_info_.enuc_D - calc_info_.enuc_A - calc_info_.enuc_B)*NA*NB);

  double **B_p_BB = get_BB_ints(1);
  double **B_p_aS = block_matrix(calc_info_.nvirB,calc_info_.nrio);
  double **C_p_RB = block_matrix(calc_info_.nvirA*calc_info_.noccB,calc_info_.nrio);
  double **vRBaS = block_matrix(calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB);
  double **taRBS = block_matrix(calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  C_DGEMM('N','N',calc_info_.nvirA,calc_info_.noccB*calc_info_.nrio,calc_info_.noccB,
    1.0,&(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
    &(B_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_RB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  psio_address next_PSIF_DF_AS = PSIO_ZERO;
  psio_address next_PSIF_tARBS = PSIO_ZERO;

  for (int a=0; a<calc_info_.noccA; a++) {
  
    if (params_.logfile) {
      fprintf(params_.logfilename,"    Starting Occupied A %3d ... ",a);
      start = time(NULL);
      fflush(params_.logfilename);
    }

    psio_read(PSIF_SAPT_AB_DF_INTS,"AS RI Integrals",(char *) &(B_p_aS[0][0]),
      sizeof(double)*calc_info_.nvirB*calc_info_.nrio,next_PSIF_DF_AS,
      &next_PSIF_DF_AS);

    for (int s=0; s<calc_info_.nvirB; s++){
      B_p_aS[s][calc_info_.nrio-3] = calc_info_.S_AB[a][s+calc_info_.noccB];
      B_p_aS[s][calc_info_.nrio-2] = NB*calc_info_.VBAB[a][s+calc_info_.noccB];
      B_p_aS[s][calc_info_.nrio-1] = enuc*calc_info_.S_AB[a][s+calc_info_.noccB];
    }

    psio_read(PSIF_SAPT_AMPS,"T2 ARBS Amplitudes",(char *) &(taRBS[0][0]),
      sizeof(double)*calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      next_PSIF_tARBS,&next_PSIF_tARBS);

    C_DGEMM('N','T',calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB,
      calc_info_.nrio,1.0,&(C_p_RB[0][0]),calc_info_.nrio,&(B_p_aS[0][0]),
      calc_info_.nrio,0.0,&(vRBaS[0][0]),calc_info_.nvirB);

    h4 -= C_DDOT(calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      &(vRBaS[0][0]),1,&(taRBS[0][0]),1);

    if (params_.logfile) {
      stop = time(NULL);
      fprintf(params_.logfilename,"finished %14ld seconds\n",stop-start);
      fflush(params_.logfilename);
    }
  }

  free_block(B_p_BB);
  free_block(B_p_aS);
  free_block(C_p_RB);
  free_block(vRBaS);
  free_block(taRBS);

  fprintf(outfile,"H4*         Energy = %18.12lf  H\n",h4);
  fflush(outfile);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  return(h4);
}

double SAPT0::exch_disp_6()
{
  time_t start;
  time_t stop;
  double q9 = 0.0;

  double **B_p_BB = get_BB_ints(1);
  double **C_p_RB = block_matrix(calc_info_.nvirA*calc_info_.noccB,calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.nvirA,calc_info_.noccB*calc_info_.nrio,calc_info_.noccB,
    1.0,&(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
    &(B_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio,0.0,&(C_p_RB[0][0]),
    calc_info_.noccB*calc_info_.nrio);

  free_block(B_p_BB);

  double **B_p_AA = get_AA_ints(1);
  double **C_p_aS = block_matrix(calc_info_.nvirB,calc_info_.nrio);
  double **vRBaS = block_matrix(calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB);
  double **taRBS = block_matrix(calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  psio_address next_PSIF_tARBS = PSIO_ZERO;

  for (int a=0; a<calc_info_.noccA; a++) {
  
    if (params_.logfile) {
      fprintf(params_.logfilename,"    Starting Occupied A %3d ... ",a);
      start = time(NULL);
      fflush(params_.logfilename);
    }

    C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.noccA,
      1.0,&(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
      &(B_p_AA[a*calc_info_.noccA][0]),calc_info_.nrio,0.0,&(C_p_aS[0][0]),
      calc_info_.nrio);

    psio_read(PSIF_SAPT_AMPS,"T2 ARBS Amplitudes",(char *) &(taRBS[0][0]),
      sizeof(double)*calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      next_PSIF_tARBS,&next_PSIF_tARBS);

    C_DGEMM('N','T',calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB,
      calc_info_.nrio,1.0,&(C_p_RB[0][0]),calc_info_.nrio,&(C_p_aS[0][0]),
      calc_info_.nrio,0.0,&(vRBaS[0][0]),calc_info_.nvirB);

    q9 += C_DDOT(calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      &(vRBaS[0][0]),1,&(taRBS[0][0]),1);

    if (params_.logfile) {
      stop = time(NULL);
      fprintf(params_.logfilename,"finished %14ld seconds\n",stop-start);
      fflush(params_.logfilename);
    }
  }

  free_block(B_p_AA);
  free_block(C_p_aS);
  free_block(C_p_RB);
  free_block(vRBaS);
  free_block(taRBS);

  fprintf(outfile,"Q9          Energy = %18.12lf  H\n",q9);
  fflush(outfile);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  return(q9);
}

double SAPT0::exch_disp_7()
{
  time_t start;
  time_t stop;
  double q12 = 0.0;

  double **B_p_AR = get_AR_ints(1);
  double **C_p_RB = block_matrix(calc_info_.nvirA*calc_info_.noccB,calc_info_.nrio);

  for (int r=0; r<calc_info_.nvirA; r++) {
    C_DGEMM('T','N',calc_info_.noccB,calc_info_.nrio,calc_info_.noccA,
      1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,&(B_p_AR[r][0]),
      calc_info_.nvirA*calc_info_.nrio,0.0,&(C_p_RB[r*calc_info_.noccB][0]),
      calc_info_.nrio);
  }

  free_block(B_p_AR);

  double **B_p_BS = get_BS_ints(1);
  double **C_p_aS = block_matrix(calc_info_.nvirB,calc_info_.nrio);
  double **vRBaS = block_matrix(calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB);
  double **taRBS = block_matrix(calc_info_.nvirA,calc_info_.noccB*calc_info_.nvirB);

  psio_address next_PSIF_tARBS = PSIO_ZERO;

  for (int a=0; a<calc_info_.noccA; a++) {
  
    if (params_.logfile) {
      fprintf(params_.logfilename,"    Starting Occupied A %3d ... ",a);
      start = time(NULL);
      fflush(params_.logfilename);
    }

    C_DGEMV('t',calc_info_.noccB,calc_info_.nvirB*calc_info_.nrio,1.0,
      &(B_p_BS[0][0]),calc_info_.nvirB*calc_info_.nrio,&(calc_info_.S_AB[a][0]),1,
      0.0,&(C_p_aS[0][0]),1);

    psio_read(PSIF_SAPT_AMPS,"T2 ARBS Amplitudes",(char *) &(taRBS[0][0]),
      sizeof(double)*calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      next_PSIF_tARBS,&next_PSIF_tARBS);

    C_DGEMM('N','T',calc_info_.nvirA*calc_info_.noccB,calc_info_.nvirB,
      calc_info_.nrio,1.0,&(C_p_RB[0][0]),calc_info_.nrio,&(C_p_aS[0][0]),
      calc_info_.nrio,0.0,&(vRBaS[0][0]),calc_info_.nvirB);

    q12 += C_DDOT(calc_info_.nvirA*calc_info_.noccB*calc_info_.nvirB,
      &(vRBaS[0][0]),1,&(taRBS[0][0]),1);

    if (params_.logfile) {
      stop = time(NULL);
      fprintf(params_.logfilename,"finished %14ld seconds\n",stop-start);
      fflush(params_.logfilename);
    }
  }

  free_block(B_p_BS);
  free_block(C_p_aS);
  free_block(C_p_RB);
  free_block(vRBaS);
  free_block(taRBS);

  fprintf(outfile,"Q12         Energy = %18.12lf  H\n",q12);
  fflush(outfile);

  if (params_.logfile) {
    fprintf(params_.logfilename,"\n");
    fflush(params_.logfilename);
  }

  return(q12);
}

void SAPT0::H1(double **Y_RB)
{
  double **B_p_RB = get_RB_ints(1);

  C_DGEMV('n',calc_info_.nvirA*calc_info_.noccB,calc_info_.nrio,2.0,&(B_p_RB[0][0]),
    calc_info_.nrio,calc_info_.diagAA,1,1.0,&(Y_RB[0][0]),1);

  free_block(B_p_RB);

  double **B_p_AR = get_AR_ints(1);
  double **B_p_AB = get_AB_ints(2);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.nvirA,calc_info_.noccB,calc_info_.nrio,-1.0,
      &(B_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,
      &(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,1.0,&(Y_RB[0][0]),
      calc_info_.noccB);
  }

  free_block(B_p_AR);
  free_block(B_p_AB);
}

void SAPT0::H3(double **Y_AS)
{
  double **B_p_AS = get_AS_ints(1);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirB,calc_info_.nrio,2.0,&(B_p_AS[0][0]),
    calc_info_.nrio,calc_info_.diagBB,1,1.0,&(Y_AS[0][0]),1);

  free_block(B_p_AS);

  double **B_p_AB = get_AB_ints(1);
  double **B_p_BS = get_BS_ints(1);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirB,calc_info_.nrio,-1.0,
      &(B_p_AB[b][0]),calc_info_.noccB*calc_info_.nrio,
      &(B_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,1.0,&(Y_AS[0][0]),
      calc_info_.nvirB);
  }

  free_block(B_p_AB);
  free_block(B_p_BS);
}

void SAPT0::Q1(double **Y_AS)
{ 
  double **B_p_AA = get_AA_ints(1);
  double **C_p_BA = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);
  
  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccA*calc_info_.nrio,calc_info_.noccA,
    1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(B_p_AA[0][0]),calc_info_.noccA*calc_info_.nrio,
    0.0,&(C_p_BA[0][0]),calc_info_.noccA*calc_info_.nrio);
  
  free_block(B_p_AA);
  
  double **B_p_BS = get_BS_ints(1);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirB,calc_info_.nrio,1.0,
      &(C_p_BA[b*calc_info_.noccA][0]),calc_info_.nrio,
      &(B_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio,1.0,&(Y_AS[0][0]),
      calc_info_.nvirB);
  }

  free_block(C_p_BA);
  free_block(B_p_BS);
}

void SAPT0::Q5(double **Y_RB)
{
  double **B_p_BB = get_BB_ints(1);
  double **C_p_AB = block_matrix(calc_info_.noccA*calc_info_.noccB,calc_info_.nrio);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.noccB*calc_info_.nrio,calc_info_.noccB,
    1.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(B_p_BB[0][0]),calc_info_.noccB*calc_info_.nrio,
    0.0,&(C_p_AB[0][0]),calc_info_.noccB*calc_info_.nrio);

  free_block(B_p_BB);

  double **B_p_AR = get_AR_ints(1);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','T',calc_info_.nvirA,calc_info_.noccB,calc_info_.nrio,1.0,
      &(B_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio,
      &(C_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,1.0,&(Y_RB[0][0]),
      calc_info_.noccB);
  }

  free_block(C_p_AB);
  free_block(B_p_AR);
}

void SAPT0::Q3(double **Y_AS)
{
  double **B_p_BS = get_BS_ints(1);
  double **X_BS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.nvirB,calc_info_.nrio,1.0,&(B_p_BS[0][0]),
    calc_info_.nrio,calc_info_.diagAA,1,1.0,&(X_BS[0][0]),1);

  free_block(B_p_BS);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirB,calc_info_.noccB,
    -2.0,&(calc_info_.S_AB[0][0]),calc_info_.nmo,
    &(X_BS[0][0]),calc_info_.nvirB,
    1.0,&(Y_AS[0][0]),calc_info_.nvirB);

  free_block(X_BS);
}

void SAPT0::Q7(double **Y_RB)
{
  double **B_p_AR = get_AR_ints(1);
  double **X_AR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.nvirA,calc_info_.nrio,1.0,&(B_p_AR[0][0]),
    calc_info_.nrio,calc_info_.diagBB,1,1.0,&(X_AR[0][0]),1);

  free_block(B_p_AR);

  C_DGEMM('T','N',calc_info_.nvirA,calc_info_.noccB,calc_info_.noccA,
    -2.0,&(X_AR[0][0]),calc_info_.nvirA,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,
    1.0,&(Y_RB[0][0]),calc_info_.noccB);

  free_block(X_AR);
}

void SAPT0::Q10(double **Y_AS)
{
  double **B_p_AA = get_AA_ints(1);
  double **X_AA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMV('n',calc_info_.noccA*calc_info_.noccA,calc_info_.nrio,1.0,&(B_p_AA[0][0]),
    calc_info_.nrio,calc_info_.diagBB,1,1.0,&(X_AA[0][0]),1);

  free_block(B_p_AA);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirB,calc_info_.noccA,
    -2.0,&(X_AA[0][0]),calc_info_.noccA,
    &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
    1.0,&(Y_AS[0][0]),calc_info_.nvirB);

  free_block(X_AA);
}

void SAPT0::Q11(double **Y_RB)
{
  double **B_p_BB = get_BB_ints(1);
  double **X_BB = block_matrix(calc_info_.noccB,calc_info_.noccB);

  C_DGEMV('n',calc_info_.noccB*calc_info_.noccB,calc_info_.nrio,1.0,&(B_p_BB[0][0]),
    calc_info_.nrio,calc_info_.diagAA,1,1.0,&(X_BB[0][0]),1);

  free_block(B_p_BB);

  C_DGEMM('N','N',calc_info_.nvirA,calc_info_.noccB,calc_info_.noccB,
    -2.0,&(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
    &(X_BB[0][0]),calc_info_.noccB,
    1.0,&(Y_RB[0][0]),calc_info_.noccB);

  free_block(X_BB);
}

void SAPT0::Q6(double **C_p_BS)
{
  double **X_BS = block_matrix(calc_info_.noccB,calc_info_.nvirB);

  C_DGEMM('T','N',calc_info_.noccB,calc_info_.nvirB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][calc_info_.noccB]),
    calc_info_.nmo,0.0,&(X_BS[0][0]),calc_info_.nvirB);

  double **B_p_BB = get_BB_ints(1);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.noccB,-2.0,
      &(X_BS[0][0]),calc_info_.nvirB,&(B_p_BB[b*calc_info_.noccB][0]),
      calc_info_.nrio,1.0,&(C_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio);
  }

  free_block(X_BS);
  free_block(B_p_BB);
}

void SAPT0::Q13(double **C_p_BS)
{
  double **X_BB = block_matrix(calc_info_.noccB,calc_info_.noccB);
  
  C_DGEMM('T','N',calc_info_.noccB,calc_info_.noccB,calc_info_.noccA,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(X_BB[0][0]),calc_info_.noccB);

  double **B_p_BS = get_BS_ints(0);

  C_DGEMM('N','N',calc_info_.noccB,calc_info_.nvirB*calc_info_.nrio,calc_info_.noccB,
    -2.0,&(X_BB[0][0]),calc_info_.noccB,&(B_p_BS[0][0]),
    calc_info_.nvirB*calc_info_.nrio,1.0,&(C_p_BS[0][0]),
    calc_info_.nvirB*calc_info_.nrio);

  free_block(X_BB);
  free_block(B_p_BS);
}

void SAPT0::Q2(double **C_p_AR)
{
  double **X_AR = block_matrix(calc_info_.noccA,calc_info_.nvirA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.nvirA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[calc_info_.noccA][0]),
    calc_info_.nmo,0.0,&(X_AR[0][0]),calc_info_.nvirA);
  
  double **B_p_AA = get_AA_ints(1);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('T','N',calc_info_.nvirA,calc_info_.nrio,calc_info_.noccA,-2.0,
      &(X_AR[0][0]),calc_info_.nvirA,&(B_p_AA[a*calc_info_.noccA][0]),
      calc_info_.nrio,1.0,&(C_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio);
  } 
    
  free_block(X_AR);
  free_block(B_p_AA);
}

void SAPT0::Q14(double **C_p_AR)
{
  double **X_AA = block_matrix(calc_info_.noccA,calc_info_.noccA);

  C_DGEMM('N','T',calc_info_.noccA,calc_info_.noccA,calc_info_.noccB,1.0,
    &(calc_info_.S_AB[0][0]),calc_info_.nmo,&(calc_info_.S_AB[0][0]),
    calc_info_.nmo,0.0,&(X_AA[0][0]),calc_info_.noccA);

  double **B_p_AR = get_AR_ints(0);

  C_DGEMM('N','N',calc_info_.noccA,calc_info_.nvirA*calc_info_.nrio,calc_info_.noccA,
    -2.0,&(X_AA[0][0]),calc_info_.noccA,&(B_p_AR[0][0]),
    calc_info_.nvirA*calc_info_.nrio,1.0,&(C_p_AR[0][0]),
    calc_info_.nvirA*calc_info_.nrio);

  free_block(X_AA);
  free_block(B_p_AR);
}

void SAPT0::H2(double **C_p_BS)
{
  double **B_p_AB = get_AB_ints(2);

  for (int b=0; b<calc_info_.noccB; b++) {
    C_DGEMM('T','N',calc_info_.nvirB,calc_info_.nrio,calc_info_.noccA,2.0,
      &(calc_info_.S_AB[0][calc_info_.noccB]),calc_info_.nmo,
      &(B_p_AB[b][0]),calc_info_.noccB*calc_info_.nrio,
      1.0,&(C_p_BS[b*calc_info_.nvirB][0]),calc_info_.nrio);
  }

  free_block(B_p_AB);
}

void SAPT0::H4(double **C_p_AR)
{
  double **B_p_AB = get_AB_ints(1);

  for (int a=0; a<calc_info_.noccA; a++) {
    C_DGEMM('N','N',calc_info_.nvirA,calc_info_.nrio,calc_info_.noccB,2.0,
      &(calc_info_.S_AB[calc_info_.noccA][0]),calc_info_.nmo,
      &(B_p_AB[a*calc_info_.noccB][0]),calc_info_.nrio,1.0,
      &(C_p_AR[a*calc_info_.nvirA][0]),calc_info_.nrio);
  }

  free_block(B_p_AB);
}

}}
