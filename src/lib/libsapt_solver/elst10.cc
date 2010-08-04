/* This function calculates the Elst10 energy */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
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

void SAPT0::elst10()
{
  double energy;

  if (params_.logfile) {
    fprintf(params_.logfilename," Elst10\n");
    fprintf(params_.logfilename,"--------\n\n");
    fflush(params_.logfilename);
  }

  if (params_.print)
    fprintf(outfile,"Begining Elst10 Calculation\n\n");

  double **B_p_A = get_diag_AA_ints(1);
  double **B_p_B = get_diag_BB_ints(1);

  calc_info_.diagAA = init_array(calc_info_.nrio);
  calc_info_.diagBB = init_array(calc_info_.nrio);

  for(int a=0; a<calc_info_.noccA; a++){
    C_DAXPY(calc_info_.nrio,1.0,&(B_p_A[a][0]),1,calc_info_.diagAA,1);
  }

  for(int b=0; b<calc_info_.noccB; b++){
    C_DAXPY(calc_info_.nrio,1.0,&(B_p_B[b][0]),1,calc_info_.diagBB,1);
  }

  energy = 4.0*C_DDOT(calc_info_.nrio,calc_info_.diagAA,1,calc_info_.diagBB,1);

  free_block(B_p_A);
  free_block(B_p_B);

  results_.elst10 = energy;

  if (params_.print) {
    fprintf(outfile,"E^{(10)}_{elst}     = %18.12lf mH\n\n",energy*1000.0);
    fflush(outfile);
  }
}

}}
