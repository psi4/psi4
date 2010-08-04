/* This function calculates the Ind20,resp energy */

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
#include "sapt0.h"

namespace psi { namespace sapt {

void SAPT0::ind20resp()
{
  if (params_.logfile) {
    fprintf(params_.logfilename," Ind20\n");
    fprintf(params_.logfilename,"-------\n\n");
    fflush(params_.logfilename);
  }

  if (params_.print)
    fprintf(outfile,"Begining Ind20,resp Calculation\n\n");

  if (params_.logfile) {
    fprintf(params_.logfilename,"  Ind20 A<-B\n\n");
    fflush(params_.logfilename);
  }

  if (params_.print) {
    fprintf(outfile,"Solving CHF equations for monomer A\n\n");
    fflush(outfile);
  }

  calc_info_.WBAR = W_ints(PSIF_SAPT_AA_DF_INTS,"AR RI Integrals",
    calc_info_.diagBB,calc_info_.VBAA,calc_info_.noccA,calc_info_.nvirA);
  calc_info_.CHFA = block_matrix(calc_info_.noccA, calc_info_.nvirA);
  results_.indrA_B = CHF(PSIF_SAPT_AA_DF_INTS,"AA RI Integrals",
    "AR RI Integrals","RR RI Integrals", calc_info_.WBAR, calc_info_.CHFA, 
    calc_info_.evalsA, calc_info_.noccA, calc_info_.nvirA);

  if (params_.logfile) {
    fprintf(params_.logfilename,"  Ind20 B<-A\n\n");
    fflush(params_.logfilename);
  }

  if (params_.print) {
    fprintf(outfile,"Solving CHF equations for monomer B\n\n");
    fflush(outfile);
  }

  calc_info_.WABS = W_ints(PSIF_SAPT_BB_DF_INTS,"BS RI Integrals",
    calc_info_.diagAA,calc_info_.VABB,calc_info_.noccB,calc_info_.nvirB);
  calc_info_.CHFB = block_matrix(calc_info_.noccB, calc_info_.nvirB);
  results_.indrB_A = CHF(PSIF_SAPT_BB_DF_INTS,"BB RI Integrals",
    "BS RI Integrals","SS RI Integrals", calc_info_.WABS, calc_info_.CHFB, 
    calc_info_.evalsB, calc_info_.noccB, calc_info_.nvirB);

  results_.indr20 = results_.indrA_B + results_.indrB_A;

  if (params_.print) {
    fprintf(outfile,"Ind,resp A<-B Energy = %18.12lf  H\n",results_.indrA_B);
    fprintf(outfile,"Ind,resp B<-A Energy = %18.12lf  H\n",results_.indrB_A);
    fprintf(outfile,"Ind20,resp    Energy = %18.12lf  H\n\n",results_.indr20);
    fflush(outfile);
  }

}

}}
