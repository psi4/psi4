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

void SAPT0::ind20()
{
  double indA_B, indB_A;

  indA_B = 2.0*C_DDOT(calc_info_.noccA*calc_info_.nvirA,calc_info_.CHFA[0],1,
    calc_info_.WBAR[0],1);
  indB_A = 2.0*C_DDOT(calc_info_.noccB*calc_info_.nvirB,calc_info_.CHFB[0],1,
    calc_info_.WABS[0],1);

  results_.ind20 = indA_B + indB_A;

  if (params_.print) {
    fprintf(outfile,"Ind,resp A<-B Energy = %18.12lf  H\n",indA_B);
    fprintf(outfile,"Ind,resp B<-A Energy = %18.12lf  H\n",indB_A);
    fprintf(outfile,"Ind20,resp    Energy = %18.12lf  H\n\n",results_.ind20);
    fflush(outfile);
  }
}

}}

