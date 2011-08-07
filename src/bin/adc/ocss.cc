
/*
 *  ocss.cc
 *  
 *
 *  Created by M.Saitow on 11/06/23.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace adc {

    double calc_ocss(char *amps, int my_irrep, int whether_d, int root, double omega)
    {
      int i, irrep;
      double f, rf, val, *mu_sq;
      dpdfile2 MU, B;
      char lbl[32];
      
      if (whether_d){
	rf    = params.rfs[my_irrep][root];
	mu_sq = params.d_mu_sqs[my_irrep][root];
      }
      else{
	rf = 1;
	mu_sq = params.mu_sqs[my_irrep][root];
	
      }
      dpd_file2_init(&B,  CC_OEI, my_irrep, 0, 1, amps);
      for(i = 0;i < 3;i++){
	if(i == 0){ 
	  irrep = moinfo.irrep_x;
	  sprintf(lbl, "MU_%1s_IA", "X");
	}
	else if(i == 1){
	  irrep = moinfo.irrep_y;
	  sprintf(lbl, "MU_%1s_IA", "Y");
	}
	else if(i == 2){
	  irrep = moinfo.irrep_z;
	  sprintf(lbl, "MU_%1s_IA", "Z");
	}
	if(irrep == my_irrep){
	  dpd_file2_init(&MU, CC_OEI, irrep, 0, 1, lbl);
	  val = dpd_file2_dot(&MU, &B);
	  mu_sq[i] = val * val / rf;
	  dpd_file2_close(&MU);
	}
      }
      f = 4 * omega * (mu_sq[0] + mu_sq[1] + mu_sq[2]) / 3; 
      dpd_file2_close(&B);
      
      return f;
    }

}}


