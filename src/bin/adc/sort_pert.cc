
/*
 *  sort_pert.cc
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

    void sort_pert(void)
    {
      int p, q, Gp, Gq, P, Q, i;
      int irrep;
      dpdfile2 f;
      double **F;
      char lbl[32];
      
      for(i = 0; i < 3; i++) {
	if(i == 0) {
	  irrep = moinfo.irrep_x;
	  F = MU_X;
	  sprintf(lbl, "MU_%1s_IA", "X");
	}
	else if(i == 1) {
	  irrep = moinfo.irrep_y;
	  F = MU_Y;
	  sprintf(lbl, "MU_%1s_IA", "Y");
	}
	else if(i == 2) {
	  irrep = moinfo.irrep_z;
	  F = MU_Z;
	  sprintf(lbl, "MU_%1s_IA", "Z");
	}
	
	dpd_file2_init(&f, CC_OEI, irrep, 0, 1, lbl);
	dpd_file2_mat_init(&f);
	for(Gp = 0; Gp < moinfo.nirreps; Gp++) {
	  Gq = irrep ^ Gp;
	  for(p = 0; p < moinfo.occpi[Gp]; p++) {
	    P = moinfo.qt2pitzer[moinfo.qt_occ[p+moinfo.occ_off[Gp]]];
	    for(q = 0; q < moinfo.virpi[Gq]; q++) {
	      Q = moinfo.qt2pitzer[moinfo.qt_vir[q+moinfo.vir_off[Gq]]];
	      f.matrix[Gp][p][q] = F[P][Q];
	    }
	  }
	}
	dpd_file2_mat_wrt(&f);
	dpd_file2_mat_close(&f);
	dpd_file2_close(&f);
	free_block(F);
      }
      
    }
    
}}
