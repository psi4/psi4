
/*
 *  ocss.cc
 *  
 *
 *  Created by M.Saitow on 11/06/23.
 *  Copyright 2010 M.Saitow. All rights reserved.
 *
 */

#include <cstdlib>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace adc {

    void transpert(void)
    {
      int nao, nso, nmo, noei_ao;
      int alpha;
      int i, j, ij;
      double *scratch, **TMP, **X, **target;
      const char *name;
      double val;
      
      nao = moinfo.nao;
      nso = moinfo.nso;
      nmo = moinfo.nmo;
      noei_ao = nao * (nao+1) / 2;
      
      TMP = block_matrix(nao, nao);
      X = block_matrix(nao, nao);
      scratch = init_array(noei_ao);
      
      MU_X = block_matrix(nmo, nmo);
      MU_Y = block_matrix(nmo, nmo);
      MU_Z = block_matrix(nmo, nmo);
      
      for(alpha = 0; alpha < 3; alpha++) {
	
	target = block_matrix(nmo,nmo);
	if     (alpha == 0) { name = PSIF_AO_MX; MU_X = target; }
	else if(alpha == 1) { name = PSIF_AO_MY; MU_Y = target; }
	else if(alpha == 2) { name = PSIF_AO_MZ; MU_Z = target; }
	
	iwl_rdone(PSIF_OEI, name, scratch, noei_ao, 0, 0, outfile);
	for(i = 0, ij = 0; i < nao; i++)
	  for(j = 0; j <= i; j++, ij++) {
	    val = scratch[ij];
	    TMP[i][j] = val;
	    TMP[j][i] = val;
	  }
	
	C_DGEMM('n', 't', nao, nso, nao, 1, &(TMP[0][0]), nao, &(moinfo.usotao[0][0]), nao,
		0, &(X[0][0]), nao);
	C_DGEMM('n', 'n', nso, nso, nao, 1, &(moinfo.usotao[0][0]), nao, &(X[0][0]), nao,
		0, &(TMP[0][0]), nao);
	
	C_DGEMM('n', 'n', nso, nmo, nso, 1, &(TMP[0][0]), nao, &(moinfo.scf[0][0]), nmo,
		0, &(X[0][0]), nao);
	C_DGEMM('t', 'n', nmo, nmo, nso, 1, &(moinfo.scf[0][0]), nmo, &(X[0][0]), nao,
		0, &(target[0][0]), nmo);
	
	zero_arr(scratch, noei_ao);
      }
      
      free(scratch);
      free_block(TMP);
      free_block(X);
    }
	
}}


