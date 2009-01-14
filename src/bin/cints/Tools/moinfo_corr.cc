/*! \file moinfo_corr.cc
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<libciomr/libciomr.h>
#include<libqt/qt.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"small_fns.h"
#include<Tools/prints.h>

namespace psi { namespace CINTS {

/*!-------------------------------------------------------------
  For correlated calculations form scf_e(vec/vals)_(u)occ, and
  compute numbers of frozen orbitals
 -------------------------------------------------------------*/
void init_moinfo_corr()
{
  int i, j;
  int irrep, size;
  int mo, mo_max;
  int pitz_offset, qts_offset;
  int errcod;
  double **mo_row;

  /*--- Read in frozen MOs ---*/
  MOInfo.frozen_docc = get_frzcpi();
  MOInfo.frozen_uocc = get_frzvpi();

  print_moinfo_corr();

  /*------------------
    Compute mo2symblk
   ------------------*/
  MOInfo.mo2symblk = init_int_array(MOInfo.num_mo);
  pitz_offset = 0;
  for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
    mo_max = MOInfo.orbspi[irrep];
    for(mo=0;mo<mo_max;mo++) {
      MOInfo.mo2symblk[mo+pitz_offset] = irrep;
    }
    pitz_offset += mo_max;
  }
  
  /*-----------------------------------------
    Isolate occupied and virtual portions of
    eigenvectors and eigenvalues, compute
    number of active orbitals
   -----------------------------------------*/
  MOInfo.nfrdocc = 0;
  MOInfo.nfruocc = 0;
  for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
    MOInfo.nfrdocc += MOInfo.frozen_docc[irrep];
    MOInfo.nfruocc += MOInfo.frozen_uocc[irrep];
  }
  MOInfo.nactdocc = MOInfo.ndocc - MOInfo.nfrdocc;
  MOInfo.nactuocc = MOInfo.nuocc - MOInfo.nfruocc;
  switch (UserOptions.reftype) {
  case rhf:
      MOInfo.scf_evec_occ[0] = block_matrix(MOInfo.ndocc,BasisSet.num_ao);
      MOInfo.scf_evec_uocc[0] = block_matrix(MOInfo.nuocc,BasisSet.num_ao);
      MOInfo.scf_evals_occ[0] = init_array(MOInfo.ndocc);
      MOInfo.scf_evals_uocc[0] = init_array(MOInfo.nuocc);
      MOInfo.mo2symblk_occ[0] = init_int_array(MOInfo.ndocc);
      MOInfo.mo2symblk_uocc[0] = init_int_array(MOInfo.nuocc);
      /*--- Frozen DOCC ---*/
      qts_offset = 0;
      if (MOInfo.nfrdocc) {
	pitz_offset = 0;
	for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
	  mo_max = MOInfo.frozen_docc[irrep];
	  for(mo=0;mo<mo_max;mo++) {
	    memcpy(MOInfo.scf_evec_occ[0][qts_offset+mo],MOInfo.scf_evec[0][pitz_offset+mo],BasisSet.num_ao*sizeof(double));
	    MOInfo.scf_evals_occ[0][qts_offset+mo] = MOInfo.scf_evals[0][pitz_offset+mo];
	    MOInfo.mo2symblk_occ[0][qts_offset+mo] = irrep;
	  }
	  pitz_offset += MOInfo.orbspi[irrep];
	  qts_offset += mo_max;
	}
      }
      /*--- Active DOCC ---*/
      pitz_offset = 0;
      for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
	pitz_offset += MOInfo.frozen_docc[irrep];
	mo_max = MOInfo.clsdpi[irrep] - MOInfo.frozen_docc[irrep];
	for(mo=0;mo<mo_max;mo++) {
	  memcpy(MOInfo.scf_evec_occ[0][qts_offset+mo],MOInfo.scf_evec[0][pitz_offset+mo],BasisSet.num_ao*sizeof(double));
	  MOInfo.scf_evals_occ[0][qts_offset+mo] = MOInfo.scf_evals[0][pitz_offset+mo];
	  MOInfo.mo2symblk_occ[0][qts_offset+mo] = irrep;
	}
	pitz_offset += MOInfo.orbspi[irrep] - MOInfo.frozen_docc[irrep];
	qts_offset += mo_max;
      }
      /*--- Active virtuals ---*/
      qts_offset = 0;
      pitz_offset = 0;
      for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
	pitz_offset += MOInfo.clsdpi[irrep];
	mo_max = MOInfo.virtpi[irrep] - MOInfo.frozen_uocc[irrep];
	for(mo=0;mo<mo_max;mo++) {
	  memcpy(MOInfo.scf_evec_uocc[0][qts_offset+mo],MOInfo.scf_evec[0][pitz_offset+mo],BasisSet.num_ao*sizeof(double));
	  MOInfo.scf_evals_uocc[0][qts_offset+mo] = MOInfo.scf_evals[0][pitz_offset+mo];
	  MOInfo.mo2symblk_uocc[0][qts_offset+mo] = irrep;
	}
	pitz_offset += MOInfo.orbspi[irrep] - MOInfo.clsdpi[irrep];
	qts_offset += mo_max;
      }
      /*--- Frozen virtuals ---*/
      if (MOInfo.nfruocc) {
	pitz_offset = 0;
	for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
	  pitz_offset += MOInfo.orbspi[irrep] - MOInfo.frozen_uocc[irrep];
	  mo_max = MOInfo.frozen_uocc[irrep];
	  for(mo=0;mo<mo_max;mo++) {
	    memcpy(MOInfo.scf_evec_uocc[0][qts_offset+mo],MOInfo.scf_evec[0][pitz_offset+mo],BasisSet.num_ao*sizeof(double));
	    MOInfo.scf_evals_uocc[0][qts_offset+mo] = MOInfo.scf_evals[0][pitz_offset+mo];
	    MOInfo.mo2symblk_uocc[0][qts_offset+mo] = irrep;
	  }
	  pitz_offset += MOInfo.frozen_uocc[irrep];
	  qts_offset += mo_max;
	}
      }
      break;

  case uhf: /*--- Not implemented yet ---*/
/*      MOInfo.nocc_alpha = MOInfo.ndocc + MOInfo.nsocc;
      MOInfo.nocc_beta = MOInfo.ndocc;
      MOInfo.nuocc_alpha = MOInfo.nuocc;
      MOInfo.nuocc_beta = MOInfo.nuocc + MOInfo.nsocc;
      MOInfo.scf_evec_occ[0] = block_matrix(MOInfo.nocc_alpha,BasisSet.num_ao);
      MOInfo.scf_evec_occ[1] = block_matrix(MOInfo.nocc_beta,BasisSet.num_ao);
      MOInfo.scf_evec_uocc[0] = block_matrix(MOInfo.nuocc_alpha,BasisSet.num_ao);
      MOInfo.scf_evec_uocc[1] = block_matrix(MOInfo.nuocc_beta,BasisSet.num_ao);
      MOInfo.scf_evals_occ[0] = init_array(MOInfo.nocc_alpha);
      MOInfo.scf_evals_occ[1] = init_array(MOInfo.nocc_beta);
      MOInfo.scf_evals_uocc[0] = init_array(MOInfo.nuocc_alpha);
      MOInfo.scf_evals_uocc[1] = init_array(MOInfo.nuocc_beta);
      MOInfo.mo2symblk_occ[0] = init_int_array(MOInfo.nocc_alpha);
      MOInfo.mo2symblk_occ[1] = init_int_array(MOInfo.nocc_beta);
      MOInfo.mo2symblk_uocc[0] = init_int_array(MOInfo.nuocc_alpha);
      MOInfo.mo2symblk_uocc[1] = init_int_array(MOInfo.nuocc_beta);*/
      /*--- Frozen DOCC ---*/
      qts_offset = 0;
      if (MOInfo.nfrdocc) {
	pitz_offset = 0;
	for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
	  mo_max = MOInfo.frozen_docc[irrep];
	  for(mo=0;mo<mo_max;mo++) {
	    memcpy(MOInfo.scf_evec_occ[0][qts_offset+mo],MOInfo.scf_evec[0][pitz_offset+mo],BasisSet.num_ao*sizeof(double));
	    memcpy(MOInfo.scf_evec_occ[1][qts_offset+mo],MOInfo.scf_evec[1][pitz_offset+mo],BasisSet.num_ao*sizeof(double));
	    MOInfo.scf_evals_occ[0][qts_offset+mo] = MOInfo.scf_evals[0][pitz_offset+mo];
	    MOInfo.scf_evals_occ[1][qts_offset+mo] = MOInfo.scf_evals[1][pitz_offset+mo];
	    MOInfo.mo2symblk_occ[0][qts_offset+mo] = irrep;
	    MOInfo.mo2symblk_occ[1][qts_offset+mo] = irrep;
	  }
	  pitz_offset += MOInfo.orbspi[irrep];
	  qts_offset += mo_max;
	}
      }
      /*--- Active DOCC ---*/
      pitz_offset = 0;
      for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
	pitz_offset += MOInfo.frozen_docc[irrep];
	mo_max = MOInfo.clsdpi[irrep] - MOInfo.frozen_docc[irrep];
	for(mo=0;mo<mo_max;mo++) {
	  memcpy(MOInfo.scf_evec_occ[0][qts_offset+mo],MOInfo.scf_evec[0][pitz_offset+mo],BasisSet.num_ao*sizeof(double));
	  MOInfo.scf_evals_occ[0][qts_offset+mo] = MOInfo.scf_evals[0][pitz_offset+mo];
	  MOInfo.mo2symblk_occ[0][qts_offset+mo] = irrep;
	}
	pitz_offset += MOInfo.orbspi[irrep] - MOInfo.frozen_docc[irrep];
	qts_offset += mo_max;
      }
      /*--- Active virtuals ---*/
      qts_offset = 0;
      pitz_offset = 0;
      for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
	pitz_offset += MOInfo.clsdpi[irrep];
	mo_max = MOInfo.virtpi[irrep] - MOInfo.frozen_uocc[irrep];
	for(mo=0;mo<mo_max;mo++) {
	  memcpy(MOInfo.scf_evec_uocc[0][qts_offset+mo],MOInfo.scf_evec[0][pitz_offset+mo],BasisSet.num_ao*sizeof(double));
	  MOInfo.scf_evals_uocc[0][qts_offset+mo] = MOInfo.scf_evals[0][pitz_offset+mo];
	  MOInfo.mo2symblk_uocc[0][qts_offset+mo] = irrep;
	}
	pitz_offset += MOInfo.orbspi[irrep] - MOInfo.clsdpi[irrep];
	qts_offset += mo_max;
      }
      /*--- Frozen virtuals ---*/
      if (MOInfo.nfruocc) {
	pitz_offset = 0;
	for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
	  pitz_offset += MOInfo.orbspi[irrep] - MOInfo.frozen_uocc[irrep];
	  mo_max = MOInfo.frozen_uocc[irrep];
	  for(mo=0;mo<mo_max;mo++) {
	    memcpy(MOInfo.scf_evec_uocc[0][qts_offset+mo],MOInfo.scf_evec[0][pitz_offset+mo],BasisSet.num_ao*sizeof(double));
	    MOInfo.scf_evals_uocc[0][qts_offset+mo] = MOInfo.scf_evals[0][pitz_offset+mo];
	    MOInfo.mo2symblk_uocc[0][qts_offset+mo] = irrep;
	  }
	  pitz_offset += MOInfo.frozen_uocc[irrep];
	  qts_offset += mo_max;
	}
      }
      break;
  }

  return;
}


void cleanup_moinfo_corr()
{
  free(MOInfo.frozen_docc);
  free(MOInfo.frozen_uocc);
  free(MOInfo.scf_evals_occ[0]);
  free(MOInfo.scf_evals_uocc[0]);
  free_block(MOInfo.scf_evec_occ[0]);
  free_block(MOInfo.scf_evec_uocc[0]);
  free(MOInfo.mo2symblk_occ[0]);
  free(MOInfo.mo2symblk_uocc[0]);
  free(MOInfo.mo2symblk);

  return;
}
}}
