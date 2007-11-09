/*! \file mkpt2.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<cstring>
#include<cstdlib>
#include<libint/libint.h>
#include<libipv1/ip_lib.h>
#include<libciomr/libciomr.h>
#include<libchkpt/chkpt.h>
#include<libpsio/psio.h>
#include<libqt/qt.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

#include"moinfo.h"
#include"moinfo_corr.h"
#include"mkpt2_ints.h"


namespace psi { namespace CINTS {

void run_mkpt2()
{
  using namespace psi::CINTS::mkpt2;

  init_moinfo();


  int nfocc = 0;
  int nfvir = 0;
  int ndocc = 0;
  int nactv = 0;

  /* Normally this would be in moinfo_corr, but it needs tweeking based on the active space */
  /* This is all ripped from PSIMRCC's moinfo, eventually it should directly use a moinfo object */
  ip_cwk_add(":MRCC");


  int *focc = new int[Symmetry.nirreps];
  int *docc = new int[Symmetry.nirreps];
  int *actv = new int[Symmetry.nirreps];
  int *fvir = new int[Symmetry.nirreps];

  for(int i=0; i<Symmetry.nirreps; i++){
    focc[i] = docc[i] = actv[i] = fvir[i] = 0;
  }

  // For single-point geometry optimizations and frequencies
  if(chkpt_exist(chkpt_build_keyword("Current Displacement Irrep"))){
    int   disp_irrep  = chkpt_rd_disp_irrep();
    char *save_prefix = chkpt_rd_prefix();
    int nirreps_ref;

    // read symmetry info and MOs for undisplaced geometry from
    // root section of checkpoint file
    chkpt_reset_prefix();
    chkpt_commit_prefix();

    char *ptgrp_ref = chkpt_rd_sym_label();

    // Lookup irrep correlation table
    int* correlation;
    correlate(ptgrp_ref, disp_irrep, nirreps_ref, Symmetry.nirreps,correlation);

    int *focc_ref    = new int[nirreps_ref];
    int *docc_ref    = new int[nirreps_ref];
    int *actv_ref    = new int[nirreps_ref];
    int *fvir_ref    = new int[nirreps_ref];

    // build orbital information for current point group
    read_mo_space(nirreps_ref,nfocc,focc_ref,"CORR_FOCC");
    read_mo_space(nirreps_ref,ndocc,docc_ref,"CORR_DOCC");
    read_mo_space(nirreps_ref,nactv,actv_ref,"CORR_ACTV");
    read_mo_space(nirreps_ref,nfvir,fvir_ref,"CORR_FVIR");
    for (int h=0; h < nirreps_ref; h++) {
      focc[ correlation[h] ] += focc_ref[h];
      docc[ correlation[h] ] += docc_ref[h];
      actv[ correlation[h] ] += actv_ref[h];
      fvir[ correlation[h] ] += fvir_ref[h];
    }
    chkpt_set_prefix(save_prefix);
    chkpt_commit_prefix();
    free(save_prefix);
    free(ptgrp_ref);
    delete [] correlation;
    delete [] focc_ref;
    delete [] docc_ref;
    delete [] actv_ref;
    delete [] fvir_ref;
  }else{
    // For a single-point only
    read_mo_space(Symmetry.nirreps,nfocc,focc,"CORR_FOCC");
    read_mo_space(Symmetry.nirreps,ndocc,docc,"CORR_DOCC");
    read_mo_space(Symmetry.nirreps,nactv,actv,"CORR_ACTV");
    read_mo_space(Symmetry.nirreps,nfvir,fvir,"CORR_FVIR");
  }

  

  int i, j;
  int irrep, size;
  int mo, mo_max;
  int errcod;
  double **mo_row;

  /*-----------------------------------------
    Isolate occupied and virtual portions of
    eigenvectors and eigenvalues, compute
    number of active orbitals
   -----------------------------------------*/
  MOInfo.nfrdocc           = nfocc;
  MOInfo.nfruocc           = nfvir;
  MOInfo.nactdocc          = nactv + ndocc;
  MOInfo.nactuocc          = MOInfo.num_mo - ndocc - nfocc - nfvir; // Actives are considered virtuals too!!
  MOInfo.ndocc             = nactv + ndocc + nfocc;
  MOInfo.nuocc             = MOInfo.nactuocc + nfvir - nactv;
  MOInfo.scf_evec_occ[0]   = block_matrix(MOInfo.ndocc,BasisSet.num_ao);
  MOInfo.scf_evec_uocc[0]  = block_matrix(MOInfo.nuocc,BasisSet.num_ao);
  MOInfo.mo2symblk_occ[0]  = init_int_array(MOInfo.ndocc);
  MOInfo.mo2symblk_uocc[0] = init_int_array(MOInfo.nuocc);
  MOInfo.mo2symblk         = init_int_array(MOInfo.num_mo);
  MOInfo.occ_to_pitzer     = init_int_array(MOInfo.ndocc);
  MOInfo.vir_to_pitzer     = init_int_array(MOInfo.nuocc);

  int occ_count   = 0;
  int pitz_offset = 0;
  for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
    for(mo=0;mo<MOInfo.orbspi[irrep];mo++) {
      MOInfo.mo2symblk[mo+pitz_offset] = irrep;
    }
    pitz_offset += MOInfo.orbspi[irrep];
  }

  /* We never want to freeze the core here, so we can build the fock matrix */
  /*--- OCC space ---*/
  pitz_offset = 0;
  for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
    for(mo=0;mo<focc[irrep]+docc[irrep]+actv[irrep];mo++) {
      memcpy(MOInfo.scf_evec_occ[0][occ_count],MOInfo.scf_evec[0][pitz_offset+mo],BasisSet.num_ao*sizeof(double));
      MOInfo.mo2symblk_occ[0][occ_count] = irrep;
      MOInfo.occ_to_pitzer[occ_count] = pitz_offset+mo;
      occ_count++;
    }
    pitz_offset += MOInfo.orbspi[irrep];
  }
  /*--- VIR space ---*/
  pitz_offset = 0;
  int vir_count   = 0;
  for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
    pitz_offset += focc[irrep] + docc[irrep] + actv[irrep];
    for(mo=0;mo<MOInfo.orbspi[irrep] - focc[irrep] - docc[irrep] - actv[irrep];mo++) {
      memcpy(MOInfo.scf_evec_uocc[0][vir_count],MOInfo.scf_evec[0][pitz_offset+mo],BasisSet.num_ao*sizeof(double));
      MOInfo.mo2symblk_uocc[0][vir_count] = irrep;
      MOInfo.vir_to_pitzer[vir_count] = pitz_offset+mo;
      vir_count++;
    }
    pitz_offset += MOInfo.orbspi[irrep] - focc[irrep] - docc[irrep] - actv[irrep];
  }

  switch(UserOptions.reftype) {
  case rhf:
      mkpt2_ints();
      break;
  case rohf:
      mkpt2_ints();
      break;
  case twocon:
      mkpt2_ints();
      break;
  default:
      throw std::domain_error("MkPT2 integrals with specified REFERENCE not implemented");
  }
  free(MOInfo.occ_to_pitzer);
  free(MOInfo.vir_to_pitzer);
  cleanup_moinfo_corr();
  cleanup_moinfo();

  return;
}

};};
