/*! \file
    \ingroup CINTS
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include <cstring>
#include <cstdlib>
#include<libint/libint.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include<libciomr/libciomr.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>

#include"moinfo.h"
#include"moinfo_corr.h"
#include"rmp2r12_energy.h"

namespace psi { namespace cints {

void mp2r12()
{
  int i,pitz_offset,mo_max,irrep,mo,qts_offset;
  double **scf_evec_so;

  switch(UserOptions.reftype) {
  case rhf:
      init_moinfo();
      init_moinfo_corr(); 
      rmp2r12_energy();
      break;

  case uhf:
      /* This should live in Tools/moinfo_corr but I don't want to interfere with Daniel's CC code */
      MOInfo.num_mo = chkpt_rd_nmo();
      MOInfo.orbspi = chkpt_rd_orbspi();
      MOInfo.clsdpi = chkpt_rd_clsdpi();
      MOInfo.openpi = chkpt_rd_openpi();
      MOInfo.frozen_docc = get_frzcpi();
      MOInfo.frozen_uocc = get_frzvpi();
      MOInfo.virtpi_alpha = init_int_array(Symmetry.nirreps);
      MOInfo.virtpi_beta  = init_int_array(Symmetry.nirreps);
      MOInfo.mo2symblk = init_int_array(MOInfo.num_mo);
 
      MOInfo.ndocc = 0;
      MOInfo.nsocc = 0;
      MOInfo.nfrdocc = 0;
      MOInfo.nfruocc = 0;
      for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
        MOInfo.virtpi_alpha[irrep] = MOInfo.orbspi[irrep] - MOInfo.clsdpi[irrep] - MOInfo.openpi[irrep];
        MOInfo.virtpi_beta[irrep]  = MOInfo.orbspi[irrep] - MOInfo.clsdpi[irrep];
        MOInfo.nfrdocc += MOInfo.frozen_docc[irrep];
        MOInfo.nfruocc += MOInfo.frozen_uocc[irrep];
        MOInfo.ndocc   += MOInfo.clsdpi[irrep];
        MOInfo.nsocc   += MOInfo.openpi[irrep];
      }
      MOInfo.alpha_act_occ = MOInfo.ndocc + MOInfo.nsocc - MOInfo.nfrdocc;
      MOInfo.beta_act_occ  = MOInfo.ndocc - MOInfo.nfrdocc;
      MOInfo.alpha_occ     = MOInfo.ndocc + MOInfo.nsocc;
      MOInfo.beta_occ      = MOInfo.ndocc;

     /* Create the full SCF eigenvectors in the AO basis*/
      /* The Alpha spin case */
      scf_evec_so = chkpt_rd_alpha_scf();
      MOInfo.scf_evec_alpha = block_matrix(MOInfo.num_mo,BasisSet.num_ao);
      mmult(Symmetry.usotao,1,scf_evec_so,0,MOInfo.scf_evec_alpha,1,BasisSet.num_ao,Symmetry.num_so,MOInfo.num_mo,0);
      /* The beta spin case */
      scf_evec_so = chkpt_rd_beta_scf();
      MOInfo.scf_evec_beta = block_matrix(MOInfo.num_mo,BasisSet.num_ao);
      mmult(Symmetry.usotao,1,scf_evec_so,0,MOInfo.scf_evec_beta,1,BasisSet.num_ao,Symmetry.num_so,MOInfo.num_mo,0);
      free_block(scf_evec_so);

      pitz_offset = 0;
      for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
        mo_max = MOInfo.orbspi[irrep];
        for(mo=0;mo<mo_max;mo++) {
          MOInfo.mo2symblk[mo+pitz_offset] = irrep;
        }
        pitz_offset += mo_max;
      }


      /* Generate the arrays to map the MO numbers to the symmetry block
         for alpha and beta occupied MOs and re-order the SCF eigenvectors*/
      MOInfo.scf_evec_occ_alpha  = block_matrix(MOInfo.alpha_occ,BasisSet.num_ao);
      MOInfo.scf_evec_occ_beta   = block_matrix(MOInfo.beta_occ,BasisSet.num_ao);
      MOInfo.mo2symblk_occ_alpha = init_int_array(MOInfo.alpha_occ);
      MOInfo.mo2symblk_occ_beta  = init_int_array(MOInfo.beta_occ);
      /*--- Frozen DOCC for alpha---*/
      qts_offset = 0;
      if (MOInfo.nfrdocc) {
        pitz_offset = 0;
        for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
          mo_max = MOInfo.frozen_docc[irrep];
          for(mo=0;mo<mo_max;mo++) {
            memcpy(MOInfo.scf_evec_occ_alpha[qts_offset+mo],MOInfo.scf_evec_alpha[pitz_offset+mo],BasisSet.num_ao*sizeof(double));
            MOInfo.mo2symblk_occ_alpha[qts_offset+mo] = irrep;
          }
          qts_offset += mo_max;
          pitz_offset += MOInfo.orbspi[irrep];
        }
      }
      /*--- Active DOCC + SOCC for alpha array---*/
      pitz_offset = 0;
      for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
        pitz_offset += MOInfo.frozen_docc[irrep];
        mo_max = MOInfo.clsdpi[irrep] + MOInfo.openpi[irrep] - MOInfo.frozen_docc[irrep];
        for(mo=0;mo<mo_max;mo++) {
          memcpy(MOInfo.scf_evec_occ_alpha[qts_offset+mo],MOInfo.scf_evec_alpha[pitz_offset+mo],BasisSet.num_ao*sizeof(double));
          MOInfo.mo2symblk_occ_alpha[qts_offset+mo] = irrep;
        }
        pitz_offset += MOInfo.orbspi[irrep] - MOInfo.frozen_docc[irrep];
        qts_offset += mo_max;
      }
      
      /*--- Frozen DOCC for beta---*/
      qts_offset = 0;
      if (MOInfo.nfrdocc) {
        pitz_offset = 0;
        for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
          mo_max = MOInfo.frozen_docc[irrep];
          for(mo=0;mo<mo_max;mo++) {
            memcpy(MOInfo.scf_evec_occ_beta[qts_offset+mo],MOInfo.scf_evec_beta[pitz_offset+mo],BasisSet.num_ao*sizeof(double));
            MOInfo.mo2symblk_occ_beta[qts_offset+mo]  = irrep;
          }
          pitz_offset += MOInfo.orbspi[irrep];
          qts_offset  += mo_max;
        }
      }
      
      /*--- Active DOCC for beta array---*/
      pitz_offset = 0;
      for(irrep=0;irrep<Symmetry.nirreps;irrep++) {
        mo_max = MOInfo.clsdpi[irrep] - MOInfo.frozen_docc[irrep];
        pitz_offset += MOInfo.frozen_docc[irrep];
        for(mo=0;mo<mo_max;mo++) {
          MOInfo.mo2symblk_occ_beta[qts_offset+mo] = irrep;
          memcpy(MOInfo.scf_evec_occ_beta[qts_offset+mo],MOInfo.scf_evec_beta[pitz_offset+mo],BasisSet.num_ao*sizeof(double));
        }
        qts_offset += mo_max;	
        pitz_offset += MOInfo.orbspi[irrep] - MOInfo.frozen_docc[irrep];
      }


      ump2r12_aa::ump2r12_energy_aa();
      ump2r12_bb::ump2r12_energy_bb();

      free(MOInfo.virtpi_alpha);
      free(MOInfo.virtpi_beta);
      free_block(MOInfo.scf_evec_alpha);
      free_block(MOInfo.scf_evec_beta);
      free_block(MOInfo.scf_evec_occ_alpha);
      free_block(MOInfo.scf_evec_occ_beta);
      free(MOInfo.mo2symblk_occ_alpha);
      free(MOInfo.mo2symblk_occ_beta);
      break;
  default:
      throw std::domain_error("MP2-R12/A energy with specified REFERENCE not implemented");
  }
  cleanup_moinfo_corr();
  cleanup_moinfo();

  return;
}
};};
