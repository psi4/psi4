/*! \file 
    \ingroup (INPUT)
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include "input.h"
#include "global.h"
#include "defines.h"

namespace psi { namespace input {

void write_scf_calc();

/*-----------------------------------------------------------------------------------------------------------------

 -----------------------------------------------------------------------------------------------------------------*/

void write_scf_to_chkpt()
{

  double *arr_double;

  chkpt_init(PSIO_OPEN_OLD);

  /*-----------------
    Update constants
   -----------------*/
  chkpt_wt_nsymhf(num_so_typs);
  chkpt_wt_iopen(iopen);
  chkpt_wt_nmo(num_mo);
  chkpt_wt_ref(ref);

  /* write the data out */
  write_scf_calc();
  
  /* Update energies */
  chkpt_wt_escf(escf);
  chkpt_wt_eref(escf);

  chkpt_close();
  return;
}

void write_scf_calc()
{
  double *zero_array;
  zero_array = init_array(num_mo*(num_mo+1)/2);
      
  if (spinrestr_ref) {
    /* SCF eigenvector */
    chkpt_wt_scf(scf_evect_so);
    
    /* SCF eigenvalues */
    chkpt_wt_evals(zero_array);
  }
  else {
    /* SCF eigenvectors */
    chkpt_wt_alpha_scf(scf_evect_so_alpha);
    chkpt_wt_beta_scf(scf_evect_so_beta);
    
    /* SCF eigenvalues */
    chkpt_wt_alpha_evals(zero_array);
    chkpt_wt_beta_evals(zero_array);
  }
  free(zero_array);

  if(scf_evect_local != NULL) chkpt_wt_local_scf(scf_evect_local);
      
  /* irrep labels for non-empty blocks */
  chkpt_wt_irr_labs(irr_labels);

  /* MOs per block */
  chkpt_wt_orbspi(orbspi);

  /* doubly-occupied MOs per block */
  chkpt_wt_clsdpi(clsdpi);

  /* singly-occupied MOs per block */
  chkpt_wt_openpi(openpi);
  
  return;
}

}} // namespace psi::input
