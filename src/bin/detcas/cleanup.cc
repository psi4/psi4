/*! \file
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
#include <cstdlib>
#include <cstdio>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "globaldefs.h"
#include "globals.h"

namespace psi { namespace detcas {

/*
** cleanup()
**
** This function frees any allocated global variables
**
*/
void cleanup(void)
{
  int i;
  
  free(CalcInfo.docc);
  free(CalcInfo.socc);
  free(CalcInfo.frozen_docc);
  free(CalcInfo.frozen_uocc);
  free(CalcInfo.rstr_docc);
  free(CalcInfo.rstr_uocc);
  free(CalcInfo.orbsym);
  free(CalcInfo.pitz2ci);
  free(CalcInfo.ci2pitz);
  free(CalcInfo.ci2relpitz);
  free(CalcInfo.first);
  free(CalcInfo.last);
  free(CalcInfo.fstact);
  free(CalcInfo.lstact);
  free(CalcInfo.active);
  free_int_matrix(CalcInfo.ras_opi);
  free_int_matrix(CalcInfo.fzc_orbs);
  free_int_matrix(CalcInfo.fzv_orbs);
  for (i=0; i<MAX_RAS_SPACES; i++) 
    free_int_matrix(CalcInfo.ras_orbs[i]);
  free(CalcInfo.ras_orbs);
  for (i=0; i<CalcInfo.nirreps; i++) 
    free(CalcInfo.labels[i]);

  for (i=0; i<CalcInfo.nirreps; i++) {
    if (CalcInfo.orbs_per_irr[i]) 
      free_block(CalcInfo.mo_coeffs[i]);
  }
  free(CalcInfo.mo_coeffs);

  free(CalcInfo.onel_ints);
  free(CalcInfo.twoel_ints);
  free_block(CalcInfo.opdm);
  free(CalcInfo.tpdm);
  free_block(CalcInfo.lag);
  free(CalcInfo.F_act);
  free(CalcInfo.mo_grad);
  if (CalcInfo.mo_hess_diag != NULL) free(CalcInfo.mo_hess_diag);
  if (CalcInfo.mo_hess != NULL) free_block(CalcInfo.mo_hess);
  free(CalcInfo.theta_cur);
  free(CalcInfo.theta_step);
  free(CalcInfo.orbs_per_irr);
}

}} // end namespace psi::detcas

