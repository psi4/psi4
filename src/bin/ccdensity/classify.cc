/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void classify(int p, int q, int r, int s, double value,
	      struct iwlbuf *ABuf, struct iwlbuf *BBuf,
	      struct iwlbuf *CBuf, struct iwlbuf *DBuf,
	      struct iwlbuf *EBuf, struct iwlbuf *FBuf)
{
  int *occ, *vir, *socc;
  int *cc_occ, *cc_vir;
  int dirac=1;
  int soccs;

  /* NB that integrals involving frozen orbitals are KEPT in this code */

  occ = frozen.occ; vir = frozen.vir; socc = frozen.socc;
  cc_occ = frozen.allcc_occ; cc_vir = frozen.allcc_vir;

  soccs = socc[p] + socc[q] + socc[r] + socc[s];
 
  /* A (oo|oo) integrals */
  if((occ[p] && occ[q] && occ[r] && occ[s]))
      iwl_buf_wrt_val(ABuf, cc_occ[p], cc_occ[q], cc_occ[r], cc_occ[s],
		      value, 0, outfile, dirac);

  /* B (vv|vv) integrals */
  if((vir[p] && vir[q] && vir[r] && vir[s]))
      iwl_buf_wrt_val(BBuf, cc_vir[p], cc_vir[q], cc_vir[r], cc_vir[s],
		      value, 0, outfile, dirac);

  /* C (oo|vv) integrals */
  if(soccs > 1) {
      if((occ[p] && occ[q] && vir[r] && vir[s]))
	  iwl_buf_wrt_val(CBuf, cc_occ[p], cc_occ[q], cc_vir[r], cc_vir[s],
			  value, 0, outfile, dirac);
      if((occ[r] && occ[s] && vir[p] && vir[q]))
	  iwl_buf_wrt_val(CBuf, cc_occ[r], cc_occ[s], cc_vir[p], cc_vir[q],
			  value, 0, outfile, dirac);
    }
  else if((occ[p] && occ[q] && vir[r] && vir[s]))
      iwl_buf_wrt_val(CBuf, cc_occ[p], cc_occ[q], cc_vir[r], cc_vir[s],
		      value, 0, outfile, dirac);
  else if((occ[r] && occ[s] && vir[p] && vir[q]))
      iwl_buf_wrt_val(CBuf, cc_occ[r], cc_occ[s], cc_vir[p], cc_vir[q],
		      value, 0, outfile, dirac);

  /* D (ov|ov) integrals */
  if(soccs > 1) {
      if((occ[p] && vir[q] && occ[r] && vir[s]))
	  iwl_buf_wrt_val(DBuf, cc_occ[p], cc_vir[q], cc_occ[r], cc_vir[s],
			  value, 0, outfile, dirac);
      if((occ[q] && vir[p] && occ[r] && vir[s]))
	  iwl_buf_wrt_val(DBuf, cc_occ[q], cc_vir[p], cc_occ[r], cc_vir[s],
			  value, 0, outfile, dirac);
      if((occ[p] && vir[q] && occ[s] && vir[r]))
	  iwl_buf_wrt_val(DBuf, cc_occ[p], cc_vir[q], cc_occ[s], cc_vir[r],
			  value, 0, outfile, dirac);
      if((occ[q] && vir[p] && occ[s] && vir[r]))
	  iwl_buf_wrt_val(DBuf, cc_occ[q], cc_vir[p], cc_occ[s], cc_vir[r],
			  value, 0, outfile, dirac);
    }
  else if((occ[p] && vir[q] && occ[r] && vir[s]))
      iwl_buf_wrt_val(DBuf, cc_occ[p], cc_vir[q], cc_occ[r], cc_vir[s],
		      value, 0, outfile, dirac);
  else if((occ[q] && vir[p] && occ[r] && vir[s]))
      iwl_buf_wrt_val(DBuf, cc_occ[q], cc_vir[p], cc_occ[r], cc_vir[s],
		      value, 0, outfile, dirac);
  else if((occ[p] && vir[q] && occ[s] && vir[r]))
      iwl_buf_wrt_val(DBuf, cc_occ[p], cc_vir[q], cc_occ[s], cc_vir[r],
		      value, 0, outfile, dirac);
  else if((occ[q] && vir[p] && occ[s] && vir[r]))
      iwl_buf_wrt_val(DBuf, cc_occ[q], cc_vir[p], cc_occ[s], cc_vir[r],
		      value, 0, outfile, dirac);

  /* E (vo|oo) integrals */
  if(soccs > 1) {
      if((vir[p] && occ[q] && occ[r] && occ[s]))
	  iwl_buf_wrt_val(EBuf, cc_vir[p], cc_occ[q], cc_occ[r], cc_occ[s],
			  value, 0, outfile, dirac);
      if((vir[q] && occ[p] && occ[r] && occ[s]))
	  iwl_buf_wrt_val(EBuf, cc_vir[q], cc_occ[p], cc_occ[r], cc_occ[s],
			  value, 0, outfile, dirac);
      if((vir[r] && occ[s] && occ[p] && occ[q]))
	  iwl_buf_wrt_val(EBuf, cc_vir[r], cc_occ[s], cc_occ[p], cc_occ[q],
			  value, 0, outfile, dirac);
      if((vir[s] && occ[r] && occ[p] && occ[q]))
	  iwl_buf_wrt_val(EBuf, cc_vir[s], cc_occ[r], cc_occ[p], cc_occ[q],
			  value, 0, outfile, dirac);
    } 
  else if((vir[p] && occ[q] && occ[r] && occ[s]))
      iwl_buf_wrt_val(EBuf, cc_vir[p], cc_occ[q], cc_occ[r], cc_occ[s],
		      value, 0, outfile, dirac);
  else if((vir[p] && occ[q] && occ[s] && occ[r]))
      iwl_buf_wrt_val(EBuf, cc_vir[q], cc_occ[p], cc_occ[r], cc_occ[s],
		      value, 0, outfile, dirac);
  else if((vir[r] && occ[s] && occ[p] && occ[q]))
      iwl_buf_wrt_val(EBuf, cc_vir[r], cc_occ[s], cc_occ[p], cc_occ[q],
		      value, 0, outfile, dirac);
  else if((vir[r] && occ[s] && occ[q] && occ[p]))
      iwl_buf_wrt_val(EBuf, cc_vir[s], cc_occ[r], cc_occ[p], cc_occ[q], 
		      value, 0, outfile, dirac);

  /* F (ov|vv) integrals */
  if(soccs > 1) {
      if((occ[p] && vir[q] && vir[r] && vir[s]))
	  iwl_buf_wrt_val(FBuf, cc_occ[p], cc_vir[q], cc_vir[r], cc_vir[s],
			  value, 0, outfile, dirac);
      if((occ[q] && vir[p] && vir[r] && vir[s]))
	  iwl_buf_wrt_val(FBuf, cc_occ[q], cc_vir[p], cc_vir[r], cc_vir[s],
			  value, 0, outfile, dirac);
      if((occ[r] && vir[s] && vir[p] && vir[q]))
	  iwl_buf_wrt_val(FBuf, cc_occ[r], cc_vir[s], cc_vir[p], cc_vir[q],
			  value, 0, outfile, dirac);
      if((occ[s] && vir[r] && vir[p] && vir[q]))
	  iwl_buf_wrt_val(FBuf, cc_occ[s], cc_vir[r], cc_vir[p], cc_vir[q],
			  value, 0, outfile, dirac);
    }
  else if((occ[p] && vir[q] && vir[r] && vir[s]))
      iwl_buf_wrt_val(FBuf, cc_occ[p], cc_vir[q], cc_vir[r], cc_vir[s],
		      value, 0, outfile, dirac);
  else if((occ[q] && vir[p] && vir[r] && vir[s]))
      iwl_buf_wrt_val(FBuf, cc_occ[q], cc_vir[p], cc_vir[r], cc_vir[s],
		      value, 0, outfile, dirac);
  else if((occ[r] && vir[s] && vir[p] && vir[q]))
      iwl_buf_wrt_val(FBuf, cc_occ[r], cc_vir[s], cc_vir[p], cc_vir[q],
		      value, 0, outfile, dirac);
  else if((occ[s] && vir[r] && vir[p] && vir[q]))
      iwl_buf_wrt_val(FBuf, cc_occ[s], cc_vir[r], cc_vir[p], cc_vir[q],
		      value, 0, outfile, dirac);
}


}} // namespace psi::ccdensity
