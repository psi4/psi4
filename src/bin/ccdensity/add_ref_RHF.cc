/*! \file 
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libiwl/iwl.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

    void add_ref_RHF(struct iwlbuf *OutBuf)
    {
      int i,j;
      int nfzc, nclsd, nopen;

      nfzc = moinfo.nfzc;
      nclsd = moinfo.nclsd;
      nopen = moinfo.nopen;

      /*** One-electron component ***/

      for(i=0; i < (nfzc + nclsd); i++)
	moinfo.opdm[i][i] += 2.0;

      for(i=nfzc + nclsd; i < (nfzc + nclsd + nopen); i++)
	moinfo.opdm[i][i] += 1.0;

      /*** Two-electron component ***/

      /* docc-docc */
      for(i=0; i < (nfzc + nclsd); i++) {
	iwl_buf_wrt_val(OutBuf, i, i, i, i, 1.0, 0, outfile, 0);
	for(j=0; j < i; j++) {
	  iwl_buf_wrt_val(OutBuf, i, i, j, j, 2.0, 0, outfile, 0);
	  iwl_buf_wrt_val(OutBuf, i, j, j, i,-1.0, 0, outfile, 0);
	}
      }

    }

  }} // namespace psi::ccdensity
