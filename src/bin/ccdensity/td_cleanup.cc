/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void td_cleanup(void)
{
  psio_close(CC_TMP,0);
  psio_close(EOM_TMP,0);
  psio_close(EOM_TMP0,0);
  psio_close(EOM_TMP1,0);
  psio_close(CC_GLG,0);
  psio_close(CC_GL,0);
  psio_close(CC_GR,0);

  psio_open(CC_TMP,PSIO_OPEN_NEW);
  psio_open(EOM_TMP,PSIO_OPEN_NEW);
  psio_open(EOM_TMP0,PSIO_OPEN_NEW);
  psio_open(EOM_TMP1,PSIO_OPEN_NEW);
  psio_open(CC_GLG,PSIO_OPEN_NEW);
  psio_open(CC_GL,PSIO_OPEN_NEW);
  psio_open(CC_GR,PSIO_OPEN_NEW);

  if((params.ref==0) || (params.ref==1)) {
    free_block(moinfo.ltd);
    free_block(moinfo.rtd);
  }
  else if(params.ref==2) {
    free_block(moinfo.ltd_a);
    free_block(moinfo.ltd_b);
    free_block(moinfo.rtd_a);
    free_block(moinfo.rtd_b);
  }

  return;
}

}} // namespace psi::ccdensity
