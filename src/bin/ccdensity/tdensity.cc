/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void ltdensity_rohf(struct TD_Params S);
void ltdensity_uhf(struct TD_Params S);
void ltdensity_intermediates(struct TD_Params S);
void sort_ltd_rohf(struct TD_Params S);
void sort_ltd_uhf(struct TD_Params S);
void rtdensity(struct TD_Params S);
void sort_rtd_rohf(struct TD_Params S);
void sort_rtd_uhf(struct TD_Params S);

void tdensity(struct TD_Params S) {

  if(params.ref == 0 || params.ref == 1) {
    ltdensity_rohf(S);
    sort_ltd_rohf(S);
    rtdensity(S);
    sort_rtd_rohf(S);
  }
  else if(params.ref == 2) {
    ltdensity_uhf(S);
    sort_ltd_uhf(S);
    rtdensity(S);
    sort_rtd_uhf(S);
  }

  return;
}


}} // namespace psi::ccdensity
