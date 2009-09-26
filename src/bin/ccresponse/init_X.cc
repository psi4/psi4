/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include <libpsio/psio.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace CCRESPONSE {

void denom1(dpdfile2 *X1, double omega);
void denom2(dpdbuf4 *X2, double omega);
void local_filter_T1(dpdfile2 *T1);
void local_filter_T2(dpdbuf4 *T2);

void init_X(const char *pert, int irrep, double omega)
{
  char lbl[32];
  dpdfile2 mu1, X1, FAE, FMI;
  dpdbuf4 X2, mu2;

  sprintf(lbl, "%sBAR_IA", pert);
  dpd_file2_init(&mu1, CC_OEI, irrep, 0, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
  if(!params.restart || !psio_tocscan(CC_OEI, lbl)) {
    dpd_file2_copy(&mu1, CC_OEI, lbl);
    dpd_file2_init(&X1, CC_OEI, irrep, 0, 1, lbl);
    if(params.local && local.filter_singles) local_filter_T1(&X1);
    else denom1(&X1, omega);
    dpd_file2_close(&X1);
  }
  else fprintf(outfile, "\tUsing existing %s amplitudes.\n", lbl);
  dpd_file2_close(&mu1);

  sprintf(lbl, "%sBAR_IjAb", pert);
  dpd_buf4_init(&mu2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  if(!params.restart || !psio_tocscan(CC_LR, lbl)) {
    dpd_buf4_copy(&mu2, CC_LR, lbl);
    dpd_buf4_init(&X2, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    if(params.local) local_filter_T2(&X2);
    else denom2(&X2, omega);
    dpd_buf4_close(&X2);
  }
  else fprintf(outfile, "\tUsing existing %s amplitudes.\n", lbl);
  dpd_buf4_close(&mu2);
}

}} // namespace psi::CCRESPONSE
