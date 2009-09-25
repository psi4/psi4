/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void Lmag(int L_irr)
{
  dpdfile2 R1, L1, LIA, Lia, RIA, Ria;
  dpdbuf4 R2, L2, LIJAB, Lijab, LIjAb, RIJAB, Rijab, RIjAb;
  double norm;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_file2_init(&LIA, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_init(&Lia, CC_LAMBDA, L_irr, 0, 1, "New Lia");
    dpd_buf4_init(&LIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&Lijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_init(&LIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    norm = dpd_file2_dot_self(&LIA);
    norm += dpd_file2_dot_self(&Lia);
    norm += dpd_buf4_dot_self(&LIJAB);
    norm += dpd_buf4_dot_self(&Lijab);
    norm += dpd_buf4_dot_self(&LIjAb);
    fprintf(outfile,"size of L <L|L>     %15.10lf\n",norm);

    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_close(&Lijab);
    dpd_buf4_close(&LIjAb);
  }
}

}} // namespace psi::cclambda
