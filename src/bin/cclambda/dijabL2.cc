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

void local_filter_T1(dpdfile2 *T1);
void local_filter_T2(dpdbuf4 *T2);

void dijabL2(int L_irr)
{
  dpdbuf4 L2, newLIJAB, newLijab, newLIjAb;
  dpdbuf4 d2, dIJAB, dijab, dIjAb;

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_copy(&newLIjAb, CC_LAMBDA, "New LIjAb Increment");
    dpd_buf4_close(&newLIjAb);

    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb Increment");
    if(params.local) local_filter_T2(&newLIjAb);
    else {
      dpd_buf4_init(&dIjAb, CC_DENOM, L_irr, 0, 5, 0, 5, 0, "dIjAb");
      dpd_buf4_dirprd(&dIjAb, &newLIjAb);
      dpd_buf4_close(&dIjAb);
    }
    dpd_buf4_close(&newLIjAb);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_copy(&L2, CC_LAMBDA, "New LIjAb");
    dpd_buf4_close(&L2);
    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb Increment");
    dpd_buf4_axpy(&L2, &newLIjAb, 1);
    dpd_buf4_close(&L2);
    /*dpd_buf4_print(&newLIjAb,outfile,1);*/
    dpd_buf4_close(&newLIjAb);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 0, 5, 1, "New LIjAb");
    dpd_buf4_copy(&L2, CC_LAMBDA, "New LIJAB");
    dpd_buf4_copy(&L2, CC_LAMBDA, "New Lijab");
    dpd_buf4_close(&L2);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_buf4_init(&newLIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_copy(&newLIJAB, CC_LAMBDA, "New LIJAB Increment");
    dpd_buf4_close(&newLIJAB);

    dpd_buf4_init(&newLIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB Increment");
    dpd_buf4_init(&dIJAB, CC_DENOM, L_irr, 1, 6, 1, 6, 0, "dIJAB");
    dpd_buf4_dirprd(&dIJAB, &newLIJAB);
    dpd_buf4_close(&dIJAB);
    dpd_buf4_close(&newLIJAB);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_copy(&L2, CC_LAMBDA, "New LIJAB");
    dpd_buf4_close(&L2);
    dpd_buf4_init(&newLIJAB, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB Increment");
    dpd_buf4_axpy(&L2, &newLIJAB, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&newLIJAB);

    dpd_buf4_init(&newLijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_copy(&newLijab, CC_LAMBDA, "New Lijab Increment");
    dpd_buf4_close(&newLijab);

    dpd_buf4_init(&newLijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab Increment");
    dpd_buf4_init(&dijab, CC_DENOM, L_irr, 1, 6, 1, 6, 0, "dijab");
    dpd_buf4_dirprd(&dijab, &newLijab);
    dpd_buf4_close(&dijab);
    dpd_buf4_close(&newLijab);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    dpd_buf4_copy(&L2, CC_LAMBDA, "New Lijab");
    dpd_buf4_close(&L2);
    dpd_buf4_init(&newLijab, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab Increment");
    dpd_buf4_axpy(&L2, &newLijab, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&newLijab);

    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_copy(&newLIjAb, CC_LAMBDA, "New LIjAb Increment");
    dpd_buf4_close(&newLIjAb);

    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb Increment");
    dpd_buf4_init(&dIjAb, CC_DENOM, L_irr, 0, 5, 0, 5, 0, "dIjAb");
    dpd_buf4_dirprd(&dIjAb, &newLIjAb);
    dpd_buf4_close(&dIjAb);
    dpd_buf4_close(&newLIjAb);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_copy(&L2, CC_LAMBDA, "New LIjAb");
    dpd_buf4_close(&L2);
    dpd_buf4_init(&newLIjAb, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb Increment");
    dpd_buf4_axpy(&L2, &newLIjAb, 1);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&newLIjAb);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&d2, CC_DENOM, L_irr, 1, 6, 1, 6, 0, "dIJAB");
    dpd_buf4_dirprd(&d2, &L2);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&d2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_init(&d2, CC_DENOM, L_irr, 11, 16, 11, 16, 0, "dijab");
    dpd_buf4_dirprd(&d2, &L2);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&d2);

    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_init(&d2, CC_DENOM, L_irr, 22, 28, 22, 28, 0, "dIjAb");
    dpd_buf4_dirprd(&d2, &L2);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&d2);

  }
}


}} // namespace psi::cclambda
