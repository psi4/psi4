/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_twopdm(void);
void uhf_twopdm(void);
void rhf_sf_twopdm(void);
void uhf_sf_twopdm(void);

void twopdm(void)
{
  if(params.ref == 0) rhf_sf_twopdm();
  else if(params.ref == 2) uhf_sf_twopdm();
}

void rhf_twopdm(void)
{
  dpdbuf4 T;

  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_copy(&T, CC_GAMMA, "GIjAb");
  dpd_buf4_close(&T);
}

void uhf_twopdm(void)
{

}

void rhf_sf_twopdm(void)
{
  dpdbuf4 T;

  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
  dpd_buf4_copy(&T, CC_GAMMA, "GIJAB");
  dpd_buf4_close(&T);

  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
  dpd_buf4_copy(&T, CC_GAMMA, "Gijab");
  dpd_buf4_close(&T);

  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_copy(&T, CC_GAMMA, "GIjAb");
  dpd_buf4_close(&T);
}

void uhf_sf_twopdm(void)
{

}

}} /* End namespaces */
