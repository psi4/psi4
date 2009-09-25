/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_sf_Iia(void);
void uhf_sf_Iia(void);

void Iia(void)
{
  if(params.ref == 0) rhf_sf_Iia();
  else if(params.ref == 2) uhf_sf_Iia();
}

void rhf_sf_Iia(void)
{
  dpdfile2 I;
  dpdbuf4 G, Eints;

  /* I'IA <-- sum_BJK <JK||IB> G(JK,AB) + 2 sum_bJk <Jk|Ib> G(Jk,Ab) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'IA");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 5, 2, 7, 0, "GIJAB");
  dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);

  /* I'ia <-- sum_bjk <jk||ib> G(jk,ab) + 2 sum_BjK <Kj|Bi> G(Kj,Ba) */
  dpd_file2_init(&I, CC_OEI, 0, 0, 1, "I'ia");

  dpd_buf4_init(&Eints, CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&G, CC_GAMMA, 0, 2, 5, 2, 7, 0, "Gijab");
  dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 0.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_sort(&G, CC_TMP0, qpsr, 0, 5, "GjIbA");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "GjIbA");
  dpd_contract442(&Eints, &G, &I, 2, 2, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Eints);

  dpd_file2_close(&I);
}

void uhf_sf_Iia(void)
{

}

}} /* End namespaces */
