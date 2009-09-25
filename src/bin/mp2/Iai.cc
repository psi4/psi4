/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <cmath>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_sf_Iai(void);
void uhf_sf_Iai(void);

void Iai(void)
{
  if(params.ref == 0) rhf_sf_Iai();
  else if(params.ref == 2) uhf_sf_Iai();
}

void rhf_sf_Iai(void)
{
  dpdfile2 D, I;
  dpdbuf4 G, Eints, Fints;

  /* I'AI <-- sum_JK <AJ||IK> (D_JK + D_KJ) + sum_jk <Aj|Ik> (D_jk + D_kj) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 0.0);
  dpd_dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);

  /* I'ai <-- sum_jk <aj||ik> (D_jk + D_kj) + sum_jk <aJ|iK> (D_JK + D_KJ) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 0.0);
  dpd_dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_dot24(&D, &Eints, &I, 0, 0, 1.0, 1.0);
  dpd_dot24(&D, &Eints, &I, 1, 0, 1.0, 1.0);
  dpd_buf4_close(&Eints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);

  /* I'AI <-- sum_BC <IC||AB> (D_BC + D_CB) + sum_bc <Ib|Ac>(D_bc + D_cb) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  dpd_dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
  dpd_dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
  dpd_dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);

  /* I'ai <-- sum_bc <ic||ab> (D_bc + D_cb) + sum_BC <iB|aC>(D_BC + D_CB) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
  dpd_dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
  dpd_dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_dot24(&D, &Fints, &I, 0, 1, 1.0, 1.0);
  dpd_dot24(&D, &Fints, &I, 1, 1, 1.0, 1.0);
  dpd_buf4_close(&Fints);
  dpd_file2_close(&D);

  dpd_file2_close(&I);

  /* I'AI <-- sum_JBC <JA||BC> G(JI,BC) + 2 sum_jbC <jA|bC> G(jI,bC) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'AI");

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 7, 2, 7, 0, "GIJAB");
  dpd_contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_sort(&G, CC_TMP0, qpsr, 0, 5, "GiJaB");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "GiJaB");
  dpd_contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_file2_close(&I);

  /* I'ai <-- sum_jbc <ja||bc> G(ji,bc) + 2 sum_JBc <Ja|Bc> G(Ji,Bc) */
  dpd_file2_init(&I, CC_OEI, 0, 1, 0, "I'ai");

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 7, 2, 7, 0, "Gijab");
  dpd_contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_contract442(&Fints, &G, &I, 1, 1, 2.0, 1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Fints);

  dpd_file2_close(&I);
}

void uhf_sf_Iai(void)
{

}

}} /* End namespaces */
