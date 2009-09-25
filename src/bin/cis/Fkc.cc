/*! \file
    \ingroup CIS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cis {

void Fkc_build(int irrep, int root, enum Spin spin)
{
  char lbl[32];
  dpdfile2 B, F, B_A, B_B;
  dpdbuf4 D;

  if(params.ref == 0) { /** RHF **/

    if(spin == singlet)
      sprintf(lbl, "BIA(%d)[%d] singlet", root, irrep);
    else 
      sprintf(lbl, "BIA(%d)[%d] triplet", root, irrep);

    dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);

    sprintf(lbl, "FKC(%d)[%d]", root, irrep);
    dpd_file2_init(&F, CC_MISC, irrep, 0, 1, lbl);

    if(spin == singlet) {
      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
      dpd_dot13(&B, &D, &F, 0, 0, 1, 0);
    }
    else {
      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      dpd_dot14(&B, &D, &F, 0, 0, -1, 0);
    }

    dpd_buf4_close(&D);

    dpd_file2_close(&F);
    dpd_file2_close(&B);
  }
  else if(params.ref == 2) { /** UHF **/

    sprintf(lbl, "BIA(%d)[%d]", root, irrep);
    dpd_file2_init(&B_A, CC_OEI, irrep, 0, 1, lbl);
    sprintf(lbl, "Bia(%d)[%d]", root, irrep);
    dpd_file2_init(&B_B, CC_OEI, irrep, 2, 3, lbl);

    sprintf(lbl, "FKC(%d)[%d]", root, irrep);
    dpd_file2_init(&F, CC_MISC, irrep, 0, 1, lbl);
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    dpd_dot13(&B_A, &D, &F, 0, 0, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_dot24(&B_B, &D, &F, 0, 0, 1, 1);
    dpd_buf4_close(&D);

    dpd_file2_close(&F);

    sprintf(lbl, "Fkc(%d)[%d]", root, irrep);
    dpd_file2_init(&F, CC_MISC, irrep, 2, 3, lbl);
    dpd_buf4_init(&D, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    dpd_dot13(&B_B, &D, &F, 0, 0, 1, 0);
    dpd_buf4_close(&D);
    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_dot13(&B_A, &D, &F, 0, 0, 1, 1);
    dpd_buf4_close(&D);
    dpd_file2_close(&F);

    dpd_file2_close(&B_A);
    dpd_file2_close(&B_B);
  }
}

}} // namespace psi::cis
