/*! \file
    \ingroup CIS
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cis {

void Fij_build(void)
{
  dpdfile2 F;
  dpdbuf4 D, T2;

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&F, CC_MISC, 0, 0, 0, "FIJ");
    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_buf4_init(&T2, CC_MISC, 0, 0, 5, 0, 5, 0, "MP2 tIjAb");
    dpd_contract442(&T2, &D, &F, 0, 0, -1, 0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);
    dpd_file2_close(&F);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&F, CC_MISC, 0, 0, 0, "FIJ");

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    dpd_buf4_init(&T2, CC_MISC, 0, 0, 7, 2, 7, 0, "MP2 tIJAB");
    dpd_contract442(&T2, &D, &F, 0, 0, -1, 0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_init(&T2, CC_MISC, 0, 22, 28, 22, 28, 0, "MP2 tIjAb");
    dpd_contract442(&T2, &D, &F, 0, 0, -1, 1);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);

    dpd_file2_close(&F);

    dpd_file2_init(&F, CC_MISC, 0, 2, 2, "Fij");

    dpd_buf4_init(&D, CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    dpd_buf4_init(&T2, CC_MISC, 0, 10, 17, 12, 17, 0, "MP2 tijab");
    dpd_contract442(&T2, &D, &F, 0, 0, -1, 0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_init(&T2, CC_MISC, 0, 22, 28, 22, 28, 0, "MP2 tIjAb");
    dpd_contract442(&T2, &D, &F, 1, 1, -1, 1);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&D);

    dpd_file2_close(&F);
  }
}

}} // namespace psi::cis
