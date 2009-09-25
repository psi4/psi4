/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void hbar_norms() {
  double tval;
  dpdfile2 FAE, Fae, FMI, Fmi, FME, Fme;
  dpdbuf4 WMBIJ, Wmbij, WMbIj, WmBiJ, W;

  fprintf(outfile,"\n");

  if ((params.eom_ref == 0) || (params.eom_ref == 1)) {
  dpd_file2_init(&FAE, CC_OEI, H_IRR, 1, 1, "FAE");
  dpd_file2_init(&Fae, CC_OEI, H_IRR, 1, 1, "Fae");
  tval = dpd_file2_dot_self(&FAE);
  tval += dpd_file2_dot_self(&Fae);
  dpd_file2_close(&Fae);
  dpd_file2_close(&FAE);
  fprintf(outfile,"Fae   dot Fae   total %15.10lf\n", tval);

  dpd_file2_init(&FMI, CC_OEI, H_IRR, 0, 0, "FMI");
  dpd_file2_init(&Fmi, CC_OEI, H_IRR, 0, 0, "Fmi");
  tval = dpd_file2_dot_self(&FMI);
  tval += dpd_file2_dot_self(&Fmi); 
  dpd_file2_close(&Fmi);
  dpd_file2_close(&FMI);
  fprintf(outfile,"Fmi   dot Fmi   total %15.10lf\n", tval);

  dpd_file2_init(&FME, CC_OEI, H_IRR, 0, 1, "FME");
  dpd_file2_init(&Fme, CC_OEI, H_IRR, 0, 1, "Fme");
  tval = dpd_file2_dot_self(&FME);
  tval += dpd_file2_dot_self(&Fme);
  dpd_file2_close(&Fme);
  dpd_file2_close(&FME);
  fprintf(outfile,"Fme   dot Fme   total %15.10lf\n", tval);


  dpd_buf4_init(&WMBIJ, CC_HBAR, H_IRR, 10, 2, 10, 2, 0, "WMBIJ");
  tval = 2 * dpd_buf4_dot_self(&WMBIJ);
  dpd_buf4_close(&WMBIJ);
  fprintf(outfile,"WMBIJ dot WMBIJ total %15.10lf\n", tval);

  dpd_buf4_init(&Wmbij, CC_HBAR, H_IRR, 10, 2, 10, 2, 0, "Wmbij");
  tval = 2 * dpd_buf4_dot_self(&Wmbij);
  dpd_buf4_close(&Wmbij);
  fprintf(outfile,"Wmbij dot Wmbij total %15.10lf\n", tval);

  dpd_buf4_init(&WMbIj, CC_HBAR, H_IRR, 10, 0, 10, 0, 0, "WMbIj");
  tval = dpd_buf4_dot_self(&WMbIj);
  dpd_buf4_close(&WMbIj);
  fprintf(outfile,"WMbIj dot WMbIj total %15.10lf\n", tval);

  dpd_buf4_init(&WmBiJ, CC_HBAR, H_IRR, 11, 0, 11, 0, 0, "WmBiJ (Bm,Ji)");
  tval = dpd_buf4_dot_self(&WmBiJ);
  dpd_buf4_close(&WmBiJ);
  fprintf(outfile,"WmBiJ dot WmBiJ total %15.10lf\n", tval);

	  if (params.full_matrix) {
      dpd_file2_init(&FME, CC_OEI, H_IRR, 0, 1, "FAI residual");
      tval = dpd_file2_dot_self(&FME);
			dpd_file2_close(&FME);
      fprintf(outfile,"FAI residual dot FAI residual %15.10lf\n", tval);
	  }
  }
  
  else if (params.eom_ref == 2) {

    dpd_file2_init(&FAE, CC_OEI, H_IRR, 1, 1, "FAE");
    dpd_file2_init(&Fae, CC_OEI, H_IRR, 3, 3, "Fae");
    tval = dpd_file2_dot_self(&FAE);
    tval += dpd_file2_dot_self(&Fae);
    dpd_file2_close(&Fae);
    dpd_file2_close(&FAE);
    fprintf(outfile,"Fae   dot Fae   total %15.10lf\n", tval);

    dpd_file2_init(&FMI, CC_OEI, H_IRR, 0, 0, "FMI");
    dpd_file2_init(&Fmi, CC_OEI, H_IRR, 2, 2, "Fmi");
    tval = dpd_file2_dot_self(&FMI);
    tval += dpd_file2_dot_self(&Fmi); 
    dpd_file2_close(&Fmi);
    dpd_file2_close(&FMI);
    fprintf(outfile,"Fmi   dot Fmi   total %15.10lf\n", tval);

    dpd_buf4_init(&W, CC_HBAR, H_IRR, 20, 30, 20, 30, 0, "WmBeJ (JB,me)"); /* (me,JB) */
    tval = dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 30, 20, 30, 20, 0, "WMbEj (jb,ME)"); /* (ME,jb) */
    tval += dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    fprintf(outfile,"WmBeJ and WMbEj dots %15.10lf\n",tval);


    /*
    dpd_file2_init(&FME, CC_OEI, H_IRR, 0, 1, "FME");
    dpd_file2_init(&Fme, CC_OEI, H_IRR, 0, 1, "Fme");
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 2, 21, 2, 21, 0, "WMNIE");
    tval = 2.0 * dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    fprintf(outfile,"WMNIE dot WMNIE total %15.10lf\n", tval);

    dpd_buf4_init(&W, CC_HBAR, H_IRR, 12, 31, 12, 31, 0, "Wmnie");
    tval += 2.0 * dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    fprintf(outfile,"Wmnie dot Wmnie total %15.10lf\n", tval);

    dpd_buf4_init(&W, CC_HBAR, H_IRR, 22, 25, 22, 25, 0, "WMnIe");
    tval += dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    fprintf(outfile,"WMnIe dot WMnIe total %15.10lf\n", tval);

    dpd_buf4_init(&W, CC_HBAR, H_IRR, 23, 26, 23, 26, 0, "WmNiE");
    tval += dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    fprintf(outfile,"WmNiE dot WmNiE total %15.10lf\n", tval);
    */
  }
  return;
}

}} // namespace psi::cceom
