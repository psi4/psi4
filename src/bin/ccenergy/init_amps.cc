/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libdpd/dpd.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void local_filter_T2(dpdbuf4 *T2);

void init_amps(void)
{
  dpdfile2 tIA, tia, fIA, fia, dIA, dia;
  dpdbuf4 tIJAB, tijab, tIjAb, D, dIJAB, dijab, dIjAb;

  if(params.ref == 0) { /** RHF **/

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    if(!params.restart || !psio_tocscan(CC_OEI, "tIA"))
      dpd_file2_scm(&tIA, 0);
    else fprintf(outfile, "\tUsing old T1 amplitudes.\n");
    dpd_file2_close(&tIA);

    if(!params.restart || !psio_tocscan(CC_TAMPS, "tIjAb")) {
      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      dpd_buf4_copy(&D, CC_TAMPS, "tIjAb");
      dpd_buf4_close(&D);

      dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      if(params.local) {
	local_filter_T2(&tIjAb);
      }
      else {
	dpd_buf4_init(&dIjAb, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
	dpd_buf4_dirprd(&dIjAb, &tIjAb);
	dpd_buf4_close(&dIjAb);
      }
      dpd_buf4_close(&tIjAb);
    }
    else fprintf(outfile, "\tUsing old T2 amplitudes.\n\n");
  }
  else if(params.ref == 1) { /*** ROHF ***/
    if(!params.restart || !psio_tocscan(CC_OEI, "tIA") ||
       !psio_tocscan(CC_OEI, "tia")) {

      dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
      dpd_file2_copy(&fIA, CC_OEI, "tIA");
      dpd_file2_close(&fIA);
      dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
      /*  dpd_oe_scm(&tIA, 0);  */
      dpd_file2_close(&tIA);
  
      dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
      dpd_file2_copy(&fia, CC_OEI, "tia");
      dpd_file2_close(&fia);
      dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
      /*  dpd_oe_scm(&tia, 0); */
      dpd_file2_close(&tia);

      dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
      dpd_file2_dirprd(&dIA, &tIA);
      dpd_file2_close(&tIA);
      dpd_file2_close(&dIA);

      dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
      dpd_file2_init(&dia, CC_OEI, 0, 0, 1, "dia");
      dpd_file2_dirprd(&dia, &tia);
      dpd_file2_close(&tia);
      dpd_file2_close(&dia);
    }
    else fprintf(outfile, "\tUsing old T1 amplitudes.\n");

    if(!params.restart || !psio_tocscan(CC_TAMPS, "tIjAb") || 
       !psio_tocscan(CC_TAMPS, "tIJAB") || !psio_tocscan(CC_TAMPS, "tijab")) {

      dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
      dpd_buf4_copy(&D, CC_TAMPS, "tIJAB");
      dpd_buf4_copy(&D, CC_TAMPS, "tijab");
      dpd_buf4_close(&D);

      dpd_buf4_init(&dIJAB, CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
      dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
      dpd_buf4_dirprd(&dIJAB, &tIJAB);
      dpd_buf4_close(&tIJAB);
      dpd_buf4_close(&dIJAB);

      dpd_buf4_init(&dijab, CC_DENOM, 0, 1, 6, 1, 6, 0, "dijab");
      dpd_buf4_init(&tijab, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
      dpd_buf4_dirprd(&dijab, &tijab);
      dpd_buf4_close(&tijab);
      dpd_buf4_close(&dijab);

      dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      dpd_buf4_copy(&D, CC_TAMPS, "tIjAb");
      dpd_buf4_close(&D);
  
      dpd_buf4_init(&dIjAb, CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
      dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      dpd_buf4_dirprd(&dIjAb, &tIjAb);
      dpd_buf4_close(&tIjAb);
      dpd_buf4_close(&dIjAb);
    }
    else 
      fprintf(outfile, "\tUsing old T2 amplitudes.\n");
  }
  else if(params.ref == 2) { /*** UHF ***/

    if(!params.restart || !psio_tocscan(CC_OEI, "tIA") ||
       !psio_tocscan(CC_OEI, "tia")) {

      dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
      dpd_file2_copy(&fIA, CC_OEI, "tIA");
      dpd_file2_close(&fIA);
      dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
      /*    dpd_file2_scm(&tIA, 0); */
      dpd_file2_close(&tIA);

      dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
      dpd_file2_copy(&fia, CC_OEI, "tia");
      dpd_file2_close(&fia);
      dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
      /*    dpd_file2_scm(&tia, 0); */
      dpd_file2_close(&tia);

      dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
      dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
      dpd_file2_dirprd(&dIA, &tIA);
      dpd_file2_close(&tIA);
      dpd_file2_close(&dIA);

      dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
      dpd_file2_init(&dia, CC_OEI, 0, 2, 3, "dia");
      dpd_file2_dirprd(&dia, &tia);
      dpd_file2_close(&tia);
      dpd_file2_close(&dia);
    }
    else fprintf(outfile, "\tUsing old T1 amplitudes.\n");

    if(!params.restart || !psio_tocscan(CC_TAMPS, "tIjAb") || 
       !psio_tocscan(CC_TAMPS, "tIJAB") || !psio_tocscan(CC_TAMPS, "tijab")) {

      dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
      dpd_buf4_copy(&D, CC_TAMPS, "tIJAB");
      dpd_buf4_close(&D);
      dpd_buf4_init(&dIJAB, CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
      dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
      dpd_buf4_dirprd(&dIJAB, &tIJAB);
      dpd_buf4_close(&tIJAB);
      dpd_buf4_close(&dIJAB);

      dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
      dpd_buf4_copy(&D, CC_TAMPS, "tijab");
      dpd_buf4_close(&D);
      dpd_buf4_init(&dIJAB, CC_DENOM, 0, 11, 16, 11, 16, 0, "dijab");
      dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
      dpd_buf4_dirprd(&dIJAB, &tIJAB);
      dpd_buf4_close(&tIJAB);
      dpd_buf4_close(&dIJAB);

      dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
      dpd_buf4_copy(&D, CC_TAMPS, "tIjAb");
      dpd_buf4_close(&D);
      dpd_buf4_init(&dIJAB, CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
      dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
      dpd_buf4_dirprd(&dIJAB, &tIJAB);
      dpd_buf4_close(&tIJAB);
      dpd_buf4_close(&dIJAB);
    }
    else 
      fprintf(outfile, "\tUsing old T2 amplitudes.\n");
  }
  else {  /*** RHF/ROHF ***/

  }
}
}} // namespace psi::ccenergy
