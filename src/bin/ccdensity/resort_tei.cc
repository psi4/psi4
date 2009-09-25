/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void distribute(void);

int file_build(dpdfile4 *File, int inputfile, double tolerance,
	       int perm_pr, int perm_qs, int perm_prqs, int keep);

void resort_tei(void)
{
  double tolerance;
  dpdfile4 A, B, C, D, E, F;

  tolerance = params.tolerance;

  distribute();

  dpd_file4_init_nocache(&A, CC_AINTS_NEW, 0, 0, 0, "A <ij|kl>");
  file_build(&A, 90, tolerance, 1, 1, 1, 0);
  dpd_file4_close(&A);

  dpd_file4_init_nocache(&B, CC_BINTS_NEW, 0, 5, 5, "B <ab|cd>");
  file_build(&B, 91, tolerance, 1, 1, 1, 0);
  dpd_file4_close(&B);

  dpd_file4_init_nocache(&C, CC_CINTS_NEW, 0, 10, 10, "C <ia|jb>");
  file_build(&C, 92, tolerance, 1, 1, 0, 0);
  dpd_file4_close(&C);

  dpd_file4_init_nocache(&D, CC_DINTS_NEW, 0, 0, 5, "D <ij|ab>");
  file_build(&D, 93, tolerance, 0, 0, 1, 0);
  dpd_file4_close(&D);

  dpd_file4_init_nocache(&E, CC_EINTS_NEW, 0, 11, 0, "E <ai|jk>");
  file_build(&E, 94, tolerance, 0, 1, 0, 0);
  dpd_file4_close(&E);

  dpd_file4_init_nocache(&F, CC_FINTS_NEW, 0, 10, 5, "F <ia|bc>");
  file_build(&F, 95, tolerance, 0, 1, 0, 0);
  dpd_file4_close(&F);

  fflush(outfile);

}

}} // namespace psi::ccdensity
