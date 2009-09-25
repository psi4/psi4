/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

#define FIRST_TMP 200

void distribute_rhf(int filenum, int first_tmp, double tolerance, int keep_input);
void distribute_uhf(const char *spin, int filenum, int first_tmp, double tolerance, int keep_input);

int file_build(dpdfile4 *File, int inputfile, double tolerance,
	       int perm_pr, int perm_qs, int perm_prqs, int keep);
int file_build_multipass(dpdfile4 *File, int inputfile, double tolerance,
	       int perm_pr, int perm_qs, int perm_prqs, int keep);

int build_abcd_packed(int, double, int);

void sort_tei(void)
{
  int keep;
  double tolerance;
  dpdfile4 A, B, C, D, E, F;

  keep = params.keep_TEIFile;
  tolerance = params.tolerance;

  if(params.ref == 2) { /*** UHF ***/
    distribute_uhf("AA", PSIF_MO_AA_TEI, FIRST_TMP, tolerance, keep);

    fflush(outfile);

    dpd_file4_init_nocache(&A, CC_AINTS, 0, 0, 0, "A <IJ|KL>");
    file_build(&A, FIRST_TMP, tolerance, 1, 1, 1, 0);
    dpd_file4_close(&A);

    if(params.make_abcd) {
      dpd_file4_init_nocache(&B, CC_BINTS, 0, 5, 5, "B <AB|CD>");
      file_build_multipass(&B, FIRST_TMP+1, tolerance, 1, 1, 1, 0);
      dpd_file4_close(&B);
    }

    dpd_file4_init_nocache(&C, CC_CINTS, 0, 20, 20, "C <IA|JB>");
    file_build(&C, FIRST_TMP+2, tolerance, 1, 1, 0, 0);
    dpd_file4_close(&C);

    dpd_file4_init_nocache(&D, CC_DINTS, 0, 0, 5, "D <IJ|AB>");
    file_build(&D, FIRST_TMP+4, tolerance, 0, 0, 1, 0);
    dpd_file4_close(&D);

    dpd_file4_init_nocache(&E, CC_EINTS, 0, 21, 0, "E <AI|JK>");
    file_build(&E, FIRST_TMP+5, tolerance, 0, 1, 0, 0);
    dpd_file4_close(&E);

    dpd_file4_init_nocache(&F, CC_FINTS, 0, 20, 5, "F <IA|BC>");
    file_build(&F, FIRST_TMP+7, tolerance, 0, 1, 0, 0);
    dpd_file4_close(&F);

    dpd_file4_init_nocache(&F, CC_FINTS, 0, 21, 5, "F <AI|BC>");
    file_build(&F, FIRST_TMP+8, tolerance, 1, 0, 0, 0);
    dpd_file4_close(&F);

    distribute_uhf("BB", PSIF_MO_BB_TEI, FIRST_TMP, tolerance, keep);

    dpd_file4_init_nocache(&A, CC_AINTS, 0, 10, 10, "A <ij|kl>");
    file_build(&A, FIRST_TMP, tolerance, 1, 1, 1, 0);
    dpd_file4_close(&A);

    if(params.make_abcd) {
      dpd_file4_init_nocache(&B, CC_BINTS, 0, 15, 15, "B <ab|cd>");
      file_build_multipass(&B, FIRST_TMP+1, tolerance, 1, 1, 1, 0);
      dpd_file4_close(&B);
    }

    dpd_file4_init_nocache(&C, CC_CINTS, 0, 30, 30, "C <ia|jb>");
    file_build(&C, FIRST_TMP+2, tolerance, 1, 1, 0, 0);
    dpd_file4_close(&C);

    dpd_file4_init_nocache(&D, CC_DINTS, 0, 10, 15, "D <ij|ab>");
    file_build(&D, FIRST_TMP+4, tolerance, 0, 0, 1, 0);
    dpd_file4_close(&D);

    dpd_file4_init_nocache(&E, CC_EINTS, 0, 31, 10, "E <ai|jk>");
    file_build(&E, FIRST_TMP+5, tolerance, 0, 1, 0, 0);
    dpd_file4_close(&E);

    dpd_file4_init_nocache(&F, CC_FINTS, 0, 30, 15, "F <ia|bc>");
    file_build(&F, FIRST_TMP+7, tolerance, 0, 1, 0, 0);
    dpd_file4_close(&F);

    dpd_file4_init_nocache(&F, CC_FINTS, 0, 31, 15, "F <ai|bc>");
    file_build(&F, FIRST_TMP+8, tolerance, 1, 0, 0, 0);
    dpd_file4_close(&F);

    distribute_uhf("AB", PSIF_MO_AB_TEI, FIRST_TMP, tolerance, keep);

    dpd_file4_init_nocache(&A, CC_AINTS, 0, 22, 22, "A <Ij|Kl>");
    file_build(&A, FIRST_TMP, tolerance, 1, 1, 0, 0);
    dpd_file4_close(&A);

    if(params.make_abcd) {
      dpd_file4_init_nocache(&B, CC_BINTS, 0, 28, 28, "B <Ab|Cd>");
      file_build_multipass(&B, FIRST_TMP+1, tolerance, 1, 1, 0, 0);
      dpd_file4_close(&B);
    }

    dpd_file4_init_nocache(&C, CC_CINTS, 0, 24, 24, "C <Ia|Jb>");
    file_build(&C, FIRST_TMP+2, tolerance, 1, 1, 0, 0);
    dpd_file4_close(&C);

    dpd_file4_init_nocache(&C, CC_CINTS, 0, 26, 26, "C <Ai|Bj>");
    file_build(&C, FIRST_TMP+3, tolerance, 1, 1, 0, 0);
    dpd_file4_close(&C);

    dpd_file4_init_nocache(&D, CC_DINTS, 0, 22, 28, "D <Ij|Ab>");
    file_build(&D, FIRST_TMP+4, tolerance, 0, 0, 0, 0);
    dpd_file4_close(&D);

    dpd_file4_init_nocache(&E, CC_EINTS, 0, 26, 22, "E <Ai|Jk>");
    file_build(&E, FIRST_TMP+5, tolerance, 0, 1, 0, 0);
    dpd_file4_close(&E);

    dpd_file4_init_nocache(&E, CC_EINTS, 0, 22, 24, "E <Ij|Ka>");
    file_build(&E, FIRST_TMP+6, tolerance, 1, 0, 0, 0);
    dpd_file4_close(&E);

    dpd_file4_init_nocache(&F, CC_FINTS, 0, 24, 28, "F <Ia|Bc>");
    file_build(&F, FIRST_TMP+7, tolerance, 0, 1, 0, 0);
    dpd_file4_close(&F);

    dpd_file4_init_nocache(&F, CC_FINTS, 0, 25, 29, "F <aI|bC>");
    file_build(&F, FIRST_TMP+8, tolerance, 1, 0, 0, 0);
    dpd_file4_close(&F);

    dpd_file4_init_nocache(&F, CC_FINTS, 0, 28, 26, "F <Ab|Ci>");
    file_build(&F, FIRST_TMP+9, tolerance, 1, 0, 0, 0);
    dpd_file4_close(&F);

    dpd_file4_init_nocache(&F, CC_FINTS, 0, 26, 28, "F <Ai|Bc>");
    file_build(&F, FIRST_TMP+10, tolerance, 1, 0, 0, 0);
    dpd_file4_close(&F);

  }
  else { /*** RHF or ROHF ***/
    distribute_rhf(PSIF_MO_TEI, FIRST_TMP, tolerance, keep);

    dpd_file4_init_nocache(&A, CC_AINTS, 0, 0, 0, "A <ij|kl>");
    file_build(&A, FIRST_TMP, tolerance, 1, 1, 1, 0);
    dpd_file4_close(&A);

    if(params.make_abcd) {
      if(params.make_unpacked_abcd) {
	dpd_file4_init_nocache(&B, CC_BINTS, 0, 5, 5, "B <ab|cd>");
	file_build_multipass(&B, FIRST_TMP+1, tolerance, 1, 1, 1, 1);
	dpd_file4_close(&B);
      }
      if(params.ref == 0) build_abcd_packed(FIRST_TMP+1, tolerance, 0);
    }

    dpd_file4_init_nocache(&C, CC_CINTS, 0, 10, 10, "C <ia|jb>");
    file_build(&C, FIRST_TMP+2, tolerance, 1, 1, 0, 0);
    dpd_file4_close(&C);

    dpd_file4_init_nocache(&D, CC_DINTS, 0, 0, 5, "D <ij|ab>");
    file_build(&D, FIRST_TMP+3, tolerance, 0, 0, 1, 0);
    dpd_file4_close(&D);

    dpd_file4_init_nocache(&E, CC_EINTS, 0, 11, 0, "E <ai|jk>");
    file_build(&E, FIRST_TMP+4, tolerance, 0, 1, 0, 0);
    dpd_file4_close(&E);

    dpd_file4_init_nocache(&F, CC_FINTS, 0, 10, 5, "F <ia|bc>");
    file_build(&F, FIRST_TMP+5, tolerance, 0, 1, 0, 0);
    dpd_file4_close(&F);

    if(params.make_aibc) {
      dpd_file4_init_nocache(&F, CC_FINTS, 0, 11, 5, "F <ai|bc>");
      file_build(&F, FIRST_TMP+6, tolerance, 1, 0, 0, 0);
      dpd_file4_close(&F);
    }
  }

  fflush(outfile);
}

}} // namespace psi::ccsort
