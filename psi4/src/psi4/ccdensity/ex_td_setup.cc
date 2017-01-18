/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void ex_td_setup(struct TD_Params S, struct TD_Params U)
{
  dpdfile2 L1,R1;
  dpdbuf4 L2,R2;

  /* Ground State LHS */
  if((params.ref == 0) || (params.ref == 1)) {
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");
    global_dpd_->file2_copy(&L1, PSIF_CC_GLG, "LIA");
    //outfile->Printf(stdout, "\t*** LAMBDA ***\n");
    //global_dpd_->file2_print(&L1, stdout);
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "Lia 0 -1");
    global_dpd_->file2_copy(&L1, PSIF_CC_GLG, "Lia");
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 2, 7, 2, 7, 0, "LIJAB 0 -1");
    global_dpd_->buf4_copy(&L2, PSIF_CC_GLG, "LIJAB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 2, 7, 2, 7, 0, "Lijab 0 -1");
    global_dpd_->buf4_copy(&L2, PSIF_CC_GLG, "Lijab");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb 0 -1");
    global_dpd_->buf4_copy(&L2, PSIF_CC_GLG, "LIjAb");
    global_dpd_->buf4_close(&L2);
  }
  else if (params.ref == 2) {
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");
    global_dpd_->file2_copy(&L1, PSIF_CC_GLG, "LIA");
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 2, 3, "Lia 0 -1");
    global_dpd_->file2_copy(&L1, PSIF_CC_GLG, "Lia");
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 2, 7, 2, 7, 0, "LIJAB 0 -1");
    global_dpd_->buf4_copy(&L2, PSIF_CC_GLG, "LIJAB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 12, 17, 12, 17, 0, "Lijab 0 -1");
    global_dpd_->buf4_copy(&L2, PSIF_CC_GLG, "Lijab");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 22, 28, 22, 28, 0, "LIjAb 0 -1");
    global_dpd_->buf4_copy(&L2, PSIF_CC_GLG, "LIjAb");
    global_dpd_->buf4_close(&L2);
  }

  /* Excited state LHS */
  if((params.ref==0) || (params.ref==1)) {
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, S.irrep, 0, 1, S.L1A_lbl);
    global_dpd_->file2_copy(&L1, PSIF_CC_GL, "LIA");
    //outfile->Printf(stdout, "\t*** LHS ***\n");
    //global_dpd_->file2_print(&L1, stdout);
    //outfile->Printf( "\t*** LHS ***\n");
    //global_dpd_->file2_print(&L1, outfile);
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, S.irrep, 0, 1, S.L1B_lbl);
    global_dpd_->file2_copy(&L1, PSIF_CC_GL, "Lia");
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, S.irrep, 2, 7, 2, 7, 0, S.L2AA_lbl);
    global_dpd_->buf4_copy(&L2, PSIF_CC_GL, "LIJAB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, S.irrep, 2, 7, 2, 7, 0, S.L2BB_lbl);
    global_dpd_->buf4_copy(&L2, PSIF_CC_GL, "Lijab");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, S.irrep, 0, 5, 0, 5, 0, S.L2AB_lbl);
    global_dpd_->buf4_copy(&L2, PSIF_CC_GL, "LIjAb");
    global_dpd_->buf4_close(&L2);
  }
  else if (params.ref == 2) {
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, S.irrep, 0, 1, S.L1A_lbl);
    global_dpd_->file2_copy(&L1, PSIF_CC_GL, "LIA");
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, S.irrep, 2, 3, S.L1B_lbl);
    global_dpd_->file2_copy(&L1, PSIF_CC_GL, "Lia");
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, S.irrep, 2, 7, 2, 7, 0, S.L2AA_lbl);
    global_dpd_->buf4_copy(&L2, PSIF_CC_GL, "LIJAB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, S.irrep, 12, 17, 12, 17, 0, S.L2BB_lbl);
    global_dpd_->buf4_copy(&L2, PSIF_CC_GL, "Lijab");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, S.irrep, 22, 28, 22, 28, 0, S.L2AB_lbl);
    global_dpd_->buf4_copy(&L2, PSIF_CC_GL, "LIjAb");
    global_dpd_->buf4_close(&L2);
  }

  /* Excited state RHS */
  if((params.ref == 0) || (params.ref == 1)) {
    global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, U.irrep, 0, 1, U.R1A_lbl);
    global_dpd_->file2_copy(&R1, PSIF_CC_GR, "RIA");
    //outfile->Printf(stdout, "\t*** RHS ***\n");
    //global_dpd_->file2_print(&R1, stdout);
    //outfile->Printf( "\t*** RHS ***\n");
    //global_dpd_->file2_print(&R1, outfile);
    global_dpd_->file2_close(&R1);

    global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, U.irrep, 0, 1, U.R1B_lbl);
    global_dpd_->file2_copy(&R1, PSIF_CC_GR, "Ria");
    global_dpd_->file2_close(&R1);

    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, U.irrep, 2, 7, 2, 7, 0, U.R2AA_lbl);
    global_dpd_->buf4_copy(&R2, PSIF_CC_GR, "RIJAB");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, U.irrep, 2, 7, 2, 7, 0, U.R2BB_lbl);
    global_dpd_->buf4_copy(&R2, PSIF_CC_GR, "Rijab");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, U.irrep, 0, 5, 0, 5, 0, U.R2AB_lbl);
    global_dpd_->buf4_copy(&R2, PSIF_CC_GR, "RIjAb");
    global_dpd_->buf4_close(&R2);
  }
  else if (params.ref == 2) {
    global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, U.irrep, 0, 1, U.R1A_lbl);
    global_dpd_->file2_copy(&R1, PSIF_CC_GR, "RIA");
    global_dpd_->file2_close(&R1);

    global_dpd_->file2_init(&R1, PSIF_CC_RAMPS, U.irrep, 2, 3, U.R1B_lbl);
    global_dpd_->file2_copy(&R1, PSIF_CC_GR, "Ria");
    global_dpd_->file2_close(&R1);

    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, U.irrep, 2, 7, 2, 7, 0, U.R2AA_lbl);
    global_dpd_->buf4_copy(&R2, PSIF_CC_GR, "RIJAB");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, U.irrep, 12, 17, 12, 17, 0, U.R2BB_lbl);
    global_dpd_->buf4_copy(&R2, PSIF_CC_GR, "Rijab");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_RAMPS, U.irrep, 22, 28, 22, 28, 0, U.R2AB_lbl);
    global_dpd_->buf4_copy(&R2, PSIF_CC_GR, "RIjAb");
    global_dpd_->buf4_close(&R2);
  }

  /* Sort Ground state LHS */
  if((params.ref==0) || (params.ref==1)) {
    if(S.irrep == 0) {
      global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
      global_dpd_->file2_scm(&L1, S.R0);
      global_dpd_->file2_close(&L1);

      global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "Lia");
      global_dpd_->file2_scm(&L1, S.R0);
      global_dpd_->file2_close(&L1);

      global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 7, 2, 7, 0, "LIJAB");
      global_dpd_->buf4_scm(&L2, S.R0);
      global_dpd_->buf4_close(&L2);

      global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 2, 7, 2, 7, 0, "Lijab");
      global_dpd_->buf4_scm(&L2, S.R0);
      global_dpd_->buf4_close(&L2);

      global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
      global_dpd_->buf4_scm(&L2, S.R0);
      global_dpd_->buf4_close(&L2);
    }

    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, qpsr, 0, 5, "LiJaB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, prqs, 10, 10, "LIAJB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 2, 7, 0, "Lijab");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, prqs, 10, 10, "Liajb");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, prqs, 10, 10, "LIAjb");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, prqs, 10, 10, "LiaJB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, psrq, 10, 10, "LIbjA");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, rqps, 10, 10, "LjAIb");
    global_dpd_->buf4_close(&L2);
  }
  else if(params.ref == 2) {
    if(S.irrep == 0) {
      global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 0, 1, "LIA");
      global_dpd_->file2_scm(&L1, S.R0);
      global_dpd_->file2_close(&L1);

      global_dpd_->file2_init(&L1, PSIF_CC_GLG, 0, 2, 3, "Lia");
      global_dpd_->file2_scm(&L1, S.R0);
      global_dpd_->file2_close(&L1);

      global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 2, 7, 0, "LIJAB");
      global_dpd_->buf4_scm(&L2, S.R0);
      global_dpd_->buf4_close(&L2);

      global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 10, 15, 12, 17, 0, "Lijab");
      global_dpd_->buf4_scm(&L2, S.R0);
      global_dpd_->buf4_close(&L2);

      global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
      global_dpd_->buf4_scm(&L2, S.R0);
      global_dpd_->buf4_close(&L2);
    }

    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, qpsr, 23, 29, "LiJaB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, prqs, 20, 20, "LIAJB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 10, 15, 12, 17, 0, "Lijab");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, prqs, 30, 30, "Liajb");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, prqs, 20, 30, "LIAjb");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, prqs, 30, 20, "LiaJB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GLG, 0, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, psrq, 24, 27, "LIbjA");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GLG, rqps, 27, 24, "LjAIb");
    global_dpd_->buf4_close(&L2);
  }

  /* Sort Excited state LHS */
  if((params.ref==0) || (params.ref==1)) { /** RHF/ROHF **/
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, qpsr, 0, 5, "LiJaB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, pqsr, 0, 5, "LIjaB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, qprs, 0, 5, "LiJAb");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, prqs, 10, 10, "LIAJB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 2, 7, 0, "Lijab");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, prqs, 10, 10, "Liajb");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, prqs, 10, 10, "LIAjb");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, prqs, 10, 10, "LiaJB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 10, 10, 10, 10, 0, "LIAjb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, psrq, 10, 10, "LIbjA");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, rqps, 10, 10, "LjAIb");
    global_dpd_->buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, qpsr, 23, 29, "LiJaB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, pqsr, 22, 29, "LIjaB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, qprs, 23, 28, "LiJAb");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 0, 5, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, prqs, 20, 20, "LIAJB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 10, 15, 12, 17, 0, "Lijab");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, prqs, 30, 30, "Liajb");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, prqs, 20, 30, "LIAjb");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, prqs, 30, 20, "LiaJB");
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&L2, PSIF_CC_GL, S.irrep, 20, 30, 20, 30, 0, "LIAjb");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, psrq, 24, 27, "LIbjA");
    global_dpd_->buf4_sort(&L2, PSIF_CC_GL, rqps, 27, 24, "LjAIb");
    global_dpd_->buf4_close(&L2);
  }

  /* Sort Excited state RHS */
  if((params.ref==0) || (params.ref==1)) { /** RHF/ROHF **/
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, U.irrep, 0, 5, 0, 5, 0, "RIjAb");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, qpsr, 0, 5, "RiJaB");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, pqsr, 0, 5, "RIjaB");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, qprs, 0, 5, "RiJAb");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_GR, U.irrep, 0, 5, 2, 7, 0, "RIJAB");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, prqs, 10, 10, "RIAJB");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_GR, U.irrep, 0, 5, 2, 7, 0, "Rijab");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, prqs, 10, 10, "Riajb");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_GR, U.irrep, 0, 5, 0, 5, 0, "RIjAb");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, prqs, 10, 10, "RIAjb");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_GR, U.irrep, 0, 5, 0, 5, 0, "RiJaB");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, prqs, 10, 10, "RiaJB");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_GR, U.irrep, 10, 10, 10, 10, 0, "RIAjb");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, psrq, 10, 10, "RIbjA");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, rqps, 10, 10, "RjAIb");
    global_dpd_->buf4_close(&R2);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->buf4_init(&R2, PSIF_CC_GR, U.irrep, 22, 28, 22, 28, 0, "RIjAb");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, qpsr, 23, 29, "RiJaB");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, pqsr, 22, 29, "RIjaB");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, qprs, 23, 28, "RiJAb");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_GR, U.irrep, 0, 5, 2, 7, 0, "RIJAB");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, prqs, 20, 20, "RIAJB");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_GR, U.irrep, 10, 15, 12, 17, 0, "Rijab");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, prqs, 30, 30, "Riajb");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_GR, U.irrep, 22, 28, 22, 28, 0, "RIjAb");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, prqs, 20, 30, "RIAjb");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_GR, U.irrep, 23, 29, 23, 29, 0, "RiJaB");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, prqs, 30, 20, "RiaJB");
    global_dpd_->buf4_close(&R2);

    global_dpd_->buf4_init(&R2, PSIF_CC_GR, U.irrep, 20, 30, 20, 30, 0, "RIAjb");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, psrq, 24, 27, "RIbjA");
    global_dpd_->buf4_sort(&R2, PSIF_CC_GR, rqps, 27, 24, "RjAIb");
    global_dpd_->buf4_close(&R2);
  }

  return;
}

}} // namespace psi::ccdensity
