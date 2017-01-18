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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
/* restart() collapses L vectors down to num vectors */

#include <cstdio>
#include <cmath>
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void restart(double **alpha, int L, int num, int C_irr, int ortho, double **alpha_old, int L_old, int use_alpha_old) {
  int i,I,j,h, A_OCC, B_OCC, A_VIR, B_VIR, cnt, L_tot;
  int AA_OCC, AA_VIR, BB_OCC, BB_VIR, AB_OCC, AB_VIR;
  char lbl[20];
  dpdfile2 C1, CME, Cme, SIA, Sia;
  double C10, CME0, S0, **alpha_tot;
  dpdbuf4 C2, CMNEF, Cmnef, CMnEf, SIJAB, Sijab, SIjAb;
  double dotval, norm;

#ifdef TIME_CCEOM
timer_on("RESTART");
#endif

  A_OCC = 0; A_VIR = 1;
  AA_OCC = 2; AA_VIR = 7;
  if (params.eom_ref <= 1) {
    B_OCC = 0; B_VIR = 1;
    BB_OCC = 2; BB_VIR = 7;
    AB_OCC = 0; AB_VIR = 5;
  }
  else if (params.eom_ref == 2) {
    B_OCC = 2; B_VIR = 3;
    BB_OCC = 12; BB_VIR = 17;
    AB_OCC = 22; AB_VIR = 28;
  }

  if (use_alpha_old)
    L_tot = L + L_old;
  else
    L_tot = L;
  alpha_tot = block_matrix(L,L_tot);

  if (use_alpha_old) {
    cnt = 0;
    for (j=0;j<num;++j,++cnt)
      for (i=0; i<L; ++i)
        alpha_tot[i][cnt] = alpha[i][j];
    for (j=0;j<num;++j,++cnt)
      for (i=0; i<L_old; ++i)
        alpha_tot[i][cnt] = alpha_old[i][j];
    for (j=num;j<L;++j,++cnt)
      for (i=0; i<L; ++i)
        alpha_tot[i][cnt] = alpha[i][j];
    for (j=num;j<L_old;++j,++cnt)
      for (i=0; i<L_old; ++i)
        alpha_tot[i][cnt] = alpha_old[i][j];
    num *= 2;
  }
  else {
    for (i=0; i<L; ++i)
      for (j=0;j<L;++j)
        alpha_tot[i][j] = alpha[i][j];
  }
  /*
  outfile->Printf("alpha\n");
  print_mat(alpha,L,L,outfile);
  outfile->Printf("alpha_old\n");
  print_mat(alpha_old,L_old,L_old,outfile);
  outfile->Printf("alpha_tot\n");
  print_mat(alpha_tot,L,L_tot,outfile);
  */

  /* Orthonormalize alpha[1] through alpha[num] */
  if ((ortho) || (use_alpha_old)) {
     /* Since the state of interest is usually the highest energy one, lets leave
     the highest energy R alone and orthonormalize the lower ones against it.
     It is hoped that this change will help excited-state CC3 computations where
     there are multiple states with high overlap of R's and the overlap between
     the EOM CCSD root and the EOM CC3 root is used to follow the right root.
     Previously, the first had been let alone.  -RAK 4-09 */
    for ( I=num-2; I>-1; --I) {
      for (i=num-1; i>I; --i) {
        dotval = 0.0;
        for (j=0;j<L;++j) {
          dotval += alpha_tot[j][i] * alpha_tot[j][I];
        }
        for (j=0; j<L; j++) alpha_tot[j][I] -= dotval * alpha_tot[j][i];
      }
      dotval = 0.0;
      for (j=0;j<L;++j) dotval += alpha_tot[j][I] * alpha_tot[j][I];
      norm = sqrt(dotval);
      for (j=0;j<L;++j) alpha_tot[j][I] = alpha_tot[j][I]/norm;
    }
/*
    for (I=1;I<num;++I) {
      for (i=0; i<I; i++) {
        dotval = 0.0;
        for (j=0;j<L;++j) {
          dotval += alpha_tot[j][i] * alpha_tot[j][I];
        }
        for (j=0; j<L; j++) alpha_tot[j][I] -= dotval * alpha_tot[j][i];
      }
      dotval = 0.0;
      for (j=0;j<L;++j) dotval += alpha_tot[j][I] * alpha_tot[j][I];
      norm = sqrt(dotval);
      for (j=0;j<L;++j) alpha_tot[j][I] = alpha_tot[j][I]/norm;
    }
*/
  }

  /* Form restart vectors Ci = Sum_j(alpha[j][i]*Cj) */
  for (i=0; i<num; ++i) {
    sprintf(lbl, "%s %d", "CME", L+i);
    global_dpd_->file2_init(&C1, PSIF_EOM_CME, C_irr, A_OCC, A_VIR, lbl);
    global_dpd_->file2_scm(&C1, 0.0);
		C10 = 0.0;
    for (j=0;j<L;++j) {
      sprintf(lbl, "%s %d", "CME", j);
      global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, A_OCC, A_VIR, lbl);
      global_dpd_->file2_axpy(&CME, &C1, alpha_tot[j][i], 0);
      global_dpd_->file2_close(&CME);
      if (params.full_matrix) {
		    sprintf(lbl, "%s %d", "C0", j);
			  psio_read_entry(PSIF_EOM_CME, lbl, (char *) &CME0, sizeof(double));
				C10 += alpha_tot[j][i] * CME0;
		  }
    }
    if (params.full_matrix) {
		  sprintf(lbl, "%s %d", "C0", L+i);
	    psio_write_entry(PSIF_EOM_CME, lbl, (char *) &C10, sizeof(double));
		}
    global_dpd_->file2_close(&C1);

    sprintf(lbl, "%s %d", "CMnEf", L+i);
    global_dpd_->buf4_init(&C2, PSIF_EOM_CMnEf, C_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, lbl);
    global_dpd_->buf4_scm(&C2, 0.0);
    for (j=0;j<L;++j) {
      sprintf(lbl, "%s %d", "CMnEf", j);
      global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, lbl);
      global_dpd_->buf4_axpy(&CMnEf, &C2, alpha_tot[j][i]);
      global_dpd_->buf4_close(&CMnEf);
    }
    global_dpd_->buf4_close(&C2);

    if (params.eom_ref > 0) {
      sprintf(lbl, "%s %d", "Cme", L+i);
      global_dpd_->file2_init(&C1, PSIF_EOM_Cme, C_irr, B_OCC, B_VIR, lbl);
      global_dpd_->file2_scm(&C1, 0.0);
      for (j=0;j<L;++j) {
        sprintf(lbl, "%s %d", "Cme", j);
        global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, B_OCC, B_VIR, lbl);
        global_dpd_->file2_axpy(&Cme, &C1, alpha_tot[j][i], 0);
        global_dpd_->file2_close(&Cme);
      }
      global_dpd_->file2_close(&C1);

      sprintf(lbl, "%s %d", "CMNEF", L+i);
      global_dpd_->buf4_init(&C2, PSIF_EOM_CMNEF, C_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, lbl);
      global_dpd_->buf4_scm(&C2, 0.0);
      for (j=0;j<L;++j) {
        sprintf(lbl, "%s %d", "CMNEF", j);
        global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, lbl);
        global_dpd_->buf4_axpy(&CMNEF, &C2, alpha_tot[j][i]);
        global_dpd_->buf4_close(&CMNEF);
      }
      global_dpd_->buf4_close(&C2);

      sprintf(lbl, "%s %d", "Cmnef", L+i);
      global_dpd_->buf4_init(&C2, PSIF_EOM_Cmnef, C_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, lbl);
      global_dpd_->buf4_scm(&C2, 0.0);
      for (j=0;j<L;++j) {
        sprintf(lbl, "%s %d", "Cmnef", j);
        global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, lbl);
        global_dpd_->buf4_axpy(&Cmnef, &C2, alpha_tot[j][i]);
        global_dpd_->buf4_close(&Cmnef);
      }
      global_dpd_->buf4_close(&C2);
    }

    sprintf(lbl, "%s %d", "SIA", L+i);
    global_dpd_->file2_init(&C1, PSIF_EOM_SIA, C_irr, A_OCC, A_VIR, lbl);
    global_dpd_->file2_scm(&C1, 0.0);
		C10 = 0.0;
    for (j=0;j<L;++j) {
      sprintf(lbl, "%s %d", "SIA", j);
      global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, A_OCC, A_VIR, lbl);
      global_dpd_->file2_axpy(&SIA, &C1, alpha_tot[j][i], 0);
      global_dpd_->file2_close(&SIA);
      if (params.full_matrix) {
		    sprintf(lbl, "%s %d", "S0", j);
			  psio_read_entry(PSIF_EOM_SIA, lbl, (char *) &S0, sizeof(double));
				C10 += alpha_tot[j][i] * S0;
		  }
    }
    if (params.full_matrix) {
		  sprintf(lbl, "%s %d", "S0", L+i);
		  psio_write_entry(PSIF_EOM_SIA, lbl, (char *) &C10, sizeof(double));
		}
    global_dpd_->file2_close(&C1);

    sprintf(lbl, "%s %d", "SIjAb", L+i);
    global_dpd_->buf4_init(&C2, PSIF_EOM_SIjAb, C_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, lbl);
    global_dpd_->buf4_scm(&C2, 0.0);
    for (j=0;j<L;++j) {
      sprintf(lbl, "%s %d", "SIjAb", j);
      global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, lbl);
      global_dpd_->buf4_axpy(&SIjAb, &C2, alpha_tot[j][i]);
      global_dpd_->buf4_close(&SIjAb);
    }
    global_dpd_->buf4_close(&C2);

    if (params.eom_ref > 0) {
      sprintf(lbl, "%s %d", "Sia", L+i);
      global_dpd_->file2_init(&C1, PSIF_EOM_Sia, C_irr, B_OCC, B_VIR, lbl);
      global_dpd_->file2_scm(&C1, 0.0);
      for (j=0;j<L;++j) {
        sprintf(lbl, "%s %d", "Sia", j);
        global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, B_OCC, B_VIR, lbl);
        global_dpd_->file2_axpy(&Sia, &C1, alpha_tot[j][i], 0);
        global_dpd_->file2_close(&Sia);
      }
      global_dpd_->file2_close(&C1);

      sprintf(lbl, "%s %d", "SIJAB", L+i);
      global_dpd_->buf4_init(&C2, PSIF_EOM_SIJAB, C_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, lbl);
      global_dpd_->buf4_scm(&C2, 0.0);
      for (j=0;j<L;++j) {
        sprintf(lbl, "%s %d", "SIJAB", j);
        global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, lbl);
        global_dpd_->buf4_axpy(&SIJAB, &C2, alpha_tot[j][i]);
        global_dpd_->buf4_close(&SIJAB);
      }
      global_dpd_->buf4_close(&C2);

      sprintf(lbl, "%s %d", "Sijab", L+i);
      global_dpd_->buf4_init(&C2, PSIF_EOM_Sijab, C_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, lbl);
      global_dpd_->buf4_scm(&C2, 0.0);
      for (j=0;j<L;++j) {
        sprintf(lbl, "%s %d", "Sijab", j);
        global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, lbl);
        global_dpd_->buf4_axpy(&Sijab, &C2, alpha_tot[j][i]);
        global_dpd_->buf4_close(&Sijab);
      }
      global_dpd_->buf4_close(&C2);
    }

  }

  /* Copy restart vectors to beginning of file */
  for (i=0; i<num; ++i) {
    sprintf(lbl, "%s %d", "CME", L+i);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, A_OCC, A_VIR, lbl);
    sprintf(lbl, "%s %d", "CME", i);
    global_dpd_->file2_copy(&CME, PSIF_EOM_CME, lbl);
    global_dpd_->file2_close(&CME);
    sprintf(lbl, "%s %d", "CMnEf", L+i);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    global_dpd_->buf4_copy(&CMnEf, PSIF_EOM_CMnEf, lbl);
    global_dpd_->buf4_close(&CMnEf);
		if (params.full_matrix) {
      sprintf(lbl, "%s %d", "C0", L+i);
			psio_read_entry(PSIF_EOM_CME, lbl, (char *) &CME0, sizeof(double));
      sprintf(lbl, "%s %d", "C0", i);
			psio_write_entry(PSIF_EOM_CME, lbl, (char *) &CME0, sizeof(double));
		}

    if (params.eom_ref > 0) {
      sprintf(lbl, "%s %d", "Cme", L+i);
      global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, B_OCC, B_VIR, lbl);
      sprintf(lbl, "%s %d", "Cme", i);
      global_dpd_->file2_copy(&Cme, PSIF_EOM_Cme, lbl);
      global_dpd_->file2_close(&Cme);
      sprintf(lbl, "%s %d", "CMNEF", L+i);
      global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, lbl);
      sprintf(lbl, "%s %d", "CMNEF", i);
      global_dpd_->buf4_copy(&CMNEF, PSIF_EOM_CMNEF, lbl);
      global_dpd_->buf4_close(&CMNEF);
      sprintf(lbl, "%s %d", "Cmnef", L+i);
      global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, lbl);
      sprintf(lbl, "%s %d", "Cmnef", i);
      global_dpd_->buf4_copy(&Cmnef, PSIF_EOM_Cmnef, lbl);
      global_dpd_->buf4_close(&Cmnef);
    }

    sprintf(lbl, "%s %d", "SIA", L+i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, A_OCC, A_VIR, lbl);
    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_copy(&SIA, PSIF_EOM_SIA, lbl);
    global_dpd_->file2_close(&SIA);
    sprintf(lbl, "%s %d", "SIjAb", L+i);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, AB_OCC, AB_VIR, AB_OCC, AB_VIR, 0, lbl);
    sprintf(lbl, "%s %d", "SIjAb", i);
    global_dpd_->buf4_copy(&SIjAb, PSIF_EOM_SIjAb, lbl);
    global_dpd_->buf4_close(&SIjAb);
		if (params.full_matrix) {
      sprintf(lbl, "%s %d", "S0", L+i);
			psio_read_entry(PSIF_EOM_SIA, lbl, (char *) &S0, sizeof(double));
      sprintf(lbl, "%s %d", "S0", i);
			psio_write_entry(PSIF_EOM_SIA, lbl, (char *) &S0, sizeof(double));
		}

    if (params.eom_ref > 0) {
      sprintf(lbl, "%s %d", "Sia", L+i);
      global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, B_OCC, B_VIR, lbl);
      sprintf(lbl, "%s %d", "Sia", i);
      global_dpd_->file2_copy(&Sia, PSIF_EOM_Sia, lbl);
      global_dpd_->file2_close(&Sia);
      sprintf(lbl, "%s %d", "SIJAB", L+i);
      global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, AA_OCC, AA_VIR, AA_OCC, AA_VIR, 0, lbl);
      sprintf(lbl, "%s %d", "SIJAB", i);
      global_dpd_->buf4_copy(&SIJAB, PSIF_EOM_SIJAB, lbl);
      global_dpd_->buf4_close(&SIJAB);
      sprintf(lbl, "%s %d", "Sijab", L+i);
      global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, BB_OCC, BB_VIR, BB_OCC, BB_VIR, 0, lbl);
      sprintf(lbl, "%s %d", "Sijab", i);
      global_dpd_->buf4_copy(&Sijab, PSIF_EOM_Sijab, lbl);
      global_dpd_->buf4_close(&Sijab);
    }
  }

  free_block(alpha_tot);

#ifdef TIME_CCEOM
timer_off("RESTART");
#endif

  return;
}


}} // namespace psi::cceom
