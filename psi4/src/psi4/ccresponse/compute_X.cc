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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void init_X(const char *pert, int irrep, double omega);
void sort_X(const char *pert, int irrep, double omega);
void cc2_sort_X(const char *pert, int irrep, double omega);
void X1_build(const char *pert, int irrep, double omega);
void X2_build(const char *pert, int irrep, double omega);
void cc2_X1_build(const char *pert, int irrep, double omega);
void cc2_X2_build(const char *pert, int irrep, double omega);
double converged(const char *pert, int irrep, double omega);
void save_X(const char *pert, int irrep, double omega);
void print_X(const char *pert, int irrep, double omega);
void update_X(const char *pert, int irrep, double omega);
void diis(int iter, const char *pert, int irrep, double omega);
double pseudopolar(const char *pert, int irrep, double omega);
void cleanup(void);
void exit_io(void);
void amp_write(const char *pert, int irrep, double omega);

void analyze(const char *pert, int irrep, double omega);

void compute_X(const char *pert, int irrep, double omega)
{
  int i, iter=0, done=0;
  double rms, polar, X2_norm;
  char lbl[32];
  dpdbuf4 X2;

  timer_on("compute_X");

  outfile->Printf( "\n\tComputing %s-Perturbed Wave Function (%5.3f E_h).\n", pert, omega);
  init_X(pert, irrep, omega);
  outfile->Printf( "\tIter   Pseudopolarizability       RMS \n");
  outfile->Printf( "\t----   --------------------   -----------\n");


  if (params.wfn == "CC2")
    cc2_sort_X(pert, irrep, omega);
  else
    sort_X(pert, irrep, omega);
  polar = -2.0*pseudopolar(pert, irrep, omega);
  outfile->Printf( "\t%4d   %20.12f\n", iter, polar);


  for(iter=1; iter <= params.maxiter; iter++) {

    if (params.wfn == "CC2") {
      cc2_sort_X(pert, irrep, omega);
      cc2_X1_build(pert, irrep, omega);
      cc2_X2_build(pert, irrep, omega);
    }
    else {
      sort_X(pert, irrep, omega);
      X1_build(pert, irrep, omega);
      X2_build(pert, irrep, omega);
    }
    update_X(pert, irrep, omega);
    rms = converged(pert, irrep, omega);
    if(rms <= params.convergence) {
      done = 1;
      save_X(pert, irrep, omega);
      if (params.wfn == "CC2")
        cc2_sort_X(pert, irrep, omega);
      else
        sort_X(pert, irrep, omega);
      outfile->Printf( "\t-----------------------------------------\n");
      outfile->Printf( "\tConverged %s-Perturbed Wfn to %4.3e\n", pert, rms);
      if(params.print & 2) {
        sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
        global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
        X2_norm = global_dpd_->buf4_dot_self(&X2);
        global_dpd_->buf4_close(&X2);
        X2_norm = sqrt(X2_norm);
        outfile->Printf( "\tNorm of the converged X2 amplitudes %20.15f\n", X2_norm);
        amp_write(pert, irrep, omega);
      }

      break;
    }
    if(params.diis) diis(iter, pert, irrep, omega);
    save_X(pert, irrep, omega);
    if (params.wfn == "CC2")
      cc2_sort_X(pert, irrep, omega);
    else
      sort_X(pert, irrep, omega);

    polar = -2.0*pseudopolar(pert, irrep, omega);
    outfile->Printf( "\t%4d   %20.12f    %4.3e\n", iter, polar, rms);


  }
  if(!done) {

    dpd_close(0);
    cleanup();
    exit_io();
    throw PsiException("Failed to converge perturbed wavefunction",__FILE__,__LINE__);
  }

  /* Clean up disk space */
  psio_close(PSIF_CC_DIIS_AMP, 0);
  psio_close(PSIF_CC_DIIS_ERR, 0);

  psio_open(PSIF_CC_DIIS_AMP, 0);
  psio_open(PSIF_CC_DIIS_ERR, 0);

  for(i=PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) {
    psio_close(i,0);
    psio_open(i,0);
  }

  if(params.analyze) analyze(pert, irrep, omega);

  /*  print_X(pert, irrep, omega); */

  timer_off("compute_X");
}

}} // namespace psi::ccresponse
