/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file    print.cc
    \ingroup optking
    \brief memory allocation
*/

#include "print.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include <cstdarg>

namespace opt {

// Argument 1 determines file for printing in Psi4.  2nd argument ignored.
// Argument 2 determines file for printing in QChem. 1st argument ignored.
void oprintf(const std::string psi_fp, const FILE *qc_fp, const char* format,...) {

  char line[256];
  va_list args;
  va_start(args, format);
  vsprintf(line, format, args);
  va_end(args);

#if defined(OPTKING_PACKAGE_PSI)
  auto mode = std::ostream::app;
  auto printer = psi_fp=="outfile" ? psi::outfile
                                   : std::make_shared<psi::PsiOutStream>(psi_fp, mode);

  printer->Printf("%s", line);
#elif defined(OPTKING_PACKAGE_QCHEM)
  fprintf(qc_fp, "%s", line);
#endif
}

void offlush_out() {
#if defined(OPTKING_PACKAGE_PSI)
  std::shared_ptr<psi::PsiOutStream> printer(psi::outfile);
#elif defined(OPTKING_PACKAGE_QCHEM)
  fflush(qc_outfile);
#endif
}

// oprintf_out is always to primary output file
void oprintf_out(const char* format,...) {
  char line[512];
  va_list args;
  va_start(args, format);
  vsprintf(line, format, args);
  va_end(args);

#if defined(OPTKING_PACKAGE_PSI)
  *(psi::outfile->stream()) << line;
#elif defined(OPTKING_PACKAGE_QCHEM)
  fprintf(qc_outfile, "%s", line);
#endif
}

void oprint_matrix(const std::string psi_fp, const FILE *qc_fp, double **A, const int nrow, const int ncol) {
  int col=0;
  const int max_col = 8;

  for (int i=0; i<nrow; ++i) {
    for (int j=0; j<ncol; ++j) {
      oprintf(psi_fp, qc_fp, "%10.6f", A[i][j]);
      ++col;
      if ((col == max_col) && (j != ncol-1)) {
        oprintf(psi_fp, qc_fp, "\n");
        col = 0;
      }
    }
    oprintf(psi_fp, qc_fp, "\n");
    col = 0;
  }
  return;
}

void oprint_matrix_out(double **A, const int nrow, const int ncol) {
  int col=0;
  const int max_col = 8;

  for (int i=0; i<nrow; ++i) {
    for (int j=0; j<ncol; ++j) {
      oprintf_out("%10.6f", A[i][j]);
      ++col;
      if ((col == max_col) && (j != ncol-1)) {
        oprintf_out("\n");
        col = 0;
      }
    }
    oprintf_out("\n");
    col = 0;
  }
  return;
}

void oprint_matrix_out_precise(double **A, const int nrow, const int ncol) {
  int col=0;
  const int max_col = 4;

  for (int i=0; i<nrow; ++i) {
    for (int j=0; j<ncol; ++j) {
      oprintf_out("%20.15f", A[i][j]);
      ++col;
      if ((col == max_col) && (j != ncol-1)) {
        oprintf_out("\n");
        col = 0;
      }
    }
    oprintf_out("\n");
    col = 0;
  }
  return;
}

void oprint_array(const std::string psi_fp, const FILE *qc_fp, double *A, const int ncol) {
  int col=0;
  const int max_col = 8;

  for (int j=0; j<ncol; ++j) {
    oprintf(psi_fp, qc_fp, "%10.6f", A[j]);
    ++col;
    if ((col == max_col) && (j != ncol-1)) {
      oprintf(psi_fp, qc_fp, "\n");
      col = 0;
    }
  }
  oprintf(psi_fp, qc_fp, "\n");
  return;
}

void oprint_array_out(double *A, const int ncol) {
  int col=0;
  const int max_col = 8;

  for (int j=0; j<ncol; ++j) {
    oprintf_out("%10.6f", A[j]);
    ++col;
    if ((col == max_col) && (j != ncol-1)) {
      oprintf_out("\n");
      col = 0;
    }
  }
  oprintf_out("\n");
  return;
}

void oprint_array_out_precise(double *A, const int ncol) {
  int col=0;
  const int max_col = 4;

  for (int j=0; j<ncol; ++j) {
    oprintf_out("%20.15f", A[j]);
    ++col;
    if ((col == max_col) && (j != ncol-1)) {
      oprintf_out("\n");
      col = 0;
    }
  }
  oprintf_out("\n");
  return;
}

}
