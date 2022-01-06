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

/*! \file    print.h
    \ingroup optking
    \brief header for print functions
*/

#ifndef _opt_print_h_
#define _opt_print_h_

#include <cstdlib>
#include <cstdio>
#include <string>

#include "package.h"

namespace opt {

// Functions to all printing in either psi or qchem.
void oprintf(std::string psi_fp, const FILE *qc_fp, const char* format,...);

void oprintf_out(const char* format,...);

void oprint_matrix(const std::string psi_fp, const FILE *qc_fp, double **A, const int x, const int y);

void oprint_matrix_out(double **A, const int x, const int y);

void oprint_matrix_out_precise(double **A, const int x, const int y);

void oprint_array(const std::string psi_fp, const FILE *qc_fp, double *A, const int x);

void oprint_array_out(double *A, const int x);

void oprint_array_out_precise(double *A, const int x);

void offlush_out();

}

#endif
