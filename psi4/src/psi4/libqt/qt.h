/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

/*!
** \file
** \brief Header file for the Quantum Trio Library
** \ingroup QT
**
** David Sherrill 1994
**
** Modifications by Daniel Crawford 1996, 1997
*/

#pragma once

#include <cstddef>
#include <string>

#include "psi4/pragma.h"
#include "psi4/psi4-dec.h"

#include "blas_level1.h"
#include "blas_level2.h"
#include "blas_level3.h"
#include "lapack.h"

namespace psi {
// I think this is forward-declaring class Options -CDS
class Options;
class Wavefunction;

void dx_write(std::shared_ptr<Wavefunction> wfn, Options& options, double** D);
void dx_read(double** V_eff, double* phi_ao, double* phi_so, int nao, int nso, double** u);
void fill_sym_matrix(double** A, int size);
double combinations(int n, int k);
double factorial(int n);
void schmidt(double** A, int rows, int cols, std::string out_fname);
PSI_API int schmidt_add(double** A, int rows, int cols, double* v);
void normalize(double** A, int rows, int cols);
double invert_matrix(double** a, double** y, int N, std::string out_fname);
void solve_2x2_pep(double** H, double S, double* evals, double** evecs);
PSI_API void reorder_qt(int* docc_in, int* socc_in, int* frozen_docc_in, int* frozen_uocc_in, int* order,
                        int* orbs_per_irrep, int nirreps);
PSI_API void reorder_qt_uhf(int* docc, int* socc, int* frozen_docc, int* frozen_uocc, int* order_alpha, int* order_beta,
                            int* orbspi, int nirreps);
// int ras_set(int nirreps, int nbfso, int freeze_core, int *orbspi,
//      int *docc, int *socc, int *frdocc, int *fruocc,
//      int **ras_opi, int *order, int ras_type);
// int ras_set2(int nirreps, int nbfso, int delete_fzdocc,
//      int delete_restrdocc, int *orbspi,
//      int *docc, int *socc, int *frdocc, int *fruocc,
//      int *restrdocc, int *restruocc, int **ras_opi, int *order,
//      int ras_type, int hoffmann, Options& options);
int ras_set3(int nirreps, int nmo, int* orbspi, int* docc, int* socc, int* frdocc, int* fruocc, int* restrdocc,
             int* restruocc, int** ras_opi, int* core_guess, int* order, int ras_type, bool is_mcscf, Options& options);
void newmm_rking(double** A, int transa, double** B, int transb, double** C, int num_rows, int num_links, int num_cols,
                 double alpha, double beta);
double dot_block(double** A, double** B, int rows, int cols, double alpha);
void dirprd_block(double** A, double** B, int rows, int cols);
int pople(double** A, double* x, int dimen, int num_vecs, double tolerance, std::string out_fname, int print_lvl);
void mat_print(double** A, int rows, int cols, std::string out_fname);

void timer_init();
void timer_done();
void timer_on(const std::string& key);
void timer_off(const std::string& key);
void parallel_timer_on(const std::string& key, int thread_rank);
void parallel_timer_off(const std::string& key, int thread_rank);
void start_skip_timers();
void stop_skip_timers();

void print_block(double*, int, int, FILE*);

int david(double** A, int N, int M, double* eps, double** v, double cutoff, int print);

int* get_frzcpi();
int* get_frzvpi();
int cc_excited(const char* wfn);
int cc_excited(std::string wfn);
void free_3d_array(double*** A, int p, int q);
double*** init_3d_array(int p, int q, int r);

#define MAX_RAS_SPACES 4
}  // namespace psi
