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

#include "dcft.h"
#include "defines.h"
#include <vector>
#include <cmath>
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libiwl/iwl.hpp"



namespace psi{ namespace dcft{

DCFTSolver::DCFTSolver(SharedWavefunction ref_wfn, Options &options):
        Wavefunction(options)
{
    reference_wavefunction_ = ref_wfn;
    shallow_copy(ref_wfn);
    Ca_ = ref_wfn->Ca()->clone();
    Cb_ = ref_wfn->Cb()->clone();
    Da_ = ref_wfn->Da()->clone();
    Db_ = ref_wfn->Db()->clone();
    Fa_ = ref_wfn->Fa()->clone();
    Fb_ = ref_wfn->Fb()->clone();

    maxiter_            = options.get_int("MAXITER");
    print_              = options.get_int("PRINT");
    maxdiis_            = options.get_int("DIIS_MAX_VECS");
    mindiisvecs_        = options.get_int("DIIS_MIN_VECS");
    regularizer_        = options.get_double("TIKHONOW_OMEGA");
    orbital_level_shift_= options.get_double("ORBITAL_LEVEL_SHIFT");
    diis_start_thresh_  = options.get_double("DIIS_START_CONVERGENCE");
    orbitals_threshold_ = options.get_double("R_CONVERGENCE");
    cumulant_threshold_ = options.get_double("R_CONVERGENCE");
    int_tolerance_      = options.get_double("INTS_TOLERANCE");
    energy_level_shift_ = options.get_double("ENERGY_LEVEL_SHIFT");

    if (!options_["E_CONVERGENCE"].has_changed())
        energy_threshold_ = options.get_double("R_CONVERGENCE");
    else
        energy_threshold_ = options.get_double("E_CONVERGENCE");

    psio_->open(PSIF_DCFT_DPD, PSIO_OPEN_OLD);

    exact_tau_ = false;
    if (options.get_str("DCFT_FUNCTIONAL") == "DC-12"
            || options.get_str("DCFT_FUNCTIONAL") == "ODC-12"
            || options.get_str("DCFT_FUNCTIONAL") == "ODC-13") exact_tau_ = true;

    orbital_optimized_ = false;
    if (options.get_str("DCFT_FUNCTIONAL") == "ODC-06"
            || options.get_str("DCFT_FUNCTIONAL") == "ODC-12"
            || options.get_str("DCFT_FUNCTIONAL") == "ODC-13") orbital_optimized_ = true;

    // Sets up the memory, and orbital info
    init();
}

/**
 * Computes A = A + alpha * B, writing the result back to A
 */
void DCFTSolver::dpd_buf4_add(dpdbuf4 *A, dpdbuf4 *B, double alpha)
{
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(A, h);
        global_dpd_->buf4_mat_irrep_init(B, h);
        global_dpd_->buf4_mat_irrep_rd(A, h);
        global_dpd_->buf4_mat_irrep_rd(B, h);

        #pragma omp parallel for
        for(int row = 0; row < A->params->rowtot[h]; ++row){
            for(int col = 0; col < A->params->coltot[h]; ++col){
                A->matrix[h][row][col] += alpha * B->matrix[h][row][col];
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(A, h);
        global_dpd_->buf4_mat_irrep_close(A, h);
        global_dpd_->buf4_mat_irrep_close(B, h);
    }
}

DCFTSolver::~DCFTSolver()
{
}

}} // Namespaces
