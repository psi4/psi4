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

#include "psi4/psi4-dec.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libtrans/integraltransform.h"

#include "psi4/libtrans/mospace.h"
#include "adc.h"

namespace psi{ namespace adc {

double
ADCWfn::rhf_init_tensors()
{
    bool do_pr;
    double ePR2, sq_norm, energy;
    dpdbuf4 K, V;
    dpdfile2 A;
    // Setting up and initialize the integraltransform object
    std::vector<std::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::occ);
    spaces.push_back(MOSpace::vir);//printf("madeok\n");
    _ints = new IntegralTransform(reference_wavefunction_, spaces, IntegralTransform::Restricted);//printf("madeok\n");
    _ints->set_keep_iwl_so_ints(true);
    _ints->set_keep_dpd_so_ints(true);
    dpd_set_default(_ints->get_dpd_id());
    // Make (OV|OV) integrals
    outfile->Printf( "\n\t==> Transforming (OV|OV) Integrals <==\n");
    _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    // Make (OO|VV) integrals
    outfile->Printf( "\n\t==> Transforming (OO|VV) Integrals <==\n");
    _ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir);
    // Make (OO|OV) integrals
    outfile->Printf( "\n\t==> Transforming (OV|OO) Integrals <==\n");
    _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::occ);
    // Make (OV|VV) integrals
    outfile->Printf( "\n\t==> Transforming (OV|VV) Integrals <==\n");
    _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::vir, MOSpace::vir);


    // Preparing MP1 amplitudes then calculating MP2 energy
    // and use of LMO is not considered in this code.
    // In ADC(2) calculation, 2 <ij|ab> - <ij|ba> typed
    // integral list is needed. So making this too.
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    psio_->open(PSIF_ADC, PSIO_OPEN_NEW);

    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&V, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    global_dpd_->buf4_sort(&V, PSIF_LIBTRANS_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "MO Ints V1243");
    global_dpd_->buf4_scmcopy(&V, PSIF_LIBTRANS_DPD, "MO Ints 2 V1234 - V1243", 2.0);
    global_dpd_->buf4_copy(&V, PSIF_ADC, "K1234");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints 2 V1234 - V1243");
    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints V1243)");
    global_dpd_->buf4_axpy(&V, &K, -1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_close(&K);

    global_dpd_->buf4_init(&K, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "K1234");

    for(int h = 0;h < nirrep_;h++){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int ij = 0;ij < K.params->rowtot[h];ij++){
            int i = K.params->roworb[h][ij][0];
            int j = K.params->roworb[h][ij][1];
            for(int ab = 0;ab < K.params->coltot[h];ab++){
                int a = K.params->colorb[h][ab][0];
                int b = K.params->colorb[h][ab][1];
                K.matrix[h][ij][ab] /= aocce_[i] + aocce_[j] - avire_[a] - avire_[b];
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&K, h);
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_sort(&K, PSIF_ADC, pqsr, ID("[O,O]"), ID("[V,V]"), "K1243");
    global_dpd_->buf4_scmcopy(&K, PSIF_ADC, "2 K1234 - K1243", 2.0);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "2 K1234 - K1243");
    global_dpd_->buf4_init(&V, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "K1243");
    global_dpd_->buf4_axpy(&V, &K, -1);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    energy = global_dpd_->buf4_dot(&V, &K);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "K1234");
    sq_norm = 1 + global_dpd_->buf4_dot(&V, &K);

    do_pr = options_.get_bool("PR");

    outfile->Printf( "\n\t==> Ground State <==\n");
    if(!do_pr) outfile->Printf( "->");
    outfile->Printf( "\tMP2 energy    = %20.14f\n", energy);
    outfile->Printf( "\t[Squared-norm of MP1 wavefunction    = %10.7f]\n", sq_norm);
    // Partially renormalized MP2 energy and the MP1 wavefunction
    // Reference: IJQC 78 (2000) 226, CPL 443 (2007) 389.
    global_dpd_->file2_init(&A, PSIF_ADC, 0, ID('O'), ID('O'), "RHO_OO");
    global_dpd_->contract442(&V, &K, &A, 0, 0, 1, 0);
    global_dpd_->buf4_close(&K);


    global_dpd_->file2_mat_init(&A);
    global_dpd_->file2_mat_rd(&A);
    global_dpd_->buf4_init(&K, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "tilde K1234");

    for(int h = 0;h < nirrep_;h++){
        global_dpd_->buf4_mat_irrep_init(&V, h);
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&V, h);
        for(int ij = 0;ij < V.params->rowtot[h];ij++){
            int i = V.params->roworb[h][ij][0];
            int j = V.params->roworb[h][ij][1];
            int isym = V.params->psym[i];
            int jsym = V.params->qsym[j];
            int I = i - V.params->poff[isym];
            int J = j - V.params->qoff[jsym];
            double Nij = 1 + (A.matrix[isym][I][I] + A.matrix[jsym][J][J]) / 2;
            for(int ab = 0;ab < V.params->coltot[h];ab++)
                K.matrix[h][ij][ab] = V.matrix[h][ij][ab] / Nij;
        }
        global_dpd_->buf4_mat_irrep_wrt(&K, h);
        global_dpd_->buf4_mat_irrep_close(&K, h);
        global_dpd_->buf4_mat_irrep_close(&V, h);
    }
    global_dpd_->buf4_close(&V);
    global_dpd_->file2_mat_close(&A);
    global_dpd_->file2_close(&A);

    global_dpd_->buf4_sort(&K, PSIF_ADC, pqsr, ID("[O,O]"), ID("[V,V]"), "tilde K1243");
    global_dpd_->buf4_scmcopy(&K, PSIF_ADC, "tilde 2 K1234 - K1243", 2.0);
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "tilde 2 K1234 - K1243");
    global_dpd_->buf4_init(&V, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "tilde K1243");
    global_dpd_->buf4_axpy(&V, &K, -1.0);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    ePR2 = global_dpd_->buf4_dot(&V, &K);
    global_dpd_->buf4_close(&V);
    global_dpd_->buf4_init(&V, PSIF_ADC, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, "tilde K1234");
    sq_norm = 1 + global_dpd_->buf4_dot(&V, &K);

#if DEBUG_
    outfile->Printf( ">> In init_tensor <<\n");
    global_dpd_->buf4_print(&K, outfile, 1);
#endif

    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&V);

    if(do_pr) outfile->Printf( "->");
    outfile->Printf( "\tPR-MP2 energy = %20.14f\n", ePR2);
    outfile->Printf( "\t[Squared-norm of PR-MP1 wavefunction = %10.7f]\n\n", sq_norm);


    if(do_pr) energy = ePR2;

    // Reordering each ERIs other than (OO|VV) type from Mulliken to Dirac notation
    // for convenience of the evaluation of the sigma tensor

    // Sort(prqs): <OV|OV> <-- (OO|VV)
    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
    global_dpd_->buf4_sort(&V, PSIF_LIBTRANS_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "MO Ints <OV|OV>");
    global_dpd_->buf4_close(&V);

    // Sort(pqrs): <OO|VO> <-- (OV|OO)
    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"), ID("[O,V]"), ID("[O>=O]+"), 0, "MO Ints (OV|OO)");
    global_dpd_->buf4_sort(&V, PSIF_LIBTRANS_DPD, prqs, ID("[O,O]"), ID("[V,O]"), "MO Ints <OO|VO>");
    global_dpd_->buf4_close(&V);

    // Sort(prqs): <OV|VV> <-- (OV|VV)
    global_dpd_->buf4_init(&V, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"), ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
    global_dpd_->buf4_sort(&V, PSIF_LIBTRANS_DPD, prqs, ID("[O,V]"), ID("[V,V]"), "MO Ints <OV|VV>");
    global_dpd_->buf4_close(&V);

    return energy;
}

}} // End Namespaces
