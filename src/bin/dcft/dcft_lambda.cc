/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include "dcft.h"
#include <libdpd/dpd.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <psifiles.h>
#include "defines.h"

#include <cmath>

namespace psi{ namespace dcft{

/**
 * Computes the residual for the lambda equations
 * R = G + F
 * @return RMS residual
 */
double
DCFTSolver::compute_cumulant_residual()
{
    dcft_timer_on("DCFTSolver::compute_lambda_residual()");

    dpdbuf4 R, G, F;
    double sumSQ = 0.0;
    double sum_F = 0.0;
    double sum_G = 0.0;
    size_t nElements = 0;

    /*
     * R_ijab = G_ijab + F_ijab
     */

    // R_IJAB = G_IJAB
    dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_->buf4_copy(&G, PSIF_DCFT_DPD, "R <OO|VV>");
    dpd_->buf4_close(&G);
    dpd_->buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");

    // R_IJAB += F_IJAB
    dpd_->buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "F <OO|VV>");
    dpd_buf4_add(&R, &F, 1.0);
    dpd_->buf4_close(&F);

    for(int h = 0; h < nirrep_; ++h)
        nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += dpd_->buf4_dot_self(&R);
    dpd_->buf4_close(&R);

    // R_IjAb = G_IjAb
    dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    dpd_->buf4_copy(&G, PSIF_DCFT_DPD, "R <Oo|Vv>");
    dpd_->buf4_close(&G);
    dpd_->buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");

    // R_IjAb += F_IjAb
    dpd_->buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "F <Oo|Vv>");
    dpd_buf4_add(&R, &F, 1.0);
    dpd_->buf4_close(&F);

    for(int h = 0; h < nirrep_; ++h)
        nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += dpd_->buf4_dot_self(&R);
    dpd_->buf4_close(&R);

    // R_ijab = G_ijab
    dpd_->buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_->buf4_copy(&G, PSIF_DCFT_DPD, "R <oo|vv>");
    dpd_->buf4_close(&G);
    dpd_->buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");

    // R_ijab += F_ijab
    dpd_->buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "F <oo|vv>");
    dpd_buf4_add(&R, &F, 1.0);
    dpd_->buf4_close(&F);

    for(int h = 0; h < nirrep_; ++h)
        nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += dpd_->buf4_dot_self(&R);
    dpd_->buf4_close(&R);

    dcft_timer_off("DCFTSolver::compute_lambda_residual()");

    if (nElements > 0) return sqrt(sumSQ / nElements);
    else return 0.0;
}


/**
 * Builds the new lambda tensor from the intermediates
 */
void
DCFTSolver::update_cumulant_jacobi()
{
    dcft_timer_on("DCFTSolver::update_lambda_from_residual()");

    dpdbuf4 L, D, R;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    /*
     * Lambda_ijab += R_ijab / D_ijab
     */
    // L_IJAB += R_IJAB / D_IJAB
    dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "D <OO|VV>");
    dpd_->buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");
    dpd_->buf4_dirprd(&D, &R);
    dpd_->buf4_close(&D);
    dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_add(&L, &R, 1.0);
    dpd_->buf4_close(&R);
    dpd_->buf4_close(&L);

    // L_IjAb += R_IjAb / D_IjAb
    dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    dpd_->buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
    dpd_->buf4_dirprd(&D, &R);
    dpd_->buf4_close(&D);
    dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_add(&L, &R, 1.0);
    dpd_->buf4_close(&R);
    dpd_->buf4_close(&L);

    // L_IJAB += R_ijab / D_ijab
    dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>=o]+"), ID("[v>=v]+"), 0, "D <oo|vv>");
    dpd_->buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");
    dpd_->buf4_dirprd(&D, &R);
    dpd_->buf4_close(&D);
    dpd_->buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_add(&L, &R, 1.0);
    dpd_->buf4_close(&R);
    dpd_->buf4_close(&L);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    dcft_timer_off("DCFTSolver::update_lambda_from_residual()");
}

}} // Namespaces


