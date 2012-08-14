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
DCFTSolver::compute_lambda_residual()
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
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_copy(&G, PSIF_DCFT_DPD, "R <OO|VV>");
    dpd_buf4_close(&G);
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");

    // R_IJAB += F_IJAB
    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "F <OO|VV>");
    dpd_buf4_add(&R, &F, 1.0);
    dpd_buf4_close(&F);

    for(int h = 0; h < nirrep_; ++h)
        nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += dpd_buf4_dot_self(&R);
    dpd_buf4_close(&R);

    // R_IjAb = G_IjAb
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    dpd_buf4_copy(&G, PSIF_DCFT_DPD, "R <Oo|Vv>");
    dpd_buf4_close(&G);
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");

    // R_IjAb += F_IjAb
    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "F <Oo|Vv>");
    dpd_buf4_add(&R, &F, 1.0);
    dpd_buf4_close(&F);

    for(int h = 0; h < nirrep_; ++h)
        nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += dpd_buf4_dot_self(&R);
    dpd_buf4_close(&R);

    // R_ijab = G_ijab
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_copy(&G, PSIF_DCFT_DPD, "R <oo|vv>");
    dpd_buf4_close(&G);
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");

    // R_ijab += F_ijab
    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "F <oo|vv>");
    dpd_buf4_add(&R, &F, 1.0);
    dpd_buf4_close(&F);

    for(int h = 0; h < nirrep_; ++h)
        nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += dpd_buf4_dot_self(&R);
    dpd_buf4_close(&R);

    dcft_timer_off("DCFTSolver::compute_lambda_residual()");


    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    sum_G += dpd_buf4_dot_self(&G);
//    dpd_buf4_print(&G, outfile, 1);
    dpd_buf4_close(&G);

    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "F <OO|VV>");
    sum_F += dpd_buf4_dot_self(&F);
//    dpd_buf4_print(&F, outfile, 1);
    dpd_buf4_close(&F);


    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    sum_G += dpd_buf4_dot_self(&G);
//    dpd_buf4_print(&G, outfile, 1);
    dpd_buf4_close(&G);

    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "F <Oo|Vv>");
    sum_F += dpd_buf4_dot_self(&F);
//    dpd_buf4_print(&F, outfile, 1);
    dpd_buf4_close(&F);

    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    sum_G += dpd_buf4_dot_self(&G);
//    dpd_buf4_print(&G, outfile, 1);
    dpd_buf4_close(&G);

    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "F <oo|vv>");
    sum_F += dpd_buf4_dot_self(&F);
//    dpd_buf4_print(&F, outfile, 1);
    dpd_buf4_close(&F);

//    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
//                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
//    dpd_buf4_print(&G, outfile, 1);
//    dpd_buf4_close(&G);

//    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
//                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
//    dpd_buf4_print(&G, outfile, 1);
//    dpd_buf4_close(&G);

//    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
//                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
//    dpd_buf4_print(&G, outfile, 1);
//    dpd_buf4_close(&G);

//    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
//    dpdfile2 FF;
//    dpd_file2_init(&FF, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
//    dpd_file2_print(&FF, outfile);
//    dpd_file2_close(&FF);
//    dpd_file2_init(&FF, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
//    dpd_file2_print(&FF, outfile);
//    dpd_file2_close(&FF);
//    psio_->close(PSIF_LIBTRANS_DPD, 1);

    aocc_tau_ = SharedMatrix(new Matrix("MO basis Tau (Alpha Occupied)", nirrep_, naoccpi_, naoccpi_));
    avir_tau_ = SharedMatrix(new Matrix("MO basis Tau (Alpha Virtual)", nirrep_, navirpi_, navirpi_));

    dpdfile2 TT_o, TT_v;
    dpd_file2_init(&TT_o, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    dpd_file2_init(&TT_v, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    dpd_file2_mat_init(&TT_o);
    dpd_file2_mat_init(&TT_v);
    dpd_file2_mat_rd(&TT_o);
    dpd_file2_mat_rd(&TT_v);

    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int j = 0; j < naoccpi_[h]; ++j){
                aocc_tau_->set(h, i, j, TT_o.matrix[h][i][j]);
            }
        }
        for(int a = 0; a < navirpi_[h]; ++a){
            for(int b = 0; b < navirpi_[h]; ++b){
                avir_tau_->set(h, a, b, TT_v.matrix[h][a][b]);
            }
        }
    }

    dpd_file2_close(&TT_o);
    dpd_file2_close(&TT_v);

//    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
//    psio_->close(PSIF_LIBTRANS_DPD, 1);

    fprintf(outfile, "Tr(Tau<O|O>) = %5.3f Tr(Tau<V|V>) = %5.3f Tr(Tau) = %5.3e \n", aocc_tau_->trace(), avir_tau_->trace(), aocc_tau_->trace() + avir_tau_->trace());
    fprintf(outfile, "Residual: F: %8.5e G: %8.5e\n", sqrt(sum_F / nElements), sqrt(sum_G / nElements));

    return sqrt(sumSQ / nElements);
}


/**
 * Builds the new lambda tensor from the intermediates
 */
void
DCFTSolver::update_lambda_from_residual()
{
    dcft_timer_on("DCFTSolver::update_lambda_from_residual()");

    dpdbuf4 L, D, R;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    /*
     * Lambda_ijab += R_ijab / D_ijab
     */
    // L_IJAB += R_IJAB / D_IJAB
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "D <OO|VV>");
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");
    dpd_buf4_dirprd(&D, &R);
    dpd_buf4_close(&D);
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_add(&L, &R, 1.0);
    dpd_buf4_close(&R);
    dpd_buf4_close(&L);

    // L_IjAb += R_IjAb / D_IjAb
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
    dpd_buf4_dirprd(&D, &R);
    dpd_buf4_close(&D);
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_add(&L, &R, 1.0);
    dpd_buf4_close(&R);
    dpd_buf4_close(&L);

    // L_IJAB += R_ijab / D_ijab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>=o]+"), ID("[v>=v]+"), 0, "D <oo|vv>");
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");
    dpd_buf4_dirprd(&D, &R);
    dpd_buf4_close(&D);
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_add(&L, &R, 1.0);
    dpd_buf4_close(&R);
    dpd_buf4_close(&L);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    dcft_timer_off("DCFTSolver::update_lambda_from_residual()");
}

}} // Namespaces


