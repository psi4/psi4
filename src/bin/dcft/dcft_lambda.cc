#include "dcft.h"
#include <libdpd/dpd.h>
#include <libtrans/integraltransform.h>
#include "psifiles.h"
#include "defines.h"

namespace psi{ namespace dcft{

/**
 * Computes the residual for the lambda equations
 * R = G + T - A
 * @return RMS residual
 */
double
DCFTSolver::compute_lambda_residual()
{
    dpdbuf4 R, G, T, A;
    double sumSQ = 0.0;
    size_t nElements = 0;

    /*
     * R_ijab = G_ijab + T_ijab - A_ijab
     */

    // R_IJAB = G_IJAB
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    dpd_buf4_copy(&G, PSIF_DCFT_DPD, "R <OO|VV>");
    dpd_buf4_close(&G);
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "R <OO|VV>");
    for(int h = 0; h < nirrep_; ++h)
        nElements += R.params->coltot[h] * R.params->rowtot[h];
    if(!options_.get_bool("IGNORE_TAU")){
        // R_IJAB += T_IJAB
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "T <OO|VV>");
        dpd_buf4_add(&R, &T, 1.0);
        dpd_buf4_close(&T);
    }
    // R_IJAB -= A_IJAB
    dpd_buf4_init(&A, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "A <OO|VV>");
    dpd_buf4_add(&R, &A, -1.0);
    dpd_buf4_close(&A);
    sumSQ += dpd_buf4_dot_self(&R);
    dpd_buf4_close(&R);

    // R_IjAb = G_IjAb
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    dpd_buf4_copy(&G, PSIF_DCFT_DPD, "R <Oo|Vv>");
    dpd_buf4_close(&G);
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
    for(int h = 0; h < nirrep_; ++h)
        nElements += R.params->coltot[h] * R.params->rowtot[h];
    if(!options_.get_bool("IGNORE_TAU")){
        // R_IjAb += T_IjAb
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "T <Oo|Vv>");
        dpd_buf4_add(&R, &T, 1.0);
        dpd_buf4_close(&T);
    }
    // R_IjAb -= A_IjAb
    dpd_buf4_init(&A, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "A <Oo|Vv>");
    dpd_buf4_add(&R, &A, -1.0);
    dpd_buf4_close(&A);
    sumSQ += dpd_buf4_dot_self(&R);
    dpd_buf4_close(&R);

    // R_ijab = G_ijab
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_buf4_copy(&G, PSIF_DCFT_DPD, "R <oo|vv>");
    dpd_buf4_close(&G);
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "R <oo|vv>");
    for(int h = 0; h < nirrep_; ++h)
        nElements += R.params->coltot[h] * R.params->rowtot[h];
    if(!options_.get_bool("IGNORE_TAU")){
        // R_ijab += T_ijab
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "T <oo|vv>");
        dpd_buf4_add(&R, &T, 1.0);
        dpd_buf4_close(&T);
    }
    // R_ijab -= A_ijab
    dpd_buf4_init(&A, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "A <oo|vv>");
    dpd_buf4_add(&R, &A, -1.0);
    dpd_buf4_close(&A);
    sumSQ += dpd_buf4_dot_self(&R);
    dpd_buf4_close(&R);

    return sqrt(sumSQ / nElements);
}


/**
 * Builds the new lambda tensor from the intermediates
 */
void
DCFTSolver::update_lambda_from_residual()
{
    dpdbuf4 L, D, R;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    /*
     * Lambda_ijab += R_ijab / D_ijab
     */
    // L_IJAB /= D_IJAB
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "R <OO|VV>");
    dpd_buf4_dirprd(&D, &R);
    dpd_buf4_close(&D);
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
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

    // L_IJAB += R_IJAB / D_IJAB
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "R <oo|vv>");
    dpd_buf4_dirprd(&D, &R);
    dpd_buf4_close(&D);
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    dpd_buf4_add(&L, &R, 1.0);
    dpd_buf4_close(&R);
    dpd_buf4_close(&L);

    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

}} // Namespaces


