#include "dcft.h"
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include "defines.h"

namespace psi{ namespace dcft{

/**
 * Uses the intermediates to compute the energy
 */
void
DCFTSolver::compute_dcft_energy()
{
    dcft_timer_on("DCFTSolver::compute_dcft_energy()");

    dpdbuf4 L, G, I;
    double eGaa, eGab, eGbb, eIaa, eIab, eIbb;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // E += 1/4 L_IJAB G_IJAB
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    eGaa = 0.25 * dpd_buf4_dot(&G, &L);
    dpd_buf4_close(&G);

    // E += 1/4 gbar_IJAB L_IJAB
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    eIaa = 0.25 * dpd_buf4_dot(&I, &L);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    // E += L_IjAb G_IjAb
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    eGab =  dpd_buf4_dot(&G, &L);
    dpd_buf4_close(&G);

    // E += gbar_IjAb L_IjAb
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    eIab = dpd_buf4_dot(&I, &L);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    // E += 1/4 L_ijab G_ijab
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    eGbb = 0.25 * dpd_buf4_dot(&G, &L);
    dpd_buf4_close(&G);


    // E += 1/4 gbar_ijab L_ijab
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    eIbb = 0.25 * dpd_buf4_dot(&I, &L);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    psio_->close(PSIF_LIBTRANS_DPD, 1);

#if PRINT_ENERGY_COMPONENTS
    fprintf(outfile, "\tAA G Energy = %20.12f\n", eGaa);
    fprintf(outfile, "\tAB G Energy = %20.12f\n", eGab);
    fprintf(outfile, "\tBB G Energy = %20.12f\n", eGbb);
    fprintf(outfile, "\tAA I Energy = %20.12f\n", eIaa);
    fprintf(outfile, "\tAB I Energy = %20.12f\n", eIab);
    fprintf(outfile, "\tBB I Energy = %20.12f\n", eIbb);
    fprintf(outfile, "\tTotal G Energy = %20.12f\n", eGaa + eGab + eGbb);
    fprintf(outfile, "\tTotal I Energy = %20.12f\n", eIaa + eIab + eIbb);
#endif

    lambda_energy_ = eGaa + eGab + eGbb + eIaa + eIab + eIbb;

    dcft_timer_off("DCFTSolver::compute_dcft_energy()");
}

}} // Namespaces



