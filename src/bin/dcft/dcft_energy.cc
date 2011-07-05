#include "dcft.h"
#include "psifiles.h"
#include <libtrans/integraltransform.h>
#include "defines.h"

namespace psi{ namespace dcft{

/**
 * Uses the intermediates to compute the energy
 */
void
DCFTSolver::compute_dcft_energy()
{
    dpdbuf4 L, G, T, A, I;
    double eGaa, eGab, eGbb, eAaa, eAab, eAbb;
    double eTaa, eTab, eTbb, eIaa, eIab, eIbb;
    old_total_energy_ = new_total_energy_;
    new_total_energy_ = scf_energy_;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // E += 1/4 L_IJAB G_IJAB
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    eGaa = 0.25 * dpd_buf4_dot(&G, &L);
    dpd_buf4_close(&G);
#if !REFACTORED
    if(options_.get_bool("IGNORE_TAU")){
        eTaa = 0.0;
    }else{
        // E += 1/8 L_IJAB T_IJAB
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "T <OO|VV>");
        eTaa = (1.0/8.0) * dpd_buf4_dot(&L, &T);
        dpd_buf4_close(&T);
    }
#endif
    // E += 1/4 gbar_IJAB L_IJAB
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    eIaa = 0.25 * dpd_buf4_dot(&I, &L);
    dpd_buf4_close(&I);
    // E -= 1/4 L_IJAB^2 D_IJAB
    dpd_buf4_init(&A, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "A <OO|VV>");
    eAaa = -0.25 * dpd_buf4_dot(&A, &L);
    dpd_buf4_close(&A);
    dpd_buf4_close(&L);

    // E += L_IjAb G_IjAb
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    eGab =  dpd_buf4_dot(&G, &L);
    dpd_buf4_close(&G);

#if !REFACTORED
    if(options_.get_bool("IGNORE_TAU")){
        eTab = 0.0;
    }else{
        // E += 1/2 L_IjAb T_IjAb
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "T <Oo|Vv>");
        eTab = 0.5 * dpd_buf4_dot(&L, &T);
        dpd_buf4_close(&T);
    }    
#endif

    // E += gbar_IjAb L_IjAb
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    eIab = dpd_buf4_dot(&I, &L);
    dpd_buf4_close(&I);
    // E -= L_IjAb^2 D_IjAb
    dpd_buf4_init(&A, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "A <Oo|Vv>");
    eAab = -dpd_buf4_dot(&A, &L);
    dpd_buf4_close(&A);
    dpd_buf4_close(&L);

    // E += 1/4 L_ijab G_ijab
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    eGbb = 0.25 * dpd_buf4_dot(&G, &L);
    dpd_buf4_close(&G);

#if !REFACTORED
    if(options_.get_bool("IGNORE_TAU")){
        eTbb = 0.0;
    }else{
        // E += 1/8 L_ijab T_ijab
        dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "T <oo|vv>");
        eTbb = (1.0/8.0) * dpd_buf4_dot(&L, &T);
        dpd_buf4_close(&T);
    }
#endif

    // E += 1/4 gbar_ijab L_ijab
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    eIbb = 0.25 * dpd_buf4_dot(&I, &L);
    dpd_buf4_close(&I);
    // E -= 1/4 L_ijab^2 D_ijab
    dpd_buf4_init(&A, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "A <oo|vv>");
    eAbb = -0.25 * dpd_buf4_dot(&A, &L);
    dpd_buf4_close(&A);
    dpd_buf4_close(&L);
    psio_->close(PSIF_LIBTRANS_DPD, 1);

#if PRINT_ENERGY_COMPONENTS
    fprintf(outfile, "\tAA G Energy = %20.12f\n", eGaa);
    fprintf(outfile, "\tAB G Energy = %20.12f\n", eGab);
    fprintf(outfile, "\tBB G Energy = %20.12f\n", eGbb);
    fprintf(outfile, "\tAA I Energy = %20.12f\n", eIaa);
    fprintf(outfile, "\tAB I Energy = %20.12f\n", eIab);
    fprintf(outfile, "\tBB I Energy = %20.12f\n", eIbb);
    fprintf(outfile, "\tAA T Energy = %20.12f\n", eTaa);
    fprintf(outfile, "\tAB T Energy = %20.12f\n", eTab);
    fprintf(outfile, "\tBB T Energy = %20.12f\n", eTbb);
    fprintf(outfile, "\tAA A Energy = %20.12f\n", eAaa);
    fprintf(outfile, "\tAB A Energy = %20.12f\n", eAab);
    fprintf(outfile, "\tBB A Energy = %20.12f\n", eAbb);
    fprintf(outfile, "\tTotal G Energy = %20.12f\n", eGaa + eGab + eGbb);
    fprintf(outfile, "\tTotal I Energy = %20.12f\n", eIaa + eIab + eIbb);
    fprintf(outfile, "\tTotal T Energy = %20.12f\n", eTaa + eTab + eTbb);
    fprintf(outfile, "\tTotal A Energy = %20.12f\n", eAaa + eAab + eAbb);
#endif

#if !REFACTORED
    new_total_energy_ += eGaa + eGab + eGbb + eAaa + eAab + eAbb
                     + eTaa + eTab + eTbb + eIaa + eIab + eIbb;
    fprintf(outfile, "\n  !LAMBDA ENERGY = %20.15f \n", new_total_energy_ - scf_energy_);

#else
    new_total_energy_ += eGaa + eGab + eGbb + eIaa + eIab + eIbb;
    fprintf(outfile, "\n  !LAMBDA ENERGY = %20.15f \n", new_total_energy_ - scf_energy_);
#endif
}

}} // Namespaces



