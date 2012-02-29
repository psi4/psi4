#include <libtrans/integraltransform.h>
#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

void
DCFTSolver::compute_gradient()
{
    bool orbResponseDone    = false;
    bool lambdaResponseDone = false;

    // Print out the header
    fprintf(outfile,"\n\n\t\t*******************************************\n");
    fprintf(outfile,    "\t\t*       DCFT Analytic Gradients Code      *\n");
    fprintf(outfile,    "\t\t*   by A.Yu. Sokolov and A.C. Simmonett   *\n");
    fprintf(outfile,    "\t\t*******************************************\n\n");

    orbital_response_guess();
//    cumulant_response_guess();

    int cycle = 0;

}

void
DCFTSolver::orbital_response_guess()
{

    dpdbuf4 L;
    dpdfile2 T;

    // Copy the converged cumulant as a guess for the cumulant response

    // Z_IJAB = L_IJAB
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_copy(&L,PSIF_DCFT_DPD,"Z <OO|VV>");
    dpd_buf4_close(&L);

    // Z_IjAb = L_IjAb
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_copy(&L,PSIF_DCFT_DPD,"Z <Oo|Vv>");
    dpd_buf4_close(&L);

    // Z_ijab = L_ijab
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    dpd_buf4_copy(&L,PSIF_DCFT_DPD,"Z <oo|vv>");
    dpd_buf4_close(&L);

    // Copy the latest tau as a guess for relaxed tau

    dpd_file2_init(&T, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('O'), _ints->DPD_ID('O'), "Tau <O|O>");
    dpd_file2_copy(&T,PSIF_DCFT_DPD,"Tau~ <O|O>");
    dpd_file2_close(&T);

    dpd_file2_init(&T, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('o'), _ints->DPD_ID('o'), "Tau <o|o>");
    dpd_file2_copy(&T,PSIF_DCFT_DPD,"Tau~ <o|o>");
    dpd_file2_close(&T);

    dpd_file2_init(&T, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('V'), _ints->DPD_ID('V'), "Tau <V|V>");
    dpd_file2_copy(&T,PSIF_DCFT_DPD,"Tau~ <V|V>");
    dpd_file2_close(&T);

    dpd_file2_init(&T, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('v'), _ints->DPD_ID('v'), "Tau <v|v>");
    dpd_file2_copy(&T,PSIF_DCFT_DPD,"Tau~ <v|v>");
    dpd_file2_close(&T);

    // Compute the generalized densities for the MO Lagrangian
    compute_density();

}

void
DCFTSolver::compute_density()
{


}


}} //End namespaces
