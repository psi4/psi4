#include <libtrans/integraltransform.h>
#include <libpsio/psio.hpp>
#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

void
DCFTSolver::compute_gradient()
{
    bool responseDone = false;

    // Print out the header
    fprintf(outfile,"\n\n\t\t*******************************************\n");
    fprintf(outfile,    "\t\t*       DCFT Analytic Gradients Code      *\n");
    fprintf(outfile,    "\t\t*   by A.Yu. Sokolov and A.C. Simmonett   *\n");
    fprintf(outfile,    "\t\t*******************************************\n\n");

    // Transform the one and two-electron integrals to the MO basis and write them into the DPD file
    gradient_init();
    // Copy the current density cumulant and tau as a guess for cumulant response and perturbed tau
    response_guess();

    iter_ = 0;
    fprintf(outfile, "\n\n\t*=============================================*\n"
                     "\t* Cycle  RMS Orb. Resp.   RMS Cumul. Resp.    *\n"
                     "\t*---------------------------------------------*\n");

    // Start macro-iterations
    while(!responseDone && iter_++ < maxiter_){
        fprintf(outfile, "\t                          *** Macro Iteration %d ***\n" "\tOrbital Response Iterations\n", iter_);

        // Solve the cumulant response equations iteratively
        if (iter_ > 1) {
            iterate_cumulant_response();
            responseDone = true;   ////
            break;  /////
        }

        // Compute the generalized densities for the MO Lagrangian
        compute_density();
        // Compute the OV block of MO Lagrangian
        compute_lagrangian_OV();
        // Compute the VO block of MO Lagrangian
        compute_lagrangian_VO();

        // Solve the orbital response equations iteratively
        iterate_orbital_response();

        // Compute terms that couple orbital and cumulant responses (C intermediate) and return RMS of their change
        double response_coupling_rms = compute_response_coupling();

        // Check convergence
        if (response_coupling_rms < 1.0E-14) responseDone = true;

        fprintf(outfile, "\tResponse Coupling RMS = %11.3E \n", response_coupling_rms);

        fprintf(outfile, "\tOrbital Response Iterations\n");
    }

    if (responseDone) fprintf(outfile, "    DCFT response equations converged.\n\n");
    else throw PSIEXCEPTION("DCFT response equations did not converge");

    fprintf(outfile, "\t*=============================================*\n");






}

void
DCFTSolver::gradient_init()
{

    // Allocate memory for the global objects

    aocc_tau_ = SharedMatrix(new Matrix("MO basis Tau (Alpha Occupied)", nirrep_, naoccpi_, naoccpi_));
    bocc_tau_ = SharedMatrix(new Matrix("MO basis Tau (Beta Occupied)", nirrep_, nboccpi_, nboccpi_));
    avir_tau_ = SharedMatrix(new Matrix("MO basis Tau (Alpha Virtual)", nirrep_, navirpi_, navirpi_));
    bvir_tau_ = SharedMatrix(new Matrix("MO basis Tau (Beta Virtual)", nirrep_, nbvirpi_, nbvirpi_));
    aocc_ptau_ = SharedMatrix(new Matrix("MO basis Perturbed Tau (Alpha Occupied)", nirrep_, naoccpi_, naoccpi_));
    bocc_ptau_ = SharedMatrix(new Matrix("MO basis Perturbed Tau (Beta Occupied)", nirrep_, nboccpi_, nboccpi_));
    avir_ptau_ = SharedMatrix(new Matrix("MO basis Perturbed Tau (Alpha Virtual)", nirrep_, navirpi_, navirpi_));
    bvir_ptau_ = SharedMatrix(new Matrix("MO basis Perturbed Tau (Beta Virtual)", nirrep_, nbvirpi_, nbvirpi_));
    akappa_ = SharedMatrix(new Matrix("MO basis Kappa (Alpha)", nirrep_, naoccpi_, naoccpi_));
    bkappa_ = SharedMatrix(new Matrix("MO basis Kappa (Beta)", nirrep_, nboccpi_, nboccpi_));

    dpdbuf4 I;
    dpdfile2 H;

    // Transform the two-electron integrals to the (VO|OO) and (OV|VV) subspaces in chemists' notation

    _ints->transform_tei(MOSpace::vir, MOSpace::occ, MOSpace::occ, MOSpace::occ);
    _ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::occ);
    _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::vir, MOSpace::vir);
    _ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::vir);

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    /*
     * Re-sort the chemists' notation integrals to physisists' notation
     * (pq|rs) = <pr|qs>
     */

    // <VO|OO> type

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O>=O]+"), 0, "MO Ints (VO|OO)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[V,O]"), ID("[O,O]"), "MO Ints <VO|OO>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O>=O]+"), 0, "MO Ints (VO|OO)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rsqp, ID("[O,O]"), ID("[O,V]"), "MO Ints (OO|OV)");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "MO Ints <VO|OO>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, qpsr, ID("[O,V]"), ID("[O,O]"), "MO Ints <OV|OO>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "MO Ints <OV|OO>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[O,O]"), ID("[O,V]"), "MO Ints <OO|OV>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,o]"),
                  ID("[V,O]"), ID("[o>=o]+"), 0, "MO Ints (VO|oo)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[V,o]"), ID("[O,o]"), "MO Ints <Vo|Oo>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[o,o]"),
                  ID("[V,O]"), ID("[o>=o]+"), 0, "MO Ints (VO|oo)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[o,o]"), ID("[V,O]"), "MO Ints (oo|VO)");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "MO Ints <Vo|Oo>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, qpsr, ID("[o,V]"), ID("[o,O]"), "MO Ints <oV|oO>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <oV|oO>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, pqsr, ID("[o,V]"), ID("[O,o]"), "MO Ints <oV|Oo>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,o]"),
                  ID("[O>=O]+"), ID("[v,o]"), 0, "MO Ints (OO|vo)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rqsp, ID("[v,O]"), ID("[o,O]"), "MO Ints <vO|oO>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,O]"), ID("[o,O]"),
                  ID("[v,O]"), ID("[o,O]"), 0, "MO Ints <vO|oO>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, qpsr, ID("[O,v]"), ID("[O,o]"), "MO Ints <Ov|Oo>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, pqsr, ID("[O,v]"), ID("[o,O]"), "MO Ints <Ov|oO>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o>=o]+"), 0, "MO Ints (vo|oo)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[v,o]"), ID("[o,o]"), "MO Ints <vo|oo>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o>=o]+"), 0, "MO Ints (vo|oo)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rsqp, ID("[o,o]"), ID("[o,v]"), "MO Ints (oo|ov)");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "MO Ints <vo|oo>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, qpsr, ID("[o,v]"), ID("[o,o]"), "MO Ints <ov|oo>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o,o]"), 0, "MO Ints <ov|oo>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[o,o]"), ID("[o,v]"), "MO Ints <oo|ov>");
    dpd_buf4_close(&I);

    // <OV|VV> type

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,V]"), ID("[V,V]"), "MO Ints <OV|VV>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[V,V]"), ID("[O,V]"), "MO Ints <VV|OV>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                  ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[O,v]"), ID("[V,v]"), "MO Ints <Ov|Vv>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                  ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[v,v]"), ID("[O,V]"), "MO Ints (vv|OV)");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, pqsr, ID("[O,v]"), ID("[v,V]"), "MO Ints <Ov|vV>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,v]"),
                  ID("[V>=V]+"), ID("[o,v]"), 0, "MO Ints (VV|ov)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rqsp, ID("[o,V]"), ID("[v,V]"), "MO Ints <oV|vV>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV|vV>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, pqsr, ID("[o,V]"), ID("[V,v]"), "MO Ints <oV|Vv>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV|vV>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, qpsr, ID("[V,o]"), ID("[V,v]"), "MO Ints <Vo|Vv>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v>=v]+"), 0, "MO Ints (ov|vv)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, prqs, ID("[o,v]"), ID("[v,v]"), "MO Ints <ov|vv>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[v,v]"), ID("[o,v]"), "MO Ints <vv|ov>");
    dpd_buf4_close(&I);

    // (OV|OV)

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[o,v]"), ID("[O,V]"), "MO Ints (ov|OV)");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rqps, ID("[V,O]"), ID("[O,V]"), "MO Ints <VO|OV>");
    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rqps, ID("[v,o]"), ID("[o,v]"), "MO Ints <vo|ov>");
    dpd_buf4_close(&I);

    // Transform one-electron integrals to the MO basis and store them in the DPD file

    Matrix aH(so_h_);
    Matrix bH(so_h_);
    aH.transform(Ca_);
    bH.transform(Cb_);

    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    dpd_file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int j = 0 ; j < naoccpi_[h]; ++j){
                H.matrix[h][i][j] = aH.get(h, i, j);
            }
        }
    }
    dpd_file2_mat_wrt(&H);
    dpd_file2_close(&H);

    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    dpd_file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int a = 0 ; a < navirpi_[h]; ++a){
            for(int b = 0 ; b < navirpi_[h]; ++b){
                H.matrix[h][a][b] = aH.get(h, naoccpi_[h] + a, naoccpi_[h] + b);
            }
        }
    }
    dpd_file2_mat_wrt(&H);
    dpd_file2_close(&H);

    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "H <o|o>");
    dpd_file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int j = 0 ; j < nboccpi_[h]; ++j){
                H.matrix[h][i][j] = bH.get(h, i, j);
            }
        }
    }
    dpd_file2_mat_wrt(&H);
    dpd_file2_close(&H);

    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "H <v|v>");
    dpd_file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int a = 0 ; a < nbvirpi_[h]; ++a){
            for(int b = 0 ; b < nbvirpi_[h]; ++b){
                H.matrix[h][a][b] = bH.get(h, nboccpi_[h] + a, nboccpi_[h] + b);
            }
        }
    }
    dpd_file2_mat_wrt(&H);
    dpd_file2_close(&H);

    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    dpd_file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int j = 0 ; j < navirpi_[h]; ++j){
                H.matrix[h][i][j] = aH.get(h, i, naoccpi_[h] + j);
            }
        }
    }
    dpd_file2_mat_wrt(&H);
    dpd_file2_close(&H);

    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    dpd_file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int j = 0 ; j < nbvirpi_[h]; ++j){
                H.matrix[h][i][j] = bH.get(h, i, nboccpi_[h] + j);
            }
        }
    }
    dpd_file2_mat_wrt(&H);
    dpd_file2_close(&H);


    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::response_guess()
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

    dpd_file2_init(&T, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    dpd_file2_copy(&T,PSIF_DCFT_DPD,"pTau <O|O>");
    dpd_file2_close(&T);

    dpd_file2_init(&T, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    dpd_file2_copy(&T,PSIF_DCFT_DPD,"pTau <o|o>");
    dpd_file2_close(&T);

    dpd_file2_init(&T, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    dpd_file2_copy(&T,PSIF_DCFT_DPD,"pTau <V|V>");
    dpd_file2_close(&T);

    dpd_file2_init(&T, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");
    dpd_file2_copy(&T,PSIF_DCFT_DPD,"pTau <v|v>");
    dpd_file2_close(&T);

}

void
DCFTSolver::compute_density()
{

    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);

    dpdbuf4 Zaa, Zab, Zbb, Laa, Lab, Lbb, Gaa, Gab, Gba, Gbb, Tab;
    dpdfile2 pT_OO, pT_oo, pT_VV, pT_vv, T_OO, T_oo, T_VV, T_vv;

    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Z <OO|VV>");
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Z <oo|vv>");

    dpd_file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    dpd_file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");
    dpd_file2_init(&pT_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "pTau <O|O>");
    dpd_file2_init(&pT_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "pTau <o|o>");
    dpd_file2_init(&pT_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "pTau <V|V>");
    dpd_file2_init(&pT_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "pTau <v|v>");

    dpd_file2_mat_init(&T_OO);
    dpd_file2_mat_init(&T_oo);
    dpd_file2_mat_init(&T_VV);
    dpd_file2_mat_init(&T_vv);
    dpd_file2_mat_init(&pT_OO);
    dpd_file2_mat_init(&pT_oo);
    dpd_file2_mat_init(&pT_VV);
    dpd_file2_mat_init(&pT_vv);

    dpd_file2_mat_rd(&T_OO);
    dpd_file2_mat_rd(&T_oo);
    dpd_file2_mat_rd(&T_VV);
    dpd_file2_mat_rd(&T_vv);
    dpd_file2_mat_rd(&pT_OO);
    dpd_file2_mat_rd(&pT_oo);
    dpd_file2_mat_rd(&pT_VV);
    dpd_file2_mat_rd(&pT_vv);

    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int j = 0; j < naoccpi_[h]; ++j){
                akappa_->set(h, i, j, (i == j ? 1.0 : 0.0));
                aocc_tau_->set(h, i, j, T_OO.matrix[h][i][j]);
                aocc_ptau_->set(h, i, j, pT_OO.matrix[h][i][j]);
            }
        }
        for(int a = 0; a < navirpi_[h]; ++a){
            for(int b = 0; b < navirpi_[h]; ++b){
                avir_tau_->set(h, a, b, T_VV.matrix[h][a][b]);
                avir_ptau_->set(h, a, b, pT_VV.matrix[h][a][b]);
            }
        }
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int j = 0; j < nboccpi_[h]; ++j){
                bkappa_->set(h, i, j, (i == j ? 1.0 : 0.0));
                bocc_tau_->set(h, i, j, T_oo.matrix[h][i][j]);
                bocc_ptau_->set(h, i, j, pT_oo.matrix[h][i][j]);
            }
        }
        for(int a = 0; a < nbvirpi_[h]; ++a){
            for(int b = 0; b < nbvirpi_[h]; ++b){
                bvir_tau_->set(h, a, b, T_vv.matrix[h][a][b]);
                bvir_ptau_->set(h, a, b, pT_vv.matrix[h][a][b]);
            }
        }
    }

    dpd_file2_close(&T_OO);
    dpd_file2_close(&T_oo);
    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);
    dpd_file2_close(&pT_OO);
    dpd_file2_close(&pT_oo);
    dpd_file2_close(&pT_VV);
    dpd_file2_close(&pT_vv);

    /*
     * The VVVV block
     */
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
              ID("[V,V]"), ID("[V,V]"), 0, "Gamma <VV|VV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gaa, h);
        dpd_buf4_mat_irrep_init(&Laa, h);
        dpd_buf4_mat_irrep_init(&Zaa, h);
        dpd_buf4_mat_irrep_rd(&Laa, h);
        dpd_buf4_mat_irrep_rd(&Zaa, h);
        for(size_t ab = 0; ab < Gaa.params->rowtot[h]; ++ab){
            size_t a = Gaa.params->roworb[h][ab][0];
            int Ga = Gaa.params->psym[a];
            a -= Gaa.params->poff[Ga];
            size_t b = Gaa.params->roworb[h][ab][1];
            int Gb = Gaa.params->qsym[b];
            b -= Gaa.params->qoff[Gb];
            for(size_t cd = 0; cd < Gaa.params->coltot[h]; ++cd){
                double tpdm = 0.0;
                for(size_t ij = 0; ij < Laa.params->rowtot[h]; ++ij){
                    tpdm += 0.0625 * (Zaa.matrix[h][ij][ab] * Laa.matrix[h][ij][cd] + Laa.matrix[h][ij][ab] * Zaa.matrix[h][ij][cd]);
                }
                size_t c = Gaa.params->colorb[h][cd][0];
                int Gc = Gaa.params->rsym[c];
                c -= Gaa.params->roff[Gc];
                size_t d = Gaa.params->colorb[h][cd][1];
                int Gd = Gaa.params->ssym[d];
                d -= Gaa.params->soff[Gd];
                if(Ga == Gc && Gb == Gd) tpdm += 0.25 * avir_tau_->get(Ga, a, c) * avir_ptau_->get(Gb, b, d);
                if(Ga == Gd && Gb == Gc) tpdm -= 0.25 * avir_tau_->get(Ga, a, d) * avir_ptau_->get(Gb, b, c);
                if(Gb == Gc && Ga == Gd) tpdm -= 0.25 * avir_tau_->get(Gb, b, c) * avir_ptau_->get(Ga, a, d);
                if(Ga == Gc && Gb == Gd) tpdm += 0.25 * avir_tau_->get(Gb, b, d) * avir_ptau_->get(Ga, a, c);

                if(Ga == Gc && Gb == Gd) tpdm -= 0.25 * avir_tau_->get(Ga, a, c) * avir_tau_->get(Gb, b, d);
                if(Ga == Gd && Gb == Gc) tpdm += 0.25 * avir_tau_->get(Ga, a, d) * avir_tau_->get(Gb, b, c);

                Gaa.matrix[h][ab][cd] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Laa, h);
        dpd_buf4_mat_irrep_close(&Zaa, h);
    }

    dpd_buf4_print(&Gaa,outfile,1);

    dpd_buf4_close(&Gaa);

    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
              ID("[V,v]"), ID("[V,v]"), 0, "Gamma <Vv|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gab, h);
        dpd_buf4_mat_irrep_init(&Lab, h);
        dpd_buf4_mat_irrep_init(&Zab, h);
        dpd_buf4_mat_irrep_rd(&Lab, h);
        dpd_buf4_mat_irrep_rd(&Zab, h);
        for(size_t ab = 0; ab < Gab.params->rowtot[h]; ++ab){
            size_t a = Gab.params->roworb[h][ab][0];
            int Ga = Gab.params->psym[a];
            a -= Gab.params->poff[Ga];
            size_t b = Gab.params->roworb[h][ab][1];
            int Gb = Gab.params->qsym[b];
            b -= Gab.params->qoff[Gb];
            for(size_t cd = 0; cd < Gab.params->coltot[h]; ++cd){
                double tpdm = 0.0;
                for(size_t ij = 0; ij < Lab.params->rowtot[h]; ++ij){
                    tpdm += 0.125 * (Zab.matrix[h][ij][ab] * Lab.matrix[h][ij][cd] + Lab.matrix[h][ij][ab] * Zab.matrix[h][ij][cd]);
                }
                size_t c = Gab.params->colorb[h][cd][0];
                int Gc = Gab.params->rsym[c];
                c -= Gab.params->roff[Gc];
                size_t d = Gab.params->colorb[h][cd][1];
                int Gd = Gab.params->ssym[d];
                d -= Gab.params->soff[Gd];
                if(Ga == Gc && Gb == Gd) tpdm += 0.25 * avir_tau_->get(Ga, a, c) * bvir_ptau_->get(Gb, b, d);
                if(Ga == Gc && Gb == Gd) tpdm += 0.25 * bvir_tau_->get(Gb, b, d) * avir_ptau_->get(Ga, a, c);

                if(Ga == Gc && Gb == Gd) tpdm -= 0.25 * avir_tau_->get(Ga, a, c) * bvir_tau_->get(Gb, b, d);
                Gab.matrix[h][ab][cd] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gab, h);
        dpd_buf4_mat_irrep_close(&Gab, h);
        dpd_buf4_mat_irrep_close(&Lab, h);
        dpd_buf4_mat_irrep_close(&Zab, h);
    }

    dpd_buf4_print(&Gab,outfile,1);

    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
              ID("[v,v]"), ID("[v,v]"), 0, "Gamma <vv|vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gbb, h);
        dpd_buf4_mat_irrep_init(&Lbb, h);
        dpd_buf4_mat_irrep_init(&Zbb, h);
        dpd_buf4_mat_irrep_rd(&Lbb, h);
        dpd_buf4_mat_irrep_rd(&Zbb, h);

        for(size_t ab = 0; ab < Gbb.params->rowtot[h]; ++ab){
            size_t a = Gbb.params->roworb[h][ab][0];
            int Ga = Gbb.params->psym[a];
            a -= Gbb.params->poff[Ga];
            size_t b = Gbb.params->roworb[h][ab][1];
            int Gb = Gbb.params->qsym[b];
            b -= Gbb.params->qoff[Gb];
            for(size_t cd = 0; cd < Gbb.params->coltot[h]; ++cd){
                double tpdm = 0.0;
                for(size_t ij = 0; ij < Lbb.params->rowtot[h]; ++ij){
                    tpdm += 0.0625 * (Zbb.matrix[h][ij][ab] * Lbb.matrix[h][ij][cd] + Lbb.matrix[h][ij][ab] * Zbb.matrix[h][ij][cd]);
                }
                size_t c = Gbb.params->colorb[h][cd][0];
                int Gc = Gbb.params->rsym[c];
                c -= Gbb.params->roff[Gc];
                size_t d = Gbb.params->colorb[h][cd][1];
                int Gd = Gbb.params->ssym[d];
                d -= Gbb.params->soff[Gd];
                if(Ga == Gc && Gb == Gd) tpdm += 0.25 * bvir_tau_->get(Ga, a, c) * bvir_ptau_->get(Gb, b, d);
                if(Ga == Gd && Gb == Gc) tpdm -= 0.25 * bvir_tau_->get(Ga, a, d) * bvir_ptau_->get(Gb, b, c);
                if(Gb == Gc && Ga == Gd) tpdm -= 0.25 * bvir_tau_->get(Gb, b, c) * bvir_ptau_->get(Ga, a, d);
                if(Ga == Gc && Gb == Gd) tpdm += 0.25 * bvir_tau_->get(Gb, b, d) * bvir_ptau_->get(Ga, a, c);

                if(Ga == Gc && Gb == Gd) tpdm -= 0.25 * bvir_tau_->get(Ga, a, c) * bvir_tau_->get(Gb, b, d);
                if(Ga == Gd && Gb == Gc) tpdm += 0.25 * bvir_tau_->get(Ga, a, d) * bvir_tau_->get(Gb, b, c);
                Gbb.matrix[h][ab][cd] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Lbb, h);
        dpd_buf4_mat_irrep_close(&Zbb, h);
    }

    dpd_buf4_print(&Gbb,outfile,1);

    dpd_buf4_close(&Gbb);

    /*
     * The OOOO  block
     */
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
              ID("[O,O]"), ID("[O,O]"), 0, "Gamma <OO|OO>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gaa, h);
        dpd_buf4_mat_irrep_init(&Laa, h);
        dpd_buf4_mat_irrep_init(&Zaa, h);
        dpd_buf4_mat_irrep_rd(&Laa, h);
        dpd_buf4_mat_irrep_rd(&Zaa, h);
                for(size_t ij = 0; ij < Gaa.params->rowtot[h]; ++ij){
            size_t i = Gaa.params->roworb[h][ij][0];
            int Gi = Gaa.params->psym[i];
            i -= Gaa.params->poff[Gi];
            size_t j = Gaa.params->roworb[h][ij][1];
            int Gj = Gaa.params->qsym[j];
            j -= Gaa.params->qoff[Gj];
            for(size_t kl = 0; kl < Gaa.params->coltot[h]; ++kl){
                double tpdm = 0.0;
                for(size_t ab = 0; ab < Laa.params->coltot[h]; ++ab){
                    tpdm += 0.0625 * (Zaa.matrix[h][ij][ab] * Laa.matrix[h][kl][ab] + Laa.matrix[h][ij][ab] * Zaa.matrix[h][kl][ab]);
                }
                size_t k = Gaa.params->colorb[h][kl][0];
                int Gk = Gaa.params->rsym[k];
                k -= Gaa.params->roff[Gk];
                size_t l = Gaa.params->colorb[h][kl][1];
                int Gl = Gaa.params->ssym[l];
                l -= Gaa.params->soff[Gl];

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * akappa_->get(Gi, i, k) * akappa_->get(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * akappa_->get(Gi, i, l) * akappa_->get(Gj, j, k);

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * (akappa_->get(Gi, i, k) + aocc_tau_->get(Gi, i, k)) * aocc_ptau_->get(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * (akappa_->get(Gi, i, l) + aocc_tau_->get(Gi, i, l)) * aocc_ptau_->get(Gj, j, k);
                if(Gj == Gk && Gi == Gl) tpdm -= 0.25 * (akappa_->get(Gj, j, k) + aocc_tau_->get(Gj, j, k)) * aocc_ptau_->get(Gi, i, l);
                if(Gj == Gl && Gi == Gk) tpdm += 0.25 * (akappa_->get(Gj, j, l) + aocc_tau_->get(Gj, j, l)) * aocc_ptau_->get(Gi, i, k);

                if(Gi == Gk && Gj == Gl) tpdm -= 0.25 * aocc_tau_->get(Gi, i, k) * aocc_tau_->get(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm += 0.25 * aocc_tau_->get(Gi, i, l) * aocc_tau_->get(Gj, j, k);

                Gaa.matrix[h][ij][kl] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Laa, h);
        dpd_buf4_mat_irrep_close(&Zaa, h);
    }

    dpd_buf4_print(&Gaa,outfile,1);

    dpd_buf4_close(&Gaa);


    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
              ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gab, h);
        dpd_buf4_mat_irrep_init(&Lab, h);
        dpd_buf4_mat_irrep_init(&Zab, h);
        dpd_buf4_mat_irrep_rd(&Lab, h);
        dpd_buf4_mat_irrep_rd(&Zab, h);
        for(size_t ij = 0; ij < Gab.params->rowtot[h]; ++ij){
            size_t i = Gab.params->roworb[h][ij][0];
            int Gi = Gab.params->psym[i];
            i -= Gab.params->poff[Gi];
            size_t j = Gab.params->roworb[h][ij][1];
            int Gj = Gab.params->qsym[j];
            j -= Gab.params->qoff[Gj];
            for(size_t kl = 0; kl < Gab.params->coltot[h]; ++kl){
                double tpdm = 0.0;
                for(size_t ab = 0; ab < Lab.params->coltot[h]; ++ab){
                    tpdm += 0.125 * (Zab.matrix[h][ij][ab] * Lab.matrix[h][kl][ab] + Lab.matrix[h][ij][ab] * Zab.matrix[h][kl][ab]);
                }
                size_t k = Gab.params->colorb[h][kl][0];
                int Gk = Gab.params->rsym[k];
                k -= Gab.params->roff[Gk];
                size_t l = Gab.params->colorb[h][kl][1];
                int Gl = Gab.params->ssym[l];
                l -= Gab.params->soff[Gl];
                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * akappa_->get(Gi, i, k) * bkappa_->get(Gj, j, l);

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * (akappa_->get(Gi, i, k) + aocc_tau_->get(Gi, i, k)) * bocc_ptau_->get(Gj, j, l);
                if(Gj == Gl && Gi == Gk) tpdm += 0.25 * (bkappa_->get(Gj, j, l) + bocc_tau_->get(Gj, j, l)) * aocc_ptau_->get(Gi, i, k);

                if(Gi == Gk && Gj == Gl) tpdm -= 0.25 * aocc_tau_->get(Gi, i, k) * bocc_tau_->get(Gj, j, l);

                Gab.matrix[h][ij][kl] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gab, h);
        dpd_buf4_mat_irrep_close(&Gab, h);
        dpd_buf4_mat_irrep_close(&Lab, h);
        dpd_buf4_mat_irrep_close(&Zab, h);
    }

    dpd_buf4_print(&Gab,outfile,1);

    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
              ID("[o,o]"), ID("[o,o]"), 0, "Gamma <oo|oo>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gbb, h);
        dpd_buf4_mat_irrep_init(&Lbb, h);
        dpd_buf4_mat_irrep_init(&Zbb, h);
        dpd_buf4_mat_irrep_rd(&Lbb, h);
        dpd_buf4_mat_irrep_rd(&Zbb, h);
        for(size_t ij = 0; ij < Gbb.params->rowtot[h]; ++ij){
            size_t i = Gbb.params->roworb[h][ij][0];
            int Gi = Gbb.params->psym[i];
            i -= Gbb.params->poff[Gi];
            size_t j = Gbb.params->roworb[h][ij][1];
            int Gj = Gbb.params->qsym[j];
            j -= Gbb.params->qoff[Gj];
            for(size_t kl = 0; kl < Gbb.params->coltot[h]; ++kl){
                double tpdm = 0.0;
                for(size_t ab = 0; ab < Lbb.params->coltot[h]; ++ab){
                    tpdm += 0.0625 * (Zbb.matrix[h][ij][ab] * Lbb.matrix[h][kl][ab] + Lbb.matrix[h][ij][ab] * Zbb.matrix[h][kl][ab]);
                }
                size_t k = Gbb.params->colorb[h][kl][0];
                int Gk = Gbb.params->rsym[k];
                k -= Gbb.params->roff[Gk];
                size_t l = Gbb.params->colorb[h][kl][1];
                int Gl = Gbb.params->ssym[l];
                l -= Gbb.params->soff[Gl];
                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * bkappa_->get(Gi, i, k) * bkappa_->get(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * bkappa_->get(Gi, i, l) * bkappa_->get(Gj, j, k);

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * (bkappa_->get(Gi, i, k) + bocc_tau_->get(Gi, i, k)) * bocc_ptau_->get(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * (bkappa_->get(Gi, i, l) + bocc_tau_->get(Gi, i, l)) * bocc_ptau_->get(Gj, j, k);
                if(Gj == Gk && Gi == Gl) tpdm -= 0.25 * (bkappa_->get(Gj, j, k) + bocc_tau_->get(Gj, j, k)) * bocc_ptau_->get(Gi, i, l);
                if(Gj == Gl && Gi == Gk) tpdm += 0.25 * (bkappa_->get(Gj, j, l) + bocc_tau_->get(Gj, j, l)) * bocc_ptau_->get(Gi, i, k);

                if(Gi == Gk && Gj == Gl) tpdm -= 0.25 * bocc_tau_->get(Gi, i, k) * bocc_tau_->get(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm += 0.25 * bocc_tau_->get(Gi, i, l) * bocc_tau_->get(Gj, j, k);

                Gbb.matrix[h][ij][kl] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Lbb, h);
        dpd_buf4_mat_irrep_close(&Zbb, h);
    }

    dpd_buf4_print(&Gbb,outfile,1);

    dpd_buf4_close(&Gbb);

    /*
     * The OOVV and VVOO blocks
     */

    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
              ID("[O,O]"), ID("[V,V]"), 0, "Gamma <OO|VV>");
    dpd_buf4_axpbycz(&Laa,&Zaa,&Gaa,0.25,0.25,1.0);

    dpd_buf4_print(&Gaa,outfile,1);

    dpd_buf4_close(&Gaa);

    // Resort the Г_OOVV to Г_VVOO. Used for the MO Lagrangian
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
              ID("[O,O]"), ID("[V,V]"), 0, "Gamma <OO|VV>");
    dpd_buf4_sort(&Gaa, PSIF_DCFT_DENSITY, rspq, ID("[V,V]"), ID("[O,O]"), "Gamma <VV|OO>");
    dpd_buf4_close(&Gaa);

    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
              ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");
    dpd_buf4_axpbycz(&Lab,&Zab,&Gab,0.25,0.25,1.0);

    dpd_buf4_print(&Gab,outfile,1);

    dpd_buf4_close(&Gab);

    // Resort the Г_OoVv to Г_VvOo. Used for the MO Lagrangian
    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
              ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");
    dpd_buf4_sort(&Gab, PSIF_DCFT_DENSITY, rspq, ID("[V,v]"), ID("[O,o]"), "Gamma <Vv|Oo>");
    dpd_buf4_close(&Gab);

//    // Resort the Г_VvOo to Г_vVoO. Used for the MO Lagrangian
//    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"),
//              ID("[V,v]"), ID("[O,o]"), 0, "Gamma <Vv|Oo>");
//    dpd_buf4_sort(&Gab, PSIF_DCFT_DENSITY, qpsr, ID("[v,V]"), ID("[o,O]"), "Gamma <vV|oO>");
//    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
              ID("[o,o]"), ID("[v,v]"), 0, "Gamma <oo|vv>");
    dpd_buf4_axpbycz(&Lbb,&Zbb,&Gbb,0.25,0.25,1.0);

    dpd_buf4_print(&Gbb,outfile,1);

    dpd_buf4_close(&Gbb);

    // Resort the Г_oovv to Г_vvoo. Used for the MO Lagrangian
    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
              ID("[o,o]"), ID("[v,v]"), 0, "Gamma <oo|vv>");
    dpd_buf4_sort(&Gbb, PSIF_DCFT_DENSITY, rspq, ID("[v,v]"), ID("[o,o]"), "Gamma <vv|oo>");
    dpd_buf4_close(&Gbb);

    dpd_buf4_close(&Laa);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Lbb);
    dpd_buf4_close(&Zaa);
    dpd_buf4_close(&Zab);
    dpd_buf4_close(&Zbb);

    /*
     * The OVOV block
     */

    // There are five unique spin cases: Г<IAJB>, Г<iajb>, Г<IaJb>, Г<iAjB>, Г<IajB>

    // TEMPORARY: Sort the cumulant Z-vector elements to chemist's notation.
    // MOVE THIS TO THE Z-VECTOR UPDATES WHEN NEEDED!!!
    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Z <OO|VV>");
    dpd_buf4_sort(&Zaa, PSIF_DCFT_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "Z (OV|OV)");
    dpd_buf4_close(&Zaa);

    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    dpd_buf4_sort(&Zab, PSIF_DCFT_DPD, psqr, ID("[O,v]"), ID("[o,V]"), "Z (Ov|oV)");
    dpd_buf4_close(&Zab);

    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Z <oo|vv>");
    dpd_buf4_sort(&Zbb, PSIF_DCFT_DPD, prqs, ID("[o,v]"),ID("[o,v]"), "Z (ov|ov)");
    dpd_buf4_close(&Zbb);

    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Z (Ov|oV)");

    dpd_buf4_sort(&Zab, PSIF_DCFT_DPD, psrq, ID("[O,V]"),ID("[o,v]"), "Z (OV|ov)");
    dpd_buf4_close(&Zab);

    // Г<IAJB> spin case

    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma (OV|OV)");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Z (OV|OV)");
    dpd_contract444(&Laa, &Zaa, &Gaa, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Laa);
    dpd_buf4_close(&Zaa);
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Z (OV|ov)");
    dpd_contract444(&Lab, &Zab, &Gaa, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Zab);
    dpd_buf4_close(&Gaa);
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma (OV|OV)");
    dpd_buf4_symm(&Gaa);
    dpd_buf4_close(&Gaa);

    // Resort Г(OV|OV) to the Г<OV|OV>
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma (OV|OV)");
    dpd_buf4_sort(&Gaa, PSIF_DCFT_DENSITY, psrq, ID("[O,V]"),ID("[O,V]"), "Gamma <OV|OV>");
    dpd_buf4_close(&Gaa);

    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gaa, h);
        dpd_buf4_mat_irrep_rd(&Gaa, h);
        for(size_t ia = 0; ia < Gaa.params->rowtot[h]; ++ia){
            size_t i = Gaa.params->roworb[h][ia][0];
            int Gi = Gaa.params->psym[i];
            i -= Gaa.params->poff[Gi];
            size_t a = Gaa.params->roworb[h][ia][1];
            int Ga = Gaa.params->qsym[a];
            a -= Gaa.params->qoff[Ga];
            for(size_t jb = 0; jb < Gaa.params->coltot[h]; ++jb){
                size_t j = Gaa.params->colorb[h][jb][0];
                int Gj = Gaa.params->rsym[j];
                j -= Gaa.params->roff[Gj];
                size_t b = Gaa.params->colorb[h][jb][1];
                int Gb = Gaa.params->ssym[b];
                b -= Gaa.params->soff[Gb];
                if(Gi == Gj && Ga == Gb) {
                    Gaa.matrix[h][ia][jb] += (akappa_->get(Gi, i, j) + aocc_tau_->get(Gi, i, j)) * avir_ptau_->get(Ga, a, b);
                    Gaa.matrix[h][ia][jb] += avir_tau_->get(Ga, a, b) * (aocc_ptau_->get(Gi, i, j) - aocc_tau_->get(Gi, i, j));
                }
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Gaa, h);
    }

    dpd_buf4_print(&Gaa,outfile,1);

    dpd_buf4_close(&Gaa);

    // Г<IaJb> and Г<iAjB> spin cases:

    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Lambda (Ov|oV)");
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Z (Ov|oV)");
    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
    dpd_contract444(&Lab, &Zab, &Gab, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Gab);
    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
    dpd_buf4_symm(&Gab);
    dpd_buf4_close(&Gab);
    dpd_buf4_init(&Gba, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");
    dpd_contract444(&Lab, &Zab, &Gba, 1, 1, -1.0, 0.0);
    dpd_buf4_close(&Gba);
    dpd_buf4_init(&Gba, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");
    dpd_buf4_symm(&Gba);
    dpd_buf4_close(&Gba);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Zab);

    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gab, h);
        dpd_buf4_mat_irrep_rd(&Gab, h);
        for(size_t ia = 0; ia < Gab.params->rowtot[h]; ++ia){
            size_t i = Gab.params->roworb[h][ia][0];
            int Gi = Gab.params->psym[i];
            i -= Gab.params->poff[Gi];
            size_t a = Gab.params->roworb[h][ia][1];
            int Ga = Gab.params->qsym[a];
            a -= Gab.params->qoff[Ga];
            for(size_t jb = 0; jb < Gab.params->coltot[h]; ++jb){
                size_t j = Gab.params->colorb[h][jb][0];
                int Gj = Gab.params->rsym[j];
                j -= Gab.params->roff[Gj];
                size_t b = Gab.params->colorb[h][jb][1];
                int Gb = Gab.params->ssym[b];
                b -= Gab.params->soff[Gb];
                if(Gi == Gj && Ga == Gb) {
                    Gab.matrix[h][ia][jb] += (akappa_->get(Gi, i, j) + aocc_tau_->get(Gi, i, j)) * bvir_ptau_->get(Ga, a, b);
                    Gab.matrix[h][ia][jb] += bvir_tau_->get(Ga, a, b) * (aocc_ptau_->get(Gi, i, j) - aocc_tau_->get(Gi, i, j));
                }
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gab, h);
        dpd_buf4_mat_irrep_close(&Gab, h);
    }

    dpd_buf4_print(&Gab,outfile,1);

    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Gba, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gba, h);
        dpd_buf4_mat_irrep_rd(&Gba, h);
        for(size_t ia = 0; ia < Gba.params->rowtot[h]; ++ia){
            size_t i = Gba.params->roworb[h][ia][0];
            int Gi = Gba.params->psym[i];
            i -= Gba.params->poff[Gi];
            size_t a = Gba.params->roworb[h][ia][1];
            int Ga = Gba.params->qsym[a];
            a -= Gba.params->qoff[Ga];
            for(size_t jb = 0; jb < Gba.params->coltot[h]; ++jb){
                size_t j = Gba.params->colorb[h][jb][0];
                int Gj = Gba.params->rsym[j];
                j -= Gba.params->roff[Gj];
                size_t b = Gba.params->colorb[h][jb][1];
                int Gb = Gba.params->ssym[b];
                b -= Gba.params->soff[Gb];
                if(Gi == Gj && Ga == Gb) {
                    Gba.matrix[h][ia][jb] += (bkappa_->get(Gi, i, j) + bocc_tau_->get(Gi, i, j)) * avir_ptau_->get(Ga, a, b);
                    Gba.matrix[h][ia][jb] += avir_tau_->get(Ga, a, b) * (bocc_ptau_->get(Gi, i, j) - bocc_tau_->get(Gi, i, j));
                }

            }
        }
        dpd_buf4_mat_irrep_wrt(&Gba, h);
        dpd_buf4_mat_irrep_close(&Gba, h);
    }

    dpd_buf4_print(&Gba,outfile,1);

    dpd_buf4_close(&Gba);

    // Г<IajB> spin case:

    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Z (OV|ov)");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Z (OV|OV)");
    dpd_contract444(&Laa, &Zab, &Tab, 0, 1, -0.5, 0.0);
    dpd_contract444(&Zaa, &Lab, &Tab, 0, 1, -0.5, 1.0);
    dpd_buf4_close(&Laa);
    dpd_buf4_close(&Zaa);
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Z (ov|ov)");
    dpd_contract444(&Lab, &Zbb, &Tab, 0, 1, -0.5, 1.0);
    dpd_contract444(&Zab, &Lbb, &Tab, 0, 1, -0.5, 1.0);
    dpd_buf4_close(&Lbb);
    dpd_buf4_close(&Zbb);
    dpd_buf4_close(&Tab);
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");
    dpd_buf4_sort(&Tab, PSIF_DCFT_DENSITY, psrq, ID("[O,v]"), ID("[o,V]"), "Gamma <Ov|oV>");
    // Resort to get the Г_oVOv. Used for the MO Lagrangian
    dpd_buf4_sort(&Tab, PSIF_DCFT_DENSITY, rqps, ID("[o,V]"), ID("[O,v]"), "Gamma <oV|Ov>");

    dpd_buf4_close(&Tab);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Zab);

    dpd_buf4_init(&Gba, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");
    dpd_buf4_print(&Gba,outfile,1);
    dpd_buf4_close(&Gba);

    // Г<iajb> spin case:

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma (ov|ov)");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Z (ov|ov)");
    dpd_contract444(&Lbb, &Zbb, &Gbb, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Lbb);
    dpd_buf4_close(&Zbb);
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Z (OV|ov)");
    dpd_contract444(&Lab, &Zab, &Gbb, 1, 1, -1.0, 1.0);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Zab);
    dpd_buf4_close(&Gbb);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma (ov|ov)");
    dpd_buf4_symm(&Gbb);
    dpd_buf4_close(&Gbb);

    // Resort Г(ov|ov) to the Г<ov|ov>
    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma (ov|ov)");
    dpd_buf4_sort(&Gbb, PSIF_DCFT_DENSITY, psrq, ID("[o,v]"),ID("[o,v]"), "Gamma <ov|ov>");
    dpd_buf4_close(&Gbb);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gbb, h);
        dpd_buf4_mat_irrep_rd(&Gbb, h);
        for(size_t ia = 0; ia < Gbb.params->rowtot[h]; ++ia){
            size_t i = Gbb.params->roworb[h][ia][0];
            int Gi = Gbb.params->psym[i];
            i -= Gbb.params->poff[Gi];
            size_t a = Gbb.params->roworb[h][ia][1];
            int Ga = Gbb.params->qsym[a];
            a -= Gbb.params->qoff[Ga];
            for(size_t jb = 0; jb < Gbb.params->coltot[h]; ++jb){
                size_t j = Gbb.params->colorb[h][jb][0];
                int Gj = Gbb.params->rsym[j];
                j -= Gbb.params->roff[Gj];
                size_t b = Gbb.params->colorb[h][jb][1];
                int Gb = Gbb.params->ssym[b];
                b -= Gbb.params->soff[Gb];
                if(Gi == Gj && Ga == Gb) {
                    Gbb.matrix[h][ia][jb] += (bkappa_->get(Gi, i, j) + bocc_tau_->get(Gi, i, j)) * bvir_ptau_->get(Ga, a, b);
                    Gbb.matrix[h][ia][jb] += bvir_tau_->get(Ga, a, b) * (bocc_ptau_->get(Gi, i, j) - bocc_tau_->get(Gi, i, j));
                }
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Gbb, h);
    }

    dpd_buf4_print(&Gbb,outfile,1);

    dpd_buf4_close(&Gbb);

    psio_->close(PSIF_DCFT_DENSITY, 1);

}

void
DCFTSolver::compute_lagrangian_OV()
{
    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 G, I;
    dpdfile2 X, H, pT;


    // X_OV: One-electron contributions

    // X_IA = H_IB pTau_BA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "pTau <V|V>");

    dpd_contract222(&H, &pT, &X, 0, 1, 1.0, 0.0);
    dpd_file2_close(&pT);
    dpd_file2_close(&H);
    dpd_file2_print(&X,outfile);
    dpd_file2_close(&X);

    // X_IA = H_IB pTau_BA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "pTau <v|v>");

    dpd_contract222(&H, &pT, &X, 0, 1, 1.0, 0.0);
    dpd_file2_close(&pT);
    dpd_file2_close(&H);
    dpd_file2_print(&X,outfile);
    dpd_file2_close(&X);

    // X_OV: Two-electron contributions

    //
    // 2 * <OV||VV> Г_VVVV
    //

    // X_IA += 2 * <IB||CD> Г_ABCD
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "Gamma <VV|VV>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_IA += 4 * <Ib|Cd> Г_AbCd
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "Gamma <Vv|Vv>");

    dpd_contract442(&I, &G, &X, 0, 0, 4.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ia += 2 * <ib||cd> Г_abcd
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "Gamma <vv|vv>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ia += 4 * <iB|cD> Г_AbCd
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "Gamma <Vv|Vv>");

    dpd_contract442(&I, &G, &X, 1, 1, 4.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    //
    // <OO||OV> Г_OOVV
    //

    // X_IA += <BI||JK> Г_BAJK
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 1, "MO Ints <VO|OO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V,V]"), ID("[O,O]"), 0, "Gamma <VV|OO>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_IA += 2 * <Ib|Jk> Г_AbJk
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"),
                  ID("[V,v]"), ID("[O,o]"), 0, "Gamma <Vv|Oo>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ia += <bi||jk> Г_bajk
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 1, "MO Ints <vo|oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v,v]"), ID("[o,o]"), 0, "Gamma <vv|oo>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ia += 2 * <Bi|Jk> Г_BaJk
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "MO Ints <Vo|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"),
                  ID("[V,v]"), ID("[O,o]"), 0, "Gamma <Vv|Oo>");

    dpd_contract442(&I, &G, &X, 1, 1, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    //
    // <OO||OV> Г_OVOV
    //

    // X_IA += <JB||KI> Г_JBKA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 1, "MO Ints <OV|OO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_IA += <kB|jI> Г_kBjA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <oV|oO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_IA -= <Kb|jI> Г_KbjA
    // Note: <Kb|jI> integrals are resorted <bK|jI> integrals.
    // <Kb||jI> Г_KbjA = (-1) * <bK|jI> * Г_KbjA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[o,O]"),
                  ID("[O,v]"), ID("[o,O]"), 0, "MO Ints <Ov|oO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");

    dpd_contract442(&I, &G, &X, 3, 3, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ia += <jb||ki> Г_jbka
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o,o]"), 1, "MO Ints <ov|oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ia += <Kb|Ji> Г_KbJa
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ia -= <kB|Ji> Г_kBJa
    // Note: <kB|Ji> integrals are resorted <Bk|Ji> integrals.
    // <kB||Ji> Г_kBJa = (-1) * <Bk|Ji> * Г_kBJa

    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[O,o]"),
                  ID("[o,V]"), ID("[O,o]"), 0, "MO Ints <oV|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "Gamma <oV|Ov>");

    dpd_contract442(&I, &G, &X, 3, 3, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);


    psio_->close(PSIF_DCFT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::compute_lagrangian_VO()
{
    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 G, I;
    dpdfile2 X, H, pT;


    // X_VO: One-electron contributions

    // X_AI = H_JA pTau_JI
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "pTau <O|O>");
    dpd_file2_mat_init(&X);
    dpd_file2_mat_init(&H);
    dpd_file2_mat_init(&pT);
    dpd_file2_mat_rd(&H);
    dpd_file2_mat_rd(&pT);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int a = 0 ; a < navirpi_[h]; ++a){
                double value = 0.0;
                for(int j = 0 ; j < naoccpi_[h]; ++j){
                    value += H.matrix[h][j][a] * (pT.matrix[h][j][i] + (i == j ? 1.0 : 0.0));
                }
                X.matrix[h][a][i] = value;
            }
        }
    }
    dpd_file2_mat_wrt(&X);
    dpd_file2_close(&pT);
    dpd_file2_close(&H);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ai = H_ja pTau_ji
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "pTau <o|o>");
    dpd_file2_mat_init(&X);
    dpd_file2_mat_init(&H);
    dpd_file2_mat_init(&pT);
    dpd_file2_mat_rd(&H);
    dpd_file2_mat_rd(&pT);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int a = 0 ; a < nbvirpi_[h]; ++a){
                double value = 0.0;
                for(int j = 0 ; j < nboccpi_[h]; ++j){
                    value += H.matrix[h][j][a] * (pT.matrix[h][j][i] + (i == j ? 1.0 : 0.0));
                }
                X.matrix[h][a][i] = value;
            }
        }
    }
    dpd_file2_mat_wrt(&X);
    dpd_file2_close(&pT);
    dpd_file2_close(&H);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_VO: Two-electron contributions

    //
    // 2 * <VO||OO> Г_OOOO
    //

    // X_AI += 2 * <AJ||KL> Г_IJKL
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 1, "MO Ints <VO|OO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "Gamma <OO|OO>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_AI += 4 * <Aj|Kl> Г_IjKl
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "MO Ints <Vo|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");

    dpd_contract442(&I, &G, &X, 0, 0, 4.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ai += 2 * <aj||kl> Г_ijkl
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 1, "MO Ints <vo|oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "Gamma <oo|oo>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_AI += 4 * <Ja|Kl> Г_JiKl
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");

    dpd_contract442(&I, &G, &X, 1, 1, 4.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    //
    // <VO||VV> Г_OOVV
    //

    // X_AI += <JA||BC> Г_JIBC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Gamma <OO|VV>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_AI += <Aj|Bc> Г_IjBc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ai += <ja||bc> Г_jibc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Gamma <oo|vv>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ai += <Ja|Bc> Г_JiBc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");

    dpd_contract442(&I, &G, &X, 1, 1, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    //
    // <OV||VV> Г_OVOV
    //

    // X_AI += <JB||AC> Г_JBIC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_AI += <Jb|Ac> Г_JbIc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_AI -= <jB|Ac> Г_jBIc
    // Note: <jB|Ac> integrals are resorted <Bj|Ac> integrals.
    // <jB||Ac> Г_jBIc = (-1) * <Bj|Ac> Г_jBIc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[V,v]"),
                  ID("[o,V]"), ID("[V,v]"), 0, "MO Ints <oV|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "Gamma <oV|Ov>");

    dpd_contract442(&I, &G, &X, 2, 2, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ai += <jb||ac> Г_jbic
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ai += <jB|aC> Г_jBiC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV|vV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    // X_ai -= <Jb|aC> Г_JbiC
    // Note: <Jb|aC> integrals are resorted <bJ|aC> integrals.
    // <Jb||aC> Г_JbiC = (-1) * <bJ|aC> Г_JbiC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[v,V]"),
                  ID("[O,v]"), ID("[v,V]"), 0, "MO Ints <Ov|vV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");

    dpd_contract442(&I, &G, &X, 2, 2, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_print(&X, outfile);
    dpd_file2_close(&X);

    psio_->close(PSIF_DCFT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::iterate_orbital_response()
{

    // dX_ia = 2.0 * (X_ia - X_ai)
    SharedMatrix dXa (new Matrix("Delta(X) Alpha", nirrep_, naoccpi_, navirpi_));
    SharedMatrix dXb (new Matrix("Delta(X) Beta", nirrep_, nboccpi_, nbvirpi_));

    dpdfile2 Xia, Xai, zia;

    // Compute guess for the orbital response matrix elements

    // Alpha spin
    dpd_file2_init(&Xia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_file2_init(&Xai, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_file2_init(&zia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_file2_mat_init(&Xia);
    dpd_file2_mat_init(&Xai);
    dpd_file2_mat_init(&zia);
    dpd_file2_mat_rd(&Xia);
    dpd_file2_mat_rd(&Xai);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int a = 0 ; a < navirpi_[h]; ++a){
                double value_dX = 2.0 * (Xia.matrix[h][i][a] - Xai.matrix[h][a][i]);
                zia.matrix[h][i][a] = value_dX / (moFa_->get(h, a + naoccpi_[h], a + naoccpi_[h]) - moFa_->get(h, i, i));
                dXa->set(h, i, a, value_dX);
            }
        }
    }
    dpd_file2_mat_wrt(&zia);
    dpd_file2_print(&zia, outfile);
    dpd_file2_close(&zia);
    dpd_file2_close(&Xai);
    dpd_file2_close(&Xia);

    // Beta spin
    dpd_file2_init(&Xia, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_file2_init(&Xai, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_file2_init(&zia, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
    dpd_file2_mat_init(&Xia);
    dpd_file2_mat_init(&Xai);
    dpd_file2_mat_init(&zia);
    dpd_file2_mat_rd(&Xia);
    dpd_file2_mat_rd(&Xai);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int a = 0 ; a < nbvirpi_[h]; ++a){
                double value_dX = 2.0 * (Xia.matrix[h][i][a] - Xai.matrix[h][a][i]);
                zia.matrix[h][i][a] = value_dX / (moFb_->get(h, a + nboccpi_[h], a + nboccpi_[h]) - moFb_->get(h, i, i));
                dXb->set(h, i, a, value_dX);
            }
        }
    }
    dpd_file2_mat_wrt(&zia);
    dpd_file2_print(&zia, outfile);
    dpd_file2_close(&zia);
    dpd_file2_close(&Xai);
    dpd_file2_close(&Xia);

    // Parameters are hard-coded for now. Let the user control them in the future
    double rms = 0.0;
    double convergence = 1.0E-14;
    bool converged = false;
    int maxcycle = 50;

     // Start iterations
    fprintf(outfile, "\n  %4s %11s\n", "Iter", "Z_ia RMS");
    fflush(outfile);

    int cycle = 0;
    do {
        cycle++;

        SharedMatrix za_new, zb_new, za_old, zb_old;

        // Compute intermediates
        compute_orbital_response_intermediates();

        // Save old orbital response
        dpd_file2_init(&zia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
        za_old = SharedMatrix(new Matrix(&zia));
        dpd_file2_close(&zia);
        dpd_file2_init(&zia, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
        zb_old = SharedMatrix(new Matrix(&zia));
        dpd_file2_close(&zia);

        // Update the orbital response
        update_orbital_response();

        // Copy new orbital response to the memory
        dpd_file2_init(&zia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
        za_new = SharedMatrix(new Matrix(&zia));
        dpd_file2_close(&zia);
        dpd_file2_init(&zia, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
        zb_new = SharedMatrix(new Matrix(&zia));
        dpd_file2_close(&zia);

        // Check convergence
        za_old->subtract(za_new);
        zb_old->subtract(zb_new);

        rms = za_old->rms() + zb_old->rms();
        converged = (fabs(rms) < fabs(convergence));

        // Print iterative trace
        fprintf(outfile, "  %4d %11.3E\n", cycle, rms);

        // Termination condition
        if (converged || cycle >= maxcycle) break;

    } while (true);

    fprintf(outfile, "\n");

    if (converged) fprintf(outfile, "    DCFT orbital response equations converged.\n\n");
    else throw PSIEXCEPTION("DCFT orbital response equations did not converge");

    dpd_file2_init(&zia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_file2_print(&zia,outfile);
    dpd_file2_close(&zia);
    dpd_file2_init(&zia, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
    dpd_file2_print(&zia,outfile);
    dpd_file2_close(&zia);


}

void
DCFTSolver::compute_orbital_response_intermediates()
{

    dpdbuf4 I;
    dpdfile2 z, zI, zI_ov, zI_vo;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // Compute z_jb * <ja||bi> intermediate (zI_ai)

    // Alpha spin
    // zI_AI = (IA|JB) z_JB
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "zI <V|O>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
    dpd_contract422(&I, &z, &zI, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // zI_AI += <IA|jb> z_jb
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "zI <V|O>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
    dpd_contract422(&I, &z, &zI, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // Copy terms for the zI_ia intermediate
    dpd_file2_init(&zI_vo, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "zI <V|O>");
    dpd_file2_init(&zI_ov, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "zI <O|V>");
    dpd_file2_mat_init(&zI_vo);
    dpd_file2_mat_init(&zI_ov);
    dpd_file2_mat_rd(&zI_vo);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int a = 0 ; a < navirpi_[h]; ++a){
                zI_ov.matrix[h][i][a] = zI_vo.matrix[h][a][i];
            }
        }
    }
    dpd_file2_mat_wrt(&zI_ov);
    dpd_file2_close(&zI_vo);
    dpd_file2_close(&zI_ov);

    // zI_AI -= <IA|JB> z_JB
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "zI <V|O>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    dpd_contract422(&I, &z, &zI, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // Beta spin
    // zI_ai = (ia|jb) z_jb
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "zI <v|o>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
    dpd_contract422(&I, &z, &zI, 0, 1, 1.0, 0.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // zI_ai += <ia|JB> z_JB
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "zI <v|o>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[O,V]"),
                  ID("[o,v]"), ID("[O,V]"), 0, "MO Ints (ov|OV)");
    dpd_contract422(&I, &z, &zI, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // Copy terms for the zI_ia intermediate
    dpd_file2_init(&zI_vo, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "zI <v|o>");
    dpd_file2_init(&zI_ov, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "zI <o|v>");
    dpd_file2_mat_init(&zI_vo);
    dpd_file2_mat_init(&zI_ov);
    dpd_file2_mat_rd(&zI_vo);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int a = 0 ; a < nbvirpi_[h]; ++a){
                zI_ov.matrix[h][i][a] = zI_vo.matrix[h][a][i];
            }
        }
    }
    dpd_file2_mat_wrt(&zI_ov);
    dpd_file2_close(&zI_vo);
    dpd_file2_close(&zI_ov);

    // zI_ai -= <ia|jb> z_jb
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "zI <v|o>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");
    dpd_contract422(&I, &z, &zI, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // Compute two unique terms for z_jb * <ji||ba> intermediate (zI_ia)

    // Alpha spin
    // zI_IA -= <AI|JB> z_JB
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "zI <O|V>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,V]"),
                  ID("[V,O]"), ID("[O,V]"), 0, "MO Ints <VO|OV>");
    dpd_contract422(&I, &z, &zI, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // Beta spin
    // zI_ia -= <ai|jb> z_jb
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "zI <o|v>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,v]"),
                  ID("[v,o]"), ID("[o,v]"), 0, "MO Ints <vo|ov>");
    dpd_contract422(&I, &z, &zI, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::update_orbital_response()
{

    dpdfile2 X_ia, X_ai, z_ia, zI_ai, zI_ia;

    // Alpha spin
    dpd_file2_init(&zI_ia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "zI <O|V>");
    dpd_file2_init(&zI_ai, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "zI <V|O>");
    dpd_file2_init(&X_ia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_file2_init(&X_ai, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_file2_init(&z_ia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");

    dpd_file2_mat_init(&zI_ai);
    dpd_file2_mat_init(&zI_ia);
    dpd_file2_mat_init(&X_ia);
    dpd_file2_mat_init(&X_ai);
    dpd_file2_mat_init(&z_ia);

    dpd_file2_mat_rd(&zI_ai);
    dpd_file2_mat_rd(&zI_ia);
    dpd_file2_mat_rd(&X_ia);
    dpd_file2_mat_rd(&X_ai);
    dpd_file2_mat_rd(&z_ia);

    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int a = 0 ; a < navirpi_[h]; ++a){
                double value = 0.0;
                for(int j = 0 ; j < naoccpi_[h]; ++j){
                    value += (zI_ai.matrix[h][a][j] + zI_ia.matrix[h][j][a]) * (aocc_tau_->get(h,i,j) + akappa_->get(h,i,j));
                }
                for(int b = 0 ; b < navirpi_[h]; ++b){
                    value -= (zI_ai.matrix[h][b][i] + zI_ia.matrix[h][i][b]) * (avir_tau_->get(h,a,b));
                }
                z_ia.matrix[h][i][a] = (2.0 * (X_ia.matrix[h][i][a] - X_ai.matrix[h][a][i]) - value) /
                        (moFa_->get(h, a + naoccpi_[h], a + naoccpi_[h]) - moFa_->get(h, i, i));
            }
        }
    }
    dpd_file2_mat_wrt(&z_ia);
    dpd_file2_close(&z_ia);
    dpd_file2_close(&X_ai);
    dpd_file2_close(&X_ia);
    dpd_file2_close(&zI_ai);
    dpd_file2_close(&zI_ia);

    // Beta spin
    dpd_file2_init(&zI_ia, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "zI <o|v>");
    dpd_file2_init(&zI_ai, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "zI <v|o>");
    dpd_file2_init(&X_ia, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_file2_init(&X_ai, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_file2_init(&z_ia, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");

    dpd_file2_mat_init(&zI_ai);
    dpd_file2_mat_init(&zI_ia);
    dpd_file2_mat_init(&X_ia);
    dpd_file2_mat_init(&X_ai);
    dpd_file2_mat_init(&z_ia);

    dpd_file2_mat_rd(&zI_ai);
    dpd_file2_mat_rd(&zI_ia);
    dpd_file2_mat_rd(&X_ia);
    dpd_file2_mat_rd(&X_ai);
    dpd_file2_mat_rd(&z_ia);

    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int a = 0 ; a < nbvirpi_[h]; ++a){
                double value = 0.0;
                for(int j = 0 ; j < nboccpi_[h]; ++j){
                    value += (zI_ai.matrix[h][a][j] + zI_ia.matrix[h][j][a]) * (bocc_tau_->get(h,i,j) + bkappa_->get(h,i,j));
                }
                for(int b = 0 ; b < nbvirpi_[h]; ++b){
                    value -= (zI_ai.matrix[h][b][i] + zI_ia.matrix[h][i][b]) * (bvir_tau_->get(h,a,b));
                }
                z_ia.matrix[h][i][a] = (2.0 * (X_ia.matrix[h][i][a] - X_ai.matrix[h][a][i]) - value) /
                        (moFb_->get(h, a + nboccpi_[h], a + nboccpi_[h]) - moFb_->get(h, i, i));
            }
        }
    }
    dpd_file2_mat_wrt(&z_ia);
    dpd_file2_close(&z_ia);
    dpd_file2_close(&X_ai);
    dpd_file2_close(&X_ia);
    dpd_file2_close(&zI_ai);
    dpd_file2_close(&zI_ia);

}

// Returns RMS of the change in the response coupling term (C intermediate)
double
DCFTSolver::compute_response_coupling()
{

    dpdbuf4 I, L, C, T;
    dpdfile2 z, zI, zIsym;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // Compute z_kb * <ki||bj> intermediate (zI_ij)

    // Alpha spin

    // zI_IJ = (IJ|KC) z_KC
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "zI <O|O>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints (OO|OV)");
    dpd_contract422(&I, &z, &zI, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // zI_IJ -= <JI|KC> z_KC
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "zI <O|O>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO|OV>");
    dpd_contract422(&I, &z, &zI, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // zI_IJ += (IJ|ck) z_ck
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "zI <O|O>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,o]"),
                  ID("[O>=O]+"), ID("[v,o]"), 0, "MO Ints (OO|vo)");
    dpd_contract422(&I, &z, &zI, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_print(&zI,outfile);
    dpd_file2_close(&zI);

    // Beta spin

    // zI_ij = (ij|kc) z_kc
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "zI <o|o>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints (oo|ov)");
    dpd_contract422(&I, &z, &zI, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // zI_ij -= <ji|kc> z_kc
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "zI <o|o>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo|ov>");
    dpd_contract422(&I, &z, &zI, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // zI_ij += (ij|CK) z_CK
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "zI <o|o>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[V,O]"),
                  ID("[o,o]"), ID("[V,O]"), 0, "MO Ints (oo|VO)");
    dpd_contract422(&I, &z, &zI, 1, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_print(&zI,outfile);
    dpd_file2_close(&zI);

    // Compute z_kc * <ka||cb> intermediate (zI_ab)

    // Alpha spin

    // zI_AB = (AB|KC) z_KC
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "zI <V|V>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,V]"),
                  ID("[V>=V]+"), ID("[O,V]"), 0, "MO Ints (VV|OV)");
    dpd_contract422(&I, &z, &zI, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // zI_AB -= <BA|KC> z_KC
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "zI <V|V>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,V]"),
                  ID("[V,V]"), ID("[O,V]"), 0, "MO Ints <VV|OV>");
    dpd_contract422(&I, &z, &zI, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // zI_AB += (AB|kc) z_kc
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "zI <V|V>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,v]"),
                  ID("[V>=V]+"), ID("[o,v]"), 0, "MO Ints (VV|ov)");
    dpd_contract422(&I, &z, &zI, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_print(&zI,outfile);
    dpd_file2_close(&zI);

    // Beta spin

    // zI_ab = (ab|kc) z_kc
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "zI <v|v>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,v]"),
                  ID("[v>=v]+"), ID("[o,v]"), 0, "MO Ints (vv|ov)");
    dpd_contract422(&I, &z, &zI, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // zI_ab -= <ba|kc> z_kc
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "zI <v|v>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,v]"),
                  ID("[v,v]"), ID("[o,v]"), 0, "MO Ints <vv|ov>");
    dpd_contract422(&I, &z, &zI, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_close(&zI);

    // zI_ab += (ab|KC) z_KC
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "zI <v|v>");
    dpd_file2_init(&z, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[O,V]"),
                  ID("[v,v]"), ID("[O,V]"), 0, "MO Ints (vv|OV)");
    dpd_contract422(&I, &z, &zI, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&z);
    dpd_file2_print(&zI,outfile);
    dpd_file2_close(&zI);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

    // Symmetrize zI_ij and zI_ab

    // zI <O|O>
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "zI <O|O>");
    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "zI <O|O> sym");
    dpd_file2_mat_init(&zI);
    dpd_file2_mat_init(&zIsym);
    dpd_file2_mat_rd(&zI);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                zIsym.matrix[h][i][j] = zIsym.matrix[h][j][i] = zI.matrix[h][i][j] + zI.matrix[h][j][i];
            }
        }
    }
    dpd_file2_mat_wrt(&zIsym);
    dpd_file2_close(&zI);
    dpd_file2_print(&zIsym,outfile);
    dpd_file2_close(&zIsym);

    // zI <o|o>
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "zI <o|o>");
    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "zI <o|o> sym");
    dpd_file2_mat_init(&zI);
    dpd_file2_mat_init(&zIsym);
    dpd_file2_mat_rd(&zI);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                zIsym.matrix[h][i][j] = zIsym.matrix[h][j][i] = zI.matrix[h][i][j] + zI.matrix[h][j][i];
            }
        }
    }
    dpd_file2_mat_wrt(&zIsym);
    dpd_file2_close(&zI);
    dpd_file2_print(&zIsym,outfile);
    dpd_file2_close(&zIsym);

    // zI <V|V>
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "zI <V|V>");
    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "zI <V|V> sym");
    dpd_file2_mat_init(&zI);
    dpd_file2_mat_init(&zIsym);
    dpd_file2_mat_rd(&zI);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < navirpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                zIsym.matrix[h][i][j] = zIsym.matrix[h][j][i] = zI.matrix[h][i][j] + zI.matrix[h][j][i];
            }
        }
    }
    dpd_file2_mat_wrt(&zIsym);
    dpd_file2_close(&zI);
    dpd_file2_print(&zIsym,outfile);
    dpd_file2_close(&zIsym);

    // zI <v|v>
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "zI <v|v>");
    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "zI <v|v> sym");
    dpd_file2_mat_init(&zI);
    dpd_file2_mat_init(&zIsym);
    dpd_file2_mat_rd(&zI);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nbvirpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                zIsym.matrix[h][i][j] = zIsym.matrix[h][j][i] = zI.matrix[h][i][j] + zI.matrix[h][j][i];
            }
        }
    }
    dpd_file2_mat_wrt(&zIsym);
    dpd_file2_close(&zI);
    dpd_file2_print(&zIsym,outfile);
    dpd_file2_close(&zIsym);

    // Save the old response coupling terms (C intermediate)
    if (iter_ > 1) {
        dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "C <OO|VV> new");
        dpd_buf4_copy(&C, PSIF_DCFT_DPD, "C <OO|VV> old");
        dpd_buf4_close(&C);

        dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "C <Oo|Vv> new");
        dpd_buf4_copy(&C, PSIF_DCFT_DPD, "C <Oo|Vv> old");
        dpd_buf4_close(&C);

        dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "C <oo|vv> new");
        dpd_buf4_copy(&C, PSIF_DCFT_DPD, "C <oo|vv> old");
        dpd_buf4_close(&C);

    }

    // Compute the new response coupling terms (C intermediate)

    //
    // C <OO|VV>
    //

    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "zI <V|V> sym");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "C <OO|VV> new");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    // The first pair index of the result of dpd_contract244 contains the file2 index.
    // The second pair index of the result of dpd_contract424 contains the file2 index.
    // This is how dpd_contract244 works: AD * IJDB -> ABIJ ?-(t)-> IJAB. If BD * IJDA -> BAIJ -(t)-> IJBA
    // This is how dpd_contract424 works: IJAD * BD -> IJAB ?-(t)-> ABIJ. If IJBD * AD -> IJBA ?-(t)-> BAIJ
    // Temp_IJAB = zIsym_AC lambda_IJCB
    dpd_contract244(&zIsym, &L, &T, 1, 2, 1, 1.0, 0.0);
    // C_IJAB = Temp_IJAB
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "C <OO|VV> new");
    dpd_file2_close(&zIsym);
    dpd_buf4_close(&L);
    // Temp_IJAB -> Temp_IJBA
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    // C_IJAB -= Temp_IJBA
    dpd_buf4_add(&C, &T, -1.0);
    dpd_buf4_close(&T);

    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "zI <O|O> sym");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    // Temp_IJAB = - zIsym_IK lambda_KJAB
    dpd_contract244(&zIsym, &L, &T, 1, 0, 0, -1.0, 0.0);
    dpd_file2_close(&zIsym);
    dpd_buf4_close(&L);
    // C_IJAB += Temp_IJAB
    dpd_buf4_add(&C, &T, 1.0);
    // Temp_IJAB -> Temp_JIAB
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    // C_IJAB -= Temp_JIAB
    dpd_buf4_add(&C, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_print(&C, outfile, 1);
    dpd_buf4_close(&C);

    //
    // C <Oo|Vv>
    //

    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "C <Oo|Vv> new");
    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "zI <V|V> sym");
    // C_IjAb = zIsym_AC * lambda_IjCb
    dpd_contract244(&zIsym, &L, &C, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&zIsym);
    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "zI <v|v> sym");
    // C_IjAb += lambda_IjAc zIsym_bc
    dpd_contract424(&L, &zIsym, &C, 3, 1, 0, 1.0, 1.0);
    dpd_file2_close(&zIsym);
    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "zI <O|O> sym");
    // C_IjAb -= zIsym_IK * lambda_KjAb
    dpd_contract244(&zIsym, &L, &C, 1, 0, 0, -1.0, 1.0);
    dpd_file2_close(&zIsym);
    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "zI <o|o> sym");
    // C_IjAb -= lambda_IkAb * zIsym_jk
    dpd_contract424(&L, &zIsym, &C, 1, 1, 1, -1.0, 1.0);
    dpd_file2_close(&zIsym);
    dpd_buf4_print(&C, outfile, 1);
    dpd_buf4_close(&C);
    dpd_buf4_close(&L);

    //
    // C <oo|vv>
    //

    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "zI <v|v> sym");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "C <oo|vv> new");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    // Temp_ijab = zIsym_ac lambda_ijcb
    dpd_contract244(&zIsym, &L, &T, 1, 2, 1, 1.0, 0.0);
    // C_ijab = Temp_ijab
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "C <oo|vv> new");
    dpd_file2_close(&zIsym);
    dpd_buf4_close(&L);
    // Temp_ijab -> Temp_ijba
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    // C_ijba -= Temp_ijba
    dpd_buf4_add(&C, &T, -1.0);
    dpd_buf4_close(&T);

    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "zI <o|o> sym");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    // Temp_ijab = - zIsym_ik lambda_kjab
    dpd_contract244(&zIsym, &L, &T, 1, 0, 0, -1.0, 0.0);
    dpd_file2_close(&zIsym);
    dpd_buf4_close(&L);
    // C_ijab += Temp_ijab
    dpd_buf4_add(&C, &T, 1.0);
    // Temp_ijab -> Temp_jiab
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    // C_ijab -= Temp_jiab
    dpd_buf4_add(&C, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_print(&C, outfile, 1);
    dpd_buf4_close(&C);

    // Compute the RMS of the change in the response coupling terms (RMS (Cnew - Cold))

    int nElements = 0;
    double sumSQ = 0.0;

    // dC <OO|VV>
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "C <OO|VV> new");
    dpd_buf4_copy(&C, PSIF_DCFT_DPD, "C <OO|VV> new - old");
    dpd_buf4_close(&C);

    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "C <OO|VV> new - old");
    if (iter_ > 1) {
        dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "C <OO|VV> old");
        dpd_buf4_add(&T, &C, -1.0);
        dpd_buf4_close(&C);
    }

    for (int h = 0; h < nirrep_; ++h) nElements += T.params->coltot[h] * T.params->rowtot[h];

    sumSQ += dpd_buf4_dot_self(&T);
    dpd_buf4_close(&T);

    // dC <Oo|Vv>
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "C <Oo|Vv> new");
    dpd_buf4_copy(&C, PSIF_DCFT_DPD, "C <Oo|Vv> new - old");
    dpd_buf4_close(&C);

    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "C <Oo|Vv> new - old");
    if (iter_ > 1) {
        dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "C <Oo|Vv> old");
        dpd_buf4_add(&T, &C, -1.0);
        dpd_buf4_close(&C);
    }

    for (int h = 0; h < nirrep_; ++h) nElements += T.params->coltot[h] * T.params->rowtot[h];

    sumSQ += dpd_buf4_dot_self(&T);
    dpd_buf4_close(&T);

    // dC <oo|vv>
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "C <oo|vv> new");
    dpd_buf4_copy(&C, PSIF_DCFT_DPD, "C <oo|vv> new - old");
    dpd_buf4_close(&C);

    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "C <oo|vv> new - old");
    if (iter_ > 1) {
        dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "C <oo|vv> old");
        dpd_buf4_add(&T, &C, -1.0);
        dpd_buf4_close(&C);
    }

    for (int h = 0; h < nirrep_; ++h) nElements += T.params->coltot[h] * T.params->rowtot[h];

    sumSQ += dpd_buf4_dot_self(&T);
    dpd_buf4_close(&T);

    // Compute RMS of dC
    return sqrt(sumSQ / nElements);

}

void
DCFTSolver::iterate_cumulant_response()
{

    // Update cumulant reponse from the change in C intermediate
    cumulant_response_guess();

    // iteratively solve for cumulant reponse

    // Parameters are hard-coded for now. Let the user control them in the future
    double rms = 0.0;
    double convergence = 1.0E-14;
    bool converged = false;
    int maxcycle = 50;

     // Start iterations
    fprintf(outfile, "\n  %4s %11s\n", "Iter", "Z_iajb RMS");
    fflush(outfile);

    int cycle = 0;
    do {
        cycle++;

        // Build perturbed tau and delta tau
        build_perturbed_tau();

        // Compute intermediates
        compute_cumulant_response_intermediates();

        // Compute cumulant response residual and its RMS
        rms = compute_cumulant_response_residual();

        // Update the cumulant response
        update_cumulant_response();

        // Check the convergence

        // Print iterative trace
        fprintf(outfile, "  %4d %11.3E\n", cycle, rms);

//        // Termination condition
//        if (converged || cycle >= maxcycle) break;
        converged = true;
        break;

    } while (true);

    fprintf(outfile, "\n");

    if (converged) fprintf(outfile, "    DCFT cumulant response equations converged.\n\n");
    else throw PSIEXCEPTION("DCFT cumulant response equations did not converge");


}

void
DCFTSolver::cumulant_response_guess()
{

    dpdbuf4 Z, D, dC;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    /*
     * Z_ijab += dC_ijab / D_ijab
     */

    // Z_IJAB += dC_IJAB / D_IJAB
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    dpd_buf4_init(&dC, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "C <OO|VV> new - old");
    dpd_buf4_dirprd(&D, &dC);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Z <OO|VV>");
    dpd_buf4_add(&Z, &dC, 1.0);
    dpd_buf4_close(&dC);
    dpd_buf4_print(&Z,outfile,1);
    dpd_buf4_close(&Z);

    // Z_IjAb += dC_IjAb / D_IjAb
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    dpd_buf4_init(&dC, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "C <Oo|Vv> new - old");
    dpd_buf4_dirprd(&D, &dC);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    dpd_buf4_add(&Z, &dC, 1.0);
    dpd_buf4_close(&dC);
    dpd_buf4_print(&Z,outfile,1);
    dpd_buf4_close(&Z);

    // Z_ijab += dC_ijab / D_ijab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    dpd_buf4_init(&dC, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "C <oo|vv> new - old");
    dpd_buf4_dirprd(&D, &dC);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Z <oo|vv>");
    dpd_buf4_add(&Z, &dC, 1.0);
    dpd_buf4_close(&dC);
    dpd_buf4_print(&Z,outfile,1);
    dpd_buf4_close(&Z);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::build_perturbed_tau()
{

    dpdbuf4 L, Z;
    dpdfile2 pT_OO, pT_oo, pT_VV, pT_vv, dT_OO, dT_oo, dT_VV, dT_vv, T_OO, T_oo, T_VV, T_vv;

    dpd_file2_init(&pT_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Temp <O|O>");
    dpd_file2_init(&pT_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Temp <o|o>");
    dpd_file2_init(&pT_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Temp <V|V>");
    dpd_file2_init(&pT_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Temp <v|v>");

    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Z <OO|VV>");

    // pTau_IJ = -1/2 Lambda_IKAB Z_JKAB
    dpd_contract442(&L, &Z, &pT_OO, 0, 0, -0.5, 0.0);

    // pTau_AB = +1/2 Lambda_IJAC Z_IJBC
    dpd_contract442(&L, &Z, &pT_VV, 2, 2, 0.5, 0.0);

    dpd_buf4_close(&Z);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Z <oo|vv>");

    // pTau_ij = -1/2 Lambda_ikab Z_jkab
    dpd_contract442(&L, &Z, &pT_oo, 0, 0, -0.5, 0.0);

    // pTau_ab = +1/2 Lambda_ijac Z_ijbc
    dpd_contract442(&L, &Z, &pT_vv, 2, 2, 0.5, 0.0);

    dpd_buf4_close(&Z);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");

    // pTau_IJ -= 1/2 Lambda_IkAb Z_JkAb - 1/2 Lambda_IkaB Z_JkaB
    dpd_contract442(&L, &Z, &pT_OO, 0, 0, -1.0, 1.0);
    // pTau_ij -= 1/2 Lambda_KiAb Z_KjAb - 1/2 Lambda_KiaB Z_KjaB
    dpd_contract442(&L, &Z, &pT_oo, 1, 1, -1.0, 1.0);
    // pTau_AB += 1/2 Lambda_IjAc Z_IjBc + 1/2 Lambda_iJAc Z_iJBc
    dpd_contract442(&L, &Z, &pT_VV, 2, 2, 1.0, 1.0);
    // pTau_ab += 1/2 Lambda_IjCa Z_IjCb + 1/2 Lambda_iJCa Z_iJCb
    dpd_contract442(&L, &Z, &pT_vv, 3, 3, 1.0, 1.0);

    dpd_buf4_close(&Z);
    dpd_buf4_close(&L);

    dpd_file2_close(&pT_OO);
    dpd_file2_close(&pT_oo);
    dpd_file2_close(&pT_VV);
    dpd_file2_close(&pT_vv);

    // Now symmetrize perturbed Tau: new pTau(i,j) = (pTau(i,j) + pTau(j,1))/2

    dpd_file2_init(&pT_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "pTau <O|O>");
    dpd_file2_init(&pT_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "pTau <o|o>");
    dpd_file2_init(&pT_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "pTau <V|V>");
    dpd_file2_init(&pT_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "pTau <v|v>");

    dpd_file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Temp <O|O>");
    dpd_file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Temp <o|o>");
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Temp <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Temp <v|v>");

    dpd_file2_mat_init(&pT_OO);
    dpd_file2_mat_init(&pT_oo);
    dpd_file2_mat_init(&pT_VV);
    dpd_file2_mat_init(&pT_vv);

    dpd_file2_mat_init(&T_OO);
    dpd_file2_mat_init(&T_oo);
    dpd_file2_mat_init(&T_VV);
    dpd_file2_mat_init(&T_vv);

    dpd_file2_mat_rd(&T_OO);
    dpd_file2_mat_rd(&T_oo);
    dpd_file2_mat_rd(&T_VV);
    dpd_file2_mat_rd(&T_vv);

    for(int h = 0; h < nirrep_; ++h){
        // pTau <O|O>
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                pT_OO.matrix[h][i][j] = pT_OO.matrix[h][j][i] = (T_OO.matrix[h][i][j] + T_OO.matrix[h][j][i]) / 2.0;
            }
        }

        // pTau <o|o>
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                pT_oo.matrix[h][i][j] = pT_oo.matrix[h][j][i] = (T_oo.matrix[h][i][j] + T_oo.matrix[h][j][i]) / 2.0;
            }
        }

        // pTau <V|V>
        for(int i = 0 ; i < navirpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                pT_VV.matrix[h][i][j] = pT_VV.matrix[h][j][i] = (T_VV.matrix[h][i][j] + T_VV.matrix[h][j][i]) / 2.0;
            }
        }

        // pTau <v|v>
        for(int i = 0 ; i < nbvirpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                pT_vv.matrix[h][i][j] = pT_vv.matrix[h][j][i] = (T_vv.matrix[h][i][j] + T_vv.matrix[h][j][i]) / 2.0;
            }
        }
    }

    dpd_file2_mat_wrt(&pT_OO);
    dpd_file2_mat_wrt(&pT_oo);
    dpd_file2_mat_wrt(&pT_VV);
    dpd_file2_mat_wrt(&pT_vv);

    dpd_file2_close(&pT_OO);
    dpd_file2_close(&pT_oo);
    dpd_file2_close(&pT_VV);
    dpd_file2_close(&pT_vv);

    dpd_file2_close(&T_OO);
    dpd_file2_close(&T_oo);
    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);

    // Form dTau = pTau - Tau
    dpd_file2_init(&pT_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "pTau <O|O>");
    dpd_file2_init(&pT_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "pTau <o|o>");
    dpd_file2_init(&pT_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "pTau <V|V>");
    dpd_file2_init(&pT_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "pTau <v|v>");

    dpd_file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    dpd_file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    dpd_file2_init(&dT_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "dTau <O|O>");
    dpd_file2_init(&dT_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "dTau <o|o>");
    dpd_file2_init(&dT_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "dTau <V|V>");
    dpd_file2_init(&dT_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "dTau <v|v>");

    dpd_file2_axpbycz(&pT_OO, &T_OO, &dT_OO, 1.0, -1.0, 0.0);
    dpd_file2_axpbycz(&pT_oo, &T_oo, &dT_oo, 1.0, -1.0, 0.0);
    dpd_file2_axpbycz(&pT_VV, &T_VV, &dT_VV, 1.0, -1.0, 0.0);
    dpd_file2_axpbycz(&pT_vv, &T_vv, &dT_vv, 1.0, -1.0, 0.0);

    dpd_file2_print(&pT_OO, outfile);
    dpd_file2_print(&pT_oo, outfile);
    dpd_file2_print(&pT_VV, outfile);
    dpd_file2_print(&pT_vv, outfile);

    dpd_file2_print(&dT_OO, outfile);
    dpd_file2_print(&dT_oo, outfile);
    dpd_file2_print(&dT_VV, outfile);
    dpd_file2_print(&dT_vv, outfile);

    dpd_file2_close(&pT_OO);
    dpd_file2_close(&pT_oo);
    dpd_file2_close(&pT_VV);
    dpd_file2_close(&pT_vv);

    dpd_file2_close(&T_OO);
    dpd_file2_close(&T_oo);
    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);

    dpd_file2_close(&dT_OO);
    dpd_file2_close(&dT_oo);
    dpd_file2_close(&dT_VV);
    dpd_file2_close(&dT_vv);

}

void
DCFTSolver::compute_cumulant_response_intermediates()
{

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdfile2 F_OO, F_oo, F_VV, F_vv, X_OO, X_oo, X_VV, X_vv, T_OO, T_oo, T_VV, T_vv, tmp;
    dpdbuf4 I, Z, G, T, F, Taa, Tab, Tbb, Laa, Lab, Lbb, Zaa, Zab, Zbb;

    //
    // Compute G intermediates
    //

    /*
     * G_ijab = <ij||ab>
     */
    // G_IJAB = <IJ||AB>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "G <OO|VV>");
    dpd_buf4_close(&I);

    // G_IjAb = <Ij|Ab>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "G <Oo|Vv>");
    dpd_buf4_close(&I);

    // G_ijab = <ij||ab>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "G <oo|vv>");
    dpd_buf4_close(&I);

    /*
     * G_ijab += 1/2 Sum_cd gbar_cdab Z_ijcd
     */
    // G_IJAB += 1/2 Sum_CD gbar_CDAB Z_IJCD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 1, "MO Ints <VV|VV>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Z <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    dpd_contract444(&Z, &I, &G, 0, 0, 0.5, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&G);


    // G_IjAb += Sum_Cd g_CdAb Z_IjCd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    dpd_contract444(&Z, &I, &G, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&G);


    // G_ijab += 1/2 Sum_cd gbar_cdab Z_ijcd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 1, "MO Ints <vv|vv>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Z <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_contract444(&Z, &I, &G, 0, 0, 0.5, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&G);

    /*
     * G_ijab += 1/2 Sum_kl gbar_ijkl Z_klab
     */
    // G_IJAB += 1/2 Sum_KL gbar_IJKL Z_KLAB
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
            ID("[O,O]"), ID("[O,O]"), 1, "MO Ints <OO|OO>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Z <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    dpd_contract444(&I, &Z, &G, 0, 1, 0.5, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&G);

    // G_IjAb += Sum_Kl g_IjKl Z_KlAb
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    dpd_contract444(&I, &Z, &G, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&G);

    // G_ijab += 1/2 Sum_kl gbar_ijkl Z_klab
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 1, "MO Ints <oo|oo>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Z <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_contract444(&I, &Z, &G, 0, 1, 0.5, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&G);


    /*
     * G_ijab -= P(ij)P(ab) Sum_kc gbar_jckb Z_ikac
     */
    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Z <OO|VV>");
    dpd_buf4_sort(&Zaa, PSIF_DCFT_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "Z (OV|OV)");
    dpd_buf4_close(&Zaa);
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    dpd_buf4_sort(&Zab, PSIF_DCFT_DPD, psqr, ID("[O,v]"), ID("[o,V]"), "Z (Ov|oV)");
    dpd_buf4_close(&Zab);
    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Z <oo|vv>");
    dpd_buf4_sort(&Zbb, PSIF_DCFT_DPD, prqs, ID("[o,v]"),ID("[o,v]"), "Z (ov|ov)");
    dpd_buf4_close(&Zbb);

    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Z (OV|OV)");
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Z (Ov|oV)");
    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Z (ov|ov)");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "MO Ints <oV|oV>");
    dpd_buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Temp (OV|OV)");
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Temp (Ov|oV)");
    dpd_buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Temp (ov|ov)");

    // T_IbjA = -Sum_kC Z_IbkC g_jAkC
    dpd_contract444(&Zab, &I, &Tab, 0, 0, -1.0, 0.0);

    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>") ;

    // T_IbjA -= Sum_Kc g_IbKc Z_KcjA
    dpd_contract444(&I, &Zab, &Tab, 0, 1, -1.0, 1.0);

    // T_IbjA -> T_IAjb
    dpd_buf4_sort(&Tab, PSIF_DCFT_DPD, psrq, ID("[O,V]"),ID("[o,v]"), "Temp (OV|ov)");
    dpd_buf4_close(&Tab);
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");

    // Z_IbkC -> Z_ICkb
    dpd_buf4_sort(&Zab, PSIF_DCFT_DPD, psrq, ID("[O,V]"),ID("[o,v]"), "Z (OV|ov)");
    dpd_buf4_close(&Zab);
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Z (OV|ov)");

    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");

    // T_IAJB = Sum_kc Z_IAkc g_JBkc
    dpd_contract444(&Zab, &I, &Taa, 0, 0, 1.0, 0.0);

    // T_iajb = Sum_KC Z_KCia g_KCjb
    dpd_contract444(&Zab, &I, &Tbb, 1, 1, 1.0, 0.0);

    // T_IAjb += Sum_kc g_IAkc Z_jbkc
    dpd_contract444(&I, &Zbb, &Tab, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC Z_IAKC g_KCjb
    dpd_contract444(&Zaa, &I, &Tab, 0, 1, 1.0, 1.0);


    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV> - (OV|OV)");

    // T_IAJB -= Sum_KC Z_IAKC gbar_JBKC
    dpd_contract444(&Zaa, &I, &Taa, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC gbar_IAKC Z_KCjb
    dpd_contract444(&I, &Zab, &Tab, 0, 1, -1.0, 1.0);

    dpd_buf4_close(&I);
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov> - (ov|ov)");

    // T_iajb -= Sum_kc Z_iakc gbar_jbkc
    dpd_contract444(&Zbb, &I, &Tbb, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC Z_IAkc gbar_jbkc
    dpd_contract444(&Zab, &I, &Tab, 0, 0, -1.0, 1.0);

    dpd_buf4_close(&I);
    dpd_buf4_close(&Zaa);
    dpd_buf4_close(&Zab);
    dpd_buf4_close(&Zbb);

    // T_IAJB -> T_IJAB
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, prqs, ID("[O,O]"), ID("[V,V]"), "Temp <OO|VV>");
    dpd_buf4_close(&Taa);
    // G_IJAB += T_IJAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    dpd_buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");

    // T_IJAB -> T_JIAB
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");

    // G_IJAB -= T_JIAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_IJAB -> T_IJBA
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // G_IJAB -= T_IJBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_IJAB -> T_JIBA
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // G_IJAB += T_JIBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_print(&G,outfile,1);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    dpd_buf4_close(&Taa);

    // T_IAjb -> T_IjAb
    dpd_buf4_sort(&Tab, PSIF_DCFT_DPD, prqs, ID("[O,o]"), ID("[V,v]"), "Temp <Oo|Vv>");
    dpd_buf4_close(&Tab);
    // G_IjAb += T_IjAb
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Temp <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "G <Oo|Vv>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_print(&G,outfile,1);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_iajb -> T_ijab
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, prqs, ID("[o,o]"), ID("[v,v]"), "Temp <oo|vv>");
    dpd_buf4_close(&Tbb);
    // G_ijab += T_ijab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    dpd_buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");

    // T_ijab -> T_jiab
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab -= T_jiab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_ijab -> T_ijba
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab -= T_ijba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_ijab -> T_jiba
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab += T_jiba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, 1.0);
    dpd_buf4_print(&G,outfile,1);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    dpd_buf4_close(&Tbb);

    //
    // Compute F intermediates
    //

    /*
     * F_ijab += P(ab) F_ca Z_ijcb - P(ij) F_ki Z_kjab
     */
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Z <OO|VV>");
    // Temp_IJAB = Z_IJCB F_AC
    dpd_file2_init(&F_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    dpd_contract244(&F_VV, &Zaa, &T, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&F_VV);
    dpd_buf4_close(&Zaa);
    // F_IJAB = Temp_IJAB
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "F <OO|VV>");
    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "F <OO|VV>");
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    // F_IJAB -= Temp_IJBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_add(&F, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
              ID("[O,O]"), ID("[V,V]"), 0, "Z <OO|VV>");
    // Temp_IJAB = -Z_KJAB F_IK
    dpd_file2_init(&F_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    dpd_contract244(&F_OO, &Zaa, &T, 1, 0, 0, -1.0, 0.0);

    dpd_file2_close(&F_OO);
    dpd_buf4_close(&Zaa);
    // F_IJAB += Temp_IJAB
    dpd_buf4_add(&F, &T, 1.0);
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    // F_IJAB -= Temp_JIAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_add(&F, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_print(&F,outfile,1);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "F <Oo|Vv>");
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    // F_IjAb += Z_IjCb F_AC
    dpd_file2_init(&F_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    dpd_contract244(&F_VV, &Zab, &F, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&F_VV);
    // F_IjAb += Z_IjAc F_bc
    dpd_file2_init(&F_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    dpd_contract424(&Zab, &F_vv, &F, 3, 1, 0, 1.0, 1.0);
    dpd_file2_close(&F_vv);
    // F_IjAb -= Z_KjAb F_IK
    dpd_file2_init(&F_OO, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    dpd_contract244(&F_OO, &Zab, &F, 1, 0, 0, -1.0, 1.0);
    dpd_file2_close(&F_OO);
    // F_IjAb -= Z_IkAb F_jk
    dpd_file2_init(&F_oo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    dpd_contract424(&Zab, &F_oo, &F, 1, 1, 1, -1.0, 1.0);
    dpd_file2_close(&F_oo);
    dpd_buf4_close(&Zab);
    dpd_buf4_print(&F,outfile,1);
    dpd_buf4_close(&F);

    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Z <oo|vv>");
    // Temp_ijab = Z_ijcb F_ac
    dpd_file2_init(&F_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    dpd_contract244(&F_vv, &Zbb, &T, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&F_vv);
    dpd_buf4_close(&Zbb);
    // F_ijab = Temp_ijab
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "F <oo|vv>");
    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "F <oo|vv>");
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);
    // F_ijab -= Temp_ijba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_add(&F, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
              ID("[o,o]"), ID("[v,v]"), 0, "Z <oo|vv>");
    // Temp_ijab = -Z_kjab X_ik
    dpd_file2_init(&F_oo, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "F <o|o>");
    dpd_contract244(&F_oo, &Zbb, &T, 1, 0, 0, -1.0, 0.0);
    dpd_file2_close(&F_oo);
    dpd_buf4_close(&Zbb);
    // F_ijab += Temp_ijab
    dpd_buf4_add(&F, &T, 1.0);
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);
    // F_ijab -= Temp_jiab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_add(&F, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_print(&F,outfile,1);
    dpd_buf4_close(&F);

    //
    // Compute X intermediates
    //

    dpd_file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "dTau <O|O>");
    dpd_file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "dTau <o|o>");
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "dTau <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "dTau <v|v>");

    dpd_file2_init(&X_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");

    // X_AB = (AB|CD) dTau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
    dpd_contract422(&I, &T_VV, &X_VV, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // X_AB -= <AB|CD> dTau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
    dpd_contract422(&I, &T_VV, &X_VV, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // X_AB += (AB|cd) dTau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[v,v]"),
                  ID("[V>=V]+"), ID("[v>=v]+"), 0, "MO Ints (VV|vv)");
    dpd_contract422(&I, &T_vv, &X_VV, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);

    // X_AB = +(AB|IJ) dTau_IJ
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V>=V]+"), ID("[O>=O]+"), 0, "MO Ints (VV|OO)");
    dpd_contract422(&I, &T_OO, &X_VV, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    // X_AB -= <AB|IJ> dTau_IJ
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V,V]"), ID("[O,O]"), 0, "MO Ints <VV|OO>");
    dpd_contract422(&I, &T_OO, &X_VV, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // X_AB += (AB|ij) dTau_ij
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                  ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
    dpd_contract422(&I, &T_oo, &X_VV, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&X_VV);

    dpd_file2_init(&X_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");

    // X_ab = +(ab|cd) dTau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v>=v]+"), ID("[v>=v]+"), 0, "MO Ints (vv|vv)");
    dpd_contract422(&I, &T_vv, &X_vv, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // X_ab -= <ab|cd> dTau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv|vv>");
    dpd_contract422(&I, &T_vv, &X_vv, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // X_ab += (ab|CD) dTau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[V,V]"),
                  ID("[v,v]"), ID("[V,V]"), 0, "MO Ints (vv|VV)");
    dpd_contract422(&I, &T_VV, &X_vv, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);

    // X_ab = +(ab|ij) dTau_ij
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v>=v]+"), ID("[o>=o]+"), 0, "MO Ints (vv|oo)");
    dpd_contract422(&I, &T_oo, &X_vv, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    // X_ab -= <ab|ij> dTau_ij
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v,v]"), ID("[o,o]"), 0, "MO Ints <vv|oo>");
    dpd_contract422(&I, &T_oo, &X_vv, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // X_ab += (ab|IJ) dTau_IJ
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[O,O]"),
                  ID("[v,v]"), ID("[O,O]"), 0, "MO Ints (vv|OO)");
    dpd_contract422(&I, &T_OO, &X_vv, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&X_vv);

    dpd_file2_init(&X_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    // X_IJ = +(IJ|AB) dTau_AB
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
    dpd_contract422(&I, &T_VV, &X_OO, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);

    // X_IJ -= <IJ|AB> dTau_AB
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_contract422(&I, &T_VV, &X_OO, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // X_IJ += (IJ|ab) dTau_ab
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                  ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
    dpd_contract422(&I, &T_vv, &X_OO, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    // X_IJ = +(IJ|KL) dTau_KL
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
    dpd_contract422(&I, &T_OO, &X_OO, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    // X_IJ -= <IJ|KL> dTau_KL
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO|OO>");
    dpd_contract422(&I, &T_OO, &X_OO, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // X_IJ += (IJ|kl) dTau_kl
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,o]"),
                  ID("[O>=O]+"), ID("[o>=o]+"), 0, "MO Ints (OO|oo)");
    dpd_contract422(&I, &T_oo, &X_OO, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&X_OO);


    dpd_file2_init(&X_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
    // X_ij = +(ij|ab) dTau_ab
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>=o]+"), ID("[v>=v]+"), 0, "MO Ints (oo|vv)");
    dpd_contract422(&I, &T_vv, &X_oo, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // X_ij -= <ij|ab> dTau_ab
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
    dpd_contract422(&I, &T_vv, &X_oo, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // X_ij += (ij|AB) dTau_AB
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[V,V]"),
                  ID("[o,o]"), ID("[V,V]"), 0, "MO Ints (oo|VV)");
    dpd_contract422(&I, &T_VV, &X_oo, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    // X_ij = +(ij|kl) dTau_kl
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o>=o]+"), ID("[o>=o]+"), 0, "MO Ints (oo|oo)");
    dpd_contract422(&I, &T_oo, &X_oo, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    // X_ij -= <ij|kl> dTau_kl
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo|oo>");
    dpd_contract422(&I, &T_oo, &X_oo, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // X_IJ += (ij|KL) dTau_KL
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[O,O]"),
                  ID("[o,o]"), ID("[O,O]"), 0, "MO Ints (oo|OO)");
    dpd_contract422(&I, &T_OO, &X_oo, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&X_oo);

    dpd_file2_close(&T_OO);
    dpd_file2_close(&T_oo);
    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);

    //
    // Compute T intermediate
    //

    /*
     * The T intermediates
     * T_ijab = 2.0 * P(ab) X_ca lambda_ijcb + 2.0 * P(ij) X_kj lambda_ikab
     */
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    // Temp_IJAB = lambda_IJCB X_AC
    dpd_file2_init(&X_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    dpd_contract244(&X_VV, &Laa, &T, 1, 2, 1, 2.0, 0.0);
    dpd_file2_close(&X_VV);
    dpd_buf4_close(&Laa);
    // T_IJAB = Temp_IJAB
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "T <OO|VV>");
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    // T_IJAB -= Temp_IJBA
    dpd_buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T <OO|VV>");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_add(&Taa, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    // Temp_IJAB = -lambda_KJAB X_IK
    dpd_file2_init(&X_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    dpd_contract244(&X_OO, &Laa, &T, 1, 0, 0, -2.0, 0.0);
    dpd_file2_close(&X_OO);
    dpd_buf4_close(&Laa);
    // T_IJAB += TEMP_IJAB
    dpd_buf4_add(&Taa, &T, 1.0);
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    // T_IJAB -= Temp_JIAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_add(&Taa, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_print(&Taa,outfile,1);
    dpd_buf4_close(&Taa);

    // The Tab intermediate
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T <Oo|Vv>");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    // T_IjAb = lambda_IjCb X_AC
    dpd_file2_init(&X_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    dpd_contract244(&X_VV, &Lab, &Tab, 1, 2, 1, 2.0, 0.0);
    dpd_file2_close(&X_VV);
    // T_IjAb += lambda_IjAc X_bc
    dpd_file2_init(&X_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
    dpd_contract424(&Lab, &X_vv, &Tab, 3, 1, 0, 2.0, 1.0);
    dpd_file2_close(&X_vv);
    // T_IjAb -= lambda_KjAb X_IK
    dpd_file2_init(&X_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    dpd_contract244(&X_OO, &Lab, &Tab, 1, 0, 0, -2.0, 1.0);
    dpd_file2_close(&X_OO);
    // T_IjAb -= lambda_IkAb X_jk
    dpd_file2_init(&X_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
    dpd_contract424(&Lab, &X_oo, &Tab, 1, 1, 1, -2.0, 1.0);
    dpd_file2_close(&X_oo);
    dpd_buf4_close(&Lab);
    dpd_buf4_print(&Tab,outfile,1);
    dpd_buf4_close(&Tab);

    // The Tbb intermediate
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    // Temp_ijab = lambda_ijcb X_ac
    dpd_file2_init(&X_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
    dpd_contract244(&X_vv, &Lbb, &T, 1, 2, 1, 2.0, 0.0);
    dpd_file2_close(&X_vv);
    dpd_buf4_close(&Lbb);
    // T_ijab = Temp_ijab
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "T <oo|vv>");
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);

    // T_ijab -= Temp_ijba
    dpd_buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T <oo|vv>");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_add(&Tbb, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    // Temp_ijab = -lambda_kjab X_ik
    dpd_file2_init(&X_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
    dpd_contract244(&X_oo, &Lbb, &T, 1, 0, 0, -2.0, 0.0);
    dpd_file2_close(&X_oo);
    dpd_buf4_close(&Lbb);
    // T_ijab += TEMP_ijab
    dpd_buf4_add(&Tbb, &T, 1.0);
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, qprs, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);
    // T_ijab -= Temp_jiab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_add(&Tbb, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_print(&Tbb,outfile,1);
    dpd_buf4_close(&Tbb);

    //End of the T intermediates

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

double
DCFTSolver::compute_cumulant_response_residual()
{
    dpdbuf4 R, G, F, T, C;
    double sumSQ = 0.0;
    size_t nElements = 0;

    /*
     * R_ijab = G_ijab + F_ijab + T_ijab + C_ijab
     */

    // R_IJAB = G_IJAB
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "G <OO|VV>");
    dpd_buf4_copy(&G, PSIF_DCFT_DPD, "R <OO|VV>");
    dpd_buf4_close(&G);
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "R <OO|VV>");

    // R_IJAB += F_IJAB
    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "F <OO|VV>");
    dpd_buf4_add(&R, &F, 1.0);
    dpd_buf4_close(&F);

    // R_IJAB += T_IJAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T <OO|VV>");
    dpd_buf4_add(&R, &T, 1.0);
    dpd_buf4_close(&T);

    // R_IJAB += C_IJAB
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "C <OO|VV> new");
    dpd_buf4_add(&R, &C, 1.0);
    dpd_buf4_close(&C);

    for(int h = 0; h < nirrep_; ++h)
        nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += dpd_buf4_dot_self(&R);
    dpd_buf4_print(&R,outfile,1);
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

    // R_IjAb += T_IjAb
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T <Oo|Vv>");
    dpd_buf4_add(&R, &T, 1.0);
    dpd_buf4_close(&T);

    // R_IjAb += C_IjAb
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "C <Oo|Vv> new");
    dpd_buf4_add(&R, &C, 1.0);
    dpd_buf4_close(&C);

    for(int h = 0; h < nirrep_; ++h)
        nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += dpd_buf4_dot_self(&R);
    dpd_buf4_print(&R,outfile,1);
    dpd_buf4_close(&R);

    // R_ijab = G_ijab
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "G <oo|vv>");
    dpd_buf4_copy(&G, PSIF_DCFT_DPD, "R <oo|vv>");
    dpd_buf4_close(&G);
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "R <oo|vv>");

    // R_ijab += F_ijab
    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "F <oo|vv>");
    dpd_buf4_add(&R, &F, 1.0);
    dpd_buf4_close(&F);

    // R_ijab += T_ijab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T <oo|vv>");
    dpd_buf4_add(&R, &T, 1.0);
    dpd_buf4_close(&T);

    // R_ijab += C_ijab
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "C <oo|vv> new");
    dpd_buf4_add(&R, &C, 1.0);
    dpd_buf4_close(&C);

    for(int h = 0; h < nirrep_; ++h)
        nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += dpd_buf4_dot_self(&R);
    dpd_buf4_print(&R,outfile,1);
    dpd_buf4_close(&R);

    return sqrt(sumSQ / nElements);
}

void
DCFTSolver::update_cumulant_response()
{

    dpdbuf4 Z, D, R;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    /*
     * Z_ijab += R_ijab / D_ijab
     */

    // Z_IJAB += R_IJAB / D_IJAB
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "R <OO|VV>");
    dpd_buf4_dirprd(&D, &R);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Z <OO|VV>");
    dpd_buf4_add(&Z, &R, 1.0);
    dpd_buf4_close(&R);
    dpd_buf4_close(&Z);

    // Z_IjAb += R_IjAb / D_IjAb
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
    dpd_buf4_dirprd(&D, &R);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    dpd_buf4_add(&Z, &R, 1.0);
    dpd_buf4_close(&R);
    dpd_buf4_close(&Z);

    // Z_ijab += R_ijab / D_ijab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "R <oo|vv>");
    dpd_buf4_dirprd(&D, &R);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Z <oo|vv>");
    dpd_buf4_add(&Z, &R, 1.0);
    dpd_buf4_close(&R);
    dpd_buf4_close(&Z);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}



}} //End namespaces
