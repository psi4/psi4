#include <libtrans/integraltransform.h>
#include <libpsio/psio.hpp>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <libdiis/diismanager.h>
#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

SharedMatrix DCFTSolver::compute_gradient()
{
    // Print out the header
    fprintf(outfile,"\n\n\t***********************************************************************************\n");
    fprintf(outfile,    "\t*                           DCFT Analytic Gradients Code                          *\n");
    fprintf(outfile,    "\t*                     by Alexander Sokolov and Andy Simmonett                     *\n");
    fprintf(outfile,    "\t***********************************************************************************\n\n");

    // Transform the one and two-electron integrals to the MO basis and write them into the DPD file
    gradient_init();

    if (!orbital_optimized_) {
        compute_gradient_dc();
    }
    else {
        compute_gradient_odc();
    }

    return SharedMatrix(new Matrix("NULL", 0, 0));

}

void
DCFTSolver::compute_gradient_dc() {

    bool responseDone = false;

    // Copy the current density cumulant and tau as a guess for cumulant response and perturbed tau
    response_guess();

    orbital_response_rms_ = 0.0;
    cumulant_response_rms_ = cumulant_convergence_;
    response_coupling_rms_ = 0.0;
    iter_ = 0;

    // Start two-step algorithm for solution of the response equations
    if(options_.get_str("RESPONSE_ALGORITHM") == "TWOSTEP"){
        fprintf(outfile, "\t*=================================================*\n"
                "\t* Cycle  RMS Orb. Resp.   RMS Cumul. Resp.   DIIS *\n"
                "\t*-------------------------------------------------*");

        // Start macro-iterations
        while(!responseDone && iter_++ < maxiter_){
            fprintf(outfile, "\n\t              *** Macro Iteration %d ***\n", iter_);

            // Solve the cumulant response equations iteratively
            if (iter_ > 1) {
                fprintf(outfile, "\t            Cumulant Response Iterations\n");
                iterate_cumulant_response();
            }

            // Compute the generalized densities for the MO Lagrangian
            compute_relaxed_density_OOOO();
            compute_relaxed_density_OOVV();
            compute_relaxed_density_OVOV();
            // Compute the OV block of MO Lagrangian
            compute_lagrangian_OV();
            // Compute the VO block of MO Lagrangian
            compute_lagrangian_VO();

            // Solve the orbital response equations iteratively
            fprintf(outfile, "\t            Orbital Response Iterations\n");
            iterate_orbital_response();

            // Compute terms that couple orbital and cumulant responses (C intermediate) and return RMS of their change
            response_coupling_rms_ = compute_response_coupling();

            // Check convergence
            if (response_coupling_rms_ < cumulant_threshold_) responseDone = true;

            fprintf(outfile, "\t   RMS of Response Coupling Change: %11.3E \n", response_coupling_rms_);
        }

        fprintf(outfile, "\t*=================================================*\n");
    }
    else {
        // Start the simultaneous algorithm for the solution of the response equations
        // Set up DIIS extrapolation
        dpdbuf4 Zaa, Zab, Zbb, Raa, Rab, Rbb;
        dpdfile2 zaa, zbb, raa, rbb;
        dpd_buf4_init(&Zaa, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");
        dpd_buf4_init(&Zab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
        dpd_buf4_init(&Zbb, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");
        dpd_file2_init(&zaa, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
        dpd_file2_init(&zbb, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
        DIISManager diisManager(maxdiis_, "DCFT DIIS orbital response vectors");
        diisManager.set_error_vector_size(5, DIISEntry::DPDFile2, &zaa,
                                          DIISEntry::DPDFile2, &zbb,
                                          DIISEntry::DPDBuf4, &Zaa,
                                          DIISEntry::DPDBuf4, &Zab,
                                          DIISEntry::DPDBuf4, &Zbb);
        diisManager.set_vector_size(5, DIISEntry::DPDFile2, &zaa,
                                    DIISEntry::DPDFile2, &zbb,
                                    DIISEntry::DPDBuf4, &Zaa,
                                    DIISEntry::DPDBuf4, &Zab,
                                    DIISEntry::DPDBuf4, &Zbb);
        dpd_buf4_close(&Zaa);
        dpd_buf4_close(&Zab);
        dpd_buf4_close(&Zbb);
        dpd_file2_close(&zaa);
        dpd_file2_close(&zbb);

        fprintf(outfile, "\t*==================================================================*\n"
                "\t* Cycle  RMS Orb. Resp.   RMS Cumul. Resp.   RMS Coupling     DIIS *\n"
                "\t*------------------------------------------------------------------*\n");

        // Start iterations
        while(!responseDone && iter_++ < maxiter_){
            std::string diisString;

            if (iter_ > 1) {
                // Update cumulant reponse from the change in C intermediate
                cumulant_response_guess();
                // Build perturbed tau and delta tau
                build_perturbed_tau();
                // Compute intermediates
                compute_cumulant_response_intermediates();
                // Compute cumulant response residual and its RMS
                cumulant_response_rms_ = compute_cumulant_response_residual();
                // Update the cumulant response
                update_cumulant_response();

                // Here's where DIIS kicks in
                if((cumulant_response_rms_ < diis_start_thresh_) && (orbital_response_rms_ < diis_start_thresh_)){
                    //Store the DIIS vectors
                    dpd_buf4_init(&Raa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                                  ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");
                    dpd_buf4_init(&Rab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                                  ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
                    dpd_buf4_init(&Rbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                                  ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");
                    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                                  ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");
                    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
                    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                                  ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");
                    dpd_file2_init(&raa, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "r <O|V>");
                    dpd_file2_init(&rbb, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "r <o|v>");
                    dpd_file2_init(&zaa, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
                    dpd_file2_init(&zbb, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");

                    if(diisManager.add_entry(10, &raa, &rbb, &Raa, &Rab, &Rbb, &zaa, &zbb, &Zaa, &Zab, &Zbb)){
                        diisString += "S";
                    }
                    if(diisManager.subspace_size() > mindiisvecs_){
                        diisString += "/E";
                        diisManager.extrapolate(5, &zaa, &zbb, &Zaa, &Zab, &Zbb);
                    }
                    dpd_file2_close(&zaa);
                    dpd_file2_close(&zbb);
                    dpd_file2_close(&raa);
                    dpd_file2_close(&rbb);
                    dpd_buf4_close(&Raa);
                    dpd_buf4_close(&Rab);
                    dpd_buf4_close(&Rbb);
                    dpd_buf4_close(&Zaa);
                    dpd_buf4_close(&Zab);
                    dpd_buf4_close(&Zbb);
                }
            }

            // Compute the generalized densities for the MO Lagrangian
            compute_relaxed_density_OOOO();
            compute_relaxed_density_OOVV();
            compute_relaxed_density_OVOV();
            // Compute the OV block of MO Lagrangian
            compute_lagrangian_OV();
            // Compute the VO block of MO Lagrangian
            compute_lagrangian_VO();

            // Compute guess for the orbital response matrix elements
            if (iter_ == 1) orbital_response_guess();

            // Compute orbital response intermediates
            compute_orbital_response_intermediates();

            // Update the orbital response
            orbital_response_rms_ = update_orbital_response();

            // Compute terms that couple orbital and cumulant responses (C intermediate) and return RMS of their change
            response_coupling_rms_ = compute_response_coupling();

            // Print iterative trace
            fprintf(outfile, "\t*%4d    %11.3E       %11.3E       %11.3E      %-4s *\n", iter_,
                    orbital_response_rms_, cumulant_response_rms_, response_coupling_rms_, diisString.c_str());

            // Check convergence
            if (orbital_response_rms_ < orbitals_threshold_ && cumulant_response_rms_ < cumulant_threshold_
                    && response_coupling_rms_ < cumulant_threshold_) responseDone = true;

        }

        fprintf(outfile, "\t*==================================================================*\n");
    }

    if (responseDone) fprintf(outfile, "\n\t   DCFT response equations converged.\n");
    else throw PSIEXCEPTION("DCFT response equations did not converge");

    // Compute the VVVV block of the relaxed TPDM
    compute_relaxed_density_VVVV();
    // Compute the OO block of MO Lagrangian
    compute_lagrangian_OO();
    // Compute the VV block of MO Lagrangian
    compute_lagrangian_VV();
    // Compute the energy-weighted density matrix
    compute_ewdm_dc();

}

void
DCFTSolver::compute_gradient_odc() {

    // Compute the VVVV block of the relaxed TPDM
    compute_unrelaxed_density_VVVV();
    fprintf(outfile, "\t Computing energy-weighted density matrix from one- and two-particle densities...\n");
    // Compute the OO block of MO Lagrangian
    compute_lagrangian_OO();
    // Compute the VV block of MO Lagrangian
    compute_lagrangian_VV();
    // Compute the energy-weighted density matrix
    compute_ewdm_odc();

}

void
DCFTSolver::gradient_init()
{

    // Allocate memory for the global objects
    aocc_ptau_ = SharedMatrix(new Matrix("MO basis Perturbed Tau (Alpha Occupied)", nirrep_, naoccpi_, naoccpi_));
    bocc_ptau_ = SharedMatrix(new Matrix("MO basis Perturbed Tau (Beta Occupied)", nirrep_, nboccpi_, nboccpi_));
    avir_ptau_ = SharedMatrix(new Matrix("MO basis Perturbed Tau (Alpha Virtual)", nirrep_, navirpi_, navirpi_));
    bvir_ptau_ = SharedMatrix(new Matrix("MO basis Perturbed Tau (Beta Virtual)", nirrep_, nbvirpi_, nbvirpi_));

    dpdbuf4 I;

    // Transform the two-electron integrals to the (VO|OO) and (OV|VV) subspaces in chemists' notation

    if ((options_.get_str("DCFT_FUNCTIONAL") == "DC-06" && options_.get_str("ALGORITHM") != "QC")
     || (options_.get_str("DCFT_FUNCTIONAL") == "DC-06" && options_.get_str("ALGORITHM") == "QC"
     && !options_.get_bool("QC_COUPLING") && options_.get_str("QC_TYPE") != "SIMULTANEOUS")) {
        _ints->transform_tei(MOSpace::vir, MOSpace::occ, MOSpace::occ, MOSpace::occ);
        _ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::occ);
        _ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::vir, MOSpace::vir);
        _ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    }

    // If the <VV|VV> integrals were not used for the energy computation (AO_BASIS = DISK) -> compute them for the gradients
    if (options_.get_str("AO_BASIS") == "DISK") _ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir);

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    /*
     * Re-sort the chemists' notation integrals to physisists' notation
     * (pq|rs) = <pr|qs>
     */

    if ((options_.get_str("DCFT_FUNCTIONAL") == "DC-06" && options_.get_str("ALGORITHM") != "QC")
            || (options_.get_str("DCFT_FUNCTIONAL") == "DC-06" && options_.get_str("ALGORITHM") == "QC"
                && !options_.get_bool("QC_COUPLING") && options_.get_str("QC_TYPE") != "SIMULTANEOUS")) {

        sort_OOOV_integrals();

        sort_OVVV_integrals();

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

    }

    // Hack for now. TODO: Implement AO_BASIS=DISK algorithm for gradients
    // (VV|VV)
    if(options_.get_str("AO_BASIS") == "DISK") sort_VVVV_integrals();

    // Transform one-electron integrals to the MO basis and store them in the DPD file

    if (!orbital_optimized_) transform_core_integrals();

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::response_guess()
{

    dpdbuf4 L;
    dpdfile2 T;

    // Copy the converged cumulant as a guess for the cumulant response

    // Z_IJAB = L_IJAB
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_copy(&L,PSIF_DCFT_DPD,"Z <OO|VV>");
    dpd_buf4_close(&L);

    // Z_IjAb = L_IjAb
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_copy(&L,PSIF_DCFT_DPD,"Z <Oo|Vv>");
    dpd_buf4_close(&L);

    // Z_ijab = L_ijab
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
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
DCFTSolver::compute_lagrangian_OV()
{
    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 G, I;
    dpdbuf4 L, W, U, Z;
    dpdfile2 X, H, pT;
    dpdfile2 T_VV, T_vv, dT_VV, dT_vv, pT_VV, pT_vv;
    dpdfile2 Y1_OV, Y1_ov, Y2_OV, Y2_ov;

    // X_OV: One-electron contributions

    // X_IA = H_IB pTau_BA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "pTau <V|V>");
    dpd_contract222(&H, &pT, &X, 0, 1, 1.0, 0.0);
    dpd_file2_close(&pT);
    dpd_file2_close(&H);
    dpd_file2_close(&X);

    // X_ia = H_ib pTau_ba
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "pTau <v|v>");
    dpd_contract222(&H, &pT, &X, 0, 1, 1.0, 0.0);
    dpd_file2_close(&pT);
    dpd_file2_close(&H);
    dpd_file2_close(&X);

    // X_OV: Two-electron contributions

    //
    // 2 * <OV||VV> Г_VVVV
    //

    // Compute contributions from VVVV density
    // 1. X_ia <-- <ib||cd> tau_ca dTau_db

    dpd_file2_init(&dT_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "dTau <V|V>");
    dpd_file2_init(&dT_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "dTau <v|v>");

    // Alpha contribution X0_IA
    dpd_file2_init(&Y1_OV, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Y1 <O|V>");

    // Y_IA = (IA|CD) dTau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
    dpd_contract422(&I, &dT_VV, &Y1_OV, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // Y_IA -= <IA|CD> dTau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    dpd_contract422(&I, &dT_VV, &Y1_OV, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_IA += (IA|cd) dTau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                  ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
    dpd_contract422(&I, &dT_vv, &Y1_OV, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);

    // X_IA = Y1_IC tau_CA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    dpd_contract222(&Y1_OV, &T_VV, &X, 0, 1, 1.0, 1.0);
    dpd_file2_close(&X);
    dpd_file2_close(&T_VV);

    dpd_file2_close(&Y1_OV);

    // Beta contribution X0_ia
    dpd_file2_init(&Y1_ov, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Y1 <o|v>");

    // Y_ib = (ib|cd) dTau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v>=v]+"), 0, "MO Ints (ov|vv)");
    dpd_contract422(&I, &dT_vv, &Y1_ov, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // Y_ib -= <ib|cd> dTau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
    dpd_contract422(&I, &dT_vv, &Y1_ov, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_ib += (ib|CD) dTau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[V>=V]+"),
                  ID("[o,v]"), ID("[V>=V]+"), 0, "MO Ints (ov|VV)");
    dpd_contract422(&I, &dT_VV, &Y1_ov, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);

    // X_ia = Y1_ic tau_ca
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");
    dpd_contract222(&Y1_ov, &T_vv, &X, 0, 1, 1.0, 1.0);
    dpd_file2_close(&X);
    dpd_file2_close(&T_vv);

    dpd_file2_close(&Y1_ov);

    dpd_file2_close(&dT_VV);
    dpd_file2_close(&dT_vv);

    // 2. X_ia <-- <ib||cd> tau_cb pTau_da

    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    // Alpha contribution X0_IA
    dpd_file2_init(&Y2_OV, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Y2 <O|V>");

    // Y2_IA = <IA|CD> Tau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    dpd_contract422(&I, &T_VV, &Y2_OV, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // Y2_IA -= (IA|CD) Tau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
    dpd_contract422(&I, &T_VV, &Y2_OV, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y2_IA -= (IA|cd) Tau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                  ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
    dpd_contract422(&I, &T_vv, &Y2_OV, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);

    // X_IA -= Y2_IC pTau_CA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_file2_init(&pT_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "pTau <V|V>");
    dpd_contract222(&Y2_OV, &pT_VV, &X, 0, 1, -1.0, 1.0);
    dpd_file2_close(&X);
    dpd_file2_close(&pT_VV);

    dpd_file2_close(&Y2_OV);

    // Beta contribution X0_ia
    dpd_file2_init(&Y2_ov, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Y2 <o|v>");

    // Y2_ib = <ib|cd> Tau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
    dpd_contract422(&I, &T_vv, &Y2_ov, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // Y2_ib -= (ib|cd) Tau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v>=v]+"), 0, "MO Ints (ov|vv)");
    dpd_contract422(&I, &T_vv, &Y2_ov, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y2_ib -= (ib|CD) Tau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[V>=V]+"),
                  ID("[o,v]"), ID("[V>=V]+"), 0, "MO Ints (ov|VV)");
    dpd_contract422(&I, &T_VV, &Y2_ov, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);

    // X_ia -= Y2_ic pTau_ca
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_file2_init(&pT_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "pTau <v|v>");
    dpd_contract222(&Y2_ov, &pT_vv, &X, 0, 1, -1.0, 1.0);
    dpd_file2_close(&X);
    dpd_file2_close(&pT_vv);

    dpd_file2_close(&Y2_ov);

    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);

    // 3. X_ia <-- 1/8 <ib||cd> Z_abkl lambda_klcd

    //  W_IBKL = <IB||CD> lambda_KLCD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V>V]-"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O>O]-"), 0, "W <OV|OO>");
    dpd_contract444(&I, &L, &W, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&W);


    //  W_KlIb = 2 lambda_KlCd <Ib|Cd>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "W <Oo|Ov>");
    dpd_contract444(&L, &I, &W, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&W);

    //  W_LkBi = 2 lambda_LkDc <Bi|Dc>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "W <Oo|Vo>");
    dpd_contract444(&L, &I, &W, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&W);

    //  W_ibkl = <ib||cd> lambda_klcd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v>v]-"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o>o]-"), 0, "W <ov|oo>");
    dpd_contract444(&I, &L, &W, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&W);


    // X_IA +=  1/8 W_IBKL Z_KLAB
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O>O]-"), 0, "W <OV|OO>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");

    dpd_contract442(&W, &Z, &X, 0, 2, 0.125, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);
    dpd_file2_close(&X);

    // X_IA +=  1/4 W_KlIb Z_KlAb
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "W <Oo|Ov>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");

    dpd_contract442(&W, &Z, &X, 2, 2, 0.25, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);
    dpd_file2_close(&X);

    // X_ia +=  1/4 W_LkBi Z_LkBa
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "W <Oo|Vo>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");

    dpd_contract442(&W, &Z, &X, 3, 3, 0.25, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);
    dpd_file2_close(&X);

    // X_ia += 1/8 W_ibkl Z_klab
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o>o]-"), 0, "W <ov|oo>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");

    dpd_contract442(&W, &Z, &X, 0, 2, 0.125, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&Z);

    dpd_file2_close(&X);

    // 4. X_ia <-- 1/8 <ib||cd> Z_cdkl lambda_klab

    //  U_IBKL = <IB||CD> Z_KLCD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V>V]-"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");
    dpd_buf4_init(&U, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O>O]-"), 0, "U <OV|OO>");
    dpd_contract444(&I, &Z, &U, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&U);


    //  U_KlIb = 2.0 Z_KlCd <Ib|Cd>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    dpd_buf4_init(&U, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "U <Oo|Ov>");
    dpd_contract444(&Z, &I, &U, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&U);

    //  U_LkBi = 2.0 Z_LkDc <Bi|Dc>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    dpd_buf4_init(&U, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "U <Oo|Vo>");
    dpd_contract444(&Z, &I, &U, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&U);

    //  U_ibkl = <ib||cd> Z_klcd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v>v]-"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");
    dpd_buf4_init(&U, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o>o]-"), 0, "U <ov|oo>");
    dpd_contract444(&I, &Z, &U, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&U);


    // X_IA +=  1/8 U_IBKL lambda_KLAB
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&U, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O>O]-"), 0, "U <OV|OO>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");

    dpd_contract442(&U, &L, &X, 0, 2, 0.125, 1.0);
    dpd_buf4_close(&U);
    dpd_buf4_close(&L);
    dpd_file2_close(&X);

    // X_IA +=  1/4 U_KlIb lambda_KlAb
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&U, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "U <Oo|Ov>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");

    dpd_contract442(&U, &L, &X, 2, 2, 0.25, 1.0);
    dpd_buf4_close(&U);
    dpd_buf4_close(&L);
    dpd_file2_close(&X);

    // X_ia +=  1/4 U_LkBi L_LkBa
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&U, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "U <Oo|Vo>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");

    dpd_contract442(&U, &L, &X, 3, 3, 0.25, 1.0);
    dpd_buf4_close(&U);
    dpd_buf4_close(&L);
    dpd_file2_close(&X);

    // X_ia +=  1/8 U_ibkl L_klab
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&U, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o>o]-"), 0, "U <ov|oo>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");

    dpd_contract442(&U, &L, &X, 0, 2, 0.125, 1.0);
    dpd_buf4_close(&U);
    dpd_buf4_close(&L);
    dpd_file2_close(&X);

    //
    // <OO||OV> Г_OOVV
    //

    // X_IA += <BI||JK> Г_BAJK
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 1, "MO Ints <VO|OO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V>V]-"), ID("[O>O]-"), 0, "Gamma <VV|OO>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
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
    dpd_file2_close(&X);

    // X_ia += <bi||jk> Г_bajk
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 1, "MO Ints <vo|oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v>v]-"), ID("[o>o]-"), 0, "Gamma <vv|oo>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
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
        #pragma omp parallel for
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
        #pragma omp parallel for
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
                  ID("[O>O]-"), ID("[O>O]-"), 0, "Gamma <OO|OO>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
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
    dpd_file2_close(&X);

    // X_ai += 2 * <aj||kl> Г_ijkl
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 1, "MO Ints <vo|oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o>o]-"), ID("[o>o]-"), 0, "Gamma <oo|oo>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
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
    dpd_file2_close(&X);

    //
    // <VO||VV> Г_OOVV
    //

    // X_AI += <JA||BC> Г_JIBC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Gamma <OO|VV>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
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
    dpd_file2_close(&X);

    // X_ai += <ja||bc> Г_jibc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Gamma <oo|vv>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
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
    dpd_file2_close(&X);

    psio_->close(PSIF_DCFT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::iterate_orbital_response()
{

    // Compute guess for the orbital response matrix elements
    if (iter_ == 1) orbital_response_guess();

    bool converged = false;

    // Initialize DIIS
    dpdfile2 zaa, zbb, raa, rbb;
    dpd_file2_init(&zaa, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_file2_init(&zbb, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");
    DIISManager ZiaDiisManager(maxdiis_, "DCFT DIIS Orbital Z",DIISManager::LargestError,DIISManager::InCore);
    ZiaDiisManager.set_error_vector_size(2, DIISEntry::DPDFile2, &zaa,
                                            DIISEntry::DPDFile2, &zbb);
    ZiaDiisManager.set_vector_size(2, DIISEntry::DPDFile2, &zaa,
                                      DIISEntry::DPDFile2, &zbb);
    dpd_file2_close(&zaa);
    dpd_file2_close(&zbb);

     // Start iterations
    int cycle = 0;
    do {
        cycle++;
        std::string diisString;

        // Compute intermediates
        compute_orbital_response_intermediates();

        // Update the orbital response
        orbital_response_rms_ = update_orbital_response();

        // Here's where DIIS kicks in
        if(orbital_response_rms_ < diis_start_thresh_){
            //Store the DIIS vectors
            dpd_file2_init(&raa, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "r <O|V>");
            dpd_file2_init(&rbb, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "r <o|v>");
            dpd_file2_init(&zaa, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
            dpd_file2_init(&zbb, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");

            if(ZiaDiisManager.add_entry(4, &raa, &rbb, &zaa, &zbb)){
                diisString += "S";
            }
            // Extrapolate orbital response
            if(ZiaDiisManager.subspace_size() >= mindiisvecs_ && maxdiis_ > 0){
                diisString += "/E";
                ZiaDiisManager.extrapolate(2, &zaa, &zbb);
            }
            dpd_file2_close(&zaa);
            dpd_file2_close(&zbb);
            dpd_file2_close(&raa);
            dpd_file2_close(&rbb);
        }

        // Check convergence
        converged = (fabs(orbital_response_rms_) < fabs(orbitals_threshold_));

        // Print iterative trace
        fprintf(outfile, "\t*%4d    %11.3E       %11.3E       %-4s *\n", cycle,
                orbital_response_rms_, cumulant_response_rms_, diisString.c_str());

        // Termination condition
        if (converged || cycle >= maxiter_) break;

    } while (true);

    if (!converged) throw PSIEXCEPTION("DCFT orbital response equations did not converge");

}

void
DCFTSolver::orbital_response_guess()
{

    dpdfile2 Xia, Xai, zia;

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
        #pragma omp parallel for
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int a = 0 ; a < navirpi_[h]; ++a){
                double value_dX = 2.0 * (Xia.matrix[h][i][a] - Xai.matrix[h][a][i]);
                zia.matrix[h][i][a] = value_dX / (moFa_->get(h, a + naoccpi_[h], a + naoccpi_[h]) - moFa_->get(h, i, i));
            }
        }
    }
    dpd_file2_mat_wrt(&zia);
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
        #pragma omp parallel for
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int a = 0 ; a < nbvirpi_[h]; ++a){
                double value_dX = 2.0 * (Xia.matrix[h][i][a] - Xai.matrix[h][a][i]);
                zia.matrix[h][i][a] = value_dX / (moFb_->get(h, a + nboccpi_[h], a + nboccpi_[h]) - moFb_->get(h, i, i));
            }
        }
    }
    dpd_file2_mat_wrt(&zia);
    dpd_file2_close(&zia);
    dpd_file2_close(&Xai);
    dpd_file2_close(&Xia);

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
        #pragma omp parallel for
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

// Returns RMS of the orbital response vector
double
DCFTSolver::update_orbital_response()
{

    dpdfile2 X_ia, X_ai, z_ia, zI_ai, zI_ia, r_ia;
    SharedMatrix a_ria (new Matrix("MO basis Orbital Response Residual (Alpha)", nirrep_, naoccpi_, navirpi_));
    SharedMatrix b_ria (new Matrix("MO basis Orbital Response Residual (Beta)", nirrep_, nboccpi_, nbvirpi_));

    // Alpha spin
    dpd_file2_init(&zI_ia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "zI <O|V>");
    dpd_file2_init(&zI_ai, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "zI <V|O>");
    dpd_file2_init(&X_ia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_file2_init(&X_ai, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_file2_init(&z_ia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");
    dpd_file2_init(&r_ia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "r <O|V>");

    dpd_file2_mat_init(&zI_ai);
    dpd_file2_mat_init(&zI_ia);
    dpd_file2_mat_init(&X_ia);
    dpd_file2_mat_init(&X_ai);
    dpd_file2_mat_init(&z_ia);
    dpd_file2_mat_init(&r_ia);

    dpd_file2_mat_rd(&zI_ai);
    dpd_file2_mat_rd(&zI_ia);
    dpd_file2_mat_rd(&X_ia);
    dpd_file2_mat_rd(&X_ai);
    dpd_file2_mat_rd(&z_ia);

    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int a = 0 ; a < navirpi_[h]; ++a){
                double value = 0.0;
                for(int j = 0 ; j < naoccpi_[h]; ++j){
                    value -= (zI_ai.matrix[h][a][j] + zI_ia.matrix[h][j][a]) * (aocc_tau_->get(h,i,j) + kappa_mo_a_->get(h,i,j));
                    value += z_ia.matrix[h][j][a] * moFa_->get(h, j, i);
                }
                for(int b = 0 ; b < navirpi_[h]; ++b){
                    value += (zI_ai.matrix[h][b][i] + zI_ia.matrix[h][i][b]) * (avir_tau_->get(h,a,b));
                    value -= z_ia.matrix[h][i][b] * moFa_->get(h, b + naoccpi_[h], a + naoccpi_[h]);
                }
                value += 2.0 * (X_ia.matrix[h][i][a] - X_ai.matrix[h][a][i]);
                a_ria->set(h, i, a, value);
                r_ia.matrix[h][i][a] = value;
                z_ia.matrix[h][i][a] += value / (moFa_->get(h, a + naoccpi_[h], a + naoccpi_[h]) - moFa_->get(h, i, i));
            }
        }
    }
    dpd_file2_mat_wrt(&z_ia);
    dpd_file2_mat_wrt(&r_ia);
    dpd_file2_close(&z_ia);
    dpd_file2_close(&r_ia);
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
    dpd_file2_init(&r_ia, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "r <o|v>");

    dpd_file2_mat_init(&zI_ai);
    dpd_file2_mat_init(&zI_ia);
    dpd_file2_mat_init(&X_ia);
    dpd_file2_mat_init(&X_ai);
    dpd_file2_mat_init(&z_ia);
    dpd_file2_mat_init(&r_ia);

    dpd_file2_mat_rd(&zI_ai);
    dpd_file2_mat_rd(&zI_ia);
    dpd_file2_mat_rd(&X_ia);
    dpd_file2_mat_rd(&X_ai);
    dpd_file2_mat_rd(&z_ia);

    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int a = 0 ; a < nbvirpi_[h]; ++a){
                double value = 0.0;
                for(int j = 0 ; j < nboccpi_[h]; ++j){
                    value -= (zI_ai.matrix[h][a][j] + zI_ia.matrix[h][j][a]) * (bocc_tau_->get(h,i,j) + kappa_mo_b_->get(h,i,j));
                    value += z_ia.matrix[h][j][a] * moFb_->get(h, j, i);
                }
                for(int b = 0 ; b < nbvirpi_[h]; ++b){
                    value += (zI_ai.matrix[h][b][i] + zI_ia.matrix[h][i][b]) * (bvir_tau_->get(h,a,b));
                    value -= z_ia.matrix[h][i][b] * moFb_->get(h, b + nboccpi_[h], a + nboccpi_[h]);
                }
                value += 2.0 * (X_ia.matrix[h][i][a] - X_ai.matrix[h][a][i]);
                b_ria->set(h, i, a, value);
                r_ia.matrix[h][i][a] = value;
                z_ia.matrix[h][i][a] += value / (moFb_->get(h, a + nboccpi_[h], a + nboccpi_[h]) - moFb_->get(h, i, i));
            }
        }
    }
    dpd_file2_mat_wrt(&z_ia);
    dpd_file2_mat_wrt(&r_ia);
    dpd_file2_close(&z_ia);
    dpd_file2_close(&r_ia);
    dpd_file2_close(&X_ai);
    dpd_file2_close(&X_ia);
    dpd_file2_close(&zI_ai);
    dpd_file2_close(&zI_ia);

    // Compute RMS
    return (a_ria->rms() + b_ria->rms());

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
        #pragma omp parallel for
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                zIsym.matrix[h][i][j] = zIsym.matrix[h][j][i] = zI.matrix[h][i][j] + zI.matrix[h][j][i];
            }
        }
    }
    dpd_file2_mat_wrt(&zIsym);
    dpd_file2_close(&zI);
    dpd_file2_close(&zIsym);

    // zI <o|o>
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "zI <o|o>");
    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "zI <o|o> sym");
    dpd_file2_mat_init(&zI);
    dpd_file2_mat_init(&zIsym);
    dpd_file2_mat_rd(&zI);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                zIsym.matrix[h][i][j] = zIsym.matrix[h][j][i] = zI.matrix[h][i][j] + zI.matrix[h][j][i];
            }
        }
    }
    dpd_file2_mat_wrt(&zIsym);
    dpd_file2_close(&zI);
    dpd_file2_close(&zIsym);

    // zI <V|V>
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "zI <V|V>");
    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "zI <V|V> sym");
    dpd_file2_mat_init(&zI);
    dpd_file2_mat_init(&zIsym);
    dpd_file2_mat_rd(&zI);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0 ; i < navirpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                zIsym.matrix[h][i][j] = zIsym.matrix[h][j][i] = zI.matrix[h][i][j] + zI.matrix[h][j][i];
            }
        }
    }
    dpd_file2_mat_wrt(&zIsym);
    dpd_file2_close(&zI);
    dpd_file2_close(&zIsym);

    // zI <v|v>
    dpd_file2_init(&zI, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "zI <v|v>");
    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "zI <v|v> sym");
    dpd_file2_mat_init(&zI);
    dpd_file2_mat_init(&zIsym);
    dpd_file2_mat_rd(&zI);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0 ; i < nbvirpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                zIsym.matrix[h][i][j] = zIsym.matrix[h][j][i] = zI.matrix[h][i][j] + zI.matrix[h][j][i];
            }
        }
    }
    dpd_file2_mat_wrt(&zIsym);
    dpd_file2_close(&zI);
    dpd_file2_close(&zIsym);

    // Save the old response coupling terms (C intermediate)
    if (iter_ > 1) {
        dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "C <OO|VV> new");
        dpd_buf4_copy(&C, PSIF_DCFT_DPD, "C <OO|VV> old");
        dpd_buf4_close(&C);

        dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "C <Oo|Vv> new");
        dpd_buf4_copy(&C, PSIF_DCFT_DPD, "C <Oo|Vv> old");
        dpd_buf4_close(&C);

        dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "C <oo|vv> new");
        dpd_buf4_copy(&C, PSIF_DCFT_DPD, "C <oo|vv> old");
        dpd_buf4_close(&C);

    }

    // Compute the new response coupling terms (C intermediate)

    //
    // C <OO|VV>
    //

    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "zI <V|V> sym");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "C <OO|VV> new");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    // The first pair index of the result of dpd_contract244 contains the file2 index.
    // The second pair index of the result of dpd_contract424 contains the file2 index.
    // This is how dpd_contract244 works: AD * IJDB -> ABIJ ?-(t)-> IJAB. If BD * IJDA -> BAIJ -(t)-> IJBA
    // This is how dpd_contract424 works: IJAD * BD -> IJAB ?-(t)-> ABIJ. If IJBD * AD -> IJBA ?-(t)-> BAIJ
    // Temp_IJAB = zIsym_AC lambda_IJCB
    dpd_contract244(&zIsym, &L, &T, 1, 2, 1, 1.0, 0.0);
    // C_IJAB = Temp_IJAB
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "C <OO|VV> new");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
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
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
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
    dpd_buf4_close(&C);
    dpd_buf4_close(&L);

    //
    // C <oo|vv>
    //

    dpd_file2_init(&zIsym, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "zI <v|v> sym");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "C <oo|vv> new");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    // Temp_ijab = zIsym_ac lambda_ijcb
    dpd_contract244(&zIsym, &L, &T, 1, 2, 1, 1.0, 0.0);
    // C_ijab = Temp_ijab
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "C <oo|vv> new");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
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
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
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
    dpd_buf4_close(&C);

    // Compute the RMS of the change in the response coupling terms (RMS (Cnew - Cold))

    int nElements = 0;
    double sumSQ = 0.0;

    // dC <OO|VV>
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "C <OO|VV> new");
    dpd_buf4_copy(&C, PSIF_DCFT_DPD, "C <OO|VV> new - old");
    dpd_buf4_close(&C);

    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "C <OO|VV> new - old");
    if (iter_ > 1) {
        dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O>O]-"), ID("[V>V]-"), 0, "C <OO|VV> old");
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
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "C <oo|vv> new");
    dpd_buf4_copy(&C, PSIF_DCFT_DPD, "C <oo|vv> new - old");
    dpd_buf4_close(&C);

    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "C <oo|vv> new - old");
    if (iter_ > 1) {
        dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o>o]-"), ID("[v>v]-"), 0, "C <oo|vv> old");
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

    // Set up DIIS extrapolation
    dpdbuf4 Zaa, Zab, Zbb, Raa, Rab, Rbb;
    dpd_buf4_init(&Zaa, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");
    dpd_buf4_init(&Zab, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    dpd_buf4_init(&Zbb, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");
    DIISManager ZDiisManager(maxdiis_, "DCFT DIIS Z",DIISManager::LargestError,DIISManager::InCore);
    ZDiisManager.set_error_vector_size(3, DIISEntry::DPDBuf4, &Zaa,
                                          DIISEntry::DPDBuf4, &Zab,
                                          DIISEntry::DPDBuf4, &Zbb);
    ZDiisManager.set_vector_size(3, DIISEntry::DPDBuf4, &Zaa,
                                    DIISEntry::DPDBuf4, &Zab,
                                    DIISEntry::DPDBuf4, &Zbb);
    dpd_buf4_close(&Zaa);
    dpd_buf4_close(&Zab);
    dpd_buf4_close(&Zbb);

    // Iteratively solve for cumulant reponse
    bool converged = false;

    int cycle = 0;
    do {
        cycle++;
        std::string diisString;

        // Build perturbed tau and delta tau
        build_perturbed_tau();

        // Compute intermediates
        compute_cumulant_response_intermediates();

        // Compute cumulant response residual and its RMS
        cumulant_response_rms_ = compute_cumulant_response_residual();

        // Update the cumulant response
        update_cumulant_response();

        // Here's where DIIS kicks in
        if(cumulant_response_rms_ < diis_start_thresh_){
            //Store the DIIS vectors
            dpd_buf4_init(&Raa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                          ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");
            dpd_buf4_init(&Rab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                          ID("[O,o]"), ID("[V,v]"), 0, "R <Oo|Vv>");
            dpd_buf4_init(&Rbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                          ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");
            dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                          ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");
            dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                          ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
            dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                          ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");

            if(ZDiisManager.add_entry(6, &Raa, &Rab, &Rbb, &Zaa, &Zab, &Zbb)){
                diisString += "S";
            }
            // Extrapolate cumulant response
            if(ZDiisManager.subspace_size() >= mindiisvecs_ && maxdiis_ > 0){
                diisString += "/E";
                ZDiisManager.extrapolate(3, &Zaa, &Zab, &Zbb);
            }
            dpd_buf4_close(&Raa);
            dpd_buf4_close(&Rab);
            dpd_buf4_close(&Rbb);
            dpd_buf4_close(&Zaa);
            dpd_buf4_close(&Zab);
            dpd_buf4_close(&Zbb);
        }

        // Check the convergence
        converged = (fabs(cumulant_response_rms_) < fabs(cumulant_threshold_));

        // Print iterative trace
        fprintf(outfile, "\t*%4d    %11.3E       %11.3E       %-4s *\n", cycle, orbital_response_rms_, cumulant_response_rms_, diisString.c_str());

        // Termination condition
        if (converged || cycle >= maxiter_) break;

    } while (true);

    if (!converged) throw PSIEXCEPTION("DCFT cumulant response equations did not converge");

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
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "D <OO|VV>");
    dpd_buf4_init(&dC, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "C <OO|VV> new - old");
    dpd_buf4_dirprd(&D, &dC);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");
    dpd_buf4_add(&Z, &dC, 1.0);
    dpd_buf4_close(&dC);
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
    dpd_buf4_close(&Z);

    // Z_ijab += dC_ijab / D_ijab
    dpd_buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>=o]+"), ID("[v>=v]+"), 0, "D <oo|vv>");
    dpd_buf4_init(&dC, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "C <oo|vv> new - old");
    dpd_buf4_dirprd(&D, &dC);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");
    dpd_buf4_add(&Z, &dC, 1.0);
    dpd_buf4_close(&dC);
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
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");

    // pTau_IJ = -1/2 Lambda_IKAB Z_JKAB
    dpd_contract442(&L, &Z, &pT_OO, 0, 0, -0.5, 0.0);

    // pTau_AB = +1/2 Lambda_IJAC Z_IJBC
    dpd_contract442(&L, &Z, &pT_VV, 2, 2, 0.5, 0.0);

    dpd_buf4_close(&Z);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");

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

    // Put perturbed tau from disk to memory
    dpd_file2_init(&pT_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "pTau <O|O>");
    dpd_file2_init(&pT_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "pTau <o|o>");
    dpd_file2_init(&pT_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "pTau <V|V>");
    dpd_file2_init(&pT_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "pTau <v|v>");

    dpd_file2_mat_init(&pT_OO);
    dpd_file2_mat_init(&pT_oo);
    dpd_file2_mat_init(&pT_VV);
    dpd_file2_mat_init(&pT_vv);

    dpd_file2_mat_rd(&pT_OO);
    dpd_file2_mat_rd(&pT_oo);
    dpd_file2_mat_rd(&pT_VV);
    dpd_file2_mat_rd(&pT_vv);

    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int j = 0; j < naoccpi_[h]; ++j){
                aocc_ptau_->set(h, i, j, pT_OO.matrix[h][i][j]);
            }
        }
        for(int a = 0; a < navirpi_[h]; ++a){
            for(int b = 0; b < navirpi_[h]; ++b){
                avir_ptau_->set(h, a, b, pT_VV.matrix[h][a][b]);
            }
        }
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int j = 0; j < nboccpi_[h]; ++j){
                bocc_ptau_->set(h, i, j, pT_oo.matrix[h][i][j]);
            }
        }
        for(int a = 0; a < nbvirpi_[h]; ++a){
            for(int b = 0; b < nbvirpi_[h]; ++b){
                bvir_ptau_->set(h, a, b, pT_vv.matrix[h][a][b]);
            }
        }
    }

    dpd_file2_close(&pT_OO);
    dpd_file2_close(&pT_oo);
    dpd_file2_close(&pT_VV);
    dpd_file2_close(&pT_vv);

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

    dpdfile2 F_OO, F_oo, F_VV, F_vv, Y_OO, Y_oo, Y_VV, Y_vv, T_OO, T_oo, T_VV, T_vv, tmp;
    dpdbuf4 I, Z, G, T, F, Taa, Tab, Tbb, Laa, Lab, Lbb, Zaa, Zab, Zbb;

    //
    // Compute G intermediates
    //

    /*
     * G_ijab = <ij||ab>
     */
    // G_IJAB = <IJ||AB>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "G <OO|VV>");
    dpd_buf4_close(&I);

    // G_IjAb = <Ij|Ab>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "G <Oo|Vv>");
    dpd_buf4_close(&I);

    // G_ijab = <ij||ab>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    dpd_buf4_copy(&I, PSIF_DCFT_DPD, "G <oo|vv>");
    dpd_buf4_close(&I);

    /*
     * G_ijab += 1/2 Sum_cd gbar_cdab Z_ijcd
     */
    // G_IJAB += 1/2 Sum_CD gbar_CDAB Z_IJCD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V>V]-"), ID("[V>V]-"),
                  ID("[V,V]"), ID("[V,V]"), 1, "MO Ints <VV|VV>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_contract444(&Z, &I, &G, 0, 0, 1.0, 1.0);
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
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v>v]-"), ID("[v>v]-"),
                  ID("[v,v]"), ID("[v,v]"), 1, "MO Ints <vv|vv>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_contract444(&Z, &I, &G, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&G);

    /*
     * G_ijab += 1/2 Sum_kl gbar_ijkl Z_klab
     */
    // G_IJAB += 1/2 Sum_KL gbar_IJKL Z_KLAB
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O>O]-"), ID("[O>O]-"),
            ID("[O,O]"), ID("[O,O]"), 1, "MO Ints <OO|OO>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_contract444(&I, &Z, &G, 0, 1, 1.0, 1.0);
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
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o>o]-"), ID("[o>o]-"),
                  ID("[o,o]"), ID("[o,o]"), 1, "MO Ints <oo|oo>");
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_contract444(&I, &Z, &G, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&G);


    /*
     * G_ijab -= P(ij)P(ab) Sum_kc gbar_jckb Z_ikac
     */
    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");
    dpd_buf4_sort(&Zaa, PSIF_DCFT_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "Z (OV|OV)");
    dpd_buf4_close(&Zaa);
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    dpd_buf4_sort(&Zab, PSIF_DCFT_DPD, psqr, ID("[O,v]"), ID("[o,V]"), "Z (Ov|oV)");
    dpd_buf4_close(&Zab);
    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");
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
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");

    // T_IAJB -= Sum_KC Z_IAKC gbar_JBKC
    dpd_contract444(&Zaa, &I, &Taa, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC gbar_IAKC Z_KCjb
    dpd_contract444(&I, &Zab, &Tab, 0, 1, -1.0, 1.0);

    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");

    // T_IAJB += Sum_KC Z_IAKC (JB|KC)
    dpd_contract444(&Zaa, &I, &Taa, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC (JB|KC) Z_KCjb
    dpd_contract444(&I, &Zab, &Tab, 0, 1, 1.0, 1.0);
    dpd_buf4_close(&I);


    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");

    // T_iajb -= Sum_kc Z_iakc gbar_jbkc
    dpd_contract444(&Zbb, &I, &Tbb, 0, 0, -1.0, 1.0);

    // T_IAjb -= Sum_KC Z_IAkc gbar_jbkc
    dpd_contract444(&Zab, &I, &Tab, 0, 0, -1.0, 1.0);

    dpd_buf4_close(&I);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");

    // T_iajb += Sum_kc Z_iakc (JB|KC)
    dpd_contract444(&Zbb, &I, &Tbb, 0, 0, 1.0, 1.0);

    // T_IAjb += Sum_KC Z_IAkc (kc|jb)
    dpd_contract444(&Zab, &I, &Tab, 0, 0, 1.0, 1.0);
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
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
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
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_IJAB -> T_IJBA
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // G_IJAB -= T_IJBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_IJAB -> T_JIBA
    dpd_buf4_sort(&Taa, PSIF_DCFT_DPD, qpsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    // G_IJAB += T_JIBA
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "G <OO|VV>");
    dpd_buf4_add(&G, &T, 1.0);
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
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_iajb -> T_ijab
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, prqs, ID("[o,o]"), ID("[v,v]"), "Temp <oo|vv>");
    dpd_buf4_close(&Tbb);
    // G_ijab += T_ijab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
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
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_ijab -> T_ijba
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab -= T_ijba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, -1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&T);

    // T_ijab -> T_jiba
    dpd_buf4_sort(&Tbb, PSIF_DCFT_DPD, qpsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    // G_ijab += T_jiba
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "G <oo|vv>");
    dpd_buf4_add(&G, &T, 1.0);
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
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");
    // Temp_IJAB = Z_IJCB F_AC
    dpd_file2_init(&F_VV, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    dpd_contract244(&F_VV, &Zaa, &T, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&F_VV);
    dpd_buf4_close(&Zaa);
    dpd_buf4_close(&T);
    // F_IJAB = Temp_IJAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "F <OO|VV>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "F <OO|VV>");
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
              ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");
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
    dpd_buf4_close(&F);

    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");
    // Temp_ijab = Z_ijcb F_ac
    dpd_file2_init(&F_vv, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "F <v|v>");
    dpd_contract244(&F_vv, &Zbb, &T, 1, 2, 1, 1.0, 0.0);
    dpd_file2_close(&F_vv);
    dpd_buf4_close(&Zbb);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    // F_ijab = Temp_ijab
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "F <oo|vv>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&F, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "F <oo|vv>");
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
              ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");
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
    dpd_buf4_close(&F);

    //
    // Compute Y intermediates
    //

    dpd_file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "dTau <O|O>");
    dpd_file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "dTau <o|o>");
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "dTau <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "dTau <v|v>");

    dpd_file2_init(&Y_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Y <V|V>");

    // Y_AB = (AB|CD) dTau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
    dpd_contract422(&I, &T_VV, &Y_VV, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // Y_AB -= <AB|CD> dTau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
    dpd_contract422(&I, &T_VV, &Y_VV, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_AB += (AB|cd) dTau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[v,v]"),
                  ID("[V>=V]+"), ID("[v>=v]+"), 0, "MO Ints (VV|vv)");
    dpd_contract422(&I, &T_vv, &Y_VV, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);

    // Y_AB = +(AB|IJ) dTau_IJ
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V>=V]+"), ID("[O>=O]+"), 0, "MO Ints (VV|OO)");
    dpd_contract422(&I, &T_OO, &Y_VV, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_AB -= <AB|IJ> dTau_IJ
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V,V]"), ID("[O,O]"), 0, "MO Ints <VV|OO>");
    dpd_contract422(&I, &T_OO, &Y_VV, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_AB += (AB|ij) dTau_ij
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                  ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
    dpd_contract422(&I, &T_oo, &Y_VV, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&Y_VV);

    dpd_file2_init(&Y_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Y <v|v>");

    // Y_ab = +(ab|cd) dTau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v>=v]+"), ID("[v>=v]+"), 0, "MO Ints (vv|vv)");
    dpd_contract422(&I, &T_vv, &Y_vv, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // Y_ab -= <ab|cd> dTau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv|vv>");
    dpd_contract422(&I, &T_vv, &Y_vv, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_ab += (ab|CD) dTau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[V,V]"),
                  ID("[v,v]"), ID("[V,V]"), 0, "MO Ints (vv|VV)");
    dpd_contract422(&I, &T_VV, &Y_vv, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);

    // Y_ab = +(ab|ij) dTau_ij
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v>=v]+"), ID("[o>=o]+"), 0, "MO Ints (vv|oo)");
    dpd_contract422(&I, &T_oo, &Y_vv, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_ab -= <ab|ij> dTau_ij
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v,v]"), ID("[o,o]"), 0, "MO Ints <vv|oo>");
    dpd_contract422(&I, &T_oo, &Y_vv, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_ab += (ab|IJ) dTau_IJ
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[O,O]"),
                  ID("[v,v]"), ID("[O,O]"), 0, "MO Ints (vv|OO)");
    dpd_contract422(&I, &T_OO, &Y_vv, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&Y_vv);

    dpd_file2_init(&Y_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Y <O|O>");
    // Y_IJ = +(IJ|AB) dTau_AB
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
    dpd_contract422(&I, &T_VV, &Y_OO, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);

    // Y_IJ -= <IJ|AB> dTau_AB
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    dpd_contract422(&I, &T_VV, &Y_OO, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_IJ += (IJ|ab) dTau_ab
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                  ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
    dpd_contract422(&I, &T_vv, &Y_OO, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_IJ = +(IJ|KL) dTau_KL
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
    dpd_contract422(&I, &T_OO, &Y_OO, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_IJ -= <IJ|KL> dTau_KL
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO|OO>");
    dpd_contract422(&I, &T_OO, &Y_OO, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_IJ += (IJ|kl) dTau_kl
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,o]"),
                  ID("[O>=O]+"), ID("[o>=o]+"), 0, "MO Ints (OO|oo)");
    dpd_contract422(&I, &T_oo, &Y_OO, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&Y_OO);


    dpd_file2_init(&Y_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Y <o|o>");
    // Y_ij = +(ij|ab) dTau_ab
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>=o]+"), ID("[v>=v]+"), 0, "MO Ints (oo|vv)");
    dpd_contract422(&I, &T_vv, &Y_oo, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // Y_ij -= <ij|ab> dTau_ab
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
    dpd_contract422(&I, &T_vv, &Y_oo, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_ij += (ij|AB) dTau_AB
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[V,V]"),
                  ID("[o,o]"), ID("[V,V]"), 0, "MO Ints (oo|VV)");
    dpd_contract422(&I, &T_VV, &Y_oo, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_ij = +(ij|kl) dTau_kl
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o>=o]+"), ID("[o>=o]+"), 0, "MO Ints (oo|oo)");
    dpd_contract422(&I, &T_oo, &Y_oo, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_ij -= <ij|kl> dTau_kl
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo|oo>");
    dpd_contract422(&I, &T_oo, &Y_oo, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y_IJ += (ij|KL) dTau_KL
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[O,O]"),
                  ID("[o,o]"), ID("[O,O]"), 0, "MO Ints (oo|OO)");
    dpd_contract422(&I, &T_OO, &Y_oo, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&I);
    dpd_file2_close(&Y_oo);

    dpd_file2_close(&T_OO);
    dpd_file2_close(&T_oo);
    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);

    //
    // Compute T intermediate
    //

    /*
     * The T intermediates
     * T_ijab = 2.0 * P(ab) Y_ca lambda_ijcb + 2.0 * P(ij) Y_kj lambda_ikab
     */
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    // Temp_IJAB = lambda_IJCB Y_AC
    dpd_file2_init(&Y_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Y <V|V>");
    dpd_contract244(&Y_VV, &Laa, &T, 1, 2, 1, 2.0, 0.0);
    dpd_file2_close(&Y_VV);
    dpd_buf4_close(&Laa);
    // T_IJAB = Temp_IJAB
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "T <OO|VV>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[O,O]"), ID("[V,V]"), "P(Temp) <OO|VV>");
    dpd_buf4_close(&T);
    // T_IJAB -= Temp_IJBA
    dpd_buf4_init(&Taa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "T <OO|VV>");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "P(Temp) <OO|VV>");
    dpd_buf4_add(&Taa, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Temp <OO|VV>");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    // Temp_IJAB = -lambda_KJAB Y_IK
    dpd_file2_init(&Y_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Y <O|O>");
    dpd_contract244(&Y_OO, &Laa, &T, 1, 0, 0, -2.0, 0.0);
    dpd_file2_close(&Y_OO);
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
    dpd_buf4_close(&Taa);

    // The Tab intermediate
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T <Oo|Vv>");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    // T_IjAb = lambda_IjCb Y_AC
    dpd_file2_init(&Y_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Y <V|V>");
    dpd_contract244(&Y_VV, &Lab, &Tab, 1, 2, 1, 2.0, 0.0);
    dpd_file2_close(&Y_VV);
    // T_IjAb += lambda_IjAc Y_bc
    dpd_file2_init(&Y_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Y <v|v>");
    dpd_contract424(&Lab, &Y_vv, &Tab, 3, 1, 0, 2.0, 1.0);
    dpd_file2_close(&Y_vv);
    // T_IjAb -= lambda_KjAb Y_IK
    dpd_file2_init(&Y_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Y <O|O>");
    dpd_contract244(&Y_OO, &Lab, &Tab, 1, 0, 0, -2.0, 1.0);
    dpd_file2_close(&Y_OO);
    // T_IjAb -= lambda_IkAb Y_jk
    dpd_file2_init(&Y_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Y <o|o>");
    dpd_contract424(&Lab, &Y_oo, &Tab, 1, 1, 1, -2.0, 1.0);
    dpd_file2_close(&Y_oo);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Tab);

    // The Tbb intermediate
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    // Temp_ijab = lambda_ijcb Y_ac
    dpd_file2_init(&Y_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Y <v|v>");
    dpd_contract244(&Y_vv, &Lbb, &T, 1, 2, 1, 2.0, 0.0);
    dpd_file2_close(&Y_vv);
    dpd_buf4_close(&Lbb);
    // T_ijab = Temp_ijab
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_copy(&T, PSIF_DCFT_DPD, "T <oo|vv>");
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_sort(&T, PSIF_DCFT_DPD, pqsr, ID("[o,o]"), ID("[v,v]"), "P(Temp) <oo|vv>");
    dpd_buf4_close(&T);

    // T_ijab -= Temp_ijba
    dpd_buf4_init(&Tbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "T <oo|vv>");
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "P(Temp) <oo|vv>");
    dpd_buf4_add(&Tbb, &T, -1.0);
    dpd_buf4_close(&T);
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Temp <oo|vv>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    // Temp_ijab = -lambda_kjab Y_ik
    dpd_file2_init(&Y_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Y <o|o>");
    dpd_contract244(&Y_oo, &Lbb, &T, 1, 0, 0, -2.0, 0.0);
    dpd_file2_close(&Y_oo);
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

    // R_IJAB += T_IJAB
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "T <OO|VV>");
    dpd_buf4_add(&R, &T, 1.0);
    dpd_buf4_close(&T);

    // R_IJAB += C_IJAB
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "C <OO|VV> new");
    dpd_buf4_add(&R, &C, 1.0);
    dpd_buf4_close(&C);

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

    // R_ijab += T_ijab
    dpd_buf4_init(&T, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "T <oo|vv>");
    dpd_buf4_add(&R, &T, 1.0);
    dpd_buf4_close(&T);

    // R_ijab += C_ijab
    dpd_buf4_init(&C, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "C <oo|vv> new");
    dpd_buf4_add(&R, &C, 1.0);
    dpd_buf4_close(&C);

    for(int h = 0; h < nirrep_; ++h)
        nElements += R.params->coltot[h] * R.params->rowtot[h];

    sumSQ += dpd_buf4_dot_self(&R);
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
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "D <OO|VV>");
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "R <OO|VV>");
    dpd_buf4_dirprd(&D, &R);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Z <OO|VV>");
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
                  ID("[o>=o]+"), ID("[v>=v]+"), 0, "D <oo|vv>");
    dpd_buf4_init(&R, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "R <oo|vv>");
    dpd_buf4_dirprd(&D, &R);
    dpd_buf4_close(&D);
    dpd_buf4_init(&Z, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Z <oo|vv>");
    dpd_buf4_add(&Z, &R, 1.0);
    dpd_buf4_close(&R);
    dpd_buf4_close(&Z);

    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::compute_lagrangian_OO()
{

    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 G, I;
    dpdfile2 X, H, pT;

    // X_OO: One-electron contributions

    if (!orbital_optimized_) {
        // X_KI = H_JK pTau_JI
        dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
        dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "H <O|O>");
        dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "pTau <O|O>");
        dpd_file2_mat_init(&X);
        dpd_file2_mat_init(&H);
        dpd_file2_mat_init(&pT);
        dpd_file2_mat_rd(&H);
        dpd_file2_mat_rd(&pT);
        for(int h = 0; h < nirrep_; ++h){
            #pragma omp parallel for
            for(int i = 0 ; i < naoccpi_[h]; ++i){
                for(int k = 0 ; k < naoccpi_[h]; ++k){
                    double value = 0.0;
                    for(int j = 0 ; j < naoccpi_[h]; ++j){
                        value += H.matrix[h][j][k] * (pT.matrix[h][j][i] + (i == j ? 1.0 : 0.0));
                    }
                    X.matrix[h][k][i] = value;
                }
            }
        }
        dpd_file2_mat_wrt(&X);
        dpd_file2_close(&pT);
        dpd_file2_close(&H);
        dpd_file2_close(&X);

        // X_ki = H_jk pTau_ji
        dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
        dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "H <o|o>");
        dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "pTau <o|o>");
        dpd_file2_mat_init(&X);
        dpd_file2_mat_init(&H);
        dpd_file2_mat_init(&pT);
        dpd_file2_mat_rd(&H);
        dpd_file2_mat_rd(&pT);
        for(int h = 0; h < nirrep_; ++h){
            #pragma omp parallel for
            for(int i = 0 ; i < nboccpi_[h]; ++i){
                for(int k = 0 ; k < nboccpi_[h]; ++k){
                    double value = 0.0;
                    for(int j = 0 ; j < nboccpi_[h]; ++j){
                        value += H.matrix[h][j][k] * (pT.matrix[h][j][i] + (i == j ? 1.0 : 0.0));
                    }
                    X.matrix[h][k][i] = value;
                }
            }
        }
        dpd_file2_mat_wrt(&X);
        dpd_file2_close(&pT);
        dpd_file2_close(&H);
        dpd_file2_close(&X);
    }
    else {
        // X_KI = H_JK Tau_JI
        dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
        dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "H <O|O>");
        dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
        dpd_file2_mat_init(&X);
        dpd_file2_mat_init(&H);
        dpd_file2_mat_init(&pT);
        dpd_file2_mat_rd(&H);
        dpd_file2_mat_rd(&pT);
        for(int h = 0; h < nirrep_; ++h){
            #pragma omp parallel for
            for(int i = 0 ; i < naoccpi_[h]; ++i){
                for(int k = 0 ; k < naoccpi_[h]; ++k){
                    double value = 0.0;
                    for(int j = 0 ; j < naoccpi_[h]; ++j){
                        value += H.matrix[h][j][k] * (pT.matrix[h][j][i] + (i == j ? 1.0 : 0.0));
                    }
                    X.matrix[h][k][i] = value;
                }
            }
        }
        dpd_file2_mat_wrt(&X);
        dpd_file2_close(&pT);
        dpd_file2_close(&H);
        dpd_file2_close(&X);

        // X_ki = H_jk Tau_ji
        dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
        dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "H <o|o>");
        dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
        dpd_file2_mat_init(&X);
        dpd_file2_mat_init(&H);
        dpd_file2_mat_init(&pT);
        dpd_file2_mat_rd(&H);
        dpd_file2_mat_rd(&pT);
        for(int h = 0; h < nirrep_; ++h){
            #pragma omp parallel for
            for(int i = 0 ; i < nboccpi_[h]; ++i){
                for(int k = 0 ; k < nboccpi_[h]; ++k){
                    double value = 0.0;
                    for(int j = 0 ; j < nboccpi_[h]; ++j){
                        value += H.matrix[h][j][k] * (pT.matrix[h][j][i] + (i == j ? 1.0 : 0.0));
                    }
                    X.matrix[h][k][i] = value;
                }
            }
        }
        dpd_file2_mat_wrt(&X);
        dpd_file2_close(&pT);
        dpd_file2_close(&H);
        dpd_file2_close(&X);
    }

    // X_OO: Two-electron contributions

    //
    // 2 * <OO||OO> Г_OOOO
    //

    // X_MI += 2 * <MJ||KL> Г_IJKL
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 1, "MO Ints <OO|OO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O>O]-"), ID("[O>O]-"), 0, "Gamma <OO|OO>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_MI += 4 * <Mj|Kl> Г_IjKl
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");

    dpd_contract442(&I, &G, &X, 0, 0, 4.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_mi += 2 * <mj||kl> Г_ijkl
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 1, "MO Ints <oo|oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o>o]-"), ID("[o>o]-"), 0, "Gamma <oo|oo>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_mi += 4 * <Jm|Kl> Г_JiKl
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");

    dpd_contract442(&I, &G, &X, 1, 1, 4.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    //
    // <OO||VV> Г_OOVV
    //

    // X_MI += <JM||BC> Г_JIBC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Gamma <OO|VV>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_MI += <Mj|Bc> Г_IjBc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_mi += <jm||bc> Г_jibc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Gamma <oo|vv>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_mi += <Jm|Bc> Г_JiBc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");

    dpd_contract442(&I, &G, &X, 1, 1, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    //
    // <OV||OV> Г_OVOV
    //

    // X_MI += <JB||MC> Г_JBIC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV> - <OV|VO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_MI += <Jb|Mc> Г_JbIc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_MI -= <Ma|jB> Г_IajB
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "MO Ints <Ov|oV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");

    dpd_contract442(&I, &G, &X, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_mi += <jb||mc> Г_jbic
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov> - <ov|vo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_mi += <jB|mC> Г_jBiC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "MO Ints <oV|oV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_mi -= <Ja|mB> Г_JaiB
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "MO Ints <Ov|oV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");

    dpd_contract442(&I, &G, &X, 2, 2, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    psio_->close(PSIF_DCFT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::compute_lagrangian_VV()
{

    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 G, I;
    dpdfile2 X, H, pT;

    // X_VV: One-electron contributions

    if (!orbital_optimized_) {
        // X_CA = H_CB pTau_BA
        dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
        dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "H <V|V>");
        dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "pTau <V|V>");

        dpd_contract222(&H, &pT, &X, 0, 1, 1.0, 0.0);
        dpd_file2_close(&pT);
        dpd_file2_close(&H);
        dpd_file2_close(&X);

        // X_ca = H_ib pTau_ba
        dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
        dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "H <v|v>");
        dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "pTau <v|v>");

        dpd_contract222(&H, &pT, &X, 0, 1, 1.0, 0.0);
        dpd_file2_close(&pT);
        dpd_file2_close(&H);
        dpd_file2_close(&X);
    }
    else {
        // X_CA = H_CB Tau_BA
        dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
        dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "H <V|V>");
        dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");

        dpd_contract222(&H, &pT, &X, 0, 1, 1.0, 0.0);
        dpd_file2_close(&pT);
        dpd_file2_close(&H);
        dpd_file2_close(&X);

        // X_ca = H_ib Tau_ba
        dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
        dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "H <v|v>");
        dpd_file2_init(&pT, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

        dpd_contract222(&H, &pT, &X, 0, 1, 1.0, 0.0);
        dpd_file2_close(&pT);
        dpd_file2_close(&H);
        dpd_file2_close(&X);
    }

    // X_OV: Two-electron contributions

    //
    // 2 * <VV||VV> Г_VVVV
    //

    // X_EA += 2 * <EB||CD> Г_ABCD
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 1, "MO Ints <VV|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V>V]-"), ID("[V>V]-"), 0, "Gamma <VV|VV>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_EA += 4 * <Eb|Cd> Г_AbCd
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "Gamma <Vv|Vv>");

    dpd_contract442(&I, &G, &X, 0, 0, 4.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ea += 2 * <ib||cd> Г_abcd
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 1, "MO Ints <vv|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v>v]-"), ID("[v>v]-"), 0, "Gamma <vv|vv>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ea += 4 * <eB|cD> Г_AbCd
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "Gamma <Vv|Vv>");

    dpd_contract442(&I, &G, &X, 1, 1, 4.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    //
    // <OO||VV> Г_OOVV
    //

    // X_CA += <BC||JK> Г_BAJK
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V,V]"), ID("[O,O]"), 1, "MO Ints <VV|OO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V>V]-"), ID("[O>O]-"), 0, "Gamma <VV|OO>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_CA += 2 * <Jk|Cb> Г_JkAb
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");

    dpd_contract442(&I, &G, &X, 2, 2, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ca += <bc||jk> Г_bajk
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v,v]"), ID("[o,o]"), 1, "MO Ints <vv|oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v>v]-"), ID("[o>o]-"), 0, "Gamma <vv|oo>");

    dpd_contract442(&I, &G, &X, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ca += 2 * <Jk|Bc> Г_JkBa
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");

    dpd_contract442(&I, &G, &X, 3, 3, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    //
    // <OO||OV> Г_OVOV
    //

    // X_CA += <JB||KC> Г_JBKA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV> - <OV|VO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_CA += <kB|jC> Г_kBjA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "MO Ints <oV|oV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_CA -= <Kb|jC> Г_KbjA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "MO Ints <Ov|oV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");

    dpd_contract442(&I, &G, &X, 3, 3, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ca += <jb||kc> Г_jbka
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov> - <ov|vo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ca += <Kb|Jc> Г_KbJa
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ca -= <kB|Jc> Г_kBJa
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "MO Ints <Ov|oV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");

    dpd_contract442(&I, &G, &X, 1, 1, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    psio_->close(PSIF_DCFT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::compute_ewdm_dc()
{

    dpdfile2 zI_OV, zI_VO, X_OV, X_VO, zI_OO, zI_VV, X_OO, X_VV, z_OV;

    Matrix aW ("Energy-weighted density matrix (Alpha)", nirrep_, nmopi_, nmopi_);
    Matrix bW ("Energy-weighted density matrix (Beta)", nirrep_, nmopi_, nmopi_);

    SharedMatrix a_opdm (new Matrix("MO basis OPDM (Alpha)", nirrep_, nmopi_, nmopi_));
    SharedMatrix b_opdm (new Matrix("MO basis OPDM (Beta)", nirrep_, nmopi_, nmopi_));
    SharedMatrix a_zia (new Matrix("MO basis Orbital Response (Alpha)", nirrep_, nmopi_, nmopi_));
    SharedMatrix b_zia (new Matrix("MO basis Orbital Response (Beta)", nirrep_, nmopi_, nmopi_));

    const int *alpha_corr_to_pitzer = _ints->alpha_corr_to_pitzer();
    int *alpha_pitzer_to_corr = new int[nmo_];
    ::memset(alpha_pitzer_to_corr, '\0', nmo_*sizeof(int));

    for(int n = 0; n < nmo_; ++n) {
        alpha_pitzer_to_corr[alpha_corr_to_pitzer[n]] = n;
    }

    const int *beta_corr_to_pitzer = _ints->beta_corr_to_pitzer();
    int *beta_pitzer_to_corr = new int[nmo_];
    ::memset(beta_pitzer_to_corr, '\0', nmo_*sizeof(int));

    for(int n = 0; n < nmo_; ++n) {
        beta_pitzer_to_corr[beta_corr_to_pitzer[n]] = n;
    }

    // Alpha spin
    dpd_file2_init(&zI_OV, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "zI <O|V>");
    dpd_file2_init(&zI_VO, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "zI <V|O>");
    dpd_file2_init(&X_OV, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_file2_init(&X_VO, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_file2_init(&zI_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "zI <O|O> sym");
    dpd_file2_init(&zI_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "zI <V|V> sym");
    dpd_file2_init(&X_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    dpd_file2_init(&X_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");
    dpd_file2_init(&z_OV, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "z <O|V>");

    dpd_file2_mat_init(&zI_VO);
    dpd_file2_mat_init(&zI_OV);
    dpd_file2_mat_init(&X_OV);
    dpd_file2_mat_init(&X_VO);
    dpd_file2_mat_init(&zI_OO);
    dpd_file2_mat_init(&zI_VV);
    dpd_file2_mat_init(&X_OO);
    dpd_file2_mat_init(&X_VV);
    dpd_file2_mat_init(&z_OV);

    dpd_file2_mat_rd(&zI_VO);
    dpd_file2_mat_rd(&zI_OV);
    dpd_file2_mat_rd(&X_OV);
    dpd_file2_mat_rd(&X_VO);
    dpd_file2_mat_rd(&zI_OO);
    dpd_file2_mat_rd(&zI_VV);
    dpd_file2_mat_rd(&X_OO);
    dpd_file2_mat_rd(&X_VV);
    dpd_file2_mat_rd(&z_OV);

    for(int h = 0; h < nirrep_; ++h){
        // O-V and V-O
        #pragma omp parallel for
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int a = 0 ; a < navirpi_[h]; ++a){
                double value = 0.0;
                for(int j = 0 ; j < naoccpi_[h]; ++j){
                    value -= 0.25 * (zI_VO.matrix[h][a][j] + zI_OV.matrix[h][j][a]) * (aocc_tau_->get(h,i,j) + kappa_mo_a_->get(h,i,j));
                    value -= 0.25 * z_OV.matrix[h][j][a] * moFa_->get(h, j, i);
                }
                for(int b = 0 ; b < navirpi_[h]; ++b){
                    value -= 0.25 * (zI_VO.matrix[h][b][i] + zI_OV.matrix[h][i][b]) * (avir_tau_->get(h,a,b));
                    value -= 0.25 * z_OV.matrix[h][i][b] * moFa_->get(h, b + naoccpi_[h], a + naoccpi_[h]);
                }
                value -= 0.5 * (X_OV.matrix[h][i][a] + X_VO.matrix[h][a][i]);
                aW.set(h, i, a + naoccpi_[h], value);
                aW.set(h, a + naoccpi_[h], i, value);
                a_zia->set(h, i, a + naoccpi_[h], z_OV.matrix[h][i][a]);
            }
        }
        // O-O
        #pragma omp parallel for
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                double value = 0.0;
                for(int k = 0 ; k < naoccpi_[h]; ++k){
                    value -= 0.25 * zI_OO.matrix[h][i][k] * (aocc_tau_->get(h,k,j) + kappa_mo_a_->get(h,k,j));
                    value -= 0.25 * zI_OO.matrix[h][j][k] * (aocc_tau_->get(h,k,i) + kappa_mo_a_->get(h,k,i));
                }
                value -= 0.5 * (X_OO.matrix[h][i][j] + X_OO.matrix[h][j][i]);
                aW.set(h, i, j, value);
                aW.set(h, j, i, value);
                a_opdm->set(h, i, j, (aocc_ptau_->get(h,i,j) + kappa_mo_a_->get(h,i,j)));
                if (i != j) a_opdm->set(h, j, i, (aocc_ptau_->get(h,i,j) + kappa_mo_a_->get(h,i,j)));
            }
        }
        // V-V
        #pragma omp parallel for
        for(int a = 0 ; a < navirpi_[h]; ++a){
            for(int b = 0 ; b <= a; ++b){
                double value = 0.0;
                for(int c = 0 ; c < navirpi_[h]; ++c){
                    value -= 0.25 * zI_VV.matrix[h][a][c] * avir_tau_->get(h,c,b);
                    value -= 0.25 * zI_VV.matrix[h][b][c] * avir_tau_->get(h,c,a);
                }
                value -= 0.5 * (X_VV.matrix[h][a][b] + X_VV.matrix[h][b][a]);
                aW.set(h, a + naoccpi_[h], b + naoccpi_[h], value);
                aW.set(h, b + naoccpi_[h], a + naoccpi_[h], value);
                a_opdm->set(h, a + naoccpi_[h], b + naoccpi_[h], avir_ptau_->get(h, a, b));
                if (a != b) a_opdm->set(h, b + naoccpi_[h], a + naoccpi_[h], avir_ptau_->get(h, a, b));
            }
        }
    }
    dpd_file2_close(&X_VO);
    dpd_file2_close(&X_OV);
    dpd_file2_close(&zI_VO);
    dpd_file2_close(&zI_OV);
    dpd_file2_close(&X_OO);
    dpd_file2_close(&X_VV);
    dpd_file2_close(&zI_OO);
    dpd_file2_close(&zI_VV);
    dpd_file2_close(&z_OV);

    // Beta spin
    dpd_file2_init(&zI_OV, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "zI <o|v>");
    dpd_file2_init(&zI_VO, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "zI <v|o>");
    dpd_file2_init(&X_OV, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_file2_init(&X_VO, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_file2_init(&zI_OO, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "zI <o|o> sym");
    dpd_file2_init(&zI_VV, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "zI <v|v> sym");
    dpd_file2_init(&X_OO, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
    dpd_file2_init(&X_VV, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");
    dpd_file2_init(&z_OV, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "z <o|v>");

    dpd_file2_mat_init(&zI_VO);
    dpd_file2_mat_init(&zI_OV);
    dpd_file2_mat_init(&X_OV);
    dpd_file2_mat_init(&X_VO);
    dpd_file2_mat_init(&zI_OO);
    dpd_file2_mat_init(&zI_VV);
    dpd_file2_mat_init(&X_OO);
    dpd_file2_mat_init(&X_VV);
    dpd_file2_mat_init(&z_OV);

    dpd_file2_mat_rd(&zI_VO);
    dpd_file2_mat_rd(&zI_OV);
    dpd_file2_mat_rd(&X_OV);
    dpd_file2_mat_rd(&X_VO);
    dpd_file2_mat_rd(&zI_OO);
    dpd_file2_mat_rd(&zI_VV);
    dpd_file2_mat_rd(&X_OO);
    dpd_file2_mat_rd(&X_VV);
    dpd_file2_mat_rd(&z_OV);

    for(int h = 0; h < nirrep_; ++h){
        // O-V and V-O
        #pragma omp parallel for
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int a = 0 ; a < nbvirpi_[h]; ++a){
                double value = 0.0;
                for(int j = 0 ; j < nboccpi_[h]; ++j){
                    value -= 0.25 * (zI_VO.matrix[h][a][j] + zI_OV.matrix[h][j][a]) * (bocc_tau_->get(h,i,j) + kappa_mo_b_->get(h,i,j));
                    value -= 0.25 * z_OV.matrix[h][j][a] * moFb_->get(h, j, i);
                }
                for(int b = 0 ; b < nbvirpi_[h]; ++b){
                    value -= 0.25 * (zI_VO.matrix[h][b][i] + zI_OV.matrix[h][i][b]) * (bvir_tau_->get(h,a,b));
                    value -= 0.25 * z_OV.matrix[h][i][b] * moFb_->get(h, b + nboccpi_[h], a + nboccpi_[h]);
                }
                value -= 0.5 * (X_OV.matrix[h][i][a] + X_VO.matrix[h][a][i]);
                b_zia->set(h, i, a + nboccpi_[h], z_OV.matrix[h][i][a]);
                bW.set(h, i, a + nboccpi_[h], value);
                bW.set(h, a + nboccpi_[h], i, value);
                b_zia->set(h, i, a + nboccpi_[h], z_OV.matrix[h][i][a]);
            }
        }
        // O-O
        #pragma omp parallel for
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                double value = 0.0;
                for(int k = 0 ; k < nboccpi_[h]; ++k){
                    value -= 0.25 * zI_OO.matrix[h][i][k] * (bocc_tau_->get(h,k,j) + kappa_mo_b_->get(h,k,j));
                    value -= 0.25 * zI_OO.matrix[h][j][k] * (bocc_tau_->get(h,k,i) + kappa_mo_b_->get(h,k,i));
                }
                value -= 0.5 * (X_OO.matrix[h][i][j] + X_OO.matrix[h][j][i]);
                bW.set(h, i, j, value);
                bW.set(h, j, i, value);
                b_opdm->set(h, i, j, (bocc_ptau_->get(h,i,j) + kappa_mo_b_->get(h,i,j)));
                if (i != j) b_opdm->set(h, j, i, (bocc_ptau_->get(h,i,j) + kappa_mo_b_->get(h,i,j)));
            }
        }
        // V-V
        #pragma omp parallel for
        for(int a = 0 ; a < nbvirpi_[h]; ++a){
            for(int b = 0 ; b <= a; ++b){
                double value = 0.0;
                for(int c = 0 ; c < nbvirpi_[h]; ++c){
                    value -= 0.25 * zI_VV.matrix[h][a][c] * bvir_tau_->get(h,c,b);
                    value -= 0.25 * zI_VV.matrix[h][b][c] * bvir_tau_->get(h,c,a);
                }
                value -= 0.5 * (X_VV.matrix[h][a][b] + X_VV.matrix[h][b][a]);
                bW.set(h, a + nboccpi_[h], b + nboccpi_[h], value);
                bW.set(h, b + nboccpi_[h], a + nboccpi_[h], value);
                b_opdm->set(h, a + nboccpi_[h], b + nboccpi_[h], bvir_ptau_->get(h, a, b));
                if (a != b) b_opdm->set(h, b + nboccpi_[h], a + nboccpi_[h], bvir_ptau_->get(h, a, b));
            }
        }
    }
    dpd_file2_close(&X_VO);
    dpd_file2_close(&X_OV);
    dpd_file2_close(&zI_VO);
    dpd_file2_close(&zI_OV);
    dpd_file2_close(&X_OO);
    dpd_file2_close(&X_VV);
    dpd_file2_close(&zI_OO);
    dpd_file2_close(&zI_VV);
    dpd_file2_close(&z_OV);

    a_opdm->add(a_zia);
    b_opdm->add(b_zia);

    // Scale the energy-weighted density matrix by -2.0 to make it the same form as in the coupled-cluster code
    aW.scale(-2.0);
    bW.scale(-2.0);

    // Reorder the energy-weighted density matrix to the QT order

    double **a_qt = block_matrix(nmo_, nmo_);
    double **b_qt = block_matrix(nmo_, nmo_);

    int offset = 0;
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = alpha_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = alpha_pitzer_to_corr[pitzer_j];
                a_qt[corr_i][corr_j] = aW.get(h, i, j);
            }
        }
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = beta_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = beta_pitzer_to_corr[pitzer_j];
                b_qt[corr_i][corr_j] = bW.get(h, i, j);
            }
        }
        offset += nmopi_[h];
    }

    // Write qt-ordered energy-weighted density matrix to the file
    psio_->open(PSIF_MO_LAG, PSIO_OPEN_OLD);
    psio_->write_entry(PSIF_MO_LAG, "MO-basis Alpha Lagrangian", (char *) a_qt[0], sizeof(double) * nmo_ * nmo_);
    psio_->write_entry(PSIF_MO_LAG, "MO-basis Beta Lagrangian", (char *) b_qt[0], sizeof(double) * nmo_ * nmo_);
    psio_->close(PSIF_MO_LAG, 1);

    // Reorder the one-particle density matrix to the QT order

    offset = 0;
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = alpha_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = alpha_pitzer_to_corr[pitzer_j];
                a_qt[corr_i][corr_j] = a_opdm->get(h, i, j);
            }
        }
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = beta_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = beta_pitzer_to_corr[pitzer_j];
                b_qt[corr_i][corr_j] = b_opdm->get(h, i, j);
            }
        }
        offset += nmopi_[h];
    }

    // Write qt-ordered OPDM to the file
    psio_->open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
    psio_->write_entry(PSIF_MO_OPDM, "MO-basis Alpha OPDM", (char *) a_qt[0], sizeof(double) * nmo_ * nmo_);
    psio_->write_entry(PSIF_MO_OPDM, "MO-basis Beta OPDM", (char *) b_qt[0], sizeof(double) * nmo_ * nmo_);
    psio_->close(PSIF_MO_OPDM, 1);

    int *aocc_qt = new int[nalpha_];
    int *bocc_qt = new int[nbeta_];
    int *avir_qt = new int[navir_];
    int *bvir_qt = new int[nbvir_];

    int aocc_count = 0;
    int bocc_count = 0;
    int avir_count = 0;
    int bvir_count = 0;
    offset = 0;

    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < naoccpi_[h]; ++i) {
            int pitzer = offset + i;
            aocc_qt[aocc_count++] = alpha_pitzer_to_corr[pitzer];
        }
        for (int i = 0; i < nboccpi_[h]; ++i) {
            int pitzer = offset + i;
            bocc_qt[bocc_count++] = beta_pitzer_to_corr[pitzer];
        }
        for (int i = naoccpi_[h]; i < (nmopi_[h] - frzvpi_[h]); ++i) {
            int pitzer = offset + i;
            avir_qt[avir_count++] = alpha_pitzer_to_corr[pitzer];
        }
        for (int i = nboccpi_[h]; i < (nmopi_[h] - frzvpi_[h]); ++i) {
            int pitzer = offset + i;
            bvir_qt[bvir_count++] = beta_pitzer_to_corr[pitzer];
        }
        offset += nmopi_[h];
    }

    dpdbuf4 G;

    struct iwlbuf AA, AB, BB;
    iwl_buf_init(&AA, PSIF_MO_AA_TPDM, 1.0E-15, 0, 0);
    iwl_buf_init(&AB, PSIF_MO_AB_TPDM, 1.0E-15, 0, 0);
    iwl_buf_init(&BB, PSIF_MO_BB_TPDM, 1.0E-15, 0, 0);

    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);

    // Compute the OOOV densities
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O>O]-"), ID("[O,V]"), 0, "Gamma <OO|OV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        #pragma omp parallel for
        for(long int ij = 0; ij < G.params->rowtot[h]; ++ij){
            size_t i = G.params->roworb[h][ij][0];
            int Gi = G.params->psym[i];
            i -= G.params->poff[Gi];
            size_t j = G.params->roworb[h][ij][1];
            int Gj = G.params->qsym[j];
            j -= G.params->qoff[Gj];
            for(size_t ka = 0; ka < G.params->coltot[h]; ++ka){
                size_t k = G.params->colorb[h][ka][0];
                int Gk = G.params->rsym[k];
                k -= G.params->roff[Gk];
                size_t a = G.params->colorb[h][ka][1];
                int Ga = G.params->ssym[a];
                a -= G.params->soff[Ga];
                if(Gi == Gk && Gj == Ga) G.matrix[h][ij][ka] = 0.5 * (kappa_mo_a_->get(Gi, i, k) + aocc_tau_->get(Gi, i, k)) * a_zia->get(Gj, j, a + naoccpi_[Gj]);
                if(Gj == Gk && Gi == Ga) G.matrix[h][ij][ka] -= 0.5 * (kappa_mo_a_->get(Gj, j, k) + aocc_tau_->get(Gj, j, k)) * a_zia->get(Gi, i, a + naoccpi_[Gi]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "Gamma <Oo|Ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        #pragma omp parallel for
        for(long int ij = 0; ij < G.params->rowtot[h]; ++ij){
            size_t i = G.params->roworb[h][ij][0];
            int Gi = G.params->psym[i];
            i -= G.params->poff[Gi];
            size_t j = G.params->roworb[h][ij][1];
            int Gj = G.params->qsym[j];
            j -= G.params->qoff[Gj];
            for(size_t ka = 0; ka < G.params->coltot[h]; ++ka){
                size_t k = G.params->colorb[h][ka][0];
                int Gk = G.params->rsym[k];
                k -= G.params->roff[Gk];
                size_t a = G.params->colorb[h][ka][1];
                int Ga = G.params->ssym[a];
                a -= G.params->soff[Ga];
                if(Gi == Gk && Gj == Ga) G.matrix[h][ij][ka] = 0.5 * (kappa_mo_a_->get(Gi, i, k) + aocc_tau_->get(Gi, i, k)) * b_zia->get(Gj, j, a + nboccpi_[Gj]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,O]"), ID("[o,V]"),
                  ID("[o,O]"), ID("[o,V]"), 0, "Gamma <oO|oV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        #pragma omp parallel for
        for(long int ij = 0; ij < G.params->rowtot[h]; ++ij){
            size_t i = G.params->roworb[h][ij][0];
            int Gi = G.params->psym[i];
            i -= G.params->poff[Gi];
            size_t j = G.params->roworb[h][ij][1];
            int Gj = G.params->qsym[j];
            j -= G.params->qoff[Gj];
            for(size_t ka = 0; ka < G.params->coltot[h]; ++ka){
                size_t k = G.params->colorb[h][ka][0];
                int Gk = G.params->rsym[k];
                k -= G.params->roff[Gk];
                size_t a = G.params->colorb[h][ka][1];
                int Ga = G.params->ssym[a];
                a -= G.params->soff[Ga];
                if(Gi == Gk && Gj == Ga) G.matrix[h][ij][ka] = 0.5 * (kappa_mo_b_->get(Gi, i, k) + bocc_tau_->get(Gi, i, k)) * a_zia->get(Gj, j, a + naoccpi_[Gj]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o>o]-"), ID("[o,v]"), 0, "Gamma <oo|ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        #pragma omp parallel for
        for(long int ij = 0; ij < G.params->rowtot[h]; ++ij){
            size_t i = G.params->roworb[h][ij][0];
            int Gi = G.params->psym[i];
            i -= G.params->poff[Gi];
            size_t j = G.params->roworb[h][ij][1];
            int Gj = G.params->qsym[j];
            j -= G.params->qoff[Gj];
            for(size_t ka = 0; ka < G.params->coltot[h]; ++ka){
                size_t k = G.params->colorb[h][ka][0];
                int Gk = G.params->rsym[k];
                k -= G.params->roff[Gk];
                size_t a = G.params->colorb[h][ka][1];
                int Ga = G.params->ssym[a];
                a -= G.params->soff[Ga];
                if(Gi == Gk && Gj == Ga) G.matrix[h][ij][ka] = 0.5 * (kappa_mo_b_->get(Gi, i, k) + bocc_tau_->get(Gi, i, k)) * b_zia->get(Gj, j, a + nboccpi_[Gj]);
                if(Gj == Gk && Gi == Ga) G.matrix[h][ij][ka] -= 0.5 * (kappa_mo_b_->get(Gj, j, k) + bocc_tau_->get(Gj, j, k)) * b_zia->get(Gi, i, a + nboccpi_[Gi]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // Compute the OVVV densities
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>V]-"), 0, "Gamma <OV|VV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        #pragma omp parallel for
        for(long int ia = 0; ia < G.params->rowtot[h]; ++ia){
            size_t i = G.params->roworb[h][ia][0];
            int Gi = G.params->psym[i];
            i -= G.params->poff[Gi];
            size_t a = G.params->roworb[h][ia][1];
            int Ga = G.params->qsym[a];
            a -= G.params->qoff[Ga];
            for(size_t bc = 0; bc < G.params->coltot[h]; ++bc){
                size_t b = G.params->colorb[h][bc][0];
                int Gb = G.params->rsym[b];
                b -= G.params->roff[Gb];
                size_t c = G.params->colorb[h][bc][1];
                int Gc = G.params->ssym[c];
                c -= G.params->soff[Gc];
                if(Gi == Gb && Ga == Gc) G.matrix[h][ia][bc] = 0.5 * avir_tau_->get(Ga, a, c) * a_zia->get(Gi, i, b + naoccpi_[Gi]);
                if(Gi == Gc && Ga == Gb) G.matrix[h][ia][bc] -= 0.5 * avir_tau_->get(Ga, a, b) * a_zia->get(Gi, i, c + naoccpi_[Gi]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "Gamma <Ov|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        #pragma omp parallel for
        for(long int ia = 0; ia < G.params->rowtot[h]; ++ia){
            size_t i = G.params->roworb[h][ia][0];
            int Gi = G.params->psym[i];
            i -= G.params->poff[Gi];
            size_t a = G.params->roworb[h][ia][1];
            int Ga = G.params->qsym[a];
            a -= G.params->qoff[Ga];
            for(size_t bc = 0; bc < G.params->coltot[h]; ++bc){
                size_t b = G.params->colorb[h][bc][0];
                int Gb = G.params->rsym[b];
                b -= G.params->roff[Gb];
                size_t c = G.params->colorb[h][bc][1];
                int Gc = G.params->ssym[c];
                c -= G.params->soff[Gc];
                if(Gi == Gb && Ga == Gc) G.matrix[h][ia][bc] = 0.5 * bvir_tau_->get(Ga, a, c) * a_zia->get(Gi, i, b + naoccpi_[Gi]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "Gamma <oV|vV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        #pragma omp parallel for
        for(long int ia = 0; ia < G.params->rowtot[h]; ++ia){
            size_t i = G.params->roworb[h][ia][0];
            int Gi = G.params->psym[i];
            i -= G.params->poff[Gi];
            size_t a = G.params->roworb[h][ia][1];
            int Ga = G.params->qsym[a];
            a -= G.params->qoff[Ga];
            for(size_t bc = 0; bc < G.params->coltot[h]; ++bc){
                size_t b = G.params->colorb[h][bc][0];
                int Gb = G.params->rsym[b];
                b -= G.params->roff[Gb];
                size_t c = G.params->colorb[h][bc][1];
                int Gc = G.params->ssym[c];
                c -= G.params->soff[Gc];
                if(Gi == Gb && Ga == Gc) G.matrix[h][ia][bc] = 0.5 * avir_tau_->get(Ga, a, c) * b_zia->get(Gi, i, b + nboccpi_[Gi]);

            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v>v]-"), 0, "Gamma <ov|vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        #pragma omp parallel for
        for(long int ia = 0; ia < G.params->rowtot[h]; ++ia){
            size_t i = G.params->roworb[h][ia][0];
            int Gi = G.params->psym[i];
            i -= G.params->poff[Gi];
            size_t a = G.params->roworb[h][ia][1];
            int Ga = G.params->qsym[a];
            a -= G.params->qoff[Ga];
            for(size_t bc = 0; bc < G.params->coltot[h]; ++bc){
                size_t b = G.params->colorb[h][bc][0];
                int Gb = G.params->rsym[b];
                b -= G.params->roff[Gb];
                size_t c = G.params->colorb[h][bc][1];
                int Gc = G.params->ssym[c];
                c -= G.params->soff[Gc];
                if(Gi == Gb && Ga == Gc) G.matrix[h][ia][bc] = 0.5 * bvir_tau_->get(Ga, a, c) * b_zia->get(Gi, i, b + nboccpi_[Gi]);
                if(Gi == Gc && Ga == Gb) G.matrix[h][ia][bc] -= 0.5 * bvir_tau_->get(Ga, a, b) * b_zia->get(Gi, i, c + nboccpi_[Gi]);
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // Dump the density to IWL
    // VVVV

    // Note (TODO): So far the contraction of TEI derivatives with Gamma is not restricted.
    // We need to restrict it by scaling Gamma and switching bk_pack in dpd_buf4_dump to 1
    // If restriction is used it seems that one needs to do a seperate dpd_buf4_sort to the chemists' notation
    // before calling dpd_buf4_dump - check that with swap23 = 0 (AYS)
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
              ID("[V>V]-"), ID("[V>V]-"), 0, "Gamma <VV|VV>");
    dpd_buf4_dump(&G, &AA, avir_qt, avir_qt, avir_qt, avir_qt, 0, 1);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
              ID("[V,v]"), ID("[V,v]"), 0, "Gamma <Vv|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ab = 0; ab < G.params->rowtot[h]; ++ab){
            int a = G.params->roworb[h][ab][0];
            int b = G.params->roworb[h][ab][1];
            int A = avir_qt[a];
            int B = bvir_qt[b];
            for(size_t cd = 0; cd < G.params->coltot[h]; ++cd){
                int c = G.params->colorb[h][cd][0];
                int d = G.params->colorb[h][cd][1];
                int C = avir_qt[c];
                int D = bvir_qt[d];
                double value = 4.0 * G.matrix[h][ab][cd];
                iwl_buf_wrt_val(&AB, A, C, B, D, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
              ID("[v>v]-"), ID("[v>v]-"), 0, "Gamma <vv|vv>");
    dpd_buf4_dump(&G, &BB, bvir_qt, bvir_qt, bvir_qt, bvir_qt, 0, 1);
    dpd_buf4_close(&G);

    // OOOO
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
              ID("[O>O]-"), ID("[O>O]-"), 0, "Gamma <OO|OO>");
    dpd_buf4_dump(&G, &AA, aocc_qt, aocc_qt, aocc_qt, aocc_qt, 0, 1);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
              ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");

    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            int I = aocc_qt[i];
            int J = bocc_qt[j];
            for(size_t kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
                int K = aocc_qt[k];
                int L = bocc_qt[l];
                double value = 4.0 * G.matrix[h][ij][kl];
                iwl_buf_wrt_val(&AB, I, K, J, L, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
              ID("[o>o]-"), ID("[o>o]-"), 0, "Gamma <oo|oo>");
    dpd_buf4_dump(&G, &BB, bocc_qt, bocc_qt, bocc_qt, bocc_qt, 0, 1);
    dpd_buf4_close(&G);

    // OOVV
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
              ID("[O>O]-"), ID("[V>V]-"), 0, "Gamma <OO|VV>");
    dpd_buf4_dump(&G, &AA, aocc_qt, aocc_qt, avir_qt, avir_qt, 0, 1);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
              ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            int I = aocc_qt[i];
            int J = bocc_qt[j];
            for(size_t ab = 0; ab < G.params->coltot[h]; ++ab){
                int a = G.params->colorb[h][ab][0];
                int b = G.params->colorb[h][ab][1];
                int A = avir_qt[a];
                int B = bvir_qt[b];
                double value = 4.0 * G.matrix[h][ij][ab];
                iwl_buf_wrt_val(&AB, I, A, J, B, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
              ID("[o>o]-"), ID("[v>v]-"), 0, "Gamma <oo|vv>");
    dpd_buf4_dump(&G, &BB, bocc_qt, bocc_qt, bvir_qt, bvir_qt, 0, 1);
    dpd_buf4_close(&G);

    // OVOV
    // Г<OV|OV> must be antisymmetrized before contracting it with TEI derivatives:
    // Г<OV|OV> contribution
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
              ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = aocc_qt[i];
            int A = avir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = aocc_qt[j];
                int B = avir_qt[b];
                double value = 0.5 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AA, I, J, A, B, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // -Г<ov|vo> contribution:
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
              ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = aocc_qt[i];
            int A = avir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = aocc_qt[j];
                int B = avir_qt[b];
                double value = -0.5 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AA, I, B, A, J, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
              ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
    dpd_buf4_dump(&G, &AB, aocc_qt, bvir_qt, aocc_qt, bvir_qt, 0, 1);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = bocc_qt[i];
            int A = avir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = bocc_qt[j];
                int B = avir_qt[b];
                double value = G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AB, A, B, I, J, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // Resort Г<Ia|jB> -> (-1.0) * Г(IB|ja) and dump it
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = aocc_qt[i];
            int A = bvir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = bocc_qt[j];
                int B = avir_qt[b];
                double value = (-2.0) * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AB, I, B, J, A, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // Г<ov|ov> must be antisymmetrized before contracting it with TEI derivatives:
    // Г<ov|ov> contribution
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
              ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = bocc_qt[i];
            int A = bvir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = bocc_qt[j];
                int B = bvir_qt[b];
                double value = 0.5 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&BB, I, J, A, B, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // -Г<ov|vo> contribution:
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = bocc_qt[i];
            int A = bvir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = bocc_qt[j];
                int B = bvir_qt[b];
                double value = -0.5 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&BB, I, B, A, J, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // OOOV
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,V]"),
              ID("[O>O]-"), ID("[O,V]"), 0, "Gamma <OO|OV>");
    dpd_buf4_dump(&G, &AA, aocc_qt, aocc_qt, aocc_qt, avir_qt, 0, 1);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,v]"),
              ID("[O,o]"), ID("[O,v]"), 0, "Gamma <Oo|Ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            int I = aocc_qt[i];
            int J = bocc_qt[j];
            for(size_t ka = 0; ka < G.params->coltot[h]; ++ka){
                int k = G.params->colorb[h][ka][0];
                int a = G.params->colorb[h][ka][1];
                int K = aocc_qt[k];
                int A = bvir_qt[a];
                double value = G.matrix[h][ij][ka];
                iwl_buf_wrt_val(&AB, I, K, J, A, value, 0, (FILE *) NULL, 0);
                iwl_buf_wrt_val(&AB, I, K, A, J, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,O]"), ID("[o,V]"),
              ID("[o,O]"), ID("[o,V]"), 0, "Gamma <oO|oV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            int I = bocc_qt[i];
            int J = aocc_qt[j];
            for(size_t ka = 0; ka < G.params->coltot[h]; ++ka){
                int k = G.params->colorb[h][ka][0];
                int a = G.params->colorb[h][ka][1];
                int K = bocc_qt[k];
                int A = avir_qt[a];
                double value = G.matrix[h][ij][ka];
                iwl_buf_wrt_val(&AB, A, J, K, I, value, 0, (FILE *) NULL, 0);
                iwl_buf_wrt_val(&AB, J, A, K, I, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,v]"),
              ID("[o>o]-"), ID("[o,v]"), 0, "Gamma <oo|ov>");
    dpd_buf4_dump(&G, &BB, bocc_qt, bocc_qt, bocc_qt, bvir_qt, 0, 1);
    dpd_buf4_close(&G);

    // OVVV
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[V,V]"),
              ID("[O,V]"), ID("[V>V]-"), 0, "Gamma <OV|VV>");
    dpd_buf4_dump(&G, &AA, aocc_qt, avir_qt, avir_qt, avir_qt, 0, 1);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[V,v]"),
              ID("[O,v]"), ID("[V,v]"), 0, "Gamma <Ov|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = aocc_qt[i];
            int A = bvir_qt[a];
            for(size_t bc = 0; bc < G.params->coltot[h]; ++bc){
                int b = G.params->colorb[h][bc][0];
                int c = G.params->colorb[h][bc][1];
                int B = avir_qt[b];
                int C = bvir_qt[c];
                double value = G.matrix[h][ia][bc];
                iwl_buf_wrt_val(&AB, I, B, A, C, value, 0, (FILE *) NULL, 0);
                iwl_buf_wrt_val(&AB, B, I, A, C, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[v,V]"),
              ID("[o,V]"), ID("[v,V]"), 0, "Gamma <oV|vV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = bocc_qt[i];
            int A = avir_qt[a];
            for(size_t bc = 0; bc < G.params->coltot[h]; ++bc){
                int b = G.params->colorb[h][bc][0];
                int c = G.params->colorb[h][bc][1];
                int B = bvir_qt[b];
                int C = avir_qt[c];
                double value = G.matrix[h][ia][bc];
                iwl_buf_wrt_val(&AB, C, A, B, I, value, 0, (FILE *) NULL, 0);
                iwl_buf_wrt_val(&AB, C, A, I, B, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[v,v]"),
              ID("[o,v]"), ID("[v>v]-"), 0, "Gamma <ov|vv>");
    dpd_buf4_dump(&G, &BB, bocc_qt, bvir_qt, bvir_qt, bvir_qt, 0, 1);
    dpd_buf4_close(&G);

    psio_->close(PSIF_DCFT_DENSITY, 1);

    iwl_buf_flush(&AA, 1);
    iwl_buf_flush(&AB, 1);
    iwl_buf_flush(&BB, 1);
    iwl_buf_close(&AA, 1);
    iwl_buf_close(&AB, 1);
    iwl_buf_close(&BB, 1);

}

void
DCFTSolver::compute_ewdm_odc()
{

    dpdfile2 X_OV, X_VO, X_OO, X_VV;

    Matrix aW ("Energy-weighted density matrix (Alpha)", nirrep_, nmopi_, nmopi_);
    Matrix bW ("Energy-weighted density matrix (Beta)", nirrep_, nmopi_, nmopi_);

    SharedMatrix a_opdm (new Matrix("MO basis OPDM (Alpha)", nirrep_, nmopi_, nmopi_));
    SharedMatrix b_opdm (new Matrix("MO basis OPDM (Beta)", nirrep_, nmopi_, nmopi_));

    const int *alpha_corr_to_pitzer = _ints->alpha_corr_to_pitzer();
    int *alpha_pitzer_to_corr = new int[nmo_];
    ::memset(alpha_pitzer_to_corr, '\0', nmo_*sizeof(int));

    for(int n = 0; n < nmo_; ++n) {
        alpha_pitzer_to_corr[alpha_corr_to_pitzer[n]] = n;
    }

    const int *beta_corr_to_pitzer = _ints->beta_corr_to_pitzer();
    int *beta_pitzer_to_corr = new int[nmo_];
    ::memset(beta_pitzer_to_corr, '\0', nmo_*sizeof(int));

    for(int n = 0; n < nmo_; ++n) {
        beta_pitzer_to_corr[beta_corr_to_pitzer[n]] = n;
    }

    // Alpha spin
    dpd_file2_init(&X_OV, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_file2_init(&X_VO, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_file2_init(&X_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "X <O|O>");
    dpd_file2_init(&X_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "X <V|V>");

    dpd_file2_mat_init(&X_OV);
    dpd_file2_mat_init(&X_VO);
    dpd_file2_mat_init(&X_OO);
    dpd_file2_mat_init(&X_VV);

    dpd_file2_mat_rd(&X_OV);
    dpd_file2_mat_rd(&X_VO);
    dpd_file2_mat_rd(&X_OO);
    dpd_file2_mat_rd(&X_VV);

    for(int h = 0; h < nirrep_; ++h){
        // O-V and V-O
        #pragma omp parallel for
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int a = 0 ; a < navirpi_[h]; ++a){
                double value = -0.5 * (X_OV.matrix[h][i][a] + X_VO.matrix[h][a][i]);
                aW.set(h, i, a + naoccpi_[h], value);
                aW.set(h, a + naoccpi_[h], i, value);
            }
        }
        // O-O
        #pragma omp parallel for
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                double value = -0.5 * (X_OO.matrix[h][i][j] + X_OO.matrix[h][j][i]);
                aW.set(h, i, j, value);
                aW.set(h, j, i, value);
                a_opdm->set(h, i, j, (aocc_tau_->get(h,i,j) + kappa_mo_a_->get(h,i,j)));
                if (i != j) a_opdm->set(h, j, i, (aocc_tau_->get(h,i,j) + kappa_mo_a_->get(h,i,j)));
            }
        }
        // V-V
        #pragma omp parallel for
        for(int a = 0 ; a < navirpi_[h]; ++a){
            for(int b = 0 ; b <= a; ++b){
                double value = -0.5 * (X_VV.matrix[h][a][b] + X_VV.matrix[h][b][a]);
                aW.set(h, a + naoccpi_[h], b + naoccpi_[h], value);
                aW.set(h, b + naoccpi_[h], a + naoccpi_[h], value);
                a_opdm->set(h, a + naoccpi_[h], b + naoccpi_[h], avir_tau_->get(h, a, b));
                if (a != b) a_opdm->set(h, b + naoccpi_[h], a + naoccpi_[h], avir_tau_->get(h, a, b));
            }
        }
    }
    dpd_file2_close(&X_VO);
    dpd_file2_close(&X_OV);
    dpd_file2_close(&X_OO);
    dpd_file2_close(&X_VV);

    // Beta spin
    dpd_file2_init(&X_OV, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_file2_init(&X_VO, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_file2_init(&X_OO, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "X <o|o>");
    dpd_file2_init(&X_VV, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "X <v|v>");

    dpd_file2_mat_init(&X_OV);
    dpd_file2_mat_init(&X_VO);
    dpd_file2_mat_init(&X_OO);
    dpd_file2_mat_init(&X_VV);

    dpd_file2_mat_rd(&X_OV);
    dpd_file2_mat_rd(&X_VO);
    dpd_file2_mat_rd(&X_OO);
    dpd_file2_mat_rd(&X_VV);

    for(int h = 0; h < nirrep_; ++h){
        // O-V and V-O
        #pragma omp parallel for
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int a = 0 ; a < nbvirpi_[h]; ++a){
                double value = -0.5 * (X_OV.matrix[h][i][a] + X_VO.matrix[h][a][i]);
                bW.set(h, i, a + nboccpi_[h], value);
                bW.set(h, a + nboccpi_[h], i, value);
            }
        }
        // O-O
        #pragma omp parallel for
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int j = 0 ; j <= i; ++j){
                double value = -0.5 * (X_OO.matrix[h][i][j] + X_OO.matrix[h][j][i]);
                bW.set(h, i, j, value);
                bW.set(h, j, i, value);
                b_opdm->set(h, i, j, (bocc_tau_->get(h,i,j) + kappa_mo_b_->get(h,i,j)));
                if (i != j) b_opdm->set(h, j, i, (bocc_tau_->get(h,i,j) + kappa_mo_b_->get(h,i,j)));
            }
        }
        // V-V
        #pragma omp parallel for
        for(int a = 0 ; a < nbvirpi_[h]; ++a){
            for(int b = 0 ; b <= a; ++b){
                double value = -0.5 * (X_VV.matrix[h][a][b] + X_VV.matrix[h][b][a]);
                bW.set(h, a + nboccpi_[h], b + nboccpi_[h], value);
                bW.set(h, b + nboccpi_[h], a + nboccpi_[h], value);
                b_opdm->set(h, a + nboccpi_[h], b + nboccpi_[h], bvir_tau_->get(h, a, b));
                if (a != b) b_opdm->set(h, b + nboccpi_[h], a + nboccpi_[h], bvir_tau_->get(h, a, b));
            }
        }
    }
    dpd_file2_close(&X_VO);
    dpd_file2_close(&X_OV);
    dpd_file2_close(&X_OO);
    dpd_file2_close(&X_VV);

    // Scale the energy-weighted density matrix by -2.0 to make it the same form as in the coupled-cluster code
    aW.scale(-2.0);
    bW.scale(-2.0);

    // Reorder the energy-weighted density matrix to the QT order
    double **a_qt = block_matrix(nmo_, nmo_);
    double **b_qt = block_matrix(nmo_, nmo_);

    int offset = 0;
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = alpha_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = alpha_pitzer_to_corr[pitzer_j];
                a_qt[corr_i][corr_j] = aW.get(h, i, j);
            }
        }
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = beta_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = beta_pitzer_to_corr[pitzer_j];
                b_qt[corr_i][corr_j] = bW.get(h, i, j);
            }
        }
        offset += nmopi_[h];
    }

    // Write qt-ordered energy-weighted density matrix to the file
    psio_->open(PSIF_MO_LAG, PSIO_OPEN_OLD);
    psio_->write_entry(PSIF_MO_LAG, "MO-basis Alpha Lagrangian", (char *) a_qt[0], sizeof(double) * nmo_ * nmo_);
    psio_->write_entry(PSIF_MO_LAG, "MO-basis Beta Lagrangian", (char *) b_qt[0], sizeof(double) * nmo_ * nmo_);
    psio_->close(PSIF_MO_LAG, 1);

    // Reorder the one-particle density matrix to the QT order
    offset = 0;
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = alpha_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = alpha_pitzer_to_corr[pitzer_j];
                a_qt[corr_i][corr_j] = a_opdm->get(h, i, j);
            }
        }
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = beta_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = beta_pitzer_to_corr[pitzer_j];
                b_qt[corr_i][corr_j] = b_opdm->get(h, i, j);
            }
        }
        offset += nmopi_[h];
    }

    // Write qt-ordered OPDM to the file
    psio_->open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
    psio_->write_entry(PSIF_MO_OPDM, "MO-basis Alpha OPDM", (char *) a_qt[0], sizeof(double) * nmo_ * nmo_);
    psio_->write_entry(PSIF_MO_OPDM, "MO-basis Beta OPDM", (char *) b_qt[0], sizeof(double) * nmo_ * nmo_);
    psio_->close(PSIF_MO_OPDM, 1);

    int *aocc_qt = new int[nalpha_];
    int *bocc_qt = new int[nbeta_];
    int *avir_qt = new int[navir_];
    int *bvir_qt = new int[nbvir_];

    int aocc_count = 0;
    int bocc_count = 0;
    int avir_count = 0;
    int bvir_count = 0;
    offset = 0;

    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < naoccpi_[h]; ++i) {
            int pitzer = offset + i;
            aocc_qt[aocc_count++] = alpha_pitzer_to_corr[pitzer];
        }
        for (int i = 0; i < nboccpi_[h]; ++i) {
            int pitzer = offset + i;
            bocc_qt[bocc_count++] = beta_pitzer_to_corr[pitzer];
        }
        for (int i = naoccpi_[h]; i < (nmopi_[h] - frzvpi_[h]); ++i) {
            int pitzer = offset + i;
            avir_qt[avir_count++] = alpha_pitzer_to_corr[pitzer];
        }
        for (int i = nboccpi_[h]; i < (nmopi_[h] - frzvpi_[h]); ++i) {
            int pitzer = offset + i;
            bvir_qt[bvir_count++] = beta_pitzer_to_corr[pitzer];
        }
        offset += nmopi_[h];
    }

    dpdbuf4 G;

    struct iwlbuf AA, AB, BB;
    iwl_buf_init(&AA, PSIF_MO_AA_TPDM, 1.0E-15, 0, 0);
    iwl_buf_init(&AB, PSIF_MO_AB_TPDM, 1.0E-15, 0, 0);
    iwl_buf_init(&BB, PSIF_MO_BB_TPDM, 1.0E-15, 0, 0);

    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);

    // Dump the density to IWL
    // VVVV

    // Note (TODO): So far the contraction of TEI derivatives with Gamma is not restricted.
    // We need to restrict it by scaling Gamma and switching bk_pack in dpd_buf4_dump to 1
    // If restriction is used it seems that one needs to do a seperate dpd_buf4_sort to the chemists' notation
    // before calling dpd_buf4_dump - check that with swap23 = 0 (AYS)
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
              ID("[V>V]-"), ID("[V>V]-"), 0, "Gamma <VV|VV>");
    dpd_buf4_dump(&G, &AA, avir_qt, avir_qt, avir_qt, avir_qt, 0, 1);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
              ID("[V,v]"), ID("[V,v]"), 0, "Gamma <Vv|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ab = 0; ab < G.params->rowtot[h]; ++ab){
            int a = G.params->roworb[h][ab][0];
            int b = G.params->roworb[h][ab][1];
            int A = avir_qt[a];
            int B = bvir_qt[b];
            for(size_t cd = 0; cd < G.params->coltot[h]; ++cd){
                int c = G.params->colorb[h][cd][0];
                int d = G.params->colorb[h][cd][1];
                int C = avir_qt[c];
                int D = bvir_qt[d];
                double value = 4.0 * G.matrix[h][ab][cd];
                iwl_buf_wrt_val(&AB, A, C, B, D, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
              ID("[v>v]-"), ID("[v>v]-"), 0, "Gamma <vv|vv>");
    dpd_buf4_dump(&G, &BB, bvir_qt, bvir_qt, bvir_qt, bvir_qt, 0, 1);
    dpd_buf4_close(&G);

    // OOOO
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
              ID("[O>O]-"), ID("[O>O]-"), 0, "Gamma <OO|OO>");
    dpd_buf4_dump(&G, &AA, aocc_qt, aocc_qt, aocc_qt, aocc_qt, 0, 1);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
              ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");

    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            int I = aocc_qt[i];
            int J = bocc_qt[j];
            for(size_t kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
                int K = aocc_qt[k];
                int L = bocc_qt[l];
                double value = 4.0 * G.matrix[h][ij][kl];
                iwl_buf_wrt_val(&AB, I, K, J, L, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
              ID("[o>o]-"), ID("[o>o]-"), 0, "Gamma <oo|oo>");
    dpd_buf4_dump(&G, &BB, bocc_qt, bocc_qt, bocc_qt, bocc_qt, 0, 1);
    dpd_buf4_close(&G);

    // OOVV
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
              ID("[O>O]-"), ID("[V>V]-"), 0, "Gamma <OO|VV>");
    dpd_buf4_dump(&G, &AA, aocc_qt, aocc_qt, avir_qt, avir_qt, 0, 1);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
              ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            int I = aocc_qt[i];
            int J = bocc_qt[j];
            for(size_t ab = 0; ab < G.params->coltot[h]; ++ab){
                int a = G.params->colorb[h][ab][0];
                int b = G.params->colorb[h][ab][1];
                int A = avir_qt[a];
                int B = bvir_qt[b];
                double value = 4.0 * G.matrix[h][ij][ab];
                iwl_buf_wrt_val(&AB, I, A, J, B, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
              ID("[o>o]-"), ID("[v>v]-"), 0, "Gamma <oo|vv>");
    dpd_buf4_dump(&G, &BB, bocc_qt, bocc_qt, bvir_qt, bvir_qt, 0, 1);
    dpd_buf4_close(&G);

    // OVOV
    // Г<OV|OV> must be antisymmetrized before contracting it with TEI derivatives:
    // Г<OV|OV> contribution
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
              ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = aocc_qt[i];
            int A = avir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = aocc_qt[j];
                int B = avir_qt[b];
                double value = 0.5 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AA, I, J, A, B, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // -Г<ov|vo> contribution:
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
              ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = aocc_qt[i];
            int A = avir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = aocc_qt[j];
                int B = avir_qt[b];
                double value = -0.5 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AA, I, B, A, J, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
              ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
    dpd_buf4_dump(&G, &AB, aocc_qt, bvir_qt, aocc_qt, bvir_qt, 0, 1);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = bocc_qt[i];
            int A = avir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = bocc_qt[j];
                int B = avir_qt[b];
                double value = G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AB, A, B, I, J, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // Resort Г<Ia|jB> -> (-1.0) * Г(IB|ja) and dump it
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = aocc_qt[i];
            int A = bvir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = bocc_qt[j];
                int B = avir_qt[b];
                double value = (-2.0) * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AB, I, B, J, A, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // Г<ov|ov> must be antisymmetrized before contracting it with TEI derivatives:
    // Г<ov|ov> contribution
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
              ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = bocc_qt[i];
            int A = bvir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = bocc_qt[j];
                int B = bvir_qt[b];
                double value = 0.5 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&BB, I, J, A, B, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // -Г<ov|vo> contribution:
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = bocc_qt[i];
            int A = bvir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = bocc_qt[j];
                int B = bvir_qt[b];
                double value = -0.5 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&BB, I, B, A, J, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    psio_->close(PSIF_DCFT_DENSITY, 1);

    iwl_buf_flush(&AA, 1);
    iwl_buf_flush(&AB, 1);
    iwl_buf_flush(&BB, 1);
    iwl_buf_close(&AA, 1);
    iwl_buf_close(&AB, 1);
    iwl_buf_close(&BB, 1);

}

}} //End namespaces
