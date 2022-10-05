/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libdiis/diismanager.h"
#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "occwave.h"
#include "defines.h"

#include <cmath>

using namespace psi;

namespace psi {
namespace occwave {

void OCCWave::occ_iterations() {
    idp();

    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    if (wfn_type_ == "OMP2")
        outfile->Printf(" ================ Performing OMP2 iterations... =============================== \n");
    else if (wfn_type_ == "OMP3")
        outfile->Printf(" ================ Performing OMP3 iterations... =============================== \n");
    else if (wfn_type_ == "OCEPA")
        outfile->Printf(" ================ Performing OLCCD iterations... ============================== \n");
    else if (wfn_type_ == "OMP2.5")
        outfile->Printf(" ================ Performing OMP2.5 iterations... ============================= \n");
    else if (wfn_type_ == "OREMP")
        outfile->Printf(" ================ Performing OREMP iterations... ============================= \n");
    outfile->Printf(" ============================================================================== \n");
    if (wfn_type_ == "OMP2")
        outfile->Printf("\t            Minimizing MP2-L Functional \n");
    else if (wfn_type_ == "OMP3")
        outfile->Printf("\t            Minimizing MP3-L Functional \n");
    else if (wfn_type_ == "OCEPA")
        outfile->Printf("\t            Minimizing LCCD-L Functional \n");
    else if (wfn_type_ == "OMP2.5")
        outfile->Printf("\t            Minimizing MP2.5-L Functional \n");
    else if (wfn_type_ == "OREMP")
        outfile->Printf("\t            Minimizing REMP-L Functional \n");
    outfile->Printf("\t            --------------------------- \n");
    outfile->Printf(" Iter       E_total           DE           RMS MO Grad      MAX MO Grad      RMS T2    \n");
    outfile->Printf(" ----    ---------------    ----------     -----------      -----------     ---------- \n");

    /********************************************************************************************/
    /************************** NR iterations **************************************************/
    /********************************************************************************************/
    itr_occ = 0;
    conver = 1;  // Assuming that the MOs will be optimized.
    mo_optimized = 0;
    itr_diis = 0;

    DIISManager orbital_diis;
    // If diis?
    // if (nooA + nooB != 1) {
    if (do_diis_ == 1) {
        orbital_diis = DIISManager(maxdiis_, "Orbital Optimized DIIS", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::OnDisk);
        std::string tensor_name = (wfn_type_ == "OCEPA" || wfn_type_ == "OREMP") ? "T2" : (wfn_type_ == "OMP2" ? "T" : "T2_1");
        if (reference_ == "RESTRICTED") {
            dpdbuf4 T;
            std::string temp1 = tensor_name + " <OO|VV>";
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0, temp1.c_str());
            if (wfn_type_ == "OMP2.5" || wfn_type_ == "OMP3") {
                orbital_diis.set_error_vector_size(kappa_bar_[SpinType::Alpha].get(), &T, &T);
                orbital_diis.set_vector_size(kappa_bar_[SpinType::Alpha].get(), &T, &T);
            } else {
                orbital_diis.set_error_vector_size(kappa_bar_[SpinType::Alpha].get(), &T);
                orbital_diis.set_vector_size(kappa_bar_[SpinType::Alpha].get(), &T);
            }
            global_dpd_->buf4_close(&T);
        } else if (reference_ == "UNRESTRICTED") {
            dpdbuf4 Taa, Tab, Tbb;
            // You're reading the below code right. The same-spin T amplitudes are stored on disk and in-memory without antisymmetry packing.
            // I don't understand either. Don't you love legacy code?
            std::string temp1 = tensor_name + " <OO|VV>";
            global_dpd_->buf4_init(&Taa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
               temp1.c_str());
            temp1 = tensor_name + " <Oo|Vv>";
            global_dpd_->buf4_init(&Tab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               temp1.c_str());
            temp1 = tensor_name + " <oo|vv>";
            global_dpd_->buf4_init(&Tbb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               temp1.c_str());
            if (wfn_type_ == "OMP2.5" || wfn_type_ == "OMP3") {
                orbital_diis.set_error_vector_size(kappa_bar_[SpinType::Alpha].get(),
                                                   kappa_bar_[SpinType::Beta].get(),
                                                   &Taa, &Tab, &Tbb, &Taa, &Tab, &Tbb);
                orbital_diis.set_vector_size(kappa_bar_[SpinType::Alpha].get(),
                                             kappa_bar_[SpinType::Beta].get(),
                                             &Taa, &Tab, &Tbb, &Taa, &Tab, &Tbb);
            } else {
                orbital_diis.set_error_vector_size(kappa_bar_[SpinType::Alpha].get(),
                                                   kappa_bar_[SpinType::Beta].get(),
                                                   &Taa, &Tab, &Tbb);
                orbital_diis.set_vector_size(kappa_bar_[SpinType::Alpha].get(),
                                             kappa_bar_[SpinType::Beta].get(),
                                             &Taa, &Tab, &Tbb);
            }
            global_dpd_->buf4_close(&Taa);
            global_dpd_->buf4_close(&Tab);
            global_dpd_->buf4_close(&Tbb);
        }
    }
    //}

    // Set up the orb-resp algorithm
    if (opt_method == "ORB_RESP") {
        // Lineq
        if (orb_resp_solver_ == "LINEQ") {
            if (reference_ == "RESTRICTED")
                Aorb = new Array2d("MO Hessian Matrix", nidpA, nidpA);
            else if (reference_ == "UNRESTRICTED") {
                nidp_tot = nidpA + nidpB;
                kappa = new Array1d("Total orb rot params vector of current step", nidp_tot);
            }
        }

        // PCG
        else if (orb_resp_solver_ == "PCG") {
            r_pcgA = new Array1d("Alpha PCG r vector", nidpA);
            S_pcgA = new Array1d("Alpha PCG search direction vector", nidpA);
            D_pcgA = new Array1d("Alpha PCG conjugate direction vector", nidpA);
            r_pcg_newA = new Array1d("Alpha New PCG r vector", nidpA);
            sigma_pcgA = new Array1d("Alpha PCG sigma vector", nidpA);
            Minv_pcgA = new Array1d("Alpha PCG inverse of M matrix", nidpA);
            r_pcgA->zero();
            S_pcgA->zero();
            sigma_pcgA->zero();
            D_pcgA->zero();
            Minv_pcgA->zero();

            if (pcg_beta_type_ == "POLAK_RIBIERE") {
                dr_pcgA = new Array1d("Alpha PCG dr vector", nidpA);
                r_pcgA->zero();
            }

            if (reference_ == "UNRESTRICTED") {
                r_pcgB = new Array1d("Beta PCG r vector", nidpB);
                S_pcgB = new Array1d("Beta PCG search direction vector", nidpB);
                D_pcgB = new Array1d("Beta PCG conjugate direction vector", nidpB);
                r_pcg_newB = new Array1d("Beta New PCG r vector", nidpB);
                sigma_pcgB = new Array1d("Beta PCG sigma vector", nidpB);
                Minv_pcgB = new Array1d("Beta PCG inverse of M matrix", nidpB);
                r_pcgB->zero();
                S_pcgB->zero();
                sigma_pcgB->zero();
                D_pcgB->zero();
                Minv_pcgB->zero();
                if (pcg_beta_type_ == "POLAK_RIBIERE") {
                    dr_pcgB = new Array1d("Alpha PCG dr vector", nidpB);
                    r_pcgB->zero();
                }
            }
        }  // pcg if
    }      // orb_resp if

    // Construct initial orbital gradient. Assume the necessary intermediates were
    // constructed earlier. (From the computation with the starting orbitals.)
    response_pdms();
    gfock();
    mograd();
    compute_orbital_step();

    /********************************************************************************************/
    /************************** Head of the Loop ************************************************/
    /********************************************************************************************/
    do {
        itr_occ++;

        /********************************************************************************************/
        /************************** update mo coefficients ******************************************/
        /********************************************************************************************/
        timer_on("update_mo");
        update_mo_spincase(SpinType::Alpha);
        if (reference_ == "UNRESTRICTED") update_mo_spincase(SpinType::Beta);
        timer_off("update_mo");

        /********************************************************************************************/
        /************************** Transform TEI from SO to MO space *******************************/
        /********************************************************************************************/
        timer_on("trans_ints");
        if (reference_ == "RESTRICTED")
            trans_ints_rhf();
        else if (reference_ == "UNRESTRICTED")
            trans_ints_uhf();
        timer_off("trans_ints");

        /********************************************************************************************/
        /************************** One-particle and two-particle density matrices ******************/
        /********************************************************************************************/
        timer_on("Response PDMs");
        response_pdms();
        timer_off("Response PDMs");

        /********************************************************************************************/
        /************************** Asymmetric Generalized-Fock matrix ******************************/
        /********************************************************************************************/
        timer_on("Generalized-Fock");
        gfock();
        timer_off("Generalized-Fock");

        /********************************************************************************************/
        /************************** Compute Lagrangian Energy ***************************************/
        /********************************************************************************************/
        if (wfn_type_ == "OMP2") {
            timer_on("MP2L Energy");
            ccl_energy();
            timer_off("MP2L Energy");
        }

        else if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            if (compute_ccl == "TRUE") {
                timer_on("MP3L Energy");
                ccl_energy();
                timer_off("MP3L Energy");
            }

            else {
                timer_on("REF Energy");
                ref_energy();
                timer_off("REF Energy");
                timer_on("MP3 Energy");
                mp3_energy();
                timer_off("MP3 Energy");
                Emp3L = Emp3;
                EcorrL = Emp3L - Escf;
                DE = Emp3L - Emp3L_old;
                Emp3L_old = Emp3L;
            }
        }

        else if (wfn_type_ == "OCEPA" || wfn_type_ == "OREMP") {
            if (compute_ccl == "TRUE" ) {
                timer_on("CEPAL Energy");
                ccl_energy();
                timer_off("CEPAL Energy");
            }

            else {
                timer_on("REF Energy");
                ref_energy();
                timer_off("REF Energy");
                timer_on("CEPA Energy");
                cepa_energy();
                timer_off("CEPA Energy");
                EcepaL = Ecepa;
                EcorrL = EcepaL - Escf;
                DE = EcepaL - EcepaL_old;
                EcepaL_old = EcepaL;
            }
        }

        /********************************************************************************************/
        /************************** new orbital gradient ********************************************/
        /********************************************************************************************/
        timer_on("MO Grad");
        mograd();
        timer_off("MO Grad");

        /********************************************************************************************/
        /************************** New orbital step ************************************************/
        /********************************************************************************************/
        timer_on("kappa orb rot");
        compute_orbital_step();
        timer_off("kappa orb rot");

        /********************************************************************************************/
        /************************** NEW amplitudes **************************************************/
        /********************************************************************************************/
        if (wfn_type_ == "OMP2") {
            timer_on("T2(1)");
            iterate_t2o1_amplitudes();
            timer_off("T2(1)");
        }

        else if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            timer_on("T2(1)");
            iterate_t2o1_amplitudes();
            timer_off("T2(1)");
            timer_on("T2(2)");
            t2_2nd_general();
            timer_off("T2(2)");
        }

        else if (wfn_type_ == "OCEPA") {
            timer_on("T2");
            t2_amps();
            timer_off("T2");
        }

        else if (wfn_type_ == "OREMP") {
            timer_on("T2");
            t2_amps_remp();
            timer_off("T2");
        }

        /********************************************************************************************/
        /************************** Print ***********************************************************/
        /********************************************************************************************/
        if (reference_ == "RESTRICTED") {
            nidp = nidpA;
            rms_wog = rms_wogA;
            biggest_mograd = biggest_mogradA;
            rms_kappa = rms_kappaA;
            biggest_kappa = biggest_kappaA;
        }

        else if (reference_ == "UNRESTRICTED") {
            nidp = MAX0(nidpA, nidpB);
            rms_wog = MAX0(rms_wogA, rms_wogB);
            biggest_mograd = MAX0(biggest_mogradA, biggest_mogradB);
            rms_kappa = MAX0(rms_kappaA, rms_kappaB);
            biggest_kappa = MAX0(biggest_kappaA, biggest_kappaB);
            rms_t2 = MAX0(rms_t2AA, rms_t2BB);
            rms_t2 = MAX0(rms_t2, rms_t2AB);
        }

        if (wfn_type_ == "OMP2")
            outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n", itr_occ, Emp2L, DE, rms_wog,
                            biggest_mograd, rms_t2);
        else if (wfn_type_ == "OMP3")
            outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n", itr_occ, Emp3L, DE, rms_wog,
                            biggest_mograd, rms_t2);
        else if (wfn_type_ == "OCEPA" || wfn_type_ == "OREMP")
            outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n", itr_occ, EcepaL, DE, rms_wog,
                            biggest_mograd, rms_t2);
        else if (wfn_type_ == "OMP2.5")
            outfile->Printf(" %3d     %12.10f  %12.2e   %12.2e     %12.2e    %12.2e \n", itr_occ, Emp3L, DE, rms_wog,
                            biggest_mograd, rms_t2);

        /********************************************************************************************/
        /********************************************************************************************/
        if (itr_occ >= mo_maxiter) {
            conver = 0;  // means MOs are NOT optimized
            break;
        }

        if (rms_wog < tol_grad && biggest_mograd < mograd_max && std::fabs(DE) < tol_Eod) break;

        if (rms_wog >= DIVERGE) {
            throw PSIEXCEPTION("OCC iterations are diverging");
        }

        // Now it's time for DIIS.
        oo_diis(orbital_diis);

        // Handle any needed resorts of amplitudes, after the DIIS.
        if (wfn_type_ == "OCEPA" || wfn_type_ == "OREMP") {
            cepa_chemist();
        } else if (wfn_type_ == "OMP2.5" || wfn_type_ == "OMP3") {
            iterative_mp_postdiis_amplitudes();
        }

    } while (true); // TODO: This is a very silly do while loop.
    // while(std::fabs(DE) >= tol_Eod || rms_wog >= tol_grad || rms_kappa >= tol_grad || biggest_mograd >= mograd_max ||
    //      biggest_kappa >= mograd_max || rms_t2 >= tol_t2);

    if (conver == 1) {
        mo_optimized = 1;
        outfile->Printf("\n");
        outfile->Printf(" ============================================================================== \n");
        if (wfn_type_ == "OMP2")
            outfile->Printf(" ======================== OMP2 ITERATIONS ARE CONVERGED ======================= \n");
        else if (wfn_type_ == "OMP3")
            outfile->Printf(" ======================== OMP3 ITERATIONS ARE CONVERGED ======================= \n");
        else if (wfn_type_ == "OCEPA")
            outfile->Printf(" ======================== OLCCD ITERATIONS ARE CONVERGED ====================== \n");
        else if (wfn_type_ == "OREMP")
            outfile->Printf(" ======================== OREMP ITERATIONS ARE CONVERGED ====================== \n");
        else if (wfn_type_ == "OMP2.5")
            outfile->Printf(" ======================== OMP2.5 ITERATIONS ARE CONVERGED ===================== \n");
        outfile->Printf(" ============================================================================== \n");

    }

    else if (conver == 0) {
        if (wfn_type_ == "OMP2")
            outfile->Printf("\n ======================== OMP2 IS NOT CONVERGED IN %2d ITERATIONS ============= \n",
                            mo_maxiter);
        else if (wfn_type_ == "OMP3")
            outfile->Printf("\n ======================== OMP3 IS NOT CONVERGED IN %2d ITERATIONS ============= \n",
                            mo_maxiter);
        else if (wfn_type_ == "OCEPA")
            outfile->Printf("\n ======================== OLCCD IS NOT CONVERGED IN %2d ITERATIONS ============ \n",
                            mo_maxiter);
        else if (wfn_type_ == "OREMP")
            outfile->Printf("\n ======================== OREMP IS NOT CONVERGED IN %2d ITERATIONS ============ \n",
                            mo_maxiter);
        else if (wfn_type_ == "OMP2.5")
            outfile->Printf("\n ======================== OMP2.5 IS NOT CONVERGED IN %2d ITERATIONS =========== \n",
                            mo_maxiter);

        throw PSIEXCEPTION("OCC iterations did not converge");
    }
    // Clean up!
    delete[] idprowA;
    delete[] idpcolA;
    delete[] idpirrA;
    delete wogA;
    delete wog_intA;
    delete kappaA;
    delete kappa_newA;

    if (reference_ == "UNRESTRICTED") {
        delete[] idprowB;
        delete[] idpcolB;
        delete[] idpirrB;
        delete wogB;
        delete wog_intB;
        delete kappaB;
        delete kappa_newB;
    }

    // Clean up the mess of ORB-RESP
    if (opt_method == "ORB_RESP") {
        if (orb_resp_solver_ == "LINEQ") {
            if (reference_ == "UNRESTRICTED") {
                delete kappa;
            }
        }

        // PCG
        else if (orb_resp_solver_ == "PCG") {
            delete r_pcgA;
            delete S_pcgA;
            delete D_pcgA;
            delete sigma_pcgA;
            delete Minv_pcgA;
            delete r_pcg_newA;
            if (pcg_beta_type_ == "POLAK_RIBIERE") delete dr_pcgA;
            if (reference_ == "UNRESTRICTED") {
                delete r_pcgB;
                delete S_pcgB;
                delete D_pcgB;
                delete sigma_pcgB;
                delete Minv_pcgB;
                delete r_pcg_newB;
                if (pcg_beta_type_ == "POLAK_RIBIERE") delete dr_pcgB;
            }
        }
    }  // end orb resp if
}

void OCCWave::response_pdms() {
    if (wfn_type_ == "OMP2")
        omp2_response_pdms();
    else if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5")
        omp3_response_pdms();
    else if (wfn_type_ == "OCEPA" || wfn_type_ == "OREMP" )
        ocepa_response_pdms();
}

void OCCWave::compute_orbital_step() {
    // Compute kappa, the update to orbital amplitudes.
    if (opt_method == "ORB_RESP") {
        if (orb_resp_solver_ == "LINEQ")
            kappa_orb_resp();
        else if (orb_resp_solver_ == "PCG")
            kappa_orb_resp_iter();
    } else if (opt_method == "MSD")
        kappa_msd();

    // Update the orbital amplitude. Trust the rest of the program to DIIS when ready.
    // DIIS is tied to both orbitals and amplitudes, so having stepped orbitals isn't enough.
    const auto kappaA_vec = std::make_shared<Vector>(idp_dimensions_[SpinType::Alpha], *kappaA);
    kappa_bar_[SpinType::Alpha]->add(*kappaA_vec);
    if (reference_ == "UNRESTRICTED") {
        const auto kappaB_vec = std::make_shared<Vector>(idp_dimensions_[SpinType::Beta], *kappaB);
        kappa_bar_[SpinType::Beta]->add(*kappaB_vec);
    }
}

void OCCWave::oo_diis(DIISManager& orbital_diis) {
    if (!do_diis_) return;

    psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);

    const auto wogA_vec = std::make_shared<Vector>(idp_dimensions_[SpinType::Alpha], *wogA);
    if (reference_ == "RESTRICTED") {
        if (wfn_type_ == "OMP2") {
            dpdbuf4 T, R;
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T <OO|VV>");
            global_dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "RT2_1 <OO|VV>");
            orbital_diis.add_entry(wogA_vec.get(), &R, kappa_bar_[SpinType::Alpha].get(), &T);
            if (orbital_diis.subspace_size() >= mindiis_) {
                orbital_diis.extrapolate(kappa_bar_[SpinType::Alpha].get(), &T);
            }
        } else if (wfn_type_ == "OMP2.5" || wfn_type_ == "OMP3") {
            dpdbuf4 T1, R1, T2, R2;
            global_dpd_->buf4_init(&T1, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_1 <OO|VV>");
            global_dpd_->buf4_init(&R1, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "RT2_1 <OO|VV>");
            global_dpd_->buf4_init(&T2, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2 <OO|VV>");
            global_dpd_->buf4_init(&R2, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "RT2_2 <OO|VV>");
            orbital_diis.add_entry(wogA_vec.get(), &R1, &R2, kappa_bar_[SpinType::Alpha].get(), &T1, &T2);
            if (orbital_diis.subspace_size() >= mindiis_) {
                orbital_diis.extrapolate(kappa_bar_[SpinType::Alpha].get(), &T1, &T2);
            }
        } else if (wfn_type_ == "OCEPA" || wfn_type_ == "OREMP") {
            dpdbuf4 T, R;
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2 <OO|VV>");
            global_dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "RT2 <OO|VV>");
            orbital_diis.add_entry(wogA_vec.get(), &R, kappa_bar_[SpinType::Alpha].get(), &T);
            if (orbital_diis.subspace_size() >= mindiis_) {
                orbital_diis.extrapolate(kappa_bar_[SpinType::Alpha].get(), &T);
            }
        }
    } else if (reference_ == "UNRESTRICTED") {
        const auto wogB_vec = std::make_shared<Vector>(idp_dimensions_[SpinType::Beta], *wogB);
        if (wfn_type_ == "OMP2") {
            dpdbuf4 T1aa, T1ab, T1bb, R1aa, R1ab, R1bb;
            global_dpd_->buf4_init(&T1aa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_1 <OO|VV>");
            global_dpd_->buf4_init(&R1aa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "RT2_1 <OO|VV>");
            global_dpd_->buf4_init(&T1ab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_1 <Oo|Vv>");
            global_dpd_->buf4_init(&R1ab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "RT2_1 <Oo|Vv>");
            global_dpd_->buf4_init(&T1bb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_1 <oo|vv>");
            global_dpd_->buf4_init(&R1bb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "RT2_1 <oo|vv>");
            orbital_diis.add_entry(wogA_vec.get(), wogB_vec.get(), &R1aa, &R1ab, &R1bb,
                    kappa_bar_[SpinType::Alpha].get(), kappa_bar_[SpinType::Beta].get(), &T1aa, &T1ab, &T1bb);
            if (orbital_diis.subspace_size() >= mindiis_) {
                orbital_diis.extrapolate(kappa_bar_[SpinType::Alpha].get(), kappa_bar_[SpinType::Beta].get(), &T1aa, &T1ab, &T1bb);
            }
            global_dpd_->buf4_close(&T1aa);
            global_dpd_->buf4_close(&T1ab);
            global_dpd_->buf4_close(&T1bb);
            global_dpd_->buf4_close(&R1aa);
            global_dpd_->buf4_close(&R1ab);
            global_dpd_->buf4_close(&R1bb);
        } else if (wfn_type_ == "OMP2.5" || wfn_type_ == "OMP3") {
            dpdbuf4 T1aa, T1ab, T1bb, R1aa, R1ab, R1bb, T2aa, T2ab, T2bb, R2aa, R2ab, R2bb;
            global_dpd_->buf4_init(&T1aa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_1 <OO|VV>");
            global_dpd_->buf4_init(&R1aa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "RT2_1 <OO|VV>");
            global_dpd_->buf4_init(&T1ab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_1 <Oo|Vv>");
            global_dpd_->buf4_init(&R1ab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "RT2_1 <Oo|Vv>");
            global_dpd_->buf4_init(&T1bb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_1 <oo|vv>");
            global_dpd_->buf4_init(&R1bb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "RT2_1 <oo|vv>");
            global_dpd_->buf4_init(&T2aa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2_2 <OO|VV>");
            global_dpd_->buf4_init(&R2aa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "RT2_2 <OO|VV>");
            global_dpd_->buf4_init(&T2ab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2_2 <Oo|Vv>");
            global_dpd_->buf4_init(&R2ab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "RT2_2 <Oo|Vv>");
            global_dpd_->buf4_init(&T2bb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2_2 <oo|vv>");
            global_dpd_->buf4_init(&R2bb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "RT2_2 <oo|vv>");
            orbital_diis.add_entry(wogA_vec.get(), wogB_vec.get(), &R1aa, &R1ab, &R1bb, &R2aa, &R2ab, &R2bb,
                    kappa_bar_[SpinType::Alpha].get(), kappa_bar_[SpinType::Beta].get(), &T1aa, &T1ab, &T1bb, &T2aa, &T2ab, &T2bb);
            if (orbital_diis.subspace_size() >= mindiis_) {
                orbital_diis.extrapolate(kappa_bar_[SpinType::Alpha].get(), kappa_bar_[SpinType::Beta].get(), &T1aa, &T1ab, &T1bb, &T2aa, &T2ab, &T2bb);
            }
            global_dpd_->buf4_close(&T1aa);
            global_dpd_->buf4_close(&T1ab);
            global_dpd_->buf4_close(&T1bb);
            global_dpd_->buf4_close(&R1aa);
            global_dpd_->buf4_close(&R1ab);
            global_dpd_->buf4_close(&R1bb);
            global_dpd_->buf4_close(&T2aa);
            global_dpd_->buf4_close(&T2ab);
            global_dpd_->buf4_close(&T2bb);
            global_dpd_->buf4_close(&R2aa);
            global_dpd_->buf4_close(&R2ab);
            global_dpd_->buf4_close(&R2bb);
        } else if (wfn_type_ == "OCEPA" || wfn_type_ == "OREMP") {
            dpdbuf4 Taa, Tab, Tbb, Raa, Rab, Rbb;
            global_dpd_->buf4_init(&Taa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2 <OO|VV>");
            global_dpd_->buf4_init(&Raa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "RT2 <OO|VV>");
            global_dpd_->buf4_init(&Tab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2 <Oo|Vv>");
            global_dpd_->buf4_init(&Rab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "RT2 <Oo|Vv>");
            global_dpd_->buf4_init(&Tbb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2 <oo|vv>");
            global_dpd_->buf4_init(&Rbb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "RT2 <oo|vv>");
            orbital_diis.add_entry(wogA_vec.get(), wogB_vec.get(), &Raa, &Rab, &Rbb,
                    kappa_bar_[SpinType::Alpha].get(), kappa_bar_[SpinType::Beta].get(), &Taa, &Tab, &Tbb);
            if (orbital_diis.subspace_size() >= mindiis_) {
                orbital_diis.extrapolate(kappa_bar_[SpinType::Alpha].get(), kappa_bar_[SpinType::Beta].get(), &Taa, &Tab, &Tbb);
            }
            global_dpd_->buf4_close(&Taa);
            global_dpd_->buf4_close(&Tab);
            global_dpd_->buf4_close(&Tbb);
            global_dpd_->buf4_close(&Raa);
            global_dpd_->buf4_close(&Rab);
            global_dpd_->buf4_close(&Rbb);
        }
    }
    psio_->close(PSIF_OCC_DPD, 1);
}
}   // namespace occwave
}  // End Namespaces
