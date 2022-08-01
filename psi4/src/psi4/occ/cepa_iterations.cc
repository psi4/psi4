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

namespace psi {
namespace occwave {

void OCCWave::cepa_iterations() {
    outfile->Printf("\n  \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ================ Performing LCCD iterations... =============================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
    outfile->Printf("  Iter    E_corr           E_total            DE           T2 RMS        \n");
    outfile->Printf("  ----   -------------    ---------------    ----------   ----------    \n");

    /********************************************************************************************/
    /************************** NR iterations **************************************************/
    /********************************************************************************************/
    itr_occ = 0;
    conver = 1;  // Assuming that the iterations will converge
                 // DIIS
    DIISManager t2_diis;
    if (nooA + nooB != 1) {
        if (reference_ == "RESTRICTED") {
            dpdbuf4 T;
            psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                                   "T2 <OO|VV>");
            t2_diis = DIISManager(maxdiis_, "CEPA DIIS T2 Amps", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::OnDisk);
            t2_diis.set_error_vector_size(&T);
            t2_diis.set_vector_size(&T);
            global_dpd_->buf4_close(&T);
            psio_->close(PSIF_OCC_DPD, 1);
        }

        else if (reference_ == "UNRESTRICTED") {
            dpdbuf4 Taa, Tbb, Tab;
            psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
            global_dpd_->buf4_init(&Taa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                                   "T2 <OO|VV>");
            global_dpd_->buf4_init(&Tbb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                                   "T2 <oo|vv>");
            global_dpd_->buf4_init(&Tab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                                   "T2 <Oo|Vv>");
            t2_diis = DIISManager(maxdiis_, "CEPA DIIS T2 Amps", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);
            t2_diis.set_error_vector_size(&Taa, &Tbb, &Tab);
            t2_diis.set_vector_size(&Taa, &Tbb, &Tab);
            global_dpd_->buf4_close(&Taa);
            global_dpd_->buf4_close(&Tbb);
            global_dpd_->buf4_close(&Tab);
            psio_->close(PSIF_OCC_DPD, 1);
        }
    }

    // head of loop
    do {
        itr_occ++;
        timer_on("T2");
        t2_amps();
        timer_off("T2");
        timer_on("CEPA Energy");
        cepa_energy();
        cepa_diis(t2_diis);
        cepa_chemist();
        timer_off("CEPA Energy");
        Ecorr = Ecepa - Escf;
        DE = Ecepa - Ecepa_old;
        Ecepa_old = Ecepa;

        if (reference_ == "UNRESTRICTED") {
            rms_t2 = MAX0(rms_t2AA, rms_t2BB);
            rms_t2 = MAX0(rms_t2, rms_t2AB);
        }

        outfile->Printf(" %3d     %12.10f    %12.10f  %12.2e %12.2e \n", itr_occ, Ecorr, Ecepa, DE, rms_t2);

        if (itr_occ >= cc_maxiter) {
            conver = 0;  // means iterations were NOT converged
            break;
        }

        if (rms_t2 >= DIVERGE) {
            throw PSIEXCEPTION("LCCD iterations are diverging");
        }

    } while (std::fabs(DE) >= (0.5 * tol_Eod) || rms_t2 >= tol_t2);
    // 0.5 scale battens down a touch tighter for spin components since tol_Eod can be satisfied by small energy increase

    if (conver == 1) {
        EcepaL = Ecepa;
        outfile->Printf("\n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ======================== LCCD ITERATIONS ARE CONVERGED ======================= \n");
        outfile->Printf(" ============================================================================== \n");

    }

    else if (conver == 0) {
        outfile->Printf("\n ======================= LCCD IS NOT CONVERGED IN %2d ITERATIONS ============ \n",
                        cc_maxiter);

        throw PSIEXCEPTION("LCCD iterations did not converge");
    }

}  // end main

void OCCWave::remp_iterations() {
    outfile->Printf("\n  \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ================ Performing REMP iterations ... ============================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
    outfile->Printf("  Iter    E_corr           E_total            DE           T2 RMS        \n");
    outfile->Printf("  ----   -------------    ---------------    ----------   ----------    \n");

    /********************************************************************************************/
    /************************** NR iterations **************************************************/
    /********************************************************************************************/
    itr_occ = 0;
    conver = 1;  // Assuming that the iterations will converge
                 // DIIS
    DIISManager t2_diis;
    if (nooA + nooB != 1) {
        if (reference_ == "RESTRICTED") {
            dpdbuf4 T;
            psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                                   "T2 <OO|VV>");
            t2_diis = DIISManager(maxdiis_, "CEPA DIIS T2 Amps", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::OnDisk);
            t2_diis.set_error_vector_size(&T);
            t2_diis.set_vector_size(&T);
            global_dpd_->buf4_close(&T);
            psio_->close(PSIF_OCC_DPD, 1);
        }

        else if (reference_ == "UNRESTRICTED") {
            dpdbuf4 Taa, Tbb, Tab;
            psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
            global_dpd_->buf4_init(&Taa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                                   "T2 <OO|VV>");
            global_dpd_->buf4_init(&Tbb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                                   "T2 <oo|vv>");
            global_dpd_->buf4_init(&Tab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                                   "T2 <Oo|Vv>");
            t2_diis = DIISManager(maxdiis_, "CEPA DIIS T2 Amps", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);
            t2_diis.set_error_vector_size(&Taa, &Tbb, &Tab);
            t2_diis.set_vector_size(&Taa, &Tbb, &Tab);
            global_dpd_->buf4_close(&Taa);
            global_dpd_->buf4_close(&Tbb);
            global_dpd_->buf4_close(&Tab);
            psio_->close(PSIF_OCC_DPD, 1);
        }
    }

    // head of loop
    do {
        itr_occ++;
        timer_on("T2");
        t2_amps_remp(); // <- the only actual modification compared to regular CEPA(0)/D
        timer_off("T2");
        timer_on("CEPA Energy");
        cepa_energy();
        cepa_diis(t2_diis);    // <- CEPA diis can be reused without modifications
        cepa_chemist();
        timer_off("CEPA Energy");
        Ecorr = Ecepa - Escf;
        DE = Ecepa - Ecepa_old;
        Ecepa_old = Ecepa;

        if (reference_ == "UNRESTRICTED") {
            rms_t2 = MAX0(rms_t2AA, rms_t2BB);
            rms_t2 = MAX0(rms_t2, rms_t2AB);
        }

        outfile->Printf(" %3d     %12.10f    %12.10f  %12.2e %12.2e \n", itr_occ, Ecorr, Ecepa, DE, rms_t2);

        if (itr_occ >= cc_maxiter) {
            conver = 0;  // means iterations were NOT converged
            break;
        }

        if (rms_t2 >= DIVERGE) {
            throw PSIEXCEPTION("REMP iterations are diverging");
        }

    } while (std::fabs(DE) >= (0.5 * tol_Eod) || rms_t2 >= tol_t2);
    // 0.5 scale battens down a touch tighter for spin components since tol_Eod can be satisfied by small energy increase

    if (conver == 1) {
        EcepaL = Ecepa;
        outfile->Printf("\n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ======================== REMP ITERATIONS ARE CONVERGED ======================= \n");
        outfile->Printf(" ============================================================================== \n");

    }

    else if (conver == 0) {
        outfile->Printf("\n ======================= REMP IS NOT CONVERGED IN %2d ITERATIONS ============ \n",
                        cc_maxiter);

        throw PSIEXCEPTION("REMP iterations did not converge");
    }

}  // end main


void OCCWave::cepa_diis(DIISManager& t2_diis) {
    psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
    if (reference_ == "RESTRICTED") {
        dpdbuf4 R, T;
        global_dpd_->buf4_init(&R, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "RT2 <OO|VV>");
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2 <OO|VV>");
        t2_diis.add_entry(&R, &T);
        if (t2_diis.subspace_size() >= mindiis_) t2_diis.extrapolate(&T);
        global_dpd_->buf4_close(&R);
        global_dpd_->buf4_close(&T);
    } else if (reference_ == "UNRESTRICTED") {
        dpdbuf4 Raa, Rbb, Rab, Taa, Tbb, Tab;
        global_dpd_->buf4_init(&Raa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "RT2 <OO|VV>");
        global_dpd_->buf4_init(&Taa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2 <OO|VV>");
        global_dpd_->buf4_init(&Rbb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "RT2 <oo|vv>");
        global_dpd_->buf4_init(&Tbb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2 <oo|vv>");
        global_dpd_->buf4_init(&Rab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "RT2 <Oo|Vv>");
        global_dpd_->buf4_init(&Tab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2 <Oo|Vv>");
        t2_diis.add_entry(&Raa, &Rbb, &Rab, &Taa, &Tbb, &Tab);
        if (t2_diis.subspace_size() >= mindiis_) t2_diis.extrapolate(&Taa, &Tbb, &Tab);
        global_dpd_->buf4_close(&Raa);
        global_dpd_->buf4_close(&Rbb);
        global_dpd_->buf4_close(&Rab);
        global_dpd_->buf4_close(&Taa);
        global_dpd_->buf4_close(&Tbb);
        global_dpd_->buf4_close(&Tab);
    }
    psio_->close(PSIF_OCC_DPD, 1);
}

void OCCWave::cepa_chemist() {
    // Build amplitudes in chemist notation
    psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
    if (reference_ == "RESTRICTED") {
        dpdbuf4 T, Tau, Ttemp, Tp;
        // Build Tau(ij,ab) = 2*T(ij,ab) - T(ji,ab)
        // Build TAA(ij,ab) = T(ij,ab) - T(ji,ab)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2 <OO|VV>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "Tau <OO|VV>");
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, "T2AA <OO|VV>");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qprs, ID("[O,O]"), ID("[V,V]"), "T2jiab <OO|VV>");
        global_dpd_->buf4_close(&T);
        global_dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Tau <OO|VV>");
        global_dpd_->buf4_init(&Tp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2AA <OO|VV>");
        global_dpd_->buf4_init(&Ttemp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2jiab <OO|VV>");
        global_dpd_->buf4_scm(&Tau, 2.0);
        global_dpd_->buf4_axpy(&Ttemp, &Tau, -1.0);  // -1.0*Ttemp + Tau -> Tau
        global_dpd_->buf4_axpy(&Ttemp, &Tp, -1.0);   // -1.0*Ttemp + Tp -> Tp
        global_dpd_->buf4_close(&Ttemp);
        global_dpd_->buf4_close(&Tp);
        global_dpd_->buf4_close(&Tau);

        // T_IJ^AB => T'(IA,JB), T"(JA,IB)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2 <OO|VV>");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "T2 (OV|OV)");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qrps, ID("[O,V]"), ID("[O,V]"), "T2pp (OV|OV)");
        global_dpd_->buf4_close(&T);

        // Tau(IJ,AB) => Tau'(IA,JB), Tau"(JA,IB)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "Tau <OO|VV>");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "Tau (OV|OV)");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qrps, ID("[O,V]"), ID("[O,V]"), "Taupp (OV|OV)");
        global_dpd_->buf4_close(&T);
    } else if (reference_ == "UNRESTRICTED") {
        dpdbuf4 T;
        // T_IJ^AB => T(IA,JB)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "T2 <OO|VV>");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "T2 (OV|OV)");
        global_dpd_->buf4_close(&T);

        // T_ij^ab => T(ia,jb)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "T2 <oo|vv>");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[o,v]"), ID("[o,v]"), "T2 (ov|ov)");
        global_dpd_->buf4_close(&T);

        // T_Ij^Ab => T(IA,jb), T(jA,Ib)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "T2 <Oo|Vv>");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,V]"), ID("[o,v]"), "T2 (OV|ov)");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qrps, ID("[o,V]"), ID("[O,v]"), "T2 (oV|Ov)");
        global_dpd_->buf4_close(&T);

        // T(IA,jb) => T(jb,IA)
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                               "T2 (OV|ov)");
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, rspq, ID("[o,v]"), ID("[O,V]"), "T2 (ov|OV)");
        global_dpd_->buf4_close(&T);
    }
    // close files
    psio_->close(PSIF_OCC_DPD, 1);
}
}
}  // End Namespaces
