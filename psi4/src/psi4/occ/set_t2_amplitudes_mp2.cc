/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "defines.h"
#include "occwave.h"

using namespace psi;

namespace psi {
namespace occwave {

// TODO: A cleaner solution would be to have the manager method pass the suffix as a variable.
// That will only be cleaner if RHF/UHF have the same variable names. That is not the case for MP2.
void OCCWave::set_t2_amplitudes_mp2() {
    //===========================================================================================
    //========================= RHF =============================================================
    //===========================================================================================
    if (reference_ == "RESTRICTED") {
        dpdbuf4 K, T, D, Tau, Ttemp, Tss;

        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);

        // Each method has its own names for the variables, so set the correct ones.
        std::string temp = (wfn_type_ == "OCEPA" || wfn_type_ == "OREMP") ? "2" : (wfn_type_ == "OMP2" ? "" : "2_1");
        std::string t = "T" + temp;
        temp = (wfn_type_ == "OCEPA" || wfn_type_ == "OREMP") ? "" : (wfn_type_ == "OMP2" ? "" : "_1");
        std::string tau = "Tau" + temp;
        std::string t_name = t + " <OO|VV>";
        std::string tau_name = tau + " <OO|VV>";
        std::string taa_name = t + "AA <OO|VV>";
        std::string jiab_name = t + "jiab <OO|VV>";
        std::string t_chem = t + " (OV|OV)";
        std::string tpp_chem = t + "pp (OV|OV)";
        std::string tau_chem = tau + " (OV|OV)";
        std::string taupp_chem = tau + "pp (OV|OV)";

        // T_ij^ab = <ij|ab>
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "MO Ints <OO|VV>");
        global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, t_name.c_str());
        global_dpd_->buf4_close(&K);

        // T_ij^ab = T_ij^ab / D_ij^ab
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "D <OO|VV>");
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               t_name.c_str());
        global_dpd_->buf4_dirprd(&D, &T);
        global_dpd_->buf4_close(&D);

        // Build Tau(ij,ab) = 2*T(ij,ab) - T(ji,ab)
        // Build TAA(ij,ab) = T(ij,ab) - T(ji,ab)
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, tau_name.c_str());
        global_dpd_->buf4_copy(&T, PSIF_OCC_DPD, taa_name.c_str());
        global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qprs, ID("[O,O]"), ID("[V,V]"), jiab_name.c_str());
        global_dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               tau_name.c_str());
        global_dpd_->buf4_init(&Tss, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               taa_name.c_str());
        global_dpd_->buf4_init(&Ttemp, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               jiab_name.c_str());
        global_dpd_->buf4_scm(&Tau, 2.0);
        global_dpd_->buf4_axpy(&Ttemp, &Tau, -1.0);  // -1.0*Ttemp + Tau -> Tau
        global_dpd_->buf4_axpy(&Ttemp, &Tss, -1.0);  // -1.0*Ttemp + Tss -> Tss
        global_dpd_->buf4_close(&Ttemp);
        global_dpd_->buf4_close(&Tau);
        global_dpd_->buf4_close(&Tss);

        if (print_ > 4) global_dpd_->buf4_print(&T, "outfile", 1);
        global_dpd_->buf4_close(&T);

        if (wfn_type_ != "OMP2") {
            // Build amplitudes in chemist notation
            // T_IJ^AB => T'(IA,JB), T"(JA,IB)
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                                   t_name.c_str());
            global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,V]"), ID("[O,V]"), t_chem.c_str());
            global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qrps, ID("[O,V]"), ID("[O,V]"), tpp_chem.c_str());
            global_dpd_->buf4_close(&T);

            // Tau(IJ,AB) => Tau'(IA,JB), Tau"(JA,IB)
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                                   tau_name.c_str());
            global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,V]"), ID("[O,V]"), tau_chem.c_str());
            global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qrps, ID("[O,V]"), ID("[O,V]"), taupp_chem.c_str());
            global_dpd_->buf4_close(&T);
        }

        psio_->close(PSIF_LIBTRANS_DPD, 1);
        psio_->close(PSIF_OCC_DPD, 1);

    }  // end if (reference_ == "RESTRICTED")

    //===========================================================================================
    //========================= UHF =============================================================
    //===========================================================================================
    else if (reference_ == "UNRESTRICTED") {
        dpdbuf4 K, T, D;

        psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
        psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);

        std::string suffix = (wfn_type_ == "OCEPA" || wfn_type_ == "OREMP") ? "" : "_1";
        std::string taa_name = "T2" + suffix + " <OO|VV>";
        std::string tab_name = "T2" + suffix + " <Oo|Vv>";
        std::string tbb_name = "T2" + suffix + " <oo|vv>";
        std::string taa_c1 = "T2" + suffix + " (OV|OV)";
        std::string tbb_c1 = "T2" + suffix + " (ov|ov)";
        std::string tab_c1 = "T2" + suffix + " (OV|ov)";
        std::string tab_c2 = "T2" + suffix + " (oV|Ov)";
        std::string tab_c3 = "T2" + suffix + " (ov|OV)";

        // Build T2AA
        // T_IJ^AB = <IJ||AB>
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "MO Ints <OO||VV>");
        global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, taa_name.c_str());
        global_dpd_->buf4_close(&K);

        // T_IJ^AB = T_IJ^AB / D_IJ^AB
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               "D <OO|VV>");
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                               taa_name.c_str());
        global_dpd_->buf4_dirprd(&D, &T);
        global_dpd_->buf4_close(&D);
        if (print_ > 1) global_dpd_->buf4_print(&T, "outfile", 1);
        global_dpd_->buf4_close(&T);

        // Build T2BB
        // T_ij^ab = <ij|ab>
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "MO Ints <oo||vv>");
        global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, tbb_name.c_str());
        global_dpd_->buf4_close(&K);

        // T_ij^ab = T_ij^ab / D_ij^ab
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               "D <oo|vv>");
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                               tbb_name.c_str());
        global_dpd_->buf4_dirprd(&D, &T);
        global_dpd_->buf4_close(&D);
        if (print_ > 1) global_dpd_->buf4_print(&T, "outfile", 1);
        global_dpd_->buf4_close(&T);

        // Build T2AB
        // T_Ij^Ab = <Ij|Ab>
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "MO Ints <Oo|Vv>");
        global_dpd_->buf4_copy(&K, PSIF_OCC_DPD, tab_name.c_str());
        global_dpd_->buf4_close(&K);

        // T_Ij^Ab = T_Ij^Ab / D_Ij^Ab
        global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               "D <Oo|Vv>");
        global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                               tab_name.c_str());
        global_dpd_->buf4_dirprd(&D, &T);
        global_dpd_->buf4_close(&D);
        if (print_ > 1) global_dpd_->buf4_print(&T, "outfile", 1);
        global_dpd_->buf4_close(&T);

        if (wfn_type_ != "OMP2") {
            // Build amplitudes in chemist notation
            // T_IJ^AB => T(IA,JB)
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"), ID("[O,O]"), ID("[V,V]"), 0,
                                   taa_name.c_str());
            global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,V]"), ID("[O,V]"), taa_c1.c_str());
            global_dpd_->buf4_close(&T);

            // T_ij^ab => T(ia,jb)
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"), ID("[o,o]"), ID("[v,v]"), 0,
                                   tbb_name.c_str());
            global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[o,v]"), ID("[o,v]"), tbb_c1.c_str());
            global_dpd_->buf4_close(&T);

            // T_Ij^Ab => T(IA,jb), T(jA,Ib)
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"), ID("[O,o]"), ID("[V,v]"), 0,
                                   tab_name.c_str());
            global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, prqs, ID("[O,V]"), ID("[o,v]"), tab_c1.c_str());
            global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, qrps, ID("[o,V]"), ID("[O,v]"), tab_c2.c_str());
            global_dpd_->buf4_close(&T);

            // T(IA,jb) => T(jb,IA)
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,V]"), ID("[o,v]"), ID("[O,V]"), ID("[o,v]"), 0,
                                   tab_c1.c_str());
            global_dpd_->buf4_sort(&T, PSIF_OCC_DPD, rspq, ID("[o,v]"), ID("[O,V]"), tab_c3.c_str());
            global_dpd_->buf4_close(&T);
        }

        psio_->close(PSIF_LIBTRANS_DPD, 1);
        psio_->close(PSIF_OCC_DPD, 1);
    }  // end if (reference_ == "UNRESTRICTED")

}  // end t2_1st_sc
}
}  // End Namespaces
