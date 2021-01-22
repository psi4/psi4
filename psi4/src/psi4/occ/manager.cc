/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#include "occwave.h"

#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/process.h"

#include <cmath>

using namespace psi;

namespace psi {
namespace occwave {

//======================================================================
//             OMP2 Manager
//======================================================================
void OCCWave::omp2_manager() {
    mo_optimized = 0;
    orbs_already_opt = 0;
    orbs_already_sc = 0;
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("REF Energy");
    ref_energy();
    timer_off("REF Energy");
    timer_on("MP2 Energy");
    mp2_energy();
    timer_off("MP2 Energy");
    Emp2L = Emp2;
    EcorrL = Emp2L - Escf;
    Emp2L_old = Emp2;
    if (ip_poles == "TRUE") omp2_ip_poles();
    if (ep_ip_poles == "TRUE") ep2_ip();

    mp2_postprocessing();

    occ_iterations();

    if (conver == 1) {
        ref_energy();
        mp2_energy();
        if (orbs_already_opt == 1) Emp2L = Emp2;

        // S2
        // if (comput_s2_ == "TRUE" && reference_ == "UNRESTRICTED") s2_response();

        // Green's function
        if (ip_poles == "TRUE") {
            if (orbs_already_sc == 0) {
                semi_canonic();
                if (reference_ == "RESTRICTED")
                    trans_ints_rhf();
                else if (reference_ == "UNRESTRICTED")
                    trans_ints_uhf();
                set_t2_amplitudes_mp2();
            }
            omp2_ip_poles();
        }

        if (ep_ip_poles == "TRUE") {
            if (orbs_already_sc == 0) {
                semi_canonic();
                if (reference_ == "RESTRICTED")
                    trans_ints_rhf();
                else if (reference_ == "UNRESTRICTED")
                    trans_ints_uhf();
            }
            ep2_ip();
        }

        // EKT
        if (ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
            if (orbs_already_sc == 1) {
                omp2_response_pdms();
                gfock();
            }
            gfock_diag();
            if (ekt_ip_ == "TRUE") ekt_ip();
            if (ekt_ea_ == "TRUE") ekt_ea();
        }

        mp2_printing();

        outfile->Printf("\n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\t================ OMP2 FINAL RESULTS ========================================== \n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tSCS-OMP2 Total Energy (a.u.)       : %20.14f\n", Escsmp2);
        outfile->Printf("\tSOS-OMP2 Total Energy (a.u.)       : %20.14f\n", Esosmp2);
        outfile->Printf("\tOMP2 Correlation Energy (a.u.)     : %20.14f\n", Emp2L - Escf);
        outfile->Printf("\tEomp2 - Eref (a.u.)                : %20.14f\n", Emp2L - Eref);
        outfile->Printf("\tOMP2 Total Energy (a.u.)           : %20.14f\n", Emp2L);
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\n");

        // Set the global variables with the energies
        variables_["OMP2 TOTAL ENERGY"] = Emp2L;
        variables_["OMP2 TOTAL ENERGY"] = Emp2L;
        variables_["SCS-OMP2 TOTAL ENERGY"] = Escsmp2;
        variables_["SOS-OMP2 TOTAL ENERGY"] = Esosmp2;

        variables_["CURRENT REFERENCE ENERGY"] = Escf;

        variables_["OMP2 CORRELATION ENERGY"] = Emp2L - Escf;
        variables_["SCS-OMP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        variables_["SOS-OMP2 CORRELATION ENERGY"] = Esosmp2 - Escf;


        variables_["CUSTOM SCS-OMP2 CORRELATION ENERGY"] = os_scale * Emp2AB + ss_scale * (Emp2AA + Emp2BB);
        variables_["CUSTOM SCS-OMP2 TOTAL ENERGY"] = Escf + variables_["CUSTOM SCS-OMP2 CORRELATION ENERGY"];

        std::map<std::string, double> spin_scale_energies = {{"SCS", Escsmp2},
                                                             {"SOS", Esosmp2},
                                                             {"CUSTOM", variables_["CUSTOM SCS-OMP2 TOTAL ENERGY"]},
                                                             {"NONE", Emp2L}};
        variables_["CURRENT ENERGY"] = spin_scale_energies[spin_scale_type_];
        energy_ = variables_["CURRENT ENERGY"];
        variables_["CURRENT CORRELATION ENERGY"] =
            variables_["CURRENT ENERGY"] - variables_["CURRENT REFERENCE ENERGY"];

        if (natorb == "TRUE") nbo();
        if (occ_orb_energy == "TRUE") semi_canonic();

        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            outfile->Printf("\tAnalytic gradient computation is starting...\n");

            coord_grad();
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }

    }  // end if (conver == 1)
}  // end omp2_manager

//======================================================================
//             MP2 Manager
//======================================================================
void OCCWave::mp2_manager() {
    time4grad = 0;  // means i will not compute the gradient
    timer_on("trans_ints");
    if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
        if (reference_ == "RESTRICTED")
            trans_ints_rhf();
        else if (reference_ == "UNRESTRICTED")
            trans_ints_uhf();
    } else {
        if (reference_ == "RESTRICTED")
            trans_ints_rmp2();
        else if (reference_ == "UNRESTRICTED" && reference == "ROHF")
            trans_ints_uhf();
        else if (reference_ == "UNRESTRICTED" && reference != "ROHF")
            trans_ints_ump2();
    }
    timer_off("trans_ints");
    // ROHF REF
    if (reference == "ROHF") {
        timer_on("T1(1)");
        t1_1st_sc();
        timer_off("T1(1)");
    }  // end if (reference == "ROHF")
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    Eref = Escf;
    timer_on("MP2 Energy");
    mp2_energy(reference == "ROHF");
    timer_off("MP2 Energy");
    Emp2L = Emp2;
    EcorrL = Emp2L - Escf;
    Emp2L_old = Emp2;
    if (ip_poles == "TRUE") omp2_ip_poles();
    if (ep_ip_poles == "TRUE") ep2_ip();

    mp2_postprocessing(reference == "ROHF");

    variables_["CURRENT REFERENCE ENERGY"] = Escf;
    std::map<std::string, double> spin_scale_energies = {
        {"SCS", Escsmp2}, {"SCSN", Escsnmp2},   {"SCSVDW", Escsmp2vdw},
        {"SOS", Esosmp2}, {"SOSPI", Esospimp2}, {"CUSTOM", variables_["CUSTOM SCS-MP2 TOTAL ENERGY"]},
        {"NONE", Emp2}};
    variables_["CURRENT ENERGY"] = spin_scale_energies[spin_scale_type_];
    energy_ = variables_["CURRENT ENERGY"];
    variables_["CURRENT CORRELATION ENERGY"] = variables_["CURRENT ENERGY"] - variables_["CURRENT REFERENCE ENERGY"];

    /* updates the wavefunction for checkpointing */
    name_ = "MP2";

    // S2
    // if (comput_s2_ == "TRUE" && reference_ == "UNRESTRICTED") s2_response();

    // Compute Analytic Gradients
    if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
        outfile->Printf("\tAnalytic gradient computation is starting...\n");
        outfile->Printf("\tComputing response density matrices...\n");

        omp2_response_pdms();
        outfile->Printf("\tComputing off-diagonal blocks of GFM...\n");

        gfock();
        outfile->Printf("\tForming independent-pairs...\n");

        idp2();
        outfile->Printf("\tComputing orbital gradient...\n");

        mograd();
        coord_grad();

        if (ekt_ip_ == "TRUE" && ekt_ea_ == "TRUE") {
            ekt_ip();
            ekt_ea();
            // outfile->Printf("\tAn EKT computation for a non-OO method requested. Analytic gradients will not be
            // computed! \n");
            // tstop();
            // exit(EXIT_SUCCESS);
        }

        else if (ekt_ip_ == "TRUE" && ekt_ea_ == "FALSE") {
            ekt_ip();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "TRUE") {
            ekt_ea();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") {
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }
    }

}  // end mp2_manager

//======================================================================
//             OMP3 Manager
//======================================================================
void OCCWave::omp3_manager() {
    mo_optimized = 0;
    orbs_already_opt = 0;
    orbs_already_sc = 0;
    time4grad = 0;  // means i will not compute the gradient
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    timer_on("REF Energy");
    ref_energy();
    timer_off("REF Energy");
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("MP2 Energy");
    mp2_energy();
    timer_off("MP2 Energy");

    mp2_postprocessing();

    timer_on("T2(2)");
    t2_2nd_sc();
    timer_off("T2(2)");
    timer_on("MP3 Energy");
    mp3_energy();
    timer_off("MP3 Energy");
    Emp3L = Emp3;
    EcorrL = Emp3L - Escf;
    Emp3L_old = Emp3;
    if (ip_poles == "TRUE") omp3_ip_poles();

    mp3_postprocessing();

    occ_iterations();

    if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod) {
        orbs_already_opt = 1;
        if (conver == 1)
            outfile->Printf("\n\tOrbitals are optimized now.\n");
        else if (conver == 0) {
            outfile->Printf("\n\tMAX MOGRAD did NOT converge, but RMS MOGRAD converged!!!\n");
            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
        }
        outfile->Printf("\tSwitching to the standard MP3 computation after semicanonicalization of the MOs... \n");

        semi_canonic();
        if (reference_ == "RESTRICTED")
            trans_ints_rhf();
        else if (reference_ == "UNRESTRICTED")
            trans_ints_uhf();
        set_t2_amplitudes_mp2();
        t2_2nd_sc();
        conver = 1;
        if (dertype == "FIRST") {
            omp3_response_pdms();
            gfock();
        }
    }

    if (conver == 1) {
        ref_energy();
        mp2_energy();
        mp3_energy();
        if (orbs_already_opt == 1) Emp3L = Emp3;

        if (ip_poles == "TRUE") {
            if (orbs_already_sc == 0) {
                semi_canonic();
                if (reference_ == "RESTRICTED")
                    trans_ints_rhf();
                else if (reference_ == "UNRESTRICTED")
                    trans_ints_uhf();
                set_t2_amplitudes_mp2();
                t2_2nd_sc();
            }
            omp3_ip_poles();
        }

        // EKT
        if (ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
            if (orbs_already_sc == 1) {
                omp3_response_pdms();
                gfock();
            }
            gfock_diag();
            if (ekt_ip_ == "TRUE") ekt_ip();
            if (ekt_ea_ == "TRUE") ekt_ea();
        }

        mp2_printing();

        mp3_printing();

        outfile->Printf("\n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\t================ OMP3 FINAL RESULTS ========================================== \n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tSCS-OMP3 Total Energy (a.u.)       : %20.14f\n", Escsmp3);
        outfile->Printf("\tSOS-OMP3 Total Energy (a.u.)       : %20.14f\n", Esosmp3);
        outfile->Printf("\tOMP3 Correlation Energy (a.u.)     : %20.14f\n", Emp3L - Escf);
        outfile->Printf("\tEomp3 - Eref (a.u.)                : %20.14f\n", Emp3L - Eref);
        outfile->Printf("\tOMP3 Total Energy (a.u.)           : %20.14f\n", Emp3L);
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\n");

        // Set the global variables with the energies
        variables_["OMP3 TOTAL ENERGY"] = Emp3L;
        variables_["SCS-OMP3 TOTAL ENERGY"] = Escsmp3;
        variables_["SOS-OMP3 TOTAL ENERGY"] = Esosmp3;

        variables_["OMP3 CORRELATION ENERGY"] = Emp3L - Escf;
        variables_["SCS-OMP3 CORRELATION ENERGY"] = Escsmp3 - Escf;
        variables_["SOS-OMP3 CORRELATION ENERGY"] = Esosmp3 - Escf;

        variables_["CUSTOM SCS-OMP3 CORRELATION ENERGY"] = os_scale * Emp3AB + ss_scale * (Emp3AA + Emp3BB);
        variables_["CUSTOM SCS-OMP3 TOTAL ENERGY"] = Escf + variables_["CUSTOM SCS-OMP3 CORRELATION ENERGY"];

        variables_["CURRENT REFERENCE ENERGY"] = Escf;
        std::map<std::string, double> spin_scale_energies = {
            {"SCS", Escsmp3}, {"SOS", Esosmp3}, {"CUSTOM", variables_["CUSTOM SCS-OMP3 TOTAL ENERGY"]}, {"NONE", Emp3}};
        variables_["CURRENT ENERGY"] = spin_scale_energies[spin_scale_type_];
        variables_["CURRENT CORRELATION ENERGY"] =
            variables_["CURRENT ENERGY"] - variables_["CURRENT REFERENCE ENERGY"];

        if (natorb == "TRUE") nbo();
        if (occ_orb_energy == "TRUE") semi_canonic();

        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
            outfile->Printf("\tAnalytic gradient computation is starting...\n");

            coord_grad();
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }

    }  // end if (conver == 1)
}  // end omp3_manager

//======================================================================
//             MP3 Manager
//======================================================================
void OCCWave::mp3_manager() {
    time4grad = 0;  // means i will not compute the gradient
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    Eref = Escf;
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("MP2 Energy");
    mp2_energy();
    timer_off("MP2 Energy");

    mp2_postprocessing();

    timer_on("T2(2)");
    t2_2nd_sc();
    timer_off("T2(2)");
    timer_on("MP3 Energy");
    mp3_energy();
    timer_off("MP3 Energy");
    Emp3L = Emp3;
    EcorrL = Emp3L - Escf;
    Emp3L_old = Emp3;
    if (ip_poles == "TRUE") omp3_ip_poles();

    mp3_postprocessing();

    variables_["CURRENT REFERENCE ENERGY"] = Escf;
    std::map<std::string, double> spin_scale_energies = {
        {"SCS", Escsmp3}, {"CUSTOM", variables_["CUSTOM SCS-MP3 TOTAL ENERGY"]}, {"NONE", Emp3}};
    variables_["CURRENT ENERGY"] = spin_scale_energies[spin_scale_type_];
    energy_ = variables_["CURRENT ENERGY"];
    variables_["CURRENT CORRELATION ENERGY"] = variables_["CURRENT ENERGY"] - variables_["CURRENT REFERENCE ENERGY"];

    // Compute Analytic Gradients
    if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
        time4grad = 1;
        outfile->Printf("\tAnalytic gradient computation is starting...\n");
        outfile->Printf("\tComputing response density matrices...\n");

        omp3_response_pdms();
        outfile->Printf("\tComputing off-diagonal blocks of GFM...\n");

        gfock();
        outfile->Printf("\tForming independent-pairs...\n");

        idp2();
        outfile->Printf("\tComputing orbital gradient...\n");

        mograd();
        coord_grad();

        if (ekt_ip_ == "TRUE" && ekt_ea_ == "TRUE") {
            ekt_ip();
            ekt_ea();
        }

        else if (ekt_ip_ == "TRUE" && ekt_ea_ == "FALSE") {
            ekt_ip();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "TRUE") {
            ekt_ea();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") {
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }
    }

}  // end mp3_manager

//======================================================================
//             OCEPA Manager
//======================================================================
void OCCWave::ocepa_manager() {
    mo_optimized = 0;  // means MOs are not optimized yet.
    orbs_already_opt = 0;
    orbs_already_sc = 0;
    time4grad = 0;  // means i will not compute the gradient
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    timer_on("REF Energy");
    ref_energy();
    timer_off("REF Energy");
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("MP2 Energy");
    mp2_energy();
    timer_off("MP2 Energy");
    Ecepa = Emp2;
    EcepaL = Ecepa;
    EcorrL = Ecorr;
    EcepaL_old = Ecepa;

    mp2_postprocessing();

    occ_iterations();

    if (conver == 1) {
        ref_energy();
        cepa_energy();
        if (orbs_already_opt == 1) EcepaL = Ecepa;

        // EKT
        if (ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
            if (orbs_already_opt == 1) {
                ocepa_response_pdms();
                gfock();
            }
            gfock_diag();
            if (ekt_ip_ == "TRUE") ekt_ip();
            if (ekt_ea_ == "TRUE") ekt_ea();
        }

        outfile->Printf("\n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\t================ OLCCD [OCEPA(0)] FINAL RESULTS ============================== \n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", EcepaAA);
        outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", EcepaAB);
        outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", EcepaBB);
        outfile->Printf("\tOLCCD Correlation Energy (a.u.)    : %20.14f\n", EcepaL - Escf);
        outfile->Printf("\tEolccd - Eref (a.u.)               : %20.14f\n", EcepaL - Eref);
        outfile->Printf("\tOLCCD Total Energy (a.u.)          : %20.14f\n", EcepaL);
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\n");

        // Set the global variables with the energies
        variables_["OLCCD TOTAL ENERGY"] = EcepaL;
        variables_["CURRENT REFERENCE ENERGY"] = Escf;

        variables_["OLCCD REFERENCE CORRECTION ENERGY"] = Eref - Escf;
        variables_["OLCCD CORRELATION ENERGY"] = EcepaL - Escf;
        variables_["OLCCD SAME-SPIN CORRELATION ENERGY"] = EcepaAA + EcepaBB;
        variables_["OLCCD OPPOSITE-SPIN CORRELATION ENERGY"] = EcepaAB;

        variables_["CUSTOM SCS-OLCCD CORRELATION ENERGY"] = os_scale * EcepaAB + ss_scale * (EcepaAA + EcepaBB);
        variables_["CUSTOM SCS-OLCCD TOTAL ENERGY"] = Escf + variables_["CUSTOM SCS-OLCCD CORRELATION ENERGY"];

        std::map<std::string, double> spin_scale_energies = {{"CUSTOM", variables_["CUSTOM SCS-OLCCD TOTAL ENERGY"]},
                                                             {"NONE", EcepaL}};
        variables_["CURRENT ENERGY"] = spin_scale_energies[spin_scale_type_];
        variables_["CURRENT CORRELATION ENERGY"] =
            variables_["CURRENT ENERGY"] - variables_["CURRENT REFERENCE ENERGY"];
        energy_ = variables_["CURRENT ENERGY"];

        // ordinary ROHF-MP2 not available in course of ROHF-OLCCD
        if (reference == "ROHF") {
            del_scalar_variable("MP2 CORRELATION ENERGY");
            del_scalar_variable("MP2 SINGLES ENERGY");
            del_scalar_variable("MP2 TOTAL ENERGY");
        }

        if (natorb == "TRUE") nbo();
        if (occ_orb_energy == "TRUE") semi_canonic();

        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
            outfile->Printf("\tAnalytic gradient computation is starting...\n");

            coord_grad();
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }

    }  // end if (conver == 1)
}  // end ocepa_manager

//======================================================================
//             CEPA Manager
//======================================================================
void OCCWave::cepa_manager() {
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    Eref = Escf;
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("MP2 Energy");
    mp2_energy();
    timer_off("MP2 Energy");
    Ecepa = Emp2;
    Ecepa_old = Emp2;

    mp2_postprocessing();

    // Perform CEPA iterations
    cepa_iterations();

    outfile->Printf("\n");
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\t================ LCCD [CEPA(0)] FINAL RESULTS ================================ \n");
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", EcepaAA);
    outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", EcepaAB);
    outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", EcepaBB);
    outfile->Printf("\tLCCD Correlation Energy (a.u.)     : %20.14f\n", Ecorr);
    outfile->Printf("\tLCCD Total Energy (a.u.)           : %20.14f\n", Ecepa);
    outfile->Printf("\t============================================================================== \n");
    outfile->Printf("\n");

    // Set the global variables with the energies
    variables_["LCCD TOTAL ENERGY"] = Ecepa;
    variables_["LCCD CORRELATION ENERGY"] = Ecorr;
    variables_["LCCD SAME-SPIN CORRELATION ENERGY"] = EcepaAA + EcepaBB;
    variables_["LCCD OPPOSITE-SPIN CORRELATION ENERGY"] = EcepaAB;
    variables_["LCCD SINGLES ENERGY"] = 0.0;  // CEPA is RHF/UHF only
    variables_["LCCD DOUBLES ENERGY"] = EcepaAB + EcepaAA + EcepaBB;
    variables_["CURRENT REFERENCE ENERGY"] = Eref;
    variables_["CUSTOM SCS-LCCD CORRELATION ENERGY"] = os_scale * EcepaAB + ss_scale * (EcepaAA + EcepaBB);
    variables_["CUSTOM SCS-LCCD TOTAL ENERGY"] = Escf + variables_["CUSTOM SCS-LCCD CORRELATION ENERGY"];
    // EcepaL = Ecepa;

    std::map<std::string, double> spin_scale_energies = {{"CUSTOM", variables_["CUSTOM SCS-LCCD TOTAL ENERGY"]},
                                                         {"NONE", Ecepa}};
    variables_["CURRENT ENERGY"] = spin_scale_energies[spin_scale_type_];
    variables_["CURRENT CORRELATION ENERGY"] = variables_["CURRENT ENERGY"] - variables_["CURRENT REFERENCE ENERGY"];
    energy_ = variables_["CURRENT ENERGY"];

    // Compute Analytic Gradients
    if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
        time4grad = 1;
        outfile->Printf("\tAnalytic gradient computation is starting...\n");
        outfile->Printf("\tComputing response density matrices...\n");

        ocepa_response_pdms();
        outfile->Printf("\tComputing off-diagonal blocks of GFM...\n");

        gfock();
        outfile->Printf("\tForming independent-pairs...\n");

        idp2();
        outfile->Printf("\tComputing orbital gradient...\n");

        mograd();
        coord_grad();

        if (ekt_ip_ == "TRUE" && ekt_ea_ == "TRUE") {
            ekt_ip();
            ekt_ea();
        }

        else if (ekt_ip_ == "TRUE" && ekt_ea_ == "FALSE") {
            ekt_ip();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "TRUE") {
            ekt_ea();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") {
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }
    }
}  // end cepa_manager

//======================================================================
//             OMP2.5 Manager
//======================================================================
void OCCWave::omp2_5_manager() {
    mo_optimized = 0;
    orbs_already_opt = 0;
    orbs_already_sc = 0;
    time4grad = 0;  // means i will not compute the gradient
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    timer_on("REF Energy");
    ref_energy();
    timer_off("REF Energy");
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("MP2 Energy");
    mp2_energy();
    timer_off("MP2 Energy");

    mp2_postprocessing();

    timer_on("T2(2)");
    t2_2nd_sc();
    timer_off("T2(2)");
    timer_on("MP3 Energy");
    mp3_energy();
    timer_off("MP3 Energy");
    Emp3L = Emp3;
    EcorrL = Emp3L - Escf;
    Emp3L_old = Emp3;
    if (ip_poles == "TRUE") omp3_ip_poles();

    mp2p5_postprocessing();

    occ_iterations();

    if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod) {
        orbs_already_opt = 1;
        if (conver == 1)
            outfile->Printf("\n\tOrbitals are optimized now.\n");
        else if (conver == 0) {
            outfile->Printf("\n\tMAX MOGRAD did NOT converge, but RMS MOGRAD converged!!!\n");
            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
        }
        outfile->Printf("\tSwitching to the standard MP2.5 computation after semicanonicalization of the MOs... \n");

        semi_canonic();
        if (reference_ == "RESTRICTED")
            trans_ints_rhf();
        else if (reference_ == "UNRESTRICTED")
            trans_ints_uhf();
        set_t2_amplitudes_mp2();
        t2_2nd_sc();
        conver = 1;
        if (dertype == "FIRST") {
            omp3_response_pdms();
            gfock();
        }
    }

    if (conver == 1) {
        ref_energy();
        mp2_energy();
        mp3_energy();
        if (orbs_already_opt == 1) Emp3L = Emp3;

        if (ip_poles == "TRUE") {
            if (orbs_already_sc == 0) {
                semi_canonic();
                if (reference_ == "RESTRICTED")
                    trans_ints_rhf();
                else if (reference_ == "UNRESTRICTED")
                    trans_ints_uhf();
                set_t2_amplitudes_mp2();
                t2_2nd_sc();
            }
            omp3_ip_poles();
        }

        // EKT
        if (ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
            if (orbs_already_sc == 1) {
                omp3_response_pdms();
                gfock();
            }
            gfock_diag();
            if (ekt_ip_ == "TRUE") ekt_ip();
            if (ekt_ea_ == "TRUE") ekt_ea();
        }

        mp2_printing();
        mp2p5_printing();

        outfile->Printf("\n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\t================ OMP2.5 FINAL RESULTS ======================================== \n");
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tOMP2.5 Correlation Energy (a.u.)   : %20.14f\n", Emp3L - Escf);
        outfile->Printf("\tEomp2.5 - Eref (a.u.)              : %20.14f\n", Emp3L - Eref);
        outfile->Printf("\tOMP2.5 Total Energy (a.u.)         : %20.14f\n", Emp3L);
        outfile->Printf("\t============================================================================== \n");
        outfile->Printf("\n");

        // Set the global variables with the energies
        variables_["OMP2.5 TOTAL ENERGY"] = Emp3L;
        variables_["OMP2.5 CORRELATION ENERGY"] = Emp3L - Escf;
        variables_["CURRENT REFERENCE ENERGY"] = Escf;

        std::map<std::string, double> spin_scale_energies = {{"CUSTOM", variables_["CUSTOM SCS-OMP2.5 TOTAL ENERGY"]},
                                                             {"NONE", Emp3L}};
        variables_["CURRENT ENERGY"] = spin_scale_energies[spin_scale_type_];
        variables_["CURRENT CORRELATION ENERGY"] =
            variables_["CURRENT ENERGY"] - variables_["CURRENT REFERENCE ENERGY"];

        if (natorb == "TRUE") nbo();
        if (occ_orb_energy == "TRUE") semi_canonic();

        // Compute Analytic Gradients
        if (dertype == "FIRST") {
            time4grad = 1;
            outfile->Printf("\tAnalytic gradient computation is starting...\n");

            coord_grad();
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }

    }  // end if (conver == 1)
}  // end omp2.5_manager

//======================================================================
//             MP2.5 Manager
//======================================================================
void OCCWave::mp2_5_manager() {
    time4grad = 0;  // means i will not compute the gradient
    timer_on("trans_ints");
    if (reference_ == "RESTRICTED")
        trans_ints_rhf();
    else if (reference_ == "UNRESTRICTED")
        trans_ints_uhf();
    timer_off("trans_ints");
    Eref = Escf;
    timer_on("T2(1)");
    set_t2_amplitudes_mp2();
    timer_off("T2(1)");
    timer_on("MP2 Energy");
    mp2_energy();
    timer_off("MP2 Energy");

    mp2_postprocessing();

    timer_on("T2(2)");
    t2_2nd_sc();
    timer_off("T2(2)");
    timer_on("MP3 Energy");
    mp3_energy();
    timer_off("MP3 Energy");
    Emp3L = Emp3;
    EcorrL = Emp3L - Escf;
    Emp3L_old = Emp3;
    if (ip_poles == "TRUE") omp3_ip_poles();

    mp2p5_postprocessing();
    variables_["CURRENT REFERENCE ENERGY"] = Eref;
    variables_["CUSTOM SCS-MP2.5 CORRELATION ENERGY"] = os_scale * Emp3AB + ss_scale * (Emp3AA + Emp3BB);
    variables_["CUSTOM SCS-MP2.5 TOTAL ENERGY"] = Escf + variables_["CUSTOM SCS-MP2.5 CORRELATION ENERGY"];
    std::map<std::string, double> spin_scale_energies = {{"CUSTOM", variables_["CUSTOM SCS-MP2.5 TOTAL ENERGY"]},
                                                         {"NONE", Emp3L}};
    variables_["CURRENT ENERGY"] = spin_scale_energies[spin_scale_type_];
    variables_["CURRENT CORRELATION ENERGY"] = variables_["CURRENT ENERGY"] - variables_["CURRENT REFERENCE ENERGY"];
    energy_ = variables_["CURRENT ENERGY"];

    // Compute Analytic Gradients
    if (dertype == "FIRST" || ekt_ip_ == "TRUE" || ekt_ea_ == "TRUE") {
        time4grad = 1;
        outfile->Printf("\tAnalytic gradient computation is starting...\n");
        outfile->Printf("\tComputing response density matrices...\n");

        omp3_response_pdms();
        outfile->Printf("\tComputing off-diagonal blocks of GFM...\n");

        gfock();
        outfile->Printf("\tForming independent-pairs...\n");

        idp2();
        outfile->Printf("\tComputing orbital gradient...\n");

        mograd();
        coord_grad();

        if (ekt_ip_ == "TRUE" && ekt_ea_ == "TRUE") {
            ekt_ip();
            ekt_ea();
        }

        else if (ekt_ip_ == "TRUE" && ekt_ea_ == "FALSE") {
            ekt_ip();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "TRUE") {
            ekt_ea();
        }

        else if (ekt_ip_ == "FALSE" && ekt_ea_ == "FALSE") {
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");
        }
    }

}  // end omp2.5_manager
}
}  // End Namespaces
