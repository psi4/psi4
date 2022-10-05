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

#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/process.h"
#include "dfocc.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

//======================================================================
//             CD-OMP2 Manager
//======================================================================
void DFOCC::cd_omp2_manager() {
    do_cd = "TRUE";
    time4grad = 0;         // means i will not compute the gradient
    mo_optimized = 0;      // means MOs are not optimized
    orbs_already_opt = 0;  // means orbitals are not optimized yet.
    orbs_already_sc = 0;   // menas orbitals are not semicanonical yet.
    timer_on("CD Integrals");
    cd_ints();
    timer_off("CD Integrals");
    timer_on("CD Trans");
    trans_cd();
    timer_off("CD Trans");

    // memalloc for density intermediates
    Jc = std::make_shared<Tensor1d>("DF_BASIS_SCF J_Q", nQ_ref);
    g1Qc = std::make_shared<Tensor1d>("DF_BASIS_SCF G1_Q", nQ_ref);
    g1Qt = std::make_shared<Tensor1d>("DF_BASIS_SCF G1t_Q", nQ_ref);
    g1Q = std::make_shared<Tensor1d>("DF_BASIS_CC G1_Q", nQ);
    g1Qt2 = std::make_shared<Tensor1d>("DF_BASIS_CC G1t_Q", nQ);

    // Fock
    fock();

    // ROHF REF
    if (reference == "ROHF") t1_1st_sc();
    t2_1st_sc();
    Emp2L = Emp2;
    EcorrL = Emp2L - Escf;
    Emp2L_old = Emp2;

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (Canonical CD-MP2)... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    energy_ = Emp2;
    variables_["CURRENT ENERGY"] = Emp2;
    variables_["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;

    variables_["CURRENT REFERENCE ENERGY"] = Escf;
    variables_["CURRENT CORRELATION ENERGY"] = Emp2 - Escf;
    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

    if (reference_ == "UNRESTRICTED") {
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        variables_["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    variables_["MP2 SINGLES ENERGY"] = Emp2_t1;

    omp2_opdm();
    omp2_tpdm();
    mp2l_energy();
    separable_tpdm();
    gfock_vo();
    gfock_ov();
    gfock_oo();
    gfock_vv();
    idp();
    mograd();
    occ_iterations();

    if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod) {
        orbs_already_opt = 1;
        if (conver == 1)
            outfile->Printf("\n\tOrbitals are optimized now.\n");
        else if (conver == 0) {
            outfile->Printf("\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
        }
        outfile->Printf("\tSwitching to the standard CD-MP2 computation after semicanonicalization of the MOs... \n");

        semi_canonic();
        trans_cd();
        fock();
        t2_1st_sc();
        conver = 1;
        if (dertype == "FIRST") {
            omp2_opdm();
            omp2_tpdm();
            separable_tpdm();
            gfock_vo();
            gfock_ov();
            gfock_oo();
            gfock_vv();
        }
    }

    if (conver == 1) {
        ref_energy();
        mp2_energy();
        if (orbs_already_opt == 1) Emp2L = Emp2;

        outfile->Printf("\n");
        outfile->Printf("\tComputing MP2 energy using optimized MOs... \n");
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
        if (reference_ == "UNRESTRICTED")
            outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
        if (reference_ == "UNRESTRICTED")
            outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
        outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Emp2 - Escf);
        outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
        outfile->Printf("\t======================================================================= \n");

        outfile->Printf("\n");
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\t================ CD-OMP2 FINAL RESULTS ================================ \n");
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-OMP2 Total Energy (a.u.)    : %20.14f\n", Escsmp2);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-OMP2 Total Energy (a.u.)    : %20.14f\n", Esosmp2);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS(N)-OMP2 Total Energy (a.u.) : %20.14f\n", Escsnmp2);
        outfile->Printf("\tCD-OMP2 Correlation Energy (a.u.)  : %20.14f\n", Emp2L - Escf);
        outfile->Printf("\tEcdomp2 - Eref (a.u.)              : %20.14f\n", Emp2L - Eref);
        outfile->Printf("\tCD-OMP2 Total Energy (a.u.)        : %20.14f\n", Emp2L);
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\n");

        // Set the global variables with the energies
        variables_["OMP2 TOTAL ENERGY"] = Emp2L;
        Process::environment.globals["SCS-OMP2 TOTAL ENERGY"] = Escsmp2;
        Process::environment.globals["SOS-OMP2 TOTAL ENERGY"] = Esosmp2;
        Process::environment.globals["SCS(N)-OMP2 TOTAL ENERGY"] = Escsnmp2;
        energy_ = Emp2L;
        variables_["CURRENT ENERGY"] = Emp2L;
        variables_["CURRENT REFERENCE ENERGY"] = Escf;
        variables_["CURRENT CORRELATION ENERGY"] = Emp2L - Escf;

        variables_["OMP2 CORRELATION ENERGY"] = Emp2L - Escf;
        variables_["OMP2 REFERENCE CORRECTION ENERGY"] = Eref - Escf;
        Process::environment.globals["SCS-OMP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
        Process::environment.globals["SOS-OMP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
        Process::environment.globals["SCS(N)-OMP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

        if (reference_ == "UNRESTRICTED") {
            variables_["OMP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
            variables_["OMP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        }

        // if scs on
        if (do_scs == "TRUE") {
            if (scs_type_ == "SCS") {
                energy_ = Escsmp2;
                variables_["CURRENT ENERGY"] = Escsmp2;
                variables_["CURRENT CORRELATION ENERGY"] = Escsmp2 - Escf;
            }

            else if (scs_type_ == "SCS(N)") {
                energy_ = Escsnmp2;
                variables_["CURRENT ENERGY"] = Escsnmp2;
                variables_["CURRENT CORRELATION ENERGY"] = Escsnmp2 - Escf;
            }
        }

        // else if sos on
        else if (do_sos == "TRUE") {
            if (sos_type_ == "SOS") {
                energy_ = Esosmp2;
                variables_["CURRENT ENERGY"] = Esosmp2;
                variables_["CURRENT CORRELATION ENERGY"] = Esosmp2 - Escf;
            }
        }

        // OEPROP
        if (oeprop_ == "TRUE") oeprop();

        // if (natorb == "TRUE") nbo();
        // if (occ_orb_energy == "TRUE") semi_canonic();

        // Compute Analytic Gradients
        /*
        if (dertype == "FIRST") {
            outfile->Printf("\tAnalytic gradient computation is starting...\n");

            coord_grad();
            outfile->Printf("\tNecessary information has been sent to DERIV, which will take care of the rest.\n");

        }
        */

    }  // end if (conver == 1)
}  // end omp2_manager

//======================================================================
//             CD-MP2 Manager
//======================================================================
void DFOCC::cd_mp2_manager() {
    do_cd = "TRUE";
    time4grad = 0;     // means i will not compute the gradient
    mo_optimized = 0;  // means MOs are not optimized
    timer_on("CD Integrals");
    cd_ints();
    if (dertype == "NONE" && ekt_ip_ == "FALSE") {
        trans_cd_mp2();
    } else
        trans_cd();
    timer_off("CD Integrals");

    // ROHF REF
    // outfile->Printf("\tI am here.\n");
    if (reference == "ROHF") t1_1st_sc();

    // t2_1st_sc();
    // mp2_energy();
    if (dertype == "NONE" && ekt_ip_ == "FALSE")
        mp2_direct();
    else {
        t2_1st_sc();
        mp2_energy();
    }
    Emp2L = Emp2;
    EcorrL = Emp2L - Escf;

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (Canonical CD-MP2)... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    energy_ = Emp2;
    variables_["CURRENT ENERGY"] = Emp2;
    variables_["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;

    variables_["CURRENT REFERENCE ENERGY"] = Escf;
    variables_["CURRENT CORRELATION ENERGY"] = Emp2 - Escf;
    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

    if (reference_ == "UNRESTRICTED") {
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        variables_["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    variables_["MP2 SINGLES ENERGY"] = Emp2_t1;

}  // end mp2_manager

//======================================================================
//             CCSD Manager
//======================================================================
void DFOCC::ccsd_manager_cd() {
    do_cd = "TRUE";
    time4grad = 0;     // means i will not compute the gradient
    mo_optimized = 0;  // means MOs are not optimized

    timer_on("CD Integrals");
    cd_ints();
    trans_cd();
    timer_off("CD Integrals");

    // Memory allocation
    T1c = std::make_shared<Tensor1d>("DF_BASIS_CC T1_Q", nQ);
    Jc = std::make_shared<Tensor1d>("DF_BASIS_SCF J_Q", nQ_ref);

    if (reference_ == "RESTRICTED") {
        t1A = std::make_shared<Tensor2d>("T1 <I|A>", naoccA, navirA);
        t1newA = std::make_shared<Tensor2d>("New T1 <I|A>", naoccA, navirA);
        FiaA = std::make_shared<Tensor2d>("Fint <I|A>", naoccA, navirA);
        FtijA = std::make_shared<Tensor2d>("Ftilde <I|J>", naoccA, naoccA);
        FtabA = std::make_shared<Tensor2d>("Ftilde <A|B>", navirA, navirA);

        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory / (1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements

        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = (long long int)nQ * (long long int)nso2_;
        cost_ampAA += (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_ampAA += (long long int)nQ * (long long int)nso_ * (long long int)navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;
        /*
        if ((cost_5amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_5amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_5amp+cost_df);
             nincore_amp = 5;
             t2_incore = true;
             df_ints_incore = true;
        }
        */
        if ((cost_4amp + cost_df) <= memory_mb) {
            outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
            outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp + cost_df);
            nincore_amp = 4;
            t2_incore = true;
            df_ints_incore = true;
        } else if ((cost_3amp + cost_df) <= memory_mb) {
            outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
            // outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
            outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
            nincore_amp = 3;
            t2_incore = false;
            df_ints_incore = false;
        } else if (cost_3amp < memory_mb && cost_df < memory_mb) {
            outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
            outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
            nincore_amp = 3;
            t2_incore = false;
            df_ints_incore = false;
        } else {
            outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
            outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp + cost_df - memory_mb);
            throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        double cost_amp1 = 0.0;
        cost_amp1 = 2.5 * (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_amp1 += (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_amp1 /= 1024.0 * 1024.0;
        cost_amp1 *= sizeof(double);
        double cost_amp2 = 0.0;
        cost_amp2 = 1.5 * (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_amp2 += 4.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_amp2 /= 1024.0 * 1024.0;
        cost_amp2 *= sizeof(double);
        double cost_amp3 = 0.0;
        cost_amp3 = 2.0 * (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_amp3 += 3.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_amp3 += 2.0 * (long long int)navirA * (long long int)navirA * (long long int)navirA;
        cost_amp3 /= 1024.0 * 1024.0;
        cost_amp3 *= sizeof(double);
        cost_amp = MAX0(cost_amp1, cost_amp2);
        cost_amp = MAX0(cost_amp, cost_amp3);
        outfile->Printf("\tMemory requirement for Wabef term (T2): %9.2lf MB \n", cost_amp);
        if (cc_lambda_ == "TRUE") {
            cost_amp1 = (long long int)navirA * (long long int)navirA * (long long int)navirA;
            cost_amp1 /= 1024.0 * 1024.0;
            cost_amp += cost_amp1;
            outfile->Printf("\tMemory requirement for Wefab term (L2): %9.2lf MB \n", cost_amp);
        }

        // cost_ppl_hm
        cost_ppl_hm = ntri_abAA / 1024.0;
        cost_ppl_hm *= cost_ppl_hm;
        cost_ppl_hm *= sizeof(double);
        cost_ppl_hm += cost_amp;
        outfile->Printf("\tMemory for high mem Wabef algorithm   : %9.2lf MB \n", cost_ppl_hm);
        if (cost_ppl_hm > memory_mb && Wabef_type_ == "AUTO") {
            do_ppl_hm = false;
            outfile->Printf("\tI will use the LOW_MEM Wabef algorithm! \n");
        } else if (cost_ppl_hm <= memory_mb && Wabef_type_ == "AUTO") {
            do_ppl_hm = true;
            outfile->Printf("\tI will use the HIGH_MEM Wabef algorithm! \n");
        }

        // Mem alloc for DF ints
        if (df_ints_incore) {
            bQijA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
            bQiaA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
            bQabA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
            bQijA->read(psio_, PSIF_DFOCC_INTS);
            bQiaA->read(psio_, PSIF_DFOCC_INTS);
            bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
        }

        //  Malloc
        if (t2_incore) {
            t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        }

    }  // end if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
        t1A = std::make_shared<Tensor2d>("T1 <I|A>", naoccA, navirA);
        t1B = std::make_shared<Tensor2d>("T1 <i|a>", naoccB, navirB);
        FiaA = std::make_shared<Tensor2d>("Fint <I|A>", naoccA, navirA);
        FiaB = std::make_shared<Tensor2d>("Fint <i|a>", naoccB, navirB);
        FtijA = std::make_shared<Tensor2d>("Ftilde <I|J>", naoccA, naoccA);
        FtabA = std::make_shared<Tensor2d>("Ftilde <A|B>", navirA, navirA);
        FtijB = std::make_shared<Tensor2d>("Ftilde <i|j>", naoccB, naoccB);
        FtabB = std::make_shared<Tensor2d>("Ftilde <a|b>", navirB, navirB);

        // memory requirements
        cost_ampAA = 0.0;
        cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_ampBB = (long long int)naocc2BB * (long long int)nvir2BB;
        cost_ampBB /= 1024.0 * 1024.0;
        cost_ampBB *= sizeof(double);
        cost_ampAB = (long long int)naocc2AB * (long long int)nvir2AB;
        cost_ampAB /= 1024.0 * 1024.0;
        cost_ampAB *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampBB);
        cost_amp = MAX0(cost_amp, cost_ampAB);
        cost_amp = 3.0 * cost_amp;
        memory = Process::environment.get_memory();
        memory_mb = (double)memory / (1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);
        outfile->Printf("\tMinimum required memory for amplitudes: %9.2lf MB \n", cost_amp);
    }  // else if (reference_ == "UNRESTRICTED")

    // memalloc for density intermediates
    if (qchf_ == "TRUE" || dertype == "FIRST") {
        g1Qc = std::make_shared<Tensor1d>("DF_BASIS_SCF G1_Q", nQ_ref);
        g1Qt = std::make_shared<Tensor1d>("DF_BASIS_SCF G1t_Q", nQ_ref);
        g1Qp = std::make_shared<Tensor1d>("DF_BASIS_SCF G1p_Q", nQ_ref);
        g1Q = std::make_shared<Tensor1d>("DF_BASIS_CC G1_Q", nQ);
        g1Qt2 = std::make_shared<Tensor1d>("DF_BASIS_CC G1t_Q", nQ);
    }

    // QCHF
    if (qchf_ == "TRUE") qchf();

    // Fock
    if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE") fock();

    // Compute MP2 energy
    if (reference == "ROHF") t1_1st_sc();
    if (t2_incore)
        ccsd_mp2();
    else
        ccsd_mp2_low();

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing DF-MP2 energy ... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;
    Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
    variables_["MP2 TOTAL ENERGY"] = Emp2;
    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        Process::environment.globals["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        Process::environment.globals["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
        variables_["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    Process::environment.globals["MP2 SINGLES ENERGY"] = Emp2_t1;
    variables_["MP2 SINGLES ENERGY"] = Emp2_t1;

    // Perform CCSD iterations
    timer_on("CCSD");
    if (t2_incore)
        ccsd_iterations();
    else
        ccsd_iterations_low();
    timer_off("CCSD");

    outfile->Printf("\n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\t================ CCSD FINAL RESULTS =================================== \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    outfile->Printf("\tCD-CCSD Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-CCSD Total Energy (a.u.)        : %20.14f\n", Eccsd);
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\n");

    Process::environment.globals["CURRENT ENERGY"] = Eccsd;
    Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
    Process::environment.globals["CURRENT CORRELATION ENERGY"] = Eccsd - Escf;
    Process::environment.globals["CCSD TOTAL ENERGY"] = Eccsd;
    Process::environment.globals["CCSD CORRELATION ENERGY"] = Eccsd - Escf;
    variables_["CURRENT CORRELATION ENERGY"] = Eccsd - Escf;
    variables_["CCSD TOTAL ENERGY"] = Eccsd;
    variables_["CCSD CORRELATION ENERGY"] = Eccsd - Escf;
    if ((reference == "RHF") or (reference == "UHF")) {
        Process::environment.globals["CCSD DOUBLES ENERGY"] = Eccsd - Escf;
        Process::environment.globals["CCSD SINGLES ENERGY"] = 0.0;
        variables_["CCSD DOUBLES ENERGY"] = Eccsd - Escf;
        variables_["CCSD SINGLES ENERGY"] = 0.0;
    }

    /* updates the wavefunction for checkpointing */
    set_energy(Eccsd);

    // CCSDL
    if (dertype == "FIRST" || cc_lambda_ == "TRUE") {
        // memalloc
        if (dertype == "FIRST") {
            GtijA = std::make_shared<Tensor2d>("Gtilde Intermediate <I|J>", naoccA, naoccA);
            GtabA = std::make_shared<Tensor2d>("Gtilde Intermediate <A|B>", navirA, navirA);
            L1c = std::make_shared<Tensor1d>("DF_BASIS_CC L1_Q", nQ);
            gQt = std::make_shared<Tensor1d>("CCSD PDM G_Qt", nQ);
        }

        timer_on("CCSDL");
        if (t2_incore) {
            tstop();
            tstart();
            lambda_title();
            outfile->Printf("\tSolving Lambda amplitude equations...\n");
            ccsdl_iterations();
        } else
            throw PSIEXCEPTION("There is NOT enough memory for Lambda equations!");
        timer_off("CCSDL");
    }

    // Compute Analytic Gradients
    if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
        tstop();
        tstart();
        pdm_title();

        // memalloc
        if (reference_ == "RESTRICTED") {
            G1c_ov = std::make_shared<Tensor2d>("Correlation OPDM <O|V>", noccA, nvirA);
            G1c_vo = std::make_shared<Tensor2d>("Correlation OPDM <V|O>", nvirA, noccA);

            outfile->Printf("\tComputing unrelaxed response density matrices...\n");
            ccsd_opdm();
            ccsd_tpdm();
        }
        ccl_energy();
        // prepare4grad();
        // if (oeprop_ == "TRUE") oeprop();
        // if (dertype == "FIRST") dfgrad();
        // if (ekt_ip_ == "TRUE") ekt_ip();
    }  // if (dertype == "FIRST" || ekt_ip_ == "TRUE")

}  // end ccsd_manager_cd

//======================================================================
//             CCSD(T) Manager
//======================================================================
void DFOCC::ccsd_t_manager_cd() {
    do_cd = "TRUE";
    time4grad = 0;     // means i will not compute the gradient
    mo_optimized = 0;  // means MOs are not optimized

    timer_on("CD Integrals");
    cd_ints();
    trans_cd();
    timer_off("CD Integrals");

    // Memory allocation
    T1c = std::make_shared<Tensor1d>("DF_BASIS_CC T1_Q", nQ);
    Jc = std::make_shared<Tensor1d>("DF_BASIS_SCF J_Q", nQ_ref);

    if (reference_ == "RESTRICTED") {
        t1A = std::make_shared<Tensor2d>("T1 <I|A>", naoccA, navirA);
        t1newA = std::make_shared<Tensor2d>("New T1 <I|A>", naoccA, navirA);
        FiaA = std::make_shared<Tensor2d>("Fint <I|A>", naoccA, navirA);
        FtijA = std::make_shared<Tensor2d>("Ftilde <I|J>", naoccA, naoccA);
        FtabA = std::make_shared<Tensor2d>("Ftilde <A|B>", navirA, navirA);

        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory / (1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements

        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = (long long int)nQ * (long long int)nso2_;
        cost_ampAA += (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_ampAA += (long long int)nQ * (long long int)nso_ * (long long int)navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;

        if ((cost_4amp + cost_df) <= memory_mb) {
            outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
            outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp + cost_df);
            nincore_amp = 4;
            t2_incore = true;
            df_ints_incore = true;
        } else if ((cost_3amp + cost_df) <= memory_mb) {
            outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
            // outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
            outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
            nincore_amp = 3;
            t2_incore = false;
            df_ints_incore = false;
        } else if (cost_3amp < memory_mb && cost_df < memory_mb) {
            outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
            outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
            nincore_amp = 3;
            t2_incore = false;
            df_ints_incore = false;
        } else {
            outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
            outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp + cost_df - memory_mb);
            throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        double cost_amp1 = 0.0;
        cost_amp1 = 2.5 * (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_amp1 += (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_amp1 /= 1024.0 * 1024.0;
        cost_amp1 *= sizeof(double);
        double cost_amp2 = 0.0;
        cost_amp2 = 1.5 * (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_amp2 += 4.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_amp2 /= 1024.0 * 1024.0;
        cost_amp2 *= sizeof(double);
        double cost_amp3 = 0.0;
        cost_amp3 = 2.0 * (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_amp3 += 3.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_amp3 += 2.0 * (long long int)navirA * (long long int)navirA * (long long int)navirA;
        cost_amp3 /= 1024.0 * 1024.0;
        cost_amp3 *= sizeof(double);
        cost_amp = MAX0(cost_amp1, cost_amp2);
        cost_amp = MAX0(cost_amp, cost_amp3);
        outfile->Printf("\tMemory requirement for Wabef term (T2): %9.2lf MB \n", cost_amp);

        if (cc_lambda_ == "TRUE") {
            cost_amp1 = (long long int)navirA * (long long int)navirA * (long long int)navirA;
            cost_amp1 /= 1024.0 * 1024.0;
            cost_amp += cost_amp1;
            outfile->Printf("\tMemory requirement for Wefab term (L2): %9.2lf MB \n", cost_amp);
        }

        // cost_ppl_hm
        cost_ppl_hm = ntri_abAA / 1024.0;
        cost_ppl_hm *= cost_ppl_hm;
        cost_ppl_hm *= sizeof(double);
        cost_ppl_hm += cost_amp;
        outfile->Printf("\tMemory for high mem Wabef algorithm   : %9.2lf MB \n", cost_ppl_hm);
        if (cost_ppl_hm > memory_mb && Wabef_type_ == "AUTO") {
            do_ppl_hm = false;
            outfile->Printf("\tI will use the LOW_MEM Wabef algorithm! \n");
        } else if (cost_ppl_hm <= memory_mb && Wabef_type_ == "AUTO") {
            do_ppl_hm = true;
            outfile->Printf("\tI will use the HIGH_MEM Wabef algorithm! \n");
        }

        // Memory for triples: 2*O^2V^2 + 5*V^3 + O^3V + V^2N + V^3/2
        cost_amp1 = 0.0;
        cost_amp1 = 2.0 * (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_amp1 += 5.0 * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_amp1 += (long long int)naoccA * (long long int)naoccA * (long long int)naoccA * (long long int)navirA;
        cost_amp1 += (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_amp1 += (long long int)navirA * (long long int)ntri_abAA;
        cost_amp1 /= 1024.0 * 1024.0;
        cost_amp1 *= sizeof(double);
        // Memory: OV^3 + 2*O^2V^2 + 2*V^3 + O^3V + V^2N
        cost_triples_iabc = 0.0;
        cost_triples_iabc = 2.0 * (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_triples_iabc += 5.0 * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_triples_iabc += (long long int)naoccA * (long long int)naoccA * (long long int)naoccA * (long long int)navirA;
        cost_triples_iabc += (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_triples_iabc /= 1024.0 * 1024.0;
        cost_amp2 = 0.0;
        cost_amp2 = (navirA * navirA) / 1024.0;
        cost_amp2 *= (naoccA * navirA) / 1024.0;
        cost_triples_iabc += cost_amp2;
        cost_triples_iabc *= sizeof(double);

        if (triples_iabc_type_ == "DISK") {
            do_triples_hm = false;
            // outfile->Printf("\n\tI will use a DISK algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_amp1);
        } else if (triples_iabc_type_ == "DIRECT") {
            do_triples_hm = false;
            outfile->Printf("\n\tI will use a DIRECT algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_amp1);
        } else if (cost_triples_iabc > memory_mb && triples_iabc_type_ == "AUTO") {
            do_triples_hm = false;
            outfile->Printf("\n\tI will use a DIRECT algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_amp1);
        } else if (cost_triples_iabc <= memory_mb && triples_iabc_type_ == "AUTO") {
            do_triples_hm = true;
            outfile->Printf("\n\tI will use an INCORE algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_triples_iabc);
        }

        // Mem alloc for DF ints
        if (df_ints_incore) {
            bQijA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
            bQiaA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
            bQabA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
            bQijA->read(psio_, PSIF_DFOCC_INTS);
            bQiaA->read(psio_, PSIF_DFOCC_INTS);
            bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
        }

        //  Malloc
        if (t2_incore) {
            t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        }

    }  // end if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
        t1A = std::make_shared<Tensor2d>("T1 <I|A>", naoccA, navirA);
        t1B = std::make_shared<Tensor2d>("T1 <i|a>", naoccB, navirB);
        FiaA = std::make_shared<Tensor2d>("Fint <I|A>", naoccA, navirA);
        FiaB = std::make_shared<Tensor2d>("Fint <i|a>", naoccB, navirB);
        FtijA = std::make_shared<Tensor2d>("Ftilde <I|J>", naoccA, naoccA);
        FtabA = std::make_shared<Tensor2d>("Ftilde <A|B>", navirA, navirA);
        FtijB = std::make_shared<Tensor2d>("Ftilde <i|j>", naoccB, naoccB);
        FtabB = std::make_shared<Tensor2d>("Ftilde <a|b>", navirB, navirB);

        // memory requirements
        cost_ampAA = 0.0;
        cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_ampBB = (long long int)naocc2BB * (long long int)nvir2BB;
        cost_ampBB /= 1024.0 * 1024.0;
        cost_ampBB *= sizeof(double);
        cost_ampAB = (long long int)naocc2AB * (long long int)nvir2AB;
        cost_ampAB /= 1024.0 * 1024.0;
        cost_ampAB *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampBB);
        cost_amp = MAX0(cost_amp, cost_ampAB);
        cost_amp = 3.0 * cost_amp;
        memory = Process::environment.get_memory();
        memory_mb = (double)memory / (1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);
        outfile->Printf("\tMinimum required memory for amplitudes: %9.2lf MB \n", cost_amp);
    }  // else if (reference_ == "UNRESTRICTED")

    // memalloc for density intermediates
    if (qchf_ == "TRUE" || dertype == "FIRST") {
        g1Qc = std::make_shared<Tensor1d>("DF_BASIS_SCF G1_Q", nQ_ref);
        g1Qt = std::make_shared<Tensor1d>("DF_BASIS_SCF G1t_Q", nQ_ref);
        g1Qp = std::make_shared<Tensor1d>("DF_BASIS_SCF G1p_Q", nQ_ref);
        g1Q = std::make_shared<Tensor1d>("DF_BASIS_CC G1_Q", nQ);
        g1Qt2 = std::make_shared<Tensor1d>("DF_BASIS_CC G1t_Q", nQ);
    }

    // QCHF
    if (qchf_ == "TRUE") qchf();

    // Fock
    if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE") fock();

    // Compute MP2 energy
    if (reference == "ROHF") t1_1st_sc();
    if (t2_incore)
        ccsd_mp2();
    else
        ccsd_mp2_low();

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing CD-MP2 energy ... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;
    Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
    variables_["MP2 TOTAL ENERGY"] = Emp2;
    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        Process::environment.globals["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        Process::environment.globals["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
        variables_["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    Process::environment.globals["MP2 SINGLES ENERGY"] = Emp2_t1;
    variables_["MP2 SINGLES ENERGY"] = Emp2_t1;

    // Perform CCSD iterations
    timer_on("CCSD");
    if (t2_incore)
        ccsd_iterations();
    else
        ccsd_iterations_low();
    timer_off("CCSD");

    outfile->Printf("\n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\t================ CCSD FINAL RESULTS =================================== \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    outfile->Printf("\tCD-CCSD Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-CCSD Total Energy (a.u.)        : %20.14f\n", Eccsd);
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\n");
    Process::environment.globals["CCSD TOTAL ENERGY"] = Eccsd;
    Process::environment.globals["CCSD CORRELATION ENERGY"] = Eccsd - Escf;
    variables_["CCSD TOTAL ENERGY"] = Eccsd;
    variables_["CCSD CORRELATION ENERGY"] = Eccsd - Escf;
    if ((reference == "RHF") or (reference == "UHF")) {
        Process::environment.globals["CCSD DOUBLES ENERGY"] = Eccsd - Escf;
        Process::environment.globals["CCSD SINGLES ENERGY"] = 0.0;
        variables_["CCSD DOUBLES ENERGY"] = Eccsd - Escf;
        variables_["CCSD SINGLES ENERGY"] = 0.0;
    }

    // CCSD(T)
    tstop();
    tstart();
    pt_title();
    outfile->Printf("\tComputing (T) correction...\n");
    timer_on("(T)");
    if (dertype == "FIRST") {
        ccsd_canonic_triples_grad();
    } else {
        if (triples_iabc_type_ == "DISK")
            ccsd_canonic_triples_disk();
        else if (triples_iabc_type_ == "AUTO") {
            if (do_triples_hm)
                ccsd_canonic_triples_hm();
            else
                ccsd_canonic_triples();
        } else if (triples_iabc_type_ == "INCORE")
            ccsd_canonic_triples_hm();
        else if (triples_iabc_type_ == "DIRECT")
            ccsd_canonic_triples();
        else if (triples_iabc_type_ == "DISK")
            ccsd_canonic_triples_disk();
    }
    timer_off("(T)");
    outfile->Printf("\t(T) Correction (a.u.)              : %20.14f\n", E_t);
    outfile->Printf("\tCD-CCSD(T) Total Energy (a.u.)     : %20.14f\n", Eccsd_t);
    if (dertype == "FIRST") {
        tstop();
        tstart();
    }

    Process::environment.globals["CURRENT ENERGY"] = Eccsd_t;
    Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
    Process::environment.globals["CURRENT CORRELATION ENERGY"] = Eccsd_t - Escf;
    Process::environment.globals["CCSD(T) CORRELATION ENERGY"] = Eccsd_t - Escf;
    Process::environment.globals["CCSD(T) TOTAL ENERGY"] = Eccsd_t;
    Process::environment.globals["(T) CORRECTION ENERGY"] = E_t;
    variables_["CURRENT ENERGY"] = Eccsd_t;
    variables_["CURRENT REFERENCE ENERGY"] = Escf;
    variables_["CURRENT CORRELATION ENERGY"] = Eccsd_t - Escf;
    variables_["CCSD(T) CORRELATION ENERGY"] = Eccsd_t - Escf;
    variables_["CCSD(T) TOTAL ENERGY"] = Eccsd_t;
    variables_["(T) CORRECTION ENERGY"] = E_t;

    /* updates the wavefunction for checkpointing */
    energy_ = Process::environment.globals["CCSD(T) TOTAL ENERGY"];
    name_ = "CD-CCSD(T)";

    // CCSDL
    if (dertype == "FIRST" || cc_lambda_ == "TRUE") {
        // memalloc
        if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
            GtijA = std::make_shared<Tensor2d>("Gtilde Intermediate <I|J>", naoccA, naoccA);
            GtabA = std::make_shared<Tensor2d>("Gtilde Intermediate <A|B>", navirA, navirA);
            L1c = std::make_shared<Tensor1d>("DF_BASIS_CC L1_Q", nQ);
            gQt = std::make_shared<Tensor1d>("CCSD PDM G_Qt", nQ);
        }

        timer_on("CCSDL");
        if (t2_incore)
            ccsdl_iterations();
        else
            throw PSIEXCEPTION("There is NOT enough memory for Lambda equations!");
        timer_off("CCSDL");
        tstop();
        tstart();
    }

    // Compute Analytic Gradients
    if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
        // memalloc
        G1c_ov = SharedTensor2d(new Tensor2d("Correlation OPDM <O|V>", noccA, nvirA));
        G1c_vo = SharedTensor2d(new Tensor2d("Correlation OPDM <V|O>", nvirA, noccA));

        outfile->Printf("\tComputing unrelaxed response density matrices...\n");
        ccsd_opdm();
        ccsd_tpdm();
        ccl_energy();
        // prepare4grad();
        // if (oeprop_ == "TRUE") oeprop();
        // if (dertype == "FIRST") dfgrad();
        // if (ekt_ip_ == "TRUE") ekt_ip();
    }  // if (dertype == "FIRST" || ekt_ip_ == "TRUE")

}  // end ccsd_t_manager_cd

//======================================================================
//             Lambda-CCSD(T) Manager
//======================================================================
void DFOCC::ccsdl_t_manager_cd() {
    do_cd = "TRUE";
    time4grad = 0;     // means i will not compute the gradient
    mo_optimized = 0;  // means MOs are not optimized

    timer_on("CD Integrals");
    cd_ints();
    trans_cd();
    timer_off("CD Integrals");

    // Memory allocation
    T1c = std::make_shared<Tensor1d>("DF_BASIS_CC T1_Q", nQ);
    Jc = std::make_shared<Tensor1d>("DF_BASIS_SCF J_Q", nQ_ref);

    if (reference_ == "RESTRICTED") {
        t1A = std::make_shared<Tensor2d>("T1 <I|A>", naoccA, navirA);
        t1newA = std::make_shared<Tensor2d>("New T1 <I|A>", naoccA, navirA);
        FiaA = std::make_shared<Tensor2d>("Fint <I|A>", naoccA, navirA);
        FtijA = std::make_shared<Tensor2d>("Ftilde <I|J>", naoccA, naoccA);
        FtabA = std::make_shared<Tensor2d>("Ftilde <A|B>", navirA, navirA);

        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory / (1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements

        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = (long long int)nQ * (long long int)nso2_;
        cost_ampAA += (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_ampAA += (long long int)nQ * (long long int)nso_ * (long long int)navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;

        if ((cost_4amp + cost_df) <= memory_mb) {
            outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
            outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp + cost_df);
            nincore_amp = 4;
            t2_incore = true;
            df_ints_incore = true;
        } else if ((cost_3amp + cost_df) <= memory_mb) {
            outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
            // outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
            outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
            nincore_amp = 3;
            t2_incore = false;
            df_ints_incore = false;
        } else if (cost_3amp < memory_mb && cost_df < memory_mb) {
            outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
            outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
            nincore_amp = 3;
            t2_incore = false;
            df_ints_incore = false;
        } else {
            outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
            outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp + cost_df - memory_mb);
            throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        double cost_amp1 = 0.0;
        cost_amp1 = 2.5 * (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_amp1 += (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_amp1 /= 1024.0 * 1024.0;
        cost_amp1 *= sizeof(double);
        double cost_amp2 = 0.0;
        cost_amp2 = 1.5 * (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_amp2 += 4.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_amp2 /= 1024.0 * 1024.0;
        cost_amp2 *= sizeof(double);
        double cost_amp3 = 0.0;
        cost_amp3 = 2.0 * (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_amp3 += 3.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_amp3 += 2.0 * (long long int)navirA * (long long int)navirA * (long long int)navirA;
        cost_amp3 /= 1024.0 * 1024.0;
        cost_amp3 *= sizeof(double);
        cost_amp = MAX0(cost_amp1, cost_amp2);
        cost_amp = MAX0(cost_amp, cost_amp3);
        outfile->Printf("\tMemory requirement for Wabef term (T2): %9.2lf MB \n", cost_amp);
        if (cc_lambda_ == "TRUE") {
            cost_amp1 = (long long int)navirA * (long long int)navirA * (long long int)navirA;
            cost_amp1 /= 1024.0 * 1024.0;
            cost_amp += cost_amp1;
            outfile->Printf("\tMemory requirement for Wefab term (L2): %9.2lf MB \n", cost_amp);
        }

        // cost_ppl_hm
        cost_ppl_hm = ntri_abAA / 1024.0;
        cost_ppl_hm *= cost_ppl_hm;
        cost_ppl_hm *= sizeof(double);
        cost_ppl_hm += cost_amp;
        outfile->Printf("\tMemory for high mem Wabef algorithm   : %9.2lf MB \n", cost_ppl_hm);
        if (cost_ppl_hm > memory_mb && Wabef_type_ == "AUTO") {
            do_ppl_hm = false;
            outfile->Printf("\tI will use the LOW_MEM Wabef algorithm! \n");
        } else if (cost_ppl_hm <= memory_mb && Wabef_type_ == "AUTO") {
            do_ppl_hm = true;
            outfile->Printf("\tI will use the HIGH_MEM Wabef algorithm! \n");
        }

        // Memory for triples: 2*O^2V^2 + 5*V^3 + O^3V + V^2N + V^3/2
        cost_amp1 = 0.0;
        cost_amp1 = 3.0 * (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_amp1 += 5.0 * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_amp1 += (long long int)naoccA * (long long int)naoccA * (long long int)naoccA * (long long int)navirA;
        cost_amp1 += (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_amp1 += (long long int)navirA * (long long int)ntri_abAA;
        cost_amp1 /= 1024.0 * 1024.0;
        cost_amp1 *= sizeof(double);
        outfile->Printf("\tMemory requirement for (AT) correction : %9.2lf MB \n", cost_amp1);

        /*
        // Memory: OV^3 + 2*O^2V^2 + 2*V^3 + O^3V + V^2N
        cost_triples_iabc = 0.0;
        cost_triples_iabc = 2.0 * naoccA * naoccA * navirA * navirA;
        cost_triples_iabc += 5.0 * naoccA * navirA * navirA;
        cost_triples_iabc += naoccA * naoccA * naoccA * navirA;
        cost_triples_iabc += nQ * navirA * navirA;
        cost_triples_iabc /= 1024.0 * 1024.0;
        cost_amp2 = 0.0;
        cost_amp2 = (navirA * navirA) / 1024.0;
        cost_amp2 *= (naoccA * navirA) / 1024.0;
        cost_triples_iabc += cost_amp2;
        cost_triples_iabc *= sizeof(double);

        if (triples_iabc_type_ == "DISK") {
            do_triples_hm = false;
            //outfile->Printf("\n\tI will use a DISK algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_amp1);
        }
        else if (triples_iabc_type_ == "DIRECT") {
            do_triples_hm = false;
            outfile->Printf("\n\tI will use a DIRECT algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_amp1);
        }
        else if (cost_triples_iabc > memory_mb && triples_iabc_type_ == "AUTO") {
            do_triples_hm = false;
            outfile->Printf("\n\tI will use a DIRECT algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_amp1);
        }
        else if (cost_triples_iabc <= memory_mb && triples_iabc_type_ == "AUTO") {
            do_triples_hm = true;
            outfile->Printf("\n\tI will use an INCORE algorithm for (ia|bc) in (T)! \n");
            outfile->Printf("\tMemory requirement for (T) correction : %9.2lf MB \n", cost_triples_iabc);
        }
        */

        // Mem alloc for DF ints
        if (df_ints_incore) {
            bQijA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
            bQiaA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
            bQabA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
            bQijA->read(psio_, PSIF_DFOCC_INTS);
            bQiaA->read(psio_, PSIF_DFOCC_INTS);
            bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
        }

        //  Malloc
        if (t2_incore) {
            t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        }

    }  // end if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
        t1A = std::make_shared<Tensor2d>("T1 <I|A>", naoccA, navirA);
        t1B = std::make_shared<Tensor2d>("T1 <i|a>", naoccB, navirB);
        FiaA = std::make_shared<Tensor2d>("Fint <I|A>", naoccA, navirA);
        FiaB = std::make_shared<Tensor2d>("Fint <i|a>", naoccB, navirB);
        FtijA = std::make_shared<Tensor2d>("Ftilde <I|J>", naoccA, naoccA);
        FtabA = std::make_shared<Tensor2d>("Ftilde <A|B>", navirA, navirA);
        FtijB = std::make_shared<Tensor2d>("Ftilde <i|j>", naoccB, naoccB);
        FtabB = std::make_shared<Tensor2d>("Ftilde <a|b>", navirB, navirB);

        // memory requirements
        cost_ampAA = 0.0;
        cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_ampBB = (long long int)naocc2BB * (long long int)nvir2BB;
        cost_ampBB /= 1024.0 * 1024.0;
        cost_ampBB *= sizeof(double);
        cost_ampAB = (long long int)nocc2AB * (long long int)nvir2AB;
        cost_ampAB /= 1024.0 * 1024.0;
        cost_ampAB *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampBB);
        cost_amp = MAX0(cost_amp, cost_ampAB);
        cost_amp = 3.0 * cost_amp;
        memory = Process::environment.get_memory();
        memory_mb = (double)memory / (1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);
        outfile->Printf("\tMinimum required memory for amplitudes: %9.2lf MB \n", cost_amp);
    }  // else if (reference_ == "UNRESTRICTED")

    // memalloc for density intermediates
    if (qchf_ == "TRUE" || dertype == "FIRST") {
        g1Qc = std::make_shared<Tensor1d>("DF_BASIS_SCF G1_Q", nQ_ref);
        g1Qt = std::make_shared<Tensor1d>("DF_BASIS_SCF G1t_Q", nQ_ref);
        g1Qp = std::make_shared<Tensor1d>("DF_BASIS_SCF G1p_Q", nQ_ref);
        g1Q = std::make_shared<Tensor1d>("DF_BASIS_CC G1_Q", nQ);
        g1Qt2 = std::make_shared<Tensor1d>("DF_BASIS_CC G1t_Q", nQ);
    }

    // QCHF
    if (qchf_ == "TRUE") qchf();

    // Fock
    if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE") fock();

    // Compute MP2 energy
    if (reference == "ROHF") t1_1st_sc();
    if (t2_incore)
        ccsd_mp2();
    else
        ccsd_mp2_low();

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing DF-MP2 energy ... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    variables_["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;
    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        variables_["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    variables_["MP2 SINGLES ENERGY"] = Emp2_t1;

    // Perform CCSD iterations
    timer_on("CCSD");
    if (t2_incore)
        ccsd_iterations();
    else
        ccsd_iterations_low();
    timer_off("CCSD");

    outfile->Printf("\n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\t================ CCSD FINAL RESULTS =================================== \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    outfile->Printf("\tCD-CCSD Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-CCSD Total Energy (a.u.)        : %20.14f\n", Eccsd);
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\n");
    variables_["CCSD TOTAL ENERGY"] = Eccsd;
    variables_["CCSD CORRELATION ENERGY"] = Eccsd - Escf;
    if ((reference == "RHF") or (reference == "UHF")) {
        variables_["CCSD DOUBLES ENERGY"] = Eccsd - Escf;
        variables_["CCSD SINGLES ENERGY"] = 0.0;
    }

    // CCSDL
    timer_on("CCSDL");
    if (t2_incore) {
        tstop();
        tstart();
        lambda_title();
        outfile->Printf("\tSolving Lambda amplitude equations...\n");
        ccsdl_iterations();
    } else
        throw PSIEXCEPTION("There is NOT enough memory for Lambda equations!");
    timer_off("CCSDL");

    // CCSD(AT)
    tstop();
    tstart();
    pat_title();
    outfile->Printf("\tComputing asymmetric triples (AT) correction...\n");
    timer_on("(AT)");
    if (reference_ == "RESTRICTED") {
        ccsdl_canonic_triples_disk();
    }
    timer_off("(AT)");
    outfile->Printf("\t(AT) Correction (a.u.)             : %20.14f\n", E_at);
    outfile->Printf("\tCD-CCSD(AT) Total Energy (a.u.)    : %20.14f\n", Eccsd_at);

    variables_["CURRENT ENERGY"] = Eccsd_at;
    variables_["CURRENT REFERENCE ENERGY"] = Escf;
    variables_["CURRENT CORRELATION ENERGY"] = Eccsd_at - Escf;
    variables_["A-CCSD(T) TOTAL ENERGY"] = Eccsd_at;
    variables_["A-CCSD(T) CORRELATION ENERGY"] = Eccsd_at - Escf;
    variables_["A-(T) CORRECTION ENERGY"] = E_at;

    /* updates the wavefunction for checkpointing */
    energy_ = Eccsd_at;
    /*
    // Compute Analytic Gradients
    if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
        tstop();
        tstart();
        pdm_title();

        // memalloc
        G1c_ov = std::make_shared<Tensor2d>("Correlation OPDM <O|V>", noccA, nvirA);
        G1c_vo = std::make_shared<Tensor2d>("Correlation OPDM <V|O>", nvirA, noccA);

        outfile->Printf("\tComputing unrelaxed response density matrices...\n");
        ccsd_opdm();
        ccsd_tpdm();
        //ccl_energy();
        prepare4grad();
        if (oeprop_ == "TRUE") oeprop();
        if (dertype == "FIRST") dfgrad();
        //if (ekt_ip_ == "TRUE") ekt_ip();
    }// if (dertype == "FIRST" || ekt_ip_ == "TRUE")
    */

}  // end ccsd_t_manager_cd

//======================================================================
//             CCD Manager
//======================================================================
void DFOCC::ccd_manager_cd() {
    do_cd = "TRUE";
    time4grad = 0;     // means i will not compute the gradient
    mo_optimized = 0;  // means MOs are not optimized

    timer_on("CD Integrals");
    cd_ints();
    trans_cd();
    timer_off("CD Integrals");

    // Memory allocation
    // T1c = SharedTensor1d(new Tensor1d("DF_BASIS_CC T1_Q", nQ));
    Jc = std::make_shared<Tensor1d>("DF_BASIS_SCF J_Q", nQ_ref);

    if (reference_ == "RESTRICTED") {
        // avaliable mem
        memory = Process::environment.get_memory();
        memory_mb = (double)memory / (1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

        // memory requirements
        // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
        cost_df = 0.0;
        cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
        cost_df *= nQ;
        cost_df /= 1024.0 * 1024.0;
        cost_df *= sizeof(double);
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);

        // Cost of Integral transform for B(Q,ab)
        cost_ampAA = 0.0;
        cost_ampAA = (long long int)nQ * (long long int)nso2_;
        cost_ampAA += (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_ampAA += (long long int)nQ * (long long int)nso_ * (long long int)navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

        // Mem for amplitudes
        cost_ampAA = 0.0;
        cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_3amp = 3.0 * cost_ampAA;
        cost_4amp = 4.0 * cost_ampAA;
        cost_5amp = 5.0 * cost_ampAA;
        /*
        if ((cost_5amp+cost_df) <= memory_mb) {
             outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_5amp);
             outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_5amp+cost_df);
             nincore_amp = 5;
             t2_incore = true;
             df_ints_incore = true;
        }
        */
        if ((cost_4amp + cost_df) <= memory_mb) {
            outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
            outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp + cost_df);
            nincore_amp = 4;
            t2_incore = true;
            df_ints_incore = true;
        } else if ((cost_3amp + cost_df) <= memory_mb) {
            outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
            // outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
            outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
            nincore_amp = 3;
            t2_incore = false;
            df_ints_incore = false;
        } else if (cost_3amp < memory_mb && cost_df < memory_mb) {
            outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
            outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
            nincore_amp = 3;
            t2_incore = false;
            df_ints_incore = false;
        } else {
            outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
            outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp + cost_df - memory_mb);
            throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
        }

        // W_abef term
        cost_ampAA = 0.0;
        cost_ampAA = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_ampAA += 2.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_ampAA += (long long int)navirA * (long long int)navirA * (long long int)navirA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        double cost_ampAA2 = 0.0;
        cost_ampAA2 = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
        cost_ampAA2 += (long long int)nQ * (long long int)navirA * (long long int)navirA;
        cost_ampAA2 += 3.0 * (long long int)navirA * (long long int)navirA * (long long int)navirA;
        cost_ampAA2 /= 1024.0 * 1024.0;
        cost_ampAA2 *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampAA2);
        outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

        // cost_ppl_hm
        cost_ppl_hm = ntri_abAA / 1024.0;
        cost_ppl_hm *= cost_ppl_hm;
        cost_ppl_hm *= sizeof(double);
        cost_ppl_hm += cost_amp;
        outfile->Printf("\tMemory for high mem Wabef algorithm   : %9.2lf MB \n", cost_ppl_hm);
        if (cost_ppl_hm > memory_mb && Wabef_type_ == "AUTO") {
            do_ppl_hm = false;
            outfile->Printf("\tI will use the LOW_MEM Wabef algorithm! \n");
        } else if (cost_ppl_hm <= memory_mb && Wabef_type_ == "AUTO") {
            do_ppl_hm = true;
            outfile->Printf("\tI will use the HIGH_MEM Wabef algorithm! \n");
        }

        // Mem alloc for DF ints
        if (df_ints_incore) {
            bQijA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA));
            bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA));
            bQabA = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA));
            bQijA->read(psio_, PSIF_DFOCC_INTS);
            bQiaA->read(psio_, PSIF_DFOCC_INTS);
            bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
        }

        //  Malloc
        if (t2_incore) {
            t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
        }

    }  // end if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
        // memory requirements
        cost_ampAA = 0.0;
        cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
        cost_ampAA /= 1024.0 * 1024.0;
        cost_ampAA *= sizeof(double);
        cost_ampBB = (long long int)naocc2BB * (long long int)nvir2BB;
        cost_ampBB /= 1024.0 * 1024.0;
        cost_ampBB *= sizeof(double);
        cost_ampAB = (long long int)nocc2AB * (long long int)nvir2AB;
        cost_ampAB /= 1024.0 * 1024.0;
        cost_ampAB *= sizeof(double);
        cost_amp = MAX0(cost_ampAA, cost_ampBB);
        cost_amp = MAX0(cost_amp, cost_ampAB);
        cost_amp = 3.0 * cost_amp;
        memory = Process::environment.get_memory();
        memory_mb = (double)memory / (1024.0 * 1024.0);
        outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);
        outfile->Printf("\tMinimum required memory for amplitudes: %9.2lf MB \n", cost_amp);
    }  // else if (reference_ == "UNRESTRICTED")

    // memalloc for density intermediates
    if (qchf_ == "TRUE" || dertype == "FIRST") {
        g1Qc = std::make_shared<Tensor1d>("DF_BASIS_SCF G1_Q", nQ_ref);
        g1Qt = std::make_shared<Tensor1d>("DF_BASIS_SCF G1t_Q", nQ_ref);
        g1Qp = std::make_shared<Tensor1d>("DF_BASIS_SCF G1p_Q", nQ_ref);
        g1Q = std::make_shared<Tensor1d>("DF_BASIS_CC G1_Q", nQ);
        g1Qt2 = std::make_shared<Tensor1d>("DF_BASIS_CC G1t_Q", nQ);
    }

    // QCHF
    if (qchf_ == "TRUE") qchf();

    // Fock
    if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE") fock();

    // Compute MP2 energy
    if (reference == "ROHF") t1_1st_sc();
    if (t2_incore)
        ccd_mp2();
    else
        ccd_mp2_low();

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing CD-MP2 energy ... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;
    Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        Process::environment.globals["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        Process::environment.globals["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    Process::environment.globals["MP2 SINGLES ENERGY"] = Emp2_t1;

    // Perform CCD iterations
    timer_on("CCD");
    if (t2_incore)
        ccd_iterations();
    else
        ccd_iterations_low();
    timer_off("CCD");

    outfile->Printf("\n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\t================ CCD FINAL RESULTS ==================================== \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    outfile->Printf("\tCD-CCD Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-CCD Total Energy (a.u.)         : %20.14f\n", Eccd);
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\n");

    Process::environment.globals["CURRENT ENERGY"] = Eccd;
    Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
    Process::environment.globals["CURRENT CORRELATION ENERGY"] = Eccd - Escf;
    Process::environment.globals["CCD TOTAL ENERGY"] = Eccd;
    Process::environment.globals["CCD CORRELATION ENERGY"] = Eccd - Escf;

    // CCDL
    if (dertype == "FIRST" || cc_lambda_ == "TRUE") {
        // memalloc
        if (dertype == "FIRST") {
            gQt = SharedTensor1d(new Tensor1d("CCSD PDM G_Qt", nQ));
        }

        timer_on("CCDL");
        if (t2_incore) {
            tstop();
            tstart();
            lambda_title();
            outfile->Printf("\tSolving Lambda amplitude equations...\n");
            ccdl_iterations();
        } else
            throw PSIEXCEPTION("There is NOT enough memory for Lambda equations!");
        timer_off("CCDL");
    }

    // Compute Analytic Gradients
    if (dertype == "FIRST" || ekt_ip_ == "TRUE") {
        tstop();
        tstart();
        pdm_title();

        // memalloc
        G1c_ov = std::make_shared<Tensor2d>("Correlation OPDM <O|V>", noccA, nvirA);
        G1c_vo = std::make_shared<Tensor2d>("Correlation OPDM <V|O>", nvirA, noccA);

        outfile->Printf("\tComputing unrelaxed response density matrices...\n");
        ccd_opdm();
        ccd_tpdm();
        // ccl_energy();
        // prepare4grad();
        // if (oeprop_ == "TRUE") oeprop();
        // if (dertype == "FIRST") dfgrad();
        // if (ekt_ip_ == "TRUE") ekt_ip();
    }  // if (dertype == "FIRST" || ekt_ip_ == "TRUE")

}  // end ccd_manager_cd

//======================================================================
//             OMP3 Manager
//======================================================================
void DFOCC::omp3_manager_cd() {
    do_cd = "TRUE";
    time4grad = 0;         // means I will not compute the gradient
    mo_optimized = 0;      // means MOs are not optimized
    orbs_already_opt = 0;  // means orbitals are not optimized yet.
    orbs_already_sc = 0;   // means orbitals are not semicanonical yet.

    timer_on("CD Integrals");
    cd_ints();
    trans_cd();
    timer_off("CD Integrals");

    // memalloc for density intermediates
    Jc = std::make_shared<Tensor1d>("DF_BASIS_SCF J_Q", nQ_ref);
    g1Qc = std::make_shared<Tensor1d>("DF_BASIS_SCF G1_Q", nQ_ref);
    g1Qt = std::make_shared<Tensor1d>("DF_BASIS_SCF G1t_Q", nQ_ref);
    g1Qp = std::make_shared<Tensor1d>("DF_BASIS_SCF G1p_Q", nQ_ref);
    g1Q = std::make_shared<Tensor1d>("DF_BASIS_CC G1_Q", nQ);
    g1Qt2 = std::make_shared<Tensor1d>("DF_BASIS_CC G1t_Q", nQ);

    // avaliable mem
    memory = Process::environment.get_memory();
    memory_mb = (double)memory / (1024.0 * 1024.0);
    outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

    // memory requirements
    // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
    cost_df = 0.0;
    cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
    cost_df *= nQ;
    cost_df /= 1024.0 * 1024.0;
    cost_df *= sizeof(double);
    if (reference_ == "RESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
    else if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0 * cost_df);

    // Cost of Integral transform for B(Q,ab)
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)nQ * (long long int)nso2_;
    cost_ampAA += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)nQ * (long long int)nso_ * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

    // Mem for amplitudes
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    cost_3amp = 3.0 * cost_ampAA;
    cost_4amp = 4.0 * cost_ampAA;
    cost_5amp = 5.0 * cost_ampAA;

    if ((cost_4amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
        outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp + cost_df);
        nincore_amp = 4;
        t2_incore = true;
        df_ints_incore = true;
    } else if ((cost_3amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        // outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else if (cost_3amp < memory_mb && cost_df < memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else {
        outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
        outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp + cost_df - memory_mb);
        throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
    }

    // W_abef term
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA += 2.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    double cost_ampAA2 = 0.0;
    cost_ampAA2 = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += 3.0 * (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 /= 1024.0 * 1024.0;
    cost_ampAA2 *= sizeof(double);
    cost_amp = MAX0(cost_ampAA, cost_ampAA2);
    outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

    // Fock
    fock();

    // QCHF
    if (qchf_ == "TRUE") qchf();

    // ROHF REF
    if (reference == "ROHF") {
        t1A = std::make_shared<Tensor2d>("T1_1 <I|A>", naoccA, navirA);
        t1B = std::make_shared<Tensor2d>("T1_1 <i|a>", naoccB, navirB);
        t1_1st_sc();
    }
    mp3_t2_1st_sc();
    Emp2L = Emp2;
    EcorrL = Emp2L - Escf;
    Emp2L_old = Emp2;

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (Canonical CD-MP2)... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    variables_["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;
    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        variables_["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    variables_["MP2 SINGLES ENERGY"] = Emp2_t1;

    // Perform MP3 iterations
    timer_on("MP3");
    t2_2nd_sc();
    timer_off("MP3");

    if (reference != "ROHF") {
        // LAB: dfocc can't compute a ROHF-MP3 explicitly, so it's unlikely to be a correct byproduct
    outfile->Printf("\n");
    outfile->Printf("\tComputing CD-MP3 energy using SCF MOs (Canonical CD-MP3)... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
    outfile->Printf("\t3rd Order Energy (a.u.)            : %20.14f\n", Emp3 - Emp2);
    outfile->Printf("\tCD-MP2.5 Correlation Energy (a.u.) : %20.14f\n", (Emp2 - Escf) + 0.5 * (Emp3 - Emp2));
    outfile->Printf("\tCD-MP2.5 Total Energy (a.u.)       : %20.14f\n", 0.5 * (Emp3 + Emp2));
    outfile->Printf("\tCD-MP3 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP3 Total Energy (a.u.)         : %20.14f\n", Emp3);
    outfile->Printf("\t======================================================================= \n");

    variables_["MP3 TOTAL ENERGY"] = Emp3;
    variables_["MP3 CORRELATION ENERGY"] = Emp3 - Escf;
        if (reference_ == "UNRESTRICTED") {
            variables_["MP3 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp3AB;
            variables_["MP3 SAME-SPIN CORRELATION ENERGY"] = Emp3AA + Emp3BB;
        }
        // RHF & UHF only
        variables_["MP3 SINGLES ENERGY"] = 0.0;
        variables_["MP3 DOUBLES ENERGY"] = Ecorr;
    }
    Emp3L = Emp3;
    EcorrL = Emp3L - Escf;
    Emp3L_old = Emp3;

    // Malloc for PDMs
    gQt = std::make_shared<Tensor1d>("CCD PDM G_Qt", nQ);
    if (reference_ == "RESTRICTED") {
        G1c_ov = std::make_shared<Tensor2d>("Correlation OPDM <O|V>", noccA, nvirA);
        G1c_vo = std::make_shared<Tensor2d>("Correlation OPDM <V|O>", nvirA, noccA);
    } else if (reference_ == "UNRESTRICTED") {
        G1c_ovA = std::make_shared<Tensor2d>("Correlation OPDM <O|V>", noccA, nvirA);
        G1c_ovB = std::make_shared<Tensor2d>("Correlation OPDM <o|v>", noccB, nvirB);
        G1c_voA = std::make_shared<Tensor2d>("Correlation OPDM <V|O>", nvirA, noccA);
        G1c_voB = std::make_shared<Tensor2d>("Correlation OPDM <v|o>", nvirB, noccB);
    }

    mp3_pdm_3index_intr();
    omp3_opdm();
    omp3_tpdm();
    // ccl_energy();
    sep_tpdm_cc();
    gfock_cc_vo();
    gfock_cc_ov();
    gfock_cc_oo();
    gfock_cc_vv();
    idp();
    mograd();
    occ_iterations();

    // main if
    if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod) {
        orbs_already_opt = 1;
        if (conver == 1)
            outfile->Printf("\n\tOrbitals are optimized now.\n");
        else if (conver == 0) {
            outfile->Printf("\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
        }
        outfile->Printf("\tTransforming MOs to the semicanonical basis... \n");
        semi_canonic();

        outfile->Printf("\tSwitching to the standard DF-MP3 computation... \n");
        trans_cd();
        fock();
        ref_energy();
        mp3_t2_1st_sc();
        t2_2nd_sc();
        conver = 1;
        if (dertype == "FIRST") {
            mp3_pdm_3index_intr();
            omp3_opdm();
            omp3_tpdm();
            sep_tpdm_cc();
            gfock_cc_vo();
            gfock_cc_ov();
            gfock_cc_oo();
            gfock_cc_vv();
        }
    }  // end main if

    else if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod && regularization == "TRUE") {
        outfile->Printf("\tOrbital gradient converged, but energy did not... \n");
        outfile->Printf("\tA tighter rms_mograd_convergence tolerance is recommended... \n");
        throw PSIEXCEPTION("A tighter rms_mograd_convergence tolerance is recommended.");
    }

    if (conver == 1) {
        if (orbs_already_opt == 1) Emp3L = Emp3;

        outfile->Printf("\n");
        outfile->Printf("\tComputing CD-MP3 energy using optimized MOs... \n");
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
        outfile->Printf("\t3rd Order Energy (a.u.)            : %20.14f\n", Emp3 - Emp2);
        outfile->Printf("\tCD-MP2.5 Correlation Energy (a.u.) : %20.14f\n", (Emp2 - Escf) + 0.5 * (Emp3 - Emp2));
        outfile->Printf("\tCD-MP2.5 Total Energy (a.u.)       : %20.14f\n", 0.5 * (Emp3 + Emp2));
        outfile->Printf("\tCD-MP3 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
        outfile->Printf("\tCD-MP3 Total Energy (a.u.)         : %20.14f\n", Emp3);
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\n");

        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\t================ CD-OMP3 FINAL RESULTS ================================ \n");
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tCD-OMP3 Correlation Energy (a.u.)  : %20.14f\n", Emp3L - Escf);
        outfile->Printf("\tEcdomp3 - Eref (a.u.)              : %20.14f\n", Emp3L - Eref);
        outfile->Printf("\tCD-OMP3 Total Energy (a.u.)        : %20.14f\n", Emp3L);
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\n");

        // Set the global variables with the energies
        // Emp3L=Emp3;
        variables_["CURRENT ENERGY"] = Emp3L;
        variables_["CURRENT REFERENCE ENERGY"] = Escf;
        variables_["CURRENT CORRELATION ENERGY"] = Emp3L - Escf;
        variables_["OMP3 TOTAL ENERGY"] = Emp3L;
        variables_["OMP3 CORRELATION ENERGY"] = Emp3L - Escf;
        variables_["OMP3 REFERENCE CORRECTION ENERGY"] = Eref - Escf;

        if (reference_ == "UNRESTRICTED") {
            variables_["OMP3 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp3AB;
            variables_["OMP3 SAME-SPIN CORRELATION ENERGY"] = Emp3AA + Emp3BB;
        }

        /* updates the wavefunction for checkpointing */
        energy_ = variables_["CURRENT ENERGY"];
        name_ = "CD-OMP3";

        response_helper();

        // Save MOs to wfn
        save_mo_to_wfn();

    }  // end if (conver == 1)

}  // end omp3_manager_cd

//======================================================================
//             MP3 Manager
//======================================================================
void DFOCC::mp3_manager_cd() {
    do_cd = "TRUE";
    time4grad = 0;     // means i will not compute the gradient
    mo_optimized = 0;  // means MOs are not optimized

    timer_on("CD Integrals");
    cd_ints();
    trans_cd();
    timer_off("CD Integrals");

    // Memory allocation
    Jc = std::make_shared<Tensor1d>("DF_BASIS_SCF J_Q", nQ_ref);

    // avaliable mem
    memory = Process::environment.get_memory();
    memory_mb = (double)memory / (1024.0 * 1024.0);
    outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

    // memory requirements
    // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
    cost_df = 0.0;
    cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
    cost_df *= nQ;
    cost_df /= 1024.0 * 1024.0;
    cost_df *= sizeof(double);
    if (reference_ == "RESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
    else if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0 * cost_df);

    // Cost of Integral transform for B(Q,ab)
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)nQ * (long long int)nso2_;
    cost_ampAA += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)nQ * (long long int)nso_ * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

    // Mem for amplitudes
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    cost_3amp = 3.0 * cost_ampAA;
    cost_4amp = 4.0 * cost_ampAA;
    cost_5amp = 5.0 * cost_ampAA;

    if ((cost_4amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
        outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp + cost_df);
        nincore_amp = 4;
        t2_incore = true;
        df_ints_incore = true;
    } else if ((cost_3amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        // outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else if (cost_3amp < memory_mb && cost_df < memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else {
        outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
        outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp + cost_df - memory_mb);
        throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
    }

    // W_abef term
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA += 2.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    double cost_ampAA2 = 0.0;
    cost_ampAA2 = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += 3.0 * (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 /= 1024.0 * 1024.0;
    cost_ampAA2 *= sizeof(double);
    cost_amp = MAX0(cost_ampAA, cost_ampAA2);
    outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

    // Mem alloc for DF ints
    /*
    if (df_ints_incore) {
        bQijA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
        bQiaA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IA)", nQ, naoccA, navirA);
        bQabA = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|AB)", nQ, navirA, navirA);
        bQijA->read(psio_, PSIF_DFOCC_INTS);
        bQiaA->read(psio_, PSIF_DFOCC_INTS);
        bQabA->read(psio_, PSIF_DFOCC_INTS, true, true);
    }
    */

    /*
    if (t2_incore) {
        t2 = std::make_shared<Tensor2d>("T2 (IA|JB)", naoccA, navirA, naoccA, navirA);
    }
    */

    // memalloc for density intermediates
    if (qchf_ == "TRUE" || dertype == "FIRST") {
        g1Qc = std::make_shared<Tensor1d>("DF_BASIS_SCF G1_Q", nQ_ref);
        g1Qt = std::make_shared<Tensor1d>("DF_BASIS_SCF G1t_Q", nQ_ref);
        g1Qp = std::make_shared<Tensor1d>("DF_BASIS_SCF G1p_Q", nQ_ref);
        g1Q = std::make_shared<Tensor1d>("DF_BASIS_CC G1_Q", nQ);
        g1Qt2 = std::make_shared<Tensor1d>("DF_BASIS_CC G1t_Q", nQ);
    }

    // QCHF
    if (qchf_ == "TRUE") qchf();

    // Compute MP2 energy
    if (reference == "ROHF") {
        t1A = std::make_shared<Tensor2d>("T1_1 <I|A>", naoccA, navirA);
        t1B = std::make_shared<Tensor2d>("T1_1 <i|a>", naoccB, navirB);
        t1_1st_sc();
    }
    mp3_t2_1st_sc();

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing CD-MP2 energy ... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;
    Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        Process::environment.globals["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        Process::environment.globals["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    Process::environment.globals["MP2 SINGLES ENERGY"] = Emp2_t1;

    // Perform MP3 iterations
    timer_on("MP3");
    t2_2nd_sc();
    timer_off("MP3");

    outfile->Printf("\n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\t================ MP3 FINAL RESULTS ==================================== \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
    outfile->Printf("\t3rd Order Energy (a.u.)            : %20.14f\n", Emp3 - Emp2);
    outfile->Printf("\tCD-MP2.5 Correlation Energy (a.u.) : %20.14f\n", (Emp2 - Escf) + 0.5 * (Emp3 - Emp2));
    outfile->Printf("\tCD-MP2.5 Total Energy (a.u.)       : %20.14f\n", 0.5 * (Emp3 + Emp2));
    outfile->Printf("\tCD-MP3 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP3 Total Energy (a.u.)         : %20.14f\n", Emp3);
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\n");

    variables_["MP2 TOTAL ENERGY"] = Emp2;
    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    variables_["MP2 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    variables_["MP2 DOUBLES ENERGY"] = Emp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
    }

    variables_["MP2.5 TOTAL ENERGY"] = 0.5 * (Emp3 + Emp2);
    variables_["MP2.5 CORRELATION ENERGY"] = (Emp2 - Escf) + 0.5 * (Emp3 - Emp2);
    variables_["MP2.5 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    variables_["MP2.5 DOUBLES ENERGY"] = (Emp2 - Escf) + 0.5 * (Emp3 - Emp2);
    if (reference_ == "UNRESTRICTED") {
        variables_["MP2.5 SAME-SPIN CORRELATION ENERGY"] = 0.5 * (Emp2AA + Emp2BB + Emp3AA + Emp3BB);
        variables_["MP2.5 OPPOSITE-SPIN CORRELATION ENERGY"] = 0.5 * (Emp2AB + Emp3AB);
    }

    variables_["CURRENT ENERGY"] = Emp3;
    variables_["CURRENT REFERENCE ENERGY"] = Escf;
    variables_["CURRENT CORRELATION ENERGY"] = Emp3 - Escf;
    variables_["MP3 TOTAL ENERGY"] = Emp3;
    variables_["MP3 CORRELATION ENERGY"] = Emp3 - Escf;
    variables_["MP3 DOUBLES ENERGY"] = Emp3 - Escf;
    variables_["MP3 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    if (reference_ == "UNRESTRICTED") {
        variables_["MP3 SAME-SPIN CORRELATION ENERGY"] = Emp3AA + Emp3BB;
        variables_["MP3 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp3AB;
    }

    Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["MP2 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    Process::environment.globals["MP2 DOUBLES ENERGY"] = Emp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
    }

    Process::environment.globals["MP2.5 TOTAL ENERGY"] = 0.5 * (Emp3 + Emp2);
    Process::environment.globals["MP2.5 CORRELATION ENERGY"] = (Emp2 - Escf) + 0.5 * (Emp3 - Emp2);
    Process::environment.globals["MP2.5 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    Process::environment.globals["MP2.5 DOUBLES ENERGY"] = (Emp2 - Escf) + 0.5 * (Emp3 - Emp2);
    if (reference_ == "UNRESTRICTED") {
        Process::environment.globals["MP2.5 SAME-SPIN CORRELATION ENERGY"] = 0.5 * (Emp2AA + Emp2BB + Emp3AA + Emp3BB);
        Process::environment.globals["MP2.5 OPPOSITE-SPIN CORRELATION ENERGY"] = 0.5 * (Emp2AB + Emp3AB);
    }

    Process::environment.globals["CURRENT ENERGY"] = Emp3;
    Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
    Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3 - Escf;
    Process::environment.globals["MP3 TOTAL ENERGY"] = Emp3;
    Process::environment.globals["MP3 CORRELATION ENERGY"] = Emp3 - Escf;
    Process::environment.globals["MP3 DOUBLES ENERGY"] = Emp3 - Escf;
    Process::environment.globals["MP3 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    if (reference_ == "UNRESTRICTED") {
        Process::environment.globals["MP3 SAME-SPIN CORRELATION ENERGY"] = Emp3AA + Emp3BB;
        Process::environment.globals["MP3 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp3AB;
    }

    /* updates the wavefunction for checkpointing */
    energy_ = Emp3;
    name_ = "CD-MP3";
    Emp3L = Emp3;


}  // end mp3_manager_cd

//======================================================================
//             OMP2.5 Manager
//======================================================================
void DFOCC::omp2_5_manager_cd() {
    do_cd = "TRUE";
    time4grad = 0;         // means I will not compute the gradient
    mo_optimized = 0;      // means MOs are not optimized
    orbs_already_opt = 0;  // means orbitals are not optimized yet.
    orbs_already_sc = 0;   // means orbitals are not semicanonical yet.

    timer_on("CD Integrals");
    cd_ints();
    trans_cd();
    timer_off("CD Integrals");

    // memalloc for density intermediates
    Jc = std::make_shared<Tensor1d>("DF_BASIS_SCF J_Q", nQ_ref);
    g1Qc = std::make_shared<Tensor1d>("DF_BASIS_SCF G1_Q", nQ_ref);
    g1Qt = std::make_shared<Tensor1d>("DF_BASIS_SCF G1t_Q", nQ_ref);
    g1Qp = std::make_shared<Tensor1d>("DF_BASIS_SCF G1p_Q", nQ_ref);
    g1Q = std::make_shared<Tensor1d>("DF_BASIS_CC G1_Q", nQ);
    g1Qt2 = std::make_shared<Tensor1d>("DF_BASIS_CC G1t_Q", nQ);

    // avaliable mem
    memory = Process::environment.get_memory();
    memory_mb = (double)memory / (1024.0 * 1024.0);
    outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

    // memory requirements
    // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
    cost_df = 0.0;
    cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
    cost_df *= nQ;
    cost_df /= 1024.0 * 1024.0;
    cost_df *= sizeof(double);
    if (reference_ == "RESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
    else if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0 * cost_df);

    // Cost of Integral transform for B(Q,ab)
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)nQ * (long long int)nso2_;
    cost_ampAA += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)nQ * (long long int)nso_ * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

    // Mem for amplitudes
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    cost_3amp = 3.0 * cost_ampAA;
    cost_4amp = 4.0 * cost_ampAA;
    cost_5amp = 5.0 * cost_ampAA;

    if ((cost_4amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
        outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp + cost_df);
        nincore_amp = 4;
        t2_incore = true;
        df_ints_incore = true;
    } else if ((cost_3amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        // outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else if (cost_3amp < memory_mb && cost_df < memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else {
        outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
        outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp + cost_df - memory_mb);
        throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
    }

    // W_abef term
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA += 2.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    double cost_ampAA2 = 0.0;
    cost_ampAA2 = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += 3.0 * (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 /= 1024.0 * 1024.0;
    cost_ampAA2 *= sizeof(double);
    cost_amp = MAX0(cost_ampAA, cost_ampAA2);
    outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

    // Fock
    fock();

    // QCHF
    if (qchf_ == "TRUE") qchf();

    // ROHF REF
    if (reference == "ROHF") {
        t1A = std::make_shared<Tensor2d>("T1_1 <I|A>", naoccA, navirA);
        t1B = std::make_shared<Tensor2d>("T1_1 <i|a>", naoccB, navirB);
        t1_1st_sc();
    }
    mp3_t2_1st_sc();
    Emp2L = Emp2;
    EcorrL = Emp2L - Escf;
    Emp2L_old = Emp2;

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (Canonical CD-MP2)... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    variables_["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;
    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        variables_["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    variables_["MP2 SINGLES ENERGY"] = Emp2_t1;

    // Perform MP3 iterations
    timer_on("MP3");
    t2_2nd_sc();
    timer_off("MP3");

    if (reference != "ROHF") {
        // LAB: dfocc can't compute a ROHF-MP3 explicitly, so it's unlikely to be a correct byproduct
    outfile->Printf("\n");
    outfile->Printf("\tComputing CD-MP2.5 energy using SCF MOs (Canonical CD-MP2.5)... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
    outfile->Printf("\tCD-MP3 Correlation Energy (a.u.)   : %20.14f\n", (Emp2 - Escf) + 2.0 * (Emp3 - Emp2));
    outfile->Printf("\tCD-MP3 Total Energy (a.u.)         : %20.14f\n", Emp2 + 2.0 * (Emp3 - Emp2));
    outfile->Printf("\tCD-MP2.5 Correlation Energy (a.u.) : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2.5 Total Energy (a.u.)       : %20.14f\n", Emp3);
    outfile->Printf("\t======================================================================= \n");

    variables_["MP2.5 TOTAL ENERGY"] = Emp3;
    variables_["MP2.5 CORRELATION ENERGY"] = Emp3 - Escf;
    variables_["MP2.5 SINGLES ENERGY"] = 0.0;
    variables_["MP2.5 DOUBLES ENERGY"] = variables_["MP2.5 CORRELATION ENERGY"];  // RHF & UHF only
    if (reference_ == "UNRESTRICTED") {
        variables_["MP2.5 SAME-SPIN CORRELATION ENERGY"] = Emp3AA + Emp3BB;
        variables_["MP2.5 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp3AB;
    }

    variables_["MP3 TOTAL ENERGY"] = Emp2 + 2.0 * (Emp3 - Emp2);
    variables_["MP3 CORRELATION ENERGY"] = Emp2 + 2.0 * (Emp3 - Emp2) - Escf;
    variables_["MP3 DOUBLES ENERGY"] = variables_["MP3 CORRELATION ENERGY"];  // RHF & UHF only
    variables_["MP3 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    if (reference_ == "UNRESTRICTED") {
        variables_["MP3 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB + 2.0 * (Emp3AA + Emp3BB - Emp2AA - Emp2BB);
        variables_["MP3 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB + 2.0 * (Emp3AB - Emp2AB);
    }
    }
    Emp3L = Emp3;
    EcorrL = Emp3L - Escf;
    Emp3L_old = Emp3;

    // Malloc for PDMs
    gQt = std::make_shared<Tensor1d>("CCD PDM G_Qt", nQ);
    if (reference_ == "RESTRICTED") {
        G1c_ov = std::make_shared<Tensor2d>("Correlation OPDM <O|V>", noccA, nvirA);
        G1c_vo = std::make_shared<Tensor2d>("Correlation OPDM <V|O>", nvirA, noccA);
    } else if (reference_ == "UNRESTRICTED") {
        G1c_ovA = std::make_shared<Tensor2d>("Correlation OPDM <O|V>", noccA, nvirA);
        G1c_ovB = std::make_shared<Tensor2d>("Correlation OPDM <o|v>", noccB, nvirB);
        G1c_voA = std::make_shared<Tensor2d>("Correlation OPDM <V|O>", nvirA, noccA);
        G1c_voB = std::make_shared<Tensor2d>("Correlation OPDM <v|o>", nvirB, noccB);
    }

    mp3_pdm_3index_intr();
    omp3_opdm();
    omp3_tpdm();
    // ccl_energy();
    sep_tpdm_cc();
    gfock_cc_vo();
    gfock_cc_ov();
    gfock_cc_oo();
    gfock_cc_vv();
    idp();
    mograd();
    occ_iterations();

    // main if
    if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod) {
        orbs_already_opt = 1;
        if (conver == 1)
            outfile->Printf("\n\tOrbitals are optimized now.\n");
        else if (conver == 0) {
            outfile->Printf("\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
        }
        outfile->Printf("\tTransforming MOs to the semicanonical basis... \n");
        semi_canonic();

        outfile->Printf("\tSwitching to the standard DF-MP3 computation... \n");
        trans_cd();
        fock();
        ref_energy();
        mp3_t2_1st_sc();
        t2_2nd_sc();
        conver = 1;
        if (dertype == "FIRST") {
            mp3_pdm_3index_intr();
            omp3_opdm();
            omp3_tpdm();
            sep_tpdm_cc();
            gfock_cc_vo();
            gfock_cc_ov();
            gfock_cc_oo();
            gfock_cc_vv();
        }
    }  // end main if

    else if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod && regularization == "TRUE") {
        outfile->Printf("\tOrbital gradient converged, but energy did not... \n");
        outfile->Printf("\tA tighter rms_mograd_convergence tolerance is recommended... \n");
        throw PSIEXCEPTION("A tighter rms_mograd_convergence tolerance is recommended.");
    }

    if (conver == 1) {
        if (orbs_already_opt == 1) Emp3L = Emp3;

        outfile->Printf("\n");
        outfile->Printf("\tComputing CD-MP2.5 energy using optimized MOs... \n");
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
        outfile->Printf("\tCD-MP3 Correlation Energy (a.u.)   : %20.14f\n", (Emp2 - Escf) + 2.0 * (Emp3 - Emp2));
        outfile->Printf("\tCD-MP3 Total Energy (a.u.)         : %20.14f\n", Emp2 + 2.0 * (Emp3 - Emp2));
        outfile->Printf("\tCD-MP2.5 Correlation Energy (a.u.) : %20.14f\n", Ecorr);
        outfile->Printf("\tCD-MP2.5 Total Energy (a.u.)       : %20.14f\n", Emp3);
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\n");

        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\t================ CD-OMP2.5 FINAL RESULTS ============================== \n");
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tCD-OMP2.5 Correlation Energy (a.u.): %20.14f\n", Emp3L - Escf);
        outfile->Printf("\tEcdomp2.5 - Eref (a.u.)            : %20.14f\n", Emp3L - Eref);
        outfile->Printf("\tCD-OMP2.5 Total Energy (a.u.)      : %20.14f\n", Emp3L);
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\n");

        // Set the global variables with the energies
        // Emp3L=Emp3;
        variables_["CURRENT ENERGY"] = Emp3L;
        variables_["CURRENT REFERENCE ENERGY"] = Escf;
        variables_["CURRENT CORRELATION ENERGY"] = Emp3L - Escf;
        variables_["OMP2.5 TOTAL ENERGY"] = Emp3L;
        variables_["OMP2.5 CORRELATION ENERGY"] = Emp3L - Escf;
        variables_["OMP2.5 REFERENCE CORRECTION ENERGY"] = Eref - Escf;
        energy_ = Emp3L;

        if (reference_ == "UNRESTRICTED") {
            variables_["OMP2.5 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp3AB;
            variables_["OMP2.5 SAME-SPIN CORRELATION ENERGY"] = Emp3AA + Emp3BB;
        }

        // Save MOs to wfn
        save_mo_to_wfn();

        response_helper();

    }  // end if (conver == 1)

}  // end omp2_5_manager_cd

//======================================================================
//             MP2.5 Manager
//======================================================================
void DFOCC::mp2_5_manager_cd() {
    do_cd = "TRUE";
    time4grad = 0;     // means i will not compute the gradient
    mo_optimized = 0;  // means MOs are not optimized

    timer_on("CD Integrals");
    cd_ints();
    trans_cd();
    timer_off("CD Integrals");

    // Memory allocation
    Jc = std::make_shared<Tensor1d>("DF_BASIS_SCF J_Q", nQ_ref);

    // avaliable mem
    memory = Process::environment.get_memory();
    memory_mb = (double)memory / (1024.0 * 1024.0);
    outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

    // memory requirements
    // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
    cost_df = 0.0;
    cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
    cost_df *= nQ;
    cost_df /= 1024.0 * 1024.0;
    cost_df *= sizeof(double);
    if (reference_ == "RESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
    else if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0 * cost_df);

    // Cost of Integral transform for B(Q,ab)
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)nQ * (long long int)nso2_;
    cost_ampAA += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)nQ * (long long int)nso_ * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

    // Mem for amplitudes
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    cost_3amp = 3.0 * cost_ampAA;
    cost_4amp = 4.0 * cost_ampAA;
    cost_5amp = 5.0 * cost_ampAA;

    if ((cost_4amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
        outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp + cost_df);
        nincore_amp = 4;
        t2_incore = true;
        df_ints_incore = true;
    } else if ((cost_3amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        // outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else if (cost_3amp < memory_mb && cost_df < memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else {
        outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
        outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp + cost_df - memory_mb);
        throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
    }

    // W_abef term
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA += 2.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    double cost_ampAA2 = 0.0;
    cost_ampAA2 = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += 3.0 * (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 /= 1024.0 * 1024.0;
    cost_ampAA2 *= sizeof(double);
    cost_amp = MAX0(cost_ampAA, cost_ampAA2);
    outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

    // memalloc for density intermediates
    if (qchf_ == "TRUE" || dertype == "FIRST") {
        g1Qc = std::make_shared<Tensor1d>("DF_BASIS_SCF G1_Q", nQ_ref);
        g1Qt = std::make_shared<Tensor1d>("DF_BASIS_SCF G1t_Q", nQ_ref);
        g1Qp = std::make_shared<Tensor1d>("DF_BASIS_SCF G1p_Q", nQ_ref);
        g1Q = std::make_shared<Tensor1d>("DF_BASIS_CC G1_Q", nQ);
        g1Qt2 = std::make_shared<Tensor1d>("DF_BASIS_CC G1t_Q", nQ);
    }

    // QCHF
    if (qchf_ == "TRUE") qchf();

    // Compute MP2 energy
    if (reference == "ROHF") {
        t1A = std::make_shared<Tensor2d>("T1_1 <I|A>", naoccA, navirA);
        t1B = std::make_shared<Tensor2d>("T1_1 <i|a>", naoccB, navirB);
        t1_1st_sc();
    }
    mp3_t2_1st_sc();

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing CD-MP2 energy ... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    Process::environment.globals["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;
    Process::environment.globals["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        Process::environment.globals["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        Process::environment.globals["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    Process::environment.globals["MP2 SINGLES ENERGY"] = Emp2_t1;

    // Perform MP3 iterations
    timer_on("MP3");
    t2_2nd_sc();
    timer_off("MP3");

    outfile->Printf("\n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\t================ MP2.5 FINAL RESULTS ================================== \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp3AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp3AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp3BB);
    outfile->Printf("\tCD-MP3 Correlation Energy (a.u.)   : %20.14f\n", (Emp2 - Escf) + 2.0 * (Emp3 - Emp2));
    outfile->Printf("\tCD-MP3 Total Energy (a.u.)         : %20.14f\n", Emp2 + 2.0 * (Emp3 - Emp2));
    outfile->Printf("\tCD-MP2.5 Correlation Energy (a.u.) : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2.5 Total Energy (a.u.)       : %20.14f\n", Emp3);
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\n");

    variables_["MP2 TOTAL ENERGY"] = Emp2;
    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    variables_["MP2 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    variables_["MP2 DOUBLES ENERGY"] = Emp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
    }

    variables_["CURRENT ENERGY"] = Emp3;
    variables_["CURRENT REFERENCE ENERGY"] = Escf;
    variables_["CURRENT CORRELATION ENERGY"] = Emp3 - Escf;
    variables_["MP2.5 TOTAL ENERGY"] = Emp3;
    variables_["MP2.5 CORRELATION ENERGY"] = Emp3 - Escf;
    variables_["MP2.5 SINGLES ENERGY"] = 0.0;
    variables_["MP2.5 DOUBLES ENERGY"] = variables_["MP2.5 CORRELATION ENERGY"];  // RHF & UHF only
    if (reference_ == "UNRESTRICTED") {
        variables_["MP2.5 SAME-SPIN CORRELATION ENERGY"] = Emp3AA + Emp3BB;
        variables_["MP2.5 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp3AB;
    }

    variables_["MP3 TOTAL ENERGY"] = Emp2 + 2.0 * (Emp3 - Emp2);
    variables_["MP3 CORRELATION ENERGY"] = Emp2 + 2.0 * (Emp3 - Emp2) - Escf;
    variables_["MP3 DOUBLES ENERGY"] = variables_["MP3 CORRELATION ENERGY"];  // RHF & UHF only
    variables_["MP3 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    if (reference_ == "UNRESTRICTED") {
        variables_["MP3 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB + 2.0 * (Emp3AA + Emp3BB - Emp2AA - Emp2BB);
        variables_["MP3 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB + 2.0 * (Emp3AB - Emp2AB);
    }

    Process::environment.globals["MP2 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    Process::environment.globals["MP2 DOUBLES ENERGY"] = Emp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
    }

    Process::environment.globals["CURRENT ENERGY"] = Emp3;
    Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
    Process::environment.globals["CURRENT CORRELATION ENERGY"] = Emp3 - Escf;
    Process::environment.globals["MP2.5 TOTAL ENERGY"] = Emp3;
    Process::environment.globals["MP2.5 CORRELATION ENERGY"] = Emp3 - Escf;
    Process::environment.globals["MP2.5 SINGLES ENERGY"] = 0.0;
    Process::environment.globals["MP2.5 DOUBLES ENERGY"] = Process::environment.globals["MP2.5 CORRELATION ENERGY"];  // RHF & UHF only
    if (reference_ == "UNRESTRICTED") {
        Process::environment.globals["MP2.5 SAME-SPIN CORRELATION ENERGY"] = Emp3AA + Emp3BB;
        Process::environment.globals["MP2.5 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp3AB;
    }

    Process::environment.globals["MP3 TOTAL ENERGY"] = Emp2 + 2.0 * (Emp3 - Emp2);
    Process::environment.globals["MP3 CORRELATION ENERGY"] = Emp2 + 2.0 * (Emp3 - Emp2) - Escf;
    Process::environment.globals["MP3 DOUBLES ENERGY"] = Process::environment.globals["MP3 CORRELATION ENERGY"];  // RHF & UHF only
    Process::environment.globals["MP3 SINGLES ENERGY"] = 0.0;  // RHF & UHF only
    if (reference_ == "UNRESTRICTED") {
        Process::environment.globals["MP3 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB + 2.0 * (Emp3AA + Emp3BB - Emp2AA - Emp2BB);
        Process::environment.globals["MP3 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB + 2.0 * (Emp3AB - Emp2AB);
    }

    Emp3L = Emp3;

    /* updates the wavefunction for checkpointing */
    energy_ = Emp3;
    name_ = "CD-MP2.5";

}  // end mp2_5_manager_cd

//======================================================================
//             OLCCD Manager
//======================================================================
void DFOCC::olccd_manager_cd() {
    // do_cd = "TRUE";
    time4grad = 0;         // means I will not compute the gradient
    mo_optimized = 0;      // means MOs are not optimized
    orbs_already_opt = 0;  // means orbitals are not optimized yet.
    orbs_already_sc = 0;   // means orbitals are not semicanonical yet.

    timer_on("CD Integrals");
    cd_ints();
    trans_cd();
    timer_off("CD Integrals");

    // memalloc for density intermediates
    Jc = std::make_shared<Tensor1d>("DF_BASIS_SCF J_Q", nQ_ref);
    g1Qc = std::make_shared<Tensor1d>("DF_BASIS_SCF G1_Q", nQ_ref);
    g1Qt = std::make_shared<Tensor1d>("DF_BASIS_SCF G1t_Q", nQ_ref);
    g1Qp = std::make_shared<Tensor1d>("DF_BASIS_SCF G1p_Q", nQ_ref);
    g1Q = std::make_shared<Tensor1d>("DF_BASIS_CC G1_Q", nQ);
    g1Qt2 = std::make_shared<Tensor1d>("DF_BASIS_CC G1t_Q", nQ);

    // avaliable mem
    memory = Process::environment.get_memory();
    memory_mb = (double)memory / (1024.0 * 1024.0);
    outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

    // memory requirements
    // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
    cost_df = 0.0;
    cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
    cost_df *= nQ;
    cost_df /= 1024.0 * 1024.0;
    cost_df *= sizeof(double);
    if (reference_ == "RESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
    else if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0 * cost_df);

    // Cost of Integral transform for B(Q,ab)
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)nQ * (long long int)nso2_;
    cost_ampAA += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)nQ * (long long int)nso_ * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

    // Mem for amplitudes
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    cost_3amp = 3.0 * cost_ampAA;
    cost_4amp = 4.0 * cost_ampAA;
    cost_5amp = 5.0 * cost_ampAA;

    if ((cost_4amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
        outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp + cost_df);
        nincore_amp = 4;
        t2_incore = true;
        df_ints_incore = true;
    } else if ((cost_3amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        // outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else if (cost_3amp < memory_mb && cost_df < memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else {
        outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
        outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp + cost_df - memory_mb);
        throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
    }

    // W_abef term
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA += 2.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    double cost_ampAA2 = 0.0;
    cost_ampAA2 = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += 3.0 * (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 /= 1024.0 * 1024.0;
    cost_ampAA2 *= sizeof(double);
    cost_amp = MAX0(cost_ampAA, cost_ampAA2);
    outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

    // Fock
    fock();

    // QCHF
    if (qchf_ == "TRUE") qchf();

    // ROHF REF
    if (reference == "ROHF") {
        t1A = std::make_shared<Tensor2d>("T1_1 <I|A>", naoccA, navirA);
        t1B = std::make_shared<Tensor2d>("T1_1 <i|a>", naoccB, navirB);
        t1_1st_sc();
    }
    lccd_t2_1st_sc();
    Elccd = Emp2;
    ElccdL = Emp2;
    EcorrL = Emp2 - Escf;
    ElccdL_old = Emp2;

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (Canonical CD-MP2)... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    variables_["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;
    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        variables_["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    variables_["MP2 SINGLES ENERGY"] = Emp2_t1;

    // Malloc for PDMs
    gQt = std::make_shared<Tensor1d>("CCD PDM G_Qt", nQ);
    if (reference_ == "RESTRICTED") {
        G1c_ov = std::make_shared<Tensor2d>("Correlation OPDM <O|V>", noccA, nvirA);
        G1c_vo = std::make_shared<Tensor2d>("Correlation OPDM <V|O>", nvirA, noccA);
    } else if (reference_ == "UNRESTRICTED") {
        G1c_ovA = std::make_shared<Tensor2d>("Correlation OPDM <O|V>", noccA, nvirA);
        G1c_ovB = std::make_shared<Tensor2d>("Correlation OPDM <o|v>", noccB, nvirB);
        G1c_voA = std::make_shared<Tensor2d>("Correlation OPDM <V|O>", nvirA, noccA);
        G1c_voB = std::make_shared<Tensor2d>("Correlation OPDM <v|o>", nvirB, noccB);
    }

    lccd_pdm_3index_intr();
    omp3_opdm();
    olccd_tpdm();
    sep_tpdm_cc();
    gfock_cc_vo();
    gfock_cc_ov();
    gfock_cc_oo();
    gfock_cc_vv();
    idp();
    mograd();
    occ_iterations();
    Elccd = ElccdL;

    // main if
    if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod) {
        orbs_already_opt = 1;
        mo_optimized = 1;
        if (conver == 1)
            outfile->Printf("\n\tOrbitals are optimized now.\n");
        else if (conver == 0) {
            outfile->Printf("\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
        }
        outfile->Printf("\tTransforming MOs to the semicanonical basis... \n");
        semi_canonic();
        outfile->Printf("\tSwitching to the standard CD-LCCD computation... \n");
        trans_cd();
        fock();
        ref_energy();
        lccd_iterations();
        conver = 1;
        if (dertype == "FIRST") {
            lccd_pdm_3index_intr();
            omp3_opdm();
            olccd_tpdm();
            sep_tpdm_cc();
            gfock_cc_vo();
            gfock_cc_ov();
            gfock_cc_oo();
            gfock_cc_vv();
        }
    }  // end main if

    else if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod && regularization == "TRUE") {
        outfile->Printf("\tOrbital gradient converged, but energy did not... \n");
        outfile->Printf("\tA tighter rms_mograd_convergence tolerance is recommended... \n");
        throw PSIEXCEPTION("A tighter rms_mograd_convergence tolerance is recommended.");
    }

    if (conver == 1) {
        if (orbs_already_opt == 1) ElccdL = Elccd;

        outfile->Printf("\n");
        outfile->Printf("\tComputing CD-LCCD energy using optimized MOs... \n");
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", ElccdAA);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", ElccdAB);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", ElccdBB);
        outfile->Printf("\tCD-LCCD Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
        outfile->Printf("\tCD-LCCD Total Energy (a.u.)        : %20.14f\n", Elccd);
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\n");

        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\t================ CD-OLCCD FINAL RESULTS =============================== \n");
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tCD-OLCCD Correlation Energy (a.u.) : %20.14f\n", ElccdL - Escf);
        outfile->Printf("\tEcdolccd - Eref (a.u.)             : %20.14f\n", ElccdL - Eref);
        outfile->Printf("\tCD-OLCCD Total Energy (a.u.)       : %20.14f\n", ElccdL);
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\n");

        // Set the global variables with the energies
        variables_["CURRENT ENERGY"] = ElccdL;
        variables_["CURRENT REFERENCE ENERGY"] = Escf;
        variables_["CURRENT CORRELATION ENERGY"] = ElccdL - Escf;
        variables_["OLCCD TOTAL ENERGY"] = ElccdL;
        variables_["OLCCD CORRELATION ENERGY"] = ElccdL - Escf;
        variables_["OLCCD REFERENCE CORRECTION ENERGY"] = Eref - Escf;

        if (reference_ == "UNRESTRICTED") {
            variables_["OLCCD OPPOSITE-SPIN CORRELATION ENERGY"] = ElccdAB;
            variables_["OLCCD SAME-SPIN CORRELATION ENERGY"] = ElccdAA + ElccdBB;
        }

        /* updates the wavefunction for checkpointing */
        energy_ = variables_["CURRENT ENERGY"];
        name_ = "CD-OLCCD";

        // OEPROP
        if (oeprop_ == "TRUE") oeprop();

        // Compute Analytic Gradients
        // if (dertype == "FIRST") dfgrad();

        // Save MOs to wfn
        save_mo_to_wfn();

    }  // end if (conver == 1)

}  // end olccd_manager_cd

//======================================================================
//             LCCD Manager
//======================================================================
void DFOCC::lccd_manager_cd() {
    // do_cd = "TRUE";
    time4grad = 0;     // means i will not compute the gradient
    mo_optimized = 0;  // means MOs are not optimized

    timer_on("CD Integrals");
    cd_ints();
    trans_cd();
    timer_off("CD Integrals");

    if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE" || qchf_ == "TRUE") {
        timer_on("DF REF Integrals");
        df_ref();
        trans_ref();
        timer_off("DF REF Integrals");
        outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
        Jc = SharedTensor1d(new Tensor1d("DF_BASIS_SCF J_Q", nQ_ref));
    }

    // avaliable mem
    memory = Process::environment.get_memory();
    memory_mb = (double)memory / (1024.0 * 1024.0);
    outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

    // memory requirements
    // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
    cost_df = 0.0;
    cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
    cost_df *= nQ;
    cost_df /= 1024.0 * 1024.0;
    cost_df *= sizeof(double);
    if (reference_ == "RESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
    else if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0 * cost_df);

    // Cost of Integral transform for B(Q,ab)
    cost_ampAA = 0.0;
    cost_ampAA = nQ * nso2_;
    cost_ampAA += nQ * navirA * navirA;
    cost_ampAA += nQ * nso_ * navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

    // Mem for amplitudes
    cost_ampAA = 0.0;
    cost_ampAA = naocc2AA * nvir2AA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    cost_3amp = 3.0 * cost_ampAA;
    cost_4amp = 4.0 * cost_ampAA;
    cost_5amp = 5.0 * cost_ampAA;

    if ((cost_4amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
        outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp + cost_df);
        nincore_amp = 4;
        t2_incore = true;
        df_ints_incore = true;
    } else if ((cost_3amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        // outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else if (cost_3amp < memory_mb && cost_df < memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else {
        outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
        outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp + cost_df - memory_mb);
        throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
    }

    // W_abef term
    cost_ampAA = 0.0;
    cost_ampAA = naoccA * naoccA * navirA * navirA;
    cost_ampAA += 2.0 * nQ * navirA * navirA;
    cost_ampAA += navirA * navirA * navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    double cost_ampAA2 = 0.0;
    cost_ampAA2 = naoccA * naoccA * navirA * navirA;
    cost_ampAA2 += nQ * navirA * navirA;
    cost_ampAA2 += 3.0 * navirA * navirA * navirA;
    cost_ampAA2 /= 1024.0 * 1024.0;
    cost_ampAA2 *= sizeof(double);
    cost_amp = MAX0(cost_ampAA, cost_ampAA2);
    outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

    // memalloc for density intermediates
    if (qchf_ == "TRUE" || dertype == "FIRST") {
        g1Qc = std::make_shared<Tensor1d>("DF_BASIS_SCF G1_Q", nQ_ref);
        g1Qt = std::make_shared<Tensor1d>("DF_BASIS_SCF G1t_Q", nQ_ref);
        g1Qp = std::make_shared<Tensor1d>("DF_BASIS_SCF G1p_Q", nQ_ref);
        g1Q = std::make_shared<Tensor1d>("DF_BASIS_CC G1_Q", nQ);
        g1Qt2 = std::make_shared<Tensor1d>("DF_BASIS_CC G1t_Q", nQ);
    }

    // QCHF
    if (qchf_ == "TRUE") qchf();

    // Compute MP2 energy
    if (reference == "ROHF") {
        t1A = std::make_shared<Tensor2d>("T1_1 <I|A>", naoccA, navirA);
        t1B = std::make_shared<Tensor2d>("T1_1 <i|a>", naoccB, navirB);
        t1_1st_sc();
    }
    lccd_t2_1st_sc();

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing CD-MP2 energy ... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    variables_["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;

    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

    if (reference_ == "UNRESTRICTED") {
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        variables_["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    variables_["MP2 SINGLES ENERGY"] = Emp2_t1;

    // Perform LCCD iterations
    timer_on("LCCD");
    lccd_iterations();
    timer_off("LCCD");

    outfile->Printf("\n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\t================ LCCD FINAL RESULTS =================================== \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", ElccdAA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", ElccdBB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", ElccdAB);
    outfile->Printf("\tCD-LCCD Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-LCCD Total Energy (a.u.)        : %20.14f\n", Elccd);
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\n");

    variables_["CURRENT ENERGY"] = Elccd;
    variables_["LCCD TOTAL ENERGY"] = Elccd;

    variables_["CURRENT REFERENCE ENERGY"] = Escf;
    variables_["CURRENT CORRELATION ENERGY"] = Elccd - Escf;
    variables_["LCCD CORRELATION ENERGY"] = Elccd - Escf;

    if (reference_ == "UNRESTRICTED") {
        variables_["LCCD OPPOSITE-SPIN CORRELATION ENERGY"] = ElccdAB;
        variables_["LCCD SAME-SPIN CORRELATION ENERGY"] = ElccdAA + ElccdBB;
        variables_["LCCD DOUBLES ENERGY"] = ElccdAB + ElccdAA + ElccdBB;
    }
    else {
        variables_["LCCD DOUBLES ENERGY"] = Ecorr;  // no ROHF
    }
    variables_["LCCD SINGLES ENERGY"] = 0.0;  // no ROHF

    Process::environment.globals["CURRENT ENERGY"] = Elccd;
    Process::environment.globals["LCCD TOTAL ENERGY"] = Elccd;

    Process::environment.globals["CURRENT REFERENCE ENERGY"] = Escf;
    Process::environment.globals["CURRENT CORRELATION ENERGY"] = Elccd - Escf;
    Process::environment.globals["LCCD CORRELATION ENERGY"] = Elccd - Escf;

    if (reference_ == "UNRESTRICTED") {
        Process::environment.globals["LCCD OPPOSITE-SPIN CORRELATION ENERGY"] = ElccdAB;
        Process::environment.globals["LCCD SAME-SPIN CORRELATION ENERGY"] = ElccdAA + ElccdBB;
        Process::environment.globals["LCCD DOUBLES ENERGY"] = ElccdAB + ElccdAA + ElccdBB;
    }
    else {
        Process::environment.globals["LCCD DOUBLES ENERGY"] = Ecorr;  // no ROHF
    }
    Process::environment.globals["LCCD SINGLES ENERGY"] = 0.0;  // no ROHF

    ElccdL = Elccd;

    /* updates the wavefunction for checkpointing */
    energy_ = Elccd;
    name_ = "CD-LCCD";

}  // end lccd_manager_cd


//======================================================================
//             OREMP Manager
//======================================================================
void DFOCC::oremp_manager_cd() {
    // do_cd = "TRUE";
    time4grad = 0;         // means I will not compute the gradient
    mo_optimized = 0;      // means MOs are not optimized
    orbs_already_opt = 0;  // means orbitals are not optimized yet.
    orbs_already_sc = 0;   // means orbitals are not semicanonical yet.

    timer_on("CD Integrals");
    cd_ints();
    trans_cd();
    timer_off("CD Integrals");

    // memalloc for density intermediates
    Jc = std::make_shared<Tensor1d>("DF_BASIS_SCF J_Q", nQ_ref); //CSB J_Q formed from the DF-REF basis, eq. (36)
    g1Qc = std::make_shared<Tensor1d>("DF_BASIS_SCF G1_Q", nQ_ref); //CSB intermediate gamma_Q, eq. (40)
    g1Qt = std::make_shared<Tensor1d>("DF_BASIS_SCF G1t_Q", nQ_ref); //CSB intermediate gamma_Q~, eq. (41)
    g1Qp = std::make_shared<Tensor1d>("DF_BASIS_SCF G1p_Q", nQ_ref); //CSB ???
    g1Q = std::make_shared<Tensor1d>("DF_BASIS_CC G1_Q", nQ);         //CSB ???
    g1Qt2 = std::make_shared<Tensor1d>("DF_BASIS_CC G1t_Q", nQ);      //CSB ???

    // avaliable mem
    memory = Process::environment.get_memory();
    memory_mb = (double)memory / (1024.0 * 1024.0);
    outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

    // memory requirements
    // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
    cost_df = 0.0;
    cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
    cost_df *= nQ;
    cost_df /= 1024.0 * 1024.0;
    cost_df *= sizeof(double);
    if (reference_ == "RESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
    else if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0 * cost_df);

    // Cost of Integral transform for B(Q,ab)
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)nQ * (long long int)nso2_;
    cost_ampAA += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)nQ * (long long int)nso_ * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

    // Mem for amplitudes
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    cost_3amp = 3.0 * cost_ampAA;
    cost_4amp = 4.0 * cost_ampAA;
    cost_5amp = 5.0 * cost_ampAA;

    if ((cost_4amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
        outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp + cost_df);
        nincore_amp = 4;
        t2_incore = true;
        df_ints_incore = true;
    } else if ((cost_3amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        // outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else if (cost_3amp < memory_mb && cost_df < memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else {
        outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
        outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp + cost_df - memory_mb);
        throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
    }

    // W_abef term
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA += 2.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    double cost_ampAA2 = 0.0;
    cost_ampAA2 = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += 3.0 * (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 /= 1024.0 * 1024.0;
    cost_ampAA2 *= sizeof(double);
    cost_amp = MAX0(cost_ampAA, cost_ampAA2);
    outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

    // Fock
    fock();

    // QCHF
    if (qchf_ == "TRUE") qchf();

    // ROHF REF
    if (reference == "ROHF") {
        t1A = std::make_shared<Tensor2d>("T1_1 <I|A>", naoccA, navirA);
        t1B = std::make_shared<Tensor2d>("T1_1 <i|a>", naoccB, navirB);
        t1_1st_sc();
    }
    lccd_t2_1st_sc();
    Eremp = Emp2;
    ErempL = Emp2;
    EcorrL = Emp2 - Escf;
    ElccdL_old = Emp2;

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing CD-MP2 energy using SCF MOs (Canonical CD-MP2)... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tDF-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    variables_["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;
    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;
    if (reference_ == "UNRESTRICTED") {
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        variables_["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    variables_["MP2 SINGLES ENERGY"] = Emp2_t1;

    // Malloc for PDMs
    gQt = std::make_shared<Tensor1d>("CCD PDM G_Qt", nQ);
    if (reference_ == "RESTRICTED") {
        G1c_ov = std::make_shared<Tensor2d>("Correlation OPDM <O|V>", noccA, nvirA);
        G1c_vo = std::make_shared<Tensor2d>("Correlation OPDM <V|O>", nvirA, noccA);
    } else if (reference_ == "UNRESTRICTED") {
        G1c_ovA = std::make_shared<Tensor2d>("Correlation OPDM <O|V>", noccA, nvirA);
        G1c_ovB = std::make_shared<Tensor2d>("Correlation OPDM <o|v>", noccB, nvirB);
        G1c_voA = std::make_shared<Tensor2d>("Correlation OPDM <V|O>", nvirA, noccA);
        G1c_voB = std::make_shared<Tensor2d>("Correlation OPDM <v|o>", nvirB, noccB);
    }

    lccd_pdm_3index_intr(); //CSB form three-index intermediates T_ia^Q, V_ij^Q, V_ij^Q', V_ab^Q, y_ia^Q
    omp3_opdm();            //CSB form the OPDMs; n o action required for REMP/OO-REMP
    oremp_tpdm();
    sep_tpdm_cc();
    gfock_cc_vo();
    gfock_cc_ov();
    gfock_cc_oo();
    gfock_cc_vv();
    idp();
    mograd();
    occ_iterations();
    Eremp = ErempL;

    // main if
    if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod) {
        orbs_already_opt = 1;
        mo_optimized = 1;
        if (conver == 1)
            outfile->Printf("\n\tOrbitals are optimized now.\n");
        else if (conver == 0) {
            outfile->Printf("\n\tMAX MOGRAD did NOT converged, but RMS MOGRAD converged!!!\n");
            outfile->Printf("\tI will consider the present orbitals as optimized.\n");
        }
        outfile->Printf("\tTransforming MOs to the semicanonical basis... \n");
        semi_canonic();
        outfile->Printf("\tSwitching to the standard CD-REMP computation... \n");
        trans_cd();
        fock();
        ref_energy();
        remp_iterations();
        conver = 1;
        if (dertype == "FIRST") {
            lccd_pdm_3index_intr();
            omp3_opdm();
            oremp_tpdm();
            sep_tpdm_cc();
            gfock_cc_vo();
            gfock_cc_ov();
            gfock_cc_oo();
            gfock_cc_vv();
        }
    }  // end main if

    else if (rms_wog <= tol_grad && std::fabs(DE) >= tol_Eod && regularization == "TRUE") {
        outfile->Printf("\tOrbital gradient converged, but energy did not... \n");
        outfile->Printf("\tA tighter rms_mograd_convergence tolerance is recommended... \n");
        throw PSIEXCEPTION("A tighter rms_mograd_convergence tolerance is recommended.");
    }

    if (conver == 1) {
        if (orbs_already_opt == 1) ErempL = Eremp;

        outfile->Printf("\n");
        outfile->Printf("\tComputing CD-REMP energy using optimized MOs... \n");
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", ErempAA);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", ErempAB);
        if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", ErempBB);
        outfile->Printf("\tCD-REMP Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
        outfile->Printf("\tCD-REMP Total Energy (a.u.)        : %20.14f\n", Eremp);
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\n");

        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\t================ CD-OREMP FINAL RESULTS =============================== \n");
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
        outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
        outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
        outfile->Printf("\tCD-OREMP Correlation Energy (a.u.) : %20.14f\n", ErempL - Escf);
        outfile->Printf("\tEcdoremp - Eref (a.u.)             : %20.14f\n", ErempL - Eref);
        outfile->Printf("\tCD-OREMP Total Energy (a.u.)       : %20.14f\n", ErempL);
        outfile->Printf("\t======================================================================= \n");
        outfile->Printf("\n");

        // Set the global variables with the energies
        variables_["CURRENT ENERGY"] = ErempL;
        variables_["CURRENT REFERENCE ENERGY"] = Escf;
        variables_["CURRENT CORRELATION ENERGY"] = ErempL - Escf;
        variables_["OREMP2 TOTAL ENERGY"] = ErempL;
        variables_["OREMP2 CORRELATION ENERGY"] = ErempL - Escf;
        variables_["OREMP2 REFERENCE CORRECTION ENERGY"] = Eref - Escf;

        /* updates the wavefunction for checkpointing */
        energy_ = variables_["CURRENT ENERGY"];
        name_ = "CD-OREMP";

        // OEPROP
        if (oeprop_ == "TRUE") oeprop();

        // Compute Analytic Gradients
        // if (dertype == "FIRST") dfgrad();

        // Save MOs to wfn
        save_mo_to_wfn();

    }  // end if (conver == 1)

}  // end oremp_manager_cd

//======================================================================
//             REMP Manager
//======================================================================
void DFOCC::remp_manager_cd() {
    // do_cd = "TRUE";
    time4grad = 0;     // means i will not compute the gradient
    mo_optimized = 0;  // means MOs are not optimized

    timer_on("CD Integrals");
    cd_ints();
    trans_cd();
    timer_off("CD Integrals");

    if (dertype == "FIRST" || oeprop_ == "TRUE" || ekt_ip_ == "TRUE" || qchf_ == "TRUE") {
        throw PSIEXCEPTION("NO FIRST DERIVATIVES FOR CD-REMP YET");
        timer_on("DF REF Integrals");
        df_ref();
        trans_ref();
        timer_off("DF REF Integrals");
        outfile->Printf("\tNumber of basis functions in the DF-HF basis: %3d\n", nQ_ref);
        Jc = std::make_shared<Tensor1d>("DF_BASIS_SCF J_Q", nQ_ref);
    }

    // avaliable mem
    memory = Process::environment.get_memory();
    memory_mb = (double)memory / (1024.0 * 1024.0);
    outfile->Printf("\n\tAvailable memory                      : %9.2lf MB \n", memory_mb);

    // memory requirements
    // DF-CC B(Q,ab) + B(Q,ia) + B(Q,ij)
    cost_df = 0.0;
    cost_df = (navirA * navirA) + (navirA * naoccA) + (naoccA * naoccA);
    cost_df *= nQ;
    cost_df /= 1024.0 * 1024.0;
    cost_df *= sizeof(double);
    if (reference_ == "RESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", cost_df);
    else if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tMemory requirement for 3-index ints   : %9.2lf MB \n", 2.0 * cost_df);

    // Cost of Integral transform for B(Q,ab)
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)nQ * (long long int)nso2_;
    cost_ampAA += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)nQ * (long long int)nso_ * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    outfile->Printf("\tMemory requirement for DF-CC int trans: %9.2lf MB \n", cost_ampAA);

    // Mem for amplitudes
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naocc2AA * (long long int)nvir2AA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    cost_3amp = 3.0 * cost_ampAA;
    cost_4amp = 4.0 * cost_ampAA;
    cost_5amp = 5.0 * cost_ampAA;

    if ((cost_4amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_4amp);
        outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_4amp + cost_df);
        nincore_amp = 4;
        t2_incore = true;
        df_ints_incore = true;
    } else if ((cost_3amp + cost_df) <= memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        // outfile->Printf("\tTotal memory requirement for DF+CC int: %9.2lf MB \n", cost_3amp+cost_df);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else if (cost_3amp < memory_mb && cost_df < memory_mb) {
        outfile->Printf("\tMemory requirement for CC contractions: %9.2lf MB \n", cost_3amp);
        outfile->Printf("\tWarning: T2 amplitudes will be stored on the disk!\n");
        nincore_amp = 3;
        t2_incore = false;
        df_ints_incore = false;
    } else {
        outfile->Printf("\tWarning: There is NOT enough memory for CC contractions!\n");
        outfile->Printf("\tIncrease memory by                    : %9.2lf MB \n", cost_3amp + cost_df - memory_mb);
        throw PSIEXCEPTION("There is NOT enough memory for CC contractions!");
    }

    // W_abef term
    cost_ampAA = 0.0;
    cost_ampAA = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA += 2.0 * (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA += (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA /= 1024.0 * 1024.0;
    cost_ampAA *= sizeof(double);
    double cost_ampAA2 = 0.0;
    cost_ampAA2 = (long long int)naoccA * (long long int)naoccA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += (long long int)nQ * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 += 3.0 * (long long int)navirA * (long long int)navirA * (long long int)navirA;
    cost_ampAA2 /= 1024.0 * 1024.0;
    cost_ampAA2 *= sizeof(double);
    cost_amp = MAX0(cost_ampAA, cost_ampAA2);
    outfile->Printf("\tMemory requirement for Wabef term     : %9.2lf MB \n", cost_amp);

    // memalloc for density intermediates
    if (qchf_ == "TRUE" || dertype == "FIRST") {
        g1Qc = std::make_shared<Tensor1d>("DF_BASIS_SCF G1_Q", nQ_ref);
        g1Qt = std::make_shared<Tensor1d>("DF_BASIS_SCF G1t_Q", nQ_ref);
        g1Qp = std::make_shared<Tensor1d>("DF_BASIS_SCF G1p_Q", nQ_ref);
        g1Q = std::make_shared<Tensor1d>("DF_BASIS_CC G1_Q", nQ);
        g1Qt2 = std::make_shared<Tensor1d>("DF_BASIS_CC G1t_Q", nQ);
    }

    // QCHF
    if (qchf_ == "TRUE") qchf();

    // Compute MP2 energy
    if (reference == "ROHF") {
        t1A = std::make_shared<Tensor2d>("T1_1 <I|A>", naoccA, navirA);
        t1B = std::make_shared<Tensor2d>("T1_1 <i|a>", naoccB, navirB);
        t1_1st_sc();
    }
    lccd_t2_1st_sc();

    outfile->Printf("\n");
    if (reference == "ROHF")
        outfile->Printf("\tComputing CD-MP2 energy (CD-ROHF-MP2)... \n");
    else
        outfile->Printf("\tComputing CD-MP2 energy ... \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tCD-HF Energy (a.u.)                : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", Emp2AA);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", Emp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", Emp2BB);
    if (reference_ == "UNRESTRICTED")
        outfile->Printf("\tScaled_SS Correlation Energy (a.u.): %20.14f\n", Escsmp2AA + Escsmp2BB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tScaled_OS Correlation Energy (a.u.): %20.14f\n", Escsmp2AB);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS-MP2 Total Energy (a.u.)     : %20.14f\n", Escsmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SOS-MP2 Total Energy (a.u.)     : %20.14f\n", Esosmp2);
    if (reference_ == "UNRESTRICTED") outfile->Printf("\tCD-SCS(N)-MP2 Total Energy (a.u.)  : %20.14f\n", Escsnmp2);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Singles Energy (a.u.)       : %20.14f\n", Emp2_t1);
    if (reference == "ROHF") outfile->Printf("\tCD-MP2 Doubles Energy (a.u.)       : %20.14f\n", Ecorr - Emp2_t1);
    outfile->Printf("\tCD-MP2 Correlation Energy (a.u.)   : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-MP2 Total Energy (a.u.)         : %20.14f\n", Emp2);
    outfile->Printf("\t======================================================================= \n");

    variables_["MP2 TOTAL ENERGY"] = Emp2;
    Process::environment.globals["SCS-MP2 TOTAL ENERGY"] = Escsmp2;
    Process::environment.globals["SOS-MP2 TOTAL ENERGY"] = Esosmp2;
    Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = Escsnmp2;

    variables_["MP2 CORRELATION ENERGY"] = Emp2 - Escf;
    Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] = Escsmp2 - Escf;
    Process::environment.globals["SOS-MP2 CORRELATION ENERGY"] = Esosmp2 - Escf;
    Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = Escsnmp2 - Escf;

    if (reference_ == "UNRESTRICTED") {
        variables_["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = Emp2AB;
        variables_["MP2 SAME-SPIN CORRELATION ENERGY"] = Emp2AA + Emp2BB;
        variables_["MP2 DOUBLES ENERGY"] = Emp2AB + Emp2AA + Emp2BB;
    }
    else {
        variables_["MP2 DOUBLES ENERGY"] = Ecorr - Emp2_t1;
    }
    variables_["MP2 SINGLES ENERGY"] = Emp2_t1;

    // Perform REMP iterations
    timer_on("LCCD");
    remp_iterations();
    timer_off("LCCD");

    outfile->Printf("\n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\t================ REMP2 FINAL RESULTS =================================== \n");
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\tNuclear Repulsion Energy (a.u.)    : %20.14f\n", Enuc);
    outfile->Printf("\tSCF Energy (a.u.)                  : %20.14f\n", Escf);
    outfile->Printf("\tREF Energy (a.u.)                  : %20.14f\n", Eref);
    if (reference_ == "UNRESTRICTED") {
        outfile->Printf("\tAlpha-Alpha Contribution (a.u.)    : %20.14f\n", ErempAA);
        outfile->Printf("\tBeta-Beta Contribution (a.u.)      : %20.14f\n", ErempBB);
        outfile->Printf("\tAlpha-Beta Contribution (a.u.)     : %20.14f\n", ErempAB);
    }
    outfile->Printf("\tCD-REMP Correlation Energy (a.u.)  : %20.14f\n", Ecorr);
    outfile->Printf("\tCD-REMP Total Energy (a.u.)        : %20.14f\n", Eremp);
    outfile->Printf("\t======================================================================= \n");
    outfile->Printf("\n");

    energy_ = Eremp;
    variables_["CURRENT ENERGY"] = Eremp;
    variables_["REMP2 TOTAL ENERGY"] = Eremp;

    variables_["CURRENT REFERENCE ENERGY"] = Escf;
    variables_["CURRENT CORRELATION ENERGY"] = Eremp - Escf;
    variables_["REMP2 CORRELATION ENERGY"] = Eremp - Escf;

    if (reference_ == "UNRESTRICTED") {
        variables_["REMP2 OPPOSITE-SPIN CORRELATION ENERGY"] = ErempAB;
        variables_["REMP2 SAME-SPIN CORRELATION ENERGY"] = ErempAA + ErempBB;
        variables_["REMP2 DOUBLES ENERGY"] = ErempAB + ErempAA + ErempBB;
    }
    else {
        variables_["REMP2 DOUBLES ENERGY"] = Ecorr;  // no ROHF
    }
    variables_["REMP2 SINGLES ENERGY"] = 0.0;  // no ROHF

    ErempL = Eremp;

    /* updates the wavefunction for checkpointing */
    name_ = "CD-REMP";

}  // end remp_manager_cd

}  // namespace dfoccwave
}  // namespace psi
