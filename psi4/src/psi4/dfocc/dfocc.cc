/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#include <fstream>
#include <cmath>
#include "psi4/libqt/qt.h"
#include "dfocc.h"
#include "defines.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"

using namespace psi;

namespace psi {
namespace dfoccwave {

DFOCC::DFOCC(SharedWavefunction ref_wfn, Options &options) : Wavefunction(options) {
    reference_wavefunction_ = ref_wfn;
    shallow_copy(ref_wfn);
    common_init();
}  //

DFOCC::~DFOCC() {}  //

void DFOCC::common_init() {
    print_ = options_.get_int("PRINT");
    if (print_ > 0) options_.print();

    module_ = "dfocc";

    cc_maxiter = options_.get_int("CC_MAXITER");
    mo_maxiter = options_.get_int("MO_MAXITER");
    num_vecs = options_.get_int("MO_DIIS_NUM_VECS");
    cc_maxdiis_ = options_.get_int("CC_DIIS_MAX_VECS");
    cc_mindiis_ = options_.get_int("CC_DIIS_MIN_VECS");
    exp_cutoff = options_.get_int("CUTOFF");
    exp_int_cutoff = options_.get_int("INTEGRAL_CUTOFF");
    pcg_maxiter = options_.get_int("PCG_MAXITER");

    step_max = options_.get_double("MO_STEP_MAX");
    lshift_parameter = options_.get_double("LEVEL_SHIFT");
    os_scale = options_.get_double("MP2_OS_SCALE");
    ss_scale = options_.get_double("MP2_SS_SCALE");
    sos_scale = options_.get_double("MP2_SOS_SCALE");
    sos_scale2 = options_.get_double("MP2_SOS_SCALE2");
    remp_a = options_.get_double("REMP_A");
    // cepa_os_scale_=options_.get_double("CEPA_OS_SCALE");
    // cepa_ss_scale_=options_.get_double("CEPA_SS_SCALE");
    // cepa_sos_scale_=options_.get_double("CEPA_SOS_SCALE");
    e3_scale = options_.get_double("E3_SCALE");
    tol_Eod = options_.get_double("E_CONVERGENCE");
    tol_t2 = options_.get_double("R_CONVERGENCE");
    reg_param = options_.get_double("REG_PARAM");
    tol_ldl = options_.get_double("CHOLESKY_TOLERANCE");

    orth_type = options_.get_str("ORTH_TYPE");
    opt_method = options_.get_str("OPT_METHOD");
    occ_orb_energy = options_.get_str("OCC_ORBS_PRINT");
    natorb = options_.get_str("NAT_ORBS");
    do_scs = options_.get_str("DO_SCS");
    do_sos = options_.get_str("DO_SOS");
    lineq = options_.get_str("LINEQ_SOLVER");
    level_shift = options_.get_str("DO_LEVEL_SHIFT");
    scs_type_ = options_.get_str("SCS_TYPE");
    sos_type_ = options_.get_str("SOS_TYPE");
    dertype = options_.get_str("DERTYPE");
    orb_resp_solver_ = options_.get_str("ORB_RESP_SOLVER");
    ekt_ip_ = options_.get_str("EKT_IP");
    reference = options_.get_str("REFERENCE");
    wfn_type_ = options_.get_str("WFN_TYPE");
    orb_opt_ = options_.get_str("ORB_OPT");
    // conv_tei_type=options_.get_str("CONV_TEI_TYPE");
    pcg_beta_type_ = options_.get_str("PCG_BETA_TYPE");
    regularization = options_.get_str("REGULARIZATION");
    read_scf_3index = options_.get_str("READ_SCF_3INDEX");
    freeze_core_ = options_.get_str("FREEZE_CORE");
    oeprop_ = options_.get_str("OEPROP");
    comput_s2_ = options_.get_str("COMPUT_S2");
    mp2_amp_type_ = options_.get_str("MP2_AMP_TYPE");
    qchf_ = options_.get_str("QCHF");
    cc_lambda_ = options_.get_str("CC_LAMBDA");
    Wabef_type_ = options_.get_str("PPL_TYPE");
    triples_iabc_type_ = options_.get_str("TRIPLES_IABC_TYPE");
    do_cd = options_.get_str("CHOLESKY");

    if ((dertype == "FIRST" || orb_opt_ == "TRUE" || wfn_type_ == "DF-CCSD(AT)") && cc_lambda_ == "FALSE") {
        cc_lambda_ = "TRUE";
        outfile->Printf("\tWarning: cc_lambda option is changed to TRUE\n");
    }

    if (!psio_) {
        throw PSIEXCEPTION("The wavefunction passed in lacks a PSIO object, crashing DFOCC. See GitHub issue #1851.");
    }

    //   Given default conjugate gradient convergence, set the criteria by what should
    //   be necessary to achive the target energy convergence.
    //   This is based solely on standard suite testing to achieve 1e-6 E & G with default convcrit.
    if (options_["PCG_CONVERGENCE"].has_changed()) {
        tol_pcg = options_.get_double("PCG_CONVERGENCE");
    } else {
        tol_pcg = 0.01 * tol_t2;
        outfile->Printf("\tFor this residual convergence, default PCG convergence is: %12.2e\n", tol_pcg);
    }

    //   Given default orbital convergence, set the criteria by what should
    //   be necessary to achieve the target energy convergence.
    //   These formulae are based on experiments and are nothing rigorous.
    //   The printed value from options_.print() will not be accurate
    //   since newly set orbital conv is not written back to options.
    //   We still want these to be the default values, after all!
    if (orb_opt_ == "TRUE") {
        if (options_["RMS_MOGRAD_CONVERGENCE"].has_changed()) {
            tol_grad = options_.get_double("RMS_MOGRAD_CONVERGENCE");
        } else {
            double temp;
            temp = (-0.9 * std::log10(tol_Eod)) - 1.6;
            if (temp < 6.5) {
                temp = 6.5;
            }
            tol_grad = pow(10.0, -temp);
            outfile->Printf("\tFor this energy convergence, default RMS orbital gradient is: %12.2e\n", tol_grad);
        }

        // Determine the MAXIMUM MOGRAD CONVERGENCE
        if (options_["MAX_MOGRAD_CONVERGENCE"].has_changed()) {
            mograd_max = options_.get_double("MAX_MOGRAD_CONVERGENCE");
        } else {
            double temp2;
            temp2 = (-0.8 * std::log10(tol_grad)) - 0.5;
            if (temp2 < 3.0) {
                temp2 = 3.0;
            }
            mograd_max = pow(10.0, -temp2 - 1);
            outfile->Printf("\tFor this energy convergence, default MAX orbital gradient is: %12.2e\n", mograd_max);
        }
    } else if (orb_opt_ == "FALSE") {
        if (options_["R_CONVERGENCE"].has_changed()) {
            tol_t2 = options_.get_double("R_CONVERGENCE");
        } else {
            double temp;
            temp = (-0.9 * std::log10(tol_Eod)) - 1.6;
            if (temp < 6.5) {
                temp = 6.5;
            }
            tol_t2 = pow(10.0, -temp);
            outfile->Printf("\tFor this energy convergence, default residual convergence is: %12.2e\n", tol_t2);
        }
    }  // end orb_opt_

    // Figure out REF
    if (reference == "RHF" || reference == "RKS")
        reference_ = "RESTRICTED";
    else if (reference == "UHF" || reference == "UKS" || reference == "ROHF" || reference == "ROKS")
        reference_ = "UNRESTRICTED";

    // title
    if (wfn_type_ == "DF-OREMP") {
        remp_title();
    } else if (!((wfn_type_ == "DF-CCSD" || wfn_type_ == "DF-CCSD(T)" || wfn_type_ == "DF-OCCD" || wfn_type_ == "DF-OCCD(T)" || wfn_type_ == "DF-OCCD(AT)") && reference_ == "UNRESTRICTED")) {
        title();
    } else {
        title_dfuccsd();
    }


    // Only ROHF-CC energy is available, not the gradients
    if (reference == "ROHF" && orb_opt_ == "FALSE" && dertype == "FIRST") {
        throw PSIEXCEPTION("ROHF DF-CC analytic gradients are not available, UHF DF-CC is recommended.");
    }

    // DIIS
    if (options_.get_str("DO_DIIS") == "TRUE")
        do_diis_ = 1;
    else if (options_.get_str("DO_DIIS") == "FALSE")
        do_diis_ = 0;

    // Figure out HESSIAN TYPE
    if (orb_opt_ == "TRUE") {
        if (options_["HESS_TYPE"].has_changed()) {
            hess_type = options_.get_str("HESS_TYPE");
        } else {
            if (reference_ == "RESTRICTED" && freeze_core_ == "FALSE") {
                hess_type = "HF";
                outfile->Printf("\tMO Hessian type is changed to 'HF'\n");
            } else if (reference_ == "RESTRICTED" && freeze_core_ == "TRUE") {
                hess_type = "APPROX_DIAG";
                outfile->Printf("\tMO Hessian type is changed to 'APPROX_DIAG'\n");
            } else if (reference_ == "UNRESTRICTED") {
                hess_type = "APPROX_DIAG";
                outfile->Printf("\tMO Hessian type is changed to 'APPROX_DIAG'\n");
            }
        }
    }  // end if (orb_opt_ == "TRUE")

    // Regularization
    if (regularization == "TRUE") {
        outfile->Printf("\n\tNOTE: A regularization procedure will be applied to the method.\n");
        outfile->Printf("\tThe regularization parameter is : %12.2f mh\n", reg_param * 1000.0);
    }

    cutoff = pow(10.0, -exp_cutoff);
    int_cutoff_ = pow(10.0, -exp_int_cutoff);
    get_moinfo();
    pair_index();

    // Frozen virtual
    if (nfrzv > 0 && dertype == "FIRST" && wfn_type_ != "DF-CCSD(T)") {
        throw PSIEXCEPTION("Frozen virtual gradients are not available for this method.");
    }
    if (nfrzv > 0 && orb_opt_ == "TRUE") {
        throw PSIEXCEPTION("Frozen virtual approximation is not available for orbital-optimized methods.");
    }

    // common malloc
    common_malloc();

    if (reference_ == "RESTRICTED") {
        outfile->Printf("\tMO spaces... \n\n");
        outfile->Printf("\t FC   OCC   VIR   FV \n");
        outfile->Printf("\t----------------------\n");
        outfile->Printf("\t%3d  %3d   %3d  %3d\n", nfrzc, naoccA, navirA, nfrzv);
    }  // end if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
        outfile->Printf("\tMO spaces... \n\n");
        outfile->Printf("\t FC   AOCC   BOCC  AVIR   BVIR   FV \n");
        outfile->Printf("\t------------------------------------------\n");
        outfile->Printf("\t%3d   %3d   %3d   %3d    %3d   %3d\n", nfrzc, naoccA, naoccB, navirA, navirB, nfrzv);
    }  // else if (reference_ == "UNRESTRICTED")
}  // end common_init

void DFOCC::title() {
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
    if (wfn_type_ == "DF-OMP2" && orb_opt_ == "TRUE" && do_cd == "FALSE")
        outfile->Printf("                      DF-OMP2 (DF-OO-MP2)   \n");
    else if (wfn_type_ == "DF-OMP2" && orb_opt_ == "FALSE" && do_cd == "FALSE")
        outfile->Printf("                       DF-MP2   \n");
    else if (wfn_type_ == "DF-CCSD" && do_cd == "FALSE")
        outfile->Printf("                       DF-CCSD   \n");
    else if (wfn_type_ == "DF-CCSD(T)" && do_cd == "FALSE")
        outfile->Printf("                       DF-CCSD   \n");
    else if (wfn_type_ == "DF-CCSD(AT)" && do_cd == "FALSE")
        outfile->Printf("                       DF-CCSD   \n");
    else if (wfn_type_ == "DF-CCD" && do_cd == "FALSE")
        outfile->Printf("                       DF-CCD   \n");
    else if (wfn_type_ == "DF-OMP3" && orb_opt_ == "TRUE" && do_cd == "FALSE")
        outfile->Printf("                     DF-OMP3 (DF-OO-MP3)   \n");
    else if (wfn_type_ == "DF-OMP3" && orb_opt_ == "FALSE" && do_cd == "FALSE")
        outfile->Printf("                     DF-MP3   \n");
    else if (wfn_type_ == "DF-OLCCD" && orb_opt_ == "TRUE" && do_cd == "FALSE")
        outfile->Printf("                     DF-OLCCD (DF-OO-LCCD)   \n");
    else if (wfn_type_ == "DF-OLCCD" && orb_opt_ == "FALSE" && do_cd == "FALSE")
        outfile->Printf("                       DF-LCCD   \n");
    else if (wfn_type_ == "DF-OMP2.5" && orb_opt_ == "TRUE" && do_cd == "FALSE")
        outfile->Printf("                    DF-OMP2.5 (DF-OO-MP2.5)   \n");
    else if (wfn_type_ == "DF-OMP2.5" && orb_opt_ == "FALSE" && do_cd == "FALSE")
        outfile->Printf("                    DF-MP2.5  \n");
    else if (wfn_type_ == "DF-OMP2" && orb_opt_ == "TRUE" && do_cd == "TRUE")
        outfile->Printf("                      CD-OMP2 (CD-OO-MP2)   \n");
    else if (wfn_type_ == "DF-OMP2" && orb_opt_ == "FALSE" && do_cd == "TRUE")
        outfile->Printf("                       CD-MP2   \n");
    else if (wfn_type_ == "DF-CCSD" && do_cd == "TRUE")
        outfile->Printf("                       CD-CCSD   \n");
    else if (wfn_type_ == "DF-CCSD(T)" && do_cd == "TRUE")
        outfile->Printf("                       CD-CCSD   \n");
    else if (wfn_type_ == "DF-CCSD(AT)" && do_cd == "TRUE")
        outfile->Printf("                       CD-CCSD   \n");
    else if (wfn_type_ == "DF-CCD" && do_cd == "TRUE")
        outfile->Printf("                       CD-CCD   \n");
    else if (wfn_type_ == "DF-OMP3" && orb_opt_ == "TRUE" && do_cd == "TRUE")
        outfile->Printf("                    CD-OMP3 (CD-OO-MP3)   \n");
    else if (wfn_type_ == "DF-OMP3" && orb_opt_ == "FALSE" && do_cd == "TRUE")
        outfile->Printf("                    CD-MP3   \n");
    else if (wfn_type_ == "DF-OMP2.5" && orb_opt_ == "TRUE" && do_cd == "TRUE")
        outfile->Printf("                   CD-OMP2.5 (CD-OO-MP2.5)   \n");
    else if (wfn_type_ == "DF-OMP2.5" && orb_opt_ == "FALSE" && do_cd == "TRUE")
        outfile->Printf("                    CD-MP2.5   \n");
    else if (wfn_type_ == "DF-OLCCD" && orb_opt_ == "TRUE" && do_cd == "TRUE")
        outfile->Printf("                    CD-OLCCD (CD-OO-LCCD)   \n");
    else if (wfn_type_ == "DF-OLCCD" && orb_opt_ == "FALSE" && do_cd == "TRUE")
        outfile->Printf("                    CD-LCCD   \n");
    else if (wfn_type_ == "QCHF")
        outfile->Printf("                      QCHF   \n");
    outfile->Printf("              Program Written by Ugur Bozkaya\n");
    outfile->Printf("              Latest Revision September 9, 2017\n");
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");

}  //

void DFOCC::lambda_title() {
    if (reference_ == "RESTRICTED") {
        outfile->Printf("\n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf("\n");
        if (wfn_type_ == "DF-CCSD")
            outfile->Printf("                       DF-CCSD-Lambda   \n");
        else if (wfn_type_ == "DF-CCSD(T)" || wfn_type_ == "DF-CCSD(AT)")
            outfile->Printf("                       DF-CCSD(T)-Lambda   \n");
        else if (wfn_type_ == "DF-CCD")
            outfile->Printf("                       DF-CCD-Lambda   \n");
        outfile->Printf("              Program Written by Ugur Bozkaya\n");
        outfile->Printf("              Latest Revision Sep 9, 2020\n");
        outfile->Printf("\n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf("\n");
    }
    else if (reference_ == "UNRESTRICTED") {
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf("\n");
        if (wfn_type_ == "DF-CCSD") {
            outfile->Printf("                              DF-CCSD-Lambda\n") ;
            outfile->Printf("              Program Written by Ugur Bozkaya and Asli Unal\n") ;
            outfile->Printf("              Latest Revision Sep 9, 2020\n") ;
        }
        outfile->Printf("\n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf(" ============================================================================== \n");
        outfile->Printf("\n");
    }
}  //

void DFOCC::pt_title() {
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
    if (wfn_type_ == "DF-CCSD(T)")
        outfile->Printf("                       DF-CCSD(T)   \n");
    else if (wfn_type_ == "DF-CCD(T)")
        outfile->Printf("                       DF-CCD(T)   \n");
    outfile->Printf("              Program Written by Ugur Bozkaya\n");
    outfile->Printf("              Latest Revision April 16, 2017\n");
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
}  //

void DFOCC::pat_title() {
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
    if (wfn_type_ == "DF-CCSD(AT)")
        outfile->Printf("                       DF-CCSD(AT)    \n");
    else if (wfn_type_ == "DF-CCD(AT)")
        outfile->Printf("                       DF-CCD(AT)  \n");
    outfile->Printf("              Program Written by Ugur Bozkaya\n");
    outfile->Printf("              Latest Revision September 9, 2015\n");
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
}  //

void DFOCC::pdm_title() {
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
    outfile->Printf("                         DFPDM   \n");
    outfile->Printf("              Particle Density Matrix Code   \n");
    outfile->Printf("               for Density-Fitted Methods       \n");
    if (reference_ == "RESTRICTED") {
        outfile->Printf("                   by Ugur Bozkaya\n");
    }
    else if (reference_ == "UNRESTRICTED") {
        if (wfn_type_ == "DF-CCSD") {
            outfile->Printf("               by Ugur Bozkaya and Asli Unal\n");
        }
        else if (wfn_type_ == "DF-CCSD(T)") {
            outfile->Printf("               by Ugur Bozkaya, Asli Unal & Yavuz Alagoz\n");
        }
        else {
            outfile->Printf("                   by Ugur Bozkaya\n");
        }
    }
    outfile->Printf("              Latest Revision Sep 5, 2020\n");
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");

}  //

double DFOCC::compute_energy() {
    // Call the appropriate manager
    // do_cd = "FALSE";
    nincore_amp = 3;
    if (wfn_type_ == "DF-OMP2" && orb_opt_ == "TRUE" && do_cd == "FALSE")
        omp2_manager();
    else if (wfn_type_ == "DF-OMP2" && orb_opt_ == "FALSE" && do_cd == "FALSE")
        mp2_manager();
    else if (wfn_type_ == "DF-CCSD" && do_cd == "FALSE")
        ccsd_manager();
    else if (wfn_type_ == "DF-CCSD(T)" && do_cd == "FALSE")
        ccsd_t_manager();
    else if (wfn_type_ == "DF-CCSD(AT)" && do_cd == "FALSE")
        ccsdl_t_manager();
    else if (wfn_type_ == "DF-CCD" && do_cd == "FALSE")
        ccd_manager();
    else if (wfn_type_ == "DF-OMP3" && orb_opt_ == "TRUE" && do_cd == "FALSE")
        omp3_manager();
    else if (wfn_type_ == "DF-OMP3" && orb_opt_ == "FALSE" && do_cd == "FALSE")
        mp3_manager();
    else if (wfn_type_ == "DF-OMP2.5" && orb_opt_ == "TRUE" && do_cd == "FALSE")
        omp2_5_manager();
    else if (wfn_type_ == "DF-OMP2.5" && orb_opt_ == "FALSE" && do_cd == "FALSE")
        mp2_5_manager();
    else if (wfn_type_ == "DF-OLCCD" && orb_opt_ == "TRUE" && do_cd == "FALSE")
        olccd_manager();
    else if (wfn_type_ == "DF-OLCCD" && orb_opt_ == "FALSE" && do_cd == "FALSE")
        lccd_manager();
    else if (wfn_type_ == "DF-OREMP" && orb_opt_ == "TRUE" && do_cd == "FALSE")
        oremp_manager();
    else if (wfn_type_ == "DF-OREMP" && orb_opt_ == "FALSE" && do_cd == "FALSE")
        remp_manager();
    else if (wfn_type_ == "DF-OMP2" && orb_opt_ == "TRUE" && do_cd == "TRUE")
        cd_omp2_manager();
    else if (wfn_type_ == "DF-OMP2" && orb_opt_ == "FALSE" && do_cd == "TRUE")
        cd_mp2_manager();
    else if (wfn_type_ == "DF-CCSD" && do_cd == "TRUE")
        ccsd_manager_cd();
    else if (wfn_type_ == "DF-CCD" && do_cd == "TRUE")
        ccd_manager_cd();
    else if (wfn_type_ == "DF-CCSD(T)" && do_cd == "TRUE")
        ccsd_t_manager_cd();
    else if (wfn_type_ == "DF-CCSD(AT)" && do_cd == "TRUE")
        ccsdl_t_manager_cd();
    else if (wfn_type_ == "DF-OMP3" && orb_opt_ == "TRUE" && do_cd == "TRUE")
        omp3_manager_cd();
    else if (wfn_type_ == "DF-OMP3" && orb_opt_ == "FALSE" && do_cd == "TRUE")
        mp3_manager_cd();
    else if (wfn_type_ == "DF-OMP2.5" && orb_opt_ == "TRUE" && do_cd == "TRUE")
        omp2_5_manager_cd();
    else if (wfn_type_ == "DF-OMP2.5" && orb_opt_ == "FALSE" && do_cd == "TRUE")
        mp2_5_manager_cd();
    else if (wfn_type_ == "DF-OLCCD" && orb_opt_ == "TRUE" && do_cd == "TRUE")
        olccd_manager_cd();
    else if (wfn_type_ == "DF-OLCCD" && orb_opt_ == "FALSE" && do_cd == "TRUE")
        lccd_manager_cd();
    else if (wfn_type_ == "DF-OREMP" && orb_opt_ == "TRUE" && do_cd == "TRUE")
        oremp_manager_cd();
    else if (wfn_type_ == "DF-OREMP" && orb_opt_ == "FALSE" && do_cd == "TRUE")
        remp_manager_cd();
    else if (wfn_type_ == "QCHF")
        qchf_manager();
    else {
        outfile->Printf("WFN_TYPE %s is not recognized\n", wfn_type_.c_str());
        throw PSIEXCEPTION("Unrecognized WFN_TYPE!");
    }

    if (wfn_type_ == "DF-OMP2")
        Etotal = Emp2L;
    else if (wfn_type_ == "DF-CCSD")
        Etotal = Eccsd;
    else if (wfn_type_ == "DF-CCSD(T)")
        Etotal = Eccsd_t;
    else if (wfn_type_ == "DF-CCSD(AT)")
        Etotal = Eccsd_at;
    else if (wfn_type_ == "DF-CCD")
        Etotal = Eccd;
    else if (wfn_type_ == "DF-OMP3")
        Etotal = Emp3L;
    else if (wfn_type_ == "DF-OMP2.5")
        Etotal = Emp3L;
    else if (wfn_type_ == "DF-OLCCD")
        Etotal = ElccdL;
    else if (wfn_type_ == "DF-OREMP")
        Etotal = ElccdL;        //CSB this might need adjustment // FIXME (behnle#1#14.06.2021): maybe replace by an own variable
    else if (wfn_type_ == "QCHF")
        Etotal = Eref;

    return Etotal;

}  // end of compute_energy

void DFOCC::remp_title() {
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
    if (wfn_type_ == "DF-OREMP" && orb_opt_ == "TRUE" && do_cd == "FALSE")
        outfile->Printf("                     DF-OREMP   \n");
    else if (wfn_type_ == "DF-OREMP" && orb_opt_ == "FALSE" && do_cd == "FALSE")
        outfile->Printf("                       DF-REMP   \n");
    else if (wfn_type_ == "DF-OREMP" && orb_opt_ == "TRUE" && do_cd == "TRUE")
        outfile->Printf("                    CD-OREMP   \n");
    else if (wfn_type_ == "DF-OREMP" && orb_opt_ == "FALSE" && do_cd == "TRUE")
        outfile->Printf("                    CD-REMP   \n");
    outfile->Printf("              Program Written by Stefan Behnle\n");
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
}

void DFOCC::common_malloc() {

    if (reference_ == "RESTRICTED") {
        // Memory allocation
        HmoA = std::make_shared<Tensor2d>("MO-basis alpha one-electron ints", nmo_, nmo_);
        FijA = std::make_shared<Tensor2d>("Fint <I|J>", naoccA, naoccA);
        FabA = std::make_shared<Tensor2d>("Fint <A|B>", navirA, navirA);
        HooA = std::make_shared<Tensor2d>("OEI <O|O>", noccA, noccA);
        HovA = std::make_shared<Tensor2d>("OEI <O|V>", noccA, nvirA);
        HvoA = std::make_shared<Tensor2d>("OEI <V|O>", nvirA, noccA);
        HvvA = std::make_shared<Tensor2d>("OEI <V|V>", nvirA, nvirA);
        eigooA = std::make_shared<Tensor1d>("epsilon <I|J>", naoccA);
        eigvvA = std::make_shared<Tensor1d>("epsilon <A|B>", navirA);

        // if we need PDMs
        if (orb_opt_ == "TRUE" || dertype != "NONE" || oeprop_ == "TRUE" || qchf_ == "TRUE" || cc_lambda_ == "TRUE" ||
            ekt_ip_ == "TRUE") {
            GijA = std::make_shared<Tensor2d>("G Intermediate <I|J>", naoccA, naoccA);
            GabA = std::make_shared<Tensor2d>("G Intermediate <A|B>", navirA, navirA);
            G1c_oo = std::make_shared<Tensor2d>("Correlation OPDM <O|O>", noccA, noccA);
            G1c_vv = std::make_shared<Tensor2d>("Correlation OPDM <V|V>", nvirA, nvirA);
            G1 = std::make_shared<Tensor2d>("MO-basis OPDM", nmo_, nmo_);
            G1ao = std::make_shared<Tensor2d>("AO-basis OPDM", nso_, nso_);
            G1c = std::make_shared<Tensor2d>("MO-basis correlation OPDM", nmo_, nmo_);
            GF = std::make_shared<Tensor2d>("MO-basis GFM", nmo_, nmo_);
            GFao = std::make_shared<Tensor2d>("AO-basis GFM", nso_, nso_);
            GFoo = std::make_shared<Tensor2d>("MO-basis GFM <O|O>", noccA, noccA);
            GFvo = std::make_shared<Tensor2d>("MO-basis GFM <V|O>", nvirA, noccA);
            GFov = std::make_shared<Tensor2d>("MO-basis GFM <O|V>", noccA, nvirA);
            GFvv = std::make_shared<Tensor2d>("MO-basis GFM <V|V>", nvirA, nvirA);
            GFtvv = std::make_shared<Tensor2d>("MO-basis Complementary GFM <V|V>", nvirA, nvirA);
            WorbA = std::make_shared<Tensor2d>("MO-basis alpha MO gradient", nmo_, nmo_);
            UorbA = std::make_shared<Tensor2d>("Alpha MO rotation matrix", nmo_, nmo_);
            KorbA = std::make_shared<Tensor2d>("Alpha K MO rotation parameters matrix", nmo_, nmo_);
            KsqrA = std::make_shared<Tensor2d>("Alpha K^2 MO rotation parameters matrix", nmo_, nmo_);
            AvoA = std::make_shared<Tensor2d>("Diagonal MO Hessian <V|O>", nvirA, noccA);
            if (orb_opt_ == "FALSE") {
                WvoA = std::make_shared<Tensor2d>("Effective MO gradient <V|O>", nvirA, noccA);
            }
            if (nfrzc > 0) AooA = std::make_shared<Tensor2d>("Diagonal MO Hessian <I|FC>", naoccA, nfrzc);
        }
    }  // end if (reference_ == "RESTRICTED")

    else if (reference_ == "UNRESTRICTED") {
        // Memory allocation
        HmoA = std::make_shared<Tensor2d>("MO-basis alpha one-electron ints", nmo_, nmo_);
        HmoB = std::make_shared<Tensor2d>("MO-basis beta one-electron ints", nmo_, nmo_);
        FijA = std::make_shared<Tensor2d>("Fint <I|J>", naoccA, naoccA);
        FijB = std::make_shared<Tensor2d>("Fint <i|j>", naoccB, naoccB);
        FabA = std::make_shared<Tensor2d>("Fint <A|B>", navirA, navirA);
        FabB = std::make_shared<Tensor2d>("Fint <a|b>", navirB, navirB);
        HooA = std::make_shared<Tensor2d>("OEI <O|O>", noccA, noccA);
        HooB = std::make_shared<Tensor2d>("OEI <o|o>", noccB, noccB);
        HovA = std::make_shared<Tensor2d>("OEI <O|V>", noccA, nvirA);
        HovB = std::make_shared<Tensor2d>("OEI <o|v>", noccB, nvirB);
        HvoA = std::make_shared<Tensor2d>("OEI <V|O>", nvirA, noccA);
        HvoB = std::make_shared<Tensor2d>("OEI <v|o>", nvirB, noccB);
        HvvA = std::make_shared<Tensor2d>("OEI <V|V>", nvirA, nvirA);
        HvvB = std::make_shared<Tensor2d>("OEI <v|v>", nvirB, nvirB);
        eigooA = std::make_shared<Tensor1d>("epsilon <I|J>", naoccA);
        eigooB = std::make_shared<Tensor1d>("epsilon <i|j>", naoccB);
        eigvvA = std::make_shared<Tensor1d>("epsilon <A|B>", navirA);
        eigvvB = std::make_shared<Tensor1d>("epsilon <a|b>", navirB);

        // if we need PDMs
        if (orb_opt_ == "TRUE" || dertype != "NONE" || oeprop_ == "TRUE" || qchf_ == "TRUE" || ekt_ip_ == "TRUE") {
            GijA = std::make_shared<Tensor2d>("G Intermediate <I|J>", naoccA, naoccA);
            GijB = std::make_shared<Tensor2d>("G Intermediate <i|j>", naoccB, naoccB);
            GabA = std::make_shared<Tensor2d>("G Intermediate <A|B>", navirA, navirA);
            GabB = std::make_shared<Tensor2d>("G Intermediate <a|b>", navirB, navirB);
            G1c_ooA = std::make_shared<Tensor2d>("Correlation OPDM <O|O>", noccA, noccA);
            G1c_ooB = std::make_shared<Tensor2d>("Correlation OPDM <o|o>", noccB, noccB);
            G1c_vvA = std::make_shared<Tensor2d>("Correlation OPDM <V|V>", nvirA, nvirA);
            G1c_vvB = std::make_shared<Tensor2d>("Correlation OPDM <v|v>", nvirB, nvirB);
            G1cA = std::make_shared<Tensor2d>("MO-basis alpha correlation OPDM", nmo_, nmo_);
            G1cB = std::make_shared<Tensor2d>("MO-basis beta correlation OPDM", nmo_, nmo_);
            G1A = std::make_shared<Tensor2d>("MO-basis alpha OPDM", nmo_, nmo_);
            G1B = std::make_shared<Tensor2d>("MO-basis beta OPDM", nmo_, nmo_);
            G1ao = std::make_shared<Tensor2d>("AO-basis OPDM", nso_, nso_);
            GFA = std::make_shared<Tensor2d>("MO-basis alpha GFM", nmo_, nmo_);
            GFB = std::make_shared<Tensor2d>("MO-basis beta GFM", nmo_, nmo_);
            GFao = std::make_shared<Tensor2d>("AO-basis GFM", nso_, nso_);
            GFooA = std::make_shared<Tensor2d>("MO-basis GFM <O|O>", noccA, noccA);
            GFooB = std::make_shared<Tensor2d>("MO-basis GFM <o|o>", noccB, noccB);
            GFvoA = std::make_shared<Tensor2d>("MO-basis GFM <V|O>", nvirA, noccA);
            GFvoB = std::make_shared<Tensor2d>("MO-basis GFM <v|o>", nvirB, noccB);
            GFovA = std::make_shared<Tensor2d>("MO-basis GFM <O|V>", noccA, nvirA);
            GFovB = std::make_shared<Tensor2d>("MO-basis GFM <o|v>", noccB, nvirB);
            GFvvA = std::make_shared<Tensor2d>("MO-basis GFM <V|V>", nvirA, nvirA);
            GFvvB = std::make_shared<Tensor2d>("MO-basis GFM <v|v>", nvirB, nvirB);
            GFtvvA = std::make_shared<Tensor2d>("MO-basis Complementary GFM <V|V>", nvirA, nvirA);
            GFtvvB = std::make_shared<Tensor2d>("MO-basis Complementary GFM <v|v>", nvirB, nvirB);
            WorbA = std::make_shared<Tensor2d>("MO-basis alpha MO gradient", nmo_, nmo_);
            WorbB = std::make_shared<Tensor2d>("MO-basis beta MO gradient", nmo_, nmo_);
            UorbA = std::make_shared<Tensor2d>("Alpha MO rotation matrix", nmo_, nmo_);
            UorbB = std::make_shared<Tensor2d>("Beta MO rotation matrix", nmo_, nmo_);
            KorbA = std::make_shared<Tensor2d>("Alpha K MO rotation parameters matrix", nmo_, nmo_);
            KorbB = std::make_shared<Tensor2d>("Beta K MO rotation parameters matrix", nmo_, nmo_);
            KsqrA = std::make_shared<Tensor2d>("Alpha K^2 MO rotation parameters matrix", nmo_, nmo_);
            KsqrB = std::make_shared<Tensor2d>("Beta K^2 MO rotation parameters matrix", nmo_, nmo_);
            AvoA = std::make_shared<Tensor2d>("Diagonal MO Hessian <V|O>", nvirA, noccA);
            AvoB = std::make_shared<Tensor2d>("Diagonal MO Hessian <v|o>", nvirB, noccB);
            if (orb_opt_ == "FALSE") {
                WvoA = std::make_shared<Tensor2d>("Effective MO gradient <V|O>", nvirA, noccA);
                WvoB = std::make_shared<Tensor2d>("Effective MO gradient <v|o>", nvirB, noccB);
            }
            if (nfrzc > 0) {
                AooA = std::make_shared<Tensor2d>("Diagonal MO Hessian <I|FC>", naoccA, nfrzc);
                AooB = std::make_shared<Tensor2d>("Diagonal MO Hessian <i|FC>", naoccB, nfrzc);
            }
        }

        // ROHF-MP2
        if (reference == "ROHF" && wfn_type_ == "DF-OMP2") {
            t1A = std::make_shared<Tensor2d>("T1_1 <I|A>", naoccA, navirA);
            t1B = std::make_shared<Tensor2d>("T1_1 <i|a>", naoccB, navirB);
            GiaA = std::make_shared<Tensor2d>("G Intermediate <I|A>", naoccA, navirA);
            GiaB = std::make_shared<Tensor2d>("G Intermediate <i|a>", naoccB, navirB);
            GaiA = std::make_shared<Tensor2d>("G Intermediate <A|I>", navirA, naoccA);
            GaiB = std::make_shared<Tensor2d>("G Intermediate <a|i>", navirB, naoccB);
            G1c_ovA = std::make_shared<Tensor2d>("Correlation OPDM <O|V>", noccA, nvirA);
            G1c_ovB = std::make_shared<Tensor2d>("Correlation OPDM <o|v>", noccB, nvirB);
            G1c_voA = std::make_shared<Tensor2d>("Correlation OPDM <V|O>", nvirA, noccA);
            G1c_voB = std::make_shared<Tensor2d>("Correlation OPDM <v|o>", nvirB, noccB);
        }
    }  // else if (reference_ == "UNRESTRICTED")

}  // end of common_malloc


void DFOCC::title_dfuccsd()
{
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
    if (wfn_type_ == "DF-CCSD")
        outfile->Printf("                       DF-CCSD   \n");
    else if ((wfn_type_ == "DF-OCCD" || wfn_type_ == "DF-OCCD(T)") && orb_opt_ == "FALSE")
        outfile->Printf("                       DF-CCD   \n");
    else if ((wfn_type_ == "DF-OCCD" || wfn_type_ == "DF-OCCD(T)" || wfn_type_ == "DF-OCCD(AT)") && orb_opt_ == "TRUE")
        outfile->Printf("                       DF-OCCD   \n");
    outfile->Printf("              Program Written by Ugur Bozkaya, Asli Unal & Yavuz Alagoz\n") ;
    outfile->Printf("              Latest Revision Sep 9, 2020\n") ;
    outfile->Printf("\n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf(" ============================================================================== \n");
    outfile->Printf("\n");
}//

}  // namespace dfoccwave
}  // namespace psi
