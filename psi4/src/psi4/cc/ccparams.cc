/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#include "ccparams.h"

#include "psi4/psi4-dec.h"

#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

namespace psi {
namespace cc {

CCParams::CCParams(Options options) {
    wfn = options.get_str("WFN");
    if (wfn == "NONE") throw PsiException("Invalid value of input keyword WFN", __FILE__, __LINE__);

    newtrips = options.get_bool("NEW_TRIPLES");

    if (wfn == "BCCD" || wfn == "BCCD_T") {
        brueckner = 1;
    } else {
        brueckner = 0;
    }

    df = (options.get_str("CC_TYPE") == "DF");

    semicanonical = false;
    ref = Reference::RHF;
    auto junk = options.get_str("REFERENCE");
    if (junk == "RHF") {
        ref = Reference::RHF;
    } else if (junk == "ROHF" &&
               (wfn == "MP2" || wfn == "CCSD_T" || wfn == "CCSD_AT" || wfn == "CC3" || wfn == "EOM_CC3" ||
                wfn == "CC2" || wfn == "EOM_CC2" || wfn == "BCCD" || wfn == "BCCD_T")) {
        ref = Reference::UHF;
        semicanonical = true;
    } else if (junk == "ROHF") {
        ref = Reference::ROHF;
    } else if (junk == "UHF") {
        ref = Reference::UHF;
    } else {
        throw PsiException("Invalid value of input keyword REFERENCE", __FILE__, __LINE__);
    }

    // Allow user to force semicanonical
    if (options["SEMICANONICAL"].has_changed()) {
        semicanonical = options.get_bool("SEMICANONICAL");
        ref = Reference::UHF;
    }

    analyze = options.get_bool("ANALYZE");

    dertype = DerivativeType::NONE;
    junk = options.get_str("DERTYPE");
    if (junk == "NONE") {
        dertype = DerivativeType::NONE;
    } else if (junk == "FIRST") {
        dertype = DerivativeType::FIRST;
    } else if (junk == "RESPONSE") {
        dertype = DerivativeType::RESPONSE; /* linear response */
    } else {
        throw PsiException("Invalid value of input keyword DERTYPE", __FILE__, __LINE__);
    }

    print = options.get_int("PRINT");
    maxiter = options.get_int("MAXITER");
    convergence = options.get_double("R_CONVERGENCE");
    e_convergence = options.get_double("E_CONVERGENCE");
    restart = options.get_bool("RESTART");

    aobasis = options.get_str("AO_BASIS");
    cachelevel = options.get_int("CACHELEVEL");

    cachetype = CacheType::LOW;
    junk = options.get_str("CACHETYPE");
    if (junk == "LOW") {
        cachetype = CacheType::LOW;
    } else if (junk == "LRU") {
        cachetype = CacheType::LRU;
    } else {
        throw PsiException("Error in input: invalid CACHETYPE", __FILE__, __LINE__);
    }

    /* No LOW cacheing yet for UHF references */
    if (ref == Reference::UHF) cachetype = CacheType::LRU;

    nthreads = Process::environment.get_n_threads();
    if (options["CC_NUM_THREADS"].has_changed()) {
        nthreads = options.get_int("CC_NUM_THREADS");
    }

    diis = options.get_bool("DIIS");
    t2_coupled = options.get_bool("T2_COUPLED");
    prop = options.get_str("PROPERTY");
    abcd = options.get_str("ABCD");
    local = options.get_bool("LOCAL");
    local_cutoff = options.get_double("LOCAL_CUTOFF");
    local_method = options.get_str("LOCAL_METHOD");
    local_weakp = options.get_str("LOCAL_WEAKP");

    local_cphf_cutoff = options.get_double("LOCAL_CPHF_CUTOFF");
    local_freeze_core = (options.get_str("FREEZE_CORE") != "FALSE");

    local_pairdef = options.get_str("LOCAL_PAIRDEF");
    if (local && dertype == DerivativeType::RESPONSE) {
        local_pairdef = "RESPONSE";
    } else if (local) {
        local_pairdef = "BP";
    }

    num_amps = options.get_int("NUM_AMPS_PRINT");
    bconv = options.get_double("BRUECKNER_ORBS_R_CONVERGENCE");

    // Tying orbital convergence to the desired e_conv,
    //   particularly important for sane numerical frequencies by energy
    if (options["BRUECKNER_ORBS_R_CONVERGENCE"].has_changed()) {
        bconv = options.get_double("BRUECKNER_ORBS_R_CONVERGENCE");
    } else {
        bconv = 100.0 * e_convergence;
    }

    print_mp2_amps = options.get_bool("MP2_AMPS_PRINT");
    print_pair_energies = options.get_bool("PAIR_ENERGIES_PRINT");
    spinadapt_energies = options.get_bool("SPINADAPT_ENERGIES");
    t3_Ws_incore = options.get_bool("T3_WS_INCORE");

    /* get parameters related to SCS-MP2 or SCS-N-MP2 */
    /* see papers by S. Grimme or J. Platz */
    scsn = options.get_bool("SCSN_MP2");
    scs = options.get_bool("SCS_MP2");
    scscc = options.get_bool("SCS_CCSD");
    scsmp2_scale_os = options.get_double("MP2_OS_SCALE");
    scsmp2_scale_ss = options.get_double("MP2_SS_SCALE");
    /* see paper by T. Takatani*/
    scscc_scale_os = options.get_double("CC_OS_SCALE");
    scscc_scale_ss = options.get_double("CC_SS_SCALE");

    if (options["MP2_OS_SCALE"].has_changed() || options["MP2_SS_SCALE"].has_changed()) {
        scs = true;
    }

    if (options["CC_OS_SCALE"].has_changed() || options["CC_SS_SCALE"].has_changed()) {
        scscc = true;
    }
}

void print_parameters(const CCParams &params, size_t memory) {
    outfile->Printf("\n    Input parameters:\n");
    outfile->Printf("    -----------------\n");
    outfile->Printf("    Wave function   =     %s\n", params.wfn.c_str());

    if (params.semicanonical) {
        outfile->Printf("    Reference wfn   =     ROHF changed to UHF for Semicanonical Orbitals\n");
    } else {
        outfile->Printf("    Reference wfn   =     %s\n",
                       (params.ref == Reference::RHF) ? "RHF" : ((params.ref == Reference::ROHF) ? "ROHF" : "UHF"));
    }
    outfile->Printf("    Brueckner       =     %s\n", params.brueckner ? "Yes" : "No");
    if (params.brueckner) outfile->Printf("    Brueckner conv. =     %3.1e\n", params.bconv);
    outfile->Printf("    Memory [GiB]    =     %.3f\n", memory * 8 / (1024 * 1024 * 1024.0));
    outfile->Printf("    Maxiter         =    %4d\n", params.maxiter);
    outfile->Printf("    R_Convergence   =     %3.1e\n", params.convergence);
    outfile->Printf("    E_Convergence   =     %3.1e\n", params.e_convergence);
    outfile->Printf("    Restart         =     %s\n", params.restart ? "Yes" : "No");
    outfile->Printf("    DIIS            =     %s\n", params.diis ? "Yes" : "No");
    outfile->Printf("    AO Basis        =     %s\n", params.aobasis.c_str());
    outfile->Printf("    ABCD            =     %s\n", params.abcd.c_str());
    outfile->Printf("    Cache Level     =     %1d\n", params.cachelevel);
    outfile->Printf("    Cache Type      =    %4s\n", (params.cachetype == CacheType::LOW) ? "LOW" : "LRU");
    outfile->Printf("    Print Level     =     %1d\n", params.print);
    outfile->Printf("    Num. of threads =     %d\n", params.nthreads);
    outfile->Printf("    # Amps to Print =     %1d\n", params.num_amps);
    outfile->Printf("    Print MP2 Amps? =     %s\n", params.print_mp2_amps ? "Yes" : "No");
    outfile->Printf("    Analyze T2 Amps =     %s\n", params.analyze ? "Yes" : "No");
    outfile->Printf("    Print Pair Ener =     %s\n", params.print_pair_energies ? "Yes" : "No");

    if (params.print_pair_energies) outfile->Printf("    Spinadapt Ener. =     %s\n", params.spinadapt_energies ? "Yes" : "No");
    outfile->Printf("    Local CC        =     %s\n", params.local ? "Yes" : "No");

    if (params.wfn == "CC3" || params.wfn == "EOM_CC3")
        outfile->Printf("    T3 Ws incore    =     %s\n", params.t3_Ws_incore ? "Yes" : "No");

    if (params.local) {
        outfile->Printf("    Local Cutoff       =     %3.1e\n", params.local_cutoff);
        outfile->Printf("    Local Method      =     %s\n", params.local_method.c_str());
        outfile->Printf("    Weak pairs        =     %s\n", params.local_weakp.c_str());
        outfile->Printf("    Local pairs       =     %s\n", params.local_pairdef.c_str());
        outfile->Printf("    Local CPHF cutoff =     %3.1e\n", params.local_cphf_cutoff);
    }
    outfile->Printf("    SCS-MP2         =     %s\n", params.scs ? "True" : "False");
    outfile->Printf("    SCSN-MP2        =     %s\n", params.scsn ? "True" : "False");
    outfile->Printf("    SCS-CCSD        =     %s\n", params.scscc ? "True" : "False");
    if (params.scs) {
        outfile->Printf("    SCS_MP2_OS_SCALE =     %.2f\n", params.scsmp2_scale_os);
        outfile->Printf("    SCS_MP2_SS_SCALE =     %.2f\n", params.scsmp2_scale_ss);
    }
    if (params.scsn) {
        outfile->Printf("    SCSN_MP2_OS_SCALE =     %.2f\n", 0.0);
        outfile->Printf("    SCSN_MP2_SS_SCALE =     %.2f\n", 1.76);
    }
    if (params.scscc) {
        outfile->Printf("    CC_OS_SCALE     =     %.2f\n", params.scscc_scale_os);
        outfile->Printf("    CC_SS_SCALE     =     %.2f\n", params.scscc_scale_ss);
    }

    outfile->Printf("\n");
}

}  // namespace cc
}  // namespace psi
