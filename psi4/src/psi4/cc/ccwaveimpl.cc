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

#include "ccwaveimpl.h"

#include <vector>

#include "psi4/psi4-dec.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

namespace psi {
namespace cc {

std::vector<int> get_pitzer2qt(std::vector<Dimension> &spaces) {
    int nirreps = spaces[0].n();

    Dimension total(nirreps);
    for (int h = 0; h < nirreps; h++)
        for (int i = 0; i < spaces.size(); i++) total[h] += spaces[i][h];
    int nmo = total.sum();

    std::vector<int> order(nmo, 0);

    Dimension offset(nirreps);
    offset[0] = 0;
    for (int h = 1; h < nirreps; h++) offset[h] = offset[h - 1] + total[h - 1];

    int count = 0;

    for (int j = 0; j < spaces.size(); j++)
        for (int h = 0; h < nirreps; h++) {
            int this_offset = offset[h];
            for (int k = 0; k < j; k++) this_offset += spaces[k][h];
            for (int i = 0; i < spaces[j][h]; i++) order[this_offset + i] = count++;
        }

    return order;
}

CCWavefunctionImpl::CCWavefunctionImpl(std::shared_ptr<Wavefunction> wfnobj, Options options) {
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

    nirreps = wfnobj->nirrep();
    nmo = wfnobj->nmo();
    nso = wfnobj->nso();
    nao = wfnobj->basisset()->nao();

    labels = wfnobj->molecule()->irrep_labels();

    enuc = wfnobj->molecule()->nuclear_repulsion_energy(wfnobj->get_dipole_field_strength());
    escf = 0.0;
    if (wfnobj->reference_wavefunction()) {
        escf = wfnobj->reference_wavefunction()->reference_energy();
    } else {
        escf = wfnobj->reference_energy();
    }
    // TODO we probably need to get the PCM polarization energy too?

    orbspi = wfnobj->nmopi();
    sopi = wfnobj->nsopi();
    frdocc = wfnobj->frzcpi();
    fruocc = wfnobj->frzvpi();
    clsdpi = wfnobj->doccpi() - frdocc;
    openpi = wfnobj->soccpi();
    uoccpi = orbspi - clsdpi - openpi - frdocc - fruocc;
    occpi = clsdpi + openpi;
    aoccpi = clsdpi + openpi;
    boccpi = clsdpi;
    virtpi = uoccpi + openpi;
    avirtpi = uoccpi;
    bvirtpi = uoccpi + openpi;
    nvirt = virtpi.sum();

    // Build Pitzer->QT and QT->Pitzer reordering arrays
    std::vector<Dimension> subspaces;
    if (ref == Reference::UHF) {
        subspaces.push_back(frdocc);
        subspaces.push_back(aoccpi);
        subspaces.push_back(avirtpi);
        subspaces.push_back(fruocc);
        pitzer2qt_a = get_pitzer2qt(subspaces);
        qt2pitzer_a.reserve(nmo);
        for (int i = 0; i < nmo; i++) qt2pitzer_a[pitzer2qt_a[i]] = i;
        subspaces.clear();
        subspaces.push_back(frdocc);
        subspaces.push_back(boccpi);
        subspaces.push_back(bvirtpi);
        subspaces.push_back(fruocc);
        pitzer2qt_b = get_pitzer2qt(subspaces);
        qt2pitzer_b.reserve(nmo);
        for (int i = 0; i < nmo; i++) qt2pitzer_b[pitzer2qt_b[i]] = i;
    } else {
        subspaces.push_back(frdocc);
        subspaces.push_back(clsdpi);
        subspaces.push_back(openpi);
        subspaces.push_back(uoccpi);
        subspaces.push_back(fruocc);
        pitzer2qt = get_pitzer2qt(subspaces);
        qt2pitzer.reserve(nmo);
        for (int i = 0; i < nmo; i++) qt2pitzer[pitzer2qt[i]] = i;
    }

    int nfzc = frdocc.sum();
    int nclsd = clsdpi.sum();
    int nopen = openpi.sum();
    int nuocc = uoccpi.sum();
    nactive = nclsd + nopen + nuocc;

    // Build QT->CC and CC->QT reordering arrays
    if (ref == Reference::UHF) {  // UHF/semicanonical
        aocc_off.push_back(0);
        avir_off.push_back(0);
        bocc_off.push_back(0);
        bvir_off.push_back(0);
        int aocount = aoccpi[0];
        int avcount = avirtpi[0];
        int bocount = boccpi[0];
        int bvcount = bvirtpi[0];
        for (int h = 1; h < nirreps; h++) {
            aocc_off.push_back(aocount);
            aocount += aoccpi[h];
            avir_off.push_back(avcount);
            avcount += avirtpi[h];
            bocc_off.push_back(bocount);
            bocount += boccpi[h];
            bvir_off.push_back(bvcount);
            bvcount += bvirtpi[h];
        }

        cc_aocc.assign(nactive, -1);
        cc_bocc.assign(nactive, -1);
        cc_avir.assign(nactive, -1);
        cc_bvir.assign(nactive, -1);
        qt_aocc.assign(nactive, -1);
        qt_bocc.assign(nactive, -1);
        qt_avir.assign(nactive, -1);
        qt_bvir.assign(nactive, -1);
        aocc_sym.assign(nactive, -1);
        bocc_sym.assign(nactive, -1);
        avir_sym.assign(nactive, -1);
        bvir_sym.assign(nactive, -1);
        for (int h = 0, count = 0, offset = 0; h < nirreps; h++) {
            if (h) offset += clsdpi[h - 1] + openpi[h - 1];
            for (int i = 0; i < clsdpi[h] + openpi[h]; i++, count++) {
                cc_aocc[offset + i] = count;
                qt_aocc[count] = nfzc + offset + i;
                aocc_sym[count] = h;
            }
        }
        for (int h = 0, count = 0, offset = 0; h < nirreps; h++) {
            if (h) offset += clsdpi[h - 1];
            for (int i = 0; i < clsdpi[h]; i++, count++) {
                cc_bocc[offset + i] = count;
                qt_bocc[count] = nfzc + offset + i;
                bocc_sym[count] = h;
            }
        }
        for (int h = 0, count = 0, offset = nclsd + nopen; h < nirreps; h++) {
            if (h) offset += uoccpi[h - 1];
            for (int i = 0; i < uoccpi[h]; i++, count++) {
                cc_avir[offset + i] = count;
                qt_avir[count] = nfzc + offset + i;
                avir_sym[count] = h;
            }
        }
        for (int h = 0, count = 0, offset = nclsd; h < nirreps; h++) {
            if (h) offset += uoccpi[h - 1] + openpi[h - 1];
            for (int i = 0; i < uoccpi[h] + openpi[h]; i++, count++) {
                cc_bvir[offset + i] = count;
                qt_bvir[count] = nfzc + offset + i;
                bvir_sym[count] = h;
            }
        }
    } else {  // RHF/ROHF
        occ_off.push_back(0);
        vir_off.push_back(0);
        int ocount = occpi[0];
        int vcount = virtpi[0];
        for (int h = 1; h < nirreps; h++) {
            occ_off.push_back(ocount);
            ocount += occpi[h];
            vir_off.push_back(vcount);
            vcount += virtpi[h];
        }

        cc_occ.assign(nactive, -1);
        cc_vir.assign(nactive, -1);
        qt_occ.assign(nactive, -1);
        qt_vir.assign(nactive, -1);
        occ_sym.assign(nactive, -1);
        vir_sym.assign(nactive, -1);
        for (int h = 0, count = 0, cl_offset = 0, op_offset = nclsd; h < nirreps; h++) {
            if (h) cl_offset += clsdpi[h - 1];
            for (int i = 0; i < clsdpi[h]; i++, count++) {
                cc_occ[cl_offset + i] = count;
                qt_occ[count] = nfzc + cl_offset + i;
                occ_sym[count] = h;
            }
            if (h) op_offset += openpi[h - 1];
            for (int i = 0; i < openpi[h]; i++, count++) {
                cc_occ[op_offset + i] = count;
                qt_occ[count] = nfzc + op_offset + i;
                occ_sym[count] = h;
            }
        }
        for (int h = 0, count = 0, vr_offset = nclsd + nopen, op_offset = nclsd; h < nirreps; h++) {
            if (h) vr_offset += uoccpi[h - 1];
            for (int i = 0; i < uoccpi[h]; i++, count++) {
                cc_vir[vr_offset + i] = count;
                qt_vir[count] = nfzc + vr_offset + i;
                vir_sym[count] = h;
            }
            if (h) op_offset += openpi[h - 1];
            for (int i = 0; i < openpi[h]; i++, count++) {
                cc_vir[op_offset + i] = count;
                qt_vir[count] = nfzc + op_offset + i;
                vir_sym[count] = h;
            }
        }
    }

    // Build sosym array (for AO-basis BT2)
    sosym.reserve(nso);
    for (int h = 0, q = 0; h < nirreps; h++)
        for (int p = 0; p < sopi[h]; p++) sosym[q++] = h;
}

void CCWavefunctionImpl::print_out(size_t memory, std::string out) const {
    auto printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));

    printer->Printf("\n    Input parameters:\n");
    printer->Printf("    -----------------\n");
    printer->Printf("    Wave function   =     %s\n", wfn.c_str());

    if (semicanonical) {
        printer->Printf("    Reference wfn   =     ROHF changed to UHF for Semicanonical Orbitals\n");
    } else {
        printer->Printf("    Reference wfn   =     %s\n",
                        (ref == Reference::RHF) ? "RHF" : ((ref == Reference::ROHF) ? "ROHF" : "UHF"));
    }
    printer->Printf("    Brueckner       =     %s\n", brueckner ? "Yes" : "No");
    if (brueckner) printer->Printf("    Brueckner conv. =     %3.1e\n", bconv);
    printer->Printf("    Memory [GiB]    =     %.3f\n", memory * 8 / (1024 * 1024 * 1024.0));
    printer->Printf("    Maxiter         =    %4d\n", maxiter);
    printer->Printf("    R_Convergence   =     %3.1e\n", convergence);
    printer->Printf("    E_Convergence   =     %3.1e\n", e_convergence);
    printer->Printf("    Restart         =     %s\n", restart ? "Yes" : "No");
    printer->Printf("    DIIS            =     %s\n", diis ? "Yes" : "No");
    printer->Printf("    AO Basis        =     %s\n", aobasis.c_str());
    printer->Printf("    ABCD            =     %s\n", abcd.c_str());
    printer->Printf("    Cache Level     =     %1d\n", cachelevel);
    printer->Printf("    Cache Type      =    %4s\n", (cachetype == CacheType::LOW) ? "LOW" : "LRU");
    printer->Printf("    Print Level     =     %1d\n", print);
    printer->Printf("    Num. of threads =     %d\n", nthreads);
    printer->Printf("    # Amps to Print =     %1d\n", num_amps);
    printer->Printf("    Print MP2 Amps? =     %s\n", print_mp2_amps ? "Yes" : "No");
    printer->Printf("    Analyze T2 Amps =     %s\n", analyze ? "Yes" : "No");
    printer->Printf("    Print Pair Ener =     %s\n", print_pair_energies ? "Yes" : "No");

    if (print_pair_energies)
        printer->Printf("    Spinadapt Ener. =     %s\n", spinadapt_energies ? "Yes" : "No");
    printer->Printf("    Local CC        =     %s\n", local ? "Yes" : "No");

    if (wfn == "CC3" || wfn == "EOM_CC3")
        printer->Printf("    T3 Ws incore    =     %s\n", t3_Ws_incore ? "Yes" : "No");

    if (local) {
        printer->Printf("    Local Cutoff       =     %3.1e\n", local_cutoff);
        printer->Printf("    Local Method      =     %s\n", local_method.c_str());
        printer->Printf("    Weak pairs        =     %s\n", local_weakp.c_str());
        printer->Printf("    Local pairs       =     %s\n", local_pairdef.c_str());
        printer->Printf("    Local CPHF cutoff =     %3.1e\n", local_cphf_cutoff);
    }
    printer->Printf("    SCS-MP2         =     %s\n", scs ? "True" : "False");
    printer->Printf("    SCSN-MP2        =     %s\n", scsn ? "True" : "False");
    printer->Printf("    SCS-CCSD        =     %s\n", scscc ? "True" : "False");
    if (scs) {
        printer->Printf("    SCS_MP2_OS_SCALE =     %.2f\n", scsmp2_scale_os);
        printer->Printf("    SCS_MP2_SS_SCALE =     %.2f\n", scsmp2_scale_ss);
    }
    if (scsn) {
        printer->Printf("    SCSN_MP2_OS_SCALE =     %.2f\n", 0.0);
        printer->Printf("    SCSN_MP2_SS_SCALE =     %.2f\n", 1.76);
    }
    if (scscc) {
        printer->Printf("    CC_OS_SCALE     =     %.2f\n", scscc_scale_os);
        printer->Printf("    CC_SS_SCALE     =     %.2f\n", scscc_scale_ss);
    }

    printer->Printf("\n");
    printer->Printf("\n\tWfn Parameters:\n");
    printer->Printf("\t--------------------\n");
    printer->Printf("\tNumber of irreps     = %d\n", nirreps);
    printer->Printf("\tNumber of MOs        = %d\n", nmo);
    printer->Printf("\tNumber of active MOs = %d\n", nactive);
    printer->Printf("\n");

    printer->Printf("\tIRREP\t# MOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
    printer->Printf("\t-----\t-----\t------\t------\t------\t------\t------\n");
    for (int i = 0; i < nirreps; i++) {
        printer->Printf("\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n", labels[i].c_str(),
                        orbspi[i], frdocc[i], clsdpi[i], openpi[i], uoccpi[i], fruocc[i]);
    }
    printer->Printf("    Nuclear Rep. energy (wfn)     = %20.15f\n", enuc);
    printer->Printf("    SCF energy          (wfn)     = %20.15f\n", escf);
}

}  // namespace cc
}  // namespace psi
