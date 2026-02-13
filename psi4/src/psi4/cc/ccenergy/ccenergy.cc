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

/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
/*
**  CCENERGY: Program to calculate coupled cluster energies.
*/

#include "Params.h"
#include "MOInfo.h"
#include "Local.h"
#include "psi4/cc/ccwave.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/psifiles.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/liboptions/liboptions.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

// Forward declaration to call cctriples
namespace psi {
namespace cctriples {
PsiReturnType cctriples(std::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
}
}  // namespace psi

namespace psi {
namespace ccenergy {

CCEnergyWavefunction::CCEnergyWavefunction(std::shared_ptr<Wavefunction> reference_wavefunction, Options &options)
    : Wavefunction(options) {
    set_reference_wavefunction(reference_wavefunction);
    init();
}

CCEnergyWavefunction::~CCEnergyWavefunction() {}

void CCEnergyWavefunction::init() {
    shallow_copy(reference_wavefunction_);
    module_ = "ccenergy";
    // shallow_copy leaves this wfn's Density as the same object as the ref's.
    // Because ccdensity/ccdensity.cc needs to modify the object at that memory
    // location, make this object's densities point to a different memory locationn
    Da_ = reference_wavefunction_->Da()->clone();
    Db_ = reference_wavefunction_->same_a_b_dens() ? Da_ : reference_wavefunction_->Db()->clone();
}

double CCEnergyWavefunction::compute_energy() {
    int done = 0;
    int **cachelist;

    moinfo_.iter = 0;

    init_io();
    title();

    get_moinfo();
    get_params(options_);

    auto cachefiles = std::vector<int>(PSIO_MAXUNIT);

    if (params_.ref == 2) { /** UHF **/
        cachelist = cacheprep_uhf(params_.cachelev, cachefiles.data());
        std::vector<int *> spaces;
        spaces.push_back(moinfo_.aoccpi);
        spaces.push_back(moinfo_.aocc_sym);
        spaces.push_back(moinfo_.avirtpi);
        spaces.push_back(moinfo_.avir_sym);
        spaces.push_back(moinfo_.boccpi);
        spaces.push_back(moinfo_.bocc_sym);
        spaces.push_back(moinfo_.bvirtpi);
        spaces.push_back(moinfo_.bvir_sym);
        delete[] dpd_list[0];
        dpd_list[0] = new DPD(0, moinfo_.nirreps, params_.memory, 0, cachefiles.data(), cachelist, nullptr, 4, spaces);
        dpd_set_default(0);

        if (params_.df) {
            form_df_ints(options_, cachelist, cachefiles.data());
        } else if (params_.aobasis != "NONE") { /* Set up new DPD's for AO-basis algorithm */
            std::vector<int *> aospaces;
            aospaces.push_back(moinfo_.aoccpi);
            aospaces.push_back(moinfo_.aocc_sym);
            aospaces.push_back(moinfo_.sopi);
            aospaces.push_back(moinfo_.sosym);
            aospaces.push_back(moinfo_.boccpi);
            aospaces.push_back(moinfo_.bocc_sym);
            aospaces.push_back(moinfo_.sopi);
            aospaces.push_back(moinfo_.sosym);
            dpd_init(1, moinfo_.nirreps, params_.memory, 0, cachefiles.data(), cachelist, nullptr, 4, aospaces);
            dpd_set_default(0);
        }

    } else { /** RHF or ROHF **/
        cachelist = cacheprep_rhf(params_.cachelev, cachefiles.data());
        init_priority_list();
        std::vector<int *> spaces;
        spaces.push_back(moinfo_.occpi);
        spaces.push_back(moinfo_.occ_sym);
        spaces.push_back(moinfo_.virtpi);
        spaces.push_back(moinfo_.vir_sym);

        dpd_init(0, moinfo_.nirreps, params_.memory, params_.cachetype, cachefiles.data(), cachelist,
                 cache_priority_list_.data(), 2, spaces);

        if (params_.df) {
            form_df_ints(options_, cachelist, cachefiles.data());
        } else if (params_.aobasis != "NONE") { /* Set up new DPD for AO-basis algorithm */
            std::vector<int *> aospaces;
            aospaces.push_back(moinfo_.occpi);
            aospaces.push_back(moinfo_.occ_sym);
            aospaces.push_back(moinfo_.sopi);
            aospaces.push_back(moinfo_.sosym);
            dpd_init(1, moinfo_.nirreps, params_.memory, 0, cachefiles.data(), cachelist, nullptr, 2, aospaces);
            dpd_set_default(0);
        }
    }

    if ((params_.just_energy) || (params_.just_residuals)) {
        one_step();
        if (params_.ref == 2)
            cachedone_uhf(cachelist);
        else
            cachedone_rhf(cachelist);
        cleanup();
        exit_io();
        return Success;
    }

    if (params_.local) {
        local_init();
        if (local_.weakp == "MP2") lmp2();
    }

    init_amps();

    /* Compute the MP2 energy while we're here */
    if (params_.ref == 0 || params_.ref == 2) {
        moinfo_.emp2 = mp2_energy();
        outfile->Printf("MP2 correlation energy %4.16f\n", moinfo_.emp2);
        psio_write_entry(PSIF_CC_INFO, "MP2 Energy", (char *)&(moinfo_.emp2), sizeof(double));
        Process::environment.globals["MP2 CORRELATION ENERGY"] = moinfo_.emp2;
        Process::environment.globals["MP2 TOTAL ENERGY"] = moinfo_.emp2 + moinfo_.eref;
        Process::environment.globals["MP2 SINGLES ENERGY"] = moinfo_.emp2_s;
        Process::environment.globals["MP2 DOUBLES ENERGY"] = moinfo_.emp2_ss + moinfo_.emp2_os;
        set_scalar_variable("MP2 CORRELATION ENERGY", moinfo_.emp2);
        set_scalar_variable("MP2 TOTAL ENERGY", moinfo_.emp2 + moinfo_.eref);
        set_scalar_variable("MP2 SINGLES ENERGY", moinfo_.emp2_s);
        set_scalar_variable("MP2 DOUBLES ENERGY", moinfo_.emp2_ss + moinfo_.emp2_os);
    }

    if (params_.print_mp2_amps) amp_write();

    tau_build();
    taut_build();
    outfile->Printf("                Solving CC Amplitude Equations\n");
    outfile->Printf("                ------------------------------\n");
    outfile->Printf("  Iter             Energy              RMS        T1Diag      D1Diag    New D1Diag    D2Diag\n");
    outfile->Printf("  ----     ---------------------    ---------   ----------  ----------  ----------   --------\n");
    moinfo_.ecc = energy();

    auto nocc_act = moinfo_.clsdpi.sum();
    std::vector<double> emp2_aa(nocc_act * (nocc_act - 1) /2);
    std::vector<double> emp2_ab(nocc_act * nocc_act);
    pair_energies(emp2_aa, emp2_ab);
    double last_energy = 0;

    moinfo_.t1diag = diagnostic();
    moinfo_.d1diag = d1diag();
    moinfo_.new_d1diag = new_d1diag();

    moinfo_.d2diag = d2diag();
    update();
    checkpoint();
    for (moinfo_.iter = 1; moinfo_.iter <= params_.maxiter; moinfo_.iter++) {
        sort_amps();

        timer_on("F build");
        Fme_build();
        Fae_build();
        Fmi_build();
        if (params_.print & 2) status("F intermediates", "outfile");
        timer_off("F build");

        t1_build();
        if (params_.print & 2) status("T1 amplitudes", "outfile");

        if (params_.wfn == "CC2" || params_.wfn == "EOM_CC2") {
            cc2_Wmnij_build();
            if (params_.print & 2) status("Wmnij", "outfile");

            timer_on("Wmbij build");
            cc2_Wmbij_build();
            if (params_.print & 2) status("Wmbij", "outfile");
            timer_off("Wmbij build");

            timer_on("Wabei build");
            cc2_Wabei_build();
            if (params_.print & 2) status("Wabei", "outfile");
            timer_off("Wabei build");

            timer_on("T2 Build");
            cc2_t2_build();
            if (params_.print & 2) status("T2 amplitudes", "outfile");
            timer_off("T2 Build");
        } else {
            timer_on("Wmbej build");
            Wmbej_build();
            if (params_.print & 2) status("Wmbej", "outfile");
            timer_off("Wmbej build");

            Z_build();
            if (params_.print & 2) status("Z", "outfile");
            Wmnij_build();
            if (params_.print & 2) status("Wmnij", "outfile");

            timer_on("T2 Build");
            t2_build();
            if (params_.print & 2) status("T2 amplitudes", "outfile");
            timer_off("T2 Build");

            if (params_.wfn == "CC3" || params_.wfn == "EOM_CC3") {
                /* step1: build cc3 intermediates, Wabei, Wmnie, Wmbij, Wamef */
                cc3_Wmnij();
                cc3_Wmbij();
                cc3_Wmnie();
                cc3_Wamef();
                cc3_Wabei();

                /* step2: loop over T3's and add contributions to T1 and T2 as you go */
                cc3();
            }
        }

        if (!params_.just_residuals) denom(); /* apply denominators to T1 and T2 */

        if (converged(last_energy - moinfo_.ecc)) {
            done = 1;

            tsave();
            tau_build();
            taut_build();
            last_energy = moinfo_.ecc;
            moinfo_.ecc = energy();
            moinfo_.t1diag = diagnostic();
            moinfo_.d1diag = d1diag();
            moinfo_.new_d1diag = new_d1diag();
            moinfo_.d2diag = d2diag();
            sort_amps();
            update();
            outfile->Printf("\n    Iterations converged.\n");

            outfile->Printf("\n");
            amp_write();
            if (params_.analyze != 0) analyze();
            break;
        }
        if (params_.diis) diis(moinfo_.iter);
        tsave();
        tau_build();
        taut_build();
        last_energy = moinfo_.ecc;
        moinfo_.ecc = energy();
        moinfo_.t1diag = diagnostic();
        moinfo_.d1diag = d1diag();
        moinfo_.new_d1diag = new_d1diag();
        moinfo_.d2diag = d2diag();
        update();
        checkpoint();
    }  // end loop over iterations

    // DGAS Edit
    Process::environment.globals["CC T1 DIAGNOSTIC"] = moinfo_.t1diag;
    Process::environment.globals["CC D1 DIAGNOSTIC"] = moinfo_.d1diag;
    Process::environment.globals["CC NEW D1 DIAGNOSTIC"] = moinfo_.new_d1diag;
    Process::environment.globals["CC D2 DIAGNOSTIC"] = moinfo_.d2diag;

    outfile->Printf("\n");
    if (!done) {
        outfile->Printf("     ** Wave function not converged to %2.1e ** \n", params_.convergence);

        if (params_.aobasis != "NONE" || params_.df) dpd_close(1);
        dpd_close(0);
        cleanup();
        exit_io();
        throw PSIEXCEPTION("Coupled Cluster wave function not converged.");
    }

    outfile->Printf("    SCF energy       (wfn)                    = %20.15f\n", moinfo_.escf);
    outfile->Printf("    Reference energy (file100)                = %20.15f\n", moinfo_.eref);

    // Process::environment.globals["SCF TOTAL ENERGY"] = moinfo.eref;

    if (params_.ref == 0 || params_.ref == 2) {
        if (params_.scs) {
            outfile->Printf("\n    OS SCS-MP2 correlation energy             = %20.15f\n",
                            moinfo_.emp2_os * params_.scsmp2_scale_os);
            outfile->Printf("    SS SCS-MP2 correlation energy             = %20.15f\n",
                            moinfo_.emp2_ss * params_.scsmp2_scale_ss);
            outfile->Printf("    SCS-MP2 correlation energy                = %20.15f\n",
                            moinfo_.emp2_os * params_.scsmp2_scale_os + moinfo_.emp2_ss * params_.scsmp2_scale_ss);
            outfile->Printf(
                "      * SCS-MP2 total energy                  = %20.15f\n",
                moinfo_.eref + moinfo_.emp2_os * params_.scsmp2_scale_os + moinfo_.emp2_ss * params_.scsmp2_scale_ss);

            Process::environment.globals["SCS-MP2 CORRELATION ENERGY"] =
                moinfo_.emp2_os * params_.scsmp2_scale_os + moinfo_.emp2_ss * params_.scsmp2_scale_ss;
            Process::environment.globals["SCS-MP2 TOTAL ENERGY"] =
                moinfo_.eref + moinfo_.emp2_os * params_.scsmp2_scale_os + moinfo_.emp2_ss * params_.scsmp2_scale_ss;
        }
        if (params_.scsn) {
            outfile->Printf("\nOS SCSN-MP2 correlation energy            = %20.15f\n", 0.0);
            outfile->Printf("    SS SCSN-MP2 correlation energy        = %20.15f\n", moinfo_.emp2_ss * 1.76);
            outfile->Printf("    SCSN-MP2 correlation energy           = %20.15f\n", moinfo_.emp2_ss * 1.76);
            outfile->Printf("      * SCSN-MP2 total energy             = %20.15f\n",
                            moinfo_.eref + moinfo_.emp2_ss * 1.76);

            Process::environment.globals["SCS(N)-MP2 CORRELATION ENERGY"] = moinfo_.emp2_ss * 1.76;
            Process::environment.globals["SCS(N)-MP2 TOTAL ENERGY"] = moinfo_.eref + moinfo_.emp2_ss * 1.76;
        }

        outfile->Printf("\n    Opposite-spin MP2 correlation energy      = %20.15f\n", moinfo_.emp2_os);
        outfile->Printf("    Same-spin MP2 correlation energy          = %20.15f\n", moinfo_.emp2_ss);
        outfile->Printf("    Singles MP2 correlation energy            = %20.15f\n", moinfo_.emp2_s);
        outfile->Printf("    MP2 correlation energy                    = %20.15f\n", moinfo_.emp2);
        outfile->Printf("      * MP2 total energy                      = %20.15f\n", moinfo_.eref + moinfo_.emp2);

        Process::environment.globals["MP2 OPPOSITE-SPIN CORRELATION ENERGY"] = moinfo_.emp2_os;
        Process::environment.globals["MP2 SAME-SPIN CORRELATION ENERGY"] = moinfo_.emp2_ss;
        Process::environment.globals["MP2 CORRELATION ENERGY"] = moinfo_.emp2;
        Process::environment.globals["MP2 TOTAL ENERGY"] = moinfo_.eref + moinfo_.emp2;
        set_scalar_variable("MP2 OPPOSITE-SPIN CORRELATION ENERGY", moinfo_.emp2_os);
        set_scalar_variable("MP2 SAME-SPIN CORRELATION ENERGY", moinfo_.emp2_ss);
        set_scalar_variable("MP2 CORRELATION ENERGY", moinfo_.emp2);
        set_scalar_variable("MP2 TOTAL ENERGY", moinfo_.eref + moinfo_.emp2);
    }
    if (params_.wfn == "CC3" || params_.wfn == "EOM_CC3") {
        outfile->Printf("    CC3 correlation energy                    = %20.15f\n", moinfo_.ecc);
        outfile->Printf("      * CC3 total energy                      = %20.15f\n", moinfo_.eref + moinfo_.ecc);

        Process::environment.globals["CC3 CORRELATION ENERGY"] = moinfo_.ecc;
        Process::environment.globals["CC3 TOTAL ENERGY"] = moinfo_.eref + moinfo_.ecc;
        set_scalar_variable("CC3 CORRELATION ENERGY", moinfo_.ecc);
        set_scalar_variable("CC3 TOTAL ENERGY", moinfo_.eref + moinfo_.ecc);
    } else if (params_.wfn == "CC2" || params_.wfn == "EOM_CC2") {
        outfile->Printf("    CC2 correlation energy                    = %20.15f\n", moinfo_.ecc);
        outfile->Printf("      * CC2 total energy                      = %20.15f\n", moinfo_.eref + moinfo_.ecc);

        Process::environment.globals["CC2 CORRELATION ENERGY"] = moinfo_.ecc;
        Process::environment.globals["CC2 TOTAL ENERGY"] = moinfo_.eref + moinfo_.ecc;
        set_scalar_variable("CC2 CORRELATION ENERGY", moinfo_.ecc);
        set_scalar_variable("CC2 TOTAL ENERGY", moinfo_.eref + moinfo_.ecc);

        if (params_.local && local_.weakp == "MP2") {
            outfile->Printf("      * LCC2 (+LMP2) total energy             = %20.15f\n",
                            moinfo_.eref + moinfo_.ecc + local_.weak_pair_energy);
            Process::environment.globals["LCC2 (+LMP2) TOTAL ENERGY"] =
                            moinfo_.eref + moinfo_.ecc + local_.weak_pair_energy;
        }
    } else {
        if (params_.scscc) {
            outfile->Printf("\n    OS SCS-CCSD correlation energy            = %20.15f\n",
                            moinfo_.ecc_os * params_.scscc_scale_os);
            outfile->Printf("    SS SCS-CCSD correlation energy            = %20.15f\n",
                            moinfo_.ecc_ss * params_.scscc_scale_ss);
            outfile->Printf("    SCS-CCSD correlation energy               = %20.15f\n",
                            moinfo_.ecc_os * params_.scscc_scale_os + moinfo_.ecc_ss * params_.scscc_scale_ss);
            outfile->Printf(
                "      * SCS-CCSD total energy                 = %20.15f\n",
                moinfo_.eref + moinfo_.ecc_os * params_.scscc_scale_os + moinfo_.ecc_ss * params_.scscc_scale_ss);

            // LAB TODO  reconsider variable names for ss/os cc
            Process::environment.globals["SCS-CCSD CORRELATION ENERGY"] =
                moinfo_.ecc_os * params_.scscc_scale_os + moinfo_.ecc_ss * params_.scscc_scale_ss;
            Process::environment.globals["SCS-CCSD TOTAL ENERGY"] =
                moinfo_.eref + moinfo_.ecc_os * params_.scscc_scale_os + moinfo_.ecc_ss * params_.scscc_scale_ss;
        }

        outfile->Printf("\n    Opposite-spin CCSD correlation energy     = %20.15f\n", moinfo_.ecc_os);
        outfile->Printf("    Same-spin CCSD correlation energy         = %20.15f\n", moinfo_.ecc_ss);
        outfile->Printf("    Singles CCSD correlation energy           = %20.15f\n", moinfo_.ecc_s);
        outfile->Printf("    CCSD correlation energy                   = %20.15f\n", moinfo_.ecc);
        outfile->Printf("      * CCSD total energy                     = %20.15f\n", moinfo_.eref + moinfo_.ecc);

        Process::environment.globals["CCSD OPPOSITE-SPIN CORRELATION ENERGY"] = moinfo_.ecc_os;
        Process::environment.globals["CCSD SAME-SPIN CORRELATION ENERGY"] = moinfo_.ecc_ss;
        Process::environment.globals["CCSD SINGLES ENERGY"] = moinfo_.ecc_s;
        Process::environment.globals["CCSD DOUBLES ENERGY"] = moinfo_.ecc_os + moinfo_.ecc_ss;
        Process::environment.globals["CCSD CORRELATION ENERGY"] = moinfo_.ecc;
        Process::environment.globals["CCSD TOTAL ENERGY"] = moinfo_.ecc + moinfo_.eref;
        Process::environment.globals["CCSD ITERATIONS"] = moinfo_.iter;
        set_scalar_variable("CCSD OPPOSITE-SPIN CORRELATION ENERGY", moinfo_.ecc_os);
        set_scalar_variable("CCSD SAME-SPIN CORRELATION ENERGY", moinfo_.ecc_ss);
        set_scalar_variable("CCSD SINGLES ENERGY", moinfo_.ecc_s);
        set_scalar_variable("CCSD DOUBLES ENERGY", moinfo_.ecc_os + moinfo_.ecc_ss);
        set_scalar_variable("CCSD CORRELATION ENERGY", moinfo_.ecc);
        set_scalar_variable("CCSD TOTAL ENERGY", moinfo_.ecc + moinfo_.eref);
        set_scalar_variable("CCSD ITERATIONS", moinfo_.iter);

        if (params_.local && local_.weakp == "MP2") {
            outfile->Printf("      * LCCSD (+LMP2) total energy            = %20.15f\n",
                            moinfo_.eref + moinfo_.ecc + local_.weak_pair_energy);
            Process::environment.globals["LCCSD (+LMP2) TOTAL ENERGY"] =
                            moinfo_.eref + moinfo_.ecc + local_.weak_pair_energy;
        }
    }
    outfile->Printf("\n");

    /* Generate the spin-adapted RHF amplitudes for later codes */
    if (params_.ref == 0) spinad_amps();

    /* Compute pair energies */
    std::vector<double> ecc_aa(nocc_act * (nocc_act - 1) /2);
    std::vector<double> ecc_ab(nocc_act * nocc_act);
    if (params_.print_pair_energies) {
        pair_energies(ecc_aa, ecc_ab);
        print_pair_energies(emp2_aa, emp2_ab, ecc_aa, ecc_ab);
    }

    if (params_.wfn == "CC2" && params_.dertype == 1) t1_ijab();

    if ((params_.wfn == "CC3" || params_.wfn == "EOM_CC3") && (params_.dertype == 1 || params_.dertype == 3) &&
        params_.ref == 0) {
        params_.ref = 1;
        /* generate the ROHF versions of the He^T1 intermediates */
        cc3_Wmnij();
        cc3_Wmbij();
        cc3_Wmnie();
        cc3_Wamef();
        cc3_Wabei();
        //    params.ref == 0;
    }

    if (params_.local) {
        /*    local_print_T1_norm(); */
        local_done();
    }

    if (params_.brueckner) Process::environment.globals["BRUECKNER CONVERGED"] = rotate();

    if (params_.aobasis != "NONE" || params_.df) dpd_close(1);
    dpd_close(0);

    if (params_.ref == 2)
        cachedone_uhf(cachelist);
    else
        cachedone_rhf(cachelist);

    cleanup();

    energy_ = moinfo_.ecc + moinfo_.eref;
    name_ = "CCSD";
    Process::environment.globals["CURRENT ENERGY"] = moinfo_.ecc + moinfo_.eref;
    Process::environment.globals["CURRENT CORRELATION ENERGY"] = moinfo_.ecc;
    set_scalar_variable("CURRENT ENERGY", moinfo_.ecc + moinfo_.eref);
    set_scalar_variable("CURRENT CORRELATION ENERGY", moinfo_.ecc);
    // Process::environment.globals["CC TOTAL ENERGY"] = moinfo.ecc+moinfo.eref;
    // Process::environment.globals["CC CORRELATION ENERGY"] = moinfo.ecc;

    exit_io();
    //  if(params.brueckner && brueckner_done)
    //     throw FeatureNotImplemented("CCENERGY", "Brueckner end loop", __FILE__, __LINE__);
    // else

    if ((options_.get_str("WFN") == "CCSD_T")) {
        // Run cctriples
        if (psi::cctriples::cctriples(reference_wavefunction_, options_) == Success) {
            energy_ = Process::environment.globals["CURRENT ENERGY"];
            set_scalar_variable("(T) CORRECTION ENERGY", reference_wavefunction_->scalar_variable("(T) CORRECTION ENERGY"));
            set_scalar_variable("CCSD(T) CORRELATION ENERGY", reference_wavefunction_->scalar_variable("CCSD(T) CORRELATION ENERGY"));
            set_scalar_variable("CCSD(T) TOTAL ENERGY", reference_wavefunction_->scalar_variable("CCSD(T) TOTAL ENERGY"));
            set_scalar_variable("CURRENT CORRELATION ENERGY", reference_wavefunction_->scalar_variable("CCSD(T) CORRELATION ENERGY"));
            set_scalar_variable("CURRENT ENERGY", reference_wavefunction_->scalar_variable("CCSD(T) TOTAL ENERGY"));
        } else {
            energy_ = 0.0;
        }
    }

    return energy_;
}

void CCEnergyWavefunction::init_io() {
    params_.just_energy = 0;
    params_.just_residuals = 0;
    timer_on("ccenergy");
    for (int i = PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i, 1);
}

void CCEnergyWavefunction::title() {
    outfile->Printf("            **************************\n");
    outfile->Printf("            *                        *\n");
    outfile->Printf("            *        CCENERGY        *\n");
    outfile->Printf("            *                        *\n");
    outfile->Printf("            **************************\n");
}

void CCEnergyWavefunction::exit_io() {
    for (int i = PSIF_CC_MIN; i < PSIF_CC_TMP; i++) psio_close(i, 1);
    for (int i = PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio_close(i, 0); /* delete CC_TMP files */
    for (int i = PSIF_CC_TMP11 + 1; i <= PSIF_CC_MAX; i++) psio_close(i, 1);
    timer_off("ccenergy");
}

void CCEnergyWavefunction::checkpoint() {
    for (int i = PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_close(i, 1);
    for (int i = PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i, 1);
}

/* just use T's on disk and don't iterate */
void CCEnergyWavefunction::one_step() {
    dpdfile2 t1;
    dpdbuf4 t2;
    double tval;

    moinfo_.ecc = energy();
    outfile->Printf("\n    Values computed from T amplitudes on disk.\n");
    outfile->Printf("Reference expectation value computed: %20.15lf\n", moinfo_.ecc);
    psio_write_entry(PSIF_CC_HBAR, "Reference expectation value", (char *)&(moinfo_.ecc), sizeof(double));

    if (params_.just_residuals) {
        Fme_build();
        Fae_build();
        Fmi_build();
        t1_build();
        Wmbej_build();
        Z_build();
        Wmnij_build();
        t2_build();
        if ((params_.ref == 0) || (params_.ref == 1)) {
            global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
            global_dpd_->file2_copy(&t1, PSIF_CC_OEI, "FAI residual");
            global_dpd_->file2_close(&t1);
            global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "FAI residual");
            tval = global_dpd_->file2_dot_self(&t1);
            global_dpd_->file2_close(&t1);
            outfile->Printf("    Norm squared of <Phi_I^A|Hbar|0> = %20.15lf\n", tval);
        }
        if (params_.ref == 1) {
            global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "New tia");
            global_dpd_->file2_copy(&t1, PSIF_CC_OEI, "Fai residual");
            global_dpd_->file2_close(&t1);
        } else if (params_.ref == 2) {
            global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 2, 3, "New tia");
            global_dpd_->file2_copy(&t1, PSIF_CC_OEI, "Fai residual");
            global_dpd_->file2_close(&t1);
        }
        if (params_.ref == 0) {
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WAbIj residual");
            global_dpd_->buf4_close(&t2);
            global_dpd_->buf4_init(&t2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WAbIj residual");
            tval = global_dpd_->buf4_dot_self(&t2);
            outfile->Printf("    Norm squared of <Phi^Ij_Ab|Hbar|0>: %20.15lf\n", tval);
            global_dpd_->buf4_close(&t2);
        } else if (params_.ref == 1) {
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WABIJ residual");
            global_dpd_->buf4_close(&t2);
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "Wabij residual");
            global_dpd_->buf4_close(&t2);
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WAbIj residual");
            global_dpd_->buf4_close(&t2);
        } else if (params_.ref == 2) {
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WABIJ residual");
            global_dpd_->buf4_close(&t2);
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "Wabij residual");
            global_dpd_->buf4_close(&t2);
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
            global_dpd_->buf4_copy(&t2, PSIF_CC_HBAR, "WAbIj residual");
            global_dpd_->buf4_close(&t2);
        }
    }
}

}  // namespace ccenergy
}  // namespace psi
