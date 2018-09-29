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

#include "ccwave.h"

#include <vector>
#include <map>

#include "psi4/psi4-dec.h"

#include "psi4/libdpd/dpd.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"

#include "ccmoinfo.h"
#include "ccparams.h"

namespace psi {
namespace cc {

void psio_on() {
    for (int i = PSIF_CC_OEI; i <= PSIF_CC_MAX; i++) psio_open(i, 1);
}

void psio_off() {
    for (int i = PSIF_CC_OEI; i < PSIF_CC_TMP; i++) psio_close(i, 1);
    for (int i = PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio_close(i, 0); /* delete CC_TMP files */
    for (int i = PSIF_CC_TMP11 + 1; i <= PSIF_CC_MAX; i++) psio_close(i, 1);
}

CCWavefunction::CCWavefunction(std::shared_ptr<Wavefunction> reference_wavefunction) : Wavefunction(Process::environment.options) {
    timer_on("ccwavefunction");
    timer_on("initialization");
    // Copy the wavefuntion then update
    shallow_copy(reference_wavefunction);
    set_reference_wavefunction(reference_wavefunction);
    common_init();
    timer_off("initialization");
}

CCWavefunction::CCWavefunction(std::shared_ptr<Wavefunction> reference_wavefunction, Options &options) : Wavefunction(options) {
    timer_on("ccwavefunction");
    // Copy the wavefuntion then update
    timer_on("initialization");
    shallow_copy(reference_wavefunction);
    set_reference_wavefunction(reference_wavefunction);
    common_init();
    timer_off("initialization");
}

CCWavefunction::~CCWavefunction() {
    // Close coupled cluster files
    psio_off();

    timer_off("ccwavefunction");
}

double CCWavefunction::compute_energy() { return 0.0; }

void CCWavefunction::title(std::string &wfn) {
    outfile->Printf("\n");
    outfile->Printf("         ---------------------------------------------------------\n");
    outfile->Printf("                          Coupled Cluster\n");
    outfile->Printf("                           %s wavefunction\n", wfn.c_str());
    outfile->Printf("\n");
    outfile->Printf("                 T. Daniel Crawford\n");
    outfile->Printf("         ---------------------------------------------------------\n");
    outfile->Printf("\n");
}

void CCWavefunction::common_init() {
    // Open coupled cluster files
    psio_on();

    // Calculation parameters
    params_ = CCParams(options_);
    title(params_.wfn);
    moinfo_ = CCMOInfo(reference_wavefunction_, params_.ref);

    // Print out information
    print_parameters(params_, memory_);
    print_ccmoinfo(moinfo_);
}

void CCWavefunction::init_dpd() {
    // cachefiles_.reserve(PSIO_MAXUNIT);
    // std::vector<int *> spaces;
    // std::vector<int *> aospaces;
    // if (params_.ref == 2) {
    //    cachelist_ = cacheprep_uhf(params_.cachelev, cachefiles_.data());
    //    spaces.push_back(moinfo_.aoccpi);
    //    spaces.push_back(moinfo_.aocc_sym);
    //    spaces.push_back(moinfo_.avirtpi);
    //    spaces.push_back(moinfo_.avir_sym);
    //    spaces.push_back(moinfo_.boccpi);
    //    spaces.push_back(moinfo_.bocc_sym);
    //    spaces.push_back(moinfo_.bvirtpi);
    //    spaces.push_back(moinfo_.bvir_sym);
    //    if (params_.aobasis != "NONE") {
    //        aospaces.push_back(moinfo_.aoccpi);
    //        aospaces.push_back(moinfo_.aocc_sym);
    //        aospaces.push_back(moinfo_.sopi);
    //        aospaces.push_back(moinfo_.sosym);
    //        aospaces.push_back(moinfo_.boccpi);
    //        aospaces.push_back(moinfo_.bocc_sym);
    //        aospaces.push_back(moinfo_.sopi);
    //        aospaces.push_back(moinfo_.sosym);
    //    }
    //} else {
    //    cachelist_ = cacheprep_rhf(params_.cachelev, cachefiles_.data());
    //    spaces.push_back(moinfo_.occpi);
    //    spaces.push_back(moinfo_.occ_sym);
    //    spaces.push_back(moinfo_.virtpi);
    //    spaces.push_back(moinfo_.vir_sym);
    //    if (params_.aobasis != "NONE") {
    //        aospaces.push_back(moinfo_.occpi);
    //        aospaces.push_back(moinfo_.occ_sym);
    //        aospaces.push_back(moinfo_.sopi);
    //        aospaces.push_back(moinfo_.sosym);
    //    }
    //}

    // dpd_["mo"].init(0, moinfo_.nirreps, params_.memory, params_.cachetype, cachefiles_.data(), cachelist_,
    //                cache_priority_list_.data(), spaces.size() / 2, spaces);

    // if (aospaces.size()) {
    //    dpd_["ao"].init(1, moinfo_.nirreps, params_.memory, 0, cachefiles_.data(), cachelist_, nullptr,
    //                    aospaces.size() / 2, aospaces);
    //}
}

void CCWavefunction::tear_down() {
    // Free up cache
    //    for (auto &&i : dpd_) {
    //        i.second.file2_cache_close();
    //        i.second.file4_cache_close();
    //    }
    //
    //    if (params_.ref == 2) {
    //        cachedone_uhf(cachelist_);
    //    } else {
    //        cachedone_rhf(cachelist_);
    //    }
}

}  // namespace cc
}  // namespace psi
