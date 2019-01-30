/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#include "psi4/cc/ccwave.h"
#include "psi4/libpsio/psio.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsi4util/process.h"
#include <memory>
namespace psi {
namespace ccenergy {
std::map<std::string, SharedMatrix> CCEnergyWavefunction::get_amplitudes() {
    // re-init dpd
    int** cachelist;
    auto cachefiles = std::vector<int>(PSIO_MAXUNIT);
    auto ref = options_.get_str("REFERENCE");
    if (ref == "ROHF" && options_.get_bool("SEMICANONICAL")){
      ref = "UHF";
    }
    std::vector<DPDMOSpace> spaces;
    std::map<std::string, SharedMatrix> amps;
    if (dpd_list[0] == nullptr) {
        if (ref == "RHF") {
            Dimension occpi_ = nalphapi_;
            Dimension virtpi_ = nmopi_ - nalphapi_;
            cachelist = cacheprep_rhf(options_.get_int("CACHELEVEL"), cachefiles.data());
            spaces = {DPDMOSpace{'o', "ijkl", occpi_}, DPDMOSpace{'v', "abcd", virtpi_}};
        } else if (ref == "ROHF") {
            Dimension doccpi_ = nalphapi_;
            Dimension soccpi_ = nbetapi_ - nalphapi_;
            Dimension occpi_ = doccpi_ + soccpi_;
            Dimension virtpi_ = (nmopi_ - doccpi_) + soccpi_;
            cachelist = cacheprep_rhf(options_.get_int("CACHELEVEL"), cachefiles.data());
            spaces = {DPDMOSpace{'o', "ijkl", occpi_}, DPDMOSpace{'v', "abcd", virtpi_}};
        } else /*UHF*/{
            Dimension aoccpi_ = nalphapi_;
            Dimension boccpi_ = nbetapi_;
            Dimension avirtpi_ = nmopi_ - nalphapi_;
            Dimension bvirtpi_ = nmopi_ - nbetapi_;
            cachelist = cacheprep_uhf(options_.get_int("CACHELEVEL"), cachefiles.data());
            spaces = {DPDMOSpace('O', "IJKL", aoccpi_), DPDMOSpace('V', "ABCD", avirtpi_),
                      DPDMOSpace('o', "ijkl", boccpi_), DPDMOSpace('v', "abcd", bvirtpi_)};
        }
        dpd_list[0] = new DPD(0, nirrep_, Process::environment.get_memory(), 0, cachefiles.data(), cachelist, nullptr,
                              (int)spaces.size(), spaces);
    }
    dpd_set_default(0);

    bool oei_open = psio_open_check(PSIF_CC_OEI);
    bool tamps_open = psio_open_check(PSIF_CC_TAMPS);
    if (!oei_open) {
        psio_open(PSIF_CC_OEI, PSIO_OPEN_OLD);
    }
    if (!tamps_open) {
        psio_open(PSIF_CC_TAMPS, PSIO_OPEN_OLD);
    }

    // grab T1
    dpdfile2 T1;
    dpdbuf4 T2;
    if (ref == "RHF") {
        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
        amps["tIA"] = std::make_shared<Matrix>(&T1);
        global_dpd_->file2_close(&T1);

        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        amps["tIjAb"] = std::make_shared<Matrix>(&T2);
        global_dpd_->buf4_close(&T2);
    } else if (ref == "ROHF") {
        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
        amps["tIA"] = std::make_shared<Matrix>(&T1);
        global_dpd_->file2_close(&T1);
        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
        amps["tia"] = std::make_shared<Matrix>(&T1);
        global_dpd_->file2_close(&T1);
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        amps["tIjAb"] = std::make_shared<Matrix>(&T2);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
        amps["tijab"] = std::make_shared<Matrix>(&T2);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
        amps["tIJAB"] = std::make_shared<Matrix>(&T2);
        global_dpd_->buf4_close(&T2);
    } else {  // if ref == UHF
        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
        amps["tIA"] = std::make_shared<Matrix>(&T1);
        global_dpd_->file2_close(&T1);
        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
        amps["tia"] = std::make_shared<Matrix>(&T1);
        global_dpd_->file2_close(&T1);
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
        amps["tIjAb"] = std::make_shared<Matrix>(&T2);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
        amps["tijab"] = std::make_shared<Matrix>(&T2);
        global_dpd_->buf4_close(&T2);
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
        amps["tIJAB"] = std::make_shared<Matrix>(&T2);
        global_dpd_->buf4_close(&T2);
    }
    dpd_close(0);
    if (ref == "RHF") {
        cachedone_rhf(cachelist);
    } else {
        cachedone_uhf(cachelist);
    }
    if (!oei_open) {
        psio_close(PSIF_CC_OEI, 1);
    }
    if (!tamps_open) {
        psio_close(PSIF_CC_TAMPS, 1);
    }
    return amps;
}
}  // namespace ccenergy
}  // namespace psi
