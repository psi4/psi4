/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#include "fcidump_helper.h"

#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <string>

#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"

#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"

namespace psi {
namespace fcidump {
void fcidump_tei_helper(int nirrep, bool restricted, std::map<std::string, int> DPD_info, double ints_tolerance,
                        std::string fname) {
    outfile->Printf("Writing TEI integrals in FCIDUMP format to " + fname + "\n");
    // Append to the file created by the fcidump function Python-side
    auto mode = std::ostream::app;
    auto intdump = std::make_shared<PsiOutStream>(fname.c_str(), mode);

    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(DPD_info["instance_id"]);

    _default_psio_lib_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    dpdbuf4 K;

    // RHF
    if (restricted) {
        /* Convert a molecular orbital index [0,1,...] to [1,2,...] (i.e. from zero-based to one-based). */
        auto mo_index = [](const int i) { return i + 1; };
        // We want only the permutationally unique integrals, see libtrans documentation for details
        // DPD_info["alpha_MO"] is DPD_ID("[A>=A]+")
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, DPD_info["alpha_MO"], DPD_info["alpha_MO"],
                               DPD_info["alpha_MO"], DPD_info["alpha_MO"], 0, "MO Ints (AA|AA)");
        detail::write_tei_to_disk(intdump, nirrep, K, ints_tolerance, mo_index, mo_index);
        global_dpd_->buf4_close(&K);
    } else {
        /* Convert an alpha spin-orbital index [0,1,...] to [1,3,...] (i.e. from
         * zero-based to one-based, with corresponding beta orbitals interwoven).
         */
        auto alpha_index = [](const int i) { return 2 * i + 1; };
        /* Convert a beta spin-orbital index [0,1,...] to [2,4,...] (i.e. from
         * zero-based to one-based, with corresponding alpha orbitals interwoven).
         */
        auto beta_index = [](const int i) { return 2 * (i + 1); };
        // We want only the permutationally unique integrals, see libtrans documentation for details
        // DPD_info["alpha_MO"] is DPD_ID("[A>=A]+"), while DPD_info["beta_MO"] is DPD_ID("[a>=a]+")
        // alpha-alpha
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, DPD_info["alpha_MO"], DPD_info["alpha_MO"],
                               DPD_info["alpha_MO"], DPD_info["alpha_MO"], 0, "MO Ints (AA|AA)");
        detail::write_tei_to_disk(intdump, nirrep, K, ints_tolerance, alpha_index, alpha_index);
        global_dpd_->buf4_close(&K);
        // beta-beta
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, DPD_info["beta_MO"], DPD_info["beta_MO"], DPD_info["beta_MO"],
                               DPD_info["beta_MO"], 0, "MO Ints (aa|aa)");
        detail::write_tei_to_disk(intdump, nirrep, K, ints_tolerance, beta_index, beta_index);
        global_dpd_->buf4_close(&K);
        // alpha-beta
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, DPD_info["alpha_MO"], DPD_info["beta_MO"],
                               DPD_info["alpha_MO"], DPD_info["beta_MO"], 0, "MO Ints (AA|aa)");
        detail::write_tei_to_disk(intdump, nirrep, K, ints_tolerance, alpha_index, beta_index);
        global_dpd_->buf4_close(&K);
    }
    _default_psio_lib_->close(PSIF_LIBTRANS_DPD, 1);
}

namespace detail {
void write_tei_to_disk(std::shared_ptr<PsiOutStream> intdump, int nirrep, dpdbuf4& K, double ints_tolerance,
                       OrbitalIndexing indx1, OrbitalIndexing indx2) {
    for (int h = 0; h < nirrep; ++h) {
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for (int pq = 0; pq < K.params->rowtot[h]; ++pq) {
            int p = K.params->roworb[h][pq][0];
            int q = K.params->roworb[h][pq][1];
            for (int rs = 0; rs < K.params->coltot[h]; ++rs) {
                int r = K.params->colorb[h][rs][0];
                int s = K.params->colorb[h][rs][1];

                if (std::abs(K.matrix[h][pq][rs]) > ints_tolerance)
                    intdump->Printf("%28.20E%4d%4d%4d%4d\n", K.matrix[h][pq][rs], indx1(p), indx1(q), indx2(r),
                                    indx2(s));
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
}
}  // End namespace detail
}  // End namespace fcidump
}  // End namespace psi
