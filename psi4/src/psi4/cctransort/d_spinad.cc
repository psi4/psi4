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

#include "psi4/libdpd/dpd.h"

namespace psi {
namespace cctransort {

void d_spinad() {
    dpdbuf4 D;

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_scmcopy(&D, PSIF_CC_DINTS, "D 2<ij|ab> - <ij|ba>", 2);
    global_dpd_->buf4_sort_axpy(&D, PSIF_CC_DINTS, pqsr, 0, 5, "D 2<ij|ab> - <ij|ba>", -1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, prqs, 10, 10, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->buf4_close(&D);
}

}  // namespace cctransort
}  // namespace psi
