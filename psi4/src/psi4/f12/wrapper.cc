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

#include "mp2.h"

namespace psi {
namespace f12 {

SharedWavefunction f12(SharedWavefunction ref_wfn, Options& options) {
    std::shared_ptr<Wavefunction> f12;
    if (options.get_str("F12_TYPE").find("DISK") != std::string::npos) {
        f12 = std::make_shared<DiskMP2F12>(ref_wfn, options);
    } else {
        f12 = std::make_shared<MP2F12>(ref_wfn, options);
    }

    return f12;
}

}  // namespace f12
}  // namespace psi
