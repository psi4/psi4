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

// Private implementation helpers shared between the C and C++ IWL APIs.
// Not part of the installed public surface; do not include from outside
// psi4/src/psi4/libiwl/.

#ifndef _psi_src_lib_libiwl_iwl_impl_h_
#define _psi_src_lib_libiwl_iwl_impl_h_

#include <climits>
#include <memory>
#include <string>
#include "config.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"

namespace psi {
namespace iwl_impl {

// Verify that all four orbital indices fit in the on-disk Label width. The
// bucket format stores each label as int16, so any larger value would be
// silently truncated and corrupt the file. Throws PSIEXCEPTION instead.
inline void check_label_fits(int p, int q, int r, int s) {
    if (p > SHRT_MAX || q > SHRT_MAX || r > SHRT_MAX || s > SHRT_MAX || p < SHRT_MIN || q < SHRT_MIN || r < SHRT_MIN ||
        s < SHRT_MIN) {
        throw PSIEXCEPTION(
            "libiwl: orbital index does not fit in 16-bit Label (max " + std::to_string(SHRT_MAX) +
            "). The on-disk IWL bucket format cannot store this case. p=" + std::to_string(p) +
            " q=" + std::to_string(q) + " r=" + std::to_string(r) + " s=" + std::to_string(s));
    }
}

// Construct a PsiOutStream once per call, lazily, only when printing is
// actually requested. The pre-refactor code constructed one per integral.
inline std::shared_ptr<PsiOutStream> make_printer(const std::string &out) {
    return (out == "outfile") ? outfile : std::make_shared<PsiOutStream>(out);
}

}  // namespace iwl_impl
}  // namespace psi

#endif
