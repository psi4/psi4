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

#include "psi4/pragma.h"
#include <memory>
#include "functional.h"
#include <xc.h>
#include "psi4/psi4-dec.h"
#include "LibXCfunctional.h"
#include "psi4/libpsi4util/exception.h"

namespace psi {

std::shared_ptr<Functional> Functional::build_base(const std::string& alias) {
    Functional* fun;

    if (xc_functional_get_number(alias.c_str()) >= 0) {
        fun = static_cast<Functional*>(new LibXCFunctional(alias, false));
    } else {
        throw PSIEXCEPTION("Functional::build_base: Unrecognized base Functional.");
    }

    return std::shared_ptr<Functional>(fun);
}
}  // namespace psi
