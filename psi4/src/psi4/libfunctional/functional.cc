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

#include "functional.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"

namespace psi {

Functional::Functional() { common_init(); }
Functional::~Functional() {}
void Functional::common_init() {
    name_ = "";
    description_ = "";
    citation_ = "";
    xclib_description_ = "";
    alpha_ = 1.0;
    omega_ = 0.0;
}
std::shared_ptr<Functional> Functional::build_worker() { throw PSIEXCEPTION("Functional: pseudo-abstract class."); }
void Functional::print(std::string out, int level) const {
    if (level < 1) return;
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    printer->Printf("   => %s Functional <=\n\n", name_.c_str());

    printer->Printf("%s", description_.c_str());
    printer->Printf("\n");

    printer->Printf("%s", citation_.c_str());
    printer->Printf("\n");

    printer->Printf("    GGA   = %14s\n", (this->is_gga() ? "TRUE" : "FALSE"));
    printer->Printf("    Meta  = %14s\n", (this->is_meta() ? "TRUE" : "FALSE"));
    printer->Printf("    LRC   = %14s\n", (this->is_lrc() ? "TRUE" : "FALSE"));
    printer->Printf("    Alpha = %14.6E\n", alpha_);
    printer->Printf("    Omega = %14.6E\n", omega_);
    printer->Printf("\n");
}
void Functional::compute_functional(const std::map<std::string, SharedVector>& in,
                                    const std::map<std::string, SharedVector>& out, int npoints, int deriv) {
    throw PSIEXCEPTION("Functional: pseudo-abstract class.");
}
void Functional::set_density_cutoff(double cut) { throw PSIEXCEPTION("Functional: pseudo-abstract class."); }
}  // namespace psi
