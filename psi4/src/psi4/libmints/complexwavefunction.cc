/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "psi4/libmints/complexwavefunction.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/process.h"

namespace psi {

void ComplexWavefunction::common_init() {
    BaseWavefunction::common_init();

    // Initialize einsums and turn off logging to stdout.
    const char* ein_argv[4] = {"psi4", "--einsums:no-profiler-report", "--einsums:log-level", "3"};
    einsums::initialize(4, ein_argv);

    nelecpi_ = Dimension(nirrep_, "Electrons per irrep");
    nelec_ = nelectron_;
}

ComplexWavefunction::ComplexWavefunction() : BaseWavefunction() {}

ComplexWavefunction::ComplexWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis,
                                         Options& options)
    : BaseWavefunction(molecule, basis, options) {
    common_init();
}

ComplexWavefunction::ComplexWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis)
    : BaseWavefunction(molecule, basis) {
    common_init();
}

ComplexWavefunction::ComplexWavefunction(Options& options) : BaseWavefunction(options) {}

ComplexWavefunction::~ComplexWavefunction() {}

SharedComplexMatrix ComplexWavefunction::C() const {
    if (!C_) {
        throw PSIEXCEPTION("ComplexWavefunction::C: Unable to obtain MO coefficients.");
    }
    return C_;
}

void ComplexWavefunction::set_gradient(SharedComplexMatrix grad) { set_array_variable("CURRENT GRADIENT", grad); }

void ComplexWavefunction::set_hessian(SharedComplexMatrix hess) { set_array_variable("CURRENT HESSIAN", hess); }

bool ComplexWavefunction::has_array_variable(const std::string& key) { return arrays_.count(to_upper_copy(key)); }

SharedComplexMatrix ComplexWavefunction::array_variable(const std::string& key) {
    std::string uc_key = to_upper_copy(key);

    auto search = arrays_.find(uc_key);
    if (search != arrays_.end()) {
        return std::make_shared<ComplexMatrix>(*search->second);
    } else {
        throw PSIEXCEPTION("ComplexWavefunction::array_variable: Requested variable " + uc_key + " was not set!\n");
    }
}

void ComplexWavefunction::set_array_variable(const std::string& key, SharedComplexMatrix val) {
    arrays_[to_upper_copy(key)] = std::make_shared<ComplexMatrix>(*val);

    if (to_upper_copy(key) == "CURRENT GRADIENT") gradient_ = std::make_shared<ComplexMatrix>(*val);
    if (to_upper_copy(key) == "CURRENT HESSIAN") hessian_ = std::make_shared<ComplexMatrix>(*val);
}

int ComplexWavefunction::del_array_variable(const std::string& key) { return arrays_.erase(to_upper_copy(key)); }

std::map<std::string, SharedComplexMatrix> ComplexWavefunction::array_variables() { return arrays_; }

}  // namespace psi
