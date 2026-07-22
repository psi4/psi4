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

#include <sstream>

namespace psi {

void ComplexWavefunction::common_init() {
    // Options-only / derived-class constructors may not have a basis yet (cf. Wavefunction(Options)).
    if (!basisset_) {
        return;
    }

    if (!molecule_->point_group()) {
        molecule_->set_point_group(molecule_->find_point_group());
    }

    integral_ = std::make_shared<IntegralFactory>(basisset_, basisset_, basisset_, basisset_);
    mintshelper_ = std::make_shared<MintsHelper>(basisset_, options_);
    sobasisset_ = std::make_shared<SOBasisSet>(basisset_, integral_);

    nsopi_ = sobasisset_->dimension();
    nsopi_.set_name("SOs per irrep");

    nirrep_ = nsopi_.n();

    // Complex overlap S_ is allocated/filled once ComplexMatrix construction from SO integrals is wired.

    nmopi_ = Dimension(nirrep_, "MOs per irrep");
    nelecpi_ = Dimension(nirrep_, "Electrons per irrep");
    frzcpi_ = Dimension(nirrep_, "Frozen core orbitals per irrep");
    frzvpi_ = Dimension(nirrep_, "Frozen virtual orbitals per irrep");

    memory_ = Process::environment.get_memory();

    nso_ = basisset_->nbf();
    nmo_ = basisset_->nbf();
    for (int k = 0; k < nirrep_; k++) {
        nmopi_[k] = 0;
        nelecpi_[k] = 0;
    }

    energy_ = 0.0;
    efzc_ = 0.0;

    debug_ = options_.get_int("DEBUG");
    print_ = options_.get_int("PRINT");

    int nelectron = 0;
    for (int i = 0; i < molecule_->natom(); ++i) {
        nelectron += (int)molecule_->Z(i);
    }
    nelectron -= molecule_->molecular_charge();

    int multiplicity = molecule_->multiplicity();
    if (multiplicity - 1 > nelectron) {
        std::ostringstream oss;
        oss << "There are not enough electrons for multiplicity = " << multiplicity << ".\n";
        oss << "Please check your input";
        throw SanityCheckError(oss.str(), __FILE__, __LINE__);
    }
    if (multiplicity % 2 == nelectron % 2) {
        std::ostringstream oss;
        oss << "A multiplicity of " << multiplicity << " with " << nelectron << " electrons is impossible.\n";
        oss << "Please check your input";
        throw SanityCheckError(oss.str(), __FILE__, __LINE__);
    }

    nelec_ = nelectron;

    perturb_h_ = options_.get_bool("PERTURB_H");
    dipole_field_type_ = nothing;
    std::fill(dipole_field_strength_.begin(), dipole_field_strength_.end(), 0.0);
    if (perturb_h_) {
        std::string perturb_with;
        if (options_["PERTURB_WITH"].has_changed()) {
            perturb_with = options_.get_str("PERTURB_WITH");
            if (perturb_with == "DIPOLE_X") {
                dipole_field_type_ = dipole_x;
                dipole_field_strength_[0] = options_.get_double("PERTURB_MAGNITUDE");
            } else if (perturb_with == "DIPOLE_Y") {
                dipole_field_type_ = dipole_y;
                dipole_field_strength_[1] = options_.get_double("PERTURB_MAGNITUDE");
            } else if (perturb_with == "DIPOLE_Z") {
                dipole_field_type_ = dipole_z;
                dipole_field_strength_[2] = options_.get_double("PERTURB_MAGNITUDE");
            } else if (perturb_with == "DIPOLE") {
                dipole_field_type_ = dipole;
                if (options_["PERTURB_DIPOLE"].size() != 3)
                    throw PSIEXCEPTION("The PERTURB dipole should have exactly three floating point numbers.");
                for (int n = 0; n < 3; ++n) dipole_field_strength_[n] = options_["PERTURB_DIPOLE"][n].to_double();
            } else if (perturb_with == "EMBPOT") {
                dipole_field_type_ = embpot;
            } else if (perturb_with == "DX") {
                dipole_field_type_ = dx;
            } else if (perturb_with == "SPHERE") {
                dipole_field_type_ = sphere;
            }
        }
    }
}

ComplexWavefunction::ComplexWavefunction() : BaseWavefunction() { common_init(); }

ComplexWavefunction::ComplexWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis,
                                         Options& options)
    : BaseWavefunction(molecule, basis, options) {
    common_init();
}

ComplexWavefunction::ComplexWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis)
    : BaseWavefunction(molecule, basis) {
    common_init();
}

ComplexWavefunction::ComplexWavefunction(Options& options) : BaseWavefunction(options) { common_init(); }

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
