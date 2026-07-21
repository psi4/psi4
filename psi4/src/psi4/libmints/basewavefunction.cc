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

#include "basewavefunction.h"
#include "wavefunction.h"

#include "psi4/libmints/mintshelper.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/process.h"

namespace psi {

void BaseWavefunction::common_init() {}

BaseWavefunction::BaseWavefunction()
    : options_(Process::environment.options), dipole_field_strength_{{0.0, 0.0, 0.0}}, PCM_enabled_(false) {}

BaseWavefunction::BaseWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis,
                                   Options& options)
    : molecule_(molecule),
      basisset_(basis),
      options_(options),
      dipole_field_strength_{{0.0, 0.0, 0.0}},
      PCM_enabled_(false) {}

BaseWavefunction::BaseWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis)
    : molecule_(molecule),
      basisset_(basis),
      options_(Process::environment.options),
      dipole_field_strength_{{0.0, 0.0, 0.0}},
      PCM_enabled_(false) {}

BaseWavefunction::BaseWavefunction(Options& options)
    : options_(options), dipole_field_strength_{{0.0, 0.0, 0.0}}, PCM_enabled_(false) {}

BaseWavefunction::~BaseWavefunction() = default;

std::map<std::string, std::shared_ptr<BasisSet>> BaseWavefunction::basissets() const {
    return mintshelper_->basissets();
}

std::shared_ptr<BasisSet> BaseWavefunction::get_basisset(std::string label) {
    return mintshelper_->get_basisset(label);
}

void BaseWavefunction::set_basisset(std::string label, std::shared_ptr<BasisSet> basis) {
    return mintshelper_->set_basisset(label, basis);
}

bool BaseWavefunction::basisset_exists(std::string label) { return mintshelper_->basisset_exists(label); }

void BaseWavefunction::set_frzvpi(const Dimension& frzvpi) {
    for (int h = 0; h < nirrep_; h++) {
        frzvpi_[h] = frzvpi[h];
    }
}

void BaseWavefunction::set_energy(double ene) { set_scalar_variable("CURRENT ENERGY", ene); }

bool BaseWavefunction::has_scalar_variable(const std::string& key) { return variables_.count(to_upper_copy(key)); }

bool BaseWavefunction::has_potential_variable(const std::string& key) {
    return potentials_.count(to_upper_copy(key));
}

double BaseWavefunction::scalar_variable(const std::string& key) {
    std::string uc_key = to_upper_copy(key);

    auto search = variables_.find(uc_key);
    if (search != variables_.end()) {
        return search->second;
    } else {
        throw PSIEXCEPTION("BaseWavefunction::scalar_variable: Requested variable " + uc_key + " was not set!\n");
    }
}

std::shared_ptr<ExternalPotential> BaseWavefunction::potential_variable(const std::string& key) {
    std::string uc_key = to_upper_copy(key);

    auto search = potentials_.find(uc_key);
    if (search != potentials_.end()) {
        return search->second;
    } else {
        throw PSIEXCEPTION("BaseWavefunction::potential_variable: Requested variable " + uc_key + " was not set!\n");
    }
}

void BaseWavefunction::set_scalar_variable(const std::string& key, double val) {
    variables_[to_upper_copy(key)] = val;

    if (to_upper_copy(key) == "CURRENT ENERGY") energy_ = val;
}

void BaseWavefunction::set_potential_variable(const std::string& key, std::shared_ptr<ExternalPotential> val) {
    potentials_[to_upper_copy(key)] = val;
}

int BaseWavefunction::del_scalar_variable(const std::string& key) { return variables_.erase(to_upper_copy(key)); }

int BaseWavefunction::del_potential_variable(const std::string& key) { return potentials_.erase(to_upper_copy(key)); }

std::map<std::string, double> BaseWavefunction::scalar_variables() { return variables_; }

std::map<std::string, std::shared_ptr<ExternalPotential>> BaseWavefunction::potential_variables() {
    return potentials_;
}

void BaseWavefunction::set_PCM(const std::shared_ptr<PCM>& pcm) {
    PCM_ = pcm;
    PCM_enabled_ = true;
}

}  // namespace psi
