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

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/sobasis.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/process.h"

#include <sstream>

namespace psi {

void BaseWavefunction::common_init() {
    // Options-only / derived-class constructors may not have a basis yet (cf. Wavefunction(Options)).
    if (!basisset_) {
        return;
    }

    // Check the point group of the molecule. If it is not set, set it.
    if (!molecule_->point_group()) {
        molecule_->set_point_group(molecule_->find_point_group());
    }

    // Create an SO basis...we need the point group for this part.
    integral_ = std::make_shared<IntegralFactory>(basisset_, basisset_, basisset_, basisset_);
    mintshelper_ = std::make_shared<MintsHelper>(basisset_, options_);
    sobasisset_ = std::make_shared<SOBasisSet>(basisset_, integral_);

    // Obtain the dimension object to initialize the factory.
    nsopi_ = sobasisset_->dimension();
    nsopi_.set_name("SOs per irrep");

    nirrep_ = nsopi_.n();

    // Initialize arrays that hold dimensionality information
    nmopi_ = Dimension(nirrep_, "MOs per irrep");
    frzcpi_ = Dimension(nirrep_, "Frozen core orbitals per irrep");
    frzvpi_ = Dimension(nirrep_, "Frozen virtual orbitals per irrep");

    // Obtain memory amount from the environment
    memory_ = Process::environment.get_memory();

    nso_ = basisset_->nbf();
    nmo_ = basisset_->nbf();
    for (int k = 0; k < nirrep_; k++) {
        nmopi_[k] = 0;
    }

    energy_ = 0.0;
    efzc_ = 0.0;

    // Read in the debug flag
    debug_ = options_.get_int("DEBUG");
    print_ = options_.get_int("PRINT");

    // Determine the number of electrons in the system
    nelectron_ = 0;
    for (int i = 0; i < molecule_->natom(); ++i) {
        nelectron_ += (int)molecule_->Z(i);
    }
    nelectron_ -= molecule_->molecular_charge();

    // Make sure that the multiplicity is reasonable
    multiplicity_ = molecule_->multiplicity();
    if (multiplicity_ - 1 > nelectron_) {
        std::ostringstream oss;
        oss << "There are not enough electrons for multiplicity = " << multiplicity_ << ".\n";
        oss << "Please check your input";
        throw SanityCheckError(oss.str(), __FILE__, __LINE__);
    }
    if (multiplicity_ % 2 == nelectron_ % 2) {
        std::ostringstream oss;
        oss << "A multiplicity of " << multiplicity_ << " with " << nelectron_ << " electrons is impossible.\n";
        oss << "Please check your input";
        throw SanityCheckError(oss.str(), __FILE__, __LINE__);
    }

    // Setup dipole field perturbation information
    perturb_h_ = options_.get_bool("PERTURB_H");
    dipole_field_type_ = nothing;
    std::fill(dipole_field_strength_.begin(), dipole_field_strength_.end(), 0.0);
    if (perturb_h_) {
        std::string perturb_with;
        if (options_["PERTURB_WITH"].has_changed()) {
            perturb_with = options_.get_str("PERTURB_WITH");
            // Do checks to see what perturb_with is.
            if (perturb_with == "DIPOLE_X") {
                dipole_field_type_ = dipole_x;
                dipole_field_strength_[0] = options_.get_double("PERTURB_MAGNITUDE");
                outfile->Printf(
                    " WARNING: the DIPOLE_X and PERTURB_MAGNITUDE keywords are deprecated."
                    "  Use DIPOLE and the PERTURB_DIPOLE array instead.");
            } else if (perturb_with == "DIPOLE_Y") {
                dipole_field_type_ = dipole_y;
                dipole_field_strength_[1] = options_.get_double("PERTURB_MAGNITUDE");
                outfile->Printf(
                    " WARNING: the DIPOLE_Y and PERTURB_MAGNITUDE keywords are deprecated."
                    "  Use DIPOLE and the PERTURB_DIPOLE array instead.");
            } else if (perturb_with == "DIPOLE_Z") {
                dipole_field_type_ = dipole_z;
                dipole_field_strength_[2] = options_.get_double("PERTURB_MAGNITUDE");
                outfile->Printf(
                    " WARNING: the DIPOLE_Z and PERTURB_MAGNITUDE keywords are deprecated."
                    "  Use DIPOLE and the PERTURB_DIPOLE array instead.");
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
                outfile->Printf(
                    "  WARNING: Option PERTURB_WITH=SPHERE is deprecated and may be removed as soon as v1.13.\n");
                dipole_field_type_ = sphere;
            } else {
                outfile->Printf("Unknown PERTURB_WITH. Applying no perturbation.\n");
            }
        } else {
            outfile->Printf("PERTURB_H is true, but PERTURB_WITH not found, applying no perturbation.\n");
        }
    }
}

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

std::shared_ptr<BasisSet> BaseWavefunction::get_basisset(std::string label) const {
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
