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

#ifndef _psi_src_lib_libmints_basewavefunction_h
#define _psi_src_lib_libmints_basewavefunction_h

#include "psi4/psi4-dec.h"
#include "psi4/libmints/dimension.h"

#include <array>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace psi {

class Molecule;
class BasisSet;
class IntegralFactory;
class MintsHelper;
class Wavefunction;
class MatrixFactory;
class Options;
class SOBasisSet;
class PCM;
class PSIO;
class ExternalPotential;

/*! \ingroup MINTS
 *  \class BaseWavefunction
 *  \brief Common polymorphic root for real and complex wavefunctions.
 */
class PSI_API BaseWavefunction {
   protected:
    /// Name of the wavefunction
    std::string name_;

    /// Module name for CURRENT ENERGY
    std::string module_;

    /// The ORBITAL basis
    std::shared_ptr<BasisSet> basisset_;

    /// Primary basis set for SO integrals
    std::shared_ptr<SOBasisSet> sobasisset_;

    /// Molecule that this wavefunction is run on
    std::shared_ptr<Molecule> molecule_;

    /// Options object
    Options& options_;

    // PSI file access variables
    std::shared_ptr<PSIO> psio_;

    /// Integral factory
    std::shared_ptr<IntegralFactory> integral_;

    /// MintsHelper
    std::shared_ptr<MintsHelper> mintshelper_;

    /// Matrix factory for creating standard sized matrices
    std::shared_ptr<MatrixFactory> factory_;

    std::shared_ptr<Wavefunction> reference_wavefunction_;

    /// How much memory you have access to.
    long int memory_;

    /// Perturb the Hamiltonian?
    int perturb_h_;
    /// With what...
    enum FieldType { nothing, dipole_x, dipole_y, dipole_z, dipole, embpot, dx, sphere };
    FieldType dipole_field_type_;

    /// How big of a field perturbation to apply
    std::array<double, 3> dipole_field_strength_;

    /// Debug flag
    size_t debug_;
    /// Print flag
    size_t print_;

    /// Total frozen core orbitals
    int nfrzc_;

    /// Number of frozen core per irrep
    Dimension frzcpi_;
    /// Number of frozen virtuals per irrep
    Dimension frzvpi_;

    /// Number of so per irrep
    Dimension nsopi_;
    /// Number of mo per irrep
    Dimension nmopi_;

    /// The energy associated with this wavefunction
    double energy_;

    /// Frozen-core energy associated with this wavefunction
    double efzc_;

    /// Total number of SOs
    int nso_;
    /// Total number of MOs
    int nmo_;
    /// Number of irreps
    int nirrep_;

    /// Should nuclear electrostatic potentials be available, they will be here
    std::shared_ptr<std::vector<double>> esp_at_nuclei_;

    /// If atomic point charges are available they will be here
    std::shared_ptr<std::vector<double>> atomic_point_charges_;

    /// Should natural orbital occupations be available, they will be here
    std::vector<std::vector<std::tuple<double, int, int>>> no_occupations_;

    // The external potential for the current wave function
    std::shared_ptr<ExternalPotential> external_pot_;

    // Collection of scalar variables
    std::map<std::string, double> variables_;

    // Collection of external potentials; this member variable is provisional and might be removed in the future.
    // This member variable is currently used for passing ExternalPotential objects to the F/I-SAPT code
    // The above defined external_pot_ member variable contains the total external potential defined for the current
    // wave function. For F/I-SAPT, we need a set of external potential that can be assigned to either the interacting
    // fragments or to the environment
    // For F/I-SAPT the keys can be A, B, or C (all optionals), where A and B signify the interacting subsystem
    // and C signify the envirnoment
    std::map<std::string, std::shared_ptr<ExternalPotential>> potentials_;

    // Polarizable continuum model
    std::shared_ptr<PCM> PCM_;
    /// Polarizable continuum model enabled?
    bool PCM_enabled_;

    /// Subclass-specific initialization after shared members are set.
    /// Must be called from the derived constructor body (not the base ctor)
    /// so virtual dispatch reaches the derived override.
    virtual void common_init();

   public:
    BaseWavefunction();

    /// Constructor for an entirely new wavefunction with an existing basis
    BaseWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis, Options& options);

    /// Constructor for an entirely new wavefunction with an existing basis and global options
    BaseWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis);

    /// Blank constructor for derived classes
    BaseWavefunction(Options& options);

    virtual ~BaseWavefunction();

    /// Returns the molecule object that pertains to this wavefunction.
    std::shared_ptr<Molecule> molecule() const { return molecule_; }
    /// Returns the basis set object that pertains to this wavefunction.
    std::shared_ptr<BasisSet> basisset() const { return basisset_; }

    /// Returns the PSIO object that pertains to this wavefunction.
    std::shared_ptr<PSIO> psio() const { return psio_; }
    Options& options() const { return options_; }

    /// An integral factory with basisset() on each center.
    std::shared_ptr<IntegralFactory> integral() const { return integral_; }
    /// An molecular integrals helper with basisset() on each center.
    std::shared_ptr<MintsHelper> mintshelper() const { return mintshelper_; }
    /// Returns the SO basis set object that pertains to this wavefunction.
    std::shared_ptr<SOBasisSet> sobasisset() const { return sobasisset_; }

    /// Getters and setters for other basis sets
    std::map<std::string, std::shared_ptr<BasisSet>> basissets() const;
    std::shared_ptr<BasisSet> get_basisset(std::string label);
    void set_basisset(std::string label, std::shared_ptr<BasisSet> basis);
    bool basisset_exists(std::string label);

    /// Returns the MatrixFactory object that pertains to this wavefunction
    std::shared_ptr<MatrixFactory> matrix_factory() const { return factory_; }
    /// Returns the reference wavefunction
    std::shared_ptr<Wavefunction> reference_wavefunction() const { return reference_wavefunction_; }
    /// Sets the reference wavefunction
    void set_reference_wavefunction(const std::shared_ptr<Wavefunction> wfn) { reference_wavefunction_ = wfn; }

    /// Returns the print level
    int get_print() const { return print_; }
    // Set the print flag level
    void set_print(size_t print) { print_ = print; }
    // Set the debug flag level
    void set_debug(size_t debug) { debug_ = debug; }

    /// Returns the number of SOs per irrep array.
    const Dimension& nsopi() const { return nsopi_; }
    /// Returns the number of MOs per irrep array.
    const Dimension& nmopi() const { return nmopi_; }
    /// Returns the frozen core orbitals per irrep array.
    const Dimension& frzcpi() const { return frzcpi_; }
    /// Returns the frozen virtual orbitals per irrep array.
    const Dimension& frzvpi() const { return frzvpi_; }
    /// Sets the frozen virtual orbitals per irrep array.
    void set_frzvpi(const Dimension& frzvpi);

    /* Return the magnitude of the dipole perturbation strength in the x,y,z direction */
    std::array<double, 3> get_dipole_field_strength() const { return dipole_field_strength_; }
    FieldType get_dipole_perturbation_type() const { return dipole_field_type_; }

    /// Return the number of frozen core orbitals
    int nfrzc() const { return nfrzc_; }
    /// Returns the number of SOs
    int nso() const { return nso_; }
    /// Returns the number of MOs
    int nmo() const { return nmo_; }
    /// Returns the number of irreps
    int nirrep() const { return nirrep_; }
    double energy() const { return energy_; }
    /// Sets the energy
    void set_energy(double ene);
    /// Returns the frozen-core energy
    double efzc() const { return efzc_; }
    /// Sets the frozen-core energy
    void set_efzc(double efzc) { efzc_ = efzc; }

    /// Returns electrostatic potentials at nuclei
    std::shared_ptr<std::vector<double>> esp_at_nuclei() const { return esp_at_nuclei_; }
    /// Sets the electrostatic potentials at nuclei
    void set_esp_at_nuclei(const std::shared_ptr<std::vector<double>>& nesps) { esp_at_nuclei_ = nesps; }

    /// Returns the atomic point charges
    std::shared_ptr<std::vector<double>> atomic_point_charges() const { return atomic_point_charges_; }
    /// Sets the atomic point charges
    void set_atomic_point_charges(const std::shared_ptr<std::vector<double>>& apcs) { atomic_point_charges_ = apcs; }

    /// Returns NO occupations
    std::vector<std::vector<std::tuple<double, int, int>>> no_occupations() const { return no_occupations_; }
    /// Sets the NO occupations
    void set_no_occupations(const std::vector<std::vector<std::tuple<double, int, int>>> no_ocs) {
        no_occupations_ = no_ocs;
    }

    /// Set the wavefunction name (e.g. "RHF", "ROHF", "UHF", "CCEnergyWavefunction")
    void set_name(const std::string& name) { name_ = name; }
    /// Returns the wavefunction name
    const std::string& name() const { return name_; }

    /// Set the module name (e.g. "OCC", "CCENERGY", "CCT3")
    void set_module(const std::string& module) { module_ = module; }
    /// Returns the module name
    const std::string& module() const { return module_; }

    // Get the external potential
    std::shared_ptr<ExternalPotential> external_pot() const { return external_pot_; }
    // Set the external potential
    void set_external_potential(std::shared_ptr<ExternalPotential> external) { external_pot_ = external; }

    /// Get and set scalar and potential variable dictionaries
    bool has_scalar_variable(const std::string& key);
    // The function below is provisional and might be removed in the future
    bool has_potential_variable(const std::string& key);
    double scalar_variable(const std::string& key);
    // The function below is provisional and might be removed in the future
    std::shared_ptr<ExternalPotential> potential_variable(const std::string& key);
    void set_scalar_variable(const std::string& key, double value);
    // The function below is provisional and might be removed in the future
    void set_potential_variable(const std::string& key, std::shared_ptr<ExternalPotential> value);
    int del_scalar_variable(const std::string& key);
    // The function below is provisional and might be removed in the future
    int del_potential_variable(const std::string& key);
    std::map<std::string, double> scalar_variables();
    // The function below is provisional and might be removed in the future
    std::map<std::string, std::shared_ptr<ExternalPotential>> potential_variables();

    /// Set PCM object
    void set_PCM(const std::shared_ptr<PCM>& pcm);
    /// Get PCM object
    std::shared_ptr<PCM> get_PCM() const { return PCM_; }
    bool PCM_enabled() const { return PCM_enabled_; }
};

}  // namespace psi

#endif
