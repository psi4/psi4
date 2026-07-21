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
};

}  // namespace psi

#endif
