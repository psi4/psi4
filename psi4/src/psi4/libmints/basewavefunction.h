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

#include <array>
#include <memory>

namespace psi {

class Molecule;
class BasisSet;
class Options;

/*! \ingroup MINTS
 *  \class BaseWavefunction
 *  \brief Common polymorphic root for real and complex wavefunctions.
 */
class PSI_API BaseWavefunction {
   protected:
    /// Molecule that this wavefunction is run on
    std::shared_ptr<Molecule> molecule_;

    /// The ORBITAL basis
    std::shared_ptr<BasisSet> basisset_;

    /// Options object
    Options& options_;

    /// How big of a field perturbation to apply
    std::array<double, 3> dipole_field_strength_;

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
