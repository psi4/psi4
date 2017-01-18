/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmints_moindexspace_h_
#define _psi_src_lib_libmints_moindexspace_h_

#include <string>
#include "typedefs.h"
#include "psi4/libmints/dimension.h"

namespace psi {
class BasisSet;
class IntegralFactory;

/**
 * @brief The OrbitalSpace class
 *
 * Describes an orbital space. Contains a unique identifier that can be used to
 * describe the space. Also contains shared_ptr references to AO->MO transformation
 * matrix and possible orbital energies. Relavent basis set and integral factory
 * are also contained.
 */
class OrbitalSpace
{
    /// Unique identifier
    std::string id_;
    /// Name of the orbital space.
    std::string name_;
    /// AO->MO transformation matrix (ao x mo) or SO->MO transformation matrix
    SharedMatrix C_;

    /// MO "eigenvalues"
    std::shared_ptr<Vector> evals_;

    /// AO basis set
    std::shared_ptr<BasisSet> basis_;
    /// Integral factory that as
    std::shared_ptr<IntegralFactory> ints_;

    /// MO Dimensionality
    Dimension dim_; // dim_.n() better equal nirrep_

    /// No default constructor
    OrbitalSpace();

public:
    OrbitalSpace(const std::string& id,
                 const std::string& name,
                 const SharedMatrix& full_C,
                 const std::shared_ptr<Vector>& evals,
                 const std::shared_ptr<BasisSet>& basis,
                 const std::shared_ptr<IntegralFactory>& ints);

    OrbitalSpace(const std::string& id,
                 const std::string& name,
                 const SharedMatrix& full_C,
                 const std::shared_ptr<BasisSet>& basis,
                 const std::shared_ptr<IntegralFactory>& ints);

    OrbitalSpace(const std::string& id,
                 const std::string& name,
                 const std::shared_ptr<Wavefunction>& wave);

    int nirrep() const;
    const std::string& id() const;
    const std::string& name() const;

    /// C - transformation matrix (AO x MO)
    const SharedMatrix& C() const;

    /// "Eigenvalues" of the C matrix
    const std::shared_ptr<Vector>& evals() const;

    /// The AO basis set used to create C
    const std::shared_ptr<BasisSet>& basisset() const;

    /// Integral factory used to create C
    const std::shared_ptr<IntegralFactory>& integral() const;

    /// MO dimensionality
    const Dimension& dim() const;

    /// Print information about the orbital space
    void print() const;

    /** Creates an OrbitalSpace from 'from' to the given basis set 'to'
      */
    static OrbitalSpace transform(const OrbitalSpace& from, const std::shared_ptr<BasisSet>& to);

    /** Returns the overlap matrix between space1 and space2.
        The matrix has dimensions of space2.C().coldim() and
        space1.C().coldim().
        Throws if the overlap cannot be computed.
      */
    static SharedMatrix overlap(const OrbitalSpace& space1, const OrbitalSpace& space2);
    /** Returns the overlap matrix between basis1 and basis2.
        Throws if the overlap cannot be computed.
      */
    static SharedMatrix overlap(const std::shared_ptr<BasisSet>& basis1,
                                const std::shared_ptr<BasisSet>& basis2);

    /** Given two spaces, it projects out one space from the other and returns the new spaces.
     * \param orb_space The space to project out. The returned space will be orthogonal to this.
     * \param ri_space The space being projected on. The returned space will be this space minus orb_space.
     * \param linear_tol The tolerance for linear dependencies.
     */
    static OrbitalSpace build_cabs_space(
            const OrbitalSpace& orb_space,
            const OrbitalSpace& ri_space,
            double linear_tol);

    /** Given two basis sets, it merges the basis sets and then constructs an orthogonalized
     * space with the same span. Linearly dependent orbitals are thrown out.
     * \param molecule molecule to construct the basis for
     * \param obs_key option keyword for orbital basis set "BASIS"
     * \param aux_key option keyword for auxiliery basis set "DF_BASIS_MP2"
     * \param lindep_tol The tolerance for linear dependencies
     */
    static OrbitalSpace build_ri_space(const std::shared_ptr<Molecule>& molecule, const std::string& obs_key, const std::string& aux_key, double lindep_tol);

    /** Given a basis set, it orthogonalizes the orbitals and returns a space with the same
     * span but orthogonal orbitals. Also, linear dependent orbitals are projected out.
     * \param aux_bs The basis to orthogonalize
     * \param lindep_tol The tolerance for linear dependencies
     */
    static OrbitalSpace build_abs_space(std::shared_ptr<BasisSet> aux_bs, std::shared_ptr<IntegralFactory> ints, double lindep_tol);

};

namespace SpaceBuilder
{
}

}

#endif // _psi_src_lib_libmints_moindexspace_h_
