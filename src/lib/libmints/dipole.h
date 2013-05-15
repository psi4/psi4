/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef _psi_src_lib_libmints_dipole_h_
#define _psi_src_lib_libmints_dipole_h_

#include <vector>
#include "typedefs.h"

namespace psi {

/*! \ingroup MINTS
 *  \class DipoleInt
 *  \brief Computes dipole integrals.
 *
 * Use an IntegralFactory to create this object. */
class DipoleInt : public OneBodyAOInt
{
    //! Obara and Saika recursion object to be used.
    ObaraSaikaTwoCenterRecursion overlap_recur_;

    //! Computes the dipole between two gaussian shells.
    void compute_pair(const GaussianShell&, const GaussianShell&);
    //! Computes the dipole derivative between two gaussian shells.
    void compute_pair_deriv1(const GaussianShell&, const GaussianShell&);

public:
    //! Constructor. Do not call directly use an IntegralFactory.
    DipoleInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>, int deriv=0);
    //! Virtual destructor
    virtual ~DipoleInt();

    //! Does the method provide first derivatives?
    bool has_deriv1() { return true; }

    /// Returns the nuclear contribution to the dipole moment
    static SharedVector nuclear_contribution(boost::shared_ptr<Molecule> mol, const Vector3 &origin);
};

}

#endif
