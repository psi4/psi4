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

#ifndef _psi_src_lib_libmints_electricfield_h_
#define _psi_src_lib_libmints_electricfield_h_

#include <vector>
#include "typedefs.h"

namespace psi {

/*! \ingroup MINTS
 *  \class ElectricFieldInt
 *  \brief Computes electric field integrals.
 *
 *  Use an IntegralFactory to create this object.
 */
class ElectricFieldInt : public OneBodyAOInt
{
    //! Obara and Saika recursion object to be used.
    ObaraSaikaTwoCenterVIDeriv2Recursion efield_recur_;

    //! Number of atoms.
    int natom_;

    //! Computes the electric field between two gaussian shells.
    void compute_pair(const GaussianShell&, const GaussianShell&);

    //! Computes the electric field gradient between two gaussian shells.
    void compute_pair_deriv1(const GaussianShell&, const GaussianShell&);

public:
    //! Constructor. Do not call directly use an IntegralFactory.
    ElectricFieldInt(std::vector<SphericalTransform>&, boost::shared_ptr<BasisSet>, boost::shared_ptr<BasisSet>, int deriv=0);
    //! Virtual destructor
    virtual ~ElectricFieldInt();

    //! Does the method provide first derivatives?
    bool has_deriv1() { return true; }

    static Vector3 nuclear_contribution(const Vector3 &origin, boost::shared_ptr<Molecule> mol);
    static SharedMatrix nuclear_contribution_to_gradient(const Vector3 &origin, boost::shared_ptr<Molecule> mol);
};

}

#endif

