/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#ifndef _psi_src_lib_libmints_giao_potential_h_
#define _psi_src_lib_libmints_giao_potential_h_

#include <vector>
#include "psi4/libmints/onebody.h" // needed because we derive from OneBodyAOInt
#include "psi4/libmints/osrecur.h"

namespace psi {

    class BasisSet;
    class GaussianShell;
    class SphericalTransform;

/*! \ingroup MINTS
 *  \class GiaoPotentialInt
 *  \brief This class computes integrals of the type 1/2 <mu|R_MN x r V|nu>, where R_MN = (R_M - R_N), 
 *   M and N are Gaussian centers and V is the potential energy operator.
 *
 *   Use an IntegralFactory to create this object.
 */
class GiaoPotentialInt : public OneBodyAOInt
{
    /// Computes integrals between two shell objects.
    void compute_pair(const GaussianShell& , const GaussianShell&);
protected:
    //! Generic Obara Saika recursion object.
    ObaraSaikaTwoCenterVIRecursion* potential_recur_;
    /// Matrix of coordinates/charges of partial charges
    SharedMatrix Zxyz_;

public:
    //! Constructor, it assumes you are computing the first derivatives by default
    //! Do not call directly use an IntegralFactory.
    GiaoPotentialInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv=0);
    //! Virtual destructor
    virtual ~GiaoPotentialInt();

    /// Set the field of charges
    void set_charge_field(SharedMatrix Zxyz) { Zxyz_ = Zxyz; }

    /// Get the field of charges
    SharedMatrix charge_field() const { return Zxyz_; }

    /// Does the method provide first derivatives?
    bool has_deriv1() { return false; }
    /// Does the method provide second derivatives?
    bool has_deriv2() { return false; }
};

}

#endif
