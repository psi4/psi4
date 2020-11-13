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

#ifndef _psi_src_lib_libmints_giao_angmom_h_
#define _psi_src_lib_libmints_giao_angmom_h_

#include <vector>
#include "psi4/libmints/onebody.h" // needed because we derive from OneBodyAOInt
#include "psi4/libmints/osrecur.h"

namespace psi {

    class BasisSet;
    class GaussianShell;
    class SphericalTransform;

/*! \ingroup MINTS
 *  \class GiaoAngmomInt
 *  \brief This class computes integrals of the type <mu|L_N|nu>,
 *   where L_N = (R - R_N) x p, and N is the center of |nu>. 
 *   
 *
 *   Use an IntegralFactory to create this object.
 */
class GiaoAngmomInt : public OneBodyAOInt
{
    //! Generic Obara Saika recursion object.
    ObaraSaikaTwoCenterRecursion overlap_recur_;

    void compute_pair(const GaussianShell& , const GaussianShell&);

public:
    //! Constructor, it assumes you are computing the first derivatives by default
    //! Do not call directly use an IntegralFactory.
    GiaoAngmomInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv=0);
    //! Virtual destructor
    virtual ~GiaoAngmomInt();

    /// Does the method provide first derivatives?
    bool has_deriv1() { return false; }
    /// Does the method provide second derivatives?
    bool has_deriv2() { return false; }
};

}

#endif
