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

#ifndef _psi_src_lib_libmints_3coverlap_h
#define _psi_src_lib_libmints_3coverlap_h

#include "psi4/libmints/osrecur.h"
#include "psi4/libmints/integral.h"
#include "psi4/libpsi4util/exception.h"


namespace psi {

class IntegralFactory;
class BasisSet;
class GaussianShell;

/** \ingroup MINTS
    \class ThreeCenterOverlapInt
    \brief Three center overlap integral.
 */
class ThreeCenterOverlapInt
{
protected:
    ObaraSaikaThreeCenterRecursion overlap_recur_;

    std::shared_ptr<BasisSet> bs1_;
    std::shared_ptr<BasisSet> bs2_;
    std::shared_ptr<BasisSet> bs3_;

    /// Buffer to hold the source integrals.
    double *buffer_;

    /// Buffer for spherical harmonics
    double *temp_;

    /// Vector of Sphericaltransforms
    std::vector<SphericalTransform> st_;

    void compute_pair(const GaussianShell& s1,
                      const GaussianShell& s2,
                      const GaussianShell& s3);

public:
    ThreeCenterOverlapInt(std::vector<SphericalTransform> st,
                          std::shared_ptr<BasisSet> bs1,
                          std::shared_ptr<BasisSet> bs2,
                          std::shared_ptr<BasisSet> bs3);

    ThreeCenterOverlapInt(std::shared_ptr<BasisSet> bs1,
                          std::shared_ptr<BasisSet> bs2,
                          std::shared_ptr<BasisSet> bs3);

    virtual ~ThreeCenterOverlapInt();

    /// Basis set on center one.
    std::shared_ptr<BasisSet> basis();
    /// Basis set on center one.
    std::shared_ptr<BasisSet> basis1();
    /// Basis set on center two.
    std::shared_ptr<BasisSet> basis2();
    /// Basis set on center three.
    std::shared_ptr<BasisSet> basis3();

    /// Buffer where the integrals are placed.
    const double *buffer() const { return buffer_; }

    /// Compute the integrals of the form (a|c|b).
    virtual void compute_shell(int, int, int);

    /// Normalize Cartesian functions based on angular momentum
    void normalize_am(const GaussianShell&,
                      const GaussianShell&,
                      const GaussianShell&);

    /// Perform pure (spherical) transform.
    void pure_transform(const GaussianShell&,
                        const GaussianShell&,
                        const GaussianShell&);
};

}

#endif
