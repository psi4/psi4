/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#ifndef _psi_src_lib_libmints_onebody_h_
#define _psi_src_lib_libmints_onebody_h_

#include <vector>
#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif
#ifdef _XOPEN_SOURCE
#undef _XOPEN_SOURCE
#endif
#include "psi4/libpsi4util/exception.h"
#include "typedefs.h"
#include "psi4/libmints/vector3.h"

#include "psi4/libpsi4util/exception.h"

namespace psi {

class BasisSet;
class GaussianShell;
class SphericalTransform;

/*! \ingroup MINTS
 *  \class OneBodyInt
 *  \brief Basis class for all one-electron integrals.
 */
class PSI_API OneBodyAOInt {
   protected:
    std::shared_ptr<BasisSet> bs1_;
    std::shared_ptr<BasisSet> bs2_;
    std::vector<SphericalTransform>& spherical_transforms_;

    Vector3 origin_;

    double* buffer_;
    double* target_;
    double* tformbuf_;

    size_t count_;
    int deriv_;
    int natom_;

    int nchunk_;

    int buffer_size_;

    OneBodyAOInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                 int deriv = 0);

    virtual void compute_pair(const GaussianShell& s1, const GaussianShell& s2) = 0;
    virtual void compute_pair_deriv1(const GaussianShell& s1, const GaussianShell& s2);
    virtual void compute_pair_deriv2(const GaussianShell& s1, const GaussianShell& s2);

    void set_chunks(int nchunk) { nchunk_ = nchunk; }
    void pure_transform(const GaussianShell&, const GaussianShell&, int = 1);

    /// Normalize Cartesian functions based on angular momentum
    void normalize_am(const GaussianShell&, const GaussianShell&, int nchunk = 1);

   public:
    virtual ~OneBodyAOInt();

    /// Basis set on center one.
    std::shared_ptr<BasisSet> basis();
    /// Basis set on center one.
    std::shared_ptr<BasisSet> basis1();
    /// Basis set on center two.
    std::shared_ptr<BasisSet> basis2();

    /// Number of chunks. Normally 1, but dipoles (3) quadrupoles (6).
    int nchunk() const { return nchunk_; }

    /// Buffer where the integrals are placed.
    const double* buffer() const;

    /// Compute the integrals between basis function in the given shell pair.
    void compute_shell(int, int);

    /*! @{
     * Computes all integrals and stores them in result
     * @param result Shared matrix object that will hold the results.
     */
    void compute(SharedMatrix& result);
    /*! @} */

    /// Computes all integrals and stores them in result by default this method throws
    virtual void compute(std::vector<SharedMatrix>& result);

    /// Does the method provide first derivatives?
    virtual bool has_deriv1() { return false; }

    /// Does the method provide second derivatives?
    virtual bool has_deriv2() { return false; }

    /// What order of derivative was requested?
    int deriv() const { return deriv_; }

    /// Computes the first derivatives and stores them in result
    virtual void compute_deriv1(std::vector<SharedMatrix>& result);

    /// Computes the second derivatives and stores them in result
    virtual void compute_deriv2(std::vector<SharedMatrix>& result);

    /// Computes the integrals between basis function in the given shell pair
    virtual void compute_shell_deriv1(int, int);
    /// Computes the integrals between basis function in the given shell pair
    virtual void compute_shell_deriv2(int, int);

    /// Return true if the clone member can be called. By default returns false.
    virtual bool cloneable() const;

    /// Returns a clone of this object. By default throws an exception.
    virtual OneBodyAOInt* clone() const;

    /// Returns the origin (useful for properties)
    Vector3 origin() const { return origin_; }

    /// Set the origin (useful for properties)
    void set_origin(const Vector3& _origin) { origin_ = _origin; }
};

typedef std::shared_ptr<OneBodyAOInt> SharedOneBodyAOInt;
}  // namespace psi

#endif
