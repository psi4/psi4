/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

namespace libint2 {
class Engine;
class Shell;
}  // namespace libint2

namespace psi {

class BasisSet;
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

    // These are scratch arrays that are used for the integral engines that do not yet use Libint2
    // under the hood.  When everything is converted to use Libint2, they can be deleted.
    double* buffer_;
    double* target_;
    double* tformbuf_;

    size_t count_;
    int deriv_;
    int natom_;

    int nchunk_;

    int buffer_size_;

    /// Libint2 engines
    std::unique_ptr<libint2::Engine> engine0_;
    std::unique_ptr<libint2::Engine> engine1_;  // derivatives
    std::unique_ptr<libint2::Engine> engine2_;  // hessians

    /// Pointers to each chunk of (derivative/multipole) integrals.  For simple integral types
    /// (e.g. overlap, kinetic) there is only one chunk of integrals.  For dipole integrals
    /// there are three such buffers in memory at a given time, pointed to by these pointers.
    /// The ordering of the dipoles is x, y, z; generally CCA library ordering is used to determine
    /// the layout of operators and derivatives.  The Libint2 wiki has a detailed description of this.
    std::vector<const double*> buffers_;

    std::vector<std::pair<int, int>> shellpairs_;

    OneBodyAOInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                 int deriv = 0);
    void set_chunks(int nchunk) { nchunk_ = nchunk; }
    void pure_transform(const libint2::Shell &s1, const libint2::Shell &s2, int nchunks = 1);

    /// Compute integrals for a given shell pair
    virtual void compute_pair(const libint2::Shell&, const libint2::Shell&);
    /// Compute first derivative integrals for a given shell pair
    virtual void compute_pair_deriv1(const libint2::Shell&, const libint2::Shell&);
    /// Compute second derivative integrals for a given shell pair
    virtual void compute_pair_deriv2(const libint2::Shell&, const libint2::Shell&);
    /// Whether the operator is antisymmetric with respect to interchange of the bra and ket
    virtual bool is_antisymmetric() const { return false; }


   public:
    virtual ~OneBodyAOInt();

    /// Basis set on center one.
    std::shared_ptr<BasisSet> basis();
    /// Basis set on center one.
    std::shared_ptr<BasisSet> basis1();
    /// Basis set on center two.
    std::shared_ptr<BasisSet> basis2();

    const auto& shellpairs() const { return shellpairs_; }

    /// Number of chunks. Normally 1, but dipoles (3) quadrupoles (6).
    int nchunk() const { return nchunk_; }

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

    /// Compute the integrals between basis function in the given shell pair.
    virtual void compute_shell(int, int);
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
    virtual void set_origin(const Vector3& _origin) { origin_ = _origin; }

    /// Buffer where each chunk of integrals is placed
    const std::vector<const double*>& buffers() const { return buffers_; }
};

typedef std::shared_ptr<OneBodyAOInt> SharedOneBodyAOInt;

/// For a pair of basis sets, provides a list of integer pairs that index all significant (i.e. whose overlap is
/// above the threshold) pairs that exist.  If the basis sets are different, the full Cartesian product is considered
/// but if they are the same, only the lower triangular significant pairs are returned.
std::vector<std::pair<int, int>> build_shell_pair_list_no_spdata(std::shared_ptr<BasisSet> bs1,
                                                                 std::shared_ptr<BasisSet> bs2, double threshold);

}  // namespace psi

#endif
