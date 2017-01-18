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

#ifndef _psi_src_lib_libmints_twobody_h
#define _psi_src_lib_libmints_twobody_h

 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif
#ifdef _XOPEN_SOURCE
#undef _XOPEN_SOURCE
#endif
#include "psi4/libpsi4util/exception.h"

namespace psi {

enum PermutedOrder { ABCD = 0, BACD = 1, ABDC = 2, BADC = 3, CDAB = 4, CDBA = 5, DCAB = 6, DCBA = 7 };


class IntegralFactory;
class AOShellCombinationsIterator;
class BasisSet;
class GaussianShell;

/*! \ingroup MINTS
 *  \class TwoBodyInt
 *  \brief Two body integral base class.
 */
class TwoBodyAOInt
{
protected:
    const IntegralFactory* integral_;

    std::shared_ptr<BasisSet> bs1_;
    std::shared_ptr<BasisSet> bs2_;
    std::shared_ptr<BasisSet> bs3_;
    std::shared_ptr<BasisSet> bs4_;

    const std::shared_ptr<BasisSet> original_bs1_;
    const std::shared_ptr<BasisSet> original_bs2_;
    const std::shared_ptr<BasisSet> original_bs3_;
    const std::shared_ptr<BasisSet> original_bs4_;

    /// Buffer to hold the final integrals.
    double *target_;
    /// Number of integrals in the current buffer
    int curr_buff_size_;
    /// Buffer to hold the transformation intermediates.
    double *tformbuf_;
    /// Buffer to hold the initially computed integrals.
    double *source_;
    /// Maximum number of unique quartets needed to compute a set of SO's
    int max_unique_quartets_;
    /// Number of atoms.
    int natom_;
    /// Derivative level.
    int deriv_;
    /// Whether to force integrals to be generated in the Cartesian (AO) basis;
    bool force_cartesian_;
    /// How the shells were reordered for libint
    PermutedOrder permuted_order_;

    void permute_target(double *s, double *t, int sh1, int sh2, int sh3, int sh4, bool p12, bool p34, bool p13p24);
    void permute_1234_to_1243(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_2134(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_2143(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_3412(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_4312(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_3421(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_4321(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);

//    TwoBodyInt(std::shared_ptr<BasisSet> bs1,
//               std::shared_ptr<BasisSet> bs2,
//               std::shared_ptr<BasisSet> bs3,
//               std::shared_ptr<BasisSet> bs4,
//               int deriv = 0);

    TwoBodyAOInt(const IntegralFactory* intsfactory, int deriv=0);

public:
    virtual ~TwoBodyAOInt();

    /// Basis set on center one
    std::shared_ptr<BasisSet> basis();
    /// Basis set on center one
    std::shared_ptr<BasisSet> basis1();
    /// Basis set on center two
    std::shared_ptr<BasisSet> basis2();
    /// Basis set on center three
    std::shared_ptr<BasisSet> basis3();
    /// Basis set on center four
    std::shared_ptr<BasisSet> basis4();

    /// Sets whether we're forcing this object to always generate Cartesian integrals
    void set_force_cartesian(bool t_f) { force_cartesian_ = t_f; }

    /// Returns the derivative level this object is setup for.
    int deriv() const { return deriv_; }

    /// Buffer where the integrals are placed
    const double *buffer() const { return target_; }

    /// Returns the integral factory used to create this object
    const IntegralFactory* integral() const { return integral_; }

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    virtual size_t compute_shell(const AOShellCombinationsIterator&) = 0;

    /// Compute the integrals
    virtual size_t compute_shell(int, int, int, int) = 0;

    /// Is the shell zero?
    virtual int shell_is_zero(int,int,int,int) { return 0; }

    /// Compute the first derivatives
    virtual size_t compute_shell_deriv1(int, int, int, int) = 0;

    /// Compute the second derivatives
    virtual size_t compute_shell_deriv2(int, int, int, int) = 0;

    /// Normalize Cartesian functions based on angular momentum
    void normalize_am(std::shared_ptr<GaussianShell>, std::shared_ptr<GaussianShell>, std::shared_ptr<GaussianShell>, std::shared_ptr<GaussianShell>, int nchunk=1);

    /// Return true if the clone member can be called. By default returns false.
    virtual bool cloneable();

    /// Returns a clone of this object. By default throws an exception
    virtual TwoBodyAOInt* clone();

    /// Results go back to buffer_
    void pure_transform(int, int, int, int, int nchunk);
};

typedef std::shared_ptr<TwoBodyAOInt> SharedTwoBodyAOInt;

}

#endif
