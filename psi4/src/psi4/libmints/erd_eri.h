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

#ifndef psi4_libmints_erd_eri_h_
#define psi4_libmints_erd_eri_h_

#ifdef USING_erd

#include "psi4/libmints/twobody.h"

typedef int F_INT;
typedef int F_BOOL;

namespace psi{


class IntegralFactory;
class AOShellCombinationsIterator;

class ERDTwoElectronInt : public TwoBodyAOInt
{


protected:
    /// The list of renormalized contraction coefficients for center 1
    double *new_cc_1_;
    /// The list of renormalized contraction coefficients for center 2
    double *new_cc_2_;
    /// The list of renormalized contraction coefficients for center 3
    double *new_cc_3_;
    /// The list of renormalized contraction coefficients for center 4
    double *new_cc_4_;
    /// All basis set 1 exponents, stored as a flat array
    double *alpha_1_;
    /// All basis set 2 exponents, stored as a flat array
    double *alpha_2_;
    /// All basis set 3 exponents, stored as a flat array
    double *alpha_3_;
    /// All basis set 4 exponents, stored as a flat array
    double *alpha_4_;
    /// The x,y,z coordinates for each shell in basis set 1
    double *xyz_1_;
    /// The x,y,z coordinates for each shell in basis set 2
    double *xyz_2_;
    /// The x,y,z coordinates for each shell in basis set 3
    double *xyz_3_;
    /// The x,y,z coordinates for each shell in basis set 4
    double *xyz_4_;
    /// The list of contraction coefficients
    double *cc_;
    /// The list of exponents
    double *alpha_;

    /// The current size of the integral buffer
    size_t d_buffer_size_;
    /// The current size of the integer scratch space
    size_t i_buffer_size_;
    /// The address of the first contraction coefficient for each shell on center 1
    int *pgto_offsets_1_;
    /// The address of the first contraction coefficient for each shell on center 2
    int *pgto_offsets_2_;
    /// The address of the first contraction coefficient for each shell on center 3
    int *pgto_offsets_3_;
    /// The address of the first contraction coefficient for each shell on center 4
    int *pgto_offsets_4_;
    /// The number of primitive GTOs per shell in basis set 1
    int *npgto_1_;
    /// The number of primitive GTOs per shell in basis set 2
    int *npgto_2_;
    /// The number of primitive GTOs per shell in basis set 3
    int *npgto_3_;
    /// The number of primitive GTOs per shell in basis set 4
    int *npgto_4_;
    /// The angular momentum of each shell in basis set 1
    int *am_1_;
    /// The angular momentum of each shell in basis set 2
    int *am_2_;
    /// The angular momentum of each shell in basis set 3
    int *am_3_;
    /// The angular momentum of each shell in basis set 4
    int *am_4_;
    /// The integer scratch space
    F_INT *iscratch_;
    /// The double scratch space, which has junk at the start, and integrals at the end
    double *dscratch_;
    /// The start address in the target integral buffer
    F_INT buffer_offset_;
    /// The first primitive in each contracted function for shells P, Q, R, and S
    F_INT ccbeg_[4];
    /// The last primitive in each contracted function for shells P, Q, R, and S
    F_INT ccend_[4];
    /// Whether to apply screening or not, within ERD
    F_BOOL screen_;
    /// Whether ERD should use spherical harmonic basis functions
    F_BOOL spheric_;
    /// Do any of the basis sets have spherical functions
    bool has_puream_;
    /// Not relating to the monotony of integral computations, but whether the basis sets are all the same
    bool same_bs_;

    void normalize_basis();
public:
    ERDTwoElectronInt(const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
    virtual ~ERDTwoElectronInt();
    void compute_scratch_size();
    virtual size_t compute_shell(const psi::AOShellCombinationsIterator&);
    virtual size_t compute_shell(int, int, int, int);
    virtual size_t compute_shell_deriv1(int, int, int, int);
    virtual size_t compute_shell_deriv2(int, int, int, int);
};

class ERDERI : public ERDTwoElectronInt
{
public:
    ERDERI(const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
    virtual ~ERDERI();
};

}//Namespace
#endif // USING_erd
#endif // header guard
