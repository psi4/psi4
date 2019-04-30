/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here
*/

/* Slater Determinant Class
**
** Based on a previous implementation by David Sherrill using bitsrings
** from the g++ library and on a symbolic Slater's rules program written
** in C by Matt Leininger
**
** February 7, 1996
**
** C. David Sherrill and Matthew L. Leininger
** Center for Computational Quantum Chemistry
** University of Georgia
** Athens, GA 30606
**
** Need to #include <cstdio> before this file
**
** Assume number of alpha electrons is greater than or equal to the number
** of beta electrons
**
** Currently matrix_element() uses static temp arrays which are never free'd.
** The arrays should always be small, so this shouldn't be a problem.  The
** malloc'ing is only done once, no matter how many times matrix_element is
** called.
*/

#ifndef _psi_src_bin_detci_slaterd_h
#define _psi_src_bin_detci_slaterd_h

#include <cstdlib>
#include <string>

namespace psi {
namespace detci {

class CIWavefunction;

class SlaterDeterminant {
    friend class CIWavefunction;

   protected:
    unsigned nalp_;
    unsigned nbet_;
    unsigned char* Occs_[2];

   public:
    SlaterDeterminant() {
        nalp_ = 0;
        nbet_ = 0;
        Occs_[0] = nullptr;
        Occs_[1] = nullptr;
    }
    ~SlaterDeterminant() {
        if (Occs_[0] != nullptr) free(Occs_[0]);
        if (Occs_[1] != nullptr) free(Occs_[1]);
    }
    void set(size_t nalp, unsigned char* alpoccs, size_t nbet, unsigned char* betoccs);
    void print();
    void print_config();
    SlaterDeterminant(const SlaterDeterminant& s);
    SlaterDeterminant& operator=(const SlaterDeterminant& s);
    friend int operator==(SlaterDeterminant& s1, SlaterDeterminant& s2);
    friend double matrix_element(SlaterDeterminant* I, SlaterDeterminant* J);
};

}  // namespace detci
}  // namespace psi

#endif  // header guard
