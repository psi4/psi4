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

#ifndef libmints_typedefs_h
#define libmints_typedefs_h

#include <memory> // for shared_ptr

// Handy mints timer macros, requires libqt to be included
#ifdef MINTS_TIMER
#   include "psi4/libqt/qt.h"
#   define mints_timer_on(a) timer_on((a));
#   define mints_timer_off(a) timer_off((a));
#else
#   define mints_timer_on(a)
#   define mints_timer_off(a)
#endif

// Forward declare psi
namespace psi {
class Matrix;
class Vector;
class Wavefunction;
using SharedMatrix = std::shared_ptr<Matrix>;
using SharedVector = std::shared_ptr<Vector>;
using SharedWavefunction =  std::shared_ptr<Wavefunction>;

// Useful when working with SO-TEIs
template<typename T>
void swap_index(T& a, T& b) {
    T temp;
    temp = b;
    b = a;
    a = temp;
}

#define SWAP_INDEX(a, b) swap_index(a ## abs, b ## abs); swap_index(a ## rel, b ## rel); swap_index(a ## irrep, b ## irrep);

}


#endif // libmints_typedefs_h
