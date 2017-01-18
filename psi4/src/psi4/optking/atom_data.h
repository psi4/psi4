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

/*!  \file element_to_Z.h
    \ingroup optking
     \brief convert element name or symbol to atomic number
*/

#ifndef _element_to_Z_h_
#define _element_to_Z_h_

#include <iostream>
#include <string>
#include <map>

namespace opt {

// atomic number to the atomic symbol
extern const char *Z_to_symbol[];

// atomic number to default mass
extern const double Z_to_mass[];

// convert symbol or name to atomic number
extern double Element_to_Z(std::string lbl);

// convert isotope name to mass
//extern double Isotope_to_mass(std::string lbl);

}

#endif