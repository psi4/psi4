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

#ifndef _findif_h_
#define _findif_h_

#include <sstream>
#include <vector>

#include "psi4/psi4-dec.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/typedefs.h"
#include "psi4/pybind11.h"
#include "psi4/libmints/vector.h"

namespace psi {
class Molecule;
class CdSalcList;

namespace findif {

// functions to generate displacements
std::vector< SharedMatrix > fd_geoms_1_0(std::shared_ptr<Molecule> mol, Options &options);
// std::vector< SharedMatrix > fd_geoms_2_0(Options &options);
std::vector< SharedMatrix > fd_geoms_freq_0(std::shared_ptr<Molecule> mol, Options &options, int irrep=-1);
std::vector< SharedMatrix > fd_geoms_freq_1(std::shared_ptr<Molecule> mol, Options &options, int irrep=-1);
std::vector< SharedMatrix > atomic_displacements(std::shared_ptr<Molecule> mol, Options &options);

// functions to carry out finite-differences
SharedMatrix fd_1_0(std::shared_ptr<Molecule> mol, Options &options, const py::list& energies);
//PsiReturnType fd_2_0(std::shared_ptr<Molecule> mol, Options &options, const py::list& energies);
SharedMatrix fd_freq_0(std::shared_ptr<Molecule> mol, Options &options,
                      const py::list& energies, int irrep=-1);
SharedMatrix fd_freq_1(std::shared_ptr<Molecule> mol, Options &options,
                      const py::list& E_list, int irrep=-1);

// for displacing along a salc
void displace_cart(std::shared_ptr<Molecule> mol, SharedMatrix geom, const CdSalcList & salclist,
  int salc_i, int disp_factor, double disp_size);

void displace_cart(std::shared_ptr<Molecule> mol, SharedMatrix geom, const CdSalcList & salclist,
  int salc_i, int salc_j, int disp_factor_i, int disp_factor_j, double disp_size);

// to massweight columns of a shared matrix
void mass_weight_columns_plus_one_half(SharedMatrix B);

// displace an atomic coordinate
void displace_atom(SharedMatrix geom, const int atom, const int coord,
                   const int sign, const double disp_size);

template <class T>
inline std::string to_string (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}

}}

#endif
