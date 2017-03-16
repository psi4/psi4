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

#ifndef _findif_h_
#define _findif_h_

#include <sstream>
#include <vector>

#include "psi4/psi4-dec.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"

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

// class to accumulate and print vibrations
class VIBRATION {
  int irrep;       // irrep
  double km;    // force constant
  double *lx;   // normal mode in mass-weighted cartesians
  double cm;    // harmonic frequency in wavenumbers

  public:
    friend PsiReturnType fd_freq_0(Options &options, const py::list& energies, int irrep);
    friend PsiReturnType fd_freq_1(Options &options, const py::list& gradients, int irrep);
    friend bool ascending(const VIBRATION *, const VIBRATION *);
    friend void print_vibrations(std::shared_ptr<Molecule> mol, std::vector<VIBRATION *> modes);

    double get_km() {return km;}
    double get_cm() {return cm;}
    double get_lx(int i) {return lx[i];}

    VIBRATION(int irrep_in, double km_in, double *lx_in) { irrep = irrep_in; km = km_in; lx = lx_in; }
    ~VIBRATION() { free(lx); }
};

// function to print vibrations
void print_vibrations(std::shared_ptr<Molecule> mol, std::vector<VIBRATION *> modes);

// to order vibrations
bool ascending(const VIBRATION *vib1, const VIBRATION *vib2);

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

// save gemetry and normal modes to files
void save_normal_modes(std::shared_ptr<Molecule> mol,
                       std::vector<VIBRATION *> modes);


template <class T>
inline std::string to_string (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}

}}

#endif
