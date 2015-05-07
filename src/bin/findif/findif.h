/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef _findif_h_
#define _findif_h_

#include <sstream>
#include <vector>
#include <libmints/mints.h>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

namespace psi { namespace findif {

// functions to generate displacements
std::vector< SharedMatrix > fd_geoms_1_0(Options &options);
std::vector< SharedMatrix > fd_geoms_2_0(Options &options);
std::vector< SharedMatrix > fd_geoms_freq_0(Options &options, int irrep=-1);
std::vector< SharedMatrix > fd_geoms_freq_1(Options &options, int irrep=-1);
std::vector< SharedMatrix > atomic_displacements(Options &options);

// functions to carry out finite-differences
PsiReturnType fd_1_0(Options &options, const boost::python::list& energies);
PsiReturnType fd_2_0(Options &options, const boost::python::list& energies);
PsiReturnType fd_freq_0(Options &options, const boost::python::list& energies, int irrep=-1);
PsiReturnType fd_freq_1(Options &options, const boost::python::list& E_list, int irrep=-1);

// class to accumulate and print vibrations
class VIBRATION {
  int irrep;       // irrep
  double km;    // force constant
  double *lx;   // normal mode in mass-weighted cartesians
  double cm;    // harmonic frequency in wavenumbers

  public:
    friend PsiReturnType fd_freq_0(Options &options, const boost::python::list& energies, int irrep);
    friend PsiReturnType fd_freq_1(Options &options, const boost::python::list& gradients, int irrep);
    friend bool ascending(const VIBRATION *, const VIBRATION *);
    friend void print_vibrations(std::vector<VIBRATION *> modes);

    VIBRATION(int irrep_in, double km_in, double *lx_in) { irrep = irrep_in; km = km_in; lx = lx_in; }
    ~VIBRATION() { free(lx); }
};

// function to print vibrations
void print_vibrations(std::vector<VIBRATION *> modes);

// to order vibrations
bool ascending(const VIBRATION *vib1, const VIBRATION *vib2);

// for displacing along a salc
void displace_cart(SharedMatrix geom, const CdSalcList & salclist,
  int salc_i, int disp_factor, double disp_size);

void displace_cart(SharedMatrix geom, const CdSalcList & salclist,
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

