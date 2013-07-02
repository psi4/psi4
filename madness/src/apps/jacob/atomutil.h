/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id$
*/


#ifndef MADNESS_ATOMUTIL_H
#define MADNESS_ATOMUTIL_H

#include <string>

/// \file atomutil.h
/// \brief Declaration of utility class and functions for atom

struct AtomicData {
    // !!! The order of declaration here must match the order in the initializer !!!

    // Nuclear info from L. Visscher and K.G. Dyall, Dirac-Fock
    // atomic electronic structure calculations using different
    // nuclear charge distributions, Atom. Data Nucl. Data Tabl., 67,
    // (1997), 207.
    //
    // http://dirac.chem.sdu.dk/doc/FiniteNuclei/FiniteNuclei.shtml
    const char* const symbol;
    const char* const symbol_lowercase;
    const unsigned int atomic_number;
    const int isotope_number;
    const double nuclear_radius;     ///< Radius of the nucleus for the finite nucleus models (in atomic units).
    const double nuclear_half_charge_radius; ///< Half charge radius in the Fermi Model (in atomic units).
    const double nuclear_gaussian_exponent; ///< Exponential parameter in the Gaussian Model (in atomic units).

    /// Covalent radii stolen without shame from NWChem
    const double covalent_radius;
};

const AtomicData& get_atomic_data(unsigned int atn);

const AtomicData& get_atomic_radius(unsigned int atn);

unsigned int symbol_to_atomic_number(const std::string& symbol);

double smoothing_parameter(double Z, double eprec);
double smoothed_potential(double r);
double dsmoothed_potential(double r);
double smoothed_density(double r);

#endif

