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
  
  $Id: mentity.h 1602 2009-12-27 19:53:06Z rjharrison $
*/
#ifndef MENTITY_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctype.h>
#include <cmath>
#include <misc/misc.h>

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

unsigned int symbol_to_atomic_number(const std::string& symbol);


class Atom {
public:
    double x, y, z, q;          ///< Coordinates and charge in atomic units
    unsigned int atomic_number; ///< Atomic number

    Atom(double x, double y, double z, double q, unsigned int atomic_number)
        : x(x), y(y), z(z), q(q), atomic_number(atomic_number)
    {}

    Atom(const Atom& a)
        : x(a.x), y(a.y), z(a.z), q(a.q), atomic_number(a.atomic_number)
    {}

    /// Default construct makes a zero charge ghost atom at origin
    Atom()
        : x(0), y(0), z(0), q(0), atomic_number(0)
    {}

    template <typename Archive>
    void serialize(Archive& ar) {ar & x & y & z & q & atomic_number;}
};

std::ostream& operator<<(std::ostream& s, const Atom& atom);

class MolecularEntity {
private:
    // If you add more fields don't forget to serialize them
    std::vector<Atom> atoms;
    std::vector<double> rcut;  // Reciprocal of the smoothing radius
    std::vector<double> rsqasymptotic;   // Value of r*r beyond which the potential is asymptotic 1/r

public:
    /// Makes a MolecularEntity with zero atoms
    MolecularEntity() : atoms() {};

    MolecularEntity(const std::string& filename, bool fractional);

    void read_file(const std::string& filename, bool fractional);

    void add_atom(double x, double y, double z, int atn, double q);

    int natom() const {return atoms.size();};

    void set_atom_coords(unsigned int i, double x, double y, double z);

    double bounding_cube() const;

    const Atom& get_atom(unsigned int i) const;

    void print() const;

    double inter_atomic_distance(unsigned int i,unsigned int j) const;

    double nuclear_repulsion_energy() const;

    double smallest_length_scale() const;

    void center();

    double total_nuclear_charge() const;

    double nuclear_attraction_potential(double x, double y, double z) const;

    double nuclear_charge_density(double x, double y, double z) const;

    template <typename Archive>
    void serialize(Archive& ar) {ar & atoms & rcut & rsqasymptotic;}
};

#define MENTITY_H_


#endif /* MENTITY_H_ */
