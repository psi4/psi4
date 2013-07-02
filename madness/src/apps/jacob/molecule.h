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
#ifndef MADNESS_MOLECULE_H
#define MADNESS_MOLECULE_H

/// \file moldft/molecule.h
/// \brief Declaration of molecule related classes and functions

#include <jacob/corepotential.h>
#include <jacob/atomutil.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <ctype.h>
#include <cmath>
#include <tensor/tensor.h>
#include <misc/misc.h>
#include <world/array.h>


class Atom {
public:
    double x, y, z, q;          ///< Coordinates and charge in atomic units
    unsigned int atomic_number; ///< Atomic number

    explicit Atom(double x, double y, double z, double q, unsigned int atomic_number)
            : x(x), y(y), z(z), q(q), atomic_number(atomic_number) {}

    Atom(const Atom& a)
            : x(a.x), y(a.y), z(a.z), q(a.q), atomic_number(a.atomic_number) {}

    /// Default construct makes a zero charge ghost atom at origin
    Atom()
            : x(0), y(0), z(0), q(0), atomic_number(0) {}

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & x & y & z & q & atomic_number;
    }
};

std::ostream& operator<<(std::ostream& s, const Atom& atom);

class Molecule {
private:
    // If you add more fields don't forget to serialize them
    std::vector<Atom> atoms;
    std::vector<double> rcut;  // Reciprocal of the smoothing radius
    std::vector<double> rsqasymptotic;// value od r*r beyond which the potential is assymptotic 1./r Jacob added
    double eprec;              // Error in energy/atom due to smoothing
    CorePotentialManager core_pot;
    madness::Tensor<double> field;

    void swapaxes(int ix, int iy);

    template <typename opT>
    bool test_for_op(double xaxis, double yaxis, double zaxis, opT op) const;

    bool test_for_c2(double xaxis, double yaxis, double zaxis) const;

    bool test_for_sigma(double xaxis, double yaxis, double zaxis) const;

    bool test_for_inverse() const;

public:
    /// Makes a molecule with zero atoms
    Molecule() : atoms(), rcut(), rsqasymptotic(), eprec(1e-4), core_pot(), field(3L) {};

    void read_file(const std::string& filename);

    void read_core_file(const std::string& filename);

    std::string guess_file() const { return core_pot.guess_file(); };

    unsigned int n_core_orb_all() const ;

    unsigned int n_core_orb(unsigned int atn) const {
        if (core_pot.is_defined(atn))
            return core_pot.n_core_orb_base(atn);
        else
            return 0;
    };

    unsigned int get_core_l(unsigned int atn, unsigned int c) const {
        return core_pot.get_core_l(atn, c);
    }

    double get_core_bc(unsigned int atn, unsigned int c) const {
        return core_pot.get_core_bc(atn, c);
    }

    double core_eval(int atom, unsigned int core, int m, double x, double y, double z) const;

    double core_derivative(int atom, int axis, unsigned int core, int m, double x, double y, double z) const;

    bool is_potential_defined(unsigned int atn) const { return core_pot.is_defined(atn); };

    bool is_potential_defined_atom(int i) const { return core_pot.is_defined(atoms[i].atomic_number); };

    void add_atom(double x, double y, double z,  double q, int atn);

    int natom() const {
        return atoms.size();
    };

    void set_atom_coords(unsigned int i, double x, double y, double z);

    madness::Tensor<double> get_all_coords() const;

    std::vector< madness::Vector<double,3> > get_all_coords_vec() const;

    std::vector<double> atomic_radii;
    
    void set_all_coords(const madness::Tensor<double>& newcoords);


    void set_eprec(double value);

    void set_rcut(double value);

    void set_core_eprec(double value) {
        core_pot.set_eprec(value);
    }

    void set_core_rcut(double value) {
        core_pot.set_rcut(value);
    }

    double get_eprec() const {
        return eprec;
    }

    double bounding_cube() const;

    const Atom& get_atom(unsigned int i) const;

    void print() const;

    double inter_atomic_distance(unsigned int i,unsigned int j) const;

    double nuclear_repulsion_energy() const;

    double nuclear_repulsion_derivative(int i, int j) const;

    double nuclear_dipole(int axis) const;
    
    double nuclear_charge_density(double x, double y, double z) const;//nuclear charge density jacob added

    double smallest_length_scale() const;

    void identify_point_group();

    void center();

    void orient();

    double total_nuclear_charge() const;

    double nuclear_attraction_potential(double x, double y, double z) const;

    double molecular_core_potential(double x, double y, double z) const;

    double core_potential_derivative(int atom, int axis, double x, double y, double z) const;

    double nuclear_attraction_potential_derivative(int atom, int axis, double x, double y, double z) const;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & atoms & rcut & rsqasymptotic & field & eprec & core_pot & atomic_radii;
    }
};


#endif
