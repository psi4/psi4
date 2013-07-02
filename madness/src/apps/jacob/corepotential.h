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

#ifndef MADNESS_COREPOTENTIAL_H
#define MADNESS_COREPOTENTIAL_H

/// \file corepotential.h
/// \brief Declaration of core potential related class

#include <madness_config.h>
#include <constants.h>
#include <moldft/atomutil.h>
#include <world/print.h>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <sstream>
using std::cout;
using std::endl;
using std::vector;

/// Represents a core potential

/// General Core Potential is able to write down as following form:
/// \f$  U(r) = \sum_k A_k r^(n_k-2) exp(-alpha_k r^2) \sum_m |Y_lm \rangle \langle Y_lm| \f$
/// CorePotential holds these parameters (l,n,A,alpha)
///
/// Note: CorePotential::eval() currently ignores `l'.
///       (It means `\f$\sum_m |Y_lm \rangle \langle Y_lm|\f$' is always `1'.)
struct CorePotential {
    std::vector<int> l; ///< Angular momentum = 0, 1, 2, ...
    std::vector<int> n;
    std::vector<double> A;
    std::vector<double> alpha;
    double eprec, rcut0, rcut;

    CorePotential() : l(), n(), A(), alpha() {};
    CorePotential(const std::vector<int>& l,
                  const std::vector<int>& n,
                  const std::vector<double>& A,
                  const std::vector<double>& alpha)
                : l(l), n(n), A(A), alpha(alpha), eprec(1e-4), rcut0(1.0), rcut(1.0) {};

    double eval(double r) const;

    double eval_derivative(double xi, double r) const;

    std::string to_string () const;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & l & n & A & alpha & eprec & rcut0 & rcut;
    }
};

struct CoreOrbital {
    double Bc;
    int type;
    vector<double> coeff, expnt;
    double rsqmax;

    CoreOrbital() : Bc(0), type(0), coeff(), expnt(), rsqmax(0.0) {}
    CoreOrbital(int type,
                const std::vector<double>& coeff,
                const std::vector<double>& expnt,
                double Bc)
                : Bc(Bc), type(type), coeff(coeff), expnt(expnt)
    {
        double minexpnt = expnt[0];
        for (unsigned int i=1; i<expnt.size(); ++i)
            minexpnt = std::min(minexpnt,expnt[i]);
        rsqmax = 18.4/minexpnt;
    };

    double eval_radial(double rsq) const;

    double eval_radial_derivative(double rsq, double xi) const;

    double eval_spherical_harmonics(int m, double x, double y, double z, double& dp, int axis) const;

    double eval(int m, double rsq, double x, double y, double z) const;

    double eval_derivative(int m, int axis, double xi, double rsq, double x, double y, double z) const;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & Bc & type & rsqmax;
        ar & coeff & expnt;
    }
};

struct AtomCore {
    unsigned int atomic_number;
    unsigned int ncore;
    std::vector<CoreOrbital> orbital;
    CorePotential potential;

    AtomCore() : atomic_number(0), ncore(0), orbital(), potential() {};

    inline unsigned int n_orbital() const { return ncore; };

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & atomic_number & ncore;
        ar & potential;
        for (std::vector<CoreOrbital>::iterator it=orbital.begin(); it != orbital.end(); ++it) {
            ar & (*it);
        }
    }
};

class CorePotentialManager {
    static const double fc;
    std::string core_type;      ///< core potential type (eg. mcp)
    std::string guess_filename; ///< filename of initial guess density data
    std::map<unsigned int,AtomCore> atom_core;
        ///< core potential data mapped atn to potential

public:
    CorePotentialManager() :
        core_type(""),
        guess_filename(""),
        atom_core() {}

    CorePotentialManager(std::string filename, double eprec) : core_type(filename) {
        read_file(filename, std::set<unsigned int>(), eprec);
    }

    inline bool is_defined(const unsigned int atn) const {
        return (atom_core.find(atn) != atom_core.end());
    }

    inline unsigned int n_core_orb(const unsigned int atn) const {
        return (*(atom_core.find(atn))).second.n_orbital();
    }

    inline unsigned int n_core_orb_base(const unsigned int atn) const {
        return (*(atom_core.find(atn))).second.orbital.size();
    }

    inline std::string guess_file() const { return guess_filename; }

    AtomCore get_atom_core(unsigned int atn) const {
        return atom_core.find(atn)->second;
    }

    CorePotential get_potential(unsigned int atn) const {
        return atom_core.find(atn)->second.potential;
    }

    inline unsigned int get_core_l(unsigned int atn, unsigned int core) const {
        return get_atom_core(atn).orbital[core].type;
    }

    inline double get_core_bc(unsigned int atn, unsigned int core) const {
        return get_atom_core(atn).orbital[core].Bc*fc/2;
    }

    inline double core_eval(unsigned int atn, unsigned int core, int m, double rsq, double x, double y, double z) const {
        return get_atom_core(atn).orbital[core].eval(m, rsq, x, y, z);
    }

    inline double core_derivative(unsigned int atn, unsigned int core, int m, int axis, double xi, double rsq, double x, double y, double z) const {
        return get_atom_core(atn).orbital[core].eval_derivative(m, axis, xi, rsq, x, y, z);
    }

    double potential(unsigned int atn, double r) const {
        AtomCore ac = (*(atom_core.find(atn))).second;
        return ac.potential.eval(r);
    }

    double potential_derivative(unsigned int atn, double xi, double r) const {
        AtomCore ac = (*(atom_core.find(atn))).second;
        return ac.potential.eval_derivative(xi, r);
    }

    void read_file(std::string filename, std::set<unsigned int> atomset, double eprec);

    void set_eprec(double value);

    void set_rcut(double value);

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & core_type;
        ar & guess_filename;
        for (std::map<unsigned int, AtomCore>::iterator it=atom_core.begin(); it != atom_core.end(); ++it) {
            ar & it->first & it->second;
        }
    }
};


#endif
