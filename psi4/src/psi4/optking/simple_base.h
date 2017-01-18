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

/*! \file simple_base.h
    \ingroup optking
    \brief simple internal coordinate base class
*/

#ifndef _opt_simple_base_h_
#define _opt_simple_base_h_

#include <cstdio>
#include <string>
#include "mem.h"

namespace opt {

enum INTCO_TYPE {min_type, stre_type, bend_type, tors_type, oofp_type, cart_type, max_type};

typedef const double * const * const GeomType;

class SIMPLE_COORDINATE {

  protected:

    INTCO_TYPE s_type; // type of simple
    int s_natom;       // # of atoms in internal definition, i.e., # of
                       // atoms with non-zero s-vectors
    int *s_atom;       // atom indices in internal definition

    bool s_frozen;     // is internal coordinate constrained frozen?

    double s_fixed_eq_val; // for adding artificial forces
    bool s_has_fixed_eq_val;

  public:

    SIMPLE_COORDINATE (INTCO_TYPE s_type_in, int s_natom_in, bool freeze_in) {
      s_type = s_type_in;
      s_natom = s_natom_in;
      s_atom = init_int_array(s_natom);
      s_frozen = freeze_in;
      s_has_fixed_eq_val = false;
    };

    virtual ~SIMPLE_COORDINATE() { // derived class destructors called first; then this one
      free_int_array(s_atom);
    }

    INTCO_TYPE g_type(void) const { return s_type; }
    int g_natom(void) const { return s_natom; }
    int g_atom(int a) const { return s_atom[a]; }

    bool is_frozen(void) { return s_frozen; }
    void freeze(void)    { s_frozen = true; }
    void unfreeze(void)  { s_frozen = false; }

    // do-nothing function overridden only by torsion class
    virtual void fix_tors_near_180(GeomType) { return; }

    // do-nothing function overridden only by oofp class
    virtual void fix_oofp_near_180(GeomType) { return; }

    // do-nothing function overridden by stretch class 
    virtual bool is_hbond(void) const { return false; }

    // do-nothing function overridden by stretch class 
    virtual void set_hbond(bool) { printf("base_hbond"); return; }

    // do-nothing function overridden by stretch class 
    virtual bool is_inverse_stre(void) const { return false; }

    // do-nothing function overridden by bend class
    virtual int g_bend_type(void) const { return -1; }

    // do-nothing function overridden by cartesian class
    virtual int g_xyz(void) const { return 0; }

    // each internal coordinate type must provide the following virtual functions:

    // function to print coordinate definitions to intco.dat
    virtual void print_intco_dat(std::string psi_fp, FILE *qc_fc, int atom_offset=0) const = 0;

    // return value of internal coordinate
    virtual double value(GeomType geom) const = 0;

    // compute s vector (dq/dx, dq/dy, dq/dz)
    virtual double ** DqDx(GeomType geom) const = 0;

    // compute second derivative (B' matrix elements)
    // dq_i dq_j / dx^2, dydx, dy^2, dzdx, dzdy, dz^2 
    virtual double ** Dq2Dx2(GeomType geom) const = 0;

    // print coordinates and value to output file
    virtual void print(std::string psi_fp, FILE *qc_fp, GeomType geom, int atom_offset=0) const = 0;

    // print coordinates and displacements 
    virtual void print_disp(std::string psi_fp, FILE *qc_fp, const double old_q, const double f_q,
      const double dq, const double new_q, int atom_offset = 0) const = 0;

    // for debugging, print s vectors to output file
    virtual void print_s(std::string psi_fp, FILE *qc_fp, GeomType geom) const = 0;

    // function to return string of coordinate definition
    virtual std::string get_definition_string(int atom_offset=0) const = 0;

    // each derived class should have an equality operator of this type that
    // first checks types with g_type() and then, if true, compares
    // the internal coordinates
    virtual bool operator==(const SIMPLE_COORDINATE & s2) const  = 0;

    // functions to save user-specified constrained eq_val
    void set_fixed_eq_val(double val) {
      s_fixed_eq_val = val;
      s_has_fixed_eq_val = true;
    }
    double fixed_eq_val(void) {
      return s_fixed_eq_val;
    }
    bool has_fixed_eq_val(void) { return s_has_fixed_eq_val; }

};

}

#endif