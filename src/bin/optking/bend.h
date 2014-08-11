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

/*! \file bend.h
    \ingroup OPTKING
    \brief BEND class declaration
*/

#ifndef _opt_bend_h_
#define _opt_bend_h_

#include "simple.h"

namespace opt {

class BEND : public SIMPLE {

    // if true, this bend is a secondary complement to another linear bend
    bool linear_bend;

  public:

    BEND(int A_in, int B_in, int C_in, bool freeze_in=false);

    ~BEND() { } // also calls ~SIMPLE()

    double value(GeomType geom) const;

    // compute and return array of first derivative (B marix elements)
    double **DqDx(GeomType geom) const;

    // compute and return array of second derivative (B' matrix elements)
    double **Dq2Dx2(GeomType geom) const;

    void print(std::string OutFileRMR, GeomType geom, int atom_offset=0) const;
    void print_intco_dat(std::string OutFileRMR, int atom_offset=0) const;
    void print_s(std::string OutFileRMR, GeomType geom) const;
    void print_disp(std::string OutFileRMR, const double old_q, const double f_q,
      const double dq, const double new_q, int atom_offset=0) const;
    bool operator==(const SIMPLE & s2) const;
    std::string get_definition_string(int atom_offset=0) const;

    void make_linear_bend(void) { linear_bend = true; }
    bool is_linear_bend(void) const { return linear_bend; }

};

}

#endif

