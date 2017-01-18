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

/*! \file tors.h
    \ingroup optking
    \brief TORS class declaration
*/

#ifndef _opt_tors_h_
#define _opt_tors_h_

#include "simple_base.h"

namespace opt {

class TORS : public SIMPLE_COORDINATE {

  private:
    int near_180;
        // +1 if positive and approaching 180
        // -1 if negative and approaching -180
        //  0 otherwise

  public:

    TORS(int A_in, int B_in, int C_in, int D_in, bool freeze_in=false);

    ~TORS() { } // also calls ~SIMPLE_COORDINATE

    double value(GeomType geom) const;
    //bool fix_tors_value_corrected(GeomType geom) const;

    // compute and return array of first derivative (B marix elements)
    double ** DqDx(GeomType geom) const;

    // compute and return array of second derivative (B' matrix elements)
    double ** Dq2Dx2(GeomType geom) const;

    void print(std::string psi_fp, FILE *qc_fp, GeomType geom, int atom_offset=0) const;
    void print_intco_dat(std::string psi_fp, FILE *qc_fp, int atom_offset=0) const;
    void print_s(std::string psi_fp, FILE *qc_fp, GeomType geom) const;
    void print_disp(std::string psi_fp, FILE *qc_fp, const double old_q, const double f_q, 
      const double dq, const double new_q, int atom_offset=0) const;
    bool operator==(const SIMPLE_COORDINATE & s2) const;
    std::string get_definition_string(int atom_offset=0) const;

    void fix_tors_near_180(GeomType geom);

    //bool check_tors_for_bad_angles(GeomType geom) const;

};

}

#endif