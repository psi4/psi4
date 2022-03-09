/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

/*! \file stre.h
    \ingroup optking
    \brief STRE class declaration
*/

#ifndef _opt_stretch_h_
#define _opt_stretch_h_

#include "simple_base.h"

namespace opt {

class STRE : public SIMPLE_COORDINATE {

    bool hbond; // whether stretch is a hydrogen bond
    bool inverse_stre; // whether stretch is really 1/R

  public:

    STRE(int A_in, int B_in, bool freeze_in=false);

    ~STRE() override { } // also calls ~SIMPLE_COORDINATE()

    double value(GeomType geom) const override;

    // compute and return array of first derivative (B matrix elements)
    // returned matrix is [atom][x,y,z]
    double **DqDx(GeomType geom) const override;

    // compute and return array of second derivative (B' matrix elements)
    // returned matrix is order 3N cart by 3N cart
    double **Dq2Dx2(GeomType geom) const override;

    void print(std::string psi_fp, FILE *qc_fp, GeomType geom, int atom_offset=0) const override;
    void print_intco_dat(std::string psi_fp, FILE *qc_fp, int atom_offset=0) const override;
    void print_s(std::string psi_fp, FILE *qc_fp, GeomType geom) const override;
    void print_disp(std::string psi_fp, FILE *qc_fp, const double old_q, const double f_q,
      const double dq, const double new_q, int atom_offset = 0) const override;
    bool operator==(const SIMPLE_COORDINATE & s2) const override;
    std::string get_definition_string(int atom_offset=0) const override;

    void set_hbond(bool val) override { hbond = val; }
    bool is_hbond() const override { return hbond; }
    void make_inverse_stre() { inverse_stre = true; }
    bool is_inverse_stre() const override { return inverse_stre; }

};

}

#endif
