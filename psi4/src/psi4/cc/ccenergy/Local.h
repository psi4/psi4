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

/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/

#ifndef _psi_src_bin_ccenergy_local_h
#define _psi_src_bin_ccenergy_local_h

#include <string>

namespace psi {
namespace ccenergy {

struct Local {
    int natom;
    int nso;
    int nocc;
    int nvir;
    int *aostart;
    int *aostop;
    int **domain;
    int **pairdomain;
    int *pairdom_len;
    int *pairdom_nrlen;
    int *weak_pairs;
    double ***V;
    double ***W;
    double *eps_occ;
    double **eps_vir;
    double cutoff;
    std::string method;
    std::string weakp;
    int filter_singles;
    double weak_pair_energy;
    double cphf_cutoff;
    int freeze_core;
    std::string pairdef;
};

}  // namespace ccenergy
}  // namespace psi

#endif  // _psi_src_bin_ccenergy_local_h
