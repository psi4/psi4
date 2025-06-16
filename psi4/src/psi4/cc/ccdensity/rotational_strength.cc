/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/cc/ccwave.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/physconst.h"
#include "ccdensity.h"
#include "Frozen.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {

#define _au2cgs 471.44353920

void transdip(const MintsHelper &mints);
void transp(const MintsHelper &mints, double sign);
void transL(const MintsHelper &mints, double sign);

void rotational_strength(ccenergy::CCEnergyWavefunction& wfn, struct TD_Params *S) {
    int i, j, k;
    int no, nv, nt;
    double rs_lx, rs_ly, rs_lz;
    double rs_rx, rs_ry, rs_rz;
    double rs_x, rs_y, rs_z;
    double rs;
    double conv;
    int nmo = moinfo.nmo;
    const auto& mints = *wfn.mintshelper();

    transdip(mints);

    outfile->Printf("\n\tLength-Gauge Rotational Strength for %d%3s\n", S->root + 1, moinfo.labels[S->irrep].c_str());
    outfile->Printf("\t                              X    \t       Y    \t       Z\n");

    rs_lx = rs_ly = rs_lz = 0.0;
    rs_rx = rs_ry = rs_rz = 0.0;
    rs_x = rs_y = rs_z = 0.0;

    auto lt_x = moinfo.ltd_mat.vector_dot(moinfo.dip[0]);
    auto lt_y = moinfo.ltd_mat.vector_dot(moinfo.dip[1]);
    auto lt_z = moinfo.ltd_mat.vector_dot(moinfo.dip[2]);

    transL(mints, +1.0);

    auto rt_x = moinfo.rtd_mat.vector_dot(moinfo.L[0]);
    auto rt_y = moinfo.rtd_mat.vector_dot(moinfo.L[1]);
    auto rt_z = moinfo.rtd_mat.vector_dot(moinfo.L[2]);

    rs_lx = lt_x * rt_x;
    rs_ly = lt_y * rt_y;
    rs_lz = lt_z * rt_z;

    outfile->Printf("\t<0|mu_e|n>              %11.8lf \t %11.8lf \t %11.8lf\n", lt_x, lt_y, lt_z);
    outfile->Printf("\t<n|mu_m|0>              %11.8lf \t %11.8lf \t %11.8lf\n", rt_x, rt_y, rt_z);

    // Complex Conjugate

    transL(mints, -1.0);

    lt_x = moinfo.ltd_mat.vector_dot(moinfo.L[0]);
    lt_y = moinfo.ltd_mat.vector_dot(moinfo.L[1]);
    lt_z = moinfo.ltd_mat.vector_dot(moinfo.L[2]);

    rt_x = moinfo.rtd_mat.vector_dot(moinfo.dip[0]);
    rt_y = moinfo.rtd_mat.vector_dot(moinfo.dip[1]);
    rt_z = moinfo.rtd_mat.vector_dot(moinfo.dip[2]);

    rs_rx = lt_x * rt_x;
    rs_ry = lt_y * rt_y;
    rs_rz = lt_z * rt_z;

    outfile->Printf("\t<0|mu_m|n>*             %11.8lf \t %11.8lf \t %11.8lf\n", lt_x, lt_y, lt_z);
    outfile->Printf("\t<n|mu_e|0>*             %11.8lf \t %11.8lf \t %11.8lf\n", rt_x, rt_y, rt_z);

    rs_x = 0.5 * (rs_lx + rs_rx);
    rs_y = 0.5 * (rs_ly + rs_ry);
    rs_z = 0.5 * (rs_lz + rs_rz);

    rs = rs_x + rs_y + rs_z;
    S->RS_length = rs;

    outfile->Printf("\n");
    outfile->Printf("\tRotational Strength (au)                 %11.8lf\n", rs);
    outfile->Printf("\tRotational Strength (10^-40 esu^2 cm^2)  %11.8lf\n", rs * _au2cgs);

    // Save rotatory strength to wfn.
    scalar_saver_ground(wfn, S, "ROTATORY STRENGTH (LEN)", rs);

    outfile->Printf("\n\tVelocity-Gauge Rotational Strength for %d%3s\n", S->root + 1, moinfo.labels[S->irrep].c_str());
    outfile->Printf("\t                              X    \t       Y    \t       Z\n");

    rs_lx = rs_ly = rs_lz = 0.0;
    rs_rx = rs_ry = rs_rz = 0.0;
    rs_x = rs_y = rs_z = 0.0;

    transp(mints, +1.0);

    lt_x = moinfo.ltd_mat.vector_dot(moinfo.nabla[0]);
    lt_y = moinfo.ltd_mat.vector_dot(moinfo.nabla[1]);
    lt_z = moinfo.ltd_mat.vector_dot(moinfo.nabla[2]);

    transL(mints, +1.0);

    rt_x = moinfo.rtd_mat.vector_dot(moinfo.L[0]);
    rt_y = moinfo.rtd_mat.vector_dot(moinfo.L[1]);
    rt_z = moinfo.rtd_mat.vector_dot(moinfo.L[2]);

    rs_lx = lt_x * rt_x;
    rs_ly = lt_y * rt_y;
    rs_lz = lt_z * rt_z;

    outfile->Printf("\t<0|mu_e|n>              %11.8lf \t %11.8lf \t %11.8lf\n", lt_x, lt_y, lt_z);
    outfile->Printf("\t<n|mu_m|0>              %11.8lf \t %11.8lf \t %11.8lf\n", rt_x, rt_y, rt_z);

    // Complex Conjugate
    transL(mints, -1.0);

    lt_x = moinfo.ltd_mat.vector_dot(moinfo.L[0]);
    lt_y = moinfo.ltd_mat.vector_dot(moinfo.L[1]);
    lt_z = moinfo.ltd_mat.vector_dot(moinfo.L[2]);

    transp(mints, -1.0);

    rt_x = moinfo.rtd_mat.vector_dot(moinfo.nabla[0]);
    rt_y = moinfo.rtd_mat.vector_dot(moinfo.nabla[1]);
    rt_z = moinfo.rtd_mat.vector_dot(moinfo.nabla[2]);

    rs_rx = lt_x * rt_x;
    rs_ry = lt_y * rt_y;
    rs_rz = lt_z * rt_z;

    rs_x = 0.5 * (rs_lx + rs_rx);
    rs_y = 0.5 * (rs_ly + rs_ry);
    rs_z = 0.5 * (rs_lz + rs_rz);

    outfile->Printf("\t<0|mu_m|n>*             %11.8lf \t %11.8lf \t %11.8lf\n", lt_x, lt_y, lt_z);
    outfile->Printf("\t<n|mu_e|0>*             %11.8lf \t %11.8lf \t %11.8lf\n", rt_x, rt_y, rt_z);

    rs_x = rs_x / S->cceom_energy;
    rs_y = rs_y / S->cceom_energy;
    rs_z = rs_z / S->cceom_energy;

    rs = rs_x + rs_y + rs_z;
    S->RS_velocity = rs;

    outfile->Printf("\n");
    outfile->Printf("\tRotational Strength (au)                 %11.8lf\n", rs);
    outfile->Printf("\tRotational Strength (10^-40 esu^2 cm^2)  %11.8lf\n", rs * _au2cgs);

    // Save rotatory strength to wfn.
    scalar_saver_ground(wfn, S, "ROTATORY STRENGTH (VEL)", rs);

    return;
}
}
}  // namespace psi
