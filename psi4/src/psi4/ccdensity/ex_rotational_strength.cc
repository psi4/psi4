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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/physconst.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#include "psi4/libmints/mintshelper.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

#define _au2cgs 471.44353920

void transdip(MintsHelper &mints);
void transp(MintsHelper &mints, double sign);
void transL(MintsHelper &mints, double sign);

void ex_rotational_strength(MintsHelper &mints, struct TD_Params *S, struct TD_Params *U, struct XTD_Params *xtd_data)
{
  int i, j, k;
  int no, nv, nt;
  double lt_x, lt_y, lt_z;
  double rt_x, rt_y, rt_z;
  double rs_lx, rs_ly, rs_lz;
  double rs_rx, rs_ry, rs_rz;
  double rs_x, rs_y, rs_z;
  double rs;
  double conv;
  double delta_ee;
  int nmo = moinfo.nmo;

  transdip(mints);

  outfile->Printf("\n\tLength-Gauge Rotational Strength for %d%3s to %d%3s Transition\n",S->root+1,
          moinfo.labels[S->irrep], U->root+1, moinfo.labels[U->irrep]);
  outfile->Printf("\t                              X    \t       Y    \t       Z\n");

  lt_x = lt_y = lt_z = 0.0;
  rt_x = rt_y = rt_z = 0.0;
  rs_lx = rs_ly = rs_lz = 0.0;
  rs_rx = rs_ry = rs_rz = 0.0;
  rs_x = rs_y = rs_z = 0.0;

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      lt_x += moinfo.ltd[i][j] * moinfo.dip[0][i][j];
      lt_y += moinfo.ltd[i][j] * moinfo.dip[1][i][j];
      lt_z += moinfo.ltd[i][j] * moinfo.dip[2][i][j];
    }

  transL(mints, +1.0);

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      rt_x += moinfo.rtd[i][j] * moinfo.L[0][i][j];
      rt_y += moinfo.rtd[i][j] * moinfo.L[1][i][j];
      rt_z += moinfo.rtd[i][j] * moinfo.L[2][i][j];
    }

  rs_lx = lt_x * rt_x;
  rs_ly = lt_y * rt_y;
  rs_lz = lt_z * rt_z;

  outfile->Printf("\t<p|mu_e|q>              %11.8lf \t %11.8lf \t %11.8lf\n",
          lt_x,lt_y,lt_z);
  outfile->Printf("\t<q|mu_m|p>              %11.8lf \t %11.8lf \t %11.8lf\n",
          rt_x,rt_y,rt_z);

  // Complex Conjugate

  lt_x = lt_y = lt_z = 0.0;
  rt_x = rt_y = rt_z = 0.0;

  for(i=0; i < 3 ; i++)
    for(j=0; j < nmo; j++)
      for(k=0; k < nmo; k++)
        moinfo.L[i][j][k] = 0.0;

  transL(mints, -1.0);

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      lt_x += moinfo.ltd[i][j] * moinfo.L[0][i][j];
      lt_y += moinfo.ltd[i][j] * moinfo.L[1][i][j];
      lt_z += moinfo.ltd[i][j] * moinfo.L[2][i][j];
    }

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      rt_x += moinfo.rtd[i][j] * moinfo.dip[0][i][j];
      rt_y += moinfo.rtd[i][j] * moinfo.dip[1][i][j];
      rt_z += moinfo.rtd[i][j] * moinfo.dip[2][i][j];
    }

  rs_rx = lt_x * rt_x;
  rs_ry = lt_y * rt_y;
  rs_rz = lt_z * rt_z;

  outfile->Printf("\t<p|mu_m|q>*             %11.8lf \t %11.8lf \t %11.8lf\n",
          lt_x,lt_y,lt_z);
  outfile->Printf("\t<q|mu_e|p>*             %11.8lf \t %11.8lf \t %11.8lf\n",
          rt_x,rt_y,rt_z);

  rs_x = 0.5 * ( rs_lx + rs_rx);
  rs_y = 0.5 * ( rs_ly + rs_ry);
  rs_z = 0.5 * ( rs_lz + rs_rz);

  rs = rs_x + rs_y + rs_z;

  //S->RS_length = rs;
  /* Fill in XTD Data */
  xtd_data->RS_length = rs;

  outfile->Printf("\n");
  outfile->Printf("\tRotational Strength (au)                 %11.8lf\n",rs);
  outfile->Printf("\tRotational Strength (10^-40 esu^2 cm^2)  %11.8lf\n",rs*_au2cgs);


  outfile->Printf("\n\tVelocity-Gauge Rotational Strength for %d%3s\n",S->root+1,
          moinfo.labels[S->irrep], U->root+1, moinfo.labels[U->irrep]);
  outfile->Printf("\t                              X    \t       Y    \t       Z\n");

  lt_x = lt_y = lt_z = 0.0;
  rt_x = rt_y = rt_z = 0.0;
  rs_lx = rs_ly = rs_lz = 0.0;
  rs_rx = rs_ry = rs_rz = 0.0;
  rs_x = rs_y = rs_z = 0.0;

  transp(mints, +1.0);

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      lt_x += moinfo.ltd[i][j] * moinfo.nabla[0][i][j];
      lt_y += moinfo.ltd[i][j] * moinfo.nabla[1][i][j];
      lt_z += moinfo.ltd[i][j] * moinfo.nabla[2][i][j];
    }

  transL(mints, +1.0);

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      rt_x += moinfo.rtd[i][j] * moinfo.L[0][i][j];
      rt_y += moinfo.rtd[i][j] * moinfo.L[1][i][j];
      rt_z += moinfo.rtd[i][j] * moinfo.L[2][i][j];
    }

  rs_lx = lt_x * rt_x;
  rs_ly = lt_y * rt_y;
  rs_lz = lt_z * rt_z;

  outfile->Printf("\t<p|mu_e|q>              %11.8lf \t %11.8lf \t %11.8lf\n",
          lt_x,lt_y,lt_z);
  outfile->Printf("\t<q|mu_m|p>              %11.8lf \t %11.8lf \t %11.8lf\n",
          rt_x,rt_y,rt_z);

  // Complex Conjugate

  lt_x = lt_y = lt_z = 0.0;
  rt_x = rt_y = rt_z = 0.0;

  for(i=0; i < 3 ; i++)
    for(j=0; j < nmo; j++)
      for(k=0; k < nmo; k++) {
        moinfo.nabla[i][j][k] = 0.0;
        moinfo.L[i][j][k] = 0.0;
      }

  transL(mints, -1.0);

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      lt_x += moinfo.ltd[i][j] * moinfo.L[0][i][j];
      lt_y += moinfo.ltd[i][j] * moinfo.L[1][i][j];
      lt_z += moinfo.ltd[i][j] * moinfo.L[2][i][j];
    }

  transp(mints, -1.0);

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      rt_x += moinfo.rtd[i][j] * moinfo.nabla[0][i][j];
      rt_y += moinfo.rtd[i][j] * moinfo.nabla[1][i][j];
      rt_z += moinfo.rtd[i][j] * moinfo.nabla[2][i][j];
    }

  rs_rx = lt_x * rt_x;
  rs_ry = lt_y * rt_y;
  rs_rz = lt_z * rt_z;

  rs_x = 0.5 * ( rs_lx + rs_rx);
  rs_y = 0.5 * ( rs_ly + rs_ry);
  rs_z = 0.5 * ( rs_lz + rs_rz);

  outfile->Printf("\t<p|mu_m|q>*             %11.8lf \t %11.8lf \t %11.8lf\n",
          lt_x,lt_y,lt_z);
  outfile->Printf("\t<q|mu_e|p>*             %11.8lf \t %11.8lf \t %11.8lf\n",
          rt_x,rt_y,rt_z);

  /* Use (w2 - w1) for rotational strengths */
     // Sign matters.  We view excitation energies as positive,
     // so we want to substract the lower state's energy from the
     // higher state's.
     // U should be the higher-energy excited state.
  delta_ee = U->cceom_energy - S->cceom_energy;

  rs_x = rs_x / delta_ee;
  rs_y = rs_y / delta_ee;
  rs_z = rs_z / delta_ee;

  rs = rs_x + rs_y + rs_z;

  /* Fill in XTD Data */
  xtd_data->RS_velocity = rs;

  outfile->Printf("\n");
  outfile->Printf("\tRotational Strength (au)                 %11.8lf\n",rs);
  outfile->Printf("\tRotational Strength (10^-40 esu^2 cm^2)  %11.8lf\n",rs*_au2cgs);


  return;
}

}} // namespace psi::ccdensity
