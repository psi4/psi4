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

/*! \file molecule.cc
    \ingroup optking
    \brief molecule class (really, molecular system class)
*/

#include "molecule.h"

#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "atom_data.h"
#include "psi4/optking/physconst.h"

#include "print.h"
#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

namespace opt {

// compute change in energy according to quadratic approximation
inline double DE_nr_energy(double step, double grad, double hess) {
  return (step * grad + 0.5 * step * step * hess);
}

void MOLECULE::nr_step(void) {
  int Nintco = Ncoord();
  double **H_inv;

  double *fq = p_Opt_data->g_forces_pointer();
  double **H = p_Opt_data->g_H_pointer();
  double *dq = p_Opt_data->g_dq_pointer();

  double *nr_u;      // unit vector in step direction
  double nr_dqnorm;   // norm of step
  double nr_g;         // gradient in step direction
  double nr_h;         // hessian in step direction
  double DE_projected; // projected energy change by quadratic approximation

  oprintf_out("\tTaking NR optimization step.\n");

  // Hinv fq = dq
  H_inv = symm_matrix_inv(H, Nintco, 1);
  opt_matrix_mult(H_inv, 0, &fq, 1, &dq, 1, Nintco, Nintco, 1, 0);
  free_matrix(H_inv);

  // Zero steps for frozen fragment
  for (std::size_t f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      oprintf_out("\tZero'ing out displacements for frozen fragment %d\n", f+1);
      for (int i=0; i<fragments[f]->Ncoord(); ++i)
        dq[ g_coord_offset(f) + i ] = 0.0;
    }
  }

  // applies maximum internal coordinate change
  apply_intrafragment_step_limit(dq);

  // get norm |q| and unit vector in the step direction
  nr_dqnorm = sqrt( array_dot(dq, dq, Nintco) );
  nr_u = init_array(Nintco);
  array_copy(dq, nr_u, Nintco);
  array_normalize(nr_u, Nintco);

  oprintf_out("\tNorm of target step-size %10.5lf\n", nr_dqnorm);

  // get gradient and hessian in step direction
  nr_g = -1 * array_dot(fq, nr_u, Nintco); // gradient, not force

  nr_h = 0;
  for (int i=0; i<Nintco; ++i)
    nr_h += nr_u[i] * array_dot(H[i], nr_u, Nintco);

  DE_projected = DE_nr_energy(nr_dqnorm, nr_g, nr_h);
  oprintf_out("\tProjected energy change by quadratic approximation: %20.10lf\n", DE_projected);

  // do displacements for each fragment separately
  for (std::size_t f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      oprintf_out("\tDisplacements for frozen fragment %d skipped.\n", f+1);
      continue;
    }
    fragments[f]->displace(&(dq[g_coord_offset(f)]), &(fq[g_coord_offset(f)]), g_atom_offset(f));
  }

  // do displacements for interfragment coordinates
  for (std::size_t I=0; I<interfragments.size(); ++I) {
    if (interfragments[I]->is_frozen() || Opt_params.freeze_interfragment) {
      oprintf_out("\tDisplacements for frozen interfragment %d skipped.\n", I+1);
      continue;
    }
    interfragments[I]->orient_fragment( &(dq[g_interfragment_coord_offset(I)]),
                                        &(fq[g_interfragment_coord_offset(I)]) );
  }

  // fix rotation matrix for rotations in QCHEM EFP code
  for (std::size_t I=0; I<fb_fragments.size(); ++I)
    fb_fragments[I]->displace( I, &(dq[g_fb_fragment_coord_offset(I)]) );

  symmetrize_geom(); // now symmetrize the geometry for next step

  // save values in step data
  p_Opt_data->save_step_info(DE_projected, nr_u, nr_dqnorm, nr_g, nr_h);

  free_array(nr_u);

}

}
