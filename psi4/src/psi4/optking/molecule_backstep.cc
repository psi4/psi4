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

// Update OPT_DATA to go back one step and try again with a smaller step.
// The current step data is replaced.
// For now we divide the last step size by 1/2 and displace from old geometry.
//
// OPT_DATA contains:
//  left unchanged : Nintco, Ncart
//  left unchanged : Hessian (throw bypasses updating of H)
//  left unchanged : rfo_eigenvector (a unit vector)
//
//  decrease by 1  : iteration (iteration if number of stored steps in opt_data
//                    - at least for now
//  increase by 1  : consecutive_backsteps
//
//  std::vector<STEP_DATA *> steps;
//   Last step data unchanged:
//    fq
//    geom
//    energy
//    unit_step
//    dq_gradient
//    dq_hessian
//
//    dq_norm = dq_norm (last) /2
//    DE_predicted  : DE_nr_energy or DE_rfo_energy (dq_norm, grad, hessian))

// compute change in energy according to quadratic approximation
inline double DE_nr_energy(double step, double grad, double hess) {
  return (step * grad + 0.5 * step * step * hess);
}

// compute change in energy according to RFO approximation
inline double DE_rfo_energy(double rfo_t, double rfo_g, double rfo_h) {
  return (rfo_t * rfo_g + 0.5 * rfo_t * rfo_t * rfo_h)/(1 + rfo_t*rfo_t);
}

void MOLECULE::backstep(void) {

  oprintf_out("\tRe-doing last optimization step - smaller this time.\n");
  oprintf_out("\tConsecutive backstep number %d.\n", p_Opt_data->g_consecutive_backsteps()+1);

  // Erase data created when opt_data was initialized for this step.
  p_Opt_data->erase_last_step();

  p_Opt_data->decrement_iteration();
  p_Opt_data->increment_consecutive_backsteps();

  int Nsteps = p_Opt_data->nsteps();
  int Nintco = Ncoord();

  // Put old cartesian geometry into molecule.
  double *x = p_Opt_data->g_geom_const_pointer(Nsteps-1);
  set_geom_array(x);

  // Compute newly desired dq.
  double *dq = p_Opt_data->g_dq_pointer(Nsteps-1);
  for (int i=0; i<Nintco; ++i)
    dq[i] /= 2;
  double dq_norm = sqrt(array_dot(dq, dq, Nintco));

  oprintf_out("\tNorm of target step-size %10.5lf\n", dq_norm);

  double *rfo_u  = p_Opt_data->g_rfo_eigenvector_pointer();
  double dq_grad = p_Opt_data->g_dq_gradient(Nsteps-1);
  double dq_hess = p_Opt_data->g_dq_hessian(Nsteps-1);

  double DE_projected = 0.0;
  if (Opt_params.step_type == OPT_PARAMS::NR)
    DE_projected = DE_nr_energy(dq_norm, dq_grad, dq_hess);
  else if (Opt_params.step_type == OPT_PARAMS::RFO)
    DE_projected = DE_rfo_energy(dq_norm, dq_grad, dq_hess);
  else if (Opt_params.step_type == OPT_PARAMS::SD)
    DE_projected = DE_nr_energy(dq_norm, dq_grad, dq_hess);

  oprintf_out( "\tNewly projected energy change : %20.10lf\n", DE_projected);

  double *fq = p_Opt_data->g_forces_pointer();

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

#if defined(OPTKING_PACKAGE_QCHEM)
  // fix rotation matrix for rotations in QCHEM EFP code
  for (std::size_t I=0; I<fb_fragments.size(); ++I)
    fb_fragments[I]->displace( I, &(dq[g_fb_fragment_coord_offset(I)]) );
#endif

  symmetrize_geom(); // now symmetrize the geometry for next step

  // save values in step data
  p_Opt_data->save_step_info(DE_projected, rfo_u, dq_norm, dq_grad, dq_hess);


} // end take RFO step

}
