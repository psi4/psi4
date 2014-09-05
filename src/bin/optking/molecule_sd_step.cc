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

/*! \file molecule_sd_step.cc
    \ingroup optking
    \brief steepest-descent step for molecule 
*/

#include "molecule.h"

#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "print.h"
#include "atom_data.h"
#include "physconst.h"
#include "psi4-dec.h"
#define EXTERN
#include "globals.h"

#if defined(OPTKING_PACKAGE_PSI)
 #include <cmath>
#elif defined (OPTKING_PACKAGE_QCHEM)
 #include "qcmath.h"
#endif

namespace opt {

// compute change in energy according to quadratic approximation
inline double DE_quadratic_energy(double step, double grad, double hess) {
  return (step * grad + 0.5 * step * step * hess);
}

void MOLECULE::sd_step(void) {
  int dim = g_nintco();
  double *fq = p_Opt_data->g_forces_pointer();
  double *dq = p_Opt_data->g_dq_pointer();

  double *last_fq = p_Opt_data->g_last_forces_pointer();
  double h = 1;
  // now obseleted by addition of cartesian coordinate type
  //bool use_cartesians = true;

  //if (use_cartesians) {
    //sd_step_cartesians();
    //return;
  //}

  if (last_fq != NULL) {
    // compute overlap of previous forces with current forces
    double *last_fq_u = init_array(dim);
    array_copy(last_fq, last_fq_u, dim);
    array_normalize(last_fq_u, dim);

    double *fq_u = init_array(dim);
    array_copy(fq, fq_u, dim);
    array_normalize(fq_u, dim);

    double fq_overlap = array_dot(last_fq_u, fq_u, dim);
    psi::outfile->Printf("\tOverlap of current forces with previous forces %8.4lf\n", fq_overlap);
    free_array(fq_u);
    free_array(last_fq_u);

    if (fq_overlap > 0.90) {
      // component of current forces in step direction (norm of fq)
      double fq_norm = sqrt( array_dot(fq, fq, dim) );
      // component of previous forces in step direction
      double last_fq_norm = array_dot(last_fq, fq, dim) / fq_norm;

      //psi::outfile->Printf( "fq_norm:        %15.10lf\n", fq_norm);
      //psi::outfile->Printf( "last_fq_norm:   %15.10lf\n", last_fq_norm);
      //psi::outfile->Printf( "g_last_dq_norm: %15.10lf\n", p_Opt_data->g_last_dq_norm());
      if (p_Opt_data->g_last_dq_norm() != 0.0)
        h = (last_fq_norm - fq_norm) / p_Opt_data->g_last_dq_norm();

      psi::outfile->Printf("\tEstimate of Hessian along step: %10.5e\n", h);
    }
  }

  for (int i=0; i<dim; ++i)
    dq[i] = fq[i] / h;

  // Zero steps for frozen fragment
  for (int f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      psi::outfile->Printf("\tZero'ing out displacements for frozen fragment %d\n", f+1);
      for (int i=0; i<fragments[f]->g_nintco(); ++i)
        dq[ g_intco_offset(f) + i ] = 0.0;
    }
  }

  apply_intrafragment_step_limit(dq);

  // norm of step
  double sd_dqnorm = sqrt( array_dot(dq, dq, dim) );

  // unit vector in step direction
  double *sd_u = init_array(dim);
  array_copy(dq, sd_u, dim);
  array_normalize(sd_u, dim);

  // gradient in step direction
  double sd_g = - sd_dqnorm;
  double sd_h = 0;

  double DE_projected = DE_quadratic_energy(sd_dqnorm, sd_g, sd_h);
  psi::outfile->Printf("\tProjected energy change: %20.10lf\n", DE_projected);

  // do displacements for each fragment separately
  for (int f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      psi::outfile->Printf("\tDisplacements for frozen fragment %d skipped.\n", f+1);
      continue;
    }
    fragments[f]->displace(&(dq[g_intco_offset(f)]), &(fq[g_intco_offset(f)]), g_atom_offset(f));
  }

  // do displacements for interfragment coordinates
  for (int I=0; I<interfragments.size(); ++I) {
    if (interfragments[I]->is_frozen() || Opt_params.freeze_interfragment) {
      psi::outfile->Printf("\tDisplacements for frozen interfragment %d skipped.\n", I+1);
      continue;
    }
    interfragments[I]->orient_fragment( &(dq[g_interfragment_intco_offset(I)]),
                                        &(fq[g_interfragment_intco_offset(I)]) );
  }

#if defined(OPTKING_PACKAGE_QCHEM)
  // fix rotation matrix for rotations in QCHEM EFP code
  for (int I=0; I<efp_fragments.size(); ++I)
    efp_fragments[I]->displace( I, &(dq[g_efp_fragment_intco_offset(I)]) );
#endif

  symmetrize_geom(); // now symmetrize the geometry for next step

  // save values in step data
  p_Opt_data->save_step_info(DE_projected, sd_u, sd_dqnorm, sd_g, sd_h);

  free_array(sd_u);
  

}

// make primitive for now
void MOLECULE::sd_step_cartesians(void) {

  double *step = g_grad_array();

  // Zero out frozen fragments
  for (int f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
     psi::outfile->Printf("\tZero'ing out displacements for frozen fragment %d\n", f+1);
      for (int i=0; i<fragments[f]->g_natom(); ++i)
        step[ g_intco_offset(f) + i ] = 0.0;
    }
  }

  // Impose step limit
  double limit = Opt_params.intrafragment_step_limit;
  int dim_xyz = 3*g_natom();

  double scale = 1;
  for (int i=0; i<dim_xyz; ++i)
    if ((scale * sqrt(array_dot(step,step,dim_xyz)))> limit)
      scale = limit / sqrt(array_dot(step,step,dim_xyz));

  if (scale != 1.0) {
   psi::outfile->Printf("\tChange in coordinate exceeds step limit of %10.5lf.\n", limit);
   psi::outfile->Printf("\tScaling displacements by %10.5lf\n", scale);
    for (int i=0; i<dim_xyz; ++i)
      step[i] *= scale;
  }

  double *x = g_geom_array();
  for (int i=0; i<dim_xyz; ++i)
    x[i] -= step[i];

  set_geom_array(x);

  double *sd_u = init_array(g_nintco());

  symmetrize_geom(); // now symmetrize the geometry for next step

  double *gx = g_grad_array();
  double DE_projected = -1.0 * array_dot(step, gx, dim_xyz);
 psi::outfile->Printf("\tProjected energy change: %20.10lf\n", DE_projected);
  free_array(gx);

  p_Opt_data->save_step_info(DE_projected, sd_u, 0.0, 0.0, 0.0);

  free_array(sd_u);
}

}

