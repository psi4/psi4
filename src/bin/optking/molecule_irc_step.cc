/*! \file molecule.cc
    \ingroup optking
    \brief molecule class (really, molecular system class)
*/

#include "molecule.h"

#include <cmath>
#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "print.h"
#include "atom_data.h"
#include "physconst.h"

#define EXTERN
#include "globals.h"

namespace opt {

// compute change in energy according to quadratic approximation
inline double DE_nr_energy(double step, double grad, double hess) {
  return (step * grad + 0.5 * step * step * hess);
}

// see molecule_nr_step.cc and molecule_rfo_step.cc for examples from basic step types

// we will assume for now at the IRC implies cartesian coordinates

// Gonzalez and Schlegel JCP 2154 (1989).

void MOLECULE::irc_step(void) {

  fprintf(outfile, "Hello World!");

  double *m = g_masses();

  // Development plan
  // 1. implement Eqn (2) scheme first:  fixed step size s ; step in direction of -1*gradient
  // 2. Muller and Brown (Figure 2)
  // 3. Gonzalez and Schlegel version (Figure 3,4 Eqn. 4-11)
  // 4. Interpolated guess at next point (Eqn. 12-15)

  // add 'step_type = IRC' user keyword
  // see opt_params.h set_params.cc and psi4/src/bin/psi4/read_options.cc (look at 'RFO')

  // modify optking() to execute irc_step if step_type=IRC

  // read gradient/forces  something like
  // double *f_q = p_Opt_data->g_forces_pointer();

  // see if we have any keyword for stepsize.  Don't think so.

  // compute displacements dq Eqn. 2

  // compute back-transformation code to get new cartesians from nr_step

  // symmetrize and save new geometry (same as in nr_step)



  // Need Hessian in mass-weighted cartesian coordinates. ?
  // Read Hessian from opt_data.
  // Mass-weight as follows: H_ixyz_jxyz' = H_ixyz_jxyz/ (sqrt(mass_i) sqrt(mass_j))
  // OPT_DATA::OPT_DATA needs modified to allocate Ncart X Ncart (or 3N by 3N) for
  // cartesian optimizations
  // other functions that read/write/update H will also need modified
  // lets create an H_dim member of optdata class that is dimension of H for functions to use

  // first step is along the eigenvector belonging to the lowest eigenvalue
  


/* copied and pasted from nr_step
  int i, f;
  int Nintco = g_nintco();
  double **H_inv;

  double *f_q = p_Opt_data->g_forces_pointer();
  double **H = p_Opt_data->g_H_pointer();
  double *dq = p_Opt_data->g_dq_pointer();

  double *nr_u;      // unit vector in step direction
  double nr_dqnorm;   // norm of step
  double nr_g;         // gradient in step direction
  double nr_h;         // hessian in step direction
  double DE_projected; // projected energy change by quadratic approximation

  // Hinv f_q = dq
  H_inv = symm_matrix_inv(H, Nintco, 1);
  opt_matrix_mult(H_inv, 0, &f_q, 1, &dq, 1, Nintco, Nintco, 1, 0);
  free_matrix(H_inv);

  // Zero steps for frozen fragment
  for (f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      fprintf(outfile,"\tZero'ing out displacements for frozen fragment %d\n", f+1);
      for (i=0; i<fragments[f]->g_nintco(); ++i)
        dq[ g_intco_offset(f) + i ] = 0.0;
    }
  }

  // applies maximum internal coordinate change
  apply_intrafragment_step_limit(dq);

  // get norm |q| and unit vector in the step direction
  nr_dqnorm = sqrt( array_dot(dq, dq, Nintco) );
  nr_u = init_array(Nintco);
  array_copy(dq, nr_u, Nintco);
  array_normalize(nr_u, Nintco);
  
  // get gradient and hessian in step direction
  nr_g = -1 * array_dot(f_q, nr_u, Nintco); // gradient, not force

  nr_h = 0;
  for (i=0; i<Nintco; ++i)
    nr_h += nr_u[i] * array_dot(H[i], nr_u, Nintco);

  DE_projected = DE_nr_energy(nr_dqnorm, nr_g, nr_h);
  fprintf(outfile,"\tProjected energy change by quadratic approximation: %20.10lf\n", DE_projected);

  // do displacements for each fragment separately
  for (f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      fprintf(outfile,"\tDisplacements for frozen fragment %d skipped.\n", f+1);
      continue;
    }
    fragments[f]->displace(&(dq[g_intco_offset(f)]), true, g_intco_offset(f));
  }

  // do displacements for interfragment coordinates
  double *q_target;
  for (int I=0; I<interfragments.size(); ++I) {

    q_target = interfragments[I]->intco_values();
    for (i=0; i<interfragments[I]->g_nintco(); ++i)
      q_target[i] += dq[g_interfragment_intco_offset(I) + i];

    interfragments[I]->orient_fragment(q_target);

    free_array(q_target);
  }

#if defined(OPTKING_PACKAGE_QCHEM)
  // fix rotation matrix for rotations in QCHEM EFP code
  for (int I=0; I<efp_fragments.size(); ++I)
    efp_fragments[I]->displace( I, &(dq[g_efp_fragment_intco_offset(I)]) );
#endif

  symmetrize_geom(); // now symmetrize the geometry for next step

  // save values in step data
  p_Opt_data->save_step_info(DE_projected, nr_u, nr_dqnorm, nr_g, nr_h);

  free_array(nr_u);
  fflush(outfile);
*/
}

}

