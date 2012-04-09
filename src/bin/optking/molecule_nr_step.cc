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

void MOLECULE::nr_step(void) {
  int i, f;
  int Nintco = g_nintco();
  double **H_inv;

  double *fq = p_Opt_data->g_forces_pointer();
  double **H = p_Opt_data->g_H_pointer();
  double *dq = p_Opt_data->g_dq_pointer();

  double *nr_u;      // unit vector in step direction
  double nr_dqnorm;   // norm of step
  double nr_g;         // gradient in step direction
  double nr_h;         // hessian in step direction
  double DE_projected; // projected energy change by quadratic approximation

  // Hinv fq = dq
  H_inv = symm_matrix_inv(H, Nintco, 1);
  opt_matrix_mult(H_inv, 0, &fq, 1, &dq, 1, Nintco, Nintco, 1, 0);
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
  nr_g = -1 * array_dot(fq, nr_u, Nintco); // gradient, not force

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
    fragments[f]->displace(&(dq[g_intco_offset(f)]), &(fq[g_intco_offset(f)]), g_atom_offset(f));
  }

  // do displacements for interfragment coordinates
  for (int I=0; I<interfragments.size(); ++I) {
    if (interfragments[I]->is_frozen() || Opt_params.freeze_interfragment) {
      fprintf(outfile,"\tDisplacements for frozen interfragment %d skipped.\n", I+1);
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
  p_Opt_data->save_step_info(DE_projected, nr_u, nr_dqnorm, nr_g, nr_h);

  free_array(nr_u);
  fflush(outfile);
}

}

