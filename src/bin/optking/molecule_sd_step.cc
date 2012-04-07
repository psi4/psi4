/*! \file molecule_sd_step.cc
    \ingroup optking
    \brief steepest-descent step for molecule 
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
inline double DE_quadratic_energy(double step, double grad, double hess) {
  return (step * grad + 0.5 * step * step * hess);
}

void MOLECULE::sd_step(void) {
  int dim = g_nintco();
  double *fq = p_Opt_data->g_forces_pointer();
  double *dq = p_Opt_data->g_dq_pointer();

  double *last_fq = p_Opt_data->g_last_forces_pointer();
  double h = 1;

  if (last_fq != NULL) {
    // compute overlap of previous forces with current forces
    double *last_fq_u = init_array(dim);
    array_copy(last_fq, last_fq_u, dim);
    array_normalize(last_fq_u, dim);

    double *fq_u = init_array(dim);
    array_copy(fq, fq_u, dim);
    array_normalize(fq_u, dim);

    double fq_overlap = array_dot(last_fq_u, fq_u, dim);
    fprintf(outfile,"\tOverlap of current forces with previous forces %8.4lf\n", fq_overlap);
    free_array(fq_u);
    free_array(last_fq_u);

    if (fq_overlap > 0.90) {
      // component of current forces in step direction (norm of fq)
      double fq_norm = sqrt( array_dot(fq, fq, dim) );
      // component of previous forces in step direction
      double last_fq_norm = array_dot(last_fq, fq, dim) / fq_norm;

      //fprintf(outfile, "fq_norm:        %15.10lf\n", fq_norm);
      //fprintf(outfile, "last_fq_norm:   %15.10lf\n", last_fq_norm);
      //fprintf(outfile, "g_last_dq_norm: %15.10lf\n", p_Opt_data->g_last_dq_norm());
      if (p_Opt_data->g_last_dq_norm() != 0.0)
        h = (last_fq_norm - fq_norm) / p_Opt_data->g_last_dq_norm();

      fprintf(outfile,"\tEstimate of Hessian along step: %10.5e\n", h);
    }
  }

  for (int i=0; i<dim; ++i)
    dq[i] = fq[i] / h;

  // Zero steps for frozen fragment
  for (int f=0; f<fragments.size(); ++f) {
    if (fragments[f]->is_frozen() || Opt_params.freeze_intrafragment) {
      fprintf(outfile,"\tZero'ing out displacements for frozen fragment %d\n", f+1);
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
  fprintf(outfile,"\tProjected energy change: %20.10lf\n", DE_projected);

  // do displacements for each fragment separately
  for (int f=0; f<fragments.size(); ++f) {
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
  p_Opt_data->save_step_info(DE_projected, sd_u, sd_dqnorm, sd_g, sd_h);

  free_array(sd_u);
  fflush(outfile);

}

}

