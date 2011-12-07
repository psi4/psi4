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
//    f_q
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

  fprintf(outfile,"\tRe-doing last optimization step - smaller this time.\n");
  fprintf(outfile,"\tConsecutive backstep number %d.\n", p_Opt_data->g_consecutive_backsteps()+1);

  // Erase data created when opt_data was initialized for this step.
  p_Opt_data->erase_last_step();

  p_Opt_data->decrement_iteration();
  p_Opt_data->increment_consecutive_backsteps();

  int Nsteps = p_Opt_data->nsteps();
  int Nintco = g_nintco();

  // Put old cartesian geometry into molecule.
  double *x = p_Opt_data->g_geom_const_pointer(Nsteps-1);
  set_geom_array(x);

  // Compute newly desired dq.
  double *dq = p_Opt_data->g_dq_pointer(Nsteps-1);
  for (int i=0; i<Nintco; ++i)
    dq[i] /= 2;
  double dq_norm = sqrt(array_dot(dq, dq, Nintco));

  double *rfo_u  = p_Opt_data->g_rfo_eigenvector_pointer();
  double dq_grad = p_Opt_data->g_dq_gradient(Nsteps-1);
  double dq_hess = p_Opt_data->g_dq_hessian(Nsteps-1);

  double DE_projected;
  if (Opt_params.step_type == OPT_PARAMS::NR)
    DE_projected = DE_nr_energy(dq_norm, dq_grad, dq_hess);
  else if (Opt_params.step_type == OPT_PARAMS::RFO)
    DE_projected = DE_rfo_energy(dq_norm, dq_grad, dq_hess);
  else if (Opt_params.step_type == OPT_PARAMS::SD)
    DE_projected = DE_nr_energy(dq_norm, dq_grad, dq_hess);

  fprintf(outfile, "\tNewly projected energy change : %20.10lf\n", DE_projected);

  // do displacements for each fragment separately
  for (int f=0; f<fragments.size(); ++f) {
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
    for (int i=0; i<interfragments[I]->g_nintco(); ++i)
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
  p_Opt_data->save_step_info(DE_projected, rfo_u, dq_norm, dq_grad, dq_hess);

  fflush(outfile);
} // end take RFO step

}
