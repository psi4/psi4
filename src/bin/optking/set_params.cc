/*! \file
    \ingroup OPT10
    \brief set optimization parameters
*/

#define EXTERN
#include "globals.h"

namespace opt {

void set_params(void) {

  Opt_params.step_type = OPT_PARAMS::RFO; // Newton-Raphson (NR) or RFO step

  Opt_params.rfo_follow_root = false; // follow initial step
  Opt_params.rfo_root = 0; // which root to follow
  // eigenvectors of augmented Hessian are divided by the last element unless it is smaller than this value
  Opt_params.rfo_normalization_min = 1.0e-8;

  // bond is assigned if interatomic distance is less than 1.3 * sum of covalent radii
  Opt_params.scale_connectivity = 1.3;

  Opt_params.fragment_mode = OPT_PARAMS::MULTI; // any of {SINGLE, MULTI}

  Opt_params.generate_intcos_only;

  Opt_params.intrafragment_H = OPT_PARAMS::SCHLEGEL;

  Opt_params.interfragment_H = OPT_PARAMS::DEFAULT;

  Opt_params.write_final_step_geometry = false; // if true, use final step geometry


  Opt_params.H_update = OPT_PARAMS::BFGS; // any of {NONE, BFGS, MS, POWELL, BOFILL}
  Opt_params.H_update_use_last = 6; // how many old Hessians to use (0=all)

  // whether to limit changes in Hessian due to update
  Opt_params.H_update_limit = true; 
  // Changes to the Hessian from the update scheme are limited to the larger of
  // (H_update_limit_scale)*(the previous value) and H_update_limit_max (in au).
  Opt_params.H_update_limit_scale = 0.50;
  Opt_params.H_update_limit_max  = 1.0;

  // maximum change in bohr or radian
  Opt_params.intrafragment_step_limit = 0.4;

  // whether to use 1/R(AB) for stretching coordinate between fragments
  Opt_params.interfragment_distance_inverse = false;

  // H-bonds are created if distance is less than this distance in au
  Opt_params.maximum_H_bond_distance = 4.3;

  // QCHEM optimization criteria
  Opt_params.conv_max_force = 3.0e-4;
  Opt_params.conv_max_DE    = 1.0e-6;
  Opt_params.conv_max_disp  = 1.2e-3;

  // tighter convergence criteria
  //Opt_params.conv_max_force = 3.0e-6;
  //Opt_params.conv_max_DE    = 1.0e-8;
  //Opt_params.conv_max_disp  = 1.2e-5;

// ** Unlikely to need modified

  // how close to pi should a torsion be to assume it may have passed through 180
  Opt_params.fix_tors_near_pi = 2.618; // ~150 degrees

  // torsional angles will not be computed if the contained bond angles are within
  // this fraction of pi from 0 or from pi
  Opt_params.tors_angle_lim = 0.01;   

  // only used for determining which atoms in a fragment are acceptable for use
  // as reference atoms.  We avoid collinear sets.
  // angle is 0/pi if the bond angle is within this fraction of pi from 0/pi
  Opt_params.interfrag_collinear_tol = 0.01;   

  // cos(torsional angle) must be this close to -1/+1 for angle to count as 0/pi
  Opt_params.tors_cos_tol = 1e-10;

  // if bend exceeds this value, then also create linear bend complement
  Opt_params.linear_bend_threshold = 3.05; // about 175 degrees

  // threshold for which entries in diagonalized redundant matrix are kept and inverted
  // while computing a generalized inverse of a matrix
  Opt_params.redundant_eval_tol = 1.0e-10;

  Opt_params.bt_max_iter = 25;
  Opt_params.bt_dx_conv = 1.0e-6;
  Opt_params.bt_dx_conv_rms_change = 1.0e-12;

  // 1=default; 2=medium; 3=lots
  Opt_params.print_lvl = 1;

  // Hessian update is avoided if the denominators (Dq*Dq) or (Dq*Dg) are smaller than this
  Opt_params.H_update_den_tol = 1e-7;
  // Hessian update is avoided if any internal coordinate has changed by more than this amount
  // in radians / au
  Opt_params.H_update_dq_tol = 0.5;

  // whether to test B matrix and derivative B matrix numerically
  Opt_params.test_B = true;
  Opt_params.test_derivative_B = false;

}

}

