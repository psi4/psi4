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

  Opt_params.scale_connectivity = 1.2;

  Opt_params.empirical_H = OPT_PARAMS::SCHLEGEL;
  //Opt_params.empirical_H = OPT_PARAMS::FISCHER;

  //Opt_params.H_update = OPT_PARAMS::BOFILL; // any of {NONE, BFGS, MS, POWELL, BOFILL}
  Opt_params.H_update = OPT_PARAMS::BFGS; // any of {NONE, BFGS, MS, POWELL, BOFILL}
  Opt_params.H_update_use_last = 0; // how many old Hessians to use (0=all)

  // maximum change in bohr or radian
  Opt_params.intrafragment_step_limit = 0.4;

  // maximum change in Hessian element allowed by Hessian update schemes is this or 50% of value
  Opt_params.H_change_limit = 0.3;

  // QCHEM optimization criteria
  Opt_params.conv_max_force = 3.0e-4;
  Opt_params.conv_max_DE    = 1.0e-6;
  Opt_params.conv_max_disp  = 1.2e-3;

// ** Unlikely to need modified

  // how close to pi should a torsion be to assume it may have passed through 180
  Opt_params.fix_tors_near_pi = 2.618; // ~150 degrees

  // torsional angles will not be computed if the contained bond angles are within
  // this fraction of pi from 0 or from pi
  Opt_params.error_tors_angle = 0.01;   

  // threshold for which entries in diagonalized redundant matrix are kept and inverted
  // while computing a generalized inverse of a matrix
  Opt_params.redundant_eval_tol = 1.0e-10;

  Opt_params.bt_max_iter = 25;
  Opt_params.bt_dq_conv_rms = 1.0e-7;
  Opt_params.bt_dq_conv_max = 1.0e-5;

  // 1=default; 2=medium; 3=lots
  Opt_params.print_lvl = 1;

}

}

