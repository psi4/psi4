/*! \file opt-params.h
    \ingroup OPT10
    \brief header for optimization parameters
*/

#ifndef _opt_opt_params_h_
#define _opt_opt_params_h_

namespace opt {

struct OPT_PARAMS {

  // convergence criteria
  double conv_max_force;
  double conv_max_DE;
  double conv_max_disp;

  // variable meanings are described in set_params.cc
  double scale_connectivity;  

  // whether to do root following
  bool rfo_follow_root;
  // default= 1
  int rfo_root;

  enum {NR, RFO} step_type; // Newton-Raphson (NR) or RFO step

  // Hessian guess

  enum INTRAFRAGMENT_HESSIAN {FISCHER, SCHLEGEL} intrafragment_H;
  enum INTERFRAGMENT_HESSIAN {DEFAULT, FISCHER_LIKE}  interfragment_H;

  enum {NONE, BFGS, MS, POWELL, BOFILL} H_update;
  int H_update_use_last;

  // related to step taken
  double intrafragment_step_limit;
  // maximum change in aJ/Ang^2 allowed by Hessian updates
  double H_change_limit;

  bool frag_dist_rho;

// ** Unlikely to need modified **

  // how close to pi should a torsion be to assume it may have passed through 180
  double fix_tors_near_pi; 

  // torsional angles will not be computed if the contained bond angles are within
  // this fraction of pi from 0 or from pi
  double error_tors_angle;  

  double collinear_lim; // as above, how close to 0 or pi counts as linear for angles

  // threshold for which entries in diagonalized redundant matrix are kept and inverted
  // while computing a generalized inverse of a matrix
  double redundant_eval_tol;

  // number of allowed iterations in backtransformation to cartesian coordinates
  double bt_max_iter;
  // required convergence for backtransformation on internal coordinates
  double bt_dq_conv_rms;
  double bt_dq_conv_max;

  //1=default; 2=medium; 3=lots
  int print_lvl;

};

}

#endif

