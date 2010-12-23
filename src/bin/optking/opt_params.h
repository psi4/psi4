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

  // variable meanings are described in more detail in set_params.cc
  double scale_connectivity;  

  enum FRAGMENT_MODE {SINGLE, MULTI} fragment_mode;

  bool generate_intcos_only;

  bool rfo_follow_root; // whether to do root following
  int rfo_root;         // which root to follow
  double rfo_normalization_min; // small threshold for rfo normalization

  enum {NR, RFO} step_type; // Newton-Raphson (NR) or RFO step

  // Hessian guess

  enum INTRAFRAGMENT_HESSIAN {FISCHER, SCHLEGEL} intrafragment_H;
  enum INTERFRAGMENT_HESSIAN {DEFAULT, FISCHER_LIKE}  interfragment_H;

  enum {NONE, BFGS, MS, POWELL, BOFILL} H_update;
  int H_update_use_last;

  // related to step taken
  double intrafragment_step_limit;

  // whether to limit changes in Hessian due to update
  bool H_update_limit;
  // changes in H are limited to H_update_limit_scale * the previous value
  // if they exceed that, then they are limited to H_update_limit_max
  double H_update_limit_scale;
  double H_update_limit_max;

  // whether to use 1/R as the distance coordinate in interfragment stretching modes
  bool interfragment_distance_inverse;

  double maximum_H_bond_distance;

// ** Unlikely to need modified **

  // how close to pi should a torsion be to assume it may have passed through 180
  double fix_tors_near_pi; 

  // torsional angles will not be computed if the contained bond angles are within
  // this fraction of pi from 0 or from pi
  double tors_angle_lim;

  // as above, how fractionally close to 0 or pi counts as linear when deciding which
  // atoms to avoid using in interfragment coordinates
  double interfrag_collinear_tol;

  double tors_cos_tol; // cos(angle) must be this close to -1/+1 for angle to count as 0/pi

  // threshold for which entries in diagonalized redundant matrix are kept and inverted
  // while computing a generalized inverse of a matrix
  double redundant_eval_tol;

  // maximum number of allowed iterations in backtransformation to cartesian coordinates
  double bt_max_iter;

  // rms and max change in cartesian coordinates in backtransformation
  double bt_dx_conv;

  // give up on backtransformation iterations if change rms from one iteration to the 
  // next is below this value
  double bt_dx_conv_rms_change;

  //1=default; 2=medium; 3=lots
  int print_lvl;

  // Hessian update is avoided if the denominators (Dq*Dq) or (Dq*Dg) are smaller than this
  double H_update_den_tol;
  // Hessian update is avoided if any internal coordinate has changed by more than this amount
  // in radians / au
  double H_update_dq_tol;


  bool test_B; // whether to test B matrices
  bool test_derivative_B; // whether to test derivative B matrices

};

}

#endif

