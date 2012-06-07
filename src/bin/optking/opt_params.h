/*! \file opt-params.h
    \ingroup optking
    \brief header for optimization parameters
      variable meanings are described in more detail in set_params.cc
*/

#ifndef _opt_opt_params_h_
#define _opt_opt_params_h_

namespace opt {

struct OPT_PARAMS {

  // convergence criteria
  double conv_max_force;
  double conv_rms_force;
  double conv_max_DE;
  double conv_max_disp;
  double conv_rms_disp;
  bool i_max_force;
  bool i_rms_force;
  bool i_max_DE;
  bool i_max_disp;
  bool i_rms_disp;
  bool i_untampered;
  std::string general_conv;

  double scale_connectivity;

  enum FRAGMENT_MODE {SINGLE, MULTI} fragment_mode;

  enum INTERFRAGMENT_MODE {FIXED, PRINCIPAL_AXES} interfragment_mode;

  bool generate_intcos_only;

  bool rfo_follow_root; // whether to do root following
  int rfo_root;         // which root to follow
  double rfo_normalization_min; // small threshold for rfo normalization

  enum OPT_TYPE {MIN, TS, IRC} opt_type;
  // Newton-Raphson (NR), rational function optimization step, steepest descent step
  enum STEP_TYPE {NR, RFO, P_RFO, SD} step_type;

  // Hessian guess
  // Note the Lindh "intrafragment" option is cartesian so it applies to all coordinates.
  enum INTRAFRAGMENT_HESSIAN {FISCHER, SCHLEGEL, SIMPLE, LINDH} intrafragment_H;
  enum INTERFRAGMENT_HESSIAN {DEFAULT, FISCHER_LIKE}  interfragment_H;

  enum H_UPDATE {NONE, BFGS, MS, POWELL, BOFILL} H_update;
  int H_update_use_last;

  enum IRC_DIRECTION {FORWARD, BACKWARD} IRC_direction;
  enum IRC_STOP {ASK, STOP, GO} IRC_stop;

  bool freeze_intrafragment; // freeze all fragments
  bool freeze_interfragment; // freeze all interfragment modes
  bool add_auxiliary_bonds;

  // related to step taken
  double intrafragment_step_limit;
  double intrafragment_step_limit_min;
  double intrafragment_step_limit_max;

  double interfragment_step_limit;

  // whether to limit changes in Hessian due to update
  bool H_update_limit;
  // changes in H are limited to H_update_limit_scale * the previous value
  // if they exceed that, then they are limited to H_update_limit_max
  double H_update_limit_scale;
  double H_update_limit_max;

  // whether to use 1/R as the distance coordinate in interfragment stretching modes
  bool interfragment_distance_inverse;

  // By default, optking prints and saves the last (previous) geometry at the end of an
  // optimization, i.e., the one at which a gradient was computed.
  // If true, then the structure obtained from the last anticipated step is printed and saved instead.
  bool write_final_step_geometry;

  double maximum_H_bond_distance;

  bool read_cartesian_H;

  bool efp_fragments;
  bool efp_fragments_only;

  int consecutive_backsteps_allowed;

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

  double linear_bend_threshold; // if bend exceeds this value, then also add linear bend complement

  // threshold for which entries in diagonalized redundant matrix are kept and inverted
  // while computing a generalized inverse of a matrix
  double redundant_eval_tol;

  // maximum number of allowed iterations in backtransformation to cartesian coordinates
  double bt_max_iter;

  double geom_maxiter;

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
  double IRC_step_size;
  bool keep_intcos; // don't delete intco.dat

  // for intcos with user-specified equilibrium values - this is the force constant
  double fixed_eq_val_force_constant;


};

}

#endif

