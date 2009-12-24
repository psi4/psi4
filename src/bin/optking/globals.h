/*! \file
    \ingroup OPTKING
    \brief globals.h : includes minimal header plus global variables
*/
#include "def.h"

#ifndef _psi3_bin_optking_globals_h_
#define _psi3_bin_optking_globals_h_

namespace psi { namespace optking {

/* Global variables */
#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

extern "C" {
  EXTERN FILE *infile, *outfile;
  EXTERN char *psi_file_prefix;
}

EXTERN FILE *fp_input, *fp_intco, *fp_fconst, *fp_opt_aux, *fp_11, *fp_fintco, *fp_energy_dat;
EXTERN int *ops_in_class;
EXTERN int nirreps, *irr;
EXTERN int num_nonzero;  /* # of non-redundant di coordinates (evects of G with nonzero eigenvalues) */
EXTERN char ptgrp[4];    /*molecular point group*/

EXTERN struct OPTInfo {
  int mode;
  int disp_num;
  int points;
  int freq_irrep;
  int points_freq_grad_ints;
  int irrep;
  int simples_present;
  int salcs_present;
  int constraints_present;
  int nconstraints;
  int *constraints;
  int test_B;

/* print options */
  int print_simples;
  int print_params;
  int print_delocalize;
  int print_symmetry;
  int print_hessian;
  int print_cartesians;
  int print_fconst;
  int print_debug_backtransformation;

/* optimization parameters */
  bool ts;  // search for a transition state? 
  bool selected_fc; // choose coordinates for which to get force constants
  int balked_last_time; // boolean to indicate we've recomputed fc's
  int cartesian;
  int optimize;
  int redundant;
  int delocalize;
  int do_only_displacements;
  int zmat;
  int zmat_simples;

  enum { NONE, BFGS, MS, POWELL, BOFILL} H_update;
  enum { FISCHER, SCHLEGEL} empirical_H;
  enum { EMPIRICAL, KEEP } nonselected_fc;
  enum { LOOSE, NORMAL, TIGHT, VERY_TIGHT, BAKER, QCHEM, G03_NORMAL, G03_TIGHT, G03_VERY_TIGHT, GENERAL } opt_conv;

  int H_update_use_last;
  bool rfo; // whether to use rfo step
  bool rfo_follow_root; // whether to follow an initial root
  int rfo_root; // which root to follow
  double step_energy_limit; // fraction error in energy prediction to tolerate
  double step_energy_limit_back; // fraction error to tolerate in guess step backward
  double step_limit; // max change in au or radians
  int dertype;
  int numerical_dertype;
  int iteration;
  int micro_iteration;
  double conv_max_force;  /* MAX force convergence */
  double conv_max_DE; /* MAX delta-E convergence */
  double conv_max_disp; /* MAX displacement convergence */
  double conv_rms_force; /* RMS force convergence */
  double conv_rms_disp; /* RMS displacement convergence */
  double ev_tol;
  double scale_connectivity;
  double disp_size;
  double step_limit_cart;
  int mix_types;
  int natom;
  int nallatom;
  int *atom_dummy;
  int *to_dummy;
  int *to_nodummy;
  int dummy_axis_1;
  char *wfn;
  char *jobtype;
  int energy_dat;
  int grad_dat;
  int external_energies; //ACS (11/07) Are we getting energy.dat from another program?
/* parameters involving fragment coordinates */
  int frag_dist_rho;
  int fix_interfragment;
  int fix_intrafragment;
  int analytic_interfragment;
  int max_consecutive_line_searches;
  double line_search_min;

/* Back-transformation parameters */
  int bt_max_iter;
  double bt_dq_conv;
  double bt_dx_conv;

/* Obscure limits in intco evaluation */
  double cos_tors_near_1_tol;
  double cos_tors_near_neg1_tol;
  double sin_phi_denominator_tol;

  int nfragment;
  int *natom_per_fragment;
  int *nallatom_per_fragment;
  int *nref_per_fragment;
  double ***fragment_coeff;
} optinfo;

EXTERN struct SYMInfo {
  char *symmetry;
  int nirreps;
  int **ct;
  char **irrep_lbls;
  char **clean_irrep_lbls;
  const char **op_lbls;
  int **ict;
  int **fict;
  int **ict_ops;
  int **ict_ops_sign;
} syminfo;

}} /* namespace psi::optking */

#endif
