
#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else
# define EXTERN
#endif

#include <cstdio>
#include <iostream>
#include <string>
#include <cmath>
#include <physconst.h>
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

extern "C" {
  EXTERN FILE *infile, *outfile;
  EXTERN char *psi_file_prefix;
}

namespace psi { namespace optking {

//#define OPT_AUX (13)
//#define SQR(A) ((A)*(A))
///*--- These are just internal type specifiers ---*/
//
//
//#define PRINT_TO_GEOM (113)
//#define PRINT_TO_30 (114)
///* Limits to hard-wired arrays */
//#define MAX_SALCS (500)
//#define MAX_ATOMS (300)
//#define MAX_ZVARS (500)
//#define MAX_SALC_LENGTH (1000)
//#define MAX(I,J) ((I>J) ? I : J)
//#define MIN(I,J) ((I>J) ? J : I)
//#define EVAL_TOL (1.0E-14)               /* tolerance for eigenvalues (used in sq_rsp() and irrep() ) */
//#define REDUNDANT_EVAL_TOL (1.0E-10)
//#define SPANNED_IRREP_TOL (0.05)         /* if character greater than this, irrep projected and kept */
//#define LABEL_LENGTH (4) // for point group and irrep labels
///* step size limits */
//#define STEP_LIMIT (0.1)     /* max step size if coord has small value */
//#define STEP_PERCENT (0.3)   /* if coord large valued, max percentage allowed for step */
//#define NONLINEAR_DIST (1.0E-4) /* designed to exclude angle for CO2 if angle exceeds 179 */
//#define MIN_DQ_STEP (1.0E-12)
//#define MIN_CART_OUT (1.0E-12)
//#define MIN_LIN_COS (1.0E-10)
///* optking running modes */
//#define MODE_DISP_NOSYMM   (10)
//#define MODE_DISP_IRREP    (11)
//#define MODE_DISP_LOAD     (12)
//#define MODE_DISP_USER     (13)
//#define MODE_LOAD_REF      (14)
//#define MODE_OPT_STEP      (15)
//#define MODE_FREQ_ENERGY   (16)
//#define MODE_GRAD_ENERGY   (17)
//#define MODE_FREQ_GRAD_NOSYMM (18)
//#define MODE_FREQ_GRAD_IRREP  (19)
//#define MODE_DISP_FREQ_GRAD_CART  (20)
//#define MODE_FREQ_GRAD_CART  (21)
//#define MODE_DISP_FREQ_ENERGY_CART  (22)
//#define MODE_FREQ_ENERGY_CART  (23)
//#define MODE_GRAD_SAVE        (24)
//#define MODE_ENERGY_SAVE      (25)
//#define MODE_RESET_PREFIX      (26)
//#define MODE_DISP_NUM_PLUS      (27)
//#define MODE_DELETE_BINARIES    (28)
//#define MODE_TEST_BMAT    (29)
//EXTERN FILE *fp_input, *fp_intco, *fp_fconst, *fp_opt_aux, *fp_11, *fp_fintco;
//EXTERN int *ops_in_class;
//EXTERN int nirreps, *irr;
//EXTERN int num_nonzero;  /* # of non-redundant di coordinates (evects of G with nonzero eigenvalues) */
//EXTERN char ptgrp[4];    /*molecular point group*/

//EXTERN void punt(char *message);
//EXTERN void open_PSIF(void);
//EXTERN void close_PSIF(void);
//EXTERN void exit_io(void);
//EXTERN void print_mat2(double **matrix, int rows, int cols, FILE *of);
//EXTERN void print_mat5(double **matrix, int rows, int cols, FILE *of);
//EXTERN void cross_product(double *u,double *v,double *out);
//EXTERN void scalar_mult(double a, double *vect, int dim);
//EXTERN void scalar_div(double a, double *vect);
//EXTERN int div_int(int big, int little);
//EXTERN double **symm_matrix_invert(double **A, int dim, int print_det, int redundant);
//EXTERN double energy_chkpt(void);
//EXTERN void dgeev_optking(int L, double **G, double *lambda, double **alpha);
//EXTERN void print_evects(double **evects, double *evals, int nrow, int ncol, FILE *out);
//EXTERN double **mass_mat(double *masses);
//EXTERN double **unit_mat(int dim);
//EXTERN void swap(int *a, int *b);
//EXTERN void swap_tors(int *a, int *b, int *c, int *d);
//EXTERN void zval_to_symbol(double zval, char *sym);
//EXTERN void symmetrize_geom(double *x);
//


struct OPTInfo {
  int natom;
  int print_lvl;
  double scale_connectivity;
//
  //int mode;
  //int disp_num;
  //int points;
  //int freq_irrep;
  //int points_freq_grad_ints;
  //int irrep;
  //int simples_present;
  //int salcs_present;
  //int constraints_present;
  //int nconstraints;
  //int *constraints;
  //int test_B;

//  int print_hessian;
//  int optimize;
  //int redundant;
  //int delocalize;
  //int do_only_displacements;
  //int zmat;
  //int zmat_simples;
  //int bfgs;
  //int bfgs_use_last;
  //int dertype;
  //int numerical_dertype;
  //int iteration;
  //int micro_iteration;
  //double conv; /* MAX force */
  //double ev_tol;
  //double scale_connectivity;
  //double disp_size;
  //int mix_types;
  //int natom;
  //int nallatom;
  //int *atom_dummy;
  //int *to_dummy;
  //int *to_nodummy;
  //int dummy_axis_1;
  //int dummy_axis_2;
  //char *wfn;
  //char *jobtype;
  //int energy_dat;
  //int grad_dat;
  //int external_energies; //ACS (11/07) Are we getting energy.dat from another program?
/* parameters involving fragment coordinates */
  //int frag_dist_rho;
  //int fix_interfragment;
  //int fix_intrafragment;
  //int freeze_intrafragment;

/* Back-transformation parameters */
  //int bt_max_iter;
  //double bt_dq_conv;
  //double bt_dx_conv;

/* Obscure limits in intco evaluation */
  //double cos_tors_near_1_tol;
  //double cos_tors_near_neg1_tol;
  //double sin_phi_denominator_tol;

  //int nfragment;
  //int *natom_per_fragment;
  //int *nallatom_per_fragment;
  //int *nref_per_fragment;
  //double ***fragment_coeff;
};

//struct SYMInfo {
  //char *symmetry;
  //int nirreps;
  //int **ct;
  //char **irrep_lbls;
  //char **clean_irrep_lbls;
  //char **op_lbls;
  //int **ict;
  //int **fict;
  //int **ict_ops;
  //int **ict_ops_sign;
//};
//
enum intco_type {stre_type, bend_type, tors_type, out_type, lin_bend_type, frag_type };

EXTERN struct OPTInfo optinfo;
//EXTERN struct SYMInfo syminfo;

}}

