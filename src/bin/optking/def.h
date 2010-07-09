/*! \file
    \ingroup OPTKING
    \brief def.h the minimal header file
      includes macros, defines, hard limits and Intco_type declaration
        basic math functions
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>

#include <physconst.h>
#include <psifiles.h>

#include <string>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>

#ifndef _psi3_bin_optking_def_h_
#define _psi3_bin_optking_def_h_

namespace psi { //namespace optking {

// very basic prototypes for everyone
void cross_product(double *u,double *v,double *out);
void scalar_mult(double a, double *vect, int dim);
void scalar_div(double a, double *vect);
void punt(const char *message);
void swap(int *a, int *b);
void swap_tors(int *a, int *b, int *c, int *d);

void opt_ffile(FILE **fptr, const char *suffix, int code);
void opt_ffile_noexit(FILE **fptr, const char *suffix, int code);
void opt_mmult(double **AF, int ta, double **BF, int tb, double **CF, int tc,
  int nr, int nl, int nc, int add);
void opt_sq_rsp(int nm, int n, double **array, double *e_vals, int matz,
  double **e_vecs, double toler);

void exit_io(void);
void print_mat(double **a, int m, int n, FILE *out);
void eivout(double **a, double *b, int m, int n, FILE *out);
void dot_array(double *a, double *b, long int n, double *value);

void intro(void);
void intro(std::string header);
void free_info(int nsimples);

// memory allocation functions in mem.cc
void zero_array(double *a, long int size);
void free_array(double *a);
void free_int_array(int *a);
double **unit_matrix(long int m);

// enumerated types

enum Intco_type {STRE, BEND, TORS, OUT, LINB, FRAG};
enum Frag_switch {FRAG_A, FRAG_B};

// macros
#define SQR(A) ((A)*(A))

// defines
#define PRINT_TO_GEOM (113)
#define MAX_SALCS (500)
#define MAX_ZVARS (500)
#define MAX_LINELENGTH (133)
#define MAX_SALC_LENGTH (1000)

#define EVAL_TOL (1.0E-14)               /* tolerance for eigenvalues (used in sq_rsp() and irrep() ) */
#define REDUNDANT_EVAL_TOL (1.0E-10)
#define SPANNED_IRREP_TOL (0.05)         /* if character greater than this, irrep projected and kept */
#define LABEL_LENGTH (4) // for point group and irrep labels

#define NONLINEAR_DIST (1.0E-4) /* designed to exclude angle for CO2 if angle exceeds 179 */
#define MIN_DQ_STEP (1.0E-12)
#define MIN_CART_OUT (1.0E-12)
#define MIN_LIN_COS (1.0E-10)
#define FIX_NEAR_180 (150.0)

// optking running modes
#define MODE_DISP_NOSYMM   (10)
#define MODE_DISP_IRREP    (11)
#define MODE_DISP_LOAD     (12)
#define MODE_DISP_USER     (13)
#define MODE_LOAD_REF      (14)
#define MODE_OPT_STEP      (15)
#define MODE_FREQ_ENERGY   (16)
#define MODE_GRAD_ENERGY   (17)
#define MODE_FREQ_GRAD_NOSYMM (18)
#define MODE_FREQ_GRAD_IRREP  (19)
#define MODE_DISP_FREQ_GRAD_CART  (20)
#define MODE_FREQ_GRAD_CART  (21)
#define MODE_DISP_FREQ_ENERGY_CART  (22)
#define MODE_FREQ_ENERGY_CART  (23)
#define MODE_GRAD_SAVE        (24)
#define MODE_ENERGY_SAVE      (25)
#define MODE_RESET_PREFIX     (26)
#define MODE_DISP_NUM_PLUS    (27)
#define MODE_DELETE_BINARIES  (28)
#define MODE_TEST_BMAT        (29)
#define MODE_OPT_REPORT       (30)
#define MODE_DISP_INTERFRAGMENT      (31)
#define MODE_FREQ_GRAD_INTERFRAGMENT (32)
#define MODE_FCONST_INIT (33)

}//} /* namespace psi::optking */

#endif
