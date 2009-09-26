/*##########################################################################*/
/*! \file
    \ingroup EXTREMA
  \brief Included header files, function declarations, and variables. */
/*##########################################################################*/

#ifndef _psi_bin_extrema_extrema_h_
#define _psi_bin_extrema_extrema_h_

#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else 
# define EXTERN
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <cctype>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <physconst.h>
#include "defines.h"
#include "math_tools.h"
#include "coord_base.h"
#include "simple.h"
#include "internals.h"
#include "zmat.h"
#include "deloc.h"

namespace psi { namespace extrema {
void punt(char*);
double **symm_matrix_invert(double**, int, int, int);

extern "C" {
EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
}

EXTERN int coord_type;
EXTERN int errcod, error;

/*this needs to be in C*/
extern "C" {
    char *gprgid();
}

/*inline function declarations*/
inline double *unit_vec(double* cart_arr, int atom1, int atom2 );
inline double vec_norm(double* cart_arr, int atom1, int atom2 );
inline double dot_pdt( double* vec1, double* vec2 );
inline double *cross_pdt( double* vec1, double* vec2 );
inline double norm(double* cart_arr, int atom1, int atom2 );
inline double compute_bond(double *car, int atm, int bnd);
inline double compute_angle(double *car, int atm, int bnd, int ang);
inline double compute_torsion(double *car, int atm, int bnd, int ang, int tor);

}} // namespace psi::extrema

#endif // header guard

