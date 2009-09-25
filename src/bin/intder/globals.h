/*! \file
    \ingroup INTDER
    \brief Enter brief description of file here 
*/
#ifndef _psi_bin_intder_globals_h_
#define _psi_bin_intder_globals_h_
#include "3dmatrix.h"

#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

namespace psi { namespace intder {

//EXTERN double*** Christoffel();
EXTERN double dot_prod(double *, double *);
EXTERN double* vect_prod(double *, double *);
EXTERN double** exp_matrix(double *);
EXTERN void mat1(double **, double *);
EXTERN void mat2(double **, double *);
EXTERN C3DMatrix* tripro();
EXTERN void fill3a(int , int , C3DMatrix *);
EXTERN void fill3b(int , int , C3DMatrix *);
EXTERN void fill4a(int , int , C4DMatrix *);
EXTERN void fill4a(int , int , C5DMatrix *);
EXTERN double** symm_matrix_invert(double **, int, int, int);
EXTERN double dot_x(double, int, double, int, int);

EXTERN FILE *fp_intco;
EXTERN FILE *fp_ider;

extern "C" {
EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
}

}} // namespace psi::intder

#define EVAL_TOL (1.0E-14)                /* tolerance for eigenvalues (used in sq_rsp() and irrep() )*/
#define REDUNDANT_EVAL_TOL (1.0E-10)
#define MAX_LINE (80)

#endif // header guard

