/*! \file
    \ingroup CPHF
    \brief Enter brief description of file here
*/

#ifndef _psi_bin_cphf_globals_h_
#define _psi_bin_cphf_globals_h_

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

namespace psi { namespace cphf {

EXTERN int *ioff;
#define IOFF_MAX 32641
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
EXTERN int X_only; /* only compute cphf coefficients, then stop */

/* setup.c */
EXTERN int natom, nmo, nso, nao, nirreps, ndocc, nuocc;
EXTERN int ntri, num_ai, num_pi, num_pq, ntei, noei, noei_ao;
EXTERN int *orbspi, *clsdpi, *openpi, *uoccpi;
EXTERN int *frdoccpi, *fruoccpi, *qtsorder;
EXTERN double *evals, *zvals, *ints, **scf, **usotao, **geom;
EXTERN int *first, *last, *ofirst, *olast, *vfirst, *vlast;
EXTERN int rottype, nnc;
EXTERN char **asymbol;
EXTERN double **dipder, **dipder_q;

EXTERN int print_lvl;
void zval_to_symbol(double zval, char *sym);

}} // namespace psi::cphf

#endif // header file
