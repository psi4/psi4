#ifndef __CSCC_H
#define __CSCC_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

#define maxatom 384
#define maxnbfn 15 * maxatom
#define mxiter 30
#define maxnnbfn maxnbfn * (maxnbfn + 1) / 2
#define pi 3.141592653589793
#define tol 0.006
#define tol2e 0.000001

//#define MAX(a,b) GA_MAX(a,b)
//#define MIN(a,b) GA_MIN(a,b)

#define MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define MIN(a,b) (((a) <= (b)) ? (a) : (b))

extern double enrep, q[maxatom], ax[maxatom], ay[maxatom], az[maxatom],
  x[maxnbfn], y[maxnbfn], z[maxnbfn], expnt[maxnbfn], rnorm[maxnbfn];
extern long long int iky[maxnbfn], nocc, nbfn, nnbfn;
extern long long int icut1,icut2,icut3,icut4; 
extern int natom; //long long int --> long 

#define ichunk 40

extern double eigv[maxnbfn];
extern int g_counter, g_dens, g_fock, g_tfock, g_schwarz, g_work, g_ident, g_orbs;
extern int g_tmp, g_proc; //temporay global array for storage major transformation

#ifdef __cplusplus
}
#endif

#endif /* __CSCC_H */
