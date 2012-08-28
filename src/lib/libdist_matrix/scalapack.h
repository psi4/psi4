#include <psiconfig.h>

#if defined(HAVE_SCALAPACK)

extern "C" {

extern void   Cblacs_pinfo( int* mypnum, int* nprocs);
extern void   Cblacs_get( int context, int request, int* value);
extern int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
extern void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
extern void   Cblacs_gridexit( int context);
extern void   Cblacs_exit( int error_code);

extern int    numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern void   descinit_( int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc,
				int *ictxt, int *lld, int *info);
extern double pdlamch_( int *ictxt , char *cmach);
extern double pdlange_( char *norm, int *m, int *n, double *A, int *ia, int *ja, int *desca, double *work);

extern void pdlacpy_( char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca,
				double *b, int *ib, int *jb, int *descb);
extern void pdgesv_( int *n, int *nrhs, double *A, int *ia, int *ja, int *desca, int* ipiv,
				double *B, int *ib, int *jb, int *descb, int *info);
extern void pdgemm_( char *TRANSA, char *TRANSB, int * M, int * N, int * K, double * ALPHA,
				double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB, int * DESCB,
				double * BETA, double * C, int * IC, int * JC, int * DESCC );
extern int  indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);

}

#endif

