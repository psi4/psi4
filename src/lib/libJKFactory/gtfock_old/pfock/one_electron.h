#ifndef __ONE_ELECTRON_H__
#define __ONE_ELECTRON_H__


#include "CInt.h"
#include "pfock.h"


void compute_S(PFock_t pfock, BasisSet_t basis,
               int startshellrow, int endshellrow,
               int startshellcol, int endshellcol,
               int ldS, double *S);

void compute_H(PFock_t pfock, BasisSet_t basis,
               int startshellrow, int endshellrow,
               int startshellcol, int endshellcol,
               int ldH, double *H);

void my_peig(int ga_A, int ga_B, int n, int nprow, int npcol, double *eval);

extern int indxg2p_(int *indxglob, int *nb, int *iproc, int *isrcproc,
                    int *nprocs);

extern int indxg2l_(int *indxglob, int *nb, int *iproc, int *isrcproc,
                    int *nprocs);

extern int indxl2g_(int *indxlocal, int *nb, int *iproc, int *isrcproc,
                    int *nprocs);

extern int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);

extern void descinit_(int *desc, int *m, int *n, int *mb, int *nb,
                      int *irsrc, int *icsrc, int *ictxt, int *lld,
                      int *info);

extern void Cblacs_pinfo(int *mypnum, int *nprocs);

extern void Cblacs_get(int context, int request, int *value);

extern int Cblacs_gridinit(int *context, char *order, int np_row,
                           int np_col);

extern void Cblacs_gridinfo(int context, int *np_row, int *np_col,
                            int *my_row, int *my_col);

extern void Cblacs_gridexit(int context);

extern void Cblacs_exit(int error_code);


#endif /* __ONE_ELECTRON_H__ */
