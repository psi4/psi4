#ifndef __PUTIL_H__
#define __PUTIL_H__

#include "PFock.h"
#include "CInt.h"
#include "purf.h"
#include <mpi.h>


#define MIN(a, b)    ((a) < (b) ? (a) : (b))
#define MAX(a, b)    ((a) > (b) ? (a) : (b))


extern"C"{ int indxg2p_ (int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);

int indxg2l_ (int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);

int indxl2g_ (int *indxlocal, int *nb, int *iproc, int *isrcproc, int *nprocs);

int numroc_ (int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);

void descinit_ (int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc,
                       int *ictxt, int *lld, int *info);

void Cblacs_pinfo (int *mypnum, int *nprocs);

void Cblacs_get (int context, int request, int* value);

int Cblacs_gridinit (int *context, char * order, int np_row, int np_col);

void Cblacs_gridinfo (int context, int * np_row, int *np_col, int *my_row, int *my_col);

void Cblacs_gridexit (int context);

void Cblacs_exit (int error_code);
}

// distributed dgemm using pipelined SUMMA with blocking
// C = alpha * A + B
// dim (A) = dim (B) = dim (C) = n x n
// bA, bB and bC are blocks of A, B and C respectively
// dim (bA) = dim (bB) = dim (bC) = nrows x ncols
// nb is the block size
// rr is the list of nrows of processors
// nc is the list of ncols of processors
void my_pdgemm (int n, int nb,
                double *A, double *B, double *C,
                int nrows, int ncols,
                int *nr, int *nc,
                MPI_Comm comm_row, MPI_Comm comm_col,
                double *work1, double *work2);

void my_peig (int ga_A, int ga_B, int n, int nprow, int npcol, double *eval);

void init_oedmat(BasisSet_t basis, PFock_t pfock, purf_t *purf, int nprow, int npcol);


#endif /* __PUTIL_H__ */
