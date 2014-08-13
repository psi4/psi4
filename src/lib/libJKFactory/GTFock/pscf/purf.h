#ifndef __PURF_H__
#define __PURF_H__


#include <mpi.h>

#include "CInt.h"


typedef struct _purf_t
{
    int runpurf;

    int nbf;
    int nobtls;
    // communicators
    MPI_Comm comm_purf;
    MPI_Comm comm_purf_row;
    MPI_Comm comm_purf_col;
    int np_purf;
    // np per row
    int nprow_purf;
    // np per col
    int npcol_purf;
    // nrows of processes in the same column
    int *nr_purf;
    // ncols of processes in the same row
    int *nc_purf;
    // positions
    int srow_purf;
    int scol_purf;
    int nrows_purf;
    int ncols_purf;
    // for matrix trace
    int tr_srow_purf;
    int tr_scol_purf;
    int tr_len_purf;
    // if contribute to trace
    int istr_purf;
    // block size for SUMMA dgemm
    int nb_purf;

    /**** local arrays ****/
    double *H_block;
    double *X_block;
    double *F_block;
    double *D_block;
    double *D2_block;
    double *FF_block;
    double *DD_block;
    // working space for purification
    double *work;
} purf_t;


purf_t *create_purf (BasisSet_t basis, int nprow_purf, int npcol_purf);

void destroy_purf (purf_t *purf);

int compute_purification (purf_t *purf, double *F_block, double *D_block);

void correction (purf_t *purf);


#endif /* __PURF_H__ */
