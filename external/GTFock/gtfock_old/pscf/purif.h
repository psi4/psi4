#ifndef __PURIF_H__
#define __PURIF_H__


#include <mpi.h>

#include "CInt.h"
#include "pfock.h"
#include "pdgemm.h"


#define MAX_DIIS 10
#define LDBMAT 128


typedef struct _purif_t
{
    int runpurif;
    int rundgemm;

    int nbf;
    int nobtls;
    // communicators
    MPI_Comm comm_purif;
    MPI_Comm comm_purif_row;
    MPI_Comm comm_purif_col;
    MPI_Comm comm_purif_grd;
    MPI_Comm comm_purif_plane;
    int np_purif;
    int nprow_purif;
    int npcol_purif;
    int npgrd_purif;
    // nrows of processes in the same column
    int *nr_purif;
    // ncols of processes in the same row
    int *nc_purif;
    // positions
    int srow_purif;
    int scol_purif;
    int nrows_purif;
    int ncols_purif;
    // for matrix trace
    int tr_srow_purif;
    int tr_scol_purif;
    int tr_len_purif;
    // if contribute to trace
    int istr_purif;
    // block size for SUMMA dgemm
    int nb_purif;

    /**** local arrays ****/
    int meshsize;
    double *H_block;
    double *X_block;
    double *S_block;
    double *F_block;
    double *D_block;
    double *D2_block;
    double *D3_block;
    tmpbuf_t tmpbuf;
    int ldx;
    // for diis
    int len_diis;
    double *diis_vecs;
    double *F_vecs;
    __declspec (align (64)) double b_mat[LDBMAT * LDBMAT]; // only on rank 0
    int bmax_id;
    double bmax; // only on rank 0

    // statistics
    double timedgemm;
    double timetr;
    double timepass;
} purif_t;


purif_t *create_purif(BasisSet_t basis, int nprow_purif,
                      int npcol_purif, int npgrd_purif);

void destroy_purif(purif_t * purif);

int compute_purification(purif_t * purif, double *F_block, double *D_block);

void compute_diis(PFock_t pfock, purif_t * purif,
                  double *D_block, double *F_block, int iter);

extern int numroc_(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);


#endif /* __PURIF_H__ */
