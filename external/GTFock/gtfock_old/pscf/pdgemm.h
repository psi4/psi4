#ifndef __PDGEMM_H__
#define __PDGEMM_H__


#include <mpi.h>


typedef struct _tmpbuf_t
{
    double *A;
    double *S;
    double *C;
    double *A_i;
    double *S_i;
    double *C_i;
} tmpbuf_t;


int pdgemm3D (int myrow, int mycol, int mygrd,
              MPI_Comm comm_row, MPI_Comm comm_col,
              MPI_Comm comm_grd, MPI_Comm comm_3D,
              int *nr, int *nc,
              int nrows, int ncols,
              double *D_, double *D2_, double *D3_,
              tmpbuf_t *tmpbuf);

void pdgemm3D_2 (int myrow, int mycol, int mygrd,
                 MPI_Comm comm_row, MPI_Comm comm_col,
                 MPI_Comm comm_grd, MPI_Comm comm_3D,
                 int *nr, int *nc, int nrows, int ncols,
                 double *A_block_, double *B_block_, double *C_block_,
                 tmpbuf_t *tmpbuf);

void allocate_tmpbuf (int nrows, int ncols, int *nr, int *nc,
                      tmpbuf_t * tmpbuf);

void dealloc_tmpbuf (tmpbuf_t * tmpbuf);


#endif /* __PDGEMM_H__ */
