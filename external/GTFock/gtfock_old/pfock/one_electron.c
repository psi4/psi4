#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ga.h>
#include <mkl.h>
#include <omp.h>
#include <mpi.h>
#include <mkl_scalapack.h>

#include "config.h"
#include "one_electron.h"


inline void matrix_block_write(double *matrix, int startrow,
                               int startcol, int ldm,
                               double *block, int nrows, int ncols)
{
    for (int k = 0; k < nrows; k++) {
        for (int l = 0; l < ncols; l++) {
            int i = startrow + k;
            int j = startcol + l;
            matrix[i * ldm + j] = block[k + nrows * l];
        }
    }
}


void compute_S(PFock_t pfock, BasisSet_t basis,
               int startshellrow, int endshellrow,
               int startshellcol, int endshellcol,
               int ldS, double *S)
{
    int nthreads = omp_get_max_threads();
    OED_t *oed = (OED_t *)malloc(sizeof(OED_t) * nthreads);
    assert(oed != NULL);
    for (int i = 0; i < nthreads; i++) {
        CInt_createOED(basis, &(oed[i]));
    }
    int start_row_id = pfock->f_startind[startshellrow];
    int start_col_id = pfock->f_startind[startshellcol];

    #pragma omp parallel
    {
        int tid = omp_get_thread_num ();
        #pragma omp for
        for (int A = startshellrow; A <= endshellrow; A++) {
            int row_id_1 = pfock->f_startind[A];
            int row_id_2 = pfock->f_startind[A + 1] - 1;
            int startrow = row_id_1 - start_row_id;
            int nrows = row_id_2 - row_id_1 + 1;
            for (int B = startshellcol; B <= endshellcol; B++) {
                int col_id_1 = pfock->f_startind[B];
                int col_id_2 = pfock->f_startind[B + 1] - 1;
                int startcol = col_id_1 - start_col_id;
                int ncols = col_id_2 - col_id_1 + 1;
                int nints;
                double *integrals;
                CInt_computePairOvl(basis, oed[tid], A, B, &integrals, &nints);
                if (nints != 0) {
                    matrix_block_write(S, startrow, startcol, ldS,
                                       integrals, nrows, ncols);
                }
            }
        }
    }

    for (int i = 0; i < nthreads; i++) {
        CInt_destroyOED(oed[i]);
    }
    free(oed);
}


void compute_H(PFock_t pfock, BasisSet_t basis,
               int startshellrow, int endshellrow,
               int startshellcol, int endshellcol,
               int ldH, double *H)
{
    int nthreads = omp_get_max_threads();
    OED_t *oed = (OED_t *)malloc(sizeof(OED_t) * nthreads);
    assert(oed != NULL);
    for (int i = 0; i < nthreads; i++) {
        CInt_createOED(basis, &(oed[i]));
    }
    
    int start_row_id = pfock->f_startind[startshellrow];
    int start_col_id = pfock->f_startind[startshellcol];
    #pragma omp parallel
    {
        int tid = omp_get_thread_num ();
        #pragma omp for
        for (int A = startshellrow; A <= endshellrow; A++) {
            int row_id_1 = pfock->f_startind[A];
            int row_id_2 = pfock->f_startind[A + 1] - 1;
            int startrow = row_id_1 - start_row_id;
            int nrows = row_id_2 - row_id_1 + 1;
            for (int B = startshellcol; B <= endshellcol; B++)
            {
                int col_id_1 = pfock->f_startind[B];
                int col_id_2 = pfock->f_startind[B + 1] - 1;
                int startcol = col_id_1 - start_col_id;
                int ncols = col_id_2 - col_id_1 + 1;
                int nints;
                double *integrals;
                CInt_computePairCoreH(basis, oed[tid],
                                      A, B, &integrals, &nints);
                if (nints != 0) {
                    matrix_block_write(H, startrow, startcol, ldH,
                                       integrals, nrows, ncols);
                }
            }
        }
    }

    for (int i = 0; i < nthreads; i++) {
        CInt_destroyOED(oed[i]);
    }
    free(oed);
}


void my_peig(int ga_A, int ga_B, int n, int nprow, int npcol, double *eval)
{
    int myrank;
    int ictxt;
    int myrow;
    int mycol;
    int nn;
    int mm;
    int izero = 0;
    int descA[9];
    int descZ[9];
    int info;
    int lo[2];
    int hi[2];
    int ld;
    int ione = 1;
#ifdef GA_NB
    ga_nbhdl_t nbnb;
#endif

    // init blacs
    int nb = MIN(n / nprow, n / npcol);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    Cblacs_pinfo(&nn, &mm);
    Cblacs_get(-1, 0, &ictxt);
    Cblacs_gridinit(&ictxt, "Row", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);

    // init matrices
    int nrows = numroc_(&n, &nb, &myrow, &izero, &nprow);
    int ncols = numroc_(&n, &nb, &mycol, &izero, &npcol);
    int itemp = nrows > 1 ? nrows : 1;
    descinit_(descA, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info);
    descinit_(descZ, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info);
    int blocksize = nrows * ncols;
    double *A = (double *)mkl_malloc(blocksize * sizeof (double), 64);
    assert (A != NULL);
    double *Z = (double *)mkl_malloc(blocksize * sizeof (double), 64);
    assert (Z != NULL);

    // distribute source matrix
    for (int i = 1; i <= nrows; i += nb) {
        lo[0] = indxl2g_(&i, &nb, &myrow, &izero, &nprow) - 1;
        hi[0] = lo[0] + nb - 1;
        hi[0] = hi[0] >= n ? n - 1 : hi[0];
        for (int j = 1; j <= ncols; j += nb) {
            lo[1] = indxl2g_(&j, &nb, &mycol, &izero, &npcol) - 1;
            hi[1] = lo[1] + nb - 1;
            hi[1] = hi[1] >= n ? n - 1 : hi[1];
            ld = ncols;
#ifdef GA_NB
            NGA_NbGet(ga_A, lo, hi, &(Z[(i - 1) * ncols + j - 1]), &ld, &nbnb);
#else
            NGA_Get(ga_A, lo, hi, &(Z[(i - 1) * ncols + j - 1]), &ld);
#endif
        }
        /* Jeff: Location of NGA_NbWait for flow-control. */
    }
#ifdef GA_NB
    /* Jeff: If one sees flow-control problems with too many
     *       outstanding NbGet operations, then move this call
     *       to the location noted above. */
    NGA_NbWait(&nbnb);
#endif
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            A[j * nrows + i] = Z[i * ncols + j];
        }
    }

    double t1 = MPI_Wtime();
    // inquire working space
    double *work = (double *)mkl_malloc(2 * sizeof (double), 64);
    assert (work != NULL);
    int lwork = -1;
#if 0
    pdsyev ("V", "U", &n, A, &ione, &ione, descA,
            eval, Z, &ione, &ione, descZ, work, &lwork, &info);
#else
    int liwork = -1;
    int *iwork = (int *)mkl_malloc(2 * sizeof (int), 64);
    assert(iwork != NULL);
    pdsyevd("V", "U", &n, A, &ione, &ione, descA,
            eval, Z, &ione, &ione, descZ,
            work, &lwork, iwork, &liwork, &info);    
#endif

    // compute eigenvalues and eigenvectors
    lwork = (int)work[0] * 2;
    mkl_free(work);
    work = (double *)mkl_malloc(lwork * sizeof (double), 64);
    assert(work != NULL);
#if 0
    pdsyev ("V", "U", &n, A, &ione, &ione, descA,
            eval, Z, &ione, &ione, descZ, work, &lwork, &info);
#else
    liwork = (int)iwork[0];
    mkl_free(iwork);
    iwork = (int *)mkl_malloc(liwork * sizeof (int), 64);
    assert(iwork != NULL);
    pdsyevd("V", "U", &n, A, &ione, &ione, descA,
            eval, Z, &ione, &ione, descZ,
            work, &lwork, iwork, &liwork, &info); 
#endif
    double t2 = MPI_Wtime();
    if (myrank == 0) {
        printf("  pdsyev_ takes %.3lf secs\n", t2 - t1);
    }

    // store desination matrix
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            A[i * ncols + j] = Z[j * nrows + i];
        }
    }
    for (int i = 1; i <= nrows; i += nb) {
        lo[0] = indxl2g_ (&i, &nb, &myrow, &izero, &nprow) - 1;
        hi[0] = lo[0] + nb - 1;
        hi[0] = hi[0] >= n ? n - 1 : hi[0];
        for (int j = 1; j <= ncols; j += nb) {
            lo[1] = indxl2g_ (&j, &nb, &mycol, &izero, &npcol) - 1;
            hi[1] = lo[1] + nb - 1;
            hi[1] = hi[1] >= n ? n - 1 : hi[1];
            ld = ncols;
#ifdef GA_NB
            NGA_NbPut(ga_B, lo, hi, &(A[(i - 1) * ncols + j - 1]), &ld, &nbnb);
#else
            NGA_Put(ga_B, lo, hi, &(A[(i - 1) * ncols + j - 1]), &ld);
#endif
        }
        /* Jeff: Location of NGA_NbWait for flow-control. */
    }
#ifdef GA_NB
    /* Jeff: If one sees flow-control problems with too many
     *       outstanding NbPut operations, then move this call
     *       to the location noted above. */
    NGA_NbWait(&nbnb);
#endif
    GA_Sync();

    mkl_free(A);
    mkl_free(Z);
    mkl_free(work);

    Cblacs_gridexit(ictxt);
}