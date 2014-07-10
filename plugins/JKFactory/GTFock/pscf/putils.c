#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <mkl.h>
#include <string.h>
#include <mkl_scalapack.h>
#include <ga.h>
#include <math.h>

#include "putils.h"


// ring broadcast (pipelined)
static void ring_bcast (double *buf, int count, MPI_Datatype type, int root, MPI_Comm comm)
{
    int me;
    int np;
    MPI_Status status;

    MPI_Comm_rank (comm, &me);
    MPI_Comm_size (comm, &np);
    if (me != root)
    {
        MPI_Recv (buf, count, type, (me - 1 + np) % np, MPI_ANY_TAG, comm,
                  &status);
    }
    if ((me + 1) % np != root)
    {
        MPI_Send (buf, count, type, (me + 1) % np, 0, comm);
    }
}


void my_pdgemm (int n, int nb,
                double *A, double *B, double *C,
                int nrows, int ncols,
                int *nr, int *nc,
                MPI_Comm comm_row, MPI_Comm comm_col,
                double *work1, double *work2)
{
    int myrank_row;
    int myrank_col;
    int myrow;
    int mycol;
    int kk;
    int iwrk;
    int icurrow;
    int icurcol;
    int ii;
    int jj;

    // get row and col communicator
    MPI_Comm_rank (comm_row, &myrank_row);
    MPI_Comm_rank (comm_col, &myrank_col);
    myrow = myrank_col;
    mycol = myrank_row;

    // zeros C
    memset (C, 0, sizeof(double) * nrows * ncols);

    icurrow = 0;
    icurcol = 0;
    ii = jj = 0;
    iwrk = 0;
    // main loop
    for (kk = 0; kk < n; kk += iwrk)
    {
        iwrk = MIN (nb, nr[icurrow] - ii);
        iwrk = MIN (iwrk, nc[icurcol] - jj);
        if (mycol == icurcol)
        {
            dlacpy_ ("General", &iwrk, &nrows, &A[jj], &ncols, work1, &iwrk);
        }
        if (myrow == icurrow)
        {
            dlacpy_ ("General", &ncols, &iwrk, &B[ii * ncols], &ncols, work2,
                     &ncols);
        }
        ring_bcast (work1, nrows * iwrk, MPI_DOUBLE, icurcol, comm_row);
        ring_bcast (work2, ncols * iwrk, MPI_DOUBLE, icurrow, comm_col);
        cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, nrows, ncols,
                     iwrk, 1.0, work1, iwrk, work2, ncols, 1.0, C, ncols);
        ii += iwrk;
        jj += iwrk;
        if (jj >= nc[icurcol])
        {
            icurcol++;
            jj = 0;
        }
        if (ii >= nr[icurrow])
        {
            icurrow++;
            ii = 0;
        }
    }
}


void my_peig (int ga_A, int ga_B, int n, int nprow, int npcol, double *eval)
{
    int myrank;
    int ictxt;
    int myrow;
    int mycol;
    int nn;
    int mm;
    int nrows;
    int ncols;
    int nb;
    int izero = 0;
    int descA[9];
    int descZ[9];
    int info;
    int itemp;
    int i;
    int blocksize;
    double *A;
    double *Z;
    double *work;
    int lo[2];
    int hi[2];
    int ld;
    int lwork;
    int ione = 1;
    int j;

    // init blacs
    nb = MIN(n/nprow, n/npcol);    
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    Cblacs_pinfo (&nn, &mm) ;
    Cblacs_get (-1, 0, &ictxt);
    Cblacs_gridinit (&ictxt, "Row", nprow, npcol);
    Cblacs_gridinfo (ictxt, &nprow, &npcol, &myrow, &mycol);

    // init matrices
    nrows = numroc_ (&n, &nb, &myrow, &izero, &nprow);
    ncols = numroc_ (&n, &nb, &mycol, &izero, &npcol);
    itemp = nrows > 1 ? nrows : 1;
    descinit_ (descA, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info);
    descinit_ (descZ, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info);    
    blocksize = nrows * ncols;
    A = (double *)malloc (blocksize * sizeof(double));
    assert (A != NULL);
    Z = (double *)malloc (blocksize * sizeof(double));
    assert (Z != NULL);

    // distribute source matrix
    for (i = 1; i <= nrows; i+=nb)
    {        
        lo[0] = indxl2g_ (&i, &nb, &myrow, &izero, &nprow) - 1;
        hi[0] = lo[0] + nb - 1;
        hi[0] = hi[0] >= n ? n - 1 : hi[0]; 
        for (j = 1; j <= ncols; j+=nb)
        {
            lo[1] = indxl2g_ (&j, &nb, &mycol, &izero, &npcol) - 1;
            hi[1] = lo[1] + nb - 1;
            hi[1] = hi[1] >= n ? n - 1 : hi[1]; 
            ld = ncols;
            NGA_Get (ga_A, lo, hi, &(Z[(i - 1) * ncols + j - 1]), &ld);
        }
    }
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            A[j * nrows + i] = Z[i * ncols + j];
        }
    }

    // inquire working space
    work = (double *)malloc (2 * sizeof(double));
    assert (work != NULL);
    lwork = -1;    
    pdsyev_ ("V", "U", &n, A, &ione, &ione, descA,
             eval, Z, &ione, &ione, descZ, work, &lwork, &info);

    // compute eigenvalues and eigenvectors
    lwork = (int)work[0];
    free (work);          
    work = (double *)malloc (lwork * sizeof(double));
    assert (work != NULL); 
    pdsyev_ ("V", "U", &n, A, &ione, &ione, descA,
             eval, Z, &ione, &ione, descZ, work, &lwork, &info);

    // store desination matrix
    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            A[i * ncols + j] = Z[j * nrows + i];
        }
    }
    for (i = 1; i <= nrows; i+=nb)
    {        
        lo[0] = indxl2g_ (&i, &nb, &myrow, &izero, &nprow) - 1;
        hi[0] = lo[0] + nb - 1;
        hi[0] = hi[0] >= n ? n - 1 : hi[0]; 
        for (j = 1; j <= ncols; j+=nb)
        {
            lo[1] = indxl2g_ (&j, &nb, &mycol, &izero, &npcol) - 1;
            hi[1] = lo[1] + nb - 1;
            hi[1] = hi[1] >= n ? n - 1 : hi[1]; 
            ld = ncols;
            NGA_Put (ga_B, lo, hi, &(A[(i - 1) * ncols + j - 1]), &ld);
        }
    }
    NGA_Sync ();

    free (A);
    free (Z);
    free (work);
    
    Cblacs_gridexit (ictxt);
}

void init_oedmat (BasisSet_t basis, PFock_t pfock, purf_t *purf, int nprow, int npcol)
{
    int srow_purf;
    int scol_purf;
    int nrows_purf;
    int ncols_purf;
    int erow_purf;
    int ecol_purf;
    int nfuncs_row;
    int nfuncs_col;
    int lo[2];
    int hi[2];
    int ld;
    double *eval;
    double lambda;
    int nbf;
    int myrank;
    double *blocktmp;
    double *blockS;
    int ga_S;
    int ga_X;
    int ga_tmp;
    int j;
    int i;
    CoreH_t hmat;
    Ovl_t smat;
    double t1;
    double t2;
    
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    if (myrank == 0)
    {
        printf ("Computing one electron matrices ...\n");
    }   
    srow_purf = purf->srow_purf;
    scol_purf = purf->scol_purf;
    nrows_purf = purf->nrows_purf;
    ncols_purf = purf->ncols_purf;
    erow_purf = srow_purf + nrows_purf - 1;
    ecol_purf = scol_purf + ncols_purf - 1;
    
    /* Compute X */
    if (myrank == 0)
    {
        printf ("  constructing X\n");
    }
    t1 = MPI_Wtime ();
    // compute S
    PFock_createOvlMat (pfock, basis, &smat);
    PFock_getOvlMatGAHandle (smat, &ga_S);
    PFock_getMatGAHandle (pfock, 0, PFOCK_MAT_TYPE_D, &ga_tmp);
#ifdef __PFOCK_SCF__
    PFock_getMatGAHandle (pfock, 0, PFOCK_MAT_TYPE_F, &ga_X);  
#else
    PFock_getMatGAHandle (pfock, 0, PFOCK_MAT_TYPE_J, &ga_X);
#endif
    // eigensolve         
    nbf = CInt_getNumFuncs (basis);
    eval = (double *)malloc (nbf * sizeof(double));
    assert (eval != NULL);
    
    //GA_Diag_std (ga_S, ga_tmp, eval);

    my_peig (ga_S, ga_tmp, nbf, nprow, npcol, eval);

    NGA_Distribution (ga_tmp, myrank, lo, hi);
    nfuncs_row = hi[0] - lo[0] + 1;
    nfuncs_col = hi[1] - lo[1] + 1;
    NGA_Access (ga_tmp, lo, hi, &blocktmp, &ld);
    NGA_Access (ga_S, lo, hi, &blockS, &ld);    
    for (j = 0; j < nfuncs_col; j++)
    {
        lambda = 1.0 / sqrt (eval[j + lo[1]]);                
        for (i = 0; i < nfuncs_row; i++)
        {
            blockS[i * nfuncs_col + j] = blocktmp[i * nfuncs_col + j] * lambda;
        }
    }
    NGA_Release (ga_tmp, lo, hi);    
    NGA_Release_update (ga_S, lo, hi); 
    // get X
    GA_Dgemm ('N', 'T', nbf, nbf, nbf, 1.0,
              ga_S, ga_tmp, 0.0, ga_X);
    if (purf->runpurf == 1)
    {
        lo[0] = srow_purf;
        hi[0] = erow_purf;
        lo[1] = scol_purf;
        hi[1] = ecol_purf;
        ld = ncols_purf;
        NGA_Get (ga_X, lo, hi, purf->X_block, &ld);
    }
    PFock_destroyOvlMat (smat);
    free (eval);

    /* Compute H */
    if (myrank == 0)
    {
        printf ("  computing H\n");
    }
    PFock_createCoreHMat (pfock, basis, &hmat);
    if (purf->runpurf == 1)
    {
        PFock_getCoreHMat (hmat, srow_purf, erow_purf,
                           scol_purf, ecol_purf,
                           purf->H_block, ncols_purf);
    }
    PFock_destroyCoreHMat (hmat);
     
    t2 = MPI_Wtime ();
    if (myrank == 0)
    {
        printf ("  takes %.3lf secs\n", t2 - t1);
        printf ("  Done\n");
    }
}
