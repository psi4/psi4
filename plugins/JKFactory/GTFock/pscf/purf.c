#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include <mkl.h>
#include <mkl_trans.h>
#include<iostream>
#include "purf.h"
#include "putils.h"

using namespace std;

static void config_purf (purf_t *purf)
{
    int nbf;
    int nrows;
    int ncols;
    int nb;
    int nprow;
    int npcol;
    int myrow;
    int mycol;
    MPI_Comm comm_col;
    MPI_Comm comm_row;
    int *nr;
    int *nc;
    int startrow;
    int endrow;
    int startcol;
    int endcol;
    int start;
    int meshsize;
    double memsize;
    int i;
    int izero = 0;

    //Get the number of basis functions
    nbf = purf->nbf;
    //Number of processors per row
    nprow = purf->nprow_purf;
    //Number of processors per column
    npcol = purf->npcol_purf;
    //Here they are figuring out the size of of the blocks and assign it to the smaller of the two
    purf->nb_purf = nb = MIN(nbf/nprow, nbf/npcol);
    cout<<"nb is:"<<nb<<endl;
    //These are communicators for the rows and columns of the matrix
    comm_col = purf->comm_purf_col;
    comm_row = purf->comm_purf_row;
    MPI_Comm_rank (comm_row, &mycol);
    MPI_Comm_rank (comm_col, &myrow);

    //Numroc returns the number of rows or columns owned by the current process
    cout<<"nbf,nb,myrow,mycol,izero,nprow,npcol "<<nbf<<" "<<nb<<" "<<myrow<<" "<<mycol<<" "<<izero<<" "<<nprow<<" "<<npcol<<endl;
    nrows = numroc_ (&nbf, &nb, &myrow, &izero, &nprow);
    ncols = numroc_ (&nbf, &nb, &mycol, &izero, &npcol);
    cout<<"nrows, ncols "<<nrows<<" "<<ncols<<endl;
    purf->nrows_purf = nrows;
    purf->ncols_purf = ncols;

    // positions of partitions
    purf->nr_purf = (int *)malloc (sizeof(int) * purf->nprow_purf);
    purf->nc_purf = (int *)malloc (sizeof(int) * purf->npcol_purf);
    assert (purf->nr_purf != NULL);
    assert (purf->nc_purf != NULL);   
    nr = purf->nr_purf;
    nc = purf->nc_purf;

    // get nr and sr
    MPI_Allgather (&nrows, 1, MPI_INT, nr, 1, MPI_INT, comm_col);
    startrow = 0;
    for (i = 0; i < myrow; i++)
    {
        startrow += nr[i];
    }
    endrow = startrow + nrows - 1;
    purf->srow_purf = startrow;


    // get nc and sc
    MPI_Allgather (&ncols, 1, MPI_INT, nc, 1, MPI_INT, comm_row);
    startcol = 0;
    for (i = 0; i < mycol; i++)
    {
        startcol += nc[i];
    }
    endcol = startcol + ncols - 1;
    purf->scol_purf = startcol;
    std::cout<<"Start "<<startrow<<","<<startcol<<" end: "<<endrow<<","<<endcol<<std::endl;
    // for matrix trace
    start = MAX (startcol, startrow);
    purf->tr_len_purf = MIN (endcol, endrow) - start + 1;
    purf->tr_scol_purf = start - startcol;
    purf->tr_srow_purf = start - startrow;
    purf->istr_purf = (purf->tr_len_purf > 0);

    // create local arrays
    meshsize = nrows * ncols;
    purf->X_block = (double *)malloc (meshsize * sizeof(double));        
    assert (purf->X_block != NULL);
    purf->H_block = (double *)malloc (meshsize * sizeof(double));
    assert (purf->H_block != NULL); 
    purf->F_block = (double *)malloc (meshsize * sizeof(double));
    assert (purf->F_block != NULL);
    purf->D_block = (double *)malloc (meshsize * sizeof(double));
    assert (purf->D_block != NULL);
    purf->FF_block = (double *)malloc (meshsize * sizeof(double));
    assert (purf->FF_block != NULL);
    purf->DD_block = (double *)malloc (meshsize * sizeof(double));
    assert (purf->DD_block != NULL);
    purf->D2_block = (double *)malloc (meshsize * sizeof(double));
    assert (purf->D2_block != NULL);
    // working space for purification
    purf->work = (double *)malloc (4 * nb * MAX (nrows, ncols) * sizeof(double));
    assert (purf->work != NULL);

    memsize = (7.0 * meshsize + 4.0 * nb * MAX (nrows, ncols)) * sizeof(double); 
    if (myrow == 0 && mycol == 0)
    {
        printf ("  using %.3lf MB\n", memsize/1024.0/1024.0);
    }
}


purf_t *create_purf (BasisSet_t basis, int nprow_purf, int npcol_purf)
{
    purf_t *purf;
    int myrank;
    int flag_purf;
    int rowid;
    int colid;
    double t1;
    double t2;
    
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);   
    if (myrank == 0)printf ("Initializing purification ...\n");
    t1 = MPI_Wtime ();
    
    // create purf object
    purf = (purf_t *)malloc (sizeof(purf_t));
    assert (purf != NULL);
    purf->nbf = CInt_getNumFuncs (basis);
    purf->nobtls = CInt_getNumOccOrb (basis);
    purf->nprow_purf = nprow_purf;
    purf->npcol_purf = npcol_purf;
    purf->np_purf = nprow_purf * npcol_purf;

    //If I am node 0 to number of processors -1 I am running the purification
    if (myrank < purf->np_purf)purf->runpurf = 1;
    else purf->runpurf = 0;

    // initialize communicators    
    flag_purf = (myrank < purf->np_purf);
    //Take the processors in the World and split them into those that are running
    //the purification (flag_purf=true), and those that are not (flag_purf=false)
    //ordered by their original rank.  These new ranks are their ranks in comm_purf
    MPI_Comm_split (MPI_COMM_WORLD, flag_purf, myrank, &(purf->comm_purf));
    if (purf->runpurf == 1)//If I am part of the processors running the purification...
    {
        //Now we have to divvy up the matrix
    	//The first operation here will give us a number between 0 and np_purf/npcol_purf
    	//Note total=np_row*np_col, so this is actually telling us which one of the np_row processors we are
    	rowid = myrank / purf->npcol_purf;
    	//This one gives a number between 0 and np_col-1, so which one of the np_col processors we are
        colid = myrank % purf->npcol_purf;
        //Now take all processors that running a given row, ordered by column and assign them to a comm
        MPI_Comm_split (purf->comm_purf, rowid, colid, &(purf->comm_purf_row));
        //Take all processors running a given column, ordered by row and assign them to a comm
        MPI_Comm_split (purf->comm_purf, colid, rowid, &(purf->comm_purf_col));

        config_purf (purf);
    }

    t2 = MPI_Wtime ();
    if (myrank == 0)
    {
        printf ("  takes %.3lf secs\n", t2 - t1);
        printf ("  Done\n");
    }
    return purf;
}


void destroy_purf (purf_t *purf)
{
    if (purf->runpurf == 1)
    {
        MPI_Comm_free (&(purf->comm_purf));
        MPI_Comm_free (&(purf->comm_purf_row));
        MPI_Comm_free (&(purf->comm_purf_col));
        free (purf->nr_purf);
        free (purf->nc_purf);
        free (purf->H_block);
        free (purf->X_block);
        free (purf->F_block);
        free (purf->D_block);
        free (purf->FF_block);
        free (purf->DD_block);
        free (purf->D2_block);
    }
    free (purf);
}


int compute_purification (purf_t *purf, double *F_block, double *D_block)
{
    // partition infor
    int nrows;
    int ncols;
    int startrow;
    int startcol;
    int *nr;
    int *nc;
    int myrow;
    int mycol;
    MPI_Comm comm_row;
    MPI_Comm comm_col;
    MPI_Comm comm_purf;
    // for matrix trace
    int starttrrow;
    int starttrcol;
    int lentr;
    // local arrays
    double *X_block;
    double *D2_block;    
    double *D3_block;
    // working arrays
    double *workm;
    double *work;
    double *work1;
    double *work2;
    double *work3;
    double *work4;
    int it;
    int i;
    int j;
    // local variables
    double mu_bar;
    double lambda;
    double hmax;
    double hmin;
    double tr;
    double tr2;
    int myrank;
    double tmp;
    int nb;
    int nbf;
    int nobtls;
    
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

    if (purf->runpurf == 1)
    {
        // initialization
        nb = purf->nb_purf;
        X_block = purf->X_block;
        D2_block = purf->D2_block;
        D3_block = purf->FF_block;
        nrows = purf->nrows_purf;
        ncols = purf->ncols_purf;
        startrow = purf->srow_purf;
        startcol = purf->scol_purf;
        starttrrow = purf->tr_srow_purf;
        starttrcol = purf->tr_scol_purf;
        lentr = purf->tr_len_purf;
        nr = purf->nr_purf;
        nc = purf->nc_purf;
        comm_row = purf->comm_purf_row;
        comm_col = purf->comm_purf_col;
        comm_purf = purf->comm_purf;
        nbf = purf->nbf;
        nobtls = purf->nobtls;

        MPI_Comm_rank (comm_purf, &myrank);
        MPI_Comm_rank (comm_row, &mycol);
        MPI_Comm_rank (comm_col, &myrow);
        
        workm = D2_block; 
        work = purf->work;    
        work1 = work;
        work2 = work + nb * nrows;
        work3 = work + ncols;
        work4 = work + 2 * ncols;
        
        // F = X'*F*X;
        my_pdgemm (nbf, nb, X_block, F_block, workm, nrows, ncols,
                   nr, nc, comm_row, comm_col, work1, work2);
        my_pdgemm (nbf, nb, workm, X_block, D_block, nrows, ncols,
                   nr, nc, comm_row, comm_col, work1, work2);
        
        // compute eigenvalue estimates using Gershgorin (F is symmetric)
        // offdF = sum(abs(F))' - abs(diag(F));
        // diagF = diag(F);
        // hmin = min(diagF - offdF);
        // hmax = max(diagF + offdF);      
        for (i = 0; i < ncols; i++)
        {
            work1[i] = 0.0;
            work3[i] = 0.0;
            for (j = 0; j < nrows; j++)
            {
                work1[i] += fabs (D_block[i + j * ncols]);
                if (j + startrow == i + startcol)
                {
                    work3[i] = D_block[i + j * ncols];
                }
            }
            work1[i] = work1[i] - fabs (work3[i]);
            tmp = work1[i] + work3[i];
            work1[i] = work3[i] - work1[i];
            work3[i] = tmp;
        }
        MPI_Reduce (work, work4, 2 * ncols, MPI_DOUBLE, MPI_SUM, 0, comm_col);

        if (myrow == 0)
        {
            work1[0] = DBL_MAX;
            work1[1] = -DBL_MAX;
            for (i = 0; i < ncols; i++)
            {
                work1[0] = work4[i] > work1[0] ? work1[0] : work4[i];
                work1[1] = work4[i + ncols] < work1[1] ? work1[1] : work4[i + ncols];
            }
            MPI_Reduce (&work1[0], &hmin, 1, MPI_DOUBLE, MPI_MIN, 0,
                        comm_row);
            MPI_Reduce (&work1[1], &hmax, 1, MPI_DOUBLE, MPI_MAX, 0,
                        comm_row);
        }
        MPI_Bcast (&hmin, 1, MPI_DOUBLE, 0, comm_purf);
        MPI_Bcast (&hmax, 1, MPI_DOUBLE, 0, comm_purf);
        // define constants, dependent on F
        // in the following:
        // 5 = no of occupied orbitals
        // 7 = no of spatial basis function (each corresponds to 2 electrons for RHF)
        // mu_bar = trace_dense_matrix(F)/7;
        work[0] = 0.0;
        for (i = 0; i < lentr; i++)
        {
            work[0] += D_block[(i + starttrrow) * ncols + i + starttrcol];
        }
        MPI_Reduce (&work[0], &tr, 1, MPI_DOUBLE, MPI_SUM, 0, comm_purf);

        if (myrank == 0)
        {
            mu_bar = tr / (double) nbf;
            // lambda = min([ 5/(hmax - mu_bar), (7-5)/(mu_bar - hmin) ]);
            lambda = MIN ((double) nobtls / (hmax - mu_bar),
                              (double) (nbf - nobtls) / (mu_bar - hmin));
        }
        MPI_Bcast (&lambda, 1, MPI_DOUBLE, 0, comm_purf);
        MPI_Bcast (&mu_bar, 1, MPI_DOUBLE, 0, comm_purf);

        // initial "guess" for density matrix
        // D = (lambda*mu_bar/7 + 5/7)*eye(7) - (lambda/7)*D;
        cblas_dscal (nrows * ncols, -lambda / (double) nbf,D_block, 1);
        
        for (i = 0; i < lentr; i++)
        {
            D_block[(i + starttrrow) * ncols + i + starttrcol] +=
                lambda * mu_bar / (double) nbf + (double) nobtls / nbf;
        }

        // McWeeny purification
        // convergence appears slow at first, before accelerating at end
        for (it = 0; it < 100; it++)
        {
            double errnorm;
            double c;
            
            // D2 = D*D;
            // D3 = D2*D;
            my_pdgemm (nbf, nb, D_block, D_block, D2_block, nrows, ncols,
                       nr, nc, comm_row, comm_col, work1, work2);

            my_pdgemm (nbf, nb, D2_block, D_block, D3_block, nrows, ncols,
                       nr, nc, comm_row, comm_col, work1, work2);
            
            // stopping criterion
            // errnorm = norm(D-D2, 'fro');
            work[0] = 0.0;
            for (i = 0; i < nrows * ncols; i++)
            {
                work[0] +=
                    (D_block[i] - D2_block[i]) * (D_block[i] - D2_block[i]);
            }
            MPI_Reduce (&work[0], &errnorm, 1, MPI_DOUBLE, MPI_SUM, 0,
                        comm_purf);
            if (myrank == 0)
            {
                errnorm = sqrt (errnorm);
            }
            MPI_Bcast (&errnorm, 1, MPI_DOUBLE, 0, comm_purf);
            if (errnorm < 1e-11)
            {
                break;
            }

            // a cheaper stopping criterion may be to check the trace of D*D
            // and stop when it is close to no. occupied orbitals (5 in this case)
            // fprintf('trace D*D //f\n', trace(D*D);
            // might be possible to "lag" the computation of c by one iteration
            // so that the global communication for traces can be overlapped.
            // Note: c appears to converge to 0.5           
            // c = trace(D2-D3) / trace(D-D2);
            work[0] = 0.0;
            work[1] = 0.0;
            for (i = 0; i < lentr; i++)
            {
                work[0] +=
                    D2_block[(i + starttrrow) * ncols + i + starttrcol] -
                    D3_block[(i + starttrrow) * ncols + i + starttrcol];
                work[1] +=
                    D_block[(i + starttrrow) * ncols + i + starttrcol] -
                    D2_block[(i + starttrrow) * ncols + i + starttrcol];
            }
            MPI_Reduce (&work[0], &tr, 1, MPI_DOUBLE, MPI_SUM, 0,
                        comm_purf);
            MPI_Reduce (&work[1], &tr2, 1, MPI_DOUBLE, MPI_SUM, 0,
                        comm_purf);
            if (myrank == 0)
            {
                c = tr / tr2;
            }
            MPI_Bcast (&c, 1, MPI_DOUBLE, 0, comm_purf);

            if (c < 0.5)
            {
                for (i = 0; i < nrows * ncols; i++)
                {
                    // D = ((1-2*c)*D + (1+c)*D2 - D3) / (1-c);
                    D_block[i] =
                        ((1.0 - 2.0 * c) * D_block[i] + 
                        (1.0 + c) * D2_block[i] -
                        D3_block[i]) / (1.0 - c);
                }
            }
            else
            {
                for (i = 0; i < nrows * ncols; i++)
                {
                    // D = ((1+c)*D2 - D3) / c;
                    D_block[i] = ((1.0 + c) * D2_block[i] - D3_block[i]) / c;
                }
            }
        }

        // change basis
        //D = X*D*X';
        my_pdgemm (nbf, nb, X_block, D_block, workm, nrows, ncols,
                   nr, nc, comm_row, comm_col, work1, work2);
        my_pdgemm (nbf, nb, workm, X_block, D_block, nrows, ncols,
                   nr, nc, comm_row, comm_col, work1, work2);
    }
    
    MPI_Barrier (MPI_COMM_WORLD);

    return it;
}


void correction (purf_t *purf)
{
    // partition infor
    int nrows;
    int ncols;
    int *nr;
    int *nc;
    int myrow;
    int mycol;
    MPI_Comm comm_row;
    MPI_Comm comm_col;
    MPI_Comm comm_purf;
    // for matrix trace
    int starttrrow;
    int starttrcol;
    int lentr;
    // working arrays
    double *work;
    double *work1;
    double *work2;
    double *workm;
    int i;
    // local variables
    double lambda;
    int myrank;
    int nb;
    int nbf;
    double *F_block;
    double *D_block;
    double *FF_block;
    double *DD_block;
    double s;
    double s0;
    double c;
    double c0;
    
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

    if (purf->runpurf == 1)
    {
        // initialization
        nb = purf->nb_purf;
        F_block = purf->F_block;
        D_block = purf->D_block;
        FF_block = purf->FF_block;
        DD_block = purf->DD_block;
        nrows = purf->nrows_purf;
        ncols = purf->ncols_purf;
        starttrrow = purf->tr_srow_purf;
        starttrcol = purf->tr_scol_purf;
        lentr = purf->tr_len_purf;
        nr = purf->nr_purf;
        nc = purf->nc_purf;
        comm_row = purf->comm_purf_row;
        comm_col = purf->comm_purf_col;
        comm_purf = purf->comm_purf;
        nbf = purf->nbf;
            
        MPI_Comm_rank (comm_purf, &myrank);
        MPI_Comm_rank (comm_row, &mycol);
        MPI_Comm_rank (comm_col, &myrow);
        
        work = purf->work;    
        work1 = work;
        work2 = work + nb * nrows;
        workm = purf->D2_block;
        
        // D = D - DD, s = Tr(F * D);    
        #pragma omp parallel for
        for (i = 0; i < nrows * ncols; i++)
        {
            D_block[i] = D_block[i] - DD_block[i];
        }
        my_pdgemm (nbf, nb, F_block, D_block, workm, nrows, ncols,
                   nr, nc, comm_row, comm_col, work1, work2);      
        s0 = 0.0;
        for (i = 0; i < lentr; i++)
        {
            s0 += workm[(i + starttrrow) * ncols + i + starttrcol];
        }
        MPI_Reduce (&s0, &s, 1, MPI_DOUBLE, MPI_SUM, 0, comm_purf);

        // F = FF - F, c = Tr (F * D)
        #pragma omp parallel for
        for (i = 0; i < nrows * ncols; i++)
        {
            FF_block[i] = FF_block[i] - F_block[i];
        }
        my_pdgemm (nbf, nb, FF_block, D_block, workm, nrows, ncols,
                   nr, nc, comm_row, comm_col, work1, work2);
        c0 = 0.0;
        for (i = 0; i < lentr; i++)
        {
            c0 += workm[(i + starttrrow) * ncols + i + starttrcol];
        }
        MPI_Reduce (&c0, &c, 1, MPI_DOUBLE, MPI_SUM, 0, comm_purf);

        // compute lambda
        if (myrank == 0)
        {
            if (c <= -s/2.0)
            {
                lambda = 1.0;
            }
            else
            {
                lambda = -s/(2.0 * c);
            }
        }
        MPI_Bcast (&lambda, 1, MPI_DOUBLE, 0, comm_purf);
         
        // DD = (1 - lambda) * DD + lambda * D
        // F = (1 - lambda) * F + lambda * FF
        #pragma omp parallel for
        for (i = 0; i < nrows * ncols; i++)
        {
            DD_block[i] = DD_block[i] + lambda * D_block[i];
            F_block[i] = F_block[i] + lambda * FF_block[i];
        }
    }

    MPI_Barrier (MPI_COMM_WORLD);
}
