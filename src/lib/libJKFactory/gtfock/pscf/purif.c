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
#include <sys/time.h>

#include "pdgemm.h"
#include "purif.h"


#define MAX_PURF_ITERS 200
#define MIN(a, b)    ((a) < (b) ? (a) : (b))
#define MAX(a, b)    ((a) > (b) ? (a) : (b))

static void config_purif (purif_t * purif, int purif_offload)
{
    int nbf;
    int nrows;
    int ncols;
    int nb;
    int nprow;
    int npcol;
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
    int izero = 0;

    nbf = purif->nbf;
    nprow = purif->nprow_purif;
    npcol = purif->npcol_purif;
    purif->nb_purif = nb = MIN (nbf / nprow, nbf / npcol);
    comm_col = purif->comm_purif_col;
    comm_row = purif->comm_purif_row;

    int coords[3];
    int myrank;
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Cart_coords (purif->comm_purif, myrank, 3, coords);
    int myrow = coords[0];
    int mycol = coords[1];
    int mygrd = coords[2];

    nrows = numroc_ (&nbf, &nb, &myrow, &izero, &nprow);
    ncols = numroc_ (&nbf, &nb, &mycol, &izero, &npcol);
    purif->nrows_purif = nrows;
    purif->ncols_purif = ncols;

    // positions of partitions
    purif->nr_purif = (int *) malloc (sizeof (int) * purif->nprow_purif);
    purif->nc_purif = (int *) malloc (sizeof (int) * purif->npcol_purif);
    assert (purif->nr_purif != NULL);
    assert (purif->nc_purif != NULL);
    nr = purif->nr_purif;
    nc = purif->nc_purif;

    // get nr and sr
    MPI_Allgather (&nrows, 1, MPI_INT, nr, 1, MPI_INT, comm_col);
    startrow = 0;
    for (int i = 0; i < myrow; i++)
    {
        startrow += nr[i];
    }
    endrow = startrow + nrows - 1;
    purif->srow_purif = startrow;

    // get nc and sc
    MPI_Allgather (&ncols, 1, MPI_INT, nc, 1, MPI_INT, comm_row);
    startcol = 0;
    for (int i = 0; i < mycol; i++)
    {
        startcol += nc[i];
    }
    endcol = startcol + ncols - 1;
    purif->scol_purif = startcol;

    // for matrix trace
    start = MAX (startcol, startrow);
    purif->tr_len_purif = MIN (endcol, endrow) - start + 1;
    purif->tr_scol_purif = start - startcol;
    purif->tr_srow_purif = start - startrow;
    purif->istr_purif = (purif->tr_len_purif > 0);

    // create local arrays
    // purif->ldx = (ncols + ALIGNSIZE - 1)/ALIGNSIZE * ALIGNSIZE;
    purif->ldx = ncols;
    meshsize = nrows * ncols;
    purif->meshsize = meshsize;
    purif->X_block = (double *) mkl_malloc (meshsize * sizeof (double), 64);
    assert (purif->X_block != NULL);
    purif->S_block = (double *) mkl_malloc (meshsize * sizeof (double), 64);
    assert (purif->S_block != NULL);
    purif->H_block = (double *) mkl_malloc (meshsize * sizeof (double), 64);
    assert (purif->H_block != NULL);
    purif->F_block = (double *) mkl_malloc (meshsize * sizeof (double), 64);
    assert (purif->F_block != NULL);
    purif->D_block = (double *) mkl_malloc (meshsize * sizeof (double), 64);
    assert (purif->D_block != NULL);
    purif->D2_block = (double *) mkl_malloc (meshsize * sizeof (double), 64);
    assert (purif->D2_block != NULL);
    purif->D3_block = (double *) mkl_malloc (meshsize * sizeof (double), 64);
    assert (purif->D3_block != NULL);
    // working space for purification
    purif->diis_vecs =
        (double *) mkl_malloc (MAX_DIIS * meshsize * sizeof (double), 64);
    assert (purif->diis_vecs != NULL);
    purif->F_vecs =
        (double *) mkl_malloc (MAX_DIIS * meshsize * sizeof (double), 64);
    assert (purif->diis_vecs != NULL);
    purif->len_diis = 0;
    purif->bmax = DBL_MIN;
    purif->bmax_id = -1;
    assert (MAX_DIIS > 1);
    for (int i = 0; i < LDBMAT; i++)
    {
        for (int j = 0; j < LDBMAT; j++)
        {
            purif->b_mat[i * LDBMAT + j] = (i != j) ? -1.0 : 0.0;
        }
    }

    #pragma omp parallel for schedule(static)
    #pragma simd
    for(int i=0; i < nrows * ncols; i++)
    {
        purif->D_block[i] = 0.0;
        purif->D2_block[i] = 0.0;
        purif->D3_block[i] = 0.0;
    }
    allocate_tmpbuf (nrows, ncols, nr, nc, &(purif->tmpbuf));
    memsize = (2.0 * MAX_DIIS + 14.0) * meshsize * sizeof (double);
    if (myrow == 0 && mycol == 0 && mygrd == 0)
    {
        printf ("  CPU uses %.3f MB\n", memsize / 1024.0 / 1024.0);
    }
}


purif_t *create_purif(BasisSet_t basis, int nprow_purif,
                      int npcol_purif, int npgrd_purif)
{
    int myrank;

    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    if (myrank == 0)
    {
        printf ("Initializing purification ...\n");
    }

    // create purif
    purif_t *purif = (purif_t *) malloc (sizeof (purif_t));
    assert (purif != NULL);
    purif->nbf = CInt_getNumFuncs (basis);
    purif->nobtls = CInt_getNumOccOrb (basis);
    purif->nprow_purif = nprow_purif;
    purif->npcol_purif = npcol_purif;
    purif->npgrd_purif = npgrd_purif;
    purif->np_purif = nprow_purif * npcol_purif * npgrd_purif;

    // set node types
    purif->runpurif = 0;
    if (myrank < purif->np_purif)
    {
        purif->rundgemm = 1;
    }
    else
    {
        purif->rundgemm = 0;
    }

    // initialize communicators
    int flag_purif = (myrank < purif->np_purif);
    MPI_Comm comm0;
    MPI_Comm_split (MPI_COMM_WORLD, flag_purif, myrank, &comm0);

    if (purif->rundgemm == 1)
    {
        int ndims = 3;
        int dim_size[3];
        int periods[3];
        int coords[3];
        dim_size[0] = nprow_purif;
        dim_size[1] = npcol_purif;
        dim_size[2] = npgrd_purif;
        periods[0] = 0;
        periods[1] = 0;
        periods[2] = 0;
        int reorder = 1;
        MPI_Cart_create (comm0, ndims, dim_size, periods, reorder,
                         &(purif->comm_purif));
        MPI_Cart_coords (purif->comm_purif, myrank, ndims, coords);
        int myrow = coords[0];
        int mycol = coords[1];
        int mygrd = coords[2];

        if (mygrd == 0)
        {
            purif->runpurif = 1;
        }
        int belongsR[3] = { 0, 1, 0 };
        MPI_Cart_sub (purif->comm_purif, belongsR, &(purif->comm_purif_row));
        int belongsC[3] = { 1, 0, 0 };
        MPI_Cart_sub (purif->comm_purif, belongsC, &(purif->comm_purif_col));
        int belongsG[3] = { 0, 0, 1 };
        MPI_Cart_sub (purif->comm_purif, belongsG, &(purif->comm_purif_grd));
        MPI_Comm_split (purif->comm_purif, mygrd, myrow * npcol_purif + mycol,
                        &(purif->comm_purif_plane));
        config_purif (purif, 0);
    }
    if (myrank == 0)
    {
        printf ("  Done\n");
    }
    return purif;
}


void destroy_purif (purif_t * purif)
{
    if (purif->rundgemm == 1)
    {
        dealloc_tmpbuf (&(purif->tmpbuf));
        
        MPI_Comm_free (&(purif->comm_purif));
        MPI_Comm_free (&(purif->comm_purif_row));
        MPI_Comm_free (&(purif->comm_purif_col));
        free (purif->nr_purif);
        free (purif->nc_purif);
        mkl_free (purif->H_block);
        mkl_free (purif->X_block);
        mkl_free (purif->S_block);
        mkl_free (purif->F_block);
        mkl_free (purif->D_block);
        mkl_free (purif->D3_block);
        mkl_free (purif->D2_block);
        mkl_free (purif->F_vecs);
        mkl_free (purif->diis_vecs);
    }
    free (purif);
}


int compute_purification(purif_t * purif, double *F_block, double *D_block)
{
    struct timeval tv1;
    struct timeval tv2;
    struct timeval tv3;
    struct timeval tv4;
    int it;
    purif->timedgemm = 0.0;
    purif->timepdgemm = 0.0;
    purif->timepass = 0.0;
    purif->timetr = 0.0;

    if (purif->rundgemm == 1) {
        gettimeofday(&tv3, NULL);
        // initialization
        int nrows = purif->nrows_purif;
        int ncols = purif->ncols_purif;
        int startrow = purif->srow_purif;
        int startcol = purif->scol_purif;
        int starttrrow = purif->tr_srow_purif;
        int starttrcol = purif->tr_scol_purif;
        int lentr = purif->tr_len_purif;
        int *nr = purif->nr_purif;
        int *nc = purif->nc_purif;
        MPI_Comm comm_row = purif->comm_purif_row;
        MPI_Comm comm_col = purif->comm_purif_col;
        MPI_Comm comm_grd = purif->comm_purif_grd;
        MPI_Comm comm_purif = purif->comm_purif_plane;
        MPI_Comm comm0 = purif->comm_purif;
        int nbf = purif->nbf;
        int nobtls = purif->nobtls;
        double *X_block = purif->X_block;
        double *D2_block = purif->D2_block;
        double *D3_block = purif->D3_block;
        double *workm = D2_block;
        int coords[3];
        int myrank;
        MPI_Comm_rank (purif->comm_purif, &myrank);
        MPI_Cart_coords (purif->comm_purif, myrank, 3, coords);
        int myrow = coords[0];
        int mycol = coords[1];
        int mygrd = coords[2];
        tmpbuf_t tmpbuf = purif->tmpbuf;
        
        if (purif->runpurif == 1) {
            // compute eigenvalue estimates using Gershgorin (F is symmetric)
            // offdF = sum(abs(F))' - abs(diag(F));
            // diagF = diag(F);
            // hmin = min(diagF - offdF);
            // hmax = max(diagF + offdF);
            double _h[2 * ncols]; // offDF, diagF
            double h[2 * ncols];  // hmin, hmax       
            for (int i = 0; i < ncols; i++) {
                _h[i] = 0.0;
                _h[i + ncols] = 0.0;
                for (int j = 0; j < nrows; j++) {
                    _h[i] += fabs(F_block[i + j * ncols]);
                    if (j + startrow == i + startcol) {
                        _h[i + ncols] = F_block[i + j * ncols];
                    }
                }
                _h[i] = _h[i] - fabs(_h[i + ncols]);
                double tmp = _h[i + ncols] + _h[i];
                _h[i] = _h[i + ncols] - _h[i];
                _h[i + ncols] = tmp;
            }
            MPI_Reduce(_h, h, 2 * ncols, MPI_DOUBLE, MPI_SUM, 0, comm_col);

            double _hmax;
            double _hmin;
            double hmax;
            double hmin;
            if (myrow == 0) {
                _hmin = DBL_MAX;
                _hmax = -DBL_MAX;
                for (int i = 0; i < ncols; i++) {
                    _hmin = h[i] > _hmin ? _hmin : h[i];
                    _hmax = h[i + ncols] < _hmax ? _hmax : h[i + ncols];
                }
                MPI_Reduce(&_hmin, &hmin, 1, MPI_DOUBLE, MPI_MIN, 0, comm_row);
                MPI_Reduce(&_hmax, &hmax, 1, MPI_DOUBLE, MPI_MAX, 0, comm_row);
            }
            MPI_Bcast(&hmin, 1, MPI_DOUBLE, 0, comm_purif);
            MPI_Bcast(&hmax, 1, MPI_DOUBLE, 0, comm_purif);
            // define constants, dependent on F
            // in the following:
            // 5 = no of occupied orbitals
            // 7 = no of spatial basis function
            // (each corresponds to 2 electrons for RHF)
            // mu_bar = trace_dense_matrix(F)/7;
            double trF;
            double _trF;
            _trF = 0.0;
            for (int i = 0; i < lentr; i++) {
                _trF += F_block[(i + starttrrow) * ncols + i + starttrcol];
            }
            MPI_Allreduce(&_trF, &trF, 1, MPI_DOUBLE, MPI_SUM, comm_purif);

            double mu_bar = trF / (double) nbf;
            // lambda = min([ 5/(hmax - mu_bar), (7-5)/(mu_bar - hmin) ]);
            double lambda = MIN((double) nobtls / (hmax - mu_bar),
                                (double) (nbf - nobtls) / (mu_bar - hmin));
            if (myrank == 0) {
                printf("mu_bar = %le, lambda = %le,"
                       " hmax = %le, hmin = %le, nobtls = %d\n",
                       mu_bar, lambda, hmax, hmin, nobtls);
            }
            
            // initial "guess" for density matrix
            // D = (lambda*mu_bar/7 + 5/7)*eye(7) - (lambda/7)*D;
            for (int i = 0; i < nrows * ncols; i++) {
                D_block[i] = F_block[i] * (-lambda / nbf);
            }
            for (int i = 0; i < lentr; i++) {
                D_block[(i + starttrrow) * ncols + i + starttrcol] +=
                    lambda * mu_bar / (double) nbf + (double) nobtls / nbf;
            }
        } /* if (purif->runpurif == 1) */

        // McWeeny purification
        // convergence appears slow at first, before accelerating at end
        for (it = 0; it < MAX_PURF_ITERS; it++) {
            gettimeofday(&tv1, NULL);
            double tmp_time;
            pdgemm3D(myrow, mycol, mygrd, comm_row, comm_col, comm_grd,
                     comm0, nr, nc, nrows, ncols, D_block, D2_block,
                     D3_block, &tmpbuf, &tmp_time);
            gettimeofday(&tv2, NULL);
            purif->timepdgemm += (tv2.tv_sec - tv1.tv_sec) +
                (tv2.tv_usec - tv1.tv_usec) / 1000.0 / 1000.0;
            purif->timedgemm += tmp_time;

            gettimeofday(&tv1, NULL);

            double errnorm;
            if (purif->runpurif == 1) {
                // stopping criterion
                // errnorm = norm(D-D2, 'fro');
                double _errnorm = 0.0;
                #pragma omp parallel for reduction(+: _errnorm)
                for (int i = 0; i < nrows * ncols; i++) {  
                    _errnorm += (D_block[i] - D2_block[i]) *
                            (D_block[i] - D2_block[i]);
                }
                MPI_Reduce(&_errnorm, &errnorm, 1,
                           MPI_DOUBLE, MPI_SUM, 0, comm_purif);
                if (myrank == 0) {
                    errnorm = sqrt(errnorm);
                }

                // a cheaper stopping criterion may be to
                // check the trace of D*D
                // and stop when it is close to no. occupied orbitals
                // (5 in this case)
                // fprintf('trace D*D //f\n', trace(D*D);
                // might be possible to "lag" the computation
                // of c by one iteration
                // so that the global communication for
                // traces can be overlapped.
                // Note: c appears to converge to 0.5           
                // c = trace(D2-D3) / trace(D-D2);
                double c;
                double _tr = 0.0;
                double _tr2 = 0.0;
                double tr;
                double tr2;
                #pragma omp parallel for reduction(+: _tr, _tr2)
                for (int i = 0; i < lentr; i++) {
                    _tr += D2_block[(i + starttrrow) * ncols + i + starttrcol] -
                        D3_block[(i + starttrrow) * ncols + i + starttrcol];
                    _tr2 += D_block[(i + starttrrow) * ncols + i + starttrcol] -
                        D2_block[(i + starttrrow) * ncols + i + starttrcol];
                }

                // Jeff: This is the result of fusion of
                //       two Reduce and one Bcast
                //       calls that were used to determine tr and tr2, hence c.
                double itmp[2], otmp[2];
                itmp[0] = _tr;
                itmp[1] = _tr2;
                MPI_Allreduce(itmp, otmp, 2, MPI_DOUBLE, MPI_SUM, comm_purif);
                tr = otmp[0];
                tr2 = otmp[1];
                c = tr / tr2;
                if (c < 0.5) {
                    #pragma omp parallel for
                    for (int i = 0; i < nrows * ncols; i++) {
                        // D = ((1-2*c)*D + (1+c)*D2 - D3) / (1-c);
                        D_block[i] = ((1.0 - 2.0 * c) * D_block[i] +
                            (1.0 + c) * D2_block[i] - D3_block[i]) / (1.0 - c);
                    }
                } else {
                    #pragma omp parallel for
                    for (int i = 0; i < nrows * ncols; i++) {
                        // D = ((1+c)*D2 - D3) / c;
                        D_block[i] =
                            ((1.0 + c) * D2_block[i] - D3_block[i]) / c;
                    }
                }
                
            }          
            MPI_Bcast(&errnorm, 1, MPI_DOUBLE, 0, comm0);
            if (errnorm < 1e-11) {
                break;
            }
            gettimeofday(&tv2, NULL);
            purif->timetr += (tv2.tv_sec - tv1.tv_sec) +
                (tv2.tv_usec - tv1.tv_usec) / 1000.0 / 1000.0;
        }

        gettimeofday(&tv1, NULL);
        double tmp_time;
        pdgemm3D_2(myrow, mycol, mygrd, comm_row, comm_col, comm_grd,
                   comm0, nr, nc, nrows, ncols,
                   X_block, D_block, workm, &tmpbuf, &tmp_time);
        purif->timedgemm += tmp_time;
        pdgemm3D_2(myrow, mycol, mygrd, comm_row, comm_col, comm_grd,
                   comm0, nr, nc, nrows, ncols,
                   workm, X_block, D_block, &tmpbuf, &tmp_time);
        purif->timedgemm += tmp_time;
        gettimeofday(&tv2, NULL);
        purif->timepdgemm += (tv2.tv_sec - tv1.tv_sec) +
            (tv2.tv_usec - tv1.tv_usec) / 1000.0 / 1000.0;

        gettimeofday(&tv4, NULL);
        purif->timepass += (tv4.tv_sec - tv3.tv_sec) +
            (tv4.tv_usec - tv3.tv_usec) / 1000.0 / 1000.0;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return it;
}


void compute_diis (PFock_t pfock, purif_t * purif,
                   double *D_block, double *F_block, int iter)
{
    int nrows = purif->nrows_purif;
    int ncols = purif->ncols_purif;
    int meshsize = purif->meshsize;
    double *X_block = purif->X_block;
    double *S_block = purif->S_block;
    double *workm = purif->D2_block;
    double *b_mat = purif->b_mat;

    double *F_vecs = purif->F_vecs;
    double *diis_vecs = purif->diis_vecs;
    int *nr = purif->nr_purif;
    int *nc = purif->nc_purif;
    int myrank;
    MPI_Comm comm_row = purif->comm_purif_row;
    MPI_Comm comm_col = purif->comm_purif_col;
    MPI_Comm comm_grd = purif->comm_purif_grd;
    MPI_Comm comm_purif = purif->comm_purif_plane;
    MPI_Comm comm0 = purif->comm_purif;
    tmpbuf_t tmpbuf = purif->tmpbuf;

    if (purif->rundgemm == 1) {    
        int coords[3];
        MPI_Comm_rank (purif->comm_purif, &myrank);
        MPI_Cart_coords (purif->comm_purif, myrank, 3, coords);
        int myrow = coords[0];
        int mycol = coords[1];
        int mygrd = coords[2];
        double *cur_F = F_block;
        int cur_idx;
        double *cur_diis;    
        if (iter > 1) {
            if (purif->len_diis < MAX_DIIS) {
                cur_idx = purif->len_diis;
                cur_diis = &(diis_vecs[cur_idx * meshsize]);
                cur_F = &(F_vecs[cur_idx * meshsize]);
                purif->len_diis++;
            } else {
                cur_idx = purif->bmax_id;
                cur_diis = &(diis_vecs[cur_idx * meshsize]);
                cur_F = &(F_vecs[cur_idx * meshsize]);
            }
            // Ctator = X*(F*D*S - S*D*F)*X;
            pdgemm3D_2(myrow, mycol, mygrd, comm_row, comm_col, comm_grd,
                       comm0, nr, nc, nrows, ncols,
                       F_block, D_block, workm, &tmpbuf, NULL);
            pdgemm3D_2(myrow, mycol, mygrd, comm_row, comm_col, comm_grd,
                       comm0, nr, nc, nrows, ncols,
                       workm, S_block, cur_diis, &tmpbuf, NULL);
            if (purif->runpurif == 1) {
                int dest;
                coords[0] = mycol;
                coords[1] = myrow;
                coords[2] = mygrd;               
                MPI_Cart_rank(purif->comm_purif, coords, &dest);
                MPI_Sendrecv(cur_diis, nrows * ncols,
                             MPI_DOUBLE, dest, dest,
                             workm, nrows * ncols,
                             MPI_DOUBLE, dest, MPI_ANY_TAG,
                             purif->comm_purif, MPI_STATUS_IGNORE);
                // F*D*S - (F*D*S)'       
                for (int i = 0; i < nrows; i++) {
                    for (int j = 0; j < ncols; j++) {
                        cur_diis[i * ncols + j] -= workm[j * nrows + i];
                    }
                }
            }
            pdgemm3D_2(myrow, mycol, mygrd, comm_row, comm_col, comm_grd,
                       comm0, nr, nc, nrows, ncols, X_block, cur_diis,
                       workm, &tmpbuf, NULL);
            pdgemm3D_2(myrow, mycol, mygrd, comm_row, comm_col, comm_grd,
                       comm0, nr, nc, nrows, ncols, workm, X_block,
                       cur_diis, &tmpbuf, NULL);
            if (purif->runpurif == 1) {
                // b_mat(i,j) = dot(vecs(:,i), vecs(:,j));
                double _dot[LDBMAT];
                for (int i = 0; i < purif->len_diis; i++) {
                    double *diis2 = &(diis_vecs[i * meshsize]);
                    _dot[i] = 0.0;
                    for (int j = 0; j < nrows; j++) {
                        for (int k = 0; k < ncols; k++) {
                            _dot[i] +=
                                cur_diis[j * ncols + k] * diis2[j * ncols + k];
                        }
                    }
                } /* end for */
                // update b_mat on rank 0          
                MPI_Reduce(_dot, &(b_mat[cur_idx * LDBMAT]),
                           purif->len_diis, MPI_DOUBLE, MPI_SUM, 0, comm_purif);
                if (myrank == 0) {
                    purif->bmax = -DBL_MAX;
                    for (int i = 0; i < purif->len_diis; i++) {
                        b_mat[i * LDBMAT + cur_idx] =
                            b_mat[cur_idx * LDBMAT + i];
                        if (purif->bmax < b_mat[i * LDBMAT + i]) {
                            purif->bmax = b_mat[i * LDBMAT + i];
                            purif->bmax_id = i;
                        }
                    }
                }
            } /* if (purif->runpurif == 1) */
            MPI_Bcast (&(purif->bmax_id), 1, MPI_DOUBLE, 0, comm0);
        } /* if (iter > 1) */

        // F = X*F*X;
        pdgemm3D_2(myrow, mycol, mygrd, comm_row, comm_col, comm_grd,
                   comm0, nr, nc, nrows, ncols,
                   X_block, F_block, workm, &tmpbuf, NULL);
        pdgemm3D_2(myrow, mycol, mygrd, comm_row, comm_col, comm_grd,
                   comm0, nr, nc, nrows, ncols, workm,
                   X_block, cur_F, &tmpbuf, NULL);
        // extrapolate
        if (iter > 1) {
            if (purif->runpurif == 1) {
                double coeffs[LDBMAT];
                // rhs = zeros(m+1,1);
                // rhs(m+1,1) = -1;
                // coeffs = inv(b_mat) * rhs;
                if (myrank == 0) {
                    int sizeb = purif->len_diis + 1;
                    __declspec(align (64)) double b_inv[LDBMAT * LDBMAT];
                    __declspec(align (64)) int ipiv[LDBMAT];
                    memcpy(b_inv, b_mat, LDBMAT * LDBMAT * sizeof (double));
                    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, sizeb, sizeb, b_inv,
                                   LDBMAT, ipiv);
                    LAPACKE_dgetri(LAPACK_ROW_MAJOR, sizeb, b_inv, LDBMAT,
                                   ipiv);
                    for (int i = 0; i < sizeb; i++) {
                        coeffs[i] = -b_inv[i * LDBMAT + sizeb - 1];
                    }
                }
                MPI_Bcast(coeffs, purif->len_diis, MPI_DOUBLE, 0, comm_purif);

                // F = 0
                // for j = 1:m
                //     F = F + coeffs(j)* F_vecs(j);
                memset(F_block, 0, sizeof (double) * meshsize);
                for (int i = 0; i < purif->len_diis; i++) {
                    double *F_vec = &(F_vecs[i * meshsize]);
                    for (int j = 0; j < meshsize; j++) {
                        F_block[j] += coeffs[i] * F_vec[j];
                    }
                }
            } /* if (purif->runpurif == 1) */
        }
    } /* if (purif->rundgemm == 1) */

    MPI_Barrier(MPI_COMM_WORLD);
}


#if 0
static void peig(int ga_A, int ga_B, int n, int nprow, int npcol, double *eval)
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
        printf("  pdsyev_ takes %.3f secs\n", t2 - t1);
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


void compute_eigensolve(int ga_tmp, purif_t * purif,
                        double *F_block, int nprow, int npcol)
{
    int nbf = purif->nbf;
    if (purif->runpurif == 1) {
        int lo[2];
        int hi[2];
        int ld;
        lo[0] = purif->srow_purif;
        hi[0] = purif->srow_purif + purif->nrows_purif - 1;
        lo[1] = purif->scol_purif;
        hi[1] = purif->scol_purif + purif->ncols_purif - 1;
        ld = purif->ldx;
        NGA_Put(ga_tmp, lo, hi, F_block, &ld);
    }

    double *eval = (double *)mkl_malloc(nbf * sizeof (double), 64);
    assert(eval != NULL);
    peig(ga_tmp, ga_tmp, nbf, nprow, npcol, eval);

    mkl_free(eval);
}

#endif