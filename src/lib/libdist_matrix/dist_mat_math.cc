#include <algorithm>
#include <libparallel/parallel.h>
#include <libmints/vector.h>
#include <libmints/matrix.h>
#include "dist_mat.h"
#include "scalapack.h"

using namespace psi;

#if defined(HAVE_MADNESS)

void Distributed_Matrix::diagonalize(Distributed_Matrix& eigvec, Vector& eigval)
{
#if defined(HAVE_SCALAPACK)
    int ictxt;   // context
    int descA[9], descB[9];
    int izero = 0, ione = 1;
    int itemp;
    int info;
    int nprow, npcol, myrow, mycol;
    int n, nb, np, nq, nqrhs, nrhs;

    Cblacs_get( 0, 0, &ictxt );

    printf("ictxt = %d\n", ictxt);
    printf("pgrid_.dim_size(0) = %d\n", pgrid_.dim_size(0));
    printf("pgrid_.dim_size(1) = %d\n", pgrid_.dim_size(1));

    Cblacs_gridinit( &ictxt, "Row", pgrid_.dim_size(0), pgrid_.dim_size(1));
    Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

    printf("nprow %d\nnpcol %d\nmyrow %d\nmycol %d\n", nprow, npcol, myrow, mycol);

    /*
     *
     * Compute the size of the local matrices (thanks to numroc)
     *
     */
    np    = numroc_( &nrows_, &tile_sz_, &myrow, &izero, &nprow );
    nq    = numroc_( &ncols_, &tile_sz_, &mycol, &izero, &npcol );
    nqrhs = numroc_( &eigvec.ncols_, &tile_sz_, &mycol, &izero, &npcol );

    SharedMatrix A(new Matrix(np,nq));
    SharedMatrix B(new Matrix(np,nq));

    A->set(1.0);
    A->set_name("A");
    B->set_name("B");

    printf("np %d  nrows_ %d  tile_sz_ %d  myrow %d  izero %d  nprow %d\n", np, nrows_, tile_sz_, myrow, izero, nprow);
    printf("nq %d  ncols_ %d  tile_sz_ %d mycol %d  izero %d  npcol %d\n", nq, ncols_, tile_sz_, mycol, izero, npcol);
    printf("nqrhs %d  ev.nrows_ %d  tile_sz_ %d mycol %d izero %d npcol %d\n", nqrhs, eigvec.nrows_, tile_sz_, mycol, izero, npcol);

    /*
     *
     * Initialize the array descriptor for the matrix A and B
     *
     */

    itemp = std::max( 1, np );
    descinit_( descA, &nrows_, &ncols_, &tile_sz_, &tile_sz_, &izero, &izero, &ictxt, &itemp, &info );
    descinit_( descB, &eigvec.nrows_, &eigvec.ncols_, &eigvec.tile_sz_, &eigvec.tile_sz_,
              &izero, &izero, &ictxt, &itemp, &info );

    /*
     **********************************************************************
     *     Call ScaLAPACK PDGESV routine
     **********************************************************************
     */
    // pivoting information
    int *ippiv = new int[np * tile_sz_];
    int lwork = 5 * nrows_ * (np + 1) +1;
    double* work = new double[lwork];
    //pdgesv_( &nrows_, &eigvec.ncols_, ptr(), &ione, &ione, descA, ippiv, eigvec.ptr(), &ione, &ione, descB, &info );
    //pdsyev_("N", "U", &nrows_, ptr(), &ione, &ione, descA, eigval.pointer(), eigvec.ptr(), &ione, &ione,
    //        descB, work, &lwork, &info);
    pdsyev_("N", "U", &eigvec.ncols_, A->pointer()[0], &ione, &ione, descA, eigval.pointer(), B->pointer()[0], &ione, &ione,
            descB, work, &lwork, &info);
    delete[] ippiv;
    delete[] work;

    Cblacs_barrier(ictxt, "A");
    Cblacs_exit(1);

    WorldComm->sync();
#else
#   warning ScaLAPACK is not available.
#endif
}

#endif

