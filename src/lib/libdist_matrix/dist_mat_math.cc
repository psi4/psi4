#include <libparallel/parallel.h>
#include <libmints/vector.h>
#include "dist_mat.h"
#include "scalapack.h"

using namespace psi;

void Distributed_Matrix::diagonalize(Distributed_Matrix& eigvec,
        Vector& eigval)
{
    int ictxt;   // context
    int descA[9];
    int izero = 0;
    int itemp;
    int info;

    Cblacs_get( -1, 0, &ictxt );
    Cblacs_gridinit( &ictxt, "Row", pgrid_.dim_size(0), pgrid_.dim_size(1));

    descinit_( descA, &nrows_, &ncols_, &tile_sz_, &tile_sz_, &izero, &izero, &ictxt, &itemp, &info );
    descinit_( descB, &n, &nrhs, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info );

}

