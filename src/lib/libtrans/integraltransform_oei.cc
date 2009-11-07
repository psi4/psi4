#include "integraltransform.h"
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

namespace psi{ namespace libtrans{

/**
* Transforms the one-electron integrals
* @param s1 - the MOSpace for the bra
* @param s2 - the MOSpace for the ket
*/
void
IntegralTransform::transform_oei(shared_ptr<MOSpace> s1, shared_ptr<MOSpace> s2,
                                 const char *label)
{
    return;
}


/**
* Transforms a packed symmetric matrix.
*
* @param m - input matrix row dimension
* @param n - output matrix row dimension
* @param input - pointer to input integrals (the lower-triangle of a symmetric matrix)
* @param pointer to output integrals (the lower-triangle of a symmetric matrix)
* @param C transformation matrix (rectangular, m X n)
* @param soOffset - the point in the full list of SOs where we want to start.  This is
*                   useful for transforming integrals one irrep at a time and in this
*                   case the offset would correspond to the position of the first
*                   orbital in the current irrep.
* @param order - a reordering array to change the order of the output
*/

void
IntegralTransform::trans_one(int m, int n, double *input, double *output, 
                             double **C, int soOffset, int* order)
{
    int dim = (m > n) ? m : n;
    int nc  = n;
    double **TMP0 = block_matrix(dim,dim);
    double **TMP1 = block_matrix(dim,dim);

    for(int p = 0; p < m; ++p){
        for(int q = 0; q <= p; ++q){
            unsigned long int pq = INDEX((p + soOffset), (q + soOffset));
            TMP0[p][q] = TMP0[q][p] = input[pq];
        }
    }

    if(m && n) {
        C_DGEMM('n', 'n', m, n, m, 1.0, TMP0[0], dim, C[0], nc, 0.0, TMP1[0], dim);
        C_DGEMM('t', 'n', n, n, m, 1.0, C[0], nc, TMP1[0], dim, 0.0, TMP0[0], dim);
    }

    for(int p = 0; p < n; ++p){
        for(int q = 0; q <= p; ++q) {
            unsigned long int pq = INDEX(order[p], order[q]);
            output[pq] = TMP0[p][q];
        }
    }

    free_block(TMP0);
    free_block(TMP1);
    
    return;
}

}} // End namespaces

