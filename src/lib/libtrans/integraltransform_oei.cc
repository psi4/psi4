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
IntegralTransform::transform_oei(shared_ptr<MOSpace> s1, shared_ptr<MOSpace> s2)
{
    return;
}


/**
* trans_one(): Transform a packed symmetric matrix.
*
* int m: input matrix row dimension
* int n: output matrix row dimension
* double *input: pointer to input integrals (the lower-triange of a symmetric matrix)
* double *output: pointer to output integrals (the lower-triangle of a symmetric matrix)
* double **C: transformation matrix (rectangular)
* int nc: column dimension of C
* int *order: reordering array for transformed indices
*
* Stolen, by ACS, from the transqt2 code written by TDC, 7/06
*/

void
IntegralTransform::trans_one(int m, int n, double *input, double *output,
                                double **C, int nc, int *order)
{
    int dim = (m > n) ? m : n;
    double **TMP0 = block_matrix(dim,dim);
    double **TMP1 = block_matrix(dim,dim);

    for(int p = 0, pq = 0; p < m; ++p){
        for(int q = 0; q <= p; ++q,++pq){
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

