#include "integraltransform.h"
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.hpp>

namespace psi{ namespace libtrans{

/**
 * Transforms the one-electron integrals.  This function is currently limited to
 * IWL input and output and Pitzer ordering, regarless of how the parameters passed
 * to the constructor.
 *
 * @param s1 - the MOSpace for the bra
 * @param s2 - the MOSpace for the ket
 *
 * N.B. This need not be called if a two-electron transformation is performed, as
 * the sort_so_tei routine will perform this transformation in addition to the
 * Fock matrix construction.
 */
void
IntegralTransform::transform_oei(shared_ptr<MOSpace> s1, shared_ptr<MOSpace> s2,
                                 const char *label)
{
    double *soInts = init_array(_nTriSo);
    double *moInts = init_array(_nTriMo);
    double *T = init_array(_nTriSo);
    if(_print>4) fprintf(outfile, "The SO basis kinetic energy integrals\n");
    IWL::read_one(_psio, PSIF_OEI, PSIF_SO_T,   T, _nTriSo, 0, _print > 4, outfile);
    if(_print>4) fprintf(outfile, "The SO basis nuclear attraction integrals\n");
    IWL::read_one(_psio, PSIF_OEI, PSIF_SO_V, soInts, _nTriSo, 0, _print > 4, outfile);

    // Add the nuclear and kinetic energy integrals
    for(int n = 0; n < _nTriSo; ++n) soInts[n] += T[n];
    free(T);
    
    int *order = init_int_array(_nmo);
    // We want to keep Pitzer ordering, so this is just an identity mapping
    for(int n = 0; n < _nmo; ++n) order[n] = n;

    if(_transformationType == Restricted){
        for(int n = 0; n < _nTriMo; ++n) moInts[n] = 0.0;
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], soInts, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis one-electron integrals\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio, PSIF_OEI, PSIF_MO_OEI, _nTriMo, moInts);
    }else{
        for(int n = 0; n < _nTriMo; ++n) moInts[n] = 0.0;
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], soInts, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis alpha one-electron integrals\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio, PSIF_OEI, PSIF_MO_A_OEI, _nTriMo, moInts);

        for(int n = 0; n < _nTriMo; ++n) moInts[n] = 0.0;
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], soInts, moInts, _Cb[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis beta one-electron integrals\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio, PSIF_OEI, PSIF_MO_B_OEI, _nTriMo, moInts);
    }
    free(order);
    free(order);
    free(moInts);
    free(soInts);
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

