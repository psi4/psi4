#include "integraltransform.h"
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.hpp>
#define EXTERN
#include "libdpd/dpd.gbl"

using namespace psi;
using namespace boost;

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
IntegralTransform::transform_oei(const shared_ptr<MOSpace> s1, const shared_ptr<MOSpace> s2,
                                 const char *label)
{
    check_initialized();

    double *soInts = init_array(_nTriSo);
    double *moInts = init_array(_nTriMo);
    double *T = init_array(_nTriSo);
    if(_print>4) fprintf(outfile, "The SO basis kinetic energy integrals\n");
    IWL::read_one(_psio.get(), PSIF_OEI, PSIF_SO_T,   T, _nTriSo, 0, _print > 4, outfile);
    if(_print>4) fprintf(outfile, "The SO basis nuclear attraction integrals\n");
    IWL::read_one(_psio.get(), PSIF_OEI, PSIF_SO_V, soInts, _nTriSo, 0, _print > 4, outfile);

    // Add the nuclear and kinetic energy integrals
    for(int n = 0; n < _nTriSo; ++n) soInts[n] += T[n];
    free(T);

    int *order = init_int_array(_nmo);
    // We want to keep Pitzer ordering, so this is just an identity mapping
    for(int n = 0; n < _nmo; ++n) order[n] = n;

    if(_transformationType == Restricted){
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], soInts, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis one-electron integrals\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_OEI, _nTriMo, moInts);
    }else{
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], soInts, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis alpha one-electron integrals\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_A_OEI, _nTriMo, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], soInts, moInts, _Cb[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis beta one-electron integrals\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_B_OEI, _nTriMo, moInts);
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
* @param offset - the point in the full list of SOs where we want to start.  This is
*                 useful for transforming integrals one irrep at a time and in this
*                 case the offset would correspond to the position of the first
*                 orbital in the current irrep.
* @param order - a reordering array to change the order of the output
* @param backtransform - whether this is a forward or backwards transformation
* @param scale - the amount of the existing output buffer to mix into the result
*/

void
IntegralTransform::trans_one(int m, int n, double *input, double *output,
                             double **C, int offset, int* order, bool backtransform, double scale)
{
    // TODO the order argument is actually not used right now.  I don't know that anybody will need it
    // so I haven't bothered so far...
    int dim = (m > n) ? m : n;
    double **TMP0 = block_matrix(dim,dim);
    double **TMP1 = block_matrix(dim,dim);

    for(int p = 0; p < m; ++p){
        for(int q = 0; q <= p; ++q){
            unsigned long int pq = INDEX((p + offset), (q + offset));
            TMP0[p][q] = TMP0[q][p] = input[pq];
        }
    }

    if(backtransform){
        int nc  = m;
        if(m && n) {
            C_DGEMM('n', 't', m, n, m, 1.0, TMP0[0], dim, C[0], nc, 0.0, TMP1[0], dim);
            C_DGEMM('n', 'n', n, n, m, 1.0, C[0], nc, TMP1[0], dim, 0.0, TMP0[0], dim);
        }
    }else{
        int nc  = n;
        if(m && n) {
            C_DGEMM('n', 'n', m, n, m, 1.0, TMP0[0], dim, C[0], nc, 0.0, TMP1[0], dim);
            C_DGEMM('t', 'n', n, n, m, 1.0, C[0], nc, TMP1[0], dim, 0.0, TMP0[0], dim);
        }
    }

    for(int p = 0; p < n; ++p){
        for(int q = 0; q <= p; ++q) {
            size_t P = order[p];
            size_t Q = order[q];
            size_t PQ = INDEX(P,Q);
            output[PQ] = scale * output[PQ] + TMP0[p][q];
        }
    }

    free_block(TMP0);
    free_block(TMP1);

    return;
}


/**
 * Generates the frozen core operator, Fock matrix and one electron integrals
 * in the MO basis by looping over the IWL SO integral file on disk.  The resulting
 * integrals are written to PSIF_SO_OEI and are labelled according to the macros in
 * psifiles.h.  Regardless of any parameters, all integrals are transformed and
 * only IWL format is used to store the results.
 */
void
IntegralTransform::generate_oei()
{
    // Set aside some memory for the frozen core density and frozen core operator
    double *aFzcD  = init_array(_nTriSo);
    double *aFzcOp = init_array(_nTriSo);
    double *aD     = init_array(_nTriSo);
    double *aFock  = init_array(_nTriSo);
    double *aoH    = init_array(_nTriSo);
    double *bFzcD  = aFzcD;
    double *bFzcOp = aFzcOp;
    double *bD     = aD;
    double *bFock  = aFock;
    if(_transformationType != Restricted){
        bFzcD  = init_array(_nTriSo);
        bFzcOp = init_array(_nTriSo);
        bD     = init_array(_nTriSo);
        bFock  = init_array(_nTriSo);
    }

    // Form the Density matrices
    for(int h = 0, soOffset = 0; h < _nirreps; ++h){
        for(int p = 0; p < _sopi[h]; ++p){
            for(int q = 0; q <= p; ++q){
                int pq = INDEX((p + soOffset), (q + soOffset));
                for(int i = 0; i < _frzcpi[h]; ++i)
                    aFzcD[pq] += _Ca[h][p][i] * _Ca[h][q][i];
                for(int i = 0; i < _clsdpi[h] + _openpi[h]; ++i)
                    aD[pq] += _Ca[h][p][i] * _Ca[h][q][i];
                if(_transformationType != Restricted){
                    for(int i = 0; i < _frzcpi[h]; ++i)
                        bFzcD[pq] += _Cb[h][p][i] * _Cb[h][q][i];
                    for(int i = 0; i < _clsdpi[h]; ++i)
                        bD[pq] += _Cb[h][p][i] * _Cb[h][q][i];
                }
            }
        }
        soOffset += _sopi[h];
    }

    double *T = init_array(_nTriSo);
    if(_print>4) fprintf(outfile, "The SO basis kinetic energy integrals\n");
    IWL::read_one(_psio.get(), PSIF_OEI, PSIF_SO_T,   T, _nTriSo, 0, _print > 4, outfile);
    if(_print>4) fprintf(outfile, "The SO basis nuclear attraction integrals\n");
    IWL::read_one(_psio.get(), PSIF_OEI, PSIF_SO_V, aoH, _nTriSo, 0, _print > 4, outfile);

    for(int pq=0; pq < _nTriSo; ++pq){
        aoH[pq] += T[pq];
        aFzcOp[pq] = aoH[pq];
        aFock[pq]  = aoH[pq];
        if(_transformationType != Restricted){
            bFock[pq]  = aoH[pq];
            bFzcOp[pq] = aoH[pq];
        }
    }
    free(T);

    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(_myDPDNum);

    int soIntFile = PSIF_SO_TEI;

    IWL *iwl = new IWL(_psio.get(), soIntFile, _tolerance, 1, 1);
    Label *lblptr = iwl->labels();
    Value *valptr = iwl->values();
    int lastbuf   = iwl->last_buffer();
    for(int index = iwl->index(); index < iwl->buffer_count(); ++index){
        int labelIndex = 4*index;
        int p = abs((int) lblptr[labelIndex++]);
        int q = (int) lblptr[labelIndex++];
        int r = (int) lblptr[labelIndex++];
        int s = (int) lblptr[labelIndex++];
        double value = (double) valptr[index];
        build_fzc_and_fock(p, q, r, s, value, aFzcD, bFzcD,
                           aFzcOp, bFzcOp, aD, bD, aFock, bFock);
    } /* end loop through current buffer */

    /* Now run through the rest of the buffers in the file */
    while(!lastbuf){
        iwl->fetch();
        lastbuf = iwl->last_buffer();
        for(int index = iwl->index(); index < iwl->buffer_count(); ++index){
            int labelIndex = 4*index;
            int p = abs((int) lblptr[labelIndex++]);
            int q = (int) lblptr[labelIndex++];
            int r = (int) lblptr[labelIndex++];
            int s = (int) lblptr[labelIndex++];
            double value = (double) valptr[index];
            build_fzc_and_fock(p, q, r, s, value, aFzcD, bFzcD,
                               aFzcOp, bFzcOp, aD, bD, aFock, bFock);
        } /* end loop through current buffer */
    } /* end loop over reading buffers */
    iwl->set_keep_flag(1);
    delete iwl;

    double *moInts = init_array(_nTriMo);
    int *order = init_int_array(_nmo);
    // We want to keep Pitzer ordering, so this is just an identity mapping
    for(int n = 0; n < _nmo; ++n) order[n] = n;
    if(_print)
        fprintf(outfile, "\tTransforming the one-electron integrals and constructing Fock matrices\n");
    if(_transformationType == Restricted){
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aoH, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis one-electron integrals\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_OEI, _nTriMo, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aFzcOp, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis frozen core operator\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_FZC, _nTriMo, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aFock, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis Fock operator\n");
            print_array(moInts, _nmo, outfile);
        }

        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_FOCK, _nTriMo, aFock);
    }else{
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aoH, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis alpha one-electron integrals\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_A_OEI, _nTriMo, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aoH, moInts, _Cb[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis beta one-electron integrals\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_B_OEI, _nTriMo, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aFzcOp, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis alpha frozen core operator\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_A_FZC, _nTriMo, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], bFzcOp, moInts, _Cb[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis beta frozen core operator\n");
            print_array(moInts, _nmo, outfile);
        }
        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], aFock, moInts, _Ca[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis alpha Fock operator\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_A_FOCK, _nTriMo, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < _nirreps; ++h){
            trans_one(_sopi[h], _mopi[h], bFock, moInts, _Cb[h], soOffset, &(order[moOffset]));
            soOffset += _sopi[h];
            moOffset += _mopi[h];
        }
        if(_print>4){
            fprintf(outfile, "The MO basis beta Fock operator\n");
            print_array(moInts, _nmo, outfile);
        }
        IWL::write_one(_psio.get(), PSIF_OEI, PSIF_MO_B_FOCK, _nTriMo, moInts);
    }
    free(order);
    free(moInts);
    free(aFzcD);
    free(aFzcOp);
    free(aD);
    free(aoH);
    free(aFock);
    if(_transformationType != Restricted){
        free(bFzcD);
        free(bFzcOp);
        free(bD);
        free(bFock);
    }

    dpd_set_default(currentActiveDPD);
}
