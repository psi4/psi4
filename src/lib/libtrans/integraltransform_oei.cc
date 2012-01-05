#include "integraltransform.h"
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libmints/matrix.h>
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

    double *soInts = init_array(nTriSo_);
    double *moInts = init_array(nTriMo_);
    double *T = init_array(nTriSo_);
    if(print_>4) fprintf(outfile, "The SO basis kinetic energy integrals\n");
    IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_T,   T, nTriSo_, 0, print_ > 4, outfile);
    if(print_>4) fprintf(outfile, "The SO basis nuclear attraction integrals\n");
    IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_V, soInts, nTriSo_, 0, print_ > 4, outfile);

    // Add the nuclear and kinetic energy integrals
    for(int n = 0; n < nTriSo_; ++n) soInts[n] += T[n];
    free(T);

    int *order = init_int_array(nmo_);
    // We want to keep Pitzer ordering, so this is just an identity mapping
    for(int n = 0; n < nmo_; ++n) order[n] = n;

    if(transformationType_ == Restricted){
        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], soInts, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The MO basis one-electron integrals\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_OEI, nTriMo_, moInts);
    }else{
        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], soInts, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The MO basis alpha one-electron integrals\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_A_OEI, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCb = Cb_->pointer(h);
            trans_one(sopi_[h], mopi_[h], soInts, moInts, pCb, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The MO basis beta one-electron integrals\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_B_OEI, nTriMo_, moInts);
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
    int nc;
    if(backtransform){
        nc  = m;
        if(m && n) {
            C_DGEMM('n', 't', m, n, m, 1.0, TMP0[0], dim, C[0], nc, 0.0, TMP1[0], dim);
            C_DGEMM('n', 'n', n, n, m, 1.0, C[0], nc, TMP1[0], dim, 0.0, TMP0[0], dim);
        }
    }else{
        nc  = n;
        if(m && n) {
            C_DGEMM('n', 'n', m, n, m, 1.0, TMP0[0], dim, C[0], nc, 0.0, TMP1[0], dim);
            C_DGEMM('t', 'n', n, n, m, 1.0, C[0], nc, TMP1[0], dim, 0.0, TMP0[0], dim);
        }
    }

    for(int p = 0; p < nc; ++p){
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
#if 0
    // Set aside some memory for the frozen core density and frozen core operator
    double *aFzcD  = init_array(nTriSo_);
    double *aFzcOp = init_array(nTriSo_);
    double *aD     = init_array(nTriSo_);
    double *aFock  = init_array(nTriSo_);
    double *aoH    = init_array(nTriSo_);
    double *bFzcD  = aFzcD;
    double *bFzcOp = aFzcOp;
    double *bD     = aD;
    double *bFock  = aFock;
    if(transformationType_ != Restricted){
        bFzcD  = init_array(nTriSo_);
        bFzcOp = init_array(nTriSo_);
        bD     = init_array(nTriSo_);
        bFock  = init_array(nTriSo_);
    }

    // Form the Density matrices
    for(int h = 0, soOffset = 0; h < nirreps_; ++h){
        double **pCa = Ca_->pointer(h);
        double **pCb = Cb_->pointer(h);
        for(int p = 0; p < sopi_[h]; ++p){
            for(int q = 0; q <= p; ++q){
                int pq = INDEX((p + soOffset), (q + soOffset));
                for(int i = 0; i < frzcpi_[h]; ++i)
                    aFzcD[pq] += pCa[p][i] * pCa[q][i];
                for(int i = 0; i < clsdpi_[h] + openpi_[h]; ++i)
                    aD[pq] += pCa[p][i] * pCa[q][i];
                if(transformationType_ != Restricted){
                    for(int i = 0; i < frzcpi_[h]; ++i)
                        bFzcD[pq] += pCb[p][i] * pCb[q][i];
                    for(int i = 0; i < clsdpi_[h]; ++i)
                        bD[pq] += pCb[p][i] * pCb[q][i];
                }
            }
        }
        soOffset += sopi_[h];
    }

    double *T = init_array(nTriSo_);
    if(print_>4) fprintf(outfile, "The SO basis kinetic energy integrals\n");
    IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_T,   T, nTriSo_, 0, print_ > 4, outfile);
    if(print_>4) fprintf(outfile, "The SO basis nuclear attraction integrals\n");
    IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_V, aoH, nTriSo_, 0, print_ > 4, outfile);

    for(int pq=0; pq < nTriSo_; ++pq){
        aoH[pq] += T[pq];
        aFzcOp[pq] = aoH[pq];
        aFock[pq]  = aoH[pq];
        if(transformationType_ != Restricted){
            bFock[pq]  = aoH[pq];
            bFzcOp[pq] = aoH[pq];
        }
    }
    free(T);

    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    int soIntFile = PSIF_SO_TEI;

    IWL *iwl = new IWL(psio_.get(), soIntFile, tolerance_, 1, 1);
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

    double *moInts = init_array(nTriMo_);
    int *order = init_int_array(nmo_);
    // We want to keep Pitzer ordering, so this is just an identity mapping
    for(int n = 0; n < nmo_; ++n) order[n] = n;
    if(print_)
        fprintf(outfile, "\tTransforming the one-electron integrals and constructing Fock matrices\n");
    if(transformationType_ == Restricted){
        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aoH, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The MO basis one-electron integrals\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_OEI, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aFzcOp, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The MO basis frozen core operator\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_FZC, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aFock, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The MO basis Fock operator\n");
            print_array(moInts, nmo_, outfile);
        }

        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_FOCK, nTriMo_, aFock);
    }else{
        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aoH, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The MO basis alpha one-electron integrals\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_A_OEI, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCb = Cb_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aoH, moInts, pCb, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The MO basis beta one-electron integrals\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_B_OEI, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aFzcOp, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The MO basis alpha frozen core operator\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_A_FZC, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCb = Cb_->pointer(h);
            trans_one(sopi_[h], mopi_[h], bFzcOp, moInts, pCb, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The MO basis beta frozen core operator\n");
            print_array(moInts, nmo_, outfile);
        }
        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCa = Ca_->pointer(h);
            trans_one(sopi_[h], mopi_[h], aFock, moInts, pCa, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The MO basis alpha Fock operator\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_A_FOCK, nTriMo_, moInts);

        for(int h = 0, moOffset = 0, soOffset = 0; h < nirreps_; ++h){
            double **pCb = Cb_->pointer(h);
            trans_one(sopi_[h], mopi_[h], bFock, moInts, pCb, soOffset, &(order[moOffset]));
            soOffset += sopi_[h];
            moOffset += mopi_[h];
        }
        if(print_>4){
            fprintf(outfile, "The MO basis beta Fock operator\n");
            print_array(moInts, nmo_, outfile);
        }
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_B_FOCK, nTriMo_, moInts);
    }
    free(order);
    free(moInts);
    free(aFzcD);
    free(aFzcOp);
    free(aD);
    free(aoH);
    free(aFock);
    if(transformationType_ != Restricted){
        free(bFzcD);
        free(bFzcOp);
        free(bD);
        free(bFock);
    }

    dpd_set_default(currentActiveDPD);
#endif
}
