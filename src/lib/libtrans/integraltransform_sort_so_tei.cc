#include "integraltransform.h"
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.hpp>
#include <libmints/matrix.h>
#include "psifiles.h"
#include "integraltransform_functors.h"
#include "mospace.h"
#define EXTERN
#include <libdpd/dpd.gbl>

using namespace psi;

/**
 * Presort the two-electron integrals into DPD buffers to prepare them for
 * the transformation.  The frozen core operator is built simultaneously.
 * If this action has already been performed, it will just load the frozen
 * core operator from disk and return.
 */
void
IntegralTransform::presort_so_tei()
{
    check_initialized();

    if(alreadyPresorted_){
        if(print_>5)
            fprintf(outfile, "\tSO integrals are already sorted, moving on...\n");
            return;
    }

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

    if(print_){
        fprintf(outfile, "\tPresorting SO-basis two-electron integrals.\n");
        fflush(outfile);
    }

    int soIntFile = PSIF_SO_TEI;

    dpdfile4 I;
    psio_->open(PSIF_SO_PRESORT, PSIO_OPEN_NEW);
    dpd_file4_init(&I, PSIF_SO_PRESORT, 0, DPD_ID("[n>=n]+"), DPD_ID("[n>=n]+"), "SO Ints (nn|nn)");

    size_t memoryd = memory_ / sizeof(double);

    int nump = 0, numq = 0;
    for(int h=0; h < nirreps_; ++h){
        nump += I.params->ppi[h];
        numq += I.params->qpi[h];
    }
    int **bucketMap = init_int_matrix(nump, numq);

    /* Room for one bucket to begin with */
    int **bucketOffset = (int **) malloc(sizeof(int *));
    bucketOffset[0] = init_int_array(nirreps_);
    int **bucketRowDim = (int **) malloc(sizeof(int *));
    bucketRowDim[0] = init_int_array(nirreps_);
    long int **bucketSize = (long int **) malloc(sizeof(long int *));
    bucketSize[0] = init_long_int_array(nirreps_);

    /* Figure out how many passes we need and where each p,q goes */
    int nBuckets = 1;
    size_t coreLeft = memoryd;
    psio_address next;
    for(int h = 0; h < nirreps_; ++h){
        size_t rowLength = (size_t) I.params->coltot[h^(I.my_irrep)];
        for(int row=0; row < I.params->rowtot[h]; ++row) {
            if(coreLeft >= rowLength){
                coreLeft -= rowLength;
                bucketRowDim[nBuckets-1][h]++;
                bucketSize[nBuckets-1][h] += rowLength;
            } else {
                nBuckets++;
                coreLeft = memoryd - rowLength;
                /* Make room for another bucket */
                bucketOffset = (int **) realloc((void *) bucketOffset,
                                             nBuckets * sizeof(int *));
                bucketOffset[nBuckets-1] = init_int_array(nirreps_);
                bucketOffset[nBuckets-1][h] = row;

                bucketRowDim = (int **) realloc((void *) bucketRowDim,
                                             nBuckets * sizeof(int *));
                bucketRowDim[nBuckets-1] = init_int_array(nirreps_);
                bucketRowDim[nBuckets-1][h] = 1;

                bucketSize = (long int **) realloc((void *) bucketSize,
                                                nBuckets * sizeof(long int *));
                bucketSize[nBuckets-1] = init_long_int_array(nirreps_);
                bucketSize[nBuckets-1][h] = rowLength;
            }
            int p = I.params->roworb[h][row][0];
            int q = I.params->roworb[h][row][1];
            bucketMap[p][q] = nBuckets - 1;
        }
    }

    if(print_) {
        fprintf(outfile, "\tSorting File: %s nbuckets = %d\n", I.label, nBuckets);
        fflush(outfile);
    }

    next = PSIO_ZERO;
    for(int n=0; n < nBuckets; ++n) { /* nbuckets = number of passes */
        /* Prepare target matrix */
        for(int h=0; h < nirreps_; h++) {
            I.matrix[h] = block_matrix(bucketRowDim[n][h], I.params->coltot[h]);
        }

        DPDFillerFunctor dpdfiller(&I,n,bucketMap,bucketOffset, false, true);
        NullFunctor null;
        IWL *iwl = new IWL(psio_.get(), PSIF_SO_TEI, tolerance_, 1, 1);
        // In the functors below, we only want to build the Fock matrix on the first pass
        if(transformationType_ == Restricted){
            FrozenCoreAndFockRestrictedFunctor fock(aD, aFzcD,aFock,aFzcOp);
            if(n)
                iwl_integrals(iwl, dpdfiller, null);
            else
                iwl_integrals(iwl, dpdfiller, fock);
        }else{
            FrozenCoreAndFockUnrestrictedFunctor fock(aD, bD, aFzcD, bFzcD,
                                                      aFock, bFock, aFzcOp, bFzcOp);
            if(n)
                iwl_integrals(iwl, dpdfiller, null);
            else
                iwl_integrals(iwl, dpdfiller, fock);
        }
        delete iwl;


        for(int h=0; h < nirreps_; ++h) {
            if(bucketSize[n][h])
                psio_->write(I.filenum, I.label, (char *) I.matrix[h][0],
                bucketSize[n][h]*((long int) sizeof(double)), next, &next);
            free_block(I.matrix[h]);
        }
    } /* end loop over buckets/passes */

    /* Get rid of the input integral file */
    psio_->open(soIntFile, PSIO_OPEN_OLD);
    psio_->close(soIntFile, keepIwlSoInts_);

    free_int_matrix(bucketMap);

    for(int n=0; n < nBuckets; ++n) {
        free(bucketOffset[n]);
        free(bucketRowDim[n]);
        free(bucketSize[n]);
    }
    free(bucketOffset);
    free(bucketRowDim);
    free(bucketSize);

    double *moInts = init_array(nTriMo_);
    int *order = init_int_array(nmo_);
    // We want to keep Pitzer ordering, so this is just an identity mapping
    for(int n = 0; n < nmo_; ++n) order[n] = n;
    if(print_)
        fprintf(outfile, "\tTransforming the one-electron integrals and constructing Fock matrices\n");
    if(transformationType_ == Restricted){

        // Compute frozen core energy
        size_t pq = 0;
        frozen_core_energy_ = 0.0;

        for(int p=0; p < nso_; p++) {
            for(int q=0; q <= p ; q++, pq++) {
                double prefact = p == q ? 1.0 : 2.0;
                frozen_core_energy_ += prefact * aFzcD[pq] * (aoH[pq] + aFzcOp[pq]);
            }
        }

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

        // Compute frozen-core energy
        size_t pq = 0;
        frozen_core_energy_ = 0.0;
        for(int p=0; p < nso_; p++) {
            for(int q=0; q <= p; q++, pq++) {
                double prefact = p == q ? 0.5 : 1.0;
                frozen_core_energy_ += prefact * aFzcD[pq] * (aoH[pq] + aFzcOp[pq]);
                frozen_core_energy_ += prefact * bFzcD[pq] * (aoH[pq] + bFzcOp[pq]);
            }
        }

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
        IWL::write_one(psio_.get(), PSIF_OEI, PSIF_MO_B_FZC, nTriMo_, moInts);
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

    alreadyPresorted_ = true;

    dpd_file4_close(&I);
    psio_->close(PSIF_SO_PRESORT, 1);
}


