#include "integraltransform.h"
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <libmints/matrix.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include "psifiles.h"
#include "ccfiles.h"
#include "mospace.h"
#define EXTERN
#include <libdpd/dpd.gbl>

using namespace psi;
using namespace boost;

void
IntegralTransform::transform_tei_second_half(const shared_ptr<MOSpace> s1, const shared_ptr<MOSpace> s2,
                                             const shared_ptr<MOSpace> s3, const shared_ptr<MOSpace> s4)
{
    check_initialized();

    char *label = new char[100];

    // Grab the transformation coefficients
    SharedMatrix c3a = aMOCoefficients_[s3->label()];
    SharedMatrix c3b = bMOCoefficients_[s3->label()];
    SharedMatrix c4a = aMOCoefficients_[s4->label()];
    SharedMatrix c4b = bMOCoefficients_[s4->label()];
    // And the number of orbitals per irrep
    int *aOrbsPI3 = aOrbsPI_[s3->label()];
    int *bOrbsPI3 = bOrbsPI_[s3->label()];
    int *aOrbsPI4 = aOrbsPI_[s4->label()];
    int *bOrbsPI4 = bOrbsPI_[s4->label()];
    // The reindexing arrays
    int *aIndex1 = aIndices_[s1->label()];
    int *bIndex1 = bIndices_[s1->label()];
    int *aIndex2 = aIndices_[s2->label()];
    int *bIndex2 = bIndices_[s2->label()];
    int *aIndex3 = aIndices_[s3->label()];
    int *bIndex3 = bIndices_[s3->label()];
    int *aIndex4 = aIndices_[s4->label()];
    int *bIndex4 = bIndices_[s4->label()];

    // Grab control of DPD for now, but store the active number to restore it later
    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    IWL *iwl;
    if(useIWL_) iwl = new IWL;
    int nBuckets;
    int thisBucketRows;
    size_t rowsPerBucket;
    size_t rowsLeft;
    size_t memFree;
    dpdbuf4 J, K;

    double **TMP = block_matrix(nso_, nso_);

    if(print_) {
        if(transformationType_ == Restricted){
            fprintf(outfile, "\tStarting second half-transformation.\n");
        }else{
            fprintf(outfile, "\tStarting AA second half-transformation.\n");
        }
        fflush(outfile);
    }

    if(useIWL_) iwl = new IWL(psio_.get(), iwlAAIntFile_, tolerance_, 0, 0);

    psio_->open(dpdIntFile_, PSIO_OPEN_OLD);
    psio_->open(aHtIntFile_, PSIO_OPEN_OLD);

    int braCore = DPD_ID(s1, s2, Alpha, true);
    int braDisk = DPD_ID(s1, s2, Alpha, true);
    int ketCore = 0;
    int ketDisk = 3;
    sprintf(label, "Half-Transformed Ints (%c%c|nn)", toupper(s1->label()), toupper(s2->label()));
    dpd_buf4_init(&J, aHtIntFile_, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
    if(print_ > 5)
        fprintf(outfile, "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
                            label, braCore, ketCore, braDisk, ketDisk);

    braCore = DPD_ID(s1, s2, Alpha, true);
    ketCore = DPD_ID(s3, s4, Alpha, false);
    braDisk = DPD_ID(s1, s2, Alpha, true);
    ketDisk = DPD_ID(s3, s4, Alpha, true);
    sprintf(label, "MO Ints (%c%c|%c%c)", toupper(s1->label()), toupper(s2->label()),
                                          toupper(s3->label()), toupper(s4->label()));
    dpd_buf4_init(&K, dpdIntFile_, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
    if(print_ > 5)
        fprintf(outfile, "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
                            label, braCore, ketCore, braDisk, ketDisk);

    for(int h=0; h < nirreps_; h++) {
        if(J.params->coltot[h] && J.params->rowtot[h]) {
            memFree = static_cast<size_t>(dpd_memfree() - J.params->coltot[h] - K.params->coltot[h]);
            rowsPerBucket = memFree/(2 * J.params->coltot[h]);
            if(rowsPerBucket > J.params->rowtot[h])
                rowsPerBucket = static_cast<size_t>(J.params->rowtot[h]);
            nBuckets = static_cast<int>(ceil(static_cast<double>(J.params->rowtot[h])/
                                        static_cast<double>(rowsPerBucket)));
            rowsLeft = static_cast<size_t>(J.params->rowtot[h] % rowsPerBucket);
        }else {
            nBuckets = 0;
            rowsPerBucket = 0;
            rowsLeft = 0;
        }

        if(print_ > 1) {
            fprintf(outfile, "\th = %d; memfree         = %lu\n", h, memFree);
            fprintf(outfile, "\th = %d; rows_per_bucket = %lu\n", h, rowsPerBucket);
            fprintf(outfile, "\th = %d; rows_left       = %lu\n", h, rowsLeft);
            fprintf(outfile, "\th = %d; nbuckets        = %d\n", h, nBuckets);
            fflush(outfile);
        }

        dpd_buf4_mat_irrep_init_block(&J, h, rowsPerBucket);
        dpd_buf4_mat_irrep_init_block(&K, h, rowsPerBucket);

        for(int n=0; n < nBuckets; n++) {
            if(nBuckets == 1)
                thisBucketRows = rowsPerBucket;
            else
                thisBucketRows = (n < nBuckets-1) ? rowsPerBucket : rowsLeft;
            dpd_buf4_mat_irrep_rd_block(&J, h, n*rowsPerBucket, thisBucketRows);
            for(int pq=0; pq < thisBucketRows; pq++) {
                for(int Gr=0; Gr < nirreps_; Gr++) {
                    // Transform ( S1 S2 | n n ) -> ( S1 S2 | n S4 )
                    int Gs = h^Gr;
                    int nrows = sopi_[Gr];
                    int ncols = aOrbsPI4[Gs];
                    int nlinks = sopi_[Gs];
                    int rs = J.col_offset[h][Gr];
                    double **pc4a = c4a->pointer(Gs);
                    if(nrows && ncols && nlinks)
                        C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, &J.matrix[h][pq][rs],
                                nlinks, pc4a[0], ncols, 0.0, TMP[0], nso_);
                    //TODO else if s4->label() == MOSPACE_NIL, copy buffer...

                    // Transform ( S1 S2 | n S4 ) -> ( S1 S2 | S3 S4 )
                    nrows = aOrbsPI3[Gr];
                    ncols = aOrbsPI4[Gs];
                    nlinks = sopi_[Gr];
                    rs = K.col_offset[h][Gr];
                    double **pc3a = c3a->pointer(Gr);
                    if(nrows && ncols && nlinks)
                        C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, pc3a[0], nrows ,
                                TMP[0], nso_, 0.0, &K.matrix[h][pq][rs], ncols);
                    //TODO else if s3->label() == MOSPACE_NIL, copy buffer...
                } /* Gr */
                if(useIWL_){
                    int p = aIndex1[K.params->roworb[h][pq+n*rowsPerBucket][0]];
                    int q = aIndex2[K.params->roworb[h][pq+n*rowsPerBucket][1]];
                    size_t PQ = INDEX(p,q);
                    for(int rs=0; rs < K.params->coltot[h]; rs++) {
                        int r = aIndex3[K.params->colorb[h][rs][0]];
                        int s = aIndex4[K.params->colorb[h][rs][1]];
                        size_t RS = INDEX(r,s);
                        if(r >= s && RS <= PQ)
                            iwl->write_value(p, q, r, s, K.matrix[h][pq][rs],
                                             printTei_, outfile, 0);
                    } /* rs */
                }
            } /* pq */
            dpd_buf4_mat_irrep_wrt_block(&K, h, n*rowsPerBucket, thisBucketRows);
        }
        dpd_buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
        dpd_buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
    }
    dpd_buf4_close(&K);
    dpd_buf4_close(&J);

    if(useIWL_){
        iwl->flush(1);
        iwl->set_keep_flag(1);
        // This closes the file too
        delete iwl;
    }

    if(transformationType_ != Restricted){
        if(print_) {
            fprintf(outfile, "\tStarting AB second half-transformation.\n");
            fflush(outfile);
        }
        if(useIWL_) iwl = new IWL(psio_.get(), iwlABIntFile_, tolerance_, 0, 0);

        braCore = braDisk = DPD_ID(s1, s2, Alpha, true);
        ketCore = 0;
        ketDisk = 3;
        sprintf(label, "Half-Transformed Ints (%c%c|nn)", toupper(s1->label()), toupper(s2->label()));
        dpd_buf4_init(&J, aHtIntFile_, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(print_ > 5)
            fprintf(outfile, "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
                                label, braCore, ketCore, braDisk, ketDisk);

        braCore = DPD_ID(s1, s2, Alpha, true);
        ketCore = DPD_ID(s3, s4, Beta,  false);
        braDisk = DPD_ID(s1, s2, Alpha, true);
        ketDisk = DPD_ID(s3, s4, Beta,  true);
        sprintf(label, "MO Ints (%c%c|%c%c)", toupper(s1->label()), toupper(s2->label()),
                                              tolower(s3->label()), tolower(s4->label()));
        dpd_buf4_init(&K, dpdIntFile_, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(print_ > 5)
            fprintf(outfile, "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
                                label, braCore, ketCore, braDisk, ketDisk);

        for(int h=0; h < nirreps_; h++) {
            if(J.params->coltot[h] && J.params->rowtot[h]) {
                memFree = static_cast<size_t>(dpd_memfree() - J.params->coltot[h] - K.params->coltot[h]);
                rowsPerBucket = memFree/(2 * J.params->coltot[h]);
                if(rowsPerBucket > J.params->rowtot[h])
                    rowsPerBucket = static_cast<size_t>(J.params->rowtot[h]);
                nBuckets = static_cast<int>(ceil(static_cast<double>(J.params->rowtot[h])/
                        static_cast<double>(rowsPerBucket)));
                rowsLeft = static_cast<size_t>(J.params->rowtot[h] % rowsPerBucket);
            }else{
                nBuckets = 0;
                rowsPerBucket = 0;
                rowsLeft = 0;
            }

            if(print_ > 1) {
                fprintf(outfile, "\th = %d; memfree         = %lu\n", h, memFree);
                fprintf(outfile, "\th = %d; rows_per_bucket = %lu\n", h, rowsPerBucket);
                fprintf(outfile, "\th = %d; rows_left       = %lu\n", h, rowsLeft);
                fprintf(outfile, "\th = %d; nbuckets        = %d\n", h, nBuckets);
                fflush(outfile);
            }

            dpd_buf4_mat_irrep_init_block(&J, h, rowsPerBucket);
            dpd_buf4_mat_irrep_init_block(&K, h, rowsPerBucket);

            for(int n=0; n < nBuckets; n++){
                if(nBuckets == 1)
                    thisBucketRows = rowsPerBucket;
                else
                    thisBucketRows = (n < nBuckets-1) ? rowsPerBucket : rowsLeft;
                dpd_buf4_mat_irrep_rd_block(&J, h, n*rowsPerBucket, thisBucketRows);
                for(int pq=0; pq < thisBucketRows; pq++) {
                    for(int Gr=0; Gr < nirreps_; Gr++) {
                        // Transform ( S1 S2 | n n ) -> ( S1 S2 | n s4 )
                        int Gs = h^Gr;
                        int nrows = sopi_[Gr];
                        int ncols = bOrbsPI4[Gs];
                        int nlinks = sopi_[Gs];
                        int rs = J.col_offset[h][Gr];
                        double **pc4b = c4b->pointer(Gs);
                        if(nrows && ncols && nlinks)
                            C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, &J.matrix[h][pq][rs],
                                    nlinks, pc4b[0], ncols, 0.0, TMP[0], nso_);
                        //TODO else if s4->label() == MOSPACE_NIL, copy buffer...

                        // Transform ( S1 S2 | n s4 ) -> ( S1 S2 | s3 s4 )
                        nrows = bOrbsPI3[Gr];
                        ncols = bOrbsPI4[Gs];
                        nlinks = sopi_[Gr];
                        rs = K.col_offset[h][Gr];
                        double **pc3b = c3b->pointer(Gr);
                        if(nrows && ncols && nlinks)
                            C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, pc3b[0], nrows,
                                    TMP[0], nso_, 0.0, &K.matrix[h][pq][rs], ncols);
                        //TODO else if s3->label() == MOSPACE_NIL, copy buffer...
                    } /* Gr */
                    if(useIWL_){
                        int p = aIndex1[K.params->roworb[h][pq+n*rowsPerBucket][0]];
                        int q = aIndex2[K.params->roworb[h][pq+n*rowsPerBucket][1]];
                        size_t PQ = INDEX(p,q);
                        for(int rs=0; rs < K.params->coltot[h]; rs++) {
                            int r = bIndex3[K.params->colorb[h][rs][0]];
                            int s = bIndex4[K.params->colorb[h][rs][1]];
                            size_t RS = INDEX(r,s);
                            if(r >= s)
                                iwl->write_value(p, q, r, s, K.matrix[h][pq][rs],
                                                 printTei_, outfile, 0);
                        } /* rs */
                    }
                } /* pq */
                dpd_buf4_mat_irrep_wrt_block(&K, h, n*rowsPerBucket, thisBucketRows);
            }
            dpd_buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
            dpd_buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
        }
        dpd_buf4_close(&K);
        dpd_buf4_close(&J);

        if(useIWL_){
            iwl->flush(1);
            iwl->set_keep_flag(1);
            // This closes the file too
            delete iwl;
        }

        /*** AA/AB two-electron integral transformation complete ***/

        if(print_) {
            fprintf(outfile, "\tStarting BB second half-transformation.\n");
            fflush(outfile);
        }
        if(useIWL_) iwl = new IWL(psio_.get(), iwlBBIntFile_, tolerance_, 0, 0);

        psio_->open(bHtIntFile_, PSIO_OPEN_OLD);

        braCore = DPD_ID(s1, s2, Beta, true);
        ketCore = 0;
        braDisk = DPD_ID(s1, s2, Beta, true);
        ketDisk = 3;
        sprintf(label, "Half-Transformed Ints (%c%c|nn)", tolower(s1->label()), tolower(s2->label()));
        dpd_buf4_init(&J, bHtIntFile_, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(print_ > 5)
            fprintf(outfile, "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
                                label, braCore, ketCore, braDisk, ketDisk);

        braCore = DPD_ID(s1, s2, Beta, true);
        ketCore = DPD_ID(s3, s4, Beta, false);
        braDisk = DPD_ID(s1, s2, Beta, true);
        ketDisk = DPD_ID(s3, s4, Beta, true);
        sprintf(label, "MO Ints (%c%c|%c%c)", tolower(s1->label()), tolower(s2->label()),
                                              tolower(s3->label()), tolower(s4->label()));
        dpd_buf4_init(&K, dpdIntFile_, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(print_ > 5)
            fprintf(outfile, "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
                                label, braCore, ketCore, braDisk, ketDisk);

        for(int h=0; h < nirreps_; h++) {
            if (J.params->coltot[h] && J.params->rowtot[h]) {
                memFree = static_cast<size_t>(dpd_memfree() - J.params->coltot[h] - K.params->coltot[h]);
                rowsPerBucket = memFree/(2 * J.params->coltot[h]);
                if(rowsPerBucket > J.params->rowtot[h])
                    rowsPerBucket = static_cast<size_t>(J.params->rowtot[h]);
                nBuckets = static_cast<int>(ceil(static_cast<double>(J.params->rowtot[h])/
                                static_cast<double>(rowsPerBucket)));
                rowsLeft = static_cast<size_t>(J.params->rowtot[h] % rowsPerBucket);
            }
            else {
                nBuckets = 0;
                rowsPerBucket = 0;
                rowsLeft = 0;
            }

            if(print_ > 1) {
                fprintf(outfile, "\th = %d; memfree         = %lu\n", h, memFree);
                fprintf(outfile, "\th = %d; rows_per_bucket = %lu\n", h, rowsPerBucket);
                fprintf(outfile, "\th = %d; rows_left       = %lu\n", h, rowsLeft);
                fprintf(outfile, "\th = %d; nbuckets        = %d\n", h, nBuckets);
                fflush(outfile);
            }

            dpd_buf4_mat_irrep_init_block(&J, h, rowsPerBucket);
            dpd_buf4_mat_irrep_init_block(&K, h, rowsPerBucket);

            for(int n=0; n < nBuckets; n++) {
                if(nBuckets == 1)
                    thisBucketRows = rowsPerBucket;
                else
                    thisBucketRows = (n < nBuckets-1) ? rowsPerBucket : rowsLeft;
                dpd_buf4_mat_irrep_rd_block(&J, h, n*rowsPerBucket, thisBucketRows);
                for(int pq=0; pq < thisBucketRows; pq++) {
                    for(int Gr=0; Gr < nirreps_; Gr++) {
                        // Transform ( s1 s2 | n n ) -> ( s1 s2 | n s4 )
                        int Gs = h^Gr;
                        int nrows = sopi_[Gr];
                        int ncols = bOrbsPI4[Gs];
                        int nlinks = sopi_[Gs];
                        int rs = J.col_offset[h][Gr];
                        double **pc4b = c4b->pointer(Gs);
                        if(nrows && ncols && nlinks)
                            C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, &J.matrix[h][pq][rs],
                                    nlinks, pc4b[0], ncols, 0.0, TMP[0], nso_);

                        // Transform ( s1 s2 | n s4 ) -> ( s1 s2 | s3 s4 )
                        nrows = bOrbsPI3[Gr];
                        ncols = bOrbsPI4[Gs];
                        nlinks = sopi_[Gr];
                        rs = K.col_offset[h][Gr];
                        double **pc3b = c3b->pointer(Gr);
                        if(nrows && ncols && nlinks)
                            C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, pc3b[0], nrows,
                                    TMP[0], nso_, 0.0, &K.matrix[h][pq][rs], ncols);
                    } /* Gr */
                    if(useIWL_){
                        int p = bIndex1[K.params->roworb[h][pq+n*rowsPerBucket][0]];
                        int q = bIndex2[K.params->roworb[h][pq+n*rowsPerBucket][1]];
                        size_t PQ = INDEX(p,q);
                        for(int rs=0; rs < K.params->coltot[h]; rs++) {
                            int r = bIndex3[K.params->colorb[h][rs][0]];
                            int s = bIndex4[K.params->colorb[h][rs][1]];
                            size_t RS = INDEX(r,s);
                            if(r >= s && RS <= PQ)
                                iwl->write_value(p, q, r, s, K.matrix[h][pq][rs],
                                                 printTei_, outfile, 0);
                        } /* rs */
                    }
                } /* pq */
                dpd_buf4_mat_irrep_wrt_block(&K, h, n*rowsPerBucket, thisBucketRows);
            }
            dpd_buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
            dpd_buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
        }
        dpd_buf4_close(&K);
        dpd_buf4_close(&J);

        psio_->close(bHtIntFile_, keepHtInts_);

        if(useIWL_){
            iwl->flush(1);
            iwl->set_keep_flag(1);
            // This closes the file too
            delete iwl;
        }
        /*** BB two-electron integral transformation complete ***/
    } // End "if not restricted transformation"


    psio_->close(dpdIntFile_, 1);
    psio_->close(aHtIntFile_, keepHtInts_);

    free_block(TMP);
    delete [] label;

    if(print_){
        fprintf(outfile, "\tTwo-electron integral transformation complete.\n");
        fflush(outfile);
    }

    // Hand DPD control back to the user
    dpd_set_default(currentActiveDPD);
}
