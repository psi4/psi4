#include "integraltransform.h"
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include "psifiles.h"
#include "ccfiles.h"
#include "mospace.h"
#include "spaceinfo.h"

namespace psi{ namespace libtrans{

void
IntegralTransform::transform_tei(shared_ptr<MOSpace> s1, shared_ptr<MOSpace> s2,
                                 shared_ptr<MOSpace> s3, shared_ptr<MOSpace> s4)
{
    // This can be safely called - it returns immediately if the SO ints are already sorted
    presort_so_tei();

    char *label = new char[100];

    // Grab the transformation coefficients
    double ***c1a = _aMOCoefficients[s1->label()];
    double ***c1b = _bMOCoefficients[s1->label()];
    double ***c2a = _aMOCoefficients[s2->label()];
    double ***c2b = _bMOCoefficients[s2->label()];
    double ***c3a = _aMOCoefficients[s3->label()];
    double ***c3b = _bMOCoefficients[s3->label()];
    double ***c4a = _aMOCoefficients[s4->label()];
    double ***c4b = _bMOCoefficients[s4->label()];
    // And the number of orbitals per irrep
    int *aOrbsPI1 = _aOrbsPI[s1->label()];
    int *bOrbsPI1 = _bOrbsPI[s1->label()];
    int *aOrbsPI2 = _aOrbsPI[s2->label()];
    int *bOrbsPI2 = _bOrbsPI[s2->label()];
    int *aOrbsPI3 = _aOrbsPI[s3->label()];
    int *bOrbsPI3 = _bOrbsPI[s3->label()];
    int *aOrbsPI4 = _aOrbsPI[s4->label()];
    int *bOrbsPI4 = _bOrbsPI[s4->label()];
    // The reindexing arrays
    int *aIndex1 = _aIndices[s1->label()];
    int *bIndex1 = _bIndices[s1->label()];
    int *aIndex2 = _aIndices[s2->label()];
    int *bIndex2 = _bIndices[s2->label()];
    int *aIndex3 = _aIndices[s3->label()];
    int *bIndex3 = _bIndices[s3->label()];
    int *aIndex4 = _aIndices[s4->label()];
    int *bIndex4 = _bIndices[s4->label()];

    struct iwlbuf MBuff;
    int nBuckets;
    int thisBucketRows;
    size_t rowsPerBucket;
    size_t rowsLeft;
    size_t memFree;
    bool doIWL = _outputType == IWLAndDPD;
    
    double **TMP = block_matrix(_nso, _nso);

    /*** AA/AB two-electron integral transformation ***/

    if(_print) {
        if(_transformationType == Restricted){
            fprintf(outfile, "\tStarting first half-transformation.\n");
        }else{
            fprintf(outfile, "\tStarting AA/AB first half-transformation.\n");
        }
        fflush(outfile);
    }

    psio_open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);
    psio_open(PSIF_HALFT0, PSIO_OPEN_NEW);
    dpdbuf4 J, K;
    dpd_buf4_init(&J, PSIF_SO_PRESORT, 0, 3, 0, 3, 3, 0, "SO Ints (nn,nn)");

    int braCore = 3;
    int ketCore = DPD_ID(s1, s2, Alpha, false);
    int braDisk = 3;
    int ketDisk = DPD_ID(s1, s2, Alpha, true);
    sprintf(label, "Half-Transformed Ints (nn,%c%c)", toupper(s1->label()), toupper(s2->label()));
    dpd_buf4_init(&K, PSIF_HALFT0, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
    if(_print > 5)
        fprintf(outfile, "Initializing %s, \nin core:(%d|%d) on disk(%d|%d)",
                            label, braCore, ketCore, braDisk, ketDisk);
    
    for(int h=0; h < _nirreps; h++) {
        if(J.params->coltot[h]) {
            memFree = static_cast<size_t>(dpd_memfree() - J.params->coltot[h] - K.params->coltot[h]);
            rowsPerBucket = memFree/(2 * J.params->coltot[h]);
            if(rowsPerBucket > J.params->rowtot[h]) rowsPerBucket = (size_t) J.params->rowtot[h];
            nBuckets = static_cast<int>(ceil(static_cast<double>(J.params->rowtot[h])/
                                        static_cast<double>(rowsPerBucket)));
            rowsLeft = static_cast<size_t>(J.params->rowtot[h] % rowsPerBucket);
        }else{
            nBuckets = 0;
            rowsPerBucket = 0;
            rowsLeft = 0;
        }

        if(_print > 1) {
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
                for(int Gr=0; Gr < _nirreps; Gr++) {
                    // Transform ( n n | n n ) -> ( n n | n S2 )
                    int Gs = h^Gr;
                    int nrows = _sopi[Gr];
                    int ncols = aOrbsPI2[Gs];
                    int nlinks = _sopi[Gs];
                    int rs = J.col_offset[h][Gr];
                    if(nrows && ncols && nlinks)
                        C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, &J.matrix[h][pq][rs],
                                nlinks, c2a[Gs][0], ncols, 0.0, TMP[0], _nso);

                    // Transform ( n n | n S2 ) -> ( n n | S1 S2 )
                    nrows = aOrbsPI1[Gr];
                    ncols = aOrbsPI2[Gs];
                    nlinks = _sopi[Gr];
                    rs = K.col_offset[h][Gr];
                    if(nrows && ncols && nlinks)
                        C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, c1a[Gr][0], nrows,
                                TMP[0], _nso, 0.0, &K.matrix[h][pq][rs], ncols);
                } /* Gr */
            } /* pq */
            dpd_buf4_mat_irrep_wrt_block(&K, h, n*rowsPerBucket, thisBucketRows);
        }
        dpd_buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
        dpd_buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
    }
    dpd_buf4_close(&K);
    dpd_buf4_close(&J);
        

    psio_close(PSIF_SO_PRESORT, 1); /* must keep the presort file for the upcoming BB transformation */

    if(_print) {
        if(_transformationType == Restricted){
            fprintf(outfile, "\tSorting half-transformed integrals.\n");
        }else{
            fprintf(outfile, "\tSorting AA/AB half-transformed integrals.\n");
        }
        fflush(outfile);
    }


    psio_open(PSIF_HALFT1, PSIO_OPEN_NEW);

    braCore = braDisk = 3;
    ketCore = braDisk = DPD_ID(s1, s2, Alpha, true);
    sprintf(label, "Half-Transformed Ints (nn,%c%c)", toupper(s1->label()), toupper(s2->label()));
    dpd_buf4_init(&K, PSIF_HALFT0, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
    if(_print > 5)
        fprintf(outfile, "Initializing %s, \nin core:(%d|%d) on disk(%d|%d)",
                            label, braCore, ketCore, braDisk, ketDisk);
    sprintf(label, "Half-Transformed Ints (%c%c,nn)", toupper(s1->label()), toupper(s2->label()));
    dpd_buf4_sort(&K, PSIF_HALFT1, rspq, ketCore, braCore, label);
    dpd_buf4_close(&K);

    psio_close(PSIF_HALFT0, 0);

    if(_print) {
        if(_transformationType == Restricted){
            fprintf(outfile, "\tStarting second half-transformation.\n");
        }else{
            fprintf(outfile, "\tStarting AA second half-transformation.\n");
        }
        fflush(outfile);
    }
    // TODO check this..
    iwl_buf_init(&MBuff, PSIF_MO_AA_TEI, _tolerance, 0, 0);

    braCore = braDisk = DPD_ID(s1, s2, Alpha, true);
    ketCore = 0;
    ketDisk = 3;
    sprintf(label, "Half-Transformed Ints (%c%c,nn)", toupper(s1->label()), toupper(s2->label()));
    dpd_buf4_init(&J, PSIF_HALFT1, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
    if(_print > 5)
        fprintf(outfile, "Initializing %s, \nin core:(%d|%d) on disk(%d|%d)",
                            label, braCore, ketCore, braDisk, ketDisk);

    braCore = DPD_ID(s1, s2, Alpha, true);
    ketCore = DPD_ID(s3, s4, Alpha, false);
    braDisk = DPD_ID(s1, s2, Alpha, true);
    ketDisk = DPD_ID(s3, s4, Alpha, true);
    sprintf(label, "MO Ints (%c%c,%c%c)", toupper(s1->label()), toupper(s2->label()),
                                          toupper(s3->label()), toupper(s4->label()));
    dpd_buf4_init(&K, CC_MISC, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
    if(_print > 5)
        fprintf(outfile, "Initializing %s, \nin core:(%d|%d) on disk(%d|%d)",
                            label, braCore, ketCore, braDisk, ketDisk);

    for(int h=0; h < _nirreps; h++) {
        if(J.params->coltot[h]) {
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

        if(_print > 1) {
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
                for(int Gr=0; Gr < _nirreps; Gr++) {
                    // Transform ( S1 S2 | n n ) -> ( S1 S2 | n S4 )
                    int Gs = h^Gr;
                    int nrows = _sopi[Gr];
                    int ncols = aOrbsPI4[Gs];
                    int nlinks = _sopi[Gs];
                    int rs = J.col_offset[h][Gr];
                    if(nrows && ncols && nlinks)
                        C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, &J.matrix[h][pq][rs],
                                nlinks, c4a[Gs][0], ncols, 0.0, TMP[0], _nso);

                    // Transform ( S1 S2 | n S4 ) -> ( S1 S2 | S3 S4 )
                    nrows = aOrbsPI3[Gr];
                    ncols = aOrbsPI4[Gs];
                    nlinks = _sopi[Gr];
                    rs = K.col_offset[h][Gr];
                    if(nrows && ncols && nlinks)
                        C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, c3a[Gr][0], nrows ,
                                TMP[0], _nso, 0.0, &K.matrix[h][pq][rs], ncols);
                } /* Gr */
                int p = aIndex1[K.params->roworb[h][pq+n*rowsPerBucket][0]];
                int q = aIndex2[K.params->roworb[h][pq+n*rowsPerBucket][1]];
                size_t PQ = INDEX(p,q);
                for(int rs=0; rs < K.params->coltot[h]; rs++) {
                    int r = aIndex3[K.params->colorb[h][rs][0]];
                    int s = aIndex4[K.params->colorb[h][rs][1]];
                    size_t RS = INDEX(r,s);
                    if(r >= s && RS <= PQ)
                        iwl_buf_wrt_val(&MBuff, p, q, r, s, K.matrix[h][pq][rs],
                                            _printTei, outfile, 0);
                } /* rs */
            } /* pq */
        }
        dpd_buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
        dpd_buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
    }
    dpd_buf4_close(&K);
    dpd_buf4_close(&J);

    iwl_buf_flush(&MBuff, 1);
    iwl_buf_close(&MBuff, 1);

    if(_transformationType == Restricted){
        // We don't need to do a second half transformation
        psio_close(PSIF_HALFT1, 0);
    }else{
        if(_print) {
            fprintf(outfile, "\tStarting AB second half-transformation.\n");
            fflush(outfile);
        }
        iwl_buf_init(&MBuff, PSIF_MO_AB_TEI, _tolerance, 0, 0);

        braCore = braDisk = DPD_ID(s1, s2, Alpha, true);
        ketCore = 0;
        ketDisk = 3;
        sprintf(label, "Half-Transformed Ints (%c%c,nn)", toupper(s1->label()), toupper(s2->label()));
        dpd_buf4_init(&J, PSIF_HALFT1, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(_print > 5)
            fprintf(outfile, "Initializing %s, \nin core:(%d|%d) on disk(%d|%d)",
                                label, braCore, ketCore, braDisk, ketDisk);

        braCore = DPD_ID(s1, s2, Alpha, true);
        ketCore = DPD_ID(s3, s4, Beta,  false);
        braDisk = DPD_ID(s1, s2, Alpha, true);
        ketDisk = DPD_ID(s3, s4, Beta,  true);
        sprintf(label, "MO Ints (%c%c,%c%c)", toupper(s1->label()), toupper(s2->label()),
                                              tolower(s3->label()), tolower(s4->label()));
        dpd_buf4_init(&K, CC_MISC, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(_print > 5)
            fprintf(outfile, "Initializing %s, \nin core:(%d|%d) on disk(%d|%d)",
                                label, braCore, ketCore, braDisk, ketDisk);
        
        for(int h=0; h < _nirreps; h++) {
            if(J.params->coltot[h]){
                static_cast<size_t>(dpd_memfree() - J.params->coltot[h] - K.params->coltot[h]);
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

            if(_print > 1) {
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
                    for(int Gr=0; Gr < _nirreps; Gr++) {
                        // Transform ( S1 S2 | n n ) -> ( S1 S2 | n s4 )
                        int Gs = h^Gr;
                        int nrows = _sopi[Gr];
                        int ncols = bOrbsPI4[Gs];
                        int nlinks = _sopi[Gs];
                        int rs = J.col_offset[h][Gr];
                        if(nrows && ncols && nlinks)
                            C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, &J.matrix[h][pq][rs],
                                    nlinks, c4b[Gs][0], ncols, 0.0, TMP[0], _nso);

                        // Transform ( S1 S2 | n s4 ) -> ( S1 S2 | s3 s4 )
                        nrows = bOrbsPI3[Gr];
                        ncols = bOrbsPI4[Gs];
                        nlinks = _sopi[Gr];
                        rs = K.col_offset[h][Gr];
                        if(nrows && ncols && nlinks)
                            C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, c3b[Gr][0], nrows,
                                    TMP[0], _nso, 0.0, &K.matrix[h][pq][rs], ncols);
                    } /* Gr */
                    int p = aIndex1[K.params->roworb[h][pq+n*rowsPerBucket][0]];
                    int q = aIndex2[K.params->roworb[h][pq+n*rowsPerBucket][1]];
                    size_t PQ = INDEX(p,q);
                    for(int rs=0; rs < K.params->coltot[h]; rs++) {
                        int r = bIndex3[K.params->colorb[h][rs][0]];
                        int s = bIndex4[K.params->colorb[h][rs][1]];
                        size_t RS = INDEX(r,s);
                        if(r >= s)
                            iwl_buf_wrt_val(&MBuff, p, q, r, s, K.matrix[h][pq][rs],
                                             _printTei, outfile, 0);
                    } /* rs */
                } /* pq */
            }
            dpd_buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
            dpd_buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
        }
        dpd_buf4_close(&K);
        dpd_buf4_close(&J);

        iwl_buf_flush(&MBuff, 1);
        iwl_buf_close(&MBuff, 1);

        psio_close(PSIF_HALFT1, 0);

        /*** AA/AB two-electron integral transformation complete ***/

        /*** BB two-electron integral transformation ***/

        if(_print) {
            fprintf(outfile, "\tStarting BB first half-transformation.\n");
            fflush(outfile);
        }

        psio_open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);
        psio_open(PSIF_HALFT0, PSIO_OPEN_NEW);

        dpd_buf4_init(&J, PSIF_SO_PRESORT, 0, 3, 0, 3, 3, 0, "SO Ints (nn,nn)");

        braCore = 3;
        ketCore = DPD_ID(s1, s2, Beta, false);
        braDisk = 3;
        ketDisk = DPD_ID(s1, s2, Beta, true);
        sprintf(label, "Half-Transformed Ints (nn,%c%c)", tolower(s1->label()), tolower(s2->label()));
        dpd_buf4_init(&K, PSIF_HALFT0, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(_print > 5)
            fprintf(outfile, "Initializing %s, \nin core:(%d|%d) on disk(%d|%d)",
                                label, braCore, ketCore, braDisk, ketDisk);

        for(int h=0; h < _nirreps; h++) {
            if(J.params->coltot[h]){
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

            if(_print > 1){
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
                    for(int Gr=0; Gr < _nirreps; Gr++) {
                        // Transform ( n n | n n ) -> ( n n | n s2 )
                        int Gs = h^Gr;
                        int nrows = _sopi[Gr];
                        int ncols = bOrbsPI2[Gs];
                        int nlinks = _sopi[Gs];
                        int rs = J.col_offset[h][Gr];
                        if(nrows && ncols && nlinks)
                            C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, &J.matrix[h][pq][rs],
                            nlinks, c2b[Gs][0], ncols, 0.0, TMP[0], _nso);

                        // Transform ( n n | n s2 ) -> ( n n | s1 s2 )
                        nrows = bOrbsPI1[Gr];
                        ncols = bOrbsPI2[Gs];
                        nlinks = _sopi[Gr];
                        rs = K.col_offset[h][Gr];
                        if(nrows && ncols && nlinks)
                            C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, c1b[Gr][0], nrows,
                                    TMP[0], _nso, 0.0, &K.matrix[h][pq][rs], ncols);
                    } /* Gr */
                } /* pq */
                dpd_buf4_mat_irrep_wrt_block(&K, h, n*rowsPerBucket, thisBucketRows);
            }
            dpd_buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
            dpd_buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
        }
        dpd_buf4_close(&K);
        dpd_buf4_close(&J);

        psio_close(PSIF_SO_PRESORT, 0);

        if(_print) {
            fprintf(outfile, "\tSorting BB half-transformed integrals.\n");
            fflush(outfile);
        }

        psio_open(PSIF_HALFT1, PSIO_OPEN_NEW);

        braCore = 3;
        ketCore = DPD_ID(s1, s2, Beta, true);
        braDisk = 3;
        ketDisk = DPD_ID(s1, s2, Beta, true);
        sprintf(label, "Half-Transformed Ints (nn,%c%c)", tolower(s1->label()), tolower(s2->label()));
        dpd_buf4_init(&K, PSIF_HALFT0, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(_print > 5)
            fprintf(outfile, "Initializing %s, \nin core:(%d|%d) on disk(%d|%d)",
                                label, braCore, ketCore, braDisk, ketDisk);

        sprintf(label, "Half-Transformed Ints (%c%c,nn)", tolower(s1->label()), tolower(s2->label()));
        dpd_buf4_sort(&K, PSIF_HALFT1, rspq, ketCore, braCore, label);
        dpd_buf4_close(&K);

        psio_close(PSIF_HALFT0, 0);

        if(_print) {
            fprintf(outfile, "\tStarting BB second half-transformation.\n");
            fflush(outfile);
        }
        iwl_buf_init(&MBuff, PSIF_MO_BB_TEI, _tolerance, 0, 0);

        
        braCore = DPD_ID(s1, s2, Beta, true);
        ketCore = 0;
        braDisk = DPD_ID(s1, s2, Beta, true);
        ketDisk = 3;
        sprintf(label, "Half-Transformed Ints (%c%c,nn)", tolower(s1->label()), tolower(s2->label()));
        dpd_buf4_init(&J, PSIF_HALFT1, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(_print > 5)
            fprintf(outfile, "Initializing %s, \nin core:(%d|%d) on disk(%d|%d)",
                                label, braCore, ketCore, braDisk, ketDisk);

        braCore = DPD_ID(s1, s2, Beta, true);
        ketCore = DPD_ID(s3, s4, Beta, false);
        braDisk = DPD_ID(s1, s2, Beta, true);
        ketDisk = DPD_ID(s3, s4, Beta, true);
        sprintf(label, "MO Ints (%c%c,%c%c)", tolower(s1->label()), tolower(s2->label()),
                                              tolower(s3->label()), tolower(s4->label()));
        dpd_buf4_init(&K, CC_MISC, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(_print > 5)
            fprintf(outfile, "Initializing %s, \nin core:(%d|%d) on disk(%d|%d)",
                                label, braCore, ketCore, braDisk, ketDisk);

        for(int h=0; h < _nirreps; h++) {
            if (J.params->coltot[h]) {
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

            if(_print > 1) {
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
                    for(int Gr=0; Gr < _nirreps; Gr++) {
                        // Transform ( s1 s2 | n n ) -> ( s1 s2 | n s4 )
                        int Gs = h^Gr;
                        int nrows = _sopi[Gr];
                        int ncols = bOrbsPI4[Gs];
                        int nlinks = _sopi[Gs];
                        int rs = J.col_offset[h][Gr];
                        if(nrows && ncols && nlinks)
                            C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, &J.matrix[h][pq][rs],
                                    nlinks, c4b[Gs][0], ncols, 0.0, TMP[0], _nso);

                        // Transform ( s1 s2 | n s4 ) -> ( s1 s2 | s3 s4 )
                        nrows = bOrbsPI3[Gr];
                        ncols = bOrbsPI4[Gs];
                        nlinks = _sopi[Gr];
                        rs = K.col_offset[h][Gr];
                        if(nrows && ncols && nlinks)
                            C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, c3b[Gr][0], nrows,
                                    TMP[0], _nso, 0.0, &K.matrix[h][pq][rs], ncols);
                    } /* Gr */
                    int p = bIndex1[K.params->roworb[h][pq+n*rowsPerBucket][0]];
                    int q = bIndex2[K.params->roworb[h][pq+n*rowsPerBucket][1]];
                    size_t PQ = INDEX(p,q);
                    for(int rs=0; rs < K.params->coltot[h]; rs++) {
                        int r = bIndex3[K.params->colorb[h][rs][0]];
                        int s = bIndex4[K.params->colorb[h][rs][1]];
                        size_t RS = INDEX(r,s);
                        if(r >= s && RS <= PQ)
                            iwl_buf_wrt_val(&MBuff, p, q, r, s, K.matrix[h][pq][rs], _printTei, outfile, 0);
                    } /* rs */
                } /* pq */
            }
            dpd_buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
            dpd_buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
        }
        dpd_buf4_close(&K);
        dpd_buf4_close(&J);

        iwl_buf_flush(&MBuff, 1);
        iwl_buf_close(&MBuff, 1);
        /*** BB two-electron integral transformation complete ***/
    } // End "if not restricted transformation"


    free_block(TMP);
    delete [] label;

    if(_print){
        fprintf(outfile, "\tTwo-electron integral transformation complete.\n");
        fflush(outfile);
    }
}

}} // End namespaces