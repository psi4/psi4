#include "integraltransform.h"
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
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
IntegralTransform::transform_tei_first_half(const shared_ptr<MOSpace> s1, const shared_ptr<MOSpace> s2)
{
    check_initialized();

    // This can be safely called - it returns immediately if the SO ints are already sorted
    presort_so_tei();

    char *label = new char[100];

    // Grab the transformation coefficients
    double ***c1a = _aMOCoefficients[s1->label()];
    double ***c1b = _bMOCoefficients[s1->label()];
    double ***c2a = _aMOCoefficients[s2->label()];
    double ***c2b = _bMOCoefficients[s2->label()];
    // And the number of orbitals per irrep
    int *aOrbsPI1 = _aOrbsPI[s1->label()];
    int *bOrbsPI1 = _bOrbsPI[s1->label()];
    int *aOrbsPI2 = _aOrbsPI[s2->label()];
    int *bOrbsPI2 = _bOrbsPI[s2->label()];
    // The reindexing arrays
    int *aIndex1 = _aIndices[s1->label()];
    int *bIndex1 = _bIndices[s1->label()];
    int *aIndex2 = _aIndices[s2->label()];
    int *bIndex2 = _bIndices[s2->label()];

    // Grab control of DPD for now, but store the active number to restore it later
    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(_myDPDNum);

    int nBuckets;
    int thisBucketRows;
    size_t rowsPerBucket;
    size_t rowsLeft;
    size_t memFree;

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

    _psio->open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);
    _psio->open(PSIF_HALFT0, PSIO_OPEN_NEW);

    dpdbuf4 J, K;
    dpd_buf4_init(&J, PSIF_SO_PRESORT, 0, 3, 0, 3, 3, 0, "SO Ints (nn|nn)");

    int braCore = 3;
    int ketCore = DPD_ID(s1, s2, Alpha, false);
    int braDisk = 3;
    int ketDisk = DPD_ID(s1, s2, Alpha, true);
    sprintf(label, "Half-Transformed Ints (nn|%c%c)", toupper(s1->label()), toupper(s2->label()));
    dpd_buf4_init(&K, PSIF_HALFT0, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
    if(_print > 5)
        fprintf(outfile, "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
                            label, braCore, ketCore, braDisk, ketDisk);

    for(int h=0; h < _nirreps; h++) {
        if(J.params->coltot[h] && J.params->rowtot[h]) {
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
                    //TODO else if s1->label() == MOSPACE_NIL, copy buffer...

                    // Transform ( n n | n S2 ) -> ( n n | S1 S2 )
                    nrows = aOrbsPI1[Gr];
                    ncols = aOrbsPI2[Gs];
                    nlinks = _sopi[Gr];
                    rs = K.col_offset[h][Gr];
                    if(nrows && ncols && nlinks)
                        C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, c1a[Gr][0], nrows,
                                TMP[0], _nso, 0.0, &K.matrix[h][pq][rs], ncols);
                    //TODO else if s2->label() == MOSPACE_NIL, copy buffer...
                } /* Gr */
            } /* pq */
            dpd_buf4_mat_irrep_wrt_block(&K, h, n*rowsPerBucket, thisBucketRows);
        }
        dpd_buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
        dpd_buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
    }
    dpd_buf4_close(&K);
    dpd_buf4_close(&J);

    if(_print) {
        if(_transformationType == Restricted){
            fprintf(outfile, "\tSorting half-transformed integrals.\n");
        }else{
            fprintf(outfile, "\tSorting AA/AB half-transformed integrals.\n");
        }
        fflush(outfile);
    }

    _psio->open(_aHtIntFile, PSIO_OPEN_NEW);

    braCore = braDisk = 3;
    ketCore = ketDisk = DPD_ID(s1, s2, Alpha, true);
    sprintf(label, "Half-Transformed Ints (nn|%c%c)", toupper(s1->label()), toupper(s2->label()));
    dpd_buf4_init(&K, PSIF_HALFT0, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
    if(_print > 5)
        fprintf(outfile, "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
                            label, braCore, ketCore, braDisk, ketDisk);
    sprintf(label, "Half-Transformed Ints (%c%c|nn)", toupper(s1->label()), toupper(s2->label()));
    dpd_buf4_sort(&K, _aHtIntFile, rspq, ketCore, braCore, label);
    dpd_buf4_close(&K);

    _psio->close(_aHtIntFile, 1);
    _psio->close(PSIF_HALFT0, 0);

    if(_transformationType != Restricted){
        /*** BB two-electron integral transformation ***/
        if(_print) {
            fprintf(outfile, "\tStarting BB first half-transformation.\n");
            fflush(outfile);
        }

        _psio->open(PSIF_HALFT0, PSIO_OPEN_NEW);

        dpd_buf4_init(&J, PSIF_SO_PRESORT, 0, 3, 0, 3, 3, 0, "SO Ints (nn|nn)");

        braCore = 3;
        ketCore = DPD_ID(s1, s2, Beta, false);
        braDisk = 3;
        ketDisk = DPD_ID(s1, s2, Beta, true);
        sprintf(label, "Half-Transformed Ints (nn|%c%c)", tolower(s1->label()), tolower(s2->label()));
        dpd_buf4_init(&K, PSIF_HALFT0, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(_print > 5)
            fprintf(outfile, "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
                                label, braCore, ketCore, braDisk, ketDisk);

        for(int h=0; h < _nirreps; h++) {
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
                        //TODO else if s2->label() == MOSPACE_NIL, copy buffer...

                        // Transform ( n n | n s2 ) -> ( n n | s1 s2 )
                        nrows = bOrbsPI1[Gr];
                        ncols = bOrbsPI2[Gs];
                        nlinks = _sopi[Gr];
                        rs = K.col_offset[h][Gr];
                        if(nrows && ncols && nlinks)
                            C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, c1b[Gr][0], nrows,
                                    TMP[0], _nso, 0.0, &K.matrix[h][pq][rs], ncols);
                        //TODO else if s1->label() == MOSPACE_NIL, copy buffer...
                    } /* Gr */
                } /* pq */
                dpd_buf4_mat_irrep_wrt_block(&K, h, n*rowsPerBucket, thisBucketRows);
            }
            dpd_buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
            dpd_buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
        }
        dpd_buf4_close(&K);
        dpd_buf4_close(&J);

        if(_print) {
            fprintf(outfile, "\tSorting BB half-transformed integrals.\n");
            fflush(outfile);
        }

        _psio->open(_bHtIntFile, PSIO_OPEN_NEW);

        braCore = 3;
        ketCore = DPD_ID(s1, s2, Beta, true);
        braDisk = 3;
        ketDisk = DPD_ID(s1, s2, Beta, true);
        sprintf(label, "Half-Transformed Ints (nn|%c%c)", tolower(s1->label()), tolower(s2->label()));
        dpd_buf4_init(&K, PSIF_HALFT0, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(_print > 5)
            fprintf(outfile, "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
                                label, braCore, ketCore, braDisk, ketDisk);

        sprintf(label, "Half-Transformed Ints (%c%c|nn)", tolower(s1->label()), tolower(s2->label()));
        dpd_buf4_sort(&K, _bHtIntFile, rspq, ketCore, braCore, label);
        dpd_buf4_close(&K);

        _psio->close(_bHtIntFile, 1);
        _psio->close(PSIF_HALFT0, 0);
    } // End "if not restricted transformation"

    _psio->close(PSIF_SO_PRESORT, _keepDpdSoInts);

    free_block(TMP);
    delete [] label;

    if(_print){
        fprintf(outfile, "\tFirst half integral transformation complete.\n");
        fflush(outfile);
    }

    // Hand DPD control back to the user
    dpd_set_default(currentActiveDPD);
}
