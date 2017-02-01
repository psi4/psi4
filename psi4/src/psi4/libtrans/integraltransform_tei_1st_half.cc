/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include "psi4/psifiles.h"
#include "mospace.h"
#define EXTERN
#include "psi4/libdpd/dpd.gbl"

using namespace psi;
;

void
IntegralTransform::transform_tei_first_half(const std::shared_ptr<MOSpace> s1, const std::shared_ptr<MOSpace> s2)
{
    check_initialized();

    // This can be safely called - it returns immediately if the SO ints are already sorted
    presort_so_tei();

    char *label = new char[100];

    // Grab the transformation coefficients
    SharedMatrix c1a = aMOCoefficients_[s1->label()];
    SharedMatrix c1b = bMOCoefficients_[s1->label()];
    SharedMatrix c2a = aMOCoefficients_[s2->label()];
    SharedMatrix c2b = bMOCoefficients_[s2->label()];
    // And the number of orbitals per irrep
    int *aOrbsPI1 = aOrbsPI_[s1->label()];
    int *bOrbsPI1 = bOrbsPI_[s1->label()];
    int *aOrbsPI2 = aOrbsPI_[s2->label()];
    int *bOrbsPI2 = bOrbsPI_[s2->label()];

    // Grab control of DPD for now, but store the active number to restore it later
    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    int nBuckets;
    int thisBucketRows;
    size_t rowsPerBucket;
    size_t rowsLeft;
    size_t memFree;

    double **TMP = block_matrix(nso_, nso_);

    /*** AA/AB two-electron integral transformation ***/

    if(print_) {
        if(transformationType_ == Restricted){
            outfile->Printf( "\tStarting first half-transformation.\n");
        }else{
            outfile->Printf( "\tStarting AA/AB first half-transformation.\n");
        }

    }

    psio_->open(PSIF_SO_PRESORT, PSIO_OPEN_OLD);
    psio_->open(PSIF_HALFT0, PSIO_OPEN_NEW);

    dpdbuf4 J, K;
    global_dpd_->buf4_init(&J, PSIF_SO_PRESORT, 0, DPD_ID("[n>=n]+"), DPD_ID("[n,n]"),
                  DPD_ID("[n>=n]+"), DPD_ID("[n>=n]+"), 0, "SO Ints (nn|nn)");

    int braCore = DPD_ID("[n>=n]+");
    int ketCore = DPD_ID(s1, s2, Alpha, false);
    int braDisk = DPD_ID("[n>=n]+");
    int ketDisk = DPD_ID(s1, s2, Alpha, true);
    sprintf(label, "Half-Transformed Ints (nn|%c%c)", toupper(s1->label()), toupper(s2->label()));
    global_dpd_->buf4_init(&K, PSIF_HALFT0, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
    if(print_ > 5)
        outfile->Printf( "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
                            label, braCore, ketCore, braDisk, ketDisk);

    for(int h=0; h < nirreps_; h++) {
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

        if(print_ > 1) {
            outfile->Printf( "\th = %d; memfree         = %lu\n", h, memFree);
            outfile->Printf( "\th = %d; rows_per_bucket = %lu\n", h, rowsPerBucket);
            outfile->Printf( "\th = %d; rows_left       = %lu\n", h, rowsLeft);
            outfile->Printf( "\th = %d; nbuckets        = %d\n", h, nBuckets);

        }

        global_dpd_->buf4_mat_irrep_init_block(&J, h, rowsPerBucket);
        global_dpd_->buf4_mat_irrep_init_block(&K, h, rowsPerBucket);

        for(int n=0; n < nBuckets; n++){
            if(nBuckets == 1)
                thisBucketRows = rowsPerBucket;
            else
                thisBucketRows = (n < nBuckets-1) ? rowsPerBucket : rowsLeft;
            global_dpd_->buf4_mat_irrep_rd_block(&J, h, n*rowsPerBucket, thisBucketRows);
            for(int pq=0; pq < thisBucketRows; pq++) {
                for(int Gr=0; Gr < nirreps_; Gr++) {
                    // Transform ( n n | n n ) -> ( n n | n S2 )
                    int Gs = h^Gr;
                    int nrows = sopi_[Gr];
                    int ncols = aOrbsPI2[Gs];
                    int nlinks = sopi_[Gs];
                    int rs = J.col_offset[h][Gr];
                    double **pc2a = c2a->pointer(Gs);
                    if(nrows && ncols && nlinks)
                        C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, &J.matrix[h][pq][rs],
                                nlinks, pc2a[0], ncols, 0.0, TMP[0], nso_);
                    //TODO else if s1->label() == MOSPACE_NIL, copy buffer...

                    // Transform ( n n | n S2 ) -> ( n n | S1 S2 )
                    nrows = aOrbsPI1[Gr];
                    ncols = aOrbsPI2[Gs];
                    nlinks = sopi_[Gr];
                    rs = K.col_offset[h][Gr];
                    double **pc1a = c1a->pointer(Gr);
                    if(nrows && ncols && nlinks)
                        C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, pc1a[0], nrows,
                                TMP[0], nso_, 0.0, &K.matrix[h][pq][rs], ncols);
                    //TODO else if s2->label() == MOSPACE_NIL, copy buffer...
                } /* Gr */
            } /* pq */
            global_dpd_->buf4_mat_irrep_wrt_block(&K, h, n*rowsPerBucket, thisBucketRows);
        }
        global_dpd_->buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
        global_dpd_->buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&J);

    if(print_) {
        if(transformationType_ == Restricted){
            outfile->Printf( "\tSorting half-transformed integrals.\n");
        }else{
            outfile->Printf( "\tSorting AA/AB half-transformed integrals.\n");
        }

    }

    psio_->open(aHtIntFile_, PSIO_OPEN_NEW);

    braCore = braDisk = DPD_ID("[n>=n]+");
    ketCore = ketDisk = DPD_ID(s1, s2, Alpha, true);
    sprintf(label, "Half-Transformed Ints (nn|%c%c)", toupper(s1->label()), toupper(s2->label()));
    global_dpd_->buf4_init(&K, PSIF_HALFT0, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
    if(print_ > 5)
        outfile->Printf( "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
                            label, braCore, ketCore, braDisk, ketDisk);
    sprintf(label, "Half-Transformed Ints (%c%c|nn)", toupper(s1->label()), toupper(s2->label()));
    global_dpd_->buf4_sort(&K, aHtIntFile_, rspq, ketCore, braCore, label);
    global_dpd_->buf4_close(&K);

    psio_->close(aHtIntFile_, 1);
    psio_->close(PSIF_HALFT0, 0);

    if(transformationType_ != Restricted){
        /*** BB two-electron integral transformation ***/
        if(print_) {
            outfile->Printf( "\tStarting BB first half-transformation.\n");

        }

        psio_->open(PSIF_HALFT0, PSIO_OPEN_NEW);

        global_dpd_->buf4_init(&J, PSIF_SO_PRESORT, 0, DPD_ID("[n>=n]+"), DPD_ID("[n,n]"),
                      DPD_ID("[n>=n]+"), DPD_ID("[n>=n]+"), 0, "SO Ints (nn|nn)");

        braCore = DPD_ID("[n>=n]+");
        ketCore = DPD_ID(s1, s2, Beta, false);
        braDisk = DPD_ID("[n>=n]+");
        ketDisk = DPD_ID(s1, s2, Beta, true);
        sprintf(label, "Half-Transformed Ints (nn|%c%c)", tolower(s1->label()), tolower(s2->label()));
        global_dpd_->buf4_init(&K, PSIF_HALFT0, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(print_ > 5)
            outfile->Printf( "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
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

            if(print_ > 1){
                outfile->Printf( "\th = %d; memfree         = %lu\n", h, memFree);
                outfile->Printf( "\th = %d; rows_per_bucket = %lu\n", h, rowsPerBucket);
                outfile->Printf( "\th = %d; rows_left       = %lu\n", h, rowsLeft);
                outfile->Printf( "\th = %d; nbuckets        = %d\n", h, nBuckets);

            }

            global_dpd_->buf4_mat_irrep_init_block(&J, h, rowsPerBucket);
            global_dpd_->buf4_mat_irrep_init_block(&K, h, rowsPerBucket);

            for(int n=0; n < nBuckets; n++) {
                if(nBuckets == 1)
                    thisBucketRows = rowsPerBucket;
                else
                    thisBucketRows = (n < nBuckets-1) ? rowsPerBucket : rowsLeft;
                global_dpd_->buf4_mat_irrep_rd_block(&J, h, n*rowsPerBucket, thisBucketRows);
                for(int pq=0; pq < thisBucketRows; pq++) {
                    for(int Gr=0; Gr < nirreps_; Gr++) {
                        // Transform ( n n | n n ) -> ( n n | n s2 )
                        int Gs = h^Gr;
                        int nrows = sopi_[Gr];
                        int ncols = bOrbsPI2[Gs];
                        int nlinks = sopi_[Gs];
                        int rs = J.col_offset[h][Gr];
                        double **pc2b = c2b->pointer(Gs);
                        if(nrows && ncols && nlinks)
                            C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, &J.matrix[h][pq][rs],
                            nlinks, pc2b[0], ncols, 0.0, TMP[0], nso_);
                        //TODO else if s2->label() == MOSPACE_NIL, copy buffer...

                        // Transform ( n n | n s2 ) -> ( n n | s1 s2 )
                        nrows = bOrbsPI1[Gr];
                        ncols = bOrbsPI2[Gs];
                        nlinks = sopi_[Gr];
                        rs = K.col_offset[h][Gr];
                        double **pc1b = c1b->pointer(Gr);
                        if(nrows && ncols && nlinks)
                            C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, pc1b[0], nrows,
                                    TMP[0], nso_, 0.0, &K.matrix[h][pq][rs], ncols);
                        //TODO else if s1->label() == MOSPACE_NIL, copy buffer...
                    } /* Gr */
                } /* pq */
                global_dpd_->buf4_mat_irrep_wrt_block(&K, h, n*rowsPerBucket, thisBucketRows);
            }
            global_dpd_->buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
            global_dpd_->buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
        }
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_close(&J);

        if(print_) {
            outfile->Printf( "\tSorting BB half-transformed integrals.\n");

        }

        psio_->open(bHtIntFile_, PSIO_OPEN_NEW);

        braCore = DPD_ID("[n>=n]+");
        ketCore = DPD_ID(s1, s2, Beta, true);
        braDisk = DPD_ID("[n>=n]+");
        ketDisk = DPD_ID(s1, s2, Beta, true);
        sprintf(label, "Half-Transformed Ints (nn|%c%c)", tolower(s1->label()), tolower(s2->label()));
        global_dpd_->buf4_init(&K, PSIF_HALFT0, 0, braCore, ketCore, braDisk, ketDisk, 0, label);
        if(print_ > 5)
            outfile->Printf( "Initializing %s, in core:(%d|%d) on disk(%d|%d)\n",
                                label, braCore, ketCore, braDisk, ketDisk);

        sprintf(label, "Half-Transformed Ints (%c%c|nn)", tolower(s1->label()), tolower(s2->label()));
        global_dpd_->buf4_sort(&K, bHtIntFile_, rspq, ketCore, braCore, label);
        global_dpd_->buf4_close(&K);

        psio_->close(bHtIntFile_, 1);
        psio_->close(PSIF_HALFT0, 0);
    } // End "if not restricted transformation"

    psio_->close(PSIF_SO_PRESORT, keepDpdSoInts_);

    free_block(TMP);
    delete [] label;

    if(print_){
        outfile->Printf( "\tFirst half integral transformation complete.\n");

    }

    // Hand DPD control back to the user
    dpd_set_default(currentActiveDPD);
}
