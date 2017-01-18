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
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/matrix.h"
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include "psi4/psifiles.h"
#include "mospace.h"
#define EXTERN
#include "psi4/libdpd/dpd.gbl"

using namespace psi;

void
IntegralTransform::backtransform_tpdm_restricted()
{
    check_initialized();

    // This can be safely called - it returns immediately if the MO TPDM is already sorted
    presort_mo_tpdm_restricted();

    // Grab the transformation coefficients
    SharedMatrix c = aMOCoefficients_[MOSPACE_ALL];

    // Grab control of DPD for now, but store the active number to restore it later
    int currentActiveDPD = psi::dpd_default;
    dpd_set_default(myDPDNum_);

    int nBuckets;
    int thisBucketRows;
    size_t rowsPerBucket;
    size_t rowsLeft;
    size_t memFree;

    double **TMP = block_matrix(nso_, nso_);

    /*** first half transformation ***/

    if(print_) {
        outfile->Printf( "\tStarting first half-transformation.\n");

    }

    psio_->open(PSIF_TPDM_PRESORT, PSIO_OPEN_OLD);
    psio_->open(PSIF_TPDM_HALFTRANS, PSIO_OPEN_NEW);

    dpdbuf4 J, K;
    global_dpd_->buf4_init(&J, PSIF_TPDM_PRESORT, 0, DPD_ID("[A>=A]+"), DPD_ID("[A,A]"),
                  DPD_ID("[A>=A]+"), DPD_ID("[A>=A]+"), 0, "MO TPDM (AA|AA)");
//    global_dpd_->buf4_print(&J, "outfile", 1);
    global_dpd_->buf4_init(&K, PSIF_TPDM_HALFTRANS, 0, DPD_ID("[A>=A]+"), DPD_ID("[n,n]"),
                  DPD_ID("[A>=A]+"), DPD_ID("[n>=n]+"), 0, "Half-Transformed TPDM (AA|nn)");

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
                    // Transform ( a a | a a ) -> ( a a | a n )
                    int Gs = h^Gr;
                    int nrows = sopi_[Gr];
                    int ncols = mopi_[Gs];
                    int nlinks = mopi_[Gs];
                    int rs = J.col_offset[h][Gr];
                    double **pc = c->pointer(Gs);
                    if(nrows && ncols && nlinks)
                        C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, &J.matrix[h][pq][rs],
                                nlinks, pc[0], ncols, 0.0, TMP[0], nso_);

                    // Transform ( a a | a n ) -> ( a a | n n )
                    nrows = sopi_[Gr];
                    ncols = sopi_[Gs];
                    nlinks = mopi_[Gr];
                    rs = K.col_offset[h][Gr];
                    pc = c->pointer(Gr);
                    if(nrows && ncols && nlinks)
                        C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, pc[0], nrows,
                                TMP[0], nso_, 0.0, &K.matrix[h][pq][rs], ncols);
                } /* Gr */
            } /* pq */
            global_dpd_->buf4_mat_irrep_wrt_block(&K, h, n*rowsPerBucket, thisBucketRows);
        }
        global_dpd_->buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
        global_dpd_->buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&J);

    psio_->close(PSIF_TPDM_PRESORT, keepDpdMoTpdm_);

    if(print_) {
        outfile->Printf( "\tSorting half-transformed TPDM.\n");

    }

    global_dpd_->buf4_init(&K, PSIF_TPDM_HALFTRANS, 0, DPD_ID("[A>=A]+"), DPD_ID("[n>=n]+"),
                  DPD_ID("[A>=A]+"), DPD_ID("[n>=n]+"), 0, "Half-Transformed TPDM (AA|nn)");
//    global_dpd_->buf4_print(&K, "outfile", 1);
    global_dpd_->buf4_sort(&K, PSIF_TPDM_HALFTRANS, rspq, DPD_ID("[n>=n]+"), DPD_ID("[A>=A]+"), "Half-Transformed TPDM (nn|AA)");
    global_dpd_->buf4_close(&K);

    if(print_){
        outfile->Printf( "\tFirst half integral transformation complete.\n");

    }

    psio_->open(PSIF_AO_TPDM, PSIO_OPEN_NEW);

    global_dpd_->buf4_init(&J, PSIF_TPDM_HALFTRANS, 0, DPD_ID("[n>=n]+"), DPD_ID("[A,A]"),
                  DPD_ID("[n>=n]+"), DPD_ID("[A>=A]+"), 0, "Half-Transformed TPDM (nn|AA)");

    global_dpd_->buf4_init(&K, PSIF_AO_TPDM, 0, DPD_ID("[n>=n]+"), DPD_ID("[n,n]"),
                  DPD_ID("[n>=n]+"), DPD_ID("[n>=n]+"), 0, "SO Basis TPDM (nn|nn)");

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
                    // Transform ( n n | a a ) -> ( n n | a n )
                    int Gs = h^Gr;
                    int nrows = sopi_[Gr];
                    int ncols = mopi_[Gs];
                    int nlinks = mopi_[Gs];
                    int rs = J.col_offset[h][Gr];
                    double **pc = c->pointer(Gs);
                    if(nrows && ncols && nlinks)
                        C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, &J.matrix[h][pq][rs],
                                nlinks, pc[0], ncols, 0.0, TMP[0], nso_);

                    // Transform ( n n | n a ) -> ( n n | n n )
                    nrows = sopi_[Gr];
                    ncols = sopi_[Gs];
                    nlinks = mopi_[Gr];
                    rs = K.col_offset[h][Gr];
                    pc = c->pointer(Gr);
                    if(nrows && ncols && nlinks)
                        C_DGEMM('n', 'n', nrows, ncols, nlinks, 1.0, pc[0], nrows,
                                TMP[0], nso_, 0.0, &K.matrix[h][pq][rs], ncols);
                } /* Gr */
            } /* pq */
            sort_so_tpdm(&K, h, n*rowsPerBucket, thisBucketRows, (h==0 && n==0));
            if(write_dpd_so_tpdm_)
                global_dpd_->buf4_mat_irrep_wrt_block(&K, h, n*rowsPerBucket, thisBucketRows);
        }
        global_dpd_->buf4_mat_irrep_close_block(&J, h, rowsPerBucket);
        global_dpd_->buf4_mat_irrep_close_block(&K, h, rowsPerBucket);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&J);

    free_block(TMP);

    psio_->close(PSIF_TPDM_HALFTRANS, keepHtTpdm_);
    psio_->close(PSIF_AO_TPDM, 1);

    // Hand DPD control back to the user
    dpd_set_default(currentActiveDPD);
}
