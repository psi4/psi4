
#include "dcft.h"
#include <libdpd/dpd.h>
#include <libmints/matrix.h>
#include <libqt/qt.h>
#include "defines.h"

namespace psi { namespace dcft {

/* half_transform(): Routine to transform the last two indices of a dpdbuf4
 * between the MO and SO bases.
 * Based on code originally written by Daniel Crawford in ccenergy/halftrans.cc
 *
 *
 * dpdbuf4 *SO:        Pointer to the (possibly half)-SO basis dpdbuf4 which has already been initialized
 * dpdbuf4 *MO:        Pointer to the (possibly half)-MO basis dpdbuf4 which has already been initialized
 * int *mospi_left:    The number of MO's per irrep for the left upper index
 * int *mospi_right:   The number of MO's per irrep for the right upper index.
 * int **mo_row:       A lookup array.  For a dpdbuf4 with MO indices (ij,ab),
 *                     given the irrep h of ij (= ab) and the irrep of orbital a, the
 *                     array returns the offset of the start of the set of b molecular
 *                     orbitals.
 * int **so_row:       Like mo_row, but for a dpdbuf4 with the last two
 *                     indices in the SO basis.
 * bool backwards:     MO --> SO if true, SO --> MO if false
 * double alpha:       Multiplicative factor for the transformation
 * double beta:        Multiplicative factor for the target
 */

void
DCFTSolver::half_transform(dpdbuf4 *SO, dpdbuf4 *MO, SharedMatrix& C1, SharedMatrix& C2,
        int *mospi_left, int *mospi_right, int **so_row, int **mo_row,
        bool backwards, double alpha, double beta)
{
    dcft_timer_on("DCFTSolver::half_transform");

    int Gc, Gd, cd, pq, ij;

    //Matrix SO_mat(_nIrreps, _soPI, _soPI);
    //Matrix MO_mat(_nIrreps, mospi_right, mospi_left);

    double **X;

        for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(SO, h);
        dpd_buf4_mat_irrep_init(MO, h);

        if(backwards){
            if(alpha != 0.0) dpd_buf4_mat_irrep_rd(MO, h);
            if(beta != 0.0) dpd_buf4_mat_irrep_rd(SO, h);
        }
        else{
            if(alpha != 0.0) dpd_buf4_mat_irrep_rd(SO, h);
            if(beta != 0.0) dpd_buf4_mat_irrep_rd(MO, h);
        }

        for(Gc=0; Gc < nirrep_; Gc++) {
            Gd = h^Gc;
            double **pC1 = C1->pointer(Gc);
            double **pC2 = C2->pointer(Gd);

            cd = mo_row[h][Gc];
            pq = so_row[h][Gc];

            if(mospi_left[Gc] && mospi_right[Gd] && nsopi_[Gc] && nsopi_[Gd]) {

                if(backwards) {
                    X = block_matrix(mospi_left[Gc], nsopi_[Gd]);

                    for(ij = 0; ij < MO->params->rowtot[h]; ij++) {

                        C_DGEMM('n','t', mospi_left[Gc], nsopi_[Gd], mospi_right[Gd], 1.0,
                                &(MO->matrix[h][ij][cd]), mospi_right[Gd], &(pC2[0][0]),
                                mospi_right[Gd], 0.0, &(X[0][0]), nsopi_[Gd]);

                        C_DGEMM('n','n', nsopi_[Gc], nsopi_[Gd], mospi_left[Gc], alpha,
                                &(pC1[0][0]), mospi_left[Gc], &(X[0][0]), nsopi_[Gd],
                                beta, &(SO->matrix[h][ij][pq]), nsopi_[Gd]);
                    }
                }
                else {
                    X = block_matrix(nsopi_[Gc],mospi_right[Gd]);

                    for(ij=0; ij < MO->params->rowtot[h]; ij++) {

                        C_DGEMM('n','n', nsopi_[Gc], mospi_right[Gd], nsopi_[Gd], 1.0,
                                &(SO->matrix[h][ij][pq]), nsopi_[Gd], &(pC2[0][0]), mospi_right[Gd],
                                0.0, &(X[0][0]), mospi_right[Gd]);

                        C_DGEMM('t','n', mospi_left[Gc], mospi_right[Gd], nsopi_[Gc], alpha,
                                &(pC1[0][0]), mospi_left[Gc], &(X[0][0]), mospi_right[Gd],
                                beta, &(MO->matrix[h][ij][cd]), mospi_right[Gd]);

                    }
                }

                free_block(X);
            }
        }

        if(!backwards) dpd_buf4_mat_irrep_wrt(MO, h);
        dpd_buf4_mat_irrep_close(MO, h);

        if(backwards) dpd_buf4_mat_irrep_wrt(SO, h);
        dpd_buf4_mat_irrep_close(SO, h);

    }
    dcft_timer_off("DCFTSolver::half_transform");
}


/* file2_transform(): Convenience routine to transform the last two indices of a dpdbuf4
 * between the MO and SO bases.
 * Based on code originally written by Daniel Crawford in ccenergy/halftrans.cc
 *
 *
 * dpdfile2 *SO:       Pointer to the SO basis dpdfile2 which has already been initialized
 * dpdfile2 *MO:       Pointer to the MO basis dpdfile2 which has already been initialized
 * Matrix *C:          Pointer to the transformation Matrix object.
 * bool backwards:     MO --> SO if true, SO --> MO if false
 */

void
DCFTSolver::file2_transform(dpdfile2 *SO, dpdfile2 *MO, SharedMatrix C, bool backwards)
{
    dcft_timer_on("DCFTSolver::file2_transform");

    if(backwards) {
        Matrix MO_mat(MO);
        MO_mat.back_transform(C);
        MO_mat.write_to_dpdfile2(SO);
    }
    else {
        Matrix SO_mat(SO);
        SO_mat.transform(C);
        SO_mat.write_to_dpdfile2(MO);
    }
    dcft_timer_off("DCFTSolver::file2_transform");
}


}} // End namespaces
