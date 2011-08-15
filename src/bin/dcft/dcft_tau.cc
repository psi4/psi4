#include "dcft.h"
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libiwl/iwl.hpp>
#include <psifiles.h>
#include <libtrans/integraltransform.h>
#include "defines.h"

using namespace std;

namespace psi{ namespace dcft{

/**
 * Forms Tau in the MO basis from the Lambda tensors, before transforming back
 * to the AO basis.
 */
void
DCFTSolver::build_tau()
{
    dpdbuf4 L1, L2;
    dpdfile2 T_OO, T_oo, T_VV, T_vv;

    dpd_file2_init(&T_OO, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('O'), _ints->DPD_ID('O'), "Tau <O|O>");
    dpd_file2_init(&T_oo, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('o'), _ints->DPD_ID('o'), "Tau <o|o>");
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('V'), _ints->DPD_ID('V'), "Tau <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('v'), _ints->DPD_ID('v'), "Tau <v|v>");

    dpd_buf4_init(&L1, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID("[O,O]"), _ints->DPD_ID("[V,V]"),
                  _ints->DPD_ID("[O,O]"), _ints->DPD_ID("[V,V]"),
                  0, "Lambda <OO|VV>");
    dpd_buf4_init(&L2, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID("[O,O]"), _ints->DPD_ID("[V,V]"),
                  _ints->DPD_ID("[O,O]"), _ints->DPD_ID("[V,V]"),
                  0, "Lambda <OO|VV>");
    /*
     * Tau_IJ = -1/2 Lambda_IKAB Lambda_JKAB
     */
    dpd_contract442(&L1, &L2, &T_OO, 0, 0, -0.5, 0.0);
    /*
     * Tau_AB = +1/2 Lambda_IJAC Lambda_IJBC
     */
    dpd_contract442(&L1, &L2, &T_VV, 2, 2, 0.5, 0.0);
    dpd_buf4_close(&L1);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L1, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID("[o,o]"), _ints->DPD_ID("[v,v]"),
                  _ints->DPD_ID("[o,o]"), _ints->DPD_ID("[v,v]"),
                  0, "Lambda <oo|vv>");
    dpd_buf4_init(&L2, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID("[o,o]"), _ints->DPD_ID("[v,v]"),
                  _ints->DPD_ID("[o,o]"), _ints->DPD_ID("[v,v]"),
                  0, "Lambda <oo|vv>");
    /*
     * Tau_ij = -1/2 Lambda_ikab Lambda_jkab
     */
    dpd_contract442(&L1, &L2, &T_oo, 0, 0, -0.5, 0.0);
    /*
     * Tau_ab = +1/2 Lambda_ijac Lambda_ijbc
     */
    dpd_contract442(&L1, &L2, &T_vv, 2, 2, 0.5, 0.0);
    dpd_buf4_close(&L1);
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L1, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID("[O,o]"), _ints->DPD_ID("[V,v]"),
                  _ints->DPD_ID("[O,o]"), _ints->DPD_ID("[V,v]"),
                  0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&L2, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID("[O,o]"), _ints->DPD_ID("[V,v]"),
                  _ints->DPD_ID("[O,o]"), _ints->DPD_ID("[V,v]"),
                  0, "Lambda <Oo|Vv>");
    /*
     * Tau_IJ -= 1/2 Lambda_IkAb Lambda_JkAb - 1/2 Lambda_IkaB Lambda_JkaB
     */
    dpd_contract442(&L1, &L2, &T_OO, 0, 0, -1.0, 1.0);
    /*
     * Tau_ij -= 1/2 Lambda_KiAb Lambda_KjAb - 1/2 Lambda_KiaB Lambda_KjaB
     */
    dpd_contract442(&L1, &L2, &T_oo, 1, 1, -1.0, 1.0);
    /*
     * Tau_AB += 1/2 Lambda_IjAc Lambda_IjBc + 1/2 Lambda_iJAc Lambda_iJBc
     */
    dpd_contract442(&L1, &L2, &T_VV, 2, 2, 1.0, 1.0);
    /*
     * Tau_ab += 1/2 Lambda_IjCa Lambda_IjCb + 1/2 Lambda_iJCa Lambda_iJCb
     */
    dpd_contract442(&L1, &L2, &T_vv, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L1);
    dpd_buf4_close(&L2);

    dpd_file2_mat_init(&T_OO);
    dpd_file2_mat_init(&T_oo);
    dpd_file2_mat_init(&T_VV);
    dpd_file2_mat_init(&T_vv);
    dpd_file2_mat_rd(&T_OO);
    dpd_file2_mat_rd(&T_oo);
    dpd_file2_mat_rd(&T_VV);
    dpd_file2_mat_rd(&T_vv);

//  Zero AO tau arrays before copmputing it in the MO basis
    a_tau_->zero();
    b_tau_->zero();

    for(int h = 0; h < nirrep_; ++h){
        if(nsopi_[h] == 0) continue;

        double **temp = block_matrix(nsopi_[h], nsopi_[h]);
        /*
         * Backtransform the Tau matrices to the AO basis
         * tauTrans = C moTau Ct
         */
        double **paOccC = aocc_c_->pointer(h);
        double **pbOccC = bocc_c_->pointer(h);
        double **paVirC = avir_c_->pointer(h);
        double **pbVirC = bvir_c_->pointer(h);
        double **pa_tau_ = a_tau_->pointer(h);
        double **pb_tau_ = b_tau_->pointer(h);

        // Alpha occupied
        if(naoccpi_[h] && nsopi_[h]){
            C_DGEMM('n', 'n', nsopi_[h], naoccpi_[h], naoccpi_[h], 1.0, paOccC[0], naoccpi_[h],
                    T_OO.matrix[h][0], naoccpi_[h], 0.0, temp[0], nsopi_[h]);
            C_DGEMM('n', 't', nsopi_[h], nsopi_[h], naoccpi_[h], 1.0, temp[0], nsopi_[h],
                    paOccC[0], naoccpi_[h], 1.0, pa_tau_[0], nsopi_[h]);
        }
        // Beta occupied
        if(nboccpi_[h] && nsopi_[h]){
            C_DGEMM('n', 'n', nsopi_[h], nboccpi_[h], nboccpi_[h], 1.0, pbOccC[0], nboccpi_[h],
                    T_oo.matrix[h][0], nboccpi_[h], 0.0, temp[0], nsopi_[h]);
            C_DGEMM('n', 't', nsopi_[h], nsopi_[h], nboccpi_[h], 1.0, temp[0], nsopi_[h],
                    pbOccC[0], nboccpi_[h], 1.0, pb_tau_[0], nsopi_[h]);
        }
        // Alpha virtual
        if(navirpi_[h] && nsopi_[h]){
            C_DGEMM('n', 'n', nsopi_[h], navirpi_[h], navirpi_[h], 1.0, paVirC[0], navirpi_[h],
                    T_VV.matrix[h][0], navirpi_[h], 0.0, temp[0], nsopi_[h]);
            C_DGEMM('n', 't', nsopi_[h], nsopi_[h], navirpi_[h], 1.0, temp[0], nsopi_[h],
                    paVirC[0], navirpi_[h], 1.0, pa_tau_[0], nsopi_[h]);
        }
        // Beta virtual
        if(nbvirpi_[h] && nsopi_[h]){
            C_DGEMM('n', 'n', nsopi_[h], nbvirpi_[h], nbvirpi_[h], 1.0, pbVirC[0], nbvirpi_[h],
                    T_vv.matrix[h][0], nbvirpi_[h], 0.0, temp[0], nsopi_[h]);
            C_DGEMM('n', 't', nsopi_[h], nsopi_[h], nbvirpi_[h], 1.0, temp[0], nsopi_[h],
                    pbVirC[0], nbvirpi_[h], 1.0, pb_tau_[0], nsopi_[h]);
        }

        free_block(temp);
    }

    dpd_file2_close(&T_OO);
    dpd_file2_close(&T_oo);
    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);

}

void
DCFTSolver::compute_tau_squared()
{

    dpdfile2 T_OO, T_oo, T_VV, T_vv;
    dpdfile2 TT_OO, TT_oo, TT_VV, TT_vv;

    dpd_file2_init(&T_OO, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('O'), _ints->DPD_ID('O'), "Tau <O|O>");
    dpd_file2_init(&T_oo, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('o'), _ints->DPD_ID('o'), "Tau <o|o>");
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('V'), _ints->DPD_ID('V'), "Tau <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('v'), _ints->DPD_ID('v'), "Tau <v|v>");

    dpd_file2_init(&TT_OO, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('O'), _ints->DPD_ID('O'), "Tau^2 <O|O>");
    dpd_file2_init(&TT_oo, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('o'), _ints->DPD_ID('o'), "Tau^2 <o|o>");
    dpd_file2_init(&TT_VV, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('V'), _ints->DPD_ID('V'), "Tau^2 <V|V>");
    dpd_file2_init(&TT_vv, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('v'), _ints->DPD_ID('v'), "Tau^2 <v|v>");

    dpd_file2_mat_init(&T_OO);
    dpd_file2_mat_init(&T_oo);
    dpd_file2_mat_init(&T_VV);
    dpd_file2_mat_init(&T_vv);
    dpd_file2_mat_init(&TT_OO);
    dpd_file2_mat_init(&TT_oo);
    dpd_file2_mat_init(&TT_VV);
    dpd_file2_mat_init(&TT_vv);
    dpd_file2_mat_rd(&T_OO);
    dpd_file2_mat_rd(&T_oo);
    dpd_file2_mat_rd(&T_VV);
    dpd_file2_mat_rd(&T_vv);

    // Compute the Tau^2 correction to Tau in the MO basis

    for(int h = 0; h < nirrep_; ++h){
        if(nsopi_[h] == 0) continue;

        // Alpha occupied
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int j = 0 ; j < naoccpi_[h]; ++j){
                TT_OO.matrix[h][i][j] = (-1.0) * T_OO.matrix[h][i][j] * T_OO.matrix[h][i][j];
            }
        }
        // Beta occupied
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int j = 0 ; j < nboccpi_[h]; ++j){
                TT_oo.matrix[h][i][j] = (-1.0) * T_oo.matrix[h][i][j] * T_oo.matrix[h][i][j];
            }
        }
        // Alpha virtual
        for(int i = 0 ; i < navirpi_[h]; ++i){
            for(int j = 0 ; j < navirpi_[h]; ++j){
                TT_VV.matrix[h][i][j] = T_VV.matrix[h][i][j] * T_VV.matrix[h][i][j];
            }
        }
        // Beta virtual
        for(int i = 0 ; i < nbvirpi_[h]; ++i){
            for(int j = 0 ; j < nbvirpi_[h]; ++j){
                TT_vv.matrix[h][i][j] = T_vv.matrix[h][i][j] * T_vv.matrix[h][i][j];
            }
        }
    }

    a_tautau_->zero();
    b_tautau_->zero();

    for(int h = 0; h < nirrep_; ++h){
        if(nsopi_[h] == 0) continue;

        double **temp = block_matrix(nsopi_[h], nsopi_[h]);
        double **paOccC = aocc_c_->pointer(h);
        double **pbOccC = bocc_c_->pointer(h);
        double **paVirC = avir_c_->pointer(h);
        double **pbVirC = bvir_c_->pointer(h);
        double **pa_tautau_ = a_tautau_->pointer(h);
        double **pb_tautau_ = b_tautau_->pointer(h);

        /*
         * Backtransform the Tau^2 correction to the AO basis
         * tau^2Trans = C moTau^2 Ct
         */

        // Alpha occupied
        if(naoccpi_[h] && nsopi_[h]){
            C_DGEMM('n', 'n', nsopi_[h], naoccpi_[h], naoccpi_[h], 1.0, paOccC[0], naoccpi_[h],
                    TT_OO.matrix[h][0], naoccpi_[h], 0.0, temp[0], nsopi_[h]);
            C_DGEMM('n', 't', nsopi_[h], nsopi_[h], naoccpi_[h], 1.0, temp[0], nsopi_[h],
                    paOccC[0], naoccpi_[h], 1.0, pa_tautau_[0], nsopi_[h]);
        }
        // Beta occupied
        if(nboccpi_[h] && nsopi_[h]){
            C_DGEMM('n', 'n', nsopi_[h], nboccpi_[h], nboccpi_[h], 1.0, pbOccC[0], nboccpi_[h],
                    TT_oo.matrix[h][0], nboccpi_[h], 0.0, temp[0], nsopi_[h]);
            C_DGEMM('n', 't', nsopi_[h], nsopi_[h], nboccpi_[h], 1.0, temp[0], nsopi_[h],
                    pbOccC[0], nboccpi_[h], 1.0, pb_tautau_[0], nsopi_[h]);
        }
        // Alpha virtual
        if(navirpi_[h] && nsopi_[h]){
            C_DGEMM('n', 'n', nsopi_[h], navirpi_[h], navirpi_[h], 1.0, paVirC[0], navirpi_[h],
                    TT_VV.matrix[h][0], navirpi_[h], 0.0, temp[0], nsopi_[h]);
            C_DGEMM('n', 't', nsopi_[h], nsopi_[h], navirpi_[h], 1.0, temp[0], nsopi_[h],
                    paVirC[0], navirpi_[h], 1.0, pa_tautau_[0], nsopi_[h]);
        }
        // Beta virtual
        if(nbvirpi_[h] && nsopi_[h]){
            C_DGEMM('n', 'n', nsopi_[h], nbvirpi_[h], nbvirpi_[h], 1.0, pbVirC[0], nbvirpi_[h],
                    TT_vv.matrix[h][0], nbvirpi_[h], 0.0, temp[0], nsopi_[h]);
            C_DGEMM('n', 't', nsopi_[h], nsopi_[h], nbvirpi_[h], 1.0, temp[0], nsopi_[h],
                    pbVirC[0], nbvirpi_[h], 1.0, pb_tautau_[0], nsopi_[h]);
        }

        free_block(temp);
    }

    dpd_file2_mat_wrt(&TT_OO);
    dpd_file2_mat_wrt(&TT_oo);
    dpd_file2_mat_wrt(&TT_VV);
    dpd_file2_mat_wrt(&TT_vv);

    dpd_file2_close(&T_OO);
    dpd_file2_close(&T_oo);
    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);

    dpd_file2_close(&TT_OO);
    dpd_file2_close(&TT_oo);
    dpd_file2_close(&TT_VV);
    dpd_file2_close(&TT_vv);
}

/**
 * Prints the occupation numbers from the OPDM
 */
void
DCFTSolver::print_opdm()
{
    dpdbuf4 L1, L2;
    dpdfile2 T_OO, T_oo, T_VV, T_vv;
    dpd_file2_init(&T_OO, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('O'), _ints->DPD_ID('O'), "Tau <O|O>");
    dpd_file2_init(&T_oo, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('o'), _ints->DPD_ID('o'), "Tau <o|o>");
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('V'), _ints->DPD_ID('V'), "Tau <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('v'), _ints->DPD_ID('v'), "Tau <v|v>");

    dpd_file2_mat_init(&T_OO);
    dpd_file2_mat_init(&T_oo);
    dpd_file2_mat_init(&T_VV);
    dpd_file2_mat_init(&T_vv);

    dpd_file2_mat_rd(&T_OO);
    dpd_file2_mat_rd(&T_oo);
    dpd_file2_mat_rd(&T_VV);
    dpd_file2_mat_rd(&T_vv);

    std::vector<std::pair<double, int> > aPairs;
    std::vector<std::pair<double, int> > bPairs;

    for(int h = 0; h < nirrep_; ++h){
        for(int row = 0; row < T_OO.params->coltot[h]; ++row)
            aPairs.push_back(std::make_pair(1.0 + T_OO.matrix[h][row][row], h));
        for(int row = 0; row < T_VV.params->coltot[h]; ++row)
            aPairs.push_back(std::make_pair(T_VV.matrix[h][row][row], h));
        for(int row = 0; row < T_oo.params->coltot[h]; ++row)
            bPairs.push_back(std::make_pair(1.0 + T_oo.matrix[h][row][row], h));
        for(int row = 0; row < T_vv.params->coltot[h]; ++row)
            bPairs.push_back(std::make_pair(T_vv.matrix[h][row][row], h));
    }
    dpd_file2_close(&T_OO);
    dpd_file2_close(&T_oo);
    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);

    sort(aPairs.begin(), aPairs.end(), greater<std::pair<double, int> >());
    sort(bPairs.begin(), bPairs.end(), greater<std::pair<double, int> >());

    int *aIrrepCount = init_int_array(nirrep_);
    int *bIrrepCount = init_int_array(nirrep_);
    char **irrepLabels = chkpt_->rd_irr_labs();

    fprintf(outfile, "\n\tOrbital occupations:\n\t\tAlpha occupied orbitals\n\t\t");
    for (int i = 0, count = 0; i < nalpha_; ++i, ++count) {
        int irrep = aPairs[i].second;
        fprintf(outfile, "%4d%-4s%11.4f  ", ++aIrrepCount[irrep], irrepLabels[irrep], aPairs[i].first);
        if (count % 4 == 3 && i != nalpha_)
            fprintf(outfile, "\n\t\t");
    }
    fprintf(outfile, "\n\n\t\tBeta occupied orbitals\n\t\t");
    for (int i = 0, count = 0; i < nbeta_; ++i, ++count) {
        int irrep = bPairs[i].second;
        fprintf(outfile, "%4d%-4s%11.4f  ", ++bIrrepCount[irrep], irrepLabels[irrep], bPairs[i].first);
        if (count % 4 == 3 && i != nbeta_)
            fprintf(outfile, "\n\t\t");
    }
    fprintf(outfile, "\n\n\t\tAlpha virtual orbitals\n\t\t");
    for (int i = nalpha_, count = 0; i < nmo_; ++i, ++count) {
        int irrep = aPairs[i].second;
        fprintf(outfile, "%4d%-4s%11.4f  ", ++aIrrepCount[irrep], irrepLabels[irrep], aPairs[i].first);
        if (count % 4 == 3 && i != nmo_)
            fprintf(outfile, "\n\t\t");
    }
    fprintf(outfile, "\n\n\t\tBeta virtual orbitals\n\t\t");
    for (int i = nbeta_, count = 0; i < nmo_; ++i, ++count) {
        int irrep = bPairs[i].second;
        fprintf(outfile, "%4d%-4s%11.4f  ", ++bIrrepCount[irrep], irrepLabels[irrep], bPairs[i].first);
        if (count % 4 == 3 && i != nmo_)
            fprintf(outfile, "\n\t\t");
    }
    fprintf(outfile, "\n\n");
    for (int h = 0; h < nirrep_; ++h)
        delete [] irrepLabels[h];
    delete[] irrepLabels;
    delete[] aIrrepCount;
    delete[] bIrrepCount;
}

}} // Namespaces


