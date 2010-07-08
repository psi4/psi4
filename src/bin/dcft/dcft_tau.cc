#include "dcft.h"
#include "libdpd/dpd.h"
#include "libqt/qt.h"
#include "libiwl/iwl.hpp"
#include "psifiles.h"
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

    for(int h = 0; h < _nIrreps; ++h){
        if(_soPI[h] == 0) continue;
        ::memset(_aTau[h][0], 0, _soPI[h]*_soPI[h]*sizeof(double));
        ::memset(_bTau[h][0], 0, _soPI[h]*_soPI[h]*sizeof(double));
        double **temp = block_matrix(_soPI[h], _soPI[h]);
        /*
         * Backtransform the Tau matrices to the AO basis
         * tauTrans = C moTau Ct
         */

        // Alpha occupied
        if(_nAOccPI[h] && _soPI[h]){
            C_DGEMM('n', 'n', _soPI[h], _nAOccPI[h], _nAOccPI[h], 1.0, _aOccC[h][0], _nAOccPI[h],
                    T_OO.matrix[h][0], _nAOccPI[h], 0.0, temp[0], _soPI[h]);
            C_DGEMM('n', 't', _soPI[h], _soPI[h], _nAOccPI[h], 1.0, temp[0], _soPI[h],
                    _aOccC[h][0], _nAOccPI[h], 1.0, _aTau[h][0], _soPI[h]);
        }
        // Beta occupied
        if(_nBOccPI[h] && _soPI[h]){
            C_DGEMM('n', 'n', _soPI[h], _nBOccPI[h], _nBOccPI[h], 1.0, _bOccC[h][0], _nBOccPI[h],
                    T_oo.matrix[h][0], _nBOccPI[h], 0.0, temp[0], _soPI[h]);
            C_DGEMM('n', 't', _soPI[h], _soPI[h], _nBOccPI[h], 1.0, temp[0], _soPI[h],
                    _bOccC[h][0], _nBOccPI[h], 1.0, _bTau[h][0], _soPI[h]);
        }
        // Alpha virtual
        if(_nAVirPI[h] && _soPI[h]){
            C_DGEMM('n', 'n', _soPI[h], _nAVirPI[h], _nAVirPI[h], 1.0, _aVirC[h][0], _nAVirPI[h],
                    T_VV.matrix[h][0], _nAVirPI[h], 0.0, temp[0], _soPI[h]);
            C_DGEMM('n', 't', _soPI[h], _soPI[h], _nAVirPI[h], 1.0, temp[0], _soPI[h],
                    _aVirC[h][0], _nAVirPI[h], 1.0, _aTau[h][0], _soPI[h]);
        }
        // Beta virtual
        if(_nBVirPI[h] && _soPI[h]){
            C_DGEMM('n', 'n', _soPI[h], _nBVirPI[h], _nBVirPI[h], 1.0, _bVirC[h][0], _nBVirPI[h],
                    T_vv.matrix[h][0], _nBVirPI[h], 0.0, temp[0], _soPI[h]);
            C_DGEMM('n', 't', _soPI[h], _soPI[h], _nBVirPI[h], 1.0, temp[0], _soPI[h],
                    _bVirC[h][0], _nBVirPI[h], 1.0, _bTau[h][0], _soPI[h]);
        }
        free_block(temp);
    }

    dpd_file2_close(&T_OO);
    dpd_file2_close(&T_oo);
    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);
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

    for(int h = 0; h < _nIrreps; ++h){
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

    int *aIrrepCount = init_int_array(_nIrreps);
    int *bIrrepCount = init_int_array(_nIrreps);
    char **irrepLabels = _chkpt->rd_irr_labs();

    fprintf(outfile, "\n\tOrbital occupations:\n\t\tAlpha occupied orbitals\n\t\t");
    for (int i = 0, count = 0; i < _nAOcc; ++i, ++count) {
        int irrep = aPairs[i].second;
        fprintf(outfile, "%4d%-4s%11.4f  ", ++aIrrepCount[irrep], irrepLabels[irrep], aPairs[i].first);
        if (count % 4 == 3 && i != _nAOcc)
            fprintf(outfile, "\n\t\t");
    }
    fprintf(outfile, "\n\n\t\tBeta occupied orbitals\n\t\t");
    for (int i = 0, count = 0; i < _nBOcc; ++i, ++count) {
        int irrep = bPairs[i].second;
        fprintf(outfile, "%4d%-4s%11.4f  ", ++bIrrepCount[irrep], irrepLabels[irrep], bPairs[i].first);
        if (count % 4 == 3 && i != _nBOcc)
            fprintf(outfile, "\n\t\t");
    }
    fprintf(outfile, "\n\n\t\tAlpha virtual orbitals\n\t\t");
    for (int i = _nAOcc, count = 0; i < _nMo; ++i, ++count) {
        int irrep = aPairs[i].second;
        fprintf(outfile, "%4d%-4s%11.4f  ", ++aIrrepCount[irrep], irrepLabels[irrep], aPairs[i].first);
        if (count % 4 == 3 && i != _nMo)
            fprintf(outfile, "\n\t\t");
    }
    fprintf(outfile, "\n\n\t\tBeta virtual orbitals\n\t\t");
    for (int i = _nBOcc, count = 0; i < _nMo; ++i, ++count) {
        int irrep = bPairs[i].second;
        fprintf(outfile, "%4d%-4s%11.4f  ", ++bIrrepCount[irrep], irrepLabels[irrep], bPairs[i].first);
        if (count % 4 == 3 && i != _nMo)
            fprintf(outfile, "\n\t\t");
    }
    fprintf(outfile, "\n\n");
    for (int h = 0; h < _nIrreps; ++h)
        delete [] irrepLabels[h];
    delete[] irrepLabels;
    delete[] aIrrepCount;
    delete[] bIrrepCount;
}

}} // Namespaces


