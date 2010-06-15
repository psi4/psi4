#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

/**
 * Checks the n-representability of the density matrix by checking the positive
 * semidefiniteness of various matrices, used as constraints by Mazziotti
 */
void
DCFTSolver::check_n_representability()
{
    // This shouldn't be used!  Just some experimentation...
    return;
    dpdbuf4 Laa, Lab, Lbb, D, Q, G;
    dpdfile2 T_OO, T_oo, T_VV, T_vv;
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
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

    for(int h = 0; h < _nIrreps; ++h){
        int nOcc = _nAOccPI[h];
        int nVir = _nAVirPI[h];
        unsigned long int dim = _moPI[h];

        dpd_buf4_mat_irrep_init(&Laa, h);
        dpd_buf4_mat_irrep_init(&Lab, h);
        dpd_buf4_mat_irrep_rd(&Laa, h);
        dpd_buf4_mat_irrep_rd(&Lab, h);

        /*
         * Alpha - Alpha and Beta - Beta are in the same matrix
         */
        double **OPDM = block_matrix(dim, dim);
        for(int i = 0; i < nOcc; ++i){
            for(int j = 0; j < nOcc; ++j){
                OPDM[i][j] = ( i == j ? 1.0 : 0.0) + T_OO.matrix[h][i][j];
            }
        }
        for(int a = nOcc; a < dim; ++a){
            for(int b = nOcc; b < dim; ++b){
                OPDM[a][b] =  T_VV.matrix[h][a-nOcc][b-nOcc];
            }
        }
        double **TPDMaa = block_matrix(dim*dim, dim*dim);
        double **TPDMab = block_matrix(dim*dim, dim*dim);
        // The OOVV and VVOO elements of the TPDM
        for(int ij = 0; ij < Laa.params->rowtot[h]; ++ij){
            int i = Laa.params->roworb[h][ij][0];
            int j = Laa.params->roworb[h][ij][1];
            unsigned long int IJ = i * dim + j;
            for(int ab = 0; ab < Laa.params->coltot[h]; ++ab){
                int a = Laa.params->colorb[h][ab][0] + nOcc;
                int b = Laa.params->colorb[h][ab][1] + nOcc;
                unsigned long int AB = a * dim + b;
                TPDMaa[AB][IJ] = TPDMaa[IJ][AB] = Laa.matrix[h][ij][ab];
            }
        }
        // The OOOO elements of the TPDM
        for(int ij = 0; ij < Laa.params->rowtot[h]; ++ij){
            int i = Laa.params->roworb[h][ij][0];
            int j = Laa.params->roworb[h][ij][1];
            unsigned long int IJ = i * dim + j;
            for(int kl = 0; kl < Laa.params->rowtot[h]; ++kl){
                int k = Laa.params->roworb[h][kl][0];
                int l = Laa.params->roworb[h][kl][1];
                unsigned long int KL = k * dim + l;
                double AAval = 0.0;
                double ABval = 0.0;
                for(int ab = 0; ab < Laa.params->coltot[h]; ++ab){
                    AAval += Laa.matrix[h][ij][ab] * Laa.matrix[h][kl][ab];
                    ABval += Lab.matrix[h][ij][ab] * Lab.matrix[h][kl][ab];
                }
                TPDMaa[IJ][KL] = 0.5 * AAval;
                TPDMab[IJ][KL] = 0.5 * ABval;
            }
        }
        // The VVVV elements of the TPDM
        for(int ab = 0; ab < Laa.params->coltot[h]; ++ab){
            int a = Laa.params->colorb[h][ab][0] + nOcc;
            int b = Laa.params->colorb[h][ab][1] + nOcc;
            unsigned long int AB = a * dim + b;
            for(int cd = 0; cd < Laa.params->coltot[h]; ++cd){
                int c = Laa.params->colorb[h][cd][0] + nOcc;
                int d = Laa.params->colorb[h][cd][1] + nOcc;
                unsigned long int CD = c * dim + d;
                double AAval = 0.0;
                double ABval = 0.0;
                for(int ij = 0; ij < Laa.params->rowtot[h]; ++ij){
                    AAval += Laa.matrix[h][ij][ab] * Laa.matrix[h][ij][cd];
                    ABval += Lab.matrix[h][ij][ab] * Lab.matrix[h][ij][cd];
                }
                TPDMaa[AB][CD] = 0.5 * AAval;
                TPDMab[AB][CD] = 0.5 * ABval;
            }
        }
        // The OVOV elements
        for(int i = 0; i < nOcc; ++i){
            int I = i;
            for(int a = 0; a < nVir; ++a){
                int A = a + nOcc;
                unsigned long int IA = I * dim + A;
                unsigned long int AI = A * dim + I;
                for(int j = 0; j < nOcc; ++j){
                    int J = j;
                    for(int b = 0; b < nVir; ++b){
                        int B = b + nOcc;
                        unsigned long int JB = J * dim + B;
                        unsigned long int BJ = B * dim + J;
                        double AAval = 0.0;
                        double ABval = 0.0;
                        for(int k = 0; k < nOcc; ++k){
                            unsigned long int ik = i * nOcc + k;
                            unsigned long int jk = j * nOcc + k;
                            for(int c = 0; c < nVir; ++c){
                                unsigned long int ac = a * nVir + c;
                                unsigned long int bc = b * nVir + c;
                                AAval += Laa.matrix[h][ik][bc] * Laa.matrix[h][jk][ac];
                                AAval += Lab.matrix[h][ik][bc] * Lab.matrix[h][jk][ac];
                                ABval += Laa.matrix[h][ik][bc] * Laa.matrix[h][jk][ac];
                                ABval += Lab.matrix[h][ik][bc] * Lab.matrix[h][jk][ac];
                            }
                        }
                        TPDMaa[IA][JB] = TPDMaa[AI][BJ] = -AAval;
                        TPDMaa[IA][BJ] = TPDMaa[AI][JB] = AAval;
                        TPDMab[IA][JB] = TPDMab[AI][BJ] = -ABval;
                        TPDMab[IA][BJ] = TPDMab[AI][JB] = ABval;
                    }
                }
            }
        }

        // Contract the integra

        print_mat(TPDMab, dim*dim, dim*dim, outfile);
        free_block(OPDM);
        free_block(TPDMaa);
        free_block(TPDMab);
    }
}

}} // Namespaces
