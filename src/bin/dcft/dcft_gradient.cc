#include <libtrans/integraltransform.h>
#include <libpsio/psio.hpp>
#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

void
DCFTSolver::compute_gradient()
{
    bool orbResponseDone    = false;
    bool lambdaResponseDone = false;

    // Print out the header
    fprintf(outfile,"\n\n\t\t*******************************************\n");
    fprintf(outfile,    "\t\t*       DCFT Analytic Gradients Code      *\n");
    fprintf(outfile,    "\t\t*   by A.Yu. Sokolov and A.C. Simmonett   *\n");
    fprintf(outfile,    "\t\t*******************************************\n\n");

    orbital_response_guess();
//    cumulant_response_guess();

    int cycle = 0;

}

void
DCFTSolver::orbital_response_guess()
{

    dpdbuf4 L;
    dpdfile2 T;

    // Copy the converged cumulant as a guess for the cumulant response

    // Z_IJAB = L_IJAB
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_copy(&L,PSIF_DCFT_DPD,"Z <OO|VV>");
    dpd_buf4_close(&L);

    // Z_IjAb = L_IjAb
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_copy(&L,PSIF_DCFT_DPD,"Z <Oo|Vv>");
    dpd_buf4_close(&L);

    // Z_ijab = L_ijab
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    dpd_buf4_copy(&L,PSIF_DCFT_DPD,"Z <oo|vv>");
    dpd_buf4_close(&L);

    // Copy the latest tau as a guess for relaxed tau

    dpd_file2_init(&T, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('O'), _ints->DPD_ID('O'), "Tau <O|O>");
    dpd_file2_copy(&T,PSIF_DCFT_DPD,"pTau <O|O>");
    dpd_file2_close(&T);

    dpd_file2_init(&T, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('o'), _ints->DPD_ID('o'), "Tau <o|o>");
    dpd_file2_copy(&T,PSIF_DCFT_DPD,"pTau <o|o>");
    dpd_file2_close(&T);

    dpd_file2_init(&T, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('V'), _ints->DPD_ID('V'), "Tau <V|V>");
    dpd_file2_copy(&T,PSIF_DCFT_DPD,"pTau <V|V>");
    dpd_file2_close(&T);

    dpd_file2_init(&T, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('v'), _ints->DPD_ID('v'), "Tau <v|v>");
    dpd_file2_copy(&T,PSIF_DCFT_DPD,"pTau <v|v>");
    dpd_file2_close(&T);

    // Compute the generalized densities for the MO Lagrangian
    compute_density();

}

void
DCFTSolver::compute_density()
{

    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 Zaa, Zab, Zbb, Laa, Lab, Lbb, Gaa, Gab, Gba, Gbb, I, L;
    dpdfile2 pT_OO, pT_oo, pT_VV, pT_vv, T_OO, T_oo, T_VV, T_vv;

    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Z <OO|VV>");
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Z <oo|vv>");
    dpd_file2_init(&T_OO, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('O'), _ints->DPD_ID('O'), "Tau <O|O>");
    dpd_file2_init(&T_oo, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('o'), _ints->DPD_ID('o'), "Tau <o|o>");
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('V'), _ints->DPD_ID('V'), "Tau <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('v'), _ints->DPD_ID('v'), "Tau <v|v>");
    dpd_file2_init(&pT_OO, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('O'), _ints->DPD_ID('O'), "pTau <O|O>");
    dpd_file2_init(&pT_oo, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('o'), _ints->DPD_ID('o'), "pTau <o|o>");
    dpd_file2_init(&pT_VV, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('V'), _ints->DPD_ID('V'), "pTau <V|V>");
    dpd_file2_init(&pT_vv, PSIF_DCFT_DPD, 0,
                  _ints->DPD_ID('v'), _ints->DPD_ID('v'), "pTau <v|v>");

    dpd_file2_mat_init(&T_OO);
    dpd_file2_mat_init(&T_oo);
    dpd_file2_mat_init(&T_VV);
    dpd_file2_mat_init(&T_vv);
    dpd_file2_mat_init(&pT_OO);
    dpd_file2_mat_init(&pT_oo);
    dpd_file2_mat_init(&pT_VV);
    dpd_file2_mat_init(&pT_vv);

    dpd_file2_mat_rd(&T_OO);
    dpd_file2_mat_rd(&T_oo);
    dpd_file2_mat_rd(&T_VV);
    dpd_file2_mat_rd(&T_vv);
    dpd_file2_mat_rd(&pT_OO);
    dpd_file2_mat_rd(&pT_oo);
    dpd_file2_mat_rd(&pT_VV);
    dpd_file2_mat_rd(&pT_vv);

    Matrix aKappa(nirrep_, naoccpi_, naoccpi_);
    Matrix bKappa(nirrep_, nboccpi_, nboccpi_);
    Matrix aOccTau(nirrep_, naoccpi_, naoccpi_);
    Matrix bOccTau(nirrep_, nboccpi_, nboccpi_);
    Matrix aVirTau(nirrep_, navirpi_, navirpi_);
    Matrix bVirTau(nirrep_, nbvirpi_, nbvirpi_);
    Matrix aOccPTau(nirrep_, naoccpi_, naoccpi_);
    Matrix bOccPTau(nirrep_, nboccpi_, nboccpi_);
    Matrix aVirPTau(nirrep_, navirpi_, navirpi_);
    Matrix bVirPTau(nirrep_, nbvirpi_, nbvirpi_);

    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int j = 0; j < naoccpi_[h]; ++j){
                aKappa.set(h, i, j, (i == j ? 1.0 : 0.0));
                aOccTau.set(h, i, j, T_OO.matrix[h][i][j]);
                aOccPTau.set(h, i, j, pT_OO.matrix[h][i][j]);
            }
        }
        for(int a = 0; a < navirpi_[h]; ++a){
            for(int b = 0; b < navirpi_[h]; ++b){
                aVirTau.set(h, a, b, T_VV.matrix[h][a][b]);
                aVirPTau.set(h, a, b, pT_VV.matrix[h][a][b]);
            }
        }
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int j = 0; j < nboccpi_[h]; ++j){
                bKappa.set(h, i, j, (i == j ? 1.0 : 0.0));
                bOccTau.set(h, i, j, T_oo.matrix[h][i][j]);
                bOccPTau.set(h, i, j, pT_oo.matrix[h][i][j]);
            }
        }
        for(int a = 0; a < nbvirpi_[h]; ++a){
            for(int b = 0; b < nbvirpi_[h]; ++b){
                bVirTau.set(h, a, b, T_vv.matrix[h][a][b]);
                bVirPTau.set(h, a, b, pT_vv.matrix[h][a][b]);
            }
        }
    }

    dpd_file2_close(&T_OO);
    dpd_file2_close(&T_oo);
    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);
    dpd_file2_close(&pT_OO);
    dpd_file2_close(&pT_oo);
    dpd_file2_close(&pT_VV);
    dpd_file2_close(&pT_vv);

    /*
     * The VVVV block
     */
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
              ID("[V,V]"), ID("[V,V]"), 0, "Gamma <VV|VV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gaa, h);
        dpd_buf4_mat_irrep_init(&Laa, h);
        dpd_buf4_mat_irrep_init(&Zaa, h);
        dpd_buf4_mat_irrep_rd(&Laa, h);
        dpd_buf4_mat_irrep_rd(&Zaa, h);
        for(size_t ab = 0; ab < Gaa.params->rowtot[h]; ++ab){
            size_t a = Gaa.params->roworb[h][ab][0];
            int Ga = Gaa.params->psym[a];
            a -= Gaa.params->poff[Ga];
            size_t b = Gaa.params->roworb[h][ab][1];
            int Gb = Gaa.params->qsym[b];
            b -= Gaa.params->qoff[Gb];
            for(size_t cd = 0; cd < Gaa.params->coltot[h]; ++cd){
                double tpdm = 0.0;
                for(size_t ij = 0; ij < Laa.params->rowtot[h]; ++ij){
                    tpdm += 0.0625 * (Zaa.matrix[h][ij][ab] * Laa.matrix[h][ij][cd] + Laa.matrix[h][ij][ab] * Zaa.matrix[h][ij][cd]);
                }
                size_t c = Gaa.params->colorb[h][cd][0];
                int Gc = Gaa.params->rsym[c];
                c -= Gaa.params->roff[Gc];
                size_t d = Gaa.params->colorb[h][cd][1];
                int Gd = Gaa.params->ssym[d];
                d -= Gaa.params->soff[Gd];
                if(Ga == Gc && Gb == Gd) tpdm += 0.25 * aVirTau(Ga, a, c) * aVirPTau(Gb, b, d);
                if(Ga == Gd && Gb == Gc) tpdm -= 0.25 * aVirTau(Ga, a, d) * aVirPTau(Gb, b, c);
                if(Gb == Gc && Ga == Gd) tpdm -= 0.25 * aVirTau(Gb, b, c) * aVirPTau(Ga, a, d);
                if(Ga == Gc && Gb == Gd) tpdm += 0.25 * aVirTau(Gb, b, d) * aVirPTau(Ga, a, c);

                if(Ga == Gc && Gb == Gd) tpdm -= 0.25 * aVirTau(Ga, a, c) * aVirTau(Gb, b, d);
                if(Ga == Gd && Gb == Gc) tpdm += 0.25 * aVirTau(Ga, a, d) * aVirTau(Gb, b, c);

                Gaa.matrix[h][ab][cd] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Laa, h);
        dpd_buf4_mat_irrep_close(&Zaa, h);
    }

    dpd_buf4_print(&Gaa,outfile,1);

    dpd_buf4_close(&Gaa);

    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
              ID("[V,v]"), ID("[V,v]"), 0, "Gamma <Vv|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gab, h);
        dpd_buf4_mat_irrep_init(&Lab, h);
        dpd_buf4_mat_irrep_init(&Zab, h);
        dpd_buf4_mat_irrep_rd(&Lab, h);
        dpd_buf4_mat_irrep_rd(&Zab, h);
        for(size_t ab = 0; ab < Gab.params->rowtot[h]; ++ab){
            size_t a = Gab.params->roworb[h][ab][0];
            int Ga = Gab.params->psym[a];
            a -= Gab.params->poff[Ga];
            size_t b = Gab.params->roworb[h][ab][1];
            int Gb = Gab.params->qsym[b];
            b -= Gab.params->qoff[Gb];
            for(size_t cd = 0; cd < Gab.params->coltot[h]; ++cd){
                double tpdm = 0.0;
                for(size_t ij = 0; ij < Lab.params->rowtot[h]; ++ij){
                    tpdm += 0.125 * (Zab.matrix[h][ij][ab] * Lab.matrix[h][ij][cd] + Lab.matrix[h][ij][ab] * Zab.matrix[h][ij][cd]);
                }
                size_t c = Gab.params->colorb[h][cd][0];
                int Gc = Gab.params->rsym[c];
                c -= Gab.params->roff[Gc];
                size_t d = Gab.params->colorb[h][cd][1];
                int Gd = Gab.params->ssym[d];
                d -= Gab.params->soff[Gd];
                if(Ga == Gc && Gb == Gd) tpdm += 0.25 * aVirTau(Ga, a, c) * bVirPTau(Gb, b, d);
                if(Ga == Gc && Gb == Gd) tpdm += 0.25 * bVirTau(Gb, b, d) * aVirPTau(Ga, a, c);

                if(Ga == Gc && Gb == Gd) tpdm -= 0.25 * aVirTau(Ga, a, c) * bVirTau(Gb, b, d);
                Gab.matrix[h][ab][cd] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gab, h);
        dpd_buf4_mat_irrep_close(&Gab, h);
        dpd_buf4_mat_irrep_close(&Lab, h);
        dpd_buf4_mat_irrep_close(&Zab, h);
    }

    dpd_buf4_print(&Gab,outfile,1);

    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
              ID("[v,v]"), ID("[v,v]"), 0, "Gamma <vv|vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gbb, h);
        dpd_buf4_mat_irrep_init(&Lbb, h);
        dpd_buf4_mat_irrep_init(&Zbb, h);
        dpd_buf4_mat_irrep_rd(&Lbb, h);
        dpd_buf4_mat_irrep_rd(&Zbb, h);

        for(size_t ab = 0; ab < Gbb.params->rowtot[h]; ++ab){
            size_t a = Gbb.params->roworb[h][ab][0];
            int Ga = Gbb.params->psym[a];
            a -= Gbb.params->poff[Ga];
            size_t b = Gbb.params->roworb[h][ab][1];
            int Gb = Gbb.params->qsym[b];
            b -= Gbb.params->qoff[Gb];
            for(size_t cd = 0; cd < Gbb.params->coltot[h]; ++cd){
                double tpdm = 0.0;
                for(size_t ij = 0; ij < Lbb.params->rowtot[h]; ++ij){
                    tpdm += 0.0625 * (Zbb.matrix[h][ij][ab] * Lbb.matrix[h][ij][cd] + Lbb.matrix[h][ij][ab] * Zbb.matrix[h][ij][cd]);
                }
                size_t c = Gbb.params->colorb[h][cd][0];
                int Gc = Gbb.params->rsym[c];
                c -= Gbb.params->roff[Gc];
                size_t d = Gbb.params->colorb[h][cd][1];
                int Gd = Gbb.params->ssym[d];
                d -= Gbb.params->soff[Gd];
                if(Ga == Gc && Gb == Gd) tpdm += 0.25 * bVirTau(Ga, a, c) * bVirPTau(Gb, b, d);
                if(Ga == Gd && Gb == Gc) tpdm -= 0.25 * bVirTau(Ga, a, d) * bVirPTau(Gb, b, c);
                if(Gb == Gc && Ga == Gd) tpdm -= 0.25 * bVirTau(Gb, b, c) * bVirPTau(Ga, a, d);
                if(Ga == Gc && Gb == Gd) tpdm += 0.25 * bVirTau(Gb, b, d) * bVirPTau(Ga, a, c);

                if(Ga == Gc && Gb == Gd) tpdm -= 0.25 * bVirTau(Ga, a, c) * bVirTau(Gb, b, d);
                if(Ga == Gd && Gb == Gc) tpdm += 0.25 * bVirTau(Ga, a, d) * bVirTau(Gb, b, c);
                Gbb.matrix[h][ab][cd] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Lbb, h);
        dpd_buf4_mat_irrep_close(&Zbb, h);
    }

    dpd_buf4_print(&Gbb,outfile,1);

    dpd_buf4_close(&Gbb);

    /*
     * The OOOO  block
     */
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
              ID("[O,O]"), ID("[O,O]"), 0, "Gamma <OO|OO>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gaa, h);
        dpd_buf4_mat_irrep_init(&Laa, h);
        dpd_buf4_mat_irrep_init(&Zaa, h);
        dpd_buf4_mat_irrep_rd(&Laa, h);
        dpd_buf4_mat_irrep_rd(&Zaa, h);
                for(size_t ij = 0; ij < Gaa.params->rowtot[h]; ++ij){
            size_t i = Gaa.params->roworb[h][ij][0];
            int Gi = Gaa.params->psym[i];
            i -= Gaa.params->poff[Gi];
            size_t j = Gaa.params->roworb[h][ij][1];
            int Gj = Gaa.params->qsym[j];
            j -= Gaa.params->qoff[Gj];
            for(size_t kl = 0; kl < Gaa.params->coltot[h]; ++kl){
                double tpdm = 0.0;
                for(size_t ab = 0; ab < Laa.params->coltot[h]; ++ab){
                    tpdm += 0.0625 * (Zaa.matrix[h][ij][ab] * Laa.matrix[h][kl][ab] + Laa.matrix[h][ij][ab] * Zaa.matrix[h][kl][ab]);
                }
                size_t k = Gaa.params->colorb[h][kl][0];
                int Gk = Gaa.params->rsym[k];
                k -= Gaa.params->roff[Gk];
                size_t l = Gaa.params->colorb[h][kl][1];
                int Gl = Gaa.params->ssym[l];
                l -= Gaa.params->soff[Gl];

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * aKappa(Gi, i, k) * aKappa(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * aKappa(Gi, i, l) * aKappa(Gj, j, k);

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * (aKappa(Gi, i, k) + aOccTau(Gi, i, k)) * aOccPTau(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * (aKappa(Gi, i, l) + aOccTau(Gi, i, l)) * aOccPTau(Gj, j, k);
                if(Gj == Gk && Gi == Gl) tpdm -= 0.25 * (aKappa(Gj, j, k) + aOccTau(Gj, j, k)) * aOccPTau(Gi, i, l);
                if(Gj == Gl && Gi == Gk) tpdm += 0.25 * (aKappa(Gj, j, l) + aOccTau(Gj, j, l)) * aOccPTau(Gi, i, k);

                if(Gi == Gk && Gj == Gl) tpdm -= 0.25 * aOccTau(Gi, i, k) * aOccTau(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm += 0.25 * aOccTau(Gi, i, l) * aOccTau(Gj, j, k);

                Gaa.matrix[h][ij][kl] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Laa, h);
        dpd_buf4_mat_irrep_close(&Zaa, h);
    }

    dpd_buf4_print(&Gaa,outfile,1);

    dpd_buf4_close(&Gaa);


    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
              ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gab, h);
        dpd_buf4_mat_irrep_init(&Lab, h);
        dpd_buf4_mat_irrep_init(&Zab, h);
        dpd_buf4_mat_irrep_rd(&Lab, h);
        dpd_buf4_mat_irrep_rd(&Zab, h);
        for(size_t ij = 0; ij < Gab.params->rowtot[h]; ++ij){
            size_t i = Gab.params->roworb[h][ij][0];
            int Gi = Gab.params->psym[i];
            i -= Gab.params->poff[Gi];
            size_t j = Gab.params->roworb[h][ij][1];
            int Gj = Gab.params->qsym[j];
            j -= Gab.params->qoff[Gj];
            for(size_t kl = 0; kl < Gab.params->coltot[h]; ++kl){
                double tpdm = 0.0;
                for(size_t ab = 0; ab < Lab.params->coltot[h]; ++ab){
                    tpdm += 0.125 * (Zab.matrix[h][ij][ab] * Lab.matrix[h][kl][ab] + Lab.matrix[h][ij][ab] * Zab.matrix[h][kl][ab]);
                }
                size_t k = Gab.params->colorb[h][kl][0];
                int Gk = Gab.params->rsym[k];
                k -= Gab.params->roff[Gk];
                size_t l = Gab.params->colorb[h][kl][1];
                int Gl = Gab.params->ssym[l];
                l -= Gab.params->soff[Gl];
                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * aKappa(Gi, i, k) * bKappa(Gj, j, l);

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * (aKappa(Gi, i, k) + aOccTau(Gi, i, k)) * bOccPTau(Gj, j, l);
                if(Gj == Gl && Gi == Gk) tpdm += 0.25 * (bKappa(Gj, j, l) + bOccTau(Gj, j, l)) * aOccPTau(Gi, i, k);

                if(Gi == Gk && Gj == Gl) tpdm -= 0.25 * aOccTau(Gi, i, k) * bOccTau(Gj, j, l);

                Gab.matrix[h][ij][kl] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gab, h);
        dpd_buf4_mat_irrep_close(&Gab, h);
        dpd_buf4_mat_irrep_close(&Lab, h);
        dpd_buf4_mat_irrep_close(&Zab, h);
    }

    dpd_buf4_print(&Gab,outfile,1);

    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
              ID("[o,o]"), ID("[o,o]"), 0, "Gamma <oo|oo>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gbb, h);
        dpd_buf4_mat_irrep_init(&Lbb, h);
        dpd_buf4_mat_irrep_init(&Zbb, h);
        dpd_buf4_mat_irrep_rd(&Lbb, h);
        dpd_buf4_mat_irrep_rd(&Zbb, h);
        for(size_t ij = 0; ij < Gbb.params->rowtot[h]; ++ij){
            size_t i = Gbb.params->roworb[h][ij][0];
            int Gi = Gbb.params->psym[i];
            i -= Gbb.params->poff[Gi];
            size_t j = Gbb.params->roworb[h][ij][1];
            int Gj = Gbb.params->qsym[j];
            j -= Gbb.params->qoff[Gj];
            for(size_t kl = 0; kl < Gbb.params->coltot[h]; ++kl){
                double tpdm = 0.0;
                for(size_t ab = 0; ab < Lbb.params->coltot[h]; ++ab){
                    tpdm += 0.0625 * (Zbb.matrix[h][ij][ab] * Lbb.matrix[h][kl][ab] + Lbb.matrix[h][ij][ab] * Zbb.matrix[h][kl][ab]);
                }
                size_t k = Gbb.params->colorb[h][kl][0];
                int Gk = Gbb.params->rsym[k];
                k -= Gbb.params->roff[Gk];
                size_t l = Gbb.params->colorb[h][kl][1];
                int Gl = Gbb.params->ssym[l];
                l -= Gbb.params->soff[Gl];
                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * bKappa(Gi, i, k) * bKappa(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * bKappa(Gi, i, l) * bKappa(Gj, j, k);

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * (bKappa(Gi, i, k) + bOccTau(Gi, i, k)) * bOccPTau(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * (bKappa(Gi, i, l) + bOccTau(Gi, i, l)) * bOccPTau(Gj, j, k);
                if(Gj == Gk && Gi == Gl) tpdm -= 0.25 * (bKappa(Gj, j, k) + bOccTau(Gj, j, k)) * bOccPTau(Gi, i, l);
                if(Gj == Gl && Gi == Gk) tpdm += 0.25 * (bKappa(Gj, j, l) + bOccTau(Gj, j, l)) * bOccPTau(Gi, i, k);

                if(Gi == Gk && Gj == Gl) tpdm -= 0.25 * bOccTau(Gi, i, k) * bOccTau(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm += 0.25 * bOccTau(Gi, i, l) * bOccTau(Gj, j, k);

                Gbb.matrix[h][ij][kl] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Lbb, h);
        dpd_buf4_mat_irrep_close(&Zbb, h);
    }

    dpd_buf4_print(&Gbb,outfile,1);

    dpd_buf4_close(&Gbb);

    /*
     * The OOVV and VVOO blocks
     */

    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
              ID("[O,O]"), ID("[V,V]"), 0, "Gamma <OO|VV>");
    dpd_buf4_axpbycz(&Laa,&Zaa,&Gaa,0.25,0.25,1.0);

    dpd_buf4_print(&Gaa,outfile,1);

    dpd_buf4_close(&Gaa);

    dpd_buf4_init(&Gab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
              ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");
    dpd_buf4_axpbycz(&Lab,&Zab,&Gab,0.25,0.25,1.0);

    dpd_buf4_print(&Gab,outfile,1);

    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
              ID("[o,o]"), ID("[v,v]"), 0, "Gamma <oo|vv>");
    dpd_buf4_axpbycz(&Lbb,&Zbb,&Gbb,0.25,0.25,1.0);

    dpd_buf4_print(&Gbb,outfile,1);

    dpd_buf4_close(&Gbb);

//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
//                      ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
//    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[V,v]"), ID("[O,o]"), "MO Ints <Vv|Oo>");
//    dpd_buf4_close(&I);

    dpd_buf4_close(&Laa);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Lbb);
    dpd_buf4_close(&Zaa);
    dpd_buf4_close(&Zab);
    dpd_buf4_close(&Zbb);

    /*
     * The OVOV block
     */

    // There are five unique spin cases: Г<IAJB>, Г<iajb>, Г<IaJb>, Г<iAjB>, Г<IajB>
    dpdbuf4 Laaaa, Laabb, Labba, Lbaab, Lbbbb, Taa, Tab, Tba, Tbb, Gaa2, Laa2, Lab2, Lbb2;

    // TEMPORARY: Sort the cumulant Z-vector elements to chemist's notation.
    // MOVE THIS TO THE Z-VECTOR UPDATES WHEN NEEDED!!!
    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Z <OO|VV>");
    dpd_buf4_sort(&Zaa, PSIF_DCFT_DPD, prqs, ID("[O,V]"), ID("[O,V]"), "Z (OV|OV)");
    dpd_buf4_close(&Zaa);
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Z <Oo|Vv>");
    dpd_buf4_sort(&Zab, PSIF_DCFT_DPD, psqr, ID("[O,v]"), ID("[o,V]"), "Z (Ov|oV)");
    dpd_buf4_close(&Zab);
    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Z <oo|vv>");
    dpd_buf4_sort(&Zbb, PSIF_DCFT_DPD, prqs, ID("[o,v]"),ID("[o,v]"), "Z (ov|ov)");
    dpd_buf4_close(&Zbb);

    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Z (Ov|oV)");

    dpd_buf4_sort(&Zab, PSIF_DCFT_DPD, psrq, ID("[O,V]"),ID("[o,v]"), "Z (OV|ov)");
    dpd_buf4_close(&Zab);

    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Z (OV|OV)");
    dpd_contract444(&Laa, &Zaa, &Gaa, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Laa);
    dpd_buf4_close(&Zaa);
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Z (OV|ov)");
    dpd_contract444(&Lab, &Zab, &Gaa, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Zab);
    dpd_buf4_symm(&Gaa);
    dpd_buf4_close(&Gaa);

    // Г<IAJB> spin case

    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gaa, h);
        dpd_buf4_mat_irrep_rd(&Gaa, h);
        for(size_t ia = 0; ia < Gaa.params->rowtot[h]; ++ia){
            size_t i = Gaa.params->roworb[h][ia][0];
            int Gi = Gaa.params->psym[i];
            i -= Gaa.params->poff[Gi];
            size_t a = Gaa.params->roworb[h][ia][1];
            int Ga = Gaa.params->qsym[a];
            a -= Gaa.params->qoff[Ga];
            for(size_t jb = 0; jb < Gaa.params->coltot[h]; ++jb){
                size_t j = Gaa.params->colorb[h][jb][0];
                int Gj = Gaa.params->rsym[j];
                j -= Gaa.params->roff[Gj];
                size_t b = Gaa.params->colorb[h][jb][1];
                int Gb = Gaa.params->ssym[b];
                b -= Gaa.params->soff[Gb];
                if(Gi == Gj && Ga == Gb) {
                    Gaa.matrix[h][ia][jb] += (aKappa(Gi, i, j) + aOccTau(Gi, i, j)) * aVirPTau(Ga, a, b);
                    Gaa.matrix[h][ia][jb] += aVirTau(Ga, a, b) * (aOccPTau(Gi, i, j) - aOccTau(Gi, i, j));
                }
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Gaa, h);
    }

    dpd_buf4_print(&Gaa,outfile,1);

    dpd_buf4_close(&Gaa);

    // Г<IaJb> and Г<iAjB> spin cases:

    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Lambda (Ov|oV)");
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Z (Ov|oV)");
    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
    dpd_contract444(&Lab, &Zab, &Gab, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Gab);
    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
    dpd_buf4_symm(&Gab);
    dpd_buf4_close(&Gab);
    dpd_buf4_init(&Gba, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");
    dpd_contract444(&Lab, &Zab, &Gba, 1, 1, -1.0, 0.0);
    dpd_buf4_close(&Gba);
    dpd_buf4_init(&Gba, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");
    dpd_buf4_symm(&Gba);
    dpd_buf4_close(&Gba);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Zab);

    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gab, h);
        dpd_buf4_mat_irrep_rd(&Gab, h);
        for(size_t ia = 0; ia < Gab.params->rowtot[h]; ++ia){
            size_t i = Gab.params->roworb[h][ia][0];
            int Gi = Gab.params->psym[i];
            i -= Gab.params->poff[Gi];
            size_t a = Gab.params->roworb[h][ia][1];
            int Ga = Gab.params->qsym[a];
            a -= Gab.params->qoff[Ga];
            for(size_t jb = 0; jb < Gab.params->coltot[h]; ++jb){
                size_t j = Gab.params->colorb[h][jb][0];
                int Gj = Gab.params->rsym[j];
                j -= Gab.params->roff[Gj];
                size_t b = Gab.params->colorb[h][jb][1];
                int Gb = Gab.params->ssym[b];
                b -= Gab.params->soff[Gb];
                if(Gi == Gj && Ga == Gb) {
                    Gab.matrix[h][ia][jb] += (aKappa(Gi, i, j) + aOccTau(Gi, i, j)) * bVirPTau(Ga, a, b);
                    Gab.matrix[h][ia][jb] += bVirTau(Ga, a, b) * (aOccPTau(Gi, i, j) - aOccTau(Gi, i, j));
                }
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gab, h);
        dpd_buf4_mat_irrep_close(&Gab, h);
    }

    dpd_buf4_print(&Gab,outfile,1);

    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Gba, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gba, h);
        dpd_buf4_mat_irrep_rd(&Gba, h);
        for(size_t ia = 0; ia < Gba.params->rowtot[h]; ++ia){
            size_t i = Gba.params->roworb[h][ia][0];
            int Gi = Gba.params->psym[i];
            i -= Gba.params->poff[Gi];
            size_t a = Gba.params->roworb[h][ia][1];
            int Ga = Gba.params->qsym[a];
            a -= Gba.params->qoff[Ga];
            for(size_t jb = 0; jb < Gba.params->coltot[h]; ++jb){
                size_t j = Gba.params->colorb[h][jb][0];
                int Gj = Gba.params->rsym[j];
                j -= Gba.params->roff[Gj];
                size_t b = Gba.params->colorb[h][jb][1];
                int Gb = Gba.params->ssym[b];
                b -= Gba.params->soff[Gb];
                if(Gi == Gj && Ga == Gb) {
                    Gba.matrix[h][ia][jb] += (bKappa(Gi, i, j) + bOccTau(Gi, i, j)) * aVirPTau(Ga, a, b);
                    Gba.matrix[h][ia][jb] += aVirTau(Ga, a, b) * (bOccPTau(Gi, i, j) - bOccTau(Gi, i, j));
                }

            }
        }
        dpd_buf4_mat_irrep_wrt(&Gba, h);
        dpd_buf4_mat_irrep_close(&Gba, h);
    }

    dpd_buf4_print(&Gba,outfile,1);

    dpd_buf4_close(&Gba);

    // Г<IajB> spin case:

    // <Ov|oV> and <oV|Ov> We will need to double it later. Check the -0.5 factor in front!!!
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Z (OV|ov)");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    dpd_buf4_init(&Zaa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Z (OV|OV)");
    dpd_contract444(&Laa, &Zab, &Tab, 0, 1, -0.5, 0.0);
    dpd_contract444(&Zaa, &Lab, &Tab, 0, 1, -0.5, 1.0);
    dpd_buf4_close(&Laa);
    dpd_buf4_close(&Zaa);
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Z (ov|ov)");
    dpd_contract444(&Lab, &Zbb, &Tab, 0, 1, -0.5, 1.0);
    dpd_contract444(&Zab, &Lbb, &Tab, 0, 1, -0.5, 1.0);
    dpd_buf4_close(&Lbb);
    dpd_buf4_close(&Zbb);
    dpd_buf4_close(&Tab);
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");
    dpd_buf4_sort(&Tab, PSIF_DCFT_DENSITY, psrq, ID("[O,v]"), ID("[o,V]"), "Gamma <Ov|oV>");
    dpd_buf4_close(&Tab);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Zab);

    dpd_buf4_init(&Gba, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");
    dpd_buf4_print(&Gba,outfile,1);
    dpd_buf4_close(&Gba);

    // Г<iajb> spin case:

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    dpd_buf4_init(&Zbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Z (ov|ov)");
    dpd_contract444(&Lbb, &Zbb, &Gbb, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Lbb);
    dpd_buf4_close(&Zbb);
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_buf4_init(&Zab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Z (OV|ov)");
    dpd_contract444(&Lab, &Zab, &Gbb, 1, 1, -1.0, 1.0);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Zab);

    dpd_buf4_symm(&Gbb);
    dpd_buf4_close(&Gbb);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gbb, h);
        dpd_buf4_mat_irrep_rd(&Gbb, h);
        for(size_t ia = 0; ia < Gbb.params->rowtot[h]; ++ia){
            size_t i = Gbb.params->roworb[h][ia][0];
            int Gi = Gbb.params->psym[i];
            i -= Gbb.params->poff[Gi];
            size_t a = Gbb.params->roworb[h][ia][1];
            int Ga = Gbb.params->qsym[a];
            a -= Gbb.params->qoff[Ga];
            for(size_t jb = 0; jb < Gbb.params->coltot[h]; ++jb){
                size_t j = Gbb.params->colorb[h][jb][0];
                int Gj = Gbb.params->rsym[j];
                j -= Gbb.params->roff[Gj];
                size_t b = Gbb.params->colorb[h][jb][1];
                int Gb = Gbb.params->ssym[b];
                b -= Gbb.params->soff[Gb];
                if(Gi == Gj && Ga == Gb) {
                    Gbb.matrix[h][ia][jb] += (bKappa(Gi, i, j) + bOccTau(Gi, i, j)) * bVirPTau(Ga, a, b);
                    Gbb.matrix[h][ia][jb] += bVirTau(Ga, a, b) * (bOccPTau(Gi, i, j) - bOccTau(Gi, i, j));
                }
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Gbb, h);
    }

    dpd_buf4_print(&Gbb,outfile,1);

    dpd_buf4_close(&Gbb);

//    /*
//     * As a sanity check, compute the energy from the density matrices.
//     */
//    Matrix aH(so_h_);
//    Matrix bH(so_h_);
//    aH.transform(Ca_);
//    bH.transform(Cb_);
//    double hCoreEnergy = 0.0;
//    for(int h = 0; h < nirrep_; ++h){
//        for(int p = 0; p < aOccOPDM.rowspi()[h]; ++p){
//            for(int q = 0; q < aOccOPDM.colspi()[h]; ++q){
//                hCoreEnergy += aH.get(h, p, q) * aOccOPDM.get(h, p, q);
//            }
//        }
//        for(int p = 0; p < bOccOPDM.rowspi()[h]; ++p){
//            for(int q = 0; q < bOccOPDM.colspi()[h]; ++q){
//                hCoreEnergy += bH.get(h, p, q) * bOccOPDM.get(h, p, q);
//            }
//        }
//        for(int p = 0; p < aVirOPDM.rowspi()[h]; ++p){
//            for(int q = 0; q < aVirOPDM.colspi()[h]; ++q){
//                hCoreEnergy += aH.get(h, p + naoccpi_[h], q + naoccpi_[h]) * aVirOPDM.get(h, p, q);
//            }
//        }
//        for(int p = 0; p < bVirOPDM.rowspi()[h]; ++p){
//            for(int q = 0; q < bVirOPDM.colspi()[h]; ++q){
//                hCoreEnergy += bH.get(h, p + nboccpi_[h], q + nboccpi_[h]) * bVirOPDM.get(h, p, q);
//            }
//        }
//    }

//    // OOVV
//    double OOVVEnergy = 0.0;
//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
//                      ID("[O,O]"), ID("[V,V]"), 0, "Gamma <OO|VV>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
//                      ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
//    OOVVEnergy += 0.25* dpd_buf4_dot(&I, &L);
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
//                      ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
//                      ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
//    OOVVEnergy += dpd_buf4_dot(&I, &L);
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
//                      ID("[o,o]"), ID("[v,v]"), 0, "Gamma <oo|vv>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
//                      ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
//    OOVVEnergy += 0.25 * dpd_buf4_dot(&I, &L);
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    // VVOO
//    double VVOOEnergy = 0.0;
//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[O,O]"),
//                      ID("[V,V]"), ID("[O,O]"), 0, "Gamma <VV|OO>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
//                      ID("[V,V]"), ID("[O,O]"), 1, "MO Ints <VV|OO>");
//    VVOOEnergy += 0.25* dpd_buf4_dot(&I, &L);
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"),
//                      ID("[V,v]"), ID("[O,o]"), 0, "Gamma <Vv|Oo>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[O,o]"),
//                      ID("[V,v]"), ID("[O,o]"), 0, "MO Ints <Vv|Oo>");
//    VVOOEnergy += dpd_buf4_dot(&I, &L);
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[o,o]"),
//                      ID("[v,v]"), ID("[o,o]"), 0, "Gamma <vv|oo>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,o]"),
//                      ID("[v,v]"), ID("[o,o]"), 1, "MO Ints <vv|oo>");
//    VVOOEnergy += 0.25 * dpd_buf4_dot(&I, &L);
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    // VVVV
//    double VVVVEnergy = 0.0;
//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
//                      ID("[V,V]"), ID("[V,V]"), 0, "Gamma <VV|VV>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
//                      ID("[V,V]"), ID("[V,V]"), 1, "MO Ints <VV|VV>");
//    VVVVEnergy += 0.25* dpd_buf4_dot(&I, &L);
////fprintf(outfile, "testaa = %16.10f\n", 0.25*dpd_buf4_dot(&I, &L));
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
//                      ID("[V,v]"), ID("[V,v]"), 0, "Gamma <Vv|Vv>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
//                      ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
//    VVVVEnergy += dpd_buf4_dot(&I, &L);
////fprintf(outfile, "testab = %16.10f\n", dpd_buf4_dot(&I, &L));
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
//                      ID("[v,v]"), ID("[v,v]"), 0, "Gamma <vv|vv>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
//                      ID("[v,v]"), ID("[v,v]"), 1, "MO Ints <vv|vv>");
//    VVVVEnergy += 0.25 * dpd_buf4_dot(&I, &L);
////fprintf(outfile, "testbb = %16.10f\n", 0.25*dpd_buf4_dot(&I, &L));
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    // OOOO
//    double OOOOEnergy = 0.0;
//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
//                      ID("[O,O]"), ID("[O,O]"), 0, "Gamma <OO|OO>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
//                      ID("[O,O]"), ID("[O,O]"), 1, "MO Ints <OO|OO>");
//    OOOOEnergy += 0.25* dpd_buf4_dot(&I, &L);
////fprintf(outfile, "testaa = %16.10f\n", 0.25*dpd_buf4_dot(&I, &L));
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
//                      ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
//                      ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
//    OOOOEnergy += dpd_buf4_dot(&I, &L);

////fprintf(outfile, "testab = %16.10f\n", dpd_buf4_dot(&I, &L));
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
//                      ID("[o,o]"), ID("[o,o]"), 0, "Gamma <oo|oo>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
//                      ID("[o,o]"), ID("[o,o]"), 1, "MO Ints <oo|oo>");
//    OOOOEnergy += 0.25 * dpd_buf4_dot(&I, &L);
////fprintf(outfile, "testbb = %16.10f\n", 0.25*dpd_buf4_dot(&I, &L));
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);


//    // OVOV: Note the different prefactor, due to 4-fold permutational symmetry
//    double OVOVEnergy = 0.0;
//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
//                      ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
//                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV> - (OV|OV)");
//    OVOVEnergy += dpd_buf4_dot(&I, &L);
////fprintf(outfile, "testaa = %16.10f\n", dpd_buf4_dot(&I, &L));
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
//                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
//                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>");
//    OVOVEnergy += dpd_buf4_dot(&L, &I);
////fprintf(outfile, "testab = %16.10f\n", dpd_buf4_dot(&I, &L));
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
//                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,V]"),
//                  ID("[o,V]"), ID("[o,V]"), 0, "MO Ints <oV|oV>");
//    OVOVEnergy += dpd_buf4_dot(&L, &I);
////fprintf(outfile, "testba = %16.10f\n", dpd_buf4_dot(&I, &L));
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
//                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[o,V]"),
//                  ID("[O,v]"), ID("[o,V]"), 0, "MO Ints <Ov|oV>");
//    OVOVEnergy += dpd_buf4_dot(&L, &I);
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
//                      ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");
//    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
//                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov> - (ov|ov)");
//    OVOVEnergy += dpd_buf4_dot(&I, &L);
////fprintf(outfile, "testbb = %16.10f\n", dpd_buf4_dot(&I, &L));
//    dpd_buf4_close(&I);
//    dpd_buf4_close(&L);

//    double twoElectronEnergy = + OOOOEnergy + VVVVEnergy + OOVVEnergy + VVOOEnergy + OVOVEnergy;

//    fprintf(outfile, "\tEnergy recomputed from the density matrices (should match the value shown above)...\n\n");
//    fprintf(outfile, "\tNuclear Repulsion energy  = %16.10f\n", enuc_);
//    fprintf(outfile, "\tOne electron energy       = %16.10f\n", hCoreEnergy);
//    fprintf(outfile, "\tTwo electron energy       = %16.10f\n", twoElectronEnergy);
//    fprintf(outfile, "\t    OOOO Energy   = %16.10f\n", OOOOEnergy);
//    fprintf(outfile, "\t    OVOV Energy   = %16.10f\n", OVOVEnergy);
//    fprintf(outfile, "\t    OOVV Energy   = %16.10f\n", OOVVEnergy);
//    fprintf(outfile, "\t    VVOO Energy   = %16.10f\n", VVOOEnergy);
//    fprintf(outfile, "\t    VVVV Energy   = %16.10f\n", VVVVEnergy);
//    fprintf(outfile, "\t--------------------------------------------\n");
//    fprintf(outfile, "\tTotal Energy              = %16.10f\n", enuc_ + hCoreEnergy + twoElectronEnergy);

    psio_->close(PSIF_DCFT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);

}


}} //End namespaces
