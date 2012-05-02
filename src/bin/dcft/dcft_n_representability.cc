#include "dcft.h"
#include "defines.h"
#include <libpsio/psio.hpp>
#include <psifiles.h>
#include <libtrans/integraltransform.h>

namespace psi{ namespace dcft{

/**
 * Dumps the MO Basis density matrix in a DPD form.
 */
void
DCFTSolver::dump_density()
{
    dcft_timer_on("DCFTSolver::dump_density()");

    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 Laa, Lab, Lbb, Gaa, Gab, Gba, Gbb, I, L;
    dpdfile2 T_OO, T_oo, T_VV, T_vv;

    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
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

    Matrix aOccOPDM(nirrep_, naoccpi_, naoccpi_);
    Matrix bOccOPDM(nirrep_, nboccpi_, nboccpi_);
    Matrix aVirOPDM(nirrep_, navirpi_, navirpi_);
    Matrix bVirOPDM(nirrep_, nbvirpi_, nbvirpi_);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int j = 0; j < naoccpi_[h]; ++j){
                aOccOPDM.set(h, i, j, (i == j ? 1.0 : 0.0) + T_OO.matrix[h][i][j]);
            }
        }
        for(int a = 0; a < navirpi_[h]; ++a){
            for(int b = 0; b < navirpi_[h]; ++b){
                aVirOPDM.set(h, a, b, T_VV.matrix[h][a][b]);
            }
        }
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int j = 0; j < nboccpi_[h]; ++j){
                bOccOPDM.set(h, i, j, (i == j ? 1.0 : 0.0) + T_oo.matrix[h][i][j]);
            }
        }
        for(int a = 0; a < nbvirpi_[h]; ++a){
            for(int b = 0; b < nbvirpi_[h]; ++b){
                bVirOPDM.set(h, a, b, T_vv.matrix[h][a][b]);
            }
        }
    }

    dpd_file2_close(&T_OO);
    dpd_file2_close(&T_oo);
    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);

    /*
     * The VVVV block
     */
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
              ID("[V,V]"), ID("[V,V]"), 0, "Gamma <VV|VV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gaa, h);
        dpd_buf4_mat_irrep_init(&Laa, h);
        dpd_buf4_mat_irrep_rd(&Laa, h);
        #pragma omp parallel for
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
                    tpdm += 0.5 * Laa.matrix[h][ij][ab] * Laa.matrix[h][ij][cd];
                }
                size_t c = Gaa.params->colorb[h][cd][0];
                int Gc = Gaa.params->rsym[c];
                c -= Gaa.params->roff[Gc];
                size_t d = Gaa.params->colorb[h][cd][1];
                int Gd = Gaa.params->ssym[d];
                d -= Gaa.params->soff[Gd];
                if(Ga == Gc && Gb == Gd) tpdm += aVirOPDM(Ga, a, c) * aVirOPDM(Gb, b, d);
                if(Ga == Gd && Gb == Gc) tpdm -= aVirOPDM(Ga, a, d) * aVirOPDM(Gb, b, c);
                Gaa.matrix[h][ab][cd] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Laa, h);
    }
    dpd_buf4_close(&Gaa);

    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
              ID("[V,v]"), ID("[V,v]"), 0, "Gamma <Vv|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gab, h);
        dpd_buf4_mat_irrep_init(&Lab, h);
        dpd_buf4_mat_irrep_rd(&Lab, h);
        #pragma omp parallel for
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
                    tpdm += Lab.matrix[h][ij][ab] * Lab.matrix[h][ij][cd];
                }
                size_t c = Gab.params->colorb[h][cd][0];
                int Gc = Gab.params->rsym[c];
                c -= Gab.params->roff[Gc];
                size_t d = Gab.params->colorb[h][cd][1];
                int Gd = Gab.params->ssym[d];
                d -= Gab.params->soff[Gd];
                if(Ga == Gc && Gb == Gd) tpdm += aVirOPDM(Ga, a, c) * bVirOPDM(Gb, b, d);
                Gab.matrix[h][ab][cd] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gab, h);
        dpd_buf4_mat_irrep_close(&Gab, h);
        dpd_buf4_mat_irrep_close(&Lab, h);
    }
    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
              ID("[v,v]"), ID("[v,v]"), 0, "Gamma <vv|vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gbb, h);
        dpd_buf4_mat_irrep_init(&Lbb, h);
        dpd_buf4_mat_irrep_rd(&Lbb, h);
        #pragma omp parallel for
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
                    tpdm += 0.5 * Lbb.matrix[h][ij][ab] * Lbb.matrix[h][ij][cd];
                }
                size_t c = Gbb.params->colorb[h][cd][0];
                int Gc = Gbb.params->rsym[c];
                c -= Gbb.params->roff[Gc];
                size_t d = Gbb.params->colorb[h][cd][1];
                int Gd = Gbb.params->ssym[d];
                d -= Gbb.params->soff[Gd];
                if(Ga == Gc && Gb == Gd) tpdm += bVirOPDM(Ga, a, c) * bVirOPDM(Gb, b, d);
                if(Ga == Gd && Gb == Gc) tpdm -= bVirOPDM(Ga, a, d) * bVirOPDM(Gb, b, c);
                Gbb.matrix[h][ab][cd] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Lbb, h);
    }
    dpd_buf4_close(&Gbb);

    /*
     * The OOOO  block
     */
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
              ID("[O,O]"), ID("[O,O]"), 0, "Gamma <OO|OO>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gaa, h);
        dpd_buf4_mat_irrep_init(&Laa, h);
        dpd_buf4_mat_irrep_rd(&Laa, h);
        #pragma omp parallel for
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
                    tpdm += 0.5 * Laa.matrix[h][ij][ab] * Laa.matrix[h][kl][ab];
                }
                size_t k = Gaa.params->colorb[h][kl][0];
                int Gk = Gaa.params->rsym[k];
                k -= Gaa.params->roff[Gk];
                size_t l = Gaa.params->colorb[h][kl][1];
                int Gl = Gaa.params->ssym[l];
                l -= Gaa.params->soff[Gl];
                if(Gi == Gk && Gj == Gl) tpdm += aOccOPDM(Gi, i, k) * aOccOPDM(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= aOccOPDM(Gi, i, l) * aOccOPDM(Gj, j, k);
                Gaa.matrix[h][ij][kl] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Laa, h);
    }
    dpd_buf4_close(&Gaa);


    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
              ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gab, h);
        dpd_buf4_mat_irrep_init(&Lab, h);
        dpd_buf4_mat_irrep_rd(&Lab, h);
        #pragma omp parallel for
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
                    tpdm += Lab.matrix[h][ij][ab] * Lab.matrix[h][kl][ab];
                }
                size_t k = Gab.params->colorb[h][kl][0];
                int Gk = Gab.params->rsym[k];
                k -= Gab.params->roff[Gk];
                size_t l = Gab.params->colorb[h][kl][1];
                int Gl = Gab.params->ssym[l];
                l -= Gab.params->soff[Gl];
                if(Gi == Gk && Gj == Gl) tpdm += aOccOPDM(Gi, i, k) * bOccOPDM(Gj, j, l);
                Gab.matrix[h][ij][kl] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gab, h);
        dpd_buf4_mat_irrep_close(&Gab, h);
        dpd_buf4_mat_irrep_close(&Lab, h);
    }
    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
              ID("[o,o]"), ID("[o,o]"), 0, "Gamma <oo|oo>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gbb, h);
        dpd_buf4_mat_irrep_init(&Lbb, h);
        dpd_buf4_mat_irrep_rd(&Lbb, h);
        #pragma omp parallel for
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
                    tpdm += 0.5 * Lbb.matrix[h][ij][ab] * Lbb.matrix[h][kl][ab];
                }
                size_t k = Gbb.params->colorb[h][kl][0];
                int Gk = Gbb.params->rsym[k];
                k -= Gbb.params->roff[Gk];
                size_t l = Gbb.params->colorb[h][kl][1];
                int Gl = Gbb.params->ssym[l];
                l -= Gbb.params->soff[Gl];
                if(Gi == Gk && Gj == Gl) tpdm += bOccOPDM(Gi, i, k) * bOccOPDM(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= bOccOPDM(Gi, i, l) * bOccOPDM(Gj, j, k);
                Gbb.matrix[h][ij][kl] = tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Lbb, h);
    }
    dpd_buf4_close(&Gbb);

    /*
     * The OOVV and VVOO blocks
     */
    dpd_buf4_init(&Gaa, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
              ID("[O,O]"), ID("[V,V]"), 0, "Lambda <OO|VV>");
    dpd_buf4_copy(&Gaa, PSIF_DCFT_DENSITY, "Gamma <OO|VV>");
    dpd_buf4_sort(&Gaa, PSIF_DCFT_DENSITY, rspq, ID("[V,V]"), ID("[O,O]"), "Gamma <VV|OO>");
    dpd_buf4_close(&Gaa);

    dpd_buf4_init(&Gab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
              ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_copy(&Gab, PSIF_DCFT_DENSITY, "Gamma <Oo|Vv>");
    dpd_buf4_sort(&Gab, PSIF_DCFT_DENSITY, rspq, ID("[V,v]"), ID("[O,o]"), "Gamma <Vv|Oo>");
    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
              ID("[o,o]"), ID("[v,v]"), 0, "Lambda <oo|vv>");
    dpd_buf4_copy(&Gbb, PSIF_DCFT_DENSITY, "Gamma <oo|vv>");
    dpd_buf4_sort(&Gbb, PSIF_DCFT_DENSITY, rspq, ID("[v,v]"), ID("[o,o]"), "Gamma <vv|oo>");
    dpd_buf4_close(&Gbb);

    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    dpd_buf4_sort(&I, PSIF_LIBTRANS_DPD, rspq, ID("[V,v]"), ID("[O,o]"), "MO Ints <Vv|Oo>");
    dpd_buf4_close(&I);

    dpd_buf4_close(&Laa);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Lbb);

    /*
     * The OVOV block
     */
    dpdbuf4 Laaaa, Laabb, Labba, Lbaab, Lbbbb, Taa, Tab, Tba, Tbb, Gaa2, Laa2, Lab2, Lbb2;

    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma (OV|OV)");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    dpd_buf4_init(&Laa2, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    dpd_contract444(&Laa, &Laa2, &Gaa, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Laa);
    dpd_buf4_close(&Laa2);
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_buf4_init(&Lab2, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_contract444(&Lab, &Lab2, &Gaa, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Lab2);
    dpd_buf4_close(&Gaa);

    // Resort Г(OV|OV) to the Г<OV|OV>
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma (OV|OV)");
    dpd_buf4_sort(&Gaa, PSIF_DCFT_DENSITY, psrq, ID("[O,V]"),ID("[O,V]"), "Gamma <OV|OV>");
    dpd_buf4_close(&Gaa);

    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gaa, h);
        dpd_buf4_mat_irrep_rd(&Gaa, h);
        #pragma omp parallel for
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
                if(Gi == Gj && Ga == Gb)
                    Gaa.matrix[h][ia][jb] += aOccOPDM(Gi, i, j) * aVirOPDM(Ga, a, b);
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Gaa, h);
    }
    dpd_buf4_close(&Gaa);

    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Lambda (Ov|oV)");
    dpd_buf4_init(&Lab2, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Lambda (Ov|oV)");
    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
    dpd_contract444(&Lab, &Lab2, &Gab, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Gab);
    dpd_buf4_init(&Gba, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");
    dpd_contract444(&Lab, &Lab2, &Gba, 1, 1, -1.0, 0.0);
    dpd_buf4_close(&Gba);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Lab2);

    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gab, h);
        dpd_buf4_mat_irrep_rd(&Gab, h);
        #pragma omp parallel for
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
                if(Gi == Gj && Ga == Gb)
                    Gab.matrix[h][ia][jb] += aOccOPDM(Gi, i, j) * bVirOPDM(Ga, a, b);
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gab, h);
        dpd_buf4_mat_irrep_close(&Gab, h);
    }
    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Gba, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gba, h);
        dpd_buf4_mat_irrep_rd(&Gba, h);
        #pragma omp parallel for
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
                if(Gi == Gj && Ga == Gb)
                    Gba.matrix[h][ia][jb] += bOccOPDM(Gi, i, j) * aVirOPDM(Ga, a, b);
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gba, h);
        dpd_buf4_mat_irrep_close(&Gba, h);
    }
    dpd_buf4_close(&Gba);

    // <Ov|oV> and <oV|Ov> (these are the same - compute 1 and double)
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");
    dpd_contract444(&Laa, &Lab, &Tab, 0, 1, 2.0, 0.0);
    dpd_buf4_close(&Laa);
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    dpd_contract444(&Lab, &Lbb, &Tab, 0, 1, 2.0, 1.0);
    dpd_buf4_close(&Lbb);
    dpd_buf4_close(&Tab);
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");
    dpd_buf4_sort(&Tab, PSIF_DCFT_DENSITY, psrq, ID("[O,v]"), ID("[o,V]"), "Gamma <Ov|oV>");
    dpd_buf4_close(&Tab);


    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma (ov|ov)");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    dpd_buf4_init(&Lbb2, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    dpd_contract444(&Lbb, &Lbb2, &Gbb, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Lbb);
    dpd_buf4_close(&Lbb2);
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_buf4_init(&Lab2, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_contract444(&Lab, &Lab2, &Gbb, 1, 1, -1.0, 1.0);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&Lab2);
    dpd_buf4_close(&Gbb);

    // Resort Г(ov|ov) to the Г<ov|ov>
    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma (ov|ov)");
    dpd_buf4_sort(&Gbb, PSIF_DCFT_DENSITY, psrq, ID("[o,v]"), ID("[o,v]"), "Gamma <ov|ov>");
    dpd_buf4_close(&Gbb);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gbb, h);
        dpd_buf4_mat_irrep_rd(&Gbb, h);
        #pragma omp parallel for
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
                if(Gi == Gj && Ga == Gb)
                    Gbb.matrix[h][ia][jb] += bOccOPDM(Gi, i, j) * bVirOPDM(Ga, a, b);
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Gbb, h);
    }
    dpd_buf4_close(&Gbb);

    /*
     * As a sanity check, compute the energy from the density matrices.
     */
    Matrix aH(so_h_);
    Matrix bH(so_h_);
    aH.transform(Ca_);
    bH.transform(Cb_);
    double hCoreEnergy = 0.0;
    for(int h = 0; h < nirrep_; ++h){
        for(int p = 0; p < aOccOPDM.rowspi()[h]; ++p){
            for(int q = 0; q < aOccOPDM.colspi()[h]; ++q){
                hCoreEnergy += aH.get(h, p, q) * aOccOPDM.get(h, p, q);
            }
        }
        for(int p = 0; p < bOccOPDM.rowspi()[h]; ++p){
            for(int q = 0; q < bOccOPDM.colspi()[h]; ++q){
                hCoreEnergy += bH.get(h, p, q) * bOccOPDM.get(h, p, q);
            }
        }
        for(int p = 0; p < aVirOPDM.rowspi()[h]; ++p){
            for(int q = 0; q < aVirOPDM.colspi()[h]; ++q){
                hCoreEnergy += aH.get(h, p + naoccpi_[h], q + naoccpi_[h]) * aVirOPDM.get(h, p, q);
            }
        }
        for(int p = 0; p < bVirOPDM.rowspi()[h]; ++p){
            for(int q = 0; q < bVirOPDM.colspi()[h]; ++q){
                hCoreEnergy += bH.get(h, p + nboccpi_[h], q + nboccpi_[h]) * bVirOPDM.get(h, p, q);
            }
        }
    }

    // OOVV
    double OOVVEnergy = 0.0;
    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 0, "Gamma <OO|VV>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                      ID("[O,O]"), ID("[V,V]"), 1, "MO Ints <OO|VV>");
    OOVVEnergy += 0.25* dpd_buf4_dot(&I, &L);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                      ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    OOVVEnergy += dpd_buf4_dot(&I, &L);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 0, "Gamma <oo|vv>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                      ID("[o,o]"), ID("[v,v]"), 1, "MO Ints <oo|vv>");
    OOVVEnergy += 0.25 * dpd_buf4_dot(&I, &L);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    // VVOO
    double VVOOEnergy = 0.0;
    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[O,O]"),
                      ID("[V,V]"), ID("[O,O]"), 0, "Gamma <VV|OO>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                      ID("[V,V]"), ID("[O,O]"), 1, "MO Ints <VV|OO>");
    VVOOEnergy += 0.25* dpd_buf4_dot(&I, &L);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"),
                      ID("[V,v]"), ID("[O,o]"), 0, "Gamma <Vv|Oo>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[O,o]"),
                      ID("[V,v]"), ID("[O,o]"), 0, "MO Ints <Vv|Oo>");
    VVOOEnergy += dpd_buf4_dot(&I, &L);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[o,o]"),
                      ID("[v,v]"), ID("[o,o]"), 0, "Gamma <vv|oo>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,o]"),
                      ID("[v,v]"), ID("[o,o]"), 1, "MO Ints <vv|oo>");
    VVOOEnergy += 0.25 * dpd_buf4_dot(&I, &L);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    // VVVV
    double VVVVEnergy = 0.0;
    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                      ID("[V,V]"), ID("[V,V]"), 0, "Gamma <VV|VV>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                      ID("[V,V]"), ID("[V,V]"), 1, "MO Ints <VV|VV>");
    VVVVEnergy += 0.25* dpd_buf4_dot(&I, &L);
//fprintf(outfile, "testaa = %16.10f\n", 0.25*dpd_buf4_dot(&I, &L));
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                      ID("[V,v]"), ID("[V,v]"), 0, "Gamma <Vv|Vv>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                      ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
    VVVVEnergy += dpd_buf4_dot(&I, &L);
//fprintf(outfile, "testab = %16.10f\n", dpd_buf4_dot(&I, &L));
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                      ID("[v,v]"), ID("[v,v]"), 0, "Gamma <vv|vv>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                      ID("[v,v]"), ID("[v,v]"), 1, "MO Ints <vv|vv>");
    VVVVEnergy += 0.25 * dpd_buf4_dot(&I, &L);
//fprintf(outfile, "testbb = %16.10f\n", 0.25*dpd_buf4_dot(&I, &L));
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    // OOOO
    double OOOOEnergy = 0.0;
    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                      ID("[O,O]"), ID("[O,O]"), 0, "Gamma <OO|OO>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                      ID("[O,O]"), ID("[O,O]"), 1, "MO Ints <OO|OO>");
    OOOOEnergy += 0.25* dpd_buf4_dot(&I, &L);
//fprintf(outfile, "testaa = %16.10f\n", 0.25*dpd_buf4_dot(&I, &L));
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                      ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                      ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
    OOOOEnergy += dpd_buf4_dot(&I, &L);

//fprintf(outfile, "testab = %16.10f\n", dpd_buf4_dot(&I, &L));
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                      ID("[o,o]"), ID("[o,o]"), 0, "Gamma <oo|oo>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                      ID("[o,o]"), ID("[o,o]"), 1, "MO Ints <oo|oo>");
    OOOOEnergy += 0.25 * dpd_buf4_dot(&I, &L);
//fprintf(outfile, "testbb = %16.10f\n", 0.25*dpd_buf4_dot(&I, &L));
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);


    // OVOV: Note the different prefactor, due to 4-fold permutational symmetry
    double OVOVEnergy = 0.0;
    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                      ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV> - <OV|VO>");
    OVOVEnergy += dpd_buf4_dot(&I, &L);
//fprintf(outfile, "testaa = %16.10f\n", dpd_buf4_dot(&I, &L));
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>");
    OVOVEnergy += dpd_buf4_dot(&L, &I);
//fprintf(outfile, "testab = %16.10f\n", dpd_buf4_dot(&I, &L));
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "MO Ints <oV|oV>");
    OVOVEnergy += dpd_buf4_dot(&L, &I);
//fprintf(outfile, "testba = %16.10f\n", dpd_buf4_dot(&I, &L));
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "MO Ints <Ov|oV>");
    OVOVEnergy += dpd_buf4_dot(&L, &I);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    dpd_buf4_init(&L, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                      ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov> - <ov|vo>");
    OVOVEnergy += dpd_buf4_dot(&I, &L);
//fprintf(outfile, "testbb = %16.10f\n", dpd_buf4_dot(&I, &L));
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);

    double twoElectronEnergy = + OOOOEnergy + VVVVEnergy + OOVVEnergy + VVOOEnergy + OVOVEnergy;

    fprintf(outfile, "\tEnergy recomputed from the density matrices (should match the value shown above)...\n\n");
    fprintf(outfile, "\tNuclear Repulsion energy  = %16.10f\n", enuc_);
    fprintf(outfile, "\tOne electron energy       = %16.10f\n", hCoreEnergy);
    fprintf(outfile, "\tTwo electron energy       = %16.10f\n", twoElectronEnergy);
    fprintf(outfile, "\t    OOOO Energy   = %16.10f\n", OOOOEnergy);
    fprintf(outfile, "\t    OVOV Energy   = %16.10f\n", OVOVEnergy);
    fprintf(outfile, "\t    OOVV Energy   = %16.10f\n", OOVVEnergy);
    fprintf(outfile, "\t    VVOO Energy   = %16.10f\n", VVOOEnergy);
    fprintf(outfile, "\t    VVVV Energy   = %16.10f\n", VVVVEnergy);
    fprintf(outfile, "\t--------------------------------------------\n");
    fprintf(outfile, "\tTotal Energy              = %16.10f\n", enuc_ + hCoreEnergy + twoElectronEnergy);

    psio_->close(PSIF_DCFT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);

    dcft_timer_off("DCFTSolver::dump_density()");
}


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

    for(int h = 0; h < nirrep_; ++h){
        int nOcc = naoccpi_[h];
        int nVir = navirpi_[h];
        unsigned long int dim = nmopi_[h];

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
        #pragma omp parallel for
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
        #pragma omp parallel for
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
        #pragma omp parallel for
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

/**
 * Prints the mulliken charges.
 */
void
DCFTSolver::mulliken_charges()
{
#if 0
    SimpleMatrix aOPDM(_aKappa.to_simple_matrix());
    SimpleMatrix bOPDM(_bKappa.to_simple_matrix());
    int offset = 0;
    for(int h = 0; h < _nIrreps; ++h){
        for(int row = 0; row < _soPI[h]; ++row){
            for(int col = 0; col < _soPI[h]; ++col){
                aOPDM.add(row + offset, col + offset, _aTau[h][row][col]);
                bOPDM.add(row + offset, col + offset, _aTau[h][row][col]);
            }
        }
        offset += _soPI[h];
    }

    SimpleMatrix aDS(_nSo, _nSo);
    SimpleMatrix bDS(_nSo, _nSo);
    SimpleMatrix S(_aoS.to_simple_matrix());
    // Form the DS matrix (well, Kappa X S really)
    aDS.gemm(0, 0, 1.0, &aOPDM, &S, 0.0);
    bDS.gemm(0, 0, 1.0, &aOPDM, &S, 0.0);
    int nAtoms         = _chkpt->rd_natom();
    int nShells        = _chkpt->rd_nshell();
    int *shellLocation = _chkpt->rd_sloc_new();
    int *shellAngMom   = _chkpt->rd_stype();
    int *shellToAtom   = _chkpt->rd_snuc();
    double *zVals      = _chkpt->rd_zvals();
    double **U         = _chkpt->rd_usotbf();
    char **labels      = _chkpt->rd_felement();
    SimpleMatrix aoADS(_nSo, _nSo);
    SimpleMatrix aoBDS(_nSo, _nSo);
    SimpleMatrix UMat(_nSo, _nSo);
    UMat.set(U);
    SimpleMatrix tmp(_nSo, _nSo);
    // Transform the SO DS matrix into the AO basis
    tmp.gemm(1, 0, 1.0, &UMat, &aDS, 0.0);
    aoADS.gemm(0, 0, 1.0, &tmp, &UMat, 0.0);
    tmp.gemm(1, 0, 1.0, &UMat, &bDS, 0.0);
    aoBDS.gemm(0, 0, 1.0, &tmp, &UMat, 0.0);
    SimpleVector aCharges(nAtoms);
    SimpleVector bCharges(nAtoms);
    for(int shell = 0; shell < nShells; ++shell){
        int atom       = shellToAtom[shell] - 1;
        int firstBF    = shellLocation[shell] - 1;
        int shellSize  = 2*shellAngMom[shell] - 1;
        for(int bf = 0; bf < shellSize; ++bf){
            aCharges[atom] += aoADS.get(bf + firstBF, bf + firstBF);
            bCharges[atom] += aoBDS.get(bf + firstBF, bf + firstBF);
        }
    }

    fprintf(outfile, "\n\t\t     Mulliken Population Analysis\n\n");
    fprintf(outfile,   "\t\t    Atom          Charge      Spin\n");
    fprintf(outfile,   "\t\t-------------------------------------\n");
    for(int atom = 0; atom < nAtoms; ++atom){
        fprintf(outfile, "\t\t  %-12s  %8.3f   %8.3f\n",
                labels[atom], zVals[atom] - aCharges[atom] - bCharges[atom],
                aCharges[atom] - bCharges[atom]);
    }
#endif
}

}} // Namespaces
