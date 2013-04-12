#include <libtrans/integraltransform.h>
#include <libpsio/psio.hpp>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <libdiis/diismanager.h>
#include "dcft.h"
#include "defines.h"

namespace psi{ namespace dcft{

void
DCFTSolver::oo_init() {

    aocc_tau_ = SharedMatrix(new Matrix("MO basis Tau (Alpha Occupied)", nirrep_, naoccpi_, naoccpi_));
    bocc_tau_ = SharedMatrix(new Matrix("MO basis Tau (Beta Occupied)", nirrep_, nboccpi_, nboccpi_));
    avir_tau_ = SharedMatrix(new Matrix("MO basis Tau (Alpha Virtual)", nirrep_, navirpi_, navirpi_));
    bvir_tau_ = SharedMatrix(new Matrix("MO basis Tau (Beta Virtual)", nirrep_, nbvirpi_, nbvirpi_));
    akappa_ = SharedMatrix(new Matrix("MO basis Kappa (Alpha)", nirrep_, naoccpi_, naoccpi_));
    bkappa_ = SharedMatrix(new Matrix("MO basis Kappa (Beta)", nirrep_, nboccpi_, nboccpi_));
    orbital_gradient_a_ = SharedMatrix(new Matrix("MO basis Orbital Gradient (Alpha)", nirrep_, nmopi_, nmopi_));
    orbital_gradient_b_ = SharedMatrix(new Matrix("MO basis Orbital Gradient (Beta)", nirrep_, nmopi_, nmopi_));
    X_a_ = SharedMatrix(new Matrix("Generator of the orbital rotations w.r.t. previous orbitals (Alpha)", nirrep_, nmopi_, nmopi_));
    X_b_ = SharedMatrix(new Matrix("Generator of the orbital rotations w.r.t. previous orbitals (Beta)", nirrep_, nmopi_, nmopi_));
    Xtotal_a_ = SharedMatrix(new Matrix("Generator of the orbital rotations w.r.t. reference orbitals (Alpha)", nirrep_, nmopi_, nmopi_));
    Xtotal_b_ = SharedMatrix(new Matrix("Generator of the orbital rotations w.r.t. reference orbitals (Beta)", nirrep_, nmopi_, nmopi_));

}

double
DCFTSolver::compute_orbital_residual() {

    dpdfile2 Xai, Xia;

    // Compute the unrelaxed densities for the orbital gradient
    compute_unrelaxed_density();

    // Compute the OV part of the orbital gradient
    compute_orbital_gradient_OV();

    // Compute the VO part of the orbital gradient
    compute_orbital_gradient_VO();

    dpd_file2_init(&Xia, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_file2_init(&Xai, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_file2_mat_init(&Xia);
    dpd_file2_mat_init(&Xai);
    dpd_file2_mat_rd(&Xia);
    dpd_file2_mat_rd(&Xai);

    double maxGradient = 0.0;
    // Alpha spin
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = 0; a < navirpi_[h]; ++a){
                double value = 2.0 * (Xia.matrix[h][i][a] - Xai.matrix[h][a][i]);
                maxGradient = (fabs(value) > maxGradient) ? fabs(value) : maxGradient;
                orbital_gradient_a_->set(h, i, a + naoccpi_[h], value);
                orbital_gradient_a_->set(h, a + naoccpi_[h], i, (-1.0) * value);
            }
        }
    }

    dpd_file2_close(&Xai);
    dpd_file2_close(&Xia);

    dpd_file2_init(&Xia, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_file2_init(&Xai, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_file2_mat_init(&Xia);
    dpd_file2_mat_init(&Xai);
    dpd_file2_mat_rd(&Xia);
    dpd_file2_mat_rd(&Xai);

    // Beta spin
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int a = 0; a < nbvirpi_[h]; ++a){
                double value = 2.0 * (Xia.matrix[h][i][a] - Xai.matrix[h][a][i]);
                maxGradient = (fabs(value) > maxGradient) ? fabs(value) : maxGradient;
                orbital_gradient_b_->set(h, i, a + nboccpi_[h], value);
                orbital_gradient_b_->set(h, a + nboccpi_[h], i, (-1.0) * value);
            }
        }
    }

    dpd_file2_close(&Xai);
    dpd_file2_close(&Xia);

    // Determine the orbital rotation step
    // Alpha spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int a = naoccpi_[h]; a < nmopi_[h]; ++a){
                double value = orbital_gradient_a_->get(h, i, a) / (2.0 * (moFa_->get(h, i, i) - moFa_->get(h, a, a)));
                X_a_->set(h, i, a, value);
                X_a_->set(h, a, i, (-1.0) * value);
            }
        }
    }

    // Beta spin
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int a = nboccpi_[h]; a < nmopi_[h]; ++a){
                double value = orbital_gradient_b_->get(h, i, a) / (2.0 * (moFb_->get(h, i, i) - moFb_->get(h, a, a)));
                X_b_->set(h, i, a, value);
                X_b_->set(h, a, i, (-1.0) * value);
            }
        }
    }

    // Determine the rotation generator with respect to the reference orbitals
    Xtotal_a_->add(X_a_);
    Xtotal_b_->add(X_b_);

    return maxGradient;
}

void
DCFTSolver::compute_unrelaxed_density() {

    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);

    dpdbuf4 LLaa, LLab, LLbb, Laa, Lab, Lbb, Gaa, Gab, Gba, Gbb, Tab;
    dpdfile2 T_OO, T_oo, T_VV, T_vv;

    // Compute the N^6 terms for Gamma OOOO

    // Gamma_ijkl = 1/16 (Lambda_ijab * Z_klab + Z_ijab * Lambda_klab)
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O>O]-"), ID("[O>O]-"),
              ID("[O>O]-"), ID("[O>O]-"), 0, "Gamma <OO|OO>");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&LLaa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_contract444(&Laa, &LLaa, &Gaa, 0, 0, 0.25, 0.0);
    dpd_buf4_close(&LLaa);
    dpd_buf4_close(&Gaa);
    dpd_buf4_close(&Laa);

    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&LLab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
              ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");
    dpd_contract444(&Lab, &LLab, &Gab, 0, 0, 0.25, 0.0);
    dpd_buf4_close(&Gab);
    dpd_buf4_close(&LLab);
    dpd_buf4_close(&Lab);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o>o]-"), ID("[o>o]-"),
              ID("[o>o]-"), ID("[o>o]-"), 0, "Gamma <oo|oo>");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&LLbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_contract444(&Lbb, &LLbb, &Gbb, 0, 0, 0.25, 0.0);
    dpd_buf4_close(&LLbb);
    dpd_buf4_close(&Gbb);
    dpd_buf4_close(&Lbb);

    // Add the terms containing one-particle densities to Gamma VVVV and OOOO

    dpd_file2_init(&T_OO, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    dpd_file2_init(&T_oo, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    dpd_file2_mat_init(&T_OO);
    dpd_file2_mat_init(&T_oo);
    dpd_file2_mat_init(&T_VV);
    dpd_file2_mat_init(&T_vv);

    dpd_file2_mat_rd(&T_OO);
    dpd_file2_mat_rd(&T_oo);
    dpd_file2_mat_rd(&T_VV);
    dpd_file2_mat_rd(&T_vv);

    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < naoccpi_[h]; ++i){
            for(int j = 0; j < naoccpi_[h]; ++j){
                akappa_->set(h, i, j, (i == j ? 1.0 : 0.0));
                aocc_tau_->set(h, i, j, T_OO.matrix[h][i][j]);
            }
        }
        for(int a = 0; a < navirpi_[h]; ++a){
            for(int b = 0; b < navirpi_[h]; ++b){
                avir_tau_->set(h, a, b, T_VV.matrix[h][a][b]);
            }
        }
        for(int i = 0; i < nboccpi_[h]; ++i){
            for(int j = 0; j < nboccpi_[h]; ++j){
                bkappa_->set(h, i, j, (i == j ? 1.0 : 0.0));
                bocc_tau_->set(h, i, j, T_oo.matrix[h][i][j]);
            }
        }
        for(int a = 0; a < nbvirpi_[h]; ++a){
            for(int b = 0; b < nbvirpi_[h]; ++b){
                bvir_tau_->set(h, a, b, T_vv.matrix[h][a][b]);
            }
        }
    }

    dpd_file2_close(&T_OO);
    dpd_file2_close(&T_oo);
    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);

    /*
     * The OOOO  block
     */
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
              ID("[O>O]-"), ID("[O>O]-"), 0, "Gamma <OO|OO>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gaa, h);
        dpd_buf4_mat_irrep_rd(&Gaa, h);

        #pragma omp parallel for
        for(long int ij = 0; ij < Gaa.params->rowtot[h]; ++ij){
            size_t i = Gaa.params->roworb[h][ij][0];
            int Gi = Gaa.params->psym[i];
            i -= Gaa.params->poff[Gi];
            size_t j = Gaa.params->roworb[h][ij][1];
            int Gj = Gaa.params->qsym[j];
            j -= Gaa.params->qoff[Gj];
            for(size_t kl = 0; kl < Gaa.params->coltot[h]; ++kl){
                double tpdm = 0.0;
                size_t k = Gaa.params->colorb[h][kl][0];
                int Gk = Gaa.params->rsym[k];
                k -= Gaa.params->roff[Gk];
                size_t l = Gaa.params->colorb[h][kl][1];
                int Gl = Gaa.params->ssym[l];
                l -= Gaa.params->soff[Gl];

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * akappa_->get(Gi, i, k) * akappa_->get(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * akappa_->get(Gi, i, l) * akappa_->get(Gj, j, k);

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * akappa_->get(Gi, i, k) * aocc_tau_->get(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * akappa_->get(Gi, i, l) * aocc_tau_->get(Gj, j, k);
                if(Gj == Gk && Gi == Gl) tpdm -= 0.25 * akappa_->get(Gj, j, k) * aocc_tau_->get(Gi, i, l);
                if(Gj == Gl && Gi == Gk) tpdm += 0.25 * akappa_->get(Gj, j, l) * aocc_tau_->get(Gi, i, k);

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * aocc_tau_->get(Gi, i, k) * aocc_tau_->get(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * aocc_tau_->get(Gi, i, l) * aocc_tau_->get(Gj, j, k);

                Gaa.matrix[h][ij][kl] += tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Gaa, h);
    }

    dpd_buf4_close(&Gaa);


    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
              ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gab, h);
        dpd_buf4_mat_irrep_rd(&Gab, h);

        #pragma omp parallel for
        for(long int ij = 0; ij < Gab.params->rowtot[h]; ++ij){
            size_t i = Gab.params->roworb[h][ij][0];
            int Gi = Gab.params->psym[i];
            i -= Gab.params->poff[Gi];
            size_t j = Gab.params->roworb[h][ij][1];
            int Gj = Gab.params->qsym[j];
            j -= Gab.params->qoff[Gj];
            for(size_t kl = 0; kl < Gab.params->coltot[h]; ++kl){
                double tpdm = 0.0;
                size_t k = Gab.params->colorb[h][kl][0];
                int Gk = Gab.params->rsym[k];
                k -= Gab.params->roff[Gk];
                size_t l = Gab.params->colorb[h][kl][1];
                int Gl = Gab.params->ssym[l];
                l -= Gab.params->soff[Gl];
                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * akappa_->get(Gi, i, k) * bkappa_->get(Gj, j, l);

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * akappa_->get(Gi, i, k) * bocc_tau_->get(Gj, j, l);
                if(Gj == Gl && Gi == Gk) tpdm += 0.25 * bkappa_->get(Gj, j, l) * aocc_tau_->get(Gi, i, k);

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * aocc_tau_->get(Gi, i, k) * bocc_tau_->get(Gj, j, l);

                Gab.matrix[h][ij][kl] += tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gab, h);
        dpd_buf4_mat_irrep_close(&Gab, h);

    }

    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
              ID("[o>o]-"), ID("[o>o]-"), 0, "Gamma <oo|oo>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gbb, h);
        dpd_buf4_mat_irrep_rd(&Gbb, h);

        #pragma omp parallel for
        for(long int ij = 0; ij < Gbb.params->rowtot[h]; ++ij){
            size_t i = Gbb.params->roworb[h][ij][0];
            int Gi = Gbb.params->psym[i];
            i -= Gbb.params->poff[Gi];
            size_t j = Gbb.params->roworb[h][ij][1];
            int Gj = Gbb.params->qsym[j];
            j -= Gbb.params->qoff[Gj];
            for(size_t kl = 0; kl < Gbb.params->coltot[h]; ++kl){
                double tpdm = 0.0;
                size_t k = Gbb.params->colorb[h][kl][0];
                int Gk = Gbb.params->rsym[k];
                k -= Gbb.params->roff[Gk];
                size_t l = Gbb.params->colorb[h][kl][1];
                int Gl = Gbb.params->ssym[l];
                l -= Gbb.params->soff[Gl];
                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * bkappa_->get(Gi, i, k) * bkappa_->get(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * bkappa_->get(Gi, i, l) * bkappa_->get(Gj, j, k);

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * bkappa_->get(Gi, i, k) * bocc_tau_->get(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * bkappa_->get(Gi, i, l) * bocc_tau_->get(Gj, j, k);
                if(Gj == Gk && Gi == Gl) tpdm -= 0.25 * bkappa_->get(Gj, j, k) * bocc_tau_->get(Gi, i, l);
                if(Gj == Gl && Gi == Gk) tpdm += 0.25 * bkappa_->get(Gj, j, l) * bocc_tau_->get(Gi, i, k);

                if(Gi == Gk && Gj == Gl) tpdm += 0.25 * bocc_tau_->get(Gi, i, k) * bocc_tau_->get(Gj, j, l);
                if(Gi == Gl && Gj == Gk) tpdm -= 0.25 * bocc_tau_->get(Gi, i, l) * bocc_tau_->get(Gj, j, k);

                Gbb.matrix[h][ij][kl] += tpdm;
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Gbb, h);

    }

    dpd_buf4_close(&Gbb);

    /*
     * The OOVV and VVOO blocks
     */

    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O>O]-"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_copy(&Laa, PSIF_DCFT_DENSITY, "Gamma <OO|VV>");
    dpd_buf4_sort(&Laa, PSIF_DCFT_DENSITY, rspq, ID("[V>V]-"), ID("[O>O]-"), "Gamma <VV|OO>");
    dpd_buf4_close(&Laa);
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O>O]-"), ID("[V>V]-"),
              ID("[O>O]-"), ID("[V>V]-"), 0, "Gamma <OO|VV>");
    dpd_buf4_scm(&Gaa, 0.5);
    dpd_buf4_close(&Gaa);
    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[V>V]-"), ID("[O>O]-"),
              ID("[V>V]-"), ID("[O>O]-"), 0, "Gamma <VV|OO>");
    dpd_buf4_scm(&Gaa, 0.5);
    dpd_buf4_close(&Gaa);

    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_copy(&Lab, PSIF_DCFT_DENSITY, "Gamma <Oo|Vv>");
    dpd_buf4_sort(&Lab, PSIF_DCFT_DENSITY, rspq, ID("[V,v]"), ID("[O,o]"), "Gamma <Vv|Oo>");
    dpd_buf4_close(&Lab);
    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
              ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");
    dpd_buf4_scm(&Gab, 0.5);
    dpd_buf4_close(&Gab);
    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"),
              ID("[V,v]"), ID("[O,o]"), 0, "Gamma <Vv|Oo>");
    dpd_buf4_scm(&Gab, 0.5);
    dpd_buf4_close(&Gab);

    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o>o]-"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_copy(&Lbb, PSIF_DCFT_DENSITY, "Gamma <oo|vv>");
    dpd_buf4_sort(&Lbb, PSIF_DCFT_DENSITY, rspq, ID("[v>v]-"), ID("[o>o]-"), "Gamma <vv|oo>");
    dpd_buf4_close(&Lbb);
    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o>o]-"), ID("[v>v]-"),
              ID("[o>o]-"), ID("[v>v]-"), 0, "Gamma <oo|vv>");
    dpd_buf4_scm(&Gbb, 0.5);
    dpd_buf4_close(&Gbb);
    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[v>v]-"), ID("[o>o]-"),
              ID("[v>v]-"), ID("[o>o]-"), 0, "Gamma <vv|oo>");
    dpd_buf4_scm(&Gbb, 0.5);
    dpd_buf4_close(&Gbb);

    /*
     * The OVOV block
     */

    // There are five unique spin cases: Г<IAJB>, Г<iajb>, Г<IaJb>, Г<iAjB>, Г<IajB>

    // Г<IAJB> spin case

    dpd_buf4_init(&Gaa, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma (OV|OV)");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    dpd_buf4_init(&LLaa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    dpd_contract444(&Laa, &LLaa, &Gaa, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Laa);
    dpd_buf4_close(&LLaa);
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_buf4_init(&LLab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_contract444(&Lab, &LLab, &Gaa, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&LLab);
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
        for(long int ia = 0; ia < Gaa.params->rowtot[h]; ++ia){
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
                    Gaa.matrix[h][ia][jb] += (akappa_->get(Gi, i, j) + aocc_tau_->get(Gi, i, j)) * avir_tau_->get(Ga, a, b);
                }
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gaa, h);
        dpd_buf4_mat_irrep_close(&Gaa, h);
    }

    dpd_buf4_close(&Gaa);

    // Г<IaJb> and Г<iAjB> spin cases:

    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Lambda (Ov|oV)");
    dpd_buf4_init(&LLab, PSIF_DCFT_DPD, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Lambda (Ov|oV)");
    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
    dpd_contract444(&Lab, &LLab, &Gab, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Gab);
    dpd_buf4_init(&Gba, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");
    dpd_contract444(&Lab, &LLab, &Gba, 1, 1, -1.0, 0.0);
    dpd_buf4_close(&Gba);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&LLab);

    dpd_buf4_init(&Gab, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gab, h);
        dpd_buf4_mat_irrep_rd(&Gab, h);

        #pragma omp parallel for
        for(long int ia = 0; ia < Gab.params->rowtot[h]; ++ia){
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
                    Gab.matrix[h][ia][jb] += (akappa_->get(Gi, i, j) + aocc_tau_->get(Gi, i, j)) * bvir_tau_->get(Ga, a, b);
                }
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
        for(long int ia = 0; ia < Gba.params->rowtot[h]; ++ia){
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
                    Gba.matrix[h][ia][jb] += (bkappa_->get(Gi, i, j) + bocc_tau_->get(Gi, i, j)) * avir_tau_->get(Ga, a, b);
                }

            }
        }
        dpd_buf4_mat_irrep_wrt(&Gba, h);
        dpd_buf4_mat_irrep_close(&Gba, h);
    }

    dpd_buf4_close(&Gba);

    // Г<IajB> spin case:

    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_buf4_init(&Laa, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Lambda (OV|OV)");
    dpd_contract444(&Laa, &Lab, &Tab, 0, 1, -1.0, 0.0);
    dpd_buf4_close(&Laa);
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    dpd_contract444(&Lab, &Lbb, &Tab, 0, 1, -1.0, 1.0);
    dpd_buf4_close(&Lbb);
    dpd_buf4_close(&Tab);
    dpd_buf4_init(&Tab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Temp (OV|ov)");
    dpd_buf4_sort(&Tab, PSIF_DCFT_DENSITY, psrq, ID("[O,v]"), ID("[o,V]"), "Gamma <Ov|oV>");
    // Resort to get the Г_oVOv. Used for the MO Lagrangian
    dpd_buf4_sort(&Tab, PSIF_DCFT_DENSITY, rqps, ID("[o,V]"), ID("[O,v]"), "Gamma <oV|Ov>");
    dpd_buf4_close(&Tab);
    dpd_buf4_close(&Lab);

    // Г<iajb> spin case:

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma (ov|ov)");
    dpd_buf4_init(&Lbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    dpd_buf4_init(&LLbb, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Lambda (ov|ov)");
    dpd_contract444(&Lbb, &LLbb, &Gbb, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&Lbb);
    dpd_buf4_close(&LLbb);
    dpd_buf4_init(&Lab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_buf4_init(&LLab, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "Lambda (OV|ov)");
    dpd_contract444(&Lab, &LLab, &Gbb, 1, 1, -1.0, 1.0);
    dpd_buf4_close(&Lab);
    dpd_buf4_close(&LLab);
    dpd_buf4_close(&Gbb);

    // Resort Г(ov|ov) to the Г<ov|ov>
    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma (ov|ov)");
    dpd_buf4_sort(&Gbb, PSIF_DCFT_DENSITY, psrq, ID("[o,v]"),ID("[o,v]"), "Gamma <ov|ov>");
    dpd_buf4_close(&Gbb);

    dpd_buf4_init(&Gbb, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&Gbb, h);
        dpd_buf4_mat_irrep_rd(&Gbb, h);

        #pragma omp parallel for
        for(long int ia = 0; ia < Gbb.params->rowtot[h]; ++ia){
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
                    Gbb.matrix[h][ia][jb] += (bkappa_->get(Gi, i, j) + bocc_tau_->get(Gi, i, j)) * bvir_tau_->get(Ga, a, b);
                }
            }
        }
        dpd_buf4_mat_irrep_wrt(&Gbb, h);
        dpd_buf4_mat_irrep_close(&Gbb, h);
    }

    dpd_buf4_close(&Gbb);

    psio_->close(PSIF_DCFT_DENSITY, 1);

}

void
DCFTSolver::compute_orbital_gradient_OV() {

    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 G, I;
    dpdbuf4 L, W, LL;
    dpdfile2 X, H, T;
    dpdfile2 T_VV, T_vv;
    dpdfile2 Y2_OV, Y2_ov;

    // X_OV: One-electron contributions

    // X_IA = H_IB Tau_BA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    dpd_file2_init(&T, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    dpd_contract222(&H, &T, &X, 0, 1, 1.0, 0.0);
    dpd_file2_close(&T);
    dpd_file2_close(&H);
    dpd_file2_close(&X);

    // X_ia = H_ib Tau_ba
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    dpd_file2_init(&T, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");
    dpd_contract222(&H, &T, &X, 0, 1, 1.0, 0.0);
    dpd_file2_close(&T);
    dpd_file2_close(&H);
    dpd_file2_close(&X);

    // X_OV: Two-electron contributions

    //
    // 2 * <OV||VV> Г_VVVV
    //

    // Compute contributions from VVVV density
    // 1. X_ia <-- <ib||cd> tau_ca tau_db
    dpd_file2_init(&T_VV, PSIF_DCFT_DPD, 0, ID('V'), ID('V'), "Tau <V|V>");
    dpd_file2_init(&T_vv, PSIF_DCFT_DPD, 0, ID('v'), ID('v'), "Tau <v|v>");

    // Alpha contribution X_IA
    dpd_file2_init(&Y2_OV, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "Y2 <O|V>");

    // Y2_IA = <IA|CD> Tau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
    dpd_contract422(&I, &T_VV, &Y2_OV, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // Y2_IA -= (IA|CD) Tau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
    dpd_contract422(&I, &T_VV, &Y2_OV, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y2_IA -= (IA|cd) Tau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                  ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
    dpd_contract422(&I, &T_vv, &Y2_OV, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);

    // X_IA -= Y2_IC Tau_CA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_contract222(&Y2_OV, &T_VV, &X, 0, 1, -1.0, 1.0);
    dpd_file2_close(&X);

    dpd_file2_close(&Y2_OV);

    // Beta contribution X_ia
    dpd_file2_init(&Y2_ov, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "Y2 <o|v>");

    // Y2_ib = <ib|cd> Tau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
    dpd_contract422(&I, &T_vv, &Y2_ov, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&I);
    // Y2_ib -= (ib|cd) Tau_cd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v>=v]+"), 0, "MO Ints (ov|vv)");
    dpd_contract422(&I, &T_vv, &Y2_ov, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);
    // Y2_ib -= (ib|CD) Tau_CD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[V>=V]+"),
                  ID("[o,v]"), ID("[V>=V]+"), 0, "MO Ints (ov|VV)");
    dpd_contract422(&I, &T_VV, &Y2_ov, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&I);

    // X_ia -= Y2_ic Tau_ca
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_contract222(&Y2_ov, &T_vv, &X, 0, 1, -1.0, 1.0);
    dpd_file2_close(&X);

    dpd_file2_close(&Y2_ov);

    dpd_file2_close(&T_VV);
    dpd_file2_close(&T_vv);

    // 2. X_ia <-- 1/4 <ib||cd> lambda_abkl lambda_klcd

    //  W_IBKL = <IB||CD> lambda_KLCD
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V>V]-"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V>V]-"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O>O]-"), 0, "W <OV|OO>");
    dpd_contract444(&I, &L, &W, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&W);


    //  W_KlIb = 2 lambda_KlCd <Ib|Cd>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "W <Oo|Ov>");
    dpd_contract444(&L, &I, &W, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&W);

    //  W_LkBi = 2 lambda_LkDc <Bi|Dc>
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "W <Oo|Vo>");
    dpd_contract444(&L, &I, &W, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&W);

    //  W_ibkl = <ib||cd> lambda_klcd
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v>v]-"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    dpd_buf4_init(&L, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v>v]-"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o>o]-"), 0, "W <ov|oo>");
    dpd_contract444(&I, &L, &W, 0, 0, 2.0, 0.0);
    dpd_buf4_close(&I);
    dpd_buf4_close(&L);
    dpd_buf4_close(&W);


    // X_IA +=  1/4 W_IBKL lambda_KLAB
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O>O]-"), 0, "W <OV|OO>");
    dpd_buf4_init(&LL, PSIF_DCFT_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Lambda <OO|VV>");

    dpd_contract442(&W, &LL, &X, 0, 2, 0.25, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&LL);
    dpd_file2_close(&X);

    // X_IA +=  1/2 W_KlIb lambda_KlAb
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "W <Oo|Ov>");
    dpd_buf4_init(&LL, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");

    dpd_contract442(&W, &LL, &X, 2, 2, 0.5, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&LL);
    dpd_file2_close(&X);

    // X_ia +=  1/2 W_LkBi lambda_LkBa
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "W <Oo|Vo>");
    dpd_buf4_init(&LL, PSIF_DCFT_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Lambda <Oo|Vv>");

    dpd_contract442(&W, &LL, &X, 3, 3, 0.5, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&LL);
    dpd_file2_close(&X);

    // X_ia += 1/4 W_ibkl lambda_klab
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&W, PSIF_DCFT_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o>o]-"), 0, "W <ov|oo>");
    dpd_buf4_init(&LL, PSIF_DCFT_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Lambda <oo|vv>");

    dpd_contract442(&W, &LL, &X, 0, 2, 0.25, 1.0);
    dpd_buf4_close(&W);
    dpd_buf4_close(&LL);
    dpd_file2_close(&X);

    //
    // <OO||OV> Г_OOVV
    //

    // X_IA += <BI||JK> Г_BAJK
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 1, "MO Ints <VO|OO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V>V]-"), ID("[O>O]-"), 0, "Gamma <VV|OO>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_IA += 2 * <Ib|Jk> Г_AbJk
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"),
                  ID("[V,v]"), ID("[O,o]"), 0, "Gamma <Vv|Oo>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ia += <bi||jk> Г_bajk
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 1, "MO Ints <vo|oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v>v]-"), ID("[o>o]-"), 0, "Gamma <vv|oo>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ia += 2 * <Bi|Jk> Г_BaJk
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "MO Ints <Vo|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[V,v]"), ID("[O,o]"),
                  ID("[V,v]"), ID("[O,o]"), 0, "Gamma <Vv|Oo>");

    dpd_contract442(&I, &G, &X, 1, 1, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    //
    // <OO||OV> Г_OVOV
    //

    // X_IA += <JB||KI> Г_JBKA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 1, "MO Ints <OV|OO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_IA += <kB|jI> Г_kBjA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[o,O]"),
                  ID("[o,V]"), ID("[o,O]"), 0, "MO Ints <oV|oO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_IA -= <Kb|jI> Г_KbjA
    // Note: <Kb|jI> integrals are resorted <bK|jI> integrals.
    // <Kb||jI> Г_KbjA = (-1) * <bK|jI> * Г_KbjA
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('O'), ID('V'), "X <O|V>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[o,O]"),
                  ID("[O,v]"), ID("[o,O]"), 0, "MO Ints <Ov|oO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");

    dpd_contract442(&I, &G, &X, 3, 3, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ia += <jb||ki> Г_jbka
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o,o]"), 1, "MO Ints <ov|oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ia += <Kb|Ji> Г_KbJa
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");

    dpd_contract442(&I, &G, &X, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ia -= <kB|Ji> Г_kBJa
    // Note: <kB|Ji> integrals are resorted <Bk|Ji> integrals.
    // <kB||Ji> Г_kBJa = (-1) * <Bk|Ji> * Г_kBJa

    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('o'), ID('v'), "X <o|v>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[O,o]"),
                  ID("[o,V]"), ID("[O,o]"), 0, "MO Ints <oV|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "Gamma <oV|Ov>");

    dpd_contract442(&I, &G, &X, 3, 3, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    psio_->close(PSIF_DCFT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);

}

void
DCFTSolver::compute_orbital_gradient_VO() {

    psio_->open(PSIF_DCFT_DENSITY, PSIO_OPEN_OLD);
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpdbuf4 G, I;
    dpdfile2 X, H, T;

    // X_VO: One-electron contributions

    // X_AI = H_JA Tau_JI
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    dpd_file2_init(&T, PSIF_DCFT_DPD, 0, ID('O'), ID('O'), "Tau <O|O>");
    dpd_file2_mat_init(&X);
    dpd_file2_mat_init(&H);
    dpd_file2_mat_init(&T);
    dpd_file2_mat_rd(&H);
    dpd_file2_mat_rd(&T);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0 ; i < naoccpi_[h]; ++i){
            for(int a = 0 ; a < navirpi_[h]; ++a){
                double value = 0.0;
                for(int j = 0 ; j < naoccpi_[h]; ++j){
                    value += H.matrix[h][j][a] * (T.matrix[h][j][i] + (i == j ? 1.0 : 0.0));
                }
                X.matrix[h][a][i] = value;
            }
        }
    }
    dpd_file2_mat_wrt(&X);
    dpd_file2_close(&T);
    dpd_file2_close(&H);
    dpd_file2_close(&X);

    // X_ai = H_ja Tau_ji
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    dpd_file2_init(&T, PSIF_DCFT_DPD, 0, ID('o'), ID('o'), "Tau <o|o>");
    dpd_file2_mat_init(&X);
    dpd_file2_mat_init(&H);
    dpd_file2_mat_init(&T);
    dpd_file2_mat_rd(&H);
    dpd_file2_mat_rd(&T);
    for(int h = 0; h < nirrep_; ++h){
        #pragma omp parallel for
        for(int i = 0 ; i < nboccpi_[h]; ++i){
            for(int a = 0 ; a < nbvirpi_[h]; ++a){
                double value = 0.0;
                for(int j = 0 ; j < nboccpi_[h]; ++j){
                    value += H.matrix[h][j][a] * (T.matrix[h][j][i] + (i == j ? 1.0 : 0.0));
                }
                X.matrix[h][a][i] = value;
            }
        }
    }
    dpd_file2_mat_wrt(&X);
    dpd_file2_close(&T);
    dpd_file2_close(&H);
    dpd_file2_close(&X);

    // X_VO: Two-electron contributions

    //
    // 2 * <VO||OO> Г_OOOO
    //

    // X_AI += 2 * <AJ||KL> Г_IJKL
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 1, "MO Ints <VO|OO>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O>O]-"), ID("[O>O]-"), 0, "Gamma <OO|OO>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_AI += 4 * <Aj|Kl> Г_IjKl
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "MO Ints <Vo|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");

    dpd_contract442(&I, &G, &X, 0, 0, 4.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ai += 2 * <aj||kl> Г_ijkl
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 1, "MO Ints <vo|oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o>o]-"), ID("[o>o]-"), 0, "Gamma <oo|oo>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_AI += 4 * <Ja|Kl> Г_JiKl
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "MO Ints <Ov|Oo>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "Gamma <Oo|Oo>");

    dpd_contract442(&I, &G, &X, 1, 1, 4.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    //
    // <VO||VV> Г_OOVV
    //

    // X_AI += <JA||BC> Г_JIBC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>O]-"), ID("[V>V]-"), 0, "Gamma <OO|VV>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_AI += <Aj|Bc> Г_IjBc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");

    dpd_contract442(&I, &G, &X, 0, 0, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ai += <ja||bc> Г_jibc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o>o]-"), ID("[v>v]-"), 0, "Gamma <oo|vv>");

    dpd_contract442(&I, &G, &X, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ai += <Ja|Bc> Г_JiBc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "Gamma <Oo|Vv>");

    dpd_contract442(&I, &G, &X, 1, 1, 2.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    //
    // <OV||VV> Г_OVOV
    //

    // X_AI += <JB||AC> Г_JBIC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 1, "MO Ints <OV|VV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "Gamma <OV|OV>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_AI += <Jb|Ac> Г_JbIc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "Gamma <Ov|Ov>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_AI -= <jB|Ac> Г_jBIc
    // Note: <jB|Ac> integrals are resorted <Bj|Ac> integrals.
    // <jB||Ac> Г_jBIc = (-1) * <Bj|Ac> Г_jBIc
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('V'), ID('O'), "X <V|O>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[V,v]"),
                  ID("[o,V]"), ID("[V,v]"), 0, "MO Ints <oV|Vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[O,v]"),
                  ID("[o,V]"), ID("[O,v]"), 0, "Gamma <oV|Ov>");

    dpd_contract442(&I, &G, &X, 2, 2, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ai += <jb||ac> Г_jbic
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 1, "MO Ints <ov|vv>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "Gamma <ov|ov>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ai += <jB|aC> Г_jBiC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[o,V]"), ID("[v,V]"),
                  ID("[o,V]"), ID("[v,V]"), 0, "MO Ints <oV|vV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[o,V]"), ID("[o,V]"),
                  ID("[o,V]"), ID("[o,V]"), 0, "Gamma <oV|oV>");

    dpd_contract442(&I, &G, &X, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    // X_ai -= <Jb|aC> Г_JbiC
    // Note: <Jb|aC> integrals are resorted <bJ|aC> integrals.
    // <Jb||aC> Г_JbiC = (-1) * <bJ|aC> Г_JbiC
    dpd_file2_init(&X, PSIF_DCFT_DPD, 0, ID('v'), ID('o'), "X <v|o>");
    dpd_buf4_init(&I, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[v,V]"),
                  ID("[O,v]"), ID("[v,V]"), 0, "MO Ints <Ov|vV>");
    dpd_buf4_init(&G, PSIF_DCFT_DENSITY, 0, ID("[O,v]"), ID("[o,V]"),
                  ID("[O,v]"), ID("[o,V]"), 0, "Gamma <Ov|oV>");

    dpd_contract442(&I, &G, &X, 2, 2, -1.0, 1.0);
    dpd_buf4_close(&G);
    dpd_buf4_close(&I);
    dpd_file2_close(&X);

    psio_->close(PSIF_DCFT_DENSITY, 1);
    psio_->close(PSIF_LIBTRANS_DPD, 1);
}

void
DCFTSolver::update_orbitals_jacobi()
{

    // Initialize the orbital rotation matrix
    SharedMatrix U_a(new Matrix("Orbital rotation matrix (Alpha)", nirrep_, nmopi_, nmopi_));
    SharedMatrix U_b(new Matrix("Orbital rotation matrix (Beta)", nirrep_, nmopi_, nmopi_));

    // Compute the orbital rotation matrix and rotate the orbitals

    // U = I
    U_a->identity();
    U_b->identity();

    // U += X
    U_a->add(Xtotal_a_);
    U_b->add(Xtotal_b_);

    // U += 0.5 * X * X
    U_a->gemm(false, false, 0.5, Xtotal_a_, Xtotal_a_, 1.0);
    U_b->gemm(false, false, 0.5, Xtotal_b_, Xtotal_b_, 1.0);

    // Orthogonalize the U vectors
    int rowA = U_a->nrow();
    int colA = U_a->ncol();

    double **U_a_block = block_matrix(rowA, colA);
    memset(U_a_block[0], 0, sizeof(double)*rowA*colA);
    U_a_block = U_a->to_block_matrix();
    schmidt(U_a_block, rowA, colA, outfile);
    U_a->set(U_a_block);
    free_block(U_a_block);

    int rowB = U_b->nrow();
    int colB = U_b->ncol();

    double **U_b_block = block_matrix(rowB, colB);
    memset(U_b_block[0], 0, sizeof(double)*rowB*colB);
    U_b_block = U_b->to_block_matrix();
    schmidt(U_b_block, rowB, colB, outfile);
    U_b->set(U_b_block);
    free_block(U_b_block);

    // Rotate the orbitals
    Ca_->gemm(false, false, 1.0, old_ca_, U_a, 0.0);
    Cb_->gemm(false, false, 1.0, old_cb_, U_b, 0.0);


}

}} //End namespaces

