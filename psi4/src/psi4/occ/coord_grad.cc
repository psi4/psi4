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

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libpsio/psio.hpp"
#include "occwave.h"
#include "defines.h"


using namespace std;

namespace psi{ namespace occwave{

void OCCWave::coord_grad()
{
      if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
          outfile->Printf("\tComputing G_abcd...\n");

          omp3_tpdm_vvvv();
      }
      else if (wfn_type_ == "OCEPA") {
          outfile->Printf("\tComputing G_abcd...\n");

          ocepa_tpdm_vvvv();
      }
      outfile->Printf("\tComputing diagonal blocks of GFM...\n");

      gfock_diag();

      // For Standard methods
      if (orb_opt_ == "FALSE" && relaxed_ == "TRUE") {
          outfile->Printf("\tSolving orbital Z-vector equations...\n");

          z_vector();
          outfile->Printf("\tForming relaxed response density matrices...\n");

          effective_pdms();
          outfile->Printf("\tForming relaxed GFM...\n");

          effective_gfock();
      }

     // OEPROP
     if (oeprop_ == "TRUE") oeprop();

      dump_ints();
      outfile->Printf("\tWriting particle density matrices and GFM to disk...\n");

      dump_pdms();
}//

//========================================================================
//         Dump Molecular Integrals
//========================================================================
void OCCWave::dump_ints()
{
    //outfile->Printf("\n dump_ints is starting... \n");
    dpdfile2 H;
    dpdbuf4 K;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    global_dpd_->file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
                H.matrix[h][i][j] = HmoA->get(h, i, j);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    global_dpd_->file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int b = 0 ; b < virtpiA[h]; ++b){
                H.matrix[h][a][b] = HmoA->get(h, a + occpiA[h], b + occpiA[h]);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    global_dpd_->file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
                H.matrix[h][i][j] = HmoA->get(h, i, j + occpiA[h]);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_close(&H);

 // Write beta-integrals if ref is UHF
 if (reference_ == "UNRESTRICTED") {
    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "H <o|o>");
    global_dpd_->file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < occpiB[h]; ++j){
                H.matrix[h][i][j] = HmoB->get(h, i, j);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "H <v|v>");
    global_dpd_->file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int a = 0 ; a < virtpiB[h]; ++a){
            for(int b = 0 ; b < virtpiB[h]; ++b){
                H.matrix[h][a][b] = HmoB->get(h, a + occpiB[h], b + occpiB[h]);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    global_dpd_->file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < virtpiB[h]; ++j){
                H.matrix[h][i][j] = HmoB->get(h, i, j + occpiB[h]);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_close(&H);
 }// end uhf

    psio_->close(PSIF_LIBTRANS_DPD, 1);
    //outfile->Printf("\n dump_ints done. \n");

}// end of dump_ints

//========================================================================
//         Dump PDMs
//========================================================================
void OCCWave::dump_pdms()
{
    //outfile->Printf("\n dump_pdms is starting... \n");

//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {
    dpdfile2 H;
    dpdbuf4 G, G2;

    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

    const int *alpha_corr_to_pitzer = ints->alpha_corr_to_pitzer();
    int *alpha_pitzer_to_corr = new int[nmo_];
    memset(alpha_pitzer_to_corr, 0, nmo_*sizeof(int));

    for(int n = 0; n < nmo_; ++n) {
        alpha_pitzer_to_corr[alpha_corr_to_pitzer[n]] = n;
    }


    // Reorder the one-particle density matrix to the QT order
    double **a_qt = block_matrix(nmo_, nmo_);

    int offset = 0;
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = alpha_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = alpha_pitzer_to_corr[pitzer_j];
                a_qt[corr_i][corr_j] = g1symm->get(h, i, j);
            }
        }
        offset += nmopi_[h];
    }
   //print_mat(a_qt, nmo_, nmo_, outfile);

   // Write qt-ordered OPDM to the file
   psio_->open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
   psio_->write_entry(PSIF_MO_OPDM, "MO-basis OPDM", (char *) a_qt[0], sizeof(double) * nmo_ * nmo_);
   psio_->close(PSIF_MO_OPDM, 1);


    // Scale the generalized-Fock matrix by 2.0 to make it the same form as in the coupled-cluster code
    GFock->scale(2.0);

    // Reorder the Generalized-Fock matrix to the QT order
    offset = 0;
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = alpha_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = alpha_pitzer_to_corr[pitzer_j];
                a_qt[corr_i][corr_j] = GFock->get(h, i, j);
            }
        }
        offset += nmopi_[h];
    }

    // Write qt-ordered Generalized-Fock matrix to the file
    psio_->open(PSIF_MO_LAG, PSIO_OPEN_OLD);
    psio_->write_entry(PSIF_MO_LAG, "MO-basis Lagrangian", (char *) a_qt[0], sizeof(double) * nmo_ * nmo_);
    psio_->close(PSIF_MO_LAG, 1);
    free_block(a_qt);

    int *aocc_qt = new int[nooA];
    int *avir_qt = new int[nvoA];

    int aocc_count = 0;
    int avir_count = 0;
    offset = 0;

    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < occpiA[h]; ++i) {
            int pitzer = offset + i;
            aocc_qt[aocc_count++] = alpha_pitzer_to_corr[pitzer];
        }
        for (int i = occpiA[h]; i < nmopi_[h]; ++i) {
            int pitzer = offset + i;
            avir_qt[avir_count++] = alpha_pitzer_to_corr[pitzer];
        }
        offset += nmopi_[h];
    }


    // Write TPDMs to psi files
    struct iwlbuf AA;
    iwl_buf_init(&AA, PSIF_MO_TPDM, 1.0E-15, 0, 0);

    // OOOO block
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
              ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(size_t ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            int I = aocc_qt[i];
            int J = aocc_qt[j];
            for(size_t kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
                int K = aocc_qt[k];
                int L = aocc_qt[l];
                int IK = index2(I,K);
                int JL = index2(J,L);
                if (I >= K && J >= L && IK >= JL) {
                   double value = 2.0 * G.matrix[h][ij][kl];
                   if (I > K) value *= 2.0;
                   if (J > L) value *= 2.0;
                   iwl_buf_wrt_val(&AA, I, K, J, L, value, 0, "NULL", 0);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

if (wfn_type_ != "OMP2") {
    // VVVV block
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
              ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(size_t ab = 0; ab < G.params->rowtot[h]; ++ab){
            int a = G.params->roworb[h][ab][0];
            int b = G.params->roworb[h][ab][1];
            int A = avir_qt[a];
            int B = avir_qt[b];
            for(size_t cd = 0; cd < G.params->coltot[h]; ++cd){
                int c = G.params->colorb[h][cd][0];
                int d = G.params->colorb[h][cd][1];
                int C = avir_qt[c];
                int D = avir_qt[d];
                int AC = index2(A,C);
                int BD = index2(B,D);
                if (A >= C && B >= D && AC >= BD) {
                   double value = 2.0 * G.matrix[h][ab][cd];
                   if (A > C) value *= 2.0;
                   if (B > D) value *= 2.0;
                   iwl_buf_wrt_val(&AA, A, C, B, D, value, 0, "NULL", 0);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);
}// end if (wfn_type_ != "OMP2") {

    // OOVV
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
              ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            int I = aocc_qt[i];
            int J = aocc_qt[j];
            for(int ab = 0; ab < G.params->coltot[h]; ++ab){
                int a = G.params->colorb[h][ab][0];
                int b = G.params->colorb[h][ab][1];
                int A = avir_qt[a];
                int B = avir_qt[b];
                int IA = ((A - nooA) * nooA)  + I;
                int JB = ((B - nooA) * nooA)  + J;
                if (IA >= JB) {
                   double value = 8.0 * G.matrix[h][ij][ab];
                   iwl_buf_wrt_val(&AA, A, I, B, J, value, 0, "NULL", 0);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // OVOV
    // Alpha-Alpha spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
              ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    global_dpd_->buf4_scm(&G, 2.0);
    global_dpd_->buf4_dump(&G, &AA, aocc_qt, avir_qt, aocc_qt, avir_qt, 0, 1);
    global_dpd_->buf4_close(&G);

   // For the standard methods I need the following contribution
   if (orb_opt_ == "FALSE") {
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "TPDM <VO|OO>");
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_dump(&G, &AA, avir_qt, aocc_qt, aocc_qt, aocc_qt, 0, 1);
    global_dpd_->buf4_close(&G);
   }// if (orb_opt_ == "FALSE")

    delete [] aocc_qt;
    delete [] avir_qt;
    delete[] alpha_pitzer_to_corr;

    psio_->close(PSIF_OCC_DENSITY, 1);

    iwl_buf_flush(&AA, 1);
    iwl_buf_close(&AA, 1);

}// end if (reference_ == "RESTRICTED")



//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference_ == "UNRESTRICTED") {

    dpdfile2 H;
    dpdbuf4 G, G2;

    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

    const int *alpha_corr_to_pitzer = ints->alpha_corr_to_pitzer();
    int *alpha_pitzer_to_corr = new int[nmo_];
    memset(alpha_pitzer_to_corr, 0, nmo_*sizeof(int));

    for(int n = 0; n < nmo_; ++n) {
        alpha_pitzer_to_corr[alpha_corr_to_pitzer[n]] = n;
    }


    const int *beta_corr_to_pitzer = ints->beta_corr_to_pitzer();
    int *beta_pitzer_to_corr = new int[nmo_];
    memset(beta_pitzer_to_corr, 0, nmo_*sizeof(int));

    for(int n = 0; n < nmo_; ++n) {
        beta_pitzer_to_corr[beta_corr_to_pitzer[n]] = n;
    }


    // Reorder the one-particle density matrix to the QT order
    double **a_qt = block_matrix(nmo_, nmo_);
    double **b_qt = block_matrix(nmo_, nmo_);

    int offset = 0;
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = alpha_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = alpha_pitzer_to_corr[pitzer_j];
                a_qt[corr_i][corr_j] = g1symmA->get(h, i, j);
            }
        }
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = beta_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = beta_pitzer_to_corr[pitzer_j];
                b_qt[corr_i][corr_j] = g1symmB->get(h, i, j);
            }
        }
        offset += nmopi_[h];
    }

   // Write qt-ordered OPDM to the file
   psio_->open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
   psio_->write_entry(PSIF_MO_OPDM, "MO-basis Alpha OPDM", (char *) a_qt[0], sizeof(double) * nmo_ * nmo_);
   psio_->write_entry(PSIF_MO_OPDM, "MO-basis Beta OPDM", (char *) b_qt[0], sizeof(double) * nmo_ * nmo_);
   psio_->close(PSIF_MO_OPDM, 1);


    // Scale the generalized-Fock matrix by 2.0 to make it the same form as in the coupled-cluster code
    GFockA->scale(2.0);
    GFockB->scale(2.0);

    // Reorder the Generalized-Fock matrix to the QT order
    offset = 0;
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = alpha_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = alpha_pitzer_to_corr[pitzer_j];
                a_qt[corr_i][corr_j] = GFockA->get(h, i, j);
            }
        }
        for(int i = 0 ; i < nmopi_[h]; ++i){
            int pitzer_i = i + offset;
            int corr_i = beta_pitzer_to_corr[pitzer_i];
            for(int j = 0 ; j < nmopi_[h]; ++j){
                int pitzer_j = j + offset;
                int corr_j = beta_pitzer_to_corr[pitzer_j];
                b_qt[corr_i][corr_j] = GFockB->get(h, i, j);
            }
        }
        offset += nmopi_[h];
    }

    // Write qt-ordered Generalized-Fock matrix to the file
    psio_->open(PSIF_MO_LAG, PSIO_OPEN_OLD);
    psio_->write_entry(PSIF_MO_LAG, "MO-basis Alpha Lagrangian", (char *) a_qt[0], sizeof(double) * nmo_ * nmo_);
    psio_->write_entry(PSIF_MO_LAG, "MO-basis Beta Lagrangian", (char *) b_qt[0], sizeof(double) * nmo_ * nmo_);
    psio_->close(PSIF_MO_LAG, 1);
    free_block(a_qt);
    free_block(b_qt);

    int *aocc_qt = new int[nooA];
    int *bocc_qt = new int[nooB];
    int *avir_qt = new int[nvoA];
    int *bvir_qt = new int[nvoB];

    int aocc_count = 0;
    int bocc_count = 0;
    int avir_count = 0;
    int bvir_count = 0;
    offset = 0;

    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < occpiA[h]; ++i) {
            int pitzer = offset + i;
            aocc_qt[aocc_count++] = alpha_pitzer_to_corr[pitzer];
        }
        for (int i = 0; i < occpiB[h]; ++i) {
            int pitzer = offset + i;
            bocc_qt[bocc_count++] = beta_pitzer_to_corr[pitzer];
        }
        for (int i = occpiA[h]; i < nmopi_[h]; ++i) {
            int pitzer = offset + i;
            avir_qt[avir_count++] = alpha_pitzer_to_corr[pitzer];
        }
        for (int i = occpiB[h]; i < nmopi_[h]; ++i) {
            int pitzer = offset + i;
            bvir_qt[bvir_count++] = beta_pitzer_to_corr[pitzer];
        }
        offset += nmopi_[h];
    }


    // Write TPDMs to psi files
    struct iwlbuf AA, AB, BB;
    iwl_buf_init(&AA, PSIF_MO_AA_TPDM, 1.0E-15, 0, 0);
    iwl_buf_init(&AB, PSIF_MO_AB_TPDM, 1.0E-15, 0, 0);
    iwl_buf_init(&BB, PSIF_MO_BB_TPDM, 1.0E-15, 0, 0);

    // OOOO: Alpha-Alpha spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
              ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
    // Dump tpdm in chemist notation, 0 means no pack, 1 means swap indices 23
    global_dpd_->buf4_dump(&G, &AA, aocc_qt, aocc_qt, aocc_qt, aocc_qt, 0, 1);
    global_dpd_->buf4_close(&G);

    // OOOO: Beta-Beta spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
              ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
    global_dpd_->buf4_dump(&G, &BB, bocc_qt, bocc_qt, bocc_qt, bocc_qt, 0, 1);
    global_dpd_->buf4_close(&G);

    // OOOO: Alpha-Beta spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
              ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(size_t ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            int I = aocc_qt[i];
            int J = bocc_qt[j];
            for(size_t kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
                int K = aocc_qt[k];
                int L = bocc_qt[l];
                double value = 4.0 * G.matrix[h][ij][kl];
                iwl_buf_wrt_val(&AB, I, K, J, L, value, 0, "NULL", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

//if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") {
if (wfn_type_ != "OMP2") {
    // VVVV: Alpha-Alpha spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
              ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
    global_dpd_->buf4_dump(&G, &AA, avir_qt, avir_qt, avir_qt, avir_qt, 0, 1);
    global_dpd_->buf4_close(&G);

    // VVVV: Beta-Beta spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
              ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
    global_dpd_->buf4_dump(&G, &BB, bvir_qt, bvir_qt, bvir_qt, bvir_qt, 0, 1);
    global_dpd_->buf4_close(&G);

    // VVVV: Alpha-Beta spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
              ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
    global_dpd_->buf4_sort(&G, PSIF_OCC_DENSITY , prqs, ID("[V,V]"), ID("[v,v]"), "TPDM (VV|vv)");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[v,v]"),
              ID("[V,V]"), ID("[v,v]"), 0, "TPDM (VV|vv)");
    global_dpd_->buf4_scm(&G, 4.0);
    global_dpd_->buf4_dump(&G, &AB, avir_qt, avir_qt, bvir_qt, bvir_qt, 0, 0);
    global_dpd_->buf4_close(&G);
}// end if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") {

    // OOVV: Alpha-Alpha spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
              ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    global_dpd_->buf4_scm(&G, 2.0);
    global_dpd_->buf4_dump(&G, &AA, aocc_qt, aocc_qt, avir_qt, avir_qt, 0, 1);
    global_dpd_->buf4_close(&G);

    // OOVV: Beta-Beta spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
              ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
    global_dpd_->buf4_scm(&G, 2.0);
    global_dpd_->buf4_dump(&G, &BB, bocc_qt, bocc_qt, bvir_qt, bvir_qt, 0, 1);
    global_dpd_->buf4_close(&G);

    // OOVV: Alpha-Beta spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
              ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(size_t ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            int I = aocc_qt[i];
            int J = bocc_qt[j];
            for(size_t ab = 0; ab < G.params->coltot[h]; ++ab){
                int a = G.params->colorb[h][ab][0];
                int b = G.params->colorb[h][ab][1];
                int A = avir_qt[a];
                int B = bvir_qt[b];
                double value = 8.0 * G.matrix[h][ij][ab];
                iwl_buf_wrt_val(&AB, I, A, J, B, value, 0, "NULL", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);


    // OVOV: Alpha-Alpha spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
              ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = aocc_qt[i];
            int A = avir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = aocc_qt[j];
                int B = avir_qt[b];
                double value = 2.0 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AA, I, J, A, B, value, 0, "NULL", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // OVOV: Beta-Beta spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
              ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = bocc_qt[i];
            int A = bvir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = bocc_qt[j];
                int B = bvir_qt[b];
                double value = 2.0 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&BB, I, J, A, B, value, 0, "NULL", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // OVOV: Alpha-Beta spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(int ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = aocc_qt[i];
            int A = bvir_qt[a];
            for(int jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = aocc_qt[j];
                int B = bvir_qt[b];
                double value = 4.0 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AB, I, J, A, B, value, 0, "NULL", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);



    // OVVO: Alpha-Alpha spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
              ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = aocc_qt[i];
            int A = avir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = aocc_qt[j];
                int B = avir_qt[b];
                double value = -2.0 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&AA, I, B, A, J, value, 0, "NULL", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);


    // OVVO: Beta-Beta spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
              ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = bocc_qt[i];
            int A = bvir_qt[a];
            for(size_t jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
                int J = bocc_qt[j];
                int B = bvir_qt[b];
                double value = -2.0 * G.matrix[h][ia][jb];
                iwl_buf_wrt_val(&BB, I, B, A, J, value, 0, "NULL", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

//if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") {
if (wfn_type_ != "OMP2") {
    // OVVO: Alpha-Beta spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
              ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(size_t ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            int I = aocc_qt[i];
            int A = bvir_qt[a];
            for(size_t bj = 0; bj < G.params->coltot[h]; ++bj){
                int b = G.params->colorb[h][bj][0];
                int j = G.params->colorb[h][bj][1];
                int B = avir_qt[b];
                int J = bocc_qt[j];
                double value = 8.0 * G.matrix[h][ia][bj];
                iwl_buf_wrt_val(&AB, I, B, A, J, value, 0, "NULL", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);
}// end if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") {


    // VOVO: Alpha-Beta spin-case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(int ai = 0; ai < G.params->rowtot[h]; ++ai){
            int a = G.params->roworb[h][ai][0];
            int i = G.params->roworb[h][ai][1];
            int A = avir_qt[a];
            int I = bocc_qt[i];
            for(int bj = 0; bj < G.params->coltot[h]; ++bj){
                int b = G.params->colorb[h][bj][0];
                int j = G.params->colorb[h][bj][1];
                int B = avir_qt[b];
                int J = bocc_qt[j];
                double value = 4.0 * G.matrix[h][ai][bj];
                iwl_buf_wrt_val(&AB, A, B, I, J, value, 0, "NULL", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

   // For the standard methods I need the following contribution
   if (orb_opt_ == "FALSE") {
    // VOOO: Alpha-Alpha Spin-Case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "TPDM <VO|OO>");
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_dump(&G, &AA, avir_qt, aocc_qt, aocc_qt, aocc_qt, 0, 1);
    global_dpd_->buf4_close(&G);

    // VOOO: Beta-Beta Spin-Case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "TPDM <vo|oo>");
    global_dpd_->buf4_scm(&G, 0.5);
    global_dpd_->buf4_dump(&G, &BB, bvir_qt, bocc_qt, bocc_qt, bocc_qt, 0, 1);
    global_dpd_->buf4_close(&G);

    // VOOO: Alpha-Beta Spin-Case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "TPDM <Vo|Oo>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(int am = 0; am < G.params->rowtot[h]; ++am){
            int a = G.params->roworb[h][am][0];
            int m = G.params->roworb[h][am][1];
            int A = avir_qt[a];
            int M = bocc_qt[m];
            for(int in = 0; in < G.params->coltot[h]; ++in){
                int i = G.params->colorb[h][in][0];
                int n = G.params->colorb[h][in][1];
                int I = aocc_qt[i];
                int N = bocc_qt[n];
                double value = G.matrix[h][am][in];
                iwl_buf_wrt_val(&AB, A, I, M, N, value, 0, "NULL", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // OVOO: Alpha-Beta Spin-Case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "TPDM <Ov|Oo>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(int ma = 0; ma < G.params->rowtot[h]; ++ma){
            int m = G.params->roworb[h][ma][0];
            int a = G.params->roworb[h][ma][1];
            int M = aocc_qt[m];
            int A = bvir_qt[a];
            for(int ni = 0; ni < G.params->coltot[h]; ++ni){
                int n = G.params->colorb[h][ni][0];
                int i = G.params->colorb[h][ni][1];
                int N = aocc_qt[n];
                int I = bocc_qt[i];
                double value = G.matrix[h][ma][ni];
                iwl_buf_wrt_val(&AB, M, N, A, I, value, 0, "NULL", 0);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);
   }// if (orb_opt_ == "FALSE")

    delete [] aocc_qt;
    delete [] bocc_qt;
    delete [] avir_qt;
    delete [] bvir_qt;
    delete [] alpha_pitzer_to_corr;
    delete [] beta_pitzer_to_corr;

    psio_->close(PSIF_OCC_DENSITY, 1);

    iwl_buf_flush(&AA, 1);
    iwl_buf_flush(&AB, 1);
    iwl_buf_flush(&BB, 1);
    iwl_buf_close(&AA, 1);
    iwl_buf_close(&AB, 1);
    iwl_buf_close(&BB, 1);

}// end if (reference_ == "UNRESTRICTED")
 //outfile->Printf("\n dump_pdms done. \n");
}// end of dump_pdms

//========================================================================
//         Form Effective PDMs
//========================================================================
void OCCWave::effective_pdms()
{
 if (reference_ == "RESTRICTED") {
     // VO and OV Blocks
     for(int x = 0; x < nidpA; x++) {
	 int a = idprowA[x];
	 int i = idpcolA[x];
	 int h = idpirrA[x];
	 g1symm->add(h, a + occpiA[h], i, 2.0 * zvectorA->get(x));
	 g1symm->add(h, i, a + occpiA[h], 2.0 * zvectorA->get(x));
     }

    // Form VOOO-block TPDM
    dpdbuf4 G, G2;
    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);
    // G_amin = 8 z_ai delta_mn
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "TPDM <VO|OO>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        #pragma omp parallel for
        for(int am = 0; am < G.params->rowtot[h]; ++am){
            int a = G.params->roworb[h][am][0];
            int m = G.params->roworb[h][am][1];
	    int ha = G.params->psym[a];
            for(int in = 0; in < G.params->coltot[h]; ++in){
                int i = G.params->colorb[h][in][0];
                int n = G.params->colorb[h][in][1];
		int hi = G.params->rsym[i];
		int aa = a - G.params->poff[ha] + occpiA[ha];
		int ii = i - G.params->roff[hi];
		if (m == n && ha == hi) G.matrix[h][am][in] = 8.0*ZmatA->get(ha,aa,ii);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // G_amni += -4 z_ai delta_mn
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "TPDM <VO|OO>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int am = 0; am < G.params->rowtot[h]; ++am){
            int a = G.params->roworb[h][am][0];
            int m = G.params->roworb[h][am][1];
	    int ha = G.params->psym[a];
            for(int ni = 0; ni < G.params->coltot[h]; ++ni){
                int n = G.params->colorb[h][ni][0];
                int i = G.params->colorb[h][ni][1];
		int hi = G.params->ssym[i];
		int aa = a - G.params->poff[ha] + occpiA[ha];
		int ii = i - G.params->soff[hi];
		if (m == n && ha == hi) G.matrix[h][am][ni] -= 4.0*ZmatA->get(ha,aa,ii);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);
    psio_->close(PSIF_OCC_DENSITY, 1);

 }// if (reference_ == "RESTRICTED") {

 else if (reference_ == "UNRESTRICTED") {
     // VO and OV Blocks
     for(int x = 0; x < nidpA; x++) {
	 int a = idprowA[x];
	 int i = idpcolA[x];
	 int h = idpirrA[x];
	 g1symmA->add(h, a + occpiA[h], i, zvectorA->get(x));
	 g1symmA->add(h, i, a + occpiA[h], zvectorA->get(x));
     }

     // vo and ov Blocks
     for(int x = 0; x < nidpB; x++) {
	 int a = idprowB[x];
	 int i = idpcolB[x];
	 int h = idpirrB[x];
	 g1symmB->add(h, a + occpiB[h], i, zvectorB->get(x));
	 g1symmB->add(h, i, a + occpiB[h], zvectorB->get(x));
     }

    // Form VOOO-block TPDM
    dpdbuf4 G, G2;
    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

    // G_AMIN = 2 * Z_AI delta_MN
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "TPDM <VO|OO>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        #pragma omp parallel for
        for(int am = 0; am < G.params->rowtot[h]; ++am){
            int a = G.params->roworb[h][am][0];
            int m = G.params->roworb[h][am][1];
	    int ha = G.params->psym[a];
            for(int in = 0; in < G.params->coltot[h]; ++in){
                int i = G.params->colorb[h][in][0];
                int n = G.params->colorb[h][in][1];
		int hi = G.params->rsym[i];
		int aa = a - G.params->poff[ha] + occpiA[ha];
		int ii = i - G.params->roff[hi];
		if (m == n && ha == hi) G.matrix[h][am][in] = 2.0*ZmatA->get(ha,aa,ii);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // G_AMNI += -2*Z_AI delta_MN
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,O]"), ID("[O,O]"),
                  ID("[V,O]"), ID("[O,O]"), 0, "TPDM <VO|OO>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int am = 0; am < G.params->rowtot[h]; ++am){
            int a = G.params->roworb[h][am][0];
            int m = G.params->roworb[h][am][1];
	    int ha = G.params->psym[a];
            for(int ni = 0; ni < G.params->coltot[h]; ++ni){
                int n = G.params->colorb[h][ni][0];
                int i = G.params->colorb[h][ni][1];
		int hi = G.params->ssym[i];
		int aa = a - G.params->poff[ha] + occpiA[ha];
		int ii = i - G.params->soff[hi];
		if (m == n && ha == hi) G.matrix[h][am][ni] -= 2.0*ZmatA->get(ha,aa,ii);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // G_amin = 2 * Z_ai delta_mn
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "TPDM <vo|oo>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        #pragma omp parallel for
        for(int am = 0; am < G.params->rowtot[h]; ++am){
            int a = G.params->roworb[h][am][0];
            int m = G.params->roworb[h][am][1];
	    int ha = G.params->psym[a];
            for(int in = 0; in < G.params->coltot[h]; ++in){
                int i = G.params->colorb[h][in][0];
                int n = G.params->colorb[h][in][1];
		int hi = G.params->rsym[i];
		int aa = a - G.params->poff[ha] + occpiB[ha];
		int ii = i - G.params->roff[hi];
		if (m == n && ha == hi) G.matrix[h][am][in] = 2.0*ZmatB->get(ha,aa,ii);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // G_amni += -2*Z_ai delta_mn
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,o]"), ID("[o,o]"),
                  ID("[v,o]"), ID("[o,o]"), 0, "TPDM <vo|oo>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int am = 0; am < G.params->rowtot[h]; ++am){
            int a = G.params->roworb[h][am][0];
            int m = G.params->roworb[h][am][1];
	    int ha = G.params->psym[a];
            for(int ni = 0; ni < G.params->coltot[h]; ++ni){
                int n = G.params->colorb[h][ni][0];
                int i = G.params->colorb[h][ni][1];
		int hi = G.params->ssym[i];
		int aa = a - G.params->poff[ha] + occpiB[ha];
		int ii = i - G.params->soff[hi];
		if (m == n && ha == hi) G.matrix[h][am][ni] -= 2.0*ZmatB->get(ha,aa,ii);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // G_AmIn = 2 * Z_AI delta_mn
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[O,o]"),
                  ID("[V,o]"), ID("[O,o]"), 0, "TPDM <Vo|Oo>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        #pragma omp parallel for
        for(int am = 0; am < G.params->rowtot[h]; ++am){
            int a = G.params->roworb[h][am][0];
            int m = G.params->roworb[h][am][1];
	    int ha = G.params->psym[a];
            for(int in = 0; in < G.params->coltot[h]; ++in){
                int i = G.params->colorb[h][in][0];
                int n = G.params->colorb[h][in][1];
		int hi = G.params->rsym[i];
		int aa = a - G.params->poff[ha] + occpiA[ha];
		int ii = i - G.params->roff[hi];
		if (m == n && ha == hi) G.matrix[h][am][in] = 2.0*ZmatA->get(ha,aa,ii);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    // G_MaNi += 2 * Z_ai delta_MN
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,o]"),
                  ID("[O,v]"), ID("[O,o]"), 0, "TPDM <Ov|Oo>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        #pragma omp parallel for
        for(int ma = 0; ma < G.params->rowtot[h]; ++ma){
            int m = G.params->roworb[h][ma][0];
            int a = G.params->roworb[h][ma][1];
	    int ha = G.params->qsym[a];
            for(int ni = 0; ni < G.params->coltot[h]; ++ni){
                int n = G.params->colorb[h][ni][0];
                int i = G.params->colorb[h][ni][1];
		int hi = G.params->ssym[i];
		int aa = a - G.params->qoff[ha] + occpiB[ha];
		int ii = i - G.params->soff[hi];
		if (m == n && ha == hi) G.matrix[h][ma][ni] = 2.0*ZmatB->get(ha,aa,ii);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);
    psio_->close(PSIF_OCC_DENSITY, 1);

 }// else if (reference_ == "UNRESTRICTED")

}// end of effective_pdms

//========================================================================
//         Form Effective Generalized-Fock Matrix
//========================================================================
void OCCWave::effective_gfock()
{
 if (reference_ == "RESTRICTED") {
        // F_ia += 2 f_ii * z_ai
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
          for(int i = 0 ; i < occpiA[h]; ++i){
	    for(int a = 0 ; a < virtpiA[h]; ++a){
                int aa = a + occpiA[h];
                GFock->add(h, i, aa, 2.0 * FockA->get(h, i, i) * ZmatA->get(h, aa, i));
            }
	  }
	}

        // F_ai += 2 f_aa * z_ai
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiA[h]; ++a){
              int aa = a + occpiA[h];
            for(int i = 0 ; i < occpiA[h]; ++i){
                GFock->add(h, aa, i, 2.0 * FockA->get(h, aa, aa) * ZmatA->get(h, aa, i));
            }
	  }
	}

    // OPEN DPD FILES
    dpdbuf4 G, K;
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // F_ai += 8 \sum_{e,m} Z_em <mi|ea>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int mi = 0; mi < K.params->rowtot[h]; ++mi){
            int m = K.params->roworb[h][mi][0];
            int i = K.params->roworb[h][mi][1];
	    int hm = K.params->psym[m];
	    int hi = K.params->qsym[i];
	    int mm = m - K.params->poff[hm];
	    int ii = i - K.params->qoff[hi];
            for(int ea = 0; ea < K.params->coltot[h]; ++ea){
                int e = K.params->colorb[h][ea][0];
                int a = K.params->colorb[h][ea][1];
	        int he = K.params->rsym[e];
	        int ha = K.params->ssym[a];
		int ee = e - K.params->roff[he] + occpiA[he];
		int aa = a - K.params->soff[ha] + occpiA[ha];
                if (he == hm && ha == hi) {
                    double value = 8.0 * ZmatA->get(he, ee, mm) * K.matrix[h][mi][ea];
                    GFock->add(ha, aa, ii, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_ai += -2 \sum_{e,m} Z_em <im|ea>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int im = 0; im < K.params->rowtot[h]; ++im){
            int i = K.params->roworb[h][im][0];
            int m = K.params->roworb[h][im][1];
	    int hi = K.params->psym[i];
	    int hm = K.params->qsym[m];
	    int ii = i - K.params->poff[hi];
	    int mm = m - K.params->qoff[hm];
            for(int ea = 0; ea < K.params->coltot[h]; ++ea){
                int e = K.params->colorb[h][ea][0];
                int a = K.params->colorb[h][ea][1];
	        int he = K.params->rsym[e];
	        int ha = K.params->ssym[a];
		int ee = e - K.params->roff[he] + occpiA[he];
		int aa = a - K.params->soff[ha] + occpiA[ha];
                if (he == hm && ha == hi) {
                    double value = -2.0 * ZmatA->get(he, ee, mm) * K.matrix[h][im][ea];
                    GFock->add(ha, aa, ii, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_ai += -2 \sum_{e,m} Z_em <ie|ma>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int ie = 0; ie < K.params->rowtot[h]; ++ie){
            int i = K.params->roworb[h][ie][0];
            int e = K.params->roworb[h][ie][1];
	    int hi = K.params->psym[i];
	    int he = K.params->qsym[e];
	    int ii = i - K.params->poff[hi];
            int ee = e - K.params->qoff[he] + occpiA[he];
            for(int ma = 0; ma < K.params->coltot[h]; ++ma){
                int m = K.params->colorb[h][ma][0];
                int a = K.params->colorb[h][ma][1];
	        int hm = K.params->rsym[m];
	        int ha = K.params->ssym[a];
	        int mm = m - K.params->roff[hm];
		int aa = a - K.params->soff[ha] + occpiA[ha];
                if (he == hm && ha == hi) {
                    double value = -2.0 * ZmatA->get(he, ee, mm) * K.matrix[h][ie][ma];
                    GFock->add(ha, aa, ii, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_ij += 8 \sum_{e,m} Z_em <jm|ie>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO|OV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int jm = 0; jm < K.params->rowtot[h]; ++jm){
            int j = K.params->roworb[h][jm][0];
            int m = K.params->roworb[h][jm][1];
	    int hj = K.params->psym[j];
	    int hm = K.params->qsym[m];
	    int jj = j - K.params->poff[hj];
	    int mm = m - K.params->qoff[hm];
            for(int ie = 0; ie < K.params->coltot[h]; ++ie){
                int i = K.params->colorb[h][ie][0];
                int e = K.params->colorb[h][ie][1];
	        int hi = K.params->rsym[i];
	        int he = K.params->ssym[e];
	        int ii = i - K.params->roff[hi];
		int ee = e - K.params->soff[he] + occpiA[he];
                if (he == hm && hi == hj) {
                    double value = 8.0 * ZmatA->get(he, ee, mm) * K.matrix[h][jm][ie];
                    GFock->add(hi, ii, jj, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_ij += -2 \sum_{e,m} Z_em <mj|ie>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO|OV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int mj = 0; mj < K.params->rowtot[h]; ++mj){
            int m = K.params->roworb[h][mj][0];
            int j = K.params->roworb[h][mj][1];
	    int hm = K.params->psym[m];
	    int hj = K.params->qsym[j];
	    int mm = m - K.params->poff[hm];
	    int jj = j - K.params->qoff[hj];
            for(int ie = 0; ie < K.params->coltot[h]; ++ie){
                int i = K.params->colorb[h][ie][0];
                int e = K.params->colorb[h][ie][1];
	        int hi = K.params->rsym[i];
	        int he = K.params->ssym[e];
	        int ii = i - K.params->roff[hi];
		int ee = e - K.params->soff[he] + occpiA[he];
                if (he == hm && hi == hj) {
                    double value = -2.0 * ZmatA->get(he, ee, mm) * K.matrix[h][mj][ie];
                    GFock->add(hi, ii, jj, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_ij += -2 \sum_{e,m} Z_em <mi|je>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO|OV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int mi = 0; mi < K.params->rowtot[h]; ++mi){
            int m = K.params->roworb[h][mi][0];
            int i = K.params->roworb[h][mi][1];
	    int hm = K.params->psym[m];
	    int hi = K.params->qsym[i];
	    int mm = m - K.params->poff[hm];
	    int ii = i - K.params->qoff[hi];
            for(int je = 0; je < K.params->coltot[h]; ++je){
                int j = K.params->colorb[h][je][0];
                int e = K.params->colorb[h][je][1];
	        int hj = K.params->rsym[j];
	        int he = K.params->ssym[e];
	        int jj = j - K.params->roff[hj];
		int ee = e - K.params->soff[he] + occpiA[he];
                if (he == hm && hi == hj) {
                    double value = -2.0 * ZmatA->get(he, ee, mm) * K.matrix[h][mi][je];
                    GFock->add(hi, ii, jj, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // close
    psio_->close(PSIF_LIBTRANS_DPD, 1);

    // clean up
    delete [] idprowA;
    delete [] idpcolA;
    delete [] idpirrA;
    delete zvectorA;
 }// if (reference_ == "RESTRICTED")


 else if (reference_ == "UNRESTRICTED") {
        // F_IA += f_II * z_AI
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
          for(int i = 0 ; i < occpiA[h]; ++i){
	    for(int a = 0 ; a < virtpiA[h]; ++a){
                int aa = a + occpiA[h];
                GFockA->add(h, i, aa, FockA->get(h, i, i) * ZmatA->get(h, aa, i));
            }
	  }
	}

        // F_ia += f_ii * z_ai
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
          for(int i = 0 ; i < occpiB[h]; ++i){
	    for(int a = 0 ; a < virtpiB[h]; ++a){
                int aa = a + occpiB[h];
                GFockB->add(h, i, aa, FockB->get(h, i, i) * ZmatB->get(h, aa, i));
            }
	  }
	}

        // F_AI += f_AA * z_AI
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiA[h]; ++a){
              int aa = a + occpiA[h];
            for(int i = 0 ; i < occpiA[h]; ++i){
                GFockA->add(h, aa, i, FockA->get(h, aa, aa) * ZmatA->get(h, aa, i));
            }
	  }
	}

        // F_ai += f_aa * z_ai
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiB[h]; ++a){
              int aa = a + occpiB[h];
            for(int i = 0 ; i < occpiB[h]; ++i){
                GFockB->add(h, aa, i, FockB->get(h, aa, aa) * ZmatB->get(h, aa, i));
            }
	  }
	}

    // OPEN DPD FILES
    dpdbuf4 G, K;
    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    // F_AI += 2 \sum_{E,M} Z_EM <MI|EA>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int mi = 0; mi < K.params->rowtot[h]; ++mi){
            int m = K.params->roworb[h][mi][0];
            int i = K.params->roworb[h][mi][1];
	    int hm = K.params->psym[m];
	    int hi = K.params->qsym[i];
	    int mm = m - K.params->poff[hm];
	    int ii = i - K.params->qoff[hi];
            for(int ea = 0; ea < K.params->coltot[h]; ++ea){
                int e = K.params->colorb[h][ea][0];
                int a = K.params->colorb[h][ea][1];
	        int he = K.params->rsym[e];
	        int ha = K.params->ssym[a];
		int ee = e - K.params->roff[he] + occpiA[he];
		int aa = a - K.params->soff[ha] + occpiA[ha];
                if (he == hm && ha == hi) {
                    double value = 2.0 * ZmatA->get(he, ee, mm) * K.matrix[h][mi][ea];
                    GFockA->add(ha, aa, ii, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_AI += -\sum_{E,M} Z_EM <IM|EA>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int im = 0; im < K.params->rowtot[h]; ++im){
            int i = K.params->roworb[h][im][0];
            int m = K.params->roworb[h][im][1];
	    int hi = K.params->psym[i];
	    int hm = K.params->qsym[m];
	    int ii = i - K.params->poff[hi];
	    int mm = m - K.params->qoff[hm];
            for(int ea = 0; ea < K.params->coltot[h]; ++ea){
                int e = K.params->colorb[h][ea][0];
                int a = K.params->colorb[h][ea][1];
	        int he = K.params->rsym[e];
	        int ha = K.params->ssym[a];
		int ee = e - K.params->roff[he] + occpiA[he];
		int aa = a - K.params->soff[ha] + occpiA[ha];
                if (he == hm && ha == hi) {
                    double value = -ZmatA->get(he, ee, mm) * K.matrix[h][im][ea];
                    GFockA->add(ha, aa, ii, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_AI += -\sum_{E,M} Z_EM <IE|MA>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV|OV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int ie = 0; ie < K.params->rowtot[h]; ++ie){
            int i = K.params->roworb[h][ie][0];
            int e = K.params->roworb[h][ie][1];
	    int hi = K.params->psym[i];
	    int he = K.params->qsym[e];
	    int ii = i - K.params->poff[hi];
            int ee = e - K.params->qoff[he] + occpiA[he];
            for(int ma = 0; ma < K.params->coltot[h]; ++ma){
                int m = K.params->colorb[h][ma][0];
                int a = K.params->colorb[h][ma][1];
	        int hm = K.params->rsym[m];
	        int ha = K.params->ssym[a];
	        int mm = m - K.params->roff[hm];
		int aa = a - K.params->soff[ha] + occpiA[ha];
                if (he == hm && ha == hi) {
                    double value = -ZmatA->get(he, ee, mm) * K.matrix[h][ie][ma];
                    GFockA->add(ha, aa, ii, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_AI += 2 \sum_{e,m} Z_em <Im|Ae>
    //outfile->Printf( "\tI am here\n");
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int im = 0; im < K.params->rowtot[h]; ++im){
            int i = K.params->roworb[h][im][0];
            int m = K.params->roworb[h][im][1];
	    int hi = K.params->psym[i];
	    int hm = K.params->qsym[m];
	    int ii = i - K.params->poff[hi];
	    int mm = m - K.params->qoff[hm];
            for(int ae = 0; ae < K.params->coltot[h]; ++ae){
                int a = K.params->colorb[h][ae][0];
                int e = K.params->colorb[h][ae][1];
	        int ha = K.params->rsym[a];
	        int he = K.params->ssym[e];
		int aa = a - K.params->roff[ha] + occpiA[ha];
		int ee = e - K.params->soff[he] + occpiB[he];
                if (he == hm && ha == hi) {
                    double value = 2.0 * ZmatB->get(he, ee, mm) * K.matrix[h][im][ae];
                    GFockA->add(ha, aa, ii, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // Beta
    // F_ai += 2 \sum_{e,m} Z_em <mi|ea>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int mi = 0; mi < K.params->rowtot[h]; ++mi){
            int m = K.params->roworb[h][mi][0];
            int i = K.params->roworb[h][mi][1];
	    int hm = K.params->psym[m];
	    int hi = K.params->qsym[i];
	    int mm = m - K.params->poff[hm];
	    int ii = i - K.params->qoff[hi];
            for(int ea = 0; ea < K.params->coltot[h]; ++ea){
                int e = K.params->colorb[h][ea][0];
                int a = K.params->colorb[h][ea][1];
	        int he = K.params->rsym[e];
	        int ha = K.params->ssym[a];
		int ee = e - K.params->roff[he] + occpiB[he];
		int aa = a - K.params->soff[ha] + occpiB[ha];
                if (he == hm && ha == hi) {
                    double value = 2.0 * ZmatB->get(he, ee, mm) * K.matrix[h][mi][ea];
                    GFockB->add(ha, aa, ii, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_ai += -\sum_{e,m} Z_em <im|ea>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int im = 0; im < K.params->rowtot[h]; ++im){
            int i = K.params->roworb[h][im][0];
            int m = K.params->roworb[h][im][1];
	    int hi = K.params->psym[i];
	    int hm = K.params->qsym[m];
	    int ii = i - K.params->poff[hi];
	    int mm = m - K.params->qoff[hm];
            for(int ea = 0; ea < K.params->coltot[h]; ++ea){
                int e = K.params->colorb[h][ea][0];
                int a = K.params->colorb[h][ea][1];
	        int he = K.params->rsym[e];
	        int ha = K.params->ssym[a];
		int ee = e - K.params->roff[he] + occpiB[he];
		int aa = a - K.params->soff[ha] + occpiB[ha];
                if (he == hm && ha == hi) {
                    double value = -ZmatB->get(he, ee, mm) * K.matrix[h][im][ea];
                    GFockB->add(ha, aa, ii, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_ai += -\sum_{e,m} Z_em <ie|ma>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov|ov>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int ie = 0; ie < K.params->rowtot[h]; ++ie){
            int i = K.params->roworb[h][ie][0];
            int e = K.params->roworb[h][ie][1];
	    int hi = K.params->psym[i];
	    int he = K.params->qsym[e];
	    int ii = i - K.params->poff[hi];
            int ee = e - K.params->qoff[he] + occpiB[he];
            for(int ma = 0; ma < K.params->coltot[h]; ++ma){
                int m = K.params->colorb[h][ma][0];
                int a = K.params->colorb[h][ma][1];
	        int hm = K.params->rsym[m];
	        int ha = K.params->ssym[a];
	        int mm = m - K.params->roff[hm];
		int aa = a - K.params->soff[ha] + occpiB[ha];
                if (he == hm && ha == hi) {
                    double value = -ZmatB->get(he, ee, mm) * K.matrix[h][ie][ma];
                    GFockB->add(ha, aa, ii, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_ai += 2 \sum_{E,M} Z_EM <Mi|Ea>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int mi = 0; mi < K.params->rowtot[h]; ++mi){
            int m = K.params->roworb[h][mi][0];
            int i = K.params->roworb[h][mi][1];
	    int hm = K.params->psym[m];
	    int hi = K.params->qsym[i];
	    int mm = m - K.params->poff[hm];
	    int ii = i - K.params->qoff[hi];
            for(int ea = 0; ea < K.params->coltot[h]; ++ea){
                int e = K.params->colorb[h][ea][0];
                int a = K.params->colorb[h][ea][1];
	        int he = K.params->rsym[e];
	        int ha = K.params->ssym[a];
		int ee = e - K.params->roff[he] + occpiA[he];
		int aa = a - K.params->soff[ha] + occpiB[ha];
                if (he == hm && ha == hi) {
                    double value = 2.0 * ZmatA->get(he, ee, mm) * K.matrix[h][mi][ea];
                    GFockB->add(ha, aa, ii, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_IJ += 2 \sum_{E,M} Z_EM <JM|IE>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO|OV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int jm = 0; jm < K.params->rowtot[h]; ++jm){
            int j = K.params->roworb[h][jm][0];
            int m = K.params->roworb[h][jm][1];
	    int hj = K.params->psym[j];
	    int hm = K.params->qsym[m];
	    int jj = j - K.params->poff[hj];
	    int mm = m - K.params->qoff[hm];
            for(int ie = 0; ie < K.params->coltot[h]; ++ie){
                int i = K.params->colorb[h][ie][0];
                int e = K.params->colorb[h][ie][1];
	        int hi = K.params->rsym[i];
	        int he = K.params->ssym[e];
	        int ii = i - K.params->roff[hi];
		int ee = e - K.params->soff[he] + occpiA[he];
                if (he == hm && hi == hj) {
                    double value = 2.0 * ZmatA->get(he, ee, mm) * K.matrix[h][jm][ie];
                    GFockA->add(hi, ii, jj, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_IJ += -\sum_{EM} Z_EM <MJ|IE>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO|OV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int mj = 0; mj < K.params->rowtot[h]; ++mj){
            int m = K.params->roworb[h][mj][0];
            int j = K.params->roworb[h][mj][1];
	    int hm = K.params->psym[m];
	    int hj = K.params->qsym[j];
	    int mm = m - K.params->poff[hm];
	    int jj = j - K.params->qoff[hj];
            for(int ie = 0; ie < K.params->coltot[h]; ++ie){
                int i = K.params->colorb[h][ie][0];
                int e = K.params->colorb[h][ie][1];
	        int hi = K.params->rsym[i];
	        int he = K.params->ssym[e];
	        int ii = i - K.params->roff[hi];
		int ee = e - K.params->soff[he] + occpiA[he];
                if (he == hm && hi == hj) {
                    double value = -ZmatA->get(he, ee, mm) * K.matrix[h][mj][ie];
                    GFockA->add(hi, ii, jj, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_IJ += -\sum_{EM} Z_EM <MI|JE>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO|OV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int mi = 0; mi < K.params->rowtot[h]; ++mi){
            int m = K.params->roworb[h][mi][0];
            int i = K.params->roworb[h][mi][1];
	    int hm = K.params->psym[m];
	    int hi = K.params->qsym[i];
	    int mm = m - K.params->poff[hm];
	    int ii = i - K.params->qoff[hi];
            for(int je = 0; je < K.params->coltot[h]; ++je){
                int j = K.params->colorb[h][je][0];
                int e = K.params->colorb[h][je][1];
	        int hj = K.params->rsym[j];
	        int he = K.params->ssym[e];
	        int jj = j - K.params->roff[hj];
		int ee = e - K.params->soff[he] + occpiA[he];
                if (he == hm && hi == hj) {
                    double value = -ZmatA->get(he, ee, mm) * K.matrix[h][mi][je];
                    GFockA->add(hi, ii, jj, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_IJ += 2 \sum_{em} Z_em <Jm|Ie>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <Oo|Ov>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int jm = 0; jm < K.params->rowtot[h]; ++jm){
            int j = K.params->roworb[h][jm][0];
            int m = K.params->roworb[h][jm][1];
	    int hj = K.params->psym[j];
	    int hm = K.params->qsym[m];
	    int jj = j - K.params->poff[hj];
	    int mm = m - K.params->qoff[hm];
            for(int ie = 0; ie < K.params->coltot[h]; ++ie){
                int i = K.params->colorb[h][ie][0];
                int e = K.params->colorb[h][ie][1];
	        int hi = K.params->rsym[i];
	        int he = K.params->ssym[e];
	        int ii = i - K.params->roff[hi];
		int ee = e - K.params->soff[he] + occpiB[he];
                if (he == hm && hi == hj) {
                    double value = 2.0 * ZmatB->get(he, ee, mm) * K.matrix[h][jm][ie];
                    GFockA->add(hi, ii, jj, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // Beta
    // F_ij += 2 \sum_{em} Z_em <jm|ie>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo|ov>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int jm = 0; jm < K.params->rowtot[h]; ++jm){
            int j = K.params->roworb[h][jm][0];
            int m = K.params->roworb[h][jm][1];
	    int hj = K.params->psym[j];
	    int hm = K.params->qsym[m];
	    int jj = j - K.params->poff[hj];
	    int mm = m - K.params->qoff[hm];
            for(int ie = 0; ie < K.params->coltot[h]; ++ie){
                int i = K.params->colorb[h][ie][0];
                int e = K.params->colorb[h][ie][1];
	        int hi = K.params->rsym[i];
	        int he = K.params->ssym[e];
	        int ii = i - K.params->roff[hi];
		int ee = e - K.params->soff[he] + occpiB[he];
                if (he == hm && hi == hj) {
                    double value = 2.0 * ZmatB->get(he, ee, mm) * K.matrix[h][jm][ie];
                    GFockB->add(hi, ii, jj, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_ij += -\sum_{em} Z_em <mj|ie>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo|ov>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int mj = 0; mj < K.params->rowtot[h]; ++mj){
            int m = K.params->roworb[h][mj][0];
            int j = K.params->roworb[h][mj][1];
	    int hm = K.params->psym[m];
	    int hj = K.params->qsym[j];
	    int mm = m - K.params->poff[hm];
	    int jj = j - K.params->qoff[hj];
            for(int ie = 0; ie < K.params->coltot[h]; ++ie){
                int i = K.params->colorb[h][ie][0];
                int e = K.params->colorb[h][ie][1];
	        int hi = K.params->rsym[i];
	        int he = K.params->ssym[e];
	        int ii = i - K.params->roff[hi];
		int ee = e - K.params->soff[he] + occpiB[he];
                if (he == hm && hi == hj) {
                    double value = -ZmatB->get(he, ee, mm) * K.matrix[h][mj][ie];
                    GFockB->add(hi, ii, jj, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_ij += -\sum_{em} Z_em <mi|je>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo|ov>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int mi = 0; mi < K.params->rowtot[h]; ++mi){
            int m = K.params->roworb[h][mi][0];
            int i = K.params->roworb[h][mi][1];
	    int hm = K.params->psym[m];
	    int hi = K.params->qsym[i];
	    int mm = m - K.params->poff[hm];
	    int ii = i - K.params->qoff[hi];
            for(int je = 0; je < K.params->coltot[h]; ++je){
                int j = K.params->colorb[h][je][0];
                int e = K.params->colorb[h][je][1];
	        int hj = K.params->rsym[j];
	        int he = K.params->ssym[e];
	        int jj = j - K.params->roff[hj];
		int ee = e - K.params->soff[he] + occpiB[he];
                if (he == hm && hi == hj) {
                    double value = -ZmatB->get(he, ee, mm) * K.matrix[h][mi][je];
                    GFockB->add(hi, ii, jj, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // F_ij += 2 \sum_{EM} Z_EM <Mj|Ei>
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "MO Ints <Oo|Vo>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int mj = 0; mj < K.params->rowtot[h]; ++mj){
            int m = K.params->roworb[h][mj][0];
            int j = K.params->roworb[h][mj][1];
	    int hm = K.params->psym[m];
	    int hj = K.params->qsym[j];
	    int mm = m - K.params->poff[hm];
	    int jj = j - K.params->qoff[hj];
            for(int ei = 0; ei < K.params->coltot[h]; ++ei){
                int e = K.params->colorb[h][ei][0];
                int i = K.params->colorb[h][ei][1];
	        int he = K.params->rsym[e];
	        int hi = K.params->ssym[i];
		int ee = e - K.params->roff[he] + occpiA[he];
	        int ii = i - K.params->soff[hi];
                if (he == hm && hi == hj) {
                    double value = 2.0 * ZmatA->get(he, ee, mm) * K.matrix[h][mj][ei];
                    GFockB->add(hi, ii, jj, value);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);

    // close
    psio_->close(PSIF_LIBTRANS_DPD, 1);

    // clean up
    delete [] idprowA;
    delete [] idpcolA;
    delete [] idpirrA;
    delete [] idprowB;
    delete [] idpcolB;
    delete [] idpirrB;
    delete zvectorA;
    delete zvectorB;
 }// else if (reference_ == "UNRESTRICTED")

}// end of effective_gfock

//=========================
// OEPROP
//=========================
void OCCWave::oeprop()
{
    outfile->Printf("\tComputing one-electron properties...\n");


    //SharedMatrix Da_ = SharedMatrix(new Matrix("MO-basis alpha OPDM", nmo_, nmo_));
    //SharedMatrix Db_ = SharedMatrix(new Matrix("MO-basis beta OPDM", nmo_, nmo_));
    SharedMatrix Da_ = SharedMatrix(new Matrix("MO-basis alpha OPDM", nirrep_, nmopi_, nmopi_));
    SharedMatrix Db_ = SharedMatrix(new Matrix("MO-basis beta OPDM", nirrep_, nmopi_, nmopi_));
    if (reference_ == "RESTRICTED") {
        Da_->copy(g1symm);
        Da_->scale(0.5);
        Db_->copy(Da_);
    }

    else if (reference_ == "UNRESTRICTED") {
        Da_->copy(g1symmA);
        Db_->copy(g1symmB);
    }

    // Compute oeprop
    std::shared_ptr<OEProp> oe(new OEProp(shared_from_this()));
    oe->set_Da_mo(Da_);
    if (reference_ == "UNRESTRICTED") oe->set_Db_mo(Db_);
    oe->add("DIPOLE");
    oe->add("QUADRUPOLE");
    oe->add("MULLIKEN_CHARGES");
    oe->add("NO_OCCUPATIONS");
    oe->set_title("DF-OMP2");
    oe->compute();
    Da_.reset();
    Db_.reset();

} // end oeprop


}} // End Namespaces
