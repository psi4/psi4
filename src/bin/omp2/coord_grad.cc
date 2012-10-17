/** Standard library includes */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector> 
 
/** Required PSI4 includes */
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libtrans/mospace.h>
#include <libtrans/integraltransform.h>


/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>

#include "omp2wave.h"
#include "defines.h"
#include "arrays.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace omp2wave{
  
void OMP2Wave::coord_grad()
{
      GFockmo_diag();
      dump_ints();
      dump_pdms();

}// 

void OMP2Wave::dump_ints()
{
    //fprintf(outfile,"\n dump_ints is starting... \n"); fflush(outfile);
    dpdfile2 H;
    dpdbuf4 K;

    psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "H <O|O>");
    dpd_file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
                H.matrix[h][i][j] = HmoA->get(h, i, j);
            }
        }
    }
    dpd_file2_mat_wrt(&H);
    dpd_file2_close(&H);

    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "H <V|V>");
    dpd_file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int b = 0 ; b < virtpiA[h]; ++b){
                H.matrix[h][a][b] = HmoA->get(h, a + occpiA[h], b + occpiA[h]);
            }
        }
    }
    dpd_file2_mat_wrt(&H);
    dpd_file2_close(&H);

    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('V'), "H <O|V>");
    dpd_file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
                H.matrix[h][i][j] = HmoA->get(h, i, j + occpiA[h]);
            }
        }
    }
    dpd_file2_mat_wrt(&H);
    dpd_file2_close(&H);

 // Write beta-integrals if ref is UHF
 if (reference == "UHF") {
    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('o'), "H <o|o>");
    dpd_file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < occpiB[h]; ++j){
                H.matrix[h][i][j] = HmoB->get(h, i, j);
            }
        }
    }
    dpd_file2_mat_wrt(&H);
    dpd_file2_close(&H);

    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('v'), ID('v'), "H <v|v>");
    dpd_file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int a = 0 ; a < virtpiB[h]; ++a){
            for(int b = 0 ; b < virtpiB[h]; ++b){
                H.matrix[h][a][b] = HmoB->get(h, a + occpiB[h], b + occpiB[h]);
            }
        }
    }
    dpd_file2_mat_wrt(&H);
    dpd_file2_close(&H);

    dpd_file2_init(&H, PSIF_LIBTRANS_DPD, 0, ID('o'), ID('v'), "H <o|v>");
    dpd_file2_mat_init(&H);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < virtpiB[h]; ++j){
                H.matrix[h][i][j] = HmoB->get(h, i, j + occpiB[h]);
            }
        }
    }
    dpd_file2_mat_wrt(&H);
    dpd_file2_close(&H);
 }// end uhf

    psio_->close(PSIF_LIBTRANS_DPD, 1);
    //fprintf(outfile,"\n dump_ints done. \n"); fflush(outfile);

}// end of dump_ints



void OMP2Wave::dump_pdms()
{
    //fprintf(outfile,"\n dump_pdms is starting... \n"); fflush(outfile);

//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference == "RHF") {
    dpdfile2 H;
    dpdbuf4 G, G2;

    psio_->open(PSIF_OMP2_DENSITY, PSIO_OPEN_OLD);

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
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
              ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
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
                   iwl_buf_wrt_val(&AA, I, K, J, L, value, 0, (FILE *) NULL, 0);
                }
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // OOVV
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
              ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
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
                   iwl_buf_wrt_val(&AA, A, I, B, J, value, 0, (FILE *) NULL, 0);
                }
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // OVOV
    // Alpha-Alpha spin-case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
              ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    dpd_buf4_scm(&G, 2.0);
    dpd_buf4_dump(&G, &AA, aocc_qt, avir_qt, aocc_qt, avir_qt, 0, 1);
    dpd_buf4_close(&G);

    delete [] aocc_qt;
    delete [] avir_qt;

    psio_->close(PSIF_OMP2_DENSITY, 1);

    iwl_buf_flush(&AA, 1);
    iwl_buf_close(&AA, 1);

}// end if (reference == "RHF") 



//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference == "UHF") {

    dpdfile2 H;
    dpdbuf4 G, G2;

    psio_->open(PSIF_OMP2_DENSITY, PSIO_OPEN_OLD);

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
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
              ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
    // Dump tpdm in chemist notation, 0 means no pack, 1 means swap indices 23
    dpd_buf4_dump(&G, &AA, aocc_qt, aocc_qt, aocc_qt, aocc_qt, 0, 1);
    dpd_buf4_close(&G);

    // OOOO: Beta-Beta spin-case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
              ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
    dpd_buf4_dump(&G, &BB, bocc_qt, bocc_qt, bocc_qt, bocc_qt, 0, 1);
    dpd_buf4_close(&G);

    // OOOO: Alpha-Beta spin-case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
              ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
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
                iwl_buf_wrt_val(&AB, I, K, J, L, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);
 
    
    // OOVV: Alpha-Alpha spin-case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
              ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    dpd_buf4_scm(&G, 2.0);
    dpd_buf4_dump(&G, &AA, aocc_qt, aocc_qt, avir_qt, avir_qt, 0, 1);
    dpd_buf4_close(&G);

    // OOVV: Beta-Beta spin-case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
              ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
    dpd_buf4_scm(&G, 2.0);
    dpd_buf4_dump(&G, &BB, bocc_qt, bocc_qt, bvir_qt, bvir_qt, 0, 1);
    dpd_buf4_close(&G);

    // OOVV: Alpha-Beta spin-case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
              ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
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
                iwl_buf_wrt_val(&AB, I, A, J, B, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);


    // OVOV: Alpha-Alpha spin-case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
              ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
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
                iwl_buf_wrt_val(&AA, I, J, A, B, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // OVOV: Beta-Beta spin-case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
              ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
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
                iwl_buf_wrt_val(&BB, I, J, A, B, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    // OVOV: Alpha-Beta spin-case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
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
                iwl_buf_wrt_val(&AB, I, J, A, B, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);



    // OVVO: Alpha-Alpha spin-case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
              ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
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
                iwl_buf_wrt_val(&AA, I, B, A, J, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);


    // OVVO: Beta-Beta spin-case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
              ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
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
                iwl_buf_wrt_val(&BB, I, B, A, J, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);



    // VOVO: Alpha-Beta spin-case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>"); 
    for(int h = 0; h < nirrep_; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        dpd_buf4_mat_irrep_rd(&G, h);
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
                iwl_buf_wrt_val(&AB, A, B, I, J, value, 0, (FILE *) NULL, 0);
            }
        }
        dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    delete [] aocc_qt;
    delete [] bocc_qt;
    delete [] avir_qt;
    delete [] bvir_qt;

    psio_->close(PSIF_OMP2_DENSITY, 1);

    iwl_buf_flush(&AA, 1);
    iwl_buf_flush(&AB, 1);
    iwl_buf_flush(&BB, 1);
    iwl_buf_close(&AA, 1);
    iwl_buf_close(&AB, 1);
    iwl_buf_close(&BB, 1);

}// end if (reference == "UHF") 
 //fprintf(outfile,"\n dump_pdms done. \n"); fflush(outfile);

}// end of dump_pdms
}} // End Namespaces

