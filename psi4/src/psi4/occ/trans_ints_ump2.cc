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

#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/matrix.h"
#include "occwave.h"
#include "defines.h"


using namespace std;

namespace psi{ namespace occwave{

void OCCWave::trans_ints_ump2()
{
    //outfile->Printf("\n trans_ints is starting... \n");
/********************************************************************************************/
/************************** Transform 2-electron int. to MO space ***************************/
/********************************************************************************************/
    ints->update_orbitals();
    ints->set_print(print_ - 2 >= 0 ? print_ - 2 : 0);
    ints->set_keep_dpd_so_ints(1);

    // Trans (OV|OV)
    timer_on("Trans (OV|OV)");
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir);
    timer_off("Trans (OV|OV)");

/********************************************************************************************/
/************************** sort chem -> phys ***********************************************/
/********************************************************************************************/
     dpdbuf4 K, G;

     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

     // Build MO ints
     timer_on("Sort chem -> phys");
     // (OV|OV) -> <OO|VV>
     timer_on("Sort (OV|OV) -> <OO|VV>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
     global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");
     global_dpd_->buf4_close(&K);

     // (ov|ov) -> <oo|vv>
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov|ov)");
     global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[o,o]"), ID("[v,v]"), "MO Ints <oo|vv>");
     global_dpd_->buf4_close(&K);

     // (OV|ov) -> <Oo|Vv>
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
     global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,o]"), ID("[V,v]"), "MO Ints <Oo|Vv>");
     global_dpd_->buf4_close(&K);

     // (OV|ov) -> <Ov|Vo>: <Ia||Bj> = <Ia|Bj> = (IB|ja)
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
     global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , psqr, ID("[O,v]"), ID("[V,o]"), "MO Ints <Ov|Vo>");
     global_dpd_->buf4_close(&K);
     timer_off("Sort (OV|OV) -> <OO|VV>");
     timer_off("Sort chem -> phys");

/********************************************************************************************/
/************************** Antisymmetrized Ints ********************************************/
/********************************************************************************************/
      timer_on("Antisymmetrize integrals");
     // <OO||VV>:  <IJ||AB> = <IJ|AB> - <IJ|BA>
     timer_on("Make <OO||VV>");
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     global_dpd_->buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <OO||VV>");
     global_dpd_->buf4_close(&K);
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
     global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
     for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
	global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            for(int ab = 0; ab < K.params->coltot[h]; ++ab){
                int a = K.params->colorb[h][ab][0];
                int b = K.params->colorb[h][ab][1];
		int ba = G.params->colidx[b][a];
		K.matrix[h][ij][ab] -= G.matrix[h][ij][ba];
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&K, h);
        global_dpd_->buf4_mat_irrep_close(&K, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);


     // <oo||vv>:  <ij||ab> = <ij|ab> - <ij|ba>
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
     global_dpd_->buf4_copy(&K, PSIF_LIBTRANS_DPD, "MO Ints <oo||vv>");
     global_dpd_->buf4_close(&K);
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
     global_dpd_->buf4_init(&G, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo|vv>");
     for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
	global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        for(int ij = 0; ij < K.params->rowtot[h]; ++ij){
            for(int ab = 0; ab < K.params->coltot[h]; ++ab){
                int a = K.params->colorb[h][ab][0];
                int b = K.params->colorb[h][ab][1];
		int ba = G.params->colidx[b][a];
		K.matrix[h][ij][ab] -= G.matrix[h][ij][ba];
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&K, h);
        global_dpd_->buf4_mat_irrep_close(&K, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_close(&G);
    timer_off("Make <OO||VV>");
    timer_off("Antisymmetrize integrals");

/********************************************************************************************/
/************************** Transform 1-electron int. to MO space ***************************/
/********************************************************************************************/
      // Trans H matrix
      timer_on("Trans OEI");
      HmoA->copy(Hso);
      HmoB->copy(Hso);
      HmoA->transform(Ca_);
      HmoB->transform(Cb_);
      timer_off("Trans OEI");

      if (print_ > 1) {
	HmoA->print();
	HmoB->print();
      }

         for(int h = 0; h < nirrep_; ++h){
             for(int i = 0; i < occpiA[h]; ++i) FockA->set(h, i, i, epsilon_a_->get(h,i));
             for(int i = 0; i < occpiB[h]; ++i) FockB->set(h, i, i, epsilon_b_->get(h,i));
             for(int a = 0; a < virtpiA[h]; ++a) FockA->set(h, a + occpiA[h], a + occpiA[h], epsilon_a_->get(h, a + occpiA[h]));
             for(int a = 0; a < virtpiB[h]; ++a) FockB->set(h, a + occpiB[h], a + occpiB[h], epsilon_b_->get(h, a + occpiB[h]));
         }

      timer_on("Build Denominators");
      denominators_ump2();
      timer_off("Build Denominators");
      psio_->close(PSIF_LIBTRANS_DPD, 1);
      //outfile->Printf("\n trans_ints done. \n");

}//


void OCCWave::denominators_ump2()
{
    //outfile->Printf("\n denominators is starting... \n");
    dpdbuf4 D;
    dpdfile2 Fo,Fv;

    double *aOccEvals = new double [nacooA];
    double *bOccEvals = new double [nacooB];
    double *aVirEvals = new double [nacvoA];
    double *bVirEvals = new double [nacvoB];
    memset(aOccEvals, 0, sizeof(double)*nacooA);
    memset(bOccEvals, 0, sizeof(double)*nacooB);
    memset(aVirEvals, 0, sizeof(double)*nacvoA);
    memset(bVirEvals, 0, sizeof(double)*nacvoB);

    // Pick out the diagonal elements of the Fock matrix, making sure that they are in the order
    // used by the DPD library, i.e. starting from zero for each space and ordering by irrep

    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0;

    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < aoccpiA[h]; ++i) aOccEvals[aOccCount++] = epsilon_a_->get(h, i + frzcpi_[h]);
	for(int i = 0; i < aoccpiB[h]; ++i) bOccEvals[bOccCount++] = epsilon_b_->get(h, i + frzcpi_[h]);
        for(int a = 0; a < avirtpiA[h]; ++a) aVirEvals[aVirCount++] = epsilon_a_->get(h, a + occpiA[h]);
	for(int a = 0; a < avirtpiB[h]; ++a) bVirEvals[bVirCount++] = epsilon_b_->get(h, a + occpiB[h]);
    }

    // Build denominators
    // The alpha-alpha spin case
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "D <OO|VV>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(aOccEvals[i] + aOccEvals[j] - aVirEvals[a] - aVirEvals[b]);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, h);
        global_dpd_->buf4_mat_irrep_close(&D, h);
    }
    if (print_ > 2) global_dpd_->buf4_print(&D, "outfile", 1);
    global_dpd_->buf4_close(&D);


    // The beta-beta spin case
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "D <oo|vv>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(bOccEvals[i] + bOccEvals[j] - bVirEvals[a] - bVirEvals[b]);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, h);
        global_dpd_->buf4_mat_irrep_close(&D, h);
    }
    if (print_ > 2) global_dpd_->buf4_print(&D, "outfile", 1);
    global_dpd_->buf4_close(&D);


    // The alpha-beta spin case
    global_dpd_->buf4_init(&D, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "D <Oo|Vv>");
    for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&D, h);
        for(int row = 0; row < D.params->rowtot[h]; ++row){
            int i = D.params->roworb[h][row][0];
            int j = D.params->roworb[h][row][1];
            for(int col = 0; col < D.params->coltot[h]; ++col){
                int a = D.params->colorb[h][col][0];
                int b = D.params->colorb[h][col][1];
                D.matrix[h][row][col] = 1.0/(aOccEvals[i] + bOccEvals[j] - aVirEvals[a] - bVirEvals[b]);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&D, h);
        global_dpd_->buf4_mat_irrep_close(&D, h);
    }
    if (print_ > 2) global_dpd_->buf4_print(&D, "outfile", 1);
    global_dpd_->buf4_close(&D);

    //Print
    if(print_ > 1){
      outfile->Printf("\n \n");
      for(int i = 0; i<nacooA; i++) {
	outfile->Printf("\taOccEvals[%1d]: %20.14f\n", i, aOccEvals[i]);

      }

      outfile->Printf("\n \n");
      for(int i = 0; i<nacooB; i++) {
	outfile->Printf("\tbOccEvals[%1d]: %20.14f\n", i, bOccEvals[i]);

      }

      outfile->Printf("\n \n");
      for(int i = 0; i<nacvoA; i++) {
	outfile->Printf("\taVirEvals[%1d]: %20.14f\n", i, aVirEvals[i]);

      }

      outfile->Printf("\n \n");
      for(int i = 0; i<nacvoB; i++) {
	outfile->Printf("\tbVirEvals[%1d]: %20.14f\n", i, bVirEvals[i]);

      }
    }

    delete [] aOccEvals;
    delete [] bOccEvals;
    delete [] aVirEvals;
    delete [] bVirEvals;

//outfile->Printf("\n denominators done. \n");
}// end denominators
}} // End Namespaces
