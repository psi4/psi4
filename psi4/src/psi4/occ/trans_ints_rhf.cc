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
#include "psi4/libiwl/iwl.hpp"
#include "psi4/psifiles.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "occwave.h"
#include "defines.h"


using namespace psi;
using namespace std;

namespace psi{ namespace occwave{

void OCCWave::trans_ints_rhf()
{
    //outfile->Printf("\n trans_ints is starting... \n");
/********************************************************************************************/
/************************** Transform 2-electron int. to MO space ***************************/
/********************************************************************************************/
    ints->update_orbitals();
    ints->set_print(print_ - 2 >= 0 ? print_ - 2 : 0);
    ints->set_keep_dpd_so_ints(1);

    // Trans (OO|OO)
    timer_on("Trans (OO|OO)");
    ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ, IntegralTransform::MakeAndKeep);
    timer_off("Trans (OO|OO)");

    // Trans (OO|OV)
    timer_on("Trans (OO|OV)");
    ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::vir, IntegralTransform::ReadAndKeep);
    timer_off("Trans (OO|OV)");

    // Trans (OO|VV)
    timer_on("Trans (OO|VV)");
    ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);
    timer_off("Trans (OO|VV)");

    // Trans (OV|OV)
    timer_on("Trans (OV|OV)");
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir, IntegralTransform::MakeAndKeep);
    timer_off("Trans (OV|OV)");

    // Trans (OV|VV)
    timer_on("Trans (OV|VV)");
    ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);
    timer_off("Trans (OV|VV)");

if (wfn_type_ != "OMP2" || ekt_ea_ == "TRUE") {
    // Trans (VV|VV)
    timer_on("Trans (VV|VV)");
    ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir);
    timer_off("Trans (VV|VV)");
}

/********************************************************************************************/
/************************** sort chem -> phys ***********************************************/
/********************************************************************************************/
     dpdbuf4 K, G;

     psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);

     // Build MO ints
     timer_on("Sort chem -> phys");
     timer_on("Sort (OO|OO) -> <OO|OO>");
     // (OO|OO) -> <OO|OO>
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
     global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[O,O]"), "MO Ints <OO|OO>");
     global_dpd_->buf4_close(&K);
     timer_off("Sort (OO|OO) -> <OO|OO>");


     timer_on("Sort (OO|OV) -> <OO|OV>");
     // (OO|OV) -> <OO|OV>
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O>=O]+"), ID("[O,V]"), 0, "MO Ints (OO|OV)");
     global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[O,V]"), "MO Ints <OO|OV>");
     global_dpd_->buf4_close(&K);
     timer_off("Sort (OO|OV) -> <OO|OV>");


     timer_on("Sort (OV|OV) -> <OO|VV>");
     // (OV|OV) -> <OO|VV>
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
     global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");
     global_dpd_->buf4_close(&K);
     timer_off("Sort (OV|OV) -> <OO|VV>");


     timer_on("Sort (OO|VV) -> <OV|OV>");
     // (OO|VV) -> <OV|OV>
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
     global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "MO Ints <OV|OV>");
     global_dpd_->buf4_close(&K);
     timer_off("Sort (OO|VV) -> <OV|OV>");


     timer_on("Sort (OV|VV) -> <OV|VV>");
     // (OV|VV) -> <OV|VV>
if (wfn_type_ == "OMP2" && incore_iabc_ == 0) {
     tei_sort_iabc();
}

else {
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
     global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,V]"), ID("[V,V]"), "MO Ints <OV|VV>");
     global_dpd_->buf4_close(&K);
}
     timer_off("Sort (OV|VV) -> <OV|VV>");


if (wfn_type_ != "OMP2" || ekt_ea_ == "TRUE") {
     timer_on("Sort (VV|VV) -> <VV|VV>");
     // (VV|VV) -> <VV|VV>
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                 ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
     global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[V,V]"), ID("[V,V]"), "MO Ints <VV|VV>");
     global_dpd_->buf4_close(&K);
     timer_off("Sort (VV|VV) -> <VV|VV>");
}// end if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") {
     timer_off("Sort chem -> phys");

/********************************************************************************************/
/************************** Transform 1-electron int. to MO space ***************************/
/********************************************************************************************/
      // Trans H matrix
      timer_on("Trans OEI");
      HmoA->copy(Hso);
      HmoA->transform(Ca_);
      timer_off("Trans OEI");

      if (print_ > 2) {
	HmoA->print();
      }

      // Trans Fock matrix
      if (orb_opt_ == "TRUE") {
          timer_on("Build Fock");
          fock_alpha();
          timer_off("Build Fock");
      }

      else if (orb_opt_ == "FALSE") {
         for(int h = 0; h < nirrep_; ++h){
             for(int i = 0; i < occpiA[h]; ++i) FockA->set(h, i, i, epsilon_a_->get(h,i));
             for(int a = 0; a < virtpiA[h]; ++a) FockA->set(h, a + occpiA[h], a + occpiA[h], epsilon_a_->get(h, a + occpiA[h]));
         }
      }

      timer_on("Build Denominators");
      if (orb_opt_ == "TRUE") denominators_rhf();
      else if (orb_opt_ == "FALSE") denominators_rmp2();
      timer_off("Build Denominators");
      psio_->close(PSIF_LIBTRANS_DPD, 1);
      //outfile->Printf("\n trans_ints done. \n");

}//


void OCCWave::denominators_rhf()
{
    //outfile->Printf("\n denominators is starting... \n");
    dpdbuf4 D;
    dpdfile2 Fo,Fv;

    double *aOccEvals = new double [nacooA];
    double *aVirEvals = new double [nacvoA];

    // Pick out the diagonal elements of the Fock matrix, making sure that they are in the order
    // used by the DPD library, i.e. starting from zero for each space and ordering by irrep

    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0;

    //Diagonal elements of the Fock matrix
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < aoccpiA[h]; ++i) aOccEvals[aOccCount++] = FockA->get(h, i + frzcpi_[h], i + frzcpi_[h]);
        for(int a = 0; a < avirtpiA[h]; ++a) aVirEvals[aVirCount++] = FockA->get(h, occpiA[h] + a, occpiA[h] + a);
    }

    // Build denominators
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


    //Print
    if(print_ > 1){
      outfile->Printf("\n \n");
      for(int i = 0; i<nacooA; i++) {
	outfile->Printf("\taOccEvals[%1d]: %20.14f\n", i, aOccEvals[i]);

      }

      outfile->Printf("\n \n");
      for(int i = 0; i<nacvoA; i++) {
	outfile->Printf("\taVirEvals[%1d]: %20.14f\n", i, aVirEvals[i]);

      }
    }

    delete [] aOccEvals;
    delete [] aVirEvals;


    // Off-diagonal elements of the Fock matrix
    // Build Occupied-Occupied block
    // The alpha-alpha spin case
    global_dpd_->file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
    global_dpd_->file2_mat_init(&Fo);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
		if (i != j) Fo.matrix[h][i][j] = FockA->get(h, i + frzcpi_[h], j + frzcpi_[h]);
		else Fo.matrix[h][i][j] = 0.0;
            }
        }
    }
    global_dpd_->file2_mat_wrt(&Fo);
    global_dpd_->file2_close(&Fo);

    if (print_ > 2) {
      global_dpd_->file2_init(&Fo, PSIF_LIBTRANS_DPD, 0, ID('O'), ID('O'), "F <O|O>");
      global_dpd_->file2_mat_init(&Fo);
      global_dpd_->file2_mat_print(&Fo, "outfile");
      global_dpd_->file2_close(&Fo);
    }

    // Build Virtual-Virtual block
    // The alpha-alpha spin case
    global_dpd_->file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
    global_dpd_->file2_mat_init(&Fv);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < avirtpiA[h]; ++i){
            for(int j = 0 ; j < avirtpiA[h]; ++j){
                if (i != j) Fv.matrix[h][i][j] = FockA->get(h, i + occpiA[h], j + occpiA[h]);
		else Fv.matrix[h][i][j] = 0.0;
            }
        }
    }
    global_dpd_->file2_mat_wrt(&Fv);
    global_dpd_->file2_close(&Fv);

    if (print_ > 2) {
      global_dpd_->file2_init(&Fv, PSIF_LIBTRANS_DPD, 0, ID('V'), ID('V'), "F <V|V>");
      global_dpd_->file2_mat_init(&Fv);
      global_dpd_->file2_mat_print(&Fv, "outfile");
      global_dpd_->file2_close(&Fv);
    }

//outfile->Printf("\n denominators done. \n");
}// end denominators
}} // End Namespaces
