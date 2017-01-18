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
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsio/psio.hpp"
#include "occwave.h"
#include "defines.h"


using namespace std;

namespace psi{ namespace occwave{

void OCCWave::trans_ints_rmp2()
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
     timer_on("Sort (OV|OV) -> <OO|VV>");
     // (OV|OV) -> <OO|VV>
     global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
     global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints <OO|VV>");
     global_dpd_->buf4_close(&K);
     timer_off("Sort (OV|OV) -> <OO|VV>");
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

         for(int h = 0; h < nirrep_; ++h){
             for(int i = 0; i < occpiA[h]; ++i) FockA->set(h, i, i, epsilon_a_->get(h,i));
             for(int a = 0; a < virtpiA[h]; ++a) FockA->set(h, a + occpiA[h], a + occpiA[h], epsilon_a_->get(h, a + occpiA[h]));
         }

      timer_on("Build Denominators");
      denominators_rmp2();
      timer_off("Build Denominators");
      psio_->close(PSIF_LIBTRANS_DPD, 1);
      //outfile->Printf("\n trans_ints done. \n");

}//


void OCCWave::denominators_rmp2()
{
    //outfile->Printf("\n denominators is starting... \n");
    dpdbuf4 D;
    dpdfile2 Fo,Fv;

    double *aOccEvals = new double [nacooA];
    double *aVirEvals = new double [nacvoA];
    memset(aOccEvals, 0, sizeof(double)*nacooA);
    memset(aVirEvals, 0, sizeof(double)*nacvoA);

    // Pick out the diagonal elements of the Fock matrix, making sure that they are in the order
    // used by the DPD library, i.e. starting from zero for each space and ordering by irrep

    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0;

    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0; i < aoccpiA[h]; ++i) aOccEvals[aOccCount++] = epsilon_a_->get(h, i + frzcpi_[h]);
        for(int a = 0; a < avirtpiA[h]; ++a) aVirEvals[aVirCount++] = epsilon_a_->get(h, a + occpiA[h]);
    }

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

    delete [] aOccEvals;
    delete [] aVirEvals;

//outfile->Printf("\n denominators done. \n");
}// end denominators
}} // End Namespaces
