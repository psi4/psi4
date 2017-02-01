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
#include "psi4/libpsio/psio.hpp"
#include "occwave.h"
#include "defines.h"


using namespace std;


namespace psi{ namespace occwave{

void OCCWave::cepa_iterations()
{

outfile->Printf("\n  \n");
outfile->Printf(" ============================================================================== \n");
outfile->Printf(" ================ Performing CEPA iterations... =============================== \n");
outfile->Printf(" ============================================================================== \n");
outfile->Printf("\n");
outfile->Printf( "  Iter    E_corr           E_total            DE           T2 RMS        \n");
outfile->Printf( "  ----   -------------    ---------------    ----------   ----------    \n");



/********************************************************************************************/
/************************** NR iterations **************************************************/
/********************************************************************************************/
      itr_occ = 0;
      conver = 1; // Assuming that the iterations will converge
 // DIIS
      if (nooA + nooB != 1) {
          if (reference_ == "RESTRICTED") {
              dpdbuf4 T;
              psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
              global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                     ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
              t2DiisManager = new DIISManager(cc_maxdiis_, "CEPA DIIS T2 Amps", DIISManager::LargestError, DIISManager::OnDisk);
              t2DiisManager->set_error_vector_size(1, DIISEntry::DPDBuf4, &T);
              t2DiisManager->set_vector_size(1, DIISEntry::DPDBuf4, &T);
              global_dpd_->buf4_close(&T);
              psio_->close(PSIF_OCC_DPD, 1);
          }

          else if (reference_ == "UNRESTRICTED") {
              dpdbuf4 Taa, Tbb, Tab;
              psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
              global_dpd_->buf4_init(&Taa, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                                     ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
              global_dpd_->buf4_init(&Tbb, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                                     ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
              global_dpd_->buf4_init(&Tab, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                                     ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
              t2DiisManager = new DIISManager(cc_maxdiis_, "CEPA DIIS T2 Amps", DIISManager::LargestError, DIISManager::InCore);
              t2DiisManager->set_error_vector_size(3, DIISEntry::DPDBuf4, &Taa,
                                                   DIISEntry::DPDBuf4, &Tbb,
                                                   DIISEntry::DPDBuf4, &Tab);
              t2DiisManager->set_vector_size(3, DIISEntry::DPDBuf4, &Taa,
                                             DIISEntry::DPDBuf4, &Tbb,
                                             DIISEntry::DPDBuf4, &Tab);
              global_dpd_->buf4_close(&Taa);
              global_dpd_->buf4_close(&Tbb);
              global_dpd_->buf4_close(&Tab);
              psio_->close(PSIF_OCC_DPD, 1);
          }
      }

// head of loop
do
{
        itr_occ++;
        timer_on("T2");
	t2_amps();
        timer_off("T2");
        timer_on("CEPA Energy");
        cepa_energy();
        timer_off("CEPA Energy");
        Ecorr = Ecepa - Escf;
        DE = Ecepa - Ecepa_old;
        Ecepa_old = Ecepa;

    if (reference_ == "UNRESTRICTED") {
	rms_t2=MAX0(rms_t2AA,rms_t2BB);
	rms_t2=MAX0(rms_t2,rms_t2AB);
    }

outfile->Printf(" %3d     %12.10f    %12.10f  %12.2e %12.2e \n", itr_occ, Ecorr, Ecepa, DE, rms_t2);


    if (itr_occ >= cc_maxiter) {
      conver = 0; // means iterations were NOT converged
      break;
    }

    if (rms_t2 >= DIVERGE) {
        throw PSIEXCEPTION("CEPA iterations are diverging");
    }

}
while(fabs(DE) >= tol_Eod || rms_t2 >= tol_t2);

//delete
delete t2DiisManager;

if (conver == 1) {
EcepaL = Ecepa;
outfile->Printf("\n");
outfile->Printf(" ============================================================================== \n");
outfile->Printf(" ======================== CEPA ITERATIONS ARE CONVERGED ======================= \n");
outfile->Printf(" ============================================================================== \n");

}

else if (conver == 0) {
  outfile->Printf("\n ======================= CEPA IS NOT CONVERGED IN %2d ITERATIONS ============ \n", cc_maxiter);

  throw PSIEXCEPTION("CEPA iterations did not converge");
}

}// end main
}} // End Namespaces
