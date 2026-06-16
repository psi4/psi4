/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libiwl/iwl.hpp"
#include "psi4/libiwl/iwl_reader.h"
#include "psi4/libiwl/iwl_writer.h"

#include "occwave.h"
#include "defines.h"

using namespace psi;

namespace psi {
namespace occwave {

void OCCWave::tei_sort_iabc() {
    // outfile->Printf("\n tei_sort_iabc is starting... \n");
    /********************************************************************************************/
    /************************** sort chem -> phys ***********************************************/
    /********************************************************************************************/
    IWLWriter out(psio_, PSIF_OCC_IABC, cutoff);

    if (print_ > 2) outfile->Printf("\n writing <ia|bc>... \n");
    IWLReader eri(psio_, PSIF_MO_TEI);
    for (const auto &integral : eri) {
        int i = std::abs(integral.p);
        int j = integral.q;
        int k = integral.r;
        int l = integral.s;
        double value = integral.value;

        // Make sure we are dealing with the (ia|bc) type integrals
        if (i < nooA && j >= nooA && k >= nooA && l >= nooA) {
            out.write(i, k, j, l, value);
            if (k > l) out.write(i, l, j, k, value);
        }
    }
    // `out` flushes the last buffer and closes (keeping the file) on scope exit.

    /*
            // Test reading
            dpdbuf4 K;
            //psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
            dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,V]"), ints->DPD_ID("[V,V]"),
                      ints->DPD_ID("[O,V]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <OV|VV>");
            dpd_buf4_print(&K, outfile, 1);
            dpd_buf4_close(&K);
            //psio_->close(PSIF_LIBTRANS_DPD, 1);

            outfile->Printf("\n reading <ia|bc>... \n");
            IWL ERIIN2(psio_.get(), PSIF_OCC_IABC, 0.0, 1, 1);
     do
     {
            ilsti = ERIIN2.last_buffer();
            nbuf = ERIIN2.buffer_count();

       fi = 0;
       for (int idx=0; idx < nbuf; idx++ )
       {

            int i = ERIIN2.labels()[fi];
                i = std::abs(i);
            int j = ERIIN2.labels()[fi+1];
            int k = ERIIN2.labels()[fi+2];
            int l = ERIIN2.labels()[fi+3];
            value = ERIIN2.values()[idx];
            fi += 4;

            outfile->Printf("%3d %3d %3d %3d %20.14f\n",i,j,k,l,value);


       }
            if(!ilsti)
              ERIIN2.fetch();

     } while(!ilsti);
     */
    // outfile->Printf("tei_sort_iabc done. \n");
}  // end sort_iabc
}
}  // End Namespaces
