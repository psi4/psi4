/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <libqt/qt.h>
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
void DFOCC::ccd_F_intr()
{

    // defs
    SharedTensor2d K, T, U, Tau;

    // read
    Tau = SharedTensor2d(new Tensor2d("T2 (Q|IA)", nQ, naoccA, navirA));
    Tau->read(psio_, PSIF_DFOCC_AMPS);

    // OO block
    // F_mi +=  \sum_{Q,e} Tau"_ie^Q b_me^Q
    FijA->zero();
    FijA->contract332(false, true, navirA, bQiaA, Tau, 1.0, 1.0);

    // VV block
    // F_ae -=  \sum_{Q,m} Tau'_ma^Q b_me^Q
    FabA->contract(true, false, navirA, navirA, nQ * naoccA, Tau, bQiaA, -1.0, 0.0);
    Tau.reset();

    //outfile->Printf("\tF int done.\n");

}// end ccd_F_intr
}} // End Namespaces
