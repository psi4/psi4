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

/** Standard library includes */
#include "psi4/libqt/qt.h"
#include "defines.h"
#include "dfocc.h"


using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{

void DFOCC::ccd_opdm()
{

    SharedTensor2d T, U, X;
    timer_on("opdm");
//if (reference_ == "RESTRICTED") {

    // G1_ij = -(G_ij + G_ji)
    T = SharedTensor2d(new Tensor2d("G Intermediate <I|J>", naoccA, naoccA));
    T->symmetrize(GijA);
    T->scale(-2.0);
    G1c_oo->set_act_oo(nfrzc, naoccA, T);
    T.reset();

    //  G1_ab = -(G_ab + G_ba)
    T = SharedTensor2d(new Tensor2d("G Intermediate <A|B>", navirA, navirA));
    T->symmetrize(GabA);
    T->scale(-2.0);
    G1c_vv->set_act_vv(T);
    T.reset();

    // set OV block
    G1c_ov->zero();

    // Build G1_ai
    G1c_vo->trans(G1c_ov);

    // Build G1c
    G1c->set_oo(G1c_oo);
    G1c->set_ov(G1c_ov);
    G1c->set_vo(G1c_vo);
    G1c->set_vv(noccA, G1c_vv);
    //G1c->print();

    // Build G1
    G1->copy(G1c);
    for (int i = 0; i < noccA; i++) G1->add(i, i, 2.0);

  if(print_ > 2) {
    G1->print();
    double trace = G1->trace();
    outfile->Printf("\t trace: %12.12f \n", trace);

  }

//}// end if (reference_ == "RESTRICTED")

//else if (reference_ == "UNRESTRICTED") {
//}// else if (reference_ == "UNRESTRICTED")
    timer_off("opdm");
} // end ccd_opdm

}} // End Namespaces
