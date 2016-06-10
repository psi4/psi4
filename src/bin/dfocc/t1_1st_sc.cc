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

/** Standard library includes */
#include <libqt/qt.h>
#include "defines.h"
#include "dfocc.h"

using namespace boost;
using namespace psi;
using namespace std;


namespace psi{ namespace dfoccwave{
  
void DFOCC::t1_1st_sc()
{   

    timer_on("1st-order T1");
    // T1A
    for(int i = 0 ; i < naoccA; ++i){
        for(int a = 0 ; a < navirA; ++a){
            double value = FockA->get(i + nfrzc, i + nfrzc) - FockA->get(a + noccA, a + noccA);
            t1A->set(i, a, FockA->get(i + nfrzc, a + noccA) / value);
        }
    }
    if (print_ > 2) t1A->print();

    // T1B
    for(int i = 0 ; i < naoccB; ++i){
        for(int a = 0 ; a < navirB; ++a){
            double value = FockB->get(i + nfrzc, i + nfrzc) - FockB->get(a + noccB, a + noccB);
            t1B->set(i, a, FockB->get(i + nfrzc, a + noccB) / value);
        }
    }

        //Singles-contribution
        Emp2_t1 = 0.0;
        //Alpha
        for(int i = 0 ; i < naoccA; ++i){
            for(int a = 0 ; a < navirA; ++a){
                Emp2_t1 += t1A->get(i, a) * FockA->get(a + noccA, i + nfrzc);
            }
        }

        // Beta
        for(int i = 0 ; i < naoccB; ++i){
            for(int a = 0 ; a < navirB; ++a){
                Emp2_t1 += t1B->get(i, a) * FockB->get(a + noccB, i + nfrzc);
            }
        }

    if (print_ > 2) t1B->print();
    timer_off("1st-order T1");
} // end t1_1st_sc

}} // End Namespaces
