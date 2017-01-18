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

#include "occwave.h"
#include "defines.h"
#include "psi4/libmints/matrix.h"


using namespace std;

namespace psi{ namespace occwave{


void OCCWave::mograd()
{

      double norm_wogA, norm_wogB;

if (reference_ == "RESTRICTED") {
      // memalloc
      WorbA->zero();

      // set W matrix
      SharedMatrix temp(GFock->transpose());
      WorbA->copy(GFock);
      WorbA->subtract(temp);
      // DO NOT scale W since Fai (spin-free) = 2 * FAI
      //WorbA->scale(2.0);

      // Form w vector
      // Alpha
      for(int x = 0; x < nidpA; x++) {
	int a = idprowA[x];
	int i = idpcolA[x];
	int h = idpirrA[x];
	wogA->set(x, WorbA->get(h, a + occpiA[h], i));
      }

    // find biggest_mograd
    biggest_mogradA=0;
    for (int i=0; i<nidpA;i++){
      if (fabs(wogA->get(i)) > biggest_mogradA)  biggest_mogradA=fabs(wogA->get(i));
    }

    // rms
    rms_wogA=0;
    for (int i=0; i<nidpA;i++) rms_wogA += wogA->get(i) * wogA->get(i);
    norm_wogA=sqrt(rms_wogA);
    rms_wogA=sqrt(rms_wogA)/nidpA;
    rms_wog=rms_wogA;

    // print
    if(print_ > 2){
      for(int i = 0; i < nidpA; i++){
        outfile->Printf("\n i, idpirrA, idprowA, idpcolA, wogA: %3d %3d %3d %3d %20.14f\n", i, idpirrA[i], idprowA[i],idpcolA[i],wogA->get(i));

      }
    }


}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
      // memalloc
      WorbA->zero();
      WorbB->zero();

      // set W matrix
      SharedMatrix tempA(GFockA->transpose());
      WorbA->copy(GFockA);
      WorbA->subtract(tempA);
      WorbA->scale(2.0);

      SharedMatrix tempB(GFockB->transpose());
      WorbB->copy(GFockB);
      WorbB->subtract(tempB);
      WorbB->scale(2.0);

      // Form w vector
      // Alpha
      for(int x = 0; x < nidpA; x++) {
	int a = idprowA[x];
	int i = idpcolA[x];
	int h = idpirrA[x];
	wogA->set(x, WorbA->get(h, a + occpiA[h], i));
      }

      // Beta
      for(int x = 0; x < nidpB; x++) {
	int a = idprowB[x];
	int i = idpcolB[x];
	int h = idpirrB[x];
	wogB->set(x, WorbB->get(h, a + occpiB[h], i));
      }


    // find biggest_mograd
    biggest_mogradA=0;
    for (int i=0; i<nidpA;i++){
      if (fabs(wogA->get(i)) > biggest_mogradA)  biggest_mogradA=fabs(wogA->get(i));
    }

    biggest_mogradB=0;
    for (int i=0; i<nidpB;i++){
      if (fabs(wogB->get(i)) > biggest_mogradB)  biggest_mogradB=fabs(wogB->get(i));
    }

    // rms
    rms_wogA=0;
    for (int i=0; i<nidpA;i++) rms_wogA += wogA->get(i) * wogA->get(i);
    norm_wogA=sqrt(rms_wogA);
    rms_wogA=sqrt(rms_wogA)/nidpA;

    rms_wogB=0;
    for (int i=0; i<nidpB;i++) rms_wogB += wogB->get(i) * wogB->get(i);
    norm_wogB=sqrt(rms_wogB);
    rms_wogB=sqrt(rms_wogB)/nidpB;
    rms_wog=MAX0(rms_wogA,rms_wogB);

    // print
    if(print_ > 2){
      for(int i = 0; i < nidpA; i++){
        outfile->Printf("\n i, idpirrA, idprowA, idpcolA, wogA: %3d %3d %3d %3d %20.14f\n", i, idpirrA[i], idprowA[i],idpcolA[i],wogA->get(i));

      }

      for(int i = 0; i < nidpB; i++){
        outfile->Printf("\n i, idpirrB, idprowB, idpcolB, wogB: %3d %3d %3d %3d %20.14f\n", i, idpirrB[i], idprowB[i],idpcolB[i],wogB->get(i));

      }
    }

}// end if (reference_ == "UNRESTRICTED")
}// end of mograd
}} // End Namespaces
