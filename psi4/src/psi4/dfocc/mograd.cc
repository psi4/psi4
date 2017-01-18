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

#include "psi4/psi4-dec.h"
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::mograd()
{
 
      double norm_wogA, norm_wogB;    

if (reference_ == "RESTRICTED") {
      // memalloc 
      WorbA->zero();

      // set W matrix
      SharedTensor2d temp(GF->transpose());
      WorbA->copy(GF);
      WorbA->subtract(temp);       
      // DO NOT scale W since Fai (spin-free) = 2 * FAI
      //WorbA->scale(2.0);
      //WorbA->print();
      
      // Form w vector
      // Alpha
      for(int x = 0; x < nidpA; x++) {
	int p = idprowA->get(x);
	int q = idpcolA->get(x);
	wogA->set(x, WorbA->get(p,q));
      }
      //wogA->print();
    
    // find biggest_mograd 
    biggest_mogradA=0;
    for (int i=0; i<nidpA;i++){ 
      if (fabs(wogA->get(i)) > biggest_mogradA)  biggest_mogradA=fabs(wogA->get(i));
    }
      
    // rms
    rms_wogA=0;
    for (int i=0; i<nidpA;i++) rms_wogA += wogA->get(i) * wogA->get(i);
    norm_wogA=sqrt(rms_wogA);
    rms_wogA = wogA->rms();
    rms_wog=rms_wogA;  
    
    // print
    if(print_ > 2){
      for(int i = 0; i < nidpA; i++){
        outfile->Printf("\ti, idprowA, idpcolA, wogA: %3d %3d %3d %20.14f\n", i, idprowA->get(i), idpcolA->get(i), wogA->get(i)); 
      }
    }


}// end if (reference_ == "RESTRICTED") 

else if (reference_ == "UNRESTRICTED") {
      // set W matrix
      SharedTensor2d tempA(GFA->transpose());
      WorbA->copy(GFA);
      WorbA->subtract(tempA);       
      WorbA->scale(2.0);
      //WorbA->print();
      
      SharedTensor2d tempB(GFB->transpose());
      WorbB->copy(GFB);
      WorbB->subtract(tempB);       
      WorbB->scale(2.0);
      //WorbB->print();

      // Form w vector
      // Alpha
      for(int x = 0; x < nidpA; x++) {
	int p = idprowA->get(x);
	int q = idpcolA->get(x);
	wogA->set(x, WorbA->get(p,q));
      }
      //wogA->print();
      
      // Beta
      for(int x = 0; x < nidpB; x++) {
	int p = idprowB->get(x);
	int q = idpcolB->get(x);
	wogB->set(x, WorbB->get(p,q));
      }
      //wogB->print();
    
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
    rms_wogA = wogA->rms();
    
    rms_wogB=0;
    for (int i=0; i<nidpB;i++) rms_wogB += wogB->get(i) * wogB->get(i);
    norm_wogB=sqrt(rms_wogB);
    rms_wogB = wogB->rms();
    rms_wog=MAX0(rms_wogA,rms_wogB);
    
    // print
    if(print_ > 2){
      for(int i = 0; i < nidpA; i++){
        outfile->Printf("\ti, idprowA, idpcolA, wogA: %3d %3d %3d %20.14f\n", i, idprowA->get(i), idpcolA->get(i), wogA->get(i)); 
	
      }
      for(int i = 0; i < nidpB; i++){
        outfile->Printf("\ti, idprowB, idpcolB, wogB: %3d %3d %3d %20.14f\n", i, idprowB->get(i), idpcolB->get(i), wogB->get(i)); 
	
      }
    }

}// end if (reference_ == "UNRESTRICTED") 

}// end of main
}} // End Namespaces
