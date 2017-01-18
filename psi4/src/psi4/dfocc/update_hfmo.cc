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
#include "defines.h"
#include "dfocc.h"

using namespace psi;
using namespace std;

namespace psi{ namespace dfoccwave{

void DFOCC::update_hfmo()
{
//fprintf(outfile,"\n update_mo is starting... \n"); fflush(outfile);
//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {

/********************************************************************************************/
/************************** initialize array ************************************************/
/********************************************************************************************/
	UorbA->zero();
	KorbA->zero();

/********************************************************************************************/
/************************** Build kappa_bar *************************************************/
/********************************************************************************************/
        kappa_barA->add(kappaA);

/********************************************************************************************/
/************************** Construct Korb **************************************************/
/********************************************************************************************/
	// alpha
        for(int x = 0; x < nidpA; x++) {
	    int p = idprowA->get(x);
	    int q = idpcolA->get(x);
	    KorbA->set(p, q, kappa_barA->get(x));
	    KorbA->set(q, p, -kappa_barA->get(x));
        }

/********************************************************************************************/
/************************** Construct Uorb **************************************************/
/********************************************************************************************/
	//set to identity
	UorbA->identity();

	// K contribution
	UorbA->add(KorbA);

	//form K^2
	KsqrA->gemm(false, false, KorbA, KorbA, 1.0, 0.0);
	KsqrA->scale(0.5);

	// 0.5*K^2 contribution
	UorbA->add(KsqrA);

/********************************************************************************************/
/************************** Orthogonalize U matrix ******************************************/
/********************************************************************************************/
        if (orth_type == "MGS") UorbA->mgs();
        else if (orth_type == "GS") UorbA->gs();

/********************************************************************************************/
/************************** Build new MO coeff. *********************************************/
/********************************************************************************************/
	CmoA->gemm(false, false, Cmo_refA, UorbA, 1.0, 0.0);

       	if (print_ > 2) {
	  KorbA->print();
	  UorbA->print();
	  CmoA->print();
	}

     // build mo coeff blocks
     mo_coeff_blocks();

}// end if (reference_ == "RESTRICTED")


//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference_ == "UNRESTRICTED") {

/********************************************************************************************/
/************************** initialize array ************************************************/
/********************************************************************************************/
	UorbA->zero();
	UorbB->zero();
	KorbA->zero();
	KorbB->zero();

/********************************************************************************************/
/************************** Build kappa_bar *************************************************/
/********************************************************************************************/
        kappa_barA->add(kappaA);
        kappa_barB->add(kappaB);

/********************************************************************************************/
/************************** Construct Korb **************************************************/
/********************************************************************************************/
	// alpha
        for(int x = 0; x < nidpA; x++) {
	    int p = idprowA->get(x);
	    int q = idpcolA->get(x);
	    KorbA->set(p, q, kappa_barA->get(x));
	    KorbA->set(q, p, -kappa_barA->get(x));
        }

	// beta
        for(int x = 0; x < nidpB; x++) {
	    int p = idprowB->get(x);
	    int q = idpcolB->get(x);
	    KorbB->set(p, q, kappa_barB->get(x));
	    KorbB->set(q, p, -kappa_barB->get(x));
        }

/********************************************************************************************/
/************************** Construct Uorb **************************************************/
/********************************************************************************************/
	//set to identity
	UorbA->identity();
	UorbB->identity();

	// K contribution
	UorbA->add(KorbA);
	UorbB->add(KorbB);

	//form K^2
	KsqrA->gemm(false, false, KorbA, KorbA, 1.0, 0.0);
	KsqrB->gemm(false, false, KorbB, KorbB, 1.0, 0.0);
	KsqrA->scale(0.5);
	KsqrB->scale(0.5);

	// 0.5*K^2 contribution
	UorbA->add(KsqrA);
	UorbB->add(KsqrB);

/********************************************************************************************/
/************************** Orthogonalize U matrix ******************************************/
/********************************************************************************************/
        if (orth_type == "MGS") {
            UorbA->mgs();
            UorbB->mgs();
        }
        else if (orth_type == "GS") {
            UorbA->gs();
            UorbB->gs();
        }

/********************************************************************************************/
/************************** Build new MO coeff. *********************************************/
/********************************************************************************************/
	CmoA->gemm(false, false, Cmo_refA, UorbA, 1.0, 0.0);
	CmoB->gemm(false, false, Cmo_refB, UorbB, 1.0, 0.0);

       	if (print_ > 2) {
	  KorbA->print();
	  KorbB->print();
	  UorbA->print();
	  UorbB->print();
	  CmoA->print();
	  CmoB->print();
	}

     // build mo coeff blocks
     mo_coeff_blocks();

}// end if (reference_ == "UNRESTRICTED")

}// end main
}} // End Namespaces
