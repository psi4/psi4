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
#include "psi4/libmints/matrix.h"
#include "occwave.h"


using namespace std;


namespace psi{ namespace occwave{

void OCCWave::update_mo()
{
//outfile->Printf("\n update_mo is starting... \n");
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
/************************ DO DIIS ***********************************************************/
/********************************************************************************************/
if (do_diis_ == 1) {

        // starting with itr = 1
        itr_diis++;

        // Form Diis Error Vector & Extrapolant Alpha Spin Case
	if (itr_diis <= num_vecs) {
	  for(int i = 0; i < nidpA; i++){
	    errvecsA->set(itr_diis-1, i, wogA->get(i));
	    vecsA->set(itr_diis-1, i, kappa_barA->get(i));
	  }
	}

	if (itr_diis > num_vecs) {
	  for(int j = 0; j < (num_vecs-1); j++){
	    for(int i = 0; i < nidpA; i++){
	      errvecsA->set(j, i, errvecsA->get(j+1, i));
	      vecsA->set(j, i, vecsA->get(j+1, i));
	    }
	  }

	  for(int i = 0; i < nidpA; i++){
	    errvecsA->set(num_vecs-1, i, wogA->get(i));
	    vecsA->set(num_vecs-1, i, kappa_barA->get(i));
	  }
	}

        // Extrapolate
        if (itr_diis >= num_vecs) {
	  diis(nidpA, vecsA, errvecsA, kappa_barA, wog_intA);
	}

}// end if (do_diis_ == 1)

/********************************************************************************************/
/************************** Construct Korb **************************************************/
/********************************************************************************************/
	// alpha
	for(int x = 0; x < nidpA; x++) {
	  int a = idprowA[x];
	  int i = idpcolA[x];
	  int h = idpirrA[x];
	  KorbA->set(h, a + occpiA[h], i, kappa_barA->get(x));
	  KorbA->set(h, i, a + occpiA[h], -kappa_barA->get(x));
	}

/********************************************************************************************/
/************************** Construct Uorb **************************************************/
/********************************************************************************************/
	//set to identity
	UorbA->identity();

	// K contribution
	UorbA->add(KorbA);

	//form K^2
	KsqrA->gemm(false, false, 1.0, KorbA, KorbA, 0.0);
	KsqrA->scale(0.5);

	// 0.5*K^2 contribution
	UorbA->add(KsqrA);

/********************************************************************************************/
/************************** Orthogonalize U matrix ******************************************/
/********************************************************************************************/
if (orth_type == "MGS") {;
    double rmgs1a,rmgs2a,rmgs1b,rmgs2b;

    // loop-over nirrep_
    for (int h=0; h<nirrep_; h++) {

      // loop-1
      for (int k = 0; k < nmopi_[h]; k++) {
	rmgs1a=0.0;

	// loop-1a
	for (int i=0; i < nmopi_[h]; i++) {
	  rmgs1a += UorbA->get(h, i, k) * UorbA->get(h, i, k);
	}// end 1a

	rmgs1a=sqrt(rmgs1a);

	// loop-1b
	for (int i=0; i < nmopi_[h]; i++) {
	  UorbA->set(h, i, k, UorbA->get(h, i, k) / rmgs1a);
	}// end 1b

	// loop-2
	for (int j=(k+1); j < nmopi_[h]; j++) {
	  rmgs2a=0;

	  // loop-2a
	  for (int i=0; i < nmopi_[h]; i++) {
	    rmgs2a += UorbA->get(h, i, k) * UorbA->get(h, i, j);
	  }// end 2a

	  // loop-2b
	  for (int i=0; i < nmopi_[h]; i++) {
	    UorbA->set(h, i, j, UorbA->get(h, i, j) - (rmgs2a * UorbA->get(h, i, k)));
	  }// end 2b

	}// end 2
      }// end 1
    }// end loop-over nirrep_
}// end main if


else if (orth_type == "GS") {
    int rowA = UorbA->nrow();
    int colA = UorbA->ncol();

    double **AdumA = block_matrix(rowA, colA);
    memset(AdumA[0], 0, sizeof(double)*rowA*colA);
    AdumA = UorbA->to_block_matrix();
    schmidt(AdumA, rowA, colA, "outfile");
    UorbA->set(AdumA);
    free_block(AdumA);
}

/********************************************************************************************/
/************************** Build new MO coeff. *********************************************/
/********************************************************************************************/
	Ca_->gemm(false, false, 1.0, Ca_ref, UorbA, 0.0);

       	if (print_ > 1) {
	  UorbA->print();
	  Ca_->print();
	}

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
/************************ DO DIIS ***********************************************************/
/********************************************************************************************/
if (do_diis_ == 1) {

        // starting with itr = 1
        itr_diis++;

        // Form Diis Error Vector & Extrapolant Alpha Spin Case
	if (itr_diis <= num_vecs) {
	  for(int i = 0; i < nidpA; i++){
	    errvecsA->set(itr_diis-1, i, wogA->get(i));
	    vecsA->set(itr_diis-1, i, kappa_barA->get(i));
	  }
	}


	if (itr_diis > num_vecs) {
	  for(int j = 0; j < (num_vecs-1); j++){
	    for(int i = 0; i < nidpA; i++){
	      errvecsA->set(j, i, errvecsA->get(j+1, i));
	      vecsA->set(j, i, vecsA->get(j+1, i));
	    }
	  }

	  for(int i = 0; i < nidpA; i++){
	    errvecsA->set(num_vecs-1, i, wogA->get(i));
	    vecsA->set(num_vecs-1, i, kappa_barA->get(i));
	  }
	}

        // Form Diis Error Vector & Extrapolant Beta Spin Case
	if (itr_diis <= num_vecs) {
	  for(int i = 0; i < nidpB; i++){
	    errvecsB->set(itr_diis-1, i, wogB->get(i));
	    vecsB->set(itr_diis-1, i, kappa_barB->get(i));
	  }
	}


	if (itr_diis > num_vecs) {
	  for(int j = 0; j < (num_vecs-1); j++){
	    for(int i = 0; i < nidpB; i++){
	      errvecsB->set(j, i, errvecsB->get(j+1, i));
	      vecsB->set(j, i, vecsB->get(j+1, i));
	    }
	  }

	  for(int i = 0; i < nidpB; i++){
	    errvecsB->set(num_vecs-1, i, wogB->get(i));
	    vecsB->set(num_vecs-1, i, kappa_barB->get(i));
	  }
	}


        // Extrapolate
        if (itr_diis >= num_vecs) {
	  diis(nidpA, vecsA, errvecsA, kappa_barA, wog_intA);
	  diis(nidpB, vecsB, errvecsB, kappa_barB, wog_intB);
	}

}// end if (do_diis_ == 1)

/********************************************************************************************/
/************************** Construct Korb **************************************************/
/********************************************************************************************/
	// alpha
	for(int x = 0; x < nidpA; x++) {
	  int a = idprowA[x];
	  int i = idpcolA[x];
	  int h = idpirrA[x];
	  KorbA->set(h, a + occpiA[h], i, kappa_barA->get(x));
	  KorbA->set(h, i, a + occpiA[h], -kappa_barA->get(x));
	}

	// beta
	for(int x = 0; x < nidpB; x++) {
	  int a = idprowB[x];
	  int i = idpcolB[x];
	  int h = idpirrB[x];
	  KorbB->set(h, a + occpiB[h], i, kappa_barB->get(x));
	  KorbB->set(h, i, a + occpiB[h], -kappa_barB->get(x));
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
	KsqrA->gemm(false, false, 1.0, KorbA, KorbA, 0.0);
	KsqrB->gemm(false, false, 1.0, KorbB, KorbB, 0.0);
	KsqrA->scale(0.5);
	KsqrB->scale(0.5);

	// 0.5*K^2 contribution
	UorbA->add(KsqrA);
	UorbB->add(KsqrB);

/********************************************************************************************/
/************************** Orthogonalize U matrix ******************************************/
/********************************************************************************************/
if (orth_type == "MGS") {;
    double rmgs1a,rmgs2a,rmgs1b,rmgs2b;

    // loop-over nirrep_
    for (int h=0; h<nirrep_; h++) {

      // loop-1
      for (int k = 0; k < nmopi_[h]; k++) {
	rmgs1a=0.0;
	rmgs1b=0.0;

	// loop-1a
	for (int i=0; i < nmopi_[h]; i++) {
	  rmgs1a += UorbA->get(h, i, k) * UorbA->get(h, i, k);
	  rmgs1b += UorbB->get(h, i, k) * UorbB->get(h, i, k);
	}// end 1a

	rmgs1a=sqrt(rmgs1a);
	rmgs1b=sqrt(rmgs1b);

	// loop-1b
	for (int i=0; i < nmopi_[h]; i++) {
	  UorbA->set(h, i, k, UorbA->get(h, i, k) / rmgs1a);
	  UorbB->set(h, i, k, UorbB->get(h, i, k) / rmgs1b);
	}// end 1b

	// loop-2
	for (int j=(k+1); j < nmopi_[h]; j++) {
	  rmgs2a=0;
	  rmgs2b=0;

	  // loop-2a
	  for (int i=0; i < nmopi_[h]; i++) {
	    rmgs2a += UorbA->get(h, i, k) * UorbA->get(h, i, j);
	    rmgs2b += UorbB->get(h, i, k) * UorbB->get(h, i, j);
	  }// end 2a

	  // loop-2b
	  for (int i=0; i < nmopi_[h]; i++) {
	    UorbA->set(h, i, j, UorbA->get(h, i, j) - (rmgs2a * UorbA->get(h, i, k)));
	    UorbB->set(h, i, j, UorbB->get(h, i, j) - (rmgs2b * UorbB->get(h, i, k)));
	  }// end 2b

	}// end 2
      }// end 1
    }// end loop-over nirrep_
}// end main if


else if (orth_type == "GS") {
    int rowA = UorbA->nrow();
    int colA = UorbA->ncol();

    double **AdumA = block_matrix(rowA, colA);
    memset(AdumA[0], 0, sizeof(double)*rowA*colA);
    AdumA = UorbA->to_block_matrix();
    schmidt(AdumA, rowA, colA, "outfile");
    UorbA->set(AdumA);
    free_block(AdumA);

    int rowB = UorbB->nrow();
    int colB = UorbB->ncol();

    double **AdumB = block_matrix(rowB, colB);
    memset(AdumB[0], 0, sizeof(double)*rowB*colB);
    AdumB = UorbB->to_block_matrix();
    schmidt(AdumB, rowB, colB, "outfile");
    UorbB->set(AdumB);
    free_block(AdumB);
}

/********************************************************************************************/
/************************** Build new MO coeff. *********************************************/
/********************************************************************************************/
	Ca_->gemm(false, false, 1.0, Ca_ref, UorbA, 0.0);
	Cb_->gemm(false, false, 1.0, Cb_ref, UorbB, 0.0);

       	if (print_ > 1) {
	  UorbA->print();
          UorbB->print();
	  Ca_->print();
	  Cb_->print();
	}

}// end if (reference_ == "UNRESTRICTED")

}// end main
}} // End Namespaces
