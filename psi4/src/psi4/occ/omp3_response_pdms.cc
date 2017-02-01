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
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/matrix.h"
#include "occwave.h"


using namespace std;

namespace psi{ namespace occwave{

void OCCWave::omp3_response_pdms()
{
        //outfile->Printf("\n response_pdms is starting... \n");

        // Build G intermediates
        timer_on("G int");
	omp3_g_int();
        timer_off("G int");

 if (reference_ == "RESTRICTED") {
        // Initialize
	gamma1corr->zero();
	g1symm->zero();

        // OPDM
        timer_on("OPDM");
	// OO-block alpha contrb.
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
                g1symm->set(h, i, j, GooA->get(h, i, j) + GooA->get(h, j, i));
            }
	  }
	}

	// VV-block alpha contrb.
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < avirtpiA[h]; ++a){
            for(int b = 0 ; b < avirtpiA[h]; ++b){
                int aa = a + occpiA[h];
                int bb = b + occpiA[h];
                g1symm->set(h, aa, bb, GvvA->get(h, a, b) + GvvA->get(h, b, a));
            }
	  }
	}

	g1symm->scale(-1.0);
	gamma1corr->copy(g1symm); // correlation opdm

        // REF contribution
	// alpha contrb.
        #pragma omp parallel for
	for(int h=0; h<nirrep_; h++) {
	  if (occpiA[h] != 0) {
	    for (int i=0; i<occpiA[h];i++) {
	      g1symm->add(h,i,i,2.0);
	    }
	  }
	}
        timer_off("OPDM");

        //print
        if (print_ > 1) {
	  g1symm->print();
        }

        // TPDM
        timer_on("V int");
        v_2nd_order();
        timer_off("V int");
        timer_on("TPDM OOVV");
	tpdm_oovv();
        timer_off("TPDM OOVV");
        timer_on("TPDM OOOO");
	tpdm_oooo();
        timer_off("TPDM OOOO");

        if (twopdm_abcd_type == "COMPUTE") {
           timer_on("TPDM VVVV");
           omp3_tpdm_vvvv();
           timer_off("TPDM VVVV");
        }

        timer_on("TPDM OVOV");
        tpdm_ovov();
        timer_off("TPDM OVOV");
        timer_on("TPDM REF");
	tpdm_ref();
        timer_off("TPDM REF");
        timer_on("TPDM CORR OPDM");
	tpdm_corr_opdm();
        timer_off("TPDM CORR OPDM");
 }// end if (reference_ == "RESTRICTED")

 else if (reference_ == "UNRESTRICTED") {
        // Initialize
	gamma1corrA->zero();
	gamma1corrB->zero();
	g1symmA->zero();
	g1symmB->zero();

        // OPDM
        timer_on("OPDM");
	// OO-block alpha contrb.
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
                g1symmA->set(h, i, j, GooA->get(h, i, j) + GooA->get(h, j, i));
            }
	  }
	}

	// OO-block beta contrb.
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int j = 0 ; j < aoccpiB[h]; ++j){
                g1symmB->set(h, i, j, GooB->get(h, i, j) + GooB->get(h, j, i));
            }
	  }
	}

	// VV-block alpha contrb.
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < avirtpiA[h]; ++a){
            for(int b = 0 ; b < avirtpiA[h]; ++b){
                int aa = a + occpiA[h];
                int bb = b + occpiA[h];
                g1symmA->set(h, aa, bb, GvvA->get(h, a, b) + GvvA->get(h, b, a));
            }
	  }
	}

        // VV-block beta contrb.
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < avirtpiB[h]; ++a){
            for(int b = 0 ; b < avirtpiB[h]; ++b){
                int aa = a + occpiB[h];
                int bb = b + occpiB[h];
                g1symmB->set(h, aa, bb, GvvB->get(h, a, b) + GvvB->get(h, b, a));
            }
	  }
	}

	g1symmA->scale(-0.5);
	g1symmB->scale(-0.5);
	gamma1corrA->copy(g1symmA); // correlation opdm
	gamma1corrB->copy(g1symmB); // correlation opdm

        // REF contribution
	// alpha contrb.
        #pragma omp parallel for
	for(int h=0; h<nirrep_; h++) {
	  if (occpiA[h] != 0) {
	    for (int i=0; i<occpiA[h];i++) {
	      g1symmA->add(h,i,i,1.0);
	    }
	  }
	}

	// beta contrb.
        #pragma omp parallel for
	for(int h=0; h<nirrep_; h++) {
	  if (occpiB[h] != 0) {
	    for (int i=0; i<occpiB[h];i++) {
	      g1symmB->add(h,i,i,1.0);
	    }
	  }
	}
        timer_off("OPDM");

        //print
        if (print_ > 1) {
	  g1symmA->print();
	  g1symmB->print();
        }

        // TPDM
        timer_on("V int");
        v_2nd_order();
        timer_off("V int");
        timer_on("TPDM OOVV");
	tpdm_oovv();
        timer_off("TPDM OOVV");
        timer_on("TPDM OOOO");
	tpdm_oooo();
        timer_off("TPDM OOOO");

        if (twopdm_abcd_type == "COMPUTE") {
           timer_on("TPDM VVVV");
           omp3_tpdm_vvvv();
           timer_off("TPDM VVVV");
        }

        timer_on("TPDM OVOV");
        tpdm_ovov();
        timer_off("TPDM OVOV");
        timer_on("TPDM VOVO");
        tpdm_vovo();
        timer_off("TPDM VOVO");
        timer_on("TPDM OVVO");
        tpdm_ovvo();
        timer_off("TPDM OVVO");
        timer_on("TPDM REF");
	tpdm_ref();
        timer_off("TPDM REF");
        timer_on("TPDM CORR OPDM");
	tpdm_corr_opdm();
        timer_off("TPDM CORR OPDM");



 }// end if (reference_ == "UNRESTRICTED")

  //outfile->Printf("\n response_pdms done... \n");

} // end of response_pdms
}} // End Namespaces


