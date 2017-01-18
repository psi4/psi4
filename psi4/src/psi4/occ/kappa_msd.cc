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
#include "psi4/libmints/matrix.h"


using namespace std;


namespace psi{ namespace occwave{

void OCCWave::kappa_msd()
{
//outfile->Printf("\n kappa_msd is starting... \n");

if (reference_ == "RESTRICTED") {
        // Get kappa
	for(int x = 0; x < nidpA; x++) {
	  int a = idprowA[x];
	  int i = idpcolA[x];
	  int h = idpirrA[x];
	  double value = FockA->get(h, a + occpiA[h], a + occpiA[h]) - FockA->get(h, i, i);
	  kappaA->set(x, -wogA->get(x) / (2*value));
	}

        // find biggest_kappa
	biggest_kappaA=0;
	for (int i=0; i<nidpA;i++) {
	    if (fabs(kappaA->get(i)) > biggest_kappaA) biggest_kappaA=fabs(kappaA->get(i));
	}

        // Scale
	if (biggest_kappaA > step_max) {
	    for (int i=0; i<nidpA;i++) kappaA->set(i, kappaA->get(i) *(step_max/biggest_kappaA));
	}


        // find biggest_kappa again
	if (biggest_kappaA > step_max)
	{
	  biggest_kappaA=0;
	  for (int i=0; i<nidpA;i++)
	  {
	      if (fabs(kappaA->get(i)) > biggest_kappaA)
	      {
		  biggest_kappaA = fabs(kappaA->get(i));
	      }
	  }
	}

        // norm
	rms_kappaA=0;
	rms_kappaA = kappaA->rms();

        // print
        if(print_ > 2) {
           kappaA->print();
           kappa_barA->print();
        }

}// end if (reference_ == "RESTRICTED")

else if (reference_ == "UNRESTRICTED") {
        // Get kappa
	// alpha
	for(int x = 0; x < nidpA; x++) {
	  int a = idprowA[x];
	  int i = idpcolA[x];
	  int h = idpirrA[x];
	  double value = FockA->get(h, a + occpiA[h], a + occpiA[h]) - FockA->get(h, i, i);
	  kappaA->set(x, -wogA->get(x) / (2*value));
	}

	// beta
	for(int x = 0; x < nidpB; x++) {
	  int a = idprowB[x];
	  int i = idpcolB[x];
	  int h = idpirrB[x];
	  double value = FockB->get(h, a + occpiB[h], a + occpiB[h]) - FockB->get(h, i, i);
	  kappaB->set(x, -wogB->get(x) / (2*value));
	}

        // find biggest_kappa
	biggest_kappaA=0;
	for (int i=0; i<nidpA;i++) {
	    if (fabs(kappaA->get(i)) > biggest_kappaA) biggest_kappaA=fabs(kappaA->get(i));
	}

	biggest_kappaB=0;
	for (int i=0; i<nidpB;i++){
	    if (fabs(kappaB->get(i)) > biggest_kappaB) biggest_kappaB=fabs(kappaB->get(i));
	}

        // Scale
	if (biggest_kappaA > step_max) {
	    for (int i=0; i<nidpA;i++) kappaA->set(i, kappaA->get(i) *(step_max/biggest_kappaA));
	}

	if (biggest_kappaB > step_max) {
	    for (int i=0; i<nidpB;i++) kappaB->set(i, kappaB->get(i) *(step_max/biggest_kappaB));
	}

        // find biggest_kappa again
	if (biggest_kappaA > step_max)
	{
	  biggest_kappaA=0;
	  for (int i=0; i<nidpA;i++)
	  {
	      if (fabs(kappaA->get(i)) > biggest_kappaA)
	      {
		  biggest_kappaA = fabs(kappaA->get(i));
	      }
	  }
	}

	if (biggest_kappaB > step_max)
	{
	  biggest_kappaB=0;
	  for (int i=0; i<nidpB;i++)
	  {
	      if (fabs(kappaB->get(i)) > biggest_kappaB)
	      {
		  biggest_kappaB=fabs(kappaB->get(i));
	      }
	  }
	}

        // norm
	rms_kappaA=0;
	rms_kappaB=0;
	rms_kappaA = kappaA->rms();
	rms_kappaB = kappaB->rms();

        // print
        if(print_ > 2){
          kappaA->print();
          kappaB->print();
        }

}// end if (reference_ == "UNRESTRICTED")
 //outfile->Printf("\n kappa_msd done. \n");

}// end main
}} // End Namespaces
