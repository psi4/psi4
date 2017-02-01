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



using namespace psi;
using namespace std;


namespace psi{ namespace occwave{

void OCCWave::diis(int dimvec, Array2d *vecs, Array2d *errvecs, Array1d *vec_new, Array1d *errvec_new)
{

/********************************************************************************************/
/************************** memalloc ********************************************************/
/********************************************************************************************/
        Array2d *Bmat = new Array2d("DIIS B Matrix", nvar, nvar);
        Array1d *Cvec = new Array1d("DIIS C Vector", nvar);
        Array1d *vrow = new Array1d(dimvec);
        Array1d *vcol = new Array1d(dimvec);
        Bmat->zero();
        Cvec->zero();
        errvec_new->zero();
        vrow->zero();
        vcol->zero();

/********************************************************************************************/
/************************** Form B matrix ***************************************************/
/********************************************************************************************/
        // num_vecs/num_vecs part
	for(int i = 0; i < num_vecs; i++){
          vrow->row_vector(errvecs, i);
	  for(int j = 0; j < num_vecs; j++){
            vcol->row_vector(errvecs, j);
            double value = vrow->dot(vcol);
            Bmat->set(i, j, value);
	  }
	}

	for(int i = 0; i < num_vecs; i++){
          Bmat->set(nvar - 1, i, -1.0);
          Bmat->set(i, nvar - 1, -1.0);
	}

        Bmat->set(nvar - 1, nvar - 1, 0.0);

	// level shift
	if (level_shift == "TRUE") {
	  //for(int i = 0; i < num_vecs; i++) Bmat->set(i, i, -lshift_parameter);// this also works
	  for(int i = 0; i < num_vecs; i++) Bmat->set(i, i, Bmat->get(i, i) * (1 + lshift_parameter) );
	}

        // Form the c vector
	Cvec->set(nvar - 1,  -1.0);

/********************************************************************************************/
/************************** Solve LINEQ *****************************************************/
/********************************************************************************************/
         if (lineq == "CDGESV") Bmat->cdgesv(Cvec);
         else if (lineq == "FLIN") {
            double det = 0.0;
            Bmat->lineq_flin(Cvec, &det);
            if (fabs(det) < DIIS_MIN_DET) {
               outfile->Printf( "Warning!!! Diis matrix is near-singular\n");
               outfile->Printf( "Determinant is %6.3E\n", det);

            }
          }
          else if (lineq == "POPLE") Bmat->lineq_pople(Cvec, num_vecs, cutoff);

/********************************************************************************************/
/************************** Extrapolate *****************************************************/
/********************************************************************************************/
	for(int i = 0; i < dimvec; i++){
          double sum1 = 0.0;
          double sum2 = 0.0;
	  for(int j = 0; j < num_vecs; j++){
	    sum1 += Cvec->get(j) * vecs->get(j,i);
            sum2 += Cvec->get(j) * errvecs->get(j,i);
	  }
          vec_new->set(i, sum1);
          errvec_new->set(i, sum2);
	}

	delete Bmat;
	delete Cvec;
	delete vrow;
	delete vcol;

/********************************************************************************************/
/********************************************************************************************/

}
}} // End Namespaces
