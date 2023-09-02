/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "defines.h"
#include "dfocc.h"
#include "psi4/psi4-dec.h"

#include <cmath>

using namespace psi;

namespace psi {
namespace dfoccwave {

void DFOCC::diis(int dimvec, SharedTensor2d &vecs, SharedTensor2d &errvecs, SharedTensor1d &vec_new,
                 SharedTensor1d &errvec_new) {
    /********************************************************************************************/
    /************************** memalloc ********************************************************/
    /********************************************************************************************/
    SharedTensor2d Bmat = std::make_shared<Tensor2d>("DIIS B Matrix", nvar, nvar);
    SharedTensor1d Cvec = std::make_shared<Tensor1d>("DIIS C Vector", nvar);
    SharedTensor1d vrow = std::make_shared<Tensor1d>("DIIS vrow", dimvec);
    SharedTensor1d vcol = std::make_shared<Tensor1d>("DIIS vcol", dimvec);

    /********************************************************************************************/
    /************************** Form B matrix ***************************************************/
    /********************************************************************************************/
    // num_vecs/num_vecs part
    for (int i = 0; i < num_vecs; i++) {
        vrow->row_vector(errvecs, i);
        for (int j = 0; j < num_vecs; j++) {
            vcol->row_vector(errvecs, j);
            double value = vrow->dot(vcol);
            Bmat->set(i, j, value);
        }
    }

    for (int i = 0; i < num_vecs; i++) {
        Bmat->set(nvar - 1, i, -1.0);
        Bmat->set(i, nvar - 1, -1.0);
    }

    Bmat->set(nvar - 1, nvar - 1, 0.0);

    // level shift
    if (level_shift == "TRUE") {
        // for(int i = 0; i < num_vecs; i++) Bmat->set(i, i, -lshift_parameter);// this also works
        for (int i = 0; i < num_vecs; i++) Bmat->set(i, i, Bmat->get(i, i) * (1 + lshift_parameter));
    }

    // Form the c vector
    Cvec->set(nvar - 1, -1.0);

    /********************************************************************************************/
    /************************** Solve LINEQ *****************************************************/
    /********************************************************************************************/
    if (lineq == "CDGESV")
        Bmat->cdgesv(Cvec);
    else if (lineq == "FLIN") {
        double det = 0.0;
        Bmat->lineq_flin(Cvec, &det);
        if (std::fabs(det) < DIIS_MIN_DET) {
            outfile->Printf("Warning!!! Diis matrix is near-singular\n");
            outfile->Printf("Determinant is %6.3E\n", det);
        }
    } else if (lineq == "POPLE")
        Bmat->lineq_pople(Cvec, num_vecs, cutoff);

    /********************************************************************************************/
    /************************** Extrapolate *****************************************************/
    /********************************************************************************************/
    for (int i = 0; i < dimvec; i++) {
        double sum1 = 0.0;
        double sum2 = 0.0;
        for (int j = 0; j < num_vecs; j++) {
            sum1 += Cvec->get(j) * vecs->get(j, i);
            sum2 += Cvec->get(j) * errvecs->get(j, i);
        }
        vec_new->set(i, sum1);
        errvec_new->set(i, sum2);
    }

    /********************************************************************************************/
    /********************************************************************************************/
}
}  // namespace dfoccwave
}  // namespace psi
