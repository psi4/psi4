/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libmints/matrix.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

/*
** DIIS: Direct inversion in the iterative subspace routine to
** accelerate convergence of the CCSD amplitude equations.
**
** Substantially improved efficiency of this routine:
** (1) Keeping at most two error vectors in core at once.
** (2) Limiting direct product (overlap) calculation to unique pairs.
** (3) Using LAPACK's linear equation solver DGESV instead of flin.
**
** -TDC  12/22/01
** -Modifications for ROHF and UHF, TDC, 6/03
**
** Condition Improvements: applying balanced, conditioned 
** pseudoinversion to prevent convergence errors
**
** -RMP 04/02/13
*/

void diis_RHF(int);
void diis_ROHF(int);
void diis_UHF(int);

void diis(int iter)
{
  if(params.ref == 0) diis_RHF(iter);
  else if(params.ref == 1) diis_ROHF(iter);
  else if(params.ref == 2) diis_UHF(iter);

  return;
}
void diis_invert_B(double** B, double* C, int dimension, double tolerance)
{
    SharedMatrix B2(new Matrix("B2", dimension, dimension));
    double** Bp = B2->pointer();
    ::memcpy((void*) Bp[0], B[0], sizeof(double) * dimension * dimension);

    double* Sp = new double[dimension];
    double* Tp = new double[dimension];

    bool is_zero = false;
    for (int i = 0; i < dimension - 1; i++) {
        if (Bp[i][i] <= 0.0) is_zero = true;
    }

    if (is_zero) {
        for (int i = 0; i < dimension; i++) {
            Sp[i] = 1.0;   
        }
    } else {
        for (int i = 0; i < dimension - 1; i++) {
            Sp[i] = pow(Bp[i][i],-1.0/2.0);
        }
        Sp[dimension - 1] = 1.0;
    }

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            Bp[i][j] *= Sp[i] * Sp[j]; 
        }
    }

    B2->power(-1.0, tolerance);
    
    C_DGEMV('N',dimension,dimension,1.0,Bp[0],dimension,C,1,0.0,Tp,1);
    
    for (int i = 0; i < dimension; i++) {
        C[i] = Sp[i] * Tp[i];
    } 
   
    delete[] Sp;
    delete[] Tp; 
}

}} // namespace psi::ccenergy
