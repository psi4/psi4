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
    \ingroup ccresponse
    \brief Compute the three tensors needed for Raman Optical Activity.

    ROA requires the following polarizability tensors:
      (1) electric-dipole/electric-dipole; 
      (2) electric-dipole/electric-quadrupole; and 
      (3) electric-dipole/magnetic-dipole.

  -TDC, August 2009
*/
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <physconst.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"
#include <vector>
#include <psi4-dec.h>
#include <libmints/mints.h>

namespace psi { namespace ccresponse {

void scatter(std::vector <SharedMatrix> dip, std::vector <SharedMatrix> rot, std::vector <SharedMatrix> quad)
{
    // COMPUTE TENSOR DERIVATIVES
    // Replicate the part of the roa.pl code that does this here in C++
    // Dipole Polarizability Tensor Gradients
    double step = 0.001;
    SharedMatrix denom(new Matrix(3,3));
    denom->set(2.0 * step);
    std::vector <SharedMatrix> dip_grad;
    for(std::vector<SharedMatrix>::iterator it=dip.begin(); it != dip.end(); ++it) {
        SharedMatrix grad_mat(new Matrix(3,3));
        grad_mat->add(*it);
        ++it;
        grad_mat->subtract(*it);
        grad_mat->apply_denominator(denom);
        grad_mat->print(stdout);
        dip_grad.push_back(grad_mat);
    }
    
    
    // COMPUTE SCATTERING
    // Update the roa.c code in ~crawdad area on cerebro to use SharedMatrix and
    // other psi4 style code
    
    
    printf("FUN-ction");

}

}} // namespace psi::ccresponse

