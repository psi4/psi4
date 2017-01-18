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

#include <cmath>
#include <algorithm>

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "algebra_interface.h"
#include "blas.h"
#include "transform.h"
#include "matrix.h"

#define CCTRANSFORM_USE_BLAS

#define MAX(i,j) ((i>j) ? i : j)
#define MIN(i,j) ((i>j) ? j : i)
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define four(i,j,k,l) INDEX(INDEX(i,j),INDEX(k,l))

#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"

namespace psi{ namespace psimrcc{
    extern MOInfo *moinfo;

using namespace std;

/*!
    \fn CCTransform::read_oei_integrals()
 */
void CCTransform::read_oei_so_integrals()
{
  // Read all the (frozen + non-frozen) OEI in Pitzer order
  allocate_oei_so();

  int nso = moinfo->get_nso();

  double* H = new double[nso*(nso+1)/2];

  // Read the kinetic energy integrals
  for(int k=0; k<nso*(nso+1)/2;++k) H[k] = 0.0;
  iwl_rdone(PSIF_OEI,const_cast<char*>(PSIF_SO_T),H,nso*(nso+1)/2,0,0,"outfile");

  for(int i=0; i < nso; i++)
    for(int j=0; j < nso; j++)
      oei_so[i][j] = H[INDEX(i,j)];

  // Read the potential energy integrals
  for(int k=0; k<nso*(nso+1)/2;++k) H[k] = 0.0;
  iwl_rdone(PSIF_OEI,const_cast<char*>(PSIF_SO_V),H,nso*(nso+1)/2,0,0,"outfile");

  for(int i=0; i < nso; i++)
    for(int j=0; j < nso; j++)
      oei_so[i][j] += H[INDEX(i,j)];

  // Read the overlap integrals
  for(int k=0; k<nso*(nso+1)/2;++k) H[k] = 0.0;
  iwl_rdone(PSIF_OEI,const_cast<char*>(PSIF_SO_S),H,nso*(nso+1)/2,0,0,"outfile");

  for(int i=0; i < nso; i++)
    for(int j=0; j < nso; j++)
      s_so[i][j] += H[INDEX(i,j)];

  delete[] H;
}

}} /* End Namespaces */
