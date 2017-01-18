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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

using namespace std;

namespace psi { namespace ccdensity {

void transL(MintsHelper &mints, double sign)
{
  int nmo, nso;
  double **scf_qt, **X;
  double **LX_MO, **LY_MO, **LZ_MO;
  double **LX_SO, **LY_SO, **LZ_SO;

  nso = moinfo.nso;
  nmo = moinfo.nmo;
  scf_qt = moinfo.scf_qt;

  /*** Transform the SO nabla integrals to the MO basis ***/
  vector<SharedMatrix> dipole = mints.so_angular_momentum();
  for(int i=0; i < 3; i++) dipole[i]->scale(-0.5 * sign);
  LX_SO = dipole[0]->to_block_matrix();
  LY_SO = dipole[1]->to_block_matrix();
  LZ_SO = dipole[2]->to_block_matrix();

  X = block_matrix(nmo,nso); /* just a temporary matrix */

  LX_MO=block_matrix(nmo,nmo);
  LY_MO=block_matrix(nmo,nmo);
  LZ_MO=block_matrix(nmo,nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(LX_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(LX_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(LY_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(LY_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(LZ_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(LZ_MO[0][0]),nmo);

  free_block(X);

  moinfo.L = (double ***) malloc(3 * sizeof(double **));
  moinfo.L[0] = LX_MO;
  moinfo.L[1] = LY_MO;
  moinfo.L[2] = LZ_MO;

  free_block(LX_SO);
  free_block(LY_SO);
  free_block(LZ_SO);

  return;
}

}} // namespace psi::ccdensity
