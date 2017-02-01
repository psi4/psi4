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

void transdip(MintsHelper &mints)
{
  int nmo, nso;
  double **scf_qt, **X;
  double **MUX_MO, **MUY_MO, **MUZ_MO;
  double **MUX_SO, **MUY_SO, **MUZ_SO;

  nso = moinfo.nso;
  nmo = moinfo.nmo;
  scf_qt = moinfo.scf_qt;

  /*** Transform the SO dipole integrals to the MO basis ***/
  vector<SharedMatrix> dipole = mints.so_dipole();
  MUX_SO = dipole[0]->to_block_matrix();
  MUY_SO = dipole[1]->to_block_matrix();
  MUZ_SO = dipole[2]->to_block_matrix();

  X = block_matrix(nmo,nso); /* just a temporary matrix */

  MUX_MO=block_matrix(nmo,nmo);
  MUY_MO=block_matrix(nmo,nmo);
  MUZ_MO=block_matrix(nmo,nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(MUX_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(MUX_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(MUY_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(MUY_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt[0][0]),nmo,&(MUZ_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(MUZ_MO[0][0]),nmo);

  free_block(X);

  moinfo.dip = (double ***) malloc(3 * sizeof(double **));
  moinfo.dip[0] = MUX_MO;
  moinfo.dip[1] = MUY_MO;
  moinfo.dip[2] = MUZ_MO;

  free_block(MUX_SO);
  free_block(MUY_SO);
  free_block(MUZ_SO);

  return;
}

}} // namespace psi::ccdensity
