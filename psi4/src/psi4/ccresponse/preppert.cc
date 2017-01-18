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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/basisset.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

using namespace std;

namespace psi { namespace ccresponse {

void sort_pert(const char *pert, double **pertints, int irrep);

/* preppert(): Prepare DPD structures for all currently needed one-electron
** property integrals in the MO basis.
**
** -TDC, 6/11
*/

void preppert(std::shared_ptr<BasisSet> primary)
{
  int i, j, ij;

  char **cartcomp = (char **) malloc(3 * sizeof(char *));
  cartcomp[0] = strdup("X");
  cartcomp[1] = strdup("Y");
  cartcomp[2] = strdup("Z");
  char lbl[32];

  MintsHelper mints(primary, Process::environment.options, 0);
  vector<SharedMatrix> dipole = mints.so_dipole();
  vector<SharedMatrix> nabla = mints.so_nabla();
  vector<SharedMatrix> angmom = mints.so_angular_momentum();
  vector<SharedMatrix> trquad = mints.so_traceless_quadrupole();

  int nso = moinfo.nso;
  int nmo = moinfo.nmo;

  double **TMP2 = block_matrix(nso,nso);

  // Electric dipole integrals
  for(i=0; i < 3; i++) {
    double **TMP1 = dipole[i]->to_block_matrix();
    double **TMP3 = block_matrix(nmo, nmo);
    C_DGEMM('n','n',nso,nmo,nso,1,TMP1[0],nso,moinfo.scf[0],nmo,0,TMP2[0],nso);
    C_DGEMM('t','n',nmo,nmo,nso,1,moinfo.scf[0],nmo,TMP2[0],nso,0,TMP3[0],nmo);
    moinfo.MU[i] = TMP3;
    sprintf(lbl, "Mu_%1s", cartcomp[i]);
    sort_pert(lbl, moinfo.MU[i], moinfo.mu_irreps[i]);
  }

  // Velocity-gauge electric dipole integrals
  for(i=0; i < 3; i++) {
    double **TMP1 = nabla[i]->to_block_matrix();
    double **TMP3 = block_matrix(nmo, nmo);
    C_DGEMM('n','n',nso,nmo,nso,1,TMP1[0],nso,moinfo.scf[0],nmo,0,TMP2[0],nso);
    C_DGEMM('t','n',nmo,nmo,nso,1,moinfo.scf[0],nmo,TMP2[0],nso,0,TMP3[0],nmo);
    moinfo.P[i] = TMP3;
    sprintf(lbl, "P_%1s", cartcomp[i]);
    sort_pert(lbl, moinfo.P[i], moinfo.mu_irreps[i]);
  }

  // Complex conjugate of velocity-gauge electric dipole integrals
  for(i=0; i < 3; i++) nabla[i]->scale(-1.0);
  for(i=0; i < 3; i++) {
    double **TMP1 = nabla[i]->to_block_matrix();
    double **TMP3 = block_matrix(nmo, nmo);
    C_DGEMM('n','n',nso,nmo,nso,1,TMP1[0],nso,moinfo.scf[0],nmo,0,TMP2[0],nso);
    C_DGEMM('t','n',nmo,nmo,nso,1,moinfo.scf[0],nmo,TMP2[0],nso,0,TMP3[0],nmo);
    moinfo.Pcc[i] = TMP3;
    sprintf(lbl, "P*_%1s", cartcomp[i]);
    sort_pert(lbl, moinfo.Pcc[i], moinfo.mu_irreps[i]);
  }

  // Magnetic dipole integrals (these require a -1/2 prefactor)
  for(i=0; i < 3; i++) {
    angmom[i]->scale(-0.5);
    double **TMP1 = angmom[i]->to_block_matrix();
    double **TMP3 = block_matrix(nmo, nmo);
    sprintf(lbl, "L_%1s", cartcomp[i]);
//    outfile->Printf( "%s Angular Momentum Integrals (SO)\n",lbl);
//    mat_print(TMP1,nmo, nmo, "outfile");
    C_DGEMM('n','n',nso,nmo,nso,1,TMP1[0],nso,moinfo.scf[0],nmo,0,TMP2[0],nso);
    C_DGEMM('t','n',nmo,nmo,nso,1,moinfo.scf[0],nmo,TMP2[0],nso,0,TMP3[0],nmo);
    moinfo.L[i] = TMP3;
    sort_pert(lbl, moinfo.L[i], moinfo.l_irreps[i]);
  }

  // Complex conjugate of magnetic dipole integrals
  for(i=0; i < 3; i++) angmom[i]->scale(-1.0);
  for(i=0; i < 3; i++) {
    double **TMP1 = angmom[i]->to_block_matrix();
    double **TMP3 = block_matrix(nmo, nmo);
    C_DGEMM('n','n',nso,nmo,nso,1,TMP1[0],nso,moinfo.scf[0],nmo,0,TMP2[0],nso);
    C_DGEMM('t','n',nmo,nmo,nso,1,moinfo.scf[0],nmo,TMP2[0],nso,0,TMP3[0],nmo);
    moinfo.Lcc[i] = TMP3;
    sprintf(lbl, "L*_%1s", cartcomp[i]);
    sort_pert(lbl, moinfo.Lcc[i], moinfo.l_irreps[i]);
  }

  // Traceless quadrupole integrals
  for(i=0,ij=0; i < 3; i++) {
    for(j=i; j < 3; j++,ij++) {
      double **TMP1 = trquad[ij]->to_block_matrix();
      double **TMP3 = block_matrix(nmo, nmo);
      C_DGEMM('n','n',nso,nmo,nso,1,TMP1[0],nso,moinfo.scf[0],nmo,0,TMP2[0],nso);
      C_DGEMM('t','n',nmo,nmo,nso,1,moinfo.scf[0],nmo,TMP2[0],nso,0,TMP3[0],nmo);
      moinfo.Q[i][j] = TMP3;
      sprintf(lbl, "Q_%1s%1s", cartcomp[i], cartcomp[j]);
      sort_pert(lbl, moinfo.Q[i][j], moinfo.mu_irreps[i]^moinfo.mu_irreps[j]);
      if(i!=j) {
        moinfo.Q[j][i] = TMP3;
        sprintf(lbl, "Q_%1s%1s", cartcomp[j], cartcomp[i]);
        sort_pert(lbl, moinfo.Q[j][i], moinfo.mu_irreps[j]^moinfo.mu_irreps[i]);
      }
    }
  }

  free_block(TMP2);
}

}} // namespace psi::ccresponse
