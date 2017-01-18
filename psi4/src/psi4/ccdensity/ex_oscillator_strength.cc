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
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/mintshelper.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

using namespace std;

namespace psi { namespace ccdensity {
#include "psi4/physconst.h"

void ex_oscillator_strength(SharedWavefunction wfn, struct TD_Params *S, struct TD_Params *U, struct XTD_Params *xtd_data)
{
  int nmo, nso, i, I, h, j, nirreps;
  int *order, *order_A, *order_B, *doccpi, *clsdpi, *openpi, *orbspi;
  double **scf_pitzer, **scf_pitzer_A, **scf_pitzer_B;
  double **scf_qt, **scf_qt_A, **scf_qt_B, **X;
  double *mu_x_ints, *mu_y_ints, *mu_z_ints;
  double **MUX_MO, **MUY_MO, **MUZ_MO;
  double **MUX_SO, **MUY_SO, **MUZ_SO;
  double **MUX_MO_A, **MUY_MO_A, **MUZ_MO_A;
  double **MUX_MO_B, **MUY_MO_B, **MUZ_MO_B;
  double lt_x, lt_y, lt_z;
  double rt_x, rt_y, rt_z;
  double ds_x, ds_y, ds_z;
  double f_x, f_y, f_z;
  double f;
  double delta_ee;
  double einstein_a, einstein_b;

  if ((params.ref == 0) || (params.ref == 1))
    scf_pitzer = wfn->Ca()->to_block_matrix();
  else if(params.ref == 2) {
    scf_pitzer_A = wfn->Ca()->to_block_matrix();
    scf_pitzer_B = wfn->Cb()->to_block_matrix();
  }

  nso = wfn->nso();
  nmo = wfn->nmo();

  lt_x = lt_y = lt_z = 0.0;
  rt_x = rt_y = rt_z = 0.0;
  ds_x = ds_y = ds_z = 0.0;
  f_x = f_y = f_z = 0.0;
  f = 0;

  doccpi = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++)
    doccpi[h] = moinfo.frdocc[h] + moinfo.clsdpi[h];

  if((params.ref == 0) || (params.ref == 1)) {
    order = init_int_array(nmo);

    reorder_qt(doccpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
               order, moinfo.orbspi, moinfo.nirreps);

    scf_qt = block_matrix(nmo, nmo);
    for(i=0; i < nmo; i++) {
      I = order[i];  /* Pitzer --> QT */
      for(j=0; j < nmo; j++) scf_qt[j][I] = scf_pitzer[j][i];
    }

    free(order);
    free_block(scf_pitzer);
  }
  else if(params.ref == 2) {
    order_A = init_int_array(nmo);
    order_B = init_int_array(nmo);

    reorder_qt_uhf(doccpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
                  order_A,order_B, moinfo.orbspi, moinfo.nirreps);

    scf_qt_A = block_matrix(nmo, nmo);
    for(i=0; i < nmo; i++) {
      I = order_A[i];  /* Pitzer --> QT */
      for(j=0; j < nmo; j++) scf_qt_A[j][I] = scf_pitzer_A[j][i];
    }

    scf_qt_B = block_matrix(nmo, nmo);
    for(i=0; i < nmo; i++) {
      I = order_B[i];  /* Pitzer --> QT */
      for(j=0; j < nmo; j++) scf_qt_B[j][I] = scf_pitzer_B[j][i];
    }

    free(order_A);
    free(order_B);
    free_block(scf_pitzer_A);
    free_block(scf_pitzer_B);
  }
  free(doccpi);

  /*** Transform the SO dipole integrals to the MO basis ***/

  MintsHelper mints(wfn->basisset(), Process::environment.options, 0);
  vector<SharedMatrix> dipole = mints.so_dipole();
  MUX_SO = dipole[0]->to_block_matrix();
  MUY_SO = dipole[1]->to_block_matrix();
  MUZ_SO = dipole[2]->to_block_matrix();

  X = block_matrix(nmo,nso); /* just a temporary matrix */

  if((params.ref == 0) || (params.ref == 1)) {

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

    free_block(scf_qt);
  }
  else if(params.ref == 2) {

    MUX_MO_A = block_matrix(nmo,nmo);
    MUY_MO_A = block_matrix(nmo,nmo);
    MUZ_MO_A = block_matrix(nmo,nmo);
    MUX_MO_B = block_matrix(nmo,nmo);
    MUY_MO_B = block_matrix(nmo,nmo);
    MUZ_MO_B = block_matrix(nmo,nmo);

    C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt_A[0][0]),nmo,&(MUX_SO[0][0]),nso,
        0,&(X[0][0]),nso);
    C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt_A[0][0]),nmo,
        0,&(MUX_MO_A[0][0]),nmo);

    C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt_B[0][0]),nmo,&(MUX_SO[0][0]),nso,
        0,&(X[0][0]),nso);
    C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt_B[0][0]),nmo,
        0,&(MUX_MO_B[0][0]),nmo);

    C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt_A[0][0]),nmo,&(MUY_SO[0][0]),nso,
        0,&(X[0][0]),nso);
    C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt_A[0][0]),nmo,
        0,&(MUY_MO_A[0][0]),nmo);

    C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt_B[0][0]),nmo,&(MUY_SO[0][0]),nso,
        0,&(X[0][0]),nso);
    C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt_B[0][0]),nmo,
        0,&(MUY_MO_B[0][0]),nmo);

    C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt_A[0][0]),nmo,&(MUZ_SO[0][0]),nso,
        0,&(X[0][0]),nso);
    C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt_A[0][0]),nmo,
        0,&(MUZ_MO_A[0][0]),nmo);

    C_DGEMM('t','n',nmo,nso,nso,1,&(scf_qt_B[0][0]),nmo,&(MUZ_SO[0][0]),nso,
        0,&(X[0][0]),nso);
    C_DGEMM('n','n',nmo,nmo,nso,1,&(X[0][0]),nso,&(scf_qt_B[0][0]),nmo,
        0,&(MUZ_MO_B[0][0]),nmo);

    free_block(scf_qt_A);
    free_block(scf_qt_B);
  }

  free_block(X);
  free_block(MUX_SO);
  free_block(MUY_SO);
  free_block(MUZ_SO);

  outfile->Printf("\n\tOscillator Strength for %d%3s to %d%3s\n",S->root+1,
          moinfo.labels[S->irrep], U->root+1, moinfo.labels[U->irrep]);
  outfile->Printf("\t                              X    \t       Y    \t       Z\n");


  if((params.ref == 0) || (params.ref == 1)) {

    for(i=0; i < nmo; i++)
      for(j=0; j < nmo; j++) {
        lt_x += MUX_MO[i][j] * moinfo.ltd[i][j];
        lt_y += MUY_MO[i][j] * moinfo.ltd[i][j];
        lt_z += MUZ_MO[i][j] * moinfo.ltd[i][j];
      }

    for(i=0; i < nmo; i++)
      for(j=0; j < nmo; j++) {
        rt_x += MUX_MO[i][j] * moinfo.rtd[i][j];
        rt_y += MUY_MO[i][j] * moinfo.rtd[i][j];
        rt_z += MUZ_MO[i][j] * moinfo.rtd[i][j];
      }

  }
  else if(params.ref == 2) {

    for(i=0; i < nmo; i++)
      for(j=0; j < nmo; j++) {
        lt_x += MUX_MO_A[i][j] * moinfo.ltd_a[i][j];
        lt_y += MUY_MO_A[i][j] * moinfo.ltd_a[i][j];
        lt_z += MUZ_MO_A[i][j] * moinfo.ltd_a[i][j];
      }

    for(i=0; i < nmo; i++)
      for(j=0; j < nmo; j++) {
        rt_x += MUX_MO_A[i][j] * moinfo.rtd_a[i][j];
        rt_y += MUY_MO_A[i][j] * moinfo.rtd_a[i][j];
        rt_z += MUZ_MO_A[i][j] * moinfo.rtd_a[i][j];
      }

    for(i=0; i < nmo; i++)
      for(j=0; j < nmo; j++) {
        lt_x += MUX_MO_B[i][j] * moinfo.ltd_b[i][j];
        lt_y += MUY_MO_B[i][j] * moinfo.ltd_b[i][j];
        lt_z += MUZ_MO_B[i][j] * moinfo.ltd_b[i][j];
      }

    for(i=0; i < nmo; i++)
      for(j=0; j < nmo; j++) {
        rt_x += MUX_MO_B[i][j] * moinfo.rtd_b[i][j];
        rt_y += MUY_MO_B[i][j] * moinfo.rtd_b[i][j];
        rt_z += MUZ_MO_B[i][j] * moinfo.rtd_b[i][j];
      }
  }

  ds_x = lt_x * rt_x;
  ds_y = lt_y * rt_y;
  ds_z = lt_z * rt_z;

  /* Use |w2 - w1| for oscillator strengths */
     // We view excitation energies as positive,
     // so we want to substract the lower state's energy from the
     // higher state's.
     // U should be the higher-energy excited state.
  delta_ee = U->cceom_energy - S->cceom_energy;

  f_x = (2*delta_ee*ds_x)/3;
  f_y = (2*delta_ee*ds_y)/3;
  f_z = (2*delta_ee*ds_z)/3;

  f = f_x + f_y + f_z;

  /* Fill in XTD_Params for this Transition */
  xtd_data->root1        = S->root;
  xtd_data->root2        = U->root;
  xtd_data->irrep1       = S->irrep;
  xtd_data->irrep2       = U->irrep;
  xtd_data->cceom_energy = delta_ee;
  xtd_data->OS           = f;

  /* Compute Einstein A,B Coefficients */
  double hartree2Hz = pc_hartree2MHz * (1.0e6);
  double hbar       = pc_h/(pc_twopi);
  /* SI Dipole Strength */
  double ds_si = (ds_x+ds_y+ds_z) * pc_dipmom_au2si * pc_dipmom_au2si;
  /* SI Transition Energy */
  double nu_si = delta_ee * hartree2Hz;
  /* Einstein Coefficients */
  einstein_b = (2.0/3.0) * (pc_pi/pow(hbar,2)) * (1.0/(4.0*pc_pi*pc_e0)) * ds_si;
  einstein_a = 8.0* pc_pi * pc_h * pow((nu_si/pc_c),3) * einstein_b;
  if(einstein_a < 1e-7) einstein_a = 0.0000000;
  xtd_data->einstein_a = einstein_a;
  xtd_data->einstein_b = einstein_b;

  outfile->Printf("\t<0|mu_e|n>              %11.8lf \t %11.8lf \t %11.8lf\n",
          lt_x,lt_y,lt_z);
  outfile->Printf("\t<n|mu_e|0>              %11.8lf \t %11.8lf \t %11.8lf\n",
          rt_x,rt_y,rt_z);
  outfile->Printf("\tDipole Strength         %11.8lf \n",ds_x+ds_y+ds_z);
  outfile->Printf("\tOscillator Strength     %11.8lf \n",f_x+f_y+f_z);
  outfile->Printf("\tEinstein A Coefficient   %11.8e \n",einstein_a);
  outfile->Printf("\tEinstein B Coefficient   %11.8e \n",einstein_b);


  if((params.ref == 0) || (params.ref == 1)) {
    free_block(MUX_MO);
    free_block(MUY_MO);
    free_block(MUZ_MO);
  }
  else if(params.ref == 2) {
    free_block(MUX_MO_A);
    free_block(MUY_MO_A);
    free_block(MUZ_MO_A);
    free_block(MUX_MO_B);
    free_block(MUY_MO_B);
    free_block(MUZ_MO_B);
  }

  return;
}

}} // namespace psi::ccdensity
