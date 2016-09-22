/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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
    \ingroup CCSORT
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

using namespace std;

namespace psi { namespace ccsort {

void local_magnetic(const char *cart, int **domain, int *domain_len,
		    int natom, int *aostart, int *aostop)
{
  int i, j, ij, a, k, max, complete, *boolean, *rank;
  int nso, nmo;
  double **TMP;
  double **L, **Z;
  double **C, **usotao;
  double mag, mag_i, mag_k, mag_k_check, value, *mag_mo, *mag_mo_check, **mag_atom;
  dpdfile2 U;
  psio_address next;

  nso = moinfo.nso;
  nmo = moinfo.nmo;

  /* grab the occupied MOs */
  next = PSIO_ZERO;
  C = block_matrix(moinfo.nso, moinfo.occpi[0]);
  psio_read(PSIF_CC_INFO, "RHF/ROHF Active Occupied Orbitals", (char *) C[0],
	    nso*moinfo.occpi[0]*sizeof(double), next, &next);

  L = block_matrix(nso,moinfo.occpi[0]);
  Z = block_matrix(nso, moinfo.occpi[0]);
  mag_mo = init_array(moinfo.occpi[0]);
  mag_mo_check = init_array(moinfo.occpi[0]);
  mag_atom = block_matrix(moinfo.occpi[0],natom);
  boolean = init_int_array(natom);
  rank = (int *) malloc(natom * sizeof(int *));

  MintsHelper mints(Process::environment.legacy_wavefunction()->basisset(),
                    Process::environment.options, 0);
  vector<SharedMatrix> angmom = mints.so_angular_momentum();

  if (!strcmp(cart, "X")) {
    angmom[0]->scale(-0.5);
    double **TMP = angmom[0]->to_block_matrix();

    C_DGEMM('n', 'n', nso, moinfo.occpi[0], nso, 1, &(TMP[0][0]), nso, &(C[0][0]), moinfo.occpi[0],
	    0, &(L[0][0]), moinfo.occpi[0]);

    global_dpd_->file2_init(&U, PSIF_CC_OEI, 0, 1, 0, "CPHF Ub_X_AI");
    global_dpd_->file2_mat_init(&U);
    global_dpd_->file2_mat_rd(&U);

    C_DGEMM('n', 'n', nso,moinfo.occpi[0],moinfo.virtpi[0], 1, &(moinfo.C[0][0][0]), moinfo.virtpi[0],
	    &(U.matrix[0][0][0]), moinfo.occpi[0], 0, &(Z[0][0]), moinfo.occpi[0]);

    mag = 0.0;
    for(i=0; i < moinfo.occpi[0]; i++) {
      mag_mo[i] = 0.0;
      mag_mo_check[i] = 0.0;
      for(k=0; k < natom; k++) {
	mag_atom[i][k] = 0.0;
	for(a=aostart[k]; a <= aostop[k]; a++) {
	  value = 4.0 * Z[a][i] * L[a][i];
	  mag += value;
	  mag_mo[i] += value;
	  mag_mo_check[i] += fabs(value);
	  mag_atom[i][k] += fabs(value);
	}
      }
    }

    for(i=0; i < moinfo.occpi[0]; i++) {
      /* Rank the atomic contributions to the orbital's magnetizability */
      for(j=0; j < natom; j++) {
	rank[j] = 0;
	boolean[j] = 0;
      }
      for(j=0,max=0; j < natom; j++) /* find the overall maximum */
	if(fabs(mag_atom[i][j]) >= fabs(mag_atom[i][max])) max = j;
      rank[0] = max;
      boolean[max] = 1;
      for(j=1; j < natom; j++) {
	max = 0;
	while(boolean[max]) max++; /* find an unused max */
	for(k=0; k < natom; k++)
	  if((fabs(mag_atom[i][k]) >= fabs(mag_atom[i][max])) && !boolean[k]) max = k;
	rank[j] = max;
	boolean[max] = 1;
      }

      /* Response domains for Alpha_xx */
      mag_i = 0.0;
      complete = 0;
      for(k=0; k < natom; k++) {
	mag_k = 0.0;
	mag_k_check = 0.0;
	for(a=aostart[rank[k]]; a <= aostop[rank[k]]; a++) {
	  value = Z[a][i] * L[a][i];
	  mag_k += value;
	  mag_k_check += fabs(value);
	}
	mag_k *= 4.0;
	mag_k_check *= 4.0;
	mag_i += mag_k_check;
	if (!complete) {
	  if(!domain[i][rank[k]]) {
	    domain[i][rank[k]] = 1;
	    domain_len[i]++;
	  }
	  if(fabs((mag_mo_check[i]-mag_i)/mag_mo_check[i]) <= local.cphf_cutoff) {
	    complete = 1;
	    if(fabs(mag_atom[i][k]-mag_atom[i][k+1]) <= 0.0001) {
	      if(!domain[i][rank[k]]) {
		domain[i][rank[k]] = 1;
		domain_len[i]++;
	      }
	    }
	  }
	}
      }
    }

    global_dpd_->file2_mat_close(&U);
    global_dpd_->file2_close(&U);
  }

  if (!strcmp(cart, "Y")) {
    angmom[1]->scale(-0.5);
    double **TMP = angmom[1]->to_block_matrix();

    C_DGEMM('n', 'n', nso, moinfo.occpi[0], nso, 1, &(TMP[0][0]), nso, &(C[0][0]), moinfo.occpi[0],
	    0, &(L[0][0]), moinfo.occpi[0]);

    global_dpd_->file2_init(&U, PSIF_CC_OEI, 0, 1, 0, "CPHF Ub_Y_AI");
    global_dpd_->file2_mat_init(&U);
    global_dpd_->file2_mat_rd(&U);

    C_DGEMM('n', 'n', nso,moinfo.occpi[0],moinfo.virtpi[0], 1, &(moinfo.C[0][0][0]), moinfo.virtpi[0],
	    &(U.matrix[0][0][0]), moinfo.occpi[0], 0, &(Z[0][0]), moinfo.occpi[0]);

    mag = 0.0;
    for(i=0; i < moinfo.occpi[0]; i++) {
      mag_mo_check[i] = 0.0;
      mag_mo[i] = 0.0;
      for(k=0; k < natom; k++) {
	mag_atom[i][k] = 0.0;
	for(a=aostart[k]; a <= aostop[k]; a++) {
	  value = 4.0 * Z[a][i] * L[a][i];
	  mag += value;
	  mag_mo[i] += value;
	  mag_mo_check[i] += fabs(value);
	  mag_atom[i][k] += fabs(value);
	}
      }
    }

    for(i=0; i < moinfo.occpi[0]; i++) {
      /* Rank the atomic contributions to the orbital's magnetizability */
      for(j=0; j < natom; j++) {
	rank[j] = 0;
	boolean[j] = 0;
      }
      for(j=0,max=0; j < natom; j++) /* find the overall maximum */
	if(fabs(mag_atom[i][j]) >= fabs(mag_atom[i][max])) max = j;
      rank[0] = max; boolean[max] = 1;
      for(j=1; j < natom; j++) {
	max = 0;
	while(boolean[max]) max++; /* find an unused max */
	for(k=0; k < natom; k++)
	  if((fabs(mag_atom[i][k]) >= fabs(mag_atom[i][max])) && !boolean[k]) max = k;
	rank[j] = max;
	boolean[max] = 1;
      }

      /* Response domains for Alpha_yy */
      mag_i = 0.0;
      complete = 0;
      for(k=0; k < natom; k++) {
	mag_k = 0.0;
	mag_k_check = 0.0;
	for(a=aostart[rank[k]]; a <= aostop[rank[k]]; a++) {
	  value = Z[a][i] * L[a][i];
	  mag_k += value;
	  mag_k_check += fabs(value);
	}
	mag_k *= 4.0;
	mag_k_check *= 4.0;
	mag_i += mag_k_check;
	if (!complete) {
	  if(!domain[i][rank[k]]) {
	    domain[i][rank[k]] = 1;
	    domain_len[i]++;
	  }
	  if(fabs((mag_mo_check[i]-mag_i)/mag_mo_check[i]) <= local.cphf_cutoff) {
	    complete = 1;
	    if(fabs(mag_atom[i][k]-mag_atom[i][k+1]) <= 0.0001) {
	      if(!domain[i][rank[k]]) {
		domain[i][rank[k]] = 1;
		domain_len[i]++;
	      }
	    }
	  }
	}
      }
    }

    global_dpd_->file2_mat_close(&U);
    global_dpd_->file2_close(&U);
  }

  if (!strcmp(cart, "Z")) {
    angmom[2]->scale(-0.5);
    double **TMP = angmom[2]->to_block_matrix();

    C_DGEMM('n', 'n', nso, moinfo.occpi[0], nso, 1, &(TMP[0][0]), nso, &(C[0][0]), moinfo.occpi[0],
	    0, &(L[0][0]), moinfo.occpi[0]);

    global_dpd_->file2_init(&U, PSIF_CC_OEI, 0, 1, 0, "CPHF Ub_Z_AI");
    global_dpd_->file2_mat_init(&U);
    global_dpd_->file2_mat_rd(&U);

    C_DGEMM('n', 'n', nso,moinfo.occpi[0],moinfo.virtpi[0], 1, &(moinfo.C[0][0][0]), moinfo.virtpi[0],
	    &(U.matrix[0][0][0]), moinfo.occpi[0], 0, &(Z[0][0]), moinfo.occpi[0]);

    mag = 0.0;
    for(i=0; i < moinfo.occpi[0]; i++) {
      mag_mo[i] = 0.0;
      mag_mo_check[i] = 0.0;
      for(k=0; k < natom; k++) {
	mag_atom[i][k] = 0.0;
	for(a=aostart[k]; a <= aostop[k]; a++) {
	  value = 4.0 * Z[a][i] * L[a][i];
	  mag += value;
	  mag_mo[i] += value;
	  mag_mo_check[i] += fabs(value);
	  mag_atom[i][k] += fabs(value);
	}
      }
    }

    for(i=0; i < moinfo.occpi[0]; i++) {
      /* Rank the atomic contributions to the orbital's magnetizability */
      for(j=0; j < natom; j++) {
	rank[j] = 0;
	boolean[j] = 0;
      }
      for(j=0,max=0; j < natom; j++) /* find the overall maximum */
	if(fabs(mag_atom[i][j]) >= fabs(mag_atom[i][max])) max = j;
      rank[0] = max; boolean[max] = 1;
      for(j=1; j < natom; j++) {
	max = 0;
	while(boolean[max]) max++; /* find an unused max */
	for(k=0; k < natom; k++)
	  if((fabs(mag_atom[i][k]) >= fabs(mag_atom[i][max])) && !boolean[k]) max = k;
	rank[j] = max;
	boolean[max] = 1;
      }

      /* Response domains for Alpha_zz */
      mag_i = 0.0;
      complete = 0;
      for(k=0; k < natom; k++) {
	mag_k = 0.0;
	mag_k_check = 0.0;
	for(a=aostart[rank[k]]; a <= aostop[rank[k]]; a++) {
	  value = Z[a][i] * L[a][i];
	  mag_k += value;
	  mag_k_check += fabs(value);
	}
	mag_k *= 4.0;
	mag_k_check *= 4.0;
	mag_i += mag_k_check;
	if (!complete) {
	  if(!domain[i][rank[k]]) {
	    domain[i][rank[k]] = 1;
	    domain_len[i]++;
	  }
	  if(fabs((mag_mo_check[i]-mag_i)/mag_mo_check[i]) <= local.cphf_cutoff) {
	    complete = 1;
	    if(fabs(mag_atom[i][k]-mag_atom[i][k+1]) <= 0.0001) {
	      if(!domain[i][rank[k]]) {
		domain[i][rank[k]] = 1;
		domain_len[i]++;
	      }
	    }
	  }
	}
      }
    }

    global_dpd_->file2_mat_close(&U);
    global_dpd_->file2_close(&U);
  }

  free_block(L);

  free_block(mag_atom);
  free_block(C);

  free(boolean);
  free(mag_mo_check);
  free(mag_mo);
  free(rank);
}

}} // namespace psi::ccsort
