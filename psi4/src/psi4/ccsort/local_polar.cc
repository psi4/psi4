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
#include <cmath>
#include <cstring>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "Local.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

using namespace std;

namespace psi { namespace ccsort {

void local_polar(const char *cart, int **domain, int *domain_len,
		 int natom, int *aostart, int *aostop)
{
  int i, j, ij, a, k, max, complete, *boolean, *rank;
  int nso, nmo;
  double **TMP;
  double **MU, **Z;
  double **C;
  double polar, polar_i, polar_k, polar_k_check, value, *polar_mo, *polar_mo_check, **polar_atom;
  dpdfile2 U;
  psio_address next;
  int stat;

  nso = moinfo.nso;
  nmo = moinfo.nmo;

  /* grab the occupied MOs */
  next = PSIO_ZERO;
  C = block_matrix(moinfo.nso, moinfo.occpi[0]);
  psio_read(PSIF_CC_INFO, "RHF/ROHF Active Occupied Orbitals", (char *) C[0],
	    nso*moinfo.occpi[0]*sizeof(double), next, &next);

  MU = block_matrix(nso,moinfo.occpi[0]);
  Z = block_matrix(nso, moinfo.occpi[0]);
  polar_mo = init_array(moinfo.occpi[0]);
  polar_mo_check = init_array(moinfo.occpi[0]);
  polar_atom = block_matrix(moinfo.occpi[0],natom);
  boolean = init_int_array(natom);
  rank = (int *) malloc(natom * sizeof(int *));

  MintsHelper mints(Process::environment.legacy_wavefunction()->basisset(),
                    Process::environment.options, 0);
  vector<SharedMatrix> dipole = mints.so_dipole();

  if (!strcmp(cart,"X")) {
    double **TMP = dipole[0]->to_block_matrix();

    C_DGEMM('n', 'n', nso, moinfo.occpi[0], nso, 1, &(TMP[0][0]), nso, &(C[0][0]), moinfo.occpi[0],
	    0, &(MU[0][0]), moinfo.occpi[0]);

    global_dpd_->file2_init(&U, PSIF_CC_OEI, 0, 1, 0, "CPHF Uf_X_AI");
    global_dpd_->file2_mat_init(&U);
    global_dpd_->file2_mat_rd(&U);

    C_DGEMM('n', 'n', nso,moinfo.occpi[0],moinfo.virtpi[0], 1, &(moinfo.C[0][0][0]), moinfo.virtpi[0],
	    &(U.matrix[0][0][0]), moinfo.occpi[0], 0, &(Z[0][0]), moinfo.occpi[0]);

    polar = 0.0;
    for(i=0; i < moinfo.occpi[0]; i++) {
      polar_mo[i] = 0.0;
      polar_mo_check[i] = 0.0;
      for(k=0; k < natom; k++) {
	polar_atom[i][k] = 0.0;
	for(a=aostart[k]; a <= aostop[k]; a++) {
	  value = 4.0 * Z[a][i] * MU[a][i];
	  polar += value;
	  polar_mo[i] += value;
	  polar_mo_check[i] += fabs(value);
	  polar_atom[i][k] += fabs(value);
	}
      }
    }

    for(i=0; i < moinfo.occpi[0]; i++) {
      /* Rank the atomic contributions to the orbital's polarizability */
      for(j=0; j < natom; j++) {
	rank[j] = 0;
	boolean[j] = 0;
      }
      for(j=0,max=0; j < natom; j++) /* find the overall maximum */
	if(fabs(polar_atom[i][j]) >= fabs(polar_atom[i][max])) max = j;
      rank[0] = max;
      boolean[max] = 1;
      for(j=1; j < natom; j++) {
	max = 0;
	while(boolean[max]) max++; /* find an unused max */
	for(k=0; k < natom; k++)
	  if((fabs(polar_atom[i][k]) >= fabs(polar_atom[i][max])) && !boolean[k]) max = k;
	rank[j] = max;
	boolean[max] = 1;
      }

      /* Response domains for Alpha_xx */
      polar_i = 0.0;
      complete = 0;
      for(k=0; k < natom; k++) {
	polar_k = 0.0;
	polar_k_check = 0.0;
	for(a=aostart[rank[k]]; a <= aostop[rank[k]]; a++) {
	  value = Z[a][i] * MU[a][i];
	  polar_k += value;
	  polar_k_check += fabs(value);
	}
	polar_k *= 4.0;
	polar_k_check *= 4.0;
	polar_i += polar_k_check;
	if (!complete) {
	  if(!domain[i][rank[k]]) {
	    domain[i][rank[k]] = 1;
	    domain_len[i]++;
	  }
	  if(fabs((polar_mo_check[i]-polar_i)/polar_mo_check[i]) <= local.cphf_cutoff) {
	    complete = 1;
	    if(fabs(polar_atom[i][k]-polar_atom[i][k+1]) <= 0.0001) {
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

  if (!strcmp(cart,"Y")) {
    double **TMP = dipole[1]->to_block_matrix();

    C_DGEMM('n', 'n', nso, moinfo.occpi[0], nso, 1, &(TMP[0][0]), nso, &(C[0][0]), moinfo.occpi[0],
	    0, &(MU[0][0]), moinfo.occpi[0]);

    global_dpd_->file2_init(&U, PSIF_CC_OEI, 0, 1, 0, "CPHF Uf_Y_AI");
    global_dpd_->file2_mat_init(&U);
    global_dpd_->file2_mat_rd(&U);

    C_DGEMM('n', 'n', nso,moinfo.occpi[0],moinfo.virtpi[0], 1, &(moinfo.C[0][0][0]), moinfo.virtpi[0],
	    &(U.matrix[0][0][0]), moinfo.occpi[0], 0, &(Z[0][0]), moinfo.occpi[0]);

    polar = 0.0;
    for(i=0; i < moinfo.occpi[0]; i++) {
      polar_mo_check[i] = 0.0;
      polar_mo[i] = 0.0;
      for(k=0; k < natom; k++) {
	polar_atom[i][k] = 0.0;
	for(a=aostart[k]; a <= aostop[k]; a++) {
	  value = 4.0 * Z[a][i] * MU[a][i];
	  polar += value;
	  polar_mo[i] += value;
	  polar_mo_check[i] += fabs(value);
	  polar_atom[i][k] += fabs(value);
	}
      }
    }

    for(i=0; i < moinfo.occpi[0]; i++) {
      /* Rank the atomic contributions to the orbital's polarizability */
      for(j=0; j < natom; j++) {
	rank[j] = 0;
	boolean[j] = 0;
      }
      for(j=0,max=0; j < natom; j++) /* find the overall maximum */
	if(fabs(polar_atom[i][j]) >= fabs(polar_atom[i][max])) max = j;
      rank[0] = max; boolean[max] = 1;
      for(j=1; j < natom; j++) {
	max = 0;
	while(boolean[max]) max++; /* find an unused max */
	for(k=0; k < natom; k++)
	  if((fabs(polar_atom[i][k]) >= fabs(polar_atom[i][max])) && !boolean[k]) max = k;
	rank[j] = max;
	boolean[max] = 1;
      }

      /* Response domains for Alpha_yy */
      polar_i = 0.0;
      complete = 0;
      for(k=0; k < natom; k++) {
	polar_k = 0.0;
	polar_k_check = 0.0;
	for(a=aostart[rank[k]]; a <= aostop[rank[k]]; a++) {
	  value = Z[a][i] * MU[a][i];
	  polar_k += value;
	  polar_k_check += fabs(value);
	}
	polar_k *= 4.0;
	polar_k_check *= 4.0;
	polar_i += polar_k_check;
	if (!complete) {
	  if(!domain[i][rank[k]]) {
	    domain[i][rank[k]] = 1;
	    domain_len[i]++;
	  }
	  if(fabs((polar_mo_check[i]-polar_i)/polar_mo_check[i]) <= local.cphf_cutoff) {
	    complete = 1;
	    if(fabs(polar_atom[i][k]-polar_atom[i][k+1]) <= 0.0001) {
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

  if (!strcmp(cart,"Z")) {
    double **TMP = dipole[2]->to_block_matrix();

    C_DGEMM('n', 'n', nso, moinfo.occpi[0], nso, 1, &(TMP[0][0]), nso, &(C[0][0]), moinfo.occpi[0],
	    0, &(MU[0][0]), moinfo.occpi[0]);

    global_dpd_->file2_init(&U, PSIF_CC_OEI, 0, 1, 0, "CPHF Uf_Z_AI");
    global_dpd_->file2_mat_init(&U);
    global_dpd_->file2_mat_rd(&U);

    C_DGEMM('n', 'n', nso,moinfo.occpi[0],moinfo.virtpi[0], 1, &(moinfo.C[0][0][0]), moinfo.virtpi[0],
	    &(U.matrix[0][0][0]), moinfo.occpi[0], 0, &(Z[0][0]), moinfo.occpi[0]);

    polar = 0.0;
    for(i=0; i < moinfo.occpi[0]; i++) {
      polar_mo[i] = 0.0;
      polar_mo_check[i] = 0.0;
      for(k=0; k < natom; k++) {
	polar_atom[i][k] = 0.0;
	for(a=aostart[k]; a <= aostop[k]; a++) {
	  value = 4.0 * Z[a][i] * MU[a][i];
	  polar += value;
	  polar_mo[i] += value;
	  polar_mo_check[i] += fabs(value);
	  polar_atom[i][k] += fabs(value);
	}
      }
    }

    for(i=0; i < moinfo.occpi[0]; i++) {
      /* Rank the atomic contributions to the orbital's polarizability */
      for(j=0; j < natom; j++) {
	rank[j] = 0;
	boolean[j] = 0;
      }
      for(j=0,max=0; j < natom; j++) /* find the overall maximum */
	if(fabs(polar_atom[i][j]) >= fabs(polar_atom[i][max])) max = j;
      rank[0] = max; boolean[max] = 1;
      for(j=1; j < natom; j++) {
	max = 0;
	while(boolean[max]) max++; /* find an unused max */
	for(k=0; k < natom; k++)
	  if((fabs(polar_atom[i][k]) >= fabs(polar_atom[i][max])) && !boolean[k]) max = k;
	rank[j] = max;
	boolean[max] = 1;
      }

      /* Response domains for Alpha_zz */
      polar_i = 0.0;
      complete = 0;
      for(k=0; k < natom; k++) {
	polar_k = 0.0;
	polar_k_check = 0.0;
	for(a=aostart[rank[k]]; a <= aostop[rank[k]]; a++) {
	  value = Z[a][i] * MU[a][i];
	  polar_k += value;
	  polar_k_check += fabs(value);
	}
	polar_k *= 4.0;
	polar_k_check *= 4.0;
	polar_i += polar_k_check;
	if (!complete) {
	  if(!domain[i][rank[k]]) {
	    domain[i][rank[k]] = 1;
	    domain_len[i]++;
	  }
	  if(fabs((polar_mo_check[i]-polar_i)/polar_mo_check[i]) <= local.cphf_cutoff) {
	    complete = 1;
	    if(fabs(polar_atom[i][k]-polar_atom[i][k+1]) <= 0.0001) {
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

  free_block(MU);

  free_block(polar_atom);
  free_block(C);

  free(boolean);
  free(polar_mo_check);
  free(polar_mo);
  free(rank);
}

}} // namespace psi::ccsort
