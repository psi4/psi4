/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void local_magnetic(const char *cart, int **domain, int *domain_len,
		    int natom, int *aostart, int *aostop)
{
  int i, j, ij, a, k, max, complete, *boolean, *rank;
  int nao, nso, nmo,noei_ao;
  double **TMP, *scratch, **X;
  double **L, **Z;
  double **C, **usotao;
  double mag, mag_i, mag_k, mag_k_check, value, *mag_mo, *mag_mo_check, **mag_atom;
  dpdfile2 U;
  psio_address next;

  nao = moinfo.nao;
  nso = moinfo.nso;
  nmo = moinfo.nmo;
  noei_ao = nao*(nao+1)/2;

  /* grab the occupied MOs */
  next = PSIO_ZERO;
  C = block_matrix(moinfo.nso, moinfo.occpi[0]);
  psio_read(CC_INFO, "RHF/ROHF Active Occupied Orbitals", (char *) C[0],
	    nso*moinfo.occpi[0]*sizeof(double), next, &next);

  /* grab the usotao matrix */
  chkpt_init(PSIO_OPEN_OLD);
  usotao = chkpt_rd_usotao();
  chkpt_close();

  L = block_matrix(nso,moinfo.occpi[0]);
  TMP = block_matrix(nao,nao);
  X = block_matrix(nao, nao);
  scratch = init_array(noei_ao);
  Z = block_matrix(nso, moinfo.occpi[0]);
  mag_mo = init_array(moinfo.occpi[0]);
  mag_mo_check = init_array(moinfo.occpi[0]);
  mag_atom = block_matrix(moinfo.occpi[0],natom);
  boolean = init_int_array(natom);
  rank = (int *) malloc(natom * sizeof(int *));

  if (!strcmp(cart, "X")) {
    iwl_rdone(PSIF_OEI, PSIF_AO_LX, scratch, noei_ao, 0, 0, outfile);
    for(i=0,ij=0; i<nao; i++)
      for(j=0; j<=i; j++,ij++) {
	TMP[i][j] = -0.5 * scratch[ij];
	TMP[j][i] = 0.5 * scratch[ij];
      }

    C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	    0,&(X[0][0]),nao);
    C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	    0,&(TMP[0][0]),nao);

    C_DGEMM('n', 'n', nso, moinfo.occpi[0], nso, 1, &(TMP[0][0]), nao, &(C[0][0]), moinfo.occpi[0],
	    0, &(L[0][0]), moinfo.occpi[0]);

    dpd_file2_init(&U, CC_OEI, 0, 1, 0, "CPHF Ub_X_AI");
    dpd_file2_mat_init(&U);
    dpd_file2_mat_rd(&U);

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
      /* Rank the atomic contributions to the orbital's magizability */
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

    dpd_file2_mat_close(&U);
    dpd_file2_close(&U);
  }

  if (!strcmp(cart, "Y")) {
    iwl_rdone(PSIF_OEI, PSIF_AO_LY, scratch, noei_ao, 0, 0, outfile);
    for(i=0,ij=0; i<nao; i++)
      for(j=0; j<=i; j++,ij++) {
	TMP[i][j] = -0.5 * scratch[ij];
	TMP[j][i] = 0.5 * scratch[ij];
      }

    C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	    0,&(X[0][0]),nao);
    C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	    0,&(TMP[0][0]),nao);

    C_DGEMM('n', 'n', nso, moinfo.occpi[0], nso, 1, &(TMP[0][0]), nao, &(C[0][0]), moinfo.occpi[0],
	    0, &(L[0][0]), moinfo.occpi[0]);

    dpd_file2_init(&U, CC_OEI, 0, 1, 0, "CPHF Ub_Y_AI");
    dpd_file2_mat_init(&U);
    dpd_file2_mat_rd(&U);

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
      /* Rank the atomic contributions to the orbital's magizability */
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

    dpd_file2_mat_close(&U);
    dpd_file2_close(&U);
  }

  if (!strcmp(cart, "Z")) {
    iwl_rdone(PSIF_OEI, PSIF_AO_LZ, scratch, noei_ao, 0, 0, outfile);
    for(i=0,ij=0; i<nao; i++)
      for(j=0; j<=i; j++,ij++) {
	TMP[i][j] = -0.5 * scratch[ij];
	TMP[j][i] = 0.5 * scratch[ij];
      }

    C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	    0,&(X[0][0]),nao);
    C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	    0,&(TMP[0][0]),nao);

    C_DGEMM('n', 'n', nso, moinfo.occpi[0], nso, 1, &(TMP[0][0]), nao, &(C[0][0]), moinfo.occpi[0],
	    0, &(L[0][0]), moinfo.occpi[0]);

    dpd_file2_init(&U, CC_OEI, 0, 1, 0, "CPHF Ub_Z_AI");
    dpd_file2_mat_init(&U);
    dpd_file2_mat_rd(&U);

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
      /* Rank the atomic contributions to the orbital's magizability */
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

    dpd_file2_mat_close(&U);
    dpd_file2_close(&U);
  }

  free_block(TMP);
  free_block(L);

  free_block(mag_atom);
  free_block(C);

  free(boolean);
  free(mag_mo_check);
  free(mag_mo);
  free(rank);
}

}} // namespace psi::ccsort
