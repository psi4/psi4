/*! \file ccinfo.cc
    \ingroup (CINTS)
    \brief Enter brief description of file here 
*/
#include<cstdio>
#include<stdlib.h>
#include<cmath>
#include<libciomr/libciomr.h>
#include<libpsio/psio.h>
#include<libdpd/dpd.h>
#include<libchkpt/chkpt.h>
#include<libint/libint.h>
#include<libqt/qt.h>
#include <psifiles.h>
#include <ccfiles.h>

#include"moinfo.h"
#include"moinfo_corr.h"

#include"defines.h"
#define EXTERN
#include"global.h"
#include <stdexcept>
#include"small_fns.h"

namespace psi { namespace CINTS {

static int *cachefiles, **cachelist;

void init_ccinfo()
{
  int h, ij, ab, row, col, row_offset, col_offset;
  int nactive, nocc, nvirt, nao;
  int *occpi, *occ_sym;
  int *virtpi, *vir_sym;
  double **A;
  double **T2_MO, **T2_AO, **T2;
  dpdbuf4 tau;

  nao = BasisSet.num_ao;

  /* grab some basic data from CC_INFO */
  psio_open(CC_INFO, PSIO_OPEN_OLD);
  psio_read_entry(CC_INFO, "No. of Active Orbitals", (char *) &(nactive),
                  sizeof(int));
  occpi = init_int_array(Symmetry.nirreps);
  virtpi = init_int_array(Symmetry.nirreps);
  occ_sym = init_int_array(nactive);
  vir_sym = init_int_array(nactive);
  psio_read_entry(CC_INFO, "Active Occ Orbs Per Irrep",
                  (char *) occpi, sizeof(int)*Symmetry.nirreps);
  psio_read_entry(CC_INFO, "Active Virt Orbs Per Irrep",
                  (char *) virtpi, sizeof(int)*Symmetry.nirreps);
  psio_read_entry(CC_INFO, "Active Occ Orb Symmetry",
                  (char *) occ_sym, sizeof(int)*Symmetry.nirreps);
  psio_read_entry(CC_INFO, "Active Virt Orb Symmetry",
                  (char *) vir_sym, sizeof(int)*Symmetry.nirreps);
  psio_close(CC_INFO, 1);

  nocc = 0; nvirt = 0;
  for(h=0; h < Symmetry.nirreps; h++) {
    nocc += occpi[h];  nvirt += virtpi[h];
  }

  cachefiles = init_int_array(PSIO_MAXUNIT);
  /* assuming no caching for DPD here */
  cachelist = init_int_matrix(12,12);

  init_moinfo();
  init_moinfo_corr();

  /* open the T-amplitude file for r/w */
  psio_open(CC_TAMPS, PSIO_OPEN_OLD);

  /*--- Initialize DPD library ---*/
  dpd_init(1, Symmetry.nirreps, UserOptions.memory, 0, cachefiles, cachelist, NULL, 
	   2, occpi, occ_sym, virtpi, vir_sym);

  /* Grab the MO-basis T2's provided by ccenergy */
  T2_MO = block_matrix(nocc*nocc,nvirt*nvirt);
  dpd_buf4_init(&tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  row_offset = col_offset = 0;
  for(h=0,row=0,col=0; h < Symmetry.nirreps; h++) {
    dpd_buf4_mat_irrep_init(&tau, h);
    dpd_buf4_mat_irrep_rd(&tau, h);

    if(h) { 
      row_offset += tau.params->rowtot[h-1];  
      col_offset += tau.params->coltot[h-1];
    }

    for(row=0; row < tau.params->rowtot[h]; row++) {
      for(col=0; col < tau.params->coltot[h]; col++) {
        T2_MO[row_offset+row][col_offset+col] = tau.matrix[h][row][col];
      }
    }

  dpd_buf4_mat_irrep_close(&tau, h);
  }
  dpd_buf4_close(&tau);

  /* half-transform the T2's */
  A = block_matrix(nvirt, nao);
  T2_AO = block_matrix(nocc*nocc, nao*nao);
  for(ij=0; ij < nocc*nocc; ij++) {


    C_DGEMM('n','n',nvirt,nao,nvirt,1.0,&(T2_MO[ij][0]),nvirt,
            &(MOInfo.scf_evec_uocc[0][0][0]),nao,0.0,A[0],nao);

    C_DGEMM('t','n',nao,nao,nvirt,1.0,&(MOInfo.scf_evec_uocc[0][0][0]),nao,
            A[0],nao,0.0,&(T2_AO[ij][0]),nao);
  }
  free_block(A);
  free_block(T2_MO);

  /* transpose T2 to AO*AO x MO*MO */
  T2 = block_matrix(nao*nao,nocc*nocc);
  for(ij=0; ij < nocc*nocc; ij++)
    for(ab=0; ab < nao*nao; ab++) 
      T2[ab][ij] = T2_AO[ij][ab];
  free_block(T2_AO);

  CCInfo.T2_s = T2;
  CCInfo.T2_t = block_matrix(nao*nao,nocc*nocc);
  CCInfo.nocc = nocc;
  CCInfo.nvirt = nvirt;

  free(occpi);  free(occ_sym);
  free(virtpi); free(vir_sym);
      
  return;
}


void cleanup_ccinfo()
{
  int h, ij, ab, row, col, row_offset, col_offset;
  int nocc, nvirt, nao;
  double **T2, **A, **T2_MO;
  dpdbuf4 tau;

  nocc = CCInfo.nocc;
  nvirt = CCInfo.nvirt;
  nao = BasisSet.num_ao;

  free_block(CCInfo.T2_s);

  /* transpose target T2 back to MO*MO x AO*AO */
  T2 = block_matrix(nocc*nocc, nao*nao); 
  for(ij=0; ij < nocc*nocc; ij++) {
    for(ab=0; ab < nao*nao; ab++) {
       T2[ij][ab] = CCInfo.T2_t[ab][ij];
    }
  }
  free_block(CCInfo.T2_t);

  /* transform T2 back to MO basis */
  A = block_matrix(nao,nvirt);
  T2_MO = block_matrix(nocc*nocc,nvirt*nvirt);
  for(ij=0; ij < nocc*nocc; ij++) {
    C_DGEMM('n','t',nao,nvirt,nao,1.0,T2[ij],nao,
            MOInfo.scf_evec_uocc[0][0],nao,0.0,A[0],nvirt);
    C_DGEMM('n','n',nvirt,nvirt,nao,1.0,MOInfo.scf_evec_uocc[0][0],nao,
            A[0],nvirt,0.0,T2_MO[ij],nvirt);
  }
  free_block(A);
  free_block(T2);

  /* Put T2 back into the DPD buffer */
  dpd_buf4_init(&tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
  row_offset = col_offset = 0;
  for(h=0; h < Symmetry.nirreps; h++) {
    dpd_buf4_mat_irrep_init(&tau, h);
    dpd_buf4_mat_irrep_rd(&tau, h);

    if(h) {
      row_offset += tau.params->rowtot[h-1];
      col_offset += tau.params->coltot[h-1];
    }

    for(row=0; row < tau.params->rowtot[h]; row++) {
      for(col=0; col < tau.params->coltot[h]; col++) {
         tau.matrix[h][row][col] += T2_MO[row+row_offset][col+col_offset];
      }
    }
    dpd_buf4_mat_irrep_wrt(&tau, h);
    dpd_buf4_mat_irrep_close(&tau, h);
  }
  dpd_buf4_close(&tau);

  free_block(T2_MO);

  free_int_matrix(cachelist);
  free(cachefiles);
  dpd_close(1);

  /* close the T-amplitude file */
  psio_close(CC_TAMPS, 1);

  return;
}

};};
