/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void transL(double sign)
{
  int nao, nso, nmo, noei_ao, nt;
  double **scfp;
  double **scf;
  double **usotao;
  double **TMP, **X, *scratch;
  int stat, i, j, ij, p, q;
  int I, h;
  int errcod;
  int *doccpi;
  int *order;
  double **LX_MO, **LY_MO, **LZ_MO;

  chkpt_init(PSIO_OPEN_OLD);
  if(params.ref == 0 || params.ref == 1) {
    scfp = chkpt_rd_scf();
  }
  usotao = chkpt_rd_usotao();
  nao = chkpt_rd_nao();
  nso = chkpt_rd_nso();
  chkpt_close();

  nmo = moinfo.nmo;
  noei_ao = nao * (nao+1)/2;

  /* doccpi array must include frozen orbitals for reorder_qt() */

  doccpi = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++)
    doccpi[h] = moinfo.frdocc[h] + moinfo.clsdpi[h];

  /*** Get the Pitzer -> QT reordering array ***/
  order = init_int_array(nmo);

  reorder_qt(doccpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc,
             order, moinfo.orbspi, moinfo.nirreps);

  /*** Reorder the SCF eigenvectors to QT ordering */
  
  scf = block_matrix(nmo, nmo);
  for(i=0; i < nmo; i++) {
    I = order[i];  
    for(j=0; j < nmo; j++) scf[j][I] = scfp[j][i];
  } 
  free(order);
  free(doccpi);
  free_block(scfp);

  scratch = init_array(noei_ao);
  TMP = block_matrix(nao, nao);
  X = block_matrix(nao, nao);

  /* NB: (a|m|b) = -1/2 (a|L|b) */
  /* NB: The angular momentum integrals are antisymmetric! */
  iwl_rdone(PSIF_OEI, PSIF_AO_LX, scratch, noei_ao, 0, 0, outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++) {
      TMP[i][j] = -0.5 * sign * scratch[ij];
      TMP[j][i] = +0.5 * sign * scratch[ij];
    }

  LX_MO = block_matrix(nmo, nmo);

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  zero_mat(TMP, nao, nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	  0,&(TMP[0][0]),nao);

  if(params.ref == 0 || params.ref == 1) {
    C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(scf[0][0]),nmo,
	    0,&(X[0][0]),nao);
    C_DGEMM('t','n',nmo,nmo,nso,1,&(scf[0][0]),nmo,&(X[0][0]),nao,
	    0,&(LX_MO[0][0]),nmo);
  }
  else if(params.ref == 2) {
  }

  zero_arr(scratch, noei_ao);
  zero_mat(TMP, nao, nao);
  zero_mat(X, nao, nao);

  iwl_rdone(PSIF_OEI, PSIF_AO_LY, scratch, noei_ao, 0, 0, outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++) {
      TMP[i][j] = -0.5 * sign * scratch[ij];
      TMP[j][i] = +0.5 * sign * scratch[ij];
    }

  LY_MO = block_matrix(nmo, nmo);

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  zero_mat(TMP, nao, nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	  0,&(TMP[0][0]),nao);

  if(params.ref == 0 || params.ref == 1) {
    C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(scf[0][0]),nmo,
	    0,&(X[0][0]),nao);
    C_DGEMM('t','n',nmo,nmo,nso,1,&(scf[0][0]),nmo,&(X[0][0]),nao,
	    0,&(LY_MO[0][0]),nmo);
  }
  else if(params.ref == 2) {
  }

  zero_arr(scratch, noei_ao);
  zero_mat(TMP, nao, nao);
  zero_mat(X, nao, nao);

  iwl_rdone(PSIF_OEI, PSIF_AO_LZ, scratch, noei_ao, 0, 0, outfile);
  for(i=0,ij=0; i < nao; i++)
    for(j=0; j <= i; j++,ij++) {
      TMP[i][j] = -0.5 * sign * scratch[ij];
      TMP[j][i] = +0.5 * sign * scratch[ij];
    }

  LZ_MO = block_matrix(nmo, nmo);

  C_DGEMM('n','t',nao,nso,nao,1,&(TMP[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(X[0][0]),nao);
  zero_mat(TMP, nao, nao);
  C_DGEMM('n','n',nso,nso,nao,1,&(usotao[0][0]),nao,&(X[0][0]),nao,
	  0,&(TMP[0][0]),nao);

  if(params.ref == 0 || params.ref == 1) {
    C_DGEMM('n','n',nso,nmo,nso,1,&(TMP[0][0]),nao,&(scf[0][0]),nmo,
	    0,&(X[0][0]),nao);
    C_DGEMM('t','n',nmo,nmo,nso,1,&(scf[0][0]),nmo,&(X[0][0]),nao,
	    0,&(LZ_MO[0][0]),nmo);
  }
  else if(params.ref == 2) {
  }

  moinfo.L = (double ***) malloc(3 * sizeof(double **));
  moinfo.L[0] = block_matrix(nmo, nmo);
  moinfo.L[1] = block_matrix(nmo, nmo);
  moinfo.L[2] = block_matrix(nmo, nmo);

  for(i=0; i<nmo; i++) 
    for(j=0; j<nmo; j++) {
      moinfo.L[0][i][j] = LX_MO[i][j];
      moinfo.L[1][i][j] = LY_MO[i][j];
      moinfo.L[2][i][j] = LZ_MO[i][j];
    }

  free(scratch);
  free_block(X);
  free_block(TMP);

  if(params.ref == 0 || params.ref == 1) {
    free_block(scf);
  }
  else if(params.ref == 2) {
  }

  free_block(usotao);
  free_block(LX_MO);
  free_block(LY_MO);
  free_block(LZ_MO);
}

}} // namespace psi::ccdensity
