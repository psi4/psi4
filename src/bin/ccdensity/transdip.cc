/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

#define IOFF_MAX 32641
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void transdip(void)
{
  int nmo, nso, nao, noei, stat, i, I, h, j;
  int *order, *doccpi, *ioff;
  double **scf_pitzer, **scfA, **scfB;
  double **scf_qt, **X, **usotao;
  double *zvals, **geom;
  double *mu_x_ints, *mu_y_ints, *mu_z_ints;
  double **MUX_AO, **MUY_AO, **MUZ_AO;
  double **MUX_MO, **MUY_MO, **MUZ_MO;

  chkpt_init(PSIO_OPEN_OLD);
  if ((params.ref == 0) || (params.ref == 1))
    scf_pitzer = chkpt_rd_scf();
  else if(params.ref == 2) {
    scfA = chkpt_rd_alpha_scf();
    scfB = chkpt_rd_beta_scf();
  }

  nao = chkpt_rd_nao();
  nso = chkpt_rd_nso();
  usotao = chkpt_rd_usotao();
  chkpt_close();

  nmo = moinfo.nmo;

  /*** Build ioff ***/
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;

  /*** Get the Pitzer -> QT reordering array ***/
  order = init_int_array(nmo);

  /* doccpi array must include frozen orbitals for reorder_qt() */
  doccpi = init_int_array(moinfo.nirreps);
  for(h=0; h < moinfo.nirreps; h++) 
    doccpi[h] = moinfo.frdocc[h] + moinfo.clsdpi[h];

  reorder_qt(doccpi, moinfo.openpi, moinfo.frdocc, moinfo.fruocc, 
             order, moinfo.orbspi, moinfo.nirreps);

  /*** Reorder the SCF eigenvectors to QT ordering */

  scf_qt = block_matrix(nmo, nmo);
  for(i=0; i < nmo; i++) {
    I = order[i];  /* Pitzer --> QT */
    for(j=0; j < nmo; j++) scf_qt[j][I] = scf_pitzer[j][i];
  }
  free(order);
  free(doccpi);
  free_block(scf_pitzer);

  /*** Read in dipole moment integrals in the AO basis ***/
  noei = nao * (nao + 1)/2;

  mu_x_ints = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MX,mu_x_ints,noei,0,0,outfile);
  mu_y_ints = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MY,mu_y_ints,noei,0,0,outfile);
  mu_z_ints = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MZ,mu_z_ints,noei,0,0,outfile);

  MUX_AO = block_matrix(nao,nao);
  MUY_AO = block_matrix(nao,nao);
  MUZ_AO = block_matrix(nao,nao);

  for(i=0; i < nao; i++)
    for(j=0; j < nao; j++) {
      MUX_AO[i][j] = mu_x_ints[INDEX(i,j)];
      MUY_AO[i][j] = mu_y_ints[INDEX(i,j)];
      MUZ_AO[i][j] = mu_z_ints[INDEX(i,j)];
    }

/*   fprintf(outfile, "MUX_AOs\n"); */
/*   print_mat(MUX_AO, nao, nao, outfile); */
/*   fprintf(outfile, "MUY_AOs\n"); */
/*   print_mat(MUY_AO, nao, nao, outfile); */
/*   fprintf(outfile, "MUZ_AOs\n"); */
/*   print_mat(MUZ_AO, nao, nao, outfile); */

  MUX_MO = block_matrix(nso,nso);
  MUY_MO = block_matrix(nso,nso);
  MUZ_MO = block_matrix(nso,nso);

  /*** Transform the AO dipole integrals to the SO basis ***/
  X = block_matrix(nso,nao); /* just a temporary matrix */

  C_DGEMM('n','n',nso,nao,nao,1,&(usotao[0][0]),nao,&(MUX_AO[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','t',nso,nso,nao,1,&(X[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(MUX_MO[0][0]),nso);

  C_DGEMM('n','n',nso,nao,nao,1,&(usotao[0][0]),nao,&(MUY_AO[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','t',nso,nso,nao,1,&(X[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(MUY_MO[0][0]),nso);

  C_DGEMM('n','n',nso,nao,nao,1,&(usotao[0][0]),nao,&(MUZ_AO[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','t',nso,nso,nao,1,&(X[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(MUZ_MO[0][0]),nso);

  free(mu_x_ints); 
  free(mu_y_ints); 
  free(mu_z_ints);
  free_block(X);
  free_block(usotao);
  free_block(MUX_AO);
  free_block(MUY_AO);
  free_block(MUZ_AO);

  /*** Transform the SO dipole integrals to the MO basis ***/

  X = block_matrix(nmo,nmo); /* just a temporary matrix */

  C_DGEMM('t','n',nmo,nmo,nmo,1,&(scf_qt[0][0]),nmo,&(MUX_MO[0][0]),nmo,
	  0,&(X[0][0]),nmo);
  C_DGEMM('n','n',nmo,nmo,nmo,1,&(X[0][0]),nmo,&(scf_qt[0][0]),nmo,
	  0,&(MUX_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nmo,nmo,1,&(scf_qt[0][0]),nmo,&(MUY_MO[0][0]),nmo,
	  0,&(X[0][0]),nmo);
  C_DGEMM('n','n',nmo,nmo,nmo,1,&(X[0][0]),nmo,&(scf_qt[0][0]),nmo,
	  0,&(MUY_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nmo,nmo,1,&(scf_qt[0][0]),nmo,&(MUZ_MO[0][0]),nmo,
	  0,&(X[0][0]),nmo);
  C_DGEMM('n','n',nmo,nmo,nmo,1,&(X[0][0]),nmo,&(scf_qt[0][0]),nmo,
	  0,&(MUZ_MO[0][0]),nmo);

  free_block(scf_qt);
  free_block(X);

  moinfo.dip = (double ***) malloc(3 * sizeof(double **));
  moinfo.dip[0] = block_matrix(nao, nao);
  moinfo.dip[1] = block_matrix(nao, nao);
  moinfo.dip[2] = block_matrix(nao, nao);

  for(i=0; i<nmo; i++) 
    for(j=0; j<nmo; j++) {
      moinfo.dip[0][i][j] = MUX_MO[i][j];
      moinfo.dip[1][i][j] = MUY_MO[i][j];
      moinfo.dip[2][i][j] = MUZ_MO[i][j];
    }
 
  free(ioff);
  free_block(MUX_MO); 
  free_block(MUY_MO); 
  free_block(MUZ_MO);

  return;
}

}} // namespace psi::ccdensity
