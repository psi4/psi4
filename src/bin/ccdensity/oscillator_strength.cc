/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
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
#include <physconst.h>

#define IOFF_MAX 32641
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

void oscillator_strength(struct TD_Params *S)
{
  int nmo, nso, nao, noei, stat, i, I, h, j, nirreps, *ioff;
  int *order, *order_A, *order_B, *doccpi, natom, *clsdpi, *openpi, *orbspi;
  double **scf_pitzer, **scf_pitzer_A, **scf_pitzer_B;
  double **scf_qt, **scf_qt_A, **scf_qt_B, **X, **usotao;
  double *zvals, **geom;
  double *mu_x_ints, *mu_y_ints, *mu_z_ints;
  double **MUX_AO, **MUY_AO, **MUZ_AO;
  double **MUX_MO, **MUY_MO, **MUZ_MO;
  double **MUX_SO, **MUY_SO, **MUZ_SO;
  double **MUX_MO_A, **MUY_MO_A, **MUZ_MO_A;
  double **MUX_MO_B, **MUY_MO_B, **MUZ_MO_B;
  double lt_x, lt_y, lt_z;
  double rt_x, rt_y, rt_z;
  double ds_x, ds_y, ds_z;
  double f_x, f_y, f_z;
  double f;

  chkpt_init(PSIO_OPEN_OLD);
  if ((params.ref == 0) || (params.ref == 1))
    scf_pitzer = chkpt_rd_scf();
  else if(params.ref == 2) {
    scf_pitzer_A = chkpt_rd_alpha_scf();
    scf_pitzer_B = chkpt_rd_beta_scf();
  }

  nso = chkpt_rd_nso();
  nao = chkpt_rd_nao();
  nmo = chkpt_rd_nmo();
  usotao = chkpt_rd_usotao();
  clsdpi = chkpt_rd_clsdpi();
  openpi = chkpt_rd_openpi();
  orbspi = chkpt_rd_orbspi();
  nirreps = chkpt_rd_nirreps();
  natom = chkpt_rd_natom();
  zvals = chkpt_rd_zvals();
  geom = chkpt_rd_geom();
  chkpt_close();

  lt_x = lt_y = lt_z = 0.0;
  rt_x = rt_y = rt_z = 0.0;
  ds_x = ds_y = ds_z = 0.0;
  f_x = f_y = f_z = 0.0;
  f = 0;

  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;

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

  MUX_SO = block_matrix(nso,nso);
  MUY_SO = block_matrix(nso,nso);
  MUZ_SO = block_matrix(nso,nso);

  /*** Transform the AO dipole integrals to the SO basis ***/
  X = block_matrix(nso,nao); /* just a temporary matrix */

  C_DGEMM('n','n',nso,nao,nao,1,&(usotao[0][0]),nao,&(MUX_AO[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','t',nso,nso,nao,1,&(X[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(MUX_SO[0][0]),nso);

  C_DGEMM('n','n',nso,nao,nao,1,&(usotao[0][0]),nao,&(MUY_AO[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','t',nso,nso,nao,1,&(X[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(MUY_SO[0][0]),nso);

  C_DGEMM('n','n',nso,nao,nao,1,&(usotao[0][0]),nao,&(MUZ_AO[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','t',nso,nso,nao,1,&(X[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(MUZ_SO[0][0]),nso);

  free(mu_x_ints); free(mu_y_ints); free(mu_z_ints);
  free_block(X);
  free_block(usotao);
  free_block(MUX_AO);
  free_block(MUY_AO);
  free_block(MUZ_AO);

  /*** Transform the SO dipole integrals to the MO basis ***/

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
  else if((params.ref == 2)) {

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

/*   fprintf(outfile, "MUX_MOs\n"); */
/*   print_mat(MUX_MO, nmo, nmo, outfile); */
/*   fprintf(outfile, "MUY_MOs\n"); */
/*   print_mat(MUY_MO, nmo, nmo, outfile); */
/*   fprintf(outfile, "MUZ_MOs\n"); */
/*   print_mat(MUZ_MO, nmo, nmo, outfile); */

  fprintf(outfile,"\n\tOscillator Strength for %d%3s\n",S->root+1,
          moinfo.labels[S->irrep]);
  fprintf(outfile,"\t                              X    \t       Y    \t       Z\n");

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
  
  f_x = (2*S->cceom_energy*ds_x)/3;
  f_y = (2*S->cceom_energy*ds_y)/3;
  f_z = (2*S->cceom_energy*ds_z)/3;
  
  f = f_x + f_y + f_z;
  S->OS = f;
  
  fprintf(outfile,"\t<0|mu_e|n>              %11.8lf \t %11.8lf \t %11.8lf\n",
          lt_x,lt_y,lt_z);
  fprintf(outfile,"\t<n|mu_e|0>              %11.8lf \t %11.8lf \t %11.8lf\n",
          rt_x,rt_y,rt_z);
  fprintf(outfile,"\tDipole Strength         %11.8lf \n",ds_x+ds_y+ds_z);
  fprintf(outfile,"\tOscillator Strength     %11.8lf \n",f_x+f_y+f_z);
  fflush(outfile);

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
