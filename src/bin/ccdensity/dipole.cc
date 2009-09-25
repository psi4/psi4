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

void dipole(void)
{
  int nmo, nso, nao, noei, stat, i, I, h, j, nirreps;
  int *order, *doccpi, *ioff, natom, *clsdpi, *openpi, *orbspi;
  double **scf_pitzer, **scfA, **scfB;
  double **scf_qt, **X, **usotao;
  double *zvals, **geom;
  double *mu_x_ints, *mu_y_ints, *mu_z_ints;
  double **MUX_AO, **MUY_AO, **MUZ_AO;
  double **MUX_SO, **MUY_SO, **MUZ_SO;
  double **MUX_MO, **MUY_MO, **MUZ_MO;
  double mu_x, mu_y, mu_z, mu;
  double mu_x_n, mu_y_n, mu_z_n;
  double mu_x_ref, mu_y_ref, mu_z_ref;
  double mu_x_tot, mu_y_tot, mu_z_tot;

  if (params.ref == 2) return; /* doesn't do UHF yet */

  chkpt_init(PSIO_OPEN_OLD);
  if ((params.ref == 0) || (params.ref == 1))
    scf_pitzer = chkpt_rd_scf();
  else if(params.ref == 2) {
    scfA = chkpt_rd_alpha_scf();
    scfB = chkpt_rd_beta_scf();
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

  mu_x_n = mu_y_n = mu_z_n = 0.0;
  mu_x_ref = mu_y_ref = mu_z_ref = 0.0;
  mu_x = mu_y = mu_z = 0.0;

  /* compute nuclear contribution to dipole moment */
  for(i=0;i<natom;i++) {
    mu_x_n += zvals[i]*geom[i][0];
    mu_y_n += zvals[i]*geom[i][1];
    mu_z_n += zvals[i]*geom[i][2]; 
  }

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

/*   fprintf(outfile, "SCF MO's\n"); */
/*   print_mat(scf_qt, nso, nmo, outfile); */

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

  MUX_MO = block_matrix(nmo,nmo);
  MUY_MO = block_matrix(nmo,nmo);
  MUZ_MO = block_matrix(nmo,nmo);


  /*** Transform the AO dipole integrals to the SO basis ***/
  X = block_matrix(nso,nao); /* just a temporary matrix */

  C_DGEMM('n','n',nso,nao,nao,1.0,&(usotao[0][0]),nao,&(MUX_AO[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','t',nso,nso,nao,1.0,&(X[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(MUX_MO[0][0]),nso);

  C_DGEMM('n','n',nso,nao,nao,1.0,&(usotao[0][0]),nao,&(MUY_AO[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','t',nso,nso,nao,1.0,&(X[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(MUY_MO[0][0]),nso);

  C_DGEMM('n','n',nso,nao,nao,1.0,&(usotao[0][0]),nao,&(MUZ_AO[0][0]),nao,
	  0,&(X[0][0]),nao);
  C_DGEMM('n','t',nso,nso,nao,1.0,&(X[0][0]),nao,&(usotao[0][0]),nao,
	  0,&(MUZ_MO[0][0]),nso);

  free(mu_x_ints); free(mu_y_ints); free(mu_z_ints);
  free_block(X);
  free_block(usotao);
  free_block(MUX_AO);
  free_block(MUY_AO);
  free_block(MUZ_AO);

  /*** Transform the SO dipole integrals to the MO basis ***/

  X = block_matrix(nmo,nso); /* just a temporary matrix */

  C_DGEMM('t','n',nmo,nso,nso,1.0,&(scf_qt[0][0]),nmo,&(MUX_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1.0,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(MUX_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nso,nso,1.0,&(scf_qt[0][0]),nmo,&(MUY_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1.0,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(MUY_MO[0][0]),nmo);

  C_DGEMM('t','n',nmo,nso,nso,1.0,&(scf_qt[0][0]),nmo,&(MUZ_SO[0][0]),nso,
	  0,&(X[0][0]),nso);
  C_DGEMM('n','n',nmo,nmo,nso,1.0,&(X[0][0]),nso,&(scf_qt[0][0]),nmo,
	  0,&(MUZ_MO[0][0]),nmo);


  free_block(scf_qt);
  free_block(MUX_SO); free_block(MUY_SO); free_block(MUZ_SO);
  free_block(X);

/*   fprintf(outfile, "MUX_MOs\n"); */
/*   print_mat(MUX_MO, nmo, nmo, outfile); */
/*   fprintf(outfile, "MUY_MOs\n"); */
/*   print_mat(MUY_MO, nmo, nmo, outfile); */
/*   fprintf(outfile, "MUZ_MOs\n"); */
/*   print_mat(MUZ_MO, nmo, nmo, outfile); */

  fprintf(outfile, "\tOne-particle density:\n");
  print_mat(moinfo.opdm,nmo,nmo,outfile);

  /*** Contract the correlated dipole moment ***/

  for(i=0; i < nmo; i++)
    for(j=0; j < nmo; j++) {
      mu_x += MUX_MO[i][j] * moinfo.opdm[i][j];
      mu_y += MUY_MO[i][j] * moinfo.opdm[i][j];
      mu_z += MUZ_MO[i][j] * moinfo.opdm[i][j];
    }

  for(i=0; i < moinfo.nclsd; i++) {
    mu_x_ref += MUX_MO[i][i] * 2.0 ;
    mu_y_ref += MUY_MO[i][i] * 2.0 ;
    mu_z_ref += MUZ_MO[i][i] * 2.0 ;
  }

  free(ioff);
  free_block(MUX_MO); free_block(MUY_MO); free_block(MUZ_MO);

  mu_x_tot = mu_x + mu_x_n + mu_x_ref ;
  mu_y_tot = mu_y + mu_y_n + mu_y_ref ;
  mu_z_tot = mu_z + mu_z_n + mu_z_ref ;

  fprintf(outfile,"\n\tReference part of electric dipole moment:\n");
  fprintf(outfile,"\t----------------------------------------------\n");
  fprintf(outfile,"\tmu(X) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_x_ref* _dipmom_au2debye, mu_x_ref* _dipmom_au2si, mu_x_ref);
  fprintf(outfile,"\tmu(Y) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_y_ref* _dipmom_au2debye, mu_y_ref* _dipmom_au2si, mu_y_ref);
  fprintf(outfile,"\tmu(Z) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_z_ref* _dipmom_au2debye, mu_z_ref* _dipmom_au2si, mu_z_ref);
  mu = sqrt(mu_x_ref*mu_x_ref+ mu_y_ref*mu_y_ref+ mu_z_ref*mu_z_ref);
  fprintf(outfile,"\t|mu|  = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu * _dipmom_au2debye, mu * _dipmom_au2si, mu);

  fprintf(outfile,"\n\tCorrelation part of electric dipole moment:\n");
  fprintf(outfile,"\t----------------------------------------------\n");
  fprintf(outfile,"\tmu(X) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_x* _dipmom_au2debye, mu_x* _dipmom_au2si, mu_x);
  fprintf(outfile,"\tmu(Y) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_y* _dipmom_au2debye, mu_y* _dipmom_au2si, mu_y);
  fprintf(outfile,"\tmu(Z) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_z* _dipmom_au2debye, mu_z* _dipmom_au2si, mu_z);
  mu = sqrt(mu_x*mu_x+ mu_y*mu_y+ mu_z*mu_z);
  fprintf(outfile,"\t|mu|  = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu * _dipmom_au2debye, mu * _dipmom_au2si, mu);

  fprintf(outfile,"\n\tNuclear part of electric dipole moment:\n");
  fprintf(outfile,"\t----------------------------------------------\n");
  fprintf(outfile,"\tmu(X) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_x_n * _dipmom_au2debye, mu_x_n * _dipmom_au2si, mu_x_n);
  fprintf(outfile,"\tmu(Y) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_y_n * _dipmom_au2debye, mu_y_n * _dipmom_au2si, mu_y_n);
  fprintf(outfile,"\tmu(Z) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_z_n * _dipmom_au2debye, mu_z_n * _dipmom_au2si, mu_z_n);
  mu = sqrt(mu_x_n*mu_x_n + mu_y_n*mu_y_n + mu_z_n*mu_z_n);
  fprintf(outfile,"\t|mu|  = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu * _dipmom_au2debye, mu * _dipmom_au2si, mu);

  fprintf(outfile,"\n\tTotal electric dipole moment:\n");
  fprintf(outfile,"\t----------------------------------------------\n");
  fprintf(outfile,"\tmu(X) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_x_tot * _dipmom_au2debye, mu_x_tot * _dipmom_au2si, mu_x_tot);
  fprintf(outfile,"\tmu(Y) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_y_tot * _dipmom_au2debye, mu_y_tot * _dipmom_au2si, mu_y_tot);
  fprintf(outfile,"\tmu(Z) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_z_tot * _dipmom_au2debye, mu_z_tot * _dipmom_au2si, mu_z_tot);
  mu = sqrt(mu_x_tot*mu_x_tot + mu_y_tot*mu_y_tot + mu_z_tot*mu_z_tot);
  fprintf(outfile,"\t|mu|  = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu * _dipmom_au2debye, mu * _dipmom_au2si, mu);

  return;
}

}} // namespace psi::ccdensity
