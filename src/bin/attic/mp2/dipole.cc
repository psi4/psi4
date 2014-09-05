/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup MP2
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
#include <physconst.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void dipole(void)
{
  int noei, stat, i, I, h, j, nirreps;
  int *order, *doccpi, natom, *clsdpi, *openpi, *orbspi;
  double **scf_pitzer;
  double **scf_qt, **X, **usotao;
  double *zvals, **geom;
  double *mu_x_ints, *mu_y_ints, *mu_z_ints;
  double **MUX_AO, **MUY_AO, **MUZ_AO;
  double **MUX_MO, **MUY_MO, **MUZ_MO;
  double mu_x, mu_y, mu_z, mu;
  double mu_x_n, mu_y_n, mu_z_n;
  double mu_x_ref, mu_y_ref, mu_z_ref;
  double mu_x_tot, mu_y_tot, mu_z_tot;

  chkpt_init(PSIO_OPEN_OLD);
  scf_pitzer = chkpt_rd_scf();
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

  /*** Get the Pitzer -> QT reordering array ***/
  order = init_int_array(mo.nmo);

  reorder_qt(mo.doccpi, mo.soccpi, mo.fzdoccpi, mo.fzvirtpi, 
             order, mo.mopi, mo.nirreps);

  /*** Reorder the SCF eigenvectors to QT ordering */

  scf_qt = block_matrix(mo.nmo,mo.nmo);
  for(i=0; i < mo.nmo; i++) {
    I = order[i];  /* Pitzer --> QT */
    for(j=0; j < mo.nmo; j++) scf_qt[j][I] = scf_pitzer[j][i];
  }
  free(order);
  free_block(scf_pitzer);

  /*** Read in dipole moment integrals in the AO basis ***/
  noei = mo.nao * (mo.nao + 1)/2;

  mu_x_ints = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MX,mu_x_ints,noei,0,0,outfile);
  mu_y_ints = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MY,mu_y_ints,noei,0,0,outfile);
  mu_z_ints = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_AO_MZ,mu_z_ints,noei,0,0,outfile);

  MUX_AO = block_matrix(mo.nao,mo.nao);
  MUY_AO = block_matrix(mo.nao,mo.nao);
  MUZ_AO = block_matrix(mo.nao,mo.nao);

  for(i=0; i < mo.nao; i++)
    for(j=0; j < mo.nao; j++) {
      MUX_AO[i][j] = mu_x_ints[INDEX(i,j)];
      MUY_AO[i][j] = mu_y_ints[INDEX(i,j)];
      MUZ_AO[i][j] = mu_z_ints[INDEX(i,j)];
    }

  /*
  outfile->Printf( "MUX_AOs\n"); 
  print_mat(MUX_AO, mo.nao, mo.nao, outfile); 
  outfile->Printf( "MUY_AOs\n"); 
  print_mat(MUY_AO, mo.nao, mo.nao, outfile); 
  outfile->Printf( "MUZ_AOs\n"); 
  print_mat(MUZ_AO, mo.nao, mo.nao, outfile); 
  
  */

  MUX_MO = block_matrix(mo.nmo,mo.nmo);
  MUY_MO = block_matrix(mo.nmo,mo.nmo);
  MUZ_MO = block_matrix(mo.nmo,mo.nmo);

  /*** Transform the AO dipole integrals to the SO basis ***/
  X = block_matrix(mo.nso,mo.nao); /* just a temporary matrix */

  C_DGEMM('n','n',mo.nso,mo.nao,mo.nao,1,&(usotao[0][0]),mo.nao,&(MUX_AO[0][0]),mo.nao,
	  0,&(X[0][0]),mo.nao);
  C_DGEMM('n','t',mo.nso,mo.nso,mo.nao,1,&(X[0][0]),mo.nao,&(usotao[0][0]),mo.nao,
	  0,&(MUX_MO[0][0]),mo.nso);

  C_DGEMM('n','n',mo.nso,mo.nao,mo.nao,1,&(usotao[0][0]),mo.nao,&(MUY_AO[0][0]),mo.nao,
	  0,&(X[0][0]),mo.nao);
  C_DGEMM('n','t',mo.nso,mo.nso,mo.nao,1,&(X[0][0]),mo.nao,&(usotao[0][0]),mo.nao,
	  0,&(MUY_MO[0][0]),mo.nso);

  C_DGEMM('n','n',mo.nso,mo.nao,mo.nao,1,&(usotao[0][0]),mo.nao,&(MUZ_AO[0][0]),mo.nao,
	  0,&(X[0][0]),mo.nao);
  C_DGEMM('n','t',mo.nso,mo.nso,mo.nao,1,&(X[0][0]),mo.nao,&(usotao[0][0]),mo.nao,
	  0,&(MUZ_MO[0][0]),mo.nso);

  free(mu_x_ints); free(mu_y_ints); free(mu_z_ints);
  free_block(X);
  free_block(usotao);
  free_block(MUX_AO);
  free_block(MUY_AO);
  free_block(MUZ_AO);

  /*** Transform the SO dipole integrals to the MO basis ***/

  X = block_matrix(mo.nmo,mo.nmo); /* just a temporary matrix */

  C_DGEMM('t','n',mo.nmo,mo.nmo,mo.nmo,1,&(scf_qt[0][0]),mo.nmo,&(MUX_MO[0][0]),mo.nmo,
	  0,&(X[0][0]),mo.nmo);
  C_DGEMM('n','n',mo.nmo,mo.nmo,mo.nmo,1,&(X[0][0]),mo.nmo,&(scf_qt[0][0]),mo.nmo,
	  0,&(MUX_MO[0][0]),mo.nmo);

  C_DGEMM('t','n',mo.nmo,mo.nmo,mo.nmo,1,&(scf_qt[0][0]),mo.nmo,&(MUY_MO[0][0]),mo.nmo,
	  0,&(X[0][0]),mo.nmo);
  C_DGEMM('n','n',mo.nmo,mo.nmo,mo.nmo,1,&(X[0][0]),mo.nmo,&(scf_qt[0][0]),mo.nmo,
	  0,&(MUY_MO[0][0]),mo.nmo);

  C_DGEMM('t','n',mo.nmo,mo.nmo,mo.nmo,1,&(scf_qt[0][0]),mo.nmo,&(MUZ_MO[0][0]),mo.nmo,
	  0,&(X[0][0]),mo.nmo);
  C_DGEMM('n','n',mo.nmo,mo.nmo,mo.nmo,1,&(X[0][0]),mo.nmo,&(scf_qt[0][0]),mo.nmo,
	  0,&(MUZ_MO[0][0]),mo.nmo);

  free_block(scf_qt);
  free_block(X);

/*   outfile->Printf( "MUX_MOs\n"); */
/*   print_mat(MUX_MO, mo.nmo, mo.nmo, outfile); */
/*   outfile->Printf( "MUY_MOs\n"); */
/*   print_mat(MUY_MO, mo.nmo, mo.nmo, outfile); */
/*   outfile->Printf( "MUZ_MOs\n"); */
/*   print_mat(MUZ_MO, mo.nmo, mo.nmo, outfile); */

  /*** Contract the correlated dipole moment ***/

  for(i=0; i < mo.nmo; i++)
    for(j=0; j < mo.nmo; j++) {
      mu_x += MUX_MO[i][j] * mo.opdm[i][j];
      mu_y += MUY_MO[i][j] * mo.opdm[i][j];
      mu_z += MUZ_MO[i][j] * mo.opdm[i][j];
    }

  for(i=0; i < mo.ndocc; i++) {
    mu_x_ref += MUX_MO[i][i] * 2.0 ;
    mu_y_ref += MUY_MO[i][i] * 2.0 ;
    mu_z_ref += MUZ_MO[i][i] * 2.0 ;
  }

  free_block(MUX_MO); free_block(MUY_MO); free_block(MUZ_MO);

  mu_x_tot = mu_x + mu_x_n + mu_x_ref ;
  mu_y_tot = mu_y + mu_y_n + mu_y_ref ;
  mu_z_tot = mu_z + mu_z_n + mu_z_ref ;

  outfile->Printf("\n\tReference part of electric dipole moment:\n");
  outfile->Printf("\t----------------------------------------------\n");
  outfile->Printf("\tmu(X) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_x_ref* pc_dipmom_au2debye, mu_x_ref* pc_dipmom_au2si, mu_x_ref);
  outfile->Printf("\tmu(Y) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_y_ref* pc_dipmom_au2debye, mu_y_ref* pc_dipmom_au2si, mu_y_ref);
  outfile->Printf("\tmu(Z) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_z_ref* pc_dipmom_au2debye, mu_z_ref* pc_dipmom_au2si, mu_z_ref);
  mu = sqrt(mu_x_ref*mu_x_ref+ mu_y_ref*mu_y_ref+ mu_z_ref*mu_z_ref);
  outfile->Printf("\t|mu|  = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu * pc_dipmom_au2debye, mu * pc_dipmom_au2si, mu);

  outfile->Printf("\n\tCorrelation part of electric dipole moment:\n");
  outfile->Printf("\t----------------------------------------------\n");
  outfile->Printf("\tmu(X) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_x* pc_dipmom_au2debye, mu_x* pc_dipmom_au2si, mu_x);
  outfile->Printf("\tmu(Y) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_y* pc_dipmom_au2debye, mu_y* pc_dipmom_au2si, mu_y);
  outfile->Printf("\tmu(Z) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_z* pc_dipmom_au2debye, mu_z* pc_dipmom_au2si, mu_z);
  mu = sqrt(mu_x*mu_x+ mu_y*mu_y+ mu_z*mu_z);
  outfile->Printf("\t|mu|  = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu * pc_dipmom_au2debye, mu * pc_dipmom_au2si, mu);

  outfile->Printf("\n\tNuclear part of electric dipole moment:\n");
  outfile->Printf("\t----------------------------------------------\n");
  outfile->Printf("\tmu(X) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_x_n * pc_dipmom_au2debye, mu_x_n * pc_dipmom_au2si, mu_x_n);
  outfile->Printf("\tmu(Y) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_y_n * pc_dipmom_au2debye, mu_y_n * pc_dipmom_au2si, mu_y_n);
  outfile->Printf("\tmu(Z) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_z_n * pc_dipmom_au2debye, mu_z_n * pc_dipmom_au2si, mu_z_n);
  mu = sqrt(mu_x_n*mu_x_n + mu_y_n*mu_y_n + mu_z_n*mu_z_n);
  outfile->Printf("\t|mu|  = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu * pc_dipmom_au2debye, mu * pc_dipmom_au2si, mu);

  outfile->Printf("\n\tTotal electric dipole moment:\n");
  outfile->Printf("\t----------------------------------------------\n");
  outfile->Printf("\tmu(X) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_x_tot * pc_dipmom_au2debye, mu_x_tot * pc_dipmom_au2si, mu_x_tot);
  outfile->Printf("\tmu(Y) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_y_tot * pc_dipmom_au2debye, mu_y_tot * pc_dipmom_au2si, mu_y_tot);
  outfile->Printf("\tmu(Z) = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu_z_tot * pc_dipmom_au2debye, mu_z_tot * pc_dipmom_au2si, mu_z_tot);
  mu = sqrt(mu_x_tot*mu_x_tot + mu_y_tot*mu_y_tot + mu_z_tot*mu_z_tot);
  outfile->Printf("\t|mu|  = %8.5lf D,  %15.8e C*m,  %11.8lf a.u.\n",
	  mu * pc_dipmom_au2debye, mu * pc_dipmom_au2si, mu);

  return;
}

}} /* End namespace */
