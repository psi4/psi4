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
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_sf_sort_I(void);
void uhf_sf_sort_I(void);

void sort_I(void)
{
  if(params.ref == 0) rhf_sf_sort_I();
  else if(params.ref == 2) uhf_sf_sort_I();
}

void rhf_sort_W(void)
{
  int h, nirreps, nmo;
  int row, col, i, j, I, J, a, b, A, B, p, q;
  int *occpi, *virpi, *occ_off, *vir_off; 
  int *occ_sym, *vir_sym;
  int *qt_occ, *qt_vir;
  double **O; 
  double value;
  dpdfile2 W;

  nmo = mo.nmo;
  nirreps = mo.nirreps;
  occpi = mo.occpi; 
  virpi = mo.virpi;
  occ_off = mo.occ_off; 
  vir_off = mo.vir_off;
  occ_sym = mo.occ_sym; 
  vir_sym = mo.vir_sym;
  qt_occ = mo.qt_occ; 
  qt_vir = mo.qt_vir;

  O = block_matrix(nmo,nmo);

  global_dpd_->file2_init(&W, PSIF_CC_OEI, 0, 0, 0, "WIJ");
  global_dpd_->file2_mat_init(&W);
  global_dpd_->file2_mat_rd(&W);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(j=0; j < occpi[h]; j++) {
        J = qt_occ[occ_off[h] + j];
        O[I][J] += W.matrix[h][i][j];
      }
    }
  }
  global_dpd_->file2_mat_close(&W);
  global_dpd_->file2_close(&W);

  global_dpd_->file2_init(&W, PSIF_CC_OEI, 0, 1, 1, "WAB");
  global_dpd_->file2_mat_init(&W);
  global_dpd_->file2_mat_rd(&W);
  for(h=0; h < nirreps; h++) {
    for(a=0; a < virpi[h]; a++) {
      A = qt_vir[vir_off[h] + a];
      for(b=0; b < virpi[h]; b++) {
        B = qt_vir[vir_off[h] + b];
        O[A][B] += W.matrix[h][a][b];
      }
    }
  }
  global_dpd_->file2_mat_close(&W);
  global_dpd_->file2_close(&W);

  global_dpd_->file2_init(&W, PSIF_CC_OEI, 0, 1, 0, "WAI");
  global_dpd_->file2_mat_init(&W);
  global_dpd_->file2_mat_rd(&W);
  for(h=0; h < nirreps; h++) {
    for(a=0; a < virpi[h]; a++) {
      A = qt_vir[vir_off[h] + a];
      for(i=0; i < occpi[h]; i++) {
        I = qt_occ[occ_off[h] + i];
        O[A][I] += W.matrix[h][a][i];
      }
    }
  }
  global_dpd_->file2_mat_close(&W);
  global_dpd_->file2_close(&W);

  global_dpd_->file2_init(&W, PSIF_CC_OEI, 0, 0, 1, "WIA");
  global_dpd_->file2_mat_init(&W);
  global_dpd_->file2_mat_rd(&W);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(a=0; a < virpi[h]; a++) {
        A = qt_vir[vir_off[h] + a];
        O[I][A] += W.matrix[h][i][a];
      }
    }
  }
  global_dpd_->file2_mat_close(&W);
  global_dpd_->file2_close(&W);

  for(p=0; p < nmo; p++) {
    for(q=0; q < p; q++) {
      value = 0.5 * (O[p][q] + O[q][p]);
      O[p][q] = O[q][p] = value;
    }
  }

  for(p=0; p < nmo; p++) {
    for(q=0; q < nmo; q++) {
      O[p][q] *= -2.0;
    }
  }

  /*
  outfile->Printf("\n\tEnergy Weighted OPDM:\n");
  print_mat(O,nmo,nmo,outfile);
  */

  psio_open(PSIF_MO_LAG, PSIO_OPEN_OLD);
  psio_write_entry(PSIF_MO_LAG, "MO-basis Lagrangian", (char *)O[0],
                   sizeof(double)*nmo*nmo);
  psio_close(PSIF_MO_LAG, 1);

  free_block(O);
}

void uhf_sort_W(void)
{

}

void rhf_sf_sort_I(void)
{
  int h, nirreps, nmo, nfzv, nfzc, nclsd, nopen;
  int row, col, i, j, I, J, a, b, A, B, p, q;
  int *occpi, *virpi, *occ_off, *vir_off; 
  int *occ_sym, *vir_sym, *openpi;
  int *qt_occ, *qt_vir;
  double **O, chksum, value;
  dpdfile2 D;

  nmo = mo.nmo;
  nirreps = mo.nirreps;
  occpi = mo.occpi; 
  virpi = mo.virtpi;
  occ_off = mo.occ_off; 
  vir_off = mo.vir_off;
  occ_sym = mo.occ_sym; 
  vir_sym = mo.vir_sym;
  qt_occ = mo.qt_occ; 
  qt_vir = mo.qt_vir;

  O = block_matrix(nmo,nmo);

  /* Sort alpha components first */
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "I(I,J)");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(j=0; j < occpi[h]; j++) {
        J = qt_occ[occ_off[h] + j];
        O[I][J] += D.matrix[h][i][j];
      }
    }
  }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "I'AB");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(a=0; a < virpi[h]; a++) {
      A = qt_vir[vir_off[h] + a];
      for(b=0; b < virpi[h]; b++) {
        B = qt_vir[vir_off[h] + b];
        O[A][B] += D.matrix[h][a][b];
      }
    }
  }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "I(I,A)");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(a=0; a < virpi[h]; a++) {
        A = qt_vir[vir_off[h] + a];
        O[A][I] += D.matrix[h][i][a];
	O[I][A] += D.matrix[h][i][a];
      }
    }
  }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  /* Sort beta components */
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 0, "I(i,j)");
  global_dpd_->file2_mat_init(&D); 
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) { 
      I = qt_occ[occ_off[h] + i];
      for(j=0; j < occpi[h]; j++) {
        J = qt_occ[occ_off[h] + j];
        O[I][J] += D.matrix[h][i][j];
      }
    }
  }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 1, "I'ab");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(a=0; a < virpi[h]; a++) {
      A = qt_vir[vir_off[h] + a];
      for(b=0; b < virpi[h]; b++) {
        B = qt_vir[vir_off[h] + b];
        O[A][B] += D.matrix[h][a][b];
      }
    }
  }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 0, 1, "I(i,a)");
  global_dpd_->file2_mat_init(&D);
  global_dpd_->file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(a=0; a < virpi[h]; a++) {
        A = qt_vir[vir_off[h] + a];
        O[A][I] += D.matrix[h][i][a];
	O[I][A] += D.matrix[h][i][a];
      }
    }
  }
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  /* Symmetrize the Lagrangian */
  for(p=0; p < nmo; p++) {
    for(q=0; q < p; q++) {
      value = 0.5*(O[p][q] + O[q][p]);
      O[p][q] = O[q][p] = value;
    }
  }

  /* Multiply the Lagrangian by -2.0 for the final energy derivative
     expression */
  for(p=0; p < nmo; p++) 
    for(q=0; q < nmo; q++) 
      O[p][q] *= -2.0;

  mo.I = O;
}

void uhf_sf_sort_I(void)
{

}

}} /* End namespaces */
