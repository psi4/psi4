/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_sf_write_data(void);
void uhf_sf_write_data(void);

void write_data(void)
{
  if(params.ref == 0) rhf_sf_write_data();
  else if(params.ref == 2) uhf_sf_write_data();
}

void rhf_sf_write_data(void)
{
  int nirreps, nmo, nfzv, nfzc;
  int nclsd, nopen;
  int *qt_occ, *qt_vir;
  int h, row, col, p, q, r, s;
  int i, j, m;
  double value;
  dpdbuf4 G;
  struct iwlbuf OutBuf;

  iwl_buf_init(&OutBuf, PSIF_MO_TPDM, 1E-14, 0, 0);

  nmo = mo.nmo;
  nfzc = mo.nfzdocc;
  nfzv = mo.nfzvirt;
  nclsd = mo.ndocc - mo.nfzdocc;
  nopen = mo.nsocc;
  qt_occ = mo.qt_occ;  
  qt_vir = mo.qt_vir;
  nirreps = mo.nirreps;

  for(p=nfzc; p < (nmo - nfzv); p++) {
    for(q=nfzc; q < (nmo - nfzv); q++) {
      value = mo.opdm[p][q];
      for(m=0; m < nfzc; m++) {
	iwl_buf_wrt_val(&OutBuf, p, q, m, m,value,0,outfile,0);
     	iwl_buf_wrt_val(&OutBuf, p, m, m, q,-0.5*value,0,outfile,0);
      }
    }
  }

  /*** One-electron component ***/

  for(i=0; i < (nfzc + nclsd); i++)
      mo.opdm[i][i] += 2.0;

  for(i=nfzc + nclsd; i < (nfzc + nclsd + nopen); i++)
      mo.opdm[i][i] += 1.0;

  /*** Two-electron component ***/

  /* docc-docc */
  for(i=0; i < (nfzc + nclsd); i++) {
    iwl_buf_wrt_val(&OutBuf, i, i, i, i, 1.0, 0, outfile, 0);
    for(j=0; j < i; j++) {
      iwl_buf_wrt_val(&OutBuf, i, i, j, j, 2.0, 0, outfile, 0);
      iwl_buf_wrt_val(&OutBuf, i, j, j, i,-1.0, 0, outfile, 0);
    }
  }

  /* socc-docc && socc-socc*/
  for(i=(nfzc + nclsd); i < (nfzc + nclsd + nopen); i++) {
    for(j=0; j < (nfzc + nclsd); j++) {
      iwl_buf_wrt_val(&OutBuf, i, i, j, j, 1.0, 0, outfile, 0);
      iwl_buf_wrt_val(&OutBuf, i, j, j, i,-0.5, 0, outfile, 0);
    }
    for(j=(nfzc + nclsd); j < i; j++) {
      iwl_buf_wrt_val(&OutBuf, i, i, j, j, 0.5, 0, outfile, 0);
      iwl_buf_wrt_val(&OutBuf, i, j, j, i,-0.5, 0, outfile, 0);
    }
  }

  psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
  psio_write_entry(PSIF_MO_OPDM, "MO-basis OPDM", (char *) mo.opdm[0],
		   sizeof(double)*nmo*nmo);
  psio_close(PSIF_MO_OPDM, 1);

  psio_open(PSIF_MO_LAG, PSIO_OPEN_OLD);
  psio_write_entry(PSIF_MO_LAG, "MO-basis Lagrangian", (char *) mo.I[0],
		   sizeof(double)*nmo*nmo);
  psio_close(PSIF_MO_LAG, 1);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
  dpd_buf4_sort(&G, CC_TMP0, prqs, 0, 0, "G(IK,JL)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 0, 0, 0, 0, 0, "G(IK,JL)");
  dpd_buf4_dump(&G, &OutBuf, qt_occ, qt_occ, qt_occ, qt_occ, 1, 0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
  dpd_buf4_sort(&G, CC_TMP0, prqs, 0, 10, "G(IK,JA)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 0, 10, 0, 10, 0, "G(IK,JA)");
  
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0];
      q = G.params->roworb[h][row][1];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0];
	s = G.params->colorb[h][col][1];

	if((qt_occ[q] == qt_vir[s]) && (p == r))
	  G.matrix[h][row][col] *= 2;
      }
    }

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }

  dpd_buf4_dump(&G, &OutBuf, qt_occ, qt_occ, qt_occ, qt_vir, 0, 0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
  dpd_buf4_sort(&G, CC_TMP9, prqs, 10, 10, "G(IA,JB)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP9, 0, 10, 10, 10, 10, 0, "G(IA,JB)");
  dpd_buf4_symm(&G);
  dpd_buf4_dump(&G, &OutBuf, qt_occ, qt_vir, qt_occ, qt_vir, 1, 0);
  dpd_buf4_close(&G);

  dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
  dpd_buf4_sort(&G, CC_TMP0, prqs, 0, 5, "G(IJ,AB)");
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, CC_TMP0, 0, 0, 5, 0, 5, 0, "G(IJ,AB)");
  dpd_buf4_scm(&G, 0.5);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);

    for(row=0; row < G.params->rowtot[h]; row++) {
      p = G.params->roworb[h][row][0];
      q = G.params->roworb[h][row][1];
      for(col=0; col < G.params->coltot[h]; col++) {
	r = G.params->colorb[h][col][0];
	s = G.params->colorb[h][col][1];

	if((qt_occ[p] == qt_vir[r]) && (qt_occ[q] == qt_vir[s]))
	  G.matrix[h][row][col] *= 2;
      }
    }

    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  
  dpd_buf4_dump(&G, &OutBuf, qt_occ, qt_occ, qt_vir, qt_vir, 0, 0);
  dpd_buf4_close(&G);

  iwl_buf_flush(&OutBuf, 1);
  iwl_buf_close(&OutBuf, 1);
}

void uhf_sf_write_data(void)
{

}

}} /* End namespaces */
