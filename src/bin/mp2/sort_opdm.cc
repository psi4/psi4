/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_sort_opdm(void);
void uhf_sort_opdm(void);
void rhf_sf_sort_opdm(void);
void uhf_sf_sort_opdm(void);

void sort_opdm(void)
{
  if(params.gradient) {
    if(params.ref == 0) rhf_sf_sort_opdm();
    else if(params.ref == 2) uhf_sf_sort_opdm();
  }
  else {
    if(params.ref == 0) rhf_sort_opdm();
    else if(params.ref == 2) uhf_sort_opdm();
  }
}

void rhf_sort_opdm(void)
{
  int h, nirreps, nmo, ndocc;
  int nfzc, nfzv;
  int row, col, i, j, I, J, a, b, A, B, p, q;
  int *occpi, *virpi, *occ_off, *vir_off; 
  int *occ_sym, *vir_sym;
  int *qt_occ, *qt_vir;
  double **OPDM;
  double value;
  dpdfile2 D;

  nmo = mo.nmo;
  nirreps = mo.nirreps;
  ndocc = mo.ndocc;
  nfzc = mo.nfzdocc;
  nfzv = mo.nfzvirt;
  occpi = mo.occpi; 
  virpi = mo.virpi;
  occ_off = mo.occ_off; 
  vir_off = mo.vir_off;
  occ_sym = mo.occ_sym; 
  vir_sym = mo.vir_sym;
  qt_occ = mo.qt_occ; 
  qt_vir = mo.qt_vir;

  OPDM = block_matrix(nmo,nmo);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(j=0; j < occpi[h]; j++) {
        J = qt_occ[occ_off[h] + j];
        OPDM[I][J] += D.matrix[h][i][j];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(a=0; a < virpi[h]; a++) {
      A = qt_vir[vir_off[h] + a];
      for(b=0; b < virpi[h]; b++) {
        B = qt_vir[vir_off[h] + b];
        OPDM[A][B] += D.matrix[h][a][b];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  if(params.relax_opdm || params.gradient) {

    dpd_file2_init(&D, CC_OEI, 0, 1, 0, "DAI");
    dpd_file2_mat_init(&D);
    dpd_file2_mat_rd(&D);
    for(h=0; h < nirreps; h++) {
      for(a=0; a < virpi[h]; a++) {
        A = qt_vir[vir_off[h] + a];
        for(i=0; i < occpi[h]; i++) {
          I = qt_occ[occ_off[h] + i];
          OPDM[A][I] += D.matrix[h][a][i];
          OPDM[I][A] += D.matrix[h][a][i];
        }
      }
    }
    dpd_file2_mat_close(&D);
    dpd_file2_close(&D);

  }

  /* Symmetrize the OPDM */

  for(p=0; p < nmo; p++) {
    for(q=0; q < p; q++) {
      value = 0.5 * (OPDM[p][q] + OPDM[q][p]);
      OPDM[p][q] = OPDM[q][p] = value;
    }
  }

  /* Add Reference Contribution */

  for(i=0; i< ndocc; i++) 
    OPDM[i][i] += 2; 

  /* Write OPDM to disk */

  psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
  psio_write_entry(PSIF_MO_OPDM, "MO-basis OPDM", (char *)OPDM[0],
                   sizeof(double)*nmo*nmo);
  psio_close(PSIF_MO_OPDM, 1);

  free_block(OPDM);
}

void uhf_sort_opdm(void)
{
  int h, nirreps, nmo, ndocc, nsocc;
  int nfzc, nfzv;
  int row, col, i, j, I, J, a, b, A, B, p, q;
  int *aoccpi, *avirpi, *aocc_off, *avir_off; 
  int *boccpi, *bvirpi, *bocc_off, *bvir_off; 
  int *aocc_sym, *avir_sym;
  int *bocc_sym, *bvir_sym;
  int *qt_aocc, *qt_avir;
  int *qt_bocc, *qt_bvir;
  double **AOPDM, **BOPDM;
  double value;
  dpdfile2 D;

  nmo = mo.nmo;
  nirreps = mo.nirreps;
  ndocc = mo.ndocc;
  nsocc = mo.nsocc;
  nfzc = mo.nfzdocc;
  nfzv = mo.nfzvirt;
  aoccpi = mo.aoccpi; 
  avirpi = mo.avirpi;
  aocc_off = mo.aocc_off; 
  avir_off = mo.avir_off;
  aocc_sym = mo.aocc_sym; 
  avir_sym = mo.avir_sym;
  qt_aocc = mo.qt_aocc; 
  qt_avir = mo.qt_avir;
  boccpi = mo.boccpi; 
  bvirpi = mo.bvirpi;
  bocc_off = mo.bocc_off; 
  bvir_off = mo.bvir_off;
  bocc_sym = mo.bocc_sym; 
  bvir_sym = mo.bvir_sym;
  qt_bocc = mo.qt_bocc; 
  qt_bvir = mo.qt_bvir;

  AOPDM = block_matrix(nmo,nmo);  
  BOPDM = block_matrix(nmo,nmo);  

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
   
  for(h=0; h < nirreps; h++) {
    for(i=0; i < aoccpi[h]; i++) {
      I = qt_aocc[aocc_off[h] + i];
      for(j=0; j < aoccpi[h]; j++) {
	J = qt_aocc[aocc_off[h] + j];
	AOPDM[I][J] += D.matrix[h][i][j];
      }
    }
  }
   
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 2, 2, "Dij");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
   
  for(h=0; h < nirreps; h++) {
    for(i=0; i < boccpi[h]; i++) {
      I = qt_bocc[bocc_off[h] + i];
      for(j=0; j < boccpi[h]; j++) {
	J = qt_bocc[bocc_off[h] + j];
	BOPDM[I][J] += D.matrix[h][i][j];
      }
    }
  }
   
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);
    
  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
    
  for(h=0; h < nirreps; h++) {
    for(a=0; a < avirpi[h]; a++) {
      A = qt_avir[avir_off[h] + a];
      for(b=0; b < avirpi[h]; b++) {
	B = qt_avir[avir_off[h] + b];
	AOPDM[A][B] += D.matrix[h][a][b];
      }
    }
  }

  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 3, 3, "Dab");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
    
  for(h=0; h < nirreps; h++) {
    for(a=0; a < bvirpi[h]; a++) {
      A = qt_bvir[bvir_off[h] + a];
      for(b=0; b < bvirpi[h]; b++) {
	B = qt_bvir[bvir_off[h] + b];
	BOPDM[A][B] += D.matrix[h][a][b];
      }
    }
  }
    
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Symmetrize AOPDM and BOPDM */
    
  for(p=0; p < nmo; p++) {
    for(q=0; q < p; q++) {
      value = 0.5 * (AOPDM[p][q] + AOPDM[q][p]);
      AOPDM[p][q] = AOPDM[q][p] = value;
    }
  }

  for(p=0; p < nmo; p++) {
    for(q=0; q < p; q++) {
      value = 0.5 * (BOPDM[p][q] + BOPDM[q][p]);
      BOPDM[p][q] = BOPDM[q][p] = value;
    }
  }

  /* Add Reference Contribution */

  for(i=0; i< (ndocc+nsocc); i++) 
    AOPDM[i][i] += 1; 

  for(i=0; i< ndocc; i++) 
    BOPDM[i][i] += 1; 

  /* Write AOPDM and BOPDM to disk */

  psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
  psio_write_entry(PSIF_MO_OPDM,"MO-basis Alpha OPDM",(char*)AOPDM[0],
                   sizeof(double)*nmo*nmo);
  psio_write_entry(PSIF_MO_OPDM,"MO-basis Beta OPDM",(char*)BOPDM[0],
                   sizeof(double)*nmo*nmo);
  psio_close(PSIF_MO_OPDM, 1);

  free_block(AOPDM);
  free_block(BOPDM);
}

void rhf_sf_sort_opdm(void)
{
  int h, nirreps, nmo, ndocc;
  int nfzc, nfzv;
  int row, col, i, j, I, J, a, b, A, B, p, q;
  int *occpi, *virpi, *occ_off, *vir_off; 
  int *occ_sym, *vir_sym;
  int *qt_occ, *qt_vir;
  double **OPDM;
  double value;
  dpdfile2 D;

  nmo = mo.nmo;
  nirreps = mo.nirreps;
  ndocc = mo.ndocc;
  nfzc = mo.nfzdocc;
  nfzv = mo.nfzvirt;
  occpi = mo.occpi; 
  virpi = mo.virpi;
  occ_off = mo.occ_off; 
  vir_off = mo.vir_off;
  occ_sym = mo.occ_sym; 
  vir_sym = mo.vir_sym;
  qt_occ = mo.qt_occ; 
  qt_vir = mo.qt_vir;

  OPDM = block_matrix(nmo,nmo);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(j=0; j < occpi[h]; j++) {
        J = qt_occ[occ_off[h] + j];
        OPDM[I][J] += D.matrix[h][i][j];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(j=0; j < occpi[h]; j++) {
        J = qt_occ[occ_off[h] + j];
        OPDM[I][J] += D.matrix[h][i][j];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(a=0; a < virpi[h]; a++) {
      A = qt_vir[vir_off[h] + a];
      for(b=0; b < virpi[h]; b++) {
        B = qt_vir[vir_off[h] + b];
        OPDM[A][B] += D.matrix[h][a][b];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(a=0; a < virpi[h]; a++) {
      A = qt_vir[vir_off[h] + a];
      for(b=0; b < virpi[h]; b++) {
        B = qt_vir[vir_off[h] + b];
        OPDM[A][B] += D.matrix[h][a][b];
      }
    }
  }
  dpd_file2_mat_close(&D);

  /* Note that this component of the density is stored occ-vir */
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DAI");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(a=0; a < virpi[h]; a++) {
        A = qt_vir[vir_off[h] + a];
        OPDM[A][I] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);
  
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(a=0; a < virpi[h]; a++) {
        A = qt_vir[vir_off[h] + a];
        OPDM[I][A] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Note that this component of the density is stored occ-vir */
  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dai");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(a=0; a < virpi[h]; a++) {
        A = qt_vir[vir_off[h] + a];
        OPDM[A][I] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(a=0; a < virpi[h]; a++) {
        A = qt_vir[vir_off[h] + a];
        OPDM[I][A] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Symmetrize the OPDM */

  for(p=0; p < nmo; p++) {
    for(q=0; q < p; q++) {
      value = 0.5 * (OPDM[p][q] + OPDM[q][p]);
      OPDM[p][q] = OPDM[q][p] = value;
    }
  }

  mo.opdm = OPDM;
}

void uhf_sf_sort_opdm(void)
{

}

}} /* End namespaces */
