/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void sort_rtd_uhf(struct TD_Params S)
{
  int h, nirreps, nmo, nfzv, nfzc, nclsd, nopen;
  int row, col, i, j, I, J, a, b, A, B, p, q;
  int *aoccpi, *avirtpi, *aocc_off, *avir_off; 
  int *boccpi, *bvirtpi, *bocc_off, *bvir_off; 
  int *aocc_sym, *avir_sym;
  int *bocc_sym, *bvir_sym;
  int *qt_aocc, *qt_avir;
  int *qt_bocc, *qt_bvir;
  double chksum, value;
  dpdfile2 D;

  nmo = moinfo.nmo;
  nfzc = moinfo.nfzc;
  nfzv = moinfo.nfzv;
  nclsd = moinfo.nclsd;
  nopen = moinfo.nopen;
  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi; avirtpi = moinfo.avirtpi;
  boccpi = moinfo.boccpi; bvirtpi = moinfo.bvirtpi;
  aocc_off = moinfo.aocc_off; avir_off = moinfo.avir_off;
  bocc_off = moinfo.bocc_off; bvir_off = moinfo.bvir_off;
  aocc_sym = moinfo.aocc_sym; avir_sym = moinfo.avir_sym;
  bocc_sym = moinfo.bocc_sym; bvir_sym = moinfo.bvir_sym;

  qt_aocc = moinfo.qt_aocc; qt_avir = moinfo.qt_avir;
  qt_bocc = moinfo.qt_bocc; qt_bvir = moinfo.qt_bvir;

  moinfo.rtd_a = block_matrix(nmo,nmo);
  moinfo.rtd_b = block_matrix(nmo,nmo);

  dpd_file2_init(&D, CC_TMP, S.irrep, 0, 0, "RTDIJ");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < aoccpi[h]; i++) {
      I = qt_aocc[aocc_off[h] + i];
      for(j=0; j < aoccpi[h^S.irrep]; j++) {
        J = qt_aocc[aocc_off[h^S.irrep] + j];
        moinfo.rtd_a[I][J] += D.matrix[h][i][j];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_TMP, S.irrep, 1, 1, "RTDAB");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(a=0; a < avirtpi[h]; a++) {
      A = qt_avir[avir_off[h] + a];
      for(b=0; b < avirtpi[h^S.irrep]; b++) {
        B = qt_avir[avir_off[h^S.irrep] + b];
        moinfo.rtd_a[A][B] += D.matrix[h][a][b];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Note that this component of the density is stored occ-vir */
  dpd_file2_init(&D, CC_TMP, S.irrep, 0, 1, "RTDAI");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < aoccpi[h]; i++) {
      I = qt_aocc[aocc_off[h] + i];
      for(a=0; a < avirtpi[h^S.irrep]; a++) {
        A = qt_avir[avir_off[h^S.irrep] + a];
        moinfo.rtd_a[A][I] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_TMP, S.irrep, 0, 1, "RTDIA");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < aoccpi[h]; i++) {
      I = qt_aocc[aocc_off[h] + i];
      for(a=0; a < avirtpi[h^S.irrep]; a++) {
        A = qt_avir[avir_off[h^S.irrep] + a];
        moinfo.rtd_a[I][A] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_TMP, S.irrep, 2, 2, "RTDij");
  dpd_file2_mat_init(&D); 
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < boccpi[h]; i++) { 
      I = qt_bocc[bocc_off[h] + i];
      for(j=0; j < boccpi[h^S.irrep]; j++) {
        J = qt_bocc[bocc_off[h^S.irrep] + j];
        moinfo.rtd_b[I][J] += D.matrix[h][i][j];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_TMP, S.irrep, 3, 3, "RTDab");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(a=0; a < bvirtpi[h]; a++) {
      A = qt_bvir[bvir_off[h] + a];
      for(b=0; b < bvirtpi[h^S.irrep]; b++) {
        B = qt_bvir[bvir_off[h^S.irrep] + b];
        moinfo.rtd_b[A][B] += D.matrix[h][a][b];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Note that this component of the density is stored occ-vir */
  dpd_file2_init(&D, CC_TMP, S.irrep, 2, 3, "RTDai");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < boccpi[h]; i++) {
      I = qt_bocc[bocc_off[h] + i];
      for(a=0; a < bvirtpi[h^S.irrep]; a++) {
        A = qt_bvir[bvir_off[h^S.irrep] + a];
        moinfo.rtd_b[A][I] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_TMP, S.irrep, 2, 3, "RTDia");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < boccpi[h]; i++) {
      I = qt_bocc[bocc_off[h] + i];
      for(a=0; a < bvirtpi[h^S.irrep]; a++) {
        A = qt_bvir[bvir_off[h^S.irrep] + a];
        moinfo.rtd_b[I][A] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /*print_mat(moinfo.rtd_a,nmo,nmo,outfile);*/
  /*print_mat(moinfo.rtd_b,nmo,nmo,outfile);*/

  return;
}

}} // namespace psi::ccdensity
