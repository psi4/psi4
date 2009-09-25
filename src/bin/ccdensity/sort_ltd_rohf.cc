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

void sort_ltd_rohf(struct TD_Params S) 
{
  int h, nirreps, nmo, nfzv, nfzc, nclsd, nopen;
  int row, col, i, j, I, J, a, b, A, B, p, q;
  int *occpi, *virtpi, *occ_off, *vir_off; 
  int *occ_sym, *vir_sym, *openpi;
  int *qt_occ, *qt_vir;
  double chksum, value;
  psio_address next;
  dpdfile2 D;

  nmo = moinfo.nmo;
  nfzc = moinfo.nfzc;
  nfzv = moinfo.nfzv;
  nclsd = moinfo.nclsd;
  nopen = moinfo.nopen;
  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;
  qt_occ = moinfo.qt_occ; qt_vir = moinfo.qt_vir;

  moinfo.ltd = block_matrix(nmo, nmo);

  dpd_file2_init(&D, CC_TMP, S.irrep, 0, 0, "LTDIJ");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(j=0; j < occpi[h^S.irrep]; j++) {
        J = qt_occ[occ_off[h^S.irrep] + j];
        moinfo.ltd[I][J] += D.matrix[h][i][j];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_TMP, S.irrep, 1, 1, "LTDAB");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(a=0; a < (virtpi[h] - openpi[h]); a++) {
      A = qt_vir[vir_off[h] + a];
      for(b=0; b < (virtpi[h^S.irrep] - openpi[h^S.irrep]); b++) {
        B = qt_vir[vir_off[h^S.irrep] + b];
        moinfo.ltd[A][B] += D.matrix[h][a][b];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Note that this component of the density is stored occ-vir */
  dpd_file2_init(&D, CC_TMP, S.irrep, 0, 1, "LTDAI");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(a=0; a < (virtpi[h^S.irrep] - openpi[h^S.irrep]); a++) {
        A = qt_vir[vir_off[h^S.irrep] + a];
        moinfo.ltd[A][I] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_TMP, S.irrep, 0, 1, "LTDIA");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++) {
      I = qt_occ[occ_off[h] + i];
      for(a=0; a < (virtpi[h^S.irrep] - openpi[h^S.irrep]); a++) {
        A = qt_vir[vir_off[h^S.irrep] + a];
        moinfo.ltd[I][A] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_TMP, S.irrep, 0, 0, "LTDij");
  dpd_file2_mat_init(&D); 
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < (occpi[h] - openpi[h]); i++) { 
      I = qt_occ[occ_off[h] + i];
      for(j=0; j < (occpi[h^S.irrep] - openpi[h^S.irrep]); j++) {
        J = qt_occ[occ_off[h^S.irrep] + j];
        moinfo.ltd[I][J] += D.matrix[h][i][j];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_TMP, S.irrep, 1, 1, "LTDab");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(a=0; a < virtpi[h]; a++) {
      A = qt_vir[vir_off[h] + a];
      for(b=0; b < virtpi[h^S.irrep]; b++) {
        B = qt_vir[vir_off[h^S.irrep] + b];
        moinfo.ltd[A][B] += D.matrix[h][a][b];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /* Note that this component of the density is stored occ-vir */
  dpd_file2_init(&D, CC_TMP, S.irrep, 0, 1, "LTDai");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < (occpi[h] - openpi[h]); i++) {
      I = qt_occ[occ_off[h] + i];
      for(a=0; a < virtpi[h^S.irrep]; a++) {
        A = qt_vir[vir_off[h^S.irrep] + a];
        moinfo.ltd[A][I] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_TMP, S.irrep, 0, 1, "LTDia");
  dpd_file2_mat_init(&D);
  dpd_file2_mat_rd(&D);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < (occpi[h] - openpi[h]); i++) {
      I = qt_occ[occ_off[h] + i];
      for(a=0; a < virtpi[h^S.irrep]; a++) {
        A = qt_vir[vir_off[h^S.irrep] + a];
        moinfo.ltd[I][A] += D.matrix[h][i][a];
      }
    }
  }
  dpd_file2_mat_close(&D);
  dpd_file2_close(&D);

  /*mat_print(moinfo.ltd,nmo,nmo,outfile);*/

  return;
}

}} // namespace psi::ccdensity
