/*! \file
    \ingroup CCRESPONSE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

/* sort_pert(): Sorts the specified MO-basis one-electron property
** integrals into CC ordering for use in building the
** similarity-transformed integrals and certain components of the
** total linear response function.
**
** NB: Some integrals are antisymmetric (e.g. L or P integrals), and
** others are symmetric (e.g. Mu integrals), so we must be careful in
** this and subsequent routines.
**
** TDC, 10/05
*/

void sort_pert(const char *pert, double **pertints, int irrep)
{
  int p, q, Gp, Gq, P, Q, i;
  dpdfile2 f;
  char prefix[32], lbl[32];

  sprintf(lbl, "%s_IJ", pert);
  dpd_file2_init(&f, CC_OEI, irrep, 0, 0, lbl);
  dpd_file2_mat_init(&f);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p+moinfo.occ_off[Gp]]];
      for(q=0; q < moinfo.occpi[Gq]; q++) {
        Q = moinfo.qt2pitzer[moinfo.qt_occ[q+moinfo.occ_off[Gq]]];
        f.matrix[Gp][p][q] = pertints[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&f);
  dpd_file2_mat_close(&f);
  dpd_file2_close(&f);

  sprintf(lbl, "%s_AB", pert);
  dpd_file2_init(&f, CC_OEI, irrep, 1, 1, lbl);
  dpd_file2_mat_init(&f);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.virtpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_vir[p+moinfo.vir_off[Gp]]];
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
        Q = moinfo.qt2pitzer[moinfo.qt_vir[q+moinfo.vir_off[Gq]]];
        f.matrix[Gp][p][q] = pertints[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&f);
  dpd_file2_mat_close(&f);
  dpd_file2_close(&f);

  sprintf(lbl, "%s_IA", pert);
  dpd_file2_init(&f, CC_OEI, irrep, 0, 1, lbl);
  dpd_file2_mat_init(&f);
  for(Gp=0; Gp < moinfo.nirreps; Gp++) { /* irrep of left-hand MO */
    Gq = irrep ^ Gp;

    for(p=0; p < moinfo.occpi[Gp]; p++) {
      P = moinfo.qt2pitzer[moinfo.qt_occ[p+moinfo.occ_off[Gp]]];
      for(q=0; q < moinfo.virtpi[Gq]; q++) {
        Q = moinfo.qt2pitzer[moinfo.qt_vir[q+moinfo.vir_off[Gq]]];
        f.matrix[Gp][p][q] = pertints[P][Q];
      }
    }
  }
  dpd_file2_mat_wrt(&f);
  dpd_file2_mat_close(&f);
  dpd_file2_close(&f);

}

}} // namespace psi::ccresponse
