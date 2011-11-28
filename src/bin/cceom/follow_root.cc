/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
/*
 root_following: returns index of vector in EOM_Cxxx with maximum overlap with
 vector "CCSD Cxxx" in CC3_MISC
 */

#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

int follow_root(int L, double **alpha, int C_irr) {
  dpdfile2 CME, Cme, RME, Rme;
  dpdbuf4 CMNEF, Cmnef, CMnEf, RMNEF, Rmnef, RMnEf;
  char lbl[32];
  double *CR_overlap, tval;
  int i,j;

  CR_overlap = init_array(L);

  /* open CCSD vector "R" from CC3_MISC */
  if (params.eom_ref == 0) {
    dpd_file2_init(&RME, CC3_MISC, C_irr, 0, 1, "CCSD CME");
    dpd_buf4_init(&RMnEf, CC3_MISC, C_irr, 0, 5, 0, 5, 0, "CCSD CMnEf");
  }
  else if (params.eom_ref == 1) {
    dpd_file2_init(&RME, CC3_MISC, C_irr, 0, 1, "CCSD CME");
    dpd_file2_init(&Rme, CC3_MISC, C_irr, 0, 1, "CCSD Cme");
    dpd_buf4_init(&RMNEF, CC3_MISC, C_irr, 2, 7, 2, 7, 0, "CCSD CMNEF");
    dpd_buf4_init(&Rmnef, CC3_MISC, C_irr, 2, 7, 2, 7, 0, "CCSD Cmnef");
    dpd_buf4_init(&RMnEf, CC3_MISC, C_irr, 0, 5, 0, 5, 0, "CCSD CMnEf");
  }
  else if (params.eom_ref == 2) {
    dpd_file2_init(&RME, CC3_MISC, C_irr, 0, 1, "CCSD CME");
    dpd_file2_init(&Rme, CC3_MISC, C_irr, 2, 3, "CCSD Cme");
    dpd_buf4_init(&RMNEF, CC3_MISC, C_irr, 2, 7, 2, 7, 0, "CCSD CMNEF");
    dpd_buf4_init(&Rmnef, CC3_MISC, C_irr, 12, 17, 12, 17, 0, "CCSD Cmnef");
    dpd_buf4_init(&RMnEf, CC3_MISC, C_irr, 22, 28, 22, 28, 0, "CCSD CMnEf");
  }

  /* loop over trial C vectors */
  for (i=0; i<L; ++i) {

    /* read C vector from EOM_Cxxx */
    if (params.eom_ref == 0) {
      sprintf(lbl, "%s %d", "CME", i);
      dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d", "CMnEf", i);
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    }
    else if (params.eom_ref == 1) {
      sprintf(lbl, "%s %d", "CME", i);
      dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d", "Cme", i);
      dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d", "CMNEF", i);
      dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d", "Cmnef", i);
      dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d", "CMnEf", i);
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    }
    else if (params.eom_ref == 2) {
      sprintf(lbl, "%s %d", "CME", i);
      dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
      sprintf(lbl, "%s %d", "Cme", i);
      dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, lbl);
      sprintf(lbl, "%s %d", "CMNEF", i);
      dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
      sprintf(lbl, "%s %d", "Cmnef", i);
      dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, lbl);
      sprintf(lbl, "%s %d", "CMnEf", i);
      dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
    }

    /* dot C vector with R vector */
    tval = 0.0;
    if (params.eom_ref == 0) {
      tval = 2.0 * dpd_file2_dot(&CME, &RME);
      tval += dpd_buf4_dot(&CMnEf, &RMnEf);
    }
    else if (params.eom_ref == 1) {
      tval = dpd_file2_dot(&CME, &RME);
      tval += dpd_file2_dot(&Cme, &Rme);
      tval += dpd_buf4_dot(&CMNEF, &RMNEF);
      tval += dpd_buf4_dot(&Cmnef, &Rmnef);
      tval += dpd_buf4_dot(&CMnEf, &RMnEf);
    }
    else if (params.eom_ref == 2) {
      tval = dpd_file2_dot(&CME, &RME);
      tval += dpd_file2_dot(&Cme, &Rme);
      tval += dpd_buf4_dot(&CMNEF, &RMNEF);
      tval += dpd_buf4_dot(&Cmnef, &Rmnef);
      tval += dpd_buf4_dot(&CMnEf, &RMnEf);
    }

    /* loop over roots and add in overlap */
    for (j=0; j<L; ++j)
      CR_overlap[j] += alpha[i][j] * tval;

    if (params.eom_ref == 0) {
      dpd_file2_close(&CME);
      dpd_buf4_close(&CMnEf);
    }
    else {
      dpd_file2_close(&CME);
      dpd_file2_close(&Cme);
      dpd_buf4_close(&CMNEF);
      dpd_buf4_close(&Cmnef);
      dpd_buf4_close(&CMnEf);
    }
  }

  if (params.eom_ref == 0) {
    dpd_file2_close(&RME);
    dpd_buf4_close(&RMnEf);
  } 
  else {
    dpd_file2_close(&RME);
    dpd_file2_close(&Rme);
    dpd_buf4_close(&RMNEF);
    dpd_buf4_close(&Rmnef);
    dpd_buf4_close(&RMnEf);
  }

  fprintf(outfile,"Overlaps of Rs with EOM CCSD eigenvector:\n");
  for(i=0;i<L;++i) {
    fprintf(outfile,"\t %d  %12.6lf\n", i, CR_overlap[i]);
  }

  /* return index with greatest overlap */
  tval = -1.0;

  for(i=0;i<L;++i) {
    if ( fabs(CR_overlap[i]) > tval) {
      tval = fabs(CR_overlap[i]);
      j = i;
    }
  }

  fprintf(outfile,"follow_root returning: %d\n", j);
  return j;
}

}} // namespace psi::cceom
