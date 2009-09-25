/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
/*
 restart_with_root: copies C's from position prop_root to position 0 
 in EOM_Cxxx files 

 also put copy in CC3_MISC for root_following
*/

#include <cstdio>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void restart_with_root(int prop_root, int C_irr) {
  dpdfile2 CME, Cme;
  dpdbuf4 CMNEF, Cmnef, CMnEf;
  char lbl[32];

  fprintf(outfile,"Copying root %d to start of EOM_Cxxx files.\n",prop_root+1);

  if (params.eom_ref == 0) {
    sprintf(lbl, "CME %d", prop_root);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    dpd_file2_copy(&CME, EOM_CME, "CME 0");
    dpd_file2_close(&CME);

    sprintf(lbl, "CMnEf %d", prop_root);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_copy(&CMnEf, EOM_CMnEf, "CMnEf 0");
    dpd_buf4_close(&CMnEf);
  }
  else if (params.eom_ref == 1) {
    sprintf(lbl, "CME %d", prop_root);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    dpd_file2_copy(&CME, EOM_CME, "CME 0");
    dpd_file2_close(&CME);

    sprintf(lbl, "Cme %d", prop_root);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
    dpd_file2_copy(&Cme, EOM_Cme, "Cme 0");
    dpd_file2_close(&Cme);

    sprintf(lbl, "CMNEF %d", prop_root);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
    dpd_buf4_copy(&CMNEF, EOM_CMNEF, "CMNEF 0");
    dpd_buf4_close(&CMNEF);

    sprintf(lbl, "Cmnef %d", prop_root);
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
    dpd_buf4_copy(&Cmnef, EOM_Cmnef, "Cmnef 0");
    dpd_buf4_close(&Cmnef);

    sprintf(lbl, "CMnEf %d", prop_root);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_copy(&CMnEf, EOM_CMnEf, "CMnEf 0");
    dpd_buf4_close(&CMnEf);
  }
  else {
    sprintf(lbl, "CME %d", prop_root);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    dpd_file2_copy(&CME, EOM_CME, "CME 0");
    dpd_file2_close(&CME);

    sprintf(lbl, "Cme %d", prop_root);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, lbl);
    dpd_file2_copy(&Cme, EOM_Cme, "Cme 0");
    dpd_file2_close(&Cme);

    sprintf(lbl, "CMNEF %d", prop_root);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
    dpd_buf4_copy(&CMNEF, EOM_CMNEF, "CMNEF 0");
    dpd_buf4_close(&CMNEF);

    sprintf(lbl, "Cmnef %d", prop_root);
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, lbl);
    dpd_buf4_copy(&Cmnef, EOM_Cmnef, "Cmnef 0");
    dpd_buf4_close(&Cmnef);

    sprintf(lbl, "CMnEf %d", prop_root);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
    dpd_buf4_copy(&CMnEf, EOM_CMnEf, "CMnEf 0");
    dpd_buf4_close(&CMnEf);
  }
  return;
}

/*
 save_C_ccsd: copies C's from position prop_root to CC3_MISC file
*/

void save_C_ccsd(int prop_root, int C_irr) {
  dpdfile2 CME, Cme;
  dpdbuf4 CMNEF, Cmnef, CMnEf;
  char lbl[32];

  fprintf(outfile,"Copying root %d to CC3_MISC file.\n",prop_root+1);

  if (params.eom_ref == 0) {
    sprintf(lbl, "CME %d", prop_root);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    dpd_file2_copy(&CME, CC3_MISC, "CCSD CME");
    dpd_file2_close(&CME);

    sprintf(lbl, "CMnEf %d", prop_root);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_copy(&CMnEf, CC3_MISC, "CCSD CMnEf");
    dpd_buf4_close(&CMnEf);
  }
  else if (params.eom_ref == 1) {
    sprintf(lbl, "CME %d", prop_root);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    dpd_file2_copy(&CME, CC3_MISC, "CCSD CME");
    dpd_file2_close(&CME);

    sprintf(lbl, "Cme %d", prop_root);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
    dpd_file2_copy(&Cme, CC3_MISC, "CCSD Cme");
    dpd_file2_close(&Cme);

    sprintf(lbl, "CMNEF %d", prop_root);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
    dpd_buf4_copy(&CMNEF, CC3_MISC, "CCSD CMNEF");
    dpd_buf4_close(&CMNEF);

    sprintf(lbl, "Cmnef %d", prop_root);
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, lbl);
    dpd_buf4_copy(&Cmnef, CC3_MISC, "CCSD Cmnef");
    dpd_buf4_close(&Cmnef);

    sprintf(lbl, "CMnEf %d", prop_root);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_copy(&CMnEf, CC3_MISC, "CCSD CMnEf");
    dpd_buf4_close(&CMnEf);
  }
  else {
    sprintf(lbl, "CME %d", prop_root);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
    dpd_file2_copy(&CME, CC3_MISC, "CCSD CME");
    dpd_file2_close(&CME);

    sprintf(lbl, "Cme %d", prop_root);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, lbl);
    dpd_file2_copy(&Cme, CC3_MISC, "CCSD Cme");
    dpd_file2_close(&Cme);

    sprintf(lbl, "CMNEF %d", prop_root);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, lbl);
    dpd_buf4_copy(&CMNEF, CC3_MISC, "CCSD CMNEF");
    dpd_buf4_close(&CMNEF);

    sprintf(lbl, "Cmnef %d", prop_root);
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, lbl);
    dpd_buf4_copy(&Cmnef, CC3_MISC, "CCSD Cmnef");
    dpd_buf4_close(&Cmnef);

    sprintf(lbl, "CMnEf %d", prop_root);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
    dpd_buf4_copy(&CMnEf, CC3_MISC, "CCSD CMnEf");
    dpd_buf4_close(&CMnEf);
  }
  return;
}

}} // namespace psi::cceom
