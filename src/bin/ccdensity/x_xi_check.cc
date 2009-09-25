/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
/*
**  X_XI_CHECK: check sum for xi
*/

#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

extern double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
        dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);
                                                                                
extern double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);

extern void c_cleanSS(dpdfile2 *CME, dpdfile2 *Cme);

void c_clean(dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);

void x_xi_check(char *term_lbl) {
  dpdfile2 Xia, XIA;
  dpdbuf4 XIJAB, Xijab, XIjAb, XIjbA;
  static double old_norm=0;
  double norm,dotval;
  char lbl[80];
  int irrep;
  irrep = params.G_irr;
  
  /*
  if (!strcmp(term_lbl,"reset"))  {
    fprintf(outfile,"resetting norm\n");
    old_norm = 0;
    return;
  }
  */
 
  if (params.ref == 0) {
    dpd_file2_init(&XIA, EOM_XI, irrep, 0, 1, "XIA");
    dpd_buf4_init(&XIjAb, EOM_XI, irrep, 0, 5, 0, 5, 0, "XIjAb");
    dpd_buf4_sort(&XIjAb, EOM_XI, pqsr, 0, 5, "XIjbA");
    dpd_buf4_init(&XIjbA, EOM_XI, irrep, 0, 5, 0, 5, 0, "XIjbA");
    
    norm = norm_C_rhf(&XIA, &XIjAb, &XIjbA);
    
    dpd_file2_close(&XIA);
    dpd_buf4_close(&XIjAb);
    dpd_buf4_close(&XIjbA);
  }
  else if (params.ref == 1) {
    dpd_file2_init(&XIA, EOM_XI, irrep, 0, 1, "XIA");
    dpd_file2_init(&Xia, EOM_XI, irrep, 0, 1, "Xia");
    dpd_buf4_init(&XIJAB, EOM_XI, irrep, 2, 7, 2, 7, 0, "XIJAB");
    dpd_buf4_init(&Xijab, EOM_XI, irrep, 2, 7, 2, 7, 0, "Xijab");
    dpd_buf4_init(&XIjAb, EOM_XI, irrep, 0, 5, 0, 5, 0, "XIjAb");

    c_clean(&XIA, &Xia, &XIJAB, &Xijab, &XIjAb);
    norm = norm_C(&XIA, &Xia, &XIJAB, &Xijab, &XIjAb);
    
    dpd_file2_close(&XIA);
    dpd_file2_close(&Xia);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_close(&Xijab);
    dpd_buf4_close(&XIjAb);
  }
  else if (params.ref == 2) {
    dpd_file2_init(&XIA, EOM_XI, irrep, 0, 1, "XIA");
    dpd_file2_init(&Xia, EOM_XI, irrep, 2, 3, "Xia");
    dpd_buf4_init(&XIJAB, EOM_XI, irrep, 2, 7, 2, 7, 0, "XIJAB");
    dpd_buf4_init(&Xijab, EOM_XI, irrep, 12, 17, 12, 17, 0, "Xijab");
    dpd_buf4_init(&XIjAb, EOM_XI, irrep, 22, 28, 22, 28, 0, "XIjAb");
    
    norm = norm_C(&XIA, &Xia, &XIJAB, &Xijab, &XIjAb);
    
    dpd_file2_close(&XIA);
    dpd_file2_close(&Xia);
    dpd_buf4_close(&XIJAB);
    dpd_buf4_close(&Xijab);
    dpd_buf4_close(&XIjAb);
  }

  fprintf(outfile,"%7s, D(norm sigma)=%15.10lf\n", term_lbl, norm - old_norm);
  fflush(outfile);
  old_norm = norm;
  return;
}

}} // namespace psi::ccdensity
