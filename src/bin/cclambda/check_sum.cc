/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <cmath>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);

double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);

void check_sum(char *term_lbl, int irrep) {
  dpdfile2 Lia, LIA;
  dpdbuf4 LIJAB, Lijab, LIjAb, LIjbA;
  static double old_norm=0;
  double norm,dotval;
  char lbl[80];

  if (!strcmp(term_lbl,"reset"))  {
    fprintf(outfile,"resetting norm\n");
  	old_norm = 0;
  	return;
  }

  if (params.ref <= 1) {
    dpd_file2_init(&LIA, CC_LAMBDA, irrep, 0, 1, "New LIA");
    dpd_file2_init(&Lia, CC_LAMBDA, irrep, 0, 1, "New Lia");
    dpd_buf4_init(&LIJAB, CC_LAMBDA, irrep, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&Lijab, CC_LAMBDA, irrep, 2, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_init(&LIjAb, CC_LAMBDA, irrep, 0, 5, 0, 5, 0, "New LIjAb");

    norm = norm_C(&LIA, &Lia, &LIJAB, &Lijab, &LIjAb); 

    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_close(&Lijab);
    dpd_buf4_close(&LIjAb);
  }
  else if (params.ref == 2) {
    dpd_file2_init(&LIA, CC_LAMBDA, irrep, 0, 1, "New LIA");
    dpd_file2_init(&Lia, CC_LAMBDA, irrep, 2, 3, "New Lia");
    dpd_buf4_init(&LIJAB, CC_LAMBDA, irrep, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&Lijab, CC_LAMBDA, irrep, 12, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_init(&LIjAb, CC_LAMBDA, irrep, 22, 28, 22, 28, 0, "New LIjAb");

    norm = norm_C(&LIA, &Lia, &LIJAB, &Lijab, &LIjAb);

    dpd_file2_close(&LIA);
    dpd_file2_close(&Lia);
    dpd_buf4_close(&LIJAB);
    dpd_buf4_close(&Lijab);
    dpd_buf4_close(&LIjAb);
  }

  fprintf(outfile,"%7s, D(norm L)=%15.10lf\n", term_lbl, norm - old_norm);
  fflush(outfile);
  old_norm = norm;
  return;
}

double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
        dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf)
{
  double norm = 0.0;
  norm += dpd_file2_dot_self(CME);
  norm += dpd_file2_dot_self(Cme);
  norm += dpd_buf4_dot_self(CMNEF);
  norm += dpd_buf4_dot_self(Cmnef);
  norm += dpd_buf4_dot_self(CMnEf);
  norm = sqrt(norm);
  return norm;
}

double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE) {
  double norm = 0.0;
  norm = 2.0 * dpd_file2_dot_self(CME);
  norm += 2.0 * dpd_buf4_dot_self(CMnEf);
  norm -= dpd_buf4_dot(CMnEf, CMnfE);
  norm = sqrt(norm);
  return norm;
} 


}} // namespace psi::cclambda
