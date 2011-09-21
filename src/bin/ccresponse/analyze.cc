/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libdpd/dpd.h>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "Params.h"
#include "MOInfo.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

double **Build_R(void);
double **Build_U(void);

void analyze(const char *pert, int irrep, double omega)
{
  int nirreps, h, i, j, a, b, ij, ab, u, v;
  int position, num_div, tot1, tot2, nvir, nso, nocc;
  double width, max, min, value, value2;
  double *amp_array;
  double **tmp, **T2trans, **T1trans;
  FILE *efile;
  dpdbuf4 I, T2, D;
  dpdfile2 T1;
  char lbl[32];

  nirreps = moinfo.nirreps;
  num_div = 500;
  max = 9;
  min = 0;
  width = (max-min) / (num_div);


  sprintf(lbl, "X_%s_%5.3f", pert, omega);
  ffile(&efile, lbl, 1);
  amp_array = init_array(num_div);

  nvir = moinfo.virtpi[0];
  nocc = moinfo.occpi[0];
  nso = moinfo.nso;

  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  dpd_buf4_init(&T2, CC_LR, 0, 0, 5, 0, 5, 0, lbl);
  dpd_buf4_mat_irrep_init(&T2, 0);
  dpd_buf4_mat_irrep_rd(&T2, 0);
  T2trans = block_matrix(nocc*nocc, nso*nso);
  tmp = block_matrix(nvir, nso);
  tot1 = 0;
  tot2 = 0;
  for(ij=0; ij<T2.params->rowtot[0]; ij++) {

    C_DGEMM('n', 't', nvir, nso, nvir, 1.0, &(T2.matrix[0][ij][0]), nvir, 
	    &(moinfo.C[0][0][0]), nvir, 0.0, &(tmp[0][0]), nso);
    C_DGEMM('n', 'n', nso, nso, nvir, 1.0, &(moinfo.C[0][0][0]), nvir,
	    tmp[0], nso, 0.0, T2trans[ij], nso);

    for(ab=0; ab<nso*nso; ab++) {
      value = fabs(log10(fabs(T2trans[ij][ab])));
      tot2++;
      if ((value >= max) && (value <= (max+width))) {
	amp_array[num_div-1]++;
	tot1++;
      }
      else if ((value <= min) && (value >= (min-width))) {
	amp_array[0]++;
	tot1++;
      }
      else if ((value < max) && (value > min)) {
	position = (int) floor((value-min)/width);
	amp_array[position]++;
	tot1++;
      }
    }
  }
  dpd_buf4_mat_irrep_close(&T2, 0);
  dpd_buf4_close(&T2);
  free_block(tmp);
  free_block(T2trans);

  value2 = 0;
  for (i = num_div-1; i >= 0; i--) {
    value = amp_array[i] / tot1;
    value2 += value;
    fprintf(efile, "%10.5lf %lf\n", -((i)*width)-min, value);
  }
  free(amp_array);
  fprintf(outfile, "Total number of converged T2 amplitudes = %d\n", tot2);
  fprintf(outfile, "Number of T2 amplitudes in analysis= %d\n", tot1);
  fclose(efile);

  num_div = 40;
  max = 2;
  min = -5;
  width = (max-min) / (num_div);

  sprintf(lbl, "X1_%s_%5.3f", pert, omega);
  ffile(&efile, lbl, 1);
  amp_array = init_array(num_div);

  sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, lbl);
  dpd_file2_print(&T1, outfile);
  dpd_file2_mat_init(&T1);
  dpd_file2_mat_rd(&T1);

  /*
  T1trans = block_matrix(nocc, nso);

  C_DGEMM('n','t', nocc, nso, nvir, 1.0, &(T1.matrix[0][0][0]), nvir,
	  &(moinfo.C[0][0][0]), nvir, 0.0, &(T1trans[0][0]), nso);
  */

  tot1 = tot2 = 0;
  for(i=0; i < nocc; i++) {
    for(a=0; a < nso; a++) {
      /*      value = fabs(log10(fabs(T1trans[i][a]))); */
      value = log10(fabs(T1.matrix[0][i][a]));
      tot2++;
      if ((value >= max) && (value <= (max+width))) {
	amp_array[num_div-1]++;
	tot1++;
      }
      else if ((value <= min) && (value >= (min-width))) {
	amp_array[0]++;
	tot1++;
      }
      else if ((value < max) && (value > min)) {
	position = (int) floor((value-min)/width);
	amp_array[position]++;
	tot1++;
      }
    }
  }
  /*  free_block(T1trans); */

  dpd_file2_mat_close(&T1);
  dpd_file2_close(&T1);

  value2 = 0;
  for (i = num_div-1; i >= 0; i--) {
    value = amp_array[i] / tot1;
    value2 += value;
    fprintf(efile, "%10.5lf %lf\n", ((i)*width)-min, value);
  }

  free(amp_array);
  fclose(efile);

}

}} // namespace psi::ccresponse
