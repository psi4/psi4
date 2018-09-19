/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "ccwave.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libciomr/libciomr.h"
namespace psi { namespace ccenergy {


void CCEnergyWavefunction::analyze(void)
{
  int nirreps, h, i, j, a, b, ij, ab, u, v;
  int position, num_div, tot1, tot2, nvir, nso, nocc;
  double width, max, min, value, value2;
  double *amp_array;
  double **tmp, **T2trans, **T1trans;

  dpdfile2 T1;
  dpdbuf4 I, T2, D;

  nirreps = moinfo_.nirreps;
  num_div = 500;
  max = 9;
  min = 0;
  width = (max-min) / (num_div);
  auto printer = std::make_shared<PsiOutStream>("tamps.dat",std::ostream::app);
  amp_array = init_array(num_div);

  nvir = moinfo_.virtpi[0];
  nocc = moinfo_.occpi[0];
  nso = moinfo_.nso;

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->buf4_mat_irrep_init(&T2, 0);
  global_dpd_->buf4_mat_irrep_rd(&T2, 0);
  T2trans = block_matrix(nocc*nocc, nso*nso);
  tmp = block_matrix(nvir, nso);
  tot1 = 0;
  tot2 = 0;
  for(ij=0; ij<T2.params->rowtot[0]; ij++) {

    C_DGEMM('n', 't', nvir, nso, nvir, 1.0, &(T2.matrix[0][ij][0]), nvir,
        &(moinfo_.Cv[0][0][0]), nvir, 0.0, &(tmp[0][0]), nso);
    C_DGEMM('n', 'n', nso, nso, nvir, 1.0, &(moinfo_.Cv[0][0][0]), nvir,
	    tmp[0], nso, 0.0, T2trans[ij], nso);

    for(ab=0; ab<nso*nso; ab++) {
      value = std::fabs(std::log10(std::fabs(T2trans[ij][ab])));
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
  global_dpd_->buf4_mat_irrep_close(&T2, 0);
  global_dpd_->buf4_close(&T2);
  free_block(tmp);
  free_block(T2trans);

  value2 = 0;
  for (i = num_div-1; i >= 0; i--) {
    value = amp_array[i] / tot1;
    value2 += value;
    printer->Printf("%10.5lf %lf\n", -((i)*width)-min, value);
  }
  free(amp_array);
  printf("Total number of converged T2 amplitudes = %d\n", tot2);
  printf("Number of T2 amplitudes in analysis= %d\n", tot1);


  num_div = 40;
  max = 2;
  min = -5;
  width = (max-min) / (num_div);
  auto printer2 = std::make_shared<PsiOutStream>("t1amps.dat",std::ostream::app);
  amp_array = init_array(num_div);

  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_print(&T1, "outfile");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);
  /*
  T1trans = block_matrix(nocc, nso);

  C_DGEMM('n','t', nocc, nso, nvir, 1.0, &(T1.matrix[0][0][0]), nvir,
	  &(moinfo.C[0][0][0]), nvir, 0.0, &(T1trans[0][0]), nso);
  */

  tot1 = tot2 = 0;
  for(i=0; i < nocc; i++) {
    for(a=0; a < nso; a++) {
      /*      value = std::fabs(log10(std::fabs(T1trans[i][a]))); */
      value = std::log10(std::fabs(T1.matrix[0][i][a]));
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

  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);

  value2 = 0;
  for (i = num_div-1; i >= 0; i--) {
    value = amp_array[i] / tot1;
    value2 += value;
    printer2->Printf("%10.5lf %lf\n", ((i)*width)-min, value);
  }

  free(amp_array);

}

}} // namespace psi::ccenergy
