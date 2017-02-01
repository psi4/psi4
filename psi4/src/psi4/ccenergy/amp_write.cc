/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
  \file
  \ingroup CCENERGY
  \brief Write the amplitudes from ccenergy
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "Params.h"
#include "ccwave.h"
#include "psi4/libparallel/ParallelPrinter.h"

namespace psi { namespace ccenergy {

struct onestack {
    double value;
    int i;
    int a;
};

struct twostack {
    double value;
    int i; int j;
    int a; int b;
};

void onestack_insert(struct onestack *stack, double value, int i, int a,
    int level, int stacklen);
void twostack_insert(struct twostack *stack, double value, int i, int j,
    int a, int b, int level, int stacklen);
void amp_write_T1(dpdfile2 *T1, int length, const char *label, std::string OutFileRMR);
void amp_write_T2(dpdbuf4 *T2, int length, const char *label, std::string OutFileRMR);

void CCEnergyWavefunction::amp_write(void)
{
  dpdfile2 T1;
  dpdbuf4 T2;

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    amp_write_T1(&T1, params_.num_amps, "\n    Largest TIA Amplitudes:\n", "outfile");
    global_dpd_->file2_close(&T1);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    amp_write_T2(&T2, params_.num_amps, "\n    Largest TIjAb Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&T2);
  }
  else if(params_.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    amp_write_T1(&T1, params_.num_amps, "\n    Largest TIA Amplitudes:\n", "outfile");
    global_dpd_->file2_close(&T1);

    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tia");
    amp_write_T1(&T1, params_.num_amps, "\n    Largest Tia Amplitudes:\n", "outfile");
    global_dpd_->file2_close(&T1);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    amp_write_T2(&T2, params_.num_amps, "\n    Largest TIJAB Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
    amp_write_T2(&T2, params_.num_amps, "\n    Largest Tijab Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    amp_write_T2(&T2, params_.num_amps, "\n    Largest TIjAb Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&T2);
  }
  else if(params_.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    amp_write_T1(&T1, params_.num_amps, "\n    Largest TIA Amplitudes:\n", "outfile");
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    amp_write_T1(&T1, params_.num_amps, "\n    Largest Tia Amplitudes:\n", "outfile");
    global_dpd_->file2_close(&T1);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    amp_write_T2(&T2, params_.num_amps, "\n    Largest TIJAB Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
    amp_write_T2(&T2, params_.num_amps, "\n    Largest Tijab Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    amp_write_T2(&T2, params_.num_amps, "\n    Largest TIjAb Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&T2);
  }
}

void amp_write_T1(dpdfile2 *T1, int length, const char *label, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   int m, h, nirreps, Gia;
  int i, I, a, A, numt1;
  int num2print=0;
  double value;
  struct onestack *t1stack;

  nirreps = T1->params->nirreps;
  Gia = T1->my_irrep;

  t1stack = (struct onestack *) malloc(length * sizeof(struct onestack));
  for(m=0; m < length; m++) { t1stack[m].value = 0; t1stack[m].i = 0; t1stack[m].a = 0; }

  global_dpd_->file2_mat_init(T1);
  global_dpd_->file2_mat_rd(T1);

  numt1 = 0;
  for(h=0; h < nirreps; h++) {

    numt1 += T1->params->rowtot[h] * T1->params->coltot[h^Gia];

    for(i=0; i < T1->params->rowtot[h]; i++) {
      I = T1->params->roworb[h][i];
      for(a=0; a < T1->params->coltot[h^Gia]; a++) {
	A = T1->params->colorb[h][a];
	value = T1->matrix[h][i][a];
	for(m=0; m < length; m++) {
	  if((fabs(value) - fabs(t1stack[m].value)) > 1e-12) {
	    onestack_insert(t1stack, value, I, A, m, length);
	    break;
	  }
	}
      }
    }
  }

  global_dpd_->file2_mat_close(T1);

  for(m=0; m < ((numt1 < length) ? numt1 : length); m++)
    if(fabs(t1stack[m].value) > 1e-8) num2print++;

  if(num2print) printer->Printf( "%s", label);

  for(m=0; m < ((numt1 < length) ? numt1 : length); m++)
    if(fabs(t1stack[m].value) > 1e-8)
      printer->Printf( "            %3d %3d %20.10f\n", t1stack[m].i, t1stack[m].a, t1stack[m].value);

  free(t1stack);
}

void onestack_insert(struct onestack *stack, double value, int i, int a, int level, int stacklen)
{
  int l;
  struct onestack temp;

  temp = stack[level];

  stack[level].value = value;
  stack[level].i = i;
  stack[level].a = a;

  value = temp.value;
  i = temp.i;
  a = temp.a;

  for(l=level; l < stacklen-1; l++) {
    temp = stack[l+1];

    stack[l+1].value = value;
    stack[l+1].i = i;
    stack[l+1].a = a;

    value = temp.value;
    i = temp.i;
    a = temp.a;
  }
}

void amp_write_T2(dpdbuf4 *T2, int length, const char *label, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   int m, h, nirreps, Gijab, numt2;
  int ij, ab, i, j, a, b;
  int num2print=0;
  double value;
  struct twostack *t2stack;

  nirreps = T2->params->nirreps;
  Gijab = T2->file.my_irrep;

  t2stack = (struct twostack *) malloc(length * sizeof(struct twostack));
  for(m=0; m < length; m++) {
    t2stack[m].value = 0;
    t2stack[m].i = 0; t2stack[m].j = 0;
    t2stack[m].a = 0; t2stack[m].b = 0;
  }

  numt2 = 0;
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(T2, h);
    global_dpd_->buf4_mat_irrep_rd(T2, h);

    numt2 += T2->params->rowtot[h] * T2->params->coltot[h^Gijab];

    for(ij=0; ij < T2->params->rowtot[h]; ij++) {
      i = T2->params->roworb[h][ij][0];
      j = T2->params->roworb[h][ij][1];
      for(ab=0; ab < T2->params->coltot[h^Gijab]; ab++) {
	a = T2->params->colorb[h^Gijab][ab][0];
	b = T2->params->colorb[h^Gijab][ab][1];

	value = T2->matrix[h][ij][ab];

	for(m=0; m < length; m++) {
	  if((fabs(value) - fabs(t2stack[m].value)) > 1e-12) {
	    twostack_insert(t2stack, value, i, j, a, b, m, length);
	    break;
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_close(T2, h);
  }

  for(m=0; m < ((numt2 < length) ? numt2 : length); m++)
    if(fabs(t2stack[m].value) > 1e-8) num2print++;

  if(num2print) printer->Printf( "%s", label);

  for(m=0; m < ((numt2 < length) ? numt2 : length); m++)
    if(fabs(t2stack[m].value) > 1e-8)
      printer->Printf( "    %3d %3d %3d %3d %20.10f\n", t2stack[m].i, t2stack[m].j,
	      t2stack[m].a, t2stack[m].b, t2stack[m].value);

  free(t2stack);
}

void twostack_insert(struct twostack *stack, double value, int i, int j, int a, int b,
		     int level, int stacklen)
{
  int l;
  struct twostack temp;

  temp = stack[level];

  stack[level].value = value;
  stack[level].i = i;
  stack[level].j = j;
  stack[level].a = a;
  stack[level].b = b;

  value = temp.value;
  i = temp.i;
  j = temp.j;
  a = temp.a;
  b = temp.b;

  for(l=level; l < stacklen-1; l++) {
    temp = stack[l+1];

    stack[l+1].value = value;
    stack[l+1].i = i;
    stack[l+1].j = j;
    stack[l+1].a = a;
    stack[l+1].b = b;

    value = temp.value;
    i = temp.i;
    j = temp.j;
    a = temp.a;
    b = temp.b;
  }
}

}} // namespace psi::ccenergy
