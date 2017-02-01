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

/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

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
void amp_write_L1(dpdfile2 *L1, int length, const char *label, std::string out);
void amp_write_L2(dpdbuf4 *L2, int length, const char *label, std::string out);

/* print largest elements in CC_LAMBDA */
void Lamp_write(struct L_Params L_params) {
  dpdfile2 L1;
  dpdbuf4 L2;
  int L_irr;
  L_irr = L_params.irrep;

  if(params.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    amp_write_L1(&L1, params.num_amps, "\n\tLargest LIA Amplitudes:\n", "outfile");
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    amp_write_L2(&L2, params.num_amps, "\n\tLargest LIjAb Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&L2);
  }
  else if(params.ref == 1) { /** ROHF **/
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    amp_write_L1(&L1, params.num_amps, "\n\tLargest LIA Amplitudes:\n", "outfile");
    global_dpd_->file2_close(&L1);

    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
    amp_write_L1(&L1, params.num_amps, "\n\tLargest Lia Amplitudes:\n", "outfile");
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    amp_write_L2(&L2, params.num_amps, "\n\tLargest LIJAB Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    amp_write_L2(&L2, params.num_amps, "\n\tLargest Lijab Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    amp_write_L2(&L2, params.num_amps, "\n\tLargest LIjAb Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    amp_write_L1(&L1, params.num_amps, "\n\tLargest LIA Amplitudes:\n", "outfile");
    global_dpd_->file2_close(&L1);
    global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");
    amp_write_L1(&L1, params.num_amps, "\n\tLargest Lia Amplitudes:\n", "outfile");
    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    amp_write_L2(&L2, params.num_amps, "\n\tLargest LIJAB Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    amp_write_L2(&L2, params.num_amps, "\n\tLargest Lijab Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    amp_write_L2(&L2, params.num_amps, "\n\tLargest LIjAb Amplitudes:\n", "outfile");
    global_dpd_->buf4_close(&L2);
  }
}

void amp_write_L1(dpdfile2 *L1, int length, const char *label, std::string out)
{
  int m, h, nirreps, Gia;
  int i, I, a, A, numt1;
  int num2print=0;
  double value;
  struct onestack *t1stack;

  nirreps = L1->params->nirreps;
  Gia = L1->my_irrep;

  t1stack = (struct onestack *) malloc(length * sizeof(struct onestack));
  for(m=0; m < length; m++) { t1stack[m].value = 0; t1stack[m].i = 0; t1stack[m].a = 0; }

  global_dpd_->file2_mat_init(L1);
  global_dpd_->file2_mat_rd(L1);

  numt1 = 0;
  for(h=0; h < nirreps; h++) {

    numt1 += L1->params->rowtot[h] * L1->params->coltot[h^Gia];

    for(i=0; i < L1->params->rowtot[h]; i++) {
      I = L1->params->roworb[h][i];
      for(a=0; a < L1->params->coltot[h^Gia]; a++) {
	A = L1->params->colorb[h^Gia][a];
	value = L1->matrix[h][i][a];
	for(m=0; m < length; m++) {
	  if((fabs(value) - fabs(t1stack[m].value)) > 1e-12) {
	    onestack_insert(t1stack, value, I, A, m, length);
	    break;
	  }
	}
      }
    }
  }

  global_dpd_->file2_mat_close(L1);

  for(m=0; m < ((numt1 < length) ? numt1 : length); m++)
    if(fabs(t1stack[m].value) > 1e-8) num2print++;

  if(num2print) outfile->Printf( "%s", label);

  for(m=0; m < ((numt1 < length) ? numt1 : length); m++)
    if(fabs(t1stack[m].value) > 1e-8)
      outfile->Printf( "\t        %3d %3d %20.10f\n", t1stack[m].i, t1stack[m].a, t1stack[m].value);

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

void amp_write_L2(dpdbuf4 *L2, int length, const char *label, std::string out)
{
  int m, h, nirreps, Gijab, numt2;
  int ij, ab, i, j, a, b;
  int num2print=0;
  double value;
  struct twostack *t2stack;

  nirreps = L2->params->nirreps;
  Gijab = L2->file.my_irrep;

  t2stack = (struct twostack *) malloc(length * sizeof(struct twostack));
  for(m=0; m < length; m++) {
    t2stack[m].value = 0;
    t2stack[m].i = 0; t2stack[m].j = 0;
    t2stack[m].a = 0; t2stack[m].b = 0;
  }

  numt2 = 0;
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(L2, h);
    global_dpd_->buf4_mat_irrep_rd(L2, h);

    numt2 += L2->params->rowtot[h] * L2->params->coltot[h^Gijab];

    for(ij=0; ij < L2->params->rowtot[h]; ij++) {
      i = L2->params->roworb[h][ij][0];
      j = L2->params->roworb[h][ij][1];
      for(ab=0; ab < L2->params->coltot[h^Gijab]; ab++) {
	a = L2->params->colorb[h^Gijab][ab][0];
	b = L2->params->colorb[h^Gijab][ab][1];

	value = L2->matrix[h][ij][ab];

	for(m=0; m < length; m++) {
	  if((fabs(value) - fabs(t2stack[m].value)) > 1e-12) {
	    twostack_insert(t2stack, value, i, j, a, b, m, length);
	    break;
	  }
	}
      }
    }

    global_dpd_->buf4_mat_irrep_close(L2, h);
  }

  for(m=0; m < ((numt2 < length) ? numt2 : length); m++)
    if(fabs(t2stack[m].value) > 1e-8) num2print++;

  if(num2print) outfile->Printf( "%s", label);

  for(m=0; m < ((numt2 < length) ? numt2 : length); m++)
    if(fabs(t2stack[m].value) > 1e-8)
      outfile->Printf( "\t%3d %3d %3d %3d %20.10f\n", t2stack[m].i, t2stack[m].j,
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


}} // namespace psi::cclambda
