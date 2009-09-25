/*! 
  \file
  \ingroup CCENERGY
  \brief Write the amplitudes from ccenergy
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libdpd/dpd.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

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
void amp_write_T1(dpdfile2 *T1, int length, const char *label, FILE *outfile);
void amp_write_T2(dpdbuf4 *T2, int length, const char *label, FILE *outfile);

void amp_write(void)
{
  dpdfile2 T1;
  dpdbuf4 T2;

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    amp_write_T1(&T1, params.num_amps, "\n\tLargest TIA Amplitudes:\n", outfile);
    dpd_file2_close(&T1);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    amp_write_T2(&T2, params.num_amps, "\n\tLargest TIjAb Amplitudes:\n", outfile);
    dpd_buf4_close(&T2);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    amp_write_T1(&T1, params.num_amps, "\n\tLargest TIA Amplitudes:\n", outfile);
    dpd_file2_close(&T1);

    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    amp_write_T1(&T1, params.num_amps, "\n\tLargest Tia Amplitudes:\n", outfile);
    dpd_file2_close(&T1);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    amp_write_T2(&T2, params.num_amps, "\n\tLargest TIJAB Amplitudes:\n", outfile);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
    amp_write_T2(&T2, params.num_amps, "\n\tLargest Tijab Amplitudes:\n", outfile);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    amp_write_T2(&T2, params.num_amps, "\n\tLargest TIjAb Amplitudes:\n", outfile);
    dpd_buf4_close(&T2);
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    amp_write_T1(&T1, params.num_amps, "\n\tLargest TIA Amplitudes:\n", outfile);
    dpd_file2_close(&T1);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    amp_write_T1(&T1, params.num_amps, "\n\tLargest Tia Amplitudes:\n", outfile);
    dpd_file2_close(&T1);

    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    amp_write_T2(&T2, params.num_amps, "\n\tLargest TIJAB Amplitudes:\n", outfile);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
    amp_write_T2(&T2, params.num_amps, "\n\tLargest Tijab Amplitudes:\n", outfile);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    amp_write_T2(&T2, params.num_amps, "\n\tLargest TIjAb Amplitudes:\n", outfile);
    dpd_buf4_close(&T2);
  }
}

void amp_write_T1(dpdfile2 *T1, int length, const char *label, FILE *outfile)
{
  int m, h, nirreps, Gia;
  int i, I, a, A, numt1;
  int num2print=0;
  double value;
  struct onestack *t1stack;

  nirreps = T1->params->nirreps;
  Gia = T1->my_irrep;

  t1stack = (struct onestack *) malloc(length * sizeof(struct onestack));
  for(m=0; m < length; m++) { t1stack[m].value = 0; t1stack[m].i = 0; t1stack[m].a = 0; }

  dpd_file2_mat_init(T1);
  dpd_file2_mat_rd(T1);

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

  dpd_file2_mat_close(T1);

  for(m=0; m < ((numt1 < length) ? numt1 : length); m++)
    if(fabs(t1stack[m].value) > 1e-8) num2print++;

  if(num2print) fprintf(outfile, "%s", label);

  for(m=0; m < ((numt1 < length) ? numt1 : length); m++)
    if(fabs(t1stack[m].value) > 1e-8)
      fprintf(outfile, "\t        %3d %3d %20.10f\n", t1stack[m].i, t1stack[m].a, t1stack[m].value);

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

void amp_write_T2(dpdbuf4 *T2, int length, const char *label, FILE *outfile)
{
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
    dpd_buf4_mat_irrep_init(T2, h);
    dpd_buf4_mat_irrep_rd(T2, h);

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

    dpd_buf4_mat_irrep_close(T2, h);
  }

  for(m=0; m < ((numt2 < length) ? numt2 : length); m++)
    if(fabs(t2stack[m].value) > 1e-8) num2print++;

  if(num2print) fprintf(outfile, "%s", label);

  for(m=0; m < ((numt2 < length) ? numt2 : length); m++)
    if(fabs(t2stack[m].value) > 1e-8)
      fprintf(outfile, "\t%3d %3d %3d %3d %20.10f\n", t2stack[m].i, t2stack[m].j, 
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
