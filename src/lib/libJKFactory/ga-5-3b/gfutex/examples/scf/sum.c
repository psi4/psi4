#include <stdlib.h>
#include <stdio.h>

#include <ga.h>

#include "sum.h"

double sum(int g_a)
{
  int me, ndims, *lo, *hi, *ld;
  double *dat, sum;
  int sz, i;

  me = GA_Nodeid();
  ndims = GA_Ndim(g_a);

  lo = malloc(sizeof(int) * ndims);
  hi = malloc(sizeof(int) * ndims);
  ld = malloc(sizeof(int) * ndims);

  NGA_Distribution(g_a, me, lo, hi);
  NGA_Access(g_a, lo, hi, &dat, ld);

  sz = 1;

  for (i = 0; i < ndims; i++) {    
    sz *= ((hi[i] - lo[i]) + 1);

    // printf("dim: %d, lo: %d, hi: %d\n", i, lo[i], hi[i]);
  }

  sum = 0.0;

  for (i = 0; i < sz; i++)
    sum += dat[i];

  NGA_Release(g_a, lo, hi);

  free(lo);
  free(hi);
  free(ld);

  return sum;
} /* sum */
