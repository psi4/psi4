
#include <stdio.h>
#include <math.h>
#include "cscc.h"
#include "ga.h"

void output(double *z, int rowlow, int rowhi, int collow, int colhi, int rowdim, int coldim,
	    int nctl)    
{
  static int kcol = 8;
  static double zero = 0.0;

  int i, j, k;
  int last, begin;

  for (i = rowlow; i < rowhi; i++) {
    /*
      for (j = collow; j < colhi; j++) 
      {
      if (z[i,j] != zero) 
      {
      goto L15;
      }
      }
    */
    if (z[i] != zero) {
      goto L15;
    }
  }

  printf(" zero matrix \n");
  goto L3;

 L15:   
  if (rowhi < rowlow) {
    goto L3;
  }
  if (colhi < collow) {
    goto L3;
  }
  
  last = MIN(colhi, collow + kcol - 1);

  for (begin = collow; begin < colhi; begin += kcol) {
    for (i = begin; i < last; i++) {
      printf("            col %3d ",i);
    }
    printf("\n");
	
    for (k = rowlow; k < rowhi; k++) {
      //for (i = begin; i < last; i++) 
      //{
      //if (z[k,i] != zero) 
      if (z[k] != zero) {
	goto L5;
      }
      //}
      // goto L1;
      continue;
    L5:
      //for (i = begin; i < last; i++) 
      //{
      //printf("row %4d %9.4f",i,z[k,i]);
      printf("row %4d %9.4f", k, z[k]);
      //}
      printf("\n");
      //L1: 
      //	   ;
    }

    last = MIN((last + kcol), colhi);
  }

 L3:  printf("\n");
 
  return;
}
