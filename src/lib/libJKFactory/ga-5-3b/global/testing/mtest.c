#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*File: mtest.c 

Anthor: Limin Zhang
purpose: testing those matrix functions developed for TAO/Global Array project.

Date: 2/28/2002

*/

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif

#include "ga.h"
#include "globalp.h"
#include "mp3.h"

#define N 100            /* dimension of matrices */
#define BLOCK_SIZE 32           /* dimension of blocks (if used) */

#define OP_SHIFT_DIAGONAL 1
#define OP_SET_DIAGONAL 2
#define OP_ADD_DIAGONAL         3
#define OP_GET_DIAGONAL 4
#define OP_NORM1    5
#define OP_NORM_INFINITY    6
#define OP_MEDIAN               7
#define OP_MEDIAN_PATCH         8
#define OP_SCALE_ROWS           9
#define OP_SCALE_COLS           10
#define OP_ZERO_DIAGONAL        11

# define THRESH 1e-5
#define MISMATCHED(x,y) GA_ABS((x)-(y))>=THRESH

/*#define BLOCK_CYCLIC*/
/*#define USE_SCALAPACK*/

#ifdef USE_SCALAPACK
#define BLOCK_CYCLIC
#endif

/* Verify that the diagonal values of a matrix have the value val and
   that the off-diagonal values are off_val */
void check_diag(int g_a, void *val, void *off_val, int *chk)
{
  int ai, bi;
  long al, bl;
  float af, bf;
  double ad, bd;
  int lo[2], hi[2];
  int i, j, ii, jj, me, ld, ioff;
  int ndim, dims[2], type;
  void *ptr;
  DoubleComplex adc, bdc;
  SingleComplex afc, bfc;

  me = GA_Nodeid();
  NGA_Inquire (g_a, &type, &ndim, dims);
  NGA_Distribution(g_a, me, lo, hi);
  NGA_Access(g_a, lo, hi, &ptr, &ld);
  *chk = 1;
  switch (type) {
    case C_INT:
      ai = *((int*)val);
      bi = *((int*)off_val);
      for (j=0; j<=hi[0]-lo[0]; j++) {
        jj = j + lo[0];
        for (i=0; i<=hi[1]-lo[1]; i++) {
          ii = i + lo[1];
          ioff = j*ld + i;
          if (ii == jj) {
            if (((int*)ptr)[ioff] != ai) *chk = 0;
          } else {
            if (((int*)ptr)[ioff] != bi) *chk = 0;
          }
        }
      }
      break;
    case C_DCPL:
      adc = *((DoubleComplex*)val);
      bdc = *((DoubleComplex*)off_val);
      for (j=0; j<=hi[0]-lo[0]; j++) {
        jj = j + lo[0];
        for (i=0; i<=hi[1]-lo[1]; i++) {
          ii = i + lo[1];
          ioff = j*ld + i;
          if (ii == jj) {
            if (((DoubleComplex*)ptr)[ioff].real != adc.real) *chk = 0;
            if (((DoubleComplex*)ptr)[ioff].imag != adc.imag) *chk = 0;
          } else {
            if (((DoubleComplex*)ptr)[ioff].real != bdc.real) *chk = 0;
            if (((DoubleComplex*)ptr)[ioff].imag != bdc.imag) *chk = 0;
          }
        }
      }
      break;
    case C_SCPL:
      afc = *((SingleComplex*)val);
      bfc = *((SingleComplex*)off_val);
      for (j=0; j<=hi[0]-lo[0]; j++) {
        jj = j + lo[0];
        for (i=0; i<=hi[1]-lo[1]; i++) {
          ii = i + lo[1];
          ioff = j*ld + i;
          if (ii == jj) {
            if (((SingleComplex*)ptr)[ioff].real != afc.real) *chk = 0;
            if (((SingleComplex*)ptr)[ioff].imag != afc.imag) *chk = 0;
          } else {
            if (((SingleComplex*)ptr)[ioff].real != bfc.real) *chk = 0;
            if (((SingleComplex*)ptr)[ioff].imag != bfc.imag) *chk = 0;
          }
        }
      }
      break;
    case C_DBL:
      ad = *((double*)val);
      bd = *((double*)off_val);
      for (j=0; j<=hi[0]-lo[0]; j++) {
        jj = j + lo[0];
        for (i=0; i<=hi[1]-lo[1]; i++) {
          ii = i + lo[1];
          ioff = j*ld + i;
          if (ii == jj) {
            if (((double*)ptr)[ioff] != ad) *chk = 0;
          } else {
            if (((double*)ptr)[ioff] != bd) *chk = 0;
          }
        }
      }
      break;
    case C_FLOAT:
      af = *((float*)val);
      bf = *((float*)off_val);
      for (j=0; j<=hi[0]-lo[0]; j++) {
        jj = j + lo[0];
        for (i=0; i<=hi[1]-lo[1]; i++) {
          ii = i + lo[1];
          ioff = j*ld + i;
          if (ii == jj) {
            if (((float*)ptr)[ioff] != af) *chk = 0;
          } else {
            if (((float*)ptr)[ioff] != bf) *chk = 0;
          }
        }
      }
      break;
    case C_LONG:
      al = *((long*)val);
      bl = *((long*)off_val);
      for (j=0; j<=hi[0]-lo[0]; j++) {
        jj = j + lo[0];
        for (i=0; i<=hi[1]-lo[1]; i++) {
          ii = i + lo[1];
          ioff = j*ld + i;
          if (ii == jj) {
            if (((long*)ptr)[ioff] != al) *chk = 0;
          } else {
            if (((long*)ptr)[ioff] != bl) *chk = 0;
          }
        }
      }
      break;
    default:
      GA_Error ("check diagonal:unknown data type.", type);
  }
  NGA_Release(g_a,lo,hi);
  GA_Igop(chk, 1, "*");
}

/* Verify that values in vector match target value val */
void check_vector(int g_v, void *val, int *chk)
{
  int ai;
  long al;
  float af;
  double ad;
  int lo, hi;
  int i, me, ld;
  int ndim, dims, type;
  void *ptr;
  DoubleComplex adc;
  SingleComplex afc;

  me = GA_Nodeid();
  NGA_Inquire (g_v, &type, &ndim, &dims);
  if (ndim != 1)
    GA_Error("Array is wrong dimension in check_vector: ",ndim);
  NGA_Distribution(g_v, me, &lo, &hi);
  NGA_Access(g_v, &lo, &hi, &ptr, &ld);
  *chk = 1;
  switch (type) {
    case C_INT:
      ai = *((int*)val);
      for (i = 0; i < hi-lo + 1; i++) {
        if (((int*)ptr)[i] != ai) *chk = 0;
      }
      break;
    case C_DCPL:
      adc = *((DoubleComplex*)val);
      for (i = 0; i < hi-lo + 1; i++) {
        if (((DoubleComplex*)ptr)[i].real != adc.real) *chk = 0;
        if (((DoubleComplex*)ptr)[i].imag != adc.imag) *chk = 0;
      }
      break;
    case C_SCPL:
      afc = *((SingleComplex*)val);
      for (i = 0; i < hi-lo + 1; i++) {
        if (((SingleComplex*)ptr)[i].real != afc.real) *chk = 0;
        if (((SingleComplex*)ptr)[i].imag != afc.imag) *chk = 0;
      }
      break;
    case C_DBL:
      ad = *((double*)val);
      for (i = 0; i < hi-lo + 1; i++) {
        if (((double*)ptr)[i] != ad) *chk = 0;
      }
      break;
    case C_FLOAT:
      af = *((float*)val);
      for (i = 0; i < hi-lo + 1; i++) {
        if (((float*)ptr)[i] != af) *chk = 0;
      }
      break;
    case C_LONG:
      al = *((long*)val);
      for (i = 0; i < hi-lo + 1; i++) {
        if (((long*)ptr)[i] != al) *chk = 0;
      }
      break;
    default:
      GA_Error ("check vector:unknown data type.", type);
  }
  NGA_Release(g_v,&lo,&hi);
  GA_Igop(chk, 1, "*");
}

void  test_scale_cols (int g_a, int g_v)
{

  int index[MAXDIM];
  void *min=NULL, *max=NULL;
  int imin, imax;
  float fmin, fmax;
  long lmin, lmax;
  double dmin, dmax;
  DoubleComplex dcmin, dcmax;
  SingleComplex fcmin, fcmax;


  void *alpha=NULL, *beta=NULL;
  int ai = 1, bi = -1;
  long al = 1, bl = -1;
  float af = 1.0, bf = -1.0;
  double ad = 1.0, bd = -1.0;
  DoubleComplex adc = { 1.0, 0.0 }, bdc =
    {
      -1.0, 0.0};

  SingleComplex afc = { 1.0, 0.0 }, bfc =
    {
      -1.0, 0.0};

  int g_b, g_c;
  int me = GA_Nodeid ();
  void *val=NULL;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval = { -2.0, 0.0 };
  SingleComplex fcval = { -2.0, 0.0 };
  void *val2=NULL;
  int ival2 = 4;
  double dval2 = 4.0;
  float fval2 = 4.0;
  long lval2 = 4;
  DoubleComplex dcval2 = { 4.0, 0.0 };
  SingleComplex fcval2 = { 4.0, 0.0 };

  int type, ndim, dims[MAXDIM];
  int vtype, vndim, vdims[MAXDIM];


  int n=0;
  /*
  int lo[2], hi[2], col, i, size;
  int vlo, vhi;
  void *buf;
  */



  NGA_Inquire (g_a, &type, &ndim, dims);
  NGA_Inquire (g_v, &vtype, &vndim, vdims);


  switch (type)
    {
    case C_INT:
      alpha = (void *) &ai;
      beta = (void *) &bi;
      break;
    case C_DCPL:
      alpha = (void *)&adc;
      beta = (void *)&bdc;
      break;

    case C_SCPL:
      alpha = (void *)&afc;
      beta = (void *)&bfc;
      break;

    case C_DBL:
      alpha = (void *) &ad;
      beta = (void *) &bd;
      break;
    case C_FLOAT:
      alpha = (void *) &af;
      beta =(void *)  &bf;
      break;
    case C_LONG:
      alpha = (void *) &al;
      beta = (void *) &bl;
      break;
    default:
      GA_Error ("test_scale_cols:wrong data type.", type);
    }

  switch (type)
    {
    case C_INT:
      val = (void *)&ival;
      val2 = (void *)&ival2;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      val2 = (void *)&dcval2;
      break;
    case C_SCPL:
      val = (void *)&fcval;
      val2 = (void *)&fcval2;
      break;
    case C_DBL:
      val = (void *)&dval;
      val2 = (void *)&dval2;
      break;
    case C_FLOAT:
      val = (void *)&fval;
      val2 = (void *)&fval2;
      break;
    case C_LONG:
      val = (void *)&lval;
      val2 = (void *)&lval2;
      break;
    default:
      GA_Error ("test_scale_cols:wrong data type.", type);
    }

  if (me == 0)
    printf ("Testing GA_Scale_cols...");


  /*
    n = vdims[0];
    printf("n=%d\n", n);
    size = GAsizeof (type);
    buf = malloc(n*size);
    if(me==0)printf("Initializing matrix A\n");
    vlo = 0;
    vhi = n-1;
    col = 0; 
    for(i=0; i<n; i++)
    switch(type){
    case C_INT:
    ((int*)buf)[i]=col+i;
    break;
    case C_LONG:
    ((long*)buf)[i]=(long)(col+i);
    break;
    case C_DBL:
    ((double*)buf)[i]=(double)(i+col);
    break;
    case C_FLOAT:
    ((float*)buf)[i]=(float)(col+i);
    break;
    case C_DCPL:
    ((DoubleComplex*)buf)[i].real=(double)(i+col);
    ((DoubleComplex*)buf)[i].imag=(double)0.0;
    break;
    case C_SCPL:
    ((SingleComplex*)buf)[i].real=(float)(i+col);
    ((SingleComplex*)buf)[i].imag=(float)0.0;
    break;
    default:
    GA_Error("test_scale:wrong data type.", type);
    }
    NGA_Put(g_v, &vlo, &vhi, buf, &n);


    n = dims[1];
    printf("dims[0]=%d \t dims[1]=%d\n", dims[0], dims[1]);
    printf("n=%d\n", n);
    size = GAsizeof (type);
    buf = malloc(n*size);
    if(me==0)printf("Initializing matrix A\n");
    lo[1]=0; hi[1]=n-1;
    for(col=me; col<dims[0]; col+= nproc){
    lo[0]=hi[0]=col;
    for(i=0; i<n; i++)
    switch(type){
    case C_INT:
    ((int*)buf)[i]=col+i;
    break;
    case C_LONG:
    ((long*)buf)[i]=(long)(col+i);
    break;
    case C_DBL:
    ((double*)buf)[i]=(double)(i+col);
    break;
    case C_FLOAT:
    ((float*)buf)[i]=(float)(col+i);
    break;
    case C_DCPL:
    ((DoubleComplex*)buf)[i].real=(double)(i+col);
    ((DoubleComplex*)buf)[i].imag=0.0;
    break;
    case C_SCPL:
    ((SingleComplex*)buf)[i].real=(float)(i+col);
    ((SingleComplex*)buf)[i].imag=0.0;
    break;
    default:
    GA_Error("test_scale_cols:wrong data type.", type);
    }
    NGA_Put(g_a, lo, hi, buf, &n);
    }

 
 
    free(buf);
    buf = NULL; 
  */

  GA_Fill (g_a, val);
  GA_Fill (g_v, val);
  GA_Scale_cols (g_a, g_v);
  /*the result is the same same as g_b filled with val2 */
  g_b = GA_Duplicate (g_a, "B");
  if (!g_b)
    GA_Error ("test_scale_cols:duplicate failed: B", n);

  g_c = GA_Duplicate (g_a, "C");
  if (!g_c)
    GA_Error ("duplicate failed: C", n);

  GA_Fill (g_b, val2);

  GA_Add (alpha, g_a, beta, g_b, g_c);

  switch (type)
    {
    case C_INT:
      max = (void *)&imax;
      min = (void *)&imin;
      break;
    case C_DCPL:
      max = (void *)&dcmax;
      min = (void *)&dcmin;
      break;
    case C_SCPL:
      max = (void *)&fcmax;
      min = (void *)&fcmin;
      break;
    case C_DBL:
      max = (void *)&dmax;
      min = (void *)&dmin;
      break;
    case C_FLOAT:
      max = (void *)&fmax;
      min = (void *)&fmin;
      break;
    case C_LONG:
      max = (void *)&lmax;
      min = (void *)&lmin;
      break;
    default:
      GA_Error ("test_scale_rows:wrong data type.", type);
    }

  NGA_Select_elem (g_c, "max", max, index);
  NGA_Select_elem (g_c, "min", min, index);


  switch (type)
    {
      double r, m;
      float rf, mf;
    case C_INT:
      if (me == 0)
    {
      if (MISMATCHED (imax, imin) || (imax != 0) || (imin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DCPL:
      r = dcmax.real - dcmin.real;
      m = dcmax.imag - dcmin.imag;
      if (me == 0)
    {
      if (MISMATCHED (dcmax.real, dcmin.real) || (dcmax.real != 0.0)
          || (dcmin.real != 0.0) || MISMATCHED (dcmax.imag, dcmin.imag)
          || (dcmax.imag != 0.0) || (dcmin.imag != 0.0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_SCPL:
      rf = fcmax.real - fcmin.real;
      mf = fcmax.imag - fcmin.imag;
      if (me == 0)
    {
      if (MISMATCHED (fcmax.real, fcmin.real) || (fcmax.real != 0.0)
          || (fcmin.real != 0.0) || MISMATCHED (fcmax.imag, fcmin.imag)
          || (fcmax.imag != 0.0) || (fcmin.imag != 0.0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DBL:
      if (me == 0)
    {
      if (MISMATCHED (dmax, dmin) || (dmax != 0) || (dmin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_FLOAT:
      if (me == 0)
    {
      if (MISMATCHED (fmax, fmin) || (fmax != 0) || (fmin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_LONG:
      if (me == 0)
    {
      if (MISMATCHED (lmax, lmin) || (lmax != 0) || (lmin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    default:
      GA_Error ("test_scale_rows:wrong data type.", type);
    }
}



void
test_scale_rows (int g_a, int g_v)
{

  int index[MAXDIM];
  void *min=NULL, *max=NULL;
  int imin, imax;
  float fmin, fmax;
  long lmin, lmax;
  double dmin, dmax;
  DoubleComplex dcmin, dcmax;
  SingleComplex fcmin, fcmax;


  void *alpha=NULL, *beta=NULL;
  int ai = 1, bi = -1;
  long al = 1, bl = -1;
  float af = 1.0, bf = -1.0;
  double ad = 1.0, bd = -1.0;
  DoubleComplex adc = { 1.0, 0.0 }, bdc =
    {
      -1.0, 0.0};

  SingleComplex afc = { 1.0, 0.0 }, bfc =
    {
      -1.0, 0.0};

  int g_b, g_c;
  int me = GA_Nodeid ();
  void *val=NULL;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval = { -2.0, 0.0 };
  SingleComplex fcval = { -2.0, 0.0 };
  void *val2=NULL;
  int ival2 = 4;
  double dval2 = 4.0;
  float fval2 = 4.0;
  long lval2 = 4;
  DoubleComplex dcval2 = { 4.0, 0.0 };
  SingleComplex fcval2 = { 4.0, 0.0 };

  int type, ndim, dims[MAXDIM];
  int vtype, vndim, vdims[MAXDIM];


  int n=0;
  /*
  int lo[2], hi[2], col, i, size;
  int vlo, vhi;
  void *buf;
  */



  NGA_Inquire (g_a, &type, &ndim, dims);
  NGA_Inquire (g_v, &vtype, &vndim, vdims);


  switch (type)
    {
    case C_INT:
      alpha = (void *) &ai;
      beta = (void *) &bi;
      break;
    case C_DCPL:
      alpha = (void *)&adc;
      beta = (void *)&bdc;
      break;

    case C_SCPL:
      alpha = (void *)&afc;
      beta = (void *)&bfc;
      break;

    case C_DBL:
      alpha = (void *) &ad;
      beta = (void *) &bd;
      break;
    case C_FLOAT:
      alpha = (void *)&af;
      beta = (void *)&bf;
      break;
    case C_LONG:
      alpha = (void *)&al;
      beta = (void *)&bl;
      break;
    default:
      GA_Error ("test_scale_rows:wrong data type.", type);
    }

  switch (type)
    {
    case C_INT:
      val = (void *) &ival;
      val2 = (void *) &ival2;
      break;
    case C_DCPL:
      val = (void *) &dcval;
      val2 = (void *) &dcval2;
      break;
    case C_SCPL:
      val = (void *) &fcval;
      val2 = (void *) &fcval2;
      break;
    case C_DBL:
      val = (void *) &dval;
      val2 = (void *) &dval2;
      break;
    case C_FLOAT:
      val = (void *) &fval;
      val2 = (void *) &fval2;
      break;
    case C_LONG:
      val = (void *) &lval;
      val2 = (void *) &lval2;
      break;
    default:
      GA_Error ("test_scale_rows:wrong data type.", type);
    }

  if (me == 0)
    printf ("Testing GA_Scale_rows...");


  /*
    n = vdims[0];
    printf("n=%d\n", n);
    size = GAsizeof (type);
    buf = malloc(n*size);
    if(me==0)printf("Initializing matrix A\n");
    vlo = 0;
    vhi = n-1;
    col = 0; 
    for(i=0; i<n; i++)
    switch(type){
    case C_INT:
    ((int*)buf)[i]=col+i;
    break;
    case C_LONG:
    ((long*)buf)[i]=(long)(col+i);
    break;
    case C_DBL:
    ((double*)buf)[i]=(double)(i+col);
    break;
    case C_FLOAT:
    ((float*)buf)[i]=(float)(col+i);
    break;
    case C_DCPL:
    ((DoubleComplex*)buf)[i].real=(double)(i+col);
    ((DoubleComplex*)buf)[i].imag=(double)0.0;
    break;
    case C_SCPL:
    ((SingleComplex*)buf)[i].real=(float)(i+col);
    ((SingleComplex*)buf)[i].imag=(float)0.0;
    break;
    default:
    GA_Error("test_scale:wrong data type.", type);
    }
    NGA_Put(g_v, &vlo, &vhi, buf, &n);


    n = dims[1];
    printf("dims[0]=%d \t dims[1]=%d\n", dims[0], dims[1]);
    printf("n=%d\n", n);
    size = GAsizeof (type);
    buf = malloc(n*size);
    if(me==0)printf("Initializing matrix A\n");
    lo[1]=0; hi[1]=n-1;
    for(col=me; col<dims[0]; col+= nproc){
    lo[0]=hi[0]=col;
    for(i=0; i<n; i++)
    switch(type){
    case C_INT:
    ((int*)buf)[i]=col+i;
    break;
    case C_LONG:
    ((long*)buf)[i]=(long)(col+i);
    break;
    case C_DBL:
    ((double*)buf)[i]=(double)(i+col);
    break;
    case C_FLOAT:
    ((float*)buf)[i]=(float)(col+i);
    break;
    case C_DCPL:
    ((DoubleComplex*)buf)[i].real=(double)(i+col);
    ((DoubleComplex*)buf)[i].imag=(double)0.0;
    break;
    case C_SCPL:
    ((SingleComplex*)buf)[i].real=(float)(i+col);
    ((SingleComplex*)buf)[i].imag=(float)0.0;
    break;
    default:
    GA_Error("test_scale:wrong data type.", type);
    }
    NGA_Put(g_a, lo, hi, buf, &n);
    }

 
 
    free(buf);
    buf = NULL; 
    GA_Print(g_v);

  */
  GA_Fill (g_a, val);
  GA_Fill (g_v, val);
  GA_Scale_rows (g_a, g_v);

  /*the result is the same same as g_b filled with val2 */
  g_b = GA_Duplicate (g_a, "B");
  if (!g_b)
    GA_Error ("duplicate failed: B", n);

  g_c = GA_Duplicate (g_a, "C");
  if (!g_c)
    GA_Error ("duplicate failed: C", n);

  GA_Fill (g_b, val2);

  GA_Add (alpha, g_a, beta, g_b, g_c);
  switch (type)
    {
    case C_INT:
      max = (void *) &imax;
      min = (void *) &imin;
      break;
    case C_DCPL:
      max = (void *) &dcmax;
      min = (void *) &dcmin;
      break;
    case C_SCPL:
      max = (void *) &fcmax;
      min = (void *) &fcmin;
      break;
    case C_DBL:
      max = (void *) &dmax;
      min = (void *) &dmin;
      break;
    case C_FLOAT:
      max = (void *) &fmax;
      min = (void *) &fmin;
      break;
    case C_LONG:
      max = (void *) &lmax;
      min = (void *) &lmin;
      break;
    default:
      GA_Error ("test_scale_rows:wrong data type.", type);
    }

  NGA_Select_elem (g_c, "max", max, index);
  NGA_Select_elem (g_c, "min", min, index);


  switch (type)
    {
      double r, m;
      float  rf, mf;
    case C_INT:
      if (me == 0)
    {
      if (MISMATCHED (imax, imin) || (imax != 0) || (imin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DCPL:
      r = dcmax.real - dcmin.real;
      m = dcmax.imag - dcmin.imag;
      if (me == 0)
    {
      if (MISMATCHED (dcmax.real, dcmin.real) || (dcmax.real != 0.0)
          || (dcmin.real != 0.0) || MISMATCHED (dcmax.imag, dcmin.imag)
          || (dcmax.imag != 0.0) || (dcmin.imag != 0.0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_SCPL:
      rf = fcmax.real - fcmin.real;
      mf = fcmax.imag - fcmin.imag;
      if (me == 0)
    {
      if (MISMATCHED (fcmax.real, fcmin.real) || (fcmax.real != 0.0)
          || (fcmin.real != 0.0) || MISMATCHED (fcmax.imag, fcmin.imag)
          || (fcmax.imag != 0.0) || (fcmin.imag != 0.0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DBL:
      if (me == 0)
    {
      if (MISMATCHED (dmax, dmin) || (dmax != 0) || (dmin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_FLOAT:
      if (me == 0)
    {
      if (MISMATCHED (fmax, fmin) || (fmax != 0) || (fmin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_LONG:
      if (me == 0)
    {
      if (MISMATCHED (lmax, lmin) || (lmax != 0) || (lmin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    default:
      GA_Error ("test_scale_rows:wrong data type.", type);
    }

}

void
test_median_patch (int g_a, int *alo, int *ahi, int g_b, int *blo, int *bhi,
           int g_c, int *clo, int *chi, int g_m, int *mlo, int *mhi)
{

  int g_e;
  int index[MAXDIM];
  void *min=NULL, *max=NULL;
  int imin, imax;
  float fmin, fmax;
  long lmin, lmax;
  double dmin, dmax;
  DoubleComplex dcmin, dcmax;
  SingleComplex fcmin, fcmax;


  void *alpha=NULL, *beta=NULL;
  int ai = 1, bi = -1;
  long al = 1, bl = -1;
  float af = 1.0, bf = -1.0;
  double ad = 1.0, bd = -1.0;
  DoubleComplex adc = { 1.0, 0.0 }, bdc =
    {
      -1.0, 0.0};
  SingleComplex afc = { 1.0, 0.0 }, bfc =
    {
      -1.0, 0.0};

  void *val=NULL;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;
  SingleComplex fcval;

  void *val2=NULL;
  int ival2 = 6;
  double dval2 = 6.0;
  float fval2 = 6.0;
  long lval2 = 6;
  DoubleComplex dcval2;
  SingleComplex fcval2;


  void *val3=NULL;
  int ival3 = 4;
  double dval3 = 4.0;
  float fval3 = 4.0;
  long lval3 = 4;
  DoubleComplex dcval3;
  SingleComplex fcval3;

  int me = GA_Nodeid ();
  int type, ndim, dims[MAXDIM];

  NGA_Inquire (g_a, &type, &ndim, dims);


  switch (type)
    {
    case C_INT:
      alpha = (void *) &ai;
      beta = (void *) &bi;
      break;
    case C_DCPL:
      alpha = (void *) &adc;
      beta = (void *) &bdc;
      break;

    case C_SCPL:
      alpha = (void *) &afc;
      beta = (void *) &bfc;
      break;

    case C_DBL:
      alpha = (void *) (void *) &ad;
      beta = (void *) (void *) &bd;
      break;
    case C_FLOAT:
      alpha = (void *) &af;
      beta = (void *) &bf;
      break;
    case C_LONG:
      alpha = (void *) &al;
      beta = (void *) &bl;
      break;
    default:
      GA_Error ("test_median:wrong data type.", type);
    }

  dcval.real = -2.0;
  dcval.imag = -0.0;

  dcval2.real = 6.0;
  dcval2.imag = 0.0;


  dcval3.real = 4.0;
  dcval3.imag = 0.0;

  fcval.real = -2.0;
  fcval.imag = -0.0;

  fcval2.real = 6.0;
  fcval2.imag = 0.0;


  fcval3.real = 4.0;
  fcval3.imag = 0.0;


  switch (type)
    {
    case C_INT:
      val = (void *) &ival;
      val2 = (void *) &ival2;
      val3 = (void *) &ival3;
      break;
    case C_DCPL:
      val = (void *) &dcval;
      val2 = (void *) &dcval2;
      val3 = (void *) &dcval3;
      break;
    case C_SCPL:
      val = (void *) &fcval;
      val2 = (void *) &fcval2;
      val3 = (void *) &fcval3;
      break;
    case C_DBL:
      val = (void *) &dval;
      val2 = (void *) &dval2;
      val3 = (void *) &dval3;
      break;
    case C_FLOAT:
      val = (void *) &fval;
      val2 = (void *) &fval2;
      val3 = (void *) &fval3;
      break;
    case C_LONG:
      val = (void *) &lval;
      val2 = (void *) &lval2;
      val3 = (void *) &lval3;
      break;
    default:
      GA_Error ("test_median:test_median:wrong data type.", type);
    }

  if (me == 0)
    printf ("Testing GA_Median_patch...");

  GA_Zero (g_a);
  GA_Zero (g_b);
  GA_Zero (g_c);
  GA_Zero (g_m);

  NGA_Fill_patch (g_a, alo, ahi, val);
  NGA_Fill_patch (g_b, blo, bhi, val2);
  NGA_Fill_patch (g_c, alo, bhi, val3);

  GA_Median_patch (g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi, g_m, mlo,
           mhi);

  /*
    The result array should        be g_c due to the value I chose: 
    val3 is the median of the three values val, val2, and val3
  */

  /* g_e = g_c - g_m */
  g_e = GA_Duplicate (g_a, "E");
  if (!g_e)
    GA_Error ("duplicate failed: E", 4);
  GA_Zero (g_e);

  NGA_Add_patch (alpha, g_c, clo, chi, beta, g_m, mlo, mhi, g_e, mlo, mhi);
  switch (type)
    {
    case C_INT:
      max = (void *)&imax;
      min = (void *)&imin;
      break;
    case C_DCPL:
      max = (void *)&dcmax;
      min = (void *)&dcmin;
      break;
    case C_SCPL:
      max = (void *)&fcmax;
      min = (void *)&fcmin;
      break;
    case C_DBL:
      max = (void *)&dmax;
      min = (void *)&dmin;
      break;
    case C_FLOAT:
      max = (void *)&fmax;
      min = (void *)&fmin;
      break;
    case C_LONG:
      max = (void *)&lmax;
      min = (void *)&lmin;
      break;
    default:
      GA_Error ("test_median:wrong data type.", type);
    }

  NGA_Select_elem (g_e, "max", max, index);
  NGA_Select_elem (g_e, "min", min, index);


  switch (type)
    {
      double r, m;
      float  rf, mf;
    case C_INT:
      if (me == 0)
    {
      if (MISMATCHED (imax, imin) || (imax != 0) || (imin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DCPL:
      r = dcmax.real - dcmin.real;
      m = dcmax.imag - dcmin.imag;
      if (me == 0)
    {
      if (MISMATCHED (dcmax.real, dcmin.real) || (dcmax.real != 0.0)
          || (dcmin.real != 0.0) || MISMATCHED (dcmax.imag, dcmin.imag)
          || (dcmax.imag != 0.0) || (dcmin.imag != 0.0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_SCPL:
      rf = fcmax.real - fcmin.real;
      mf = fcmax.imag - fcmin.imag;
      if (me == 0)
    {
      if (MISMATCHED (fcmax.real, fcmin.real) || (fcmax.real != 0.0)
          || (fcmin.real != 0.0) || MISMATCHED (fcmax.imag, fcmin.imag)
          || (fcmax.imag != 0.0) || (fcmin.imag != 0.0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DBL:
      if (me == 0)
    {
      if (MISMATCHED (dmax, dmin) || (dmax != 0) || (dmin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_FLOAT:
      if (me == 0)
    {
      if (MISMATCHED (fmax, fmin) || (fmax != 0) || (fmin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_LONG:
      if (me == 0)
    {
      if (MISMATCHED (lmax, lmin) || (lmax != 0) || (lmin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    default:
      GA_Error ("test_median:wrong data type.", type);
    }


}


void
test_median (int g_a, int g_b, int g_c, int g_m)
{

  int g_e;
  int index[MAXDIM];
  void *min=NULL, *max=NULL;
  int imin, imax;
  float fmin, fmax;
  long lmin, lmax;
  double dmin, dmax;
  DoubleComplex dcmin, dcmax;
  SingleComplex fcmin, fcmax;


  void *alpha=NULL, *beta=NULL;
  int ai = 1, bi = -1;
  long al = 1, bl = -1;
  float af = 1.0, bf = -1.0;
  double ad = 1.0, bd = -1.0;
  DoubleComplex adc = { 1.0, 0.0 }, bdc =
    {
      -1.0, 0.0};
  SingleComplex afc = { 1.0, 0.0 }, bfc =
    {
      -1.0, 0.0};

  void *val=NULL;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;
  SingleComplex fcval;

  void *val2=NULL;
  int ival2 = 6;
  double dval2 = 6.0;
  float fval2 = 6.0;
  long lval2 = 6;
  DoubleComplex dcval2;
  SingleComplex fcval2;


  void *val3=NULL;
  int ival3 = 4;
  double dval3 = 4.0;
  float fval3 = 4.0;
  long lval3 = 4;
  DoubleComplex dcval3;
  SingleComplex fcval3;

  int me = GA_Nodeid ();
  int type, ndim, dims[MAXDIM];

  NGA_Inquire (g_a, &type, &ndim, dims);

  switch (type)
    {
    case C_INT:
      alpha = (void *) &ai;
      beta = (void *) &bi;
      break;
    case C_DCPL:
      alpha = (void *)&adc;
      beta =(void *) &bdc;
      break;
    case C_SCPL:
      alpha = (void *)&afc;
      beta =(void *) &bfc;
      break;

    case C_DBL:
      alpha = (void *) &ad;
      beta = (void *) &bd;
      break;
    case C_FLOAT:
      alpha = (void *)&af;
      beta =(void *) &bf;
      break;
    case C_LONG:
      alpha = (void *)&al;
      beta = (void *)&bl;
      break;
    default:
      GA_Error ("test_median:wrong data type.", type);
    }

  dcval.real = -2.0;
  dcval.imag = -0.0;

  dcval2.real = 6.0;
  dcval2.imag = 0.0;


  dcval3.real = 4.0;
  dcval3.imag = 0.0;

  fcval.real = -2.0;
  fcval.imag = -0.0;

  fcval2.real = 6.0;
  fcval2.imag = 0.0;


  fcval3.real = 4.0;
  fcval3.imag = 0.0;


  switch (type)
    {
    case C_INT:
      val = (void *)&ival;
      val2 = (void *)&ival2;
      val3 = (void *)&ival3;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      val2 = (void *)&dcval2;
      val3 = (void *)&dcval3;
      break;
    case C_SCPL:
      val = (void *)&fcval;
      val2 = (void *)&fcval2;
      val3 = (void *)&fcval3;
      break;
    case C_DBL:
      val = (void *)&dval;
      val2 = (void *)&dval2;
      val3 = (void *)&dval3;
      break;
    case C_FLOAT:
      val = (void *)&fval;
      val2 = (void *)&fval2;
      val3 = (void *)&fval3;
      break;
    case C_LONG:
      val = (void *)&lval;
      val2 = (void *)&lval2;
      val3 = (void *)&lval3;
      break;
    default:
      GA_Error ("test_median:test_median:wrong data type.", type);
    }

  if (me == 0)
    printf ("Testing GA_Median...");

  GA_Zero (g_a);
  GA_Zero (g_b);
  GA_Zero (g_c);
  GA_Zero (g_m);

  GA_Fill (g_a, val);
  GA_Fill (g_b, val2);
  GA_Fill (g_c, val3);

  GA_Median (g_a, g_b, g_c, g_m);

  /*
  if (type == C_INT) {
    int achk[100000];
    int ii,jj,ndim,atype;
    int lo[2],hi[2],ld[2],dims[2];
    NGA_Inquire(g_a,&atype,&ndim,dims);
    for (ii=0; ii<ndim; ii++) {
      lo[ii] = 0;
      hi[ii] = dims[ii]-1;
      ld[ii] = dims[ii];
    }
    NGA_Get(g_a,lo,hi,achk,&ld[1]);
    if (me == 0) {
      for (ii=0; ii<dims[0]; ii++) {
        printf("\n");
        for (jj=0; jj<dims[1]; jj++) {
          printf("%8d",achk[ii*ld[1]+jj]);
        }
      }
      printf("\n");
    }
    printf("\n");
    NGA_Inquire(g_b,&atype,&ndim,dims);
    for (ii=0; ii<ndim; ii++) {
      lo[ii] = 0;
      hi[ii] = dims[ii]-1;
      ld[ii] = dims[ii];
    }
    NGA_Get(g_b,lo,hi,achk,&ld[1]);
    if (me == 0) {
      for (ii=0; ii<dims[0]; ii++) {
        printf("\n");
        for (jj=0; jj<dims[1]; jj++) {
          printf("%8d",achk[ii*ld[1]+jj]);
        }
      }
      printf("\n");
    }
    printf("\n");
    NGA_Inquire(g_c,&atype,&ndim,dims);
    for (ii=0; ii<ndim; ii++) {
      lo[ii] = 0;
      hi[ii] = dims[ii]-1;
      ld[ii] = dims[ii];
    }
    NGA_Get(g_c,lo,hi,achk,&ld[1]);
    if (me == 0) {
      for (ii=0; ii<dims[0]; ii++) {
        printf("\n");
        for (jj=0; jj<dims[1]; jj++) {
          printf("%8d",achk[ii*ld[1]+jj]);
        }
      }
      printf("\n");
    }
    printf("\n");
    NGA_Inquire(g_m,&atype,&ndim,dims);
    for (ii=0; ii<ndim; ii++) {
      lo[ii] = 0;
      hi[ii] = dims[ii]-1;
      ld[ii] = dims[ii];
    }
    NGA_Get(g_m,lo,hi,achk,&ld[1]);
    if (me == 0) {
      for (ii=0; ii<dims[0]; ii++) {
        printf("\n");
        for (jj=0; jj<dims[1]; jj++) {
          printf("%8d",achk[ii*ld[1]+jj]);
        }
      }
      printf("\n");
    }
  }
  */

  /*
    The result array should        be g_c due to the value I chose: 
    val3 is the median of the three values val, val2, and val3
  */

  /* g_e = g_c - g_m */
  g_e = GA_Duplicate (g_a, "E");
  if (!g_e)
    GA_Error ("duplicate failed: E", 4);

  GA_Add (alpha, g_c, beta, g_m, g_e);
  switch (type)
    {
    case C_INT:
      max = (void *)&imax;
      min = (void *)&imin;
      break;
    case C_DCPL:
      max = (void *)&dcmax;
      min = (void *)&dcmin;
      break;
    case C_SCPL:
      max = (void *)&fcmax;
      min = (void *)&fcmin;
      break;
    case C_DBL:
      max = (void *)&dmax;
      min = (void *)&dmin;
      break;
    case C_FLOAT:
      max = (void *)&fmax;
      min = (void *)&fmin;
      break;
    case C_LONG:
      max = (void *)&lmax;
      min = (void *)&lmin;
      break;
    default:
      GA_Error ("test_median:wrong data type.", type);
    }

  NGA_Select_elem (g_e, "max", max, index);
  NGA_Select_elem (g_e, "min", min, index);


  switch (type)
    {
      double r, m;
      float rf, mf;
    case C_INT:
      if (me == 0)
    {
      if (MISMATCHED (imax, imin) || (imax != 0) || (imin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DCPL:
      r = dcmax.real - dcmin.real;
      m = dcmax.imag - dcmin.imag;
      if (me == 0)
    {
      if (MISMATCHED (dcmax.real, dcmin.real) || (dcmax.real != 0.0)
          || (dcmin.real != 0.0) || MISMATCHED (dcmax.imag, dcmin.imag)
          || (dcmax.imag != 0.0) || (dcmin.imag != 0.0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_SCPL:
      rf = fcmax.real - fcmin.real;
      mf = fcmax.imag - fcmin.imag;
      if (me == 0)
    {
      if (MISMATCHED (fcmax.real, fcmin.real) || (fcmax.real != 0.0)
          || (fcmin.real != 0.0) || MISMATCHED (fcmax.imag, fcmin.imag)
          || (fcmax.imag != 0.0) || (fcmin.imag != 0.0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DBL:
      if (me == 0)
    {
      if (MISMATCHED (dmax, dmin) || (dmax != 0) || (dmin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_FLOAT:
      if (me == 0)
    {
      if (MISMATCHED (fmax, fmin) || (fmax != 0) || (fmin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_LONG:
      if (me == 0)
    {
      if (MISMATCHED (lmax, lmin) || (lmax != 0) || (lmin != 0))
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    default:
      GA_Error ("test_median:wrong data type.", type);
    }


}


void
test_norm_infinity (int g_a)
{

  void *val=NULL;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;
  SingleComplex fcval;

  double norm_infinity = -1.0, result = -1.0;

  int me = GA_Nodeid ();
  int type, ndim, dims[MAXDIM];

  NGA_Inquire (g_a, &type, &ndim, dims);
  dcval.real = -2.0;
  dcval.imag = -0.0;

  fcval.real = -2.0;
  fcval.imag = -0.0;

  switch (type)
    {
    case C_INT:
      val = (void *)&ival;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      break;
    case C_SCPL:
      val = (void *)&fcval;
      break;
    case C_DBL:
      val = (void *)&dval;
      break;
    case C_FLOAT:
      val = (void *)&fval;
      break;
    case C_LONG:
      val = (void *)&lval;
      break;
    default:
      GA_Error ("test_norm_infinity:wrong data type.", type);
    }

  if (me == 0)
    printf ("Testing GA_Norm_infinity...");
  GA_Fill (g_a, val);
  GA_Norm_infinity (g_a, &norm_infinity);
  /* GA_Print(g_a);
     printf("norm_infinity = %lf\n",norm_infinity);*/
  switch (type)
    {
    case C_INT:
      result = (double) GA_ABS (ival);
      break;
    case C_LONG:
      result = (double) GA_ABS (lval);
      break;
    case C_FLOAT:
      result = (double) GA_ABS (fval);
      break;

    case C_DBL:
      result = GA_ABS (dval);
      break;

    case C_DCPL:
      result = sqrt (dcval.real * dcval.real + dcval.imag * dcval.imag);
      break;
    case C_SCPL:
      result = sqrt (fcval.real * fcval.real + fcval.imag * fcval.imag);
      break;
    default:
      GA_Error ("test_norm_infinity: wrong data type.\n", type);
    }
  result = result * dims[0];
  if (me == 0)
    {
      if (MISMATCHED (result, norm_infinity))
    printf ("not ok.\n");
      else
    printf ("ok.\n");
    }
}


void
test_norm1 (int g_a)
{

  void *val=NULL;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;
  SingleComplex fcval;

  double norm1 = 0.0, result = -1.0;

  int me = GA_Nodeid ();
  int type, ndim, dims[MAXDIM];

  NGA_Inquire (g_a, &type, &ndim, dims);
  dcval.real = -2.0;
  dcval.imag = -0.0;

  fcval.real = -2.0;
  fcval.imag = -0.0;

  switch (type)
    {
    case C_INT:
      val = (void *) &ival;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      break;
    case C_SCPL:
      val = (void *)&fcval;
      break;
    case C_DBL:
      val = (void *)&dval;
      break;
    case C_FLOAT:
      val = (void *)&fval;
      break;
    case C_LONG:
      val = (void *)&lval;
      break;
    default:
      GA_Error ("test_norm1:wrong data type.", type);
    }

  if (me == 0)
    printf ("Testing GA_Norm1...");
  GA_Fill (g_a, val);
  GA_Norm1 (g_a, &norm1);
  /* GA_Print(g_a);
     printf("norm1=%lf\n", norm1);*/
  switch (type)
    {
    case C_INT:
      result = (double) GA_ABS (ival);
      break;
    case C_LONG:
      result = (double) GA_ABS (lval);
      break;
    case C_FLOAT:
      result = (double) GA_ABS (fval);
      break;

    case C_DBL:
      result = GA_ABS (dval);
      break;

    case C_DCPL:
      result = sqrt (dcval.real * dcval.real + dcval.imag * dcval.imag);
      break;
    case C_SCPL:
      result = sqrt (fcval.real * fcval.real + fcval.imag * fcval.imag);
      break;
    default:
      GA_Error ("test_norm1: wrong data type.\n", type);
    }
  result = result * dims[1];
  if (me == 0)
    {
      if (MISMATCHED (result, norm1))
    printf ("not ok.\n");
      else
    printf ("ok.\n");
    }
}


void
test_get_diagonal (int g_a, int g_v)
{

  int me = GA_Nodeid ();
  void *val=NULL;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;
  SingleComplex fcval;

  int type, ndim, dims[MAXDIM];
  int vtype, vndim, vdims[MAXDIM];
  int chk;

  NGA_Inquire (g_a, &type, &ndim, dims);
  NGA_Inquire (g_v, &vtype, &vndim, vdims);
  dcval.real = -2.0;
  dcval.imag = -0.0;

  fcval.real = -2.0;
  fcval.imag = -0.0;

  switch (type)
    {
    case C_INT:
      val = (void *)&ival;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      break;
    case C_SCPL:
      val = (void *)&fcval;
      break;
    case C_DBL:
      val = (void *)&dval;
      break;
    case C_FLOAT:
      val = (void *)&fval;
      break;
    case C_LONG:
      val = (void *)&lval;
      break;
    default:
      GA_Error ("test_get_diagonal:wrong data type.", type);
    }

  if (me == 0)
    printf ("Testing GA_Get_diag...");
  GA_Zero (g_v);
  GA_Fill (g_a, val);
  GA_Get_diag (g_a, g_v);
  switch (type)
    {
    case C_INT:
      check_vector(g_v, (void*)&ival, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_LONG:
      check_vector(g_v, (void*)&lval, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_FLOAT:
      check_vector(g_v, (void*)&fval, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DBL:
      check_vector(g_v, (void*)&dval, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DCPL:
      check_vector(g_v, (void*)&dcval, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_SCPL:
      check_vector(g_v, (void*)&fcval, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    default:
      GA_Error ("test_get_diagonal:wrong data type:", type);
    }
}


void
test_add_diagonal (int g_a, int g_v)
{


  int me = GA_Nodeid ();
  void *val=NULL;
  int ival = -2, izero = 0;
  double dval = -2.0, dzero = 0.0;
  float fval = -2.0, fzero = 0.0;
  long lval = -2, lzero = 0;
  DoubleComplex dcval, dczero = {0.0,0.0};
  SingleComplex fcval, fczero = {0.0,0.0};

  int type, ndim, dims[MAXDIM];
  int vtype, vndim, vdims[MAXDIM];
  int chk;


  NGA_Inquire (g_a, &type, &ndim, dims);
  dcval.real = -2.0;
  dcval.imag = -0.0;

  fcval.real = -2.0;
  fcval.imag = -0.0;

  NGA_Inquire (g_v, &vtype, &vndim, vdims);

  switch (type)
    {
    case C_INT:
      val = (void *)&ival;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      break;
    case C_SCPL:
      val = (void *)&fcval;
      break;
    case C_DBL:
      val =(void *) &dval;
      break;
    case C_FLOAT:
      val = (void *)&fval;
      break;
    case C_LONG:
      val = (void *)&lval;
      break;
    default:
      GA_Error ("test_add_diagonal:wrong data type.", type);
    }


  if (me == 0)
    printf ("Testing GA_Add_diagonal...");
  GA_Zero (g_a);
  GA_Fill (g_v, val);
  GA_Set_diagonal (g_a, g_v);

  /*reassign value to val */
  ival = 3;
  dval = 3.0;
  fval = 3.0;
  lval = 3;
  dcval.real = 3.0;
  dcval.imag = -0.0;
  fcval.real = 3.0;
  fcval.imag = -0.0;

  /*refile the global array g_v */
  GA_Fill (g_v, val);

  /*Add g_v to the diagonal of g_a */
  GA_Add_diagonal (g_a, g_v);
  /*after this line, the g_a should only have 1 on the diagonal and zeros every where else */

  switch (type)
    {
    case C_INT:
      ival = 1;
      check_diag(g_a, (void*)&ival, (void*)&izero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_LONG:
      lval = 1;
      check_diag(g_a, (void*)&lval, (void*)&lzero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_FLOAT:
      fval = 1.0;
      check_diag(g_a, (void*)&fval, (void*)&fzero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DBL:
      dval = 1.0;
      check_diag(g_a, (void*)&dval, (void*)&dzero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DCPL:
      dcval.real = 1.0;
      dcval.imag = 0.0;
      check_diag(g_a, (void*)&dcval, (void*)&dczero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_SCPL:
      fcval.real = 1.0;
      fcval.imag = 0.0;
      check_diag(g_a, (void*)&fcval, (void*)&fczero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    default:
      GA_Error ("test_add_diagonal:wrong data type:", type);
    }





}

void
test_set_diagonal (int g_a, int g_v)
{


  int me = GA_Nodeid ();
  void *val=NULL;
  int ival = -2, izero = 0;
  double dval = -2.0, dzero = 0.0;
  float fval = -2.0, fzero = 0.0;
  long lval = -2, lzero = 0;
  DoubleComplex dcval, dczero = {0.0,0.0};
  SingleComplex fcval, fczero = {0.0,0.0};

  int type, ndim, dims[MAXDIM];
  int vtype, vndim, vdims[MAXDIM];
  int chk;

  NGA_Inquire (g_a, &type, &ndim, dims);
  NGA_Inquire (g_v, &vtype, &vndim, vdims);
  dcval.real = -2.0;
  dcval.imag = -0.0;

  fcval.real = -2.0;
  fcval.imag = -0.0;

  switch (type)
    {
    case C_INT:
      val = (void *)&ival;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      break;
    case C_SCPL:
      val = (void *)&fcval;
      break;
    case C_DBL:
      val = (void *)&dval;
      break;
    case C_FLOAT:
      val = (void *)&fval;
      break;
    case C_LONG:
      val = (void *)&lval;
      break;
    default:
      GA_Error ("test_set_diagonal:wrong data type.", type);
    }


  if (me == 0)
    printf ("Testing GA_Set_diagonal...");
  GA_Zero (g_a);
  GA_Fill (g_v, val);
  GA_Set_diagonal (g_a, g_v);
  switch (type)
    {
    case C_INT:
      check_diag(g_a, (void*)&ival, (void*)&izero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_LONG:
      check_diag(g_a, (void*)&lval, (void*)&lzero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_FLOAT:
      check_diag(g_a, (void*)&fval, (void*)&fzero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DBL:
      check_diag(g_a, (void*)&dval, (void*)&dzero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DCPL:
      check_diag(g_a, (void*)&dcval, (void*)&dczero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_SCPL:
      check_diag(g_a, (void*)&fcval, (void*)&fczero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    default:
      GA_Error ("test_set_diagonal:wrong data type:", type);
    }

}

void
test_zero_diagonal (int g_a)
{


  int me = GA_Nodeid ();
  void *val=NULL;
  int ival = -2, izero = 0;
  double dval = -2.0, dzero = 0.0;
  float fval = -2.0, fzero = 0.0;
  long lval = -2, lzero = 0;
  DoubleComplex dcval, dczero = {0.0,0.0};
  SingleComplex fcval, fczero = {0.0,0.0};

  int g_b;
  int ialpha, ibeta;
  long lalpha, lbeta;
  float falpha, fbeta;
  double dalpha, dbeta;
  DoubleComplex zalpha, zbeta;
  SingleComplex calpha, cbeta;
  int vdims;
  void *alpha, *beta;

  int type, ndim, dims[MAXDIM];
  int chk;

  NGA_Inquire (g_a, &type, &ndim, dims);
  dcval.real = -2.0;
  dcval.imag = -0.0;

  fcval.real = -2.0;
  fcval.imag = -0.0;

  vdims = GA_MIN(dims[0],dims[1]);

  switch (type)
  {
    case C_INT:
      val = (void *)&ival;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      break;
    case C_SCPL:
      val = (void *)&fcval;
      break;
    case C_DBL:
      val = (void *)&dval;
      break;
    case C_FLOAT:
      val = (void *)&fval;
      break;
    case C_LONG:
      val = (void *)&lval;
      break;
    default:
      GA_Error ("test_zero_diagonal:wrong data type.", type);
  }


  if (me == 0)
    printf ("Testing GA_Zero_diagonal...");
  GA_Fill (g_a, val);
  GA_Zero_diagonal (g_a);
  g_b = GA_Duplicate(g_a, "tmp_array");
  GA_Fill (g_b, val);

  switch (type)
  {
    case C_INT:
      ialpha = -1;
      ibeta = 1;
      alpha = (void*)&ialpha;
      beta = (void*)&ibeta;
      GA_Add(alpha,g_a,beta,g_b,g_a);
      check_diag(g_a, (void*)&ival, (void*)&izero, &chk);
      if (me == 0)
      {
        if (chk != 1)
          printf ("not ok.\n");
        else
          printf ("ok.\n");
      }
      break;
    case C_LONG:
      lalpha = -1;
      lbeta = 1;
      alpha = (void*)&lalpha;
      beta = (void*)&lbeta;
      GA_Add(alpha,g_a,beta,g_b,g_a);
      check_diag(g_a, (void*)&lval, (void*)&lzero, &chk);
      if (me == 0)
      {
        if (chk != 1)
          printf ("not ok.\n");
        else
          printf ("ok.\n");
      }
      break;
    case C_FLOAT:
      falpha = -1.0;
      fbeta = 1.0;
      alpha = (void*)&falpha;
      beta = (void*)&fbeta;
      GA_Add(alpha,g_a,beta,g_b,g_a);
      check_diag(g_a, (void*)&fval, (void*)&fzero, &chk);
      if (me == 0)
      {
        if (chk != 1)
          printf ("not ok.\n");
        else
          printf ("ok.\n");
      }
      break;
    case C_DBL:
      dalpha = -1.0;
      dbeta = 1.0;
      alpha = (void*)&dalpha;
      beta = (void*)&dbeta;
      GA_Add(alpha,g_a,beta,g_b,g_a);
      check_diag(g_a, (void*)&dval, (void*)&dzero, &chk);
      if (me == 0)
      {
        if (chk != 1)
          printf ("not ok.\n");
        else
          printf ("ok.\n");
      }
      break;
    case C_DCPL:
      zalpha.real = -1.0;
      zalpha.imag = 0.0;
      zbeta.real = 1.0;
      zbeta.imag = 0.0;
      alpha = (void*)&zalpha;
      beta = (void*)&zbeta;
      GA_Add(alpha,g_a,beta,g_b,g_a);
      check_diag(g_a, (void*)&dcval, (void*)&dczero, &chk);
      if (me == 0)
      {
        if (chk != 1)
          printf ("not ok.\n");
        else
          printf ("ok.\n");
      }
      break;
    case C_SCPL:
      calpha.real = -1.0;
      calpha.imag = 0.0;
      cbeta.real = 1.0;
      cbeta.imag = 0.0;
      alpha = (void*)&calpha;
      beta = (void*)&cbeta;
      GA_Add(alpha,g_a,beta,g_b,g_a);
      check_diag(g_a, (void*)&fcval, (void*)&fczero, &chk);
      if (me == 0)
      {
        if (chk != 1)
          printf ("not ok.\n");
        else
          printf ("ok.\n");
      }
      break;
    default:
      GA_Error ("test_zero_diagonal:wrong data type:", type);
  }
  GA_Destroy(g_b);
}

void
test_shift_diagonal (int g_a)
{

  int me = GA_Nodeid ();
  void *val=NULL;
  int ival = -2, izero = 0;
  double dval = -2.0, dzero = 0.0;
  float fval = -2.0, fzero = 0.0;
  long lval = -2, lzero = 0;
  DoubleComplex dcval, dczero = {0.0,0.0};
  SingleComplex fcval, fczero = {0.0,0.0};

  int type, ndim, dims[MAXDIM];
  int dim;            /*the length of the diagonal */
  int chk;


  NGA_Inquire (g_a, &type, &ndim, dims);

  dim = GA_MIN (dims[0], dims[1]);

  dcval.real = -2.0;
  dcval.imag = -0.0;
  fcval.real = -2.0;
  fcval.imag = -0.0;
  switch (type)
    {
    case C_INT:
      val = (void *)&ival;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      break;
    case C_SCPL:
      val = (void *)&fcval;
      break;
    case C_DBL:
      val = (void *)&dval;
      break;
    case C_FLOAT:
      val = (void *)&fval;
      break;
    case C_LONG:
      val = (void *)&lval;
      break;
    default:
      GA_Error ("test_shift_diagonal:wrong data type.", type);
    }


  if (me == 0)
     printf ("Testing GA_Shift_diagonal...");
  
  GA_Zero (g_a);
  GA_Shift_diagonal (g_a, val);

  switch (type)
    {
    case C_INT:
      check_diag(g_a, (void*)&ival, (void*)&izero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_LONG:
      check_diag(g_a, (void*)&lval, (void*)&lzero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_FLOAT:
      check_diag(g_a, (void*)&fval, (void*)&fzero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DBL:
      check_diag(g_a, (void*)&dval, (void*)&dzero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_DCPL:
      check_diag(g_a, (void*)&dcval, (void*)&dczero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    case C_SCPL:
      check_diag(g_a, (void*)&fcval, (void*)&fczero, &chk);
      if (me == 0)
    {
      if (chk != 1)
        printf ("not ok.\n");
      else
        printf ("ok.\n");
    }
      break;
    default:
      GA_Error ("test_shift_diagonal:wrong data type: ", type);
    }



}

void
do_work (int type, int op)
{
  int g_a, g_b, g_c, g_m, g_v;
  int n = N;
  int nproc = GA_Nnodes ();
  int dims[2] = { N,        /*N columns */
          N + 2            /*N+2 rows */
  };
  int vdim, i;
  int lo[2], hi[2];


  int block_size[2], proc_grid[2], proc_cnt;

  proc_cnt = 1;
  for (i = 0; i<2; i++) {
    block_size[i] = BLOCK_SIZE;
  }

  if (nproc%2 == 0) {
    proc_grid[0] = 2;
    proc_grid[1] = nproc/2;
  } else {
    proc_grid[0] = nproc;
    proc_grid[1] = 1;
  }

  lo[0] = 1;
  hi[0] = dims[0] - 1;
  lo[1] = 1;
  hi[1] = dims[1] - 1;


  switch (op)
    {

    case OP_SHIFT_DIAGONAL:
#ifndef BLOCK_CYCLIC
      g_a = NGA_Create (type, 2, dims, "A", NULL);
      if (!g_a)
        GA_Error ("create failed: A", n);
#else
      g_a = GA_Create_handle();
      GA_Set_data(g_a, 2, dims, type);
      GA_Set_array_name(g_a,"A");
#ifdef USE_SCALAPACK
      GA_Set_block_cyclic_proc_grid(g_a, block_size, proc_grid);
#else
      GA_Set_block_cyclic(g_a, block_size);
#endif
      GA_Allocate(g_a);
#endif
      test_shift_diagonal (g_a);
      GA_Destroy (g_a);
      break;
    case OP_SET_DIAGONAL:
#ifndef BLOCK_CYCLIC
      g_a = NGA_Create (type, 2, dims, "A", NULL);
      if (!g_a)
        GA_Error ("create failed: A", n);
#else
      g_a = GA_Create_handle();
      GA_Set_data(g_a, 2, dims, type);
      GA_Set_array_name(g_a,"A");
#ifdef USE_SCALAPACK
      GA_Set_block_cyclic_proc_grid(g_a, block_size, proc_grid);
#else
      GA_Set_block_cyclic(g_a, block_size);
#endif
      GA_Allocate(g_a);
#endif
      /*find out the diagonal length of the matrix A */
      vdim = GA_MIN (dims[0], dims[1]);
      g_v = NGA_Create (type, 1, &vdim, "V", NULL);
      if (!g_v)
        GA_Error ("create failed:V", n);
      test_set_diagonal (g_a, g_v);
      GA_Destroy (g_a);
      GA_Destroy (g_v);
      break;
    case OP_ZERO_DIAGONAL:
#ifndef BLOCK_CYCLIC
      g_a = NGA_Create (type, 2, dims, "A", NULL);
      if (!g_a)
        GA_Error ("create failed: A", n);
#else
      g_a = GA_Create_handle();
      GA_Set_data(g_a, 2, dims, type);
      GA_Set_array_name(g_a,"A");
#ifdef USE_SCALAPACK
      GA_Set_block_cyclic_proc_grid(g_a, block_size, proc_grid);
#else
      GA_Set_block_cyclic(g_a, block_size);
#endif
      GA_Allocate(g_a);
#endif
      test_zero_diagonal (g_a);
      GA_Destroy (g_a);
      break;
    case OP_ADD_DIAGONAL:
#ifndef BLOCK_CYCLIC
      g_a = NGA_Create (type, 2, dims, "A", NULL);
      if (!g_a)
        GA_Error ("create failed: A", n);
#else
      g_a = GA_Create_handle();
      GA_Set_data(g_a, 2, dims, type);
      GA_Set_array_name(g_a,"A");
#ifdef USE_SCALAPACK
      GA_Set_block_cyclic_proc_grid(g_a, block_size, proc_grid);
#else
      GA_Set_block_cyclic(g_a, block_size);
#endif
      GA_Allocate(g_a);
#endif
      /*find out the diagonal length of the matrix A */
      vdim = GA_MIN (dims[0], dims[1]);
      g_v = NGA_Create (type, 1, &vdim, "V", NULL);
      if (!g_v)
        GA_Error ("create failed:V", n);
      vdim = GA_MIN (dims[0], dims[1]);
      test_add_diagonal (g_a, g_v);
      GA_Destroy (g_a);
      GA_Destroy (g_v);
      break;
    case OP_GET_DIAGONAL:
#ifndef BLOCK_CYCLIC
      g_a = NGA_Create (type, 2, dims, "A", NULL);
      if (!g_a)
        GA_Error ("create failed: A", n);
#else
      g_a = GA_Create_handle();
      GA_Set_data(g_a, 2, dims, type);
      GA_Set_array_name(g_a,"A");
#ifdef USE_SCALAPACK
      GA_Set_block_cyclic_proc_grid(g_a, block_size, proc_grid);
#else
      GA_Set_block_cyclic(g_a, block_size);
#endif
      GA_Allocate(g_a);
#endif
      /*find out the diagonal length of the matrix A */
      vdim = GA_MIN (dims[0], dims[1]);
      g_v = NGA_Create (type, 1, &vdim, "V", NULL);
      if (!g_v)
        GA_Error ("create failed:V", n);
      test_get_diagonal (g_a, g_v);
      GA_Destroy (g_a);
      GA_Destroy (g_v);
      break;
    case OP_NORM1:
#ifndef BLOCK_CYCLIC
      g_a = NGA_Create (type, 2, dims, "A", NULL);
      if (!g_a)
        GA_Error ("create failed: A", n);
#else
      g_a = GA_Create_handle();
      GA_Set_data(g_a, 2, dims, type);
      GA_Set_array_name(g_a,"A");
#ifdef USE_SCALAPACK
      GA_Set_block_cyclic_proc_grid(g_a, block_size, proc_grid);
#else
      GA_Set_block_cyclic(g_a, block_size);
#endif
      GA_Allocate(g_a);
#endif
      test_norm1 (g_a);
      GA_Destroy (g_a);
      break;

    case OP_NORM_INFINITY:
#ifndef BLOCK_CYCLIC
      g_a = NGA_Create (type, 2, dims, "A", NULL);
      if (!g_a)
        GA_Error ("create failed: A", n);
#else
      g_a = GA_Create_handle();
      GA_Set_data(g_a, 2, dims, type);
      GA_Set_array_name(g_a,"A");
#ifdef USE_SCALAPACK
      GA_Set_block_cyclic_proc_grid(g_a, block_size, proc_grid);
#else
      GA_Set_block_cyclic(g_a, block_size);
#endif
      GA_Allocate(g_a);
#endif
      test_norm_infinity (g_a);
      GA_Destroy (g_a);
      break;

    case OP_MEDIAN:
#ifndef BLOCK_CYCLIC
      g_a = NGA_Create (type, 2, dims, "A", NULL);
      if (!g_a)
        GA_Error ("create failed: A", n);
#else
      g_a = GA_Create_handle();
      GA_Set_data(g_a, 2, dims, type);
      GA_Set_array_name(g_a,"A");
#ifdef USE_SCALAPACK
      GA_Set_block_cyclic_proc_grid(g_a, block_size, proc_grid);
#else
      GA_Set_block_cyclic(g_a, block_size);
#endif
      GA_Allocate(g_a);
#endif
      /*duplicate g_a */
      g_b = GA_Duplicate (g_a, "B");
      if (!g_b)
        GA_Error ("duplicate failed: B", n);

      g_c = GA_Duplicate (g_a, "C");
      if (!g_c)
        GA_Error ("duplicate failed: C", n);
#if 0 /* test g_m is different from g_a, g_b, amd g_c */
      g_m = GA_Duplicate (g_a, "M");
      if (!g_m)
        GA_Error ("duplicate failed: M", n);
      test_median (g_a, g_b, g_c, g_m);
#else /* test g_m = g_c */
      test_median (g_a, g_b, g_c, g_a);
#endif
      GA_Destroy (g_a);
      GA_Destroy (g_b);
      GA_Destroy (g_c);
#if 0  /* test g_m is different from g_a, g_b, g_c */
      GA_Destroy (g_m);
#endif
      break;

    case OP_MEDIAN_PATCH:
#ifndef BLOCK_CYCLIC
      g_a = NGA_Create (type, 2, dims, "A", NULL);
      if (!g_a)
        GA_Error ("create failed: A", n);
#else
      g_a = GA_Create_handle();
      GA_Set_data(g_a, 2, dims, type);
      GA_Set_array_name(g_a,"A");
#ifdef USE_SCALAPACK
      GA_Set_block_cyclic_proc_grid(g_a, block_size, proc_grid);
#else
      GA_Set_block_cyclic(g_a, block_size);
#endif
      GA_Allocate(g_a);
#endif
      /*duplicate g_a */
      g_b = GA_Duplicate (g_a, "B");
      if (!g_b)
        GA_Error ("duplicate failed: B", n);

      g_c = GA_Duplicate (g_a, "C");
      if (!g_c)
        GA_Error ("duplicate failed: C", n);

      g_m = GA_Duplicate (g_a, "M");
      if (!g_m)
        GA_Error ("duplicate failed: M", n);

      test_median_patch (g_a, lo, hi, g_b, lo, hi, g_c, lo, hi, g_m, lo, hi);
      GA_Destroy (g_a);
      GA_Destroy (g_b);
      GA_Destroy (g_c);
      GA_Destroy (g_m);
      break;
    case OP_SCALE_ROWS:
#ifndef BLOCK_CYCLIC
      g_a = NGA_Create (type, 2, dims, "A", NULL);
      if (!g_a)
    GA_Error ("create failed: A", n);
#else
      g_a = GA_Create_handle();
      GA_Set_data(g_a, 2, dims, type);
      GA_Set_array_name(g_a,"A");
#ifdef USE_SCALAPACK
      GA_Set_block_cyclic_proc_grid(g_a, block_size, proc_grid);
#else
      GA_Set_block_cyclic(g_a, block_size);
#endif
      GA_Allocate(g_a);
#endif
      /*find out the diagonal length of the matrix A */
      vdim = dims[1];
      g_v = NGA_Create (type, 1, &vdim, "V", NULL);
      if (!g_v)
        GA_Error ("create failed:V", n);
      test_scale_rows (g_a, g_v);
      GA_Destroy (g_a);
      GA_Destroy (g_v);
      break;
    case OP_SCALE_COLS:
#ifndef BLOCK_CYCLIC
      g_a = NGA_Create (type, 2, dims, "A", NULL);
      if (!g_a)
        GA_Error ("create failed: A", n);
#else
      g_a = GA_Create_handle();
      GA_Set_data(g_a, 2, dims, type);
      GA_Set_array_name(g_a,"A");
#ifdef USE_SCALAPACK
      GA_Set_block_cyclic_proc_grid(g_a, block_size, proc_grid);
#else
      GA_Set_block_cyclic(g_a, block_size);
#endif
      GA_Allocate(g_a);
#endif
      /*find out the diagonal length of the matrix A */
      vdim = dims[0];
      g_v = NGA_Create (type, 1, &vdim, "V", NULL);
      if (!g_v)
        GA_Error ("create failed:V", n);
      test_scale_cols (g_a, g_v);
      GA_Destroy (g_a);
      GA_Destroy (g_v);
      break;

    default:
      GA_Error ("test_function: wrong operation.", op);
    }


}



int
main (argc, argv)
     int argc;
     char **argv;
{
  int heap = 20000, stack = 20000;
  int me, nproc;
  int type, op;

  type = C_DBL;

  MP_INIT(argc,argv);

  GA_INIT(argc,argv);        /* initialize GA */
  me = GA_Nodeid ();
  nproc = GA_Nnodes ();
  if (me == 0)
    {
      if (GA_Uses_fapi ())
    GA_Error ("Program runs with C array API only", 1);
      printf ("Using %ld processes\n", (long) nproc);
      fflush (stdout);
    }

  heap /= nproc;
  stack /= nproc;
  if (!MA_init (MT_F_DBL, stack, heap))
    GA_Error ("MA_init failed", stack + heap);    /* initialize memory allocator */


  for (op = 1; op < 12; op++)
    {
      if(me == 0) printf ("\n\n");
      if (me == 0)
    printf ("type = C_INT \t ");
      do_work (C_INT, op);
      if (me == 0)
    printf ("type = C_LONG \t ");
      do_work (C_LONG, op);
      if (me == 0)
    printf ("type = C_FLOAT \t ");
      do_work (C_FLOAT, op);
      if (me == 0)
    printf ("type = C_DBL \t ");
      do_work (C_DBL, op);
      if (me == 0)
    printf ("type = C_DCPL \t ");
      do_work (C_DCPL, op);
      if (me == 0)
    printf ("type = C_SCPL \t ");
      do_work (C_SCPL, op);
    }

  if(me==0) printf("\nSuccess\n");

  GA_Terminate ();
  
  MP_FINALIZE();

  return 0;
}
