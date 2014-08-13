#if HAVE_CONFIG_H
#   include "config.h"
#endif

/**\
File: elempatch.c
Purpose: test the interfaces:

GA_Abs_value_patch(g_a)
GA_Add_constant_patch(g_a, alpha)
GA_Recip_patch_patch(g_a)
GA_Elem_multiply_patch(g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi)
GA_Elem_divide_patch(g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi)
GA_Elem_maximum_patch(g_a, alo, ahi, g_b, blo, bhi, g_c, clo, ch
GA_Elem_minimum_patch(g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi)

that are for TAO/Global Array Project

Author:

Limin Zhang, Ph.D.
Mathematics Department
Columbia Basin College
Pasco, WA 99301

Mentor:

Jarek Nieplocha
Pacific Northwest National Laboratory

Date: Jauary 30, 2002
Revised on February 26, 2002.

\**/

#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STDIO_H
#   include <stdio.h>
#endif

#include "ga.h"
#include "macdecls.h"
#include "mp3.h"
#include "../src/globalp.h"
#include "ga-papi.h"

#define BLOCK_CYCLIC 0

#if BLOCK_CYCLIC
#define USE_SCALAPACK 1
#endif

#ifndef GA_HALF_MAX_INT 
#define GA_HALF_MAX_INT ((((int)1) << ((int)(8*sizeof(int))-2)) - 1)
#endif

#ifndef GA_INFINITY_I
#define GA_INFINITY_I (GA_HALF_MAX_INT + GA_HALF_MAX_INT + 1)
/* 
   Original value below.
   Seemed too small arbitrarily.
#define GA_INFINITY_I 100000
 */
#endif

#ifndef GA_NEGATIVE_INFINITY_I
#define GA_NEGATIVE_INFINITY_I (- GA_INFINITY_I)


/* 
   Original value below. 
   Seemed too small arbitrarily.
#define GA_NEGATIVE_INFINITY_I -100000
 */
#endif

#ifndef GA_HALF_MAX_LONG
#define GA_HALF_MAX_LONG ((((long)1) << ((int)(8*sizeof(long))-2)) - 1)
#endif

#ifndef GA_INFINITY_L
#define GA_INFINITY_L (GA_HALF_MAX_LONG + GA_HALF_MAX_LONG + 1)
/* Original value was
#define GA_INFINITY_L 100000
 */
#endif

#ifndef GA_NEGATIVE_INFINITY_L
#define GA_NEGATIVE_INFINITY_L (- GA_INFINITY_L)
#endif
/* 
   Original value was:
#define GA_NEGATIVE_INFINITY_L -100000
 */

/* 
   Modified by Doug Baxter 01/24/04 to distinguish between
   Double inifinity and float infinity.
#ifndef GA_INFINITY
#define GA_INFINITY 1.0e20
#endif

#ifndef GA_NEGATIVE_INFINITY
#define GA_NEGATIVE_INFINITY -1.0e20
#endif
 */
#ifndef GA_INFINITY_F
#define GA_INFINITY_F 1.0e37
#endif
/*
Original value below.
#define GA_INFINITY_F 1.0e20
 */
#ifndef GA_NEGATIVE_INFINITY_F
#define GA_NEGATIVE_INFINITY_F -1.0e37
#endif
/*
Original value below.
#define GA_NEGATIVE_INFINITY_F -1.0e20
 */
#ifndef GA_INFINITY_D
#define GA_INFINITY_D 1.0e307
#endif
/*
Original value below.
#define GA_INFINITY_D 1.0e20
 */
#ifndef GA_NEGATIVE_INFINITY_D
#define GA_NEGATIVE_INFINITY_D -1.0e307
#endif


#define THRESH 1e-5
#define MISMATCHED(x,y) GA_ABS((x)-(y))>=THRESH

#define N 100
#define BLOCK_SIZE 20
#define OP_ELEM_MULT 0
#define OP_ELEM_DIV 1
#define OP_ELEM_MAX 2
#define OP_ELEM_MIN 3
#define OP_ABS 4
#define OP_ADD_CONST 5
#define OP_RECIP 6
#define OP_STEP_MAX 7
#define OP_STEP_BOUND_INFO 8
#define MY_TYPE 2002

Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM], _ga_work[MAXDIM];
#  define COPYINDEX_C2F(carr, farr, n){\
  int i; for(i=0; i< (n); i++)(farr)[n-i-1]=(Integer)(carr)[i]+1;}

void nga_vfill_patch(Integer *g_a, Integer *lo, Integer *hi);
void nga_pnfill_patch(Integer *g_a, Integer *lo, Integer *hi);

void NGA_Vfill_patch(int g_a, int lo[], int hi[])
{
  Integer a=(Integer)g_a;
  Integer ndim = pnga_ndim(a);
  COPYINDEX_C2F(lo,_ga_lo, ndim);
  COPYINDEX_C2F(hi,_ga_hi, ndim);

  nga_vfill_patch(&a, _ga_lo, _ga_hi);
}


void NGA_Pnfill_patch(int g_a, int lo[], int hi[])
{
  Integer a=(Integer)g_a;
  Integer ndim = pnga_ndim(a);
  COPYINDEX_C2F(lo,_ga_lo, ndim);
  COPYINDEX_C2F(hi,_ga_hi, ndim);

  nga_pnfill_patch(&a, _ga_lo, _ga_hi);
}

int
ifun (int k)
{
  int result;
  result = -k - 1;
  result = -2;
  return result;
}

int
ifun2 (int k)
{
  int result;
  result = k + 1;
  result = -3;
  return result;
}

void
fill_func (int nelem, int type, void *buf)
{
  int i;


  switch (type)
  {
    case C_FLOAT:
      for (i = 0; i < nelem; i++)
        ((float *) buf)[i] = (float) ifun (i);
      break;
    case C_LONG:
      for (i = 0; i < nelem; i++)
        ((long *) buf)[i] = (long) ifun (i);
      break;
    case C_DBL:
      for (i = 0; i < nelem; i++)
        ((double *) buf)[i] = (double) ifun (i);
      break;
    case C_DCPL:
      for (i = 0; i < 2 * nelem; i++)
        ((double *) buf)[i] = (double) ifun (i);
      break;
    case C_SCPL:
      for (i = 0; i < 2 * nelem; i++)
        ((float *) buf)[i] = (float) ifun (i);
      break;
    case C_INT:
      for (i = 0; i < nelem; i++)
        ((int *) buf)[i] = ifun (i);
      break;
    default:
      GA_Error (" wrong data type ", type);

  }
}

void
fill_func2 (int nelem, int type, void *buf)
{
  /* int i,size=MA_sizeof(MT_CHAR,type,1);*/

  int i;

  switch (type)
  {
    case C_FLOAT:
      for (i = 0; i < nelem; i++)
        ((float *) buf)[i] = (float) ifun2 (i);
      break;
    case C_LONG:
      for (i = 0; i < nelem; i++)
        ((long *) buf)[i] = (long) ifun2 (i);
      break;
    case C_DBL:
      for (i = 0; i < nelem; i++)
        ((double *) buf)[i] = (double) ifun2 (i);
      break;
    case C_DCPL:
      for (i = 0; i < 2 * nelem; i++)
        ((double *) buf)[i] = (double) ifun2 (i);
      break;
    case C_SCPL:
      for (i = 0; i < 2 * nelem; i++)
        ((float *) buf)[i] = (float) ifun2 (i);
      break;
    case C_INT:
      for (i = 0; i < nelem; i++)
        ((int *) buf)[i] = ifun2 (i);
      break;
    default:
      GA_Error (" wrong data type ", type);

  }
}

void
fill_func3 (int nelem, int type, void *buf)
  /*taking the absolute of the ifun() */
{
  /*int i,size=MA_sizeof(MT_CHAR,type,1);*/

  int i;

  switch (type)
  {
    case C_FLOAT:
      for (i = 0; i < nelem; i++)
        ((float *) buf)[i] = (float) GA_ABS (ifun (i));
      break;
    case C_LONG:
      for (i = 0; i < nelem; i++)
        ((long *) buf)[i] = (long) GA_ABS (ifun (i));
      break;
    case C_DBL:
      for (i = 0; i < nelem; i++)
        ((double *) buf)[i] = (double) GA_ABS (ifun (i));
      break;
    case C_DCPL:
      for (i = 0; i < 2 * nelem - 1; i = i + 2)
      {
        ((double *) buf)[i] =
          sqrt ((double)
              (ifun (i) * ifun (i) + ifun (i + 1) * ifun (i + 1)));
        ((double *) buf)[i + 1] = 0.0;
      }
      break;
    case C_SCPL:
      for (i = 0; i < 2 * nelem - 1; i = i + 2)
      {
        ((float *) buf)[i] =
          sqrt ((float)
              (ifun (i) * ifun (i) + ifun (i + 1) * ifun (i + 1)));
        ((float *) buf)[i + 1] = 0.0;
      }
      break;
    case C_INT:
      for (i = 0; i < nelem; i++)
        ((int *) buf)[i] = GA_ABS (ifun (i));
      break;
    default:
      GA_Error (" wrong data type ", type);

  }
}






int
test_fun (int type, int dim, int OP)
{
  DoubleComplex *dcptr;
  int ld[7],locp[7],hicp[7];
  void *boundminx=NULL,*boundmaxx=NULL,*wolfeminx=NULL;
  double boundmind=0,boundmaxd=0,wolfemind=0,aboundmind=0,aboundmaxd=0,awolfemind=0;
  long boundminl=0,boundmaxl=0,wolfeminl=0,aboundminl=0,aboundmaxl=0,awolfeminl=0;
  float boundminf=0,boundmaxf=0,wolfeminf=0,aboundminf=0,aboundmaxf=0,awolfeminf=0;
  int boundmini=0,boundmaxi=0,wolfemini=0,aboundmini=0,aboundmaxi=0,awolfemini=0;

  int g_a, g_b, g_c, g_d, g_e;
  int g_f, g_g, g_h, g_i, g_j;
  int g_k, g_l, g_m, g_n;
  int n = N;
  int me = GA_Nodeid ();
  int i;
  int dims[MAXDIM];
  int lo[MAXDIM], hi[MAXDIM];
  int index[MAXDIM];
  int index2[MAXDIM];
  int index3[MAXDIM];
  int block_size[MAXDIM], proc_grid[MAXDIM], proc_cnt;
  int needs_scaled_result;
  void *val=NULL;
  void *val3=NULL;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;
  SingleComplex fcval;
  void *val2=NULL;
  void *val4=NULL;
  int ival2 = -3;
  double dval2 = -3.0;
  float fval2 = -3.0;
  long lval2 = -3;
  void *val5=NULL;
  int ival5 = 5;
  long lval5 = 5;
  double dval5 = 5.0;
  float fval5 = 5.0;
  DoubleComplex dcval2;
  DoubleComplex dcval3;
  DoubleComplex dcval4;
  DoubleComplex dcval5;
  DoubleComplex dcval6;
  DoubleComplex dcval7;
  SingleComplex fcval2;
  SingleComplex fcval3;
  SingleComplex fcval4;
  SingleComplex fcval5;
  SingleComplex fcval6;
  SingleComplex fcval7;
  void *vresult=NULL;
  int ivresult;
  double dvresult;
  float fvresult;
  long lvresult;
  DoubleComplex dcvresult;
  SingleComplex fcvresult;
  void *bvresult=NULL;
  DoubleComplex dcbvresult;
  SingleComplex fcbvresult;
  int ok = 1;
  int result=0,result2=0,result3=0;
  void *max=NULL;
  int imax;
  float fmax;
  long lmax;
  double dmax;
  DoubleComplex dcmax;
  SingleComplex fcmax;
  void *max2=NULL;
  DoubleComplex dcmax2;
  SingleComplex fcmax2;
  void *max3=NULL;
  DoubleComplex dcmax3;
  SingleComplex fcmax3;

  void *alpha=NULL, *beta=NULL;
  int ai = 1, bi = -1;
  long al = 1, bl = -1;
  float af = 1.0, bf = -1.0;
  double ad = 1.0, bd = -1.0;
  DoubleComplex adc, bdc;
  SingleComplex afc, bfc;
  double x1, x2, x3, x4;
  float fx1, fx2, fx3, fx4;
  void    *resultx=NULL;
  long    resultl=0,aresultl=0;
  double  resultd=0,aresultd=0;
  float   resultf=0,aresultf=0;
  int resulti=0,aresulti=0;

  adc.real = 1.0;
  adc.imag = 0.0;
  bdc.real = -1.0;
  bdc.imag = 0.0;

  afc.real = 1.0;
  afc.imag = 0.0;
  bfc.real = -1.0;
  bfc.imag = 0.0;

  needs_scaled_result = 0;

  dcval.real = -sin (3.0);
  dcval.imag = -cos (3.0);
  dcval2.real = 2 * sin (3.0);
  dcval2.imag = 2 * cos (3.0);
  dcval3.real = dcval.real*1.0e2;
  dcval3.imag = dcval.imag*1.0e2;
  dcval4.real = dcval2.real*1.0e2;
  dcval4.imag = dcval2.imag*1.0e2;
  dcval5.real = 5.0;
  dcval5.imag = 0.0;
  dcval6.real = dcval3.imag;
  dcval6.imag = dcval3.real;
  dcval7.real = dcval4.imag;
  dcval7.imag = dcval4.real;

  fcval.real = -sin (3.0);
  fcval.imag = -cos (3.0);
  fcval2.real = 2 * sin (3.0);
  fcval2.imag = 2 * cos (3.0);
  fcval3.real = fcval.real*1.0e2;
  fcval3.imag = fcval.imag*1.0e2;
  fcval4.real = fcval2.real*1.0e2;
  fcval4.imag = fcval2.imag*1.0e2;
  fcval5.real = 5.0;
  fcval5.imag = 0.0;
  fcval6.real = fcval3.imag;
  fcval6.imag = fcval3.real;
  fcval7.real = fcval4.imag;
  fcval7.imag = fcval4.real;


  proc_cnt=1;
  for (i = 0; i < dim; i++) {
    dims[i] = N;
    block_size[i] = BLOCK_SIZE;
    if (i<dim-1 && proc_cnt < GA_Nnodes()) {
      proc_grid[i] = 2;
      proc_cnt *= 2;
    } else if (proc_cnt >= GA_Nnodes()) {
      proc_grid[i] = 1;
    } else {
      proc_grid[i] = GA_Nnodes()/proc_cnt;
    }
  }

  for (i = 0; i < dim; i++)
  {
    lo[i] = 0;
    hi[i] = N - 1;
  }
#if BLOCK_CYCLIC
  g_a = GA_Create_handle();
  GA_Set_data(g_a,dim,dims,type);
  GA_Set_array_name(g_a,"A");
# if USE_SCALAPACK
  GA_Set_block_cyclic_proc_grid(g_a,block_size,proc_grid);
# else
  GA_Set_block_cyclic(g_a,block_size);
# endif
  GA_Allocate(g_a);
#else
  g_a = NGA_Create (type, dim, dims, "A", NULL);
  if (!g_a)
    GA_Error ("create failed: A", n);
#endif

  g_b = GA_Duplicate (g_a, "B");
  if (!g_b)
    GA_Error ("duplicate failed: B", n);

  g_c = GA_Duplicate (g_a, "C");
  if (!g_c)
    GA_Error ("duplicate failed: C", n);

  g_d = GA_Duplicate (g_a, "D");
  if (!g_d)
    GA_Error ("duplicate failed: D", n);

  g_e = GA_Duplicate (g_a, "E");
  if (!g_e)
    GA_Error ("duplicate failed: E", n);

  g_f = GA_Duplicate (g_a, "F");
  if (!g_f)
    GA_Error ("duplicate failed: F", n);

  g_g = GA_Duplicate (g_a, "G");
  if (!g_g)
    GA_Error ("duplicate failed: G", n);

  g_h = GA_Duplicate (g_a, "H");
  if (!g_h)
    GA_Error ("duplicate failed: H", n);

  g_i = GA_Duplicate (g_a, "I");
  if (!g_i)
    GA_Error ("duplicate failed: I", n);

  g_j = GA_Duplicate (g_a, "J");
  if (!g_j)
    GA_Error ("duplicate failed: J", n);

  g_k = GA_Duplicate (g_a, "K");
  if (!g_k)
    GA_Error ("duplicate failed: K", n);

  g_l = GA_Duplicate (g_a, "L");
  if (!g_l)
    GA_Error ("duplicate failed: L", n);

  g_m = GA_Duplicate (g_a, "M");
  if (!g_m)
    GA_Error ("duplicate failed: M", n);

  g_n = GA_Duplicate (g_a, "N");
  if (!g_m)
    GA_Error ("duplicate failed: N", n);
  /*initialize  with zero */
  GA_Zero (g_a);
  GA_Zero (g_b);
  GA_Zero (g_c);
  GA_Zero (g_d);
  GA_Zero (g_e);
  GA_Zero (g_f);
  GA_Zero (g_g);
  GA_Zero (g_h);
  GA_Zero (g_i);
  GA_Zero (g_j);
  GA_Zero (g_k);
  GA_Zero (g_l);
  GA_Zero (g_m);
  GA_Zero (g_n);

  switch (type)
  {
    case C_INT:
      val = &ival;
      val2 = &ival2;
      val5 = &ival5;
      vresult = &ivresult;
      resultx = &resulti;
      boundminx = &boundmini;
      boundmaxx = &boundmaxi;
      wolfeminx = &wolfemini;
      break;
    case C_DCPL:
      val = &dcval;
      val2 = &dcval2;
      val3 = &dcval3;
      val4 = &dcval4;
      val5 = &dcval5;
      vresult = &dcvresult;
      bvresult = &dcbvresult;
      break;
    case C_SCPL:
      val = &fcval;
      val2 = &fcval2;
      val3 = &fcval3;
      val4 = &fcval4;
      val5 = &fcval5;
      vresult = &fcvresult;
      bvresult = &fcbvresult;
      break;

    case C_DBL:
      val = &dval;
      val2 = &dval2;
      val5 = &dval5;
      vresult = &dvresult;
      resultx = &resultd;
      boundminx = &boundmind;
      boundmaxx = &boundmaxd;
      wolfeminx = &wolfemind;
      break;
    case C_FLOAT:
      val = &fval;
      val2 = &fval2;
      val5 = &fval5;
      vresult = &fvresult;
      resultx = &resultf;
      boundminx = &boundminf;
      boundmaxx = &boundmaxf;
      wolfeminx = &wolfeminf;
      break;
    case C_LONG:
      val = &lval;
      val2 = &lval2;
      val5 = &lval5;
      vresult = &lvresult;
      resultx = &resultl;
      boundminx = &boundminl;
      boundmaxx = &boundmaxl;
      wolfeminx = &wolfeminl;
      break;
    default:
      GA_Error ("wrong data type.", type);
  }


  NGA_Fill_patch (g_a, lo, hi, val);
  NGA_Fill_patch (g_b, lo, hi, val2);
  NGA_Pnfill_patch (g_j, lo, hi);
  switch (OP)
  {
    double tmp, tmp2;
    float  tmpf, tmp2f;
    case OP_ABS:
    if (me == 0)
      printf ("Testing GA_Abs_value...");
    GA_Abs_value_patch (g_a, lo, hi);
    ivresult = GA_ABS (ival);
    dvresult = GA_ABS (dval);
    fvresult = GA_ABS (fval);
    lvresult = GA_ABS (lval);

#if 0
    if (GA_ABS(dcval.real) >= GA_ABS(dcval.imag)) {
      if (dcval.real == (double)0.0) {
        dcvresult.real = (double)0.0;
      } else {
        x1 = dcval.imag/dcval.real;
        dcvresult.real = GA_ABS(dcval.real)*sqrt(((double)1.0)+(x1*x1));
      }
    } else {
      x1 = dcval.real/dcval.imag;
      dcvresult.real = GA_ABS(dcval.imag)*sqrt(((double)1.0)+(x1*x1));
    }
#endif
    dcvresult.real = sqrt(dcval.real*dcval.real+dcval.imag*dcval.imag);
    dcvresult.imag = 0.0;
    NGA_Fill_patch (g_d, lo, hi, vresult);

    if (type == C_DCPL) {
      needs_scaled_result = 1;
      NGA_Fill_patch(g_f,lo,hi,val3);
      GA_Abs_value_patch (g_f, lo, hi);
#if 0
      if (GA_ABS(dcval3.real) >= GA_ABS(dcval3.imag)) {
        if (dcval3.real == (double)0.0) {
          dcbvresult.real = (double)0.0;
        } else {
          x1 = dcval3.imag/dcval3.real;
          dcbvresult.real = GA_ABS(dcval3.real)*sqrt(((double)1.0)+(x1*x1));
        }
      } else {
        x1 = dcval3.real/dcval3.imag;
        dcbvresult.real = GA_ABS(dcval3.imag)*sqrt(((double)1.0)+(x1*x1));
      }
#endif
      dcbvresult.real = sqrt(dcval3.real*dcval3.real+dcval3.imag*dcval3.imag);
      dcbvresult.imag = (double)0.0;
      NGA_Fill_patch (g_i, lo, hi, bvresult);

      NGA_Fill_patch(g_k,lo,hi,&dcval6);
      GA_Abs_value_patch (g_k, lo, hi);
#if 0
      if (GA_ABS(dcval6.real) >= GA_ABS(dcval6.imag)) {
        if (dcval6.real == (double)0.0) {
          dcbvresult.real = (double)0.0;
        } else {
          x1 = dcval6.imag/dcval6.real;
          dcbvresult.real = GA_ABS(dcval6.real)*sqrt(((double)1.0)+(x1*x1));
        }
      } else {
        x1 = dcval6.real/dcval6.imag;
        dcbvresult.real = GA_ABS(dcval6.imag)*sqrt(((double)1.0)+(x1*x1));
      }
#endif
      dcbvresult.real = sqrt(dcval6.real*dcval6.real+dcval6.imag*dcval6.imag);
      dcbvresult.imag = (double)0.0;
      NGA_Fill_patch (g_n, lo, hi, bvresult);
    }
    if (type == C_SCPL) {
      needs_scaled_result = 1;
      NGA_Fill_patch(g_f,lo,hi,val3);
      GA_Abs_value_patch (g_f, lo, hi);
#if 0
      if (GA_ABS(fcval3.real) >= GA_ABS(fcval3.imag)) {
        if (fcval3.real == (float )0.0) {
          fcbvresult.real = (float )0.0;
        } else {
          fx1 = fcval3.imag/fcval3.real;
          fcbvresult.real = GA_ABS(fcval3.real)*sqrt(((float)1.0)+(fx1*fx1));
        }
      } else {
        fx1 = fcval3.real/fcval3.imag;
        fcbvresult.real = GA_ABS(fcval3.imag)*sqrt(((float)1.0)+(fx1*fx1));
      }
#endif
      fcbvresult.real = sqrt(fcval3.real*fcval3.real+fcval3.imag*fcval3.imag);
      fcbvresult.imag = (double)0.0;
      NGA_Fill_patch (g_i, lo, hi, bvresult);
      NGA_Fill_patch(g_k,lo,hi,&dcval6);
      GA_Abs_value_patch (g_k, lo, hi);
#if 0
      if (GA_ABS(fcval6.real) >= GA_ABS(fcval6.imag)) {
        if (fcval6.real == (float)0.0) {
          fcbvresult.real = (float)0.0;
        } else {
          fx1 = fcval6.imag/fcval6.real;
          fcbvresult.real = GA_ABS(fcval6.real)*sqrt(((float)1.0)+(fx1*fx1));
        }
      } else {
        fx1 = fcval6.real/fcval6.imag;
        fcbvresult.real = GA_ABS(fcval6.imag)*sqrt(((float)1.0)+(fx1*fx1));
      }
#endif
      fcbvresult.real = sqrt(fcval6.real*fcval6.real+fcval6.imag*fcval6.imag);
      fcbvresult.imag = (double)0.0;
      NGA_Fill_patch (g_n, lo, hi, bvresult);
    }
    break;
    case OP_ADD_CONST:
    if (me == 0)
      printf ("Testing GA_Add_const...");
    GA_Add_constant_patch (g_a, lo, hi, val2);
    ivresult = ival + ival2;
    dvresult = dval + dval2;
    fvresult = fval + fval2;
    lvresult = lval + lval2;
    dcvresult.real = dcval.real + dcval2.real;
    dcvresult.imag = dcval.imag + dcval2.imag;
    fcvresult.real = fcval.real + fcval2.real;
    fcvresult.imag = fcval.imag + fcval2.imag;
    NGA_Fill_patch (g_d, lo, hi, vresult);
    break;
    case OP_RECIP:
    if (me == 0)
      printf ("Testing GA_Recip...");
    GA_Recip_patch (g_a, lo, hi);
    ivresult = ((int)1) / ival;
    dvresult = ((double)1.0) / dval;
    fvresult = ((float)1.0) / fval;
    lvresult = ((long)1) / lval;

    if (GA_ABS(dcval.real) >= GA_ABS(dcval.imag)) {
      if (dcval.real != (double)0.0) {
        tmp = dcval.imag/dcval.real;
        tmp2 = ((double)1.0)/((((double)1.0)+(tmp*tmp))*dcval.real);
        dcvresult.real = tmp2;
        dcvresult.imag = -tmp * tmp2;
      } else {
        printf("Error in testing GA_Recip dcval = 0.0\n");
      }
    } else {
      tmp = dcval.real/dcval.imag;
      tmp2 = ((double)1.0)/((((double)1.0)+(tmp*tmp))*dcval.imag);
      dcvresult.real = tmp * tmp2;
      dcvresult.imag = -tmp2;
    }
    NGA_Fill_patch (g_d, lo, hi, vresult);

    if (type == C_DCPL) {
      needs_scaled_result = 1;
      NGA_Fill_patch (g_f, lo, hi, val3);
      GA_Recip_patch (g_f, lo, hi);
      if (GA_ABS(dcval3.real) >= GA_ABS(dcval3.imag)) {
        if (dcval3.real == (double)0.0) {
          printf("Error testing GA_Recip, dcval3.real = 0.0\n");
        } else {
          tmp = dcval3.imag/dcval3.real;
          tmp2 = ((double)1.0)/((((double)1.0)+(tmp*tmp))*dcval3.real);
          dcbvresult.real = tmp2;
          dcbvresult.imag = -tmp * tmp2;
        }
      } else {
        tmp = dcval3.real/dcval3.imag;
        tmp2 = ((double)1.0)/((((double)1.0)+(tmp*tmp))*dcval3.imag);
        dcbvresult.real = tmp * tmp2;
        dcbvresult.imag = -tmp2;
      }
      NGA_Fill_patch (g_i, lo, hi, bvresult);
      NGA_Fill_patch(g_k,lo,hi,&dcval6);
      GA_Recip_patch (g_k, lo, hi);
      if (GA_ABS(dcval6.real) >= GA_ABS(dcval6.imag)) {
        if (dcval6.real == (double)0.0) {
          printf("Error testing GA_Recip, dcval6.real = 0.0\n");
        } else {
          tmp = dcval6.imag/dcval6.real;
          tmp2 = ((double)1.0)/((((double)1.0)+(tmp*tmp))*dcval6.real);
          dcbvresult.real = tmp2;
          dcbvresult.imag = -tmp * tmp2;
        }
      } else {
        tmp = dcval6.real/dcval6.imag;
        tmp2 = ((double)1.0)/((((double)1.0)+(tmp*tmp))*dcval6.imag);
        dcbvresult.real = tmp * tmp2;
        dcbvresult.imag = -tmp2;
      }
      NGA_Fill_patch (g_n, lo, hi, bvresult);
    }
    if (type == C_SCPL) {
      needs_scaled_result = 1;
      NGA_Fill_patch (g_f, lo, hi, val3);
      GA_Recip_patch (g_f, lo, hi);
      if (GA_ABS(fcval3.real) >= GA_ABS(fcval3.imag)) {
        if (fcval3.real == (float )0.0) {
          printf("Error testing GA_Recip, fcval3.real = 0.0\n");
        } else {
          tmpf = fcval3.imag/fcval3.real;
          tmp2f = ((float)1.0)/((((float)1.0)+(tmpf*tmpf))*fcval3.real);
          fcbvresult.real = tmp2f;
          fcbvresult.imag = -tmpf * tmp2f;
        }
      } else {
        tmpf = fcval3.real/fcval3.imag;
        tmp2f = ((float)1.0)/((((float)1.0)+(tmpf*tmpf))*fcval3.imag);
        fcbvresult.real = tmpf * tmp2f;
        fcbvresult.imag = -tmp2f;
      }
      NGA_Fill_patch (g_i, lo, hi, bvresult);
      NGA_Fill_patch(g_k,lo,hi,&dcval6);
      GA_Recip_patch (g_k, lo, hi);
      if (GA_ABS(fcval6.real) >= GA_ABS(fcval6.imag)) {
        if (fcval6.real == (float)0.0) {
          printf("Error testing GA_Recip, fcval6.real = 0.0\n");
        } else {
          tmpf = fcval6.imag/fcval6.real;
          tmp2f = ((float)1.0)/((((float)1.0)+(tmpf*tmpf))*fcval6.real);
          fcbvresult.real = tmp2f;
          fcbvresult.imag = -tmpf * tmp2f;
        }
      } else {
        tmpf = fcval6.real/fcval6.imag;
        tmp2f = ((float)1.0)/((((float)1.0)+(tmpf*tmpf))*fcval6.imag);
        fcbvresult.real = tmpf * tmp2f;
        fcbvresult.imag = -tmp2f;
      }
      NGA_Fill_patch (g_n, lo, hi, bvresult);
    }
    break;
    case OP_ELEM_MULT:
    if (me == 0)
      printf ("Testing GA_Elem_multiply...");
    NGA_Fill_patch (g_b, lo, hi, val2);
    /* g_c is different from g_a or g_b*/
    GA_Elem_multiply_patch (g_a, lo, hi, g_b, lo, hi, g_c, lo, hi);

    ivresult = ival * ival2;
    dvresult = dval * dval2;
    fvresult = fval * fval2;
    lvresult = lval * lval2;
    dcvresult.real = dcval.real * dcval2.real - dcval.imag * dcval2.imag;
    dcvresult.imag = dcval.real * dcval2.imag + dcval2.real * dcval.imag;
    NGA_Fill_patch (g_d, lo, hi, vresult);
    break;
    case OP_ELEM_DIV:
    if (me == 0)
      printf ("Testing GA_Elem_divide...");
    NGA_Fill_patch (g_b, lo, hi, val2);
    GA_Elem_divide_patch (g_a, lo, hi, g_b, lo, hi, g_c, lo, hi);
    ivresult = ival / ival2;
    dvresult = dval / dval2;
    fvresult = fval / fval2;
    lvresult = lval / lval2;
    dcvresult.real = 0.0;
    dcvresult.imag = 0.0;

    if (GA_ABS(dcval2.real) >= GA_ABS(dcval2.imag)) {
      if (dcval2.real != (double)0.0) {
        tmp = dcval2.imag/dcval2.real;
        tmp2 = ((double)1.0)/(dcval2.real*(((double)1.0)+(tmp*tmp)));
        dcvresult.real = (dcval.real + dcval.imag*tmp)*tmp2;
        dcvresult.imag = (dcval.imag - dcval.real*tmp)*tmp2;
      } else {
        printf("Error in testing GA_Elem_divide dcval = 0.0\n");
      } 
    } else {
      tmp = dcval2.real/dcval2.imag;
      tmp2 = 1.0/(dcval2.imag*(1.0+(tmp*tmp)));
      dcvresult.real = (dcval.real*tmp + dcval.imag)*tmp2;
      dcvresult.imag = (dcval.imag*tmp - dcval.real)*tmp2;
    }
    NGA_Fill_patch (g_d, lo, hi, vresult);

    if (type == C_DCPL) {
      needs_scaled_result = 1;
      NGA_Fill_patch (g_f, lo, hi, val3);
      NGA_Fill_patch (g_g, lo, hi, val4);
      GA_Elem_divide_patch (g_f, lo, hi, g_g, lo, hi, g_h, lo, hi);
      dcbvresult.real = (double)0.0;
      dcbvresult.imag = (double)0.0;
      if (GA_ABS(dcval4.real) >= GA_ABS(dcval4.imag)) {
        if (dcval4.real != (double)0.0) {
          tmp = dcval4.imag/dcval4.real;
          tmp2 = ((double)1.0)/(dcval4.real*(((double)1.0)+(tmp*tmp)));
          dcbvresult.real = (dcval3.real + dcval3.imag*tmp)*tmp2;
          dcbvresult.imag = (dcval3.imag - dcval3.real*tmp)*tmp2;
        } else {
          printf("Error in testing GA_Elem_divide dcval4 = 0.0\n");
        } 
      } else {
        tmp = dcval4.real/dcval4.imag;
        tmp2 = ((double)1.0)/(dcval4.imag*(((double)1.0)+(tmp*tmp)));
        dcbvresult.real = (dcval3.real*tmp + dcval3.imag)*tmp2;
        dcbvresult.imag = (dcval3.imag*tmp - dcval3.real)*tmp2;
      }
      NGA_Fill_patch (g_i, lo, hi, bvresult);
      NGA_Fill_patch (g_k, lo, hi, &dcval6);
      NGA_Fill_patch (g_l, lo, hi, &dcval7);
      GA_Elem_divide_patch (g_k, lo, hi, g_l, lo, hi, g_m, lo, hi);
      dcbvresult.real = (double)0.0;
      dcbvresult.imag = (double)0.0;
      if (GA_ABS(dcval7.real) >= GA_ABS(dcval7.imag)) {
        if (dcval7.real != (double)0.0) {
          tmp = dcval7.imag/dcval7.real;
          tmp2 = ((double)1.0)/(dcval7.real*(((double)1.0)+(tmp*tmp)));
          dcbvresult.real = (dcval6.real + dcval6.imag*tmp)*tmp2;
          dcbvresult.imag = (dcval6.imag - dcval6.real*tmp)*tmp2;
        } else {
          printf("Error in testing GA_Elem_divide dcval7 = 0.0\n");
        } 
      } else {
        tmp = dcval7.real/dcval7.imag;
        tmp2 = ((double)1.0)/(dcval7.imag*(((double)1.0)+(tmp*tmp)));
        dcbvresult.real = (dcval6.real*tmp + dcval6.imag)*tmp2;
        dcbvresult.imag = (dcval6.imag*tmp - dcval6.real)*tmp2;
      }
      NGA_Fill_patch (g_n, lo, hi, bvresult);
    }
    if (type == C_SCPL) {
      needs_scaled_result = 1;
      NGA_Fill_patch (g_f, lo, hi, val3);
      NGA_Fill_patch (g_g, lo, hi, val4);
      GA_Elem_divide_patch (g_f, lo, hi, g_g, lo, hi, g_h, lo, hi);
      fcbvresult.real = (float)0.0;
      fcbvresult.imag = (float)0.0;
      if (GA_ABS(fcval4.real) >= GA_ABS(fcval4.imag)) {
        if (fcval4.real != (float)0.0) {
          tmpf = fcval4.imag/fcval4.real;
          tmp2f = ((float)1.0)/(fcval4.real*(((float)1.0)+(tmpf*tmpf)));
          fcbvresult.real = (fcval3.real + fcval3.imag*tmpf)*tmp2f;
          fcbvresult.imag = (fcval3.imag - fcval3.real*tmpf)*tmp2f;
        } else {
          printf("Error in testing GA_Elem_divide fcval4 = 0.0\n");
        } 
      } else {
        tmpf = fcval4.real/fcval4.imag;
        tmp2f = ((float)1.0)/(fcval4.imag*(((float)1.0)+(tmpf*tmpf)));
        fcbvresult.real = (fcval3.real*tmpf + fcval3.imag)*tmp2f;
        fcbvresult.imag = (fcval3.imag*tmpf - fcval3.real)*tmp2f;
      }
      NGA_Fill_patch (g_i, lo, hi, bvresult);
      NGA_Fill_patch (g_k, lo, hi, &dcval6);
      NGA_Fill_patch (g_l, lo, hi, &dcval7);
      GA_Elem_divide_patch (g_k, lo, hi, g_l, lo, hi, g_m, lo, hi);
      fcbvresult.real = (float)0.0;
      fcbvresult.imag = (float)0.0;
      if (GA_ABS(fcval7.real) >= GA_ABS(fcval7.imag)) {
        if (fcval7.real != (float)0.0) {
          tmpf = fcval7.imag/fcval7.real;
          tmp2f = ((float)1.0)/(fcval7.real*(((float)1.0)+(tmpf*tmpf)));
          fcbvresult.real = (fcval6.real + fcval6.imag*tmpf)*tmp2f;
          fcbvresult.imag = (fcval6.imag - fcval6.real*tmpf)*tmp2f;
        } else {
          printf("Error in testing GA_Elem_divide fcval7 = 0.0\n");
        } 
      } else {
        tmpf = fcval7.real/fcval7.imag;
        tmp2f = ((float)1.0)/(fcval7.imag*(((float)1.0)+(tmpf*tmpf)));
        fcbvresult.real = (fcval6.real*tmpf + fcval6.imag)*tmp2f;
        fcbvresult.imag = (fcval6.imag*tmpf - fcval6.real)*tmp2f;
      }
      NGA_Fill_patch (g_n, lo, hi, bvresult);
    }
    break;

    case OP_ELEM_MAX:
    if (me == 0)
      printf ("Testing GA_Elem_maximum...");
    /*NGA_Fill_patch (g_b, lo, hi, val2);*/
    GA_Elem_maximum_patch (g_a, lo, hi, g_b, lo, hi, g_c, lo, hi);
    ivresult = GA_MAX (ival, ival2);
    dvresult = GA_MAX (dval, dval2);
    fvresult = GA_MAX (fval, fval2);
    lvresult = GA_MAX (lval, lval2);
    tmp  = GA_MAX(GA_ABS(dcval.real),GA_ABS(dcval.imag));
    tmp2 = GA_MAX(GA_ABS(dcval2.real),GA_ABS(dcval2.imag));
    tmp  = GA_MAX(tmp,tmp2);
    dcvresult.real = dcval.real;
    dcvresult.imag = dcval.imag;
    if (tmp != 0.0) {
      tmp = ((double)1.0)/tmp;
      x1 = dcval.real*tmp;
      x2 = dcval.imag*tmp;
      x3 = dcval2.real*tmp;
      x4 = dcval2.imag*tmp;
      tmp = x1*x1 + x2*x2;
      tmp2 = x3*x3 + x4*x4;
      if (tmp2 > tmp) {
        dcvresult.real = dcval2.real;
        dcvresult.imag = dcval2.imag;
      }
    }
    NGA_Fill_patch (g_d, lo, hi, vresult);
    if (type == C_DCPL) {
      needs_scaled_result = 1;
      NGA_Fill_patch (g_f, lo, hi, val3);
      NGA_Fill_patch (g_g, lo, hi, val4);
      GA_Elem_maximum_patch (g_f, lo, hi, g_g, lo, hi, g_h, lo, hi);
      tmp  = GA_MAX(GA_ABS(dcval3.real),GA_ABS(dcval3.imag));
      tmp2 = GA_MAX(GA_ABS(dcval4.real),GA_ABS(dcval4.imag));
      tmp  = GA_MAX(tmp,tmp2);
      dcvresult.real = dcval3.real;
      dcvresult.imag = dcval3.imag;
      if (tmp != 0.0) {
        tmp = ((double)1.0)/tmp;
        x1 = dcval3.real*tmp;
        x2 = dcval3.imag*tmp;
        x3 = dcval4.real*tmp;
        x4 = dcval4.imag*tmp;
        tmp = x1*x1 + x2*x2;
        tmp2 = x3*x3 + x4*x4;
        if (tmp2 > tmp) {
          dcvresult.real = dcval4.real;
          dcvresult.imag = dcval4.imag;
        }
      }
      NGA_Fill_patch (g_i, lo, hi, vresult);
      NGA_Fill_patch (g_k, lo, hi, &dcval6);
      NGA_Fill_patch (g_l, lo, hi, &dcval7);
      GA_Elem_maximum_patch (g_k, lo, hi, g_l, lo, hi, g_m, lo, hi);
      tmp  = GA_MAX(GA_ABS(dcval6.real),GA_ABS(dcval6.imag));
      tmp2 = GA_MAX(GA_ABS(dcval7.real),GA_ABS(dcval7.imag));
      tmp  = GA_MAX(tmp,tmp2);
      dcvresult.real = dcval6.real;
      dcvresult.imag = dcval6.imag;
      if (tmp != 0.0) {
        tmp = ((double)1.0)/tmp;
        x1 = dcval6.real*tmp;
        x2 = dcval6.imag*tmp;
        x3 = dcval7.real*tmp;
        x4 = dcval7.imag*tmp;
        tmp = x1*x1 + x2*x2;
        tmp2 = x3*x3 + x4*x4;
        if (tmp2 > tmp) {
          dcvresult.real = dcval7.real;
          dcvresult.imag = dcval7.imag;
        }
      }
      NGA_Fill_patch (g_n, lo, hi, vresult);
    }
    break;
    case OP_ELEM_MIN:
    if (me == 0)
      printf ("Testing GA_Elem_minimum...");
    NGA_Fill_patch (g_b, lo, hi, val2);
    GA_Elem_minimum_patch (g_a, lo, hi, g_b, lo, hi, g_c, lo, hi);
    ivresult = GA_MIN (ival, ival2);
    dvresult = GA_MIN (dval, dval2);
    fvresult = GA_MIN (fval, fval2);
    lvresult = GA_MIN (lval, lval2);
    tmp  = GA_MAX(GA_ABS(dcval.real),GA_ABS(dcval.imag));
    tmp2 = GA_MAX(GA_ABS(dcval2.real),GA_ABS(dcval2.imag));
    tmp  = GA_MAX(tmp,tmp2);
    dcvresult.real = dcval.real;
    dcvresult.imag = dcval.imag;
    if (tmp != 0.0) {
      tmp = 1.0/tmp;
      x1 = dcval.real*tmp;
      x2 = dcval.imag*tmp;
      x3 = dcval2.real*tmp;
      x4 = dcval2.imag*tmp;
      tmp = x1*x1 + x2*x2;
      tmp2 = x3*x3 + x4*x4;
      if (tmp2 < tmp) {
        dcvresult.real = dcval2.real;
        dcvresult.imag = dcval2.imag;
      }
    }
    NGA_Fill_patch (g_d, lo, hi, vresult);
    if (type == C_SCPL) {
      needs_scaled_result = 1;
      NGA_Fill_patch (g_f, lo, hi, val3);
      NGA_Fill_patch (g_g, lo, hi, val4);
      GA_Elem_minimum_patch (g_f, lo, hi, g_g, lo, hi, g_h, lo, hi);
      tmpf  = GA_MAX(GA_ABS(fcval3.real),GA_ABS(fcval3.imag));
      tmp2f = GA_MAX(GA_ABS(fcval4.real),GA_ABS(fcval4.imag));
      tmpf  = GA_MAX(tmpf,tmp2f);
      fcvresult.real = fcval3.real;
      fcvresult.imag = fcval3.imag;
      if (tmpf != 0.0) {
        tmpf = ((float)1.0)/tmpf;
        fx1 = fcval3.real*tmpf;
        fx2 = fcval3.imag*tmpf;
        fx3 = fcval4.real*tmpf;
        fx4 = fcval4.imag*tmpf;
        tmpf = fx1*fx1 + fx2*fx2;
        tmp2f = fx3*fx3 + fx4*fx4;
        if (tmp2f < tmpf) {
          fcvresult.real = fcval4.real;
          fcvresult.imag = fcval4.imag;
        }
      }
      NGA_Fill_patch (g_i, lo, hi, vresult);
      NGA_Fill_patch (g_k, lo, hi, &dcval6);
      NGA_Fill_patch (g_l, lo, hi, &dcval7);
      GA_Elem_minimum_patch (g_k, lo, hi, g_l, lo, hi, g_m, lo, hi);
      tmpf  = GA_MAX(GA_ABS(fcval6.real),GA_ABS(fcval6.imag));
      tmp2f = GA_MAX(GA_ABS(fcval7.real),GA_ABS(fcval7.imag));
      tmpf  = GA_MAX(tmpf,tmp2f);
      fcvresult.real = fcval6.real;
      fcvresult.imag = fcval6.imag;
      if (tmpf != 0.0) {
        tmpf = ((float)1.0)/tmpf;
        fx1 = fcval6.real*tmpf;
        fx2 = fcval6.imag*tmpf;
        fx3 = fcval7.real*tmpf;
        fx4 = fcval7.imag*tmpf;
        tmpf = fx1*fx1 + fx2*fx2;
        tmp2f = fx3*fx3 + fx4*fx4;
        if (tmp2f < tmpf) {
          fcvresult.real = fcval7.real;
          fcvresult.imag = fcval7.imag;
        }
      }
      NGA_Fill_patch (g_n, lo, hi, vresult);
    }
    break;
    case OP_STEP_MAX:
    if (me == 0)
      printf ("Testing GA_Step_max...");
    if (type != C_DCPL || type != C_SCPL) {
      /*NGA_Fill_patch (g_b, lo, hi, val2);*/
      GA_Abs_value_patch (g_b, lo, hi);
      GA_Step_max_patch (g_b, lo, hi, g_j, lo, hi, resultx);
      /*
      printf(" GA_Stepmax_patch type = %d, resultx = %le\n",type,resultx);
      fflush(stdout);
       */
      /* 
         It would be more robust to use GA_Elem_min_patch 
         here to determine the minimum g_j value, but for
         now we set it to -2.
       */
      aresulti = ((int)(GA_ABS(ival2)/GA_ABS(ival))) - resulti;
      aresultd = GA_ABS(dval2/dval) - resultd;
      aresultf = ((float)GA_ABS(fval2/fval)) - resultf;
      aresultl = ((long)(GA_ABS(lval2)/GA_ABS(lval))) - resultl;
    }
    break;

    case OP_STEP_BOUND_INFO:
    if (me == 0)
      printf ("Testing GA_Step_bound_info...");
    if (type != C_DCPL || type != C_SCPL) {
      /*NGA_Fill_patch (g_b, lo, hi, val2);*/
      GA_Abs_value_patch (g_b, lo, hi);
      GA_Abs_value_patch (g_a, lo, hi);
      /*GA_Abs_value_patch (g_j, lo, hi);*/
      NGA_Fill_patch(g_c, lo, hi, val5);
      GA_Step_bound_info_patch (g_b,lo,hi, g_j,lo,hi, g_a,lo,hi, g_c,lo,hi, boundminx,wolfeminx,boundmaxx);
      /*
      printf(" GA_Stepmax2_patch type = %d, resultx = %le\n",type,resultx);
      fflush(stdout);
       */
      /* 
         This is currently hardwired. would need to change if 
         val, val2 or val5 change.
       */
      switch (type)
      {
        case C_INT:
          awolfemini = ((int)(((int)1)/((int)2))) - wolfemini;
          aboundmini = ((int)(((int)1)/((int)2))) - boundmini;
          aboundmaxi = (int)GA_INFINITY_I - boundmaxi;
          break;
        case C_DBL:
          awolfemind = (((double)1.0)/((double)2.0)) - wolfemind;
          aboundmind = (((double)1.0)/((double)2.0)) - boundmind;
          aboundmaxd = (double)GA_INFINITY_D - boundmaxd;
          break;
        case C_FLOAT:
          awolfeminf = (((float)1.0)/((float)2.0))   - wolfeminf;
          aboundminf = (((float)1.0)/((float)2.0))   - boundminf;
          aboundmaxf = (float)GA_INFINITY_F - boundmaxf;
          break;
        case C_LONG:
          awolfeminl = ((long)(((long)1)/((long)2))) - wolfeminl; 
          aboundminl = ((long)(((long)1)/((long)2))) - boundminl; 
          aboundmaxl = (long)GA_INFINITY_L - boundmaxl;
          break;
        default:
          GA_Error ("GA_step_bound_info wrong data type.", type);
      }
    }
    break;
    default:
    GA_Error ("test_function: wrong operation.", OP);
  }
  switch (type)
  {
    case C_INT:
      alpha = &ai;
      beta = &bi;
      break;
    case C_DCPL:
      alpha = &adc;
      beta = &bdc;
      break;

    case C_SCPL:
      alpha = &afc;
      beta = &bfc;
      break;

    case C_DBL:
      alpha = &ad;
      beta = &bd;
      break;
    case C_FLOAT:
      alpha = &af;
      beta = &bf;
      break;
    case C_LONG:
      alpha = &al;
      beta = &bl;
      break;
    default:
      GA_Error ("wrong data type.", type);
  }

  if (OP < 4) {
    /* 
       Binary operation. 
     */
    NGA_Add_patch (alpha, g_c, lo, hi, beta, g_d, lo, hi, g_e, lo, hi);
    if (needs_scaled_result == 1) {
      NGA_Add_patch (alpha, g_h, lo, hi, beta, g_i, lo, hi, g_j, lo, hi);
      NGA_Add_patch (alpha, g_m, lo, hi, beta, g_n, lo, hi, g_n, lo, hi);
    }
  } else {

    /*
    Unary operation.
     */
    if (OP < 7) {
      NGA_Add_patch (alpha, g_a, lo, hi, beta, g_d, lo, hi, g_e, lo, hi);
      if (needs_scaled_result == 1) {
        NGA_Add_patch (alpha, g_f, lo, hi, beta, g_i, lo, hi, g_j, lo, hi);
        NGA_Add_patch (alpha, g_k, lo, hi, beta, g_n, lo, hi, g_n, lo, hi);
      }
    }
    /* 
       Else it was a reduction operation (one of the step_max functions).
     */
  }

  switch (type)
  {
    case C_INT:
      max = &imax;
      break;
    case C_DCPL:
      max = &dcmax;
      max2 = &dcmax2;
      max3 = &dcmax3;
      break;
    case C_SCPL:
      max = &fcmax;
      max2 = &fcmax2;
      max3 = &fcmax3;
      break;
    case C_DBL:
      max = &dmax;
      break;
    case C_FLOAT:
      max = &fmax;
      break;
    case C_LONG:
      max = &lmax;
      break;
    default:
      GA_Error ("wrong data type.", type);
  }

  /*  
      for unary and binary operators extract the maximum difference between
      computed and correct solutions.
   */
  if (OP < 7) {
    NGA_Select_elem (g_e, "max", max, index);
    if (needs_scaled_result == 1) {
      NGA_Select_elem (g_j, "max", max2, index2);
      NGA_Select_elem (g_n, "max", max3, index3);
    }
  }
  /*  NGA_Select_elem (g_e, "min", min, index);*/

  if (OP < 7) {
    /* 
       Binary or Unary operators.
     */
    switch (type)
    {
      double r, im;
      float rf, imf;
      case C_INT:
      /*      result = (int)(imax - imin);*/
      result = imax;
      if (result != (int)0) result = 1;
      break;
      case C_DCPL:
      /*
      r = dcmax.real - dcmin.real;
      im = dcmax.imag - dcmin.imag;
       */
      r = dcmax.real;
      im = dcmax.imag;
      if ((GA_ABS(r) + GA_ABS(im)) < THRESH) {
        result = 0;
      } else {
        result = 1;
      }
      if (needs_scaled_result == 1) {
        result2 = 0;
        r = dcmax2.real;
        im = dcmax2.imag;
        if ((GA_ABS(r) + GA_ABS(im)) < THRESH) {
          result2 = 0;
        } else {
          result2 = 1;
        }
        r = dcmax3.real;
        im = dcmax3.imag;
        if ((GA_ABS(r) + GA_ABS(im)) < THRESH) {
          result3 = 0;
        } else {
          result3 = 1;
        }
        result = result | result2 | result3;
      }
      break;
      case C_SCPL:
      /*
      rf = fcmax.real - fcmin.real;
      imf = fcmax.imag - fcmin.imag;
       */
      rf = fcmax.real;
      imf = fcmax.imag;
      if ((GA_ABS(rf) + GA_ABS(imf)) == (float)0.0) {
        result = 0;
      } else {
        result = 1;
      }
      if (needs_scaled_result == 1) {
        result2 = 0;
        rf = fcmax2.real;
        imf = fcmax2.imag;
        if ((GA_ABS(rf) + GA_ABS(imf)) == (float)0.0) {
          result2 = 0;
        } else {
          result2 = 1;
        }
        rf = fcmax3.real;
        imf = fcmax3.imag;
        if ((GA_ABS(rf) + GA_ABS(imf)) == (float)0.0) {
          result3 = 0;
        } else {
          result3 = 1;
        }
        result = result | result2 | result3;
      }
      break;
      case C_DBL:
      if (dmax == (double)0.0) {
        result = 0;
      } else {
        result = 1;
      }
      break;
      case C_FLOAT:
      if (fmax == (float)0.0) {
        result = 0;
      } else {
        result = 1;
      }
      break;
      case C_LONG:
      if (lmax == (long)0) {
        result = 0;
      } else {
        result = 1;
      }
      break;
      default:
      GA_Error ("wrong data type.", type);
    }
  } else {
    /*
    A reduction operation, Step_max or Step_bound_info.
     */
    if (type == C_DCPL || type == C_SCPL) {
      result = 0;
    } else {
      if (OP == OP_STEP_MAX) {
        /* Step_max */
        switch (type) 
        {
          case C_INT:
            if (aresulti == 0) {
              result = 0;
            } else {
              result = 1;
            }
            break;
          case C_DBL:
            if (aresultd == (double)0.0) {
              result = 0;
            } else {
              result = 1;
            }
            break;
          case C_FLOAT:
            if (aresultf == (float)0.0) {
              result = 0;
            } else {
              result = 1;
            }
            break;
          case C_LONG:
            if (aresultl == (long)0) {
              result = 0;
            } else {
              result = 1;
            }
            break;
          default:
            GA_Error ("Stepmax op, wrong data type.", type);
        }
      } else {
        /* OP = 8 so Step_bound_info */
        switch (type) 
        {
          case C_INT:
            if (awolfemini == 0) {
              result = 0;
            } else {
              result = 1;
            }
            if (aboundmini == 0) {
              result2 = 0;
            } else {
              result2 = 1;
            }
            if (aboundmaxi == 0) {
              result3 = 0;
            } else {
              result3 = 1;
            }
            break;
          case C_DBL:
            if (awolfemind == ((double)0.0)) {
              result = 0;
            } else {
              result = 1;
            }
            if (aboundmind == ((double)0.0)) {
              result2 = 0;
            } else {
              result2 = 1;
            }
            if (aboundmaxd == ((double)0.0)) {
              result3 = 0;
            } else {
              result3 = 1;
            }
            break;
          case C_FLOAT:
            if (awolfeminf == ((float)0.0)) {
              result = 0;
            } else {
              result = 1;
            }
            if (aboundminf == ((float)0.0)) {
              result2 = 0;
            } else {
              result2 = 1;
            }
            if (aboundmaxf == ((float)0.0)) {
              result3 = 0;
            } else {
              result3 = 1;
            }
            break;
          case C_LONG:
            if (awolfeminl == ((long)0)) {
              result = 0;
            } else {
              result = 1;
            }
            if (aboundminl == ((long)0)) {
              result2 = 0;
            } else {
              result2 = 1;
            }
            if (aboundmaxl == ((long)0)) {
              result3 = 0;
            } else {
              result3 = 1;
            }
            break;
          default:
            GA_Error ("Stepmax op, wrong data type.", type);
        }
        result = result | result2 | result3;
      }
    }
  }
  if (me == 0)
  {
    if (MISMATCHED (result, 0)) {
      printf ("is not ok\n");
      GA_Error("aborting", 1);
    } else {
      printf ("is ok.\n");
    }
  }

  /*
  NGA_Print_patch(g_a, lo, hi, 1);
  NGA_Print_patch(g_d, lo, hi, 1);
  NGA_Print_patch(g_e, lo, hi, 1);
   */

  GA_Destroy (g_a);
  GA_Destroy (g_b);
  GA_Destroy (g_c);
  GA_Destroy (g_d);
  GA_Destroy (g_e);
  GA_Destroy (g_f);
  GA_Destroy (g_g);
  GA_Destroy (g_h);
  GA_Destroy (g_i);
  GA_Destroy (g_j);
  GA_Destroy (g_k);
  GA_Destroy (g_l);
  GA_Destroy (g_m);
  GA_Destroy (g_n);

  return ok;
}

int
main (argc, argv)
  int argc;
  char **argv;
{
  int heap = 20000, stack = 20000;
  int me, nproc;
  int d, op;
  int ok = 1;

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
  if (!MA_init (C_DBL, stack, heap))
    GA_Error ("MA_init failed", stack + heap);    /* initialize memory allocator */


  /* op = 8;*/
  for (op = 0; op < 9; op++)
  {
    /*for (d = 1; d < 2; d++)*/
    for (d = 1; d < 4; d++)
    {
      if (me == 0) 
        printf ("\n\ndim =%d\n\n", d);
      if (me == 0) 
        printf ("\ndata type: INT\t\t");
      ok = test_fun (C_INT, d, op);
      if (me == 0) 
        printf ("\ndata type: double\t");
      ok = test_fun (C_DBL, d, op);
      if (me == 0)
        printf ("\ndata type: float\t");
      ok = test_fun (C_FLOAT, d, op);
      if (me == 0)
        printf ("\ndata type: long\t\t");
      ok = test_fun (C_LONG, d, op);
      if (op < 7) {
        if (me == 0)
          printf ("\ndata type: double complex\t");
        ok = test_fun (C_DCPL, d, op);
       }
    }
  }

  if (me==0) printf("\nAll tests successful\n");

  GA_Terminate();

  MP_FINALIZE();

  return 0;
}

/*\ FILL IN ARRAY WITH Varying VALUEs. (from 0 to product of dims-1).
  For complex arrays make the real and imaginary parts equal.
  \*/
void nga_vfill_patch(Integer *g_a, Integer *lo, Integer *hi)
{
  Integer i, j;
  Integer ndim, dims[MAXDIM], type;
  Integer loA[MAXDIM], hiA[MAXDIM], ld[MAXDIM];
  void *data_ptr;
  Integer idx, n1dim;
  Integer bvalue[MAXDIM], bunit[MAXDIM], baseld[MAXDIM];
  Integer me= pnga_nodeid();
  int local_sync_begin,local_sync_end;

  local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
  _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
  if(local_sync_begin)GA_Sync(); 

  GA_PUSH_NAME("nga_vfill_patch");

  pnga_inquire(*g_a,  &type, &ndim, dims);

  /* get limits of VISIBLE patch */ 
  pnga_distribution(*g_a, me, loA, hiA);

  /*  determine subset of my local patch to access  */
  /*  Output is in loA and hiA */
  if(pnga_patch_intersect(lo, hi, loA, hiA, ndim)){

    /* get data_ptr to corner of patch */
    /* ld are leading dimensions INCLUDING ghost cells */
    pnga_access_ptr(*g_a, loA, hiA, &data_ptr, ld);

    /* number of n-element of the first dimension */
    n1dim = 1; for(i=1; i<ndim; i++) n1dim *= (hiA[i] - loA[i] + 1);

    /* calculate the destination indices */
    bvalue[0] = 0; bvalue[1] = 0; bunit[0] = 1; bunit[1] = 1;
    /* baseld[0] = ld[0]
     * baseld[1] = ld[0] * ld[1]
     * baseld[2] = ld[0] * ld[1] * ld[2] .....
     */
    baseld[0] = ld[0]; baseld[1] = baseld[0] *ld[1];
    for(i=2; i<ndim; i++) {
      bvalue[i] = 0;
      bunit[i] = bunit[i-1] * (hiA[i-1] - loA[i-1] + 1);
      baseld[i] = baseld[i-1] * ld[i];
    }

    switch (type){
      case C_INT:
        for(i=0; i<n1dim; i++) {
          idx = 0;
          for(j=1; j<ndim; j++) {
            idx += bvalue[j] * baseld[j-1];
            if(((i+1) % bunit[j]) == 0) bvalue[j]++;
            if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
          }
          for(j=0; j<(hiA[0]-loA[0]+1); j++)
            ((int *)data_ptr)[idx+j] = (int)(idx+j);
        }
        break;
      case C_DCPL:
        for(i=0; i<n1dim; i++) {
          idx = 0;
          for(j=1; j<ndim; j++) {
            idx += bvalue[j] * baseld[j-1];
            if(((i+1) % bunit[j]) == 0) bvalue[j]++;
            if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
          }

          for(j=0; j<(hiA[0]-loA[0]+1); j++) {
            ((DoubleComplex *)data_ptr)[idx+j].real = (double)(idx+j);
            ((DoubleComplex *)data_ptr)[idx+j].imag = (double)(idx+j);
          }
        }

        break;
      case C_SCPL:
        for(i=0; i<n1dim; i++) {
          idx = 0;
          for(j=1; j<ndim; j++) {
            idx += bvalue[j] * baseld[j-1];
            if(((i+1) % bunit[j]) == 0) bvalue[j]++;
            if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
          }

          for(j=0; j<(hiA[0]-loA[0]+1); j++) {
            ((SingleComplex *)data_ptr)[idx+j].real = (float)(idx+j);
            ((SingleComplex *)data_ptr)[idx+j].imag = (float)(idx+j);
          }
        }

        break;
      case C_DBL:
        for(i=0; i<n1dim; i++) {
          idx = 0;
          for(j=1; j<ndim; j++) {
            idx += bvalue[j] * baseld[j-1];
            if(((i+1) % bunit[j]) == 0) bvalue[j]++;
            if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
          }

          for(j=0; j<(hiA[0]-loA[0]+1); j++) 
            ((double*)data_ptr)[idx+j] = (double)(idx+j);
        }
        break;
      case C_FLOAT:
        for(i=0; i<n1dim; i++) {
          idx = 0;
          for(j=1; j<ndim; j++) {
            idx += bvalue[j] * baseld[j-1];
            if(((i+1) % bunit[j]) == 0) bvalue[j]++;
            if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
          }
          for(j=0; j<(hiA[0]-loA[0]+1); j++)
            ((float *)data_ptr)[idx+j] = (float)(idx+j);
        }
        break;     
      case C_LONG:
        for(i=0; i<n1dim; i++) {
          idx = 0;
          for(j=1; j<ndim; j++) {
            idx += bvalue[j] * baseld[j-1];
            if(((i+1) % bunit[j]) == 0) bvalue[j]++;
            if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
          }
          for(j=0; j<(hiA[0]-loA[0]+1); j++)
            ((long *)data_ptr)[idx+j] = (long)(idx+j);
        } 
        break;                          
      default: GA_Error(" wrong data type ",type);
    }

    /* release access to the data */
    pnga_release_update(*g_a, loA, hiA);
  }
  GA_POP_NAME;
  if(local_sync_end)GA_Sync();
}
/*\ Utility function to actually set positive/negative values
  \*/
void ngai_do_pnfill_patch(Integer type, Integer ndim, Integer *loA, Integer *hiA,
    Integer *ld, void *data_ptr)
{
  Integer i, j;
  Integer idx, n1dim;
  Integer bvalue[MAXDIM], bunit[MAXDIM], baseld[MAXDIM];
  /* number of n-element of the first dimension */
  n1dim = 1; for(i=1; i<ndim; i++) n1dim *= (hiA[i] - loA[i] + 1);

  /* calculate the destination indices */
  bvalue[0] = 0; bvalue[1] = 0; bunit[0] = 1; bunit[1] = 1;
  /* baseld[0] = ld[0]
   * baseld[1] = ld[0] * ld[1]
   * baseld[2] = ld[0] * ld[1] * ld[2] .....
   */
  baseld[0] = ld[0]; baseld[1] = baseld[0] *ld[1];
  for(i=2; i<ndim; i++) {
    bvalue[i] = 0;
    bunit[i] = bunit[i-1] * (hiA[i-1] - loA[i-1] + 1);
    baseld[i] = baseld[i-1] * ld[i];
  }

  switch (type){
    case C_INT:
      for(i=0; i<n1dim; i++) {
        idx = 0;
        for(j=1; j<ndim; j++) {
          idx += bvalue[j] * baseld[j-1];
          if(((i+1) % bunit[j]) == 0) bvalue[j]++;
          if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
        }

        for(j=0; j<(hiA[0]-loA[0]+1); j++)
          ((int *)data_ptr)[idx+j] = (int)(((idx+j)&3)-2);
      }
      break;
    case C_DCPL:
      for(i=0; i<n1dim; i++) {
        idx = 0;
        for(j=1; j<ndim; j++) {
          idx += bvalue[j] * baseld[j-1];
          if(((i+1) % bunit[j]) == 0) bvalue[j]++;
          if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
        }

        for(j=0; j<(hiA[0]-loA[0]+1); j++) {
          ((DoubleComplex *)data_ptr)[idx+j].real = (double)(((idx+j)&3)-2);
          ((DoubleComplex *)data_ptr)[idx+j].imag = (double)(((idx+j)&3)-2);
        }
      }
      break;
    case C_SCPL:
      for(i=0; i<n1dim; i++) {
        idx = 0;
        for(j=1; j<ndim; j++) {
          idx += bvalue[j] * baseld[j-1];
          if(((i+1) % bunit[j]) == 0) bvalue[j]++;
          if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
        }

        for(j=0; j<(hiA[0]-loA[0]+1); j++) {
          ((SingleComplex *)data_ptr)[idx+j].real = (float)(((idx+j)&3)-2);
          ((SingleComplex *)data_ptr)[idx+j].imag = (float)(((idx+j)&3)-2);
        }
      }
      break;
    case C_DBL:
      for(i=0; i<n1dim; i++) {
        idx = 0;
        for(j=1; j<ndim; j++) {
          idx += bvalue[j] * baseld[j-1];
          if(((i+1) % bunit[j]) == 0) bvalue[j]++;
          if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
        }

        for(j=0; j<(hiA[0]-loA[0]+1); j++) 
          ((double*)data_ptr)[idx+j] = (double)(((idx+j)&3)-2);
      }
      break;
    case C_FLOAT:
      for(i=0; i<n1dim; i++) {
        idx = 0;
        for(j=1; j<ndim; j++) {
          idx += bvalue[j] * baseld[j-1];
          if(((i+1) % bunit[j]) == 0) bvalue[j]++;
          if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
        }
        for(j=0; j<(hiA[0]-loA[0]+1); j++)
          ((float *)data_ptr)[idx+j] = (float)(((idx+j)&3)-2);
      }
      break;     
    case C_LONG:
      for(i=0; i<n1dim; i++) {
        idx = 0;
        for(j=1; j<ndim; j++) {
          idx += bvalue[j] * baseld[j-1];
          if(((i+1) % bunit[j]) == 0) bvalue[j]++;
          if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
        }
        for(j=0; j<(hiA[0]-loA[0]+1); j++)
          ((long *)data_ptr)[idx+j] = (long)(((idx+j)&3)-2);
      } 
      break;                          
    default: GA_Error(" wrong data type ",type);
  }

}

/*\ FILL IN ARRAY WITH Varying positive and negative VALUEs. 
  (from -2 to 1).
  For complex arrays make the real and imaginary parts equal.
  \*/
void nga_pnfill_patch(Integer *g_a, Integer *lo, Integer *hi)
{
  Integer i;
  Integer ndim, dims[MAXDIM], type;
  Integer loA[MAXDIM], hiA[MAXDIM], ld[MAXDIM];
  void *data_ptr;
  Integer me= pnga_nodeid();
  Integer num_blocks;
  int local_sync_begin,local_sync_end;

  local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
  _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
  if(local_sync_begin)GA_Sync(); 

  GA_PUSH_NAME("nga_pnfill_patch");

  pnga_inquire(*g_a,  &type, &ndim, dims);

  num_blocks = pnga_total_blocks(*g_a);

  if (num_blocks < 0) {
    /* get limits of VISIBLE patch */ 
    pnga_distribution(*g_a, me, loA, hiA);

    /*  determine subset of my local patch to access  */
    /*  Output is in loA and hiA */
    if(pnga_patch_intersect(lo, hi, loA, hiA, ndim)){

      /* get data_ptr to corner of patch */
      /* ld are leading dimensions INCLUDING ghost cells */
      pnga_access_ptr(*g_a, loA, hiA, &data_ptr, ld);

      ngai_do_pnfill_patch(type, ndim, loA, hiA, ld, data_ptr);

      /* release access to the data */
      pnga_release_update(*g_a, loA, hiA);
    }
  } else {
    Integer offset, j, jtmp, chk;
    Integer loS[MAXDIM];
    Integer nproc = pnga_nnodes();
    /* using simple block-cyclic data distribution */
    if (!pnga_uses_proc_grid(*g_a)){
      for (i=me; i<num_blocks; i += nproc) {
        /* get limits of patch */
        pnga_distribution(*g_a, i, loA, hiA);

        /* loA is changed by pnga_patch_intersect, so
           save a copy */
        for (j=0; j<ndim; j++) {
          loS[j] = loA[j];
        }

        /*  determine subset of my local patch to access  */
        /*  Output is in loA and hiA */
        if(pnga_patch_intersect(lo, hi, loA, hiA, ndim)){

          /* get data_ptr to corner of patch */
          /* ld are leading dimensions for block */
          pnga_access_block_ptr(*g_a, i, &data_ptr, ld);

          /* Check for partial overlap */
          chk = 1;
          for (j=0; j<ndim; j++) {
            if (loS[j] < loA[j]) {
              chk=0;
              break;
            }
          }
          if (!chk) {
            /* Evaluate additional offset for pointer */
            offset = 0;
            jtmp = 1;
            for (j=0; j<ndim-1; j++) {
              offset += (loA[j]-loS[j])*jtmp;
              jtmp *= ld[j];
            }
            offset += (loA[ndim-1]-loS[ndim-1])*jtmp;
            switch (type){
              case C_INT:
                data_ptr = (void*)((int*)data_ptr + offset);
                break;
              case C_DCPL:
                data_ptr = (void*)((double*)data_ptr + 2*offset);
                break;
              case C_SCPL:
                data_ptr = (void*)((float*)data_ptr + 2*offset);
                break;
              case C_DBL:
                data_ptr = (void*)((double*)data_ptr + offset);
                break;
              case C_FLOAT:
                data_ptr = (void*)((float*)data_ptr + offset);
                break;
              case C_LONG:
                data_ptr = (void*)((long*)data_ptr + offset);
                break;
              default: GA_Error(" wrong data type ",type);
            }
          }
          /* fill in patch */
          ngai_do_pnfill_patch(type, ndim, loA, hiA, ld, data_ptr);

          /* release access to the data */
          pnga_release_update_block(*g_a, i);
        }
      }
    } else {
      /* using scalapack block-cyclic data distribution */
      Integer proc_index[MAXDIM], index[MAXDIM];
      Integer topology[MAXDIM];
      Integer blocks[MAXDIM], block_dims[MAXDIM];
      pnga_get_proc_index(*g_a, me, proc_index);
      pnga_get_proc_index(*g_a, me, index);
      pnga_get_block_info(*g_a, blocks, block_dims);
      pnga_get_proc_grid(*g_a, topology);
      while (index[ndim-1] < blocks[ndim-1]) {
        /* find bounding coordinates of block */
        chk = 1;
        for (i = 0; i < ndim; i++) {
          loA[i] = index[i]*block_dims[i]+1;
          hiA[i] = (index[i] + 1)*block_dims[i];
          if (hiA[i] > dims[i]) hiA[i] = dims[i];
          if (hiA[i] < loA[i]) chk = 0;
        }

        /* loA is changed by pnga_patch_intersect, so
           save a copy */
        for (j=0; j<ndim; j++) {
          loS[j] = loA[j];
        }

        /*  determine subset of my local patch to access  */
        /*  Output is in loA and hiA */
        if(pnga_patch_intersect(lo, hi, loA, hiA, ndim)){

          /* get data_ptr to corner of patch */
          /* ld are leading dimensions for block */
          pnga_access_block_grid_ptr(*g_a, index, &data_ptr, ld);

          /* Check for partial overlap */
          chk = 1;
          for (j=0; j<ndim; j++) {
            if (loS[j] < loA[j]) {
              chk=0;
              break;
            }
          }
          if (!chk) {
            /* Evaluate additional offset for pointer */
            offset = 0;
            jtmp = 1;
            for (j=0; j<ndim-1; j++) {
              offset += (loA[j]-loS[j])*jtmp;
              jtmp *= ld[j];
            }
            offset += (loA[ndim-1]-loS[ndim-1])*jtmp;
            switch (type){
              case C_INT:
                data_ptr = (void*)((int*)data_ptr + offset);
                break;
              case C_DCPL:
                data_ptr = (void*)((double*)data_ptr + 2*offset);
                break;
              case C_SCPL:
                data_ptr = (void*)((float*)data_ptr + 2*offset);
                break;
              case C_DBL:
                data_ptr = (void*)((double*)data_ptr + offset);
                break;
              case C_FLOAT:
                data_ptr = (void*)((float*)data_ptr + offset);
                break;
              case C_LONG:
                data_ptr = (void*)((long*)data_ptr + offset);
                break;
              default: GA_Error(" wrong data type ",type);
            }
          }
          /* fill in patch */
          ngai_do_pnfill_patch(type, ndim, loA, hiA, ld, data_ptr);

          /* release access to the data */
          pnga_release_update_block_grid(*g_a, index);
        }
        /* increment index to get next block on processor */
        index[0] += topology[0];
        for (i = 0; i < ndim; i++) {
          if (index[i] >= blocks[i] && i<ndim-1) {
            index[i] = proc_index[i];
            index[i+1] += topology[i+1];
          }
        }
      }
    }
  }
  GA_POP_NAME;
  if(local_sync_end)GA_Sync();
}


