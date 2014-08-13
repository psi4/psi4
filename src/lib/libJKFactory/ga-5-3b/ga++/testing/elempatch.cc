#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <iostream>
#include <cstdio>
#include <cmath>
using namespace std;
#include "ga++.h"








#define N     10    // First dimension
#define NDIM  4     // Number of dimensions
#define BASE  0
#define PERMUTE_

#define GA_DATA_TYPE MT_F_REAL

#define MAXDIM GA_MAX_DIM
#define C_DBL MT_C_DBL
#define C_INT MT_C_INT
#define C_FLOAT MT_C_FLOAT
#define C_DCPL MT_C_DCPL
#define C_LONG MT_C_LONGINT
#define C_SCPL MT_C_SCPL
#define GA_ABS(a)   (((a) >= 0) ? (a) : (-(a)))
#define GA_MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define GA_MIN(a,b) (((a) <= (b)) ? (a) : (b))

#define THRESH 1e-5
#define MISMATCHED(x,y) GA_ABS((x)-(y))>=THRESH

#define OP_ELEM_MULT 0
#define OP_ELEM_DIV 1
#define OP_ELEM_MAX 2
#define OP_ELEM_MIN 3
#define OP_ABS 4
#define OP_ADD_CONST 5
#define OP_RECIP 6
#define MY_TYPE 2002

/* Integer _ga_lo[MAXDIM], _ga_hi[MAXDIM], _ga_work[MAXDIM];*/
#  define COPYINDEX_C2F(carr, farr, n){\
   int i; for(i=0; i< (n); i++)(farr)[n-i-1]=(Integer)(carr)[i]+1;}


int
test_fun (int type, int dim, int OP) {
  
  GA::GlobalArray *g_a, *g_b, *g_c, *g_d, *g_e;
  int me = GA_Nodeid ();
  int i;
  int dims[MAXDIM];
  int lo[MAXDIM], hi[MAXDIM];
  int index[MAXDIM];
  void *val;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;
  void *val2;
  int ival2 = -3;
  double dval2 = -3.0;
  float fval2 = -3.0;
  long lval2 = -3;
  DoubleComplex dcval2;
  int ok = 1;
  int result;
  void *min, *max;
  float fmin, fmax;
  long lmin, lmax;
  double dmin, dmax;
  DoubleComplex dcmin, dcmax;


  void *alpha, *beta;
  int ai = 1, bi = -1;
  long al = 1, bl = -1;
  float af = 1.0, bf = -1.0;
  double ad = 1.0, bd = -1.0;
  DoubleComplex adc, bdc;

  char name_A[] = "A";
  char name_B[] = "B";
  char name_C[] = "C";
  char name_D[] = "D";
  char name_E[] = "E";
  char name_max[] = "max";
  char name_min[] = "min";

  adc.real = 1.0;
  adc.imag = 0.0;
  bdc.real = -1.0;
  bdc.imag = 0.0;


  dcval.real = -sin (3.0);
  dcval.imag = -cos (3.0);
  dcval2.real = 2 * sin (3.0);
  dcval2.imag = 2 * cos (3.0);

  for (i = 0; i < dim; i++)  dims[i] = N;

  for (i = 0; i < dim; i++) {
    lo[i] = 0;
    hi[i] = N - 1;
  }
  
  g_a = GA::SERVICES.createGA (type, dim, dims, name_A, NULL);
  g_b = GA::SERVICES.createGA (g_a, name_B);
  g_c = GA::SERVICES.createGA (g_a, name_C);
  g_d = GA::SERVICES.createGA (g_a, name_D);
  g_e = GA::SERVICES.createGA (g_a, name_E);

  /*initialize  with zero */
  g_a->zero ();
  g_b->zero ();
  g_c->zero ();
  g_d->zero ();
  g_e->zero ();
  
  switch (type)
    {
    case C_INT:
      val  = (void *)&ival;
      val2 = (void *)&ival2;
      break;
    case C_DCPL:
      val  = (void *)&dcval;
      val2 = (void *)&dcval2;
      break;

    case C_DBL:
      val  = (void *)&dval;
      val2 = (void *)&dval2;
      break;
    case C_FLOAT:
      val  = (void *)&fval;
      val2 = (void *)&fval2;
      break;
    case C_LONG:
      val  = (void *)&lval;
      val2 = (void *)&lval2;
      break;
    default:
      GA::SERVICES.error ("wrong data type.", type);
    }


  g_a->fillPatch (lo, hi, val);

  switch (OP)
    {
      double tmp, tmp2;
      DoubleComplex dctemp;
    case OP_ABS:
      if (me == 0)
	printf ("Testing GA_Abs_value...");
      g_a->absValuePatch (lo, hi);
      ival = GA_ABS (ival);
      dval = GA_ABS (dval);
      fval = GA_ABS (fval);
      lval = GA_ABS (lval);
      dcval.real = dcval.real * dcval.real + dcval.imag * dcval.imag;
      dcval.imag = 0.0;
      g_d->fillPatch (lo, hi, val);
      break;
    case OP_ADD_CONST:
      if (me == 0)
	printf ("Testing GA_Add_const...");
      g_a->addConstantPatch (lo, hi, val2);
      ival = ival + ival2;
      dval = dval + dval2;
      fval = fval + fval2;
      lval = lval + lval2;
      dcval.real = dcval.real + dcval2.real;
      dcval.imag = dcval.imag + dcval2.imag;
      g_d->fillPatch (lo, hi, val);
      break;
    case OP_RECIP:
      if (me == 0)
	printf ("Testing GA_Recip...");
      g_a->recipPatch (lo, hi);
      ival = 1 / ival;
      dval = 1.0 / dval;
      fval = 1.0 / fval;
      lval = 1 / lval;
      tmp = dcval.real * dcval.real + dcval.imag * dcval.imag;
      dcval.real = dcval.real / tmp;
      dcval.imag = -dcval.imag / tmp;
      g_d->fillPatch (lo, hi, val);
      break;
    case OP_ELEM_MULT:
      if (me == 0)
	printf ("Testin GA_Elem_multiply...");
      g_b->fillPatch (lo, hi, val);
#if 0
      //g_c is different from g_a or g_b
      g_c->elemMultiplyPatch (g_a, lo, hi, g_b, lo, hi, lo, hi);
#else
      //g_c is g_b 
      g_b->elemMultiplyPatch (g_a, lo, hi, g_b, lo, hi, lo, hi);
#endif
      ival = ival * ival2;
      dval = dval * dval2;
      fval = fval * fval2;
      lval = lval * lval2;
      dctemp.real = dcval.real * dcval2.real - dcval.imag * dcval2.imag;
      dctemp.imag = dcval.real * dcval2.imag + dcval2.real * dcval.imag;
      dcval = dctemp;
      g_d->fillPatch (lo, hi, val);
      break;
    case OP_ELEM_DIV:
      if (me == 0)
	printf ("Testin GA_Elem_divide...");
      g_b->fillPatch (lo, hi, val2);
      g_c->elemDividePatch (g_a, lo, hi, g_b, lo, hi, lo, hi);
      ival = ival / ival2;
      dval = dval / dval2;
      fval = fval / fval2;
      lval = lval / lval2;
      tmp = dcval2.real * dcval2.real + dcval2.imag * dcval2.imag;
      dctemp.real =
	(dcval.real * dcval2.real + dcval.imag * dcval2.imag) / tmp;
      dctemp.imag =
	(-dcval.real * dcval2.imag + dcval2.real * dcval.imag) / tmp;
      dcval = dctemp;
      g_d->fillPatch (lo, hi, val);
      break;

    case OP_ELEM_MAX:
      if (me == 0)
	printf ("Testin GA_Elem_maximum...");
      g_b->fillPatch (lo, hi, val2);
      g_c->elemMaximumPatch (g_a, lo, hi, g_b, lo, hi, lo, hi);
      ival = GA_MAX (ival, ival2);
      dval = GA_MAX (dval, dval2);
      fval = GA_MAX (fval, fval2);
      lval = GA_MAX (lval, lval2);
      tmp = dcval.real * dcval.real + dcval.imag * dcval.imag;
      tmp2 = dcval2.real * dcval2.real + dcval2.imag * dcval2.imag;
      if (tmp2 > tmp)
	dcval = dcval2;
      g_d->fillPatch (lo, hi, val);
      break;
    case OP_ELEM_MIN:
      if (me == 0)
	printf ("Testin GA_Elem_minimum...");
      g_b->fillPatch (lo, hi, val2);
      g_c->elemMinimumPatch (g_a, lo, hi, g_b, lo, hi, lo, hi);
      ival = GA_MIN (ival, ival2);
      dval = GA_MIN (dval, dval2);
      fval = GA_MIN (fval, fval2);
      lval = GA_MIN (lval, lval2);
      tmp = dcval.real * dcval.real + dcval.imag * dcval.imag;
      tmp2 = dcval2.real * dcval2.real + dcval2.imag * dcval2.imag;
      if (tmp2 < tmp)
	dcval = dcval2;
      g_d->fillPatch (lo, hi, val);
      break;
    default:
      GA::SERVICES.error ("test_function: wrong operation.", OP);

    }
  switch (type)
    {
    case C_INT:
      alpha = (void *)&ai;
      beta  = (void *)&bi;
      break;
    case C_DCPL:
      alpha = (void *)&adc;
      beta  = (void *)&bdc;
      break;

    case C_DBL:
      alpha = (void *)&ad;
      beta  = (void *)&bd;
      break;
    case C_FLOAT:
      alpha = (void *)&af;
      beta  =(void *) &bf;
      break;
    case C_LONG:
      alpha = (void *)&al;
      beta  = (void *)&bl;
      break;
    default:
      GA::SERVICES.error ("wrong data type.", type);
    }

  if (OP < 4) 
    g_e->addPatch (alpha, g_c, lo, hi, beta, g_d, lo, hi, lo, hi);
  else
    g_e->addPatch (alpha, g_a, lo, hi, beta, g_d, lo, hi, lo, hi);

  switch (type)
    {
    case C_INT:
      max =  (void *)&lmax;
      min = (void *)&lmin;
      break;
    case C_DCPL:
      max = (void *)&dcmax;
      min = (void *)&dcmin;
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
      GA::SERVICES.error ("wrong data type.", type);
    }
  
  g_e->selectElem (name_max, max, index);
  g_e->selectElem (name_min, min, index);

  switch (type)
    {
      double r, im;
    case C_INT:
      result = lmax - lmin;
      break;
    case C_DCPL:
      r = dcmax.real - dcmin.real;
      im = dcmax.imag - dcmin.imag;
      result = (int) (GA_ABS (r) + GA_ABS (im));
      break;
    case C_DBL:
      result = (int) (dmax - dmin);
      break;
    case C_FLOAT:
      result = (int) (fmax - fmin);
      break;
    case C_LONG:
      result = (int) (lmax - lmin);
      break;
    default:
      GA::SERVICES.error ("wrong data type.", type);
    }


  if (me == 0)
    {
      if (MISMATCHED (result, 0))
	printf ("is not ok\n");
      else
	printf ("is ok.\n");
    }
  
  /*
    g_a->printPatch(lo, hi, 1);
    g_d->printPatch(lo, hi, 1);
    g_a->printPatch(lo, hi, 1);
  */

  g_a->destroy ();
  g_b->destroy ();
  g_c->destroy ();
  g_d->destroy ();
  g_e->destroy ();
  
  return ok;
}


int
main(int argc, char *argv[]) {
 
  int me, nproc;
  int heap  = 200000, stack = 200000;
  int d, op, ok = 1;

  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, 0);
  me=GA_Nodeid(); nproc=GA_Nnodes();
  cout << "Rank = " << me << " : Size = " << nproc << "\n";

  cout << "After Initialize()\n";
  
  for (op = 0; op < 7; op++) {
    for (d = 1; d < 4; d++) {
      if (me == 0)  printf ("\n\ndim =%d\n\n", d);
  
      if (me == 0)  printf ("data type: int\t\t");
      ok = test_fun (C_INT, d, op);
      
      if (me == 0)  printf ("data type: double\t");
      ok = test_fun (C_DBL, d, op);
      
      if (me == 0)  printf ("data type: float\t");
      ok = test_fun (C_FLOAT, d, op);
      
      if (me == 0)  printf ("data type: long\t\t");
      ok = test_fun (C_LONG, d, op);
      
      if (me == 0)  printf ("data type: complex\t");
      test_fun (C_DCPL, d, op);
    }
  }
  
  if(!me) cout << "Terminating\n";
  GA::Terminate();
}
