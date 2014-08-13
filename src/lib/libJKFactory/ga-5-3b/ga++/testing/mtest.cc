#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <iostream>
#include <cstdio>
#include <cmath>
using namespace std;
#include "ga++.h"


#define GA_DATA_TYPE MT_F_REAL

#define N 4			/* dimension of matrices */

#define MAXDIM GA_MAX_DIM

#define OP_SHIFT_DIAGONAL 1
#define OP_SET_DIAGONAL 2
#define OP_ADD_DIAGONAL         3
#define OP_GET_DIAGONAL 4
#define OP_NORM1	5
#define OP_NORM_INFINITY	6
#define OP_MEDIAN               7
#define OP_MEDIAN_PATCH         8
#define OP_SCALE_ROWS           9
#define OP_SCALE_COLS           10

#define GA_ABS(a)   (((a) >= 0) ? (a) : (-(a)))
#define GA_MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define GA_MIN(a,b) (((a) <= (b)) ? (a) : (b))

# define THRESH 1e-5
#define MISMATCHED(x,y) GA_ABS((x)-(y))>=THRESH


void  
test_scale_cols (GA::GlobalArray *g_a, 
		 GA::GlobalArray *g_v) {
  
  int index[MAXDIM];
  void *min, *max;
  int imin, imax;
  float fmin, fmax;
  long lmin, lmax;
  double dmin, dmax;
  DoubleComplex dcmin, dcmax;


  void *alpha, *beta;
  int ai = 1, bi = -1;
  long al = 1, bl = -1;
  float af = 1.0, bf = -1.0;
  double ad = 1.0, bd = -1.0;
  DoubleComplex adc = { 1.0, 0.0 }, bdc =
    {
      -1.0, 0.0};

  GA::GlobalArray * g_b, *g_c;
  int me = GA_Nodeid ();
  void *val;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval = { -2.0, 0.0 };
  void *val2;
  int ival2 = 4;
  double dval2 = 4.0;
  float fval2 = 4.0;
  long lval2 = 4;
  DoubleComplex dcval2 = { 4.0, 0.0 };

  int type, ndim, dims[MAXDIM];
  int vtype, vndim, vdims[MAXDIM];

  g_a->inquire (&type, &ndim, dims);
  g_v->inquire (&vtype, &vndim, vdims);

  switch (type)
    {
    case C_INT:
      alpha = (void *)&ai;
      beta = (void *)&bi;
      break;
    case C_DCPL:
      alpha = (void *)&adc;
      beta = (void *)&bdc;
      break;

    case C_DBL:
      alpha = (void *)&ad;
      beta =(void *) &bd;
      break;
    case C_FLOAT:
      alpha =(void *) &af;
      beta = (void *)&bf;
      break;
    case C_LONG:
      alpha = (void *)&al;
      beta = (void *)&bl;
      break;
    default:
      GA::SERVICES.error ((char *)"test_scale_cols:wrong data type.", type);
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
      GA::SERVICES.error ((char *)"test_scale_cols:wrong data type.", type);
    }

  if (me == 0)
    printf ("Testing GA_Scale_cols...");


  g_a->fill(val);
  g_v->fill(val);
  
  g_a->scaleCols (g_v);
  /*the result is the same same as g_b filled with val2 */
  g_b = GA::SERVICES.createGA (g_a, (char *)"B");
  g_c = GA::SERVICES.createGA (g_a, (char *)"C");

  g_b->fill(val2);

  g_c->add(alpha, g_a, beta, g_b);

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
    case C_DBL:
      max = (void *)&dmax;
      min =(void *) &dmin;
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
      GA::SERVICES.error ((char *)"test_scale_rows:wrong data type.", type);
    }

  g_c->selectElem ((char *)"max", max, index);
  g_c->selectElem ((char *)"min", min, index);


  switch (type)
    {
      double r, m;
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
      GA::SERVICES.error ((char *)"test_scale_rows:wrong data type.", type);
    }
}



void
test_scale_rows (GA::GlobalArray *g_a, 
		 GA::GlobalArray *g_v) {
  
  int index[MAXDIM];
  void *min, *max;
  int imin, imax;
  float fmin, fmax;
  long lmin, lmax;
  double dmin, dmax;
  DoubleComplex dcmin, dcmax;


  void *alpha, *beta;
  int ai = 1, bi = -1;
  long al = 1, bl = -1;
  float af = 1.0, bf = -1.0;
  double ad = 1.0, bd = -1.0;
  DoubleComplex adc = { 1.0, 0.0 };
  DoubleComplex bdc = {-1.0, 0.0 };

  GA::GlobalArray *g_b, *g_c;
  int me = GA_Nodeid ();
  void *val;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval = { -2.0, 0.0 };
  void *val2;
  int ival2 = 4;
  double dval2 = 4.0;
  float fval2 = 4.0;
  long lval2 = 4;
  DoubleComplex dcval2 = { 4.0, 0.0 };

  int type, ndim, dims[MAXDIM];
  int vtype, vndim, vdims[MAXDIM];

  g_a->inquire (&type, &ndim, dims);
  g_v->inquire (&vtype, &vndim, vdims);

  switch (type)
    {
    case C_INT:
      alpha = (void *)&ai;
      beta = (void *)&bi;
      break;
    case C_DCPL:
      alpha = (void *)&adc;
      beta = (void *)&bdc;
      break;

    case C_DBL:
      alpha = (void *)&ad;
      beta =(void *) &bd;
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
      GA::SERVICES.error ((char *)"test_scale_rows:wrong data type.", type);
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
    case C_DBL:
      val = (void *)&dval;
      val2 = (void *)&dval2;
      break;
    case C_FLOAT:
      val =(void *) &fval;
      val2 = (void *)&fval2;
      break;
    case C_LONG:
      val = (void *)&lval;
      val2 = (void *)&lval2;
      break;
    default:
      GA::SERVICES.error ((char *)"test_scale_rows:wrong data type.", type);
    }

  if (me == 0)
    printf ("Testing GA_Scale_rows...");



  g_a->fill (val);
  g_v->fill (val);

  g_a->scaleRows (g_v);
  /*the result is the same same as g_b filled with val2 */
  g_b = GA::SERVICES.createGA (g_a, (char *)"B");
  g_c = GA::SERVICES.createGA (g_a, (char *)"C");

  g_b->fill (val2);

  g_c->add (alpha, g_a, beta, g_b);
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
      min =(void *) &lmin;
      break;
    default:
      GA::SERVICES.error ((char *)"test_scale_rows:wrong data type.", type);
    }

  g_c->selectElem ((char *)"max", max, index);
  g_c->selectElem ((char *)"min", min, index);


  switch (type)
    {
      double r, m;
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
      GA::SERVICES.error ((char *)"test_scale_rows:wrong data type.", type);
    }

}

void 
test_median_patch (GA::GlobalArray * g_a, int *alo, int *ahi, 
		   GA::GlobalArray * g_b, int *blo, int *bhi,
		   GA::GlobalArray * g_c, int *clo, int *chi, 
		   GA::GlobalArray * g_m, int *mlo, int *mhi) {

  GA::GlobalArray * g_e;
  int index[MAXDIM];
  void *min, *max;
  int imin, imax;
  float fmin, fmax;
  long lmin, lmax;
  double dmin, dmax;
  DoubleComplex dcmin, dcmax;


  void *alpha, *beta;
  int ai = 1, bi = -1;
  long al = 1, bl = -1;
  float af = 1.0, bf = -1.0;
  double ad = 1.0, bd = -1.0;
  DoubleComplex adc = { 1.0, 0.0 };
  DoubleComplex bdc = {-1.0, 0.0 };

  void *val;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;

  void *val2;
  int ival2 = 6;
  double dval2 = 6.0;
  float fval2 = 6.0;
  long lval2 = 6;
  DoubleComplex dcval2;

  void *val3;
  int ival3 = 4;
  double dval3 = 4.0;
  float fval3 = 4.0;
  long lval3 = 4;
  DoubleComplex dcval3;

  int me = GA_Nodeid ();
  int type, ndim, dims[MAXDIM];

  g_a->inquire (&type, &ndim, dims);


  switch (type)
    {
    case C_INT:
      alpha = (void *)&ai;
      beta =(void *) &bi;
      break;
    case C_DCPL:
      alpha =(void *) &adc;
      beta = (void *)&bdc;
      break;

    case C_DBL:
      alpha = (void *)&ad;
      beta = (void *)&bd;
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
      GA::SERVICES.error ((char *)"test_median:wrong data type.", type);
    }

  dcval.real = -2.0;
  dcval.imag = -0.0;

  dcval2.real = 6.0;
  dcval2.imag = 0.0;


  dcval3.real = 4.0;
  dcval3.imag = 0.0;


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
      GA::SERVICES.error ((char *)"test_median:test_median:wrong data type.", 
			  type);
    }

  if (me == 0)
    printf ("Testing GA_Median_patch...");

  g_a->zero ();
  g_b->zero ();
  g_c->zero ();
  g_m->zero ();

  g_a->fillPatch (alo, ahi, val);
  g_b->fillPatch (blo, bhi, val2);
  g_c->fillPatch (alo, bhi, val3);

  g_m->medianPatch (g_a, alo, ahi, g_b, blo, bhi, g_c, clo, chi, mlo, mhi);

  /*
    The result array should        be g_c due to the value I chose: 
    val3 is the median of the three values val, val2, and val3
  */

  /* g_e = g_c - g_m */
  g_e = GA::SERVICES.createGA(g_a, (char *)"E");
  g_e->zero ();
  g_e->addPatch (alpha, g_c, clo, chi, beta, g_m, mlo, mhi, alo, ahi);

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
      GA_Error ((char *)"test_median:wrong data type.", type);
    }

  g_e->selectElem ((char *)"max", max, index);
  g_e->selectElem ((char *)"min", min, index);


  switch (type)
    {
      double r, m;
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
      GA::SERVICES.error ((char *)"test_median:wrong data type.", type);
    }


}

void 
test_median (GA::GlobalArray * g_a, GA::GlobalArray * g_b,
	     GA::GlobalArray * g_c, GA::GlobalArray * g_m)
{
  
  GA::GlobalArray *g_e;
  int index[MAXDIM];
  void *min, *max;
  int imin, imax;
  float fmin, fmax;
  long lmin, lmax;
  double dmin, dmax;
  DoubleComplex dcmin, dcmax;


  void *alpha, *beta;
  int ai = 1, bi = -1;
  long al = 1, bl = -1;
  float af = 1.0, bf = -1.0;
  double ad = 1.0, bd = -1.0;
  DoubleComplex adc = { 1.0, 0.0 };
  DoubleComplex bdc = {-1.0, 0.0 };

  void *val;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;

  void *val2;
  int ival2 = 6;
  double dval2 = 6.0;
  float fval2 = 6.0;
  long lval2 = 6;
  DoubleComplex dcval2;

  void *val3;
  int ival3 = 4;
  double dval3 = 4.0;
  float fval3 = 4.0;
  long lval3 = 4;
  DoubleComplex dcval3;

  int me = GA_Nodeid ();
  int type, ndim, dims[MAXDIM];

  g_a->inquire (&type, &ndim, dims);

  switch (type)
    {
    case C_INT:
      alpha = (void *)&ai;
      beta = (void *)&bi;
      break;
    case C_DCPL:
      alpha = (void *)&adc;
      beta =(void *) &bdc;
      break;

    case C_DBL:
      alpha = (void *)&ad;
      beta = (void *)&bd;
      break;
    case C_FLOAT:
      alpha =(void *) &af;
      beta = (void *)&bf;
      break;
    case C_LONG:
      alpha = (void *)&al;
      beta = (void *)&bl;
      break;
    default:
      GA::SERVICES.error ((char *)"test_median:wrong data type.", type);
    }

  dcval.real = -2.0;
  dcval.imag = -0.0;

  dcval2.real = 6.0;
  dcval2.imag = 0.0;


  dcval3.real = 4.0;
  dcval3.imag = 0.0;


  switch (type)
    {
    case C_INT:
      val = (void *)&ival;
      val2 =(void *) &ival2;
      val3 = (void *)&ival3;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      val2 =(void *) &dcval2;
      val3 =(void *) &dcval3;
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
      GA_Error ((char *)"test_median:test_median:wrong data type.", type);
    }

  if (me == 0)
    printf ("Testing GA_Median...");

  g_a->zero();
  g_b->zero();
  g_c->zero();
  g_m->zero();

  g_a->fill (val);
  g_b->fill (val2);
  g_c->fill (val3);

  g_m->median (g_a, g_b, g_c);
  
  /*
    The result array should        be g_c due to the value I chose: 
    val3 is the median of the three values val, val2, and val3
  */

  /* g_e = g_c - g_m */
  g_e = GA::SERVICES.createGA (g_a, (char *)"E");
  g_e->add (alpha, g_c, beta, g_m);
  
  switch (type)
    {
    case C_INT:
      max =(void *) &imax;
      min = (void *)&imin;
      break;
    case C_DCPL:
      max =(void *) &dcmax;
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
      GA::SERVICES.error ((char *)"test_median:wrong data type.", type);
    }
  
  g_e->selectElem ((char *)"max", max, index);
  g_e->selectElem ((char *)"min", min, index);

  switch (type)
    {
      double r, m;
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
      GA::SERVICES.error ((char *)"test_median:wrong data type.", type);
    }


}


void
test_norm_infinity (GA::GlobalArray * g_a) {

  void *val;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;

  double norm_infinity = -1.0, result = -1.0;

  int me = GA_Nodeid ();
  int type, ndim, dims[MAXDIM];

  g_a->inquire (&type, &ndim, dims);
  dcval.real = -2.0;
  dcval.imag = -0.0;

  switch (type)
    {
    case C_INT:
      val = (void *)&ival;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      break;
    case C_DBL:
      val = (void *)&dval;
      break;
    case C_FLOAT:
      val =(void *) &fval;
      break;
    case C_LONG:
      val = (void *)&lval;
      break;
    default:
      GA::SERVICES.error ((char *)"test_norm_infinity:wrong data type.", type);
    }

  if (me == 0)
    printf ("Testing GA_Norm_infinity...");
  g_a->fill (val);
  g_a->normInfinity (&norm_infinity);

  // GA_Print(g_a);
  //printf("norm_infinity = %lf\n",norm_infinity);
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
    default:
      GA::SERVICES.error ((char *)"test_norm_infinity: wrong data type.\n", type);
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
test_norm1 (GA::GlobalArray * g_a)
{

  void *val;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;

  double norm1 = 0.0, result = -1.0;

  int me = GA_Nodeid ();
  int type, ndim, dims[MAXDIM];

  g_a->inquire (&type, &ndim, dims);
  dcval.real = -2.0;
  dcval.imag = -0.0;

  switch (type)
    {
    case C_INT:
      val = (void *)&ival;
      break;
    case C_DCPL:
      val =(void *) &dcval;
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
      GA::SERVICES.error ((char *)"test_norm1:wrong data type.", type);
    }

  if (me == 0)
    printf ("Testing GA_Norm1...");
  g_a->fill (val);
  g_a->norm1 (&norm1);
  // GA_Print(g_a);
  //printf("norm1=%lf\n", norm1);
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
    default:
      GA::SERVICES.error ((char *)"test_norm1: wrong data type.\n", type);
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
test_get_diagonal (GA::GlobalArray * g_a, 
		   GA::GlobalArray * g_v) {
  
  int me = GA_Nodeid ();
  void *val;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;

  int idot, iresult, ldot, lresult;
  double fdot, ddot, fresult, dresult;
  DoubleComplex zdot, zresult;

  int type, ndim, dims[MAXDIM];
  int vtype, vndim, vdims[MAXDIM];

  g_a->inquire (&type, &ndim, dims);
  g_v->inquire (&vtype, &vndim, vdims);
  dcval.real = -2.0;
  dcval.imag = -0.0;

  switch (type)
    {
    case C_INT:
      val = (void *)&ival;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      break;
    case C_DBL:
      val = (void *)&dval;
      break;
    case C_FLOAT:
      val =(void *) &fval;
      break;
    case C_LONG:
      val =(void *) &lval;
      break;
    default:
      GA::SERVICES.error ((char *)"test_get_diagonal:wrong data type.", type);
    }

  if (me == 0)
    printf ("Testing GA_Get_diagonal...");
  g_v->zero ();
  g_a->fill (val);
  g_v->getDiagonal (g_a);
  switch (type)
    {
    case C_INT:
      idot = vdims[0] * ival * ival;
      iresult = g_v->idot (g_v);
      if (me == 0)
	{
	  if (MISMATCHED (idot, iresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_LONG:
      ldot = ((long) vdims[0]) * lval * lval;
      lresult = g_v->ldot (g_v);
      if (me == 0)
	{
	  if (MISMATCHED (ldot, lresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_FLOAT:
      fdot = ((float) vdims[0]) * fval * fval;
      fresult = g_v->fdot (g_v);
      if (me == 0)
	{
	  if (MISMATCHED (fdot, fresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_DBL:
      ddot = ((double) vdims[0]) * dval * dval;
      dresult = g_v->ddot (g_v);
      if (me == 0)
	{
	  if (MISMATCHED (ddot, dresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_DCPL:
      zdot.real =
	((double) vdims[0]) * (dcval.real * dcval.real -
			       dcval.imag * dcval.imag);
      zdot.imag = ((double) vdims[0]) * (2.0 * dcval.real * dcval.imag);
      zresult = g_v->zdot (g_v);
      if (me == 0)
	{
	  if (MISMATCHED (zdot.real, zresult.real)
	      || MISMATCHED (zdot.imag, zresult.imag))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    default:
      GA::SERVICES.error ((char *)"test_get_diagonal:wrong data type:", type);
    }




}


void
test_add_diagonal (GA::GlobalArray * g_a, 
		   GA::GlobalArray * g_v)
{

  int me = GA_Nodeid ();
  void *val;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;

  int idot, iresult, ldot, lresult;
  double fdot, ddot, fresult, dresult;
  DoubleComplex zdot, zresult;

  int type, ndim, dims[MAXDIM];
  int vtype, vndim, vdims[MAXDIM];


  g_a->inquire (&type, &ndim, dims);
  dcval.real = -2.0;
  dcval.imag = -0.0;

  g_v->inquire (&vtype, &vndim, vdims);

  switch (type)
    {
    case C_INT:
      val =(void *) &ival;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      break;
    case C_DBL:
      val =(void *) &dval;
      break;
    case C_FLOAT:
      val =(void *) &fval;
      break;
    case C_LONG:
      val = (void *)&lval;
      break;
    default:
      GA::SERVICES.error ((char *)"test_add_diagonal:wrong data type.", type);
    }


  if (me == 0)
    printf ("Testing GA_Add_diagonal...");
  g_a->zero ();
  g_v->fill (val);
  g_a->setDiagonal (g_v);

  /*reassign value to val */
  ival = 3;
  dval = 3.0;
  fval = 3.0;
  lval = 3;
  dcval.real = 3.0;
  dcval.imag = -0.0;

  /*refile the global array g_v */
  g_v->fill (val);

  /*Add g_v to the diagonal of g_a */
  g_a->addDiagonal (g_v);
  /*after this line, the g_a should only have 1 on the diagonal and zeros every where else */

  switch (type)
    {
    case C_INT:
      idot = vdims[0];
      iresult = g_a->idot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (idot, iresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_LONG:
      ldot = ((long) vdims[0]);
      lresult = g_a->ldot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (ldot, lresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_FLOAT:
      fdot = ((float) vdims[0]);
      fresult = g_a->fdot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (fdot, fresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_DBL:
      ddot = (double) vdims[0];
      dresult = g_a->ddot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (ddot, dresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_DCPL:
      zdot.real = ((double) vdims[0]);
      zdot.imag = 0.0;
      zresult = g_a->zdot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (zdot.real, zresult.real)
	      || MISMATCHED (zdot.imag, zresult.imag))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    default:
      GA::SERVICES.error ((char *)"test_add_diagonal:wrong data type:", type);
    }

}

void 
test_set_diagonal (GA::GlobalArray * g_a, 
		   GA::GlobalArray * g_v)  {
  

  int me = GA_Nodeid ();
  void *val;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;

  int idot, iresult, ldot, lresult;
  double fdot, ddot, fresult, dresult;
  DoubleComplex zdot, zresult;

  int type, ndim, dims[MAXDIM];
  int vtype, vndim, vdims[MAXDIM];

  g_a->inquire (&type, &ndim, dims);
  g_v->inquire (&vtype, &vndim, vdims);
  dcval.real = -2.0;
  dcval.imag = -0.0;

  switch (type)
    {
    case C_INT:
      val = (void *)&ival;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      break;
    case C_DBL:
      val = (void *)&dval;
      break;
    case C_FLOAT:
      val = (void *)&fval;
      break;
    case C_LONG:
      val =(void *) &lval;
      break;
    default:
      GA::SERVICES.error ((char *)"test_set_diagonal:wrong data type.", type);
    }


  if (me == 0)
    printf ("Testing GA_Set_diagonal...");
  g_a->zero ();
  g_v->fill (val);
  g_a->setDiagonal (g_v);

  switch (type)
    {
    case C_INT:
      idot = vdims[0] * ival * ival;
      iresult = g_a->idot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (idot, iresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_LONG:
      ldot = ((long) vdims[0]) * lval * lval;
      lresult = g_a->ldot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (ldot, lresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_FLOAT:
      fdot = ((float) vdims[0]) * fval * fval;
      fresult = g_a->fdot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (fdot, fresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_DBL:
      ddot = ((double) vdims[0]) * dval * dval;
      dresult = g_a->ddot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (ddot, dresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_DCPL:
      zdot.real =
	((double) vdims[0]) * (dcval.real * dcval.real -
			       dcval.imag * dcval.imag);
      zdot.imag = ((double) dims[0]) * (2.0 * dcval.real * dcval.imag);
      zresult = g_a->zdot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (zdot.real, zresult.real)
	      || MISMATCHED (zdot.imag, zresult.imag))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    default:
      GA::SERVICES.error ((char *)"test_set_diagonal:wrong data type:", type);
    }

}

void
test_shift_diagonal (GA::GlobalArray *g_a) {
  
  int me = GA_Nodeid ();
  void *val;
  int ival = -2;
  double dval = -2.0;
  float fval = -2.0;
  long lval = -2;
  DoubleComplex dcval;

  int idot, iresult, ldot, lresult;
  double fdot, ddot, fresult, dresult;
  DoubleComplex zdot, zresult;
  int type, ndim, dims[MAXDIM];
  int dim;			/*the length of the diagonal */


  g_a->inquire (&type, &ndim, dims);

  dim = GA_MIN (dims[0], dims[1]);

  dcval.real = -2.0;
  dcval.imag = -0.0;
  switch (type)
    {
    case C_INT:
      val = (void *)&ival;
      break;
    case C_DCPL:
      val = (void *)&dcval;
      break;
    case C_DBL:
      val = (void *)&dval;
      break;
    case C_FLOAT:
      val =(void *) &fval;
      break;
    case C_LONG:
      val =(void *) &lval;
      break;
    default:
      GA::SERVICES.error ((char *)"test_shift_diagonal:wrong data type.", 
			  type);
    }


  if (me == 0)
    printf ("Testing GA_Shift_diagonal...");
  g_a->zero ();
  g_a->shiftDiagonal (val);

  switch (type)
    {
    case C_INT:
      idot = dim * ival * ival;
      iresult = g_a->idot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (idot, iresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_LONG:
      ldot = ((long) dim) * lval * lval;
      lresult = g_a->ldot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (ldot, lresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_FLOAT:
      fdot = ((float) dim) * fval * fval;
      fresult = g_a->fdot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (fdot, fresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_DBL:
      ddot = ((double) dim) * dval * dval;
      dresult = g_a->ddot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (ddot, dresult))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    case C_DCPL:
      zdot.real =
	((double) dim) * (dcval.real * dcval.real - dcval.imag * dcval.imag);
      zdot.imag = ((double) dim) * (2.0 * dcval.real * dcval.imag);
      zresult = g_a->zdot (g_a);
      if (me == 0)
	{
	  if (MISMATCHED (zdot.real, zresult.real)
	      || MISMATCHED (zdot.imag, zresult.imag))
	    printf ("not ok.\n");
	  else
	    printf ("ok.\n");
	}
      break;
    default:
      GA::SERVICES.error ((char *)"test_shift_diagonal:wrong data type: ", 
			  type);
    }

}



void
do_work (int type, int op) {
  GA::GlobalArray *g_a, *g_b, *g_c, *g_m, *g_v;
  int n = N;
  int dims[2] = { N,		/*N columns */
		  N + 2		/*N+2 rows */
  };
  int vdim;
  int lo[2], hi[2];

  lo[0] = 1;
  hi[0] = dims[0] - 1;
  lo[1] = 1;
  hi[1] = dims[1] - 1;

  switch (op)
    {

    case OP_SHIFT_DIAGONAL:
      g_a = GA::SERVICES.createGA (type, 2, dims, (char *)"A", NULL);
      test_shift_diagonal (g_a);
      g_a->destroy();
      break;
    case OP_SET_DIAGONAL:
      g_a = GA::SERVICES.createGA (type, 2, dims, (char *)"A", NULL);
      /*find out the diagonal length of the matrix A */
      vdim = GA_MIN (dims[0], dims[1]);
      g_v = GA::SERVICES.createGA (type, 1, &vdim, (char *)"V", NULL);
      test_set_diagonal (g_a, g_v);
      g_a->destroy ();
      g_v->destroy ();
      break;
    case OP_ADD_DIAGONAL:
      g_a = GA::SERVICES.createGA (type, 2, dims, (char *)"A", NULL);
      /*find out the diagonal length of the matrix A */
      vdim = GA_MIN (dims[0], dims[1]);
      g_v = GA::SERVICES.createGA (type, 1, &vdim, (char *)"V", NULL);
      test_add_diagonal (g_a, g_v);
      g_a->destroy ();
      g_v->destroy ();
      break;
    case OP_GET_DIAGONAL:
      g_a = GA::SERVICES.createGA (type, 2, dims, (char *)"A", NULL);
      /*find out the diagonal length of the matrix A */
      vdim = GA_MIN (dims[0], dims[1]);
      g_v = GA::SERVICES.createGA (type, 1, &vdim, (char *)"V", NULL);
      test_get_diagonal (g_a, g_v);
      g_a->destroy ();
      g_v->destroy ();
      break;
    case OP_NORM1:
      g_a = GA::SERVICES.createGA (type, 2, dims, (char *)"A", NULL);
      if (!g_a)
	GA_Error ((char *)"create failed: A", n);
      test_norm1 (g_a);
      g_a->destroy ();
      break;

    case OP_NORM_INFINITY:
      g_a = GA::SERVICES.createGA (type, 2, dims,(char *) "A", NULL);
      test_norm_infinity (g_a);
      g_a->destroy ();
      break;

    case OP_MEDIAN:
      g_a = GA::SERVICES.createGA (type, 2, dims, (char *)"A", NULL);
      /*duplicate g_a */
      g_b = GA::SERVICES.createGA (g_a, (char *)"B");      
      g_c = GA::SERVICES.createGA (g_a, (char *)"C");
#if 0 //test g_m is different from g_a, g_b, amd g_c
      g_m = GA::SERVICES.createGA (g_a, (char *)"M");
      test_median (g_a, g_b, g_c, g_m);
#else //test g_m = g_c
      test_median (g_a, g_b, g_c, g_a);
#endif
      g_a->destroy ();
      g_b->destroy ();
      g_c->destroy ();
#if 0 //test g_m is different from g_a, g_b, g_c
      g_m->destroy ();
#endif
      break;

    case OP_MEDIAN_PATCH:
      g_a = GA::SERVICES.createGA (type, 2, dims, (char *)"A", NULL);
      /*duplicate g_a */
      g_b = GA::SERVICES.createGA (g_a, (char *)"B");
      g_c = GA::SERVICES.createGA (g_a, (char *)"C");
      g_m = GA::SERVICES.createGA (g_a, (char *)"M");
      test_median_patch (g_a, lo, hi, g_b, lo, hi, g_c, lo, hi, g_m, lo, hi);
      g_a->destroy ();
      g_b->destroy ();
      g_c->destroy ();
      g_m->destroy ();
      break;
    case OP_SCALE_ROWS:
      g_a = GA::SERVICES.createGA (type, 2, dims, (char *)"A", NULL);
      /*find out the diagonal length of the matrix A */
      vdim = dims[1];
      g_v = GA::SERVICES.createGA (type, 1, &vdim, (char *)"V", NULL);
      test_scale_rows (g_a, g_v);
      g_a->destroy ();
      g_v->destroy ();
      break;
    case OP_SCALE_COLS:
      g_a = GA::SERVICES.createGA (type, 2, dims, (char *)"A", NULL);
      /*find out the diagonal length of the matrix A */
      vdim = dims[0];
      g_v = GA::SERVICES.createGA (type, 1, &vdim, (char *)"V", NULL);
      test_scale_cols (g_a, g_v);
      g_a->destroy ();
      g_v->destroy ();
      break;

    default:
      GA::SERVICES.error ((char *)"test_function: wrong operation.", op);
    }
}




int 
main(int argc, char *argv[]) {
 
  int me, nproc;
  int heap  = 200000, stack = 200000;
  int op;

  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, 0);
  me=GA_Nodeid(); nproc=GA_Nnodes();
  if(!me) cout << "Using " << nproc << " processes\n";
  for (op = 1; op < 11; op++) {
    if(me == 0) printf ("\n\n");
    if (me == 0)  printf ("type = C_INT \t ");
    do_work (C_INT, op);
    
    if (me == 0)  printf ("type = C_LONG \t ");
    do_work (C_LONG, op);
    
    if (me == 0)  printf ("type = C_FLOAT \t ");
    do_work (C_FLOAT, op);
    
    if (me == 0)  printf ("type = C_DBL \t ");
    do_work (C_DBL, op);
    
    if (me == 0)  printf ("type = C_DCPL \t ");
    do_work (C_DCPL, op);
  }
  
  if(!me) cout << "Terminating\n";
  GA::Terminate();
}
