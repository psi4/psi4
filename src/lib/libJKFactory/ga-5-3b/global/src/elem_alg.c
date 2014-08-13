#if HAVE_CONFIG_H
#   include "config.h"
#endif

/**************************************************************
File: elem.alg.c

Elementwise operations on patches and whole arrays

Author: Limin Zhang, Ph.D.
	Mathematics Department
        Columbia Basin College
        Pasco, WA 99301
        Limin.Zhang@cbc2.org

Mentor: Jarek Nieplocha.
  	Environmental Molecular Science Laboratory
        Pacific Northwest National Laboratory
	Richland, WA 99352

Date: 1/18/2002
    
Purpose:
      to design and implement some interfaces between TAO and
      global arrays.

Modified 3/2004 By Doug Baxter to increase robustness.

**************************************************************/

#include "globalp.h"
#if HAVE_MATH_H
#   include <math.h>
#endif
#include "abstract_ops.h"
#include "ga-papi.h"
#include "ga-wapi.h"

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
/*
  Original value below.
#define GA_NEGATIVE_INFINITY_D -1.0e20
*/
/*
  End of 01/24/04 Modification.
  Perhaps it would be more appropriate to have GA_INFINITY_D BE 1.0e307
  These ranges make assumptions about the data.
*/
#define OP_ABS 0
#define OP_ADD_CONST 1
#define OP_RECIP 2
#define OP_ELEM_MULT 3
#define OP_ELEM_DIV 4
#define OP_ELEM_MAX 5
#define OP_ELEM_MIN 6
#define OP_STEPMAX 7
#define OP_STEPBOUNDINFO 8
#define OP_ELEM_SDIV 9
#define OP_ELEM_SDIV2 10
#define OP_STEP_MASK 11
#define OP_FILL 100 /*The OP_FILL is not currently in use */

int debug_gai_oper_elem = 1;

static void do_stepboundinfo(void *ptr, int nelem, int type)
/*look at elements one by one and replace the positive infinity with negative infinity */ 
{
    int i;
    switch (type){
         int *ia;
         double *da;
         float *fa;
         long *la;

  	 case C_DBL:
                /*Only double data type will be handled for TAO/GA project*/ 
              da = (double *) ptr;
              for(i=0;i<nelem;i++)
		/* Modified 01/24/04 to add _D ending to constants. */
		if(da[i]>=GA_INFINITY_D) da[i]=-GA_INFINITY_D;
              break;
         case C_INT:
	   /* This block added 01/24/04 */
	   ia = (int *) ptr;
	   for (i=0;i<nelem;i++)
	     if (ia[i]>= GA_INFINITY_I) ia[i] = GA_NEGATIVE_INFINITY_I;
	   break;
         case C_DCPL: 
         case C_SCPL: 
	   /* This operation is not well defined for complex
	      numbers . This statement added when drop through
	      behavior changed by adding code for C_FLOAT and C_LONG
	      cases below. 01/24/04
	   */
	   pnga_error("do_stepboundinfo:wrong data type",type);
         case C_FLOAT:
	   /* This case added 01/24/04 */
	   fa = (float *) ptr;
	   for (i=0;i<nelem;i++)
	     if (fa[i]>= GA_INFINITY_F) fa[i] = GA_NEGATIVE_INFINITY_F;
	   break;
 	case C_LONG:
	   /* This case added 01/24/04 */
	   la = (long *) ptr;
	   for (i=0;i<nelem;i++)
	     if (la[i]>= GA_INFINITY_L) la[i] = GA_NEGATIVE_INFINITY_L;
	   break;
         default: pnga_error("do_stepboundinfo:wrong data type",type);
    }
}
static void do_stepmax(void *ptr, int nelem, int type)
/*
  Look at elements one by one and replace the positive with negative infinity.
*/ 
{
    int i;
    switch (type){
         int *ia,i_0;
         double *da,d_0;
         float *fa,f_0;
	 long *la,l_0;


  	 case C_DBL:
                /*Only double data type will be handled for TAO/GA project*/ 
              da = (double *) ptr;
	      d_0 = (double) 0.0;
	      /* Modified 01/24/04 to use _D ending. */
              for(i=0;i<nelem;i++) 
		if(da[i]>d_0) da[i]=-GA_INFINITY_D;
              break;
         case C_INT:
	   /* Thix case added 01/24/04*/
              ia = (int *) ptr;
	      i_0 = (int)0;
              for(i=0;i<nelem;i++) 
		if(ia[i]>i_0)ia[i]=-GA_INFINITY_I;
              break;
         case C_DCPL:
         case C_SCPL:
	   /* This operation is not well defined for complex
	      numbers . This statement added when drop through
	      behavior changed by adding code for C_FLOAT and C_LONG
	      cases below. 01/24/04
	   */
	   pnga_error("do_stepmax:wrong data type",type);
         case C_FLOAT:
	   /* Thix case added 01/24/04*/
              fa = (float *) ptr;
	      f_0 = (float) 0.0;
              for(i=0;i<nelem;i++) 
		if(fa[i]>f_0) fa[i]=-GA_INFINITY_F;
              break;
         case C_LONG:
	   /* Thix case added 01/24/04*/
              la = (long *) ptr;
	      l_0 = (long)0;
              for(i=0;i<nelem;i++)
		if(la[i]>l_0) la[i]=-GA_INFINITY_L;
              break;
         default: pnga_error("do_stepmax:wrong data type",type);
    }
}




static void do_abs(void *ptr, int nelem, int type)
{
    int i;
    double x2;
    float sx2;
    switch (type){
         int *ia;
         double *da;
         float *fa;
         DoubleComplex *ca,val;
         SingleComplex *cfa,cval;
         long *la;

         case C_INT:
              ia = (int *)ptr; 
              for(i=0;i<nelem;i++)
                  ia[i]= GA_ABS(ia[i]);
              break; 
         case C_DCPL:
              ca = (DoubleComplex *) ptr;
              for(i=0;i<nelem;i++){
#if HAVE_HYPOT
                  ca[i].real = hypot(ca[i].real, ca[i].imag);
#else
                  val = ca[i];
                  /* DJB: This algorithm can lead to overflows when
                     and underflows when not necessary.
                  ca[i].real = sqrt(val.real * val.real + val.imag *val.imag);
                  ca[i].imag = 0.0;
                     Better (but slower) is: */
                  if (GA_ABS(val.real) >= GA_ABS(val.imag)) {
                      if (val.real == (double)0.0) {
                          ca[i].real = (double)0.0;
                      } else {
                          x2 = val.imag/val.real;
                          ca[i].real = GA_ABS(val.real)*sqrt(1.0+(x2*x2));
                      }
                  } else {
                      x2 = val.real/val.imag;
                      ca[i].real = GA_ABS(val.imag)*sqrt(1.0+(x2*x2));
                  }
#endif
                  ca[i].imag=(double)0.0;
              }
              break;
         case C_SCPL:
              cfa = (SingleComplex *) ptr;
              for(i=0;i<nelem;i++){
#if HAVE_HYPOT
                  cfa[i].real = hypotf(cfa[i].real, cfa[i].imag);
#else
                  cval = cfa[i];
                  /* DJB: This algorithm can lead to overflows when
                     and underflows when not necessary.
                  cfa[i].real = sqrt(cval.real * cval.real + cval.imag *cval.imag);
                  cfa[i].imag = 0.0;
                     Better (but slower) is: */
                  if (GA_ABS(cval.real) >= GA_ABS(cval.imag)) {
                      if (cval.real == 0.0f) {
                          cfa[i].real = 0.0f;
                      } else {
                          sx2 = cval.imag/cval.real;
                          cfa[i].real = GA_ABS(cval.real)*sqrt(1.0f+(sx2*sx2));
                      }
                  } else {
                      sx2 = cval.real/cval.imag;
                      cfa[i].real = GA_ABS(cval.imag)*sqrt(1.0f+(sx2*sx2));
                  }
#endif
                  cfa[i].imag=0.0f;
              }
              break;
  	 case C_DBL:
              da = (double *) ptr;
              for(i=0;i<nelem;i++)
                  da[i]= GA_ABS(da[i]);
              break;
         case C_FLOAT:
              fa = (float *)ptr;
              for(i=0;i<nelem;i++)
                  fa[i]= GA_ABS(fa[i]);
              break;
 	case C_LONG:
              la = (long *)ptr;
              for(i=0;i<nelem;i++)
                  la[i]= GA_ABS(la[i]);
              break;

         default: pnga_error("wrong data type",type);
    }
} 

static void do_recip(void *ptr, int nelem, int type)
{
  /*
    DJB general comment, as I found this routine, it
    would return some form of infinity when having 
    a zero denominator. I find that technically incorrect
    and the commented out error message to be preferrable.
    If returning infinity is the behavior desired I would recommend
    having a separate specialized reciprocal, where it's
    in the GA documentation, what should happen on a divide
    by zero. In IEEE standard, the default is to throw
    a floating point exception (FPE) when division by zero
    occurs. The ga_error message seems to be the
    moral equivalent of that.
    I have commented out the Infinity returns and returned 
    to the ga_error calls. Also on the commented out 
    infinity value returns I have been more specific in
    the INFINITY type returned (the trailing _*).
  */
  double magi, magr, x1, x2, c, d;
  float smagi, smagr, sx1, sx2, sc, sd;
    int i;
    switch (type){
         int *ia;
         double *da; /*, temp; */
         float *fa;
         DoubleComplex *ca;
         SingleComplex *cfa;
         long *la; 

         case C_INT:
              ia = (int *)ptr;
              for(i=0;i<nelem;i++)
                  if(ia[i]!=0) ia[i]= 1/ia[i];
                     else
                   pnga_error("zero value at index",i);
		       /*
			 ia[i] = GA_INFINITY_I;
		       */
              break;
         case C_DCPL:
              ca = (DoubleComplex *) ptr;
              for(i=0;i<nelem;i++){
		/* 
		  Again, as for absolute value the following
		  algorithm can lead to unecessary overflow/underflow
		   
                  temp = ca[i].real*ca[i].real + ca[i].imag*ca[i].imag;
                  if( temp!=0.0){
                   ca[i].real =ca[i].real/temp;
                   ca[i].imag =-ca[i].imag/temp;
                  }
                  else{
 		     pnga_error("zero value at index",i); 
                     OR
		       ca[i].real = GA_INFINITY_D;
		       ca[i].imag = GA_INFINITY_D;

                 }
		 Better (but slower) is: 
	       */		     
		x1 = ca[i].real;
		x2 = ca[i].imag;
		/*
		printf(" do_recip i = %d, x1 = %le, x2 = %le\n",
		       i,x1,x2);
		*/
		magr = GA_ABS(x1);
		magi = GA_ABS(x2);
		/*
		printf(" do_recip i = %d, magr = %le, magi = %le\n",
		       i,magr,magi);
		*/
		if (magr >= magi) {
		  if (magr != ((double)0.0)) {
		    c = x2/x1;
		    d = ((double)1.0)/((((double)1.0) + (c*c))*x1);
		    ca[i].real = d;
		    ca[i].imag = -c*d;
		  } else {
		    pnga_error("zero value at index",i); 
		  }
		} else {
		  c = x1/x2;
		  d = ((double)1.0)/((((double)1.0) + (c*c))*x2);
		  ca[i].real = c*d;
		  ca[i].imag = -d;
		}
                /*
		printf(" do_recip ca[%d].real = %le, ca[%d].imag = %le\n",
		       i,ca[i].real,i,ca[i].imag);
		*/

	      }
              break;
         case C_SCPL:
              cfa = (SingleComplex *) ptr;
              for(i=0;i<nelem;i++){
		/* 
		  Again, as for absolute value the following
		  algorithm can lead to unecessary overflow/underflow
		   
                  temp = cfa[i].real*cfa[i].real + cfa[i].imag*cfa[i].imag;
                  if( temp!=0.0){
                   cfa[i].real =cfa[i].real/temp;
                   cfa[i].imag =-cfa[i].imag/temp;
                  }
                  else{
 		     pnga_error("zero value at index",i); 
                     OR
		       cfa[i].real = GA_INFINITY_D;
		       cfa[i].imag = GA_INFINITY_D;

                 }
		 Better (but slower) is: 
	       */		     
		sx1 = cfa[i].real;
		sx2 = cfa[i].imag;
		/*
		printf(" do_recip i = %d, x1 = %le, x2 = %le\n",
		       i,x1,x2);
		*/
		smagr = GA_ABS(sx1);
		smagi = GA_ABS(sx2);
		/*
		printf(" do_recip i = %d, magr = %le, magi = %le\n",
		       i,magr,magi);
		*/
		if (smagr >= smagi) {
		  if (smagr != ((float)0.0)) {
		    sc = sx2/sx1;
		    sd = ((float)1.0)/((((float)1.0) + (sc*sc))*sx1);
		    cfa[i].real = sd;
		    cfa[i].imag = -sc*sd;
		  } else {
		    pnga_error("zero value at index",i); 
		  }
		} else {
		  sc = sx1/sx2;
		  sd = ((float)1.0)/((((float)1.0) + (sc*sc))*sx2);
		  cfa[i].real = sc*sd;
		  cfa[i].imag = -sd;
		}
                /*
		printf(" do_recip ca[%d].real = %le, ca[%d].imag = %le\n",
		       i,ca[i].real,i,ca[i].imag);
		*/

	      }
              break;
         case C_DBL:
              da = (double *) ptr;
              for(i=0;i<nelem;i++)
                  if(da[i]!=(double)0.0) da[i]= ((double)1.0)/da[i];
  		     else
		   pnga_error("zero value at index",i); 
		    /* 
		       da[i] = GA_INFINITY_D;
		    */
              break;
         case C_FLOAT:
              fa = (float *)ptr;
              for(i=0;i<nelem;i++)
                  if(fa[i]!=(float)0.0) fa[i]= ((float)1.0)/fa[i];
                     else
		   pnga_error("zero value at index",i); 
         	   /*
		     fa[i] = GA_INFINITY_F
		   */;
              break;
	case C_LONG:
              la = (long *)ptr;
              for(i=0;i<nelem;i++)
                  if(la[i]!=(long)0) la[i]= ((long)1)/la[i];
                     else
                  pnga_error("zero value at index",i); 
	          /*
		    la[i] = GA_INFINITY_I;
		  */
              break;


         default: pnga_error("wrong data type",type);
    }
} 

static void do_add_const(void *ptr, int nelem, int type, void *alpha)
{
    int i;
    switch (type){
         int *ia;
         double *da;
         float *fa;
         DoubleComplex *ca,val;
         SingleComplex *cfa,cval;
	 long *la;

         case C_INT:
              ia = (int *)ptr;
              for(i=0;i<nelem;i++)
                  ia[i] += *(int *)alpha;
              break;
         case C_DCPL:
              ca = (DoubleComplex *) ptr;
              for(i=0;i<nelem;i++){
                  val = *(DoubleComplex*)alpha;
                  ca[i].real += val.real;
                  ca[i].imag += val.imag;
              }
              break;
         case C_SCPL:
              cfa = (SingleComplex *) ptr;
              for(i=0;i<nelem;i++){
                  cval = *(SingleComplex*)alpha;
                  cfa[i].real += cval.real;
                  cfa[i].imag += cval.imag;
              }
              break;
         case C_DBL:
              da = (double *) ptr;
              for(i=0;i<nelem;i++)
                  da[i] += *(double*)alpha;
              break;
         case C_FLOAT:
              fa = (float *)ptr;
              for(i=0;i<nelem;i++)
                  fa[i] += *(float*)alpha;
              break;
	 case C_LONG:
              la = (long *)ptr;
              for(i=0;i<nelem;i++)
                  la[i] += *(long *)alpha;
              break;

         default: pnga_error("wrong data type",type);
    }
} 

/*
void do_fill(void *ptr, int nelem, int type, void *alpha)
{
    int i;
    switch (type){
         int *ia;
         double *da;
         float *fa;
         DoubleComplex *ca,val;
         long *la;

         case C_INT:
              ia = (int *)ptr;
              for(i=0;i<nelem;i++)
                  ia[i] = *(int *)alpha;
              break;
         case C_DCPL:
              ca = (DoubleComplex *) ptr;
              for(i=0;i<nelem;i++){
                  val = *(DoubleComplex*)alpha;
                  ca[i].real = val.real;
                  ca[i].imag = val.imag;
              }
              break;
         case C_SCPL:
              ca = (SingleComplex *) ptr;
              for(i=0;i<nelem;i++){
                  val = *(SingleComplex*)alpha;
                  ca[i].real = val.real;
                  ca[i].imag = val.imag;
              }
              break;
         case C_DBL:
              da = (double *) ptr;
              for(i=0;i<nelem;i++)
                  da[i] = *(double*)alpha;
              break;
         case C_FLOAT:
              fa = (float *)ptr;
              for(i=0;i<nelem;i++)
                  fa[i] = *(float*)alpha;
              break;
	 case C_LONG:
              la = (long *)ptr;
              for(i=0;i<nelem;i++)
                  la[i] = *(long *)alpha;
              break;

         default: pnga_error("wrong data type",type);
    }
} 
*/

/*
Input Parameters

int g_a -- the global array handle
int *lo, *hi--the integer arrays that define the patch of the global array
void *scalar -- the pointer that points to the data to pass. When it is NULL, no scalar will be passed.
int op -- the operations to handle. For example op can be  

OP_ABS for pointwise taking absolute function 
OP_ADD_CONSTANT 2 for pointwise adding the same constant
OP_RECIP for pointwise taking reciprocal
OP_FILL for pointwise filling value 

Output Parameters

None

*/

/*\ Internal utility function to do operation.
\*/
static
void ngai_do_oper_elem(Integer type, Integer ndim, Integer *loA, Integer *hiA,
                       Integer *ld, void *data_ptr, void *scalar, Integer op)
{
  Integer i, j;    
  Integer bvalue[MAXDIM], bunit[MAXDIM], baseld[MAXDIM];
  void *temp = NULL;
  Integer idx, n1dim;
  /* number of n-element of the first dimension */
  n1dim = 1; for(i=1; i<ndim; i++) n1dim *= (hiA[i] - loA[i] + 1);

  /* calculate the destination indices */
  bvalue[0] = 0; bvalue[1] = 0; bunit[0] = 1; bunit[1] = 1;
  /* baseld[0] = ld[0]
   * baseld[1] = ld[0] * ld[1]
   * baseld[2] = ld[0] * ld[1] * ld[2] ....
   */
  baseld[0] = ld[0]; baseld[1] = baseld[0] *ld[1];
  for(i=2; i<ndim; i++) {
    bvalue[i] = 0;
    bunit[i] = bunit[i-1] * (hiA[i-1] - loA[i-1] + 1);
    baseld[i] = baseld[i-1] * ld[i];
  }

  for(i=0; i<n1dim; i++) {
    idx = 0;
    for(j=1; j<ndim; j++) {
      idx += bvalue[j] * baseld[j-1];
      if(((i+1) % bunit[j]) == 0) bvalue[j]++;
      if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
    }

    switch(type){
      case C_INT: 
        temp=((int*)data_ptr)+idx; 
        break;
      case C_DCPL: 
        temp=((DoubleComplex*)data_ptr)+idx; 
        break;
      case C_SCPL: 
        temp=((SingleComplex*)data_ptr)+idx; 
        break;
      case C_DBL: 
        temp=((double*)data_ptr)+idx; 
        break;
      case C_FLOAT:
        temp=((float*)data_ptr)+idx;
        break;
      case C_LONG:
        temp=((long *)data_ptr)+idx;
        break;
      default: pnga_error("wrong data type.",type);	

    }

    switch(op){
      case OP_ABS:
        do_abs(temp ,hiA[0] -loA[0] +1, type); break;
        break;
      case OP_ADD_CONST:
        do_add_const(temp ,hiA[0] -loA[0] +1, type, scalar); 
        break;
      case OP_RECIP:
        do_recip(temp ,hiA[0] -loA[0] +1, type); break;
        break;
      default: pnga_error("bad operation",op);
    }
  }

}

static void gai_oper_elem(Integer g_a, Integer *lo, Integer *hi, void *scalar, Integer op)
{

  Integer ndim, dims[MAXDIM], type;
  Integer loA[MAXDIM], hiA[MAXDIM], ld[MAXDIM];
  void /* *temp,*/ *data_ptr;
  Integer me= pnga_nodeid();
  Integer num_blocks;
  int local_sync_begin,local_sync_end;

  local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
  _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
  if(local_sync_begin)pnga_sync();

  pnga_check_handle(g_a, "gai_oper_elem");

  GA_PUSH_NAME("gai_oper_elem");

  pnga_inquire(g_a,  &type, &ndim, dims);
  num_blocks = pnga_total_blocks(g_a);

  if (num_blocks < 0) {
    /* get limits of VISIBLE patch */
    pnga_distribution(g_a, me, loA, hiA);

    /*  determine subset of my local patch to access  */
    /*  Output is in loA and hiA */
    if(pnga_patch_intersect(lo, hi, loA, hiA, ndim)){

      /* get data_ptr to corner of patch */
      /* ld are leading dimensions INCLUDING ghost cells */
      pnga_access_ptr(g_a, loA, hiA, &data_ptr, ld);

      /* perform operation on all elements in local patch */
      ngai_do_oper_elem(type, ndim, loA, hiA, ld, data_ptr, scalar, op);

      /* release access to the data */
      pnga_release_update(g_a, loA, hiA);
    }
  } else {
    Integer offset, i, j, jtmp, chk;
    Integer loS[MAXDIM];
    Integer nproc = pnga_nnodes();
    /* using simple block-cyclic data distribution */
    if (!pnga_uses_proc_grid(g_a)){
      for (i=me; i<num_blocks; i += nproc) {

        /* get limits of patch */
        pnga_distribution(g_a, i, loA, hiA); 

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
          pnga_access_block_ptr(g_a, i, &data_ptr, ld);

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
              default: pnga_error(" wrong data type ",type);
            }
          }
          /* perform operation on all elements in local patch */
          ngai_do_oper_elem(type, ndim, loA, hiA, ld, data_ptr, scalar, op);

          /* release access to the data */
          pnga_release_update_block(g_a, i);
        }
      }
    } else {
      /* using scalapack block-cyclic data distribution */
      Integer proc_index[MAXDIM], index[MAXDIM];
      Integer topology[MAXDIM];
      Integer blocks[MAXDIM], block_dims[MAXDIM];
      pnga_get_proc_index(g_a, me, proc_index);
      pnga_get_proc_index(g_a, me, index);
      pnga_get_block_info(g_a, blocks, block_dims);
      pnga_get_proc_grid(g_a, topology);
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
          pnga_access_block_grid_ptr(g_a, index, &data_ptr, ld);

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
              default: pnga_error(" wrong data type ",type);
            }
          }

          /* perform operation on all elements in local patch */
          ngai_do_oper_elem(type, ndim, loA, hiA, ld, data_ptr, scalar, op);

          /* release access to the data */
          pnga_release_update_block_grid(g_a, index);
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
  if(local_sync_end)pnga_sync();
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_abs_value_patch = pnga_abs_value_patch
#endif
void pnga_abs_value_patch(Integer g_a, Integer *lo, Integer *hi)
{
    gai_oper_elem(g_a, lo, hi, NULL, OP_ABS);
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_recip_patch = pnga_recip_patch
#endif
void pnga_recip_patch(Integer g_a, Integer *lo, Integer *hi)
{
    gai_oper_elem(g_a, lo, hi, NULL, OP_RECIP);

}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_add_constant_patch = pnga_add_constant_patch
#endif
void pnga_add_constant_patch(Integer g_a, Integer *lo, Integer *hi, void *alpha)
{
    gai_oper_elem(g_a, lo, hi, alpha, OP_ADD_CONST);

}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_abs_value = pnga_abs_value
#endif
void pnga_abs_value(Integer g_a)
{
   Integer type, ndim;
   Integer lo[MAXDIM],hi[MAXDIM];

    pnga_inquire(g_a,  &type, &ndim, hi);
    while(ndim){
        lo[ndim-1]=1;
        ndim--;
    }
    _ga_sync_begin = 1; /*just to be on the safe side*/
    gai_oper_elem(g_a, lo, hi, NULL, OP_ABS);
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_add_constant = pnga_add_constant
#endif
void pnga_add_constant(Integer g_a, void *alpha)
{
   Integer type, ndim;
   Integer lo[MAXDIM],hi[MAXDIM];

    pnga_inquire(g_a,  &type, &ndim, hi);
    while(ndim){
        lo[ndim-1]=1;
        ndim--;
    }
    _ga_sync_begin = 1; /*just to be on the safe side*/
    gai_oper_elem(g_a, lo, hi, alpha, OP_ADD_CONST);
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_recip = pnga_recip
#endif
void pnga_recip(Integer g_a)
{        
   Integer type, ndim;
   Integer lo[MAXDIM],hi[MAXDIM];
        
    pnga_inquire(g_a,  &type, &ndim, hi);
    while(ndim){
        lo[ndim-1]=1; 
        ndim--;
    }
    _ga_sync_begin = 1; /*just to be on the safe side*/
    gai_oper_elem(g_a, lo, hi, NULL, OP_RECIP);
}    



static void do_multiply(void *pA, void *pB, void *pC, Integer nelems, Integer type){
#if 1
    int i;
    switch (type) {
#define TYPE_CASE(MT,T,AT)                                                  \
        case MT:                                                            \
            {                                                               \
                const T* const restrict aptr = (const T* const restrict)pA; \
                const T* const restrict bptr = (const T* const restrict)pB; \
                T* const restrict cptr = (T* const restrict)pC;             \
                for (i=0; i<nelems; i++) {                                  \
                    assign_mul_##AT(cptr[i],aptr[i],bptr[i]);               \
                }                                                           \
                break;                                                      \
            }
#include "types.xh"
#undef TYPE_CASE
        default: pnga_error("do_multiply:wrong data type ",type);
    }
#else
  Integer i;
  
  switch(type){
    double aReal, aImag, bReal, bImag;
    
  case C_DBL:
    for(i = 0; i<nelems; i++)
      ((double*)pC)[i]= ((double*)pA)[i]*((double*)pB)[i]; 
    break;
  case C_DCPL:
    for(i = 0; i<nelems; i++) {
      aReal = ((DoubleComplex*)pA)[i].real; 
      bReal = ((DoubleComplex*)pB)[i].real; 
      aImag = ((DoubleComplex*)pA)[i].imag; 
      bImag = ((DoubleComplex*)pB)[i].imag; 
      ((DoubleComplex*)pC)[i].real = aReal*bReal-aImag*bImag;
      ((DoubleComplex*)pC)[i].imag = aReal*bImag+aImag*bReal;
    }
    break;
  case C_SCPL:
    for(i = 0; i<nelems; i++) {
      aReal = ((SingleComplex*)pA)[i].real; 
      bReal = ((SingleComplex*)pB)[i].real; 
      aImag = ((SingleComplex*)pA)[i].imag; 
      bImag = ((SingleComplex*)pB)[i].imag; 
      ((SingleComplex*)pC)[i].real = aReal*bReal-aImag*bImag;
      ((SingleComplex*)pC)[i].imag = aReal*bImag+aImag*bReal;
    }
    break;
  case C_INT:
    for(i = 0; i<nelems; i++)
      ((int*)pC)[i] = ((int*)pA)[i]* ((int*)pB)[i];
    break;
  case C_FLOAT:
    for(i = 0; i<nelems; i++)
      ((float*)pC)[i]=  ((float*)pA)[i]*((float*)pB)[i];
    break;
  case C_LONG:
    for(i = 0; i<nelems; i++)
      ((long *)pC)[i]= ((long *)pA)[i]* ((long *)pB)[i];
    break;
    
  default: pnga_error(" wrong data type ",type);
  }
#endif
}


static void do_divide(void *pA, void *pB, void *pC, Integer nelems, Integer type){
  /* 
    Modified by Doug Baxter to call ga_error on divide by 0.
    A do_step_divide is added below to return properly signed
    infinity for the step_max functions.
  */
  Integer i;
  double aReal, aImag, bReal, bImag;
  double x1,x2;

  switch(type){
  
  case C_DBL:
    for(i = 0; i<nelems; i++) {
      if(((double*)pB)[i]!=(double)0.0)
	((double*)pC)[i]=  ((double*)pA)[i]/((double*)pB)[i];
      else{
	pnga_error("zero divisor ",((double*)pB)[i]); 
      }
    }
    break;
  case C_DCPL:
    for(i = 0; i<nelems; i++) {
      aReal = ((DoubleComplex*)pA)[i].real;
      bReal = ((DoubleComplex*)pB)[i].real;
      aImag = ((DoubleComplex*)pA)[i].imag;
      bImag = ((DoubleComplex*)pB)[i].imag;
      /* 
        The following original algorithm overflows
	when it need not.
      temp = bReal*bReal+bImag*bImag;
      if(temp!=0.0){
	((DoubleComplex*)pC)[i].real
	  =(aReal*bReal+aImag*bImag)/temp;
	((DoubleComplex*)pC)[i].imag
	  =(aImag*bReal-aReal*bImag)/temp;
      }
      else{
	pnga_error("zero divisor ",temp); 
      }
      */
      if (GA_ABS(bReal) >= GA_ABS(bImag)) {
	if (bReal != (double)0.0) {
	  x1 = bImag/bReal;
          /* So x1 <= 1 */
	  x2 = ((double)1.0)/(bReal*(((double)1.0)+(x1*x1)));
	  ((DoubleComplex*)pC)[i].real = (aReal + aImag*x1)*x2;
	  ((DoubleComplex*)pC)[i].imag = (aImag - aReal*x1)*x2;
	}
	else{
	  pnga_error("zero divisor ",bReal); 
	}
      } else {
	x1 = bReal/bImag;
        /* So x1 <= 1 */
	x2 = ((double)1.0)/(bImag*(((double)1.0)+(x1*x1)));
	((DoubleComplex*)pC)[i].real = (aReal*x1 + aImag)*x2;
	((DoubleComplex*)pC)[i].imag = (aImag*x1 - aReal)*x2;
      }
    }
    break;
  case C_SCPL:
    for(i = 0; i<nelems; i++) {
      aReal = ((SingleComplex*)pA)[i].real;
      bReal = ((SingleComplex*)pB)[i].real;
      aImag = ((SingleComplex*)pA)[i].imag;
      bImag = ((SingleComplex*)pB)[i].imag;
      /* 
        The following original algorithm overflows
	when it need not.
      temp = bReal*bReal+bImag*bImag;
      if(temp!=0.0){
	((SingleComplex*)pC)[i].real
	  =(aReal*bReal+aImag*bImag)/temp;
	((SingleComplex*)pC)[i].imag
	  =(aImag*bReal-aReal*bImag)/temp;
      }
      else{
	pnga_error("zero divisor ",temp); 
      }
      */
      if (GA_ABS(bReal) >= GA_ABS(bImag)) {
	if (bReal != (float)0.0) {
	  x1 = bImag/bReal;
          /* So x1 <= 1 */
	  x2 = ((float)1.0)/(bReal*(((float)1.0)+(x1*x1)));
	  ((SingleComplex*)pC)[i].real = (aReal + aImag*x1)*x2;
	  ((SingleComplex*)pC)[i].imag = (aImag - aReal*x1)*x2;
	}
	else{
	  pnga_error("zero divisor ",bReal); 
	}
      } else {
	x1 = bReal/bImag;
        /* So x1 <= 1 */
	x2 = ((float)1.0)/(bImag*(((float)1.0)+(x1*x1)));
	((SingleComplex*)pC)[i].real = (aReal*x1 + aImag)*x2;
	((SingleComplex*)pC)[i].imag = (aImag*x1 - aReal)*x2;
      }
    }
    break;
  case C_INT:
    for(i = 0; i<nelems; i++){
      if(((int*)pB)[i]!=0)
	((int*)pC)[i] = ((int*)pA)[i]/((int*)pB)[i];
      else{
	pnga_error("zero divisor ",((int*)pB)[i]); 
      } 
    }
    break;
  case C_FLOAT:
    for(i = 0; i<nelems; i++){
      if(((float*)pB)[i]!=(float)0.0) 
	((float*)pC)[i]=  ((float*)pA)[i]/((float*)pB)[i];
      else{
	pnga_error("zero divisor ",((float*)pB)[i]); 
      }
    }
    break;
  case C_LONG:
    for(i = 0; i<nelems; i++){
      if(((long *)pB)[i]!=0)
	((long *)pC)[i]=  ((long *)pA)[i]/((long *)pB)[i];
      else{
	pnga_error("zero divisor ",((long*)pB)[i]); 
      }
    }
    break;		
  default: pnga_error(" wrong data type ",type);
  }
}
 

static void do_step_divide(void *pA, void *pB, void *pC, Integer nelems, Integer type){
  /* Elementwise divide, not aborting on a zero denominator, but
     returning an infinity. If an element in the numerator vector (PA)
     is zero, then infinity is returned if the corresponding denominator
     element is non-negative, else zero is returned for that element.
     DJB 4/02/04
  */
  Integer i;
  double d_0;
  long l_0;
  float f_0;
  Integer i_0;

  d_0 = (double)0.0;
  l_0 = (long)0;
  f_0 = (float)0.0;
  i_0 = (int)0;
  switch(type){
  
  case C_DBL:
    for(i = 0; i<nelems; i++) {
      if(((double*)pA)[i] == d_0) {
	if(((double*)pB)[i]>=d_0) {
	  ((double*)pC)[i]=  GA_INFINITY_D;
	} else {
	  ((double*)pC)[i]=  d_0;
	}
      } else {
	if(((double*)pB)[i]!=d_0)
	  ((double*)pC)[i]=  ((double*)pA)[i]/((double*)pB)[i];
	else{
	  /* if b is zero an infinite number could be added without
	     changing the sign of a.
	  */
	  ((double*)pC)[i]=  GA_INFINITY_D;
	}
      }
    }
    break;
  case C_DCPL:
    pnga_error(" do_step_divide called with type C_DCPL",C_DCPL);
    break;
  case C_SCPL:
    pnga_error(" do_step_divide called with type C_SCPL",C_SCPL);
    break;
  case C_INT:
    i_0 = (int)0;
    for(i = 0; i<nelems; i++){
      if(((int*)pA)[i]==i_0) {
	if(((int*)pB)[i]>=i_0) {
	  ((int*)pC)[i]=GA_INFINITY_I;
	} else {
	  ((int*)pC)[i]=i_0;
	}
      } else {
	if(((int*)pB)[i]!=i_0)
	  ((int*)pC)[i] = ((int*)pA)[i]/((int*)pB)[i];
	else{
	  ((int*)pC)[i]=GA_INFINITY_I;
	}
      } 
    }
    break;
  case C_FLOAT:
    f_0 = (float)0.0;
    for(i = 0; i<nelems; i++){
      if(((float*)pA)[i]==f_0) {
	if(((float*)pB)[i]>=f_0) {
	  ((float*)pC)[i]= GA_INFINITY_F;
	} else {
	  ((float*)pC)[i]= f_0;
	}
      } else {
	if(((float*)pB)[i]!=f_0) {
	  ((float*)pC)[i]=  ((float*)pA)[i]/((float*)pB)[i];
	} else {
	/* _F added 01/24/04 */
	  ((float*)pC)[i]= GA_INFINITY_F;
	}
      }
    }
    break;
  case C_LONG:
    l_0 = (long)0;
    for(i = 0; i<nelems; i++){
      if(((long *)pA)[i]==l_0) {
	if(((long *)pB)[i]>=l_0) {
	  ((long *)pC)[i] = GA_INFINITY_L;
	} else {
	  ((long *)pC)[i] = l_0;
	}
      } else {
	if(((long *)pB)[i]!=l_0)
	  ((long *)pC)[i]=  ((long *)pA)[i]/((long *)pB)[i];
	else{
	  ((long *)pC)[i] = GA_INFINITY_L;
	}
      }
    }
    break;		
  default: pnga_error(" wrong data type ",type);
  }
}

static void do_stepb_divide(void *pA, void *pB, void *pC, Integer nelems, Integer type){
  /* Elementwise divide, not aborting on a zero denominator, but
     returning an infinity. If an element in the numerator vector (PA)
     is zero, then infinity is returned if the corresponding denominator
     element is non-negative, else zero is returned for that element.
     DJB 4/02/04
  */
  Integer i;
  double d_0;
  long l_0;
  float f_0;
  Integer i_0;

  d_0 = (double)0.0;
  l_0 = (long)0;
  f_0 = (float)0.0;
  i_0 = (int)0;
  switch(type){
  
  case C_DBL:
    for(i = 0; i<nelems; i++) {
      if(((double*)pA)[i] == d_0) {
	if(((double*)pB)[i]>d_0) {
	  ((double*)pC)[i]=  d_0;
	} else {
	  ((double*)pC)[i]=  GA_INFINITY_D;
	}
      } else {
	if(((double*)pB)[i]!=d_0)
	  ((double*)pC)[i]=  ((double*)pA)[i]/((double*)pB)[i];
	else{
	  /* if b is zero an infinite number could be added without
	     changing the sign of a.
	  */
	  ((double*)pC)[i]=  GA_INFINITY_D;
	}
      }
    }
    break;
  case C_DCPL:
    pnga_error(" do_stepb_divide called with type C_DCPL",C_DCPL);
    break;
  case C_SCPL:
    pnga_error(" do_stepb_divide called with type C_SCPL",C_SCPL);
    break;
  case C_INT:
    i_0 = (int)0;
    for(i = 0; i<nelems; i++){
      if(((int*)pA)[i]==i_0) {
	if(((int*)pB)[i]>i_0) {
	  ((int*)pC)[i]=i_0;
	} else {
	  ((int*)pC)[i]=GA_INFINITY_I;
	}
      } else {
	if(((int*)pB)[i]!=i_0)
	  ((int*)pC)[i] = ((int*)pA)[i]/((int*)pB)[i];
	else{
	  ((int*)pC)[i]=GA_INFINITY_I;
	}
      } 
    }
    break;
  case C_FLOAT:
    f_0 = (float)0.0;
    for(i = 0; i<nelems; i++){
      if(((float*)pA)[i]==f_0) {
	if(((float*)pB)[i]>f_0) {
	  ((float*)pC)[i]= f_0;
	} else {
	  ((float*)pC)[i]= GA_INFINITY_F;
	}
      } else {
	if(((float*)pB)[i]!=f_0) {
	  ((float*)pC)[i]=  ((float*)pA)[i]/((float*)pB)[i];
	} else {
	/* _F added 01/24/04 */
	  ((float*)pC)[i]= GA_INFINITY_F;
	}
      }
    }
    break;
  case C_LONG:
    l_0 = (long)0;
    for(i = 0; i<nelems; i++){
      if(((long *)pA)[i]==l_0) {
	if(((long *)pB)[i]>l_0) {
	  ((long *)pC)[i] = l_0;
	} else {
	  ((long *)pC)[i] = GA_INFINITY_L;
	}
      } else {
	if(((long *)pB)[i]!=l_0)
	  ((long *)pC)[i]=  ((long *)pA)[i]/((long *)pB)[i];
	else{
	  ((long *)pC)[i] = GA_INFINITY_L;
	}
      }
    }
    break;		
  default: pnga_error(" do_stepb_divide: wrong data type ",type);
  }
}

static void do_step_mask(void *pA, void *pB, void *pC, Integer nelems, Integer type){
  /* 
    Set vector C to vector B wherever vector A is nonzero,
    and to zero wherever vector A is zero.
  */
  Integer i;

  switch(type){
  
  case C_DBL:
    for(i = 0; i<nelems; i++) {
      if(((double*)pA)[i]!=(double)0.0) {
	((double*)pC)[i]=  ((double*)pB)[i];
      }else {
	((double*)pC)[i]=  (double)0.0;
      }
    }
    break;
  case C_DCPL:
    pnga_error(" do_step_mask called with type C_DCPL",C_DCPL);
    break;
  case C_SCPL:
    pnga_error(" do_step_mask called with type C_SCPL",C_SCPL);
    break;
  case C_INT:
    for(i = 0; i<nelems; i++){
      if(((int*)pA)[i]!=(int)0){
	((int*)pC)[i] = ((int*)pB)[i];
      }else{
	((int*)pC)[i] = (int)0;
      } 
    }
    break;
  case C_FLOAT:
    for(i = 0; i<nelems; i++){
      if(((float*)pA)[i]!=(float)0.0) { 
	((float*)pC)[i]=  ((float*)pB)[i];
      }else{
	((float*)pC)[i]=  (float)0.0;
      }
    }
    break;
  case C_LONG:
    for(i = 0; i<nelems; i++){
      if(((long *)pA)[i]!=(long)0){
	((long *)pC)[i]=  ((long *)pB)[i];
      }else{
	((long *)pC)[i]=(long)0;  
      }
    }
    break;		
  default: pnga_error(" do_step_mask: wrong data type ",type);
  }
}


static void do_maximum(void *pA, void *pB, void *pC, Integer nelems, Integer type){
  /*
    This routine was modified by DJB to scale components
    so as not to unecessarily overflow.
  */
  Integer i;
  double aReal, aImag, bReal, bImag, temp1, temp2;
  double x1,x2;
  switch(type){
    
  case C_DBL:
    for(i = 0; i<nelems; i++)
      ((double*)pC)[i] = GA_MAX(((double*)pA)[i],((double*)pB)[i]);
    break;
  case C_DCPL:
    for(i = 0; i<nelems; i++) {
      aReal = ((DoubleComplex*)pA)[i].real;
      bReal = ((DoubleComplex*)pB)[i].real;
      aImag = ((DoubleComplex*)pA)[i].imag;
      bImag = ((DoubleComplex*)pB)[i].imag;
      x1    = GA_MAX(GA_ABS(aReal),GA_ABS(aImag));
      x2    = GA_MAX(GA_ABS(bReal),GA_ABS(bImag));
      x1    = GA_MAX(x1,x2);
      if (x1 == (double)0.0) {
	((DoubleComplex*)pC)[i].real=((DoubleComplex*)pA)[i].real;
	((DoubleComplex*)pC)[i].imag=((DoubleComplex*)pA)[i].imag;
      } else {
	x1 = ((double)1.0)/x1;
	aReal = aReal*x1;
	aImag = aImag*x1;
	bReal = bReal*x1;
	bImag = bImag*x1;
	temp1 = (aReal*aReal)+(aImag*aImag);
	temp2 = (bReal*bReal)+(bImag*bImag);
	if(temp1>temp2){
	  ((DoubleComplex*)pC)[i].real=((DoubleComplex*)pA)[i].real;
	  ((DoubleComplex*)pC)[i].imag=((DoubleComplex*)pA)[i].imag;
	}
	else{
	  ((DoubleComplex*)pC)[i].real=((DoubleComplex*)pB)[i].real;
	  ((DoubleComplex*)pC)[i].imag=((DoubleComplex*)pB)[i].imag;
	}
      }
    }
    break;
  case C_SCPL:
    for(i = 0; i<nelems; i++) {
      aReal = ((SingleComplex*)pA)[i].real;
      bReal = ((SingleComplex*)pB)[i].real;
      aImag = ((SingleComplex*)pA)[i].imag;
      bImag = ((SingleComplex*)pB)[i].imag;
      x1    = GA_MAX(GA_ABS(aReal),GA_ABS(aImag));
      x2    = GA_MAX(GA_ABS(bReal),GA_ABS(bImag));
      x1    = GA_MAX(x1,x2);
      if (x1 == (double)0.0) {
	((SingleComplex*)pC)[i].real=((SingleComplex*)pA)[i].real;
	((SingleComplex*)pC)[i].imag=((SingleComplex*)pA)[i].imag;
      } else {
	x1 = ((double)1.0)/x1;
	aReal = aReal*x1;
	aImag = aImag*x1;
	bReal = bReal*x1;
	bImag = bImag*x1;
	temp1 = (aReal*aReal)+(aImag*aImag);
	temp2 = (bReal*bReal)+(bImag*bImag);
	if(temp1>temp2){
	  ((SingleComplex*)pC)[i].real=((SingleComplex*)pA)[i].real;
	  ((SingleComplex*)pC)[i].imag=((SingleComplex*)pA)[i].imag;
	}
	else{
	  ((SingleComplex*)pC)[i].real=((SingleComplex*)pB)[i].real;
	  ((SingleComplex*)pC)[i].imag=((SingleComplex*)pB)[i].imag;
	}
      }
    }
    break;
  case C_INT:
    for(i = 0; i<nelems; i++)
      ((int*)pC)[i] =GA_MAX(((int*)pA)[i],((int*)pB)[i]);
    break;
  case C_FLOAT:
    for(i = 0; i<nelems; i++)
      ((float*)pC)[i]=GA_MAX(((float*)pA)[i],((float*)pB)[i]);
    break;
    
  case C_LONG:
    for(i = 0; i<nelems; i++)
      ((long *)pC)[i]=GA_MAX(((long *)pA)[i],((long *)pB)[i]);
    break;
    
  default: pnga_error(" wrong data type ",type);
  }
}


static void do_minimum(void *pA, void *pB, void *pC, Integer nelems, Integer type){
  /*
    This routine was modified by DJB to scale components
    so as not to unecessarily overflow.
  */
  Integer i;
  double x1,x2;

  switch(type){
    double aReal, aImag, bReal, bImag, temp1, temp2;
    
  case C_DBL:
    for(i = 0; i<nelems; i++)
      ((double*)pC)[i] = GA_MIN(((double*)pA)[i],((double*)pB)[i]);
    break;
  case C_DCPL:
    for(i = 0; i<nelems; i++) {
      aReal = ((DoubleComplex*)pA)[i].real;
      bReal = ((DoubleComplex*)pB)[i].real;
      aImag = ((DoubleComplex*)pA)[i].imag;
      bImag = ((DoubleComplex*)pB)[i].imag;
      x1    = GA_MAX(GA_ABS(aReal),GA_ABS(aImag));
      x2    = GA_MAX(GA_ABS(bReal),GA_ABS(bImag));
      x1    = GA_MAX(x1,x2);
      if (x1 == (double)0.0) {
	((DoubleComplex*)pC)[i].real=((DoubleComplex*)pA)[i].real;
	((DoubleComplex*)pC)[i].imag=((DoubleComplex*)pA)[i].imag;
      } else {
	x1 = ((double)1.0)/x1;
	aReal = aReal*x1;
	aImag = aImag*x1;
	bReal = bReal*x1;
	bImag = bImag*x1;
	temp1 = aReal*aReal+aImag*aImag;
	temp2 = bReal*bReal+bImag*bImag;
	if(temp1<temp2){ 
	  ((DoubleComplex*)pC)[i].real=((DoubleComplex*)pA)[i].real; 
	  ((DoubleComplex*)pC)[i].imag=((DoubleComplex*)pA)[i].imag; 
	} 
	else{ 
	  ((DoubleComplex*)pC)[i].real=((DoubleComplex*)pB)[i].real; 
	  ((DoubleComplex*)pC)[i].imag=((DoubleComplex*)pB)[i].imag; 
	}
      }
    }
    break;
  case C_SCPL:
    for(i = 0; i<nelems; i++) {
      aReal = ((SingleComplex*)pA)[i].real;
      bReal = ((SingleComplex*)pB)[i].real;
      aImag = ((SingleComplex*)pA)[i].imag;
      bImag = ((SingleComplex*)pB)[i].imag;
      x1    = GA_MAX(GA_ABS(aReal),GA_ABS(aImag));
      x2    = GA_MAX(GA_ABS(bReal),GA_ABS(bImag));
      x1    = GA_MAX(x1,x2);
      if (x1 == (double)0.0) {
	((SingleComplex*)pC)[i].real=((SingleComplex*)pA)[i].real;
	((SingleComplex*)pC)[i].imag=((SingleComplex*)pA)[i].imag;
      } else {
	x1 = ((double)1.0)/x1;
	aReal = aReal*x1;
	aImag = aImag*x1;
	bReal = bReal*x1;
	bImag = bImag*x1;
	temp1 = aReal*aReal+aImag*aImag;
	temp2 = bReal*bReal+bImag*bImag;
	if(temp1<temp2){ 
	  ((SingleComplex*)pC)[i].real=((SingleComplex*)pA)[i].real; 
	  ((SingleComplex*)pC)[i].imag=((SingleComplex*)pA)[i].imag; 
	} 
	else{ 
	  ((SingleComplex*)pC)[i].real=((SingleComplex*)pB)[i].real; 
	  ((SingleComplex*)pC)[i].imag=((SingleComplex*)pB)[i].imag; 
	}
      }
    }
    break;
  case C_INT:
    for(i = 0; i<nelems; i++)
      ((int*)pC)[i] =GA_MIN(((int*)pA)[i],((int*)pB)[i]);
    break;
  case C_FLOAT:
    for(i = 0; i<nelems; i++)
      ((float*)pC)[i]=GA_MIN(((float*)pA)[i],((float*)pB)[i]);
    break;
  case C_LONG:
    for(i = 0; i<nelems; i++)
      ((long *)pC)[i]=GA_MIN(((long *)pA)[i],((long *)pB)[i]);
    break;
    
  default: pnga_error(" wrong data type ",type);
  }
} 

static
void ngai_do_elem2_oper(Integer atype, Integer cndim, Integer *loC, Integer *hiC,
                        Integer *ldC, void *A_ptr, void *B_ptr, void *C_ptr, int op)
{
  Integer i, j;
  Integer bvalue[MAXDIM], bunit[MAXDIM], baseldC[MAXDIM];
  void *tempA = NULL, *tempB = NULL, *tempC = NULL;
  Integer idx, n1dim;
  /* compute "local" operation accoording to op */

  /* number of n-element of the first dimension */
  n1dim = 1; for(i=1; i<cndim; i++) n1dim *= (hiC[i] - loC[i] + 1);

  /* calculate the destination indices */
  bvalue[0] = 0; bvalue[1] = 0; bunit[0] = 1; bunit[1] = 1;
  /* baseld[0] = ld[0]
   * baseld[1] = ld[0] * ld[1]
   * baseld[2] = ld[0] * ld[1] * ld[2] .....
   */
  baseldC[0] = ldC[0]; baseldC[1] = baseldC[0] *ldC[1];
  for(i=2; i<cndim; i++) {
    bvalue[i] = 0;
    bunit[i] = bunit[i-1] * (hiC[i-1] - loC[i-1] + 1);
    baseldC[i] = baseldC[i-1] * ldC[i];
  }


  for(i=0; i<n1dim; i++) {
    idx = 0;
    for(j=1; j<cndim; j++) {
      idx += bvalue[j] * baseldC[j-1];
      if(((i+1) % bunit[j]) == 0) bvalue[j]++;
      if(bvalue[j] > (hiC[j]-loC[j])) bvalue[j] = 0;
    }

    switch(atype){
      case C_DBL:
        tempA=((double*)A_ptr)+idx;
        tempB=((double*)B_ptr)+idx;
        tempC=((double*)C_ptr)+idx;
        break;
      case C_DCPL:
        tempA=((DoubleComplex*)A_ptr)+idx;
        tempB=((DoubleComplex*)B_ptr)+idx;
        tempC=((DoubleComplex*)C_ptr)+idx;
        break;
      case C_SCPL:
        tempA=((SingleComplex*)A_ptr)+idx;
        tempB=((SingleComplex*)B_ptr)+idx;
        tempC=((SingleComplex*)C_ptr)+idx;
        break;
      case C_INT:
        tempA=((int*)A_ptr)+idx;
        tempB=((int*)B_ptr)+idx;
        tempC=((int*)C_ptr)+idx;
        break;
      case C_FLOAT:
        tempA=((float*)A_ptr)+idx;
        tempB=((float*)B_ptr)+idx;
        tempC=((float*)C_ptr)+idx;
        break;
      case C_LONG:
        tempA=((long *)A_ptr)+idx;
        tempB=((long *)B_ptr)+idx;
        tempC=((long *)C_ptr)+idx;
        break;

      default: pnga_error(" wrong data type ",atype);
    }   
    switch((int)op)
    {
      case OP_ELEM_MULT:
        do_multiply(tempA,tempB,tempC,hiC[0]-loC[0]+1,atype);
        break;
      case OP_ELEM_DIV:
        do_divide(tempA,tempB,tempC,hiC[0]-loC[0]+1,atype);
        break;
      case OP_ELEM_SDIV:
        do_step_divide(tempA,tempB,tempC,hiC[0]-loC[0]+1,atype);
        break;
      case OP_ELEM_SDIV2:
        do_stepb_divide(tempA,tempB,tempC,hiC[0]-loC[0]+1,atype);
        break;
      case OP_STEP_MASK:
        do_step_mask(tempA,tempB,tempC,hiC[0]-loC[0]+1,atype);
        break;
      case  OP_ELEM_MAX:
        do_maximum(tempA,tempB,tempC,hiC[0]-loC[0]+1,atype);
        break;
      case  OP_ELEM_MIN:
        do_minimum(tempA,tempB,tempC,hiC[0]-loC[0]+1,atype);
        break;
      default: 
        printf("op : OP_ELEM_MULT = %d:%d\n", op, OP_ELEM_MULT);
        pnga_error(" wrong operation ",op);
    }
  }
}

/*\  generic operation of two patches
\*/
static void ngai_elem2_patch_(g_a, alo, ahi, g_b, blo, bhi,
                              g_c, clo, chi, op)
Integer g_a, *alo, *ahi;    /* patch of g_a */
Integer g_b, *blo, *bhi;    /* patch of g_b */
Integer g_c, *clo, *chi;    /* patch of g_c */
int op; /* operation to be perform between g_a and g_b */
{
  Integer i, j;
  Integer compatible;
  Integer atype, btype, ctype;
  Integer andim, adims[MAXDIM], bndim, bdims[MAXDIM], cndim, cdims[MAXDIM];
  Integer loA[MAXDIM], hiA[MAXDIM], ldA[MAXDIM];
  Integer loB[MAXDIM], hiB[MAXDIM], ldB[MAXDIM];
  Integer loC[MAXDIM], hiC[MAXDIM], ldC[MAXDIM];
  void *A_ptr, *B_ptr, *C_ptr;
  Integer idx, n1dim;
  Integer atotal, btotal;
  Integer g_A = g_a, g_B = g_b;
  Integer num_blocks_a, num_blocks_b, num_blocks_c;
  Integer me= pnga_nodeid(), A_created=0, B_created=0;
  char *tempname = "temp", notrans='n';
  int local_sync_begin,local_sync_end;

  local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
  _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
  if(local_sync_begin)pnga_sync();
  pnga_check_handle(g_a, "gai_elem2_patch_");
  GA_PUSH_NAME("ngai_elem2_patch_");

  pnga_inquire(g_a, &atype, &andim, adims);
  pnga_inquire(g_b, &btype, &bndim, bdims);
  pnga_inquire(g_c, &ctype, &cndim, cdims);

  if(atype != btype || atype != ctype ) pnga_error(" types mismatch ", 0L); 

  /* check if patch indices and dims match */
  for(i=0; i<andim; i++)
    if(alo[i] <= 0 || ahi[i] > adims[i])
      pnga_error("g_a indices out of range ", g_a);
  for(i=0; i<bndim; i++)
    if(blo[i] <= 0 || bhi[i] > bdims[i])
      pnga_error("g_b indices out of range ", g_b);
  for(i=0; i<cndim; i++)
    if(clo[i] <= 0 || chi[i] > cdims[i])
      pnga_error("g_c indices out of range ", g_c);

  /* check if numbers of elements in patches match each other */
  n1dim = 1; for(i=0; i<cndim; i++) n1dim *= (chi[i] - clo[i] + 1);
  atotal = 1; for(i=0; i<andim; i++) atotal *= (ahi[i] - alo[i] + 1);
  btotal = 1; for(i=0; i<bndim; i++) btotal *= (bhi[i] - blo[i] + 1);

  if((atotal != n1dim) || (btotal != n1dim))
    pnga_error("  capacities of patches do not match ", 0L);

  num_blocks_a = pnga_total_blocks(g_a);
  num_blocks_b = pnga_total_blocks(g_b);
  num_blocks_c = pnga_total_blocks(g_c);

  if (num_blocks_a < 0 && num_blocks_b < 0 && num_blocks_c < 0) {
    /* find out coordinates of patches of g_a, g_b and g_c that I own */
    pnga_distribution(g_A, me, loA, hiA);
    pnga_distribution(g_B, me, loB, hiB);
    pnga_distribution(g_c, me, loC, hiC);

    /* test if the local portion of patches matches */
    if(pnga_comp_patch(andim, loA, hiA, cndim, loC, hiC) &&
        pnga_comp_patch(andim, alo, ahi, cndim, clo, chi)) compatible = 1;
    else compatible = 0;
    pnga_gop(pnga_type_f2c(MT_F_INT), &compatible, 1, "*");
    if(!compatible) {
      /* either patches or distributions do not match:
       *        - create a temp array that matches distribution of g_c
       *        - do C<= A
       */
      if(g_b != g_c) {
        pnga_copy_patch(&notrans, g_a, alo, ahi, g_c, clo, chi);
        andim = cndim;
        g_A = g_c;
        pnga_distribution(g_A, me, loA, hiA);
      }
      else {
        if (!pnga_duplicate(g_c, &g_A, tempname))
          pnga_error("ga_dadd_patch: dup failed", 0L);
        pnga_copy_patch(&notrans, g_a, alo, ahi, g_A, clo, chi);
        andim = cndim;
        A_created = 1;
        pnga_distribution(g_A, me, loA, hiA);
      }
    }

    /* test if the local portion of patches matches */
    if(pnga_comp_patch(bndim, loB, hiB, cndim, loC, hiC) &&
        pnga_comp_patch(bndim, blo, bhi, cndim, clo, chi)) compatible = 1;
    else compatible = 0;
    pnga_gop(pnga_type_f2c(MT_F_INT), &compatible, 1, "*");
    if(!compatible) {
      /* either patches or distributions do not match:
       *        - create a temp array that matches distribution of g_c
       *        - copy & reshape patch of g_b into g_B
       */
      if (!pnga_duplicate(g_c, &g_B, tempname))
        pnga_error("ga_dadd_patch: dup failed", 0L);
      pnga_copy_patch(&notrans, g_b, blo, bhi, g_B, clo, chi);
      bndim = cndim;
      B_created = 1;
      pnga_distribution(g_B, me, loB, hiB);
    }        

    if(andim > bndim) cndim = bndim;
    if(andim < bndim) cndim = andim;

    if(!pnga_comp_patch(andim, loA, hiA, cndim, loC, hiC))
      pnga_error(" A patch mismatch ", g_A); 
    if(!pnga_comp_patch(bndim, loB, hiB, cndim, loC, hiC))
      pnga_error(" B patch mismatch ", g_B);

    /*  determine subsets of my patches to access  */
    if (pnga_patch_intersect(clo, chi, loC, hiC, cndim)){
      pnga_access_ptr(g_A, loC, hiC, &A_ptr, ldA);
      pnga_access_ptr(g_B, loC, hiC, &B_ptr, ldB);
      pnga_access_ptr(g_c, loC, hiC, &C_ptr, ldC);

      /* compute "local" operation accoording to op */
      ngai_do_elem2_oper(atype, cndim, loC, hiC, ldC, A_ptr, B_ptr, C_ptr, op);

      /* release access to the data */
      pnga_release       (g_A, loC, hiC);
      pnga_release       (g_B, loC, hiC); 
      pnga_release_update(g_c, loC, hiC); 

    }
  } else {
    /* create copies of arrays A and B that are identically distributed
       as C*/
    if (!pnga_duplicate(g_c, &g_A, tempname))
      pnga_error("ga_dadd_patch: dup failed", 0L);
    pnga_copy_patch(&notrans, g_a, alo, ahi, g_A, clo, chi);
    andim = cndim;
    A_created = 1;

    if (!pnga_duplicate(g_c, &g_B, tempname))
      pnga_error("ga_dadd_patch: dup failed", 0L);
    pnga_copy_patch(&notrans, g_b, blo, bhi, g_B, clo, chi);
    bndim = cndim;
    B_created = 1;

    /* C is normally distributed so just add copies together for regular
       arrays */
    if (num_blocks_c < 0) {
      pnga_distribution(g_c, me, loC, hiC);
      if(andim > bndim) cndim = bndim;
      if(andim < bndim) cndim = andim;
      if (pnga_patch_intersect(clo, chi, loC, hiC, cndim)){
        pnga_access_ptr(g_A, loC, hiC, &A_ptr, ldA);
        pnga_access_ptr(g_B, loC, hiC, &B_ptr, ldB);
        pnga_access_ptr(g_c, loC, hiC, &C_ptr, ldC);

        /* compute "local" operation accoording to op */
        ngai_do_elem2_oper(atype, cndim, loC, hiC, ldC, A_ptr, B_ptr, C_ptr, op);

        /* release access to the data */
        pnga_release       (g_A, loC, hiC);
        pnga_release       (g_B, loC, hiC);
        pnga_release_update(g_c, loC, hiC);
      }
    } else {
      Integer lod[MAXDIM];
      /* Integer hid[MAXDIM]; */
      /* Integer chk; */
      Integer offset, last, jtot;
      if (!pnga_uses_proc_grid(g_c)) {
        Integer nproc = pnga_nnodes();
        for (idx = me; idx < num_blocks_c; idx += nproc) {

          pnga_distribution(g_c, idx, loC, hiC);
          /* make temporary copies of loC and hiC since pnga_patch_intersect
             destroys original versions */
          for (j=0; j<cndim; j++) {
            lod[j] = loC[j];
            /* hid[j] = hiC[j]; */
          }

          if (pnga_patch_intersect(clo, chi, loC, hiC, cndim)) {
            pnga_access_block_ptr(g_A, idx, &A_ptr, ldA);
            pnga_access_block_ptr(g_B, idx, &B_ptr, ldB);
            pnga_access_block_ptr(g_c, idx, &C_ptr, ldC);

            /* evaluate offsets for system */
            offset = 0;
            last = cndim - 1;
            jtot = 1;
            for (j=0; j<last; j++) {
              offset += (loC[j] - lod[j])*jtot;
              jtot *= ldC[j];
            }
            offset += (loC[last]-lod[last])*jtot;
            switch(ctype) {
              case C_DBL:
                A_ptr = (void*)((double*)(A_ptr) + offset);
                B_ptr = (void*)((double*)(B_ptr) + offset);
                C_ptr = (void*)((double*)(C_ptr) + offset);
                break;
              case C_INT:
                A_ptr = (void*)((int*)(A_ptr) + offset);
                B_ptr = (void*)((int*)(B_ptr) + offset);
                C_ptr = (void*)((int*)(C_ptr) + offset);
                break;
              case C_DCPL:
                A_ptr = (void*)((DoubleComplex*)(A_ptr) + offset);
                B_ptr = (void*)((DoubleComplex*)(B_ptr) + offset);
                C_ptr = (void*)((DoubleComplex*)(C_ptr) + offset);
                break;
              case C_SCPL:
                A_ptr = (void*)((SingleComplex*)(A_ptr) + offset);
                B_ptr = (void*)((SingleComplex*)(B_ptr) + offset);
                C_ptr = (void*)((SingleComplex*)(C_ptr) + offset);
                break;
              case C_FLOAT:
                A_ptr = (void*)((float*)(A_ptr) + offset);
                B_ptr = (void*)((float*)(B_ptr) + offset);
                C_ptr = (void*)((float*)(C_ptr) + offset);
                break;
              case C_LONG:
                A_ptr = (void*)((long*)(A_ptr) + offset);
                B_ptr = (void*)((long*)(B_ptr) + offset);
                C_ptr = (void*)((long*)(C_ptr) + offset);
                break;
              default:
                break;
            }

            /* compute "local" operation accoording to op */
            ngai_do_elem2_oper(atype, cndim, loC, hiC, ldC, A_ptr, B_ptr, C_ptr, op);

            /* release access to the data */
            pnga_release_block       (g_A, idx);
            pnga_release_block       (g_B, idx);
            pnga_release_update_block(g_c, idx);
          }
        }
      } else {
        /* Uses scalapack block-cyclic data distribution */
        Integer proc_index[MAXDIM], index[MAXDIM];
        Integer topology[MAXDIM];
        Integer blocks[MAXDIM], block_dims[MAXDIM];
        pnga_get_proc_index(g_c, me, proc_index);
        pnga_get_proc_index(g_c, me, index);
        pnga_get_block_info(g_c, blocks, block_dims);
        pnga_get_proc_grid(g_c, topology);
        while (index[cndim-1] < blocks[cndim-1]) {
          /* find bounding coordinates of block */
          /* chk = 1; */
          for (i = 0; i < cndim; i++) {
            loC[i] = index[i]*block_dims[i]+1;
            hiC[i] = (index[i] + 1)*block_dims[i];
            if (hiC[i] > cdims[i]) hiC[i] = cdims[i];
            /* if (hiC[i] < loC[i]) chk = 0; */
          }
          /* make temporary copies of loC and hiC since pnga_patch_intersect
             destroys original versions */
          for (j=0; j<cndim; j++) {
            lod[j] = loC[j];
            /* hid[j] = hiC[j]; */
          }

          if (pnga_patch_intersect(clo, chi, loC, hiC, cndim)) {
            pnga_access_block_grid_ptr(g_A, index, &A_ptr, ldA);
            pnga_access_block_grid_ptr(g_B, index, &B_ptr, ldB);
            pnga_access_block_grid_ptr(g_c, index, &C_ptr, ldC);

            /* evaluate offsets for system */
            offset = 0;
            last = cndim - 1;
            jtot = 1;
            for (j=0; j<last; j++) {
              offset += (loC[j] - lod[j])*jtot;
              jtot *= ldC[j];
            }
            offset += (loC[last]-lod[last])*jtot;
            switch(ctype) {
              case C_DBL:
                A_ptr = (void*)((double*)(A_ptr) + offset);
                B_ptr = (void*)((double*)(B_ptr) + offset);
                C_ptr = (void*)((double*)(C_ptr) + offset);
                break;
              case C_INT:
                A_ptr = (void*)((int*)(A_ptr) + offset);
                B_ptr = (void*)((int*)(B_ptr) + offset);
                C_ptr = (void*)((int*)(C_ptr) + offset);
                break;
              case C_DCPL:
                A_ptr = (void*)((DoubleComplex*)(A_ptr) + offset);
                B_ptr = (void*)((DoubleComplex*)(B_ptr) + offset);
                C_ptr = (void*)((DoubleComplex*)(C_ptr) + offset);
                break;
              case C_SCPL:
                A_ptr = (void*)((SingleComplex*)(A_ptr) + offset);
                B_ptr = (void*)((SingleComplex*)(B_ptr) + offset);
                C_ptr = (void*)((SingleComplex*)(C_ptr) + offset);
                break;
              case C_FLOAT:
                A_ptr = (void*)((float*)(A_ptr) + offset);
                B_ptr = (void*)((float*)(B_ptr) + offset);
                C_ptr = (void*)((float*)(C_ptr) + offset);
                break;
              case C_LONG:
                A_ptr = (void*)((long*)(A_ptr) + offset);
                B_ptr = (void*)((long*)(B_ptr) + offset);
                C_ptr = (void*)((long*)(C_ptr) + offset);
                break;
              default:
                break;
            }

            /* compute "local" operation accoording to op */
            ngai_do_elem2_oper(atype, cndim, loC, hiC, ldC, A_ptr, B_ptr, C_ptr, op);

            /* release access to the data */
            pnga_release_block_grid       (g_A, index);
            pnga_release_block_grid       (g_B, index);
            pnga_release_update_block_grid(g_c, index);
          }
          /* increment index to get next block on processor */
          index[0] += topology[0];
          for (i = 0; i < cndim; i++) {
            if (index[i] >= blocks[i] && i<cndim-1) {
              index[i] = proc_index[i];
              index[i+1] += topology[i+1];
            }
          }
        }
      }
    }
  }

  if(A_created) pnga_destroy(g_A);
  if(B_created) pnga_destroy(g_B);

  GA_POP_NAME;
  if(local_sync_end)pnga_sync();
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_elem_multiply = pnga_elem_multiply
#endif
void pnga_elem_multiply(Integer g_a, Integer g_b, Integer g_c){
 
   Integer atype, andim;
   Integer btype, bndim;
   Integer ctype, cndim;
   Integer alo[MAXDIM],ahi[MAXDIM];
   Integer blo[MAXDIM],bhi[MAXDIM];
   Integer clo[MAXDIM],chi[MAXDIM];
 
    pnga_inquire(g_a,  &atype, &andim, ahi);
    pnga_inquire(g_b,  &btype, &bndim, bhi);
    pnga_inquire(g_c,  &ctype, &cndim, chi);
    if((andim!=bndim)||(andim!=cndim))
	pnga_error("global arrays have different dimmensions.", andim);
    while(andim){
        alo[andim-1]=1;
        blo[bndim-1]=1;
        clo[cndim-1]=1;
        andim--;
        bndim--;
        cndim--;
    }
    _ga_sync_begin = 1; /*just to be on the safe side*/
    ngai_elem2_patch_(g_a, alo, ahi, g_b, blo, bhi,g_c,clo,chi,OP_ELEM_MULT);

}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_elem_divide = pnga_elem_divide
#endif
void pnga_elem_divide(Integer g_a, Integer g_b, Integer g_c){
 
   Integer atype, andim;
   Integer btype, bndim;
   Integer ctype, cndim;
   Integer alo[MAXDIM],ahi[MAXDIM];
   Integer blo[MAXDIM],bhi[MAXDIM];
   Integer clo[MAXDIM],chi[MAXDIM];
 
    pnga_inquire(g_a,  &atype, &andim, ahi);
    pnga_inquire(g_b,  &btype, &bndim, bhi);
    pnga_inquire(g_c,  &ctype, &cndim, chi);
    if((andim!=bndim)||(andim!=cndim))
        pnga_error("global arrays have different dimmensions.", andim);
    while(andim){
        alo[andim-1]=1;
        blo[bndim-1]=1;
        clo[cndim-1]=1;
        andim--;
        bndim--;
        cndim--;
    }

    _ga_sync_begin = 1; /*just to be on the safe side*/
  ngai_elem2_patch_(g_a, alo, ahi, g_b, blo, bhi,g_c,clo,chi,OP_ELEM_DIV);
 
}

 


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_elem_maximum = pnga_elem_maximum
#endif
void pnga_elem_maximum(Integer g_a, Integer g_b, Integer g_c){

   Integer atype, andim;
   Integer btype, bndim;
   Integer ctype, cndim;
   Integer alo[MAXDIM],ahi[MAXDIM];
   Integer blo[MAXDIM],bhi[MAXDIM];
   Integer clo[MAXDIM],chi[MAXDIM];

    pnga_inquire(g_a,  &atype, &andim, ahi);
    pnga_inquire(g_b,  &btype, &bndim, bhi);
    pnga_inquire(g_c,  &ctype, &cndim, chi);
    if((andim!=bndim)||(andim!=cndim))
        pnga_error("global arrays have different dimmensions.", andim);
    while(andim){
        alo[andim-1]=1;
        blo[bndim-1]=1;
        clo[cndim-1]=1;
        andim--;
        bndim--;
        cndim--;
    }

    _ga_sync_begin = 1; /*just to be on the safe side*/
    ngai_elem2_patch_(g_a, alo, ahi, g_b, blo, bhi,g_c,clo,chi,OP_ELEM_MAX);

}

 
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_elem_minimum = pnga_elem_minimum
#endif
void pnga_elem_minimum(Integer g_a, Integer g_b, Integer g_c){
 
   Integer atype, andim;
   Integer btype, bndim;
   Integer ctype, cndim;
   Integer alo[MAXDIM],ahi[MAXDIM];
   Integer blo[MAXDIM],bhi[MAXDIM];
   Integer clo[MAXDIM],chi[MAXDIM];
 
    pnga_inquire(g_a,  &atype, &andim, ahi);
    pnga_inquire(g_b,  &btype, &bndim, bhi);
    pnga_inquire(g_c,  &ctype, &cndim, chi);
    if((andim!=bndim)||(andim!=cndim))
        pnga_error("global arrays have different dimmensions.", andim);
    while(andim){
        alo[andim-1]=1;
        blo[bndim-1]=1;
        clo[cndim-1]=1;
        andim--;
        bndim--;
        cndim--;
    }
 
    _ga_sync_begin = 1; /*just to be on the safe side*/
    ngai_elem2_patch_(g_a, alo, ahi, g_b, blo, bhi,g_c,clo,chi,OP_ELEM_MIN);
 
}
 
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_elem_multiply_patch = pnga_elem_multiply_patch
#endif
void pnga_elem_multiply_patch(Integer g_a,Integer *alo,Integer *ahi,Integer g_b,Integer *blo,Integer *bhi,Integer g_c,Integer *clo,Integer *chi){

    ngai_elem2_patch_(g_a, alo, ahi, g_b, blo, bhi,g_c,clo,chi,OP_ELEM_MULT);

}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_elem_divide_patch = pnga_elem_divide_patch
#endif
void pnga_elem_divide_patch(Integer g_a,Integer *alo,Integer *ahi,
Integer g_b,Integer *blo,Integer *bhi,Integer g_c, Integer *clo,Integer *chi){

    ngai_elem2_patch_(g_a, alo, ahi, g_b, blo, bhi,g_c,clo,chi,OP_ELEM_DIV);

}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_elem_step_divide_patch = pnga_elem_step_divide_patch
#endif
void pnga_elem_step_divide_patch(Integer g_a,Integer *alo,Integer *ahi,
Integer g_b,Integer *blo,Integer *bhi,Integer g_c, Integer *clo,Integer *chi){

    ngai_elem2_patch_(g_a, alo, ahi, g_b, blo, bhi,g_c,clo,chi,OP_ELEM_SDIV);

}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_elem_stepb_divide_patch = pnga_elem_stepb_divide_patch
#endif
void pnga_elem_stepb_divide_patch(Integer g_a,Integer *alo,Integer *ahi,
Integer g_b,Integer *blo,Integer *bhi,Integer g_c, Integer *clo,Integer *chi){

    ngai_elem2_patch_(g_a, alo, ahi, g_b, blo, bhi,g_c,clo,chi,OP_ELEM_SDIV2);

}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_step_mask_patch = pnga_step_mask_patch
#endif
void pnga_step_mask_patch(Integer g_a,Integer *alo,Integer *ahi,
Integer g_b,Integer *blo,Integer *bhi,Integer g_c, Integer *clo,Integer *chi){

    ngai_elem2_patch_(g_a, alo, ahi, g_b, blo, bhi,g_c,clo,chi,OP_STEP_MASK);

}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_elem_maximum_patch = pnga_elem_maximum_patch
#endif
void pnga_elem_maximum_patch(Integer g_a,Integer *alo,Integer *ahi,
Integer g_b,Integer *blo,Integer *bhi,Integer g_c,Integer *clo,Integer *chi){

    ngai_elem2_patch_(g_a, alo, ahi, g_b, blo, bhi,g_c,clo,chi,OP_ELEM_MAX);

}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_elem_minimum_patch = pnga_elem_minimum_patch
#endif
void pnga_elem_minimum_patch(Integer g_a,Integer *alo,Integer *ahi,
Integer g_b,Integer *blo,Integer *bhi,Integer g_c,Integer *clo,Integer *chi){

    ngai_elem2_patch_(g_a, alo, ahi, g_b, blo, bhi,g_c,clo,chi,OP_ELEM_MIN);

}

static
void ngai_do_elem3_patch(Integer atype, Integer andim, Integer *loA, Integer *hiA,
                         Integer *ldA, void *A_ptr, Integer op)
{
  Integer i, j;
  void *tempA = NULL;
  Integer bvalue[MAXDIM], bunit[MAXDIM], baseldA[MAXDIM];
  Integer idx, n1dim;

  /* number of n-element of the first dimension */
  n1dim = 1; for(i=1; i<andim; i++) n1dim *= (hiA[i] - loA[i] + 1);

  /* calculate the destination indices */
  bvalue[0] = 0; bvalue[1] = 0; bunit[0] = 1; bunit[1] = 1;
  /* baseld[0] = ld[0]
   * baseld[1] = ld[0] * ld[1]
   * baseld[2] = ld[0] * ld[1] * ld[2] .....
   */
  baseldA[0] = ldA[0]; baseldA[1] = baseldA[0] *ldA[1];
  for(i=2; i<andim; i++) {
    bvalue[i] = 0;
    bunit[i] = bunit[i-1] * (hiA[i-1] - loA[i-1] + 1);
    baseldA[i] = baseldA[i-1] * ldA[i];
  }

  for(i=0; i<n1dim; i++) {
    idx = 0;
    for(j=1; j<andim; j++) {
      idx += bvalue[j] * baseldA[j-1];
      if(((i+1) % bunit[j]) == 0) bvalue[j]++;
      if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
    }

    switch(atype){
      case C_DBL:
        tempA=((double*)A_ptr)+idx;
        break;
      case C_DCPL:
      case C_SCPL:
        pnga_error(" ngai_elem3_patch_: wrong data type ",atype);
        break;
      case C_INT:
        tempA=((int*)A_ptr)+idx;
        break;
      case C_FLOAT:
        tempA=((float*)A_ptr)+idx;
        break;
      case C_LONG:
        tempA=((long *)A_ptr)+idx;
        break;

      default: pnga_error(" ngai_elem3_patch_: wrong data type ",atype);
    }

    switch(op){
      case  OP_STEPMAX:
        do_stepmax(tempA,hiA[0]-loA[0]+1, atype);
        break;
      case  OP_STEPBOUNDINFO:
        do_stepboundinfo(tempA,hiA[0]-loA[0]+1, atype);
        break;
      default: pnga_error(" wrong operation ",op);
    }
  }
}

static void ngai_elem3_patch_(Integer g_a, Integer *alo, Integer *ahi, int op)
  /*do some preprocess jobs for stepMax and stepMax2*/
{
  Integer i;
  Integer atype;
  Integer andim, adims[MAXDIM];
  Integer loA[MAXDIM], hiA[MAXDIM], ldA[MAXDIM];
  void *A_ptr;
  Integer me= pnga_nodeid();
  Integer num_blocks;
  int local_sync_begin,local_sync_end;

  local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
  _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
  if(local_sync_begin)pnga_sync();

  pnga_check_handle(g_a, "gai_elem3_patch_");
  GA_PUSH_NAME("ngai_elem3_patch_");

  pnga_inquire(g_a, &atype, &andim, adims);
  num_blocks = pnga_total_blocks(g_a);

  /* check if patch indices and dims match */
  for(i=0; i<andim; i++)
    if(alo[i] <= 0 || ahi[i] > adims[i])
      pnga_error("g_a indices out of range ", g_a);

  if (num_blocks < 0) {
  /* find out coordinates of patches of g_a, g_b and g_c that I own */
    pnga_distribution(g_a, me, loA, hiA);

    /*  determine subsets of my patches to access  */
    if (pnga_patch_intersect(alo, ahi, loA, hiA, andim)){
      pnga_access_ptr(g_a, loA, hiA, &A_ptr, ldA);

      /* compute "local" operation accoording to op */
      ngai_do_elem3_patch(atype, andim, loA, hiA, ldA, A_ptr, op);

      /* release access to the data */
      pnga_release (g_a, loA, hiA);
    }
  } else {
    Integer offset, j, jtmp, chk;
    Integer loS[MAXDIM], nproc;
    nproc = pnga_nnodes();
    /* using simple block-cyclic data distribution */
    if (!pnga_uses_proc_grid(g_a)){
      for (i=me; i<num_blocks; i += nproc) {
        /* get limits of patch */
        pnga_distribution(g_a, i, loA, hiA);

        /* loA is changed by pnga_patch_intersect, so
         *            save a copy */
        for (j=0; j<andim; j++) {
          loS[j] = loA[j];
        }

        /*  determine subset of my local patch to access  */
        /*  Output is in loA and hiA */
        if(pnga_patch_intersect(alo, ahi, loA, hiA, andim)){

          /* get data_ptr to corner of patch */
          /* ld are leading dimensions for block */
          pnga_access_block_ptr(g_a, i, &A_ptr, ldA);

          /* Check for partial overlap */
          chk = 1;
          for (j=0; j<andim; j++) {
            if (loS[j] < loA[j]) {
              chk=0;
              break;
            }
          }
          if (!chk) {
            /* Evaluate additional offset for pointer */
            offset = 0;
            jtmp = 1;
            for (j=0; j<andim-1; j++) {
              offset += (loA[j]-loS[j])*jtmp;
              jtmp *= ldA[j];
            }
            offset += (loA[andim-1]-loS[andim-1])*jtmp;
            switch (atype){
              case C_INT:
                A_ptr = (void*)((int*)A_ptr + offset);
                break;
              case C_DCPL:
                A_ptr = (void*)((double*)A_ptr + 2*offset);
                break;
              case C_SCPL:
                A_ptr = (void*)((float*)A_ptr + 2*offset);
                break;
              case C_DBL:
                A_ptr = (void*)((double*)A_ptr + offset);
                break;
              case C_FLOAT:
                A_ptr = (void*)((float*)A_ptr + offset);
                break;
              case C_LONG:
                A_ptr = (void*)((long*)A_ptr + offset);
                break;
              default: pnga_error(" wrong data type ",atype);
            }
          }

          /* compute "local" operation accoording to op */
          ngai_do_elem3_patch(atype, andim, loA, hiA, ldA, A_ptr, op);

          /* release access to the data */
          pnga_release_update_block(g_a, i);
        }
      }
    } else {
      /* using scalapack block-cyclic data distribution */
      Integer proc_index[MAXDIM], index[MAXDIM];
      Integer topology[MAXDIM];
      Integer blocks[MAXDIM], block_dims[MAXDIM];
      pnga_get_proc_index(g_a, me, proc_index);
      pnga_get_proc_index(g_a, me, index);
      pnga_get_block_info(g_a, blocks, block_dims);
      pnga_get_proc_grid(g_a, topology);
      while (index[andim-1] < blocks[andim-1]) {
        /* find bounding coordinates of block */
        for (i = 0; i < andim; i++) {
          loA[i] = index[i]*block_dims[i]+1;
          hiA[i] = (index[i] + 1)*block_dims[i];
          if (hiA[i] > adims[i]) hiA[i] = adims[i];
        }
        /* loA is changed by pnga_patch_intersect, so
         *            save a copy */
        for (j=0; j<andim; j++) {
          loS[j] = loA[j];
        }

        /*  determine subset of my local patch to access  */
        /*  Output is in loA and hiA */
        if(pnga_patch_intersect(alo, ahi, loA, hiA, andim)){

          /* get data_ptr to corner of patch */
          /* ld are leading dimensions for block */
          pnga_access_block_grid_ptr(g_a, index, &A_ptr, ldA);

          /* Check for partial overlap */
          chk = 1;
          for (j=0; j<andim; j++) {
            if (loS[j] < loA[j]) {
              chk=0;
              break;
            }
          }
          if (!chk) {
            /* Evaluate additional offset for pointer */
            offset = 0;
            jtmp = 1;
            for (j=0; j<andim-1; j++) {
              offset += (loA[j]-loS[j])*jtmp;
              jtmp *= ldA[j];
            }
            offset += (loA[andim-1]-loS[andim-1])*jtmp;
            switch (atype){
              case C_INT:
                A_ptr = (void*)((int*)A_ptr + offset);
                break;
              case C_DCPL:
                A_ptr = (void*)((double*)A_ptr + 2*offset);
                break;
              case C_SCPL:
                A_ptr = (void*)((float*)A_ptr + 2*offset);
                break;
              case C_DBL:
                A_ptr = (void*)((double*)A_ptr + offset);
                break;
              case C_FLOAT:
                A_ptr = (void*)((float*)A_ptr + offset);
                break;
              case C_LONG:
                A_ptr = (void*)((long*)A_ptr + offset);
                break;
              default: pnga_error(" wrong data type ",atype);
            }
          }

          /* compute "local" operation accoording to op */
          ngai_do_elem3_patch(atype, andim, loA, hiA, ldA, A_ptr, op);

          /* release access to the data */
          pnga_release_update_block_grid(g_a, index);
        }
        /* increment index to get next block on processor */
        index[0] += topology[0];
        for (i = 0; i < andim; i++) {
          if (index[i] >= blocks[i] && i<andim-1) {
            index[i] = proc_index[i];
            index[i+1] += topology[i+1];
          }
        }
      }
    }
  }

  GA_POP_NAME;
  if(local_sync_end)pnga_sync();
}

static void ngai_has_negative_element(Integer atype, Integer andim, Integer *loA, Integer *hiA, Integer *ldA, void *A_ptr, Integer *iretval)
{
  Integer i, j;
  Integer bvalue[MAXDIM], bunit[MAXDIM], baseldA[MAXDIM];
  Integer idx, n1dim;
  double *tempA;
  int    *itempA;
  long   *ltempA;
  float  *ftempA;

  /* number of n-element of the first dimension */
  n1dim = 1; for(i=1; i<andim; i++) n1dim *= (hiA[i] - loA[i] + 1);

  /* calculate the destination indices */
  bvalue[0] = 0; bvalue[1] = 0; bunit[0] = 1; bunit[1] = 1;
  /* baseld[0] = ld[0]
   * baseld[1] = ld[0] * ld[1]
   * baseld[2] = ld[0] * ld[1] * ld[2] .....
   */
  baseldA[0] = ldA[0]; baseldA[1] = baseldA[0] *ldA[1];
  for(i=2; i<andim; i++) {
    bvalue[i] = 0;
    bunit[i] = bunit[i-1] * (hiA[i-1] - loA[i-1] + 1);
    baseldA[i] = baseldA[i-1] * ldA[i];
  }

  for(i=0; i<n1dim; i++) {
    idx = 0;
    for(j=1; j<andim; j++) {
      idx += bvalue[j] * baseldA[j-1];
      if(((i+1) % bunit[j]) == 0) bvalue[j]++;
      if(bvalue[j] > (hiA[j]-loA[j])) bvalue[j] = 0;
    }

    switch(atype){
      case C_DBL:
        /*double is the only type that is handled for Tao/GA project*/
        /* 
           DJB modification to add types int, float and long.
           This operation does not make sense for complex.
         */
        tempA=((double*)A_ptr)+idx;
        for(j=0;j<hiA[0]-loA[0]+1;j++)
          if(tempA[j]<(double)0.0) *iretval=1;
        break;
      case C_DCPL:
      case C_SCPL:
        pnga_error(" has_negative_elem: wrong data type ",
            atype);
        break;
      case C_INT:
        itempA=((int*)A_ptr)+idx;
        for(j=0;j<hiA[0]-loA[0]+1;j++)
          if(itempA[j]<(int)0) *iretval=1;
        break;
      case C_FLOAT:
        ftempA=((float*)A_ptr)+idx;
        for(j=0;j<hiA[0]-loA[0]+1;j++)
          if(ftempA[j]<(float)0.0) *iretval=1;
        break;
      case C_LONG:
        ltempA=((long*)A_ptr)+idx;
        for(j=0;j<hiA[0]-loA[0]+1;j++)
          if(ltempA[j]<(long)0) *iretval=1;
        break;

      default: pnga_error(" has_negative_elem: wrong data type ",
                   atype);
    }

  }
}

static Integer has_negative_elem(g_a, alo, ahi)
Integer g_a, *alo, *ahi;    /* patch of g_a */
/*returned value: 1=found; 0 = not found*/
{
  Integer i;
  Integer atype;
  Integer andim, adims[MAXDIM];
  Integer loA[MAXDIM], hiA[MAXDIM], ldA[MAXDIM];
  void *A_ptr; 
  Integer iretval;
  Integer num_blocks;
  Integer me= pnga_nodeid();


  pnga_sync();
  pnga_check_handle(g_a, "has_negative_elem");
  GA_PUSH_NAME("has_negative_elem");

  pnga_inquire(g_a, &atype, &andim, adims);
  num_blocks = pnga_total_blocks(g_a);

  /* check if patch indices and dims match */
  for(i=0; i<andim; i++)
    if(alo[i] <= 0 || ahi[i] > adims[i])
      pnga_error("g_a indices out of range ", g_a);

  if (num_blocks < 0) {
    /* find out coordinates of patches of g_a, g_b and g_c that I own */
    pnga_distribution(g_a, me, loA, hiA);
    iretval = 0;
    /*  determine subsets of my patches to access  */
    if (pnga_patch_intersect(alo, ahi, loA, hiA, andim)){
      pnga_access_ptr(g_a, loA, hiA, &A_ptr, ldA);

      ngai_has_negative_element(atype, andim, loA, hiA, ldA, A_ptr, &iretval);

      /* release access to the data */
      pnga_release (g_a, loA, hiA);
    }
  } else {
    Integer offset, j, jtmp, chk;
    Integer loS[MAXDIM], nproc;
    nproc = pnga_nnodes();
    /* using simple block-cyclic data distribution */
    if (!pnga_uses_proc_grid(g_a)){
      for (i=me; i<num_blocks; i += nproc) {
        /* get limits of patch */
        pnga_distribution(g_a, i, loA, hiA);

        /* loA is changed by pnga_patch_intersect, so
         *            save a copy */
        for (j=0; j<andim; j++) {
          loS[j] = loA[j];
        }

        /*  determine subset of my local patch to access  */
        /*  Output is in loA and hiA */
        if(pnga_patch_intersect(alo, ahi, loA, hiA, andim)){

          /* get data_ptr to corner of patch */
          /* ld are leading dimensions for block */
          pnga_access_block_ptr(g_a, i, &A_ptr, ldA);

          /* Check for partial overlap */
          chk = 1;
          for (j=0; j<andim; j++) {
            if (loS[j] < loA[j]) {
              chk=0;
              break;
            }
          }
          if (!chk) {
            /* Evaluate additional offset for pointer */
            offset = 0;
            jtmp = 1;
            for (j=0; j<andim-1; j++) {
              offset += (loA[j]-loS[j])*jtmp;
              jtmp *= ldA[j];
            }
            offset += (loA[andim-1]-loS[andim-1])*jtmp;
            switch (atype){
              case C_INT:
                A_ptr = (void*)((int*)A_ptr + offset);
                break;
              case C_DCPL:
                A_ptr = (void*)((double*)A_ptr + 2*offset);
                break;
              case C_SCPL:
                A_ptr = (void*)((float*)A_ptr + 2*offset);
                break;
              case C_DBL:
                A_ptr = (void*)((double*)A_ptr + offset);
                break;
              case C_FLOAT:
                A_ptr = (void*)((float*)A_ptr + offset);
                break;
              case C_LONG:
                A_ptr = (void*)((long*)A_ptr + offset);
                break;
              default: pnga_error(" wrong data type ",atype);
            }
          }

          /* check all values in patch */
          ngai_has_negative_element(atype, andim, loA, hiA, ldA, A_ptr, &iretval);

          /* release access to the data */
          pnga_release_update_block(g_a, i);
        }
      }
    } else {
      /* using scalapack block-cyclic data distribution */
      Integer proc_index[MAXDIM], index[MAXDIM];
      Integer topology[MAXDIM];
      Integer blocks[MAXDIM], block_dims[MAXDIM];
      pnga_get_proc_index(g_a, me, proc_index);
      pnga_get_proc_index(g_a, me, index);
      pnga_get_block_info(g_a, blocks, block_dims);
      pnga_get_proc_grid(g_a, topology);
      while (index[andim-1] < blocks[andim-1]) {
        /* find bounding coordinates of block */
        for (i = 0; i < andim; i++) {
          loA[i] = index[i]*block_dims[i]+1;
          hiA[i] = (index[i] + 1)*block_dims[i];
          if (hiA[i] > adims[i]) hiA[i] = adims[i];
        }
        /* loA is changed by pnga_patch_intersect, so
         *            save a copy */
        for (j=0; j<andim; j++) {
          loS[j] = loA[j];
        }

        /*  determine subset of my local patch to access  */
        /*  Output is in loA and hiA */
        if(pnga_patch_intersect(alo, ahi, loA, hiA, andim)){

          /* get data_ptr to corner of patch */
          /* ld are leading dimensions for block */
          pnga_access_block_grid_ptr(g_a, index, &A_ptr, ldA);

          /* Check for partial overlap */
          chk = 1;
          for (j=0; j<andim; j++) {
            if (loS[j] < loA[j]) {
              chk=0;
              break;
            }
          }
          if (!chk) {
            /* Evaluate additional offset for pointer */
            offset = 0;
            jtmp = 1;
            for (j=0; j<andim-1; j++) {
              offset += (loA[j]-loS[j])*jtmp;
              jtmp *= ldA[j];
            }
            offset += (loA[andim-1]-loS[andim-1])*jtmp;
            switch (atype){
              case C_INT:
                A_ptr = (void*)((int*)A_ptr + offset);
                break;
              case C_DCPL:
                A_ptr = (void*)((double*)A_ptr + 2*offset);
                break;
              case C_SCPL:
                A_ptr = (void*)((float*)A_ptr + 2*offset);
                break;
              case C_DBL:
                A_ptr = (void*)((double*)A_ptr + offset);
                break;
              case C_FLOAT:
                A_ptr = (void*)((float*)A_ptr + offset);
                break;
              case C_LONG:
                A_ptr = (void*)((long*)A_ptr + offset);
                break;
              default: pnga_error(" wrong data type ",atype);
            }
          }

          /* check all values in patch */
          ngai_has_negative_element(atype, andim, loA, hiA, ldA, A_ptr, &iretval);

          /* release access to the data */
          pnga_release_update_block_grid(g_a, index);
        }

        /* increment index to get next block on processor */
        index[0] += topology[0];
        for (i = 0; i < andim; i++) {
          if (index[i] >= blocks[i] && i<andim-1) {
            index[i] = proc_index[i];
            index[i+1] += topology[i+1];
          }
        }
      }
    }
  }

  GA_POP_NAME;
  pnga_sync();
  return iretval; /*negative element is not found in g_a*/
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_step_bound_info_patch = pnga_step_bound_info_patch
#endif
void pnga_step_bound_info_patch(
     Integer g_xx, Integer *xxlo, Integer *xxhi,    /* patch of g_xx */
     Integer g_vv, Integer *vvlo, Integer *vvhi,    /* patch of g_vv */
     Integer g_xxll, Integer *xxlllo, Integer *xxllhi,    /* patch of g_xxll */
     Integer g_xxuu, Integer *xxuulo, Integer *xxuuhi,    /* patch of g_xxuu */
     void *boundmin, void* wolfemin, void *boundmax)
{
     /*double  result1,result2;*/
     double  dresult,dresult2;
     long    lresult,lresult2;

     Integer index[MAXDIM];
     Integer xxtype;
     Integer xxndim, xxdims[MAXDIM];
     Integer loXX[MAXDIM], hiXX[MAXDIM];
     Integer vvtype;
     Integer vvndim, vvdims[MAXDIM];
     Integer loVV[MAXDIM], hiVV[MAXDIM];
     Integer xxtotal,vvtotal;
     Integer xltype;
     Integer xlndim, xldims[MAXDIM];
     Integer loXL[MAXDIM], hiXL[MAXDIM];
     Integer xutype;
     Integer xundim, xudims[MAXDIM];
     Integer loXU[MAXDIM], hiXU[MAXDIM];
     Integer xltotal,xutotal;
     Integer me= pnga_nodeid();
     Integer g_Q;
     Integer g_R;
     Integer g_S;
     Integer g_T;
     double dalpha = (double)1.0, dbeta = (double)(-1.0);
     long   lalpha = (long)1, lbeta = (long)(-1);
     int ialpha = (int)1, ibeta = (int)(-1);
     float   falpha = (float)1.0, fbeta = (float)(-1.0);
     int iresult,iresult2;
     float   fresult,fresult2;
     Integer compatible;
     Integer compatible2;
     Integer compatible3;
     void *sresult = NULL;
     void *sresult2 = NULL;
     void *alpha = NULL,*beta = NULL;
     int local_sync_begin,local_sync_end;
     int i;

     local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
     _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
     if(local_sync_begin)pnga_sync();

     /* Check for valid ga handles. */

     pnga_check_handle(g_xx, "pnga_step_bound_info_patch");
     pnga_check_handle(g_vv, "pnga_step_bound_info_patch");
     pnga_check_handle(g_xxll, "pnga_step_bound_info_patch");
     pnga_check_handle(g_xxuu, "pnga_step_bound_info_patch");

     GA_PUSH_NAME("pnga_step_bound_info_patch");

     /* get chaacteristics of the input ga patches */

     pnga_inquire(g_xx, &xxtype, &xxndim, xxdims);
     pnga_inquire(g_vv, &vvtype, &vvndim, vvdims);
     pnga_inquire(g_xxll, &xltype, &xlndim, xldims);
     pnga_inquire(g_xxuu, &xutype, &xundim, xudims);

     /* Check for matching types. */

     if(xxtype != vvtype) pnga_error(" pnga_step_bound_info_patch: types mismatch ", 0L); 
     if(xxtype != xltype) pnga_error(" pnga_step_bound_info_patch: types mismatch ", 0L); 
     if(xxtype != xutype) pnga_error(" pnga_step_bound_info_patch: types mismatch ", 0L); 

     /* check if patch indices and dims match */
     for(i=0; i<xxndim; i++)
       if(xxlo[i] <= 0 || xxhi[i] > xxdims[i])
	 pnga_error("ga_elem_step_bound_info_patch: g_a indices out of range ", g_xx);

     for(i=0; i<vvndim; i++)
       if(vvlo[i] <= 0 || vvhi[i] > vvdims[i])
	 pnga_error("ga_elem_step_bound_info_patch: g_a indices out of range ", g_vv);

     for(i=0; i<xlndim; i++)
       if(xxlllo[i] <= 0 || xxllhi[i] > xldims[i])
	 pnga_error("ga_elem_step_bound_info_patch: g_a indices out of range ", g_xxll);
     for(i=0; i<xundim; i++)
       if(xxuulo[i] <= 0 || xxuuhi[i] > xudims[i])
	 pnga_error("ga_elem_step_bound_info_patch: g_a indices out of range ", g_xxuu);
     
     /* check if numbers of elements in patches match each other */
     xxtotal = 1; for(i=0; i<xxndim; i++) xxtotal *= (xxhi[i] - xxlo[i] + 1);
     vvtotal = 1; for(i=0; i<vvndim; i++) vvtotal *= (vvhi[i] - vvlo[i] + 1);
     xltotal = 1; for(i=0; i<xlndim; i++) xltotal *= (xxllhi[i] - xxlllo[i] + 1);
     xutotal = 1; for(i=0; i<xundim; i++) xutotal *= (xxuuhi[i] - xxuulo[i] + 1);
 
     if(xxtotal != vvtotal)
        pnga_error(" pnga_step_bound_info_patch capacities of patches do not match ", 0L);
     if(xxtotal != xltotal)
        pnga_error(" pnga_step_bound_info_patch capacities of patches do not match ", 0L);
     if(xxtotal != xutotal)
        pnga_error(" pnga_step_bound_info_patch capacities of patches do not match ", 0L);
     /* find out coordinates of patches of g_a, and g_b that I own */
     pnga_distribution(g_xx, me, loXX, hiXX);
     pnga_distribution(g_vv, me, loVV, hiVV);
     pnga_distribution(g_xxll, me, loXL, hiXL);
     pnga_distribution(g_xxuu, me, loXU, hiXU);
     

     /* test if the local portion of patches matches */
     if(pnga_comp_patch(xxndim, loXX, hiXX, vvndim, loVV, hiVV) &&
	pnga_comp_patch(xxndim, xxlo, xxhi, vvndim, vvlo, vvhi)) {
       compatible = 1;
     }
     else {
       compatible = 0;
     }
     if(pnga_comp_patch(xxndim, loXX, hiXX, xlndim, loXL, hiXL) &&
	pnga_comp_patch(xxndim, xxlo, xxhi, xlndim, xxlllo, xxllhi)) {
       compatible2 = 1;
     }
     else {
       compatible2 = 0;
     }
     if(pnga_comp_patch(xxndim, loXX, hiXX, xundim, loXU, hiXU) &&
	pnga_comp_patch(xxndim, xxlo, xxhi, xundim, xxuulo, xxuuhi)) {
       compatible3 = 1;
     }
     else {
       compatible3 = 0;
     }
     compatible = compatible * compatible2 * compatible3;
     pnga_gop(pnga_type_f2c(MT_F_INT), &compatible, 1, "*");
     if(!compatible) {
       pnga_error(" pnga_step_bound_info_patch mismatched patchs ",0);
     }
     switch (xxtype)
       {
       case C_INT:
	 /* This should point to iresult but we use lresult
	    due to the strange implementation if pnga_select_elem.
	 */
	 sresult = &iresult;
	 sresult2 = &iresult2;
	 alpha    = &ialpha;
	 beta     = &ibeta;
	 break;
       case C_DCPL:
       case C_SCPL:
	 pnga_error("Ga_step_bound_info_patch_: unavalable for complex datatype.", 
		   xxtype);
	 break;
       case C_DBL:
	 sresult = &dresult;
	 sresult2 = &dresult2;
	 alpha    = &dalpha;
	 beta     = &dbeta;
	 break;
       case C_FLOAT:
	 sresult = &fresult;
	 sresult2 = &fresult2;
	 alpha    = &falpha;
	 beta     = &fbeta;
	 break;
       case C_LONG:
	 sresult = &lresult;
	 sresult2 = &lresult2;
	 alpha    = &lalpha;
	 beta     = &lbeta;
	 break;
       default:
	 pnga_error("Ga_step_max_patch_: alpha/beta set wrong data type.", xxtype);
       }

     /*duplicatecate an array Q to hold the temparary result */
     pnga_duplicate(g_xx, &g_Q, "TempQ");
     if(g_Q==0)
       pnga_error("pnga_step_bound_info_patch:fail to duplicate array Q", g_Q);
     
     /*duplicatecate an array R to hold the temparary result */
     pnga_duplicate(g_xx, &g_R, "TempR");
     if(g_R==0)
       pnga_error("pnga_step_bound_info_patch:fail to duplicate array R", g_R);

     /*duplicatecate an array s to hold the temparary result */
     pnga_duplicate(g_xx, &g_S, "TempS");
     if(g_S==0)
       pnga_error("pnga_step_bound_info_patch:fail to duplicate array S", g_S);
     
     /*duplicatecate an array T to hold the temparary result */
     pnga_duplicate(g_xx, &g_T, "TempT");
     if(g_T==0)
       pnga_error("pnga_step_bound_info_patch:fail to duplicate array T", g_T);

     /*First, compute xu - xx */
     pnga_add_patch(alpha, g_xxuu, xxuulo, xxuuhi, beta, g_xx, xxlo, xxhi, g_S, xxlo, xxhi); 

     /*Check for negative elements in g_s, if it has any then xxuu was
       not an upper bound, exit with error message.
     */
     if(has_negative_elem(g_S, xxlo, xxhi) == 1)
       pnga_error("pnga_step_bound_info_patch: Upper bound is not > xx.", -1);

     /* Then compute t = positve elements of vv */
     pnga_zero(g_T);
     pnga_elem_maximum(g_vv,g_T,g_T);

     /* Then, compute (xu-xx)/vv */
     pnga_elem_stepb_divide_patch(g_S, xxlo, xxhi, g_T, vvlo, vvhi, g_T, xxlo, xxhi); 

     /* Then, we will select the minimum of the array g_t*/ 
     pnga_select_elem(g_T, "min", sresult, &index[0]); 

     switch (xxtype)
       {
       case C_INT:
	 /* This should be iresult but is lresult because of
	    the strange implementation of nga_select_elem.
	 */
           /* result1 = (double)(iresult); */
           break;
       case C_DCPL:
       case C_SCPL:
	 pnga_error("Ga_step_bound_info_patch_: unavalable for complex datatype.", 
		   xxtype);
	 break;
       case C_DBL:
	 /* result1 = dresult; */
	 break;
       case C_FLOAT:
	 /* result1 = (double)fresult; */
	 break;
       case C_LONG:
	 /* result1 = (double)lresult; */
	 break;
       default:
	 pnga_error("Ga_step_bound_info_patch_: result set: wrong data type.", xxtype);
       }

     /*Now doing the same thing to get (xx-xxll)/dv */
     /*First, compute xl - xx */
     pnga_add_patch(alpha, g_xx, xxlo, xxhi, beta, g_xxll, xxlllo, xxllhi, g_Q, xxlo, xxhi); 
     /*Check for negative elements in g_s, if it has any then xxll was
       not a lower bound, exit with error message.
     */
     if(has_negative_elem(g_Q, xxlo, xxhi) == 1)
       pnga_error("pnga_step_bound_info_patch: Lower bound is not < xx.", -1);

     /* Then compute r = negative elements of vv */
     pnga_zero(g_R);
     pnga_elem_minimum(g_vv,g_R,g_R);
     pnga_abs_value(g_R);

     /* Then, compute (xx-xl)/vv */
     pnga_elem_stepb_divide_patch(g_Q, xxlo, xxhi, g_R, vvlo, vvhi, g_R, xxlo, xxhi); 
     /* Then, we will select the minimum of the array g_t*/ 
     pnga_select_elem(g_R, "min", sresult2, &index[0]); 
     switch (xxtype)
       {
       case C_INT:
	 *(int*)wolfemin = GA_ABS(GA_MIN(iresult,iresult2));
	 break;
       case C_DCPL:
       case C_SCPL:
	 pnga_error("Ga_step_bound_info_patch_: unavalable for complex datatype.", 
		   xxtype);
	 break;
       case C_DBL:
	 *(double*)wolfemin = GA_ABS(GA_MIN(dresult,dresult2));
	 break;
       case C_FLOAT:
	 *(float*)wolfemin = GA_ABS(GA_MIN(fresult,fresult2));
	 break;
       case C_LONG:
	 *(long*)wolfemin =  GA_ABS(GA_MIN(lresult,lresult2));
	 break;
       default:
	 pnga_error("Ga_step_bound_info_patch_: result2 set: wrong data type.", xxtype);
       }
     /* 
       Now set T to be the elementwise minimum of R and T. 
       So, T is infinity only where ever g_vv is zero.
     */
     pnga_elem_minimum(g_R,g_T,g_T);
     /*
       Now we want to set T to be zero whenever g_vv was zero
       and gxx coincides with either boundary vector.
       Set S to be the element-wise product of S and Q.
       It will be zero when either of them is zero.
     */
     pnga_elem_multiply(g_Q,g_S,g_S);
     /*
       Set Q to the |vv|.
     */
     pnga_copy(g_vv,g_Q);
     pnga_abs_value(g_Q);
     /* 
       Now add q and s to get a vector that is zero only
       where g_vv was zero and g_xx meets one of the
       boundary vectors.
     */
     pnga_add_patch(alpha, g_Q, xxlo, xxhi, alpha, g_S, xxlo, xxhi, g_S, xxlo, xxhi); 
     /* 
       Then use that vector as a mask to set certain
       elements of T to be zero (so we have a collection
       of the a_i and c_i elements as per the TAO StepBoundInfo
       function).
     */
     pnga_step_mask_patch(g_S,xxlo,xxhi,g_T,xxlo,xxhi,g_T,xxlo,xxhi);

     /* 
       Then, we will select the minimum of the array g_t, that will
       be boundmin .
     */ 
     pnga_select_elem(g_T, "min", sresult, &index[0]); 
     switch (xxtype)
       {
       case C_INT:
	 /* This should be iresult but is lresult because of
	    the strange implementation of nga_select_elem.
	 */
           *(int*)boundmin = iresult;
           break;
       case C_DCPL:
       case C_SCPL:
	 pnga_error("Ga_step_bound_info_patch_: unavalable for complex datatype.", 
		   xxtype);
	 break;
       case C_DBL:
	 *(double*)boundmin = dresult;
	 break;
       case C_FLOAT:
	 *(float*)boundmin = fresult;
	 break;
       case C_LONG:
	 *(long*)boundmin = lresult;
	 break;
       default:
	 pnga_error("Ga_step_bound_info_patch_: result set: wrong data type.", xxtype);
       }
     /* 
       Then, we will select the maximum of the array g_t, that will
       be boundmax .
     */ 
     pnga_select_elem(g_T, "max", sresult, &index[0]); 
     switch (xxtype)
       {
       case C_INT:
	 /* This should be iresult but is lresult because of
	    the strange implementation of nga_select_elem.
	 */
           *(int*)boundmax = iresult;
           break;
       case C_DCPL:
       case C_SCPL:
	 pnga_error("Ga_step_bound_info_patch_: unavalable for complex datatype.", 
		   xxtype);
	 break;
       case C_DBL:
	 *(double*)boundmax = dresult;
	 break;
       case C_FLOAT:
	 *(float*)boundmax = fresult;
	 break;
       case C_LONG:
	 *(long*)boundmax = lresult;
	 break;
       default:
	 pnga_error("Ga_step_bound_info_patch_: result set: wrong data type.", xxtype);
       }
     pnga_destroy(g_Q); 
     pnga_destroy(g_R); 
     pnga_destroy(g_S); 
     pnga_destroy(g_T); 
     GA_POP_NAME;
     if(local_sync_end)pnga_sync();
}

/*\ generic  routine for element wise operation between two array
\*/
#if 0 /* I want to delete op parameter */
void ga_step_max_patch_(g_a,  alo, ahi, g_b,  blo, bhi, result, op) 
#else
#endif

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_step_max_patch = pnga_step_max_patch
#endif
void pnga_step_max_patch(g_a,  alo, ahi, g_b,  blo, bhi, result) 
     Integer g_a, *alo, *ahi;    /* patch of g_a */
     Integer g_b, *blo, *bhi;    /* patch of g_b */
     void *result;
#if 0
     Integer op; /* operations */
#endif

{
  double  dresult;
  long    lresult;
  Integer atype;
  Integer andim, adims[MAXDIM];
  Integer btype;
  Integer bndim, bdims[MAXDIM];
  Integer index[MAXDIM];
  /* Integer num_blocks_a, num_blocks_b; */
  /* double result = -1; */
  Integer g_c;
  int iresult;
  Integer atotal,btotal;
  float   fresult;
  int local_sync_begin,local_sync_end;
  int i;
  Integer compatible;
  void *sresult = NULL;

  local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
  _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
  if(local_sync_begin)pnga_sync();

  /* Check for valid ga handles. */

  pnga_check_handle(g_a, "ga_step_max_patch_");
  pnga_check_handle(g_b, "ga_step_max_patch_");

  GA_PUSH_NAME("ga_step_max_patch_");

  /* get chacteristics of the input ga patches */

  pnga_inquire(g_a, &atype, &andim, adims);
  pnga_inquire(g_b, &btype, &bndim, bdims);
  /* num_blocks_a = pnga_total_blocks(g_a); */
  /* num_blocks_b = pnga_total_blocks(g_b); */

  /* Check for matching types. */
  if(atype != btype) pnga_error(" ga_step_max_patch_: types mismatch ", 0L); 

  /* check if patch indices and dims match */
  for(i=0; i<andim; i++)
    if(alo[i] <= 0 || ahi[i] > adims[i])
      pnga_error("g_a indices out of range ", g_a);
  for(i=0; i<bndim; i++)
    if(blo[i] <= 0 || bhi[i] > bdims[i])
      pnga_error("g_b indices out of range ", g_b);

  /* check if numbers of elements in patches match each other */
  atotal = 1; for(i=0; i<andim; i++) atotal *= (ahi[i] - alo[i] + 1);
  btotal = 1; for(i=0; i<bndim; i++) btotal *= (bhi[i] - blo[i] + 1);

  if(btotal != atotal)
    pnga_error(" ga_step_max_patch_ capacities of patches do not match ", 0L);

  /* test if patches match */
  if(pnga_comp_patch(andim, alo, ahi, bndim, blo, bhi)) compatible = 1;
  else compatible = 0;
  pnga_gop(pnga_type_f2c(MT_F_INT), &compatible, 1, "*");
  if(!compatible) {
    pnga_error(" ga_step_max_patch_ mismatched patchs ",0);
  }

  switch (atype)
  {
    case C_INT:
      sresult = &iresult;
      break;
    case C_DCPL:
    case C_SCPL:
      pnga_error("Ga_step_max_patch_: unavalable for complex datatype.", 
          atype);
      break;
    case C_DBL:
      sresult = &dresult;
      break;
    case C_FLOAT:
      sresult = &fresult;
      break;
    case C_LONG:
      sresult = &lresult;
      break;
    default:
      pnga_error("Ga_step_max_patch_: wrong data type.", atype);
  }

  if(g_a == g_b) {
    /* It used to say 1, but if ga and gb are the same, and 
       ga is nonnegative then any number of multiples of gb
       can be added to ga still leaving it nonnegative.
     *result = (double)1.0;
     */
    switch (atype)
    {
      case C_INT:
        *(int*)result = GA_INFINITY_I;
        break;
      case C_DCPL:
      case C_SCPL:
        pnga_error("Ga_step_max_patch_: unavailable for complex datatype.", 
            atype);
        break;
      case C_DBL:
        *(double*)result = GA_INFINITY_D;
        break;
      case C_FLOAT:
        *(float*)result = GA_INFINITY_F;
        break;
      case C_LONG:
        *(long*)result = GA_INFINITY_L;
        break;
      default:
        pnga_error("Ga_step_max_patch_: wrong data type.", atype);
    }
  } else {
    /*Now look at each element of the array g_a. 
      If an element of g_a is negative, then simply return */ 
    if(has_negative_elem(g_a, alo, ahi) == 1)
      pnga_error("ga_step_max_patch_: g_a has negative element.", -1);

    /*duplicate an array c to hold the temparate result = g_a/g_b; */
    pnga_duplicate(g_a, &g_c, "Temp");
    if(g_c==0)
      pnga_error("ga_step_max_patch_:fail to duplicate array c", g_a);

    /*
       pnga_elem_divide_patch(g_a, alo, ahi, g_b, blo, bhi, g_c, alo, ahi);
     */
    pnga_elem_step_divide_patch(g_a, alo, ahi, g_b, blo, bhi, 
        g_c, alo, ahi);

    /*Now look at each element of the array g_c. If an element of g_c is positive,
      then replace it with -GA_INFINITY */ 
    ngai_elem3_patch_(g_c, alo, ahi, OP_STEPMAX);  
    /*Then, we will select the maximum of the array g_c*/ 
    pnga_select_elem(g_c, "max", sresult, index); 
    switch (atype)
    {
      case C_INT:
        *(int*)result = GA_ABS(iresult);
        break;
      case C_DCPL:
      case C_SCPL:
        pnga_error("Ga_step_max_patch_: unavailable for complex datatype.", 
            atype);
        break;
      case C_DBL:
        *(double*)result = GA_ABS(dresult);
        break;
      case C_FLOAT:
        *(float*)result = GA_ABS(fresult);
        break;
      case C_LONG:
        *(long*)result = GA_ABS(lresult);
        break;
      default:
        pnga_error("Ga_step_max_patch_: wrong data type.", atype);
    }
    pnga_destroy(g_c);
  }
  GA_POP_NAME;
  if(local_sync_end)pnga_sync();
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_step_max = pnga_step_max
#endif
void pnga_step_max(Integer g_a, Integer g_b, void *retval)
{
   Integer atype, andim;
   Integer btype, bndim;
   Integer alo[MAXDIM],ahi[MAXDIM];
   Integer blo[MAXDIM],bhi[MAXDIM];
 
    pnga_inquire(g_a,  &atype, &andim, ahi);
    pnga_inquire(g_b,  &btype, &bndim, bhi);
    while(andim){
        alo[andim-1]=1;
        andim--;
        blo[bndim-1]=1;
        bndim--;
    }
    
#if 0
    pnga_step_max_patch(g_a, alo, ahi, g_b, blo, bhi, retval, OP_STEPMAX);
#else
    pnga_step_max_patch(g_a, alo, ahi, g_b, blo, bhi, retval);
#endif
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_step_bound_info = pnga_step_bound_info
#endif
void pnga_step_bound_info(Integer g_xx, Integer g_vv, Integer g_xxll,
                          Integer g_xxuu,  void *boundmin, void *wolfemin,
                          void *boundmax)
{
  Integer xxtype, xxndim;
  Integer vvtype, vvndim;
  Integer xxlltype, xxllndim;
  Integer xxuutype, xxuundim;
  Integer xxlo[MAXDIM],xxhi[MAXDIM];
  Integer vvlo[MAXDIM],vvhi[MAXDIM];
  Integer xxlllo[MAXDIM],xxllhi[MAXDIM];
  Integer xxuulo[MAXDIM],xxuuhi[MAXDIM];

  pnga_inquire(g_xx,  &xxtype, &xxndim, xxhi);
  pnga_inquire(g_vv,  &vvtype, &vvndim, vvhi);
  pnga_inquire(g_xxll,  &xxlltype, &xxllndim, xxllhi);
  pnga_inquire(g_xxuu,  &xxuutype, &xxuundim, xxuuhi);
  while(xxndim){
    xxlo[xxndim-1]=1;
    xxndim--;
    vvlo[vvndim-1]=1;
    vvndim--;
    xxlllo[xxllndim-1]=1;
    xxllndim--;
    xxuulo[xxuundim-1]=1;
    xxuundim--;
  }

  pnga_step_bound_info_patch(g_xx,xxlo,xxhi, g_vv,vvlo,vvhi, g_xxll,xxlllo,xxllhi, g_xxuu,xxuulo,xxuuhi, boundmin, wolfemin, boundmax);
}

