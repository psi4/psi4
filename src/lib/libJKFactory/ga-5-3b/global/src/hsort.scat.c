#if HAVE_CONFIG_H
#   include "config.h"
#endif

/**************************************************************************\
 
 Sort routines for scatter and gather.
 scatter requires sorting of index and value arrays.
 gather requires sorting of index arrays only.

\**************************************************************************/

  
#include "typesf2c.h"
#include "macommon.h"
#include "ga-papi.h"
#include "ga-wapi.h"

#define GT(a,b) (*(a) > *(b))
#define GE(a,b) (*(a) >= *(b))


#define INDEX_SORT(base,pn,SWAP){\
  unsigned long gap, g;\
  Integer *p, *q, *hi, n=*pn;\
  Integer  *base0=base - 1;\
\
  gap = n >>1;\
  hi = base0 + gap + gap;\
  if (n & 1) hi ++;\
\
  for ( ; gap != 1; gap--) {\
    for (p = base0 + (g = gap) ; (q = p + g) <= hi ; p = q) {\
      g += g;\
      if (q != hi && GT(q+1, q)) {\
        q++;\
        g++;\
      }\
      if (GE(p,q)) break;\
      SWAP(p , q);\
    }\
  }\
\
  for ( ; hi != base ; hi--) {\
    p = base;\
    for (g = 1 ; (q = p + g) <= hi ; p = q) {\
      g += g;\
      if (q != hi && GT(q+1,q)) {\
        q++;\
        g++;\
      }\
      if (GE(p,q)) break;\
      SWAP(p, q);\
    }\
    SWAP(base, hi);\
  }\
}




static void ga_sort_scat_dcpl_(pn, v, i, j, base)
     Integer *pn;
     DoubleComplex *v;
     Integer *i;
     Integer *j;
     Integer *base;
{

  if (*pn < 2) return;

#  undef SWAP  
#  define SWAP(a,b) { \
    Integer ltmp; \
    DoubleComplex dtmp; \
    long ia = a - base; \
    long ib = b - base; \
    ltmp=*a; *a=*b; *b=ltmp; \
    dtmp=v[ia]; v[ia]=v[ib]; v[ib]=dtmp; \
    ltmp=i[ia]; i[ia]=i[ib]; i[ib]=ltmp; \
    ltmp=j[ia]; j[ia]=j[ib]; j[ib]=ltmp; \
  }
  INDEX_SORT(base,pn,SWAP);
}

static void ga_sort_scat_scpl_(pn, v, i, j, base)
     Integer *pn;
     SingleComplex *v;
     Integer *i;
     Integer *j;
     Integer *base;
{

  if (*pn < 2) return;

#  undef SWAP  
#  define SWAP(a,b) { \
    Integer ltmp; \
    SingleComplex dtmp; \
    long ia = a - base; \
    long ib = b - base; \
    ltmp=*a; *a=*b; *b=ltmp; \
    dtmp=v[ia]; v[ia]=v[ib]; v[ib]=dtmp; \
    ltmp=i[ia]; i[ia]=i[ib]; i[ib]=ltmp; \
    ltmp=j[ia]; j[ia]=j[ib]; j[ib]=ltmp; \
  }
  INDEX_SORT(base,pn,SWAP);
}

void ga_sort_permutation(pn, index, base)
     Integer *pn;
     Integer *index;
     Integer *base;
{
  if (*pn < 2) return;
#  undef SWAP  
#  define SWAP(a,b) { \
    Integer ltmp; \
    Integer itmp;\
    long ia = a - base; \
    long ib = b - base; \
    ltmp=*a; *a=*b; *b=ltmp; \
    itmp=index[ia]; index[ia]=index[ib]; index[ib] = itmp;\
   }
  INDEX_SORT(base,pn,SWAP);
}

     


static void ga_sort_scat_dbl_(pn, v, i, j, base)
     Integer *pn;
     DoublePrecision *v;
     Integer *i;
     Integer *j;
     Integer *base;
{
  
  if (*pn < 2) return;

#  undef SWAP  
#  define SWAP(a,b) { \
    Integer ltmp; \
    DoublePrecision dtmp; \
    long ia = a - base; \
    long ib = b - base; \
    ltmp=*a; *a=*b; *b=ltmp; \
    dtmp=v[ia]; v[ia]=v[ib]; v[ib]=dtmp; \
    ltmp=i[ia]; i[ia]=i[ib]; i[ib]=ltmp; \
    ltmp=j[ia]; j[ia]=j[ib]; j[ib]=ltmp; \
  }
  INDEX_SORT(base,pn,SWAP);
}


static void ga_sort_scat_int_(pn, v, i, j, base)
     Integer *pn;
     int *v;
     Integer *i;
     Integer *j;
     Integer *base;
{

  if (*pn < 2) return;

#  undef SWAP  
#  define SWAP(a,b) { \
    int ltmp; \
    int dtmp; \
    long ia = a - base; \
    long ib = b - base; \
    ltmp=*a; *a=*b; *b=ltmp; \
    dtmp=v[ia]; v[ia]=v[ib]; v[ib]=dtmp; \
    ltmp=i[ia]; i[ia]=i[ib]; i[ib]=ltmp; \
    ltmp=j[ia]; j[ia]=j[ib]; j[ib]=ltmp; \
  }
  INDEX_SORT(base,pn,SWAP);
}




static void ga_sort_scat_long_(pn, v, i, j, base)
     Integer *pn;
     long *v;
     Integer *i;
     Integer *j;
     Integer *base;
{
 
  if (*pn < 2) return;
 
#  undef SWAP
#  define SWAP(a,b) { \
    long ltmp; \
    long dtmp; \
    long ia = a - base; \
    long ib = b - base; \
    ltmp=*a; *a=*b; *b=ltmp; \
    dtmp=v[ia]; v[ia]=v[ib]; v[ib]=dtmp; \
    ltmp=i[ia]; i[ia]=i[ib]; i[ib]=ltmp; \
    ltmp=j[ia]; j[ia]=j[ib]; j[ib]=ltmp; \
  }
  INDEX_SORT(base,pn,SWAP);
}

static void ga_sort_scat_flt_(pn, v, i, j, base)
     Integer *pn;
     float   *v;
     Integer *i;
     Integer *j;
     Integer *base;
{
 
  if (*pn < 2) return;
 
#  undef SWAP
#  define SWAP(a,b) { \
    Integer ltmp; \
    float dtmp; \
    int ia = a - base; \
    int ib = b - base; \
    ltmp=*a; *a=*b; *b=ltmp; \
    dtmp=v[ia]; v[ia]=v[ib]; v[ib]=dtmp; \
    ltmp=i[ia]; i[ia]=i[ib]; i[ib]=ltmp; \
    ltmp=j[ia]; j[ia]=j[ib]; j[ib]=ltmp; \
  }
  INDEX_SORT(base,pn,SWAP);
}                                   

void ga_sort_scat(pn, v, i, j, base, type)
     Integer *pn;
     void    *v;
     Integer *i;
     Integer *j;
     Integer *base;
     Integer type;
{ 
   switch (type){
     case C_DBL:  ga_sort_scat_dbl_(pn, (double*)v, i,j,base);break;
     case C_DCPL: ga_sort_scat_dcpl_(pn, (DoubleComplex*)v, i,j,base); break;
     case C_SCPL: ga_sort_scat_scpl_(pn, (SingleComplex*)v, i,j,base); break;
     case C_INT:  ga_sort_scat_int_(pn, (int*)v, i, j, base); break;
     case C_FLOAT:  ga_sort_scat_flt_(pn, (float*)v, i, j, base); break; 
     case C_LONG:  ga_sort_scat_long_(pn, (long*)v, i, j, base); break;
     default:        pnga_error("ERROR:ga_sort_scat: wrong type",type);
   } 
}


void ga_sort_gath(pn, i, j, base)
     Integer *pn;
     Integer *i;
     Integer *j;
     Integer *base;
{

  if (*pn < 2) return;
  
#  undef SWAP  
#  define SWAP(a,b) { \
    Integer ltmp; \
    long ia = a - base; \
    long ib = b - base; \
    ltmp=*a; *a=*b; *b=ltmp; \
    ltmp=i[ia]; i[ia]=i[ib]; i[ib]=ltmp; \
    ltmp=j[ia]; j[ia]=j[ib]; j[ib]=ltmp; \
  }
  INDEX_SORT(base,pn,SWAP);
}


void gai_hsort(Integer *list, int num)
{
  if(num<2) return;
# undef SWAP
# define SWAP(a,b) { Integer ltmp; ltmp=*a; *a=*b; *b=ltmp;}
  INDEX_SORT(list,&num,SWAP);
}
  
  

