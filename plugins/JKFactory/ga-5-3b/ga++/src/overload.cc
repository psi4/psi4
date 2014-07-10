#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "ga++.h"

#ifdef FALSE
#undef FALSE
#endif
#ifdef TRUE
#undef TRUE
#endif
#define FALSE 0
#define TRUE  1

/**
 * More operator overloading stuff (a lot!!) to come.
 */

GA::GlobalArray&
GA::GlobalArray::operator=(const GA::GlobalArray &g_a) {

  if(this != &g_a) { 
    GA_Destroy(mHandle);
    
    mHandle = GA_Duplicate(g_a.mHandle, g_a.inquireName());
    if(!mHandle)  GA_Error((char *)" GA creation failed",0);
    
    GA_Copy(g_a.mHandle, mHandle); 
  }
  return *this;
}

int
GA::GlobalArray::operator==(const GA::GlobalArray &g_a) const {

  long isEqual = TRUE;
  
  int i, type1, type2, ndim1, ndim2, dims1[GA_MAX_DIM], dims2[GA_MAX_DIM];
  int alo[GA_MAX_DIM], ahi[GA_MAX_DIM], blo[GA_MAX_DIM], bhi[GA_MAX_DIM];

  NGA_Inquire(mHandle, &type1, &ndim1, dims1);
  NGA_Inquire(g_a.mHandle, &type2, &ndim2, dims2);
  if(type1 != type2) isEqual = FALSE; // check type
  if(GA_Compare_distr(mHandle, g_a.mHandle)) isEqual = FALSE; 
  NGA_Distribution(mHandle, GA_Nodeid(), alo, ahi);
  NGA_Distribution(g_a.mHandle, GA_Nodeid(), blo, bhi);
  if(ahi[0] != bhi[0]) isEqual = FALSE; // check process owns data?
  
  if(ahi[0] >= 0) { // true => process owns data 
    void *ptr1 = NULL, *ptr2 = NULL;
    int ld1[GA_MAX_DIM];
    int ld2[GA_MAX_DIM];
    int num = 0;

    NGA_Access(mHandle, alo, ahi, &ptr1, ld1);
    NGA_Access(g_a.mHandle, blo, bhi, &ptr2, ld2);

    // number of elements I own.
    for(i=0; i<ndim1; i++) num += ahi[i] - alo[i] + 1;
    
    for(i=0; i<num; i++)
      switch (type1) {
      case C_INT:
	if(((int *)ptr1)[i] != ((int *)ptr2)[i]) isEqual = FALSE;
	break;
      case C_LONG:
	if(((long *)ptr1)[i] != ((long *)ptr2)[i]) isEqual = FALSE;
	break;
      case C_FLOAT:
	if(((float *)ptr1)[i] != ((float *)ptr2)[i]) isEqual = FALSE;
	break;
      case C_DBL:
	if(((double *)ptr1)[i] != ((double *)ptr2)[i]) isEqual = FALSE;
	break;
      case C_DCPL: 
	if(((DoubleComplex *)ptr1)[i].real !=
	   ((DoubleComplex *)ptr2)[i].real) isEqual = FALSE;
	if(((DoubleComplex *)ptr1)[i].imag !=
	   ((DoubleComplex *)ptr2)[i].imag) isEqual = FALSE;
      }
    NGA_Release(mHandle, alo, ahi);
    NGA_Release(g_a.mHandle, blo, bhi);
  }
  
  GA_Lgop(&isEqual, 1, (char *)"*");
  if(isEqual == TRUE) return TRUE;
  else return FALSE;
}

int
GA::GlobalArray::operator!=(const GA::GlobalArray &g_a) const {
  if(*this == g_a) return FALSE;
  return TRUE;
}
