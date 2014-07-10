#if HAVE_CONFIG_H
#   include "config.h"
#endif


/**
 * Symmetrizes matrix A:  A := .5 * (A+A`)
 * diag(A) remains unchanged
 * 
 */

#include "globalp.h"
#include "macdecls.h"
#include "ga-papi.h"
#include "ga-wapi.h"

static void gai_add(
        Integer *lo, Integer *hi, void *a, void *b, DoublePrecision alpha,
        Integer type, Integer nelem, Integer ndim)
{
  Integer i, j, m=0;
  Integer nrow, ncol, indexA=0, indexB=0;
  DoublePrecision *A = (DoublePrecision *)a, *B = (DoublePrecision*)b;
  Integer offset1=1, offset2=1;

  nrow = hi[ndim-2] - lo[ndim-2] + 1;
  ncol = hi[ndim-1] - lo[ndim-1] + 1;
  
  for(i=0; i<ndim-2; i++) {
    offset1 *= hi[i] - lo[i] + 1;
    offset2 *= hi[i] - lo[i] + 1;
  }
  offset1 *= nrow;
  
  for(j=0; j<offset2; ++j,indexA=j,indexB=j,m=0) {
    for(i=0; i<nrow*ncol; i++, indexA += offset1, indexB += offset2) {
      if(indexA >= nelem) indexA = j + ++m*offset2;
      A[indexA] = alpha *(A[indexA] + B[indexB]);
    }
  }
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_symmetrize = pnga_symmetrize
#endif
void pnga_symmetrize(Integer g_a) {
  
  DoublePrecision alpha = 0.5;
  Integer i, me = pnga_nodeid();
  Integer alo[GA_MAX_DIM], ahi[GA_MAX_DIM], lda[GA_MAX_DIM], nelem=1;
  Integer blo[GA_MAX_DIM], bhi[GA_MAX_DIM], ldb[GA_MAX_DIM];
  Integer ndim, dims[GA_MAX_DIM], type;
  Logical have_data;
  Integer g_b; /* temporary global array (b = A') */
  Integer num_blocks_a;
  void *a_ptr=NULL, *b_ptr=NULL;
  int local_sync_begin,local_sync_end;
  char *tempB = "A_transpose";

  local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
  _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
  if(local_sync_begin)pnga_sync();

  GA_PUSH_NAME("ga_symmetrize");
  
  num_blocks_a = pnga_total_blocks(g_a);

  pnga_inquire(g_a, &type, &ndim, dims);

  if (type != C_DBL)
    pnga_error("ga_symmetrize: only implemented for double precision",0);

  if (num_blocks_a < 0) {

    if (dims[ndim-1] != dims[ndim-2]) 
      pnga_error("ga_sym: can only sym square matrix", 0L);

    /* Find the local distribution */
    pnga_distribution(g_a, me, alo, ahi);


    have_data = ahi[0]>0;
    for(i=1; i<ndim; i++) have_data = have_data && ahi[i]>0;

    if(have_data) {
      pnga_access_ptr(g_a, alo, ahi, &a_ptr, lda); 

      for(i=0; i<ndim; i++) nelem *= ahi[i]-alo[i] +1;
      b_ptr = (void *) ga_malloc(nelem, MT_F_DBL, "v");

      for(i=0; i<ndim-2; i++) {bhi[i]=ahi[i]; blo[i]=alo[i]; }

      /* switch rows and cols */
      blo[ndim-1]=alo[ndim-2];
      bhi[ndim-1]=ahi[ndim-2];
      blo[ndim-2]=alo[ndim-1];
      bhi[ndim-2]=ahi[ndim-1];

      for (i=0; i < ndim-1; i++) 
        ldb[i] = bhi[i] - blo[i] + 1; 
      pnga_get(g_a, blo, bhi, b_ptr, ldb);
    }
    pnga_sync(); 

    if(have_data) {
      gai_add(alo, ahi, a_ptr, b_ptr, alpha, type, nelem, ndim);
      pnga_release_update(g_a, alo, ahi);
      ga_free(b_ptr);
    }
  } else {
    /* For block-cyclic data, probably most efficient solution is to
       create duplicate copy, transpose it and add the results together */
    DoublePrecision half = 0.5;
    if (!pnga_duplicate(g_a, &g_b, tempB))
      pnga_error("ga_symmetrize: duplicate failed", 0L);
    pnga_transpose(g_a, g_b);
    pnga_add(&half, g_a, &half, g_b, g_a);
    pnga_destroy(g_b);
  }
  GA_POP_NAME;
  if(local_sync_end)pnga_sync();
}
