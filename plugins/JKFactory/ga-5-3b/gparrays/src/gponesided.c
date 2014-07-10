#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif

#include "gacommon.h"
#include "typesf2c.h"
#include "gpbase.h"
#include "armci.h"
#include "message.h"
#include "ga-papi.h"
#include "gp-papi.h"
#include "gp-wapi.h"

#define gpm_GetRangeFromMap(p, ndim, plo, phi){                           \
  Integer   _mloc = p* ndim *2;                                           \
            *plo  = (Integer*)_gp_map + _mloc;                            \
            *phi  = *plo + ndim;                                          \
}

/**
 *  Utility arrays used in onesided communication
 */
Integer *_gp_map;
Integer *_gp_proclist;

/**
 *  Initialize utility arrays
 */
void gpi_onesided_init()
{
  Integer nproc;
  nproc = pnga_pgroup_nnodes(pnga_pgroup_get_world());
  _gp_proclist = (Integer*)malloc((size_t)nproc*sizeof(Integer));
  _gp_map = (Integer*)malloc((size_t)(nproc*2*GP_MAX_DIM)*sizeof(Integer));
}

/**
 * Clean utility arrays
 */
void gpi_onesided_clean()
{
  free(_gp_map);
  free(_gp_proclist);
}

/**
 * Get sizes of element in GP array and return them to a local buffer. Also
 * evaluate the total size of the data to be copied and return that in the
 * variable size. intsize is an internal variable that can be used to
 * distinguish between 4 and 8 byte integers
 * @param[in] g_p                pointer array handle
 * @param[in] lo[ndim]           lower corner of pointer array block
 * @param[in] hi[ndim]           upper corner of pointer array block
 * @param[out] buf               buffer that holds element sizes
 * @param[in] ld[ndim-1]         physical dimensions of buffer
 * @param[out] size              total size of requested data
 * @param[in] intsize            parameter to distinguish between 4 and 8
 *                               byte integers
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_get_size = pgp_get_size
#endif

void pgp_get_size(Integer g_p, Integer *lo, Integer *hi,
                  Integer *size, Integer intsize)
{
  Integer handle, ndim, i, nelems;
  Integer block_ld[GP_MAX_DIM-1];
  int *int_ptr;
  int *long_ptr;
  handle = g_p + GP_OFFSET;
  if (!GP[handle].active) {
    pnga_error("gp_get_size: inactive array handle specified", g_p);
  }
  ndim = pnga_ndim(GP[handle].g_ptr_array);
  for (i=0; i<ndim; i++) {
    if (GP[handle].lo[i] > GP[handle].hi[i])
      pnga_error("gp_get_size: illegal block size specified", g_p);
  }

  /* Find total number of elements in block and get strides of requested block */
  ndim = GP[handle].ndim;
  nelems = 1;
  for (i=0; i<ndim; i++) {
    nelems *= (hi[i] - lo[i] + 1);
    if (i<ndim-1) {
      block_ld[i] = (hi[i] - lo[i] + 1);
    }
  }


  *size = 0;
  /* Find total size of elements in block */
  if (intsize == 4) {
    int_ptr = (int*)malloc(nelems*sizeof(int));
    pnga_get(GP[handle].g_size_array, lo, hi, int_ptr, block_ld);
    for (i=0; i<nelems; i++) {
      *size += (Integer)int_ptr[i];
    }
    free(int_ptr);
  } else {
    long_ptr = (int*)malloc(nelems*sizeof(int64_t));
    pnga_get(GP[handle].g_size_array, lo, hi, long_ptr, block_ld);
    for (i=0; i<nelems; i++) {
      *size += (Integer)long_ptr[i];
    }
    free(long_ptr);
  }

}

/**
 * Get data from a GP array and return it to a local buffer. Also return
 * an array of pointers to data in local buffer.
 * @param[in] g_p                pointer array handle
 * @param[in] lo[ndim]           lower corner of pointer array block
 * @param[in] hi[ndim]           upper corner of pointer array block
 * @param[out] buf               buffer that holds data
 * @param[out] buf_ptr           buffer that holds pointers to data
 * @param[in] ld[ndim-1]         physical dimensions of buf_ptr
 * @param[out] buf_size          buffer that holds size of data
 * @param[in] ld_sz[ndim-1]      physical dimensions of buf_size
 * @param[out] size              total size of requested data
 * @param[in] intsize            parameter to distinguish between 4 and 8
 *                               byte integers
 * @param[in] setbuf             controls how information about local
 *                               buffers and element sizes is used.
 *                                  if setbuf = 0, buf is assumed to be
 *                                  a valid pointer and buf_ptr and buf_size
 *                                  are assigned on output. Otherwise,
 *                                  buf is ignored and buf_ptr is assumed to
 *                                  contain pointers to memory locations big
 *                                  enough to hold incoming data. The size of
 *                                  these locations are in buf_size
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_get = pgp_get
#endif

void pgp_get(Integer g_p, Integer *lo, Integer *hi, void *buf,
             void **buf_ptr, Integer *ld, void *buf_size, Integer *ld_sz,
             Integer *size, Integer intsize, Integer setbuf)
{
  Integer handle, ndim, i, j, d, itmp, offset_sz, np;
  Integer idx, offset_d, offset_ptr, offset_rem;
  Integer nelems, index[GP_MAX_DIM];
  Integer block_ld[GP_MAX_DIM], block_ld_loc[GP_MAX_DIM];
  Integer me = (Integer)armci_msg_me();
  void ***src_array, ***dst_array;
  armci_meminfo_t *rem_ptr;
  int rc, bytes;
  armci_giov_t *desc;
  handle = g_p + GP_OFFSET;
  if (!GP[handle].active) {
    pnga_error("gp_get: inactive array handle specified", g_p);
  }
  ndim = pnga_ndim(GP[handle].g_ptr_array);
  for (i=0; i<ndim; i++) {
    if (GP[handle].lo[i] > GP[handle].hi[i])
      pnga_error("gp_get: illegal block size specified", g_p);
  }
  /*bjp
    printf("p[%d] (gp_get) lo[0]: %d hi[0]: %d lo[1]: %d hi[1]: %d\n",
    me,lo[0],hi[0],lo[1],hi[1]);
   */

  if (!setbuf) {
    pnga_get(GP[handle].g_size_array, lo, hi, buf_size, ld_sz);
  }

  /* Get strides of requested block */
  ndim = GP[handle].ndim;
  nelems = 1;
  for (i=0; i<ndim; i++) {
    block_ld[i] = (hi[i] - lo[i] + 1);
    nelems *= block_ld[i];
  }

  /* Based on sizes, construct buf_ptr array */
  idx = 0;
  offset_ptr = 0;

  /*bjp
    printf("p[%d] lo[0]: %d hi[0]: %d lo[1]: %d hi[1]: %d\n",me,lo[0],hi[0],lo[1],hi[1]);
    */
  /*bjp
    printf("p[%d] ld[0]: %d ld_sz[0]: %d\n",me,ld[0],ld_sz[0]);
   */
  while(idx<nelems && !setbuf) {
    /* find local index for idx in the requested block */
    itmp = idx;
    for (j=0; j<ndim-1; j++) {
      index[j] = itmp%block_ld[j];
      itmp = (itmp - index[j])/block_ld[j];
    }
    index[ndim-1] = itmp;
    /* use index to evaluate offset in size array and
       buf_ptr array */
    offset_sz = index[ndim-1];
    offset_d = index[ndim-1];
    for (d=ndim-2; d>=0; d--) {
      offset_sz = offset_sz*ld_sz[d] + index[d];
      offset_d = offset_d*ld[d] + index[d];
    }
    /* evaluate offset in data buffer */
    buf_ptr[offset_d] = (void*)(((char*)buf)+offset_ptr);
    if (intsize == 4) {
      offset_ptr += ((int*)buf_size)[offset_sz];
    } else {
      offset_ptr += ((int64_t*)buf_size)[offset_sz];
    }
    /*bjp
      printf("p[%d] (gp_get) buf_size[%d]: %d\n",me,offset_sz,((int*)buf_size)[offset_sz]);
     */
    idx++;
  }
  *size = offset_ptr;

  /* locate the processors containing some portion of the patch represented by
   * lo and hi and return the results in _gp_map, gp_proclist, and np.
   * _gp_proclist contains a list of processors containing some portion of the
   * patch, _gp_map contains the lower and upper indices of the portion of the
   * patch held by a given processor and np contains the number of processors
   * that contain some portion of the patch.
   */
  pnga_locate_region(GP[handle].g_size_array, lo, hi, _gp_map, _gp_proclist,
      &np);

  /* Loop over processors containing data */
  for (idx=0; idx<np; idx++) {
    Integer p = _gp_proclist[idx];
    Integer *plo, *phi;
    Integer jcnt;
    int *test_ptr;
    gpm_GetRangeFromMap(idx, ndim, &plo, &phi);
    nelems = 1;
    /* Find out how big patch is */
    for (i=0; i<ndim; i++) {
      nelems *= (phi[i]-plo[i] + 1);
      if (i<ndim-1) {
        block_ld_loc[i] = phi[i]-plo[i] + 1;
      }
    }

    /* allocate src and destination arrays */
    src_array = (void***)malloc(sizeof(void*)*nelems);
    dst_array = (void***)malloc(sizeof(void*)*nelems);

    /* Allocate a buffer to hold remote pointers for patch and and
       array of descriptors for GetV operation */
    rem_ptr = (armci_meminfo_t*)malloc((size_t)(nelems)*sizeof(armci_meminfo_t));
    desc = (armci_giov_t*)malloc((size_t)(nelems)*sizeof(armci_giov_t));

    /* Get remote pointers */
    /*bjp
      printf("p[%d] plo[0]: %d phi[0]: %d plo[1]: %d phi[1]: %d\n",
      me,plo[0],phi[0],plo[1],phi[1]);
      */
    pnga_get(GP[handle].g_ptr_array, plo, phi, rem_ptr, block_ld_loc);
    /* Construct descriptors */
    jcnt = 0;
    for (j=0; j<nelems; j++) {
      itmp = j;
      for (d=0; d<ndim-1; d++) {
        index[d] = itmp%block_ld_loc[d];
        itmp = (itmp - index[d])/block_ld_loc[d];
      }
      index[ndim-1] = itmp;
      src_array[j]=(void**)malloc(sizeof(void*));
      dst_array[j]=(void**)malloc(sizeof(void*));

      /* evaluate local offsets */
      offset_rem = index[ndim-1];
      offset_sz = index[ndim-1] + plo[ndim-1] - lo[ndim-1];
      offset_d = index[ndim-1] + plo[ndim-1] - lo[ndim-1];
      for (d=ndim-2; d>=0; d--) {
        offset_rem = offset_rem*block_ld_loc[d] + index[d];
        offset_sz = offset_sz*ld_sz[d] + index[d] + plo[d] - lo[d];
        offset_d = offset_d*ld[d] + index[d] + plo[d] - lo[d];
      }
      /*bjp
        printf("p[%d] j: %d offset_rem: %d offset_sz: %d offset_d: %d\n",
            me, j, offset_rem, offset_sz, offset_d);
            */
      if (intsize == 4) {
        bytes = (int)((int*)buf_size)[offset_sz];
      } else {
        bytes = (int)((int64_t*)buf_size)[offset_sz];
      }
      /*bjp
        printf("p[%d] bytes: %d\n",me,bytes);
       */
      if (bytes > 0) {
        if (rem_ptr[offset_rem].cpid == me) {
          (src_array[j])[0] = ((void*)(rem_ptr[offset_rem].addr));
          /* handle remote and SMP case */
        } else if (pnga_cluster_proc_nodeid(me) ==
            pnga_cluster_proc_nodeid(rem_ptr[offset_rem].cpid)) {
          (src_array[j])[0] = ARMCI_Memat(&rem_ptr[offset_rem],
              sizeof(armci_meminfo_t));
        } else {
          (src_array[j])[0] = (void*)rem_ptr[offset_rem].armci_addr;
        }
        (dst_array[j])[0] = (void*)buf_ptr[offset_d];
#define DBG
        desc[jcnt].src_ptr_array = src_array[j];
        desc[jcnt].dst_ptr_array = dst_array[j];
        desc[jcnt].bytes = bytes;
        desc[jcnt].ptr_array_len = 1;
        /*bjp
          printf("p[%ld] jcnt: %d nelems: %ld index[%ld,%ld] bytes: %d src_ptr: %ld p: %d dst_ptr: %ld\n",
          (long)pnga_nodeid(), jcnt, (long)nelems,
          (long)(index[0]+plo[0]), (long)(index[1]+plo[1]),
          desc[jcnt].bytes, (long)desc[jcnt].src_ptr_array[0], (int)p,
          (long)desc[jcnt].dst_ptr_array[0]);
          */
        jcnt++;
      } else {
        /*bjp
          printf("p[%ld] null pointer at i: %ld j: %ld\n", (long)pnga_nodeid(),
          (long)(index[0]+plo[0]), (long)(index[1]+plo[1]));
         */
      }
    }
    /*bjp
    printf("p[%ld] (gp_get) jcnt: %d p: %d\n",(long)pnga_nodeid(),jcnt,p);
    */
#ifdef XDBG
    for (j=0; j<jcnt; j++) {
      ARMCI_Get(desc[j].src_ptr_array[0], desc[j].dst_ptr_array[0],
          desc[j].bytes, p);
    }
#else
    if (jcnt > 0) {
      rc = ARMCI_GetV(desc, (int)jcnt, (int)p);
      if (rc) pnga_error("ARMCI_GetV failure in gp_get",rc);
    }
#endif
    for (j=0; j<jcnt; j++) {
      free(src_array[j]);
      free(dst_array[j]);
    }
    /* Free temporary buffers */
    free(rem_ptr);
    free(desc);
    free(src_array);
    free(dst_array);
  }
}

/**
 * Put data from a set of local buffers into a GP array.
 * @param[in] g_p                pointer array handle
 * @param[in] lo[ndim]           lower corner of pointer array block
 * @param[in] hi[ndim]           upper corner of pointer array block
 * @param[in] buf_ptr            buffer that holds pointers to local data
 * @param[in] ld[ndim-1]         physical dimensions of buf_ptr
 * @param[in] buf_size           buffer that holds size of local data
 * @param[in] ld_sz[ndim-1]      physical dimensions of buf_size
 * @param[out] size              total size of transmitted data
 * @param[in] checksize          check that sizes in buf_size are OK
 * @param[in] intsize            parameter to distinguish between 4 and 8
 *                               byte integers
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_put = pgp_put
#endif

void pgp_put(Integer g_p, Integer *lo, Integer *hi, void **buf_ptr,
             Integer *ld, void *buf_size, Integer *ld_sz, Integer *size,
             Integer checksize, Integer intsize)
{
  Integer handle, ndim, i, j, d, itmp, offset_sz, np;
  Integer idx, offset_d, offset_ptr, offset_rem;
  Integer nelems, index[GP_MAX_DIM];
  Integer block_ld[GP_MAX_DIM], block_ld_loc[GP_MAX_DIM];
  Integer me = (Integer)armci_msg_me();
  void ***src_array, ***dst_array;
  armci_meminfo_t *rem_ptr;
  int rc, bytes;
  int *tmpsize32;
  int64_t *tmpsize64;
  armci_giov_t *desc;
  handle = g_p + GP_OFFSET;
  if (!GP[handle].active) {
    pnga_error("gp_put: inactive array handle specified", g_p);
  }
  ndim = pnga_ndim(GP[handle].g_ptr_array);
  for (i=0; i<ndim; i++) {
    if (GP[handle].lo[i] > GP[handle].hi[i])
      pnga_error("gp_put: illegal block size specified", g_p);
  }

  /* Get strides of target block */
  ndim = GP[handle].ndim;
  nelems = 1;
  for (i=0; i<ndim; i++) {
    block_ld[i] = (hi[i] - lo[i] + 1);
    nelems *= block_ld[i];
  }

  /* Check size of elements in buf_size against size of elements in GP size
   * array. Throw an error if there is a mismatch */
  if (checksize) {
    if (intsize == 4) {
      tmpsize32 = (int*)malloc(nelems*sizeof(int));
      pnga_get(GP[handle].g_size_array, lo, hi, tmpsize32, block_ld);
    } else {
      tmpsize64 = (int64_t*)malloc(nelems*sizeof(int64_t));
      pnga_get(GP[handle].g_size_array, lo, hi, tmpsize64, block_ld);
    }
    if (intsize == 4) {
      for (i=0; i<nelems; i++) {
        if (tmpsize32[i] != (int)((int*)buf_size)[i])
          pnga_error("gp_put: mismatch in element sizes", i);
      }
    } else {
      for (i=0; i<nelems; i++) {
        if (tmpsize64[i] != (int64_t)((int64_t*)buf_size)[i])
          pnga_error("gp_put: mismatch in element sizes", i);
      }
    }
    if (intsize == 4) {
      free(tmpsize32);
    } else {
      free(tmpsize64);
    }
  }

  idx = 0;
  offset_ptr = 0;

  /* locate the processors containing some portion of the patch represented by
   * lo and hi and return the results in _gp_map, gp_proclist, and np.
   * _gp_proclist contains a list of processors containing some portion of the
   * patch, _gp_map contains the lower and upper indices of the portion of the
   * patch held by a given processor and np contains the number of processors
   * that contain some portion of the patch.
   */
  pnga_locate_region(GP[handle].g_size_array, lo, hi, _gp_map, _gp_proclist,
      &np);
  /* Loop over processors containing data */
  for (idx=0; idx<np; idx++) {
    Integer p = _gp_proclist[idx];
    Integer *plo, *phi;
    Integer jcnt;
    gpm_GetRangeFromMap(idx, ndim, &plo, &phi);
    nelems = 1;
    /* Find out how big patch is */
    for (i=0; i<ndim; i++) {
      nelems *= (phi[i]-plo[i] + 1);
      if (i<ndim-1) {
        block_ld_loc[i] = phi[i]-plo[i] + 1;
      }
    }

    /* allocate src and destination arrays */
    src_array = (void***)malloc(sizeof(void*)*nelems);
    dst_array = (void***)malloc(sizeof(void*)*nelems);

    /* Allocate a buffer to hold remote pointers for patch and and
       array of descriptors for GetV operation */
    rem_ptr = (armci_meminfo_t*)malloc((size_t)(nelems)*sizeof(armci_meminfo_t));
    desc = (armci_giov_t*)malloc((size_t)(nelems)*sizeof(armci_giov_t));

    /* Get remote pointers */
    pnga_get(GP[handle].g_ptr_array, plo, phi, rem_ptr, block_ld_loc);
    /* Construct descriptors */
    jcnt = 0;
    for (j=0; j<nelems; j++) {
      itmp = j;
      for (d=0; d<ndim-1; d++) {
        index[d] = itmp%block_ld_loc[d];
        itmp = (itmp - index[d])/block_ld_loc[d];
      }
      index[ndim-1] = itmp;
      src_array[j]=(void**)malloc(sizeof(void*));
      dst_array[j]=(void**)malloc(sizeof(void*));

      /* evaluate local offsets */
      offset_rem = index[ndim-1];
      offset_sz = index[ndim-1] + plo[ndim-1] - lo[ndim-1];
      offset_d = index[ndim-1] + plo[ndim-1] - lo[ndim-1];
      for (d=ndim-2; d>=0; d--) {
        offset_rem = offset_rem*block_ld_loc[d] + index[d];
        offset_sz = offset_sz*ld_sz[d] + index[d] + plo[d] - lo[d];
        offset_d = offset_d*ld[d] + index[d] + plo[d] - lo[d];
      }
      if (intsize == 4) {
        bytes = (int)((int*)buf_size)[offset_sz];
      } else {
        bytes = (int)((int64_t*)buf_size)[offset_sz];
      }
      if (bytes > 0) {
        if (rem_ptr[offset_rem].cpid == me) {
          (dst_array[j])[0] = ((void*)(rem_ptr[offset_rem].addr));
          /* handle remote and SMP case */
        } else if (pnga_cluster_proc_nodeid(me) ==
            pnga_cluster_proc_nodeid(rem_ptr[offset_rem].cpid)) {
          (dst_array[j])[0] = ARMCI_Memat(&rem_ptr[offset_rem],
              sizeof(armci_meminfo_t));
        } else {
          (dst_array[j])[0] = (void*)rem_ptr[offset_rem].armci_addr;
        }
        (src_array[j])[0] = (void*)buf_ptr[offset_d];
#define DBG
        desc[jcnt].src_ptr_array = src_array[j];
        desc[jcnt].dst_ptr_array = dst_array[j];
        desc[jcnt].bytes = bytes;
        desc[jcnt].ptr_array_len = 1;
        jcnt++;
      } else {
        /*bjp
          printf("p[%ld] null pointer at i: %ld j: %ld\n", (long)pnga_nodeid(),
          (long)(index[0]+plo[0]), (long)(index[1]+plo[1]));
         */
      }
    }
#ifdef XDBG
    for (j=0; j<jcnt; j++) {
      ARMCI_Put(desc[j].src_ptr_array[0], desc[j].dst_ptr_array[0],
          desc[j].bytes, p);
    }
#else
    if (jcnt > 0) {
      rc = ARMCI_PutV(desc, (int)jcnt, (int)p);
      if (rc) pnga_error("ARMCI_PutV failure in gp_put",rc);
    }
#endif
    for (j=0; j<jcnt; j++) {
      free(src_array[j]);
      free(dst_array[j]);
    }
    /* Free temporary buffers */
    free(rem_ptr);
    free(desc);
  }
  free(src_array);
  free(dst_array);
}

/**
 * Return a pointer to a locally held data element, along with the size of the
 * data element.
 * @param g_p[in]        pointer array handle
 * @param subscript[in]  array containing indices of element
 * @param  ptr[out]      pointer to contents of array element
 * @param size[out]      size of array element (in bytes)
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_access_element = pgp_access_element
#endif

void pgp_access_element(Integer g_p, Integer *subscript, void *ptr, Integer *size)
{
  armci_meminfo_t *gp_ptr;
  Integer handle, ld[GP_MAX_DIM-1], i;
  Integer one=1;
  Integer *lo, *hi;

  handle = g_p + GP_OFFSET;
  lo = GP[handle].lo;
  hi = GP[handle].hi;

  /* check to make sure that element is located in local block
   * of GP array */
  for (i=0; i<GP[handle].ndim; i++) {
    if (subscript[i]<lo[i] || subscript[i]>hi[i]) {
      pnga_error("gp_access_element: subscript out of bounds", i);
    }
    if (i<GP[handle].ndim-1) {
      ld[i] = one;
    }
  }
  pnga_get(GP[handle].g_size_array, subscript, subscript, size, ld);
  pnga_access_ptr(GP[handle].g_ptr_array,subscript,subscript, &gp_ptr, ld);
  *(char**)ptr = (char*)(*gp_ptr).addr;
}

/**
 * Release access to a data element
 * @param g_p[in]        pointer array handle
 * @param subscript[in]  array containing indices of element
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_release_element = pgp_release_element
#endif

void pgp_release_element(Integer g_p, Integer *subscript)
{
  Integer handle = g_p + GP_OFFSET;
  pnga_release(GP[handle].g_ptr_array, subscript, subscript);
}

/**
 * Release access and update contents of a data element
 * @param g_p[in]        pointer array handle
 * @param subscript[in]  array containing indices of element
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_release_update_element = pgp_release_update_element
#endif

void pgp_release_update_element(Integer g_p, Integer *subscript)
{
  Integer handle = g_p + GP_OFFSET;
  pnga_release_update(GP[handle].g_ptr_array, subscript, subscript);
}

/**
 * Determine the total size of a random set of elements from a GP array
 * @param g_p[in]        pointer array handle
 * @param nv[in]         number of elements being requested
 * @param subscript[in]  array containing element indices
 * @param size[out]      total size (in bytes) of requested data
 * @param[in] intsize    parameter to distinguish between 4 and 8
 *                       byte integers
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_gather_size = pgp_gather_size
#endif

void pgp_gather_size(Integer g_p, Integer nv, Integer *subscript, Integer *size,
                     Integer intsize)
{
  void *size_buf;
  Integer i, asize, handle;

  handle = g_p + GP_OFFSET;
  if (intsize == 4) {
    size_buf = (void*)malloc((int)nv*sizeof(int));
    /* BJP initialize to value other than 0 
    for (i=0; i<nv; i++) {
      ((int*)size_buf)[i] = -1;
    }
    */
  } else {
    size_buf = (void*)malloc((int)nv*sizeof(int64_t));
  }
  pnga_gather(GP[handle].g_size_array, size_buf, subscript, 0, nv);

  /* sum up all sizes */
  asize = 0;
  if (intsize == 4) {
    for (i=0; i<nv; i++) {
      asize += (Integer)((int*)size_buf)[i];
    /* BJP
    printf("p[%d] Test1 size: %d\n",pnga_nodeid(),((int*)size_buf)[i]);
    */
    }
  } else {
    for (i=0; i<nv; i++) {
      asize += (Integer)((int64_t*)size_buf)[i];
    }
  }
  free(size_buf);
  *size = asize;
}

/**
 * Gather a list of random elements from a GP array and store them in local
 * buffers
 * @param g_p[in]        pointer array handle
 * @param nv[in]         number of elements being requested
 * @param subscript[in]  array containing element indices
 * @param buf[in]        pointer to local buffer that will contain results
 * @param buf_ptr[out]   pointers to local copies of elements
 * @param buf_size[out]  array with element sizes
 * @param size[out]      total size (in bytes) of requested data
 * @param[in] intsize    parameter to distinguish between 4 and 8
 *                       byte integers
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_gather = pgp_gather
#endif

void pgp_gather(Integer g_p, Integer nv, Integer *subscript, void *buf,
                void **buf_ptr, void *buf_size, Integer *size, Integer intsize,
                Integer setbuf)
{
  Integer handle;
  Integer *header, *list, *nelems;
  Integer me, iproc, nproc, i, j, idx;
  void *l_ptr;
  armci_meminfo_t *info_buf;
  armci_giov_t *desc;
  void ***src_array, ***dst_array;
  int rc, bytes;

  handle = g_p + GP_OFFSET;
  nproc = pnga_nnodes();
  me = pnga_nodeid();

  info_buf = (armci_meminfo_t*)malloc((int)nv*sizeof(armci_meminfo_t)); 
  pnga_gather(GP[handle].g_ptr_array, info_buf, subscript, 0, nv);
  if (!setbuf) {
    pnga_gather(GP[handle].g_size_array, buf_size, subscript, 0, nv);
  }

  /* create link list arrays and other utility arrays
   * header[iproc]: first element in list for processor iproc
   * list[i]: pointer to next element in list
   * nelems[iproc]: total number of elements on processor iproc
   */
  header = (Integer*)malloc((int)nproc*sizeof(Integer));
  list = (Integer*)malloc((int)nv*sizeof(Integer));
  nelems = (Integer*)malloc((int)nproc*sizeof(Integer));

  for (i=0; i<nproc; i++) {
    nelems[i] = 0;
    header[i] = -1;
  }
  for (i=0; i<nv; i++) {
    list[i] = 0;
  }

  l_ptr = buf;
  if (intsize == 4) {
    for (i=0; i<nv; i++) {
      idx = (Integer)info_buf[i].cpid;
      if (!setbuf) {
        buf_ptr[i] = l_ptr;
      }
      l_ptr = (void*)((char*)l_ptr+(int)((int*)buf_size)[i]);
      nelems[idx]++;
      j = header[idx];
      header[idx] = i;
      list[i] = j;
    }
  } else {
    for (i=0; i<nv; i++) {
      idx = (Integer)info_buf[i].cpid;
      if (!setbuf) {
        buf_ptr[i] = l_ptr;
      }
      l_ptr = (void*)((char*)l_ptr+(int)((int64_t*)buf_size)[i]);
      nelems[idx]++;
      j = header[idx];
      header[idx] = i;
      list[i] = j;
    }
  }
  
  /* scan through linked list and get data from each processor */
  for (iproc=0; iproc<nproc; iproc++) {
    if (nelems[iproc] > 0) {
      idx = header[iproc];
      /* allocate descriptor array for this vector call */
      desc = (armci_giov_t*)malloc((int)nelems[iproc]*sizeof(armci_giov_t));
      src_array = (void***)malloc((int)nelems[iproc]*sizeof(void**));
      dst_array = (void***)malloc((int)nelems[iproc]*sizeof(void**));
      j = 0;
      while (idx > -1) {
        if (intsize == 4) {
          bytes = (int)((int*)buf_size)[idx];
        } else {
          bytes = (int)((int64_t*)buf_size)[idx];
        }
        if (bytes>0) {
          src_array[j] = (void**)malloc(sizeof(void*));
          dst_array[j] = (void**)malloc(sizeof(void*));
          if (iproc == me) {
            (src_array[j])[0] = ((void*)(info_buf[idx].addr));
          } else if (pnga_cluster_proc_nodeid(me) ==
              pnga_cluster_proc_nodeid(iproc)) {
            (src_array[j])[0] = ARMCI_Memat(&info_buf[idx],sizeof(armci_meminfo_t));
          } else { 
            (src_array[j])[0] = (void*)(info_buf[idx].armci_addr);
          }
          (dst_array[j])[0] = (void*)(buf_ptr[idx]);
          if (intsize == 4) {
            desc[j].bytes = (int)((int*)buf_size)[idx];
          } else {
            desc[j].bytes = (int)((int64_t*)buf_size)[idx];
          }
          desc[j].src_ptr_array = src_array[j];
          desc[j].dst_ptr_array = dst_array[j];
          desc[j].ptr_array_len = 1;

          j++;
        }
        idx = list[idx];
      }

      /* gather data from remote locations */
#ifdef XDBG
    for (idx=0; idx<j; idx++) {
      ARMCI_Get(desc[idx].src_ptr_array[0], desc[idx].dst_ptr_array[0],
          desc[idx].bytes, iproc);
    }
#else
      if (j > 0) {
        rc = ARMCI_GetV(desc, (int)j, (int)iproc);
        if (rc) pnga_error("ARMCI_GetV failure in gp_gather",rc);
      }
#endif

      /* free arrays */
      for (j=0; j<nelems[iproc]; j++) {
        free(src_array[j]);
        free(dst_array[j]);
      }
      free(src_array);
      free(dst_array);
      free(desc);
    }
  }

  /* free remaining temporary arrays */
  free(info_buf);
  free(header);
  free(list);
  free(nelems);
  
  *size = 0;
  for (i=0; i<nv; i++) {
    *size += ((Integer*)buf_size)[i];
  }
}

/**
 * Scatter a list of random elements from local buffers and store them in
 * a GP array
 * @param g_p[in]        pointer array handle
 * @param nv[in]         number of elements being requested
 * @param subscript[in]  array containing element indices
 * @param buf_ptr[out]   pointers to local copies of elements
 * @param buf_size[out]  array with element sizes
 * @param size[out]      total size (in bytes) of data to be moved
 * @param checksize[in]  check size of data
 * @param[in] intsize    parameter to distinguish between 4 and 8
 *                       byte integers
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_scatter = pgp_scatter
#endif

void pgp_scatter(Integer g_p, Integer nv, Integer *subscript, void **buf_ptr,
                 void *buf_size, Integer *size, Integer checksize, Integer intsize)
{
  Integer handle;
  Integer *header, *list, *nelems;
  Integer me, iproc, nproc, i, j, idx;
  void *l_ptr;
  armci_meminfo_t *info_buf;
  armci_giov_t *desc;
  void ***src_array, ***dst_array;
  int rc, bytes;
  int *tmpsize32;
  int64_t *tmpsize64;

  handle = g_p + GP_OFFSET;
  nproc = pnga_nnodes();
  me = pnga_nodeid();

  /* Check size of elements in buf_size against size of elements in GP size
   * array. Throw an error if there is a mismatch */
  if (checksize) {
    if (intsize == 4) {
      tmpsize32 = (int*)malloc(nv*sizeof(int));
      pnga_gather(GP[handle].g_size_array, tmpsize32, subscript, 0, nv);
    } else {
      tmpsize64 = (int64_t*)malloc(nv*sizeof(int64_t));
      pnga_gather(GP[handle].g_size_array, tmpsize64, subscript, 0, nv);
    }
    if (intsize == 4) {
      for (i=0; i<nv; i++) {
        if (tmpsize32[i] != (int)((int*)buf_size)[i]) {
           printf("p[%ld] tmpsize[%ld]: %d buf_size[%ld]: %ld\n",
                   (long)me,(long)i,tmpsize32[i],
                   (long)i,(long)((int*)buf_size)[i]);
       /*   pnga_error("gp_scatter: mismatch in element sizes", i); */
        }
      }
    } else {
      for (i=0; i<nv; i++) {
        if (tmpsize64[i] != (int64_t)((int64_t*)buf_size)[i]) {
        /*  pnga_error("gp_scatter: mismatch in element sizes", i); */
        }
      }
    }
    if (intsize == 4) {
      free(tmpsize32);
    } else {
      free(tmpsize64);
    }
  }
  info_buf = (armci_meminfo_t*)malloc((int)nv*sizeof(armci_meminfo_t)); 
  pnga_gather(GP[handle].g_ptr_array, info_buf, subscript, 0, nv);

  /* create link list arrays and other utility arrays
   * header[iproc]: first element in list for processor iproc
   * list[i]: pointer to next element in list
   * nelems[iproc]: total number of elements on processor iproc
   */
  header = (Integer*)malloc((int)nproc*sizeof(Integer));
  list = (Integer*)malloc((int)nv*sizeof(Integer));
  nelems = (Integer*)malloc((int)nproc*sizeof(Integer));

  for (i=0; i<nproc; i++) {
    nelems[i] = 0;
    header[i] = -1;
  }
  for (i=0; i<nv; i++) {
    list[i] = 0;
  }

  for (i=0; i<nv; i++) {
    idx = (Integer)info_buf[i].cpid;
    nelems[idx]++;
    j = header[idx];
    header[idx] = i;
    list[i] = j;
  }
  
  /* scan through linked list and get data from each processor */
  for (iproc=0; iproc<nproc; iproc++) {
    if (nelems[iproc] > 0) {
      idx = header[iproc];
      /* allocate descriptor array for this vector call */
      desc = (armci_giov_t*)malloc((int)nelems[iproc]*sizeof(armci_giov_t));
      src_array = (void***)malloc((int)nelems[iproc]*sizeof(void**));
      dst_array = (void***)malloc((int)nelems[iproc]*sizeof(void**));
      j = 0;
      while (idx > -1) {
        if (intsize == 4) {
          bytes = (int)((int*)buf_size)[idx];
        } else {
          bytes = (int)((int64_t*)buf_size)[idx];
        }
        if (bytes>0) {
          src_array[j] = (void**)malloc(sizeof(void*));
          dst_array[j] = (void**)malloc(sizeof(void*));
          if (iproc == me) {
            (dst_array[j])[0] = ((void*)(info_buf[idx].addr));
          } else if (pnga_cluster_proc_nodeid(me) ==
              pnga_cluster_proc_nodeid(iproc)) {
            (dst_array[j])[0] = ARMCI_Memat(&info_buf[idx],sizeof(armci_meminfo_t));
          } else { 
            (dst_array[j])[0] = (void*)(info_buf[idx].armci_addr);
          }
          (src_array[j])[0] = (void*)(buf_ptr[idx]);
          if (intsize == 4) {
            desc[j].bytes = (int)((int*)buf_size)[idx];
          } else {
            desc[j].bytes = (int)((int64_t*)buf_size)[idx];
          }
          desc[j].src_ptr_array = src_array[j];
          desc[j].dst_ptr_array = dst_array[j];
          desc[j].ptr_array_len = 1;

          j++;
        }
        idx = list[idx];
      }

      /* gather data from remote locations */
#ifdef XDBG
    for (idx=0; idx<j; idx++) {
      ARMCI_Put(desc[idx].src_ptr_array[0], desc[idx].dst_ptr_array[0],
          desc[idx].bytes, iproc);
    }
#else
      if (j > 0) {
        rc = ARMCI_PutV(desc, (int)j, (int)iproc);
        if (rc) pnga_error("ARMCI_PutV failure in gp_gather",rc);
      }
#endif

      /* free arrays */
      for (j=0; j<nelems[iproc]; j++) {
        free(src_array[j]);
        free(dst_array[j]);
      }
      free(src_array);
      free(dst_array);
      free(desc);
    }
  }

  /* free remaining temporary arrays */
  free(info_buf);
  free(header);
  free(list);
  free(nelems);
  
  *size = 0;
  for (i=0; i<nv; i++) {
    *size += (Integer)((Integer*)buf_size)[i];
  }
}
