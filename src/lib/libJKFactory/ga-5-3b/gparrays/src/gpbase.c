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
#include "ga-papi.h"
#include "gp-papi.h"
#include "gp-wapi.h"
#include "gpbase.h"
#include "armci.h"

extern void gpi_onesided_init();
extern void gpi_onesided_clean();

gp_array_t *GP;
int GP_pointer_type;

/**
 *  Initialize internal library structures for Global Pointer Arrays
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_initialize = pgp_initialize
#endif

void pgp_initialize()
{
  Integer i;
  GP = (gp_array_t*)malloc(sizeof(gp_array_t)*GP_MAX_ARRAYS);
  GP_pointer_type = pnga_register_type(sizeof(armci_meminfo_t));
  if (!GP) {
    pnga_error("gp_initialize: malloc GP failed",0);
  }
  for (i=0; i<GP_MAX_ARRAYS; i++) {
    GP[i].active = 0;
  }
  gpi_onesided_init();
}

/**
 *  Deallocate all arrays and clean up library
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_terminate = pgp_terminate
#endif

void pgp_terminate()
{
  Integer i;
  for (i=0; i<GP_MAX_ARRAYS; i++) {
    if (GP[i].active) {
      pnga_destroy(GP[i].g_size_array);
      pnga_destroy(GP[i].g_ptr_array);
      GP[i].active = 0;
    }
  }
  pnga_deregister_type(GP_pointer_type);
  gpi_onesided_clean();
}

/**
 *  Special malloc() for Global Pointer Arrays
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_malloc = pgp_malloc
#endif

void* pgp_malloc(size_t size)
{
  armci_meminfo_t meminfo, *meminfo_ptr=NULL;
    size_t meminfo_sz = sizeof(armci_meminfo_t);
    
    ARMCI_Memget(size+meminfo_sz, &meminfo, 0);
    /*bjp
    printf("%d: GP_Malloc: ptr = %p\n", pnga_nodeid(), meminfo.addr);fflush(stdout);
    */

    /* store the meminfo handle at the beginning of segment */
    memcpy( meminfo.addr, &meminfo, meminfo_sz);
    meminfo_ptr = (armci_meminfo_t*)meminfo.addr;

    /* update the meminfo structure */
    meminfo_ptr->armci_addr = ((char*)meminfo_ptr->armci_addr) + meminfo_sz;
    meminfo_ptr->addr       = ((char*)meminfo_ptr->addr) + meminfo_sz;
    meminfo_ptr->size      -= meminfo_sz;
    /*bjp
    printf("p[%d]: armci_addr = %ld\n", pnga_nodeid(), (long)meminfo_ptr->armci_addr);
    */
    
    return meminfo_ptr->addr;
}

/**
 *  Special free() for Global Pointer Arrays
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_free = pgp_free
#endif

void pgp_free(void* ptr)
{
    armci_meminfo_t meminfo;
    size_t meminfo_sz = sizeof(armci_meminfo_t);
    
    if(!ptr) pnga_error("gp_free: Invalid pointer",0);
    
    memcpy( &meminfo, ((char*)ptr)-meminfo_sz, meminfo_sz);
    
    /* update the meminfo structure */
    meminfo.armci_addr = ((char*)meminfo.armci_addr) - meminfo_sz;
    meminfo.addr       = ((char*)meminfo.addr) - meminfo_sz;
    meminfo.size       += meminfo_sz;
    ARMCI_Memctl(&meminfo);
}

/**
 *  Create a handle for a GP array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_create_handle = pgp_create_handle
#endif

Integer pgp_create_handle()
{
  Integer i, handle=-GP_OFFSET-1;
  for (i=0; i<GP_MAX_ARRAYS; i++) {
    if (!GP[i].active) {
      handle = i-GP_OFFSET;
      GP[i].g_size_array = pnga_create_handle();
      GP[i].g_ptr_array = pnga_create_handle();
      GP[i].active = 1;
      GP[i].ndim = -1;
      break;
    }
  }
  return handle;
}

/**
 *  Set array dimensions
 *  @param[in] g_p         pointer array handle
 *  @param[in] ndim        dimension of pointer array
 *  @param[in] dims[ndim]  dimension of array axes
 *  @param[in] intsize     size of integers in calling program
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_set_dimensions = pgp_set_dimensions
#endif

void pgp_set_dimensions(Integer g_p, Integer ndim, Integer *dims,
                        Integer intsize)
{
  Integer handle, i, type;
  handle = g_p + GP_OFFSET;

  /* Do some basic checks on parameters */
  if (!GP[handle].active) {
    pnga_error("gp_set_dimensions: Global Pointer handle is not active", 0);
  }
  if (ndim < 0 || ndim > GP_MAX_DIM) {
    pnga_error("gp_set_dimensions: dimension is not valid", ndim);
  }
  for (i=0; i<ndim; i++) {
    if (dims[i] < 0) {
      pnga_error("gp_set_dimensions: invalid dimension found", dims[i]);
    }
  }

  if (intsize == 4) {
    type = C_INT;
  } else {
    type = MT_F_INT;
  }
  pnga_set_data(GP[handle].g_size_array, ndim, dims, type);
  type = GP_pointer_type;
  pnga_set_data(GP[handle].g_ptr_array, ndim, dims, type);
  GP[handle].ndim = ndim;
  for (i=0; i<ndim; i++) {
    GP[handle].dims[i] = dims[i];
  }
}

/**
 *  Determine decomposition of array accross processors
 *  @param[in] g_p          pointer array handle
 *  @param[in] mapc         array giving first index of each
 *                          block for each axis
 *  @param[in] nblock[ndim] number of blocks along each dimension
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_set_irreg_distr = pgp_set_irreg_distr
#endif

void pgp_set_irreg_distr(Integer g_p, Integer *mapc, Integer *nblock)
{
  Integer handle, i, ichk, ndim;
  handle = g_p + GP_OFFSET;

  /* Do some basic checks on parameters */
  if (!GP[handle].active) {
    pnga_error("gp_set_irreg_distr: Global Pointer handle is not active", 0);
  }
  ndim = GP[handle].ndim;
  for (i=0; i<ndim; i++) {
    if (GP[handle].dims[i]<(Integer)nblock[i]) {
      pnga_error("gp_set_irreg_distr: number of blocks <= corresponding dimension", i);
    }
  }

  pnga_set_irreg_distr(GP[handle].g_size_array, mapc, nblock);
  pnga_set_irreg_distr(GP[handle].g_ptr_array, mapc, nblock);
}

/**
 * Return the dimension of the pointer array. This is mostly useful for setting
 * up the C interface.
 * @param[in] g_p         pointer array handle
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_get_dimension = pgp_get_dimension
#endif

Integer pgp_get_dimension(Integer g_p)
{
  Integer handle;
  handle = g_p + GP_OFFSET;
  return (Integer)GP[handle].ndim;
}

/**
 *  Set chunk array dimensions. This determines the minimum dimension of a
 *  local block of data
 *  @param[in] g_p         pointer array handle
 *  @param[in] chunk[ndim] minimum dimensions of array blocks
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_set_chunk = pgp_set_chunk
#endif

void pgp_set_chunk(Integer g_p, Integer *chunk)
{
  Integer handle;
  handle = g_p + GP_OFFSET;
  pnga_set_chunk(GP[handle].g_size_array, chunk);
  pnga_set_chunk(GP[handle].g_ptr_array, chunk);
}

/**
 *  Allocate memory for a Global Pointer array.
 *  @param[in] g_p         pointer array handle
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_allocate = pgp_allocate
#endif

logical pgp_allocate(Integer g_p)
{
  logical status;
  Integer handle, me, i;
  handle = g_p + GP_OFFSET;
  status = pnga_allocate(GP[handle].g_size_array);
  status = status && pnga_allocate(GP[handle].g_ptr_array);
  if (!status) {
    pnga_error("gp_allocate: unable to allocate GP array", 0);
  }
  pnga_zero(GP[handle].g_size_array);
  pnga_zero(GP[handle].g_ptr_array);
  me = pnga_nodeid();
  pnga_distribution(GP[handle].g_ptr_array, me, GP[handle].lo,
          GP[handle].hi);
  GP[handle].active = 1;
  for (i=0; i<GP[handle].ndim-1; i++) {
    GP[handle].ld[i] = GP[handle].hi[i] - GP[handle].lo[i] + 1;
  }
  return status;
}

/**
 *  Destroy Global Pointer array and free memory for reuse.
 *  @param[in] g_p         pointer array handle
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_destroy = pgp_destroy
#endif

logical pgp_destroy(Integer g_p)
{
  logical status;
  Integer handle;
  handle = g_p + GP_OFFSET;
  status = pnga_destroy(GP[handle].g_size_array);
  status = status && pnga_destroy(GP[handle].g_ptr_array);
  if (!status) {
    pnga_error("gp_destroy: unable to destroy GP array", 0);
  }
  GP[handle].active = 0;
  return status;
}

/**
 *  Shell that can be used to insert debug code
 *  @param[in] g_p                pointer array handle
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_debug = pgp_debug
#endif

void pgp_debug(Integer g_p, Integer intsize)
{
  Integer handle;
  void *ptr;
  Integer lo[2],hi[2],ld;
  Integer i, j, idim, jdim, size;
  handle = g_p + GP_OFFSET;
  lo[0] = 1;
  lo[1] = 1;
  idim = GP[handle].dims[0];
  jdim = GP[handle].dims[1];
  hi[0] = idim;
  hi[1] = jdim;
  ld = idim;
  if (pnga_nodeid() == 0) {
    size = idim*jdim;
    if (intsize == 4) {
      ptr = (int*)malloc(size*sizeof(int));
    } else {
      ptr = (int*)malloc(size*sizeof(int64_t));
    }
    pnga_get(GP[handle].g_size_array, lo, hi, ptr, &ld);
    size = 0;
    for (i=0; i<idim; i++) {
      for (j=0; j<jdim; j++) {
        if (intsize == 4) {
          printf("  %5d",((int*)ptr)[j*idim+i]);
          size += ((int*)ptr)[j*idim+i];
        } else {
          printf("  %5d",(int)((int64_t*)ptr)[j*idim+i]);
          size += (Integer)((int64_t*)ptr)[j*idim+i];
        }
      }
      printf("\n");
    }
    printf("total size of array: %ld\n",(long)size);
    free(ptr);
  }
}

/**
 *  Return coordinates of a GP patch associated with processor proc
 *  @param[in] g_p                pointer array handle
 *  @param[in] proc               processor for which patch coordinates
 *                                are being requested
 *  @param[out] lo[ndim],hi[ndim] bounding indices of patch
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_distribution = pgp_distribution
#endif

void pgp_distribution(Integer g_p, Integer proc, Integer *lo, Integer *hi)
{
  Integer handle, ndim, i;
  handle = g_p + GP_OFFSET;
  if (pnga_nodeid() == proc) {
    ndim = pnga_ndim(GP[handle].g_ptr_array);
    for (i=0; i<ndim; i++) {
      lo[i] = GP[handle].lo[i];
      hi[i] = GP[handle].hi[i];
    }
  } else {
    pnga_distribution(GP[handle].g_ptr_array, proc, lo, hi);
  }
}

/**
 *  Assign data object to a pointer array element. Pointer array element
 *  must be on the same processor as the data object.
 *  @param[in] g_p             pointer array handle
 *  @param[in] subscript[ndim] location of element in pointer array
 *  @param[in] ptr             ptr to local data ojbect
 *  @param[in] size            size of local data ojbect
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_assign_local_element = pgp_assign_local_element
#endif

void pgp_assign_local_element(Integer g_p, Integer *subscript, void *ptr,
                              Integer size, Integer intsize)
{
  void *gp_ptr;
  Integer handle, ld[GP_MAX_DIM-1], i;
  handle = g_p + GP_OFFSET;
  /* check to make sure that element is located in local block of GP array */
  for (i=0; i<GP[handle].ndim; i++) {
    if (subscript[i]<GP[handle].lo[i] || subscript[i]>GP[handle].hi[i]) {
      /*bjp
      printf("p[%d] subscript[%d]: %d\n",pnga_nodeid(),i,subscript[i]);
      printf("p[%d] lo[%d]: %d hi[%d]: %d\n",pnga_nodeid(),i,GP[handle].lo[i],i,
             GP[handle].hi[i]);
             */
      /*
      printf("p[%d] subscript[%d]: %d lo[%d]: %d hi[%d]: %d\n",pnga_nodeid(),
        i, subscript[i], i, GP[handle].lo[i], i, GP[handle].hi[i]);
        */
      pnga_error("gp_assign_local_element: subscript out of bounds", i);
    }
  }
  pnga_access_ptr(GP[handle].g_size_array,subscript,subscript,&gp_ptr,ld);
  if (intsize == 4) {
    *((int*)gp_ptr) = (int)size;
  } else {
    *((int64_t*)gp_ptr) = (int64_t)size;
  }
  /*bjp
  printf("p[%ld] (internal) size %d at location [%ld:%ld]\n",
          (long)pnga_nodeid(), *((int*)gp_ptr),
          (long)subscript[0],(long)subscript[1]);
          */
  pnga_release_update(GP[handle].g_size_array, subscript, subscript);
  pnga_access_ptr(GP[handle].g_ptr_array,subscript,subscript,&gp_ptr,ld);
  *((armci_meminfo_t*)gp_ptr) =
    *((armci_meminfo_t*)(((char*)ptr)-sizeof(armci_meminfo_t)));
  pnga_release_update(GP[handle].g_ptr_array, subscript, subscript);
}

/**
 * Free local data element using access via the Global Pointer array.
 *  @param[in] g_p             pointer array handle
 *  @param[in] subscript[ndim] location of element in pointer array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_free_local_element = pgp_free_local_element
#endif

void* pgp_free_local_element(Integer g_p, Integer *subscript)
{
  armci_meminfo_t *gp_ptr;
  void *ptr;
  Integer handle, ld[GP_MAX_DIM-1], i;
  GP_Int buf;
  handle = g_p + GP_OFFSET;
  /* check to make sure that element is located in local block of GP array */
  for (i=0; i<GP[handle].ndim; i++) {
    if (subscript[i]<GP[handle].lo[i] || subscript[i]>GP[handle].hi[i]) {
      pnga_error("gp_free_local_element: subscript out of bounds", i);
    }
  }
  pnga_access_ptr(GP[handle].g_ptr_array,subscript,subscript,&gp_ptr,ld);
  ptr = (*gp_ptr).addr;
  memset((void*)gp_ptr,0,sizeof(armci_meminfo_t));
  pnga_release_update(GP[handle].g_ptr_array, subscript, subscript);

  /* set corresponding element of size array to zero */
  buf = 0;
  for (i=0; i<GP[handle].ndim-1; i++) {
    ld[i] = 1;
  }
  pnga_put(GP[handle].g_size_array, subscript, subscript, &buf, ld);
  return ptr;
}

/**
 * Zero all allocated bits in a Global Pointer array.
 *  @param[in] g_p             pointer array handle
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_memzero = pgp_memzero
#endif

void pgp_memzero(Integer g_p, Integer intsize)
{
  void *gp_ptr, *size_array;
  Integer handle, i, me, nelems, ndim;
  Integer lo[GP_MAX_DIM], hi[GP_MAX_DIM], ld[GP_MAX_DIM-1];
  Integer j;

  pnga_sync();
  handle = g_p + GP_OFFSET;
  me = pnga_nodeid();
  ndim = GP[handle].ndim;

  /* Determine number of elements held locally */
  pgp_distribution(g_p, me, lo, hi);
  nelems = 1;
  for (i=0; i<ndim; i++) {
    nelems *= (hi[i]-lo[i]+1);
  }

  /* Get pointers to local data elements and their sizes */
  pnga_access_ptr(GP[handle].g_ptr_array,lo,hi,&gp_ptr,ld);
  pnga_access_ptr(GP[handle].g_size_array,lo,hi,&size_array,ld);

  /* Zero bits in data elements */
  /*bjp
  printf("p[%d] nelems: %d ld[0]: %d\n",me,nelems,ld[0]);
  */
  for (i=0; i<nelems; i++) {
    /*bjp
    printf("p[%d] gp_ptr[%d].addr: %p size_array[%d]: %d\n",me,
        i,(void*)((armci_meminfo_t*)gp_ptr)[i].addr,i,
        (int)((int*)size_array)[i]);
    */
    if (intsize == 4) {
      memset((void*)((armci_meminfo_t*)gp_ptr)[i].addr, 0,
          (size_t)((int*)size_array)[i]);
#if 0
      for (j=0; j<((int*)size_array)[i]; j++) {
        if (((char*)((armci_meminfo_t*)gp_ptr)[i].addr)[j] != 0) {
          printf("p[%d] mismatch for i: %d j: %d\n",me,i,j);
        }
      }
#endif
    } else {
      memset((void*)((armci_meminfo_t*)gp_ptr)[i].addr, 0,
          (size_t)((int64_t*)size_array)[i]);
#if 0
      for (j=0; j<((int64_t*)size_array)[i]; j++) {
        if (((char*)((armci_meminfo_t*)gp_ptr)[i].addr)[j] != 0) {
          printf("p[%d] mismatch for i: %d j: %d\n",me,i,j);
        }
      }
#endif
    }
  }
  pnga_release_update(GP[handle].g_ptr_array,lo,hi);
  pnga_release_update(GP[handle].g_size_array,lo,hi);
  pnga_sync();
}

/**
 * Synchronize system and flush all outstanding communicaion.
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wgp_sync = pgp_sync
#endif

void pgp_sync()
{
  pnga_sync();
}
