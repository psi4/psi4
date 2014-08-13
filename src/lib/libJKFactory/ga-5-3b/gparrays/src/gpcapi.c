/**
 * @file gpcapi.c
 *
 * Implements the C interface.
 * These calls forward to the (possibly) weak symbols of the internal
 * implementations.
 */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "gp.h"
#include "gpbase.h"
#include "gp-papi.h"

#if ENABLE_PROFILING
#   include "gp-wapi.h"
#else
#   include "gp-wapidefs.h"
#endif

#ifdef USE_FAPI
#  define COPYC2F(carr, farr, n){\
   int i; for(i=0; i< (n); i++)(farr)[i]=(Integer)(carr)[i];}
#  define COPYF2C(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[i]=(int)(farr)[i];}
#  define COPYINDEX_F2C     COPYF2C
#else
#  define COPYC2F(carr, farr, n){\
   int i; for(i=0; i< (n); i++)(farr)[n-i-1]=(Integer)(carr)[i];}
#  define COPYINDEX_C2F(carr, farr, n){\
   int i; for(i=0; i< (n); i++)(farr)[n-i-1]=(Integer)(carr)[i]+1;}
#  define COPYINDEX_F2C(farr, carr, n){\
   int i; for(i=0; i< (n); i++)(carr)[n-i-1]=(int)(farr)[i] -1;}
#endif

static Integer* gp_copy_map(int block[], int block_ndim, int map[]);

void GP_Access_element(int g_p, int *subscript, void *ptr, int *size)
{
  Integer ag_p = (Integer)g_p;
  Integer _gp_idx[GP_MAX_DIM];
  Integer asize;
  int ndim = wgp_get_dimension(ag_p);
  COPYINDEX_C2F(subscript,_gp_idx,ndim);
  wgp_access_element(ag_p, _gp_idx, ptr, &asize);
  *size = (int)asize;
}

int GP_Allocate(int g_p)
{
  Integer ag_p;
  ag_p = (Integer)g_p;
  return (int)wgp_allocate(ag_p);
}

void GP_Assign_local_element(int g_p, int *subscript, void *ptr, int size)
{
  Integer ag_p = (Integer)g_p;
  Integer _gp_idx[GP_MAX_DIM];
  Integer asize = (Integer)size;
  int ndim = wgp_get_dimension(ag_p);
  COPYINDEX_C2F(subscript,_gp_idx,ndim);
  wgp_assign_local_element(ag_p, _gp_idx, ptr, asize, 4);
}

int GP_Create_handle()
{
  return (int)wgp_create_handle();
}

void GP_Debug(int g_p)
{
  Integer ag_p;
  ag_p = (Integer)g_p;
  wgp_debug(ag_p, 4);
}

int GP_Destroy(int g_p)
{
  Integer ag_p;
  ag_p = (Integer)g_p;
  return (int)wgp_destroy(ag_p);
}

void GP_Distribution(int g_p, int proc, int *lo, int *hi)
{
  Integer ag_p = (Integer)g_p;
  Integer p = (Integer)proc;
  int ndim = (int)wgp_get_dimension(ag_p);
  Integer _gp_lo[GP_MAX_DIM], _gp_hi[GP_MAX_DIM];
  wgp_distribution(ag_p, p, _gp_lo, _gp_hi);
  COPYINDEX_F2C(_gp_lo, lo, ndim);
  COPYINDEX_F2C(_gp_hi, hi, ndim);
}

void GP_Free(void* ptr)
{
  wgp_free(ptr);
}

void* GP_Free_local_element(int g_p, int *subscript)
{
  Integer ag_p = (Integer)g_p;
  int ndim = wgp_get_dimension(ag_p);
  Integer _gp_idx[GP_MAX_DIM];
  COPYINDEX_C2F(subscript, _gp_idx, ndim);
  return wgp_free_local_element(ag_p, _gp_idx);
}

int GP_Get_dimension(int g_p)
{
  return wgp_get_dimension(g_p);
}

void GP_Gather_size(int g_p, int nv, int *subscript, int *size)
{
  int idx, i;
  Integer ag_p = (Integer)g_p;
  Integer anv = (Integer)nv;
  Integer asize;
  Integer *asubscript;
  int ndim = wgp_get_dimension(ag_p);
  asubscript = (Integer*)malloc((int)ndim*nv*sizeof(Integer));
  if (asubscript == NULL)
    GA_Error("Memory allocation in GP_Gather_size failed",0);

  /* adjust the indices for fortran interface */
  for (idx=0; idx<nv; idx++) {
    for (i=0; i<ndim; i++) {
      asubscript[idx*ndim +(ndim-i-1)] = subscript[idx*ndim+i] + 1;
    }
  }

  wgp_gather_size(ag_p, anv, asubscript, &asize, 4);
  
  *size = (int)asize;

  free(asubscript);
}

void GP_Gather(int g_p, int nv, int *subscript, void *buf, void **buf_ptr,
               int *buf_size, int *size, int setbuf)
{
  int idx, i;
  Integer ag_p = (Integer)g_p;
  Integer anv = (Integer)nv;
  Integer aset = (Integer)setbuf;
  Integer asize;
  Integer *asubscript;
  int ndim = wgp_get_dimension(ag_p);
  asubscript = (Integer*)malloc((int)ndim*nv*sizeof(Integer));
  if (asubscript == NULL)
    GA_Error("Memory allocation in GP_Gather failed",0);

  /* adjust the indices for fortran interface */
  for (idx=0; idx<nv; idx++) 
    for (i=0; i<ndim; i++)
      asubscript[idx*ndim +(ndim-i-1)] = subscript[idx*ndim+i] + 1;

  wgp_gather(ag_p, anv, asubscript, buf, buf_ptr, buf_size, &asize, 4, aset);
  
  *size = (int)asize;

  free(asubscript);
}

void GP_Get_size(int g_p, int *lo, int *hi, int *size)
{
  Integer ag_p = (Integer)g_p;
  int ndim = wgp_get_dimension(ag_p);
  Integer asize;
  Integer _gp_lo[GP_MAX_DIM], _gp_hi[GP_MAX_DIM];
  COPYINDEX_C2F(lo, _gp_lo, ndim);
  COPYINDEX_C2F(hi, _gp_hi, ndim);
  wgp_get_size(ag_p, _gp_lo, _gp_hi, &asize, 4);
  *size = (int)asize;
}

void GP_Get(int g_p, int *lo, int *hi, void *buf, void **buf_ptr, int *ld,
            void *buf_size, int *ld_sz, int *size, int setbuf)
{
  Integer ag_p = (Integer)g_p;
  int ndim = wgp_get_dimension(ag_p);
  Integer asize;
  Integer asetbuf = (Integer)setbuf;
  Integer _gp_lo[GP_MAX_DIM], _gp_hi[GP_MAX_DIM];
  Integer _gp_ld[GP_MAX_DIM], _gp_ld_sz[GP_MAX_DIM];
  COPYINDEX_C2F(lo, _gp_lo, ndim);
  COPYINDEX_C2F(hi, _gp_hi, ndim);
  COPYC2F(ld, _gp_ld, ndim-1);
  COPYC2F(ld_sz, _gp_ld_sz, ndim-1);
  wgp_get(ag_p, _gp_lo, _gp_hi, buf, buf_ptr, _gp_ld,
          buf_size, _gp_ld_sz, &asize, 4, asetbuf);
  *size = (int)asize;
}

void GP_Initialize()
{
  wgp_initialize();
}

void* GP_Malloc(size_t size)
{
  return wgp_malloc(size);
}

void GP_Memzero(int g_p)
{
  Integer ag_p;
  ag_p = (Integer)g_p;
  wgp_memzero(ag_p, 4);
}

void GP_Put(int g_p, int *lo, int *hi, void **buf_ptr, int *ld,
            void *buf_size, int *ld_sz, int *size, int checksize)
{
  Integer ag_p = (Integer)g_p;
  int ndim = wgp_get_dimension(ag_p);
  Integer achecksize = (Integer)checksize;
  Integer asize;
  Integer _gp_lo[GP_MAX_DIM], _gp_hi[GP_MAX_DIM];
  Integer _gp_ld[GP_MAX_DIM], _gp_ld_sz[GP_MAX_DIM];
  COPYINDEX_C2F(lo, _gp_lo, ndim);
  COPYINDEX_C2F(hi, _gp_hi, ndim);
  COPYC2F(ld, _gp_ld, ndim-1);
  COPYC2F(ld_sz, _gp_ld_sz, ndim-1);
  wgp_put(ag_p, _gp_lo, _gp_hi, buf_ptr, _gp_ld,
          buf_size, _gp_ld_sz, &asize, achecksize, 4);
  *size = (int)asize;
}

void GP_Release_element(int g_p, int *subscript)
{
  Integer ag_p = (Integer)g_p;
  Integer _gp_idx[GP_MAX_DIM];
  int ndim = wgp_get_dimension(ag_p);
  COPYINDEX_C2F(subscript,_gp_idx,ndim);
  wgp_release_element(ag_p, _gp_idx);
}

void GP_Release_update_element(int g_p, int *subscript)
{
  Integer ag_p = (Integer)g_p;
  Integer _gp_idx[GP_MAX_DIM];
  int ndim = wgp_get_dimension(ag_p);
  COPYINDEX_C2F(subscript,_gp_idx,ndim);
  wgp_release_update_element(ag_p, _gp_idx);
}

void GP_Scatter(int g_p, int nv, int *subscript, void **buf_ptr,
               int *buf_size, int *size, int checksize)
{
  int idx, i;
  Integer ag_p = (Integer)g_p;
  Integer anv = (Integer)nv;
  Integer asize;
  Integer acheck = (Integer)checksize;
  Integer *asubscript;
  int ndim = wgp_get_dimension(ag_p);
  asubscript = (Integer*)malloc((int)ndim*nv*sizeof(Integer));
  if (asubscript == NULL)
    GA_Error("Memory allocation in GP_Scatter failed",0);

  /* adjust the indices for fortran interface */
  for (idx=0; idx<nv; idx++) 
    for (i=0; i<ndim; i++)
      asubscript[idx*ndim +(ndim-i-1)] = subscript[idx*ndim+i] + 1;

  wgp_scatter(ag_p, anv, asubscript, buf_ptr, buf_size, &asize, acheck, 4);
  
  *size = (int)asize;

  free(asubscript);
}

void GP_Set_chunk(int g_p, int *chunk)
{
  Integer ag_p;
  Integer _gp_chunk[GP_MAX_DIM];
  int ndim;
  ag_p = (Integer)g_p;
  ndim = (int)wgp_get_dimension(ag_p);
  COPYC2F(chunk,_gp_chunk,ndim);
  wgp_set_chunk(ag_p, _gp_chunk);
}

void GP_Set_dimensions(int g_p, int ndim, int *dims)
{
  Integer ag_p, andim;
  Integer _gp_dims[GP_MAX_DIM];
  COPYC2F(dims,_gp_dims,ndim);
  ag_p = (Integer)g_p;
  andim = (Integer)ndim;
  wgp_set_dimensions(ag_p, andim, _gp_dims, 4);
}

void GP_Set_irreg_distr(int g_p, int *mapc, int *nblock)
{
  Integer ag_p, andim;
  Integer ablock[GP_MAX_DIM];
  Integer *amapc;
  ag_p = (Integer)g_p;
  andim = (int)wgp_get_dimension(ag_p);
  COPYC2F(nblock,ablock,andim);
  amapc = gp_copy_map(nblock, (int)andim, mapc);
  wgp_set_irreg_distr(ag_p, amapc, ablock);
  free(amapc);
}

void GP_Sync()
{
  wgp_sync();
}

void GP_Terminate()
{
  wgp_terminate();
}

/* Internal utility functions */

static Integer* gp_copy_map(int block[], int block_ndim, int map[])
{
  int d;
  int i,sum=0,capi_offset=0,map_offset=0;
  Integer *_ga_map_capi;

  for (d=0; d<block_ndim; d++) {
    sum += block[d];
  }

  _ga_map_capi = (Integer*)malloc(sum * sizeof(Integer));

  capi_offset = sum;
  for (d=0; d<block_ndim; d++) {
    capi_offset -= block[d];
    for (i=0; i<block[d]; i++) {
      _ga_map_capi[capi_offset+i] = map[map_offset+i] + 1;
    }
    map_offset += block[d];
  }

  return _ga_map_capi;
}
