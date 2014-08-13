/**
 * @file fapi.c
 *
 * Implements the Fortran interface.
 * These calls forward to the (possibly) weak symbols of the internal
 * implementations.
 */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#include "cnames.h"
#include "farg.h"
#include "globalp.h"
#include "macommon.h"
#include "matmul.h"

int _ga_initialize_f=0;

#include "ga-papi.h"
#if ENABLE_PROFILING
#   include "ga-wapi.h"
#else
#   include "ga-wapidefs.h"
#endif

#define FNAM 31
#define FMSG 256

/* Routines from base.c */

logical FATR ga_allocate_(Integer *g_a)
{
  return wnga_allocate(*g_a);
}

logical FATR nga_allocate_(Integer *g_a)
{
  return wnga_allocate(*g_a);
}

logical FATR ga_compare_distr_(Integer *g_a, Integer *g_b)
{
  return wnga_compare_distr(*g_a, *g_b);
}

logical FATR nga_compare_distr_(Integer *g_a, Integer *g_b)
{
  return wnga_compare_distr(*g_a, *g_b);
}

logical FATR ga_create_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        type, dim1, dim2, array_name, chunk1, chunk2, g_a, slen
#else
        type, dim1, dim2, array_name, slen, chunk1, chunk2, g_a
#endif
        )
Integer *type, *dim1, *dim2, *chunk1, *chunk2, *g_a;
int slen;
char* array_name;
{
  char buf[FNAM];
  Integer ndim, dims[2], chunk[2];
  ga_f2cstring(array_name ,slen, buf, FNAM);
  dims[0] = *dim1;
  dims[1] = *dim2;
  ndim = 2;
  chunk[0] = (*chunk1==0)? -1 : *chunk1;
  chunk[1] = (*chunk2==0)? -1 : *chunk2;

  return(wnga_create(*type, ndim, dims, buf, chunk, g_a));
}

logical FATR nga_create_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
    Integer *type, Integer *ndim,
    Integer *dims, char* array_name, Integer *chunk,
    Integer *g_a, int slen
#else
    Integer *type, Integer *ndim,
    Integer *dims, char* array_name, int slen,
    Integer *p_handle, Integer *g_a
#endif
    )
{
  char buf[FNAM];
  ga_f2cstring(array_name ,slen, buf, FNAM);

  return (wnga_create(*type, *ndim,  dims, buf, chunk, g_a));
}

logical FATR nga_create_config_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
    Integer *type, Integer *ndim,
    Integer *dims, char* array_name, Integer *chunk,
    Integer *p_handle, Integer *g_a, int
    slen
#else
    Integer *type, Integer *ndim,
    Integer *dims, char* array_name, int slen,
    Integer *chunk, Integer *p_handle,
    Integer *g_a
#endif
    )
{
  char buf[FNAM];
  ga_f2cstring(array_name ,slen, buf, FNAM);

  return (wnga_create_config(*type, *ndim,  dims, buf, chunk, *p_handle,
                             g_a));
}

logical FATR nga_create_ghosts_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
    Integer *type, Integer *ndim, Integer *dims,
    Integer *width, char* array_name, Integer *chunk, Integer *g_a,
    int slen
#else
    Integer *type, Integer *ndim, Integer *dims,
    Integer *width, char* array_name, int slen,
    Integer *chunk, Integer *g_a
#endif
    )
{
  char buf[FNAM];
  ga_f2cstring(array_name ,slen, buf, FNAM);

  return (wnga_create_ghosts(*type, *ndim,  dims, width, buf, chunk, g_a));
}

logical FATR nga_create_ghosts_config_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
    Integer *type, Integer *ndim,
    Integer *dims, Integer *width, char* array_name,
    Integer *chunk, Integer *p_handle,
    Integer *g_a,
    int slen
#else
    Integer *type, Integer *ndim,
    Integer *dims, Integer *width, char* array_name,
    int slen,
    Integer *chunk,
    Integer *p_handle,
    Integer *g_a
#endif
    )
{
  char buf[FNAM];
  ga_f2cstring(array_name ,slen, buf, FNAM);

  return (wnga_create_ghosts_config(*type, *ndim,  dims, width, buf, chunk,
                                    *p_handle, g_a));
}

logical FATR nga_create_ghosts_irreg_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
    Integer *type, Integer *ndim,
    Integer *dims, Integer width[], char* array_name, Integer map[],
    Integer block[], Integer *g_a, int slen
#else
    Integer *type, Integer *ndim,
    Integer *dims, Integer width[], char* array_name, int slen,
    Integer map[], Integer block[], Integer *g_a
#endif
    )
{
  char buf[FNAM];
  Integer st;
  ga_f2cstring(array_name ,slen, buf, FNAM);

  st = wnga_create_ghosts_irreg(*type, *ndim,  dims, width,
                                buf, map, block, g_a);
  return st;
}

logical FATR nga_create_ghosts_irreg_config_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
    Integer *type,
    Integer *ndim, Integer *dims, Integer width[], char* array_name,
    Integer map[], Integer block[], Integer *p_handle, Integer *g_a,
    int slen
#else
    Integer *type,
    Integer *ndim, Integer *dims, Integer width[], char* array_name,
    int slen, Integer map[], Integer block[],
    Integer *p_handle, Integer *g_a
#endif
    )
{
  char buf[FNAM];
  Integer st;
  ga_f2cstring(array_name ,slen, buf, FNAM);

  st = wnga_create_ghosts_irreg_config(*type, *ndim,  dims,
      width, buf, map, block, *p_handle, g_a);
  return st;
}

logical FATR ga_create_irreg_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
    Integer *type, Integer *dim1, Integer *dim2, char *array_name,
    Integer *map1, Integer *nblock1, Integer *map2, Integer *nblock2,
    Integer *g_a, int slen
#else
    Integer *type, Integer *dim1, Integer *dim2, char *array_name, int
    slen, Integer *map1, Integer *nblock1, Integer *map2, Integer
    *nblock2, Integer *g_a
#endif
    )
{
  char buf[FNAM];
  Integer i, ndim, dims[2], block[2], *map;
  Integer status;
  ga_f2cstring(array_name ,slen, buf, FNAM);
  dims[0] = *dim1;
  dims[1] = *dim2;
  block[0] = *nblock1;
  block[1] = *nblock2;
  ndim = 2;
  map = (Integer*)malloc((wnga_nnodes()+1)*sizeof(Integer));
  for(i=0; i<*nblock1; i++) map[i] = map1[i];
  for(i=0; i<*nblock2; i++) map[i+ *nblock1] = map2[i];
  status = wnga_create_irreg(*type, ndim, dims, buf, map, block, g_a);
  free(map);
  return status;
}

logical FATR nga_create_irreg_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
    Integer *type, Integer *ndim, Integer *dims,
    char* array_name, Integer map[], Integer block[],
    Integer *g_a, int slen
#else
    Integer *type, Integer *ndim, Integer *dims,
    char* array_name, int slen,
    Integer map[], Integer block[], Integer *g_a
#endif
    )
{
  char buf[FNAM];
  Integer st;
  ga_f2cstring(array_name ,slen, buf, FNAM);

  st = wnga_create_irreg(*type, *ndim,  dims, buf, map, block, g_a);
  return st;
}

logical FATR nga_create_irreg_config_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
    Integer *type, Integer *ndim,
    Integer *dims, char* array_name, Integer map[],
    Integer block[], Integer *p_handle, Integer *g_a,
    int slen
#else
    Integer *type, Integer *ndim,
    Integer *dims, char* array_name, int slen, Integer map[],
    Integer block[], Integer *p_handle, Integer *g_a
#endif
    )
{
  char buf[FNAM];
  Integer st;
  ga_f2cstring(array_name ,slen, buf, FNAM);

  st = wnga_create_irreg_config(*type, *ndim,  dims, buf, map,
      block, *p_handle, g_a);
  return st;
}

Integer FATR ga_create_handle_()
{
  return wnga_create_handle();
}

Integer FATR nga_create_handle_()
{
  return wnga_create_handle();
}

logical FATR ga_create_mutexes_(Integer *num)
{
  return wnga_create_mutexes(*num);
}

logical FATR nga_create_mutexes_(Integer *num)
{
  return wnga_create_mutexes(*num);
}

logical FATR ga_destroy_(Integer *g_a)
{
  return wnga_destroy(*g_a);
}

logical FATR nga_destroy_(Integer *g_a)
{
  return wnga_destroy(*g_a);
}

logical FATR ga_destroy_mutexes_()
{
  return wnga_destroy_mutexes();
}

logical FATR nga_destroy_mutexes_()
{
  return wnga_destroy_mutexes();
}

void FATR ga_distribution_(Integer *g_a, Integer *proc, Integer *ilo,
                           Integer *ihi, Integer *jlo, Integer *jhi)
{
  Integer lo[2], hi[2];
  wnga_distribution(*g_a, *proc, lo, hi);
  *ilo = lo[0];
  *jlo = lo[1];
  *ihi = hi[0];
  *jhi = hi[1];
}

void FATR nga_distribution_(Integer *g_a, Integer *proc, Integer *lo, Integer *hi)
{
  wnga_distribution(*g_a, *proc, lo, hi);
}

logical FATR ga_duplicate_( Integer *g_a, Integer *g_b, char *array_name, int slen)
{
  char buf[FNAM];

  ga_f2cstring(array_name ,slen, buf, FNAM);
  return(wnga_duplicate(*g_a, g_b, buf));
}

logical FATR nga_duplicate_( Integer *g_a, Integer *g_b, char *array_name, int slen)
{
  char buf[FNAM];

  ga_f2cstring(array_name ,slen, buf, FNAM);
  return(wnga_duplicate(*g_a, g_b, buf));
}

void FATR ga_fill_(Integer *g_a, void* val)
{
  wnga_fill(*g_a, val);
}

void FATR nga_fill_(Integer *g_a, void* val)
{
  wnga_fill(*g_a, val);
}

void FATR ga_get_block_info_(Integer *g_a, Integer *num_blocks,
                             Integer *block_dims)
{
  wnga_get_block_info(*g_a, num_blocks, block_dims);
}

void FATR nga_get_block_info_(Integer *g_a, Integer *num_blocks,
                             Integer *block_dims)
{
  wnga_get_block_info(*g_a, num_blocks, block_dims);
}

logical FATR ga_get_debug_()
{
  return wnga_get_debug();
}

logical FATR nga_get_debug_()
{
  return wnga_get_debug();
}

Integer FATR ga_get_dimension_(Integer *g_a)
{
  return wnga_get_dimension(*g_a);
}

Integer FATR nga_get_dimension_(Integer *g_a)
{
  return wnga_get_dimension(*g_a);
}

Integer FATR ga_get_pgroup_(Integer *g_a)
{
  return wnga_get_pgroup(*g_a);
}

Integer FATR nga_get_pgroup_(Integer *g_a)
{
  return wnga_get_pgroup(*g_a);
}

Integer FATR ga_get_pgroup_size_(Integer *grp_id)
{
  return wnga_get_pgroup_size(*grp_id);
}

Integer FATR nga_get_pgroup_size_(Integer *grp_id)
{
  return wnga_get_pgroup_size(*grp_id);
}

void FATR ga_get_proc_grid_(Integer *g_a, Integer *dims)
{
  wnga_get_proc_grid(*g_a, dims);
}

void FATR ga_get_proc_index_(Integer *g_a, Integer *iproc, Integer *index)
{
  wnga_get_proc_index(*g_a, *iproc, index);
}

void FATR nga_get_proc_index_(Integer *g_a, Integer *iproc, Integer *index)
{
  wnga_get_proc_index(*g_a, *iproc, index);
}

void FATR nga_get_proc_grid_(Integer *g_a, Integer *dims)
{
  wnga_get_proc_grid(*g_a, dims);
}

logical FATR ga_has_ghosts_(Integer *g_a)
{
  return wnga_has_ghosts(*g_a);
}

logical FATR nga_has_ghosts_(Integer *g_a)
{
  return wnga_has_ghosts(*g_a);
}

void FATR ga_initialize_()
{
  _ga_initialize_f=1;
  wnga_initialize();
}

void FATR nga_initialize_()
{
  _ga_initialize_f=1;
  wnga_initialize();
}

void FATR ga_initialize_ltd_(Integer *limit)
{
  _ga_initialize_f=1;
  wnga_initialize_ltd(*limit);
}

void FATR nga_initialize_ltd_(Integer *limit)
{
  _ga_initialize_f=1;
  wnga_initialize_ltd(*limit);
}

logical FATR ga_initialized_()
{
  return wnga_initialized();
}

logical FATR nga_initialized_()
{
  return wnga_initialized();
}

void FATR ga_inquire_(Integer *g_a, Integer *type, Integer *dim1, Integer *dim2)
{
  Integer dims[2], ndim;
  wnga_inquire(*g_a, type, &ndim, dims);
  if (ndim != 2) wnga_error("Wrong array dimension in ga_inquire",ndim);
  *type = pnga_type_c2f(*type);
  *dim1 = dims[0];
  *dim2 = dims[1];
}

void FATR nga_inquire_(Integer *g_a, Integer *type, Integer *ndim, Integer *dims)
{
  wnga_inquire(*g_a, type, ndim, dims);
  *type = pnga_type_c2f(*type);
}

Integer FATR ga_inquire_memory_()
{
  return wnga_inquire_memory();
}

Integer FATR nga_inquire_memory_()
{
  return wnga_inquire_memory();
}

void FATR ga_inquire_name_(Integer *g_a, char *array_name, int len)
{
  char *c_name;
  wnga_inquire_name(*g_a, &c_name);
  ga_c2fstring(c_name, array_name, len);
}

void FATR nga_inquire_name_(Integer *g_a, char *array_name, int len)
{
  char *c_name;
  wnga_inquire_name(*g_a, &c_name);
  ga_c2fstring(c_name, array_name, len);
}

logical FATR ga_is_mirrored_(Integer *g_a)
{
  return wnga_is_mirrored(*g_a);
}

logical FATR nga_is_mirrored_(Integer *g_a)
{
  return wnga_is_mirrored(*g_a);
}

void FATR ga_list_nodeid_(Integer *list, Integer *nprocs)
{
  wnga_list_nodeid(list, *nprocs);
}

void FATR nga_list_nodeid_(Integer *list, Integer *nprocs)
{
  wnga_list_nodeid(list, *nprocs);
}

Integer FATR nga_locate_num_blocks_(Integer *g_a, Integer *lo, Integer *hi)
{
  return wnga_locate_num_blocks(*g_a,lo,hi);
}

logical FATR ga_locate_nnodes_( Integer *g_a,
                                Integer *ilo,
                                Integer *ihi,
                                Integer *jlo,
                                Integer *jhi,
                                Integer *np)
{
  Integer lo[2], hi[2];
  lo[0] = *ilo;
  lo[1] = *jlo;
  hi[0] = *ihi;
  hi[1] = *jhi;
  return wnga_locate_nnodes(*g_a, lo, hi, np);
}

logical FATR nga_locate_nnodes_(Integer *g_a, Integer *lo, Integer *hi, Integer *np)
{
    return wnga_locate_nnodes(*g_a, lo, hi, np);
}

logical FATR ga_locate_region_( Integer *g_a,
                                Integer *ilo,
                                Integer *ihi,
                                Integer *jlo,
                                Integer *jhi,
                                Integer map[][5],
                                Integer *np)
{
  logical status = FALSE;
  Integer lo[2], hi[2], p;
  Integer *mapl, *proclist;
  proclist = (Integer*)malloc(wnga_nnodes()*sizeof(Integer));
  mapl = (Integer*)malloc(5*wnga_nnodes()*sizeof(Integer));
  lo[0] = *ilo;
  lo[1] = *jlo;
  hi[0] = *ihi;
  hi[1] = *jhi;
  if (wnga_locate_num_blocks(*g_a, lo, hi) == -1) {
    status = wnga_locate_region(*g_a, lo, hi, mapl, proclist, np);
    /* need to swap elements (ilo,jlo,ihi,jhi) -> (ilo,ihi,jlo,jhi) */
    for(p = 0; p< *np; p++){
      map[p][0] = mapl[4*p];
      map[p][1] = mapl[4*p + 2];
      map[p][2] = mapl[4*p + 1];
      map[p][3] = mapl[4*p + 3];
      map[p][4] = proclist[p];
    }
  } else {
    wnga_error("Must call nga_locate_region on block-cyclic data distribution",0);
  }
  free(proclist);
  free(mapl);
  return status;
}

logical FATR nga_locate_region_( Integer *g_a,
                                 Integer *lo,
                                 Integer *hi,
                                 Integer *map,
                                 Integer *proclist,
                                 Integer *np)
{
  return wnga_locate_region(*g_a, lo, hi, map, proclist, np);
}

void FATR ga_lock_(Integer *mutex)
{
  wnga_lock(*mutex);
}

void FATR nga_lock_(Integer *mutex)
{
  wnga_lock(*mutex);
}

logical FATR ga_locate_(Integer *g_a, Integer *i, Integer *j, Integer *owner)
{
  Integer subscript[2];
  subscript[0] = *i;
  subscript[1] = *j;
  return wnga_locate(*g_a, subscript, owner);
}

logical FATR nga_locate_(Integer *g_a, Integer *subscript, Integer *owner)
{
  return wnga_locate(*g_a, subscript, owner);
}

void FATR ga_mask_sync_(Integer *begin, Integer *end)
{
  wnga_mask_sync(*begin, *end);
}

void FATR nga_mask_sync_(Integer *begin, Integer *end)
{
  wnga_mask_sync(*begin, *end);
}

Integer FATR ga_memory_avail_()
{
  return wnga_memory_avail();
}

Integer FATR nga_memory_avail_()
{
  return wnga_memory_avail();
}

logical FATR ga_memory_limited_()
{
  return wnga_memory_limited();
}

logical FATR nga_memory_limited_()
{
  return wnga_memory_limited();
}

void FATR nga_merge_distr_patch_(Integer *g_a, Integer *alo, Integer *ahi,
                                 Integer *g_b, Integer *blo, Integer *bhi)
{
  wnga_merge_distr_patch(*g_a, alo, ahi, *g_b, blo, bhi);
}

void FATR ga_merge_mirrored_(Integer *g_a)
{
  wnga_merge_mirrored(*g_a);
}

void FATR nga_merge_mirrored_(Integer *g_a)
{
  wnga_merge_mirrored(*g_a);
}

void FATR ga_nblock_(Integer *g_a, Integer *nblock)
{
  wnga_nblock(*g_a, nblock);
}

void FATR nga_nblock_(Integer *g_a, Integer *nblock)
{
  wnga_nblock(*g_a, nblock);
}

Integer FATR ga_ndim_(Integer *g_a)
{
  return wnga_ndim(*g_a);
}

Integer FATR nga_ndim_(Integer *g_a)
{
  return wnga_ndim(*g_a);
}

Integer FATR ga_nnodes_()
{
  return wnga_nnodes();
}

Integer FATR nga_nnodes_()
{
  return wnga_nnodes();
}

Integer FATR ga_nodeid_()
{
  return wnga_nodeid();
}

Integer FATR nga_nodeid_()
{
  return wnga_nodeid();
}

Integer FATR ga_pgroup_absolute_id_(Integer *grp, Integer *pid)
{
  return wnga_pgroup_absolute_id(*grp, *pid);
}

Integer FATR nga_pgroup_absolute_id_(Integer *grp, Integer *pid)
{
  return wnga_pgroup_absolute_id(*grp, *pid);
}

Integer FATR ga_pgroup_create_(Integer *list, Integer *count)
{
  return wnga_pgroup_create(list, *count);
}

Integer FATR nga_pgroup_create_(Integer *list, Integer *count)
{
  return wnga_pgroup_create(list, *count);
}

logical FATR ga_pgroup_destroy_(Integer *grp)
{
  return wnga_pgroup_destroy(*grp);
}

logical FATR nga_pgroup_destroy_(Integer *grp)
{
  return wnga_pgroup_destroy(*grp);
}

Integer FATR ga_pgroup_get_default_()
{
  return wnga_pgroup_get_default();
}

Integer FATR nga_pgroup_get_default_()
{
  return wnga_pgroup_get_default();
}

Integer FATR ga_pgroup_get_mirror_()
{
  return wnga_pgroup_get_mirror();
}

Integer FATR nga_pgroup_get_mirror_()
{
  return wnga_pgroup_get_mirror();
}

Integer FATR ga_pgroup_get_world_()
{
  return wnga_pgroup_get_world();
}

Integer FATR nga_pgroup_get_world_()
{
  return wnga_pgroup_get_world();
}

void FATR ga_pgroup_set_default_(Integer *grp)
{
  wnga_pgroup_set_default(*grp);
}

void FATR nga_pgroup_set_default_(Integer *grp)
{
  wnga_pgroup_set_default(*grp);
}

Integer FATR ga_pgroup_split_(Integer *grp, Integer *grp_num)
{
  return wnga_pgroup_split(*grp, *grp_num);
}

Integer FATR nga_pgroup_split_(Integer *grp, Integer *grp_num)
{
  return wnga_pgroup_split(*grp, *grp_num);
}

Integer FATR ga_pgroup_split_irreg_(Integer *grp, Integer *mycolor)
{
  return wnga_pgroup_split_irreg(*grp, *mycolor);
}

Integer FATR nga_pgroup_split_irreg_(Integer *grp, Integer *mycolor)
{
  return wnga_pgroup_split_irreg(*grp, *mycolor);
}

Integer FATR ga_pgroup_nnodes_(Integer *grp)
{
  return wnga_pgroup_nnodes(*grp);
}

Integer FATR nga_pgroup_nnodes_(Integer *grp)
{
  return wnga_pgroup_nnodes(*grp);
}

Integer FATR ga_pgroup_nodeid_(Integer *grp)
{
  return wnga_pgroup_nodeid(*grp);
}

Integer FATR nga_pgroup_nodeid_(Integer *grp)
{
  return wnga_pgroup_nodeid(*grp);
}

void FATR ga_proc_topology_(Integer* g_a, Integer* proc, Integer* pr,
                            Integer *pc)
{
  Integer subscript[2];
  wnga_proc_topology(*g_a, *proc, subscript);
  *pr = subscript[0];
  *pc = subscript[1];
}

void FATR nga_proc_topology_(Integer* g_a, Integer* proc, Integer* subscript)
{
  wnga_proc_topology(*g_a, *proc, subscript);
}

void FATR ga_randomize_(Integer *g_a, void* val)
{
  wnga_randomize(*g_a, val);
}

void FATR nga_randomize_(Integer *g_a, void* val)
{
  wnga_randomize(*g_a, val);
}

void FATR ga_set_array_name_(Integer *g_a, char *array_name, int slen)
{
  char buf[FNAM];
  ga_f2cstring(array_name ,slen, buf, FNAM);
  wnga_set_array_name(*g_a, buf);
}

void FATR nga_set_array_name_(Integer *g_a, char *array_name, int slen)
{
  char buf[FNAM];
  ga_f2cstring(array_name ,slen, buf, FNAM);
  wnga_set_array_name(*g_a, buf);
}

void FATR ga_set_block_cyclic_(Integer *g_a, Integer *dims)
{
  wnga_set_block_cyclic(*g_a, dims);
}

void FATR nga_set_block_cyclic_(Integer *g_a, Integer *dims)
{
  wnga_set_block_cyclic(*g_a, dims);
}

void FATR ga_set_block_cyclic_proc_grid_(Integer *g_a, Integer *dims,
                                         Integer *proc_grid)
{
  wnga_set_block_cyclic_proc_grid(*g_a, dims, proc_grid);
}

void FATR nga_set_block_cyclic_proc_grid_(Integer *g_a, Integer *dims,
                                          Integer *proc_grid)
{
  wnga_set_block_cyclic_proc_grid(*g_a, dims, proc_grid);
}

void FATR ga_set_chunk_(Integer *g_a, Integer *chunk)
{
  wnga_set_chunk(*g_a, chunk);
}

void FATR nga_set_chunk_(Integer *g_a, Integer *chunk)
{
  wnga_set_chunk(*g_a, chunk);
}

void FATR ga_set_data_(Integer *g_a, Integer *ndim, Integer *dims,
                       Integer *type)
{
  wnga_set_data(*g_a, *ndim, dims, *type);
}

void FATR nga_set_data_(Integer *g_a, Integer *ndim, Integer *dims,
                        Integer *type)
{
  wnga_set_data(*g_a, *ndim, dims, *type);
}

void FATR ga_set_debug_(logical *flag)
{
  wnga_set_debug(*flag);
}

void FATR nga_set_debug_(logical *flag)
{
  wnga_set_debug(*flag);
}

void FATR ga_set_ghosts_(Integer *g_a, Integer *width)
{
  wnga_set_ghosts(*g_a,width);
}

void FATR nga_set_ghosts_(Integer *g_a, Integer *width)
{
  wnga_set_ghosts(*g_a,width);
}

void FATR ga_set_irreg_distr_(Integer *g_a, Integer *mapc, Integer *nblock)
{
  wnga_set_irreg_distr(*g_a, mapc, nblock);
}

void FATR nga_set_irreg_distr_(Integer *g_a, Integer *mapc, Integer *nblock)
{
  wnga_set_irreg_distr(*g_a, mapc, nblock);
}

void FATR ga_set_irreg_flag_(Integer *g_a, logical *flag)
{
  wnga_set_irreg_flag(*g_a, *flag);
}

void FATR nga_set_irreg_flag_(Integer *g_a, logical *flag)
{
  wnga_set_irreg_flag(*g_a, *flag);
}

void FATR ga_set_memory_limit_(Integer *mem_limit)
{
  wnga_set_memory_limit(*mem_limit);
}

void FATR nga_set_memory_limit_(Integer *mem_limit)
{
  wnga_set_memory_limit(*mem_limit);
}

void FATR ga_set_pgroup_(Integer *g_a, Integer *p_handle)
{
  wnga_set_pgroup(*g_a, *p_handle);
}

void FATR nga_set_pgroup_(Integer *g_a, Integer *p_handle)
{
  wnga_set_pgroup(*g_a, *p_handle);
}

void FATR ga_set_restricted_(Integer *g_a, Integer *list, Integer *size)
{
  wnga_set_restricted(*g_a, list, *size);
}

void FATR nga_set_restricted_(Integer *g_a, Integer *list, Integer *size)
{
  wnga_set_restricted(*g_a, list, *size);
}

void FATR ga_set_restricted_range_(Integer *g_a, Integer *lo_proc, Integer *hi_proc)
{
  wnga_set_restricted_range(*g_a, *lo_proc, *hi_proc);
}

void FATR nga_set_restricted_range_(Integer *g_a, Integer *lo_proc, Integer *hi_proc)
{
  wnga_set_restricted_range(*g_a, *lo_proc, *hi_proc);
}

void FATR  ga_terminate_()
{
  wnga_terminate();
  _ga_initialize_f=0;
}

void FATR  nga_terminate_()
{
  wnga_terminate();
  _ga_initialize_f=0;
}

Integer FATR ga_total_blocks_(Integer *g_a)
{
  return wnga_total_blocks(*g_a);
}

Integer FATR nga_total_blocks_(Integer *g_a)
{
  return wnga_total_blocks(*g_a);
}

void FATR ga_unlock_(Integer *mutex)
{
  wnga_unlock(*mutex);
}

void FATR nga_unlock_(Integer *mutex)
{
  wnga_unlock(*mutex);
}

logical FATR ga_uses_ma_()
{
  return wnga_uses_ma();
}

logical FATR nga_uses_ma_()
{
  return wnga_uses_ma();
}

logical FATR ga_uses_proc_grid_(Integer *g_a)
{
  return wnga_uses_proc_grid(*g_a);
}

logical FATR nga_uses_proc_grid_(Integer *g_a)
{
  return wnga_uses_proc_grid(*g_a);
}

logical FATR ga_valid_handle_(Integer *g_a)
{
  return wnga_valid_handle(*g_a);
}

logical FATR nga_valid_handle_(Integer *g_a)
{
  return wnga_valid_handle(*g_a);
}

Integer FATR ga_verify_handle_(Integer *g_a)
{
  return wnga_verify_handle(*g_a);
}

Integer FATR nga_verify_handle_(Integer *g_a)
{
  return wnga_verify_handle(*g_a);
}

void FATR ga_check_handle_(Integer *g_a, char *fstring, int slen)
{
    char buf[FMSG];

    ga_f2cstring(fstring , slen, buf, FMSG);
    wnga_check_handle(*g_a, buf);
}

void FATR nga_check_handle_(Integer *g_a, char *fstring, int slen)
{
    char buf[FMSG];

    ga_f2cstring(fstring , slen, buf, FMSG);
    wnga_check_handle(*g_a, buf);
}

/* Routines from onesided.c */

void FATR ga_acc_(Integer *g_a, Integer *ilo, Integer *ihi,
                  Integer *jlo, Integer *jhi, void *buf, Integer *ld,
                  void *alpha)
{
    Integer lo[2], hi[2];
    lo[0]=*ilo;
    lo[1]=*jlo;
    hi[0]=*ihi;
    hi[1]=*jhi;
    wnga_acc(*g_a, lo, hi, buf, ld, alpha);
}

void FATR nga_acc_(Integer *g_a, Integer *lo, Integer *hi,
                   void *buf, Integer *ld, void *alpha)
{
    wnga_acc(*g_a, lo, hi, buf, ld, alpha);
}

void FATR ga_access_(Integer *g_a, Integer *ilo, Integer *ihi,
                     Integer *jlo, Integer *jhi, AccessIndex* index,
                     Integer *ld)
{
    Integer lo[2], hi[2];
    lo[0]=*ilo;
    lo[1]=*jlo;
    hi[0]=*ihi;
    hi[1]=*jhi;
    wnga_access_idx(*g_a, lo, hi, index, ld);
}

void FATR nga_access_(Integer* g_a, Integer lo[], Integer hi[],
                      AccessIndex* index, Integer ld[])
{
    wnga_access_idx(*g_a, lo, hi, index, ld);
}

void FATR nga_access_block_(Integer* g_a, Integer* idx,
                            AccessIndex* index, Integer *ld)
{
  wnga_access_block_idx(*g_a, *idx, index, ld);
}

void FATR nga_access_block_grid_(Integer* g_a, Integer* subscript,
                                 AccessIndex *index, Integer *ld)
{
  wnga_access_block_grid_idx(*g_a, subscript, index, ld);
}

void FATR nga_access_block_segment_(Integer* g_a, Integer *proc,
                                    AccessIndex* index, Integer *len)
{
  wnga_access_block_segment_idx(*g_a, *proc, index, len);
}

void FATR nga_alloc_gatscat_buf_(Integer *nelems)
{
  wnga_alloc_gatscat_buf(*nelems);
}

void FATR ga_fence_()
{
  wnga_fence();
}

void FATR nga_fence_()
{
  wnga_fence();
}

void FATR nga_free_gatscat_buf_()
{
  wnga_free_gatscat_buf();
}

void FATR  ga_gather_(Integer *g_a, void *v, Integer *i, Integer *j,
                      Integer *nv)
{
  wnga_gather2d(*g_a, v, i, j, *nv);
}

void FATR nga_gather_(Integer *g_a, void* v, Integer subscript[], Integer *nv)
{
  wnga_gather(*g_a, v, subscript, 0, *nv);
}

void FATR ga_get_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo,
                  Integer *jhi, void *buf, Integer *ld)
{
  Integer lo[2], hi[2];
  lo[0] = *ilo;
  lo[1] = *jlo;
  hi[0] = *ihi;
  hi[1] = *jhi;
  wnga_get(*g_a, lo, hi, buf, ld);
}


void FATR nga_get_(Integer *g_a, Integer *lo, Integer *hi,
                   void *buf, Integer *ld)
{
  wnga_get(*g_a, lo, hi, buf, ld);
}

void FATR ga_init_fence_()
{
  wnga_init_fence();
}

void FATR nga_init_fence_()
{
  wnga_init_fence();
}

void FATR ga_nbacc_(Integer *g_a, Integer *ilo, Integer *ihi,
                    Integer *jlo, Integer *jhi, void *buf, Integer *ld,
                    void *alpha, Integer *nbhandle)
{
  Integer lo[2], hi[2];
  lo[0]=*ilo;
  lo[1]=*jlo;
  hi[0]=*ihi;
  hi[1]=*jhi;
  wnga_nbacc(*g_a, lo, hi, buf, ld, alpha, nbhandle);
}

void FATR nga_nbacc_(Integer *g_a, Integer *lo, Integer *hi,
                     void *buf, Integer *ld, void *alpha, Integer *nbhandle)
{
  wnga_nbacc(*g_a, lo, hi, buf, ld, alpha, nbhandle);
}

void FATR ga_nbget_(Integer *g_a, Integer *ilo, Integer *ihi,
                    Integer *jlo, Integer *jhi, void *buf,
                    Integer *ld, Integer *nbhandle)
{
  Integer lo[2], hi[2];
  lo[0]=*ilo;
  lo[1]=*jlo;
  hi[0]=*ihi;
  hi[1]=*jhi;
  wnga_nbget(*g_a, lo, hi, buf, ld, nbhandle);
}

void FATR nga_nbget_(Integer *g_a, Integer *lo,
                     Integer *hi, void *buf, Integer *ld,
                     Integer *nbhandle)
{
  wnga_nbget(*g_a, lo, hi, buf, ld, nbhandle);
}

void FATR ga_nbput_(Integer *g_a, Integer *ilo, Integer *ihi,
                    Integer *jlo, Integer *jhi, void *buf,
                    Integer *ld, Integer *nbhandle)
{
  Integer lo[2], hi[2];
  lo[0]=*ilo;
  lo[1]=*jlo;
  hi[0]=*ihi;
  hi[1]=*jhi;
  wnga_nbput(*g_a, lo, hi, buf, ld, nbhandle);
}

void FATR nga_nbput_(Integer *g_a, Integer *lo,
                     Integer *hi, void *buf, Integer *ld,
                     Integer *nbhandle)
{
  wnga_nbput(*g_a, lo, hi, buf, ld, nbhandle);
}

void FATR nga_nbput_notify_(Integer *g_a, Integer *lo, Integer *hi, void *buf, Integer *ld,
			    Integer *g_b, Integer *ecoords, void *bufn, Integer *nbhandle)
{
  wnga_nbput_notify(*g_a, lo, hi, buf, ld, *g_b, ecoords, bufn, nbhandle);
} /* nga_nbput_notify_ */

void FATR nga_nbwait_notify_(Integer *nbhandle)
{
  wnga_nbwait_notify(nbhandle);
} /* nga_nbwait_notify_ */

Integer FATR ga_nbtest_(Integer *nbhandle)
{
  return wnga_nbtest(nbhandle);
}

Integer FATR nga_nbtest_(Integer *nbhandle)
{
  return wnga_nbtest(nbhandle);
}

void FATR ga_nbwait_(Integer *nbhandle)
{
  wnga_nbwait(nbhandle);
}

void FATR nga_nbwait_(Integer *nbhandle)
{
  wnga_nbwait(nbhandle);
}

void FATR ga_pgroup_sync_(Integer *grp_id)
{
  wnga_pgroup_sync(*grp_id);
}

void FATR nga_pgroup_sync_(Integer *grp_id)
{
  wnga_pgroup_sync(*grp_id);
}

void FATR ga_put_(Integer *g_a, Integer *ilo, Integer *ihi,
                  Integer *jlo, Integer *jhi, void *buf, Integer *ld)
{
  Integer lo[2], hi[2];
  lo[0]=*ilo;
  lo[1]=*jlo;
  hi[0]=*ihi;
  hi[1]=*jhi;
  wnga_put(*g_a, lo, hi, buf, ld);
}

void FATR nga_put_(Integer *g_a, Integer *lo,
                   Integer *hi, void *buf, Integer *ld)
{
  wnga_put(*g_a, lo, hi, buf, ld);
}

Integer FATR ga_read_inc_(Integer *g_a, Integer *i, Integer *j,
                          Integer *inc)
{
  Integer subscript[2];
  subscript[0] = *i;
  subscript[1] = *j;
  return wnga_read_inc(*g_a, subscript, *inc);
}

Integer FATR nga_read_inc_(Integer *g_a, Integer *subscript,
                           Integer *inc)
{
  return wnga_read_inc(*g_a, subscript, *inc);
}

void FATR ga_release_(Integer *g_a, Integer *ilo, Integer *ihi,
                      Integer *jlo, Integer *jhi)
{
  Integer lo[2], hi[2];
  lo[0]=*ilo;
  lo[1]=*jlo;
  hi[0]=*ihi;
  hi[1]=*jhi;
  wnga_release(*g_a, lo, hi);
}

void FATR nga_release_(Integer *g_a, Integer *lo, Integer *hi)
{
  wnga_release(*g_a, lo, hi);
}

void FATR nga_release_block_(Integer *g_a, Integer *iblock)
{
  wnga_release_block(*g_a, *iblock);
}

void FATR nga_release_block_grid_(Integer *g_a, Integer *index)
{
  wnga_release_block_grid(*g_a, index);
}

void FATR nga_release_block_segment_(Integer *g_a, Integer *iproc)
{
  wnga_release_block_segment(*g_a, *iproc);
}

void FATR ga_release_update_(Integer *g_a, Integer *ilo, Integer *ihi,
                             Integer *jlo, Integer *jhi)
{
  Integer lo[2], hi[2];
  lo[0]=*ilo;
  lo[1]=*jlo;
  hi[0]=*ihi;
  hi[1]=*jhi;
  wnga_release_update(*g_a, lo, hi);
}

void FATR nga_release_update_(Integer *g_a, Integer *lo, Integer *hi)
{
  wnga_release_update(*g_a, lo, hi);
}

void FATR nga_release_update_block_(Integer *g_a, Integer *iblock)
{
  wnga_release_update_block(*g_a, *iblock);
}

void FATR nga_release_update_block_grid_(Integer *g_a, Integer *index)
{
  wnga_release_update_block_grid(*g_a, index);
}

void FATR nga_release_update_block_segment_(Integer *g_a, Integer *iproc)
{
  wnga_release_update_block_segment(*g_a, *iproc);
}

void FATR ga_scatter_(Integer *g_a, void *v, Integer *i, Integer *j, Integer *nv)
{
  wnga_scatter2d(*g_a, v, i, j, *nv);
}

void FATR nga_scatter_(Integer *g_a, void* v, Integer subscript[], Integer *nv)
{
  wnga_scatter(*g_a, v, subscript, 0, *nv);
}

void FATR ga_scatter_acc_(Integer *g_a, void *v, Integer *i, Integer *j,
                          Integer *nv, void *alpha)
{
  wnga_scatter_acc2d(*g_a, v, i, j, *nv, alpha);
}

void FATR nga_scatter_acc_(Integer *g_a, void* v, Integer subscript[],
                           Integer *nv, void *alpha)
{
  wnga_scatter_acc(*g_a, v, subscript, 0, *nv, alpha);
}

void FATR nga_strided_acc_(Integer *g_a, Integer *lo, Integer *hi,
                           Integer *skip, void *buf, Integer *ld, void *alpha)
{
  wnga_strided_acc(*g_a, lo, hi, skip, buf, ld, alpha);
}

void FATR nga_strided_get_(Integer *g_a, Integer *lo, Integer *hi,
                           Integer *skip, void *buf, Integer *ld)
{
  wnga_strided_get(*g_a, lo, hi, skip, buf, ld);
}

void FATR nga_strided_put_(Integer *g_a, Integer *lo, Integer *hi,
                           Integer *skip, void *buf, Integer *ld)
{
  wnga_strided_put(*g_a, lo, hi, skip, buf, ld);
}

void FATR ga_sync_()
{
  wnga_sync();
}

void FATR nga_sync_()
{
  wnga_sync();
}

/* Routines from global.util.c */

void FATR ga_print_stats_()
{
    wnga_print_stats();
}

void FATR nga_print_stats_()
{
    wnga_print_stats();
}

void FATR ga_error_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
    char *string, Integer *icode, int slen
#else
    char *string, int slen, Integer *icode
#endif
    )
{
  char buf[FMSG];
  ga_f2cstring(string,slen, buf, FMSG);
  wnga_error(buf,*icode);
}

void FATR nga_error_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
    char *string, Integer *icode, int slen
#else
    char *string, int slen, Integer *icode
#endif
    )
{
  char buf[FMSG];
  ga_f2cstring(string,slen, buf, FMSG);
  wnga_error(buf,*icode);
}

Integer FATR ga_cluster_nodeid_()
{
    return wnga_cluster_nodeid();
}

Integer FATR nga_cluster_nodeid_()
{
    return wnga_cluster_nodeid();
}

Integer FATR ga_cluster_nprocs_(Integer *node)
{
    return wnga_cluster_nprocs(*node);
}

Integer FATR nga_cluster_nprocs_(Integer *node)
{
    return wnga_cluster_nprocs(*node);
}

Integer FATR ga_cluster_procid_(Integer *node, Integer *loc_proc_id)
{
    return wnga_cluster_procid(*node, *loc_proc_id);
}

Integer FATR nga_cluster_procid_(Integer *node, Integer *loc_proc_id)
{
    return wnga_cluster_procid(*node, *loc_proc_id);
}

Integer FATR ga_cluster_nnodes_()
{
    return wnga_cluster_nnodes();
}

Integer FATR nga_cluster_nnodes_()
{
    return wnga_cluster_nnodes();
}

Integer FATR ga_cluster_proc_nodeid_(Integer *proc)
{
    return wnga_cluster_proc_nodeid(*proc);
}

Integer FATR nga_cluster_proc_nodeid_(Integer *proc)
{
    return wnga_cluster_proc_nodeid(*proc);
}

void FATR ga_print_distribution_(Integer* g_a)
{
    wnga_print_distribution(1, *g_a);
}

void FATR nga_print_distribution_(Integer* g_a)
{
    wnga_print_distribution(1, *g_a);
}

DoublePrecision FATR ga_wtime_()
{
    return wnga_wtime();
}

DoublePrecision FATR nga_wtime_()
{
    return wnga_wtime();
}

/* Routines from collect.c */

void FATR ga_brdcst_(
        Integer *type, void *buf, Integer *len, Integer *originator)
{
    wnga_brdcst(*type, buf, *len, *originator);
}

void FATR nga_brdcst_(
        Integer *type, void *buf, Integer *len, Integer *originator)
{
    wnga_brdcst(*type, buf, *len, *originator);
}

void FATR ga_pgroup_brdcst_(Integer *grp_id, Integer *type, void *buf, Integer *len, Integer *originator)
{
    wnga_pgroup_brdcst(*grp_id, *type, buf, *len, *originator);
}

void FATR ga_pgroup_gop_(Integer *grp, Integer *type, void *x, Integer *n, char *op, int len)
{
    wnga_pgroup_gop(*grp, pnga_type_f2c(*type), x, *n, op);
}

void FATR nga_pgroup_gop_(Integer *grp, Integer *type, void *x, Integer *n, char *op, int len)
{
    wnga_pgroup_gop(*grp, pnga_type_f2c(*type), x, *n, op);
}

void FATR ga_pgroup_igop_(Integer *grp, Integer *type, Integer *x, Integer *n, char *op, int len)
{
    wnga_pgroup_gop(*grp, pnga_type_f2c(MT_F_INT), x, *n, op);
}

void FATR nga_pgroup_igop_(Integer *grp, Integer *type, Integer *x, Integer *n, char *op, int len)
{
    wnga_pgroup_gop(*grp, pnga_type_f2c(MT_F_INT), x, *n, op);
}

void FATR ga_pgroup_sgop_(Integer *grp, Integer *type, Real *x, Integer *n, char *op, int len)
{
    wnga_pgroup_gop(*grp, pnga_type_f2c(MT_F_REAL), x, *n, op);
}

void FATR nga_pgroup_sgop_(Integer *grp, Integer *type, Real *x, Integer *n, char *op, int len)
{
    wnga_pgroup_gop(*grp, pnga_type_f2c(MT_F_REAL), x, *n, op);
}

void FATR ga_pgroup_dgop_(Integer *grp, Integer *type, DoublePrecision *x, Integer *n, char *op, int len)
{
    wnga_pgroup_gop(*grp, pnga_type_f2c(MT_F_DBL), x, *n, op);
}

void FATR nga_pgroup_dgop_(Integer *grp, Integer *type, DoublePrecision *x, Integer *n, char *op, int len)
{
    wnga_pgroup_gop(*grp, pnga_type_f2c(MT_F_DBL), x, *n, op);
}

void FATR ga_pgroup_cgop_(Integer *grp, Integer *type, SingleComplex *x, Integer *n, char *op, int len)
{
    wnga_pgroup_gop(*grp, pnga_type_f2c(MT_F_SCPL), x, *n, op);
}

void FATR nga_pgroup_cgop_(Integer *grp, Integer *type, SingleComplex *x, Integer *n, char *op, int len)
{
    wnga_pgroup_gop(*grp, pnga_type_f2c(MT_F_SCPL), x, *n, op);
}

void FATR ga_pgroup_zgop_(Integer *grp, Integer *type, DoubleComplex *x, Integer *n, char *op, int len)
{
    wnga_pgroup_gop(*grp, pnga_type_f2c(MT_F_DCPL), x, *n, op);
}

void FATR nga_pgroup_zgop_(Integer *grp, Integer *type, DoubleComplex *x, Integer *n, char *op, int len)
{
    wnga_pgroup_gop(*grp, pnga_type_f2c(MT_F_DCPL), x, *n, op);
}

void FATR ga_gop_(Integer *type, void *x, Integer *n, char *op, int len)
{
    wnga_gop(pnga_type_f2c(*type), x, *n, op);
}

void FATR nga_gop_(Integer *type, void *x, Integer *n, char *op, int len)
{
    wnga_gop(pnga_type_f2c(*type), x, *n, op);
}

void FATR ga_igop_(Integer *type, Integer *x, Integer *n, char *op, int len)
{
    wnga_gop(pnga_type_f2c(MT_F_INT), x, *n, op);
}

void FATR nga_igop_(Integer *type, Integer *x, Integer *n, char *op, int len)
{
    wnga_gop(pnga_type_f2c(MT_F_INT), x, *n, op);
}

void FATR ga_sgop_(Integer *type, Real *x, Integer *n, char *op, int len)
{
    wnga_gop(pnga_type_f2c(MT_F_REAL), x, *n, op);
}

void FATR nga_sgop_(Integer *type, Real *x, Integer *n, char *op, int len)
{
    wnga_gop(pnga_type_f2c(MT_F_REAL), x, *n, op);
}

void FATR ga_dgop_(Integer *type, DoublePrecision *x, Integer *n, char *op, int len)
{
    wnga_gop(pnga_type_f2c(MT_F_DBL), x, *n, op);
}

void FATR nga_dgop_(Integer *type, DoublePrecision *x, Integer *n, char *op, int len)
{
    wnga_gop(pnga_type_f2c(MT_F_DBL), x, *n, op);
}

void FATR ga_cgop_(Integer *type, SingleComplex *x, Integer *n, char *op, int len)
{
    wnga_gop(pnga_type_f2c(MT_F_SCPL), x, *n, op);
}

void FATR nga_cgop_(Integer *type, SingleComplex *x, Integer *n, char *op, int len)
{
    wnga_gop(pnga_type_f2c(MT_F_SCPL), x, *n, op);
}

void FATR ga_zgop_(Integer *type, DoubleComplex *x, Integer *n, char *op, int len)
{
    wnga_gop(pnga_type_f2c(MT_F_DCPL), x, *n, op);
}

void FATR nga_zgop_(Integer *type, DoubleComplex *x, Integer *n, char *op, int len)
{
    wnga_gop(pnga_type_f2c(MT_F_DCPL), x, *n, op);
}

#ifdef MPI
#   include "ga-mpi.h"
#   define ga_mpi_comm_ F77_FUNC_(ga_mpi_comm,GA_MPI_COMM)
void FATR ga_mpi_comm_(Integer *fcomm)
{
    MPI_Comm ccomm;
    GA_MPI_Comm(&ccomm);
    *fcomm = (Integer)(MPI_Comm_c2f(ccomm));
}
#define ga_mpi_comm_pgroup_ F77_FUNC_(ga_mpi_comm_pgroup,GA_MPI_COMM_PGROUP)
void FATR ga_mpi_comm_pgroup_(Integer *fcomm, Integer *pgroup)
{
    MPI_Comm ccomm;
    ccomm = GA_MPI_Comm_pgroup((int)(*pgroup));
    *fcomm = (Integer)(MPI_Comm_c2f(ccomm));
}
#define ga_mpi_comm_pgroup_default_ F77_FUNC_(ga_mpi_comm_pgroup_default,GA_MPI_COMM_PGROUP_DEFAULT)
void FATR ga_mpi_comm_pgroup_default_(Integer *fcomm)
{
    MPI_Comm ccomm;
    ccomm = GA_MPI_Comm_pgroup_default();
    *fcomm = (Integer)(MPI_Comm_c2f(ccomm));
}
#endif

/* Routines from elem_alg.c */

void FATR ga_abs_value_patch_(Integer *g_a, Integer *lo, Integer *hi)
{
    wnga_abs_value_patch(*g_a, lo, hi);
}

void FATR nga_abs_value_patch_(Integer *g_a, Integer *lo, Integer *hi)
{
    wnga_abs_value_patch(*g_a, lo, hi);
}

void FATR ga_recip_patch_(Integer *g_a, Integer *lo, Integer *hi)
{
    wnga_recip_patch(*g_a, lo, hi);
}

void FATR nga_recip_patch_(Integer *g_a, Integer *lo, Integer *hi)
{
    wnga_recip_patch(*g_a, lo, hi);
}

void FATR ga_add_constant_patch_(Integer *g_a, Integer *lo, Integer *hi, void *alpha)
{
    wnga_add_constant_patch(*g_a, lo, hi, alpha);
}

void FATR nga_add_constant_patch_(Integer *g_a, Integer *lo, Integer *hi, void *alpha)
{
    wnga_add_constant_patch(*g_a, lo, hi, alpha);
}

void FATR ga_abs_value_(Integer *g_a)
{
    wnga_abs_value(*g_a);
}

void FATR nga_abs_value_(Integer *g_a)
{
    wnga_abs_value(*g_a);
}

void FATR ga_add_constant_(Integer *g_a, void *alpha)
{
    wnga_add_constant(*g_a, alpha);
}

void FATR nga_add_constant_(Integer *g_a, void *alpha)
{
    wnga_add_constant(*g_a, alpha);
}

void FATR ga_recip_(Integer *g_a)
{
    wnga_recip(*g_a);
}

void FATR ga_elem_multiply_(Integer *g_a, Integer *g_b, Integer *g_c)
{
    wnga_elem_multiply(*g_a, *g_b, *g_c);
}

void FATR nga_elem_multiply_(Integer *g_a, Integer *g_b, Integer *g_c)
{
    wnga_elem_multiply(*g_a, *g_b, *g_c);
}

void FATR ga_elem_divide_(Integer *g_a, Integer *g_b, Integer *g_c)
{
    wnga_elem_divide(*g_a, *g_b, *g_c);
}

void FATR nga_elem_divide_(Integer *g_a, Integer *g_b, Integer *g_c)
{
    wnga_elem_divide(*g_a, *g_b, *g_c);
}

void FATR ga_elem_maximum_(Integer *g_a, Integer *g_b, Integer *g_c)
{
    wnga_elem_maximum(*g_a, *g_b, *g_c);
}

void FATR nga_elem_maximum_(Integer *g_a, Integer *g_b, Integer *g_c)
{
    wnga_elem_maximum(*g_a, *g_b, *g_c);
}

void FATR ga_elem_minimum_(Integer *g_a, Integer *g_b, Integer *g_c)
{
    wnga_elem_minimum(*g_a, *g_b, *g_c);
}

void FATR nga_elem_minimum_(Integer *g_a, Integer *g_b, Integer *g_c)
{
    wnga_elem_minimum(*g_a, *g_b, *g_c);
}

void FATR ga_elem_multiply_patch_(Integer *g_a,Integer *alo,Integer *ahi,Integer *g_b,Integer *blo,Integer *bhi,Integer *g_c,Integer *clo,Integer *chi)
{
    wnga_elem_multiply_patch(*g_a,alo,ahi,*g_b,blo,bhi,*g_c,clo,chi);
}

void FATR nga_elem_multiply_patch_(Integer *g_a,Integer *alo,Integer *ahi,Integer *g_b,Integer *blo,Integer *bhi,Integer *g_c,Integer *clo,Integer *chi)
{
    wnga_elem_multiply_patch(*g_a,alo,ahi,*g_b,blo,bhi,*g_c,clo,chi);
}

void FATR ga_elem_divide_patch_(Integer *g_a,Integer *alo,Integer *ahi, Integer *g_b,Integer *blo,Integer *bhi,Integer *g_c, Integer *clo,Integer *chi)
{
    wnga_elem_divide_patch(*g_a,alo,ahi, *g_b,blo,bhi,*g_c, clo,chi);
}

void FATR nga_elem_divide_patch_(Integer *g_a,Integer *alo,Integer *ahi, Integer *g_b,Integer *blo,Integer *bhi,Integer *g_c, Integer *clo,Integer *chi)
{
    wnga_elem_divide_patch(*g_a,alo,ahi, *g_b,blo,bhi,*g_c, clo,chi);
}

void FATR ga_elem_step_divide_patch_(Integer *g_a,Integer *alo,Integer *ahi, Integer *g_b,Integer *blo,Integer *bhi,Integer *g_c, Integer *clo,Integer *chi)
{
    wnga_elem_step_divide_patch(*g_a,alo,ahi, *g_b,blo,bhi,*g_c, clo,chi);
}

void FATR nga_elem_step_divide_patch_(Integer *g_a,Integer *alo,Integer *ahi, Integer *g_b,Integer *blo,Integer *bhi,Integer *g_c, Integer *clo,Integer *chi)
{
    wnga_elem_step_divide_patch(*g_a,alo,ahi, *g_b,blo,bhi,*g_c, clo,chi);
}

void FATR ga_elem_stepb_divide_patch_(Integer *g_a,Integer *alo,Integer *ahi, Integer *g_b,Integer *blo,Integer *bhi,Integer *g_c, Integer *clo,Integer *chi)
{
    wnga_elem_stepb_divide_patch(*g_a,alo,ahi, *g_b,blo,bhi,*g_c, clo,chi);
}

void FATR nga_elem_stepb_divide_patch_(Integer *g_a,Integer *alo,Integer *ahi, Integer *g_b,Integer *blo,Integer *bhi,Integer *g_c, Integer *clo,Integer *chi)
{
    wnga_elem_stepb_divide_patch(*g_a,alo,ahi, *g_b,blo,bhi,*g_c, clo,chi);
}

void FATR ga_elem_maximum_patch_(Integer *g_a,Integer *alo,Integer *ahi, Integer *g_b,Integer *blo,Integer *bhi,Integer *g_c,Integer *clo,Integer *chi)
{
    wnga_elem_maximum_patch(*g_a,alo,ahi, *g_b,blo,bhi,*g_c,clo,chi);
}

void FATR nga_elem_maximum_patch_(Integer *g_a,Integer *alo,Integer *ahi, Integer *g_b,Integer *blo,Integer *bhi,Integer *g_c,Integer *clo,Integer *chi)
{
    wnga_elem_maximum_patch(*g_a,alo,ahi, *g_b,blo,bhi,*g_c,clo,chi);
}

void FATR ga_elem_minimum_patch_(Integer *g_a,Integer *alo,Integer *ahi, Integer *g_b,Integer *blo,Integer *bhi,Integer *g_c,Integer *clo,Integer *chi)
{
    wnga_elem_minimum_patch(*g_a,alo,ahi, *g_b,blo,bhi,*g_c,clo,chi);
}

void FATR nga_elem_minimum_patch_(Integer *g_a,Integer *alo,Integer *ahi, Integer *g_b,Integer *blo,Integer *bhi,Integer *g_c,Integer *clo,Integer *chi)
{
    wnga_elem_minimum_patch(*g_a,alo,ahi, *g_b,blo,bhi,*g_c,clo,chi);
}

void FATR ga_step_bound_info_patch_(Integer *g_xx, Integer *xxlo, Integer *xxhi, Integer *g_vv, Integer *vvlo, Integer *vvhi, Integer *g_xxll, Integer *xxlllo, Integer *xxllhi, Integer *g_xxuu, Integer *xxuulo, Integer *xxuuhi, void *boundmin, void* wolfemin, void *boundmax)
{
    wnga_step_bound_info_patch(*g_xx, xxlo, xxhi, *g_vv, vvlo, vvhi, *g_xxll, xxlllo, xxllhi, *g_xxuu, xxuulo, xxuuhi, boundmin, wolfemin, boundmax);
}

void FATR nga_step_bound_info_patch_(Integer *g_xx, Integer *xxlo, Integer *xxhi, Integer *g_vv, Integer *vvlo, Integer *vvhi, Integer *g_xxll, Integer *xxlllo, Integer *xxllhi, Integer *g_xxuu, Integer *xxuulo, Integer *xxuuhi, void *boundmin, void* wolfemin, void *boundmax)
{
    wnga_step_bound_info_patch(*g_xx, xxlo, xxhi, *g_vv, vvlo, vvhi, *g_xxll, xxlllo, xxllhi, *g_xxuu, xxuulo, xxuuhi, boundmin, wolfemin, boundmax);
}

void FATR ga_step_max_patch_(Integer *g_a,  Integer *alo, Integer *ahi, Integer *g_b,  Integer *blo, Integer *bhi, void *result)
{
    wnga_step_max_patch(*g_a, alo, ahi, *g_b, blo, bhi, result);
}

void FATR nga_step_max_patch_(Integer *g_a,  Integer *alo, Integer *ahi, Integer *g_b,  Integer *blo, Integer *bhi, void *result)
{
    wnga_step_max_patch(*g_a, alo, ahi, *g_b, blo, bhi, result);
}

void FATR ga_step_max_(Integer *g_a, Integer *g_b, void *retval)
{
    wnga_step_max(*g_a, *g_b, retval);
}

void FATR nga_step_max_(Integer *g_a, Integer *g_b, void *retval)
{
    wnga_step_max(*g_a, *g_b, retval);
}

void FATR ga_step_bound_info_(Integer *g_xx, Integer *g_vv, Integer *g_xxll, Integer *g_xxuu,  void *boundmin, void *wolfemin, void *boundmax)
{
    wnga_step_bound_info(*g_xx, *g_vv, *g_xxll, *g_xxuu, boundmin, wolfemin, boundmax);
}

void FATR nga_step_bound_info_(Integer *g_xx, Integer *g_vv, Integer *g_xxll, Integer *g_xxuu,  void *boundmin, void *wolfemin, void *boundmax)
{
    wnga_step_bound_info(*g_xx, *g_vv, *g_xxll, *g_xxuu, boundmin, wolfemin, boundmax);
}

/* Routines from ga_solve_seq.c */

void FATR ga_lu_solve_seq_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *trans, Integer *g_a, Integer *g_b, int len
#else
        char *trans, int len, Integer *g_a, Integer *g_b
#endif
        )
{
    char buf[FNAM];
    ga_f2cstring(trans, len, buf, FNAM);
    wnga_lu_solve_seq(buf, *g_a, *g_b);
}

/* Routines from global.util.c */

void FATR ga_print_(Integer *g_a)
{
    wnga_print(*g_a);
}

void FATR nga_print_(Integer *g_a)
{
    wnga_print(*g_a);
}

void FATR ga_print_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, Integer *pretty)
{
    wnga_print_patch2d(*g_a, *ilo, *ihi, *jlo, *jhi, *pretty);
}

void FATR nga_print_patch_(Integer *g_a, Integer *lo, Integer *hi, Integer *pretty)
{
    wnga_print_patch(*g_a, lo, hi, *pretty);
}

void FATR ga_summarize_(Integer *verbose)
{
    wnga_summarize(*verbose);
}

void FATR nga_summarize_(Integer *verbose)
{
    wnga_summarize(*verbose);
}

/* Routines from ghosts.c */

void FATR nga_access_ghost_element_(Integer* g_a, AccessIndex* index, Integer subscript[], Integer ld[])
{
    wnga_access_ghost_element(*g_a, index, subscript, ld);
}

void FATR nga_access_ghosts_(Integer* g_a, Integer dims[], AccessIndex* index, Integer ld[])
{
    wnga_access_ghosts(*g_a, dims, index, ld);
}

void FATR nga_release_ghost_element_(Integer* g_a, Integer subscript[])
{
    wnga_release_ghost_element(*g_a, subscript);
}

void FATR nga_release_update_ghost_element_(Integer* g_a, Integer subscript[])
{
    wnga_release_update_ghost_element(*g_a, subscript);
}

void FATR nga_release_ghosts_(Integer* g_a)        
{
    wnga_release_ghosts(*g_a);
}

void FATR nga_release_update_ghosts_(Integer* g_a)
{
    wnga_release_update_ghosts(*g_a);
}

void FATR nga_get_ghost_block_(Integer *g_a, Integer *lo, Integer *hi, void *buf, Integer *ld) 
{
    wnga_get_ghost_block(*g_a, lo, hi, buf, ld);
}

void FATR ga_update1_ghosts_(Integer *g_a)
{
    wnga_update1_ghosts(*g_a);
}

void FATR nga_update1_ghosts_(Integer *g_a)
{
    wnga_update1_ghosts(*g_a);
}

logical FATR ga_update2_ghosts_(Integer *g_a)
{
    return wnga_update2_ghosts(*g_a);
}

logical FATR nga_update2_ghosts_(Integer *g_a)
{
    return wnga_update2_ghosts(*g_a);
}

logical FATR ga_update3_ghosts_(Integer *g_a)
{
    return wnga_update3_ghosts(*g_a);
}

logical FATR nga_update3_ghosts_(Integer *g_a)
{
    return wnga_update3_ghosts(*g_a);
}

logical FATR ga_set_update4_info_(Integer *g_a)
{
    return wnga_set_update4_info(*g_a);
}

logical FATR nga_set_update4_info_(Integer *g_a)
{
    return wnga_set_update4_info(*g_a);
}

logical FATR ga_update4_ghosts_(Integer *g_a)
{
    return wnga_update4_ghosts(*g_a);
}

logical FATR nga_update4_ghosts_(Integer *g_a)
{
    return wnga_update4_ghosts(*g_a);
}

logical FATR ga_update44_ghosts_(Integer *g_a)
{
    return wnga_update44_ghosts(*g_a);
}

logical FATR nga_update44_ghosts_(Integer *g_a)
{
    return wnga_update44_ghosts(*g_a);
}

logical FATR ga_update55_ghosts_(Integer *g_a)
{
    return wnga_update55_ghosts(*g_a);
}

logical FATR nga_update55_ghosts_(Integer *g_a)
{
    return wnga_update55_ghosts(*g_a);
}

logical FATR nga_update_ghost_dir_(Integer *g_a, Integer *pdim, Integer *pdir, logical *pflag)
{
    return wnga_update_ghost_dir(*g_a, *pdim, *pdir, *pflag);
}

logical FATR ga_update5_ghosts_(Integer *g_a)
{
    return wnga_update5_ghosts(*g_a);
}

logical FATR nga_update5_ghosts_(Integer *g_a)
{
    return wnga_update5_ghosts(*g_a);
}

logical FATR ga_set_update5_info_(Integer *g_a)
{
    return wnga_set_update5_info(*g_a);
}

logical FATR nga_set_update5_info_(Integer *g_a)
{
    return wnga_set_update5_info(*g_a);
}

void FATR ga_update_ghosts_(Integer *g_a)
{
    wnga_update_ghosts(*g_a);
}

void FATR nga_update_ghosts_(Integer *g_a)
{
    wnga_update_ghosts(*g_a);
}

void FATR nga_update_ghosts_nb_(Integer *g_a, Integer *nb)
{
    wnga_update_ghosts_nb(*g_a, nb);
}

logical FATR ga_update6_ghosts_(Integer *g_a)
{
    return wnga_update6_ghosts(*g_a);
}

logical FATR nga_update6_ghosts_(Integer *g_a)
{
    return wnga_update6_ghosts(*g_a);
}

logical FATR ga_update7_ghosts_(Integer *g_a)
{
    return wnga_update7_ghosts(*g_a);
}

logical FATR nga_update7_ghosts_(Integer *g_a)
{
    return wnga_update7_ghosts(*g_a);
}

void FATR ga_ghost_barrier_()
{
    wnga_ghost_barrier();
}

void FATR nga_ghost_barrier_()
{
    wnga_ghost_barrier();
}

void FATR nga_nbget_ghost_dir_(Integer *g_a, Integer *mask, Integer *nbhandle)
{
    wnga_nbget_ghost_dir(*g_a, mask, nbhandle);
}

void FATR ga_set_ghost_corner_flag_(Integer *g_a, logical *flag)
{
    wnga_set_ghost_corner_flag(*g_a, *flag);
}

void FATR nga_set_ghost_corner_flag_(Integer *g_a, logical *flag)
{
    wnga_set_ghost_corner_flag(*g_a, *flag);
}

logical FATR ga_set_ghost_info_(Integer *g_a)
{
    return wnga_set_ghost_info(*g_a);
}

logical FATR nga_set_ghost_info_(Integer *g_a)
{
    return wnga_set_ghost_info(*g_a);
}

/* Routines from global.nalg.c */

void FATR ga_zero_(Integer *g_a)
{
    wnga_zero(*g_a);
}

void FATR nga_zero_(Integer *g_a)
{
    wnga_zero(*g_a);
}

void FATR ga_copy_(Integer *g_a, Integer *g_b)
{
    wnga_copy(*g_a, *g_b);
}

void FATR nga_copy_(Integer *g_a, Integer *g_b)
{
    wnga_copy(*g_a, *g_b);
}

Integer FATR ga_idot_(Integer *g_a, Integer *g_b)
{
    Integer sum;
    wnga_dot(pnga_type_f2c(MT_F_INT), *g_a, *g_b, &sum);
    return sum;
}

Integer FATR nga_idot_(Integer *g_a, Integer *g_b)
{
    Integer sum;
    wnga_dot(pnga_type_f2c(MT_F_INT), *g_a, *g_b, &sum);
    return sum;
}

DoublePrecision FATR ga_ddot_(Integer *g_a, Integer *g_b)
{
    DoublePrecision sum;
    wnga_dot(pnga_type_f2c(MT_F_DBL), *g_a, *g_b, &sum);
    return sum;
}

DoublePrecision FATR nga_ddot_(Integer *g_a, Integer *g_b)
{
    DoublePrecision sum;
    wnga_dot(pnga_type_f2c(MT_F_DBL), *g_a, *g_b, &sum);
    return sum;
}

Real FATR ga_sdot_(Integer *g_a, Integer *g_b)
{
    Real sum;
    wnga_dot(pnga_type_f2c(MT_F_REAL), *g_a, *g_b, &sum);
    return sum;
}            

Real FATR nga_sdot_(Integer *g_a, Integer *g_b)
{
    Real sum;
    wnga_dot(pnga_type_f2c(MT_F_REAL), *g_a, *g_b, &sum);
    return sum;
}            

void FATR gai_zdot_(Integer *g_a, Integer *g_b, DoubleComplex *sum)
{
    wnga_dot(pnga_type_f2c(MT_F_DCPL), *g_a, *g_b, sum);
}

void FATR ngai_zdot_(Integer *g_a, Integer *g_b, DoubleComplex *sum)
{
    wnga_dot(pnga_type_f2c(MT_F_DCPL), *g_a, *g_b, sum);
}

void gai_cdot_(Integer *g_a, Integer *g_b, SingleComplex *sum)
{
    wnga_dot(pnga_type_f2c(MT_F_SCPL), *g_a, *g_b, sum);
}

void ngai_cdot_(Integer *g_a, Integer *g_b, SingleComplex *sum)
{
    wnga_dot(pnga_type_f2c(MT_F_SCPL), *g_a, *g_b, sum);
}

void FATR ga_scale_(Integer *g_a, void* alpha)
{
    wnga_scale(*g_a, alpha);
}

void FATR ga_cscale_(Integer *g_a, SingleComplex* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_SCPL) pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR ga_cscal_(Integer *g_a, SingleComplex* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_SCPL) pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR ga_dscale_(Integer *g_a, DoublePrecision* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_DBL) pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR ga_dscal_(Integer *g_a, DoublePrecision* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_DBL) pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR ga_iscale_(Integer *g_a, Integer* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_INT || atype != C_LONG || atype != C_LONGLONG)
        pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR ga_iscal_(Integer *g_a, Integer* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_INT || atype != C_LONG || atype != C_LONGLONG)
        pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR ga_sscale_(Integer *g_a, Real* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_FLOAT) pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR ga_sscal_(Integer *g_a, Real* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_FLOAT) pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR ga_zscale_(Integer *g_a, DoubleComplex* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_DCPL) pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR ga_zscal_(Integer *g_a, DoubleComplex* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_DCPL) pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR nga_scale_(Integer *g_a, void* alpha)
{
    wnga_scale(*g_a, alpha);
}

void FATR nga_cscale_(Integer *g_a, SingleComplex* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_SCPL) pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR nga_dscale_(Integer *g_a, DoublePrecision* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_DBL) pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR nga_iscale_(Integer *g_a, Integer* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_INT || atype != C_LONG || atype != C_LONGLONG)
        pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR nga_sscale_(Integer *g_a, Real* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_FLOAT) pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR nga_zscale_(Integer *g_a, DoubleComplex* alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_DCPL) pnga_error(" wrong type ", 0L);
    wnga_scale(*g_a, alpha);
}

void FATR ga_add_(void *alpha, Integer* g_a, void* beta, Integer* g_b, Integer* g_c)
{
    wnga_add(alpha, *g_a, beta, *g_b, *g_c);
}

void FATR ga_cadd_(SingleComplex *alpha, Integer *g_a, SingleComplex *beta,
        Integer *g_b, Integer *g_c)
{
    wnga_add(alpha, *g_a, beta, *g_b, *g_c);
}

void FATR ga_dadd_(DoublePrecision *alpha, Integer *g_a, DoublePrecision *beta,
        Integer *g_b, Integer *g_c)
{
    wnga_add(alpha, *g_a, beta, *g_b, *g_c);
}

void FATR ga_iadd_(Integer *alpha, Integer *g_a, Integer *beta,
        Integer *g_b, Integer *g_c)
{
    wnga_add(alpha, *g_a, beta, *g_b, *g_c);
}

void FATR ga_sadd_(Real *alpha, Integer *g_a, Real *beta,
        Integer *g_b, Integer *g_c)
{
    wnga_add(alpha, *g_a, beta, *g_b, *g_c);
}

void FATR ga_zadd_(DoubleComplex *alpha, Integer *g_a, DoubleComplex *beta,
        Integer *g_b, Integer *g_c)
{
    wnga_add(alpha, *g_a, beta, *g_b, *g_c);
}

void FATR nga_add_(void *alpha, Integer* g_a, void* beta, Integer* g_b, Integer* g_c)
{
    wnga_add(alpha, *g_a, beta, *g_b, *g_c);
}

void FATR nga_cadd_(SingleComplex *alpha, Integer *g_a, SingleComplex *beta,
        Integer *g_b, Integer *g_c)
{
    wnga_add(alpha, *g_a, beta, *g_b, *g_c);
}

void FATR nga_dadd_(DoublePrecision *alpha, Integer *g_a, DoublePrecision *beta,
        Integer *g_b, Integer *g_c)
{
    wnga_add(alpha, *g_a, beta, *g_b, *g_c);
}

void FATR nga_iadd_(Integer *alpha, Integer *g_a, Integer *beta,
        Integer *g_b, Integer *g_c)
{
    wnga_add(alpha, *g_a, beta, *g_b, *g_c);
}

void FATR nga_sadd_(Real *alpha, Integer *g_a, Real *beta,
        Integer *g_b, Integer *g_c)
{
    wnga_add(alpha, *g_a, beta, *g_b, *g_c);
}

void FATR nga_zadd_(DoubleComplex *alpha, Integer *g_a, DoubleComplex *beta,
        Integer *g_b, Integer *g_c)
{
    wnga_add(alpha, *g_a, beta, *g_b, *g_c);
}

void FATR ga_transpose_(Integer *g_a, Integer *g_b)
{
    wnga_transpose(*g_a, *g_b);
}

void FATR nga_transpose_(Integer *g_a, Integer *g_b)
{
    wnga_transpose(*g_a, *g_b);
}

/* Routines from global.npatch.c */

void FATR ga_copy_patch_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *trans, Integer *g_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, int len
#else
        char *trans, int len, Integer *g_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi
#endif
        )
{
    Integer alo[2], ahi[2], blo[2], bhi[2];

    alo[0] = *ailo; alo[1] = *ajlo;
    ahi[0] = *aihi; ahi[1] = *ajhi;
    blo[0] = *bilo; blo[1] = *bjlo;
    bhi[0] = *bihi; bhi[1] = *bjhi;

    wnga_copy_patch(trans, *g_a, alo, ahi, *g_b, blo, bhi);
}

void FATR nga_copy_patch_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *trans, Integer *g_a, Integer *alo, Integer *ahi, Integer *g_b, Integer *blo, Integer *bhi, int len
#else
        char *trans, int len, Integer *g_a, Integer *alo, Integer *ahi, Integer *g_b, Integer *blo, Integer *bhi
#endif
        )
{
    wnga_copy_patch(trans,*g_a,alo,ahi,*g_b,blo,bhi);
}

void FATR ga_zero_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi)
{
    Integer lo[2], hi[2];

    lo[0] = *ilo; lo[1] = *jlo;
    hi[0] = *ihi; lo[1] = *jhi;
    wnga_zero_patch(*g_a, lo, hi);
}

void FATR nga_zero_patch_(Integer *g_a, Integer *lo, Integer *hi)
{
    wnga_zero_patch(*g_a, lo, hi);
}

static void sga_dot_patch(Integer g_a, char *t_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer g_b, char *t_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, void *retval)
{
    Integer alo[2], ahi[2], blo[2], bhi[2];

    alo[0] = *ailo; alo[1] = *ajlo;
    ahi[0] = *aihi; ahi[1] = *ajhi;
    blo[0] = *bilo; blo[1] = *bjlo;
    bhi[0] = *bihi; bhi[1] = *bjhi;
    wnga_dot_patch(g_a, t_a, alo, ahi, g_b, t_b, blo, bhi, retval);
}

#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
DoublePrecision ga_ddot_patch_(Integer *g_a, char *t_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, char *t_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, int alen, int blen)
#else
DoublePrecision ga_ddot_patch_(Integer *g_a, char *t_a, int alen, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, char *t_b, int blen, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi)
#endif
{
    DoublePrecision retval;
    Integer atype, btype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    if (atype != btype || atype != C_DBL) pnga_error(" wrong types ", 0L);
    sga_dot_patch(*g_a, t_a, ailo, aihi, ajlo, ajhi, *g_b, t_b, bilo, bihi, bjlo, bjhi, &retval);

    return retval;
}

#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
Integer ga_idot_patch_(Integer *g_a, char *t_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, char *t_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, int alen, int blen)
#else
Integer ga_idot_patch_(Integer *g_a, char *t_a, int alen, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, char *t_b, int blen, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi)
#endif
{
    Integer retval;
    Integer atype, btype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    if (atype != btype
            || (atype != C_INT && atype != C_LONG && atype != C_LONGLONG))
        pnga_error(" wrong types ", 0L);
    sga_dot_patch(*g_a, t_a, ailo, aihi, ajlo, ajhi, *g_b, t_b, bilo, bihi, bjlo, bjhi, &retval);

    return retval;
}

#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
Real ga_sdot_patch_(Integer *g_a, char *t_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, char *t_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, int alen, int blen)
#else
Real ga_sdot_patch_(Integer *g_a, char *t_a, int alen, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, char *t_b, int blen, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi)
#endif
{
    Real retval;
    Integer atype, btype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    if (atype != btype || atype != C_FLOAT) pnga_error(" wrong types ", 0L);
    sga_dot_patch(*g_a, t_a, ailo, aihi, ajlo, ajhi, *g_b, t_b, bilo, bihi, bjlo, bjhi, &retval);

    return retval;
}

#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
void gai_zdot_patch_(Integer *g_a, char *t_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, char *t_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, DoubleComplex *retval, int alen, int blen)
#else
void gai_zdot_patch_(Integer *g_a, char *t_a, int alen, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, char *t_b, int blen, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, DoubleComplex *retval)
#endif
{
    Integer atype, btype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    if (atype != btype || atype != C_DCPL) pnga_error(" wrong types ", 0L);
    sga_dot_patch(*g_a, t_a, ailo, aihi, ajlo, ajhi, *g_b, t_b, bilo, bihi, bjlo, bjhi, retval);
}

#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
void gai_cdot_patch_(Integer *g_a, char *t_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, char *t_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, SingleComplex *retval, int alen, int blen)
#else
void gai_cdot_patch_(Integer *g_a, char *t_a, int alen, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, char *t_b, int blen, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, SingleComplex *retval)
#endif
{
    Integer atype, btype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    if (atype != btype || atype != C_SCPL) pnga_error(" wrong types ", 0L);
    sga_dot_patch(*g_a, t_a, ailo, aihi, ajlo, ajhi, *g_b, t_b, bilo, bihi, bjlo, bjhi, retval);
}

#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
DoublePrecision nga_ddot_patch_(Integer *g_a, char *t_a, Integer *alo, Integer *ahi, Integer *g_b, char *t_b, Integer *blo, Integer *bhi, int alen, int blen)
#else
DoublePrecision nga_ddot_patch_(Integer *g_a, char *t_a, int alen, Integer *alo, Integer *ahi, Integer *g_b, char *t_b, int blen, Integer *blo, Integer *bhi)
#endif
{
    DoublePrecision retval;
    Integer atype, btype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    if (atype != btype || atype != C_DBL) pnga_error(" wrong types ", 0L);
    wnga_dot_patch(*g_a, t_a, alo, ahi, *g_b, t_b, blo, bhi, &retval);

    return retval;
}

#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
Integer nga_idot_patch_(Integer *g_a, char *t_a, Integer *alo, Integer *ahi, Integer *g_b, char *t_b, Integer *blo, Integer *bhi, int alen, int blen)
#else
Integer nga_idot_patch_(Integer *g_a, char *t_a, int alen, Integer *alo, Integer *ahi, Integer *g_b, char *t_b, int blen, Integer *blo, Integer *bhi)
#endif
{
    Integer retval;
    Integer atype, btype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    if (atype != btype
            || (atype != C_INT && atype != C_LONG && atype != C_LONGLONG))
        pnga_error(" wrong types ", 0L);
    wnga_dot_patch(*g_a, t_a, alo, ahi, *g_b, t_b, blo, bhi, &retval);

    return retval;
}

#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
Real nga_sdot_patch_(Integer *g_a, char *t_a, Integer *alo, Integer *ahi, Integer *g_b, char *t_b, Integer *blo, Integer *bhi, int alen, int blen)
#else
Real nga_sdot_patch_(Integer *g_a, char *t_a, int alen, Integer *alo, Integer *ahi, Integer *g_b, char *t_b, int blen, Integer *blo, Integer *bhi)
#endif
{
    Real retval;
    Integer atype, btype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    if (atype != btype || atype != C_FLOAT) pnga_error(" wrong types ", 0L);
    wnga_dot_patch(*g_a, t_a, alo, ahi, *g_b, t_b, blo, bhi, &retval);

    return retval;
}

#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
void ngai_cdot_patch_(Integer *g_a, char *t_a, Integer *alo, Integer *ahi, Integer *g_b, char *t_b, Integer *blo, Integer *bhi, SingleComplex *retval, int alen, int blen)
#else
void ngai_cdot_patch_(Integer *g_a, char *t_a, int alen, Integer *alo, Integer *ahi, Integer *g_b, char *t_b, int blen, Integer *blo, Integer *bhi, SingleComplex *retval)
#endif
{
    Integer atype, btype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    if (atype != btype || atype != C_SCPL) pnga_error(" wrong types ", 0L);
    wnga_dot_patch(*g_a, t_a, alo, ahi, *g_b, t_b, blo, bhi, retval);
}

#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
void ngai_zdot_patch_(Integer *g_a, char *t_a, Integer *alo, Integer *ahi, Integer *g_b, char *t_b, Integer *blo, Integer *bhi, DoubleComplex *retval, int alen, int blen)
#else
void ngai_zdot_patch_(Integer *g_a, char *t_a, int alen, Integer *alo, Integer *ahi, Integer *g_b, char *t_b, int blen, Integer *blo, Integer *bhi, DoubleComplex *retval)
#endif
{
    Integer atype, btype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    if (atype != btype || atype != C_DCPL) pnga_error(" wrong types ", 0L);
    wnga_dot_patch(*g_a, t_a, alo, ahi, *g_b, t_b, blo, bhi, retval);
}

static void sga_fill_patch(Integer g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, void* val)
{
    Integer lo[2], hi[2];

    lo[0] = *ilo; lo[1] = *jlo;
    hi[0] = *ihi; hi[1] = *jhi;
    wnga_fill_patch(g_a, lo, hi, val);
}

void FATR ga_fill_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, void* val)
{
    sga_fill_patch(*g_a, ilo, ihi, jlo, jhi, val);
}

void FATR ga_cfill_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, SingleComplex *val)
{
    sga_fill_patch(*g_a, ilo, ihi, jlo, jhi, val);
}

void FATR ga_dfill_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, DoublePrecision *val)
{
    sga_fill_patch(*g_a, ilo, ihi, jlo, jhi, val);
}

void FATR ga_ifill_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, Integer *val)
{
    sga_fill_patch(*g_a, ilo, ihi, jlo, jhi, val);
}

void FATR ga_sfill_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, Real *val)
{
    sga_fill_patch(*g_a, ilo, ihi, jlo, jhi, val);
}

void FATR ga_zfill_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, DoubleComplex *val)
{
    sga_fill_patch(*g_a, ilo, ihi, jlo, jhi, val);
}

void FATR nga_fill_patch_(Integer *g_a, Integer *lo, Integer *hi, void* val)
{
    wnga_fill_patch(*g_a, lo, hi, val);
}

void FATR nga_cfill_patch_(Integer *g_a, Integer *lo, Integer *hi, SingleComplex* val)
{
    wnga_fill_patch(*g_a, lo, hi, val);
}

void FATR nga_dfill_patch_(Integer *g_a, Integer *lo, Integer *hi, DoublePrecision* val)
{
    wnga_fill_patch(*g_a, lo, hi, val);
}

void FATR nga_ifill_patch_(Integer *g_a, Integer *lo, Integer *hi, Integer* val)
{
    wnga_fill_patch(*g_a, lo, hi, val);
}

void FATR nga_sfill_patch_(Integer *g_a, Integer *lo, Integer *hi, SingleComplex* val)
{
    wnga_fill_patch(*g_a, lo, hi, val);
}

void FATR nga_zfill_patch_(Integer *g_a, Integer *lo, Integer *hi, DoubleComplex* val)
{
    wnga_fill_patch(*g_a, lo, hi, val);
}

static void sga_scale_patch(Integer g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, void *alpha)
{
    Integer lo[2], hi[2];

    lo[0] = *ilo; lo[1] = *jlo;
    hi[0] = *ihi; hi[1] = *jhi;
    wnga_scale_patch(g_a, lo, hi, alpha);
}

void FATR ga_scale_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, void *alpha)
{
    sga_scale_patch(*g_a, ilo, ihi, jlo, jhi, alpha);
}

void FATR ga_cscal_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, SingleComplex *alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_SCPL) pnga_error(" wrong types ", 0L);
    sga_scale_patch(*g_a, ilo, ihi, jlo, jhi, alpha);
}

void FATR ga_dscal_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, DoublePrecision *alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_DBL) pnga_error(" wrong types ", 0L);
    sga_scale_patch(*g_a, ilo, ihi, jlo, jhi, alpha);
}

void FATR ga_iscal_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, Integer *alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_INT || atype != C_LONG || atype != C_LONGLONG)
        pnga_error(" wrong types ", 0L);
    sga_scale_patch(*g_a, ilo, ihi, jlo, jhi, alpha);
}

void FATR ga_sscal_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, Real *alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_FLOAT) pnga_error(" wrong types ", 0L);
    sga_scale_patch(*g_a, ilo, ihi, jlo, jhi, alpha);
}

void FATR ga_zscal_patch_(Integer *g_a, Integer *ilo, Integer *ihi, Integer *jlo, Integer *jhi, DoubleComplex *alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_DCPL) pnga_error(" wrong types ", 0L);
    sga_scale_patch(*g_a, ilo, ihi, jlo, jhi, alpha);
}

void FATR nga_scale_patch_(Integer *g_a, Integer *lo, Integer *hi, void *alpha)
{
    wnga_scale_patch(*g_a, lo, hi, alpha);
}

void FATR nga_cscale_patch_(Integer *g_a, Integer *lo, Integer *hi, SingleComplex *alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_SCPL) pnga_error(" wrong types ", 0L);
    wnga_scale_patch(*g_a, lo, hi, alpha);
}

void FATR nga_dscale_patch_(Integer *g_a, Integer *lo, Integer *hi, DoublePrecision *alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_DBL) pnga_error(" wrong types ", 0L);
    wnga_scale_patch(*g_a, lo, hi, alpha);
}

void FATR nga_iscale_patch_(Integer *g_a, Integer *lo, Integer *hi, Integer *alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_INT || atype != C_LONG || atype != C_LONGLONG)
        pnga_error(" wrong types ", 0L);
    wnga_scale_patch(*g_a, lo, hi, alpha);
}

void FATR nga_sscale_patch_(Integer *g_a, Integer *lo, Integer *hi, Real *alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_FLOAT) pnga_error(" wrong types ", 0L);
    wnga_scale_patch(*g_a, lo, hi, alpha);
}

void FATR nga_zscale_patch_(Integer *g_a, Integer *lo, Integer *hi, DoubleComplex *alpha)
{
    Integer atype;

    pnga_inquire_type(*g_a, &atype);
    if (atype != C_DCPL) pnga_error(" wrong types ", 0L);
    wnga_scale_patch(*g_a, lo, hi, alpha);
}

static void sga_add_patch(void *alpha, Integer g_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, void *beta, Integer g_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, Integer g_c, Integer *cilo, Integer *cihi, Integer *cjlo, Integer *cjhi)
{
    Integer alo[2], ahi[2], blo[2], bhi[2], clo[2], chi[2];

    alo[0] = *ailo; alo[1] = *ajlo;
    ahi[0] = *aihi; ahi[1] = *ajhi;
    blo[0] = *bilo; blo[1] = *bjlo;
    bhi[0] = *bihi; bhi[1] = *bjhi;
    clo[0] = *cilo; clo[1] = *cjlo;
    chi[0] = *cihi; chi[1] = *cjhi;
    wnga_add_patch(alpha, g_a, alo, ahi, beta, g_b, blo, bhi, g_c, clo, chi);
}

void FATR ga_add_patch_(void *alpha, Integer *g_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, void *beta, Integer *g_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, Integer *g_c, Integer *cilo, Integer *cihi, Integer *cjlo, Integer *cjhi)
{
    sga_add_patch(alpha, *g_a, ailo, aihi, ajlo, ajhi, beta, *g_b, bilo, bihi, bjlo, bjhi, *g_c, cilo, cihi, cjlo, cjhi);
}

void FATR ga_cadd_patch_(SingleComplex *alpha, Integer *g_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, SingleComplex *beta, Integer *g_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, Integer *g_c, Integer *cilo, Integer *cihi, Integer *cjlo, Integer *cjhi)
{
    Integer atype, btype, ctype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    pnga_inquire_type(*g_c, &ctype);
    if (atype != btype || atype != ctype || atype != C_SCPL)
        pnga_error(" wrong types ", 0L);
    sga_add_patch(alpha, *g_a, ailo, aihi, ajlo, ajhi, beta, *g_b, bilo, bihi, bjlo, bjhi, *g_c, cilo, cihi, cjlo, cjhi);
}


void FATR ga_dadd_patch_(DoublePrecision *alpha, Integer *g_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, DoublePrecision *beta, Integer *g_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, Integer *g_c, Integer *cilo, Integer *cihi, Integer *cjlo, Integer *cjhi)
{
    Integer atype, btype, ctype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    pnga_inquire_type(*g_c, &ctype);
    if (atype != btype || atype != ctype || atype != C_DBL)
        pnga_error(" wrong types ", 0L);
    sga_add_patch(alpha, *g_a, ailo, aihi, ajlo, ajhi, beta, *g_b, bilo, bihi, bjlo, bjhi, *g_c, cilo, cihi, cjlo, cjhi);
}


void FATR ga_iadd_patch_(Integer *alpha, Integer *g_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *beta, Integer *g_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, Integer *g_c, Integer *cilo, Integer *cihi, Integer *cjlo, Integer *cjhi)
{
    Integer atype, btype, ctype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    pnga_inquire_type(*g_c, &ctype);
    if (atype != btype || atype != ctype
            || (atype != C_INT && atype != C_LONG && atype != C_LONGLONG))
        pnga_error(" wrong types ", 0L);
    sga_add_patch(alpha, *g_a, ailo, aihi, ajlo, ajhi, beta, *g_b, bilo, bihi, bjlo, bjhi, *g_c, cilo, cihi, cjlo, cjhi);
}


void FATR ga_sadd_patch_(Real *alpha, Integer *g_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Real *beta, Integer *g_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, Integer *g_c, Integer *cilo, Integer *cihi, Integer *cjlo, Integer *cjhi)
{
    Integer atype, btype, ctype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    pnga_inquire_type(*g_c, &ctype);
    if (atype != btype || atype != ctype || atype != C_FLOAT)
        pnga_error(" wrong types ", 0L);
    sga_add_patch(alpha, *g_a, ailo, aihi, ajlo, ajhi, beta, *g_b, bilo, bihi, bjlo, bjhi, *g_c, cilo, cihi, cjlo, cjhi);
}


void FATR ga_zadd_patch_(DoubleComplex *alpha, Integer *g_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, DoubleComplex *beta, Integer *g_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, Integer *g_c, Integer *cilo, Integer *cihi, Integer *cjlo, Integer *cjhi)
{
    Integer atype, btype, ctype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    pnga_inquire_type(*g_c, &ctype);
    if (atype != btype || atype != ctype || atype != C_DCPL)
        pnga_error(" wrong types ", 0L);
    sga_add_patch(alpha, *g_a, ailo, aihi, ajlo, ajhi, beta, *g_b, bilo, bihi, bjlo, bjhi, *g_c, cilo, cihi, cjlo, cjhi);
}

void FATR nga_add_patch_(void *alpha, Integer *g_a, Integer *alo, Integer *ahi, void *beta, Integer *g_b, Integer *blo, Integer *bhi, Integer *g_c, Integer *clo, Integer *chi)
{
    wnga_add_patch(alpha, *g_a, alo, ahi, beta, *g_b, blo, bhi, *g_c, clo, chi);
}

void FATR nga_cadd_patch_(SingleComplex *alpha, Integer *g_a, Integer *alo, Integer *ahi, SingleComplex *beta, Integer *g_b, Integer *blo, Integer *bhi, Integer *g_c, Integer *clo, Integer *chi)
{
    Integer atype, btype, ctype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    pnga_inquire_type(*g_c, &ctype);
    if (atype != btype || atype != ctype || atype != C_SCPL)
        pnga_error(" wrong types ", 0L);
    wnga_add_patch(alpha, *g_a, alo, ahi, beta, *g_b, blo, bhi, *g_c, clo, chi);
}

void FATR nga_dadd_patch_(DoublePrecision *alpha, Integer *g_a, Integer *alo, Integer *ahi, DoublePrecision *beta, Integer *g_b, Integer *blo, Integer *bhi, Integer *g_c, Integer *clo, Integer *chi)
{
    Integer atype, btype, ctype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    pnga_inquire_type(*g_c, &ctype);
    if (atype != btype || atype != ctype || atype != C_DBL)
        pnga_error(" wrong types ", 0L);
    wnga_add_patch(alpha, *g_a, alo, ahi, beta, *g_b, blo, bhi, *g_c, clo, chi);
}

void FATR nga_iadd_patch_(Integer *alpha, Integer *g_a, Integer *alo, Integer *ahi, Integer *beta, Integer *g_b, Integer *blo, Integer *bhi, Integer *g_c, Integer *clo, Integer *chi)
{
    Integer atype, btype, ctype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    pnga_inquire_type(*g_c, &ctype);
    if (atype != btype || atype != ctype
            || (atype != C_INT && atype != C_LONG && atype != C_LONGLONG))
        pnga_error(" wrong types ", 0L);
    wnga_add_patch(alpha, *g_a, alo, ahi, beta, *g_b, blo, bhi, *g_c, clo, chi);
}

void FATR nga_sadd_patch_(Real *alpha, Integer *g_a, Integer *alo, Integer *ahi, Real *beta, Integer *g_b, Integer *blo, Integer *bhi, Integer *g_c, Integer *clo, Integer *chi)
{
    Integer atype, btype, ctype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    pnga_inquire_type(*g_c, &ctype);
    if (atype != btype || atype != ctype || atype != C_FLOAT)
        pnga_error(" wrong types ", 0L);
    wnga_add_patch(alpha, *g_a, alo, ahi, beta, *g_b, blo, bhi, *g_c, clo, chi);
}

void FATR nga_zadd_patch_(DoubleComplex *alpha, Integer *g_a, Integer *alo, Integer *ahi, DoubleComplex *beta, Integer *g_b, Integer *blo, Integer *bhi, Integer *g_c, Integer *clo, Integer *chi)
{
    Integer atype, btype, ctype;

    pnga_inquire_type(*g_a, &atype);
    pnga_inquire_type(*g_b, &btype);
    pnga_inquire_type(*g_c, &ctype);
    if (atype != btype || atype != ctype || atype != C_DCPL)
        pnga_error(" wrong types ", 0L);
    wnga_add_patch(alpha, *g_a, alo, ahi, beta, *g_b, blo, bhi, *g_c, clo, chi);
}

/* Routines from select.c */

void FATR nga_select_elem_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
    Integer *g_a, char* op, void* val, Integer *subscript, int oplen
#else
    Integer *g_a, char* op, int oplen, void* val, Integer *subscript
#endif
    )
{
    wnga_select_elem(*g_a, op, val, subscript);
}

/* Routines from sparse.c */

void FATR ga_patch_enum_(Integer* g_a, Integer* lo, Integer* hi, void* start, void* stride)
{
    wnga_patch_enum(*g_a, *lo, *hi, start, stride);
}

void FATR nga_patch_enum_(Integer* g_a, Integer* lo, Integer* hi, void* start, void* stride)
{
    wnga_patch_enum(*g_a, *lo, *hi, start, stride);
}

void FATR ga_scan_copy_(Integer* g_a, Integer* g_b, Integer* g_sbit, Integer* lo, Integer* hi)
{
    wnga_scan_copy(*g_a, *g_b, *g_sbit, *lo, *hi);
}

void FATR nga_scan_copy_(Integer* g_a, Integer* g_b, Integer* g_sbit, Integer* lo, Integer* hi)
{
    wnga_scan_copy(*g_a, *g_b, *g_sbit, *lo, *hi);
}

void FATR ga_scan_add_(Integer* g_a, Integer* g_b, Integer* g_sbit, Integer* lo, Integer* hi, Integer* excl)
{
    wnga_scan_add(*g_a, *g_b, *g_sbit, *lo, *hi, *excl);
}

void FATR nga_scan_add_(Integer* g_a, Integer* g_b, Integer* g_sbit, Integer* lo, Integer* hi, Integer* excl)
{
    wnga_scan_add(*g_a, *g_b, *g_sbit, *lo, *hi, *excl);
}

void FATR ga_pack_(Integer* g_a, Integer* g_b, Integer* g_sbit, Integer* lo, Integer* hi, Integer* icount)
{
    wnga_pack(*g_a, *g_b, *g_sbit, *lo, *hi, icount);
}

void FATR nga_pack_(Integer* g_a, Integer* g_b, Integer* g_sbit, Integer* lo, Integer* hi, Integer* icount)
{
    wnga_pack(*g_a, *g_b, *g_sbit, *lo, *hi, icount);
}

void FATR ga_unpack_(Integer* g_a, Integer* g_b, Integer* g_sbit, Integer* lo, Integer* hi, Integer* icount)
{
    wnga_unpack(*g_a, *g_b, *g_sbit, *lo, *hi, icount);
}

void FATR nga_unpack_(Integer* g_a, Integer* g_b, Integer* g_sbit, Integer* lo, Integer* hi, Integer* icount)
{
    wnga_unpack(*g_a, *g_b, *g_sbit, *lo, *hi, icount);
}

logical FATR ga_create_bin_range_(Integer *g_bin, Integer *g_cnt, Integer *g_off, Integer *g_range)
{
    return wnga_create_bin_range(*g_bin, *g_cnt, *g_off, g_range);
}

logical FATR nga_create_bin_range_(Integer *g_bin, Integer *g_cnt, Integer *g_off, Integer *g_range)
{
    return wnga_create_bin_range(*g_bin, *g_cnt, *g_off, g_range);
}

void FATR ga_bin_sorter_(Integer *g_bin, Integer *g_cnt, Integer *g_off)
{
    wnga_bin_sorter(*g_bin, *g_cnt, *g_off);
}

void FATR nga_bin_sorter_(Integer *g_bin, Integer *g_cnt, Integer *g_off)
{
    wnga_bin_sorter(*g_bin, *g_cnt, *g_off);
}

void FATR ga_bin_index_(Integer *g_bin, Integer *g_cnt, Integer *g_off, Integer *values, Integer *subs, Integer *n, Integer *sortit)
{
    wnga_bin_index(*g_bin, *g_cnt, *g_off, values, subs, *n, *sortit);
}

void FATR nga_bin_index_(Integer *g_bin, Integer *g_cnt, Integer *g_off, Integer *values, Integer *subs, Integer *n, Integer *sortit)
{
    wnga_bin_index(*g_bin, *g_cnt, *g_off, values, subs, *n, *sortit);
}

/* Routines from matrix.c */

void FATR ga_median_patch_(Integer *g_a, Integer *alo, Integer *ahi, Integer *g_b, Integer *blo, Integer *bhi, Integer *g_c, Integer *clo, Integer *chi, Integer *g_m, Integer *mlo, Integer *mhi)
{
    wnga_median_patch(*g_a, alo, ahi, *g_b, blo, bhi, *g_c, clo, chi, *g_m, mlo, mhi);
}

void FATR nga_median_patch_(Integer *g_a, Integer *alo, Integer *ahi, Integer *g_b, Integer *blo, Integer *bhi, Integer *g_c, Integer *clo, Integer *chi, Integer *g_m, Integer *mlo, Integer *mhi)
{
    wnga_median_patch(*g_a, alo, ahi, *g_b, blo, bhi, *g_c, clo, chi, *g_m, mlo, mhi);
}

void FATR ga_median_(Integer * g_a, Integer * g_b, Integer * g_c, Integer * g_m){
    wnga_median(*g_a, *g_b, *g_c, *g_m);
}

void FATR nga_median_(Integer * g_a, Integer * g_b, Integer * g_c, Integer * g_m){
    wnga_median(*g_a, *g_b, *g_c, *g_m);
}

void FATR ga_norm_infinity_(Integer * g_a, double *nm)
{
    wnga_norm_infinity(*g_a, nm);
}

void FATR nga_norm_infinity_(Integer * g_a, double *nm)
{
    wnga_norm_infinity(*g_a, nm);
}

void FATR ga_norm1_(Integer * g_a, double *nm)
{
    wnga_norm1(*g_a, nm);
}

void FATR nga_norm1_(Integer * g_a, double *nm)
{
    wnga_norm1(*g_a, nm);
}

void FATR ga_get_diag_(Integer * g_a, Integer * g_v)
{
    wnga_get_diag(*g_a, *g_v);
}

void FATR nga_get_diag_(Integer * g_a, Integer * g_v)
{
    wnga_get_diag(*g_a, *g_v);
}

void FATR ga_add_diagonal_(Integer * g_a, Integer * g_v)
{
    wnga_add_diagonal(*g_a, *g_v);
}

void FATR nga_add_diagonal_(Integer * g_a, Integer * g_v)
{
    wnga_add_diagonal(*g_a, *g_v);
}

void FATR ga_set_diagonal_(Integer * g_a, Integer * g_v)
{
    wnga_set_diagonal(*g_a, *g_v);
}

void FATR nga_set_diagonal_(Integer * g_a, Integer * g_v)
{
    wnga_set_diagonal(*g_a, *g_v);
}

void FATR ga_shift_diagonal_(Integer * g_a, void *c)
{
    wnga_shift_diagonal(*g_a, c);
}

void FATR nga_shift_diagonal_(Integer * g_a, void *c)
{
    wnga_shift_diagonal(*g_a, c);
}

void FATR ga_zero_diagonal_(Integer * g_a)
{
    wnga_zero_diagonal(*g_a);
}

void FATR nga_zero_diagonal_(Integer * g_a)
{
    wnga_zero_diagonal(*g_a);
}

void FATR ga_scale_rows_(Integer *g_a, Integer *g_v)
{
    wnga_scale_rows(*g_a, *g_v);
}

void FATR nga_scale_rows_(Integer *g_a, Integer *g_v)
{
    wnga_scale_rows(*g_a, *g_v);
}

void FATR ga_scale_cols_(Integer *g_a, Integer *g_v)
{
    wnga_scale_cols(*g_a, *g_v);
}

void FATR nga_scale_cols_(Integer *g_a, Integer *g_v)
{
    wnga_scale_cols(*g_a, *g_v);
}

/* Routines from ga_symmetr.c */

void FATR ga_symmetrize_(Integer *g_a)
{
    wnga_symmetrize(*g_a);
}

void FATR nga_symmetrize_(Integer *g_a)
{
    wnga_symmetrize(*g_a);
}

/* Routines from global.periodic.c */

void FATR nga_periodic_get_(Integer *g_a, Integer *lo, Integer *hi,
                            void *buf, Integer *ld)
{
    wnga_periodic(*g_a, lo, hi, buf, ld, NULL, PERIODIC_GET);
}

void FATR nga_periodic_put_(Integer *g_a, Integer *lo, Integer *hi,
                            void *buf, Integer *ld)
{
    wnga_periodic(*g_a, lo, hi, buf, ld, NULL, PERIODIC_PUT);
}

void FATR nga_periodic_acc_(Integer *g_a, Integer *lo, Integer *hi,
                            void *buf, Integer *ld, void *alpha)
{
    wnga_periodic(*g_a, lo, hi, buf, ld, alpha, PERIODIC_ACC);
}

/* Routines from matmul.c */

void FATR ga_matmul_mirrored_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        transa, transb, alpha, beta, g_a, ailo, aihi, ajlo, ajhi, g_b, bilo, bihi, bjlo, bjhi, g_c, cilo, cihi, cjlo, cjhi, alen, blen
#else
        transa, alen, transb, blen, alpha, beta, g_a, ailo, aihi, ajlo, ajhi, g_b, bilo, bihi, bjlo, bjhi, g_c, cilo, cihi, cjlo, cjhi
#endif
        )
Integer *g_a, *ailo, *aihi, *ajlo, *ajhi;    /* patch of g_a */
Integer *g_b, *bilo, *bihi, *bjlo, *bjhi;    /* patch of g_b */
Integer *g_c, *cilo, *cihi, *cjlo, *cjhi;    /* patch of g_c */
void    *alpha, *beta;
char    *transa, *transb;
int      alen, blen;
{
    wnga_matmul_mirrored(transa, transb, alpha, beta, *g_a, *ailo, *aihi, *ajlo, *ajhi, *g_b, *bilo, *bihi, *bjlo, *bjhi, *g_c, *cilo, *cihi, *cjlo, *cjhi);
}

void FATR nga_matmul_patch_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *transa, char *transb, void *alpha, void *beta, Integer *g_a, Integer alo[], Integer ahi[], Integer *g_b, Integer blo[], Integer bhi[], Integer *g_c, Integer clo[], Integer chi[], int alen, int blen
#else
        char *transa, int alen, char *transb, int blen, void *alpha, void *beta, Integer *g_a, Integer alo[], Integer ahi[], Integer *g_b, Integer blo[], Integer bhi[], Integer *g_c, Integer clo[], Integer chi[]
#endif
        )
{
    wnga_matmul_patch(transa, transb, alpha, beta, *g_a, alo, ahi, *g_b, blo, bhi, *g_c, clo, chi);
}

void FATR ga_matmul_patch_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *transa, char *transb, DoublePrecision *alpha, DoublePrecision *beta, Integer *g_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, Integer *g_c, Integer *cilo, Integer *cihi, Integer *cjlo, Integer *cjhi, int alen, int blen
#else
        char *transa, int alen, char *transb, int blen, DoublePrecision *alpha, DoublePrecision *beta, Integer *g_a, Integer *ailo, Integer *aihi, Integer *ajlo, Integer *ajhi, Integer *g_b, Integer *bilo, Integer *bihi, Integer *bjlo, Integer *bjhi, Integer *g_c, Integer *cilo, Integer *cihi, Integer *cjlo, Integer *cjhi
#endif
        )
{
#if 0
Integer alo[2], ahi[2]; 
Integer blo[2], bhi[2];
Integer clo[2], chi[2];
        alo[0]=*ailo; ahi[0]=*aihi; alo[1]=*ajlo; ahi[1]=*ajhi;
        blo[0]=*bilo; bhi[0]=*bihi; blo[1]=*bjlo; bhi[1]=*bjhi;
        clo[0]=*cilo; chi[0]=*cihi; clo[1]=*cjlo; chi[1]=*cjhi;
    pnga_matmul_patch(transa, transb, alpha, beta, g_a, alo, ahi,
                         g_b, blo, bhi, g_c, clo, chi);
#else
    if(pnga_is_mirrored(*g_a))
       wnga_matmul_mirrored(transa, transb, (void*)alpha, (void*)beta,
                  *g_a, *ailo, *aihi, *ajlo, *ajhi,
                  *g_b, *bilo, *bihi, *bjlo, *bjhi,
                  *g_c, *cilo, *cihi, *cjlo, *cjhi);
    else {
       gai_matmul_patch_flag(SET);
       wnga_matmul(transa, transb, (void*)alpha, (void*)beta,
             *g_a, *ailo, *aihi, *ajlo, *ajhi,
             *g_b, *bilo, *bihi, *bjlo, *bjhi,
             *g_c, *cilo, *cihi, *cjlo, *cjhi);
       gai_matmul_patch_flag(UNSET);
    }
#endif
}

void FATR ga_matmul_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        transa, transb, alpha, beta, g_a, ailo, aihi, ajlo, ajhi, g_b, bilo, bihi, bjlo, bjhi, g_c, cilo, cihi, cjlo, cjhi, alen, blen
#else
        transa, alen, transb, blen, alpha, beta, g_a, ailo, aihi, ajlo, ajhi, g_b, bilo, bihi, bjlo, bjhi, g_c, cilo, cihi, cjlo, cjhi
#endif
        )
Integer *g_a, *ailo, *aihi, *ajlo, *ajhi;    /* patch of g_a */
Integer *g_b, *bilo, *bihi, *bjlo, *bjhi;    /* patch of g_b */
Integer *g_c, *cilo, *cihi, *cjlo, *cjhi;    /* patch of g_c */
void    *alpha, *beta;
char    *transa, *transb;
int      alen, blen;
{
    wnga_matmul(transa, transb, alpha, beta, *g_a, *ailo, *aihi, *ajlo, *ajhi, *g_b, *bilo, *bihi, *bjlo, *bjhi, *g_c, *cilo, *cihi, *cjlo, *cjhi);
}

/* use ga_dgemm in ga_dgemmf.F as accumulate is sloooow in CRAY_XT */
#ifdef CRAY_XT
#   define GA_DGEMM ga_dgemm_DISABLE 
#else
#   define GA_DGEMM ga_dgemm_
#endif

#define  SET_GEMM_INDICES\
      Integer ailo = 1;\
  Integer aihi = *m;\
  Integer ajlo = 1;\
  Integer ajhi = *k;\
\
  Integer bilo = 1;\
  Integer bihi = *k;\
  Integer bjlo = 1;\
  Integer bjhi = *n;\
\
  Integer cilo = 1;\
  Integer cihi = *m;\
  Integer cjlo = 1;\
  Integer cjhi = *n

void FATR GA_DGEMM(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *transa, char *transb,
        Integer *m, Integer *n, Integer *k,
        void *alpha, Integer *g_a, Integer *g_b,
        void *beta, Integer *g_c, int talen, int tblen
#else
        char *transa, int talen, char *transb, int tblen,
        Integer *m, Integer *n, Integer *k,
        void *alpha, Integer *g_a, Integer *g_b,
        void *beta, Integer *g_c
#endif
        )
{
SET_GEMM_INDICES;

 wnga_matmul(transa, transb, alpha, beta,
       *g_a, ailo, aihi, ajlo, ajhi,
       *g_b, bilo, bihi, bjlo, bjhi,
       *g_c, cilo, cihi, cjlo, cjhi);
}

void FATR ga_cgemm_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *transa, char *transb,
        Integer *m, Integer *n, Integer *k,
        void *alpha, Integer *g_a, Integer *g_b,
        void *beta, Integer *g_c, int talen, int tblen
#else
        char *transa, int talen, char *transb, int tblen,
        Integer *m, Integer *n, Integer *k,
        void *alpha, Integer *g_a, Integer *g_b,
        void *beta, Integer *g_c
#endif
        )
{
SET_GEMM_INDICES;

  wnga_matmul (transa, transb, alpha, beta,
         *g_a, ailo, aihi, ajlo, ajhi,
         *g_b, bilo, bihi, bjlo, bjhi,
         *g_c, cilo, cihi, cjlo, cjhi);
}

void FATR ga_sgemm_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *transa, char *transb,
        Integer *m, Integer *n, Integer *k,
        void *alpha, Integer *g_a, Integer *g_b,
        void *beta, Integer *g_c, int talen, int tblen
#else
        char *transa, int talen, char *transb, int tblen,
        Integer *m, Integer *n, Integer *k,
        void *alpha, Integer *g_a, Integer *g_b,
        void *beta, Integer *g_c
#endif
        )
{
SET_GEMM_INDICES;

  wnga_matmul (transa, transb, alpha, beta,
         *g_a, ailo, aihi, ajlo, ajhi,
         *g_b, bilo, bihi, bjlo, bjhi,
         *g_c, cilo, cihi, cjlo, cjhi);
}

void FATR ga_zgemm_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *transa, char *transb,
        Integer *m, Integer *n, Integer *k,
        void *alpha, Integer *g_a, Integer *g_b,
        void *beta, Integer *g_c, int talen, int tblen
#else
        char *transa, int talen, char *transb, int tblen,
        Integer *m, Integer *n, Integer *k,
        void *alpha, Integer *g_a, Integer *g_b,
        void *beta, Integer *g_c
#endif
        )
{
SET_GEMM_INDICES;

  wnga_matmul (transa, transb, alpha, beta,
         *g_a, ailo, aihi, ajlo, ajhi,
         *g_b, bilo, bihi, bjlo, bjhi,
         *g_c, cilo, cihi, cjlo, cjhi);
}

/* Routines from ga_diag_seqc.c */

void FATR ga_diag_seq_(Integer *g_a, Integer *g_s, Integer *g_v, DoublePrecision *eval)
{
    wnga_diag_seq(*g_a, *g_s, *g_v, eval);
}

void FATR ga_diag_std_seq_(Integer * g_a, Integer * g_v, DoublePrecision *eval)
{
    wnga_diag_std_seq(*g_a, *g_v, eval);
}

/* Routines from peigstubs.c */

void FATR ga_diag_(Integer * g_a, Integer * g_s, Integer * g_v, DoublePrecision *eval)
{
    wnga_diag(*g_a, *g_s, *g_v, eval);
}

void FATR ga_diag_std_(Integer * g_a, Integer * g_v, DoublePrecision *eval)
{
    wnga_diag_std(*g_a, *g_v, eval);
}

void FATR ga_diag_reuse_(Integer * reuse, Integer * g_a, Integer * g_s,
           Integer * g_v, DoublePrecision *eval)
{
    wnga_diag_reuse(*reuse, *g_a, *g_s, *g_v, eval);
}

/* Routines from sclstubs.c */

void FATR ga_lu_solve_alt_(Integer *tran, Integer * g_a, Integer * g_b)
{
    wnga_lu_solve_alt(*tran, *g_a, *g_b);
}

void FATR ga_lu_solve_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *tran, Integer * g_a, Integer * g_b, int len
#else
        char *tran, int len, Integer * g_a, Integer * g_b
#endif
        )
{
    wnga_lu_solve(tran, *g_a, *g_b);
}

Integer FATR ga_llt_solve_(Integer * g_a, Integer * g_b)
{
    return wnga_llt_solve(*g_a, *g_b);
}

Integer FATR nga_llt_solve_(Integer * g_a, Integer * g_b)
{
    return wnga_llt_solve(*g_a, *g_b);
}

Integer FATR ga_solve_(Integer * g_a, Integer * g_b)
{
    return wnga_solve(*g_a, *g_b);
}

Integer FATR nga_solve_(Integer * g_a, Integer * g_b)
{
    return wnga_solve(*g_a, *g_b);
}

Integer FATR ga_spd_invert_(Integer * g_a)
{
    return wnga_spd_invert(*g_a);
}

Integer FATR nga_spd_invert_(Integer * g_a)
{
    return wnga_spd_invert(*g_a);
}

/* Routines from DP.c */

void FATR ga_copy_patch_dp_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        trans, g_a, ailo, aihi, ajlo, ajhi, g_b, bilo, bihi, bjlo, bjhi, translen
#else
        trans, translen, g_a, ailo, aihi, ajlo, ajhi, g_b, bilo, bihi, bjlo, bjhi
#endif
        )
Integer *g_a, *ailo, *aihi, *ajlo, *ajhi;
Integer *g_b, *bilo, *bihi, *bjlo, *bjhi;
char *trans;
int translen;
{
    wnga_copy_patch_dp(trans,*g_a,*ailo,*aihi,*ajlo,*ajhi,*g_b,*bilo,*bihi,*bjlo,*bjhi);
}

DoublePrecision FATR ga_ddot_patch_dp_(
#if F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        g_a, t_a, ailo, aihi, ajlo, ajhi, g_b, t_b, bilo, bihi, bjlo, bjhi, alen, blen
#else
        g_a, t_a, alen, ailo, aihi, ajlo, ajhi, g_b, t_b, blen, bilo, bihi, bjlo, bjhi
#endif
        )
Integer *g_a, *ailo, *aihi, *ajlo, *ajhi;    /* patch of g_a */
Integer *g_b, *bilo, *bihi, *bjlo, *bjhi;    /* patch of g_b */
char    *t_a, *t_b;                          /* transpose operators */
int alen, blen;
{
    return wnga_ddot_patch_dp(*g_a, t_a, *ailo, *aihi, *ajlo, *ajhi, *g_b, t_b, *bilo, *bihi, *bjlo, *bjhi);
}

/* Routines from ga_trace.c */

double FATR ga_timer_()
{
    return wnga_timer();
}

double FATR nga_timer_()
{
    return wnga_timer();
}

Integer nga_register_type_(Integer *size) {
  return wnga_register_type((size_t)*size);
}

Integer nga_deregister_type_(Integer *type) {
  return wnga_deregister_type((int)*type);
}

void nga_get_field_(Integer *g_a, Integer *lo, Integer *hi, Integer *foff, Integer *fsize,
		   void *buf, Integer *ld) {
  wnga_get_field(*g_a, lo, hi, *foff, *fsize, buf, ld);
}

void nga_nbget_field_(Integer *g_a, Integer *lo, Integer *hi, Integer *foff, Integer *fsize,
		     void *buf, Integer *ld, Integer *nbhandle) {
  wnga_nbget_field(*g_a, lo, hi, *foff, *fsize, buf, ld, nbhandle);
}

void nga_nbput_field_(Integer *g_a, Integer *lo, Integer *hi, Integer *foff, Integer *fsize,
		     void *buf, Integer *ld, Integer *nbhandle) {
  wnga_nbput_field(*g_a, lo, hi, *foff, *fsize, buf, ld, nbhandle);
}

void nga_put_field_(Integer *g_a, Integer *lo, Integer *hi, Integer *foff, Integer *fsize,
		   void *buf, Integer *ld) {
  wnga_put_field(*g_a, lo, hi, *foff, *fsize, buf, ld);
}
