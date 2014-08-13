/* $Id: onesided.c,v 1.80.2.18 2007/12/18 22:22:27 d3g293 Exp $ */
/* 
 * module: onesided.c
 * author: Jarek Nieplocha
 * description: implements GA primitive communication operations --
 *              accumulate, scatter, gather, read&increment & synchronization 
 * 
 * DISCLAIMER
 *
 * This material was prepared as an account of work sponsored by an
 * agency of the United States Government.  Neither the United States
 * Government nor the United States Department of Energy, nor Battelle,
 * nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
 * COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
 * SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
 * INFRINGE PRIVATELY OWNED RIGHTS.
 *
 *
 * ACKNOWLEDGMENT
 *
 * This software and its documentation were produced with United States
 * Government support under Contract Number DE-AC06-76RLO-1830 awarded by
 * the United States Department of Energy.  The United States Government
 * retains a paid-up non-exclusive, irrevocable worldwide license to
 * reproduce, prepare derivative works, perform publicly and display
 * publicly by or for the US Government, including the right to
 * distribute to other US Government contractors.
 */
#if HAVE_CONFIG_H
#   include "config.h"
#endif
 
/*#define PERMUTE_PIDS */
/*#define USE_GATSCAT_NEW 1 */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STDINT_H
#   include <stdint.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif
#if HAVE_STDDEF_H
#include <stddef.h>
#endif

#include "global.h"
#include "globalp.h"
#include "base.h"
#include "armci.h"
#include "macdecls.h"
#include "ga-papi.h"
#include "ga-wapi.h"

#define DEBUG 0
#define USE_MALLOC 1
#define INVALID_MA_HANDLE -1 
#define NEAR_INT(x) (x)< 0.0 ? ceil( (x) - 0.5) : floor((x) + 0.5)

#if !defined(CRAY_YMP)
#define BYTE_ADDRESSABLE_MEMORY
#endif

#ifdef PROFILE_OLD
#include "ga_profile.h"
#endif


#define DISABLE_NBOPT /* disables Non-Blocking OPTimization */

/*uncomment line below to verify consistency of MA in every sync */
/*#define CHECK_MA yes */

char *fence_array;
static int GA_fence_set=0;
Integer *_ga_map;       /* used in get/put/acc */
int *ProcListPerm;

static int GA_prealloc_gatscat = 0;
static Integer *GA_header;
static Integer *GA_list;
static Integer *GA_elems;

extern void armci_read_strided(void*, int, int*, int*, char*);
extern void armci_write_strided(void*, int, int*, int*, char*);
extern armci_hdl_t* get_armci_nbhandle(Integer *);

/***************************************************************************/

/**
 *  Synchronize all processes on group
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_pgroup_sync = pnga_pgroup_sync
#endif

void pnga_pgroup_sync(Integer grp_id)
{
#ifdef CHECK_MA
    Integer status;
#endif
 
/*    printf("p[%d] calling ga_pgroup_sync on group: %d\n",GAme,*grp_id); */
    if (grp_id > 0) {
#   ifdef MPI
       /* fence on all processes in this group */
       { int p; for(p=0;p<PGRP_LIST[grp_id].map_nproc;p++)
          ARMCI_Fence(ARMCI_Absolute_id(&PGRP_LIST[grp_id].group, p)); }
#ifdef BGML 
       ARMCI_Barrier();
#endif
       pnga_msg_pgroup_sync(grp_id);
       if(GA_fence_set)bzero(fence_array,(int)GAnproc);
       GA_fence_set=0;
#   else
       pnga_error("ga_pgroup_sync_(): MPI not defined. ga_msg_pgroup_sync_()  can be called only if GA is built with MPI", 0);
#   endif
    } else {
      /* printf("p[%d] calling regular sync in ga_pgroup_sync\n",GAme); */
       ARMCI_AllFence();
#ifdef BGML 
       ARMCI_Barrier();
#endif
       pnga_msg_pgroup_sync(grp_id);
       if(GA_fence_set)bzero(fence_array,(int)GAnproc);
       GA_fence_set=0;
    }
#ifdef CHECK_MA
    status = MA_verify_allocator_stuff();
#endif
}

/**
 *  Synchronize all processes and complete all outstanding communications
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_sync = pnga_sync
#endif

void pnga_sync()
{
#ifdef CHECK_MA
Integer status;
#endif
       
       if (GA_Default_Proc_Group == -1) {
	  ARMCI_AllFence();
#ifdef BGML
          ARMCI_Barrier();
#endif
	  pnga_msg_sync();
	  if(GA_fence_set)bzero(fence_array,(int)GAnproc);
	  GA_fence_set=0;
       } else {
         Integer grp_id = (Integer)GA_Default_Proc_Group;
#ifdef BGML
          ARMCI_Barrier();
#endif
         pnga_pgroup_sync(grp_id);
       }
#ifdef CHECK_MA
       status = MA_verify_allocator_stuff();
#endif
}


/**
 *  Wait until all requests initiated by calling process are completed
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_fence = pnga_fence
#endif

void pnga_fence()
{
    int proc;
    if(GA_fence_set<1)pnga_error("ga_fence: fence not initialized",0);
    GA_fence_set--;
    for(proc=0;proc<GAnproc;proc++)if(fence_array[proc])ARMCI_Fence(proc);
    bzero(fence_array,(int)GAnproc);
}

/**
 *  Initialize tracing of request completion
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_init_fence = pnga_init_fence
#endif

void pnga_init_fence()
{
    GA_fence_set++;
}

void gai_init_onesided()
{
    _ga_map = (Integer*)malloc((size_t)(GAnproc*2*MAXDIM +1)*sizeof(Integer));
    if(!_ga_map) pnga_error("ga_init:malloc failed (_ga_map)",0);
    fence_array = calloc((size_t)GAnproc,1);
    if(!fence_array) pnga_error("ga_init:calloc failed",0);
}


void gai_finalize_onesided()
{
    free(_ga_map);
    _ga_map = NULL;
    free(fence_array);
    fence_array = NULL;
}


/*\ prepare permuted list of processes for remote ops
\*/
#define gaPermuteProcList(nproc)\
{\
  if((nproc) ==1) ProcListPerm[0]=0; \
  else{\
    int _i, iswap, temp;\
    if((nproc) > GAnproc) pnga_error("permute_proc: error ", (nproc));\
\
    /* every process generates different random sequence */\
    (void)srand((unsigned)GAme); \
\
    /* initialize list */\
    for(_i=0; _i< (nproc); _i++) ProcListPerm[_i]=_i;\
\
    /* list permutation generated by random swapping */\
    for(_i=0; _i< (nproc); _i++){ \
      iswap = (int)(rand() % (nproc));  \
      temp = ProcListPerm[iswap]; \
      ProcListPerm[iswap] = ProcListPerm[_i]; \
      ProcListPerm[_i] = temp; \
    } \
  }\
}
     




/*\ internal malloc that bypasses MA and uses internal buf when possible
\*/
#define MBUFLEN 256
#define MBUF_LEN MBUFLEN+2
static double ga_int_malloc_buf[MBUF_LEN];
static int mbuf_used=0;
#define MBUF_GUARD -1998.1998
void *gai_malloc(int bytes)
{
    void *ptr;
    if(!mbuf_used && bytes <= MBUF_LEN){
       if(DEBUG){
          ga_int_malloc_buf[0]= MBUF_GUARD;
          ga_int_malloc_buf[MBUFLEN]= MBUF_GUARD;
       }
       ptr = ga_int_malloc_buf+1;
       mbuf_used++;
    }else{
        Integer elems = (bytes+sizeof(double)-1)/sizeof(double)+1;
	ptr=ga_malloc(elems, MT_DBL, "GA malloc temp");
	/* *((Integer*)ptr)= handle;
	   ptr = ((double*)ptr)+ 1;*/ /*needs sizeof(double)>=sizeof(Integer)*/
    }
    return ptr;
}

void gai_free(void *ptr)
{
    if(ptr == (ga_int_malloc_buf+1)){
        if(DEBUG){
          assert(ga_int_malloc_buf[0]== MBUF_GUARD);
          assert(ga_int_malloc_buf[MBUFLEN]== MBUF_GUARD);
          assert(mbuf_used ==1);
        }
        mbuf_used =0;
    }else{
        /* Integer handle= *( (Integer*) (-1 + (double*)ptr)); */
      ga_free(ptr);
    }
}

        


#define gaShmemLocation(proc, g_a, _i, _j, ptr_loc, _pld)                      \
{                                                                              \
  Integer _ilo, _ihi, _jlo, _jhi, offset, proc_place, g_handle=(g_a)+GA_OFFSET;\
  Integer _lo[2], _hi[2], _p_handle;                                           \
  Integer _iw = GA[g_handle].width[0];                                         \
  Integer _jw = GA[g_handle].width[1];                                         \
  Integer _num_rstrctd = GA[g_handle].num_rstrctd;                             \
                                                                               \
  ga_ownsM(g_handle, (proc),_lo,_hi);                                          \
  _p_handle = GA[g_handle].p_handle;                                           \
  if (_p_handle != 0) {                                                        \
    proc_place =  proc;                                                        \
    if (_num_rstrctd != 0) {                                                   \
      proc_place = GA[g_handle].rstrctd_list[proc_place];                      \
    }                                                                          \
  } else {                                                                     \
    proc_place = PGRP_LIST[_p_handle].inv_map_proc_list[proc];                 \
  }                                                                            \
  _ilo = _lo[0]; _ihi=_hi[0];                                                  \
  _jlo = _lo[1]; _jhi=_hi[1];                                                  \
                                                                               \
  if((_i)<_ilo || (_i)>_ihi || (_j)<_jlo || (_j)>_jhi){                        \
    char err_string[ERR_STR_LEN];                                              \
    sprintf(err_string,"%s:p=%ld invalid i/j (%ld,%ld)><(%ld:%ld,%ld:%ld)",    \
        "gaShmemLocation", (long)proc, (long)(_i),(long)(_j),                  \
        (long)_ilo, (long)_ihi, (long)_jlo, (long)_jhi);                       \
    pnga_error(err_string, g_a );                                                \
  }                                                                            \
  offset = ((_i)-_ilo+_iw) + (_ihi-_ilo+1+2*_iw)*((_j)-_jlo+_jw);              \
                                                                               \
  /* find location of the proc in current cluster pointer array */             \
  *(ptr_loc) = GA[g_handle].ptr[proc_place] +                                  \
  offset*GAsizeofM(GA[g_handle].type);                                         \
  *(_pld) = _ihi-_ilo+1+2*_iw;                                                 \
}

/*\ Return pointer (ptr_loc) to location in memory of element with subscripts
 *  (subscript). Also return physical dimensions of array in memory in ld.
\*/
#define gam_Location(proc, g_handle,  subscript, ptr_loc, ld)                  \
{                                                                              \
Integer _offset=0, _d, _w, _factor=1, _last=GA[g_handle].ndim-1;               \
Integer _lo[MAXDIM], _hi[MAXDIM], _pinv, _p_handle;                            \
                                                                               \
      ga_ownsM(g_handle, proc, _lo, _hi);                                      \
      gaCheckSubscriptM(subscript, _lo, _hi, GA[g_handle].ndim);               \
      if(_last==0) ld[0]=_hi[0]- _lo[0]+1+2*(Integer)GA[g_handle].width[0];    \
      __CRAYX1_PRAGMA("_CRI shortloop");                                       \
      for(_d=0; _d < _last; _d++)            {                                 \
          _w = (Integer)GA[g_handle].width[_d];                                \
          _offset += (subscript[_d]-_lo[_d]+_w) * _factor;                     \
          ld[_d] = _hi[_d] - _lo[_d] + 1 + 2*_w;                               \
          _factor *= ld[_d];                                                   \
      }                                                                        \
      _offset += (subscript[_last]-_lo[_last]                                  \
               + (Integer)GA[g_handle].width[_last])                           \
               * _factor;                                                      \
      _p_handle = GA[g_handle].p_handle;                                       \
      if (_p_handle != 0) {                                                    \
        if (GA[g_handle].num_rstrctd == 0) {                                   \
          _pinv = proc;                                                        \
        } else {                                                               \
          _pinv = GA[g_handle].rstrctd_list[proc];                             \
        }                                                                      \
      } else {                                                                 \
        _pinv = PGRP_LIST[_p_handle].inv_map_proc_list[proc];                  \
      }                                                                        \
      *(ptr_loc) =  GA[g_handle].ptr[_pinv]+_offset*GA[g_handle].elemsize;     \
}


#define gam_GetRangeFromMap(p, ndim, plo, phi){\
Integer   _mloc = p* ndim *2;\
          *plo  = (Integer*)_ga_map + _mloc;\
          *phi  = *plo + ndim;\
}

/* compute index of point subscripted by plo relative to point
   subscripted by lo, for a block with dimensions dims */
#define gam_ComputePatchIndex(ndim, lo, plo, dims, pidx){                      \
Integer _d, _factor;                                                           \
          *pidx = plo[0] -lo[0];                                               \
          __CRAYX1_PRAGMA("_CRI shortloop");                                   \
          for(_d= 0,_factor=1; _d< ndim -1; _d++){                             \
             _factor *= (dims[_d]);                                            \
             *pidx += _factor * (plo[_d+1]-lo[_d+1]);                          \
          }                                                                    \
}

#define gam_GetBlockPatch(plo,phi,lo,hi,blo,bhi,ndim) {                        \
  Integer _d;                                                                  \
  for (_d=0; _d<ndim; _d++) {                                                  \
    if (lo[_d] <= phi[_d] && lo[_d] >= plo[_d]) blo[_d] = lo[_d];              \
    else blo[_d] = plo[_d];                                                    \
    if (hi[_d] <= phi[_d] && hi[_d] >= plo[_d]) bhi[_d] = hi[_d];              \
    else bhi[_d] = phi[_d];                                                    \
  }                                                                            \
}

/**
 *  A routine to test for a non-blocking call
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_nbtest = pnga_nbtest
#endif

Integer pnga_nbtest(Integer *nbhandle) 
{
    return nga_test_internal((Integer *)nbhandle);
} 

/**
 *  A routine to wait for a non-blocking call to complete
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_nbwait = pnga_nbwait
#endif

void pnga_nbwait(Integer *nbhandle) 
{
    nga_wait_internal((Integer *)nbhandle);
} 

static void ngai_puts(char *loc_base_ptr, char *pbuf, int *stride_loc, char *prem, int *stride_rem,
		      int *count, int nstrides, int proc, int field_off, 
		      int field_size, int type_size) {
  if(field_size<0 || field_size == type_size) {
    ARMCI_PutS(pbuf,stride_loc,prem,stride_rem,count,nstrides,proc);
  }
  else {
    int i;
    count -= 1;
    stride_loc -= 1;
    stride_rem -= 1;
    nstrides += 1;

    for(i=1; i<nstrides; i++) {
      stride_loc[i] /= type_size;
      stride_loc[i] *= field_size;
    }

    pbuf = loc_base_ptr + (pbuf - loc_base_ptr)/type_size*field_size;
    prem += field_off;

    count[1] /= type_size; 
    ARMCI_PutS(pbuf,stride_loc,prem,stride_rem,count,nstrides,proc);
    count[1] *= type_size; /*restore*/
    for(i=1; i<nstrides; i++) {
      stride_loc[i] /= field_size;
      stride_loc[i] *= type_size;
    }
  }
}


static void ngai_nbputs(char *loc_base_ptr, char *pbuf,int *stride_loc, char *prem, int *stride_rem,
		      int *count, int nstrides, int proc, int field_off, 
			int field_size, int type_size, armci_hdl_t *nbhandle) {
  if(field_size<0 || field_size == type_size) {
    ARMCI_NbPutS(pbuf,stride_loc,prem,stride_rem,count,nstrides,proc,nbhandle);
  }
  else {
    int i;
    count -= 1;
    stride_loc -= 1;
    stride_rem -= 1;
    nstrides += 1;

    for(i=1; i<nstrides; i++) {
      stride_loc[i] /= type_size;
      stride_loc[i] *= field_size;
    }

    pbuf = loc_base_ptr + (pbuf - loc_base_ptr)/type_size*field_size;
    prem += field_off;

    count[1] /= type_size; 
    ARMCI_NbPutS(pbuf,stride_loc,prem,stride_rem,count,nstrides,proc, nbhandle);
    count[1] *= type_size; /*restore*/
    for(i=1; i<nstrides; i++) {
      stride_loc[i] /= field_size;
      stride_loc[i] *= type_size;
    }
  }
}

static void ngai_nbgets(char *loc_base_ptr, char *prem,int *stride_rem, char *pbuf, int *stride_loc,
			int *count, int nstrides, int proc, int field_off, 
			int field_size, int type_size, armci_hdl_t *nbhandle) {
  if(field_size<0 || field_size == type_size) {
    ARMCI_NbGetS(prem,stride_rem,pbuf,stride_loc,count,nstrides,proc,nbhandle);
  }
  else {
    int i;
    count -= 1;
    stride_loc -= 1;
    stride_rem -= 1;

    pbuf = loc_base_ptr + (pbuf - loc_base_ptr)/type_size*field_size;
    prem += field_off;

    count[1] /= type_size; 
    nstrides += 1;

    for(i=1; i<nstrides; i++) {
      stride_loc[i] /= type_size;
      stride_loc[i] *= field_size;
    }
/*     { */
/*       int i; */
/*       printf("%s: calling armci_nbgets: field_off=%d field_size=%d type_size=%d proc=%d nstrides=%d\n",  */
/* 	     __FUNCTION__, (int)field_off, (int)field_size, (int)type_size, proc, nstrides); */
/*       printf("me=%d loc_ptr=%p rem_ptr=%p count=[", pnga_nodeid(), pbuf, prem); */
/*       for(i=0; i<=nstrides; i++) { */
/* 	printf("%d ",count[i]); */
/*       } */
/*       printf("]\n"); */
/*       printf("src_stride_arr=["); */
/*       for(i=0; i<nstrides; i++) { */
/* 	printf("%d ",stride_rem[i]); */
/*       } */
/*       printf("]\n"); */
/*       printf("dst_stride_arr=["); */
/*       for(i=0; i<nstrides; i++) { */
/* 	printf("%d ",stride_loc[i]); */
/*       } */
/*       printf("]\n"); */
/*     } */
    ARMCI_NbGetS(prem,stride_rem,pbuf,stride_loc,count,nstrides,proc,nbhandle);
    count[1] *= type_size; /*restore*/
    for(i=1; i<nstrides; i++) {
      stride_loc[i] /= field_size;
      stride_loc[i] *= type_size;
    }
  }
}

static void ngai_gets(char *loc_base_ptr, char *prem,int *stride_rem, char *pbuf, int *stride_loc,
		      int *count, int nstrides, int proc, int field_off, 
		      int field_size, int type_size) {
#if 1
  armci_hdl_t nbhandle;
  ARMCI_INIT_HANDLE(&nbhandle);
  ngai_nbgets(loc_base_ptr, prem, stride_rem, pbuf, stride_loc, count, nstrides, proc, 
	      field_off, field_size, type_size, &nbhandle);
  ARMCI_Wait(&nbhandle);
#else
  if(field_size<0 || field_size == type_size) {
    ARMCI_GetS(prem,stride_rem,pbuf,stride_loc,count,nstrides,proc);
  } else {
    int i;
    count -= 1;
    stride_loc -= 1;
    stride_rem -= 1;

    pbuf = loc_base_ptr + (pbuf - loc_base_ptr)/type_size*field_size;
    prem += field_off;

    count[1] /= type_size; 
    nstrides += 1;

    for(i=1; i<nstrides; i++) {
      stride_loc[i] /= type_size;
      stride_loc[i] *= field_size;
    }
    ARMCI_GetS(prem,stride_rem,pbuf,stride_loc,count,nstrides,proc);
    count[1] *= type_size; /*restore*/
    for(i=1; i<nstrides; i++) {
      stride_loc[i] /= field_size;
      stride_loc[i] *= type_size;
    }
  }
#endif
}

/**
 *  A common routine called by both non-blocking and blocking GA put calls.
 */
#ifdef __crayx1
#pragma _CRI inline pnga_locate_region
#endif
void ngai_put_common(Integer g_a, 
                   Integer *lo,
                   Integer *hi,
                   void    *buf,
                   Integer *ld,
		     Integer field_off,
		     Integer field_size,
		     Integer *nbhandle) 
{
  Integer  p, np, handle=GA_OFFSET + g_a;
  Integer  idx, elems, size, p_handle;
  int proc, ndim, loop, cond;
  int num_loops=2; /* 1st loop for remote procs; 2nd loop for local procs */
  Integer use_blocks;
  Integer n_rstrctd;
  Integer *rank_rstrctd;
#if defined(__crayx1) || defined(DISABLE_NBOPT)
#else
  Integer ga_nbhandle;
  int counter=0;
#endif

  int _stride_rem[MAXDIM+1], _stride_loc[MAXDIM+1], _count[MAXDIM+1];
  int *stride_rem=&_stride_rem[1], *stride_loc=&_stride_loc[1], *count=&_count[1];

  GA_PUSH_NAME("ngai_put_common");

  ga_check_handleM(g_a, "ngai_put_common");

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  use_blocks = GA[handle].block_flag;
  p_handle = GA[handle].p_handle;
  n_rstrctd = (Integer)GA[handle].num_rstrctd;
  rank_rstrctd = GA[handle].rank_rstrctd;

  /*initial  stride portion for field-wise operations*/
  _stride_rem[0] = size;
  _stride_loc[0] = field_size;
  _count[0] = field_size;

  if (!use_blocks) {

    /* Locate the processors containing some portion of the patch
       specified by lo and hi and return the results in _ga_map,
       GA_proclist, and np. GA_proclist contains a list of processors
       containing some portion of the patch, _ga_map contains
       the lower and upper indices of the portion of the patch held
       by a given processor, and np contains the total number of
       processors that contain some portion of the patch.
     */
    if(!pnga_locate_region(g_a, lo, hi, _ga_map, GA_proclist, &np ))
      ga_RegionError(pnga_ndim(g_a), lo, hi, g_a);

#ifndef NO_GA_STATS
    gam_CountElems(ndim, lo, hi, &elems);
    GAbytes.puttot += (double)size*elems;
    GAstat.numput++;
    GAstat.numput_procs += np;
#endif

    if(nbhandle)ga_init_nbhandle(nbhandle);
#if !defined(__crayx1) && !defined(DISABLE_NBOPT)
    else ga_init_nbhandle(&ga_nbhandle);
#endif

#ifdef PROFILE_OLD
    ga_profile_start((int)handle, (long)size*elems, ndim, lo, hi, 
        ENABLE_PROFILE_PUT);
#endif

    gaPermuteProcList(np);

#if !defined(__crayx1) && !defined(DISABLE_NBOPT)
    for(loop=0; loop<num_loops; loop++) {
      __CRAYX1_PRAGMA("_CRI novector");
#endif
      for(idx=0; idx<np; idx++){
        Integer ldrem[MAXDIM];
        Integer idx_buf, *plo, *phi;
        char *pbuf, *prem;

        p = (Integer)ProcListPerm[idx];
        proc = (int)GA_proclist[p];
        if (p_handle >= 0) {
          proc = PGRP_LIST[p_handle].inv_map_proc_list[proc];
        }
#ifdef PERMUTE_PIDS
        if(GA_Proc_list) proc = GA_inv_Proc_list[proc];
#endif

        /* check if it is local to SMP */

#if !defined(__crayx1) && !defined(DISABLE_NBOPT)
        cond = armci_domain_same_id(ARMCI_DOMAIN_SMP,(int)proc);
        if(loop==0) cond = !cond;
        if(cond) {
#endif
          /* Find  visible portion of patch held by processor p and
             return the result in plo and phi. Also get actual processor
             index corresponding to p and store the result in proc. */
          gam_GetRangeFromMap(p, ndim, &plo, &phi);
          proc = (int)GA_proclist[p];

          if (n_rstrctd == 0) {
            gam_Location(proc,handle, plo, &prem, ldrem);
          } else {
            gam_Location(rank_rstrctd[proc],handle, plo, &prem, ldrem);
          }

          /* find the right spot in the user buffer */
          gam_ComputePatchIndex(ndim, lo, plo, ld, &idx_buf);
          pbuf = size*idx_buf + (char*)buf;        

          gam_ComputeCount(ndim, plo, phi, count); 

          /* scale number of rows by element size */
          count[0] *= size;
          gam_setstride(ndim, size, ld, ldrem, stride_rem, stride_loc);

          if (p_handle >= 0) {
            proc = (int)GA_proclist[p];
            /* BJP */
            proc = PGRP_LIST[p_handle].inv_map_proc_list[proc];
          }
          if(GA_fence_set)fence_array[proc]=1;

#ifdef PERMUTE_PIDS
          if(GA_Proc_list) proc = GA_inv_Proc_list[proc];
#endif

#ifndef NO_GA_STATS	    
          if(proc == GAme){
            gam_CountElems(ndim, plo, phi, &elems);
            GAbytes.putloc += (double)size*elems;
          }
#endif

          /*casting what ganb_get_armci_handle function returns to armci_hdl is 
            very crucial here as on 64 bit platforms, pointer is 64 bits where 
            as temporary is only 32 bits*/ 
#if defined(__crayx1) || defined(DISABLE_NBOPT)
/*           ARMCI_PutS(pbuf,stride_loc,prem,stride_rem,count,ndim-1,proc); */
	  ngai_puts(buf, pbuf,stride_loc,prem,stride_rem,count,ndim-1,proc, field_off, field_size, size);
#else
          if(nbhandle)  {
/*             ARMCI_NbPutS(pbuf, stride_loc, prem, stride_rem, count, ndim -1, */
/*                 proc,(armci_hdl_t*)get_armci_nbhandle(nbhandle)); */
            ngai_nbputs(buf,pbuf, stride_loc, prem, stride_rem, count, ndim -1,
			proc,field_off, field_size, size, 
			(armci_hdl_t*)get_armci_nbhandle(nbhandle));
	  }
          else {
            /* do blocking put for local processes. If all processes
               are remote processes then do blocking put for the last one */
            if((loop==0 && counter==(int)np-1) || loop==1) {
/*               ARMCI_PutS(pbuf,stride_loc,prem,stride_rem,count,ndim-1,proc); */
              ngai_puts(buf, pbuf,stride_loc,prem,stride_rem,count,ndim-1,proc, field_off, field_size, size);
	    }
            else {
              ++counter;
/*               ARMCI_NbPutS(pbuf,stride_loc,prem,stride_rem,count, ndim-1, */
/*                   proc,(armci_hdl_t*)get_armci_nbhandle(&ga_nbhandle)); */
              ngai_nbputs(buf,pbuf,stride_loc,prem,stride_rem,count, ndim-1,
			  proc, field_off, field_size, size,
			  (armci_hdl_t*)get_armci_nbhandle(&ga_nbhandle));
            }
          }
        } /* end if(cond) */
      }
#endif
    }
#if !defined(__crayx1) && !defined(DISABLE_NBOPT)
    if(!nbhandle) nga_wait_internal(&ga_nbhandle);  
#endif
  } else {
    Integer offset, l_offset, last, pinv;
    Integer blo[MAXDIM],bhi[MAXDIM];
    Integer plo[MAXDIM],phi[MAXDIM];
    Integer idx, j, jtot, chk, iproc;
    Integer idx_buf, ldrem[MAXDIM];
    Integer blk_tot = GA[handle].block_total;
    int check1, check2;
    char *pbuf, *prem;

    /* GA uses simple block cyclic data distribution */
    if (GA[handle].block_sl_flag == 0) {
      for(loop=0; loop<num_loops; loop++) {
        /* Loop through all blocks owned by this processor. Decompose
           this loop into a loop over all processors and then a loop
           over all blocks owned by this processor. */
        for (iproc = 0; iproc<GAnproc; iproc++) {

#ifndef __crayx1
          cond = armci_domain_same_id(ARMCI_DOMAIN_SMP,(int)iproc);
          if(loop==0) cond = !cond;
          if(cond) {
#endif
            /* Initialize offset for each processor to zero */
            offset = 0;
            for (idx=iproc; idx<blk_tot; idx += GAnproc) {

              /* get the block corresponding to the virtual processor proc */
              ga_ownsM(handle, idx, blo, bhi);

              /* check to see if this block overlaps with requested block
               * defined by lo and hi */
              chk = 1;
              for (j=0; j<ndim; j++) {
                /* check to see if at least one end point of the interval
                 * represented by blo and bhi falls in the interval
                 * represented by lo and hi */
                check1 = ((blo[j] >= lo[j] && blo[j] <= hi[j]) ||
                          (bhi[j] >= lo[j] && bhi[j] <= hi[j]));
                /* check to see if interval represented by lo and hi
                 * falls entirely within interval represented by blo and bhi */
                check2 = ((lo[j] >= blo[j] && lo[j] <= bhi[j]) &&
                          (hi[j] >= blo[j] && hi[j] <= bhi[j]));
                if (!check1 && !check2) {
                  chk = 0;
                }
              }

              if (chk) {

                /* get the patch of block that overlaps requested region */
                gam_GetBlockPatch(blo,bhi,lo,hi,plo,phi,ndim);

                /* evaluate offset within block */
                last = ndim - 1;
                jtot = 1;
                if (last == 0) ldrem[0] = bhi[0] - blo[0] + 1;
                l_offset = 0;
                for (j=0; j<last; j++) {
                  l_offset += (plo[j]-blo[j])*jtot;
                  ldrem[j] = bhi[j]-blo[j]+1;
                  jtot *= ldrem[j];
                }
                l_offset += (plo[last]-blo[last])*jtot;
                l_offset += offset;

                /* get pointer to data on remote block */
                pinv = idx%GAnproc;
                if (p_handle > 0) {
                  pinv = PGRP_LIST[p_handle].inv_map_proc_list[pinv];
                }
                prem =  GA[handle].ptr[pinv]+l_offset*GA[handle].elemsize;

                gam_ComputePatchIndex(ndim, lo, plo, ld, &idx_buf);
                pbuf = size*idx_buf + (char*)buf;        

                gam_ComputeCount(ndim, plo, phi, count); 
                /* scale number of rows by element size */
                count[0] *= size; 
                gam_setstride(ndim, size, ld, ldrem, stride_rem, stride_loc);

                proc = pinv;
                if(GA_fence_set)fence_array[proc]=1;

#ifndef NO_GA_STATS	    
                if(proc == GAme){
                  gam_CountElems(ndim, plo, phi, &elems);
                  GAbytes.putloc += (double)size*elems;
                }
#endif

                /*casting what ganb_get_armci_handle function returns to armci_hdl is 
                  very crucial here as on 64 bit platforms, pointer is 64 bits where 
                  as temporary is only 32 bits*/ 
#ifdef __crayx1
/*                 ARMCI_PutS(pbuf,stride_loc,prem,stride_rem,count,ndim-1,proc); */
		ngai_puts(buf,pbuf,stride_loc,prem,stride_rem,count,ndim-1,proc, field_off, field_size, size);
#else
                if(nbhandle) {
/*                   ARMCI_NbPutS(pbuf, stride_loc, prem, stride_rem, count, ndim -1, */
/*                       proc,(armci_hdl_t*)get_armci_nbhandle(nbhandle)); */
		  ngai_nbputs(buf,pbuf, stride_loc, prem, stride_rem, count, ndim -1,
			      proc,field_off, field_size, size, 
			      (armci_hdl_t*)get_armci_nbhandle(nbhandle));
		}
                else {
                  /* do blocking put for local processes. If all processes
                     are remote processes then do blocking put for the last one */
                     /*
                  if(loop==1)
                    */
/*                     ARMCI_PutS(pbuf,stride_loc,prem,stride_rem,count,ndim-1,proc); */
		  ngai_puts(buf,pbuf,stride_loc,prem,stride_rem,count,ndim-1,proc, field_off, field_size, size);
                     /*
                  else {
                    ARMCI_NbPutS(pbuf,stride_loc,prem,stride_rem,count, ndim-1,
                        proc,(armci_hdl_t*)get_armci_nbhandle(&ga_nbhandle));
                  }
                    */
                }
#endif
              }
              /* evaluate offset for block idx */
              jtot = 1;
              for (j=0; j<ndim; j++) {
                jtot *= bhi[j]-blo[j]+1;
              }
              offset += jtot;
            }
          }
        }
      }
    } else {
    /* GA uses ScaLAPACK block cyclic data distribution */
      Integer proc_index[MAXDIM], index[MAXDIM];
      Integer itmp;
#if COMPACT_SCALAPACK
#else
      Integer /*blk_size[MAXDIM],*/ blk_num[MAXDIM], blk_dim[MAXDIM];
      Integer blk_inc[MAXDIM], blk_jinc;
      Integer blk_ld[MAXDIM],hlf_blk[MAXDIM];
      /*C_Integer *num_blocks;*/
      C_Integer *block_dims;
      int *proc_grid;

      /* Calculate some properties associated with data distribution */
      proc_grid = GA[handle].nblock;
      /*num_blocks = GA[handle].num_blocks;*/
      block_dims = GA[handle].block_dims;
      for (j=0; j<ndim; j++)  {
        blk_dim[j] = block_dims[j]*proc_grid[j];
        blk_num[j] = GA[handle].dims[j]/blk_dim[j];
        /*blk_size[j] = block_dims[j]*blk_num[j];*/
        blk_inc[j] = GA[handle].dims[j]-blk_num[j]*blk_dim[j];
        blk_ld[j] = blk_num[j]*block_dims[j];
        hlf_blk[j] = blk_inc[j]/block_dims[j];
      }
#endif
      for(loop=0; loop<num_loops; loop++) {
        /* Loop through all blocks owned by this processor. Decompose
           this loop into a loop over all processors and then a loop
           over all blocks owned by this processor. */
        for (iproc = 0; iproc<GAnproc; iproc++) {

#ifndef __crayx1
          cond = armci_domain_same_id(ARMCI_DOMAIN_SMP,(int)iproc);
          if(loop==0) cond = !cond;
          if(cond) {
#endif
            gam_find_proc_indices(handle, iproc, proc_index);
            gam_find_proc_indices(handle, iproc, index);

            /* Initialize offset for each processor to zero */
            offset = 0;
            while (index[ndim-1] < GA[handle].num_blocks[ndim-1]) {

              /* get bounds for current block */
              for (idx = 0; idx < ndim; idx++) {
                blo[idx] = GA[handle].block_dims[idx]*index[idx]+1;
                bhi[idx] = GA[handle].block_dims[idx]*(index[idx]+1);
                if (bhi[idx] > GA[handle].dims[idx]) bhi[idx] = GA[handle].dims[idx];
              }

              /* check to see if this block overlaps with requested block
               * defined by lo and hi */
              chk = 1;
              for (j=0; j<ndim; j++) {
                /* check to see if at least one end point of the interval
                 * represented by blo and bhi falls in the interval
                 * represented by lo and hi */
                check1 = ((blo[j] >= lo[j] && blo[j] <= hi[j]) ||
                          (bhi[j] >= lo[j] && bhi[j] <= hi[j]));
                /* check to see if interval represented by lo and hi
                 * falls entirely within interval represented by blo and bhi */
                check2 = ((lo[j] >= blo[j] && lo[j] <= bhi[j]) &&
                          (hi[j] >= blo[j] && hi[j] <= bhi[j]));
                if (!check1 && !check2) {
                  chk = 0;
                }
              }
              if (chk) {

                /* get the patch of block that overlaps requested region */
                gam_GetBlockPatch(blo,bhi,lo,hi,plo,phi,ndim);

                /* evaluate offset within block */
                last = ndim - 1;
#if COMPACT_SCALAPACK
                jtot = 1;
                if (last == 0) ldrem[0] = bhi[0] - blo[0] + 1;
                l_offset = 0;
                for (j=0; j<last; j++) {
                  l_offset += (plo[j]-blo[j])*jtot;
                  ldrem[j] = bhi[j]-blo[j]+1;
                  jtot *= ldrem[j];
                }
                l_offset += (plo[last]-blo[last])*jtot;
                l_offset += offset;
#else
                l_offset = 0;
                jtot = 1;
                for (j=0; j<last; j++)  {
                  ldrem[j] = blk_ld[j];
                  blk_jinc = GA[handle].dims[j]%block_dims[j];
                  if (blk_inc[j] > 0) {
                    if (proc_index[j]<hlf_blk[j]) {
                      blk_jinc = block_dims[j];
                    } else if (proc_index[j] == hlf_blk[j]) {
                      blk_jinc = blk_inc[j]%block_dims[j];
                      /*
                      if (blk_jinc == 0) {
                        blk_jinc = block_dims[j];
                      }
                      */
                    } else {
                      blk_jinc = 0;
                    }
                  }
                  ldrem[j] += blk_jinc;
                  l_offset += (plo[j]-blo[j]
                            + ((blo[j]-1)/blk_dim[j])*block_dims[j])*jtot;
                  jtot *= ldrem[j];
                }
                l_offset += (plo[last]-blo[last]
                    + ((blo[last]-1)/blk_dim[j])*block_dims[last])*jtot;
#endif

                /* get pointer to data on remote block */
                pinv = iproc;
                if (p_handle > 0) {
                  pinv = PGRP_LIST[p_handle].inv_map_proc_list[pinv];
                }
                prem =  GA[handle].ptr[pinv]+l_offset*GA[handle].elemsize;

                gam_ComputePatchIndex(ndim, lo, plo, ld, &idx_buf);
                pbuf = size*idx_buf + (char*)buf;        

                gam_ComputeCount(ndim, plo, phi, count); 
                /* scale number of rows by element size */
                count[0] *= size; 
                gam_setstride(ndim, size, ld, ldrem, stride_rem, stride_loc);

                proc = pinv;
                if(GA_fence_set)fence_array[proc]=1;

#ifndef NO_GA_STATS	    
                if(proc == GAme){
                  gam_CountElems(ndim, plo, phi, &elems);
                  GAbytes.putloc += (double)size*elems;
                }
#endif

                /*casting what ganb_get_armci_handle function returns to armci_hdl is 
                  very crucial here as on 64 bit platforms, pointer is 64 bits where 
                  as temporary is only 32 bits*/ 
#ifdef __crayx1
/*                 ARMCI_PutS(pbuf,stride_loc,prem,stride_rem,count,ndim-1,proc); */
		ngai_puts(buf,pbuf,stride_loc,prem,stride_rem,count,ndim-1,proc, field_off, field_size, size);
#else
                if(nbhandle) {
/*                   ARMCI_NbPutS(pbuf, stride_loc, prem, stride_rem, count, ndim -1, */
/*                       proc,(armci_hdl_t*)get_armci_nbhandle(nbhandle)); */
		  ngai_nbputs(buf,pbuf, stride_loc, prem, stride_rem, count, ndim -1,
			      proc,field_off, field_size, size, 
			      (armci_hdl_t*)get_armci_nbhandle(nbhandle));
		}
                else {
                  /* do blocking put for local processes. If all processes
                     are remote processes then do blocking put for the last one */
                     /*
                  if(loop==1)
                    */
/*                     ARMCI_PutS(pbuf,stride_loc,prem,stride_rem,count,ndim-1,proc); */
		  ngai_puts(buf,pbuf,stride_loc,prem,stride_rem,count,ndim-1,proc, field_off, field_size, size);
                     /*
                  else {
                    ARMCI_NbPutS(pbuf,stride_loc,prem,stride_rem,count, ndim-1,
                        proc,(armci_hdl_t*)get_armci_nbhandle(&ga_nbhandle));
                  }
                    */
                }
#endif
              }

              /* increment offset to account for all elements on this block */
              itmp = 1;
              for (idx = 0; idx < ndim; idx++) {
                itmp *= (bhi[idx] - blo[idx] + 1);
              }
              offset += itmp;
              
              /* increment block indices to get the next block on processor iproc */
              index[0] += GA[handle].nblock[0];
              for (idx= 0; idx < ndim; idx++) {
                if (index[idx] >= GA[handle].num_blocks[idx] && idx < ndim-1) {
                  index[idx] = proc_index[idx];
                  index[idx+1] += GA[handle].nblock[idx+1];
                }
              }
            }
#ifndef __crayx1
          }
#endif
        }
      }
    }
#if !defined(__crayx1) && !defined(DISABLE_NBOPT)
    if(!nbhandle) nga_wait_internal(&ga_nbhandle);  
#endif
  }

  GA_POP_NAME;
#ifdef PROFILE_OLD
  ga_profile_stop();
#endif
}


/**
 * (Non-blocking) Put an N-dimensional patch of data into a Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_nbput = pnga_nbput
#endif

void pnga_nbput(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld, Integer *nbhandle)
{
  ngai_put_common(g_a,lo,hi,buf,ld,0,-1,nbhandle); 
}

/**
 * (Non-blocking) Put an N-dimensional patch of data into a Global Array and notify the other
                  side with information on another Global Array
 */

#define HANDLES_OUTSTANDING 100
/* Maximum number of outstanding put/notify handles */

typedef struct {
  Integer *orighdl;
  Integer firsthdl;
  Integer elementhdl;
  void *elem_copy;
} gai_putn_hdl_t;

static int putn_handles_initted = 0;
static gai_putn_hdl_t putn_handles[HANDLES_OUTSTANDING];

static int putn_check_single_elem(Integer g_a, Integer *lon, Integer *hin)
{
  int ndims, i;

  ndims = pnga_ndim(g_a);
  for (i = 0; i < ndims; i++)
    if (lon[i] != hin[i])
      return 0;

  return 1;
} /* putn_check_single_elem */

static int putn_find_empty_slot(void)
{
  int i;

  for (i = 0; i < HANDLES_OUTSTANDING; i++)
    if (!putn_handles[i].orighdl)
      return i;

  return -1;
} /* putn_find_empty_slot */

static int putn_intersect_coords(Integer g_a, Integer *lo, Integer *hi, Integer *ecoords)
{
  int ndims, i;

  ndims = pnga_ndim(g_a);

  for (i = 0; i < ndims; i++)
    if ((ecoords[i] < lo[i]) || (ecoords[i] > hi[i]))
      return 0;

  return 1;
} /* putn_intersect_coords */

static int putn_verify_element_in_buf(Integer g_a, Integer *lo, Integer *hi, void *buf,
				      Integer *ld, Integer *ecoords, void *bufn,
				      Integer elemSize)
{
  int i, ndims;
#ifdef HAVE_STDDEF_H
  ptrdiff_t off = (char *)bufn - (char *)buf;
#else
  Integer off = (char *)bufn - (char *)buf;
#endif
  Integer eoff = 0;

  off /= elemSize; /* Offset in terms of elements */

  ndims = pnga_ndim(g_a);
  eoff = ecoords[0] - lo[0];

  /* Check in Fortran ordering */
  for (i = 1; i < ndims; i++)
    eoff += (ecoords[i] - lo[i]) * ld[i - 1];

  return (eoff == (Integer)off); /* Must be the same for a correct notify buffer */
} /* putn_verify_element_in_buf */

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_nbput_notify = pnga_nbput_notify
#endif

void pnga_nbput_notify(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld, Integer g_b, Integer *ecoords, void *bufn, Integer *nbhandle)
{
  Integer ldn[MAXDIM] = { 1 };
  int pos, intersect;

  /* Make sure everything has been initialized */
  if (!putn_handles_initted) {
    memset(putn_handles, 0, sizeof(putn_handles));
    putn_handles_initted = 1;
  }

  pos = putn_find_empty_slot();
  if (pos == -1) /* no empty handles available */
    pnga_error("Too many outstanding put/notify operations!", 0);

  putn_handles[pos].orighdl = nbhandle; /* Store original handle for nbwait_notify */

  if (g_a == g_b)
    intersect = putn_intersect_coords(g_a, lo, hi, ecoords);
  else
    intersect = 0;

  if (!intersect) { /* Simpler case */
    ngai_put_common(g_a, lo, hi, buf, ld, 0, -1, &putn_handles[pos].firsthdl);
    ngai_put_common(g_b, ecoords, ecoords, bufn, ldn, 0, -1, &putn_handles[pos].elementhdl);

    putn_handles[pos].elem_copy = NULL;
  }
  else {
    int ret, i;
    Integer handle = GA_OFFSET + g_a, size;
    void *elem_copy;
    char *elem;

    size = GA[handle].elemsize;
    ret = putn_verify_element_in_buf(g_a, lo, hi, buf, ld, ecoords, bufn, size);

    if (!ret)
      pnga_error("Intersecting buffers, but notify element is not in buffer!", 0);

    elem_copy = malloc(size);
    memcpy(elem_copy, bufn, size);

    elem = bufn;
    for (i = 0; i < size; i++)
      elem[i] += 1; /* Increment each byte by one, safe? */

    putn_handles[pos].elem_copy = elem_copy;

    ngai_put_common(g_a, lo, hi, buf, ld, 0, -1, &putn_handles[pos].firsthdl);
    ngai_put_common(g_a, ecoords, ecoords, elem_copy, ldn, 0, -1,
		    &putn_handles[pos].elementhdl);
  }
} /* pnga_nbput_notify */

/**
 *  Wait for a non-blocking put/notify to complete
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_nbwait_notify = pnga_nbwait_notify
#endif

void pnga_nbwait_notify(Integer *nbhandle)
{
  int i;

  for (i = 0; i < HANDLES_OUTSTANDING; i++)
    if (putn_handles[i].orighdl == nbhandle)
      break;

  if (i >= HANDLES_OUTSTANDING)
    return; /* Incorrect handle used or maybe wait was called multiple times? */

  nga_wait_internal(&putn_handles[i].firsthdl);
  nga_wait_internal(&putn_handles[i].elementhdl);

  if (putn_handles[i].elem_copy) {
    free(putn_handles[i].elem_copy);
    putn_handles[i].elem_copy = NULL;
  }

  putn_handles[i].orighdl = NULL;
} /* pnga_nbwait_notify */

/**
 * Put an N-dimensional patch of data into a Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_put = pnga_put
#endif

void pnga_put(Integer g_a, Integer *lo, Integer *hi, void *buf, Integer *ld)
{
  ngai_put_common(g_a,lo,hi,buf,ld,0,-1,NULL); 
}

/**
 * (Non-blocking) Put a field in a an N-dimensional patch of data into a Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_nbput_field = pnga_nbput_field
#endif

void pnga_nbput_field(Integer g_a, Integer *lo, Integer *hi, Integer foff, Integer fsize, void *buf, Integer *ld, Integer *nbhandle)
{
  ngai_put_common(g_a,lo,hi,buf,ld,foff, fsize, nbhandle); 
}


/**
 * Put a field in a an N-dimensional patch of data into a Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_put_field = pnga_put_field
#endif

void pnga_put_field(Integer g_a, Integer *lo, Integer *hi, Integer foff, Integer fsize, void *buf, Integer *ld)
{
  ngai_put_common(g_a,lo,hi,buf,ld,foff, fsize, NULL); 
}



/*\ A common routine called by both non-blocking and blocking GA Get calls.
\*/
void ngai_get_common(Integer g_a,
                   Integer *lo,
                   Integer *hi,
                   void    *buf,
                   Integer *ld,
		     Integer field_off,
		     Integer field_size,
                   Integer *nbhandle)
{
                   /* g_a:   Global array handle
                      lo[]:  Array of lower indices of patch of global array
                      hi[]:  Array of upper indices of patch of global array
                      buf[]: Local buffer that array patch will be copied into
                      ld[]:  Array of physical ndim-1 dimensions of local buffer */

  Integer  p, np, handle=GA_OFFSET + g_a;
  Integer  idx, elems, size, p_handle;
  int proc, ndim, loop, cond;
  int num_loops=2; /* 1st loop for remote procs; 2nd loop for local procs */
  Integer use_blocks;
  Integer n_rstrctd;
  Integer *rank_rstrctd;
#if defined(__crayx1) || defined(DISABLE_NBOPT)
#else
  Integer ga_nbhandle;
  int counter=0;
#endif

  int _stride_rem[MAXDIM+1], _stride_loc[MAXDIM+1], _count[MAXDIM+1];
  int *stride_rem=&_stride_rem[1], *stride_loc=&_stride_loc[1], *count=&_count[1];

  GA_PUSH_NAME("nga_get_common");

  ga_check_handleM(g_a, "nga_get_common");

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  use_blocks = GA[handle].block_flag;
  p_handle = (Integer)GA[handle].p_handle;
  n_rstrctd = (Integer)GA[handle].num_rstrctd;
  rank_rstrctd = GA[handle].rank_rstrctd;

  /*initial  stride portion for field-wise operations*/
  _stride_rem[0] = size;
  _stride_loc[0] = field_size;
  _count[0] = field_size;

  if (!use_blocks) {

    /* Locate the processors containing some portion of the patch
       specified by lo and hi and return the results in _ga_map,
       GA_proclist, and np. GA_proclist contains a list of processors
       containing some portion of the patch, _ga_map contains
       the lower and upper indices of the portion of the patch held
       by a given processor, and np contains the total number of
       processors that contain some portion of the patch.
     */
    if(!pnga_locate_region(g_a, lo, hi, _ga_map, GA_proclist, &np ))
      ga_RegionError(pnga_ndim(g_a), lo, hi, g_a);

    /* get total size of patch */
#ifndef NO_GA_STATS
    gam_CountElems(ndim, lo, hi, &elems);
    GAbytes.gettot += (double)size*elems;
    GAstat.numget++;
    GAstat.numget_procs += np;
#endif

    if(nbhandle)ga_init_nbhandle(nbhandle);
#if !defined(__crayx1) && !defined(DISABLE_NBOPT)
    else ga_init_nbhandle(&ga_nbhandle);
#endif

#ifdef PROFILE_OLD
    ga_profile_start((int)handle, (long)size*elems, ndim, lo, hi,
        ENABLE_PROFILE_GET);
#endif

    gaPermuteProcList(np);

#if !defined(__crayx1) && !defined(DISABLE_NBOPT)
    for(loop=0; loop<num_loops; loop++) {
      __CRAYX1_PRAGMA("_CRI novector");
#endif
      for(idx=0; idx< np; idx++){
        Integer ldrem[MAXDIM];
        Integer idx_buf, *plo, *phi;
        char *pbuf, *prem;

        p = (Integer)ProcListPerm[idx];
        proc = (int)GA_proclist[p];
        if (p_handle >= 0) {
          proc = PGRP_LIST[p_handle].inv_map_proc_list[proc];
        }
#ifdef PERMUTE_PIDS
        if(GA_Proc_list) proc = GA_inv_Proc_list[proc];
#endif

        /* check if it is local to SMP */
#if !defined(__crayx1) && !defined(DISABLE_NBOPT)
        cond = armci_domain_same_id(ARMCI_DOMAIN_SMP,(int)proc);
        if(loop==0) cond = !cond;
        if(cond) {
#endif

          /* Find  visible portion of patch held by processor p and
             return the result in plo and phi. Also get actual processor
             index corresponding to p and store the result in proc. */
          gam_GetRangeFromMap(p, ndim, &plo, &phi);
          proc = (int)GA_proclist[p];

          /* get pointer prem to location indexed by plo. Also get
             leading physical dimensions in memory in ldrem */
          if (n_rstrctd == 0) {
            gam_Location(proc,handle, plo, &prem, ldrem);
          } else {
            gam_Location(rank_rstrctd[proc],handle, plo, &prem, ldrem);
          }

          /* find the right spot in the user buffer for the point
             subscripted by plo given that the corner of the user
             buffer is subscripted by lo */
          gam_ComputePatchIndex(ndim, lo, plo, ld, &idx_buf);
          pbuf = size*idx_buf + (char*)buf;

          /* compute number of elements in each dimension and store the
             result in count */
          gam_ComputeCount(ndim, plo, phi, count);

          /* Scale first element in count by element size. The ARMCI_GetS
             routine uses this convention to figure out memory sizes.*/
          count[0] *= size; 

          /* Return strides for memory containing global array on remote
             processor indexed by proc (stride_rem) and for local buffer
             buf (stride_loc) */
          gam_setstride(ndim, size, ld, ldrem, stride_rem, stride_loc);

          if (p_handle >= 0) {
            proc = (int)GA_proclist[p];
            /* BJP */
            proc = (int)PGRP_LIST[p_handle].inv_map_proc_list[proc];
          }
#ifdef PERMUTE_PIDS
          if(GA_Proc_list) proc = GA_inv_Proc_list[proc];
#endif

#ifndef NO_GA_STATS	    
          if(proc == GAme){
            gam_CountElems(ndim, plo, phi, &elems);
            GAbytes.getloc += (double)size*elems;
          }
#endif
#if defined(__crayx1) || defined(DISABLE_NBOPT)
/*           ARMCI_GetS(prem,stride_rem,pbuf,stride_loc,count,ndim-1,proc); */
          ngai_gets(buf,prem,stride_rem,pbuf,stride_loc,count,ndim-1,proc, field_off, field_size, size);
#else
          if(nbhandle)  {
/*             ARMCI_NbGetS(prem, stride_rem, pbuf, stride_loc, count, ndim -1, */
/*                 proc,(armci_hdl_t*)get_armci_nbhandle(nbhandle)); */
            ngai_nbgets(buf,prem, stride_rem, pbuf, stride_loc, count, ndim -1,
			proc,field_off, field_size, size,
			(armci_hdl_t*)get_armci_nbhandle(nbhandle));
	  } else {
            if((loop==0 && counter==(int)np-1) || loop==1) {
/*               ARMCI_GetS(prem,stride_rem,pbuf,stride_loc,count,ndim-1,proc); */
	      ngai_gets(buf,prem,stride_rem,pbuf,stride_loc,count,ndim-1,proc, field_off, field_size, size);
	    } else {
              ++counter;
/*               ARMCI_NbGetS(prem,stride_rem,pbuf,stride_loc,count,ndim-1, */
/*                   proc,(armci_hdl_t*)get_armci_nbhandle(&ga_nbhandle)); */
	      ngai_nbgets(buf,prem, stride_rem, pbuf, stride_loc, count, ndim -1,
			  proc,field_off, field_size, size,
			  (armci_hdl_t*)get_armci_nbhandle(nbhandle));
            }
          }
        } /* end if(cond) */
      }
#endif
    }

#if !defined(__crayx1) && !defined(DISABLE_NBOPT)
    if(!nbhandle) nga_wait_internal(&ga_nbhandle);  
#endif
  } else {
    Integer offset, l_offset, last, pinv;
    Integer blo[MAXDIM],bhi[MAXDIM];
    Integer plo[MAXDIM],phi[MAXDIM];
    Integer idx, j, jtot, chk, iproc;
    Integer idx_buf, ldrem[MAXDIM];
    Integer blk_tot = GA[handle].block_total;
    int check1, check2;
    char *pbuf, *prem;

    /* GA uses simple block cyclic data distribution */
    if (GA[handle].block_sl_flag == 0) {
      for(loop=0; loop<num_loops; loop++) {
        /* Loop through all blocks owned by this processor. Decompose
           this loop into a loop over all processors and then a loop
           over all blocks owned by a processor */
        for (iproc = 0; iproc<GAnproc; iproc++) {

#ifndef __crayx1
          cond = armci_domain_same_id(ARMCI_DOMAIN_SMP,(int)iproc);
          if(loop==0) cond = !cond;
          if(cond) {
#endif
            /* Initialize offset for each processor to zero */
            offset = 0;
            for (idx=iproc; idx<blk_tot; idx += GAnproc) {

              /* get the block corresponding to the virtual processor proc */
              ga_ownsM(handle, idx, blo, bhi);

              /* check to see if this block overlaps with requested block
               * defined by lo and hi */
              chk = 1;
              for (j=0; j<ndim; j++) {
                /* check to see if at least one end point of the interval
                 * represented by blo and bhi falls in the interval
                 * represented by lo and hi */
                check1 = ((blo[j] >= lo[j] && blo[j] <= hi[j]) ||
                          (bhi[j] >= lo[j] && bhi[j] <= hi[j]));
                /* check to see if interval represented by lo and hi
                 * falls entirely within interval represented by blo and bhi */
                check2 = ((lo[j] >= blo[j] && lo[j] <= bhi[j]) &&
                          (hi[j] >= blo[j] && hi[j] <= bhi[j]));
                if (!check1 && !check2) {
                  chk = 0;
                }
              }
              if (chk) {
                /* get the patch of block that overlaps requested region */
                gam_GetBlockPatch(blo,bhi,lo,hi,plo,phi,ndim);

                /* evaluate offset within block */
                last = ndim - 1;
                jtot = 1;
                if (last == 0) ldrem[0] = bhi[0] - blo[0] + 1;
                l_offset = 0;
                for (j=0; j<last; j++) {
                  l_offset += (plo[j]-blo[j])*jtot;
                  ldrem[j] = bhi[j]-blo[j]+1;
                  jtot *= ldrem[j];
                }
                l_offset += (plo[last]-blo[last])*jtot;
                l_offset += offset;

                /* get pointer to data on remote block */
                pinv = idx%GAnproc;
                if (p_handle > 0) {
                  pinv = PGRP_LIST[p_handle].inv_map_proc_list[pinv];
                }
                prem =  GA[handle].ptr[pinv]+l_offset*GA[handle].elemsize;

                gam_ComputePatchIndex(ndim, lo, plo, ld, &idx_buf);
                pbuf = size*idx_buf + (char*)buf;        

                gam_ComputeCount(ndim, plo, phi, count); 
                /* scale number of rows by element size */
                count[0] *= size; 
                gam_setstride(ndim, size, ld, ldrem, stride_rem, stride_loc);

                proc = pinv;
                if(GA_fence_set)fence_array[proc]=1;

#ifndef NO_GA_STATS	    
                if(proc == GAme){
                  gam_CountElems(ndim, plo, phi, &elems);
                  GAbytes.putloc += (double)size*elems;
                }
#endif

                /*casting what ganb_get_armci_handle function returns to armci_hdl is 
                  very crucial here as on 64 bit platforms, pointer is 64 bits where 
                  as temporary is only 32 bits*/ 
#ifdef __crayx1
/*                 ARMCI_GetS(prem,stride_rem,pbuf,stride_loc,count,ndim-1,proc); */
		ngai_gets(buf,prem,stride_rem,pbuf,stride_loc,count,ndim-1,proc, field_off, field_size, size);
#else
                if(nbhandle) {
/*                   ARMCI_NbGetS(prem, stride_rem, pbuf, stride_loc, count, ndim -1, */
/*                       proc,(armci_hdl_t*)get_armci_nbhandle(nbhandle)); */
		  ngai_nbgets(buf,prem, stride_rem, pbuf, stride_loc, count, ndim -1,
			      proc,field_off, field_size, size,
			      (armci_hdl_t*)get_armci_nbhandle(nbhandle));
                } else {
                  /* do blocking put for local processes. If all processes
                     are remote processes then do blocking put for the last one */
                     /*
                  if(loop==1) 
                    */
/*                     ARMCI_GetS(prem,stride_rem,pbuf,stride_loc,count,ndim-1,proc); */
		  ngai_gets(buf,prem,stride_rem,pbuf,stride_loc,count,ndim-1,proc, field_off, field_size, size);
                     /*
                  else {
                    ARMCI_NbGetS(prem,stride_rem,pbuf,stride_loc,count, ndim-1,
                        proc,(armci_hdl_t*)get_armci_nbhandle(&ga_nbhandle));
                  }
                    */
                }
#endif
              }
              /* evaluate size of  block idx and use it to increment offset */
              jtot = 1;
              for (j=0; j<ndim; j++) {
                jtot *= bhi[j]-blo[j]+1;
              }
              offset += jtot;
            }
#ifndef __crayx1
          }
#endif
        }
      }
    } else {
      /* GA uses ScaLAPACK block cyclic data distribution */
      Integer proc_index[MAXDIM], index[MAXDIM];
#if COMPACT_SCALAPACK
      Integer itmp;
#else
      Integer /*blk_size[MAXDIM],*/ blk_num[MAXDIM], blk_dim[MAXDIM];
      Integer blk_inc[MAXDIM], blk_jinc;
      Integer blk_ld[MAXDIM], hlf_blk[MAXDIM];
      /*C_Integer *num_blocks;*/
      C_Integer *block_dims;
      int *proc_grid;

      /* Calculate some properties associated with data distribution */
      proc_grid = GA[handle].nblock;
      /*num_blocks = GA[handle].num_blocks;*/
      block_dims = GA[handle].block_dims;
      for (j=0; j<ndim; j++)  {
        blk_dim[j] = block_dims[j]*proc_grid[j];
        blk_num[j] = GA[handle].dims[j]/blk_dim[j];
        /*blk_size[j] = block_dims[j]*blk_num[j];*/
        blk_inc[j] = GA[handle].dims[j] - blk_num[j]*blk_dim[j];
        blk_ld[j] = blk_num[j]*block_dims[j];
        hlf_blk[j] = blk_inc[j]/block_dims[j];
      }
#endif
      for(loop=0; loop<num_loops; loop++) {
        /* Loop through all blocks owned by this processor. Decompose
           this loop into a loop over all processors and then a loop
           over all blocks owned by this processor. */
        for (iproc = 0; iproc<GAnproc; iproc++) {

#ifndef __crayx1
          cond = armci_domain_same_id(ARMCI_DOMAIN_SMP,(int)iproc);
          if(loop==0) cond = !cond;
          if(cond) {
#endif
            gam_find_proc_indices(handle, iproc, proc_index);
            gam_find_proc_indices(handle, iproc, index);

            /* Initialize offset for each processor to zero */
            offset = 0;
            while (index[ndim-1] < GA[handle].num_blocks[ndim-1]) {

              /* get bounds for current block */
              for (idx = 0; idx < ndim; idx++) {
                blo[idx] = GA[handle].block_dims[idx]*index[idx]+1;
                bhi[idx] = GA[handle].block_dims[idx]*(index[idx]+1);
                if (bhi[idx] > GA[handle].dims[idx]) bhi[idx] = GA[handle].dims[idx];
              }

              /* check to see if this block overlaps with requested block
               * defined by lo and hi */
              chk = 1;
              for (j=0; j<ndim; j++) {
                /* check to see if at least one end point of the interval
                 * represented by blo and bhi falls in the interval
                 * represented by lo and hi */
                check1 = ((blo[j] >= lo[j] && blo[j] <= hi[j]) ||
                          (bhi[j] >= lo[j] && bhi[j] <= hi[j]));
                /* check to see if interval represented by lo and hi
                 * falls entirely within interval represented by blo and bhi */
                check2 = ((lo[j] >= blo[j] && lo[j] <= bhi[j]) &&
                          (hi[j] >= blo[j] && hi[j] <= bhi[j]));
                if (!check1 && !check2) {
                  chk = 0;
                }
              }
              if (chk) {
                /* get the patch of block that overlaps requested region */
                gam_GetBlockPatch(blo,bhi,lo,hi,plo,phi,ndim);

                /* evaluate offset within block */
                last = ndim - 1;
#if COMPACT_SCALAPACK
                jtot = 1;
                if (last == 0) ldrem[0] = bhi[0] - blo[0] + 1;
                l_offset = 0;
                for (j=0; j<last; j++) {
                  l_offset += (plo[j]-blo[j])*jtot;
                  ldrem[j] = bhi[j]-blo[j]+1;
                  jtot *= ldrem[j];
                }
                l_offset += (plo[last]-blo[last])*jtot;
                l_offset += offset;
#else
                l_offset = 0;
                jtot = 1;
                for (j=0; j<last; j++)  {
                  ldrem[j] = blk_ld[j];
                  blk_jinc = GA[handle].dims[j]%block_dims[j];
                  if (blk_inc[j] > 0) {
                    if (proc_index[j]<hlf_blk[j]) {
                      blk_jinc = block_dims[j];
                    } else if (proc_index[j] == hlf_blk[j]) {
                      blk_jinc = blk_inc[j]%block_dims[j];
                      /*
                      if (blk_jinc == 0) {
                        blk_jinc = block_dims[j];
                      }
                      */
                    } else {
                      blk_jinc = 0;
                    }
                  }
                  ldrem[j] += blk_jinc;
                  l_offset += (plo[j]-blo[j]
                            + ((blo[j]-1)/blk_dim[j])*block_dims[j])*jtot;
                  jtot *= ldrem[j];
                }
                l_offset += (plo[last]-blo[last]
                    + ((blo[last]-1)/blk_dim[j])*block_dims[last])*jtot;
#endif

                /* get pointer to data on remote block */
                pinv = iproc;
                if (p_handle > 0) {
                  pinv = PGRP_LIST[p_handle].inv_map_proc_list[pinv];
                }
                prem =  GA[handle].ptr[pinv]+l_offset*GA[handle].elemsize;

                gam_ComputePatchIndex(ndim, lo, plo, ld, &idx_buf);
                pbuf = size*idx_buf + (char*)buf;        

                gam_ComputeCount(ndim, plo, phi, count); 
                /* scale number of rows by element size */
                count[0] *= size; 
                gam_setstride(ndim, size, ld, ldrem, stride_rem, stride_loc);

                proc = pinv;
                if(GA_fence_set)fence_array[proc]=1;

#ifndef NO_GA_STATS	    
                if(proc == GAme){
                  gam_CountElems(ndim, plo, phi, &elems);
                  GAbytes.putloc += (double)size*elems;
                }
#endif

                /*casting what ganb_get_armci_handle function returns to armci_hdl is 
                  very crucial here as on 64 bit platforms, pointer is 64 bits where 
                  as temporary is only 32 bits*/ 
#ifdef __crayx1
/*                 ARMCI_GetS(prem,stride_rem,pbuf,stride_loc,count,ndim-1,proc); */
		ngai_gets(buf,prem,stride_rem,pbuf,stride_loc,count,ndim-1,proc, field_off, field_size, size);
#else
                if(nbhandle)  {
/*                   ARMCI_NbGetS(prem, stride_rem, pbuf, stride_loc, count, ndim -1, */
/*                       proc,(armci_hdl_t*)get_armci_nbhandle(nbhandle)); */
		  ngai_nbgets(buf,prem, stride_rem, pbuf, stride_loc, count, ndim -1,
			      proc,field_off, field_size, size,
			      (armci_hdl_t*)get_armci_nbhandle(nbhandle));
                } else {
                  /* do blocking put for local processes. If all processes
                     are remote processes then do blocking put for the last one */
                     /*
                  if(loop==1) 
                    */
/*                     ARMCI_GetS(prem,stride_rem,pbuf,stride_loc,count,ndim-1,proc); */
		  ngai_gets(buf,prem,stride_rem,pbuf,stride_loc,count,ndim-1,proc, field_off, field_size, size);
                     /*
                  else {
                    ARMCI_NbGetS(prem,stride_rem,pbuf,stride_loc,count, ndim-1,
                        proc,(armci_hdl_t*)get_armci_nbhandle(&ga_nbhandle));
                  }
                    */
                }
#endif
              }
#if COMPACT_SCALAPACK
              /* increment offset to account for all elements on this block */
              itmp = 1;
              for (idx = 0; idx < ndim; idx++) {
                itmp *= (bhi[idx] - blo[idx] + 1);
              }
              offset += itmp;
#endif
              
              /* increment block indices to get the next block on processor iproc */
              index[0] += GA[handle].nblock[0];
              for (idx= 0; idx < ndim; idx++) {
                if (index[idx] >= GA[handle].num_blocks[idx] && idx < ndim-1) {
                  index[idx] = proc_index[idx];
                  index[idx+1] += GA[handle].nblock[idx+1];
                }
              }
            }
#ifndef __crayx1
          }
#endif
        }
      }

    }
#if !defined(__crayx1) && !defined(DISABLE_NBOPT)
    if(!nbhandle) nga_wait_internal(&ga_nbhandle);  
#endif
  }

  GA_POP_NAME;
#ifdef PROFILE_OLD
  ga_profile_stop();
#endif
}

/**
 * Get an N-dimensional patch of data from a Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_get = pnga_get
#endif

void pnga_get(Integer g_a, Integer *lo, Integer *hi,
              void *buf, Integer *ld)
{
  ngai_get_common(g_a,lo,hi,buf,ld,0,-1,(Integer *)NULL);
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_nbget = pnga_nbget
#endif

void pnga_nbget(Integer g_a, Integer *lo, Integer *hi,
               void *buf, Integer *ld, Integer *nbhandle)
{
  ngai_get_common(g_a,lo,hi,buf,ld,0,-1,nbhandle);
}

/**
 * Get a field in an N-dimensional patch of data from a Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_get_field = pnga_get_field
#endif

void pnga_get_field(Integer g_a, Integer *lo, Integer *hi,Integer foff, Integer fsize,
              void *buf, Integer *ld)
{
  ngai_get_common(g_a,lo,hi,buf,ld,foff,fsize,(Integer *)NULL);
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_nbget_field = pnga_nbget_field
#endif

void pnga_nbget_field(Integer g_a, Integer *lo, Integer *hi,Integer foff, Integer fsize,
               void *buf, Integer *ld, Integer *nbhandle)
{
  ngai_get_common(g_a,lo,hi,buf,ld,foff,fsize,nbhandle);
}

#ifdef __crayx1 
#  pragma _CRI inline ga_get_
#  pragma _CRI inline ngai_get_common
#endif

/**
 *  A common routine called by both non-blocking and blocking GA acc calls.
 */
void ngai_acc_common(Integer g_a,
                   Integer *lo,
                   Integer *hi,
                   void    *buf,
                   Integer *ld,
                   void    *alpha,
                   Integer *nbhandle)
{
  Integer  p, np, handle=GA_OFFSET + g_a;
  Integer  idx, elems, size, type, p_handle, ga_nbhandle;
  int optype=-1, proc, loop, ndim, cond;
  int num_loops=2; /* 1st loop for remote procs; 2nd loop for local procs */
  Integer use_blocks;
  Integer n_rstrctd;
  Integer *rank_rstrctd;


  GA_PUSH_NAME("nga_acc_common");

  ga_check_handleM(g_a, "nga_acc_common");

  size = GA[handle].elemsize;
  type = GA[handle].type;
  ndim = GA[handle].ndim;
  use_blocks = GA[handle].block_flag;
  p_handle = GA[handle].p_handle;
  n_rstrctd = (Integer)GA[handle].num_rstrctd;
  rank_rstrctd = GA[handle].rank_rstrctd;

  if(type==C_DBL) optype= ARMCI_ACC_DBL;
  else if(type==C_FLOAT) optype= ARMCI_ACC_FLT;
  else if(type==C_DCPL)optype= ARMCI_ACC_DCP;
  else if(type==C_SCPL)optype= ARMCI_ACC_CPL;
  else if(type==C_INT)optype= ARMCI_ACC_INT;
  else if(type==C_LONG)optype= ARMCI_ACC_LNG;
  else pnga_error("type not supported",type);

  if (!use_blocks) {
    /* Locate the processors containing some portion of the patch
       specified by lo and hi and return the results in _ga_map,
       GA_proclist, and np. GA_proclist contains a list of processors
       containing some portion of the patch, _ga_map contains
       the lower and upper indices of the portion of the patch held
       by a given processor, and np contains the total number of
       processors that contain some portion of the patch.
     */
    if(!pnga_locate_region(g_a, lo, hi, _ga_map, GA_proclist, &np ))
      ga_RegionError(pnga_ndim(g_a), lo, hi, g_a);

#ifndef NO_GA_STATS
    gam_CountElems(ndim, lo, hi, &elems);
    GAbytes.acctot += (double)size*elems;
    GAstat.numacc++;
    GAstat.numacc_procs += np;
#endif

    if(nbhandle)ga_init_nbhandle(nbhandle);
    else ga_init_nbhandle(&ga_nbhandle);

#ifdef PROFILE_OLD
    ga_profile_start((int)handle, (long)size*elems, ndim, lo, hi,
        ENABLE_PROFILE_ACC);
#endif

    gaPermuteProcList(np);

    for(loop=0; loop<num_loops; loop++) {
      for(idx=0; idx< np; idx++){
        Integer ldrem[MAXDIM];
        int stride_rem[MAXDIM], stride_loc[MAXDIM], count[MAXDIM];
        Integer idx_buf, *plo, *phi;
        char *pbuf, *prem;

        p = (Integer)ProcListPerm[idx];
        proc = (int)GA_proclist[p];
        if (p_handle >= 0) {
          proc = PGRP_LIST[p_handle].inv_map_proc_list[proc];
        }
#ifdef PERMUTE_PIDS
        if(GA_Proc_list) proc = GA_inv_Proc_list[proc];
#endif

        /* check if it is local to SMP */
        cond = armci_domain_same_id(ARMCI_DOMAIN_SMP,(int)proc);
        if(loop==0) cond = !cond;

        if(cond) {
          /* Find  visible portion of patch held by processor p and
             return the result in plo and phi. Also get actual processor
             index corresponding to p and store the result in proc. */
          gam_GetRangeFromMap(p, ndim, &plo, &phi);
          proc = (int)GA_proclist[p];

          /* get pointer prem to location indexed by plo. Also get
             leading physical dimensions in memory in ldrem */
          if (n_rstrctd == 0) {
            gam_Location(proc,handle, plo, &prem, ldrem);
          } else {
            gam_Location(rank_rstrctd[proc],handle, plo, &prem, ldrem);
          }

          /* find the right spot in the user buffer */
          gam_ComputePatchIndex(ndim,lo, plo, ld, &idx_buf);
          pbuf = size*idx_buf + (char*)buf;

          gam_ComputeCount(ndim, plo, phi, count);

          /* scale number of rows by element size */
          count[0] *= size;
          gam_setstride(ndim, size, ld, ldrem, stride_rem, stride_loc);

          if (p_handle >= 0) {
            proc = (int)GA_proclist[p];
            /* BJP */
            proc = (int)PGRP_LIST[p_handle].inv_map_proc_list[proc];
          }

          if(GA_fence_set)fence_array[proc]=1;

#ifdef PERMUTE_PIDS
          if(GA_Proc_list) proc = GA_inv_Proc_list[proc];
#endif
#ifndef NO_GA_STATS
          if(proc == GAme){
            gam_CountElems(ndim, plo, phi, &elems);
            GAbytes.accloc += (double)size*elems;
          }
#endif

          if(nbhandle) 
            ARMCI_NbAccS(optype, alpha, pbuf, stride_loc, prem,
                stride_rem, count, ndim-1, proc,
                (armci_hdl_t*)get_armci_nbhandle(nbhandle));
          else {
#  if 0 /* disabled, as nbacc fails in quadrics */
            if((loop==0 && counter==(int)np-1) || loop==1)
              ARMCI_AccS(optype, alpha, pbuf, stride_loc, prem, stride_rem, 
                  count, ndim-1, proc);
            else {
              ++counter;
              ARMCI_NbAccS(optype, alpha, pbuf, stride_loc, prem, 
                  stride_rem, count, ndim-1, proc,
                  (armci_hdl_t*)get_armci_nbhandle(&ga_nbhandle));
            }
#  else
            ARMCI_AccS(optype, alpha, pbuf, stride_loc, prem, stride_rem,
                count, ndim-1, proc);
#  endif
          }
        } /* end if(cond) */
      }
    }
#if 0
  if(!nbhandle) nga_wait_internal(&ga_nbhandle);
#endif
  } else {
    Integer offset, l_offset, last, pinv;
    Integer blo[MAXDIM],bhi[MAXDIM];
    Integer plo[MAXDIM],phi[MAXDIM];
    Integer idx, j, jtot, chk, iproc;
    Integer idx_buf, ldrem[MAXDIM];
    Integer blk_tot = GA[handle].block_total;
    int check1, check2;
    int stride_rem[MAXDIM], stride_loc[MAXDIM], count[MAXDIM];
    char *pbuf, *prem;

    /* GA uses simple block cyclic data distribution */
    if (GA[handle].block_sl_flag == 0) {
      for(loop=0; loop<num_loops; loop++) {
        /* Loop through all blocks owned by this processor. Decompose
           this loop into a loop over all processors and then a loop
           over all blocks owned by a processor */
        for (iproc = 0; iproc<GAnproc; iproc++) {

          /* check if it is local to SMP */
          cond = armci_domain_same_id(ARMCI_DOMAIN_SMP,(int)iproc);
          if(loop==0) cond = !cond;
          if(cond) {

            /* Initialize offset for each processor to zero */
            offset = 0;
            for (idx=iproc; idx<blk_tot; idx += GAnproc) {

              /* get the block corresponding to the virtual processor proc */
              ga_ownsM(handle, idx, blo, bhi);

              /* check to see if this block overlaps with requested block
               * defined by lo and hi */
              chk = 1;
              for (j=0; j<ndim; j++) {
                /* check to see if at least one end point of the interval
                 * represented by blo and bhi falls in the interval
                 * represented by lo and hi */
                check1 = ((blo[j] >= lo[j] && blo[j] <= hi[j]) ||
                          (bhi[j] >= lo[j] && bhi[j] <= hi[j]));
                /* check to see if interval represented by lo and hi
                 * falls entirely within interval represented by blo and bhi */
                check2 = ((lo[j] >= blo[j] && lo[j] <= bhi[j]) &&
                          (hi[j] >= blo[j] && hi[j] <= bhi[j]));
                if (!check1 && !check2) {
                  chk = 0;
                }
              }
              if (chk) {

                /* get the patch of block that overlaps requested region */
                gam_GetBlockPatch(blo,bhi,lo,hi,plo,phi,ndim);

                /* evaluate offset within block */
                last = ndim - 1;
                jtot = 1;
                if (last == 0) ldrem[0] = bhi[0] - blo[0] + 1;
                l_offset = 0;
                for (j=0; j<last; j++) {
                  l_offset += (plo[j]-blo[j])*jtot;
                  ldrem[j] = bhi[j]-blo[j]+1;
                  jtot *= ldrem[j];
                }
                l_offset += (plo[last]-blo[last])*jtot;
                l_offset += offset;

                /* get pointer to data on remote block */
                pinv = idx%GAnproc;
                if (p_handle > 0) {
                  pinv = PGRP_LIST[p_handle].inv_map_proc_list[pinv];
                }
                prem =  GA[handle].ptr[pinv]+l_offset*GA[handle].elemsize;

                gam_ComputePatchIndex(ndim, lo, plo, ld, &idx_buf);
                pbuf = size*idx_buf + (char*)buf;        

                gam_ComputeCount(ndim, plo, phi, count); 
                /* scale number of rows by element size */
                count[0] *= size; 
                gam_setstride(ndim, size, ld, ldrem, stride_rem, stride_loc);

                proc = pinv;
                if(GA_fence_set)fence_array[proc]=1;

#ifndef NO_GA_STATS	    
                if(proc == GAme){
                  gam_CountElems(ndim, plo, phi, &elems);
                  GAbytes.putloc += (double)size*elems;
                }
#endif
                if(nbhandle) 
                  ARMCI_NbAccS(optype, alpha, pbuf, stride_loc, prem,
                      stride_rem, count, ndim-1, proc,
                      (armci_hdl_t*)get_armci_nbhandle(nbhandle));
                else {
#  if 0 /* disabled, as nbacc fails in quadrics */
                  if((loop==0 && counter==(int)np-1) || loop==1)
                    ARMCI_AccS(optype, alpha, pbuf, stride_loc, prem, stride_rem, 
                        count, ndim-1, proc);
                  else {
                    ++counter;
                    ARMCI_NbAccS(optype, alpha, pbuf, stride_loc, prem, 
                        stride_rem, count, ndim-1, proc,
                        (armci_hdl_t*)get_armci_nbhandle(&ga_nbhandle));
                  }
#  else
                  ARMCI_AccS(optype, alpha, pbuf, stride_loc, prem, stride_rem,
                      count, ndim-1, proc);
#  endif
                }
              }
              /* evaluate offset for block idx */
              jtot = 1;
              for (j=0; j<ndim; j++) {
                jtot *= bhi[j]-blo[j]+1;
              }
              offset += jtot;
            }
          }
        }
      }
    } else {
      /* GA uses ScaLAPACK block cyclic data distribution */
      Integer proc_index[MAXDIM], index[MAXDIM];
#if COMPACT_SCALAPACK
      Integer itmp;
#else
      Integer /*blk_size[MAXDIM],*/ blk_num[MAXDIM], blk_dim[MAXDIM];
      Integer blk_inc[MAXDIM], blk_jinc;
      Integer blk_ld[MAXDIM],hlf_blk[MAXDIM];
      /*C_Integer *num_blocks;*/
      C_Integer *block_dims;
      int *proc_grid;

      /* Calculate some properties associated with data distribution */
      proc_grid = GA[handle].nblock;
      /*num_blocks = GA[handle].num_blocks;*/
      block_dims = GA[handle].block_dims;
      for (j=0; j<ndim; j++)  {
        blk_dim[j] = block_dims[j]*proc_grid[j];
        blk_num[j] = GA[handle].dims[j]/blk_dim[j];
        /*blk_size[j] = block_dims[j]*blk_num[j];*/
        blk_inc[j] = GA[handle].dims[j]-blk_num[j]*blk_dim[j];
        blk_ld[j] = blk_num[j]*block_dims[j];
        hlf_blk[j] = blk_inc[j]/block_dims[j];
      }
#endif
      for(loop=0; loop<num_loops; loop++) {
        /* Loop through all blocks owned by this processor. Decompose
           this loop into a loop over all processors and then a loop
           over all blocks owned by this processor. */
        for (iproc = 0; iproc<GAnproc; iproc++) {

#ifndef __crayx1
          cond = armci_domain_same_id(ARMCI_DOMAIN_SMP,(int)iproc);
          if(loop==0) cond = !cond;
          if(cond) {
#endif
            gam_find_proc_indices(handle, iproc, proc_index);
            gam_find_proc_indices(handle, iproc, index);

            /* Initialize offset for each processor to zero */
            offset = 0;
            while (index[ndim-1] < GA[handle].num_blocks[ndim-1]) {

              /* get bounds for current block */
              for (idx = 0; idx < ndim; idx++) {
                blo[idx] = GA[handle].block_dims[idx]*index[idx]+1;
                bhi[idx] = GA[handle].block_dims[idx]*(index[idx]+1);
                if (bhi[idx] > GA[handle].dims[idx]) bhi[idx] = GA[handle].dims[idx];
              }

              /* check to see if this block overlaps with requested block
               * defined by lo and hi */
              chk = 1;
              for (j=0; j<ndim; j++) {
                /* check to see if at least one end point of the interval
                 * represented by blo and bhi falls in the interval
                 * represented by lo and hi */
                check1 = ((blo[j] >= lo[j] && blo[j] <= hi[j]) ||
                          (bhi[j] >= lo[j] && bhi[j] <= hi[j]));
                /* check to see if interval represented by lo and hi
                 * falls entirely within interval represented by blo and bhi */
                check2 = ((lo[j] >= blo[j] && lo[j] <= bhi[j]) &&
                          (hi[j] >= blo[j] && hi[j] <= bhi[j]));
                if (!check1 && !check2) {
                  chk = 0;
                }
              }
              if (chk) {

                /* get the patch of block that overlaps requested region */
                gam_GetBlockPatch(blo,bhi,lo,hi,plo,phi,ndim);

                /* evaluate offset within block */
                last = ndim - 1;
#if COMPACT_SCALAPACK
                jtot = 1;
                if (last == 0) ldrem[0] = bhi[0] - blo[0] + 1;
                l_offset = 0;
                for (j=0; j<last; j++) {
                  l_offset += (plo[j]-blo[j])*jtot;
                  ldrem[j] = bhi[j]-blo[j]+1;
                  jtot *= ldrem[j];
                }
                l_offset += (plo[last]-blo[last])*jtot;
                l_offset += offset;
#else
                l_offset = 0;
                jtot = 1;
                for (j=0; j<last; j++)  {
                  ldrem[j] = blk_ld[j];
                  blk_jinc = GA[handle].dims[j]%block_dims[j];
                  if (blk_inc[j] > 0) {
                    if (proc_index[j]<hlf_blk[j]) {
                      blk_jinc = block_dims[j];
                    } else if (proc_index[j] == hlf_blk[j]) {
                      blk_jinc = blk_inc[j]%block_dims[j];
                      /*
                      if (blk_jinc == 0) {
                        blk_jinc = block_dims[j];
                      }
                      */
                    } else {
                      blk_jinc = 0;
                    }
                  }
                  ldrem[j] += blk_jinc;
                  l_offset += (plo[j]-blo[j]
                            + ((blo[j]-1)/blk_dim[j])*block_dims[j])*jtot;
                  jtot *= ldrem[j];
                }
                l_offset += (plo[last]-blo[last]
                    + ((blo[last]-1)/blk_dim[j])*block_dims[last])*jtot;
#endif

                /* get pointer to data on remote block */
                pinv = iproc;
                if (p_handle > 0) {
                  pinv = PGRP_LIST[p_handle].inv_map_proc_list[pinv];
                }
                prem =  GA[handle].ptr[pinv]+l_offset*GA[handle].elemsize;

                gam_ComputePatchIndex(ndim, lo, plo, ld, &idx_buf);
                pbuf = size*idx_buf + (char*)buf;        

                gam_ComputeCount(ndim, plo, phi, count); 
                /* scale number of rows by element size */
                count[0] *= size; 
                gam_setstride(ndim, size, ld, ldrem, stride_rem, stride_loc);

                proc = pinv;
                if(GA_fence_set)fence_array[proc]=1;

#ifndef NO_GA_STATS	    
                if(proc == GAme){
                  gam_CountElems(ndim, plo, phi, &elems);
                  GAbytes.putloc += (double)size*elems;
                }
#endif
                if(nbhandle) 
                  ARMCI_NbAccS(optype, alpha, pbuf, stride_loc, prem,
                      stride_rem, count, ndim-1, proc,
                      (armci_hdl_t*)get_armci_nbhandle(nbhandle));
                else {
#  if 0 /* disabled, as nbacc fails in quadrics */
                  if((loop==0 && counter==(int)np-1) || loop==1)
                    ARMCI_AccS(optype, alpha, pbuf, stride_loc, prem, stride_rem, 
                        count, ndim-1, proc);
                  else {
                    ++counter;
                    ARMCI_NbAccS(optype, alpha, pbuf, stride_loc, prem, 
                        stride_rem, count, ndim-1, proc,
                        (armci_hdl_t*)get_armci_nbhandle(&ga_nbhandle));
                  }
#  else
                  ARMCI_AccS(optype, alpha, pbuf, stride_loc, prem, stride_rem,
                      count, ndim-1, proc);
#  endif
                }
              }

              /* increment offset to account for all elements on this block */
#if COMPACT_SCALAPACK
              itmp = 1;
              for (idx = 0; idx < ndim; idx++) {
                itmp *= (bhi[idx] - blo[idx] + 1);
              }
              offset += itmp;
#endif
              
              /* increment block indices to get the next block on processor iproc */
              index[0] += GA[handle].nblock[0];
              for (idx= 0; idx < ndim; idx++) {
                if (index[idx] >= GA[handle].num_blocks[idx] && idx < ndim-1) {
                  index[idx] = proc_index[idx];
                  index[idx+1] += GA[handle].nblock[idx+1];
                }
              }
            }
#ifndef __crayx1
          }
#endif
        }
      }
    }
#if 0
  if(!nbhandle) nga_wait_internal(&ga_nbhandle);
#endif
  }

  GA_POP_NAME;
#ifdef PROFILE_OLD
  ga_profile_stop();
#endif
}

/**
 *  Accumulate operation for an N-dimensional patch of a Global Array
 *       g_a += alpha * patch
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_acc = pnga_acc
#endif

void pnga_acc(Integer g_a,
              Integer *lo,
              Integer *hi,
              void    *buf,
              Integer *ld,
              void    *alpha)
{
    ngai_acc_common(g_a,lo,hi,buf,ld,alpha,NULL);
}

/**
 *  (Non-blocking) Accumulate operation for an N-dimensional patch of a Global Array
 *       g_a += alpha * patch
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_nbacc = pnga_nbacc
#endif

void pnga_nbacc(Integer g_a,
                Integer *lo,
                Integer *hi,
                void    *buf,
                Integer *ld,
                void    *alpha,
                Integer *nbhndl)
{
    ngai_acc_common(g_a,lo,hi,buf,ld,alpha,nbhndl);
}

/**
 * Return a pointer to local data in a Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_access_ptr = pnga_access_ptr
#endif

void pnga_access_ptr(Integer g_a, Integer lo[], Integer hi[],
                      void* ptr, Integer ld[])
{
char *lptr;
Integer  handle = GA_OFFSET + g_a;
Integer  ow,i,p_handle;

   GA_PUSH_NAME("nga_access_ptr");
   p_handle = GA[handle].p_handle;
   if (!pnga_locate(g_a,lo,&ow)) pnga_error("locate top failed",0);
   if (p_handle != -1)
      ow = PGRP_LIST[p_handle].inv_map_proc_list[ow];
   if ((armci_domain_id(ARMCI_DOMAIN_SMP, ow) != armci_domain_my_id(ARMCI_DOMAIN_SMP)) && (ow != GAme)) 
      pnga_error("cannot access top of the patch",ow);
   if (!pnga_locate(g_a,hi, &ow)) pnga_error("locate bottom failed",0);
   if (p_handle != -1)
      ow = PGRP_LIST[p_handle].inv_map_proc_list[ow];
   if ((armci_domain_id(ARMCI_DOMAIN_SMP, ow) != armci_domain_my_id(ARMCI_DOMAIN_SMP)) && (ow != GAme))
      pnga_error("cannot access bottom of the patch",ow);

   for (i=0; i<GA[handle].ndim; i++)
       if(lo[i]>hi[i]) {
           ga_RegionError(GA[handle].ndim, lo, hi, g_a);
       }

   if (p_handle >= 0) {
     ow = PGRP_LIST[p_handle].map_proc_list[ow];
   }
   if (GA[handle].num_rstrctd > 0) {
     ow = GA[handle].rank_rstrctd[ow];
   }
   gam_Location(ow,handle, lo, &lptr, ld);
   *(char**)ptr = lptr; 
   GA_POP_NAME;
}

/*\ RETURN A POINTER TO BEGINNING OF LOCAL DATA BLOCK
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_access_block_grid_ptr = pnga_access_block_grid_ptr
#endif

void pnga_access_block_grid_ptr(Integer g_a, Integer *index, void* ptr, Integer *ld)
            /* g_a: array handle [input]
             * index: subscript of a particular block  [input]
             * ptr: pointer to data in block [output]
             * ld:  array of strides for block data [output]
             */
{
  char *lptr;
  Integer  handle = GA_OFFSET + g_a;
  Integer  i/*, p_handle*/, offset, factor, inode;
  Integer ndim;
  C_Integer *num_blocks, *block_dims;
  int *proc_grid;
  Integer *dims;
  Integer last;
  Integer proc_index[MAXDIM];
  Integer /*blk_size[MAXDIM],*/ blk_num[MAXDIM], blk_dim[MAXDIM], blk_inc[MAXDIM];
  Integer blk_ld[MAXDIM],hlf_blk[MAXDIM],blk_jinc;
#if COMPACT_SCALAPACK
  Integer j, lo, hi;
  Integer lld[MAXDIM], block_count[MAXDIM], loc_block_dims[MAXDIM];
  Integer ldims[MAXDIM];
#endif

  GA_PUSH_NAME("nga_access_block_grid_ptr");
  /*p_handle = GA[handle].p_handle;*/
  if (!GA[handle].block_sl_flag) {
    pnga_error("Array is not using ScaLAPACK data distribution",0);
  }
  proc_grid = GA[handle].nblock;
  num_blocks = GA[handle].num_blocks;
  block_dims = GA[handle].block_dims;
  dims = GA[handle].dims;
  ndim = GA[handle].ndim;
  for (i=0; i<ndim; i++) {
    if (index[i] < 0 || index[i] >= num_blocks[i])
      pnga_error("block index outside allowed values",index[i]);
  }

  /* find out what processor block is located on */
  gam_find_proc_from_sl_indices(handle, inode, index);

  /* get proc indices of processor that owns block */
  gam_find_proc_indices(handle,inode,proc_index)
  last = ndim-1;
 
  /* Find strides of requested block */
#if COMPACT_SCALAPACK
  for (i=0; i<last; i++)  {
    lo = index[i]*block_dims[i]+1;
    hi = (index[i]+1)*block_dims[i];
    if (hi > dims[i]) hi = dims[i]; 
    ld[i] = (hi - lo + 1);
  }
#else
  for (i=0; i<last; i++)  {
    blk_dim[i] = block_dims[i]*proc_grid[i];
    blk_num[i] = GA[handle].dims[i]/blk_dim[i];
    /*blk_size[i] = block_dims[i]*blk_num[i];*/
    blk_inc[i] = GA[handle].dims[i] - blk_num[i]*blk_dim[i];
    blk_ld[i] = blk_num[i]*block_dims[i];
    hlf_blk[i] = blk_inc[i]/block_dims[i];
    ld[i] = blk_ld[i];
    blk_jinc = dims[i]%block_dims[i];
    if (blk_inc[i] > 0) {
      if (proc_index[i]<hlf_blk[i]) {
        blk_jinc = block_dims[i];
      } else if (proc_index[i] == hlf_blk[i]) {
        blk_jinc = blk_inc[i]%block_dims[i];
        /*
        if (blk_jinc == 0) {
          blk_jinc = block_dims[i];
        }
        */
      } else {
        blk_jinc = 0;
      }
    }
    ld[i] += blk_jinc;
  }
#endif
 
#if COMPACT_SCALAPACK
  /* Find dimensions of local block grid stored on processor inode
     and store results in loc_block_dims. Also find the local grid
     index of block relative to local block grid and store result
     in block_count.
     Find physical dimensions of locally held data and store in 
     lld and set values in ldim, which is used to evaluate the
     offset for the requested block. */
  for (i=0; i<ndim; i++) {
    block_count[i] = 0;
    loc_block_dims[i] = 0;
    lld[i] = 0;
    lo = 0;
    hi = -1;
    for (j=proc_index[i]; j<num_blocks[i]; j += proc_grid[i]) {
      lo = j*block_dims[i] + 1;
      hi = (j+1)*block_dims[i];
      if (hi > dims[i]) hi = dims[i]; 
      if (i<last) lld[i] += (hi - lo + 1);
      if (j<index[i]) block_count[i]++;
      loc_block_dims[i]++;
    }
    if (index[i] < num_blocks[i] - 1 || i == 0) {
      ldims[i] = block_dims[i];
    } else {
      lo = index[i]*block_dims[i] + 1;
      hi = (index[i]+1)*block_dims[i];
      if (hi > dims[i]) hi = dims[i];
      ldims[i] = hi - lo + 1;
    }
  }

  /* Evaluate offset for requested block. This algorithm has only been tested in
     2D and is otherwise completely incomprehensible. */
  offset = 0;
  for (i=0; i<ndim; i++) {
    factor = 1;
    for (j=0; j<ndim; j++) {
      if (j<i) {
        factor *= lld[j];
      } else {
        factor *= ldims[j];
        ldims[j] = block_dims[j];
      }
    }
    offset += block_count[i]*factor;
  }
#else
  /* Evalauate offset for block */
  offset = 0;
  factor = 1;
  for (i = 0; i<ndim; i++) {
    offset += ((index[i]-proc_index[i])/proc_grid[i])*block_dims[i]*factor;
    if (i<ndim-1) factor *= ld[i];
  }
#endif

  lptr = GA[handle].ptr[inode]+offset*GA[handle].elemsize;

  *(char**)ptr = lptr; 
  GA_POP_NAME;
}

/**
 *  Return a pointer to the beginning of a local data block in a block-cyclic
 *  data distribution
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_access_block_ptr = pnga_access_block_ptr
#endif

void pnga_access_block_ptr(Integer g_a, Integer idx, void* ptr, Integer *ld)
            /* g_a: array handle [input]
             * idx: block index  [input]
             * ptr: pointer to data in block [output]
             * ld:  array of strides for block data [output]
             */
{
  char *lptr;
  Integer  handle = GA_OFFSET + g_a;
  Integer  i, j/*, p_handle*/, nblocks, offset, tsum, inode;
  Integer ndim, lo[MAXDIM], hi[MAXDIM], index;

  GA_PUSH_NAME("nga_access_block_ptr");
  /*p_handle = GA[handle].p_handle;*/
  nblocks = GA[handle].block_total;
  ndim = GA[handle].ndim;
  index = idx;
  if (index < 0 || index >= nblocks)
    pnga_error("block index outside allowed values",index);

  if (GA[handle].block_sl_flag == 0) {
    offset = 0;
    inode = index%GAnproc;
    for (i=inode; i<index; i += GAnproc) {
      ga_ownsM(handle,i,lo,hi); 
      tsum = 1;
      for (j=0; j<ndim; j++) {
        tsum *= (hi[j]-lo[j]+1);
      }
      offset += tsum;
    }
    lptr = GA[handle].ptr[inode]+offset*GA[handle].elemsize;

    ga_ownsM(handle,index,lo,hi); 
    for (i=0; i<ndim-1; i++) {
      ld[i] = hi[i]-lo[i]+1;
    }
  } else {
    Integer indices[MAXDIM];
    /* find block indices */
    gam_find_block_indices(handle,index,indices);
    /* find pointer */
    pnga_access_block_grid_ptr(g_a, indices, &lptr, ld);
  }
  *(char**)ptr = lptr; 

  GA_POP_NAME;
}

/**
 *  Return a pointer to the beginning of local data on a processor containing
 *  block-cyclic data
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_access_block_segment_ptr = pnga_access_block_segment_ptr
#endif

void pnga_access_block_segment_ptr(Integer g_a, Integer proc, void* ptr, Integer *len)
            /* g_a:  array handle [input]
             * proc: processor for data [input]
             * ptr:  pointer to data start of data on processor [output]
             * len:  length of data contained on processor [output]
             */
{
  char *lptr;
  Integer  handle = GA_OFFSET + g_a;
  /*Integer  p_handle, nblocks;*/
  Integer /*ndim,*/ index;

  GA_PUSH_NAME("ga_access_block_segment_ptr");
  /*p_handle = GA[handle].p_handle;*/
  /*nblocks = GA[handle].block_total;*/
  /*ndim = GA[handle].ndim;*/
  index = proc;
  if (index < 0 || index >= GAnproc)
    pnga_error("processor index outside allowed values",index);

  if (index != GAme)
    pnga_error("Only get accurate number of elements for processor making request",0);
  lptr = GA[handle].ptr[index];

  *len = GA[handle].size/GA[handle].elemsize;
  *(char**)ptr = lptr; 
  GA_POP_NAME;
}

/**
 * Provide access to a patch of a Global Array using an index
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_access_idx = pnga_access_idx
#endif

void pnga_access_idx(Integer g_a, Integer lo[], Integer hi[],
                     AccessIndex* index, Integer ld[])
{
char     *ptr;
Integer  handle = GA_OFFSET + g_a;
Integer  ow,i,p_handle;
unsigned long    elemsize;
unsigned long    lref=0, lptr;

   GA_PUSH_NAME("nga_access");
   p_handle = GA[handle].p_handle;
   if(!pnga_locate(g_a,lo,&ow))pnga_error("locate top failed",0);
   if (p_handle != -1)
      ow = PGRP_LIST[p_handle].inv_map_proc_list[ow];
   if ((armci_domain_id(ARMCI_DOMAIN_SMP, ow) != armci_domain_my_id(ARMCI_DOMAIN_SMP)) && (ow != GAme)) 
      pnga_error("cannot access top of the patch",ow);
   if(!pnga_locate(g_a,hi, &ow))pnga_error("locate bottom failed",0);
   if (p_handle != -1)
      ow = PGRP_LIST[p_handle].inv_map_proc_list[ow];
   if ((armci_domain_id(ARMCI_DOMAIN_SMP, ow) != armci_domain_my_id(ARMCI_DOMAIN_SMP)) && (ow != GAme)) 
      pnga_error("cannot access bottom of the patch",ow);

   for (i=0; i<GA[handle].ndim; i++)
       if(lo[i]>hi[i]) {
           ga_RegionError(GA[handle].ndim, lo, hi, g_a);
       }

   
   if (p_handle != -1)
      ow = PGRP_LIST[p_handle].map_proc_list[ow];

   gam_Location(ow,handle, lo, &ptr, ld);

   /*
    * return patch address as the distance elements from the reference address
    *
    * .in Fortran we need only the index to the type array: dbl_mb or int_mb
    *  that are elements of COMMON in the the mafdecls.h include file
    * .in C we need both the index and the pointer
    */

   elemsize = (unsigned long)GA[handle].elemsize;

   /* compute index and check if it is correct */
   switch (pnga_type_c2f(GA[handle].type)){
     case MT_F_DBL:
        *index = (AccessIndex) ((DoublePrecision*)ptr - DBL_MB);
        lref = (unsigned long)DBL_MB;
        break;

     case MT_F_DCPL:
        *index = (AccessIndex) ((DoubleComplex*)ptr - DCPL_MB);
        lref = (unsigned long)DCPL_MB;
        break;

     case MT_F_SCPL:
        *index = (AccessIndex) ((SingleComplex*)ptr - SCPL_MB);
        lref = (unsigned long)SCPL_MB;
        break;

     case MT_F_INT:
        *index = (AccessIndex) ((Integer*)ptr - INT_MB);
        lref = (unsigned long)INT_MB;
        break;

     case MT_F_REAL:
        *index = (AccessIndex) ((float*)ptr - FLT_MB);
        lref = (unsigned long)FLT_MB;
        break;        
   }

#ifdef BYTE_ADDRESSABLE_MEMORY
   /* check the allignment */
   lptr = (unsigned long)ptr;
   if( lptr%elemsize != lref%elemsize ){ 
       printf("%d: lptr=%lu(%lu) lref=%lu(%lu)\n",(int)GAme,lptr,lptr%elemsize,
                                                    lref,lref%elemsize);
       pnga_error("nga_access: MA addressing problem: base address misallignment",
                 handle);
   }
#endif

   /* adjust index for Fortran addressing */
   (*index) ++ ;
   FLUSH_CACHE;

   GA_POP_NAME;
}

/*\ PROVIDE ACCESS TO AN INDIVIDUAL DATA BLOCK OF A GLOBAL ARRAY
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_access_block_idx = pnga_access_block_idx
#endif

void pnga_access_block_idx(Integer g_a, Integer idx, AccessIndex* index, Integer *ld)
{
char     *ptr;
Integer  handle = GA_OFFSET + g_a;
Integer  /*p_handle,*/ iblock;
unsigned long    elemsize;
unsigned long    lref=0, lptr;

   GA_PUSH_NAME("nga_access_block");
   /*p_handle = GA[handle].p_handle;*/
   iblock = idx;
   if (iblock < 0 || iblock >= GA[handle].block_total)
     pnga_error("block index outside allowed values",iblock);

   pnga_access_block_ptr(g_a,iblock,&ptr,ld);
   /*
    * return patch address as the distance elements from the reference address
    *
    * .in Fortran we need only the index to the type array: dbl_mb or int_mb
    *  that are elements of COMMON in the the mafdecls.h include file
    * .in C we need both the index and the pointer
    */

   elemsize = (unsigned long)GA[handle].elemsize;

   /* compute index and check if it is correct */
   switch (pnga_type_c2f(GA[handle].type)){
     case MT_F_DBL:
        *index = (AccessIndex) ((DoublePrecision*)ptr - DBL_MB);
        lref = (unsigned long)DBL_MB;
        break;

     case MT_F_DCPL:
        *index = (AccessIndex) ((DoubleComplex*)ptr - DCPL_MB);
        lref = (unsigned long)DCPL_MB;
        break;

     case MT_F_SCPL:
        *index = (AccessIndex) ((SingleComplex*)ptr - SCPL_MB);
        lref = (unsigned long)SCPL_MB;
        break;

     case MT_F_INT:
        *index = (AccessIndex) ((Integer*)ptr - INT_MB);
        lref = (unsigned long)INT_MB;
        break;

     case MT_F_REAL:
        *index = (AccessIndex) ((float*)ptr - FLT_MB);
        lref = (unsigned long)FLT_MB;
        break;        
   }

#ifdef BYTE_ADDRESSABLE_MEMORY
   /* check the allignment */
   lptr = (unsigned long)ptr;
   if( lptr%elemsize != lref%elemsize ){ 
       printf("%d: lptr=%lu(%lu) lref=%lu(%lu)\n",(int)GAme,lptr,lptr%elemsize,
                                                    lref,lref%elemsize);
       pnga_error("nga_access: MA addressing problem: base address misallignment",
                 handle);
   }
#endif

   /* adjust index for Fortran addressing */
   (*index) ++ ;
   FLUSH_CACHE;

   GA_POP_NAME;
}

/**
 *  Provide access to an individual data block of a Global Array
 *  with a SCALAPACK type block-cyclic data distribution
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_access_block_grid_idx = pnga_access_block_grid_idx
#endif

void pnga_access_block_grid_idx(Integer g_a, Integer* subscript,
                                AccessIndex *index, Integer *ld)
{
char     *ptr;
Integer  handle = GA_OFFSET + g_a;
Integer  i,ndim/*,p_handle*/;
unsigned long    elemsize;
unsigned long    lref=0, lptr;

   GA_PUSH_NAME("nga_access_block_grid");
   /*p_handle = GA[handle].p_handle;*/
   ndim = GA[handle].ndim;
   for (i=0; i<ndim; i++) 
     if (subscript[i]<0 || subscript[i] >= GA[handle].num_blocks[i]) 
       pnga_error("index outside allowed values",subscript[i]);

   pnga_access_block_grid_ptr(g_a,subscript,&ptr,ld);
   /*
    * return patch address as the distance elements from the reference address
    *
    * .in Fortran we need only the index to the type array: dbl_mb or int_mb
    *  that are elements of COMMON in the the mafdecls.h include file
    * .in C we need both the index and the pointer
    */

   elemsize = (unsigned long)GA[handle].elemsize;

   /* compute index and check if it is correct */
   switch (pnga_type_c2f(GA[handle].type)){
     case MT_F_DBL:
        *index = (AccessIndex) ((DoublePrecision*)ptr - DBL_MB);
        lref = (unsigned long)DBL_MB;
        break;

     case MT_F_DCPL:
        *index = (AccessIndex) ((DoubleComplex*)ptr - DCPL_MB);
        lref = (unsigned long)DCPL_MB;
        break;

     case MT_F_SCPL:
        *index = (AccessIndex) ((SingleComplex*)ptr - SCPL_MB);
        lref = (unsigned long)SCPL_MB;
        break;

     case MT_F_INT:
        *index = (AccessIndex) ((Integer*)ptr - INT_MB);
        lref = (unsigned long)INT_MB;
        break;

     case MT_F_REAL:
        *index = (AccessIndex) ((float*)ptr - FLT_MB);
        lref = (unsigned long)FLT_MB;
        break;        
   }

#ifdef BYTE_ADDRESSABLE_MEMORY
   /* check the allignment */
   lptr = (unsigned long)ptr;
   if( lptr%elemsize != lref%elemsize ){ 
       printf("%d: lptr=%lu(%lu) lref=%lu(%lu)\n",(int)GAme,lptr,lptr%elemsize,
                                                    lref,lref%elemsize);
       pnga_error("nga_access: MA addressing problem: base address misallignment",
                 handle);
   }
#endif

   /* adjust index for Fortran addressing */
   (*index) ++ ;
   FLUSH_CACHE;

   GA_POP_NAME;
}

/*\ PROVIDE ACCESS TO A PATCH OF A GLOBAL ARRAY
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_access_block_segment_idx = pnga_access_block_segment_idx
#endif

void pnga_access_block_segment_idx(Integer g_a, Integer proc,
                                        AccessIndex* index, Integer *len)
{
char     *ptr;
Integer  handle = GA_OFFSET + g_a;
/*Integer  p_handle;*/
unsigned long    elemsize;
unsigned long    lref=0, lptr;

   GA_PUSH_NAME("nga_access_block_segment");
   /*p_handle = GA[handle].p_handle;*/

   /*
    * return patch address as the distance elements from the reference address
    *
    * .in Fortran we need only the index to the type array: dbl_mb or int_mb
    *  that are elements of COMMON in the the mafdecls.h include file
    * .in C we need both the index and the pointer
    */
   pnga_access_block_segment_ptr(g_a, proc, &ptr, len);

   elemsize = (unsigned long)GA[handle].elemsize;

   /* compute index and check if it is correct */
   switch (pnga_type_c2f(GA[handle].type)){
     case MT_F_DBL:
        *index = (AccessIndex) ((DoublePrecision*)ptr - DBL_MB);
        lref = (unsigned long)DBL_MB;
        break;

     case MT_F_DCPL:
        *index = (AccessIndex) ((DoubleComplex*)ptr - DCPL_MB);
        lref = (unsigned long)DCPL_MB;
        break;

     case MT_F_SCPL:
        *index = (AccessIndex) ((SingleComplex*)ptr - SCPL_MB);
        lref = (unsigned long)SCPL_MB;
        break;

     case MT_F_INT:
        *index = (AccessIndex) ((Integer*)ptr - INT_MB);
        lref = (unsigned long)INT_MB;
        break;

     case MT_F_REAL:
        *index = (AccessIndex) ((float*)ptr - FLT_MB);
        lref = (unsigned long)FLT_MB;
        break;        
   }

#ifdef BYTE_ADDRESSABLE_MEMORY
   /* check the allignment */
   lptr = (unsigned long)ptr;
   if( lptr%elemsize != lref%elemsize ){ 
       printf("%d: lptr=%lu(%lu) lref=%lu(%lu)\n",(int)GAme,lptr,lptr%elemsize,
                                                    lref,lref%elemsize);
       pnga_error("nga_access_block_segment: MA addressing problem: base address misallignment",
                 handle);
   }
#endif

   /* adjust index for Fortran addressing */
   (*index) ++ ;
   FLUSH_CACHE;

   GA_POP_NAME;
}

/**
 *  Release access to a patch of a Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_release = pnga_release
#endif

void pnga_release(Integer g_a, Integer *lo, Integer *hi)
{}

/**
 *  Release access and update a patch of a Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_release_update = pnga_release_update
#endif

void pnga_release_update(Integer g_a, Integer *lo, Integer *hi)
{}

/**
 *  Release access to a block in a block-cyclic Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_release_block = pnga_release_block
#endif

void pnga_release_block(Integer g_a, Integer iblock)
{
  /*
  Integer i;
  Integer handle = GA_OFFSET + g_a;
  for (i=0; i<GA[handle].num_rstrctd; i++) {
    printf("p[%d] location: %d rstrctd_list[%d]: %d\n",GAme,*iblock,
        i,GA[handle].rstrctd_list[i]);
  }
        */
}

/**
 *  Release access and update a block in a block-cyclic Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_release_update_block = pnga_release_update_block
#endif

void pnga_release_update_block(Integer g_a, Integer iblock)
{}

/**
 *  Release access to a block in a block-cyclic Global Array with proc grid
 *  (SCALAPack) layout
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_release_block_grid = pnga_release_block_grid
#endif

void pnga_release_block_grid(Integer g_a, Integer *index)
{}

/**
 *  Release access and update a block in a block-cyclic Global Array with
 *  proc grid (SCALAPack) layout
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_release_update_block_grid = pnga_release_update_block_grid
#endif

void pnga_release_update_block_grid(Integer g_a, Integer *index)
{}

/**
 *  Release access to data segment in a block-cyclic Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_release_block_segment = pnga_release_block_segment
#endif

void pnga_release_block_segment(Integer g_a, Integer iproc)
{}

/**
 *  Release access and update a block in a block-cyclic Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_release_update_block_segment = pnga_release_update_block_segment
#endif

void pnga_release_update_block_segment(Integer g_a, Integer iproc)
{}

void gai_scatter_acc_local(Integer g_a, void *v,Integer *i,Integer *j,
                          Integer nv, void* alpha, Integer proc) 
{
void **ptr_src, **ptr_dst;
char *ptr_ref;
Integer ldp, item_size, ilo, ihi, jlo, jhi, type;
Integer lo[2], hi[2];
Integer handle,p_handle,iproc;
armci_giov_t desc;
register Integer k, offset;
int rc=0;
int use_blocks;

  if (nv < 1) return;

  GA_PUSH_NAME("ga_scatter_local");
  handle = GA_OFFSET + g_a;
  p_handle = GA[handle].p_handle;
  use_blocks = GA[handle].block_flag;

  pnga_distribution(g_a, proc, lo, hi);
  ilo = lo[0];
  jlo = lo[1];
  ihi = hi[0];
  jhi = hi[1];

  iproc = proc;
  if (!use_blocks) {
    if (GA[handle].num_rstrctd > 0)
      iproc = GA[handle].rank_rstrctd[iproc];
    gaShmemLocation(iproc, g_a, ilo, jlo, &ptr_ref, &ldp);
  } else {
    /*
    Integer lo[2];
    lo[0] = ilo;
    lo[1] = jlo;
    */
    pnga_access_block_ptr(g_a, iproc, &ptr_ref, &ldp);
    pnga_release_block(g_a, iproc);
    if (GA[handle].block_sl_flag == 0) {
      proc = proc%pnga_nnodes();
    } else {
      Integer index[2];
      gam_find_block_indices(handle, proc, index);
      index[0] = index[0]%GA[handle].nblock[0];
      index[1] = index[1]%GA[handle].nblock[1];
      gam_find_proc_from_sl_indices(handle,proc,index);
    }

  }

  type = GA[handle].type;
  item_size = GAsizeofM(type);

  ptr_src = gai_malloc((int)nv*2*sizeof(void*));
  if(ptr_src==NULL)pnga_error("gai_malloc failed",nv);
  ptr_dst=ptr_src+ nv;

  for(k=0; k< nv; k++){
     if(i[k] < ilo || i[k] > ihi  || j[k] < jlo || j[k] > jhi){
       char err_string[ERR_STR_LEN];
       sprintf(err_string,"proc=%d invalid i/j=(%ld,%ld)>< [%ld:%ld,%ld:%ld]",
               (int)proc, (long)i[k], (long)j[k], (long)ilo, 
               (long)ihi, (long)jlo, (long)jhi); 
       pnga_error(err_string,g_a);
     }

     offset  = (j[k] - jlo)* ldp + i[k] - ilo;
     ptr_dst[k] = ptr_ref + item_size * offset;
     ptr_src[k] = ((char*)v) + k*item_size;
  }
  desc.bytes = (int)item_size;
  desc.src_ptr_array = ptr_src;
  desc.dst_ptr_array = ptr_dst;
  desc.ptr_array_len = (int)nv;

  if(GA_fence_set)fence_array[proc]=1;

#ifdef PERMUTE_PIDS
    if(GA_Proc_list) proc = GA_inv_Proc_list[proc];
#endif

  if (p_handle >= 0) {
    proc = PGRP_LIST[p_handle].inv_map_proc_list[proc];
  }
  if(alpha != NULL) {
    int optype=-1;
    if(type==C_DBL) optype= ARMCI_ACC_DBL;
    else if(type==C_DCPL)optype= ARMCI_ACC_DCP;
    else if(type==C_SCPL)optype= ARMCI_ACC_CPL;
    else if(type==C_INT)optype= ARMCI_ACC_INT;
    else if(type==C_LONG)optype= ARMCI_ACC_LNG;
    else if(type==C_FLOAT)optype= ARMCI_ACC_FLT;  
    else pnga_error("type not supported",type);
    rc= ARMCI_AccV(optype, alpha, &desc, 1, (int)proc);
  }

  if(rc) pnga_error("scatter/_acc failed in armci",rc);

  gai_free(ptr_src);

  GA_POP_NAME;
}

/**
 *  Scatter nv elements of v into a Global Array at locations specified
 *  by arrays i and j
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_scatter2d = pnga_scatter2d
#endif

void pnga_scatter2d(Integer g_a, void *v, Integer *i, Integer *j, Integer nv)
{
    register Integer k;
    Integer kk;
    Integer item_size;
    Integer proc, type=GA[GA_OFFSET + g_a].type;
    Integer nproc, p_handle, iproc;
    Integer subscrpt[2];

    Integer *aproc, naproc; /* active processes and numbers */
    Integer *map;           /* map the active processes to allocated space */
    char *buf1, *buf2;
    Integer handle = g_a + GA_OFFSET;
    
    Integer *count;   /* counters for each process */
    Integer *nelem;   /* number of elements for each process */
    /* source and destination pointers for each process */
    void ***ptr_src, ***ptr_dst; 
    void **ptr_org; /* the entire pointer array */
    armci_giov_t desc;
    Integer *ilo, *ihi, *jlo, *jhi, *ldp, *owner;
    Integer lo[2], hi[2];
    char **ptr_ref;
    int use_blocks;
    Integer num_blocks=0;
    
    if (nv < 1) return;
    
    ga_check_handleM(g_a, "ga_scatter");
    GA_PUSH_NAME("ga_scatter");
    GAstat.numsca++;
    /* determine how many processors are associated with array */
    p_handle = GA[handle].p_handle;
    if (p_handle < 0) {
      nproc = GAnproc;
      if (GA[handle].num_rstrctd > 0) {
        nproc = GA[handle].num_rstrctd;
      }
    } else {
      nproc = PGRP_LIST[p_handle].map_nproc;
    }
    use_blocks = GA[handle].block_flag;
    
    /* allocate temp memory */
    if (!use_blocks) {
      buf1 = gai_malloc((int) (nproc *4 +nv)* (sizeof(Integer)));
      if(buf1 == NULL) pnga_error("gai_malloc failed", 3*nproc);
    } else {
      num_blocks = GA[handle].block_total;
      buf1 = gai_malloc((int) (num_blocks *4 +nv)* (sizeof(Integer)));
      if(buf1 == NULL) pnga_error("gai_malloc failed", 3*num_blocks);
    }
   
    owner = (Integer *)buf1;  
    count = owner+ nv;  
    if (!use_blocks) {
      nelem =  count + nproc;  
      aproc = count + 2 * nproc;  
      map = count + 3 * nproc;  
    } else {
      nelem =  count + num_blocks;  
      aproc = count + 2 * num_blocks;  
      map = count + 3 * num_blocks;  
    }

    /* initialize the counters and nelem */
    if (!use_blocks) {
      for(kk=0; kk<nproc; kk++) {
        count[kk] = 0; nelem[kk] = 0;
      }
    } else {
      for(kk=0; kk<num_blocks; kk++) {
        count[kk] = 0; nelem[kk] = 0;
      }
    }
    
    /* find proc that owns the (i,j) element; store it in temp:  */
    if (GA[handle].num_rstrctd == 0) {
      for(k=0; k< nv; k++) {
        subscrpt[0] = *(i+k);
        subscrpt[1] = *(j+k);
        if(! pnga_locate(g_a, subscrpt, owner+k)){
          char err_string[ERR_STR_LEN];
          sprintf(err_string,"invalid i/j=(%ld,%ld)", (long)i[k], (long)j[k]);
          pnga_error(err_string,g_a);
        }
        iproc = owner[k];
        nelem[iproc]++;
      }
    } else {
      for(k=0; k< nv; k++) {
        subscrpt[0] = *(i+k);
        subscrpt[1] = *(j+k);
        if(! pnga_locate(g_a, subscrpt, owner+k)){
          char err_string[ERR_STR_LEN];
          sprintf(err_string,"invalid i/j=(%ld,%ld)", (long)i[k], (long)j[k]);
          pnga_error(err_string,g_a);
        }
        iproc = GA[handle].rank_rstrctd[owner[k]];
        nelem[iproc]++;
      }
    }

    naproc = 0;
    if (!use_blocks) {
      for(k=0; k<nproc; k++) if(nelem[k] > 0) {
        aproc[naproc] = k;
        map[k] = naproc;
        naproc ++;
      }
    } else {
      for(k=0; k<num_blocks; k++) if(nelem[k] > 0) {
        aproc[naproc] = k;
        map[k] = naproc;
        naproc ++;
      }
    }
    
    GAstat.numsca_procs += naproc;

    buf2 = gai_malloc((int)(2*naproc*sizeof(void **) + 2*nv*sizeof(void *) +
                      5*naproc*sizeof(Integer) + naproc*sizeof(char*)));
    if(buf2 == NULL) pnga_error("gai_malloc failed", naproc);
 
    ptr_src = (void ***)buf2;
    ptr_dst = (void ***)(buf2 + naproc*sizeof(void **));
    ptr_org = (void **)(buf2 + 2*naproc*sizeof(void **));
    ptr_ref = (char **)(buf2+2*naproc*sizeof(void **)+2*nv*sizeof(void *));
    ilo = (Integer *)(((char*)ptr_ref) + naproc*sizeof(char*));
    ihi = ilo + naproc;
    jlo = ihi + naproc;
    jhi = jlo + naproc;
    ldp = jhi + naproc;

    if (!use_blocks) {
      for(kk=0; kk<naproc; kk++) {
        iproc = aproc[kk];
        if (GA[handle].num_rstrctd > 0)
                    iproc = GA[handle].rstrctd_list[iproc];
        pnga_distribution(g_a, iproc, lo, hi);
        ilo[kk] = lo[0];
        jlo[kk] = lo[1];
        ihi[kk] = hi[0];
        jhi[kk] = hi[1];

        /* get address of the first element owned by proc */
        gaShmemLocation(aproc[kk], g_a, ilo[kk], jlo[kk], &(ptr_ref[kk]),
            &(ldp[kk]));
      }
    } else {
      for(kk=0; kk<naproc; kk++) {
        iproc = aproc[kk];
        pnga_distribution(g_a, iproc, lo, hi);
        ilo[kk] = lo[0];
        jlo[kk] = lo[1];
        ihi[kk] = hi[0];
        jhi[kk] = hi[1];

        /* get address of the first element owned by proc */
        pnga_access_block_ptr(g_a, iproc, &(ptr_ref[kk]), &(ldp[kk]));
        pnga_release_block(g_a, iproc);
      }
    }
    
    /* determine limit for message size --  v,i, & j will travel together */
    item_size = GAsizeofM(type);
    GAbytes.scatot += (double)item_size*nv ;
    iproc = owner[GAme];
    GAbytes.scaloc += (double)item_size* nelem[iproc];
    ptr_src[0] = ptr_org; ptr_dst[0] = ptr_org + nv;
    for(k=1; k<naproc; k++) {
        ptr_src[k] = ptr_src[k-1] + nelem[aproc[k-1]];
        ptr_dst[k] = ptr_dst[k-1] + nelem[aproc[k-1]];
    }
    
    for(k=0; k<nv; k++){
        Integer this_count;
        proc = owner[k];
        if (GA[handle].num_rstrctd > 0)
                    proc = GA[handle].rank_rstrctd[owner[k]];
        this_count = count[proc]; 
        count[proc]++;
        proc = map[proc];
        ptr_src[proc][this_count] = ((char*)v) + k * item_size;

        if(i[k] < ilo[proc] || i[k] > ihi[proc]  ||
           j[k] < jlo[proc] || j[k] > jhi[proc]){
          char err_string[ERR_STR_LEN];
          sprintf(err_string,"proc=%d invalid i/j=(%ld,%ld)><[%ld:%ld,%ld:%ld]",
             (int)proc, (long)i[k], (long)j[k], (long)ilo[proc], 
             (long)ihi[proc], (long)jlo[proc], (long)jhi[proc]);
            pnga_error(err_string, g_a);
        }
        ptr_dst[proc][this_count] = ptr_ref[proc] + item_size *
            ((j[k] - jlo[proc])* ldp[proc] + i[k] - ilo[proc]);
    }
    
    /* source and destination pointers are ready for all processes */
    if (!use_blocks) {
      for(k=0; k<naproc; k++) {
        int rc;

        desc.bytes = (int)item_size;
        desc.src_ptr_array = ptr_src[k];
        desc.dst_ptr_array = ptr_dst[k];
        desc.ptr_array_len = (int)nelem[aproc[k]];
        if (p_handle < 0) {
          iproc = aproc[k];
          if (GA[handle].num_rstrctd > 0)
            iproc = GA[handle].rstrctd_list[iproc];
        } else {
          iproc = PGRP_LIST[p_handle].inv_map_proc_list[aproc[k]];
        }

        rc = ARMCI_PutV(&desc, 1, (int)iproc);
        if(rc) pnga_error("scatter failed in armci",rc);
      }
    } else {
      if (GA[handle].block_sl_flag == 0) {
        for(k=0; k<naproc; k++) {
          int rc;

          desc.bytes = (int)item_size;
          desc.src_ptr_array = ptr_src[k];
          desc.dst_ptr_array = ptr_dst[k];
          desc.ptr_array_len = (int)nelem[aproc[k]];
          iproc = aproc[k];
          iproc = iproc%nproc;
          if (p_handle >= 0) {
            iproc = PGRP_LIST[p_handle].inv_map_proc_list[iproc];
          }

          rc = ARMCI_PutV(&desc, 1, (int)iproc);
          if(rc) pnga_error("scatter failed in armci",rc);
        }
      } else {
        Integer index[MAXDIM];
        for(k=0; k<naproc; k++) {
          int rc;

          desc.bytes = (int)item_size;
          desc.src_ptr_array = ptr_src[k];
          desc.dst_ptr_array = ptr_dst[k];
          desc.ptr_array_len = (int)nelem[aproc[k]];
          iproc = aproc[k];
          gam_find_block_indices(handle, iproc, index);
          index[0] = index[0]%GA[handle].nblock[0];
          index[1] = index[1]%GA[handle].nblock[1];
          gam_find_proc_from_sl_indices(handle,iproc,index);
          if (p_handle >= 0) {
            iproc = PGRP_LIST[p_handle].inv_map_proc_list[iproc];
          }

          rc = ARMCI_PutV(&desc, 1, (int)iproc);
          if(rc) pnga_error("scatter failed in armci",rc);
        }
      }
    }

    gai_free(buf2);
    gai_free(buf1);

    GA_POP_NAME;
}

/**
 *  Add elements of v to contents at random locations in Global Array
 *  specified by vectors i and j
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_scatter_acc2d = pnga_scatter_acc2d
#endif

void pnga_scatter_acc2d(g_a, v, i, j, nv, alpha)
     Integer g_a, nv, *i, *j;
     void *v, *alpha;
{
register Integer k;
Integer item_size;
Integer first, nelem, proc, type=GA[GA_OFFSET + g_a].type;
Integer *int_ptr;
Integer subscrpt[2];

  if (nv < 1) return;

  ga_check_handleM(g_a, "ga_scatter_acc");
  GA_PUSH_NAME("ga_scatter_acc");
  GAstat.numsca++;

  int_ptr = (Integer*) ga_malloc(nv, MT_F_INT, "ga_scatter_acc--p");

  /* find proc that owns the (i,j) element; store it in temp: int_ptr */
  for(k=0; k< nv; k++) {
    subscrpt[0] = *(i+k);
    subscrpt[1] = *(j+k);
    if(! pnga_locate(g_a, subscrpt, int_ptr+k)){
         char err_string[ERR_STR_LEN];
         sprintf(err_string,"invalid i/j=(%ld,%ld)", (long)i[k], (long)j[k]);
         pnga_error(err_string,g_a);
    }
  }

  /* determine limit for message size --  v,i, & j will travel together */
  item_size = GAsizeofM(type);
  GAbytes.scatot += (double)item_size*nv ;

  /* Sort the entries by processor */
  ga_sort_scat(&nv, v, i, j, int_ptr, type );

  /* go through the list again executing scatter for each processor */

  first = 0;
  do {
      proc = int_ptr[first];
      nelem = 0;

      /* count entries for proc from "first" to last */
      for(k=first; k< nv; k++){
        if(proc == int_ptr[k]) nelem++;
        else break;
      }

      if(proc == GAme){
             GAbytes.scaloc += (double)item_size* nelem ;
      }

      gai_scatter_acc_local(g_a, ((char*)v)+item_size*first, i+first,
                           j+first, nelem, alpha, proc);
      first += nelem;

  }while (first< nv);

  ga_free(int_ptr);

  GA_POP_NAME;
}

#define SCATTER -99
#define GATHER -98
#define SCATTER_ACC -97

/*\ GATHER OPERATION elements from the global array into v
\*/
void gai_gatscat(int op, Integer g_a, void* v, Integer subscript[],
                 Integer nv, double *locbytes, double* totbytes, void *alpha)
{
    Integer k, handle=g_a+GA_OFFSET;
    int  ndim, item_size, type;
    Integer *proc;
    Integer nproc, p_handle, iproc;

    Integer *aproc, naproc; /* active processes and numbers */
    Integer *map;           /* map the active processes to allocated space */
    char *buf1, *buf2;
    
    Integer *count;        /* counters for each process */
    Integer *nelem;        /* number of elements for each process */
    /* source and destination pointers for each process */
    void ***ptr_src, ***ptr_dst; 
    void **ptr_org; /* the entire pointer array */
    armci_giov_t desc;
    int use_blocks;
    Integer num_blocks=0;
    
    GA_PUSH_NAME("gai_gatscat");

    proc=(Integer *)ga_malloc(nv, MT_F_INT, "ga_gat-p");

    ndim = GA[handle].ndim;
    type = GA[handle].type;
    item_size = GA[handle].elemsize;
    *totbytes += (double)item_size*nv;

    /* determine how many processors are associated with array */
    p_handle = GA[handle].p_handle;
    if (p_handle < 0) {
      nproc = GAnproc;
      if (GA[handle].num_rstrctd > 0) nproc = GA[handle].num_rstrctd;
    } else {
      nproc = PGRP_LIST[p_handle].map_nproc;
    }
    use_blocks = GA[handle].block_flag;

    /* allocate temp memory */
    if (!use_blocks) {
      buf1 = gai_malloc((int) nproc * 4 * (sizeof(Integer)));
      if(buf1 == NULL) pnga_error("gai_malloc failed", 3*nproc);
    } else {
      num_blocks = GA[handle].block_total;
      buf1 = gai_malloc((int) num_blocks * 4 * (sizeof(Integer)));
      if(buf1 == NULL) pnga_error("gai_malloc failed", 3*num_blocks);
    }
    
    count = (Integer *)buf1;
    if (!use_blocks) {
      nelem = (Integer *)(buf1 + nproc * sizeof(Integer));
      aproc = (Integer *)(buf1 + 2 * nproc * sizeof(Integer));
      map = (Integer *)(buf1 + 3 * nproc * sizeof(Integer));
    } else {
      num_blocks = GA[handle].block_total;
      nelem = (Integer *)(buf1 + num_blocks * sizeof(Integer));
      aproc = (Integer *)(buf1 + 2 * num_blocks * sizeof(Integer));
      map = (Integer *)(buf1 + 3 * num_blocks * sizeof(Integer));
    }
    
    /* initialize the counters and nelem */
    if (!use_blocks) {
      for(k=0; k<nproc; k++) count[k] = 0; 
      for(k=0; k<nproc; k++) nelem[k] = 0;
    } else {
      for(k=0; k<num_blocks; k++) count[k] = 0; 
      for(k=0; k<num_blocks; k++) nelem[k] = 0;
    }

    /* get the process id that the element should go and count the
     * number of elements for each process
     */
    if (GA[handle].num_rstrctd == 0) {
      for(k=0; k<nv; k++) {
        if(!pnga_locate(g_a, subscript+k*ndim, proc+k)) {
          gai_print_subscript("invalid subscript",ndim, subscript+k*ndim,"\n");
          pnga_error("failed -element:",k);
        }
        iproc = proc[k];
        nelem[iproc]++;
      }
    } else {
      for(k=0; k<nv; k++) {
        if(!pnga_locate(g_a, subscript+k*ndim, proc+k)) {
          gai_print_subscript("invalid subscript",ndim, subscript+k*ndim,"\n");
          pnga_error("failed -element:",k);
        }
        iproc = GA[handle].rank_rstrctd[proc[k]];
        nelem[iproc]++;
      }
    }

    /* find the active processes (with which transfer data) */
    naproc = 0;
    if (!use_blocks) {
      for(k=0; k<nproc; k++) if(nelem[k] > 0) {
        aproc[naproc] = k;
        map[k] = naproc;
        naproc ++;
      }
    } else {
      for(k=0; k<num_blocks; k++) if(nelem[k] > 0) {
        aproc[naproc] = k;
        map[k] = naproc;
        naproc ++;
      }
    }

    buf2 = gai_malloc((int)(2*naproc*sizeof(void **) + 2*nv*sizeof(void *)));
    if(buf2 == NULL) pnga_error("gai_malloc failed", 2*naproc);
    
    ptr_src = (void ***)buf2;
    ptr_dst = (void ***)(buf2 + naproc * sizeof(void **));
    ptr_org = (void **)(buf2 + 2 * naproc * sizeof(void **));
    
    /* set the pointers as
     *    P0            P1                  P0          P1        
     * ptr_src[0]   ptr_src[1] ...       ptr_dst[0]  ptr_dst[1] ...
     *        \          \                    \          \
     * ptr_org |-------------------------------|---------------------------|
     *         |                nv             |              nv           |
     *         | nelem[0] | nelem[1] |...      | nelem[0] | nelem[1] |...
     */
    ptr_src[0] = ptr_org; ptr_dst[0] = ptr_org + nv;
    for(k=1; k<naproc; k++) {
        ptr_src[k] = ptr_src[k-1] + nelem[aproc[k-1]];
        ptr_dst[k] = ptr_dst[k-1] + nelem[aproc[k-1]];
    }

    *locbytes += (double)item_size* nelem[GAme];
    
/*
#ifdef PERMUTE_PIDS
    if(GA_Proc_list) p = GA_inv_Proc_list[p];
#endif
*/    
    switch(op) { 
      case GATHER:
        /* go through all the elements
         * for process 0: ptr_src[0][0, 1, ...] = subscript + offset0...
         *                ptr_dst[0][0, 1, ...] = v + offset0...
         * for process 1: ptr_src[1][...] ...
         *                ptr_dst[1][...] ...
         */  
        if (!use_blocks) {
          for(k=0; k<nv; k++){
            iproc = proc[k];
            if (GA[handle].num_rstrctd > 0)
              iproc = GA[handle].rank_rstrctd[iproc];
            ptr_dst[map[iproc]][count[iproc]] = ((char*)v) + k * item_size;
            if (p_handle != 0) {
              gam_Loc_ptr(iproc, handle,  (subscript+k*ndim),
                  ptr_src[map[iproc]]+count[iproc]);
            } else {
              gam_Loc_ptr(proc[k], handle,  (subscript+k*ndim),
                  ptr_src[map[iproc]]+count[iproc]);
            }
            count[iproc]++;
          }
        } else {
          Integer lo[MAXDIM], hi[MAXDIM], ld[MAXDIM];
          Integer j, jtot, last, offset;
          for(k=0; k<nv; k++){
            iproc = proc[k];
            ptr_dst[map[iproc]][count[iproc]] = ((char*)v) + k * item_size;
            pnga_distribution(g_a, iproc, lo, hi);
            pnga_access_block_ptr(g_a, iproc, &(ptr_src[map[iproc]][count[iproc]]), ld);
            pnga_release_block(g_a, iproc);
            /* calculate remaining offset */
            offset = 0;
            last = ndim - 1;
            jtot = 1;
            for (j=0; j<last; j++) {
              offset += ((subscript+k*ndim)[j]-lo[j])*jtot;
              jtot *= ld[j];
            }
            offset += ((subscript+k*ndim)[last]-lo[last])*jtot;
            ptr_src[map[iproc]][count[iproc]]=((char *)ptr_src[map[iproc]][count[iproc]])+offset*GA[handle].elemsize;
            count[iproc]++;
          }
        }
        
        /* source and destination pointers are ready for all processes */
        if (!use_blocks) {
          for(k=0; k<naproc; k++) {
            int rc;

            desc.bytes = (int)item_size;
            desc.src_ptr_array = ptr_src[k];
            desc.dst_ptr_array = ptr_dst[k];
            desc.ptr_array_len = (int)nelem[aproc[k]];

            if (p_handle < 0) {
              iproc = aproc[k];
              if (GA[handle].num_rstrctd > 0)
                iproc = GA[handle].rstrctd_list[iproc];
            } else {
              iproc = PGRP_LIST[p_handle].inv_map_proc_list[aproc[k]];
            }
            rc=ARMCI_GetV(&desc, 1, (int)iproc);
            if(rc) pnga_error("gather failed in armci",rc);
          }
        } else {
          if (GA[handle].block_sl_flag == 0) {
            for(k=0; k<naproc; k++) {
              int rc;

              desc.bytes = (int)item_size;
              desc.src_ptr_array = ptr_src[k];
              desc.dst_ptr_array = ptr_dst[k];
              desc.ptr_array_len = (int)nelem[aproc[k]];
              iproc = aproc[k];
              iproc = iproc%nproc;
              if (p_handle >= 0) {
                iproc = PGRP_LIST[p_handle].inv_map_proc_list[iproc];
              }
              rc=ARMCI_GetV(&desc, 1, (int)iproc);
              if(rc) pnga_error("gather failed in armci",rc);
            }
          } else {
            Integer j, index[MAXDIM];
            for(k=0; k<naproc; k++) {
              int rc;

              desc.bytes = (int)item_size;
              desc.src_ptr_array = ptr_src[k];
              desc.dst_ptr_array = ptr_dst[k];
              desc.ptr_array_len = (int)nelem[aproc[k]];
              iproc = aproc[k];
              gam_find_block_indices(handle, iproc, index);
              for (j=0; j<ndim; j++) {
                index[j] = index[j]%GA[handle].nblock[j];
              }
              gam_find_proc_from_sl_indices(handle,iproc,index);
              if (p_handle >= 0) {
                iproc = PGRP_LIST[p_handle].inv_map_proc_list[iproc];
              }
              rc=ARMCI_GetV(&desc, 1, (int)iproc);
              if(rc) pnga_error("gather failed in armci",rc);
            }
          }
        }
        GAstat.numgat_procs += naproc;
        break;
      case SCATTER:
        /* go through all the elements
         * for process 0: ptr_src[0][0, 1, ...] = v + offset0...
         *                ptr_dst[0][0, 1, ...] = subscript + offset0...
         * for process 1: ptr_src[1][...] ...
         *                ptr_dst[1][...] ...
         */
        if (!use_blocks) {
          for(k=0; k<nv; k++){
            iproc = proc[k];
            if (GA[handle].num_rstrctd > 0)
              iproc = GA[handle].rank_rstrctd[iproc];
            ptr_src[map[iproc]][count[iproc]] = ((char*)v) + k * item_size;
            if (p_handle != 0) {
              gam_Loc_ptr(iproc, handle,  (subscript+k*ndim),
                  ptr_dst[map[iproc]]+count[iproc]);
            } else {
              gam_Loc_ptr(proc[k], handle,  (subscript+k*ndim),
                  ptr_dst[map[iproc]]+count[iproc]);
            }
            count[iproc]++;
          }
        } else {
          Integer lo[MAXDIM], hi[MAXDIM], ld[MAXDIM];
          Integer j, jtot, last, offset;
          for(k=0; k<nv; k++){
            iproc = proc[k];
            ptr_src[map[iproc]][count[iproc]] = ((char*)v) + k * item_size;
            pnga_distribution(g_a, iproc, lo, hi);
            pnga_access_block_ptr(g_a, iproc, &(ptr_dst[map[iproc]][count[iproc]]), ld);
            pnga_release_block(g_a, iproc);
            /* calculate remaining offset */
            offset = 0;
            last = ndim - 1;
            jtot = 1;
            for (j=0; j<last; j++) {
              offset += ((subscript+k*ndim)[j]-lo[j])*jtot;
              jtot *= ld[j];
            }
            offset += ((subscript+k*ndim)[last]-lo[last])*jtot;
            ptr_dst[map[iproc]][count[iproc]]=((char*)ptr_dst[map[iproc]][count[iproc]])+offset*GA[handle].elemsize;
            count[iproc]++;
          }
        }

        /* source and destination pointers are ready for all processes */
        if (!use_blocks) {
          for(k=0; k<naproc; k++) {
            int rc;

            desc.bytes = (int)item_size;
            desc.src_ptr_array = ptr_src[k];
            desc.dst_ptr_array = ptr_dst[k];
            desc.ptr_array_len = (int)nelem[aproc[k]];


            if (p_handle < 0) {
              iproc = aproc[k];
              if (GA[handle].num_rstrctd > 0)
                iproc = GA[handle].rstrctd_list[iproc];
            } else {
              iproc = PGRP_LIST[p_handle].inv_map_proc_list[aproc[k]];
            }
            if(GA_fence_set) fence_array[iproc]=1;

            rc=ARMCI_PutV(&desc, 1, (int)iproc);
            if(rc) pnga_error("scatter failed in armci",rc);
          }
        } else {
          if (GA[handle].block_sl_flag == 0) {
            for(k=0; k<naproc; k++) {
              int rc;

              desc.bytes = (int)item_size;
              desc.src_ptr_array = ptr_src[k];
              desc.dst_ptr_array = ptr_dst[k];
              desc.ptr_array_len = (int)nelem[aproc[k]];

              iproc = aproc[k];
              iproc = iproc%nproc;
              if (p_handle >= 0) {
                iproc = PGRP_LIST[p_handle].inv_map_proc_list[iproc];
              }
              if(GA_fence_set) fence_array[iproc]=1;

              rc=ARMCI_PutV(&desc, 1, (int)iproc);
              if(rc) pnga_error("scatter failed in armci",rc);
            }
          } else {
            Integer j, index[MAXDIM];
            for(k=0; k<naproc; k++) {
              int rc;

              desc.bytes = (int)item_size;
              desc.src_ptr_array = ptr_src[k];
              desc.dst_ptr_array = ptr_dst[k];
              desc.ptr_array_len = (int)nelem[aproc[k]];

              iproc = aproc[k];
              gam_find_block_indices(handle, iproc, index);
              for (j=0; j<ndim; j++) {
                index[j] = index[j]%GA[handle].nblock[j];
              }
              gam_find_proc_from_sl_indices(handle,iproc,index);
              if (p_handle >= 0) {
                iproc = PGRP_LIST[p_handle].inv_map_proc_list[iproc];
              }
              if(GA_fence_set) fence_array[iproc]=1;

              rc=ARMCI_PutV(&desc, 1, (int)iproc);
              if(rc) pnga_error("scatter failed in armci",rc);
            }
          }
        }
        GAstat.numsca_procs += naproc;
        break;
      case SCATTER_ACC:
        /* go through all the elements
         * for process 0: ptr_src[0][0, 1, ...] = v + offset0...
         *                ptr_dst[0][0, 1, ...] = subscript + offset0...
         * for process 1: ptr_src[1][...] ...
         *                ptr_dst[1][...] ...
         */
        if (!use_blocks) {
          for(k=0; k<nv; k++){
            iproc = proc[k];
            if (GA[handle].num_rstrctd > 0)
              iproc = GA[handle].rank_rstrctd[iproc];
            ptr_src[map[iproc]][count[iproc]] = ((char*)v) + k * item_size;
            if (p_handle != 0) {
              gam_Loc_ptr(iproc, handle,  (subscript+k*ndim),
                  ptr_dst[map[iproc]]+count[iproc]);
            } else {
              gam_Loc_ptr(proc[k], handle,  (subscript+k*ndim),
                  ptr_dst[map[iproc]]+count[iproc]);
            }
            count[iproc]++;
          }
        } else {
          Integer lo[MAXDIM], hi[MAXDIM], ld[MAXDIM];
          Integer j, jtot, last, offset;
          for(k=0; k<nv; k++){
            iproc = proc[k];
            ptr_src[map[iproc]][count[iproc]] = ((char*)v) + k * item_size;
            pnga_distribution(g_a, iproc, lo, hi);
            pnga_access_block_ptr(g_a, iproc, &(ptr_dst[map[iproc]][count[iproc]]), ld);
            pnga_release_block(g_a, iproc);
            /* calculate remaining offset */
            offset = 0;
            last = ndim - 1;
            jtot = 1;
            for (j=0; j<last; j++) {
              offset += ((subscript+k*ndim)[j]-lo[j])*jtot;
              jtot *= ld[j];
            }
            offset += ((subscript+k*ndim)[last]-lo[last])*jtot;
            ptr_dst[map[iproc]][count[iproc]]=((char*)ptr_dst[map[iproc]][count[iproc]])+offset*GA[handle].elemsize;
            count[iproc]++;
          }
        }

        /* source and destination pointers are ready for all processes */
        if (!use_blocks) {
          for(k=0; k<naproc; k++) {
            int rc=0;

            desc.bytes = (int)item_size;
            desc.src_ptr_array = ptr_src[k];
            desc.dst_ptr_array = ptr_dst[k];
            desc.ptr_array_len = (int)nelem[aproc[k]];

            if (p_handle < 0) {
              iproc = aproc[k];
              if (GA[handle].num_rstrctd > 0)
                iproc = GA[handle].rstrctd_list[iproc];
            } else {
              iproc = PGRP_LIST[p_handle].inv_map_proc_list[aproc[k]];
            }
            if(GA_fence_set) fence_array[iproc]=1;

            if(alpha != NULL) {
              int optype=-1;
              if(type==C_DBL) optype= ARMCI_ACC_DBL;
              else if(type==C_DCPL)optype= ARMCI_ACC_DCP;
              else if(type==C_SCPL)optype= ARMCI_ACC_CPL;
              else if(type==C_INT)optype= ARMCI_ACC_INT;
              else if(type==C_LONG)optype= ARMCI_ACC_LNG;
              else if(type==C_FLOAT)optype= ARMCI_ACC_FLT; 
              else pnga_error("type not supported",type);
              rc= ARMCI_AccV(optype, alpha, &desc, 1, (int)iproc);
            }
            if(rc) pnga_error("scatter_acc failed in armci",rc);
          }
        } else {
          if (GA[handle].block_sl_flag == 0) {
            for(k=0; k<naproc; k++) {
              int rc=0;

              desc.bytes = (int)item_size;
              desc.src_ptr_array = ptr_src[k];
              desc.dst_ptr_array = ptr_dst[k];
              desc.ptr_array_len = (int)nelem[aproc[k]];

              iproc = aproc[k];
              iproc = iproc%nproc;
              if (p_handle >= 0) {
                iproc = PGRP_LIST[p_handle].inv_map_proc_list[iproc];
              }
              if(GA_fence_set) fence_array[iproc]=1;

              if(alpha != NULL) {
                int optype=-1;
                if(type==C_DBL) optype= ARMCI_ACC_DBL;
                else if(type==C_DCPL)optype= ARMCI_ACC_DCP;
                else if(type==C_SCPL)optype= ARMCI_ACC_CPL;
                else if(type==C_INT)optype= ARMCI_ACC_INT;
                else if(type==C_LONG)optype= ARMCI_ACC_LNG;
                else if(type==C_FLOAT)optype= ARMCI_ACC_FLT; 
                else pnga_error("type not supported",type);
                rc= ARMCI_AccV(optype, alpha, &desc, 1, (int)iproc);
              }
              if(rc) pnga_error("scatter_acc failed in armci",rc);
            }
          } else {
            Integer j, index[MAXDIM];
            for(k=0; k<naproc; k++) {
              int rc=0;

              desc.bytes = (int)item_size;
              desc.src_ptr_array = ptr_src[k];
              desc.dst_ptr_array = ptr_dst[k];
              desc.ptr_array_len = (int)nelem[aproc[k]];

              iproc = aproc[k];
              gam_find_block_indices(handle, iproc, index);
              for (j=0; j<ndim; j++) {
                index[j] = index[j]%GA[handle].nblock[j];
              }
              gam_find_proc_from_sl_indices(handle,iproc,index);
              if (p_handle >= 0) {
                iproc = PGRP_LIST[p_handle].inv_map_proc_list[iproc];
              }
              if(GA_fence_set) fence_array[iproc]=1;

              if(alpha != NULL) {
                int optype=-1;
                if(type==C_DBL) optype= ARMCI_ACC_DBL;
                else if(type==C_DCPL)optype= ARMCI_ACC_DCP;
                else if(type==C_SCPL)optype= ARMCI_ACC_CPL;
                else if(type==C_INT)optype= ARMCI_ACC_INT;
                else if(type==C_LONG)optype= ARMCI_ACC_LNG;
                else if(type==C_FLOAT)optype= ARMCI_ACC_FLT; 
                else pnga_error("type not supported",type);
                rc= ARMCI_AccV(optype, alpha, &desc, 1, (int)iproc);
              }
              if(rc) pnga_error("scatter_acc failed in armci",rc);
            }
          }
        }
        break;        
      default: pnga_error("operation not supported",op);
    }

    gai_free(buf2); gai_free(buf1);
    
    ga_free(proc);
    GA_POP_NAME;
}

/**
 *  Preallocate data for internal linked-list arrays for gather/scatter calls.
 *  @param nelems: the maximum number of elements that will be moved in
 *  gather/scatter call using these arrays.
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_alloc_gatscat_buf = pnga_alloc_gatscat_buf
#endif
void pnga_alloc_gatscat_buf(Integer nelems)
{
  Integer nprocs = pnga_nnodes();
  if (GA_prealloc_gatscat)
    pnga_error("Gather/scatter buffers already allocated",nelems);
  GA_prealloc_gatscat = nelems;
  GA_header =(Integer *)ga_malloc(nprocs, MT_F_INT, "ga_gat_header");
  GA_list =(Integer *)ga_malloc(nelems, MT_F_INT, "ga_gat_list");
  GA_elems =(Integer *)ga_malloc(nprocs, MT_F_INT, "ga_gat_nelems");
}

/**
 *  Free up preallocate data for internal linked-list arrays for gather/scatter calls.
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_free_gatscat_buf = pnga_free_gatscat_buf
#endif
void pnga_free_gatscat_buf()
{
  if (!GA_prealloc_gatscat)
    pnga_error("Gather/scatter buffers not allocated",0);
  GA_prealloc_gatscat = 0;
  ga_free(GA_elems);
  ga_free(GA_list);
  ga_free(GA_header);
}

#define gam_c2f_index(index_c, index_f, ndim)        \
{                                                    \
  Integer _i;                                        \
  for (_i = 0; _i < (ndim); _i++) {                  \
    index_f[_i] = (Integer)(index_c[(ndim)-1-_i]+1); \
  }                                                  \
}

/*\ GATHER OPERATION elements from the global array into v
\*/
void gai_gatscat_new(int op, Integer g_a, void* v, void *subscript,
                     Integer c_flag, Integer nv, double *locbytes,
                     double* totbytes, void *alpha)
{
    Integer handle=g_a+GA_OFFSET;
    int  ndim, i, j, k, type, item_size;
    Integer idx, p_handle, num_rstrctd;
    Integer *proc;
    Integer nprocs, me, iproc, tproc, index[MAXDIM];
    Integer lo[MAXDIM], hi[MAXDIM], ld[MAXDIM-1];
    Integer jtot, last, offset;
    Integer *subscript_ptr;

    Integer *header, *list, *nelems;
    char *buf;
    void **ptr_rem, **ptr_loc; 
    int rc, maxlen, use_blocks, use_sl_blocks;
    int *nblock;
    armci_giov_t desc;

    GA_PUSH_NAME("gai_gatscat_new");

    me = pnga_nodeid();
    num_rstrctd = GA[handle].num_rstrctd;
    use_blocks = GA[handle].block_flag;
    use_sl_blocks = GA[handle].block_sl_flag;
    nblock = GA[handle].nblock;

    /* determine how many processors are associated with array */
    p_handle = GA[handle].p_handle;
    if (p_handle < 0) {
      nprocs = GAnproc;
      if (num_rstrctd > 0) nprocs = num_rstrctd;
    } else {
      nprocs = PGRP_LIST[p_handle].map_nproc;
    }

    if (!GA_prealloc_gatscat) {
      header =(Integer *)ga_malloc(nprocs, MT_F_INT, "ga_gat_header");
      list =(Integer *)ga_malloc(nv, MT_F_INT, "ga_gat_list");
      nelems =(Integer *)ga_malloc(nprocs, MT_F_INT, "ga_gat_nelems");
    } else {
      if (GA_prealloc_gatscat < nv)
        pnga_error("Gather/scatter vector exceeds allocation length ",
                   GA_prealloc_gatscat);
      header = (Integer*)GA_header;
      list = (Integer*)GA_list;
      nelems = (Integer*)GA_elems;
    }

    ndim = GA[handle].ndim;
    type = GA[handle].type;
    item_size = GA[handle].elemsize;
    *totbytes = item_size * nv;

    /* Initialize linked list data structures */
    for (i=0; i<nprocs; i++) {
      nelems[i] = 0;
      header[i] = -1;
    }
    for (i=0; i<nv; i++) {
      list[i] = 0;
    }

    /* Set up linked list that partitions elements defined in subscripts array
     * into groups corresponding to the elements home processor
     */
    maxlen = 1;
    for (i=0; i<nv; i++) {
      if (c_flag) {
        gam_c2f_index(((int**)subscript)[i], index, ndim);
        subscript_ptr = index;
      } else {
        subscript_ptr = ((Integer*)subscript)+i*ndim;
      }
      if(!pnga_locate(g_a, subscript_ptr, &idx)) {
        gai_print_subscript("invalid subscript",ndim, subscript_ptr,"\n");
        pnga_error("failed -element:",i);
      }
      if (num_rstrctd > 0) idx = GA[handle].rank_rstrctd[idx];
      /* map from block to processor, if necessary */
      if (use_blocks) {
        if (use_sl_blocks == 0) {
          idx = idx%nprocs;
        } else {
          gam_find_block_indices(handle,idx,index);
          for (j=0; j<ndim; j++) {
            index[j] = index[j]%nblock[j];
          }
          gam_find_proc_from_sl_indices(handle,idx,index);
        }
      }
      nelems[idx]++;
      if (maxlen<nelems[idx]) maxlen=nelems[idx];
      j = header[idx];
      header[idx] = i;
      list[i] = j;
    }
    *locbytes = item_size * nelems[me];

    /* allocate buffers for individual vector calls */
    buf = (char*)gai_malloc(2*maxlen*sizeof(void*));
    ptr_loc = (void**)buf;
    ptr_rem = (void**)(buf+maxlen*sizeof(void*));

    for (iproc=0; iproc<nprocs; iproc++) {
      if (nelems[iproc] > 0) {
        /* loop through linked list to find all data elements associated with
         * the remote processor iproc
         */
        idx = header[iproc];
        j = 0;
        while (idx > -1) {
          if (c_flag) {
            gam_c2f_index(((int**)subscript)[idx], index, ndim);
            subscript_ptr = index;
          } else {
            subscript_ptr = ((Integer*)subscript)+idx*ndim;
          }
          if (!use_blocks) {
          /* gam_Loc_ptr modifies the value of the processor variable for
           * restricted arrays or processor groups, so make a temporary copy
           */
            tproc = iproc;
            gam_Loc_ptr(tproc, handle,  (subscript_ptr),
                ptr_rem+j);
          } else {
            /* TODO: Figure out how to get correct pointers for block cyclic
             * distributions */
            /* find block index from subscript */
            if(!pnga_locate(g_a, subscript_ptr, &tproc)) {
              gai_print_subscript("invalid subscript for block-cyclic array",
                                   ndim, subscript_ptr,"\n");
              pnga_error("failed -element:",idx);
            }
            pnga_distribution(g_a, tproc, lo, hi);
            pnga_access_block_ptr(g_a, tproc, ptr_rem+j, ld);
            pnga_release_block(g_a, tproc);
            offset = 0;
            last = ndim -1;
            jtot = 1;
            for (k=0; k<last; k++) {
              offset += ((subscript_ptr)[k]-lo[k])*jtot;
              jtot *= ld[k];
            }
            offset += ((subscript_ptr)[last]-lo[last])*jtot;
            ptr_rem[j] = (void*)(((char*)ptr_rem[j])+offset*item_size);
          }
          ptr_loc[j] = (void*)(((char*)v) + idx * item_size);
          idx = list[idx];
          j++;
        }
        /* correct remote proc if restricted arrays or processor groups are
         * being used
         */
        if (num_rstrctd > 0) {
          tproc = GA[handle].rstrctd_list[iproc];
        } else {
          if (p_handle < 0) {
            tproc = iproc;
          } else {
            tproc = PGRP_LIST[p_handle].inv_map_proc_list[iproc];
          }
        }
        /* perform vector operation */
        switch(op) { 
          case GATHER:
            desc.bytes = (int)item_size;
            desc.src_ptr_array = ptr_rem;
            desc.dst_ptr_array = ptr_loc;
            desc.ptr_array_len = (int)nelems[iproc];
            rc=ARMCI_GetV(&desc, 1, (int)tproc);
            if(rc) pnga_error("gather failed in armci",rc);
            break;
          case SCATTER:
            desc.bytes = (int)item_size;
            desc.src_ptr_array = ptr_loc;
            desc.dst_ptr_array = ptr_rem;
            desc.ptr_array_len = (int)nelems[iproc];
            rc=ARMCI_PutV(&desc, 1, (int)tproc);
            if(rc) pnga_error("scatter failed in armci",rc);
            break;
          case SCATTER_ACC:
            desc.bytes = (int)item_size;
            desc.src_ptr_array = ptr_loc;
            desc.dst_ptr_array = ptr_rem;
            desc.ptr_array_len = (int)nelems[iproc];
            if(alpha != NULL) {
              int optype=-1;
              if(type==C_DBL) optype= ARMCI_ACC_DBL;
              else if(type==C_DCPL)optype= ARMCI_ACC_DCP;
              else if(type==C_SCPL)optype= ARMCI_ACC_CPL;
              else if(type==C_INT)optype= ARMCI_ACC_INT;
              else if(type==C_LONG)optype= ARMCI_ACC_LNG;
              else if(type==C_FLOAT)optype= ARMCI_ACC_FLT; 
              else pnga_error("type not supported",type);
              rc= ARMCI_AccV(optype, alpha, &desc, 1, (int)tproc);
            }
            if(rc) pnga_error("scatter_acc failed in armci",rc);
            break;
          default: pnga_error("operation not supported",op);
        }
      }
    }
    gai_free(buf);
    if (!GA_prealloc_gatscat) {
      ga_free(nelems);
      ga_free(list);
      ga_free(header);
    }

    GA_POP_NAME;
}

/**
 *  Gather random elements from a global array into local buffer v
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_gather = pnga_gather
#endif

void pnga_gather(Integer g_a, void* v, void *subscript, Integer c_flag, Integer nv)
{

  if (nv < 1) return;
  ga_check_handleM(g_a, "nga_gather");
  GA_PUSH_NAME("nga_gather");
  GAstat.numgat++;

#ifdef USE_GATSCAT_NEW
  gai_gatscat_new(GATHER,g_a,v,subscript,c_flag,nv,&GAbytes.gattot,&GAbytes.gatloc, NULL);
#else
  gai_gatscat(GATHER,g_a,v,subscript,nv,&GAbytes.gattot,&GAbytes.gatloc, NULL);
#endif

  GA_POP_NAME;
}

/**
 *  Scatter elements from a local buffer v into random locations in a Global
 *  Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_scatter = pnga_scatter
#endif

void pnga_scatter(Integer g_a, void* v, void *subscript, Integer c_flag, Integer nv)
{

  if (nv < 1) return;
  ga_check_handleM(g_a, "nga_scatter");
  GA_PUSH_NAME("nga_scatter");
  GAstat.numsca++;

#ifdef USE_GATSCAT_NEW
  gai_gatscat_new(SCATTER,g_a,v,subscript,c_flag,nv,&GAbytes.scatot,&GAbytes.scaloc, NULL);
#else
  gai_gatscat(SCATTER,g_a,v,subscript,nv,&GAbytes.scatot,&GAbytes.scaloc, NULL);
#endif

  GA_POP_NAME;
}

/**
 *  Add elements of v to contents at random locations in Global Array
 *  specified by array subscript
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_scatter_acc = pnga_scatter_acc
#endif

void pnga_scatter_acc(Integer g_a, void* v, void *subscript, Integer c_flag,
                      Integer nv, void *alpha)
{

  if (nv < 1) return;
  ga_check_handleM(g_a, "nga_scatter_acc");
  GA_PUSH_NAME("nga_scatter_acc");
  GAstat.numsca++;

#ifdef USE_GATSCAT_NEW
  gai_gatscat_new(SCATTER_ACC, g_a, v, subscript, c_flag, nv, &GAbytes.scatot,
              &GAbytes.scaloc, alpha);
#else
  gai_gatscat(SCATTER_ACC, g_a, v, subscript, nv, &GAbytes.scatot,
              &GAbytes.scaloc, alpha);
#endif

  GA_POP_NAME;
}

#if 000
void FATR  ga_gather000_(g_a, v, i, j, nv)
     Integer g_a, nv, *i, *j;
     void *v;
{
int k;
Integer *sbar = (Integer*)malloc(2*sizeof(Integer)*  (int)nv);
     if(!sbar) pnga_error("gather:malloc failed",nv);
     for(k=0;k<nv;k++){
          sbar[2*k] = i[k];
          sbar[2*k+1] = j[k];
     }
     nga_gather_(g_a,v,sbar,nv);
     free(sbar);
}
  


/*\ SCATTER OPERATION elements of v into the global array
\*/
void FATR  ga_scatter000_(g_a, v, i, j, nv)
     Integer g_a, nv, *i, *j;
     void *v;
{
int k;
Integer *sbar = (Integer*)malloc(2*sizeof(Integer)* (int) nv);
     if(!sbar) pnga_error("scatter:malloc failed",nv);
     for(k=0;k<nv;k++){
          sbar[2*k] = i[k];
          sbar[2*k+1] = j[k];
     }
     nga_scatter_(g_a,v,sbar,nv);
     free(sbar);
}

#endif


/**
 *  Gather random elements from 2D global array into local buffer v
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_gather2d = pnga_gather2d
#endif

void pnga_gather2d(Integer g_a, void *v, Integer *i, Integer *j,
                         Integer nv)
{
    Integer k, kk, proc, item_size;
    Integer *aproc, naproc; /* active processes and numbers */
    Integer *map;           /* map the active processes to allocated space */
    char *buf1, *buf2;
    Integer handle = g_a + GA_OFFSET;
    Integer nproc, p_handle, iproc;
    
    Integer *count;   /* counters for each process */
    Integer *nelem;   /* number of elements for each process */
    /* source and destination pointers for each process */
    void ***ptr_src, ***ptr_dst; 
    void **ptr_org; /* the entire pointer array */
    armci_giov_t desc;
    Integer *ilo, *ihi, *jlo, *jhi, *ldp, *owner;
    Integer lo[2], hi[2];
    char **ptr_ref;
    int use_blocks;
    Integer num_blocks=0;
    Integer subscrpt[2];
    
    if (nv < 1) return;

    ga_check_handleM(g_a, "ga_gather");
    GA_PUSH_NAME("ga_gather");
    GAstat.numgat++;

    /* determine how many processors are associated with array */
    p_handle = GA[handle].p_handle;
    if (p_handle < 0) {
      nproc = GAnproc;
      if (GA[handle].num_rstrctd > 0) {
        nproc = GA[handle].num_rstrctd;
      }
    } else {
      nproc = PGRP_LIST[p_handle].map_nproc;
    }
    use_blocks = GA[handle].block_flag;

    /* allocate temp memory */
    if (!use_blocks) {
      buf1 = gai_malloc((int)(nproc *4  + nv)*  (sizeof(Integer)));
      if(buf1 == NULL) pnga_error("gai_malloc failed", 3*nproc);
    } else {
      num_blocks = GA[handle].block_total;
      buf1 = gai_malloc((int)(num_blocks *4  + nv)*  (sizeof(Integer)));
      if(buf1 == NULL) pnga_error("gai_malloc failed", 3*num_blocks);
    }
    
    owner = (Integer *)buf1; 
    count = owner+ nv;
    if (!use_blocks) {
      nelem = count + nproc;
      aproc = count + 2 * nproc;
      map =   count + 3 * nproc;
    } else {
      nelem = count + num_blocks;
      aproc = count + 2 * num_blocks;
      map =   count + 3 * num_blocks;
    }
   
    if (!use_blocks) {
      /* initialize the counters and nelem */
      for(kk=0; kk<nproc; kk++) {
        count[kk] = 0; nelem[kk] = 0;
      }
    } else {
      for(kk=0; kk<num_blocks; kk++) {
        count[kk] = 0; nelem[kk] = 0;
      }
    }

    /* find proc or block that owns the (i,j) element; store it in temp: */
    if (GA[handle].num_rstrctd == 0) {
      for(k=0; k< nv; k++) {
        subscrpt[0] = *(i+k);
        subscrpt[1] = *(j+k);
        if(! pnga_locate(g_a, subscrpt, owner+k)){
          char err_string[ERR_STR_LEN];
          sprintf(err_string,"invalid i/j=(%ld,%ld)", (long)i[k], (long)j[k]);
          pnga_error(err_string, g_a);
        }
        iproc = owner[k];
        nelem[iproc]++;
      }
    } else {
      for(k=0; k< nv; k++) {
        subscrpt[0] = *(i+k);
        subscrpt[1] = *(j+k);
        if(! pnga_locate(g_a, subscrpt, owner+k)){
          char err_string[ERR_STR_LEN];
          sprintf(err_string,"invalid i/j=(%ld,%ld)", (long)i[k], (long)j[k]);
          pnga_error(err_string, g_a);
        }
        iproc = GA[handle].rank_rstrctd[owner[k]];
        nelem[iproc]++;
      }
    }

    naproc = 0;
    if (!use_blocks) {
      for(k=0; k<nproc; k++) if(nelem[k] > 0) {
        aproc[naproc] = k;
        map[k] = naproc;
        naproc ++;
      }
    } else {
      for(k=0; k<num_blocks; k++) if(nelem[k] > 0) {
        aproc[naproc] = k;
        map[k] = naproc;
        naproc ++;
      }
    }
    GAstat.numgat_procs += naproc;
    
    buf2 = gai_malloc((int)(2*naproc*sizeof(void **) + 2*nv*sizeof(void *) +
                      5*naproc*sizeof(Integer) + naproc*sizeof(char*)));
    if(buf2 == NULL) pnga_error("gai_malloc failed", naproc);
 
    ptr_src = (void ***)buf2;
    ptr_dst = (void ***)(buf2 + naproc*sizeof(void **));
    ptr_org = (void **)(buf2 + 2*naproc*sizeof(void **));
    ptr_ref = (char **)(buf2+2*naproc*sizeof(void **)+2*nv*sizeof(void *));
    ilo = (Integer *)(((char*)ptr_ref) + naproc*sizeof(char*));
    ihi = ilo + naproc;
    jlo = ihi + naproc;
    jhi = jlo + naproc;
    ldp = jhi + naproc;

    if (!use_blocks) {
      for(kk=0; kk<naproc; kk++) {
        iproc = aproc[kk];
        if (GA[handle].num_rstrctd > 0)
          iproc = GA[handle].rstrctd_list[iproc];
        pnga_distribution(g_a, iproc, lo, hi);
        ilo[kk] = lo[0];
        jlo[kk] = lo[1];
        ihi[kk] = hi[0];
        jhi[kk] = hi[1];

        /* get address of the first element owned by proc */
        gaShmemLocation(aproc[kk], g_a, ilo[kk], jlo[kk], &(ptr_ref[kk]),
            &(ldp[kk]));
      }
    } else {
      for(kk=0; kk<naproc; kk++) {
        iproc = aproc[kk];
        pnga_distribution(g_a, iproc, lo, hi);
        ilo[kk] = lo[0];
        jlo[kk] = lo[1];
        ihi[kk] = hi[0];
        jhi[kk] = hi[1];

        /* get address of the first element owned by proc */
        pnga_access_block_ptr(g_a, iproc, &(ptr_ref[kk]), &(ldp[kk]));
        pnga_release_block(g_a, iproc);
      }
    }
    
    item_size = GA[GA_OFFSET + g_a].elemsize;
    GAbytes.gattot += (double)item_size*nv;
    /*This stuff is probably completely wrong. Doesn't affect performance,
     * just statistics. */
    iproc = owner[GAme];
    GAbytes.gatloc += (double)item_size * nelem[iproc];

    ptr_src[0] = ptr_org; ptr_dst[0] = ptr_org + nv;
    for(k=1; k<naproc; k++) {
        ptr_src[k] = ptr_src[k-1] + nelem[aproc[k-1]];
        ptr_dst[k] = ptr_dst[k-1] + nelem[aproc[k-1]];
    }
    
    for(k=0; k<nv; k++){
        Integer this_count;
        proc = owner[k];
        if (GA[handle].num_rstrctd > 0)
          proc = GA[handle].rank_rstrctd[owner[k]];
        this_count = count[proc]; 
        count[proc]++;
        proc = map[proc]; 
        ptr_dst[proc][this_count] = ((char*)v) + k * item_size;

        if(i[k] < ilo[proc] || i[k] > ihi[proc]  ||
           j[k] < jlo[proc] || j[k] > jhi[proc]){
          char err_string[ERR_STR_LEN];
          sprintf(err_string,"proc=%d invalid i/j=(%ld,%ld)><[%ld:%ld,%ld:%ld]",
                 (int)proc,(long)i[k],(long)j[k],
                 (long)ilo[proc],(long)ihi[proc],
                 (long)jlo[proc], (long)jhi[proc]);
            pnga_error(err_string, g_a);
        }
        ptr_src[proc][this_count] = ptr_ref[proc] + item_size *
            ((j[k] - jlo[proc])* ldp[proc] + i[k] - ilo[proc]);
    }
    /* source and destination pointers are ready for all processes */
    if (!use_blocks) {
      for(k=0; k<naproc; k++) {
        int rc;

        desc.bytes = (int)item_size;
        desc.src_ptr_array = ptr_src[k];
        desc.dst_ptr_array = ptr_dst[k];
        desc.ptr_array_len = (int)nelem[aproc[k]];
        if (p_handle < 0) {
          iproc = aproc[k];
          if (GA[handle].num_rstrctd > 0)
            iproc = GA[handle].rstrctd_list[iproc];
        } else {
          iproc = PGRP_LIST[p_handle].inv_map_proc_list[aproc[k]]; 
        }
        rc=ARMCI_GetV(&desc, 1, (int)iproc);
        if(rc) pnga_error("gather failed in armci",rc);
      }
    } else {
      if (GA[handle].block_sl_flag == 0) {
        for(k=0; k<naproc; k++) {
          int rc;

          desc.bytes = (int)item_size;
          desc.src_ptr_array = ptr_src[k];
          desc.dst_ptr_array = ptr_dst[k];
          desc.ptr_array_len = (int)nelem[aproc[k]];
          iproc = aproc[k];
          iproc = iproc%nproc;
          if (p_handle >= 0) {
            iproc = PGRP_LIST[p_handle].inv_map_proc_list[iproc]; 
          }
          rc=ARMCI_GetV(&desc, 1, (int)iproc);
          if(rc) pnga_error("gather failed in armci",rc);
        }
      } else {
        Integer index[MAXDIM];
        for(k=0; k<naproc; k++) {
          int rc;

          desc.bytes = (int)item_size;
          desc.src_ptr_array = ptr_src[k];
          desc.dst_ptr_array = ptr_dst[k];
          desc.ptr_array_len = (int)nelem[aproc[k]];
          iproc = aproc[k];
          gam_find_block_indices(handle, iproc, index);
          index[0] = index[0]%GA[handle].nblock[0];
          index[1] = index[1]%GA[handle].nblock[1];
          gam_find_proc_from_sl_indices(handle,iproc,index);
          if (p_handle >= 0) {
            iproc = PGRP_LIST[p_handle].inv_map_proc_list[iproc]; 
          }
          rc=ARMCI_GetV(&desc, 1, (int)iproc);
          if(rc) pnga_error("gather failed in armci",rc);
        }
      }
    }

    gai_free(buf2);
    gai_free(buf1);
    GA_POP_NAME;
}

/**
 *  Read and increment an element of a Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_read_inc = pnga_read_inc
#endif

Integer pnga_read_inc(Integer g_a, Integer* subscript, Integer inc)
{
char *ptr;
Integer ldp[MAXDIM], proc, handle=GA_OFFSET+g_a, p_handle, ndim;
int optype,ivalue;
long lvalue;
void *pval;

    ga_check_handleM(g_a, "nga_read_inc");
    GA_PUSH_NAME("nga_read_inc");
    /* BJP printf("p[%d] g_a: %d subscript: %d inc: %d\n",GAme, g_a, subscript[0], inc); */

    if(GA[handle].type!=C_INT && GA[handle].type!=C_LONG &&
       GA[handle].type!=C_LONGLONG)
       pnga_error("type must be integer",GA[handle].type);

    GAstat.numrdi++;
    GAbytes.rditot += (double)sizeof(Integer);
    p_handle = GA[handle].p_handle;
    ndim = GA[handle].ndim;

    /* find out who owns it */
    pnga_locate(g_a, subscript, &proc);

    /* get an address of the g_a(subscript) element */
    if (!GA[handle].block_flag) {
      if (GA[handle].num_rstrctd == 0) {
        gam_Location(proc, handle,  subscript, (char**)&ptr, ldp);
      } else {
        gam_Location(GA[handle].rank_rstrctd[proc], handle,  subscript, (char**)&ptr, ldp);
      }
    } else {
      Integer lo[MAXDIM], hi[MAXDIM];
      Integer j, jtot, last, offset;
      pnga_distribution(g_a, proc, lo, hi);
      pnga_access_block_ptr(g_a, proc, (char**)&ptr, ldp);
      pnga_release_block(g_a, proc);
      offset = 0;
      last = ndim - 1;
      jtot = 1;
      for (j=0; j<last; j++) {
        offset += (subscript[j]-lo[j])*jtot;
        jtot *= ldp[j];
      }
      offset += (subscript[last]-lo[last])*jtot;
      ptr += offset*GA[handle].elemsize;
    }

    if(GA[handle].type==C_INT){
       optype = ARMCI_FETCH_AND_ADD;
       pval = &ivalue;
    }else{
       optype = ARMCI_FETCH_AND_ADD_LONG;
       pval = &lvalue;
    }

    if(GAme == proc)GAbytes.rdiloc += (double)sizeof(Integer);

#ifdef PERMUTE_PIDS
    if(GA_Proc_list) proc = GA_inv_Proc_list[proc];
#endif
    if (GA[handle].block_flag) {
      if (GA[handle].block_sl_flag == 0) {
        proc = proc%pnga_nnodes();
      } else {
        Integer j, index[MAXDIM];
        gam_find_block_indices(handle, proc, index);
        for (j=0; j<ndim; j++) {
          index[j] = index[j]%GA[handle].nblock[j];
        }
        gam_find_proc_from_sl_indices(handle,proc,index);
      }
    }
    if (p_handle != -1) {   
       proc=PGRP_LIST[p_handle].inv_map_proc_list[proc];
       /*printf("\n%d:proc=%d",GAme,proc);fflush(stdout);*/
    }

    /**
     * On BGL, when datatype is C_LONGLONG or if USE_INTEGER8 is defined,
     * (i.e. default Fortran integer is C_LONGLONG), there is no support for
     * 64-bit integers in ARMCI_Rmw. Therefore, offset by 4 bytes as BGL uses
     * Big-endian format, and increment only the 32-bit value.
     * NOTE: On BGL, ARMCI/BGML doesnot support 64-bit integer increment.
     */
#ifdef BGML
    if(GA[handle].type==C_LONGLONG) ptr = ((char*)ptr) + 4 ;
#endif
    
    ARMCI_Rmw(optype, pval, (int*)ptr, (int)inc, (int)proc);

   GA_POP_NAME;

    if(GA[handle].type==C_INT)
         return (Integer) ivalue;
    else
         return (Integer) lvalue;
}

/**
 *  Utility function to correct patch indices for presence of
 *  skips. If no data is left on patch it returns false.
 */
Integer gai_correct_strided_patch(Integer ndim,
                                  Integer *lo,
                                  Integer *skip,
                                  Integer *plo,
                                  Integer *phi)
{
  Integer i, delta;
  for (i=0; i<ndim; i++) {
    delta = plo[i]-lo[i];
    if (delta%skip[i] != 0) {
      plo[i] = plo[i] - delta%skip[i] + skip[i];
    }
    delta = phi[i]-lo[i];
    if (delta%skip[i] != 0) {
      phi[i] = phi[i] - delta%skip[i];
    }
    if (phi[i]<plo[i]) return 0;
  }
  return 1;
}

/**
 *  Utility function to calculate the number of elements in each
 *  stride direction
 */

int gai_ComputeCountWithSkip(Integer ndim, Integer *lo, Integer *hi,
                             Integer *skip, int *count, Integer *nstride)
{
#if 1
  Integer idx;
  int i, istride = 0;
  /*
  if (GAme==0) {
    printf("p[%d] ComputeCount lo[0]: %d hi[0]: %d lo[1]: %d hi[1]: %d\n",GAme,
        lo[0],hi[0],lo[1],hi[1]);
  }
  */
  count[0]=1;
  istride++;
  for (i=0; i<(int)ndim; i++) {
    idx = hi[i] - lo[i];
    if (idx < 0) return 0;
    if (skip[i] > 1) {
      count[istride] = (int)(idx/skip[i]+1);
      istride++;
    } else {
      count[istride] = (int)idx+1;
      istride++;
    }
  }
  *nstride = istride;
  return 1;
#else
  Integer idx;
  int i, istride = 0;
  /*
  if (GAme==0) {
    printf("p[%d] ComputeCount lo[0]: %d hi[0]: %d lo[1]: %d hi[1]: %d\n",GAme,
        lo[0],hi[0],lo[1],hi[1]);
  }
  */
  for (i=0; i<(int)ndim; i++) {
    idx = hi[i] - lo[i];
    if (idx < 0) return 0;
    if (skip[i] > 1) {
      count[istride] = 1;
      istride++;
      count[istride] = (int)(idx/skip[i]+1);
      istride++;
    } else {
      count[istride] = (int)idx+1;
      istride++;
    }
  }
  *nstride = istride;
  return 1;
#endif
}

/**
 *  Utility function to calculate stride lengths on local and
 *  remote processors, taking into account the presence of skips
 */
void gai_SetStrideWithSkip(Integer ndim, Integer size, Integer *ld,
                          Integer *ldrem, int *stride_rem,
                          int *stride_loc, Integer *skip)
{
#if 1
  int i;
  stride_rem[0] = stride_loc[0] = (int)size;
  for (i=0; i<ndim; i++) {
    stride_rem[i+1] = stride_rem[i];
    stride_rem[i] *= skip[i];
    stride_rem[i+1] *= (int)ldrem[i];
    stride_loc[i+1] = stride_loc[i];
    stride_loc[i+1] *= (int)ld[i];
  }
#else
  int i, istride;
  stride_rem[0] = stride_loc[0] = (int)size;
  istride = 0;
  for (i=0; i<ndim; i++) {
    if (skip[i] > 1) {
      stride_rem[istride] *= (int)skip[i];
      stride_rem[istride+1] = stride_rem[istride]/((int)skip[i]);
      stride_loc[istride+1] = stride_loc[istride];
      istride++;
      if (i<ndim-1) {
        stride_rem[istride] *= (int)ldrem[i];
        stride_rem[istride+1] = stride_rem[istride];
        stride_loc[istride] *= (int)ld[i];
        stride_loc[istride+1] = stride_loc[istride];
      }
    } else {
      if (i<ndim-1) {
        stride_rem[istride] *= (int)ldrem[i];
        stride_rem[istride+1] = stride_rem[istride];
        stride_loc[istride] *= ld[i];
        stride_loc[istride+1] = stride_loc[istride];
      }
    }
    istride++;
  }
#endif
}

void gai_ComputePatchIndexWithSkip(Integer ndim, Integer *lo, Integer *plo,
                                   Integer *skip, Integer *ld, Integer *idx_buf)
{
  Integer i, delta, inc, factor;
  delta = plo[0] - lo[0];
  inc = delta%skip[0];
  delta -= inc;
  delta /=  skip[0];
  *idx_buf = delta;
  for (i=0; i<ndim-1; i++) {
    factor = ld[i];
    delta = plo[i+1]-lo[i+1];
    inc = delta%skip[i+1];
    delta -= inc;
    delta /=  skip[i+1];
    *idx_buf += factor*delta;
  }
}

/**
 *  Put an N-dimensional patch of strided data into a Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_strided_put = pnga_strided_put
#endif

void pnga_strided_put(Integer g_a, Integer *lo, Integer *hi, Integer *skip,
                     void    *buf, Integer *ld)
{
  /* g_a:    Global Array handle
     lo[]:   Array of lower indices of patch of global array
     hi[]:   Array of upper indices of patch of global array
     skip[]: Array of skips for each dimension
     buf[]:  Local buffer that patch will be copied from
     ld[]:   ndim-1 physical dimensions of local buffer */
  Integer p, np, handle = GA_OFFSET + g_a;
  Integer idx, size, nstride, p_handle, nproc;
  int i, proc, ndim;
  int use_blocks;

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  use_blocks = GA[handle].block_flag;
  nproc = pnga_nnodes();
  p_handle = GA[handle].p_handle;

  /* check values of skips to make sure they are legitimate */
  for (i = 0; i<ndim; i++) {
    if (skip[i]<1) {
      pnga_error("nga_strided_put: Invalid value of skip along coordinate ",i);
    }
  }

  GA_PUSH_NAME("nga_strided_put");

  if (!use_blocks) {
    /* Locate the processors containing some portion of the patch
       specified by lo and hi and return the results in _ga_map,
       GA_proclist, and np. GA_proclist contains the list of processors
       containing some portion of the patch, _ga_map contains the
       lower and upper indices of the portion of the total patch held by
       a given processor, and np contains the total number of processors
       that contain some portion of the patch. */
    if (!pnga_locate_region(g_a, lo, hi, _ga_map, GA_proclist, &np))
      ga_RegionError(pnga_ndim(g_a), lo, hi, g_a);

    /* Loop over all processors containing a portion of patch */
    gaPermuteProcList(np);
    for (idx=0; idx<np; idx++) {
      Integer ldrem[MAXDIM];
      int stride_rem[2*MAXDIM], stride_loc[2*MAXDIM], count[2*MAXDIM];
      Integer idx_buf, *plo, *phi;
      char *pbuf, *prem;

      p = (Integer)ProcListPerm[(int)idx];
      /* find visible portion of patch held by processor p and return
         the result in plo and phi. Also, get actual processor index
         corresponding to p and store the result in proc. */
      gam_GetRangeFromMap(p, ndim, &plo, &phi);
      proc = (int)GA_proclist[(int)p];

      /* Correct ranges to account for skips in original patch. If no
         data is left in patch jump to next processor in loop. */
      if (!gai_correct_strided_patch((Integer)ndim, lo, skip, plo, phi))
        continue;

      /* get pointer prem in remote buffer to location indexed by plo.
         Also get leading physical dimensions of remote buffer in memory
         in ldrem */
      gam_Location(proc, handle, plo, &prem, ldrem);

      /* get pointer in local buffer to point indexed by plo given that
         the corner of the buffer corresponds to the point indexed by lo */
      gai_ComputePatchIndexWithSkip(ndim, lo, plo, skip, ld, &idx_buf);
      pbuf = size*idx_buf + (char*)buf;

      /* Compute number of elements in each stride region and compute the
         number of stride regions. Store the results in count and nstride */
      if (!gai_ComputeCountWithSkip(ndim, plo, phi, skip, count, &nstride))
        continue;

      /* Scale first element in count by element size. The ARMCI_PutS routine
         uses this convention to figure out memory sizes. */
      count[0] *= size;

      /* Calculate strides in memory for remote processor indexed by proc and
         local buffer */ 
      gai_SetStrideWithSkip(ndim, size, ld, ldrem, stride_rem, stride_loc,
          skip);

      /* BJP */
      if (p_handle != -1) {
        proc = PGRP_LIST[p_handle].inv_map_proc_list[proc];
      }
#ifdef PERMUTE_PIDS
      if (GA_Proc_list) proc = GA_inv_Proc_list[proc];
#endif
      ARMCI_PutS(pbuf, stride_loc, prem, stride_rem, count, nstride-1, proc);
    }
  } else {
    Integer offset, l_offset, last, pinv;
    Integer blo[MAXDIM],bhi[MAXDIM];
    Integer plo[MAXDIM],phi[MAXDIM];
    Integer idx, j, jtot, chk, iproc;
    Integer idx_buf, ldrem[MAXDIM];
    Integer blk_tot = GA[handle].block_total;
    int check1, check2;
    int stride_rem[2*MAXDIM], stride_loc[2*MAXDIM], count[2*MAXDIM];
    char *pbuf, *prem;

    /* GA uses simple block cyclic data distribution */
    if (GA[handle].block_sl_flag == 0) {

      /* loop over all processors */
      for (iproc = 0; iproc < nproc; iproc++) {
        /* loop over all blocks on the processor */
        offset = 0;
        for (idx=iproc; idx < blk_tot; idx += nproc) {
          /* get the block corresponding to the virtual processor idx */
          ga_ownsM(handle, idx, blo, bhi);

          /* check to see if this block overlaps requested block */
          chk = 1;
          for (j=0; j<ndim; j++) {
            /* check to see if at least one end point of the interval
             * represented by blo and bhi falls in the interval
             * represented by lo and hi */
            check1 = ((blo[j] >= lo[j] && blo[j] <= hi[j]) ||
                (bhi[j] >= lo[j] && bhi[j] <= hi[j]));
            /* check to see if interval represented by lo and hi
             * falls entirely within interval represented by blo and bhi */
            check2 = ((lo[j] >= blo[j] && lo[j] <= bhi[j]) &&
                (hi[j] >= blo[j] && hi[j] <= bhi[j]));
            if (!check1 && !check2) {
              chk = 0;
            }
          }
          if (chk) {
            /* get the patch of block that overlaps requested region */
            gam_GetBlockPatch(blo,bhi,lo,hi,plo,phi,ndim);

            /* Correct ranges to account for skips in original patch. If no
               data is left in patch jump to next processor in loop. */
            if (!gai_correct_strided_patch((Integer)ndim, lo, skip, plo, phi))
              continue;

            /* evaluate offset within block */
            last = ndim - 1;
            jtot = 1;
            if (last == 0) ldrem[0] = bhi[0] - blo[0] + 1;
            l_offset = 0;
            for (j=0; j<last; j++) {
              l_offset += (plo[j]-blo[j])*jtot;
              ldrem[j] = bhi[j]-blo[j]+1;
              jtot *= ldrem[j];
            }
            l_offset += (plo[last]-blo[last])*jtot;
            l_offset += offset;

            /* get pointer to data on remote block */
            pinv = idx%nproc;
            if (p_handle > 0) {
              pinv = PGRP_LIST[p_handle].inv_map_proc_list[pinv];
            }
            prem =  GA[handle].ptr[pinv]+l_offset*GA[handle].elemsize;

            gai_ComputePatchIndexWithSkip(ndim, lo, plo, skip, ld, &idx_buf);
            pbuf = size*idx_buf + (char*)buf;

            /* Compute number of elements in each stride region and compute the
               number of stride regions. Store the results in count and nstride */
            if (!gai_ComputeCountWithSkip(ndim, plo, phi, skip, count, &nstride))
              continue;

            /* scale number of rows by element size */
            count[0] *= size;

            /* Calculate strides in memory for remote processor indexed by proc and
               local buffer */
            gai_SetStrideWithSkip(ndim, size, ld, ldrem, stride_rem, stride_loc,
                skip);

            proc = pinv;

            ARMCI_PutS(pbuf,stride_loc,prem,stride_rem,count,nstride-1,proc);

          }
          /* evaluate size of  block idx and use it to increment offset */
          jtot = 1;
          for (j=0; j<ndim; j++) {
            jtot *= bhi[j]-blo[j]+1;
          }
          offset += jtot;
        }
      }
    } else {
      /* GA uses ScaLAPACK block cyclic data distribution */
      Integer proc_index[MAXDIM], index[MAXDIM];
#if COMPACT_SCALAPACK
      Integer itmp;
#else
      Integer /*blk_size[MAXDIM],*/ blk_num[MAXDIM], blk_dim[MAXDIM];
      Integer blk_inc[MAXDIM], blk_jinc;
      Integer blk_ld[MAXDIM],hlf_blk[MAXDIM];
      /*C_Integer *num_blocks;*/
      C_Integer *block_dims;
      int *proc_grid;

      /* Calculate some properties associated with data distribution */
      proc_grid = GA[handle].nblock;
      /*num_blocks = GA[handle].num_blocks;*/
      block_dims = GA[handle].block_dims;
      for (j=0; j<ndim; j++)  {
        blk_dim[j] = block_dims[j]*proc_grid[j];
        blk_num[j] = GA[handle].dims[j]/blk_dim[j];
        /*blk_size[j] = block_dims[j]*blk_num[j];*/
        blk_inc[j] = GA[handle].dims[j]-blk_num[j]*blk_dim[j];
        blk_ld[j] = blk_num[j]*block_dims[j];
        hlf_blk[j] = blk_inc[j]/block_dims[j];
      }
#endif
        /* Loop through all blocks owned by this processor. Decompose
           this loop into a loop over all processors and then a loop
           over all blocks owned by this processor. */
      for (iproc = 0; iproc<GAnproc; iproc++) {
        gam_find_proc_indices(handle, iproc, proc_index);
        gam_find_proc_indices(handle, iproc, index);

        /* Initialize offset for each processor to zero */
        offset = 0;
        while (index[ndim-1] < GA[handle].num_blocks[ndim-1]) {

          /* get bounds for current block */
          for (idx = 0; idx < ndim; idx++) {
            blo[idx] = GA[handle].block_dims[idx]*index[idx]+1;
            bhi[idx] = GA[handle].block_dims[idx]*(index[idx]+1);
            if (bhi[idx] > GA[handle].dims[idx]) bhi[idx] = GA[handle].dims[idx];
          }

          /* check to see if this block overlaps with requested block
           * defined by lo and hi */
          chk = 1;
          for (j=0; j<ndim; j++) {
            /* check to see if at least one end point of the interval
             * represented by blo and bhi falls in the interval
             * represented by lo and hi */
            check1 = ((blo[j] >= lo[j] && blo[j] <= hi[j]) ||
                (bhi[j] >= lo[j] && bhi[j] <= hi[j]));
            /* check to see if interval represented by lo and hi
             * falls entirely within interval represented by blo and bhi */
            check2 = ((lo[j] >= blo[j] && lo[j] <= bhi[j]) &&
                (hi[j] >= blo[j] && hi[j] <= bhi[j]));
            if (!check1 && !check2) {
              chk = 0;
            }
          }
          if (chk) {
            /* get the patch of block that overlaps requested region */
            gam_GetBlockPatch(blo,bhi,lo,hi,plo,phi,ndim);

            /* Correct ranges to account for skips in original patch. If no
               data is left in patch jump to next processor in loop. */
            if (!gai_correct_strided_patch((Integer)ndim, lo, skip, plo, phi))
              continue;

            /* evaluate offset within block */
            last = ndim - 1;
#if COMPACT_SCALAPACK
            jtot = 1;
            if (last == 0) ldrem[0] = bhi[0] - blo[0] + 1;
            l_offset = 0;
            for (j=0; j<last; j++) {
              l_offset += (plo[j]-blo[j])*jtot;
              ldrem[j] = bhi[j]-blo[j]+1;
              jtot *= ldrem[j];
            }
            l_offset += (plo[last]-blo[last])*jtot;
            l_offset += offset;
#else
            l_offset = 0;
            jtot = 1;
            for (j=0; j<last; j++)  {
              ldrem[j] = blk_ld[j];
              blk_jinc = GA[handle].dims[j]%block_dims[j];
              if (blk_inc[j] > 0) {
                if (proc_index[j]<hlf_blk[j]) {
                  blk_jinc = block_dims[j];
                } else if (proc_index[j] == hlf_blk[j]) {
                  blk_jinc = blk_inc[j]%block_dims[j];
                  /*
                  if (blk_jinc == 0) {
                    blk_jinc = block_dims[j];
                  }
                  */
                } else {
                  blk_jinc = 0;
                }
              }
              ldrem[j] += blk_jinc;
              l_offset += (plo[j]-blo[j]
                  + ((blo[j]-1)/blk_dim[j])*block_dims[j])*jtot;
              jtot *= ldrem[j];
            }
            l_offset += (plo[last]-blo[last]
                + ((blo[last]-1)/blk_dim[j])*block_dims[last])*jtot;
#endif

            /* get pointer to data on remote block */
            pinv = iproc;
            if (p_handle > 0) {
              pinv = PGRP_LIST[p_handle].inv_map_proc_list[pinv];
            }
            prem =  GA[handle].ptr[pinv]+l_offset*GA[handle].elemsize;

            gai_ComputePatchIndexWithSkip(ndim, lo, plo, skip, ld, &idx_buf);
            pbuf = size*idx_buf + (char*)buf;        

            /* Compute number of elements in each stride region and compute the
               number of stride regions. Store the results in count and nstride */
            if (!gai_ComputeCountWithSkip(ndim, plo, phi, skip, count, &nstride))
              continue;

            /* scale number of rows by element size */
            count[0] *= size; 

            /* Calculate strides in memory for remote processor indexed by proc and
               local buffer */
            gai_SetStrideWithSkip(ndim, size, ld, ldrem, stride_rem, stride_loc,
                skip);

            proc = pinv;

            ARMCI_PutS(pbuf,stride_loc,prem,stride_rem,count,nstride-1,proc);

          }
          /* increment offset to account for all elements on this block */
#if COMPACT_SCALAPACK
          itmp = 1;
          for (idx = 0; idx < ndim; idx++) {
            itmp *= (bhi[idx] - blo[idx] + 1);
          }
          offset += itmp;
#endif

          /* increment block indices to get the next block on processor iproc */
          index[0] += GA[handle].nblock[0];
          for (idx= 0; idx < ndim; idx++) {
            if (index[idx] >= GA[handle].num_blocks[idx] && idx < ndim-1) {
              index[idx] = proc_index[idx];
              index[idx+1] += GA[handle].nblock[idx+1];
            }
          }
        }
      }
    }
  }
  GA_POP_NAME;
}

/**
 *  Get an N-dimensional patch of strided data from a Global Array
 */
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_strided_get = pnga_strided_get
#endif

void pnga_strided_get(Integer g_a, Integer *lo, Integer *hi, Integer *skip,
                   void    *buf, Integer *ld)
{
  /* g_a:    Global Array handle
     lo[]:   Array of lower indices of patch of global array
     hi[]:   Array of upper indices of patch of global array
     skip[]: Array of skips for each dimension
     buf[]:  Local buffer that patch will be copied from
     ld[]:   ndim-1 physical dimensions of local buffer */
  Integer p, np, handle = GA_OFFSET + g_a;
  Integer idx, size, nstride, p_handle, nproc;
  int i, proc, ndim;
  int use_blocks;

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  use_blocks = GA[handle].block_flag;
  nproc = pnga_nnodes();
  p_handle = GA[handle].p_handle;

  /* check values of skips to make sure they are legitimate */
  for (i = 0; i<ndim; i++) {
    if (skip[i]<1) {
      pnga_error("nga_strided_get: Invalid value of skip along coordinate ",i);
    }
  }

  GA_PUSH_NAME("nga_strided_get");

  if (!use_blocks) {
    /* Locate the processors containing some portion of the patch
       specified by lo and hi and return the results in _ga_map,
       GA_proclist, and np. GA_proclist contains the list of processors
       containing some portion of the patch, _ga_map contains the
       lower and upper indices of the portion of the total patch held by
       a given processor, and np contains the total number of processors
       that contain some portion of the patch. */
    if (!pnga_locate_region(g_a, lo, hi, _ga_map, GA_proclist, &np))
      ga_RegionError(pnga_ndim(g_a), lo, hi, g_a);

    /* Loop over all processors containing a portion of patch */
    gaPermuteProcList(np);
    for (idx=0; idx<np; idx++) {
      Integer ldrem[MAXDIM];
      int stride_rem[2*MAXDIM], stride_loc[2*MAXDIM], count[2*MAXDIM];
      Integer idx_buf, *plo, *phi;
      char *pbuf, *prem;

      p = (Integer)ProcListPerm[(int)idx];
      /* find visible portion of patch held by processor p and return
         the result in plo and phi. Also, get actual processor index
         corresponding to p and store the result in proc. */
      gam_GetRangeFromMap(p, ndim, &plo, &phi);
      proc = (int)GA_proclist[(int)p];

      /* Correct ranges to account for skips in original patch. If no
         data is left in patch jump to next processor in loop. */
      if (!gai_correct_strided_patch((Integer)ndim, lo, skip, plo, phi))
        continue;

      /* get pointer prem in remote buffer to location indexed by plo.
         Also get leading physical dimensions of remote buffer in memory
         in ldrem */
      gam_Location(proc, handle, plo, &prem, ldrem);

      /* get pointer in local buffer to point indexed by plo given that
         the corner of the buffer corresponds to the point indexed by lo */
      gai_ComputePatchIndexWithSkip(ndim, lo, plo, skip, ld, &idx_buf);
      /*
      if (GAme==0) {
        printf("p[%d] Get patch at plo[0]: %d plo[1]: %d idx_buf: %d\n",
               GAme,plo[0],plo[1],idx_buf);
      }
      */
      pbuf = size*idx_buf + (char*)buf;

      /* Compute number of elements in each stride region and compute the
         number of stride regions. Store the results in count and nstride */
      if (!gai_ComputeCountWithSkip(ndim, plo, phi, skip, count, &nstride))
        continue;
      /*
      if (GAme==0) {
        printf("p[%d] count[0]: %d count[1]: %d count[2]: %d count[3]: %d nstride: %d\n",
            GAme,count[0],count[1],count[2],count[3],nstride);
      }
      */

      /* Scale first element in count by element size. The ARMCI_PutS routine
         uses this convention to figure out memory sizes. */
      count[0] *= size;

      /* Calculate strides in memory for remote processor indexed by proc and
         local buffer */ 
      gai_SetStrideWithSkip(ndim, size, ld, ldrem, stride_rem, stride_loc,
          skip);
      /*
      if (GAme==0) {
        printf("p[%d] stride_loc[0]: %d stride_loc[1]: %d stride_loc[2]: %d nstride: %d\n",
            GAme,stride_loc[0],stride_loc[1],stride_loc[2],nstride);
        printf("p[%d] stride_rem[0]: %d stride_rem[1]: %d stride_rem[2]: %d nstride: %d\n",
            GAme,stride_rem[0],stride_rem[1],stride_rem[2],nstride);
      }
      */
      

      /* BJP */
      if (p_handle != -1) {
        proc = PGRP_LIST[p_handle].inv_map_proc_list[proc];
      }
#ifdef PERMUTE_PIDS
      if (GA_Proc_list) proc = GA_inv_Proc_list[proc];
#endif
      ARMCI_GetS(prem, stride_rem, pbuf, stride_loc, count, nstride-1, proc);
    }
  } else {
    Integer offset, l_offset, last, pinv;
    Integer blo[MAXDIM],bhi[MAXDIM];
    Integer plo[MAXDIM],phi[MAXDIM];
    Integer idx, j, jtot, chk, iproc;
    Integer idx_buf, ldrem[MAXDIM];
    Integer blk_tot = GA[handle].block_total;
    int check1, check2;
    int stride_rem[2*MAXDIM], stride_loc[2*MAXDIM], count[2*MAXDIM];
    char *pbuf, *prem;

    /* GA uses simple block cyclic data distribution */
    if (GA[handle].block_sl_flag == 0) {

      /* loop over all processors */
      for (iproc = 0; iproc < nproc; iproc++) {
        /* loop over all blocks on the processor */
        offset = 0;
        for (idx=iproc; idx < blk_tot; idx += nproc) {
          /* get the block corresponding to the virtual processor idx */
          ga_ownsM(handle, idx, blo, bhi);

          /* check to see if this block overlaps requested block */
          chk = 1;
          for (j=0; j<ndim; j++) {
            /* check to see if at least one end point of the interval
             * represented by blo and bhi falls in the interval
             * represented by lo and hi */
            check1 = ((blo[j] >= lo[j] && blo[j] <= hi[j]) ||
                (bhi[j] >= lo[j] && bhi[j] <= hi[j]));
            /* check to see if interval represented by lo and hi
             * falls entirely within interval represented by blo and bhi */
            check2 = ((lo[j] >= blo[j] && lo[j] <= bhi[j]) &&
                (hi[j] >= blo[j] && hi[j] <= bhi[j]));
            if (!check1 && !check2) {
              chk = 0;
            }
          }
          if (chk) {
            /* get the patch of block that overlaps requested region */
            gam_GetBlockPatch(blo,bhi,lo,hi,plo,phi,ndim);

            /* Correct ranges to account for skips in original patch. If no
               data is left in patch jump to next processor in loop. */
            if (!gai_correct_strided_patch((Integer)ndim, lo, skip, plo, phi))
              continue;

            /* evaluate offset within block */
            last = ndim - 1;
            jtot = 1;
            if (last == 0) ldrem[0] = bhi[0] - blo[0] + 1;
            l_offset = 0;
            for (j=0; j<last; j++) {
              l_offset += (plo[j]-blo[j])*jtot;
              ldrem[j] = bhi[j]-blo[j]+1;
              jtot *= ldrem[j];
            }
            l_offset += (plo[last]-blo[last])*jtot;
            l_offset += offset;

            /* get pointer to data on remote block */
            pinv = idx%nproc;
            if (p_handle > 0) {
              pinv = PGRP_LIST[p_handle].inv_map_proc_list[pinv];
            }
            prem =  GA[handle].ptr[pinv]+l_offset*GA[handle].elemsize;

            gai_ComputePatchIndexWithSkip(ndim, lo, plo, skip, ld, &idx_buf);
            pbuf = size*idx_buf + (char*)buf;

            /* Compute number of elements in each stride region and compute the
               number of stride regions. Store the results in count and nstride */
            if (!gai_ComputeCountWithSkip(ndim, plo, phi, skip, count, &nstride))
              continue;

            /* scale number of rows by element size */
            count[0] *= size;

            /* Calculate strides in memory for remote processor indexed by proc and
               local buffer */
            gai_SetStrideWithSkip(ndim, size, ld, ldrem, stride_rem, stride_loc,
                skip);

            proc = pinv;

            ARMCI_GetS(prem,stride_rem,pbuf,stride_loc,count,nstride-1,proc);

          }
          /* evaluate size of  block idx and use it to increment offset */
          jtot = 1;
          for (j=0; j<ndim; j++) {
            jtot *= bhi[j]-blo[j]+1;
          }
          offset += jtot;
        }
      }
    } else {
      /* GA uses ScaLAPACK block cyclic data distribution */
      Integer proc_index[MAXDIM], index[MAXDIM];
#if COMPACT_SCALAPACK
      Integer itmp;
#else
      Integer /*blk_size[MAXDIM],*/ blk_num[MAXDIM], blk_dim[MAXDIM];
      Integer blk_inc[MAXDIM], blk_jinc;
      Integer blk_ld[MAXDIM],hlf_blk[MAXDIM];
      /*C_Integer *num_blocks;*/
      C_Integer *block_dims;
      int *proc_grid;

      /* Calculate some properties associated with data distribution */
      proc_grid = GA[handle].nblock;
      /*num_blocks = GA[handle].num_blocks;*/
      block_dims = GA[handle].block_dims;
      for (j=0; j<ndim; j++)  {
        blk_dim[j] = block_dims[j]*proc_grid[j];
        blk_num[j] = GA[handle].dims[j]/blk_dim[j];
        /*blk_size[j] = block_dims[j]*blk_num[j];*/
        blk_inc[j] = GA[handle].dims[j]-blk_num[j]*blk_dim[j];
        blk_ld[j] = blk_num[j]*block_dims[j];
        hlf_blk[j] = blk_inc[j]/block_dims[j];
      }
#endif
        /* Loop through all blocks owned by this processor. Decompose
           this loop into a loop over all processors and then a loop
           over all blocks owned by this processor. */
      for (iproc = 0; iproc<GAnproc; iproc++) {
        gam_find_proc_indices(handle, iproc, proc_index);
        gam_find_proc_indices(handle, iproc, index);

        /* Initialize offset for each processor to zero */
        offset = 0;
        while (index[ndim-1] < GA[handle].num_blocks[ndim-1]) {

          /* get bounds for current block */
          for (idx = 0; idx < ndim; idx++) {
            blo[idx] = GA[handle].block_dims[idx]*index[idx]+1;
            bhi[idx] = GA[handle].block_dims[idx]*(index[idx]+1);
            if (bhi[idx] > GA[handle].dims[idx]) bhi[idx] = GA[handle].dims[idx];
          }

          /* check to see if this block overlaps with requested block
           * defined by lo and hi */
          chk = 1;
          for (j=0; j<ndim; j++) {
            /* check to see if at least one end point of the interval
             * represented by blo and bhi falls in the interval
             * represented by lo and hi */
            check1 = ((blo[j] >= lo[j] && blo[j] <= hi[j]) ||
                (bhi[j] >= lo[j] && bhi[j] <= hi[j]));
            /* check to see if interval represented by lo and hi
             * falls entirely within interval represented by blo and bhi */
            check2 = ((lo[j] >= blo[j] && lo[j] <= bhi[j]) &&
                (hi[j] >= blo[j] && hi[j] <= bhi[j]));
            if (!check1 && !check2) {
              chk = 0;
            }
          }
          if (chk) {
            /* get the patch of block that overlaps requested region */
            gam_GetBlockPatch(blo,bhi,lo,hi,plo,phi,ndim);

            /* Correct ranges to account for skips in original patch. If no
               data is left in patch jump to next processor in loop. */
            if (!gai_correct_strided_patch((Integer)ndim, lo, skip, plo, phi))
              continue;

            /* evaluate offset within block */
            last = ndim - 1;
#if COMPACT_SCALAPACK
            jtot = 1;
            if (last == 0) ldrem[0] = bhi[0] - blo[0] + 1;
            l_offset = 0;
            for (j=0; j<last; j++) {
              l_offset += (plo[j]-blo[j])*jtot;
              ldrem[j] = bhi[j]-blo[j]+1;
              jtot *= ldrem[j];
            }
            l_offset += (plo[last]-blo[last])*jtot;
            l_offset += offset;
#else
            l_offset = 0;
            jtot = 1;
            for (j=0; j<last; j++)  {
              ldrem[j] = blk_ld[j];
              blk_jinc = GA[handle].dims[j]%block_dims[j];
              if (blk_inc[j] > 0) {
                if (proc_index[j]<hlf_blk[j]) {
                  blk_jinc = block_dims[j];
                } else if (proc_index[j] == hlf_blk[j]) {
                  blk_jinc = blk_inc[j]%block_dims[j];
                  /*
                  if (blk_jinc == 0) {
                    blk_jinc = block_dims[j];
                  }
                  */
                } else {
                  blk_jinc = 0;
                }
              }
              ldrem[j] += blk_jinc;
              l_offset += (plo[j]-blo[j]
                  + ((blo[j]-1)/blk_dim[j])*block_dims[j])*jtot;
              jtot *= ldrem[j];
            }
            l_offset += (plo[last]-blo[last]
                + ((blo[last]-1)/blk_dim[j])*block_dims[last])*jtot;
#endif

            /* get pointer to data on remote block */
            pinv = iproc;
            if (p_handle > 0) {
              pinv = PGRP_LIST[p_handle].inv_map_proc_list[pinv];
            }
            prem =  GA[handle].ptr[pinv]+l_offset*GA[handle].elemsize;

            gai_ComputePatchIndexWithSkip(ndim, lo, plo, skip, ld, &idx_buf);
            pbuf = size*idx_buf + (char*)buf;        

            /* Compute number of elements in each stride region and compute the
               number of stride regions. Store the results in count and nstride */
            if (!gai_ComputeCountWithSkip(ndim, plo, phi, skip, count, &nstride))
              continue;

            /* scale number of rows by element size */
            count[0] *= size; 

            /* Calculate strides in memory for remote processor indexed by proc and
               local buffer */
            gai_SetStrideWithSkip(ndim, size, ld, ldrem, stride_rem, stride_loc,
                skip);

            proc = pinv;

            ARMCI_GetS(prem,stride_rem,pbuf,stride_loc,count,nstride-1,proc);

          }
          /* increment offset to account for all elements on this block */
#if COMPACT_SCALAPACK
          itmp = 1;
          for (idx = 0; idx < ndim; idx++) {
            itmp *= (bhi[idx] - blo[idx] + 1);
          }
          offset += itmp;
#endif

          /* increment block indices to get the next block on processor iproc */
          index[0] += GA[handle].nblock[0];
          for (idx= 0; idx < ndim; idx++) {
            if (index[idx] >= GA[handle].num_blocks[idx] && idx < ndim-1) {
              index[idx] = proc_index[idx];
              index[idx+1] += GA[handle].nblock[idx+1];
            }
          }
        }
      }
    }
  }
  GA_POP_NAME;
}

/**
 *  Accumulate operation for an N-dimensional patch of strided data of
 *  a Global Array
 *       g_a += alpha * patch
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_strided_acc = pnga_strided_acc
#endif

void pnga_strided_acc(Integer g_a, Integer *lo, Integer *hi, Integer *skip,
                      void *buf, Integer *ld, void *alpha)
{
  /* g_a:    Global Array handle
     lo[]:   Array of lower indices of patch of global array
     hi[]:   Array of upper indices of patch of global array
     skip[]: Array of skips for each dimension
     buf[]:  Local buffer that patch will be copied from
     ld[]:   ndim-1 physical dimensions of local buffer
     alpha:  muliplicative scale factor */
  Integer p, np, handle = GA_OFFSET + g_a;
  Integer idx, size, nstride, type, p_handle, nproc;
  int i, optype=-1, proc, ndim;
  int use_blocks;

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  type = GA[handle].type;
  use_blocks = GA[handle].block_flag;
  nproc = pnga_nnodes();
  p_handle = GA[handle].p_handle;

  if (type == C_DBL) optype = ARMCI_ACC_DBL;
  else if (type == C_FLOAT) optype = ARMCI_ACC_FLT;
  else if (type == C_DCPL) optype = ARMCI_ACC_DCP;
  else if (type == C_SCPL) optype = ARMCI_ACC_CPL;
  else if (type == C_INT) optype = ARMCI_ACC_INT;
  else if (type == C_LONG) optype = ARMCI_ACC_LNG;
  else pnga_error("nga_strided_acc: type not supported",type);

  /* check values of skips to make sure they are legitimate */
  for (i = 0; i<ndim; i++) {
    if (skip[i]<1) {
      pnga_error("nga_strided_acc: Invalid value of skip along coordinate ",i);
    }
  }

  GA_PUSH_NAME("nga_strided_acc");

  if (!use_blocks) {
    /* Locate the processors containing some portion of the patch
       specified by lo and hi and return the results in _ga_map,
       GA_proclist, and np. GA_proclist contains the list of processors
       containing some portion of the patch, _ga_map contains the
       lower and upper indices of the portion of the total patch held by
       a given processor, and np contains the total number of processors
       that contain some portion of the patch. */
    if (!pnga_locate_region(g_a, lo, hi, _ga_map, GA_proclist, &np))
      ga_RegionError(pnga_ndim(g_a), lo, hi, g_a);

    /* Loop over all processors containing a portion of patch */
    gaPermuteProcList(np);
    for (idx=0; idx<np; idx++) {
      Integer ldrem[MAXDIM];
      int stride_rem[2*MAXDIM], stride_loc[2*MAXDIM], count[2*MAXDIM];
      Integer idx_buf, *plo, *phi;
      char *pbuf, *prem;

      p = (Integer)ProcListPerm[(int)idx];
      /* find visible portion of patch held by processor p and return
         the result in plo and phi. Also, get actual processor index
         corresponding to p and store the result in proc. */
      gam_GetRangeFromMap(p, ndim, &plo, &phi);
      proc = (int)GA_proclist[(int)p];

      /* Correct ranges to account for skips in original patch. If no
         data is left in patch jump to next processor in loop. */
      if (!gai_correct_strided_patch((Integer)ndim, lo, skip, plo, phi))
        continue;

      /* get pointer prem in remote buffer to location indexed by plo.
         Also get leading physical dimensions of remote buffer in memory
         in ldrem */
      gam_Location(proc, handle, plo, &prem, ldrem);

      /* get pointer in local buffer to point indexed by plo given that
         the corner of the buffer corresponds to the point indexed by lo */
      gai_ComputePatchIndexWithSkip(ndim, lo, plo, skip, ld, &idx_buf);
      pbuf = size*idx_buf + (char*)buf;

      /* Compute number of elements in each stride region and compute the
         number of stride regions. Store the results in count and nstride */
      if (!gai_ComputeCountWithSkip(ndim, plo, phi, skip, count, &nstride))
        continue;

      /* Scale first element in count by element size. The ARMCI_AccS routine
         uses this convention to figure out memory sizes. */
      count[0] *= size;

      /* Calculate strides in memory for remote processor indexed by proc and
         local buffer */ 
      gai_SetStrideWithSkip(ndim, size, ld, ldrem, stride_rem, stride_loc,
          skip);

      /* BJP */
      if (p_handle != -1) {
        proc = PGRP_LIST[p_handle].inv_map_proc_list[proc];
      }
#ifdef PERMUTE_PIDS
      if (GA_Proc_list) proc = GA_inv_Proc_list[proc];
#endif
      ARMCI_AccS(optype, alpha, pbuf, stride_loc, prem, stride_rem, count,
          nstride-1, proc);
    }
  } else {
    Integer offset, l_offset, last, pinv;
    Integer blo[MAXDIM],bhi[MAXDIM];
    Integer plo[MAXDIM],phi[MAXDIM];
    Integer idx, j, jtot, chk, iproc;
    Integer idx_buf, ldrem[MAXDIM];
    Integer blk_tot = GA[handle].block_total;
    int check1, check2;
    int stride_rem[2*MAXDIM], stride_loc[2*MAXDIM], count[2*MAXDIM];
    char *pbuf, *prem;

    /* GA uses simple block cyclic data distribution */
    if (GA[handle].block_sl_flag == 0) {

      /* loop over all processors */
      for (iproc = 0; iproc < nproc; iproc++) {
        /* loop over all blocks on the processor */
        offset = 0;
        for (idx=iproc; idx < blk_tot; idx += nproc) {
          /* get the block corresponding to the virtual processor idx */
          ga_ownsM(handle, idx, blo, bhi);

          /* check to see if this block overlaps requested block */
          chk = 1;
          for (j=0; j<ndim; j++) {
            /* check to see if at least one end point of the interval
             * represented by blo and bhi falls in the interval
             * represented by lo and hi */
            check1 = ((blo[j] >= lo[j] && blo[j] <= hi[j]) ||
                (bhi[j] >= lo[j] && bhi[j] <= hi[j]));
            /* check to see if interval represented by lo and hi
             * falls entirely within interval represented by blo and bhi */
            check2 = ((lo[j] >= blo[j] && lo[j] <= bhi[j]) &&
                (hi[j] >= blo[j] && hi[j] <= bhi[j]));
            if (!check1 && !check2) {
              chk = 0;
            }
          }
          if (chk) {
            /* get the patch of block that overlaps requested region */
            gam_GetBlockPatch(blo,bhi,lo,hi,plo,phi,ndim);

            /* Correct ranges to account for skips in original patch. If no
               data is left in patch jump to next processor in loop. */
            if (!gai_correct_strided_patch((Integer)ndim, lo, skip, plo, phi))
              continue;

            /* evaluate offset within block */
            last = ndim - 1;
            jtot = 1;
            if (last == 0) ldrem[0] = bhi[0] - blo[0] + 1;
            l_offset = 0;
            for (j=0; j<last; j++) {
              l_offset += (plo[j]-blo[j])*jtot;
              ldrem[j] = bhi[j]-blo[j]+1;
              jtot *= ldrem[j];
            }
            l_offset += (plo[last]-blo[last])*jtot;
            l_offset += offset;

            /* get pointer to data on remote block */
            pinv = idx%nproc;
            if (p_handle > 0) {
              pinv = PGRP_LIST[p_handle].inv_map_proc_list[pinv];
            }
            prem =  GA[handle].ptr[pinv]+l_offset*GA[handle].elemsize;

            gai_ComputePatchIndexWithSkip(ndim, lo, plo, skip, ld, &idx_buf);
            pbuf = size*idx_buf + (char*)buf;

            /* Compute number of elements in each stride region and compute the
               number of stride regions. Store the results in count and nstride */
            if (!gai_ComputeCountWithSkip(ndim, plo, phi, skip, count, &nstride))
              continue;

            /* scale number of rows by element size */
            count[0] *= size;

            /* Calculate strides in memory for remote processor indexed by proc and
               local buffer */
            gai_SetStrideWithSkip(ndim, size, ld, ldrem, stride_rem, stride_loc,
                skip);

            proc = pinv;

            ARMCI_AccS(optype, alpha, pbuf, stride_loc, prem, stride_rem, count,
                       nstride-1, proc);

          }
          /* evaluate size of  block idx and use it to increment offset */
          jtot = 1;
          for (j=0; j<ndim; j++) {
            jtot *= bhi[j]-blo[j]+1;
          }
          offset += jtot;
        }
      }
    } else {
      /* GA uses ScaLAPACK block cyclic data distribution */
      Integer proc_index[MAXDIM], index[MAXDIM];
#if COMPACT_SCALAPACK
      Integer itmp;
#else
      Integer /*blk_size[MAXDIM],*/ blk_num[MAXDIM], blk_dim[MAXDIM];
      Integer blk_inc[MAXDIM], blk_jinc;
      Integer blk_ld[MAXDIM],hlf_blk[MAXDIM];
      /*C_Integer *num_blocks;*/
      C_Integer *block_dims;
      int *proc_grid;

      /* Calculate some properties associated with data distribution */
      proc_grid = GA[handle].nblock;
      /*num_blocks = GA[handle].num_blocks;*/
      block_dims = GA[handle].block_dims;
      for (j=0; j<ndim; j++)  {
        blk_dim[j] = block_dims[j]*proc_grid[j];
        blk_num[j] = GA[handle].dims[j]/blk_dim[j];
        /*blk_size[j] = block_dims[j]*blk_num[j];*/
        blk_inc[j] = GA[handle].dims[j]-blk_num[j]*blk_dim[j];
        blk_ld[j] = blk_num[j]*block_dims[j];
        hlf_blk[j] = blk_inc[j]/block_dims[j];
      }
#endif
        /* Loop through all blocks owned by this processor. Decompose
           this loop into a loop over all processors and then a loop
           over all blocks owned by this processor. */
      for (iproc = 0; iproc<GAnproc; iproc++) {
        gam_find_proc_indices(handle, iproc, proc_index);
        gam_find_proc_indices(handle, iproc, index);

        /* Initialize offset for each processor to zero */
        offset = 0;
        while (index[ndim-1] < GA[handle].num_blocks[ndim-1]) {

          /* get bounds for current block */
          for (idx = 0; idx < ndim; idx++) {
            blo[idx] = GA[handle].block_dims[idx]*index[idx]+1;
            bhi[idx] = GA[handle].block_dims[idx]*(index[idx]+1);
            if (bhi[idx] > GA[handle].dims[idx]) bhi[idx] = GA[handle].dims[idx];
          }

          /* check to see if this block overlaps with requested block
           * defined by lo and hi */
          chk = 1;
          for (j=0; j<ndim; j++) {
            /* check to see if at least one end point of the interval
             * represented by blo and bhi falls in the interval
             * represented by lo and hi */
            check1 = ((blo[j] >= lo[j] && blo[j] <= hi[j]) ||
                (bhi[j] >= lo[j] && bhi[j] <= hi[j]));
            /* check to see if interval represented by lo and hi
             * falls entirely within interval represented by blo and bhi */
            check2 = ((lo[j] >= blo[j] && lo[j] <= bhi[j]) &&
                (hi[j] >= blo[j] && hi[j] <= bhi[j]));
            if (!check1 && !check2) {
              chk = 0;
            }
          }
          if (chk) {
            /* get the patch of block that overlaps requested region */
            gam_GetBlockPatch(blo,bhi,lo,hi,plo,phi,ndim);

            /* Correct ranges to account for skips in original patch. If no
               data is left in patch jump to next processor in loop. */
            if (!gai_correct_strided_patch((Integer)ndim, lo, skip, plo, phi))
              continue;

            /* evaluate offset within block */
            last = ndim - 1;
#if COMPACT_SCALAPACK
            jtot = 1;
            if (last == 0) ldrem[0] = bhi[0] - blo[0] + 1;
            l_offset = 0;
            for (j=0; j<last; j++) {
              l_offset += (plo[j]-blo[j])*jtot;
              ldrem[j] = bhi[j]-blo[j]+1;
              jtot *= ldrem[j];
            }
            l_offset += (plo[last]-blo[last])*jtot;
            l_offset += offset;
#else
            l_offset = 0;
            jtot = 1;
            for (j=0; j<last; j++)  {
              ldrem[j] = blk_ld[j];
              blk_jinc = GA[handle].dims[j]%block_dims[j];
              if (blk_inc[j] > 0) {
                if (proc_index[j]<hlf_blk[j]) {
                  blk_jinc = block_dims[j];
                } else if (proc_index[j] == hlf_blk[j]) {
                  blk_jinc = blk_inc[j]%block_dims[j];
                  /*
                  if (blk_jinc == 0) {
                    blk_jinc = block_dims[j];
                  }
                  */
                } else {
                  blk_jinc = 0;
                }
              }
              ldrem[j] += blk_jinc;
              l_offset += (plo[j]-blo[j]
                  + ((blo[j]-1)/blk_dim[j])*block_dims[j])*jtot;
              jtot *= ldrem[j];
            }
            l_offset += (plo[last]-blo[last]
                + ((blo[last]-1)/blk_dim[j])*block_dims[last])*jtot;
#endif

            /* get pointer to data on remote block */
            pinv = iproc;
            if (p_handle > 0) {
              pinv = PGRP_LIST[p_handle].inv_map_proc_list[pinv];
            }
            prem =  GA[handle].ptr[pinv]+l_offset*GA[handle].elemsize;

            gai_ComputePatchIndexWithSkip(ndim, lo, plo, skip, ld, &idx_buf);
            pbuf = size*idx_buf + (char*)buf;        

            /* Compute number of elements in each stride region and compute the
               number of stride regions. Store the results in count and nstride */
            if (!gai_ComputeCountWithSkip(ndim, plo, phi, skip, count, &nstride))
              continue;

            /* scale number of rows by element size */
            count[0] *= size; 

            /* Calculate strides in memory for remote processor indexed by proc and
               local buffer */
            gai_SetStrideWithSkip(ndim, size, ld, ldrem, stride_rem, stride_loc,
                skip);

            proc = pinv;

            ARMCI_AccS(optype, alpha, pbuf, stride_loc, prem, stride_rem, count,
                       nstride-1, proc);

          }
          /* increment offset to account for all elements on this block */
#if COMPACT_SCALAPACK
          itmp = 1;
          for (idx = 0; idx < ndim; idx++) {
            itmp *= (bhi[idx] - blo[idx] + 1);
          }
          offset += itmp;
#endif

          /* increment block indices to get the next block on processor iproc */
          index[0] += GA[handle].nblock[0];
          for (idx= 0; idx < ndim; idx++) {
            if (index[idx] >= GA[handle].num_blocks[idx] && idx < ndim-1) {
              index[idx] = proc_index[idx];
              index[idx+1] += GA[handle].nblock[idx+1];
            }
          }
        }
      }
    }
  }
  GA_POP_NAME;
}
