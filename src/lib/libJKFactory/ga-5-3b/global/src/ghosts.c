#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: ghosts.c,v 1.47.4.2 2007-05-02 16:23:39 d3g293 Exp $ */
/* 
 * module: ghosts.c
 * author: Bruce Palmer
 * description: implements GA collective communication operations to
 * update ghost cell regions.
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

 
/*#define PERMUTE_PIDS */

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
#include "global.h"
#include "globalp.h"
#include "base.h"
#include "armci.h"
#include "message.h"
#include "macdecls.h"
#include "ga-papi.h"
#include "ga-wapi.h"

/* from armcip.h, but armcip.h is private so we should not include it */
extern void armci_write_strided(void *ptr, int stride_levels, int stride_arr[], int count[], char *buf);
extern void armci_read_strided(void *ptr, int stride_levels, int stride_arr[], int count[], char *buf);
extern armci_hdl_t* get_armci_nbhandle(Integer *);

#define USE_MALLOC 1
#define INVALID_MA_HANDLE -1 
#define NEAR_INT(x) (x)< 0.0 ? ceil( (x) - 0.5) : floor((x) + 0.5)

#if !defined(CRAY_YMP)
#define BYTE_ADDRESSABLE_MEMORY
#endif

/*uncomment line below to verify consistency of MA in every sync */
/*#define CHECK_MA yes */

/***************************************************************************/

/*\ Return a pointer to the location indicated by subscript and and an array
 * of leading dimensions (ld). Assume that subscript refers to a set of local
 * coordinates relative to the origin of the array and account for the
 * presence of ghost cells.
\*/
#define gam_LocationWithGhosts(proc, handle, subscript, ptr_loc, ld)           \
{                                                                              \
Integer _d, _factor = 1, _last=GA[handle].ndim - 1, _offset=0;                 \
Integer _lo[MAXDIM], _hi[MAXDIM];                                              \
  ga_ownsM(handle, proc, _lo, _hi);                                            \
  if (_last == 0) ld[0] = _hi[0] - _lo[0] + 1 + 2*(Integer)GA[handle].width[0];\
  for (_d = 0; _d < _last; _d++) {                                             \
    _offset += subscript[_d] * _factor;                                        \
    ld[_d] = _hi[_d] - _lo[_d] + 1 + 2*(Integer)GA[handle].width[_d];          \
    _factor *= ld[_d];                                                         \
  }                                                                            \
  _offset += subscript[_last] * _factor;                                       \
  *(ptr_loc) = GA[handle].ptr[proc] + _offset*GA[handle].elemsize;             \
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_access_ghost_ptr = pnga_access_ghost_ptr
#endif
void pnga_access_ghost_ptr(Integer g_a, Integer dims[],
                      void* ptr, Integer ld[])

{
char *lptr;
Integer  handle = GA_OFFSET + g_a;
Integer  i, lo[MAXDIM], hi[MAXDIM];
Integer ndim = GA[handle].ndim;
Integer me = pnga_nodeid();

   GA_PUSH_NAME("pnga_access_ghost_ptr");

   pnga_distribution(g_a, me, lo, hi);

   for (i=0; i < ndim; i++) {
     dims[i] = 0;
   }

   gam_LocationWithGhosts(me, handle, dims, &lptr, ld);
   *(char**)ptr = lptr; 
   for (i=0; i < ndim; i++)
     dims[i] = hi[i] - lo[i] + 1 + 2*(Integer)GA[handle].width[i];
   GA_POP_NAME;
}

/*\  PROVIDE INDEX TO LOCALLY HELD DATA, ACCOUNTING FOR
 *   PRESENCE OF GHOST CELLS
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_access_ghost_element = pnga_access_ghost_element
#endif
void pnga_access_ghost_element(Integer g_a, AccessIndex* index,
                        Integer subscript[], Integer ld[])
{
char *ptr=NULL;
Integer  handle = GA_OFFSET + g_a;
Integer i=0;
Integer tmp_sub[MAXDIM];
unsigned long    elemsize=0;
unsigned long    lref=0, lptr=0;
Integer me = pnga_nodeid();
   GA_PUSH_NAME("nga_access_ghost_element");
   /* Indices conform to Fortran convention. Shift them down 1 so that
      gam_LocationWithGhosts works. */
   for (i=0; i<GA[handle].ndim; i++) tmp_sub[i] = subscript[i] - 1;
   gam_LocationWithGhosts(me, handle, tmp_sub, &ptr, ld);
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

/*\  PROVIDE POINTER TO LOCALLY HELD DATA, ACCOUNTING FOR 
 *   PRESENCE OF GHOST CELLS 
\*/ 
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_access_ghost_element_ptr = pnga_access_ghost_element_ptr
#endif
void pnga_access_ghost_element_ptr(Integer g_a, void *ptr, 
                        Integer subscript[], Integer ld[]) 
{ 
  char *lptr; 
  Integer  handle = GA_OFFSET + g_a; 
  Integer i; 
  Integer tmp_sub[MAXDIM]; 
  Integer me = pnga_nodeid(); 
  GA_PUSH_NAME("nga_access_ghost_element_ptr"); 
  /* Indices conform to Fortran convention. Shift them down 1 so that 
     gam_LocationWithGhosts works. */ 
  for (i=0; i<GA[handle].ndim; i++) tmp_sub[i] = subscript[i] - 1; 
  gam_LocationWithGhosts(me, handle, tmp_sub, &lptr, ld); 
 
  *(char**)ptr = lptr; 
  GA_POP_NAME; 
} 
 
/*\ PROVIDE ACCESS TO LOCAL PATCH OF A GLOBAL ARRAY WITH GHOST CELLS
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_access_ghosts = pnga_access_ghosts
#endif
void pnga_access_ghosts(Integer g_a, Integer dims[],
                      AccessIndex* index, Integer ld[])
{
char     *ptr=NULL;
Integer  handle = GA_OFFSET + g_a;
unsigned long    elemsize=0;
unsigned long    lref=0, lptr=0;

   GA_PUSH_NAME("nga_access_ghosts");
   pnga_access_ghost_ptr(g_a, dims, &ptr, ld);

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

/*\ RELEASE ACCESS TO A GHOST ELEMENT
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_release_ghost_element = pnga_release_ghost_element
#endif
void pnga_release_ghost_element(Integer g_a, Integer subscript[])
{
}

/*\ RELEASE ACCESS & UPDATE A GHOST ELEMENT
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_release_update_ghost_element = pnga_release_update_ghost_element
#endif
void pnga_release_update_ghost_element(Integer g_a, Integer subscript[])
{
}

/*\ RELEASE ACCESS TO A GHOST BLOCK
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_release_ghosts = pnga_release_ghosts
#endif
void pnga_release_ghosts(Integer g_a)
{
}

/*\ RELEASE ACCESS & UPDATE A GHOST BLOCK
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_release_update_ghosts = pnga_release_update_ghosts
#endif
void pnga_release_update_ghosts(Integer g_a)
{
}

/*\ GET DATA FROM LOCAL BLOCK
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_get_ghost_block = pnga_get_ghost_block
#endif
void pnga_get_ghost_block(Integer g_a,
                               Integer *lo,
                               Integer *hi,
                               void *buf,
                               Integer *ld)
{
  /* g_a:      Global array handle
     lo[]:     Array of lower indices of patch of global array
     hi[]:     Array of upper indices of patch of global array
     buf[]:    Local buffer that array patch will be copied into
     ld[]:     Array of physical ndim-1 dimensions of local buffer */
  Integer handle=GA_OFFSET + g_a, ndim;
  Integer i, glo[MAXDIM], ghi[MAXDIM], ichk, me, grp_id;
  Integer llo[MAXDIM];
  int  stride_rem[MAXDIM], stride_loc[MAXDIM], count[MAXDIM];
  Integer ldrem[MAXDIM];
  Integer offset, factor, size;
  char *ptr;

  me = GAme;
  grp_id = (Integer)GA[handle].p_handle;
  if (grp_id>0) me = PGRP_LIST[grp_id].map_proc_list[me];
  ndim = GA[handle].ndim;

  /* Figure out whether or not lo and hi can be accessed completely
     from local data */
  pnga_distribution(g_a, me, glo, ghi);
  ichk = 1;
  for (i=0; i<ndim; i++) {
    if (lo[i] < glo[i]-(Integer)GA[handle].width[i]) ichk = 0;
    if (hi[i] > ghi[i]+(Integer)GA[handle].width[i]) ichk = 0;
    llo[i] = glo[i] - (Integer)GA[handle].width[i];
    if (i<ndim-1) ldrem[i] = ghi[i] - glo[i] + 1
      + 2*(Integer)GA[handle].width[i];
  }

  /* Get data. Use local copy if possible, otherwise use a periodic get */
  if (ichk) {
    offset = 0;
    factor = 1;
    size = GA[handle].elemsize;
    for (i=0; i<ndim-1; i++) {
      offset += (lo[i]-llo[i])*factor;
      factor *= ghi[i] - glo[i] + 1 + 2*(Integer)GA[handle].width[i];
    }
    offset += (lo[ndim-1]-llo[ndim-1])*factor;
    ptr = GA[handle].ptr[me] + size*offset;
    /* compute number of elements in each dimension and store result in count */
    gam_ComputeCount(ndim, lo, hi, count);

    /* scale first element in count by element size. The ARMCI_GetS
       routine uses this convention to figure out memory sizes.*/
    count[0] *= size;

    /* Return strides for memory containing global array on remote
       processor indexed by proc (stride_rem) and for local buffer
       buf (stride_loc) */
    gam_setstride(ndim, size, ld, ldrem, stride_rem, stride_loc);
    ARMCI_GetS(ptr,stride_rem,buf,stride_loc,count,ndim-1,me);
  } else {
    pnga_periodic(g_a,lo,hi,buf,ld,NULL,PERIODIC_GET);
  }
}

/*\ UPDATE GHOST CELLS OF GLOBAL ARRAY USING SHIFT ALGORITHM
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_update1_ghosts = pnga_update1_ghosts
#endif
void pnga_update1_ghosts(Integer g_a)
{
  Integer idx, ipx, inx, i, np, handle=GA_OFFSET + g_a, proc_rem;
  Integer size, ndim, nwidth, offset, slice, increment[MAXDIM];
  Integer width[MAXDIM];
  Integer dims[MAXDIM], imax=0;
  Integer lo_loc[MAXDIM], hi_loc[MAXDIM];
  Integer plo_loc[MAXDIM]/*, phi_loc[MAXDIM]*/;
  Integer lo_rem[MAXDIM], hi_rem[MAXDIM];
  Integer tlo_rem[MAXDIM], thi_rem[MAXDIM];
  Integer slo_rem[MAXDIM], shi_rem[MAXDIM];
  Integer plo_rem[MAXDIM], phi_rem[MAXDIM];
  Integer ld_loc[MAXDIM], ld_rem[MAXDIM];
  int corner_flag;
  int stride_loc[MAXDIM], stride_rem[MAXDIM], count[MAXDIM];
  char *ptr_loc, *ptr_rem;
  logical hasData = TRUE;
  Integer me = pnga_nodeid();
  Integer p_handle;

  /* This routine makes use of the shift algorithm to update data in the
   * ghost cells bounding the local block of visible data. The shift
   * algorithm starts by updating the blocks of data along the first
   * dimension by grabbing a block of data that is width[0] deep but
   * otherwise matches the  dimensions of the data residing on the
   * calling processor. The update of the second dimension, however,
   * grabs a block that is width[1] deep in the second dimension but is
   * ldim0 + 2*width[0] in the first dimensions where ldim0 is the
   * size of the visible data along the first dimension. The remaining
   * dimensions are left the same. For the next update, the width of the
   * second dimension is also increased by 2*width[1] and so on. This
   * algorith makes use of the fact that data for the dimensions that
   * have already been updated is available on each processor and can be
   * used in the updates of subsequent dimensions. The total number of
   * separate updates is 2*ndim, an update in the negative and positive
   * directions for each dimension.
   *
   * To perform the update, this routine makes use of several copies of
   * indices marking the upper and lower limits of data. Indices
   * beginning with the character "p" are relative indices marking the
   * location of the data set relative to the origin the local patch of
   * the global array, all other indices are in absolute coordinates and
   * mark locations in the total global array. The indices used by this
   * routine are described below.
   *
   *       lo_loc[], hi_loc[]: The lower and upper indices of the visible
   *       block of data held by the calling processor.
   *
   *       lo_rem[], hi_rem[]: The lower and upper indices of the block
   *       of data on a remote processor or processors that is needed to
   *       fill in the calling processors ghost cells. These indices are
   *       NOT corrected for wrap-around (periodic) boundary conditions
   *       so they can be negative or greater than the array dimension
   *       values held in dims[].
   *
   *       slo_rem[], shi_rem[]: Similar to lo_rem[] and hi_rem[], except
   *       that these indices have been corrected for wrap-around
   *       boundary conditions. If lo_rem[] and hi_rem[] cross a global
   *        array boundary, as opposed to being entirely located on one
   *       side or the other of the array, then two sets of slo_rem[] and
   *       shi_rem[] will be created. One set will correspond to the
   *       block of data on one side of the global array boundary and the
   *       other set will correspond to the remaining block. This
   *       situation will only occur if the value of the ghost cell width
   *       is greater than the dimension of the visible global array
   *       data on a single processor.
   *
   *       thi_rem[], thi_rem[]: The lower and upper indices of the visible
   *       data on a remote processor.
   *
   *       plo_loc[], phi_loc[]: The indices of the local data patch that
   *       is going to be updated.
   *
   *       plo_rem[], phi_rem[]: The indices of the data patch on the
   *       remote processor that will be used to update the data on the
   *       calling processor. Note that the dimensions of the patches
   *       represented by plo_loc[], plo_rem[] and plo_loc[], phi_loc[]
   *       must be the same.
   *
   * For the case where the width of the ghost cells is more than the
   * width of the visible data held on a processor, special problems
   * arise. It now takes several updates to fill in one block of boundary
   * data and it is now necessary to keep track of where each of these
   * blocks of data go in the ghost cell region. To do this two extra
   * variables are needed. These are offset and slice. Slice is equal to
   * the width of the visible data along the dimension being updated
   * minus one coming from the remote processor. Offset is the amount
   * that this data must be moved inward from the lower boundary of the
   * ghost cell region. Another variable that is also used to handle
   * this case is imax. If this variable is set to 2, then this means
   * that the block of data that is needed to update the ghost cells
   * crosses a global array boundary and the block needs to be broken
   * up into two pieces. */

  /* if global array has no ghost cells, just return */
  if (!pnga_has_ghosts(g_a)) return;

  GA_PUSH_NAME("ga_update1_ghosts");

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  corner_flag = GA[handle].corner_flag;
  p_handle = GA[handle].p_handle;

  /* Get pointer to local memory */
  ptr_loc = GA[handle].ptr[me];
  /* obtain range of data that is held by local processor */
  pnga_distribution(g_a,me,lo_loc,hi_loc);
  /* initialize range increments and get array dimensions */
  for (idx=0; idx < ndim; idx++) {
    increment[idx] = 0;
    width[idx] = (Integer)GA[handle].width[idx];
    dims[idx] = (Integer)GA[handle].dims[idx];
    if (lo_loc[idx] == 0 && hi_loc[idx] == -1) hasData = FALSE;
  }

  /* loop over dimensions for sequential update using shift algorithm */
  for (idx=0; idx < ndim; idx++) {
    nwidth = (Integer)width[idx];

    /* Do not bother with update if nwidth is zero or processor has
       no data */
    if (nwidth != 0 && hasData) {

      /* Perform update in negative direction. Start by getting rough
         estimate of block of needed data*/
      for (i = 0; i < ndim; i++) {
        if (i == idx) {
          lo_rem[i] = lo_loc[i] - nwidth;
          hi_rem[i] = lo_loc[i] - 1;
          /* Check to see if we will need to update ghost cells using
             one or two major patches of the global array. */
          if (lo_rem[i] < 1) {
            if (hi_rem[i] > 0) {
              imax = 2;
            } else {
              imax = 1;
            }
          } else {
            imax = 1;
          }
        } else {
          lo_rem[i] = lo_loc[i];
          hi_rem[i] = hi_loc[i];
        }
      }

      for (inx = 0; inx < imax; inx++) {
        /* Check to see if boundary is being updated in one patch or two,
           adjust lower boundary accordingly. */
        for (i=0; i<ndim; i++) {
          if (imax == 2 && i == idx) {
            if (inx == 0) {
              slo_rem[i] = 1;
              shi_rem[i] = hi_rem[i];
            } else {
              slo_rem[i] = lo_rem[i] + dims[i];
              shi_rem[i] = dims[i];
            }
          } else if (i == idx) {
            if (lo_rem[i] < 1) {
              slo_rem[i] = dims[i] - nwidth + 1;
              shi_rem[i] = dims[i];
            } else {
              slo_rem[i] = lo_rem[i];
              shi_rem[i] = hi_rem[i];
            }
          } else {
            slo_rem[i] = lo_rem[i];
            shi_rem[i] = hi_rem[i];
          }
        }
        /* locate processor with this data */
        if (!pnga_locate_region(g_a, slo_rem, shi_rem, _ga_map,
            GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
            slo_rem, shi_rem, g_a);

        for (ipx = 0; ipx < np; ipx++) {
          /* Get actual coordinates of desired chunk of remote
             data as well as the actual coordinates of the local chunk
             of data that will receive the remote data (these
             coordinates take into account the presence of ghost
             cells). Start by finding out what data is actually held by
             remote processor. */
          proc_rem = GA_proclist[ipx];
          pnga_distribution(g_a, proc_rem, tlo_rem, thi_rem);
          for (i = 0; i < ndim; i++) {
            if (increment[i] == 0) {
              if (i == idx) {
                if (np == 1 && imax == 1) {
                  plo_rem[i] = thi_rem[i] - tlo_rem[i] + 1;
                  phi_rem[i] = thi_rem[i] - tlo_rem[i] + width[i];
                  plo_loc[i] = 0;
                  /*phi_loc[i] = width[i] - 1;*/
                } else {
                  if (tlo_rem[i] >= slo_rem[i]) {
                    offset = tlo_rem[i] - lo_rem[i];
                    slice = thi_rem[i] - tlo_rem[i];
                  } else {
                    offset = 0;
                    slice = thi_rem[i] - slo_rem[i];
                  }
                  if (offset < 0) offset = offset + dims[i];
                  if (offset >= dims[i]) offset = offset - dims[i];
                  plo_rem[i] = thi_rem[i] - tlo_rem[i] + width[i] - slice;
                  phi_rem[i] = thi_rem[i] - tlo_rem[i] + width[i];
                  plo_loc[i] = offset;
                  /*phi_loc[i] = offset + slice;*/
                }
              } else {
                plo_rem[i] = width[i];
                phi_rem[i] = thi_rem[i] - tlo_rem[i] + width[i];
                plo_loc[i] = width[i];
                /*phi_loc[i] = hi_loc[i] - lo_loc[i] + width[i];*/
              }
            } else {
              plo_rem[i] = 0;
              phi_rem[i] = thi_rem[i] - tlo_rem[i] + increment[i];
              plo_loc[i] = 0;
              /*phi_loc[i] = hi_loc[i] - lo_loc[i] + increment[i];*/
            }
          }

          /* Get pointer to local data buffer and remote data
             buffer as well as lists of leading dimenstions */
          gam_LocationWithGhosts(me, handle, plo_loc, &ptr_loc, ld_loc);
          gam_LocationWithGhosts(proc_rem, handle, plo_rem, &ptr_rem, ld_rem);

          /* Evaluate strides on local and remote processors */
          gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
              stride_loc);

          /* Compute the number of elements in each dimension and store
             result in count. Scale the first element in count by the
             element size. */
          gam_ComputeCount(ndim, plo_rem, phi_rem, count);
          count[0] *= size;
 
          /* get remote data */
          if (p_handle >= 0) {
            proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem];
          }
          ARMCI_GetS(ptr_rem, stride_rem, ptr_loc, stride_loc, count,
              (int)(ndim - 1), (int)proc_rem);
        }
      }

      /* Perform update in positive direction. Start by getting rough
         estimate of block of needed data*/
      for (i = 0; i < ndim; i++) {
        if (i == idx) {
          lo_rem[i] = hi_loc[i] + 1;
          hi_rem[i] = hi_loc[i] + nwidth;
          /* Check to see if we will need to update ghost cells using
             one or two major patches of the global array. */
          if (hi_rem[i] > dims[i]) {
            if (lo_rem[i] <= dims[i]) {
              imax = 2;
            } else {
              imax = 1;
            }
          } else {
            imax = 1;
          }
        } else {
          lo_rem[i] = lo_loc[i];
          hi_rem[i] = hi_loc[i];
        }
      }

      for (inx = 0; inx < imax; inx++) {
        /* Check to see if boundary is being updated in one patch or two,
           adjust lower boundary accordingly. */
        for (i=0; i<ndim; i++) {
          if (imax == 2 && i == idx) {
            if (inx == 0) {
              slo_rem[i] = lo_rem[i];
              shi_rem[i] = dims[i];
            } else {
              slo_rem[i] = 1;
              shi_rem[i] = hi_rem[i] - dims[i];
            }
          } else if (i == idx) {
            if (hi_rem[i] > dims[i]) {
              slo_rem[i] = 1;
              shi_rem[i] = nwidth;
            } else {
              slo_rem[i] = lo_rem[i];
              shi_rem[i] = hi_rem[i];
            }
          } else {
            slo_rem[i] = lo_rem[i];
            shi_rem[i] = hi_rem[i];
          }
        }
        /* locate processor with this data */
        if (!pnga_locate_region(g_a, slo_rem, shi_rem, _ga_map,
            GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
            slo_rem, shi_rem, g_a);

        for (ipx = 0; ipx < np; ipx++) {
          /* Get actual coordinates of desired chunk of remote
             data as well as the actual coordinates of the local chunk
             of data that will receive the remote data (these
             coordinates take into account the presence of ghost
             cells). Start by finding out what data is actually held by
             remote processor. */
          proc_rem = GA_proclist[ipx];
          pnga_distribution(g_a, proc_rem, tlo_rem, thi_rem);
          for (i = 0; i < ndim; i++) {
            if (increment[i] == 0) {
              if (i == idx) {
                if (np == 1 && imax == 1) {
                  plo_rem[i] = width[i];
                  phi_rem[i] = 2*width[i] - 1;
                  plo_loc[i] = hi_loc[i] - lo_loc[i] + 1 + width[i];
                  /*phi_loc[i] = hi_loc[i] - lo_loc[i] + 2*width[i];*/
                } else {
                  offset = tlo_rem[i] - hi_loc[i] - 1;
                  if (thi_rem[i] <= shi_rem[i]) {
                    slice = thi_rem[i] - tlo_rem[i];
                  } else {
                    slice = shi_rem[i] - tlo_rem[i];
                  }
                  if (offset < 0) offset = offset + dims[i];
                  if (offset >= dims[i]) offset = offset - dims[i];
                  plo_rem[i] = width[i];
                  phi_rem[i] = width[i] + slice;
                  plo_loc[i] = hi_loc[i] - lo_loc[i] + width[i] + 1 + offset;
                  /*phi_loc[i] = hi_loc[i] - lo_loc[i] + width[i] + 1
                    + offset + slice;*/
                }
              } else {
                plo_rem[i] = width[i];
                phi_rem[i] = thi_rem[i] - tlo_rem[i] + width[i];
                plo_loc[i] = width[i];
                /*phi_loc[i] = hi_loc[i] - lo_loc[i] + width[i];*/
              }
            } else {
              plo_rem[i] = 0;
              phi_rem[i] = thi_rem[i] - tlo_rem[i] + increment[i];
              plo_loc[i] = 0;
              /*phi_loc[i] = hi_loc[i] - lo_loc[i] + increment[i];*/
            }
          }

          /* Get pointer to local data buffer and remote data
             buffer as well as lists of leading dimenstions */
          gam_LocationWithGhosts(me, handle, plo_loc, &ptr_loc, ld_loc);
          gam_LocationWithGhosts(proc_rem, handle, plo_rem, &ptr_rem, ld_rem);

          /* Evaluate strides on local and remote processors */
          gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
              stride_loc);

          /* Compute the number of elements in each dimension and store
             result in count. Scale the first element in count by the
             element size. */
          gam_ComputeCount(ndim, plo_rem, phi_rem, count);
          count[0] *= size;
 
          /* get remote data */
          if (p_handle >= 0) {
            proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem];
          }
          ARMCI_GetS(ptr_rem, stride_rem, ptr_loc, stride_loc, count,
              (int)(ndim - 1), (int)proc_rem);
        }
      }
    }
    /* synchronize all processors and update increment array */
    if (idx < ndim-1) pnga_pgroup_sync(p_handle);
    if (corner_flag)
      increment[idx] = 2*nwidth;
  }

  GA_POP_NAME;
}

/*\ UTILITY FUNCTION TO MAKE SURE GHOST CELLS WIDTHS ARE
 *  LESS THAN VISIBLE DATA WIDTHS
\*/
static logical gai_check_ghost_distr(Integer g_a)
{
  Integer handle=GA_OFFSET + g_a;
  Integer idx, ndim, np, ipx;
  ndim = GA[handle].ndim;
  ipx = 0;
  for (idx = 0; idx < ndim; idx++) {
    for (np = 0; np < GA[handle].nblock[idx]; np++) {
      if (np < GA[handle].nblock[idx] - 1) {
        if (GA[handle].mapc[ipx+1]-GA[handle].mapc[ipx]+1
            <GA[handle].width[idx]) {
          return FALSE;
        }
      } else {
        if (GA[handle].dims[idx]-GA[handle].mapc[ipx]+1
            <GA[handle].width[idx]) {
          return FALSE;
        }
      }
      ipx++;
    }
  }
  return TRUE;
}

/*\ UPDATE GHOST CELLS OF GLOBAL ARRAY USING PUT CALLS
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_update2_ghosts = pnga_update2_ghosts
#endif
logical pnga_update2_ghosts(Integer g_a)
{
  Integer idx, ipx, np, handle=GA_OFFSET + g_a, proc_rem;
  Integer ntot, mask[MAXDIM];
  Integer size, ndim, i, itmp;
  Integer width[MAXDIM], dims[MAXDIM];
  Integer lo_loc[MAXDIM], hi_loc[MAXDIM];
  /*Integer tlo_loc[MAXDIM], thi_loc[MAXDIM];*/
  Integer plo_loc[MAXDIM], phi_loc[MAXDIM];
  Integer tlo_rem[MAXDIM], thi_rem[MAXDIM];
  Integer plo_rem[MAXDIM];
  Integer ld_loc[MAXDIM], ld_rem[MAXDIM];
  logical mask0;
  int stride_loc[MAXDIM], stride_rem[MAXDIM],count[MAXDIM];
  char *ptr_loc, *ptr_rem;
  Integer me = pnga_nodeid();
  Integer p_handle;

  /* if global array has no ghost cells, just return */
  if (!pnga_has_ghosts(g_a)) {
    return TRUE;
  }

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  p_handle = GA[handle].p_handle;
  /* initialize ghost cell widths and get array dimensions */
  for (idx=0; idx < ndim; idx++) {
    width[idx] = (Integer)GA[handle].width[idx];
    dims[idx] = (Integer)GA[handle].dims[idx];
  }

  /* Check to make sure that global array is well-behaved (all processors
     have data and the width of the data in each dimension is greater than
     the corresponding value in width[]). */
  if (!gai_check_ghost_distr(g_a)) return FALSE;

  GA_PUSH_NAME("ga_update2_ghosts");
  /* Get pointer to local memory */
  ptr_loc = GA[handle].ptr[me];
  /* obtain range of data that is held by local processor */
  pnga_distribution(g_a,me,lo_loc,hi_loc);

  /* evaluate total number of PUT operations that will be required */
  ntot = 1;
  for (idx=0; idx < ndim; idx++) ntot *= 3;

  /* Loop over all PUT operations. The operation corresponding to the
     mask of all zeros is left out. */
  for (ipx=0; ipx < ntot; ipx++) {
    /* Convert ipx to corresponding mask values */
    itmp = ipx;
    mask0 = TRUE;
    for (idx = 0; idx < ndim; idx++) {
      i = itmp%3;
      mask[idx] = i-1;
      if (mask[idx] != 0) mask0 = FALSE;
      itmp = (itmp-i)/3;
    }
    if (mask0) continue;

    /* check to see if ghost cell block has zero elements*/
    mask0 = FALSE;
    itmp = 0;
    for (idx = 0; idx < ndim; idx++) {
      if (mask[idx] != 0 && width[idx] == 0) mask0 = TRUE;
      if (mask[idx] != 0) itmp++;
    }
    /*if (itmp>1) mask0 = TRUE; */
    if (mask0) continue;
    /* Now that mask has been determined, find data that is to be moved
     * and identify processor to which it is going. Wrap boundaries
     * around, if necessary */
    for (idx = 0; idx < ndim; idx++) {
      if (mask[idx] == 0) {
        /*tlo_loc[idx] = lo_loc[idx];*/
        /*thi_loc[idx] = hi_loc[idx];*/
        tlo_rem[idx] = lo_loc[idx];
        thi_rem[idx] = hi_loc[idx];
      } else if (mask[idx] == -1) {
        /*tlo_loc[idx] = lo_loc[idx];*/
        /*thi_loc[idx] = lo_loc[idx]+width[idx]-1;*/
        if (lo_loc[idx] > 1) {
          tlo_rem[idx] = lo_loc[idx]-width[idx];
          thi_rem[idx] = lo_loc[idx]-1;
        } else {
          tlo_rem[idx] = dims[idx]-width[idx]+1;
          thi_rem[idx] = dims[idx];
        }
      } else if (mask[idx] == 1) {
        /*tlo_loc[idx] = hi_loc[idx]-width[idx]+1;*/
        /*thi_loc[idx] = hi_loc[idx];*/
        if (hi_loc[idx] < dims[idx]) {
          tlo_rem[idx] = hi_loc[idx] + 1;
          thi_rem[idx] = hi_loc[idx] + width[idx];
        } else {
          tlo_rem[idx] = 1;
          thi_rem[idx] = width[idx];
        }
      } else {
        fprintf(stderr,"Illegal mask value found\n");
      }
    }
    /* Locate remote processor to which data must be sent */
    if (!pnga_locate_region(g_a, tlo_rem, thi_rem, _ga_map,
       GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
       tlo_rem, thi_rem, g_a);
    if (np > 1) {
      fprintf(stderr,"More than one remote processor found\n");
    }
    /* Remote processor has been identified, now get ready to send
       data to it. Start by getting distribution on remote
       processor.*/
    proc_rem = GA_proclist[0];
    pnga_distribution(g_a, proc_rem, tlo_rem, thi_rem);
    for (idx = 0; idx < ndim; idx++) {
      if (mask[idx] == 0) {
        plo_loc[idx] = width[idx];
        phi_loc[idx] = hi_loc[idx]-lo_loc[idx]+width[idx];
        plo_rem[idx] = plo_loc[idx];
      } else if (mask[idx] == -1) {
        plo_loc[idx] = width[idx];
        phi_loc[idx] = 2*width[idx]-1;
        plo_rem[idx] = thi_rem[idx]-tlo_rem[idx]+width[idx]+1;
      } else if (mask[idx] == 1) {
        plo_loc[idx] = hi_loc[idx]-lo_loc[idx]+1;
        phi_loc[idx] = hi_loc[idx]-lo_loc[idx]+width[idx];
        plo_rem[idx] = 0;
      }
    }
    /* Get pointer to local data buffer and remote data
       buffer as well as lists of leading dimenstions */
    gam_LocationWithGhosts(me, handle, plo_loc, &ptr_loc, ld_loc);
    gam_LocationWithGhosts(proc_rem, handle, plo_rem, &ptr_rem, ld_rem);

    /* Evaluate strides on local and remote processors */
    gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
                  stride_loc);

    /* Compute the number of elements in each dimension and store
       result in count. Scale the first element in count by the
       element size. */
    gam_ComputeCount(ndim, plo_loc, phi_loc, count);
    count[0] *= size;
 
    /* put data on remote processor */
    /*ARMCI_PutS(ptr_loc, stride_loc, ptr_rem, stride_rem, count,
          (int)(ndim - 1), (int)proc_rem);*/
    if (p_handle >= 0) {
      proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem];
    }
    ARMCI_NbPutS(ptr_loc, stride_loc, ptr_rem, stride_rem, count,
          (int)(ndim - 1), (int)proc_rem, NULL); 
  }

  ARMCI_WaitAll();
  GA_POP_NAME;
  return TRUE;
}

/*\ GET INDICES ON REMOTE BLOCK IN NEGATIVE DIRECTION FOR UPDATE
\*/
static void get_remote_block_neg(Integer idx, Integer ndim, Integer *lo_loc,
                          Integer *hi_loc, Integer *slo_rem, Integer *shi_rem,
                          Integer *dims, Integer *width)
{
  Integer i, lo_rem[MAXDIM], hi_rem[MAXDIM];
  /*Start by getting rough idea of where data needs to go. */
  for (i = 0; i < ndim; i++) {
    if (i == idx) {
      lo_rem[i] = lo_loc[i] - width[i];
      hi_rem[i] = lo_loc[i] - 1;
    } else {
      lo_rem[i] = lo_loc[i];
      hi_rem[i] = hi_loc[i];
    }
  }

  /* Account for boundaries, if necessary. */
  for (i=0; i<ndim; i++) {
    if (i == idx) {
      if (lo_rem[i] < 1) {
        slo_rem[i] = dims[i] - width[i] + 1;
        shi_rem[i] = dims[i];
      } else {
        slo_rem[i] = lo_rem[i];
        shi_rem[i] = hi_rem[i];
      }
    } else {
      slo_rem[i] = lo_rem[i];
      shi_rem[i] = hi_rem[i];
    }
  }
}

/*\ GET INDICES ON REMOTE BLOCK IN POSITIVE DIRECTION FOR UPDATE
\*/
static void get_remote_block_pos(Integer idx, Integer ndim, Integer *lo_loc,
                          Integer *hi_loc, Integer *slo_rem, Integer *shi_rem,
                          Integer *dims, Integer *width)
{
  Integer i, lo_rem[MAXDIM], hi_rem[MAXDIM];
  /* Start by getting rough idea of where data needs to go. */
  for (i = 0; i < ndim; i++) {
    if (i == idx) {
      lo_rem[i] = hi_loc[i] + 1;
      hi_rem[i] = hi_loc[i] + width[i];
    } else {
      lo_rem[i] = lo_loc[i];
      hi_rem[i] = hi_loc[i];
    }
  }

  /* Account for boundaries, if necessary. */
  for (i=0; i<ndim; i++) {
    if (i == idx) {
      if (hi_rem[i] > dims[i]) {
        slo_rem[i] = 1;
        shi_rem[i] = width[i];
      } else {
        slo_rem[i] = lo_rem[i];
        shi_rem[i] = hi_rem[i];
      }
    } else {
      slo_rem[i] = lo_rem[i];
      shi_rem[i] = hi_rem[i];
    }
  }
}

/*\ UPDATE GHOST CELLS OF GLOBAL ARRAY USING SHIFT ALGORITHM AND PUT CALLS
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_update3_ghosts = pnga_update3_ghosts
#endif
logical pnga_update3_ghosts(Integer g_a)
{
  Integer idx, i, np, handle=GA_OFFSET + g_a, proc_rem;
  Integer size, ndim, nwidth, increment[MAXDIM];
  Integer width[MAXDIM];
  Integer dims[MAXDIM];
  Integer lo_loc[MAXDIM], hi_loc[MAXDIM];
  Integer plo_loc[MAXDIM]/*, phi_loc[MAXDIM]*/;
  Integer tlo_rem[MAXDIM], thi_rem[MAXDIM];
  Integer slo_rem[MAXDIM], shi_rem[MAXDIM];
  Integer plo_rem[MAXDIM], phi_rem[MAXDIM];
  Integer ld_loc[MAXDIM], ld_rem[MAXDIM];
  int stride_loc[MAXDIM], stride_rem[MAXDIM],count[MAXDIM];
  char *ptr_loc, *ptr_rem;
  Integer me = pnga_nodeid();
  Integer p_handle;

  /* This routine makes use of the shift algorithm to update data in the
   * ghost cells bounding the local block of visible data. The shift
   * algorithm starts by updating the blocks of data along the first
   * dimension by grabbing a block of data that is width[0] deep but
   * otherwise matches the  dimensions of the data residing on the
   * calling processor. The update of the second dimension, however,
   * grabs a block that is width[1] deep in the second dimension but is
   * ldim0 + 2*width[0] in the first dimensions where ldim0 is the
   * size of the visible data along the first dimension. The remaining
   * dimensions are left the same. For the next update, the width of the
   * second dimension is also increased by 2*width[1] and so on. This
   * algorith makes use of the fact that data for the dimensions that
   * have already been updated is available on each processor and can be
   * used in the updates of subsequent dimensions. The total number of
   * separate updates is 2*ndim, an update in the negative and positive
   * directions for each dimension. This implementation uses simple put
   * operations to perform the updates along each dimension with an
   * intervening synchronization call being used to make sure that the
   * necessary data is available on each processor before starting the
   * update along the next dimension.
   *
   * To perform the update, this routine makes use of several copies of
   * indices marking the upper and lower limits of data. Indices
   * beginning with the character "p" are relative indices marking the
   * location of the data set relative to the origin the local patch of
   * the global array, all other indices are in absolute coordinates and
   * mark locations in the total global array. The indices used by this
   * routine are described below.
   *
   *       lo_loc[], hi_loc[]: The lower and upper indices of the visible
   *       block of data held by the calling processor.
   *
   *       lo_rem[], hi_rem[]: The lower and upper indices of the block
   *       of data on a remote processor or processors that is needed to
   *       fill in the calling processors ghost cells. These indices are
   *       NOT corrected for wrap-around (periodic) boundary conditions
   *       so they can be negative or greater than the array dimension
   *       values held in dims[].
   *
   *       slo_rem[], shi_rem[]: Similar to lo_rem[] and hi_rem[], except
   *       that these indices have been corrected for wrap-around
   *       boundary conditions. 
   *
   *       thi_rem[], thi_rem[]: The lower and upper indices of the visible
   *       data on a remote processor.
   *
   *       plo_loc[], phi_loc[]: The indices of the local data patch that
   *       is going to be updated.
   *
   *       plo_rem[], phi_rem[]: The indices of the data patch on the
   *       remote processor that will be used to update the data on the
   *       calling processor. Note that the dimensions of the patches
   *       represented by plo_loc[], plo_rem[] and plo_loc[], phi_loc[]
   *       must be the same.
   */

  /* if global array has no ghost cells, just return */
  if (!pnga_has_ghosts(g_a)) return TRUE;

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  p_handle = GA[handle].p_handle;

  /* obtain range of data that is held by local processor */
  pnga_distribution(g_a,me,lo_loc,hi_loc);

  /* initialize range increments and get array dimensions */
  for (idx=0; idx < ndim; idx++) {
    increment[idx] = 0;
    width[idx] = (Integer)GA[handle].width[idx];
    dims[idx] = (Integer)GA[handle].dims[idx];
    if (lo_loc[idx] == 0 && hi_loc[idx] == -1) return FALSE;
  }

  /* Check to make sure that global array is well-behaved (all processors
     have data and the width of the data in each dimension is greater
     than the corresponding value in width[]. */
  if (!gai_check_ghost_distr(g_a)) return FALSE;

  GA_PUSH_NAME("ga_update3_ghosts");

  /* Get pointer to local memory */
  ptr_loc = GA[handle].ptr[me];

  /* loop over dimensions for sequential update using shift algorithm */
  for (idx=0; idx < ndim; idx++) {
    nwidth = width[idx];

    /* Do not bother with update if nwidth is zero */
    if (nwidth != 0) {

      /* Perform update in negative direction. */
      get_remote_block_neg(idx, ndim, lo_loc, hi_loc, slo_rem, shi_rem,
                           dims, width);
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rem, shi_rem, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rem, shi_rem, g_a);

      /* Get actual coordinates of desired location of remote
         data as well as the actual coordinates of the local chunk
         of data that will be sent to remote processor (these
         coordinates take into account the presence of ghost
         cells). Start by finding out what data is actually held by
         remote processor. */
      proc_rem = GA_proclist[0];
      pnga_distribution(g_a, proc_rem, tlo_rem, thi_rem);
      for (i = 0; i < ndim; i++) {
        if (increment[i] == 0) {
          if (i == idx) {
            plo_rem[i] = thi_rem[i] - tlo_rem[i] + width[i] + 1;
            phi_rem[i] = thi_rem[i] - tlo_rem[i] + 2*width[i];
            plo_loc[i] = width[i];
            /*phi_loc[i] = 2*width[i] - 1;*/
          } else {
            plo_rem[i] = width[i];
            phi_rem[i] = thi_rem[i] - tlo_rem[i] + width[i];
            plo_loc[i] = width[i];
            /*phi_loc[i] = hi_loc[i] - lo_loc[i] + width[i];*/
          }
        } else {
          plo_rem[i] = 0;
          phi_rem[i] = thi_rem[i] - tlo_rem[i] + increment[i];
          plo_loc[i] = 0;
          /*phi_loc[i] = hi_loc[i] - lo_loc[i] + increment[i];*/
        }
      }

      /* Get pointer to local data buffer and remote data
         buffer as well as lists of leading dimenstions */
      gam_LocationWithGhosts(me, handle, plo_loc, &ptr_loc, ld_loc);
      gam_LocationWithGhosts(proc_rem, handle, plo_rem, &ptr_rem, ld_rem);

      /* Evaluate strides on local and remote processors */
      gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
          stride_loc);

      /* Compute the number of elements in each dimension and store
         result in count. Scale the first element in count by the
         element size. */
      gam_ComputeCount(ndim, plo_rem, phi_rem, count);
      count[0] *= size;

      /* Put local data on remote processor */
      if (p_handle >= 0) {
        proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem];
      }
      ARMCI_PutS(ptr_loc, stride_loc, ptr_rem, stride_rem, count,
          (int)(ndim - 1), (int)proc_rem);

      /* Perform update in positive direction */
      get_remote_block_pos(idx, ndim, lo_loc, hi_loc, slo_rem, shi_rem,
                           dims, width);
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rem, shi_rem, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rem, shi_rem, g_a);

      /* Get actual coordinates of desired chunk of remote
         data as well as the actual coordinates of the local chunk
         of data that will receive the remote data (these
         coordinates take into account the presence of ghost
         cells). Start by finding out what data is actually held by
         remote processor. */
      proc_rem = GA_proclist[0];
      pnga_distribution(g_a, proc_rem, tlo_rem, thi_rem);
      for (i = 0; i < ndim; i++) {
        if (increment[i] == 0) {
          if (i == idx) {
            plo_rem[i] = 0;
            phi_rem[i] = width[i] - 1;
            plo_loc[i] = hi_loc[i] - lo_loc[i] + width[i] - 1;
            /*phi_loc[i] = hi_loc[i] - lo_loc[i] + 2*width[i] - 1;*/
          } else {
            plo_rem[i] = width[i];
            phi_rem[i] = thi_rem[i] - tlo_rem[i] + width[i];
            plo_loc[i] = width[i];
            /*phi_loc[i] = hi_loc[i] - lo_loc[i] + width[i];*/
          }
        } else {
          plo_rem[i] = 0;
          phi_rem[i] = thi_rem[i] - tlo_rem[i] + increment[i];
          plo_loc[i] = 0;
          /*phi_loc[i] = hi_loc[i] - lo_loc[i] + increment[i];*/
        }
      }

      /* Get pointer to local data buffer and remote data
         buffer as well as lists of leading dimenstions */
      gam_LocationWithGhosts(me, handle, plo_loc, &ptr_loc, ld_loc);
      gam_LocationWithGhosts(proc_rem, handle, plo_rem, &ptr_rem, ld_rem);

      /* Evaluate strides on local and remote processors */
      gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
          stride_loc);

      /* Compute the number of elements in each dimension and store
         result in count. Scale the first element in count by the
         element size. */
      gam_ComputeCount(ndim, plo_rem, phi_rem, count);
      count[0] *= size;

      /* get remote data */
      if (p_handle >= 0) {
        proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem];
      }
      ARMCI_PutS(ptr_loc, stride_loc, ptr_rem, stride_rem, count,
          (int)(ndim - 1), (int)proc_rem);
    }
    /* synchronize all processors and update increment array */
    if (idx < ndim-1) pnga_pgroup_sync(p_handle);
    increment[idx] = 2*nwidth;
  }

  GA_POP_NAME;
  return TRUE;
}

#define GHOST_PRINT 0
/*\ UPDATE GHOST CELLS OF GLOBAL ARRAY USING SHIFT ALGORITHM AND
 *  MESSAGE PASSING
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_set_update4_info = pnga_set_update4_info
#endif
logical pnga_set_update4_info(Integer g_a)
{
  Integer idx, idir, i, np, handle=GA_OFFSET + g_a;
  Integer size, buflen, buftot, *bufsize, ndim, increment[MAXDIM];
  Integer *proc_rem_snd, *proc_rem_rcv;
  Integer *length;
  Integer width[MAXDIM], dims[MAXDIM], index[MAXDIM];
  Integer lo_loc[MAXDIM], hi_loc[MAXDIM];
  Integer plo_snd[MAXDIM], phi_snd[MAXDIM];
  Integer lo_rcv[MAXDIM], hi_rcv[MAXDIM];
  Integer slo_rcv[MAXDIM], shi_rcv[MAXDIM];
  Integer plo_rcv[MAXDIM], phi_rcv[MAXDIM];
  Integer ld_loc[MAXDIM];
  int *stride_snd, *stride_rcv, *count, cache_size;
  /*int corner_flag;*/
  char **ptr_snd, **ptr_rcv, *cache;
  char *current;
  Integer me = pnga_nodeid();
  Integer p_handle;

  /* This routine sets the arrays that are used to transfer data using
   * the update4. To perform the update, this routine makes use of several
   * copies of indices marking the upper and lower limits of data. Indices
   * beginning with the character "p" are relative indices marking the
   * location of the data set relative to the origin the local patch of
   * the global array, all other indices are in absolute coordinates and
   * mark locations in the total global array. The indices used by this
   * routine are described below.
   *
   *       lo_loc[], hi_loc[]: The lower and upper indices of the visible
   *       block of data held by the calling processor.
   *
   *       lo_rcv[], hi_rcv[]: The lower and upper indices of the blocks
   *       of data that will be either sent to or received from a remote
   *       processor. These indices are NOT corrected for wrap-around
   *       (periodic) boundary conditions so they can be negative or greater
   *       than the array dimension values held in dims[].
   *
   *       slo_rcv[], shi_rcv[]: Similar to lo_rcv[] and hi_rcv[], except
   *       that these indices have been corrected for wrap-around
   *       boundary conditions.
   *
   *       plo_rcv[], phi_rcv[]: The local indices of the local data patch
   *       that receive that message from the remote processor.
   *
   *       plo_snd[], phi_snd[]: The local indices of the data patch
   *       that will be sent to the remote processor. Note that the
   *       dimensions of the patches represented by plo_rec[], plo_rec[] and
   *       plo_snd[], phi_snd[] must be the same.
   */

  /* if global array has no ghost cells, just return */
  if (!pnga_has_ghosts(g_a)) return TRUE;

  /* Check to make sure that global array is well-behaved (all processors
     have data and the width of the data in each dimension is greater
     than the corresponding value in width[]. */
  if (!gai_check_ghost_distr(g_a)) return FALSE;

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  p_handle = GA[handle].p_handle;
  cache_size = 2*sizeof(char *)+3*ndim*sizeof(int)+3*sizeof(Integer);
  cache_size = 2* ndim *((cache_size/8) + 1) + 1;
  GA[handle].cache = (double *)malloc(sizeof(double)*cache_size);
  cache = (char *)GA[handle].cache;
  /*corner_flag = GA[handle].corner_flag;*/
#if GHOST_PRINT
      printf("p[%d]a cache_size: %d\n",GAme,cache_size);
#endif

  /* initialize range increments and get array dimensions */

  pnga_distribution(g_a,me,lo_loc,hi_loc);
  for (idx=0; idx < ndim; idx++) {
    increment[idx] = 0;
    width[idx] = (Integer)GA[handle].width[idx];
    dims[idx] = (Integer)GA[handle].dims[idx];
    if (lo_loc[idx] == 0 && hi_loc[idx] == -1) {
      *(char**)cache = NULL;
      return FALSE;
    }
  }

  /* Get indices of processor in virtual grid */
  pnga_proc_topology(g_a, me, index);

  /* Try to find maximum size of message that will be sent during
   * update operations and use this to allocate memory for message
   * passing buffers. */
  buftot = 1;
  for (i=0; i<ndim; i++) {
    buftot *= (hi_loc[i]-lo_loc[i] + 1 + 2*width[i]);
  }
  buflen = 1;
  for (i = 0; i < ndim; i++) {
    idir =  hi_loc[i] - lo_loc[i] + 1;
    if (buflen < (buftot/(idir + 2*width[i]))*width[i]) {
      buflen = (buftot/(idir + 2*width[i]))*width[i];
    }
  }
  bufsize = (Integer*)cache;
#if GHOST_PRINT
      printf("p[%d]a initial pointer: %d\n",GAme,(Integer)bufsize);
      fflush(stdout);
#endif
  current = (char*)(bufsize+1);

  *bufsize = size*buflen;
#if GHOST_PRINT
      printf("p[%d]a buflen: %d size: %d bufsize: %d\n",GAme,buflen,size,*bufsize);
      fflush(stdout);
#endif

  /* loop over dimensions for sequential update using shift algorithm */
  for (idx=0; idx < ndim; idx++) {

    /* Do not bother with update if nwidth is zero */
    if (width[idx] != 0) {

      ptr_snd = (char**)current;
      ptr_rcv = (char**)(ptr_snd+1);
      proc_rem_snd = (Integer*)(ptr_rcv+1);
      proc_rem_rcv = (Integer*)(proc_rem_snd+1);
      stride_snd = (int*)(proc_rem_rcv+1);
      stride_rcv = (int*)(stride_snd+ndim);
      length = (Integer*)(stride_rcv+ndim);
      count = (int*)(length+1);
      current = (char*)(count+ndim);

      /* Perform update in negative direction. */
      get_remote_block_neg(idx, ndim, lo_loc, hi_loc, slo_rcv, shi_rcv,
                           dims, width);
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rcv, shi_rcv, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rcv, shi_rcv, g_a);
      *proc_rem_snd = GA_proclist[0];
      if (p_handle >= 0) {
        *proc_rem_snd = PGRP_LIST[p_handle].inv_map_proc_list[*proc_rem_snd];
      }

      /* Find processor from which data will be recieved */
      for (i = 0; i < ndim; i++) {
        if (i == idx) {
          lo_rcv[i] = hi_loc[i] + 1;
          hi_rcv[i] = hi_loc[i] + width[i];
        } else {
          lo_rcv[i] = lo_loc[i];
          hi_rcv[i] = hi_loc[i];
        }
      }

      /* Account for boundaries, if necessary. */
      for (i=0; i<ndim; i++) {
        if (i == idx) {
          if (hi_rcv[i] > dims[i]) {
            slo_rcv[i] = 1;
            shi_rcv[i] = width[i];
          } else {
            slo_rcv[i] = lo_rcv[i];
            shi_rcv[i] = hi_rcv[i];
          }
        } else {
          slo_rcv[i] = lo_rcv[i];
          shi_rcv[i] = hi_rcv[i];
        }
      }
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rcv, shi_rcv, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rcv, shi_rcv, g_a);
      *proc_rem_rcv = GA_proclist[0];
      if (p_handle >= 0) {
        *proc_rem_rcv = PGRP_LIST[p_handle].inv_map_proc_list[*proc_rem_rcv];
      }

      /* Get actual coordinates of chunk of data that will be sent to
       * remote processor as well as coordinates of the array space that
       * will receive data from remote processor. */
      for (i = 0; i < ndim; i++) {
        if (increment[i] == 0) {
          if (i == idx) {
            plo_snd[i] = width[i];
            phi_snd[i] = 2*width[i] - 1;
            plo_rcv[i] = hi_loc[i] - lo_loc[i] + width[i] + 1;
            phi_rcv[i] = hi_loc[i] - lo_loc[i] + 2*width[i];
          } else {
            plo_snd[i] = width[i];
            phi_snd[i] = hi_loc[i] - lo_loc[i] + width[i];
            plo_rcv[i] = width[i];
            phi_rcv[i] = hi_loc[i] - lo_loc[i] + width[i];
          }
        } else {
          plo_rcv[i] = 0;
          phi_rcv[i] = hi_loc[i] - lo_loc[i] + increment[i];
          plo_snd[i] = 0;
          phi_snd[i] = hi_loc[i] - lo_loc[i] + increment[i];
        }
      }

      /* Get pointer to local data buffer and remote data
         buffer as well as lists of leading dimenstions */
      gam_LocationWithGhosts(me, handle, plo_snd, ptr_snd, ld_loc);
#if GHOST_PRINT
      printf("p[%d]a 1: plo_snd[0]: %d plo_snd[1]: %d ptr_snd: %d\n",
          GAme, plo_snd[0], plo_snd[1], (Integer)*ptr_snd);
      fflush(stdout);
#endif
      gam_LocationWithGhosts(me, handle, plo_rcv, ptr_rcv, ld_loc);
#if GHOST_PRINT
      printf("p[%d]a 1: plo_rcv[0]: %d plo_rcv[1]: %d ptr_rcv: %d\n",
          GAme, plo_rcv[0], plo_rcv[1], (Integer)*ptr_rcv);
      fflush(stdout);
#endif

      /* Evaluate strides for send and recieve */
      gam_setstride(ndim, size, ld_loc, ld_loc, stride_rcv,
          stride_snd);

      /* Compute the number of elements in each dimension and store
         result in count. Scale the first element in count by the
         element size. */
      gam_ComputeCount(ndim, plo_rcv, phi_rcv, count);
      gam_CountElems(ndim, plo_snd, phi_snd, length);
      *length *= size;
      count[0] *= size;

#if GHOST_PRINT
      printf("p[%d]a 1: length: %d bufsize: %d proc_rem_snd: %d proc_rem_rcv: %d\n",
          GAme, *length, *bufsize, (int)*proc_rem_snd, (int)*proc_rem_rcv);
      printf("p[%d]a 1: count[0]: %d stride_rcv[0]: %d stride_rcv[1]: %d\n",
          GAme, count[0], stride_rcv[0],stride_rcv[1]);
      printf("p[%d]a 1: count[1]: %d stride_rcv[2]: %d stride_rcv[3]: %d\n",
          GAme, count[1], stride_rcv[2],stride_rcv[3]);
      printf("p[%d]a 1: count[2]: %d stride_snd[0]: %d stride_snd[1]: %d\n",
          GAme, count[2], stride_snd[0],stride_snd[1]);
      printf("p[%d]a 1: count[3]: %d stride_snd[2]: %d stride_snd[3]: %d\n",
          GAme, count[3], stride_snd[2],stride_snd[3]);
      fflush(stdout);
#endif

      ptr_snd = (char**)current;
      ptr_rcv = (char**)(ptr_snd+1);
      proc_rem_snd = (Integer*)(ptr_rcv+1);
      proc_rem_rcv = (Integer*)(proc_rem_snd+1);
      stride_snd = (int*)(proc_rem_rcv+1);
      stride_rcv = (int*)(stride_snd+ndim);
      length = (Integer*)(stride_rcv+ndim);
      count = (int*)(length+1);
      current = (char*)(count+ndim);

      /* Find parameters for message in positive direction. */
      get_remote_block_pos(idx, ndim, lo_loc, hi_loc, slo_rcv, shi_rcv,
                           dims, width);
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rcv, shi_rcv, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rcv, shi_rcv, g_a);
      *proc_rem_snd = GA_proclist[0];
      if (p_handle >= 0) {
        *proc_rem_snd = PGRP_LIST[p_handle].inv_map_proc_list[*proc_rem_snd];
      }

      /* Find processor from which data will be recieved */
      for (i = 0; i < ndim; i++) {
        if (i == idx) {
          lo_rcv[i] = lo_loc[i] - width[i];
          hi_rcv[i] = lo_loc[i] - 1;
        } else {
          lo_rcv[i] = lo_loc[i];
          hi_rcv[i] = hi_loc[i];
        }
      }

      /* Account for boundaries, if necessary. */
      for (i=0; i<ndim; i++) {
        if (i == idx) {
          if (hi_rcv[i] < 1) {
            slo_rcv[i] = dims[i] - width[i] + 1;
            shi_rcv[i] = dims[i];
          } else {
            slo_rcv[i] = lo_rcv[i];
            shi_rcv[i] = hi_rcv[i];
          }
        } else {
          slo_rcv[i] = lo_rcv[i];
          shi_rcv[i] = hi_rcv[i];
        }
      }
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rcv, shi_rcv, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rcv, shi_rcv, g_a);
      *proc_rem_rcv = GA_proclist[0];
      if (p_handle >= 0) {
        *proc_rem_rcv = PGRP_LIST[p_handle].inv_map_proc_list[*proc_rem_rcv];
      }
      /* Get actual coordinates of chunk of data that will be sent to
       * remote processor as well as coordinates of the array space that
       * will receive data from remote processor. */
      for (i = 0; i < ndim; i++) {
        if (increment[i] == 0) {
          if (i == idx) {
            plo_snd[i] = hi_loc[i] - lo_loc[i] + 1;
            phi_snd[i] = hi_loc[i] - lo_loc[i] + width[i];
            plo_rcv[i] = 0;
            phi_rcv[i] = width[i] - 1;
          } else {
            plo_snd[i] = width[i];
            phi_snd[i] = hi_loc[i] - lo_loc[i] + width[i];
            plo_rcv[i] = width[i];
            phi_rcv[i] = hi_loc[i] - lo_loc[i] + width[i];
          }
        } else {
          plo_rcv[i] = 0;
          phi_rcv[i] = hi_loc[i] - lo_loc[i] + increment[i];
          plo_snd[i] = 0;
          phi_snd[i] = hi_loc[i] - lo_loc[i] + increment[i];
        }
      }

      /* Get pointer to local data buffer and remote data
         buffer as well as lists of leading dimenstions */
      gam_LocationWithGhosts(me, handle, plo_snd, ptr_snd, ld_loc);
#if GHOST_PRINT
      printf("p[%d]a 2: plo_snd[0]: %d plo_snd[1]: %d ptr_snd: %d\n",
          GAme, plo_snd[0], plo_snd[1], (Integer)*ptr_snd);
      fflush(stdout);
#endif
      gam_LocationWithGhosts(me, handle, plo_rcv, ptr_rcv, ld_loc);
#if GHOST_PRINT
      printf("p[%d]a 2: plo_rcv[0]: %d plo_rcv[1]: %d ptr_rcv: %d\n",
          GAme, plo_rcv[0], plo_rcv[1], (Integer)*ptr_rcv);
      fflush(stdout);
#endif

      /* Evaluate strides for send and recieve */
      gam_setstride(ndim, size, ld_loc, ld_loc, stride_rcv,
          stride_snd);

      /* Compute the number of elements in each dimension and store
         result in count. Scale the first element in count by the
         element size. */
      gam_ComputeCount(ndim, plo_rcv, phi_rcv, count);
      gam_CountElems(ndim, plo_snd, phi_snd, length);
      *length *= size;
      count[0] *= size;
#if GHOST_PRINT
      printf("p[%d]a 2: length: %d bufsize: %d proc_rem_snd: %d proc_rem_rcv: %d\n",
          GAme, *length, *bufsize, (int)*proc_rem_snd, (int)*proc_rem_rcv);
      printf("p[%d]a 2: count[0]: %d stride_rcv[0]: %d stride_rcv[1]: %d\n",
          GAme, count[0], stride_rcv[0],stride_rcv[1]);
      printf("p[%d]a 2: count[1]: %d stride_rcv[2]: %d stride_rcv[3]: %d\n",
          GAme, count[1], stride_rcv[2],stride_rcv[3]);
      printf("p[%d]a 2: count[2]: %d stride_snd[0]: %d stride_snd[1]: %d\n",
          GAme, count[2], stride_snd[0],stride_snd[1]);
      printf("p[%d]a 2: count[3]: %d stride_snd[2]: %d stride_snd[3]: %d\n",
          GAme, count[3], stride_snd[2],stride_snd[3]);
      fflush(stdout);
#endif

    }
    if (GA[handle].corner_flag)
      increment[idx] = 2*width[idx];
  }
#if GHOST_PRINT
      printf("p[%d]a final pointer: %d\n",GAme,(Integer)(Integer*)current);
      fflush(stdout);
#endif
  return TRUE;
}

/*\ UPDATE GHOST CELLS OF GLOBAL ARRAY USING SHIFT ALGORITHM AND
 *  MESSAGE PASSING
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_update4_ghosts = pnga_update4_ghosts
#endif
logical pnga_update4_ghosts(Integer g_a)
{
  Integer idx, i, handle=GA_OFFSET + g_a;
  Integer *size, bufsize, buflen, ndim, elemsize;
  Integer *proc_rem_snd, *proc_rem_rcv, pmax;
  Integer msgcnt, *length;
  Integer index[MAXDIM], width[MAXDIM];
  int *stride_snd, *stride_rcv, *count, msglen;
  char **ptr_snd, **ptr_rcv, *cache, *current;
  char send_name[32], rcv_name[32];
  void *snd_ptr, *rcv_ptr, *snd_ptr_orig, *rcv_ptr_orig;
  Integer me = pnga_nodeid();
  /*Integer p_handle;*/

  /* This routine makes use of the shift algorithm to update data in the
   * ghost cells bounding the local block of visible data. The shift
   * algorithm starts by updating the blocks of data along the first
   * dimension by grabbing a block of data that is width[0] deep but
   * otherwise matches the  dimensions of the data residing on the
   * calling processor. The update of the second dimension, however,
   * grabs a block that is width[1] deep in the second dimension but is
   * ldim0 + 2*width[0] in the first dimensions where ldim0 is the
   * size of the visible data along the first dimension. The remaining
   * dimensions are left the same. For the next update, the width of the
   * second dimension is also increased by 2*width[1] and so on. This
   * algorith makes use of the fact that data for the dimensions that
   * have already been updated is available on each processor and can be
   * used in the updates of subsequent dimensions. The total number of
   * separate updates is 2*ndim, an update in the negative and positive
   * directions for each dimension.
   *
   * This implementation make use of explicit message passing to perform
   * the update. Separate message types for the updates in each coordinate
   * direction are used to maintain synchronization locally and to
   * guarantee that the data is present before the updates in a new
   * coordinate direction take place.
   */

  /* if global array has no ghost cells, just return */
  if (!pnga_has_ghosts(g_a)) return TRUE;

  ndim = GA[handle].ndim;
  /*p_handle = GA[handle].p_handle;*/
  cache = (char *)GA[handle].cache;
  elemsize = GA[handle].elemsize;
  for (i=0; i<ndim; i++) {
    width[i] = (Integer)GA[handle].width[i];
  }

  GA_PUSH_NAME("ga_update4_ghosts");
  msgcnt = 0;

  /* Get indices of processor in virtual grid */
  pnga_proc_topology(g_a, me, index);

  size = (Integer*)cache;
  current = (char*)(size+1);
  bufsize = *size;
  buflen = bufsize/elemsize;
#if GHOST_PRINT
  printf("p[%d] bufsize: %d buflen: %d\n", GAme, bufsize, buflen);
  fflush(stdout);
#endif

  strcpy(send_name,"send_buffer");
  strcpy(rcv_name,"receive_buffer");
  snd_ptr_orig = snd_ptr = ga_malloc(buflen, GA[handle].type, send_name);
  rcv_ptr_orig = rcv_ptr = ga_malloc(buflen, GA[handle].type, rcv_name);

  /* loop over dimensions for sequential update using shift algorithm */
  for (idx=0; idx < ndim; idx++) {

    /* Do not bother with update if nwidth is zero */
    if (width[idx] != 0) {

      /* send messages in negative direction */
      snd_ptr = snd_ptr_orig;
      rcv_ptr = rcv_ptr_orig;

      ptr_snd = (char**)current;
      ptr_rcv = (char**)(ptr_snd+1);
      proc_rem_snd = (Integer*)(ptr_rcv+1);
      proc_rem_rcv = (Integer*)(proc_rem_snd+1);
      stride_snd = (int*)(proc_rem_rcv+1);
      stride_rcv = (int*)(stride_snd+ndim);
      length = (Integer*)(stride_rcv+ndim);
      count = (int*)(length+1);
      current = (char*)(count+ndim);

#if GHOST_PRINT
      printf("p[%d] 1: ptr_snd: %d ptr_rcv: %d\n", GAme, (Integer)*ptr_snd,
              (Integer)*ptr_rcv);
      printf("p[%d] 1: length: %d proc_rem_snd: %d proc_rem_rcv: %d\n",
          GAme, (int)*length, (int)*proc_rem_snd, (int)*proc_rem_rcv);
      printf("p[%d] 1: count[0]: %d stride_rcv[0]: %d stride_rcv[1]: %d\n",
          GAme, count[0], stride_rcv[0],stride_rcv[1]);
      printf("p[%d] 1: count[1]: %d stride_rcv[2]: %d stride_rcv[3]: %d\n",
          GAme, count[1], stride_rcv[2],stride_rcv[3]);
      printf("p[%d] 1: count[2]: %d stride_snd[0]: %d stride_snd[1]: %d\n",
          GAme, count[2], stride_snd[0],stride_snd[1]);
      printf("p[%d] 1: count[3]: %d stride_snd[2]: %d stride_snd[3]: %d\n",
          GAme, count[3], stride_snd[2],stride_snd[3]);
      printf("p[%d] 1: snd_ptr: %d rcv_ptr: %d\n", GAme, (Integer)snd_ptr,
          (Integer)rcv_ptr);
      fflush(stdout);
#endif
      /* Fill send buffer with data. */
      armci_write_strided(*ptr_snd, (int)ndim-1, stride_snd, count, snd_ptr);
#if GHOST_PRINT
      printf("p[%d] completed armci_write_strided\n",GAme);
      fflush(stdout);
#endif

      /* Send Messages. If processor has odd index in direction idx, it
       * sends message first, if processor has even index it receives
       * message first. Then process is reversed. Also need to account
       * for whether or not there are an odd number of processors along
       * update direction. */

      if (GAme != *proc_rem_snd) {
        if (GA[handle].nblock[idx]%2 == 0) {
          if (index[idx]%2 != 0) {
            armci_msg_snd(msgcnt, snd_ptr, *length, *proc_rem_snd);
          } else {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, *proc_rem_rcv);
          }
          if (index[idx]%2 != 0) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, *proc_rem_rcv);
          } else {
            armci_msg_snd(msgcnt, snd_ptr, *length, *proc_rem_snd);
          }
        } else {
          pmax = GA[handle].nblock[idx] - 1;
          if (index[idx]%2 != 0) {
            armci_msg_snd(msgcnt, snd_ptr, *length, *proc_rem_snd);
          } else if (index[idx] != pmax) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, *proc_rem_rcv);
          }
          if (index[idx]%2 != 0) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, *proc_rem_rcv);
          } else if (index[idx] != 0) {
            armci_msg_snd(msgcnt, snd_ptr, *length, *proc_rem_snd);
          }
          /* make up for odd processor at end of string */
          if (index[idx] == 0) {
            armci_msg_snd(msgcnt, snd_ptr, *length, *proc_rem_snd);
          }
          if (index[idx] == pmax) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, *proc_rem_rcv);
          }
        }
      } else {
        rcv_ptr = snd_ptr;
      }
      msgcnt++;
      /* copy data back into global array */
      armci_read_strided(*ptr_rcv, (int)ndim-1, stride_rcv, count, rcv_ptr);
#if GHOST_PRINT
      printf("p[%d] completed armci_read_strided\n",GAme);
      fflush(stdout);
#endif

      /* send messages in positive direction */
      snd_ptr = snd_ptr_orig;
      rcv_ptr = rcv_ptr_orig;

      ptr_snd = (char**)current;
      ptr_rcv = (char**)(ptr_snd+1);
      proc_rem_snd = (Integer*)(ptr_rcv+1);
      proc_rem_rcv = (Integer*)(proc_rem_snd+1);
      stride_snd = (int*)(proc_rem_rcv+1);
      stride_rcv = (int*)(stride_snd+ndim);
      length = (Integer*)(stride_rcv+ndim);
      count = (int*)(length+1);
      current = (char*)(count+ndim);

#if GHOST_PRINT
      printf("p[%d] 2: ptr_snd: %d ptr_rcv: %d\n", GAme, (Integer)*ptr_snd,
              (Integer)*ptr_rcv);
      printf("p[%d] 2: length: %d proc_rem_snd: %d proc_rem_rcv: %d\n",
          GAme, (int)*length, (int)*proc_rem_snd, (int)*proc_rem_rcv);
      printf("p[%d] 2: count[0]: %d stride_rcv[0]: %d stride_rcv[1]: %d\n",
          GAme, count[0], stride_rcv[0],stride_rcv[1]);
      printf("p[%d] 2: count[1]: %d stride_rcv[2]: %d stride_rcv[3]: %d\n",
          GAme, count[1], stride_rcv[2],stride_rcv[3]);
      printf("p[%d] 2: count[2]: %d stride_snd[0]: %d stride_snd[1]: %d\n",
          GAme, count[2], stride_snd[0],stride_snd[1]);
      printf("p[%d] 2: count[3]: %d stride_snd[2]: %d stride_snd[3]: %d\n",
          GAme, count[3], stride_snd[2],stride_snd[3]);
      printf("p[%d] 2: snd_ptr: %d rcv_ptr: %d\n", GAme, (Integer)snd_ptr,
          (Integer)rcv_ptr);
      fflush(stdout);
#endif
      /* Fill send buffer with data. */
      armci_write_strided(*ptr_snd, (int)ndim-1, stride_snd, count, snd_ptr);

      /* Send Messages. If processor has odd index in direction idx, it
       * sends message first, if processor has even index it receives
       * message first. Then process is reversed. Also need to account
       * for whether or not there are an odd number of processors along
       * update direction. */

      if (GAme != *proc_rem_rcv) {
        if (GA[handle].nblock[idx]%2 == 0) {
          if (index[idx]%2 != 0) {
            armci_msg_snd(msgcnt, snd_ptr, *length, *proc_rem_snd);
          } else {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, *proc_rem_rcv);
          }
          if (index[idx]%2 != 0) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, *proc_rem_rcv);
          } else {
            armci_msg_snd(msgcnt, snd_ptr, *length, *proc_rem_snd);
          }
        } else {
          pmax = GA[handle].nblock[idx] - 1;
          if (index[idx]%2 != 0) {
            armci_msg_snd(msgcnt, snd_ptr, *length, *proc_rem_snd);
          } else if (index[idx] != 0) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, *proc_rem_rcv);
          }
          if (index[idx]%2 != 0) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, *proc_rem_rcv);
          } else if (index[idx] != pmax) {
            armci_msg_snd(msgcnt, snd_ptr, *length, *proc_rem_snd);
          }
          /* make up for odd processor at end of string */
          if (index[idx] == pmax) {
            armci_msg_snd(msgcnt, snd_ptr, *length, *proc_rem_snd);
          }
          if (index[idx] == 0) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, *proc_rem_rcv);
          }
        }
      } else {
        rcv_ptr = snd_ptr;
      }
      /* copy data back into global array */
      armci_read_strided(*ptr_rcv, (int)ndim-1, stride_rcv, count, rcv_ptr);
      msgcnt++;
    }
  }

  ga_free(rcv_ptr_orig);
  ga_free(snd_ptr_orig);
  GA_POP_NAME;
  return TRUE;
}

/*\ UPDATE GHOST CELLS OF GLOBAL ARRAY USING SHIFT ALGORITHM AND
 *  MESSAGE PASSING
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_update44_ghosts = pnga_update44_ghosts
#endif
logical pnga_update44_ghosts(Integer g_a)
{
  Integer idx, idir, i, np, handle=GA_OFFSET + g_a;
  Integer size, buflen, buftot, bufsize, ndim, increment[MAXDIM];
  Integer proc_rem_snd, proc_rem_rcv, pmax;
  Integer msgcnt, length;
  Integer width[MAXDIM], dims[MAXDIM], index[MAXDIM];
  Integer lo_loc[MAXDIM], hi_loc[MAXDIM];
  Integer plo_snd[MAXDIM], phi_snd[MAXDIM];
  Integer lo_rcv[MAXDIM], hi_rcv[MAXDIM];
  Integer slo_rcv[MAXDIM], shi_rcv[MAXDIM];
  Integer plo_rcv[MAXDIM], phi_rcv[MAXDIM];
  Integer ld_loc[MAXDIM];
  int msglen;
  int stride_snd[MAXDIM], stride_rcv[MAXDIM],count[MAXDIM];
  char *ptr_snd, *ptr_rcv;
  char send_name[32], rcv_name[32];
  void *snd_ptr, *rcv_ptr, *snd_ptr_orig, *rcv_ptr_orig;
  Integer me = pnga_nodeid();
  Integer p_handle;

  /* This routine makes use of the shift algorithm to update data in the
   * ghost cells bounding the local block of visible data. The shift
   * algorithm starts by updating the blocks of data along the first
   * dimension by grabbing a block of data that is width[0] deep but
   * otherwise matches the  dimensions of the data residing on the
   * calling processor. The update of the second dimension, however,
   * grabs a block that is width[1] deep in the second dimension but is
   * ldim0 + 2*width[0] in the first dimensions where ldim0 is the
   * size of the visible data along the first dimension. The remaining
   * dimensions are left the same. For the next update, the width of the
   * second dimension is also increased by 2*width[1] and so on. This
   * algorith makes use of the fact that data for the dimensions that
   * have already been updated is available on each processor and can be
   * used in the updates of subsequent dimensions. The total number of
   * separate updates is 2*ndim, an update in the negative and positive
   * directions for each dimension.
   *
   * This implementation make use of explicit message passing to perform
   * the update. Separate message types for the updates in each coordinate
   * direction are used to maintain synchronization locally and to
   * guarantee that the data is present before the updates in a new
   * coordinate direction take place.
   *
   * To perform the update, this routine makes use of several copies of
   * indices marking the upper and lower limits of data. Indices
   * beginning with the character "p" are relative indices marking the
   * location of the data set relative to the origin the local patch of
   * the global array, all other indices are in absolute coordinates and
   * mark locations in the total global array. The indices used by this
   * routine are described below.
   *
   *       lo_loc[], hi_loc[]: The lower and upper indices of the visible
   *       block of data held by the calling processor.
   *
   *       lo_rcv[], hi_rcv[]: The lower and upper indices of the blocks
   *       of data that will be either sent to or received from a remote
   *       processor. These indices are NOT corrected for wrap-around
   *       (periodic) boundary conditions so they can be negative or greater
   *       than the array dimension values held in dims[].
   *
   *       slo_rcv[], shi_rcv[]: Similar to lo_rcv[] and hi_rcv[], except
   *       that these indices have been corrected for wrap-around
   *       boundary conditions.
   *
   *       plo_rcv[], phi_rcv[]: The local indices of the local data patch
   *       that receive that message from the remote processor.
   *
   *       plo_snd[], phi_snd[]: The local indices of the data patch
   *       that will be sent to the remote processor. Note that the
   *       dimensions of the patches represented by plo_rec[], plo_rec[] and
   *       plo_snd[], phi_snd[] must be the same.
   */

  /* if global array has no ghost cells, just return */
  if (!pnga_has_ghosts(g_a)) return TRUE;

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  p_handle = GA[handle].p_handle;

  /* initialize range increments and get array dimensions */
  for (idx=0; idx < ndim; idx++) {
    increment[idx] = 0;
    width[idx] = (Integer)GA[handle].width[idx];
    dims[idx] = (Integer)GA[handle].dims[idx];
  }

  /* Check to make sure that global array is well-behaved (all processors
     have data and the width of the data in each dimension is greater
     than the corresponding value in width[]. */
  if (!gai_check_ghost_distr(g_a)) return FALSE;

  GA_PUSH_NAME("ga_update4_ghosts");
  msgcnt = 0;

  /* obtain range of data that is held by local processor */
  pnga_distribution(g_a,me,lo_loc,hi_loc);
  /* Get indices of processor in virtual grid */
  pnga_proc_topology(g_a, me, index);

  /* Try to find maximum size of message that will be sent during
   * update operations and use this to allocate memory for message
   * passing buffers. */
  buftot = 1;
  for (i=0; i<ndim; i++) {
    buftot *= (hi_loc[i]-lo_loc[i] + 1 + 2*width[i]);
  }
  buflen = 1;
  for (i = 0; i < ndim; i++) {
    idir =  hi_loc[i] - lo_loc[i] + 1;
    if (buflen < (buftot/(idir + 2*width[i]))*width[i]) {
      buflen = (buftot/(idir + 2*width[i]))*width[i];
    }
  }
  bufsize = size*buflen;
  strcpy(send_name,"send_buffer");
  strcpy(rcv_name,"receive_buffer");
  snd_ptr_orig = snd_ptr = ga_malloc(buflen, GA[handle].type, send_name);
  rcv_ptr_orig = rcv_ptr = ga_malloc(buflen, GA[handle].type, rcv_name);

  /* loop over dimensions for sequential update using shift algorithm */
  for (idx=0; idx < ndim; idx++) {

    /* Do not bother with update if nwidth is zero */
    if (width[idx] != 0) {

      /* Perform update in negative direction. */
      get_remote_block_neg(idx, ndim, lo_loc, hi_loc, slo_rcv, shi_rcv,
                           dims, width);
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rcv, shi_rcv, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rcv, shi_rcv, g_a);
      proc_rem_snd = GA_proclist[0];
      if (p_handle >= 0) {
        proc_rem_snd = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem_snd];
      }

      /* Find processor from which data will be recieved */
      for (i = 0; i < ndim; i++) {
        if (i == idx) {
          lo_rcv[i] = hi_loc[i] + 1;
          hi_rcv[i] = hi_loc[i] + width[i];
        } else {
          lo_rcv[i] = lo_loc[i];
          hi_rcv[i] = hi_loc[i];
        }
      }

      /* Account for boundaries, if necessary. */
      for (i=0; i<ndim; i++) {
        if (i == idx) {
          if (hi_rcv[i] > dims[i]) {
            slo_rcv[i] = 1;
            shi_rcv[i] = width[i];
          } else {
            slo_rcv[i] = lo_rcv[i];
            shi_rcv[i] = hi_rcv[i];
          }
        } else {
          slo_rcv[i] = lo_rcv[i];
          shi_rcv[i] = hi_rcv[i];
        }
      }
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rcv, shi_rcv, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rcv, shi_rcv, g_a);
      proc_rem_rcv = GA_proclist[0];
      if (p_handle >= 0) {
        proc_rem_rcv = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem_rcv];
      }

      /* Get actual coordinates of chunk of data that will be sent to
       * remote processor as well as coordinates of the array space that
       * will receive data from remote processor. */
      for (i = 0; i < ndim; i++) {
        if (increment[i] == 0) {
          if (i == idx) {
            plo_snd[i] = width[i];
            phi_snd[i] = 2*width[i] - 1;
            plo_rcv[i] = hi_loc[i] - lo_loc[i] + width[i] + 1;
            phi_rcv[i] = hi_loc[i] - lo_loc[i] + 2*width[i];
          } else {
            plo_snd[i] = width[i];
            phi_snd[i] = hi_loc[i] - lo_loc[i] + width[i];
            plo_rcv[i] = width[i];
            phi_rcv[i] = hi_loc[i] - lo_loc[i] + width[i];
          }
        } else {
          plo_rcv[i] = 0;
          phi_rcv[i] = hi_loc[i] - lo_loc[i] + increment[i];
          plo_snd[i] = 0;
          phi_snd[i] = hi_loc[i] - lo_loc[i] + increment[i];
        }
      }

      /* Get pointer to local data buffer and remote data
         buffer as well as lists of leading dimenstions */
      gam_LocationWithGhosts(me, handle, plo_snd, &ptr_snd, ld_loc);
#if GHOST_PRINT
      printf("p[%d] 1: plo_snd[0]: %d plo_snd[1]: %d ptr_snd: %d\n",
          GAme, plo_snd[0], plo_snd[1], (Integer)ptr_snd);
      fflush(stdout);
#endif
      gam_LocationWithGhosts(me, handle, plo_rcv, &ptr_rcv, ld_loc);
#if GHOST_PRINT
      printf("p[%d] 1: plo_rcv[0]: %d plo_rcv[1]: %d ptr_rcv: %d\n",
          GAme, plo_rcv[0], plo_rcv[1], (Integer)ptr_rcv);
      fflush(stdout);
#endif

      /* Evaluate strides for send and recieve */
      gam_setstride(ndim, size, ld_loc, ld_loc, stride_rcv,
          stride_snd);

      /* Compute the number of elements in each dimension and store
         result in count. Scale the first element in count by the
         element size. */
      gam_ComputeCount(ndim, plo_rcv, phi_rcv, count);
      gam_CountElems(ndim, plo_snd, phi_snd, &length);
      length *= size;
      count[0] *= size;

      /* Fill send buffer with data. */
#if GHOST_PRINT
      printf("p[%d]b 1: ptr_snd: %d ptr_rcv: %d\n", GAme, (Integer)ptr_snd,
              (Integer)ptr_rcv);
      printf("p[%d]b 1: length: %d proc_rem_snd: %d proc_rem_rcv: %d\n",
          GAme, (int)length, (int)proc_rem_snd, (int)proc_rem_rcv);
      printf("p[%d]b 1: count[0]: %d stride_rcv[0]: %d stride_rcv[1]: %d\n",
          GAme, count[0], stride_rcv[0],stride_rcv[1]);
      printf("p[%d]b 1: count[1]: %d stride_rcv[2]: %d stride_rcv[3]: %d\n",
          GAme, count[1], stride_rcv[2],stride_rcv[3]);
      printf("p[%d]b 1: count[2]: %d stride_snd[0]: %d stride_snd[1]: %d\n",
          GAme, count[2], stride_snd[0],stride_snd[1]);
      printf("p[%d]b 1: count[3]: %d stride_snd[2]: %d stride_snd[3]: %d\n",
          GAme, count[3], stride_snd[2],stride_snd[3]);
      printf("p[%d]b 1: snd_ptr: %d rcv_ptr: %d\n", GAme, (Integer)snd_ptr,
          (Integer)rcv_ptr);
      fflush(stdout);
#endif
      armci_write_strided(ptr_snd, (int)ndim-1, stride_snd, count, snd_ptr);

      /* Send Messages. If processor has odd index in direction idx, it
       * sends message first, if processor has even index it receives
       * message first. Then process is reversed. Also need to account
       * for whether or not there are an odd number of processors along
       * update direction. */

#if GHOST_PRINT
      printf("p[%d] 1: msgcnt: %d length: %d bufsize: %d proc_rem_snd: %d proc_rem_rcv: %d\n",
          GAme, msgcnt, length, bufsize, (int)proc_rem_snd, (int)proc_rem_rcv);
      fflush(stdout);
#endif
      snd_ptr = snd_ptr_orig;
      rcv_ptr = rcv_ptr_orig;
      if (GAme != proc_rem_snd) {
        if (GA[handle].nblock[idx]%2 == 0) {
          if (index[idx]%2 != 0) {
            armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
          } else {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
          }
          if (index[idx]%2 != 0) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
          } else {
            armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
          }
        } else {
          pmax = GA[handle].nblock[idx] - 1;
          if (index[idx]%2 != 0) {
            armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
          } else if (index[idx] != pmax) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
          }
          if (index[idx]%2 != 0) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
          } else if (index[idx] != 0) {
            armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
          }
          /* make up for odd processor at end of string */
          if (index[idx] == 0) {
            armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
          }
          if (index[idx] == pmax) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
          }
        }
      } else {
        rcv_ptr = snd_ptr;
      }
      msgcnt++;
      /* copy data back into global array */
      armci_read_strided(ptr_rcv, (int)ndim-1, stride_rcv, count, rcv_ptr);

      /* Find parameters for message in positive direction. */
      get_remote_block_pos(idx, ndim, lo_loc, hi_loc, slo_rcv, shi_rcv,
                           dims, width);
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rcv, shi_rcv, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rcv, shi_rcv, g_a);
      proc_rem_snd = GA_proclist[0];
      if (p_handle >= 0) {
        proc_rem_snd = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem_snd];
      }

      /* Find processor from which data will be recieved */
      for (i = 0; i < ndim; i++) {
        if (i == idx) {
          lo_rcv[i] = lo_loc[i] - width[i];
          hi_rcv[i] = lo_loc[i] - 1;
        } else {
          lo_rcv[i] = lo_loc[i];
          hi_rcv[i] = hi_loc[i];
        }
      }

      /* Account for boundaries, if necessary. */
      for (i=0; i<ndim; i++) {
        if (i == idx) {
          if (hi_rcv[i] < 1) {
            slo_rcv[i] = dims[i] - width[i] + 1;
            shi_rcv[i] = dims[i];
          } else {
            slo_rcv[i] = lo_rcv[i];
            shi_rcv[i] = hi_rcv[i];
          }
        } else {
          slo_rcv[i] = lo_rcv[i];
          shi_rcv[i] = hi_rcv[i];
        }
      }
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rcv, shi_rcv, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rcv, shi_rcv, g_a);
      proc_rem_rcv = GA_proclist[0];
      if (p_handle >= 0) {
        proc_rem_rcv = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem_rcv];
      }
      /* Get actual coordinates of chunk of data that will be sent to
       * remote processor as well as coordinates of the array space that
       * will receive data from remote processor. */
      for (i = 0; i < ndim; i++) {
        if (increment[i] == 0) {
          if (i == idx) {
            plo_snd[i] = hi_loc[i] - lo_loc[i] + 1;
            phi_snd[i] = hi_loc[i] - lo_loc[i] + width[i];
            plo_rcv[i] = 0;
            phi_rcv[i] = width[i] - 1;
          } else {
            plo_snd[i] = width[i];
            phi_snd[i] = hi_loc[i] - lo_loc[i] + width[i];
            plo_rcv[i] = width[i];
            phi_rcv[i] = hi_loc[i] - lo_loc[i] + width[i];
          }
        } else {
          plo_rcv[i] = 0;
          phi_rcv[i] = hi_loc[i] - lo_loc[i] + increment[i];
          plo_snd[i] = 0;
          phi_snd[i] = hi_loc[i] - lo_loc[i] + increment[i];
        }
      }

      /* Get pointer to local data buffer and remote data
         buffer as well as lists of leading dimenstions */
      gam_LocationWithGhosts(me, handle, plo_snd, &ptr_snd, ld_loc);
#if GHOST_PRINT
      printf("p[%d] 2: plo_snd[0]: %d plo_snd[1]: %d ptr_snd: %d\n",
          GAme, plo_snd[0], plo_snd[1], (Integer)ptr_snd);
      fflush(stdout);
#endif
      gam_LocationWithGhosts(me, handle, plo_rcv, &ptr_rcv, ld_loc);
#if GHOST_PRINT
      printf("p[%d] 2: plo_rcv[0]: %d plo_rcv[1]: %d ptr_rcv: %d\n",
          GAme, plo_rcv[0], plo_rcv[1], (Integer)ptr_rcv);
      fflush(stdout);
#endif

      /* Evaluate strides for send and recieve */
      gam_setstride(ndim, size, ld_loc, ld_loc, stride_rcv,
          stride_snd);

      /* Compute the number of elements in each dimension and store
         result in count. Scale the first element in count by the
         element size. */
      gam_ComputeCount(ndim, plo_rcv, phi_rcv, count);
      gam_CountElems(ndim, plo_snd, phi_snd, &length);
      length *= size;
      count[0] *= size;

      /* Need to reallocate memory if length > buflen */
      /* TO DO */

      /* Fill send buffer with data. */
#if GHOST_PRINT
      printf("p[%d]b 2: ptr_snd: %d ptr_rcv: %d\n", GAme, (Integer)ptr_snd,
              (Integer)ptr_rcv);
      printf("p[%d]b 2: length: %d proc_rem_snd: %d proc_rem_rcv: %d\n",
          GAme, (int)length, (int)proc_rem_snd, (int)proc_rem_rcv);
      printf("p[%d]b 2: count[0]: %d stride_rcv[0]: %d stride_rcv[1]: %d\n",
          GAme, count[0], stride_rcv[0],stride_rcv[1]);
      printf("p[%d]b 2: count[1]: %d stride_rcv[2]: %d stride_rcv[3]: %d\n",
          GAme, count[1], stride_rcv[2],stride_rcv[3]);
      printf("p[%d]b 2: count[2]: %d stride_snd[0]: %d stride_snd[1]: %d\n",
          GAme, count[2], stride_snd[0],stride_snd[1]);
      printf("p[%d]b 2: count[3]: %d stride_snd[2]: %d stride_snd[3]: %d\n",
          GAme, count[3], stride_snd[2],stride_snd[3]);
      printf("p[%d]b 2: snd_ptr: %d rcv_ptr: %d\n", GAme, (Integer)snd_ptr,
          (Integer)rcv_ptr);
      fflush(stdout);
#endif
      armci_write_strided(ptr_snd, (int)ndim-1, stride_snd, count, snd_ptr);

      /* Send Messages. If processor has odd index in direction idx, it
       * sends message first, if processor has even index it receives
       * message first. Then process is reversed. Also need to account
       * for whether or not there are an odd number of processors along
       * update direction. */

#if GHOST_PRINT
      printf("p[%d] 2: msgcnt: %d length: %d bufsize: %d proc_rem_snd: %d proc_rem_rcv: %d\n",
          GAme, msgcnt, length, bufsize, (int)proc_rem_snd, (int)proc_rem_rcv);
      fflush(stdout);
#endif
      snd_ptr = snd_ptr_orig;
      rcv_ptr = rcv_ptr_orig;
      if (GAme != proc_rem_rcv) {
        if (GA[handle].nblock[idx]%2 == 0) {
          if (index[idx]%2 != 0) {
            armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
          } else {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
          }
          if (index[idx]%2 != 0) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
          } else {
            armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
          }
        } else {
          pmax = GA[handle].nblock[idx] - 1;
          if (index[idx]%2 != 0) {
            armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
          } else if (index[idx] != 0) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
          }
          if (index[idx]%2 != 0) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
          } else if (index[idx] != pmax) {
            armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
          }
          /* make up for odd processor at end of string */
          if (index[idx] == pmax) {
            armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
          }
          if (index[idx] == 0) {
            armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
          }
        }
      } else {
        rcv_ptr = snd_ptr;
      }
      /* copy data back into global array */
      armci_read_strided(ptr_rcv, (int)ndim-1, stride_rcv, count, rcv_ptr);
      msgcnt++;
    }
    if (GA[handle].corner_flag)
      increment[idx] = 2*width[idx];
  }

  ga_free(rcv_ptr_orig);
  ga_free(snd_ptr_orig);
  GA_POP_NAME;
  return TRUE;
}

/* Utility function for ga_update5_ghosts routine */
static inline double waitforflags (int *ptr1, int *ptr2)
{
  int i = 1;
  double val = 0;
  while (*ptr1 ==  0 || *ptr2 == 0) {
    val = exp(-(double)i++);
  }
#if 0
  printf("%d: flags set at %p and %p\n",GAme,ptr1,ptr2); fflush(stdout);
#endif
  return(val);
}

#if 0
/* Stub in new ARMCI_PutS_flag call until actual implementation is
   available */
int ARMCI_PutS_flag__(
      void* src_ptr,        /* pointer to 1st segment at source */
      int src_stride_arr[], /* array of strides at source */
      void* dst_ptr,        /* pointer to 1st segment at destination */
      int dst_stride_arr[], /* array of strides at destination */
      int count[],          /* number of units at each stride level,
                               count[0] = #bytes */
      int stride_levels,    /* number of stride levels */
      int *flag,            /* pointer to remote flag */
      int *val,             /* pointer to value to set flag upon completion of
                               data transfer */
      int proc              /* remote process(or) ID */
      )
{
  int bytes;
  /* Put local data on remote processor */
  ARMCI_PutS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr,
             count, stride_levels, proc);

  /* Send signal to remote processor that data transfer has
   * been completed. */
  bytes = sizeof(int);
  ARMCI_Put(val, flag, bytes, proc);
  return 1;
}
#endif

/*\ UPDATE GHOST CELLS OF GLOBAL ARRAY USING SHIFT ALGORITHM AND PUT CALLS
 *  WITHOUT ANY BARRIERS
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_update55_ghosts = pnga_update55_ghosts
#endif
logical pnga_update55_ghosts(Integer g_a)
{
  Integer idx, i, np, handle=GA_OFFSET + g_a, proc_rem;
  Integer size, ndim, nwidth, increment[MAXDIM];
  Integer width[MAXDIM];
  Integer dims[MAXDIM];
  Integer lo_loc[MAXDIM], hi_loc[MAXDIM];
  Integer plo_loc[MAXDIM]/*, phi_loc[MAXDIM]*/;
  Integer tlo_rem[MAXDIM], thi_rem[MAXDIM];
  Integer slo_rem[MAXDIM], shi_rem[MAXDIM];
  Integer plo_rem[MAXDIM], phi_rem[MAXDIM];
  Integer ld_loc[MAXDIM], ld_rem[MAXDIM];
  int stride_loc[MAXDIM], stride_rem[MAXDIM],count[MAXDIM];
  int msgcnt;
  char *ptr_loc, *ptr_rem;
  Integer me = pnga_nodeid();
  Integer p_handle;

  /* This routine makes use of the shift algorithm to update data in the
   * ghost cells bounding the local block of visible data. The shift
   * algorithm starts by updating the blocks of data along the first
   * dimension by grabbing a block of data that is width[0] deep but
   * otherwise matches the  dimensions of the data residing on the
   * calling processor. The update of the second dimension, however,
   * grabs a block that is width[1] deep in the second dimension but is
   * ldim0 + 2*width[0] in the first dimensions where ldim0 is the
   * size of the visible data along the first dimension. The remaining
   * dimensions are left the same. For the next update, the width of the
   * second dimension is also increased by 2*width[1] and so on. This
   * algorith makes use of the fact that data for the dimensions that
   * have already been updated is available on each processor and can be
   * used in the updates of subsequent dimensions. The total number of
   * separate updates is 2*ndim, an update in the negative and positive
   * directions for each dimension.
   *
   * This operation is implemented using put calls to place the
   * appropriate data on remote processors. To signal the remote
   * processor that it has received the data, a second put call
   * consisting of a single integer is sent after the first put call and
   * used to update a signal buffer on the remote processor. Each
   * processor can determine how much data it has received by checking
   * its signal buffer. 
   *
   * To perform the update, this routine makes use of several copies of
   * indices marking the upper and lower limits of data. Indices
   * beginning with the character "p" are relative indices marking the
   * location of the data set relative to the origin the local patch of
   * the global array, all other indices are in absolute coordinates and
   * mark locations in the total global array. The indices used by this
   * routine are described below.
   *
   *       lo_loc[], hi_loc[]: The lower and upper indices of the visible
   *       block of data held by the calling processor.
   *
   *       lo_rem[], hi_rem[]: The lower and upper indices of the block
   *       of data on a remote processor or processors that is needed to
   *       fill in the calling processors ghost cells. These indices are
   *       NOT corrected for wrap-around (periodic) boundary conditions
   *       so they can be negative or greater than the array dimension
   *       values held in dims[].
   *
   *       slo_rem[], shi_rem[]: Similar to lo_rem[] and hi_rem[], except
   *       that these indices have been corrected for wrap-around
   *       boundary conditions. 
   *
   *       thi_rem[], thi_rem[]: The lower and upper indices of the visible
   *       data on a remote processor.
   *
   *       plo_loc[], phi_loc[]: The indices of the local data patch that
   *       is going to be updated.
   *
   *       plo_rem[], phi_rem[]: The indices of the data patch on the
   *       remote processor that will be used to update the data on the
   *       calling processor. Note that the dimensions of the patches
   *       represented by plo_loc[], plo_rem[] and plo_loc[], phi_loc[]
   *       must be the same.
   */

  /* if global array has no ghost cells, just return */
  if (!pnga_has_ghosts(g_a)) return TRUE;

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  p_handle = GA[handle].p_handle;

  /* initialize range increments and get array dimensions */
  for (idx=0; idx < ndim; idx++) {
    increment[idx] = 0;
    width[idx] = (Integer)GA[handle].width[idx];
    dims[idx] = (Integer)GA[handle].dims[idx];
    if (lo_loc[idx] == 0 && hi_loc[idx] == -1) return FALSE;
  }

  /* Check to make sure that global array is well-behaved (all processors
     have data and the width of the data in each dimension is greater
     than the corresponding value in width[]. */
  if (!gai_check_ghost_distr(g_a)) return FALSE;

  GA_PUSH_NAME("ga_update55_ghosts");

  /* Get pointer to local memory */
  ptr_loc = GA[handle].ptr[GAme];
  /* obtain range of data that is held by local processor */
  pnga_distribution(g_a,me,lo_loc,hi_loc);

  /* loop over dimensions for sequential update using shift algorithm */
  msgcnt = 0;
  (*GA_Update_Signal) = 1;
  for (idx=0; idx < ndim; idx++) {
    nwidth = width[idx];

    /* Do not bother with update if nwidth is zero */
    if (nwidth != 0) {

      /* Perform update in negative direction. */
      get_remote_block_neg(idx, ndim, lo_loc, hi_loc, slo_rem, shi_rem,
                           dims, width);
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rem, shi_rem, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rem, shi_rem, g_a);

      /* Get actual coordinates of desired location of remote
         data as well as the actual coordinates of the local chunk
         of data that will be sent to remote processor (these
         coordinates take into account the presence of ghost
         cells). Start by finding out what data is actually held by
         remote processor. */
      proc_rem = GA_proclist[0];
      pnga_distribution(g_a, proc_rem, tlo_rem, thi_rem);
      for (i = 0; i < ndim; i++) {
        if (increment[i] == 0) {
          if (i == idx) {
            plo_rem[i] = thi_rem[i] - tlo_rem[i] + width[i] + 1;
            phi_rem[i] = thi_rem[i] - tlo_rem[i] + 2*width[i];
            plo_loc[i] = width[i];
            /*phi_loc[i] = 2*width[i] - 1;*/
          } else {
            plo_rem[i] = width[i];
            phi_rem[i] = thi_rem[i] - tlo_rem[i] + width[i];
            plo_loc[i] = width[i];
            /*phi_loc[i] = hi_loc[i] - lo_loc[i] + width[i];*/
          }
        } else {
          plo_rem[i] = 0;
          phi_rem[i] = thi_rem[i] - tlo_rem[i] + increment[i];
          plo_loc[i] = 0;
          /*phi_loc[i] = hi_loc[i] - lo_loc[i] + increment[i];*/
        }
      }

      /* Get pointer to local data buffer and remote data
         buffer as well as lists of leading dimenstions */
      gam_LocationWithGhosts(me, handle, plo_loc, &ptr_loc, ld_loc);
      gam_LocationWithGhosts(proc_rem, handle, plo_rem, &ptr_rem, ld_rem);

      /* Evaluate strides on local and remote processors */
      gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
          stride_loc);

      /* Compute the number of elements in each dimension and store
         result in count. Scale the first element in count by the
         element size. */
      gam_ComputeCount(ndim, plo_rem, phi_rem, count);
      count[0] *= size;

      /* Put local data on remote processor */
      if (p_handle >= 0) {
        proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem];
      }
#if 0
      ARMCI_PutS(ptr_loc, stride_loc, ptr_rem, stride_rem, count, ndim- 1, proc_rem);
      /* Send signal to remote processor that data transfer has been completed. */
      bytes = sizeof(int);
      ARMCI_Put(GA_Update_Signal, GA_Update_Flags[proc_rem]+msgcnt, bytes, proc_rem);
#else
      ARMCI_PutS_flag(ptr_loc, stride_loc, ptr_rem, stride_rem, count,
          (int)(ndim - 1), GA_Update_Flags[proc_rem]+msgcnt,
          *GA_Update_Signal, (int)proc_rem);
#endif
      msgcnt++;

      /* Perform update in positive direction. */
      get_remote_block_pos(idx, ndim, lo_loc, hi_loc, slo_rem, shi_rem,
                           dims, width);
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rem, shi_rem, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rem, shi_rem, g_a);

      /* Get actual coordinates of desired chunk of remote
         data as well as the actual coordinates of the local chunk
         of data that will receive the remote data (these
         coordinates take into account the presence of ghost
         cells). Start by finding out what data is actually held by
         remote processor. */
      proc_rem = GA_proclist[0];
      pnga_distribution(g_a, proc_rem, tlo_rem, thi_rem);
      for (i = 0; i < ndim; i++) {
        if (increment[i] == 0) {
          if (i == idx) {
            plo_rem[i] = 0;
            phi_rem[i] = width[i] - 1;
            plo_loc[i] = hi_loc[i] - lo_loc[i] + width[i] - 1;
            /*phi_loc[i] = hi_loc[i] - lo_loc[i] + 2*width[i] - 1;*/
          } else {
            plo_rem[i] = width[i];
            phi_rem[i] = thi_rem[i] - tlo_rem[i] + width[i];
            plo_loc[i] = width[i];
            /*phi_loc[i] = hi_loc[i] - lo_loc[i] + width[i];*/
          }
        } else {
          plo_rem[i] = 0;
          phi_rem[i] = thi_rem[i] - tlo_rem[i] + increment[i];
          plo_loc[i] = 0;
          /*phi_loc[i] = hi_loc[i] - lo_loc[i] + increment[i];*/
        }
      }

      /* Get pointer to local data buffer and remote data
         buffer as well as lists of leading dimenstions */
      gam_LocationWithGhosts(GAme, handle, plo_loc, &ptr_loc, ld_loc);
      gam_LocationWithGhosts(proc_rem, handle, plo_rem, &ptr_rem, ld_rem);

      /* Evaluate strides on local and remote processors */
      gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
          stride_loc);

      /* Compute the number of elements in each dimension and store
         result in count. Scale the first element in count by the
         element size. */
      gam_ComputeCount(ndim, plo_rem, phi_rem, count);
      count[0] *= size;

      /* Put local data on remote processor */
      if (p_handle >= 0) {
        proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem];
      }
#if 0
      ARMCI_PutS(ptr_loc, stride_loc, ptr_rem, stride_rem, count, ndim- 1, proc_rem);
      /* Send signal to remote processor that data transfer has been completed. */
      bytes = sizeof(int);
      ARMCI_Put(GA_Update_Signal, GA_Update_Flags[proc_rem]+msgcnt, bytes, proc_rem);

#else
      ARMCI_PutS_flag(ptr_loc, stride_loc, ptr_rem, stride_rem, count,
          (int)(ndim - 1), GA_Update_Flags[proc_rem]+msgcnt,
          *GA_Update_Signal, (int)proc_rem);
#endif
      msgcnt++;
    }
    /* check to make sure that all messages have been recieved before
       starting update along new dimension */
    waitforflags((GA_Update_Flags[GAme]+msgcnt-2),
        (GA_Update_Flags[GAme]+msgcnt-1));
    /* update increment array */
    increment[idx] = 2*nwidth;
  }

  /* set GA_Update_Flags array to zero for next update operation. */
  for (idx=0; idx < 2*ndim; idx++) {
    GA_Update_Flags[GAme][idx] = 0;
  }

  GA_POP_NAME;
  return TRUE;
}

/*\ UPDATE GHOST CELLS OF GLOBAL ARRAY USING NON-BLOCKING GET CALLS AND RETURN
 *  A NON-BLOCKING HANDLE
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_update_ghosts_nb = pnga_update_ghosts_nb
#endif
void pnga_update_ghosts_nb(Integer g_a, Integer *nbhandle)
{
  Integer idx, ipx, np, handle=GA_OFFSET + g_a, proc_rem;
  Integer ntot, mask[MAXDIM];
  Integer size, ndim, i, itmp;
  Integer width[MAXDIM], dims[MAXDIM];
  Integer lo_loc[MAXDIM], hi_loc[MAXDIM];
  /*Integer tlo_loc[MAXDIM], thi_loc[MAXDIM];*/
  Integer plo_loc[MAXDIM], phi_loc[MAXDIM];
  Integer tlo_rem[MAXDIM], thi_rem[MAXDIM];
  Integer plo_rem[MAXDIM];
  Integer ld_loc[MAXDIM], ld_rem[MAXDIM];
  logical mask0;
  int stride_loc[MAXDIM], stride_rem[MAXDIM],count[MAXDIM];
  char *ptr_loc, *ptr_rem;
  Integer me = pnga_nodeid();
  Integer p_handle;

  /* if global array has no ghost cells, just return */
  if (!pnga_has_ghosts(g_a)) {
    return;
  }

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  p_handle = GA[handle].p_handle;
  /* initialize ghost cell widths and get array dimensions */
  for (idx=0; idx < ndim; idx++) {
    width[idx] = (Integer)GA[handle].width[idx];
    dims[idx] = (Integer)GA[handle].dims[idx];
  }

  /* Check to make sure that global array is well-behaved (all processors
     have data and the width of the data in each dimension is greater than
     the corresponding value in width[]). */
  if (!gai_check_ghost_distr(g_a)) return;

  /* Create non-blocking handle */
  ga_init_nbhandle(nbhandle);

  GA_PUSH_NAME("ga_update_ghosts_nb");
  /* Get pointer to local memory */
  ptr_loc = GA[handle].ptr[me];
  /* obtain range of data that is held by local processor */
  pnga_distribution(g_a,me,lo_loc,hi_loc);

  /* evaluate total number of PUT operations that will be required */
  ntot = 1;
  for (idx=0; idx < ndim; idx++) ntot *= 3;

  /* Loop over all GET operations. The operation corresponding to the
     mask of all zeros is left out. */
  for (ipx=0; ipx < ntot; ipx++) {
    /* Convert ipx to corresponding mask values */
    itmp = ipx;
    mask0 = TRUE;
    for (idx = 0; idx < ndim; idx++) {
      i = itmp%3;
      mask[idx] = i-1;
      if (mask[idx] != 0) mask0 = FALSE;
      itmp = (itmp-i)/3;
    }
    if (mask0) continue;

    /* check to see if ghost cell block has zero elements*/
    mask0 = FALSE;
    itmp = 0;
    for (idx = 0; idx < ndim; idx++) {
      if (mask[idx] != 0 && width[idx] == 0) mask0 = TRUE;
      if (mask[idx] != 0) itmp++;
    }
    if (mask0) continue;
    /* Now that mask has been determined, find data that is to be moved
     * and identify processor to which it is going. Wrap boundaries
     * around, if necessary */
    for (idx = 0; idx < ndim; idx++) {
      if (mask[idx] == 0) {
        tlo_rem[idx] = lo_loc[idx];
        thi_rem[idx] = hi_loc[idx];
      } else if (mask[idx] == -1) {
        if (lo_loc[idx] > 1) {
          tlo_rem[idx] = lo_loc[idx]-width[idx];
          thi_rem[idx] = lo_loc[idx]-1;
        } else {
          tlo_rem[idx] = dims[idx]-width[idx]+1;
          thi_rem[idx] = dims[idx];
        }
      } else if (mask[idx] == 1) {
        if (hi_loc[idx] < dims[idx]) {
          tlo_rem[idx] = hi_loc[idx] + 1;
          thi_rem[idx] = hi_loc[idx] + width[idx];
        } else {
          tlo_rem[idx] = 1;
          thi_rem[idx] = width[idx];
        }
      } else {
        fprintf(stderr,"Illegal mask value found\n");
      }
    }
    /* Locate remote processor from which data must be retrieved */
    if (!pnga_locate_region(g_a, tlo_rem, thi_rem, _ga_map,
       GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
       tlo_rem, thi_rem, g_a);
    if (np > 1) {
      fprintf(stderr,"More than one remote processor found\n");
    }
    /* Remote processor has been identified, now get ready to get
       data from it. Start by getting distribution on remote
       processor.*/
    proc_rem = GA_proclist[0];
    pnga_distribution(g_a, proc_rem, tlo_rem, thi_rem);
    for (idx = 0; idx < ndim; idx++) {
      if (mask[idx] == 0) {
        plo_loc[idx] = width[idx];
        phi_loc[idx] = hi_loc[idx]-lo_loc[idx]+width[idx];
        plo_rem[idx] = plo_loc[idx];
      } else if (mask[idx] == -1) {
        plo_loc[idx] = 0;
        phi_loc[idx] = width[idx]-1;
        plo_rem[idx] = thi_rem[idx]-tlo_rem[idx]+1;
      } else if (mask[idx] == 1) {
        plo_loc[idx] = hi_loc[idx]-lo_loc[idx]+width[idx]+1;
        phi_loc[idx] = hi_loc[idx]-lo_loc[idx]+2*width[idx];
        plo_rem[idx] = width[idx];
      }
    }
    /* Get pointer to local data buffer and remote data
       buffer as well as lists of leading dimenstions */
    gam_LocationWithGhosts(me, handle, plo_loc, &ptr_loc, ld_loc);
    gam_LocationWithGhosts(proc_rem, handle, plo_rem, &ptr_rem, ld_rem);

    /* Evaluate strides on local and remote processors */
    gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
                  stride_loc);

    /* Compute the number of elements in each dimension and store
       result in count. Scale the first element in count by the
       element size. */
    gam_ComputeCount(ndim, plo_loc, phi_loc, count);
    count[0] *= size;
 
    /* get data from remote processor */
    if (p_handle >= 0) {
      proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem];
    }
    ARMCI_NbGetS(ptr_rem, stride_rem, ptr_loc, stride_loc, count,
        (int)(ndim - 1), (int)proc_rem, 
        (armci_hdl_t*)get_armci_nbhandle(nbhandle));
  }

  GA_POP_NAME;
  return;
}

/*\ UPDATE GHOST CELLS OF GLOBAL ARRAY ALONG ONE SIDE OF ARRAY
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_update_ghost_dir = pnga_update_ghost_dir
#endif
logical pnga_update_ghost_dir(Integer g_a,    /* GA handle */
                                   Integer pdim,   /* Dimension of update */
                                   Integer pdir,   /* Direction of update (+/-1) */
                                   logical pflag)  /* include corner cells */
{
  Integer idx, ipx, inx, np, handle=GA_OFFSET + g_a, proc_rem;
  Integer ntot, mask[MAXDIM],lmask[MAXDIM];
  Integer size, ndim, i, itmp, idim, idir;
  Integer width[MAXDIM], dims[MAXDIM];
  Integer lo_loc[MAXDIM], hi_loc[MAXDIM];
  Integer plo_loc[MAXDIM], phi_loc[MAXDIM];
  Integer tlo_rem[MAXDIM], thi_rem[MAXDIM];
  Integer plo_rem[MAXDIM]/*, phi_rem[MAXDIM]*/;
  Integer ld_loc[MAXDIM], ld_rem[MAXDIM];
  logical flag;
  int stride_loc[MAXDIM], stride_rem[MAXDIM],count[MAXDIM];
  char *ptr_loc, *ptr_rem;
  Integer me = pnga_nodeid();
  Integer p_handle;

  int local_sync_begin,local_sync_end;

  local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
  _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/

  /* if global array has no ghost cells, just return */
  if (!pnga_has_ghosts(g_a)) 
    return TRUE;
  
  p_handle = GA[handle].p_handle;
  if(local_sync_begin)pnga_pgroup_sync(p_handle);
  idim = pdim;
  idir = pdir;
  flag = pflag;

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  /* initialize ghost cell widths and get array dimensions */
  for (idx=0; idx < ndim; idx++) {
    width[idx] = (Integer)GA[handle].width[idx];
    dims[idx] = (Integer)GA[handle].dims[idx];
  }

  /* Check to make sure that global array is well-behaved (all processors
     have data and the width of the data in each dimension is greater than
     the corresponding value in width[]). */
  ipx = 0;
  for (idx = 0; idx < ndim; idx++) {
    for (np = 0; np < GA[handle].nblock[idx]; np++) {
      if (np < GA[handle].nblock[idx] - 1) {
        if (GA[handle].mapc[ipx+1]-GA[handle].mapc[ipx]+1<width[idx]) {
          return FALSE;
        }
      } else {
        if (GA[handle].dims[idx]-GA[handle].mapc[ipx]+1<width[idx]) {
          return FALSE;
        }
      }
      ipx++;
    }
  }

  GA_PUSH_NAME("nga_update_ghost_dir");
  /* Get pointer to local memory */
  ptr_loc = GA[handle].ptr[GAme];
  /* obtain range of data that is held by local processor */
  pnga_distribution(g_a,me,lo_loc,hi_loc);

  /* evaluate total number of GET operations */
  ntot = 1;
  if (flag) {
    for (idx=0; idx < ndim-1; idx++) ntot *= 3;
  }

  /* Loop over all GET operations. */
  for (ipx=0; ipx < ntot; ipx++) {
    /* Convert ipx to corresponding mask values */
    if (flag) {
      itmp = ipx;
      for (idx = 0; idx < ndim-1; idx++) {
        i = itmp%3;
        lmask[idx] = i-1;
        itmp = (itmp-i)/3;
      }
    } else {
      for (idx = 0; idx < ndim-1; idx++) lmask[idx] = 0;
    }
    inx = 0;
    for (idx = 0; idx < ndim; idx++) {
      if (idx == idim-1) {
        mask[idx] = idir;
      } else {
        mask[idx] = lmask[inx];
        inx++;
      }
    }
    /* Now that mask has been determined, find processor that contains
     * data needed by the corresponding block of ghost cells */
    for (idx = 0; idx < ndim; idx++) {
      if (mask[idx] == 0) {
        tlo_rem[idx] = lo_loc[idx];
        thi_rem[idx] = hi_loc[idx];
      } else if (mask[idx] == -1) {
        if (lo_loc[idx] > 1) {
          tlo_rem[idx] = lo_loc[idx]-width[idx];
          thi_rem[idx] = lo_loc[idx]-1;
        } else {
          tlo_rem[idx] = dims[idx]-width[idx]+1;
          thi_rem[idx] = dims[idx];
        }
      } else if (mask[idx] == 1) {
        if (hi_loc[idx] < dims[idx]) {
          tlo_rem[idx] = hi_loc[idx] + 1;
          thi_rem[idx] = hi_loc[idx] + width[idx];
        } else {
          tlo_rem[idx] = 1;
          thi_rem[idx] = width[idx];
        }
      } else {
        fprintf(stderr,"Illegal mask value found\n");
      }
    }
    /* Locate remote processor to which data must be sent */
    if (!pnga_locate_region(g_a, tlo_rem, thi_rem, _ga_map,
       GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
       tlo_rem, thi_rem, g_a);
    if (np > 1) {
      fprintf(stderr,"More than one remote processor found\n");
    }
    /* Remote processor has been identified, now get ready to get
       data from it. Start by getting distribution on remote
       processor.*/
    proc_rem = GA_proclist[0];
    pnga_distribution(g_a, proc_rem, tlo_rem, thi_rem);
    for (idx = 0; idx < ndim; idx++) {
      if (mask[idx] == 0) {
        plo_loc[idx] = width[idx];
        phi_loc[idx] = hi_loc[idx]-lo_loc[idx]+width[idx];
        plo_rem[idx] = plo_loc[idx];
        /*phi_rem[idx] = phi_loc[idx];*/
      } else if (mask[idx] == -1) {
        plo_loc[idx] = 0;
        phi_loc[idx] = width[idx]-1;
        plo_rem[idx] = thi_rem[idx]-tlo_rem[idx]+1;
        /*phi_rem[idx] = thi_rem[idx]-tlo_rem[idx]+width[idx];*/
      } else if (mask[idx] == 1) {
        plo_loc[idx] = hi_loc[idx]-lo_loc[idx]+width[idx]+1;
        phi_loc[idx] = hi_loc[idx]-lo_loc[idx]+2*width[idx];
        plo_rem[idx] = width[idx];
        /*phi_rem[idx] = 2*width[idx]-1;*/
      }
    }
    /* Get pointer to local data buffer and remote data
       buffer as well as lists of leading dimenstions */
    gam_LocationWithGhosts(me, handle, plo_loc, &ptr_loc, ld_loc);
    gam_LocationWithGhosts(proc_rem, handle, plo_rem, &ptr_rem, ld_rem);

    /* Evaluate strides on local and remote processors */
    gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
                  stride_loc);

    /* Compute the number of elements in each dimension and store
       result in count. Scale the first element in count by the
       element size. */
    gam_ComputeCount(ndim, plo_loc, phi_loc, count);
    count[0] *= size;
 
    /* get data from remote processor */
    if (p_handle >= 0) {
      proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem];
    }
    ARMCI_GetS(ptr_rem, stride_rem, ptr_loc, stride_loc, count,
          (int)(ndim - 1), (int)proc_rem);
  }

  GA_POP_NAME;
  if(local_sync_end)pnga_pgroup_sync(p_handle);
  return TRUE;
}

/*uncomment for using message passing sendrecv in north south direction */
/*#define USE_MP_NORTHSOUTH */


/*\ UPDATE GHOST CELLS OF GLOBAL ARRAY USING PUT CALLS WITHOUT CORNERS AND
 *  WITHOUT ANY BARRIERS
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_update5_ghosts = pnga_update5_ghosts
#endif
logical pnga_update5_ghosts(Integer g_a)
{
  Integer idx, i, handle=GA_OFFSET + g_a;
  Integer /*size,*/ ndim, nwidth;
  Integer width[MAXDIM];
  Integer* proc_rem_ptr;
  int *stride_loc, *stride_rem,*count;
  int msgcnt, corner_flag, proc_rem;
  /* int bytes; */
  char *ptr_loc, *ptr_rem,*cache;
  int local_sync_begin,local_sync_end;
  Integer p_handle;
#ifdef USE_MP_NORTHSOUTH
  char send_name[32], rcv_name[32];
  void *snd_ptr, *rcv_ptr;
#endif
  /* This routine makes use of the shift algorithm to update data in the
   * ghost cells bounding the local block of visible data. The shift
   * algorithm starts by updating the blocks of data along the first
   * dimension by grabbing a block of data that is width[0] deep but
   * otherwise matches the  dimensions of the data residing on the
   * calling processor. The update of the second dimension, however,
   * grabs a block that is width[1] deep in the second dimension but is
   * ldim0 + 2*width[0] in the first dimensions where ldim0 is the
   * size of the visible data along the first dimension. The remaining
   * dimensions are left the same. For the next update, the width of the
   * second dimension is also increased by 2*width[1] and so on. This
   * algorith makes use of the fact that data for the dimensions that
   * have already been updated is available on each processor and can be
   * used in the updates of subsequent dimensions. The total number of
   * separate updates is 2*ndim, an update in the negative and positive
   * directions for each dimension.
   *
   * This operation is implemented using put calls to place the
   * appropriate data on remote processors. To signal the remote
   * processor that it has received the data, a second put call
   * consisting of a single integer is sent after the first put call and
   * used to update a signal buffer on the remote processor. Each
   * processor can determine how much data it has received by checking
   * its signal buffer. 
   */

  local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
  _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
  p_handle = GA[handle].p_handle;
  if(local_sync_begin)pnga_pgroup_sync(p_handle);

#ifdef USE_MP_NORTHSOUTH
  strcpy(send_name,"send_buffer");
  strcpy(rcv_name,"receive_buffer");

  snd_ptr = ga_malloc(buflen, GA[handle].type, send_name);
  rcv_ptr = ga_malloc(buflen, GA[handle].type, rcv_name);
#endif
 
  cache = (char *)GA[handle].cache;
  /* if global array has no ghost cells, just return */
  if (!pnga_has_ghosts(g_a)) return TRUE;

  /*size = GA[handle].elemsize;*/
  ndim = GA[handle].ndim;
  for (i=0; i<ndim; i++) {
    width[i] = (Integer)GA[handle].width[i];
  }

  if (!gai_check_ghost_distr(g_a)) return FALSE;

  GA_PUSH_NAME("pnga_update5_ghosts");

  /* loop over dimensions for sequential update using shift algorithm */
  msgcnt = 0;
  corner_flag = GA[handle].corner_flag;
  (*GA_Update_Signal) = 1;
  for (idx=0; idx < ndim; idx++) {
    nwidth = width[idx];
    if (nwidth != 0) {

      /* Perform update in negative direction. */
      ptr_rem = *(char **)(cache);
      if(ptr_rem==NULL) return FALSE;
      ptr_loc = *(char **)(cache+sizeof(char *));
      stride_loc = (int *)(cache+2*sizeof(char *));
      stride_rem = (int *)(stride_loc+ndim);
      count = (int *)(stride_rem+ndim);
      proc_rem_ptr = (Integer *)(count+ndim);
      proc_rem = (int)(*proc_rem_ptr);
      cache = (char *)(proc_rem_ptr+1);
          
      if (p_handle >= 0) {
        proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem];
      }
      if(count[0]>1000000){
        /*tries to use armci direct put when possible */
        ARMCI_PutS_flag(ptr_loc, stride_loc, ptr_rem, stride_rem, count,
            (int)(ndim - 1), GA_Update_Flags[proc_rem]+msgcnt,
            *GA_Update_Signal, proc_rem);
      }
      else{
#ifndef USE_MP_NORTHSOUTH
        ARMCI_PutS_flag(ptr_loc, stride_loc, ptr_rem, stride_rem, count,
            (int)(ndim - 1), GA_Update_Flags[proc_rem]+msgcnt,
            *GA_Update_Signal, proc_rem);
#else
#endif
      }

      msgcnt++;

      /* Perform update in positive direction. */
      ptr_rem = *(char **)(cache);
      ptr_loc = *(char **)(cache+sizeof(char *));
      stride_loc = (int *)(cache+2*sizeof(char *));
      stride_rem = (int *)(stride_loc+ndim);
      count = (int *)(stride_rem+ndim);
      proc_rem_ptr = (Integer *)(count+ndim);
      proc_rem = (int)(*proc_rem_ptr);
      cache = (char *)(proc_rem_ptr+1);

      if (p_handle >= 0) {
        proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem];
      }
      if(count[0]>1000000){
        /*tries to use armci direct put when possible */
        ARMCI_PutS_flag(ptr_loc, stride_loc, ptr_rem, stride_rem, count,
            (int)(ndim - 1), GA_Update_Flags[proc_rem]+msgcnt,
            *GA_Update_Signal, proc_rem);
      }
      else{
#ifndef USE_MP_NORTHSOUTH
        ARMCI_PutS_flag(ptr_loc, stride_loc, ptr_rem, stride_rem, count,
            (int)(ndim - 1), GA_Update_Flags[proc_rem]+msgcnt,
            *GA_Update_Signal, proc_rem);
#else
#endif

      }

      msgcnt++;

      if (corner_flag){
        /* check to make sure that last two messages have been recieved
           before starting update along a new dimension */
        waitforflags((GA_Update_Flags[GAme]+msgcnt-2),
          (GA_Update_Flags[GAme]+msgcnt-1));
        GA_Update_Flags[GAme][msgcnt-1]=0;
        GA_Update_Flags[GAme][msgcnt-2]=0;
      }
    }
  }
#if 1
  if (!corner_flag) {
    /* check to make sure that all messages have been recieved */
    while(msgcnt){
      waitforflags((GA_Update_Flags[GAme]+msgcnt-1),
          (GA_Update_Flags[GAme]+msgcnt-2));
      GA_Update_Flags[GAme][msgcnt-1]=0;
      GA_Update_Flags[GAme][msgcnt-2]=0;
      msgcnt-=2;
    }
  }
#endif 
  GA_POP_NAME;
  if(local_sync_end)pnga_pgroup_sync(p_handle);
  return TRUE;
}

/*#define UPDATE_SAMENODE_GHOSTS_FIRST*/

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_set_update5_info = pnga_set_update5_info
#endif
logical pnga_set_update5_info(Integer g_a)
{
  int i;
  Integer *proc_rem;
  Integer size, ndim, nwidth, increment[MAXDIM],np;
  Integer width[MAXDIM];
  Integer dims[MAXDIM];
  Integer lo_loc[MAXDIM], hi_loc[MAXDIM];
  Integer plo_loc[MAXDIM]/*, phi_loc[MAXDIM]*/;
  Integer tlo_rem[MAXDIM], thi_rem[MAXDIM];
  Integer slo_rem[MAXDIM], shi_rem[MAXDIM];
  Integer plo_rem[MAXDIM], phi_rem[MAXDIM];
  Integer ld_loc[MAXDIM], ld_rem[MAXDIM];
  int *stride_loc, *stride_rem,*count;
  int idx, corner_flag;
  char **ptr_loc, **ptr_rem,*cache;
  Integer handle = GA_OFFSET + g_a;
  int cache_size;
#ifdef UPDATE_SAMENODE_GHOSTS_FIRST
  int scope;
#endif
  Integer me = pnga_nodeid();
  Integer p_handle;

  /* This routine sets up the arrays that are used to transfer data
   * using the update5 algorithm. The arrays begining with the character
   * "p" represent relative indices marking the location of the data set
   * relative to the origin the local patch of the global array, all
   * other indices are in absolute coordinates and mark locations in the
   * total global array. The indices used by this routine are described
   * below.
   *
   *       lo_loc[], hi_loc[]: The lower and upper indices of the visible
   *       block of data held by the calling processor.
   *
   *       lo_rem[], hi_rem[]: The lower and upper indices of the block
   *       of data on a remote processor or processors that is needed to
   *       fill in the calling processors ghost cells. These indices are
   *       NOT corrected for wrap-around (periodic) boundary conditions
   *       so they can be negative or greater than the array dimension
   *       values held in dims[].
   *
   *       slo_rem[], shi_rem[]: Similar to lo_rem[] and hi_rem[], except
   *       that these indices have been corrected for wrap-around
   *       boundary conditions. 
   *
   *       thi_rem[], thi_rem[]: The lower and upper indices of the visible
   *       data on a remote processor.
   *
   *       plo_loc[], phi_loc[]: The indices of the local data patch that
   *       is going to be updated.
   *
   *       plo_rem[], phi_rem[]: The indices of the data patch on the
   *       remote processor that will be used to update the data on the
   *       calling processor. Note that the dimensions of the patches
   *       represented by plo_loc[], plo_rem[] and plo_loc[], phi_loc[]
   *       must be the same.
   */

  /* if global array has no ghost cells, just return */
  if (!pnga_has_ghosts(g_a)) return TRUE;

  /* Check to make sure that global array is well-behaved (all processors
     have data and the width of the data in each dimension is greater
     than the corresponding value in width[]. */
  if (!gai_check_ghost_distr(g_a)) return FALSE;

  ndim = GA[handle].ndim;
  p_handle = GA[handle].p_handle;
  size = GA[handle].elemsize;
  cache_size = 2*sizeof(char *)+3*sizeof(int)+sizeof(Integer);
  cache_size = 2*ndim*((cache_size/sizeof(double)) + 1);
  GA[handle].cache = (double *)malloc(sizeof(double)*cache_size);
  cache = (char *)GA[handle].cache;
  corner_flag = GA[handle].corner_flag;

  pnga_distribution(g_a,me,lo_loc,hi_loc); 
  for (idx=0; idx < ndim; idx++) {
    increment[idx] = 0;
    width[idx] = (Integer)GA[handle].width[idx];
    dims[idx] = (Integer)GA[handle].dims[idx];
    if (lo_loc[idx] == 0 && hi_loc[idx] == -1){
      *(char **)cache = NULL; 
      return FALSE;
    }
  } 
#ifdef UPDATE_SAMENODE_GHOSTS_FIRST
  for(scope=0;scope < 2; scope ++)
#endif
    for (idx=0; idx < ndim; idx++) {
      nwidth = width[idx];
      if (nwidth != 0) {  
      
        ptr_rem = (char **)cache;
        ptr_loc = (char **)(cache+sizeof(char *));
        stride_loc = (int *)(cache+2*sizeof(char *));
        stride_rem = (int *)(stride_loc+ndim);
        count = (int *)(stride_rem+ndim);
        proc_rem = (Integer *)(count+ndim);

        get_remote_block_neg(idx, ndim, lo_loc, hi_loc, slo_rem, shi_rem,
                             dims, width);
        if (!pnga_locate_region(g_a, slo_rem, shi_rem, _ga_map,
            GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
            slo_rem, shi_rem, g_a);

        *proc_rem = (Integer)GA_proclist[0];
        if (p_handle >= 0) {
          *proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[*proc_rem];
        }

#ifdef UPDATE_SAMENODE_GHOSTS_FIRST
        if(scope == 0 && ARMCI_Same_node(*proc_rem))
          goto do_negative;
#endif

        cache = (char *)(proc_rem+1);

        pnga_distribution(g_a, *proc_rem, tlo_rem, thi_rem);
        

        for (i = 0; i < ndim; i++) {
          if (increment[i] == 0) {
            if (i == idx) {
              plo_rem[i] = thi_rem[i] - tlo_rem[i] + width[i] + 1;
              phi_rem[i] = thi_rem[i] - tlo_rem[i] + 2*width[i];
              plo_loc[i] = width[i];
              /*phi_loc[i] = 2*width[i] - 1;*/
            } else {
              plo_rem[i] = width[i];
              phi_rem[i] = thi_rem[i] - tlo_rem[i] + width[i];
              plo_loc[i] = width[i];
              /*phi_loc[i] = hi_loc[i] - lo_loc[i] + width[i];*/
            }
          } else {
            plo_rem[i] = 0;
            phi_rem[i] = thi_rem[i] - tlo_rem[i] + increment[i];
            plo_loc[i] = 0;
            /*phi_loc[i] = hi_loc[i] - lo_loc[i] + increment[i];*/
          }
        }
        gam_LocationWithGhosts(me, handle, plo_loc, ptr_loc, ld_loc);
        gam_LocationWithGhosts(*proc_rem, handle, plo_rem, ptr_rem, ld_rem);
 
        /* Evaluate strides on local and remote processors */
        gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
            stride_loc);
        gam_ComputeCount(ndim, plo_rem, phi_rem, count);
        count[0] *= size;
        if (p_handle >= 0) {
          *proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[*proc_rem];
        }

#ifdef UPDATE_SAMENODE_GHOSTS_FIRST
        do_negative:
#endif

       /*BJP proc_rem++; */
        ptr_rem = (char **)cache;
        ptr_loc = (char **)(cache+sizeof(char *));
        stride_loc = (int *)(cache+2*sizeof(char *));
        stride_rem = (int *)(stride_loc+ndim);
        count = (int *)(stride_rem+ndim);
        proc_rem = (Integer *)(count+ndim);

        get_remote_block_pos(idx, ndim, lo_loc, hi_loc, slo_rem, shi_rem,
                             dims, width);

        if (!pnga_locate_region(g_a, slo_rem, shi_rem, _ga_map,
            GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
            slo_rem, shi_rem, g_a);

        *proc_rem = (Integer)GA_proclist[0];
        if (p_handle >= 0) {
          *proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[*proc_rem];
        }

#ifdef UPDATE_SAMENODE_GHOSTS_FIRST
        if(scope == 0 && ARMCI_Same_node(*proc_rem))
          continue;
#endif

        cache = (char *)(proc_rem+1);

        pnga_distribution(g_a, *proc_rem, tlo_rem, thi_rem);



        for (i = 0; i < ndim; i++) {
          if (increment[i] == 0) {
            if (i == idx) {
              plo_rem[i] = 0;
              phi_rem[i] = width[i] - 1;
              plo_loc[i] = hi_loc[i] - lo_loc[i] + width[i] - 1;
              /*phi_loc[i] = hi_loc[i] - lo_loc[i] + 2*width[i] - 1;*/
            } else {
              plo_rem[i] = width[i];
              phi_rem[i] = thi_rem[i] - tlo_rem[i] + width[i];
              plo_loc[i] = width[i];
              /*phi_loc[i] = hi_loc[i] - lo_loc[i] + width[i];*/
            }
          } else {
            plo_rem[i] = 0;
            phi_rem[i] = thi_rem[i] - tlo_rem[i] + increment[i];
            plo_loc[i] = 0;
            /*phi_loc[i] = hi_loc[i] - lo_loc[i] + increment[i];*/
          }
        }


        gam_LocationWithGhosts(GAme, handle, plo_loc, ptr_loc, ld_loc);
        gam_LocationWithGhosts(*proc_rem, handle, plo_rem, ptr_rem, ld_rem);

        gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
            stride_loc);

        gam_ComputeCount(ndim, plo_rem, phi_rem, count);
        count[0] *= size;
        if (p_handle >= 0) {
          *proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[*proc_rem];
        }

        if (corner_flag)
          increment[idx] = 2*nwidth;
      }
    }
    return TRUE;
}

/*\ UPDATE GHOST CELLS OF GLOBAL ARRAY USING SHIFT ALGORITHM
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_update_ghosts = pnga_update_ghosts
#endif
void pnga_update_ghosts(Integer g_a)
{
  /* Wrapper program for ghost cell update operations. If optimized
     update operation fails then use slow but robust version of
     update operation */
   int local_sync_begin,local_sync_end;
   Integer handle = GA_OFFSET + g_a;

   local_sync_begin = _ga_sync_begin; local_sync_end = _ga_sync_end;
   _ga_sync_begin = 1; _ga_sync_end=1; /*remove any previous masking*/
   if(local_sync_begin)pnga_pgroup_sync(GA[handle].p_handle);

#ifdef CRAY_T3D
   if (!pnga_update5_ghosts(g_a))
#else
   if (!pnga_update4_ghosts(g_a))
#endif
   {
     pnga_update1_ghosts(g_a);
   }

   if(local_sync_end)pnga_pgroup_sync(GA[handle].p_handle);
}

/* Utility function for ga_update6_ghosts routine */
static double waitformixedflags (int flag1, int flag2, int *ptr1, int *ptr2) {
  int i = 1;
  double val = 0;
  while ((flag1 && *ptr1 ==  0) ||
         (flag2 && *ptr2 == 0)) {
    val = exp(-(double)i++);
  }
  return(val);
}

/*\ UPDATE GHOST CELLS OF GLOBAL ARRAY USING SHIFT ALGORITHM AND
 *  MESSAGE PASSING
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_update6_ghosts = pnga_update6_ghosts
#endif
logical pnga_update6_ghosts(Integer g_a)
{
  Integer idx, idir, i, np, handle=GA_OFFSET + g_a;
  Integer size, buflen, buftot, bufsize, ndim, increment[MAXDIM];
  Integer proc_rem_snd, proc_rem_rcv, pmax;
  Integer msgcnt, length;
  Integer width[MAXDIM], dims[MAXDIM], index[MAXDIM];
  Integer lo_loc[MAXDIM], hi_loc[MAXDIM];
  Integer plo_rem[MAXDIM]/*, phi_rem[MAXDIM]*/;
  Integer tlo_rem[MAXDIM], thi_rem[MAXDIM];
  Integer plo_snd[MAXDIM], phi_snd[MAXDIM];
  Integer lo_rcv[MAXDIM], hi_rcv[MAXDIM];
  Integer slo_rcv[MAXDIM], shi_rcv[MAXDIM];
  Integer plo_rcv[MAXDIM], phi_rcv[MAXDIM];
  Integer ld_loc[MAXDIM], ld_rem[MAXDIM];
  int msglen;
  int stride_snd[MAXDIM], stride_rcv[MAXDIM],count[MAXDIM];
  int stride_rem[MAXDIM];
  int flag1=0, flag2=0, sprocflag, rprocflag;
  char *ptr_snd, *ptr_rcv;
  char /* *ptr_loc,*/ *ptr_rem;
  char send_name[32], rcv_name[32];
  void *snd_ptr, *rcv_ptr, *snd_ptr_orig, *rcv_ptr_orig;
  Integer me = pnga_nodeid();
  Integer p_handle, wproc;

  /* This routine makes use of the shift algorithm to update data in the
   * ghost cells bounding the local block of visible data. The shift
   * algorithm starts by updating the blocks of data along the first
   * dimension by grabbing a block of data that is width[0] deep but
   * otherwise matches the  dimensions of the data residing on the
   * calling processor. The update of the second dimension, however,
   * grabs a block that is width[1] deep in the second dimension but is
   * ldim0 + 2*width[0] in the first dimensions where ldim0 is the
   * size of the visible data along the first dimension. The remaining
   * dimensions are left the same. For the next update, the width of the
   * second dimension is also increased by 2*width[1] and so on. This
   * algorith makes use of the fact that data for the dimensions that
   * have already been updated is available on each processor and can be
   * used in the updates of subsequent dimensions. The total number of
   * separate updates is 2*ndim, an update in the negative and positive
   * directions for each dimension.
   *
   * This implementation make use of a combination of explicit message
   * passing between processors on different nodes and shared memory
   * copies with an additional flag between processors on the same node
   * to perform the update. Separate message types for the messages and
   * the use of the additional flag are for the updates in each
   * coordinate direction are used to maintain synchronization locally
   * and to guarantee that the data is present before the updates in a
   * new coordinate direction take place.
   *
   * To perform the update, this routine makes use of several copies of
   * indices marking the upper and lower limits of data. Indices
   * beginning with the character "p" are relative indices marking the
   * location of the data set relative to the origin the local patch of
   * the global array, all other indices are in absolute coordinates and
   * mark locations in the total global array. The indices used by this
   * routine are described below.
   *
   *       lo_loc[], hi_loc[]: The lower and upper indices of the visible
   *       block of data held by the calling processor.
   *
   *       lo_rcv[], hi_rcv[]: The lower and upper indices of the blocks
   *       of data that will be either sent to or received from a remote
   *       processor. These indices are NOT corrected for wrap-around
   *       (periodic) boundary conditions so they can be negative or greater
   *       than the array dimension values held in dims[].
   *
   *       slo_rcv[], shi_rcv[]: Similar to lo_rcv[] and hi_rcv[], except
   *       that these indices have been corrected for wrap-around
   *       boundary conditions.
   *
   *       plo_rcv[], phi_rcv[]: The local indices of the local data patch
   *       that receive that message from the remote processor.
   *
   *       plo_snd[], phi_snd[]: The local indices of the data patch
   *       that will be sent to the remote processor. Note that the
   *       dimensions of the patches represented by plo_rec[], plo_rec[] and
   *       plo_snd[], phi_snd[] must be the same.
   *
   *       tlo_rem[], thi_rem[]: The indices of the locally held visible
   *       portion of the global array on the remote processor that will be
   *       receiving the data using a shared memory copy.
   *
   *       plo_rem[], phi_rem[]: The local indices of the coordinate patch
   *       that will be put on the remote processor using a shared memory
   *       copy.
   */

  /* if global array has no ghost cells, just return */
  if (!pnga_has_ghosts(g_a)) return TRUE;

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  p_handle = GA[handle].p_handle;

  /* initialize range increments and get array dimensions */
  for (idx=0; idx < ndim; idx++) {
    increment[idx] = 0;
    width[idx] = (Integer)GA[handle].width[idx];
    dims[idx] = (Integer)GA[handle].dims[idx];
  }

  /* Check to make sure that global array is well-behaved (all processors
     have data and the width of the data in each dimension is greater
     than the corresponding value in width[]. */
  if (!gai_check_ghost_distr(g_a)) return FALSE;

  GA_PUSH_NAME("ga_update6_ghosts");
  msgcnt = 0;

  /* Get pointer to local memory */
  /*ptr_loc = GA[handle].ptr[me];*/
  /* obtain range of data that is held by local processor */
  pnga_distribution(g_a,me,lo_loc,hi_loc);
  /* Get indices of processor in virtual grid */
  pnga_proc_topology(g_a, me, index);

  /* Try to find maximum size of message that will be sent during
   * update operations and use this to allocate memory for message
   * passing buffers. */
  buftot = 1;
  for (i=0; i<ndim; i++) {
    buftot *= (hi_loc[i]-lo_loc[i] + 1 + 2*width[i]);
  }
  buflen = 1;
  for (i = 0; i < ndim; i++) {
    idir =  hi_loc[i] - lo_loc[i] + 1;
    if (buflen < (buftot/(idir + 2*width[i]))*width[i]) {
      buflen = (buftot/(idir + 2*width[i]))*width[i];
    }
  }
  bufsize = size*buflen;
  strcpy(send_name,"send_buffer");
  strcpy(rcv_name,"receive_buffer");
  snd_ptr_orig = snd_ptr = ga_malloc(buflen, GA[handle].type, send_name);
  rcv_ptr_orig = rcv_ptr = ga_malloc(buflen, GA[handle].type, rcv_name);

  /* loop over dimensions for sequential update using shift algorithm */
  msgcnt = 0;
  (*GA_Update_Signal) = 1;
  for (idx=0; idx < ndim; idx++) {

    /* Do not bother with update if nwidth is zero */
    if (width[idx] != 0) {

      /* Perform update in negative direction. */
      get_remote_block_neg(idx, ndim, lo_loc, hi_loc, slo_rcv, shi_rcv,
                           dims, width);
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rcv, shi_rcv, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rcv, shi_rcv, g_a);
      /* find out if this processor is on the same node */
      wproc = GA_proclist[0];
      if (p_handle >= 0) {
        wproc = PGRP_LIST[p_handle].inv_map_proc_list[wproc];
      }
      rprocflag = ARMCI_Same_node(wproc);
      proc_rem_snd = GA_proclist[0];

      /* Find processor from which data will be received */
      for (i = 0; i < ndim; i++) {
        if (i == idx) {
          lo_rcv[i] = hi_loc[i] + 1;
          hi_rcv[i] = hi_loc[i] + width[i];
        } else {
          lo_rcv[i] = lo_loc[i];
          hi_rcv[i] = hi_loc[i];
        }
      }

      /* Account for boundaries, if necessary. */
      for (i=0; i<ndim; i++) {
        if (i == idx) {
          if (hi_rcv[i] > dims[i]) {
            slo_rcv[i] = 1;
            shi_rcv[i] = width[i];
          } else {
            slo_rcv[i] = lo_rcv[i];
            shi_rcv[i] = hi_rcv[i];
          }
        } else {
          slo_rcv[i] = lo_rcv[i];
          shi_rcv[i] = hi_rcv[i];
        }
      }
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rcv, shi_rcv, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rcv, shi_rcv, g_a);
      wproc = GA_proclist[0];
      if (p_handle >= 0) {
        wproc = PGRP_LIST[p_handle].inv_map_proc_list[wproc];
      }
      sprocflag = ARMCI_Same_node(wproc);
      proc_rem_rcv = GA_proclist[0];
      pnga_distribution(g_a, proc_rem_rcv, tlo_rem, thi_rem);

      /* Get actual coordinates of chunk of data that will be sent to
       * remote processor as well as coordinates of the array space that
       * will receive data from remote processor. */
      for (i = 0; i < ndim; i++) {
        if (increment[i] == 0) {
          if (i == idx) {
            plo_snd[i] = width[i];
            phi_snd[i] = 2*width[i] - 1;
            plo_rcv[i] = hi_loc[i] - lo_loc[i] + width[i] + 1;
            phi_rcv[i] = hi_loc[i] - lo_loc[i] + 2*width[i];
            plo_rem[i] = thi_rem[i] - tlo_rem[i] + width[i] + 1;
            /*phi_rem[i] = thi_rem[i] - tlo_rem[i] + 2*width[i];*/
          } else {
            plo_snd[i] = width[i];
            phi_snd[i] = hi_loc[i] - lo_loc[i] + width[i];
            plo_rcv[i] = width[i];
            phi_rcv[i] = hi_loc[i] - lo_loc[i] + width[i];
            plo_rem[i] = width[i];
            /*phi_rem[i] = thi_rem[i] - tlo_rem[i] + width[i];*/
          }
        } else {
          plo_snd[i] = 0;
          phi_snd[i] = hi_loc[i] - lo_loc[i] + increment[i];
          plo_rcv[i] = 0;
          phi_rcv[i] = hi_loc[i] - lo_loc[i] + increment[i];
          plo_rem[i] = 0;
          /*phi_rem[i] = thi_rem[i] - tlo_rem[i] + increment[i];*/
        }
      }

      /* Get pointer to local data buffer and remote data
         buffer as well as lists of leading dimenstions */
      gam_LocationWithGhosts(me, handle, plo_snd, &ptr_snd, ld_loc);
      gam_LocationWithGhosts(me, handle, plo_rcv, &ptr_rcv, ld_loc);
      gam_LocationWithGhosts(proc_rem_snd, handle, plo_rem, &ptr_rem, ld_rem);

      /* Evaluate strides for send and receive */
      gam_setstride(ndim, size, ld_loc, ld_loc, stride_rcv,
          stride_snd);
      gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
          stride_snd);

      /* Compute the number of elements in each dimension and store
         result in count. Scale the first element in count by the
         element size. */
      gam_ComputeCount(ndim, plo_rcv, phi_rcv, count);
      gam_CountElems(ndim, plo_snd, phi_snd, &length);
      length *= size;
      count[0] *= size;

      /* If we are sending data to another node, then use message passing */
      if (!rprocflag) {
        /* Fill send buffer with data. */
        armci_write_strided(ptr_snd, (int)ndim-1, stride_snd, count, snd_ptr);
      }

      /* Send Messages. If processor has odd index in direction idx, it
       * sends message first, if processor has even index it receives
       * message first. Then process is reversed. Also need to account
       * for whether or not there are an odd number of processors along
       * update direction. */

      if (p_handle >= 0) {
        proc_rem_snd = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem_snd];
        proc_rem_rcv = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem_rcv];
      }
      if (GA[handle].nblock[idx]%2 == 0) {
        if (index[idx]%2 != 0 && !rprocflag) {
          armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
        } else if (index[idx]%2 == 0 && !sprocflag) {
          armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
        } 
        if (rprocflag) {
#if 0
          ARMCI_PutS(ptr_snd, stride_snd, ptr_rem, stride_rem, count, ndim- 1,
                     proc_rem_snd);
          /* Send signal to remote processor that data transfer has been completed. */
          bytes = sizeof(int);
          ARMCI_Put(GA_Update_Signal, GA_Update_Flags[proc_rem_snd]+msgcnt, bytes,
                    proc_rem_snd);
#else
          ARMCI_PutS_flag(ptr_snd, stride_snd, ptr_rem, stride_rem, count,
                          (int)(ndim-1), GA_Update_Flags[proc_rem_snd]+msgcnt,
                          *GA_Update_Signal, (int)proc_rem_snd);
#endif
        }
        if (index[idx]%2 != 0 && !sprocflag) {
          armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
        } else if (index[idx]%2 == 0 && !rprocflag) {
          armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
        }
      } else {
        /* account for wrap-around boundary condition, if necessary */
        pmax = GA[handle].nblock[idx] - 1;
        if (index[idx]%2 != 0 && !rprocflag) {
          armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
        } else if (index[idx]%2 == 0 && index[idx] != pmax && !sprocflag) {
          armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
        }
        if (rprocflag) {
#if 0
          ARMCI_PutS(ptr_snd, stride_snd, ptr_rem, stride_rem, count, ndim- 1,
                     proc_rem_snd);
          /* Send signal to remote processor that data transfer has been completed. */
          bytes = sizeof(int);
          ARMCI_Put(GA_Update_Signal, GA_Update_Flags[proc_rem_snd]+msgcnt, bytes,
                    proc_rem_snd);
#else
          ARMCI_PutS_flag(ptr_snd, stride_snd, ptr_rem, stride_rem, count,
                          (int)(ndim-1), GA_Update_Flags[proc_rem_snd]+msgcnt,
                          *GA_Update_Signal, (int)proc_rem_snd);
#endif
        }
        if (index[idx]%2 != 0 && !sprocflag) {
          armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
        } else if (index[idx]%2 == 0 && index[idx] != 0 && !rprocflag) {
          armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
        }
        /* make up for odd processor at end of string */
        if (index[idx] == 0 && !rprocflag) {
          armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
        }
        if (index[idx] == pmax && !sprocflag) {
          armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
        }
      }
      if (sprocflag) {
        flag1 = 1;
      } else {
        flag1 = 0;
      }
      msgcnt++;
      /* copy data back into global array */
      if (!sprocflag) {
        armci_read_strided(ptr_rcv, (int)ndim-1, stride_rcv, count, rcv_ptr);
      }

      /* Find parameters for message in positive direction. */
      get_remote_block_pos(idx, ndim, lo_loc, hi_loc, slo_rcv, shi_rcv,
                           dims, width);
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rcv, shi_rcv, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rcv, shi_rcv, g_a);
      wproc = GA_proclist[0];
      if (p_handle >= 0) {
        wproc = PGRP_LIST[p_handle].inv_map_proc_list[wproc];
      }
      rprocflag = ARMCI_Same_node(wproc);
      proc_rem_snd = GA_proclist[0];

      /* Find processor from which data will be recieved */
      for (i = 0; i < ndim; i++) {
        if (i == idx) {
          lo_rcv[i] = lo_loc[i] - width[i];
          hi_rcv[i] = lo_loc[i] - 1;
        } else {
          lo_rcv[i] = lo_loc[i];
          hi_rcv[i] = hi_loc[i];
        }
      }

      /* Account for boundaries, if necessary. */
      for (i=0; i<ndim; i++) {
        if (i == idx) {
          if (hi_rcv[i] < 1) {
            slo_rcv[i] = dims[i] - width[i] + 1;
            shi_rcv[i] = dims[i];
          } else {
            slo_rcv[i] = lo_rcv[i];
            shi_rcv[i] = hi_rcv[i];
          }
        } else {
          slo_rcv[i] = lo_rcv[i];
          shi_rcv[i] = hi_rcv[i];
        }
      }
      /* locate processor with this data */
      if (!pnga_locate_region(g_a, slo_rcv, shi_rcv, _ga_map,
          GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
          slo_rcv, shi_rcv, g_a);
      wproc = GA_proclist[0];
      if (p_handle >= 0) {
        wproc = PGRP_LIST[p_handle].inv_map_proc_list[wproc];
      }
      sprocflag = ARMCI_Same_node(wproc);
      proc_rem_rcv = GA_proclist[0];
      pnga_distribution(g_a, proc_rem_rcv, tlo_rem, thi_rem);

      /* Get actual coordinates of chunk of data that will be sent to
       * remote processor as well as coordinates of the array space that
       * will receive data from remote processor. */
      for (i = 0; i < ndim; i++) {
        if (increment[i] == 0) {
          if (i == idx) {
            plo_snd[i] = hi_loc[i] - lo_loc[i] + 1;
            phi_snd[i] = hi_loc[i] - lo_loc[i] + width[i];
            plo_rcv[i] = 0;
            phi_rcv[i] = width[i] - 1;
            plo_rem[i] = 0;
            /*phi_rem[i] = width[i] - 1;*/
          } else {
            plo_snd[i] = width[i];
            phi_snd[i] = hi_loc[i] - lo_loc[i] + width[i];
            plo_rcv[i] = width[i];
            phi_rcv[i] = hi_loc[i] - lo_loc[i] + width[i];
            plo_rem[i] = width[i];
            /*phi_rem[i] = thi_rem[i] - tlo_rem[i] + width[i];*/
          }
        } else {
          plo_snd[i] = 0;
          phi_snd[i] = hi_loc[i] - lo_loc[i] + increment[i];
          plo_rcv[i] = 0;
          phi_rcv[i] = hi_loc[i] - lo_loc[i] + increment[i];
          plo_rem[i] = 0;
          /*phi_rem[i] = thi_rem[i] - tlo_rem[i] + increment[i];*/
        }
      }

      /* Get pointer to local data buffer and remote data
         buffer as well as lists of leading dimenstions */
      gam_LocationWithGhosts(me, handle, plo_snd, &ptr_snd, ld_loc);
      gam_LocationWithGhosts(me, handle, plo_rcv, &ptr_rcv, ld_loc);
      gam_LocationWithGhosts(proc_rem_snd, handle, plo_rem, &ptr_rem, ld_rem);

      /* Evaluate strides for send and recieve */
      gam_setstride(ndim, size, ld_loc, ld_loc, stride_rcv,
          stride_snd);
      gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
          stride_snd);

      /* Compute the number of elements in each dimension and store
         result in count. Scale the first element in count by the
         element size. */
      gam_ComputeCount(ndim, plo_rcv, phi_rcv, count);
      gam_CountElems(ndim, plo_snd, phi_snd, &length);
      length *= size;
      count[0] *= size;

      /* if we are sending data to another node, use message passing */
      if (!rprocflag) {
        /* Fill send buffer with data. */
        armci_write_strided(ptr_snd, (int)ndim-1, stride_snd, count, snd_ptr);
      }

      /* Send Messages. If processor has odd index in direction idx, it
       * sends message first, if processor has even index it receives
       * message first. Then process is reversed. Also need to account
       * for whether or not there are an odd number of processors along
       * update direction. */

      if (p_handle >= 0) {
        proc_rem_snd = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem_snd];
        proc_rem_rcv = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem_rcv];
      }
      if (GA[handle].nblock[idx]%2 == 0) {
        if (index[idx]%2 != 0 && !rprocflag) {
          armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
        } else if (index[idx]%2 == 0 && !sprocflag) {
          armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
        } 
        if (rprocflag) {
#if 0
          ARMCI_PutS(ptr_snd, stride_snd, ptr_rem, stride_rem, count, ndim- 1,
                     proc_rem_snd);
          /* Send signal to remote processor that data transfer has been completed. */
          bytes = sizeof(int);
          ARMCI_Put(GA_Update_Signal, GA_Update_Flags[proc_rem_snd]+msgcnt, bytes,
                    proc_rem_snd);
#else
          ARMCI_PutS_flag(ptr_snd, stride_snd, ptr_rem, stride_rem, count,
                          (int)(ndim-1), GA_Update_Flags[proc_rem_snd]+msgcnt,
                          *GA_Update_Signal, (int)proc_rem_snd);
#endif
        }
        if (index[idx]%2 != 0 && !sprocflag) {
          armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
        } else if (index[idx]%2 == 0 && !rprocflag) {
          armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
        }
      } else {
        /* account for wrap-around boundary condition, if necessary */
        pmax = GA[handle].nblock[idx] - 1;
        if (index[idx]%2 != 0 && !rprocflag) {
          armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
        } else if (index[idx]%2 == 0 && index[idx] != 0 && !sprocflag) {
          armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
        }
        if (rprocflag) {
#if 0
          ARMCI_PutS(ptr_snd, stride_snd, ptr_rem, stride_rem, count, ndim- 1,
                     proc_rem_snd);
          /* Send signal to remote processor that data transfer has been completed. */
          bytes = sizeof(int);
          ARMCI_Put(GA_Update_Signal, GA_Update_Flags[proc_rem_snd]+msgcnt, bytes,
                    proc_rem_snd);
#else
          ARMCI_PutS_flag(ptr_snd, stride_snd, ptr_rem, stride_rem, count,
                          (int)(ndim-1), GA_Update_Flags[proc_rem_snd]+msgcnt,
                          *GA_Update_Signal, (int)proc_rem_snd);
#endif
        }
        if (index[idx]%2 != 0 && !sprocflag) {
          armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
        } else if (index[idx]%2 == 0 && index[idx] != pmax && !rprocflag) {
          armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
        }
        /* make up for odd processor at end of string */
        if (index[idx] == pmax && !rprocflag) {
          armci_msg_snd(msgcnt, snd_ptr, length, proc_rem_snd);
        }
        if (index[idx] == 0 && !sprocflag) {
          armci_msg_rcv(msgcnt, rcv_ptr, bufsize, &msglen, proc_rem_rcv);
        }
      }
      /* copy data back into global array */
      if (!sprocflag) {
        armci_read_strided(ptr_rcv, (int)ndim-1, stride_rcv, count, rcv_ptr);
      }
      if (sprocflag) {
        flag2 = 1;
      } else {
        flag2 = 0;
      }
      msgcnt++;
    }
    /* check to make sure any outstanding puts have showed up */
    waitformixedflags(flag1, flag2, GA_Update_Flags[GAme]+msgcnt-2,
                      GA_Update_Flags[GAme]+msgcnt-1);
    /* update increment array */
    increment[idx] = 2*width[idx];
  }

  ga_free(rcv_ptr_orig);
  ga_free(snd_ptr_orig);
  /* set update flags to zero for next operation */
  for (idx=0; idx < 2*ndim; idx++) {
    GA_Update_Flags[GAme][idx] = 0;
  }

  GA_POP_NAME;
  return TRUE;
}

/*\ UPDATE GHOST CELLS OF GLOBAL ARRAY USING GET CALLS
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_update7_ghosts = pnga_update7_ghosts
#endif
logical pnga_update7_ghosts(Integer g_a)
{
  Integer idx, ipx, np, handle=GA_OFFSET + g_a, proc_rem;
  Integer ntot, mask[MAXDIM];
  Integer size, ndim, i, itmp;
  Integer width[MAXDIM], dims[MAXDIM];
  Integer lo_loc[MAXDIM], hi_loc[MAXDIM];
  Integer plo_loc[MAXDIM], phi_loc[MAXDIM];
  Integer tlo_rem[MAXDIM], thi_rem[MAXDIM];
  Integer plo_rem[MAXDIM];
  Integer ld_loc[MAXDIM], ld_rem[MAXDIM];
  logical mask0;
  int stride_loc[MAXDIM], stride_rem[MAXDIM],count[MAXDIM];
  char *ptr_loc, *ptr_rem;
  Integer me = pnga_nodeid();
  Integer p_handle;

  /* if global array has no ghost cells, just return */
  if (!pnga_has_ghosts(g_a)) {
    return TRUE;
  }

  size = GA[handle].elemsize;
  ndim = GA[handle].ndim;
  p_handle = GA[handle].p_handle;
  /* initialize ghost cell widths and get array dimensions */
  for (idx=0; idx < ndim; idx++) {
    width[idx] = (Integer)GA[handle].width[idx];
    dims[idx] = (Integer)GA[handle].dims[idx];
  }

  /* Check to make sure that global array is well-behaved (all processors
     have data and the width of the data in each dimension is greater than
     the corresponding value in width[]). */
  if (!gai_check_ghost_distr(g_a)) return FALSE;

  GA_PUSH_NAME("ga_update7_ghosts");
  /* Get pointer to local memory */
  ptr_loc = GA[handle].ptr[me];
  /* obtain range of data that is held by local processor */
  pnga_distribution(g_a,me,lo_loc,hi_loc);

  /* evaluate total number of GET operations that will be required */
  ntot = 1;
  for (idx=0; idx < ndim; idx++) ntot *= 3;

  /* Loop over all GET operations. The operation corresponding to the
     mask of all zeros is left out. */
  for (ipx=0; ipx < ntot; ipx++) {
    /* Convert ipx to corresponding mask values */
    itmp = ipx;
    mask0 = TRUE;
    for (idx = 0; idx < ndim; idx++) {
      i = itmp%3;
      mask[idx] = i-1;
      if (mask[idx] != 0) mask0 = FALSE;
      itmp = (itmp-i)/3;
    }
    if (mask0) continue;

    /* check to see if ghost cell block has zero elements*/
    mask0 = FALSE;
    itmp = 0;
    for (idx = 0; idx < ndim; idx++) {
      if (mask[idx] != 0 && width[idx] == 0) mask0 = TRUE;
      if (mask[idx] != 0) itmp++;
    }
    /*if (itmp>1) mask0 = TRUE; */
    if (mask0) continue;
    /* Now that mask has been determined, find data that is to be moved
     * and identify processor from which it is coming. Wrap boundaries
     * around, if necessary */
    for (idx = 0; idx < ndim; idx++) {
      if (mask[idx] == 0) {
        tlo_rem[idx] = lo_loc[idx];
        thi_rem[idx] = hi_loc[idx];
      } else if (mask[idx] == -1) {
        if (lo_loc[idx] > 1) {
          tlo_rem[idx] = lo_loc[idx]-width[idx];
          thi_rem[idx] = lo_loc[idx]-1;
        } else {
          tlo_rem[idx] = dims[idx]-width[idx]+1;
          thi_rem[idx] = dims[idx];
        }
      } else if (mask[idx] == 1) {
        if (hi_loc[idx] < dims[idx]) {
          tlo_rem[idx] = hi_loc[idx] + 1;
          thi_rem[idx] = hi_loc[idx] + width[idx];
        } else {
          tlo_rem[idx] = 1;
          thi_rem[idx] = width[idx];
        }
      } else {
        fprintf(stderr,"Illegal mask value found\n");
      }
    }
    /* Locate remote processor to which data must be sent */
    if (!pnga_locate_region(g_a, tlo_rem, thi_rem, _ga_map,
       GA_proclist, &np)) ga_RegionError(pnga_ndim(g_a),
       tlo_rem, thi_rem, g_a);
    if (np > 1) {
      fprintf(stderr,"More than one remote processor found\n");
    }
    /* Remote processor has been identified, now get ready to send
       data to it. Start by getting distribution on remote
       processor.*/
    proc_rem = GA_proclist[0];
    pnga_distribution(g_a, proc_rem, tlo_rem, thi_rem);
    for (idx = 0; idx < ndim; idx++) {
      if (mask[idx] == 0) {
        plo_loc[idx] = width[idx];
        phi_loc[idx] = hi_loc[idx]-lo_loc[idx]+width[idx];
        plo_rem[idx] = plo_loc[idx];
      } else if (mask[idx] == -1) {
        plo_loc[idx] = 0;
        phi_loc[idx] = width[idx]-1;
        plo_rem[idx] = thi_rem[idx]-tlo_rem[idx]+1;
      } else if (mask[idx] == 1) {
        plo_loc[idx] = hi_loc[idx]-lo_loc[idx]+width[idx]+1;
        phi_loc[idx] = hi_loc[idx]-lo_loc[idx]+2*width[idx];
        plo_rem[idx] = width[idx];
      }
    }
    /* Get pointer to local data buffer and remote data
       buffer as well as lists of leading dimenstions */
    gam_LocationWithGhosts(me, handle, plo_loc, &ptr_loc, ld_loc);
    gam_LocationWithGhosts(proc_rem, handle, plo_rem, &ptr_rem, ld_rem);

    /* Evaluate strides on local and remote processors */
    gam_setstride(ndim, size, ld_loc, ld_rem, stride_rem,
                  stride_loc);

    /* Compute the number of elements in each dimension and store
       result in count. Scale the first element in count by the
       element size. */
    gam_ComputeCount(ndim, plo_loc, phi_loc, count);
    count[0] *= size;

    if (p_handle >= 0) {
      proc_rem = PGRP_LIST[p_handle].inv_map_proc_list[proc_rem];
    }
    /* put data on remote processor */
/*    ARMCI_GetS(ptr_rem, stride_rem, ptr_loc, stride_loc, count,
          (int)(ndim - 1), (int)proc_rem); */
    ARMCI_NbGetS(ptr_rem, stride_rem, ptr_loc, stride_loc, count,
          (int)(ndim - 1), (int)proc_rem, NULL);
  }

  ARMCI_WaitAll();
  GA_POP_NAME;
  return TRUE;
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_ghost_barrier = pnga_ghost_barrier
#endif
void pnga_ghost_barrier()
{
#ifdef LAPI
  int signal = 1, n = 1;
  int *ptr;
  ptr = &signal;
  armci_msg_igop(ptr,n,"+");
#else
  armci_msg_barrier();
#endif
}
/*\ UPDATE THE GHOST CELLS ON A PROCESSOR IN A SPECIFIC DIRECTION
 *  USING NON-BLOCKING GET CALLS
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_nbget_ghost_dir = pnga_nbget_ghost_dir
#endif
void pnga_nbget_ghost_dir(Integer g_a,
                               Integer *mask,
                               Integer *nbhandle)
{
  Integer handle = GA_OFFSET + g_a;
  Integer lo_loc[MAXDIM], hi_loc[MAXDIM], lo_rem[MAXDIM], hi_rem[MAXDIM];
  Integer subscript[MAXDIM], ld[MAXDIM];
  Integer i, ndim, dim, width;
  char *ptr_loc;
  Integer me = pnga_nodeid();
  /*Integer p_handle;*/
  GA_PUSH_NAME("nga_nbget_ghost_dir");
  ndim = GA[handle].ndim;
  /*p_handle = GA[handle].p_handle;*/
  /* check mask to see that it corresponds to a valid direction */
  for (i=0; i<ndim; i++) {
    if (abs(mask[i]) != 0 && abs(mask[i]) != 1)
      pnga_error("nga_nbget_ghost_dir: invalid mask entry", mask[i]);
  }

  /* get range of data on local processor */
  pnga_distribution(g_a,me,lo_loc,hi_loc);

  /* locate data on remote processor */
  for (i=0; i<ndim; i++) {
    dim = (Integer)GA[handle].dims[i];
    width = (Integer)GA[handle].width[i];
    if (mask[i] == 1) {
      if (hi_loc[i] == dim) {
        lo_rem[i] = 1;
        hi_rem[i] = width;
      } else {
        lo_rem[i] = hi_loc[i]+1;
        hi_rem[i] = hi_loc[i]+width;
      }
    } else if (mask[i] == -1) {
      if (lo_loc[i] == 1) {
        lo_rem[i] = dim - width + 1;
        hi_rem[i] = dim;
      } else {
        lo_rem[i] = lo_loc[i] - width;
        hi_rem[i] = lo_loc[i] - 1;
      }
    } else {
      lo_rem[i] = lo_loc[i];
      hi_rem[i] = hi_loc[i];
    }
  }
  
  /* Get pointer to data destination on local block. Start by
     by finding subscript to origin of destination block */
  for (i=0; i<ndim; i++) {
    if (mask[i] == 1) {
      subscript[i] = hi_loc[i]-lo_loc[i]+1+GA[handle].width[i];
    } else if (mask[i] == -1) {
      subscript[i] = 0;
    } else {
      subscript[i] = GA[handle].width[i];
    }
    ld[i] = hi_loc[i]-lo_loc[i]+1+2*GA[handle].width[i];
  }
  gam_LocationWithGhosts(me, handle, subscript, &ptr_loc, ld);
  /* get data */
  pnga_nbget(g_a,lo_rem,hi_rem,ptr_loc,ld,nbhandle);  
  GA_POP_NAME;
}

/*\ SET PRECOMPUTED INFO FOR UPDATING GHOST CELLS
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_set_ghost_info = pnga_set_ghost_info
#endif
logical pnga_set_ghost_info(Integer g_a)
{
  Integer handle = g_a + GA_OFFSET;
  if (GA[handle].cache != NULL)
    free(GA[handle].cache);
  GA[handle].cache = NULL;
  if (GA[handle].actv == 1) {
#ifdef CRAY_T3D
    return pnga_set_update5_info(g_a);
#else
    return pnga_set_update4_info(g_a);
#endif
  }
  return TRUE;
}

/*\ SET FLAG ON WHETHER OR NOT TO UPDATE GHOST CELL CORNER DATA
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_set_ghost_corner_flag = pnga_set_ghost_corner_flag
#endif
void pnga_set_ghost_corner_flag(Integer g_a, logical flag)
{
  Integer handle = g_a + GA_OFFSET;
  GA[handle].corner_flag = (int)flag;
  if (GA[handle].actv == 1) {
    (void)pnga_set_ghost_info(g_a);
  }
}

