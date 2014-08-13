#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: noncont.c,v 1.3.2.2 2007-05-04 16:43:35 d3p687 Exp $
 * noncont.c
 *
 * Developed by Andriy Kot <andriy.kot@pnl.gov>
 * Copyright (c) 2006 Pacific Northwest National Laboratory
 *
 * Alternative version of non-contiguous calls using non-blocking ones
 *
 * Changelog:
 * 2006-09-08 - created
 *
 */

#include "armcip.h"
#include "copy.h"
#include "acc.h"
#include "memlock.h"
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#ifdef PORTALS
#include "armci_portals.h"
#endif


#if 0
#   define PRN_DBG_MSG3(m,a1,a2,a3) \
           fprintf(stderr,"DBG %d: " m,armci_me,a1,a2,a3);fflush(stderr)
#   define PRN_DBG_MSG(m) PRN_DBG_MSG3(m,0,0,0)
#   define PRN_DBG_MSG1(m,a1) PRN_DBG_MSG3(m,a1,0,0)
#   define PRN_DBG_MSG2(m,a1,a2) PRN_DBG_MSG3(m,a1,a2,0)
#else
#   define PRN_DBG_MSG(m)
#   define PRN_DBG_MSG1(m,a1)
#   define PRN_DBG_MSG2(m,a1,a2)
#   define PRN_DBG_MSG3(m,a1,a2,a3)
#endif

#if 0
#   define CALL_IN(_func)  { if (armci_me == 0) printf("ENTERED %s\n", _func); fflush(stdout); }
#   define CALL_OUT(_func) { if (armci_me == 0) printf("EXITING %s\n", _func); fflush(stdout); }
#else
#   define CALL_IN(_func)
#   define CALL_OUT(_func)
#endif

#ifdef  NB_NONCONT

#if   defined(QUADRICS)
typedef ELAN_EVENT *HTYPE;
#define SHMEM_HANDLE_SUPPORTED
#elif defined(CRAY_SHMEM)
typedef void *HTYPE;
#else
typedef armci_ireq_t HTYPE;
#endif

#define MAX_SLOTS_LL    64
#define MIN_OUTSTANDING 6
static int max_pending = 16; /* throttle number of outstanding nb calls */

/* might have to use MAX_SLOTS_LL < MAX_PENDING due to throttling problem */
#define MAX_PENDING 6
#define ZR (HTYPE)0

static HTYPE put_dscr[MAX_SLOTS_LL];
static HTYPE get_dscr[MAX_SLOTS_LL];
/* static variables alreay initialize to 0 (?)
static HTYPE put_dscr[MAX_SLOTS_LL]= {
ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,
ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR};

static HTYPE get_dscr[MAX_SLOTS_LL] = {
ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,
ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR,ZR};
*/

#if defined(PORTALS)
extern ARMCI_MEMHDL_T *mhloc;
extern ARMCI_MEMHDL_T *mhrem;
#   define INI_HDL(_hdl, _op, _proc) {      \
            (_hdl).tag = GET_NEXT_NBTAG();  \
            (_hdl).op = _op;                \
            (_hdl).proc = _proc;            \
            (_hdl).bufid = NB_NONE;         \
           }
#   define CLR_HDL(_hdl) ((_hdl).tag = 0)
#   define CHK_HDL(_hdl) (_hdl.tag)
#else
#   define CLR_HDL(_hdl) ((_hdl) = ZR)
#   define CHK_HDL(_hdl) (_hdl)
#   define INI_HDL(_hdl, _op, _proc)
#endif

static int cur_get=0;
static int cur_put=0;
static int pending_get=0;
static int pending_put=0;

/* strided put, nonblocking */
void armcill_put2D(int proc, int bytes, int count, void* src_ptr,int src_stride,
                                                   void* dst_ptr,int dst_stride)
{
    CALL_IN("armcill_put2D");

    int _j, i, batch, issued=0;
    char *ps=src_ptr, *pd=dst_ptr;

    for (_j = 0;  _j < count;  ){
      /* how big a batch of requests can we issue */
      batch = (count - _j )<max_pending ? count - _j : max_pending;
      _j += batch;
#ifdef SHMEM_HANDLE_SUPPORTED
      for(i=0; i< batch; i++){
        if (CHK_HDL(put_dscr[cur_put])) armcill_nb_wait(put_dscr[cur_put]);
        else pending_put++;
        INI_HDL(put_dscr[cur_put], PUT, proc);
        armcill_nb_put(pd,ps,bytes,proc,put_dscr[cur_put]);
#else
      shmem_quiet();
      for(i=0; i< batch; i++){
        HTYPE dummy;
        armcill_nb_put(pd,ps,bytes,proc,dummy);
#endif
        issued++;
        ps += src_stride;
        pd += dst_stride;
        cur_put++;
        if(cur_put>=max_pending)cur_put=0;
      }
    }

    if(issued != count)
       armci_die2("armcill_put2D: mismatch %d %d \n", count,issued);

    CALL_OUT("armcill_put2D");
}


/* blocking vector put */
void armcill_putv(int proc, int bytes, int count, void* src[], void* dst[])
{
    int _j, i, batch, issued=0;
    void *ps, *pd;

    for (_j = 0;  _j < count;  ){
      /* how big a batch of requests can we issue */
      batch = (count - _j )<max_pending ? count - _j : max_pending;
      _j += batch;
#ifdef SHMEM_HANDLE_SUPPORTED
      for(i=0; i< batch; i++){
        if (CHK_HDL(put_dscr[cur_put])) armcill_nb_wait(put_dscr[cur_put]);
        else pending_put++;
        ps = src[issued];
        pd = dst[issued];
        INI_HDL(put_dscr[cur_put], PUT, proc);
        armcill_nb_put(pd,ps,bytes,proc,put_dscr[cur_put]);
#else
      shmem_quiet();
      for(i=0; i< batch; i++){
        HTYPE dummy;
        armcill_nb_put(pd,ps,bytes,proc,dummy);
#endif
        issued++;
        cur_put++;
        if(cur_put>=max_pending)cur_put=0;
      }
    }
    if(issued != count)
       armci_die2("armcill_putv: mismatch\n", count,issued);

#ifdef SHMEM_HANDLE_SUPPORTED
    for(i=0; i<max_pending; i++) if(CHK_HDL(put_dscr[i])){
        armcill_nb_wait(put_dscr[i]);
        CLR_HDL(put_dscr[i]);
    }
#else
    shmem_quiet();
#endif
}




/* strided get, nonblocking */
void armcill_get2D(int proc, int bytes, int count, void* src_ptr,int src_stride,
                                                   void* dst_ptr,int dst_stride)
{
    CALL_IN("armcill_get2D");
    PRN_DBG_MSG3("armcill_get2D: proc=%d, bytes=%d, count=%d\n", proc, bytes, count);

    int _j, i, batch, issued=0;
    char *ps=src_ptr, *pd=dst_ptr;

    for (_j = 0;  _j < count;  ){
      /* how big a batch of requests can we issue */
      batch = (count - _j )<max_pending ? count - _j : max_pending;
      _j += batch;
#ifdef SHMEM_HANDLE_SUPPORTED
      for(i=0; i< batch; i++){
        PRN_DBG_MSG2("inner loop: cur_ptr=%d, tag=%d\n", cur_get, get_dscr[cur_get].tag);
        if (CHK_HDL(get_dscr[cur_get])) armcill_nb_wait(get_dscr[cur_get]);
        else pending_get++;
        PRN_DBG_MSG1("inner loop: pending_get=%d\n", pending_get);
        INI_HDL(get_dscr[cur_get], GET, proc);
        armcill_nb_get(pd,ps,bytes,proc,get_dscr[cur_get]);
        PRN_DBG_MSG("inner loop: after get\n");
#else
      shmem_quiet();
      for(i=0; i< batch; i++){
        HTYPE dummy;
        armcill_nb_get(pd,ps,bytes,proc,dummy);
#endif
        issued++;
        ps += src_stride;
        pd += dst_stride;
        cur_get++;
        if(cur_get>=max_pending)cur_get=0;
      }
    }

    if(issued != count)
       armci_die2("armcill_get2D: mismatch %d %d \n", count,issued);

    CALL_OUT("armcill_get2D");
}


/* blocking vector get */
void armcill_getv(int proc, int bytes, int count, void* src[], void* dst[])
{
    int _j, i, batch, issued=0;
    void *ps, *pd;

    for (_j = 0;  _j < count;  ){
      /* how big a batch of requests can we issue */
      batch = (count - _j )<max_pending ? count - _j : max_pending;
      _j += batch;
#ifdef SHMEM_HANDLE_SUPPORTED
      for(i=0; i< batch; i++){
        if (CHK_HDL(get_dscr[cur_get])) armcill_nb_wait(get_dscr[cur_get]);
        else pending_get++;
        ps = src[issued];
        pd = dst[issued];
        INI_HDL(get_dscr[cur_get], GET, proc);
        armcill_nb_get(pd,ps,bytes,proc,get_dscr[cur_get]);
#else
      shmem_quiet();
      for(i=0; i< batch; i++){
        HTYPE dummy;
        armcill_nb_get(pd,ps,bytes,proc,dummy);
#endif
        issued++;
        cur_get++;
        if(cur_get>=max_pending)cur_get=0;
      }
    }
    if(issued != count)
       armci_die2("armcill_getv: mismatch %d %d \n", count,issued);

#ifdef SHMEM_HANDLE_SUPPORTED
    for(i=0; i<max_pending; i++) if(CHK_HDL(get_dscr[i])){
        armcill_nb_wait(get_dscr[i]);
        CLR_HDL(get_dscr[i]);
    }
#else
    shmem_quiet();
#endif
}


void armcill_wait_get()
{
    CALL_IN("armcill_wait_get");
#ifdef SHMEM_HANDLE_SUPPORTED
    int i;
    if(!pending_get)return;
    else pending_get=0;
    for(i=0; i<max_pending; i++) if(CHK_HDL(get_dscr[i])){
        armcill_nb_wait(get_dscr[i]);
        CLR_HDL(get_dscr[i]);
    }
#else
    shmem_quiet();
#endif
    CALL_OUT("armcill_wait_get");
}


void armcill_wait_put()
{
    CALL_IN("armcill_wait_put");
#ifdef SHMEM_HANDLE_SUPPORTED
    int i;
    if(!pending_put)return;
    else pending_put=0;
    for(i=0; i<max_pending; i++) if(CHK_HDL(put_dscr[i])){
        armcill_nb_wait(put_dscr[i]);
        CLR_HDL(put_dscr[i]);
    }
#else
    shmem_quiet();
#endif
    CALL_OUT("armcill_wait_put");
}


#endif/*NB_NONCONT*/
