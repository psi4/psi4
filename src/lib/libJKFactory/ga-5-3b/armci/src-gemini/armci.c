#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: armci.c,v 1.114.2.17 2007-08-30 22:58:18 manoj Exp $ */

/* DISCLAIMER
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

#define  EXTERN
/*#define  PRINT_BT*/
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#if defined(CRAY) && !defined(__crayx1)
#  include <sys/category.h>
#  include <sys/resource.h>
#  include <unistd.h>
#endif
#ifdef LAPI
#  include "lapidefs.h"
#endif
#include <errno.h>
#include "armcip.h"
#include "copy.h"
#include "memlock.h"
#include "shmem.h"
#include "signaltrap.h"

#ifdef ARMCIX
#include "x/armcix.h"
#endif
#ifdef BGML
#include "bgml.h"
#include <assert.h>
#include "bgmldefs.h"
extern void armci_msg_barrier(void);
#endif

#ifdef CRAY_SHMEM
#  ifdef CRAY_XT
#    include <mpp/shmem.h>
#  else
#    include <shmem.h>
#  endif
#endif

#include <malloc.h>

/* global variables */
int armci_me, armci_Sme, armci_nproc;
int armci_clus_me, armci_nclus, armci_master;
int armci_clus_first, armci_clus_last;
int *_armci_argc=NULL;
char ***_armci_argv=NULL;
int _armci_initialized_args=0;
int _armci_initialized=0;
int _armci_terminating =0;
thread_id_t armci_usr_tid;
armci_ireq_t armci_inb_handle[ARMCI_MAX_IMPLICIT];/*implicit non-blocking handle*/
#ifndef HITACHI
double armci_internal_buffer[BUFSIZE_DBL];
#endif
#if defined(SYSV) || defined(WIN32) || defined(MMAP) || defined(HITACHI) || defined(CATAMOUNT) || defined(BGML)
#   include "locks.h"
    lockset_t lockid;
#endif

int* armci_prot_switch_fence=NULL;
int armci_prot_switch_preproc = -1;
int armci_prot_switch_preop = -1;

#ifdef BGML
/*   void armci_allocate_locks(); */
   void armci_init_memlock();
#endif

#ifdef LIBELAN_ATOMICS
ELAN_ATOMIC *a;
#warning "Enabling new atomics"
#endif

typedef struct{
  int sent;
  int received;
  int waited;
}armci_notify_t;

armci_notify_t **_armci_notify_arr;

void ARMCI_Cleanup()
{
#if defined(DATA_SERVER)
#if defined(LIBONESIDED)
    dsTurnOff();
#else
    if(armci_nclus >1){
        armci_wait_for_server();
    }
#endif
#endif

#if (defined(SYSV) || defined(WIN32) || defined(MMAP))&& !defined(HITACHI) 
    Delete_All_Regions();
    if(armci_nproc>1)
#if !defined(LAPI) 
       DeleteLocks(lockid);
#endif

#ifndef WIN32
    ARMCI_RestoreSignals();
#endif

#endif
  armci_transport_cleanup();

}

int armci_getbufsize()
{
        return(BUFSIZE);
}

void armci_notify_init()
{
  int rc,bytes=sizeof(armci_notify_t)*armci_nproc;

#ifdef DOELAN4
  armci_elan_notify_init();
  return;
#endif

  _armci_notify_arr=
        (armci_notify_t**)malloc(armci_nproc*sizeof(armci_notify_t*));
  if(!_armci_notify_arr)armci_die("armci_notify_ini:malloc failed",armci_nproc);

  if((rc=PARMCI_Malloc((void **)_armci_notify_arr, bytes))) 
        armci_die(" armci_notify_init: armci_malloc failed",bytes); 
  bzero(_armci_notify_arr[armci_me], bytes);
}


static void armci_perror_msg()
{
    char perr_str[80];
    if(!errno)return;
    sprintf(perr_str,"Last System Error Message from Task %d:",armci_me);
    perror(perr_str);
}

static void armci_abort(int code)
{
    abort();
#if !defined(BGML)
    armci_perror_msg();
#endif
    ARMCI_Cleanup();
    /* data server process cannot use message-passing library to abort
     * it simply exits, parent will get SIGCHLD and abort the program
     */
#if defined(DATA_SERVER)
    if(armci_me<0)_exit(1);
    else
#endif
    armci_msg_abort(code);
}


void armci_die(char *msg, int code)
{
    void *bt[100];

    if(_armci_terminating)return;
    else _armci_terminating=1;

    if(SERVER_CONTEXT){
       fprintf(stdout,"%d(s):%s: %d\n",armci_me, msg, code); fflush(stdout);
    // fprintf(stderr,"%d(s):%s: %d\n",armci_me, msg, code);
    }else{
      fprintf(stdout,"%d:%s: %d\n",armci_me, msg, code); fflush(stdout);
    //fprintf(stderr,"%d:%s: %d\n",armci_me, msg, code);
    }

#ifdef PRINT_BT
    backtrace_symbols_fd(bt, backtrace(bt, 100), 2);
#endif

    armci_abort(code);
}


void armci_die2(char *msg, int code1, int code2)
{
    void *bt[100];

    if(_armci_terminating)return;
    else _armci_terminating=1;

    if(SERVER_CONTEXT){
      fprintf(stdout,"%d(s):%s: (%d,%d)\n",armci_me,msg,code1,code2);
      fflush(stdout);
      fprintf(stderr,"%d(s):%s: (%d,%d)\n",armci_me,msg,code1,code2);
    }else{
      fprintf(stdout,"%d:%s: (%d,%d)\n",armci_me,msg,code1,code2);
      fflush(stdout);
      fprintf(stderr,"%d:%s: (%d,%d)\n",armci_me,msg,code1,code2);
    }
#ifdef PRINT_BT
        backtrace_symbols_fd(bt, backtrace(bt, 100), 2);
#endif
    armci_abort(code1);
}


void ARMCI_Error(char *msg, int code)
{
    armci_die(msg,code);
}


void armci_allocate_locks()
{
    /* note that if ELAN_ACC is defined the scope of locks is limited to SMP */
#if !defined(CRAY_SHMEM) && (defined(HITACHI) || defined(CATAMOUNT) || \
    (defined(QUADRICS) && defined(_ELAN_LOCK_H) && !defined(ELAN_ACC)))
       armcill_allocate_locks(NUM_LOCKS);
#elif (defined(SYSV) || defined(WIN32) || defined(MMAP)) && !defined(HITACHI)
       if(armci_nproc == 1)return;
#  if defined(SPINLOCK) || defined(PMUTEX) || defined(PSPIN)
       CreateInitLocks(NUM_LOCKS, &lockid);
#  else
       if(armci_master==armci_me)CreateInitLocks(NUM_LOCKS, &lockid);
       armci_msg_clus_brdcst(&lockid, sizeof(lockid));
       if(armci_master != armci_me)InitLocks(NUM_LOCKS, lockid);
#  endif
#endif
}


void ARMCI_Set_shm_limit(unsigned long shmemlimit)
{
#if (defined(SYSV) || defined(WIN32)  || defined(MMAP)) && !defined(HITACHI)
#define EXTRASHM  1024   /* extra shmem used internally in ARMCI */
unsigned long limit;
    limit = shmemlimit + EXTRASHM;
    armci_set_shmem_limit_per_core(limit);
#endif
}


 
/*\ allocate and initialize memory locking data structure
\*/
void armci_init_memlock()
{
    int bytes = MAX_SLOTS*sizeof(memlock_t);
    int rc, msize_per_proc=bytes;
    
#ifdef MEMLOCK_SHMEM_FLAG    
    /* last proc on node allocates memlock flag in shmem */
    if(armci_clus_last == armci_me) bytes += sizeof(int);
#endif

    memlock_table_array = malloc(armci_nproc*sizeof(void*));
    if(!memlock_table_array) armci_die("malloc failed for ARMCI lock array",0);

    rc = PARMCI_Malloc(memlock_table_array, bytes);
    if(rc) armci_die("failed to allocate ARMCI memlock array",rc);

    armci_msg_barrier();

    bzero(memlock_table_array[armci_me],bytes);

#ifdef BGML
    bgml_init_locks ((void *) memlock_table_array[armci_me]);
#elif ARMCIX
    ARMCIX_init_memlock ((memlock_t *) memlock_table_array[armci_me]);
#endif


#ifdef MEMLOCK_SHMEM_FLAG    
    /* armci_use_memlock_table is a pointer to local memory variable=1
     * we overwrite the pointer with address of shared memory variable 
     * armci_use_memlock_table and initialize it >0
     */
    armci_use_memlock_table = (int*) (msize_per_proc + 
                              (char*) memlock_table_array[armci_clus_last]);  
                              
    /* printf("%d: last=%d bytes=%d ptr =(%d, %d)\n",
           armci_me,armci_clus_last,bytes,armci_use_memlock_table,  
           memlock_table_array[armci_clus_last]); fflush(stdout); */

    if(armci_clus_last == armci_me) *armci_use_memlock_table =1+armci_me;

#endif

    *armci_use_memlock_table = 0;
    armci_msg_barrier();
}


#if defined(SYSV) || defined(WIN32)
static void armci_check_shmmax()
{
  long mylimit, limit;
  mylimit = limit = (long) armci_max_region();
  armci_msg_bcast_scope(SCOPE_MASTERS, &limit, sizeof(long), 0);
  if(mylimit != limit){
     printf("%d:Shared mem limit in ARMCI is %ld bytes on node %s vs %ld on %s\n",
            armci_me,mylimit<<10,armci_clus_info[armci_clus_me].hostname,
            limit<<10, armci_clus_info[0].hostname);
     fflush(stdout); sleep(1);
     armci_die("All nodes must have the same SHMMAX limit if NO_SHM is not defined",0);
  }
}
#endif

extern void armci_region_shm_malloc(void *ptr_arr[], size_t bytes);


void ARMCI_NetInit()
{
  /*armci_portals_net_init();*/
}

int PARMCI_Init_args(int *argc, char ***argv)
{
    armci_msg_init(argc,argv);

    _armci_argc = argc;
    _armci_argv = argv;
    _armci_initialized_args=1;
    PARMCI_Init();
}


extern void *sbrk(intptr_t);
extern void code_summary();

int PARMCI_Init()
{
    caddr_t atbeginbrval = (caddr_t)sbrk(0);
    if(_armci_initialized>0) return 0;
#ifdef NEW_MALLOC
    mallopt(M_MMAP_MAX, 0);
    mallopt(M_TRIM_THRESHOLD, -1);
#endif

    armci_msg_init(NULL, NULL);

    armci_nproc = armci_msg_nproc();
    armci_me = armci_msg_me();
    armci_usr_tid = THREAD_ID_SELF(); /*remember the main user thread id */
    armci_init_clusinfo();
    armci_prot_switch_fence = malloc(sizeof(int*)*armci_nproc);
    assert(armci_prot_switch_fence !=NULL);
  # ifdef LIBONESIDED
    armci_onesided_init();
  # endif
#ifdef MPI
    armci_group_init();
#endif
#ifndef NEW_MALLOC
    armci_krmalloc_init_localmem();
#endif
#if defined(SYSV) || defined(WIN32) || defined(MMAP)
    if(ARMCI_Uses_shm() ) {
      armci_shmem_init();
    }
#endif
    armci_allocate_locks();
    armci_init_fence();
#if ARMCI_ENABLE_GPC_CALLS
    gpc_init_signals();
#endif
    armci_msg_barrier();
    armci_init_memlock(); /* allocate data struct for locking memory areas */
    armci_msg_barrier();
    //if(armci_me == 0) code_summary();
    armci_msg_barrier();
    armci_msg_gop_init();
    _armci_initialized++;
    return 0;
}


void PARMCI_Finalize()
{
    if(!_armci_initialized)return;
    _armci_initialized--;
    if(_armci_initialized)return;

    _armci_terminating =1;
    armci_msg_barrier();
    if(armci_me==armci_master) ARMCI_ParentRestoreSignals();

#ifdef PORTALS
    request_header_t msg;
    portals_ds_req_t req;
    ptl_process_id_t dsid = portals_id_map[armci_me];
    msg.operation = QUIT;

    if(armci_me == armci_master) {
       portalsBlockingRemoteOperationToNode(&msg,sizeof(request_header_t),armci_clus_me);
    }

    armci_msg_barrier();
    portals_cp_finalize();

#else

    ARMCI_Cleanup();
    armci_msg_barrier();
    armci_group_finalize();
    free(armci_prot_switch_fence);
#endif
#ifdef MPI
    MPI_Comm_free(&ARMCI_COMM_WORLD); /*JD: free at last*/
#endif
}


/* Indicates whether ARMCI_Init or ARMCI_Init_args has been called. */
int PARMCI_Initialized()
{
    return (_armci_initialized > 0) ? 1 : 0;
}


#if !(defined(SYSV) || defined(WIN32))
void ARMCI_Set_shmem_limit(unsigned long shmemlimit)
{
   /* not applicable here
    * aborting would  make user's life harder
    */
}
#endif



void ARMCI_Copy(void *src, void *dst, int n)
{
 armci_copy(src,dst,n);
}

extern void cpu_yield();
void armci_util_wait_int(volatile int *p, int val, int maxspin)
{
int count=0;
extern void cpu_yield();
       while(*p != val)
            if((++count)<maxspin);
            else{
               /*printf("\n%d:flag=%d val=%d",armci_me,*p,val);*/
               cpu_yield();
               count =0;
#if defined(MACX) && defined(__ppc__) && defined(__GNUC__)
               __asm__ __volatile__ ("sync" ::: "memory");
#endif
               __asm__ __volatile__ ("mfence" ::: "memory");
               __asm__ __volatile__ ("sfence" ::: "memory");
            }
}
void armci_util_wait_long(volatile long *p, long val, int maxspin)
 {
int count=0;
extern void cpu_yield();
       while(*p != val)
            if((++count)<maxspin);
            else{
               /*printf("\n%d:flag=%d val=%d",armci_me,*p,val);*/
               cpu_yield();
               count =0;
#if defined(MACX) && defined(__ppc__) && defined(__GNUC__)
               __asm__ __volatile__ ("sync" ::: "memory");
#endif
               __asm__ __volatile__ ("mfence" ::: "memory");
               __asm__ __volatile__ ("sfence" ::: "memory");
            }
} 

/*\ returns 1 if specified process resides on the same smp node as calling task 
\*/
int ARMCI_Same_node(int proc)
{
   int direct=SAMECLUSNODE(proc);
   return direct;
}

/*\ blocks the calling process until a nonblocking operation represented
 *  by the user handle completes
\*/
int PARMCI_Wait(armci_hdl_t* usr_hdl)
{
        armci_ihdl_t nb_handle = (armci_ihdl_t)usr_hdl;
        int i,success=0;
        int direct=SAMECLUSNODE(nb_handle->proc);

        if(direct) {
           return(success);
        }

        if(nb_handle) {

           if(nb_handle->onesided_direct) {
              for(i=0; i<MAX_OUTSTANDING_ONESIDED_GETS; i++) {
                 if(nb_handle->comm_desc[i].state) {
                    onesided_wait(&nb_handle->comm_desc[i]);
                    cpMemDeregister(&nb_handle->comm_desc[i].local_mdesc);
                 }
              }
              __asm__ __volatile__ ("mfence" ::: "memory");
              __asm__ __volatile__ ("sfence" ::: "memory");
              ARMCI_INIT_HANDLE(nb_handle);
              return(success);
           }

           if(nb_handle->agg_flag) {
              armci_agg_complete(nb_handle, UNSET);
              return (success);
           }
      
           if(nb_handle->tag!=0 && nb_handle->bufid==NB_NONE) {
              ARMCI_NB_WAIT(nb_handle->cmpl_info);
              __asm__ __volatile__ ("mfence" ::: "memory");
              __asm__ __volatile__ ("sfence" ::: "memory");
              return(success);
           }
     
         # ifdef COMPLETE_HANDLE
           COMPLETE_HANDLE(nb_handle->bufid,nb_handle->tag,(&success));
         # endif
        }

        __asm__ __volatile__ ("mfence" ::: "memory");
        __asm__ __volatile__ ("sfence" ::: "memory");
        return(success);
}

/** 
 * implicit handle 
 */
static char hdl_flag[ARMCI_MAX_IMPLICIT];
static int impcount=0;
armci_ihdl_t armci_set_implicit_handle (int op, int proc) {
 
  int i=impcount%ARMCI_MAX_IMPLICIT;
  if(hdl_flag[i]=='1')
    PARMCI_Wait((armci_hdl_t*)&armci_inb_handle[i]);

#ifdef BGML
   armci_inb_handle[i].count=0;
#endif
  armci_inb_handle[i].tag   = GET_NEXT_NBTAG();
  armci_inb_handle[i].op    = op;
  armci_inb_handle[i].proc  = proc;
  armci_inb_handle[i].bufid = NB_NONE;
  armci_inb_handle[i].agg_flag = 0;
  hdl_flag[i]='1';
  ++impcount;
  return &armci_inb_handle[i];
}
 
 
/* wait for all non-blocking operations to finish */
int PARMCI_WaitAll (void) {
#ifdef BGML
  BGML_WaitAll();
#elif ARMCIX
  ARMCIX_WaitAll ();
#else
  int i;
  if(impcount) {
    for(i=0; i<ARMCI_MAX_IMPLICIT; i++) {
      if(hdl_flag[i] == '1') {
        PARMCI_Wait((armci_hdl_t*)&armci_inb_handle[i]);
        hdl_flag[i]='0';
      }
    }
  }
  impcount=0;
#endif
  return 0;
}
 
/* wait for all non-blocking operations to a particular process to finish */
int PARMCI_WaitProc (int proc) {
#ifdef BGML
  BGML_WaitProc(proc);
#elif ARMCIX
  ARMCIX_WaitProc (proc);
#else
  int i;
  if(impcount) {
    for(i=0; i<ARMCI_MAX_IMPLICIT; i++) {
      if(hdl_flag[i]=='1' && armci_inb_handle[i].proc==proc) {
        PARMCI_Wait((armci_hdl_t*)&armci_inb_handle[i]);
        hdl_flag[i]='0';
      }
    }
  }
#endif
  return 0;
}

static unsigned int _armci_nb_tag=0;
unsigned int _armci_get_next_tag(){
    return((++_armci_nb_tag));
}

void ARMCI_SET_AGGREGATE_HANDLE(armci_hdl_t* nb_handle) { 
      ((armci_ihdl_t)(nb_handle))->agg_flag = 1;
      ((armci_ihdl_t)(nb_handle))->proc = -1;
}
 
void ARMCI_UNSET_AGGREGATE_HANDLE(armci_hdl_t* nb_handle) {
      ((armci_ihdl_t)(nb_handle))->agg_flag = 0;
      ((armci_ihdl_t)(nb_handle))->proc = -1;
}

int parmci_notify(int proc)
{
#ifdef DOELAN4
  if(proc==armci_me){
    return 0;
  }
#endif
#if defined(GM) || (defined(DOELAN4) && defined(ELAN_ACC))
  {
    extern int armci_inotify_proc(int);
    return(armci_inotify_proc(proc));
  }
#else
   armci_notify_t *pnotify = _armci_notify_arr[armci_me]+proc;
   pnotify->sent++;
# ifdef MEM_FENCE
   if(SAMECLUSNODE(proc)) MEM_FENCE;
# endif
   PARMCI_Put(&pnotify->sent,&(_armci_notify_arr[proc]+armci_me)->received, 
             sizeof(pnotify->sent),proc);
   return(pnotify->sent);
#endif
}


/*\ blocks until received count becomes >= waited count
 *  return received count and store waited count in *pval
\*/
int parmci_notify_wait(int proc,int *pval)
{
  int retval;
#ifdef DOELAN4
  if(proc==armci_me){
#ifdef MEM_FENCE
       MEM_FENCE;
#endif
    return 0;
  }
#endif

#if defined(GM) || (defined(DOELAN4) && defined(ELAN_ACC))
  {
     extern int armci_inotify_wait(int,int*);
     retval=armci_inotify_wait(proc,pval);
  }
#else
  {
     long loop=0;
     armci_notify_t *pnotify = _armci_notify_arr[armci_me]+proc;
     pnotify->waited++;
     while( pnotify->waited > pnotify->received) {
         if(++loop == 1000) { loop=0;cpu_yield(); }
         armci_util_spin(loop, pnotify);
     }
     *pval = pnotify->waited;
     retval=pnotify->received;
  }
#endif

  return retval;
}

long armci_util_long_getval(long* p)
{
   return *p;
}

int armci_util_int_getval(int* p)
{
   return *p;
}


int PARMCI_Test(armci_hdl_t *usr_hdl)
{
armci_ihdl_t nb_handle = (armci_ihdl_t)usr_hdl;
int success=0;
#ifdef BGML
   success=(int)nb_handle->count;
#else
int direct=SAMECLUSNODE(nb_handle->proc);
   if(direct)return(success);
    if(nb_handle) {
      if(nb_handle->agg_flag) {
         armci_die("test for aggregate handle not yet implemented\n",0);
      }
    }
    if(nb_handle){
#     ifdef ARMCI_NB_TEST
        if(nb_handle->tag==0){
              ARMCI_NB_TEST(nb_handle->cmpl_info,&success);
              return(success);
        }
#       ifdef LAPI
         if(nb_handle->tag!=0 && nb_handle->bufid==NB_NONE){
               ARMCI_NB_TEST(nb_handle->cmpl_info,&success);
               return(success);
         }
#       endif
#     endif
#     ifdef TEST_HANDLE
       TEST_HANDLE(nb_handle->bufid,nb_handle->tag,(&success));
#     endif
    }
#endif
    return(success);
}

#ifdef DO_CKPT
void ARMCI_Ckpt_create_ds(armci_ckpt_ds_t *ckptds, int count)
{
    armci_create_ckptds(ckptds,count);
}

int ARMCI_Ckpt_init(char *filename, ARMCI_Group *grp, int savestack, int saveheap, armci_ckpt_ds_t *ckptds)
{
int rid;
    rid = armci_icheckpoint_init(filename,grp,savestack,saveheap,ckptds);
    return(rid);
}

int ARMCI_Ckpt(int rid)
{
    return(armci_icheckpoint(rid));
}

void ARMCI_Ckpt_Recover(int rid, int iamreplacement)
{
    armci_irecover(rid, iamreplacement);
}
void ARMCI_Ckpt_finalize(int rid)
{
    armci_icheckpoint_finalize(rid);
}
#endif
#if ARMCI_ENABLE_GPC_CALLS
int armci_gpc(int hndl, int proc, void  *hdr, int hlen,  void *data,  int dlen,
              void *rhdr, int rhlen, void *rdata, int rdlen,
              armci_hdl_t* nbh) {
armci_ihdl_t nb_handle = (armci_ihdl_t)nbh;
armci_giov_t darr[2]; /* = {{&rhdr, &rhdr, 1, rhlen}, {&rdata, &rdata, 1, rdlen}};*/
gpc_send_t send;
char *ptr;

    /* initialize giov */
    darr[0].src_ptr_array = &rhdr;
    darr[0].dst_ptr_array = &rhdr;
    darr[0].ptr_array_len = 1;
    darr[0].bytes         = rhlen;

    darr[1].src_ptr_array = &rdata;
    darr[1].dst_ptr_array = &rdata;
    darr[1].ptr_array_len = 1;
    darr[1].bytes         = rdlen;

  
/*    if(hlen<0 || hlen>=ARMCI_Gpc_get_hlen()) */
/*      return FAIL2; */
/*    if(rhlen<0 || rhlen>=ARMCI_Gpc_get_hlen()) */
/*      return FAIL2; */
/*    if(dlen<0 || dlen>=ARMCI_Gpc_get_dlen())  */
/*      return FAIL2; */
/*    if(rdlen<0 || rdlen>=ARMCI_Gpc_get_dlen())  */
/*      return FAIL2; */

    if(hlen>0 && hdr==NULL) 
      return FAIL3;
    if(rhlen>0 && rhdr==NULL) 
      return FAIL3;
    if(dlen>0 && data==NULL) 
      return FAIL3;
    if(rdlen>0 && rdata==NULL) 
      return FAIL3;

    if(proc<0 || proc >= armci_nproc)
      return FAIL4;

    send.hndl = hndl;
    send.hlen = hlen;
    send.dlen = dlen;
    send.hdr = hdr;
    send.data = data;

    if(nb_handle){
      nb_handle->tag = GET_NEXT_NBTAG();
      nb_handle->op  = GET;
      nb_handle->proc= proc;
      nb_handle->bufid=NB_NONE;
    }
    else {
      ORDER(GET,proc); /*ensure ordering */      
      nb_handle = NULL;
    }  

#if defined(LAPI) || defined(GM) || defined(VAPI) || defined(QUADRICS)
    if(armci_rem_gpc(GET, darr, 2, &send, proc, 1, nb_handle))
#endif
      return FAIL2;
    return 0;
}

int armci_sameclusnode(int proc) {
  return SAMECLUSNODE(proc);
}
#endif

void _armci_init_handle(armci_hdl_t *hdl)
{
    ((double *)((hdl)->data))[0]=0;
    ((double *)((hdl)->data))[1]=0;
}

static inline int val_to_char(int v)
{
    if (v >= 0 && v < 10)
      return '0' + v;
    else if (v >= 10 && v < 16)
      return ('a' - 10) + v;
    else
      return -1;
}
static const char *nexttoken(const char *q, int sep)
{
    if (q)
      q = strchr(q, sep);
    if (q)
      q++;
    return q;
}

#ifdef PORTALS_UNRESOLVED
int cstr_to_cpuset(cpu_set_t * mask, const char *str)
{
const char     *p, *q;
q = str;
    CPU_ZERO(mask);

    while (p = q, q = nexttoken(q, ','), p) {
    unsigned int    a;	/* beginning of range */
    unsigned int    b;	/* end of range */
    unsigned int    s;	/* stride */
    const char     *c1, *c2;
      if (sscanf(p, "%u", &a) < 1)
        return 1;
      b = a;
      s = 1;
      c1 = nexttoken(p, '-');
      c2 = nexttoken(p, ',');
      if (c1 != NULL && (c2 == NULL || c1 < c2)) {
        if (sscanf(c1, "%u", &b) < 1)
          return 1;
        c1 = nexttoken(c1, ':');
        if (c1 != NULL && (c2 == NULL || c1 < c2))
          if (sscanf(c1, "%u", &s) < 1) {
            return 1;
          }
      }
      if (!(a <= b))
        return 1;
      while (a <= b) {
        CPU_SET(a, mask);
        a += s;
      }
    }
    return 0;
}

char *cpuset_to_cstr(cpu_set_t * mask, char *str)
{
int             i;
char           *ptr = str;
int             entry_made = 0;
    for (i = 0; i < CPU_SETSIZE; i++) {
      if (CPU_ISSET(i, mask)) {
      int             j;
      int             run = 0;
        entry_made = 1;
        for (j = i + 1; j < CPU_SETSIZE; j++) {
          if (CPU_ISSET(j, mask))
            run++;
          else
            break;
        }
        if (!run)
          sprintf(ptr, "%d,", i);
        else if (run == 1) {
          sprintf(ptr, "%d,%d,", i, i + 1);
          i++;
        } else {
          sprintf(ptr, "%d-%d,", i, i + run);
          i += run;
        }
        while (*ptr != 0)
          ptr++;
      }
    }
    ptr -= entry_made;
    *ptr = 0;
    return str;
}

char *cpuset_to_str(cpu_set_t * mask, char *str)
{
int             base;
char           *ptr = str;
char           *ret = 0;
    for (base = CPU_SETSIZE - 4; base >= 0; base -= 4) {
    char    val = 0;
      if (CPU_ISSET(base, mask))
        val |= 1;
      if (CPU_ISSET(base + 1, mask))
        val |= 2;
      if (CPU_ISSET(base + 2, mask))
        val |= 4;
      if (CPU_ISSET(base + 3, mask))
        val |= 8;
      if (!ret && val)
        ret = ptr;
      *ptr++ = val_to_char(val);
    }
    *ptr = 0;
    return ret ? ret : ptr - 1;
}
#endif


long armci_cksm_copy(char *src, char *dst, size_t bytes)
{
long sum = 0;
size_t count=bytes;
  while( count > 1 )  {
    sum += * (unsigned int *) src++;
    count -= 4;
  }

  if( count > 0 ){
          printf("\nblistering barnicles");
    sum += * (unsigned char *) src;
  }

  while (sum>>32)
    sum = (sum & 0xffffffff) + (sum >> 32);
  return(~sum);
}

void code_summary() {
  printf("\nActive #defines that could affect ARMCI");
  printf("\n----------------------------------------");
# ifdef ORNL_USE_DS_FOR_REMOTE_GETS
  printf("\n#define ORNL_USE_DS_FOR_REMOTE_GETS");
# endif

# ifdef PORTALS_USE_RENDEZ_VOUS
  printf("\n#define PORTALS_USE_RENDEZ_VOUS");
# endif  

# ifdef PORTALS_LIMIT_REMOTE_REQUESTS_BY_NODE
  printf("\n#define PORTALS_LIMIT_REMOTE_REQUESTS_BY_NODE");
# endif

# ifdef PORTALS_AFFINITY
  printf("\n#define PORTALS_AFFINITY");
# endif

/* 
# ifdef CRAY_USE_MDMD_COPY
  printf("\n#define CRAY_USE_MDMD_COPY");
# endif
*/
  printf("\n----------------------------------------");
  printf("\nInfo @ armci/src/code_options.h");
  printf("\n----------------------------------------\n");

# ifdef PORTALS
  portals_print_summary();
# endif
}
