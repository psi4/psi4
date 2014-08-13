#if HAVE_CONFIG_H
#   include "config.h"
#endif

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

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STDARG_H
#   include <stdarg.h>
#endif
#if defined(CRAY) && !defined(__crayx1)
#  include <sys/category.h>
#  include <sys/resource.h>
#  if HAVE_UNISTD_H
#   include <unistd.h>
#  endif
#endif
#ifdef LAPI
#  include "lapidefs.h"
#endif
#if HAVE_ERRNO_H
#   include <errno.h>
#endif
#include "armcip.h"
#include "copy.h"
#include "memlock.h"
#include "shmem.h"
#include "signaltrap.h"

#ifdef ARMCIX
#include "armcix.h"
#endif
#ifdef BGML
#include "bgml.h"
#if HAVE_ASSERT_H
#   include <assert.h>
#endif
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

/* global variables -- Initialized in PARMCI_Init() and never modified*/
int armci_me, armci_nproc;
int armci_clus_me, armci_nclus, armci_master;
int armci_clus_first, armci_clus_last;
int _armci_initialized=0;
int _armci_initialized_args=0;
int _armci_terminating =0;
int *_armci_argc=NULL;
char ***_armci_argv=NULL;
thread_id_t armci_usr_tid;

#if !defined(HITACHI) && !defined(THREAD_SAFE)
double armci_internal_buffer[BUFSIZE_DBL];
#endif
#if defined(SYSV) || defined(WIN32) || defined(MMAP) || defined(HITACHI) || defined(CATAMOUNT) || defined(BGML)
#   include "locks.h"
    lockset_t lockid;
#endif

#ifdef ALLOW_PIN
int* armci_prot_switch_fence=NULL;
int armci_prot_switch_preproc = -1;
int armci_prot_switch_preop = -1;
#endif

#ifdef BGML
/*   void armci_allocate_locks(); */
   void armci_init_memlock();
#endif


typedef struct{
  int sent;
  int received;
  int waited;
}armci_notify_t;

armci_notify_t **_armci_notify_arr;

#ifdef CRAY_XT
int _armci_malloc_local_region;
#endif

void ARMCI_Cleanup()
{
#if (defined(SYSV) || defined(WIN32) || defined(MMAP))&& !defined(HITACHI) 
    Delete_All_Regions();
    if(armci_nproc>1)
#if !defined(LAPI) 
       DeleteLocks(lockid);
#endif

    /* in case of an error notify server that it is time to quit */
#if defined(DATA_SERVER)
    if(armci_nclus >1){
        /* send quit request to server unless it is already dead */
        armci_wait_for_server();
        armci_transport_cleanup();
    }
#endif
    armci_finalize_fence();
#ifndef WIN32
    ARMCI_RestoreSignals();
#endif
#endif
}


void armci_notify_init()
{
  int rc, bytes= sizeof(armci_notify_t)*armci_nproc;

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
    if(!errno)
        return;
    sprintf(perr_str,"Last System Error Message from Task %d:",armci_me);
    perror(perr_str);
}


#if defined(IBM) || defined(IBM64)
int AR_caught_sigint;
int AR_caught_sigterm;
#else
extern int AR_caught_sigint;
extern int AR_caught_sigterm;
#endif

void armci_abort(int code)
{
#if !defined(BGML)
    armci_perror_msg();
#endif
    ARMCI_Cleanup();

    /* data server process cannot use message-passing library to abort
     * it simply exits, parent will get SIGCHLD and abort the program
     */
#if defined(IBM) || defined(IBM64)
     /* hack for a problem in POE signal handlers in non-LAPI MPI  */
     if(AR_caught_sigint || AR_caught_sigterm) 
         _exit(1);
#endif

#if defined(DATA_SERVER)
    if(armci_me<0)
        _exit(1);
    else
#endif
    armci_msg_abort(code);
}


/*For now, until no code requires a function pointer to ARMCI_Error
  (used by GA now).*/
void ARMCI_Error(char *msg, int code)
{
    armci_die(msg,code);
}



void armci_allocate_locks()
{
    /* note that if ELAN_ACC is defined the scope of locks is limited to SMP */
#if !defined(CRAY_SHMEM) && \
    ( defined(HITACHI) || defined(CATAMOUNT) || \
      (defined(QUADRICS) && defined(_ELAN_LOCK_H) && !defined(ELAN_ACC)) )
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
    limit = shmemlimit+EXTRASHM;
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

    armci_msg_barrier();
}


#if defined(SYSV) || defined(WIN32) || defined(MMAP)
#   if defined(QUADRICS) && !defined(NO_SHM)
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
#   endif
#endif

extern void armci_region_shm_malloc(void *ptr_arr[], size_t bytes);

#ifdef ENABLE_CHECKPOINT
int armci_ft_spare_procs;
void armci_set_spare_procs(int spare)
{
    armci_ft_spare_procs = spare;
}
ARMCI_Group armci_ft_group;
ARMCI_Group *ARMCI_Get_ft_group()
{
    return(&armci_ft_group);
}
void armci_create_ft_group()
{
    int i, list[MAX_PROC];
    for(i=0;i<armci_nproc-armci_ft_spare_procs;i++)list[i] = i;
    printf("\n%d:here ok\n",armci_me);fflush(stdout);
    ARMCI_Group_create((armci_nproc-armci_ft_spare_procs),list,&armci_ft_group);
    printf("\n%d:done with group create\n",armci_me);fflush(stdout);
}
#endif

int PARMCI_Init_args(int *argc, char ***argv) 
{
    armci_msg_init(argc,argv);

#ifdef MPI_SPAWN
    /* If this is data server process, then it should call
     * armci_mpi2_server_init() instead of PARMCI_Init(). PARMCI_Init() should
     * only be called by clients */
    {
       MPI_Comm parent_comm;
       MPI_Comm_get_parent(&parent_comm);
       if (parent_comm != MPI_COMM_NULL) 
           armci_mpi2_server();
    }
#endif

    _armci_argc = argc;
    _armci_argv = argv;
    _armci_initialized_args=1;
    return PARMCI_Init();
}

void _armci_test_connections() 
{
  int i;
  int nprocs = armci_msg_nproc();
  int **bufs = (int **)malloc(nprocs*sizeof(int*));
  int *val;
  dassert(1, bufs);

  PARMCI_Malloc((void **)bufs, sizeof(int));
  dassert(1, bufs[armci_me]);
  val = PARMCI_Malloc_local(sizeof(int));
  dassert(1, val);

  for(i=0; i<nprocs; i++) {
    dassert(1, bufs[i]);
    PARMCI_Put(val, bufs[i], sizeof(int), i);
  }
  PARMCI_AllFence();
  PARMCI_Barrier();
  PARMCI_Free(bufs[armci_me]);
  PARMCI_Free_local(val);
  free(bufs);
  if(armci_me==0) {
    printf("All connections between all procs tested: SUCCESS\n");
    fflush(stdout);
  }
}

int PARMCI_Init()
{
    char *uval;
#if defined(THREAD_SAFE)
    int th_idx;
#endif

    if(_armci_initialized>0) return 0;
    dassertp(1,sizeof(armci_ireq_t) <= sizeof(armci_hdl_t),
	     ("nb handle sizes: internal(%d) should be <= external(%d)\n",
	      sizeof(armci_ireq_t), sizeof(armci_hdl_t)));

    /* let's hope that the message passing environment was initialized outside
     * of ARMCI such that passing NULL for argc/argv here is okay */
    armci_msg_init(NULL, NULL);

#ifdef MPI_SPAWN
    if(!_armci_initialized_args)
       armci_die("ARMCI is built w/ ARMCI_NETWORK=MPI-SPAWN. For this network "
                 "setting, ARMCI must be initialized with PARMCI_Init_args() "
                 "instead of PARMCI_Init(). Please replace PARMCI_Init() "
                 " with PARMCI_Init_args(&argc, &argv) as in the API docs", 0L);
#endif
#if defined(MPI_MT) || defined(DCMF)
    {
        int provided;
        MPI_Query_thread(&provided);
        if (provided == MPI_THREAD_SINGLE) {
            armci_die("ARMCI is built w/ ARMCI_NETWORK=MPI_MT but the "
                    "provided MPI threading level is MPI_THREAD_SINGLE "
                    " not MPI_THREAD_MULTIPLE", 1);
        }
        else if (provided == MPI_THREAD_FUNNELED) {
            armci_die("ARMCI is built w/ ARMCI_NETWORK=MPI_MT but the "
                    "provided MPI threading level is MPI_THREAD_FUNNELED "
                    " not MPI_THREAD_MULTIPLE", 1);
        }
        else if (provided == MPI_THREAD_SERIALIZED) {
            armci_die("ARMCI is built w/ ARMCI_NETWORK=MPI_MT but the "
                    "provided MPI threading level is MPI_THREAD_SERIALIZED "
                    " not MPI_THREAD_MULTIPLE", 1);
        }
        else if (provided == MPI_THREAD_MULTIPLE) {
        }
    }
#endif
    
#ifdef BGML
    BGML_Messager_Init();
    BG1S_Configuration_t config;
    config=BG1S_Configure(NULL);
    config.consistency= BG1S_ConsistencyModel_Weak;
    BG1S_Configure(&config);

    unsigned long long available = BGML_Messager_available();
    if (available & BGML_MESSAGER_GI)
      bgml_barrier = (BGML_Barrier) BGGI_Barrier;
    else
      bgml_barrier = (BGML_Barrier) BGTr_Barrier;
#endif
#ifdef ARMCIX
    ARMCIX_Init ();
#endif
    armci_nproc = armci_msg_nproc();
    armci_me = armci_msg_me();
    armci_usr_tid = THREAD_ID_SELF(); /*remember the main user thread id */

#if defined(THREAD_SAFE)
    armci_init_threads();
    th_idx = ARMCI_THREAD_IDX;
    if (th_idx)
        printf("WARNING: PARMCI_Init is called from thread %d, should be 0\n",th_idx);
#endif

#ifdef _CRAYMPP
    cmpl_proc=-1;
#endif
#ifdef LAPI
#   ifdef AIX
    {
       char *tmp1 = getenv("RT_GRQ"), *tmp2 = getenv("AIXTHREAD_SCOPE");
       if(tmp1 == NULL || strcmp((const char *)tmp1,"ON")) 
	  armci_die("Armci_Init: environment variable RT_GRQ not set. It should be set as RT_GRQ=ON, to restore original thread scheduling LAPI relies upon",0);
       if(tmp2 == NULL || strcmp((const char *)tmp2,"S")) 
	  armci_die("Armci_Init: environment variable AIXTHREAD_SCOPE=S should be set to assure correct operation of LAPI", 0);
    }
#   endif
    armci_init_lapi();
#endif

#ifdef PORTALS
    armci_init_portals();
    shmem_init();
#endif

#ifdef CRAY_SHMEM
    shmem_init();
#endif

    armci_init_clusinfo();

#ifdef MPI
    armci_group_init();
#endif
    armci_krmalloc_init_localmem();

#ifndef BLRTS
    /* trap signals to cleanup ARMCI system resources in case of crash */
    if(armci_me == armci_master) {
        ARMCI_ParentTrapSignals();
    }
    ARMCI_ChildrenTrapSignals();
#endif

#if defined(SYSV) || defined(WIN32) || defined(MMAP)
    /* init shared/K&R memory */
    if(ARMCI_Uses_shm() ) {
#      ifdef SGIALTIX
          armci_altix_shm_init();
#      else
          armci_shmem_init();
#      endif
    }

#   if defined(QUADRICS) && !defined(NO_SHM)
       if(armci_me == armci_master)armci_check_shmmax();
#   endif
#endif

#ifdef REGION_ALLOC
       {
	 void* test_ptr_arr = malloc(sizeof(void *)*MAX_PROC);
	 dassert(1,test_ptr_arr);
	 PARMCI_Malloc(test_ptr_arr,256*1024*1024);
       PARMCI_Free(test_ptr_arr[armci_me]);
       free(test_ptr_arr);
       }
#endif

#ifdef MULTI_CTX
    /* this is a hack for the Elan-3 multi-tiled memory (qsnetlibs v 1.4.10) 
     * we need to allocate and then free memory to satisfy libelan requirements
     * for symmetric memory addresses
     */ 
    if(armci_nclus >1){ 
       int segments, segsize, seg;
       void **addr;
       armci_nattach_preallocate_info(&segments, &segsize);

       segsize -= 1024*1024; /* leave some for the K&RM headers */
       if(armci_me!=armci_master)segsize=0; /* only one allocates mem on node*/

       addr = (void*) malloc(segments*armci_nproc*sizeof(void*));
       if(!addr)armci_die("armci_init:addr malloc failed",segments*armci_nproc);

       for(seg=0; seg< segments; seg++) /* allocate segments */
          if(PARMCI_Malloc(addr+armci_nproc*seg,segsize))
             armci_die("problem in Elan-3 mem preallocation",seg);
       
       for(seg=0; seg< segments; seg++) /* return to free pool */
         if(armci_me==armci_master)
           if(PARMCI_Free(*(addr+armci_nproc*seg+armci_me)))
              armci_die("problem in Elan-3 mem preallocation - free stage",seg);
       free(addr);

#if 0
       if(armci_me==armci_master){
          printf("%d:preallocated %d segments %d bytes each\n",armci_me,
                 segments, segsize); fflush(stdout);
       }
#endif

    }
#endif

    /* allocate locks: we need to do it before server is started */
    armci_allocate_locks();
    armci_init_fence();

#if ARMCI_ENABLE_GPC_CALLS
    gpc_init_signals();
#endif

#ifdef ALLOW_PIN
    armci_prot_switch_fence = malloc(sizeof(int*)*armci_nproc);
    armci_prot_switch_preproc = -1;
    armci_prot_switch_preop = -1;
#endif

    /* NOTE: FOR PROCESS-BASED DATA SERVER WE CANNOT call PARMCI_Malloc yet */

#   if defined(DATA_SERVER) 
       if(armci_nclus >1) 
           armci_start_server();
#   endif
#if defined(GM) || defined(VAPI) || defined(PORTALS) || (defined(LAPI) && defined(LAPI_RDMA))
    /* initialize registration of memory */
    armci_region_init();
#endif

    armci_msg_barrier();
    armci_init_memlock(); /* allocate data struct for locking memory areas */
#if !defined(GM) 
    armci_notify_init();
#endif
    armci_msg_barrier();
    armci_msg_gop_init();

    _armci_initialized=1;
#ifdef ENABLE_CHECKPOINT
    armci_init_checkpoint(armci_ft_spare_procs);
#endif

#ifdef MPI_MT
    _armci_test_connections();
#else
    uval = getenv("ARMCI_TEST_CONNECTIONS"); 
    if(uval!=NULL) {
      _armci_test_connections();
    }
#endif

#if MSG_COMMS_TCGMSGMPI
    install_nxtval(NULL, NULL);
#endif
    return 0;
}

/* ARMCI Finalize is called multiple times, if both GA and TCGMSG are used
 * */
void PARMCI_Finalize()
{
    if(_armci_initialized <= 0 ) {
        return;
    }

    _armci_initialized = 0;
    _armci_terminating =1;
    _armci_initialized_args=0;
    _armci_argc = NULL;
    _armci_argv = NULL;

    armci_msg_barrier();
    if(armci_me==armci_master) 
        ARMCI_ParentRestoreSignals();

#if defined(DATA_SERVER)
    if(armci_nclus >1){
        armci_wait_for_server();
        armci_msg_barrier(); 
    }
#endif

#ifdef PORTALS
    armci_fini_portals();
#endif
#ifdef LAPI
    armci_term_lapi();
#endif
#ifdef ALLOW_PIN
    free(armci_prot_switch_fence);
#endif
    armci_msg_gop_finalize();
    ARMCI_Cleanup();
    armci_msg_barrier();
#ifdef MPI
    armci_group_finalize();
#endif
#ifdef ARMCIX
    ARMCIX_Finalize ();
#endif
#ifdef MPI
    MPI_Comm_free(&ARMCI_COMM_WORLD); /*SK: free at last*/
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
       while(*p != val)
            if((++count)<maxspin) armci_util_spin(count,(int *)p);
            else{
               cpu_yield();
               count =0;
#if defined(MACX) && defined(__ppc__) && defined(__GNUC__)
               __asm__ __volatile__ ("sync" ::: "memory");
#endif
            }
}
  

int ARMCI_Same_node(int proc)
{
   int direct = SAMECLUSNODE(proc);
   return direct;
}

unsigned int _armci_get_next_tag(){
  static unsigned int _armci_nb_tag=0;
  unsigned int rval;
  THREAD_LOCK(armci_user_threads.lock);
  rval = ++_armci_nb_tag;
  THREAD_UNLOCK(armci_user_threads.lock);
  return rval;
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
   armci_notify_t *pnotify = _armci_notify_arr[armci_me]+proc;
   pnotify->sent++;
# ifdef MEM_FENCE
   if(SAMECLUSNODE(proc)) MEM_FENCE;
# endif
#ifdef OPENIB
   /* IB will optimze a simple Put by using RDMA.  This can bypass non-RDMA
    * Puts and lead to incorrect behavour.  Avoid that by using PutV, which
    * presently does not optimize to RDMA.
    * This workaround is sub-optimal for two reasons:
    * 1. This adds more overhead when there is may be no need.
    * 2. There is no guarantee that PutV will always be un-optimized.
    */
   void *sp = &pnotify->sent;
   void *dp = &(_armci_notify_arr[proc]+armci_me)->received;
   armci_giov_t gv;

   gv.src_ptr_array = &sp;
   gv.dst_ptr_array = &dp;
   gv.ptr_array_len = 1;
   gv.bytes = sizeof(pnotify->sent);

   PARMCI_PutV(&gv, 1, proc);
#else
   PARMCI_Put(&pnotify->sent,&(_armci_notify_arr[proc]+armci_me)->received, 
             sizeof(pnotify->sent),proc);
#endif /* OPENIB */
   return(pnotify->sent);
}


/* blocks until received count becomes >= waited count
 *  return received count and store waited count in *pval
 */
int parmci_notify_wait(int proc,int *pval)
{
  int retval;
     long loop=0;
     armci_notify_t *pnotify = _armci_notify_arr[armci_me]+proc;
     pnotify->waited++;
     while( pnotify->waited > pnotify->received) {
         if(++loop == 1000) { loop=0;cpu_yield(); }
         armci_util_spin(loop, pnotify);
     }
     *pval = pnotify->waited;
     retval=pnotify->received;
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

#ifdef CHANGE_SERVER_AFFINITY
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

static int in_error_cleanup=0;

void derr_printf(const char *format, ...) {
    
  if(!in_error_cleanup) {
#ifdef SYSV
    if((!AR_caught_sigterm && !AR_caught_sigint) || armci_me==0) 
#endif
    {
      va_list ap;
      va_start(ap, format);
      vprintf(format, ap);
      va_end(ap);
    }
  }
}


int dassertp_fail(const char *cond_string, const char *file, 
		  const char *func, unsigned int line, int code) {
  if(!in_error_cleanup) {
    /* JAD 02/23/2012 for applications, an exit/error code of 0 indicates
     * success, it is therefore wrong to call dassertp_fail with a zero value */
    if (0 == code) {
        code = -1;
    }
    in_error_cleanup=1;
#ifdef SYSV
    if((!AR_caught_sigterm && !AR_caught_sigint) || armci_me==0)
#endif
    {
      printf("(rank:%d hostname:%s pid:%d):ARMCI DASSERT fail. %s:%s():%d cond:%s\n",
	     armci_me,armci_clus_info[armci_clus_me].hostname, 
	     getpid(), file,func,line,cond_string);
#if defined(PRINT_BT)
      backtrace_symbols_fd(bt, backtrace(bt, 100), 2);
#endif
    }
    armci_abort(code);
  }
  return code;
}
