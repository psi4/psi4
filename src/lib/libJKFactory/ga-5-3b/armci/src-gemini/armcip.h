/* $Id: armcip.h,v 1.82.2.9 2007-08-29 17:32:31 manoj Exp $ */
/* armci private header file */
#ifndef _ARMCI_P_H

#define _ARMCI_P_H
#include <stdlib.h>
#include "armci.h"
#include "message.h"
// #include "code_options.h"
#if 0
#define ARMCI_PR_DBG(__ARMCI_ST,__ARMCI_NU) \
        printf("\n%d:%s:%d:%s:%s:%d",armci_me,__FILE__,__LINE__,__FUNCTION__,__ARMCI_ST,__ARMCI_NU);fflush(stdout)
#define ARMCI_PR_SDBG(__ARMCI_ST,__ARMCI_NU) \
        printf("\n(%d):%s:%d:%s:%s:%d",armci_me,__FILE__,__LINE__,__FUNCTION__,__ARMCI_ST,__ARMCI_NU);fflush(stdout)
#else
#define ARMCI_PR_DBG(__ARMCI_ST,__ARMCI_NU)
#define ARMCI_PR_SDBG(__ARMCI_ST,__ARMCI_NU)
#endif

#ifdef LIBONESIDED
#include "armci-onesided.h"
#endif

#define DATA_SERVER
#define SERVER_THREAD
/*#define ARMCI_CHECK_STATE*/

#define ARMCI_STAMP 11214
#define ARMCI_TAIL 31121
#ifdef QUADRICS
#include <elan/elan.h>
#ifdef QSNETLIBS_VERSION_CODE
#ifndef DECOSF
#  define ELAN_ACC
#  define PENDING_OPER(x) ARMCI_ACC_INT
#endif

#  if QSNETLIBS_VERSION_CODE > QSNETLIBS_VERSION(1,5,0)
#     define LIBELAN_ATOMICS
#  endif

#endif
extern void armci_elan_fence(int p);
#endif

/* we got problems on IA64/Linux64 with Elan if inlining is used */
#if defined(__GNUC__) && !defined(QUADRICS)
#   define INLINE inline
#else
#   define INLINE
#endif

#ifdef WIN32
#include <windows.h>
#define sleep(x) Sleep(100*(x))
#else
#include <unistd.h>
#endif

#if (defined(SYSV) || defined(WIN32)|| defined(MMAP)) && !defined(NO_SHM) && !defined(HITACHI) && !defined(CATAMOUNT)
#define CLUSTER

#ifdef SERVER_THREAD
#  define SERVER_NODE(c) (armci_clus_info[(c)].master);
#  define NODE_SERVER(c) (c);
#else
#  define SOFFSET -1000000
#  define SERVER_NODE(c) ((int)(SOFFSET -armci_clus_info[(c)].master));
#  define NODE_SERVER(c) ((int)(SOFFSET - c))
#endif

#endif


/*\GPC call stuff
\*/
typedef struct {
  int hndl, hlen, dlen;
  void *hdr, *data;
}gpc_send_t;


/*\ Stuff for non-blocking API
\*/ 
#define NB_MULTI -1 /*more than one armci buffer(buffers.c) used for nbcall*/
#define NB_NONE  -2 /*no armci buffer(buffers.c) used for nbcall*/
extern unsigned int _armci_get_next_tag();
#define GET_NEXT_NBTAG _armci_get_next_tag
#define ARMCI_MAX_IMPLICIT 15

typedef struct{
  int len;
  int last;
  void *exthdr;
} ext_header_t;

typedef struct{
int val;
void *ptr;
} armci_flag_t;


#if defined(LAPI) || defined(PTHREADS) || defined(POSIX_THREADS)
# include <pthread.h>
  typedef pthread_t thread_id_t;
# define  THREAD_ID_SELF pthread_self
#elif defined(WIN32)
# include <windows.h>
  typedef DWORD thread_id_t;
# define  THREAD_ID_SELF GetCurrentThreadId  
#else
  typedef int thread_id_t;
# define  THREAD_ID_SELF() 1  
#endif

extern thread_id_t armci_usr_tid;
#ifdef SERVER_THREAD
#  define SERVER_CONTEXT (armci_usr_tid != THREAD_ID_SELF())
#else
#  define SERVER_CONTEXT (armci_me<0)
#endif

#if defined(LAPI) || defined(CLUSTER) || defined(CRAY) \
        || defined(CRAY_SHMEM) || defined(BGML) || defined(DCMF)
#  include "request.h"
#endif

/* ------------------------ ARMCI threads support ------------------------- */
#define ARMCI_THREADS_LIMIT 32

#include "utils.h"
#if defined(THREAD_SAFE)
typedef struct {
    int max;                /* max # of threads per proc */
    int avail;              /* next available position */
    thread_id_t *ids;       /* list of threads' ids */
    thread_lock_t lock;     /* general case lock */
    thread_lock_t buf_lock; /* lock for buffer access */
    thread_lock_t net_lock; /* lock for network accees */
} armci_user_threads_t;

extern armci_user_threads_t armci_user_threads;

extern void armci_init_threads();
extern void armci_finalize_threads();
extern int armci_thread_idx();
extern INLINE int armci_register_thread(thread_id_t id);

#define ARMCI_THREAD_IDX armci_thread_idx() /* needs to be optimized */

#else
#   define ARMCI_THREAD_IDX 0
#endif

/* ------------------------------------------------------------------------ */

/* min amount of data in strided request to be sent in single TCP/IP message*/
#if defined(SOCKETS) || defined(MPI_SPAWN_ZEROCOPY)
#  define TCP_PAYLOAD 128 
#  define LONG_GET_THRESHOLD  TCP_PAYLOAD  
#  define LONG_GET_THRESHOLD_STRIDED LONG_GET_THRESHOLD
#  define LONG_PUT_THRESHOLD 128  
#endif

#ifdef WIN32
#  define bzero(a,len){\
     int _i;\
     char *_c = (char*)(a);\
     for(_i=0; _i< (int)(len); _i++)_c[_i]=(char)0;\
   }
#  define bcopy(a,b,len) memcpy(b,a,len)
#else
# include <strings.h>
#endif

/*#define ACC_COPY*/
#if defined(CRAY_T3E) || defined(FUJITSU)\
       || defined(HITACHI) || (defined(QUADRICS) && !defined(ELAN_ACC))
#define ACC_COPY
#endif

#ifndef FATR
# ifdef WIN32
#   define FATR __stdcall
# else
#   define FATR 
# endif
#endif

#define MAX_PROC 8096
#define MAX_STRIDE_LEVEL ARMCI_MAX_STRIDE_LEVEL

/* msg tag ARMCI uses in collective ops */
#define ARMCI_TAG 30000

#ifndef EXTRA_MSG_BUFLEN_DBL
#   define RESERVED_BUFLEN ((sizeof(request_header_t)>>3)+3*MAX_STRIDE_LEVEL)
#else
#   define RESERVED_BUFLEN ((sizeof(request_header_t)>>3)+3*MAX_STRIDE_LEVEL +\
                           EXTRA_MSG_BUFLEN_DBL)
#endif

#if defined(HITACHI)
#  define BUFSIZE  ((0x50000) * sizeof(double))
#else   
   /* packing algorithm for double complex numbers requires even number */
#  ifdef MSG_BUFLEN_DBL
#    define BUFSIZE_DBL (MSG_BUFLEN_DBL - RESERVED_BUFLEN)
#  else
#    define BUFSIZE_DBL 32768
#  endif
#  define BUFSIZE  (BUFSIZE_DBL * sizeof(double))
#endif

/* note opcodes must be lower than ARMCI_ACC_OFF !!! */
#define PUT 1
#define GET 2
#define RMW 3
#define LOCK   4
#define UNLOCK 5
#define ACK 6 
#define STATE 11214

/* must fit in two bits, see msginfo->format in request.h */
#define STRIDED 1
#define VECTOR  2

extern  int armci_me, armci_nproc;
extern int _armci_initialized;
#ifdef HITACHI
   extern int sr8k_server_ready;
   extern  double *armci_internal_buffer;
#else
   extern  double armci_internal_buffer[BUFSIZE_DBL];
#endif
extern int armci_getbufsize();
extern void armci_shmem_init();
extern void armci_krmalloc_init_localmem();
extern void armci_die(char *msg, int code);
extern void armci_die2(char *msg, int code1, int code2);
extern void armci_write_strided(void *ptr, int stride_levels, 
                                int stride_arr[], int count[], char *buf);
extern void armci_read_strided(void *ptr, int stride_levels, 
                               int stride_arr[], int count[], char *buf);
extern int armci_op_strided(int op, void* scale, int proc,void *src_ptr, 
			int src_stride_arr[],  void* dst_ptr, 
                        int dst_stride_arr[], int count[],  
                        int stride_levels, int lockit,armci_ihdl_t nb_handle);
extern int armci_copy_vector(int op, /* operation code */
                armci_giov_t darr[], /* descriptor array */
                int len,  /* length of descriptor array */
                int proc /* remote process(or) ID */
              );

extern int armci_acc_vector(int op, /* operation code */
                void *scale,        /* scale factor */
                armci_giov_t darr[],/* descriptor array */
                int len,  /* length of descriptor array */
                int proc  /* remote process(or) ID */
              );

extern int armci_pack_strided(int op, void* scale, int proc,
                       void *src_ptr, int src_stride_arr[],
                       void* dst_ptr, int dst_stride_arr[],
                       int count[], int stride_levels, ext_header_t *hdr,
                       int fit_level, int nb, int last,armci_ihdl_t nb_handle);

extern int armci_pack_vector(int op, void *scale, 
                    armci_giov_t darr[],int len,int proc,armci_ihdl_t nb_handle);

extern void armci_lockmem(void *pstart, void* pend, int proc);
extern void armci_unlockmem(int proc);

extern int armci_acc_copy_strided(int optype, void* scale, int proc,
                                  void* src_ptr, int src_stride_arr[],  
		                  void* dst_ptr, int dst_stride_arr[], 
                                  int count[], int stride_levels);

extern void armci_vector_to_buf(armci_giov_t darr[], int len, void* buf);
extern void armci_vector_from_buf(armci_giov_t darr[], int len, void* buf);
extern void armci_init_fence();

#ifdef SOCKETS
#ifdef SERVER_THREAD
  extern void armci_create_server_thread ( void* (* func)(void*) );
  extern void armci_terminate_server_thread();
#else
  extern void armci_create_server_process ( void* (* func)(void*) );
  extern void armci_wait_server_process();
  extern void RestoreSigChldDfl();
#endif
#endif

#ifdef MPI_SPAWN
  extern void armci_create_server_MPIprocess ();
#endif


#define ARMCI_MAX(a,b) (((a)>(b))?(a):(b))
#define ARMCI_MIN(a,b) (((a)<(b))?(a):(b))
#define ARMCI_ABS(a)   (((a) >= 0) ? (a) : (-(a)))
#define ARMCI_ACC(op)  ((((int)(op))-ARMCI_ACC_INT)>=0)


#ifdef CLUSTER
   extern char *_armci_fence_arr;
#  define FENCE_ARR(p_) (_armci_fence_arr[p_])

#  define SAMECLUSNODE(p)\
     ( ((p) <= armci_clus_last) && ((p) >= armci_clus_first) )
#elif defined(__crayx1)
#  define SAMECLUSNODE(p) 1
#elif defined(ARMCIX)
#  define SAMECLUSNODE(p) 0
#else
#  define SAMECLUSNODE(p) ((p)==armci_me) 
#endif


#if defined(LAPI) || defined(ELAN_ACC)
#  define ORDER(op,proc)\
        if( proc == armci_me || ( ARMCI_ACC(op) && ARMCI_ACC(PENDING_OPER(proc))) );\
        else  FENCE_NODE(proc)
#  define UPDATE_FENCE_INFO(proc_)
#elif defined(CLUSTER) && !defined(QUADRICS) && !defined(HITACHI)\
        && !defined(CRAY_SHMEM)
#  define ORDER(op_,proc_)\
          if(!SAMECLUSNODE(proc_) && op_ != GET )FENCE_ARR(proc_)=1
#  define UPDATE_FENCE_INFO(proc_) if(!SAMECLUSNODE(proc_))FENCE_ARR(proc_)=1
#else
#  if defined(GM) && defined(ACK_FENCE) 
#   define ORDER(op,proc)
#  else
#   define ORDER(op,proc) if(proc != armci_me) FENCE_NODE(proc)
#  endif 
#  define UPDATE_FENCE_INFO(proc_)
#endif
        
typedef struct {
    int  ptr_array_len;
    int bytes;
    void **ptr_array;
} armci_riov_t;

/*\ consider up to HOSTNAME_LEN characters in host name 
 *  we can truncate names of the SP nodes since it is not used
 *  to establish socket communication like on the networks of workstations
 *  SP node names must be distinct within first HOSTNAME_LEN characters
\*/
#if defined(LAPI) && defined(AIX)
#  define HOSTNAME_TRUNCATE 
#  define HOSTNAME_LEN 12
#else
#  define HOSTNAME_LEN 64
#endif

typedef struct {
  int master;
  int nslave;
  char hostname[HOSTNAME_LEN];
} armci_clus_t;

extern armci_clus_t *armci_clus_info;
extern int armci_nclus, armci_clus_me, armci_master;
extern int armci_clus_first, armci_clus_last;
extern int armci_clus_id(int p);
extern void armci_init_clusinfo();
extern void armci_set_mem_offset(void *ptr);
extern int _armci_terminating;
extern void armci_acc_2D(int op, void* scale, int proc, void *src_ptr, 
                         void *dst_ptr, int bytes, int cols, int src_stride, 
                         int dst_stride, int lockit);
extern void armci_lockmem_scatter(void *ptr_array[], int len, int bytes, int p);
extern void armci_generic_rmw(int op, void *ploc, void *prem, int extra, int p);
extern unsigned long armci_max_region();
extern void armci_dispatch_strided(void *ptr, int stride_arr[], int count[],
                            int strides, int fit_level, int nb, int bufsize, 
                            void (*fun)(void*,int*,int*,int,void*), void *arg);
extern void armci_msg_gop_init();
extern void armci_util_spin(int n, void *notused);

#if defined(SYSV) || defined(WIN32)
extern void armci_shmem_init();
extern void armci_set_shmem_limit_per_core(unsigned long shmemlimit);
extern void armci_set_shmem_limit_per_node(int nslaves);
extern void armci_set_shmem_limit(unsigned long shmemlimit);
#endif

#define ALIGN_PTR_LONG(type, x) if( ((long)(x)) % sizeof(long)) { long _y = (long)(x);\
              if(sizeof(long)==8){_y>>=3; _y<<=3; }\
                            else { _y>>=2; _y<<=2; }\
              _y += sizeof(long); (x) = (type*)_y; }

#define SIXTYFOUR 64
#define ALIGN64ADD(buf) (SIXTYFOUR-(((ssize_t)(buf))%SIXTYFOUR))
#define ALIGNLONGADD(buf) ((((ssize_t)(buf))%sizeof(long))?(sizeof(long)-(((ssize_t)(buf))%sizeof(long))):0)

#define SET   1
#define UNSET 0

extern int armci_agg_save_strided_descriptor(void *src_ptr, 
					     int src_stride_arr[],
					     void* dst_ptr, 
					     int dst_stride_arr[],
					     int count[], 
					     int stride_levels, int proc,
					     int op, armci_ihdl_t nb_handle);
     
extern int armci_agg_save_giov_descriptor(armci_giov_t darr[], int len, 
					   int proc, int op, 
					   armci_ihdl_t nb_handle);

extern int armci_agg_save_descriptor(void *src, void *dst, int bytes, 
				      int proc, int op, int is_registered_put,
				      armci_ihdl_t nb_handle);

extern void armci_agg_complete(armci_ihdl_t nb_handle, int condition);

extern armci_ihdl_t armci_set_implicit_handle (int op, int proc);

extern int armci_getnumcpus(void);
extern long armci_util_long_getval(long* p);
extern int armci_util_int_getval(int* p);
extern void armci_region_register_shm(void *start, long size);
extern void armci_region_register_loc(void *start, long size);
extern void armci_region_clus_record(int node, void *start, long size);
extern void armci_region_init();
extern int armci_region_clus_found(int node, void *start, int size);
extern int armci_region_loc_found(void *start, int size);
extern int armci_region_both_found(void *loc, void *rem, int size, int node);
#ifdef REGIONS_REQUIRE_MEMHDL
extern int get_armci_region_local_hndl(void *loc, int node, ARMCI_MEMHDL_T **loc_memhdl);
#endif
extern void armci_region_exchange(void *start, long size);
extern void cpu_yield();

#ifdef ALLOW_PIN
extern void armci_global_region_exchange(void *, long); 
#endif


/* -------------------- ARMCI Groups ---------------------- */
/* data structure that caches a group's attribute */
#ifdef BGML
#define   PCLASS 3
#endif
#ifdef MPI
typedef int ARMCI_Datatype;

extern int ATTR_KEY; /* attribute key */

/* #define ARMCI_GROUP /\*Generic ARMCI implementation*\/ */
 
typedef struct {
  armci_clus_t *grp_clus_info;
  int grp_me;              /* my group id */
  int grp_nclus;           /* number of cluster nodes */
  int grp_clus_me;         /* my cluster node id */
  int mem_offset;          /* memory offset */
#ifdef ARMCI_GROUP
  int nproc;               /* #procs in this group*/
  int *proc_list;          /* Ids of procs in this group
			      (w.r.t. MPI_COMM_WORLD)*/
#endif
}armci_grp_attr_t;
 
#include "mpi.h"
 
/**dup of MPI_COMM_WORLD for internal MPI communication*/
extern MPI_Comm ARMCI_COMM_WORLD;

#ifdef PORTALS
#include "portals.h"
#endif

typedef MPI_Comm ARMCI_Comm;
typedef struct {
#ifndef ARMCI_GROUP
  MPI_Comm icomm;
  MPI_Group igroup;
#endif
  armci_grp_attr_t grp_attr;
}ARMCI_iGroup;
 
armci_grp_attr_t *ARMCI_Group_getattr(ARMCI_Group *grp);

extern void armci_group_init();
extern void armci_group_finalize();
extern ARMCI_iGroup* armci_get_igroup_from_group(ARMCI_Group *group);
 
#endif /* ifdef MPI */ 
/* -------------------------------------------------------- */

/* ------------ ARMCI Chekcpointing/Recovery -------------- */
#ifdef DO_CKPT
extern int armci_init_checkpoint();
extern void armci_create_ckptds(armci_ckpt_ds_t *ckptds, int count);
extern int armci_icheckpoint_init(char *filename, ARMCI_Group *grp,
                                  int savestack, int saveheap,
                                  armci_ckpt_ds_t *ckptds);
extern int armci_icheckpoint(int rid);
extern int armci_irecover(int rid,int iamreplacement);
extern void armci_icheckpoint_finalize(int rid);


#endif /* ifdef DO_CKPT */
/* -------------------------------------------------------- */

#endif
