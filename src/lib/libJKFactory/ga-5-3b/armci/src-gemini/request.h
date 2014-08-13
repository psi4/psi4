#ifndef _REQUEST_H_
#define _REQUEST_H_


/********  client buffer managment ops ****************************/
extern void  _armci_buf_init();
extern char* _armci_buf_get(int size, int operation, int to);
extern void  _armci_buf_release(void *buf);
extern int   _armci_buf_to_index(void *buf);
extern char* _armci_buf_ptr_from_id(int id);
extern void  _armci_buf_ensure_one_outstanding_op_per_node(void *buf, int node);
#if defined(SERV_QUEUE)
extern void  _armci_buf_ensure_pend_outstanding_op_per_node(void *buf, int node);
#endif
extern void _armci_buf_complete_nb_request(int bufid,unsigned int tag, int *retcode);
extern void _armci_buf_test_nb_request(int bufid,unsigned int tag, int *retcode);
extern void _armci_buf_set_tag(void *bufptr,unsigned int tag,short int protocol);
extern void _armci_buf_clear_all();
extern void x_buf_send_complete(void *);

extern INLINE char *_armci_buf_get_clear_busy(int size, int operation, int to);
extern INLINE void _armci_buf_set_busy(void *buf, int state);
extern INLINE void _armci_buf_set_busy_idx(int tbl_idx, int state);
extern INLINE int  _armci_buf_cmpld(int bufid);
extern INLINE void _armci_buf_set_cmpld(void *buf, int state);
extern INLINE void _armci_buf_set_cmpld_idx(int idx, int state);

#ifdef LAPI
#  include "lapidefs.h"
#elif LIBONESIDED
   typedef armci_onesided_msg_tag_t msg_tag_t;
#elif PORTALS
#  include "armci_portals.h"
#elif defined(GM)
#  include "myrinet.h"
#elif defined(DOELAN4)
#  include "elandefs.h"
#elif defined(QUADRICS)
#  include <elan/elan.h>
   typedef void* msg_tag_t; 
#  ifdef _ELAN_PUTGET_H
#    define NB_CMPL_T ELAN_EVENT*
#  endif
#elif defined(VIA)
#  include "via.h"
   typedef void* msg_tag_t;
#elif defined(VAPI)
#  include "armci-vapi.h"
#elif defined(SOCKETS)
#  include "sockets.h"
   typedef long msg_tag_t;
   typedef unsigned short msg_id_t;
#   define DTAG_ ((1<<(sizeof(msg_id_t)*8))-1)
#   define NB_SOCKETS_ /* define NB_SOCKETS to allow non-blocking path */
#elif defined(HITACHI)
#  include "sr8k.h"
#elif defined(BGML)
#  include "bgml.h"
#  include "bgmldefs.h"
#  define NB_CMPL_T BG1S_t  
    typedef long msg_tag_t;
#elif defined(MPI_SPAWN)
#  include "mpi2.h"
#  define MSG_BUFLEN_DBL 500000
   typedef long msg_tag_t;
#else
   typedef long msg_tag_t;
#endif

#ifndef CLEAR_HNDL_FIELD 
#   define CLEAR_HNDL_FIELD(_x) 
#endif

#define ACK_QUIT 0
#define QUIT 33
#define ATTACH 34
#define REGISTER 35
   
/*\ the internal request structure for non-blocking api. 
\*/
typedef struct{
   unsigned int tag;
   int bufid;
   int agg_flag;
   int op;
   int proc;

#ifdef NB_CMPL_T
   NB_CMPL_T cmpl_info;
#endif

   int onesided_direct;
   cos_desc_t comm_desc[MAX_OUTSTANDING_ONESIDED_GETS];
} armci_ireq_t;
/*\ the internal request structure for non-blocking api. 
\*/
typedef armci_ireq_t* armci_ihdl_t;
extern void armci_set_nbhandle_bufid(armci_ihdl_t nb_handle, char *buf, int val);
extern void set_nbhandle(armci_ihdl_t *nbh, armci_hdl_t *nb_handle,
                                int op, int proc);

typedef struct {
   int   to;            /* message recipient */
   int from;            /* message sender */
   int   operation;   /* operation code */
   int   format;      /* data format used */
   int   bytes;      /* number of bytes requested */
   int   datalen;       /* >0 in lapi means that data is included */
   int   ehlen;       /* size of extra header and the end of descr */
   int   dscrlen;    /* >0 in lapi means that descriptor is included */
   msg_tag_t tag;       /* message tag for response to this request, MUST BE LAST */
}request_header_t;


typedef struct _buf_ackresp{
   long val,valc;
   cos_request_t req;
   struct _buf_ackresp *next, *previous;
} _buf_ackresp_t;

/*******gpc call strctures*************/
#include <signal.h>
#define MAX_GPC_REQ 1
#define MAX_GPC_REPLY_LEN (64*1024)
#define MAX_GPC_SEND_LEN (64*1024)
#define GPC_COMPLETION_SIGNAL SIGUSR1

typedef struct {
  int hndl;
  int hlen, dlen;
  void *hdr, *data;
  int rhlen, rdlen;
  void *rhdr, *rdata;
} gpc_call_t;

typedef struct {
  int active;
/*    int zombie; */
  request_header_t  msginfo;
  gpc_call_t call;
  char send[MAX_GPC_SEND_LEN];
  char reply[MAX_GPC_REPLY_LEN];
} gpc_buf_t;

/*  gpc_buf_t *gpc_req; */
extern gpc_buf_t *gpc_req;

extern void block_pthread_signal(int signo);
extern void unblock_pthread_signal(int signo);

/*******structures copied from async.c for storing cmpl dscr for nb req*******/
#define UBUF_LEN 112

typedef struct {
  unsigned int tag;             /* request id*/
  _buf_ackresp_t ar;
  short int bufid;              /* communication buffer id */
  short int protocol;           /* what does this buf hold?*/
  union {                 
        void *dscrbuf;          /*in case dscr below is not enough, do a*/
        double pad;             /*malloc, save pointer in dscrbuf and use it*/
  }ptr;
  char dscr[UBUF_LEN];          /*place to store the dscr*/
}_buf_info_t;

#define BUF_INFO_T _buf_info_t
extern BUF_INFO_T *_armci_buf_to_bufinfo(void *buf);
#define BUF_TO_BUFINFO _armci_buf_to_bufinfo

void armci_complete_req_buf(BUF_INFO_T *info, void *buffer);
extern INLINE BUF_INFO_T *_armci_id_to_bufinfo(int bufid);

#ifndef MAX_BUFS
#define MAX_BUFS 8
#error "MAX_BUFS set to 8"
#endif

#ifndef MAX_SMALL_BUFS
#define MAX_SMALL_BUFS 16
#error "MAX_SMALL_BUFS set to 16"
#endif


/* tracks sockets used for receiving responces from data server (GET) */
typedef struct {
  int socks[MAX_BUFS+MAX_SMALL_BUFS]; /* sock # or -1 if not used */
  int ready[MAX_BUFS+MAX_SMALL_BUFS]; /* 1 - ready, 0 - not */
} active_socks_t;



/*valid values for the element protocol in BUF_INFO_T*/
#define SDSCR_IN_PLACE 1 /*indicated that strided descriptor is in place*/
#define VDSCR_IN_PLACE 2 /*indicated that vector descriptor is in place*/
#define VDSCR_IN_PTR   3 /*indicates that the vector descriptor in allocated 
                           and pointer stored in dscrbuf */
/****************************************************************************/

/* this effects: buf_ext_t, portalsEagerMessageSendSize, portals ds buffer size */
/* note: MSG_BUFLEN_DBL is being defined earlier in armci-portals.h */
#ifndef MSG_BUFLEN_DBL
# error "MSG_BUFLEN_DBL not yet defined"
# if defined(HITACHI)
#  define MSG_BUFLEN_DBL 0x50000
# else
#  ifdef PORTALS_USE_RENDEZ_VOUS
#    define MSG_BUFLEN_DBL 50000 /* for rendez-vous, this can go bigger i think */ 
#  else
#    define MSG_BUFLEN_DBL 8192  /* this is smaller when rendez-vous is off */
#  endif
# endif
#endif

#define MSG_BUFLEN  sizeof(double)*MSG_BUFLEN_DBL
extern  char* MessageRcvBuffer;
extern  char* MessageSndBuffer;

#ifdef LAPI
#  define GET_SEND_BUFFER_(_size)(MessageSndBuffer+sizeof(lapi_cmpl_t));\
          CLEAR_COUNTER(*((lapi_cmpl_t*)MessageSndBuffer));\
          SET_COUNTER(*((lapi_cmpl_t*)MessageSndBuffer),1);
#  define GET_SEND_BUFFER _armci_buf_get
#  define GA_SEND_REPLY armci_lapi_send
#else
#  ifdef SOCKETS
#    define GA_SEND_REPLY(tag, buf, len, p) armci_sock_send(p,buf,len)
#  else
#    define GA_SEND_REPLY(tag, buf, len, p)  
#  endif
#endif

#ifdef QUADRICS_
#  define GET_SEND_BUFFER(_size,_op,_to) MessageSndBuffer;\
                    while(((request_header_t*)MessageSndBuffer)->tag)\
                    armci_util_spin(100, MessageSndBuffer)
#  define FREE_SEND_BUFFER(_ptr) ((request_header_t*)MessageSndBuffer)->tag = (void*)0 
#endif

#ifndef GET_SEND_BUFFER
#  define GET_SEND_BUFFER(_size,_op,_to) MessageSndBuffer
#endif

#ifndef FREE_SEND_BUFFER
#define FREE_SEND_BUFFER(_ptr)  
#endif

#ifndef INIT_SENDBUF_INFO
#define INIT_SENDBUF_INFO(_hdl,_buf,_op,_proc) 
#endif

typedef struct {
           char *buf; char* buf_posted; int count; int proc; int op; int extra;
} buf_arg_t;

/*includes for SERVER_LOCK*/
#if defined(SERVER_THREAD) && !defined(VIA)
   extern void armci_rem_lock(int mutex, int proc, int *ticket);
   extern void armci_rem_unlock(int mutex, int proc, int ticket);
   extern void armci_unlock_waiting_process(msg_tag_t tag,int proc, int ticket);
#endif


#ifdef PIPE_BUFSIZE 
   extern void armcill_pipe_post_bufs(void *ptr, int stride_arr[], int count[],
                                      int strides, void* argvoid);
   extern void armcill_pipe_extract_data(void *ptr,int stride_arr[],int count[],
                                         int strides, void* argvoid);
   extern void armcill_pipe_send_chunk(void *data, int stride_arr[],int count[],
                                       int strides, void* argvoid);
#endif

extern void armci_send_strided(int proc, request_header_t *msginfo, char *bdata,
                         void *ptr, int strides, int stride_arr[], int count[],int tag);

extern void armci_rcv_hdlr(request_header_t* msginfo);

extern char *armci_rcv_data(int proc, request_header_t *msginfo, int rcvlen);
extern void armci_rcv_strided_data_bypass(int proc, request_header_t *msginfo,
                                          void *ptr, int stride_levels);
extern void armci_send_strided_data_bypass(int proc, request_header_t *msginfo,
            void *loc_buf, int msg_buflen, void *loc_ptr, int *loc_stride_arr,
            void *rem_ptr, int *rem_stride_arr, int *count, int stride_levels);

extern void armci_rcv_strided_data(int proc, request_header_t* msginfo, 
                  int datalen, void *ptr, int strides,int stride_arr[],int count[]);
extern void armci_send_strided_data(int proc,  request_header_t *msginfo, 
            char *bdata, void *ptr, int strides, int stride_arr[], int count[]);
extern void armci_send_req(int proc, request_header_t* msginfo, int len,int tag);
extern void armci_server_rmw(request_header_t* msginfo,void* ptr, void* pextra);
extern int armci_rem_vector(int op, void *scale, armci_giov_t darr[],int len,
                            int proc,int flag,armci_ihdl_t nb_handle);
extern int armci_rem_strided(int op, void* scale, int proc,
                       void *src_ptr, int src_stride_arr[],
                       void* dst_ptr, int dst_stride_arr[],
                       int count[], int stride_levels, 
                       ext_header_t *h, int lockit,armci_ihdl_t nb_handle);

extern void armci_rem_rmw(int op, void *ploc, void *prem, int extra, int proc);
extern void armci_rem_ack(int clus);
extern void armci_server(request_header_t *msginfo, char *dscr, char* buf, 
                         int buflen);
extern void armci_server_vector(request_header_t *msginfo,
                                char *dscr, char* buf, int buflen);
extern void *armci_server_ptr(int);
extern void armci_serv_attach_req(void *info, int ilen, long size,
                                  void* resp,int rlen);
extern void armci_server_lock(request_header_t *msginfo);
extern void armci_server_unlock(request_header_t *msginfo, char* dscr);
extern void armci_create_server_thread ( void* (* func)(void*) );
extern int armci_server_lock_mutex(int mutex, int proc, msg_tag_t tag);
extern void armci_send_data(request_header_t* msginfo, void *data);
extern int armci_server_unlock_mutex(int mutex, int p, int tkt, msg_tag_t* tag);
extern void armci_rcv_vector_data(int p, request_header_t* msginfo, armci_giov_t dr[], int len);

#if !defined(LAPI) 
extern void armci_wait_for_server();
extern void armci_start_server();
extern void armci_transport_cleanup();
extern int armci_send_req_msg(int proc, void *buf, int bytes,int tag);
extern void armci_WriteToDirect(int proc, request_header_t* msginfo, void *buf);
extern char *armci_ReadFromDirect(int proc, request_header_t *msginfo, int len);
extern void armci_init_connections();
extern void *armci_server_code(void *data);
extern void armci_rcv_req(void *mesg, void *phdr, void *pdescr, 
                          void *pdata, int *buflen);
extern void armci_client_connect_to_servers();
extern void armci_data_server(void *mesg);
extern void armci_server_initial_connection();
extern void armci_call_data_server();
#endif
#ifdef SOCKETS
extern void armci_ReadStridedFromDirect(int proc, request_header_t* msginfo,
                  void *ptr, int strides, int stride_arr[], int count[]);
extern void armci_WriteStridedToDirect(int proc, request_header_t* msginfo,
                         void *ptr, int strides, int stride_arr[], int count[]);
extern void armci_serv_quit();
extern int armci_send_req_msg_strided(int proc, request_header_t *msginfo,
                          char *ptr, int strides, int stride_arr[],int count[]);
extern void armci_server_goodbye(request_header_t* msginfo);
#endif
#ifdef MPI_SPAWN
extern void armci_serv_quit();
extern void armci_server_goodbye(request_header_t* msginfo);
#endif
#ifdef HITACHI
extern void armci_server_goodbye(request_header_t* msginfo);
extern void armci_serv_quit();
#endif
extern void armci_server_ipc(request_header_t* msginfo, void* descr,
                             void* buffer, int buflen);

#ifdef PIPE_BUFSIZE
extern void armci_pipe_prep_receive_strided(request_header_t *msginfo,char *buf,
                       int strides, int stride_arr[], int count[], int bufsize);
extern void armci_pipe_receive_strided(request_header_t* msginfo, void *ptr,
                                int stride_arr[], int count[], int strides);
extern void armci_pipe_send_req(int proc, void *buf, int bytes);
#endif

extern void armci_rcv_strided_data_bypass_both(int, request_header_t*,void*, int*, int);
extern int armci_rem_get(int proc, void *src_ptr, int src_stride_arr[],
                  void* dst_ptr, int dst_stride_arr[], int count[], int stride_levels,
                  armci_ihdl_t nb_handle,void *mhloc,void *mhrem);

#if defined(ALLOW_PIN) && defined(VAPI)
extern int armci_two_phase_send(int proc,void *src_ptr,int src_stride_arr[],
                    void *dst_ptr,int dst_stride_arr[],int count[],
                    int stride_levels,void ** context_ptr,armci_ihdl_t nbhandle,
                    ARMCI_MEMHDL_T *mhloc);
extern int armci_two_phase_get(int proc, void*src_ptr, int src_stride_arr[],
                    void*dst_ptr,int dst_stride_arr[], int count[],
                    int stride_levels, void**context_ptr,
                    armci_ihdl_t nbhandle, ARMCI_MEMHDL_T *mhloc);
#endif
#endif
