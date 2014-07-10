#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: openib.c,v 1.4.2.9 2007-10-18 06:08:03 d3h325 Exp $
 *
 * File organized as follows
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRINGS_H
#   include <strings.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#include <mpi.h>

#include "cbuf.h"
#include "armcip.h"
#include "copy.h"
#include "request.h"
#include "armci-vapi.h"
#include "iterator.h"
#define DEBUG_INIT 0
#define DEBUG_FINALIZE 0
#define DEBUG_SERVER 0
#define DEBUG_CLN 0
#define TIME_INIT 0
#  define VAPIDEV_NAME "InfiniHost0"
#  define INVAL_HNDL 0xFFFFFFFF
#define RNR_TIMER 12

/*Debug macros used to tune what is being tested -- mostly openib calls*/
#define DBG_INIT  1
#define DBG_POLL  1
#define DBG_ALL   1

#define QP_INACTIVE 5
#define QP_REQ_SENT 2
#define QP_ACK_RCVD 3
#define QP_ACTIVE 4

u_int32_t armci_max_num_sg_ent;
u_int32_t armci_max_qp_ous_swr;
u_int32_t armci_max_qp_ous_rwr;

typedef struct {
   struct ibv_qp *qp;
   uint32_t sqpnum;                /*we need to exchng qp nums,arr for that*/
   uint16_t lid;
   uint16_t state;
   void *next; 
} armci_connect_t;

armci_connect_t *CLN_con, *SRV_con;
static uint32_t *SRV_rqpnums, *CLN_rqpnums; /*relevant rqp num arrs, to connect to svr and client*/
static uint32_t *CLN_rqpnumtmpbuf=NULL; /*temporary buf used during connection setup*/
/*\
 * datastrucure for infinihost NIC
\*/
typedef struct {
  uint16_t *lid_arr;                /*we need to exchange lids, arr for that*/
  struct ibv_context *handle;       /*device context/handle*/
  int maxtransfersize;
  struct ibv_device_attr attr;      /*device properties*/
  struct ibv_port_attr hca_port;    /*mostly for getting lid*/
  uint8_t active_port;
  struct ibv_pd *ptag;              /*protection tag*/
  const char *vendor;
  struct ibv_cq *scq;               /*send completion queue*/
  struct ibv_cq *rcq;               /*recv completion queue*/
  struct ibv_comp_channel *sch;     /*send completion channel*/
  struct ibv_comp_channel *rch;     /*recv completion channel*/
  void *scq_cntx;                   /*send context for completion queue*/
  void *rcq_cntx;                   /*recv context for completion queue*/
  int scv;                          /*send completion vector*/
  int rcv;                          /*recv completion vector*/
} vapi_nic_t;

typedef struct {
  armci_vapi_memhndl_t *prem_handle; /*address server to store memory handle*/
  armci_vapi_memhndl_t handle;
}ack_t;

armci_vapi_memhndl_t *CLN_handle;
armci_vapi_memhndl_t serv_memhandle, client_memhandle;
armci_vapi_memhndl_t *handle_array;
armci_vapi_memhndl_t *pinned_handle;

static vapi_nic_t nic_arr[3];
static vapi_nic_t *SRV_nic= nic_arr;
static vapi_nic_t *CLN_nic= nic_arr+1;
static int armci_server_terminating;

#define NONE -1
static int armci_ack_proc=NONE;

static int armci_vapi_server_ready;
static int armci_vapi_server_stage1=0;
static int armci_vapi_client_stage1=0;
static int armci_vapi_server_stage2=0;
static int armci_vapi_client_ready;
int _s=-1,_c=-1;
static int armci_vapi_max_inline_size=-1;
#define CLIENT_STAMP 101
#define SERV_STAMP 99
#define MAX_PROC_INLINE_SIZE 2048

static char * client_tail;
static char * serv_tail;
static ack_t *SRV_ack;

#if defined(PEND_BUFS)
typedef immbuf_t vapibuf_t;
typedef pendbuf_t vapibuf_pend_t;
#else
typedef struct {
    struct ibv_recv_wr  dscr;
    struct ibv_sge      sg_entry;
    char buf[CBUF_DLEN];
} vapibuf_t;
#endif

typedef struct {
    struct ibv_send_wr  snd_dscr;
    struct ibv_sge      ssg_entry;
    struct ibv_recv_wr  rcv_dscr;
    struct ibv_sge      rsg_entry;
  char buf[VBUF_DLEN];
} vapibuf_ext_t;

typedef struct {
    struct ibv_send_wr  rmw_dscr;
    struct ibv_sge      rmw_entry;
} vapirmw_t;

unsigned int armci_use_odcm = 0;
unsigned int armci_use_lazy_break = 0;
unsigned int armci_use_apm = 0;
unsigned int armci_use_apm_test = 0;
unsigned int armci_use_srq = 0;     
unsigned int armci_use_snft = 0;    
unsigned int armci_use_affinity = 0;
  
unsigned int armci_srq_size = 4096; 
  
pthread_t armci_async_thread[4];

void async_thread_hca_events(void *ctx);
void async_thread_ud_events(void *ctx);
  
void init_apm_lock(void);
static struct ibv_srq *create_srq(vapi_nic_t *nic);

void setup_ud_channel(void);

void process_recv_completion_from_server(armci_ud_rank *h, cbuf *v);
void process_recv_completion_from_client(armci_ud_rank *h, cbuf *v);

struct ibv_srq      *CLN_srq_hndl;
struct ibv_srq      *SRV_srq_hndl;
void post_recv(void);
struct Remote_Buf
{
    char            **buf;
    uint32_t        *qp_num;
    uint16_t        *lid;
    uint32_t        *rkey;
};

struct HCA 
{
    struct ibv_device   *ib_dev;
    struct ibv_context  *context;
    struct ibv_pd       *pd;
    struct ibv_cq       *cq;
    struct ibv_srq      *srq_hndl;
    struct ibv_comp_channel *comp_channel;
};

struct RC_Conn
{
    struct ibv_qp       **qp;
    uint16_t            *lid;
    int                 *status;
    uint32_t            *qp_num;
    struct ibv_ah       **ud_ah;
};

struct Remote_Buf rbuf;

struct RC_Conn      conn;
    
struct HCA          hca;
    
void handle_network_fault(struct ibv_wc *pdscr);

int process_recv_completion_from_client_flag;

int total_active_conn_to_server, total_active_conn_to_client, total_breaks;

void check_state_of_ib_connection(int a, int b);

static vapibuf_t **serv_buf_arr;
#if !defined(PEND_BUFS)
/*These are typically used as spare buffers for communication. Since
  we do not wait on completion anymore, we need to ensure things work
  fine when these have in-flight messages. Disabled for now.*/
static vapibuf_t *spare_serv_buf, *spare_serv_bufptr;
static vapibuf_ext_t *serv_buf;
#endif

static vapirmw_t rmw[64];

static int *flag_arr; /* flag indicates its receiving scatter data */
#define SERV 2
#define CLN 1

#define MAX_DESCR 2
typedef struct {
    int avail;
    struct ibv_qp *qp;
    struct ibv_recv_wr *descr;
} descr_pool_t;

static int* _gtmparr;
static void* test_ptr;
static int test_stride_arr[1];
static int test_count[2];
static int test_stride_levels;
char *MessageRcvBuffer;

extern void armci_util_wait_int(volatile int *,int,int);
void armci_send_data_to_client(int proc, void *buf,int bytes,void *dbuf);
void armci_server_register_region(void *,long,ARMCI_MEMHDL_T *);
static descr_pool_t serv_descr_pool = {MAX_DESCR,NULL,NULL};
static descr_pool_t client_descr_pool = {MAX_DESCR,NULL,NULL};

/**Buffer (long[1] used to set msginfo->tag.ack_ptr in
   client-side. See usage in SERVER_SEND_ACK macro*/
static long *ack_buf;

#define GET_DATA_PTR(buf) (sizeof(request_header_t) + (char*)buf)

#define BUF_TO_SDESCR(buf) ((struct ibv_send_wr *)(&((armci_vapi_field_t *)((char *)(buf) - sizeof(armci_vapi_field_t)))->sdscr))

#define BUF_TO_RDESCR(buf) ((struct ibv_recv_wr *)(&((armci_vapi_field_t *)((char *)(buf) - sizeof(armci_vapi_field_t)))->rdscr))

#define BUF_TO_SSGLST(buf) ((struct ibv_sge *)(&((armci_vapi_field_t *)((char *)(buf) - sizeof(armci_vapi_field_t)))->ssg_entry))

#define BUF_TO_RSGLST(buf) ((struct ibv_sge *)(&((armci_vapi_field_t *)((char *)(buf) - sizeof(armci_vapi_field_t)))->rsg_entry))

#define BUF_TO_ECBUF(buf) (vapibuf_ext_t*)(((char*)buf) - (sizeof(struct ibv_send_wr)+sizeof(struct ibv_recv_wr)+2*sizeof(struct ibv_sge)))

#define SERVER_SEND_ACK(p) do {            \
    assert(*ack_buf == ARMCI_STAMP);       \
    assert((p)>=0);                        \
    armci_send_data_to_client((p),ack_buf, \
      sizeof(long),msginfo->tag.ack_ptr);  \
  } while(0)
/* #define SERVER_SEND_ACK(p) {assert(serv_buf!=NULL);assert(msginfo->from==(p));*((long *)serv_buf->buf)=ARMCI_STAMP;armci_send_data_to_client((p),serv_buf->buf,sizeof(long),msginfo->tag.ack_ptr);} */

#define SERVER_SEND_DATA(_SS_proc,_SS_src,_SS_dst,_SS_size) {armci_send_data_to_client(_SS_proc,_SS_src,_SS_size,_SS_dst);}
#define SERVER_GET_DATA(_SG_proc,_SG_src,_SG_dst,_SG_size) {armci_get_data_from_client(_SG_proc,_SG_src,_SG_size,_SG_dst);}


/*\ descriptors will have unique ID's for the wait on descriptor routine to
 * complete a descriptor and know where it came from
\*/

#define NUMOFBUFFERS (MAX_BUFS+MAX_SMALL_BUFS)
#define DSCRID_FROMBUFS 1
#define DSCRID_FROMBUFS_END (DSCRID_FROMBUFS+NUMOFBUFFERS)

#define DSCRID_NBDSCR 10000
#define DSCRID_NBDSCR_END (10000+MAX_PENDING)

#define DSCRID_SCATGAT 20000
#define DSCRID_SCATGAT_END 20000+MAX_PENDING

#define DSCRID_RMW 30000
#define DSCRID_RMW_END 30000+9999

#if defined(PEND_BUFS)
#define DSCRID_PENDBUF (40000)
#define DSCRID_PENDBUF_END (DSCRID_PENDBUF + 2*PENDING_BUF_NUM+1)

#define DSCRID_IMMBUF_RECV     (200000)
#define DSCRID_IMMBUF_RECV_END (600000)

#define DSCRID_IMMBUF_RESP     (600000)
#define DSCRID_IMMBUF_RESP_END (1000000)
#endif

extern double MPI_Wtime();
static double inittime0=0,inittime1=0,inittime2=0,inittime3=0,inittime4=0;

static int mark_buf_send_complete[NUMOFBUFFERS+1];
static sr_descr_t armci_vapi_client_nbsdscr_array[MAX_PENDING];
static sr_descr_t armci_vapi_client_nbrdscr_array[MAX_PENDING];
static sr_descr_t armci_vapi_serv_nbsdscr_array[MAX_PENDING];
static sr_descr_t armci_vapi_serv_nbrdscr_array[MAX_PENDING];

void armci_server_transport_cleanup();
/********************FUNCTIONS TO CHECK OPENIB RETURN STATUS*******************/
void armci_check_status(int debug, int rc,char *msg)
{
  dassertp(debug,rc==0,("%d: %s, rc=%d\n",armci_me,msg,rc));
/*     if(debug)printf("%d:%s, rc = %d\n", armci_me,msg, rc); */
/*     if(rc!=0)armci_die(msg,rc); */
}

void armci_vapi_check_return(int debug, int ret, const char *ss)
{
}

void armci_vapi_print_dscr_info(struct ibv_send_wr *sr, struct ibv_recv_wr *rr)
{
int i;
    if(rr){
       printf("\n%d:print_dscr rr id=%ld sg_lst_len=%d",
              armci_me, rr->wr_id, rr->num_sge);
       for (i = 0; i < rr->num_sge; i++) {
         printf("\n\t:sg_entry=%d addr=%p len=%d",
                i, rr->sg_list[i].addr, rr->sg_list[i].length);
       }
       fflush(stdout);
    }
    if(sr){
       printf("\n%d:print_dscr sr id=%d opcode=%d sg_lst_len=%d",
              armci_me, sr->wr_id, sr->opcode, sr->num_sge);
       for (i = 0; i < sr->num_sge; i++) {
         printf("\n\t:sg_entry=%d addr=%p len=%d",
                i, sr->sg_list[i].addr, sr->sg_list[i].length);
       }
       fflush(stdout);
    }
}

/*****************END FUNCTIONS TO CHECK VAPI RETURN STATUS********************/

void armci_recv_complete(struct ibv_recv_wr *rcv_dscr, 
        char *from, int numofrecvs) 
{
int rc=0;
struct ibv_wc pdscr1;
struct ibv_wc *pdscr = &pdscr1;
sr_descr_t *rdscr_arr;
vapi_nic_t *nic;
int debug,i,done=0;

    if(SERVER_CONTEXT){
       rdscr_arr = armci_vapi_serv_nbrdscr_array;
       nic=CLN_nic;
       debug = DEBUG_SERVER;
    }
    else{
       rdscr_arr = armci_vapi_client_nbrdscr_array;
       nic=SRV_nic;
       debug = DEBUG_CLN;
    }
    if(debug){
       printf("\n%d%s:recv_complete called from %s id=%ld\n",armci_me,
               ((SERVER_CONTEXT)?"(s)":" "),from,rcv_dscr->wr_id);fflush(stdout);
    }
    for(i=0;i<numofrecvs;i++){
    do{
      while(rc == 0) {
         rc = ibv_poll_cq(nic->rcq, 1, pdscr);
      }
      dassertp(DBG_POLL|DBG_ALL,rc>=0,
	       ("%d: rc=%d id=%d status=%d (%d/%d)\n",
		armci_me,rc,pdscr->wr_id,pdscr->status,i,numofrecvs));
      dassert1(1,pdscr->status==IBV_WC_SUCCESS,pdscr->status);
       if(debug){
         if(pdscr->wr_id >= DSCRID_SCATGAT && pdscr->wr_id < DSCRID_SCATGAT_END)
           printf("\n%d:recv from %s complete id=%d num=%d",armci_me,
             from,pdscr->wr_id,rdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].numofrecvs);
       }
       if(pdscr->wr_id >= DSCRID_SCATGAT && pdscr->wr_id < DSCRID_SCATGAT_END){
         rdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].numofrecvs--;
         if(rdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].numofrecvs==0)
           rdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].tag=0;
       }
       else if(pdscr->wr_id == (DSCRID_SCATGAT + MAX_PENDING)){
               /*this was from a blocking call, do nothing*/
         continue;
       }
       else {
         armci_die("\nclient should be posting only one kind of recv",armci_me);
       }
       rc = 0;
   }while(pdscr->wr_id!=rcv_dscr->wr_id);
   rc = 0;
   }

}

void armci_vapi_set_mark_buf_send_complete(int id)
{
    mark_buf_send_complete[id]=0;
}

void armci_send_complete(struct ibv_send_wr *snd_dscr, char *from,int numoftimes)
{
    int rc=0;
    struct ibv_wc pdscr1;
    struct ibv_wc *pdscr = &pdscr1;
    sr_descr_t *sdscr_arr;
    vapi_nic_t *nic;
    int debug,i;

    pdscr1.status = IBV_WC_SUCCESS;
    /*  bzero(&pdscr1, sizeof(pdscr1)); */
    /* printf("%d: Waiting for send with wr_id=%d to complete\n", armci_me, snd_dscr->wr_id); */
    /* fflush(stdout); */

    if(SERVER_CONTEXT){
        sdscr_arr = armci_vapi_serv_nbsdscr_array;
        nic=CLN_nic;
        debug = DEBUG_SERVER;
    }
    else{
        sdscr_arr = armci_vapi_client_nbsdscr_array;
        nic=SRV_nic;
        debug = DEBUG_CLN;
    }

    if(debug) {
        printf("\n%d%s:send_complete called from %s id=%ld nt=%d\n",armci_me,
                ((SERVER_CONTEXT)?"(s)":" "),from,snd_dscr->wr_id,numoftimes);
        fflush(stdout);
    }
    for(i=0;i<numoftimes;i++){
        do{
            while(rc == 0){  
#if defined(PEND_BUFS) 
                if(SERVER_CONTEXT)
                    rc = ibv_poll_cq(nic->rcq,1,pdscr);
                else
#endif
                    rc = ibv_poll_cq(nic->scq,1, pdscr);
            }  
            dassertp(DBG_POLL|DBG_ALL,rc>=0,
                    ("%d:rc=%d status=%d id=%d (%d/%d)",armci_me,
                     rc,pdscr->status,(int)pdscr->wr_id,i,numoftimes));
            dassert1(1,pdscr->status==IBV_WC_SUCCESS,pdscr->status);
            /*       printf("%d: Obtained completion of wr_id=%d\n", armci_me, pdscr->wr_id); */
            /*       fflush(stdout); */
            if(SERVER_CONTEXT){
                if(debug)printf("%d:completed id %d i=%d\n",armci_me,pdscr->wr_id,i);
                if(pdscr->wr_id >=DSCRID_SCATGAT && pdscr->wr_id < DSCRID_SCATGAT_END){
                    sdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].numofsends--;
                    if(sdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].numofsends==0)
                        sdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].tag=0;
                }
                else if(pdscr->wr_id >=armci_nproc && pdscr->wr_id < 2*armci_nproc){
                    /*its coming from send_data_to_client just return*/
                }
#if defined(PEND_BUFS)
                else if(pdscr->wr_id >= DSCRID_IMMBUF_RESP && pdscr->wr_id>DSCRID_IMMBUF_RESP_END) {
                    /*send from server to client completed*/
                }
#endif
                else armci_die("server send complete got weird id",pdscr->wr_id);
            }
            else{
                if(debug)printf("%d:completed id %d i=%d\n",armci_me,pdscr->wr_id,i);
                if(pdscr->wr_id >=DSCRID_FROMBUFS && pdscr->wr_id < DSCRID_FROMBUFS_END) {
                    /*	   printf("%d: marking send buffer %d as complete\n", armci_me, pdscr->wr_id);*/
                    mark_buf_send_complete[pdscr->wr_id]=1;
                }
                else if(pdscr->wr_id >=DSCRID_NBDSCR && pdscr->wr_id < DSCRID_NBDSCR_END){
                    sdscr_arr[pdscr->wr_id-DSCRID_NBDSCR].numofsends--;
                    if(sdscr_arr[pdscr->wr_id-DSCRID_NBDSCR].numofsends==0)
                        sdscr_arr[pdscr->wr_id-DSCRID_NBDSCR].tag=0;
                }
                else if(pdscr->wr_id >=DSCRID_SCATGAT && pdscr->wr_id < DSCRID_SCATGAT_END){
                    sdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].numofsends--;
                    if(sdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].numofsends==0)
                        sdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].tag=0;
                }
                else if(pdscr->wr_id == (DSCRID_SCATGAT + MAX_PENDING)){
                    /* 	   printf("%d: completed a blocking scatgat descriptor\n", armci_me); */
                    /*this was from a blocking call, do nothing*/
                    continue;
                }
                else armci_die("client send complete got weird id",pdscr->wr_id);
            }
            rc = 0;
        }while(pdscr->wr_id!=snd_dscr->wr_id);
        rc = 0;
    }
}


void armci_dscrlist_recv_complete(int tag, char* from,sr_descr_t *dscr)
{
    int i,nr,j;
    sr_descr_t *retdscr,*rdscr_arr;
    
    if(dscr == NULL){
        if(SERVER_CONTEXT)
            rdscr_arr = armci_vapi_serv_nbrdscr_array;
        else
            rdscr_arr = armci_vapi_client_nbrdscr_array;

        for(i=0;i<MAX_PENDING;i++){
            if(rdscr_arr[i].tag==tag)
                break;
        }

        if(i==MAX_PENDING)return;
        retdscr = &rdscr_arr[i];
    }
    else
        retdscr=dscr;

    nr = retdscr->numofrecvs;
    armci_recv_complete(&(retdscr->rdescr),"(s)list_send_complete",nr);
}


void armci_dscrlist_send_complete(int tag,char *from, sr_descr_t *dscr)
{
    int i,ns,j;
    sr_descr_t *retdscr,*sdscr_arr;
    if(dscr==NULL){
        if(SERVER_CONTEXT)
            sdscr_arr = armci_vapi_serv_nbsdscr_array;
        else
            sdscr_arr = armci_vapi_client_nbsdscr_array;

        for(i=0;i<MAX_PENDING;i++){
            if(sdscr_arr[i].tag==tag)
                break;
        }
        if(i==MAX_PENDING)return;
        retdscr=&sdscr_arr[i];
    }
    else
        retdscr=dscr;

    ns = retdscr->numofsends;

    armci_send_complete(&(retdscr->sdescr),"dscrlist_send_complete",ns);

}

void armci_client_nbcall_complete(sr_descr_t *dscr, int tag, int op)
{
    if(tag != dscr->tag)
        return;

	THREAD_LOCK(armci_user_threads.net_lock);

    if(op == GET){
       if(dscr->issg){
         if(dscr->numofrecvs>0)
           armci_dscrlist_recv_complete(tag,"armci_client_nbcall_complete recv",
                           dscr);
       }
       else{
         if(dscr->numofsends>0)
           armci_dscrlist_send_complete(tag,"armci_client_nbcall_complete send",
                           dscr);
       }
    }
    if(op == PUT){
       if(dscr->numofsends>0)
         armci_dscrlist_send_complete(tag,"armci_client_nbcall_complete send",
                         dscr);
    }

	THREAD_UNLOCK(armci_user_threads.net_lock);
}


static int cur_serv_pend_descr;
static int cur_client_pend_descr;

sr_descr_t *armci_vapi_get_next_rdescr(int nbtag,int sg)
{
static int serverthreadavail=-1; /*client thread can't touch this*/
static int clientthreadavail=-1; /*server thread can't touch this*/
int avail,newavail;
sr_descr_t *retdscr,*rdscr_arr;

    if(SERVER_CONTEXT){
       rdscr_arr = armci_vapi_serv_nbrdscr_array;
       avail = serverthreadavail;
       /*printf("\n%d:serv thread avail=%d",armci_me,serverthreadavail);*/
    }
    else{
       rdscr_arr = armci_vapi_client_nbrdscr_array;
       avail = clientthreadavail;
    }
    if(avail==-1){
       int i;
       for(i=0;i<MAX_PENDING;i++){
         rdscr_arr[i].tag=0;
         bzero(&rdscr_arr[i].rdescr,sizeof(struct ibv_recv_wr)); 
         if(sg)
           rdscr_arr[i].rdescr.wr_id = DSCRID_SCATGAT + i;
         else
           rdscr_arr[i].rdescr.wr_id = DSCRID_NBDSCR + i; 
       }
       avail=0;
    }

    if(rdscr_arr[avail].tag!=0){
       armci_dscrlist_recv_complete(rdscr_arr[avail].tag,
                         "armci_vapi_get_next_rdescr",&rdscr_arr[avail]);
    }

    rdscr_arr[avail].tag=nbtag;
    rdscr_arr[avail].issg=sg;
    retdscr= (rdscr_arr+avail);

    memset(&retdscr->rdescr,0,sizeof(struct ibv_recv_wr));

    if(sg)
       retdscr->rdescr.wr_id = DSCRID_SCATGAT + avail;
    else{
       retdscr->rdescr.wr_id = DSCRID_NBDSCR + avail; 
       retdscr->numofrecvs=1;
    }

    newavail = (avail+1)%MAX_PENDING;

    if(SERVER_CONTEXT){
      cur_serv_pend_descr = avail;
      serverthreadavail=newavail;
    }
    else{
      cur_client_pend_descr = avail;
      clientthreadavail=newavail;
    }

    return(retdscr);

}

sr_descr_t *armci_vapi_get_next_sdescr(int nbtag,int sg)
{
    static int serverthreadavail=-1; /*client thread can't touch this*/
    static int clientthreadavail=-1; /*server thread can't touch this*/
    int avail,newavail;
    sr_descr_t *retdscr,*sdscr_arr;

    if(SERVER_CONTEXT){
        sdscr_arr = armci_vapi_serv_nbsdscr_array;
        avail = serverthreadavail;
    }
    else{
        sdscr_arr = armci_vapi_client_nbsdscr_array;
        avail = clientthreadavail;
    }

    if(avail==-1){ /*first call*/
        int i;
        for(i=0;i<MAX_PENDING;i++){
            sdscr_arr[i].tag=0;
            bzero(&sdscr_arr[i].sdescr,sizeof(struct ibv_send_wr));
            if(sg)
                sdscr_arr[i].sdescr.wr_id = DSCRID_SCATGAT+i;
            else
                sdscr_arr[i].sdescr.wr_id = DSCRID_NBDSCR + i;
        }
        avail=0;
    }

    if(sdscr_arr[avail].tag!=0){
        armci_dscrlist_send_complete(sdscr_arr[avail].tag,
                "armci_vapi_get_next_sdescr",&sdscr_arr[avail]);
    }

    sdscr_arr[avail].tag=nbtag;
    sdscr_arr[avail].issg=sg;
    retdscr= (sdscr_arr+avail);

    memset(&retdscr->sdescr,0,sizeof(struct ibv_recv_wr));

    if(sg)
        retdscr->sdescr.wr_id = DSCRID_SCATGAT + avail;
    else{
        retdscr->sdescr.wr_id = DSCRID_NBDSCR + avail;
        retdscr->numofsends=1;
    }

    newavail = (avail+1)%MAX_PENDING;

    if(SERVER_CONTEXT){
        cur_serv_pend_descr = avail;
        serverthreadavail=newavail;
    }
    else{
        cur_client_pend_descr = avail;
        clientthreadavail=newavail;
    }
    return(retdscr);
}

void armci_wait_for_server()
{
    armci_server_terminating = 1;
}


/* ibv_create_qp does not use separate structure to return properties,
   seems it is all inside ibv_qp */
static void armci_create_qp(vapi_nic_t *nic, struct ibv_qp **qp)
{
    struct ibv_qp_init_attr initattr;

    bzero(&initattr, sizeof(struct ibv_qp_init_attr));

    *qp = NULL;

    initattr.cap.max_send_wr = armci_max_qp_ous_swr;
    initattr.cap.max_recv_wr = armci_max_qp_ous_rwr;
    initattr.cap.max_recv_sge = armci_max_num_sg_ent;
    initattr.cap.max_send_sge = armci_max_num_sg_ent;
#if defined(PEND_BUFS)
    if(nic==CLN_nic) {
        initattr.send_cq = nic->rcq;
        initattr.recv_cq = nic->rcq;
    }
    else 
#endif
    {
        initattr.send_cq = nic->scq;
        initattr.recv_cq = nic->rcq;      
    }
    initattr.qp_type = IBV_QPT_RC;

    *qp = ibv_create_qp(nic->ptag, &initattr);
    dassert(1,*qp!=NULL);

    /* The value of inline size should be dependent on number of processes
     * */
    if (armci_nproc >= MAX_PROC_INLINE_SIZE) {
        armci_vapi_max_inline_size = -1;
    }
    else {
        armci_vapi_max_inline_size = initattr.cap.max_inline_data;
    }
}

int armci_openib_sl;
int armci_openib_server_poll;
void armci_openib_env_init()
{
    char *value;

    if ((value = getenv("ARMCI_OPENIB_USE_SL")) != NULL){
        armci_openib_sl = atoi(value);
    } 
    else {
        armci_openib_sl = 0;
    }
        
    /* Don't enable server polling by default */
    if ((value = getenv("ARMCI_OPENIB_SERVER_POLL")) != NULL){
        armci_openib_server_poll = atoi(value);
    } 
    else {
        armci_openib_server_poll = 0;
    }
}

static void armci_init_nic(vapi_nic_t *nic, int scq_entries, int
        rcq_entries)
{
    int rc, ndevs, i;
    struct ibv_device **devs=NULL;
    struct ibv_context *cxt;

    if (nic == SRV_nic) {
        /* Initialize OpenIB runtime variables only once*/
        armci_openib_env_init();
    }

    bzero(nic,sizeof(vapi_nic_t));
    nic->lid_arr = (uint16_t *)calloc(armci_nproc,sizeof(uint16_t));
    dassert(1,nic->lid_arr!=NULL);

    devs = ibv_get_device_list(&ndevs);

    nic->handle = ibv_open_device(*devs); 

    nic->maxtransfersize = MAX_RDMA_SIZE;

    nic->vendor = ibv_get_device_name(*devs);

    rc = ibv_query_device(nic->handle, &nic->attr);

    int down_port_count_check = 0;
    for (i = 1; i <= 2; i++) {
        rc = ibv_query_port(nic->handle, (uint8_t)i, &nic->hca_port);
        if (IBV_PORT_ACTIVE == nic->hca_port.state) {
            nic->active_port = i;
            break;
        } 
        else {
            down_port_count_check++;
        }
    }

    /* Assert that the number of inactive ports is not equal to the number
     * of down ports on any adapter */
    assert(down_port_count_check != 2);
        
    /*save the lid for doing a global exchange later */
    nic->lid_arr[armci_me] = nic->hca_port.lid;

    /*allocate tag (protection domain) */
    nic->ptag = ibv_alloc_pd(nic->handle);

    /* properties of scq and rcq required for the cq number, this also needs
     * to be globally exchanged
     */
    nic->scv = 1;
    nic->rcv = 2;
    nic->scq = nic->rcq = NULL; 
    
    if(scq_entries) {
        nic->sch = ibv_create_comp_channel(nic->handle);
        nic->scq = ibv_create_cq(nic->handle, 16000,
                nic->scq_cntx,nic->sch, 0);
    }
    
    if(rcq_entries) {
        nic->rch = ibv_create_comp_channel(nic->handle);
        nic->rcq = ibv_create_cq(nic->handle, 32768,
                nic->rcq_cntx,nic->rch, 0);
    }
    
    ibv_free_device_list(devs);

    armci_max_num_sg_ent = 29; 
    armci_max_qp_ous_swr = 100;
    armci_max_qp_ous_rwr = 50;

    char *value;
    if ((value = getenv("ARMCI_USE_ODCM")) != NULL){
        armci_use_odcm = atoi(value);
    } else {
        armci_use_odcm = 0;
    }

    armci_use_lazy_break = 0;

    if(armci_max_qp_ous_rwr + armci_max_qp_ous_swr>nic->attr.max_qp_wr){
        armci_max_qp_ous_swr = nic->attr.max_qp_wr/16;
        armci_max_qp_ous_rwr = nic->attr.max_qp_wr - armci_max_qp_ous_swr;
    }
    if(armci_max_num_sg_ent >= nic->attr.max_sge){
        armci_max_num_sg_ent = nic->attr.max_sge - 1;
    }
    
}


void armci_setaffinity(char *cpu_mapping) {
    long N_CPUs_online;
    cpu_set_t affinity_mask;
    unsigned long affinity_mask_len = sizeof(affinity_mask);
    char *tp;
    char *cp;
    char tp_str[8];
    int cpu, i, j;
               
    
    if (!armci_use_affinity)
        return;

    /*Get number of CPU on machine */
    if ((N_CPUs_online = sysconf(_SC_NPROCESSORS_ONLN)) < 1) {
        perror("sysconf");
    }
    
    if (cpu_mapping) {
        tp = cpu_mapping;
        j = 0;
        while (*tp != '\0') {
            i = 0;
            cp = tp;
            while (*cp != '\0' && *cp != ',' && *cp != ':') {
                cp++;
                i++;
            }           
            strncpy(tp_str, tp, i);
            tp_str[i] = '\0';
            cpu = atoi(tp_str);

            if (j == armci_me - armci_master) {
                CPU_ZERO(&affinity_mask);
                CPU_SET(cpu, &affinity_mask);
                if (sched_setaffinity(0,
                    affinity_mask_len, &affinity_mask)<0 ) {
                    perror("sched_setaffinity");
                }
                break;
            }
            tp = cp + 1;
            j++;
        }

        free(cpu_mapping);
    } else {
        CPU_ZERO(&affinity_mask);
        CPU_SET(((armci_me) - armci_master) %N_CPUs_online, &affinity_mask);

        if (sched_setaffinity(0,affinity_mask_len,&affinity_mask)<0 ) {
            perror("sched_setaffinity");
        }
    }
}

/****************MEMORY ALLOCATION REGISTRATION DEREGISTRATION****************/
static char * serv_malloc_buf_base;
#if ARMCI_ENABLE_GPC_CALLS
extern gpc_buf_t *gpc_req;
#endif
void armci_server_alloc_bufs()
{
    int rc;
    int mod, bytes, total, extra =sizeof(struct ibv_recv_wr)*MAX_DESCR+SIXTYFOUR;
    int mhsize = armci_nproc*sizeof(armci_vapi_memhndl_t); /* ack */
    char *tmp, *tmp0;
    int i, j=0;
#if defined(PEND_BUFS)
    int clients = (IMM_BUF_NUM+1)*armci_nproc;
#else
    int clients = armci_nproc;
#endif

    /* allocate memory for the recv buffers-must be alligned on 64byte bnd */
    /* note we add extra one to repost it for the client we are received req */
    bytes = (clients+1)*sizeof(vapibuf_t)+sizeof(vapibuf_ext_t) + extra+ mhsize
#if ARMCI_ENABLE_GPC_CALLS
      + MAX_GPC_REQ * sizeof(gpc_buf_t)
#endif
#if defined(PEND_BUFS)
      + (clients+1)*IMM_BUF_LEN
      + PENDING_BUF_NUM*(sizeof(vapibuf_pend_t)+PENDING_BUF_LEN)
#endif
      + sizeof(long)
      + 7*SIXTYFOUR;
    total = bytes + SIXTYFOUR;
    if(total%4096!=0)
       total = total - (total%4096) + 4096;
    tmp0=tmp = malloc(total);
    serv_malloc_buf_base = tmp0;

    dassert1(1,tmp!=NULL,(int)total);
    /* stamp the last byte */
    serv_tail= tmp + bytes+SIXTYFOUR-1;
    *serv_tail=SERV_STAMP;
    /* allocate memory for client memory handle to support put response
     *         in dynamic memory registration protocols */
    CLN_handle = (armci_vapi_memhndl_t*)tmp;
    memset(CLN_handle,0,mhsize); /* set it to zero */
    tmp += mhsize;

#if ARMCI_ENABLE_GPC_CALLS
    /* gpc_req memory*/
    tmp += SIXTYFOUR - ((ssize_t)tmp % SIXTYFOUR);
    gpc_req = (gpc_buf_t *)tmp;
    tmp += MAX_GPC_REQ * sizeof(gpc_buf_t);
#endif

    /* setup descriptor memory */
    tmp += SIXTYFOUR - ((ssize_t)tmp % SIXTYFOUR);
    serv_descr_pool.descr= (struct ibv_recv_wr *)(tmp);
    tmp += extra;

    /* setup ack buffer*/
    tmp += SIXTYFOUR - ((ssize_t)tmp % SIXTYFOUR);
    ack_buf = (long *)(tmp);
    *ack_buf=ARMCI_STAMP;
    tmp += sizeof(long);

    /* setup buffer pointers */
    tmp += SIXTYFOUR - ((ssize_t)tmp % SIXTYFOUR);
    serv_buf_arr = (vapibuf_t **)malloc(sizeof(vapibuf_t*)*clients);
    for(i=0;i<clients;i++){
      serv_buf_arr[i] = (vapibuf_t*)(tmp) + i;
    }
    tmp = (char *)(serv_buf_arr[0]+clients);

#if defined(PEND_BUFS)
    /*setup buffers in immediate buffers*/
    tmp += SIXTYFOUR - ((ssize_t)tmp % SIXTYFOUR);
    for(i=0; i<clients; i++) {
      serv_buf_arr[i]->buf = tmp + i*IMM_BUF_LEN;
    }
    tmp += clients*IMM_BUF_LEN;

    /*setup pending buffers*/
    tmp += SIXTYFOUR - ((ssize_t)tmp % SIXTYFOUR);
    serv_pendbuf_arr = (vapibuf_pend_t *)(tmp);
    tmp=(char *)(serv_pendbuf_arr+PENDING_BUF_NUM);
    tmp += SIXTYFOUR - ((ssize_t)tmp % SIXTYFOUR);
    for(i=0; i<PENDING_BUF_NUM; i++) {
      serv_pendbuf_arr[i].buf = tmp+i*PENDING_BUF_LEN;
      assert(serv_pendbuf_arr[i].buf != NULL);
    }
    tmp += PENDING_BUF_NUM*PENDING_BUF_LEN;
    MessageRcvBuffer = NULL;    
#else
    tmp += SIXTYFOUR - ((ssize_t)tmp % SIXTYFOUR);
    spare_serv_buf = (vapibuf_t *)tmp; /* spare buffer is at the end */
    spare_serv_bufptr = spare_serv_buf;    /* save the pointer for later */
    serv_buf =(vapibuf_ext_t*)(spare_serv_buf+1);
    tmp = (char *)(serv_buf+1);

    MessageRcvBuffer = serv_buf->buf;
#endif

   flag_arr = (int *)malloc(sizeof(int)*armci_nproc);
   for (i =0; i<armci_nproc; i++) flag_arr[i] = 9999;

    if(DEBUG_SERVER){
      printf("\n%d(s):registering mem %p %dbytes ptag=%ld handle=%d\n",
             armci_me, tmp0,total,CLN_nic->ptag,CLN_nic->handle);fflush(stdout);
    }

    serv_memhandle.memhndl = ibv_reg_mr(CLN_nic->ptag, tmp0, total,
                                        IBV_ACCESS_LOCAL_WRITE |
                                        IBV_ACCESS_REMOTE_WRITE |
                                        IBV_ACCESS_REMOTE_READ);
    dassert1(1,serv_memhandle.memhndl!=NULL,total);
    serv_memhandle.lkey=serv_memhandle.memhndl->lkey;
    serv_memhandle.rkey=serv_memhandle.memhndl->rkey;

    /* exchange address of ack/memhandle flag on servers */
    if(DEBUG_SERVER){
       printf("%d(s):registered mem %p %dbytes mhandle=%d mharr starts%p\n",
              armci_me, tmp0, total, serv_memhandle.memhndl,CLN_handle);
       fflush(stdout);
    }
}

static char * client_malloc_buf_base;
char * armci_vapi_client_mem_alloc(int size)
{
    int rc;
    int mod, total;
    int extra = MAX_DESCR*sizeof(struct ibv_recv_wr)+SIXTYFOUR;
    char *tmp,*tmp0;

    /*we use the size passed by the armci_init_bufs routine instead of bytes*/

    total = size + extra + 2*SIXTYFOUR;

    if(total%4096!=0)
       total = total - (total%4096) + 4096;
    tmp0  = tmp = malloc(total);
    dassert1(1,tmp!=NULL,total);
    client_malloc_buf_base = tmp;
#if 0
    /*SK: could this lead to a problem at ibv_reg_mr() because of unfixed 'total'?*/
    if(ALIGN64ADD(tmp0))tmp0+=ALIGN64ADD(tmp0);
#endif
    /* stamp the last byte */
    client_tail= tmp + extra+ size +2*SIXTYFOUR-1;
    *client_tail=CLIENT_STAMP;

    /* we also have a place to store memhandle for zero-copy get */
    pinned_handle =(armci_vapi_memhndl_t *) (tmp + extra+ size +SIXTYFOUR-16);

    mod = ((ssize_t)tmp)%SIXTYFOUR;
    client_descr_pool.descr= (struct ibv_recv_wr*)(tmp+SIXTYFOUR-mod);
    tmp += extra;

    client_memhandle.memhndl = ibv_reg_mr(SRV_nic->ptag, tmp0, total,
                                          IBV_ACCESS_LOCAL_WRITE |
                                          IBV_ACCESS_REMOTE_WRITE |
                                          IBV_ACCESS_REMOTE_READ);
    dassert(1,client_memhandle.memhndl!=NULL);
    
    client_memhandle.lkey = client_memhandle.memhndl->lkey;
    client_memhandle.rkey = client_memhandle.memhndl->rkey;
    handle_array[armci_me].lkey = client_memhandle.lkey;
    handle_array[armci_me].rkey = client_memhandle.rkey;
  
    handle_array[armci_me].memhndl = client_memhandle.memhndl;

    if(DEBUG_INIT){
       printf("%d: registered client memory %p %dsize tmp=%p \n",
               armci_me,tmp0, total, tmp);
       fflush(stdout);
    }
    /*now that we have the handle array, we get every body elses RDMA handle*/
    total = (sizeof(armci_vapi_memhndl_t)*armci_nproc)/sizeof(int);
    armci_msg_gop_scope(SCOPE_ALL,handle_array,total,"+",ARMCI_INT);

    return(tmp);
}


void armci_server_register_region(void *ptr,long bytes, ARMCI_MEMHDL_T *memhdl)
{
    bzero(memhdl,sizeof(ARMCI_MEMHDL_T));

    memhdl->memhndl = ibv_reg_mr(CLN_nic->ptag, ptr, bytes,
               IBV_ACCESS_LOCAL_WRITE |
               IBV_ACCESS_REMOTE_WRITE |
               IBV_ACCESS_REMOTE_READ);
    dassert(1,memhdl->memhndl!=NULL);

    memhdl->lkey=memhdl->memhndl->lkey;
    memhdl->rkey=memhdl->memhndl->rkey;

    if(DEBUG_SERVER){
       printf("\n%d(s):registered lkey=%d rkey=%d ptr=%p end=%p %p\n",armci_me,
               memhdl->lkey,memhdl->rkey,ptr,(char *)ptr+bytes,memhdl);
       fflush(stdout);
    }
}

int armci_pin_contig_hndl(void *ptr, size_t bytes, ARMCI_MEMHDL_T *memhdl)
{
    memhdl->memhndl = ibv_reg_mr(SRV_nic->ptag, ptr, bytes,
               IBV_ACCESS_LOCAL_WRITE |
               IBV_ACCESS_REMOTE_WRITE |
               IBV_ACCESS_REMOTE_READ);
    dassert(1,memhdl->memhndl!=NULL);
    memhdl->lkey=memhdl->memhndl->lkey;
    memhdl->rkey=memhdl->memhndl->rkey;
    if(0){
       printf("\n%d:registered lkey=%d rkey=%d ptr=%p end=%p\n",armci_me,
               memhdl->lkey,memhdl->rkey,ptr,(char *)ptr+bytes);fflush(stdout);
    }
    return 1;
}

#if 1
void armci_network_client_deregister_memory(ARMCI_MEMHDL_T *mh)
{
    int rc;
    rc = ibv_dereg_mr(mh->memhndl);
    dassert1(1,rc==0,rc);
    armci_vapi_check_return(DEBUG_FINALIZE,rc,
                        "armci_network_client_deregister_memory:deregister_mr");
}
void armci_network_server_deregister_memory(ARMCI_MEMHDL_T *mh)
{
    int rc;
return; /* ??? why ??? */
    printf("\n%d:deregister ptr=%p",armci_me,mh);fflush(stdout);
    rc = ibv_dereg_mr(mh->memhndl);
    dassert1(1,rc==0,rc);
    armci_vapi_check_return(DEBUG_FINALIZE,rc,
                        "armci_network_server_deregister_memory:deregister_mr");
}
#else
#   define armci_network_client_deregister_memory(mh)           \
           armci_vapi_check_return(DEBUG_FINALIZE,              \
                                   ibv_dereg_mr(mh->memhndl),   \
                                   "armci_network_client_deregister_memory:deregister_mr")
#   define armci_network_server_deregister_memory(mh)           \
           armci_vapi_check_return(DEBUG_FINALIZE,              \
                                   ibv_dereg_mr(mh->memhndl),   \
                                   "armci_network_server_deregister_memory:deregister_mr")
#endif

void armci_set_serv_mh()
{
int s, ratio = sizeof(ack_t)/sizeof(int);
    /* first collect addrresses on all masters */
    if(armci_me == armci_master){
       SRV_ack[armci_clus_me].prem_handle=CLN_handle;
       SRV_ack[armci_clus_me].handle =serv_memhandle;
       armci_msg_gop_scope(SCOPE_MASTERS,SRV_ack,ratio*armci_nclus,"+",
                           ARMCI_INT);
    }
    /* next master broadcasts the addresses within its node */
    armci_msg_bcast_scope(SCOPE_NODE,SRV_ack,armci_nclus*sizeof(ack_t),
                          armci_master);

    /* Finally save address corresponding to my id on each server */
    for(s=0; s< armci_nclus; s++){
       SRV_ack[s].prem_handle += armci_me;
    }

}
/**********END MEMORY ALLOCATION REGISTRATION AND DEREGISTRATION**************/

/*\
 * init_connections, client_connect_to_servers -- client code
 * server_initial_connection, all_data_server -- server code 
\*/ 
void armci_init_connections()
{
    int c,s;
    int sz;
    uint32_t *tmpbuf;
    int *tmparr;
    if(TIME_INIT)inittime0 = MPI_Wtime(); 

#if defined(PEND_BUFS)
    armci_pbuf_init_buffer_env();
#endif

    armci_setaffinity(NULL);

    /* initialize nic connection for qp numbers and lid's */
    armci_init_nic(SRV_nic,1,1);
    for(c=0; c < NUMOFBUFFERS+1; c++) {
        mark_buf_send_complete[c]=1;
    }
    _gtmparr = (int *)calloc(armci_nproc,sizeof(int)); 

    /*qp_numbers and lids need to be exchanged globally*/
    tmparr = (int *)calloc(armci_nproc,sizeof(int));
    tmparr[armci_me] = SRV_nic->lid_arr[armci_me];
    sz = armci_nproc;
    armci_msg_gop_scope(SCOPE_ALL,tmparr,sz,"+",ARMCI_INT);
    for(c=0;c<armci_nproc;c++){
        SRV_nic->lid_arr[c]=tmparr[c];
        tmparr[c]=0;
    }
    /*SRV_con is for client to connect to servers */
    SRV_con=(armci_connect_t *)malloc(sizeof(armci_connect_t)*armci_nclus);
    dassert1(1,SRV_con!=NULL,sizeof(armci_connect_t)*armci_nclus);
    bzero(SRV_con,sizeof(armci_connect_t)*armci_nclus);

    CLN_con=(armci_connect_t*)malloc(sizeof(armci_connect_t)*armci_nproc);
    dassert1(1,CLN_con!=NULL,sizeof(armci_connect_t)*armci_nproc);
    bzero(CLN_con,sizeof(armci_connect_t)*armci_nproc);

    /*every client creates a qp with every server other than the one on itself*/
    SRV_rqpnums = (uint32_t*)malloc(sizeof(uint32_t)*armci_nproc);
    dassert(1,SRV_rqpnums);
    tmpbuf = (uint32_t*)calloc(armci_nproc,sizeof(uint32_t));
    dassert(1,tmpbuf);

    sz = armci_nproc*(sizeof(uint32_t)/sizeof(int));
    armci_vapi_max_inline_size = 0;


    if (!armci_use_odcm) {
        for(s = 0; s < armci_nclus; s++){
            armci_connect_t *con = SRV_con + s;
            {
                armci_create_qp(SRV_nic,&con->qp);
                con->sqpnum  = con->qp->qp_num;
                tmpbuf[armci_clus_info[s].master] = con->qp->qp_num;
                con->lid = SRV_nic->lid_arr[s];
            }
        }
        MPI_Alltoall(tmpbuf,sizeof(uint32_t),MPI_CHAR,SRV_rqpnums,
                sizeof(uint32_t),MPI_CHAR,ARMCI_COMM_WORLD);
        free(tmpbuf);

        if(armci_me != armci_master) {
            free(SRV_rqpnums);
            SRV_rqpnums = NULL;
        }
    }
    else {
        for(s = 0; s < armci_nclus; s++){
            armci_connect_t *con = SRV_con + s;
            con->state = QP_INACTIVE;
        }
    }

    SRV_ack = (ack_t*)calloc(armci_nclus, sizeof(ack_t));
    dassert1(1,SRV_ack!=NULL,armci_nclus*sizeof(ack_t));

    handle_array = (armci_vapi_memhndl_t *)calloc(sizeof(armci_vapi_memhndl_t),
            armci_nproc);
    dassert1(1,handle_array!=NULL,sizeof(armci_vapi_memhndl_t)*armci_nproc);
   
    if (armci_use_odcm) {
        setup_ud_channel();
    }
}

static void vapi_connect_client()
{
    int i, start, sz=0, c, rc;
    struct ibv_qp_attr qp_attr;
    struct ibv_qp_cap qp_cap;
    enum ibv_qp_attr_mask qp_attr_mask;

    if (TIME_INIT) inittime0 = MPI_Wtime();
    if (armci_me == armci_master)
        armci_util_wait_int(&armci_vapi_server_stage1, 1, 10);
    if (TIME_INIT) printf("\n%d:wait for server to get to stage 1 time for "
                          "vapi_connect_client is %f",
                          armci_me, (inittime1 = MPI_Wtime()) - inittime0);
    sz = armci_nproc;
    if (armci_me == armci_master) {
       armci_msg_gop_scope(SCOPE_MASTERS, _gtmparr, sz, "+", ARMCI_INT);
       for (c=0; c<armci_nproc; c++) {
         CLN_nic->lid_arr[c] = _gtmparr[c];
         _gtmparr[c] = 0;
       }
       if (DEBUG_CLN) {
         printf("\n%d(svc): mylid = %d",armci_me,CLN_nic->lid_arr[armci_me]);
         fflush(stdout);
       }
    }

    armci_vapi_client_stage1 = 1;

    /* allocate and initialize connection structs */
    sz = armci_nproc*sizeof(uint32_t)/sizeof(int);

    if (armci_me == armci_master)
       armci_util_wait_int(&armci_vapi_server_stage2, 1, 10);
#if 0
    for (c = 0; c < armci_nproc; c++){
       armci_connect_t *con = CLN_con + c;
       if (armci_me != armci_master) {
         char *ptrr;
         int extra;
         ptrr = malloc(8 + sizeof(uint32_t) * armci_nproc);
         extra = ALIGNLONGADD(ptrr);
         ptrr = ptrr + extra;
         con->rqpnum = (uint32_t *)ptrr;
         bzero(con->rqpnum, sizeof(uint32_t) * armci_nproc);
       }
       armci_msg_gop_scope(SCOPE_ALL, con->rqpnum, sz, "+", ARMCI_INT);
    }
#else
    CLN_rqpnums = (uint32_t*)malloc(sizeof(uint32_t)*armci_nproc);

    if (!armci_use_odcm) {
        if(armci_me != armci_master) {
            /*just has junk*/
            CLN_rqpnumtmpbuf = (uint32_t*)malloc(sizeof(uint32_t)*armci_nproc);
        }
        dassert(1, CLN_rqpnumtmpbuf);
        MPI_Alltoall(CLN_rqpnumtmpbuf, sizeof(uint32_t), MPI_CHAR,
                CLN_rqpnums, sizeof(uint32_t), MPI_CHAR, ARMCI_COMM_WORLD);
        free(CLN_rqpnumtmpbuf);
        CLN_rqpnumtmpbuf=NULL;
#endif

        if (TIME_INIT) printf("\n%d:wait for server tog et to stage 2 time for "
                "vapi_connect_client is %f",
                armci_me, (inittime2 = MPI_Wtime()) - inittime1);
        /*armci_set_serv_mh();*/

        if (DEBUG_CLN) {
            printf("%d:all connections ready\n", armci_me);
            fflush(stdout);
        }

        /* For sanity */
        memset(&qp_attr, 0, sizeof qp_attr);
        /* Modifying  QP to INIT */
        qp_attr_mask = IBV_QP_STATE
            | IBV_QP_PKEY_INDEX
            | IBV_QP_PORT
            | IBV_QP_ACCESS_FLAGS;

        qp_attr.qp_state = IBV_QPS_INIT;
        qp_attr.pkey_index = DEFAULT_PKEY_IX;
        qp_attr.port_num = SRV_nic->active_port;
        qp_attr.qp_access_flags = IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_REMOTE_READ;

        /* start from from server on my_node -1 */
        start = (armci_clus_me == 0) ? armci_nclus - 1 : armci_clus_me - 1;
        for (i = 0; i < armci_nclus; i++) {
            armci_connect_t *con;
            con = SRV_con + i;
            rc = ibv_modify_qp(con->qp, &qp_attr, qp_attr_mask);
            dassertp(1,!rc,("%d: client RST->INIT i=%d rc=%d\n",armci_me,i,rc));
        }

        if (TIME_INIT) printf("\n%d:to init time for vapi_connect_client is %f",
                armci_me, (inittime1 = MPI_Wtime()) - inittime2);
        qp_attr_mask = IBV_QP_STATE
            | IBV_QP_MAX_DEST_RD_ATOMIC
            | IBV_QP_PATH_MTU
            | IBV_QP_RQ_PSN
            | IBV_QP_MIN_RNR_TIMER;
        memset(&qp_attr, 0, sizeof qp_attr);

        qp_attr.qp_state        = IBV_QPS_RTR;
        qp_attr.max_dest_rd_atomic   = 4;
        qp_attr.path_mtu        = IBV_MTU_1024;
        qp_attr.rq_psn          = 0;
        qp_attr.min_rnr_timer   = RNR_TIMER;

        /* AV: Adding the service level parameter */
        qp_attr.ah_attr.sl      = armci_openib_sl;

        start = (armci_clus_me == 0) ? armci_nclus - 1 : armci_clus_me - 1;
        for (i = 0; i < armci_nclus; i++) {
            armci_connect_t *con;
            armci_connect_t *conS;
            con = SRV_con + i;
#if 0
            conS = CLN_con + armci_me;
#endif
            qp_attr_mask |= IBV_QP_AV | IBV_QP_DEST_QPN;
#if 0
            qp_attr.dest_qp_num = conS->rqpnum[armci_clus_info[i].master];
#else
            qp_attr.dest_qp_num = CLN_rqpnums[armci_clus_info[i].master];
#endif
            qp_attr.ah_attr.dlid = SRV_nic->lid_arr[armci_clus_info[i].master];
            qp_attr.ah_attr.port_num = SRV_nic->active_port;

            qp_attr.ah_attr.sl = armci_openib_sl;

            rc = ibv_modify_qp(con->qp, &qp_attr, qp_attr_mask);
            dassertp(1,!rc,("%d: INIT->RTR client i=%d rc=%d\n",armci_me,i,rc));
        }
    }
    /*to to to RTS, other side must be in RTR*/

    armci_msg_barrier();
    if (TIME_INIT) printf("\n%d:init to rtr time for vapi_connect_client is %f",
                          armci_me, (inittime2 = MPI_Wtime()) - inittime1);
    armci_vapi_client_ready=1;
    
    if (!armci_use_odcm) {

        qp_attr_mask = IBV_QP_STATE
            | IBV_QP_SQ_PSN
            | IBV_QP_TIMEOUT
            | IBV_QP_RETRY_CNT
            | IBV_QP_RNR_RETRY
            | IBV_QP_MAX_QP_RD_ATOMIC;

        memset(&qp_attr, 0, sizeof qp_attr);

        qp_attr.qp_state            = IBV_QPS_RTS;
        qp_attr.sq_psn              = 0;
        qp_attr.timeout             = 18;
        qp_attr.retry_cnt           = 7;
        qp_attr.rnr_retry           = 7;
        qp_attr.max_rd_atomic  = 4;

        start = (armci_clus_me == 0) ? armci_nclus - 1 : armci_clus_me - 1;
        for (i = 0; i < armci_nclus; i++){
            armci_connect_t *con;
            con = SRV_con + i;
            rc = ibv_modify_qp(con->qp, &qp_attr, qp_attr_mask);
            dassertp(1,!rc,("%d: client RTR->RTS i=%d rc=%d\n",armci_me,i,rc));
        }
        if (TIME_INIT) printf("\n%d:rtr to rts time for vapi_connect_client is %f",
                armci_me, (inittime1 = MPI_Wtime()) - inittime2);
        free(CLN_rqpnums);
        CLN_rqpnums=NULL;
    }

}


void armci_client_connect_to_servers()
{
    extern void armci_util_wait_int(volatile int *,int,int);
    if (TIME_INIT) inittime0 = MPI_Wtime();
    _armci_buf_init();

    vapi_connect_client();
    if (armci_me == armci_master) 
       armci_util_wait_int(&armci_vapi_server_ready,1,10);
    armci_msg_barrier();
    if (DEBUG_CLN && armci_me == armci_master) {
       printf("\n%d:server_ready=%d\n",armci_me,armci_vapi_server_ready);
       fflush(stdout);
    }
    if (TIME_INIT) printf("\n%d:time for client_connect_to_s is %f",
                          armci_me,MPI_Wtime()-inittime0);
}


void armci_init_vapibuf_recv(struct ibv_recv_wr *rd, struct ibv_sge *sg_entry,
                             char *buf, int len, armci_vapi_memhndl_t *mhandle)
{
     memset(rd,0,sizeof(struct ibv_recv_wr));
     rd->next = NULL;
     rd->num_sge    = 1;
     rd->sg_list    = sg_entry;
     rd->wr_id      = 0;

     sg_entry->lkey     = mhandle->lkey;
     sg_entry->addr     = (uint64_t)buf;
     sg_entry->length   = len;
}


void armci_init_vapibuf_send(struct ibv_send_wr *sd, struct ibv_sge *sg_entry,
                             char *buf, int len, armci_vapi_memhndl_t *mhandle)
{
     sd->opcode = IBV_WR_SEND;
     sd->next = NULL;
     sd->send_flags = IBV_SEND_SIGNALED;
     sd->num_sge            = 1;
     sd->sg_list            = sg_entry;

     sg_entry->lkey     = mhandle->lkey;
     sg_entry->addr     = (uint64_t)buf;
     sg_entry->length   = len;
}


static void armci_init_cbuf_srdma(struct ibv_send_wr *sd, struct ibv_sge *sg_entry,
                                  char *lbuf, char *rbuf, int len,
                                  armci_vapi_memhndl_t *lhandle,
                                  armci_vapi_memhndl_t *rhandle)
{
     /* NOTE: sd->wr is a union, sr->wr.ud might conflict with sr->wr.rdma */
     sd->opcode = IBV_WR_RDMA_WRITE;
     sd->send_flags = IBV_SEND_SIGNALED;
     sd->next = NULL;
     sd->num_sge                    = 1;
     sd->sg_list                    = sg_entry;
     if (rhandle) sd->wr.rdma.rkey  = rhandle->rkey;
     sd->wr.rdma.remote_addr        = (uint64_t)rbuf;

     if (lhandle) sg_entry->lkey    = lhandle->lkey;
     sg_entry->addr                 = (uint64_t)lbuf;
     sg_entry->length               = len;
}


static void armci_init_cbuf_rrdma(struct ibv_send_wr *sd, struct ibv_sge
        *sg_entry, char *lbuf, char *rbuf, int len, armci_vapi_memhndl_t
        *lhandle, armci_vapi_memhndl_t *rhandle)
{
     sd->opcode = IBV_WR_RDMA_READ;
     sd->next = NULL;
     sd->send_flags = IBV_SEND_SIGNALED;
     sd->num_sge                    = 1;
     sd->sg_list                    = sg_entry;
     sd->wr.ud.remote_qkey          = 0;
     if (rhandle) sd->wr.rdma.rkey  = rhandle->rkey;
     sd->wr.rdma.remote_addr        = (uint64_t)rbuf;

     if (lhandle) sg_entry->lkey    = lhandle->lkey;
     sg_entry->addr                 = (uint64_t)lbuf;
     sg_entry->length               = len;
     /* sd->wr is a union, sr->wr.ud might conflict with sr->wr.rdma */
}


void armci_server_initial_connection()
{
  int c, rc, i, j;
    struct ibv_qp_attr qp_attr;
    struct ibv_qp_init_attr qp_init_attr;
    struct ibv_qp_cap qp_cap;
    enum ibv_qp_attr_mask qp_attr_mask;
    char *enval;
    struct ibv_recv_wr *bad_wr;

    if (TIME_INIT) 
        inittime0 = MPI_Wtime();

    if (DEBUG_SERVER) {
        printf("in server after fork %d (%d)\n",armci_me,getpid());
        fflush(stdout);
    }

#if defined(PEND_BUFS) && !defined(SERVER_THREAD)
    armci_pbuf_init_buffer_env();
#endif
    armci_init_nic(CLN_nic,1,1);
    if (!armci_openib_server_poll) {
	/*
	 * Start a notify event request immediately after creation so
	 * nothing is missed.
	 */
	rc = ibv_req_notify_cq(CLN_nic->rcq, 0);
	dassert1(1,rc==0,rc);
    }

    _gtmparr[armci_me] = CLN_nic->lid_arr[armci_me];
    armci_vapi_server_stage1 = 1;
    armci_util_wait_int(&armci_vapi_client_stage1, 1, 10);

    CLN_rqpnumtmpbuf = (uint32_t*)malloc(sizeof(uint32_t)*armci_nproc);
    dassert(1, CLN_rqpnumtmpbuf);

    if (!armci_use_odcm) {
        for (c = 0; c < armci_nproc; c++) {
            char *ptrr;           
            int extra;            
            armci_connect_t *con = CLN_con + c;
            armci_create_qp(CLN_nic, &con->qp);
            con->sqpnum = con->qp->qp_num;
            con->lid    = CLN_nic->lid_arr[c];
            CLN_rqpnumtmpbuf[c] = con->qp->qp_num;
        }
    }
    else {
        for (c = 0; c < armci_nproc; c++) {
            armci_connect_t *con = CLN_con + c;
            con->state = QP_INACTIVE;
        }            
    }

    armci_vapi_server_stage2 = 1;

    if (!armci_use_odcm) {
    qp_attr_mask = IBV_QP_STATE
                 | IBV_QP_PKEY_INDEX
                 | IBV_QP_PORT
                 | IBV_QP_ACCESS_FLAGS;

    memset(&qp_attr, 0, sizeof qp_attr);
    qp_attr.qp_state        = IBV_QPS_INIT;
    qp_attr.pkey_index      = DEFAULT_PKEY_IX;
    qp_attr.port_num        = CLN_nic->active_port;
    qp_attr.qp_access_flags = IBV_ACCESS_LOCAL_WRITE |
        IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_REMOTE_READ;

    for (c = 0; c < armci_nproc; c++) {
       armci_connect_t *con = CLN_con + c;
       rc = ibv_modify_qp(con->qp, &qp_attr, qp_attr_mask);
       dassertp(1,!rc,("%d: RTS->INIT server c=%d rc=%d\n",armci_me,c,rc));
    }
    
    memset(&qp_attr, 0, sizeof qp_attr);
    qp_attr_mask = IBV_QP_STATE
                 | IBV_QP_MAX_DEST_RD_ATOMIC
                 | IBV_QP_PATH_MTU
                 | IBV_QP_RQ_PSN
                 | IBV_QP_MIN_RNR_TIMER;
    qp_attr.qp_state           = IBV_QPS_RTR;
    qp_attr.path_mtu           = IBV_MTU_1024;          
    qp_attr.max_dest_rd_atomic = 4;
    qp_attr.min_rnr_timer      = RNR_TIMER;
    qp_attr.rq_psn             = 0;

    for(c = 0; c < armci_nproc; c++) {
       armci_connect_t *con = CLN_con + c;
       qp_attr_mask |= IBV_QP_DEST_QPN | IBV_QP_AV;
       qp_attr.dest_qp_num  = SRV_rqpnums[c];
       qp_attr.ah_attr.dlid = SRV_nic->lid_arr[c];
       qp_attr.ah_attr.port_num = CLN_nic->active_port;

       rc = ibv_modify_qp(con->qp, &qp_attr, qp_attr_mask); 
       dassertp(1,!rc,("%d: INIT->RTR server cln=%d rc=%d\n",armci_me,c,rc));
    }
    }
    armci_util_wait_int(&armci_vapi_client_ready,1,10);
    memset(&qp_attr, 0, sizeof qp_attr);

    if (!armci_use_odcm) {
    qp_attr_mask = IBV_QP_STATE
                 | IBV_QP_SQ_PSN
                 | IBV_QP_TIMEOUT
                 | IBV_QP_RETRY_CNT
                 | IBV_QP_RNR_RETRY
                 | IBV_QP_MAX_QP_RD_ATOMIC;

    qp_attr.qp_state            = IBV_QPS_RTS;
    qp_attr.sq_psn              = 0;
    qp_attr.timeout             = 18;
    qp_attr.retry_cnt           = 7;
    qp_attr.rnr_retry           = 7;
    qp_attr.max_rd_atomic  = 4;

    for (c = 0; c < armci_nproc; c++) {
       armci_connect_t *con = CLN_con + c;
       rc = ibv_modify_qp(con->qp, &qp_attr,qp_attr_mask);
       dassertp(1,!rc,("%d: server RTR->RTS cln=%d rc=%d\n",armci_me,c,rc));
    }
    free(SRV_rqpnums);
    SRV_rqpnums = NULL;
    }

    armci_server_alloc_bufs();

    if (!armci_use_odcm) {
    /* setup descriptors and post nonblocking receives */
#if defined(PEND_BUFS)
    assert(armci_nproc*(IMM_BUF_NUM+1)<DSCRID_IMMBUF_RECV_END-DSCRID_IMMBUF_RECV);
    for(i =  0; i < armci_nproc; i++) {
      for(j=0; j<IMM_BUF_NUM+1; j++) {
	vapibuf_t *cbuf;
	cbuf = serv_buf_arr[i*(IMM_BUF_NUM+1)+j];
	armci_init_vapibuf_recv(&cbuf->dscr, &cbuf->sg_entry, cbuf->buf,
				IMM_BUF_LEN, &serv_memhandle);
	/* we use index of the buffer to identify the buffer, this index is
	 * returned with a call to ibv_poll_cq inside the ibv_wr */
	cbuf->dscr.wr_id = i*(IMM_BUF_NUM+1)+j + DSCRID_IMMBUF_RECV;
	if (DEBUG_SERVER) {
	  printf("\n%d(s):posted rr with lkey=%d",armci_me,cbuf->sg_entry.lkey);
	  fflush(stdout);
	}
	rc = ibv_post_recv((CLN_con+i)->qp, &cbuf->dscr, &bad_wr);
	dassert1(1,rc==0,rc);
      }
    }
#else
    for(i =  0; i < armci_nproc; i++) {
      vapibuf_t *cbuf;
      cbuf = serv_buf_arr[i];
      armci_init_vapibuf_recv(&cbuf->dscr, &cbuf->sg_entry, cbuf->buf,
			      CBUF_DLEN, &serv_memhandle);
      /* we use index of the buffer to identify the buffer, this index is
       * returned with a call to ibv_poll_cq inside the ibv_wr */
      cbuf->dscr.wr_id = i+armci_nproc;
      if (DEBUG_SERVER) {
	printf("\n%d(s):posted rr with lkey=%d",armci_me,cbuf->sg_entry.lkey);
	fflush(stdout);
      }
      rc = ibv_post_recv((CLN_con+i)->qp, &cbuf->dscr, &bad_wr);
      dassert1(1,rc==0,rc);
    }
#endif
    }
    if (TIME_INIT) printf("\n%d:post time for server_initial_conn is %f",
                          armci_me, MPI_Wtime() - inittime4);

    armci_vapi_server_ready=1;

    if (DEBUG_SERVER) {
       printf("%d: server connected to all clients\n",armci_me); fflush(stdout);
    }
    if (TIME_INIT) printf("\n%d:time for server_initial_conn is %f",
                          armci_me, MPI_Wtime() - inittime0);
}

static void armci_finalize_nic(vapi_nic_t *nic)
{
    int ret;

    ret = ibv_destroy_cq(nic->scq);
    dassert1(1,ret==0,ret);
    armci_vapi_check_return(DEBUG_FINALIZE,ret,"armci_finalize_nic:destroy_scq");

    ret = ibv_destroy_comp_channel(nic->sch);
    dassert1(1,ret==0,ret);
    armci_vapi_check_return(DEBUG_FINALIZE,ret,"armci_finalize_nic:destroy_sch");

    ret = ibv_destroy_cq(nic->rcq);
    dassert1(1,ret==0,ret);
    armci_vapi_check_return(DEBUG_FINALIZE,ret,"armci_finalize_nic:destroy_rcq");

    ret = ibv_destroy_comp_channel(nic->rch);
    dassert1(1,ret==0,ret);
    armci_vapi_check_return(DEBUG_FINALIZE,ret,"armci_finalize_nic:destroy_rch");

    ret = ibv_close_device(nic->handle);
    dassert1(1,ret==0,ret);

    armci_vapi_check_return(DEBUG_FINALIZE,ret,"armci_finalize_nic:release_hca");

}


void armci_server_transport_cleanup()
{
    int s;
    int rc;

    /*first we have empty send/recv queues TBD*/
    if(serv_malloc_buf_base){
        rc = ibv_dereg_mr(serv_memhandle.memhndl);
	dassert1(1,rc==0,rc);
        armci_vapi_check_return(DEBUG_FINALIZE,rc,
                                "armci_server_transport_cleanup:deregister_mr");
       /*now free it*/
       free(serv_malloc_buf_base);
    }
    /*now deregister all my regions from regionskk.c*/
    armci_server_region_destroy();
    if (CLN_con) {
        for (s = 0; s < armci_nproc; s++) {
            armci_connect_t *con = CLN_con + s;
            if (con->qp) {
                rc = ibv_destroy_qp(con->qp);
                armci_vapi_check_return(DEBUG_FINALIZE,rc,
                                        "armci_server_transport_cleanup:destroy_qp");
            }
#if 0
            free(con->rqpnum);
#endif
        }
        free(CLN_con);
    }
    armci_finalize_nic(CLN_nic);
}

void armci_transport_cleanup()
{
    int s;
    int rc;

    /*first deregister buffers memory */
    if (client_malloc_buf_base) {
        rc = ibv_dereg_mr(client_memhandle.memhndl);
	dassert1(1,rc==0,rc);
        armci_vapi_check_return(DEBUG_FINALIZE,rc,"armci_client_transport_cleanup:deregister_mr");
        /*now free it*/
        free(client_malloc_buf_base);
    }
    /*now deregister all my regions from regions.c*/
    armci_region_destroy();
    if (SRV_con) {
        for (s = 0; s < armci_nclus; s++) {
            armci_connect_t *con = SRV_con + s;
            if (con->qp) {
                rc = ibv_destroy_qp(con->qp);
		dassert1(1,rc==0,rc);
                armci_vapi_check_return(DEBUG_FINALIZE,rc,"armci_client_transport_cleanup:destroy_qp");
            }
#if 0
            free(con->rqpnum);
#endif
        }
        free(SRV_con);
    }
    armci_finalize_nic(SRV_nic);
}

/** Post an immediate buffer back for the client to send.
 */
static void _armci_pendbuf_post_immbuf(vapibuf_t *cbuf, int to) {
  int rc; 
  struct ibv_recv_wr *bad_wr;
#if defined(PEND_BUFS)
  assert(cbuf->dscr.wr_id == cbuf-serv_buf_arr[0]+DSCRID_IMMBUF_RECV);
#endif
  rc = ibv_post_recv((CLN_con+to)->qp, &(cbuf->dscr), &bad_wr);
  dassert1(1,rc==0,rc);
}

#if defined(PEND_BUFS)
#define DSCRID_TO_IMMBUFID(x) (x-DSCRID_IMMBUF_RECV)
#else
#define DSCRID_TO_IMMBUFID(x) ((x)-armci_nproc)
#endif

#if defined(PEND_BUFS)

/**Obtain a message receive buffer to receive a message. Used in place
 *  of MessageRcvBuffer. Should not be used.
 */
char *armci_openib_get_msg_rcv_buf(int proc)
{
  armci_die("PEND_BUFS in OPENIB: MessageRcvBuffer not available. Should use the in-place buffers to receive data", proc);
  return NULL;
}

/** Check that the data is in a server allocated buffer. This is
 *   guaranteed to be pinned. Ideally, this should always be true. Any
 *   operation that request alternative support will have to fix this
 *   function and possibly @armci_openib_get_msg_rcv_buf().
 * @param br IN Buffer pointer being checked
 * @return 1 if it is a server-allocated buffer. 0 otherwise.
 */
int armci_data_in_serv_buf(void *br)
{
  if(br>=(void *)serv_malloc_buf_base && br<(void *)serv_tail)
    return 1;
  if(DEBUG_SERVER) {
    printf("%d:: serv_bufs=%p<->%p. br=%p out of range\n",
	   armci_me, serv_malloc_buf_base, serv_tail, br); 
    fflush(stdout);
  }
  return 0;
}

#define PBUF_BUFID_TO_PUT_WRID(_pbufid) (DSCRID_PENDBUF+(_pbufid)*2)
#define PBUF_BUFID_TO_GET_WRID(_pbufid) (DSCRID_PENDBUF+(_pbufid)*2+1)
#define PBUF_WRID_TO_PBUFID(_id) (((_id)-DSCRID_PENDBUF)/2)
#define PBUF_IS_GET_WRID(_id) (((_id)-DSCRID_PENDBUF)&1)
#define PBUF_IS_PUT_WRID(_id) (!(((_id)-DSCRID_PENDBUF)&1))

/**Complete processing this immediate buffer. Parameters is void *,
 * since vapibuf_t*|immbuf_t* is not available in armci-vapi.h
 */
void armci_complete_immbuf(void *buf) {
  vapibuf_t *cbuf = (vapibuf_t*)buf;
  request_header_t *msginfo=(request_header_t*)cbuf->buf;
  
#if SRI_CORRECT
#error
  cbuf->send_pending = 0;
#else
    _armci_pendbuf_post_immbuf(cbuf,msginfo->from);
#endif
  armci_data_server(cbuf);
  if(msginfo->operation==PUT || ARMCI_ACC(msginfo->operation)) {
    SERVER_SEND_ACK(msginfo->from);
  }  
#if SRI_CORRECT
  if(!cbuf->send_pending) {
    _armci_pendbuf_post_immbuf(cbuf,msginfo->from);
  }
#endif
}

/**Complete processing this pending buffer. Parameters is void *,
 * since vapibuf_t*|immbuf_t* is not available in armci-vapi.h. Note
 * that the pending buffer may not yet be available for reuse. This
 * will depend on the state of the pending buffer (which might have to
 * wait for a communication innitiated by armci_data_server() to
 * complete.
 */
void armci_complete_pendbuf(void *buf) {
  vapibuf_pend_t *pbuf = (vapibuf_pend_t *)buf;
  request_header_t *msginfo=(request_header_t*)pbuf->buf;

  assert(pbuf->vbuf);
#if SRI_CORRECT
  pbuf->cbuf->send_pending=0;
#else
  _armci_pendbuf_post_immbuf(pbuf->vbuf,msginfo->from);
#endif
  armci_data_server(pbuf);
  if(msginfo->operation==PUT || ARMCI_ACC(msginfo->operation)) {
    SERVER_SEND_ACK(msginfo->from);
  }
#if SRI_CORRECT
#error  
 assert(!pbuf->cbuf->send_pending);
  _armci_pendbuf_post_immbuf(pbuf->cbuf,msginfo->from);
#endif
}

void _armci_get_data_from_client(int proc, struct ibv_send_wr *sdscr, 
				 int dscrid, struct ibv_sge *ssg_entry, 
				 void *rbuf, void *lbuf, int bytes) ;
void _armci_send_data_to_client_pbuf(int proc, struct ibv_send_wr *sdscr, 
				     int dscrid, struct ibv_sge *ssg_entry, 
				     void *rbuf, void *lbuf, int bytes);

int no_srv_copy_nsegs_ulimit() {
  return armci_max_qp_ous_swr*armci_max_num_sg_ent/10;
}

/** Initiate a get operation to progress a pending buffer.
 * @param msginfo Request header for any additional processing
 * @param src Pointer to src of data (remote for GET)
 * @param dst Pointer to dst
 * @param bytes #bytes to transfer
 * @param proc proc to transfer from(for get)/to(for put)
 * @param pbufid Index of pending buffer
 */
void armci_pbuf_start_get(void *msg_info, void *src, void *dst, 
			  int bytes, int proc, int pbufid) {
  struct ibv_send_wr sdscr;
  struct ibv_sge sg_entry;
  int wrid = PBUF_BUFID_TO_GET_WRID(pbufid);
  request_header_t *msginfo=(request_header_t *)msg_info;
  void armci_server_rdma_contig_to_strided(char *src_ptr, int proc,
					   char *dst_ptr, 
					   int dst_stride_arr[],
					   int seg_count[],
					   int stride_levels,
					   request_header_t *msginfo);


#if defined(PUT_NO_SRV_COPY)
  if(msginfo->operation==PUT && msginfo->format==STRIDED 
     && !msginfo->pinned && src==msginfo->tag.data_ptr)   {
    char *loc_ptr, *rem_ptr;
    int stride_levels, *count;
    int *loc_stride_arr;
    char *dscr = (char *)(msginfo+1);
    ARMCI_MEMHDL_T *mhloc=NULL;
    int nsegs, i;

    /* unpack descriptor record */
    loc_ptr = *(void**)dscr;           dscr += sizeof(void*);
    stride_levels = *(int*)dscr;       dscr += sizeof(int);
    loc_stride_arr = (int*)dscr;       dscr += stride_levels*sizeof(int);
    count = (int*)dscr;

    rem_ptr = msginfo->tag.data_ptr;

    nsegs = 1;
    for(i=0; i<stride_levels; i++) 
      nsegs *= count[i+1];    

    dassert(1,proc==msginfo->from);
    if(nsegs<no_srv_copy_nsegs_ulimit() &&
       get_armci_region_local_hndl(loc_ptr,armci_clus_id(armci_me),&mhloc)) {
/*       printf("%d(s): direct rdma from client buffers to server-side memory\n",armci_me); */
/*       fflush(stdout); */
   
      armci_server_rdma_contig_to_strided(rem_ptr, proc,
					  loc_ptr,loc_stride_arr,
					  count, stride_levels,
					  msginfo);
    return;
   }
  }
#endif
/*   printf("%d(s): rdma from client buffers to pending buffers\n",armci_me); */
/*   fflush(stdout);   */
  _armci_get_data_from_client(proc,&sdscr,wrid,&sg_entry,src,dst,bytes);
}

/** Initiate a put operation to progress a pending buffer.
 * @param src Pointer to src of data (local for PUT)
 * @param dst Pointer to dst
 * @param bytes #bytes to transfer
 * @param proc proc to transfer from(for get)/to(for put)
 * @param pbufid Index of pending buffer
 */
void armci_pbuf_start_put(void *src, void *dst, int bytes, int proc, 
			  int pbufid) {
  struct ibv_send_wr sdscr;
  struct ibv_sge sg_entry;
  int wrid = PBUF_BUFID_TO_PUT_WRID(pbufid);

  _armci_send_data_to_client_pbuf(proc,&sdscr,wrid,&sg_entry,src,dst,bytes);
}

/**
  * function to get data from remote client called by data
  * server. Note that this is only called for pending buffers.
  * @param proc IN the id of remote client
  * @param sdscr IN/OUT Descriptor to be used to post the get
  * @param dscrid IN ID to be used for the descriptor
  * @param ssg_entry IN Scatter/gather list
  * @param rbuf IN the remote buffer to get from
  * @param lbuf IN local buf to get the data into, this is the queue buffer for SERVER_QUEUE path
  * @param bytes IN the size of get
  * @see SERVER_QUEUE
  * @see armci_send_data_to_client
  */
/*static*/ void _armci_get_data_from_client(int proc, struct ibv_send_wr *sdscr, 
				int dscrid, struct ibv_sge *ssg_entry, 
				void *rbuf, void *lbuf, int bytes) 
{
    int rc = 0;

    if(DEBUG_SERVER){
       printf("\n%d(s):sending data to client %d at %p flag = %p bytes=%d\n",
               armci_me,
               proc,lbuf,(char *)lbuf+bytes-sizeof(int),bytes);fflush(stdout);
    }

    memset(sdscr,0,sizeof(struct ibv_send_wr));
    armci_init_cbuf_rrdma(sdscr,ssg_entry,lbuf,rbuf,bytes,
                          &serv_memhandle,(handle_array+proc));

    if(DEBUG_SERVER){
       printf("\n%d(s):handle_array[%d]=%p lbuf=%p flag=%p bytes=%d\n",armci_me,
              proc,&handle_array[proc],(char *)lbuf,
              (char *)lbuf+bytes-sizeof(int),bytes);
       fflush(stdout);
    }

    assert(sizeof(request_header_t)+bytes<PENDING_BUF_LEN);

    sdscr->wr_id = dscrid;
    struct ibv_send_wr *bad_wr;
    rc = ibv_post_send((CLN_con+proc)->qp, sdscr, &bad_wr);
    dassert1(1,rc==0,rc);
}

void _armci_send_data_to_client_pbuf(int proc, struct ibv_send_wr *sdscr, 
				     int dscrid, struct ibv_sge *ssg_entry, 
				     void *rbuf, void *lbuf, int bytes)  {
    int rc = 0;

    if(DEBUG_SERVER) {
       printf("\n%d(s):sending data to client %d at %p flag = %p bytes=%d\n",
               armci_me,
               proc,rbuf,(char *)rbuf+bytes-sizeof(int),bytes);fflush(stdout);
    }
    memset(sdscr,0,sizeof(struct ibv_send_wr));
    armci_init_cbuf_srdma(sdscr,ssg_entry,lbuf,rbuf,bytes,
                          &serv_memhandle,(handle_array+proc));
    if(DEBUG_SERVER){
       printf("\n%d(s):handle_array[%d]=%p dbuf=%p flag=%p bytes=%d\n",armci_me,
              proc,&handle_array[proc],(char *)rbuf,
              (char *)rbuf+bytes-sizeof(int),bytes);
       fflush(stdout);
    }
    sdscr->wr_id = dscrid;
    struct ibv_send_wr *bad_wr;
    rc = ibv_post_send((CLN_con+proc)->qp, sdscr, &bad_wr);
    dassert1(1,rc==0,rc);
}
#endif

#define DATA_SERVER_YIELD_CPU
void armci_call_data_server()
{
int rc = 0;
int rc1 = 0;
vapibuf_t *cbuf,*cbufs;
request_header_t *msginfo,*msg;
int c,i,need_ack,pollcount;
static int mytag=1;
int rrr,serverwcount=0;

#ifdef CHANGE_SERVER_AFFINITY
cpu_set_t mycpuid,new_mask;
char str[CPU_SETSIZE];
char cid[8];
extern char * cpuset_to_cstr(cpu_set_t *mask, char *str);
int nslave=armci_clus_info[armci_clus_me].nslave;
    rrr=sched_getaffinity(0, sizeof(mycpuid), &mycpuid);
#endif

#if ARMCI_ENABLE_GPC_CALLS
    unblock_thread_signal(GPC_COMPLETION_SIGNAL);
#endif
#if defined(PEND_BUFS)
    armci_pendbuf_init(); 
#endif

    for (;;) {
      struct ibv_wc *pdscr=NULL;
      struct ibv_wc pdscr1;
      pdscr = &pdscr1;
      pdscr->status = IBV_WC_SUCCESS;
      rc = 0;
#ifdef CHANGE_SERVER_AFFINITY
      static int ccc;
      serverwcount++;
      if(serverwcount==100){
        serverwcount=0;
        ccc=(ccc+1)%nslave;
        sprintf (cid, "%d", ccc);
        rrr = cstr_to_cpuset(&new_mask,cid);
        if (sched_setaffinity(0, sizeof (new_mask), &new_mask)) {
          perror("sched_setaffinity");
          printf("failed to set pid %d's affinity.\n", getpid());
        }
        rrr=sched_getaffinity(0, sizeof(mycpuid), &mycpuid);
        if(rrr)perror("sched_getaffinity");
      }
#else
#ifdef DATA_SERVER_YIELD_CPU_
      serverwcount++;
      if(serverwcount==50){
        serverwcount=0;usleep(1);
      }
#endif
#endif

#if ARMCI_ENABLE_GPC_CALLS
      block_thread_signal(GPC_COMPLETION_SIGNAL);
#endif
      bzero(pdscr, sizeof(*pdscr));
      do {
        rc = ibv_poll_cq(CLN_nic->rcq, 1, pdscr);
	if (armci_server_terminating) {
	  /* server is interrupted when clients terminate connections */
	  armci_server_transport_cleanup();
	  sleep(1);
	  _exit(0);
	}
	if (rc == 0 && !armci_openib_server_poll) {
	  /* wait for a notify event */
          rc1 = ibv_get_cq_event(CLN_nic->rch,&CLN_nic->rcq,&CLN_nic->rcq_cntx);
          dassert1(1,rc1==0,rc1);
          ibv_ack_cq_events(CLN_nic->rcq, 1);
	  /* re-arm notify event */
          rc1 = ibv_req_notify_cq(CLN_nic->rcq, 0);
          dassert1(1,rc1==0,rc1);
	  /* note: an event receive does not guarantee an actual completion */
	  continue;
	}
      } while (rc == 0);

      if(DEBUG_SERVER) {
        printf("\n%d:pdscr=%p %p %d %d %d %d\n",armci_me,pdscr,&pdscr1,
                           pdscr->status,pdscr->opcode,pdscr->vendor_err,
                           pdscr->src_qp);
        fflush(stdout);
      }
      dassertp(1,rc>=0,("%d: rc=%d id=%d status=%d",
			armci_me,rc,(int)pdscr->wr_id,pdscr->status));
      dassert1(1,pdscr->status==IBV_WC_SUCCESS,pdscr->status);
                            
       if (DEBUG_SERVER) {
         printf("%d(s) : NEW MESSAGE bytelen %d \n",armci_me,pdscr->byte_len);
         printf("%d(s) : NEW MESSAGE id is %ld \n",armci_me,pdscr->wr_id);
         fflush(stdout);
       }
#if defined(PEND_BUFS)
      if(pdscr->wr_id>=DSCRID_IMMBUF_RESP && pdscr->wr_id<DSCRID_IMMBUF_RESP_END) {
/* 	fprintf(stderr, "%d(s) : Got server response msg completion\n", armci_me); */
#if SRI_CORRECT
	int id = pdscr->wr_id - DSCRID_IMMBUF_RESP;
	if(id>=0 && id<armci_nproc*(IMM_BUF_NUM+1)) {
	  int dest = id/(IMM_BUF_NUM+1);
	  dassert(1,serv_buf_arr[id]->send_pending==1);
	  serv_buf_arr[id]->send_pending = 0;
	  _armci_pendbuf_post_immbuf(serv_buf_arr[id],dest);
	}
#endif
	continue;
      }
       if (pdscr->wr_id>=DSCRID_PENDBUF && pdscr->wr_id<DSCRID_PENDBUF_END) {
	 int pbufid = PBUF_WRID_TO_PBUFID(pdscr->wr_id);
/* 	 printf("%d(s) : Progressing pending msg (something completed) pbufid=%d id=%ld byte_len=%d status=%d\n", armci_me, pbufid,pdscr->wr_id,pdscr->byte_len,done_status); */
/* 	 fflush(stdout); */
	 if(PBUF_IS_GET_WRID(pdscr->wr_id))
	   armci_pendbuf_done_get(pbufid);
	 else if(PBUF_IS_PUT_WRID(pdscr->wr_id))
	   armci_pendbuf_done_put(pbufid);
	 else
	   armci_die("Pending buffer op completed. But not PUT or GET!",pdscr->wr_id);
	 continue;
       }
#endif
       if (pdscr->wr_id >= DSCRID_SCATGAT && pdscr->wr_id < DSCRID_SCATGAT_END) {
	 sr_descr_t *sdscr_arr, *rdscr_arr;
         if (DEBUG_SERVER) {
           printf("%d(s) : received SCATGAT DATA id = %ld, length = %d\n",
                  armci_me,pdscr->wr_id, pdscr->byte_len);
           fflush(stdout);
	 }
#if defined(PEND_BUFS)
	 sdscr_arr = armci_vapi_serv_nbsdscr_array;
	 assert(sdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].numofsends>0);
	 sdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].numofsends--;
	 if(sdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].numofsends==0)
	     sdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].tag=0;
#else
	 rdscr_arr;
	 rdscr_arr = armci_vapi_serv_nbrdscr_array;
	 rdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].numofrecvs--;
	 if(rdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].numofrecvs==0)
	   rdscr_arr[pdscr->wr_id-DSCRID_SCATGAT].tag=0;
#endif
         continue;
       }

#if defined(PEND_BUFS)
       assert(pdscr->wr_id>=DSCRID_IMMBUF_RECV && pdscr->wr_id<DSCRID_IMMBUF_RECV_END);
#endif
       cbuf = serv_buf_arr[DSCRID_TO_IMMBUFID(pdscr->wr_id)];
       assert(cbuf->dscr.wr_id == pdscr->wr_id);

       msginfo = (request_header_t*)cbuf->buf;
       armci_ack_proc = c = msginfo->from;

       if (DEBUG_SERVER) {
         printf("%d(s) : request id is %ld operation is %d, length is %d from=%d cbuf->dscr.wr_id=%d\n",
		armci_me,pdscr->wr_id,msginfo->operation,pdscr->byte_len,msginfo->from, (int)cbuf->dscr.wr_id);
         fflush(stdout);
       }

#if defined(PEND_BUFS)
       cbufs = cbuf;
       armci_init_vapibuf_recv(&cbufs->dscr, &cbufs->sg_entry,cbufs->buf,
			       IMM_BUF_LEN, &serv_memhandle);
       cbufs->dscr.wr_id = pdscr->wr_id;
#else
       cbufs = serv_buf_arr[pdscr->wr_id - armci_nproc] = spare_serv_buf;
       armci_init_vapibuf_recv(&cbufs->dscr, &cbufs->sg_entry,cbufs->buf,
			       CBUF_DLEN, &serv_memhandle);
       cbufs->dscr.wr_id = c + armci_nproc;

       spare_serv_buf = cbuf;
#endif

       if(DEBUG_SERVER) {
	 printf("%d(s):Came out of poll id=%ld\n",armci_me,pdscr->wr_id);
	 fflush(stdout);
       }

       if(msginfo->operation == PUT &&msginfo->pinned == 1){
	 int found, num;
	 int stride_arr[MAX_STRIDE_LEVEL]; /*should be MAX_STRIDE_LEVELS*/
	 int count[MAX_STRIDE_LEVEL];
	 void *dest_ptr;
	 int stride_levels;
	 ARMCI_MEMHDL_T *loc_memhandle;
	 void armci_post_scatter(void *,int *,int *,int, armci_vapi_memhndl_t *,int,int,int,sr_descr_t **);

	 /*unpack decsriptor_record : should call a function instead */
	 msg = msginfo + 1;
	 test_ptr = dest_ptr = *(void**)msg;
	 msg = (request_header_t *) ((char*)msg + sizeof(void*));
	 test_stride_levels=stride_levels = *(int*)msg;
	 msg = (request_header_t *) ((char*)msg + sizeof(int));
	 for(i =0; i<stride_levels; i++){
	   test_stride_arr[i] = stride_arr[i] = *(int*)msg;
	   msg = (request_header_t*) ((int*)msg + 1);
	 }
	 for(i=0; i<stride_levels+1; i++){
	   test_count[i] = count[i] = *(int*)msg;
	   msg = (request_header_t*) ((int*)msg + 1);
	 }

	 if (DEBUG_SERVER) {
	   printf(" server:the dest_ptr is %p\n", dest_ptr);
	   for(i =0; i<stride_levels; i++)
	     printf("stride_arr[i] is %d,value of count[i] is %d\n",
		    stride_arr[i], count[i]);
	   printf("the value of stride_levels is %d\n", stride_levels);
	   fflush(stdout);
	 }

	 found =get_armci_region_local_hndl(dest_ptr,armci_me, &loc_memhandle);
	 dassertp(1,found!=0,("%d:SERVER : local region not found id=%d",
			      armci_me,pdscr->wr_id));

	 if(DEBUG_SERVER) {
	   printf("%d(s) : about to call armci_post_scatter\n",armci_me);
	   fflush(stdout);
	 }

	 armci_post_scatter(dest_ptr, stride_arr, count, stride_levels,
			    loc_memhandle,msginfo->from, mytag, SERV,NULL );

	 mytag = (mytag+1)%(MAX_PENDING);
	 if(mytag==0)mytag=1;

	 if(DEBUG_SERVER) {
	   printf("%d(s) : finished posting %d scatter\n",armci_me,num);
	   fflush(stdout);
	 }
	 _armci_pendbuf_post_immbuf(cbufs, msginfo->from);
	 SERVER_SEND_ACK(msginfo->from);
	 need_ack = 0;
       }
       else if(msginfo->operation == REGISTER){
	 if (DEBUG_SERVER) {
            printf("%d(s) : Register_op id is %d, comp_dscr_id is  %ld\n",
                     armci_me,msginfo->operation,pdscr->wr_id);
            fflush(stdout);
          }

          armci_server_register_region(*((void **)(msginfo+1)),
                           *((long *)((char *)(msginfo+1)+sizeof(void *))),
                           (ARMCI_MEMHDL_T *)(msginfo->tag.data_ptr));
	  _armci_pendbuf_post_immbuf(cbufs, msginfo->from);
          *(long *)(msginfo->tag.ack_ptr) = ARMCI_STAMP;
          continue;
       }
       else {
         if(DEBUG_SERVER) {
	   printf("%d(s) : request is %ld about to call armci_data_server\n",
		  armci_me, pdscr->wr_id);
	   fflush(stdout);
         }
#if defined(PEND_BUFS)
	 armci_pendbuf_service_req(cbuf);
#else
	 _armci_pendbuf_post_immbuf(cbufs, msginfo->from);
	 armci_data_server(cbuf);
	 
	 if((msginfo->operation == PUT) || ARMCI_ACC(msginfo->operation)) { 
	   /* for operations that do not send data back we can send ACK now */
	   SERVER_SEND_ACK(msginfo->from);
	   need_ack=0;
	   if(DEBUG_SERVER){
	     printf("%d(s) : posted ack\n\n",armci_me);
	     fflush(stdout);
	   }
	 } else need_ack=1;
#endif
       }
       if (0) {
	 printf("%d(s):Done processed request\n\n",armci_me);
	 fflush(stdout);
       }
       
#if ARMCI_ENABLE_GPC_CALLS
       unblock_thread_signal(GPC_COMPLETION_SIGNAL);
#endif
   }/* end of for */
}


void armci_vapi_complete_buf(armci_vapi_field_t *field,int snd,int rcv,int to,int op) {
  struct ibv_send_wr *snd_dscr;

  BUF_INFO_T *info;
  info = (BUF_INFO_T *)((char *)field-sizeof(BUF_INFO_T));

  if(info->tag && op==GET)return;

  if(snd){
    request_header_t *msginfo = (request_header_t *)(field+1);
    snd_dscr=&(field->sdscr);
    if(mark_buf_send_complete[snd_dscr->wr_id]==0)
      armci_send_complete(snd_dscr,"armci_vapi_complete_buf",1);
  }

  if(rcv){
    int *last;
    long *flag;
    int loop = 0;
    request_header_t *msginfo = (request_header_t *)(field+1);
    flag = (long *)&msginfo->tag.ack;

  
    if(op==PUT || ARMCI_ACC(op)){
      if(msginfo->bypass && msginfo->pinned && msginfo->format == STRIDED &&
	 op == PUT);
      else{
	while(armci_util_long_getval(flag) != ARMCI_STAMP) {
	  loop++;
	  loop %=100000;
	  if(loop==0){
	  }
	}
      }
      /* 	 printf("%d: client complete_buf. op=%d loop=%d till *flag=ARMCI_STAMP\n", armci_me,op,loop); */
      /* 	 fflush(stdout); */
      *flag = 0L;
    }
    else{
      /*SK: I think we get here only for GET with result directly
	going to client's pinned memory. (info.tag==0 && op==GET)*/
      last = (int *)((char *)msginfo+msginfo->datalen-sizeof(int));
      while(armci_util_int_getval(last) == ARMCI_STAMP &&
	    armci_util_long_getval(flag)  != ARMCI_STAMP){
	loop++;
	loop %=100000;
	if(loop==0){
	  if(DEBUG_CLN){
	    printf("%d: client last(%p)=%d flag(%p)=%ld off=%d\n",
		   armci_me,last,*last,flag,*flag,msginfo->datalen);
	    fflush(stdout);
	  }
	}
      }
    }
  }
}

void armci_vapi_test_buf(armci_vapi_field_t *field,int snd,int rcv,int to,int op, int *retval) {
  struct ibv_send_wr *snd_dscr;

  BUF_INFO_T *info;
  info = (BUF_INFO_T *)((char *)field-sizeof(BUF_INFO_T));

  *retval = 0;

  if(info->tag && op==GET)return;

  if(snd){
    request_header_t *msginfo = (request_header_t *)(field+1);
    snd_dscr=&(field->sdscr);
    if(mark_buf_send_complete[snd_dscr->wr_id]==0) {
/*       printf("%d: test buf. send not complete\n",armci_me); */
/*       fflush(stdout); */
      return;
    }
  }

  if(rcv){
    int *last;
    long *flag;
    int loop = 0;
    request_header_t *msginfo = (request_header_t *)(field+1);
    flag = (long *)&msginfo->tag.ack;
  
    if(op==PUT || ARMCI_ACC(op)){
      if(msginfo->bypass && msginfo->pinned && msginfo->format == STRIDED &&
	 op == PUT) 
	*retval=1;
      else{
	if(armci_util_long_getval(flag) == ARMCI_STAMP) {
	  *retval = 1;
	}
      }
      return;
    }
    else{
      /*SK: I think we get here only for GET with result directly
	going to client's pinned memory. (info.tag==0 && op==GET)*/
      last = (int *)((char *)msginfo+msginfo->datalen-sizeof(int));
      if(armci_util_int_getval(last) != ARMCI_STAMP ||
	    armci_util_long_getval(flag)  == ARMCI_STAMP){
	*retval=1;
      }
      return;
    }
  }
}


static inline void armci_vapi_post_send(int isclient,int con_offset,
                                        struct ibv_send_wr *snd_dscr,char *from)
{
    int rc = 0;
    vapi_nic_t *nic;
    armci_connect_t *con;
    int total = 0;

    if(!isclient){
       nic = CLN_nic;
       con = CLN_con+con_offset;
    }
    else{
       nic = SRV_nic;
       con = SRV_con+con_offset;
    }

    if(DEBUG_CLN){
       printf("vapi_post_send: snd_dscr->num_sge=%d, snd_dscr->sg_list->length=%d\n",
              snd_dscr->num_sge, snd_dscr->sg_list->length);
       fflush(stdout);
    }


    /* find the total length of all the segments */
    total = snd_dscr->sg_list->length * snd_dscr->num_sge;
    if(DEBUG_CLN){
       printf("%d(c) : total is %d\t, max_size is %d\n",armci_me,total,
                    armci_vapi_max_inline_size);
    }

    struct ibv_send_wr *bad_wr;
    if (total > armci_vapi_max_inline_size) {
        rc = ibv_post_send(con->qp, snd_dscr, &bad_wr);
    } else {
        rc = ibv_post_send(con->qp, snd_dscr, &bad_wr);
        /* no corresponding call, using ibv_post_send
       rc = EVAPI_post_inline_sr(nic->handle,con->qp,snd_dscr);*/
    }
    dassert1(1,rc==0,rc);
}

/** Send request to server. 
  */
int armci_send_req_msg(int proc, void *buf, int bytes)
{
  int cluster = armci_clus_id(proc), i;
    request_header_t *msginfo = (request_header_t *)buf;
    struct ibv_send_wr *snd_dscr;
    struct ibv_sge *ssg_lst;

    THREAD_LOCK(armci_user_threads.net_lock);   

    check_state_of_ib_connection(proc, 0);

    snd_dscr = BUF_TO_SDESCR((char *)buf);
    ssg_lst  = BUF_TO_SSGLST((char *)buf);

    /*Stamp end of buffers as needed*/
    if(msginfo->operation == GET && !msginfo->pinned) {
      const int dscrlen = msginfo->dscrlen;
      const int datalen = msginfo->datalen;
      int *last;
      if(dscrlen < (datalen - sizeof(int)))
	last = (int*)(((char*)(msginfo+1))+(datalen-sizeof(int)));
      else
	last = (int*)(((char*)(msginfo+1))+(dscrlen+datalen-sizeof(int)));
      *last = ARMCI_STAMP;
#ifdef GET_STRIDED_COPY_PIPELINED
      if(msginfo->format == STRIDED) {
	const int ssize = GET_STRIDED_COPY_PIPELINED_SIZE/sizeof(int);
	int *sfirst = (int*)(dscrlen+(char*)(msginfo+1))+ssize; /*stamping
							    can start here*/
	int *slast = last, *ptr;
	for(ptr=sfirst; ptr<slast; ptr+=ssize) {
	  *ptr = ARMCI_STAMP;
	}
      }
#endif
    }
    if(msginfo->operation == ACK) {
      *(int *)(msginfo +1) = ARMCI_STAMP+1;
      *(((int *)(msginfo +1))+1) = ARMCI_STAMP+1;
    }
    

#if defined(PEND_BUFS)
    if((msginfo->operation==PUT || ARMCI_ACC(msginfo->operation)) 
       && bytes > IMM_BUF_LEN) {
      msginfo->tag.imm_msg=0;
      assert(sizeof(request_header_t)<IMM_BUF_LEN); /*sanity check*/
      bytes = ARMCI_MIN(bytes-msginfo->datalen, IMM_BUF_LEN);
      assert(bytes==IMM_BUF_LEN||(bytes==sizeof(*msginfo)+msginfo->dscrlen));
    }
    else if(msginfo->operation==GET
	    && !(msginfo->datalen+sizeof(request_header_t)+msginfo->dscrlen<IMM_BUF_LEN)) {
      assert(sizeof(request_header_t) < IMM_BUF_LEN);
      msginfo->tag.imm_msg=0;
      bytes = ARMCI_MIN(sizeof(request_header_t)+msginfo->dscrlen, IMM_BUF_LEN);
    }
#if defined(PUT_NO_SRV_COPY) && 0 /*SK:disabled. Imm msgs are sent inline
				    for latency reasons*/
    else if(msginfo->operation==PUT && !msginfo->pinned && msginfo->format==STRIDED && msginfo->tag.data_len>=2048) {
      msginfo->tag.imm_msg = 0;
      assert(sizeof(request_header_t)<IMM_BUF_LEN); /*sanity check*/
      bytes = ARMCI_MIN(bytes-msginfo->datalen, IMM_BUF_LEN);
      assert(bytes==IMM_BUF_LEN||(bytes==sizeof(*msginfo)+msginfo->dscrlen));
    }
#endif
    else{
      msginfo->tag.imm_msg=1;
    }
/*    printf("%d: send_req: op=%d bytes=%d data_len=%d imm=%d\n",*/
/*	   armci_me, msginfo->operation, bytes, msginfo->datalen,msginfo->tag.imm_msg);*/
/*    fflush(stdout);*/
    if(bytes<0 || bytes>IMM_BUF_LEN) {
      printf("%d(pid=%d): Trying to send too large a mesg. op=%d bytes=%d(max=%d) to=%d\n", armci_me, getpid(),msginfo->operation,bytes,IMM_BUF_LEN, proc);
      fflush(stdout);
      pause();
      assert(bytes>=0);
      assert(bytes <= IMM_BUF_LEN);
    }
    _armci_buf_ensure_pend_outstanding_op_per_node(buf,cluster);
/*     printf("%d: send_req. ensured pend os per node. to=%d op=%d\n", armci_me, msginfo->to,msginfo->operation); */
/*     fflush(stdout); */
#else
    _armci_buf_ensure_one_outstanding_op_per_node(buf,cluster);
#endif

    if(msginfo->operation == PUT || ARMCI_ACC(msginfo->operation)){
#if defined(PEND_BUFS)
      if(!msginfo->tag.imm_msg){
        msginfo->tag.data_ptr = (char *)(msginfo+1)+msginfo->dscrlen;
        msginfo->tag.data_len = msginfo->datalen;
      }
      else
	msginfo->tag.data_ptr = NULL;
#else
      {
          msginfo->tag.data_ptr = (void *)&msginfo->tag.ack;
      }
#endif
    }
    else {
       if(msginfo->operation == GET && !msginfo->bypass && msginfo->dscrlen
                       >= (msginfo->datalen-sizeof(int)))
         msginfo->tag.data_ptr = (char *)(msginfo+1)+msginfo->dscrlen;
       else
         msginfo->tag.data_ptr = GET_DATA_PTR(buf);
    }

    /*this has to be reset so that we can wait on it
      see ReadFromDirect*/
    msginfo->tag.ack = 0;
    msginfo->tag.ack_ptr = &(msginfo->tag.ack);

    if(DEBUG_CLN){
       printf("%d:the ack_ptr is initialised to %p, ack->value is %ld\n",
                 armci_me,msginfo->tag.ack_ptr,msginfo->tag.ack);fflush(stdout);
    }

    armci_init_vapibuf_send(snd_dscr, ssg_lst,buf, 
                            bytes, &client_memhandle);

/*    printf("%d: Sending req wr_id=%d to=%d\n",armci_me,snd_dscr->wr_id,proc);*/
/*    fflush(stdout);*/
    armci_vapi_post_send(1,cluster,snd_dscr,"send_req_msg:post_send");

    THREAD_UNLOCK(armci_user_threads.net_lock);

    if(DEBUG_CLN){
       printf("%d:client sent REQ=%d %d bytes serv=%d qp=%ld id =%ld lkey=%d\n",
               armci_me,msginfo->operation,bytes,cluster,
               (SRV_con+cluster)->qp,snd_dscr->wr_id,ssg_lst->lkey);
       fflush(stdout);
    }
    return(0);
}


/*\
 *  client waits for first phase ack before posting gather desr
\*/
void armci_wait_ack(char *buffer)
{
   long *flag;
   request_header_t *msginfo = (request_header_t *)(buffer);
   flag = (long*)&msginfo->tag.ack;

   while(armci_util_long_getval(flag) != ARMCI_STAMP);
   flag = 0;
}




void armci_client_direct_send(int p,void *src_buf, void *dst_buf, int len,void** contextptr,int nbtag,ARMCI_MEMHDL_T *lochdl,ARMCI_MEMHDL_T *remhdl)
{
sr_descr_t *dirdscr;
int clus = armci_clus_id(p);

    check_state_of_ib_connection(p, 0);

    THREAD_LOCK(armci_user_threads.net_lock);

    /*ID for the desr that comes from get_next_descr is already set*/
    dirdscr = armci_vapi_get_next_sdescr(nbtag,0);
    if(nbtag)*contextptr = dirdscr;

    armci_init_cbuf_srdma(&dirdscr->sdescr,dirdscr->sg_entry,src_buf,dst_buf,
                          len,lochdl,remhdl);

    armci_vapi_post_send(1,clus,&(dirdscr->sdescr),
                         "client_direct_send:post_send");

    /* the following unlock/lock ensures fairness (in case other threads are waiting
       on the lock) not required to work */
#if 1 
    THREAD_UNLOCK(armci_user_threads.net_lock);
    THREAD_LOCK(armci_user_threads.net_lock);
#endif

    if(nbtag==0)
       armci_send_complete(&(dirdscr->sdescr),"armci_client_direct_send",1);

    THREAD_UNLOCK(armci_user_threads.net_lock);
}

/*\ RDMA get
\*/
void armci_client_direct_get(int p, void *src_buf, void *dst_buf, int len,
                             void** cptr,int nbtag,ARMCI_MEMHDL_T *lochdl,
                             ARMCI_MEMHDL_T *remhdl)
{
int rc = 0;
sr_descr_t *dirdscr;
int clus = armci_clus_id(p);
    check_state_of_ib_connection(p, 0);
struct ibv_send_wr *bad_wr;

    THREAD_LOCK(armci_user_threads.net_lock);

    /*ID for the desr that comes from get_next_descr is already set*/
    dirdscr = armci_vapi_get_next_sdescr(nbtag,0);
    if(nbtag)*cptr = dirdscr;

    if(DEBUG_CLN){
      printf("\n%d: in direct get lkey=%d rkey=%d\n",armci_me,lochdl->lkey,
               remhdl->rkey);fflush(stdout);
    }

    armci_init_cbuf_rrdma(&dirdscr->sdescr,dirdscr->sg_entry,dst_buf,src_buf,
                          len,lochdl,remhdl);
    rc = ibv_post_send((SRV_con+clus)->qp, &(dirdscr->sdescr), &bad_wr);
    dassert1(1,rc==0,rc);

    /* unlock/lock to ensure fairness: allows others thread post before
       waiting for completion */
    /*VT?check to see if this should be UNLOCK followed by lock*/
#if 1
    THREAD_UNLOCK(armci_user_threads.net_lock);
    THREAD_LOCK(armci_user_threads.net_lock);
#endif

    if(!nbtag){
       armci_send_complete(&(dirdscr->sdescr),"armci_client_direct_get",1);
    }

    THREAD_UNLOCK(armci_user_threads.net_lock);
}

#define WQE_LIST_LENGTH 32
#define WQE_LIST_COUNT  1

/** Direct put into remote processor memory. Assumes that (and invoked
 *  only when) the source buffers in user memory are pinned as well.
 * @param operation PUT/GET
 * @param src_ptr Source pointer for data
 * @param src_stride_arr Strides on the source array
 * @param dst_ptr Destination pointer to start writing to
 * @param seq_count[stride_levels+1] #els in each stride
 * level. seg_count[0] is contiguous bytes
 * @param proc Destimation process
 * @param cptr OUT Pointer to store the descriptor to wait on for completion
 * @param nbtag IN Non-blocking tag (non-blocking op if nbtag!=0)
 * @param lochdl IN Local memory handle/key (registered memory stuff)
 * @param remhdl IN Remote memory handle/key
 * 
 */
#if 0
void armci_client_direct_rdma_strided(int operation, int proc,
				      char *src_ptr, int src_stride_arr[],
				      char *dst_ptr, int dst_stride_arr[],
				      int seg_count[],
				      int stride_levels,
				      void **cptr, int nbtag,
				      ARMCI_MEMHDL_T *lochdl,
				      ARMCI_MEMHDL_T *remhdl) {
  
  int rc;
  sr_descr_t *dirdscr;
  const int clus = armci_clus_id(proc);
  struct ibv_send_wr *bad_wr;
  struct ibv_send_wr sdscr[WQE_LIST_COUNT][WQE_LIST_LENGTH];
  struct ibv_sge     sg_entry[WQE_LIST_COUNT][WQE_LIST_LENGTH];
  int busy[WQE_LIST_COUNT], wait_count[WQE_LIST_COUNT],clst;
  int i, j, c, numposts;
  int idx[MAX_STRIDE_LEVEL];

  THREAD_LOCK(armci_user_threads.net_lock);

  assert(stride_levels >= 0);
  assert(stride_levels<=MAX_STRIDE_LEVEL);
  /*ID for the desr that comes from get_next_descr is already set*/
  dirdscr = armci_vapi_get_next_sdescr(nbtag,0);
  if(nbtag)*cptr = dirdscr;
  assert(dirdscr->tag == nbtag);

  if(DEBUG_CLN) {
    printf("\n%d: in direct rdma strided id=%d lkey=%ld rkey=%ld\n",
	   armci_me,dirdscr->sdescr.wr_id,lochdl->lkey,remhdl->rkey);fflush(stdout);
  }

  for(c=0; c<WQE_LIST_COUNT; c++) {
    busy[c]=0;
  }
  /*initialize fixed values for descriptors*/
  bzero(sdscr, WQE_LIST_COUNT*WQE_LIST_LENGTH*sizeof(struct ibv_send_wr));
  bzero(sg_entry, WQE_LIST_COUNT*WQE_LIST_LENGTH*sizeof(struct ibv_sge));
  for(j=0; j<WQE_LIST_COUNT; j++) {
    for(i=0; i<WQE_LIST_LENGTH; i++) {
      if(operation == PUT) 
	armci_init_cbuf_srdma(&sdscr[j][i],&sg_entry[j][i],NULL,NULL,seg_count[0],lochdl,remhdl);
      else if(operation == GET) 
	armci_init_cbuf_rrdma(&sdscr[j][i],&sg_entry[j][i],NULL,NULL,seg_count[0],lochdl,remhdl);
      else
	armci_die("rdma_strided: unsupported operation",operation);
      sdscr[j][i].wr_id = dirdscr->sdescr.wr_id;
      sdscr[j][i].send_flags = 0; /*non-signalled*/
      if(i<WQE_LIST_LENGTH-1)
	sdscr[j][i].next = &sdscr[j][i+1];
    }
  }
  /*post requests in a loop*/
  numposts=1;
  for(i=1; i<=stride_levels; i++) {
    numposts *= seg_count[i];
  }
/*   printf("%d: client rdma op=%d numposts=%d\n",armci_me,operation,numposts); */
  
  dirdscr->numofsends=0;
  bzero(idx, stride_levels*sizeof(int));
  int count = (numposts%WQE_LIST_LENGTH) ? (numposts%WQE_LIST_LENGTH):WQE_LIST_LENGTH;
  assert(count == ARMCI_MIN(count, numposts));
  clst=0;
  for(i=0; i<numposts; ) {
    for(j=i; j<i+count; j++) {
      int src_offset=0, dst_offset=0;
      for(c=0; c<stride_levels; c++) {
	src_offset += idx[c]*src_stride_arr[c];
	dst_offset += idx[c]*dst_stride_arr[c];
      }

/*       armci_client_direct_send(proc,src_ptr+src_offset,  */
/* 			       dst_ptr+dst_offset, seg_count[0], */
/* 			       NULL,0,lochdl,remhdl); */
      if(busy[clst]) {
	assert(wait_count[clst]>0);
	armci_send_complete(&dirdscr->sdescr,"client_direct_rdma_strided",wait_count[clst]);
	dirdscr->numofsends -= wait_count[clst];
	busy[clst]=0;
	wait_count[clst]=0;
      }

      if(operation == PUT) {
	sg_entry[clst][j-i].addr        = (uint64_t)(src_ptr + src_offset);
	sdscr[clst][j-i].wr.rdma.remote_addr = (uint64_t)(dst_ptr + dst_offset);
      }
      else if (operation == GET) {
	sg_entry[clst][j-i].addr        = (uint64_t)(dst_ptr + dst_offset);
	sdscr[clst][j-i].wr.rdma.remote_addr = (uint64_t)(src_ptr + src_offset);
      }
      assert(sg_entry[clst][j-i].length == seg_count[0]);

      idx[0] += 1;
      for(c=0;c<stride_levels-1 && idx[c]==seg_count[c+1]; c++) {
	idx[c]=0; idx[c+1]++;
      }
    }
    sdscr[clst][count-1].next=NULL;
    sdscr[clst][count-1].send_flags=IBV_SEND_SIGNALED; /*only the last one*/
    for(c=0; c<count-1; c++) {
      assert(sdscr[clst][c].next == &sdscr[clst][c+1]);
    }
    rc = ibv_post_send(SRV_con[clus].qp, sdscr[clst], &bad_wr);
    dassert1(1,rc==0,rc);
    dirdscr->numofsends += 1;
    wait_count[clst] = 1;
/*     armci_send_complete(&dirdscr->sdescr,"armci_client_direct_rdma_strided",count); */

    if(count < WQE_LIST_LENGTH) {
      sdscr[clst][count-1].next=&sdscr[clst][count]; /*reset it*/ 
    }
    sdscr[clst][count-1].send_flags=0; /*reset it*/
    i += count;
    count = ARMCI_MIN(WQE_LIST_LENGTH,numposts-i);
    assert(count==0 || count==WQE_LIST_LENGTH);
    clst = (clst+1)%WQE_LIST_COUNT;
  }

  if(!nbtag) {
    armci_send_complete(&dirdscr->sdescr,"armci_client_direct_get",dirdscr->numofsends);
    dirdscr->numofsends = 0;
    dirdscr->tag = 0;
  }  
  THREAD_UNLOCK(armci_user_threads.net_lock);
}
#else
void armci_client_direct_rdma_strided(int operation, int proc,
				      char *src_ptr, int src_stride_arr[],
				      char *dst_ptr, int dst_stride_arr[],
				      int seg_count[],
				      int stride_levels,
				      void **cptr, int nbtag,
				      ARMCI_MEMHDL_T *lochdl,
				      ARMCI_MEMHDL_T *remhdl) {
  int rc, i, j, c, busy[WQE_LIST_COUNT], clst, ctr;
  sr_descr_t *dirdscr;
  const int clus = armci_clus_id(proc);
  struct ibv_send_wr *bad_wr;
  struct ibv_send_wr sdscr[WQE_LIST_COUNT][WQE_LIST_LENGTH];
  struct ibv_sge     sg_entry[WQE_LIST_COUNT][WQE_LIST_LENGTH];
  stride_info_t sinfo, dinfo;

  THREAD_LOCK(armci_user_threads.net_lock);

  assert(stride_levels >= 0);
  assert(stride_levels<=MAX_STRIDE_LEVEL);
  /*ID for the desr that comes from get_next_descr is already set*/
  dirdscr = armci_vapi_get_next_sdescr(nbtag,0);
  if(nbtag)*cptr = dirdscr;
  assert(dirdscr->tag == nbtag);

  if(DEBUG_CLN) {
    printf("\n%d: in direct rdma strided id=%d lkey=%ld rkey=%ld\n",
	   armci_me,dirdscr->sdescr.wr_id,lochdl->lkey,remhdl->rkey);fflush(stdout);
  }

  /*initialize fixed values for descriptors*/
  bzero(sdscr, WQE_LIST_COUNT*WQE_LIST_LENGTH*sizeof(struct ibv_send_wr));
  bzero(sg_entry, WQE_LIST_COUNT*WQE_LIST_LENGTH*sizeof(struct ibv_sge));
  for(j=0; j<WQE_LIST_COUNT; j++) {
    for(i=0; i<WQE_LIST_LENGTH; i++) {
      if(operation == PUT) 
	armci_init_cbuf_srdma(&sdscr[j][i],&sg_entry[j][i],NULL,NULL,seg_count[0],lochdl,remhdl);
      else if(operation == GET) 
	armci_init_cbuf_rrdma(&sdscr[j][i],&sg_entry[j][i],NULL,NULL,seg_count[0],lochdl,remhdl);
      else
	armci_die("rdma_strided: unsupported operation",operation);
      sdscr[j][i].wr_id = dirdscr->sdescr.wr_id;
      sdscr[j][i].send_flags = 0; /*non-signalled*/
      if(i<WQE_LIST_LENGTH-1)
	sdscr[j][i].next = &sdscr[j][i+1];
    }
  }

  /*post requests in a loop*/
  armci_stride_info_init(&sinfo,src_ptr,stride_levels,src_stride_arr,seg_count);
  armci_stride_info_init(&dinfo,dst_ptr,stride_levels,dst_stride_arr,seg_count);
  assert(armci_stride_info_size(&sinfo)==armci_stride_info_size(&dinfo));
  
  dirdscr->numofsends=0;
  clst=ctr=0;
  bzero(busy, sizeof(int)*WQE_LIST_COUNT);
  while(armci_stride_info_has_more(&sinfo)) {
    assert(armci_stride_info_has_more(&dinfo));
    uint64_t saddr = (uint64_t)armci_stride_info_seg_ptr(&sinfo);
    uint64_t daddr = (uint64_t)armci_stride_info_seg_ptr(&dinfo);
    if(operation == PUT) {
      sg_entry[clst][ctr].addr = saddr;
      sdscr[clst][ctr].wr.rdma.remote_addr = daddr;
    }
    else if (operation == GET) {
      sg_entry[clst][ctr].addr = daddr;
      sdscr[clst][ctr].wr.rdma.remote_addr = saddr;
    }
    assert(sg_entry[clst][ctr].length == seg_count[0]);

    ctr+=1;
    armci_stride_info_next(&sinfo);
    armci_stride_info_next(&dinfo);
    if(ctr == WQE_LIST_LENGTH || !armci_stride_info_has_more(&sinfo)) {
      sdscr[clst][ctr-1].next=NULL;
      sdscr[clst][ctr-1].send_flags=IBV_SEND_SIGNALED; /*only the last one*/
      for(c=0; c<ctr-1; c++) {
	assert(sdscr[clst][c].next == &sdscr[clst][c+1]);
      }

      check_state_of_ib_connection(armci_clus_info[clus].master, 0);
      rc = ibv_post_send(SRV_con[clus].qp, sdscr[clst], &bad_wr);
      dassert1(1,rc==0,rc);
      busy[clst] = 1;
      dirdscr->numofsends += 1;
      if(ctr<WQE_LIST_LENGTH) 
	sdscr[clst][ctr-1].next = &sdscr[clst][ctr];
      sdscr[clst][ctr-1].send_flags = 0;

      ctr=0;
      clst = (clst+1)%WQE_LIST_COUNT;
      if(busy[clst]) {
	armci_send_complete(&dirdscr->sdescr,"client_direct_rdma_strided",1);
	dirdscr->numofsends -= 1;
	busy[clst]=0;	
      }
    }
  }
  armci_stride_info_destroy(&sinfo);
  armci_stride_info_destroy(&dinfo);

  if(!nbtag) {
    armci_send_complete(&dirdscr->sdescr,"armci_client_direct_get",dirdscr->numofsends);
    dirdscr->numofsends = 0;
    dirdscr->tag = 0;
  }  
  THREAD_UNLOCK(armci_user_threads.net_lock);
}
#endif

#if defined(PEND_BUFS)
int armci_server_msginfo_to_pbuf_index(request_header_t *msginfo) {
  int index=-1, i;
  vapibuf_pend_t *pbuf=NULL;
  
  assert(!msginfo->tag.imm_msg);
  for(i = 0; i<PENDING_BUF_NUM; i++) {
    if(serv_pendbuf_arr[i].buf == (char *)msginfo) {
      pbuf = &serv_pendbuf_arr[i];
      index = i;
      break;
    }
  }
  return index;
}


/** Routine for server to RDMA strided data to the client-side buffers
 * (allocated through buffers.c). This is to be used instead of
 * copying the data to immediate or pending buffers when possible.
 */
#if 0
void armci_server_rdma_strided_to_contig(char *src_ptr, int src_stride_arr[],
					 int seg_count[],
					 int stride_levels,
					 char *dst_ptr, int proc,
					 request_header_t *msginfo) {
  int rc, i, j, c, busy[WQE_LIST_COUNT], clst, ctr, wr_id;
  sr_descr_t *dirdscr;
  struct ibv_send_wr *bad_wr, sdscr1;
  struct ibv_send_wr sdscr[WQE_LIST_COUNT][WQE_LIST_LENGTH];
  struct ibv_sge     sg_entry[WQE_LIST_COUNT][WQE_LIST_LENGTH];
  stride_info_t sinfo;
  uint64_t daddr;
  ARMCI_MEMHDL_T *loc_memhdl;
  ARMCI_MEMHDL_T *rem_memhdl = &handle_array[proc];
  
  THREAD_LOCK(armci_user_threads.net_lock);

  assert(msginfo->operation == GET);
  assert(stride_levels >= 0);
  assert(stride_levels<=MAX_STRIDE_LEVEL);

  if(!get_armci_region_local_hndl(src_ptr,armci_clus_id(armci_me), &loc_memhdl)) {
    armci_die("rdma_strided_to_contig: failed to get local handle\n",0);
  }

  if(!msginfo->tag.imm_msg) {
    int index = armci_server_msginfo_to_pbuf_index(msginfo);
    assert(index>=0);
    wr_id = PBUF_BUFID_TO_PUT_WRID(index);
  }
  else {
    wr_id = DSCRID_IMMBUF_RESP_END-1-proc;
  }
  bzero(&sdscr1, sizeof(sdscr1));
  sdscr1.wr_id = wr_id;

  if(DEBUG_CLN) {
    printf("\n%d: in rdma strided to contig id=%d lkey=%ld rkey=%ld\n",
	   armci_me,wr_id,loc_memhdl->lkey,rem_memhdl->rkey);
    fflush(stdout);
  }

  /*initialize fixed values for descriptors*/
  bzero(sdscr, WQE_LIST_COUNT*WQE_LIST_LENGTH*sizeof(struct ibv_send_wr));
  bzero(sg_entry, WQE_LIST_COUNT*WQE_LIST_LENGTH*sizeof(struct ibv_sge));
  for(j=0; j<WQE_LIST_COUNT; j++) {
    for(i=0; i<WQE_LIST_LENGTH; i++) {
      armci_init_cbuf_srdma(&sdscr[j][i],&sg_entry[j][i],NULL,NULL,seg_count[0],loc_memhdl,rem_memhdl);
      sdscr[j][i].wr_id = wr_id;
      sdscr[j][i].send_flags = 0; /*non-signalled*/
/*       sdscr[j][i].send_flags = IBV_SEND_SIGNALED; /\*signalled*\/ */
      if(i<WQE_LIST_LENGTH-1)
	sdscr[j][i].next = &sdscr[j][i+1];
    }
  }

  /*post requests in a loop*/
  sinfo = armci_stride_info_init(src_ptr,stride_levels,src_stride_arr,seg_count);
  
  clst=ctr=0;
  bzero(busy, sizeof(int)*WQE_LIST_COUNT);
  daddr = (uint64_t)dst_ptr;
  while(armci_stride_info_has_more(sinfo)) {
    uint64_t saddr = (uint64_t)armci_stride_info_seg_ptr(sinfo);
    sg_entry[clst][ctr].addr = saddr;
    sdscr[clst][ctr].wr.rdma.remote_addr = daddr;
    assert(sg_entry[clst][ctr].length == seg_count[0]);

    ctr+=1;
    daddr += seg_count[0];
    armci_stride_info_next(sinfo);
    if(ctr == WQE_LIST_LENGTH || !armci_stride_info_has_more(sinfo)) {
      sdscr[clst][ctr-1].next=NULL;
      if(!armci_stride_info_has_more(sinfo)) {
	sdscr[clst][ctr-1].send_flags=IBV_SEND_SIGNALED; /*only the last one*/
      }
      for(c=0; c<ctr-1; c++) {
	assert(sdscr[clst][c].next == &sdscr[clst][c+1]);
      }
      rc = ibv_post_send(CLN_con[proc].qp, sdscr[clst], &bad_wr);
      dassert1(1,rc==0,rc);
      busy[clst] = 1;
#if 0
      armci_send_complete(&sdscr1,"serv_rdma_to_contig",ctr);
      busy[clst] = 0;
#endif
      if(ctr<WQE_LIST_LENGTH) 
	sdscr[clst][ctr-1].next = &sdscr[clst][ctr];
      sdscr[clst][ctr-1].send_flags = 0;

      ctr=0;
      clst = (clst+1)%WQE_LIST_COUNT;
#if 0
      if(busy[clst]) {
	armci_send_complete(&sdscr1,"client_direct_rdma_strided",1);
	busy[clst]=0;	
      }
#endif
    }
  }
  armci_stride_info_destroy(&sinfo);
  assert(proc == msginfo->from);
  THREAD_UNLOCK(armci_user_threads.net_lock);
}
#else

#define MAX_NUM_SGE 64

/*same as above, but uses gather rdma writes*/
void armci_server_rdma_strided_to_contig(char *src_ptr, int src_stride_arr[],
					 int seg_count[],
					 int stride_levels,
					 char *dst_ptr, int proc,
					 request_header_t *msginfo) {
  int rc, ctr, wr_id, bytes;
  struct ibv_send_wr *bad_wr, sdscr1, sdscr;
  struct ibv_sge     sg_entry[MAX_NUM_SGE];
  stride_info_t sinfo;
  uint64_t daddr;
  ARMCI_MEMHDL_T *loc_memhdl;
  ARMCI_MEMHDL_T *rem_memhdl = &handle_array[proc];
  const int max_num_sge = ARMCI_MIN(MAX_NUM_SGE, armci_max_num_sg_ent);
  int numposts=0, numsegs=0;
  
  THREAD_LOCK(armci_user_threads.net_lock);

  assert(msginfo->operation == GET);
  assert(stride_levels >= 0);
  assert(stride_levels<=MAX_STRIDE_LEVEL);

  if(!get_armci_region_local_hndl(src_ptr,armci_clus_id(armci_me), &loc_memhdl)) {
    armci_die("rdma_strided_to_contig: failed to get local handle\n",0);
  }

  if(!msginfo->tag.imm_msg) {
    int index = armci_server_msginfo_to_pbuf_index(msginfo);
    assert(index>=0);
    wr_id = PBUF_BUFID_TO_PUT_WRID(index);
  }
  else {
    wr_id = DSCRID_IMMBUF_RESP_END-1-proc;
  }
  bzero(&sdscr1, sizeof(sdscr1));
  sdscr1.wr_id = wr_id;

  if(DEBUG_CLN) {
    printf("\n%d: in rdma strided to contig id=%d lkey=%ld rkey=%ld\n",
	   armci_me,wr_id,loc_memhdl->lkey,rem_memhdl->rkey);
    fflush(stdout);
  }

  /*initialize fixed values for descriptors*/
  bzero(&sdscr, sizeof(sdscr));
  bzero(sg_entry, max_num_sge*sizeof(struct ibv_sge));
  armci_init_cbuf_srdma(&sdscr,&sg_entry[0],NULL,NULL,seg_count[0],loc_memhdl,rem_memhdl);
  sdscr.send_flags = 0; /*non-signalled*/
  sdscr.num_sge    = 0; /*set below in the loop*/
  sdscr.wr_id      = wr_id;

  for(ctr=0; ctr<max_num_sge; ctr++) {
    sg_entry[ctr].length = seg_count[0];
    sg_entry[ctr].lkey = loc_memhdl->lkey;
  }
  
  /*post requests in a loop*/
  armci_stride_info_init(&sinfo,src_ptr,stride_levels,src_stride_arr,seg_count);
  
  numposts = numsegs = 0;
  ctr=0;
  daddr = (uint64_t)dst_ptr;
  bytes=0;
  while(armci_stride_info_has_more(&sinfo)) {
    sg_entry[ctr].addr = (uint64_t)armci_stride_info_seg_ptr(&sinfo);
    assert(sg_entry[ctr].length == seg_count[0]);

    sdscr.num_sge += 1;
    bytes += seg_count[0];
    ctr+=1;
    numsegs += 1;
    armci_stride_info_next(&sinfo);
    if(ctr == max_num_sge || !armci_stride_info_has_more(&sinfo)) {
      sdscr.wr.rdma.remote_addr = daddr;
      if(!armci_stride_info_has_more(&sinfo)) {
	sdscr.send_flags=IBV_SEND_SIGNALED; /*only the last one*/
      }
      else {
	assert(sdscr.send_flags == 0);
      }
      assert(ctr == sdscr.num_sge);
      rc = ibv_post_send(CLN_con[proc].qp, &sdscr, &bad_wr);
      dassert1(1,rc==0,rc);

      numposts += 1;
      ctr=0;
      sdscr.num_sge = 0;
      daddr += bytes;
      bytes = 0;
    }
  }
/*   printf("%d(s): scatgat write numposts=%d numsegs=%d\n",armci_me,numposts,numsegs); */
  armci_stride_info_destroy(&sinfo);
  assert(proc == msginfo->from);
  THREAD_UNLOCK(armci_user_threads.net_lock);
}

/*Directly read data from client buffers into remote memory. Data is
  contiguous in client-side. */
void armci_server_rdma_contig_to_strided(char *src_ptr, int proc,
					 char *dst_ptr, 
					 int dst_stride_arr[],
					 int seg_count[],
					 int stride_levels,
					 request_header_t *msginfo) {
  int rc, ctr, wr_id, bytes;
  struct ibv_send_wr *bad_wr, sdscr1, sdscr;
  struct ibv_sge     sg_entry[MAX_NUM_SGE];
  stride_info_t dinfo;
  uint64_t saddr;
  ARMCI_MEMHDL_T *loc_memhdl;
  ARMCI_MEMHDL_T *rem_memhdl = &handle_array[proc];
  const int max_num_sge = ARMCI_MIN(MAX_NUM_SGE, armci_max_num_sg_ent);
  int numposts=0, numsegs=0;
  
  THREAD_LOCK(armci_user_threads.net_lock);

  assert(msginfo->operation == PUT);
  assert(stride_levels >= 0);
  assert(stride_levels<=MAX_STRIDE_LEVEL);

  if(!get_armci_region_local_hndl(dst_ptr,armci_clus_id(armci_me), &loc_memhdl)) {
    armci_die("rdma_strided_to_contig: failed to get local handle\n",0);
  }

  if(!msginfo->tag.imm_msg) {
    int index = armci_server_msginfo_to_pbuf_index(msginfo);
    assert(index>=0);
    wr_id = PBUF_BUFID_TO_GET_WRID(index);
  }
  else {
    wr_id = DSCRID_IMMBUF_RESP_END-1-proc;
  }
  bzero(&sdscr1, sizeof(sdscr1));
  sdscr1.wr_id = wr_id;

  if(DEBUG_CLN) {
    printf("\n%d: in rdma strided to contig id=%d lkey=%ld rkey=%ld\n",
	   armci_me,wr_id,loc_memhdl->lkey,rem_memhdl->rkey);
    fflush(stdout);
  }

  /*initialize fixed values for descriptors*/
  bzero(&sdscr, sizeof(sdscr));
  bzero(sg_entry, max_num_sge*sizeof(struct ibv_sge));
  armci_init_cbuf_rrdma(&sdscr,&sg_entry[0],NULL,NULL,seg_count[0],loc_memhdl,rem_memhdl);
  sdscr.send_flags = 0; /*non-signalled*/
  sdscr.num_sge    = 0; /*set below in the loop*/
  sdscr.wr_id      = wr_id;

  for(ctr=0; ctr<max_num_sge; ctr++) {
    sg_entry[ctr].length = seg_count[0];
    sg_entry[ctr].lkey = loc_memhdl->lkey;
  }
  
  /*post requests in a loop*/
  armci_stride_info_init(&dinfo,dst_ptr,stride_levels,dst_stride_arr,seg_count);
  
  numposts = numsegs = 0;
  ctr=0;
  saddr = (uint64_t)src_ptr;
  bytes=0;
  while(armci_stride_info_has_more(&dinfo)) {
    sg_entry[ctr].addr = (uint64_t)armci_stride_info_seg_ptr(&dinfo);
    assert(sg_entry[ctr].length == seg_count[0]);

    sdscr.num_sge += 1;
    bytes += seg_count[0];
    ctr+=1;
    numsegs += 1;
    armci_stride_info_next(&dinfo);
    if(ctr == max_num_sge || !armci_stride_info_has_more(&dinfo)) {
      sdscr.wr.rdma.remote_addr = saddr;
      if(!armci_stride_info_has_more(&dinfo)) {
	sdscr.send_flags=IBV_SEND_SIGNALED; /*only the last one*/
      }
      else {
	assert(sdscr.send_flags == 0);
      }
      assert(ctr == sdscr.num_sge);
      rc = ibv_post_send(CLN_con[proc].qp, &sdscr, &bad_wr);
      dassert1(1,rc==0,rc);

      numposts += 1;
      ctr=0;
      sdscr.num_sge = 0;
      saddr += bytes;
      bytes = 0;
    }
  }
/*   printf("%d(s): scatgat write numposts=%d numsegs=%d\n",armci_me,numposts,numsegs); */
  armci_stride_info_destroy(&dinfo);
  assert(proc == msginfo->from);
  THREAD_UNLOCK(armci_user_threads.net_lock);
}

#endif
#endif

char *armci_ReadFromDirect(int proc, request_header_t *msginfo, int len)
{
int cluster = armci_clus_id(proc);
vapibuf_ext_t* ecbuf=BUF_TO_ECBUF(msginfo);
char *dataptr = GET_DATA_PTR(ecbuf->buf);
extern void armci_util_wait_int(volatile int *,int,int);

    if(DEBUG_CLN){ printf("%d(c):read direct %d qp=%p\n",armci_me,
                len,&(SRV_con+cluster)->qp); fflush(stdout);
    }

    if(mark_buf_send_complete[ecbuf->snd_dscr.wr_id]==0)
       armci_send_complete(&(ecbuf->snd_dscr),"armci_ReadFromDirect",1); 

    if(!msginfo->bypass){
       long *flag;
       int *last;
       int loop = 0;
       flag = &(msginfo->tag.ack);
       if(msginfo->operation==GET){
         last = (int *)(dataptr+len-sizeof(int));
         if(msginfo->dscrlen >= (len-sizeof(int))){
           last = (int *)(dataptr+len+msginfo->dscrlen-sizeof(int));
           dataptr+=msginfo->dscrlen;
         }

         if(DEBUG_CLN){
           printf("\n%d:flagval=%d at ptr=%p ack=%ld dist=%d\n",armci_me,*last,
                   last,*flag,len);fflush(stdout);
         }

         while(armci_util_int_getval(last) == ARMCI_STAMP &&
               armci_util_long_getval(flag)  != ARMCI_STAMP){
           loop++;
           loop %=100000;
           if(loop==0){
             if(DEBUG_CLN){
               printf("%d: client last(%p)=%d flag(%p)=%ld off=%d\n",
                      armci_me,last,*last,flag,*flag,msginfo->datalen);
               fflush(stdout);
             }
           }
         }
         *flag = 0L;
       }
       else if(msginfo->operation == REGISTER){
         while(armci_util_long_getval(flag)  != ARMCI_STAMP){
           loop++;
           loop %=100000;
           if(loop==0){
             if(DEBUG_CLN){
               printf("%d: client flag(%p)=%ld off=%d\n",
                      armci_me,flag,*flag,msginfo->datalen);
               fflush(stdout);
             }
           }
         }
       }
       else{
         int *flg = (int *)(dataptr+len);
         while(armci_util_int_getval(flg) != ARMCI_STAMP){
           loop++;
           loop %=100000;
           if(loop==0){
             if(DEBUG_CLN){
               printf("%d: client waiting (%p)=%d off=%d\n",
                      armci_me,flg,*flg,len);
               fflush(stdout);
             }
           }
         }
       }
    }
    return dataptr;
}


#ifdef GET_STRIDED_COPY_PIPELINED
/**Same as armci_ReadFromDirect, except reads partial segments
 *  (identify by stamping done in armci_send_req_msg() and
 *  returns. Note that the return value is the starting pointer of the
 *  buffer containig the data. It is the same for all the segments
 *  read for a message. 
 * @param proc IN Read data corresponding to an earlier req to this proc
 * @param msginfo IN The request for which we are reading now
 * @param len IN #bytes in the total response
 * @param bytes_done OUT @bytes of the total response read so far (monotonic)
 * @return Starting pointer to the buffer containing the data
 */
char *armci_ReadFromDirectSegment(int proc, request_header_t *msginfo, int len, int *bytes_done) {
  int cluster = armci_clus_id(proc);
  vapibuf_ext_t* ecbuf=BUF_TO_ECBUF(msginfo);
  char *dataptr = GET_DATA_PTR(ecbuf->buf);
  extern void armci_util_wait_int(volatile int *,int,int);

  if(DEBUG_CLN){ printf("%d(c):read direct %d qp=%p\n",armci_me,
			len,&(SRV_con+cluster)->qp); fflush(stdout);
  }

  if(mark_buf_send_complete[ecbuf->snd_dscr.wr_id]==0)
    armci_send_complete(&(ecbuf->snd_dscr),"armci_ReadFromDirect",1); 

  if(!msginfo->bypass){
    long *flag;
    int *last, *mid1, *mid2, third;
    int loop = 0;
    flag = &(msginfo->tag.ack);
    if(msginfo->operation==GET){
      last = (int *)(dataptr+len-sizeof(int));
      if(msginfo->dscrlen >= (len-sizeof(int))){
	last = (int *)(dataptr+len+msginfo->dscrlen-sizeof(int));
	dataptr+=msginfo->dscrlen;
      }
      third = (last-(int*)(msginfo->dscrlen+(char*)(msginfo+1)))/3;
      mid2 = (last - third);
      mid1 = mid2 - third;

      if(DEBUG_CLN){
	printf("\n%d:flagval=%d at ptr=%p ack=%ld dist=%d\n",armci_me,*last,
	       last,*flag,len);fflush(stdout);
      }

      while(armci_util_int_getval(last) == ARMCI_STAMP &&
	    armci_util_long_getval(flag)  != ARMCI_STAMP){
	loop++;
	loop %=100000;
	if(loop==0){
	  if(DEBUG_CLN){
	    printf("%d: client last(%p)=%d flag(%p)=%ld off=%d\n",
		   armci_me,last,*last,flag,*flag,msginfo->datalen);
	    fflush(stdout);
	  }
	}

	{
	  int ssize = GET_STRIDED_COPY_PIPELINED_SIZE/sizeof(int);
	  int *sfirst = (int*)(msginfo->dscrlen+(char*)(msginfo+1))+ssize; /*stamping
									     can start here*/
	  int *slast = last;
	  int off = (((int *)(dataptr+*bytes_done)-sfirst+ssize)/ssize)*ssize;
	  int *ptr = sfirst+off;
	  dassert(1,off>=0);
	  dassert(1,(void *)sfirst>dataptr);
	  dassert(1,(void *)ptr>dataptr);
	  if(ptr<=slast && armci_util_int_getval(ptr)!=ARMCI_STAMP) {
	    *bytes_done = ((char*)ptr)-dataptr;
	    return dataptr;
	  }
	}
      }
      *flag = 0L;
      *bytes_done = len;
      return dataptr;
    }
    else if(msginfo->operation == REGISTER){
      while(armci_util_long_getval(flag)  != ARMCI_STAMP){
	loop++;
	loop %=100000;
	if(loop==0){
	  if(DEBUG_CLN){
	    printf("%d: client flag(%p)=%ld off=%d\n",
		   armci_me,flag,*flag,msginfo->datalen);
	    fflush(stdout);
	  }
	}
      }
    }
    else{
      int *flg = (int *)(dataptr+len);
      while(armci_util_int_getval(flg) != ARMCI_STAMP){
	loop++;
	loop %=100000;
	if(loop==0){
	  if(DEBUG_CLN){
	    printf("%d: client waiting (%p)=%d off=%d\n",
		   armci_me,flg,*flg,len);
	    fflush(stdout);
	  }
	}
      }
    }
  }
  *bytes_done = len;
  return dataptr;
}
#endif

/**
  * @param proc IN id of remote client to put to
  * @param buf IN local buf (has to be registered)
 */
void armci_send_data_to_client(int proc, void *buf, int bytes,void *dbuf)
{
  int i, rc = 0;
    struct ibv_send_wr *bad_wr;
    struct ibv_send_wr sdscr;
    struct ibv_sge ssg_entry;

    if(DEBUG_SERVER){
       printf("\n%d(s):sending data to client %d at %p flag = %p bytes=%d\n",
               armci_me,
	      proc,dbuf,(char *)dbuf+bytes-sizeof(int),bytes);fflush(stdout);
    }

    memset(&sdscr,0,sizeof(struct ibv_send_wr));
    memset(&ssg_entry,0,sizeof(ssg_entry));
    armci_init_cbuf_srdma(&sdscr,&ssg_entry,buf,dbuf,bytes,
                          &serv_memhandle,(handle_array+proc));

    if(DEBUG_SERVER){
       printf("\n%d(s):handle_array[%d]=%p dbuf=%p flag=%p bytes=%d\n",armci_me,
              proc,&handle_array[proc],(char *)dbuf,
              (char *)dbuf+bytes-sizeof(int),bytes);
       fflush(stdout);
    }

#if defined(PEND_BUFS)
    for(i=proc*(IMM_BUF_NUM+1); i<(proc+1)*(IMM_BUF_NUM+1); i++) {
      if((char*)buf>= serv_buf_arr[i]->buf && 
	 (char*)buf<IMM_BUF_LEN+(char*)serv_buf_arr[i]->buf)
	break;
    }

#if SRI_CORRECT
     if(i<(proc+1)*(IMM_BUF_NUM+1)) {
      /*Message from an immediate buffer*/
     assert(serv_buf_arr[i]->send_pending==0);
      serv_buf_arr[i]->send_pending=1;
      sdscr.wr_id = DSCRID_IMMBUF_RESP+i;
    }
    else 
#endif
      {
	sdscr.wr_id = DSCRID_IMMBUF_RESP+armci_nproc*(IMM_BUF_NUM+1)+1;
      }
/* #endif */

/* #if defined(PEND_BUFS) */
/*     { */
/*       static uint64_t ctr=DSCRID_IMMBUF_RESP; */
/*       sdscr.wr_id = ctr; */
/*       ctr = (ctr+1-DSCRID_IMMBUF_RESP)%(DSCRID_IMMBUF_RESP_END-DSCRID_IMMBUF_RESP)+DSCRID_IMMBUF_RESP; */
/*     } */
#else
    sdscr.wr_id = proc+armci_nproc;
#endif
    rc = ibv_post_send((CLN_con+proc)->qp, &sdscr, &bad_wr);
    dassert1(1,rc==0,rc);

#if !defined(PEND_BUFS)
    armci_send_complete(&sdscr,"armci_send_data_to_client",1);
#endif
}

void armci_WriteToDirect(int proc, request_header_t* msginfo, void *buf)
{
int bytes;
int *last;
    ARMCI_PR_DBG("enter",0);
    bytes = (int)msginfo->datalen;
    if(DEBUG_SERVER){
      printf("%d(s):write to direct sent %d to %d at %p\n",armci_me,
             bytes,proc,(char *)msginfo->tag.data_ptr);
      fflush(stdout);
    }
    if(msginfo->operation!=GET){
       *(int *)((char *)buf+bytes)=ARMCI_STAMP;
       bytes+=sizeof(int);
    }
#if defined(PEND_BUFS)
    if(!msginfo->tag.imm_msg) {
      int i;
/*       fprintf(stderr, "%d:: Not immediate mesg operated on\n", armci_me); */
      assert(msginfo->operation == GET); /*nothing else uses this for now*/
      /**This is a pending buf*/
      vapibuf_pend_t *pbuf=NULL;
      int index;
      for(i = 0; i<PENDING_BUF_NUM; i++) {
	if(serv_pendbuf_arr[i].buf == (char *)msginfo) {
	  pbuf = &serv_pendbuf_arr[i];
	  index = i;
	  break;
	}
      }
      assert(pbuf != NULL);
      assert(sizeof(request_header_t)+msginfo->dscrlen+bytes<PENDING_BUF_LEN);
      _armci_send_data_to_client_pbuf(proc, &pbuf->sdscr,
				      PBUF_BUFID_TO_PUT_WRID(index),
				      &pbuf->sg_entry,
				      msginfo->tag.data_ptr, buf,
				      bytes);
    }
    else 
#endif
    {
      armci_send_data_to_client(proc,buf,bytes,msginfo->tag.data_ptr);
    }
    /*if(msginfo->dscrlen >= (bytes-sizeof(int)))
       last = (int*)(((char*)(buf)) + (msginfo->dscrlen+bytes - sizeof(int)));
    else*/
       last = (int*)(((char*)(buf)) + (bytes - sizeof(int)));

    if(msginfo->operation==GET && *last == ARMCI_STAMP){
       SERVER_SEND_ACK(msginfo->from);
    }
    armci_ack_proc=NONE;
    ARMCI_PR_DBG("exit",0);
}


#if defined(PEND_BUFS)
void armci_rcv_req(void *mesg,void *phdr,void *pdescr,void *pdata,int *buflen)
{
  request_header_t *msginfo = *(request_header_t**)mesg;
  *(void **)phdr = msginfo;

  if(msginfo->tag.imm_msg) 
    *buflen = IMM_BUF_LEN - sizeof(request_header_t) - msginfo->dscrlen;
  else 
    *buflen = PENDING_BUF_LEN - sizeof(request_header_t) - msginfo->dscrlen;
  
  *(void **)pdata = msginfo->dscrlen + (char *)(msginfo+1);
  if(msginfo->bytes)
    *(void **)pdescr = msginfo+1;
  else
    *(void **)pdescr = NULL;
}
#else
void armci_rcv_req(void *mesg,void *phdr,void *pdescr,void *pdata,int *buflen)
{
  vapibuf_t *cbuf = (vapibuf_t*)mesg;
  request_header_t *msginfo = (request_header_t *)cbuf->buf;
  *(void **)phdr = msginfo;

  ARMCI_PR_DBG("enter",msginfo->operation);
  if(DEBUG_SERVER){
    printf("%d(server): got %d req (dscrlen=%d datalen=%d) from %d\n",
	   armci_me, msginfo->operation, msginfo->dscrlen,
	   msginfo->datalen, msginfo->from); fflush(stdout);
  }

  /* we leave room for msginfo on the client side */
  *buflen = MSG_BUFLEN - sizeof(request_header_t);
  
  if(msginfo->bytes) {
    *(void **)pdescr = msginfo+1;
    if(msginfo->operation == GET)
      *(void **)pdata = MessageRcvBuffer;
    else
      *(void **)pdata = msginfo->dscrlen + (char*)(msginfo+1);
  }else {
    *(void**)pdescr = NULL;
    *(void**)pdata = MessageRcvBuffer;
  }
  ARMCI_PR_DBG("exit",msginfo->operation);
}
#endif

static void posts_scatter_desc(sr_descr_t *pend_dscr,int proc,int type)
{
int rc;
int cluster = armci_clus_id(proc);
struct ibv_recv_wr *scat_dscr;
struct ibv_recv_wr *bad_wr;

    scat_dscr = &pend_dscr->rdescr;

    /*armci_vapi_print_dscr_info(NULL,scat_dscr);*/
    if((type==SERV && DEBUG_SERVER) || (type==CLN && DEBUG_CLN)){
       printf("%d(%d) : inside posts scatter dscr, id is %d\n",
	      armci_me,type,scat_dscr->wr_id);
       fflush(stdout);
    }

    if(type == SERV)
        rc = ibv_post_recv((CLN_con + proc)->qp, scat_dscr, &bad_wr);
    else
        rc = ibv_post_recv((SRV_con+cluster)->qp, scat_dscr, &bad_wr);
    dassert1(1,rc==0,rc);

    if((type==SERV && DEBUG_SERVER) || (type==CLN && DEBUG_CLN) ) {
       printf("\n%d: list_length is %d, id is %ld\n",
	      armci_me,scat_dscr->num_sge,scat_dscr->wr_id);
       fflush(stdout);
    }
}


/*\
 *  client calls from request.c
 *  server calls from ds-shared.c
\*/
static sr_descr_t serv_blocking_scatter_dscr;
static sr_descr_t client_blocking_scatter_dscr;
void armci_post_scatter(void *dest_ptr, int dest_stride_arr[], int count[],
     int stride_levels, armci_vapi_memhndl_t *mhandle,
     int proc, int nbtag, int type, sr_descr_t **srd)
{
    int i;
    int total_size = 0;
    int total_of_2D = 1;
    int index[MAX_STRIDE_LEVEL], unit[MAX_STRIDE_LEVEL];
    int j,k,y,z;
    int num_dscr = 0;
    int num_xmit = 0, num_seg, max_seg, rem_seg,vecind;
    char* src, *src1;
    sr_descr_t *pend_dscr;
    struct ibv_sge *scat_sglist;
    struct ibv_recv_wr *scat_dscr;

    if((type==SERV && DEBUG_SERVER) || (type==CLN && DEBUG_CLN) ){
       printf("%d(%d)  : inside post_scatter %d\n",armci_me,type,nbtag);
       fflush(stdout);
    }

    max_seg =  armci_max_num_sg_ent;

    THREAD_LOCK(armci_user_threads.net_lock);

    if(nbtag){
       pend_dscr = armci_vapi_get_next_rdescr(nbtag,1);
       if(srd!=NULL)*srd=pend_dscr;
    }
    else{
       pend_dscr = &client_blocking_scatter_dscr;
       pend_dscr->rdescr.wr_id=DSCRID_SCATGAT + MAX_PENDING;
    }

    /*pend_dscr->proc = proc;*/
    pend_dscr->numofrecvs=0;

    scat_dscr = &pend_dscr->rdescr;
    scat_sglist = pend_dscr->sg_entry;
    /* scat_dscr->opcode = VAPI_RECEIVE; no ->opcode in ibv_recv_wr */
    /* scat_dscr->comp_type = VAPI_SIGNALED; no ->comp_type in ibv_recv_wr */
    scat_dscr->sg_list = scat_sglist;
    scat_dscr->num_sge = 0;

    index[2] = 0; unit[2] = 1;
    if(stride_levels > 1){
       total_of_2D = count[2];
       for(j=3; j<=stride_levels; j++){
         index[j] = 0; unit[j] = unit[j-1]*count[j-1];
         total_of_2D*=count[j];
       }
    }

    num_xmit = total_of_2D*count[1]/max_seg;
    rem_seg = (total_of_2D*count[1])%max_seg;
    if(num_xmit == 0) num_xmit = 1;
    else if(rem_seg!= 0)num_xmit++;


    if ((type==SERV && DEBUG_SERVER) || (type==CLN && DEBUG_CLN) ) {
       printf("%d(%d):armci_post_scatter num_xmit = %d\t, rem_seg = %d\n",
               armci_me,type,num_xmit,rem_seg);
       fflush(stdout);
    }

    k=0; vecind = 0;
    if(rem_seg!=0 && k==(num_xmit-1))num_seg = rem_seg;
    else num_seg = max_seg;

    y=0,z=0;
    for(i=0;i<total_of_2D;i++){
       src = (char *)dest_ptr;
       for(j=2;j<=stride_levels;j++){
         src+= index[j]*dest_stride_arr[j-1];
         if(((i+1)%unit[j]) == 0) index[j]++;
         if(index[j] >= count[j]) index[j] =0;
       }
       src1 = src;

       for(j=0; j<count[1]; j++, vecind++){
         if(vecind == num_seg) {
           posts_scatter_desc(pend_dscr,proc,type);
           pend_dscr->numofrecvs++;

           /* the previous one has been posted, start off new*/
           scat_dscr->num_sge = 0;
           y = 0; /* reuse the same scatter descriptor */
           vecind=0;total_size=0;k++;
           if(rem_seg!=0 && k==(num_xmit-1))num_seg = rem_seg;
         }
         /* fill the scatter descriptor */
         scat_sglist[y].addr = (uint64_t)src1;
         scat_sglist[y].lkey = mhandle->lkey;
         scat_sglist[y].length = count[0];
         scat_dscr->num_sge++;
         src1 += dest_stride_arr[0];
         y++;

       }

       if(vecind == num_seg){
         posts_scatter_desc(pend_dscr,proc,type);
         pend_dscr->numofrecvs++;

         /* the previous one has been posted, start off new*/
         scat_dscr->num_sge = 0;
         y =0 ;
         vecind = 0; total_size=0; k++;
         if(rem_seg!=0 && k==(num_xmit-1))num_seg=rem_seg;
         else num_seg = max_seg;
       }

    }

    THREAD_UNLOCK(armci_user_threads.net_lock);

/*     printf("%d(s): num scatters posted=%d\n", armci_me,pend_dscr->numofrecvs); */
    if(!nbtag){
       /*if blocking call wait_for_blocking_scatter to complete*/
    }
    return;
}

void armci_wait_for_blocking_scatter()
{
sr_descr_t *pend_dscr=&client_blocking_scatter_dscr;
int i;
    armci_recv_complete(&pend_dscr->rdescr,"armci_post_scatter",pend_dscr->numofrecvs);
}


/*\
 *  function used by armci_post_gather to actually post the sctter list
\*/
static void posts_gather_desc(sr_descr_t *pend_dscr,int proc,int type)
{
    int rc;
    int cluster = armci_clus_id(proc);
    struct ibv_send_wr *gat_dscr;
    struct ibv_send_wr *bad_wr;

    THREAD_LOCK(armci_user_threads.net_lock);

    gat_dscr = &pend_dscr->sdescr;
    /*armci_vapi_print_dscr_info(gat_dscr,NULL);*/
    if((type==SERV && DEBUG_SERVER) || (type==CLN && DEBUG_CLN)){
       printf("%d: type(client=1)=%d inside posts gather dscr, id is %d\n",
	      armci_me,type,gat_dscr->wr_id);
       fflush(stdout);
    }

    rc = 0;
    if(type == CLN){
       rc = ibv_post_send((SRV_con+cluster)->qp, gat_dscr, &bad_wr);
    }
    else{
        rc = ibv_post_send((CLN_con + proc)->qp, gat_dscr, &bad_wr);
    }
    dassert1(1,rc==0,rc);

    THREAD_UNLOCK(armci_user_threads.net_lock);

}

/*\
 *  posts a bunch of gather descriptors
\*/ 
static sr_descr_t serv_blocking_gather_dscr;
static sr_descr_t client_blocking_gather_dscr;
void armci_post_gather(void *src_ptr, int src_stride_arr[], int count[],
      int stride_levels, armci_vapi_memhndl_t *mhandle,
      int proc,int nbtag, int type, sr_descr_t **srd)
{
    int i;
    int total_of_2D = 1;
    int total_size = 0;
    int index[MAX_STRIDE_LEVEL], unit[MAX_STRIDE_LEVEL];
    int j,k,y,z;
    int num_posted = 0;
    char *src, *src1;
    int num_xmit = 0, num_seg, max_seg, rem_seg,vecind;
    sr_descr_t *pend_dscr;

    struct ibv_sge *gat_sglist;
    struct ibv_send_wr *gat_dscr;

    if((type==SERV && DEBUG_SERVER) || (type==CLN && DEBUG_CLN)){
      printf("%d(%d)  : inside post_gather\n",armci_me,type);
      fflush(stdout);
    }

    max_seg =  armci_max_num_sg_ent;
    if(nbtag){
       pend_dscr = armci_vapi_get_next_sdescr(nbtag,1);
       if(srd!=NULL)*srd=pend_dscr;
    }
    else{
       pend_dscr = &client_blocking_gather_dscr;
       pend_dscr->sdescr.wr_id=DSCRID_SCATGAT + MAX_PENDING;
    }
    pend_dscr->numofsends=0;

    gat_dscr = &pend_dscr->sdescr;
    gat_sglist = pend_dscr->sg_entry;
    gat_dscr->opcode = IBV_WR_SEND;
    gat_dscr->send_flags = IBV_SEND_SIGNALED;
    gat_dscr->sg_list = gat_sglist;
    gat_dscr->num_sge = 0;
/*     gat_dscr->send_flags = 0; */

    index[2] = 0; unit[2] = 1;
    if(stride_levels > 1){
      total_of_2D = count[2];
      for(j=3; j<=stride_levels; j++){
        index[j] = 0; unit[j] = unit[j-1]*count[j-1];
        total_of_2D*=count[j];
      }
    }

    num_xmit = total_of_2D*count[1]/max_seg;
    rem_seg = (total_of_2D*count[1])%max_seg;
    if(num_xmit == 0) num_xmit = 1;
    else if(rem_seg!= 0)num_xmit++;

    if((type==SERV && DEBUG_SERVER) || (type==CLN && DEBUG_CLN) ){ 
       printf("%d(%d):armci_post_gather total_2D=%d, num_xmit=%d, rem_seg =%d, count[1] = %d\n",armci_me,type,total_of_2D, num_xmit,rem_seg,count[1]);
      fflush(stdout);
    }

    k=0; vecind = 0;
    if(rem_seg!=0 && k==(num_xmit-1))num_seg = rem_seg;
    else num_seg = max_seg;

    y=0,z=0;
    for(i=0;i<total_of_2D;i++){
       src = (char *)src_ptr;
       for(j=2;j<=stride_levels;j++){
         src+= index[j]*src_stride_arr[j-1];
         if(((i+1)%unit[j]) == 0) index[j]++;
         if(index[j] >= count[j]) index[j] =0;
       }
       src1 = src;

       for(j=0; j<count[1]; j++, vecind++){
         if(vecind == num_seg){
           posts_gather_desc(pend_dscr,proc,type);
           pend_dscr->numofsends++;

           /* the previous one has been posted, start off new*/
           gat_dscr->num_sge = 0;
           y = 0;
           vecind=0;total_size=0;k++;
           if(rem_seg!=0 && k==(num_xmit-1))num_seg = rem_seg;
         }

         /* fill the gather descriptor */
         gat_sglist[y].addr = (uint64_t)src1;
         gat_sglist[y].lkey = mhandle->lkey;
         gat_sglist[y].length = count[0];
         gat_dscr->num_sge++;
         src1 += src_stride_arr[0];
         y++;

       }

       if(vecind == num_seg){
         posts_gather_desc(pend_dscr,proc,type);
         pend_dscr->numofsends++;
         if((type==SERV && DEBUG_SERVER) || (type==CLN && DEBUG_CLN) ){
           printf("%d(%d)posts_gather_desc done\n",armci_me,type);
           fflush(stdout);
         }

         /* the previous one has been posted, start off new*/
         gat_dscr->num_sge = 0;
         y = 0;
         vecind = 0; total_size=0; k++;
         if(rem_seg!=0 && k==(num_xmit-1))num_seg=rem_seg;
         else num_seg = max_seg;
       }
    }
/*     printf("%d: num gathers posted =%d\n",armci_me,pend_dscr->numofsends); */
    if(!nbtag){
       /*complete here*/
       armci_send_complete(&pend_dscr->sdescr,"armci_post_gather",pend_dscr->numofsends);
    }
    return;
}
/***********************END SCATTER GATHER STUFF******************************/



/***********************SPECIAL SEND/RECV*************************************/
void armci_server_direct_send(int dst, char *src_buf, char *dst_buf, int len,
                              uint32_t *lkey, uint32_t *rkey)
{
    int rc = 0;
    struct ibv_wc *pdscr=NULL;
    struct ibv_wc pdscr1;
    struct ibv_send_wr sdscr;
    struct ibv_sge ssg_entry;

    pdscr = &pdscr1;

    if(DEBUG_SERVER){
       printf("\n%d(s):sending dir data to client %d at %p bytes=%d last=%p\n",
                armci_me,dst,dst_buf,len,(dst_buf+len-4));fflush(stdout);
    }

    memset(&sdscr,0,sizeof(struct ibv_send_wr));
    armci_init_cbuf_srdma(&sdscr,&ssg_entry,src_buf,dst_buf,len,NULL,NULL);
    sdscr.wr.rdma.rkey = *rkey;
    ssg_entry.lkey = *lkey;

    sdscr.wr_id = dst+armci_nproc;
    struct ibv_send_wr *bad_wr;
    rc = ibv_post_send((CLN_con+dst)->qp, &sdscr, &bad_wr);
    dassert1(1,rc==0,rc);

    while (rc == 0) {
       rc = ibv_poll_cq(CLN_nic->scq, 1, pdscr);
    }
    dassertp(1,rc>=0,("%d: rc=%d id=%d status=%d\n",
		      armci_me,rc,(int)pdscr->wr_id,pdscr->status));
    dassert1(1,pdscr->status==IBV_WC_SUCCESS,pdscr->status);;
}



void armci_send_contig_bypass(int proc, request_header_t *msginfo,
                              void *src_ptr, void *rem_ptr, int bytes)
{
    int *last;
    uint32_t *lkey=NULL;
    uint32_t *rkey;    
    int dscrlen = msginfo->dscrlen;

    last = (int*)(((char*)(src_ptr)) + (bytes - sizeof(int)));
    dassertp(1,msginfo->pinned,("%d: not pinned proc=%d",armci_me,proc));

    rkey = (uint32_t *)((char *)(msginfo+1)+dscrlen-(sizeof(uint32_t)+sizeof(uint32_t)));

    if(DEBUG_SERVER){
       printf("%d(server): sending data bypass to %d (%p,%p) %d %d\n", armci_me,
               msginfo->from,src_ptr, rem_ptr,*lkey,*rkey);
       fflush(stdout);
    }
    armci_server_direct_send(msginfo->from,src_ptr,rem_ptr,bytes,lkey,rkey);

    if(*last == ARMCI_STAMP){
       SERVER_SEND_ACK(msginfo->from);
    }
}

void armci_rcv_strided_data_bypass_both(int proc, request_header_t *msginfo,
                                       void *ptr, int *count, int stride_levels)
{
int datalen = msginfo->datalen;
int *last;
long *ack;
int loop=0;

    if(DEBUG_CLN){ printf("%d:rcv_strided_data_both bypass from %d\n",
                armci_me,  proc); fflush(stdout);
    }
    if(!stride_levels){
      last = (int*)(((char*)(ptr)) + (count[0] -sizeof(int)));
      ack  = (long *)&msginfo->tag;
      while(armci_util_int_getval(last) == ARMCI_STAMP &&
            armci_util_long_getval(ack)  != ARMCI_STAMP){
        loop++;
        loop %=1000000;
        if(loop==0){
          if(DEBUG_CLN){
            printf("%d: client last(%p)=%d ack(%p)=%ld off=%d\n",
                  armci_me,last,*last,ack,*ack,(int)((char*)last - (char*)ptr));
            fflush(stdout);
          }
        }
      }
    }
    else {
      printf("\n%d:rcv_strided_data called, it should never be called\n",armci_me);
      armci_dscrlist_recv_complete(0,"armci_rcv_strided_data_bypass_both",NULL);
    }

    if(DEBUG_CLN){printf("%d:rcv_strided_data bypass both: %d bytes from %d\n",
                          armci_me, datalen, proc); fflush(stdout);
    }
}


int armci_pin_memory(void *ptr, int stride_arr[], int count[], int strides)
{
    fprintf(stderr, "[%d]:armci_pin_memory not implemented\n",armci_me);
    fflush(stderr);
    return 0;
}


void armci_client_send_ack(int proc, int n)
{
    printf("\n%d:client_send_ack not implemented",armci_me);fflush(stdout);
}


void armci_rcv_strided_data_bypass(int proc, request_header_t* msginfo,
                                   void *ptr, int stride_levels)
{
    printf("\n%d:armci_rcv_strided_data_bypass not implemented",armci_me);
    fflush(stdout);
}


void armci_unpin_memory(void *ptr, int stride_arr[], int count[], int strides)
{
    printf("\n%d:armci_unpin_memory not implemented",armci_me);fflush(stdout);
}


int armcill_server_wait_ack(int proc, int n)
{
    printf("\n%d:armcill_server_wait_ack not implemented",armci_me);
    fflush(stdout);
    return(0);
}


void armcill_server_put(int proc, void* s, void *d, int len)
{
    printf("\n%d:armcill_server_put not implemented",armci_me);fflush(stdout);
}


/*\
 *  initialising the atomic send descriptor
\*/
void armci_init_vapibuf_atomic(struct ibv_send_wr *sd, struct ibv_sge *sg,
                   int op, int*ploc,int *prem, int extra,
                   int id,ARMCI_MEMHDL_T *lhandle,
                   ARMCI_MEMHDL_T *rhandle)
{
    if (1) {
       printf("%d(c) : entered armci_init_vapibuf_atomic\n",armci_me);
       fflush(stdout);
    }
    memset(sd,0,sizeof(struct ibv_send_wr));
    if (op == ARMCI_FETCH_AND_ADD_LONG ) {
       printf("%d(c) :setting opcode for snd dscr to FETCH_AND_ADD\n",armci_me);
       sd->opcode = IBV_WR_ATOMIC_FETCH_AND_ADD;
       sd->wr.atomic.compare_add = (uint64_t)extra;
    } else if(op == ARMCI_SWAP_LONG){
       sd->opcode = IBV_WR_ATOMIC_CMP_AND_SWP;
       sd->wr.atomic.swap = (uint64_t)extra;
    }
    sd->send_flags = IBV_SEND_SIGNALED;
    sg->length = 8; /* 64 bit atomic*/
    printf("--------\n");
    sg->addr= (uint64_t)(void *)ploc;
    if(lhandle)
    sg->lkey = lhandle->lkey;
    sd->sg_list = sg;
    sd->num_sge = 1;
    sd->wr.atomic.remote_addr = (uint64_t)(void *)prem;
    if(rhandle)
       sd->wr.atomic.rkey = rhandle->rkey; /* how do we get the remote key  */
    sd->wr_id = DSCRID_RMW + armci_me;

    if(1){
       printf("%d(c) : finished initialising atomic send desc id is %ld,armci_ime = %d\n",armci_me,sd->wr_id,armci_me);
       fflush(stdout);
    }   
}
/*\
 *   using vapi remote atomic operations
\*/
void client_rmw_complete(struct ibv_send_wr *snd_dscr, char *from)
{
    int rc = 0;
    struct ibv_wc pdscr1;
    struct ibv_wc *pdscr=&pdscr1;

  printf("%d(c) : inside client_rmw_complete\n",armci_me);
  do {
      while (rc == 0) {
	rc =  ibv_poll_cq(CLN_nic->scq, 1, pdscr);
      }
      dassertp(DBG_POLL|DBG_ALL,rc>=0,
	       ("%d: rc=%d id=%d status=%d\n",
		armci_me,rc,pdscr->wr_id,pdscr->status));
      dassert1(1,pdscr->status==IBV_WC_SUCCESS,pdscr->status);
      rc = 0;
    } while(pdscr->wr_id != snd_dscr->wr_id);
}


void armci_direct_rmw(int op, int*ploc, int *prem, int extra, int proc,
                      ARMCI_MEMHDL_T *lhandle, ARMCI_MEMHDL_T *rhandle)
{
    int rc = 0;
    struct ibv_send_wr *sd;
    struct ibv_sge *sg;
    vapi_nic_t *nic;
    armci_connect_t *con;

    nic = SRV_nic;
    con = CLN_con+proc;

    sd = &(rmw[armci_me].rmw_dscr);
    sg = &(rmw[armci_me].rmw_entry);

    if (1) {
        printf("%d(c) : about to call armci_init_vapibuf_atomic\n",armci_me);
        fflush(stdout);
    }

  armci_init_vapibuf_atomic(sd, sg, op,ploc,prem,extra,proc,lhandle,rhandle);

  if (1) {
     printf("%d(c) : finished armci_init_vapibuf_atomic\n",armci_me);
     fflush(stdout);
  }

  struct ibv_send_wr * bad_wr;
  rc = ibv_post_send(con->qp, sd, &bad_wr);
  dassert1(1,rc==0,rc);

  if (1) {
     printf("%d(c) : finished posting desc\n",armci_me);
     fflush(stdout);
  }

  /*armci_send_complete(sd,"send_remote_atomic");*/
  client_rmw_complete(sd,"send_remote_atomic");

  return;
}

struct node *dto_q = NULL;

void process_con_break_from_client(armci_ud_rank *h, cbuf *v)
{

    struct ibv_qp_attr qp_attr;
    struct ibv_qp_init_attr qp_init_attr;
    struct ibv_qp_cap qp_cap;
    enum ibv_qp_attr_mask qp_attr_mask;
    char *enval;
    struct ibv_recv_wr *bad_wr;

    int rc, j;
    armci_connect_t *con = CLN_con + h->src_rank;

    cbuf *v1 = get_cbuf();
    
    assert( v1 != NULL);
    
    v1->desc.u.sr.wr.ud.remote_qpn = rbuf.qp_num[h->src_rank];
    v1->desc.u.sr.wr.ud.remote_qkey = 0;
    v1->desc.u.sr.wr.ud.ah = conn.ud_ah[h->src_rank];
    
    armci_ud_rank *h1 = (armci_ud_rank *)CBUF_BUFFER_START(v1);
    h1->src_rank = armci_me;
    h1->qpnum = con->sqpnum;
    h1->msg_type = QP_CON_BREAK_FROM_SERVER;
    cbuf_init_send(v1, sizeof(armci_ud_rank));
    
    /* Release the receiving cbuf */
    release_cbuf(v);
    
    struct ibv_send_wr *send_wr; 
    if(ibv_post_send(conn.qp[0], &(v1->desc.u.sr), &send_wr)) {
        fprintf(stderr, "Error posting send\n");
    }
}


void process_con_break_from_server(armci_ud_rank *h, cbuf *v)
{

    int rc, j;
    armci_connect_t *con = SRV_con + armci_clus_id(h->src_rank);
    con->state = QP_INACTIVE;
    release_cbuf(v);
}

void process_recv_completion_from_client(armci_ud_rank *h, cbuf *v)
{
    
    struct ibv_qp_attr qp_attr;
    struct ibv_qp_init_attr qp_init_attr;
    struct ibv_qp_cap qp_cap;
    enum ibv_qp_attr_mask qp_attr_mask;
    char *enval;
    struct ibv_recv_wr *bad_wr;
    
    int rc, j;
    
    static int qp_flag = 1;
    
    if (qp_flag == 1) {
        total_active_conn_to_client = 0;
        qp_flag = 0;
    }
    
    total_active_conn_to_client++;
    double t_init_start, t_init_end, t_rtr_start, t_rts_end;
    static int flag = 0;
    
    t_init_start = MPI_Wtime();
    armci_connect_t *con = CLN_con + h->src_rank;
    
    if (con->state == QP_ACTIVE) {
        
        qp_attr_mask = IBV_QP_STATE;
        
        memset(&qp_attr, 0, sizeof qp_attr);
        qp_attr.qp_state        = IBV_QPS_SQD;
        
        if (ibv_modify_qp(con->qp, &qp_attr, qp_attr_mask)) {
            fprintf(stdout," Error modifying QP\n");
            fflush(stdout);
        }
        
        if (ibv_destroy_qp(con->qp)) {
            printf("Error destroying QP\n");
        }
        
        total_active_conn_to_client--;
    }
    armci_create_qp(CLN_nic, &con->qp);
    con->sqpnum = con->qp->qp_num;
    con->lid    = CLN_nic->lid_arr[h->src_rank];
    CLN_rqpnumtmpbuf[h->src_rank] = con->qp->qp_num;
    qp_attr_mask = IBV_QP_STATE
        | IBV_QP_PKEY_INDEX
        | IBV_QP_PORT
        | IBV_QP_ACCESS_FLAGS;

    memset(&qp_attr, 0, sizeof qp_attr);
    qp_attr.qp_state        = IBV_QPS_INIT;
    qp_attr.pkey_index      = DEFAULT_PKEY_IX;
    qp_attr.port_num        = CLN_nic->active_port;
    qp_attr.qp_access_flags = IBV_ACCESS_LOCAL_WRITE|
        IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_REMOTE_READ;

    rc = ibv_modify_qp(con->qp, &qp_attr, qp_attr_mask);

    t_init_end = MPI_Wtime();

    memset(&qp_attr, 0, sizeof qp_attr);
    qp_attr_mask = IBV_QP_STATE
        | IBV_QP_MAX_DEST_RD_ATOMIC
        | IBV_QP_PATH_MTU
        | IBV_QP_RQ_PSN
        | IBV_QP_MIN_RNR_TIMER;
    qp_attr.qp_state           = IBV_QPS_RTR;
    qp_attr.path_mtu           = IBV_MTU_1024;          /*MTU*/
    qp_attr.max_dest_rd_atomic = 4;
    qp_attr.min_rnr_timer      = RNR_TIMER;
    qp_attr.rq_psn             = 0;

    /* Fill in the information from the header */
    SRV_rqpnums[h->src_rank] = h->qpnum;

    qp_attr_mask |= IBV_QP_DEST_QPN | IBV_QP_AV;
    qp_attr.dest_qp_num  = SRV_rqpnums[h->src_rank];
    qp_attr.ah_attr.dlid = SRV_nic->lid_arr[h->src_rank];
    qp_attr.ah_attr.port_num = CLN_nic->active_port;

    rc = ibv_modify_qp(con->qp, &qp_attr, qp_attr_mask);

    memset(&qp_attr, 0, sizeof qp_attr);

    qp_attr_mask = IBV_QP_STATE
        | IBV_QP_SQ_PSN
        | IBV_QP_TIMEOUT
        | IBV_QP_RETRY_CNT
        | IBV_QP_RNR_RETRY
        | IBV_QP_MAX_QP_RD_ATOMIC;

    qp_attr.qp_state            = IBV_QPS_RTS;
    qp_attr.sq_psn              = 0;
    qp_attr.timeout             = 18;
    qp_attr.retry_cnt           = 7;
    qp_attr.rnr_retry           = 7;
    qp_attr.max_rd_atomic  = 4;
    rc = ibv_modify_qp(con->qp, &qp_attr,qp_attr_mask);
    t_rts_end = MPI_Wtime();

#if defined(PEND_BUFS)
    assert(armci_nproc*(IMM_BUF_NUM+1)<DSCRID_IMMBUF_RECV_END-DSCRID_IMMBUF_RECV);
    for(j=0; j<IMM_BUF_NUM+1; j++) {
        vapibuf_t *cbuf;
        cbuf = serv_buf_arr[h->src_rank*(IMM_BUF_NUM+1)+j];
        armci_init_vapibuf_recv(&cbuf->dscr, &cbuf->sg_entry, cbuf->buf,
                IMM_BUF_LEN, &serv_memhandle);
        cbuf->dscr.wr_id = h->src_rank*(IMM_BUF_NUM+1)+j + DSCRID_IMMBUF_RECV;
        if (DEBUG_SERVER) {
            printf("\n%d(s):posted rr with lkey=%d",armci_me,cbuf->sg_entry.lkey);
            fflush(stdout);
        }
        if (armci_use_srq) {
            rc = ibv_post_srq_recv(CLN_srq_hndl, &cbuf->dscr, &bad_wr);
        }
        else {
            rc = ibv_post_recv((CLN_con + h->src_rank)->qp, &cbuf->dscr, &bad_wr);
        }

        dassert1(1,rc==0,rc);
    }
#else
    int i;
    for(i =  0; i < armci_nproc; i++) {
        vapibuf_t *cbuf;
        cbuf = serv_buf_arr[h->src_rank];
        armci_init_vapibuf_recv(&cbuf->dscr, &cbuf->sg_entry, cbuf->buf,
                CBUF_DLEN, &serv_memhandle);
        cbuf->dscr.wr_id = h->src_rank+armci_nproc;
        if (DEBUG_SERVER) {
            printf("\n%d(s):posted rr with lkey=%d",armci_me,cbuf->sg_entry.lkey);
            fflush(stdout);
        }
        if (armci_use_srq) {
            rc = ibv_post_srq_recv(CLN_srq_hndl, &cbuf->dscr, &bad_wr);
        }
        else {
            rc = ibv_post_recv((CLN_con+h->src_rank)->qp, &cbuf->dscr, &bad_wr);
        }
        dassert1(1,rc==0,rc);
    }
#endif

    /* Now send back the information */

    struct cbuf *v1 = get_cbuf();

    assert( v1 != NULL);

    v1->desc.u.sr.wr.ud.remote_qpn = rbuf.qp_num[h->src_rank];
    v1->desc.u.sr.wr.ud.remote_qkey = 0;
    v1->desc.u.sr.wr.ud.ah = conn.ud_ah[h->src_rank];

    armci_ud_rank *h1 = (armci_ud_rank *)CBUF_BUFFER_START(v1);
    h1->src_rank = armci_me;
    h1->qpnum = con->sqpnum;
    h1->msg_type = 2;
    cbuf_init_send(v1, sizeof(armci_ud_rank));

    /* Release the receiving cbuf */
    release_cbuf(v);

    struct ibv_send_wr *send_wr;
    if(ibv_post_send(conn.qp[0], &(v1->desc.u.sr), &send_wr)) {
        fprintf(stderr, "Error posting send\n");
    }

    con->state = QP_ACTIVE;
}

void progress_engine()
{
    struct ibv_cq *ev_cq;
    void *ev_ctx;
    struct ibv_wc *sc;
    int ne;
    struct ibv_wc wc;
    int send_comp = 0;
    int recv_comp = 0;
    int count;

    void *cbuf_addr;

    cbuf *v;
    int ret;

    do {
        ne = ibv_poll_cq(hca.cq, 1, &wc);
    } while (ne < 1);

    /* Okay, got an entry, check for errors */
    if(ne < 0) {
        fprintf(stderr,"Error Polling CQ\n");
    }

    if(wc.status != IBV_WC_SUCCESS) {
        fprintf(stderr, "[%d] Failed status %d\n",
                armci_me, wc.status);
    }
    cbuf_addr = (void *) ((aint_t) wc.wr_id);
    assert(cbuf_addr != NULL);
    v = (cbuf *)cbuf_addr;
    armci_ud_rank *h = (armci_ud_rank *)(CBUF_BUFFER_START(v) + 40);

    if(IBV_WC_SEND == wc.opcode
            || IBV_WC_RDMA_WRITE == wc.opcode) {
        release_cbuf(v);
        /* Do nothing, just release the cbuf */
        /* Send Completion */
    } else if (IBV_WC_RECV == wc.opcode) {
        /* Recv completion */
        post_recv();

        /* Check if the message received is from a data server */
        assert((h->msg_type == QP_CON_REQ)
                || (h->msg_type == QP_CON_ACK)
                || (h->msg_type == QP_CON_BREAK_FROM_CLIENT)
                || (h->msg_type == QP_CON_BREAK_FROM_SERVER));
        if (h->msg_type == QP_CON_REQ) {
            process_recv_completion_from_client(h, v);
        }
        else if (h->msg_type == QP_CON_BREAK_FROM_CLIENT) {
            process_con_break_from_client(h, v);
        }
        else if (h->msg_type == QP_CON_BREAK_FROM_SERVER) {
            process_con_break_from_server(h, v);
        }
        else {
            process_recv_completion_from_server(h, v);
        }
    } else {
        fprintf(stderr, "Unknown opcode recv'd\n");
    }

}


void async_thread_ud_events(void *context)
{
    struct ibv_cq *ev_cq;
    void *ev_ctx;
    struct ibv_wc *sc;
    int ne;
    struct ibv_wc wc;
    int send_comp = 0;
    int recv_comp = 0;
    int count;

    void *cbuf_addr;

    cbuf *v;
    int ret;
    while (1) {
        progress_engine();
    }
}


static struct ibv_srq *create_srq(vapi_nic_t *nic)
{
    struct ibv_srq_init_attr srq_init_attr;
    struct ibv_srq *srq_ptr = NULL;

    memset(&srq_init_attr, 0, sizeof(srq_init_attr));

    srq_init_attr.srq_context = nic->handle;
#ifdef PEND_BUFS
    srq_init_attr.attr.max_wr = armci_nproc * (IMM_BUF_NUM + 1) + 200;
#else
    srq_init_attr.attr.max_wr = armci_nproc + 200;
#endif
    srq_init_attr.attr.max_sge = 1;
    /* The limit value should be ignored during SRQ create */
    srq_init_attr.attr.srq_limit = 30;

    srq_ptr = ibv_create_srq(nic->ptag, &srq_init_attr);
    return srq_ptr;
}

int init_params(void)
{
    conn.qp = (struct ibv_qp **) malloc(armci_nproc * sizeof(struct ibv_qp *));
    conn.lid = (uint16_t *) malloc(armci_nproc * sizeof(int));
    conn.qp_num = (uint32_t *) malloc(armci_nproc * sizeof(int));
    conn.ud_ah = (struct ibv_ah **) malloc (armci_nproc * sizeof (struct ibv_ah *));
    conn.status = (int *) malloc(armci_nproc * sizeof(int));

    memset((void *)conn.status, 0, sizeof(int) * armci_nproc);

    rbuf.qp_num = (uint32_t *) malloc(armci_nproc * sizeof(int));
    rbuf.lid = (uint16_t *) malloc(armci_nproc * sizeof(int));
    rbuf.rkey = (uint32_t *) malloc(armci_nproc * sizeof(int));
    rbuf.buf = (char **) malloc(armci_nproc * sizeof(char *));

    assert(conn.qp && conn.lid && rbuf.qp_num && rbuf.lid && rbuf.rkey && rbuf.buf);
    return 0;
}


int open_hca(void)
{
    struct ibv_device **dev_list;
    struct ibv_device *ib_dev;

    int num_hcas;
    dev_list = ibv_get_device_list(&num_hcas);


    hca.ib_dev = dev_list[0];

    hca.context = ibv_open_device(hca.ib_dev);

    if(!hca.context) {
        fprintf(stderr,"Couldn't get context %s\n",
                ibv_get_device_name(hca.ib_dev));
        return 1;
    }

    hca.pd = ibv_alloc_pd(hca.context);

    if(!hca.pd) {
        fprintf(stderr,"Couldn't get pd %s\n",
                ibv_get_device_name(hca.ib_dev));
        return 1;
    }
    return 0;
}

int create_cq(void)
{
    hca.cq = ibv_create_cq(hca.context, 16000, NULL,
                    NULL, 0);
    if(!hca.cq) {
        fprintf(stderr, "Couldn't create CQ\n");
        return 1;
    }
    if (armci_me == armci_master) {
        pthread_create(&armci_async_thread[1], NULL,
                (void *) async_thread_ud_events, (void *) hca.context);
    }

    return 0;
}

int get_lid(void)
{
    struct ibv_port_attr port_attr[2];
    int i, j, count;
    int active_port_search_count;

    for (j = 1; j <= 1; j++) {
        if (!ibv_query_port(hca.context, j, &port_attr[j - 1])
                && (port_attr[j - 1].state == IBV_PORT_ACTIVE)) {
            for (i = 0; i < armci_nproc; i++)
                conn.lid[i] = port_attr[j - 1].lid;
            return 0;
        }
    }

    return 1;
}

int exch_addr(void)
{
    MPI_Status status;
    int i;

    MPI_Allgather((void *)conn.qp_num, sizeof(uint32_t), MPI_BYTE,
            (void *)rbuf.qp_num, sizeof(uint32_t), MPI_BYTE, ARMCI_COMM_WORLD);
    MPI_Alltoall((void *)conn.lid, sizeof(uint16_t), MPI_BYTE,
            (void *)rbuf.lid, sizeof(uint16_t), MPI_BYTE, ARMCI_COMM_WORLD);

#ifdef DEBUG
    for (i = 0; i < nprocs; i++) {
        if (me == i)
            continue;
        fprintf(stdout,"[%d] Remote QP %d, Remote LID %u, Rkey %u, Lkey %u\n"
                " LBuf %p, RBuf %p\n",
                me, rbuf.qp_num[i], rbuf.lid[i], rbuf.rkey[i], lbuf.mr->lkey,
                lbuf.buf, rbuf.buf[i]);
        fflush(stdout);
    }
#endif

    return 0;
}

int create_qp(void)
{
    struct ibv_qp_attr qp_attr;
    int i;

    memset(&qp_attr, 0, sizeof qp_attr);

    struct ibv_qp_init_attr attr = {
        .send_cq = hca.cq,
        .recv_cq = hca.cq,
        .cap     = {
            .max_send_wr  = 8192,
            .max_recv_wr  = 8192,
            .max_send_sge = 1,
            .max_recv_sge = 1,
            .max_inline_data = 0
        },
        .qp_type = IBV_QPT_UD
    };

    conn.qp[0] = ibv_create_qp(hca.pd, &attr);
    if(!conn.qp[0]) {
        fprintf(stderr,"Couldn't create QP\n");
        return 1;
    }

    conn.qp_num[0] = conn.qp[0]->qp_num;
    qp_attr.qp_state = IBV_QPS_INIT;
    qp_attr.pkey_index = 0;
    qp_attr.port_num   = 1;
    qp_attr.qkey = 0;

    if(ibv_modify_qp(conn.qp[0], &qp_attr,
                IBV_QP_STATE              |
                IBV_QP_PKEY_INDEX         |
                IBV_QP_PORT               |
                IBV_QP_QKEY)) {
        fprintf(stderr,"Could not modify QP to INIT\n");
        return 1;
    }
#ifdef DEBUG
    fprintf(stdout,"[%d] Created QP %d, LID %d\n", me,
            conn.qp_num[0], conn.lid[0]);
    fflush(stdout);
#endif

    return 0;
}
struct ibv_mr* armci_register_memory(void* buf, int len)
{
    return (ibv_reg_mr(hca.pd, buf, len,
                IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_WRITE |
                IBV_ACCESS_REMOTE_READ));
}

int connect_qp(void)
{
    struct ibv_qp_attr attr;
    int i;
    memset(&attr, 0 , sizeof attr);

    attr.qp_state       = IBV_QPS_RTR;
    if (ibv_modify_qp(conn.qp[0], &attr,
                IBV_QP_STATE)) {
        fprintf(stderr, "Failed to modify QP to RTR\n");
        return 1;
    }

    attr.qp_state       = IBV_QPS_RTS;
    attr.sq_psn         = 0;
    if (ibv_modify_qp(conn.qp[0], &attr,
                IBV_QP_STATE              |
                IBV_QP_SQ_PSN)) {
        fprintf(stderr, "Failed to modify QP to RTS\n");
        return 1;
    }
    for (i = 0; i < armci_nproc; i++) {
        struct ibv_ah_attr ah_attr;
        memset(&ah_attr, 0, sizeof(ah_attr));
        ah_attr.is_global = 0;
        ah_attr.dlid = rbuf.lid[i];
        ah_attr.sl = 0;
        ah_attr.src_path_bits = 0;
        ah_attr.port_num = 1;

        conn.ud_ah[i] = ibv_create_ah(hca.pd, &ah_attr);

        if (!conn.ud_ah[i]) {
            fprintf(stderr, "Error creating address handles\n");
        }
    }

    return 0;
}

int total_ud_recv_buffers = 4096;


void post_recv()
{
    cbuf *v = get_cbuf();
    cbuf_init_recv(v, CBUF_BUFFER_SIZE);
    struct ibv_recv_wr *bad_wr;

    /* set id to be the global_rank */
    v->grank = -1;

    if(ibv_post_recv(conn.qp[0], &(v->desc.u.rr), &bad_wr)) {
        fprintf(stderr," Error posting UD recv\n");
        fflush(stderr);
    }
}


int post_buffers()
{
    init_cbuf_lock();
    allocate_cbufs(8192);
    int i;
    for (i = 0; i < total_ud_recv_buffers; i++) {
        post_recv();
    }
    return 0;
}

void test_connectivity(void)
{
    int i;
    struct ibv_send_wr *bad_wr;

    cbuf *v = NULL;
    for (i = 0; i < armci_nproc;i++) {
        if (armci_me == i)
            continue;

        int j;

        for (j = 0 ; j < 2; j++) {
            v = get_cbuf();
            assert(v != NULL);
            v->desc.u.sr.wr.ud.remote_qpn = rbuf.qp_num[i];
            v->desc.u.sr.wr.ud.remote_qkey = 0;
            v->desc.u.sr.wr.ud.ah = conn.ud_ah[i];

            armci_ud_rank *h = (armci_ud_rank *)CBUF_BUFFER_START(v);
            h->src_rank = armci_me;
            cbuf_init_send(v, sizeof(armci_ud_rank));
            if(ibv_post_send(conn.qp[0], &(v->desc.u.sr), &bad_wr)) {
                fprintf(stderr, "Error posting send\n");
            }
        }
    }
}

int recycle_dead_qp()
{
}

void handle_network_fault(struct ibv_wc *pdscr)
{

    recycle_dead_qp();
}


void setup_ud_channel()
{
    if(init_params()) {
    }

    if(open_hca()) {
    }

    if(create_cq()) {
    }

    if(get_lid()) {
    }

    if(create_qp()) {
    }

    if (exch_addr()) {
    }

    if(connect_qp()) {
    }

    MPI_Barrier(ARMCI_COMM_WORLD);

    if(post_buffers()) {
    }

    MPI_Barrier(ARMCI_COMM_WORLD);
}


void process_recv_completion_from_server(armci_ud_rank *h, cbuf *v)
{

    struct ibv_qp_attr qp_attr;
    struct ibv_qp_init_attr qp_init_attr;
    struct ibv_qp_cap qp_cap;
    enum ibv_qp_attr_mask qp_attr_mask;
    char *enval;
    struct ibv_recv_wr *bad_wr;

    int rc, j;

    armci_connect_t *con = SRV_con + armci_clus_id(h->src_rank);

    qp_attr_mask = IBV_QP_STATE
        | IBV_QP_PKEY_INDEX
        | IBV_QP_PORT
        | IBV_QP_ACCESS_FLAGS;

    memset(&qp_attr, 0, sizeof qp_attr);
    qp_attr.qp_state        = IBV_QPS_INIT;
    qp_attr.pkey_index      = DEFAULT_PKEY_IX;
    qp_attr.port_num        = SRV_nic->active_port;
    qp_attr.qp_access_flags = IBV_ACCESS_LOCAL_WRITE|
        IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_REMOTE_READ;

    rc = ibv_modify_qp(con->qp, &qp_attr, qp_attr_mask);

    memset(&qp_attr, 0, sizeof qp_attr);
    qp_attr_mask = IBV_QP_STATE
        | IBV_QP_MAX_DEST_RD_ATOMIC
        | IBV_QP_PATH_MTU
        | IBV_QP_RQ_PSN
        | IBV_QP_MIN_RNR_TIMER;
    qp_attr.qp_state           = IBV_QPS_RTR;
    qp_attr.path_mtu           = IBV_MTU_1024;          /*MTU*/
    qp_attr.max_dest_rd_atomic = 4;
    qp_attr.min_rnr_timer      = RNR_TIMER;
    qp_attr.rq_psn             = 0;

    /* Fill in the information from the header */
    CLN_rqpnums[h->src_rank] = h->qpnum;

    qp_attr_mask |= IBV_QP_DEST_QPN | IBV_QP_AV;
    qp_attr.dest_qp_num  = CLN_rqpnums[h->src_rank];
    qp_attr.ah_attr.dlid = SRV_nic->lid_arr[h->src_rank];
    qp_attr.ah_attr.port_num = SRV_nic->active_port;
    rc = ibv_modify_qp(con->qp, &qp_attr, qp_attr_mask);

    memset(&qp_attr, 0, sizeof qp_attr);

    qp_attr_mask = IBV_QP_STATE
        | IBV_QP_SQ_PSN
        | IBV_QP_TIMEOUT
        | IBV_QP_RETRY_CNT
        | IBV_QP_RNR_RETRY
        | IBV_QP_MAX_QP_RD_ATOMIC;

    qp_attr.qp_state            = IBV_QPS_RTS;
    qp_attr.sq_psn              = 0;
    qp_attr.timeout             = 18;
    qp_attr.retry_cnt           = 7;
    qp_attr.rnr_retry           = 7;
    qp_attr.max_rd_atomic  = 4;

    rc = ibv_modify_qp(con->qp, &qp_attr,qp_attr_mask);
   release_cbuf(v);

    con->state = QP_ACTIVE;
}

struct con_q_t {
    struct armci_connect_t *head;
    void *next;
};

struct con_q_t *con_q;

#define MAX_CLIENT_TO_SERVER_CONN 2


armci_connect_t * dequeue_conn()
{
    return NULL;
}

int get_max_client_to_server_conn()
{
    static int env_read_flag = 0;

    static int max_client_to_server_conn = 0;
    if (!env_read_flag) {
        char *value;
        if ((value = getenv("ARMCI_MAX_CLIENT_TO_SERVER_FACTOR")) != NULL){
            max_client_to_server_conn = armci_nclus/atoi(value);
            /* We need a minimum of 4 connections */
            if (max_client_to_server_conn <= 1)
                max_client_to_server_conn = 2;

            if (armci_me == 0) {
                fprintf(stdout, "max_client_to_server_conn[%d]\n",
                        max_client_to_server_conn);
            }
        } else {
            max_client_to_server_conn = armci_nclus / 2;
            if (max_client_to_server_conn <= 1)
                max_client_to_server_conn = 2;
        }
    }

    return max_client_to_server_conn;

}

int get_the_victim_connection()
{
    return 0;
}

void break_a_connection_if_needed()
{
    assert(!SERVER_CONTEXT);

    if (!armci_use_lazy_break)
        return;

    int victim, proc;

    armci_connect_t *con;
    int max_client_to_server_conn = get_max_client_to_server_conn();

    if ((total_active_conn_to_server >=
                max_client_to_server_conn) && (armci_me != armci_master)) {

        int proc, clus_id;
        do {
            proc = get_the_victim_connection();
            /* Not enough on the queue */
            if (proc == -1)
                return;
            clus_id = armci_clus_id(proc);
        } while ((armci_clus_me == clus_id) ||
                ((SRV_con + clus_id)->state != QP_ACTIVE));


        double t_fence_start, t_fence_end;

        t_fence_start = MPI_Wtime();

        ARMCI_WaitAll();

        ARMCI_Fence(proc);

        struct ibv_qp_attr qp_attr;
        struct ibv_qp_init_attr qp_init_attr;
        struct ibv_qp_cap qp_cap;
        enum ibv_qp_attr_mask qp_attr_mask;
        char *enval;

        double t_destroy_start, t_destroy_end;

        t_destroy_start = MPI_Wtime();
        armci_connect_t *con = SRV_con + clus_id;
        int rc, j;

        dassert(1, con->state == QP_ACTIVE);
        qp_attr_mask = IBV_QP_STATE;

        memset(&qp_attr, 0, sizeof qp_attr);
        qp_attr.qp_state        = IBV_QPS_SQD;

        if (ibv_modify_qp(con->qp, &qp_attr, qp_attr_mask)) {
            fprintf(stdout," Error modifying QP\n");
            fflush(stdout);
        }
        cbuf *v = get_cbuf();
        assert(v != NULL);
        v->desc.u.sr.wr.ud.remote_qpn =
            rbuf.qp_num[armci_clus_info[clus_id].master];
        v->desc.u.sr.wr.ud.remote_qkey = 0;
        v->desc.u.sr.wr.ud.ah = conn.ud_ah[armci_clus_info[clus_id].master];

        armci_ud_rank *h = (armci_ud_rank *)CBUF_BUFFER_START(v);
        h->src_rank = armci_me;
        h->qpnum = con->sqpnum;

        h->msg_type = QP_CON_BREAK_FROM_CLIENT;

        struct ibv_send_wr *bad_wr;

        cbuf_init_send(v, sizeof(armci_ud_rank));

        if(ibv_post_send(conn.qp[0], &(v->desc.u.sr), &bad_wr)) {
            printf("Error posting send\n");
        }


        while (con->state != QP_INACTIVE) {
            if (armci_me != armci_master)
                progress_engine();
        }
        if (ibv_destroy_qp(con->qp)) {
            printf("Error destroying QP\n");
        }
    }
}


void check_state_of_ib_connection(int proc, int force)
{
    int clus_id;
    armci_connect_t *con;

    /* return if ARMCI does not use on demand connection management */
    if (!armci_use_odcm)
        return;

    static int flag = 1;

    static int last = 0;

    if (flag == 1) {
        total_active_conn_to_server = 0;
        total_breaks = 0;
        flag = 0;
    }
    /* Check clus id */
    clus_id = armci_clus_id(proc);

    con = SRV_con + clus_id;

    assert(!SERVER_CONTEXT);

    if (con->state != QP_ACTIVE) {

        if (!force)
            break_a_connection_if_needed();

        total_active_conn_to_server++;

        double t_create_start, t_create_end;

        t_create_start = MPI_Wtime();
        armci_create_qp(SRV_nic, &con->qp);
        con->sqpnum  = con->qp->qp_num;
        con->lid = SRV_nic->lid_arr[clus_id];

        cbuf *v = get_cbuf();
        assert(v != NULL);
        v->desc.u.sr.wr.ud.remote_qpn =
            rbuf.qp_num[armci_clus_info[clus_id].master];
        v->desc.u.sr.wr.ud.remote_qkey = 0;
        v->desc.u.sr.wr.ud.ah = conn.ud_ah[armci_clus_info[clus_id].master];

        armci_ud_rank *h = (armci_ud_rank *)CBUF_BUFFER_START(v);
        h->src_rank = armci_me;
        h->qpnum = con->sqpnum;

        h->msg_type = QP_CON_REQ;

        struct ibv_send_wr *bad_wr;

        cbuf_init_send(v, sizeof(armci_ud_rank));

        if(ibv_post_send(conn.qp[0], &(v->desc.u.sr), &bad_wr)) {
            printf("Error posting send\n");
        }

        if (armci_me != armci_master) {
            while (con->state != QP_ACTIVE) {
                progress_engine();
            }
        } else {
            while (con->state != QP_ACTIVE) {
                usleep(1);
            }
        }
        t_create_end = MPI_Wtime();
    }

}


