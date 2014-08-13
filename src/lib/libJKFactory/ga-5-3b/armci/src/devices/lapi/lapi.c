#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* initialization of data structures and setup of lapi internal parameters */ 
 
#include <pthread.h>
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRINGS_H
#   include <strings.h>
#endif
#include "lapidefs.h"
#include "armcip.h"
#include "copy.h"
#ifdef AIX
#include <sys/atomic_op.h>
#endif
#if (defined(PPC) || defined(__PPC__) || defined(__PPC))
#   include "asm-ppc.h"
#endif

#define DEBUG_ 0
#define ERROR(str,val) armci_die((str),(val))
#define BUF_TO_EVBUF(buf) ((lapi_cmpl_t*)(((char*)buf) - sizeof(lapi_cmpl_t)))

char err_msg_buf[LAPI_MAX_ERR_STRING]; /* for error msg returned by LAPI */

/*
** macro to check return code of function calls. keeps return
** code checking logic from needing to be in main logic
*/
#define CHECK(func_and_args)                                      \
{                                                                 \
    int rc;                                                       \
    if ((rc = (func_and_args)) != LAPI_SUCCESS) {                 \
        LAPI_Msg_string(rc, err_msg_buf);                         \
        fprintf(stderr,                                           \
                "LAPI ERROR: %s, rc = %d\n", err_msg_buf, rc);    \
        armci_die("LAPI Error", 0);                               \
    }                                                             \
}

#if ARMCI_ENABLE_GPC_CALLS
    extern gpc_buf_t *gpc_req;
#endif
int lapi_max_uhdr_data_sz; /* max  data payload */
lapi_cmpl_t *cmpl_arr;     /* completion state array, dim=NPROC */
lapi_cmpl_t  hdr_cntr;     /* AM header buffer counter  */
lapi_cmpl_t  buf_cntr;     /* AM data buffer counter    */
lapi_cmpl_t*  ack_cntr;    /* ACK counter used in handshaking protocols
                              between origin and target */
lapi_cmpl_t*  get_cntr;    /* counter used with lapi_get  */
lapi_user_cxt_t *lapi_remote_cxt; /* Remote context for RDMA call */

int intr_status;
lapi_info_t     lapi_info;
#ifndef TCGMSG
lapi_handle_t   lapi_handle;
#endif
pthread_mutex_t _armci_mutex_thread=PTHREAD_MUTEX_INITIALIZER;


double _armci_rcv_buf[MSG_BUFLEN_DBL];
char* MessageRcvBuffer = (char*)_armci_rcv_buf;
char* MessageSndBuffer = (char*)0;

extern void armci_waitsome(int factor);

/************* LAPI Active Message handlers *******************************/

volatile static int hndlcnt=0, header_cnt=0;
static int hhnum=0;
static long num_malloc=0;   /* trace and limit the number malloc calls in HH */
#define MAX_NUM_MALLOC 100

/* trace state of accumulate lock */
int kevin_ok=1; /* "1" indicates that no other thread is holding the lock */


void armci_completion_handler(lapi_handle_t *t_hndl, void *save)
{
    lapi_handle_t hndl = *t_hndl;
    int need_data;
    void *message;
    int whofrom, msglen;
    request_header_t *msginfo = (request_header_t *)save;
    char *descr= (char*)(msginfo+1), *buf=MessageRcvBuffer;
    int buflen=MSG_BUFLEN;
#if ARMCI_ENABLE_GPC_CALLS
    extern pthread_t data_server;
    data_server = pthread_self();
#endif

    if(DEBUG_)
        fprintf(stderr,"%d:CH:op=%d from=%d datalen=%d dscrlen=%d\n", armci_me,
                msginfo->operation, msginfo->from,msginfo->datalen,msginfo->dscrlen);

    /*** assure that descriptor and data are in the right format and place ***/
    if(  msginfo->dscrlen < 0 || msginfo->datalen <0 ){
        /* for large put/acc/scatter need to get the data */
        int rc;
        lapi_cntr_t req_cntr;    
        int bytes=0;
        char *origin_ptr = msginfo->tag.buf;

        if (msginfo->dscrlen<0) {
            descr =MessageRcvBuffer;
            msginfo->dscrlen = -msginfo->dscrlen;
            buf = descr + msginfo->dscrlen;
            buflen += msginfo->dscrlen;
            bytes += msginfo->dscrlen;

        }
        if (msginfo->datalen <0){
            msginfo->datalen = -msginfo->datalen;
            bytes += msginfo->datalen;
        }

        if(rc=LAPI_Setcntr(hndl, &req_cntr, 0)) ERROR("CH:setcntr failed",rc);
        if(rc=LAPI_Get(hndl, (uint)msginfo->from, bytes,
                    origin_ptr, MessageRcvBuffer,
                    msginfo->tag.cntr,&req_cntr))ERROR("CH:LAPI_Get failed",rc);

        if(rc=LAPI_Waitcntr(hndl, &req_cntr,1,NULL))ERROR("CH:Waitcntr failed",rc);


    } else{

        /* desc is in save, data could be but not for GET */
        if(msginfo->operation !=GET)buf = descr + msginfo->dscrlen;
        buflen = MSG_BUFLEN;
    }

    /*   fprintf(stderr,"CH: val=%lf\n",*(double*)(buf+msginfo->datalen -8));*/


    /*** dispatch request to the appropriate handler function ***/
    switch(msginfo->operation){
        case LOCK:   armci_server_lock(msginfo); 
                     break;
        case UNLOCK: armci_server_unlock(msginfo, descr); 
                     break;
        default:
                     if(msginfo->format == STRIDED)
                         armci_server(msginfo, descr, buf, buflen);
                     else
                         armci_server_vector(msginfo, descr, buf, buflen);
    }

    free(msginfo);
#ifdef LINUX
    (void)fetch_and_add(&num_malloc, (long)-1);
#else
    (void)fetch_and_addlp(&num_malloc, (long)-1);
#endif
}




void* armci_header_handler(lapi_handle_t *t_hndl, void *uhdr, uint *t_uhdrlen,
        uint *msglen, compl_hndlr_t **handler, void** psave)
{
    lapi_handle_t hndl = *t_hndl;
    uint uhdrlen = *t_uhdrlen;
    request_header_t *msginfo = (request_header_t *)uhdr;

    if(DEBUG_)
        fprintf(stderr,"%d:HH: op=%d from %d\n",armci_me,msginfo->operation,
                msginfo->from);
    if(msginfo->to != armci_me)
        armci_die("wrong message delivered",msginfo->to);

    /* process small requests that do not require comms in header handler */
    if(msginfo->datalen >0 && msginfo->dscrlen>0 && msginfo->operation != GET 
            && msginfo->operation != LOCK && msginfo->operation != UNLOCK){

        /* If another thread is in accumulate use compl. handler path:
         * Try to avoid blocking inside HH which degrades Lapi performance.
         * The completion handler path requires malloc to save request info.
         * Only up to approx. MAX_NUM_MALLOC requests can be rescheduled to
         * run in CH instead of HH. 
         * MAX_NUM_MALLOC is a soft limit to avoid cost of locking when reading 
         */

        if( msginfo->operation==PUT || num_malloc>MAX_NUM_MALLOC || kevin_ok){

            char *descr = (char*)(msginfo+1);
            char *buf   = descr + msginfo->dscrlen;
            int buflen = uhdrlen - sizeof(request_header_t) - msginfo->dscrlen;

            if(DEBUG_)
                fprintf(stderr,"%d:HH: buf =%lf\n",armci_me,*(double*)buf);
            if(msginfo->format == STRIDED)
                armci_server(msginfo, descr, buf, buflen);
            else
                armci_server_vector(msginfo, descr, buf, buflen);

            /*             fprintf(stderr,"%d:HH: getting out of server\n",armci_me);*/
            *psave = NULL;
            *handler = NULL;
            return(NULL);
        }
    }

#ifdef LINUX
    (void)fetch_and_add(&num_malloc, (long)1);
#else
    (void)fetch_and_addlp(&num_malloc, (long)1); /* AIX atomic increment */
#endif

    msginfo  = (request_header_t*) malloc(uhdrlen); /* recycle pointer */
    if(!msginfo) ERROR("HH: malloc failed in header handler",num_malloc);

    /* save the request info for processing in compl. handler */
    memcpy((char*)msginfo, uhdr, uhdrlen);
    *psave = msginfo;
    *handler = armci_completion_handler;

    return(NULL);
}


void armci_send_req(int proc, request_header_t* msginfo, int len)
{
    int msglen = sizeof(request_header_t);
    lapi_cntr_t *pcmpl_cntr, *pcntr = &(BUF_TO_EVBUF(msginfo)->cntr);
    int rc;

    msginfo->tag.cntr= pcntr;
#if ARMCI_ENABLE_GPC_CALLS
    if(msginfo->operation==GET && msginfo->format==VECTOR && msginfo->ehlen){ 
        msginfo->tag.buf = (char *)(msginfo+1)+msginfo->dscrlen;
    }
    else 
#endif
        msginfo->tag.buf = msginfo+1;

    if(msginfo->operation==GET || msginfo->operation==LOCK){

        SET_COUNTER(*(lapi_cmpl_t*)pcntr,1);/*dataarrive in same buf*/
        /*The GPC case. Note that we don't use the parameter len*/
        if(msginfo->format==VECTOR && msginfo->ehlen > 0) 
            msglen += msginfo->datalen;
        if(lapi_max_uhdr_data_sz < msginfo->dscrlen){

            msginfo->dscrlen = -msginfo->dscrlen; /* no room for descriptor */
            pcntr = NULL; /* GET(descr) from CH will increment buf cntr */

        }else msglen += msginfo->dscrlen;

        /*
           we should send the mutex, too. When op==LOCK, Value of len parameter
           is already sizeof(reqest_header_t)+sizeof(int), since we dont use 
           len but construct our own msglen, we need to add sizeof(int).
           */
        if(msginfo->operation==LOCK) msglen += sizeof(int);

        pcmpl_cntr=NULL; /* don't trace completion status for load ops */

    }else if (msginfo->operation==UNLOCK){

        msglen += msginfo->dscrlen;
        pcmpl_cntr=NULL; /* don't trace completion status for unlock */

    }else{

        if(lapi_max_uhdr_data_sz < (msginfo->datalen + msginfo->dscrlen)){

            msginfo->datalen = -msginfo->datalen;
            msginfo->dscrlen = -msginfo->dscrlen;
            pcntr = NULL; /* GET/LOCK from CH will increment buf cntr */

        }else msglen += msginfo->dscrlen+msginfo->datalen;

        /* trace completion of store ops */
        pcmpl_cntr = &cmpl_arr[msginfo->to].cntr; 

    }

    if(msginfo->operation==PUT || ARMCI_ACC(msginfo->operation)) 
        UPDATE_FENCE_STATE(msginfo->to, msginfo->operation, 1);

    if((rc=LAPI_Amsend(lapi_handle,(uint)msginfo->to,
                    (void*)armci_header_handler, msginfo, msglen, NULL, 0,
                    NULL, pcntr, pcmpl_cntr))) armci_die("AM failed",rc);

    if(DEBUG_) fprintf(stderr,"%d sending req=%d to %d\n",
            armci_me, msginfo->operation, proc);
}



/*\ client sends strided data + request to server
  \*/
void armci_send_strided(int proc, request_header_t *msginfo, char *bdata,
        void *ptr, int strides, int stride_arr[], int count[])
{

    armci_write_strided(ptr, strides, stride_arr, count, bdata);
    armci_send_req(proc,msginfo,msginfo->bytes + sizeof(request_header_t));
}


/*\ server sends data back to client
  \*/
void armci_send_data(request_header_t* msginfo, void *data)
{
    armci_lapi_send(msginfo->tag, data, msginfo->datalen, msginfo->from);
}


/*\ server sends strided data back to client
  \*/
void armci_send_strided_data(int proc,  request_header_t *msginfo, char *bdata,
        void *ptr, int strides, int stride_arr[], int count[])
{
    armci_write_strided(ptr, strides, stride_arr, count, bdata);
    armci_lapi_send(msginfo->tag, bdata, msginfo->datalen, msginfo->from);
}


char* armci_rcv_data(int proc, request_header_t *msginfo)
{
    lapi_cmpl_t *pcntr=BUF_TO_EVBUF(msginfo);
    CLEAR_COUNTER((*pcntr));
#if ARMCI_ENABLE_GPC_CALLS
    if(msginfo->operation==GET && msginfo->format==VECTOR && msginfo->ehlen){
        return((char *)(msginfo+1)+msginfo->dscrlen);
    }
    else
#endif
        return (char*)(msginfo+1);
}



/*\ client receives strided data from server
  \*/
void armci_rcv_strided_data(int proc, request_header_t* msginfo, int datalen,
        void *ptr, int strides, int stride_arr[], int count[])
{
    lapi_cmpl_t *pcntr=BUF_TO_EVBUF(msginfo);
    CLEAR_COUNTER((*pcntr));
    armci_read_strided(ptr, strides, stride_arr, count, (char*)(msginfo+1));
}



/*\ client receives vector data from server to buffer and unpacks it            
  \*/
void armci_rcv_vector_data(int proc, request_header_t* msginfo, 
        armci_giov_t darr[], int len)
{
    char *buf;
    buf = armci_rcv_data(proc,msginfo);
    armci_vector_from_buf(darr, len, buf);
}


/*\ initialization of LAPI related data structures
  \*/
void armci_init_lapi()
{
    int rc, p;
    int lapi_max_uhdr_sz;
    lapi_cmpl_t *pcntr;
    lapi_remote_cxt_t util_cxt;  /* For call to obtain rCxt */

#ifndef TCGMSG
    rc = LAPI_Init(&lapi_handle, &lapi_info);
    if(rc) ERROR("lapi_init failed",rc);
#endif

    /* set the max limit for AM header data length */
    rc = LAPI_Qenv(lapi_handle,MAX_UHDR_SZ, &lapi_max_uhdr_sz);
    if(rc) ERROR("armci_init_lapi:  LAPI_Qenv failed", rc); 

    /*     fprintf(stderr,"max header size = %d\n",lapi_max_uhdr_sz);*/

    /* how much data can fit into AM header ? */
    lapi_max_uhdr_data_sz = lapi_max_uhdr_sz - sizeof(request_header_t);

    /* allocate memory for completion state array */
    cmpl_arr = (lapi_cmpl_t*)malloc(armci_nproc*sizeof(lapi_cmpl_t));
    if(cmpl_arr==NULL) ERROR("armci_init_lapi:malloc for cmpl_arr failed",0);

    /* allocate memory for ack and get counters, 1 if not thread safe */
#ifdef THREAD_SAFE
    ack_cntr = calloc(armci_user_threads.max, sizeof(lapi_cmpl_t));
    get_cntr = calloc(armci_user_threads.max, sizeof(lapi_cmpl_t));
#else
    ack_cntr = calloc(1, sizeof(lapi_cmpl_t));
    get_cntr = calloc(1, sizeof(lapi_cmpl_t));
#endif
    if (!(ack_cntr && get_cntr))
        ERROR("armci_init_lapi:calloc for ack or get counters failed",0);

    /* initialize completion state array */
    for(p = 0; p< armci_nproc; p++){
        rc = LAPI_Setcntr(lapi_handle, &cmpl_arr[p].cntr, 0);
        if(rc) ERROR("armci_init_lapi: LAPI_Setcntr failed (arr)",rc);
        cmpl_arr[p].oper = -1;
        cmpl_arr[p].val = 0;
    }

    /* initialize ack/buf/hdr counters */
#ifdef THREAD_SAFE
#   define N armci_user_threads.max
#else
#   define N 1
#endif
    for (p = 0; p < N; p++) {
        rc = LAPI_Setcntr(lapi_handle, &(ack_cntr[p].cntr), 0);
        if(rc) ERROR("armci_init_lapi: LAPI_Setcntr failed (ack)",rc);
        ack_cntr[p].val = 0;

        rc = LAPI_Setcntr(lapi_handle, &(get_cntr[p].cntr), 0);
        if(rc) ERROR("armci_init_lapi: LAPI_Setcntr failed (get)",rc);
        get_cntr[p].val = 0;
    }
    rc = LAPI_Setcntr(lapi_handle, &hdr_cntr.cntr, 0);
    if(rc) ERROR("armci_init_lapi: LAPI_Setcntr failed (hdr)",rc);
    hdr_cntr.val = 0;
    rc = LAPI_Setcntr(lapi_handle, &buf_cntr.cntr, 0);
    if(rc) ERROR("armci_init_lapi: LAPI_Setcntr failed (buf)",rc);
    buf_cntr.val = 0;
#if 0
    pcntr = (lapi_cmpl_t*)MessageSndBuffer;
    rc = LAPI_Setcntr(lapi_handle, &pcntr->cntr, 0);
    if(rc) ERROR("armci_init_lapi: LAPI_Setcntr failed (bufcntr)",rc);
    pcntr->val = 0;
#endif

#ifdef LAPI_RDMA
    /* allocate rCxt */
    lapi_remote_cxt = (lapi_user_cxt_t*)malloc(armci_nproc *
            sizeof(lapi_user_cxt_t));
    if(lapi_remote_cxt==NULL) ERROR("armci_init_lapi: rCxt malloc failed",0);

    /* obtain remote context "rCxt" for RDMA Operation of all procs */
    for(p = 0; p< armci_nproc; p++){
        if(p==armci_me) continue;
        util_cxt.Util_type   = LAPI_REMOTE_RCXT;
        util_cxt.operation   = LAPI_RDMA_ACQUIRE;
        util_cxt.dest        = p;
        CHECK(LAPI_Util(lapi_handle, (lapi_util_t *) &util_cxt));
        lapi_remote_cxt[p]   =  util_cxt.usr_rcxt;
    }
#endif

#if  !defined(LAPI2)

    /* for high performance, disable LAPI internal error checking */
    LAPI_Senv(lapi_handle, ERROR_CHK, 0);

#endif

    /* make sure that interrupt mode is on */
    LAPI_Senv(lapi_handle, INTERRUPT_SET, 1);

    /* initialize buffer managment module */
    _armci_buf_init();

#ifdef LAPI_RDMA
    CHECK((LAPI_Gfence(lapi_handle)));
#endif
#if ARMCI_ENABLE_GPC_CALLS
    gpc_req = (gpc_buf_t *)malloc(sizeof(gpc_buf_t)*MAX_GPC_REQ);
    if(gpc_req==NULL)armci_die("malloc for gpc failed",sizeof(gpc_buf_t));
    gpc_init();
#endif
}


void armci_term_lapi()
{
    int p;
    lapi_remote_cxt_t util_cxt;  /* For call to obtain rCxt */

#ifdef LAPI_RDMA 
    CHECK((LAPI_Gfence(lapi_handle)));

    /* release remote context "rCxt" for RDMA Operation of all procs */
    for(p = 0; p< armci_nproc; p++){
        if(p==armci_me) continue;
        util_cxt.Util_type   = LAPI_REMOTE_RCXT;
        util_cxt.operation   = LAPI_RDMA_RELEASE;
        util_cxt.dest        = p;
        util_cxt.usr_rcxt    = lapi_remote_cxt[p];
        CHECK(LAPI_Util(lapi_handle, (lapi_util_t *) &util_cxt));
    }
    free(lapi_remote_cxt);
#endif

#ifndef TCGMSG
    CHECK((LAPI_Term(lapi_handle))); /* terminate the LAPI handle */
#endif
    free(cmpl_arr);
    free(ack_cntr);
    free(get_cntr);
}

/* primitive pseudo message-passing on top of lapi */ 

/* send data to remote process using p specified message tag */
/* tag contains address of receive buffer guarded by cntr at process p */
void armci_lapi_send(msg_tag_t tag, void* data, int len, int p)
{
    int rc;
    lapi_cntr_t org_cntr;
    void *buf = tag.buf;
    lapi_cntr_t *cntr = tag.cntr;
    if(!buf)ERROR("armci_lapi_send: NULL tag(buf) error",0);
    if(!cntr)ERROR("armci_lapi_send:  NULL tag(cntr) error",0);

    rc=LAPI_Setcntr(lapi_handle, &org_cntr, 0);
    if(rc) ERROR("armci_lapi_send:setcntr failed",rc);
    rc=LAPI_Put(lapi_handle, (uint)p, (uint)len, buf, data, 
            cntr, &org_cntr, NULL);
    if(rc) ERROR("armci_lapi_send:put failed",rc);
    rc+=LAPI_Waitcntr(lapi_handle, &org_cntr, 1, NULL);
    if(rc) ERROR("armci_lapi_send:waitcntr failed",rc);
}

/* subroutine versions of macros disabling and enabling interrupts */
void intr_off_()
{
    INTR_OFF;
}

void intr_on_()
{
    INTR_ON;
}


void print_counters_()
{
    int i;
    printf("bufcntr: val =%d cntr=%d\n", buf_cntr.val, buf_cntr.cntr);
    for(i=0; i< armci_nproc;i++){
        printf("cmpl_arr: val=%d cntr=%d oper=%d\n",cmpl_arr[i].val,
                cmpl_arr[i].cntr, cmpl_arr[i].oper);
    }
    fflush(stdout);
}

#ifdef LAPI_RDMA 
/* LAPI Put RDMA */
void armci_client_direct_send(int p, void *src_buf, void *dst_buf,
        int len, void** contextptr, int nbtag,
        ARMCI_MEMHDL_T *lochdl,ARMCI_MEMHDL_T *remhdl) {

    lapi_xfer_t      xfer_struct;   /* Data structure for the xfer call */
    lapi_rdma_tag_t  lapi_rdma_tag; /* RDMA notification tag */
    uint src_offset, tgt_offset;
    int  val, rc;

    /* can be any number that fits in ushort */
    lapi_rdma_tag = 22;

    /* CHECK: offset problem. what if client and server attached (shmat) at
       diff address */ 
    src_offset = (char *)src_buf- (char *)lochdl->start;
    tgt_offset = (char *)dst_buf - (char *)remhdl->start;

#if DEBUG_
    printf("%d: Doing LAPI_Xfer (RDMA Put): dst=%d srchdl_start=%p remhdl_start=%p (bytes=%ld src_off=%d tgt_off=%d)\n", armci_me, p, lochdl->start, remhdl->start, len, src_offset, tgt_offset); fflush(stdout);
#endif

    bzero(&xfer_struct, sizeof(xfer_struct));
    xfer_struct.HwXfer.Xfer_type    = LAPI_RDMA_XFER;
    xfer_struct.HwXfer.tgt          = p;
    /*xfer_struct.HwXfer.op         = LAPI_RDMA_PUT|LAPI_RCNTR_UPDATE;*/
    xfer_struct.HwXfer.op           = LAPI_RDMA_PUT;
    xfer_struct.HwXfer.rdma_tag     = lapi_rdma_tag;
    xfer_struct.HwXfer.remote_cxt   = lapi_remote_cxt[p];
    xfer_struct.HwXfer.src_pvo      = lochdl->pvo;
    xfer_struct.HwXfer.tgt_pvo      = remhdl->pvo;
    xfer_struct.HwXfer.src_offset   = src_offset;
    xfer_struct.HwXfer.tgt_offset   = tgt_offset;
    xfer_struct.HwXfer.len          = (ulong) (len);  
    xfer_struct.HwXfer.shdlr        = (scompl_hndlr_t *) NULL;
    xfer_struct.HwXfer.sinfo        = (void *) NULL;
    xfer_struct.HwXfer.org_cntr     = &(ack_cntr->cntr);

    /* Initiate RDMA Xfer */
    if((rc = LAPI_Xfer(lapi_handle, &xfer_struct)) != LAPI_SUCCESS) {
        LAPI_Msg_string(rc, err_msg_buf);
        fprintf(stderr, "LAPI ERROR: %s, rc = %d\n", err_msg_buf, rc);
        armci_die("LAPI_Xfer (RDMA Put) failed", 0);
    }

    /* wait for RDMA completion */
    rc = LAPI_Waitcntr(lapi_handle, &(ack_cntr->cntr),1,&val);
    if(rc != LAPI_SUCCESS) {
        LAPI_Msg_string(rc, err_msg_buf);
        fprintf(stderr, "LAPI ERROR: %s, rc = %d\n", err_msg_buf, rc);
        armci_die("LAPI_Waitcntr (RDMA Put) failed", 0);
    }

    /* CHECK((LAPI_Fence(lapi_handle))); */

#if DEBUG_
    printf("%d: Completed LAPI_Xfer RDMA (Put): dst=%d\n", armci_me, p);
#endif
}

/* LAPI Get RDMA */
void armci_client_direct_get(int p, void *src_buf, void *dst_buf, 
        int len, void** cptr, int nbtag,
        ARMCI_MEMHDL_T *lochdl, ARMCI_MEMHDL_T *remhdl) {

    lapi_xfer_t      xfer_struct;   /* Data structure for the xfer call */
    lapi_rdma_tag_t  lapi_rdma_tag; /* RDMA notification tag */
    uint src_offset, tgt_offset;
    int  val, rc;

    /* can be any number that fits in ushort */
    lapi_rdma_tag = 21;

    /* CHECK: offset problem. what if client and server attached (shmat) at
       diff address */ 
    src_offset = (char *)dst_buf- (char *)lochdl->start;
    tgt_offset = (char *)src_buf - (char *)remhdl->start;

#if DEBUG_
    printf("%d: Doing LAPI_Xfer (RDMA Get): dst=%d srchdl_start=%p remhdl_start=%p (bytes=%ld src_off=%d tgt_off=%d)\n", armci_me, p, lochdl->start, remhdl->start, len, src_offset, tgt_offset); fflush(stdout);
#endif
    bzero(&xfer_struct, sizeof(xfer_struct));
    xfer_struct.HwXfer.Xfer_type    = LAPI_RDMA_XFER;
    xfer_struct.HwXfer.tgt          = p;
    /*xfer_struct.HwXfer.op         = LAPI_RDMA_GET|LAPI_RCNTR_UPDATE;*/
    xfer_struct.HwXfer.op           = LAPI_RDMA_GET;
    xfer_struct.HwXfer.rdma_tag     = lapi_rdma_tag;
    xfer_struct.HwXfer.remote_cxt   = lapi_remote_cxt[p];
    xfer_struct.HwXfer.src_pvo      = lochdl->pvo;
    xfer_struct.HwXfer.tgt_pvo      = remhdl->pvo;
    xfer_struct.HwXfer.src_offset   = src_offset;
    xfer_struct.HwXfer.tgt_offset   = tgt_offset;
    xfer_struct.HwXfer.len          = (ulong) (len);  
    xfer_struct.HwXfer.shdlr        = (scompl_hndlr_t *) NULL;
    xfer_struct.HwXfer.sinfo        = (void *) NULL;
    xfer_struct.HwXfer.org_cntr     = &(get_cntr->cntr);

    /* Initiate RDMA Xfer */
    if((rc = LAPI_Xfer(lapi_handle, &xfer_struct)) != LAPI_SUCCESS) {
        LAPI_Msg_string(rc, err_msg_buf);
        fprintf(stderr, "LAPI ERROR: %s, rc = %d\n", err_msg_buf, rc);
        armci_die("LAPI_Xfer (RDMA Get) failed", 0);
    }

    /* wait for RDMA completion */
    rc = LAPI_Waitcntr(lapi_handle, &(get_cntr->cntr),1,&val);
    if(rc != LAPI_SUCCESS) {
        LAPI_Msg_string(rc, err_msg_buf);
        fprintf(stderr, "LAPI ERROR: %s, rc = %d\n", err_msg_buf, rc);
        armci_die("LAPI_Waitcntr (RDMA Get) failed", 0);
    }

#if DEBUG_
    printf("%d: Completed LAPI_Xfer (RDMA Get): dst=%d\n", armci_me, p);
#endif

}

int armci_pin_contig_hndl(void *ptr, int bytes, ARMCI_MEMHDL_T *memhdl)
{

    lapi_get_pvo_t util_pvo;     /* For call to obtain PVO */
    int rc;

    /* translate and pin the buffer to the adapter */
    util_pvo.Util_type = LAPI_XLATE_ADDRESS;
    util_pvo.length    = bytes;
    util_pvo.usr_pvo   = 0;
    util_pvo.address   = ptr;
    util_pvo.operation = LAPI_RDMA_ACQUIRE;
    /*bzero(ptr, bytes);*/ /* CHECK: Is touching the entire shmem sgement feasible */
    if((rc=LAPI_Util(lapi_handle, (lapi_util_t *) &util_pvo)) != LAPI_SUCCESS) {
        return 0;
    }

    memhdl->pvo   = util_pvo.usr_pvo;
    memhdl->start = ptr;

#if DEBUG_
    printf("\n%d:armci_pin_contig_hndl(): memhdl(pvo)=%ld ptr=%p bytes=%ld\n",
            armci_me, (long)memhdl->pvo, ptr, bytes);fflush(stdout);
#endif

    return 1;
}

void armci_network_client_deregister_memory(ARMCI_MEMHDL_T *mh)
{
    fprintf(stderr," armci_network_client_deregister_memory not available\n");
    fflush(stderr);
}


void armci_network_server_deregister_memory(ARMCI_MEMHDL_T *mh)
{
    fprintf(stderr," armci_network_client_deregister_memory not available\n");
    fflush(stderr);
}
#endif /* LAPI_RDMA */

#ifdef AIX

#define LOCKED 1
void armci_lapi_lock(int *lock)
{
    atomic_p word_addr = (atomic_p)lock;
    int spin = 1;


    while(1){

        if(_check_lock(word_addr, 0, LOCKED) == FALSE )
            break; /* we got the lock */ 

        if(spin){
            armci_waitsome(1);
            spin = 0;
        }else{

            /* yield processor to another thread */
            /* cannot yield w/o affecting thread priority - better sleep */
            /* yield(); */

            /* call usleep to notify scheduler */
            (void)usleep(5);
        }
    }
}


void armci_lapi_unlock(int *lock)
{
    atomic_p word_addr = (atomic_p)lock;

    if(_check_lock(word_addr, LOCKED, 0) == TRUE ) 
        armci_die("somebody else unlocked",0);
}
#endif

#ifdef LAPI2
#include "lapi2.c"
#endif

