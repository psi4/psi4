#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "armcip.h"
#include "request.h"
#include "message.h"
#include "memlock.h"
#include "copy.h"
#include "gpc.h"
#include "iterator.h"
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif
#if HAVE_PROCESS_H
#   include <process.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif

#define DEBUG_ 0
#define DEBUG1 0

#ifndef SERV
#     define SERV 2
#endif

#ifdef SOCKETS
#   define EQ_TAGS(a_, b_) ((a_) == (b_))
#else
#   define EQ_TAGS(a_, b_) !memcmp(&(a_), &(b_), sizeof(a_))
#endif

int _armci_server_started=0;

#if defined(SOCKETS)
extern active_socks_t *_armci_active_socks;
extern void armci_sock_send(int to, void *data, int len);
#endif

/**************************** pipelining for medium size msg ***********/
#ifdef PIPE_BUFSIZE


static int pack_size(int len)
{
int oldlen = len;
#define PIPE_ROUNDUP 512 
#define PIPE_SHORT_ROUNDUP (1024) 
int n;
 if(len <4*PIPE_BUFSIZE){  
   len /=2;
    n = len%PIPE_SHORT_ROUNDUP;
   if(n)len += (PIPE_SHORT_ROUNDUP-n);
 } 
#if defined(VIA) || defined(VAPI)
 else if(len <25*PIPE_BUFSIZE){
   len /=4;
   n = len%PIPE_SHORT_ROUNDUP;
   if(n)len += (PIPE_SHORT_ROUNDUP-n);
 }
else if(len <41*PIPE_BUFSIZE){
   len /=8;
   n = len%PIPE_SHORT_ROUNDUP;
   if(n)len += (PIPE_SHORT_ROUNDUP-n);
 }
#else
 else if(len <32*PIPE_BUFSIZE){
   len /=8;
   n = len%PIPE_SHORT_ROUNDUP;
   if(n)len += (PIPE_SHORT_ROUNDUP-n);
 }
#endif
else
#if defined(VIA) || defined(VAPI)
   len = 8*4096;
#elif defined(HITACHI)
   len = 128*1024-128;
#else
   len = 64*1024-128;
#endif
#ifdef MAX_PIPELINE_CHUNKS

 if(oldlen/len > MAX_PIPELINE_CHUNKS-1){
  len = oldlen/MAX_PIPELINE_CHUNKS;
  n = len%PIPE_SHORT_ROUNDUP;
  if(n)len += (PIPE_SHORT_ROUNDUP-n);
 }
#endif
 return len;
}

#define PACK_SIZE1(_len) ((_len)<PIPE_BUFSIZE)?PIPE_MIN_BUFSIZE:PIPE_BUFSIZE;
#define PACK_SIZE(_len) pack_size(_len)

void armci_pipe_prep_receive_strided(request_header_t *msginfo, char *buf,
                        int strides, int stride_arr[], int count[], int bufsize)
{
buf_arg_t arg;
int  packsize = PACK_SIZE(msginfo->datalen);

     arg.buf_posted = arg.buf   = buf;
#ifdef HITACHI
     arg.count = 0;
#else
     arg.count = bufsize;
#endif
     arg.proc  = (msginfo->operation==GET)?msginfo->to:msginfo->from;
     arg.op    = msginfo->operation;

     armci_dispatch_strided(buf, stride_arr, count, strides, -1, -1,
                            packsize, armcill_pipe_post_bufs,&arg);
}

void armci_pipe_receive_strided(request_header_t* msginfo, void *ptr,
                                int stride_arr[], int count[], int strides)
{
buf_arg_t arg;
int  packsize = PACK_SIZE(msginfo->datalen);
#if defined(GM)
     arg.buf_posted   = msginfo->tag.data_ptr;
#endif
#if (defined(VIA) && defined(VIA_USES_RDMA)) || defined(VAPI)
     arg.buf_posted   = msginfo->tag;
#endif

     arg.buf   = ptr;
     arg.count = 0;
     arg.proc  = (msginfo->operation==GET)?msginfo->to:msginfo->from;
     arg.op    = msginfo->operation;

     armci_dispatch_strided(ptr, stride_arr, count, strides, -1, -1,
                            packsize, armcill_pipe_extract_data, &arg);
}

void armci_pipe_send_strided(request_header_t *msginfo, void *buf, int buflen,
                             void *ptr, int *stride_arr,int count[],int strides)
{
buf_arg_t arg;
int  packsize = PACK_SIZE(msginfo->datalen);

#if defined(GM) || defined(HITACHI)
     arg.buf_posted   = msginfo->tag.data_ptr;
#endif
#if (defined(VIA) && defined(VIA_USES_RDMA)) || defined(VAPI)
     arg.buf_posted   = msginfo->tag;
#endif

     arg.buf   = buf;
     arg.count = 0;
     arg.proc  = (msginfo->operation==GET)?msginfo->from:msginfo->to;
     arg.op    = msginfo->operation;

     armci_dispatch_strided(ptr, stride_arr, count, strides, -1, -1,
                            packsize, armcill_pipe_send_chunk, &arg);
#ifdef GM
     armci_serv_send_nonblocking_complete(0);
#endif
}
#endif
/**************************** end of pipelining for medium size msg ***********/


#if defined(CLIENT_BUF_BYPASS) && !defined(GM)
/**************** NOTE: for now this code can only handle contiguous data *****/
void armci_send_strided_data_bypass(int proc, request_header_t *msginfo,
                                    void *loc_buf, int msg_buflen,
                                    void *loc_ptr, int *loc_stride_arr,
                                    void *rem_ptr, int *rem_stride_arr,
                                    int *count, int stride_levels)
{
    int armcill_server_wait_ack(int,int);
    if(DEBUG_){
      printf("%d(s): strided(%d) get bypass from %d\n",armci_me,stride_levels,
             msginfo->from);
      fflush(stdout);
    }

#ifdef VAPI
    if(stride_levels==0 && msginfo->pinned){
      armci_send_contig_bypass(proc,msginfo,loc_ptr,rem_ptr,count[0]);
      return;
    }
    else {
      armci_die("***Contact Developers with machine/network info at hpctools@emsl.pnl.gov: bypass path wrongly invoked***",0);
    }
#endif

    armci_pin_memory(loc_ptr, loc_stride_arr,count, stride_levels);
    /*wait until client ready*/
    if(!armcill_server_wait_ack(msginfo->from,1)){
       /*client was not able to pin memory, it will revert to default protocol
         hence, unpin the memory and leave. */
       armci_unpin_memory(loc_ptr, loc_stride_arr,count, stride_levels);
       return;
    }

    armcill_server_put(msginfo->from,loc_ptr,rem_ptr,count[0]);
    armci_unpin_memory(loc_ptr, loc_stride_arr,count, stride_levels);
    if(DEBUG_){
      printf("%d(s): strided(%d) get bypass done \n",armci_me,stride_levels);
      fflush(stdout);
    }
}
#endif

/*\ client initialization
\*/
void armci_client_code()
{
   if(DEBUG_){
      printf("in client after fork %d(%d)\n",armci_me,getpid()); fflush(stdout);
   }

   armci_client_connect_to_servers();
   armci_msg_barrier();

   if(DEBUG_){
      printf("%d client connected to all %d servers\n",armci_me, armci_nclus-1);
      fflush(stdout);
   }
}


/*\ client sends request to server
\*/
void armci_send_req(int proc, request_header_t* msginfo, int len)
{
int hdrlen = sizeof(request_header_t);
int bytes;

    if(msginfo->operation == GET) {
        if(msginfo->format==VECTOR && msginfo->ehlen > 0)
            bytes = msginfo->dscrlen + hdrlen + msginfo->datalen;
        else
            bytes = msginfo->dscrlen + hdrlen;
    } else bytes = msginfo->bytes + hdrlen;

    if(DEBUG_){printf("%d: sending req %d (len=%d dscr=%d data=%d) to %d \n",
               armci_me, msginfo->operation, bytes,msginfo->dscrlen,
               msginfo->datalen,proc); fflush(stdout);
    }
    if(bytes > len)armci_die2("armci_send_req:buffer overflow",bytes,len);

#ifdef PIPE_BUFSIZE
    if(
#   ifdef CLIENT_BUF_BYPASS 
     (!msginfo->bypass) &&
#   endif
     (msginfo->datalen>2*PIPE_MIN_BUFSIZE) && (msginfo->operation == GET)
                                        && (msginfo->format == STRIDED)){
      char *buf = sizeof(void*) + (char*)(msginfo+1);
      int  *ibuf = (int*)buf;
      int  *strides =ibuf;
      int  *stride_arr= ibuf +1;
      int  *count = stride_arr + *strides;
      armci_pipe_prep_receive_strided(msginfo, buf, *strides, stride_arr, count,
                                      len-2**strides*sizeof(int)-sizeof(void*));
      armci_pipe_send_req(proc,msginfo, bytes);
    }else
#endif

       armci_send_req_msg(proc,msginfo, bytes);
}


/*\ client sends strided data + request to server
\*/
void armci_send_strided(int proc, request_header_t *msginfo, char *bdata,
                        void *ptr, int strides, int stride_arr[], int count[])
{
    int hdrlen = sizeof(request_header_t);
    int dscrlen = msginfo->dscrlen;
    int datalen = msginfo->datalen;
    int cluster = armci_clus_id(proc);
    int bytes;

    bytes = msginfo->bytes + hdrlen;

    if(DEBUG_){
       printf("%d:sending strided %d to(%d,%d,%d) bytes=%d dslen=%d dlen=%d,\n",
                armci_me, msginfo->operation, msginfo->to,
                cluster, proc, bytes, dscrlen, datalen); fflush(stdout);
    }

#if defined(SOCKETS)
    /* zero-copy optimization for large requests */
    if(count[0] >  TCP_PAYLOAD){
       if(armci_send_req_msg_strided(proc, msginfo,ptr,strides,
          stride_arr, count))armci_die("armci_send_strided_req long: failed",0);
       return; /************** done **************/
    }
#elif defined(MPI_SPAWN_ZEROCOPY)
    /* zero-copy optimization for large requests */
    if(msginfo->operation==PUT && msginfo->datalen==0 && count[0]>TCP_PAYLOAD){
       if(armci_send_req_msg_strided(proc, msginfo,ptr,strides,
          stride_arr, count))armci_die("armci_send_strided_req long: failed",0);
       return; /************** done **************/
    }
#elif defined(PIPE_BUFSIZE___)
#warning Network resource is only locked inside armci_send_req_msg, no common lock
    if((msginfo->datalen>2*PIPE_MIN_BUFSIZE) && (msginfo->operation == PUT)){
       msginfo->bytes =0; /*** this tells server that we use pipelined put ****/
       armci_send_req_msg(proc,msginfo, hdrlen+dscrlen);
       armci_pipe_send_strided(msginfo, bdata, datalen,
                               ptr, stride_arr, count, strides);
       return; /************** done **************/
    }
#endif
    /*  copy into a buffer before sending */

#  ifdef SERV_BUF_IDX_T
    msginfo->inbuf = armcill_getbidx((msginfo->datalen+msginfo->dscrlen), proc, &msginfo->tag.ack);
    msginfo->tag.ack_ptr = &msginfo->tag.ack;
#  endif
    armci_write_strided(ptr, strides, stride_arr, count, bdata);
    if(armci_send_req_msg(proc,msginfo, bytes))
       armci_die("armci_send_strided_req: failed",0);
}

#ifdef SOCKETS
/* main handler to process responses from dataserver */
void armci_rcv_hdlr(request_header_t* msginfo)
{
    int not_rcvd, my_id, nready, i, n, rc;
    msg_tag_t rcvd_id;
    BUF_INFO_T *info;

    n = MAX_BUFS + MAX_SMALL_BUFS;
    my_id = _armci_buf_to_bufinfo(msginfo)->bufid;
    not_rcvd = 1;

    while (not_rcvd) {
        THREAD_LOCK(armci_user_threads.net_lock);

        if (!_armci_buf_cmpld(my_id)) {
            /* my buffer has not been completed yet */
            nready=armci_WaitSock(_armci_active_socks->socks,n,_armci_active_socks->ready);
            if (nready) {
                for (i = 0; i < n && nready; i++) {
                    if (!_armci_active_socks->ready[i]) continue; /* not a ready sock */
                    nready--;
                    /* receive data from socket _armci_active_socks->ready[i]
                     * Note: socks[i] is the socket with incoming data HOWEVER
                     * i is NOT necessarily the index of the associated buffer.
                     * This is because armci_WaitSock will mark ALL entries for
                     * the same socket as ready
                     */
                    rc = armci_ReadFromSocket(_armci_active_socks->socks[i],
                                            &rcvd_id, sizeof(rcvd_id));
                    if(rc<0)armci_die("armci_rcv_strided_data: read tag failed",rc);


                    /* receive response and process it */
                    msginfo = (request_header_t *)_armci_buf_ptr_from_id(rcvd_id);
                    switch (msginfo->operation) {
                        case PUT:
                        case UNLOCK:
                            armci_die("armci_rcv_hdlr: unexpected op",msginfo->operation);
                            break;

                        case GET:
                            info = _armci_id_to_bufinfo(rcvd_id);
                            armci_complete_req_buf(info, msginfo);
                            break;

                        case LOCK:
                        case RMW:
                        case ARMCI_SWAP:
                        case ARMCI_SWAP_LONG:
                        case ARMCI_FETCH_AND_ADD:
                        case ARMCI_FETCH_AND_ADD_LONG:
                            armci_rcv_data(msginfo->to, msginfo);
                            break;

                        case ACK:
                            armci_rcv_data(NODE_SERVER(msginfo->to), msginfo);
                            break;

                        default:
                            armci_die("armci_rcv_hdlr: unrecognized op",msginfo->operation);
                    }

                    /* clear this socket in active sockets */
                    _armci_active_socks->socks[i] = -1;

                    THREAD_UNLOCK(armci_user_threads.net_lock);

                    /* check if the data we received were sent to us */
                    if (rcvd_id == my_id) not_rcvd = 0;
                    /* mark received buffer as completed */
                    _armci_buf_set_cmpld_idx(rcvd_id, 1);
                }
                if(nready)armci_die("armci_rcv_hdlr:nready in not consistent",nready);
            } else { /* timed out in select */
                THREAD_UNLOCK(armci_user_threads.net_lock);
                cpu_yield();
            }
        } else { /* buffer was completed by another thread */
            THREAD_UNLOCK(armci_user_threads.net_lock);
            not_rcvd = 0;
        }
    }
}


#if 0
/* receives plain(contiguous) data from dataserver */
void armci_rcv_data_hdlr(int bufid)
{
    request_header_t* msginfo;
    int proc, datalen;
    char *buf;

    /* obtain buffer and buffer info associated with this receive */
    msginfo = (request_header_t *)_armci_buf_ptr_from_id(bufid);
    proc = msginfo->to;
    datalen = msginfo->datalen;

    if(datalen == 0) armci_die("armci_rcv_data_hdlr: no data to receive",datalen);
    if(datalen > (MSG_BUFLEN-sizeof(request_header_t)-sizeof(long)))
        armci_die("armci_rcv_data_hdlr: data overflowing rcv buffer",datalen);

    /* fills msginfo buffer */
    buf = armci_ReadFromDirect(proc, msginfo, datalen);

    if ((char *)(msginfo+1) != buf)
        armci_die("armci_rcv_data_hdlr: buf != msginfo+1",datalen);
}

/* received strided data from dataserver */
void armci_rcv_strided_data_hdlr(int bufid)
{
    BUF_INFO_T *info;
    char *dscr;
    request_header_t *msginfo;
    void *ptr;
    int proc, strides, *stride_arr, *count;
    char *databuf;

    /* obtain buffer and buffer info associated with this receive */
    info = _armci_id_to_bufinfo(bufid);
    dscr = info->dscr;
    msginfo = (request_header_t *)_armci_buf_ptr_from_id(bufid);
    proc = msginfo->to;

    /* ptr, strides, stride_arr and count should be extracted from buf_info */
    ptr = *(void**)dscr;       dscr += sizeof(void*);
    strides = *(int*)dscr;     dscr += sizeof(int);
    stride_arr = (int*)dscr;   dscr += strides*sizeof(int);
    count = (int*)dscr;

    /* actual rcv: copied from old armci_rcv_strided_data */
    /* zero-copy optimization for large requests */
    if(count[0] >  TCP_PAYLOAD){
       armci_ReadStridedFromDirect(proc,msginfo,ptr,strides,stride_arr, count);
       return; /*********************** done ************************/
    }

    databuf = armci_ReadFromDirect(proc,msginfo,msginfo->datalen);
    armci_read_strided(ptr, strides, stride_arr, count, databuf);
}

/* receives vector data from dataserver */
void armci_rcv_vector_data_hdlr(int bufid)
{
    request_header_t* msginfo;
    int proc, datalen;
    char *buf;
    armci_giov_t *darr;

    /* obtain buffer and buffer info associated with this receive */
    msginfo = (request_header_t *)_armci_buf_ptr_from_id(bufid);
    proc = msginfo->to;
    /*datalen = msginfo->datalen;*/
    buf = (char *)(msginfo + 1);

    /* receive vector as cont block, data is in buf */
    armci_rcv_data_hdlr(bufid);

    /* unpack vector */
    /* armci_giov_t darr[], int len) */
     armci_vector_from_buf(darr, len, buf);
}
#endif
#endif

/*\ client receives data from server
\*/
char *armci_rcv_data(int proc, request_header_t* msginfo)
{
int datalen = msginfo->datalen;
char *buf;
    if(DEBUG_) {
        printf("%d:armci_rcv_data:  bytes= %d \n", armci_me, datalen);
        fflush(stdout);
    }

    if(datalen == 0) armci_die("armci_rcv_data: no data to receive",datalen);
    if(datalen > (((int)MSG_BUFLEN)-((int)sizeof(request_header_t))-((int)sizeof(long))))
        armci_die("armci_rcv_data:data overflowing rcv buffer",datalen);

    buf = armci_ReadFromDirect(proc, msginfo, datalen);

    if(DEBUG_){
        printf("%d:armci_rcv_data: got %d bytes \n",armci_me,datalen);
        fflush(stdout);
    }
    return(buf);
}

/*\ client receives vector data from server and unpacks to the right loc
\*/
void armci_rcv_vector_data(int proc, request_header_t* msginfo, armci_giov_t darr[], int len)
{
    char *buf = armci_rcv_data(proc, msginfo);
    armci_vector_from_buf(darr, len, buf);
}

/*\ client receives strided data from server
\*/
#if 0
void armci_rcv_strided_data(int proc, request_header_t* msginfo, int datalen, 
                            void *ptr, int strides,int stride_arr[],int count[])
{
extern BUF_INFO_T *_armci_tag_to_bufinfo(msg_tag_t tag);
extern BUF_INFO_T *_armci_buf_to_bufinfo(void *buf);
    int not_received = 1;
    int sel, idx, rc, n=MAX_BUFS+MAX_SMALL_BUFS;
    char *databuf;
    msg_tag_t tag;
    BUF_INFO_T *info;
    char *dscr;
    void *loc_ptr;
    int stride_levels, *loc_stride_arr;
    request_header_t* buf;

    while (not_received) {
         THREAD_LOCK(armci_user_threads.net_lock);

         if (!_armci_buf_cmpld(msginfo)) {
             /* buffer not completed */
#ifdef SOCKETS
             sel=armci_WaitSock(_armci_active_socks->socks,n,_armci_active_socks->ready);
#endif
             if (sel > 0) {
                 /* pick a socket (should I check if sel > 1?) */
                 for(idx=0;idx<n;idx++)if(_armci_active_socks->ready[idx])break;

                 /* socks[idx] is the socket with incoming data HOWEVER idx is
                  * NOT necessarily the index of the associated buffer. It is
                  * because armci_WaitSock will mark ALL entries for the same
                  * socket as ready */

                 /* read tag */
#ifdef SOCKETS
                 rc=armci_ReadFromSocket(_armci_active_socks->socks[idx],&tag,sizeof(tag));
                 if(rc<0)armci_die("armci_rcv_strided_data: read tag failed",rc);
#if 0 || defined(DTAG)
                 idx = tag & DTAG;
                 printf("DAG RCV: dtag=%ld,idx=%d,",tag,idx);
                 tag >>= (sizeof(msg_id_t) * 8);
                 printf("tag=%d,",tag);

                 /* find proper buffer idx */
                 info = _armci_tag_to_bufinfo(tag);
                 printf("idx(tag)=%d\n",info->bufid); fflush(stdout);
                 if(info->bufid!=idx)armci_die("armci_rcv_strided_data: bad tag",tag);
#else
                 idx = tag;
#endif
#endif
                 info = _armci_id_to_bufinfo(idx);
                 dscr = info->dscr;

                 /* network complete -- old armci_rcv_strided_data
                  * ptr, strides, stride_arr and count should be extracted from buf_info */
                 ptr = *(void**)dscr;       dscr += sizeof(void*);
                 strides = *(int*)dscr;     dscr += sizeof(int);
                 stride_arr = (int*)dscr;   dscr += strides*sizeof(int);
                 count = (int*)dscr;

                 /* find appropriate msginfo for received response */
                 buf = (request_header_t *)_armci_buf_ptr_from_id(idx);
                 proc = buf->to;

#ifdef CLIENT_BUF_BYPASS
#error THIS PATH IS NOT UPDATED
                 if(msginfo->bypass){
                     /* zero-copy protocol: get ACK and then unpin user buffer */
                     armci_rcv_strided_data_bypass(proc, msginfo, ptr, strides);
                     armci_unpin_memory(ptr, stride_arr, count, strides);
                     return; /* we are done */
                 }
#endif

#ifdef SOCKETS
                 /* zero-copy optimization for large requests */
                 if (count[0] > TCP_PAYLOAD)
                    armci_ReadStridedFromDirect(proc,buf,ptr,strides,stride_arr,count);
                 else
#elif defined(PIPE_BUFSIZE)
                 if (buf->datalen > 2*PIPE_MIN_BUFSIZE)
                    armci_pipe_receive_strided(buf, ptr, stride_arr, count, strides);
                 else
#endif
                 {
                    databuf = armci_ReadFromDirect(proc,buf,datalen);
                    armci_read_strided(ptr, strides, stride_arr, count, databuf);
                 }

                 /* update active sockets */
                 _armci_active_socks->socks[idx] = -1;

                 THREAD_UNLOCK(armci_user_threads.net_lock);

                 /* check if the response we received was for our request */
                 if (EQ_TAGS(tag, _armci_buf_to_bufinfo(msginfo)->tag)) {
/*                 if (tag == _armci_buf_to_bufinfo(buf)->tag) {
                     _armci_buf_release(msginfo); released in armci_rem_strided */
                     not_received = 0;
                 } else {
                     _armci_buf_set_cmpld_idx(idx, 1); /* completed */
                 }
             } else {
                 THREAD_UNLOCK(armci_user_threads.net_lock);
#if 0
                 if(sel)armci_die("armci_rcv_strided_data: error in select",errno);
                 else
#endif
                     cpu_yield(); /* timed out in select */
             }
         } else { /* buffer was completed by another thread */
             THREAD_UNLOCK(armci_user_threads.net_lock);

/*             _armci_buf_release(msginfo); released in armci_rem_strided */
             not_received = 0;
         }
    }
}
#else
void armci_rcv_strided_data(int proc, request_header_t* msginfo, int datalen,
                            void *ptr, int strides,int stride_arr[],int count[])
{
    char *databuf;

    if(DEBUG_){
        printf("%d: armci_rcv_strided_data: expecting datalen %d from %d\n",
                armci_me, datalen, proc); fflush(stdout);
    }

#ifdef CLIENT_BUF_BYPASS
    if(msginfo->bypass){
       /* zero-copy protocol: get ACK and then unpin user buffer */
       armci_rcv_strided_data_bypass(proc, msginfo, ptr, strides);
       armci_unpin_memory(ptr, stride_arr, count, strides);
       return; /* we are done */
    }
#endif

#if defined(SOCKETS) || defined(MPI_SPAWN_ZEROCOPY)
    /* zero-copy optimization for large requests */
    if(count[0] >  TCP_PAYLOAD){
       armci_ReadStridedFromDirect(proc,msginfo,ptr,strides,stride_arr, count);
       return; /*********************** done ************************/
    }
#elif defined(PIPE_BUFSIZE)
    if(msginfo->datalen>2*PIPE_MIN_BUFSIZE){
       armci_pipe_receive_strided(msginfo, ptr, stride_arr, count, strides);
       return; /*********************** done ************************/
    }
#endif

#if !defined(GET_STRIDED_COPY_PIPELINED)
    databuf = armci_ReadFromDirect(proc,msginfo,datalen);
    armci_read_strided(ptr, strides, stride_arr, count, databuf);
#else
    {
      int bytes_buf = 0, bytes_usr = 0, seg_off=0;
      int ctr=0;
      stride_info_t sinfo;
      char *armci_ReadFromDirectSegment(int proc,request_header_t *msginfo,
					int datalen, int *bytes_buf);

      armci_stride_info_init(&sinfo,ptr,strides,stride_arr,count);
      do {
	databuf = armci_ReadFromDirectSegment(proc,msginfo,datalen,&bytes_buf);
	bytes_usr += armci_read_strided_inc(&sinfo,&databuf[bytes_usr],bytes_buf-bytes_usr, &seg_off);
      } while(bytes_buf<datalen);
      dassert(1,bytes_buf == bytes_usr);
      armci_stride_info_destroy(&sinfo);
    }
#endif
}
#endif



/*\ get ACK from server
\*/
void armci_rem_ack(int clus)
{
    int bufsize = sizeof(request_header_t)+sizeof(int);
    int destproc = 0;
    request_header_t *msginfo;
    destproc = SERVER_NODE(clus);
    msginfo = (request_header_t *)GET_SEND_BUFFER(bufsize,ACK,destproc);

    bzero(msginfo, sizeof(request_header_t));
    msginfo->dscrlen = 0;
    msginfo->from  = armci_me;
    msginfo->to    = SERVER_NODE(clus);
    msginfo->operation = ACK;
    msginfo->bytes   =0;
    msginfo->datalen =sizeof(int);
#ifdef SOCKETS
    msginfo->tag = BUF_TO_BUFINFO(msginfo)->bufid;
#endif

    if(DEBUG_){
        printf("%d(c):sending ACKreq to %d clus=%d\n",armci_me,msginfo->to,clus);
        fflush(stdout);
    }

    armci_send_req(armci_clus_info[clus].master, msginfo, bufsize);
#ifdef SOCKETS
    armci_rcv_hdlr(msginfo);
#else
    armci_rcv_data(armci_clus_info[clus].master, msginfo);  /* receive ACK */
#endif
    assert(*(int*)(msginfo+1) == ACK);
#ifdef VAPI
    assert(*(((int *)(msginfo+1))+1) == ARMCI_STAMP);
#endif
    FREE_SEND_BUFFER(msginfo);
}



/***************************** server side *********************************/

static void armci_check_req(request_header_t *msginfo, int buflen)
{

    if((msginfo->to != armci_me && msginfo->to < armci_master) ||
       msginfo->to >= armci_master + armci_clus_info[armci_clus_me].nslave)
        armci_die("armci_rcv_req: invalid to", msginfo->to);
#if 0
    /* should be done in recv_req */
    if(msginfo->operation != GET && msginfo->bytes > buflen)
        armci_die2("armci_rcv_req: message overflowing rcv buffer",
                  msginfo->bytes,MSG_BUFLEN);
#endif

    if(msginfo->dscrlen < 0)
        armci_die("armci_rcv_req: dscrlen < 0", msginfo->dscrlen);
    if(msginfo->datalen < 0)
        armci_die("armci_rcv_req: datalen < 0", msginfo->datalen);
#ifndef PIPE_BUFSIZE
    if(msginfo->dscrlen > (int)msginfo->bytes)
        armci_die2("armci_rcv_req: dsclen > bytes", msginfo->dscrlen,
                   msginfo->bytes);
#endif
}


/*\ server response - send data to client
\*/
void armci_send_data(request_header_t* msginfo, void *data)
{
    int to = msginfo->from;

#if defined(VIA) || defined(GM) || defined(VAPI)
    /* if the data is in the pinned buffer: MessageRcvBuffer */
#if defined(PEND_BUFS)
    extern int armci_data_in_serv_buf(void *);
    if(armci_data_in_serv_buf(data))
#else
    if((data > (void *)MessageRcvBuffer) &&
       (data < (void *)(MessageRcvBuffer + MSG_BUFLEN)))
#endif
        /* write the message to the client */
        armci_WriteToDirect(to, msginfo, data);
    else {
        /* copy the data to the MessageRcvBuffer */
#ifdef GM
        /* leave space for header ack */
        char *buf = MessageRcvBuffer + sizeof(long);
#else
        char *buf = MessageRcvBuffer;
# if defined(PEND_BUFS)
	fprintf(stderr, "%d:: op=%d len=%d ptr=%p working on unpinned memory. aborting!\n", armci_me, msginfo->operation,msginfo->datalen, data);
	assert(0);
	buf = NULL;
/*         extern char *armci_openib_get_msg_rcv_buf(int); */
/* 	buf = armci_openib_get_msg_rcv_buf(msginfo->from); */
# endif
#endif
	assert(buf != NULL);
        armci_copy(data, buf, msginfo->datalen);
        armci_WriteToDirect(to, msginfo, buf);
    }
#else
#ifdef DOELAN4
        /*this is because WriteToDirect is a no-op in elan4.c so we have
         * to do a put. This will not cause problems anywhere else in the
         * code and this part on elan4 will only be invoked in a GPC
         */
        PARMCI_Put(data,msginfo->tag.data_ptr,msginfo->datalen,to);
#else
        armci_WriteToDirect(to, msginfo, data);
#endif
#endif
}


/*\ server sends strided data back to client
\*/
void armci_send_strided_data(int proc,  request_header_t *msginfo,
                             char *bdata, void *ptr, int strides,
                             int stride_arr[], int count[])
{

    int to = msginfo->from;

    if(DEBUG_){ printf("%d(server): sending datalen = %d to %d %p\n",
                armci_me, msginfo->datalen, to,ptr); fflush(stdout); }
 
#if defined(SOCKETS) || defined(MPI_SPAWN_ZEROCOPY)
    /* zero-copy optimization for large requests */
    if(count[0] >  TCP_PAYLOAD){
       armci_WriteStridedToDirect(to,msginfo,ptr, strides, stride_arr, count);
       return; /*********************** done ************************/
    }
#elif defined(PIPE_BUFSIZE)
    if(msginfo->datalen>2*PIPE_MIN_BUFSIZE) {
        armci_pipe_send_strided(msginfo, bdata, msginfo->datalen,
                                ptr, stride_arr, count, strides);
        return;
    }
#endif

#if defined(GET_NO_SRV_COPY)
    {
      ARMCI_MEMHDL_T *mhloc=NULL, *mhrem=NULL;
      int nsegs, i;
/*       printf("%d(s): TRYING to use rdma contig to strided\n",armci_me); */
/*       fflush(stdout); */
      nsegs = 1;
      for(i=0;i<strides; i++) 
	nsegs *= count[i+1];
      if(nsegs<no_srv_copy_nsegs_ulimit() &&
	 msginfo->operation==GET && !msginfo->pinned && strides>=0 
	 && get_armci_region_local_hndl(ptr,armci_clus_id(armci_me),&mhloc)) {
/* 	printf("%d(s): using rdma contig to strided\n",armci_me); */
/* 	fflush(stdout); */
	armci_server_rdma_strided_to_contig(ptr, stride_arr,
					    count, strides,
					    msginfo->tag.data_ptr, to,
					    msginfo);
	return;
      }
      else {
/* 	printf("%d(s): not taking rdma to contig path. mhloc=%p mhrem=%p\n",armci_me,mhloc,mhrem); */
      }
    }
#endif
    /* for small contiguous blocks copy into a buffer before sending */
    armci_write_strided(ptr, strides, stride_arr, count, bdata);
    /* write the message to the client */
    armci_WriteToDirect(to, msginfo, bdata);

    if(DEBUG_){
        printf("%d(serv):sent len=%d to %d\n",armci_me,msginfo->datalen,to);
        fflush(stdout);
    }
}


/*\ server sends ACK to client
\*/
void armci_server_ack(request_header_t* msginfo)
{
     int ack=ACK;
     if(DEBUG_){
        printf("%d server: sending ACK to %d\n",armci_me,msginfo->from);
        fflush(stdout);
     }

     if(msginfo->datalen != sizeof(int))
        armci_die("armci_server_ack: bad datalen=",msginfo->datalen);
#if defined(PEND_BUFS)
     {
       /*Send from server known memory -- avoid extra buffers and
	 copying on server*/ 
       int *ack1 = (int *)(msginfo+1); /*msginfo is in some server
					buffer. I can overwrite the descriptor*/
       *ack1 = ACK;
       assert(sizeof(request_header_t)+2*sizeof(int)<IMM_BUF_LEN); 
       armci_send_data(msginfo, ack1);
     }
#else
     armci_send_data(msginfo, &ack);
#endif
}


/*  main routine for data server process in a cluster environment
 *  the process is blocked until message arrives from
 *  the clients and services the requests
 */
void armci_data_server(void *mesg)
{
    /* message */
    request_header_t *msginfo;
    void *descr;
    void *buffer;
    int buflen;
    int from;
#if defined(VAPI)
    static int mytag=1;
#endif

    /* read header, descriptor, data, and buffer length */
    armci_rcv_req(mesg, &msginfo, &descr, &buffer, &buflen );

    /* check what we got */
    armci_check_req(msginfo,buflen);
    from = msginfo->from;

    if(DEBUG_){ 
       printf("%d(serv):got %d request from %d\n",armci_me,msginfo->operation,
               from);
       fflush(stdout);
    }

/*if(msginfo->operation==GET)fprintf(stderr,"GET request received with tag: %d\n",msginfo->tag);*/

    switch(msginfo->operation){
      case ACK:
          if(DEBUG_) {
              fprintf(stdout, "%d(server): got ACK request from %d\n",
                      armci_me, msginfo->from); fflush(stdout);
          }
#ifdef SOCKETS
          armci_sock_send(msginfo->from, &(msginfo->tag), sizeof(msg_tag_t));
#endif
          armci_server_ack(msginfo);
          break;

      case ATTACH: 
          if(DEBUG_){
             printf("%d(serv):got ATTACH request from%d\n",armci_me, from);
             fflush(stdout);
          }
          armci_server_ipc(msginfo, descr, buffer, buflen);
          break;
#if defined(SOCKETS) || defined(HITACHI) || defined(MPI_SPAWN) || defined(MPI_MT)
      case QUIT:   
          if(DEBUG_){ 
             printf("%d(serv):got QUIT request from %d\n",armci_me, from);
             fflush(stdout);
          }
          armci_server_goodbye(msginfo);
          break;
#endif

      case ARMCI_SWAP:
      case ARMCI_SWAP_LONG:
      case ARMCI_FETCH_AND_ADD:
      case ARMCI_FETCH_AND_ADD_LONG:
          armci_server_rmw(msginfo,descr,buffer);
          break;

      case LOCK:
          armci_server_lock(msginfo);
          break;

      case UNLOCK:
          armci_server_unlock(msginfo, descr);
          break;

      default:
          if(msginfo->format ==VECTOR)
              armci_server_vector(msginfo, descr, buffer, buflen);
          else if(msginfo->format ==STRIDED){
#if defined(VAPI) /* buffer bypass protocol */
              if(msginfo->pinned == 1){
                  int armci_post_gather(void *, int *, int *,int, 
                                  armci_vapi_memhndl_t *,int,int,int,void *);
                  void * src_ptr;
                  int stride_levels;
                  int count[MAX_STRIDE_LEVEL];
                  int src_stride_arr[MAX_STRIDE_LEVEL];    
                  int found;
                  ARMCI_MEMHDL_T *mhandle;
                  int i,num,id;
                  
                  if(DEBUG1){
                     printf("%d(s) : unpacking dscr\n",armci_me);
                     fflush(stdout);
                  }
                  
                  src_ptr = *(void**)descr;
                  descr = (char*)descr + sizeof(void*);
                  stride_levels = *(int*)descr;
                  descr = (char*)descr + sizeof(int);
                  for(i =0; i<stride_levels;i++){
                      src_stride_arr[i] = *(int *)descr;
                      descr = (int *)descr + 1;
                  }
                  for(i =0;i<stride_levels+1;i++){
                      count[i] = *(int*)descr;
                      descr = (int*)descr + 1;   
                  
                  }

                  found = get_armci_region_local_hndl(src_ptr, armci_me,
                                 &mhandle);
                  if(!found){
                     armci_die("SERVER : local region not found",id);
                  }
                   
                  num =  armci_post_gather(src_ptr,src_stride_arr,
                                  count,stride_levels, mhandle,
                                  msginfo->from,mytag,SERV,NULL );
                  mytag =  (mytag+1)%MAX_PENDING;
                  if(mytag==0)mytag=1;
                  if(DEBUG1){
                     printf("%d(s) : finished posting %d gather\n", 
                                     armci_me,num);
                     fflush(stdout);
                  }     
                 
              }
              else        
#endif
                armci_server(msginfo, descr, buffer, buflen);
          }
          else
              armci_die2("armci_data_serv: unknown format code",
                         msginfo->format, msginfo->from);
    }
}


/*\ initialize connection and start server thread/processes
\*/
void armci_start_server()
{
    armci_init_connections();

#if defined(MPI_SPAWN) 
    
    /* For MPI_SPAWN, this should be called by all processes */
    armci_create_server_MPIprocess( );
    
#else
    
    if(armci_me == armci_master) {  
# ifdef SERVER_THREAD
       armci_create_server_thread( armci_server_code );
# else
       armci_create_server_process( armci_server_code );
# endif
    }
    
#endif
    
    armci_client_code();
    _armci_server_started=1;
}




void *armci_server_code(void *data)
{
#ifdef SERVER_THREAD
#if (defined(GM) || defined(VAPI) || defined(QUADRICS)) && ARMCI_ENABLE_GPC_CALLS
#  ifdef PTHREADS
  extern pthread_t data_server;
  data_server = pthread_self();
#  else  
  armci_die("armci_server_code: threaded data servers not using pthreads not supported by gpc", 0);
#  endif
#endif
#endif

    if(DEBUG_)
        printf("%d: in server after creating thread.\n",armci_me);

    /* make initial contact with all the computing process */
    armci_server_initial_connection();

    if(DEBUG_) {
        printf("%d(server): connected to all computing processes\n",armci_me);
        fflush(stdout);
    }
#if ARMCI_ENABLE_GPC_CALLS
    gpc_init();
#endif
    armci_call_data_server();

    armci_transport_cleanup();

    return(NULL);
}



/*\ request to QUIT sent by client
\*/
void armci_serv_quit()
{
int bufsize = sizeof(request_header_t)+sizeof(int);
int destproc;
request_header_t *msginfo;
destproc = SERVER_NODE(armci_clus_me);  
msginfo = (request_header_t*)GET_SEND_BUFFER(bufsize,QUIT,destproc);

    if(DEBUG_){ printf("%d master: sending quit request to server\n",armci_me);
        fflush(stdout);
    }

    msginfo->dscrlen = 0;
    msginfo->from  = armci_me;
    msginfo->to    = SERVER_NODE(armci_clus_me);
    msginfo->operation = QUIT;
    if(ACK_QUIT)
       msginfo->bytes   = msginfo->datalen = sizeof(int); /* ACK */
    else
       msginfo->bytes   = msginfo->datalen = 0; /* no ACK */

    armci_send_req(armci_master, msginfo, bufsize);

    if(ACK_QUIT){
       int stat;
       stat = *(int*)armci_rcv_data(armci_master,msginfo);  /* receive ACK */
       if(stat  != QUIT)
            armci_die("armci_serv_quit: wrong response from server", stat);
       FREE_SEND_BUFFER(msginfo);
    }
}


/*\ server action triggered by request to quit
\*/
void armci_server_goodbye(request_header_t* msginfo)
{
     int ack=QUIT;
     if(DEBUG_){
        printf("%d server: terminating request by %d\n",armci_me,msginfo->from);
        fflush(stdout);
     }

     if(msginfo->datalen){
       msginfo->datalen = -msginfo->datalen;
       if(msginfo->datalen != sizeof(int))
          armci_die("armci_server_goodbye: bad datalen=",msginfo->datalen);

       armci_send_data(msginfo, &ack);
     }

     armci_transport_cleanup();

     /* Finalizing data server process w.r.t. MPI is not portable
      */
     _exit(0);
}

