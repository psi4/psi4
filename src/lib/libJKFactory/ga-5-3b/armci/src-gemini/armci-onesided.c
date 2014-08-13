#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "armcip.h"
int armci_onesided_ds_handler(void *);


#ifdef ARMCI_REGISTER_SHMEM
typedef struct {
       void *base_ptr;
       void *serv_ptr;
       size_t size;
       int islocal;
       int valid;
} aptl_reginfo_t;

typedef struct {
       aptl_reginfo_t reginfo[MAX_MEM_REGIONS];
       int reg_count;
} rem_meminfo_t;

static rem_meminfo_t *_rem_meminfo;
static aptl_reginfo_t *_tmp_rem_reginfo;
#define IN_REGION(_ptr__,_reg__) ((_reg__.valid) && (_ptr__)>=(_reg__.serv_ptr) \
                && (_ptr__) <= ( (char *)(_reg__.serv_ptr)+_reg__.size))
#endif

static cos_mdesc_t _send_mdesc, _recv_mdesc;
static cos_mdesc_t *send_mdesc = NULL;
static cos_mdesc_t *recv_mdesc = NULL;

int armci_onesided_direct_get_enabled = 1;
int armci_onesided_direct_put_enabled = 1;

cos_desc_t __global_1sided_direct_comm_desc;
cos_desc_t __global_1sided_direct_get_comm_desc;

// linked-list to hold mdh arrays for all ARMCI_Malloc calls
remote_mdh_node_t *remote_mdh_base_node = NULL;

char **client_buf_ptrs;

int 
armci_onesided_init()
{
        int i;
        cos_parameters_t cos_params;

        cos_params.options        = ONESIDED_DS_PER_NUMA;
        cos_params.nDataServers   = 1;
        cos_params.maxDescriptors = ARMCI_MAX_DESCRIPTORS*10;
        cos_params.maxRequestSize = ARMCI_MAX_REQUEST_SIZE;
        cos_params.dsHandlerFunc  = armci_onesided_ds_handler;

        bzero(&__global_1sided_direct_comm_desc,sizeof(cos_desc_t));

     // check to make sure things are properly sized
        if(armci_me == 0) {
        // ARMCI_ONESIDED_SIZEOF_IREQ is defined in armci.h
           if(sizeof(armci_ireq_t) != ARMCI_ONESIDED_SIZEOF_IREQ) {
              printf("ARMCI_ONESIDED_SIZEOF_IREQ is not sized correctly.\n");
              printf("ARMCI_ONESIDED_SIZEOF_IREQ = %d\nsizeof(armci_ireq_t) = %d\n",
                      ARMCI_ONESIDED_SIZEOF_IREQ,sizeof(armci_ireq_t));
              abort();
           }
        }

     // initialize libonesided
        COS_Init( &cos_params );

     // initialize armci memory
      # ifdef ARMCI_REGISTER_SHMEM
        _rem_meminfo = (rem_meminfo_t *)calloc(armci_nproc,sizeof(rem_meminfo_t));
        _tmp_rem_reginfo = (aptl_reginfo_t *)malloc(sizeof(aptl_reginfo_t)*armci_nproc);
        if( _rem_meminfo==NULL || _tmp_rem_reginfo ==NULL) {
           armci_die("malloc failed in init_portals",0);
        }
        if(armci_me == 0) {
           printf("sizeof(rem_meminfo_t)=%ld\n",sizeof(rem_meminfo_t));
        }
      # endif
        client_buf_ptrs = (char **) calloc(armci_nproc,sizeof(char *));
        assert(client_buf_ptrs);
        armci_msg_barrier();
        _armci_buf_init();

     // each armci buffer has a cos_request_t associated with it
     // initialize that cos_request_t now
     // moved into the above _armci_buf_init routine
     // for(i=0; i<MAX_BUFS; i++) cpReqCreate(&_armci_buffers[i].id.ar.req);
     // for(i=0; i<MAX_SMALL_BUFS; i++) cpReqCreate(&_armci_smbuffers[i].id.ar.req);

        return 0;
}



void 
armci_transport_cleanup()
{
    /*for i=0tomaxpendingclean*/
    ARMCI_PR_DBG("enter",0);
    free(client_buf_ptrs);
    ARMCI_PR_DBG("exit",0);
} 



static void 
armci_onesided_send(void *buffer, request_header_t *msginfo, int remote_node, cos_request_t *req)
{
        size_t length = sizeof(request_header_t) + msginfo->dscrlen + msginfo->datalen;

     // print_data(msginfo);
        cpReqInit(remote_node, req);
        cpPrePostRecv(buffer, length, req);
        cpCopyLocalDataMDesc(req, &msginfo->tag.response_mdesc);
        if(length > ARMCI_MAX_REQUEST_SIZE) length = sizeof(request_header_t); 
        cpReqSend(msginfo, length, req);
     // cpReqWait(req); // required until a new fence operation is created
}



void
print_data(void* buf)
{
        request_header_t *msginfo = (request_header_t *) buf;
        char *buffer = (char *) buf;
        buffer += sizeof(request_header_t) + msginfo->dscrlen;

        int ndouble = msginfo->datalen/8;
        double *data = (double *) buffer;

        printf("%d: [0]=%lf; [%d]=%lf; from=%d; to=%d\n",armci_me,data[0],ndouble-1,data[ndouble-1],msginfo->from, msginfo->to);
}



static void 
armci_onesided_recv(void* buffer, request_header_t *msginfo, int remote_node, cos_request_t *req)
{
        size_t length = sizeof(request_header_t) + msginfo->dscrlen;
        size_t reg_len = length + msginfo->datalen;

        cpReqInit(remote_node, req);
        cpPrePostRecv(buffer, reg_len, req);
        cpCopyLocalDataMDesc(req, &msginfo->tag.response_mdesc);
        if(length > ARMCI_MAX_REQUEST_SIZE) length = sizeof(request_header_t); 
        cpReqSend(msginfo, length, req);
}
        


static void
armci_onesided_oper(void* buffer, request_header_t *msginfo, int remote_node, cos_request_t *req)
{
        size_t length = sizeof(request_header_t);

        cpReqInit(remote_node, req);
        cpPrePostRecv(buffer, length, req);
        cpCopyLocalDataMDesc(req, &msginfo->tag.response_mdesc);
        cpReqSend(msginfo, length, req);
}



static void
armci_onesided_rmw(void *buffer, request_header_t *msginfo, int remote_node, cos_request_t *req)
{
        size_t length = sizeof(request_header_t) + msginfo->dscrlen + msginfo->datalen;

        cpReqInit(remote_node, req);
        cpPrePostRecv(buffer, msginfo->datalen, req);
        cpCopyLocalDataMDesc(req, &msginfo->tag.response_mdesc);
        cpReqSend(msginfo, length, req);
}

extern _buf_ackresp_t *_buf_ackresp_first,*_buf_ackresp_cur;


#if defined CRAY_REGISTER_ARMCI_MALLOC && HAVE_ONESIDED_FADD
void
armci_onesided_fadd(void *ploc, void *prem, int extra, int proc)
{
        onesided_hnd_t cp_hnd;
        cos_desc_t comm_desc;
        cos_mdesc_t local_mdh, remote_mdh, *mdh = NULL;

        cpGetOnesidedHandle(&cp_hnd);
        armci_onesided_search_remote_mdh_list(prem, proc, &remote_mdh);
        onesided_mem_register(cp_hnd, ploc, sizeof(long), NULL, &local_mdh);
        onesided_desc_init(cp_hnd, &local_mdh, &remote_mdh, 0, &comm_desc);
        onesided_fadd(extra, &comm_desc);
        onesided_wait(&comm_desc);
}
#endif

int
armci_send_req_msg(int proc, void *buf, int bytes, int tag)
{
        int cluster = armci_clus_id(proc);
        int serv    = armci_clus_info[cluster].master;
        char *buffer = (char *) buf;
        request_header_t *msginfo = (request_header_t *) buf;

      # ifdef ARMCI_LIMIT_REMOTE_REQUESTS_BY_NODE
        _armci_buf_ensure_one_outstanding_op_per_node(buf,cluster);
      # endif

      # ifdef SPECIAL_PUT_OPERATION_BROKEN_WHEN_INITIATED_FROM_USER_BUFFER
     // ensure any outstanding onesided direct operations have finished
        int state = __global_1sided_direct_comm_desc.state;
        onesided_wait(&__global_1sided_direct_comm_desc);
        if(state) cpMemDeregister(&__global_1sided_direct_comm_desc.local_mdesc);
      # endif

        BUF_INFO_T *bufinfo=_armci_buf_to_bufinfo(msginfo);
        _buf_ackresp_t *ar = &bufinfo->ar;
        cos_request_t *req = &ar->req;

        if(msginfo->operation == PUT || ARMCI_ACC(msginfo->operation)) {
           armci_onesided_send(buffer, msginfo, cluster, req);
        }

        else if(msginfo->operation == GET) {
        // move the buffer shift into the data server handler 
        // buffer = (char *) buf; 
        // buffer += sizeof(request_header_t);
        // buffer += msginfo->dscrlen;
           armci_onesided_recv(buffer, msginfo, cluster, req);
        }

        else if(msginfo->operation == ACK) {
           armci_onesided_oper(buffer, msginfo, cluster, req);
#if HAVE_ONESIDED_MEM_HTFLUSH
           onesided_mem_htflush(cluster);
#endif
        }

        else if(msginfo->operation == ARMCI_SWAP || msginfo->operation == ARMCI_SWAP_LONG ||
                msginfo->operation == ARMCI_FETCH_AND_ADD || 
                msginfo->operation == ARMCI_FETCH_AND_ADD_LONG) {
           buffer = (char *) buf;
           buffer += sizeof(request_header_t);
           buffer += msginfo->dscrlen; 
           armci_onesided_rmw(buffer, msginfo, cluster, req);
        }

        else {
           cosError("armci_send_req_msg: operation not supported",msginfo->operation);
        }


     // this had to be included in the portals version or shit would go down ... not sure y!
     // for now, we'll leave it in and see what happens later when we remote it
      # if 1
        ar->val = ar->valc = 0;
        if(ar==_buf_ackresp_first)_buf_ackresp_first=ar->next;
        if(ar->next!=NULL){
          ar->next->previous=ar->previous;
        }
        if(ar->previous!=NULL){
          ar->previous->next=ar->next;
          if(_buf_ackresp_cur==ar)_buf_ackresp_cur=ar->previous;
        }
        if(_buf_ackresp_cur==ar)_buf_ackresp_cur=NULL;
        ar->previous=ar->next=NULL;
      # endif

        return 0;
}



char *
armci_ReadFromDirect(int proc, request_header_t *msginfo, int len)
{
     // this is a CP funciton
        BUF_INFO_T *bufinfo = _armci_buf_to_bufinfo(msginfo);
        cos_request_t *req = &bufinfo->ar.req;
        cpReqWait(req);
        
     // return pointer to data
        char *ret = (char *) msginfo;
        ret += sizeof(request_header_t);
        ret += msginfo->dscrlen;
        return ret;
}



void
armci_WriteToDirect(int proc, request_header_t *msginfo, void *buf)
{
     // this is a DS function
        cos_desc_t resp_desc;
        cos_mdesc_t *resp_mdesc = &msginfo->tag.response_mdesc;
        dsDescInit(resp_mdesc, &resp_desc);
        resp_desc.event_type = EVENT_LOCAL | EVENT_REMOTE;
        if(send_mdesc == NULL) {
           send_mdesc = &_send_mdesc;
           dsMemRegister(MessageSndBuffer, sizeof(double)*MSG_BUFLEN_DBL, send_mdesc);
        }
        memcpy(&resp_desc.local_mdesc, send_mdesc, sizeof(cos_mdesc_t));
        resp_desc.local_mdesc.addr   = (uint64_t) buf;
        resp_desc.local_mdesc.length = (uint64_t) msginfo->datalen;
     // cosPut(buf, msginfo->datalen, &resp_desc);
        cosPutWithDesc(&resp_desc);
        dsDescWait(&resp_desc);
}



int armci_onesided_ds_handler(void *buffer)
{
        size_t length = 0;
        cos_desc_t get_desc;
        cos_mdesc_t *mdesc = NULL;
        void *buffer_to_data_server = buffer;
        request_header_t *request = (request_header_t *) buffer;
        if(request->operation == PUT || ARMCI_ACC(request->operation)) {
           length = sizeof(request_header_t) + request->dscrlen + request->datalen;
           if(length > ARMCI_MAX_REQUEST_SIZE) {
              char *get_buffer = (char *) MessageRcvBuffer;
              if(recv_mdesc == NULL) {
                 recv_mdesc = &_recv_mdesc;
                 dsMemRegister(MessageRcvBuffer, sizeof(double)*MSG_BUFLEN_DBL, recv_mdesc);
              }
              mdesc = &request->tag.response_mdesc;
              dsDescInit(mdesc, &get_desc);
              get_desc.event_type = EVENT_LOCAL;
              memcpy(&get_desc.local_mdesc, recv_mdesc, sizeof(cos_mdesc_t));
              get_desc.local_mdesc.length = length;
              assert(length <= sizeof(double)*MSG_BUFLEN_DBL);
              cosGetWithDesc(&get_desc);
              dsDescWait(&get_desc);
              buffer_to_data_server = (void *) get_buffer;
           }
        }
        else if(request->operation == GET) {
           length = sizeof(request_header_t) + request->dscrlen;
           if(length > ARMCI_MAX_REQUEST_SIZE) {
              // printf("[ds %d]: boom - rz fetch of get dscr\n",armci_me);
              char *get_buffer = (char *) MessageRcvBuffer;
              if(recv_mdesc == NULL) {
                 recv_mdesc = &_recv_mdesc;
                 dsMemRegister(MessageRcvBuffer, sizeof(double)*MSG_BUFLEN_DBL, recv_mdesc);
              }
              mdesc = &request->tag.response_mdesc;
              dsDescInit(mdesc, &get_desc);
              get_desc.event_type = EVENT_LOCAL;
              memcpy(&get_desc.local_mdesc, recv_mdesc, sizeof(cos_mdesc_t));
              get_desc.local_mdesc.length = length;
              assert(length <= sizeof(double)*MSG_BUFLEN_DBL);
              cosGetWithDesc(&get_desc);
              dsDescWait(&get_desc);
              buffer_to_data_server = (void *) get_buffer;
           }
           // regardless of rendez-vous or eager protocols
           // we have to shift the buffer and data lengths in the response_mdesc tag
           request = (request_header_t *) buffer_to_data_server;
           char *rbuf = (char *) request->tag.response_mdesc.addr;
           rbuf += length;
           request->tag.response_mdesc.addr = (uint64_t) rbuf;
           request->tag.response_mdesc.length -= length;
        } 

        if(request->operation == 0) {
           printf("%d [ds] possible zeroed buffer problem\n",armci_me);
           abort();
        }

        armci_data_server(buffer_to_data_server);
}



void 
armci_rcv_req(void *mesg,void *phdr,void *pdescr,void *pdata,int *buflen)
{
int i,na;
char *a;
double *tmp;

    request_header_t *msginfo = (request_header_t *)mesg;

    ARMCI_PR_SDBG("enter",msginfo->operation);
    *(void **) phdr = msginfo;

    if(0) {     
        printf("%d [ds]: got %d req (hdrlen=%d dscrlen=%d datalen=%d %d) from %d\n",
               armci_me, msginfo->operation, sizeof(request_header_t), msginfo->dscrlen,
               msginfo->datalen, msginfo->bytes,msginfo->from);
               fflush(stdout);
    }
    /* we leave room for msginfo on the client side */
    *buflen = MSG_BUFLEN - sizeof(request_header_t);
     

    if(send_mdesc == NULL) {
       send_mdesc = &_send_mdesc;
       dsMemRegister(MessageSndBuffer, sizeof(double)*MSG_BUFLEN_DBL, send_mdesc);
    } 

    // printf("%d [ds] oper=%d; bytes=%d\n",armci_me,msginfo->operation,msginfo->bytes);
    if(msginfo->bytes) {
       *(void **) pdescr = msginfo+1;
       *(void **) pdata  = msginfo->dscrlen + (char*)(msginfo+1);
          
       if(msginfo->operation == GET) {
          // the descriptor will exists after the request header
          // but there will be no data buffer
          // use the MessageRcvBuffer
          *(void**) pdata = MessageSndBuffer;
//        printf("%s (server) overriding pdata in rcv_req\n",Portals_ID());
          if(send_mdesc == NULL) {
             send_mdesc = &_send_mdesc;
             dsMemRegister(MessageSndBuffer, sizeof(double)*MSG_BUFLEN_DBL, send_mdesc);
          // printf("send_mdesc registered\n");
          // fflush(stdout);
          }
       }     
    }
    else {
    // printf("%d [ds]: hit this\n",armci_me);
       *(void**) pdescr = NULL;
       *(void**) pdata = MessageRcvBuffer;
       if(recv_mdesc == NULL) {
          recv_mdesc = &_recv_mdesc;
          dsMemRegister(MessageRcvBuffer, sizeof(double)*MSG_BUFLEN_DBL, recv_mdesc);
       }
    }
    ARMCI_PR_SDBG("exit",msginfo->operation);
}



void
armci_server_send_ack(request_header_t *msginfo)
{
     // this is a DS function
        cos_desc_t resp_desc;
        cos_mdesc_t *resp_mdesc = &msginfo->tag.response_mdesc;
        dsDescInit(resp_mdesc, &resp_desc);
        resp_desc.event_type = EVENT_LOCAL | EVENT_REMOTE;
        cosPut(NULL, 0, &resp_desc);
        dsDescWait(&resp_desc);
}



void
x_buf_wait_ack(request_header_t *msginfo, BUF_INFO_T *bufinfo)
{
        armci_die("x_buf_wait_ack not implemented",911);
}



void
x_net_send_ack(request_header_t *msginfo, int proc, void *dst, void *src)
{
        armci_die("x_net_send_ack not implemented",911);
}



long 
x_net_offset(char *buf, int proc)
{
        armci_die("x_net_offset not implemented",911);
      # if 0
        ARMCI_PR_DBG("enter",_rem_meminfo[proc].reg_count);
        if(DEBUG_COMM) { 
           printf("\n%d:%s:buf=%p",armci_me,__FUNCTION__,buf);fflush(stdout); 
        }
        for(i=0;i<_rem_meminfo[proc].reg_count;i++) {
            if(IN_REGION(buf,_rem_meminfo[proc].reginfo[i])) {
               return((long)((char *)_rem_meminfo[proc].reginfo[i].serv_ptr-(char *)_rem_meminfo[proc].reginfo[i].base_ptr));
            }
        }
        ARMCI_PR_DBG("exit",0);
      # endif
        return 0;
}


// currently our list of remote mdhs appears that it can get several entries with various
// lengths.  we should scan the mdh list first to see if an entry exists in the list
// if so, that could be an indication that the remote list entry function is not working
// properly, or that a different type of armci_free call is being used to by pass the
// removal of the mdh entry.  either way, we need to examine these occurences.
void
armci_onesided_append_remote_mdh_list(void* tgt_ptr, int proc, cos_mdesc_t *ret_mdh)
{

}

void
armci_onesided_search_remote_mdh_list(void* tgt_ptr, int proc, cos_mdesc_t *ret_mdh) 
{
        int node = armci_clus_id(proc);
        uint64_t length;
        uint64_t rem_addr;
        uint64_t tgt_addr = (uint64_t) tgt_ptr;
        remote_mdh_node_t *ll = remote_mdh_base_node;
        const cos_mdesc_t *mdh = NULL;

     // search the link-list for remote address and return the 
        while(ll) {
        // if we are in this routine, we are doing a direct onesided operations on a chuck of local
        // memory that was registered by the master process on this node.  typically, an armci operation
        // would work directly off the virtual address of that data as attached by the current process; 
        // however, because we are going to do a UGNI operation targetted at the MDH registered by the
        // armci_master rank on this node, we have to translate the virtual address on this rank to the
        // virtual address on armci_master.  this means we have to find the mdh by searching the ptrs
        // array and not the mdhs[*].addr values
           if(SAMECLUSNODE(proc) && armci_me != armci_master) {
              rem_addr = (uint64_t) ll->ptrs[proc];
           } else {
              rem_addr = (uint64_t) ll->mdhs[proc].addr;
           }
           length   = ll->mdhs[proc].length;
           if(tgt_addr >= rem_addr && tgt_addr < (rem_addr+length) /* check length of msg */) {
             mdh = &ll->mdhs[proc];
             break;
           }
           ll = ll->next;
        }

     // if remote mdh not found
        if(mdh == NULL) {
           printf("[cp %d]: warning - could not locate remote mdh for a direct put.\n",armci_me);
           printf("[cp %d]: searching for tgt_ptr=%p on node=%d / proc=%d\n",armci_me,tgt_ptr,node,proc);
           ll = remote_mdh_base_node;
           while(ll) {
              rem_addr = (uint64_t) ll->ptrs[proc];
              length   = ll->mdhs[proc].length;
              printf("[cp %d]: ll->ptrs[proc]=%p; ll->mdhs[node].length=%ld\n",armci_me,ll->ptrs[proc], length);
              ll = ll->next;
           }
           abort();
        }

     // setup return mdh
     // on the remote side the node master is the only rank that registers the "shared" memmory.  however,
     // shmat doesn't guarantee that all ranks on the node share the same starting virtual address.  that
     // is why we have to calculate the offset from the starting address on the node master based on the
     // actual virutal addresses on the remote rank.
        memcpy(ret_mdh, mdh, sizeof(cos_mdesc_t));
     // ret_mdh->addr += (tgt_addr-rem_addr);
     // if(ret_mdh->addr != tgt_addr) {
     //    printf("%d: ret_mdh->addr=%ld; tgt_addr=%ld\n",armci_me,ret_mdh->addr, tgt_addr);
     //    fflush(stdout);
     // }

     // if we are targeting a rank on the node for a direct operation, we need to translate the address
     // if not, then we can use the tgt_addr as passed in
        if(SAMECLUSNODE(proc) && armci_me != armci_master) {
           ret_mdh->addr += (tgt_addr-rem_addr);
        } else {
           ret_mdh->addr = tgt_addr;
        }
}

void
armci_onesided_remove_from_remote_mdh_list(void *tgt_ptr)
{
        cos_comm_t info;
        cos_mdesc_t *mdh = NULL;
        onesided_hnd_t cp_hnd;
        int node = armci_clus_id(armci_me);
        long total_bytes;
        remote_mdh_node_t *rm_ll, *ll = remote_mdh_base_node;

        NTK_MPI_GetComm(MPI_COMM_WORLD, &info);

     // get the onesided v2.0 api handle for the compute process
        cpGetOnesidedHandle(&cp_hnd);

     // find mdh
        while(ll) {
           if(tgt_ptr == ll->ptrs[armci_me]) {
              mdh = &ll->mdhs[armci_me];
              break;
           }
           ll = ll->next;
        }

     // ensure we have a valid mdh
        if(mdh == NULL) abort();

     // sum the total bytes allocated on the node
        MPI_Allreduce(&mdh->length, &total_bytes, 1, MPI_LONG, MPI_SUM, info.numa_comm);

     // node master only  
        if(info.numa_me == 0 && total_bytes) {

        // deregister memory
           onesided_mem_deregister(cp_hnd, mdh);
        // cpMemDeregister(mdh);
        }

     // free mdhs
        free(ll->mdhs);
        ll->mdhs = NULL;

     // update linked-list
        rm_ll = ll;
        if(rm_ll == remote_mdh_base_node) remote_mdh_base_node = rm_ll->next;
        else {
          ll = remote_mdh_base_node;
          while(ll->next != rm_ll) ll = ll->next;
          assert(ll->next == rm_ll);
          ll->next = rm_ll->next;
        }
        free(rm_ll);
}


void ARMCI_INIT_HANDLE(void *hdl)
{
        bzero(hdl, ARMCI_ONESIDED_SIZEOF_IREQ);
}


void armci_direct_on()
{
        armci_onesided_direct_get_enabled = 1;
        armci_onesided_direct_put_enabled = 1;
}

void armci_direct_off()
{
        armci_onesided_direct_get_enabled = 0;
        armci_onesided_direct_put_enabled = 0;
}

void armci_direct_on_() { armci_direct_on(); }
void armci_direct_off_() { armci_direct_off(); }
