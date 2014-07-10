#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* ---------------------------------------------------------------------------------------------- *\
   portals_cp.c  -- compute process portals calls
   author: ryan olson
   email:  ryan@cray.com
\* ---------------------------------------------------------------------------------------------- */
 # include "armcip.h"
 # include <assert.h>
 # include <stdio.h>
 # include <sched.h>
 # include <unistd.h>

/* ---------------------------------------------------------------------------------------------- *\
\* ---------------------------------------------------------------------------------------------- */
   static ptl_handle_ni_t cp_nih;
   static ptl_handle_eq_t cp_eqh;
   static ptl_handle_eq_t cp_tx_eqh;

   static void  *portals_eager_send_buffer  = NULL;
   static size_t portals_unique_msg_counter = 373;

   static int portals_smp_sem = -1;
   static int *active_requests_by_node = NULL;

/* ---------------------------------------------------------------------------------------------- *\
\* ---------------------------------------------------------------------------------------------- */
   int portals_cp_finished = 0;


/* ---------------------------------------------------------------------------------------------- *\
   Implementation
\* ---------------------------------------------------------------------------------------------- */

int 
portals_cp_init(void) 
{
        int rc;
        int me;
        ptl_process_id_t id;

        rc = portals_init(&cp_nih);
        if(rc != PTL_OK) {
           printf("error in portals_init: err %d\n",rc);
           Fatal_error(rc);
        }

        rc = portals_create_eq(cp_nih,10*PORTALS_MAX_DESCRIPTORS,&cp_eqh);
        if(rc != PTL_OK) {
           printf("failed to create cp event queue; err %d\n",rc);
           Fatal_error(911);
        }

        rc = portals_create_eq(cp_nih,30,&cp_tx_eqh);
        if(rc != PTL_OK) {
           printf("failed to create cp_tx event queue; err %d\n",rc);
           Fatal_error(911);
        }

        rc = portals_cp_getid(&id);
        if(rc != PTL_OK) {
           printf("failed to get the portals id; err %d\n",rc);
           Fatal_error(rc);
        }

     /* creating an smp/intra-node communicator */
        MPI_Comm_rank(ARMCI_COMM_WORLD,&me);
        MPI_Comm_split(ARMCI_COMM_WORLD,id.nid,me,&portals_smp_comm);

     /* set affinity */
      # ifdef PORTALS_AFFINITY
        int smp_np, smp_me;
        unsigned long mask;
        unsigned int len = sizeof(mask);
        unsigned long ncpus;
        unsigned int nsockets, siblings;
        int cores_per_socket, cps_per_socket;
        int verbose = 0;

        MPI_Comm_size(portals_smp_comm,&smp_np);
        MPI_Comm_rank(portals_smp_comm,&smp_me);


        if((ncpus = sysconf(_SC_NPROCESSORS_ONLN)) < 0) {
           printf("%d [cp] sysconf(_SC_NPROCESSORS_ONLN) failed; err=%d\n", ncpus);
           armci_die("sysconf in init_throttle",911);
        }

           
        if(sched_getaffinity(0, len, &mask) < 0) {
           perror("sched_getaffinity");
           armci_die("getaffinity error in ds_init",911);
        }

        if(armci_clus_me == 0 && /* verbose */ 0 ) {
           printf("%d [cp]: old affinity = 0x%x, ncpus = %d\n", armci_me, mask, ncpus);
        }

        if(smp_me == 0) {
           mask = 1 << (ncpus-1);
           if(sched_setaffinity(0, len, (cpu_set_t *) &mask) < 0) {
              perror("sched_setaffinity to probe the socket count");
              armci_die("setaffinity error in ds_init",911);
           }
           siblings = cpuid_ebx(1) >>16 & 0xff;
           nsockets = ncpus / siblings;
        } 
        MPI_Bcast(&nsockets,1,MPI_INT,0,portals_smp_comm);
           
        cores_per_socket = ncpus/nsockets;
        cps_per_socket  = (smp_np / nsockets);
        cps_per_socket += (smp_np % nsockets);
        if(nsockets > 2) {
           armci_die("nsockets > 2 not supported",911);
        }
        if(smp_me < cps_per_socket) {
           mask = 1 << smp_me;
        } else {
           mask = 1 << (smp_me + (cores_per_socket - cps_per_socket));
        }

        if(sched_setaffinity(0, len, (cpu_set_t *) &mask) < 0) {
           perror("sched_setaffinity");
           armci_die("setaffinity error in ds_init",911);
        }

        if(sched_getaffinity(0, len, &mask) < 0) {
           perror("sched_getaffinity");
           armci_die("getaffinity error (#2) in ds_init",911);
        }

        if(armci_clus_me == 0 && verbose) {
           printf("%d [cp]: new affinity = 0x%x, ncpus = %d\n", armci_me, mask, ncpus);
        }
      # endif

        return PTL_OK;
}
   

int
portals_cp_finalize()
{
        int rc;

      # ifdef PORTALS_LIMIT_REMOTE_REQUESTS_BY_NODE
        armci_semrm(portals_smp_sem);
      # endif

        rc = portals_free_eq(cp_eqh);
        if (rc != PTL_OK) {
           printf("error freeing cp_eqh; err %d\n",rc);
        }

        MPI_Barrier(ARMCI_COMM_WORLD);
        MPI_Finalize();

        portals_cp_finished = 1;
        exit(0);

        return PTL_OK;
//      return portals_finalize(cp_nih);
}


int
portals_cp_getid(ptl_process_id_t *id)
{
        return portals_getid(cp_nih, id);
}


static size_t
portals_get_unique_msg_id(void) {
   size_t val = armci_me*1000;
   portals_unique_msg_counter++;
   if(portals_unique_msg_counter == 1000) portals_unique_msg_counter=1;
   val += portals_unique_msg_counter;
   return val;
}


static void
portals_req_clear(portals_ds_req_t *req)
{
        req->active         = 0;
        req->unique_msg_id  = 0;

        req->req_desc.done  = 1;
        req->req_desc.state = 0;
        req->req_desc.eqh   = cp_tx_eqh;

        req->ack_desc.done  = 1;
        req->ack_desc.state = 0;
        req->ack_desc.eqh   = cp_eqh;

        req->data_desc.done  = 1;
        req->data_desc.state = 0;
        req->data_desc.eqh   = cp_eqh;

        req->remote_node     = -1;
}


static ptl_process_id_t
portals_get_dsid_from_node(int remote_node)
{
        int rank = armci_clus_info[remote_node].master;
        if(portals_cloned_id_map) return portals_cloned_id_map[rank];
        else                      return portals_id_map[rank];
}


static ptl_process_id_t
portals_get_dsid_from_rank(int remote_id)
{
        if(portals_cloned_id_map) return portals_cloned_id_map[remote_id];
        else                      return portals_id_map[remote_id];
}

void
portals_req_nbsend(void *buffer, size_t size, portals_ds_req_t *req)
{
        int rc;
        portals_desc_t *desc = &req->req_desc;

        assert(req->unique_msg_id);
        assert(size < portalsMaxEagerMessageSize);
        assert(req->remote_node >= 0);

     /* ---------------------------------------------------------------------------- *\
        if we get here, we can guarantee that where are no outstanding requests from
        this PE to the remote node; however, we can not guarantee that other PEs on
        this node aren't talking to the intended data server ... so now we wait on
        value in the "shared" array.
     \* ---------------------------------------------------------------------------- */
      # ifdef PORTALS_LIMIT_REMOTE_REQUESTS_BY_NODE
        int got_lock = 0;
        while(!got_lock) {
           portalsSpinLockOnInt(&active_requests_by_node[req->remote_node],0,1000);
           semaphoreAcquire(portals_smp_sem,1,PORTALS_WRITE_ACCESS);
           if(active_requests_by_node[req->remote_node] == 0) {
             active_requests_by_node[req->remote_node] = 1;
             got_lock = 1;
           }
           semaphoreRelease(portals_smp_sem,1,PORTALS_WRITE_ACCESS);
        }
      # endif

        desc->buffer = buffer;
        desc->length = size;
        desc->id     = req->dsid;
        desc->mbits  = MATCH_ALL_MBITS;
        desc->hdr    = req->unique_msg_id;
        desc->state  = 0;
        desc->eqh    = cp_tx_eqh;
        desc->nih    = cp_nih;

        rc = portals_put(desc);
        if(rc != PTL_OK) {
                printf("portals_put err %d\n",rc);
                Fatal_error(rc);
        }
}

void
portals_req_send(void *buffer, size_t size, portals_ds_req_t *req)
{
        int rc;
        portals_desc_t *desc = &req->req_desc;

        portals_req_nbsend(buffer,size,req);

        rc = portals_wait(desc);
        if(rc != PTL_OK) {
                printf("portals_wait err %d\n",rc);
                Fatal_error(rc);
        }
}


static inline void
portals_req_wait(portals_ds_req_t *req) 
{
        int rc;

        if(req->req_desc.state) {
           rc = portals_wait( &(req->req_desc) );
           if(rc != PTL_OK) {
              printf("portals wait error on req_desc in req_wait; err=%d\n",rc);
              Fatal_error(rc);
           }
        }

        if(req->ack_desc.state) {
           rc = portals_wait( &(req->ack_desc) );
           if(rc != PTL_OK) {
              printf("portals wait error on ack_desc in req_wait; err=%d\n",rc);
              Fatal_error(rc);
           }
        }
        if(req->data_desc.state) {
           rc = portals_wait( &(req->data_desc) );
           if(rc != PTL_OK) {
              printf("portals wait error on data_desc in req_wait; err=%d\n",rc);
              Fatal_error(rc);
           }
        }

        req->active = 0;
        return;
}


void
portalsWaitOnRequest(portals_ds_req_t *req) {
        portals_req_wait(req);
}


static int
portals_prepost_ack_from_ds(portals_ds_req_t *req)
{ 
        int rc;
        ptl_md_t md;
        portals_desc_t *desc = &req->ack_desc;
        unsigned long mbits = req->unique_msg_id;

        assert(req->unique_msg_id);
        assert(req->remote_node >= 0);

      # ifdef PORTALS_LIMIT_REMOTE_REQUESTS_BY_NODE
        desc->buffer = &active_requests_by_node[req->remote_node];
        desc->length = sizeof(int);
      # else
        desc->buffer = NULL;
        desc->length = 0;
      # endif
        desc->id     = req->dsid;
        desc->mbits  = mbits | DS_RESPONSE_ACK;
        desc->hdr    = mbits;
        desc->eqh    = cp_eqh;

        rc = portals_me_attach(cp_nih,desc->id,desc->mbits,0,&desc->meh);
        if(rc != PTL_OK) {
                printf("me failed in prepost ack\n");
                Fatal_error(rc);
        }

        md.start     = desc->buffer;
        md.length    = desc->length;
        md.threshold = 1;
        md.options   = PTL_MD_OP_PUT | PTL_MD_EVENT_START_DISABLE;
        md.user_ptr  = desc;
        md.eq_handle = cp_eqh;

        rc = portals_md_attach(desc->meh,md,PTL_UNLINK,&desc->mdh);
        if(rc != PTL_OK) {      
                printf("md failed in prepost ack\n");
                Fatal_error(rc);
        } 

        // desc->state = STATE_PUT_END;
        // |= needed for rendez-vous gets; put and get using the same descriptor
        desc->state |= STATE_PUT_END;
        desc->done  = 0;
}


static int
portals_prepost_put_from_ds(void *buffer, size_t size, portals_ds_req_t *req) 
{
        int rc;
        int nputs;
        ptl_md_t md;
        portals_desc_t *desc = &req->data_desc;
        unsigned long mbits = req->unique_msg_id;

        assert(req->unique_msg_id);

        desc->buffer = buffer;
        desc->length = size;
        desc->id     = req->dsid;
        desc->mbits  = mbits | DS_RESPONSE_PUT;
        desc->hdr    = mbits;
        desc->eqh    = cp_eqh;

        rc = portals_me_attach(cp_nih,desc->id,desc->mbits,0,&desc->meh);
        if(rc != PTL_OK) {
                printf("me failed in prepost put\n");
                Fatal_error(rc);
        }

        md.start     = buffer;
        md.length    = size;
        md.threshold = desc->noperations;
        md.options   = PTL_MD_OP_PUT
                     | PTL_MD_EVENT_AUTO_UNLINK_ENABLE
                     | PTL_MD_EVENT_START_DISABLE 
                     | PTL_MD_EVENT_END_DISABLE;
        md.user_ptr  = (void *) desc;
        md.eq_handle = cp_eqh;

        rc = portals_md_attach(desc->meh,md,PTL_UNLINK,&desc->mdh);
        if(rc != PTL_OK) {      
                printf("md failed in prepost put\n");
                Fatal_error(rc);
        } 

        // desc->state = STATE_UNLINK;
        // |= needed for rendez-vous gets; put and get using the same descriptor
        desc->state |= STATE_UNLINK;
        desc->done  = 0;
}


static int
portals_prepost_get_from_ds(void *buffer, size_t size, portals_ds_req_t *req) {
        
        int rc;
        ptl_md_t md;
        portals_desc_t *desc = &req->data_desc;
        unsigned long mbits = req->unique_msg_id;

        assert(req->unique_msg_id);

        desc->buffer = buffer;
        desc->length = size;
        desc->id     = req->dsid;
        desc->mbits  = mbits | DS_RESPONSE_GET;
        desc->hdr    = mbits;
        desc->eqh    = cp_eqh;

        rc = portals_me_attach(cp_nih,desc->id,desc->mbits,0,&desc->meh);
        if(rc != PTL_OK) {
                printf("me failed in prepost get\n");
                Fatal_error(rc);
        }

        md.start     = buffer;
        md.length    = size;
        md.threshold = desc->noperations;
        md.options   = PTL_MD_OP_GET
                     | PTL_MD_EVENT_START_DISABLE;
                  // | PTL_MD_EVENT_AUTO_UNLINK_ENABLE
                  // | PTL_MD_EVENT_START_DISABLE 
                  // | PTL_MD_EVENT_END_DISABLE;
        md.user_ptr  = (void *) desc;
        md.eq_handle = cp_eqh;

        rc = portals_md_attach(desc->meh,md,PTL_UNLINK,&desc->mdh);
        if(rc != PTL_OK) {      
                printf("md failed in prepost get\n");
                Fatal_error(rc);
        } 

     // printf("%d: preposted get of lenght=%ld\n",armci_me,size);
     // desc->state = STATE_UNLINK;
        // desc->state = STATE_GET_END;
        // |= needed for rendez-vous gets; put and get using the same descriptor
        desc->state |= STATE_GET_END;
        desc->done  = 0;
}


void portalsBlockingRemoteOperationToNode(void *buffer, size_t length, int remote_node) {
        portals_ds_req_t req;
        portals_req_clear(&req);        
        portalsRemoteOperationToNode(buffer,length,remote_node,&req);
        portalsWaitOnRequest(&req);
}


void portalsRemoteOperationToNode(void *buffer, size_t length, int remote_node, portals_ds_req_t *req)
{
        ptl_process_id_t id = portals_get_dsid_from_node(remote_node);
        req->remote_node = remote_node;
        portalsRemoteOperation(buffer,length,id,req);
}


/*
void portalsRemoteOperationToRank(void *buffer, size_t length, int remote_rank, portals_ds_req_t *req) {
        ptl_process_id_t id = portals_get_dsid_from_rank(remote_rank);
        portalsRemoteOperation(buffer,length,id,req);
}
*/


void
portalsRemoteOperation(void *buffer, size_t length, ptl_process_id_t dsid, portals_ds_req_t *req)
{
     /* --------------------------------------------------------------------- *\
        initialize the data server request
     \* --------------------------------------------------------------------- */
     // portals_req_clear(req);
        req->active = 1;
        req->unique_msg_id = portals_get_unique_msg_id();
        req->dsid = dsid;

     /* --------------------------------------------------------------------- *\
        the only response from the ds will be a 0-byte ack coming in as a put
     \* --------------------------------------------------------------------- */
        portals_prepost_ack_from_ds(req);

     /* --------------------------------------------------------------------- *\
        send data request; this is a completely blocking req
     \* --------------------------------------------------------------------- */
        portals_req_send(buffer,length,req);
}


void
portals_send_oper(int remote_node,int val, portals_ds_req_t *req)
{
        int rc;
        request_header_t msg;

     /* --------------------------------------------------------------------- *\
        initialize the data server request
     \* --------------------------------------------------------------------- */
        portals_req_clear(req);
        req->active = 1;
        req->unique_msg_id = portals_get_unique_msg_id();
        req->dsid = portals_get_dsid_from_node(remote_node);
        req->remote_node = remote_node;

     /* --------------------------------------------------------------------- *\
        the only response from the ds will be a 0-byte ack coming in as a put
     \* --------------------------------------------------------------------- */
        portals_prepost_ack_from_ds(req);

     /* --------------------------------------------------------------------- *\
        prepare data request and send it; this is a completely blocking req
     \* --------------------------------------------------------------------- */
        msg.operation = val;
        portals_req_send(&msg,sizeof(request_header_t),req);
        return; 
}


void
portals_send_QUIT(int remote_node)
{
        portals_ds_req_t req;
        portals_send_oper(remote_node,QUIT,&req);
        portals_req_wait(&req);
}


static int
portals_determine_remote_op_count(request_header_t *msg)
{
#ifdef DDI
        int nr,nc,np;
        int datatype_extent = sizeof(double);

     /* --------------------------------------------------------------------- *\
        previously we have worked with words, but to provide support for 
        other data types, we must work with bytes.  note to developers:
        datatype_extent = the size in bytes of the stored datatype
     \* --------------------------------------------------------------------- */
        if(msg->size*datatype_extent <= MAX_DS_MSG_SIZE) return 1;

     /* --------------------------------------------------------------------- *\
        the data must be moved in segments; determine patch dimensions
     \* --------------------------------------------------------------------- */
        nr = msg->ihi - msg->ilo + 1;
        nc = msg->jhi - msg->jlo + 1;
     
     /* --------------------------------------------------------------------- *\
        each column individually is too long to fit in the buffer
     \* --------------------------------------------------------------------- */
        if(nr*datatype_extent < MAX_DS_MSG_SIZE) {

        /* ------------------------------------------------------------------ *\
           np the number of "evenly" sized passed needed to send a column
        \* ------------------------------------------------------------------ */
           np = 2;
           while(((nr/np)+((nr%np)?1:0)*datatype_extent)>MAX_DS_MSG_SIZE) np++;

        /* ------------------------------------------------------------------ *\
           noperations is np times the number of columns to be sent
        \* ------------------------------------------------------------------ */
           return np*nc;

        }

     /* --------------------------------------------------------------------- *\
        determine the number of full columns that can be sent in one pass
        break down the subpatch on this metric
     \* --------------------------------------------------------------------- */
        else {

        /* ------------------------------------------------------------------ *\
           np is the number of passes needed to send the full patch which
           is broken down into "evenly" sized sets of columns that fit in 
           the allocated buffer region
        \* ------------------------------------------------------------------ */
           np = 2;
           while(nr*((nc/np)+((nc%np)?1:0))*datatype_extent>MAX_DS_MSG_SIZE) np++;

        /* ------------------------------------------------------------------ *\
           noperations is np 
        \* ------------------------------------------------------------------ */
           return np;

        }

        assert(0); // should not happen
        return -1;
#else
        return 1;
#endif
}

void 
portals_remote_rmw(void *buffer, request_header_t *msginfo, int remote_node, portals_ds_req_t *req)
{
        ptl_size_t length;

     /* --------------------------------------------------------------------- *\
        initialize the data server request
     \* --------------------------------------------------------------------- */
        portals_req_clear(req);
        req->active = 1;
        req->unique_msg_id = portals_get_unique_msg_id();
        req->dsid = portals_get_dsid_from_node(remote_node);
        req->remote_node = remote_node;

     /* --------------------------------------------------------------------- *\
        prepare the buffer into which the ds will put data
     \* --------------------------------------------------------------------- */
        req->data_desc.noperations=portals_determine_remote_op_count(msginfo);
        portals_prepost_put_from_ds(buffer,msginfo->datalen,req);

      # ifdef PORTALS_LIMIT_REMOTE_REQUESTS_BY_NODE
        portals_prepost_ack_from_ds(req);
      # endif

     /* --------------------------------------------------------------------- *\
        send data request
        note: from armci_send_req - if get, the value of bytes (local: length)
        is msginfo->dscrlen + (hdrlen=sizeof(request_header_t) ... this is
        the size of the "data server request message" to be sent
     \* --------------------------------------------------------------------- */
        length = sizeof(request_header_t) + msginfo->dscrlen + msginfo->datalen;
        portals_req_send(msginfo,length,req);
}

void 
portals_remote_get(void *buffer, request_header_t *msginfo, int remote_node)
{
        portals_ds_req_t req;
        portals_remote_nbget(buffer,msginfo,remote_node,&req);
        portals_req_wait(&req);
}
        
void 
portals_remote_nbget(void *buffer, request_header_t *msginfo, int remote_node, portals_ds_req_t *req)
{
        ptl_size_t length;

     /* --------------------------------------------------------------------- *\
        initialize the data server request
     \* --------------------------------------------------------------------- */
        portals_req_clear(req);
        req->active = 1;
        req->unique_msg_id = portals_get_unique_msg_id();
        req->dsid = portals_get_dsid_from_node(remote_node);
        req->remote_node = remote_node;

     /* --------------------------------------------------------------------- *\
        prepare the buffer into which the ds will put data
     \* --------------------------------------------------------------------- */
        req->data_desc.noperations=portals_determine_remote_op_count(msginfo);
        portals_prepost_put_from_ds(buffer,msginfo->datalen,req);

      # ifdef PORTALS_LIMIT_REMOTE_REQUESTS_BY_NODE
        portals_prepost_ack_from_ds(req);
      # endif

     /* --------------------------------------------------------------------- *\
        send data request
        note: from armci_send_req - if get, the value of bytes (local: length)
        is msginfo->dscrlen + (hdrlen=sizeof(request_header_t) ... this is
        the size of the "data server request message" to be sent
     \* --------------------------------------------------------------------- */
        length = sizeof(request_header_t) + msginfo->dscrlen;

      # if defined(PORTALS_USE_RENDEZ_VOUS)
        if(length < portalsMaxEagerMessageSize) portals_req_send(msginfo,length,req);
        else {
           req->data_desc.noperations = 1;
           portals_prepost_get_from_ds(msginfo,length,req);
        
        /* ------------------------------------------------------------------ *\
           send data request: branch here for eager vs. rendez-vous
        \* ------------------------------------------------------------------ */
           assert(length <= PORTALS_BUF_SIZE);
           portals_req_send(msginfo,sizeof(request_header_t),req);
        }
      # else   
        portals_req_send(msginfo,length,req);
      # endif
}


void 
portals_remote_put(void *buffer, request_header_t *msginfo, int remote_node)
{
        portals_ds_req_t req;
        portals_remote_nbput(buffer,msginfo,remote_node,&req);
        portals_req_wait(&req);
}


void 
portals_remote_nbput(void *buffer, request_header_t *msginfo, int remote_node, portals_ds_req_t *req)
{
        char *eagerBuffer = NULL;
        size_t eagerSendSize = 0;

     /* --------------------------------------------------------------------- *\
        initialize the data server request
     \* --------------------------------------------------------------------- */
        portals_req_clear(req);
        req->active = 1;
        req->unique_msg_id = portals_get_unique_msg_id();
        req->dsid = portals_get_dsid_from_node(remote_node);
        req->remote_node = remote_node;

     /* --------------------------------------------------------------------- *\
        prepost ack response from the data server
     \* --------------------------------------------------------------------- */
        portals_prepost_ack_from_ds(req);

     /* --------------------------------------------------------------------- *\
        eager vs. rendez-vous messaging
        eager: pack and send the message immediate (only for small messages)
        developers note: since portals_eager_send_buffer only exists once,
        this has to be a blocking send (ie the data is on the wire when
        req_send has finished and the buffer can be reused.  for greater
        overlap, create a set of eager send buffers ... however they have to
        be managed ... probably best to do it in a ring. 

        note: armci put/acc buffer is prepacked.
     \* --------------------------------------------------------------------- */
        eagerSendSize = sizeof(request_header_t) + msginfo->dscrlen + msginfo->datalen;
        if(eagerSendSize < portalsMaxEagerMessageSize) {
//         printf("sending eager message\n");
         # if 0 /* armci prepacked */
           eagerBuffer = (char *) portals_eager_send_buffer;
           memcpy(eagerBuffer,msginfo,sizeof(request_header_t));
           eagerBuffer += sizeof(request_header_t);
           memcpy(eagerBuffer,buffer,msginfo->bytes);
         # endif
           eagerBuffer = (char *) msginfo;  /* buffer == msginfo for armci */
           portals_req_send(eagerBuffer,eagerSendSize,req);
        }

     /* --------------------------------------------------------------------- *\
        rendez-vous: send the ds a request; ds will "get/pull" data 
     \* --------------------------------------------------------------------- */
        else {
         # ifdef PORTALS_USE_RENDEZ_VOUS
        /* ------------------------------------------------------------------ *\
           prepare the buffer into which the ds will put data
        \* ------------------------------------------------------------------ */
        // req->data_desc.noperations=portals_determine_remote_op_count(msginfo);
           req->data_desc.noperations = 1;
           portals_prepost_get_from_ds(msginfo,eagerSendSize,req);

        /* ------------------------------------------------------------------ *\
           send data request: branch here for eager vs. rendez-vous
        \* ------------------------------------------------------------------ */
           assert(eagerSendSize <= PORTALS_BUF_SIZE);
           portals_req_send(msginfo,sizeof(request_header_t),req);

         # else
           printf("%d [cp]: rendez-vous messaging not supported\n",armci_me);
           abort();
         # endif
        }
          
          
}


#if 0
void 
portals_remote_acc(void *buffer, request_header_t *msginfo, int remote_node)
{
        portals_ds_req_t req;
        portals_remote_nbacc(buffer,msginfo,remote_node,&req);
        portals_req_wait(&req);
}


void 
portals_remote_nbacc(void *buffer, request_header_t *msginfo, int remote_node, portals_ds_req_t *req)
{
        char *eagerBuffer = NULL;
        size_t eagerSendSize = 0;

        assert(msginfo->bytes);

     /* --------------------------------------------------------------------- *\
        initialize the data server request
     \* --------------------------------------------------------------------- */
        portals_req_clear(req);
        req->active = 1;
        req->unique_msg_id = portals_get_unique_msg_id();
        req->dsid = portals_get_dsid_from_node(remote_node);

     /* --------------------------------------------------------------------- *\
        eager vs. rendez-vous messaging
        eager: pack and send the message immediate (only for small messages)
     \* --------------------------------------------------------------------- */
        eagerSendSize = msginfo->bytes + sizeof(request_header_t);
        if(eagerSendSize < portalsMaxEagerMessageSize) {

        /* ------------------------------------------------------------------ *\
           prepost ack response from the data server
           developers note: if you globally fence an array with a collective
           operation prior to a section of code and defence it after, then you
           don't need to micro manage the fence on a per request basis in that
           section; this eliminates the need for a DS ack
        \* ------------------------------------------------------------------ */
           portals_prepost_ack_from_ds(req);

        /* ------------------------------------------------------------------ *\
           pack and send eager data request
           blocking for now, since portals_eager_send_buffer only exists once
           create multiple eager buffers for greater overlap
        \* ------------------------------------------------------------------ */
           eagerBuffer = (char *) portals_eager_send_buffer;
           memcpy(eagerBuffer,msginfo,sizeof(request_header_t));
           eagerBuffer += sizeof(request_header_t);
           memcpy(eagerBuffer,buffer,msginfo->bytes);
           eagerBuffer = (char *) portals_eager_send_buffer;
           portals_req_send(eagerBuffer,eagerSendSize,req);
        }

     /* --------------------------------------------------------------------- *\
        rendez-vous: send the ds a request; ds will "get/pull" data 
        developers note: a ds ack is not required for a rendez-vous pull,
        this is because the ds will not start the pull until a local fence
        has been raised (if needed - see note above)
     \* --------------------------------------------------------------------- */
        else {
        /* ------------------------------------------------------------------ *\
           prepare the buffer from which the ds will pull data
        \* ------------------------------------------------------------------ */
           req->data_desc.noperations=portals_determine_remote_op_count(msginfo);
           portals_prepost_get_from_ds(buffer,msginfo->bytes,req);

        /* ------------------------------------------------------------------ *\
           send data request
        \* ------------------------------------------------------------------ */
           portals_req_send(msginfo,sizeof(request_header_t),req);
           portalsWaitOnRequest(req);
        }
}
#endif

extern int armci_shmget(size_t,char*);
extern int armci_semget(int);
extern void *shmat(int,int,int);

void    
portals_cp_init_throttle(int nnodes) 
{         
        int i, shmid, smp_np, smp_me;
        size_t size = nnodes*sizeof(int);
        char *buf = NULL;
        
        MPI_Comm_size(portals_smp_comm,&smp_np);
        MPI_Comm_rank(portals_smp_comm,&smp_me);


      # ifdef PORTALS_LIMIT_REMOTE_REQUESTS_BY_NODE 
        if(armci_me == armci_master) {
           if(smp_me != 0) armci_die("smp_me and armci_master are different",911);
        }  
           
        if(smp_me == 0) {
           shmid = armci_shmget(size,"portals_cp_init_throttle");
           active_requests_by_node = (int *) shmat(shmid,0,0);
           if(active_requests_by_node == (void *) -1) {
              printf("%d [cp] shmat failed for shmid %d\n",armci_me,shmid);
              armci_die("badness",911);
           }
           armci_shmrm(shmid);
           for(i=0; i<nnodes; i++) active_requests_by_node[i] = 0;
           portals_smp_sem = armci_semget(2);
           semaphoreOperation(portals_smp_sem,1,PORTALS_WRITE_ACCESS);
        }
        
        MPI_Bcast(&shmid,1,MPI_INT,0,portals_smp_comm);
        MPI_Bcast(&portals_smp_sem,1,MPI_INT,0,portals_smp_comm);
        
        if(smp_me != 0) {
           active_requests_by_node = (int *) shmat(shmid,0,0);
        } else {
        }

        MPI_Barrier(portals_smp_comm);
      # endif
}

