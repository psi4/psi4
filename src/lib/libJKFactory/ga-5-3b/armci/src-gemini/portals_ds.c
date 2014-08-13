#if HAVE_CONFIG_H
#   include "config.h"
#endif

 # include "armcip.h"
 # include <assert.h>
 # include <stdio.h>
 # include <sched.h>
 # include <unistd.h>

static ptl_handle_ni_t ds_nih;
static ptl_handle_eq_t ds_eqh;
static ptl_handle_eq_t request_eqh;
static ptl_handle_me_t matchall_meh;

static int      request_buffer_cur_block;
static ptl_md_t request_buffer_md[PORTALS_NREQUEST_BUFFERS];
static ptl_handle_me_t request_buffer_meh[PORTALS_NREQUEST_BUFFERS];

int portals_ds_ready = 0;

// void *portals_ds_working_buffer = NULL;

void*
portals_ds_thread(void* args)
{
        portals_ds_init();
        portals_ds();
        portals_ds_finalize();
        portalsSpinLockOnInt(&portals_cp_finished,1,1000);
        exit(0);
        return NULL;
}


int
portals_ds_init()
{
        int i,rc;
        size_t bufferSize;
        float  warningSize;

        portals_ds_ready = 0;

     /* --------------------------------------------------------------------- *\
        unhook set affinity ... data servers can roam
     \* --------------------------------------------------------------------- */
      # ifdef PORTALS_AFFINITY
        int smp_np, smp_me;
        unsigned long mask;
        unsigned int len = sizeof(mask);
        unsigned long ncpus;
        int verbose = 0;

        MPI_Comm_size(portals_smp_comm,&smp_np);
        MPI_Comm_rank(portals_smp_comm,&smp_me);

        if((ncpus = sysconf(_SC_NPROCESSORS_ONLN)) < 0) {
           printf("%d [ds] sysconf(_SC_NPROCESSORS_ONLN) failed; err=%d\n", armci_me, ncpus);
           armci_die("sysconf in init_throttle",911);
        } 

        if(sched_getaffinity(0, len, &mask) < 0) {
           perror("sched_getaffinity");
           armci_die("getaffinity error in ds_init",911);
        }

        if(armci_clus_me == 0 && /* verbose */ 0 ) {
           printf("%d [ds]: old affinity = 0x%x, ncpus = %d\n", armci_me, mask, ncpus);
        }

        if(smp_np == ncpus) { 
           mask = (1 << ncpus) - 1; /* let the data server roam over all cores */
        } else {
           mask = 1 << (ncpus - 1); /* pin the ds to the last core on the node */   
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
           printf("%d [ds]: new affinity = 0x%x, ncpus = %d\n", armci_me, mask, ncpus);
        }
      # endif

     /* --------------------------------------------------------------------- *\
        initialize the network interface
     \* --------------------------------------------------------------------- */
        rc = portals_init(&ds_nih);
        if (rc != PTL_OK) {
                printf("failed to initialize portals on ds; err %d\n",rc);
                Fatal_error(rc);
        }

     /* --------------------------------------------------------------------- *\
        used for responding to data requests; this keeps the response events
        in a separate queue from the multitude of incoming data requests
     \* --------------------------------------------------------------------- */
        rc = portals_create_eq(ds_nih, 200, &ds_eqh);

     /* --------------------------------------------------------------------- *\
        used to process incoming data requests.  at very large scale we will
        have to do some sort of messaging by node group to reduce the worst
        case scenario off all to one type operations.  use the data server
        to message forward from node groups.
     \* --------------------------------------------------------------------- */
        i  = ARMCI_MAX(6*PORTALS_MAX_DESCRIPTORS*armci_nproc,200);
        i  = ARMCI_MAX(6*armci_nproc,200);
        rc = portals_create_eq(ds_nih, i, &request_eqh);
        if (rc != PTL_OK) {
                printf("failed to create request event queue");
                Fatal_error(rc);
        }

     /* --------------------------------------------------------------------- *\
        create ME list that matches all incoming data requests
        this will be a dead ME with no MD ... it will only be used as a
        place holder in which the "active" me/md will be placed in front of.
     \* --------------------------------------------------------------------- */
        rc = portals_create_matchall_me(&matchall_meh);
        if (rc != PTL_OK) {
                printf("failed to create matchall ME\n");
                Fatal_error(rc);
        }

     /* --------------------------------------------------------------------- *\
        create buffer space for the ds buffer
     \* --------------------------------------------------------------------- */
        assert(portalsMaxEagerMessageSize > sizeof(request_header_t));
        bufferSize = portalsMaxEagerMessageSize;

      # ifdef PORTALS_LIMIT_REMOTE_REQUESTS_BY_NODE
        bufferSize *= armci_nclus;
      # else
        bufferSize *= armci_nproc;
      # endif
        bufferSize = bufferSize/(PORTALS_NREQUEST_BUFFERS-2);
        bufferSize = ARMCI_MAX(bufferSize,portalsMaxEagerMessageSize);

        // if(armci_me == 0) printf("%s: bufferSize=%ld\n",Portals_ID(),bufferSize);
/*        
        if(bufferSize*PORTALS_NREQUEST_BUFFERS > PORTALS_REQUEST_BUFFER_SIZE_WARNING) {
           warningSize = (float) bufferSize * PORTALS_NREQUEST_BUFFERS;
           warningSize /= ONE_MB;
           printf("[data server]: internal request buffer is %.2f MB\n",warningSize);
        }
*/
        for(i=0; i<PORTALS_NREQUEST_BUFFERS; i++) {
                request_buffer_md[i].start     = malloc(bufferSize);
                request_buffer_md[i].length    = bufferSize;
                request_buffer_md[i].threshold = PTL_MD_THRESH_INF;
                request_buffer_md[i].max_size  = portalsMaxEagerMessageSize;
                request_buffer_md[i].options   = PTL_MD_OP_PUT 
                                               | PTL_MD_MAX_SIZE
                                               | PTL_MD_EVENT_AUTO_UNLINK_ENABLE
                                               | PTL_MD_EVENT_START_DISABLE;
                request_buffer_md[i].user_ptr  = (void *) (long) i;
                request_buffer_md[i].eq_handle = request_eqh;
        }

     /* --------------------------------------------------------------------- *\
        attach request buffers' ME/MD to the matchall ME placeholder
     \* --------------------------------------------------------------------- */
        for(i=0; i<PORTALS_NREQUEST_BUFFERS; i++) {
                portals_ds_requeue_md(i);
        }

        request_buffer_cur_block = 0;

     /* --------------------------------------------------------------------- *\
        sync with the compute processes
     \* --------------------------------------------------------------------- */
        portals_ds_ready = 1; 
        // printf("%s: ds ready\n",Portals_ID());

        return PTL_OK;
}



int
portals_ds(void)
{
        int rc = 0;
        int active = 1;
        size_t buffersize = 0;
        ptl_event_t ev;
        ptl_process_id_t from;
        request_header_t *request = NULL;
        char *buffer = NULL;

        do { 
        /* --------------------------------------------------------------------- *\
           wait for a portals event
        \* --------------------------------------------------------------------- */
           rc = portals_eqwait(request_eqh, &ev);
           if(rc != PTL_OK) {
              printf("eqwait failed in data_server\n");
              Fatal_error(911);
           }
           
        /* --------------------------------------------------------------------- *\
           process the event
        \* --------------------------------------------------------------------- */
           switch(ev.type) {
              case PTL_EVENT_PUT_START:
                   break;

              case PTL_EVENT_PUT_END:
              //   portals_print_event_details(&ev);
                /* ------------------------------------------------------------- *\
                   get the location of the request buffer
                \* ------------------------------------------------------------- */
                   buffer = (char *) ev.md.start;
                   buffer += ev.offset;
                   request = (request_header_t *) buffer;
                   request->tag.user_ptr = (void *) &ev;

                   if(request->operation == PUT || ARMCI_ACC(request->operation)) {
                      buffersize = sizeof(request_header_t) + request->dscrlen + request->datalen;
                      if(buffersize >= portalsMaxEagerMessageSize) {
                         buffer = (char *) MessageRcvBuffer;
                         portals_ds_get_from_cp(buffer,buffersize,ev.initiator,ev.hdr_data);
                         request = (request_header_t *) buffer;
                         request->tag.user_ptr = (void *) &ev;
                         armci_data_server(buffer);
                      // printf("%d: FINISHED RENDEZ-VOUS!\n",armci_me);
                         break;
                      }
                   }

                   if(request->operation == GET) {
                      buffersize = sizeof(request_header_t) + request->dscrlen;
                      if(buffersize >= portalsMaxEagerMessageSize) {
                         buffer = (char *) MessageRcvBuffer;
                         portals_ds_get_from_cp(buffer,buffersize,ev.initiator,ev.hdr_data);
                         request = (request_header_t *) buffer;
                         request->tag.user_ptr = (void *) &ev;
                         armci_data_server(buffer);
                      // printf("%d: FINISHED RENDEZ-VOUS!\n",armci_me);
                         break;
                      }
                   }

                /* ------------------------------------------------------------- *\
                   process request
                \* ------------------------------------------------------------- */
                   armci_data_server(buffer);
                   if(request->operation == QUIT) active = 0;
                   break;
  
              case PTL_EVENT_UNLINK:
//                 printf("captured an unlink event!!\n");
//                 portals_print_event_details(&ev);
                /* 
                   if((long) ev.md.user_ptr != request_buffer_cur_block) {
                      printf("sanity check failed: user_ptr=%ld; cur_block=%ld\n",(long) ev.md.user_ptr, request_buffer_cur_block);
                      armci_die("hummm ... unlink issue?",911);
                   }
                */
                   portals_ds_requeue_md((long) ev.md.user_ptr);
                   break;

              default:
                   printf("unexpected event type %d in recvany\n");
                   Fatal_error(911);
                   break;
           }

        } while(active);

// flush out event q; the only thing that should remain is possibly 1 unlink event;
        while( (rc=PtlEQGet(request_eqh, &ev)) != PTL_EQ_EMPTY) {
           if(rc == PTL_OK) {
              if(ev.type != PTL_EVENT_UNLINK) {
                 printf("%s: flushing request_eqh: event type=%d\n",Portals_ID(),ev.type);
              } else {
                 portals_ds_requeue_md((long) ev.md.user_ptr);
              }
           }
           else if(rc == PTL_EQ_DROPPED) {
              printf("%s: eq dropped\n",Portals_ID());
           } 
           else {
              printf("%s: some error in PtlEQGet; err=%d\n",Portals_ID(),rc);
              Fatal_error(rc);
           }
        }

        return PTL_OK;
}




int
portals_ds_finalize()
{
        int i,rc;

        // unlink and request buffers
        for(i=0; i<PORTALS_NREQUEST_BUFFERS; i++) {
           rc = portals_me_unlink(request_buffer_meh[i]);
           if(rc != PTL_OK) {
              printf("error unlinking me #%d (current working me = %d)\n",i,request_buffer_cur_block);
           }
           free(request_buffer_md[i].start);
        }

        rc = portals_free_eq(ds_eqh);
        if (rc != PTL_OK) {
           printf("error freeing ds_eqh; err %d\n",rc);
        }

        rc = portals_free_eq(request_eqh);
        if (rc != PTL_OK) {

        }

        return PTL_OK;
    //  return portals_finalize(ds_nih); // thread can't shutdown the nih before the cp is finished
}


int
portals_ds_getid(ptl_process_id_t *id)
{
        return portals_getid(ds_nih, id);
}


void
portals_print_event_details(ptl_event_t *ev)
{
        printf("%s [ds] event type=%d; offset=%d; mlength=%d; hdr_data=%lx; md.user_ptr=%ld\n",
                Portals_ID(), ev->type, ev->offset, ev->mlength, ev->hdr_data, (long) ev->md.user_ptr);
        fflush(stdout);
}


int
portals_ds_requeue_md(int i) 
{
        int rc;
        ptl_handle_me_t  meh;
        ptl_handle_md_t  mdh;
        ptl_process_id_t match_id;
        ptl_match_bits_t match_bits  = MATCH_ALL_MBITS;
        ptl_match_bits_t ignore_bits = MATCH_ALL_IBITS;

        match_id.nid = PTL_NID_ANY;
        match_id.pid = PTL_PID_ANY;

        rc = portals_me_insert(matchall_meh,match_id,match_bits,ignore_bits,&meh);
        if(rc != PTL_OK) {
                printf("me insert failed in ds requeue md; err %d\n",rc);
                Fatal_error(rc);
        }

        rc = portals_md_attach(meh,request_buffer_md[i],PTL_UNLINK,&mdh);
        if(rc != PTL_OK) {
                printf("md attach failed in ds requeue md; err %d\n",rc);
                Fatal_error(rc);
        }

        request_buffer_meh[i] = meh;
        request_buffer_cur_block++;
        if(request_buffer_cur_block == PORTALS_NREQUEST_BUFFERS) request_buffer_cur_block=0;

        return PTL_OK;
}


int
portals_create_matchall_me(ptl_handle_me_t* me_handle)
{
        int rc;
        ptl_process_id_t match_id;
        ptl_match_bits_t match_bits  = MATCH_ALL_MBITS;
        ptl_match_bits_t ignore_bits = MATCH_ALL_IBITS;

        match_id.nid = PTL_NID_ANY;
        match_id.pid = PTL_PID_ANY;

        rc = portals_me_attach(ds_nih,match_id,match_bits,ignore_bits,&matchall_meh);
            
        if (rc != PTL_OK) {
                printf("PtlMEAttachAny err %d in portals_create_melist\n",rc);
                return rc;
        }

        return rc;
}


void
portals_ds_send_ack(ptl_process_id_t id, ptl_match_bits_t mbits)
{
        portals_desc_t desc;
      # ifdef PORTALS_LIMIT_REMOTE_REQUESTS_BY_NODE
        static int ack = 0;
        desc.buffer = &ack;
        desc.length = sizeof(int);
      # else
        desc.buffer = NULL;
        desc.length = 0;
      # endif
        desc.id     = id;
        desc.mbits  = mbits | DS_RESPONSE_ACK;
        desc.hdr    = mbits;
        desc.state  = 0;
        desc.eqh    = ds_eqh;
        desc.nih    = ds_nih;
        portals_put(&desc);
        portals_wait(&desc);
}


void
portals_ds_send_put(void *buffer, ptl_size_t length, ptl_process_id_t id, ptl_match_bits_t mbits)
{
        portals_desc_t desc;
        desc.buffer = buffer;
        desc.length = length;
        desc.id     = id;
        desc.mbits  = mbits | DS_RESPONSE_PUT;
        desc.hdr    = mbits;
        desc.state  = 0;
        desc.eqh    = ds_eqh;
        desc.nih    = ds_nih;
        portals_put(&desc);
        portals_wait(&desc);
}


void
portals_ds_get_from_cp(void *buffer, ptl_size_t length, ptl_process_id_t id, ptl_match_bits_t mbits)
{
        portals_desc_t desc;
        desc.buffer = buffer;
        desc.length = length;
        desc.id     = id;
        desc.mbits  = mbits | DS_RESPONSE_GET;
        desc.hdr    = mbits;
        desc.state  = 0;
        desc.eqh    = ds_eqh;
        desc.nih    = ds_nih;
        portals_get(&desc);
        portals_wait(&desc);
}


#ifdef DDI
static void
ds_handler(DDI_Patch *request, ptl_process_id_t from)
{
        int i,j,nr,nc;
        long array[10],*a;
        size_t size;
        char *data_ptr;
        portals_desc_t desc;
        ptl_event_t *ev = (ptl_event_t *) request->user_ptr;

        switch(request->oper) {

        case DDI_GET:
//              printf("%s received DDI_GET request of size %d\n",Portals_ID(),request->size);
                nr = request->ihi - request->ilo + 1;
                nc = request->jhi - request->jlo + 1;
                if(nr < 0 || nc < 0 || nr > 10 || nc > 1) {
                   printf("test get dimension problem\n");
                   abort();
                }

                if(nr*sizeof(long) != request->size) {
                   printf("test get request size does not match\n");
                   abort();
                }

                for(i=0,j=317; i<nr; i++) array[i] = (long) j++;
                
                desc.buffer = &array[0];
                desc.length = nr*sizeof(long);
                desc.id     = ev->initiator;
                desc.mbits  = ev->hdr_data | DS_RESPONSE_PUT;
                desc.hdr    = ev->hdr_data;
                desc.state  = 0;
                desc.eqh    = ds_eqh;
                desc.nih    = ds_nih;
                portals_put(&desc);
                portals_wait(&desc);
                break;

        case DDI_PUT:
                nr = request->ihi - request->ilo + 1;
                nc = request->jhi - request->jlo + 1;

                data_ptr = NULL;
                if(ev->mlength > sizeof(DDI_Patch)) {
                   printf("recv'ed eager put - size %d\n",ev->mlength-sizeof(DDI_Patch));
                   data_ptr = (char *) request;
                   data_ptr += sizeof(DDI_Patch);
                }

                if(request->size != ev->mlength-sizeof(DDI_Patch)) {
                   printf("eager msg buffer length does not match request size %d\n",request->size);
                   abort();
                }

                a = (long *) data_ptr;
                for(i=0; i<nr; i++) printf("put[%d] = %ld\n",i,a[i]);
                portals_ds_send_ack(ev->initiator,ev->hdr_data);
                break;

        case DDI_QUIT:
//              printf("%s received DDI_QUIT request\n",Portals_ID());
                portals_ds_send_ack(ev->initiator,ev->hdr_data);
/*
                desc.buffer = NULL;
                desc.length = 0;
                desc.id     = ev->initiator;
                desc.mbits  = ev->hdr_data | DS_RESPONSE_ACK;
                desc.hdr    = ev->hdr_data;
                desc.state  = 0;
                portals_put(&desc);
                portals_wait(&desc);
*/
                break;

        case DDI_MEMORY:
                DDI_Memory_server(request->size);
                portals_ds_send_ack(ev->initiator,ev->hdr_data);
                break;

        default:
                printf("%s unknown operation in request=%d\n",Portals_ID(),request->oper);
                abort();
                break;
        }

        return;
}
#endif
