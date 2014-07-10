#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* ---------------------------------------------------------------------------------------------- *\
   portals.c  -- wrapper for commonly used portals calls
   author: ryan olson
   email:  ryan@cray.com
\* ---------------------------------------------------------------------------------------------- */
 # include "armcip.h"

/* ---------------------------------------------------------------------------------------------- *\
   global variables
\* ---------------------------------------------------------------------------------------------- */
   ptl_process_id_t *portals_id_map = NULL;
   ptl_process_id_t *portals_cloned_id_map = NULL;

   size_t portalsMaxEagerMessageSize;

   MPI_Comm portals_smp_comm;

/* ---------------------------------------------------------------------------------------------- *\
   static variables for this object
\* ---------------------------------------------------------------------------------------------- */
   static int portals_verbose = 0;


/* ---------------------------------------------------------------------------------------------- *\
   portals wrappers
\* ---------------------------------------------------------------------------------------------- */


int
portals_init(ptl_handle_ni_t *nih)
{
        int num_interfaces = 0;
        int rc;

        rc = PtlInit(&num_interfaces);
        if (rc != PTL_OK) {
                printf("PtlInit err %d\n", rc);
                return rc;
        }

        rc = PtlNIInit(CRAY_UK_SSNAL, PTL_PID_ANY, NULL, NULL, nih);
        if (rc != PTL_OK && rc != PTL_IFACE_DUP) {
                printf("PtlNIInit err %d\n", rc);
                return rc;
        }

        portalsMaxEagerMessageSize = PORTALS_MAX_EAGER_MESSAGE_SIZE;

        return PTL_OK;
}


int
portals_finalize(ptl_handle_ni_t nih)
{
        PtlNIFini(nih);
        PtlFini();
        return PTL_OK;
}


int
portals_getid(ptl_handle_ni_t nih, ptl_process_id_t *id)
{
        int rc;

        rc = PtlGetId(nih, id);
        if(rc != PTL_OK) {
                printf("PtlGetId err %d\n",rc);
                return rc;
        }

        return PTL_OK;
}


int
portals_create_eq(ptl_handle_ni_t nih, ptl_size_t count, ptl_handle_eq_t *eq_handle) 
{
        int rc;

        rc = PtlEQAlloc(nih, count, PTL_EQ_HANDLER_NONE, eq_handle);
        if (rc != PTL_OK) {
                printf("PtlEQAlloc err %d\n", rc);
                return rc;
        }

        return PTL_OK;
}


int
portals_free_eq(ptl_handle_eq_t eq)
{
        int rc;

        rc = PtlEQFree(eq);
        if (rc != PTL_OK) {
                printf("PtlEQFree err %d\n",rc);
                return rc;
        }

        return PTL_OK;
}

/* 
   permanent buffers - such as unexpected receive buffers or data requests
   buffers should not be unlinked.  client side buffers, such as large puts/accs
   would create a ME in front of the MATCH ALL unexpected buffer/data req ME. 
   on the client side, the MATCH ALL ME should catch the ACKs 
*/


int
portals_me_attach(ptl_handle_ni_t  nih,
                  ptl_process_id_t match_id,
                  ptl_match_bits_t match_bits,
                  ptl_match_bits_t ignore_bits,
                  ptl_handle_me_t *me_handle)
{
        int rc = PtlMEAttach(nih,PORTALS_INDEX,match_id,match_bits,ignore_bits,
                             PTL_UNLINK,PTL_INS_BEFORE,me_handle);
        if (rc != PTL_OK) {
                printf("PtlAttach err %d in me_attach\n",rc);
                return rc;
        }

        return PTL_OK;
}

int
portals_me_insert(ptl_handle_me_t  base,
                  ptl_process_id_t pe_match_id, 
                  ptl_match_bits_t match_bits,
                  ptl_match_bits_t ignore_bits,
                  ptl_handle_me_t *me_handle)
{
        int rc = PtlMEInsert(base,pe_match_id,match_bits,ignore_bits,
                             PTL_UNLINK,PTL_INS_BEFORE,me_handle);
        if (rc != PTL_OK) {
                printf("PtlME err %d in portals_push_me\n",rc);
                return rc;
        }

        return rc;
}


int
portals_me_unlink(ptl_handle_me_t meh)
{
        int rc = PtlMEUnlink(meh);

        if(rc != PTL_OK) {
           printf("PtlMEUnlink err %d in me_unlink\n",rc);
        }

        return rc;
}


int
portals_md_attach(ptl_handle_me_t me_handle,
                  ptl_md_t md,
                  ptl_unlink_t unlink_op,
                  ptl_handle_md_t *md_handle)
{
        int rc = PtlMDAttach(me_handle, md, unlink_op, md_handle);
        if (rc != PTL_OK) {
                printf("PtlMDAttach err %d\n",rc);
                return rc;
        }

        return PTL_OK;
}


int 
portals_md_bind(ptl_handle_ni_t nih,
                ptl_md_t md,
                ptl_unlink_t unlink_op,
                ptl_handle_md_t *md_handle)
{
        int rc = PtlMDBind(nih, md, unlink_op, md_handle);
        if (rc != PTL_OK) {
                printf("PtlMDBind err %d\n",rc);
                return rc;
        }

        return rc;
}


int
portals_eqwait(ptl_handle_eq_t eqh, ptl_event_t *ev) 
{
        int rc = PtlEQWait(eqh, ev);
        if (rc != PTL_OK) {
                printf("PtlEQWait err %d\n",rc);
                return rc;
        }

        return  PTL_OK;
}


static int
notify(portals_desc_t *desc, int state, char *name) {
        if(desc->state & state) {
                desc->state &= ~state;
                if(desc->state == 0) desc->done = 1;
                return 1;
        } else {
                printf("event: %s with desc state %x not %x\n",name,desc->state,state);
                abort();
                return 0;
        }
}


int
portals_wait(portals_desc_t *wait_on_desc) {

        int rc;
        ptl_event_t ev;
        portals_desc_t *desc = NULL;

        while(wait_on_desc->state) {

                rc = portals_eqwait(wait_on_desc->eqh, &ev);
                if (rc != PTL_OK) {
                        printf("eq wait error in portals_wait\n");
                        abort();
                }
                
                desc = (portals_desc_t *) ev.md.user_ptr;

                switch(ev.type) {

                case PTL_EVENT_SEND_START:
                         if (portals_verbose) printf("%s event: send start\n",Portals_ID());
                         notify(desc, STATE_SEND_START, "send start");
                         break;
 
                case PTL_EVENT_SEND_END:
                         if (portals_verbose) printf("%s event: send end\n",Portals_ID());
                         notify(desc, STATE_SEND_END, "send end");
                         break;
 
                case PTL_EVENT_REPLY_START:
                         if (portals_verbose) printf("%s event: reply start\n",Portals_ID());
                         notify(desc, STATE_REPLY_START, "reply start");
                         break;
 
                case PTL_EVENT_REPLY_END:
                         if (portals_verbose) printf("%s event: reply end\n",Portals_ID());
                         notify(desc, STATE_REPLY_END, "reply end");
                         break;
 
                case PTL_EVENT_ACK:
                         if (portals_verbose) printf("%s event: ack\n",Portals_ID());
                         printf("%s event ack: md.threshold=%d\n",Portals_ID(),ev.md.threshold);
                         notify(desc, STATE_ACK, "ack");
                         break;
 
                case PTL_EVENT_PUT_START:
                         if (portals_verbose) printf("%s event: put start\n",Portals_ID());
                         notify(desc, STATE_PUT_START, "put start");
                         break;
 
                case PTL_EVENT_PUT_END:
                         if (portals_verbose) printf("%s event: put end\n",Portals_ID());
                         if (notify(desc, STATE_PUT_END, "put end")) {
                         //      desc->len = ev.mlength;
                         //      desc->off = ev.offset;
                         }
                         break;

                case PTL_EVENT_GET_START:
                        if (portals_verbose) printf("%s event: get start\n",Portals_ID());
                        notify(desc, STATE_GET_START, "get start");
                        break;

                case PTL_EVENT_GET_END:
                        if (portals_verbose) printf("%s event: get end\n",Portals_ID());
                        notify(desc, STATE_GET_END, "get end");
                        break;

                case PTL_EVENT_UNLINK:
                        if (portals_verbose) printf("%s event: unlink\n",Portals_ID());
                        notify(desc, STATE_UNLINK, "unlink");
                        break;

                default:
                        printf("%s event: %d\n",Portals_ID(), ev.type);
                        break;
                }

        }

        return PTL_OK;
}


int
portals_put(portals_desc_t *desc)
{
        int rc;
        int threshold = 1;
        ptl_md_t md = { 0 };
        ptl_handle_md_t md_handle;
        ptl_ack_req_t ack_req = PTL_NOACK_REQ; 

      # ifdef PORTALS_PUT_USE_ACK
        ack_req = PTL_ACK_REQ;
        threshold++;
      # endif

        md.start     = desc->buffer;
        md.length    = desc->length;
        md.threshold = threshold;
        md.options   = 0;
      # ifndef PORTALS_PUT_USE_START
        md.options  |= PTL_MD_EVENT_START_DISABLE;
      # endif
        md.user_ptr  = desc;
        md.eq_handle = desc->eqh;

        rc = portals_md_bind(desc->nih, md, PTL_UNLINK, &md_handle);
        if (rc != PTL_OK) {
                printf("failed to bind local md in put; err %d",rc);
                Fatal_error(rc);
        }

        rc = PtlPut(md_handle,
                    ack_req,
                    desc->id,
                    PORTALS_INDEX,
                    0,
                    desc->mbits,
                    0,
                    desc->hdr);
        if (rc != PTL_OK) {
                printf("PtlPut err %d\n",rc);
                return rc;
        }

        desc->done  = 0;
        desc->state = STATE_SEND_END;

      # ifdef PORTALS_PUT_USE_START
        desc->state |= STATE_SEND_START;
      # endif

      # ifdef PORTALS_PUT_USE_ACK
        desc->state |= STATE_ACK;
      # endif

        return PTL_OK;
}


int
portals_get(portals_desc_t* desc)
{
        int rc;
        ptl_md_t md = { 0 };
        ptl_handle_md_t md_handle;

        md.start     = desc->buffer;
        md.length    = desc->length;
        md.threshold = 2;
        md.options   = 0;
      # ifndef PORTALS_GET_USE_START
        md.options  |= PTL_MD_EVENT_START_DISABLE;
      # endif
        md.user_ptr  = desc;
        md.eq_handle = desc->eqh;

        rc = portals_md_bind(desc->nih, md, PTL_UNLINK, &md_handle);
        if (rc != PTL_OK) {
                printf("failed to bind local md in get; err %d\n",rc);
                Fatal_error(rc);
        }

        rc = PtlGet(md_handle,
                    desc->id,
                    PORTALS_INDEX,
                    0,
                    desc->mbits,
                    0);
        if (rc != PTL_OK) {
                printf("PtlGet err %d\n",rc);
                Fatal_error(rc);
        }

        desc->done  = 0;
        desc->state = STATE_REPLY_END | STATE_SEND_END;

      # ifdef PORTALS_GET_USE_START
        desc->state |= STATE_REPLY_START;
      # endif

        return PTL_OK;
}


/*
portals_desc_t*
portals_get_free_desc(void)
{
        int i,rc;
        portals_desc_t *desc = NULL;

        while(desc == NULL) {
           for(i=0; i<PORTALS_MAX_DESCRIPTORS; i++) {
              if(portals_desc_list[i].done) {
                 desc = &portals_desc_list[i];
                 break;
              }
           }
           if(desc) break;
           portals_wait_any(&portals_desc_list,PORTALS_MAX_DESCRIPTORS);
        }

        bzero(desc,sizeof(portals_desc_t));
        return desc;
}
*/

void * 
portalsCloneDataServer( void * (*func)(void *) )
{
      char *stack = malloc(256*ONE_KB);
      char *stack_top = &stack[256*ONE_KB-1];
      pid_t pid = clone(func, (void *) stack_top,
                        CLONE_THREAD|CLONE_SIGHAND|CLONE_VM, NULL);

      if(pid == -1) {
         printf("clone failed in systemCreateClone\n");
         Fatal_error(911);
      }

      // data_server_pid = pid;
      return;
}


void
portalsSpinLockOnInt(volatile int *ptr, int val, int maxspin)
{
        int count = 0;
        extern void cpu_yield();

        while(*ptr != val) {

           if(++count < maxspin);
           else{
//             cpu_yield();
               count =0;
             # if defined(MACOSX) && defined(__ppc__) && defined(__GNUC__)
               __asm__ __volatile__ ("sync" ::: "memory");
             # endif
               __asm__ __volatile__ ("mfence" ::: "memory");
               __asm__ __volatile__ ("sfence" ::: "memory");
            }
        }
}


void Fatal_error(int rc) {
        armci_die("rmo fatal error",rc);
}

const char *Portals_ID() {
        static char string[128];
        sprintf(string,"[%i]:",armci_me);
        return string;
}


void
hex_print(char* data, int length)
{
        int ptr = 0;
        for(;ptr < length;ptr++)
        {
                printf("0x%02x ",(unsigned char)*(data+ptr));
        }
        printf("\n");
}

void
bit_print(const char* data, int length)
{
        unsigned char mask = 0x01;
        int ptr = 0;
        int bit = 0;
        for(;ptr < length;ptr++)
        {
                for(bit = 7;bit >= 0;bit--)
                {
                        if ((mask << bit) & (unsigned char)*(data+ptr))
                        {
                                printf("1");
                        }
                        else
                        {
                                printf("0");
                        }
                }
                printf(" ");
        }
        printf("\n");
}


void
portals_print_summary()
{
        printf("PORTALS_MAX_DESCRIPTORS        = %d\n",PORTALS_MAX_DESCRIPTORS);
        printf("PORTALS_MAX_BUFS               = %d\n",PORTALS_MAX_BUFS);
        printf("PORTALS_MAX_SMALL_BUFS         = %d\n",PORTALS_MAX_SMALL_BUFS);
        printf("PORTALS_BUF_SIZE               = %d\n",PORTALS_BUF_SIZE);
        printf("PORTALS_SMALL_BUF_SIZE         = %d\n",PORTALS_SMALL_BUF_SIZE);
        printf("PORTALS_NREQUEST_BUFFERS       = %d\n",PORTALS_NREQUEST_BUFFERS);
        printf("PORTALS_MAX_EAGER_MESSAGE_SIZE = %d\n",PORTALS_MAX_EAGER_MESSAGE_SIZE);
        return;
}
