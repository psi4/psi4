/* ---------------------------------------------------------------------------------------------- *\
   portals.h header
\* ---------------------------------------------------------------------------------------------- */
 # ifndef _PORTALS_H_
 # define _PORTALS_H_

 # define PORTALS_INDEX    1

 # define ONE_KB 1024
 # define ONE_MB 1048576

 # define MAX_DS_MSG_SIZE ONE_MB

 # define PORTALS_MAX_DESCRIPTORS (MAX_BUFS+MAX_SMALL_BUFS)
 # define PORTALS_MAX_BUFS        MAX_BUFS
 # define PORTALS_MAX_SMALL_BUFS  MAX_SMALL_BUFS
 # define PORTALS_BUF_SIZE        MSG_BUFLEN                 /* defined in requesh.h */

/* define small buf length here - formerly request.h */
 # ifdef PORTALS_USE_RENDEZ_VOUS
 #  define PORTALS_SMALL_BUF_SIZE   1024 /* for use with nwchem only -- will not pass armci test.x */
 #  define PORTALS_MAX_EAGER_MESSAGE_SIZE PORTALS_SMALL_BUF_SIZE
 # else
 #  define PORTALS_SMALL_BUF_SIZE  1024
 #  define PORTALS_MAX_EAGER_MESSAGE_SIZE PORTALS_BUF_SIZE
 # endif

 # define PORTALS_NREQUEST_BUFFERS 40
 # define PORTALS_REQUEST_BUFFER_SIZE_WARNING (128*ONE_MB)

 # define PORTALS_READ_ACCESS     1
 # define PORTALS_WRITE_ACCESS 1000

 # define MATCH_ALL_MBITS     0x8000000000000000 /* should be set for all data requests */
 # define MATCH_ALL_IBITS     ~MATCH_ALL_MBITS   /* used to mask out all other bits, but MATCH_ALL */


 # define STATE_SEND_START        0x1
 # define STATE_SEND_END          0x2
 # define STATE_REPLY_START       0x4
 # define STATE_REPLY_END         0x8
 # define STATE_ACK               0x10
 # define STATE_PUT_START         0x20
 # define STATE_PUT_END           0x40
 # define STATE_GET_START         0x80
 # define STATE_GET_END           0x100
 # define STATE_UNLINK            0x200


 # define DS_RESPONSE_ACK         0x100000000000000
 # define DS_RESPONSE_PUT         0x200000000000000
 # define DS_RESPONSE_GET         0x400000000000000

 # define PORTALS_ALLOW_NBGETS
 # define PORTALS_USE_ARMCI_CLIENT_BUFFERS

 # define PORTALS_PUT_USE_ACK_TURNED_OFF
 # define PORTALS_PUT_USE_START_TURNED_OFF
 # define PORTALS_GET_USE_START_TURNED_OFF

/* ---------------------------------------------------------------------------------------------- *\
   portals types
\* ---------------------------------------------------------------------------------------------- */
   typedef struct portals_desc_s {
      void*             buffer; // used for the md
      ptl_size_t        length; // used for the md
      ptl_process_id_t  id;     // on whom the operation is acting on
      ptl_match_bits_t  mbits;  // operations destination mbits
      ptl_hdr_data_t    hdr;    // used for puts/unique counter value

      ptl_handle_ni_t   nih;    // network interface handle
      ptl_handle_eq_t   eqh;    // event handler
      ptl_handle_me_t   meh;    // me handle (if necessary)
      ptl_handle_md_t   mdh;    // md handle (if necessary)

      int               state;  // track outstanding events remaining on the descriptor
      int               done;   // flag to test whether all work on the descriptor is finished
      int               noperations; // the number of remote operations allowed on buffer
                                     // this is only used when preposting/pinning CP memory
                                     // for remote operations initiated by the data server
   } portals_desc_t;


   typedef struct portals_ds_req_s {
      portals_desc_t req_desc;
      portals_desc_t ack_desc;
      portals_desc_t data_desc;
      ptl_process_id_t dsid;
      size_t unique_msg_id;
      int active;
      int remote_node;
   } portals_ds_req_t;


/* ---------------------------------------------------------------------------------------------- *\
   portals global variables
\* ---------------------------------------------------------------------------------------------- */
   ptl_handle_ni_t   cp_nih;
   ptl_handle_ni_t   ds_nih;
   ptl_handle_eq_t   cp_eqh;
   ptl_handle_eq_t   ds_eqh;
   ptl_process_id_t *portals_id_map;
   ptl_process_id_t *portals_cloned_id_map;

   int portals_ds_ready;
   int portals_cp_finished;

   size_t portalsMaxEagerMessageSize;


/* ---------------------------------------------------------------------------------------------- *\
   portals prototypes
\* ---------------------------------------------------------------------------------------------- */
   int portals_init(ptl_handle_ni_t*);
   int portals_finalize(ptl_handle_ni_t);
   int portals_getid(ptl_handle_ni_t,ptl_process_id_t *);
   int portals_free_eq(ptl_handle_eq_t);
   int portals_create_eq(ptl_handle_ni_t, ptl_size_t, ptl_handle_eq_t*);
   int portals_create_matchall_me(ptl_handle_me_t*);
   int portals_me_attach(ptl_handle_ni_t,ptl_process_id_t,ptl_match_bits_t,ptl_match_bits_t,ptl_handle_me_t*);
   int portals_me_insert(ptl_handle_me_t,ptl_process_id_t,ptl_match_bits_t,ptl_match_bits_t,ptl_handle_me_t*);
   int portals_me_unlink(ptl_handle_me_t);
   int portals_md_attach(ptl_handle_me_t,ptl_md_t,ptl_unlink_t,ptl_handle_md_t*);
   int portals_md_bind(ptl_handle_ni_t,ptl_md_t,ptl_unlink_t,ptl_handle_md_t*);
 
   int portals_eqwait(ptl_handle_eq_t,ptl_event_t*);
   int portals_put(portals_desc_t*);
   int portals_get(portals_desc_t*);
   int portals_wait(portals_desc_t*);


   void* portalsCloneDataServer(void *);
   void  portalsSpinLockOnInt(volatile int*, int, int);

   void portals_print_event_details(ptl_event_t *ev);

   void Fatal_error(int);
   const char *Portals_ID();

   void bit_print(const char *,int);
   void hex_print(const char *,int);
   void portals_print_summary();

/* ---------------------------------------------------------------------------------------------- *\
   portals data server prototypes
\* ---------------------------------------------------------------------------------------------- */
   void* portals_ds_thread(void* args);
   int portals_ds_init(void);
   int portals_ds(void);
   int portal_send_test_ack(int to,int val);
   int portals_ds_requeue_md(int);

   void portals_ds_get_from_cp(void*,ptl_size_t,ptl_process_id_t,ptl_match_bits_t);

   //void ds_handler(DDI_Patch*,ptl_process_id_t);

/* ---------------------------------------------------------------------------------------------- *\
   portals compute process prototypes
\* ---------------------------------------------------------------------------------------------- */
   int portals_cp_init(void);
   int portals_cp_getid(ptl_process_id_t *id);

   void portals_req_send(void *buffer, size_t size, portals_ds_req_t *req);
   void portals_req_nbsend(void *buffer, size_t size, portals_ds_req_t *req);
   void portals_req_wait(portals_ds_req_t *req);

   void portals_remote_get(void *buffer, request_header_t *msginfo, int remote_node);
   void portals_remote_put(void *buffer, request_header_t *msginfo, int remote_node);
   void portals_remote_acc(void *buffer, request_header_t *msginfo, int remote_node);
   void portals_remote_rmw(void *buffer, request_header_t *msginfo, int remote_node, portals_ds_req_t *req);
   void portals_remote_nbget(void *buffer, request_header_t *msginfo, int remote_node, portals_ds_req_t *req);
   void portals_remote_nbput(void *buffer, request_header_t *msginfo, int remote_node, portals_ds_req_t *req);
   void portals_remote_nbacc(void *buffer, request_header_t *msginfo, int remote_node, portals_ds_req_t *req);


   void portalsRemoteOperation(void*,size_t,ptl_process_id_t,portals_ds_req_t*);
   void portalsRemoteOperationToRank(void*,size_t,int,portals_ds_req_t*);
   void portalsRemoteOperationToNode(void*,size_t,int,portals_ds_req_t*);

   void portalsBlockingRemoteOperationToNode(void*,size_t,int);


static inline unsigned int cpuid_ebx(unsigned int op)
{
        unsigned int eax, ebx;

        __asm__("cpuid"
                : "=a" (eax), "=b" (ebx)
                : "0" (op)
                : "cx", "dx" );
        return ebx;
}

 # endif
