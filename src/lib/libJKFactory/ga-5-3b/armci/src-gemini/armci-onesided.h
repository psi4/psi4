#ifndef __ARMCI_ONESIDED_H__
#define __ARMCI_ONESIDED_H__

#include "onesided.h"

#define NUM_SERV_BUFS           1
#define MAX_MEM_REGIONS         30

#define ARMCI_BUF_SIZE          262144
#define ARMCI_SMALL_BUF_SIZE    2048

#define ARMCI_MAX_BUFS          4
#define ARMCI_MAX_SMALL_BUFS    8

#define ARMCI_MAX_DESCRIPTORS   (ARMCI_MAX_BUFS+ARMCI_MAX_SMALL_BUFS)
#define ARMCI_MAX_REQUEST_SIZE  ARMCI_SMALL_BUF_SIZE

/*
 There is a problem with ga_transpose when CRAY_REGISTER_ARMCI_MALLOC
 is defined.

 The fix is a special hook in ga_transpose to turn off the direct puts
 during the transpose.  This may indicate a race condition in the 
 transpose code or a problem with the direct fencing.
*/
#define CRAY_REGISTER_ARMCI_MALLOC
#define ARMCI_LIMIT_REMOTE_REQUESTS_BY_NODE_TURNED_OFF
#define MAX_OUTSTANDING_ONESIDED_GETS 64

#define ARMCI_ONESIDED_GETS_USES_NBGETS

/* typedefs */

typedef struct armci_onesided_msg_tag_s {
        int msgid;
        cos_mdesc_t response_mdesc;
} armci_onesided_msg_tag_t;


typedef struct remote_mdh_node {
        void **ptrs;
        cos_mdesc_t *mdhs;
        struct remote_mdh_node *next;
} remote_mdh_node_t;

// linked-list of remote mdhs
// a new node is created on each ARMCI_Malloc operation
// not an ideal scenario -- perhaps use Abhinav's new dreg routines
// to store this data ... it would require two pieces of info ...
// remote target and remote virtual addr ... and return the mdh
// for now: manually look up the mdh in ARMCI_GetS and ARMCI_PutS
extern remote_mdh_node_t *remote_mdh_base_node;

/* functions */
int armci_onesided_init();
void armci_transport_cleanup();
void armci_rcv_req(void *,void *,void *,void *,int *);

void print_data(void *);

void armci_onesided_search_remote_mdh_list(void* tgt_ptr, int proc, cos_mdesc_t *mdh);
void armci_onesided_remove_from_remote_mdh_list(void *tgt_ptr);

#if defined CRAY_REGISTER_ARMCI_MALLOC && HAVE_ONESIDED_FADD
void armci_onesided_fadd(void *ploc, void *prem, int extra, int proc);
#endif

extern int armci_onesided_direct_get_enabled;
extern int armci_onesided_direct_put_enabled;
extern cos_desc_t __global_1sided_direct_comm_desc;
extern cos_desc_t __global_1sided_direct_get_comm_desc;


/* set up internals */

#ifdef MAX_BUFS
#error "MAX_BUFS should not be defined yet"
#else
#define MAX_BUFS                ARMCI_MAX_BUFS
#endif

#ifdef MAX_SMALL_BUFS
#error "MAX_SMALL_BUFS should not be defined yet"
#else
#define MAX_SMALL_BUFS          ARMCI_MAX_SMALL_BUFS
#endif

#ifdef MSG_BUFLEN_DBL
#error "MSG_BUFLEN_DBL should not be defined yet"
#else
#define MSG_BUFLEN_DBL          ARMCI_BUF_SIZE
#endif


/* for buffers */

extern char **client_buf_ptrs;
#define BUF_ALLOCATE armci_portals_client_buf_allocate
//define BUF_EXTRA_FIELD_T comp_desc* 
//define INIT_SEND_BUF(_field,_snd,_rcv) _snd=1;_rcv=1;_field=NULL
#define GET_SEND_BUFFER _armci_buf_get
#define FREE_SEND_BUFFER _armci_buf_release

//define CLEAR_SEND_BUF_FIELD(_field,_snd,_rcv,_to,_op) if((_op==UNLOCK || _op==PUT || ARMCI_ACC(_op)) && _field!=NULL)x_buf_wait_ack((request_header_t *)((void **)&(_field)+1),((char *)&(_field)-sizeof(BUF_INFO_T)));_field=NULL;
//define TEST_SEND_BUF_FIELD(_field,_snd,_rcv,_to,_op,_ret)

#define CLEAR_SEND_BUF_FIELD(_field,_snd,_rcv,_to,_op)
#define TEST_SEND_BUF_FIELD(_field,_snd,_rcv,_to,_op,_ret)

#define COMPLETE_HANDLE _armci_buf_complete_nb_request

//define NB_CMPL_T comp_desc*
#if 0
#define ARMCI_NB_WAIT(_cntr) if(_cntr){\
        int rc;\
        if(nb_handle->tag)\
          if(nb_handle->tag==_cntr->tag)\
          rc = armci_client_complete(0,nb_handle->proc,nb_handle->tag,_cntr);\
} else{\
printf("\n%d:wait null ctr\n",armci_me);}
#endif
#define ARMCI_NB_WAIT(_cntr)



#endif
