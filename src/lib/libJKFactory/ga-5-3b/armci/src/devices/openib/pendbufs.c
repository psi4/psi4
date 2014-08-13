#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if defined(PEND_BUFS)

#include "pendbufs.h"
#include "armcip.h"
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif

#define DEBUG_SERVER 0

/*-------------------Attributes-------------------------*/

/**Attributes to control buffer count and sizes. Implement this way to
   hide the global variables, and provide get/set methods.*/

#define NUM_ATTRIBUTES 4
#define ATTRIB_IMMBUF_LEN 0
#define ATTRIB_IMMBUF_NUM 1
#define ATTRIB_PNDBUF_LEN 2
#define ATTRIB_PNDBUF_NUM 3

/** List of hidden attributes and their operations. 
 * @param attid IN Attribute id. Choose from the list above
 * @param gs IN Get(=0)/Set(=1)
 * @param v IN Value (used only when gs==1)
 * @return Value of the attribute on return
 */
static int att_ops(int attid, int gs, int v) {
  static not_first[NUM_ATTRIBUTES]; /*auto-init to zero*/
  static val[NUM_ATTRIBUTES];
  assert(attid>=0 && attid<NUM_ATTRIBUTES);
  assert(gs==0 || gs==1);
  if(!((gs && !not_first[attid]) || (!gs && not_first[attid]))) {
    printf("%d(%d) pausing...\n", armci_me, (int)getpid());
    pause();
  }
  if(gs) {
    not_first[attid]=1;
    val[attid]=v;
  }
  return val[attid];
}

int get_immbuf_len() { return att_ops(ATTRIB_IMMBUF_LEN,0,-1);}
void set_immbuf_len(int val) { att_ops(ATTRIB_IMMBUF_LEN,1,val);}
int get_immbuf_num() { return att_ops(ATTRIB_IMMBUF_NUM,0,-1);}
void set_immbuf_num(int val) { att_ops(ATTRIB_IMMBUF_NUM,1,val);}
int get_pendbuf_len() { return att_ops(ATTRIB_PNDBUF_LEN,0,-1);}
void set_pendbuf_len(int val) { att_ops(ATTRIB_PNDBUF_LEN,1,val);}
int get_pendbuf_num() { return att_ops(ATTRIB_PNDBUF_NUM,0,-1);}
void set_pendbuf_num(int val) { att_ops(ATTRIB_PNDBUF_NUM,1,val);}

#undef NUM_ATTRIBUTES
#undef ATTRIB_IMMBUF_LEN
#undef ATTRIB_IMMBUF_NUM
#undef ATTRIB_PNDBUF_LEN
#undef ATTRIB_PNDBUF_NUM

void armci_pbuf_init_buffer_env() {
  set_immbuf_len(IMM_BUF_LEN_DEFAULT);
  set_immbuf_num(IMM_BUF_NUM_DEFAULT);
  set_pendbuf_len(PENDING_BUF_LEN_DEFAULT);
  set_pendbuf_num(PENDING_BUF_NUM_DEFAULT);
}

/*---------------Waiting lists for ordering--------------------*/

/** Waiting procs list. This is a structure to maintain a list of
 *  processes with any waiting messages (immediate or non-immediate). 
 * Invariants: 
 * - !prev equivals !next (boolean type)
 * - !immbuf_wlist_head equivals !immbuf_wlist_tail (boolean type)
 * - !prev equivals !immbuf_wlist_head (boolean type)
 * - !immbuf_waitlist_head equivals (n_pending!=0) (boolean type)
 * - !order_head equivals !order_tail (boolean type)
 * In particular, note that order_head and order_tail are set and
 * reset independent of any immbufs waiting to be processed.
 */
typedef struct proc_waitlist_t {
  struct proc_waitlist_t *prev; /*<Previous proc in list of procs with
				  waiting immbufs (circular dlist)*/ 
  struct proc_waitlist_t *next; /*<Next proc in list of procs with
				  waiting immbufs*/  
  immbuf_t *immbuf_wlist_head; /**<Waiting immbufs for this proc - head*/
  immbuf_t *immbuf_wlist_tail; /**<Waiting immbufs for this proc - tail*/
  pendbuf_t *order_head,*order_tail; /**<List of pending buffers processing requests from this proc*/
/*   vapibuf_pend_t *waiting_on;   /\**<This pendbuf is for req from the */
/* 				   same proc*\/ */
  int n_pending; /**<#pending reqs fromt his proc*/
} proc_waitlist_t;


/** #pending buffers being used.
 */
static int _nPendBufsUsed=0;

/** The ordering semantics used. Need to provide accessors later if
 *  necessory.  
 */
static enum PendBufOrderingRule _pbufOrder = GET_GET_REORDER;

/*static*/ pendbuf_t *serv_pendbuf_arr;   /**<Array of pending buffers*/

static proc_waitlist_t *pbuf_ordering_plist_head; /**<Head of list of procs
						     with waiting immbufs*/
static proc_waitlist_t *pbuf_proc_list_info; /**<Array of structs for
						all the procs*/ 

/*-------------Pending states---------------------------*/

/**Status of pending buffers. Used in pending making progress and
 *  implementing teh corresponding finite-state-machine.
 */
enum PendBufStatus { 
  RECV_DATA_DONE=4, /**<Done recv-ing data. Used for PUT/ACCs*/
  SEND_DATA_DONE=5, /**<Done sending data. For GETs*/
  RECV_DSCR_DONE=6, /**<Done recv-ing descriptor. When
		       request_header_t+dscrlen>IMM_BUF_LEN*/ 
  INIT=7,             /**<Initial state for buffer*/
  RECV_DATA_PENDING=0,/**<Data recv posted and not completed*/
  SEND_DATA_PENDING=1,/**<Data send posted and not completed*/
  RECV_DSCR_PENDING=2 /**<Descriptor recv posted & not completed*/
};




/**Get a buffer to process an incoming request with non-immediate
 *  data. 
 * @return Pointer to available pending buffer
 */
static void* _armci_serv_pendbuf_getbuf(){
  int i;
  for(i=0; i<PENDING_BUF_NUM; i++) {
    if(serv_pendbuf_arr[i].avail) {
      serv_pendbuf_arr[i].status = INIT;
      serv_pendbuf_arr[i].avail=0;
      serv_pendbuf_arr[i].order_prev = NULL;
      serv_pendbuf_arr[i].order_next = NULL;
      serv_pendbuf_arr[i].commit_me = 0;
      /*fprintf(stderr, "%d:: getbuf returns idx=%d pbuf=%p\n",
	armci_me, i, &serv_pendbuf_arr[i]);*/
      _nPendBufsUsed += 1;
      assert(_nPendBufsUsed<=PENDING_BUF_NUM);
      ARMCI_PR_DBG("exit",0);
      return &serv_pendbuf_arr[i];        
    }
  }
  return NULL;
}

/**Assign a pending buffer for this immediate buffer. If a buffer is
 * available, it initializes the pending buffer appropriately.
 * @param vbuf IN Immediate buffer to be assign a pending buffer
 * @return Pointer to pending buffer, it assigned. NULL otherwise.
 */
static pendbuf_t* _armci_serv_pendbuf_assignbuf(immbuf_t *vbuf) {
  pendbuf_t *pbuf;
  const request_header_t *msginfo = (request_header_t *)vbuf->buf;
  proc_waitlist_t *info = &pbuf_proc_list_info[msginfo->from];
  assert(msginfo->tag.imm_msg == 0);
  pbuf = _armci_serv_pendbuf_getbuf();
  if(pbuf) {
    pbuf->status = INIT;
    pbuf->avail  = 0;
    pbuf->vbuf = vbuf;
    memcpy(pbuf->buf, vbuf->buf, sizeof(request_header_t)+msginfo->dscrlen);
/*     pbuf_proc_list_info[msginfo->from].waiting_on=pbuf; */
    pbuf->order_prev = info->order_tail;
    if(info->order_tail) info->order_tail->order_next = pbuf;
    info->order_tail = pbuf;
    if(!info->order_head) info->order_head = pbuf;
  }
  return pbuf;
}


/**Free a pending buffer
 * @param pbuf IN Pointer to Pending buffer to be freed
 * @return none
 */
static void _armci_serv_pendbuf_freebuf(pendbuf_t *pbuf){
  const request_header_t *msginfo = (request_header_t *)pbuf->buf;
  proc_waitlist_t *info = &pbuf_proc_list_info[msginfo->from];
  ARMCI_PR_DBG("enter",0);
  assert(pbuf != NULL);
  pbuf->avail=1;
  pbuf->status = -1;
  pbuf->vbuf = NULL;
/*   assert(info->waiting_on == pbuf); */
/*   info->waiting_on = NULL; */
  if(pbuf->order_prev) 
    pbuf->order_prev->order_next = pbuf->order_next;
  if(pbuf->order_next)
    pbuf->order_next->order_prev = pbuf->order_prev;
  if(info->order_head == pbuf) {
    assert(pbuf->order_prev == NULL);
    info->order_head = pbuf->order_next;
  }
  if(info->order_tail == pbuf) {
    assert(pbuf->order_next == NULL);
    info->order_tail = pbuf->order_prev;
  }
  pbuf->order_prev = pbuf->order_next = NULL; /*not necessary here*/
  
  _nPendBufsUsed -= 1;
  assert(_nPendBufsUsed>=0);
  ARMCI_PR_DBG("exit",0);
}

#if 0
/*Messages are processed in-place in immediate buffers or issued
  into pending buffers for progress in order (like
  ONE_PBUF_PER_MESG). This rule relaxes ONE_PBUF_PER_MESG by
  allowing ACCs to be processed in-place/issued
  without waiting for the prior reqs to complete*/
static int _can_progress_accnoorder(immbuf_t *vbuf) {
  const request_header_t *msginfo=(request_header_t*)vbuf->buf;
  const int proc = msginfo->from;
  const proc_waitlist_t *info = &pbuf_proc_list_info[proc];
  int i, nwaiting_on, nacc;
  pendbuf_t *ptr;
  
  assert(_pbufOrder == ACC_NO_ORDER);
  if(!IS_IMM_MSG(*msginfo) && _nPendBufsUsed==PENDING_BUF_NUM) {
    /*       printf("%d(s): op=%d from=%d datalen=%d waiting for pending buffers\n",armci_me,msginfo->operation,msginfo->from,msginfo->tag.data_len); */
    /*       fflush(stdout); */
    return 0; /*This buffer needs a free pending buffer*/
  }
  if(IS_IMM_MSG(*msginfo) && ARMCI_ACC(msginfo->operation)) {
    return 1;
  }
  if(info->immbuf_wlist_head && info->immbuf_wlist_head!=vbuf) {
    /*       printf("%d(s): op=%d from=%d datalen=%d not queue head\n",armci_me,msginfo->operation,msginfo->from,msginfo->tag.data_len); */
    /*       fflush(stdout); */
    return 0; /*in order issue*/
  }

  if(!ARMCI_ACC(msginfo->operation)) {
    if(info->order_head)
      return 0;
    return 1;
  }
  
  assert(ARMCI_ACC(msginfo->operation));
  for(ptr=info->order_head; ptr!=NULL; ptr=ptr->order_next) {
    request_header_t *m = (request_header_t *)ptr->buf;
    assert(m->from == msginfo->from);
    if(!ARMCI_ACC(m->operation)) 
      break;
  }
  if(ptr != NULL) 
    return 0;
  return 1;
}

static int _can_progress_putaccsplitorder(immbuf_t*vbuf) {
  if(!IS_IMM_MSG(*msginfo) && _nPendBufsUsed==PENDING_BUF_NUM) {
    return 0; /*This buffer needs a free pending buffer*/
  }
  if(IS_IMM_MSG(*msginfo) && ARMCI_ACC(msginfo->operation)) {
    return 1;
  }
  if(info->immbuf_wlist_head && info->immbuf_wlist_head!=vbuf) {
    return 0;
  }
  if(msginfo->operation!=PUT && !ARMCI_ACC(msginfo->operation)) {
    if(info->order_head)
      return 0;
    return 1;
  }
  if(IS_IMM_MSG(*msginfo) && info->order_head)
    return 0;
  return 1;  
}
#endif

/** Implement ordering between messages. This function needs to be
 * implemented in conjunction with @_armci_serv_pendbuf_promote to
 * ensure ordered processing of messages. 
 * @param vbuf IN Message in immediate buffer being checked 
 * @return 1 if the message can be progressed (either in-place or
 * after copying to a pending buffer). 0 therwise.
 */
static int _armci_serv_pendbuf_can_progress(immbuf_t *vbuf) {
  const request_header_t *msginfo=(request_header_t*)vbuf->buf;
  const int proc = msginfo->from;
  const proc_waitlist_t *info = &pbuf_proc_list_info[proc];

  if(_pbufOrder == ONE_PBUF_MESG) {
    /*Only one pending buffer used at any time*/
    if(_nPendBufsUsed>0) 
      return 0;
    return 1;
  }
  if(_pbufOrder == ONE_PBUF_MESG_PER_PROC) {
    /*Only one non-immediate mesg can be assigned to the pending
      buffers at any time*/
    if(info->order_head 
       || (info->immbuf_wlist_head && info->immbuf_wlist_head!=vbuf)) {
      return 0;/*other requests from this process remain*/
    }
    if(!IS_IMM_MSG(*msginfo) && _nPendBufsUsed==PENDING_BUF_NUM) {
      return 0; /*This buffer needs a free pending buffer*/
    }
    assert(info->n_pending == 0 || info->immbuf_wlist_head==vbuf);
    return 1;
  }
  if(_pbufOrder == ACC_NO_ORDER) {
    /*Messages are processed in-place in immediate buffers or issued
      into pending buffers for progress in order (like
      ONE_PBUF_PER_MESG). This rule relaxes ONE_PBUF_PER_MESG by
      allowing a sequence of ACCs to be processed in-place/issued
      without waiting for the prior ones to complete*/
    int i, nwaiting_on, nacc;
    pendbuf_t *ptr;
    if(!IS_IMM_MSG(*msginfo) && _nPendBufsUsed==PENDING_BUF_NUM) {
/*       printf("%d(s): op=%d from=%d datalen=%d waiting for pending buffers\n",armci_me,msginfo->operation,msginfo->from,msginfo->tag.data_len); */
/*       fflush(stdout); */
      return 0; /*This buffer needs a free pending buffer*/
    }
#if 1 /*commented for now: it does work*/
    if(IS_IMM_MSG(*msginfo) && ARMCI_ACC(msginfo->operation)) {
      return 1;
    }
#endif
    if(info->immbuf_wlist_head && info->immbuf_wlist_head!=vbuf) {
/*       printf("%d(s): op=%d from=%d datalen=%d not queue head\n",armci_me,msginfo->operation,msginfo->from,msginfo->tag.data_len); */
/*       fflush(stdout); */
      return 0; /*in order issue*/
    }

    if(!ARMCI_ACC(msginfo->operation)) {
      if(info->order_head)
	return 0;
      return 1;
    }

    assert(ARMCI_ACC(msginfo->operation));
    for(ptr=info->order_head; ptr!=NULL; ptr=ptr->order_next) {
      request_header_t *m = (request_header_t *)ptr->buf;
      assert(m->from == msginfo->from);
      if(!ARMCI_ACC(m->operation)) 
	break;
    }
    if(ptr != NULL) 
      return 0;
    return 1;
  }
  if(_pbufOrder == PUTACC_SPLIT_ORDER) {
    if(!IS_IMM_MSG(*msginfo) && _nPendBufsUsed==PENDING_BUF_NUM) {
      return 0; /*This buffer needs a free pending buffer*/
    }
    if(info->immbuf_wlist_head && info->immbuf_wlist_head!=vbuf) {
      return 0;
    }
    if(msginfo->operation!=PUT && !ARMCI_ACC(msginfo->operation)) {
      if(info->order_head)
	return 0;
      return 1;
    }
    if(IS_IMM_MSG(*msginfo) && ARMCI_ACC(msginfo->operation)) {
      return 1;
    }
    if(IS_IMM_MSG(*msginfo) && info->order_head)
      return 0;
    return 1;
  }
  if(_pbufOrder == GET_GET_REORDER) {
    if(!IS_IMM_MSG(*msginfo) && _nPendBufsUsed==PENDING_BUF_NUM) {
      return 0; /*This buffer needs a free pending buffer*/
    }
    if(IS_IMM_MSG(*msginfo) && ARMCI_ACC(msginfo->operation)) {
      return 1;
    }
    if(info->immbuf_wlist_head && info->immbuf_wlist_head!=vbuf) {
      return 0;
    }
    if(msginfo->operation!=PUT && !ARMCI_ACC(msginfo->operation)) {
      if(info->order_tail) {
	request_header_t *m=(request_header_t*)info->order_tail->buf;
	if(msginfo->operation==GET && m->operation == GET) {
/* 	  printf("%d: Get Get progressing\n", armci_me); */
	  return 1;
	}
	return 0;
      }
      return 1;
    }
    if(IS_IMM_MSG(*msginfo) && info->order_head)
      return 0;
    return 1;
  }
  armci_die("Unknown pbuf ordering rule",_pbufOrder);
  return 0;
}

/** Goes through the set of immediate buffers waiting to be processed
 *  and completed, and identifies a buffer that can be processed
 *  now. Removes it from the list and returns it. Promote also
 *  considers availability of pending buffers if need be.
 * @return Pointer to buffer that can be processed now. NULL if none exists.
 */
static immbuf_t* _armci_serv_pendbuf_promote() {
  immbuf_t *immbuf = NULL;
  proc_waitlist_t *info;

  ARMCI_PR_DBG("enter",0);

  assert(_nPendBufsUsed>=0);
  if(!pbuf_ordering_plist_head) {
    return NULL; /*nothing to promote*/
  }

  info = pbuf_ordering_plist_head;
  do {
    if(info->immbuf_wlist_head==NULL) {    
      printf("%d(s): Why is info->immbuf_wlist_head NULL\n", armci_me);
      fflush(stdout);
      pause();
    }
    assert(info->immbuf_wlist_head!=NULL);
    assert(info->n_pending>0);
    if(_armci_serv_pendbuf_can_progress(info->immbuf_wlist_head)) {
      immbuf = info->immbuf_wlist_head;
      info->immbuf_wlist_head = immbuf->immbuf_list_next;
      info->n_pending -= 1;
      immbuf->immbuf_list_next = NULL;
      if(!info->immbuf_wlist_head) {
         assert(info->immbuf_wlist_tail == immbuf);
         info->immbuf_wlist_tail = NULL;
	/*remove this proc from proc list*/
         info->prev->next = info->next;
         info->next->prev = info->prev;
         if(pbuf_ordering_plist_head == info) {
            pbuf_ordering_plist_head = (info->next==info)?NULL:info->next;
         }
         info->prev = info->next = NULL;
      }
      break; 
    }
    info = info->next;
  } while(info != pbuf_ordering_plist_head);

  if(DEBUG_SERVER) if(immbuf) {
    request_header_t *msginfo=(request_header_t*)immbuf->buf;
    printf("%d:: promoting a buffer immbuf=%p op=%d from=%d n_pending=%d\n", 
	   armci_me,immbuf,msginfo->operation,msginfo->from,info->n_pending);
    fflush(stdout);
  }
  ARMCI_PR_DBG("exit",0);
  return immbuf;
}

/** Enqueue a message. It could be an immediate message that cannot
 *  make progress or a non-immediate message that cannot make progress
 *  either due to ordering constraints or lack of pending buffers. 
 *  @param vbuf IN Immediate buffer to be enqueud
 *  @return Pending buffer into which the message was enqueued. NULL
 *  if no pending buffer was allocated (which is always the case for
 *  immediate messages)
 */
static pendbuf_t* _armci_serv_pendbuf_enqueue(immbuf_t *vbuf) {
  request_header_t *msginfo=(request_header_t *)vbuf->buf;
  int from = msginfo->from;
  pendbuf_t *pbuf;
  immbuf_t *ibuf;
  proc_waitlist_t *info = &pbuf_proc_list_info[msginfo->from];
  ARMCI_PR_DBG("enter",0);

/*   printf("%d: Entered serv_pbuf_enqueue\n", armci_me); */

  pbuf=NULL;
  if(msginfo->tag.imm_msg) {
    assert(!_armci_serv_pendbuf_can_progress(vbuf));
  }
  else if(_armci_serv_pendbuf_can_progress(vbuf)) {
    pbuf = _armci_serv_pendbuf_assignbuf(vbuf);
    assert(pbuf != NULL); /*can_progress() should ensure this*/
  }
  if(pbuf == NULL) {    
/*     printf("%d(s):: Enqueing op=%d imm=%d from %d. n_pending=%d\n", armci_me, msginfo->operation, msginfo->tag.imm_msg, msginfo->from,info->n_pending); */
/*     fflush(stdout); */
    vbuf->immbuf_list_next = NULL;
    assert(info->n_pending < IMM_BUF_NUM); /*How another message now?*/
    for(ibuf=info->immbuf_wlist_head;ibuf!=NULL;ibuf=ibuf->immbuf_list_next) {
       assert(vbuf != ibuf); /*enqueueing the same buffer twice*/
    }
    
    info->n_pending += 1;

    if(!info->immbuf_wlist_head) {
      assert(!info->immbuf_wlist_tail);
      assert(!info->prev && !info->next);
      /*insert proc into proc list*/
      if(!pbuf_ordering_plist_head) {
	pbuf_ordering_plist_head=info->next=info->prev=info;
      }
      else {
	info->next = pbuf_ordering_plist_head;
	info->prev = pbuf_ordering_plist_head->prev;
	pbuf_ordering_plist_head->prev->next = info;
	pbuf_ordering_plist_head->prev = info;
      }
    }
    /*insert vbuf into immbuf list for this proc*/
    if(info->immbuf_wlist_tail)
      info->immbuf_wlist_tail->immbuf_list_next=vbuf;
    info->immbuf_wlist_tail = vbuf;
    if(!info->immbuf_wlist_head)
      info->immbuf_wlist_head = vbuf;
  }
/*   printf("%d: Leaving serv_pbuf_enqueue\n", armci_me); */
  ARMCI_PR_DBG("exit",0);
  return pbuf;
}

/** Progress GET requests. 
 * @param pbuf IN Pending buffer containing the GET request
 * @return none
 */
static void _armci_serv_pendbuf_progress_get(pendbuf_t *pbuf) {
  int index = (pbuf - serv_pendbuf_arr);
  request_header_t *msginfo = (request_header_t *)pbuf->buf;
  void *buffer =((char *)(msginfo+1)+msginfo->dscrlen);
  int *status = &pbuf->status;  

  assert(sizeof(request_header_t)+msginfo->dscrlen+msginfo->datalen<PENDING_BUF_LEN);
  switch(*status) {
  case INIT:
    if(sizeof(request_header_t)+msginfo->dscrlen <= IMM_BUF_LEN) {
      /*Have the header and descriptor; go process*/
      armci_complete_pendbuf(pbuf);
      if(msginfo->pinned) {
	*status = SEND_DATA_DONE;
	_armci_serv_pendbuf_freebuf(pbuf);
      }
      else {
	*status = SEND_DATA_PENDING;
      }
    }
    else { /*Need to get rest of descriptor*/
      const int bytes = sizeof(request_header_t)+msginfo->dscrlen-IMM_BUF_LEN;
#warning "PEND_BUFS: Abusing msginfo->tag.ack_ptr for GETS with large descriptors!"
      assert(msginfo->tag.ack_ptr != NULL); /*sanity check. Should point to tag.ack on the client side*/
      void *lptr = ((char *)msginfo)+IMM_BUF_LEN;
      void *rptr = ((char *)msginfo->tag.ack_ptr) - (int)(&((request_header_t *)0)->tag.ack) + IMM_BUF_LEN;
/*       printf("%d(s):: GET getting rest of descriptor index=%d bytes=%d ptr=%p from=%d\n", */
/* 	     armci_me,index,bytes,rptr,msginfo->from); */
/*       fflush(stdout); */
      assert(IMM_BUF_LEN+bytes < PENDING_BUF_LEN);
      armci_pbuf_start_get(msginfo,rptr,lptr,bytes,msginfo->from,index);
      *status = RECV_DSCR_PENDING;
    }
    break;
  case RECV_DSCR_PENDING:
    armci_die("call_data_server should set status to RECV_DSCR_DONE before calling progress",*status);
    break;
  case SEND_DATA_PENDING:
    armci_die("call_data_server should set status to SEND_DATA_DONE before calling progress",*status);
    break;
  case RECV_DSCR_DONE:
/*     printf("%d(s):: GET. Done recving descriptor index=%d op=%d datalen=%d from=%d\n", */
/* 	   armci_me,index,msginfo->operation,msginfo->datalen,msginfo->from); */
/*     fflush(stdout); */
    armci_complete_pendbuf(pbuf);
    if(msginfo->pinned) {
      *status = SEND_DATA_DONE;
      _armci_serv_pendbuf_freebuf(pbuf);
    }
    else {
      *status = SEND_DATA_PENDING;
    }
    break;
  case SEND_DATA_DONE:
    _armci_serv_pendbuf_freebuf(pbuf);
    break;
  case RECV_DATA_PENDING:
  case RECV_DATA_DONE:
  default:
    armci_die("pendbuf_progress_get: invalid status", *status);
  }
}

/** Progress PUT/ACC requests.
 * @param pbuf IN Pending buffer containing the PUT/ACC request
 * @return none
 */
static void _armci_serv_pendbuf_progress_putacc(pendbuf_t *pbuf) {
  int index = (pbuf - serv_pendbuf_arr);
  request_header_t *msginfo = (request_header_t *)pbuf->buf;
  void *buffer =((char *)(msginfo+1))+msginfo->dscrlen;
  int *status = &pbuf->status;  

  assert(msginfo->operation==PUT || ARMCI_ACC(msginfo->operation));
  assert(sizeof(request_header_t)+msginfo->dscrlen+msginfo->datalen<PENDING_BUF_LEN);
  switch(*status) {
  case INIT:
/*     printf("%d(s): progressing new msg. index=%d op=%d from=%d\n", armci_me,index,msginfo->operation,msginfo->from); */
/*     fflush(stdout); */
    if(sizeof(request_header_t)+msginfo->dscrlen <= IMM_BUF_LEN) {
      /*Have the header and descriptor; go process*/
      assert(sizeof(request_header_t)+msginfo->dscrlen+msginfo->tag.data_len < PENDING_BUF_LEN);
      armci_pbuf_start_get(msginfo,msginfo->tag.data_ptr,buffer,msginfo->tag.data_len,
			   msginfo->from, index);
      /*       printf("%d(s): PUT/ACC getting data. pbuf_num=%d data_ptr=%p data_len=%d bytes=%d\n", armci_me,index,msginfo->tag.data_ptr, msginfo->tag.data_len,msginfo->bytes); */
      *status = RECV_DATA_PENDING;
    }
    else { /*Need to get rest of descriptor*/
      const int bytes = sizeof(request_header_t)+msginfo->dscrlen-IMM_BUF_LEN;
#warning "PEND_BUFS: Abusing msginfo->tag.ack_ptr for GETS with large descriptors!"
      assert(msginfo->tag.ack_ptr != NULL); /*sanity check. Should point to tag.ack on the client side*/
      void *lptr = ((char *)msginfo)+IMM_BUF_LEN;
      void *rptr = ((char *)msginfo->tag.ack_ptr) - (int)(&((request_header_t *)0)->tag.ack) + IMM_BUF_LEN;
/*       printf("%d(s):: PUT getting rest of descriptor index=%d bytes=%d ptr=%p from=%d\n", */
/* 	     armci_me,index,bytes,rptr,msginfo->from); */
/*       fflush(stdout); */
      assert(IMM_BUF_LEN+bytes < PENDING_BUF_LEN);
      armci_pbuf_start_get(msginfo,rptr,lptr,bytes,msginfo->from,index);
      *status = RECV_DSCR_PENDING;      
    }
    break;
  case RECV_DSCR_PENDING:
    armci_die("call_data_server should set status to RECV_DSCR_DONE before calling progress",*status);
    break;
  case RECV_DATA_PENDING:
    armci_die("call_data_server should set status to RECV_DONE before calling progress",*status);
    break;
  case RECV_DSCR_DONE:
      assert(sizeof(request_header_t)+msginfo->dscrlen+msginfo->tag.data_len < PENDING_BUF_LEN);
      armci_pbuf_start_get(msginfo,msginfo->tag.data_ptr,buffer,msginfo->tag.data_len,
			   msginfo->from, index);
/*     printf("%d(s): PUT/ACC getting data. pbuf_num=%d data_ptr=%p data_len=%d bytes=%d\n", armci_me,index,msginfo->tag.data_ptr, msginfo->tag.data_len,msginfo->bytes); */
    *status = RECV_DATA_PENDING;
    break;
  case RECV_DATA_DONE:
/*     printf("%d(s):: Done PUT/ACC with buf index=%d op=%d datalen=%d from=%d\n", */
/* 	   armci_me,index,msginfo->operation,msginfo->datalen,msginfo->from); */
/*     fflush(stdout); */
    if(msginfo->operation == PUT && pbuf->order_prev!=NULL) {
      assert(pbuf->commit_me == 0); /*Why called so many times in thie
				      state?*/
      pbuf->commit_me = 1;
      break;
    }
    pbuf->commit_me = 0;
    armci_complete_pendbuf(pbuf);
    _armci_serv_pendbuf_freebuf(pbuf);
    break;
  case SEND_DATA_PENDING:
  case SEND_DATA_DONE:
  default:
    armci_die("pendbuf_progress_putacc: invalid status", *status);
  }
}


/** Make progress on processing a pending buffer. This function, also
 * ensures any other waiting messages get processed if they can
 * be. Thus, progress and eventual termination is guaranteed by this
 * function. 
 * @param _pbuf IN Pending buffer to make progress on
 * @return none
 */
static void _armci_serv_pendbuf_progress(pendbuf_t *_pbuf){
  request_header_t *msginfo = (request_header_t *)_pbuf->buf;
  immbuf_t *vbuf = _pbuf->vbuf;
  pendbuf_t *pbuf = _pbuf;

  assert(pbuf->vbuf!=NULL);
  do {
    if(vbuf && !IS_IMM_MSG(*msginfo)) { assert(pbuf->vbuf == vbuf); }
/*     printf("%d(s):: progressing op=%d imm=%d from=%d datalen=%d pbuf=%p vbuf=%p n_pending=%d\n", armci_me, */
/* 	   msginfo->operation,msginfo->tag.imm_msg,msginfo->from,msginfo->datalen, pbuf,vbuf,pbuf_proc_list_info[msginfo->from].n_pending); */
/*     fflush(stdout); */
    if(IS_IMM_MSG(*msginfo)) {
      armci_complete_immbuf(vbuf);
    }
    else { /*non-immediate message*/
      proc_waitlist_t* info = &pbuf_proc_list_info[msginfo->from];
      
      do {
	assert(pbuf->vbuf == vbuf);
	if(msginfo->operation == PUT || ARMCI_ACC(msginfo->operation)) {
	  _armci_serv_pendbuf_progress_putacc(pbuf);
	}
	else if (msginfo->operation == GET) {
	  _armci_serv_pendbuf_progress_get(pbuf);
	}
	else {
	  armci_die("pending buffer processing for this op not yet implemented", msginfo->operation);
	}
	pbuf = info->order_head;
	vbuf = pbuf ? pbuf->vbuf : NULL;
      } while(info->order_head && info->order_head->commit_me);
    }
/*     sleep(2); */
    vbuf = _armci_serv_pendbuf_promote();
    if(vbuf) {
      msginfo = (request_header_t *)vbuf->buf;
      if(!msginfo->tag.imm_msg) {
	pbuf = _armci_serv_pendbuf_assignbuf(vbuf);
	assert(pbuf != NULL);
      }
    }
  } while(vbuf != NULL);
}



/*----------------External functions--------------------*/


/** Initialize array of pending buffers
 *  @return none
 */
void armci_pendbuf_init() {
  int i;

  ARMCI_PR_DBG("enter",0);

/*   bzero(serv_pendbuf_arr, sizeof(pendbuf_t)*PENDING_BUF_NUM); */
  for(i=0; i<PENDING_BUF_NUM; i++) {
    pendbuf_t *pbuf = &serv_pendbuf_arr[i];
    char *buf = pbuf->buf;
    bzero(pbuf, sizeof(pendbuf_t));
    pbuf->buf = buf;
    pbuf->avail=1;
  }
  
  pbuf_ordering_plist_head=NULL;
  pbuf_proc_list_info = (proc_waitlist_t *)malloc(sizeof(proc_waitlist_t)*armci_nproc);
  assert(pbuf_proc_list_info != NULL);
  bzero(pbuf_proc_list_info, sizeof(proc_waitlist_t)*armci_nproc);
  ARMCI_PR_DBG("exit",0);
}

void armci_pendbuf_service_req(immbuf_t *immbuf) {
  pendbuf_t *pbuf;
  request_header_t *msginfo=(request_header_t*)immbuf->buf;

  if(IS_IMM_MSG(*msginfo) && _armci_serv_pendbuf_can_progress(immbuf)) {
    /* 	   printf("%d: msg vbuf=%p op=%d from=%d imm=%d datalen=%d bytes=%d data_ptr=%p can progress. Completing it now!\n", */
    /* 		  armci_me, vbuf, msginfo->operation, msginfo->from, msginfo->tag.imm_msg,msginfo->datalen,msginfo->bytes,msginfo->tag.data_ptr); */
    /* 	   fflush(stdout); */
    armci_complete_immbuf(immbuf);
  }
  else if(pbuf = _armci_serv_pendbuf_enqueue(immbuf)) {
    /* 	   printf("%d: msg vbuf=%p op=%d from=%d imm=%d datalen=%d bytes=%d data_ptr=%p got pending buf. Progressing it!\n", */
    /* 		  armci_me, vbuf, msginfo->operation, msginfo->from, msginfo->tag.imm_msg,msginfo->datalen,msginfo->bytes,msginfo->tag.data_ptr); */
    /* 	   fflush(stdout); */
    _armci_serv_pendbuf_progress(pbuf);
  }
  else {	  
    /*  	   printf("%d: msg vbuf=%p op=%d from=%d imm=%d datalen=%d bytes=%d data_ptr=%p in waitlist!\n", armci_me, vbuf, msginfo->operation, msginfo->from, msginfo->tag.imm_msg,msginfo->datalen,msginfo->bytes,msginfo->tag.data_ptr); */
    /* 	   fflush(stdout); */
  } 
}

/**Network layer reporting to split buffers code that a put completed
 *  on a pending buffer.
 * @param pbufid IN Pending buffer id (specified when starting a
 * put).
 * @return void
 */
void armci_pendbuf_done_put(int pbufid) {
  assert(pbufid>=0 && pbufid<PENDING_BUF_NUM);
  assert(serv_pendbuf_arr[pbufid].status == SEND_DATA_PENDING);
  serv_pendbuf_arr[pbufid].status = SEND_DATA_DONE;
  _armci_serv_pendbuf_progress(&serv_pendbuf_arr[pbufid]);
}

/**Network layer reporting to split buffers code that a get completed
 *  on a pending buffer.
 * @param pbufid IN Pending buffer id (specified when starting a
 * get).
 * @return void
 */
void armci_pendbuf_done_get(int pbufid) {
  int done_status;
  pendbuf_t *pbuf;
  assert(pbufid>=0 && pbufid<PENDING_BUF_NUM);
  pbuf = &serv_pendbuf_arr[pbufid];

  switch(pbuf->status) {
  case RECV_DSCR_PENDING:
    pbuf->status = RECV_DSCR_DONE;
    break;
  case RECV_DATA_PENDING:
    pbuf->status = RECV_DATA_DONE;
    break;
  default:
    armci_die("Reporting get done on buf with inappropriate status",pbufid);
  }
  _armci_serv_pendbuf_progress(pbuf);
}


#endif /*PEND_BUFS*/

