#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: buffers.c,v 1.29.6.9 2007-07-02 05:16:50 d3p687 Exp $    **/

#define SIXTYFOUR 64
#define DEBUG_  0
#define DEBUG2_ 0
#define EXTRA_ERR_CHECK

/**********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "armcip.h"
#include "request.h"
#ifdef WIN32
#  include <windows.h>
   typedef unsigned long ssize_t;
#else
#  include <unistd.h>
#endif

#   define EQ_TAGS(a_, b_) !memcmp(&(a_), &(b_), sizeof(a_))

#define ALIGN64ADD(buf) (SIXTYFOUR-(((ssize_t)(buf))%SIXTYFOUR))
/* the following symbols should be defined if needed in protocol specific
   header file:  BUF_EXTRA_FIELD, BUF_ALLOCATE
*/

#ifndef BUF_ALLOCATE
#   define BUF_ALLOCATE malloc
#endif
#if defined PORTALS
# define SMALL_BUF_LEN PORTALS_SMALL_BUF_SIZE
#else
# if defined(SERV_QUEUE)
#  define SMALL_BUF_LEN 4096
# else
#  define SMALL_BUF_LEN 2048
# endif
#endif

#ifndef MSG_BUFLEN_SMALL
# define MSG_BUFLEN_SMALL (MSG_BUFLEN >>0)
#endif

#define LEFT_GUARD  11.11e11
#define RIGHT_GUARD 22.22e22
#define CLEAR_TABLE_SLOT(idx) *((int*)(_armci_buf_state->table+(idx))) =0

#ifndef BUF_NET_INIT
#define BUF_NET_INIT(x,xX,Xx)
#endif
_buf_ackresp_t *_buf_ackresp_first,*_buf_ackresp_cur;
/* we allow multiple buffers (up to 15) per single request
 * adjacent buffers can be coalesced into a large one
 */
typedef struct {
  int op;     /* pending operation code */
  int snd;    /* if 1 then buffer is used for sending request */
  int rcv;    /* if 1 then buffer is used for receiving data */
  int async;  /* if 1 then request is nonblocking */
  int first;  /* id of the 1st buffer in the set in same request */
  int count;  /* count is not used and is always 1 (or 0???) */
  /*unsigned int count:4;  \* how many buffers used for this request 8 possible */
  int busy;   /* if 1 buffer is used and cannot be completed */
  int cmpl;   /* set to 1 if buffer was completed and can be released */
  int to;    /* serv/proc to which request was sent, 8172 possible */
}buf_state_t;


#ifndef BUFID_PAD_T
#define BUFID_PAD_T BUF_INFO_T
#endif

/* message send buffer data structure */
typedef struct {
  BUF_INFO_T id;
# ifdef BUF_EXTRA_FIELD_T
        BUF_EXTRA_FIELD_T field;
# endif
  char buffer[MSG_BUFLEN_SMALL];
} buf_ext_t;

/* message send buffer data structure */
typedef struct {
  BUF_INFO_T id;
# ifdef BUF_EXTRA_FIELD_T
        BUF_EXTRA_FIELD_T field;
# endif
  char buffer[SMALL_BUF_LEN];
} buf_smext_t;

/* we keep table and buffer pointer together for better locality */
typedef struct {
  double left_guard;        /* stamp to verify if array was corrupted */
  buf_state_t table[MAX_BUFS+MAX_SMALL_BUFS]; /*array with state of buffer */
  buf_ext_t *buf;           /* address of buffer pool */
  buf_smext_t *smallbuf;      /* address of the large buffer pool */
  int avail;
  int smavail;
  int pad;
  double right_guard;       /* stamp to verify if array was corrupted */

  unsigned buf_bitmap;  /* bitmaps to track available buffers: */
  unsigned smbuf_bitmap;/*  1 - available, 0 - not available   */
} reqbuf_pool_t;

#ifndef  BUF_EXTRA_FIELD_T
#  define        SIZE_BUF_EXTRA_FIELD 0 
#  define BUF_TO_EBUF(buf) (buf_ext_t*)(((char*)buf) - sizeof(BUFID_PAD_T) -\
                                      SIZE_BUF_EXTRA_FIELD)
#  define BUF_TO_SMEBUF(buf) (buf_smext_t*)(((char*)buf)- sizeof(BUFID_PAD_T) -\
                                      SIZE_BUF_EXTRA_FIELD)
#else
#  define BUF_TO_EBUF(buf) (buf_ext_t*)(((char*)buf) - sizeof(BUFID_PAD_T) -\
				      sizeof(BUF_EXTRA_FIELD_T))
#  define BUF_TO_SMEBUF(buf) (buf_smext_t*)(((char*)buf)- sizeof(BUFID_PAD_T) -\
				      sizeof(BUF_EXTRA_FIELD_T))
#endif

#define BUF_TO_BUFINDEX(buf) (BUF_TO_EBUF((buf)))->id.bufid
#define BUF_TO_SMBUFINDEX(buf) (BUF_TO_SMEBUF((buf)))->id.bufid


buf_ext_t *_armci_buffers;        /* these are the actual buffers */
buf_smext_t *_armci_smbuffers;    /* no, these are the actual buffers */
reqbuf_pool_t* _armci_buf_state;  /* array that describes state of each buf */

extern active_socks_t *_armci_active_socks;

/* returns bufinfo, given bufid */
INLINE BUF_INFO_T *_armci_id_to_bufinfo(int bufid) {
  if (bufid < 0 || bufid >= (MAX_BUFS+MAX_SMALL_BUFS))
      armci_die2("_armci_id_to_bufinfo: bad id",bufid,MAX_BUFS);

  return bufid < MAX_BUFS ? &(_armci_buf_state->buf[bufid].id) :
                            &(_armci_buf_state->smallbuf[bufid-MAX_BUFS].id);
}



/*\ we allocate alligned buffer space
 *  this operation can be implemented in platform specific files
\*/ 
void _armci_buf_init()
{
char *tmp;
int  extra=0;
int smallbuf_size = sizeof(buf_smext_t)*(MAX_SMALL_BUFS);
 //  tmp = (char *) BUF_ALLOCATE((MAX_BUFS*sizeof(buf_ext_t) + 64 + smallbuf_size + 64));
     tmp = (char *) malloc((MAX_BUFS*sizeof(buf_ext_t) + 64 + smallbuf_size + 64));
     bzero(tmp,MAX_BUFS*sizeof(buf_ext_t) + 64 + smallbuf_size + 64);
     extra= ALIGN64ADD(tmp);

     _armci_buffers = (buf_ext_t *) (tmp + extra); 

     tmp = (char *)(_armci_buffers + MAX_BUFS);
     extra = ALIGN64ADD(tmp);
     _armci_smbuffers = (buf_smext_t *) (tmp + extra); 
     

     if(DEBUG2_){
	printf("%d:armci_init_bufs: pointer %p, before align ptr=%p bufptr=%p end of region is %p  size=%d extra=%d\n",
               armci_me,_armci_buffers,tmp,_armci_buffers->buffer,(_armci_buffers+MAX_BUFS),
               MAX_BUFS*sizeof(buf_ext_t),extra);
	fflush(stdout);
     }

     /* now allocate state array */
     tmp  = malloc(sizeof(reqbuf_pool_t) + 64);
     bzero(tmp,sizeof(reqbuf_pool_t) + 64);
     if(!tmp)armci_die("_armci_buf_init calloc failed",0);
     extra= ALIGN64ADD(tmp);
     _armci_buf_state = (reqbuf_pool_t*)(tmp + extra); 

     /* initialize it */
     _armci_buf_state->left_guard  = LEFT_GUARD;
     _armci_buf_state->right_guard = RIGHT_GUARD;
     _armci_buf_state->avail =0;
     _armci_buf_state->smavail =MAX_BUFS;
     _armci_buf_state->buf = _armci_buffers;
     _armci_buf_state->smallbuf = _armci_smbuffers;

     _buf_ackresp_first=_buf_ackresp_cur=NULL;
     
     if(BUF_TO_EBUF(_armci_buf_state->buf[0].buffer)!=_armci_buf_state->buf)
        armci_die("buffers.c, internal structure alignment problem",0);
}


/*\ convert buffer pointer to index (in state array)
\*/
int _armci_buf_to_index(void *buf)
{
int index;
char *ptr = (char*)buf;

   if(DEBUG2_){
     printf("%d: in _armci_buf_to_index %p\n",armci_me, buf);
     fflush(stdout);
   }
   if(buf > (void *)_armci_buffers && buf < (void *)(_armci_buffers+MAX_BUFS)){
      index = BUF_TO_BUFINDEX(ptr);
      if((index >= MAX_BUFS)|| (index<0)) 
        armci_die2("armci_buf_to_index: bad index:",index,MAX_BUFS);
      return(index);
   }
   else if(buf > (void *)_armci_smbuffers && buf < (void *)(_armci_smbuffers+MAX_SMALL_BUFS)){
      index = BUF_TO_SMBUFINDEX(ptr);
      if((index >= MAX_BUFS+MAX_SMALL_BUFS)|| (index<MAX_BUFS)) 
        armci_die2("armci_buf_to_ind:indexwrong",index,MAX_BUFS+MAX_SMALL_BUFS);
      return(index);
   } 
   else {
        armci_die("armci_buf_to_index: bad pointer",0);
        return(0);
   }
}


void x_buf_send_complete(void *buf)
{
int index = _armci_buf_to_index(buf);
buf_state_t *buf_state = _armci_buf_state->table + index;
    ARMCI_PR_DBG("enter",0);
    if(index>=MAX_BUFS){
    int relidx;
      relidx = index-MAX_BUFS;
      CLEAR_SEND_BUF_FIELD(_armci_buf_state->smallbuf[relidx].field,buf_state->snd,buf_state->rcv,buf_state->to,buf_state->op);
    }
    else 
      CLEAR_SEND_BUF_FIELD(_armci_buf_state->buf[index].field,buf_state->snd,buf_state->rcv,buf_state->to,buf_state->op);
    ARMCI_PR_DBG("exit",0);
}

/*\  complete outstanding operation that uses the specified buffer
\*/
void _armci_buf_complete_index(int idx, int called)
{
int count;
buf_state_t *buf_state = _armci_buf_state->table +idx;
portals_ds_req_t *req  = NULL; 

    count = buf_state->count;
    if(DEBUG_ || 0) {
       printf("%d:buf_complete_index:%d op=%d first=%d count=%d called=%d\n",
              armci_me,idx,buf_state->op,buf_state->first,buf_state->count,
              called); 
       fflush(stdout);
    }

    if(buf_state->first != (unsigned int)idx){ 
      armci_die2("complete_buf_index:inconsistent Index:",idx,buf_state->first);
    }

    /* need to call platform specific function */
    if(idx>=MAX_BUFS){
      int relidx,rr;
      relidx = idx-MAX_BUFS; 
    //printf("\n%d:in clear idx=%d %d",armci_me,idx,_armci_buf_state->smallbuf[relidx].id.tag);fflush(stdout);    
   /* ------------------------------------------------------------------------------------------- *\
      active buffers need to be completed
   \* ------------------------------------------------------------------------------------------- */
    # ifdef PORTALS
      req = &_armci_buf_state->smallbuf[relidx].id.ar.req;
      if(req->active) {
     //  printf("%s [cp buf_complete_index] waiting on request %p\n",Portals_ID(),req);
         portalsWaitOnRequest(req);
     //  printf("%s [cp buf_complete_index] request %p completed\n",Portals_ID(),req);
      } else {
     //  printf("%s [cp buf_complete_index] request %p already completed\n",Portals_ID(),req);
      } 

    # else

      if(_armci_buf_state->smallbuf[relidx].id.tag && (_armci_buf_state->smallbuf[relidx].field)->tag>0) { 
        printf("%s [cp] calling armci_client_complete\n",Portals_ID()); 
        rr=armci_client_complete(0,buf_state->to,_armci_buf_state->smallbuf[relidx].id.tag,_armci_buf_state->smallbuf[relidx].field);
      } 
      CLEAR_SEND_BUF_FIELD(_armci_buf_state->smallbuf[relidx].field,buf_state->snd,buf_state->rcv,buf_state->to,buf_state->op);

    # endif

      /*later, we might just need to do this for all operations, not just get*/
    # ifdef PORTALS_ALLOW_NBGETS
      if(_armci_buf_state->smallbuf[relidx].id.tag!=0 &&(buf_state->op == GET)){
        armci_complete_req_buf(&(_armci_buf_state->smallbuf[relidx].id),
                                _armci_buf_state->smallbuf[relidx].buffer);
      }
    # endif
      _armci_buf_state->smallbuf[relidx].id.tag=0;
    }
    else {
      int rr;

   /* ------------------------------------------------------------------------------------------- *\
      active buffers need to be completed
   \* ------------------------------------------------------------------------------------------- */
    # ifdef PORTALS
      req = &_armci_buf_state->buf[idx].id.ar.req;
      if(req->active) portalsWaitOnRequest(req);
   
    # else

      if(_armci_buf_state->buf[idx].id.tag  && (_armci_buf_state->buf[idx].field)->tag>0 )
        rr=armci_client_complete(0,buf_state->to,_armci_buf_state->buf[idx].id.tag,_armci_buf_state->buf[idx].field);
      CLEAR_SEND_BUF_FIELD(_armci_buf_state->buf[idx].field,buf_state->snd,buf_state->rcv,buf_state->to,buf_state->op);
    //printf("\n%d:in clear large idx=%d %d",armci_me,idx,_armci_buf_state->buf[idx].id.tag);fflush(stdout);     
    # endif

      /*later, we might just need to do this for all operations, not just get*/
    # ifdef PORTALS_ALLOW_NBGETS
      if(_armci_buf_state->buf[idx].id.tag!=0 &&(buf_state->op == GET)){
        armci_complete_req_buf(&(_armci_buf_state->buf[idx].id),
                                _armci_buf_state->buf[idx].buffer);
      }
    # endif
      _armci_buf_state->buf[idx].id.tag=0;
    }
    /* clear table slots for all the buffers in the set for this request */
    for(; count; count--, buf_state++) *(int*)buf_state = 0;
}


/*\  test outstanding operation that uses the specified buffer for complete
 *   It is important not to change the state of the buffer, the buffer has
 *   to remain as it was, only completion has to be indicated
\*/
int _armci_buf_test_index(int idx, int called)
{
int count,retval=0;
buf_state_t *buf_state = _armci_buf_state->table +idx;
    count = buf_state->count;
    if(DEBUG_ ){
       printf("%d:buf_test_index:%d op=%d first=%d count=%d called=%d\n",
              armci_me,idx,buf_state->op,buf_state->first,buf_state->count,
              called); 
       fflush(stdout);
    }
    if(buf_state->first != (unsigned int)idx){ 
      armci_die2("_buf_test_index:inconsistent index:",idx,buf_state->first);
    }
#   ifdef BUF_EXTRA_FIELD_T
    /* need to call platform specific function */
    if(idx>=MAX_BUFS){
       int relidx;
       relidx = idx-MAX_BUFS; 
       /*printf("\n%d:relidx=%d \n",armci_me,relidx);fflush(stdout);*/
       TEST_SEND_BUF_FIELD(_armci_buf_state->smallbuf[relidx].field,buf_state->snd,buf_state->rcv,buf_state->to,buf_state->op,&retval);

    }
    else {
       TEST_SEND_BUF_FIELD(_armci_buf_state->buf[idx].field,buf_state->snd,buf_state->rcv,buf_state->to,buf_state->op,&retval);
    }
#   endif
    if(DEBUG_ ){
       printf("%d:buf_test_index:%d op=%d first=%d count=%d called=%d ret=%d\n",
              armci_me,idx,buf_state->op,buf_state->first,buf_state->count,
              called,retval); 
       fflush(stdout);
    }
    return(retval);
}

/**
an addition to the below operation to allow for multiple outstanding operations
per server node
*/
void _armci_buf_ensure_pend_outstanding_op_per_node(void *buf, int node)
{
int i;
int index =_armci_buf_to_index(buf); 
int this = _armci_buf_state->table[index].first;
int nfirst, nlast;
void _armci_buf_release_index(int i);
int buf_pend_count=0;
int changeid=0;
    nfirst=armci_clus_info[node].master;
    nlast = nfirst+armci_clus_info[node].nslave-1;
    if(_armci_buf_state->table[index].to<0){
      _armci_buf_state->table[index].to = 0-1e6-_armci_buf_state->table[index].to;
      changeid=1;
    }

    if((_armci_buf_state->table[index].to<(unsigned int) nfirst) || 
      (_armci_buf_state->table[index].to>(unsigned int) nlast))
        armci_die2("_armci_buf_ensure_pend_outstanding_op_per_node: bad to",node,
                        (int)_armci_buf_state->table[index].to);

    buf_pend_count=0;
    for(i=0;i<MAX_BUFS+MAX_SMALL_BUFS;i++){
      buf_state_t *buf_state = _armci_buf_state->table +i;
      if((buf_state->to >= nfirst) && (buf_state->to<= (unsigned int) nlast))
        if( (buf_state->first != (unsigned int) this) && (buf_state->first==(unsigned int) i) && buf_state->op){
          buf_pend_count++;
          if(buf_pend_count == NUM_SERV_BUFS){
            _armci_buf_complete_index(i,0);
	          _armci_buf_release_index(i);
            break;
          }
        }
    }
    if(changeid)_armci_buf_state->table[index].to = 0-1e6-_armci_buf_state->table[index].to;
}

/*\ make sure that there are no other pending operations to that smp node
 *  this operation is called from platforms specific routine that sends
 *  request
 *  we could have accomplished the same in armci_buf_get but as Vinod
 *  is pointing out, it is better to delay completing outstanding
 *  calls to overlap memcpy for the current buffer with communication
\*/
void _armci_buf_ensure_one_outstanding_op_per_node(void *buf, int node)
{
    int i;
    int index =_armci_buf_to_index(buf); 
    int this = _armci_buf_state->table[index].first;
    int nfirst, nlast;
    void _armci_buf_release_index(int i);

    nfirst=armci_clus_info[node].master;
    nlast = nfirst+armci_clus_info[node].nslave-1;
    if((_armci_buf_state->table[index].to<(unsigned int) nfirst) || 
       (_armci_buf_state->table[index].to>(unsigned int) nlast))
        armci_die2("_armci_buf_ensure_one_outstanding_op_per_node: bad to",node,
                (int)_armci_buf_state->table[index].to);

    for(i=0;i<MAX_BUFS+MAX_SMALL_BUFS;i++){
        buf_state_t *buf_state = _armci_buf_state->table +i;
        if((buf_state->to >= nfirst) && (buf_state->to<= (unsigned int) nlast)) {
          if((buf_state->first != (unsigned int) this)&&(buf_state->first==(unsigned int) i) && buf_state->op) {
	    _armci_buf_complete_index(i,0);
	    _armci_buf_release_index(i);
	  }
	}
    }
}

/*\ same as above but for process
\*/
void _armci_buf_ensure_one_outstanding_op_per_proc(void *buf, int proc)
{
    int i;
    int index = _armci_buf_to_index(buf); 
    int this = _armci_buf_state->table[index].first;
    void _armci_buf_release_index(int i);

    if(_armci_buf_state->table[index].to !=(unsigned int)  proc )
       armci_die2("_armci_buf_ensure_one_outstanding_op_per_proc: bad to", proc,
                (int)_armci_buf_state->table[index].to);

    for(i=0;i<MAX_BUFS+MAX_SMALL_BUFS;i++){
        buf_state_t *buf_state = _armci_buf_state->table +i;
        if(buf_state->to == (unsigned int) proc) {
          if((buf_state->first != (unsigned int) this)&&(buf_state->first==(unsigned int) i) && buf_state->op) {
	    _armci_buf_complete_index(i,0);
	    _armci_buf_release_index(i);
	  }
	}
    }
}


#define HISTORY__ 
#ifdef HISTORY
typedef struct{ int size; int op; int count; int id; } history_t;
history_t history[100];
int h=0;

void print_history()
{
int i;
    fflush(stdout);
    printf("%d records\n",h);
    for(i=0; i<h;i++) printf("size=%d id=%d ptr=%p count=%d op=%d\n",
        history[i].size, history[i].id,
       _armci_buf_state->buf[history[i].id].buffer, history[i].count,
        history[i].op);

    fflush(stdout);
}
#endif

/*\  call corresponding to GET_SEND_BUF
\*/
char *_armci_buf_get_small(int size, int operation, int to)
{
int avail=_armci_buf_state->smavail,i;
_buf_ackresp_t *ar;
    if(_armci_buf_state->table[avail].op || 
       _armci_buf_state->table[avail].first ||
       _armci_buf_state->smallbuf[avail-MAX_BUFS].id.ar.req.active) {
       
       for(i=MAX_BUFS;i<MAX_BUFS+MAX_SMALL_BUFS;i++){
          if(!_armci_buf_state->table[i].op &&
             !_armci_buf_state->table[i].first &&
             !_armci_buf_state->smallbuf[i-MAX_BUFS].id.ar.req.active)
            break;
       }
       if(i<(MAX_SMALL_BUFS+MAX_BUFS))avail = i;
       else {
          _armci_buf_complete_index(avail,1);
       }
    }
    _armci_buf_state->table[avail].op = operation;
    _armci_buf_state->table[avail].to = to;
    _armci_buf_state->table[avail].count=  1;
    _armci_buf_state->table[avail].first = avail;
    _armci_buf_state->smallbuf[avail-MAX_BUFS].id.tag=0;
    _armci_buf_state->smallbuf[avail-MAX_BUFS].id.bufid= avail; 
    _armci_buf_state->smallbuf[avail-MAX_BUFS].id.protocol=0;
    ar=&_armci_buf_state->smallbuf[avail-MAX_BUFS].id.ar;
    assert(ar->val==0);assert(ar->next==NULL);assert(ar->previous==NULL);
  # ifdef PORTALS
    assert(ar->req.active == 0);
  # endif
    ar->req.active = 1;
    if(_buf_ackresp_cur!=NULL)
      _buf_ackresp_cur->next=ar;
    if(_buf_ackresp_first==NULL)
      _buf_ackresp_first=ar;
    ar->previous=_buf_ackresp_cur;
    ar->next=NULL;
    _buf_ackresp_cur=ar;

    if(DEBUG_ || 0) {
      printf("%d:buf_get_sm1:size=%d max=%d got %d ptr=%p op=%d to=%d count=%d first=%d\n",
             armci_me,size,SMALL_BUF_LEN,avail,
	     _armci_buf_state->smallbuf[avail-MAX_BUFS].buffer,operation,to,
	     (int)_armci_buf_state->table[avail].count,(int)_armci_buf_state->table[avail].first);
      fflush(stdout);
    }

# ifdef BUF_EXTRA_FIELD_T
    INIT_SEND_BUF(_armci_buf_state->smallbuf[avail-MAX_BUFS].field,_armci_buf_state->table[avail].snd,_armci_buf_state->table[avail].rcv);
#endif
     
    _armci_buf_state->smavail = (avail+1-MAX_BUFS)%MAX_SMALL_BUFS + MAX_BUFS;

    if(DEBUG_ || 0) {
      printf("%d:buf_get_sm:size=%d max=%d got %d ptr=%p op=%d to=%d count=%d first=%d\n",
             armci_me,size,SMALL_BUF_LEN,avail,
	     _armci_buf_state->smallbuf[avail-MAX_BUFS].buffer,operation,to,
	     _armci_buf_state->table[avail].count,_armci_buf_state->table[avail].first);
      fflush(stdout);
    }
    
    return(_armci_buf_state->smallbuf[avail-MAX_BUFS].buffer); 

}

/*\  call corresponding to GET_SEND_BUF
\*/
static char *rmo_buffer = NULL;

char *_armci_buf_get(int size, int operation, int to)
{
#ifndef PORTALS_USE_ARMCI_CLIENT_BUFFERS
       if(rmo_buffer) return rmo_buffer;
       rmo_buffer = (char *) valloc(MSG_BUFLEN);
       return rmo_buffer;
#else
int avail=_armci_buf_state->avail;
int count=1, i;
_buf_ackresp_t *ar;

    /*if small buffer, we go to another routine that gets smallbuf*/
    if(size<SMALL_BUF_LEN) return(_armci_buf_get_small(size,operation,to));
    /* compute number of buffers needed (count) to satisfy the request */
    if((size > MSG_BUFLEN_SMALL) ){ 
       double val = (double)size;  /* use double due to a bug in gcc */
       val /= MSG_BUFLEN_SMALL;
       count=(int)val;
       if(size%MSG_BUFLEN_SMALL) count++; 
       assert(0);
    }
    /* start from 0 if there is not enough bufs available from here */
    if((avail+count) > MAX_BUFS)avail = 0;

    /* avail should never point to buffer in a middle of a set of used bufs */
    if(_armci_buf_state->table[avail].op && 
      (_armci_buf_state->table[avail].first != (unsigned int) avail)){ sleep(1); 
      printf("%d: inconsistent first. avail=%d count=%d first=%d size=%d\n",
	     armci_me, avail, count, _armci_buf_state->table[avail].first, size);
      armci_die2("armci_buf_get: inconsistent first", avail,
		 _armci_buf_state->table[avail].first);
    }
    
    /* we need complete "count" number of buffers */
    for(i=0;i<count;i++){
        int cur = i +avail;
        if((_armci_buf_state->table[cur].op &&
           _armci_buf_state->table[cur].first==(unsigned int) cur) ||
           _armci_buf_state->buf[cur].id.ar.req.active) {
              _armci_buf_complete_index(cur,1);
        }
    }

    for(i=0; i<count; i++){
       _armci_buf_state->table[avail+i].op = operation;
       _armci_buf_state->table[avail+i].to = to;
       _armci_buf_state->table[avail+i].count=  count;
       _armci_buf_state->table[avail+i].first = avail;
    }

    _armci_buf_state->buf[avail].id.tag=0;
    _armci_buf_state->buf[avail].id.bufid=avail; 
    _armci_buf_state->buf[avail].id.protocol=0;
    ar=&_armci_buf_state->buf[avail].id.ar;

    assert(ar->val==0);assert(ar->next==NULL);assert(ar->previous==NULL);
    assert(ar->req.active == 0);

    ar->req.active = 1;
    
    if(_buf_ackresp_cur!=NULL)
      _buf_ackresp_cur->next=ar;
    if(_buf_ackresp_first==NULL)
      _buf_ackresp_first=ar;
    ar->previous=_buf_ackresp_cur;
    ar->next=NULL;
    _buf_ackresp_cur = ar;

# ifdef BUF_EXTRA_FIELD_T
    INIT_SEND_BUF(_armci_buf_state->buf[avail].field,_armci_buf_state->table[avail].snd,_armci_buf_state->table[avail].rcv);
#endif

#ifdef HISTORY
    history[h].size=size;
    history[h].op=operation;
    history[h].count=count;
    history[h].id = avail;
    h++;
#endif

    if(DEBUG_ || 0) {
      printf("%d:buf_get:size=%d max=%d got %d ptr=%p count=%d op=%d to=%d\n",
             armci_me,size,MSG_BUFLEN_SMALL,avail,
            _armci_buf_state->buf[avail].buffer, count,operation,to);
      fflush(stdout);
    }

    /* select candidate buffer for next allocation request */
    _armci_buf_state->avail = avail+count;
    _armci_buf_state->avail %= MAX_BUFS;

    return(_armci_buf_state->buf[avail].buffer); 
#endif
}


void _armci_buf_release_index(int index) {
  int count;
  buf_state_t *buf_state = _armci_buf_state->table +index;
  char *_armci_buf_ptr_from_id(int id);  

  if((index >= MAX_BUFS+MAX_SMALL_BUFS)|| (index<0))
    armci_die2("armci_buf_release: bad index:",index,MAX_BUFS);
  
  count =  _armci_buf_state->table[index].count;
  
  if(DEBUG_ || 0) {
    printf("%d:_armci_buf_release_index %d ptr=%p count=%d op=%d smavail=%d\n",
	   armci_me,index,_armci_buf_ptr_from_id(index),count, _armci_buf_state->table[index].op,_armci_buf_state->smavail);
    fflush(stdout);
  }
  
  /* clear table slots for all the buffers in the set for this request */
  for(; count; count--, buf_state++) *(int*)buf_state = 0;
  if(index >= MAX_BUFS){
    _armci_buf_state->smallbuf[index-MAX_BUFS].id.tag=0;
    //_armci_buf_state->smavail = index;
  }
  else{
    _armci_buf_state->buf[index].id.tag=0;
   // _armci_buf_state->avail = index;
  }
  /* the current buffer is prime candidate to satisfy next buffer request */
}


/*\ release buffer when it becomes free
\*/
void _armci_buf_release(void *buf) {
#ifdef PORTALS_USE_ARMCI_CLIENT_BUFFERS
  _armci_buf_release_index(_armci_buf_to_index(buf));
#endif
}


/*\ return pointer to buffer number id
\*/
char *_armci_buf_ptr_from_id(int id)
{
  if(id <0 || id >=(MAX_BUFS+MAX_SMALL_BUFS)) 
              armci_die2("armci_buf_ptr_from_id: bad id",id,MAX_BUFS);
  if(id >=MAX_BUFS)return(_armci_buf_state->smallbuf[id-MAX_BUFS].buffer);
  return(_armci_buf_state->buf[id].buffer);
}



/*\function called from PARMCI_Wait to wait for non-blocking ops
\*/
void _armci_buf_complete_nb_request(int bufid,unsigned int tag, int *retcode) 
{
int i=0;
#if 0
    printf("\n%d:wait called with bufid=%d tag=%d \n",armci_me,bufid,tag);
    fflush(stdout);
#endif
 
    if(bufid == NB_NONE) *retcode=0;
    else if(bufid == NB_MULTI) {
       for(i=0;i<MAX_BUFS;i++){ 
         if(tag && tag==_armci_buf_state->buf[i].id.tag)
           _armci_buf_complete_index(i,1); 
       }
       for(i=0;i<MAX_SMALL_BUFS;i++){ 
         if(tag && tag==_armci_buf_state->smallbuf[i].id.tag)
           _armci_buf_complete_index(i+MAX_BUFS,1); 
       }
       *retcode=0;
    }
    else {
       if(bufid<MAX_BUFS){
         if(tag && tag==_armci_buf_state->buf[bufid].id.tag)
           _armci_buf_complete_index(bufid,1);
       }
       else{
         if(tag && tag==_armci_buf_state->smallbuf[bufid-MAX_BUFS].id.tag)
           _armci_buf_complete_index(bufid,1);
       }
       *retcode=0;
    } 
}


/*\function called from PARMCI_Test to test completion of non-blocking ops
\*/
void _armci_buf_test_nb_request(int bufid,unsigned int tag, int *retcode) 
{
int i;
    if(bufid == NB_NONE) *retcode=0;
    else if(bufid == NB_MULTI) {
       for(i=0;i<MAX_BUFS;i++){ 
         if(tag && tag==_armci_buf_state->buf[i].id.tag){
           if(_armci_buf_test_index(i,1)){
             *retcode=1;
	     break;
	   }
	 }
       }
       for(i=0;i<MAX_SMALL_BUFS;i++){ 
         if(tag && tag==_armci_buf_state->smallbuf[i].id.tag)
           if(_armci_buf_test_index(i+MAX_BUFS,1)){
             *retcode=1;
	     break;
	   }
       }
    }
    else {
       if(bufid<MAX_BUFS){
         if(tag && tag==_armci_buf_state->buf[bufid].id.tag)
           *retcode = _armci_buf_test_index(bufid,1);
       }
       else{
         if(tag && tag==_armci_buf_state->smallbuf[bufid-MAX_BUFS].id.tag)
           *retcode = _armci_buf_test_index(bufid,1);
       }
    }
}

/*\function to set the buffer tag and the protocol
\*/
void _armci_buf_set_tag(void *bufptr,unsigned int tag,short int protocol)
{
int  index = _armci_buf_to_index(bufptr);
   /*_armci_buf_state->table[index].async=1;*/
   if(index<MAX_BUFS){
      _armci_buf_state->buf[index].id.tag=tag;
      _armci_buf_state->buf[index].id.protocol=protocol;
   }
   else{
      _armci_buf_state->smallbuf[index-MAX_BUFS].id.tag=tag;
      _armci_buf_state->smallbuf[index-MAX_BUFS].id.protocol=protocol;
   }
}

int _armci_buf_get_tag(void *bufptr)
{
int  index = _armci_buf_to_index(bufptr);
    if(index<MAX_BUFS)
      return(_armci_buf_state->buf[index].id.tag);
    else
      return(_armci_buf_state->smallbuf[index-MAX_BUFS].id.tag);
}

/*\function to return bufinfo, given buf ptr
\*/
BUF_INFO_T *_armci_buf_to_bufinfo(void *buf){
    if(buf > (void *)_armci_buffers && buf < (void *)(_armci_buffers+MAX_BUFS)){
       return(&((BUF_TO_EBUF(buf))->id));
    }
   else if(buf > (void *)_armci_smbuffers && buf < (void *)(_armci_smbuffers+MAX_SMALL_BUFS)){
       return(&((BUF_TO_SMEBUF(buf))->id));
   }
   else {
        armci_die("armci_buf_to_index: bad pointer",0);
        return(0);
   }
}

/*\function to clear all buffers
\*/
void _armci_buf_clear_all()
{
int i; 
    for(i=0;i<MAX_BUFS;i++){
# ifdef BUF_EXTRA_FIELD_T
       if(_armci_buf_state->table[i].op || _armci_buf_state->table[i].first)
         CLEAR_SEND_BUF_FIELD(_armci_buf_state->buf[i].field,_armci_buf_state->table[i].snd,_armci_buf_state->table[i].rcv,_armci_buf_state->table[i].to,_armci_buf_state->table[i].op);
#endif
    }
    for(i=MAX_BUFS;i<MAX_BUFS+MAX_SMALL_BUFS;i++){
# ifdef BUF_EXTRA_FIELD_T
       if(_armci_buf_state->table[i].op || _armci_buf_state->table[i].first)
         CLEAR_SEND_BUF_FIELD(_armci_buf_state->smallbuf[i-MAX_BUFS].field,_armci_buf_state->table[i].snd,_armci_buf_state->table[i].rcv,_armci_buf_state->table[i].to,_armci_buf_state->table[i].op);
#endif
    }
}

/* function to return bufinfo, given buf tag */
BUF_INFO_T *_armci_tag_to_bufinfo(msg_tag_t tag) {
    int idx;

    for (idx=0; idx < MAX_BUFS; idx++)
        if (EQ_TAGS(_armci_buffers[idx].id.tag, tag)) break;

    if (idx == MAX_BUFS) {/* not found is regular buffers */
        for (idx = 0; idx < MAX_SMALL_BUFS; idx++)
            if (EQ_TAGS(_armci_smbuffers[idx].id.tag, tag)) break;
        if (idx == MAX_SMALL_BUFS) /* not found at all */
            armci_die("_armci_tag_to_bufinfo: bad tag",0);

        return &(_armci_smbuffers[idx].id);
    } else return &(_armci_buffers[idx].id);
}


/* inline primitives for buffer state management */
INLINE char *_armci_buf_get_clear_busy(int size, int operation, int to) {
    char *buf = _armci_buf_get(size, operation, to);
    _armci_buf_set_busy(buf, 0);
    return buf;
}

INLINE void _armci_buf_set_busy(void *buf, int state) {
        _armci_buf_state->table[_armci_buf_to_index(buf)].busy = state;
}

INLINE void _armci_buf_set_busy_idx(int idx, int state) {
    _armci_buf_state->table[idx].busy = state;
}

#if 0
INLINE int _armci_buf_cmpld(void *buf) {
    return _armci_buf_state->table[_armci_buf_to_index(buf)].cmpl;
}
#else
INLINE int _armci_buf_cmpld(int bufid) {
        return _armci_buf_state->table[bufid].cmpl;
}
#endif


INLINE void _armci_buf_set_cmpld(void *buf, int state) {
        _armci_buf_state->table[_armci_buf_to_index(buf)].cmpl = state;
}

INLINE void _armci_buf_set_cmpld_idx(int idx, int state) {
    _armci_buf_state->table[idx].cmpl = state;
}


