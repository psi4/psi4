#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: buffers.c,v 1.29.6.9 2007-07-02 05:16:50 d3p687 Exp $    **/

#define SIXTYFOUR 64
#define DEBUG_  0
#define DEBUG2_ 0
#define EXTRA_ERR_CHECK

/**********************************************************************/
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif
#include "armcip.h"
#include "request.h"

#if HAVE_UNISTD_H
#   include <unistd.h>
#elif HAVE_WINDOWS_H
#   include <windows.h>
    typedef unsigned long ssize_t;
#endif

#ifdef SOCKETS
#   define EQ_TAGS(a_, b_) ((a_) == (b_))
#else
#   define EQ_TAGS(a_, b_) !memcmp(&(a_), &(b_), sizeof(a_))
#endif

#define ALIGN64ADD(buf) (SIXTYFOUR-(((ssize_t)(buf))%SIXTYFOUR))
/* the following symbols should be defined if needed in protocol specific
   header file:  BUF_EXTRA_FIELD, BUF_ALLOCATE
*/

#ifndef BUF_ALLOCATE
#   define BUF_ALLOCATE malloc
#endif
#if defined(PEND_BUFS)
#define SMALL_BUF_LEN 8192
#else
#define SMALL_BUF_LEN 8192
#endif
#ifndef MSG_BUFLEN_SMALL
#define MSG_BUFLEN_SMALL (MSG_BUFLEN >>0)
#endif
#define LEFT_GUARD  11.11e11
#define RIGHT_GUARD 22.22e22
/* #define CLEAR_TABLE_SLOT(idx) *((int*)(_armci_buf_state->table+(idx))) =0 */
#define CLEAR_TABLE_SLOT(idx) (memset(_armci_buf_state->table+(idx),'\0',sizeof(buf_state_t))


/* we allow multiple buffers (up to 15) per single request
 * adjacent buffers can be coalesced into a large one
 */
#if 0
typedef struct {
  unsigned int op:8;     /* pending operation code */
  unsigned int snd:1;    /* if 1 then buffer is used for sending request */
  unsigned int rcv:1;    /* if 1 then buffer is used for receiving data */
  unsigned int async:1;  /* if 1 then request is nonblocking */
  unsigned int first:5;  /* id of the 1st buffer in the set in same request */
  unsigned int count:1;  /* count is not used and is always 1 (or 0???) */
  /*unsigned int count:4;  \* how many buffers used for this request 8 possible */
  unsigned int busy:1;   /* if 1 buffer is used and cannot be completed */
  unsigned int cmpl:1;   /* set to 1 if buffer was completed and can be released */
  unsigned int to:13;    /* serv/proc to which request was sent, 8172 possible */
}buf_state_t;
#else
typedef struct {
  unsigned int op:8;     /* pending operation code */
  unsigned int snd:1;    /* if 1 then buffer is used for sending request */
  unsigned int rcv:1;    /* if 1 then buffer is used for receiving data */
  unsigned int async:1;  /* if 1 then request is nonblocking */
  unsigned int first:20;  /* id of the 1st buffer in the set in same request */
  unsigned int count:1;  /* count is not used and is always 1 (or 0???) */
  /*unsigned int count:4;  \* how many buffers used for this request 8 possible */
  unsigned int busy:1;   /* if 1 buffer is used and cannot be completed */
  unsigned int cmpl:1;   /* set to 1 if buffer was completed and can be released */
  unsigned int to:30;    /* serv/proc to which request was sent, can handle pretty large counts (can be used for fields in the future) */
}buf_state_t;
#endif

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


#if (defined(THREAD_SAFE) || defined(SOCKETS))

/* check if buffer was completed and can be released */
int armci_test_network_complete() {
    int idx;
    for (idx=0;idx<MAX_BUFS+MAX_SMALL_BUFS;idx++)
        if (_armci_buf_state->table[idx].cmpl) break;
    return (idx == MAX_BUFS+MAX_SMALL_BUFS ? -1 : idx);
}


/*\ we allocate alligned buffer space
 *  this operation can be implemented in platform specific files
\*/
void _armci_buf_init()
{
char *tmp;
int extra=0;
int smallbuf_size = sizeof(buf_smext_t)*(MAX_SMALL_BUFS);
     tmp = (char *)BUF_ALLOCATE((MAX_BUFS*sizeof(buf_ext_t) + 64 + smallbuf_size));
     extra= ALIGN64ADD(tmp);
/*      if(sizeof(buf_state_t) != sizeof(int)) */
/*         armci_die("armci_buf_init size buf_state_t!=int",sizeof(buf_state_t)); */
     dassert(1,MAX_BUFS<sizeof(int)*8); /*should fit in the bitmap*/
     dassert(1,MAX_SMALL_BUFS<sizeof(int)*8); /*should fit in the bitmap*/

     _armci_buffers = (buf_ext_t *) (tmp + extra);

     tmp = (char *)(_armci_buffers + MAX_BUFS);
     extra = ALIGN64ADD(tmp);
     _armci_smbuffers = (buf_smext_t *) (tmp + extra);

     if(DEBUG2_){
	printf("%d:armci_init_bufs: pointer %p, before align ptr=%p bufptr=%p end of region is %p  size=%lu extra=%d\n",
               armci_me, (void*)_armci_buffers, tmp, _armci_buffers->buffer,
               (void*)(_armci_buffers+MAX_BUFS),
               (long unsigned)MAX_BUFS*sizeof(buf_ext_t), extra);
	fflush(stdout);
     }

     /* now allocate state array */
     tmp  = calloc(1, sizeof(reqbuf_pool_t) + 64);
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

     /* initialize bitmaps */
     /* should not be more than sizeof(unsigned int) buffers */
     if (MAX_BUFS > sizeof(unsigned)*8 && MAX_SMALL_BUFS > sizeof(unsigned)*8)
         armci_die("_armci_buf_init: cannot allocate this many buffers",0);

     _armci_buf_state->buf_bitmap = (1 << MAX_BUFS) - 1;
     _armci_buf_state->smbuf_bitmap = (1 << MAX_SMALL_BUFS) - 1;

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
        armci_die2("armci_buf_to_index: index wrong",index,MAX_BUFS+MAX_SMALL_BUFS);
      return(index);
   } 
   else {
        armci_die("armci_buf_to_index: bad pointer",0);
        return(0);
   }
}



/*\  complete outstanding operation that uses the specified buffer
\*/
void _armci_buf_complete_index(int idx, int called)
{
int count;
buf_state_t *buf_state = _armci_buf_state->table +idx;

/*  fprintf(stderr, "%d:: entered %s. called=%d\n", armci_me, FUNCTION_NAME); */

    count = buf_state->count;
    if(DEBUG_ ) {
       printf("%d:buf_complete_index:%d op=%d first=%d count=%d called=%d\n",
              armci_me,idx,buf_state->op,buf_state->first,buf_state->count,
              called); 
       fflush(stdout);
    }

    if(buf_state->first != (unsigned int)idx){ 
      armci_die2("complete_buf_index:Inconsistent index:",idx,buf_state->first);
    }

    if(buf_state->async){
      /* completion of strided get should release that buffer */
      if(buf_state->op == GET);
      else
         armci_die2("buf_complete_index: async mode not avail for this op",
                     buf_state->op,idx);
    }
#   ifdef BUF_EXTRA_FIELD_T
    else{
       /* need to call platform specific function */
       if(idx>=MAX_BUFS){
         int relidx;
         relidx = idx-MAX_BUFS; 
         
         CLEAR_SEND_BUF_FIELD(_armci_buf_state->smallbuf[relidx].field,buf_state->snd,buf_state->rcv,buf_state->to,buf_state->op);

       /*later, we might just need to do this for all operations, not just get*/
         if(_armci_buf_state->smallbuf[relidx].id.tag!=0 &&(buf_state->op == GET)){
          armci_complete_req_buf(&(_armci_buf_state->smallbuf[relidx].id),
                                _armci_buf_state->smallbuf[relidx].buffer);
         }
         _armci_buf_state->smallbuf[relidx].id.tag=0;
       }
       else {
         CLEAR_SEND_BUF_FIELD(_armci_buf_state->buf[idx].field,buf_state->snd,buf_state->rcv,buf_state->to,buf_state->op);

       /*later, we might just need to do this for all operations, not just get*/
         if(_armci_buf_state->buf[idx].id.tag!=0 &&(buf_state->op == GET)){
           armci_complete_req_buf(&(_armci_buf_state->buf[idx].id),
                                _armci_buf_state->buf[idx].buffer);
         }
         _armci_buf_state->buf[idx].id.tag=0;
       }
    }
#   endif

    /* clear table slots for all the buffers in the set for this request */
    for(; count; count--, buf_state++) {
/*       *(int*)buf_state = 0; */
      memset(buf_state,'\0',sizeof(buf_state_t));
    }
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
        if((buf_state->to >= nfirst) && (buf_state->to<= (unsigned int) nlast))
          if((buf_state->first != (unsigned int) this)&&(buf_state->first==(unsigned int) i) && buf_state->op){
                _armci_buf_complete_index(i,0);
                _armci_buf_release_index(i);
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

    if(_armci_buf_state->table[index].to !=(unsigned int)  proc )
       armci_die2("_armci_buf_ensure_one_outstanding_op_per_proc: bad to", proc,
                (int)_armci_buf_state->table[index].to);

    for(i=0;i<MAX_BUFS+MAX_SMALL_BUFS;i++){
        buf_state_t *buf_state = _armci_buf_state->table +i;
        if(buf_state->to == (unsigned int) proc)
          if((buf_state->first != (unsigned int) this)&&(buf_state->first==(unsigned int) i) && buf_state->op)
                _armci_buf_complete_index(i,0);
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

enum {CALL_GET, CALL_RELEASE} last_call_ = CALL_RELEASE;
#if DEBUG3_
static int get_count_ = 0;
static int release_count_ = 0;
#endif


/* release buffer, update free buffers bitmaps */
void _armci_buf_release_index(int tbl_idx) {
    int idx;

#if DEBUG3_
     release_count_++;
     if (last_call_ != CALL_GET) {
         printf("same call: trying to release %s buffer\n"
                "           %d release calls to date\n",
                tbl_idx >= MAX_BUFS ? "small" : "regular", release_count_);
     } else last_call_ = CALL_RELEASE;

     int mark = !(tbl_idx < MAX_BUFS ? _armci_buf_state->buf_bitmap
                                     : _armci_buf_state->smbuf_bitmap);
#endif

    if ((tbl_idx >= MAX_BUFS+MAX_SMALL_BUFS) || (tbl_idx < 0))
      armci_die2("armci_buf_release: bad index:", tbl_idx, MAX_BUFS);

    THREAD_LOCK(armci_user_threads.buf_lock);

    if (tbl_idx < MAX_BUFS) {
        idx = tbl_idx;
        _armci_buf_state->buf_bitmap |= 1 << idx;
        _armci_buf_state->buf[idx].id.tag = 0;
    } else {
        idx = tbl_idx - MAX_BUFS;
        _armci_buf_state->smbuf_bitmap |= 1 << idx;
        _armci_buf_state->smallbuf[idx].id.tag = 0;
    }
    _armci_buf_state->table[tbl_idx].busy = 0;

    THREAD_UNLOCK(armci_user_threads.buf_lock);

#if DEBUG3_
    if (mark) printf("release of empty buffer pool, after: %d\n",
                     tbl_idx < MAX_BUFS ? _armci_buf_state->buf_bitmap
                                        : _armci_buf_state->smbuf_bitmap);
#endif
}

/*\ release buffer when it becomes free
\*/
INLINE void _armci_buf_release(void *buf)
{
    _armci_buf_release_index(_armci_buf_to_index(buf));
}

/*  call corresponding to GET_SEND_BUF */
char *_armci_buf_get(int size, int operation, int to)
{
    int avail; /* needed by INIT_SEND_BUF for MELLANOX */
    unsigned bitmap;
    int small = size < SMALL_BUF_LEN;
    int max_bufs = small ? MAX_SMALL_BUFS : MAX_BUFS;
    int idx = 0; /* same type buffers index: 0..{MAX_BUFS|MAX_SMALL_BUFS}-1 */
    int tbl_idx; /* global index in table: 0..MAX_BUFS+MAX_SMALL_BUFS-1 */
    int not_ready = 1;

#if DEBUG3_
     get_count_++;
     if (last_call_ != CALL_RELEASE) {
         printf("same call: trying to get %s buffer\n"
                "           %d get calls to date, %d release calls to date, bitmap: %d\n",
                small ? "small" : "regular", get_count_, release_count_,
                small ? _armci_buf_state->smbuf_bitmap : _armci_buf_state->buf_bitmap);
     } else last_call_ = CALL_GET;
     /* debug hook: no buffers left */
     if (!(small ? _armci_buf_state->smbuf_bitmap : _armci_buf_state->buf_bitmap))
         bitmap = small ? _armci_buf_state->smbuf_bitmap : _armci_buf_state->buf_bitmap;
#endif

    while (not_ready) {
        THREAD_LOCK(armci_user_threads.buf_lock);

        bitmap = small ? _armci_buf_state->smbuf_bitmap
                       : _armci_buf_state->buf_bitmap;
        /* check if there are available buffers */
        if (bitmap) {
            /* find available buffer in the bitmap */
            for (idx = 0; idx < max_bufs; idx++) {
                if (bitmap & (1 << idx)) break;
            }
            if (idx >= max_bufs)
                armci_die("_armci_buf_get: buffer idx is out of the range",idx);

            /* mark buffer as taken in the bitmap and busy */
            bitmap &= ~((unsigned)(1 << idx));
            if (small) {
                tbl_idx = idx + MAX_BUFS;
                _armci_buf_state->smbuf_bitmap = bitmap;
            } else {
                tbl_idx = idx;
                _armci_buf_state->buf_bitmap = bitmap;
            }
#if 0
            _armci_buf_state->table[tbl_idx].busy = 1;
#endif

            THREAD_UNLOCK(armci_user_threads.buf_lock);
            not_ready = 0;
        } else {
            THREAD_UNLOCK(armci_user_threads.buf_lock);

            /* try network complete */
#if defined(SOCKETS) || defined(MELLANOX)
            tbl_idx = armci_test_network_complete();
#else /* all network should eventually use armci_test_network_complete */
	    tbl_idx = small ? _armci_buf_state->smavail : _armci_buf_state->avail;
#endif
            avail = tbl_idx;
            if ((tbl_idx >= MAX_BUFS+MAX_SMALL_BUFS) || (tbl_idx < 0 && tbl_idx != -1))
                armci_die2("_armci_buf_get: bad idx:", tbl_idx, MAX_BUFS);

            if (tbl_idx < 0) {
                /*printf("armci_test_network_complete returned -1\n");fflush(stdout);*/
                cpu_yield();
            } else { /* could complete a buffer */
                /* ignore a busy buffer */
#if DEBUG3_
                printf("armci_test_network_complete, gets:%d, releases:%d\n",
                       get_count_, release_count_);
#endif
                if (_armci_buf_state->table[tbl_idx].busy) {
                    printf("BUFFER BUSY 1\n");
                    continue;
                }

                /* complete buffer */
                _armci_buf_complete_index(tbl_idx, 0);

                /* tbl_idx < MAX_BUFS ^ small - 1 if completed compatible buffer */
                if ((tbl_idx < MAX_BUFS) ^ small) {
                    THREAD_LOCK(armci_user_threads.buf_lock);

                    /* is this check really necessary ??? */
                    if (!_armci_buf_state->table[tbl_idx].busy) {
#if 0
                        _armci_buf_state->table[tbl_idx].busy = 1;
#endif
                        not_ready = 0;
#if DEBUG3_
                        release_count_++;
                        printf("released in get, gets:%d, releases:%d\n",
                                get_count_, release_count_);
#endif
                    } else {
                        printf("BUFFER BUSY 2\n");
                    }

                    THREAD_UNLOCK(armci_user_threads.buf_lock);
                    idx = small ? tbl_idx - MAX_BUFS : tbl_idx;
                } else
                    _armci_buf_release_index(tbl_idx);
            }
        }
    }

    /* initialize buffer */
    _armci_buf_state->table[tbl_idx].op = operation;
    _armci_buf_state->table[tbl_idx].to = to;
    _armci_buf_state->table[tbl_idx].count = 1;
    _armci_buf_state->table[tbl_idx].first = tbl_idx;
    _armci_buf_state->table[tbl_idx].cmpl = 0;

    /* Note: tbl_idx is used in vapi vesrion of INIT_SEND_BUF */
    if (small) {
        _armci_buf_state->smallbuf[idx].id.tag = 0;
        _armci_buf_state->smallbuf[idx].id.bufid = tbl_idx;
        _armci_buf_state->smallbuf[idx].id.protocol = 0;
# ifdef BUF_EXTRA_FIELD_T
         INIT_SEND_BUF(_armci_buf_state->smallbuf[idx].field,
                      _armci_buf_state->table[tbl_idx].snd,
                      _armci_buf_state->table[tbl_idx].rcv);
#endif
    } else {
        _armci_buf_state->buf[idx].id.tag = 0;
        _armci_buf_state->buf[idx].id.bufid = tbl_idx;
        _armci_buf_state->buf[idx].id.protocol = 0;
# ifdef BUF_EXTRA_FIELD_T
        INIT_SEND_BUF(_armci_buf_state->buf[idx].field,
                      _armci_buf_state->table[idx].snd,
                      _armci_buf_state->table[idx].rcv);
#endif
    }
    return small ? _armci_buf_state->smallbuf[idx].buffer
                 : _armci_buf_state->buf[idx].buffer;
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

#ifdef VAPI
/* this will work for vapi as there is no get pipeline enabled in vapi
 * with get pipeline, this will break very badly
 */
void _armci_buf_update_scatter_count(int id)
{
int i,num,last,first;
    for(i=0;i<MAX_BUFS;i++){
# ifdef BUF_EXTRA_FIELD_T
       if(_armci_buf_state->table[i].op==GET){
         request_header_t *msginfo;
	 msginfo = (request_header_t*)_armci_buf_state->buf[i].buffer; 
         if(msginfo->pinned && msginfo->bypass && msginfo->format == STRIDED){
           num = *(int *)((char *)msginfo+msginfo->bytes); 
           last = *(int *)((char *)msginfo+msginfo->bytes+sizeof(int));
           first = last - num+1;
           if(first < 0 )first+=DSCRID_SCATTERCLIENT_END-DSCRID_SCATTERCLIENT-1;
           if(id == first && num!=0){
             *(int *)((char *)msginfo+msginfo->bytes) = (--num);
             return;
           }
         }
       }
# endif
    }
    for(i=MAX_BUFS;i<MAX_BUFS+MAX_SMALL_BUFS;i++){
# ifdef BUF_EXTRA_FIELD_T
       if(_armci_buf_state->table[i].op==GET){
         request_header_t *msginfo = (request_header_t*)_armci_buf_state->smallbuf[i-MAX_BUFS].buffer;
         if(msginfo->pinned && msginfo->bypass && msginfo->format == STRIDED){
           num = *(int *)((char *)msginfo+msginfo->bytes); 
           last = *(int *)((char *)msginfo+msginfo->bytes+sizeof(int));
           first = last - num+1;
           if(first < 0 )first+=DSCRID_SCATTERCLIENT_END-DSCRID_SCATTERCLIENT-1;
           if(id == first && num!=0){
             *(int *)((char *)msginfo+msginfo->bytes) = (--num);
             return;
           }
         }
       }
# endif
    }

}
#endif


#else /* (defined(THREAD_SAFE) || defined(SOCKETS)) */


/*\ we allocate alligned buffer space
 *  this operation can be implemented in platform specific files
\*/ 
void _armci_buf_init()
{
char *tmp;
int  extra=0;
int smallbuf_size = sizeof(buf_smext_t)*(MAX_SMALL_BUFS);
     tmp = (char *)BUF_ALLOCATE((MAX_BUFS*sizeof(buf_ext_t) + 64 + smallbuf_size));
     extra= ALIGN64ADD(tmp);
#if 0
     if(sizeof(buf_state_t) != sizeof(int)) 
        armci_die("armci_buf_init size buf_state_t!=int",sizeof(buf_state_t));
#endif                   
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
     tmp  = calloc(1, sizeof(reqbuf_pool_t) + 64);
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



/*\  complete outstanding operation that uses the specified buffer
\*/
void _armci_buf_complete_index(int idx, int called)
{
int count;
buf_state_t *buf_state = _armci_buf_state->table +idx;

/* printf("%d: buf_complete_index. idx=%d\n", armci_me, idx);*/
/* fflush(stdout);*/

/*  if(buf_state->op==GET) { */
/*    printf("%d: %s(): op is get\n",armci_me,FUNCTION_NAME); */
/*  } */

    count = buf_state->count;
    if(DEBUG_ ) {
       printf("%d:buf_complete_index:%d op=%d first=%d count=%d called=%d\n",
              armci_me,idx,buf_state->op,buf_state->first,buf_state->count,
              called); 
       fflush(stdout);
    }

    if(buf_state->first != (unsigned int)idx){ 
      armci_die2("complete_buf_index:inconsistent Index:",idx,buf_state->first);
    }

    if(buf_state->async){
      /* completion of strided get should release that buffer */
      if(buf_state->op == GET);
      else
         armci_die2("buf_complete_index: async mode not avail for this op",
                     buf_state->op,idx);
    }
#   ifdef BUF_EXTRA_FIELD_T
    else{
       /* need to call platform specific function */
       if(idx>=MAX_BUFS){
         int relidx;
         relidx = idx-MAX_BUFS; 

/* 	 printf("%d:%s(): Calling clear_send_buf_field\n",armci_me,FUNCTION_NAME); */
         CLEAR_SEND_BUF_FIELD(_armci_buf_state->smallbuf[relidx].field,buf_state->snd,buf_state->rcv,buf_state->to,buf_state->op);

       /*later, we might just need to do this for all operations, not just get*/
         if(_armci_buf_state->smallbuf[relidx].id.tag!=0 &&(buf_state->op == GET)){
          armci_complete_req_buf(&(_armci_buf_state->smallbuf[relidx].id),
                                _armci_buf_state->smallbuf[relidx].buffer);
         }
         _armci_buf_state->smallbuf[relidx].id.tag=0;
       }
       else {
         CLEAR_SEND_BUF_FIELD(_armci_buf_state->buf[idx].field,buf_state->snd,buf_state->rcv,buf_state->to,buf_state->op);

       /*later, we might just need to do this for all operations, not just get*/
         if(_armci_buf_state->buf[idx].id.tag!=0 &&(buf_state->op == GET)){
           armci_complete_req_buf(&(_armci_buf_state->buf[idx].id),
                                _armci_buf_state->buf[idx].buffer);
         }
         _armci_buf_state->buf[idx].id.tag=0;
       }
    }
#   endif

    /* clear table slots for all the buffers in the set for this request */
    assert(count==1);
    for(; count; count--, buf_state++) {
/*       *(int*)buf_state = 0; */
      memset(buf_state, '\0', sizeof(buf_state_t));
    }
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
#if defined(PEND_BUFS)
/** Need to implement effective mechanisms to ensure a given number of outstanding operations. The key is to not wait for a recently sent message. Another issue might be priority between the large and small messages. This function will be improved as needed.  X_BUFS+
 */
void _armci_buf_ensure_pend_outstanding_op_per_node(void *buf, int node)
{
unsigned int i;
int index =_armci_buf_to_index(buf); 
int this = _armci_buf_state->table[index].first;
unsigned int nfirst, nlast;
void _armci_buf_release_index(int i);
int buf_pend_count=0;
 int max_pend_count=IMM_BUF_NUM; 

/* printf("%d: ensure_pend_os_per_node. idx=%d\n", armci_me, index);*/
/* fflush(stdout);*/

    nfirst=armci_clus_info[node].master;
    nlast = nfirst+armci_clus_info[node].nslave-1;

    if((_armci_buf_state->table[index].to<(unsigned int) nfirst) || 
      (_armci_buf_state->table[index].to>(unsigned int) nlast))
        armci_die2("_armci_buf_ensure_pend_outstanding_op_per_node: bad to",node,
                        (int)_armci_buf_state->table[index].to);

#if 0
    for(i=0;i<MAX_BUFS+MAX_SMALL_BUFS;i++){
      buf_state_t *buf_state = _armci_buf_state->table +i;
      if((buf_state->to >= nfirst) && (buf_state->to<= (unsigned int) nlast))
        if((buf_state->first != (unsigned int) this)&&(buf_state->first==(unsigned int) i) && buf_state->op){
          buf_pend_count++;
          if(buf_pend_count > IMM_BUF_NUM-1){
/* 	    printf("%d: client. pend_os_per_node. completing buf index=%d\n", armci_me, i); */
/* 	    fflush(stdout); */
            _armci_buf_complete_index(i,0);
	    _armci_buf_release_index(i);
/* 	    printf("%d: client. pend_os_per_node. done completing buf index=%d\n", armci_me, i); */
	    fflush(stdout);
            buf_pend_count--;
          }
        }
    }
#else
    const int smavail=_armci_buf_state->smavail;
    const int avail = _armci_buf_state->avail;
    const int startsmall = (smavail-MAX_BUFS-1+MAX_SMALL_BUFS)%MAX_SMALL_BUFS+MAX_BUFS;
    const int start =(avail-1+MAX_BUFS)%MAX_BUFS;
    buf_pend_count=0;
    buf_state_t *bs;

    i = start;
    bs = &_armci_buf_state->table[i];
    do {
      if((bs->to>=nfirst) && (bs->to<=nlast)) {
	if((bs->first!=this) && (bs->first==i) && bs->op) {
#if defined(OPENIB) /*SK: not tested on other platforms*/
	  if(!_armci_buf_test_index(i,1)) {
	    buf_pend_count++;
	  }
#else
	  buf_pend_count++;
#endif
	  if(buf_pend_count > max_pend_count-1) {
/* 	    printf("%d:%s():complete largebuf %d\n",armci_me,FUNCTION_NAME,i); */
	    _armci_buf_complete_index(i,0);
	    _armci_buf_release_index(i);
	    buf_pend_count--;
	  }
	}
      }
      i = (i-1+MAX_BUFS)%MAX_BUFS;
      bs = &_armci_buf_state->table[i];
    } while(i != start);

/*     printf("%d: largebuf pend=%d\n",armci_me,buf_pend_count); */

    i = startsmall;
    bs = &_armci_buf_state->table[i];
    do {
      if((bs->to>=nfirst) && (bs->to<=nlast)) {
	if((bs->first!=this) && (bs->first==i) && bs->op) {
#if defined(OPENIB) /*SK:not tested on other platforms*/
	  if(!_armci_buf_test_index(i,1)) {
	    buf_pend_count++;
	  }
#else
	  buf_pend_count++;
#endif
	  if(buf_pend_count > max_pend_count-1) {
/* 	    printf("%d:%s():complete smallbuf %d\n",armci_me,FUNCTION_NAME,i); */
	    _armci_buf_complete_index(i,0);
	    _armci_buf_release_index(i);
	    buf_pend_count--;
	  }
	}
      }
      i = (i-MAX_BUFS-1+MAX_SMALL_BUFS)%MAX_SMALL_BUFS+MAX_BUFS;
      bs = &_armci_buf_state->table[i];
    } while(i != startsmall);


    buf_pend_count=0;
    for(i=0; i<MAX_BUFS+MAX_SMALL_BUFS; i++) {
      if((bs->to>=nfirst) && (bs->to<=nlast)) {
	if((bs->first!=this) && (bs->first==i) && bs->op) {
	  buf_pend_count += 1;
	}
      }      
    }
    assert(buf_pend_count <= max_pend_count);
#endif
}
void _armci_buf_ensure_pend_outstanding_op_per_proc(void *buf, int proc)
{
int i;
int index = _armci_buf_to_index(buf); 
int this = _armci_buf_state->table[index].first;
void _armci_buf_release_index(int i);
int buf_pend_count=0;

    if(_armci_buf_state->table[index].to !=(unsigned int)  proc )
       armci_die2("_armci_buf_ensure_one_outstanding_op_per_proc: bad to", proc,
                (int)_armci_buf_state->table[index].to);

    for(i=0;i<MAX_BUFS+MAX_SMALL_BUFS;i++){
      buf_state_t *buf_state = _armci_buf_state->table +i;
      if(buf_state->to == (unsigned int) proc) {
	if((buf_state->first != (unsigned int) this)&&(buf_state->first==(unsigned int) i) && buf_state->op){
          buf_pend_count++;
          if(buf_pend_count > IMM_BUF_NUM-1){
	    _armci_buf_complete_index(i,0);
	    _armci_buf_release_index(i);
	    buf_pend_count--;
          }
	}
      }
    }
}
#endif

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
/*  int ous_in=0,ous_out; */
/*  for(i=MAX_BUFS; i<MAX_BUFS+MAX_SMALL_BUFS; i++) { */
/*    if(_armci_buf_state->table[i].first) { */
/*      assert(_armci_buf_state->table[i].op); */
/*    } */
/*    if(_armci_buf_state->table[i].op) { */
/*      assert(_armci_buf_state->table[i].first==i); */
/*      ous_in++; */
/*    } */
/*  } */
    if((_armci_buf_state->table[avail].op || 
       _armci_buf_state->table[avail].first)) {
       for(i=MAX_BUFS;i<MAX_BUFS+MAX_SMALL_BUFS;i++){
          if(!_armci_buf_state->table[i].op &&
             !_armci_buf_state->table[i].first)
            break;
#if defined(OPENIB) /*not tested on other platforms*/
	  if(_armci_buf_test_index(i,1)) {
/* 	    ous_in-=1; */
	    break;
	  }
#endif
       }
       if(i<(MAX_SMALL_BUFS+MAX_BUFS))avail = i;
       else {
/* 	 printf("%d: no free smallbuf.  complete index %d\n",armci_me,avail); */
          _armci_buf_complete_index(avail,1);
/* 	  ous_in-=1; */
       }
    }
    _armci_buf_state->table[avail].op = operation;
    _armci_buf_state->table[avail].to = to;
    _armci_buf_state->table[avail].count=  1;
    _armci_buf_state->table[avail].first = avail;
    _armci_buf_state->smallbuf[avail-MAX_BUFS].id.tag=0;
    _armci_buf_state->smallbuf[avail-MAX_BUFS].id.bufid= avail; 
    _armci_buf_state->smallbuf[avail-MAX_BUFS].id.protocol=0;
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

/*     ous_out=0; */
/*     for(i=MAX_BUFS; i<MAX_BUFS+MAX_SMALL_BUFS; i++) { */
/*       if(_armci_buf_state->table[i].first) { */
/* 	assert(_armci_buf_state->table[i].op); */
/*       } */
/*       if(_armci_buf_state->table[i].op) { */
/* 	assert(_armci_buf_state->table[i].first==i); */
/* 	ous_out++; */
/*       } */
/*     } */
/*     assert(ous_in+1 == ous_out); */

    return(_armci_buf_state->smallbuf[avail-MAX_BUFS].buffer); 
}

/*\  call corresponding to GET_SEND_BUF
\*/
char *_armci_buf_get(int size, int operation, int to)
{
int avail=_armci_buf_state->avail;
int count=1, i;
/*  int ous_in, ous_out; */
    /*if small buffer, we go to another routine that gets smallbuf*/
    if(size<SMALL_BUF_LEN) return(_armci_buf_get_small(size,operation,to));
    /* compute number of buffers needed (count) to satisfy the request */
    if((size > MSG_BUFLEN_SMALL) ){ 
       double val = (double)size;  /* use double due to a bug in gcc */
       val /= MSG_BUFLEN_SMALL;
       count=(int)val;
       if(size%MSG_BUFLEN_SMALL) count++; 
       assert(0); /*FOR NOW:Should never require multiple buffers; what else is pack.c for?*/
    }
    /* start from 0 if there is not enough bufs available from here */
    if((avail+count) > MAX_BUFS)avail = 0;

/*     ous_in=0; */
/*     for(i=0; i<MAX_BUFS; i++) { */
/*       if(_armci_buf_state->table[i].first) { */
/* 	assert(_armci_buf_state->table[i].op); */
/*       } */
/*       if(_armci_buf_state->table[i].op) { */
/* 	assert(_armci_buf_state->table[i].first==i); */
/* 	ous_in++; */
/*       } */
/*     } */


    /* avail should never point to buffer in a middle of a set of used bufs */
    if(_armci_buf_state->table[avail].op && 
      (_armci_buf_state->table[avail].first != (unsigned int) avail)){ sleep(1); 
      printf("%d: inconsistent first. avail=%d count=%d first=%d size=%d\n",
	     armci_me, avail, count, _armci_buf_state->table[avail].first, size);
      armci_die2("armci_buf_get: inconsistent first", avail,
		 _armci_buf_state->table[avail].first);
    }

#if 0
    /* we need complete "count" number of buffers */
    for(i=0;i<count;i++){
        int cur = i +avail;
        if(_armci_buf_state->table[cur].op &&
           _armci_buf_state->table[cur].first==(unsigned int) cur)
                                   _armci_buf_complete_index(cur,1);
    }
#else
    assert(count == 1);
    if(_armci_buf_state->table[avail].op ||
       _armci_buf_state->table[avail].first) {
      for(i=0; i<MAX_BUFS; i++) {
	if(!_armci_buf_state->table[i].op &&
	   !_armci_buf_state->table[i].first)
	  break;
#if defined(OPENIB)
	if(_armci_buf_test_index(i,1)) {
/* 	  ous_in -= 1; */
	  break;
	}
#endif
      }
      if(i<MAX_BUFS) avail=i;
      else {
/* 	printf("%d: no free large buf. completing index=%d\n",armci_me,avail); */
	_armci_buf_complete_index(avail,1);
/* 	ous_in -= 1; */
      }
    }    
#endif

    for(i=0; i<count; i++){
       _armci_buf_state->table[avail+i].op = operation;
       _armci_buf_state->table[avail+i].to = to;
       _armci_buf_state->table[avail+i].count=  count;
       _armci_buf_state->table[avail+i].first = avail;
    }

    _armci_buf_state->buf[avail].id.tag=0;
    _armci_buf_state->buf[avail].id.bufid=avail; 
    _armci_buf_state->buf[avail].id.protocol=0;

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

/*     ous_out=0; */
/*     for(i=0; i<MAX_BUFS; i++) { */
/*       if(_armci_buf_state->table[i].first) { */
/* 	assert(_armci_buf_state->table[i].op); */
/*       } */
/*       if(_armci_buf_state->table[i].op) { */
/* 	assert(_armci_buf_state->table[i].first==i); */
/* 	ous_out++; */
/*       } */
/*     } */
/*     assert(ous_in+1 == ous_out); */
    return(_armci_buf_state->buf[avail].buffer); 
}


void _armci_buf_release_index(int index) {
  int count;
  buf_state_t *buf_state = _armci_buf_state->table +index;
  char *_armci_buf_ptr_from_id(int id);  

  if((index >= MAX_BUFS+MAX_SMALL_BUFS)|| (index<0))
    armci_die2("armci_buf_release: bad index:",index,MAX_BUFS);
  
  count =  _armci_buf_state->table[index].count;
  
  if(DEBUG_) {
    printf("%d:_armci_buf_release_index %d ptr=%p count=%d op=%d smavail=%d\n",
	   armci_me,index,_armci_buf_ptr_from_id(index),count, _armci_buf_state->table[index].op,_armci_buf_state->smavail);
    fflush(stdout);
  }
  
  /* clear table slots for all the buffers in the set for this request */
  for(; count; count--, buf_state++) {
    memset(buf_state, '\0', sizeof(buf_state_t));
/*     *(int*)buf_state = 0; */
  }
  if(index >= MAX_BUFS){
    _armci_buf_state->smallbuf[index-MAX_BUFS].id.tag=0;
#ifndef VAPI
    _armci_buf_state->smavail = index;
#endif
  }
  else{
    _armci_buf_state->buf[index].id.tag=0;
#ifndef VAPI
    _armci_buf_state->avail = index;
#endif
  }
  /* the current buffer is prime candidate to satisfy next buffer request */
}


/*\ release buffer when it becomes free
\*/
void _armci_buf_release(void *buf) {
  _armci_buf_release_index(_armci_buf_to_index(buf));
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
#if 0
       for(i=0;i<MAX_BUFS;i++){ 
         if(tag && tag==_armci_buf_state->buf[i].id.tag)
           _armci_buf_complete_index(i,1); 
       }
       for(i=0;i<MAX_SMALL_BUFS;i++){ 
         if(tag && tag==_armci_buf_state->smallbuf[i].id.tag)
           _armci_buf_complete_index(i+MAX_BUFS,1); 
       }
#else
       {
	 const int smavail=_armci_buf_state->smavail;
	 const int avail = _armci_buf_state->avail;
	 const int startsmall = (smavail-MAX_BUFS+1+MAX_SMALL_BUFS)%MAX_SMALL_BUFS+MAX_BUFS;
	 const int start =(avail+1+MAX_BUFS)%MAX_BUFS;
	 i = start;
	 do {
	   if(tag && tag==_armci_buf_state->buf[i].id.tag)
	     _armci_buf_complete_index(i,1);
	   i = (i+1+MAX_BUFS)%MAX_BUFS;
	 } while(i != start);

	 i = startsmall;
	 do {
	   if(tag && tag==_armci_buf_state->smallbuf[i-MAX_BUFS].id.tag)
	     _armci_buf_complete_index(i,1);
	   i = (i-MAX_BUFS+1+MAX_SMALL_BUFS)%MAX_SMALL_BUFS+MAX_BUFS;
	 } while(i != startsmall);
       }
#endif
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

#ifdef VAPI
/* this will work for vapi as there is no get pipeline enabled in vapi
 * with get pipeline, this will break very badly
 */
void _armci_buf_update_scatter_count(int id)
{
int i,num,last,first;
    for(i=0;i<MAX_BUFS;i++){
# ifdef BUF_EXTRA_FIELD_T
       if(_armci_buf_state->table[i].op==GET){
         request_header_t *msginfo;
	 msginfo = (request_header_t*)_armci_buf_state->buf[i].buffer; 
         if(msginfo->pinned && msginfo->bypass && msginfo->format == STRIDED){
           num = *(int *)((char *)msginfo+msginfo->bytes); 
           last = *(int *)((char *)msginfo+msginfo->bytes+sizeof(int));
           first = last - num+1;
           if(first < 0 )first+=DSCRID_SCATTERCLIENT_END-DSCRID_SCATTERCLIENT-1;
           if(id == first && num!=0){
             *(int *)((char *)msginfo+msginfo->bytes) = (--num);
             return;
           }
         }
       }
# endif
    }
    for(i=MAX_BUFS;i<MAX_BUFS+MAX_SMALL_BUFS;i++){
# ifdef BUF_EXTRA_FIELD_T
       if(_armci_buf_state->table[i].op==GET){
         request_header_t *msginfo = (request_header_t*)_armci_buf_state->smallbuf[i-MAX_BUFS].buffer;
         if(msginfo->pinned && msginfo->bypass && msginfo->format == STRIDED){
           num = *(int *)((char *)msginfo+msginfo->bytes); 
           last = *(int *)((char *)msginfo+msginfo->bytes+sizeof(int));
           first = last - num+1;
           if(first < 0 )first+=DSCRID_SCATTERCLIENT_END-DSCRID_SCATTERCLIENT-1;
           if(id == first && num!=0){
             *(int *)((char *)msginfo+msginfo->bytes) = (--num);
             return;
           }
         }
       }
# endif
    }

}
#endif

#endif

/* function to return bufinfo, given buf tag */
BUF_INFO_T *_armci_tag_to_bufinfo(msg_tag_t tag) {
    int idx;

    for (idx=0; idx < MAX_BUFS; idx++)
        if (EQ_TAGS(_armci_buffers[idx].id.tag, tag)) break;

    if (idx == MAX_BUFS) {/* not found is regular buffers */
        for (idx = 0; idx < MAX_SMALL_BUFS; idx++)
            if (EQ_TAGS(_armci_smbuffers[idx].id.tag, tag)) break;
        if (idx == MAX_SMALL_BUFS) /* not found at all */
#ifdef SOCKETS
            armci_die("_armci_tag_to_bufinfo: bad tag",tag);
#else
            armci_die("_armci_tag_to_bufinfo: bad tag",0);
#endif

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


