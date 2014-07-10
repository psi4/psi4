#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: request.c,v 1.74.2.11 2007-10-18 06:09:37 d3h325 Exp $ */
#include "armcip.h"
#include "request.h"
#include "memlock.h"
#include "shmem.h"
#include "copy.h"
#include "gpc.h"
#include <stdio.h>
#include <signal.h>

#define DEBUG_ 0
#define DEBUG_MEM 0

#if 0
#   define MARK_ENTER(func_) { fprintf(stdout, "ENTERING %s\n", func_); fflush(stdout); }
#   define MARK_EXIT(func_) { fprintf(stdout, "EXITING %s\n", func_); fflush(stdout); }
#else
#   define MARK_ENTER(func_)
#   define MARK_EXIT(func_)
#endif

#if 0
#   define PRNDBG3(m,a1,a2,a3) \
               fprintf(stderr,"DBG %d: " m,armci_me,a1,a2,a3);fflush(stderr)
#   define PRNDBG(m) PRNDBG3(m,0,0,0)
#   define PRNDBG1(m,a1) PRNDBG3(m,a1,0,0)
#   define PRNDBG2(m,a1,a2) PRNDBG3(m,a1,a2,0)
#else
#   define PRNDBG(m)
#   define PRNDBG1(m,a1)
#   define PRNDBG2(m,a1,a2)
#   define PRNDBG3(m,a1,a2,a3)
#endif


#if !defined(GM) && !defined(VIA) && !defined(LAPI) &&!defined(VAPI)
  double _armci_rcv_buf[MSG_BUFLEN_DBL];
  double _armci_snd_buf[MSG_BUFLEN_DBL]; 
  char* MessageSndBuffer = (char*)_armci_snd_buf;
  char* MessageRcvBuffer = (char*)_armci_rcv_buf;
#endif


#define MAX_EHLEN 248
#define ADDBUF(buf,type,val) *(type*)(buf) = (val); (buf) += sizeof(type)
#define GETBUF(buf,type,var) (var) = *(type*)(buf); (buf) += sizeof(type)

#define ALLIGN8(buf){size_t _adr=(size_t)(buf); \
                    _adr>>=3; _adr<<=3; _adr+=8; (buf) = (char*)_adr; }

#ifndef CLN
#   define CLN 1
#endif
#ifndef SERV
#   define SERV 2
#endif

/*******************Routines to handle completion descriptor******************/
/*\
 *Following the the routines to fill a completion descriptor, if necessary
 *copy the data to destination based on completion descriptor
 *NOTE, THE FOLLOWING ROUTINES ARE FOR CLIENTS ONLY
\*/


/*\Routine to complete a vector request, data is in buf and descriptor in dscr
\*/
extern int armci_direct_vector_get(request_header_t *msginfo , armci_giov_t darr[], int len, int proc);
static void armci_complete_vector_get(armci_giov_t darr[],int len,void *buf)
{
int proc;
request_header_t *msginfo = (request_header_t*) buf;
    proc = msginfo->to;
#if defined(USE_SOCKET_VECTOR_API)
    armci_direct_vector_get(msginfo, darr, len, proc);
#else
    armci_rcv_vector_data(proc, msginfo, darr, len);
#endif
    FREE_SEND_BUFFER(buf);
}






/*\ Routine called from buffers.c to complete a request for which the buffer was
 *  used for, so that the buffer can be reused.
\*/
void armci_complete_req_buf(BUF_INFO_T *info, void *buffer)
{
request_header_t *msginfo = (request_header_t*) buffer;
    ARMCI_PR_DBG("enter",0);
    if(info->protocol==0)return;
    else if(info->protocol==SDSCR_IN_PLACE){
       char *dscr = info->dscr;
       void *loc_ptr;
       int stride_levels;
       int *loc_stride_arr,*count;

       loc_ptr = *(void**)dscr;           dscr += sizeof(void*);
       stride_levels = *(int*)dscr;       dscr += sizeof(int);
       loc_stride_arr = (int*)dscr;       dscr += stride_levels*sizeof(int);
       count = (int*)dscr;
       if(0 || DEBUG_){ 
         if(armci_me==0){
           printf("\n%d:extracted loc_ptr=%p, stridelevels=%d\n",armci_me,
                  loc_ptr,stride_levels);
           fflush(stdout);
         }
       }
       armci_rcv_strided_data(msginfo->to, msginfo, msginfo->datalen, loc_ptr,
                                stride_levels,loc_stride_arr,count);
       FREE_SEND_BUFFER(msginfo);
    }
    else if(info->protocol==VDSCR_IN_PLACE || info->protocol==VDSCR_IN_PTR){
       char *dscr;
       int len,i;
       if(info->protocol==VDSCR_IN_PLACE){
               dscr = info->dscr;
               //printf("\n%d:vdscr in place\n",armci_me);
       }
       else {
               dscr = info->ptr.dscrbuf;
               //printf("\n%d:vdscr in buf\n",armci_me);
       }
       GETBUF(dscr, long ,len);
       {
         armci_giov_t *darr;
         darr = (armci_giov_t *)malloc(sizeof(armci_giov_t)*len);
         if(!darr)armci_die("malloc in complete_req_buf failed",len);
         for(i = 0; i< len; i++){
           int parlen, bytes;
           GETBUF(dscr, int, parlen);
           GETBUF(dscr, int, bytes);
           darr[i].ptr_array_len = parlen;
           darr[i].bytes = bytes;
           if(msginfo->operation==GET)darr[i].dst_ptr_array=(void **)dscr;
           else darr[i].src_ptr_array=(void **)dscr;
           dscr+=sizeof(void *)*parlen;
         }
         if (msginfo->operation==GET) armci_complete_vector_get(darr,len,buffer);
       }
    }
    else
       armci_die("armci_complete_req_buf,protocol val invalid",info->protocol);
    ARMCI_PR_DBG("exit",0);
}

extern long x_net_offset(void *,int);
/*\ save a part of strided descriptor needed to complete request

rmo: it seems as if save_

\*/
void armci_save_strided_dscr(char **bptr, void *rem_ptr,int rem_stride_arr[],
                             int count[], int stride_levels,int is_nb,int proc)
{
int i;
char *bufptr=*bptr;
BUF_INFO_T *info=NULL;
long network_offset,tmpoffset;
    ARMCI_PR_DBG("enter",0);

  # ifdef PORTALS_UNRESOLVED
    if(!is_nb){
      network_offset=x_net_offset(rem_ptr,proc);
      if(DEBUG_){printf("\n%d:rem_ptr=%p offset=%d newrem=%p",armci_me,rem_ptr,network_offset,(char *)rem_ptr+network_offset);fflush(stdout);}
      rem_ptr = (char *)rem_ptr+network_offset;
    }
  # endif

    if(is_nb){
       info=BUF_TO_BUFINFO(*bptr);
       bufptr = (info->dscr);
    }
    *(void**)bufptr = rem_ptr;         bufptr += sizeof(void*);
    *(int*)bufptr = stride_levels;     bufptr += sizeof(int);
    for(i=0;i<stride_levels;i++)((int*)bufptr)[i] = rem_stride_arr[i];
    bufptr += stride_levels*sizeof(int);
    for(i=0;i< stride_levels+1;i++)((int*)bufptr)[i] = count[i];
    bufptr += (1+stride_levels)*sizeof(int);
    if((0 || DEBUG_) && is_nb){
      bufptr = (info->dscr);
      if(armci_me==0)
        printf("\n%d:rem_ptr %p=%p stride_levels %d=%d\n",armci_me,
                *(void**)bufptr,rem_ptr,
                *(int*)(bufptr + sizeof(void*)),stride_levels);
    }
    /*remote_strided expects the pointer to point to the end of descr hence..*/
    if(is_nb)
       info->protocol=SDSCR_IN_PLACE;
    else
       *bptr=bufptr;
    ARMCI_PR_DBG("exit",0);

}


/*\ save a part of vector descriptor needed to complete request
\*/
void armci_save_vector_dscr(char **bptr,armci_giov_t darr[],int len,
                            int op,int is_nb, int proc)
{
int i,size=sizeof(int);
BUF_INFO_T *info;
char *buf,*bufptr=*bptr;
void *rem_ptr;
long offst;
    ARMCI_PR_DBG("enter",0);
    if(is_nb){    
       for(i=0;i<len;i++){
         size+=(2*sizeof(int)+darr[i].ptr_array_len * sizeof(void*));
       }
       info=BUF_TO_BUFINFO(bufptr);
       /*if descr fits in available buffer, use it else do malloc */
       if(size<=UBUF_LEN){
         buf = info->dscr;
         info->protocol=VDSCR_IN_PLACE;
       }
       else {
         info->ptr.dscrbuf = (void *)malloc(size);
         buf = (char *)info->ptr.dscrbuf;
         info->protocol=VDSCR_IN_PTR;
       }
    }
    else
       buf=bufptr;
      
    ADDBUF(buf,long,len); /* number of sets */
    for(i=0;i<len;i++){
        int j;
        ADDBUF(buf,int,darr[i].ptr_array_len); /* number of elements */
        ADDBUF(buf,int,darr[i].bytes);         /* sizeof element */
        if(op==GET) {
          if(is_nb){
            rem_ptr = darr[i].dst_ptr_array;
          }
          else {
          # ifdef PORTALS_UNRESOLVED
            for(j=0;j<darr[i].ptr_array_len;j++){
              offst=x_net_offset(darr[i].src_ptr_array[j],proc);
              darr[i].src_ptr_array[j]= (char*)darr[i].src_ptr_array[j]+offst;
            }
          # endif
            rem_ptr = darr[i].src_ptr_array;
          }
        }
        else {
        # ifdef PORTALS_UNRESOLVED
          for(j=0;j<darr[i].ptr_array_len;j++){
            offst=x_net_offset(darr[i].dst_ptr_array[j],proc);
            darr[i].dst_ptr_array[j]= (char*)darr[i].dst_ptr_array[j]+offst;
          } 
        # endif
          rem_ptr = darr[i].dst_ptr_array;
        }
        armci_copy(rem_ptr,buf, darr[i].ptr_array_len * sizeof(void*));
        buf += darr[i].ptr_array_len*sizeof(void*);
    }
    if(!is_nb)
       *bptr=buf;
    ARMCI_PR_DBG("exit",0);
}

/*\
 * If buf==null, set handle->bufid to val, else set it to the id of the buf
\*/
void armci_set_nbhandle_bufid(armci_ihdl_t nb_handle,char *buf,int val)
{
BUF_INFO_T *info;
    if(buf){
       info = BUF_TO_BUFINFO(buf);
       val = info->bufid;
    }
    nb_handle->bufid = val; 
} 

/**************End--Routines to handle completion descriptor******************/


/*\ send request to server to LOCK MUTEX
\*/
void armci_rem_lock(int mutex, int proc, int *ticket)
{
request_header_t *msginfo;
int *ibuf;
int bufsize = sizeof(request_header_t)+sizeof(int);

    msginfo = (request_header_t*)GET_SEND_BUFFER(bufsize,LOCK,proc);
    bzero(msginfo,sizeof(request_header_t));

    msginfo->datalen = sizeof(int);
    msginfo->dscrlen = 0;
    msginfo->from  = armci_me;
    msginfo->to    = proc;
    msginfo->operation = LOCK;
    msginfo->format  = mutex;
    msginfo->bytes = msginfo->datalen + msginfo->dscrlen;

    ibuf = (int*)(msginfo+1);
    *ibuf = mutex;

    armci_send_req(proc, msginfo, bufsize, 0);

    /* receive ticket from server */
    *ticket = *(int*)armci_rcv_data(proc,msginfo,0);
    FREE_SEND_BUFFER(msginfo);

    if(DEBUG_)fprintf(stderr,"%d receiving ticket %d\n",armci_me, *ticket);
}




void armci_server_lock(request_header_t *msginfo)
{
int *ibuf = (int*)(msginfo+1);
int proc  = msginfo->from;
int mutex;
int ticket;
    ARMCI_PR_DBG("enter",0);

    mutex = *(int*)ibuf;

    /* acquire lock on behalf of requesting process */
    ticket = armci_server_lock_mutex(mutex, proc, msginfo->tag);

    if(ticket >-1){
       /* got lock */
       msginfo->datalen = sizeof(int);
       armci_send_data(msginfo, &ticket);
    }
    ARMCI_PR_DBG("exit",0);
}


/*\ send request to server to UNLOCK MUTEX
\*/
void armci_rem_unlock(int mutex, int proc, int ticket)
{
request_header_t *msginfo;
int *ibuf;
int bufsize = sizeof(request_header_t)+sizeof(ticket);

    msginfo = (request_header_t*)GET_SEND_BUFFER(bufsize,UNLOCK,proc);
    bzero(msginfo,sizeof(request_header_t));

    msginfo->dscrlen = msginfo->bytes = sizeof(ticket); 
    msginfo->datalen = 0; 
    msginfo->from  = armci_me;
    msginfo->to    = proc;
    msginfo->operation = UNLOCK;
    msginfo->format  = mutex;
    ibuf = (int*)(msginfo+1);
    *ibuf = ticket;

    if(DEBUG_)fprintf(stderr,"%d sending unlock\n",armci_me);
    armci_send_req(proc, msginfo, bufsize,0);
}
    


/*\ server unlocks mutex and passes lock to the next waiting process
\*/
void armci_server_unlock(request_header_t *msginfo, char* dscr)
{
    int ticket = *(int*)dscr;
    int mutex  = msginfo->format;
    int proc   = msginfo->to;
    int waiting;
    
    waiting = armci_server_unlock_mutex(mutex,proc,ticket,&msginfo->tag);

    if(waiting >-1){ /* -1 means that nobody is waiting */

       ticket++;
       /* pass ticket to the waiting process */
       msginfo->from = waiting;
       msginfo->datalen = sizeof(ticket);
       armci_send_data(msginfo, &ticket);

    }
}

void armci_unlock_waiting_process(msg_tag_t tag, int proc, int ticket)
{
request_header_t header;
request_header_t *msginfo = &header;

  msginfo->datalen = sizeof(int);
  msginfo->tag     = tag;
  msginfo->from      = proc;
  msginfo->to    = armci_me;
  armci_send_data(msginfo, &ticket); 
}

void * armci_server_ptr(int id){
char *buf;
int bufsize = sizeof(int);
request_header_t *msginfo = (request_header_t*)GET_SEND_BUFFER(bufsize,ATTACH,armci_me);
  bzero(msginfo,sizeof(request_header_t));
  msginfo->from  = armci_me;
  msginfo->to    = SERVER_NODE(armci_clus_me);
  msginfo->dscrlen   = 0;
  msginfo->datalen = sizeof(int);
  msginfo->operation =  ATTACH;
  msginfo->bytes = msginfo->dscrlen+ msginfo->datalen;
  armci_copy(&id, msginfo +1, sizeof(int));
  if(DEBUG_MEM){
    printf("\n%d:attach req:sending id %d \n",armci_me,id);fflush(stdout);
  }
  armci_send_req(armci_master, msginfo, bufsize,0);
  buf= armci_rcv_data(armci_master,msginfo,sizeof(void *));/* receive response */
  if(DEBUG_MEM){
    printf("\n%d:attach req:got %p \n",armci_me,buf);fflush(stdout);
  }
  FREE_SEND_BUFFER(msginfo);
  ARMCI_PR_DBG("exit",0);
  return (void *)buf;

}

/*\ control message to the server, e.g.: ATTACH to shmem, return ptr etc.
\*/
void armci_serv_attach_req(void *info, int ilen, long size, void* resp,int rlen)
{
char *buf;
    ARMCI_PR_DBG("enter",0);
int bufsize = 2*sizeof(request_header_t)+ilen + sizeof(long)+sizeof(rlen);
long *idlist=(long *)info;
request_header_t *msginfo = (request_header_t*)GET_SEND_BUFFER(bufsize,ATTACH,armci_me);
    bzero(msginfo,sizeof(request_header_t));

    msginfo->from  = armci_me;
    msginfo->to    = SERVER_NODE(armci_clus_me);
    msginfo->dscrlen   = ilen;
    msginfo->datalen = sizeof(long)+sizeof(int);
    msginfo->operation =  ATTACH;
    msginfo->bytes = msginfo->dscrlen+ msginfo->datalen;

    armci_copy(info, msginfo +1, ilen);
    if(DEBUG_MEM){printf("\n%d:sending idlist+1 %d, size %d, idlist[0] %d, idlist[1] %d\n",armci_me,idlist+1,size,idlist[0],idlist[1]);}
    buf = ((char*)msginfo) + ilen + sizeof(request_header_t);
    *((long*)buf) =size;
    *(int*)(buf+ sizeof(long)) = rlen;
    armci_send_req(armci_master, msginfo, bufsize,0);
    if(rlen){
      buf= armci_rcv_data(armci_master, msginfo,rlen);  /* receive response */
      bcopy(buf, resp, rlen);
      FREE_SEND_BUFFER(msginfo);

      if(DEBUG_MEM){printf("%d:client attaching got ptr=%p %d bytes\n",armci_me,buf,rlen);
         fflush(stdout);
      }
    }
    ARMCI_PR_DBG("exit",0);
}


/*\ server initializes its copy of the memory lock data structures
\*/
static void server_alloc_memlock(void *ptr_myclus)
{
int i;

    /* for protection, set pointers for processes outside local node NULL */
    memlock_table_array = calloc(armci_nproc,sizeof(void*));
    if(!memlock_table_array) armci_die("malloc failed for ARMCI lock array",0);

    /* set pointers for processes on local cluster node
     * ptr_myclus - corresponds to the master process
     */
    for(i=0; i< armci_clus_info[armci_clus_me].nslave; i++){
        memlock_table_array[armci_master +i] = ((char*)ptr_myclus)
                + MAX_SLOTS*sizeof(memlock_t)*i;
    }

    /* set pointer to the use flag */
#ifdef MEMLOCK_SHMEM_FLAG
    armci_use_memlock_table = (int*) (MAX_SLOTS*sizeof(memlock_t) +
                      (char*) memlock_table_array[armci_clus_last]);
    
    if(DEBUG_)
      fprintf(stderr,"server initialized memlock %p\n",armci_use_memlock_table);
#endif
}


static int allocate_memlock=1;

/*\ server actions triggered by client request to ATTACH
\*/
void armci_server_ipc(request_header_t* msginfo, void* descr,
                      void* buffer, int buflen)
{
double *ptr;
long *idlist = (long*)descr;
long size = *(long*)buffer;
int rlen = *(int*)(sizeof(long)+(char*)buffer);
extern int **_armci_int_mutexes;
    ARMCI_PR_DBG("enter",0);
    if(size<0) armci_die("armci_server_ipc: size<0",(int)size);
    if(DEBUG_MEM)printf("\n%d:got idlist+1 %p, size %d, idlist[0] %d, idlist[1] %d",armci_me,idlist+1,size,idlist[0],idlist[1]);
    ptr=(double*)Attach_Shared_Region(idlist+1,size,idlist[0]);
    if(!ptr)armci_die("armci_server_ipc: failed to attach",0);
    /* provide data server with access to the memory lock data structures */
    if(allocate_memlock){
      allocate_memlock = 0;
      server_alloc_memlock(ptr);
    }
    if(_armci_int_mutexes==NULL){
      printf("unresolved portals external\n");
      abort();
    # ifdef PORTALS_UNRESOLVED
      extern int _armci_server_mutex_ready;
      extern void *_armci_server_mutex_ptr;
      if(_armci_server_mutex_ready){
        _armci_int_mutexes=(int **)_armci_server_mutex_ptr;
      }
    # endif
    }
   if(size>0)armci_set_mem_offset(ptr);

   if(msginfo->datalen != sizeof(long)+sizeof(int))
      armci_die("armci_server_ipc: bad msginfo->datalen ",msginfo->datalen);

   if(rlen==sizeof(ptr)){
     msginfo->datalen = rlen;
     armci_send_data(msginfo, &ptr);
   }
   else armci_die("armci_server_ipc: bad rlen",rlen);
   ARMCI_PR_DBG("exit",0);
}


/*\ send RMW request to server
\*/
void armci_rem_rmw(int op, void *ploc, void *prem, int extra, int proc)
{
request_header_t *msginfo;
char *buf;
void *buffer;
int bufsize = sizeof(request_header_t)+sizeof(long)+sizeof(void*);
long offst;

    ARMCI_PR_DBG("enter",0);
    msginfo = (request_header_t*)GET_SEND_BUFFER(bufsize,op,proc);
    bzero(msginfo,sizeof(request_header_t));

    msginfo->dscrlen = sizeof(void*);
    msginfo->from  = armci_me;
    msginfo->to    = proc;
    msginfo->operation = op;
    msginfo->datalen = sizeof(long);
  # ifdef PORTALS_UNRESOLVED
    offst=x_net_offset(prem,proc);
    prem = ((char *)prem+offst);
  # endif
    buf = (char*)(msginfo+1);
    ADDBUF(buf, void*, prem); /* pointer is shipped as descriptor */

    /* data field: extra argument in fetch&add and local value in swap */
    if(op==ARMCI_SWAP){
       ADDBUF(buf, int, *((int*)ploc));
    }else if(op==ARMCI_SWAP_LONG) {
       ADDBUF(buf, long, *((long*)ploc) );
       msginfo->datalen = sizeof(long);
    }else {
       ADDBUF(buf, int, extra);
    }

    msginfo->bytes   = msginfo->datalen+msginfo->dscrlen ;

    if(DEBUG_){
        printf("%d sending RMW request %d to %d\n",armci_me,op,proc);
        fflush(stdout);
    }
    armci_send_req(proc, msginfo, bufsize,0);
    buffer = armci_rcv_data(proc,msginfo,0);  /* receive response */

    if(op==ARMCI_FETCH_AND_ADD || op== ARMCI_SWAP)
        *(int*)ploc = *(int*)buffer;
    else
        *(long*)ploc = *(long*)buffer;

    FREE_SEND_BUFFER(msginfo);
    ARMCI_PR_DBG("exit",0);
}


/*\ server response to RMW 
\*/
void armci_server_rmw(request_header_t* msginfo,void* ptr, void* pextra)
{
long lold;
int iold;
void *pold=0;
int op = msginfo->operation;

    ARMCI_PR_DBG("enter",0);
     if(DEBUG_){
       printf("%d server: executing RMW from %d. op=%d pextra=%p\n",armci_me,msginfo->from, op, pextra);
        fflush(stdout);
     }
     if(msginfo->datalen != sizeof(long))
          armci_die2("armci_server_rmw: bad datalen=",msginfo->datalen,op);

     /* for swap operations *pextra has the  value to swap
      * for fetc&add it carries the increment argument
      */
     switch(op){
     case ARMCI_SWAP:
        iold = *(int*) pextra;
     case ARMCI_FETCH_AND_ADD:
        pold = &iold;
        break;

     case ARMCI_SWAP_LONG:
        lold = *(long*) pextra;
     case ARMCI_FETCH_AND_ADD_LONG:
        pold = &lold;
        break;

     default:
          armci_die("armci_server_rmw: bad operation code=",op);
     }

     armci_generic_rmw(op, pold, *(int**)ptr, *(int*) pextra, msginfo->to);

     armci_send_data(msginfo, pold);
    ARMCI_PR_DBG("exit",0);
}

extern int armci_direct_vector_snd(request_header_t *msginfo , armci_giov_t darr[], int len, int proc);
extern int armci_direct_vector(request_header_t *msginfo , armci_giov_t darr[], int len, int proc);
int armci_rem_vector(int op, void *scale, armci_giov_t darr[],int len,int proc,int flag, armci_ihdl_t nb_handle)
{
    char *buf,*buf0;
    request_header_t *msginfo;
    int bytes =0, s, slen=0;
    size_t adr;
    int bufsize = sizeof(request_header_t);
    int tag=0;

    if(nb_handle)tag=nb_handle->tag;

    /* compute size of the buffer needed */
    for(s=0; s<len; s++){
        bytes   += darr[s].ptr_array_len * darr[s].bytes; /* data */
        bufsize += darr[s].ptr_array_len *sizeof(void*)+2*sizeof(int); /*descr*/
    }

    bufsize += bytes + sizeof(long) +2*sizeof(double) +8; /*+scale+allignment*/

    buf = buf0= GET_SEND_BUFFER(bufsize,op,proc);
    msginfo = (request_header_t*)buf;
    bzero(msginfo,sizeof(request_header_t));

/*     printf("%d:: rem_vector. len=%d. ptr_len[len-1]=%d bytes[len-1]=%d bufsize=%d\n",  */
/* 	   armci_me, len, darr[len-1].ptr_array_len, darr[len-1].bytes,bufsize); */
/*     fflush(stdout); */


    if(nb_handle){
   /*   INIT_SENDBUF_INFO(nb_handle,buf,op,proc); redundant -- see armci_rem_strided */
      _armci_buf_set_tag(buf,nb_handle->tag,0);
      if(nb_handle->bufid == NB_NONE)
        armci_set_nbhandle_bufid(nb_handle,buf,0);
    }

    buf += sizeof(request_header_t);

    /* fill vector descriptor */
    armci_save_vector_dscr(&buf,darr,len,op,0,proc);

    /* align buf for doubles (8-bytes) before copying data */
    adr = (size_t)buf;
    adr >>=3;
    adr <<=3;
    adr +=8;
    buf = (char*)adr;

    msginfo->ehlen = 0;

    /* fill message header */
    msginfo->dscrlen = buf - buf0 - sizeof(request_header_t);
    msginfo->from  = armci_me;
    msginfo->to    = proc;
    msginfo->operation  = op;
    msginfo->format  = VECTOR;
    msginfo->datalen = bytes;

    /* put scale for accumulate */
    switch(op){
    case ARMCI_ACC_INT:
               *(int*)buf = *(int*)scale; slen= sizeof(int); break;
    case ARMCI_ACC_DCP:
               ((double*)buf)[0] = ((double*)scale)[0];
               ((double*)buf)[1] = ((double*)scale)[1];
               slen=2*sizeof(double);break;
    case ARMCI_ACC_DBL:
               *(double*)buf = *(double*)scale; slen = sizeof(double); break;
    case ARMCI_ACC_CPL:
               ((float*)buf)[0] = ((float*)scale)[0];
               ((float*)buf)[1] = ((float*)scale)[1];
               slen=2*sizeof(float);break;
    case ARMCI_ACC_FLT:
               *(float*)buf = *(float*)scale; slen = sizeof(float); break;
    default: slen=0;
    }
    buf += slen;
    msginfo->datalen += slen;
    msginfo->bytes = msginfo->datalen+msginfo->dscrlen;


    /* for put and accumulate copy data into buffer */
    if(op != GET){
/*       fprintf(stderr,"sending %lf\n",*(double*)darr[0].src_ptr_array[0]);*/
       armci_vector_to_buf(darr, len, buf);
    }

    armci_send_req(proc, msginfo, bufsize,tag);
    /*x_buf_send_complete(buf0);*/

    if(nb_handle && op==GET) armci_save_vector_dscr(&buf0,darr,len,op,1,proc);
    if(op == GET&& !nb_handle){
       armci_complete_vector_get(darr,len,msginfo);
    }

    return 0;
}

#define CHUN_ (8*8096)
#define CHUN 200000

/*\ client version of remote strided operation
\*/
int armci_rem_strided(int op, void* scale, int proc,
                       void *src_ptr, int src_stride_arr[],
                       void* dst_ptr, int dst_stride_arr[],
                       int count[], int stride_levels,
                       ext_header_t *h, int flag,armci_ihdl_t nb_handle)
{
    char *buf, *buf0;
    request_header_t *msginfo;
    int  i, slen=0, bytes;
    void *rem_ptr;
    int  *rem_stride_arr;
    int bufsize = sizeof(request_header_t);
    int ehlen =0;
    msg_tag_t msg_tag;
    int tag=0;

    /* we send ext header only for last chunk */
#if 0
    if(h)  ehlen = h->len;
#else
    if(h) if(h->last)  ehlen = h->len;
#endif
    if(ehlen>MAX_EHLEN || ehlen <0) 
       armci_die2("armci_rem_strided ehlen out of range",MAX_EHLEN,ehlen);
    /* calculate size of the buffer needed */
    for(i=0, bytes=1;i<=stride_levels;i++)bytes*=count[i];
    bufsize += bytes+sizeof(void*)+2*sizeof(int)*(stride_levels+1) +ehlen
               +2*sizeof(double) + 16; /* +scale+alignment */

    if (flag){
      printf("%d: flag=%d\n",armci_me,flag);
      if(op==GET)bufsize -=bytes;
    }

    buf = buf0= GET_SEND_BUFFER((bufsize),op,proc);
    msginfo = (request_header_t*)buf;
    bzero(msginfo,sizeof(request_header_t));


    if(nb_handle)
#ifdef ACC_SMP
	 if(!ARMCI_ACC(op))
#endif
    {
    // printf("%s: non-blocking ops not yet supported\n",Portals_ID());
    // abort();
/*    INIT_SENDBUF_INFO(nb_handle,buf,op,proc); same as _armci_buf_set_tag, why here? */
     _armci_buf_set_tag(buf,nb_handle->tag,0);
     if(nb_handle->bufid == NB_NONE)
        armci_set_nbhandle_bufid(nb_handle,buf,0);
     tag = nb_handle->tag;
    }

    if(op == GET){
       rem_ptr = src_ptr;
       rem_stride_arr = src_stride_arr;
    }else{
       rem_ptr = dst_ptr;
       rem_stride_arr = dst_stride_arr;
    }

    msginfo->datalen=bytes;

    /* fill strided descriptor */
    buf += sizeof(request_header_t);
    /*this function fills the dscr into buf and also moves the buf ptr to the
      end of the dscr*/
    armci_save_strided_dscr(&buf,rem_ptr,rem_stride_arr,count,stride_levels,0,proc);

    /* align buf for doubles (8-bytes) before copying data */
    ALLIGN8(buf);

    /* fill message header */
    msginfo->from   = armci_me;
    msginfo->to     = proc;
    msginfo->format = STRIDED;
    msginfo->operation  = op;

    /* put scale for accumulate */
    switch(op){
    case ARMCI_ACC_INT:
               *(int*)buf = *(int*)scale; slen= sizeof(int); break;
    case ARMCI_ACC_DCP:
               ((double*)buf)[0] = ((double*)scale)[0];
               ((double*)buf)[1] = ((double*)scale)[1];
               slen=2*sizeof(double);break;
    case ARMCI_ACC_DBL:
               *(double*)buf = *(double*)scale; slen = sizeof(double); break;
    case ARMCI_ACC_CPL:
               ((float*)buf)[0] = ((float*)scale)[0];
               ((float*)buf)[1] = ((float*)scale)[1];
               slen=2*sizeof(float);break;
    case ARMCI_ACC_FLT:
               *(float*)buf = *(float*)scale; slen = sizeof(float); break;
    case ARMCI_ACC_LNG:
               *(long*)buf = *(long*)scale; slen = sizeof(long); break;
    default: slen=0;
    }

    /*
	if(ARMCI_ACC(op))printf("%d client len=%d alpha=%lf data=%lf,%lf\n",
	     armci_me, buf-(char*)msginfo,((double*)buf)[0],*((double*)src_ptr),             ((double*)buf)[1]);
    */

    buf += slen;

    /**** add extended header *******/
    if(ehlen){
       bcopy(h->exthdr,buf,ehlen);
       i = ehlen%8; ehlen += (8-i); /* make sure buffer is still alligned */
       buf += ehlen;
    }

    msginfo->ehlen  = ehlen;
    msginfo->dscrlen = buf - buf0 - sizeof(request_header_t);
    msginfo->bytes = msginfo->datalen+msginfo->dscrlen;

    if(op == GET){
    /*
      if(nb_handle) {
         printf("%s rem_strided: nb gets not yet available\n",Portals_ID());
         abort();
      }
    */
      armci_send_req(proc, msginfo, bufsize,tag);
      armci_save_strided_dscr(&buf0,dst_ptr,dst_stride_arr,count,
                                 stride_levels,1,proc);
      
      if(!nb_handle){
        armci_rcv_strided_data(proc, msginfo, msginfo->datalen,
                              dst_ptr, stride_levels, dst_stride_arr, count);
        FREE_SEND_BUFFER(msginfo);
      }
    } else {
       /* for put and accumulate send data */
       armci_send_strided(proc,msginfo, buf,
                          src_ptr, stride_levels, src_stride_arr, count,tag);
    }

    return 0;
}


void armci_process_extheader(request_header_t *msginfo, char *dscr, char* buf, int buflen)
{
 armci_flag_t *h;
 int *flag;

   h = (armci_flag_t*)(dscr + msginfo->dscrlen - msginfo->ehlen);
#if 0
   if(msginfo->ehlen)printf("%d:server from=%d len=%d: ptr=%p val=%d\n",armci_me,msginfo->from, msginfo->ehlen,h->ptr,h->val);
   fflush(stdout);
#endif
   flag = (int*)(h->ptr);
   *flag = h->val;
}

void armci_server(request_header_t *msginfo, char *dscr, char* buf, int buflen)
{
int  buf_stride_arr[MAX_STRIDE_LEVEL+1];
int  *loc_stride_arr,slen; 
int  *count, stride_levels;
void *buf_ptr, *loc_ptr;
void *scale;
char *dscr_save = dscr;
int  rc, i,proc;
int stat;
      
    ARMCI_PR_DBG("enter",msginfo->datalen);fflush(stdout);
    /*return if using readv/socket for put*/
    if(msginfo->operation==PUT && msginfo->datalen==0){
      if(msginfo->ehlen) /* process extra header if available */
         armci_process_extheader(msginfo, dscr, buf, buflen);
      return;
    }
    
    /* unpack descriptor record */
    loc_ptr = *(void**)dscr;           dscr += sizeof(void*);
    stride_levels = *(int*)dscr;       dscr += sizeof(int);
    loc_stride_arr = (int*)dscr;       dscr += stride_levels*sizeof(int);
    count = (int*)dscr;

    /* compute stride array for buffer */
    buf_stride_arr[0]=count[0];
    for(i=0; i< stride_levels; i++)
        buf_stride_arr[i+1]= buf_stride_arr[i]*count[i+1];

    /* get scale for accumulate, adjust buf to point to data */
    switch(msginfo->operation){
    case ARMCI_ACC_INT:     slen = sizeof(int); break;
    case ARMCI_ACC_DCP:     slen = 2*sizeof(double); break;
    case ARMCI_ACC_DBL:     slen = sizeof(double); break;
    case ARMCI_ACC_CPL:     slen = 2*sizeof(float); break;
    case ARMCI_ACC_FLT:     slen = sizeof(float); break;
    case ARMCI_ACC_LNG:     slen = sizeof(long); break;
	default:				slen=0;
    }

    scale = dscr_save+ (msginfo->dscrlen - slen -msginfo->ehlen);
/*
    if(ARMCI_ACC(msginfo->operation))
      fprintf(stderr,"%d in server len=%d slen=%d alpha=%lf data=%lf\n", 
               armci_me, msginfo->dscrlen, slen, *(double*)scale,*(double*)buf);
*/

    buf_ptr = buf; /*  data in buffer */

    proc = msginfo->to;

    if(msginfo->operation == GET){
       armci_send_strided_data(proc, msginfo, buf,
                               loc_ptr, stride_levels, loc_stride_arr, count);
       /* fprintf(stderr, "GET response sent with tag: %d\n, msginfo->tag",
          msginfo->tag); */
    } else{
       if((rc = armci_op_strided(msginfo->operation, scale, proc,
               buf_ptr, buf_stride_arr, loc_ptr, loc_stride_arr,
               count, stride_levels, 1,NULL)))
               armci_die("server_strided: op from buf failed",rc);
    }

    if(msginfo->ehlen) /* process extra header if available */
         armci_process_extheader(msginfo, dscr_save, buf, buflen);
    ARMCI_PR_DBG("exit",0);
}


void armci_server_vector( request_header_t *msginfo, 
                          char *dscr, char* buf, int buflen)
{
    int  proc;
    long  len;
    void *scale;
    int  i,s;
    char *sbuf = buf;
    if(msginfo->operation==PUT && msginfo->datalen==0)return;/*return if using readv/socket for put*/
    /* unpack descriptor record */
    GETBUF(dscr, long ,len);
    
    /* get scale for accumulate, adjust buf to point to data */
    scale = buf;
    switch(msginfo->operation){
    case ARMCI_ACC_INT:     buf += sizeof(int); break;
    case ARMCI_ACC_DCP:     buf += 2*sizeof(double); break;
    case ARMCI_ACC_DBL:     buf += sizeof(double); break;
    case ARMCI_ACC_CPL:     buf += 2*sizeof(float); break;
    case ARMCI_ACC_FLT:     buf += sizeof(float); break;
    }

    proc = msginfo->to;

    /*fprintf(stderr,"scale=%lf\n",*(double*)scale);*/
    /* execute the operation */

    switch(msginfo->operation) {
    case GET:
/*        fprintf(stderr, "%d:: Got a vector message!!\n", armci_me); */
      if(msginfo->ehlen) {
	armci_die("Unexpected vector message with non-zero ehlen. GPC call?",
		   msginfo->ehlen);
      }
      else {
	for(i = 0; i< len; i++){
	  int parlen, bytes;
	  void **ptr;
	  GETBUF(dscr, int, parlen);
	  GETBUF(dscr, int, bytes);
	  /*        fprintf(stderr,"len=%d bytes=%d parlen=%d\n",len,bytes,parlen);*/
	  ptr = (void**)dscr; dscr += parlen*sizeof(char*);
	  for(s=0; s< parlen; s++){
	    armci_copy(ptr[s], buf, bytes);
	    buf += bytes;
	  }
	}
/*     fprintf(stderr,"%d:: VECTOR GET. server sending buffer %p datalen=%d\n",armci_me, sbuf, msginfo->datalen); */
	armci_send_data(msginfo, sbuf);
      }
      break;

    case PUT:

/*    fprintf(stderr,"received in buffer %lf\n",*(double*)buf);*/
      for(i = 0; i< len; i++){
        int parlen, bytes;
        void **ptr;
        GETBUF(dscr, int, parlen);
        GETBUF(dscr, int, bytes);
        ptr = (void**)dscr; dscr += parlen*sizeof(char*);
        for(s=0; s< parlen; s++){
/*
          armci_copy(buf, ptr[s], bytes);
*/
          bcopy(buf, ptr[s], (size_t)bytes);
          buf += bytes;
        }
      }
      break;

     default:

      /* this should be accumulate */
      if(!ARMCI_ACC(msginfo->operation))
               armci_die("v server: wrong op code",msginfo->operation);

/*      fprintf(stderr,"received first=%lf last =%lf in buffer\n",*/
/*                     *((double*)buf),((double*)buf)[99]);*/

      for(i = 0; i< len; i++){
        int parlen, bytes;
        void **ptr;
        GETBUF(dscr, int, parlen);
        GETBUF(dscr, int, bytes);
        ptr = (void**)dscr; dscr += parlen*sizeof(char*);
        armci_lockmem_scatter(ptr, parlen, bytes, proc); 
        for(s=0; s< parlen; s++){
          armci_acc_2D(msginfo->operation, scale, proc, buf, ptr[s],
                       bytes, 1, bytes, bytes, 0);
          buf += bytes;
        }
        ARMCI_UNLOCKMEM(proc);
      }
    }
}
