#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "armcip.h"
#include "request.h"
#include "message.h"
#include "memlock.h"
#include "copy.h"
#include "gpc.h"
#include <stdio.h>
#include <assert.h>
#ifdef WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

#define DEBUG_ 0
#define DEBUG1 0

#ifndef SERV
#     define SERV 2
#endif

#ifdef SOCKETS
#   define EQ_TAGS(a_, b_) ((a_) == (b_))
#else
#   define EQ_TAGS(a_, b_) !memcmp(&(a_), &(b_), sizeof(a_))
#endif

int _armci_server_started=0;

extern active_socks_t *_armci_active_socks;

#ifdef ARMCI_CHECK_STATE
typedef struct sarns{
  int data;
  long data1;
  struct sarns *next; 
} sarnode;

sarnode **sarn_np=NULL;

sarnode * sarlist_add(int pr, int i,long j)
{
sarnode **p = &sarn_np[pr];
  sarnode *n = (sarnode *)malloc(sizeof(sarnode));
  assert(n != NULL);
   
  n->next = *p; 
  *p = n; 
  n->data = i;
  n->data1 = j;
  return *p;
}
 
void sarlist_remove(sarnode **p)
{
  if(*p != NULL){
    sarnode *n = *p;
    *p = (*p)->next;
    free(n);
  }
}
 
sarnode **sarlist_search(sarnode **n, long i)
{
  while (*n != NULL){
    if ((*n)->data == i){
      return n;
    }
    n = &(*n)->next;
  }
  return NULL;
}
 
void sarlist_print(int proc)
{
  sarnode *n =sarn_np[proc]; 
  if (n == NULL){
    /*printf("sarlist is empty\n");*/
  }
  while (n != NULL){
    printf("(%d):%d %d next=%d\n", armci_me,n->data,n->data1,(n->next==NULL)?0:1);
    n = n->next;
  }
}
#endif 

/*\ client sends request to server
\*/
void armci_send_req(int proc, request_header_t* msginfo, int len,int tag)
{
int hdrlen = sizeof(request_header_t);
int bytes;

    ARMCI_PR_DBG("enter",0);
    if(msginfo->operation == GET) {
        if(msginfo->format==VECTOR && msginfo->ehlen > 0) {
            printf("%d [cp] unhandled condition in send_req for VECTOR and ehlen\n",armci_me);
            abort();
            bytes = msginfo->dscrlen + hdrlen + msginfo->datalen;
        } else {
            bytes = msginfo->dscrlen + hdrlen;
        }
    } else bytes = msginfo->bytes + hdrlen;

    if(0 || DEBUG_) {
       printf("%d: sending req %d (len=%d dscr=%d data=%d) to %d \n",
         armci_me, msginfo->operation, bytes,msginfo->dscrlen,
         msginfo->datalen,proc); fflush(stdout);
    }
    if(bytes > len) armci_die2("armci_send_req:buffer overflow",bytes,len);
 // msginfo->tag.data_ptr = (msginfo+1); // not really data, but dscr ptr
    armci_send_req_msg(proc,msginfo, bytes,tag);
    ARMCI_PR_DBG("exit",0);
}


/*\ client sends strided data + request to server
\*/
void armci_send_strided(int proc, request_header_t *msginfo, char *bdata,
                        void *ptr, int strides, int stride_arr[], int count[],int tag)
{
int hdrlen = sizeof(request_header_t);
int dscrlen = msginfo->dscrlen;
int datalen = msginfo->datalen;
int cluster = armci_clus_id(proc);
int bytes;
int i,na;
char *a;
double *tmp;
    ARMCI_PR_DBG("enter",0);
    bytes = msginfo->bytes + hdrlen;
    if(0){
       printf("%d:sending strided %d to(%d,%d,%d) bytes=%d dslen=%d dlen=%d,\n",
                armci_me, msginfo->operation, msginfo->to,
                cluster, proc, bytes, dscrlen, datalen); fflush(stdout);
    }
    armci_write_strided(ptr, strides, stride_arr, count, bdata);
 // msginfo->tag.data_ptr = (msginfo+1);
#ifdef RMO_DEBUG_
    a  = (char *) (msginfo + 1);
    a += msginfo->dscrlen;
    tmp = (double *) a;
    na = msginfo->datalen/sizeof(double);
    for(i=0; i<na; i++) {
       printf("%s [ds] tmp[%d]=%lf\n",Portals_ID(),i,tmp[i]);
    }
    fflush(stdout);
#endif
    if(armci_send_req_msg(proc,msginfo, bytes,tag))
       armci_die("armci_send_strided_req: failed",0);
    ARMCI_PR_DBG("exit",0);
}

/*\ client receives data from server
\*/
char *armci_rcv_data(int proc, request_header_t* msginfo,int rcvlen)
{
int datalen = msginfo->datalen;
char *buf;
    ARMCI_PR_DBG("enter",0);
    if(rcvlen)datalen=rcvlen;
    if(DEBUG_) {
        printf("%d:armci_rcv_data:  bytes= %d \n", armci_me, datalen);
        fflush(stdout);
    }

    if(datalen == 0) armci_die("armci_rcv_data: no data to receive",datalen);

    buf = armci_ReadFromDirect(proc, msginfo, datalen);

    if(DEBUG_){
        printf("%d:armci_rcv_data: got %d bytes \n",armci_me,datalen);
        fflush(stdout);
    }
    ARMCI_PR_DBG("exit",0);
    return(buf);
}

/*\ client receives vector data from server and unpacks to the right loc
\*/
void armci_rcv_vector_data(int proc, request_header_t* msginfo, armci_giov_t darr[], int len)
{
    ARMCI_PR_DBG("enter",0);
    char *buf = armci_rcv_data(proc, msginfo, 0);
    armci_vector_from_buf(darr, len, buf);
    ARMCI_PR_DBG("exit",0);
}

/*\ client receives strided data from server
\*/
void armci_rcv_strided_data(int proc, request_header_t* msginfo, int datalen,
                            void *ptr, int strides,int stride_arr[],int count[])
{
    char *databuf;
    ARMCI_PR_DBG("enter",0);
    if(DEBUG_){
        printf("%d: armci_rcv_strided_data: expecting datalen %d from %d\n",
                armci_me, datalen, proc); fflush(stdout);
    }
    databuf = armci_ReadFromDirect(proc,msginfo,0);
    armci_read_strided(ptr, strides, stride_arr, count, databuf);
    ARMCI_PR_DBG("exit",0);
}

void armci_rem_state(int clus)
{
int bufsize = sizeof(request_header_t)+sizeof(int);
int destproc = 0;
request_header_t *msginfo;
destproc = SERVER_NODE(clus);
msginfo = (request_header_t *)GET_SEND_BUFFER(bufsize,STATE,destproc);
int tag=0;

    ARMCI_PR_DBG("enter",0);
    msginfo->dscrlen = 0;
    msginfo->from  = armci_me;
    msginfo->to    = SERVER_NODE(clus);
    msginfo->operation = STATE;
    msginfo->bytes   =0;
    msginfo->datalen =sizeof(int);
 // msginfo->tag.data_ptr = (msginfo+1);

    if(DEBUG_){
       printf("%d(c):sending ACKreq to %d clus=%d\n",armci_me,msginfo->to,clus);
        fflush(stdout);
    }

    armci_send_req(armci_clus_info[clus].master, msginfo, bufsize,tag);
    armci_rcv_data(armci_clus_info[clus].master, msginfo,0);  /* receive  */
    FREE_SEND_BUFFER(msginfo);
    ARMCI_PR_DBG("exit",0);
}


/*\ get ACK from server
\*/
void armci_rem_ack(int clus)
{
int bufsize = sizeof(request_header_t)+sizeof(int);
int destproc = 0;
request_header_t *msginfo;
destproc = SERVER_NODE(clus);
msginfo = (request_header_t *) GET_SEND_BUFFER(bufsize,ACK,destproc);
int tag=0;

    ARMCI_PR_DBG("enter",0);
    msginfo->dscrlen = 0;
    msginfo->from  = armci_me;
    msginfo->to    = SERVER_NODE(clus);
    msginfo->operation = ACK;
    msginfo->bytes   =0;
    msginfo->datalen =sizeof(int);
 // msginfo->tag.data_ptr = (msginfo+1);

    if(DEBUG_){
       printf("%d(c):sending ACKreq to %d clus=%d\n",armci_me,msginfo->to,clus);
        fflush(stdout);
    }

    armci_send_req(armci_clus_info[clus].master, msginfo, bufsize,tag);
    armci_rcv_data(armci_clus_info[clus].master, msginfo,0);  /* receive ACK */
    FREE_SEND_BUFFER(msginfo);
    ARMCI_PR_DBG("exit",0);
}


/*\ request to QUIT sent by client
\*/
void armci_serv_quit()
{
int bufsize = sizeof(request_header_t)+sizeof(int);
int destproc;
request_header_t *msginfo;
destproc = SERVER_NODE(armci_clus_me);  
msginfo = (request_header_t*)GET_SEND_BUFFER(bufsize,QUIT,destproc);
int tag=0;

    ARMCI_PR_DBG("enter",0);
    if(DEBUG_){ printf("%d master: sending quit request to server\n",armci_me);
        fflush(stdout);
    }

    msginfo->dscrlen = 0;
    msginfo->from  = armci_me;
    msginfo->to    = SERVER_NODE(armci_clus_me);
    msginfo->operation = QUIT;
    if(ACK_QUIT)
       msginfo->bytes   = msginfo->datalen = sizeof(int); /* ACK */
    else
       msginfo->bytes   = msginfo->datalen = 0; /* no ACK */

    armci_send_req(armci_master, msginfo, bufsize,tag);

    if(ACK_QUIT){
       int stat;
       stat = *(int*)armci_rcv_data(armci_master,msginfo,0);  /* receive ACK */
       if(stat  != QUIT)
            armci_die("armci_serv_quit: wrong response from server", stat);
       FREE_SEND_BUFFER(msginfo);
    }
    ARMCI_PR_SDBG("exit",0);
}


/***************************** server side *********************************/

static void armci_check_req(request_header_t *msginfo, int buflen)
{

    ARMCI_PR_SDBG("enter",msginfo->operation);
    if((msginfo->to != armci_me && msginfo->to < armci_master) ||
       msginfo->to >= armci_master + armci_clus_info[armci_clus_me].nslave)
        /*armci_die("armci_check_req: invalid to", msginfo->to);*/
            printf("\n%d:got following to %d",armci_me,msginfo->to);
    if(msginfo->dscrlen < 0)
        armci_die("armci_check_req: dscrlen < 0", msginfo->dscrlen);
    if(msginfo->datalen < 0)
        armci_die("armci_check_req: datalen < 0", msginfo->datalen);
    if(msginfo->dscrlen > (int)msginfo->bytes)
        armci_die2("armci_check_req: dsclen > bytes", msginfo->dscrlen,
                   msginfo->bytes);
    ARMCI_PR_SDBG("exit",0);
}


/*\ server response - send data to client
\*/
void armci_send_data(request_header_t* msginfo, void *data)
{
    int to = msginfo->from;
    ARMCI_PR_SDBG("enter",0);
    armci_WriteToDirect(to, msginfo, data);
    ARMCI_PR_SDBG("exit",0);
}


/*\ server sends strided data back to client
\*/
void armci_send_strided_data(int proc,  request_header_t *msginfo,
                             char *bdata, void *ptr, int strides,
                             int stride_arr[], int count[])
{
    int i,na;
    double *a = NULL;
    int to = msginfo->from;
    ARMCI_PR_SDBG("enter",0);

    if(DEBUG_){ printf("%d(server): sending datalen = %d to %d %p\n",
                armci_me, msginfo->datalen, to,ptr); fflush(stdout); }
 
    /* for small contiguous blocks copy into a buffer before sending */
    armci_write_strided(ptr, strides, stride_arr, count, bdata);

#ifdef RMO_PORTALS_DEBUG_GET    
    a = (double *) bdata;
    na = msginfo->datalen/sizeof(double);
    for(i=0; i<na; i++) {
       printf("%s [ds] a[%d]= ",Portals_ID(),i); bit_print(a+i,sizeof(double));
       fflush(stdout);
    }
#endif

    /* write the message to the client */
    armci_WriteToDirect(to, msginfo, bdata);

    if(DEBUG_){
        printf("%d(serv):sent len=%d to %d\n",armci_me,msginfo->datalen,to);
        fflush(stdout);
    }
    ARMCI_PR_SDBG("exit",0);
}


/*\ server sends ACK to client
ptl_event_t *ev = (ptl_event_t *) msginfo->tag.user_ptr;

  ARMCI_PR_SDBG("enter",0);
    if(DEBUG_){
        printf("%d server: terminating request by %d\n",armci_me,msginfo->from);
            fflush(stdout);
              }

                portals_ds_send_ack(ev->initiator,ev->hdr_data);

\*/
void armci_server_ack(request_header_t* msginfo)
{
int *ack=(int*)(msginfo+1);

  ARMCI_PR_SDBG("enter",0);
  if(DEBUG_){
    printf("%d server: sending ACK to %d\n",armci_me,msginfo->from);
    fflush(stdout);
  }

  *ack = ACK;
  if(msginfo->datalen != sizeof(int))
    armci_die("armci_server_ack: bad datalen=",msginfo->datalen);
  armci_send_data(msginfo,ack);

  ARMCI_PR_SDBG("exit",0);
}



/*\ server action triggered by request to quit
\*/
void armci_server_goodbye(request_header_t* msginfo)
{
#ifdef LIBONESIDED

#else
  ptl_event_t *ev = (ptl_event_t *) msginfo->tag.user_ptr;

  ARMCI_PR_SDBG("enter",0);
  if(DEBUG_){
    printf("%d server: terminating request by %d\n",armci_me,msginfo->from);
    fflush(stdout);
  }

  portals_ds_send_ack(ev->initiator,ev->hdr_data);

#ifdef ARMCI_CHECK_STATE_
  for(int i=0;i<armci_nproc;i++)
    sarlist_print(i);
#endif

  ARMCI_PR_SDBG("exit",0);
  /* Finalizing data server process w.r.t. MPI is not portable
  */
  //_exit(0);
#endif
}


/*  main routine for data server process in a cluster environment
 *  the process is blocked until message arrives from
 *  the clients and services the requests
 */
void armci_data_server(void *mesg)
{
    /* message */
request_header_t non_volatile_mesg;
request_header_t *msginfo;
void *descr;
void *buffer;
int buflen;
int from, i;
static int mytag=1;
#   ifdef ARMCI_CHECK_STATE
    if(sarn_np==NULL)sarn_np=calloc(armci_nproc,sizeof(sarnode));
#   endif
    ARMCI_PR_SDBG("enter",0);

 // read header, descriptor, data, and buffer length 
    armci_rcv_req(mesg, &msginfo, &descr, &buffer, &buflen );

 // create a non volatile copy of the request_header_t
 // as soon as data is returned, the remote CP can overwrite mesg
 // we might not be completely done with this routine when that happens
 // this was THE cause for the original portals race condition
    memcpy(&non_volatile_mesg, msginfo, sizeof(request_header_t));
    msginfo = &non_volatile_mesg;

//  check what we got 
//  armci_check_req(msginfo,buflen);
    from = msginfo->from;

    if(DEBUG_){ 
       printf("%d(serv):got %d request from %d\n",armci_me,msginfo->operation,
               from);
       fflush(stdout);
    }

/*if(msginfo->operation==GET)fprintf(stderr,"GET request received with tag: %d\n",msginfo->tag);*/

    switch(msginfo->operation){
#    ifdef ARMCI_CHECK_STATE
      case STATE:
          printf("[ds %d]: operation=%d not supported yet\n",armci_me,msginfo->operation);
          abort();
          if(DEBUG_){printf("\n%d:state request\n",armci_me);fflush(stdout);}
          sarlist_print(msginfo->from);
          armci_WriteToDirect(msginfo->from, msginfo, (msginfo+1));
          break;
#    endif
      case QUIT:
          if(DEBUG_){
            printf("%d(serv):got QUIT request from %d\n",armci_me, from);
            fflush(stdout);
          }
          armci_server_goodbye(msginfo);
          break; /*pessimism?*/

      case ACK:
      //  printf("[ds %d]: operation=%d not supported yet\n",armci_me,msginfo->operation);
      //  abort();
          if(DEBUG_) {
              fprintf(stdout, "%d(server): got ACK request from %d\n",
                      armci_me, msginfo->from); fflush(stdout);
          }
          armci_server_ack(msginfo);
          break;

      case ATTACH: 
          printf("[ds %d]: operation=%d not supported yet\n",armci_me,msginfo->operation);
          abort();
       // if(DEBUG_){
       //    printf("%d(serv):got ATTACH request from%d\n",armci_me, from);
       //    fflush(stdout);
       // }
       // armci_server_ipc(msginfo, descr, buffer, buflen);
          break;
      case ARMCI_SWAP:
      case ARMCI_SWAP_LONG:
      case ARMCI_FETCH_AND_ADD:
      case ARMCI_FETCH_AND_ADD_LONG:
          armci_server_rmw(msginfo,descr,buffer);
          break;

      case LOCK:
          printf("[ds %d]: operation=%d not supported yet\n",armci_me,msginfo->operation);
          abort();
       // armci_server_lock(msginfo);
          break;

      case UNLOCK:
          printf("[ds %d]: operation=%d not supported yet\n",armci_me,msginfo->operation);
          abort();
       // armci_server_unlock(msginfo, descr);
       // msginfo->tag.ack=ARMCI_STAMP;
       // x_net_send_ack(msginfo,msginfo->from,msginfo->tag.ack_ptr,&msginfo->tag.ack);
          break;

      default:
          if(msginfo->format == VECTOR){
              armci_server_vector(msginfo, descr, buffer, buflen);      // point 1
              if(msginfo->operation==PUT || ARMCI_ACC(msginfo->operation)) {  // point 2
                 armci_server_send_ack(msginfo);
              }
           // the obove if clause and the similar cause below for a strided operation
           // was the reason for the race condition in the original portals code.
           // if the original request was a get, it could return it's data to the CP
           // once the the data is returned, the CP could fire off a new request which
           // would overwrite the 'now' volatile msginfo ... in which case, after returning
           // from armci_server_vector (having finished the get); the operation could now be
           // a put, in which case, it would repy back that it has also finished the put, 
           // with out actually doing it. Msginfo could be different at points 1 and 2 if 
           // at point 1 the operation was a get.
          }
          else if(msginfo->format == STRIDED){
           // if(msginfo->operation != PUT && msginfo->operation != GET && !ARMCI_ACC(msginfo->operation)) {
           //    printf("[ds %d]: operation=%d (format==STRIDED) not supported yet\n",armci_me,msginfo->operation);
           //    abort();
           // }
              armci_server(msginfo, descr, buffer, buflen);             // point 1
              if(msginfo->operation==PUT || ARMCI_ACC(msginfo->operation)){   // point 2
                 armci_server_send_ack(msginfo);
              }
          }
          else
              armci_die2("armci_data_serv: unknown format code",
                         msginfo->format, msginfo->from);
    }
    ARMCI_PR_SDBG("exit",0);
}



