#if HAVE_CONFIG_H
#   include "config.h"
#endif

/**
 * mpi2_server.c: MPI_SPAWN Server Code
 * Manojkumar Krishnan
 */

#if HAVE_STDARG_H
#   include <stdarg.h>
#endif
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#include "mpi.h"

#include "armcip.h"
#include "mpi2.h"
#include "kr_malloc.h"
#include "locks.h"

/* Inter-communicators for communicating with clients */
MPI_Comm MPI_COMM_SERVER2CLIENT=MPI_COMM_NULL;

static int armci_server_me=-1, armci_nserver=-1;
static int armci_client_first=-1, armci_nclients=-1;

extern Header *armci_get_shmem_ptr(int shmid, long shmoffset, size_t shmsize);

/* ====================== MUltiple Buffers Code ======================== */
#ifdef MULTIPLE_BUFS

#define MPI2_MAX_BUFS 10

typedef struct req_waitlist 
{
   int reqid;
   int p;   
   int tag;
   struct req_waitlist *next;   
}
req_waitlist_t;

static double      _mpi2_rcv_buf[MPI2_MAX_BUFS][MSG_BUFLEN_DBL];
static MPI_Request _mpi_request[MPI2_MAX_BUFS];
static int *_next_tag=NULL;
static int _reqid_ready=0;
static req_waitlist_t *_req_waitlist_head = NULL;
static req_waitlist_t *_req_waitlist_tail = NULL;

/* increment the tag value to maintain flowcontrol */
#define INCR_TAG(p) \
{                                                         \
    if(++_next_tag[p] > ARMCI_MPI_SPAWN_TAG_END)          \
            _next_tag[p] = ARMCI_MPI_SPAWN_TAG_BEGIN;     \
}


static int wlist_get_req(int *p, int *tag, int *reqid) {

    req_waitlist_t *itr, *prev=NULL;
    
    for(itr=_req_waitlist_head; itr != NULL; itr=itr->next)
    {
       if(itr->tag == _next_tag[itr->p]) 
       {
          /* mark the request id as ready to be processed, and update the
           * tag order and waitlist */
          *p     = itr->p;
          *tag   = itr->tag;
          *reqid = itr->reqid;
          INCR_TAG(*p);
          
          /* remove this request from the waiting list */
          if(itr==_req_waitlist_head) _req_waitlist_head = itr->next;
          if(itr==_req_waitlist_tail) _req_waitlist_tail = prev;
          if(prev != NULL) prev->next = itr->next;
          free(itr);
          
          return 1;
       }       
       prev = itr;
    }
    
    return 0;
}

static void wlist_add_req(int reqid, int p, int tag) 
{
    req_waitlist_t *node = (req_waitlist_t *) malloc(sizeof(req_waitlist_t));
    
    node->reqid = reqid;
    node->p     = p;
    node->tag   = tag;
    node->next  = NULL;
    
    if(_req_waitlist_head == NULL)
    {
       _req_waitlist_head = node;
       _req_waitlist_tail = node;
       return;
    }

    /* append the new request at the end of the list */
    _req_waitlist_tail->next = node;
    _req_waitlist_tail = node;
}
#endif

/* ==================== END: Multiple Buffers Code ====================== */

#if MPI_SPAWN_DEBUG
void armci_mpi2_server_debug(int rank, const char *format, ...) 
{
    va_list arg;
    
    if(rank == armci_server_me) {
       va_start(arg, format);
       printf("**** Server %d: ", armci_server_me);
       vprintf(format, arg);
       va_end(arg);
       fflush(stdout);
    }
}
#else
# define armci_mpi2_server_debug(x, ...)
#endif

#if MPI_SPAWN_DEBUG
static inline int MPI_Check (int status)
{
    if(status != MPI_SUCCESS) 
    {
       armci_mpi2_server_debug(armci_me, "MPI Check failed.\n");
       armci_die("MPI_Check failed.", 0);
    }
}
#else
# define MPI_Check(x) x
#endif


/**************************************************************************
 * Platform specific server code as required by the ARMCI s/w layer. (BEGIN)
 */

/* establish connections with client (i.e compute) processes */
void armci_server_initial_connection()
{
    armci_mpi2_server_debug(0, "armci_server_initial_connection\n");   
}


/* close all open connections, called before terminating/aborting */
void armci_transport_cleanup()
{
    /* armci_transport_cleanup is called by all procs (clients and servers).
       Therefore, only in server case we need to finalize MPI before exit. */
    if(MPI_COMM_SERVER2CLIENT != MPI_COMM_NULL) 
    {
       armci_mpi2_server_debug(0, "Calling MPI_Finalize\n");
       MPI_Finalize();
       exit(EXIT_SUCCESS); /* server termination */
    }
}

static void armci_mpi_rcv_strided_data(request_header_t *msginfo,
                                       void *vdscr, int from) 
{    
    int bytes;
    void *ptr;
    char *dscr;
    int stride_levels, *stride_arr, *count;
    
    bytes = msginfo->dscrlen;
    dscr  = (char*)(msginfo + 1);
    *(void**)vdscr = (void *)dscr;
    
    ptr = *(void**)dscr;           dscr += sizeof(void*);
    stride_levels = *(int*)dscr;   dscr += sizeof(int);
    stride_arr = (int*)dscr;       dscr += stride_levels*sizeof(int);
    count = (int*)dscr;            dscr += (stride_levels+1)*sizeof(int);

#ifdef MPI_USER_DEF_DATATYPE
    if(stride_levels>0) 
    {
       armci_mpi_strided2(RECV, ptr, stride_levels, stride_arr, count, from,
                          MPI_COMM_SERVER2CLIENT);
    }
    else
#endif
    {
       armci_mpi_strided(RECV, ptr, stride_levels, stride_arr, count, from,
                         MPI_COMM_SERVER2CLIENT);
    }
}

static void armci_mpi_rcv_vector_data(request_header_t *msginfo,
                                      void *vdscr, int proc) 
{
       armci_die("armci_mpi_rcv_vector_data(): Not yet implemented!", 0);    
}

/* server receives request */
void armci_rcv_req (void *mesg, void *phdr, void *pdescr,
                    void *pdata, int *buflen)
{
    request_header_t *msginfo = NULL;
    int hdrlen = sizeof(request_header_t);
    int p=-1;
    int bytes;

#if !defined(MULTIPLE_BUFS)
    MPI_Status status;
    msginfo = (request_header_t*) MessageRcvBuffer;
    p = * (int *) mesg;
    
    MPI_Check(
       MPI_Recv(MessageRcvBuffer, MSG_BUFLEN, MPI_BYTE, p, ARMCI_MPI_SPAWN_TAG,
                MPI_COMM_SERVER2CLIENT, &status)
       );
#else
    int reqid = _reqid_ready;;/*get request id that is ready to be processed */
    
    msginfo = (request_header_t*) _mpi2_rcv_buf[reqid];
    p = * (int *) mesg;
    if(p != msginfo->from)
       armci_die("armci_rcv_req: invalid client", p);
#endif
    
    * (void **) phdr = msginfo;    

    if( !(p >= 0 && p < armci_nproc) )
       armci_die("armci_rcv_req: request from invalid client", p);
    
    armci_mpi2_server_debug(armci_server_me,
                            "armci_rcv_req: op=%d mesg=%p, phdr=%p "
                            "pdata=%p, buflen=%p, p=%d\n", msginfo->operation,
                            mesg, phdr, pdata, buflen, p, MSG_BUFLEN);
    
#ifdef MPI_SPAWN_ZEROCOPY
    if(msginfo->operation==PUT && msginfo->datalen==0)
    {
       if(msginfo->format==STRIDED) 
       {
          armci_mpi_rcv_strided_data(msginfo, pdescr, p);
       }
       if(msginfo->format==VECTOR)
       {
          armci_mpi_rcv_vector_data(msginfo, pdescr, p);
       }
       return;
    }
#endif
    
    *buflen = MSG_BUFLEN - hdrlen;
    if (msginfo->operation == GET)
    {
       bytes = msginfo->dscrlen;
    }
    else
    {
       bytes = msginfo->bytes;
       if (bytes > *buflen)
          armci_die2("armci_rcv_req: message overflowing rcv buf",
                     msginfo->bytes, *buflen);
    }

#if MPI_SPAWN_DEBUG && !defined(MPI_SPAWN_ZEROCOPY) && 0
    {
       int count;
       MPI_Get_count(&status, MPI_BYTE, &count);
       if (count != (bytes + hdrlen))
       {
          armci_mpi2_server_debug(armci_server_me, "armci_rcv_req: "
                                  "got %d bytes, expected %d bytes\n",
                                  count, bytes + hdrlen);
          printf("%d: armci_rcv_req: got %d bytes, expected %d bytes (%d)\n",
                 armci_me, count, bytes + hdrlen,  msginfo->datalen);
          armci_die("armci_rcv_req: count check failed.\n", 0);
       }
    }
#endif
    
    if (msginfo->bytes)
    {
       * (void **) pdescr = msginfo + 1;
       * (void **) pdata  = msginfo->dscrlen + (char *) (msginfo+1); 
       *buflen -= msginfo->dscrlen; 
       
       if (msginfo->operation != GET && msginfo->datalen)
       {
          *buflen -= msginfo->datalen;
       }
    }
    else
    {
       * (void**) pdata  = msginfo + 1;
       * (void**) pdescr = NULL;
    }
    
    if (msginfo->datalen > 0 && msginfo->operation != GET)
    {
       if (msginfo->datalen > (MSG_BUFLEN - hdrlen - msginfo->dscrlen)) {
          armci_die2("armci_rcv_req:data overflowing buffer",
                     msginfo->dscrlen, msginfo->datalen);
       }
       *buflen -= msginfo->datalen;
    }
}

/* server sends data back to client */
void armci_WriteToDirect (int to, request_header_t *msginfo, void *data)
{
    armci_mpi2_server_debug(armci_server_me, "armci_WriteToDirect: "
                            "to=%d, msginfo=%p, data=%p, bytes=%d\n",
                            to, msginfo, data, msginfo->datalen);

    if( !(to >= 0 && to < armci_nproc) )
       armci_die("armci_WriteToDirect: send request to invalid client", to);
    
    MPI_Check(
       MPI_Send(data, msginfo->datalen, MPI_BYTE, to,
                ARMCI_MPI_SPAWN_TAG, MPI_COMM_SERVER2CLIENT)
       );
}

/*\ server sends strided data back to client 
\*/
void armci_WriteStridedToDirect(int to, request_header_t* msginfo,
                                void *ptr, int strides, int stride_arr[],
                                int count[])
{
    armci_mpi2_server_debug(armci_server_me, "armci_WriteStridedToDirect: "
                            "to=%d, stride_levels=%d, bytes=%d\n", to, strides,
                            msginfo->datalen);
 
#ifdef MPI_USER_DEF_DATATYPE
    if(strides>0) 
    {
       armci_mpi_strided2(SEND, ptr, strides, stride_arr, count, to,
                          MPI_COMM_SERVER2CLIENT);
    }
    else
#endif
    {
       armci_mpi_strided(SEND, ptr, strides, stride_arr, count, to,
                         MPI_COMM_SERVER2CLIENT);
    }
    
}


void armci_call_data_server() 
{
    int p=-1;
    MPI_Status status;
    
    armci_mpi2_server_debug(0, "armci_call_data_server(): Server main loop\n");
    
#if !defined(MULTIPLE_BUFS)
    /* server main loop; wait for and service requests until QUIT requested */
    for(;;)
    {       
       MPI_Check(
          MPI_Probe(MPI_ANY_SOURCE, ARMCI_MPI_SPAWN_TAG,MPI_COMM_SERVER2CLIENT,
                    &status)
          );

       p = status.MPI_SOURCE;
       armci_mpi2_server_debug(armci_server_me,
                               "Processing message from client %d\n", p);
       
       armci_data_server(&p);  
    }
#else

    int i, tag, reqid, do_waitlist=0;

    /* server multiple bufs setup */
    _req_waitlist_head = NULL;
    _req_waitlist_tail = NULL;
    /* Initialize "next tag" array, which manages flow control */
    if( (_next_tag = (int*) malloc(armci_nproc*sizeof(int)) ) == NULL)
       armci_die("mpi2_server: _next_tag malloc failed", 0);
    for(i=0; i<armci_nproc; i++) _next_tag[i] = ARMCI_MPI_SPAWN_TAG_BEGIN;
    
    
    /* server posts multiple receive buffers in advance */
    for(i=0; i<MPI2_MAX_BUFS; i++) 
    {
       MPI_Check(
          MPI_Irecv(_mpi2_rcv_buf[i], MSG_BUFLEN, MPI_BYTE, MPI_ANY_SOURCE,
                    ARMCI_MPI_SPAWN_TAG, MPI_COMM_SERVER2CLIENT,
                    &_mpi_request[i])
          );
    }
    
    for(;;)
    {
       /* process wait-listed requests, if any */
       do_waitlist = 0;       
       if(_req_waitlist_head != NULL) 
       {
          do_waitlist = wlist_get_req(&p, &tag, &reqid);
       }
       
       if(!do_waitlist) 
       {  
          /* process the first completed incoming request */
          MPI_Check(
             MPI_Waitany(MPI2_MAX_BUFS, _mpi_request, &reqid, &status)
             );
          p   = status.MPI_SOURCE;
          /* tag = status.MPI_TAG; */
          tag = ((request_header_t*) _mpi2_rcv_buf[reqid])->tag;
          
          /* check if it is in or out of order request */
          if(tag == _next_tag[p]) { INCR_TAG(p); }
          else 
          {
             /* out of order req - enforce ordering by waitlisting this req */
             wlist_add_req(reqid, p, tag);
             continue;
          }
       }
       
       /* mark the request id that is ready to processed */
       _reqid_ready = reqid;
       
       /* server process the incoming (or waitlisted) request */
       armci_data_server(&p);

       /* After completing the request (which also frees a buffer), server
        * posts a receive using this buffer */
       MPI_Check(
          MPI_Irecv(_mpi2_rcv_buf[reqid], MSG_BUFLEN, MPI_BYTE, MPI_ANY_SOURCE,
                    ARMCI_MPI_SPAWN_TAG, MPI_COMM_SERVER2CLIENT,
                    &_mpi_request[reqid])
          );
    }
#endif

}
/**
 * Platform specific server code ENDs here.
 **************************************************************************/

static void emulate_armci_init_clusinfo() 
{
    int psize;
    
    MPI_Comm_remote_size(MPI_COMM_SERVER2CLIENT, &psize);
    
    /* server id (i.e. server's armci_me) is derived from node master's id.
       Similar to armci_create_server_process() to set SERVER_CONTEXT */
    armci_me        = SOFFSET - armci_client_first; 
    armci_nproc     = psize;
    armci_usr_tid   = THREAD_ID_SELF(); /*remember the main user thread id */
    armci_master    = armci_client_first;
    
    /* ***** emulate armci_init_clusinfo() ***** */
    armci_clus_me    = armci_server_me;
    armci_nclus      = armci_nserver;
    armci_clus_first = armci_clus_info[armci_clus_me].master;
    armci_clus_last  = (armci_clus_first +
                        armci_clus_info[armci_clus_me].nslave - 1);
    
    if(armci_clus_first != armci_client_first ||
       armci_nclients   != armci_clus_info[armci_clus_me].nslave) 
    {
       armci_mpi2_server_debug(armci_server_me,
                               "armci_clus_first=%d, armci_clus_last=%d\n",
                               armci_clus_first, armci_clus_last);
       armci_die("mpi2_server: armci_clus_info is incorrect.", 0);
    }
}

static void emulate_armci_allocate_locks(long *shm_info)
{
    int shmid      = (int) shm_info[0];
    long shmoffset = shm_info[1];
    size_t shmsize = (size_t) shm_info[2];
    
    _armci_int_mutexes = (PAD_LOCK_T*) armci_get_shmem_ptr(shmid, shmoffset,
                                                           shmsize);
}


void armci_mpi2_server_init() 
{

    int namelen, version, subversion;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    long shm_info[3];
    MPI_Status status;
    
    MPI_Comm_rank(ARMCI_COMM_WORLD, &armci_server_me);
    MPI_Comm_size(ARMCI_COMM_WORLD, &armci_nserver);
    MPI_Get_processor_name(processor_name, &namelen);
    MPI_Get_version(&version, &subversion);
    
    armci_mpi2_server_debug(armci_server_me,
                            "I'm %d of %d SERVERS running on %s (MPI %d.%d)\n",
                            armci_server_me, armci_nserver, processor_name,
                            version, subversion);
    
    /* get parent's groupinfo */
    MPI_Comm_get_parent(&MPI_COMM_SERVER2CLIENT);
    if (MPI_COMM_SERVER2CLIENT == MPI_COMM_NULL) 
    {
       armci_die("mpi2_server: Invalid spawn. No parent process found.\n",0);
    }

    /* receive my clients info */
    {
       int msg[3];
       MPI_Recv(msg, 3, MPI_INT, MPI_ANY_SOURCE, ARMCI_MPI_SPAWN_INIT_TAG,
                MPI_COMM_SERVER2CLIENT, &status);
       if(msg[0]-ARMCI_MPI_SPAWN_INIT_TAG != armci_server_me) 
       {
          armci_die("mpi2_server: Recv failed", msg[0]);
       }

       armci_client_first = msg[1];
       armci_nclients     = msg[2];
       armci_mpi2_server_debug(armci_server_me,
                               "My clients are [%d-%d]\n", armci_client_first,
                               armci_client_first+armci_nclients-1);
    }

    /**********************************************************************
     * Emulate PARMCI_Init().
     * Spawned Data server processes emulate PARMCI_Init() to complete the
     * ARMCI Initalization process similar to clients
     */
    
    armci_clus_info =(armci_clus_t*)malloc(armci_nserver*sizeof(armci_clus_t));
    if(armci_clus_info == NULL)
    {
       armci_die("mpi2_server: armci_clus_info malloc failed", 0);
    }
    
    /* receive and emulate clus info, lock info */
    MPI_Recv(armci_clus_info, armci_nserver*sizeof(armci_clus_t), MPI_BYTE,
             armci_client_first, ARMCI_MPI_SPAWN_INIT_TAG,
             MPI_COMM_SERVER2CLIENT, &status);
    
    MPI_Recv(shm_info, 3, MPI_LONG, armci_client_first,
             ARMCI_MPI_SPAWN_INIT_TAG, MPI_COMM_SERVER2CLIENT, &status);
    
    
    /* server setup clusinfo&locks, exactly as this node's armci master */
    emulate_armci_init_clusinfo();  /* armci_init_clusinfo()  in PARMCI_Init */
    emulate_armci_allocate_locks(shm_info); /* armci_allocate_locks() */

    /* Fence data structures should be initialized by server too. see
     * armci_generic_rmw (called by armci_server_rmw) */
    armci_init_fence();
    
    /**
     * End of PARMCI_Init() emulation.
     * *******************************************************************/

    MPI_Barrier(ARMCI_COMM_WORLD);
}

void armci_mpi2_server() 
{

    armci_mpi2_server_init();
    
    armci_server_code(NULL);

    /* Should never come here */
    armci_die("mpi2_server: server died in an unexpected manner", 0);
}

