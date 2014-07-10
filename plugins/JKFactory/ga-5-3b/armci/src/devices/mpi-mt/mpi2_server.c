#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <assert.h>
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

/* Abhinav Vishnu */
void armci_mpi2_server_init_twosided();

MPI_Comm MPI_COMM_SERVER2CLIENT = MPI_COMM_NULL;

static int armci_server_me=-1, armci_nserver=-1;
static int armci_client_first=-1, armci_nclients=-1;

extern Header *armci_get_shmem_ptr(int shmid, long shmoffset, size_t shmsize);


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

#if 1
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

void armci_mpi_strided_s2c(int op, void *ptr, int stride_levels, int stride_arr[],
                       int count[], int proc, MPI_Comm comm)
{
    int i, j;
    long idx;    /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int bvalue[MAX_STRIDE_LEVEL], bunit[MAX_STRIDE_LEVEL];
    MPI_Status status;
    
    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++)
        n1dim *= count[i];
  
    /* calculate the destination indices */
    bvalue[0] = 0; bvalue[1] = 0; bunit[0] = 1; bunit[1] = 1;
    for(i=2; i<=stride_levels; i++)
    {
        bvalue[i] = 0;
        bunit[i] = bunit[i-1] * count[i-1];
    }
     
    for(i=0; i<n1dim; i++)
    {
        idx = 0;
        for(j=1; j<=stride_levels; j++)
        {     
            idx += bvalue[j] * stride_arr[j-1];
            if((i+1) % bunit[j] == 0) bvalue[j]++;
            if(bvalue[j] > (count[j]-1)) bvalue[j] = 0;
        }          

        if(op == SEND)
        {
           MPI_Check(
              MPI_Send(((char*)ptr)+idx, count[0], MPI_BYTE, proc,
                       ARMCI_MPI_SERVER2CLIENT_TAG, comm)
              );
        }                      
        else /* ( op == RECV) */
        {
           MPI_Check(
              MPI_Recv(((char*)ptr)+idx, count[0], MPI_BYTE, proc,
                       ARMCI_MPI_CLIENT2SERVER_TAG, comm, &status)
              );
        }
    }
}   

void check_comm()
{
    int result;
    assert(SERVER_CONTEXT);
    assert(ARMCI_COMM_WORLD != MPI_COMM_NULL);
}

/**************************************************************************
 * Platform specific server code as required by the ARMCI s/w layer. (BEGIN)
 */

/* establish connections with client (i.e compute) processes */
void armci_server_initial_connection()
{
    armci_mpi2_server_init_twosided();
    armci_mpi2_server_debug(0, "armci_server_initial_connection\n");   
}


/* close all open connections, called before terminating/aborting */
void armci_transport_cleanup()
{
#if 0
    /* armci_transport_cleanup is called by all procs (clients and servers).
       Therefore, only in server case we need to finalize MPI before exit. */
    if(ARMCI_COMM_WORLD != MPI_COMM_NULL) 
    {
       armci_mpi2_server_debug(0, "Calling MPI_Finalize\n");
       MPI_Finalize();
       exit(EXIT_SUCCESS); /* server termination */
    }
#endif
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

    check_comm();
    {
       armci_mpi_strided_s2c(RECV, ptr, stride_levels, stride_arr, count, from,
                         ARMCI_COMM_WORLD);
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

    check_comm();
    MPI_Status status;
    msginfo = (request_header_t*) MessageRcvBuffer;
    p = * (int *) mesg;
    
    MPI_Check(
       MPI_Recv(MessageRcvBuffer, MSG_BUFLEN, MPI_BYTE, p, ARMCI_MPI_CLIENT2SERVER_TAG,
                ARMCI_COMM_WORLD, &status)
       );
    
    * (void **) phdr = msginfo;    

    if( !(p >= 0 && p < armci_nproc) )
       armci_die("armci_rcv_req: request from invalid client", p);
    
    armci_mpi2_server_debug(armci_server_me,
                            "armci_rcv_req: op=%d mesg=%p, phdr=%p "
                            "pdata=%p, buflen=%p, p=%d\n", msginfo->operation,
                            mesg, phdr, pdata, buflen, p, MSG_BUFLEN);
    
#ifdef MPI_SPAWN_ZEROCOPY
    assert(0);
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

    check_comm();
    MPI_Check(
            MPI_Send(data, msginfo->datalen, MPI_BYTE, to,
                ARMCI_MPI_SERVER2CLIENT_TAG, ARMCI_COMM_WORLD)
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
    check_comm();

    {
       armci_mpi_strided_s2c(SEND, ptr, strides, stride_arr, count, to,
                         ARMCI_COMM_WORLD);
    }
    
}


void armci_call_data_server() 
{
    int p=-1;
    MPI_Status status;
    
    armci_mpi2_server_debug(0, "armci_call_data_server(): Server main loop\n");
   
    int result;

    check_comm();
    /* server main loop; wait for and service requests until QUIT requested */
    int flag = 0;
    for(;;)
    {       
#if 1
       while (!flag) { 
            MPI_Check(
                    MPI_Iprobe(MPI_ANY_SOURCE, ARMCI_MPI_CLIENT2SERVER_TAG,
                        ARMCI_COMM_WORLD,
                        &flag, &status)
                    );
       }
#else
            MPI_Check(
                    MPI_Probe(MPI_ANY_SOURCE, ARMCI_MPI_CLIENT2SERVER_TAG,
                        ARMCI_COMM_WORLD,
                        &status)
                    );
#endif
       flag = 0;
       p = status.MPI_SOURCE;
       assert((p>= 0) && (p < armci_nproc));
       armci_mpi2_server_debug(armci_server_me,
                               "Processing message from client %d\n", p);
       
       armci_data_server(&p);  
    }

}
/**
 * Platform specific server code ENDs here.
 **************************************************************************/

static void emulate_armci_init_clusinfo() 
{
    assert(0);
}

static void emulate_armci_allocate_locks(long *shm_info)
{
    assert(0);
}

/* Abhinav Vishnu */
void armci_mpi2_server_init_twosided()
{
    int namelen, version, subversion;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    long shm_info[3];
    MPI_Status status;

    assert(ARMCI_COMM_WORLD != MPI_COMM_NULL);
    MPI_Check(MPI_Comm_rank(ARMCI_COMM_WORLD, &armci_server_me));
    armci_nserver = armci_nclus;
    MPI_Check(MPI_Get_processor_name(processor_name, &namelen));
    MPI_Check(MPI_Get_version(&version, &subversion));

}

void armci_mpi2_server_init() 
{
    assert(0);
}

void armci_mpi2_server() 
{
    assert(0);
}

void armci_comm_dup_server(MPI_Comm * comm)
{
    /* Do nothing */
}

