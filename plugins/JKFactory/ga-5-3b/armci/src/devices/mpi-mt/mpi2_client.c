#if HAVE_CONFIG_H
#   include "config.h"
#endif

/**
 * MPI_SPAWN: ARMCI on top of MPI Multithreaded
 * Abhinav Vishnu 
 */

#if HAVE_STDARG_H
#   include <stdarg.h>
#endif
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#include "mpi.h"

#include "mpi2.h"
#include "armcip.h"
#include "request.h"
#include "shmem.h"
#include "locks.h"

#include <assert.h>

#define ARMCI_ROOT 0 /* root process */

/* Inter-communicators for communicating between clients and data servers */
MPI_Comm MPI_COMM_CLIENT2SERVER=MPI_COMM_NULL;

static int armci_nserver=-1;
static int *_armci_mpi_tag=NULL;

extern char ***_armci_argv;
extern int armci_get_shmem_info(char *addrp,  int* shmid, long *shmoffset,
                                size_t *shmsize);

#if MPI_SPAWN_DEBUG
void armci_mpi2_debug(int rank, const char *format, ...) 
{
    va_list arg;
    if(rank == armci_me) 
    {
       va_start(arg, format);
       printf("%d: ", rank);
       vprintf(format, arg);
       va_end(arg);
       fflush(stdout);
    }
}
#else
#define armci_mpi2_debug(x, ...)
#endif

#if MPI_SPAWN_DEBUG
static inline int MPI_Check (int status)
{
    if(status != MPI_SUCCESS) 
    {
       armci_mpi2_debug(armci_me, "MPI Check failed.\n");
       armci_die("MPI_Check failed.", 0);
    }
}
#else
# define MPI_Check(x) x
#endif


/**************************************************************************
 * Platform specific server code as required by the ARMCI s/w layer. (BEGIN)
 */

/* Create connections between clients and servers */
void armci_init_connections()
{
    armci_mpi2_debug(0, "armci_init_connections\n");
    _armci_buf_init();    /* CHECK: Is this correct ? */
    MPI_Check(MPI_Barrier(ARMCI_COMM_WORLD));
    /* Abhinav Vishnu */
    armci_create_server_MPIprocess();

    armci_mpi2_debug(0, "armci_init_connections completed\n");
}

void armci_wait_for_server()
{
    armci_mpi2_debug(0, "armci_wait_for_server: wait for server to quit\n");
    if (armci_me == armci_master)
    {
       armci_serv_quit();  
    }
}

void armci_client_connect_to_servers()
{

    /* Abhinav Vishnu */
}

/* NOTE: armci_mpi_strided and armci_mpi_strided2 are the only 2 functions
 * that are common to client and server part  */
void armci_mpi_strided_c2s(int op, void *ptr, int stride_levels, int stride_arr[], 
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
                       ARMCI_MPI_CLIENT2SERVER_TAG, comm)
              );
        }
        else /* ( op == RECV) */
        {
           MPI_Check(
              MPI_Recv(((char*)ptr)+idx, count[0], MPI_BYTE, proc,
                       ARMCI_MPI_SERVER2CLIENT_TAG, comm, &status)
              );
        }
    }
}

/* This is the only function that is common to client and server part */
void armci_mpi_strided2(int op, void *ptr, int stride_levels, int stride_arr[],
                        int count[], int proc, MPI_Comm comm)
{
    /* Not supported */
    assert(0);
}

/*\ client sends request message to server
\*/
int armci_send_req_msg (int proc, void *buf, int bytes)
{
  int clus_id = armci_clus_id(proc);
  int server ;

  /* Abhinav Vishnu */
  server = armci_clus_info[clus_id].master;

  armci_mpi2_debug(armci_me, "armci_send_req_msg(): proc=%d, server=%d, "
                   "buf=%p, bytes=%d\n", proc, server, buf, bytes);

  MPI_Check(
     MPI_Send(buf, bytes, MPI_BYTE, server, ARMCI_MPI_CLIENT2SERVER_TAG,
              ARMCI_COMM_WORLD)
     );
  armci_mpi2_debug(armci_me, "armci_send_req_msg(): send msg to server(%d), to"
                   "fwd to client %d\n", server, proc);

  return 0;
}

/*\ client sends strided data + request to server
\*/
int armci_send_req_msg_strided(int proc, request_header_t *msginfo,char *ptr,
                               int strides, int stride_arr[], int count[])
{
    int server;
    int clus_id = armci_clus_id(proc);
    int bytes;

    /* Abhinav Vishnu */
    server = armci_clus_info[clus_id].master;

    armci_mpi2_debug(armci_me, "armci_send_req_msg_strided: proc=%d server=%d "
                     "bytes=%d (op=%d)\n", proc, server, msginfo->datalen,
                     msginfo->operation);


    /* we write header + descriptor of strided data  */
    bytes = sizeof(request_header_t) + msginfo->dscrlen;
    armci_send_req_msg(proc, msginfo, bytes);
    
    {
       /* for larger blocks write directly thus avoiding memcopy */
       armci_mpi_strided_c2s(SEND, ptr, strides, stride_arr, count, server,
                         ARMCI_COMM_WORLD);
    }
       

    armci_mpi2_debug(armci_me, "armci_send_req_msg_strided(): send msg to "
                     "server(%d), to fwd to client %d\n", server, proc);

    return 0;
}

/*\ client receives data from server
\*/
char *armci_ReadFromDirect (int proc, request_header_t *msginfo, int len)
{

    int server;
    int clus_id = armci_clus_id(proc);
    MPI_Status status;

    server = armci_clus_info[clus_id].master;

    armci_mpi2_debug(armci_me, "armci_ReadFromDirect: proc=%d, server=%d, "
                     "msginfo=%p, bytes=%d (op=%d)\n", proc, server, msginfo,
                     len, msginfo->operation);
    MPI_Check(
       MPI_Recv(msginfo + 1, len, MPI_BYTE, server, ARMCI_MPI_SERVER2CLIENT_TAG,
                ARMCI_COMM_WORLD, &status)
       );

    
    armci_mpi2_debug(armci_me, "recv msg from server(%d), fwd by client %d\n",
                     server, proc);

    {
        int count;
        MPI_Get_count(&status, MPI_BYTE, &count);
        if (count != len) 
        {
            armci_mpi2_debug(armci_me, "armci_ReadFromDirect: got %d bytes, "
                    "expected %d bytes\n", count, len);
            armci_die("armci_ReadFromDirect: MPI_Recv failed.", count);
        }
    }
    
 return (char *) (msginfo+1);
}

/*\ client receives strided data from server
\*/
void armci_ReadStridedFromDirect(int proc, request_header_t* msginfo,
                                 void *ptr, int strides, int stride_arr[],
                                 int count[])
{

    int server;
    int clus_id = armci_clus_id(proc);

    /* Abhinav Vishnu */
    server = armci_clus_info[clus_id].master;
    
    armci_mpi2_debug(armci_me, "armci_ReadStridedFromDirect: proc=%d "
                     "stride_levels=%d, server=%d bytes=%d (op=%d)\n",
                     proc, strides, server, msginfo->datalen,
                     msginfo->operation);

    {
       armci_mpi_strided_c2s(RECV, ptr, strides, stride_arr, count, server,
                         ARMCI_COMM_WORLD);
    }
}



/**
 * Platform specific server code ENDs here. (END)
 **************************************************************************/

static void armci_gather_hostnames(char **hostname_arr) 
{
    /* Code must not flow through here */
    assert(0);
}



/**
 * Create server processes. This is called in armci_start_server.
 * Must be called after armci_init_clusinfo().
 */
void armci_create_server_MPIprocess ()
{
    int rank, size, flag, i;

    MPI_Check(MPI_Initialized(&flag));
    if (flag == 0)
       armci_die("ARMCI error: MPI_Init must be called before PARMCI_Init()",0);
    
    MPI_Check(MPI_Comm_rank(ARMCI_COMM_WORLD, &rank));
    MPI_Check(MPI_Comm_size(ARMCI_COMM_WORLD, &size));

    armci_nserver = armci_nclus;
    
    /* makesure all processes sync here. CHECK: does it ensure global sync ? */
    MPI_Check(MPI_Barrier(ARMCI_COMM_WORLD));
 
    armci_mpi2_debug(0, "armci_create_server_MPIprocess: Servers spawned!\n");
}

  
