#if HAVE_CONFIG_H
#   include "config.h"
#endif

/**
 * MPI_SPAWN: ARMCI on top of MPI.
 * Manojkumar Krishnan
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
    armci_mpi2_debug(0, "armci_client_connect_to_servers\n");   
}

/* NOTE: armci_mpi_strided and armci_mpi_strided2 are the only 2 functions
 * that are common to client and server part  */
void armci_mpi_strided(int op, void *ptr, int stride_levels, int stride_arr[], 
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
                       ARMCI_MPI_SPAWN_DATA_TAG, comm)
              );
        }
        else /* ( op == RECV) */
        {
           MPI_Check(
              MPI_Recv(((char*)ptr)+idx, count[0], MPI_BYTE, proc,
                       ARMCI_MPI_SPAWN_DATA_TAG, comm, &status)
              );
        }
    }
}

/* This is the only function that is common to client and server part */
void armci_mpi_strided2(int op, void *ptr, int stride_levels, int stride_arr[],
                        int count[], int proc, MPI_Comm comm)
{
    int i, stride=1;
    MPI_Status status;
    MPI_Datatype type[MAX_STRIDE_LEVEL];

    if(stride_levels == 0) 
    {
        armci_mpi_strided(op, ptr, stride_levels, stride_arr, count, proc,
                          comm);
        return;
    }

    /* convert stided data desciption to MPI type */
    type[0] = MPI_BYTE;
    for(i=1; i<=stride_levels; i++) 
    {
       stride *= stride_arr[i-1];
       MPI_Check( MPI_Type_hvector(count[i], count[i-1], stride,
                                  type[i-1], &type[i]) );
    }
    MPI_Check( MPI_Type_commit(&type[stride_levels]) );
    
    if(op == SEND) 
    {
       MPI_Check( MPI_Send(ptr, 1, type[stride_levels], proc,
                           ARMCI_MPI_SPAWN_VDATA_TAG, comm) );
    }
    else /* ( op == RECV) */
    {
       MPI_Check( MPI_Recv(ptr, 1, type[stride_levels], proc,
                           ARMCI_MPI_SPAWN_VDATA_TAG, comm, &status) );
    }
}

/*\ client sends request message to server
\*/
int armci_send_req_msg (int proc, void *buf, int bytes)
{
  int server = armci_clus_id(proc);

  armci_mpi2_debug(armci_me, "armci_send_req_msg(): proc=%d, server=%d, "
                   "buf=%p, bytes=%d\n", proc, server, buf, bytes);
  
  if( !(server >= 0 && server < armci_nserver) )
     armci_die("armci_send_req_msg: Invalid server.", 0);

#ifdef MULTIPLE_BUFS
  /**
   * Sequentially ordered tags to ensure flow control at the server side.
   * For example, a put followed by get from a client should be processed in
   * ORDER at the server side. If we don't have the flow control, the server
   * might process the get request first instead of put (and thus violating
   * ARMCI's ordering semantics.
   */
  ((request_header_t*)buf)->tag = _armci_mpi_tag[server];
  MPI_Check(
     MPI_Send(buf, bytes, MPI_BYTE, server, ARMCI_MPI_SPAWN_TAG,
              MPI_COMM_CLIENT2SERVER)
     );

  _armci_mpi_tag[server]++;
  if(_armci_mpi_tag[server] > ARMCI_MPI_SPAWN_TAG_END) 
     _armci_mpi_tag[server] = ARMCI_MPI_SPAWN_TAG_BEGIN;
  
#else
  MPI_Check(
     MPI_Send(buf, bytes, MPI_BYTE, server, ARMCI_MPI_SPAWN_TAG,
              MPI_COMM_CLIENT2SERVER)
     );
#endif
  armci_mpi2_debug(armci_me, "armci_send_req_msg(): send msg to server(%d), to"
                   "fwd to client %d\n", server, proc);

  return 0;
}

/*\ client sends strided data + request to server
\*/
int armci_send_req_msg_strided(int proc, request_header_t *msginfo,char *ptr,
                               int strides, int stride_arr[], int count[])
{
    int server = armci_clus_id(proc);
    int bytes;

    armci_mpi2_debug(armci_me, "armci_send_req_msg_strided: proc=%d server=%d "
                     "bytes=%d (op=%d)\n", proc, server, msginfo->datalen,
                     msginfo->operation);

    THREAD_LOCK(armci_user_threads.net_lock);

    /* we write header + descriptor of strided data  */
    bytes = sizeof(request_header_t) + msginfo->dscrlen;
    armci_send_req_msg(proc, msginfo, bytes);
    
#ifdef MPI_USER_DEF_DATATYPE
    if(strides>0) 
    {
       armci_mpi_strided2(SEND, ptr, strides, stride_arr, count, server,
                          MPI_COMM_CLIENT2SERVER);
    }
    else
#endif
    {
       /* for larger blocks write directly thus avoiding memcopy */
       armci_mpi_strided(SEND, ptr, strides, stride_arr, count, server,
                         MPI_COMM_CLIENT2SERVER);
    }
       
    THREAD_UNLOCK(armci_user_threads.net_lock);

    armci_mpi2_debug(armci_me, "armci_send_req_msg_strided(): send msg to "
                     "server(%d), to fwd to client %d\n", server, proc);

    return 0;
}

/*\ client receives data from server
\*/
char *armci_ReadFromDirect (int proc, request_header_t *msginfo, int len)
{

    int server = armci_clus_id(proc);
    MPI_Status status;

    armci_mpi2_debug(armci_me, "armci_ReadFromDirect: proc=%d, server=%d, "
                     "msginfo=%p, bytes=%d (op=%d)\n", proc, server, msginfo,
                     len, msginfo->operation);
    
    if( !(server >= 0 && server < armci_nserver) )
       armci_die("armci_ReadFromDirect: Invalid server.", 0);
    
    MPI_Check(
       MPI_Recv(msginfo + 1, len, MPI_BYTE, server, ARMCI_MPI_SPAWN_TAG,
                MPI_COMM_CLIENT2SERVER, &status)
       );

    
    armci_mpi2_debug(armci_me, "recv msg from server(%d), fwd by client %d\n",
                     server, proc);

#if MPI_SPAWN_DEBUG
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
#endif
    
 return (char *) (msginfo+1);
}

/*\ client receives strided data from server
\*/
void armci_ReadStridedFromDirect(int proc, request_header_t* msginfo,
                                 void *ptr, int strides, int stride_arr[],
                                 int count[])
{

    int server=armci_clus_id(proc);
    
    armci_mpi2_debug(armci_me, "armci_ReadStridedFromDirect: proc=%d "
                     "stride_levels=%d, server=%d bytes=%d (op=%d)\n",
                     proc, strides, server, msginfo->datalen,
                     msginfo->operation);

    
    if( !(server >= 0 && server < armci_nserver) )
       armci_die("armci_ReadStridedFromDirect: Invalid server.", 0);

#ifdef MPI_USER_DEF_DATATYPE
    if(strides > 0) 
    {
       armci_mpi_strided2(RECV, ptr, strides, stride_arr, count, server,
                          MPI_COMM_CLIENT2SERVER);
    }
    else
#endif
    {
       armci_mpi_strided(RECV, ptr, strides, stride_arr, count, server,
                         MPI_COMM_CLIENT2SERVER);
    }
}



/**
 * Platform specific server code ENDs here. (END)
 **************************************************************************/

static void armci_gather_hostnames(char **hostname_arr) 
{
    int i, j, k, namelen, is_master;
    char hostname[MPI_MAX_PROCESSOR_NAME], *hostnames=NULL;
    int *master_arr=NULL;
    
    master_arr = (int*)  malloc(armci_nproc * sizeof(int));
    hostnames  = (char*) malloc(armci_nproc * MPI_MAX_PROCESSOR_NAME *
                                   sizeof(char));
    
    if(hostnames==NULL || master_arr==NULL) 
    {
       armci_die("armci_gather_hostnames: malloc failed.", 0);
    }
    
    
    MPI_Get_processor_name(hostname, &namelen);
    MPI_Check(
       MPI_Allgather(hostname,  MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
                     hostnames, MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
                     ARMCI_COMM_WORLD)
       );

    if(armci_me == armci_master) 
    {
       is_master = 1;
    }
    else 
    {
       is_master = 0;
    }

    MPI_Check(MPI_Allgather(&is_master, 1, MPI_INT, master_arr, 1, MPI_INT,
                            ARMCI_COMM_WORLD));

    {
       /* get only the hostname of armci master processes */
       for(i=0,j=0,k=0; i<armci_nproc; i++) 
       {
          if(master_arr[i] == 1) 
          {
             if(j>=armci_nserver)
                armci_die("armci_gather_hostnames: Invalid masters.",0);
             strncpy(hostname_arr[j++], &hostnames[k], MPI_MAX_PROCESSOR_NAME);
          }
          k += MPI_MAX_PROCESSOR_NAME;
       }
    }
    
    free(hostnames);
    free(master_arr);
}


void select_server_program(char *server_program,
                           int num_procs)
{

    strcpy(server_program, (*_armci_argv)[0]);

    return;
}


static void armci_mpi2_spawn() 
{

    int i;
    char server_program[100];
    char **command_arr=NULL, **hostname_arr=NULL, **nid_arr=NULL;
    int *size_arr=NULL;
    MPI_Info *info_arr;
    
    /* we need to start 1 data server process on each node. So a total of
       "armci_nclus" data servers */
    armci_nserver = armci_nclus;
    select_server_program(server_program, armci_nserver);
    
    armci_mpi2_debug(0, "armci_mpi2_init(): Spawning %d data server processes "
                     "running %s\n", armci_nserver, server_program);

    /* allocate necessary data structures */
    {
       command_arr  = (char**)    malloc(armci_nserver * sizeof(char*));
       size_arr     = (int*)      malloc(armci_nserver * sizeof(int));
       info_arr     = (MPI_Info*) malloc(armci_nserver * sizeof(MPI_Info));
       hostname_arr = (char**)    malloc(armci_nserver * sizeof(char*));
#ifdef SPAWN_CRAY_XT
       nid_arr      = (char**)    malloc(armci_nserver * sizeof(char*));;
#endif
       for(i=0; i<armci_nserver; i++) 
       {
          hostname_arr[i] = (char*)malloc(MPI_MAX_PROCESSOR_NAME*sizeof(char));
       }

       if(command_arr==NULL || size_arr==NULL || info_arr==NULL ||
          hostname_arr==NULL) 
       {
          armci_die("armci_mpi2_spawn: malloc failed.", 0);
       }
    }
    
    /**
     * 1. root process collects hostnames (i.e. machine names) of where to
     * spawn dataservers. ARMCI masters of respective node will return their
     * hostnames. 
     */
    armci_gather_hostnames(hostname_arr);
    
       
    /** 2. initialize MPI_Comm_spawn_multiple() arguments */
    {   
       for(i=0; i<armci_nserver; i++)
       {
          command_arr[i] = (*_armci_argv)[0];  /*CHECK: path needs fix */
          size_arr[i]    = 1;                /* 1 data server in each node */
          MPI_Info_create(&info_arr[i]);
#ifdef SPAWN_CRAY_XT
          asprintf(&nid_arr[i], "%d", atoi((hostname_arr[i] + 3)));
          MPI_Info_set(info_arr[i], "host", nid_arr[i]); /*portability? */
#else
          MPI_Info_set(info_arr[i], "host", hostname_arr[i]); /*portability? */
#endif
       }
    }

    
    /**
     * 3. MPI_Comm_spawn_multiple(): This is a collective call.
     * Intercommunicator "ds_intercomm" contains only new dataserver processes.
     */
    MPI_Check(
       MPI_Comm_spawn_multiple(armci_nserver, command_arr, MPI_ARGVS_NULL,
                               size_arr, info_arr, ARMCI_ROOT, ARMCI_COMM_WORLD,
                               &MPI_COMM_CLIENT2SERVER, MPI_ERRCODES_IGNORE)
       );


    {  
       for(i=0; i<armci_nserver; i++)  free(hostname_arr[i]);
       
       free(command_arr);
       free(size_arr);
       free(info_arr);
       free(hostname_arr);
#ifdef SPAWN_CRAY_XT
       free(nid_arr);
#endif
    }
}
       
/**
 * Create server processes. This is called in armci_start_server.
 * Must be called after armci_init_clusinfo().
 */
void armci_create_server_MPIprocess ()
{
    int rank, size, flag, i;

    MPI_Initialized(&flag);
    if (flag == 0)
       armci_die("ARMCI error: MPI_Init must be called before PARMCI_Init()",0);
    
    MPI_Comm_rank(ARMCI_COMM_WORLD, &rank);
    MPI_Comm_size(ARMCI_COMM_WORLD, &size);

    /* spawn one data server process (i.e. additional MPI proc) on each node */
    armci_mpi2_spawn();
    
    /**
     * Armci masters send the following info to their corresponding server as
     * the server was not part of the initialization step in PARMCI_Init()
     *    1. cluster info ( i.e. armci_init_clusinfo() )
     *    2. lock info    ( i.e.armci_allocate_locks() )
     */
    
    if(armci_me == armci_master) {
       int msg[3];
       long shm_info[3], shmoffset;
       int shmid;
       size_t shmsize;
       
       /**
        * 1. Cluster info
        */
       msg[0] = ARMCI_MPI_SPAWN_INIT_TAG + armci_clus_me; /* for validation */
       msg[1] = armci_me;
       msg[2] = armci_clus_info[armci_clus_me].nslave;
       MPI_Send(msg, 3, MPI_INT, armci_clus_me, ARMCI_MPI_SPAWN_INIT_TAG,
                MPI_COMM_CLIENT2SERVER);

       /* send the entire clus info to its data server */
       MPI_Send(armci_clus_info, armci_nclus*sizeof(armci_clus_t), MPI_BYTE,
                armci_clus_me, ARMCI_MPI_SPAWN_INIT_TAG,
                MPI_COMM_CLIENT2SERVER);
       
       /**
        * 2. lock info
        */
       armci_get_shmem_info((char*)_armci_int_mutexes, &shmid, &shmoffset,
                            &shmsize);
       shm_info[0] = (long) shmid;
       shm_info[1] = (long) shmoffset;
       shm_info[2] = (long) shmsize;
       
       MPI_Send(shm_info, 3, MPI_LONG, armci_clus_me, ARMCI_MPI_SPAWN_INIT_TAG,
                MPI_COMM_CLIENT2SERVER);       
    }
     

    /* initialize tags for flow control */
    _armci_mpi_tag = (int*) malloc(armci_nserver*sizeof(int));
    for(i=0; i<armci_nserver; i++)
       _armci_mpi_tag[i]=ARMCI_MPI_SPAWN_TAG_BEGIN;
    
    /* makesure all processes sync here. CHECK: does it ensure global sync ? */
    MPI_Barrier(ARMCI_COMM_WORLD);
    
    armci_mpi2_debug(0, "armci_create_server_MPIprocess: Servers spawned!\n");
}

/*
  CHECK:
  =====
  1. MPI_Info portability issue. e.g. MPI_Info_set

  2. For sockets with server process option (i.e. not thread, but an actual
  process as data server), all clients call armci_init_connections() in
  armci_start_server(), which ic called by PARMCI_Init(). Should this init
  connections be done for MPI_SPAWN as well
  NOTE: armci_init_connections call _armci_buf_init(). In our MPI_SPAWN case we
  never call _armci_buf_init() either in clients or in (spawned) dataservers.

  3. Implement non-blocking in future for sure. For example:
  armci_save_strided_dscr() is disabled in armci_rem_strided (request.c) and
  armci_complete_vector_get is enabled in armci_rev_vector (request.c)

*/
  
