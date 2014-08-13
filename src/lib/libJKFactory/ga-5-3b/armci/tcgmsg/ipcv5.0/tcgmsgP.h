#ifndef TCGMSGP_H_
#define TCGMSGP_H_

#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif
#if HAVE_SIGNAL_H
#   include <signal.h>
#endif

#ifdef LAPI
#   include <lapi.h>
#endif

#include "tcgshmem.h"
#include "sndrcv.h"
#include "srftoc.h"

/* TODO autoconf way to detect this?? */
#define MAX_PROC 512
/* #define MAX_PROC 16 */
/* under Cygnus we got only serial execution */
/* #define MAX_PROC 1 */

#define INTERNAL_SYNC_TYPE 33333
#define MAX_N_OUTSTANDING_MSG 64

extern void USleep(long);
#ifndef LAPI
extern long *nxtval_shmem;
#endif

extern long DEBUG_;
extern long TCGMSG_nodeid;        /**> The id of this process */
extern long TCGMSG_nnodes;        /**> Total no. of processes */
extern char *  TCGMSG_shmem;         /**> Pointer to shared-memory segment */
extern long TCGMSG_shmem_id;      /**> ID of shared-memory segment */
extern long TCGMSG_shmem_size;    /**> Size of shared-memory segment */
extern long TCGMSG_caught_sigint; /**> True if SIGINT was trapped */

/* Structure defines shared memory buffer ... each process has
   one for every process that can send to it via shared memory.

   Adjust SHMEM_BUF_SIZE so that sizeof(ShmemBuf) is an integer
   multiple of page sizes. Structure of this buffer is exploited
   in T3D code. */

#ifdef NOTIFY_SENDER
#  ifdef LAPI
#    define RESERVED (6*sizeof(long) + sizeof(lapi_cntr_t))
#  else
#    define RESERVED 6*sizeof(long)
#  endif
#else
#  define RESERVED 4*sizeof(long)
#endif

#if defined(MACX)
#   define WHOLE_BUF_SIZE 2*65536
#elif defined(LAPI)
#   define WHOLE_BUF_SIZE (3*4096)
#else
#   define WHOLE_BUF_SIZE (16*8192)
#endif

#define SHMEM_BUF_SIZE (WHOLE_BUF_SIZE - RESERVED)

#ifdef  LAPI
#   define SND_RESERVED (4*sizeof(long) + sizeof(lapi_cntr_t) + sizeof(void*))
#   define SEND_BUF_SIZE (WHOLE_BUF_SIZE - SND_RESERVED) 
#   define SENDBUF_NUM 2
typedef struct {
    lapi_cntr_t cntr;
    void *next;
    long info[4];
    char buf[SEND_BUF_SIZE];
} sendbuf_t;
sendbuf_t *sendbuf_arr, *localbuf;
#endif

typedef struct {
    long info[4];          /**< 0=type, 1=length, 2=tag, 3=full */
    char buf[SHMEM_BUF_SIZE]; /**< Message buffer */
#ifdef NOTIFY_SENDER
    long stamp;
#   ifdef LAPI
    lapi_cntr_t cntr;
#   endif
    long flag;             /**< JN: used by receiver to signal sender */
#endif
} ShmemBuf;

/* Structure defines an entry in the send q */

typedef struct {
    long msgid;         /**< Message id for msg_status */
    long type;          /**< Message type */
    long node;          /**< Destination node */
    long tag;           /**< Message tag */
    char *buf;             /**< User or internally malloc'd buffer */
    long lenbuf;        /**< Length of user buffer in bytes */
    long written;       /**< Amount already sent */
    long buffer_number; /**< No. of buffers alread sent */
    long free_buf_on_completion; /* Boolean true if free buffer using free */
    void *next;            /**< Pointer to next entry in linked list */
    void *next_in_ring;    /**< Pointer to next entry in ring of free entries */
    long active;        /**< 0/1 if free/allocated */
} SendQEntry;

/* This structure holds basically all process specific information */

#define COMM_MODE_NONE  0
#define COMM_MODE_SHMEM 1
#define COMM_MODE_SOCK  2

typedef struct {
    ShmemBuf *sendbuf; /**< Shared-memory buffer for sending to node*/
    ShmemBuf *recvbuf; /**< Shared-memory buffer for receiving from */
    int sock;          /**< Socket for send/receive */
    int comm_mode;     /**< Defines communication info */
    pid_t pid;         /**< Unix process id (or 0 if unknown) */
    long tag_rcv;   /**< Expected tag from next rcv() */
    long n_snd;     /**< No. of messages sent from this process */
    long n_rcv;     /**< No. of messages recv from this process */
    SendQEntry *sendq; /**< Queue of messages to be sent */
} ProcInfo;

extern ProcInfo *TCGMSG_proc_info; /**< Will point to array of structures */

extern SendQEntry *TCGMSG_sendq_ring; /**< Circular ring of SendQEntry
                                        structures for fast allocation/free */

#endif /* TCGMSGP_H_ */
