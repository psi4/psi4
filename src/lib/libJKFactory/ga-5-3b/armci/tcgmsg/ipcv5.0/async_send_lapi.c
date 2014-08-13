#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif

#include "tcgmsgP.h"

/*#define DEBUG 1*/
/*#define DEBUG2 1*/
static const long false = 0;
static const long true  = 1;

typedef struct {
    int from:16;
    int to:16;
} nodepair_t;

typedef union{
    long fromto;
    nodepair_t n;
} pair_t;


extern void Busy(int);

/* All data movement to/from shared memory is done using the
   COPY_TO/FROM_SHMEM macros */
extern lapi_handle_t lapi_handle;
/* ShmemBuf *localbuf = &tmp_snd_buf; */
extern void lapi_put_c(void* dest, void* src, long bytes, long node, lapi_cntr_t *cntr);
extern void lapi_put(void* dest, void* src, long bytes, long node);
extern void lapi_get(void* dest, void* src, long bytes, long node);
#define COPY_TO_LOCAL(src, dest, n) (void) memcpy(dest, src, (long) n)
#define COPY_FROM_LOCAL(src, dest, n) (void)memcpy(dest, src, (long) n)
#define COPY_TO_REMOTE(src,dest,n,node) lapi_put(dest, src, (long) n, node)
#define COPY_FROM_REMOTE(src,dest,n,node)lapi_get(dest, src, (long) n,node)
/* #define COPY_TO_REMOTE_CNTR(src, dest, n, node, pcntr) lapi_put_c(dest, src, (long) n, node, pcntr) */
#define COPY_TO_REMOTE_CNTR(localbuf, dest, n, node, pcntr) do { \
    if (LAPI_Put(lapi_handle,(uint)node, (uint)n, dest,localbuf->info, pcntr, &localbuf->cntr, NULL)) { \
        Error("TCG:lapi_put_c failed",0); \
    } \
} while (0)
#define NEXT_LOC_BUF(localbuf) localbuf = (sendbuf_t*)localbuf->next;
#define GET_LOC_BUF(localbuf) do { \
    if(LAPI_Waitcntr(lapi_handle, &localbuf->cntr, 1, NULL)) { \
        Error("TCG:LAPI_Waitcntr failed",0); \
    } \
} while (0)
#ifndef FLUSH_CACHE 
#   define FLUSH_CACHE          
#endif
#ifndef FLUSH_CACHE_LINE
#   define FLUSH_CACHE_LINE(x) 
#endif

#define TCG_ABS(a) (((a) >= 0) ? (a) : (-(a)))


/**
 * Return the value of a volatile variable in shared memory
 * that is REMOTE to this processor
 */
static long remote_flag(long *p, long node)
{
    long tmp;

    /*  FLUSH_CACHE;*/ /* no need to flush for one word only*/
    COPY_FROM_REMOTE(p, &tmp, sizeof(tmp), node);
    return tmp;
}


/**
 * Return the value of a volatile variable in shared memory
 * that is LOCAL to this processor
 */
static long local_flag(void *p)
{
    long val;
    FLUSH_CACHE_LINE(p);  
    val = *(long*)p;
    return(val);
}


void set_local_flag(void *p, long val)
{
    *(long*)p = val;
}


void set_remote_flag(void *p, long val, long node)
{
    COPY_TO_REMOTE(&val, p, sizeof(long), node); 
}


/**
 * Wait on Lapi counter for data to appear
 * check if *p == value
 */
static void lapi_await(long *p, long value, lapi_cntr_t* cntr)
{
    int val;
    long pval;

    if(LAPI_Waitcntr(lapi_handle, cntr, 1, &val))
        Error("lapi_await: error",-1);

#if 0
    if ( (pval = local_flag(p)) != value) {
        fprintf(stdout,"%2ld: invalid value=%ld, local_flag=%lx %ld\n",
                TCGMSG_nodeid, value, (unsigned long)p, pval);
        fflush(stdout);
        Error("lapi_await: exiting..",-1);;
    }
#endif
}


/**
 * Wait for (*p == value)
 */
static void local_await(long *p, long value)
{
    long pval;
    long nspin = 0;
    long spinlim = 100000000;
    long waittim = 100000;
    extern void flush_send_q(void);

    while ((pval = local_flag(p)) != value) {

        if (pval && (pval != value)) {
            fprintf(stdout,"%2ld: invalid value=%ld, local_flag=%lx %ld\n", 
                    TCGMSG_nodeid, value, (unsigned long)p, pval);
            fflush(stdout);
            exit(1);
        }
        nspin++;
        if((nspin&7)==0)flush_send_q();
        if (nspin < spinlim)
            Busy(100);
        else 
            USleep(waittim);
    }
}


/**
 * Entry points to info about a message ... determine which
 * transport mechanism is appropriate and send as much as
 * possible without blocking.
 *
 * Right now just shared memory ... when sockets are working this
 * routine will become async_shmem_send.
 *
 * Shared-memory protocol aims for low latency. Each process has
 * one buffer for every other process.  Thus, to send a message U
 * merely have to determine if the receivers buffer for you is empty
 * and copy directly into the receivers buffer.
 *
 * Return 0 data has not been sent, 1 if the send is complete.
 */
long async_send(SendQEntry *entry)
{
    long node = entry->node;
    ShmemBuf *sendbuf= TCGMSG_proc_info[node].sendbuf;
#ifdef NOTIFY_SENDER
    void *busy_flag = &TCGMSG_proc_info[node].recvbuf->flag;
#endif
    long ncopy, complete;
    long pval;
    long info[4];
    pair_t pair;

#ifdef DEBUG2
    (void) fprintf(stdout,"%2ld: sending to %ld buf=%lx len=%ld\n",
                   TCGMSG_nodeid, node, entry->buf, entry->lenbuf); 
    (void) fprintf(stdout,"%2ld: sendbuf=%lx\n", TCGMSG_nodeid, sendbuf);
    (void) fflush(stdout);
#endif

    /* return if the receiver buffer is not available */
#ifdef NOTIFY_SENDER
    pval = local_flag(busy_flag);
#else
    pval = remote_flag(&sendbuf->info[3], node);
#endif
    if (pval) {
#ifdef DEBUG
        {
            long info[4];
            FLUSH_CACHE;
            COPY_FROM_REMOTE(sendbuf->info, info, sizeof(info), node);
            fprintf(stdout,"%2ld: snd info after full = %ld %ld %ld\n",
                    TCGMSG_nodeid, info[0], info[1], info[2]);
            fflush(stdout);
            sleep(1);
        }
#endif

        return 0;
    }

    /* if data has been written already and we are here, operation is complete */
    if(entry->written) return 1L;

#ifdef NOTIFY_SENDER
    set_local_flag(busy_flag,true);
#endif

    info[0] = entry->type; info[1] = entry->lenbuf; info[2] = entry->tag;
#if 0
    entry->buffer_number++;
    info[3] = entry->buffer_number;
#else
    pair.n.from = TCGMSG_nodeid;
    pair.n.to   = node;
    info[3] =  pair.fromto;
#endif

    /* Copy over the message if it fits in the receiver buffer */
    ncopy = (long) (( entry->lenbuf <= SHMEM_BUF_SIZE) ? entry->lenbuf : 0 );

    GET_LOC_BUF(localbuf);

    if (ncopy) {
#ifdef DEBUG
        printf("%ld:snd:copying data node=%ld adr=%lx %ld bytes\n",
                TCGMSG_nodeid, node, sendbuf->buf, ncopy);
        fflush(stdout);
#endif
        COPY_TO_LOCAL(entry->buf+entry->written, localbuf->buf, ncopy);
        complete = 1;
    } else {
#ifdef DEBUG
        printf("%ld:snd:copying addr node=%ld adr=%lx %ld bytes\n",
                TCGMSG_nodeid, node, sendbuf->buf, ncopy);
        fflush(stdout);
#endif
        /* copy address of the user buffer to the send buffer */
        COPY_TO_LOCAL(&(entry->buf), localbuf->buf, sizeof(char*));
        ncopy = sizeof(char*);
        complete = 0;  /* sent is complete only when receiver gets the data */
        entry->written = 1; 
    }

#ifdef DEBUG
    printf("%ld:snd:copying info to node=%ld adr=%lx %ld bytes\n",
            TCGMSG_nodeid, node, sendbuf->info, sizeof(info));
    fflush(stdout);
#endif

    COPY_TO_LOCAL(info, localbuf->info, sizeof(info));
    COPY_TO_REMOTE_CNTR(localbuf,sendbuf,sizeof(info)+ncopy,node,&sendbuf->cntr); 

    /* advance to next buf */
    NEXT_LOC_BUF(localbuf);

    return complete;
}


/**
 * Receive a message of given type from the specified node, returning
 * the message and length of the message.
 * 
 * Right now just shared memory ... when sockets are working this
 * routine will become msg_shmem_rcv
 *
 * Shared-memory protocol aims for low latency. Each process has
 * one buffer for every other process.  Thus, to send a message U
 * merely have to determine if the receivers buffer for you is empty
 * and copy directly into the receivers buffer.
 */
void msg_rcv(long type, char *buf, long lenbuf, long *lenmes, long node)
{
    long me = TCGMSG_nodeid;
    ShmemBuf *recvbuf;        /* Points to receving buffer */
    long nleft;
    long msg_type, msg_tag, msg_len;
    long buffer_number = 1;
    long expected_tag = TCGMSG_proc_info[node].tag_rcv++;
#ifdef NOTIFY_SENDER
    void *busy_flag= &TCGMSG_proc_info[node].sendbuf->flag;
#endif

    if (node<0 || node>=TCGMSG_nnodes)
        Error("msg_rcv: node is out of range", node);

    recvbuf = TCGMSG_proc_info[node].recvbuf;  

    /* Wait for first part message to be written */

#ifdef DEBUG
    (void) fprintf(stdout,"%2ld: receiving from %ld buf=%lx len=%ld\n",
                   me, node, recvbuf,lenbuf); 
    (void) fprintf(stdout,"%2ld: user buf=%lx len=%ld\n", me, buf, lenbuf);
    (void) fflush(stdout);
#endif


#ifdef LAPI
    lapi_await(&recvbuf->info[3], buffer_number, &recvbuf->cntr);
#else
    local_await(&recvbuf->info[3], buffer_number);
#endif

    /* Copy over the header information */

    msg_type = recvbuf->info[0]; 
    msg_len  = recvbuf->info[1];
    msg_tag  = recvbuf->info[2];

#ifdef DEBUG
    (void) fprintf(stdout,"%2ld: received msg from %ld len=%ld\n",
                   me, node, msg_len); 
    (void) fflush(stdout);
#endif

    /* Check type and size information */
    if(msg_tag != expected_tag) {
        pair_t pair;
        pair.fromto = recvbuf->info[3];
        fprintf(stdout,
                "rcv: me=%ld from=%ld type=%ld expectedtag=%ld lenbuf=%ld\ngot: to=%d from=%d type=%ld msg_tag=%ld msg_len=%ld info[3]=%ld\n",
                me, node, type, expected_tag, lenbuf,
                (int)pair.n.to, (int)pair.n.from, msg_type, msg_tag, msg_len,
                recvbuf->info[3]);
        fflush(stdout);
        Error("msg_rcv: tag mismatch ... transport layer failed????", 0L);
    }

    if (msg_type != type) {
        (void) fprintf(stderr,
                       "rcv: me=%ld from=%ld type=(%ld != %ld) tag=%ld len=%ld\n",
                       me, node, type, msg_type, msg_tag, msg_len);
        Error("msg_rcv: type mismatch ... strong typing enforced\n", 0L);
    }

    if (msg_len > lenbuf) {
        (void) fprintf(stderr,
                       "rcv: me=%ld from=%ld type=%ld tag=%ld len=(%ld > %ld)\n",
                       me, node, type, msg_tag, msg_len, lenbuf);
        Error("msg_rcv: message too long for buffer\n", 0L);
    }

    nleft = *lenmes = msg_len;

    if (nleft) {
        long ncopy = nleft;

        /* for short messages data is in local buffer, for long in remote buffer */

        if(nleft <= SHMEM_BUF_SIZE) { 

            FLUSH_CACHE;

            COPY_FROM_LOCAL(recvbuf->buf, buf, ncopy);

        }else {

            char *addr = *((char**)recvbuf->buf);

            COPY_FROM_REMOTE(addr, buf, nleft, node);

        }
    }

    recvbuf->info[3] = false;
#ifdef NOTIFY_SENDER
    /* confirm that data has been transfered */
    set_remote_flag(busy_flag,false,node);
#endif

}


long MatchShmMessage(long node, long type)
{
    ShmemBuf *recvbuf;
    long  msg_type;

    recvbuf = TCGMSG_proc_info[node].recvbuf;

    if(recvbuf->info[3] == false) return (0); /* no message to receive */

    /* we have a message but let's see if want it */

    FLUSH_CACHE_LINE(recvbuf->info);
    COPY_FROM_LOCAL(recvbuf->info, &msg_type, sizeof(long));
    if(type == msg_type) return (1);
    return (0);
}
