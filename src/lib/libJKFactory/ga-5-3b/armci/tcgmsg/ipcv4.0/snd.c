#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/snd.c,v 1.21 2004-04-01 02:04:57 manoj Exp $ */

#include <stdio.h>
#include <stdlib.h>
#ifdef SEQUENT
#include <strings.h>
#else
#include <string.h>
#endif

#ifdef AIX
#include <sys/select.h>
#endif

#include <sys/types.h>
#include <sys/time.h>

#if defined(SHMEM) || defined(SYSV)
#   if (defined(SGI_N32) || defined(SGITFP))
#       define PARTIALSPIN
#   else
#       define NOSPIN
#   endif
#endif

#if (defined(SUN) && !defined(SOLARIS))
    extern char *sprintf();
#endif

extern void Error();

#include "sndrcv.h"
#include "sndrcvP.h"
#include "tcgsockets.h"

#ifdef GOTXDR
#include "xdrstuff.h"
#endif

#if defined(SHMEM) || defined(SYSV)
#if !defined(SEQUENT) && !defined(CONVEX)
#include <memory.h>
#endif
#include "sema.h"
#include "tcgshmem.h"
#if defined(USE_SRMOVER)
extern void SRmover();
#else
#define SRmover(a,b,n) memcpy(a,b,n)
#endif
#endif 

#ifdef EVENTLOG
#include "evlog.h"
#endif

extern void ListenOnSock(int sock);
extern int AcceptConnection(int sock);

void PrintProcInfo()
/*
  Print out the SR_proc_info structure array for this process
*/
{
  long i;

  (void) fprintf(stderr,"Process info for node %ld: \n",NODEID_());

  for (i=0; i<NNODES_(); i++)
    (void) fprintf(stderr,"[%ld] = {\n\
     clusid = %-8ld    slaveid = %-8ld      local = %-8ld\n\
       sock = %-8d      shmem = %-8p shmem_size = %-8ld\n\
   shmem_id = %-8ld     buffer = %-8p     buflen = %-8ld\n\
     header = %-8p      semid = %-8ld   sem_read = %-8ld\n\
sem_written = %-8ld      n_rcv = %-8ld     nb_rcv = %-8ld\n\
      t_rcv = %-8ld      n_snd = %-8ld     nb_snd = %-8ld\n\
      t_snd = %-8ld,     peeked = %-8ld}\n",
		   i,
		   SR_proc_info[i].clusid,
		   SR_proc_info[i].slaveid,
		   SR_proc_info[i].local,
		   SR_proc_info[i].sock,
		   SR_proc_info[i].shmem,
		   SR_proc_info[i].shmem_size,
		   SR_proc_info[i].shmem_id,
		   SR_proc_info[i].buffer,
		   SR_proc_info[i].buflen,
		   SR_proc_info[i].header,
		   SR_proc_info[i].semid,
		   SR_proc_info[i].sem_read,
		   SR_proc_info[i].sem_written,
		   SR_proc_info[i].n_rcv,
   	   (long)  SR_proc_info[i].nb_rcv,
	   (long)  SR_proc_info[i].t_rcv,
		   SR_proc_info[i].n_snd,
	   (long)  SR_proc_info[i].nb_snd,
	   (long)  SR_proc_info[i].t_snd,
		   SR_proc_info[i].peeked);
		   
  (void) fflush(stderr);
}

static void PrintMessageHeader(info, header)
     char *info;
     MessageHeader *header;
/*
  Print out the contents of a message header along with info message
*/
{
  (void) printf("%2ld:%s: type=%ld, from=%ld, to=%ld, len=%ld, tag=%ld\n",
		NODEID_(),info, header->type, header->nodefrom, 
		header->nodeto, header->length, header->tag);
  (void) fflush(stdout);
}


#if defined(SHMEM) || defined(SYSV)

static int DummyRoutine()
{int i, sum=0; for(i=0; i<10; i++) sum += i; return sum;}

static long flag(p)
     long *p;
{
#if defined(CONVEX) && defined(HPUX)
  asm("sync");
#endif

  return *p;
}

static void Await(p, value)
     long *p;
     long value;
/*
  Wait until the value pointed to by p equals value.
  Since *ptr is volatile but cannot usually declare this
  include another level of procedure call to protect
  against compiler optimization.
*/
{
  int nspin = 0;
  if (DEBUG_) {
    printf("%2ld: Await p=%p, value=%ld\n", NODEID_(), p, value);
    fflush(stdout);
  }

  for (; flag(p) != value; nspin++) {
#if defined(NOSPIN) && !defined(PARTIALSPIN)
    if (nspin < 100)
      (void) DummyRoutine();
    else 
      USleep((long) 10000);
#else
    if (nspin < 10000000)
      (void) DummyRoutine();
    else {
/*      printf("%2ld: Await sleeping\n", NODEID_()); fflush(stdout); */
      USleep((long) 100000);
    }
#endif
  }
}

static void rcv_local(type, buf, lenbuf, lenmes, nodeselect, nodefrom)
     long *type;
     char *buf;
     long *lenbuf;
     long *lenmes;
     long *nodeselect;
     long *nodefrom;
{
  long me = NODEID_();
  long node = *nodeselect;
  MessageHeader *head = SR_proc_info[node].header;
  long buflen = SR_proc_info[node].buflen;
  char *buffer = SR_proc_info[node].buffer;
  long nodeto, len;
#ifdef NOSPIN
  long semid = SR_proc_info[node].semid;
  long sem_read = SR_proc_info[node].sem_read;
  long sem_written = SR_proc_info[node].sem_written;
  long semid_to = SR_proc_info[me].semid;
  long sem_pend = SR_proc_info[me].sem_pend;
#endif
#if !defined(NOSPIN) || defined(PARTIALSPIN)
  long *buffer_full = SR_proc_info[node].buffer_full;
#endif
  
  /* Error checking */
  
  if ( (buffer == (char *) NULL) || (head == (MessageHeader *) NULL) )
    Error("rcv_local: invalid shared memory", (long) node);

#ifdef NOSPIN
  if ( (semid < 0) || (sem_read < 0) || (sem_written < 0) || 
       (semid_to < 0) || (sem_pend < 0) )
    Error("rcv_local: invalid semaphore set", (long) node);
#endif

#ifdef NOSPIN
  SemWait(semid_to, sem_pend);
#endif

  Await(&head->nodeto, me);	/* Still have this possible spin */

#if defined(NOSPIN) && !defined(PARTIALSPIN)
  SemWait(semid, sem_written);
#else
  Await(buffer_full, (long) TRUE);
#endif

  /* Now have a message for me ... check the header info and
     copy the first block of the message. */
  
  if (DEBUG_)
    PrintMessageHeader("rcv_local ",head);

  nodeto = head->nodeto;	/* Always me ... history here */
  head->nodeto = -1;

  *nodefrom = head->nodefrom;

  if (head->type != *type) {
    PrintMessageHeader("rcv_local ",head);
/*    printf("rcv_local: type mismatch ... strong typing enforced\n"); */
/*    abort(); */
    Error("rcv_local: type mismatch ... strong typing enforced", (long) *type);
  }

  *lenmes = len = head->length;

   if ( *lenmes > *lenbuf )
     Error("rcv_local: message too long for buffer", (long) *lenmes);
   if (nodeto != me)
     Error("rcv_local: message meant for someone else?", (long) nodeto);
     
  if (len)
    (void) SRmover(buf, buffer, (len > buflen) ? buflen : len);

#if defined(NOSPIN) && !defined(PARTIALSPIN)
  SemPost(semid, sem_read);
#else
  *buffer_full = FALSE;
# if defined(CONVEX) && defined(HPUX)
     asm("sync");
# endif

#endif

  len -= buflen;
  buf += buflen;

  /* Copy the remainder of the message */
  
  while (len > 0) {
#if defined(NOSPIN) && !defined(PARTIALSPIN)
    SemWait(semid, sem_written);
#else    
    Await(buffer_full, (long) TRUE);
#endif
    (void) SRmover(buf, buffer, (len > buflen) ? buflen : len);
#if defined(NOSPIN) && !defined(PARTIALSPIN)
    SemPost(semid, sem_read);
#else
    *buffer_full = FALSE;
#endif
    len -= buflen;
    buf += buflen;
  }
}

static void snd_local(type, buf, lenbuf, node)
     long *type;
     char *buf;
     long *lenbuf;
     long *node;
{
  long me = NODEID_();
  MessageHeader *head = SR_proc_info[me].header;
  long buflen = SR_proc_info[me].buflen;
  long len = *lenbuf;
  char *buffer = SR_proc_info[me].buffer;
  long tag = SR_proc_info[*node].n_snd;
#ifdef NOSPIN
  long semid = SR_proc_info[me].semid;
  long sem_read = SR_proc_info[me].sem_read;
  long sem_written = SR_proc_info[me].sem_written;
  long semid_to = SR_proc_info[*node].semid;
  long sem_pend = SR_proc_info[*node].sem_pend;
#endif
#if !defined(NOSPIN) || defined(PARTIALSPIN)
  long *buffer_full = SR_proc_info[me].buffer_full;
#endif

  /* Error checking */
  
  if ( (buffer == (char *) NULL) || (head == (MessageHeader *) NULL) )
    Error("snd_local: invalid shared memory", (long) *node);

#ifdef NOSPIN
  if ( (semid < 0) || (semid_to < 0) || (sem_read < 0) || (sem_written < 0) )
    Error("snd_local: invalid semaphore set", (long) *node);
#endif

  /* Check that final segment of last message has been consumed */

#if defined(NOSPIN) && !defined(PARTIALSPIN)
  SemWait(semid, sem_read);
#else
  Await(buffer_full, (long) FALSE);
#endif

  /* Fill in message header */
  
  head->nodefrom = (char) me;
  head->type = *type;
  head->length = *lenbuf;
  head->tag = tag;
  head->nodeto = (char) *node;
#if defined(CONVEX) && defined(HPUX)
  asm("sync");
#endif

  if (DEBUG_) {
    PrintMessageHeader("snd_local ",head);
    (void) fflush(stdout);
  }

  /* Copy the first piece of the message so that send along with
     header to minimize use of semaphores. Also need to send header
     even for messages of zero length */

  if (len)
    (void) SRmover(buffer, buf, (len > buflen) ? buflen : len);

#if defined(NOSPIN) && !defined(PARTIALSPIN)
  SemPost(semid, sem_written);
#else
  *buffer_full = TRUE;
# if defined(CONVEX) && defined(HPUX)
   asm("sync");
# endif
#endif
#ifdef NOSPIN
  SemPost(semid_to, sem_pend);
#endif

  len -= buflen;
  buf += buflen;

  while (len > 0) {
#if defined(NOSPIN) && !defined(PARTIALSPIN)
    SemWait(semid, sem_read);
#else
    Await(buffer_full, (long) FALSE);
#endif
    (void) SRmover(buffer, buf, (len > buflen) ? buflen : len);
#if defined(NOSPIN) && !defined(PARTIALSPIN)
    SemPost(semid, sem_written);
#else
    *buffer_full = TRUE;
#   if defined(CONVEX) && defined(HPUX)
       asm("sync");
#   endif
#endif
    len -= buflen;
    buf += buflen;
  }
}    
#endif

static void snd_remote(type, buf, lenbuf, node)
     long *type;
     char *buf;
     long *lenbuf;
     long *node;
/*
  synchronous send to remote process

  long *type     = user defined integer message type (input)
  char *buf      = data buffer (input)
  long *lenbuf   = length of buffer in bytes (input)
  long *node     = node to send to (input)

  for zero length messages only the header is sent
*/
{
#define SHORT_MSG_BUF_SIZE (2048 + 40)
  static char fudge[SHORT_MSG_BUF_SIZE]; 
  MessageHeader header;
  long me=NODEID_();
  int sock=SR_proc_info[*node].sock;
  long len;
#ifdef SOCK_FULL_SYNC
  char sync=0;
#endif
 
  if ( sock < 0 )
    Error("snd_remote: sending to process without socket", (long) *node);

  header.nodefrom = me;
  header.nodeto   = *node;
  header.type     = *type;
  header.length   = *lenbuf;
  header.tag      = SR_proc_info[*node].n_snd;

  /* header.length is the no. of items if XDR is used or just the
     number of bytes */

#ifdef GOTXDR
  if ( *type & MSGDBL )
    header.length = *lenbuf / sizeof(double);
  else if ( *type & MSGINT ) 
    header.length = *lenbuf / sizeof(long);
  else if ( *type & MSGCHR )
    header.length = *lenbuf / sizeof(char);
  else
    header.length = *lenbuf;
#else
  header.length = *lenbuf;
#endif

  if (DEBUG_)
    PrintMessageHeader("snd_remote",&header);

#ifndef GOTXDR
  /* Combine header and messages less than a certain size to avoid
   * performance problem on (older?) linuxes */
  if ((*lenbuf + sizeof(header)) <= sizeof(fudge)) {
    memcpy(fudge,(char *) &header, sizeof(header));
    memcpy(fudge+sizeof(header), buf, *lenbuf);
    if ( (len = WriteToSocket(sock, fudge, sizeof(header)+*lenbuf)) != 
	                                  ((long)sizeof(header)+*lenbuf))
      Error("snd_remote: writing message to socket",
            (long) (len+100000*(sock + 1000* *node)));
    return;
  }
#endif

#ifdef GOTXDR
  (void) WriteXdrLong(sock, (long *) &header, 
		      (long) (sizeof(header)/sizeof(long)));
#else
  if ( (len = WriteToSocket(sock, (char *) &header, (long) sizeof(header)))
                                                        != sizeof(header) )
      Error("snd_remote: writing header to socket", len);
#endif

  if (*lenbuf)  {
#ifdef GOTXDR
    if ( *type & MSGDBL )
      (void) WriteXdrDouble(sock, (double *) buf, header.length);
    else if ( *type & MSGINT )
      (void) WriteXdrLong(sock, (long *) buf, header.length);
    else if ( *type & MSGCHR )
      (void) WriteXdrChar(sock, (char *) buf, header.length);
    else if ( (len = WriteToSocket(sock, buf, header.length)) != 
                                                      header.length)
      Error("snd_remote: writing message to socket",
            (long) (len+100000*(sock + 1000* *node)));
#else
    if ( (len = WriteToSocket(sock, buf, header.length)) != 
	                                 header.length)
      Error("snd_remote: writing message to socket",
            (long) (len+100000*(sock + 1000* *node)));
#endif
  }

#ifdef SOCK_FULL_SYNC
  /* this read (and write in rcv_remote) of an acknowledgment 
     forces synchronous */

  if ( ReadFromSocket(sock, &sync, (long) 1) != 1)
    Error("snd_remote: reading acknowledgement",
	  (long) (len+100000*(sock + 1000* *node)));
#endif
}

/*ARGSUSED*/
void SND_(type, buf, lenbuf, node, sync)
     long *type;
     void *buf;
     long *lenbuf;
     long *node;
     long *sync;
/*
  mostly syncrhonous send

  long *type     = user defined integer message type (input)
  void *buf     = data buffer (input)
  long *lenbuf   = length of buffer in bytes (input)
  long *node     = node to send to (input)
  long *sync    = flag for sync/async ... IGNORED

  for zero length messages only the header is sent
*/
{
  long me=NODEID_();
  long nproc=NNODES_();
#ifdef TIMINGS
  double start;
#endif

  /* Error checking */

  if (*node == me)
    Error("SND_: cannot send message to self", (long) me);

  if ( (*node < 0) || (*node > nproc) )
    Error("SND_: out of range node requested", (long) *node);

  if ( (*lenbuf < 0) || (*lenbuf > (long)BIG_MESSAGE) )
    Error("SND_: message length out of range", (long) *lenbuf);

#ifdef EVENTLOG
  evlog(EVKEY_BEGIN,     EVENT_SND,
	EVKEY_MSG_LEN,  (int) *lenbuf,
	EVKEY_MSG_FROM, (int)  me,
	EVKEY_MSG_TO,   (int) *node,
	EVKEY_MSG_TYPE, (int) *type,
	EVKEY_MSG_SYNC, (int) *sync,
	EVKEY_LAST_ARG);
#endif

  /* Send via shared memory or sockets */

#ifdef TIMINGS
  start = TCGTIME_();
#endif

#if defined(SHMEM) || defined(SYSV)
  if (SR_proc_info[*node].local){
#ifdef KSR_NATIVE
    KSR_snd_local(type, buf, lenbuf, node);
#else
    snd_local(type, buf, lenbuf, node);
#endif
  } else {
#endif
      snd_remote(type, buf, lenbuf, node);
#if defined(SHMEM) || defined(SYSV)
  }
#endif

  /* Collect statistics */

  SR_proc_info[*node].n_snd += 1;
  SR_proc_info[*node].nb_snd += *lenbuf;

#ifdef TIMINGS
  SR_proc_info[*node].t_snd += TCGTIME_() - start;
#endif

#ifdef EVENTLOG
  evlog(EVKEY_END, EVENT_SND, EVKEY_LAST_ARG);
#endif
}    
    
static long MatchMessage(header, me, type)
     MessageHeader *header;
     long me, type;
/*
  Wrapper round check on if header is to me and of required
  type so that compiler does not optimize out fetching
  header info from shared memory.
*/
{
  return (long) ((header->nodeto == me) && (header->type == type));
}

static long NextReadyNode(type)
     long type;
/*
  Select a node from which input is pending ... also match the
  desired type.

  next_node is maintained as the last node that NextReadyNode chose
  plus one modulo NNODES_(). This aids in ensuring fairness.

  First use select to get info about the sockets and then loop
  through processes looking either at the bit in the fd_set for
  the socket (remote process) or the message header in the shared
  memory buffer (local process).

  This may be an expensive operation but fairness seems important.

  If only sockets are in use, just block in select until data is
  available.  
*/
{
  static long  next_node = 0;

  long  nproc = NNODES_();
  long  me = NODEID_();
  int i, nspin = 0;

  if (!SR_using_shmem) {
    int list[MAX_PROCESS];
    int nready;
    nready = WaitForSockets(SR_nsock,SR_socks,list);
    if (nready == 0) 
      Error("NextReadyNode: nready = 0\n", 0);

    /* Insert here type checking logic ... not yet done */

    return SR_socks_proc[list[0]];
  }

  /* With both local and remote processes end up with a busy wait
     as no way to wait for both a semaphore and a socket.
     Moderate this slightly by having short timeout in select */

  while (1) {
    
    for(i=0; i<nproc; i++, next_node = (next_node + 1) % nproc) {

      if (next_node == me) {
        ;  /* can't receive from self */
      }
      else if (SR_proc_info[next_node].local) {
	/* Look for local message */

#ifdef KSR_NATIVE
	if (KSR_MatchMessage(next_node, me, type))
#else
	if (MatchMessage(SR_proc_info[next_node].header, me, type))
#endif
	  break;
      }
      else if (SR_proc_info[next_node].sock >= 0) {
	/* Look for message over socket */

	int sock = SR_proc_info[next_node].sock;

	/* Have we already peeked at this socket? */

	if (SR_proc_info[next_node].peeked) {
	  if (SR_proc_info[next_node].head_peek.type == type)
	    break;
	}
	else if (PollSocket(sock)) {
	  /* Data is available ... let's peek at it */
#ifdef GOTXDR
	  (void) ReadXdrLong(sock, 
			     (long *) &SR_proc_info[next_node].head_peek,
			     (long) (sizeof(MessageHeader)/sizeof(long)));
#else
	  if (ReadFromSocket(sock, 
			     (char *) &SR_proc_info[next_node].head_peek,
			     (long) sizeof(MessageHeader))
	      != sizeof(MessageHeader) )
	    Error("NextReadyNode: reading header from socket", next_node);
#endif
	  SR_proc_info[next_node].peeked = TRUE;
	  if (DEBUG_)
	    PrintMessageHeader("peeked_at ",
			       &SR_proc_info[next_node].head_peek);

	  if (SR_proc_info[next_node].head_peek.type == type)
	    break;
	}
      }
    }
    if (i < nproc)       /* If found a node skip out of the while loop */
      break;

    nspin++;		 /* Compromise between low latency and low cpu use */
    if (nspin < 10)
      continue;
    else if (nspin < 100)
      USleep((long) 1000);
    else if (nspin < 600)
      USleep((long) 10000);
    else
      USleep((long) 100000);
  }

  i = next_node;
  next_node = (next_node + 1) % nproc;
  
  return (long) i;
}

long PROBE_(type, node)
     long *type, *node;
     /*
       Return 1/0 (TRUE/FALSE) if a message of the given type is available
       from the given node.  If the node is specified as -1, then all nodes
       will be examined.  Some attempt is made at ensuring fairness.
       
       First use select to get info about the sockets and then loop
       through processes looking either at the bit in the fd_set for
       the socket (remote process) or the message header in the shared
       memory buffer (local process).
       
       This may be an expensive operation but fairness seems important.
       */
{
  long  nproc = NNODES_();
  long  me = NODEID_();
  int i, proclo, prochi;

  if (*node == me)
    Error("PROBE_ : cannot recv message from self, msgtype=", *type);
      
      if (*node == -1) {		/* match anyone */
	proclo = 0;
	prochi = nproc-1;
      }
      else
	proclo = prochi = *node;
  
  for(i=proclo; i<=prochi; i++) {
    
    if (i == me) {
      ;  /* can't receive from self */
    }
    else if (SR_proc_info[i].local) {
      /* Look for local message */
      
#ifdef KSR_NATIVE
      if (KSR_MatchMessage(i, me, type))
#else
      if (MatchMessage(SR_proc_info[i].header, me, *type))
#endif
	break;
    }
    else if (SR_proc_info[i].sock >= 0) {
      /* Look for message over socket */
      
      int sock = SR_proc_info[i].sock;
      
      /* Have we already peeked at this socket? */
      
      if (SR_proc_info[i].peeked) {
	if (SR_proc_info[i].head_peek.type == *type)
	  break;
      }
      else if (PollSocket(sock)) {
	/* Data is available ... let's peek at it */
#ifdef GOTXDR
	(void) ReadXdrLong(sock, 
			   (long *) &SR_proc_info[i].head_peek,
			   (long) (sizeof(MessageHeader)/sizeof(long)));
#else
	if (ReadFromSocket(sock, 
			   (char *) &SR_proc_info[i].head_peek,
			   (long) sizeof(MessageHeader))
	    != sizeof(MessageHeader) )
	  Error("NextReadyNode: reading header from socket", (long) i);
#endif
	SR_proc_info[i].peeked = TRUE;
	if (DEBUG_)
	  PrintMessageHeader("peeked_at ",
			     &SR_proc_info[i].head_peek);
	
	if (SR_proc_info[i].head_peek.type == *type)
	  break;
      }
    }
  }

  if (i <= prochi)
    return 1;
  else
    return 0;
}


static void rcv_remote(type, buf, lenbuf, lenmes, nodeselect, nodefrom)
     long *type;
     char *buf;
     long *lenbuf;
     long *lenmes;
     long *nodeselect;
     long *nodefrom;
/*
  synchronous receive of data

  long *type        = user defined type of received message (input)
  char *buf        = data buffer (output)
  long *lenbuf      = length of buffer in bytes (input)
  long *lenmes      = length of received message in bytes (output)
                     (exceeding receive buffer is hard error)
  long *nodeselect  = node to receive from (input)
                     -1 implies that any pending message may be received
                       
  long *nodefrom    = node message is received from (output)
*/
{
  long me = NODEID_();
  long node = *nodeselect;
  int sock = SR_proc_info[node].sock;
  long len;
  MessageHeader header;
#ifdef SOCK_FULL_SYNC
  char sync = 0;
#endif

  if ( sock < 0 )
    Error("rcv_remote: receiving from process without socket", (long) node);

  /* read the message header and check contents */

  if (SR_proc_info[node].peeked) {
    /* Have peeked at this socket ... get message header from buffer */

    if (DEBUG_)
      printf("%2ld: rcv_remote message has been peeked at\n", me);

    (void) memcpy((char *) &header, (char *) &SR_proc_info[node].head_peek,
		  sizeof(MessageHeader));
    SR_proc_info[node].peeked = FALSE;
  }
  else {
#ifdef GOTXDR
    (void) ReadXdrLong(sock, (long *) &header,
		       (long) (sizeof(header)/sizeof(long)));
#else
    if ( (len = ReadFromSocket(sock, (char *) &header, (long) sizeof(header)))
	!= sizeof(header) )
      Error("rcv_remote: reading header from socket", len);
#endif
  }

  if (DEBUG_)
    PrintMessageHeader("rcv_remote",&header);

  if (header.nodeto != me) {
    PrintMessageHeader("rcv_remote",&header);
    Error("rcv_remote: got message meant for someone else",
	  (long) header.nodeto);
  }

  *nodefrom = header.nodefrom;
  if (*nodefrom != node)
    Error("rcv_remote: got message from someone on incorrect socket",
	  (long) *nodefrom);

  if (header.type != *type) {
    PrintMessageHeader("rcv_remote",&header);
    printf("rcv_remote: type mismatch ... strong typing enforced\n");
    abort();
    Error("rcv_remote: type mismatch ... strong typing enforced", (long) *type);
  }

#ifdef GOTXDR
  if ( *type & MSGDBL )
    *lenmes = header.length * sizeof(double);
  else if ( *type & MSGINT )
    *lenmes = header.length * sizeof(long);
  else if ( *type & MSGCHR )
    *lenmes = header.length * sizeof(char);
  else
    *lenmes = header.length; 
#else
  *lenmes = header.length; 
#endif
  
  if ( (*lenmes < 0) || (*lenmes > (long)BIG_MESSAGE) || (*lenmes > *lenbuf) ) {
    PrintMessageHeader("rcv_remote",&header);
    (void) fprintf(stderr, "rcv_remote err: lenbuf=%ld\n",*lenbuf);
    Error("rcv_remote: message length out of range",(long) *lenmes);
  }

  if (*lenmes > 0) {
#ifdef GOTXDR
    if ( *type & MSGDBL )
      (void) ReadXdrDouble(sock, (double *) buf, header.length);
    else if ( *type & MSGINT ) 
      (void) ReadXdrLong(sock, (long *) buf, header.length);
    else if ( *type & MSGCHR )
      (void) ReadXdrChar(sock, (char *) buf, header.length);
    else if ( (len = ReadFromSocket(sock, buf, *lenmes)) != *lenmes)
      Error("rcv_remote: reading message from socket",
            (long) (len+100000*(sock+ 1000* *nodefrom)));
#else
    if ( (len = ReadFromSocket(sock, buf, *lenmes)) != *lenmes)
      Error("rcv_remote: reading message from socket",
            (long) (len+100000*(sock+ 1000* *nodefrom)));
#endif
  }

  /* this write (and read in snd_remote) makes the link synchronous */

#ifdef SOCK_FULL_SYNC
  if ( WriteToSocket(sock, &sync, (long) 1) != 1)
    Error("rcv_remote: writing sync to socket", (long) node);
#endif

}

/*ARGSUSED*/
void RCV_(type, buf, lenbuf, lenmes, nodeselect, nodefrom, sync)
     long *type;
     void *buf;
     long *lenbuf;
     long *lenmes;
     long *nodeselect;
     long *nodefrom;
     long *sync;
/*
  long *type        = user defined type of received message (input)
  void *buf        = data buffer (output)
  long *lenbuf      = length of buffer in bytes (input)
  long *lenmes      = length of received message in bytes (output)
                     (exceeding receive buffer is hard error)
  long *nodeselect  = node to receive from (input)
                     -1 implies that any pending message may be received
                       
  long *nodefrom    = node message is received from (output)
  long *sync        = 0 for asynchronous, 1 for synchronous (NOT USED)
*/
{
  long me = NODEID_();
  long nproc = NNODES_();
  long node;
#ifdef TIMINGS
  double start;
#endif

#ifdef EVENTLOG
  evlog(EVKEY_BEGIN,     EVENT_RCV,
	EVKEY_MSG_FROM, (int) *nodeselect,
	EVKEY_MSG_TO,   (int)  me,
	EVKEY_MSG_TYPE, (int) *type,
	EVKEY_MSG_SYNC, (int) *sync,
	EVKEY_LAST_ARG);
#endif

  /* Assign the desired node or the next ready node */

#ifdef TIMINGS
  start = TCGTIME_();
#endif

  if (*nodeselect == -1)
    node = NextReadyNode(*type);
  else
    node = *nodeselect;

  /* Check for some errors ... need more checking here ...
     note that the overall master process has id nproc */

  if (node == me)
    Error("RCV_: cannot receive message from self", (long) me);

  if ( (node < 0) || (node > nproc) )
    Error("RCV_: out of range node requested", (long) node);

  /* Receive the message ... use shared memory, switch or socket */

#if defined(SHMEM) || defined(SYSV)
  if (SR_proc_info[node].local){
#ifdef KSR_NATIVE
    KSR_rcv_local(type, buf, lenbuf, lenmes, &node, nodefrom);
#else
    rcv_local(type, buf, lenbuf, lenmes, &node, nodefrom);
#endif
  } else {
#endif
      rcv_remote(type, buf, lenbuf, lenmes, &node, nodefrom);
#if defined(SHMEM) || defined(SYSV)
  }
#endif

  /* Collect statistics */

  SR_proc_info[node].n_rcv += 1;
  SR_proc_info[node].nb_rcv += *lenmes;

#ifdef TIMINGS
  SR_proc_info[node].t_rcv += TCGTIME_() - start;
#endif

#ifdef EVENTLOG
  evlog(EVKEY_END, EVENT_RCV,
	EVKEY_MSG_FROM, (int) node,
	EVKEY_MSG_LEN, (int) *lenmes,
	EVKEY_LAST_ARG);
#endif
}    
  
void RemoteConnect(a, b, c)
     long a, b, c;
/*
  Make a socket connection between processes a and b via the
  process c to which both are already connected.
*/
{
  long me = NODEID_();
  long nproc = NNODES_();
  long type = TYPE_CONNECT;  /* Overriden below */
  char cport[8];
  long tmp, lenmes, nodefrom, clusid, lenbuf, sync=1;
  int sock, port;
  long lport;

  if ((a == b) || (a == c) || (b == c) )
    return;        /* Gracefully ignore redundant connections */

  if ( (me != a) && (me != b) && (me != c) )
    return;        /* I'm not involved in this connection */
    

  if (a < b) {
    tmp = a; a = b; b = tmp;
  }

  type = (a + nproc*b) | MSGINT;  /* Create a unique type */

  if (DEBUG_) {
    (void) printf("RC a=%ld, b=%ld, c=%ld, me=%ld\n",a,b,c,me);
    (void) fflush(stdout);
  }

  if (a == me) {
    CreateSocketAndBind(&sock, &port);  /* Create port */
    if (DEBUG_) {
      (void) printf("RC node=%ld, sock=%d, port=%d\n",me, sock, port);
      (void) fflush(stdout);
    }
    lport = port;
    lenbuf = sizeof lport;
    ListenOnSock(sock);
    SND_(&type, (char *) &lport, &lenbuf, &c, &sync); /* Port to intermediate */
    SR_proc_info[b].sock = AcceptConnection(sock); /* Accept connection
						     and save socket info */
  }
  else if (b == me) {
    clusid = SR_proc_info[a].clusid;
    lenbuf = sizeof lport;
    RCV_(&type, (char *) &lport, &lenbuf, &lenmes, &c, &nodefrom, &sync);
    port = lport;
    (void) sprintf(cport,"%d",port);
    lenbuf = strlen(cport) + 1;
    if (lenbuf > (long)sizeof(cport))
      Error("RemoteConnect: cport too small", (long) lenbuf);
    SR_proc_info[a].sock = 
      CreateSocketAndConnect(SR_clus_info[clusid].hostname, cport); 
  }
  else if (c == me) {
    lenbuf = sizeof lport;
    RCV_(&type, (char *) &lport, &lenbuf, &lenmes, &a, &nodefrom, &sync);
    SND_(&type, (char *) &lport, &lenbuf, &b, &sync);
  }
}
