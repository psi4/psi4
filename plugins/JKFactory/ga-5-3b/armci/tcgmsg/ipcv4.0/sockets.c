#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/sockets.c,v 1.12 2005-04-08 16:55:04 vinodtipparaju Exp $ */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef SEQUENT
#include <strings.h>
#else
#include <string.h>
#endif
#if defined(SUN) || defined(ALLIANT) || defined(ENCORE) || \
                    defined(SEQUENT) || defined(AIX)    || \
                    defined(NEXT)    || defined(LINUX)
#include <sys/wait.h>
#endif

#ifdef AIX
#include <sys/select.h>
#endif
#ifdef CONVEX
#include <errno.h>
#else
#include <sys/errno.h>
#endif
#include <sys/time.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <netdb.h>

#ifdef CRAY
#include <memory.h>
#include <errno.h>
#else
extern int errno;
#endif

#include "sndrcv.h"
#include "sndrcvP.h"

long WaitForSockets(int nsock, int *socks, int *list)
/*
  Wait until one or more sockets are ready or have an exception.

  Returns the number of ready sockets and sets corresponding 
  numbers in list.  I.e., list[i]=k meaning sock[k] is ready.
*/
{
  fd_set ready;
  int i;
  long nready;
  int sockmax = 0;

again:
  FD_ZERO(&ready);
  for (i=0; i<nsock; i++) {
    FD_SET(socks[i], &ready);
    if (socks[i] > sockmax) sockmax = socks[i];
  }
  nready = (long) select(sockmax+1, &ready, (fd_set *) NULL, (fd_set *) NULL,
			 (struct timeval *) NULL);
  if (nready < 0) {
    if (errno == EINTR) {
      /*fprintf(stderr,"wait in sockets got interrupted\n");*/
      goto again;
    }
    else {
      Error("WaitForSockets: error from select",  0L);
    }
  }
  else {
    int n = 0;
    for (i=0; i<nsock; i++) {
      if (FD_ISSET(socks[i],&ready)) list[n++] = i;
    }
  }	

  return nready;
}


long PollSocket(sock)
     int sock;
/*
  Poll the socket for available input.

  Return 1 if data is available, 0 otherwise.
*/
{
  fd_set ready;
  struct timeval timelimit;
  int nready;

  if (sock < 0)
    return 0;

again:
  FD_ZERO(&ready);
  FD_SET(sock, &ready);
  timelimit.tv_sec = 0;
  timelimit.tv_usec = 0;

  nready = (long) select(sock+1, &ready, (fd_set *) NULL, (fd_set *) NULL,
			 &timelimit);
  if (nready < 0) {
    if (errno == EINTR)
      goto again;
    else
      Error("PollSocket: error from select",  (long) sock);
  }

  return nready;
}

void TcpNoDelay(sock)
  int sock;
/*
  Turn off waiting for more input to improve buffering 
  by TCP layer ... improves performance for small messages by
  a factor of 30 or more. Slightly degrades performance for
  large messages.
*/
{
  int status, level, value=1;
#ifdef AIX
  struct protoent *proto = getprotobyname("tcp");
#else
  struct protoent *proto = getprotobyname("TCP");
#endif

#if defined(APOLLO)
  if (value)
    return;
#endif

  if (proto == (struct protoent *) NULL)
    Error("TcpNoDelay: getprotobyname on TCP failed!", (long) -1);

  level = proto->p_proto;

  status = setsockopt(sock, level, TCP_NODELAY, &value, sizeof(int));

  if (status != 0)
    Error("TcpNoDelay: setsockopt failed", (long) status);
}

void ShutdownAll()
/* 
   close all sockets discarding any pending data in either direction.
*/
{
   int i;

   for (i=0; i<NNODES_(); i++)
      if (SR_proc_info[i].sock >= 0) {
         (void) shutdown(SR_proc_info[i].sock, 2);
         (void) close(SR_proc_info[i].sock);
      }
}

int ReadFromSocket(sock, buf, lenbuf)
     int sock;
     char *buf;
     long lenbuf;
/*
   Read from the socket until we get all we want.
*/
{
   int nread, status;

   status = lenbuf;
   while (lenbuf > 0) {
again:
     if ( (nread = recv(sock, buf, (int) lenbuf, 0)) < 0) {
       if (errno == EINTR)
         goto again;
       else {
         (void) fprintf(stderr,"sock=%d, pid=%ld, nread=%d, len=%ld\n",
                               sock, NODEID_(), nread, lenbuf);
         (void) fflush(stderr);
         status = -1;
         break;
       }
     }
     buf += nread;
     lenbuf -= nread;
   }
   
   return status;
}

int WriteToSocket(sock, buf, lenbuf)
     int sock;
     char *buf;
     long lenbuf;
/*
  Write to the socket in packets of PACKET_SIZE bytes
*/
{
  int status = lenbuf;
  int nsent, len;
  
  while (lenbuf > 0) {
    
    len = (lenbuf > PACKET_SIZE) ? PACKET_SIZE : lenbuf;
    nsent = send(sock, buf, (int) len, 0);
    
    if (nsent < 0) { /* This is bad news */
      (void) fprintf(stderr,"sock=%d, pid=%ld, nsent=%d, len=%ld\n",
		     sock, NODEID_(), nsent, lenbuf);
      (void) fflush(stderr);
      status = -1; break;
    }

    buf += nsent;
    lenbuf -= nsent;
  }
  
  return status;
}

void CreateSocketAndBind(sock, port)
     int *sock;
     int *port;
/*
  Create a socket, bind it to a wildcard internet name and return
  the info so that its port number may be advertised
*/
{
  unsigned int length;
  struct sockaddr_in server;
  int size = SR_SOCK_BUF_SIZE;
  int on = 1;
#if defined(LINUX) && defined(__powerpc__)
  int dupsock;
#endif

  length = sizeof (struct sockaddr_in);

  /* Create socket */

  if ( (*sock = socket(AF_INET, SOCK_STREAM, 0)) < 0)
    Error("CreateSocketAndBind: socket creation failed", (long) *sock);

#if defined(LINUX) && defined(__powerpc__)
  if(*sock==0)
    dupsock = dup(*sock);
  *sock = dupsock;
#endif

  if(setsockopt(*sock, SOL_SOCKET, SO_REUSEADDR, 
		(char *) &on, sizeof on) == -1)
	Error("CreateSocketAndBind: error from setsockopt", (long) -1);

  /* Increase size of socket buffers to improve long message
     performance and increase size of message that goes asynchronously */

  if(setsockopt(*sock, SOL_SOCKET, SO_RCVBUF, (char *) &size, sizeof size))
    Error("CreateSocketAndBind: error setting SO_RCVBUF", (long) size);
  if(setsockopt(*sock, SOL_SOCKET, SO_SNDBUF, (char *) &size, sizeof size))
    Error("CreateSocketAndBind: error setting SO_SNDBUF", (long) size);

#ifndef ARDENT
  TcpNoDelay(*sock);
#endif

  /* Name socket with wildcards */

  server.sin_family = AF_INET;
  server.sin_addr.s_addr = INADDR_ANY;
  server.sin_port = 0;
  if (bind(*sock, (struct sockaddr *) &server, length) < 0)
    Error("CreateSocketAndBind: bind failed", (long) 0);

  /* Find out port number etc. */

  if (getsockname(*sock, (struct sockaddr *) &server, &length) < 0)
    Error("CreateSocketAndBind: getsockname failed", (long) 0);

  *port = ntohs(server.sin_port);

}

void ListenOnSock(sock)
  int sock;
/*
  Listen for a connection on the specified socket
  which was created with CreateSocketAndBind
*/
{
againlist:
  if (listen(sock, 1) < 0) {
    if (errno == EINTR)
      goto againlist;
    else
      Error("ListenAndAccept: listen failed", (long) 0);
  }

  if (DEBUG_) {
    (void) printf("process %ld out of listen on socket %d\n",NODEID_(),sock);
    (void) fflush(stdout);
  }
}

int AcceptConnection(sock)
  int sock;
/*
  Accept a connection on the specified socket
  which was created with CreateSocketAndBind and
  listen has been called.
*/
{
  fd_set ready;
  struct timeval timelimit;
  int msgsock, nready;
  int size = SR_SOCK_BUF_SIZE;
  
  /* Use select to wait for someone to try and establish a connection
     so that we can add a short timeout to avoid hangs */

againsel:
  FD_ZERO(&ready);
  FD_SET(sock, &ready);

  timelimit.tv_sec = TIMEOUT_ACCEPT;
  timelimit.tv_usec = 0;
  nready = select(sock+1, &ready, (fd_set *) NULL, (fd_set *) NULL,
		  &timelimit);
  if ( (nready <= 0) && (errno == EINTR) )
    goto againsel;
  else if (nready < 0)
    Error("ListenAndAccept: error from select",  (long) nready);
  else if (nready == 0)
    Error("ListenAndAccept: timeout waiting for connection", 
          (long) nready);

  if (!FD_ISSET(sock, &ready))
    Error("ListenAndAccept: out of select but not ready!", (long) nready);

againacc:
  msgsock = accept(sock, (struct sockaddr *) NULL, (unsigned int *) NULL);
  if (msgsock == -1) {
    if (errno == EINTR)
      goto againacc;
    else
      Error("ListenAndAccept: accept failed", (long) msgsock);
  }

  if (DEBUG_) {
    (void) printf("process %ld out of accept on socket %d\n",
		  NODEID_(),msgsock);
    (void) fflush(stdout);
  }

  /* Increase size of socket buffers to improve long message
     performance and increase size of message that goes asynchronously */

  if(setsockopt(msgsock, SOL_SOCKET, SO_RCVBUF, (char *) &size, sizeof size))
    Error("ListenAndAccept: error setting SO_RCVBUF", (long) size);
  if(setsockopt(msgsock, SOL_SOCKET, SO_SNDBUF, (char *) &size, sizeof size))
    Error("ListenAndAccept: error setting SO_SNDBUF", (long) size);

#ifndef ARDENT
  TcpNoDelay(msgsock);
#endif

  (void) close(sock); /* will not be needing this again */
  return msgsock;
}

int ListenAndAccept(sock)
  int sock;
/*
  Listen and accept a connection on the specified socket
  which was created with CreateSocketAndBind
*/
{
  fd_set ready;
  struct timeval timelimit;
  int msgsock, nready;
  int size = SR_SOCK_BUF_SIZE;
  
againlist:
  if (listen(sock, 1) < 0) {
    if (errno == EINTR)
      goto againlist;
    else
      Error("ListenAndAccept: listen failed", (long) 0);
  }

  if (DEBUG_) {
    (void) printf("process %ld out of listen on socket %d\n",NODEID_(),sock);
    (void) fflush(stdout);
  }

  /* Use select to wait for someone to try and establish a connection
     so that we can add a short timeout to avoid hangs */

againsel:
  FD_ZERO(&ready);
  FD_SET(sock, &ready);

  timelimit.tv_sec = TIMEOUT_ACCEPT;
  timelimit.tv_usec = 0;
  nready = select(sock+1, &ready, (fd_set *) NULL, (fd_set *) NULL,
		  &timelimit);
  if ( (nready <= 0) && (errno == EINTR) )
    goto againsel;
  else if (nready < 0)
    Error("ListenAndAccept: error from select",  (long) nready);
  else if (nready == 0)
    Error("ListenAndAccept: timeout waiting for connection", 
          (long) nready);

  if (!FD_ISSET(sock, &ready))
    Error("ListenAndAccept: out of select but not ready!", (long) nready);

againacc:
  msgsock = accept(sock, (struct sockaddr *) NULL, (unsigned int *) NULL);
  if (msgsock == -1) {
    if (errno == EINTR)
      goto againacc;
    else
      Error("ListenAndAccept: accept failed", (long) msgsock);
  }

  if (DEBUG_) {
    (void) printf("process %ld out of accept on socket %d\n",
		  NODEID_(),msgsock);
    (void) fflush(stdout);
  }

  /* Increase size of socket buffers to improve long message
     performance and increase size of message that goes asynchronously */

  if(setsockopt(msgsock, SOL_SOCKET, SO_RCVBUF, (char *) &size, sizeof size))
    Error("ListenAndAccept: error setting SO_RCVBUF", (long) size);
  if(setsockopt(msgsock, SOL_SOCKET, SO_SNDBUF, (char *) &size, sizeof size))
    Error("ListenAndAccept: error setting SO_SNDBUF", (long) size);

#ifndef ARDENT
  TcpNoDelay(msgsock);
#endif

  (void) close(sock); /* will not be needing this again */
  return msgsock;
}

int CreateSocketAndConnect(hostname, cport)
     char *hostname;
     char *cport;
/*
  Return the file descriptor of the socket which connects me to the
  remote process on hostname at port in string cport

  hostname = hostname of the remote process
  cport    = asci string containing port number of remote socket
*/
{
  int sock, status;
  struct sockaddr_in server;
  struct hostent *hp;
  int on = 1;
  int size = SR_SOCK_BUF_SIZE;
#ifndef SGI
  struct hostent *gethostbyname();
#endif
#if defined(LINUX) && defined(__powerpc__)
  int dupsock;
#endif

  /* Create socket */

  if ( (sock = socket(AF_INET, SOCK_STREAM, 0)) < 0 ) {
    (void) fprintf(stderr,"trying to connect to host=%s, port=%s\n",
                   hostname, cport);
    Error("CreateSocketAndConnect: socket failed",  (long) sock);
  }

#if defined(LINUX) && defined(__powerpc__)
  if(sock==0)
    dupsock = dup(sock);
  sock = dupsock;
#endif

  if (setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, 
		 (char *) &on, sizeof on) == -1)
	Error("CreateSocketAndConnect: error setting REUSEADDR", (long) -1);

  /* Increase size of socket buffers to improve long message
     performance and increase size of message that goes asynchronously */

  if(setsockopt(sock, SOL_SOCKET, SO_RCVBUF, (char *) &size, sizeof size))
    Error("CreateSocketAndConnect: error setting SO_RCVBUF", (long) size);
  if(setsockopt(sock, SOL_SOCKET, SO_SNDBUF, (char *) &size, sizeof size))
    Error("CreateSocketAndConnect: error setting SO_SNDBUF", (long) size);

#ifndef ARDENT
  TcpNoDelay(sock);
#endif

  /* Connect socket */

  server.sin_family = AF_INET;
  hp = gethostbyname(hostname);
  if (hp == 0) {
    (void) fprintf(stderr,"trying to connect to host=%s, port=%s\n",
                   hostname, cport);
    Error("CreateSocketAndConnect: gethostbyname failed", (long) 0);
  }

  bcopy((char *) hp->h_addr, (char *) &server.sin_addr, hp->h_length);
  server.sin_port = htons((ushort) atoi(cport));

againcon:
  if ((status = 
     connect(sock, (struct sockaddr *) &server, sizeof server)) < 0) {
    if (errno == EINTR)
      goto againcon;
    else {
      (void) fprintf(stderr,"trying to connect to host=%s, port=%s\n",
                   hostname, cport);
      Error("CreateSocketAndConnect: connect failed", (long) status);
    }
  }
  
  return sock;
}
