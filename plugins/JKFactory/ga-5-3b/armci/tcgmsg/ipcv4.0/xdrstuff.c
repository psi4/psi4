#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/xdrstuff.c,v 1.5 2004-04-01 02:04:57 manoj Exp $ */

#ifdef GOTXDR

#include <rpc/types.h>
#include <rpc/xdr.h>

#ifdef CRAY
extern bool_t xdr_char();
extern char *malloc();
#define NULL 0
#endif

#ifdef SEQUENT
static bool_t xdr_char();  /* below from sun distribution tape */
#define NULL 0
#endif

#ifdef HPUX
#define NULL 0
#endif

#if defined(ULTRIX) || defined(SGI) || defined(NEXT) || defined(HPUX)|| \
    defined(AIX) || defined(KSR) || defined(DECOSF)
extern void *malloc();
#else
extern char *malloc();
#endif

#define XDR_BUF_LEN 4096        /* Size of XDR buffer in bytes */
#define XDR_DOUBLE_LEN 8        /* Size of XDR double in bytes */
#define XDR_LONG_LEN 4          /* Size of XDR long in bytes */
#define XDR_CHAR_LEN 4          /* Size of XDR char in bytes */

static char *xdrbuf_decode;
static char *xdrbuf_encode;
static XDR xdr_decode;
static XDR xdr_encode;

static int xdr_buf_allocated = 0; /* =1 if buffers allocated, 0 otherwise */

extern void Error();

void CreateXdrBuf()
/*
  Call at start to allocate the XDR buffers
*/
{
  if (!xdr_buf_allocated) {

    /* Malloc the buffer space */
  
    if ( (xdrbuf_decode = malloc((unsigned) XDR_BUF_LEN)) == (char *) NULL)
      Error("CreateXdrBuf: malloc of xdrbuf_decode failed",
            (long) XDR_BUF_LEN);

    if ( (xdrbuf_encode = malloc((unsigned) XDR_BUF_LEN)) == (char *) NULL)
      Error("CreateXdrBuf: malloc of xdrbuf_encode failed",
            (long) XDR_BUF_LEN);

    /* Associate the xdr memory streams with the buffers */

    xdrmem_create(&xdr_decode, xdrbuf_decode, XDR_BUF_LEN, XDR_DECODE);

    xdrmem_create(&xdr_encode, xdrbuf_encode, XDR_BUF_LEN, XDR_ENCODE);

    xdr_buf_allocated = 1;
  }
}

void DestroyXdrBuf()
/*
  Call to free the xdr buffers
*/
{
  if (xdr_buf_allocated) {

    /* Destroy the buffers and free the space */

    xdr_destroy(&xdr_encode);
    xdr_destroy(&xdr_decode);
    (void) free(xdrbuf_encode);
    (void) free(xdrbuf_decode);

    xdr_buf_allocated = 0;
  }
}

int WriteXdrDouble(sock, x, n_double)
    int sock;
    double *x;
    long n_double;
/*
  Write double x[n_double] to the socket translating to XDR representation.

  Returned is the number of bytes written to the socket.

  All errors are treated as fatal.
*/
{
  int nd_per_buf = (XDR_BUF_LEN-4)/XDR_DOUBLE_LEN;
                                              /* No. of XDR doubles per buf */
  int status, nb=0;
  u_int len;
  long lenb;

  if (!xdr_buf_allocated)
    CreateXdrBuf();

  /* Loop thru buffer loads */
    
  while (n_double > 0) {

    len = (n_double > nd_per_buf) ? nd_per_buf : n_double;

    /* Position the xdr buffer to the beginning */

    if (!xdr_setpos(&xdr_encode, (u_int) 0))
      Error("WriteXdrDouble: xdr_setpos failed", (long) -1);

    /* Translate the buffer and then write it to the socket */

    if (!xdr_array(&xdr_encode, (char **) &x, &len, (u_int) XDR_BUF_LEN,
                   (u_int) sizeof(double), xdr_double))
       Error("WriteXdrDouble: xdr_array failed", (long) -1);

    lenb = xdr_getpos(&xdr_encode);

    if ((status = WriteToSocket(sock, xdrbuf_encode, lenb)) != lenb)
      Error("WriteXdrDouble: WriteToSocket failed", (long) status);

    nb += lenb;
    n_double -= len;
    x += len;
  }

  return nb;
}

int ReadXdrDouble(sock, x, n_double)
    int sock;
    double *x;
    long n_double;
/*
  Read double x[n_double] from the socket translating from XDR representation.

  Returned is the number of bytes read from the socket.

  All errors are treated as fatal.
*/
{
  int nd_per_buf = (XDR_BUF_LEN-4)/XDR_DOUBLE_LEN; 
                                              /* No. of XDR doubles per buf */
  int status, nb=0;
  u_int len;
  long lenb;

  if (!xdr_buf_allocated)
    CreateXdrBuf();

  /* Loop thru buffer loads */
    
  while (n_double > 0) {

    len = (n_double > nd_per_buf) ? nd_per_buf : n_double;
    lenb = 4 + len * XDR_DOUBLE_LEN;

    /* Position the xdr buffer to the beginning */

    if (!xdr_setpos(&xdr_decode, (u_int) 0))
      Error("ReadXdrDouble: xdr_setpos failed", (long) -1);

    /* Read from the socket and then translate the buffer */

    if ((status = ReadFromSocket(sock, xdrbuf_decode, lenb)) != lenb)
      Error("ReadXdrDouble: ReadFromSocket failed", (long) status);

    if (!xdr_array(&xdr_decode, (char **) &x, &len, (u_int) XDR_BUF_LEN, 
                   (u_int) sizeof(double), xdr_double))
       Error("ReadXdrDouble: xdr_array failed", (long) -1);

    nb += lenb;
    n_double -= len;
    x += len;
  }

  return nb;
}

int WriteXdrLong(sock, x, n_long)
    int sock;
    long *x;
    long n_long;
/*
  Write long x[n_long] to the socket translating to XDR representation.

  Returned is the number of bytes written to the socket.

  All errors are treated as fatal.
*/
{
  int nd_per_buf = (XDR_BUF_LEN-4)/XDR_LONG_LEN;
                                              /* No. of XDR longs per buf */
  int status, nb=0;
  u_int len;
  long lenb;

  if (!xdr_buf_allocated)
    CreateXdrBuf();

  /* Loop thru buffer loads */
    
  while (n_long > 0) {

    len = (n_long > nd_per_buf) ? nd_per_buf : n_long;

    /* Position the xdr buffer to the beginning */

    if (!xdr_setpos(&xdr_encode, (u_int) 0))
      Error("WriteXdrLong: xdr_setpos failed", (long) -1);

    /* Translate the buffer and then write it to the socket */

    if (!xdr_array(&xdr_encode, (char **) &x, &len, (u_int) XDR_BUF_LEN,
                   (u_int) sizeof(long), xdr_long))
       Error("WriteXdrLong: xdr_array failed", (long) -1);

    lenb = xdr_getpos(&xdr_encode);

    if ((status = WriteToSocket(sock, xdrbuf_encode, lenb)) != lenb)
      Error("WriteXdrLong: WriteToSocket failed", (long) status);

    nb += lenb;
    n_long -= len;
    x += len;
  }

  return nb;
}

int ReadXdrLong(sock, x, n_long)
    int sock;
    long *x;
    long n_long;
/*
  Read long x[n_long] from the socket translating from XDR representation.

  Returned is the number of bytes read from the socket.

  All errors are treated as fatal.
*/
{
  int nd_per_buf = (XDR_BUF_LEN-4)/XDR_LONG_LEN; 
                                              /* No. of XDR longs per buf */
  int status, nb=0;
  u_int len;
  long lenb;

  if (!xdr_buf_allocated)
    CreateXdrBuf();

  /* Loop thru buffer loads */
    
  while (n_long > 0) {

    len = (n_long > nd_per_buf) ? nd_per_buf : n_long;
    lenb = 4 + len * XDR_LONG_LEN;

    /* Position the xdr buffer to the beginning */

    if (!xdr_setpos(&xdr_decode, (u_int) 0))
      Error("ReadXdrLong: xdr_setpos failed", (long) -1);

    /* Read from the socket and then translate the buffer */

    if ((status = ReadFromSocket(sock, xdrbuf_decode, lenb)) != lenb)
      Error("ReadXdrLong: ReadFromSocket failed", (long) status);

    if (!xdr_array(&xdr_decode, (char **) &x, &len, (u_int) XDR_BUF_LEN, 
                   (u_int) sizeof(long), xdr_long))
       Error("ReadXdrLong: xdr_array failed", (long) -1);

    nb += lenb;
    n_long -= len;
    x += len;
  }

  return nb;
}

int WriteXdrChar(sock, x, n_char)
    int sock;
    char *x;
    long n_char;
/*
  Write char x[n_char] to the socket translating to XDR representation.

  Returned is the number of bytes written to the socket.

  All errors are treated as fatal.
*/
{
  int nc_per_buf = (XDR_BUF_LEN-4)/XDR_CHAR_LEN;
                                              /* No. of XDR chars per buf */
  int status, nb=0;
  u_int len;
  long lenb;

  if (!xdr_buf_allocated)
    CreateXdrBuf();

  /* Loop thru buffer loads */
    
  while (n_char > 0) {

    len = (n_char > nc_per_buf) ? nc_per_buf : n_char;

    /* Position the xdr buffer to the beginning */

    if (!xdr_setpos(&xdr_encode, (u_int) 0))
      Error("WriteXdrChar: xdr_setpos failed", (long) -1);

    /* Translate the buffer and then write it to the socket */

    if (!xdr_array(&xdr_encode, (char **) &x, &len, (u_int) XDR_BUF_LEN,
                   (u_int) sizeof(char), xdr_char))
       Error("WriteXdrChar: xdr_array failed", (long) -1);

    lenb = xdr_getpos(&xdr_encode);

    if ((status = WriteToSocket(sock, xdrbuf_encode, lenb)) != lenb)
      Error("WriteXdrChar: WriteToSocket failed", (long) status);

    nb += lenb;
    n_char -= len;
    x += len;
  }

  return nb;
}

int ReadXdrChar(sock, x, n_char)
    int sock;
    char *x;
    long n_char;
/*
  Read char x[n_char] from the socket translating from XDR representation.

  Returned is the number of bytes read from the socket.

  All errors are treated as fatal.
*/
{
  int nc_per_buf = (XDR_BUF_LEN-4)/XDR_CHAR_LEN; 
                                              /* No. of XDR chars per buf */
  int status, nb=0;
  u_int len;
  long lenb;

  if (!xdr_buf_allocated)
    CreateXdrBuf();

  /* Loop thru buffer loads */
    
  while (n_char > 0) {

    len = (n_char > nc_per_buf) ? nc_per_buf : n_char;
    lenb = 4 + len * XDR_CHAR_LEN;

    /* Position the xdr buffer to the beginning */

    if (!xdr_setpos(&xdr_decode, (u_int) 0))
      Error("ReadXdrChar: xdr_setpos failed", (long) -1);

    /* Read from the socket and then translate the buffer */

    if ((status = ReadFromSocket(sock, xdrbuf_decode, lenb)) != lenb)
      Error("ReadXdrChar: ReadFromSocket failed", (long) status);

    if (!xdr_array(&xdr_decode, (char **) &x, &len, (u_int) XDR_BUF_LEN, 
                   (u_int) sizeof(char), xdr_char))
       Error("ReadXdrChar: xdr_array failed", (long) -1);

    nb += lenb;
    n_char -= len;
    x += len;
  }

  return nb;
}

#ifdef SEQUENT
/*
 * XDR a char
 */
static bool_t xdr_char(xdrs, cp)
        XDR *xdrs;
        char *cp;
{
        int i;

        i = (*cp);
        if (!xdr_int(xdrs, &i)) {
                return (FALSE);
        }
        *cp = i;
        return (TRUE);
}
#endif

#else
/* dummy function to make this source file legitimate */
#include <stdio.h>
#include <stdlib.h>
void _dummy_ZefP_() {printf("XDR:Illegal function call\n"); exit(1);}
#endif
