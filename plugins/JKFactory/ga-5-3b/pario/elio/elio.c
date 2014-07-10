#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 **********************************************************************\
 ELementary I/O (ELIO) disk operations for parallel I/O libraries   
 Authors: Jarek Nieplocha (PNNL) and Jace Mogill (ANL)
 *
 * DISCLAIMER
 *
 * This material was prepared as an account of work sponsored by an
 * agency of the United States Government.  Neither the United States
 * Government nor the United States Department of Energy, nor Battelle,
 * nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
 * COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
 * SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
 * INFRINGE PRIVATELY OWNED RIGHTS.
 *
 * ACKNOWLEDGMENT
 *
 * This software and its documentation were produced with United States
 * Government support under Contract Number DE-AC06-76RLO-1830 awarded by
 * the United States Department of Energy.  The United States Government
 * retains a paid-up non-exclusive, irrevocable worldwide license to
 * reproduce, prepare derivative works, perform publicly and display
 * publicly by or for the US Government, including the right to
 * distribute to other US Government contractors.
 */
#ifdef USE_LUSTRE
#include <lustre/lustre_user.h> /* for O_LOV_DELAY_CREATE, LL_IOC_LOV_SETSTRIPE */
#include <linux/lustre_idl.h> /* for struct lov_mds_md, LOV_MAGIC */
#include <sys/ioctl.h> /* for ioctl */
#endif

#include "eliop.h"

#include "../sf/coms.h"

#if defined(CRAY) && defined(__crayx1)
#undef CRAY
#endif

#if  defined(AIX) || defined(DECOSF) || defined(SGITFP) || defined(SGI64) || defined(SGI_N32) || defined(CRAY) || defined(LINUXAIO)
     /* systems with Asynchronous I/O */
#else
#    ifndef NOAIO
#      define NOAIO
#    endif
#endif

/****************** Internal Constants and Parameters **********************/

#define  MAX_AIO_REQ  4
#define  NULL_AIO    -123456
#define  FOPEN_MODE 0644
#define  MAX_ATTEMPTS 10


#ifndef NOAIO
#   define AIO 1
#endif


#ifdef FFIO
#  define WRITE  ffwrite
#  define WRITEA ffwritea
#  define READ   ffread
#  define READA  ffreada
#  define CLOSE  ffclose
#  define SEEK   ffseek
#  define OPEN   ffopens
#  define DEFARG FULL
#else
#  define WRITE  write
#  define WRITEA writea
#  define READ   read
#  define READA  reada
#  define CLOSE  close
#  define SEEK   lseek
#  define OPEN   open
#  define DEFARG 0
#endif


#ifdef WIN32
#define ELIO_FSYNC _commit
#else
#include <unistd.h>
#define ELIO_FSYNC fsync
#endif

/* structure to emulate control block in Posix AIO */
#if defined (CRAY)
#   if defined(FFIO)
       typedef struct { struct ffsw stat; int filedes; }io_status_t;
#   else 
#      include <sys/iosw.h>
       typedef struct { struct iosw stat; int filedes; }io_status_t;
#   endif
    io_status_t cb_fout[MAX_AIO_REQ];
    io_status_t *cb_fout_arr[MAX_AIO_REQ];

#elif defined(AIO)
#   include <aio.h>
#   if defined(AIX)
#      define INPROGRESS EINPROG
#   else 
#      define INPROGRESS EINPROGRESS
#   endif
    struct aiocb          cb_fout[MAX_AIO_REQ];
#ifndef AIX
    const
#endif
           struct aiocb   *cb_fout_arr[MAX_AIO_REQ];
#endif

#ifndef INPROGRESS
#   define INPROGRESS 1
#endif

static long           aio_req[MAX_AIO_REQ]; /* array for AIO requests */
static int            first_elio_init = 1;  /* intialization status */
int                   _elio_Errors_Fatal=0; /* sets mode of handling errors */


/****************************** Internal Macros *****************************/
#if defined(AIO)
#  define AIO_LOOKUP(aio_i) {\
      aio_i = 0;\
      while(aio_req[aio_i] != NULL_AIO && aio_i < MAX_AIO_REQ) aio_i++;\
}
#else
#  define AIO_LOOKUP(aio_i) aio_i = MAX_AIO_REQ
#endif

#define SYNC_EMULATE(op) *req_id = ELIO_DONE; \
  if((stat= elio_ ## op (fd, offset, buf, bytes)) != bytes ){ \
       ELIO_ERROR(stat,0);  \
  }else \
       stat    = 0; 

#ifndef MIN 
#define PARIO_MIN(a,b) (((a) <= (b)) ? (a) : (b))
#endif

/* 
 * Offsets bigger than ABSURDLY_LARGE generate a SEEKFAIL.
 * The maximum no. of extents permitted for a file is MAX_EXTENT.
 */

#if defined(_LARGE_FILES) || defined(_LARGEFILE_SOURCE) || defined(_LARGEFILE64_SOURCE) || _FILE_OFFSET_BITS+0 == 64 || SIZEOF_VOIDP == 8
#   define LARGE_FILES
#endif

#define MAX_EXTENT 127
#ifdef LARGE_FILES
#define ABSURDLY_LARGE 1e14
#else
#define ABSURDLY_LARGE (MAX_EXTENT*2147483648.0)
#endif

/*****************************************************************************/

static Off_t elio_max_file_size(Fd_t fd)
     /* 
      * Return the maximum size permitted for this PHYSICAL file.
      * Presently not file dependent.
      */
{
#ifdef LARGE_FILES
  return ABSURDLY_LARGE;
#else
  return (2047.0*1024.0*1024.0);        /* 2 GB - 1 MB */  
#endif
}

static Fd_t elio_get_next_extent(Fd_t fd)
     /*
      * Return a pointer to the file descriptor that forms
      * the next extent of this file.  If the extension file
      * does not exist then it is opened.  If the open fails
      * then the usual error condition of elio_open is returned.
      */
{
  Fd_t next_fd = (Fd_t) fd->next;
  if (!next_fd) {
    /* Eventually need to replace this with user controllable naming 
     * and combine with similar logic in delete routine.
     */
    char fname[ELIO_FILENAME_MAX];
    int len;
    if (fd->extent >= MAX_EXTENT)
      return 0;
    strcpy(fname, fd->name);
    len = strlen(fname);
    if (fd->extent) len -= 4;
    sprintf(fname+len,"x%3.3d",fd->extent+1);
    /*printf("Opening extent %d with name '%s'\n",fd->extent+1,fname);*/
    if ((next_fd = elio_open(fname, fd->type, fd->mode))) {
      next_fd->extent = fd->extent + 1;
      fd->next = (struct fd_struct *) next_fd;
    }
  }
  return next_fd;
}

void elio_errors_fatal(int onoff)
{
    _elio_Errors_Fatal = onoff;
}
 

/*\ Blocking Write 
 *    - returns number of bytes written or error code (<0) if failed
\*/
Size_t elio_write(Fd_t fd, Off_t  doffset, const void* buf, Size_t bytes)
{
  off_t offset;
  Size_t stat, bytes_to_write = bytes;
  Size_t nextbytes;

  if (doffset >= ABSURDLY_LARGE) 
    ELIO_ERROR(SEEKFAIL,0);

  /* Follow the linked list of extents down until we hit the file
     that contains the offset */
  if (doffset >= elio_max_file_size(fd)) {
    Fd_t next_fd = elio_get_next_extent(fd);
    if (!next_fd) ELIO_ERROR(OPENFAIL,0);
    doffset -= elio_max_file_size(fd);
    return elio_write(next_fd, doffset, buf, bytes);
  }

  /* Figure out if the write continues onto the next extent */
  offset = (off_t) doffset;
  nextbytes = 0;
  if ((doffset+bytes_to_write) >= elio_max_file_size(fd)) {
    nextbytes = bytes_to_write;
    bytes_to_write = (Size_t) (elio_max_file_size(fd)-doffset);
    nextbytes -= bytes_to_write;
  }
  /*printf("TRYING TO WRITE AT doffset=%f offset=%lu bw=%lu nb=%lu\n", doffset, offset,
    bytes_to_write, nextbytes);*/

  /* Write to this extent */

#ifdef PABLO
  int pablo_code = PABLO_elio_write;
  PABLO_start( pablo_code );
#endif
  
  if(offset != SEEK(fd->fd,offset,SEEK_SET)) ELIO_ERROR(SEEKFAIL,0);
  
  while (bytes_to_write) {
    stat = WRITE(fd->fd, buf, bytes_to_write);
    if ((stat == -1) && ((errno == EINTR) || (errno == EAGAIN))) {
      ; /* interrupted write should be restarted */
    } else if (stat > 0) {
      bytes_to_write -= stat;
      buf = stat + (char*)buf; /*advance pointer by # bytes written*/
    } else {
      ELIO_ERROR(WRITFAIL, stat);
    }
  }

  /* Only get here if all has gone OK */
  
#ifdef PABLO
  PABLO_end(pablo_code);
#endif

  /* Write to next extent(s) ... relies on incrementing of buf */
  if (nextbytes) {
    Fd_t next_fd = elio_get_next_extent(fd);
    if (!next_fd) ELIO_ERROR(OPENFAIL,0);
    stat = elio_write(next_fd, (Off_t) 0, buf, nextbytes);
    if (stat != nextbytes)
      ELIO_ERROR(WRITFAIL, stat);
  }
  
  return bytes;
}

int elio_set_cb(Fd_t fd, Off_t doffset, int reqn, void *buf, Size_t bytes)
{
#if defined(AIO)
    off_t offset = (off_t) doffset;
#   if defined(CRAY)
       if(offset != SEEK(fd->fd, offset, SEEK_SET))return (SEEKFAIL);
       cb_fout_arr[reqn] = cb_fout+reqn;
       cb_fout[reqn].filedes    = fd->fd;
#   else
       cb_fout[reqn].aio_offset = offset;
       cb_fout_arr[reqn] = cb_fout+reqn;
         cb_fout[reqn].aio_buf    = buf;
         cb_fout[reqn].aio_nbytes = bytes;
#        if defined(AIX)
           cb_fout[reqn].aio_whence = SEEK_SET;
#        else
           cb_fout[reqn].aio_sigevent.sigev_notify = SIGEV_NONE;
           cb_fout[reqn].aio_fildes    = fd->fd;
#        endif
#   endif
#endif
    return ELIO_OK;
}


/*\ Asynchronous Write: returns 0 if succeded or err code if failed
\*/
int elio_awrite(Fd_t fd, Off_t doffset, const void* buf, Size_t bytes, io_request_t * req_id)
{
  off_t offset;
  Size_t stat;
#ifdef AIO
  int    aio_i;
#endif

  if (doffset >= ABSURDLY_LARGE) 
    ELIO_ERROR(SEEKFAIL,0);

  /* Follow the linked list of extents down until we hit the file
     that contains the offset */
  if (doffset >= elio_max_file_size(fd)) {
    Fd_t next_fd = elio_get_next_extent(fd);
    if (!next_fd) ELIO_ERROR(OPENFAIL,0);
    doffset -= elio_max_file_size(fd);
    return elio_awrite(next_fd, doffset, buf, bytes, req_id);
  }

  /* Figure out if the write continues onto the next extent 
   * ... if so then force the entire request to be done synchronously
   * so that we don't have to manage multiple async requests */

  if ((doffset+((Off_t) bytes)) >= elio_max_file_size(fd)) {
    *req_id = ELIO_DONE;
    if (elio_write(fd, doffset, buf, bytes) != bytes)
      return -1;
    else
      return 0;
  }

  offset = (off_t) doffset;

#ifdef PABLO
  int pablo_code = PABLO_elio_awrite;
  PABLO_start( pablo_code );
#endif

  *req_id = ELIO_DONE;

#ifdef AIO
   AIO_LOOKUP(aio_i);

   /* blocking io when request table is full */
   if(aio_i >= MAX_AIO_REQ){
#     if defined(DEBUG) && defined(ASYNC)
         fprintf(stderr, "elio_awrite: Warning- asynch overflow\n");
#     endif
      SYNC_EMULATE(write);
   } else {
      int rc;
      *req_id = (io_request_t) aio_i;
      if((rc=elio_set_cb(fd, offset, aio_i, (void*) buf, bytes)))
                                                 ELIO_ERROR(rc,0);

#    if defined(CRAY)
       rc = WRITEA(fd->fd, (char*)buf, bytes, &cb_fout[aio_i].stat, DEFARG);
       stat = (rc < 0)? -1 : 0; 
#    elif defined(AIX) 
#       if !defined(AIX52) && !defined(_AIO_AIX_SOURCE)
       stat = aio_write(fd->fd, cb_fout + aio_i);
#       endif
#    else
       stat = aio_write(cb_fout+aio_i);
#    endif
     aio_req[aio_i] = *req_id;
  }

#else
      /* call blocking write when AIO not available */
      SYNC_EMULATE(write);
#endif

  if(stat ==-1) ELIO_ERROR(AWRITFAIL, 0);

#ifdef PABLO
  PABLO_end(pablo_code);
#endif

  return((int)stat);
}


/*\ Truncate the file at the specified length.
\*/
int elio_truncate(Fd_t fd, Off_t dlength)
{
  off_t length = (off_t) dlength;
#ifdef WIN32
#   define ftruncate _chsize 
#endif

#ifdef PABLO
    int pablo_code = PABLO_elio_truncate;
    PABLO_start( pablo_code );
#endif
    if(dlength >= elio_max_file_size(fd)){
      Fd_t next_fd = elio_get_next_extent(fd);
      dlength -= elio_max_file_size(fd);
#       if defined(DEBUG)
      printf(stderr," calling ftruncate with length = %f \n", dlength);
#endif
      return elio_truncate(next_fd, dlength);
    }
    (void) SEEK(fd->fd, 0L, SEEK_SET);
    if (ftruncate(fd->fd, length))
    return TRUNFAIL;
    else {
    return ELIO_OK;
    }
#ifdef PABLO
    PABLO_end(pablo_code);
#endif
}


/*\ Return in length the length of the file
\*/
int elio_length(Fd_t fd, Off_t *dlength)
{
  off_t length;
  int status;

  /* Add up the lengths of any extents */
  if (fd->next) {
    status = elio_length((Fd_t) fd->next, dlength);
    *dlength += elio_max_file_size(fd);
    return status;
  }
  else {
#ifdef PABLO
    int pablo_code = PABLO_elio_length;
    PABLO_start( pablo_code );
#endif
    
    if ((length = SEEK(fd->fd, (off_t) 0, SEEK_END)) != -1)
      status = ELIO_OK;
    else
      status = SEEKFAIL;
    
#ifdef PABLO
    PABLO_end(pablo_code);
#endif
    
    *dlength = (Off_t) length;
    return status;
  }
}


/*\ Blocking Read 
 *      - returns number of bytes read or error code (<0) if failed
\*/
Size_t elio_read(Fd_t fd, Off_t doffset, void* buf, Size_t bytes)
{
off_t offset;
Size_t stat, bytes_to_read = bytes;
Size_t nextbytes;
int    attempt=0;

 if (doffset >= ABSURDLY_LARGE) 
    ELIO_ERROR(SEEKFAIL,0);

  /* Follow the linked list of extents down until we hit the file
     that contains the offset */
  if (doffset >= elio_max_file_size(fd)) {
    Fd_t next_fd = elio_get_next_extent(fd);
    if (!next_fd) ELIO_ERROR(OPENFAIL,0);
    doffset -= elio_max_file_size(fd);
    return elio_read(next_fd, doffset, buf, bytes);
  }

  /* Figure out if the read continues onto the next extent */
  offset = (off_t) doffset;
  nextbytes = 0;
  if ((doffset+bytes_to_read) >= elio_max_file_size(fd)) {
    nextbytes = bytes_to_read;
    bytes_to_read = (Size_t) (elio_max_file_size(fd)-doffset);
    nextbytes -= bytes_to_read;
  }


  /* Read from this physical file */

#ifdef PABLO
  int pablo_code = PABLO_elio_read;
  PABLO_start( pablo_code );
#endif

  if(offset != SEEK(fd->fd,offset,SEEK_SET)) ELIO_ERROR(SEEKFAIL,0);
  
  while (bytes_to_read) {
    stat = READ(fd->fd, buf, bytes_to_read);
    if(stat==0){
      ELIO_ERROR(EOFFAIL, stat);
    } else if ((stat == -1) && ((errno == EINTR) || (errno == EAGAIN))) {
      ; /* interrupted read should be restarted */
    } else if (stat > 0) {
      bytes_to_read -= stat;
      buf = stat + (char*)buf; /*advance pointer by # bytes read*/
    } else {
      ELIO_ERROR(READFAIL, stat);
    }
    attempt++;
  }
  
  /* Only get here if all went OK */
  
#ifdef PABLO
  PABLO_end(pablo_code);
#endif
  
  /* Read from next extent(s) ... relies on incrementing of buf */
  if (nextbytes) {
    Fd_t next_fd = elio_get_next_extent(fd);
    if (!next_fd) ELIO_ERROR(OPENFAIL,0);
    stat = elio_read(next_fd, (Off_t) 0, buf, nextbytes);
    if (stat != nextbytes)
      ELIO_ERROR(READFAIL, stat);
  }


  return bytes;
}



/*\ Asynchronous Read: returns 0 if succeded or -1 if failed
\*/
int elio_aread(Fd_t fd, Off_t doffset, void* buf, Size_t bytes, io_request_t * req_id)
{
  off_t offset = (off_t) doffset;
  Size_t stat;
#ifdef AIO
  int    aio_i;
#endif
#ifdef CRAY
  int rc;
#endif

  if (doffset >= ABSURDLY_LARGE) 
    ELIO_ERROR(SEEKFAIL,0);

  /* Follow the linked list of extents down until we hit the file
     that contains the offset */
  if (doffset >= elio_max_file_size(fd)) {
    Fd_t next_fd = elio_get_next_extent(fd);
    if (!next_fd) ELIO_ERROR(OPENFAIL,0);
    doffset -= elio_max_file_size(fd);
    return elio_aread(next_fd, doffset, buf, bytes, req_id);
  }

  /* Figure out if the read continues onto the next extent 
   * ... if so then force the entire request to be done synchronously
   * so that we don't have to manage multiple async requests */

  if ((doffset+((Off_t) bytes)) >= elio_max_file_size(fd)) {
    *req_id = ELIO_DONE;
    if (elio_read(fd, doffset, buf, bytes) != bytes)
      return -1;
    else
      return 0;
  }

  offset = (off_t) doffset;

#ifdef PABLO
  int pablo_code = PABLO_elio_aread;
  PABLO_start( pablo_code );
#endif

  *req_id = ELIO_DONE;

#ifdef AIO
    AIO_LOOKUP(aio_i);

    /* blocking io when request table is full */
    if(aio_i >= MAX_AIO_REQ){
#       if defined(DEBUG)
           fprintf(stderr, "elio_read: Warning- asynch overflow\n");
#       endif
        SYNC_EMULATE(read);

    } else {

       *req_id = (io_request_t) aio_i;
        if((stat=elio_set_cb(fd, offset, aio_i, (void*) buf, bytes)))
                                                 ELIO_ERROR((int)stat,0);
#       if defined(CRAY)
          rc = READA(fd->fd, buf, bytes, &cb_fout[aio_i].stat, DEFARG);
          stat = (rc < 0)? -1 : 0;
#       elif defined(AIX)
#if    !defined(AIX52) && !defined(_AIO_AIX_SOURCE)
          stat = aio_read(fd->fd, cb_fout+aio_i);
#endif
#       else
          stat = aio_read(cb_fout+aio_i);
#       endif
        aio_req[aio_i] = *req_id;
    }
#else

    /* call blocking write when AIO not available */
    SYNC_EMULATE(read);

#endif

    if(stat ==-1) ELIO_ERROR(AWRITFAIL, 0);

#ifdef PABLO
    PABLO_end(pablo_code);
#endif

    return((int)stat);
}


/*\ Wait for asynchronous I/O operation to complete. Invalidate id.
\*/
int elio_wait(io_request_t *req_id)
{
  int  aio_i=0;
  int  rc;

  rc=0; /* just to remove the compiler warning */
#ifdef PABLO
  int pablo_code = PABLO_elio_wait;
  PABLO_start( pablo_code );
#endif

  if(*req_id != ELIO_DONE ) { 

#    ifdef AIO
#      if defined(CRAY)

#        if defined(FFIO)
         {
            struct ffsw dumstat, *prdstat=&(cb_fout[*req_id].stat);
            fffcntl(cb_fout[*req_id].filedes, FC_RECALL, prdstat, &dumstat);
            if (FFSTAT(*prdstat) == FFERR) ELIO_ERROR(SUSPFAIL,0);
         }
#        else
         {
            struct iosw *statlist[1];
            statlist[0] = &(cb_fout[*req_id].stat);
            recall(cb_fout[*req_id].filedes, 1, statlist); 
         }
#        endif

#      elif defined(AIX) 
#         if    !defined(AIX52) && !defined(_AIO_AIX_SOURCE)
              do {    /* I/O can be interrupted on SP through rcvncall ! */
                   rc =(int)aio_suspend(1, cb_fout_arr+(int)*req_id);
              } while(rc == -1 && errno == EINTR); 
#         endif

#  else
      if((int)aio_suspend((const struct aiocb *const*)(cb_fout_arr+(int)*req_id), 1, NULL) != 0) rc =-1;
#  endif
      if(rc ==-1) ELIO_ERROR(SUSPFAIL,0);

#  if defined(DECOSF)
      /* on DEC aio_return is required to clean internal data structures */
      if(aio_return(cb_fout+(int)*req_id) == -1) ELIO_ERROR(RETUFAIL,0);
#  endif
#endif

      while(aio_req[aio_i] != *req_id && aio_i < MAX_AIO_REQ) aio_i++;
      if(aio_i >= MAX_AIO_REQ) ELIO_ERROR(HANDFAIL, aio_i);

      aio_req[aio_i] = NULL_AIO;
      *req_id = ELIO_DONE;
   }

#ifdef PABLO
   PABLO_end(pablo_code);
#endif

   return ELIO_OK;
}



/*\ Check if asynchronous I/O operation completed. If yes, invalidate id.
\*/
int elio_probe(io_request_t *req_id, int* status)
{
  int    errval=-1;
  int    aio_i = 0;
     
#ifdef PABLO
  int pablo_code = PABLO_elio_probe;
  PABLO_start( pablo_code );
#endif

  if(*req_id == ELIO_DONE){
      *status = ELIO_DONE;
  } else {
      
#ifdef AIO
#    if defined(CRAY)

#     if defined(FFIO)
      {
         struct ffsw dumstat, *prdstat=&(cb_fout[*req_id].stat);
         fffcntl(cb_fout[*req_id].filedes, FC_ASPOLL, prdstat, &dumstat);
         errval = (FFSTAT(*prdstat) == 0) ? INPROGRESS: 0;
      }
#     else

         errval = ( IO_DONE(cb_fout[*req_id].stat) == 0)? INPROGRESS: 0;

#     endif

#   elif defined(AIX)
      errval = aio_error(cb_fout[(int)*req_id].aio_handle);
#   else
      errval = aio_error(cb_fout+(int)*req_id);
#   endif
#endif
      switch (errval) {
      case 0: 
          while(aio_req[aio_i] != *req_id && aio_i < MAX_AIO_REQ) aio_i++;
          if(aio_i >= MAX_AIO_REQ) ELIO_ERROR(HANDFAIL, aio_i);

      *req_id = ELIO_DONE; 
      *status = ELIO_DONE;
      aio_req[aio_i] = NULL_AIO;
      break;
      case INPROGRESS:
      *status = ELIO_PENDING; 
      break;
      default:
          return PROBFAIL;
      }
  }

#ifdef PABLO
  PABLO_end(pablo_code);
#endif

  return ELIO_OK;
}


#if defined(CRAY) && defined(FFIO)
static int cray_part_info(char *dirname,long *pparts,long *sparts)
{
  struct statfs stats;
  long temp,count=0;

  if(statfs(dirname, &stats, sizeof(struct statfs), 0) == -1) return -1;

  temp = stats.f_priparts;
  while(temp != 0){
      count++;
      temp <<= 1;
  }
 *pparts = count;

 if(stats.f_secparts != 0){

    temp = (stats.f_secparts << count);
    count = 0;
    while(temp != 0){
           count++;
           temp <<= 1;
    }
    *sparts = count;
 }
 return ELIO_OK;

}

#endif


/*\ Noncollective File Open
\*/
Fd_t  elio_open(const char* fname, int type, int mode)
{
  Fd_t fd=NULL;
  stat_t statinfo;
  int ptype=0, rc;
  char dirname[ELIO_FILENAME_MAX];

  /*
    Create a file for writing to in lustre with
    a specified pagesize and stripe.
    pagesize = 1048576;
    lustre_stripe_count = 32;
    are good choices.
  */
#ifdef USE_LUSTRE
  struct lov_mds_md stripecfg;
  int    lustre_file;
  int  lustre_stripe_count;
  int  pagesize;
  pagesize = 1048576;
  lustre_stripe_count = 32;
#endif

#ifdef PABLO
  int pablo_code = PABLO_elio_open;
  PABLO_start( pablo_code );
#endif

  if(first_elio_init) elio_init();

   switch(type){
     case ELIO_W:  ptype = O_CREAT | O_TRUNC | O_WRONLY;
                   break;
     case ELIO_R:  ptype = O_RDONLY;
                   break;
     case ELIO_RW: ptype = O_CREAT | O_RDWR;
                   break;
     default:      
                   ELIO_ERROR_NULL(MODEFAIL, type);
   }

#if defined(WIN32) || defined(CYGNUS)
   ptype |= O_BINARY;
#endif

  if((fd = (Fd_t ) malloc(sizeof(fd_struct)) ) == NULL) 
                   ELIO_ERROR_NULL(ALOCFAIL, 0);

  if( (rc = elio_dirname(fname, dirname, ELIO_FILENAME_MAX)) != ELIO_OK) {
                   free(fd);
                   ELIO_ERROR_NULL(rc, 0);
  }

  if( (rc = elio_stat(dirname, &statinfo)) != ELIO_OK) {
                   free(fd);
                   ELIO_ERROR_NULL(rc, 0);
  }

  fd->fs = statinfo.fs;
  fd->mode = mode;
  fd->type = type;
  fd->extent = 0;
  fd->next = NULL;
  
#ifdef USE_LUSTRE
  lustre_file = (strncmp(fname,"/dtemp",6) == 0) && (access(fname, F_OK) != 0) && (ME() == 0);
  if (lustre_file) {
    ptype = ptype | O_LOV_DELAY_CREATE ;
  }
#endif

#if defined(CRAY) && defined(FFIO)
  {
    struct ffsw ffstat;
    long pparts, sparts, cbits, cblocks;
    extern long _MPP_MY_PE;
    char *ffio_str="cache:256"; /*  intern I/O buffer/cache 256*4096 bytes */ 
                                /*  JN: we do not want read-ahead write-behind*/

    if(cray_part_info(dirname,&pparts,&sparts) != ELIO_OK){
                   free(fd);
                   ELIO_ERROR_NULL(STATFAIL, 0);
    }

    ptype |= ( O_BIG | O_PLACE | O_RAW );
    cbits = (sparts != 0) ? 1 : 0;

    if( sparts != 0) {

      /* stripe is set so we only select secondary partitions with cbits */
      if(mode == ELIO_SHARED){
         cbits = ~((~0L)<<PARIO_MIN(32,sparts)); /* use all secondary partitions */
         cblocks = 100;
      }else{
         cbits = 1 << (_MPP_MY_PE%sparts);  /* round robin over s part */
      }

      cbits <<= pparts;        /* move us out of the primary partitions */

     }

     
/*     printf ("parts=%d cbits = %X\n",sparts,cbits);*/

     if(mode == ELIO_SHARED)
      fd->fd = OPEN(fname, ptype, FOPEN_MODE, cbits, cblocks, &ffstat, NULL);
     else
      fd->fd = OPEN(fname, ptype, FOPEN_MODE, 0L   , 0      , &ffstat, ffio_str);

  }
#else
  fd->fd = OPEN(fname, ptype, FOPEN_MODE );
#endif

  if( (int)fd->fd == -1) {
                   free(fd);
                   ELIO_ERROR_NULL(OPENFAIL, 0);
  }
  
  fd->name = strdup(fname);

#ifdef USE_LUSTRE
    if (lustre_file) {
      stripecfg.lmm_magic = LOV_MAGIC;
      stripecfg.lmm_pattern = 0; /* Only available option for now. */
      stripecfg.lmm_stripe_size = pagesize; /* Stripe size in bytes. */
      stripecfg.lmm_stripe_count  = lustre_stripe_count;
      if (ioctl((int)fd->fd, LL_IOC_LOV_SETSTRIPE, &stripecfg) < 0) {
        fprintf(stderr,
              "fp_create_out_filefp: Error: unable to stripe %s file.\n"
              "error was %s\n",
              fname,strerror(errno));
        fflush(stderr);
        free(fd);
        ELIO_ERROR_NULL(OPENFAIL, 0);
      }
    } /* end if (luster_file) (is in /dtemp) */
#endif

#ifdef PABLO
  PABLO_end(pablo_code);
#endif

  return(fd);
}

/*\ Close File
\*/
int elio_close(Fd_t fd)
{
    int status = ELIO_OK;
#ifdef PABLO
    pablo_code = PABLO_elio_close;
    PABLO_start( pablo_code );
#endif

    if (fd->next)
      status = elio_close((Fd_t) fd->next);

    /*printf("Closing extent %d name %s\n", fd->extent, fd->name);*/
    if(CLOSE(fd->fd)==-1 || (status != ELIO_OK)) 
      ELIO_ERROR(CLOSFAIL, 0);

    free(fd->name);
    free(fd);

#ifdef PABLO
    PABLO_end(pablo_code);
#endif
    return ELIO_OK;
}



/*\ Close File
\*/
int elio_fsync(Fd_t fd)
{
    int status = ELIO_OK;

#ifdef ELIO_FSYNC
    if (fd->next)
      status = elio_fsync((Fd_t) fd->next);

    /* printf("syncing extent %d name %s\n", fd->extent, fd->name); */
    /*   if(ELIO_FSYNC(fd->fd)==-1 || (status != ELIO_OK)) */
#ifndef WIN32
#if !defined(__INTERIX) 
#ifdef CATAMOUNT
      status = fsync((Fd_t) fd->next);
#else
    sync();
#endif
#endif
#endif
    if(ELIO_FSYNC(fd->fd)==-1 )
      ELIO_ERROR(FSYNCFAIL, 0);
#endif

    return ELIO_OK;
}


/*\ Delete File
\*/
int elio_delete(const char* filename)
{
    int rc;

    if (access(filename, F_OK) != 0) /* Succeed if the file does not exist */
      return ELIO_OK;

#ifdef PABLO
    int pablo_code = PABLO_elio_delete;
    PABLO_start( pablo_code );
#endif

    rc = unlink(filename);

    /* Remeber the first rc ... now delete possible extents until
       one fails */

    {
      int extent;
      for (extent=1; extent<MAX_EXTENT; extent++) {
    char fname[ELIO_FILENAME_MAX];
    sprintf(fname,"%sx%3.3d",filename,extent);
    /*printf("Deleting extent %d with name '%s'\n",extent,fname);*/
    if (unlink(fname)) break;
      }
    }
    
    if(rc ==-1) ELIO_ERROR(DELFAIL,0);

#ifdef PABLO
    PABLO_end(pablo_code);
#endif
    return(ELIO_OK);
}



/*\ Initialize ELIO
\*/
void elio_init(void)
{
  if(first_elio_init) {
#     if defined(ASYNC)
           int i;
           for(i=0; i < MAX_AIO_REQ; i++)
         aio_req[i] = NULL_AIO;
#     endif
      first_elio_init = 0;
  }
}


/*\ Return Error String Associated with Given Error Code 
\*/
void elio_errmsg(int code, char *msg)
{
     if(code==ELIO_OK){
         (void) strcpy(msg, ">OK");
         return;
     }
     else if(code == ELIO_PENDING_ERR) code = elio_pending_error;

     if(code<OFFSET || code >OFFSET+ERRLEN) *msg=(char)0;
     else (void) strcpy(msg, errtable[-OFFSET + code]);
}
      

int elio_pending_error=UNKNFAIL;

char *errtable[ERRLEN] ={
">Unable to Seek",
">Write Failed",
">Asynchronous Write Failed",
">Read Failed",
">Asynchronous Read Failed",
">Suspend Failed",
">I/O Request Handle not in Table",
">Incorrect File Mode",
">Unable to Determine Directory",
">Stat For Specified File or Directory Failed",
">Open Failed",
">Unable To Allocate Internal Data Structure",
">Unsupported Feature",
">Unlink Failed",
">Close Failed",
">Operation Interrupted Too Many Times",
">AIO Return Failed",
">Name String too Long",
">Unable to Determine Filesystem Type",
">Numeric Conversion Error", 
">Incorrect Filesystem/Device Type",
">Error in Probe",
">Unable to Truncate",
">End of File",
">Fsync Failed",
""};
