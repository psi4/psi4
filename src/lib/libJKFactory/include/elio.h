/* file name: elio.h  to be included by all apps that use ELIO */


/*include file that contains some common constants, also included by fortran */
#include "chemio.h" 
#include <sys/types.h>
#include <stdio.h>

#define   ELIO_UFS         0    /* Unix filesystem type */
#define   ELIO_PFS           1    /* PFS Intel parallel filesystem type */
#define   ELIO_PIOFS         2    /* IBM SP parallel filesystem type */
#define   ELIO_PENDING_ERR -44  /* error code for failing elio_(g)open */
#define   ELIO_SHARED       77
#define   ELIO_PRIVATE      88


/*********************** type definitions for ELIO interface *****************/
typedef long Size_t;         /* size of I/O request type */ 
typedef double Off_t;         /* size of offset type - double = 56 bit integer*/
typedef struct {
  int   fd;            /* OS handle */
  int   fs;            /* ??? */
  int   mode;            /* ??? */
  int   type;            /* ??? */
  char  *name;                /* Name of physical file */
  int   extent;            /* Counts extents of logical files */
  struct fd_struct *next;    /* Next extent */
} fd_struct;            /* file descriptor type definition */
typedef fd_struct* Fd_t;
#if defined(IBM) || defined(SOLARIS) || defined(HPUX)
typedef unsigned long long avail_t;
#else
typedef unsigned long avail_t;
#endif

typedef struct{
  int   fs;
  avail_t  avail;
} stat_t;
typedef long io_request_t;   /* asynchronous I/O request type */


/********************** prototypes for elio functions ***********************/
extern Size_t elio_read(Fd_t fd, Off_t offset, void *buf, Size_t bytes); 
extern int    elio_aread(Fd_t fd, Off_t offset, void *buf,
                         Size_t bytes, io_request_t *req_id);
extern Size_t elio_write(Fd_t fd, Off_t offset, const void *buf, 
                         Size_t bytes); 
extern int    elio_awrite(Fd_t fd, Off_t offset, const void *buf,
                          Size_t bytes, io_request_t *req_id);
extern int    elio_wait(io_request_t *id);
extern int    elio_probe(io_request_t *id, int* status);
extern int    elio_delete(const char *filename);
extern Fd_t   elio_open(const char *fname, int type, int mode);
extern int    elio_close(Fd_t fd);
extern int    elio_stat(char *fname, stat_t *statinfo);
extern int    elio_dirname(const char *fname, char *statinfo, int len);
extern int    elio_truncate(Fd_t fd, Off_t length);
extern int    elio_length(Fd_t fd, Off_t *length);
extern void   elio_errmsg(int code, char *msg);
extern int    elio_fsync(Fd_t fd);

