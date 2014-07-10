#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#include "eaf.h"
#include "eafP.h"
#include "typesf2c.h"
#include "farg.h"

#if defined(CRAY) && defined(__crayx1)
#undef CRAY
#endif

#define eaf_aread_       F77_FUNC_(eaf_aread,EAF_AREAD)
#define eaf_awrite_      F77_FUNC_(eaf_awrite,EAF_AWRITE)
#define eaf_close_       F77_FUNC_(eaf_close,EAF_CLOSE)
#define eaf_delete_      F77_FUNC_(eaf_delete,EAF_DELETE)
#define eaf_eof_         F77_FUNC_(eaf_eof,EAF_EOF)
#define eaf_errmsg_      F77_FUNC_(eaf_errmsg,EAF_ERRMSG)
#define eaf_error_       F77_FUNC_(eaf_error,EAF_ERROR)
#define eaf_length_      F77_FUNC_(eaf_length,EAF_LENGTH)
#define eaf_open_        F77_FUNC_(eaf_open,EAF_OPEN)
#define eaf_print_stats_ F77_FUNC_(eaf_print_stats,EAF_PRINT_STATS)
#define eaf_probe_       F77_FUNC_(eaf_probe,EAF_PROBE)
#define eaf_read_        F77_FUNC_(eaf_read,EAF_READ)
#define eaf_stat_        F77_FUNC_(eaf_stat,EAF_STAT)
#define eaf_truncate_    F77_FUNC_(eaf_truncate,EAF_TRUNCATE)
#define eaf_util_random_ F77_FUNC_(eaf_util_random,EAF_UTIL_RANDOM)
#define eaf_util_szint_  F77_FUNC_(eaf_util_szint,EAF_UTIL_SZINT)
#define eaf_wait_        F77_FUNC_(eaf_wait,EAF_WAIT)
#define eaf_write_       F77_FUNC_(eaf_write,EAF_WRITE)


static int valid_offset(double offset)
{
    return ((offset - ((double) ((eaf_off_t) offset))) == 0.0);
}

    
Integer FATR eaf_write_(Integer *fd, double *offset, const void *buf, 
           Integer *bytes)
{
    if (!valid_offset(*offset)) return EAF_ERR_NONINTEGER_OFFSET;
    return (Integer) EAF_Write((int) *fd, (eaf_off_t) *offset, buf, 
                   (size_t) *bytes);
}


Integer FATR eaf_awrite_(Integer *fd, double *offset, const void *buf, 
            Integer *bytes, Integer *req_id)
{
    int req, status;

    if (!valid_offset(*offset)) return EAF_ERR_NONINTEGER_OFFSET;
    status = EAF_Awrite((int) *fd, (eaf_off_t) *offset, buf, 
            (size_t) *bytes, &req);
    *req_id = (Integer) req;
    return (Integer) status;
}


Integer FATR eaf_read_(Integer *fd, double *offset, void *buf, Integer *bytes)
{
    if (!valid_offset(*offset)) return EAF_ERR_NONINTEGER_OFFSET;
    return (Integer) EAF_Read((int) *fd, (eaf_off_t) *offset, buf, 
                  (size_t) *bytes);
}


Integer FATR eaf_aread_(Integer *fd, double *offset, void *buf, 
            Integer *bytes, Integer *req_id)
{
    int req, status;

    if (!valid_offset(*offset)) return EAF_ERR_NONINTEGER_OFFSET;
    status = EAF_Aread((int) *fd, (eaf_off_t) *offset, buf, 
               (size_t) *bytes, &req);
    *req_id = (Integer) req;
    return (Integer) status;
}


Integer FATR eaf_wait_(Integer *fd, Integer *id)
{
    return (Integer) EAF_Wait((int) *fd, (int) *id);
}


void FATR eaf_print_stats_(Integer *fd)
{
    EAF_Print_stats((int) *fd);
}


Integer FATR eaf_truncate_(Integer *fd, double *length)
{
    return (Integer) EAF_Truncate((int) *fd, (eaf_off_t) *length);
}


Integer FATR eaf_probe_(Integer *id, Integer *status)
{
    int s, code;

    code = EAF_Probe((int) *id, &s);
    *status = (Integer) s;

    return (Integer) code;
}


Integer FATR eaf_close_(Integer *fd)
{
    return (Integer) EAF_Close((int) *fd);
}


Integer FATR eaf_length_(Integer *fd, double *length)
{
    eaf_off_t len;
    int code;

    code = EAF_Length((int) *fd, &len);
    if (!code) *length = (double) len;

    return code;
}


logical eaf_eof_(Integer *code)
{
    return (logical) EAF_Eof((int) *code);
}


Integer FATR eaf_open_(
#ifdef F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *fname,
        Integer *type,
        Integer *fd,
        int flen
#else
        char *fname,
        int flen,
        Integer *type,
        Integer *fd
#endif
        )
{
    char buf[1024];
    int code, tmp;

    ga_f2cstring(fname, flen, buf, sizeof(buf));
    /* return (Integer) EAF_ERR_TOO_LONG; TODO errcheck? */

    code = EAF_Open(buf, (int) *type, &tmp);
    *fd = (Integer) tmp;

    return (Integer)code;
}


Integer FATR eaf_delete_(char *fname, int flen)
{
    char buf[1024];

    ga_f2cstring(fname, flen, buf, sizeof(buf));
    /* return (Integer) EAF_ERR_TOO_LONG; TODO errcheck? */

    return (Integer) EAF_Delete(buf);
}


Integer FATR eaf_stat_(
#ifdef F2C_HIDDEN_STRING_LENGTH_AFTER_ARGS
        char *path,
        Integer *avail_kb,
        char *fstype,
        int pathlen,
        int fslen
#else
        char *path,
        int pathlen,
        Integer *avail_kb,
        char *fstype,
        int fslen
#endif
        )
{
    char pbuf[1024];
    char fbuf[32];

    int code, kb;

    ga_f2cstring(path, pathlen, pbuf, sizeof(pbuf));
    /* return (Integer) EAF_ERR_TOO_LONG; TODO errcheck? */

    code = EAF_Stat(pbuf, &kb, fbuf, sizeof(fbuf));

    if (!code) {
        ga_c2fstring(fbuf, fstype, fslen);
        /* return (Integer) EAF_ERR_TOO_SHORT; TODO errcheck? */
        *avail_kb = (double) kb;
    }

    return code;
}
    

void FATR eaf_errmsg_(Integer *code,  char *msg, int msglen)
{
    char buf[80];

    EAF_Errmsg((int) *code, buf);

    (void) ga_c2fstring(buf, msg, msglen);
}


double FATR eaf_util_random_(Integer* seed)
{
#ifdef NEC
  if(*seed) srand((unsigned) *seed);
  return ((double) rand())*4.6566128752458e-10;
#else
  if(*seed) srandom((unsigned) *seed);
  return ((double) random())*4.6566128752458e-10;
#endif
}


Integer FATR eaf_util_szint_()
{
  return (Integer)sizeof(Integer);
}
