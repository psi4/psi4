#ifndef  ELIOP_H
#define  ELIOP_H

#ifdef WIN32
#include <io.h>
#include "winutil.h"
#define F_OK 00
#endif

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "typesf2c.h"

#define PRINT_AND_ABORT(msg, val) GA_Error(msg, (int)val)
#ifndef GLOBAL_H
extern void GA_Error(char*, int);
#endif

#if (defined(SP) || defined(SP1))
#define PIOFS 1
#endif


#if (defined(CRAY) && !defined(__crayx1)) || defined(NEC)
#        include <sys/statfs.h>
#        define  STATVFS statfs
#elif defined(__FreeBSD__) || defined(MACX)
#        include <sys/param.h>
#        include <sys/mount.h>
#        define  STATVFS statfs
#        define NO_F_FRSIZE 
#elif defined(WIN32)
#        define  STATVFS _stat 
#        define  S_ISDIR(mode) ((mode&S_IFMT) == S_IFDIR)
#        define  S_ISREG(mode) ((mode&S_IFMT) == S_IFREG)
#elif defined(CYGNUS) ||  defined(LINUX)  ||  defined(CYGWIN) || defined(BGL) || defined(BGP) || defined(BGQ) || defined(HPUX)
#        include <sys/vfs.h>
#        define  STATVFS statfs
#        define NO_F_FRSIZE 
#else
#        include <sys/statvfs.h>
#        define  STATVFS statvfs
#endif

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <fcntl.h>

#if (defined(CRAY) && defined(FFIO))
#        include <ffio.h>
#        include <sys/fstyp.h>
#        include <sys/fsid.h>
#endif


#include "elio.h"
#include "pablo.h"

extern int                   _elio_Errors_Fatal;
extern void                  elio_init(void);
extern int                   elio_pending_error;


#if !defined(PRINT_AND_ABORT)
#   if defined(SUN) && !defined(SOLARIS)
      extern int fprintf();
      extern void fflush();
#   endif
#   define PRINT_AND_ABORT(msg, val){\
     fprintf(stderr, "ELIO fatal error: %s %ld\n", msg,  val);\
     fprintf(stdout, "ELIO fatal error: %s %ld\n", msg,  val);\
     fflush(stdout);\
     exit(val);\
   }
#endif

/**************************** Error Macro ******************************/
/* ELIO defines error macro called in case of error
 * the macro can also use user-provided error routine PRINT_AND_ABORT
 * defined as macro to do some cleanup in the application before
 * aborting
 * The requirement is that PRINT_AND_ABORT is defined before
 * including ELIO header file - this file
 */

#define ELIO_ERROR_NULL(code, val){\
 PABLO_end(pablo_code);\
 if(! _elio_Errors_Fatal){\
     elio_pending_error= code;\
     return NULL;\
 }\
 if( _elio_Errors_Fatal)\
     PRINT_AND_ABORT(errtable[code-OFFSET], val);\
}

#define ELIO_ERROR(code, val) { \
 PABLO_end(pablo_code);\
 if(! _elio_Errors_Fatal) return(code);\
 else PRINT_AND_ABORT(errtable[code-OFFSET], val);\
}


/* error codes and messages */

#define ERRLEN 26
#define OFFSET    (-2000)
#define SEEKFAIL  (OFFSET + 0)
#define WRITFAIL  (OFFSET + 1)
#define AWRITFAIL (OFFSET + 2)
#define READFAIL  (OFFSET + 3)
#define AREADFAIL (OFFSET + 4)
#define SUSPFAIL  (OFFSET + 5)
#define HANDFAIL  (OFFSET + 6)
#define MODEFAIL  (OFFSET + 7)
#define DIRFAIL   (OFFSET + 8)
#define STATFAIL  (OFFSET + 9)
#define OPENFAIL  (OFFSET + 10)
#define ALOCFAIL  (OFFSET + 11)
#define UNSUPFAIL (OFFSET + 12)
#define DELFAIL   (OFFSET + 13)
#define CLOSFAIL  (OFFSET + 14)
#define INTRFAIL  (OFFSET + 15)
#define RETUFAIL  (OFFSET + 16)
#define LONGFAIL  (OFFSET + 17)
#define FTYPFAIL  (OFFSET + 18)
#define CONVFAIL  (OFFSET + 19)
#define TYPEFAIL  (OFFSET + 20)
#define PROBFAIL  (OFFSET + 21)
#define TRUNFAIL  (OFFSET + 22)
#define EOFFAIL   (OFFSET + 23)
#define FSYNCFAIL (OFFSET + 24)
#define UNKNFAIL  (OFFSET + 25)

extern  char *errtable[ERRLEN];

#define ELIO_FILENAME_MAX 1024
#define SDIRS_INIT_SIZE 1024

#endif
