#ifndef TCGMSGP_H_
#define TCGMSGP_H_

#include <mpi.h>

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif

#include "sndrcv.h"

#if SIZEOF_F77_INTEGER == SIZEOF_SHORT
#   define TCG_INT MPI_SHORT
#   define FMT_INT "%hd"
#elif SIZEOF_F77_INTEGER == SIZEOF_INT
#   define TCG_INT MPI_INT
#   define FMT_INT "%d"
#elif SIZEOF_F77_INTEGER == SIZEOF_LONG
#   define TCG_INT MPI_LONG
#   define FMT_INT "%ld"
#elif SIZEOF_F77_INTEGER == SIZEOF_LONG_LONG
#   define TCG_INT MPI_LONG_LONG_INT
#   define FMT_INT "%lld"
#else
#   error Cannot determine TCG_INT
#endif

#if SIZEOF_F77_DOUBLE_PRECISION == SIZEOF_FLOAT
#   define TCG_DBL MPI_FLOAT
#   define FMT_DBL "%f"
#elif SIZEOF_F77_DOUBLE_PRECISION == SIZEOF_DOUBLE
#   define TCG_DBL MPI_DOUBLE
#   define FMT_DBL "%f"
#elif SIZEOF_F77_DOUBLE_PRECISION == SIZEOF_LONG_DOUBLE
#   define TCG_DBL MPI_LONG_DOUBLE
#   define FMT_DBL "%Lf"
#else
#   error Cannot determine TCG_DBL
#endif

#define MAX_PROCESS 100000
#define TYPE_NXTVAL 33333

extern MPI_Comm TCGMSG_Comm;
extern int      SR_parallel;
extern int      SR_single_cluster;
extern long     DEBUG_;
extern int      _tcg_initialized;

#define TCG_MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define TCG_MIN(a,b) (((a) <= (b)) ? (a) : (b))
#define TCG_ABS(a)   (((a) >= 0)   ? (a) : (-(a)))

#define TCG_ERR_LEN 80
#define ERR_STR_LEN TCG_ERR_LEN + MPI_MAX_ERROR_STRING
extern char tcgmsg_err_string[ERR_STR_LEN];

#define tcgmsg_test_statusM(_str, _status)\
{\
  if( _status != MPI_SUCCESS){\
      int _tot_len, _len = TCG_MIN(ERR_STR_LEN, strlen(_str));\
      strncpy(tcgmsg_err_string, _str, _len);\
      MPI_Error_string( _status, tcgmsg_err_string + _len, &_tot_len);\
      Error(tcgmsg_err_string, (int)_status);\
  }\
}

extern void finalize_nxtval();
extern void install_nxtval(int *argc, char **argv[]);

#endif /* TCGMSGP_H_ */
