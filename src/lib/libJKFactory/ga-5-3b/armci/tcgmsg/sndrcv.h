/** @file
 * This header file declares stubs and show prototypes of the 
 * public sndrcv calls
 *
 * srftoc.h contains macros which define the names of c routines
 * accessible from FORTRAN and vice versa
 */
#ifndef SNDRCV_H_
#define SNDRCV_H_

#include "msgtypesc.h"
#include "srftoc.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void   BRDCST_(long *type, void *buf, long *lenbuf, long *originator);
extern void   DGOP_(long *type, double *x, long *n, char *op, int oplen);
extern double DRAND48_();
extern void   IGOP_(long *type, long *x, long *n, char *op, int oplen);
extern void   LLOG_();
extern long   MDTOB_(long *n);
extern long   MDTOI_(long *n);
extern long   MITOB_(long *n);
extern long   MITOD_(long *n);
extern long   MTIME_();
extern long   NICEFTN_(long *ival);
extern long   NNODES_();
extern long   NODEID_();
extern long   NXTVAL_(long *mproc);
extern void   PARERR_(long *code);
extern void   PBEGINF_();
extern void   PBGINF_();
extern void   PEND_();
extern void   PFCOPY_(long *type, long *node0, char *filename, int flen);
extern void   PFILECOPY_(long *type, long *node0, char *filename);
extern long   PROBE_(long *type, long *node);
extern void   RCV_(long *type, void *buf, long *lenbuf, long *lenmes, long *nodeselect, long *nodefrom, long *sync);
extern void   SETDBG_(long *value);
extern void   SND_(long *type, void *buf, long *lenbuf, long *node, long *sync);
extern void   SRAND48_(long *seed);
extern void   STATS_();
extern void   SYNCH_(long *type);
extern long   TCGREADY_();
extern double TCGTIME_();
extern void   WAITCOM_(long *node);

/*
  Miscellaneous routines for internal use only?
*/

extern void Error(char *string, long integer);
extern void MtimeReset();
extern void PrintProcInfo();
extern void RemoteConnect(long a, long b, long c);
extern void tcgi_pbegin(int argc, char **argv);
extern void USleep(long us);

#ifdef __cplusplus
}
#endif

#endif /* SNDRCV_H_ */
