/** @file
 * This header file declares the public C API to tcgmsg.
 */
#ifndef TCGMSG_H_
#define TCGMSG_H_

#ifdef __cplusplus
extern "C" {
#endif

extern void   tcg_alt_pbegin();
extern void   tcg_brdcst(long type, void *buf, long lenbuf, long originator);
extern void   tcg_dgop(long type, double *x, long n, char *op);
extern double tcg_drand48();
extern void   tcg_error(char *msg, long code);
extern void   tcg_igop(long type, long *x, long n, char *op);
extern void   tcg_llog();
extern long   tcg_mdtob(long n);
extern long   tcg_mdtoi(long n);
extern long   tcg_mitob(long n);
extern long   tcg_mitod(long n);
extern long   tcg_mtime();
extern long   tcg_niceftn(long ival);
extern long   tcg_nnodes();
extern long   tcg_nodeid();
extern long   tcg_nxtval(long mproc);
extern void   tcg_parerr(long code);
extern void   tcg_pbeginf();
extern void   tcg_pbegin(int argc, char **argv);
extern void   tcg_pbginf();
extern void   tcg_pend();
extern void   tcg_pfcopy(long type, long node0, char *fname);
extern long   tcg_probe(long type, long node);
extern void   tcg_rcv(long type, void *buf, long lenbuf, long *lenmes, long nodeselect, long *nodefrom, long sync);
extern long   tcg_ready();
extern void   tcg_setdbg(long value);
extern void   tcg_snd(long type, void *buf, long lenbuf, long node, long sync);
extern void   tcg_srand48(long seed);
extern void   tcg_stats();
extern void   tcg_synch(long type);
extern double tcg_time();
extern void   tcg_waitcom(long node);

#ifdef __cplusplus
}
#endif

#endif /* TCGMSG_H_ */
