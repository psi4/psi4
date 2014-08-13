#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 * This header file declares the public C API to tcgmsg.
 */

#include <string.h>
#include <stdlib.h>

#include "sndrcv.h"

void tcg_alt_pbegin(int *argc, char **argv[])
{
    tcgi_pbegin(*argc, *argv);
}


void tcg_brdcst(long type, void *buf, long lenbuf, long originator)
{
    long atype = type;
    long alenbuf = lenbuf;
    long aoriginator = originator;

    BRDCST_(&atype, buf, &alenbuf, &aoriginator);
}


void tcg_dgop(long type, double *x, long n, char *op)
{
    long atype = type;
    long an = n;
    double *ax;
    int i;

    ax = (double*)malloc(n * sizeof(double));
    for (i=0; i<n; ++i) {
        ax[i] = (double)x[i];
    }
    DGOP_(&atype, ax, &an, op, strlen(op));
    for (i=0; i<n; ++i) {
        x[i] = (double)ax[i];
    }
    free(ax);
}


double tcg_drand48()
{
    return ( (double) random() ) * 4.6566128752458e-10;
}


void tcg_igop(long type, long *x, long n, char *op)
{
    long atype = type;
    long an = n;
    long *ax;
    int i;

    ax = (long*)malloc(n * sizeof(long));
    for (i=0; i<n; ++i) {
        ax[i] = (long)x[i];
    }
    IGOP_(&atype, ax, &an, op, strlen(op));
    for (i=0; i<n; ++i) {
        x[i] = (long)ax[i];
    }
    free(ax);
}


void tcg_llog()
{
    LLOG_();
}


long tcg_mdtob(long n)
{
    long an = n;

    return MDTOB_(&an);
}


long tcg_mdtoi(long n)
{
    long an = n;

    return MDTOI_(&an);
}


long tcg_mitob(long n)
{
    long an = n;

    return MITOB_(&an);
}


long tcg_mitod(long n)
{
    long an = n;

    return MITOD_(&an);
}


long tcg_mtime()
{
    return MTIME_();
}


long tcg_niceftn(long ival)
{
    long aival = ival;

    return NICEFTN_(&aival);
}


long tcg_nnodes()
{
    return NNODES_();
}


long tcg_nodeid()
{
    return NODEID_();
}


long tcg_nxtval(long mproc)
{
    long amproc = mproc;

    return NXTVAL_(&amproc);
}


void tcg_error(char *msg, long code)
{
    long acode = code;

    Error(msg, acode);
}


void tcg_pbegin(int argc, char **argv)
{
    tcgi_pbegin(argc,argv);
}


void tcg_pbeginf()
{
    PBEGINF_();
}


void tcg_pbginf()
{
    PBGINF_();
}


void tcg_pend()
{
    PEND_();
}


void tcg_pfcopy(long type, long node0, char *fname)
{
    long atype = type;
    long anode0 = node0;

    PFILECOPY_(&atype, &anode0, fname);
}


long tcg_probe(long type, long node)
{
    long atype = type;
    long anode = node;

    return PROBE_(&atype, &anode);
}


void tcg_rcv(long type, void *buf, long lenbuf, long *lenmes,
        long nodeselect, long *nodefrom, long sync)
{
    long atype = type;
    long alenbuf = lenbuf;
    long alenmes = *lenmes;
    long anodeselect = nodeselect;
    long anodefrom = *nodefrom;
    long async = sync;

    RCV_(&atype, buf, &alenbuf, &alenmes, &anodeselect, &anodefrom, &async);

    *lenmes = alenmes;
    *nodefrom = anodefrom;
}


void tcg_setdbg(long value)
{
    long avalue = value;

    SETDBG_(&avalue);
}


void tcg_snd(long type, void *buf, long lenbuf, long node, long sync)
{
    long atype = type;
    long alenbuf = lenbuf;
    long anode = node;
    long async = sync;

    SND_(&atype, buf, &alenbuf, &anode, &async);
}


void tcg_srand48(long seed)
{
    unsigned int aseed = seed;

    srandom(aseed);
}


void tcg_stats()
{
    STATS_();
}


void tcg_synch(long type)
{
    long atype = type;

    SYNCH_(&atype);
}


long tcg_ready()
{
    return TCGREADY_();
}


double tcg_time()
{
    return TCGTIME_();
}


void tcg_waitcom(long node)
{
    long anode = node;

    WAITCOM_(&anode);
}
