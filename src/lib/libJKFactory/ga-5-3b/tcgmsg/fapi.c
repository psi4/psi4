/**
 * @file fapi.c
 *
 * Implements the Fortran interface to TCGMSG.
 */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <string.h>
#include <stdlib.h>

#include "farg.h"
#include "typesf2c.h"
#include "srftoc.h"
#include "sndrcv.h"
#define LEN 255

#define _BRDCST_     F77_FUNC(brdcst, BRDCST)
#define _DGOP_       F77_FUNC(dgop, DGOP)
#define _DRAND48_    F77_FUNC(drand48, DRAND48)
#define _IGOP_       F77_FUNC(igop, IGOP)
#define _LLOG_       F77_FUNC(llog, LLOG)
#define _MDTOB_      F77_FUNC(mdtob, MDTOB)
#define _MDTOI_      F77_FUNC(mdtoi, MDTOI)
#define _MITOB_      F77_FUNC(mitob, MITOB)
#define _MITOD_      F77_FUNC(mitod, MITOD)
#define _MTIME_      F77_FUNC(mtime, MTIME)
#define _NICEFTN_    F77_FUNC(niceftn, NICEFTN)
#define _NNODES_     F77_FUNC(nnodes, NNODES)
#define _NODEID_     F77_FUNC(nodeid, NODEID)
#define _NXTVAL_     F77_FUNC(nxtval, NXTVAL)
#define _PARERR_     F77_FUNC(parerr, PARERR)
#define _PBEGINF_    F77_FUNC(pbeginf, PBEGINF)
#define _PBGINF_     F77_FUNC(pbginf, PBGINF)
#define _PEND_       F77_FUNC(pend, PEND)
#define _PFCOPY_     F77_FUNC(pfcopy, PFCOPY)
#define _PFILECOPY_  F77_FUNC(pfilecopy, PFILECOPY)
#define _PROBE_      F77_FUNC(probe, PROBE)
#define _RCV_        F77_FUNC(rcv, RCV)
#define _SETDBG_     F77_FUNC(setdbg, SETDBG)
#define _SND_        F77_FUNC(snd, SND)
#define _SRAND48_    F77_FUNC(srand48, SRAND48)
#define _STATS_      F77_FUNC(stats, STATS)
#define _SYNCH_      F77_FUNC(synch, SYNCH)
#define _TCGREADY_   F77_FUNC(tcgready, TCGREADY)
#define _TCGTIME_    F77_FUNC(tcgtime, TCGTIME)
#define _WAITCOM_    F77_FUNC(waitcom, WAITCOM)


void FATR _BRDCST_(Integer *type, void *buf, Integer *lenbuf, Integer *originator)
{
    long atype = *type;
    long alenbuf = *lenbuf;
    long aoriginator = *originator;

    BRDCST_(&atype, buf, &alenbuf, &aoriginator);
}


void FATR _DGOP_(Integer *type, DoublePrecision *x, Integer *n, char *op, int oplen)
{
    long atype = *type;
    long an = *n;
    double *ax;
    int i;

    ax = (double*)malloc(an * sizeof(double));
    for (i=0; i<an; ++i) {
        ax[i] = (double)x[i];
    }
    DGOP_(&atype, ax, &an, op, strlen(op));
    for (i=0; i<an; ++i) {
        x[i] = (DoublePrecision)ax[i];
    }
    free(ax);
}


DoublePrecision  FATR _DRAND48_()
{
    return ( (double) random() ) * 4.6566128752458e-10;
}


void FATR _IGOP_(Integer *type, Integer *x, Integer *n, char *op, int oplen)
{
    long atype = *type;
    long an = *n;
    long *ax;
    int i;

    ax = (long*)malloc(an * sizeof(long));
    for (i=0; i<an; ++i) {
        ax[i] = (long)x[i];
    }
    IGOP_(&atype, ax, &an, op, strlen(op));
    for (i=0; i<an; ++i) {
        x[i] = (Integer)ax[i];
    }
    free(ax);
}


void FATR _LLOG_()
{
    LLOG_();
}


Integer  FATR _MDTOB_(Integer *n)
{
    long an = *n;

    return MDTOB_(&an);
}


Integer  FATR _MDTOI_(Integer *n)
{
    long an = *n;

    return MDTOI_(&an);
}


Integer  FATR _MITOB_(Integer *n)
{
    long an = *n;

    return MITOB_(&an);
}


Integer  FATR _MITOD_(Integer *n)
{
    long an = *n;

    return MITOD_(&an);
}


Integer  FATR _MTIME_()
{
    return MTIME_();
}


Integer  FATR _NICEFTN_(Integer *ival)
{
    long aival = *ival;

    return NICEFTN_(&aival);
}


Integer  FATR _NNODES_()
{
    return NNODES_();
}


Integer  FATR _NODEID_()
{
    return NODEID_();
}


Integer  FATR _NXTVAL_(Integer *mproc)
{
    long amproc = *mproc;

    return NXTVAL_(&amproc);
}


void FATR _PARERR_(Integer *code)
{
    long acode = *code;

    Error("User detected error in FORTRAN", acode);
}


void FATR _PBEGINF_()
{
    Integer argc = F2C_IARGC();
    Integer i, len;
    char *argv[LEN], arg[LEN];

    for (i=0; i<argc; i++) {
        F2C_GETARG(&i, arg, LEN);
        for(len = LEN-2; len && (arg[len] == ' '); len--);
        len++;
        arg[len] = '\0';
        //printf("%10s, len=%ld\n", arg, (long)len);  fflush(stdout);
        argv[i] = strdup(arg);
    }

    tcgi_pbegin(argc, argv);
}


void FATR _PBGINF_()
{
    PBGINF_();
}


void FATR _PEND_()
{
    PEND_();
}


void FATR _PFCOPY_(Integer *type, Integer *node0, char *filename, int flen)
{
    long atype = *type;
    long anode0 = *node0;

    PFILECOPY_(&atype, &anode0, filename);
}


void FATR _PFILECOPY_(Integer *type, Integer *node0, char *filename)
{
    long atype = *type;
    long anode0 = *node0;

    PFILECOPY_(&atype, &anode0, filename);
}


Integer  FATR _PROBE_(Integer *type, Integer *node)
{
    long atype = *type;
    long anode = *node;

    return PROBE_(&atype, &anode);
}


void FATR _RCV_(Integer *type, void *buf, Integer *lenbuf, Integer *lenmes, Integer *nodeselect, Integer *nodefrom, Integer *sync)
{
    long atype = *type;
    long alenbuf = *lenbuf;
    long alenmes = *lenmes;
    long anodeselect = *nodeselect;
    long anodefrom = *nodefrom;
    long async = *sync;

    RCV_(&atype, buf, &alenbuf, &alenmes, &anodeselect, &anodefrom, &async);

    *lenmes = alenmes;
    *nodefrom = anodefrom;
}


void FATR _SETDBG_(Integer *value)
{
    long avalue = *value;

    SETDBG_(&avalue);
}


void FATR _SND_(Integer *type, void *buf, Integer *lenbuf, Integer *node, Integer *sync)
{
    long atype = *type;
    long alenbuf = *lenbuf;
    long anode = *node;
    long async = *sync;

    SND_(&atype, buf, &alenbuf, &anode, &async);
}


void FATR _SRAND48_(Integer *seed)
{
    unsigned int aseed = *seed;

    srandom(aseed);
}


void FATR _STATS_()
{
    STATS_();
}


void FATR _SYNCH_(Integer *type)
{
    long atype = *type;

    SYNCH_(&atype);
}


Integer  FATR _TCGREADY_()
{
    return TCGREADY_();
}


DoublePrecision  FATR _TCGTIME_()
{
    return TCGTIME_();
}


void FATR _WAITCOM_(Integer *node)
{
    long anode = *node;

    WAITCOM_(&anode);
}

