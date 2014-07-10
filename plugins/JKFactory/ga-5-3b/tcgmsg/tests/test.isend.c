#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv5.0/test.isend.c,v 1.2 1994-12-30 20:56:23 d3h325 Exp $ */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_MEMORY_H
#   include <memory.h>
#endif

#include "evlog.h"
#include "tcgmsg.h"

extern unsigned char CheckByte();


/**
 * Process 0 sleeps for 20 seconds and then sends each
 * process a message.  The other processes sleep for periods
 * of 2 seconds and probes until it gets a message.  All processes
 * respond to process 0 which recieves using a wildcard probe.
 */
static void TestProbe()
{
    long type_syn = 32;
    long type_msg = 33;
    long type_ack = 34;
    long me = NODEID_();
    char buf;
    long lenbuf = sizeof buf;
    long sync = 1;


    if (me == 0) {
        (void) printf("Probe test ... processes should sleep for 20s only\n");
        (void) printf("----------\n\n");
        (void) fflush(stdout);
    }

    SYNCH_(&type_syn);

    if (me == 0) {
        long nproc = NNODES_();
        long anyone = -1;
        long ngot = 0;
        long node;

        (void) sleep((unsigned) 20);

        for (node=1; node<nproc; node++) {
            SND_(&type_msg, &buf, &lenbuf, &node, &sync);
            (void) printf("    Sent message to %ld\n", node);
            (void) fflush(stdout);
        }

        while (ngot < (nproc-1))
            if (PROBE_(&type_ack, &anyone)) {
                RCV_(&type_ack, &buf, &lenbuf, &lenbuf, &anyone, &node, &sync);
                (void) printf("    Got response from %ld\n", node);
                (void) fflush(stdout);
                ngot++;
            }
    }
    else {
        long node = 0;
        while (!PROBE_(&type_msg, &node)) {
            (void) printf("    Node %ld sleeping\n", me);
            (void) fflush(stdout);
            (void) sleep((unsigned) 2);
        }
        (void) printf("    Node %ld detected available message\n", me);
        (void) fflush(stdout);
        RCV_(&type_msg, &buf, &lenbuf, &lenbuf, &node, &node, &sync);
        (void) printf("    Node %ld sending response\n", me);
        (void) fflush(stdout);
        SND_(&type_ack, &buf, &lenbuf, &node, &sync);
    }

    SYNCH_(&type_syn);
}
    

static void TestGlobals()
{
#define MAXLENG 256*1024
    double *dtest;
    long *itest;
    long len;
    long me = NODEID_(), nproc = NNODES_(), from=NNODES_()-1;
    long itype=3+MSGINT, dtype=4+MSGDBL;

    if (me == 0) {
        (void) printf("Global test ... test brodcast, igop and dgop\n----------\n\n");
        (void) fflush(stdout);
    }

    if (!(dtest = (double *) malloc((unsigned) (MAXLENG*sizeof(double)))))
        Error("TestGlobals: failed to allocated dtest", (long) MAXLENG);
    if (!(itest = (long *) malloc((unsigned) (MAXLENG*sizeof(long)))))
        Error("TestGlobals: failed to allocated itest", (long) MAXLENG);

    for (len=1; len<MAXLENG; len*=2) {
        long ilen = len*sizeof(long);
        long dlen = len*sizeof(double);
        long i;

        if (me == 0) {
            printf("Test length = %d ... ", len);
            fflush(stdout);
        }

        /* Test broadcast */

        if (me == (nproc-1)) {
            for (i=0; i<len; i++) {
                itest[i] = i;
                dtest[i] = (double) itest[i];
            }
        }
        else {
            for (i=0; i<len; i++) {
                itest[i] = 0;
                dtest[i] = 0.0;
            }
        }
        BRDCST_(&itype, (char *) itest, &ilen, &from);
        BRDCST_(&dtype, (char *) dtest, &dlen, &from);

        for (i=0; i<len; i++)
            if (itest[i] != i || dtest[i] != (double) i)
                Error("TestGlobal: broadcast failed", (long) i);

        if (me == 0) {
            printf("broadcast OK ...");
            fflush(stdout);
        }

        /* Test global sum */

        for (i=0; i<len; i++) {
            itest[i] = i*me;
            dtest[i] = (double) itest[i];
        }

        IGOP_(&itype, itest, &len, "+");
        DGOP_(&dtype, dtest, &len, "+");

        for (i=0; i<len; i++) {
            long iresult = i*nproc*(nproc-1)/2;
            if (itest[i] != iresult || dtest[i] != (double) iresult)
                Error("TestGlobals: global sum failed", (long) i);
        }

        if (me == 0) {
            printf("global sums OK\n");
            fflush(stdout);
        }
    }

    free((char *) itest);
    free((char *) dtest);
}
  

/**
 * Everyone says hi to everyone else
 */
static void Hello()
{
    char buf[30];
    long lenbuf = sizeof buf;
    long type=19 | MSGCHR;
    long node, kode, nodefrom, lenmes;
    long sync = 1;

    if (NODEID_() == 0) {
        (void) printf("Hello test ... show network integrity\n----------\n\n");
        (void) fflush(stdout);
    }

    for (node = 0; node<NNODES_(); node++) {
        if (node == NODEID_()) {
            for (kode = 0; kode<NNODES_(); kode++) {
                (void) sprintf(buf, "Hello to %ld from %ld", kode, NODEID_());
                if (node != kode)
                    SND_(&type, buf, &lenbuf, &kode, &sync);
            }
        }
        else {
            RCV_(&type, buf, &lenbuf, &lenmes, &node, &nodefrom, &sync);
            (void) printf("me=%ld, from=%ld: %s\n",NODEID_(), node, buf);
            (void) fflush(stdout);
        }
    }

}


/**
 * Fill list with n random integers between lo & hi inclusively
 */
static void RandList(long lo, long hi, long *list, long n)
{
    long i, ran;
    double dran;

    for (i=0; i<n; i++) {
        dran = DRAND48_();
        ran = lo + (long) (dran * (double) (hi-lo+1));
        if (ran < lo)
            ran = lo;
        if (ran > hi)
            ran = hi;
        list[i] = ran;
    }
}


/**
 * Stress the system by passing messages between a ranomly selected
 * list of nodes
 */
void Stress()
{
    long me = NODEID_();
    long nproc = NNODES_();
    long type, lenbuf, node, lenmes, nodefrom, i, j, from, to;
    long *list_i, *list_j, *list_n;
#define N_LEN 11
    static long len[N_LEN] = {0,1,2,4,8,4095,4096,4097,16367,16368,16369};
    char *buf1, *buf2;
    long n_stress, mod;
    long sync = 1;

    from = 0;
    lenbuf = sizeof(long);

    if (me == 0) {
        (void) printf("Stress test ... randomly exchange messages\n-----------");
        (void) printf("\n\nInput no. of messages: ");
        (void) fflush(stdout);
        if (scanf("%ld",&n_stress) != 1)
            Error("Stress: error reading n_stress",(long) -1);
        if ( (n_stress <= 0) || (n_stress > 100000) )
            n_stress = 100;
    }
    type = 13 | MSGINT;
    BRDCST_(&type, (char *) &n_stress, &lenbuf, &from);
    type++;

    lenbuf = n_stress * sizeof(long);

#if HAVE_MEMALIGN
    list_i = (long *) memalign(sizeof(long), (unsigned) lenbuf);
#else
    list_i = (long *) malloc((unsigned) lenbuf);
#endif
    if ( list_i == (long *) NULL ) {
        Error("Stress: failed to allocate list_i",n_stress);
    }

#if HAVE_MEMALIGN 
    list_j = (long *) memalign(sizeof(long), (unsigned) lenbuf);
#else
    list_j = (long *) malloc((unsigned) lenbuf);
#endif
    if ( list_j == (long *) NULL ) {
        Error("Stress: failed to allocate list_j",n_stress);
    }
    
#if HAVE_MEMALIGN
    list_n = (long *) memalign(sizeof(long), (unsigned) lenbuf);
#else   
    list_n = (long *) malloc((unsigned) lenbuf);
#endif
    if ( list_n  == (long *) NULL ) {
        Error("Stress: failed to allocate list_n",n_stress);
    }

    if ( (buf1 = malloc((unsigned) 16376)) == (char *) NULL )
        Error("Stress: failed to allocate buf1", (long) 16376);

    if ( (buf2 = malloc((unsigned) 16376)) == (char *) NULL )
        Error("Stress: failed to allocate buf2", (long) 16376);


    if (me == 0) { /* Make random list of node pairs and message lengths */

        RandList((long) 0, (long) (NNODES_()-1), list_i, n_stress);
        RandList((long) 0, (long) (NNODES_()-1), list_j, n_stress);
        RandList((long) 0, (long) (N_LEN-1), list_n, n_stress);
        for (i=0; i<n_stress; i++)
            list_n[i] = len[list_n[i]];
    }

    node = 0;
    BRDCST_(&type, (char *) list_i, &lenbuf, &node);
    type++;
    BRDCST_(&type, (char *) list_j, &lenbuf, &node);
    type++;
    BRDCST_(&type, (char *) list_n, &lenbuf, &node);
    type++;

    type = 8;

    for (j=0; j<16370; j++)
        buf1[j] = (char) (j%127);

    j = 0;
    mod = (n_stress-1)/10 + 1;
    for (i=0; i < n_stress; i++) {

        from   = list_i[i];
        to     = list_j[i];
        lenbuf = list_n[i];
        type++;

        if ( (from < 0) || (from >= nproc) )
            Error("Stress: from is out of range", from);
        if ( (to < 0) || (to >= nproc) )
            Error("Stress: to is out of range", to);

        if (from == to)
            continue;

        if ( (me == 0) && (j%mod == 0) ) {
            (void) printf("Stress: test=%ld: from=%ld, to=%ld, len=%ld\n",
                          i, from, to, lenbuf);
            (void) fflush(stdout);
        }

        j++;

        if (from == me)
            SND_(&type, buf1, &lenbuf, &to, &sync);
        else if (to == me) {
            (void) bzero(buf2, (int) lenbuf);    /* Initialize the receive buffer */
            buf2[lenbuf] = '+';

            RCV_(&type, buf2, &lenbuf, &lenmes, &from, &nodefrom, &sync);

            if (buf2[lenbuf] != '+')
                Error("Stress: overran buffer on receive",lenbuf);
            if (CheckByte((unsigned char *) buf1, lenbuf) != 
                    CheckByte((unsigned char *) buf2, lenbuf))
                Error("Stress: invalid checksum on receive",lenbuf);
            if (lenmes != lenbuf)
                Error("Stress: invalid message length on receive",lenbuf);
        }
    }

    (void) free(buf2);
    (void) free(buf1);
    (void) free((char *) list_n);
    (void) free((char *) list_j);
    (void) free((char *) list_i);
}


void IsendTest()
{
#define N_ISEND_MSG 5
    long me = NODEID_();
    long type = 4;
    long left = (me + NNODES_() - 1) % NNODES_();
    long right = (me + 1) % NNODES_();
    char *source[N_ISEND_MSG], *dest[N_ISEND_MSG];
    long lenbuf, lenmes, from;
    long sync = 1;
    long max_len;
    long ids[N_ISEND_MSG];
    long i;


    if (me == 0) {
        (void) printf("Isend test ... test async. comms\n---------\n\n");
        (void) printf("Input maximum message size: ");
        (void) fflush(stdout);
        if (scanf("%ld", &max_len) != 1)
            Error("IsendTest: error reading max_len", (long) -1);
        if ( (max_len <= 0) || (max_len >= 4*1024*1024) )
            max_len = 256*1024;
    }

    type = 4 | MSGINT;
    lenbuf = sizeof(long);
    i = 0;
    BRDCST_(&type, (char *) &max_len, &lenbuf, &i);

    for (i=0; i<N_ISEND_MSG; i++) {
        if ( (source[i] = malloc((unsigned) max_len)) == (char *) NULL)
            Error("failed to allocate buf",max_len);
        if ( (dest[i] = malloc((unsigned) max_len)) == (char *) NULL)
            Error("failed to allocate buf",max_len);
        /* Fill in source with something verifiable here */
    }

    type = 5;
    SYNCH_(&type);

    type = 6;
    for (i=0; i<N_ISEND_MSG; i++) {
        ids[i] = msg_async_snd(type, source[i], max_len, left);
        type++;
    }

    type = 6;
    for (i=0; i<N_ISEND_MSG; i++) {
        RCV_(&type, dest[i], &max_len, &lenmes, &right, &from, &sync);
        /* Verify buffer here */
        type++;
    }

    for (i=0; i<N_ISEND_MSG; i++)
        msg_wait(ids[i]);

    for (i=0; i<N_ISEND_MSG; i++) {
        free(source[i]);
        free(dest[i]);
    }

    SYNCH_(&type);
}
       

/**
 * Time passing a message round a ring.
 */
void RingTest()
{
    long me = NODEID_();
    long type = 4;
    long left = (me + NNODES_() - 1) % NNODES_();
    long right = (me + 1) % NNODES_();
    char *buf, *buf2;
    unsigned char sum, sum2;
    long lenbuf, lenmes, nodefrom;
    double start, used, rate;
    long max_len;
    long i;
    long sync = 1;

    i = 0;
    lenbuf = sizeof(long);

    if (me == 0) {
        (void) printf("Ring test ... time network performance\n---------\n\n");
        (void) printf("Input maximum message size: ");
        (void) fflush(stdout);
        if (scanf("%ld", &max_len) != 1)
            Error("RingTest: error reading max_len",(long) -1);
        if ( (max_len <= 0) || (max_len >= 4*1024*1024) )
            max_len = 256*1024;
    }
    type = 4 | MSGINT;
    BRDCST_(&type, (char *) &max_len, &lenbuf, &i);

    if ( (buf = malloc((unsigned) max_len)) == (char *) NULL)
        Error("failed to allocate buf",max_len);

    if (me == 0) {
        if ( (buf2 = malloc((unsigned) max_len)) == (char *) NULL)
            Error("failed to allocate buf2",max_len);

        for (i=0; i<max_len; i++)
            buf[i] = (char) (i%127);
    }

    type = 5;
    lenbuf = 1;
    while (lenbuf <= max_len) {
        int nloops = 10 + 1000/lenbuf;
        /*    int nloops = 1; */
        int loop = nloops;
        if (me == 0) {
            sum = CheckByte((unsigned char *) buf, lenbuf);
            (void) bzero(buf2, (int) lenbuf);
            start = TCGTIME_();
            while (loop--) {
                SND_(&type, buf, &lenbuf, &left, &sync);
                RCV_(&type, buf2, &lenbuf, &lenmes, &right, &nodefrom, &sync);
            }
            used = TCGTIME_() - start;
            sum2 = CheckByte((unsigned char *) buf2, lenbuf);
            if (used > 0)
                rate = 1.0e-6 * (double) (NNODES_() * lenbuf) / (double) used;
            else
                rate = 0.0;
            rate = rate * nloops;
            printf("len=%6ld bytes, nloop=%4ld, used=%8.4f s, rate=%8.4f Mb/s (0x%x, 0x%x)\n",
                    lenbuf, nloops, used, rate, sum, sum2);
            (void) fflush(stdout);
        }
        else {
            while (loop--) {
                RCV_(&type, buf, &lenbuf, &lenmes, &right, &nodefrom, &sync);
                SND_(&type, buf, &lenbuf, &left, &sync);
            }
        }
        lenbuf *= 2;
    }

    if (me == 0)
        (void) free(buf2);

    (void) free(buf);
}


/**
 * Test receiveing a message from any node
 */
void RcvAnyTest()
{
    long me = NODEID_();
    long type = 337 | MSGINT;
    char buf[8];
    long i, j, node, lenbuf, lenmes, nodefrom, receiver, n_msg;
    long sync = 1;

    lenbuf = sizeof(long);

    if (me == 0) {
        (void) printf("RCV any test ... check is working!\n-----------\n\n");
        (void) printf("Input node to receive : ");
        (void) fflush(stdout);
        if (scanf("%ld", &receiver) != 1)
            Error("RcvAnyTest: error reading receiver",(long) -1);
        if ( (receiver < 0) || (receiver >= NNODES_()) )
            receiver = NNODES_()-1;
        (void) printf("Input number of messages : ");
        (void) fflush(stdout);
        if (scanf("%ld", &n_msg) != 1)
            Error("RcvAnyTest: error reading n_msg",(long) -1);
        if ( (n_msg <= 0) || (n_msg > 10) )
            n_msg = 5;
    }

    node = 0;
    BRDCST_(&type, (char *) &receiver, &lenbuf, &node);
    type++;
    BRDCST_(&type, (char *) &n_msg, &lenbuf, &node);
    type++;

    lenbuf = 0;

    type = 321;
    for (i=0; i<n_msg; i++) {

        if (me == receiver) {
            for (j = 0; j<NNODES_(); j++)
                if (j !=  me) {
                    node = -1;
                    RCV_(&type, buf, &lenbuf, &lenmes, &node, &nodefrom, &sync);
                    (void) printf("RcvAnyTest: message received from %ld\n",nodefrom);
                    (void) fflush(stdout);
                }
        }
        else
            SND_(&type, buf, &lenbuf, &receiver, &sync);
    }

}


/**
 * Test the load balancing mechanism
 */
void NextValueTest()
{
    long nproc = NNODES_();
    long me = NODEID_();
    long type = 51 | MSGINT;
    long i, node, lenbuf, n_val, next;
    long ngot, ntimes;
    double start, used, rate;

    lenbuf = sizeof(long);

    if (me == 0) {
        (void) printf("Next value test ... time overhead!\n---------------\n\n");
        (void) printf("Input iteration count : ");
        (void) fflush(stdout);
        if (scanf("%ld", &n_val) != 1)
            Error("NextValueTest: error reading n_val",(long) -1);
        if ( (n_val < 0) || (n_val >= 10000) )
            n_val = 100;
    }
    node = 0;
    BRDCST_(&type, (char *) &n_val, &lenbuf, &node);

    /* Loop thru a few values to visually show it is working */

    next = -1;
    for (i=0; i<10; i++) {

        if (i > next)
            next = NXTVAL_(&nproc);

        if (i == next) {
            (void) printf("node %ld got value %ld\n",me, i);
            (void) fflush(stdout);
        }
    }
    nproc = -nproc;
    next = NXTVAL_(&nproc);
    nproc = -nproc;

    /* Now time it for real .. twice*/

    for (ntimes=0; ntimes<2; ntimes++) {
        if (me == 0)
            start = TCGTIME_();

        next = -1;
        ngot = 0;

        for (i=0; i<n_val; i++) {
            if (i > next)
                next = NXTVAL_(&nproc);
            if (i == next)
                ngot++;
        }

        nproc = -nproc;
        next = NXTVAL_(&nproc);
        nproc = -nproc;

        if (me == 0) {
            used =  TCGTIME_() - start;
            if(ngot){
                rate = used / ngot;
                (void) printf("node 0: From %ld busy iters did %ld, used=%fs per call\n",
                              n_val, ngot, rate);
            }else printf("node 0: From %ld busy iters did 0\n",n_val);
            (void) fflush(stdout);
        }

        type++;
        SYNCH_(&type);
    }
}


void ToggleDebug()
{
    static long on = 0;
    long me = NODEID_();
    long type = 666 | MSGINT;
    long lenbuf = sizeof(long);
    long from=0;
    long node;

    if (me == 0) {
        (void) printf("\nInput node to debug (-1 = all) : ");
        (void) fflush(stdout);
        if (scanf("%ld", &node) != 1)
            Error("ToggleDebug: error reading node",(long) -1);
    }
    BRDCST_(&type, (char *) &node, &lenbuf, &from);

    if ((node < 0) || (node == me)) {
        on = (on + 1)%2;
        SETDBG_(&on);
    }
}


int main(int argc, char **argv)
{
    long type;
    long lenbuf;
    long node, opt;

    tcg_pbegin(argc, argv);

    (void) printf("In process %ld\n", NODEID_());
    (void) fflush(stdout);

    /* Read user input for action */

    lenbuf = sizeof(long);
    node = 0;

    while (1) {

        (void) fflush(stdout);
        if (NODEID_() == 0)
            (void) sleep(1);
        type = 999;
        SYNCH_(&type);
        (void) sleep(1);

        if (NODEID_() == 0) {
again:
            (void) printf("\n\
                          0=quit\n\
                          1=Ring             5=NxtVal\n\
                          2=Stress           6=Global\n\
                          3=Hello            7=Debug\n\
                          4=RcvAny           8=Probe\n\
                          9=Isend\n\n\
                          Enter test number : ");

            (void) fflush(stdout);
            if (scanf("%ld", &opt) != 1)
                Error("test: input of option failed",(long) -1);
            (void) printf("\n");
            (void) fflush(stdout);
            if ( (opt < 0) || (opt > 9) )
                goto again;

        }
        type = 2 | MSGINT;
        BRDCST_(&type, (char *) &opt, &lenbuf, &node);

        switch (opt) {
            case 0:
                if (NODEID_() == 0)
                    STATS_();
                PEND_();
                return 0;

            case 1:
                RingTest();
                break;

            case 2:
                Stress();
                break;

            case 3:
                Hello();
                break;

            case 4:
                RcvAnyTest();
                break;

            case 5:
                NextValueTest();
                break;

            case 6:
                TestGlobals();
                break;

            case 7:
                ToggleDebug();
                break;

            case 8:
                TestProbe();
                break;

            case 9:
                IsendTest();
                break;

            default:
                Error("test: invalid option", opt);
                break;
        }
    }
}
