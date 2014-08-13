#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_MALLOC_H
#   include <malloc.h>
#endif
#if HAVE_MEMORY_H
#   include <memory.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#elif HAVE_WINDOWS_H
#   include <windows.h>
#   define sleep(x) Sleep(1000*(x))
#endif

#define MAXLENG 256*1024
#define N_LEN 11

#include "typesf2c.h"
#include <msgtypesc.h>
#include <tcgmsg.h>

extern unsigned char CheckByte();
extern void Error(const char *, Integer);


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
    long me = tcg_nodeid();
    char buf;
    long lenbuf = sizeof buf;
    long sync = 1;

    if (me == 0) {
        (void) printf("Probe test ... processes should sleep for 20s only\n");
        (void) printf("----------\n\n");
        (void) fflush(stdout);
    }

    tcg_synch(type_syn);

    if (me == 0) {
        long nproc = tcg_nnodes();
        long anyone = -1;
        long ngot = 0;
        long node;

        (void) sleep((unsigned) 20);

        for (node=1; node<nproc; node++) {
            tcg_snd(type_msg, &buf, lenbuf, node, sync);
            (void) printf("    Sent message to %ld\n", node);
            (void) fflush(stdout);
        }

        while (ngot < (nproc-1)) {
            if (tcg_probe(type_ack, anyone)) {
                tcg_rcv(type_ack, &buf, lenbuf, &lenbuf, anyone, &node, sync);
                (void) printf("    Got response from %ld\n", node);
                (void) fflush(stdout);
                ngot++;
            }
        }
    }
    else {
        long node = 0;
        while (!tcg_probe(type_msg, node)) {
            (void) printf("    Node %ld sleeping\n", me);
            (void) fflush(stdout);
            (void) sleep((unsigned) 2);
        }
        tcg_rcv(type_msg, &buf, lenbuf, &lenbuf, node, &node, sync);
        tcg_snd(type_ack, &buf, lenbuf, node, sync);
    }

    tcg_synch(type_syn);
}
    

static void TestGlobals()
{
    double *dtest;
    long *itest;
    long len;
    long me = tcg_nodeid(), nproc = tcg_nnodes(), from=tcg_nnodes()-1;
    long itype=3+MSGINT, dtype=4+MSGDBL;

    if (me == 0) {
        (void) printf("Global test ... test brodcast, igop and dgop\n----------\n\n");
        (void) fflush(stdout);
    }

    if (!(dtest = (double *) malloc((unsigned) (MAXLENG*sizeof(double))))) {
        Error("TestGlobals: failed to allocated dtest", (long) MAXLENG);
    }
    if (!(itest = (long *) malloc((unsigned) (MAXLENG*sizeof(long))))) {
        Error("TestGlobals: failed to allocated itest", (long) MAXLENG);
    }

    for (len=1; len<MAXLENG; len*=2) {
        long ilen = len*sizeof(long);
        long dlen = len*sizeof(double);
        long i;

        if (me == 0) {
            printf("Test length = %ld ... ", len);
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
        tcg_brdcst(itype, (void *) itest, ilen, from);
        tcg_brdcst(dtype, (void *) dtest, dlen, from);

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

        tcg_igop(itype, itest, len, "+");
        tcg_dgop(dtype, dtest, len, "+");

        for (i=0; i<len; i++) {
            long iresult = i*nproc*(nproc-1)/2;
            if (itest[i] != iresult || dtest[i] != (double) iresult){
                printf(" dt %f it %ld ir %ld \n",
                        dtest[i],itest[i],iresult);
                Error("TestGlobals: global sum failed", (long) i);
            }
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

    if (tcg_nodeid() == 0) {
        (void) printf("Hello test ... show network integrity\n----------\n\n");
        (void) fflush(stdout);
    }

    for (node = 0; node<tcg_nnodes(); node++) {
        if (node == tcg_nodeid()) {
            for (kode = 0; kode<tcg_nnodes(); kode++) { 
                (void) sprintf(buf, "Hello to %ld from %ld",
                               kode, tcg_nodeid());
                if (node != kode) {
                    tcg_snd(type, buf, lenbuf, kode, sync);
                }
            }
        }
        else {
            tcg_rcv(type, buf, lenbuf, &lenmes, node, &nodefrom, sync);
            (void) printf("me=%ld, from=%ld: %s\n",
                          tcg_nodeid(), node, buf);
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
        dran = tcg_drand48();
        ran = lo + (long) (dran * (double) (hi-lo+1));
        if (ran < lo) {
            ran = lo;
        }
        if (ran > hi) {
            ran = hi;
        }
        list[i] = ran;
    }
}


/**
 * Stress the system by passing messages between a ranomly selected
 * list of nodes
 */
void Stress()
{
    long me = tcg_nodeid();
    long nproc = tcg_nnodes();
    long type, lenbuf, node, lenmes, nodefrom, i, j, from, to;
    long *list_i, *list_j, *list_n;
    static long len[N_LEN] = {0,1,2,4,8,4095,4096,4097,16367,16368,16369};
    char *buf1, *buf2;
    long n_stress, mod;
    long sync = 1;

    from = 0;
    lenbuf = sizeof(long);

    if (me == 0) {
        (void) printf("Stress test ... randomly exchange messages\n---------");
        (void) printf("\n\nInput no. of messages: ");
        (void) fflush(stdout);
        if (scanf("%ld",&n_stress) != 1) {
            Error("Stress: error reading n_stress",(long) -1);
        }
        if ( (n_stress <= 0) || (n_stress > 100000) ) {
            n_stress = 100;
        }
    }
    type = 13 | MSGINT;
    tcg_brdcst(type, (void *) &n_stress, lenbuf, from);
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
    if ( list_n == (long *) NULL ) {
        Error("Stress: failed to allocate list_n",n_stress);
    }

    if ( (buf1 = malloc((unsigned) 16376)) == (char *) NULL ) {
        Error("Stress: failed to allocate buf1", (long) 16376);
    }

    if ( (buf2 = malloc((unsigned) 16376)) == (char *) NULL ) {
        Error("Stress: failed to allocate buf2", (long) 16376);
    }

    if (me == 0) { /* Make random list of node pairs and message lengths */
        RandList((long) 0, (long) (tcg_nnodes()-1), list_i, n_stress);
        RandList((long) 0, (long) (tcg_nnodes()-1), list_j, n_stress);
        RandList((long) 0, (long) (N_LEN-1), list_n, n_stress);
        for (i=0; i<n_stress; i++)
            list_n[i] = len[list_n[i]];
    }

    node = 0;
    tcg_brdcst(type, (void *) list_i, lenbuf, node);
    type++;
    tcg_brdcst(type, (void *) list_j, lenbuf, node);
    type++;
    tcg_brdcst(type, (void *) list_n, lenbuf, node);
    type++;

    type = 8;

    for (j=0; j<16370; j++) {
        buf1[j] = (char) (j%127);
    }

    j = 0;
    mod = (n_stress-1)/10 + 1;
    for (i=0; i < n_stress; i++) {
        from   = list_i[i];
        to     = list_j[i];
        lenbuf = list_n[i];
        type++;

        if ( (from < 0) || (from >= nproc) ) {
            Error("Stress: from is out of range", from);
        }
        if ( (to < 0) || (to >= nproc) ) {
            Error("Stress: to is out of range", to);
        }

        if (from == to) {
            continue;
        }

        if ( (me == 0) && (j%mod == 0) ) {
            (void) printf("Stress: test=%ld: from=%ld, to=%ld, len=%ld\n",
                    i, from, to, lenbuf);
            (void) fflush(stdout);
        }

        j++;

        if (from == me) {
            tcg_snd(type, buf1, lenbuf, to, sync);
        }
        else if (to == me) {
            (void) bzero(buf2, (int) lenbuf); /* Initialize the rcv buffer */
            buf2[lenbuf] = '+';

            tcg_rcv(type, buf2, lenbuf, &lenmes, from, &nodefrom, sync);

            if (buf2[lenbuf] != '+') {
                Error("Stress: overran buffer on receive",lenbuf);
            }
            if (CheckByte((unsigned char *) buf1, lenbuf) != 
                    CheckByte((unsigned char *) buf2, lenbuf)) {
                Error("Stress: invalid checksum on receive",lenbuf);
            }
            if (lenmes != lenbuf) {
                Error("Stress: invalid message length on receive",lenbuf);
            }
        }
    }

    (void) free(buf2);
    (void) free(buf1);
    (void) free((char *) list_n);
    (void) free((char *) list_j);
    (void) free((char *) list_i);
}


/** Time passing a message round a ring */
void RingTest()
{
    long me = tcg_nodeid();
    long type = 4;
    long left = (me + tcg_nnodes() - 1) % tcg_nnodes();
    long right = (me + 1) % tcg_nnodes();
    char *buf=NULL, *buf2=NULL;
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
        if (scanf("%ld", &max_len) != 1) {
            Error("RingTest: error reading max_len",(long) -1);
        }
        if ( (max_len <= 0) || (max_len >= 4*1024*1024) ) {
            max_len = 256*1024;
        }
    }
    type = 4 | MSGINT;
    tcg_brdcst(type, (void *) &max_len, lenbuf, i);

    if ( (buf = malloc((unsigned) max_len)) == (char *) NULL) {
        Error("failed to allocate buf",max_len);
    }

    if (me == 0) {
        if ( (buf2 = malloc((unsigned) max_len)) == (char *) NULL) {
            Error("failed to allocate buf2",max_len);
        }

        for (i=0; i<max_len; i++) {
            buf[i] = (char) (i%127);
        }
    }

    type = 5;
    lenbuf = 1;
    while (lenbuf <= max_len) {
        int nloops = 10 + 1000/lenbuf;
        int loop = nloops;
        if (me == 0) {
            sum = CheckByte((unsigned char *) buf, lenbuf);
            (void) bzero(buf2, (int) lenbuf);
            start = tcg_time();
            while (loop--) {
                tcg_snd(type, buf, lenbuf, left, sync);
                tcg_rcv(type, buf2, lenbuf, &lenmes, right, &nodefrom, sync);
            }
            used = tcg_time() - start;
            sum2 = CheckByte((unsigned char *) buf2, lenbuf);
            if (used > 0) {
                rate = 1.0e-6 * (double) (tcg_nnodes() * lenbuf) / (double) used;
            }
            else {
                rate = 0.0;
            }
            rate = rate * nloops;
            printf("len=%ld bytes, nloop=%d, used=%8.4f s, rate=%8.4f Mb/s (0x%x, 0x%x)\n",
                    lenbuf, nloops, used, rate, sum, sum2);
            (void) fflush(stdout);
        }
        else {
            while (loop--) {
                tcg_rcv(type, buf, lenbuf, &lenmes, right, &nodefrom, sync);
                tcg_snd(type, buf, lenbuf, left, sync);
            }
        }
        lenbuf *= 2;
    }

    if (me == 0) {
        (void) free(buf2);
    }

    (void) free(buf);
}


/** Test receiveing a message from any node */
void RcvAnyTest()
{
    long me = tcg_nodeid();
    long type = 337 | MSGINT;
    char buf[8];
    long i, j, node, lenbuf, lenmes, nodefrom, receiver, n_msg;
    long sync = 1;

    lenbuf = sizeof(long);

    if (me == 0) {
        (void) printf("RCV any test ... check is working!\n-----------\n\n");
        (void) printf("Input node to receive : ");
        (void) fflush(stdout);
        if (scanf("%ld", &receiver) != 1) {
            Error("RcvAnyTest: error reading receiver",(long) -1);
        }
        if ( (receiver < 0) || (receiver >= tcg_nnodes()) ) {
            receiver = tcg_nnodes()-1;
        }
        (void) printf("Input number of messages : ");
        (void) fflush(stdout);
        if (scanf("%ld", &n_msg) != 1) {
            Error("RcvAnyTest: error reading n_msg",(long) -1);
        }
        if ( (n_msg <= 0) || (n_msg > 10) ) {
            n_msg = 5;
        }
    }

    node = 0;
    tcg_brdcst(type, (void *) &receiver, lenbuf, node);
    type++;
    tcg_brdcst(type, (void *) &n_msg, lenbuf, node);
    type++;

    lenbuf = 0;

    type = 321;
    for (i=0; i<n_msg; i++) {
        if (me == receiver) {
            for (j = 0; j<tcg_nnodes(); j++) {
                if (j !=  me) {
                    node = -1;
                    tcg_rcv(type, buf, lenbuf, &lenmes, node, &nodefrom, sync);
                    (void) printf("RcvAnyTest: message received from %ld\n",
                                  nodefrom);
                    (void) fflush(stdout);
                }
            }
        }
        else {
            tcg_snd(type, buf, lenbuf, receiver, sync);
        }
    }
}


/** Test the load balancing mechanism */
void NextValueTest()
{
    long nproc = tcg_nnodes();
    long me = tcg_nodeid();
    long type = 51 | MSGINT;
    long i, node, lenbuf, n_val, next;
    long ngot, ntimes;
    double start=0, used=0, rate=0;

    lenbuf = sizeof(long);

    if (me == 0) {
        (void) printf("Next value test ... time overhead!\n---------------\n\n");
        (void) printf("Input iteration count : ");
        (void) fflush(stdout);
        if (scanf("%ld", &n_val) != 1) {
            Error("NextValueTest: error reading n_val",(long) -1);
        }
        if ( (n_val < 0) || (n_val >= 10000) ) {
            n_val = 100;
        }
    }
    node = 0;
    tcg_brdcst(type, (void *) &n_val, lenbuf, node);

    /* Loop thru a few values to visually show it is working */

    next = -1;
    for (i=0; i<10; i++) {
        if (i > next) {
            next = tcg_nxtval(nproc);
        }
        sleep(1);
        if (i == next) {
            (void) printf("node %ld got value %ld\n", me, i);
            (void) fflush(stdout);
        }
    }
    nproc = -nproc;
    next = tcg_nxtval(nproc);
    nproc = -nproc;

    /* Now time it for real .. twice*/

    for (ntimes=0; ntimes<2; ntimes++) {
        if (me == 0) {
            start = tcg_time();
        }

        next = -1;
        ngot = 0;

        for (i=0; i<n_val; i++) {
            if (i > next) {
                next = tcg_nxtval(nproc);
            }
            if (i == next) {
                ngot++;
            }
        }

        nproc = -nproc;
        next = tcg_nxtval(nproc);
        nproc = -nproc;

        if (me == 0) {
            used =  tcg_time() - start;
            rate = ngot ? used / ngot: 0.;
            printf("node 0: From %ld busy iters did %ld, used=%lfs per call\n",
                    n_val, ngot, rate);
            fflush(stdout);
        }

        type++;
        tcg_synch(type);
    }
}


void ToggleDebug()
{
    static long on = 0;
    long me = tcg_nodeid();
    long type = 666 | MSGINT;
    long lenbuf = sizeof(long);
    long from=0;
    long node;

    if (me == 0) {
        (void) printf("\nInput node to debug (-1 = all) : ");
        (void) fflush(stdout);
        if (scanf("%ld", &node) != 1) {
            Error("ToggleDebug: error reading node",(long) -1);
        }
    }
    tcg_brdcst(type, (void *) &node, lenbuf, from);

    if ((node < 0) || (node == me)) {
        on = (on + 1)%2;
        tcg_setdbg(on);
    }
}


int main(int argc, char **argv)
{
    long type;
    long lenbuf;
    long node, opt;

    tcg_alt_pbegin(&argc, &argv);

    (void) printf("In process %ld\n", tcg_nodeid());
    (void) fflush(stdout);

    /* Read user input for action */

    lenbuf = sizeof(long);
    node = 0;

    while (1) {

        (void) fflush(stdout);
        if (tcg_nodeid() == 0) {
            (void) sleep(1);
        }
        type = 999;
        tcg_synch(type);
        (void) sleep(1);

        if (tcg_nodeid() == 0) {
again:
            (void) printf("\n\
                          0=quit\n\
                          1=Ring             5=NxtVal\n\
                          2=Stress           6=Global\n\
                          3=Hello            7=Debug\n\
                          4=RcvAny           8=Probe\n\n\
                          Enter test number : ");

            (void) fflush(stdout);
            if (scanf("%ld", &opt) != 1) {
                Error("test: input of option failed",(long) -1);
            }
            (void) printf("\n");
            (void) fflush(stdout);
            if ( (opt < 0) || (opt > 8) ) {
                goto again;
            }
        }
        type = 2 | MSGINT;
        tcg_brdcst(type, (void *) &opt, lenbuf, node);

        switch (opt) {
            case 0:
                if (tcg_nodeid() == 0) {
                    tcg_stats();
                }
                tcg_pend();
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

            default:
                Error("test: invalid option", opt);
                break;
        }
    }
}
