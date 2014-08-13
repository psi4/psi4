#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif

#include <pvm3.h>

#define MAX_PROC 128

#define TCGTIME_ TCGTIME

double  TCGTIME_();

#define ENCODING PvmDataRaw


long NODEID_()
{
    return((long) pvm_get_PE(pvm_mytid()));
}


long NNODES_()
{
    return((long)pvm_gsize(0));
}


/**
 * Error handler
 */
void Error(string, code)
    char *string;
    long code;
{

    (void) fflush(stdout);
    (void) fflush(stderr);

    (void) fprintf(stderr, "%3d:%s %ld(%x)\n", NODEID_(), string, code, code);
    (void) perror("system message");

    (void) fflush(stdout);
    (void) fflush(stderr);

    globalexit(1);
}


void SND_(type, buf, lenbuf, node, sync)
    long *type;
    char *buf;
    long *lenbuf;
    long *node;
    long *sync;
{
    long tid=pvm_gettid("", *node);

    pvm_psend(tid, *type, buf, *lenbuf, PVM_BYTE); 
}


void RCV_(type, buf, lenbuf, lenmes, nodeselect, nodefrom, sync)
    long *type;
    char *buf;
    long *lenbuf;
    long *lenmes;
    long *nodeselect;
    long *nodefrom;
    long *sync;
{
    int tid=*nodeselect, tidfrom;

    if(tid >-1) tid=pvm_gettid("", *nodeselect);
    pvm_precv(tid, *type, buf, *lenbuf, PVM_BYTE, &tidfrom, 0, 0);
    *nodefrom = pvm_get_PE(tidfrom);
}


/**
 * Time passing a message round a ring
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
    char *malloc();

    i = 0;
    lenbuf = sizeof(long);

    if (me == 0) {
        (void) printf("Ring test ... time network performance\n---------\n\n");
        /*
           (void) printf("Input maximum message size: ");
           (void) fflush(stdout);
           if (scanf("%ld", &max_len) != 1)
           Error("RingTest: error reading max_len",(long) -1);
           if ( (max_len <= 0) || (max_len >= 4*1024*1024) )
           max_len = 256*1024;
           */
    }
    max_len = 512*1024;
    /*  type = 4 | MSGINT;*/
    /*  BRDCST_(&type, (char *) &max_len, &lenbuf, &i);*/

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


int main(int argc, char **argv)
{
    long me, nproc;
    long buf[10];
    long type=1, len=sizeof(buf), node, sync=1, i;


    me = NODEID_();
    nproc = NNODES_();
    /*
       if(me==0){
       for(i=0;i<10;i++)buf[i]=i;

       node=1;
       SND_(&type, buf, &len, &node, &sync);
       printf("me=%d nproc = %d\n", me, nproc);
       fflush(stdout);
       }else{
       node=0;
       RCV_(&type, buf, &len,  &len, &node, &node, &sync);
       printf("me=%d nproc = %d\n", me, nproc);
       fflush(stdout);
       for(i=0;i<10;i++)printf("%d ",buf[i]);;
       }
       */
    RingTest();
}
