#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "srftoc.h"
#include "tcgmsgP.h"

#define BUF_SIZE  10000
#define IBUF_SIZE (BUF_SIZE * sizeof(double)/sizeof(long)) 
double _gops_work[BUF_SIZE];

long one=1;

#define TCG_MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define TCG_MIN(a,b) (((a) <= (b)) ? (a) : (b))
#define TCG_ABS(a)   (((a) >= 0) ? (a) : (-(a)))


void BRDCST_(long *type, void *buf, long *len, long *originator)
{
    long me=NODEID_(), nproc=NNODES_(), lenmes, from, root=0;
    long up, left, right;

    /* determine location in the binary tree */
    up    = (me-1)/2;    if(up >= nproc)       up = -1;
    left  =  2* me + 1;  if(left >= nproc)   left = -1;
    right =  2* me + 2;  if(right >= nproc) right = -1;

    /*  originator sends data to root */
    if (*originator != root ){
        if(me == *originator) SND_(type, buf, len, &root, &one); 
        if(me == root) RCV_(type, buf, len, &lenmes, originator, &from, &one); 
    }

    if (me != root) RCV_(type, buf, len, &lenmes, &up, &from, &one);
    if (left > -1)  SND_(type, buf, len, &left, &one);
    if (right > -1) SND_(type, buf, len, &right, &one);
}


/**
 * implements x = op(x,work) for integer datatype
 *  x[n], work[n] -  arrays of n integers
 */ 
static void idoop(long n, char *op, long *x, long *work)
{
    if (strncmp(op,"+",1) == 0)
        while(n--)
            *x++ += *work++;
    else if (strncmp(op,"*",1) == 0)
        while(n--)
            *x++ *= *work++;
    else if (strncmp(op,"max",3) == 0)
        while(n--) {
            *x = TCG_MAX(*x, *work);
            x++; work++;
        }
    else if (strncmp(op,"min",3) == 0)
        while(n--) {
            *x = TCG_MIN(*x, *work);
            x++; work++;
        }
    else if (strncmp(op,"absmax",6) == 0)
        while(n--) {
            register long x1 = TCG_ABS(*x), x2 = TCG_ABS(*work);
            *x = TCG_MAX(x1, x2);
            x++; work++;
        }
    else if (strncmp(op,"absmin",6) == 0)
        while(n--) {
            register long x1 = TCG_ABS(*x), x2 = TCG_ABS(*work);
            *x = TCG_MIN(x1, x2);
            x++; work++;
        }
    else if (strncmp(op,"or",2) == 0)
        while(n--) {
            *x |= *work;
            x++; work++;
        }
    else
        Error("idoop: unknown operation requested", (long) n);
}


/**
 * implements x = op(x,work) for double datatype
 *  x[n], work[n] -  arrays of n doubles
 */ 
static void ddoop(long n, char *op, double *x, double *work)
{
    if (strncmp(op,"+",1) == 0)
        while(n--)
            *x++ += *work++;
    else if (strncmp(op,"*",1) == 0)
        while(n--)
            *x++ *= *work++;
    else if (strncmp(op,"max",3) == 0)
        while(n--) {
            *x = TCG_MAX(*x, *work);
            x++; work++;
        }
    else if (strncmp(op,"min",3) == 0)
        while(n--) {
            *x = TCG_MIN(*x, *work);
            x++; work++;
        }
    else if (strncmp(op,"absmax",6) == 0)
        while(n--) {
            register double x1 = TCG_ABS(*x), x2 = TCG_ABS(*work);
            *x = TCG_MAX(x1, x2);
            x++; work++;
        }
    else if (strncmp(op,"absmin",6) == 0)
        while(n--) {
            register double x1 = TCG_ABS(*x), x2 = TCG_ABS(*work);
            *x = TCG_MIN(x1, x2);
            x++; work++;
        }
    else
        Error("ddoop: unknown operation requested", (long) n);
}


void DGOP_(
        long *type, double *x, long *n, char *op, int oplen)
{
    long me=NODEID_(), nproc=NNODES_(), len, lenmes, from, root=0;
    double *work = _gops_work, *origx = x;
    long ndo, up, left, right, np=*n, orign = *n;

    /* determine location in the binary tree */
    up    = (me-1)/2;    if(up >= nproc)       up = -1;
    left  =  2* me + 1;  if(left >= nproc)   left = -1;
    right =  2* me + 2;  if(right >= nproc) right = -1;

    while ((ndo = (np <= BUF_SIZE) ? np : BUF_SIZE)) {
        len = lenmes = ndo*sizeof(double);

        if (left > -1) {
            RCV_(type, (char *) work, &len, &lenmes, &left, &from, &one);
            ddoop(ndo, op, x, work);
        }
        if (right > -1) {
            RCV_(type, (char *) work, &len, &lenmes, &right, &from, &one);
            ddoop(ndo, op, x, work);
        }
        if (me != root) SND_(type, x, &len, &up, &one); 

        np -=ndo;
        x  +=ndo;
    }

    /* Now, root broadcasts the result down the binary tree */
    len = orign*sizeof(double);
    BRDCST_(type, (char *) origx, &len, &root);
}


void IGOP_(long *type, long *x, long *n, char *op, int oplen)
{
    long me=NODEID_(), nproc=NNODES_(), len, lenmes, from, root=0;
    long *work = (long*)_gops_work;
    long *origx = x;
    long ndo, up, left, right, np=*n, orign =*n;

    /* determine location in the binary tree */
    up    = (me-1)/2;    if(up >= nproc)       up = -1;
    left  =  2* me + 1;  if(left >= nproc)   left = -1;
    right =  2* me + 2;  if(right >= nproc) right = -1;

    while ((ndo = (np<=IBUF_SIZE) ? np : IBUF_SIZE)) {
        len = lenmes = ndo*sizeof(long);

        if (left > -1) {
            RCV_(type, (char *) work, &len, &lenmes, &left, &from, &one);
            idoop(ndo, op, x, work);
        }
        if (right > -1) {
            RCV_(type, (char *) work, &len, &lenmes, &right, &from, &one);
            idoop(ndo, op, x, work);
        }
        if (me != root) SND_(type, x, &len, &up, &one); 

        np -=ndo;
        x  +=ndo;
    }

    /* Now, root broadcasts the result down the binary tree */
    len = orign*sizeof(long);
    BRDCST_(type, (char *) origx, &len, &root);
}
