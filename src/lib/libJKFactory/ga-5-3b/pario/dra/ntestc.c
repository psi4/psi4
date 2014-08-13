#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#define BASE_NAME "dra.file"
#ifdef  HPIODIR
#   define FNAME HPIODIR/*BASE_NAME*/
#else
#   define FNAME BASE_NAME
#endif

#define NDIM 3
#define SIZE 20
#define NSIZE 8000
#define LSIZE 216000
#define TRUE (logical)1
#define FALSE (logical)0
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

#include "dra.h"
#include "ga.h"
#include "macdecls.h"
#include "mp3.h"

#define MAXDIM GA_MAX_DIM

#define COPY_INT_DRAT(src,dst) {\
int i; \
for (i=0; i<NDIM; ++i) { \
    dst[i] = (dra_size_t)src[i]; \
} \
}

void init_char(char* str, Integer length, char a)
{
    Integer i;
    for (i=0; i<length; i++) str[i] = a;
}


void test_io_int()
{
    int n, nsize;
    int a[NSIZE];
    int g_a, g_b, d_a;
    int plus, minus;
    int i, j, index, req, err, type, op, ndim=NDIM;
    int me, nproc;
    int ga_dims[MAXDIM], ga_chunk[MAXDIM], lo[MAXDIM], hi[MAXDIM], ld[MAXDIM];
    dra_size_t dims[MAXDIM], reqdim[MAXDIM], chunk[MAXDIM];
    char filename[200], fname[200];
    char name[80];

    nproc = GA_Nnodes();
    me    = GA_Nodeid();
    n = pow(NSIZE,1.0/ndim)+0.5;
    nsize = pow(n,ndim);

    init_char(name,100, ' ');
    init_char(filename,200, ' ');

    /* a() is a local copy of what the l array should start as */

    for (i=0; i<nsize; i++) {
        a[i] = i;
    }

    if (me == 0) printf("Creating global arrays\n");
    GA_Sync();
    for (i=0; i<ndim; i++) {
        dims[i] = n;
        ga_dims[i] = n;
        chunk[i] = 1;
        ga_chunk[i] = 1;
        ld[i] = n;
    }
    type = MT_INT;
    g_a = NGA_Create(type, ndim, ga_dims, "a", ga_chunk);
    if (!g_a) GA_Error("NGA_Create failed: a", 0);
    g_b = NGA_Create(type, ndim, ga_dims, "b", ga_chunk);
    if (!g_b) GA_Error("NGA_Create failed: b", 0);

    for (j=me; j<n; j += nproc) {
        lo[0] = j;
        hi[0] = j;
        for (i=1; i<ndim; i++) {
            lo[i] = 0;
            hi[i] = n-1;
        }
        index = j;
        for (i=0; i<ndim-1;i++) index *= SIZE;
        NGA_Put(g_a, lo, hi, &a[index], ld);
    }

    if (me == 0) {
        printf("Creating Disk Array %d",n);
        for (i=1; i<ndim; i++) printf(" x %d",n);
        printf("\n");
    }
    dims[0] = n;
    reqdim[0] = 1;
    for (i=1; i<ndim; i++) {
        dims[i] = n;
        reqdim[i] = n;
    }
    type = MT_INT;
    op = DRA_RW;
    /*bjp sprintf(fname,"%s%d",FNAME,me); */
    sprintf(fname,"%s%d",FNAME,0);
    if (NDRA_Create(type, ndim, dims, "array A", 
                fname, op, reqdim, &d_a) != 0)
        GA_Error("NDRA_Create failed: ",0);
    if (me == 0) printf("OK\n\n");

    if (me == 0) printf("Writing Global Array to Disk Array\n");
    if (NDRA_Write(g_a, d_a, &req) != 0)
        GA_Error("NDRA_Write failed:",0);
    if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed: ",req);
    if (me == 0) printf("OK\n\n");
    if (me == 0) printf("Closing Disk Array\n");
    if (DRA_Close(d_a) != 0) GA_Error("DRA_Close failed: ",d_a);
    if (me == 0) printf("OK\n\n");

    if (me == 0) printf("Opening Existing Disk Array\n");
    op = DRA_R;
    if (DRA_Open(fname ,op, &d_a) != 0)
        GA_Error("DRA_Open failed",0);

    if (NDRA_Inquire(d_a, &type, &ndim, dims, name, filename) != 0)
        GA_Error("NDRA_Inquire failed",0);
    for (i=0; i<ndim; i++) {
        if (dims[i] != n) GA_Error("dims error",dims[i]);
    }
    if (type != MT_INT) GA_Error("type error",type);
    if (me == 0) printf("array name read from disk is: %s\n",name);
    GA_Sync();
    if (me == 0) printf("OK\n\n");

    if (me == 0) printf("Checking NDRA_Read\n");
    if (NDRA_Read(g_b, d_a, &req) != 0)
        GA_Error("NDRA_Read failed:",0);
    fflush(stdout);
    if(DRA_Wait(req) != 0) GA_Error("DRA_Wait failed: ",req);

    /*     error checking: (g_a - g_b)^2 */
    plus = 1;
    minus = -1;
    GA_Add(&plus, g_a, &minus, g_b, g_b);
    err = GA_Idot(g_b, g_b); 

    if (err != 0) {
        if( me == 0) GA_Error("NDRA_Read comparison failed", err);
    } else {
        if (me == 0) printf("OK\n");
    }
    if (me == 0) printf("\n");

    if (me == 0) printf("Checking DRA_Delete\n");
    if (DRA_Delete(d_a) != 0) GA_Error("DRA_Delete failed",0);
    if (me == 0) printf("OK\n\n");
    GA_Destroy(g_a);
    GA_Destroy(g_b);
}


float ran0(long *idum)
{
    long k;
    float ans;

    *idum ^= MASK;
    k=(*idum)/IQ;
    *idum = IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    ans=AM*(*idum);
    *idum ^= MASK;
    if(ans<0)return -ans;
    else return ans;
}


int iran(int i, long *idum)
{
    float fr = ran0(idum);
    if(fr<0)printf("problem %f\n",fr);
#if 1
    return (int) (i * fr);
#else
    return (int) ((i-1) * fr)+1;
#endif
}


void swap(int *i,int *j)
{
    int a;
    a = *i;
    *i = *j;
    *j = a;
    return;
}


void test_io_dbl()
{
    int n, m, nsize, msize;
    double a[NSIZE],  err, plus, minus;
    int g_a, g_b, d_a;
    int i, j, req, loop, index, ndim=NDIM;
    int dlo[NDIM],dhi[NDIM];
    int glo[NDIM],ghi[NDIM];
    int elem, type, op;
    int me, nproc;
    int ga_dims[NDIM], ga_chunk[NDIM], ld[NDIM];
    dra_size_t dims[NDIM], reqdim[NDIM], chunk[NDIM];
    dra_size_t dlo_d[NDIM], dhi_d[NDIM];
    char fname[200];
    long idum;

    loop  = 30;
    nproc = GA_Nnodes();
    me    = GA_Nodeid();
    n = pow(NSIZE,1.0/ndim)+0.5;
    m = pow(LSIZE,1.0/ndim)+0.5;
    nsize = pow(n,ndim);
    msize = pow(m,ndim);

    /*     a() is a local copy of what the l array should start as */

    for (i=0; i<nsize; i++) {
        a[i] = i;
    }

    GA_Sync();
    for (i=0; i<ndim; i++) {
        dims[i] = n;
        ga_dims[i] = n;
        chunk[i] = 1;
        ga_chunk[i] = 1;
        ld[i] = n;
    }
    if (me == 0) printf("Creating global arrays\n\n");
    g_a = NGA_Create(MT_DBL, ndim, ga_dims, "a", ga_chunk);
    if (g_a == 0) GA_Error("GA_Create failed: a", 0);
    g_b = NGA_Create(MT_DBL, ndim, ga_dims, "b", ga_chunk);
    if (g_b == 0) GA_Error("GA_Create failed: b", 0);

    if (me == 0) printf("Zeroing global arrays\n\n");
    GA_Zero(g_a);
    GA_Zero(g_b);

    for (j=me; j<n; j += nproc) {
        dlo[0] = j;
        dhi[0] = j;
        for (i=1; i<ndim; i++) {
            dlo[i] = 0;
            dhi[i] = n-1;
        }
        index = j;
        for (i=0; i<ndim-1;i++) index *= SIZE;
        NGA_Put(g_a, dlo, dhi, &a[index], ld);
    }

    if (me == 0) {
        printf("Creating Disk Array %d",n);
        for (i=1; i<ndim; i++) printf(" x %d",n);
        printf("\n");
    }
    reqdim[0] = 3;
    for (i=1; i<ndim; i++) reqdim[i] = n;

    type = MT_DBL;
    op = DRA_RW;
    /*bjp sprintf(fname,"%s%d",FNAME,me);*/
    sprintf(fname,"%s%d",FNAME,1);
    if (NDRA_Create(type, ndim, dims, "A", 
                fname, op, reqdim, &d_a) != 0)
        GA_Error("NDRA_Create failed: ",0);

    if (me == 0) printf("Writing Global Array to Disk Array\n");
    if (NDRA_Write(g_a, d_a, &req) != 0)
        GA_Error("NDRA_Write failed:",0);
    if (me == 0) printf("\n");
    if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed: ",req);

    if (DRA_Close(d_a) != 0) GA_Error("DRA_Close failed: ",d_a);

    if (me == 0) printf("Checking NDRA_Read\n");
    op = DRA_R;
    if (DRA_Open(fname, op, &d_a) != 0) GA_Error("DRA_Open failed",0);
    if (NDRA_Read(g_b, d_a, &req) != 0) GA_Error("NDRA_Read failed:",0);
    if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed: ",req);

    /*     error checking: (g_a - g_b)^2 */

    plus = 1.0;
    minus = -1.0;
    GA_Add(&plus, g_a, &minus, g_b, g_b);
    err = GA_Ddot(g_b, g_b);
    if (err != 0) {
        if (me == 0) printf("NDRA_Read error = %f\n\n",err);
    } else {
        if (me == 0) printf("OK\n\n");
    }

    if (me == 0) printf("Checking NDRA_Read_section\n");

    idum = 7889;
    ran0(&idum);
    GA_Zero(g_b);
    for (j=0; j<loop; j++) {
        for (i=0; i<ndim; i++) {
            dlo[i] = iran(n,&idum);
            dhi[i] = iran(n,&idum);
            if (dlo[i] > dhi[i]) swap(&dlo[i],&dhi[i]);
            elem = dhi[i] - dlo[i] + 1;
            glo[i] = iran(n-elem,&idum);
            ghi[i] = glo[i] + elem - 1;
        }

        if (me ==  0) {
            printf(" read global[");
            printf("%4d:%4d",glo[0],ghi[0]);
            for (i=1; i<ndim; i++) printf(",%4d:%4d",glo[i],ghi[i]);
            printf("] from disk[");
            printf("%4d:%4d",dlo[0],dhi[0]);
            for (i=1; i<ndim; i++) printf(",%4d:%4d",dlo[i],dhi[i]);
            printf("]\n");
            fflush(stdout);
        }

        COPY_INT_DRAT(dlo,dlo_d);
        COPY_INT_DRAT(dhi,dhi_d);
        if (NDRA_Read_section(FALSE, g_b, glo, ghi, d_a, dlo_d, dhi_d, &req) != 0)
            GA_Error("NDRA_Read_section failed:",0);
        if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed:",req);

        NGA_Add_patch(&plus, g_a, dlo, dhi, &minus, g_b, glo, ghi,
                g_b, glo, ghi);
        err = NGA_Ddot_patch(g_b, 'n', glo, ghi, g_b, 'n', glo, ghi);
        if (err !=0.0 && me == 0) {
            printf("error = %f\n",err);
            GA_Error("failed",0);
        }
    }
    if (me == 0) printf("OK\n\n");
    if (DRA_Delete(d_a) != 0)
        GA_Error("DRA_Delete failed",0);

    /*  now d_a is 2^NDIM times larger than g_a */

    if (me == 0) {
        printf("Creating New Disk Array %d",m);
        for (i=1; i<ndim; i++) printf(" x %d",m);
        printf("\n");
    }
    dims[0] = m;
    reqdim[0] = 2;
    for (i=1; i<ndim; i++) {
        dims[i] = m;
        reqdim[i] = n;
    }
    type = MT_DBL;
    op = DRA_RW;
    sprintf(fname,"%s%d",FNAME,me);
    if (NDRA_Create(type, ndim, dims, "A", 
                fname, op, reqdim, &d_a) != 0)
        GA_Error("NDRA_Create failed: ",0);
    if (me == 0) printf("OK\n\n");

    if (me == 0) printf("Testing NDRA_Write_section\n");
    for (j=0; j<loop; j++) {
        for (i=0; i<ndim; i++) {
            glo[i] = iran(n,&idum);
            if (glo[i] > ghi[i]) swap(&glo[i],&ghi[i]);
            elem = ghi[i] - glo[i] +1;
            dlo[i] = iran(m-elem,&idum);
            dhi[i] = dlo[i]+elem-1;
        }

        if (me ==  0) {
            printf(" writing global[");
            printf("%4d:%4d",glo[0],ghi[0]);
            for (i=1; i<ndim; i++) printf(",%4d:%4d",glo[i],ghi[i]);
            printf("] to disk[");
            printf("%4d:%4d",dlo[0],dhi[0]);
            for (i=1; i<ndim; i++) printf(",%4d:%4d",dlo[i],dhi[i]);
            printf("]\n");
            fflush(stdout);
        }

        COPY_INT_DRAT(dlo,dlo_d);
        COPY_INT_DRAT(dhi,dhi_d);
        if (NDRA_Write_section(FALSE, g_a, glo, ghi, d_a, dlo_d, dhi_d, &req) != 0)
            GA_Error("NDRA_Write_section failed:",0);
        if(DRA_Wait(req) != 0) GA_Error("DRA_Wait failed:",req);

        /*     NDRA_Read was tested already and we use it for testing
               NDRA_Write_section */

        COPY_INT_DRAT(dlo,dlo_d);
        COPY_INT_DRAT(dhi,dhi_d);
        if (NDRA_Read_section(FALSE, g_b, glo, ghi, d_a, dlo_d, dhi_d, &req) != 0)
            GA_Error("NDRA_Read_section failed:",0);
        if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed:",req);

        NGA_Add_patch(&plus, g_a, glo, ghi, &minus, g_b, glo, ghi,
                g_b, glo, ghi);
        err = NGA_Ddot_patch(g_b, 'n', glo, ghi, g_b, 'n', glo, ghi);
        if (err != 0.0 && me == 0) {
            printf("error = %f",err);
            GA_Error("error in NDRA_Write_section",0);
        }
    }
    if (me == 0) printf("OK\n");

    if (DRA_Delete(d_a) != 0) GA_Error("DRA_Delete failed",0);
    GA_Destroy(g_a);
    GA_Destroy(g_b);
}


int main(int argc, char **argv)
{
    int status, me;
    int max_arrays = 10;
    int stack=80000, heap=80000;
    double max_sz=100000000.0, max_disk=10000000000.0, max_mem=1000000.0;

    MP_INIT(argc,argv);
    if(MA_init(MT_DBL, stack, heap) ) {
        GA_Initialize();
        me    = GA_Nodeid();
        if (DRA_Init(max_arrays, max_sz, max_disk, max_mem) != 0)
            GA_Error("DRA_Init failed: ",0);
        if (me == 0) printf("\n  TESTING INTEGERS\n\n");
        test_io_int();
        if (me == 0) printf("\n  TESTING DOUBLES\n\n");
        test_io_dbl();
        status = DRA_Terminate();
        GA_Terminate();
    } else {
        printf("MA_init failed\n");
    }
    MP_FINALIZE();

    return 0;
}
