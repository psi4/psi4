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
#if HAVE_STRING_H
#   include <string.h>
#endif

#define FNAME      "/scratch/da.try"
#define FNAME_ALT  "/tmp/da.try"
#define FNAME1     "/scratch/da1.try"
#define FNAME1_ALT "/tmp/da1.try"
#define FNAME2     "/scratch/da2.try"
#define FNAME2_ALT "/tmp/da2.try"

#include "dra.h"
#include "ga.h"
#include "macdecls.h"
#include "mp3.h"

#ifndef MAXDIM
#   define MAXDIM GA_MAX_DIM
#endif
#ifndef TRUE
#   define TRUE (Logical)1
#endif
#ifndef FALSE
#   define FALSE (Logical)0
#endif

/* If USER_CONFIG set to:
 *   0: Use default configuration
 *   1: Number of files is 1, number of I/O procs equals
 *      the number of nodes
 *   2: Number of files and number of I/O procs equals
 *      the number of nodes
 *   3: Number of files and number of I/O procs equals 1
 * These configurations only apply if files are not located
 * on local scratch disk. For local scratch, system defaults
 * to 1 I/O proc per node and 1 file per node. Note that
 * USER_CONFIG=1 will only work if system has parallel I/O.
 *
 * Memory and Disk Usage:
 * The aggregate memory required to run this test is approximately
 * 2*SIZE**NDIM*sizeof(double) bytes. The amount of disk space
 * required is approximately 1+2**NDIM times this amount.
 */
#define USER_CONFIG 0
#define TEST_TRANSPOSE 0

#define NDIM 3
#define SIZE 250

/*
#define NDIM 2
#define SIZE 4000

#define NDIM 1
#define SIZE 16000000
*/

#define SWITCH 0
/*#define MAXDIM 7*/
/*#define TRUE (logical)1*/
/*#define FALSE (logical)0*/

#define MULTFILES 0

#ifdef SOLARIS
#   if MULTFILES
#       define USEMULTFILES 1
#   endif
#else
#  define USEMULTFILES 0
#endif

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876


void filename_check(char *result, const char *fname, const char *fname_alt)
{
    FILE *fd;
    strcpy(result, fname);
    if (! (fd = fopen(result, "w"))) {
        strcpy(result, fname_alt);
        if (! (fd = fopen(result, "w"))) {
            GA_Error("Could not open file", 0);
        }
    }
    fclose(fd);
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
    return ans;
}


void fill_random(double *a, int isize)
{
    long *idum;
    long i, j;

    j = 38282;
    idum = &j;
    a[0] = (double)ran0(idum);
    for (i=0; i<(long)isize; i++) {
        a[i] = (double)(10000.0*ran0(idum));
    }
}


void test_io_dbl()
{
    int ndim = NDIM;
    double err, tt0, tt1, mbytes;
    int g_a, g_b, d_a, d_b;
#if TEST_TRANSPOSE
    int g_c, g_d, d_c;
#endif
    int i, req, loop;
    dra_size_t dlo[MAXDIM],dhi[MAXDIM];
    dra_size_t n, m;
    dra_size_t ddims[MAXDIM], reqdims[MAXDIM];
    int glo[MAXDIM],ghi[MAXDIM];
    int dims[MAXDIM];
    int me, nproc, isize, numfiles, nioprocs;
    double plus, minus;
    double *index;
    int ld[MAXDIM], chunk[MAXDIM];
#if USEMULTFILES
    int ilen;
#endif
    char filename[80], filename1[80];

    n = SIZE;
    m = 2*SIZE;

    loop  = 30;
    req = -1;
    nproc = GA_Nnodes();
    me    = GA_Nodeid();
    nioprocs = GA_Cluster_nnodes();
    numfiles = nioprocs;

    if (me == 0) {
        printf("Creating temporary global arrays %ld",(long)n);
        for (i=1; i<ndim; i++) {
            printf(" x %ld",(long)n);
        }
        printf("\n");
    }
    if (me == 0) fflush(stdout);
    GA_Sync();
    for (i=0; i<ndim; i++) {
        dims[i] = n;
        chunk[i] = 1;
    }

    g_a = NGA_Create(MT_DBL, ndim, dims, "a", chunk);
    if (!g_a) GA_Error("NGA_Create failed: a", 0);
    g_b = NGA_Create(MT_DBL, ndim, dims, "b", chunk);
    if (!g_b) GA_Error("NGA_Create failed: b", 0);
    if (me == 0) printf("done\n");
    if (me == 0) fflush(stdout);

    /*     initialize g_a, g_b with random values
           ... use ga_access to avoid allocating local buffers for ga_put */

    GA_Sync();
    NGA_Distribution(g_a, me, glo, ghi);
    NGA_Access(g_a, glo, ghi, &index, ld);
    isize = 1;
    for (i=0; i<ndim; i++) isize *= (ghi[i]-glo[i]+1);
    fill_random(index, isize);
    GA_Sync();
    GA_Zero(g_b);

    /*.......................................................................*/
#if SWITCH
    if (me == 0) {
        printf("Creating Disk array %d",m);
        for (i=1; i<ndim; i++) {
            printf(" x %d",m);
        }
        printf("\n");
    }
    if (me == 0) fflush(stdout);
    for (i=0; i<ndim; i++) {
        ddims[i] = m;
        reqdims[i] = n;
    }
    GA_Sync();
    filename_check(filename1, FNAME1, FNAME1_ALT);
    if (NDRA_Create(MT_DBL, ndim, ddims, "B", filename1, DRA_RW,
                reqdims, &d_b) != 0) {
        GA_Error("NDRA_Create failed(d_b): ",0);
    }

    if (me == 0) printf("non alligned blocking write\n");
    if (me == 0) fflush(stdout);

    for (i=0; i<ndim; i++) {
        glo[i] = 0;
        ghi[i] = n-1;
        dlo[i] = 1;
        dhi[i] = n;
    }
    GA_Sync();
    tt0 = MP_TIMER();
    if (NDRA_Write_section(FALSE, g_a, glo, ghi,
                d_b, dlo, dhi, &req) != 0)
        GA_Error("ndra_write_section failed:",0);

    if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed(d_b): ",req);
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    mbytes = 1.e-6*(double)(pow(n,ndim)*sizeof(double));
    if (me == 0) {
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }

    if (DRA_Close(d_b) != 0) GA_Error("DRA_Close failed(d_b): ",d_b);
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    if (me == 0) {
        printf("Time including DRA_Close\n");
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }

    if (me == 0) printf("\n");
    if (me == 0) printf("disk array closed\n");
    if (me == 0) fflush(stdout);

    GA_Sync();
    if (me == 0) {
        printf("Creating Disk array %d",n);
        for (i=1; i<ndim; i++) {
            printf(" x %d",n);
        }
        printf("\n");
    }
    for (i=0; i<ndim; i++) {
        ddims[i] = n;
        reqdims[i] = n;
    }
    filename_check(filename, FNAME, FNAME_ALT);
    if (NDRA_Create(MT_DBL, ndim, ddims, "A", filename, DRA_RW,
                reqdims, &d_a) != 0)
    {
        GA_Error("NDRA_Create failed(d_a): ",0);
    }
    if (me == 0) printf("alligned blocking write\n");
    fflush(stdout);
    tt0 = MP_TIMER();
    if (NDRA_Write(g_a, d_a, &req) != 0) GA_Error("NDRA_Write failed(d_a):",0);
    if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed(d_a): ",req);
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    mbytes = 1.e-6 * (double)(pow(n,ndim)*sizeof(double));
    if (me == 0) {
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }

    if (DRA_Close(d_a) != 0) GA_Error("DRA_Close failed(d_a): ",d_a);
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    if (me == 0) {
        printf("Time including DRA_Close\n");
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }

    if (me == 0) printf("\n");
    if (me == 0) printf("disk array closed\n");
    if (me == 0) fflush(stdout);
#else /* SWITCH */
    if (me == 0) {
        printf("Creating Disk array %ld",(long)n);
        for (i=1; i<ndim; i++) {
            printf(" x %ld",(long)n);
        }
        printf("\n");
    }
    if (me == 0) fflush(stdout);
    for (i=0; i<ndim; i++) {
        ddims[i] = n;
        reqdims[i] = n;
    }
    GA_Sync();
    filename_check(filename, FNAME, FNAME_ALT);
    if (NDRA_Create(MT_DBL, ndim, ddims, "A", filename, DRA_RW,
                reqdims, &d_a) != 0) {
        GA_Error("NDRA_Create failed(d_a): ",0);
    }
    if (me == 0) printf("alligned blocking write\n");
    fflush(stdout);
    tt0 = MP_TIMER();
    if (NDRA_Write(g_a, d_a, &req) != 0) GA_Error("NDRA_Write failed(d_a):",0);
    if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed(d_a): ",req);
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    mbytes = 1.e-6 * (double)(pow(n,ndim)*sizeof(double));
    if (me == 0) {
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }

    if (DRA_Close(d_a) != 0) GA_Error("DRA_Close failed(d_a): ",d_a);
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    if (me == 0) {
        printf("Time including DRA_Close\n");
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }

    if (me == 0) printf("\n");
    if (me == 0) printf("disk array closed\n");
    if (me == 0) fflush(stdout);
    /*.......................................................................*/

    if (me == 0) {
        printf("Creating Disk array %ld",(long)m);
        for (i=1; i<ndim; i++) {
            printf(" x %ld",(long)m);
        }
        printf("\n");
    }
    for (i=0; i<ndim; i++) {
        ddims[i] = m;
        reqdims[i] = n;
    }
    filename_check(filename1, FNAME1, FNAME1_ALT);
    if (NDRA_Create(MT_DBL, ndim, ddims, "B", filename1, DRA_RW,
                reqdims, &d_b) != 0) {
        GA_Error("NDRA_Create failed(d_b): ",0);
    }

    if (me == 0) printf("non alligned blocking write\n");
    if (me == 0) fflush(stdout);

    for (i=0; i<ndim; i++) {
        glo[i] = 0;
        ghi[i] = n-1;
        dlo[i] = 1;
        dhi[i] = n;
    }
    tt0 = MP_TIMER();
    if (NDRA_Write_section(FALSE, g_a, glo, ghi,
                d_b, dlo, dhi, &req) != 0)
        GA_Error("ndra_write_section failed:",0);

    if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed(d_b): ",req);
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    mbytes = 1.e-6*(double)(pow(n,ndim)*sizeof(double));
    if (me == 0) {
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }

    if (DRA_Close(d_b) != 0) GA_Error("DRA_Close failed(d_b): ",d_b);
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    if (me == 0) {
        printf("Time including DRA_Close\n");
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }

    if (me == 0) printf("\n");
    if (me == 0) printf("disk array closed\n");
    if (me == 0) fflush(stdout);
#endif /* SWITCH */
    /*.......................................................................*/

#if SWITCH
    if (me == 0) printf("\n");
    if (me == 0) printf("opening disk array\n");
    if (DRA_Open(filename1, DRA_R, &d_b) != 0) GA_Error("DRA_Open failed",0);
    if (me == 0) printf("non alligned blocking read\n");
    if (me == 0) fflush(stdout);
    tt0 = MP_TIMER();
    if (NDRA_Read_section(FALSE, g_b, glo, ghi, d_b, dlo, dhi, &req) != 0)
        GA_Error("NDRA_Read_section failed:",0);
    if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed: ",req);
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    if (me == 0) {
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }
    plus = 1.0;
    minus = -1.0;
    GA_Add(&plus, g_a, &minus, g_b, g_b);
    err = GA_Ddot(g_b, g_b);
    if (err != 0) {
        if (me == 0) printf("BTW, we have error = %f\n",err);
    } else {
        if (me == 0) printf("OK\n");
    }
    if (DRA_Delete(d_b) != 0) GA_Error("DRA_Delete failed",0);
    /*.......................................................................*/
    if (me == 0) printf("\n");
    if (me == 0) printf("opening disk array\n");
    if (DRA_Open(filename, DRA_R, &d_a) != 0) GA_Error("DRA_Open failed",0);
    if (me == 0) printf("alligned blocking read\n");
    if (me == 0) fflush(stdout);
    tt0 = MP_TIMER();
    if (NDRA_Read(g_b, d_a, &req) != 0) GA_Error("NDRA_Read failed:",0);
    if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed: ",req);
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    if (me == 0) {
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }
    GA_Add(&plus, g_a, &minus, g_b, g_b);
    err = GA_Ddot(g_b, g_b);
    if (err != 0) {
        if (me == 0) printf("BTW, we have error = %f\n",err);
    } else {
        if (me == 0) printf("OK\n");
    }
    if (DRA_Delete(d_a) != 0) GA_Error("DRA_Delete failed",0);
#else /* SWITCH */
    if (me == 0) printf("\n");
    if (me == 0) printf("opening disk array\n");
    if (DRA_Open(filename, DRA_R, &d_a) != 0) GA_Error("DRA_Open failed",0);
    if (me == 0) printf("alligned blocking read\n");
    if (me == 0) fflush(stdout);
    tt0 = MP_TIMER();
    if (NDRA_Read(g_b, d_a, &req) != 0) GA_Error("NDRA_Read failed:",0);
    if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed: ",req);
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    if (me == 0) {
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }
    plus = 1.0;
    minus = -1.0;
    GA_Add(&plus, g_a, &minus, g_b, g_b);
    err = GA_Ddot(g_b, g_b);
    if (err != 0) {
        if (me == 0) printf("BTW, we have error = %f\n",err);
    } else {
        if (me == 0) printf("OK\n");
    }
    if (DRA_Delete(d_a) != 0) GA_Error("DRA_Delete failed",0);
    /*.......................................................................*/

    if (me == 0) printf("\n");
    if (me == 0) printf("opening disk array\n");
    if (DRA_Open(filename1, DRA_R, &d_b) != 0) GA_Error("DRA_Open failed",0);
    if (me == 0) printf("non alligned blocking read\n");
    if (me == 0) fflush(stdout);
    tt0 = MP_TIMER();
    if (NDRA_Read_section(FALSE, g_b, glo, ghi, d_b, dlo, dhi, &req) != 0)
        GA_Error("NDRA_Read_section failed:",0);
    if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed: ",req);
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    if (me == 0) {
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }
    GA_Add(&plus, g_a, &minus, g_b, g_b);
    err = GA_Ddot(g_b, g_b);
    if (err != 0) {
        if (me == 0) printf("BTW, we have error = %f\n",err);
    } else {
        if (me == 0) printf("OK\n");
    }
    if (DRA_Delete(d_b) != 0) GA_Error("DRA_Delete failed",0);
#endif /* SWITCH */
    /*.......................................................................*/
    GA_Destroy(g_a);
    GA_Destroy(g_b);
    /*.......................................................................*/
#if TEST_TRANSPOSE
    /* Test transpose function for DRAs */
    dims[0] = n;
    for (i=1; i<ndim; i++) dims[i] = n/2;
    for (i=0; i<ndim; i++) chunk[i] = 1;
    if (me == 0) printf("Creating asymmetric arrays to test transpose\n\n");
    g_c = NGA_Create(MT_DBL, ndim, dims, "c", chunk);
    if (!g_c) GA_Error("NGA_Create failed: c", 0);
    g_d = NGA_Create(MT_DBL, ndim, dims, "d", chunk);
    if (!g_d) GA_Error("NGA_Create failed: d", 0);
    if (me == 0) printf("done\n");
    if (me == 0) fflush(stdout);

    /*     initialize g_a, g_b with random values
           ... use ga_access to avoid allocating local buffers for ga_put */

    GA_Sync();
    NGA_Distribution(g_c, me, glo, ghi);
    NGA_Access(g_c, glo, ghi, &index, ld);
    isize = 1;
    for (i=0; i<ndim; i++) isize *= (ghi[i]-glo[i]+1);
    fill_random(index, isize);
    GA_Sync();
    GA_Zero(g_c);
    GA_Copy(g_c,g_d);

    for (i=0; i<ndim; i++) {
        ddims[i] = m;
        reqdims[i] = n;
    }
    filename_check(filename2, FNAME2, FNAME2_ALT);
    if (me == 0) printf("Creating DRA for transpose test\n");
    if (NDRA_Create(MT_DBL, ndim, ddims, "C", filename2, DRA_RW,
                reqdims, &d_c) != 0) {
        GA_Error("NDRA_Create failed(d_c): ",0);
    }
    if (me == 0) printf("done\n");
    if (me == 0) fflush(stdout);
    GA_Sync();
    for (i=0; i<ndim-1; i++) {
        dlo[i] = 1;
        dhi[i] = n/2;
    }
    dlo[ndim-1] = 1;
    dhi[ndim-1] = n;
    glo[0] = 0;
    ghi[0] = n-1;
    for (i=1; i<ndim; i++) {
        glo[i] = 0;
        ghi[i] = n/2-1;
    }
    if (me == 0) printf("non-aligned blocking write with transpose\n");
    tt0 = MP_TIMER();
    if (NDRA_Write_section(TRUE,g_c,glo,ghi,d_c,dlo,dhi,&req) != 0)
        GA_Error("NDRA_Write_section (transpose) failed: ",0);
    if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed: ",req);
    isize = 1;
    for (i=0; i<ndim; i++) isize *= (ghi[i]-glo[i]+1);
    mbytes = 1.e-6 * (double)(isize*sizeof(double));
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    if (me == 0) {
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }
    if (DRA_Close(d_c) != 0) GA_Error("DRA_Close failed(d_c): ",d_c);
    if (me == 0) printf("\n");
    if (me == 0) printf("disk array closed\n");
    if (me == 0) fflush(stdout);

    if (me == 0) printf("\n");
    if (me == 0) printf("opening disk array\n");
    if (DRA_Open(filename2, DRA_R, &d_c) != 0) GA_Error("DRA_Open failed",0);

    GA_Zero(g_c);
    if (me == 0) printf("non-aligned blocking read with transpose\n");
    tt0 = MP_TIMER();
    if (NDRA_Read_section(TRUE,g_c,glo,ghi,d_c,dlo,dhi,&req) != 0)
        GA_Error("NDRA_Read_section (transpose) failed: ",0);
    if (DRA_Wait(req) != 0) GA_Error("DRA_Wait failed: ",req);
    tt1 = MP_TIMER() - tt0;
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    if (me == 0) {
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }
    GA_Add(&plus, g_c, &minus, g_d, g_d);
    err = GA_Ddot(g_d, g_d);
    if (err != 0) {
        if (me == 0) printf("BTW, we have error = %f\n",err);
    } else {
        if (me == 0) printf("OK\n");
    }
    if (DRA_Delete(d_c) != 0) GA_Error("DRA_Delete failed",0);
    /*.......................................................................*/
    GA_Destroy(g_c);
    GA_Destroy(g_d);
#endif /* TEST_TRANSPOSE */
}


int main(int argc, char **argv)
{
    int status, me;
    int max_arrays = 10;
    double max_sz = 1e8, max_disk = 1e10, max_mem = 1e6;
    int numfiles, numioprocs;
#if defined(IBM)
    int stack = 9000000, heap = 4000000;
#else
    int stack = 12000000, heap = 8000000;
#endif

    MP_INIT(argc,argv);
    GA_Initialize();
    if (!GA_Uses_ma()) {
        if (GA_Nodeid() == 0) printf("GA not using MA\n");
        stack = 100000;
        heap = 2*sizeof(double)*pow(SIZE,NDIM)/GA_Nnodes()+10000000;
    }

    me = GA_Nodeid();
    if (MA_init(MT_F_DBL, stack, heap) ) {
        if (DRA_Init(max_arrays, max_sz, max_disk, max_mem) != 0)
            GA_Error("DRA_Init failed: ",0);
        if (USER_CONFIG == 0) {
            numfiles = -1;
            numioprocs = -1;
        } else if (USER_CONFIG == 1) {
            numfiles = 1;
            numioprocs = GA_Cluster_nnodes();
        } else if (USER_CONFIG == 2) {
            numfiles = GA_Cluster_nnodes();
            numioprocs = GA_Cluster_nnodes();
        } else {
            numfiles = 1;
            numioprocs = 1;
        }
        if (me==0) {
            printf("Disk resident arrays configured as:\n");
            printf("    Number of files: %d\n",numfiles);
            printf("    Number of I/O processors: %d\n",numioprocs);
        }
        DRA_Set_default_config(numfiles,numioprocs);
        if (me == 0) printf("\n");
        if (me == 0) printf("TESTING PERFORMANCE OF DISK ARRAYS\n");
        if (me == 0) printf("\n");
        test_io_dbl();
        status = DRA_Terminate();
        GA_Terminate();
    } else {
        printf("MA_init failed\n");
    }
    if(me == 0) printf("all done ...\n");
    MP_FINALIZE();
    return 0;
}
