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

#define FNAME     "/scratch/da.try"
#define FNAME_ALT "/tmp/da.try"

#include "dra.h"
#include "eaf.h"
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
 * required is approximately NFACTOR**NDIM times this amount.
 */
#define USER_CONFIG 0

#define NDIM 3
#define SIZE 250
#define NFACTOR 7

/*
#define NDIM 2
#define SIZE 4000
#define NFACTOR 3

#define NDIM 1
#define SIZE 12500000
#define NFACTOR 640
*/

/*#define MAXDIM 7*/
/*#define TRUE (logical)1*/
/*#define FALSE (logical)0*/

#define MULTFILES 0

#ifdef SOLARIS
#   if MULTFILES
#       define USEMULTFILES 1
#   endif
#else
#   define USEMULTFILES 0
#endif

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876


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
    int n, ndim = NDIM;
    double err, tt0, tt1, mbytes;
    int g_a, g_b, d_a;
    int i, itmp, j, req, loop;
    int glo[MAXDIM],ghi[MAXDIM];
    dra_size_t dlo[MAXDIM],dhi[MAXDIM];
    dra_size_t ddims[MAXDIM],reqdims[MAXDIM];
    dra_size_t m;
    int index[MAXDIM], dims[MAXDIM];
    int me, nproc, isize;
    double *ptr;
    double plus, minus;
    int ld[MAXDIM], chunk[MAXDIM];
    char filename[80];
    FILE *fd;

    n = SIZE;
    m = ((dra_size_t)NFACTOR)*((dra_size_t)SIZE);

    loop  = 1;
    for (i=0; i<ndim; i++) loop *= NFACTOR;
    req = -1;
    nproc = GA_Nnodes();
    me    = GA_Nodeid();

    if (me == 0) {
        printf("Creating temporary global arrays %d",n);
        for (i=1; i<ndim; i++) {
            printf(" x %d",n);
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
    NGA_Access(g_a, glo, ghi, &ptr, ld);
    isize = 1;
    for (i=0; i<ndim; i++) isize *= (ghi[i]-glo[i]+1);
    fill_random(ptr, isize);
    GA_Sync();
    GA_Zero(g_b);


    /*.......................................................................*/
    if (me == 0) {
        printf("Creating Disk array %ld",m);
        for (i=1; i<ndim; i++) {
            printf(" x %ld",m);
        }
        printf("\n");
    }
    if (me == 0) fflush(stdout);
    for (i=0; i<ndim; i++) {
        ddims[i] = m;
        reqdims[i] = (dra_size_t)n;
    }
    GA_Sync();
    strcpy(filename,FNAME);
    if (! (fd = fopen(filename, "w"))) {
        strcpy(filename,FNAME_ALT);
        if (! (fd = fopen(filename, "w"))) {
            GA_Error("open failed",0);
        }
    }
    fclose(fd);
    if (NDRA_Create(MT_DBL, ndim, ddims, "A", filename, DRA_RW,
                reqdims, &d_a) != 0) {
        GA_Error("NDRA_Create failed(d_a): ",0);
    }
    if (me == 0) printf("testing write\n");
    fflush(stdout);
    tt1 = 0.0;
    for (i=0; i<loop; i++) {
        itmp=i;
        for (j=0; j<ndim; j++) {
            index[j] = itmp%NFACTOR;
            itmp = (itmp - index[j])/NFACTOR;
        }
        for (j=0; j<ndim; j++) {
            glo[j] = 0;
            ghi[j] = SIZE - 1;
            dlo[j] = ((dra_size_t)index[j])*((dra_size_t)SIZE);
            dhi[j] = (((dra_size_t)index[j])+(dra_size_t)1)
                * ((dra_size_t)SIZE) - (dra_size_t)1;
        }
        tt0 = MP_TIMER();
        if (NDRA_Write_section(FALSE, g_a, glo, ghi,
                    d_a, dlo, dhi, &req) != 0) {
            GA_Error("ndra_write_section failed:",0);
        }
        if (DRA_Wait(req) != 0) {
            GA_Error("DRA_Wait failed(d_a): ",req);
        }
        tt1 += (MP_TIMER() - tt0);
    }
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    mbytes = 1.e-6 * (double)(pow(m,ndim)*sizeof(double));
    if (me == 0) {
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }

    if (DRA_Close(d_a) != 0) {
        GA_Error("DRA_Close failed(d_a): ",d_a);
    }

    if (me == 0) printf("\n");
    if (me == 0) printf("disk array closed\n");
    if (me == 0) fflush(stdout);

    /*..........................................................*/

    if (me == 0) printf("\n");
    if (me == 0) printf("opening disk array\n");
    if (DRA_Open(filename, DRA_R, &d_a) != 0) {
        GA_Error("DRA_Open failed",0);
    }
    if (me == 0) printf("testing read\n");
    /*  printf("testing read on proc %d\n",me); */
    if (me == 0) fflush(stdout);
    tt1 = 0.0;
    for (i=0; i<loop; i++) {
        itmp=i;
        for (j=0; j<ndim; j++) {
            index[j] = itmp%NFACTOR;
            itmp = (itmp - index[j])/NFACTOR;
        }
        for (j=0; j<ndim; j++) {
            glo[j] = 0;
            ghi[j] = SIZE - 1;
            dlo[j] = ((dra_size_t)index[j])*((dra_size_t)SIZE);
            dhi[j] = (((dra_size_t)index[j])+(dra_size_t)1)
                * ((dra_size_t)SIZE) - (dra_size_t)1;
        }
        tt0 = MP_TIMER();
        if (NDRA_Read_section(FALSE, g_b, glo, ghi,
                    d_a, dlo, dhi, &req) != 0) {
            GA_Error("ndra_read_section failed:",0);
        }
        if (DRA_Wait(req) != 0) {
            GA_Error("DRA_Wait failed(d_a): ",req);
        }
        tt1 += (MP_TIMER() - tt0);
        plus = 1.0;
        minus = -1.0;
        GA_Add(&plus, g_a, &minus, g_b, g_b);
        err = GA_Ddot(g_b, g_b);
        if (err != 0) {
            if (me == 0) {
                printf("BTW, we have error = %f on loop value %d\n",
                        err,i);
            }
            GA_Error(" bye",0);
        }
    }
    GA_Dgop(&tt1,1,"+");
    tt1 = tt1/((double)nproc);
    if (me == 0) {
        printf("%11.2f MB  time = %11.2f rate = %11.3f MB/s\n",
                mbytes,tt1,mbytes/tt1);
    }
    if (DRA_Delete(d_a) != 0) GA_Error("DRA_Delete failed",0);
    /*.......................................................................*/
    GA_Destroy(g_a);
    GA_Destroy(g_b);
}


int main(int argc, char **argv)
{
    int status, me;
    int max_arrays = 10;
    double max_sz = 1e8, max_disk = 1e10, max_mem = 1e6;
    int stack = 120000, heap = 3200000;
    int numfiles, numioprocs;
    int total, nproc;

    MP_INIT(argc,argv);
    GA_Initialize();
    me = GA_Nodeid();
    nproc = GA_Nnodes();
    total = pow(SIZE,NDIM)*sizeof(double);
    if (!GA_Uses_ma()) {
        if (GA_Nodeid() == 0) {
            printf("GA is not using MA\n");
        }
        stack = 100000;
        heap  = (int)(2.2*(float)(total));
    }

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
