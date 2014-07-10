#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#include "dra.h"
#include "ga.h"
#include "macommon.h"
#include "macdecls.h"
#include "mp3.h"

#define ERROR(msg,code){printf("ERROR:%s\n",(msg)); fflush(stdout); exit(1);}

#define BUFSIZE 8000000
int main(int argc, char **argv)
{
    int heap=400000, stack=400000;
    int me, nproc;
    int max_arrays=2;
    double max_sz=1e8, max_disk=2e8, max_mem=1e6;
    int d_a, mode=DRA_R;
    int g_a,rows,cols;
    char name[1024], fname[1024];
    logical transp=0;
    int reqid;
    int i,ilo,ihi,jlo,jhi,type,ndim,glo[2],ghi[2],gld[1],gdims[2];
    dra_size_t dlo[2],dhi[2],ddims[2];
    size_t size, nitems;
    void *ptr;

    if(argc!=6){
        printf("Usage: dra2arviz <dra_filename> <ilo> <ihi> <jlo> <jhi>\n");
        printf("       dra_filename is the meta-file name for disk resident array\n");
        printf("       [ilo:ihi, jlo:jhi]  array section to read\n\n\n");
        return(1);
    }

    MP_INIT(argc,argv);
    MP_MYID(&me);
    MP_PROCS(&nproc);

    heap /= nproc;
    stack /= nproc;
    if(! MA_init((Integer)MT_F_DBL, stack, heap))
        GA_Error("MA_init failed",stack+heap); /* initialize memory allocator*/
    GA_Initialize();                         /* initialize GA */

    if(nproc != 1)ERROR("Error: does not run in parallel",nproc);

    if(DRA_Init(max_arrays, max_sz, max_disk, max_mem))
        ERROR("initialization failed",0);

    if(DRA_Open(argv[1], mode, &d_a)) ERROR("dra_open failed",0);

    ilo = atoi(argv[2]); ihi = atoi(argv[3]);
    jlo = atoi(argv[4]); jhi = atoi(argv[5]);
    dlo[0] = ilo;
    dlo[1] = jlo;
    dhi[0] = ihi;
    dhi[1] = jhi;
    rows = ihi - ilo +1; 
    cols = jhi - jlo +1; 
    glo[0] = 0;
    glo[1] = 0;
    ghi[0] = rows-1;
    ghi[1] = cols-1;
    gdims[0] = rows;
    gdims[1] = cols;

    if(NDRA_Inquire(d_a, &type, &ndim, ddims, name, fname))
        ERROR("dra_inquire failed",0);

    switch (type) {
        case  MT_F_INT:  size = sizeof(Integer); break;
        case  MT_F_DBL:  size = sizeof(DoublePrecision); break;
        case  MT_F_DCPL:  size = sizeof(DoubleComplex); break;
        default: ERROR("type not supported",type);
    }

    g_a = NGA_Create(type, 2, gdims, "temp", NULL);

    if(NDRA_Read_section(transp, g_a, glo, ghi,
                d_a, dlo, dhi, &reqid))
        ERROR("dra_read_section failed",0);

    if(DRA_Wait(reqid)) ERROR("dra_wait failed",0);

    NGA_Access(g_a, glo, ghi, &ptr, gld);

    if(gld[0] != rows) ERROR("ld != rows",gld[0]); 

    fwrite("OK\n",1,3,stdout);
    nitems = (size_t)rows;
    /* write data by columns */
    for(i=0; i<cols; i++){
        if(type == MT_F_DBL)
            fwrite(ptr,size,nitems,stdout);
        else if(type == MT_F_DCPL)
            fwrite(ptr,size,nitems,stdout);
        else
            fwrite(ptr,size,nitems,stdout);
    }
    fflush(stdout);

    GA_Destroy(g_a);

    GA_Terminate();
    MP_FINALIZE();
    return 0;
}
