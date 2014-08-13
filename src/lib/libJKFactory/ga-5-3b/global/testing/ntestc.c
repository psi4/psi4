#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#include "macdecls.h"
#include "ga.h"
#include "mp3.h"

/* utilities for GA test programs */
#include "testutil.h"

#define N 10            /* first dimension  */
#define NDIM 4          /* number of dimensions */
#define BASE 0
#define PERMUTE_ 

/*#define NEW_API*/


/*\ fill n-dimensional array section with value
\*/
void fill_patch(double *ptr, int dim[], int ld[], int ndim, double val)
{
     int i, j, stride=1;

     switch (ndim){
     case 0: GA_Error("fill_patch: error",ndim);
     case 1: for(i=0;i <dim[0];i++)ptr[i]=val;
             break;
     case 2: for(i=0; i< dim[0]; i++){
               double *arr = ptr + i*ld[0];
               for(j=0; j< dim[1]; j++)arr[j]=val; 
             }
             break;
     default: 
             for(i=0; i<ndim-1; i++)stride *=ld[i];
             for(i=0; i<dim[0]; i++){
                 double *arr = ptr + stride*i;
                 fill_patch(arr, dim+1, ld+1, ndim-1, val);
             } 
     }
}
        


void do_work()
{
int g_a, g_b;
int me=GA_Nodeid(), nproc=GA_Nnodes(), proc, loop;
int dims[NDIM], lo[NDIM], hi[NDIM], block[NDIM], ld[NDIM-1];
int64_t lo64[NDIM], hi64[NDIM];
int i,d,*proclist, offset;
int adims[NDIM], ndim,type;
typedef struct {
       int lo[NDIM];
       int hi[NDIM];
} patch_t;
patch_t *regions;
int *map;
double *buf;

     /***** create array A with default distribution  *****/
     if(me==0){printf("Creating array A\n"); fflush(stdout);}
     for(i = 0; i<NDIM; i++)dims[i] = N*(i+1);
#ifdef NEW_API
     g_a = GA_Create_handle();
     GA_Set_data(g_a,NDIM,dims,MT_F_DBL);
     GA_Set_array_name(g_a,"array A");
     (void)GA_Allocate(g_a);
#else
     g_a = NGA_Create(MT_F_DBL, NDIM, dims, "array A", NULL);
#endif
     if(!g_a) GA_Error("create failed: A",1); 
     if(me==0)printf("OK\n\n");

     /* print info about array we got */
     NGA_Inquire(g_a, &type, &ndim, adims);
     GA_Print_distribution(g_a);

     GA_Sync();
     /* duplicate array A with ga_create irreg rather than ga_duplicate
      * -- want to show distribution control 
      * -- with ga_duplicate it would be g_b=GA_Duplicate(g_a,name) 
      */
     if(me==0)printf("\nReconstructing distribution description for A\n");

     /* get memory for arrays describing distribution */
     proclist = (int*)malloc(nproc*sizeof(int));
     if(!proclist)GA_Error("malloc failed for proclist",1);
     regions = (patch_t*)malloc(nproc*sizeof(patch_t));
     if(!regions)GA_Error("malloc failed for regions",1);
     map = (int*)malloc((nproc+ndim)*sizeof(int)); /* ubound= nproc+mdim */
     if(!map)GA_Error("malloc failed for map",1);

     /* first find out how array g_a is distributed */
     for(i=0;i<ndim;i++)lo[i]=BASE;
     for(i=0;i<ndim;i++)hi[i]=adims[i] -1 + BASE;
     for(i=0;i<ndim;i++)lo64[i]=BASE;
     for(i=0;i<ndim;i++)hi64[i]=adims[i] -1 + BASE;
     proc = NGA_Locate_nnodes(g_a, lo, hi);
     if(proc<1) GA_Error("error in NGA_Locate_nnodes",proc);
     proc = NGA_Locate_nnodes64(g_a, lo64, hi64);
     if(proc<1) GA_Error("error in NGA_Locate_nnodes",proc);
     proc = NGA_Locate_region(g_a, lo, hi, (int*)regions, proclist);
     if(proc<1) GA_Error("error in NGA_Locate_region",proc);

     /* determine blocking for each dimension */
     for(i=0;i<ndim;i++)block[i]=0;
     for(i=0;i<ndim;i++)adims[i]=0;
     
     offset =0;
     for(d=0; d<ndim; d++)
         for(i=0;i<proc;i++)
             if( regions[i].hi[d]>adims[d] ){
                map[offset] = regions[i].lo[d];
                offset++;
                block[d]++;
                adims[d]= regions[i].hi[d];
             }
            
     if(me==0){
        printf("Distribution map contains %d elements\n",offset); 
        print_subscript("number of blocks for each dimension",ndim,block,"\n");
        print_subscript("distribution map",offset,map,"\n\n");
        fflush(stdout);
     }
     
     if(me==0)printf("Creating array B applying distribution of A\n");

#    ifdef USE_DUPLICATE
       g_b = GA_Duplicate(g_a,"array B");
#    else
       g_b = NGA_Create_irreg(MT_F_DBL, NDIM, dims, "array B", block,map);
#    endif
     if(!g_b) GA_Error("create failed: B",1); 
     if(me==0)printf("OK\n\n");
     free(proclist); free(regions); free(map);
     
     GA_Print_distribution(g_b);

     GA_Sync();

     if(me==0){
       printf("\nCompare distributions of A and B\n");
       if(GA_Compare_distr(g_a,g_b))
          printf("Failure: distributions NOT identical\n");
       else 
          printf("Success: distributions identical\n");
       fflush(stdout);
     }
       

     if(me==0){
        printf("\nAccessing local elements of A: set them to the owner process id\n");
        fflush(stdout);
     }
     GA_Sync();

     NGA_Distribution(g_a,me,lo,hi);

     if(hi[0]>=0){/* -1 means no elements stored on this processor */
         double *ptr;
         int locdim[NDIM];
         NGA_Access(g_a, lo,hi, &ptr, ld);
         for(i=0;i<ndim;i++)locdim[i]=hi[i]-lo[i]+1;
         fill_patch(ptr, locdim, ld, ndim,(double)me);
     }

     for(i=0;i<nproc; i++){
       if(me==i && hi[0]>=0){
         char msg[100];
         sprintf(msg,"%d: leading dimensions",me);
         print_subscript(msg,ndim-1,ld,"\n");
         fflush(stdout);
       }
       GA_Sync();
     }
     
     GA_Sync();
     if(me==0)printf("\nRandomly checking the update using ga_get on array sections\n");
     GA_Sync();

     /* show ga_get working and verify array updates 
      * every process does N random gets
      * for simplicity get only a single row at a time
      */
     srand(me); /* different seed for every process */
     hi[ndim-1]=adims[ndim-1] -1 + BASE;
     for(i=1;i<ndim-1; i++)ld[i]=1; ld[ndim-2]=adims[ndim-1] -1 + BASE;

     /* get buffer memory */
     buf = (double*)malloc(adims[ndim-1]*sizeof(double));
     if(!buf)GA_Error("malloc failed for buf",1);

     /* half of the processes check the result */
     if(me<=nproc/2) 
     for(loop = 0; loop< N; loop++){ /* task parallel loop */
         lo[ndim-1]=BASE;
         for (i= 0; i < ndim -1; i ++){
              lo[i] = hi[i] = rand()%adims[i]+BASE; 
         }

         /* print_subscript("getting",ndim,lo,"\n");*/
         NGA_Get(g_a,lo,hi,buf,ld); 
         
         /* check values */
         for(i=0;i<adims[ndim-1]; i++){
             int p = NGA_Locate(g_a, lo);
             if((double)p != buf[i]) {
                char msg[100];
                sprintf(msg,"%d: wrong value: %d != %f a",me, p, buf[i]);
                print_subscript(msg,ndim,lo,"\n");
                GA_Error("Error - bye",i);  
             }
             lo[ndim-1]++;
          }
     }
             
     free(buf);
     GA_Sync();
           
     if(me==0)printf("OK\n");
     
     GA_Destroy(g_a);
     GA_Destroy(g_b);
}
     

int main(argc, argv)
int argc;
char **argv;
{
Integer heap=300000, stack=300000;
int me, nproc;

    MP_INIT(argc,argv);

    GA_INIT(argc,argv);                           /* initialize GA */

    nproc = GA_Nnodes();
    me = GA_Nodeid();

    if(me==0) printf("Using %d processes\n\n",nproc);

    if(!MA_init((Integer)MT_F_DBL, stack/nproc, heap/nproc))
       GA_Error("MA_init failed bytes= %d",stack+heap);   

#ifdef PERMUTE
      {
        int i, *list = (int*)malloc(nproc*sizeof(int));
        if(!list)GA_Error("malloc failed",nproc);

        for(i=0; i<nproc;i++)list[i]=nproc-1-i;

        GA_Register_proclist(list, nproc);
        free(list);
      }
#endif

    if(GA_Uses_fapi())GA_Error("Program runs with C API only",1);
    
    do_work();

    if(me==0)printf("\nAll tests successful\n\n");
    GA_Terminate();

    MP_FINALIZE();

    return 0;
}

