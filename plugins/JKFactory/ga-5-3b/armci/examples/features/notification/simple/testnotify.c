#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: testnotify.c,v 1.1.2.1 2007-06-20 17:42:09 vinod Exp $ */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRINGS_H
#   include <strings.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif

#define DEBUG__ 

#include "armci.h"
#include "message.h"

#define DIM1 5
#define DIM2 3
#ifdef __sun
/* Solaris has shared memory shortages in the default system configuration */
# define DIM3 6
# define DIM4 5
# define DIM5 4
#elif defined(__alpha__)
# define DIM3 8
# define DIM4 5
# define DIM5 6
#else
# define DIM3 8
# define DIM4 9
# define DIM5 7
#endif
#define DIM6 3
#define DIM7 2


#define OFF 1
#define EDIM1 (DIM1+OFF)
#define EDIM2 (DIM2+OFF)
#define EDIM3 (DIM3+OFF)
#define EDIM4 (DIM4+OFF)
#define EDIM5 (DIM5+OFF)
#define EDIM6 (DIM6+OFF)
#define EDIM7 (DIM7+OFF)

#define DIMS 4
#define MAXDIMS 7
#define MAX_DIM_VAL 50 
#define LOOP 200

#define BASE 100.
#define MAXPROC 128
#define TIMES 100

#ifdef CRAY
# define ELEMS 800
#else
# define ELEMS 200
#endif



/***************************** macros ************************/
#define COPY(src, dst, bytes) memcpy((dst),(src),(bytes))
#define ARMCI_MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define ARMCI_MIN(a,b) (((a) <= (b)) ? (a) : (b))
#define ARMCI_ABS(a) (((a) <0) ? -(a) : (a))

/***************************** global data *******************/
int me, nproc;
void* work[MAXPROC]; /* work array for propagating addresses */



/*\ generate random range for a section of multidimensional array 
\*/
void get_range(int ndim, int dims[], int lo[], int hi[]) 
{
    int dim;
    for(dim=0; dim <ndim;dim++){
        int toss1, toss2;
        toss1 = rand()%dims[dim];
        toss2 = rand()%dims[dim];
        if(toss1<toss2){
            lo[dim]=toss1;
            hi[dim]=toss2;
        }else {
              hi[dim]=toss1;
            lo[dim]=toss2;
        }
    }
}



/*\ generates a new random range similar to the input range for an array with specified dimensions
\*/
void new_range(int ndim, int dims[], int lo[], int hi[],int new_lo[], int new_hi[])
{
    int dim;
    for(dim=0; dim <ndim;dim++){
        int toss, range;
        int diff = hi[dim] -lo[dim]+1;
        assert(diff <= dims[dim]);
                range = dims[dim]-diff;
                toss = (range > 0)? rand()%range : lo[dim];
        new_lo[dim] = toss;
        new_hi[dim] = toss + diff -1;
        assert(new_hi[dim] < dims[dim]);
        assert(diff == (new_hi[dim] -new_lo[dim]+1));
    }
}


/*\ print range of ndim dimensional array with two strings before and after
\*/
void print_range(char *pre,int ndim, int lo[], int hi[], char* post)
{
    int i;

    printf("%s[",pre);
    for(i=0;i<ndim;i++){
        printf("%d:%d",lo[i],hi[i]);
        if(i==ndim-1)printf("] %s",post);
        else printf(",");
    }
}

/*\ print subscript of ndim dimensional array with two strings before and after
\*/
void print_subscript(char *pre,int ndim, int subscript[], char* post)
{
    int i;

    printf("%s [",pre);
    for(i=0;i<ndim;i++){
        printf("%d",subscript[i]);
        if(i==ndim-1)printf("] %s",post);
        else printf(",");
    }
}


/*\ initialize array: a[i,j,k,..]=i+100*j+10000*k+ ... 
\*/
void init(double *a, int ndim, int elems, int dims[])
{
  int idx[MAXDIMS];
  int i,dim;

     for(i=0; i<elems; i++){
        int Index = i;
        double field, val;
        
        for(dim = 0; dim < ndim; dim++){
            idx[dim] = Index%dims[dim];
            Index /= dims[dim];
        }
        
                field=1.; val=0.;
        for(dim=0; dim< ndim;dim++){
            val += field*idx[dim];
            field *= BASE;
        }
        a[i] = val;
        /* printf("(%d,%d,%d)=%6.0f",idx[0],idx[1],idx[2],val); */
    }
}


/*\ compute Index from subscript
 *  assume that first subscript component changes first
\*/
int Index(int ndim, int subscript[], int dims[])
{
    int idx = 0, i, factor=1;
    for(i=0;i<ndim;i++){
        idx += subscript[i]*factor;
        factor *= dims[i];
    }
    return idx;
}


void update_subscript(int ndim, int subscript[], int lo[], int hi[], int dims[])
{
    int i;
    for(i=0;i<ndim;i++){
        if(subscript[i] < hi[i]) { subscript[i]++; return; }
        subscript[i] = lo[i];
    }
}

    

void compare_patches(double eps, int ndim, double *patch1, int lo1[], int hi1[],
                     int dims1[],double *patch2, int lo2[], int hi2[], 
                     int dims2[])
                               
{
int i,j, elems=1;    
int subscr1[MAXDIMS], subscr2[MAXDIMS];
double diff,max;

    for(i=0;i<ndim;i++){
    int diff = hi1[i]-lo1[i];
      assert(diff == (hi2[i]-lo2[i]));
      assert(diff < dims1[i]);
      assert(diff < dims2[i]);
      elems *= diff+1;
      subscr1[i]= lo1[i];
      subscr2[i]=lo2[i];
    }
    for(j=0; j< elems; j++){ 
    int idx1=0, idx2=0, offset1=0, offset2=0;
      idx1 = Index(ndim, subscr1, dims1);
      idx2 = Index(ndim, subscr2, dims2);
      if(j==0){
        offset1 =idx1;
    offset2 =idx2;
      }
      idx1 -= offset1;
      idx2 -= offset2;
      diff = patch1[idx1] - patch2[idx2];
      max  = ARMCI_MAX(ARMCI_ABS(patch1[idx1]),ARMCI_ABS(patch2[idx2]));
      if(max == 0. || max <eps) max = 1.; 
      if(eps < ARMCI_ABS(diff)/max){
       char msg[48];
         sprintf(msg,"(proc=%d):%f",me,patch1[idx1]);
     print_subscript("ERROR: a",ndim,subscr1,msg);
     sprintf(msg,"%f\n",patch2[idx2]);
     print_subscript(" b",ndim,subscr2,msg);
         fflush(stdout);
         sleep(1);
         ARMCI_Error("Bailing out",0);
      }
      update_subscript(ndim, subscr1, lo1,hi1, dims1);
      update_subscript(ndim, subscr2, lo2,hi2, dims2);
    }
}


void create_array(void *a[], int elem_size, int ndim, int dims[])
{
armci_size_t bytes=elem_size;
int i, rc;

    assert(ndim<=MAXDIMS);
    for(i=0;i<ndim;i++)bytes*=dims[i];

    rc = ARMCI_Malloc(a, bytes);
    assert(rc==0);
     
#ifdef DEBUG_
    printf("%d after malloc ndim=%d b=%d ptr=%p\n",me,ndim,(int) bytes,a[me]);
    fflush(stdout);
#endif

    assert(a[me]);
    bzero(a[me],bytes);
}

void destroy_array(void *ptr[])
{
    armci_msg_barrier();
    assert(!ARMCI_Free(ptr[me]));
}


int get_elems(int ndim, int *stride, int *dims, int size)
{
int elems=1,i;
    stride[0]=size;
    for(i=0;i<ndim;i++){
      stride[i] *= dims[i];
      if(i<ndim-1){
        stride[i+1] = stride[i];
      }
      elems *= dims[i];
    }
    return elems;
}


int dimsB[MAXDIMS]={30,EDIM2,EDIM3,EDIM4,EDIM5,EDIM6,EDIM7};
#define BATCH 5
    
void test_notify(int ndim)
{
int lo[MAXDIMS], hi[MAXDIMS], count[MAXDIMS];
int stride[MAXDIMS];
int dim,elems;
int i,Idx=1,idx=0;
void *b[MAXPROC], *a[MAXPROC];
int left = (me+nproc-1) % nproc;
int right = (me+1) % nproc;
int loopcnt=1, less=2, strl; /* less>1 takes a partial plane */


    /* create shared and local arrays */
    create_array(b, sizeof(double),ndim,dimsB);
    create_array(a, sizeof(double),ndim,dimsB);

    elems = get_elems(ndim, stride, dimsB, sizeof(double));
    init((double*)a[me], ndim, elems, dimsB);

    for(i=0; i<ndim; i++){
       lo[i]=0; 
       hi[i]= (less > dimsB[i]) ? dimsB[i]-1: dimsB[i]-less;
       count[i]=hi[i]-lo[i]+1;
    }
    count[0]*=sizeof(double);

    for(i=0; i<ndim-1; i++)Idx *= dimsB[i];

    ARMCI_Barrier();
    if(me==0){
       printf("--------array[%d",dimsB[0]);
       for(dim=1;dim<ndim;dim++)printf(",%d",dimsB[dim]);
       printf("]--------\n");
       fflush(stdout);
    }

    ARMCI_Barrier();
    loopcnt = (ndim>1)? dimsB[ndim-1] : 1;
    strl    = (ndim>1)? ndim-2: 0; /* strides of the subpatch to transfer */

    for(i=0;i<loopcnt;i++){
        int wc;

        if(me==0){

          ARMCI_PutS((double*)a[me]+idx, stride, 
                     (double*)b[left]+idx, stride, count, strl, left);
#if DEBUG_ 
          printf("%d-%d: ps=%p pd=%p i=%d idx=%d count=%d\n",me,left,(double*)
              a[me]+idx, (double*)b[left]+idx,i, idx,count[0]); fflush(stdout); 
#endif
          (void)armci_notify(left);
          (void)armci_notify_wait(right,&wc); 

        } else{


          (void)armci_notify_wait(right,&wc); 
          ARMCI_PutS((double*)b[me]+idx, stride, 
                     (double*)b[left]+idx, stride, count, strl, left);
#if DEBUG_ 
          printf("%d: ps=%p pd=%p i=%d idx=%d count=%d\n",me,(double*)b[me]+idx,
                      (double*)b[left]+idx,i, idx,count[0]); fflush(stdout); 
#endif
          (void)armci_notify(left);
        }

        idx += Idx; /* advance to the next slab */
    }

    ARMCI_Barrier();

    if(me==0){
       compare_patches(0.,ndim,(double*)a[0],lo,hi,dimsB,
                               (double*)b[0],lo,hi,dimsB);
       printf("OK\n");
    }

    ARMCI_Barrier();
    destroy_array(b);
    destroy_array(a);
}


int main(int argc, char* argv[])
{
int ndim;

    armci_msg_init(&argc, &argv);
    nproc = armci_msg_nproc();
    me = armci_msg_me();

    ARMCI_Init();

    armci_msg_barrier();
    if(me==0){
      printf("\nTesting armci_notify\n");
      fflush(stdout);
      sleep(1);
    }
    armci_msg_barrier();
        
    for(ndim=1; ndim<=MAXDIMS; ndim++) test_notify(ndim);
    armci_msg_barrier();

    ARMCI_Finalize();
    armci_msg_finalize();
    return(0);
}
