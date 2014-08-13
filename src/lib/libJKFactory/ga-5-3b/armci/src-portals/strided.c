#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <assert.h>

#include "armcip.h"
#include "copy.h"
#include "acc.h"
#include "memlock.h"
#include <stdio.h>
#include <assert.h>

#define DATA_SERVER_ 1

#ifdef ORNL_USE_DS_FOR_REMOTE_GETS
#define DATA_SERVER_GET_ 1
#else
#define DATA_SERVER_GET_ 0
#endif

#define ARMCI_OP_2D(op, scale, proc, src, dst, bytes, count, src_stride, dst_stride,lockit)\
if(op == GET || op ==PUT)\
      armci_copy_2D(op, proc, src, dst, bytes, count, src_stride,dst_stride);\
else if(count==1) armci_acc_1D(op, scale, proc, src, dst, bytes,lockit);\
else\
      armci_acc_2D(op, scale, proc, src, dst, bytes, count, src_stride,dst_stride,lockit) 

/* macro supports run-time selection of request sending scheme */
#if defined(CLIENT_BUF_BYPASS)
#define CAN_REQUEST_DIRECTLY _armci_bypass
#else
#  if defined(HITACHI)
#    define CAN_REQUEST_DIRECTLY 0
#  else
#    define CAN_REQUEST_DIRECTLY 1
#  endif
#endif

#define PREPROCESS_STRIDED(tmp_count) {\
 tmp_count=0;\
 if(stride_levels) \
    for(;stride_levels;stride_levels--)if(count[stride_levels]>1)break;\
 if(stride_levels&&(count[0]==src_stride_arr[0]&&count[0]==dst_stride_arr[0])){\
      tmp_count=seg_count[1];\
      count = seg_count+1;\
      seg_count[1] = seg_count[0] * seg_count[1];\
      stride_levels --;\
      src_stride_arr ++;  dst_stride_arr++ ;\
 }\
}
#define POSTPROCESS_STRIDED(tmp_count) if(tmp_count)seg_count[1]=tmp_count

#define SERVER_GET 1
#define SERVER_NBGET 2
#define DIRECT_GET 3
#define DIRECT_NBGET 4
#define SERVER_PUT 5
#define SERVER_NBPUT 6
#define DIRECT_PUT 7
#define DIRECT_NBPUT 8


#  define DO_FENCE(__proc,__prot) if(__prot==SERVER_GET);\
        else if(__prot==SERVER_PUT);\
        else if(__prot==DIRECT_GET || __prot==DIRECT_NBGET){\
          if(armci_prot_switch_fence[__proc]==SERVER_PUT)\
            ARMCI_DoFence(__proc);\
        }\
        else if(__prot==DIRECT_PUT || __prot==DIRECT_NBPUT){\
          if(armci_prot_switch_fence[__proc]==SERVER_PUT)\
            ARMCI_DoFence(__proc);\
        }\
        else;\
        armci_prot_switch_fence[__proc]=__prot

#ifndef REGIONS_REQUIRE_MEMHDL 
#  define ARMCI_MEMHDL_T void
#endif

ARMCI_MEMHDL_T *mhloc=NULL,*mhrem=NULL; 

#ifdef REGIONS_REQUIRE_MEMHDL 
   int armci_region_both_found_hndl(void *loc, void *rem, int size, int node,
                 ARMCI_MEMHDL_T **loc_memhdl,ARMCI_MEMHDL_T **rem_memhdl);
#  define ARMCI_REGION_BOTH_FOUND(_s,_d,_b,_p) \
    armci_region_both_found_hndl((_s),(_d),(_b),(_p),&mhloc,&mhrem)
#else
#  define ARMCI_REGION_BOTH_FOUND(_s,_d,_b,_p) \
    armci_region_both_found((_s),(_d),(_b),(_p))
#endif

#ifdef HAS_RDMA_GET
        
#  ifdef REGIONS_REQUIRE_MEMHDL 
   void armci_client_direct_get(int p, void *src_buf, void *dst_buf, int len,
         void** cptr,int nbtag,ARMCI_MEMHDL_T *lochdl,ARMCI_MEMHDL_T *remhdl);
#  else
   void armci_client_direct_get(int p, void *src_buf, void *dst_buf, int len,
                    void** contextptr,int nbtag,void *mhdl,void *mhdl1);
#  endif
#  define ARMCI_NBREM_GET(_p,_s,_sst,_d,_dst,_cou,_lev,_hdl) \
    armci_client_direct_get((_p),(_s),(_d),(_cou)[0],&((_hdl)->cmpl_info),(_hdl)->tag,(void *)mhloc,(void *)mhrem); \

#  define ARMCI_REM_GET(_p,_s,_sst,_d,_dst,_cou,_lev,_hdl) \
    armci_client_direct_get((_p),(_s),(_d),(_cou)[0],NULL,0,(void *)mhloc,(void *)mhrem); \

#else

#  define ARMCI_REM_GET(_p,_s,_sst,_d,_dst,_cou,_lev,_hdl) \
    armci_rem_get((_p),(_s),(_sst),(_d),(_dst),(_cou),(_lev),(_hdl),(void *)mhloc,(void *)mhrem)
#  define ARMCI_NBREM_GET ARMCI_REM_GET
        
#endif

 extern int* armci_prot_switch_fence;
 extern int armci_prot_switch_preproc;
 extern int armci_prot_switch_preop;

        
int armci_iwork[MAX_STRIDE_LEVEL];

/*\ 2-dimensional array copy
\*/
static void armci_copy_2D(int op, int proc, void *src_ptr, void *dst_ptr, 
                          int bytes, int count, int src_stride, int dst_stride)
{
    int armci_th_idx = ARMCI_THREAD_IDX;
    
#ifdef LAPI2__
#  define COUNT 1
#else
#  define COUNT count
#endif

#ifdef __crayx1
    int shmem = 1;
#else
    int shmem = SAMECLUSNODE(proc);
#endif
    if(shmem) {
      /* data is in local/shared memory -- can use memcpy */
//    printf("%s: shmem==true; count==%d\n",Portals_ID(),count);
      if(count==1){
        armci_copy(src_ptr, dst_ptr, bytes); 
//      printf("%s: shmem==true; finished\n",Portals_ID(),count);
      }else {
        char *ps=(char*)src_ptr;
        char *pd=(char*)dst_ptr;
        int j;
        for (j = 0;  j < count;  j++){
          bcopy(ps,pd,bytes);
          ps += src_stride;
          pd += dst_stride;
        }
      }
    } else {
        
        /* data not in local/shared memory-access through global address space*/
        
        if(op==PUT){ 

            printf("%s: pre UPDATE_FENCE_STATE\n",Portals_ID());
            UPDATE_FENCE_STATE(proc, PUT, COUNT);
            printf("%s: post UPDATE_FENCE_STATE\n",Portals_ID());
#ifdef LAPI
            SET_COUNTER(ack_cntr[armci_th_idx],COUNT);
#endif
            if(count==1){
                armci_put(src_ptr, dst_ptr, bytes, proc);
            }else{
                armci_put2D(proc, bytes, count, src_ptr, src_stride,
                            dst_ptr, dst_stride);
            }
            
        }else{
            
#ifdef LAPI
            SET_COUNTER(get_cntr[armci_th_idx], COUNT);
#endif
            if(count==1){
                armci_get(src_ptr, dst_ptr, bytes, proc);
            }else{
                armci_get2D(proc, bytes, count, src_ptr, src_stride,
                            dst_ptr, dst_stride);
            }
        }
    }
}


#if (defined(CRAY) && !defined(__crayx1)) || defined(FUJITSU)
#ifdef CRAY
#  define DAXPY  SAXPY
#else
#  define DAXPY  daxpy_
#endif

static int ONE=1;
#define THRESH_ACC 32

static void daxpy_2d_(void* alpha, int *rows, int *cols, void *a, int *ald,
               void* b, int *bld)
{
   int c,r;   
   double *A = (double*)a;
   double *B = (double*)b;
   double Alpha = *(double*)alpha;

   if(*rows < THRESH_ACC)
      for(c=0;c<*cols;c++)
         for(r=0;r<*rows;r++)
           A[c* *ald+ r] += Alpha * B[c* *bld+r];
   else for(c=0;c<*cols;c++)
         DAXPY(rows, alpha, B + c* *bld, &ONE, A + c* *ald, &ONE);
}
#endif


void armci_acc_1D(int op, void *scale, int proc, void *src, void *dst, int bytes, int lockit)
{
int rows;
void (ATR *func)(void*, void*, void*, int*);
    ARMCI_PR_DBG("enter",0);
      switch (op){
      case ARMCI_ACC_INT:
          rows = bytes/sizeof(int);
          func = I_ACCUMULATE_1D;
          break;
      case ARMCI_ACC_LNG:
          rows = bytes/sizeof(long);
          func = L_ACCUMULATE_1D;
          break;
      case ARMCI_ACC_DBL:
          rows = bytes/sizeof(double);
          func = D_ACCUMULATE_1D;
          break;
      case ARMCI_ACC_DCP:
          rows = bytes/(2*sizeof(double));
          func = Z_ACCUMULATE_1D;
          break;
      case ARMCI_ACC_CPL:
          rows = bytes/(2*sizeof(float));
          func = C_ACCUMULATE_1D;
          break;
      case ARMCI_ACC_FLT:
          rows = bytes/sizeof(float);
          func = F_ACCUMULATE_1D;
          break;
      default: armci_die("ARMCI accumulate: operation not supported",op);
          func = F_ACCUMULATE_1D; /*avoid compiler whining */
      }


      if(lockit){
          ARMCI_LOCKMEM(dst, bytes + (char*)dst, proc);
      }
      func(scale, dst, src, &rows);
      if(lockit)ARMCI_UNLOCKMEM(proc);
    ARMCI_PR_DBG("exit",0);
}

/*\ 2-dimensional accumulate
\*/
void armci_acc_2D(int op, void* scale, int proc, void *src_ptr, void *dst_ptr,
                  int bytes, int cols, int src_stride, int dst_stride, int lockit)
{
int   rows, lds, ldd, span;
void (ATR *func)(void*, int*, int*, void*, int*, void*, int*);

    ARMCI_PR_DBG("enter",0);

/*
      if((long)src_ptr%ALIGN)armci_die("src not aligned",(long)src_ptr);
      if((long)dst_ptr%ALIGN)armci_die("src not aligned",(long)dst_ptr);
*/

      switch (op){
      case ARMCI_ACC_INT:
          rows = bytes/sizeof(int);
          ldd  = dst_stride/sizeof(int);
          lds  = src_stride/sizeof(int);
          func = I_ACCUMULATE_2D;
          break;
      case ARMCI_ACC_LNG:
          rows = bytes/sizeof(long);
          ldd  = dst_stride/sizeof(long);
          lds  = src_stride/sizeof(long);
          func = L_ACCUMULATE_2D;
          break;
      case ARMCI_ACC_DBL:
          rows = bytes/sizeof(double);
          ldd  = dst_stride/sizeof(double);
          lds  = src_stride/sizeof(double);
          func = D_ACCUMULATE_2D;
          break;
      case ARMCI_ACC_DCP:
          rows = bytes/(2*sizeof(double));
          ldd  = dst_stride/(2*sizeof(double));
          lds  = src_stride/(2*sizeof(double));
          func = Z_ACCUMULATE_2D;
          break;
      case ARMCI_ACC_CPL:
          rows = bytes/(2*sizeof(float));
          ldd  = dst_stride/(2*sizeof(float));
          lds  = src_stride/(2*sizeof(float));
          func = C_ACCUMULATE_2D;
          break;
      case ARMCI_ACC_FLT:
          rows = bytes/sizeof(float);
          ldd  = dst_stride/sizeof(float);
          lds  = src_stride/sizeof(float);
          func = F_ACCUMULATE_2D;
          break;
      case ARMCI_ACC_RA:
          rows = bytes/sizeof(long);
          ldd  = dst_stride/sizeof(long);
          lds  = src_stride/sizeof(long);
          func = RA_ACCUMULATE_2D;
          break;
      default: armci_die("ARMCI accumulate: operation not supported",op);
          func = F_ACCUMULATE_2D; /*avoid compiler whining */
      }

             
      if(lockit){ 
          span = cols*dst_stride;
          ARMCI_LOCKMEM(dst_ptr, span + (char*)dst_ptr, proc);
      }
      func(scale, &rows, &cols, dst_ptr, &ldd, src_ptr, &lds);
      if(lockit)ARMCI_UNLOCKMEM(proc);
    ARMCI_PR_DBG("exit",0);

}


/*\ compute range of strided data AND lock it
\*/
static void 
armci_lockmem_patch(void* dst_ptr, int dst_stride_arr[], int count[], int stride_levels, int proc)
{
    long span = count[stride_levels];
    ARMCI_PR_DBG("enter",0);
    span *= dst_stride_arr[stride_levels-1];

    /* lock region of remote memory */
    ARMCI_LOCKMEM(dst_ptr, span + (char*)dst_ptr, proc);
    ARMCI_PR_DBG("exit",0);
}


/*\ strided accumulate on top of remote memory copy:
 *  copies remote data to local buffer, accumulates, puts it back 
 *  Note: if we are here then remote patch must fit in the ARMCI buffer
\*/
int armci_acc_copy_strided(int optype, void* scale, int proc,
                                  void* src_ptr, int src_stride_arr[],  
		                  void* dst_ptr, int dst_stride_arr[], 
                                  int count[], int stride_levels)
{
    void *buf_ptr = armci_internal_buffer;
    int  rc, i, *buf_stride_arr = armci_iwork;
    ARMCI_PR_DBG("enter",0);
    armci_lockmem_patch(dst_ptr,dst_stride_arr, count, stride_levels, proc);

    /* setup stride array for internal buffer */
    buf_stride_arr[0]=count[0];
    for(i=0; i< stride_levels; i++) {
         buf_stride_arr[i+1]= buf_stride_arr[i]*count[i+1];
    }

    /* get remote data to local buffer */
    rc = armci_op_strided(GET, scale, proc, dst_ptr, dst_stride_arr, buf_ptr, 
                          buf_stride_arr, count, stride_levels, 0,NULL);

    if(rc) { ARMCI_UNLOCKMEM(proc); return(rc); }

    /* call local accumulate with lockit=0 (we locked it already) and proc=me */
    rc = armci_op_strided(optype, scale, armci_me, src_ptr, src_stride_arr, 
                          buf_ptr,buf_stride_arr, count, stride_levels,0,NULL);
    if(rc) { ARMCI_UNLOCKMEM(proc); return(rc); }

    /* put data back from the buffer to remote location */
    rc = armci_op_strided(PUT, scale, proc, buf_ptr, buf_stride_arr, dst_ptr, 
                          dst_stride_arr, count, stride_levels,0,NULL);

    FENCE_NODE(proc); /* make sure put completes before unlocking */
    ARMCI_UNLOCKMEM(proc);    /* release memory lock */
    ARMCI_PR_DBG("exit",0);

    return(rc);
}



/*\ Strided  operation
\*/
int armci_op_strided(int op, void* scale, int proc,void *src_ptr, 
                     int src_stride_arr[], void* dst_ptr, int dst_stride_arr[], 
                     int count[], int stride_levels, int lockit,
                     armci_ihdl_t nb_handle)
{
char *src = (char*)src_ptr, *dst=(char*)dst_ptr;
int s2, s3, i,j, unlockit=0;
int total_of_2D;
int index[MAX_STRIDE_LEVEL], unit[MAX_STRIDE_LEVEL];
    ARMCI_PR_DBG("enter",op);
#   if defined(ACC_COPY)
      
#      ifdef ACC_SMP
         if(ARMCI_ACC(op) && !(SAMECLUSNODE(proc)) )
#      else
         if ( ARMCI_ACC(op) && proc!=armci_me)
#      endif
             /* copy remote data, accumulate, copy back*/
             return (armci_acc_copy_strided(op,scale, proc, src_ptr, src_stride_arr,
                                       dst_ptr, dst_stride_arr, count, stride_levels));

         else; /* do it directly through shared/local memory */
#   endif


      if(ARMCI_ACC(op) && (stride_levels>2) && lockit){
         /* we need one lock operation only - must be done outside 2d acc */
         armci_lockmem_patch(dst_ptr,dst_stride_arr, count, stride_levels, proc);
         unlockit=1;
         lockit =0;
      }
/*    if(proc!=armci_me) INTR_OFF;*/
    if(armci_me>=0 && !SAMECLUSNODE(proc)) {
       printf("%s network_strided not supported (in op_strided)\n",Portals_ID());
       abort();
       armci_network_strided(op,scale,proc,src_ptr,src_stride_arr,dst_ptr,
                         dst_stride_arr,count,stride_levels,nb_handle);
    }
      else {
//     printf("%s in large switch stmt in op_strided (stride_levels=%d)\n",Portals_ID(),stride_levels);
    switch (stride_levels) {
      case 0: /* 1D copy */ 

          ARMCI_OP_2D(op, scale, proc, src_ptr, dst_ptr, count[0], 1, 
                      count[0], count[0], lockit); 
          
          break;
          
      case 1: /* 2D op */
          ARMCI_OP_2D(op, scale, proc, src_ptr, dst_ptr, count[0], count[1], 
                      src_stride_arr[0], dst_stride_arr[0], lockit);
          break;

      case 2: /* 3D op */
          for (s2= 0; s2  < count[2]; s2++){ /* 2D copy */
              ARMCI_OP_2D(op, scale, proc, src+s2*src_stride_arr[1], 
                          dst+s2*dst_stride_arr[1], count[0], count[1], 
                          src_stride_arr[0], dst_stride_arr[0], lockit );
          }
          break;
          
      case 3: /* 4D op */
          for(s3=0; s3< count[3]; s3++){
              src = (char*)src_ptr + src_stride_arr[2]*s3;
              dst = (char*)dst_ptr + dst_stride_arr[2]*s3;
              for (s2= 0; s2  < count[2]; s2++){ /* 3D copy */
                  ARMCI_OP_2D(op, scale, proc, src+s2*src_stride_arr[1],
                              dst+s2*dst_stride_arr[1],
                              count[0], count[1],src_stride_arr[0],
                              dst_stride_arr[0],lockit);
              }
          }
          break;
          
      default: /* N-dimensional */ 
      {
	  /* stride_levels is not the same as ndim. it is ndim-1
	   * For example a 10x10x10... array, suppose the datatype is byte
	   * the stride_arr is 10, 10x10, 10x10x10 ....
	   */
          index[2] = 0; unit[2] = 1; total_of_2D = count[2];
          for(j=3; j<=stride_levels; j++) {
              index[j] = 0; unit[j] = unit[j-1] * count[j-1];
              total_of_2D *= count[j];
          }

          for(i=0; i<total_of_2D; i++) {
              src = (char *)src_ptr; dst = (char *)dst_ptr;
              for(j=2; j<=stride_levels; j++) {
                  src += index[j] * src_stride_arr[j-1];
                  dst += index[j] * dst_stride_arr[j-1];
                  
                  if(((i+1) % unit[j]) == 0) index[j]++;
                  if(index[j] >= count[j]) index[j] = 0;
              }
              
              ARMCI_OP_2D(op, scale, proc, src, dst, count[0], count[1], 
                          src_stride_arr[0], dst_stride_arr[0], lockit);
          }
          
      }
    }
    } // ends else block

//  printf("%s after switch stmt; prior to fence/lock\n",Portals_ID());

    if(unlockit){
#      if defined(ACC_COPY)
          FENCE_NODE(proc); 
#      endif
       ARMCI_UNLOCKMEM(proc);    /* release memory lock */
    }

//  printf("%s after fence/lock; leaving op_strided\n",Portals_ID());
    ARMCI_PR_DBG("exit",op);
    return 0;
}


int PARMCI_PutS( void *src_ptr,        /* pointer to 1st segment at source*/ 
		int src_stride_arr[], /* array of strides at source */
		void* dst_ptr,        /* pointer to 1st segment at destination*/
		int dst_stride_arr[], /* array of strides at destination */
		int seg_count[],      /* number of segments at each stride 
                                         levels: count[0]=bytes*/
		int stride_levels,    /* number of stride levels */
                int proc              /* remote process(or) ID */
                )
{
int rc=0, direct=1;
int *count=seg_count, tmp_count=0;

    ARMCI_PR_DBG("enter",proc);
    if(src_ptr == NULL || dst_ptr == NULL) return FAIL;
    if(count[0]<0)return FAIL3;
    if(stride_levels <0 || stride_levels > MAX_STRIDE_LEVEL) return FAIL4;
    if(proc<0)return FAIL5;

    ORDER(PUT,proc); /* ensure ordering */
    PREPROCESS_STRIDED(tmp_count);

#if DATA_SERVER_
    if(stride_levels) direct=SAMECLUSNODE(proc);
    direct=SAMECLUSNODE(proc);
#endif

    // printf("%s direct=%d, proc=%d\n",Portals_ID(),direct,proc);

    if(!direct){
       DO_FENCE(proc,SERVER_PUT);
//     printf("%s calling pack_strided in PARMCI_PutS\n",Portals_ID());
       rc = armci_pack_strided(PUT, NULL, proc, src_ptr, src_stride_arr,dst_ptr,
                  dst_stride_arr, count, stride_levels, NULL, -1, -1, -1,NULL);
    }
    else
    {
       if(!SAMECLUSNODE(proc))DO_FENCE(proc,DIRECT_PUT);
//     printf("%s calling op_strided in PARMCI_PutS\n",Portals_ID());
       rc = armci_op_strided( PUT, NULL, proc, src_ptr, src_stride_arr, 
				 dst_ptr, dst_stride_arr,count,stride_levels, 
				 0,NULL);
    }
    POSTPROCESS_STRIDED(tmp_count);

    ARMCI_PR_DBG("exit",proc);
    if(rc) return FAIL6;
    else return 0;

}

int PARMCI_PutS_flag(
      void* src_ptr,        /* pointer to 1st segment at source */
      int src_stride_arr[], /* array of strides at source */
      void* dst_ptr,        /* pointer to 1st segment at destination */
      int dst_stride_arr[], /* array of strides at destination */
      int count[],          /* number of units at each stride level,
                               count[0] = #bytes */
      int stride_levels,    /* number of stride levels */
      int *flag,            /* pointer to remote flag */
      int val,              /* value to set flag upon completion of
                               data transfer */
      int proc              /* remote process(or) ID */
      )
{
  int bytes;
  /* Put local data on remote processor */
  PARMCI_PutS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr,
             count, stride_levels, proc);

  /* Send signal to remote processor that data transfer has
   * been completed. */
  bytes = sizeof(int);
  ARMCI_Put(&val, flag, bytes, proc);
  return 1;
}


int PARMCI_Put_flag(void *src, void* dst,int bytes,int *f,int v,int proc) {
  return  PARMCI_PutS_flag(src, NULL, dst, NULL, &bytes, 0, f, v, proc);
}


int PARMCI_PutS_flag_dir(void *src_ptr,   int src_stride_arr[],
            void* dst_ptr,   int dst_stride_arr[],
            int seg_count[], int stride_levels,
            int *flag, int val, int proc) {
  return PARMCI_PutS_flag(src_ptr, src_stride_arr,dst_ptr,dst_stride_arr,
             seg_count, stride_levels, flag, val, proc);
}


int PARMCI_GetS( void *src_ptr,  	/* pointer to 1st segment at source*/ 
		int src_stride_arr[],   /* array of strides at source */
		void* dst_ptr,          /* 1st segment at destination*/
		int dst_stride_arr[],   /* array of strides at destination */
		int seg_count[],       /* number of segments at each stride 
					   levels: count[0]=bytes*/
		int stride_levels,      /* number of stride levels */
                int proc                /* remote process(or) ID */
                )
{
int rc,direct=1;
int *count=seg_count, tmp_count=0;
    ARMCI_PR_DBG("enter",proc);

    if(src_ptr == NULL || dst_ptr == NULL) return FAIL;
    if(seg_count[0]<0)return FAIL3;
    if(stride_levels <0 || stride_levels > MAX_STRIDE_LEVEL) return FAIL4;
    if(proc<0||proc>=armci_nproc){printf("\n%d:%s:proc=%d",armci_me,FUNCTION_NAME,proc);fflush(stdout);return FAIL5;}
    
    ORDER(GET,proc); /* ensure ordering */
    PREPROCESS_STRIDED(tmp_count);

#if DATA_SERVER_GET_
    if(stride_levels)direct=SAMECLUSNODE(proc);
    direct=SAMECLUSNODE(proc);
#endif
    if(!direct){
       DO_FENCE(proc,SERVER_GET);
       rc = armci_pack_strided(GET, NULL, proc, src_ptr, src_stride_arr,
                                 dst_ptr,dst_stride_arr,count,stride_levels,
                                 NULL,-1,-1,-1,NULL);
               
    }else{
       if(!SAMECLUSNODE(proc))DO_FENCE(proc,DIRECT_GET);
       rc = armci_op_strided(GET, NULL, proc, src_ptr, src_stride_arr, dst_ptr,
                             dst_stride_arr,count, stride_levels,0,NULL);
    }

    POSTPROCESS_STRIDED(tmp_count);
    ARMCI_PR_DBG("exit",proc);
    if(rc) return FAIL6;
    else return 0;
}




int PARMCI_AccS( int  optype,            /* operation */
                void *scale,            /* scale factor x += scale*y */
                void *src_ptr,          /* pointer to 1st segment at source*/ 
		int src_stride_arr[],   /* array of strides at source */
		void* dst_ptr,          /* 1st segment at destination*/
		int dst_stride_arr[],   /* array of strides at destination */
		int seg_count[],        /* number of segments at each stride 
                                           levels: count[0]=bytes*/
		int stride_levels,      /* number of stride levels */
                int proc                /* remote process(or) ID */
                )
{
int rc, direct=1;
int *count=seg_count, tmp_count=0;

    ARMCI_PR_DBG("enter",proc);
    if(src_ptr == NULL || dst_ptr == NULL) return FAIL;
    if(src_stride_arr == NULL || dst_stride_arr ==NULL) return FAIL2;
    if(count[0]<0)return FAIL3;
    if(stride_levels <0 || stride_levels > MAX_STRIDE_LEVEL) return FAIL4;
    if(proc<0)return FAIL5;

    ORDER(optype,proc); /* ensure ordering */
    PREPROCESS_STRIDED(tmp_count);

#if DATA_SERVER_
    direct=SAMECLUSNODE(proc);
#endif

#   if defined(ACC_COPY) && !defined(ACC_SMP)
       if(armci_me != proc) direct=0;
#   endif
    if(direct)
      rc = armci_op_strided(optype,scale, proc, src_ptr, src_stride_arr,dst_ptr,
                           dst_stride_arr, count, stride_levels,1,NULL);
    else{
      DO_FENCE(proc,SERVER_PUT);
      rc = armci_pack_strided(optype,scale,proc,src_ptr, src_stride_arr,dst_ptr,
                      dst_stride_arr,count,stride_levels,NULL,-1,-1,-1,NULL);
    }
    POSTPROCESS_STRIDED(tmp_count);
    ARMCI_PR_DBG("exit",proc);
    if(rc) return FAIL6;
    else return 0;
}


/* 
   whatever original put and get functions were here have been
   replaced with the proper ones from the main armci branch.
   the old functions were entirely responsible for causing the
   test_vector_acc test to fail in test.x
*/
    
int PARMCI_Put(void *src, void* dst, int bytes, int proc) {
  int rc=0;
//ARMCI_PROFILE_START_STRIDED(&bytes, 0, proc, ARMCI_PROF_PUT);
  rc = PARMCI_PutS(src, NULL, dst, NULL, &bytes, 0, proc);
//ARMCI_PROFILE_STOP_STRIDED(ARMCI_PROF_PUT);
  assert(rc==0);
  return rc;
}

int PARMCI_Get(void *src, void* dst, int bytes, int proc) {
  int rc=0;
//ARMCI_PROFILE_START_STRIDED(&bytes, 0, proc, ARMCI_PROF_GET);

#ifdef __crayx1
  memcpy(dst,src,bytes);
#else
  rc = PARMCI_GetS(src, NULL, dst, NULL, &bytes, 0, proc);
#endif
//ARMCI_PROFILE_STOP_STRIDED(ARMCI_PROF_GET);
//dassert(1,rc==0);
  assert(rc==0);
  return rc;
}

int PARMCI_Acc(int optype, void *scale, void *src, void* dst, int bytes, int proc) {
  int rc=0;
  rc = PARMCI_AccS(optype, scale, src, NULL, dst, NULL, &bytes, 0, proc);
  return rc;
}

#define PACK1D 1

#if PACK1D 
#  define armci_read_strided1  armci_read_strided
#  define armci_write_strided1 armci_write_strided
#else
#  define armci_read_strided2  armci_read_strided
#  define armci_write_strided2 armci_write_strided
#endif

void armci_write_strided1(void *ptr, int stride_levels, int stride_arr[],
                   int count[], char *buf)
{
    int i, j;
    long idx;    /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int bvalue[MAX_STRIDE_LEVEL], bunit[MAX_STRIDE_LEVEL];
    int bytes = count[0];
    ARMCI_PR_DBG("enter",stride_levels);

    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++)
        n1dim *= count[i];

    /* calculate the destination indices */
    bvalue[0] = 0; bvalue[1] = 0; bunit[0] = 1; bunit[1] = 1;
    for(i=2; i<=stride_levels; i++) {
        bvalue[i] = 0;
        bunit[i] = bunit[i-1] * count[i-1];
    }

    for(i=0; i<n1dim; i++) {
        idx = 0;
        for(j=1; j<=stride_levels; j++) {
            idx += bvalue[j] * stride_arr[j-1];
            if((i+1) % bunit[j] == 0) bvalue[j]++;
            if(bvalue[j] > (count[j]-1)) bvalue[j] = 0;
        }

        armci_copy( ((char*)ptr)+idx, buf, bytes);
        buf += count[0];
    }
    ARMCI_PR_DBG("exit",stride_levels);
}


void armci_write_strided2(void *ptr, int stride_levels, int stride_arr[],
                          int count[], char *buf)
{                  
    int i, j;
    int total;   /* number of 2 dim block */
    int index[MAX_STRIDE_LEVEL], unit[MAX_STRIDE_LEVEL];
    ARMCI_PR_DBG("enter",stride_levels);
    
    if(stride_levels == 0){
       armci_copy( ptr, buf, count[0]);
    }else if (count[0]%ALIGN_SIZE || (unsigned long)ptr%ALIGN_SIZE ) 
              armci_write_strided1(ptr,stride_levels, stride_arr,count,buf);
          else {
             int rows, ld, idx, ldd;
             char *src;
             rows = count[0]/8;
             ld   = stride_arr[0]/8;
             switch(stride_levels){
             case 1: 
                     DCOPY21(&rows, count+1, ptr, &ld, buf, &idx);
                     break;
             case 2: 
                     ldd = stride_arr[1]/stride_arr[0];
                     DCOPY31(&rows, count+1, count+2, ptr, &ld, &ldd, buf,&idx);

                     break;
             default: 
                     index[2] = 0; unit[2] = 1; total = count[2];
                     for(j=3; j<=stride_levels; j++) {
                         index[j] = 0; unit[j] = unit[j-1] * count[j-1];
                         total *= count[j];
                     }
                     for(i=0; i<total; i++) {
                         src = (char *)ptr; 
                         for(j=2; j<=stride_levels; j++) {
                             src += index[j] * stride_arr[j-1];
                             if(((i+1) % unit[j]) == 0) index[j]++;
                             if(index[j] >= count[j]) index[j] = 0;
                         }
                         DCOPY21(&rows, count+1,src, &ld, buf, &idx); 
                         buf = (char*) ((double*)buf + idx);
                     }
            } /*switch */
         } /*else */
    ARMCI_PR_DBG("exit",stride_levels);
}


void armci_read_strided1(void *ptr, int stride_levels, int stride_arr[],
                        int count[], char *buf)
{
    int i, j;
    long idx;    /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int bvalue[MAX_STRIDE_LEVEL], bunit[MAX_STRIDE_LEVEL];
    int bytes = count[0];

    ARMCI_PR_DBG("enter",stride_levels);
    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++)
        n1dim *= count[i];

    /* calculate the destination indices */
    bvalue[0] = 0; bvalue[1] = 0; bunit[0] = 1; bunit[1] = 1;
    for(i=2; i<=stride_levels; i++) {
        bvalue[i] = 0;
        bunit[i] = bunit[i-1] * count[i-1];
    }

    for(i=0; i<n1dim; i++) {
        idx = 0;
        for(j=1; j<=stride_levels; j++) {
            idx += bvalue[j] * stride_arr[j-1];
            if((i+1) % bunit[j] == 0) bvalue[j]++;
            if(bvalue[j] > (count[j]-1)) bvalue[j] = 0;
        }

        armci_copy(buf, ((char*)ptr)+idx,bytes);
        buf += count[0];
    }
    ARMCI_PR_DBG("exit",stride_levels);
}


void armci_read_strided2(void *ptr, int stride_levels, int stride_arr[],
                         int count[], char *buf)
{                  
    int i, j;
    int total;   /* number of 2 dim block */
    int index[MAX_STRIDE_LEVEL], unit[MAX_STRIDE_LEVEL];
   
    ARMCI_PR_DBG("enter",stride_levels);
    if(stride_levels == 0){
       armci_copy( buf, ptr, count[0]);
    }else if (count[0]%ALIGN_SIZE || (unsigned long)ptr%ALIGN_SIZE) 
              armci_read_strided1(ptr,stride_levels, stride_arr,count,buf);
          else {
             int rows, ld, idx, ldd;
             char *src;
             rows = count[0]/8;
             ld   = stride_arr[0]/8;
             switch(stride_levels){
             case 1: 
                     DCOPY12(&rows, count+1, ptr, &ld, buf, &idx);
                     break;
             case 2:
                     ldd = stride_arr[1]/stride_arr[0];   
                     DCOPY13(&rows, count+1, count+2, ptr, &ld, &ldd, buf,&idx);
                     break;
             default:
                     index[2] = 0; unit[2] = 1; total = count[2];
                     for(j=3; j<=stride_levels; j++) {
                         index[j] = 0; unit[j] = unit[j-1] * count[j-1];
                         total *= count[j];
                     }
                     for(i=0; i<total; i++) {
                         src = (char *)ptr; 
                         for(j=2; j<=stride_levels; j++) {
                             src += index[j] * stride_arr[j-1];
                             if(((i+1) % unit[j]) == 0) index[j]++;
                             if(index[j] >= count[j]) index[j] = 0;
                         }
                         DCOPY12(&rows, count+1,src, &ld, buf, &idx);
                         buf = (char*) ((double*)buf + idx);
                     }
            } /*switch */
         } /*else */
    ARMCI_PR_DBG("exit",stride_levels);
}

/*\Non-Blocking API
\*/
int PARMCI_NbPutS( void *src_ptr,        /* pointer to 1st segment at source*/ 
		int src_stride_arr[], /* array of strides at source */
		void* dst_ptr,        /* pointer to 1st segment at destination*/
		int dst_stride_arr[], /* array of strides at destination */
		int seg_count[],      /* number of segments at each stride 
                                         levels: count[0]=bytes*/
		int stride_levels,    /* number of stride levels */
                int proc,             /* remote process(or) ID */
                armci_hdl_t* usr_hdl  /* armci non-blocking call handle*/
                )
{
armci_ihdl_t nb_handle = (armci_ihdl_t)usr_hdl;
int *count=seg_count, tmp_count=0;
int rc=0, direct=1;
    ARMCI_PR_DBG("enter",proc);
    if(src_ptr == NULL || dst_ptr == NULL) return FAIL;
    if(count[0]<0)return FAIL3;
    if(stride_levels <0 || stride_levels > MAX_STRIDE_LEVEL) return FAIL4;
    if(proc<0)return FAIL5;

    PREPROCESS_STRIDED(tmp_count);

#if DATA_SERVER_
    if(stride_levels)direct=SAMECLUSNODE(proc);
    direct=SAMECLUSNODE(proc);
#endif

    /* aggregate put */
    if(nb_handle && nb_handle->agg_flag == SET) {
      if(!direct){ 
	rc= armci_agg_save_strided_descriptor(src_ptr, src_stride_arr, 
						 dst_ptr, dst_stride_arr, 
						 count, stride_levels, proc, 
						 PUT, nb_handle);
        POSTPROCESS_STRIDED(tmp_count);
        return(rc);
      }
    } 
    else {
      UPDATE_FENCE_INFO(proc);
      
      /*set tag and op in the nb handle*/
      if(nb_handle){
	nb_handle->tag = GET_NEXT_NBTAG();
	nb_handle->op  = PUT;
	nb_handle->proc= proc;
	nb_handle->bufid=NB_NONE;
      }
      else
        nb_handle = armci_set_implicit_handle(PUT, proc);
    }

    if(!direct){
      DO_FENCE(proc,SERVER_NBPUT);
      rc = armci_pack_strided(PUT, NULL, proc, src_ptr, src_stride_arr,dst_ptr,
                  dst_stride_arr, count, stride_levels,NULL,-1,-1,-1,nb_handle);
    }
    else{
      if(!SAMECLUSNODE(proc))DO_FENCE(proc,DIRECT_NBPUT);
      rc = armci_op_strided( PUT, NULL, proc, src_ptr, src_stride_arr,
                      dst_ptr,dst_stride_arr,count,stride_levels, 0,nb_handle);
    }
    
    POSTPROCESS_STRIDED(tmp_count);
    ARMCI_PR_DBG("exit",proc);
    if(rc) return FAIL6;
    else return 0;
}

int PARMCI_NbGetS( void *src_ptr,  	/* pointer to 1st segment at source*/ 
		int src_stride_arr[],   /* array of strides at source */
		void* dst_ptr,          /* 1st segment at destination*/
		int dst_stride_arr[],   /* array of strides at destination */
		int seg_count[],        /* number of segments at each stride 
                                           levels: byte_count[0]=bytes*/
		int stride_levels,      /* number of stride levels */
                int proc,               /* remote process(or) ID */
                armci_hdl_t* usr_hdl  /* armci non-blocking call handle*/
                )
{
armci_ihdl_t nb_handle = (armci_ihdl_t)usr_hdl;
int rc=0,direct=1;
int *count=seg_count, tmp_count=0;

    ARMCI_PR_DBG("enter",proc);
    if(src_ptr == NULL || dst_ptr == NULL) return FAIL;
    if(seg_count[0]<0)return FAIL3;
    if(stride_levels <0 || stride_levels > MAX_STRIDE_LEVEL) return FAIL4;
    if(proc<0)return FAIL5;

#if DATA_SERVER_GET_
    if(stride_levels)direct=SAMECLUSNODE(proc);
    direct=SAMECLUSNODE(proc);
#endif

    PREPROCESS_STRIDED(tmp_count);

    if(nb_handle && nb_handle->agg_flag == SET) {
      if(!direct){ 
	rc= armci_agg_save_strided_descriptor(src_ptr, src_stride_arr,
					 dst_ptr, dst_stride_arr, 
					 count, stride_levels, proc, 
					 GET, nb_handle);
        POSTPROCESS_STRIDED(tmp_count);
        return(rc);
      }
    } 
    else {
      /* ORDER(GET,proc); ensure ordering */
      /*set tag and op in the nb handle*/
      if(nb_handle){
	nb_handle->tag = GET_NEXT_NBTAG();
	nb_handle->op  = GET;
	nb_handle->proc= proc;
	nb_handle->bufid=NB_NONE;
      }
      else
        nb_handle = armci_set_implicit_handle(GET, proc);
    }
    
    if(!direct){
       DO_FENCE(proc,SERVER_NBGET);
          rc = armci_pack_strided(GET, NULL, proc, src_ptr, src_stride_arr,
                                 dst_ptr,dst_stride_arr,count,stride_levels,
                                 NULL,-1,-1,-1,nb_handle);
    }
    else{
       if(!SAMECLUSNODE(proc))DO_FENCE(proc,DIRECT_GET);
       rc = armci_op_strided(GET, NULL, proc, src_ptr, src_stride_arr, dst_ptr,
                             dst_stride_arr,count, stride_levels,0,nb_handle);
    }

    POSTPROCESS_STRIDED(tmp_count);

    ARMCI_PR_DBG("exit",proc);
    if(rc) return FAIL6;
    else return 0;
}


int PARMCI_NbAccS( int  optype,            /* operation */
                void *scale,            /* scale factor x += scale*y */
                void *src_ptr,          /* pointer to 1st segment at source*/ 
		int src_stride_arr[],   /* array of strides at source */
		void* dst_ptr,          /* 1st segment at destination*/
		int dst_stride_arr[],   /* array of strides at destination */
		int seg_count[],        /* number of segments at each stride 
                                           levels: count[0]=bytes*/
		int stride_levels,      /* number of stride levels */
                int proc,               /* remote process(or) ID */
                armci_hdl_t* usr_hdl  /* armci non-blocking call handle*/
                )
{
armci_ihdl_t nb_handle = (armci_ihdl_t)usr_hdl;
int *count=seg_count, tmp_count=0;
int rc, direct=1;

    ARMCI_PR_DBG("enter",proc);
    if(src_ptr == NULL || dst_ptr == NULL) return FAIL;
    if(src_stride_arr == NULL || dst_stride_arr ==NULL) return FAIL2;
    if(count[0]<0)return FAIL3;
    if(stride_levels <0 || stride_levels > MAX_STRIDE_LEVEL) return FAIL4;
    if(proc<0)return FAIL5;

    UPDATE_FENCE_INFO(proc);
    PREPROCESS_STRIDED(tmp_count);
    
#if DATA_SERVER_
    direct=SAMECLUSNODE(proc);
#endif

#   if defined(ACC_COPY) && !defined(ACC_SMP)
       if(armci_me != proc) direct=0;
#   endif
 
    /*set tag and op in the nb handle*/
    if(nb_handle){
      nb_handle->tag = GET_NEXT_NBTAG();
      nb_handle->op  = optype;
      nb_handle->proc= proc;
      nb_handle->bufid=NB_NONE;
    }
    else
      nb_handle = armci_set_implicit_handle(optype, proc);


    if(direct){
      rc = armci_op_strided(optype,scale, proc, src_ptr, src_stride_arr,dst_ptr,
			    dst_stride_arr, count, stride_levels,1,NULL);
    }
    else{
      DO_FENCE(proc,SERVER_NBPUT);
      rc = armci_pack_strided(optype,scale,proc,src_ptr, src_stride_arr,dst_ptr,
                    dst_stride_arr,count,stride_levels,NULL,-1,-1,-1,nb_handle);
    }

    POSTPROCESS_STRIDED(tmp_count);

    ARMCI_PR_DBG("exit",proc);
    if(rc) return FAIL6;
    else return 0;
}


#if !defined(ACC_COPY)&&!defined(CRAY_YMP)&&!defined(CYGNUS)&&!defined(CYGWIN) &&!defined(BGML)
#   define REMOTE_OP
#endif

#define INIT_NB_HANDLE(nb,o,p) if(nb){\
             (nb)->tag = 0;\
             (nb)->op  = (o); (nb)->proc= (p);\
             (nb)->bufid=NB_NONE;}\
             else { (nb)=armci_set_implicit_handle(o, p); (nb)->tag=0; }

void set_nbhandle(armci_ihdl_t *nbh, armci_hdl_t *nb_handle, int op,
                         int proc)
{
   if(nb_handle)
   {
      *nbh=(armci_ihdl_t)nb_handle;
   }
   else
   {
      *nbh=armci_set_implicit_handle(op, proc);
   }
}


int PARMCI_NbPut(void *src, void* dst, int bytes, int proc,armci_hdl_t* uhandle)
{

int rc=0, direct=0;
armci_ihdl_t nb_handle = (armci_ihdl_t)uhandle;
    ARMCI_PR_DBG("enter",proc);
    
    if(src == NULL || dst == NULL) return FAIL;

    direct =SAMECLUSNODE(proc);

    /* aggregate put */
    if(nb_handle && nb_handle->agg_flag == SET) {
      if(direct) { armci_copy(src,dst,bytes); rc=0; }
      else
	rc=armci_agg_save_descriptor(src,dst,bytes,proc,PUT,0,nb_handle); 
      return rc;
    }

    if(direct) {
      /*armci_wait needs proc to compute direct*/
      INIT_NB_HANDLE(nb_handle,PUT,proc);
      armci_copy(src,dst,bytes);
    }
    else{
    # ifdef PORTALS
      rc=PARMCI_NbPutS(src, NULL,dst,NULL, &bytes,0,proc,uhandle);
    # else
#     ifdef ARMCI_NB_PUT
      INIT_NB_HANDLE(nb_handle,PUT,proc);
      UPDATE_FENCE_STATE(proc, PUT, 1);
      ARMCI_NB_PUT(src, dst, bytes, proc, &nb_handle->cmpl_info);
#     else
      rc=PARMCI_NbPutS(src, NULL,dst,NULL, &bytes,0,proc,uhandle);
#     endif
    # endif
    }

    ARMCI_PR_DBG("exit",proc);
    return(rc);
}


int PARMCI_NbGet(void *src, void* dst, int bytes, int proc,armci_hdl_t* uhandle)
{

int rc=0, direct=0;
armci_ihdl_t nb_handle = (armci_ihdl_t)uhandle;
    ARMCI_PR_DBG("enter",proc);
    
    if(src == NULL || dst == NULL) return FAIL;

    direct =SAMECLUSNODE(proc);

    if(nb_handle && nb_handle->agg_flag == SET) {
      if(direct) { armci_copy(src,dst,bytes); rc=0; }
      else
	rc=armci_agg_save_descriptor(src,dst,bytes,proc,GET,0,nb_handle);
      return rc;
    }

    if(direct) {
      /*armci_wait needs proc to compute direct*/
      INIT_NB_HANDLE(nb_handle,PUT,proc);
      armci_copy(src,dst,bytes);
    }else{
    
    # ifdef PORTALS
      rc=PARMCI_NbGetS(src, NULL,dst,NULL, &bytes,0,proc,uhandle);
    # else
#     ifdef ARMCI_NB_GET
      /*set tag and op in the nb handle*/
      INIT_NB_HANDLE(nb_handle,GET,proc);
      
      ARMCI_NB_GET(src, dst, bytes, proc, &nb_handle->cmpl_info);
#     else
      rc=PARMCI_NbGetS(src, NULL,dst,NULL, &bytes,0,proc,uhandle);
#     endif
    # endif
    }
    ARMCI_PR_DBG("exit",proc);
    return(rc);
}


static void _armci_rem_value(int op, void *src, void *dst, int proc, 
			     int bytes) {  
    int rc=0;
    int armci_th_idx = ARMCI_THREAD_IDX;

    ORDER(op,proc); /* ensure ordering */

#if defined(REMOTE_OP) && !defined(QUADRICS)
    rc = armci_rem_strided(op, NULL, proc, src, NULL, dst, NULL,
			   &bytes, 0, NULL, 0, NULL);
    if(rc) armci_die("ARMCI_Value: armci_rem_strided incomplete", FAIL6);
#else

    if(op==PUT) {
      UPDATE_FENCE_STATE(proc, PUT, 1);
#     ifdef LAPI
      SET_COUNTER(ack_cntr[armci_th_idx], 1);
#     endif
#if defined(BGML)
      /* fprintf(stderr,"bytes: %d\n",bytes); */
      /* this call is blocking, so local count is fine */
      BG1S_t req;
      unsigned count=1;
      BGML_Callback_t cb_wait={wait_callback, &count};
      BG1S_Memput(&req, proc, src, 0, dst, bytes, &cb_wait, 1);
      BGML_Wait(&count);
#else

      armci_put(src, dst, bytes, proc);
#endif
    }
    else {
#     ifdef LAPI
      SET_COUNTER(get_cntr[armci_th_idx], 1);
#     endif
#if defined(BGML)
      /* fprintf(stderr,"before memget\n"); */
   BG1S_t req;
   unsigned count=1;
   BGML_Callback_t cb_wait={wait_callback, &count};
   BG1S_Memget(&req, proc, dst, 0, src, bytes, &cb_wait, 1);
   BGML_Wait(&count);

#else
      armci_get(src, dst, bytes, proc);
#endif
    }
    
    /* deal with non-blocking loads and stores */
# if defined(LAPI) || defined(_ELAN_PUTGET_H)
    if(proc != armci_me){
      if(op == GET){
	WAIT_FOR_GETS; /* wait for data arrival */
      }else {
	WAIT_FOR_PUTS; /* data must be copied out*/
      }
    }
#endif
#endif
}

/* non-blocking remote value put/get operation */
static void _armci_nb_rem_value(int op, void *src, void *dst, int proc, 
				int bytes, armci_ihdl_t nb_handle) {  
    int rc=0, pv=0;
    int armci_th_idx = ARMCI_THREAD_IDX;

    if(nb_handle && nb_handle->agg_flag == SET) {
      if(op==PUT) pv = 1;
      (void)armci_agg_save_descriptor(src,dst,bytes,proc,op,pv,nb_handle);
      return;
    }
    else {
      if(op==PUT) UPDATE_FENCE_INFO(proc); 
      
      /*set tag and op in the nb handle*/
      if(nb_handle){
	nb_handle->tag = GET_NEXT_NBTAG();
	nb_handle->op  = op;
	nb_handle->proc= proc;
	nb_handle->bufid=NB_NONE;
      }
      else 
	nb_handle = armci_set_implicit_handle(op, proc);
    }

#if defined(REMOTE_OP) && !defined(QUADRICS)
    rc = armci_rem_strided(op, NULL, proc, src, NULL, dst, NULL,
			   &bytes, 0, NULL, 0, nb_handle);
    if(rc) armci_die("ARMCI_Value: armci_rem_strided incomplete", FAIL6);
#else
    
    if(op==PUT) {
      UPDATE_FENCE_STATE(proc, PUT, 1);
#     ifdef LAPI
      SET_COUNTER(ack_cntr[armci_th_idx], 1);
#     endif
      armci_put(src, dst, bytes, proc);
    }
    else {
#     ifdef LAPI
      SET_COUNTER(get_cntr[armci_th_idx], 1);
#     endif
      armci_get(src, dst, bytes, proc);
    }
    
    /* deal with non-blocking loads and stores */
# if defined(LAPI) || defined(_ELAN_PUTGET_H)
#   ifdef LAPI
    if(!nb_handle)
#   endif
      {
	if(proc != armci_me){
          if(op == GET){
	    WAIT_FOR_GETS; /* wait for data arrival */
          }else {
            WAIT_FOR_PUTS; /* data must be copied out*/
          }
	}
      }
# endif
#endif
}


#define CHK_ERR(dst, proc)       \
    if(dst==NULL) armci_die("PARMCI_PutValue: NULL pointer passed",FAIL);  \
    if(proc<0) armci_die("PARMCI_PutValue: Invalid process rank", proc);

#define CHK_ERR_GET(src, dst, proc, bytes)       \
    if(src==NULL || dst==NULL) armci_die("PARMCI_GetValue: NULL pointer passed",FAIL);  \
    if(proc<0) armci_die("PARMCI_GetValue: Invalid process rank", proc); \
    if(bytes<0) armci_die("PARMCI_GetValue: Invalid size", bytes);

/** 
 * Register-Originated Put.
 */
int PARMCI_PutValueInt(int src, void *dst, int proc) 
{
    CHK_ERR(dst, proc);
    if( SAMECLUSNODE(proc) ) *(int *)dst = src;
    else _armci_rem_value(PUT, &src, dst, proc, sizeof(int));
    return 0;
}

int PARMCI_PutValueLong(long src, void *dst, int proc) 
{
    CHK_ERR(dst, proc);
    if( SAMECLUSNODE(proc) ) *(long *)dst = src;
    else _armci_rem_value(PUT, &src, dst, proc, sizeof(long));
    return 0;
}

int PARMCI_PutValueFloat(float src, void *dst, int proc) 
{
    CHK_ERR(dst, proc);
    if( SAMECLUSNODE(proc) ) *(float *)dst = src;
    else _armci_rem_value(PUT, &src, dst, proc, sizeof(float));
    return 0;
}

int PARMCI_PutValueDouble(double src, void *dst, int proc) 
{
    CHK_ERR(dst, proc);
    if( SAMECLUSNODE(proc) ) *(double *)dst = src;
    else _armci_rem_value(PUT, &src, dst, proc, sizeof(double));
    return 0;
}

/**
 * Non-Blocking register-originated put.
 */
int PARMCI_NbPutValueInt(int src, void *dst, int proc, armci_hdl_t* usr_hdl) 
{
    CHK_ERR(dst, proc);
    if( SAMECLUSNODE(proc) ) *(int *)dst = src;
    else _armci_nb_rem_value(PUT,&src,dst,proc,sizeof(int),(armci_ihdl_t)usr_hdl);
    return 0;
}

int PARMCI_NbPutValueLong(long src, void *dst, int proc, armci_hdl_t* usr_hdl) 
{
    CHK_ERR(dst, proc);
    if( SAMECLUSNODE(proc) ) *(long *)dst = src;
    else _armci_nb_rem_value(PUT,&src,dst,proc,sizeof(long),(armci_ihdl_t)usr_hdl);
    return 0;
}

int PARMCI_NbPutValueFloat(float src, void *dst, int proc, armci_hdl_t* usr_hdl)
{
    CHK_ERR(dst, proc);
    if( SAMECLUSNODE(proc) ) *(float *)dst = src;
    else  _armci_nb_rem_value(PUT,&src,dst,proc,sizeof(float),(armci_ihdl_t)usr_hdl);
    return 0;
}

int PARMCI_NbPutValueDouble(double src, void *dst, int proc, armci_hdl_t* usr_hdl)
{
    CHK_ERR(dst, proc);
    if( SAMECLUSNODE(proc) ) *(double *)dst = src;
    else  _armci_nb_rem_value(PUT,&src,dst,proc,sizeof(double),(armci_ihdl_t)usr_hdl);
     return 0;
 }

#if 1
/** 
 * Register-Originated Get.
 */
int PARMCI_GetValueInt(void *src, int proc) 
{
    int dst;
    if( SAMECLUSNODE(proc) ) return *(int *)src;
    else _armci_rem_value(GET, src, &dst, proc, sizeof(int));
    return dst;
}

long PARMCI_GetValueLong(void *src, int proc) 
{
    long dst;
    if( SAMECLUSNODE(proc) ) return *(long *)src;
    else _armci_rem_value(GET, src, &dst, proc, sizeof(long));
    return dst;
}

float PARMCI_GetValueFloat(void *src, int proc) 
{
    float dst;
    if( SAMECLUSNODE(proc) ) return *(float *)src;
    else _armci_rem_value(GET, src, &dst, proc, sizeof(float));
    return dst;
}

double PARMCI_GetValueDouble(void *src, int proc) 
{
    double dst;
    if( SAMECLUSNODE(proc) ) return *(double *)src;
    else _armci_rem_value(GET, src, &dst, proc, sizeof(double));
    return dst;
}

#endif

#if 0
/**
 * Register-Originated Get.
 */
int PARMCI_GetValue(void *src, void *dst, int proc, int bytes) 
{
    CHK_ERR_GET(src, dst, proc, bytes);
    if( SAMECLUSNODE(proc) ) { armci_copy(src, dst, bytes); }
    else _armci_rem_value(GET, src, dst, proc, bytes);
    return 0;
}

/**
 * Non-Blocking register-originated get.
 */
int PARMCI_NbGetValue(void *src, void *dst, int proc, int bytes, armci_hdl_t* usr_hdl) 
{
    CHK_ERR_GET(src, dst, proc, bytes);
    if( SAMECLUSNODE(proc) ) { armci_copy(src, dst, bytes); }
    else _armci_nb_rem_value(GET, src, dst, proc, bytes, (armci_ihdl_t)usr_hdl);
    return 0;
}
#endif

