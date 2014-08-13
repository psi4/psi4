#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: strided.c,v 1.117.2.6 2007-08-29 17:46:40 manoj Exp $ */
#include "armcip.h"
#include "copy.h"
#include "acc.h"
#include "memlock.h"
#include "armci.h"
#include "iterator.h"
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif

#define ARMCI_OP_2D(op, scale, proc, src, dst, bytes, count, src_stride, dst_stride,lockit) \
  if(op == GET || op ==PUT)						\
    armci_copy_2D(op, proc, src, dst, bytes, count, src_stride,dst_stride); \
  else if(count==1) armci_acc_1D(op, scale, proc, src, dst, bytes,lockit); \
  else									\
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

#if defined(BGML) || defined(ARMCIX)
#define PREPROCESS_STRIDED(tmp_count)
#define POSTPROCESS_STRIDED(tmp_count)
#else
#define PREPROCESS_STRIDED(tmp_count) {					\
    tmp_count=0;							\
    if(stride_levels)							\
      for(;stride_levels;stride_levels--)if(count[stride_levels]>1)break; \
    if(stride_levels&&(count[0]==src_stride_arr[0]&&count[0]==dst_stride_arr[0])){ \
      tmp_count=seg_count[1];						\
      count = seg_count+1;						\
      seg_count[1] = seg_count[0] * seg_count[1];			\
      stride_levels --;							\
      src_stride_arr ++;  dst_stride_arr++ ;				\
    }									\
  }
#define POSTPROCESS_STRIDED(tmp_count) if(tmp_count)seg_count[1]=tmp_count
#endif

#define SERVER_GET 1
#define SERVER_NBGET 2
#define DIRECT_GET 3
#define DIRECT_NBGET 4
#define SERVER_PUT 5
#define SERVER_NBPUT 6
#define DIRECT_PUT 7
#define DIRECT_NBPUT 8


#ifdef ALLOW_PIN
#  define DO_FENCE(__proc,__prot) if(__prot==SERVER_GET);	\
  else if(__prot==SERVER_PUT);					\
  else if(__prot==DIRECT_GET || __prot==DIRECT_NBGET){		\
    if(armci_prot_switch_fence[__proc]==SERVER_PUT)		\
      PARMCI_Fence(__proc);					\
  }								\
  else if(__prot==DIRECT_PUT || __prot==DIRECT_NBPUT){		\
    if(armci_prot_switch_fence[__proc]==SERVER_PUT)		\
      PARMCI_Fence(__proc);					\
  }								\
  else;								\
  armci_prot_switch_fence[__proc]=__prot
#else

#  define DO_FENCE(__proc,__prot)  
        
#endif


#ifndef REGIONS_REQUIRE_MEMHDL 
#  define ARMCI_MEMHDL_T void
#endif

ARMCI_MEMHDL_T *mhloc=NULL,*mhrem=NULL; 

#ifdef REGIONS_REQUIRE_MEMHDL 
int armci_region_both_found_hndl(void *loc, void *rem, int size, int node,
				 ARMCI_MEMHDL_T **loc_memhdl,ARMCI_MEMHDL_T **rem_memhdl);
#  define ARMCI_REGION_BOTH_FOUND(_s,_d,_b,_p)				\
  armci_region_both_found_hndl((_s),(_d),(_b),(_p),&mhloc,&mhrem)
#else
#  define ARMCI_REGION_BOTH_FOUND(_s,_d,_b,_p)	\
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
  armci_client_direct_get((_p),(_s),(_d),(_cou)[0],NULL,0,(void *)mhloc,(void *)mhrem) \

#else

#  define ARMCI_REM_GET(_p,_s,_sst,_d,_dst,_cou,_lev,_hdl)		\
  armci_rem_get((_p),(_s),(_sst),(_d),(_dst),(_cou),(_lev),(_hdl),(void *)mhloc,(void *)mhrem)
#  define ARMCI_NBREM_GET ARMCI_REM_GET
        
#endif

#ifdef ALLOW_PIN
extern int* armci_prot_switch_fence;
extern int armci_prot_switch_preproc;
extern int armci_prot_switch_preop;
#endif

        
int armci_iwork[MAX_STRIDE_LEVEL];

/*\ 2-dimensional array copy
  \*/
static void armci_copy_2D(int op, int proc, void *src_ptr, void *dst_ptr, 
                          int bytes, int count, int src_stride, int dst_stride)
{
#ifdef LAPI
  int armci_th_idx = ARMCI_THREAD_IDX;
#endif
    
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

    if(count==1 && bytes <THRESH1D){
	  
      armci_copy(src_ptr, dst_ptr, bytes); 

    }else {
            
      if(bytes < THRESH){ /* low-latency copy for small data segments */        
#if defined(__crayx1)
	if( !(bytes%sizeof(float)) ) {
	  float *ps=(float*)src_ptr;
	  float *pd=(float*)dst_ptr;
	  long fsstride = src_stride/sizeof(float);
	  long fdstride = dst_stride/sizeof(float);
	  int j;
                
	  for (j = 0;  j < count;  j++){
	    int i;
#pragma _CRI concurrent
	    for(i=0;i<bytes/sizeof(float);i++) pd[i] = ps[i];
	    ps += fsstride;
	    pd += fdstride;
	  }
	} else
#endif
	  {
	    char *ps=(char*)src_ptr;
	    char *pd=(char*)dst_ptr;
	    int j;
		  
	    for (j = 0;  j < count;  j++){
	      int i;
	      for(i=0;i<bytes;i++) pd[i] = ps[i];
	      ps += src_stride;
	      pd += dst_stride;
	    }
	  }
      } else if(bytes %ALIGN_SIZE
		|| dst_stride % ALIGN_SIZE
		|| src_stride % ALIGN_SIZE
#ifdef PTR_ALIGN
		|| (unsigned long)src_ptr%ALIGN_SIZE
		|| (unsigned long)dst_ptr%ALIGN_SIZE
#endif
		){ 

	/* size/address not alligned */
	ByteCopy2D(bytes, count, src_ptr, src_stride, dst_ptr, dst_stride);
                
      }else { /* size aligned -- should be the most efficient copy */
                
	DCopy2D(bytes/ALIGN_SIZE, count,src_ptr, src_stride/ALIGN_SIZE, 
		dst_ptr, dst_stride/ALIGN_SIZE);
      }
    }
        
  } else {
        
    /* data not in local/shared memory-access through global address space*/
        
    if(op==PUT){ 
            
      UPDATE_FENCE_STATE(proc, PUT, COUNT);
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
  switch (op){
  case ARMCI_ACC_INT:
    rows = bytes/sizeof(int);
    break;
  case ARMCI_ACC_LNG:
    rows = bytes/sizeof(long);
    break;
  case ARMCI_ACC_DBL:
    rows = bytes/sizeof(double);
    break;
  case ARMCI_ACC_DCP:
    rows = bytes/(2*sizeof(double));
    break;
  case ARMCI_ACC_CPL:
    rows = bytes/(2*sizeof(float));
    break;
  case ARMCI_ACC_FLT:
    rows = bytes/sizeof(float);
    break;
  default:
    armci_die("ARMCI accumulate: operation not supported",op);
  }
  if(lockit)ARMCI_LOCKMEM(dst, bytes + (char*)dst, proc);
  switch (op){
  case ARMCI_ACC_INT:
    I_ACCUMULATE_1D(scale, dst, src, &rows);
    break;
  case ARMCI_ACC_LNG:
    L_ACCUMULATE_1D(scale, dst, src, &rows);
    break;
  case ARMCI_ACC_DBL:
    D_ACCUMULATE_1D(scale, dst, src, &rows);
    break;
  case ARMCI_ACC_DCP:
    Z_ACCUMULATE_1D(scale, dst, src, &rows);
    break;
  case ARMCI_ACC_CPL:
    C_ACCUMULATE_1D(scale, dst, src, &rows);
    break;
  case ARMCI_ACC_FLT:
    F_ACCUMULATE_1D(scale, dst, src, &rows);
    break;
  default:
    break;
  }
  if(lockit)ARMCI_UNLOCKMEM(proc);
}

/*\ 2-dimensional accumulate
  \*/
  void armci_acc_2D(int op, void* scale, int proc, void *src_ptr, void *dst_ptr,
		    int bytes, int cols, int src_stride, int dst_stride, int lockit)
{
  int   rows, lds, ldd, span;

  /*
    if((long)src_ptr%ALIGN)armci_die("src not aligned",(long)src_ptr);
    if((long)dst_ptr%ALIGN)armci_die("src not aligned",(long)dst_ptr);
  */

  switch (op){
  case ARMCI_ACC_INT:
    rows = bytes/sizeof(int);
    ldd  = dst_stride/sizeof(int);
    lds  = src_stride/sizeof(int);
    break;
  case ARMCI_ACC_LNG:
    rows = bytes/sizeof(long);
    ldd  = dst_stride/sizeof(long);
    lds  = src_stride/sizeof(long);
    break;
  case ARMCI_ACC_DBL:
    rows = bytes/sizeof(double);
    ldd  = dst_stride/sizeof(double);
    lds  = src_stride/sizeof(double);
    break;
  case ARMCI_ACC_DCP:
    rows = bytes/(2*sizeof(double));
    ldd  = dst_stride/(2*sizeof(double));
    lds  = src_stride/(2*sizeof(double));
    break;
  case ARMCI_ACC_CPL:
    rows = bytes/(2*sizeof(float));
    ldd  = dst_stride/(2*sizeof(float));
    lds  = src_stride/(2*sizeof(float));
    break;
  case ARMCI_ACC_FLT:
    rows = bytes/sizeof(float);
    ldd  = dst_stride/sizeof(float);
    lds  = src_stride/sizeof(float);
    break;
  default:
    armci_die("ARMCI accumulate: operation not supported",op);
  }
  if(lockit){ 
    span = cols*dst_stride;
    ARMCI_LOCKMEM(dst_ptr, span + (char*)dst_ptr, proc);
  }
  switch (op){
  case ARMCI_ACC_INT:
    I_ACCUMULATE_2D(scale, &rows, &cols, dst_ptr, &ldd, src_ptr, &lds);
    break;
  case ARMCI_ACC_LNG:
    L_ACCUMULATE_2D(scale, &rows, &cols, dst_ptr, &ldd, src_ptr, &lds);
    break;
  case ARMCI_ACC_DBL:
    D_ACCUMULATE_2D(scale, &rows, &cols, dst_ptr, &ldd, src_ptr, &lds);
    break;
  case ARMCI_ACC_DCP:
    Z_ACCUMULATE_2D(scale, &rows, &cols, dst_ptr, &ldd, src_ptr, &lds);
    break;
  case ARMCI_ACC_CPL:
    C_ACCUMULATE_2D(scale, &rows, &cols, dst_ptr, &ldd, src_ptr, &lds);
    break;
  case ARMCI_ACC_FLT:
    F_ACCUMULATE_2D(scale, &rows, &cols, dst_ptr, &ldd, src_ptr, &lds);
    break;
  default:
    break;
  }
  if(lockit)ARMCI_UNLOCKMEM(proc);
}


/*\ compute range of strided data AND lock it
  \*/
  static void 
  armci_lockmem_patch(void* dst_ptr, int dst_stride_arr[], int count[], int stride_levels, int proc)
{
  long span = count[stride_levels];
  span *= dst_stride_arr[stride_levels-1];

  /* lock region of remote memory */
  ARMCI_LOCKMEM(dst_ptr, span + (char*)dst_ptr, proc);
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
  void *buf_ptr;
  int  rc, i, *buf_stride_arr = armci_iwork;

  dassert(1, !SERVER_CONTEXT);
  buf_ptr = malloc(sizeof(double)*BUFSIZE_DBL);
  dassert(1,buf_ptr);
  armci_lockmem_patch(dst_ptr,dst_stride_arr, count, stride_levels, proc);

  /* setup stride array for internal buffer */
  buf_stride_arr[0]=count[0];
  for(i=0; i< stride_levels; i++) {
    buf_stride_arr[i+1]= buf_stride_arr[i]*count[i+1];
  }

  /* get remote data to local buffer */
  rc = armci_op_strided(GET, scale, proc, dst_ptr, dst_stride_arr, buf_ptr, 
			buf_stride_arr, count, stride_levels, 0,NULL);

  if(rc) { ARMCI_UNLOCKMEM(proc); free(buf_ptr); return(rc); }

  /* call local accumulate with lockit=0 (we locked it already) and proc=me */
  rc = armci_op_strided(optype, scale, armci_me, src_ptr, src_stride_arr, 
			buf_ptr,buf_stride_arr, count, stride_levels,0,NULL);
  if(rc) { ARMCI_UNLOCKMEM(proc); free(buf_ptr); return(rc); }

  /* put data back from the buffer to remote location */
  rc = armci_op_strided(PUT, scale, proc, buf_ptr, buf_stride_arr, dst_ptr, 
			dst_stride_arr, count, stride_levels,0,NULL);

  FENCE_NODE(proc); /* make sure put completes before unlocking */
  ARMCI_UNLOCKMEM(proc);    /* release memory lock */
  free(buf_ptr);
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
  long s2, s3, i,j;
  int unlockit=0;
  int total_of_2D;
  int index[MAX_STRIDE_LEVEL], unit[MAX_STRIDE_LEVEL];

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

#  if defined(LAPI2) || defined(PORTALS) /*|| defined(DOELAN4) && !defined(NB_NONCONT)*/
  /*even 1D armci_nbput has to use different origin counters for 1D */
#   if defined(LAPI2)
  if(!ARMCI_ACC(op) && !SAMECLUSNODE(proc) && (nb_handle || 
					 (!nb_handle && stride_levels>=1 && count[0]<=LONG_PUT_THRESHOLD))) 
#   elif defined(DOELAN4) && !defined(NB_NONCONT)
    /*if(!ARMCI_ACC(op) && !SAMECLUSNODE(proc) && nb_handle && stride_levels<2)*/
    if(!ARMCI_ACC(op) && !SAMECLUSNODE(proc) && stride_levels<2)
#   else
      if(!SAMECLUSNODE(proc))
#   endif
	armci_network_strided(op,scale,proc,src_ptr,src_stride_arr,dst_ptr,
			      dst_stride_arr,count,stride_levels,nb_handle);
      else
#  endif
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
    
  /* deal with non-blocking loads and stores */
#if defined(LAPI) || defined(_ELAN_PUTGET_H) || defined(NB_NONCONT)
#   if defined(LAPI)
  if(!nb_handle)
#   endif
    {
      if(!(SAMECLUSNODE(proc))){
	if(op == GET){
	  WAIT_FOR_GETS; /* wait for data arrival */
	}else { 
	  WAIT_FOR_PUTS; /* data must be copied out*/ 
	}
      }
    }
#endif

  /*    if(proc!=armci_me) INTR_ON;*/

  /*      __asm__ __volatile__ ("sfence":::"memory"); */

  if(unlockit){
#      if defined(ACC_COPY)
    FENCE_NODE(proc); 
#      endif
    ARMCI_UNLOCKMEM(proc);    /* release memory lock */
  }

  return 0;
}

/**Internal puts function. Combines both blocking and non-blocking
 * variants. Any use of implicit handles should be outside of this
 * function. 
 * @param src_ptr  pointer to 1st segment at source
 * @param src_stride_arr array of strides at source 
 * @param dst_ptr   pointer to 1st segment at destination
 * @param dst_stride_arr  array of strides at destination
 * @param seg_count number of segments at each stride levels: count[0]=bytes
 * @param stride_levels number of stride levels 
 * @param proc remote process(or) ID
 * @param nbh non-blocking handle (NULL implies blocking call)
 * @param put_flag Flag to set after the PUT is remote complete (if any)
 * @return 
 */
static int _armci_puts(void *src_ptr,
		       int src_stride_arr[],
		       void* dst_ptr,       
		       int dst_stride_arr[],
		       int seg_count[],     
		       int stride_levels,  
		       int proc,              
		       armci_ihdl_t nbh,
		       armci_flag_t *put_flag) {
  int *count=seg_count, tmp_count=0;
  int rc=0, direct=1;
  
  if(src_ptr == NULL || dst_ptr == NULL) return FAIL;
  if(count[0]<0)return FAIL3;
  if(stride_levels <0 || stride_levels > MAX_STRIDE_LEVEL) return FAIL4;
  if(proc<0)return FAIL5;

#ifdef __crayx1
  if(!stride_levels) {
    memcpy(dst_ptr, src_ptr,count[0]);
    return 0;
  }
#endif

  PREPROCESS_STRIDED(tmp_count);
#  if (!defined(QUADRICS) || defined(PACKPUT))
  direct=SAMECLUSNODE(proc);
#  endif /*(!QUADRICS||!PACKPUT)&&!PORTALS*/

  if(put_flag) dassert(1,nbh==NULL);

  if(!nbh) {
    ORDER(PUT,proc);/* ensure ordering */
  }
  else {/* aggregate put */
    if(nbh->agg_flag == SET) {
      if(!direct){ 
	rc= armci_agg_save_strided_descriptor(src_ptr, src_stride_arr, 
					      dst_ptr, dst_stride_arr, 
					      count, stride_levels, proc, 
					      PUT, nbh);
	POSTPROCESS_STRIDED(tmp_count);
	return(rc);
      }
    } else {
      /*ORDER(PUT,proc);  ensure ordering */
      UPDATE_FENCE_INFO(proc);
	
      /*set tag and op in the nb handle*/
      nbh->tag = GET_NEXT_NBTAG();
      nbh->op  = PUT;
      nbh->proc= proc;
      nbh->bufid=NB_NONE;
    }
  }
   
#ifdef BGML
  if(nbh) {
    nbh->count = 1;
    BGML_Callback_t cb_wait={wait_callback, &nbh->count};
    BG1S_MemputS (&nbh->cmpl_info, proc,
		  src_ptr, src_stride_arr,
		  dst_ptr, dst_stride_arr,
		  seg_count, stride_levels,
		  0, &cb_wait, 1);
  }
  else if(!stride_levels) {
    unsigned temp_count=1;
    BGML_Callback_t cb_wait={wait_callback, &temp_count};
    BG1S_t request;
    BGML_CriticalSection_enter();
    BG1S_Memput(&request, proc, src_ptr, 0, dst_ptr, count[0], &cb_wait, 1);
    /*BGML_Wait(&count);*/
    while (temp_count) BGML_Messager_advance();
    BGML_CriticalSection_exit();
  }
  else {
    armci_hdl_t nb_handle;
    ARMCI_INIT_HANDLE(&nb_handle);
    PARMCI_NbPutS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count,
		 stride_levels, proc, &nb_handle);
    PARMCI_Wait(&nb_handle);
  }
  if(put_flag) { /*=>!nbh*/
    PARMCI_Fence(proc);
    PARMCI_Put(&put_flag->val,put_flag->ptr,sizeof(int),proc);
  }
#elif ARMCIX
  if(nbh) 
    ARMCIX_NbPutS (src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, nbh);
  else if(!stride_levels) {
    ARMCIX_Put(src_ptr, dst_ptr, count[0], proc);
  }
  else {
    ARMCIX_PutS (src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc);
  }
  if(put_flag) { /*=>!nbh*/
    PARMCI_Fence(proc);
    PARMCI_Put(&put_flag->val,put_flag->ptr,sizeof(int),proc);
  }
#else /*BGML*/

  /* use direct protocol for remote access when performance is better */
#  if defined(LAPI) || defined(DOELAN4)
  if(!direct) {
    switch(stride_levels) {
    case 0:
#      ifndef LAPI_RDMA
       direct =1;
#      endif
       break;
    case 1:  if((count[1]<PACKPUT)||count[0]>LONG_PUT_THRESHOLD) direct =1; break;
    default: if(count[0]> LONG_PUT_THRESHOLD )direct=1; break;
    }
  }
#  endif /*LAPI||DOELAN4*/
#  ifdef PORTALS
     if(stride_levels) direct=1;
#  endif
  
#  if !defined(LAPI2) || defined(LAPI_RDMA)
  if(!direct){
#    ifdef ALLOW_PIN /*if we can pin, we do*/
    if(!stride_levels && 
       ARMCI_REGION_BOTH_FOUND(src_ptr,dst_ptr,count[0],armci_clus_id(proc))){
      if(nbh) { 
	DO_FENCE(proc,DIRECT_NBPUT); 
	armci_client_direct_send(proc, src_ptr, dst_ptr, count[0],
				 (void **)(&nbh->cmpl_info),
				 nbh->tag,mhloc,mhrem);
      }
      else { 
	DO_FENCE(proc,DIRECT_PUT); 
	armci_client_direct_send(proc,src_ptr,dst_ptr,count[0],NULL,0,mhloc,mhrem);	
      }
      POSTPROCESS_STRIDED(tmp_count);
      if(put_flag) {
	PARMCI_Fence(proc);
	PARMCI_Put(&put_flag->val,put_flag->ptr,sizeof(int),proc);
      }
      return 0;
    }
#      if  0 && defined(VAPI)
#        if !defined(PEND_BUFS)
    if(stride_levels==1 && count[0]>VAPI_SGPUT_MIN_COLUMN &&
       ARMCI_REGION_BOTH_FOUND(src_ptr,dst_ptr,count[0],armci_clus_id(proc))){
      if(nbh) { DO_FENCE(proc,DIRECT_NBPUT); }
      else    { DO_FENCE(proc,DIRECT_PUT); }
      /* 	   printf("%d:Calling two phase send\n", armci_me); */
      armci_two_phase_send(proc, src_ptr, src_stride_arr, dst_ptr,
			   dst_stride_arr,count,stride_levels,NULL,nbh,mhloc);
      if(put_flag) {
	PARMCI_Fence(proc);
	PARMCI_Put(&put_flag->val,put_flag->ptr,sizeof(int),proc);
      }
      return 0;  
    }
#        else /*!PEND_BUFS*/
    {
      int i, off;
      for(i=0; i<stride_levels; i++) {
	if(i==0) assert(src_stride_arr[0]>0);
	if(i!=0) assert(src_stride_arr[i]>=src_stride_arr[i-1]*count[i]);
      }
      off = (stride_levels>0)?
	src_stride_arr[stride_levels-1]*count[stride_levels]
	: count[0];
	
      mhloc=mhrem=NULL;
      if(ARMCI_REGION_BOTH_FOUND(src_ptr,dst_ptr,off,armci_clus_id(proc))) {
	assert(mhloc != NULL);
	assert(mhrem != NULL);
	if(nbh) {
	  DO_FENCE(proc, DIRECT_NBPUT);
	  armci_client_direct_rdma_strided(PUT,proc,src_ptr,src_stride_arr,
					   dst_ptr,dst_stride_arr,
					   count, stride_levels,
					   (void**)&nbh->cmpl_info,nbh->tag,
					   mhloc,mhrem);
	}
	else {
	  DO_FENCE(proc, DIRECT_PUT);
	  armci_client_direct_rdma_strided(PUT,proc,src_ptr,src_stride_arr,
					   dst_ptr,dst_stride_arr,
					   count, stride_levels,NULL,0,
					   mhloc,mhrem);
	}
	if(put_flag) {
	  PARMCI_Fence(proc);
	  PARMCI_Put(&put_flag->val,put_flag->ptr,sizeof(int),proc);
	}
	return 0;
      }
    }
#        endif /*!PEND_BUFS*/
#      endif /*VAPI*/
#    endif /*ALLOW_PIN*/
  }
#endif /* !LAPI2||LAPI_RDMA */
  
#  ifndef LAPI2
  if(!direct){
    if(nbh) { DO_FENCE(proc,SERVER_PUT); }
    else    { DO_FENCE(proc,SERVER_NBPUT); }

#    if defined(DATA_SERVER) && defined(SOCKETS) && defined(USE_SOCKET_VECTOR_API)
    if(count[0]> LONG_PUT_THRESHOLD && stride_levels>0){
      ext_header_t h, *hdr;
      h.exthdr = put_flag;
      h.len = sizeof(armci_flag_t);
      hdr = put_flag?&h:NULL;
      rc = armci_rem_strided(PUT, NULL, proc, src_ptr, src_stride_arr,
			     dst_ptr, dst_stride_arr, count, stride_levels,
			     hdr,1, nbh);
    }
    else
#    endif /*DATA_SERVER && SOCKETS && USE_SOCKET_VECTOR_API*/
      {
	ext_header_t h,*hdr;
	h.exthdr = put_flag;
	h.len = sizeof(armci_flag_t);
	hdr = put_flag?&h:NULL;
	if(nbh) {
	  nbh->tag =0; /* packed request is completed locally */ 
	  CLEAR_HNDL_FIELD(nbh->cmpl_info);
	}
	rc = armci_pack_strided(PUT,NULL,proc,src_ptr,src_stride_arr,
				dst_ptr,dst_stride_arr, 
				count, stride_levels,hdr,-1,-1,-1,NULL);
      }
  }
  else
#  endif /*!LAPI*/
    {
      if(!nbh && stride_levels == 0) {
	armci_copy_2D(PUT, proc, src_ptr, dst_ptr, count[0], 1, count[0],
		      count[0]);
#  if defined(LAPI) || defined(_ELAN_PUTGET_H)
	if(proc != armci_me) { WAIT_FOR_PUTS; }
#  endif /*LAPI||_ELAN_PUTGET_H*/
      }
      else {
	rc = armci_op_strided( PUT, NULL, proc, src_ptr, src_stride_arr, 
			       dst_ptr, dst_stride_arr,count,stride_levels, 
			       0,nbh);
      }
      if(put_flag) { /*=>!nbh*/
	PARMCI_Fence(proc);
	PARMCI_Put(&put_flag->val,put_flag->ptr,sizeof(int),proc);
      }
    }
#endif /*BGML*/
  POSTPROCESS_STRIDED(tmp_count);
  if(rc) return FAIL6;
  else return 0;
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
#if 1
  return _armci_puts(src_ptr,src_stride_arr,dst_ptr,dst_stride_arr,
		     seg_count,stride_levels,proc,NULL,NULL);
#else
  armci_hdl_t nbh;
  ARMCI_INIT_HANDLE(&nbh);
  PARMCI_NbPutS(src_ptr,src_stride_arr,dst_ptr,dst_stride_arr,seg_count,stride_levels,proc,&nbh);
  PARMCI_Wait(&nbh);
  return 0;
#endif
}


int PARMCI_PutS_flag_dir(void *src_ptr,   int src_stride_arr[],
			void* dst_ptr,   int dst_stride_arr[],
			int seg_count[], int stride_levels,
			int *flag, int val, int proc) {
  return PARMCI_PutS_flag(src_ptr, src_stride_arr,dst_ptr,dst_stride_arr,
			 seg_count, stride_levels, flag, val, proc);
}

int PARMCI_PutS_flag(void *src_ptr,     int src_stride_arr[],
		    void* dst_ptr,     int dst_stride_arr[],
		    int seg_count[],   int stride_levels,   
		    int *flag,  int val, int proc) {
  armci_flag_t put_flag;
  put_flag.val=val;
  put_flag.ptr=flag;
  return _armci_puts(src_ptr,src_stride_arr,dst_ptr,dst_stride_arr,seg_count,stride_levels,proc,NULL,&put_flag);
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
  armci_hdl_t nbh;
  
  ORDER(GET,proc);
  ARMCI_INIT_HANDLE(&nbh);
  PARMCI_NbGetS(src_ptr,src_stride_arr,dst_ptr,dst_stride_arr,seg_count,stride_levels,proc,&nbh);
  PARMCI_Wait(&nbh);
  return 0;
}

/**Internal strided accumulate. Implicit handles should be used
 * outsise this function.
 * @param optype            operation 
 * @param scale             scale factor x += scale*y 
 * @param src_ptr           pointer to 1st segment at source 
 * @param src_stride_arr[]  array of strides at source 
 * @param dst_ptr           1st segment at destination
 * @param dst_stride_arr[]  array of strides at destination 
 * @param seg_count[]       number of segments at each stride
 *                           levels: count[0]=bytes 
 * @param stride_levels     number of stride levels 
 * @param proc              remote process(or) ID 
 * @param nbh               armci non-blocking call handle 
 * @return
 */
static int _armci_accs( int  optype,    void *scale,         
			void *src_ptr,	int src_stride_arr[],
			void* dst_ptr,	int dst_stride_arr[],
			int seg_count[],int stride_levels,   
			int proc,       armci_ihdl_t nbh) {
  int rc, direct=1;
  int *count=seg_count, tmp_count=0;

  if(src_ptr == NULL || dst_ptr == NULL) return FAIL;
  if(src_stride_arr == NULL || dst_stride_arr ==NULL) return FAIL2;
  if(count[0]<0)return FAIL3;
  if(stride_levels <0 || stride_levels > MAX_STRIDE_LEVEL) return FAIL4;
  if(proc<0)return FAIL5;

  if(!nbh) { ORDER(optype,proc); }
  else { 
    UPDATE_FENCE_INFO(proc); 
    nbh->tag = GET_NEXT_NBTAG();
    nbh->op  = optype;
    nbh->proc= proc;
    nbh->bufid=NB_NONE;
  }

  PREPROCESS_STRIDED(tmp_count);
#ifdef BGML
  armci_ihdl_t inbh;
  armci_hdl_t tmp_hdl;
  if(nbh) inbh = nbh;
  else {
    ARMCI_INIT_HANDLE(&tmp_hdl);
    inbh = (armci_ihdl_t)&tmp_hdl;
  }
  inbh->count=1;
  BGML_Callback_t cb_wait={wait_callback, &inbh->count};
    
  BGML_Op oper1=BGML_PROD;
  BGML_Op oper2=BGML_SUM;
  BGML_Dt dt;
  switch(optype) {
  case ARMCI_ACC_INT:
  case ARMCI_ACC_LNG:
    dt=BGML_SIGNED_INT;
    break;
#if 0
  case ARMCI_ACC_LNG:
    dt=BGML_SIGNED_LONG;
    break;
#endif
  case ARMCI_ACC_DBL:
    dt=BGML_DOUBLE;
    break;
  case ARMCI_ACC_CPL:
    dt=BGML_SINGLE_COMPLEX;
    break;
  case ARMCI_ACC_DCP:
    dt=BGML_DOUBLE_COMPLEX;
    break;
  case ARMCI_ACC_FLT:
    dt=BGML_FLOAT;
    break;
  default:
    assert(0);
  }
    
  BG1S_AccumulateS (&inbh->cmpl_info, proc,
		    src_ptr, src_stride_arr,
		    dst_ptr, dst_stride_arr,
		    seg_count, stride_levels,
		    scale, 0,
		    dt, oper1, oper2,
		    &cb_wait, 1);

  if(!nbh) PARMCI_Wait(&tmp_hdl);
#elif ARMCIX
  if(!nbh)
    ARMCIX_AccS (optype, scale, src_ptr, src_stride_arr, dst_ptr,
		 dst_stride_arr, count, stride_levels, proc);
  else
    ARMCIX_NbAccS (optype, scale, src_ptr, src_stride_arr, dst_ptr,
		   dst_stride_arr, count, stride_levels, proc, nbh);
#else 

  direct=SAMECLUSNODE(proc);

#   if defined(ACC_COPY) && !defined(ACC_SMP)
  if(armci_me != proc) direct=0;
#   endif /*ACC_COPY && !ACC_SMP*/
       
  if(direct)
    rc = armci_op_strided(optype,scale, proc, src_ptr, src_stride_arr,dst_ptr,
			  dst_stride_arr, count, stride_levels,1,NULL);
  else{
    if(nbh) { DO_FENCE(proc,SERVER_NBPUT); }
    else { DO_FENCE(proc,SERVER_PUT); }
    rc = armci_pack_strided(optype,scale,proc,src_ptr, src_stride_arr,dst_ptr,
			    dst_stride_arr,count,stride_levels,NULL,-1,-1,-1,nbh);
  }
#endif /*BGML*/
  POSTPROCESS_STRIDED(tmp_count);
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
                ) {
  return _armci_accs(optype,scale,src_ptr,src_stride_arr,dst_ptr,
		     dst_stride_arr,seg_count,stride_levels,proc,NULL);
}



int PARMCI_Put(void *src, void* dst, int bytes, int proc) {
  int rc=0;
  rc = PARMCI_PutS(src, NULL, dst, NULL, &bytes, 0, proc);
  return rc;
}

int PARMCI_Acc(int optype, void *scale, void *src, void* dst, int bytes, int proc) {
  int rc=0;
  rc = PARMCI_AccS(optype, scale, src, NULL, dst, NULL, &bytes, 0, proc);
  return rc;
}


int PARMCI_Put_flag(void *src, void* dst,int bytes,int *f,int v,int proc) {
  return  PARMCI_PutS_flag(src, NULL, dst, NULL, &bytes, 0, f, v, proc);
}

int PARMCI_Get(void *src, void* dst, int bytes, int proc) {
  int rc=0;
  
#ifdef __crayx1
  memcpy(dst,src,bytes);   
#else
  rc = PARMCI_GetS(src, NULL, dst, NULL, &bytes, 0, proc);
#endif
  
  dassert(1,rc==0);
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
			  int count[], char *buf) {
  const int seg_size = count[0];
  int off=0;
  stride_info_t sinfo;
  armci_stride_info_init(&sinfo, ptr,stride_levels,stride_arr,count);
  while(armci_stride_info_has_more(&sinfo)) {
    char *sptr = armci_stride_info_seg_ptr(&sinfo);
    armci_copy(sptr,&buf[off],seg_size);
    off += seg_size;
    armci_stride_info_next(&sinfo);
  }
  armci_stride_info_destroy(&sinfo);
}


void armci_write_strided2(void *ptr, int stride_levels, int stride_arr[],
                          int count[], char *buf)
{                  
  int i, j;
  int total;   /* number of 2 dim block */
  int index[MAX_STRIDE_LEVEL], unit[MAX_STRIDE_LEVEL];
    
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
      DCOPY21(&rows, count+1, ptr, &ld, (void*)buf, &idx);
      break;
    case 2: 
#if 0
      for(i=0; i< count[2]; i++){ 
	DCOPY21(&rows, count+1, ptr, &ld, buf, &idx);
	ptr = ((char*)ptr)+stride_arr[1];
	buf = (char*) ((double*)buf + idx);
      }
#endif
      ldd = stride_arr[1]/stride_arr[0];
      DCOPY31(&rows, count+1, count+2, ptr, &ld, &ldd, (void*)buf, &idx);

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
	DCOPY21(&rows, count+1, (void*)src, &ld, (void*)buf, &idx); 
	buf = (char*) ((double*)buf + idx);
      }
    } /*switch */
  } /*else */
}


void armci_read_strided1(void *ptr, int stride_levels, int stride_arr[],
			 int count[], char *buf) {
  const int seg_size = count[0];
  int off=0;
  stride_info_t sinfo;
  armci_stride_info_init(&sinfo,ptr,stride_levels,stride_arr,count);
  while(armci_stride_info_has_more(&sinfo)) {
    char *dptr = armci_stride_info_seg_ptr(&sinfo);
    armci_copy(&buf[off],dptr,seg_size);
    off += seg_size;
    armci_stride_info_next(&sinfo);
  }
  armci_stride_info_destroy(&sinfo);
}


void armci_read_strided2(void *ptr, int stride_levels, int stride_arr[],
                         int count[], char *buf)
{                  
  int i, j;
  int total;   /* number of 2 dim block */
  int index[MAX_STRIDE_LEVEL], unit[MAX_STRIDE_LEVEL];
   
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
      DCOPY12(&rows, count+1, ptr, &ld, (void*)buf, &idx);
      break;
    case 2:
#if 0
      for(i=0; i< count[2]; i++){
	DCOPY12(&rows, count+1, ptr, &ld, buf, &idx);
	ptr = ((char*)ptr)+stride_arr[1];
	buf = (char*) ((double*)buf + idx);
      }
#endif
      ldd = stride_arr[1]/stride_arr[0];   
      DCOPY13(&rows, count+1, count+2, ptr, &ld, &ldd, (void*)buf, &idx);
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
	DCOPY12(&rows, count+1, (void*)src, &ld, (void*)buf, &idx);
	buf = (char*) ((double*)buf + idx);
      }
    } /*switch */
  } /*else */
}

/**Read data from buffer into the locations pointed to by the stride
 * iterator. The reading happens incrementally. The stride iterator is
 * traversed to copy as much data as possible in the buffer. When all
 * the data in buf is consumed the function returns with the number of
 * bytes consumed from the buffer.
 * @param sinfo Stride iterator
 * @param buf IN Pointer to data to be read into user memory
 * @param bytes IN #bytes available in buf for reading
 * @param seg_off INOUT Bytes of the current segment written in the
 * last call (on partial segment write). On return, this parameter
 * contains the bytes of the last segment written if it was partial. 
 * @return #bytes read from buf into user memory.
 */
int armci_read_strided_inc(stride_info_t *sinfo, const char *buf,int bytes, int *seg_off) {
  int off=0;
  const int seg_size = armci_stride_info_seg_size(sinfo);

  dassert1(1,bytes>0,bytes);
  off=0;
  if(*seg_off) {
    char *sptr = (char*) &buf[off];
    char *dptr = ((char*)armci_stride_info_seg_ptr(sinfo))+*seg_off;
    int size = ARMCI_MIN(seg_size-*seg_off,bytes);
    /*     printf("%d:%s(): seg_size=%d,seg_off=%d,bytes=%d\n",armci_me,FUNCTION_NAME,seg_size,*seg_off,bytes); */
    dassert(1,armci_stride_info_has_more(sinfo));
    armci_copy(sptr,dptr,size);
    off += size;
    if(*seg_off+size == seg_size) {
      armci_stride_info_next(sinfo);
    }
  }
  while(bytes>off) {
    int size = ARMCI_MIN(seg_size, bytes-off);
    dassert(1,armci_stride_info_has_more(sinfo));
    armci_copy(&buf[off],armci_stride_info_seg_ptr(sinfo),size);
    if(size==seg_size) armci_stride_info_next(sinfo);
    off += size;
  }
  dassertp(1,off==bytes,("%d:off=%d bytes=%d",armci_me,off,bytes));
  *seg_off = (bytes + *seg_off) % seg_size;
  return bytes;
}

#define INIT_NB_HANDLE(nb,o,p) if(nb){				\
    (nb)->tag = 0;						\
    (nb)->op  = (o); (nb)->proc= (p);				\
    (nb)->bufid=NB_NONE;}					\
  else { (nb)=(armci_ihdl_t)armci_set_implicit_handle(o, p); (nb)->tag=0; }



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
  if(!usr_hdl) usr_hdl = armci_set_implicit_handle(PUT, proc);
  return _armci_puts(src_ptr, src_stride_arr,dst_ptr,dst_stride_arr,
		     seg_count,stride_levels,proc,(armci_ihdl_t)usr_hdl,NULL);
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

  if(src_ptr == NULL || dst_ptr == NULL) return FAIL;
  if(seg_count[0]<0)return FAIL3;
  if(stride_levels <0 || stride_levels > MAX_STRIDE_LEVEL) return FAIL4;
  if(proc<0)return FAIL5;

#ifdef BGML
  armci_ihdl_t nbh;
  set_nbhandle(&nbh, usr_hdl, PUT, proc);
  nbh->count=1;
  BGML_Callback_t cb_wait={wait_callback, &nbh->count};

  BG1S_MemgetS (&nbh->cmpl_info, proc,
		src_ptr, src_stride_arr,
		dst_ptr, dst_stride_arr,
		seg_count, stride_levels,
		0, &cb_wait, 1);
#else

#if !defined(QUADRICS)
  direct=SAMECLUSNODE(proc);
#endif
  PREPROCESS_STRIDED(tmp_count);

  /* aggregate get */
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
    UPDATE_FENCE_INFO(proc);
    /*set tag and op in the nb handle*/
    if(nb_handle){
      nb_handle->tag = GET_NEXT_NBTAG();
      nb_handle->op  = GET;
      nb_handle->proc= proc;
      nb_handle->bufid=NB_NONE;
    }
    else
      nb_handle = (armci_ihdl_t)armci_set_implicit_handle(GET, proc);
  }

#ifdef LAPI_RDMA
  if(stride_levels == 0 || count[0] > LONG_GET_THRESHOLD) 
      direct=0;
#endif

#ifdef PORTALS
  if(stride_levels) 
      direct=1;
#endif
  
#if !defined(LAPI2) || defined(LAPI_RDMA)
  if(!direct){
#     ifdef ALLOW_PIN
#if defined(VAPI)
    extern int armci_max_num_sg_ent;
#endif
    if(!stride_levels && 
       ARMCI_REGION_BOTH_FOUND(dst_ptr,src_ptr,count[0],armci_clus_id(proc))){
      DO_FENCE(proc,DIRECT_NBGET);
      ARMCI_NBREM_GET(proc, src_ptr,NULL,dst_ptr,NULL,count, 0, nb_handle);
      POSTPROCESS_STRIDED(tmp_count);
      return 0;
    }
#     endif
  }
#endif /*!LAPI||LAPI_RDMA */
  
#ifndef LAPI2
  if(!direct){
    DO_FENCE(proc,SERVER_NBGET);
#if defined(DATA_SERVER) && (defined(SOCKETS) || defined(CLIENT_BUF_BYPASS) )
    /* for larger strided or 1D reqests buffering can be avoided to send data
     * we can try to bypass the packetization step and send request directly
     */
    if(CAN_REQUEST_DIRECTLY && ((count[0]> LONG_GET_THRESHOLD) ||
				(stride_levels && count[0]>LONG_GET_THRESHOLD_STRIDED) ) ) {

      int nobuf =1; /* tells the sending routine not to buffer */
      rc = armci_rem_strided(GET, NULL, proc,src_ptr,src_stride_arr,dst_ptr,
			     dst_stride_arr, count, stride_levels,
			     (ext_header_t*)0,nobuf,nb_handle);
      if(rc) goto DefaultPath; /* attempt to avoid buffering failed */ 

    }else
    DefaultPath: /* standard buffered path */
#endif
#ifdef ARMCIX
      rc = ARMCIX_NbGetS (src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, nb_handle);
#else
    rc = armci_pack_strided(GET, NULL, proc, src_ptr, src_stride_arr,
			    dst_ptr,dst_stride_arr,count,stride_levels,
			    NULL,-1,-1,-1,nb_handle);
#endif
  }else
#else
    /* avoid LAPI_GetV */
    if(stride_levels==1 && count[0]>320 && !direct) 
      ARMCI_REM_GET(proc,src_ptr,src_stride_arr,dst_ptr,
		    dst_stride_arr, count, stride_levels, nb_handle);
    else
#endif
      rc = armci_op_strided(GET, NULL, proc, src_ptr, src_stride_arr, dst_ptr,
			    dst_stride_arr,count, stride_levels,0,nb_handle);

  POSTPROCESS_STRIDED(tmp_count);
#endif /*bgml*/

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
		  ) {
  if(!usr_hdl) usr_hdl = armci_set_implicit_handle(optype, proc);
  return _armci_accs(optype,scale,src_ptr,src_stride_arr,dst_ptr,
		     dst_stride_arr,seg_count,stride_levels,proc,(armci_ihdl_t)usr_hdl);
}


#if !defined(ACC_COPY)&&!defined(CRAY_YMP)&&!defined(CYGNUS)&&!defined(CYGWIN) &&!defined(BGML)&&!defined(DCMF)
#   define REMOTE_OP
#endif

void set_nbhandle(armci_ihdl_t *nbh, armci_hdl_t *nb_handle, int op,
		  int proc) {
  if(nb_handle) {
    *nbh=(armci_ihdl_t)nb_handle;
  }
  else  {
    *nbh=(armci_ihdl_t)armci_set_implicit_handle(op, proc);
  }
}


int PARMCI_NbPut(void *src, void* dst, int bytes, int proc,armci_hdl_t* uhandle) {
  int rc;
  rc = PARMCI_NbPutS(src,NULL,dst,NULL,&bytes,0,proc,uhandle);
  return(rc);
}


int PARMCI_NbGet(void *src, void* dst, int bytes, int proc,armci_hdl_t* uhandle) {
  int rc;
  rc=PARMCI_NbGetS(src,NULL,dst,NULL,&bytes,0,proc,uhandle);
  return(rc);
}


static void _armci_op_value(int op, void *src, void *dst, int proc, 
			    int bytes, armci_hdl_t *usr_hdl) {
  int rc=0,pv=0;
#ifdef LAPI
  int armci_th_idx = ARMCI_THREAD_IDX;
#endif
  armci_ihdl_t nbh = (armci_ihdl_t)usr_hdl;

  if(!nbh) {
    ORDER(op,proc); /* ensure ordering */
  }else {
    if(nbh->agg_flag == SET) {
      if(op==PUT) pv = 1;
      (void)armci_agg_save_descriptor(src,dst,bytes,proc,op,pv,nbh);
      return;
    }
    else {
      if(op==PUT) UPDATE_FENCE_INFO(proc); 
	
      /*set tag and op in the nb handle*/
      nbh->tag = GET_NEXT_NBTAG();
      nbh->op  = op;
      nbh->proc= proc;
      nbh->bufid=NB_NONE;
    }
  }
#if defined(REMOTE_OP) && !defined(QUADRICS)
  rc = armci_rem_strided(op, NULL, proc, src, NULL, dst, NULL,
			 &bytes, 0, NULL, 0, nbh);
  if(rc) armci_die("ARMCI_Value: armci_rem_strided incomplete", FAIL6);
#else
  if(op==PUT) {
    UPDATE_FENCE_STATE(proc, PUT, 1);
#  ifdef LAPI
    SET_COUNTER(ack_cntr[armci_th_idx], 1);
#  endif
#  if defined(BGML) || defined(ARMCIX)
    if(usr_hdl) PARMCI_NbPut(src,dst,bytes,proc,usr_hdl);
    else PARMCI_Put(src,dst,bytes,proc);
#  else
    armci_put(src, dst, bytes, proc);
#  endif
  }
  else {
#  ifdef LAPI
    SET_COUNTER(get_cntr[armci_th_idx], 1);
#  endif
#  if defined(BGML) || defined(ARMCIX)
    if(usr_hdl) PARMCI_NbGet(src,dst,bytes,proc,usr_hdl);
    else PARMCI_Get(src,dst,bytes,proc);
#  else
    armci_get(src, dst, bytes, proc);
#  endif
  }
    
  /* deal with non-blocking loads and stores */
#  if defined(LAPI) || defined(_ELAN_PUTGET_H)
#    ifdef LAPI
  if(!nbh)
#    endif
    {
      if(proc != armci_me){
	if(op == GET){
	  WAIT_FOR_GETS; /* wait for data arrival */
	}else {
	  WAIT_FOR_PUTS; /* data must be copied out*/
	}
      }
    }
#  endif
#endif
}

static void _armci_rem_value(int op, void *src, void *dst, int proc, 
			     int bytes) {  
  _armci_op_value(op,src,dst,proc,bytes,NULL);
}

/* non-blocking remote value put/get operation */
static void _armci_nb_rem_value(int op, void *src, void *dst, int proc, 
				int bytes, armci_hdl_t *usr_hdl) {  
  if(!usr_hdl) usr_hdl = (armci_hdl_t*)armci_set_implicit_handle(op,proc);
  _armci_op_value(op,src,dst,proc,bytes,usr_hdl);
}


#define CHK_ERR(dst, proc)						\
  if(dst==NULL) armci_die("ARMCI_PutValue: NULL pointer passed",FAIL);  \
  if(proc<0) armci_die("ARMCI_PutValue: Invalid process rank", proc);

#define CHK_ERR_GET(src, dst, proc, bytes)				\
  if(src==NULL || dst==NULL) armci_die("ARMCI_GetValue: NULL pointer passed",FAIL); \
  if(proc<0) armci_die("ARMCI_GetValue: Invalid process rank", proc);	\
  if(bytes<0) armci_die("ARMCI_GetValue: Invalid size", bytes);

/** 
 * Register-Originated Put.
 */
int PARMCI_PutValueInt(int src, void *dst, int proc) {
  CHK_ERR(dst, proc);
  if( SAMECLUSNODE(proc) ) *(int *)dst = src;
  else _armci_rem_value(PUT, &src, dst, proc, sizeof(int));
  return 0;
}

int PARMCI_PutValueLong(long src, void *dst, int proc) {
  CHK_ERR(dst, proc);
  if( SAMECLUSNODE(proc) ) *(long *)dst = src;
  else _armci_rem_value(PUT, &src, dst, proc, sizeof(long));
  return 0;
}

int PARMCI_PutValueFloat(float src, void *dst, int proc) {
  CHK_ERR(dst, proc);
  if( SAMECLUSNODE(proc) ) *(float *)dst = src;
  else _armci_rem_value(PUT, &src, dst, proc, sizeof(float));
  return 0;
}

int PARMCI_PutValueDouble(double src, void *dst, int proc) {
  CHK_ERR(dst, proc);
  if( SAMECLUSNODE(proc) ) *(double *)dst = src;
  else _armci_rem_value(PUT, &src, dst, proc, sizeof(double));
  return 0;
}

/**
 * Non-Blocking register-originated put.
 */
int PARMCI_NbPutValueInt(int src, void *dst, int proc, armci_hdl_t* usr_hdl) {
  CHK_ERR(dst, proc);
  if( SAMECLUSNODE(proc) ) *(int *)dst = src;
  else _armci_nb_rem_value(PUT,&src,dst,proc,sizeof(int),usr_hdl);
  return 0;
}

int PARMCI_NbPutValueLong(long src, void *dst, int proc, armci_hdl_t* usr_hdl) {
  CHK_ERR(dst, proc);
  if( SAMECLUSNODE(proc) ) *(long *)dst = src;
  else _armci_nb_rem_value(PUT,&src,dst,proc,sizeof(long),usr_hdl);
  return 0;
}

int PARMCI_NbPutValueFloat(float src, void *dst, int proc, armci_hdl_t* usr_hdl) {
  CHK_ERR(dst, proc);
  if( SAMECLUSNODE(proc) ) *(float *)dst = src;
  else  _armci_nb_rem_value(PUT,&src,dst,proc,sizeof(float),usr_hdl);
  return 0;
}

int PARMCI_NbPutValueDouble(double src, void *dst, int proc, armci_hdl_t* usr_hdl) {
  CHK_ERR(dst, proc);
  if( SAMECLUSNODE(proc) ) *(double *)dst = src;
  else  _armci_nb_rem_value(PUT,&src,dst,proc,sizeof(double),usr_hdl);
  return 0;
}

/** 
 * Register-Originated Get.
 */
int PARMCI_GetValueInt(void *src, int proc) 
{
    int dst;

    if (SAMECLUSNODE(proc)) 
        return *(int *)src;
    else 
        _armci_rem_value(GET, src, &dst, proc, sizeof(int));

    return dst;
}

long PARMCI_GetValueLong(void *src, int proc) 
{
    long dst;

    if (SAMECLUSNODE(proc)) 
        return *(long *)src;
    else 
        _armci_rem_value(GET, src, &dst, proc, sizeof(long));

    return dst;
}

float PARMCI_GetValueFloat(void *src, int proc) 
{
    float dst;

    if (SAMECLUSNODE(proc)) 
        return *(float *)src;
    else 
        _armci_rem_value(GET, src, &dst, proc, sizeof(float));

    return dst;
}

double PARMCI_GetValueDouble(void *src, int proc) 
{
  double dst;
 
  if(SAMECLUSNODE(proc)) 
      return *(double *)src;
  else 
      _armci_rem_value(GET, src, &dst, proc, sizeof(double));
  
  return dst;
}



