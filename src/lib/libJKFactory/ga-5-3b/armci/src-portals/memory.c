#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: memory.c,v 1.56.2.3 2007-04-25 23:49:55 d3p687 Exp $ */
#ifndef NEW_MALLOC
#include <stdio.h>
#include <assert.h>
#include "armcip.h"
#include "message.h"
#include "kr_malloc.h"

#define DEBUG_ 0
#define USE_MALLOC
#define USE_SHMEM_
#define SHM_UNIT 1024

static context_t ctx_localmem;
/*
static context_t ctx_mlocalmem;
*/
#if defined(SYSV) || defined(WIN32) || defined(MMAP) || defined(HITACHI)
#include "shmem.h"

#if !defined(USE_SHMEM) && (defined(HITACHI) || defined(MULTI_CTX))
#    define USE_SHMEM
#endif

#if !(defined(LAPI)||defined(QUADRICS)||defined(SERVER_THREAD)) ||\
    defined(USE_SHMEM)
#define RMA_NEEDS_SHMEM
#endif

void kr_check_local()
{
#if 0
kr_malloc_print_stats(&ctx_localmem);
#endif
kr_malloc_verify(&ctx_localmem);
}

void  armci_print_ptr(void **ptr_arr, int bytes, int size, void* myptr, int off)
{
int i;
int nproc = armci_clus_info[armci_clus_me].nslave;

    ARMCI_PR_DBG("enter",0);
    for(i=0; i< armci_nproc; i++){
      int j;
      if(armci_me ==i){
        printf("%d master =%d nproc=%d off=%d\n",armci_me, 
               armci_master,nproc, off);
        printf("%d:bytes=%d mptr=%p s=%d ",armci_me, bytes, myptr,size);
        for(j = 0; j< armci_nproc; j++)printf(" %p",ptr_arr[j]);
        printf("\n"); fflush(stdout);
      }
      armci_msg_barrier();
    }
    ARMCI_PR_DBG("exit",0);
}


/*\ master exports its address of shmem region at the beggining of that region
\*/
static void armci_master_exp_attached_ptr(void* ptr)
{
    ARMCI_PR_DBG("enter",0);
    if(!ptr) armci_die("armci_master_exp_att_ptr: null ptr",0);
    *(volatile void**)ptr = ptr;
    ARMCI_PR_DBG("exit",0);
}


/*\ Collective Memory Allocation on shared memory systems
\*/
void armci_shmem_malloc(void *ptr_arr[], armci_size_t bytes)
{
    void *myptr=NULL, *ptr=NULL;
    long idlist[SHMIDLEN];
    long size=0, offset=0;
    long *size_arr;
    void **ptr_ref_arr;
    int  i,cn, len;
    int  nproc = armci_clus_info[armci_clus_me].nslave;
    ARMCI_PR_DBG("enter",0);
    bzero((char*)ptr_arr,armci_nproc*sizeof(void*));

    /* allocate work arrays */
    size_arr = (long*)calloc(armci_nproc,sizeof(long));
    if(!size_arr)armci_die("armci_malloc:calloc failed",armci_nproc);
    /* allocate arrays for cluster address translations */

    ptr_ref_arr = calloc(armci_nclus,sizeof(void*)); /* must be zero */
    if(!ptr_ref_arr)armci_die("armci_malloc:calloc 2 failed",armci_nclus);

    /* combine all memory requests into size_arr  */
    size_arr[armci_me] = bytes;
    armci_msg_lgop(size_arr, armci_nproc, "+");

    /* determine aggregate request size on the cluster node */
    for(i=0, size=0; i< nproc; i++) size += size_arr[i+armci_master];

    /* master process creates shmem region and then others attach to it */
    if(armci_me == armci_master ){
       /* can malloc if there is no data server process and has 1 process/node*/
#      ifndef RMA_NEEDS_SHMEM
             if(nproc == 1)
                myptr = kr_malloc(size, &ctx_localmem, 0, NULL, NULL);
             else
#      endif
                myptr = Create_Shared_Region(idlist+1,size,idlist);
       if(!myptr && size>0 )armci_die("armci_malloc: could not create", (int)(size>>10));

       /* place its address at begining of attached region for others to see */
       if(size)armci_master_exp_attached_ptr(myptr);

       if(DEBUG_){
         printf("%d:armci_malloc addr mptr=%p size=%ld\n",armci_me,myptr,size);
         fflush(stdout);
       }
    }

    /* broadcast shmem id to other processes on the same cluster node */
    armci_msg_clus_brdcst(idlist, SHMIDLEN*sizeof(long));

    if(armci_me != armci_master){
       myptr=(double*)Attach_Shared_Region(idlist+1,size,idlist[0]);
       if(!myptr)armci_die("armci_malloc: could not attach", (int)(size>>10));

       /* now every process in a SMP node needs to find out its offset
        * w.r.t. master - this offset is necessary to use memlock table
        */
       if(size) armci_set_mem_offset(myptr);
       if(DEBUG_){
          printf("%d:armci_malloc attached addr mptr=%p ref=%p size=%ld\n",
                 armci_me,myptr, *(void**)myptr,size); fflush(stdout);
       }
    }
#   if defined(DATA_SERVER)
       /* get server reference address for every cluster node to perform
        * remote address translation for global address space */
       if(armci_nclus>1){
          if(armci_me == armci_master){
#            ifdef SERVER_THREAD
               ptr_ref_arr[armci_clus_me]=myptr;
#            else
               {
                  extern int _armci_server_started;
                  if(_armci_server_started) {
                     armci_serv_attach_req(idlist, SHMIDLEN*sizeof(long), size,
                                           &ptr, sizeof(void*));
                     ptr_ref_arr[armci_clus_me]= myptr; /* from server*/
                  }
                  else /* server not yet started */
                     ptr_ref_arr[armci_clus_me]=myptr;
               }

               if(DEBUG_){
                 printf("%d:addresses server=%p myptr=%p\n",armci_me,ptr,myptr);
                 fflush(stdout);
               }
#            endif
          }
          /* exchange ref addr of shared memory region on every cluster node*/
          armci_exchange_address(ptr_ref_arr, armci_nclus);
        # ifdef ARMCI_REGISTER_SHMEM
          armci_register_shmem(myptr,size,idlist+1,idlist[0],ptr_ref_arr[armci_clus_me]);
        # endif
       }else {
          ptr_ref_arr[armci_master] = myptr;
       }
       /* translate addresses for all cluster nodes */
       for(cn = 0; cn < armci_nclus; cn++){
         int master = armci_clus_info[cn].master;
         offset = 0;
         /* on local cluster node use myptr directly */
         ptr = (armci_clus_me == cn) ? myptr: ptr_ref_arr[cn];
         /* compute addresses pointing to the memory regions on cluster node*/
         for(i=0; i< armci_clus_info[cn].nslave; i++){
           /* NULL if request size is 0*/
           ptr_arr[i+master] = (size_arr[i+master])? ((char*)ptr)+offset : NULL;
           offset += size_arr[i+master];
         }
       }
#   else
      /* compute addresses for local cluster node */
      offset =0;
      for(i=0; i< nproc; i++) {
        ptr_ref_arr[i] = (size_arr[i+armci_master])? ((char*)myptr)+offset : 0L;
        offset += size_arr[i+armci_master];
      }
      /* exchange addreses with all other processes */
      ptr_arr[armci_me] = (char*)ptr_ref_arr[armci_me-armci_master]; 
      armci_exchange_address(ptr_arr, armci_nproc);
      /* overwrite entries for local cluster node with ptr_ref_arr */
      bcopy((char*)ptr_ref_arr, (char*)(ptr_arr+armci_master), nproc*sizeof(void*)); 
     /*  armci_print_ptr(ptr_arr, bytes, size, myptr, offset);*/
#   endif

    armci_msg_barrier();

    /* free work arrays */
    free(ptr_ref_arr);
    free(size_arr);
    ARMCI_PR_DBG("exit",0);

}

/********************************************************************
 * Non-collective Memory Allocation on shared memory systems
\*/
void armci_shmem_memget(armci_meminfo_t *meminfo, size_t size) {
    void *myptr=NULL;
    void *armci_ptr=NULL; /* legal ARCMIptr used in ARMCI data xfer ops */
    long idlist[SHMIDLEN];
    
    /* can malloc if there is no data server process & has 1 process/node*/
#ifndef RMA_NEEDS_SHMEM
    if( armci_clus_info[armci_clus_me].nslave == 1)
       myptr = kr_malloc(size, &ctx_localmem, 0, NULL, NULL);
    else
#endif
       myptr = Create_Shared_Region(idlist+1,size,idlist);
    
    if(!myptr && size>0 )
       armci_die("armci_shmem_memget: create failed", (int)(size>>10));
    
    if(DEBUG_)
    {
       printf("%d: armci_shmem_memget: addr=%p size=%ld %ld %ld \n", armci_me,
              myptr, size, idlist[0], idlist[1]);
       fflush(stdout);
    }

    armci_ptr = myptr;
    
#if defined(DATA_SERVER)
    
    /* get server reference address to perform
     * remote address translation for global address space */
    if(armci_nclus>1)
    {
#   ifdef SERVER_THREAD

       /* data server thread runs on master process */
       if(armci_me != armci_master) {
          armci_serv_attach_req(idlist, SHMIDLEN*sizeof(long), size,
                                &armci_ptr, sizeof(void*));
       }
       
#   else
       /* ask dataserver process to attach to region and get ptr*/
       {
          extern int _armci_server_started;
          if(_armci_server_started) {
             armci_serv_attach_req(idlist, SHMIDLEN*sizeof(long), size,
                                   &armci_ptr, sizeof(void*));
          }
       }      
#   endif
    }
#endif

    /* fill the meminfo structure */
    meminfo->armci_addr = armci_ptr;
    meminfo->addr       = myptr;
    meminfo->size       = size;
    meminfo->cpid       = armci_me;
    bcopy(idlist, meminfo->idlist, SHMIDLEN*sizeof(long));

}

void* armci_shmem_memat(armci_meminfo_t *meminfo) {
    void *ptr=NULL;
    long size    = (long)  meminfo->size;
    long *idlist = (long*) meminfo->idlist;
  
    if(SAMECLUSNODE(meminfo->cpid))
    {
       /* Attach to the shared memory segment */
       ptr=(double*)Attach_Shared_Region(idlist+1,size,idlist[0]);
       if(!ptr)armci_die("ARMCi_Memat: could not attach", (int)(size>>10));

       /* CHECK: now every process in a SMP node needs to find out its offset
        * w.r.t. master - this offset is necessary to use memlock table
        */
       if(size) armci_set_mem_offset(ptr);
    }
    else
    {
       ptr = meminfo->armci_addr; /* remote address  */
    }

    return ptr;
}

void armci_shmem_memctl(armci_meminfo_t *meminfo) {

    /* only the creator can delete the segment */
    if(meminfo->cpid == armci_me) {
       void *ptr = meminfo->addr;
    
#ifdef RMA_NEEDS_SHMEM
       Free_Shmem_Ptr(0,0,ptr);
#else
       if(armci_clus_info[armci_clus_me].nslave>1)
          Free_Shmem_Ptr(0,0,ptr);
       else kr_free(ptr, &ctx_localmem);
#endif
    }
}

/****** End: Non-collective memory allocation on shared memory systems *****/

#ifdef MPI
/********************************************************************
 * Group Memory Allocation on shared memory systems for ARMCI Groups
\*/
void armci_shmem_malloc_group(void *ptr_arr[], armci_size_t bytes,
                              ARMCI_Group *group)
{
    void *myptr=NULL, *ptr=NULL;
    long idlist[SHMIDLEN];
    long size=0, offset=0;
    long *size_arr;
    void **ptr_ref_arr;
    int  i,cn, len;
    /* int  nproc = armci_clus_info[armci_clus_me].nslave; ? change ? */
    int grp_me, grp_nproc, grp_nclus, grp_master, grp_clus_nproc, grp_clus_me;
    armci_grp_attr_t *grp_attr=ARMCI_Group_getattr(group);
    ARMCI_PR_DBG("enter",0);

    /* Get the group info: group size & group rank */
    ARMCI_Group_size(group, &grp_nproc);
    ARMCI_Group_rank(group, &grp_me);
    if(grp_me == MPI_UNDEFINED) { /* check if the process is in this group */
       armci_die("armci_malloc_group: process is not a member in this group",
                 armci_me);
    }

    grp_nclus      = grp_attr->grp_nclus;
    grp_clus_me    = grp_attr->grp_clus_me;
    grp_master     = grp_attr->grp_clus_info[grp_clus_me].master;
    grp_clus_nproc = grp_attr->grp_clus_info[grp_clus_me].nslave;

    bzero((char*)ptr_arr,grp_nproc*sizeof(void*));

    /* allocate work arrays */
    size_arr = (long*)calloc(grp_nproc,sizeof(long));
    if(!size_arr)armci_die("armci_malloc_group:calloc failed",grp_nproc);

    /* allocate arrays for cluster address translations */
#   if defined(DATA_SERVER)
        len = grp_nclus;
#   else
        len = grp_clus_nproc;
#   endif

    ptr_ref_arr = calloc(len,sizeof(void*)); /* must be zero */
    if(!ptr_ref_arr)armci_die("armci_malloc_group:calloc 2 failed",len);

    /* combine all memory requests into size_arr  */
    size_arr[grp_me] = bytes;
    armci_msg_group_gop_scope(SCOPE_ALL, size_arr, grp_nproc, "+", ARMCI_LONG,
                              group);

    /* determine aggregate request size on the cluster node */
    for(i=0, size=0; i< grp_clus_nproc; i++) size += size_arr[i+grp_master];

    /* master process creates shmem region and then others attach to it */
    if(grp_me == grp_master ){


       /* can malloc if there is no data server process and has 1 process/node*/
#     ifndef RMA_NEEDS_SHMEM
       if( armci_clus_info[armci_clus_me].nslave == 1)
         myptr = kr_malloc(size, &ctx_localmem, 0, NULL, NULL);
       else
#     endif
         myptr = Create_Shared_Region(idlist+1,size,idlist);
       if(!myptr && size>0 )
         armci_die("armci_malloc_group: could not create", (int)(size>>10));

       /* place its address at begining of attached region for others to see */
       if(size)armci_master_exp_attached_ptr(myptr);

       if(DEBUG_){
         printf("%d:armci_malloc_group addr mptr=%p ref=%p size=%ld %ld %ld \n",armci_me,myptr,*(void**)myptr, size,idlist[0],idlist[1]);
         fflush(stdout);
       }
    }

    /* broadcast shmem id to other processes (in the same group) on the 
       same cluster node */
    armci_grp_clus_brdcst(idlist, SHMIDLEN*sizeof(long), grp_master, 
                          grp_clus_nproc, group);

    if(grp_me != grp_master){
       myptr=(double*)Attach_Shared_Region(idlist+1,size,idlist[0]);
       if(!myptr)armci_die("armci_malloc_group: could not attach", (int)(size>>10));

       /* now every process in a SMP node needs to find out its offset
        * w.r.t. master - this offset is necessary to use memlock table
        */
       if(size) armci_set_mem_offset(myptr);
       if(DEBUG_){
          printf("%d:armci_malloc_group attached addr mptr=%p ref=%p size=%ld\n",
                 armci_me,myptr, *(void**)myptr,size); fflush(stdout);
       }
    }
    
#   if defined(DATA_SERVER)
 
    /* get server reference address for every cluster node in the group 
     * to perform remote address translation for global address space */
    if(grp_nclus>1){
       if(grp_me == grp_master){
 
#            ifdef SERVER_THREAD
 
          /* data server thread runs on master process */
          if(ARMCI_Absolute_id(group,grp_master)!=armci_master){
            /*printf("\n%d: grp_master=%d %ld %ld \n",armci_me,ARMCI_Absolute_id(group,grp_master),idlist[0],idlist[1]);*/
            armci_serv_attach_req(idlist, SHMIDLEN*sizeof(long), size,
                                  &ptr, sizeof(void*));
            ptr_ref_arr[grp_clus_me]= ptr; /* from server*/
          }
          else
            ptr_ref_arr[grp_clus_me]=myptr;
          
#            else
          /* ask data server process to attach to the region and get ptr */
          {
             extern int _armci_server_started;
             if(_armci_server_started) {
                armci_serv_attach_req(idlist, SHMIDLEN*sizeof(long), size,
                                      &ptr, sizeof(void*));
                ptr_ref_arr[grp_clus_me]= ptr; /* from server*/
             }
             else /* server not yet started */
                ptr_ref_arr[grp_clus_me]=myptr;
          }
 
          if(DEBUG_){
             printf("%d:addresses server=%p myptr=%p\n",grp_me,ptr,myptr);
               fflush(stdout);
          }
#            endif
       }
       /* exchange ref addr of shared memory region on every cluster node*/
       {
          int ratio = sizeof(void*)/sizeof(int);
          if(DEBUG_)printf("%d: exchanging %ld ratio=%d\n",armci_me,
                           (long)ptr_arr[grp_me], ratio);
          armci_msg_group_gop_scope(SCOPE_ALL, ptr_ref_arr, grp_nclus*ratio,
                                    "+", ARMCI_INT, group);
        # ifdef ARMCI_REGISTER_SHMEM
          armci_register_shmem_grp(myptr,size,idlist+1,idlist[0],ptr_ref_arr[armci_clus_me],group);
        # endif
       }
    }else {
       
       ptr_ref_arr[grp_master] = myptr;
       
    }
    
    /* translate addresses for all cluster nodes */
    for(cn = 0; cn < grp_nclus; cn++){
       
       int master = grp_attr->grp_clus_info[cn].master;
       offset = 0;
 
       /* on local cluster node use myptr directly */
       ptr = (grp_clus_me == cn) ? myptr: ptr_ref_arr[cn];

       /* compute addresses pointing to the memory regions on cluster node*/
       for(i=0; i< grp_attr->grp_clus_info[cn].nslave; i++){
 
          /* NULL if request size is 0*/
          ptr_arr[i+master] =(size_arr[i+master])? ((char*)ptr)+offset: NULL;
            offset += size_arr[i+master];
       }
    }
 
#   else

    /* compute addresses for local cluster node */
    offset =0;
    for(i=0; i< grp_clus_nproc; i++) {
       
       ptr_ref_arr[i] = (size_arr[i+grp_master])? ((char*)myptr)+offset : 0L;
       offset += size_arr[i+grp_master];
       
    }
    
    /* exchange addreses with all other processes */
    ptr_arr[grp_me] = (char*)ptr_ref_arr[grp_me-grp_master]; 
    armci_exchange_address_grp(ptr_arr, grp_nproc, group);

    /* overwrite entries for local cluster node with ptr_ref_arr */
    bcopy((char*)ptr_ref_arr, (char*)(ptr_arr+grp_master), grp_clus_nproc*sizeof(void*)); 
     
#   endif

    /*  armci_print_ptr(ptr_arr, bytes, size, myptr, offset);*/
       
    armci_msg_group_barrier(group);

    /* free work arrays */
    free(ptr_ref_arr);
    free(size_arr);
    ARMCI_PR_DBG("exit",0);
}
#endif /* ifdef MPI */

#else

void armci_shmem_malloc(void* ptr_arr[], int bytes)
{
  armci_die("armci_shmem_malloc should never be called on this system",0);
}
void armci_shmem_memget(armci_meminfo_t *meminfo, size_t size) {
  armci_die("armci_shmem_memget should never be called on this system",0);
}
void* armci_shmem_memat(armci_meminfo_t *meminfo) {
  armci_die("armci_shmem_memat should never be called on this system",0);
}
void armci_shmem_memctl(armci_meminfo_t *meminfo) {
  armci_die("armci_shmem_memctl should never be called on this system",0);  
}
# ifdef MPI
  void armci_shmem_malloc_group(void *ptr_arr[], armci_size_t bytes,
                                ARMCI_Group *group) {
      armci_die("armci_shmem_malloc_group should never be called on this system",0);
  }
# endif

#endif


/* public constructor to initialize the kr_malloc context */
void armci_krmalloc_init_localmem() {
#if defined(ALLOW_PIN)
    kr_malloc_init(0, 0, 0, reg_malloc, 0, &ctx_localmem);
    kr_malloc_init(0, 0, 0, malloc, 0, &ctx_mlocalmem);
    ctx_mlocalmem.ctx_type = KR_CTX_LOCALMEM;
#elif defined(CRAY_SHMEM) && defined(CRAY_XT)
#   ifdef CATAMOUNT
    int units_avail = (cnos_shmem_size() - 1024 * 1024) / SHM_UNIT;
#   else
    extern size_t get_xt_heapsize();
    int units_avail = (get_xt_heapsize() - 1024 * 1024) / SHM_UNIT;
#   endif

    if(DEBUG_) 
    {
       fprintf(stderr,"%d:krmalloc_init_localmem: symheap=%llu,units(%d)=%d\n",
               armci_me, SHM_UNIT*units_avail, SHM_UNIT, units_avail);
    }
    kr_malloc_init(SHM_UNIT, units_avail, units_avail, shmalloc, 0,
                   &ctx_localmem);
    armci_shmalloc_exchange_offsets(&ctx_localmem);
#else

    kr_malloc_init(0, 0, 0, malloc, 0, &ctx_localmem);

#endif

    ctx_localmem.ctx_type = KR_CTX_LOCALMEM;
}

/**
 * Local Memory Allocation and Free
 */
void *PARMCI_Malloc_local(armci_size_t bytes) {
    void *rptr;
    ARMCI_PR_DBG("enter",0);
    ARMCI_PR_DBG("exit",0);
    rptr = (void *)kr_malloc((size_t)bytes, &ctx_localmem, 0, NULL, NULL);
    //printf("\n%d:%s:%d:%p\n",armci_me,FUNCTION_NAME,bytes,rptr);
    return rptr;
}

int PARMCI_Free_local(void *ptr) {
    ARMCI_PR_DBG("enter",0);
    kr_free((char *)ptr, &ctx_localmem);
    ARMCI_PR_DBG("exit",0);
    return 0;
}


/*\ Collective Memory Allocation
 *  returns array of pointers to blocks of memory allocated by everybody
 *  Note: as the same shared memory region can be mapped at different locations
 *        in each process address space, the array might hold different values
 *        on every process. However, the addresses are legitimate
 *        and can be used in the ARMCI data transfer operations.
 *        ptr_arr[nproc]
\*/
int PARMCI_Malloc(void *ptr_arr[], armci_size_t bytes)
{
    void *ptr;
    char *new_base;
    size_t new_size=0;
    ARMCI_PR_DBG("enter",0);
    if(DEBUG_){ 
       fprintf(stderr,"%d bytes in armci_malloc %d\n",armci_me, (int)bytes);
       fflush(stderr);
       armci_msg_barrier();
    }

#  ifdef USE_MALLOC
    if(armci_nproc == 1) {
      ptr = kr_malloc((size_t) bytes, &ctx_localmem, 0, NULL, NULL);
      if(bytes) if(!ptr) armci_die("armci_malloc:malloc 1 failed",(int)bytes);
      ptr_arr[armci_me] = ptr;
      ARMCI_PR_DBG("exit",0);
      return (0);
    }
#  endif

   // static int one_time = 0;

    if( ARMCI_Uses_shm() ) {
    //  if(one_time++ == 0 && armci_me==0) printf("%d: ARMCI_Uses_shm = true\n",armci_me);
      armci_shmem_malloc(ptr_arr,bytes);
    } else {
      /* on distributed-memory systems just malloc & collect all addresses */
      ptr = kr_malloc(bytes, &ctx_localmem, 1, &new_base, &new_size);
      if(bytes) if(!ptr) armci_die("armci_malloc:malloc 2 failed",bytes);

      bzero((char*)ptr_arr,armci_nproc*sizeof(void*));
      ptr_arr[armci_me] = ptr;
      
      /* now combine individual addresses into a single array */
      armci_exchange_address(ptr_arr, armci_nproc);
    # ifdef ARMCI_REGISTER_SHMEM
      if(new_size)
        armci_register_shmem(new_base,new_size,NULL,0,new_base);
      else
        armci_register_shmem(ptr,bytes,NULL,0,ptr);
    # endif
    }
    ARMCI_PR_DBG("exit",0);
    //printf("\n%d:%s:%d:%p\n",armci_me,FUNCTION_NAME,bytes,ptr_arr[armci_me]);
    return(0);
}



/*\ shared memory is released to kr_malloc only on process 0
 *  with data server malloc cannot be used
\*/
int PARMCI_Free(void *ptr)
{
    ARMCI_PR_DBG("enter",0);
    if(!ptr)return 1;

#    if (defined(SYSV) || defined(WIN32) || defined(MMAP)) && !defined(NO_SHM)
#       ifdef USE_MALLOC
          if(armci_nproc > 1)
#       endif
             if(ARMCI_Uses_shm()){
                if(armci_me==armci_master){
#               ifdef RMA_NEEDS_SHMEM
                   Free_Shmem_Ptr(0,0,ptr);
#               else
                   if(armci_clus_info[armci_clus_me].nslave>1)
                      Free_Shmem_Ptr(0,0,ptr);
                   else kr_free(ptr, &ctx_localmem);
#               endif
                }
                ptr = NULL;
                return 0;
             }
#    endif
	  kr_free(ptr, &ctx_localmem);
    //armci_unregister_shmem(ptr,0);
    ptr = NULL;
    ARMCI_PR_DBG("exit",0);
    return 0;
}


int ARMCI_Uses_shm()
{
    int uses=0;

#if (defined(SYSV) || defined(WIN32) || defined(MMAP) ||defined(HITACHI)) \
    && !defined(NO_SHM)
#   ifdef RMA_NEEDS_SHMEM
      if(armci_nproc >1) uses= 1; /* always unless serial mode */
#   else
      if(armci_nproc != armci_nclus)uses= 1; /* only when > 1 node used */
#   endif
#endif
    if(DEBUG_) fprintf(stderr,"%d:uses shmem %d\n",armci_me, uses);
    return uses;
}
#ifdef MPI

int ARMCI_Uses_shm_grp(ARMCI_Group *group) 
{    
    int uses=0, grp_me, grp_nproc, grp_nclus;
    ARMCI_PR_DBG("enter",0);
    armci_grp_attr_t *grp_attr=ARMCI_Group_getattr(group);

    ARMCI_Group_size(group, &grp_nproc);
    ARMCI_Group_rank(group, &grp_me);
    grp_nclus = grp_attr->grp_nclus;
    
#if (defined(SYSV) || defined(WIN32) || defined(MMAP) ||defined(HITACHI)) \
    && !defined(NO_SHM)
#   ifdef RMA_NEEDS_SHMEM
      if(grp_nproc >1) uses= 1; /* always unless serial mode */
#   else
      if(grp_nproc != grp_nclus)uses= 1; /* only when > 1 node used */
#   endif
#endif
    if(DEBUG_) fprintf(stderr,"%d (grp_id=%d):uses shmem %d\n",armci_me, grp_me, uses);
    ARMCI_PR_DBG("exit",0);
    return uses;
}

/*\ ************** Begin Group Collective Memory Allocation ******************
 *  returns array of pointers to blocks of memory allocated by everybody
 *  Note: as the same shared memory region can be mapped at different locations
 *        in each process address space, the array might hold different values
 *        on every process. However, the addresses are legitimate
 *        and can be used in the ARMCI data transfer operations.
 *        ptr_arr[nproc]
\*/
int ARMCI_Malloc_group(void *ptr_arr[], armci_size_t bytes,
                       ARMCI_Group *group)
{
    void *ptr;
    int grp_me, grp_nproc;
    ARMCI_PR_DBG("enter",0);
    ARMCI_Group_size(group, &grp_nproc);
    ARMCI_Group_rank(group, &grp_me);
    if(DEBUG_)fprintf(stderr,"%d (grp_id=%d) bytes in armci_malloc_group %d\n",
                      armci_me, grp_me, (int)bytes);
#ifdef USE_MALLOC
    if(grp_nproc == 1) {
       ptr = kr_malloc((size_t) bytes, &ctx_localmem, 0, NULL, NULL);
       if(bytes) if(!ptr) armci_die("armci_malloc_group:malloc 1 failed",(int)bytes);
       ptr_arr[grp_me] = ptr;
       ARMCI_PR_DBG("exit",0);
       return (0);
    }
#endif
    
    if( ARMCI_Uses_shm_grp(group) ) {
#      ifdef SGIALTIX
          armci_altix_shm_malloc_group(ptr_arr,bytes,group);
#      else   
          armci_shmem_malloc_group(ptr_arr,bytes,group);
#      endif
    }
    else {
       void *new_base=NULL;
       size_t new_size=NULL;
       ptr = kr_malloc(bytes, &ctx_localmem, 1, &new_base, &new_size);
       if(bytes) if(!ptr) armci_die("armci_malloc:malloc 2 failed",bytes);
       
       bzero((char*)ptr_arr,grp_nproc*sizeof(void*));
       ptr_arr[grp_me] = ptr;
       
       /* now combine individual addresses into a single array */
       armci_exchange_address_grp(ptr_arr, grp_nproc, group);
      
    }
    ARMCI_PR_DBG("exit",0);
    return(0);
}


/*\ shared memory is released to kr_malloc only on process 0
 *  with data server malloc cannot be used
 \*/
int ARMCI_Free_group(void *ptr, ARMCI_Group *group)
{
    int grp_me, grp_nproc, grp_master, grp_clus_me;
    armci_grp_attr_t *grp_attr=ARMCI_Group_getattr(group);
    ARMCI_PR_DBG("enter",0);
    
    if(!ptr)return 1;

    ARMCI_Group_size(group, &grp_nproc);
    ARMCI_Group_rank(group, &grp_me);
    if(grp_me == MPI_UNDEFINED) { /* check if the process is in this group */
       armci_die("armci_malloc_group: process is not a member in this group",
                 armci_me);
    }
    /* get the group cluster info */
    grp_clus_me    = grp_attr->grp_clus_me;
    grp_master     = grp_attr->grp_clus_info[grp_clus_me].master;

#   if (defined(SYSV) || defined(WIN32) || defined(MMAP)) && !defined(NO_SHM)
#      ifdef USE_MALLOC
         if(grp_nproc > 1)
#      endif
       if(ARMCI_Uses_shm_grp(group)){
          if(grp_me == grp_master) {
#            ifdef RMA_NEEDS_SHMEM
             Free_Shmem_Ptr(0,0,ptr);
#            else
             if(armci_clus_info[armci_clus_me].nslave>1) Free_Shmem_Ptr(0,0,ptr);
             else kr_free(ptr, &ctx_localmem);
#            endif
          }
          ptr = NULL;
          ARMCI_PR_DBG("exit",0);
          return 0;
       }
#   endif
    kr_free(ptr, &ctx_localmem);

    ptr = NULL;
    ARMCI_PR_DBG("exit",0);
    return 0;
}
/* ***************** End Group Collective Memory Allocation ******************/

/* ************** Begin Non-Collective Memory Allocation ******************
 * Prototype similar to SysV shared memory.
 */

/**
 * CHECK: On Altix we are forced to use SysV as shmalloc is collective. We
 * may use a preallocated shmalloc memory, however, it may NOT still solve
 * our problem...
 * NOTE: "int memflg" option for future optimiztions.
 */
void PARMCI_Memget(size_t bytes, armci_meminfo_t *meminfo, int memflg) {

    void *myptr=NULL;
    void *armci_ptr=NULL; /* legal ARCMI ptr used in ARMCI data xfer ops*/
    size_t size = bytes;
    
    if(size<=0) armci_die("ARMCI_Memget: size must be > 0", (int)size);
    if(meminfo==NULL) armci_die("ARMCI_Memget: Invalid arg #2 (NULL ptr)",0);
    if(memflg!=0) armci_die("ARMCI_Memget: Invalid memflg", memflg);

    if( !ARMCI_Uses_shm() )
    {
       armci_ptr = myptr = kr_malloc(size, &ctx_localmem, 0, NULL, NULL);
       if(size) if(!myptr) armci_die("ARMCI_Memget failed", (int)size);

       /* fill the meminfo structure */
       meminfo->armci_addr = armci_ptr;
       meminfo->addr       = myptr;
       meminfo->size       = size;
       meminfo->cpid       = armci_me;
       /* meminfo->attr       = NULL; */
    }
    else
    {
       armci_shmem_memget(meminfo, size);
    }
    
    if(DEBUG_){
       printf("%d: ARMCI_Memget: addresses server=%p myptr=%p bytes=%ld\n",
              armci_me, meminfo->armci_addr, meminfo->addr, bytes);
       fflush(stdout);
    }    
}

void* PARMCI_Memat(armci_meminfo_t *meminfo, long offset) {
    void *ptr=NULL;
    
    if(meminfo==NULL) armci_die("ARMCI_Memget: Invalid arg #2 (NULL ptr)",0);
    if(offset!=0) armci_die("ARMCI_Memget: Invalid offset", offset);

    if(meminfo->cpid==armci_me) { ptr = meminfo->addr; return ptr; }

    if( !ARMCI_Uses_shm())
    {
       ptr = meminfo->addr;
    }
    else
    {
       ptr = armci_shmem_memat(meminfo);
    }
    
    if(DEBUG_)
    {
       printf("%d:ARMCI_Memat: attached addr mptr=%p size=%ld\n",
              armci_me, ptr, meminfo->size); fflush(stdout);
    }
    
    return ptr;
}

void ARMCI_Memdt(armci_meminfo_t *meminfo, long offset) {
  /**
   * Do nothing. May be we need to have reference counting in future. This
   * is to avoid the case of dangling pointers when the creator of shm
   * segment calls Memctl and other processes are still attached to this
   * segment
   */
}

void ARMCI_Memctl(armci_meminfo_t *meminfo) {

    if(meminfo==NULL) armci_die("ARMCI_Memget: Invalid arg #2 (NULL ptr)",0);

    /* only the creator can delete the segment */
    if(meminfo->cpid == armci_me)
    {
       if( !ARMCI_Uses_shm() )
       {
          void *ptr = meminfo->addr;
          kr_free(ptr, &ctx_localmem);
       }
       else
       {
          armci_shmem_memctl(meminfo);
       }
    }

    meminfo->addr       = NULL;
    meminfo->armci_addr = NULL;
    /* if(meminfo->attr!=NULL) free(meminfo->attr); */
}

/* ***************** End Non-Collective Memory Allocation ******************/

#endif
#endif
