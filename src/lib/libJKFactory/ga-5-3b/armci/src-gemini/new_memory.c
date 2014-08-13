#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/param.h>
#include "armcip.h"
#include "message.h"

#define DEBUG_ 0
#define USE_SHMEM_
#define SHM_UNIT 1024

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


/********************************************************************
 * Non-collective Memory Allocation on shared memory systems
\*/
void armci_shmem_memget(armci_meminfo_t *meminfo, size_t size) {
    void *myptr=NULL;
    void *armci_ptr=NULL; /* legal ARCMIptr used in ARMCI data xfer ops */
    
    /* can malloc if there is no data server process & has 1 process/node*/
}

void* armci_shmem_memat(armci_meminfo_t *meminfo) {
    return NULL;
}

void armci_shmem_memctl(armci_meminfo_t *meminfo) {

}

/****** End: Non-collective memory allocation on shared memory systems *****/

/**
 * Local Memory Allocation and Free
 */
void *PARMCI_Malloc_local(armci_size_t bytes) {
    void *rptr;
    ARMCI_PR_DBG("enter",0);
    ARMCI_PR_DBG("exit",0);
		return malloc(bytes);
}

int PARMCI_Free_local(void *ptr) {
    ARMCI_PR_DBG("enter",0);
		free(ptr);
    ARMCI_PR_DBG("exit",0);
    return 0;
}

/*\ A wrapper to shmget. Just to be sure that ID is not 0.
\*/
static int armci_shmget(size_t size,char *from)
{
int id;

    id = shmget(IPC_PRIVATE, size, (IPC_CREAT | 00600));

    /*attaching with id 0 somehow fails (Seen on pentium4+linux24+gm163)
     *so if id=0, shmget again. */
    while(id==0){
       /* free id=0 and get a new one */
       if(shmctl((int)id,IPC_RMID,(struct shmid_ds *)NULL)) {
         fprintf(stderr,"id=%d \n",id);
         armci_die("allocate: failed to _delete_ shared region ",id);
       }
       id = shmget(IPC_PRIVATE, size, (IPC_CREAT | 00600));
    }
    if(DEBUG_){
       printf("\n%d:armci_shmget sz=%ld caller=%s id=%d\n",armci_me,(long)size,
               from,id);
       fflush(stdout);
    }
    return(id);
}


/*\ Collective Memory Allocation
 *  returns array of pointers to blocks of memory allocated by everybody
 *  Note: as the same shared memory region can be mapped at different locations
 *        in each process address space, the array might hold different values
 *        on every process. However, the addresses are legitimate
 *        and can be used in the ARMCI data transfer operations.
 *        ptr_arr[nproc]
\*/
#define CLEANUP_CMD(command) sprintf(command,"/usr/bin/ipcrm shm %d",id);
int PARMCI_Malloc(void *ptr_arr[], armci_size_t bytes)
{
int mynslave = armci_clus_info[armci_clus_me].nslave;
void *servptr,*mynodeptrs[mynslave];
int id,nodeids[mynslave],mynodeid=armci_me-armci_master;

    ARMCI_PR_DBG("enter",0);
#ifdef DEBUG_MEM
    fprintf(stderr,"%d bytes in armci_malloc %d\n",armci_me, (int)bytes);
    fflush(stderr);
    armci_msg_barrier();
#endif
		if(bytes>0){
		  if(mynslave>1){

#ifdef DEBUG_MEM
        printf("\n%d:%s:mynslave is %d",armci_me,__FUNCTION__,mynslave);fflush(stdout);
#endif
        bzero((void *)nodeids,sizeof(int)*mynslave);
        id =nodeids[mynodeid]= armci_shmget(bytes,"PARMCI_Malloc");
        armci_msg_gop_scope(SCOPE_NODE,nodeids,mynslave,"+",ARMCI_INT);
        for(int i=0;i<mynslave;i++){
          if((long)((mynodeptrs[i] = shmat(nodeids[i],mynodeptrs[i],0))) == -1L){
            char command[64];
            CLEANUP_CMD(command);
            if(system(command) == -1) 
              printf("clean shared memory (id=%d): see man ipcrm\n",nodeids[i]);
            armci_die("allocate: failed to attach to shared id=",nodeids[i]);
          }
        }

        /* mark my id for rem so OS cleans when all attached processes vanish*/
        if(shmctl( id, IPC_RMID, (struct shmid_ds *)NULL))
          fprintf(stderr,"failed to remove shm id=%d\n",id);

#ifdef DEBUG_MEM
          printf("%d:attach:id=%d paddr=%p size=%ld\n",armci_me,id,mynodeptrs[mynodeid],bytes);
          fflush(stdout);
#endif

			  if(armci_nclus>1){
          servptr = armci_server_ptr(id);
			  }
        else servptr = mynodeptrs[mynodeid];

		  }
			else{
#ifdef DEBUG_MEM
        printf("\n%d:%s:mynslave is %d, doing malloc",armci_me,__FUNCTION__,mynslave);fflush(stdout);
#endif
        mynodeptrs[mynodeid] = servptr = malloc(bytes);
			}
		}
		else{
        mynodeptrs[mynodeid] = servptr = NULL;    
		}

    bzero((char*)ptr_arr,armci_nproc*sizeof(void*));
    /*ptr_arr[armci_me] = servptr;*/
    ptr_arr[armci_me] = mynodeptrs[mynodeid];
    armci_exchange_address(ptr_arr,armci_nproc);

    if(mynslave>1)for(int i=0;i<mynslave;i++){
      ptr_arr[armci_master+i] = mynodeptrs[i];
    }

    if(armci_nclus>1){
      armci_portals_memsetup((long)servptr-(long)ptr_arr[armci_me]);
    }

    ARMCI_PR_DBG("exit",0);
    return(0);

}



int PARMCI_Free(void *ptr)
{
    ARMCI_PR_DBG("enter",0);
    if(!ptr)return 1;

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
    
    ARMCI_PR_DBG("exit",0);
    return(0);
}


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
    
    if(size<=0) armci_die("PARMCI_Memget: size must be > 0", (int)size);
    if(meminfo==NULL) armci_die("PARMCI_Memget: Invalid arg #2 (NULL ptr)",0);
    if(memflg!=0) armci_die("PARMCI_Memget: Invalid memflg", memflg);

    if( !ARMCI_Uses_shm() )
    {

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
       printf("%d: PARMCI_Memget: addresses server=%p myptr=%p bytes=%ld\n",
              armci_me, meminfo->armci_addr, meminfo->addr, bytes);
       fflush(stdout);
    }    
}

void* PARMCI_Memat(armci_meminfo_t *meminfo, long memflg) {
    void *ptr=NULL;
    
    if(meminfo==NULL) armci_die("PARMCI_Memget: Invalid arg #2 (NULL ptr)",0);
    if(memflg!=0) armci_die("PARMCI_Memget: Invalid memflg", memflg);

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
       printf("%d:PARMCI_Memat: attached addr mptr=%p size=%ld\n",
              armci_me, ptr, meminfo->size); fflush(stdout);
    }
    
    return ptr;
}

void ARMCI_Memdt(armci_meminfo_t *meminfo, int memflg) {
  /**
   * Do nothing. May be we need to have reference counting in future. This
   * is to avoid the case of dangling pointers when the creator of shm
   * segment calls Memctl and other processes are still attached to this
   * segment
   */
}

void ARMCI_Memctl(armci_meminfo_t *meminfo) {

    if(meminfo==NULL) armci_die("PARMCI_Memget: Invalid arg #2 (NULL ptr)",0);

    /* only the creator can delete the segment */
    if(meminfo->cpid == armci_me)
    {
       if( !ARMCI_Uses_shm() )
       {
          void *ptr = meminfo->addr;
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
