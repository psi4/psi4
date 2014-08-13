#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: shmem.c,v 1.87.2.2 2007-09-10 23:31:32 manoj Exp $ */
/* System V shared memory allocation and managment
 *
 * Interface:
 * ~~~~~~~~~
 *  char *Create_Shared_Region(long *idlist, long size, long *offset)
 *       . to be called by just one process. 
 *       . calls kr_malloc,  malloc-like memory allocator from the K&R book. 
 *         kr_malloc inturn calls armci_allocate() that does shmget() and shmat(). 
 *       . idlist might be just a pointer to integer or a true array in the
 *         MULTIPLE_REGIONS versions (calling routine has to take care of it) 
 *  char *Attach_Shared_Region(long *idlist, long size, long offset)
 *       . called by any other process to attach to existing shmem region or
 *         if attached just calculate the address based on the offset
 *       . has to be called after shmem region was created
 *  void  Free_Shmem_Ptr(long id, long size, char* addr)
 *       . called ONLY by the process that created shmem region (id) to return
 *         pointer to kr_malloc (shmem is not destroyed)
 *  void  Delete_All_Regions()
 *       . destroys all shared memory regions
 *       . can be called by any process assuming that all processes attached
 *         to alllocated shmem regions 
 *       . needs to by called by cleanup procedure(s)
 *
 * Jarek Nieplocha, 06.13.94
 * 
 */

#ifdef SYSV
 
 
#define DEBUG_ 0
#define DEBUG1 0
#define DEBUG2_ 0

/* For debugging purposes at the beginning of the shared memory region
 * creator process can write a stamp which then is read by attaching processes
 * NOTE: on clusters we cannot use it anymore since ARMCI node master
 * uses it since Nov 99 to write the value of address it attached at
 * This feature is used in the ARMCI memlock table.
 */
#define STAMP 0


#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/sem.h>
#include <sys/param.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include "shmem.h"
#include "kr_malloc.h"
#include "shmlimit.h"
#include "message.h"
#include "armcip.h"

#ifdef   ALLOC_MUNMAP
#include <sys/mman.h>
#include <unistd.h>
static  size_t pagesize=0;
static  int logpagesize=0;
/* allow only that big shared memory segment (in MB)- incresed from 128 11/02 */
#define MAX_ALLOC_MUNMAP 128
#define MAX_ALLOC_MUNMAP_ 368
static long max_alloc_munmap=MAX_ALLOC_MUNMAP;
#endif

#if defined(SUN)
  extern char *shmat();
#endif

#define SHM_UNIT (1024)


/* Need to determine the max shmem segment size. There are 2 alternatives:
 * 1. use predefined SHMMAX if available or set some reasonable values, or
 * 2. trial-and-error search for a max value (default)
 *    case a) fork a process to determine shmmax size (more accurate)
 *    case b) search w/o forking until success (less accurate)
 */

/* under Myrinet GM, we cannot fork */
#if defined(GM) || defined(VAPI)
#   define SHMMAX_SEARCH_NO_FORK 
#endif
#if defined(LAPI) || defined(AIX) || defined(SHMMAX_SEARCH_NO_FORK) || defined(CRAY_XT) || defined(CRAY_UGNI)
#   define NO_SHMMAX_SEARCH
#endif

/* Limits for the largest shmem segment are in Kilobytes to avoid passing
 * Gigavalues to kr_malloc
 * the limit for the KSR is lower than SHMMAX in sys/param.h because
 * shmat would fail -- SHMMAX cannot be trusted (a bug)
 */
#define _SHMMAX 4*1024

#if defined(SUN)||defined(SOLARIS)
#  undef _SHMMAX
#  define _SHMMAX (1024)  /* memory in KB */
#elif defined(SGI64) || defined(AIX) || defined(CONVEX)
#  undef _SHMMAX
#  define _SHMMAX ((unsigned long)512*1024)
#elif defined(SGI) && !defined(SGI64)
#  undef _SHMMAX
#  define _SHMMAX ((unsigned long)128*1024)
#elif defined(KSR)
#  undef _SHMMAX
#  define _SHMMAX ((unsigned long)512*1024)
#elif defined(HPUX)
#  undef _SHMMAX
#  define _SHMMAX ((unsigned long)64*1024)
#elif defined(__FreeBSD__)
#  undef _SHMMAX
#  define _SHMMAX ((unsigned long)3*1024)
#elif defined(LINUX) 
#  if !defined(SHMMAX) /* Red Hat does not define SHMMAX */
#     undef _SHMMAX
#     if defined(__sparc__) || defined(__powerpc__) 
#       define _SHMMAX ((unsigned long)16*1024)
#     elif defined(__alpha__)
#       define _SHMMAX ((unsigned long)4072)
#     else
        /* Intel */
#       define _SHMMAX ((unsigned long)32*1024)
#     endif
#  endif
#elif defined(SHMMAX)
#  undef _SHMMAX
#  define _SHMMAX (((unsigned long)SHMMAX)>>10)
#endif

static  unsigned long MinShmem_per_core = 0;
static  unsigned long MaxShmem_per_core = 0;
static  unsigned long MinShmem = _SHMMAX;  
static  unsigned long MaxShmem = MAX_REGIONS*_SHMMAX;
static  context_t ctx_shmem; /* kr_malloc context */ 
static  context_t *ctx_shmem_global; /* kr_malloc context stored in shmem */
static  int create_call=0;

#ifdef  SHMMAX_SEARCH_NO_FORK
static  char *ptr_search_no_fork = (char*)0;
static  int id_search_no_fork=0;
#endif


#ifdef LINUX
#define CLEANUP_CMD(command) sprintf(command,"/usr/bin/ipcrm shm %d",id);
#elif  defined(SOLARIS) 
#define CLEANUP_CMD(command) sprintf(command,"/bin/ipcrm -m %d",id);
#elif  defined(SGI) 
#define CLEANUP_CMD(command) sprintf(command,"/usr/sbin/ipcrm -m %d",id);
#else
#define CLEANUP_CMD(command) sprintf(command,"/usr/bin/ipcrm -m %d",id);
#endif


#ifdef   ALLOC_MUNMAP
#ifdef QUADRICS
#  include <elan/elan.h>
#  include <elan3/elan3.h>
   static  char *armci_elan_starting_address = (char*)0;

#  ifdef __ia64__
#    define ALLOC_MUNMAP_ALIGN 1024*1024
#  else
#    define ALLOC_MUNMAP_ALIGN 64*1024
#  endif

#  define ALGN_MALLOC(s,a) elan_allocMain(elan_base->state, (a), (s))
#else 
#  define ALGN_MALLOC(s,a) malloc((s))
#endif

static char* alloc_munmap(size_t size)
{
char *tmp;
unsigned long iptr;
size_t bytes = size+pagesize-1;
 
    if(armci_elan_starting_address){
       tmp = armci_elan_starting_address;
       armci_elan_starting_address += size; 
#      ifdef ALLOC_MUNMAP_ALIGN
         armci_elan_starting_address += ALLOC_MUNMAP_ALIGN;
#      endif
       if(DEBUG_) {printf("%d: address for shm attachment is %p size=%ld\n",
                         armci_me,tmp,(long)size); fflush(stdout); }
    } else {
      tmp = ALGN_MALLOC(bytes, getpagesize());
      if(tmp){
        iptr = (unsigned long)tmp + pagesize-1;
        iptr >>= logpagesize; iptr <<= logpagesize;
        if(DEBUG_) printf("%d:unmap ptr=%p->%p size=%d pagesize=%d\n",armci_me, 
                          tmp,(char*)iptr,(int)size,pagesize);
        tmp = (char*)iptr;
        if(munmap(tmp, size) == -1) armci_die("munmap failed",0);
        if(DEBUG_){printf("%d: unmap OK\n",armci_me); fflush(stdout);}
      }else armci_die("alloc_munmap: malloc failed",(int)size);
    }
    return tmp;
}
#endif

/*\ A wrapper to shmget. Just to be sure that ID is not 0.
\*/
int armci_shmget(size_t size,char *from)
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

static int 
Semget(key_t key,int nsems,int semflg) {
      int ret;

      if((ret = semget(key,nsems,semflg)) == -1) {
         switch(errno) {
           case EACCES: fprintf(stdout," semget errno=EACCES.\n"); break;
           case EINVAL: fprintf(stdout," semget errno=EINVAL.\n"); break;
           case ENOENT: fprintf(stdout," semget errno=ENOENT.\n"); break;
           case ENOSPC: fprintf(stdout," semget errno=ENOSPC -- check system limit for sysv semaphores.\n"); break;
           case ENOMEM: fprintf(stdout," semget errno=ENOMEM.\n"); break;
           case EEXIST: fprintf(stdout," semget errno=EEXIST.\n"); break;
           default:
             fprintf(stdout," unusual semget errno=%d\n",errno); break;
         }
         armci_die("semget failed",errno);
      }

      return ret;
}


int armci_semget(int count) 
{
        int id;
        id = Semget(IPC_PRIVATE,2,0600);
        return id;
}

int armci_semrm(int id) 
{
        semctl(id,0,IPC_RMID);
}

int armci_shmrm(int id)
{
        int ret;
        if((ret = shmctl(id,IPC_RMID,NULL)) != 0) {
          fprintf(stdout,"[cp]: shmctl return an error.\n");
          switch(errno) {
            case EINVAL:
              fprintf(stdout," Error EINVAL: shmid is not a valid shared memory segment.\n");
              break;
            case EFAULT:
              fprintf(stdout," Error EFAULT: argument 3 is not a valid struct shmid_ds.\n");
              break;
            case EPERM:
              fprintf(stdout," Error EPREM: permission to access/change shared mem segment denied.\n");
              break;
            default:
              fprintf(stdout," unusual shmctl errno=%d\n",errno); break;
              break;
            }
            armci_die("error deleting shmid",id);
        }
        return ret;
}

/*\ test is a shared memory region of a specified size can be allocated
 *  return 0 (no) or 1 (yes)
\*/
int armci_test_allocate(long size)
{
   char *ptr;
   int id = armci_shmget((size_t)size,"armci_test_allocate");
   if (id <0) return 0;

   /* attach to segment */
   ptr =  shmat(id, (char *) NULL, 0);

   /* delete segment id */
   if(shmctl( id, IPC_RMID, (struct shmid_ds *)NULL))
      fprintf(stderr,"failed to remove shm id=%d\n",id);

   /* test pointer */
   if (((long)ptr) == -1L) return 0;
   else return 1;
}


/*\ try to allocate a shared memory region of a specified size; return pointer
\*/
static int armci_shmalloc_try(long size)
{
#ifdef  SHMMAX_SEARCH_NO_FORK
   char *ptr;
   int id = armci_shmget((size_t) size,"armci_shmalloc_try");
   if (id <0) return 0;

   /* attach to segment */
   ptr =  shmat(id, (char *) NULL, 0);

   /* test pointer */
   if (((long)ptr) == -1L) return 0;

   ptr_search_no_fork = ptr;
   id_search_no_fork = id;
#endif
   return 1;
}




/* parameters that define range and granularity of search for shm segment size
 * UBOUND is chosen to be < 2GB to avoid overflowing on 32-bit systems
 * smaller PAGE gives more accurate results but with more search steps
 * LBOUND  is set to minimum amount for our purposes
 * change UBOUND=512MB if you need larger arrays than 512 MB
 */
#define PAGE (16*65536L)
#define LBOUND  1048576L
#if defined(MULTI_CTX) && defined(QUADRICS)
#define UBOUND (256*LBOUND)
#else
#define UBOUND (512*LBOUND)
#endif

static long get_user_shmmax()
{
char *uval;
long x=0;
     uval = getenv("ARMCI_DEFAULT_SHMMAX"); 
     if(uval != NULL){
       sscanf(uval,"%ld",&x);
       if(x<1L || x> 2048L){ 
          fprintf(stderr,"incorrect ARMCI_DEFAULT_SHMMAX should be <1,2048>mb and 2^N Found=%ld\n",x);
          x=0;
       }
     }
     return x*1048576; /* return value in bytes */
}

/*\ determine the max shmem segment size using bisection
\*/
int armci_shmem_test()
{
long x;
int  i,rc;
long upper_bound=UBOUND;
long lower_bound=0;

     x = get_user_shmmax();
     if(!x) x = upper_bound;
     else upper_bound =x;
     
     if(DEBUG_){printf("%d: x = %ld upper_bound=%ld\n",armci_me, x, upper_bound); fflush(stdout);}

     for(i=1;;i++){
        long step;
        rc = armci_test_allocate(x);
        if(DEBUG_) 
           printf("%d:test %d size=%ld bytes status=%d\n",armci_me,i,x,rc);
        if(rc){
          lower_bound = x;
          step = (upper_bound -x)>>1;
          if(step < PAGE) break;
          x += step;
        }else{
          upper_bound = x;
          step = (x-lower_bound)>>1;
          if(step<PAGE) break;
          x -= step;
        }
        /* round it up to a full base-2 MB */
        x += 1048576L -1L;
        x >>=20;
        x <<=20; 
      }

      if(!lower_bound){
          /* try if can get LBOUND - necessary if search starts from UBOUND */
          lower_bound=LBOUND;
          rc = armci_test_allocate(lower_bound);
          if(!rc) return(0);
      }

      if(DEBUG_) printf("%ld bytes segment size, %d calls \n",lower_bound,i);
      return (int)( lower_bound>>20); /* return shmmax in mb */
}


#ifdef SHMMAX_SEARCH_NO_FORK
/*\ determine the max shmem segment size by halving
\*/
static int armci_shmem_test_no_fork()                          
{
long x;                                                     
int  i,rc;
long lower_bound=_SHMMAX*SHM_UNIT;
#define UBOUND_SEARCH_NO_FORK (256*SHM_UNIT*SHM_UNIT)

     x = get_user_shmmax();
     if(!x) x = UBOUND_SEARCH_NO_FORK;

     for(i=1;;i++){

        rc = armci_shmalloc_try(x);
        if(DEBUG_)
           printf("%d:test by halving size=%ld bytes rc=%d\n",armci_me,x,rc);

        if(rc){
          lower_bound = x;
          break;
        }else{
          x >>= 1 ;
          if(x<lower_bound) break;
        }
     }

     if(DEBUG_) printf("%ld: shmax test no fork: bytes segment size, %d calls \n",lower_bound,i);
     return (int)( lower_bound>>20); /* return shmmax in mb */
}
#endif


#ifdef MULTI_CTX
void armci_nattach_preallocate_info(int* segments, int *segsize)
{
     int x;
     char *uval;
     uval = getenv("LIBELAN_NATTACH");
     if(uval != NULL){
        sscanf(uval,"%d",&x);
        if(x<2 || x>8) armci_die("Error in LIBELAN_NATTACH <8, >1 ",(int)x);
     }else
        armci_die("Inconsistent configuration: ARMCI needs LIBELAN_NATTACH",0);
     *segments =x;
     *segsize = (int) (SHM_UNIT * MinShmem);

}
#endif
        
/* Create shared region to store kr_malloc context in shared memory */
void armci_krmalloc_init_ctxshmem() {
    void *myptr=NULL;
    long idlist[SHMIDLEN];
    long size; 
    int offset = sizeof(void*)/sizeof(int);

    /* to store shared memory context and  myptr */
    size = SHMEM_CTX_MEM;
    
    if(armci_me == armci_master ){
       myptr = Create_Shared_Region(idlist+1,size,idlist);
       if(!myptr && size>0 ) armci_die("armci_krmalloc_init_ctxshmem: could not create", (int)(size>>10));
       if(size) *(volatile void**)myptr = myptr;
       if(DEBUG_){
	  printf("%d:armci_krmalloc_init_ctxshmem addr mptr=%p ref=%p size=%ld\n", armci_me, myptr, *(void**)myptr, size);
	  fflush(stdout);
       }
       
       /* Bootstrapping: allocate storage for ctx_shmem_global. NOTE:there is 
	  offset,as master places its address at begining for others to see */
       ctx_shmem_global = (context_t*) ( ((int*)myptr)+offset );
       *ctx_shmem_global = ctx_shmem; /*master copies ctx into shared region */
    }

    /* broadcast shmem id to other processes on the same cluster node */
    armci_msg_clus_brdcst(idlist, SHMIDLEN*sizeof(long));

    if(armci_me != armci_master){
       myptr=(double*)Attach_Shared_Region(idlist+1,size,idlist[0]);
       if(!myptr)armci_die("armci_krmalloc_init_ctxshmem: could not attach", (int)(size>>10));
       
       /* now every process in a SMP node needs to find out its offset
        * w.r.t. master - this offset is necessary to use memlock table
        */
       if(size) armci_set_mem_offset(myptr);
       if(DEBUG_){
          printf("%d:armci_krmalloc_init_ctxshmem attached addr mptr=%p ref=%p size=%ld\n", armci_me,myptr, *(void**)myptr,size); fflush(stdout);
       }
       /* store context info */
       ctx_shmem_global = (context_t*) ( ((int*)myptr)+offset );
       if(DEBUG_){
	  printf("%d:armci_krmalloc_init_ctxshmem: shmid=%d off=%ld size=%ld\n", armci_me, ctx_shmem_global->shmid, ctx_shmem_global->shmoffset,
		 (long)ctx_shmem_global->shmsize);
	  fflush(stdout);
       }
    }
}

void armci_shmem_init()
{

#ifdef ALLOC_MUNMAP

#if defined(QUADRICS) 
#   if (defined(__ia64__) || defined(__alpha)) && !defined(DECOSF) 

      /* this is to determine size of Elan Main memory allocator for munmap */
      long x;
      char *uval;
      uval = getenv("LIBELAN_ALLOC_SIZE");
      if(uval != NULL){
        sscanf(uval,"%ld",&x);
        if((x>80000000) && (x< 4*1024*1024*1024L)){ 
          max_alloc_munmap = (x>>20) - 72;
          if(DEBUG_){
            printf("%d: max_alloc_munmap is %ld\n",armci_me,max_alloc_munmap);
            fflush(stdout);
          }
        }
      }

      /* an alternative approach is to use MMAP area where we get
         the address from the Elan environment variable in qsnetlibs 1.4+  */
      uval = getenv("LIBELAN3_MMAPBASE");
      if(uval != NULL){
         sscanf(uval,"%p",&armci_elan_starting_address);
      }

#   endif
#   if defined(__ia64__)
       /* need aligment on 1MB boundary rather than the actual pagesize */
       pagesize = 1024*1024;
       logpagesize = 20;
#   else
       /* determine log2(pagesize) needed for address alignment */
       int tp=512;
       logpagesize = 9;
       pagesize = getpagesize();
       if(tp>pagesize)armci_die("armci_shmem_init:pagesize",pagesize);

       while(tp<pagesize){
         tp <<= 1;
         logpagesize++;
       }
       if(tp!=pagesize)armci_die("armci_shmem_init:pagesize pow 2",pagesize);
#   endif

   if(DEBUG_) {
     printf("page size =%d log=%d\n",pagesize,logpagesize); fflush(stdout); }

#endif
#endif

   if(armci_me == armci_master){
#if !defined(NO_SHMMAX_SEARCH) || defined(SHMMAX_SEARCH_NO_FORK)
#       ifdef SHMMAX_SEARCH_NO_FORK
          int x = armci_shmem_test_no_fork();
#       else
          int x = armci_child_shmem_init();
#       endif

        if(x<1)
          armci_die("no usable amount of shared memory available: only got \n",
          (int)LBOUND);

#       if defined(ALLOC_MUNMAP)
          /* cap down for special memory allocator unless ARMCI_DEFAULT_SHMMAX
             not set - the user knows what is doing*/
#         if !defined(REGION_ALLOC)
          if(!getenv("ARMCI_DEFAULT_SHMMAX"))
            if(x>max_alloc_munmap && !armci_elan_starting_address) x=max_alloc_munmap;
#         else
	    x = 10; /* mb */
#         endif
#       endif

        if(DEBUG_){
           printf("%d:shmem_init: %d mbytes max segment size\n",armci_me,x);fflush(stdout);}

        MinShmem = (long)(x<<10); /* make sure it is in kb: mb <<10 */ 
        MaxShmem = MAX_REGIONS*MinShmem;
#       ifdef REPORT_SHMMAX
              printf("%d using x=%d SHMMAX=%ldKB\n", armci_me,x, MinShmem);
              fflush(stdout);
#       endif
#else

      /* nothing to do here - limits were given */

#endif
    }

    armci_krmalloc_init_ctxshmem();
    if(DEBUG_)printf("%d: out of shmem_init\n",armci_me);
}

void armci_set_shmem_limit_per_node(int nslaves)
{
     if (MaxShmem_per_core > 0) MaxShmem = nslaves*MaxShmem_per_core;
     if (MinShmem_per_core > 0) MinShmem = nslaves*MinShmem_per_core;
}

void armci_set_shmem_limit_per_core(unsigned long shmemlimit)
{
     MaxShmem_per_core = (shmemlimit + SHM_UNIT - 1)/SHM_UNIT;
     MinShmem_per_core = (shmemlimit + SHM_UNIT - 1)/SHM_UNIT;
}
/*\ application can reset the upper limit (bytes) for memory allocation
\*/
void armci_set_shmem_limit(unsigned long shmemlimit)
{
     unsigned long kbytes;
     kbytes = (shmemlimit + SHM_UNIT -1)/SHM_UNIT;
     if(MaxShmem > kbytes) MaxShmem = kbytes;
     if(MinShmem > kbytes) MinShmem = kbytes;
}


static void shmem_errmsg(size_t size)
{
long sz=(long)size;
    printf("******************* ARMCI INFO ************************\n");
    printf("The application attempted to allocate a shared memory segment ");
    printf("of %ld bytes in size. This might be in addition to segments ",sz);
    printf("that were allocated succesfully previously. ");
    printf("The current system configuration does not allow enough ");
    printf("shared memory to be allocated to the application.\n");
    printf("This is most often caused by:\n1) system parameter SHMMAX ");
    printf("(largest shared memory segment) being too small or\n");
    printf("2) insufficient swap space.\n");
    printf("Please ask your system administrator to verify if SHMMAX ");
    printf("matches the amount of memory needed by your application and ");
    printf("the system has sufficient amount of swap space. ");
    printf("Most UNIX systems can be easily reconfigured ");
    printf("to allow larger shared memory segments,\n");
    printf("see http://www.emsl.pnl.gov/docs/global/support.shtml\n");
    printf("In some cases, the problem might be caused by insufficient swap space.\n");
    printf("*******************************************************\n");
}


static struct shm_region_list{
   char     *addr;
   long     id;
   long     sz;
   long     attached;
}region_list[MAX_REGIONS];
static int alloc_regions=0;
static long occup_blocks=0;

/* Terminology
 *   region - actual piece of shared memory allocated from OS
 *   block  - a part of allocated shmem that is given to the requesting process
 */


static int last_allocated=-1;


unsigned long armci_max_region()
{
  return MinShmem;
}


int find_regions(char *addrp,  long* id, int *region)
{
int nreg, reg;

    if(last_allocated!=-1){
       reg=last_allocated;
       last_allocated = -1;
    } else{
       
       for(reg=-1,nreg=0;nreg<alloc_regions; nreg++)
       {
          if(addrp >= region_list[nreg].addr &&
             addrp < (region_list[nreg].addr + region_list[nreg].sz))
          {
             reg = nreg;
             break;
          }
       }
       
       if(reg == -1)
          armci_die("find_regions: failed to locate shared region", 0L);
    }

    *region = reg;
    *id = region_list[reg].id;

    return 1;
}

/* returns the shmem info based on the addr */
int armci_get_shmem_info(char *addrp,  int* shmid, long *shmoffset, 
			 size_t *shmsize) 
{    
    int region; long id;

    find_regions(addrp, &id, &region);
    *shmid     = id;
    *shmoffset = (long)(addrp - region_list[region].addr);
    *shmsize   = region_list[region].sz;

    return 1;
}

long armci_shm_reg_size(int i, long id)
{
     if(i<0 || i>= MAX_REGIONS)armci_die("armci_shmem_reg_size: bad i",i); 
     return region_list[i].sz;
}

void* armci_shm_reg_ptr(int i)
{
     if(i<0 || i>= MAX_REGIONS)armci_die("armci_shmem_reg_ptr: bad i",i); 
     return (void *)region_list[i].addr;
}

Header *armci_get_shmem_ptr(int shmid, long shmoffset, size_t shmsize) 
{
/* returns, address of the shared memory region based on shmid, offset.
 * (i.e. return_addr = stating address of shmid + offset)*/
    long idlist[SHMIDLEN];
    Header *p = NULL;

    idlist[1] = (long)shmid;
    idlist[0] = shmoffset;
    idlist[IDLOC+1] = shmsize; /* CHECK : idlist in CreateShmem????*/

    if(!(p=(Header*)Attach_Shared_Region(idlist+1, shmsize, idlist[0])))
       armci_die("kr_malloc:could not attach",(int)(p->s.shmsize>>10));
#if DEBUG_
    printf("%d: armci_get_shmem_ptr: %d %ld %ld %p\n",
           armci_me, idlist[1], idlist[0], shmsize, p);
    fflush(stdout);    
#endif
    return p;
}


char *Attach_Shared_Region(id, size, offset)
     long *id, offset, size;
{
int reg, found, shmflag=0;
static char *temp;

  if(alloc_regions>=MAX_REGIONS)
       armci_die("Attach_Shared_Region: to many regions ",0);

  if(DEBUG_){
      printf("%d:AttachSharedRegion %d:size=%ld id=%ld\n",
             armci_me, create_call++, size,*id);
      fflush(stdout);
  }


  /* under Linux we can get valid id=0 */
#ifndef LINUX
  if(!*id) armci_die("Attach_Shared_Region: shmem ID=0 ",(int)*id);
#endif

  /* first time needs to initialize region_list structure */
  if(!alloc_regions){
      for(reg=0;reg<MAX_REGIONS;reg++){
        region_list[reg].addr=(char*)0;
        region_list[reg].attached=0;
        region_list[reg].id=0;
      }
      MinShmem= id[SHMIDLEN-2];
      if(DEBUG2_){
         printf("%d:attach: allocation unit: %ldK\n",armci_me,MinShmem);
         fflush(stdout);
      }
  }

  /* search region_list for the current shmem id */
  for(found = 0, reg=0; reg < MAX_REGIONS;reg++)
      if((found=(region_list[reg].id == *id)))break;

  if(!found){
     reg = alloc_regions;
     region_list[reg].id =*id;
     alloc_regions++;
  }

  /* we need to use the actual shared memory segment not user req size */ 
  size = id[IDLOC];

  /* attach if not attached yet */
  if(!region_list[reg].attached){

#   ifdef ALLOC_MUNMAP
       char *pref_addr = alloc_munmap((size_t) (size));
#   else
       char *pref_addr = (char*)0;
#   endif
    if ( (long) (temp = shmat((int) *id, pref_addr, shmflag)) == -1L){
       fprintf(stderr,"%d:attach error:id=%ld off=%ld seg=%ld\n",armci_me,*id,offset,MinShmem);
       shmem_errmsg((size_t)MinShmem*1024);
       armci_die("Attach_Shared_Region:failed to attach to segment id=",(int)*id);
    }
    if(DEBUG_){
        printf("%d:attached: id=%d address=%p\n",armci_me,(int)*id, temp);
        fflush(stdout);
    }
    POST_ALLOC_CHECK(temp,size);
    region_list[reg].addr = temp; 
    region_list[reg].attached = 1;
    region_list[reg].sz= size;
    if(DEBUG_){
    printf("\n%d:%s:new memory region\n",armci_me,__FUNCTION__);fflush(stdout);
    }

  }

  if(STAMP)
  /* check stamp to make sure that we are attached in the right place */
  if(*((int*)(region_list[reg].addr+ offset))!= alloc_regions-1)
      armci_die("Attach_Shared_Region: wrong stamp value !", 
                *((int*)(region_list[reg].addr+ offset)));
  occup_blocks++;
  return (region_list[reg].addr+ offset);
}

#ifdef ALLOW_PIN
extern void armci_region_register_shm(void *start, long size);
#endif

/*\ allocates shmem, to be called by krmalloc that is called by process that
 *  creates shmem region
\*/
void *armci_allocate(long size)
{
char * temp;
int id,shmflag=0;
size_t sz = (size_t)size;
#ifdef ALLOC_MUNMAP
       char *pref_addr = alloc_munmap((size_t) (MinShmem*SHM_UNIT));
#else
       char *pref_addr = (char*)0;
#endif
#if defined(SGI_N32) && defined(SHM_SGI_ANYADDR)
  shmflag= SHM_SGI_ANYADDR;
#endif

    if(DEBUG1){
       printf("%d:allocate: Shmem allocate size %ld bytes\n",armci_me,size); 
       fflush(stdout);
    }

    if( alloc_regions >= MAX_REGIONS)
       armci_die("Create_Shared_Region:allocate:too many regions allocated ",0);

    last_allocated = alloc_regions;

#ifdef SHMMAX_SEARCH_NO_FORK
    if (ptr_search_no_fork){
        temp = ptr_search_no_fork;
        id   = id_search_no_fork;
        ptr_search_no_fork = (char*)0; /* do not look at it again */
    }else 
#endif
    {
       if ( (id = armci_shmget(sz,"armci_allocate")) < 0 ) {
          fprintf(stderr,"id=%d size=%ld\n",id, size);
          shmem_errmsg(sz);
          armci_die("allocate: failed to create shared region ",id);
       }

       if ( (long)( (temp = shmat(id, pref_addr, shmflag))) == -1L){
          char command[64];
          CLEANUP_CMD(command);
          if(system(command) == -1) 
            printf("Please clean shared memory (id=%d): see man ipcrm\n",id);
          armci_die("allocate: failed to attach to shared region id=",id);
       }
       if(DEBUG_){
         printf("%d:allocate:attach:id=%d paddr=%p size=%ld\n",armci_me,id,temp,size);
         fflush(stdout);
       }
#if 1
       /* delete segment id so that OS cleans it when all attached processes are gone */
       if(shmctl( id, IPC_RMID, (struct shmid_ds *)NULL))
          fprintf(stderr,"failed to remove shm id=%d\n",id);
#endif

    }
    POST_ALLOC_CHECK(temp,sz);

    region_list[alloc_regions].addr = temp;
    region_list[alloc_regions].id = id;
    region_list[alloc_regions].attached=1;
    region_list[alloc_regions].sz=sz;
    alloc_regions++;

    if(DEBUG2_){
      printf("%d:allocate:id=%d addr=%p size=%ld\n",armci_me,id,temp,size);
      fflush(stdout);
    }

    return (void*) (temp);
}
    
/******************** common code for the two versions *********************/


/*\ Allocate a block of shared memory - called by master process
\*/
char *Create_Shared_Region(long *id, long size, long *offset)
{
  char *temp;  
int  reg, refreg=0,nreg;
  
    if(alloc_regions>=MAX_REGIONS)
       armci_die("Create_Shared_Region: to many regions ",0);

    if(DEBUG_){
      printf("%d:CreateSharedRegion %d:size=%ld\n",armci_me,create_call++,size);
      fflush(stdout);
    }

    /*initialization: 1st allocation request */
    if(!alloc_regions){
       for(reg=0;reg<MAX_REGIONS;reg++){
          region_list[reg].addr=(char*)0;
          region_list[reg].attached=0;
          region_list[reg].id=0;
       }
       if(DEBUG_){
          printf("%d:1st CreateSharedRegion: allocation unit:%ldK,shmax:%ldK\n",
                 armci_me,MinShmem,MaxShmem);
          fflush(stdout);
       }

       kr_malloc_init(SHM_UNIT, (size_t)MinShmem, (size_t)MaxShmem, 
		      armci_allocate, 0, &ctx_shmem);
       ctx_shmem.ctx_type = KR_CTX_SHMEM;
       id[SHMIDLEN-2]=MinShmem;
    }

    if(!alloc_regions)  temp = kr_malloc((size_t)size, &ctx_shmem,0, NULL, NULL);
    else temp = kr_malloc((size_t)size, ctx_shmem_global,0,NULL,NULL);

    if(temp == (char*)0 )
       armci_die("CreateSharedRegion:kr_malloc failed KB=",(int)size>>10);
    
    if(!(nreg=find_regions(temp,id,&reg)))
        armci_die("CreateSharedRegion: allocation inconsitent",0);

#ifndef MULTIPLE_REGIONS
    refreg = reg;
#endif

    if(STAMP) *((int*)temp) = alloc_regions-1;
    *offset = (long) (temp - region_list[refreg].addr);
    id[IDLOC]=region_list[reg].sz; /* elan post check */
    occup_blocks++;
  
    if(DEBUG_){ 
      printf("%d:CreateShmReg:reg=%d id=%ld off=%ld ptr=%p adr=%p s=%d n=%d sz=%ld\n",
           armci_me,reg,region_list[reg].id,*offset,region_list[reg].addr,
           temp,(int)size,nreg,id[IDLOC]);
      fflush(stdout);
    }

    return temp;
}



/*\ only process that created shared region returns the pointer to kr_malloc 
\*/
void Free_Shmem_Ptr( id, size, addr)
     long id, size;
     char *addr;
{
  kr_free(addr, ctx_shmem_global);
}


void Delete_All_Regions()
{
int reg;
int code=0;
extern int armci_me;

  for(reg = 0; reg < MAX_REGIONS; reg++){
    if(region_list[reg].addr != (char*)0){
      code += shmctl((int)region_list[reg].id,IPC_RMID,(struct shmid_ds *)NULL);
      region_list[reg].addr = (char*)0;
      region_list[reg].attached = 0;
      if(DEBUG_)
         fprintf(stderr,"%d Delete_All_Regions id=%d code=%d\n",armci_me, 
                (int)region_list[reg].id, code);
    }
  }
}


static int
Semop(int semid,struct sembuf *opers,size_t nops) {
      int ret;

      if((ret = semop(semid,opers,nops)) == -1) {
         printf("%d: semop error performing %i operation(s) on semid %i.\n",armci_me,(int) nops,semid);
         switch(errno) {
           case EFBIG:  printf(" semop errno=EFBIG.\n"); break;
           case E2BIG:  printf(" semop errno=E2BIG.\n"); break;
           case EINTR:  printf(" semop errno=EINTR.\n"); break;
           case EINVAL: printf(" semop errno=EINVAL.\n"); break;
           case EACCES: printf(" semop errno=EACCES.\n"); break;
           case EAGAIN: printf(" semop errno=EAGAIN.\n"); break;
           case ENOSPC: printf(" semop errno=ENOSPC.\n"); break;
           case ERANGE: printf(" semop errno=ERANGE.\n"); break;
           case EFAULT: printf(" semop errno=EFAULT.\n"); break;
           default:
              printf(" unusual semop errno=%d\n",errno); break;
         }
         armci_die("semaphore error",911);
      }

      return ret;
}


int
semaphoreOperation(int semid,int semnum,int semop) {
      struct sembuf op;
      op.sem_op  = semop;
      op.sem_num = semnum;
      op.sem_flg = 0;
      return Semop(semid,&op,1);
}


void
semaphoreAcquire(int semid, int semnum,int access) {
     struct sembuf op;
      op.sem_op = -access;
      op.sem_num = semnum;
      op.sem_flg = 0;
      Semop(semid,&op,1);
}


void
semaphoreRelease(int semid,int semnum,int access) {
      struct sembuf op;
      op.sem_op  = access;
      op.sem_num = semnum;
      op.sem_flg = 0;
      Semop(semid,&op,1);
}

#else
 what are doing here ?
#endif
