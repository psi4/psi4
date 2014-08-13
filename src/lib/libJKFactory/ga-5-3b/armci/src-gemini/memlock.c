#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: memlock.c,v 1.24.2.3 2007-08-29 17:32:32 manoj Exp $ */
#include "armcip.h"
#include "locks.h"
#include "copy.h"
#include "memlock.h"
#include <stdio.h>

#define DEBUG_ 0
#define INVALID_VAL -9999999

#ifdef DATA_SERVER
#  define CORRECT_PTR 
#endif
size_t armci_mem_offset=0;

/* We start by  using table: assign address of local variable set to 1
 * On shmem systems, this addres is overwritten by a shared memory location
 * when memlock array is allocated in armci_init 
 * Therefore, any process within shmem node can reset armci_use_memlock_table
 * to "not used" when offset changes. Since the variable is in shmem, everybody
 * on that SMP node will see the change and use the same locking functions
 */ 
int init_use_memlock_table=1;
int *armci_use_memlock_table=&init_use_memlock_table;

static int locked_slot=INVALID_VAL;

volatile double armci_dummy_work=0.;
void **memlock_table_array;

/* constants for cache line alignment */
#  define CALGN 64
#  define LOG_CALGN 6

#define ALIGN_ADDRESS(x) (char*)((((unsigned long)x) >> LOG_CALGN) << LOG_CALGN) 
static memlock_t table[MAX_SLOTS];


/*\ simple locking scheme that ignores addresses
\*/
void armci_lockmem_(void *pstart, void *pend, int proc)
{

#if defined(CLUSTER) && !defined(SGIALTIX)
    int lock = (proc-armci_clus_info[armci_clus_id(proc)].master)%NUM_LOCKS;
#else
    int lock = 0;
#endif

    if(DEBUG_){
      printf("%d: armci_lockmem_ proc=%d lock=%d\n",armci_me,proc,lock);
      fflush(stdout);
    }

    NATIVE_LOCK(lock,proc);
    if(DEBUG_){
      printf("%d: armci_lockmem_ done\n",armci_me);
      fflush(stdout);
    }
}

void armci_unlockmem_(int proc)
{

#if defined(CLUSTER) && !defined(SGIALTIX) 
    int lock = (proc-armci_clus_info[armci_clus_id(proc)].master)%NUM_LOCKS;
#else
    int lock = 0;
#endif
    if(DEBUG_){
      printf("%d: armci_unlockmem_ proc=%d lock=%d\n",armci_me,proc,lock);
      fflush(stdout);
    }
     NATIVE_UNLOCK(lock,proc);
}



/*\ idle for a time proportional to factor 
\*/
void armci_waitsome(int factor)
{
int i=factor*100000;

   if(factor <= 1) armci_dummy_work =0.;
   if(factor < 1) return;
   while(--i){
      armci_dummy_work = armci_dummy_work + 1./(double)i;  
   }
}
   
/*\ acquire exclusive LOCK to MEMORY area <pstart,pend> owned by process "proc"
 *   . only one area can be locked at a time by the calling process
 *   . must unlock it with armci_unlockmem
\*/
void armci_lockmem(void *start, void *end, int proc)
{
     register void* pstart, *pend;
     register  int slot, avail=0;
     int turn=0, conflict=0;
     memlock_t *memlock_table;
#if defined(CLUSTER) && !defined(SGIALTIX)
    int lock = (proc-armci_clus_info[armci_clus_id(proc)].master)%NUM_LOCKS;
#else
    int lock = 0;
#endif

#ifdef CORRECT_PTR
     if(! *armci_use_memlock_table){
       /* if offset invalid, use dumb locking scheme ignoring addresses */
       armci_lockmem_(start, end, proc); 
       return;
     }

#  ifndef SGIALTIX
     /* when processes are attached to a shmem region at different addresses,
      * addresses written to memlock table must be adjusted to the node master
      */
     if(armci_mem_offset){
        start = armci_mem_offset + (char*)start;
        end   = armci_mem_offset + (char*)end;
     }
#  endif
#endif

     if(DEBUG_){
       printf("%d: calling armci_lockmem for %d range %p -%p\n",
              armci_me, proc, start,end);
       fflush(stdout);
     }
     memlock_table = (memlock_t*)memlock_table_array[proc];


#ifdef ALIGN_ADDRESS
     /* align address range on cache line boundary to avoid false sharing */
     pstart = ALIGN_ADDRESS(start);
     pend = CALGN -1 + ALIGN_ADDRESS(end);
#else
     pstart=start;
     pend =end;
#endif

#ifdef CRAY_SHMEM
     { /* adjust according the remote process raw address */
        long bytes = (long) ((char*)pend-(char*)pstart);
        extern void* armci_shmalloc_remote_addr(void *ptr, int proc);
        pstart = armci_shmalloc_remote_addr(pstart, proc);
        pend   = (char*)pstart + bytes;
     }
#endif
     while(1){
        NATIVE_LOCK(lock,proc);

        armci_get(memlock_table, table, sizeof(table), proc);
/*        armci_copy(memlock_table, table, sizeof(table));*/
        
        /* inspect the table */
        conflict = 0; avail =-1;
        for(slot = 0; slot < MAX_SLOTS; slot ++){

            /* nonzero starting address means the slot is occupied */ 
            if(table[slot].start == NULL){

              /* remember a free slot to store address range */
              avail = slot;  
          
            }else{
              /*check for conflict: overlap between stored and current range*/
              if(  (pstart >= table[slot].start && pstart <= table[slot].end)
                 || (pend >= table[slot].start && pend <= table[slot].end) ){

                  conflict = 1;
                  break;

              }
              /*
              printf("%d: locking %ld-%ld (%d) conflict\n",
                  armci_me,  */
            }
       }
        
       if(avail != -1 && !conflict) break;

       NATIVE_UNLOCK(lock,proc);
       armci_waitsome( ++turn );

     }

     /* we got the memory lock: enter address into the table */
     table[avail].start = pstart;
     table[avail].end = pend;
     armci_put(table+avail,memlock_table+avail,sizeof(memlock_t),proc);

     FENCE_NODE(proc);

     NATIVE_UNLOCK(lock,proc);
     locked_slot = avail;

}
        

/*\ release lock to the memory area locked by previous call to armci_lockemem
\*/
void armci_unlockmem(int proc)
{
     void *null[2] = {NULL,NULL};
     memlock_t *memlock_table;

#ifdef CORRECT_PTR
     if(! *armci_use_memlock_table){
       /* if offset invalid, use dumb locking scheme ignoring addresses */
       armci_unlockmem_(proc);               
       return;
     }
#endif

#ifdef DEBUG
     if(locked_slot == INVALID_VAL) armci_die("armci_unlock: empty",0);
     if(locked_slot >= MAX_SLOTS || locked_slot <0) 
        armci_die("armci_unlock: corrupted slot?",locked_slot);
#endif

     memlock_table = (memlock_t*)memlock_table_array[proc];
     armci_put(null,&memlock_table[locked_slot].start,2*sizeof(void*),proc);

}



/*\ based on address for set by master, determine correction for
 *  memory addresses set in memlock table
 *  if the correction/offset ever changes stop using memlock table locking
\*/ 
void armci_set_mem_offset(void *ptr)
{
   extern size_t armci_mem_offset;
   size_t off;
   static int first_time=1;
   volatile void *ref_ptr;

    ARMCI_PR_DBG("enter",0);
   /* do not care if memlock not used */
   if(! *armci_use_memlock_table) return;

   if(!ptr) armci_die("armci_set_mem_offset : null ptr",0);
   ref_ptr = *(void**)ptr;
   off = (size_t)((char*)ref_ptr - (char*)ptr);

   if(first_time){

      armci_mem_offset =off;
      first_time =0;
      if(DEBUG_){
        printf("%d memlock offset=%ld ref=%p ptr=%p\n",armci_me,
                  (long)armci_mem_offset, ref_ptr, ptr); fflush(stdout);
      }

   }else{
      if(armci_mem_offset != off){
         *armci_use_memlock_table =0;
         fprintf(stderr, "%d: WARNING:armci_set_mem_offset: offset changed %ld to %ld\n",
                 armci_me, (long)armci_mem_offset, (long)off); fflush(stdout);
      }
   }
}
