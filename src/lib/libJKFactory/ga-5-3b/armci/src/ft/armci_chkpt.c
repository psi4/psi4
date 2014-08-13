#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*interfaces for checkpointing */

/* TODO
 * work on the case if pagenum==firstpage or lastpage when writing pages
 */
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_SETJMP_H
#   include <setjmp.h>
#endif
#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_SYSCALL_H
#   include <sys/syscall.h>
#endif
#if HAVE_SYS_MMAN_H
#   include <sys/mman.h>
#endif
#if HAVE_SYS_PARAM_H
#   include <sys/param.h>
#endif
#if HAVE_SYS_WAIT_H
#   include <sys/wait.h>
#endif
#if HAVE_SIGNAL_H
#   include <signal.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif
#if HAVE_FCNTL_H
#   include <fcntl.h>
#endif
#if HAVE_ERRNO_H
#   include <errno.h>
#endif
#if HAVE_DIRENT_H
#   include <dirent.h>
#endif
#if HAVE_STDARG_H
#   include <stdarg.h>
#endif
#include <asm/page.h>
#include "armcip.h"
#include "message.h"
#include "armci_chkpt.h"

#define DEBUG 0
#define DEBUG_ 1
/*\
 * ----------CORE FUNCTIONS -----------
 * armci_init_checkpoint() - the checkpointing code is initialized 
 * armci_icheckpoint_init() - called when with first checkpoint
 * armci_icheckpoint_finalize() - called when done with chkpt
 * armci_icheckpoint() - called every time we checkpoint
 * armci_recover() - called to recoved
 * ----------SUPPORT FUNCTIONS -----------
 * armci_ckpt_pgfh() - pagefault handler, to set the list of dirty pages
 * armci_monitor_addr() - monitors address to look for and pages to set readonly
 * armci_create_record() - create the record for local part of data being stored
 * armci_protect_pages() - to protect pages in the dirty page list
 * armci_storage_read() - reads record from storage 
 * armci_storage_write() - writes into storage
 * armci_storage_fopen() - opens a file in storage
 *
\*/

/*\ global variables
\*/
static armci_storage_record_t armci_storage_record[1001];
static int isheap;
static int number_of_records=1;
static int next_available_rid=0;
static int mypagesize; 
int **armci_rec_ind;
static armci_page_info_t armci_dpage_info;
static int checkpointing_initialized=0;
static int armci_recovering=0;
static char *armci_ckpt_bspBottom=NULL;

/* return the current position (index) of the ckpt file */
#define CURR_FILE_POS(rid)  armci_storage_record[rid].fileinfo.startindex
#define UPDATE_FILE_POS(rid,size) \
    armci_storage_record[rid].fileinfo.startindex += (size)

void ARMCI_Ckpt_create_ds(armci_ckpt_ds_t *ckptds, int count)
{   
    armci_create_ckptds(ckptds,count);
}   

int ARMCI_Ckpt_init(char *filename, ARMCI_Group *grp, int savestack,
        int saveheap, armci_ckpt_ds_t *ckptds)
{     
int rid;
    rid = armci_icheckpoint_init(filename,grp,savestack,saveheap,ckptds);
    return(rid);
}   
      
int ARMCI_Ckpt(int rid)
{   
    return(armci_icheckpoint(rid));
}
    
void ARMCI_Ckpt_Recover(int rid, int iamreplacement)
{     
    armci_irecover(rid, iamreplacement);
}
void ARMCI_Ckpt_finalize(int rid)
{
    armci_icheckpoint_finalize(rid);
}

/* ----------SUPPORT FUNCTIONS ----------- */

/* This function is called from ARMCI sigv and sigbus handler */
int armci_ckpt_pgfh(void *addr, int errno, int fd)
{
    char *paddr;
    unsigned long pagenum;
    /*find the page number and the corresponding page aligned address*/
    pagenum = (unsigned long)((long)addr/mypagesize);
    paddr = (char*)(pagenum*mypagesize); 

    if(DEBUG)printf("%d:paddr=%p addr=%p %lu\n",armci_me,paddr,addr,pagenum);

    /*page is being touched change page permission to READ/WRITE*/
    mprotect(paddr, mypagesize, PROT_READ | PROT_WRITE);
    
    /*mark pagenumber dirty in dirty page array*/
    armci_dpage_info.touched_page_arr[armci_dpage_info.num_touched_pages++] =
            pagenum;
    if(pagenum<armci_dpage_info.firstpage)
            armci_dpage_info.firstpage = pagenum;
    if(pagenum>armci_dpage_info.lastpage)
            armci_dpage_info.lastpage = pagenum;
#if 1
    printf("%d: armci_ckpt_pgfh(): num_touched_pages = %ld\n",armci_me,armci_dpage_info.num_touched_pages);fflush(stdout);
#endif
    return(0);
}

#if 0
printf("%d:pagenum=%d first=%d last=%d\n",armci_me,pagenum,armci_storage_record[i].firstpage,(armci_storage_record[i].firstpage+armci_storage_record[i].totalpages));fflush(stdout);
#endif

static int armci_create_record(ARMCI_Group *group, int count)
{
    int recind;
    int relprocid;
    int rc;

    rc = ARMCI_Group_rank(group,&relprocid);
    /*create and broadcast new index in the records data structure*/
    next_available_rid = 0;
    if(next_available_rid==0){
       if(relprocid==0)
         ARMCI_Rmw(ARMCI_FETCH_AND_ADD,&recind,armci_rec_ind[0],1,0);
    }
    else recind=next_available_rid;
    armci_msg_group_bcast_scope(SCOPE_ALL,&recind,sizeof(int),0,group);

    if(recind>1001) armci_die("create_record, failure",recind);
    armci_storage_record[recind].pid = armci_me;
    armci_storage_record[recind].rid = recind;
    armci_storage_record[recind].rel_pid = relprocid;
    memcpy(&armci_storage_record[recind].group,group,sizeof(ARMCI_Group));
    if(count!=0)
       armci_storage_record[recind].user_addr = (armci_monitor_address_t *)malloc(sizeof(armci_monitor_address_t)*count);
    armci_storage_record[recind].user_addr_count=count;
    
    if(next_available_rid!=0)
       next_available_rid = 0;
    else
       number_of_records++;

    return(recind);
}


static void armci_protect_pages(unsigned long startpagenum,unsigned long numpages)
{
    char *addr;
    addr =(char *)((unsigned long)(startpagenum*mypagesize));
    mprotect(addr, mypagesize*numpages,PROT_READ);
    if(DEBUG)printf("%d:protecting address %p %ld\n",armci_me,addr,mypagesize*numpages);
}

/* CHECK: This is a temporary function - remove later. Then make sure remove
 * the ifdef CHECKPOINT2 in armci_init_checkpoint(). Note this should be
 * called inside main(int argc, char **argv), I guess ...  */
void armci_init_checkpoint2()
{
    printf("%d:in armci init checkpoint2\n",armci_me);fflush(stdout);
#ifdef __ia64
    /* get backing store bottom */
    asm("mov %0=ar.bsp": "=r"(armci_ckpt_bspBottom));
    printf("%d: armci_ckpt_bspBottom=%p\n", armci_me, armci_ckpt_bspBottom);
#endif
}

/*\ ----------CORE FUNCTIONS -----------
\*/
/*called in armci init*/
int armci_init_checkpoint(int spare)
{
    int rc;
    extern void ARMCI_Register_Signal_Handler(int sig, void  (*func)());
    extern void armci_create_ft_group();
    
#ifdef CHECKPOINT2
    printf("%d:in armci init checkpoint\n",armci_me);fflush(stdout);
#ifdef __ia64
    /* get backing store bottom */
    asm("mov %0=ar.bsp": "=r"(armci_ckpt_bspBottom));
    printf("%d: armci_ckpt_bspBottom=%p\n", armci_me, armci_ckpt_bspBottom);
#endif
#endif
    mypagesize = getpagesize();
    if(checkpointing_initialized)return(0);

    /* malloc for record index */
    armci_rec_ind = (int **)malloc(sizeof(int *)*armci_nproc);
    if(armci_me==0){
       rc = ARMCI_Malloc((void **)armci_rec_ind, 2*sizeof(int));
       armci_rec_ind[armci_me][0]=armci_rec_ind[armci_me][1]=1;
    }
    else
       rc = ARMCI_Malloc((void **)armci_rec_ind, 0);
    assert(rc==0);

    ARMCI_Register_Signal_Handler(SIGSEGV,(void *)armci_ckpt_pgfh);
    armci_dpage_info.touched_page_arr = (unsigned long *)malloc(sizeof(unsigned long)*100000);
    armci_dpage_info.num_touched_pages=armci_dpage_info.lastpage=0;
    armci_dpage_info.firstpage = 99999999;
    armci_create_ft_group(spare);

    checkpointing_initialized = 1;
    return(0);
}

void armci_create_ckptds(armci_ckpt_ds_t *ckptds, int count)
{
    printf("%d:in armci_create_ckptds with count=%d\n",armci_me,count);fflush(stdout);
    ckptds->count=count;
    ckptds->ptr_arr=(void **)malloc(sizeof(void *)*(count+1));
    ckptds->sz=(size_t *)malloc(sizeof(size_t)*(count+1));
    ckptds->saveonce=(int *)calloc((count+1),sizeof(int));
    if( ckptds->saveonce==NULL || ckptds->ptr_arr==NULL || ckptds->sz == NULL )
      armci_die("malloc failed in armci_create_ckptds",count);
}

void armci_free_ckptds(armci_ckpt_ds_t *ckptds)
{
    free(ckptds->ptr_arr);
    free(ckptds->sz);
}

static void armci_create_protect_pages(armci_monitor_address_t *addrds, void *ptr, unsigned long bytes,int callprotect)
{
    unsigned long laddr;
    unsigned long totalpages;
    unsigned long j;
    addrds->ptr=(void *)(ptr);
    addrds->bytes = bytes;
    laddr = (unsigned long)(addrds->ptr);
    addrds->firstpage = (unsigned long)((unsigned long)laddr/mypagesize);

    if(laddr%mypagesize ==0){
       totalpages = (int)(bytes/mypagesize);
       if(bytes%mypagesize)totalpages++;
    }
    else {
       int shift;
       shift = mypagesize - laddr%mypagesize;
       if(DEBUG);{
         printf("%d:shift=%d bytes=%ld\n",armci_me,shift,bytes);
         fflush(stdout);
       }
       if(bytes<shift)totalpages=1;
       else{
         totalpages = 1+(bytes-shift)/mypagesize;
         if((bytes-shift)%mypagesize)totalpages++;
       }
    }
    addrds->totalpages = totalpages;
    addrds->num_touched_pages = totalpages;
    addrds->touched_page_arr = malloc(totalpages*sizeof(unsigned long));
    if(addrds->touched_page_arr==NULL)
       armci_die("malloc failed in armci_icheckpoint_init",totalpages);
    addrds->touched_page_arr[0]=addrds->firstpage;
    for(j=1;j<totalpages;j++){
       addrds->touched_page_arr[j]=addrds->touched_page_arr[j-1]+1;
    }
    if(DEBUG);{
       printf("%d: addrds=%p first=%lu total=%lu laddr=%lu (%p)\n",armci_me, addrds, addrds->firstpage, addrds->totalpages, laddr, (void*)laddr);
       fflush(stdout);
    }
    if(callprotect){
       if(isheap)
         armci_protect_pages(addrds->firstpage+16,addrds->totalpages);
       else
         armci_protect_pages(addrds->firstpage,addrds->totalpages);
    }
}

/*called everytime a new checkpoint record is created*/
int armci_icheckpoint_init(char *filename,ARMCI_Group *grp, int savestack, 
                int saveheap, armci_ckpt_ds_t *ckptds)
{
    int rid;
    int i=0;
    char line[256],tmpfilename[100];
    unsigned long databottom;
    FILE *fp;
    printf("%d:in armci ckpt init\n",armci_me);fflush(stdout);
    if(DEBUG && ckptds!=NULL)
            printf("%d:ckptdscount=%d\n",armci_me,ckptds->count);

    /*create the record*/
    if(ckptds!=NULL) rid = armci_create_record(grp,ckptds->count);
    else rid = armci_create_record(grp,0);
    if(DEBUG) printf("%d:got rid = %d\n",armci_me,rid);fflush(stdout);

    armci_storage_record[rid].ckpt_heap = saveheap;
    armci_storage_record[rid].ckpt_stack = savestack;

    /*open the file for reading and writing*/
    if(filename == NULL){
      filename = (char *)malloc(sizeof(char)*(11+1+6+1+4));
      if(filename==NULL)armci_die("alloc for filename failed",11+1+6+1+4);
      sprintf(filename,"%s","armci_chkpt_");
      sprintf((filename+strlen(filename)),"%d",armci_me);
      sprintf((filename+strlen(filename)),"%s","_");
      sprintf(filename+strlen(filename),"%d",rid);
    }
    armci_storage_record[rid].fileinfo.filename = malloc(sizeof(char)*strlen(filename));
    if(NULL==armci_storage_record[rid].fileinfo.filename)
      armci_die("malloc failed for filename in ga_icheckpoint_init",0);
    strcpy(armci_storage_record[rid].fileinfo.filename,filename);
    armci_storage_record[rid].fileinfo.fd = armci_storage_fopen(filename);
    armci_storage_record[rid].fileinfo.startindex = 0; /* initialize */
    if(DEBUG){printf("filename=%s\n",filename);fflush(stdout);}

    sprintf(tmpfilename, "/proc/%d/maps",getpid());
    fp=fopen(tmpfilename, "r");
    if (fp == NULL){
      armci_die("couldnt find dseg base",getpid());
    }

    if(armci_storage_record[rid].ckpt_stack){
       unsigned long stackbottom;
       char *start,*end,*tmp;
       armci_monitor_address_t *addrds =&armci_storage_record[rid].stack_mon;
       do {
         fgets(line, 255, fp);
       } while ( line[0] != '6' );
       sscanf(line, "%p", (void**)&databottom);
       printf("%p databot\n",(char *)databottom);
#   ifdef KERNEL_2_4       
       do {
         start = fgets(line, 255, fp);
       } while(start != NULL);
#   else /* KERNEL_2_6 */
       {
         char tmpline[256];
         do {
           strncpy(line, tmpline, 256);
           start = fgets(tmpline, 255, fp);
         } while(tmpline[0] == '6');
       }
#   endif       
       start = strstr(line, "-") + 1;
       if(DEBUG);{printf("stack top=%s\n",start);fflush(stdout);}
       end   = strstr(line, " ");
       *end = 0;
       sscanf(start, "%p", (void**)&stackbottom);
       addrds->ptr = (void *)(stackbottom);
    }
    if(armci_storage_record[rid].ckpt_heap){
       char *datatop;
       armci_monitor_address_t *addrds =&armci_storage_record[rid].heap_mon;
       datatop=(char *)sbrk(0);
       if(!armci_storage_record[rid].ckpt_stack){
         do {
           fgets(line, 255, fp);
         } while ( line[0] != '6' );
         sscanf(line, "%p", (void**)&databottom);
       }
       printf("%d:databot=%p datatop=%p %ld %p %p\n",armci_me,(void*)databottom,datatop,(unsigned long)((char *)datatop-(char *)databottom),&armci_storage_record[rid].jmp, &addrds);
       isheap = 1;
       printf("I'm here 1\n"); fflush(stdout);
       armci_create_protect_pages(addrds,(void *)databottom,(unsigned long)((char *)datatop-(char *)databottom),1);
       printf("I'm here 1a\n"); fflush(stdout);
       isheap = 0;
       mprotect(armci_storage_record,sizeof(armci_storage_record_t)*1000, PROT_READ | PROT_WRITE);
       mprotect(&armci_recovering,sizeof(int), PROT_READ | PROT_WRITE);
       mprotect(&armci_dpage_info,sizeof(armci_page_info_t), PROT_READ | PROT_WRITE);
       mprotect(armci_dpage_info.touched_page_arr,sizeof(unsigned long)*99999999, PROT_READ | PROT_WRITE);
       printf("I'm here 2\n"); fflush(stdout);
    }
    else {
       if(ckptds!=NULL)
          for(i=0;i<ckptds->count;i++){
             armci_monitor_address_t *addrds =&armci_storage_record[rid].user_addr[i];
             addrds->saveonce = ckptds->saveonce[i];
             if(addrds->saveonce)
                armci_create_protect_pages(addrds,ckptds->ptr_arr[i],ckptds->sz[i],0);
             else
                armci_create_protect_pages(addrds,ckptds->ptr_arr[i],ckptds->sz[i],1);
          }
    }
    if(DEBUG){printf("%d:completed init\n",armci_me);fflush(stdout);}
    return(rid);
}

int armci_create_touchedpagearray(unsigned long *tpa,unsigned long firstpage, unsigned long totalpages){
unsigned long i,j=0;
    printf("In armci_create_touchedpagearray(): armci_dpage_info.num_touched_page=%ld\n", armci_dpage_info.num_touched_pages); 
    tpa = (unsigned long *)malloc(sizeof(unsigned long)*armci_dpage_info.num_touched_pages);
    for(i=0;i<armci_dpage_info.num_touched_pages;i++){
       if(DEBUG){
         printf("%d: armci_create_tpa %lu %lu %lu\n",armci_me,
                         armci_dpage_info.touched_page_arr[i],firstpage,
                         (firstpage+totalpages));
         fflush(stdout);
       }
       if(armci_dpage_info.touched_page_arr[i]>=firstpage || armci_dpage_info.touched_page_arr[i] < (firstpage+totalpages)){
         tpa[j]=armci_dpage_info.touched_page_arr[i];
         j++;
       }
    }

    armci_dpage_info.num_touched_pages=0;
    return(j);

}

/* CHECK: put all the defines here in the armci_chpt.h header or in some
 * meaningful header*/
/* returns the stack head address. SP register contains this address. We save
   some EXTRA_SPACE as there may be some space used below the current SP */
#define EXTRA_STACK_SPACE   (1 << 11)
#define TMP_STACK_SIZE      (1 << 12)  /* CHECK: 4K is enough, I guess */
#ifndef JB_SP
#define JB_SP 0
#endif

#define EST_OFFSET          100
#define ARMCI_STACK_VERIFY  1234

#define STACK_TOP ((unsigned long)((unsigned long)(&dummy_first) - EST_OFFSET))
#if 0
static char* stack_head_addr() 
{
    jmp_buf tmp_j;
    setjmp(tmp_j);
    return ((char*) ((tmp_j->__jmpbuf[JB_SP] - EXTRA_STACK_SPACE) & (~PAGE_SIZE)) );
}

/* to get the stack pointer */
unsigned long get_esp()
{
    __asm__(" movl %esp,%eax ");
}
#endif

void what_is_going_on() {
    int a;
    printf("what_is_going_on(): a=%p\n", &a);
}

/*
  In IA64, there is a seperate stack called register stack engine (RSE,
  contains 96 registers) to manage across functions calls (e.g. these
  registers stores the function return address, etc.. In order to save this
  info, flush the stack registers to backing store and save backing store.
  NOTE: backing store is a cache to register stack. 
 */
#if defined(__ia64)
static void armci_ckpt_write_backstore(int rid) 
{
    char *bspTop; /* in IA64 only, back store pointer (bsp) */
    off_t ofs;
    
    /* flush the register stack */
    asm("flushrs");
    
    /* getting back store pointer (bsp). BSP is similar to stack pointer,
     * which points to the top of backing store */
    asm("mov %0=ar.bsp": "=r"(bspTop));
    printf("BSP Pointer (ar.bsp) = %p\n", bspTop);

    armci_storage_record[rid].bsp_mon.ptr = armci_ckpt_bspBottom;
    armci_storage_record[rid].bsp_mon.bytes=((unsigned long)(bspTop) - (unsigned long)(armci_ckpt_bspBottom));

    ofs=CURR_FILE_POS(rid);
    UPDATE_FILE_POS(rid, armci_storage_record[rid].bsp_mon.bytes);
    armci_storage_record[rid].bsp_mon.fileoffset = ofs;
    
    printf("%d: Save Backing store: %p to %p (bytes=%ld: off=%ld)\n\n",armci_me, armci_ckpt_bspBottom, armci_ckpt_bspBottom+armci_storage_record[rid].bsp_mon.bytes, armci_storage_record[rid].bsp_mon.bytes, ofs);fflush(stdout);

    armci_storage_write_ptr(armci_storage_record[rid].fileinfo.fd,
                            armci_storage_record[rid].bsp_mon.ptr,
                            armci_storage_record[rid].bsp_mon.bytes,
                            armci_storage_record[rid].bsp_mon.fileoffset);
    
}
#endif

static void armci_ckpt_write_stack(int rid)
{
    int dummy_first=ARMCI_STACK_VERIFY;
    char *top=NULL;
    off_t ofs = 0;
    int dummy_last=ARMCI_STACK_VERIFY;

    /* top = stack_head_addr(); */
    top = (char*)STACK_TOP;
    
    armci_storage_record[rid].stack_mon.bytes=(unsigned long)(armci_storage_record[rid].stack_mon.ptr)-((unsigned long)(top));

    ofs=CURR_FILE_POS(rid);
    UPDATE_FILE_POS(rid, armci_storage_record[rid].stack_mon.bytes);
    armci_storage_record[rid].stack_mon.fileoffset = ofs;

    printf("%d: ptr = %p\n",armci_me,armci_storage_record[rid].stack_mon.ptr);
    
    printf("%d: Save stack: %p to %p (bytes=%ld : off=%ld)\n\n",armci_me, top, top+armci_storage_record[rid].stack_mon.bytes, armci_storage_record[rid].stack_mon.bytes, ofs);fflush(stdout);
    armci_storage_write_ptr(armci_storage_record[rid].fileinfo.fd,top,
                            armci_storage_record[rid].stack_mon.bytes,
                            armci_storage_record[rid].stack_mon.fileoffset);

#if defined(__ia64)
    /* In IA64, write Backing Store, as it is the cache for Stack Registers */
    armci_ckpt_write_backstore(rid);    
#endif    
}

static void armci_ckpt_write_heap(int rid)
{
    int j;
    off_t ofs = 0;
    char *addr=NULL;
    
    armci_monitor_address_t *addrds =&armci_storage_record[rid].heap_mon;
    addr = sbrk(0);
    ofs=CURR_FILE_POS(rid); 
    /* UPDATE_FILE_POS(rid, armci_storage_record[rid].heap_mon.bytes); */
    armci_storage_record[rid].heap_mon.fileoffset = ofs;
    
    if(addr > (char *)armci_storage_record[rid].heap_mon.ptr){ 
       /*this means change in data segment - save what ever is new and reset size*/
       void *tmpaddr = addrds->ptr;
       void *tmpaddr1 = (void *)((char *)addrds->ptr+addrds->bytes);
       unsigned long firstpage = addrds->firstpage;
       unsigned long totalpages = addrds->totalpages;
       /*first save new pages*/
       ofs=(off_t)(armci_storage_record[rid].heap_mon.fileoffset+totalpages*mypagesize);
       /*problem here - remove mallocs*/
       isheap = 1;
       armci_create_protect_pages(addrds,tmpaddr1,(addr-(char *)tmpaddr1),1);
       isheap = 0;
       printf("%d: Hello 1\n", armci_me);
       armci_storage_write_pages(armci_storage_record[rid].fileinfo.fd,addrds->firstpage,addrds->touched_page_arr,addrds->num_touched_pages,mypagesize,ofs);
       /*now write the touched pages*/
       addrds->ptr = tmpaddr;
       addrds->bytes = addr-(char *)tmpaddr;
       addrds->firstpage = firstpage;
       /*problem here if last data seg addr is not a page boundary*/
       addrds->totalpages+=totalpages;
    }
    
    /*write touched pages since*/
    ofs=(off_t)(armci_storage_record[rid].heap_mon.fileoffset);
    addrds->num_touched_pages = armci_create_touchedpagearray(addrds->touched_page_arr,addrds->firstpage,addrds->totalpages);
    printf("%d: Hello 2\n", armci_me);
    /*problem here - remove mallocs*/
    if(addrds->num_touched_pages!=0)
       addrds->num_touched_pages = armci_create_touchedpagearray(addrds->touched_page_arr,addrds->firstpage,addrds->totalpages);
    printf("%d: Hello 3: %ld %p %ld %d %ld\n", armci_me,addrds->firstpage,addrds->touched_page_arr,addrds->num_touched_pages,mypagesize,ofs);
    armci_storage_write_pages(armci_storage_record[rid].fileinfo.fd,addrds->firstpage,addrds->touched_page_arr,addrds->num_touched_pages,mypagesize,ofs);
    printf("%d: Hello 3a\n", armci_me);
    for(j=0;j<addrds->num_touched_pages;j++){
       addr =(char *)(addrds->touched_page_arr[j]*mypagesize);
       mprotect(addr, mypagesize,PROT_READ);
    }
    printf("%d: Hello 4\n", armci_me);
}

static void armci_ckpt_write_data(int rid)
{
    int i,j;
    off_t ofs=0;
    char *addr=NULL;

    printf("%d: armci_ckpt_write_data(): Saving data ...\n", armci_me);

    for(i=0;i<armci_storage_record[rid].user_addr_count;i++){
       armci_monitor_address_t *addrds = &armci_storage_record[rid].user_addr[i];
       /* get the file offset */
       ofs = CURR_FILE_POS(rid);
       addrds->fileoffset = ofs;
       
       /*if(addrds->num_touched_pages!=0) CHECK this... */
          addrds->num_touched_pages = armci_create_touchedpagearray(addrds->touched_page_arr,addrds->firstpage,addrds->totalpages);
       
       printf("%d: DATA:[i=%d] addrds=%p off=%ld size=%ld (#ofPages=%ld total=%ld pagesize=%d)\n",
              armci_me, i, addrds, addrds->fileoffset,
              addrds->num_touched_pages*mypagesize, addrds->num_touched_pages,
              addrds->totalpages, mypagesize);
       
       armci_storage_write_pages(armci_storage_record[rid].fileinfo.fd,addrds->firstpage,addrds->touched_page_arr,addrds->num_touched_pages,mypagesize,ofs);
       
       for(j=0;j<addrds->num_touched_pages;j++){
          addr =(char *)(addrds->touched_page_arr[j]*mypagesize);
          mprotect(addr, mypagesize,PROT_READ);
       }

       /* CHECK: should I uncomment this ? check it out!
       bzero(addrds->touched_page_arr,
             sizeof(unsigned long)*addrds->num_touched_pages);
       addrds->num_touched_pages = 0;
       */
       
       /* update the file offset */
       UPDATE_FILE_POS(rid, addrds->bytes);
    }
}

static void armci_ckpt_write_rid(int rid)
{
    
}

/*get the list of changed pages from touched_page_array and rewrite the 
 * changed pages*/
int armci_icheckpoint(int rid)
{
    int rc=ARMCI_CKPT;
    off_t ofs;

    if(DEBUG);{
       printf("%d: in checkpoint rid=%d %p\n",armci_me,rid,&armci_recovering);fflush(stdout);
    }

    if(armci_storage_record[rid].ckpt_stack ||
       armci_storage_record[rid].ckpt_heap) {

#if defined(__ia64)
       {
          char *tmp_bsp;
          /* flush the register stack */
          asm("flushrs");
          /* get the top of backing store */
          asm("mov %0=ar.bsp": "=r"(tmp_bsp));
          printf("tmp: ar.bsp = %p\n", tmp_bsp);
       }
#endif
       if((armci_recovering=setjmp(armci_storage_record[rid].jmp))==0){

          /* 1. file offsets */
          armci_storage_record[rid].fileinfo.startindex = 0;
          ofs=CURR_FILE_POS(rid); UPDATE_FILE_POS(rid,sizeof(jmp_buf));
    
          /* 1a. save jmp buffer env */
          printf("%d: Save jmp_buf: %p to %p (bytes=%ld : off=%ld)\n\n",armci_me, &armci_storage_record[rid].jmp, &armci_storage_record[rid].jmp+armci_storage_record[rid].stack_mon.bytes, sizeof(jmp_buf), ofs);fflush(stdout);
          armci_storage_write_ptr(armci_storage_record[rid].fileinfo.fd,
                                  &armci_storage_record[rid].jmp,
                                  sizeof(jmp_buf), ofs);
          
          /* 2. save stack */
          printf("%d: Saving stack\n", armci_me);
          if(armci_storage_record[rid].ckpt_stack){
             armci_ckpt_write_stack(rid);
          }
          
          /* 3. save data segment (entire heap (or) user specified data ) */
          if(armci_storage_record[rid].ckpt_heap) armci_ckpt_write_heap(rid);
          else armci_ckpt_write_data(rid);
          
          /* 4. CHECK: sync file system, thus data is flushed to disk */
          /* armci_storage_fsync(armci_storage_record[rid].fileinfo.fd); */
          
          /* 5. TODO: save the record index in the file.
             Caution: there are mallocs in the structure. beware. */
          armci_ckpt_write_rid(rid);
          
       }
       else { /*long jump brings us here */

          /* CHECK: open the ckpt files*/
          
          printf("long jump brought us here. Performed Recovery..\n");
          printf("address(rid)=%p address(rc)=%p\n", &rid, &rc);
          what_is_going_on();
          printf("rid=%d address(rid)=%p\n", rid, &rid);
          rc = ARMCI_RESTART;
       }
    }
    else{
       armci_ckpt_write_data(rid);      
    }
         
    armci_msg_group_barrier(&armci_storage_record[rid].group);
    printf("%d: After Barrier\n", armci_me);
    return(rc);
}

/**
 * Recover Backing Store.
 */
#if defined(__ia64)
static void armci_recover_backstore(int rid) 
{
    off_t offset = armci_storage_record[rid].bsp_mon.fileoffset;
    size_t size  = armci_storage_record[rid].bsp_mon.bytes;
    char *bspTop = (char*)((unsigned long)(armci_storage_record[rid].bsp_mon.ptr) + size);
    char *bsp;
    
    asm("flushrs");
    asm("mov %0=ar.bsp": "=r"(bsp));

    /* CHECK: expand the backing store so that the current backing store is
       replaced with saved backing store (CHECK: register stack can be as
       large as 96 registers, so 96*8 bytes) */
    if( (unsigned long)bsp < (unsigned long)(bspTop + 96*8 + EST_OFFSET) ) {
       armci_recover_backstore(rid);
    }
    else{
       printf("%d: armci_recover_backstore(): size=%ld offset=%ld backing store: %p to %p\n", armci_me, size, offset, bspTop, armci_storage_record[rid].bsp_mon.ptr);
       armci_storage_read_ptr(armci_storage_record[rid].fileinfo.fd, bspTop, size, offset);

       printf("%d: armci_recover_backstore(): rid=%d\n", armci_me, rid);
       
       /* CHECK: Is there a way to verify backing store recovery
          (similar to stack) */
    }
    /**
     * CHECK: Do nothing here. Recursive function in action.
     */
}
#endif


/**
 * Recover stack: restore a saved stack by overwriting the current stack
 * of this process . The idea of restoring the stack is, we are going to
 * replace the contents of current stack, so that longjmp is legitimate.
 */
#if 0
static void armci_recover_stack(int rid) 
{   
    off_t offset = sizeof(jmp_buf)+4*sizeof(int);
    size_t size  = armci_storage_record[rid].stack_mon.bytes;
    char *stacktop = (char*)((unsigned long)(armci_storage_record[rid].stack_mon.ptr) - size);
    int dummy;
    printf("check=%p %p; rid=%d\n", &dummy, &offset, rid);

    /* CHECK: check whether current stack frame is above the old (saved)
       stack. If so, the recover the stack, else call thus recursively
       until the current stack is above the old stack */
    if( (unsigned long)&dummy >= (unsigned long)(stacktop-EST_OFFSET) ) {
       armci_recover_stack(rid);
    }
    else {
       printf("%d: armci_recover_stack(): size=%ld offset=%ld stack: %p to %p\n", armci_me, size, offset, stacktop, armci_storage_record[rid].stack_mon.ptr);
       armci_storage_read_ptr(armci_storage_record[rid].fileinfo.fd, stacktop, size, offset);
       
       { /* verify stack recovery */
          int dummy = *((int*)(stacktop+EST_OFFSET));
          if(dummy != ARMCI_STACK_VERIFY) {
             printf("WARNING: armci_recover_stack FAILED: %d", dummy);
             armci_die("armci_recover_stack FAILED", dummy);
          }
          else if(DEBUG_)
             printf("%d: armci_recover_stack SUCCESS (%d)\n", armci_me, dummy);
       }
       
#ifdef __ia64
       /* recover the backing store (BSP) */
       armci_recover_backstore(rid);
#endif

    }
    /**
     * CHECK: Do nothing here...recursive function in action here..
     */
}
#endif

static void armci_recover_memory(int rid) 
{   
    int dummy;
    off_t ofs;
    size_t stacksize  = armci_storage_record[rid].stack_mon.bytes;
    char *stacktop = (char*)((unsigned long)(armci_storage_record[rid].stack_mon.ptr) - stacksize);

#ifdef __ia64
    size_t bspsize = armci_storage_record[rid].bsp_mon.bytes;
    char *bspTop   = (char*)((unsigned long)(armci_storage_record[rid].bsp_mon.ptr) + bspsize);
    char *bsp;
#endif
    
    printf("armci_recover_stack(): check=%p ; rid=%d\n", &dummy, rid);
    /* call recursively until current stack is above saved stack */
    if( (unsigned long)&dummy >= (unsigned long)(stacktop-EST_OFFSET) )
       armci_recover_memory(rid);
    
#ifdef __ia64
    asm("flushrs");
    asm("mov %0=ar.bsp": "=r"(bsp));

    printf("armci_recover_bsp(): check=%p ; rid=%d\n", &dummy, rid);
    /* similarly, call recursively until current backing store expands
       (register stack can be as large as 96 registers) */
    if( (unsigned long)bsp < (unsigned long)(bspTop + 97*8) )
       armci_recover_memory(rid);
#endif
    
    /* ------------------ recover stack segment ------------------- */
    printf("%d: armci_recover_stack(): fp=%p size=%ld off=%ld stack: %p to %p\n", armci_me, armci_storage_record[rid].fileinfo.fd, stacksize, armci_storage_record[rid].stack_mon.fileoffset, stacktop, armci_storage_record[rid].stack_mon.ptr);
    armci_storage_read_ptr(armci_storage_record[rid].fileinfo.fd, stacktop, stacksize, armci_storage_record[rid].stack_mon.fileoffset);
    { /* verify stack recovery */
       int dummy = *((int*)(stacktop+EST_OFFSET));
       if(dummy != ARMCI_STACK_VERIFY) {
          printf("WARNING: armci_recover_stack FAILED: %d", dummy);
          armci_die("armci_recover_stack FAILED", dummy);
       }
       else if(DEBUG_)
          printf("%d: armci_recover_stack SUCCESS (%d)\n", armci_me, dummy);
    }
    
    /* -------- recover register stack (RSE) segment (IA64 only) -------- */
#ifdef __ia64
    {
       size_t bspsize = armci_storage_record[rid].bsp_mon.bytes;
       char *bspTop   = (char*)((unsigned long)(armci_storage_record[rid].bsp_mon.ptr) + bspsize);

       bsp = (char*)armci_storage_record[rid].bsp_mon.ptr; /* CHECK: */
       
       printf("%d: armci_recover_backstore(): size=%ld off=%ld backing store: %p to %p\n", armci_me, bspsize, armci_storage_record[rid].bsp_mon.fileoffset, armci_storage_record[rid].bsp_mon.ptr, bspTop);
       armci_storage_read_ptr(armci_storage_record[rid].fileinfo.fd, bsp, bspsize, armci_storage_record[rid].bsp_mon.fileoffset);
       printf("%d: armci_recover_backstore(): rid=%d\n", armci_me, rid);
       
       /* CHECK: Is there a way to verify backing store recovery
          (similar to stack) */
    }
#endif

    ofs=0; /* jmp_buf is the first one to be stored in ckpt file, so ofs=0 */
    printf("%d: armci_recover jmp_buf(): size=%ld off=%ld (%p to %p)\n", armci_me, sizeof(jmp_buf), ofs, &armci_storage_record[rid].jmp, (char*)(&armci_storage_record[rid].jmp)+sizeof(jmp_buf));
    armci_storage_read_ptr(armci_storage_record[rid].fileinfo.fd, &armci_storage_record[rid].jmp, sizeof(jmp_buf), ofs);
    armci_msg_group_barrier(&armci_storage_record[rid].group);
    printf("%d: restoring original stack starts @ %p\n", armci_me, (void*)armci_storage_record[rid].jmp->__jmpbuf[JB_SP]);

    longjmp(armci_storage_record[rid].jmp,1);/*goto the restored stack*/
}

static void armci_recover_stack(int rid, char *filename) 
{
    armci_recover_memory(rid);
}

static void armci_recover_heap(int rid, char *filename) 
{
    /* TODO: restore heap */
}

static void armci_recover_data(int rid, char *filename) 
{
    int i,j;
    printf("%d: armci_recover_data(): rid=%d\n", armci_me, rid);
    
    for(i=0;i<armci_storage_record[rid].user_addr_count;i++){
       armci_monitor_address_t *addrds =
         &armci_storage_record[rid].user_addr[i];

       /* CHECK: beware of malloc'ed memory in addrds */
       printf("%d: armci_recover_data(): [i=%d] from=%p size=%ld off=%ld\n",
              armci_me, i, addrds->ptr, addrds->bytes, addrds->fileoffset);
       armci_storage_read_pages(armci_storage_record[rid].fileinfo.fd,
                                addrds->firstpage, addrds->touched_page_arr,
                                addrds->num_touched_pages, mypagesize,
                                addrds->fileoffset);
    }
}

int armci_irecover_NEW(int rid,int iamreplacement)
{
    /* CHECK: do you need to fille rest of fileinfo like name, startindex? */
    char *filename = armci_storage_record[rid].fileinfo.filename;
    sprintf(filename,"%s","armci_chkpt_");
    sprintf((filename+strlen(filename)),"%d",armci_me);
    sprintf((filename+strlen(filename)),"%s","_");
    sprintf(filename+strlen(filename),"%d",rid);
    
    printf("\n%d: Starting recovery...\n", armci_me);
    printf("%d: filename = %s\n", armci_me, filename);

    if(armci_storage_record[rid].ckpt_heap) {
       armci_recover_heap(rid,filename);
    } else {
      armci_recover_data(rid, filename);
    }

    /* stack should be the last thing recovered, as it calls longjmp() */
    if(armci_storage_record[rid].ckpt_stack) armci_recover_stack(rid,filename);

    /*we should never come here things are hosed */
    armci_die("recovery hosed",0);
    return(1);   
}

int armci_irecover(int rid,int iamreplacement) 
{
    armci_irecover_NEW(rid, iamreplacement);
    return 1;
}

#if 0
static int tmpStack[TMP_STACK_SIZE];
int armci_irecover_OLD(int rid,int iamreplacement)
{
    int rc;
    jmp_buf jmp;
    
    /* Save "rid" and "iamreplacement" in a global variable as we are going
       to replace the contents of the current stack. */
    RID = rid; /* CHECK: save rid in a file or somewhere instead of
                * global variable*/
    tmp_iamreplacement = iamreplacement;
    
#if 0
    /* create a temporary stack */
    rc = _setjmp(jmp);
    
    if (rc == 0) {
       /* Goto a temporary stack as we still running on the original stack. To
          do this, update Stack Pointer (SP) to be in a temp stack area. */
       jmp->__jmpbuf[JB_SP] = ((long)((char *)(tmpStack + TMP_STACK_SIZE) - EXTRA_STACK_SPACE) & ~0xf);
       printf("%d: temporary stack starts @ %p\n", armci_me, jmp->__jmpbuf[JB_SP]);
       
       /* CHECK: make this TMP_STACK_SIZE dynamic, by measuring the size of
          the stack from file */
       
       /*
        * Jump back ...
        * But with new 'jmp'
        */
       _longjmp(jmp, 1);
    }
    else
#endif
    {

       /**
        * Now we are on temporary stack. So it is safe to recover stack.
        */
       armci_recover_stack(RID);
       
       /**
       * go to the restored stack by calling longjmp(). Read jmpbuf from file 
       */
       if(tmp_iamreplacement){ /* CHECK: what is iamreplacement */
          rc=armci_storage_read_ptr(armci_storage_record[RID].fileinfo.fd,&armci_storage_record[RID].jmp,sizeof(jmp_buf),4*sizeof(int));
       }
       armci_msg_group_barrier(&armci_storage_record[RID].group);
       printf("%d: restoring original stack starts @ %p\n", armci_me, armci_storage_record[RID].jmp->__jmpbuf[JB_SP]);
       longjmp(armci_storage_record[RID].jmp,1);/*goto the restored stack*/
    }

    /*we should never come here things are hosed */
    armci_die2("recovery hosed",RID,iamreplacement);
    return(1);   
}
#endif

void armci_icheckpoint_finalize(int rid)
{
    int i;
    armci_msg_group_barrier(&armci_storage_record[rid].group);
    for(i=0;i<armci_storage_record[rid].user_addr_count;i++){
       armci_monitor_address_t *addrds=&armci_storage_record[rid].user_addr[i];
       free(addrds->touched_page_arr);
    }
    free(armci_storage_record[rid].user_addr);
    free(armci_storage_record[rid].fileinfo.filename);
    armci_storage_fclose(armci_storage_record[rid].fileinfo.fd);
    next_available_rid = rid;
}

/*
  TODO:
  1. organize all the $ifdef __ia64's properly..They are scattered all
  over and it is difficult to track down and potentially buggy
  2. checkpoint shared memory and mmap regions
  3. I/O file open/close, signals and other system specific stuff ???
  4. memory leaks due to malloc()....free'em
*/
