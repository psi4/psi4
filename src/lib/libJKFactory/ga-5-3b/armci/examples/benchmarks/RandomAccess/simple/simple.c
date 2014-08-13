#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <mpi.h>

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_STAT_H
#   include <sys/stat.h>
#endif
#if HAVE_FCNTL_H
#   include <fcntl.h>
#endif

#include "armci.h"
#include "message.h"
#include "../ra_common.h"

#define DEBUG 0

static int me, nproc;
static u64Int procnumupdates,myglobalstart,globaltablelen,mytablelen,mintablesize,*procglobalstart; 
static int bigtables,Remainder;
static int my_free_handle;
FILE *fd;
u64Int **HPCC_Table;

struct vector_dscr_t{
  int active;
  int dstproc;
  u64Int **xmitbuffer;
  struct vector_dscr_t *next;
  armci_giov_t a_v;
};
struct vector_dscr_t *vec_dscr;
static struct vector_dscr_t curdscr_val;
struct vector_dscr_t *curdscr_ptr = &curdscr_val;

static armci_hdl_t* _get_next_handle()
{
    my_free_handle++;
    my_free_handle %= MAX_OUTSTANDING_HANDLES;
}
/* from hpcc RandomAccess/utility.c */
static u64Int HPCC_starts(s64Int n)
{
int i,j;
u64Int temp,ran,m2[64];

    while (n < 0) n += PERIOD;
    while (n > PERIOD) n -= PERIOD;
    if (n == 0) return 0x1;

    temp = 0x1;
    for (i=0; i<64; i++) {
      m2[i] = temp;
      temp = (temp << 1) ^ ((s64Int) temp < 0 ? POLY : 0);
      temp = (temp << 1) ^ ((s64Int) temp < 0 ? POLY : 0);
    }

    for (i=62; i>=0; i--)
      if ((n >> i) & 1)
        break;

    ran = 0x2;
    while (i > 0) {
      temp = 0;
      for (j=0; j<64; j++)
        if ((ran >> j) & 1)
          temp ^= m2[j];
      ran = temp;
      i -= 1;
      if ((n >> i) & 1)
      ran = (ran << 1) ^ ((s64Int) ran < 0 ? POLY : 0);
    }
    return ran;

}

static void xmitvector(){
struct vector_dscr_t *tmp=curdscr_val.next, *tmp1;
u64Int myscale=0;
    curdscr_ptr = &curdscr_val;
    while(tmp!=NULL){
      if(tmp->dstproc!=me)
        ARMCI_NbAccV(ARMCI_ACC_RA,&myscale,&tmp->a_v,1,tmp->dstproc,NULL);
      tmp=tmp->next;
    }
    tmp=curdscr_val.next;
    while(tmp!=NULL){
      if(tmp->dstproc==me)
        ARMCI_AccV(ARMCI_ACC_RA,&myscale,&tmp->a_v,1,tmp->dstproc);
      tmp=tmp->next;
    }
    tmp=curdscr_val.next;
    while(tmp!=NULL){
      tmp->active=0;
      tmp->a_v.ptr_array_len=0;
      tmp1=tmp->next;
      tmp->next=NULL;
      tmp = tmp1;
    }
}

static void addtovector(u64Int ran){
u64Int offset;
int proc;
struct vector_dscr_t *tmp;
    offset = ran & (globaltablelen-1);
    if(offset < bigtables)
      proc=offset/(mintablesize+1);
    else
      proc=(offset-Remainder)/mintablesize;
    tmp =&vec_dscr[proc]; 
    if(!tmp->active){
      tmp->active = 1;
      tmp->next = NULL;
      curdscr_ptr->next = tmp;
      curdscr_ptr = tmp;
    }
    *(u64Int *)tmp->a_v.src_ptr_array[tmp->a_v.ptr_array_len]=ran;
    tmp->a_v.dst_ptr_array[tmp->a_v.ptr_array_len]=(void *)(HPCC_Table[proc]+(offset-procglobalstart[proc]));
    tmp->a_v.ptr_array_len+=1;
}

void HPCCRandom_Access()
{
int i;
u64Int ran;
    ran = HPCC_starts (4 * myglobalstart);
    for(i=0;i<procnumupdates;i++){
      ran = (ran << 1) ^ ((s64Int)ran < ZERO64B ? POLY : ZERO64B);
      addtovector(ran);
      if(i && i%MAX_TOTAL_PENDING_UPDATES==0)xmitvector();
    }
    if(i%MAX_TOTAL_PENDING_UPDATES)xmitvector();
}

static void initialize_tables()
{
int rc,i,j;
armci_domain_t d;
    curdscr_val.next=NULL;
    curdscr_ptr = &curdscr_val;
    procglobalstart = (u64Int *)calloc(nproc,sizeof(u64Int)); 
    procglobalstart[me]=myglobalstart;
#if   SIZEOF_LONG == 8
    armci_msg_lgop(procglobalstart,nproc,"+");
#elif SIZEOF_LONG_LONG == 8
    armci_msg_llgop(procglobalstart,nproc,"+");
#else
#   error could not determine 64 bit int type
#endif
    HPCC_Table = (u64Int **)malloc(sizeof(u64Int *)*nproc);    
    if(HPCC_Table == NULL)
      ARMCI_Error("initialize_tables:Table pointer malloc failed",(mytablelen*sizeof(u64Int)));
    if(rc=ARMCI_Malloc((void **)HPCC_Table,mytablelen*sizeof(u64Int)))
      ARMCI_Error("initialize_tables:Global Table malloc failed",(mytablelen*sizeof(u64Int)));
    for(i=0;i<mytablelen;i++){
      HPCC_Table[me][i] = i + myglobalstart;
    }
    vec_dscr = (struct vector_dscr_t *)malloc(sizeof(struct vector_dscr_t)*nproc);
    if(vec_dscr == NULL)
      ARMCI_Error("initialize_tables:vec_dscr malloc fail",sizeof(struct vector_dscr_t)*nproc);
    for(i=0;i<nproc;i++){
      vec_dscr[i].next=NULL;
      vec_dscr[i].active=0;
      vec_dscr[i].xmitbuffer = (u64Int **)malloc(sizeof(u64Int *)*nproc);
      if(vec_dscr[i].xmitbuffer == NULL)
         ARMCI_Error("initialize_tables:xmitbuffer malloc failed",sizeof(u64Int *)*nproc);
      if(rc=ARMCI_Malloc((void **)vec_dscr[i].xmitbuffer,MAX_TOTAL_PENDING_UPDATES*sizeof(u64Int))) 
        ARMCI_Error("initialize_tables:xmitbuffer armci_malloc failed",sizeof(u64Int)*MAX_TOTAL_PENDING_UPDATES);
      vec_dscr[i].dstproc = i;
      vec_dscr[i].a_v.src_ptr_array=(void **)malloc(sizeof(void *)*MAX_TOTAL_PENDING_UPDATES);
      for(j=0;j<MAX_TOTAL_PENDING_UPDATES;j++)
        vec_dscr[i].a_v.src_ptr_array[j]=(void *)(vec_dscr[i].xmitbuffer[me]+j);
      
      vec_dscr[i].a_v.dst_ptr_array=(void **)malloc(sizeof(void *)*MAX_TOTAL_PENDING_UPDATES);
      if(vec_dscr[i].a_v.src_ptr_array==NULL || vec_dscr[i].a_v.dst_ptr_array==NULL)
        ARMCI_Error("initialize_tables:.src_ptr_array malloc fail",sizeof(void *)*MAX_TOTAL_PENDING_UPDATES);
      vec_dscr[i].a_v.ptr_array_len = 0;
      vec_dscr[i].a_v.bytes = 8;
    }

}

static void finalize_tables()
{
    ARMCI_Free(HPCC_Table[me]);
    free(HPCC_Table);
}     

int main(argc, argv)
int argc;
char **argv;
{
s64Int i;
int log2nproc;
double CPUTime;        /* CPU  time to update table */
double RealTime;       /* Real time to update table */

double TotalMem;
int PowerofTwo;

u64Int NumUpdates;     /* actual number of updates to table */
s64Int ProcNumUpdates; /* number of updates per processor */

FILE *outFile = NULL;
double *GUPs;
double max_time,min_time,avg_time,time_start,time_stop,total_time;
    armci_msg_init(&argc,&argv);
    nproc = armci_msg_nproc();
    me = armci_msg_me();
    ARMCI_Init();      /* initialize ARMCI */

    if(me==0)printf("\n                          RANDOM ACCESS EXAMPLE\n");
    if(argc<2){
       if(me==0){
         printf(" CORRECT USAGE IS:");
         printf("\n\n <launch commands> simple.x inpfile\n");
         fflush(stdout);
       }
       ARMCI_Finalize();
       armci_msg_finalize();
       return 0;
    }
    globaltablelen = atoi(argv[1]);

    mintablesize = globaltablelen/nproc;
    Remainder = globaltablelen - mintablesize*nproc;
    bigtables = (mintablesize+1)*Remainder;
    if(me<Remainder){
      mytablelen = mintablesize+1;
      myglobalstart = mytablelen*me; 
    }
    else{
      mytablelen = mintablesize;
      myglobalstart = mytablelen*me+Remainder; 
    }
    procnumupdates = 4*mytablelen;
#if DEBUG
    printf("\n%d:%d is totaltable, mintablesize=%d rem=%d big=%d glosta=%d tablelen=%d numup=%d\n",me, globaltablelen, mintablesize,Remainder,bigtables,myglobalstart,mytablelen,procnumupdates);
#endif
    armci_msg_barrier();

    initialize_tables();

    if(me==0)printf("\n\nStarting Random Access....");
    armci_msg_barrier();
    time_start=armci_timer();
    HPCCRandom_Access();
    time_stop=armci_timer();
    armci_msg_barrier();
    total_time=(time_stop-time_start);
    max_time=total_time;
    min_time=total_time;
    avg_time=total_time;
    armci_msg_dgop(&max_time,1,"max");
    armci_msg_dgop(&min_time,1,"min");
    armci_msg_dgop(&avg_time,1,"+");
    avg_time/=nproc;
    if(me==0)printf("\nGUPs = %.9f %.9f %.9f Billion(10^9) Updates/PE  per second [GUP/s]\n",1e-9*procnumupdates/max_time,1e-9*procnumupdates/min_time,1e-9*procnumupdates/avg_time);

    finalize_tables();
    if(me==0)printf("Terminating..\n");
    ARMCI_Finalize();
    armci_msg_finalize();
    return 0;
}
