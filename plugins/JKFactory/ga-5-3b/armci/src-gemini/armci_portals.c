#if HAVE_CONFIG_H
#   include "config.h"
#endif

        /*$id:$*/
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include "armcip.h"
#include "message.h"
#include <assert.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <limits.h>
#include <unistd.h>

#define DEBUG_COMM 0
#define DEBUG_INIT 0
#define DEBUG_SERV 0
#define PUT_LOCAL_ONLY_COMPLETION__

typedef struct arminfo{
  caddr_t ptr[MAX_DS];
  size_t size[MAX_DS];
  long serv_offs[MAX_DS];
  int cur_ds;
}rm_info_t;

static rm_info_t *all_meminfo;

static int client_md_count=0,serv_md_count=0;

typedef struct arns{
  long data;
  long data1;
  struct arns *next; 
} arnode;
 
#ifdef ARMCI_CHECK_STATE
arnode * arlist_add(arnode **p, long i,long j)
{
  arnode *n = (arnode *)malloc(sizeof(arnode));
  if(n == NULL)
    return NULL;
  n->next = *p; 
  *p = n; 
  n->data = i;
  n->data1 = j;
  return *p;
}
 
void arlist_remove(arnode **p)
{
  if(*p != NULL){
    arnode *n = *p;
    *p = (*p)->next;
    free(n);
  }
}
 
arnode **arlist_search(arnode **n, long i)
{
  while (*n != NULL){
    if ((*n)->data == i){
      return n;
    }
    n = &(*n)->next;
  }
  return NULL;
}
 
void arlist_print(arnode *n)
{
  if (n == NULL){
    /*printf("arlist is empty\n");*/
  }
  while (n != NULL){
    printf("%d:%d %d next=%d\n", armci_me,n->data,n->data1,(n->next==NULL)?0:1);
    n = n->next;
  }
}
#endif 

extern void armci_util_wait_int(volatile int *, int , int );
extern void armci_util_wait_long(volatile long *, long, int );

int _armci_portals_server_ready=0;
int _armci_portals_client_ready=0;
int _armci_server_mutex_ready=0;
void *_armci_server_mutex_ptr = NULL;

#ifdef ARMCI_REGISTER_SHMEM
typedef struct {
       void *base_ptr;
       void *serv_ptr;
       size_t size;
       int islocal;
       int valid;
}aptl_reginfo_t;

typedef struct {
       aptl_reginfo_t reginfo[MAX_MEM_REGIONS];  
       int reg_count;
} rem_meminfo_t;
#endif

typedef struct serv_buf_t{
  ptl_handle_md_t md_h;
  ptl_handle_me_t me_h;
  ptl_md_t        md;
  char            *buf;
  char            *bufend;
} serv_buf_t;

char **client_buf_ptrs;
static int armci_server_terminating=0;

serv_buf_t *serv_bufs;
long servackval=ARMCI_STAMP,*serv_ack_ptr=&servackval;
ptl_handle_md_t serv_ack_md_h,serv_response_md_h;

static armci_portals_proc_t _armci_portals_proc_struct;
static armci_portals_serv_t _armci_portals_serv_struct;
static armci_portals_proc_t *portals = &_armci_portals_proc_struct;
static armci_portals_serv_t *serv_portals = &_armci_portals_serv_struct;
/*static */comp_desc _compdesc_array[NUM_COMP_DSCR];

static arnode *arn = NULL;

#ifdef ARMCI_REGISTER_SHMEM
static rem_meminfo_t *_rem_meminfo;
static aptl_reginfo_t *_tmp_rem_reginfo;

#define IN_REGION(_ptr__,_reg__) ((_reg__.valid) && (_ptr__)>=(_reg__.serv_ptr) \
                && (_ptr__) <= ( (char *)(_reg__.serv_ptr)+_reg__.size))
#endif

static int ptl_initialized = 0;
extern pid_t server_pid;

ptl_ni_limits_t armci_ptl_nilimits;
ptl_ni_limits_t armci_ptl_Snilimits;

void armci_portals_init_ptl()
{
int rc;
int npes,i;
    ARMCI_PR_DBG("enter",0);

    /*initialize data structures*/
    portals->ptl = ARMCI_PORTALS_PTL_NUMBER; /* our own ptl number */

    rc=PtlNIInit(IFACE_FROM_BRIDGE_AND_NALID(PTL_BRIDGE_UK,PTL_IFACE_SS), 
                    PTL_PID_ANY, NULL, &armci_ptl_nilimits, &(portals->ni_h));
    switch(rc) {
      case PTL_OK:
        /*printf("\n%d:ok for nii\n",armci_me);*/
        break;
      case PTL_IFACE_DUP:
        /*printf("\n%d:dup for nii\n",armci_me);*/
        break;
      default:
        printf( "PtlNIInit() failed %d error=%s\n",rc,ARMCI_NET_ERRTOSTR(rc) );
        exit(1);
    }

    if((rc=PtlGetId(portals->ni_h,&portals->rank)) !=PTL_OK) {
      printf("%s: PtlGetId failed: %d(%d)\n",__FUNCTION__, rc, server_pid);
      exit(1);
    }
    ARMCI_PR_DBG("exit",0);
}

static inline void init_serv_buf(serv_buf_t *tmp)
{
int rc;
ptl_match_bits_t ignbits = 0xFFFFFFFFF00000FF;
ptl_match_bits_t mbits;
ptl_process_id_t match_id;
ptl_md_t *md_ptr,md;

    ARMCI_PR_DBG("enter",0);
    tmp->md.user_ptr=tmp;
    tmp->md.start=tmp->buf;
    tmp->md.length=armci_nproc*NUM_SERV_BUFS*VBUF_DLEN;
    tmp->md.eq_handle=portals->Seq_h;
    tmp->md.max_size=0;
    tmp->md.threshold=PTL_MD_THRESH_INF;
    tmp->md.options=PTL_MD_OP_GET | PTL_MD_OP_PUT | PTL_MD_EVENT_START_DISABLE;
    {
      match_id.nid = PTL_NID_ANY;
      match_id.pid = PTL_PID_ANY;
      mbits = 16<<8;

      rc = PtlMEAttach(portals->Sni_h,portals->ptl,match_id,mbits,ignbits,
                      PTL_RETAIN,PTL_INS_AFTER,&(tmp->me_h)); 
      if (rc != PTL_OK) {
        printf("(%d):PtlMEAttach: %s\n", portals->Srank,ARMCI_NET_ERRTOSTR(rc));
        armci_die("portals attach error isb",rc);
      }
      tmp->md.options=PTL_MD_OP_GET | PTL_MD_OP_PUT | PTL_MD_EVENT_START_DISABLE | PTL_MD_MANAGE_REMOTE;
 
      rc = PtlMDAttach((tmp->me_h),tmp->md,PTL_RETAIN,&(tmp->md_h));
      if (rc != PTL_OK) {
        printf("%d:PtlMDAttach: %s %d\n", portals->Srank, ARMCI_NET_ERRTOSTR(rc),(serv_md_count+client_md_count) );
        exit(1);
      }
      serv_md_count++;
    }
    /*set up for sending acks */
    md_ptr            = &(md);
    md_ptr->start     = serv_ack_ptr;
    md_ptr->length    = sizeof(long);
    md_ptr->threshold = PTL_MD_THRESH_INF;
    md_ptr->options   =  PTL_MD_OP_PUT | PTL_MD_EVENT_START_DISABLE;
    md_ptr->user_ptr  = NULL;
    md_ptr->max_size  = sizeof(long);
    md_ptr->eq_handle = portals->Seq_h;

    rc = PtlMDBind(portals->Sni_h,md,PTL_RETAIN,&serv_ack_md_h);
    if (rc != PTL_OK){
      fprintf(stderr, "%d:PtlMDBindxn: %s %d\n", portals->Srank.nid, 
                      ARMCI_NET_ERRTOSTR(rc),(serv_md_count+client_md_count)); 
      armci_die("ptlmdbind failed",0);
    }
    serv_md_count++;
    /*set up for sending response */
    md_ptr            = &(md);
    md_ptr->start     = tmp->buf;
    md_ptr->length    = tmp->md.length;
    md_ptr->threshold = PTL_MD_THRESH_INF;
    md_ptr->options   =  PTL_MD_OP_PUT | PTL_MD_EVENT_START_DISABLE;
    md_ptr->user_ptr  = NULL;
    md_ptr->max_size  = tmp->md.length;
    md_ptr->eq_handle = portals->Seq_h;

    rc = PtlMDBind(portals->Sni_h, md, PTL_RETAIN, &serv_response_md_h);
    if (rc != PTL_OK){
      fprintf(stderr, "%d:PtlMDBindxn: %s %d\n", portals->Srank.nid, 
                      ARMCI_NET_ERRTOSTR(rc),(serv_md_count+client_md_count)); 
      armci_die("ptlmdbind failed",0);
    }
    serv_md_count++;
    ARMCI_PR_DBG("exit",0);
}

void armci_portals_wait_for_client()
{
int rc;  
int *procidinfo;
extern armci_clus_t *armci_clus_info;
ptl_process_id_t *tmp;
    ARMCI_PR_SDBG("enter",0);
    //printf(" ");
    armci_util_wait_int(&_armci_portals_client_ready,1,1000);
    if((armci_me)!=armci_master){
      exit(0);
    }
    else{
      if(DEBUG_SERV){
        printf("\n%d:chosen one nid,pid=%d,%d\n",armci_me,portals->Srank.nid,portals->Srank.pid);
      }
    }
    ARMCI_PR_SDBG("exit",0);
}


void armci_portals_prepare_server()
{
int rc,i,j;
    ARMCI_PR_SDBG("enter",0);
    serv_bufs=(serv_buf_t *)malloc(sizeof(serv_buf_t));
    bzero(serv_bufs,sizeof(serv_buf_t));
    assert(serv_bufs);
    serv_bufs->buf=(char *)malloc((NUM_SERV_BUFS*armci_nproc*VBUF_DLEN));
    bzero(serv_bufs->buf,(NUM_SERV_BUFS*armci_nproc*VBUF_DLEN));
    assert(serv_bufs->buf);
    serv_bufs->bufend=(char *)serv_bufs->buf+(NUM_SERV_BUFS*armci_nproc*VBUF_DLEN);
    rc = PtlEQAlloc(portals->Sni_h,4*(NUM_SERV_BUFS*armci_nproc),NULL, &(portals->Seq_h));
    if (rc != PTL_OK) {
      printf("(%d):Ptleaalloc() failed: %s %d (%d)\n",portals->Srank, 
                      ARMCI_NET_ERRTOSTR(rc),(NUM_SERV_BUFS*armci_nproc),rc);
      armci_die("EQ Alloc failed",rc);
    }
    init_serv_buf(serv_bufs);
    _armci_portals_server_ready=1;
    ARMCI_PR_SDBG("exit",0);
}


void *armci_server_code(void *data)
{
int rc,num_interface;
    ARMCI_PR_SDBG("enter",0);
    if(DEBUG_INIT)
        printf("%d: in server after creating thread.\n",armci_me);

    rc = PtlInit(&num_interface);
    if (rc != PTL_OK) {
       printf("PtlInit() failed %d %s\n",rc, ARMCI_NET_ERRTOSTR(rc) );
       exit(1);
    }

    rc=PtlNIInit(IFACE_FROM_BRIDGE_AND_NALID(PTL_BRIDGE_UK,PTL_IFACE_SS), 
                    PTL_PID_ANY, NULL, &armci_ptl_Snilimits, &(portals->Sni_h));
    switch(rc) {
      case PTL_OK:
        //printf("\n(%d):ok for serv nii\n",armci_me);
        break;
      case PTL_IFACE_DUP:
        //printf("\n(%d):dup for serv nii\n",armci_me);
        break;
      default:
        printf( "PtlNIInit() serv failed %d error=%s\n",rc,ARMCI_NET_ERRTOSTR(rc) );
        exit(1);
    }

    if((rc=PtlGetId(portals->Sni_h,&portals->Srank)) !=PTL_OK) {
      printf("%s: PtlGetId failed: %d(%d)\n",__FUNCTION__, rc, server_pid);
      exit(1);
    }
    /*printf("\n(%d):server nid=%d pid=%d\n",armci_me,portals->Srank.nid,portals->Srank.pid);*/

    armci_portals_wait_for_client();
    armci_portals_prepare_server();

    if(DEBUG_INIT) {
        printf("(%d): connected to all computing processes\n",armci_me);
        fflush(stdout);
    }
    armci_call_data_server();

    armci_transport_cleanup();
    ARMCI_PR_SDBG("exit",0);
    return(NULL);
}


void armci_client_connect_to_servers()
{
int rc;
ptl_size_t offset_local = 0, offset_remote=0;
ptl_md_t md_local;
ptl_handle_md_t md_hdl_local;
ptl_process_id_t *tmp;
int *procidinfo;
int c_info;
int *flag,shmid;
void *addr;
char *buf;
extern int _armci_server_started;
    ARMCI_PR_DBG("enter",0);

    _armci_portals_client_ready=1;
    if(armci_me==armci_master){
      armci_util_wait_int(&_armci_portals_server_ready,1,1000);
    }

    armci_msg_barrier();
    _armci_server_started=1;

    if(armci_me==armci_master){
      portals->servid_map[armci_clus_me].pid=portals->Srank.pid;
      portals->servid_map[armci_clus_me].nid=portals->Srank.nid;
    }

    armci_msg_gop_scope(SCOPE_ALL,portals->servid_map,(sizeof(ptl_process_id_t)*armci_nclus)/sizeof(int),"+",ARMCI_INT);

    ARMCI_PR_DBG("exit",0);
}

static int check_meminfo(void *ptr, long size, int proc)
{
  for(int i=0;i<all_meminfo[proc].cur_ds;i++){
    long left = (caddr_t)ptr - all_meminfo[proc].ptr[i];
    long right= all_meminfo[proc].size[i]-left;
#ifdef DEBUG_MEM
    printf("\n%d:%s:proc=%d curds=%d/%d ptr=%p base=%p top=%p left=%ld right=%ld size=%ld",armci_me,__FUNCTION__,proc,i,all_meminfo[proc].cur_ds,ptr,all_meminfo[proc].ptr[i],(all_meminfo[proc].ptr[i]+all_meminfo[proc].size[i]),left,right,size);fflush(stdout);
#endif
    if((left>=0) && (right>=size))
      return(i+1);
  }
  return 0;
}


static void add_meminfo(void *ptr, size_t size, int proc)
{
  if(check_meminfo(ptr,(long)size,proc)!=0)armci_die("repeat add request for dss",proc);
  all_meminfo[proc].cur_ds++;
  all_meminfo[proc].ptr[all_meminfo[proc].cur_ds]=ptr;
  all_meminfo[proc].size[all_meminfo[proc].cur_ds]=size;
#ifdef DEBUG_MEM
  printf("\n%d:%s:adding %p %ld %d at %d",armci_me,__FUNCTION__,ptr,size,proc,all_meminfo[proc].cur_ds);
#endif
}


typedef struct{
        void *ptr;
        size_t size;
        size_t serv_offs;
} meminfo_t;


void armci_exchange_meminfo(void *ptr, size_t size,size_t off)
{
  static meminfo_t exng[armci_nproc];
  bzero(exng,sizeof(meminfo_t)*armci_nproc);
  exng[armci_me].ptr=ptr; exng[armci_me].size=size; exng[armci_me].serv_offs = off;
  armci_msg_gop_scope(SCOPE_ALL,exng,(sizeof(meminfo_t)*armci_nproc)/sizeof(int),"+",ARMCI_INT);
  for(int i=0;i<armci_nproc;i++){
    if(exng[i].size!=0){
      add_meminfo(exng[i].ptr,exng[i].size,i);
    }
  }
}


static caddr_t get_heap_bottom_addr()                                                                                       
{                                                                                                                           
    extern caddr_t _end;
    printf("\n%d:%s:_end=%p addr_end=%p",armci_me,__FUNCTION__,_end,&_end);fflush(stdout);
    return((caddr_t)&_end);                                                                                                 
}       


void armci_portals_memsetup(long serv_offset)
{
  void *br_val,*ptr=NULL;
  size_t size=0;
  br_val = sbrk(0);
  int myerr;

  if(br_val!=portals->brval[portals->cur_ds]){
  ptl_md_t *md_ptr;
  ptl_match_bits_t ignbits = 0xFFFFFFFFFFFFFF00;
  ptl_process_id_t match_id;
  int rc,cds=++portals->cur_ds;
    if(cds>=MAX_DS)armci_die("increase MAX_CDS",cds);
    portals->dsbase[cds]=portals->brval[cds-1];
    //portals->dsbase[cds]=sbrk(0);
    ptr = portals->brval[cds] = br_val;
    size = portals->dssizes[cds]=((caddr_t)portals->brval[cds] - portals->dsbase[cds]);
    portals->serv_offs[cds] = serv_offset;
    printf("\n%d:%s:base=%p brval=%p dslen=%ld %p end=%p",armci_me,__FUNCTION__,portals->dsbase[cds],br_val,
                      portals->dssizes[cds],portals->brval[cds],get_heap_bottom_addr());

    md_ptr            = &(portals->heap_md[cds]);
    md_ptr->start     = portals->dsbase[cds];
    md_ptr->length    = portals->dssizes[cds];
    md_ptr->threshold = PTL_MD_THRESH_INF;
    md_ptr->options   =  PTL_MD_OP_PUT | PTL_MD_OP_GET | PTL_MD_MANAGE_REMOTE;
    md_ptr->user_ptr  = NULL;
    md_ptr->max_size  = 0;
    md_ptr->eq_handle = PTL_EQ_NONE;

    portals->heap_mb[cds]=cds+1;
 
    match_id.nid = PTL_NID_ANY;
    match_id.pid = PTL_PID_ANY; 
    rc = PtlMEAttach(portals->ni_h,portals->ptl,match_id,
                    portals->heap_mb[cds],
                    ignbits,
                    PTL_RETAIN,PTL_INS_AFTER,
                    &(portals->heap_me_h[cds])); 
    if (rc != PTL_OK) {
      printf("%d:PtlMEAttach: %s\n", portals->rank, ARMCI_NET_ERRTOSTR(rc) );
      armci_die("portals attach error reg",rc);
    } 

    rc = PtlMDAttach((portals->heap_me_h[cds]),
                    *md_ptr,PTL_RETAIN,
                    &(portals->heap_md_h[cds]));
    if (rc != PTL_OK) {
      printf("%d:PtlMDAttach: %s %d\n", portals->rank, ARMCI_NET_ERRTOSTR(rc),(client_md_count+serv_md_count));
      armci_die("portals attach error reg",rc);
    }
  }
  else{
#ifdef DEBUG_MEM_
    extern caddr_t _end;
    printf("\n%d:%s:curds=%d brvalin=%p curbrval=%p _end=%p &_end=%p",armci_me,__FUNCTION__,portals->cur_ds,portals->brval[portals->cur_ds],br_val,_end,&_end);
#endif
  }
  armci_exchange_meminfo(ptr,size,serv_offset);
}


#ifndef PMI_SUCCESS
#define PMI_SUCCESS 0
#endif

static int *_client_servbuf_count;
int armci_init_portals(caddr_t atbeginbrval)
{
#ifndef OLD_PORTALS_CODE
        int i,rc,np,me;
        ptl_process_id_t id;
        ptl_process_id_t clone_id;

        MPI_Comm_size(ARMCI_COMM_WORLD,&np);
        MPI_Comm_rank(ARMCI_COMM_WORLD,&me);

        if(armci_me != me) {
           printf("[mpi %d]: armci_me=%d ... this is a problem\n",me,armci_me);
           armci_die("mpi rank mismatch",911);
        }
        
        portals_cp_init();

        MPI_Barrier(ARMCI_COMM_WORLD);

        portals_ds_ready = 0;
        if(armci_me == armci_master) {
           portalsCloneDataServer( portals_ds_thread );
           portalsSpinLockOnInt( &portals_ds_ready,1,10000 );
        }
        MPI_Barrier(ARMCI_COMM_WORLD);

        i=0;
        if((rc=PMI_Initialized(&i))!=PMI_SUCCESS){
          printf("PMI_Initialized failed\n");
        }
          
        if(i==0){
          if((rc==PMI_Init(&i))!=PMI_SUCCESS){
            printf("MPI_Init failed (npes=%d)\n", armci_nproc);
          }
        }
         
        if((rc=PMI_CNOS_Get_nidpid_map(&portals_id_map))!=PMI_SUCCESS){
          printf("Getting proc map failed (npes=%d)\n", armci_nproc);
        }

     /* create intra-node communicator */
        MPI_Barrier(ARMCI_COMM_WORLD);
        portals_cp_init_throttle(armci_nclus);
        MPI_Barrier(ARMCI_COMM_WORLD);

   /* stuff from old code ... */
    bzero(portals,sizeof(armci_portals_proc_t));
    // note: i got rid of this rem_meminfo stuff with the gemini version
    // see that code to see how to remove it here
  # ifdef ARMCI_REGISTER_SHMEM
    _rem_meminfo = (rem_meminfo_t *)calloc(armci_nproc,sizeof(rem_meminfo_t));
    _tmp_rem_reginfo = (aptl_reginfo_t *)malloc(sizeof(aptl_reginfo_t)*armci_nproc);
    if( _rem_meminfo==NULL || _tmp_rem_reginfo ==NULL)
      armci_die("malloc failed in init_portals",0);
    //if(armci_me == 0) {
    //  printf("sizeof(rem_meminfo_t)=%ld\n",sizeof(rem_meminfo_t));
    //}
  # endif
  # ifdef CRAY_USE_ARMCI_CLIENT_BUFFERS
    client_buf_ptrs = (char **) calloc(armci_nproc,sizeof(char *));
    assert(client_buf_ptrs);
    armci_msg_barrier();
    _armci_buf_init();
  # endif
   /* end old stuff */

        return 0;

#else

int num_interface;
int rc;
int npes,i;
    ARMCI_PR_DBG("enter",0);
   
    bzero(portals,sizeof(armci_portals_proc_t));

    _rem_meminfo = (rem_meminfo_t *)calloc(armci_nproc,sizeof(rem_meminfo_t));
    _tmp_rem_reginfo = (aptl_reginfo_t *)malloc(sizeof(aptl_reginfo_t)*armci_nproc);
    if( _rem_meminfo==NULL || _tmp_rem_reginfo ==NULL)
      armci_die("malloc failed in init_portals",0);

    portals->servid_map=(ptl_process_id_t *)calloc(armci_nclus,sizeof(ptl_process_id_t));
    if(portals->servid_map==NULL)armci_die("calloc of servidmap failed",0);

    rc = PtlInit(&num_interface);
    if (rc != PTL_OK) {
       printf("PtlInit() failed %d %s\n",rc, ARMCI_NET_ERRTOSTR(rc) );
       exit(1);
    }
    armci_portals_init_ptl();

#if 1
    i=0;
    if((rc=PMI_Initialized(&i))!=PMI_SUCCESS){
      printf("PMI_Initialized failed\n");
    }

    if(i==0){
      if((rc==PMI_Init(&i))!=PMI_SUCCESS){
        printf("MPI_Init failed (npes=%d)\n", armci_nproc);
      }
    }

    if((rc=PMI_CNOS_Get_nidpid_map(&portals->procid_map))!=PMI_SUCCESS){
      printf("Getting proc map failed (npes=%d)\n", armci_nproc);
    }
    //printf(" ");
# else
    portals->procid_map = (ptl_process_id_t *) calloc(armci_nproc,sizeof(ptl_process_id_t));
    portals->procid_map[armci_me]=portals->rank;
    armci_msg_gop_scope(SCOPE_ALL,portals->procid_map,(sizeof(ptl_process_id_t)*armci_nproc)/sizeof(int),"+",ARMCI_INT);
    //printf(" ");
#endif

    client_buf_ptrs = (char **) calloc(armci_nproc,sizeof(char *));
    assert(client_buf_ptrs);
    armci_msg_barrier();

    if(armci_me==armci_master)armci_create_server_process( armci_server_code );

    rc = PtlEQAlloc(portals->ni_h,16*NUM_COMP_DSCR,NULL, &(portals->eq_h));
    if (rc != PTL_OK) {
      printf("%d:Ptleqalloc() failed: %s (%d)\n",
                      armci_me, ARMCI_NET_ERRTOSTR(rc) , rc);
      armci_die("EQ Alloc failed",rc);
    }
    /*printf("\n%d:client nid=%d pid=%d\n",armci_me,portals->rank.nid,portals->rank.pid);*/

    _armci_buf_init();

    for(i=0;i<NUM_COMP_DSCR;i++){
      _compdesc_array[i].active=0;
      _compdesc_array[i].tag=-1;
      _compdesc_array[i].dest_id=-1;
      _compdesc_array[i].mem_dsc.eq_handle=portals->eq_h;
      _compdesc_array[i].mem_dsc.max_size=0;
      /*_compdesc_array[i].mem_dsc.threshold=PTL_MD_THRESH_INF;*/
      _compdesc_array[i].mem_dsc.threshold=2;
      _compdesc_array[i].mem_dsc.options=PTL_MD_OP_GET | PTL_MD_OP_PUT | PTL_MD_EVENT_START_DISABLE;
    }

    ptl_initialized = 1;
    portals->free_comp_desc_index=0;
    /*for(i=0;i<armci_nproc;i++)printf("\t%d:%d,%d",i,portals->procid_map[i].nid,portals->procid_map[i].pid);*/
    _client_servbuf_count = calloc(armci_nclus,sizeof(int));
    armci_msg_barrier();
    armci_client_connect_to_servers();
    armci_msg_barrier();
    if(DEBUG_COMM){
    cpu_set_t mycpuid,new_mask;
    char cid[8],*cidptr;
    int rrr;
    extern char * cpuset_to_cstr(cpu_set_t *mask, char *str);
      rrr=sched_getaffinity(0, sizeof(mycpuid), &mycpuid);
      if(rrr)perror("sched_getaffinity");
      cidptr = cpuset_to_cstr(&mycpuid,cid);
      printf("%d:my affinity is to %s\n",armci_me,cid);
    }
#ifdef NEW_MALLOC
    /*post entire heap wildcard for direct communication*/
    {
    ptl_md_t *md_ptr;
    ptl_match_bits_t ignbits = 0xFFFFFFFFFFFFFF00;
    ptl_process_id_t match_id;
      portals->cur_ds = 0; 
      portals->dsbase[0]=get_heap_bottom_addr();
      //portals->brval[0] = sbrk(0);
      portals->brval[0] = atbeginbrval;
      portals->dssizes[0]=((caddr_t)portals->brval[0] - portals->dsbase[0]);
      printf("\n%d:base=%p dslen=%ld %p",armci_me,portals->dsbase[0],
                      portals->dssizes[0],portals->brval[0]);

      md_ptr            = &(portals->heap_md[0]);
      md_ptr->start     = portals->dsbase[0];
      md_ptr->length    = portals->dssizes[0];
      md_ptr->threshold = PTL_MD_THRESH_INF;
      md_ptr->options   =  PTL_MD_OP_PUT | PTL_MD_OP_GET | PTL_MD_MANAGE_REMOTE;
      md_ptr->user_ptr  = NULL;
      md_ptr->max_size  = 0;
      md_ptr->eq_handle = PTL_EQ_NONE;

      portals->heap_mb[0]=1;
 
      match_id.nid = PTL_NID_ANY;
      match_id.pid = PTL_PID_ANY; 
      rc = PtlMEAttach(portals->ni_h,portals->ptl,match_id,
                    portals->heap_mb[0],
                    ignbits,
                    PTL_RETAIN,PTL_INS_AFTER,
                    &(portals->heap_me_h[0])); 
      if (rc != PTL_OK) {
        printf("%d:PtlMEAttach: %s\n", portals->rank, ARMCI_NET_ERRTOSTR(rc) );
        armci_die("portals attach error reg",rc);
      } 

      rc = PtlMDAttach((portals->heap_me_h[0]),
                    *md_ptr,PTL_RETAIN,
                    &(portals->heap_md_h[0]));
      if (rc != PTL_OK) {
        printf("%d:PtlMDAttach: %s %d\n", portals->rank, ARMCI_NET_ERRTOSTR(rc),(client_md_count+serv_md_count));
        armci_die("portals attach error reg",rc);
      }
      all_meminfo = (rm_info_t *)malloc(sizeof(rm_info_t)*armci_nproc);
      all_meminfo[armci_me].cur_ds = -1;
      armci_exchange_meminfo(portals->dsbase[0],portals->dssizes[0],0);
    }
#endif
    ARMCI_PR_DBG("exit",0);
    return 0;   
#endif
}


void armci_fini_portals()
{
    ARMCI_PR_DBG("enter",0);
    if(DEBUG_INIT){
      printf("ENTERING ARMCI_FINI_PORTALS\n");fflush(stdout);
    }    
#ifdef ARMCI_CHECK_STATE
    arlist_print(arn);
#endif
    PtlNIFini(portals->ni_h);
    /*PtlFini();*/
    if(DEBUG_INIT){
      printf("LEAVING ARMCI_FINI_PORTALS\n");fflush(stdout);    
    }
    ARMCI_PR_DBG("exit",0);
}


void armci_pin_contig1(void *start,size_t bytes)
{

}

#ifdef ARMCI_REGISTER_SHMEM
#ifndef NEW_MALLOC
void armci_register_req(void *start,int bytes, int ID)
{
int rc;
ptl_md_t *md_ptr;
ptl_match_bits_t ignbits = 0xFFFFFFFFFFFFFF00;
ptl_process_id_t match_id;

    ARMCI_PR_DBG("enter",serv_portals->reg_count);

#ifdef DEBUG_MEM
      printf("\n(%d):armci_register_req start=%p bytes=%d\n",
                      armci_me,start,bytes);fflush(stdout);
#endif

    md_ptr            = &(serv_portals->meminfo[serv_portals->reg_count].md);
    md_ptr->start     = start;
    md_ptr->length    = bytes;
    md_ptr->threshold = PTL_MD_THRESH_INF;
    md_ptr->options   =  PTL_MD_OP_PUT | PTL_MD_OP_GET | PTL_MD_MANAGE_REMOTE;
    md_ptr->user_ptr  = NULL;
    md_ptr->max_size  = 0;
    md_ptr->eq_handle = PTL_EQ_NONE;

    serv_portals->meminfo[serv_portals->reg_count].mb=serv_portals->reg_count+1;
 
    match_id.nid = PTL_NID_ANY;
    match_id.pid = PTL_PID_ANY; 
    rc = PtlMEAttach(portals->ni_h,portals->ptl,match_id,
                    serv_portals->meminfo[serv_portals->reg_count].mb,
                    ignbits,
                    PTL_RETAIN,PTL_INS_AFTER,
                    &(serv_portals->meminfo[serv_portals->reg_count].me_h)); 
    if (rc != PTL_OK) {
      printf("%d:PtlMEAttach: %s\n", portals->rank, ARMCI_NET_ERRTOSTR(rc) );
      armci_die("portals attach error reg",rc);
    }

    rc = PtlMDAttach((serv_portals->meminfo[serv_portals->reg_count].me_h),
                    *md_ptr,PTL_RETAIN,
                    &serv_portals->meminfo[serv_portals->reg_count].md_h);
    if (rc != PTL_OK) {
      printf("%d:PtlMDAttach: %s %d\n", portals->rank, ARMCI_NET_ERRTOSTR(rc),(client_md_count+serv_md_count));
      armci_die("portals attach error reg",rc);
    }
    client_md_count++;
    serv_portals->reg_count++;     

    ARMCI_PR_DBG("exit",serv_portals->reg_count);
}
#endif
#endif

int armci_must_remotecomplete=1;
extern _buf_ackresp_t *_buf_ackresp_first,*_buf_ackresp_cur;
int x_net_wait_ackresp(_buf_ackresp_t *ar)
{
int rc;  
ptl_event_t ev_t;
ptl_event_t *ev=&ev_t;
comp_desc *temp_comp = NULL;
int loop=1;
int temp_proc;
    ARMCI_PR_DBG("enter",0);
    while(ar->val){
      ev->type=0;
      if((rc = PtlEQWait(portals->eq_h, ev)) != PTL_OK){
        printf("%d:PtlEQWait(): %d %s\n", portals->rank,rc,
                        ARMCI_NET_ERRTOSTR(rc) ); 
        armci_die("EQWait problem",rc);
      }
      if (ev->ni_fail_type != PTL_NI_OK) {
        temp_comp = (comp_desc *)ev->md.user_ptr;
        printf("%d:NI sent %d in event %d,%d.\n",
         armci_me,portals->rank.nid, portals->rank.pid,  ev->ni_fail_type); 
        armci_die("event failure problem",temp_comp->dest_id);
      }
      if(DEBUG_COMM){
        printf("\n%d:net_wait_ackresp:done waiting type=%d\n",armci_me,
                        ev->type);
        fflush(stdout);
      }
      if (ev->type == PTL_EVENT_SEND_END){
        if(DEBUG_COMM){
          printf("\n%d:net_wait_ackresp:event send end\n",armci_me);
          fflush(stdout);
        }
        temp_comp = (comp_desc *)ev->md.user_ptr;
        if(temp_comp->type==ARMCI_PORTALS_GETPUT || temp_comp->type==ARMCI_PORTALS_NBGETPUT){
          temp_comp->active=0;
          temp_comp->tag=-1;
          continue;
        }
        if(!armci_must_remotecomplete){
          if(temp_comp->type==ARMCI_PORTALS_PUT || temp_comp->type==ARMCI_PORTALS_NBPUT){
            temp_comp->active=0;
            temp_comp->tag=-1;
          }
          else
            continue;
        }
        else{
          temp_comp->active++;
          continue;
        }


      }

      else if (ev->type == PTL_EVENT_REPLY_END){
        temp_comp = (comp_desc *)ev->md.user_ptr;
        if(DEBUG_COMM){
          printf("\n%d:net_wait_ackresp:reply end tag=%d\n",armci_me,temp_comp->tag);
          fflush(stdout);
        }
        temp_comp->active = 0; /*this was a get request, so we are done*/
        temp_comp->tag=-1;
        continue;
      }
      else if (ev->type == PTL_EVENT_ACK){
        temp_comp = (comp_desc *)ev->md.user_ptr;
        if(DEBUG_COMM){
          printf("\n%d:net_wait_ackresp:event ack tag=%d\n",armci_me,temp_comp->tag);
          fflush(stdout);
        }
        temp_comp->active=0;
        temp_comp->tag=-1;
        portals->outstanding_puts--; 
      }
      else if (ev->type==PTL_EVENT_PUT_END){
        _buf_ackresp_t *sweep=_buf_ackresp_first;
        if(DEBUG_COMM){printf("\n%d:put end offset=%d",armci_me,ev->offset);fflush(stdout);}
        if(ar->val==ev->offset){
          /*bingo!*/
          ar->val=0;
        }
        else{
          while(sweep!=NULL){
            if(sweep->val==ev->offset){
              sweep->val=0;
              break;
            }
            sweep=sweep->next;
          }
          /*if(sweep==NULL)armci_die("server wrote data at unexpected offset",ev->offset);*/
          if(sweep==NULL){
            int y;
            printf("%d:server wrote data at unexpected offset %d",armci_me,ev->offset);fflush(stdout);
            abort();
#          ifdef ARMCI_CHECK_STATE
            for(y=0;y<armci_nclus;y++)
              if(portals->servid_map[y].pid==ev->initiator.pid && portals->servid_map[y].nid==ev->initiator.nid)break;
            assert(y!=armci_nclus); 
            arlist_print(arn);
            armci_rem_state(y);
#          endif
          }
        }
      }
      else 
        armci_die("in net_wait_ackresp unknown event",ev->type);
    }

#  ifdef ARMCI_CHECK_STATE               
    arlist_remove(arlist_search(&arn, ar->valc));
#  endif
    ar->valc=0;
    if(ar==_buf_ackresp_first)_buf_ackresp_first=ar->next;
    if(ar->next!=NULL){
      ar->next->previous=ar->previous;
    }
    if(ar->previous!=NULL){
      /*printf("\n%d:prev=%p %d %p %p\n",armci_me,ar->previous, ar->val,ar->next,ar);fflush(stdout);*/
      ar->previous->next=ar->next;
      if(_buf_ackresp_cur==ar)_buf_ackresp_cur=ar->previous;
    }
    if(_buf_ackresp_cur==ar)_buf_ackresp_cur=NULL;

    ar->previous=ar->next=NULL;

    ARMCI_PR_DBG("exit",0);

    return rc; 
}
int armci_client_complete(ptl_event_kind_t evt,int proc_id, int nb_tag,
                comp_desc *cdesc)
{
int rc;  
ptl_event_t ev_t;
ptl_event_t *ev=&ev_t;
comp_desc *temp_comp = NULL;
int loop=1;
int temp_proc;
    ARMCI_PR_DBG("enter",0);
    if(DEBUG_COMM){
      printf("\n%d:enter:client_complete active=%d tag=%d %d\n",armci_me,
                      cdesc->active,cdesc->tag,nb_tag);fflush(stdout);
    }
    if(nb_tag>0){
      if(cdesc->tag!=nb_tag)return 0;
    }
    while(cdesc->active!=0){
      ev->type=0;
      if((rc = PtlEQWait(portals->eq_h, ev)) != PTL_OK){
        printf("%d:PtlEQWait(): %d %s\n", portals->rank,rc,
                        ARMCI_NET_ERRTOSTR(rc) ); 
        armci_die("EQWait problem",rc);
      }
      if (ev->ni_fail_type != PTL_NI_OK) {
        temp_comp = (comp_desc *)ev->md.user_ptr;
        printf("%d:NI sent %d in event %d,%d.\n",
         armci_me,portals->rank.nid, portals->rank.pid,  ev->ni_fail_type); 
        armci_die("event failure problem",temp_comp->dest_id);
      }
      if(DEBUG_COMM){
        printf("\n%d:armci_client_complete:done waiting type=%d\n",armci_me,
                        ev->type);
        fflush(stdout);
      }
      if(cdesc!=ev->md.user_ptr){
        /*printf("\n%d:expecting desc %p completing %p\n",armci_me,cdesc,ev->md.user_ptr);*/
      }
      if (ev->type == PTL_EVENT_SEND_END){
        if(DEBUG_COMM){
          printf("\n%d:armci_client_complete:event send end\n",armci_me);
          fflush(stdout);
        }
        temp_comp = (comp_desc *)ev->md.user_ptr;
        if(temp_comp->type==ARMCI_PORTALS_GETPUT || temp_comp->type==ARMCI_PORTALS_NBGETPUT){
          temp_comp->active=0;
          temp_comp->tag=-1;
          continue;
        }
        if(!armci_must_remotecomplete){
          if(temp_comp->type==ARMCI_PORTALS_PUT || temp_comp->type==ARMCI_PORTALS_NBPUT){
            temp_comp->active=0;
            temp_comp->tag=-1;
          }
          else
            continue;
        }
        else{
          temp_comp->active++;
          continue;
        }
      }

      else if (ev->type == PTL_EVENT_REPLY_END){
        temp_comp = (comp_desc *)ev->md.user_ptr;
        if(DEBUG_COMM){
          printf("\n%d:client_send_complete:reply end tag=%d\n",armci_me,temp_comp->tag);
          fflush(stdout);
        }
        temp_comp->active = 0; /*this was a get request, so we are done*/
        temp_comp->tag=-1;
        continue;
      }
      else if (ev->type == PTL_EVENT_ACK){
        temp_comp = (comp_desc *)ev->md.user_ptr;
        if(DEBUG_COMM){
          printf("\n%d:client_send_complete:event ack tag=%d\n",armci_me,temp_comp->tag);
          fflush(stdout);
        }
        temp_comp->active=0;
        temp_comp->tag=-1;
        portals->outstanding_puts--; 
      }
      else if (ev->type==PTL_EVENT_PUT_END){
        _buf_ackresp_t *ar=_buf_ackresp_first;
        while(ar!=NULL){
          if(ar->val==ev->offset){
            ar->val=0;
            break;
          }
          ar=ar->next;
        }
        if(ar==NULL)armci_die("server wrote data at unexpected offset",ev->offset);      
        if(DEBUG_COMM){printf("\n%d:put end offset=%d",armci_me,ev->offset);fflush(stdout);}
      }
      else 
        armci_die("in client_complete unknown event",ev->type);
    }
    if(DEBUG_COMM){
      printf("\n%d:exit:client_complete active=%d tag=%d %d\n",armci_me,
                      cdesc->active,cdesc->tag,nb_tag);fflush(stdout);
    }

    ARMCI_PR_DBG("exit",0);

    return rc; 
}


comp_desc * get_free_comp_desc(int * comp_id)
{
comp_desc * c;     
int rc = PTL_OK;

    ARMCI_PR_DBG("enter",0);

    c = &(_compdesc_array[portals->free_comp_desc_index]);
    if(c->active!=0 && c->tag>0)
      armci_client_complete(0,c->dest_id,c->tag,c);
    else{
      /*
      if(c->active!=0)
        printf("\n%d:potential problem:active completion descriptor but tag=%d",armci_me,c->tag);
      else
        printf("\n%d:potential problem:active completion descriptor with tag=%d",armci_me,c->tag);
        */
    }
    if(!armci_must_remotecomplete){
      do{
           rc = PtlMDUnlink(c->mem_dsc_hndl);
      }while(rc==PTL_MD_IN_USE);
    }
    
    *comp_id = portals->free_comp_desc_index;
    if(DEBUG_COMM){
      printf("\nthe value of comp_desc_id is %d\n",*comp_id);
      fflush(stdout);
    }
    portals->free_comp_desc_index = (portals->free_comp_desc_index+1) % NUM_COMP_DSCR;

    ARMCI_PR_DBG("exit",0);

    return c;
}


void print_mem_desc(ptl_md_t * md)
{
    printf("%d:%p:start=%p length=%d threshold=%d max_size=%d options=%d eq_handle=%d\n",armci_me,md,md->start, md->length,md->threshold,md->max_size,md->options,md->eq_handle);
    fflush(stdout);
}


#ifndef NEW_MALLOC
#if 0
void armci_unregister_shmem(void *my_ptr, long size)
{
int i=0,dst,found=0;
long id ;
long reg_size=0;
int reg_num = _rem_meminfo[armci_me].reg_count;
void *tptr;
     
    ARMCI_PR_DBG("enter",reg_num);
#ifdef DEBUG_MEM
    printf("%d:%s:got size=%ld myptr %p\n",armci_me,__FUNCTION__,size,my_ptr);
    fflush(stdout);
#endif
    bzero(_tmp_rem_reginfo,sizeof(aptl_reginfo_t)*armci_nproc);
    if(reg_num>=MAX_MEM_REGIONS)
        armci_die("reg_num corrupted",reg_num);
    for(i=0;i<reg_num;i++)
      if(IN_REGION(sptr,_rem_meminfo[armci_me].reginfo[i])){
        found=1;
        break;
      } 
    if(found){

    }
    else
      armci_die("region dissappeared, cannot free it",reg_num);
}
#endif

#ifdef ARMCI_REGISTER_SHMEM
void armci_register_shmem(void *my_ptr, long size, long *idlist, long off,
       void *sptr)
{
int i=0,dst,found=0;
long id=0;
long reg_size=0;
int reg_num = _rem_meminfo[armci_me].reg_count;
extern void *armci_shm_reg_ptr(int);
void *tptr;
ARMCI_Group def_group;
     
    ARMCI_PR_DBG("enter",reg_num);
#ifdef DEBUG_MEM
      printf("%d:%s:got id=%ld size=%ld myptr %p, sptr %p offset=%d\n",armci_me,__FUNCTION__,id,size,my_ptr,sptr,idlist,off);
      fflush(stdout);
#endif
    bzero(_tmp_rem_reginfo,sizeof(aptl_reginfo_t)*armci_nproc);
    if(size){
      if(reg_num>=MAX_MEM_REGIONS)
        armci_die("reg_num corrupted",reg_num);
      for(i=0;i<reg_num;i++)
        if(IN_REGION(sptr,_rem_meminfo[armci_me].reginfo[i])){
          found=1;
          break;
        } 
      if(!found){ 
        /* record new region id */
        if(idlist==NULL){     
        _tmp_rem_reginfo[armci_me].base_ptr = sptr;
        _tmp_rem_reginfo[armci_me].serv_ptr = tptr = (char *)sptr-off;
        _tmp_rem_reginfo[armci_me].size = size;
        }
        else{
        _tmp_rem_reginfo[armci_me].base_ptr = tptr= armci_shm_reg_ptr(reg_num); /* this gives my shmem ptr*/
        _tmp_rem_reginfo[armci_me].serv_ptr = (char *)sptr-off;
        _tmp_rem_reginfo[armci_me].size = armci_shm_reg_size(reg_num,id);
        _tmp_rem_reginfo[armci_me].base_ptr = _tmp_rem_reginfo[armci_me].serv_ptr;
        }
        armci_register_req(tptr,_tmp_rem_reginfo[armci_me].size,0);
        _tmp_rem_reginfo[armci_me].islocal = 0;
#ifdef DEBUG_MEM
        {printf("\n%d:%s:not found serv_ptr=%p base=%p size=%d",armci_me,__FUNCTION__,((char *)sptr-off),_tmp_rem_reginfo[armci_me].base_ptr,_tmp_rem_reginfo[armci_me].size);}
#endif
      }
    }
    ARMCI_Group_get_default(&def_group);
    armci_msg_group_gop_scope(SCOPE_ALL,_tmp_rem_reginfo,(sizeof(aptl_reginfo_t)*armci_nproc/sizeof(int)),"+",ARMCI_INT,&def_group);
    for(i=0;i<armci_nproc;i++)
      if(_tmp_rem_reginfo[i].size){
        reg_num = _rem_meminfo[i].reg_count;
        _rem_meminfo[i].reginfo[reg_num].base_ptr = _tmp_rem_reginfo[i].base_ptr;
        _rem_meminfo[i].reginfo[reg_num].serv_ptr = _tmp_rem_reginfo[i].serv_ptr;
        _rem_meminfo[i].reginfo[reg_num].size = _tmp_rem_reginfo[i].size;
        _rem_meminfo[i].reginfo[reg_num].islocal = _tmp_rem_reginfo[i].islocal;
        _rem_meminfo[i].reginfo[reg_num].valid = 1;
#ifdef DEBUG_MEM
        {printf("\n%d:%s:new region %d, proc=%d base=%p sptr=%p size=%d\n",armci_me,__FUNCTION__,reg_num,i,_tmp_rem_reginfo[i].base_ptr,_tmp_rem_reginfo[i].serv_ptr,_tmp_rem_reginfo[i].size);fflush(stdout);}
#endif
        _rem_meminfo[i].reg_count++;
        if(_rem_meminfo[i].reg_count>=MAX_MEM_REGIONS-1){
          printf("\n%d:more than expected regions -- %d, increase MAX_MEM_REGIONS",armci_me,_rem_meminfo[i].reg_count++);fflush(stdout);
          armci_die2("more than expected regions",_rem_meminfo[i].reg_count,MAX_MEM_REGIONS);
        }
      }
#ifdef DEBUG_MEM
      printf("%d: regist id=%ld found=%d size=%ld reg_num=%d\n",
                      armci_me,id,found,reg_size,reg_num);
      fflush(stdout);
#endif
    ARMCI_PR_DBG("exit",0);
}

void armci_register_shmem_grp(void *my_ptr, long size, long *idlist, long off,
       void *sptr,ARMCI_Group *group)
{
ARMCI_Group orig_group;
    ARMCI_PR_DBG("enter",0);
    ARMCI_Group_get_default(&orig_group);
    ARMCI_Group_set_default(group);
    armci_register_shmem(my_ptr,size,idlist,off,sptr);
    ARMCI_Group_set_default(&orig_group);
    ARMCI_PR_DBG("enter",0);
}
#endif
#endif // end #ifdef ARMCI_REGISTER_SHMEM

static int _get_rem_servinfo(int serv,size_t bytes, size_t* offset)
{
int i;
    ARMCI_PR_DBG("enter",0);
    i = 16<<8;
    *offset=(armci_me*NUM_SERV_BUFS+_client_servbuf_count[serv])*VBUF_DLEN;
    _client_servbuf_count[serv] = (_client_servbuf_count[serv]+1)%NUM_SERV_BUFS;
    ARMCI_PR_DBG("exit",i);
    return i;
}

static int _get_rem_info(int proc, void *ptr,size_t bytes, size_t* offset)
{
#ifdef ARMCI_REGISTER_SHMEM
int i;
    ARMCI_PR_DBG("enter",0);
#ifdef NEW_MALLOC
    i = check_meminfo(ptr,(long)bytes,proc);
    if(i==0){
      printf("\n%d:ptr=%p bytes=%d proc=%d",armci_me,ptr,bytes,proc);
      armci_die("region not found",proc);
    }
    *offset = (size_t)((caddr_t)ptr-(caddr_t)portals->dsbase[i-1]);
    printf("\n%d:ptr=%p dsbase[0]=%p offs=%ld",armci_me,ptr,portals->dsbase[0],*offset);fflush(stdout);
    if(*offset>=0){
      ARMCI_PR_DBG("exit A",(i+1));
      return(i);
    }
#else
rem_meminfo_t *mem_info=&(_rem_meminfo[proc]);
aptl_reginfo_t *memreg = mem_info->reginfo;
    for(i=0;i<mem_info->reg_count;i++){
      /*for now size is not verified*/
      if(DEBUG_COMM){
        printf("\n%d:proc=%d regcount=%d reg=%d base=%p size=%d end=%p checkptr=%p\n",armci_me,proc,mem_info->reg_count,i,memreg[i].base_ptr,memreg[i].size, ((char *)memreg[i].base_ptr+memreg[i].size), ptr);fflush(stdout);
      }
      if((memreg[i].valid) && ptr>= memreg[i].base_ptr && 
                      ptr< ((char *)memreg[i].base_ptr+memreg[i].size)){
        *offset = ((char *)ptr-(char *)memreg[i].base_ptr);
        ARMCI_PR_DBG("exit A",(i+1));
        return (i+1);
      }
    }
#endif
    ARMCI_PR_DBG("exit B",i);
    armci_die("_get_rem_info, rem memory region not found",bytes);
#else
 printf("_get_rem_info called ... this shouldn't happen"); abort();
#endif
}

void armci_client_direct_get(ptl_process_id_t dest_proc,
                ptl_size_t offset_remote, ptl_match_bits_t mb, size_t bytes,
                ptl_md_t *md_local, 
                ptl_handle_md_t *md_hdl_local)
{
int rc;
ptl_size_t offset_local = 0;

    ARMCI_PR_DBG("enter",0);

    if(DEBUG_COMM){
      printf("\n%d:armci_client_direct_get:BYTES = %d\n",armci_me,bytes);
      printf("\n%d:offr=%ld offl=%ld\n",armci_me,offset_remote,offset_local);
      fflush(stdout);
    }

    rc = PtlMDBind(portals->ni_h,*md_local, PTL_UNLINK, md_hdl_local);
    if (rc != PTL_OK){
      printf("%d:PtlMDBind: %s\n", portals->rank, ARMCI_NET_ERRTOSTR(rc) );
      armci_die("ptlmdbind get failed",0);
    }

#ifdef CRAY_USE_MDMD_COPY
    if (dest_proc.nid == portals->rank.nid) {
        rc = PtlMDMDCopy(*md_hdl_local, dest_proc,
                   portals->ptl,
                   0,
                   mb,
                   offset_remote);
    } else {
#endif    
        rc = PtlGetRegion(*md_hdl_local,offset_local,bytes,dest_proc,
                   portals->ptl,
                   0,
                   mb,
                   offset_remote);
#ifdef CRAY_USE_MDMD_COPY
    }
#endif

    if (rc != PTL_OK){
      printf("%d:PtlGetRegion: %s\n", portals->rank, ARMCI_NET_ERRTOSTR(rc) );
      armci_die("PtlGetRegion failed",0); 
    }

    if(DEBUG_COMM){
      printf("\n%d:issued get\n",armci_me);fflush(stdout);
    }

    ARMCI_PR_DBG("exit",0);
}

void armci_portals_get(int proc, void *src_buf, void *dst_buf, int bytes,
                             void** cptr,int tag)
{
int rc;
ptl_size_t offset_local = 0, offset_remote=0;
ptl_md_t *md_local;
ptl_handle_md_t *md_hdl_local;
int rem_info;
comp_desc *cdesc;
ptl_process_id_t dest_proc;
int c_info;
int cluster = armci_clus_id(proc);

    ARMCI_PR_DBG("enter",0);

    /*first remote process information*/
    /*dest_proc.nid = portals->procid_map[proc].nid;
    dest_proc.pid = portals->procid_map[proc].pid;*/
    dest_proc.nid = portals->servid_map[cluster].nid;
    dest_proc.pid = portals->servid_map[cluster].pid;

    /*create local xfer info*/
    cdesc = get_free_comp_desc(&c_info);
    md_local = &cdesc->mem_dsc;
    md_hdl_local = &cdesc->mem_dsc_hndl; 
    md_local->length=bytes;
    md_local->start=dst_buf;
    md_local->user_ptr = (void *)cdesc;
    md_local->options =  PTL_MD_OP_GET | PTL_MD_EVENT_START_DISABLE;

    /*get remote info*/
    rem_info = _get_rem_info(proc,src_buf,bytes,&offset_remote);

    cdesc->dest_id = proc;
    if (tag){
      *((comp_desc **)cptr) = cdesc;
      cdesc->tag = tag;
      cdesc->type = ARMCI_PORTALS_NBGET;
      /*printf("\n%d:get tag=%d c_info=%d
       * %p",armci_me,tag,c_info,cdesc);fflush(stdout);*/
    }
    else{
      cdesc->tag = 0;
      cdesc->type = ARMCI_PORTALS_GET; 
    }

    cdesc->active = 1;
    armci_client_direct_get(dest_proc,offset_remote,(ptl_match_bits_t)rem_info,
                bytes,md_local,md_hdl_local);

    if(!tag){ 
       armci_client_complete(0,proc,0,cdesc); /* check this later */
    }

    ARMCI_PR_DBG("exit",0);
}


void armci_client_nb_get(int proc, void *src_buf, int *src_stride_arr, 
                             void *dst_buf, int *dst_stride_arr, int bytes,
                             void** cptr,int tag)
{
}

void armci_client_direct_send(ptl_process_id_t dest_proc,
                ptl_size_t offset_remote, ptl_match_bits_t mb, size_t bytes,
                ptl_md_t *md_local, 
                ptl_handle_md_t *md_hdl_local)
{
int rc;
ptl_size_t offset_local = 0;

    ARMCI_PR_DBG("enter",0);

    if(DEBUG_COMM){
      printf("%d:armci_client_direct_send:BYTES = %d\n",armci_me,bytes);
      printf("\n%d:offr=%ld offl=%ld\n",armci_me,offset_remote,offset_local);
      fflush(stdout);
    }
    /*print_mem_desc(md_local);*/
    rc = PtlMDBind(portals->ni_h,*md_local, PTL_UNLINK, md_hdl_local);
    if (rc != PTL_OK){
      fprintf(stderr, "%d:PtlMDBind: %s\n", portals->rank, 
                      ARMCI_NET_ERRTOSTR(rc)); 
      armci_die("ptlmdbind send failed",0);
    }
    if(armci_must_remotecomplete){ 
    rc = PtlPutRegion(*md_hdl_local,offset_local,bytes,
                    PTL_ACK_REQ,
                    dest_proc,portals->ptl,0, mb,offset_remote, 0);
    }
    else{
    rc = PtlPutRegion(*md_hdl_local,offset_local,bytes,
                    PTL_NOACK_REQ,
                    dest_proc,portals->ptl,0, mb,offset_remote, 0);
    }

    if (rc != PTL_OK){
      fprintf(stderr, "%d:PtlPutRegion: %s\n", portals->rank, 
                      ARMCI_NET_ERRTOSTR(rc) );
      armci_die("PtlPutRegion failed",0);
    }

    ARMCI_PR_DBG("exit",0);
}


void armci_portals_put(int proc, void *src_buf, void *dst_buf, int bytes,
                             void** cptr,int tag)
{
int rc;
ptl_size_t offset_local = 0, offset_remote=0;
ptl_md_t *md_local;
ptl_handle_md_t *md_hdl_local;
int rem_info;
comp_desc *cdesc;
ptl_process_id_t dest_proc;
int c_info;
int cluster = armci_clus_id(proc);

    ARMCI_PR_DBG("enter",0);

    /*first process information*/
    dest_proc.nid = portals->servid_map[cluster].nid;
    dest_proc.pid = portals->servid_map[cluster].pid;
    /*dest_proc.nid = portals->procid_map[proc].nid;
    dest_proc.pid = portals->procid_map[proc].pid;*/

    /*create local xfer info*/
    cdesc = get_free_comp_desc(&c_info);
    md_local = &cdesc->mem_dsc;
    md_hdl_local = &cdesc->mem_dsc_hndl; 
    md_local->length=bytes;
    md_local->start=src_buf;
    md_local->user_ptr = (void *)cdesc;
    md_local->options =  PTL_MD_OP_PUT | PTL_MD_EVENT_START_DISABLE;
    
    /*get remote info*/
    rem_info = _get_rem_info(proc,dst_buf,bytes,&offset_remote);
                    

    if(DEBUG_COMM){
      printf("\n%d:offr=%ld offl=%ld\n",armci_me,offset_remote,offset_local);
    }

    cdesc->dest_id = proc;
    if (tag){
      *((comp_desc **)cptr) = cdesc;
      cdesc->tag = tag;
      cdesc->type = ARMCI_PORTALS_NBPUT;
      /*printf("\n%d:put tag=%d c_info=%d
       * %p",armci_me,tag,c_info,cdesc);fflush(stdout);*/
    }
    else{
      cdesc->tag = 0;
      cdesc->type = ARMCI_PORTALS_PUT; 
    }
    
    cdesc->active = 1;

    armci_client_direct_send(dest_proc,offset_remote,(ptl_match_bits_t)rem_info,
                bytes,md_local,md_hdl_local);

    if(!tag){ 
       armci_client_complete(0,proc,0,cdesc); /* check this later */
    }
    else
       portals->outstanding_puts++;


    ARMCI_PR_DBG("exit",0);

}

void armci_client_nb_send(int proc, void *src_buf, int *src_stride_arr, 
                             void *dst_buf, int *dst_stride_arr, int bytes,
                             void** cptr,int tag)
                             
{
}

/*using non-blocking for multiple 1ds inside a 2d*/
void armci_network_strided(int op, void* scale, int proc,void *src_ptr,
                int src_stride_arr[], void* dst_ptr, int dst_stride_arr[],
                int count[], int stride_levels, armci_ihdl_t nb_handle)
{
int i, j,tag=0;
long idxs,idxd;    /* index offset of current block position to ptr */
int n1dim;  /* number of 1 dim block */
int bvalue_s[MAX_STRIDE_LEVEL], bunit[MAX_STRIDE_LEVEL];
int bvalue_d[MAX_STRIDE_LEVEL];
int bytes = count[0];
void *sptr,*dptr;
NB_CMPL_T cptr;
ptl_process_id_t dest_proc;
ptl_size_t offset_remote;
comp_desc *cdesc;
int c_info; 
ptl_md_t *md_local;
int rem_info;
int cluster = armci_clus_id(proc);

    ARMCI_PR_DBG("enter",0);

       printf("%s calling abort ... network_strided not implemented\n",Portals_ID());
       abort();

    if(nb_handle)tag=nb_handle->tag;

    /*first remote process information*/
    dest_proc.nid = portals->servid_map[cluster].nid;
    dest_proc.pid = portals->servid_map[cluster].pid;
    /*dest_proc.nid = portals->procid_map[proc].nid;
    dest_proc.pid = portals->procid_map[proc].pid;*/

    rem_info = _get_rem_info(proc,(op==GET)?src_ptr:dst_ptr,bytes,&offset_remote);

    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++)
        n1dim *= count[i];

    /* calculate the destination indices */
    bvalue_s[0] = 0; bvalue_s[1] = 0; bunit[0] = 1; 
    bvalue_d[0] = 0; bvalue_d[1] = 0; bunit[1] = 1;
    for(i=2; i<=stride_levels; i++) {
        bvalue_s[i] = bvalue_d[i] = 0;
        bunit[i] = bunit[i-1] * count[i-1];
    }

    if(ARMCI_ACC(op)){ /*for now die for acc*/
      /*lock here*/
#          ifdef ARMCI_CHECK_STATE
            arlist_print(arn);
            armci_rem_state(armci_clus_info[proc].master%armci_clus_info[0].nslave);
#          endif
      printf("\nSHOULD NOT DO NETWORK_STRIDED FOR ACCS \n",armci_me);
      fflush(stdout);
      armci_die("network_strided called for acc",proc);
    }

    /*loop over #contig chunks*/
    for(i=0; i<n1dim; i++) {
    ptl_handle_md_t *md_hdl_local;
      tag = GET_NEXT_NBTAG();      
      idxs = 0;
      idxd = 0;
      for(j=1; j<=stride_levels; j++) {
        idxs += bvalue_s[j] * src_stride_arr[j-1];
        idxd += bvalue_d[j] * dst_stride_arr[j-1];
        if((i+1) % bunit[j] == 0) {bvalue_s[j]++;bvalue_d[j]++;}
        if(bvalue_s[j] > (count[j]-1)) bvalue_s[j] = 0;
        if(bvalue_d[j] > (count[j]-1)) bvalue_d[j] = 0;
      }
      sptr = ((char *)src_ptr)+idxs;
      dptr = ((char *)dst_ptr)+idxd;
      cdesc = get_free_comp_desc(&c_info);
      md_local = &cdesc->mem_dsc;
      md_hdl_local = &cdesc->mem_dsc_hndl;
      md_local->length=bytes;
      md_local->start=(op==GET)?dptr:sptr;
      md_local->user_ptr = (void *)cdesc;
      cdesc->dest_id = proc;
      cdesc->tag = tag;
      
      if(op==GET){
        md_local->options =  PTL_MD_OP_GET | PTL_MD_EVENT_START_DISABLE;
        cdesc->active = 1;
        cdesc->type = ARMCI_PORTALS_NBGET;
        /*
        printf("\n%d:reminfo=%d off=%d idxs=%d idxd=%d",armci_me, rem_info,
                        offset_remote, idxs, idxd);
                        */
        armci_client_direct_get( dest_proc,offset_remote+idxs,rem_info,
                        bytes,md_local,md_hdl_local);
      }
      else if(op==PUT){
        cdesc->active = 1;
        cdesc->type = ARMCI_PORTALS_NBPUT;
        md_local->options =  PTL_MD_OP_PUT | PTL_MD_EVENT_START_DISABLE;
        armci_client_direct_send(dest_proc,offset_remote+idxd,rem_info,
                bytes,md_local,md_hdl_local);
        if(op==PUT)portals->outstanding_puts++;
      }
      else if(ARMCI_ACC(op)){
        assert(0);
      }
      else{
        ARMCI_PR_DBG("exit",0);
        armci_die("in network_strided unknown opcode",op);
      }
      armci_client_complete(0,proc,tag,cdesc);
    }

    if(ARMCI_ACC(op)){
    /*unlock here*/
    }

    if(nb_handle){
      /* completing the last call is sufficient, given ordering semantics*/
      nb_handle->tag=tag;
      nb_handle->cmpl_info=cdesc;
    }
    else{
      /*completing the last call ensures everything before it is complete this
       * is one of the main reasons why dataserver is necessary*/
      /*armci_client_complete(0,proc,tag,cdesc);*/
    }
    ARMCI_PR_DBG("exit",0);
}

void armci_client_direct_getput(ptl_process_id_t dest_proc,
                ptl_size_t offset_remote, ptl_match_bits_t mb, size_t bytes,
                ptl_md_t *md_local_get,ptl_md_t *md_local_put, 
                ptl_handle_md_t *md_hdl_local_get, ptl_handle_md_t
                *md_hdl_local_put)
{
int rc;
ptl_size_t offset_get = 0;
ptl_size_t offset_put = 0;

    ARMCI_PR_DBG("enter",0);

    if(DEBUG_COMM){
      printf("%d:armci_client_direct_getput:BYTES = %d\n",armci_me,bytes);
      printf("\n%d:offr=%ld\n",armci_me,offset_remote);fflush(stdout);
    }

    rc = PtlGetPutRegion(*md_hdl_local_get, offset_get, *md_hdl_local_put,
                    offset_put,bytes,dest_proc, portals->ptl,0,mb,
                    offset_remote,0);
    if (rc != PTL_OK){
      printf("%d:PtlGetPutRegion: %s\n", portals->rank, ARMCI_NET_ERRTOSTR(rc) );
      armci_die("PtlGetPutRegion failed",0);
    }
    
    ARMCI_PR_DBG("exit",0);

}


long a_p_putfrom;
long a_p_getinto;


int armci_portals_rmw_(int op, int *ploc, int *prem, int extra, int proc)
{
        printf("error rmw");
    return(0);
}

void armci_portals_shmalloc_allocate_mem(int num_lks)
{
void **ptr_arr;
void *ptr;
armci_size_t bytes = 128;
int i;    
    
    ARMCI_PR_DBG("enter",0);
    ptr_arr    = (void**)malloc(armci_nproc*sizeof(void*));
    if(!ptr_arr) armci_die("armci_shmalloc_get_offsets: malloc failed", 0);
    bzero((char*)ptr_arr,armci_nproc*sizeof(void*));

    PARMCI_Malloc(ptr_arr,bytes);
    ARMCI_PR_DBG("exit",0);
    
    return;
}


void armci_wait_for_server()
{
    ARMCI_PR_DBG("enter",0);
    armci_server_terminating=1;
    armci_serv_quit();
    ARMCI_PR_DBG("exit",0);
}

/*client buffers info*/
void armci_portals_client_buf_info(char *buf, ptl_match_bits_t *mb, ptl_size_t *offset,int proc)
{
    ARMCI_PR_DBG("enter",0);
    *mb = (1<<30); 
    *offset = buf-client_buf_ptrs[proc];
    if(DEBUG_SERV){printf("\n(%d):serv writing to ofset %d on %d\n",armci_me,*offset,proc);fflush(stdout);}
    ARMCI_PR_DBG("exit",0);
}

/*memory for client buffers*/
char *armci_portals_client_buf_allocate(int bytes)
{
void *ptr;
ptl_match_bits_t ignbits = 0xFFFFFFFF0FFFFFFF;
ptl_match_bits_t mbits = 1;
ptl_md_t *md_ptr,md;
ptl_process_id_t match_id;
ptl_handle_me_t me_h;
ptl_handle_md_t md_h;
int rc;
    ARMCI_PR_DBG("enter",sizeof(ptl_match_bits_t));
    ptr = malloc(bytes);
    bzero(ptr,bytes);
    assert(ptr);

    mbits = (1<<30);
    md_ptr            = &(md);
    md_ptr->start     = ptr;
    md_ptr->length    = bytes;
    md_ptr->threshold = PTL_MD_THRESH_INF;
    md_ptr->options   =  PTL_MD_OP_PUT | PTL_MD_OP_GET | PTL_MD_MANAGE_REMOTE | PTL_MD_EVENT_START_DISABLE;
    md_ptr->user_ptr  = NULL;
    md_ptr->max_size  = 0;
    /*logic that says, eq_h is now recieving data for the buffers, including acks! */
    md_ptr->eq_handle = portals->eq_h;
    match_id.nid = PTL_NID_ANY;
    match_id.pid = PTL_PID_ANY; 
    rc = PtlMEAttach(portals->ni_h,portals->ptl,match_id,
                    mbits,ignbits,PTL_RETAIN,PTL_INS_AFTER,&(me_h)); 
    if (rc != PTL_OK){
      printf("%d:PtlMEAttach: %s\n", portals->rank, ARMCI_NET_ERRTOSTR(rc) );
      armci_die("portals attach error2",rc);
    }
    rc = PtlMDAttach(me_h,md,PTL_RETAIN,&md_h);
    if (rc != PTL_OK) {
      printf("%d:PtlMDAttach: %s %d\n", portals->rank, ARMCI_NET_ERRTOSTR(rc),(client_md_count+serv_md_count) );
      armci_die("portals attach error CBA",rc);
    }
    client_md_count++;

    client_buf_ptrs[armci_me]=ptr;
    armci_msg_barrier();
    armci_exchange_address(client_buf_ptrs,armci_nproc);

    ARMCI_PR_DBG("exit",0);
    return(ptr);
}

void armci_transport_cleanup()
{
    /*for i=0tomaxpendingclean*/
    ARMCI_PR_DBG("enter",0);
    free(client_buf_ptrs);
    ARMCI_PR_DBG("exit",0);
}

void free_serv_bufs()
{
    if(serv_bufs) free(serv_bufs);
}


int armci_send_req_msg(int proc, void *buf, int bytes,int tag)
{
#ifndef OLD_PORTALS_CODE
        int cluster = armci_clus_id(proc);
        int serv    = armci_clus_info[cluster].master;
        char *buffer = NULL;
        request_header_t *msginfo = (request_header_t *) buf;

//    # ifdef PORTALS_LIMIT_REMOTE_REQUESTS_BY_NODE
        _armci_buf_ensure_one_outstanding_op_per_node(buf,cluster);
//    # endif

      # ifdef PORTALS_USE_ARMCI_CLIENT_BUFFERS
        BUF_INFO_T *bufinfo=_armci_buf_to_bufinfo(msginfo);
        _buf_ackresp_t *ar = &bufinfo->ar;
        portals_ds_req_t *req = &ar->req;
      # endif

        if(msginfo->operation == PUT || ARMCI_ACC(msginfo->operation)) {
        // printf("%s cp: sending packed put\n",Portals_ID());
         # ifdef PORTALS_USE_ARMCI_CLIENT_BUFFERS
           portals_remote_nbput(buf, buf, cluster, req);
        // portalsWaitOnRequest(req);
         # else
           portals_remote_put(buf, buf, cluster);
         # endif
        // printf("%s cp: finished packed put\n",Portals_ID());
        }

        else if(msginfo->operation == GET) {
           buffer = (char *) buf;
           buffer += sizeof(request_header_t);
           buffer += msginfo->dscrlen;
        // printf("%s cp: sending blocking get request\n",Portals_ID());
         # ifdef PORTALS_USE_ARMCI_CLIENT_BUFFERS
           portals_remote_nbget(buffer, msginfo, cluster, req);
        // portalsWaitOnRequest(req);
         # else
           portals_remote_get(buffer, msginfo, cluster); 
         # endif
        // printf("%s cp: get request finished\n",Portals_ID());
        }

        else if(msginfo->operation == ACK) {
         # ifdef PORTALS_USE_ARMCI_CLIENT_BUFFERS
           portalsRemoteOperationToNode(buf, bytes, cluster, req);
        // portalsWaitOnRequest(req);
         # else
           portalsBlockingRemoteOperationToNode(buf, bytes, cluster);
         # endif
        }

        else if(msginfo->operation == ARMCI_SWAP || msginfo->operation == ARMCI_SWAP_LONG ||
                msginfo->operation == ARMCI_FETCH_AND_ADD || msginfo->operation == ARMCI_FETCH_AND_ADD_LONG) {
           buffer = (char *) buf;
           buffer += sizeof(request_header_t);
           buffer += msginfo->dscrlen;
           portals_remote_rmw(buffer, msginfo, cluster, req);
         # ifndef PORTALS_USE_ARMCI_CLIENT_BUFFERS
           portalsWaitOnOperation(req);
         # endif 
        }

        else {
           printf("%s cp: msginfo->operation=%d not supported yet\n",Portals_ID(),msginfo->operation);
           abort();
        }

      # ifdef PORTALS_USE_ARMCI_CLIENT_BUFFERS
/* for now, clear the ackresp structure because the call had to have been blocking 
   later, we will allow a modified x_net_wait_ackresp clear it */
        ar->val = ar->valc = 0;
        if(ar==_buf_ackresp_first)_buf_ackresp_first=ar->next;
        if(ar->next!=NULL){
          ar->next->previous=ar->previous;
        }
        if(ar->previous!=NULL){
          ar->previous->next=ar->next;
          if(_buf_ackresp_cur==ar)_buf_ackresp_cur=ar->previous;
        }
        if(_buf_ackresp_cur==ar)_buf_ackresp_cur=NULL;
        ar->previous=ar->next=NULL;
      # endif

        return 0;

#else

int rc;
ptl_size_t offset_local = 0, offset_remote=0;
ptl_md_t *md_local;
ptl_handle_md_t *md_hdl_local;
int rem_info;
comp_desc *cdesc;
void *cptr;
ptl_process_id_t dest_proc;
int c_info;
int cluster = armci_clus_id(proc);
int serv = armci_clus_info[cluster].master;
request_header_t *msginfo = (request_header_t *)buf;

    ARMCI_PR_DBG("enter",0);
    if(msginfo->operation==GET && msginfo->dscrlen<=msginfo->datalen){
      *(long *)((char *)(msginfo+1)+msginfo->datalen)=0;
    }
     
    /*badbadbad*/
    msginfo->tag.ack_ptr=&(msginfo->tag.ack);
    cptr = (void *)((double *)buf-1);
    /*first process information*/
    dest_proc.nid = portals->servid_map[cluster].nid;
    dest_proc.pid = portals->servid_map[cluster].pid;
    /*create local xfer info*/
    cdesc = get_free_comp_desc(&c_info);
    md_local = &cdesc->mem_dsc;
    md_hdl_local = &cdesc->mem_dsc_hndl; 
    md_local->length=bytes;
    md_local->start=buf;
    md_local->user_ptr = (void *)cdesc;
    md_local->options =  PTL_MD_OP_PUT | PTL_MD_EVENT_START_DISABLE;
    
    /*get remote info*/
    rem_info = _get_rem_servinfo(cluster,(size_t)bytes,&offset_remote);
                    
    if(DEBUG_COMM){
      printf("\n%d:offr=%ld offl=%ld\n",armci_me,offset_remote,offset_local);
    }

    cdesc->dest_id = serv;
    *((comp_desc **)cptr) = cdesc;
    if(tag==0)tag=GET_NEXT_NBTAG();
    cdesc->tag = tag;
    cdesc->type = ARMCI_PORTALS_NBPUT;
    /*printf("\n%d:put tag=%d c_info=%d
    * %p",armci_me,tag,c_info,cdesc);fflush(stdout);*/
    cdesc->active = 1;

    if(msginfo->operation==PUT || msginfo->operation == UNLOCK || ARMCI_ACC(msginfo->operation)){
      _buf_ackresp_cur->valc = _buf_ackresp_cur->val = (char *)msginfo->tag.ack_ptr-client_buf_ptrs[armci_me];
#    ifdef ARMCI_CHECK_STATE
      arlist_add(&arn,_buf_ackresp_cur->val,msginfo->operation);
#    endif
    }
    else {
      _buf_ackresp_cur->valc = _buf_ackresp_cur->val = (char *)msginfo->tag.data_ptr-client_buf_ptrs[armci_me];
#    ifdef ARMCI_CHECK_STATE
      arlist_add(&arn,_buf_ackresp_cur->val,msginfo->operation);
#    endif
    }

    if(DEBUG_COMM){printf("\n%d:registered %d in val to %d at %d %d\n",armci_me,_buf_ackresp_cur->val,serv,offset_remote,msginfo->operation);fflush(stdout);}
    _armci_buf_ensure_pend_outstanding_op_per_node(buf,cluster);
    armci_client_direct_send(dest_proc,offset_remote,(ptl_match_bits_t)rem_info,
                bytes,md_local,md_hdl_local);
    /*if(msginfo->operation==GET){
       BUF_INFO_T *info=((char *)msginfo-sizeof(BUF_EXTRA_FIELD_T)-sizeof(BUF_INFO_T));
       armci_client_complete(0,proc,cdesc->tag,cdesc); 
    }*/
    /*armci_client_complete(0,proc,cdesc->tag,cdesc);*/ 

    portals->outstanding_puts++;

    ARMCI_PR_DBG("exit",0);
    return 0;
#endif

}


char *armci_ReadFromDirect(int proc, request_header_t *msginfo, int len)
{
#ifndef OLD_PORTALS_CODE
      # ifdef PORTALS_USE_ARMCI_CLIENT_BUFFERS
        BUF_INFO_T *bufinfo = _armci_buf_to_bufinfo(msginfo);
        portals_ds_req_t *req = &bufinfo->ar.req;
        portalsWaitOnRequest(req);
      # endif
        char *ret = (char *) msginfo;
        ret += sizeof(request_header_t);
        ret += msginfo->dscrlen;
        return ret;
#else
long *flag;
int loop;
BUF_INFO_T *bufinfo=_armci_buf_to_bufinfo(msginfo);

    ARMCI_PR_DBG("enter",0);
    if(len)
      flag = (long *)((char *)(msginfo+1)+len);
    else
      flag = (long *)((char *)(msginfo+1)+msginfo->datalen);
    x_net_wait_ackresp(&(bufinfo->ar));
                   
    while(armci_util_long_getval(flag)  != ARMCI_TAIL){
      loop++;
      loop %=100000;
      if(loop==0){
        if(DEBUG_COMM){
          printf("%d: client flag(%p)=%ld off=%d %d\n",
                          armci_me,flag,*flag,msginfo->datalen,*((int*)(msginfo+1)));
          fflush(stdout);
        }
      }
    }
    *flag=0;
    ARMCI_PR_DBG("exit",0);
    return (msginfo+1);
#endif
}

#ifdef ARMCI_CHECK_STATE
extern void sarlist_add(int,int,long);
#endif

void armci_WriteToDirect(int proc, request_header_t* msginfo, void *buf)
{
#ifndef OLD_PORTALS_CODE
        ptl_size_t bytes = (ptl_size_t) msginfo->datalen;
        ptl_event_t *ev = (ptl_event_t *) msginfo->tag.user_ptr;
        portals_ds_send_put(buf, msginfo->datalen, ev->initiator, ev->hdr_data);
        // you could do an assertion that the portals_id_map of proc == ev->initiator
#else
long *tail;
int bytes;
void *dst_addr = msginfo->tag.data_ptr;
ptl_match_bits_t ignbits = 0xFFFFFFFF0FFFFFFF;
ptl_match_bits_t mbits = 1;
ptl_md_t *md_ptr,md;
ptl_process_id_t match_id;
ptl_handle_me_t me_h;
ptl_size_t offst,localoffset;
int rc;
   
    /* set tail ack, make sure it is alligned */
    ARMCI_PR_SDBG("enter",0);
    bytes = msginfo->datalen+sizeof(long);
    if(!(buf>=serv_bufs->buf && buf<serv_bufs->bufend)){
      bcopy(buf,(msginfo+1),bytes);
      buf=(msginfo+1);
    }
    tail = (long*)(buf + msginfo->datalen);
    *tail = ARMCI_TAIL;
 
    armci_portals_client_buf_info((char *)dst_addr,&mbits,&offst,proc); 

#  ifdef ARMCI_CHECK_STATE
    sarlist_add(proc,msginfo->operation,offst);
#  endif

    match_id.nid = portals->procid_map[proc].nid;
    match_id.pid = portals->procid_map[proc].pid;
    localoffset=(char *)buf-(char *)serv_bufs->buf;
    if(DEBUG_COMM){
    printf("\n(%d):dst=%p,mbits=%d,localoffset=%d,offst=%d,proc=%d,nid=%d,pid=%d len=%d\n",armci_me,
                    dst_addr,mbits,localoffset,offst,proc,portals->procid_map[proc].nid,
                    portals->procid_map[proc].pid,bytes);fflush(stdout);
    }
    rc = PtlPutRegion(serv_response_md_h,localoffset,bytes,PTL_NOACK_REQ,
                    match_id,portals->ptl,0,mbits,offst,0);
    if (rc != PTL_OK){
      fprintf(stderr, "%d:PtlPutRegion: %s\n", portals->Srank, 
                      ARMCI_NET_ERRTOSTR(rc) );
      armci_die("PtlPutRegion failed",0);
    }
    ARMCI_PR_SDBG("exit",0);
#endif
}


void armci_rcv_req(void *mesg,void *phdr,void *pdescr,void *pdata,int *buflen)
{
int i,na;
char *a;
double *tmp;
request_header_t *msginfo = (request_header_t *)mesg;
    ARMCI_PR_SDBG("enter",msginfo->operation);
    *(void **)phdr = msginfo;
    if(0){
        printf("%s [ds]: got %d req (hdrlen=%d dscrlen=%d datalen=%d %d) from %d\n",
               Portals_ID(), msginfo->operation, sizeof(request_header_t), msginfo->dscrlen,
               msginfo->datalen, msginfo->bytes,msginfo->from);
               fflush(stdout);
    }
    /* we leave room for msginfo on the client side */
    *buflen = MSG_BUFLEN - sizeof(request_header_t);

    if(msginfo->bytes) {
      *(void **)pdescr = msginfo+1;
      *(void **)pdata = msginfo->dscrlen + (char*)(msginfo+1);

      if(msginfo->operation == GET) {
         // the descriptor will exists after the request header
         // but there will be no data buffer
         // use the MessageRcvBuffer
         *(void**) pdata = MessageSndBuffer;
//       printf("%s (server) overriding pdata in rcv_req\n",Portals_ID());
      }
      printf("%s [ds] oper=%d; bytes=%d\n",armci_me,msginfo->operation,msginfo->bytes);
    }
    else {
      printf("%s [ds] bytes=%d\n",armci_me,msginfo->bytes);
      *(void**)pdescr = NULL;
      *(void**)pdata = MessageRcvBuffer;
    }
    ARMCI_PR_SDBG("exit",msginfo->operation);
}

void armci_call_data_server()
{
int rc;  
ptl_event_t ev_t;
ptl_event_t *ev=&ev_t;
serv_buf_t *compbuf = NULL;
int loop=1;
int temp_proc;
int ccc=2,rrr;
cpu_set_t mycpuid,new_mask;
char str[CPU_SETSIZE];
char ncid[8],*cidptr,cid[8];
extern char * cpuset_to_cstr(cpu_set_t *mask, char *str);
    ARMCI_PR_SDBG("enter",0);
    //if(armci_me==0)unsetenv("CRAY_PORTALS_USE_BLOCKING_POLL");
    sprintf (cid, "%d", ccc);
    rrr = cstr_to_cpuset(&new_mask,cid);

/* ------------------------------------------------------------ *\
   Change affinity for the data server
\* ------------------------------------------------------------ */
    if(sched_setaffinity(0, sizeof (new_mask), &new_mask)) {
      perror("sched_setaffinity");
      printf("failed to set pid %d's affinity.\n", getpid());
    }
    if(DEBUG_SERV){
      rrr=sched_getaffinity(0, sizeof(mycpuid), &mycpuid);
      if(rrr)perror("sched_getaffinity");
      cidptr = cpuset_to_cstr(&mycpuid,ncid);
      printf("(%d):my affinity is to %s\n",armci_me,ncid);
      fflush(stdout);
    }

/* ------------------------------------- *\
   Main data server loop
\* ------------------------------------- */
   while(armci_server_terminating==0){

   /* ------------------------------------------------------------ *\
      check event queue for incoming data requests from remote CPs
   \* ------------------------------------------------------------ */
      ev->type=0;
      if((rc = PtlEQWait(portals->Seq_h, ev)) != PTL_OK){
        printf("(%d):PtlEQWait(): %d %s\n", armci_me,rc,ARMCI_NET_ERRTOSTR(rc) );
        armci_die("EQWait problem",rc);
      }
      if (ev->ni_fail_type != PTL_NI_OK) {
        printf("(%d)%d,%d:NI sent %d in event.\n",
          armci_me,portals->Srank.nid, portals->Srank.pid,ev->ni_fail_type);
        fflush(stdout);
        armci_die2("event failure problem",ev->initiator.nid,ev->initiator.pid);
      }
      if(DEBUG_SERV){
        printf("\n(%d):armci_call_data_server: ptl event detected=%d\n",armci_me,ev->type);
        fflush(stdout);
      }

   /* ------------------------------------------------------------ *\
      PTL_EVENT_SEND_END: is ignored.  This event is triggered as
      the DS returns data to a remote CP via a PtlPut.  This event
      signals that that PtlPut has complete.
   \* ------------------------------------------------------------ */
      if(ev->type == PTL_EVENT_SEND_END) continue;


   /* ------------------------------------------------------------ *\
      PTL_EVENT_PUT_END: this is the key portals event for the DS.
      PUT_END signifies that a remote data request has come in
      from a remote CP.  This data request will be handled by the
      data server: armci_data_server
   \* ------------------------------------------------------------ */
      else if(ev->type == PTL_EVENT_PUT_END) {
         if(DEBUG_SERV) {
            printf("\n(%d):ev->offset=%d from %d%d",armci_me,ev->offset,
                                                    ev->initiator.pid,ev->initiator.nid);
            fflush(stdout);
         }
         armci_data_server(((char *)serv_bufs->buf+ev->offset));
      }

   /* ------------------------------------------------------------ *\
      Unexpected Portals Event -- Panic!
   \* ------------------------------------------------------------ */
      else { 
         armci_die("unexpected event in data server",ev->type);
      }
   }
   ARMCI_PR_SDBG("exit",0);
}

void x_buf_wait_ack(request_header_t *msginfo, BUF_INFO_T *bufinfo)
{
  ARMCI_PR_DBG("enter",bufinfo->ar.val);
  if(DEBUG_COMM){printf("\n%d:waiting for ack at %p",armci_me,&(msginfo->tag.ack));fflush(stdout);}
  x_net_wait_ackresp(&(bufinfo->ar));
  armci_util_wait_long(&(msginfo->tag.ack),ARMCI_STAMP,10000);
  if(DEBUG_COMM){printf("\n%d:done waiting for ack at %p",armci_me,&(msginfo->tag.ack));fflush(stdout);}
  msginfo->tag.ack=0;
  ARMCI_PR_DBG("exit",0);
}

void x_net_send_ack(request_header_t *msginfo, int proc,void *dst,void *src)
{
long *tail;
int bytes=sizeof(long);
ptl_size_t offst;
ptl_match_bits_t ignbits = 0xFFFFFFFF0FFFFFFF;
ptl_match_bits_t mbits = 1;
ptl_process_id_t match_id;
int rc;
   
    /* set tail ack, make sure it is alligned */
    ARMCI_PR_SDBG("enter",0);

 
    armci_portals_client_buf_info((char *)dst,&mbits,&offst,proc); 

#  ifdef ARMCI_CHECK_STATE
    sarlist_add(proc,msginfo->operation,offst);
#  endif

    match_id.nid = portals->procid_map[proc].nid;
    match_id.pid = portals->procid_map[proc].pid;
    if(DEBUG_SERV){
      printf("\n(%d):dst=%p,mbits=%d,offst=%d,proc=%d,nid=%d,pid=%d len=%d\n",armci_me,
                    dst,mbits,offst,proc,portals->procid_map[proc].nid,
                    portals->procid_map[proc].pid,bytes);fflush(stdout);
    }

    rc = PtlPutRegion(serv_ack_md_h,0,bytes,PTL_NOACK_REQ,
                    match_id,portals->ptl,0,mbits,offst,0);
    if (rc != PTL_OK){
      fprintf(stderr, "%d:PtlPutRegion: %s\n", portals->Srank, 
                      ARMCI_NET_ERRTOSTR(rc) );
      armci_die("PtlPutRegion failed",0);
    }
    ARMCI_PR_SDBG("exit",0);
}

long x_net_offset(char *buf,int proc)
{
#ifdef ARMCI_REGISTER_SHMEM
int i;
#if NEW_MALLOC
    if((i=check_meminfo(buf,1,proc))==0)
      armci_die("x_net_offset,reg not found",proc);
    return(all_meminfo[proc].serv_offs[i]);
#else
    ARMCI_PR_DBG("enter",_rem_meminfo[proc].reg_count);
    if(DEBUG_COMM){printf("\n%d:%s:buf=%p",armci_me,__FUNCTION__,buf);fflush(stdout); }
    for(i=0;i<_rem_meminfo[proc].reg_count;i++){
      if(IN_REGION(buf,_rem_meminfo[proc].reginfo[i])){
#ifdef DEBUG_MEM
        {printf("\n%d:found it in reg=%d  (%p,%d) for proc=%d",armci_me,i,_rem_meminfo[proc].reginfo[i].base_ptr,_rem_meminfo[proc].reginfo[i].size,proc);}
#endif
        return((long)((char *)_rem_meminfo[proc].reginfo[i].serv_ptr-(char *)_rem_meminfo[proc].reginfo[i].base_ptr));
      }
    }
#endif
    ARMCI_PR_DBG("exit",0);
#else
 printf("x_net_offset called; this shouldn't happen ...\n"); abort();
#endif
}

void armci_set_serv_mutex_arr(void *ptr)
{
int i;
long offset;
    ARMCI_PR_DBG("enter",0);
    offset=x_net_offset(ptr,armci_me);

    _armci_server_mutex_ready=1;
    _armci_server_mutex_ptr = (char *)ptr+offset;
    ARMCI_PR_DBG("exit",0);

}

