/* Author: Abhinav Vishnu */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <mpi.h>

#include <assert.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include <infiniband/verbs.h>

#include "comex.h"
#include "comex_impl.h"
#include "openib.h"
#include "reg_cache.h"

/* Global vars */
int me;
int nprocs;

double timecalc[20000];
struct RC_Conn      conn;

struct User_Opts    opts;
struct HCA          hca;
struct Local_Buf    lbuf;
struct Remote_Buf   rbuf;

struct ibv_recv_wr  rr_desc;
struct ibv_sge      rr_sg_entry;

pthread_t           test_async_thread;

void openib_init_envs()
{
    char *value;
    if ((value = getenv("COMEX_OPENIB_USE_DREG")) != NULL){
        l_state.comex_openib_use_dreg = (atoi(value));
    }
    else {
        l_state.comex_openib_use_dreg = 0;
    }


}

void init_params(void)
{
    opts.msg_size = DEFAULT_MSG_SIZE;
    opts.iterations = DEFAULT_ITERATIONS;

    opts.rdma_flag = 0;
    opts.sendrecv_flag = 0;
    opts.all_msgs = 1;
    opts.memperf = 0;	
    opts.align = ALIGNMENT;

    opts.ib_devname = NULL;
    opts.num_cqe = DEFAULT_NUM_CQE;
    opts.send_wr = DEFAULT_NUM_OUST_SEND;
    opts.recv_wr = DEFAULT_NUM_OUST_RECV;
    opts.num_sge = DEFAULT_NUM_SGE;
    opts.latency = 0;
    opts.bandwidth = 0;
    opts.bibandwidth = 0;
    opts.window = DEFAULT_WINDOW_SIZE;
    opts.num_ports = 1;
    opts.striping_threshold = STRIPING_THRESHOLD;
    opts.num_qp_per_port = 1;
    opts.default_port = DEFAULT_PORT;
    hca.ib_dev = NULL;
    hca.context = NULL;
    hca.pd = NULL;

    lbuf.buf_original = NULL;
    lbuf.buf = NULL;
    lbuf.mr = NULL;

    conn.qp = (struct ibv_qp **) malloc(nprocs * sizeof(struct ibv_qp *));
    conn.lid = (uint16_t *) malloc(nprocs * sizeof(uint16_t));
    conn.qp_num = (uint32_t *) malloc(nprocs * sizeof(uint32_t));

    rbuf.qp_num = (uint32_t *) malloc(nprocs * sizeof(uint32_t));
    rbuf.lid = (uint16_t *) malloc(nprocs * sizeof(uint16_t));
    rbuf.rkey = (uint32_t *) malloc(nprocs * sizeof(uint32_t));
    rbuf.buf = (char **) malloc(nprocs * sizeof(char *));

    assert(conn.qp && conn.lid && rbuf.qp_num && rbuf.lid && rbuf.rkey && rbuf.buf);
}


double my_wtime()
{
    double t;

    struct timeval tv;
    static unsigned long initialized = 0;
    static unsigned long  sec_base;
    gettimeofday(&tv, NULL);
    if (!initialized) {
        sec_base = tv.tv_sec;
        initialized = 1;
    }
    t = (double) (tv.tv_sec - sec_base) + 
        (double) tv.tv_usec * 1.0e-6;

    return t;
}

void release_resources(void)
{
    int i;
   
    // Destroy the registration cache

    reg_cache_destroy(nprocs);

    for(i = 0; i < nprocs ; i++) {
        if (me != i) {
            if(ibv_destroy_qp(conn.qp[i])) {
                printf("Exiting\n");
                exit(1);
            }
        }
    }
#if 0
    if(ibv_destroy_cq(hca.cq)) {
        fprintf(stderr,"Couldn't destroy cq %s\n",
                ibv_get_device_name(hca.ib_dev));
    }
    if(ibv_dealloc_pd(hca.pd)) {
        fprintf(stderr,"Couldn't free pd %s\n",
                ibv_get_device_name(hca.ib_dev));
    }
#endif

    free(conn.qp);
    free(conn.lid);
    free(conn.qp_num);

    free(rbuf.qp_num);
    free(rbuf.lid);
#if 0
    free(rbuf.rkey); 
    free(rbuf.buf);
#endif
}

int open_hca(void)
{
    struct ibv_device **dev_list=NULL;
    int num_hcas;

    dev_list = ibv_get_device_list(&num_hcas);


    // Assume that the first device has an ACTIVE port
    // if it does not we do not handle this situation for now
    
    hca.ib_dev = dev_list[0];

    hca.context = ibv_open_device(hca.ib_dev);

    if(!hca.context) {
        fprintf(stderr,"Couldn't get context %s\n",
                ibv_get_device_name(hca.ib_dev));
        return 1;
    }

    hca.pd = ibv_alloc_pd(hca.context);

    assert(hca.pd != NULL);
    
    if(!hca.pd) {
        fprintf(stderr,"Couldn't get pd %s\n",
                ibv_get_device_name(hca.ib_dev));
        return 1;
    }

    return 0;
}


int create_cq(void)
{
    hca.cq = ibv_create_cq(hca.context, opts.num_cqe, NULL, NULL, 0);
    assert(hca.cq);    

    return 0;
}

int exch_addr(void)
{
	int rc;

    rc = MPI_Alltoall((void *)conn.qp_num, sizeof(uint32_t), MPI_BYTE, 
            (void *)rbuf.qp_num, sizeof(uint32_t), MPI_BYTE, l_state.world_comm);
   
    assert(!rc); 
    rc = MPI_Alltoall((void *)conn.lid, sizeof(uint16_t), MPI_BYTE, 
            (void *)rbuf.lid, sizeof(uint16_t), MPI_BYTE, l_state.world_comm);
    assert(!rc); 

#ifdef DEBUG
    for (i = 0; i < nprocs; i++) {
        if (me == i)
            continue;
        fprintf(stdout,"[%d] Remote QP %d, Remote LID %u, Rkey %u, Lkey %u\n"
                " LBuf %p, RBuf %p\n", 
                me, rbuf.qp_num[i], rbuf.lid[i], rbuf.rkey[i], lbuf.mr->lkey,
                lbuf.buf, rbuf.buf[i]);
        fflush(stdout);
    }
#endif

    return 0;
}

int create_qp(void)
{
    struct ibv_qp_attr qp_attr;
	int i;
    
    memset(&qp_attr, 0, sizeof qp_attr);

    struct ibv_qp_init_attr attr;

    memset(&attr, 0, sizeof attr);
    
    attr.send_cq = hca.cq;
    attr.recv_cq = hca.cq;
    attr.cap.max_send_wr  = opts.send_wr;
    attr.cap.max_recv_wr  = opts.recv_wr;
    attr.cap.max_send_sge = opts.num_sge;
    attr.cap.max_recv_sge = opts.num_sge;
    attr.cap.max_inline_data = 1;
    attr.qp_type = IBV_QPT_RC;

    // Create a connection to yourself 
    for(i = 0; i < nprocs; i++) {
        conn.qp[i] = ibv_create_qp(hca.pd, &attr);
        if(!conn.qp[i]) {
            fprintf(stderr,"Couldn't create QP\n");
            return 1;
        }

        conn.qp_num[i] = conn.qp[i]->qp_num;
        qp_attr.qp_state = IBV_QPS_INIT;
        qp_attr.pkey_index = 0;
        qp_attr.port_num   = 1;

        qp_attr.qp_access_flags = IBV_ACCESS_LOCAL_WRITE| 
            IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_REMOTE_READ |
            IBV_ACCESS_REMOTE_ATOMIC;

        if(ibv_modify_qp(conn.qp[i], &qp_attr,
                    IBV_QP_STATE              |
                    IBV_QP_PKEY_INDEX         |
                    IBV_QP_PORT               |
                    IBV_QP_ACCESS_FLAGS)) {
            fprintf(stderr,"Could not modify QP to INIT\n");
            return 1;
        }
#ifdef DEBUG
        fprintf(stdout,"[%d] Created QP %d, LID %d\n", me, 
                conn.qp_num[i], conn.lid[i]);
        fflush(stdout);
#endif
    }

    return 0;
}

int get_lid(void)
{
    struct ibv_port_attr port_attr[MAX_PORTS];
    int i, j;

    int active_port_found = 0;

    for (j = 1; j <= opts.num_ports; j++) {
        if (!ibv_query_port(hca.context, j, &port_attr[j - 1]) 
                && (port_attr[j - 1].state == IBV_PORT_ACTIVE)) {
            active_port_found = 1;
            for (i = 0; i < nprocs; i++)
                conn.lid[i] = port_attr[j - 1].lid;
            return 0;
        } 
    }

    assert(active_port_found);
    
    return 1;
}

int connect_qp(void)
{
    struct ibv_qp_attr attr;
	int i;

    memset(&attr, 0 , sizeof attr);

    for(i = 0; i < nprocs; i++){

	    attr.qp_state       = IBV_QPS_RTR;
    	attr.path_mtu       = IBV_MTU_2048;
	    attr.dest_qp_num    = rbuf.qp_num[i];
	    attr.rq_psn         = 0;
    	attr.max_dest_rd_atomic = 10;
	    attr.min_rnr_timer  = 20;
	    attr.ah_attr.is_global  = 0;
	    attr.ah_attr.dlid       = rbuf.lid[i];
	    attr.ah_attr.sl         = 0;
    	attr.ah_attr.src_path_bits = 0;
		attr.ah_attr.port_num   = 1;
		attr.ah_attr.static_rate = 0;
	    if (ibv_modify_qp(conn.qp[i], &attr,
    	            IBV_QP_STATE              |
        	        IBV_QP_AV                 |
            	    IBV_QP_PATH_MTU           |
                	IBV_QP_DEST_QPN           |
	                IBV_QP_RQ_PSN             |
    	            IBV_QP_MAX_DEST_RD_ATOMIC |
        	        IBV_QP_MIN_RNR_TIMER)) {
	        fprintf(stderr, "Failed to modify QP to RTR\n");
    	    return 1;
    	}

	    attr.qp_state       = IBV_QPS_RTS;
    	attr.timeout        = 14;
	    attr.retry_cnt      = 7;
    	attr.rnr_retry      = 7;
	    attr.sq_psn         = 0;
    	attr.max_rd_atomic  = 0;
	    if (ibv_modify_qp(conn.qp[i], &attr,
	                IBV_QP_STATE              |
    	            IBV_QP_TIMEOUT            |
        	        IBV_QP_RETRY_CNT          |
	                IBV_QP_RNR_RETRY          |
    	            IBV_QP_SQ_PSN             |
        	        IBV_QP_MAX_QP_RD_ATOMIC)) {
	        fprintf(stderr, "Failed to modify QP to RTS\n");
    	    return 1;
    	}
    }
    return 0;
}

struct ibv_send_wr  sr_desc;
struct ibv_sge      sr_sg_entry;

// All data transfers must go through this function
void prepare_and_post_send_desc(void *src, void *dst,
        int dest, int len, int lkey, int rkey, int type, int lock_or_unlock)
{
    sr_desc.send_flags = IBV_SEND_SIGNALED;
    sr_desc.next = NULL;
    sr_desc.opcode = type;
    sr_desc.wr_id = 0;
    sr_desc.num_sge = 1;

    if(IBV_WR_RDMA_WRITE == type) { 
        sr_desc.wr.rdma.remote_addr = (uintptr_t) (dst);
        sr_sg_entry.addr = (uintptr_t) (src);
        sr_sg_entry.length = len;
        sr_desc.wr.rdma.rkey = rkey;
    }

    if (IBV_WR_RDMA_READ == type) {
        sr_desc.wr.rdma.remote_addr = (uintptr_t) (src);
        sr_sg_entry.addr = (uintptr_t) (dst);
        sr_sg_entry.length = len;
        sr_desc.wr.rdma.rkey = rkey;
    }

    if (IBV_WR_ATOMIC_CMP_AND_SWP == type) {
        sr_desc.wr.atomic.remote_addr = (uintptr_t) (dst);
        sr_desc.wr.atomic.rkey = rkey;
        sr_sg_entry.addr = (uintptr_t) (src);
        sr_sg_entry.length = sizeof(long);

        if (lock_or_unlock == OPENIB_LOCK) {
            sr_desc.wr.atomic.compare_add = 0;
            sr_desc.wr.atomic.swap = l_state.rank + 1;
        }
        else if (lock_or_unlock == OPENIB_UNLOCK){
            sr_desc.wr.atomic.compare_add = l_state.rank + 1;
            sr_desc.wr.atomic.swap = 0;
        }
        else {
            assert(0);
        }
    }
    sr_sg_entry.lkey = lkey;

    sr_desc.sg_list = &(sr_sg_entry);
    struct ibv_send_wr *bad_wr;
    
    if(ibv_post_send(conn.qp[dest], &sr_desc, &bad_wr)) {
        fprintf(stderr,"[%d] Error posting send\n",
                me);
        fflush(stderr);
    }
   
    // Increment outstanding and check whether we need to make progress 
    increment_outstanding();
}


/* Poll CQ till n_sends are complete and n_recvs
 * are complete
 */
int poll_cq(int n_sends, int n_recvs)
{
    int ne;
    struct ibv_wc wc;
    int send_comp = 0;
    int recv_comp = 0;

    // return if poll cq is called with no outstanding messages
    if (!n_sends)
        return 0;

    // Make sure that no recv requests are called
    assert(!n_recvs);

    while (1) {
        do {
            ne = ibv_poll_cq(hca.cq, 1, &wc);
        } while (ne < 1);

        /* Okay, got an entry, check for errors */
        if(ne < 0) {
            fprintf(stderr,"Error Polling CQ\n");
            return 1;
        }

        if(wc.status != IBV_WC_SUCCESS) {
            fprintf(stderr, "[%d] Failed status %d\n", 
                    me, wc.status);
            return 1;
        }
        /* What type is it? */
        /* In this test suite, we are using either send descriptors
         * or RDMA, hence the following is OK */
        if(IBV_WC_RDMA_WRITE == wc.opcode || IBV_WC_RDMA_READ == wc.opcode
                || IBV_WC_COMP_SWAP == wc.opcode) {
            /* Send Completion */
            if(n_sends) {
                send_comp++;
            }
        } else {
            fprintf(stderr, "Unknown opcode recv'd\n");
            assert(0);
        }

        if((n_sends == send_comp) &&
                (n_recvs == recv_comp)) {
            break;
        }
    }

    return 0;
}

void openib_waitall()
{
    // Poll the CQ for all outstanding data transfer
    int rc = poll_cq(l_state.num_outstanding, 0);
    assert(!rc);
    
    l_state.num_outstanding = 0;
}

void increment_outstanding()
{
    l_state.num_outstanding++;

    if ((l_state.num_outstanding == 
            DEFAULT_NUM_OUST_SEND)) {
        openib_waitall();
    }
}

void * openib_register_memory(void *buf, int len)
{
    struct ibv_mr *mr;
    mr = ibv_reg_mr(hca.pd, buf, len,
            IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_WRITE
            | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_ATOMIC);

    if (mr) {
        // Insert in the registration cache
        reg_cache_insert(l_state.rank, buf, len, 
                mr->lkey, mr->rkey, mr);
    }

    return (void *)mr;
}

int openib_deregister_memory(void *buf)
{
    struct _reg_entry_t *reg = reg_cache_find(l_state.rank, buf, 0);
    assert(reg);

    if (ibv_dereg_mr(reg->mr)) {
        assert(0);
    }

    reg_cache_delete(l_state.rank, buf);
    return 0;
}

int openib_put_nbi(void *src, void *dst, int bytes, int proc)
{
    void *src_ptr;
    int local_reg_failure = 0;

    struct _reg_entry_t *local_reg, *remote_reg;

    // Search for local key
    local_reg = reg_cache_find(l_state.rank, src, bytes);

    // Search for remote key
    remote_reg = reg_cache_find(proc, dst, bytes);

    if (!local_reg && l_state.comex_openib_use_dreg) {
        openib_register_memory(src, bytes);
        local_reg = reg_cache_find(l_state.rank, src, bytes);
    }

    if (local_reg) {
        src_ptr = src;
    }
    else {
        local_reg_failure = 1;
        openib_waitall();
        src_ptr = l_state.put_buf;
        assert(bytes <= l_state.put_buf_len);
        memcpy(src_ptr, src, bytes);
        local_reg = reg_cache_find(l_state.rank, src_ptr, bytes);
        assert(local_reg);
    }
    // Ensure that the registration entries are valid
    assert(remote_reg);

    // Prepare a send decriptor and post
    prepare_and_post_send_desc(src_ptr, dst, proc, bytes,
            local_reg->lkey, remote_reg->rkey, IBV_WR_RDMA_WRITE, -1);

    if (local_reg_failure) {
        openib_waitall();
    }

    return 0;
}

int openib_get_nbi(void *src, void *dst, int bytes, int proc)
{
    void *dst_ptr;
    int local_reg_failure = 0;

    // Due to location consistency semantics, we need to call waitall here
    //openib_waitall();

    // The source is the buffer from which the data is read
    struct _reg_entry_t *local_reg, *remote_reg;

    // Search for local key
    local_reg = reg_cache_find(l_state.rank, dst, bytes);

    // Search for remote key
    remote_reg = reg_cache_find(proc, src, bytes);
    assert(remote_reg);

    if (!local_reg && l_state.comex_openib_use_dreg) {
        openib_register_memory(dst, bytes);
        local_reg = reg_cache_find(l_state.rank, dst, bytes);
    }

    // Ensure that the registration entries are valid
    if (local_reg) {
        dst_ptr = dst;
    }
    else {
        local_reg_failure = 1;
        assert(bytes <= l_state.get_buf_len);
        dst_ptr = l_state.get_buf;
        local_reg = reg_cache_find(l_state.rank, dst_ptr, bytes);
        assert(local_reg);
        openib_waitall();
    }

    assert(remote_reg);

    // Prepare a send decriptor and post
    prepare_and_post_send_desc(src, dst_ptr, proc, bytes,
            local_reg->lkey, remote_reg->rkey, IBV_WR_RDMA_READ, -1);

    if (local_reg_failure) {
        openib_waitall();
        memcpy(dst, dst_ptr, bytes);
    }

    return 0;
}

void openib_network_lock(int proc)
{
    void *src = l_state.local_lock_buf;
    void *dst = l_state.atomic_lock_buf[proc];

    assert(src && dst);
    int bytes = sizeof(long);
    // The source is the buffer from which the data is read
    struct _reg_entry_t *local_reg, *remote_reg;

    // Search for local key
    local_reg = reg_cache_find(l_state.rank, 
            src, sizeof(long));

    // Search for remote key
    remote_reg = reg_cache_find(proc, dst, 
            sizeof(long));

    // Ensure that the registration entries are valid
    assert(local_reg && remote_reg);

    // Prepare a send decriptor and post
    //

    do {
        prepare_and_post_send_desc(src, dst, proc, bytes, local_reg->lkey,
                remote_reg->rkey, IBV_WR_ATOMIC_CMP_AND_SWP, OPENIB_LOCK);
        openib_waitall();
    } while (*(long *)(src) != 0);
}

void openib_network_unlock(int proc)
{
    void *src = l_state.local_lock_buf;
    void *dst = l_state.atomic_lock_buf[proc];

    assert(src && dst);

    int bytes = sizeof(long);
    // The source is the buffer from which the data is read
    struct _reg_entry_t *local_reg, *remote_reg;

    // Search for local key
    local_reg = reg_cache_find(l_state.rank, 
            src, sizeof(long));

    // Search for remote key
    remote_reg = reg_cache_find(proc, dst, 
            sizeof(long));

    // Ensure that the registration entries are valid
    assert(local_reg);
    assert(remote_reg);

    // Prepare a send decriptor and post
    do {
        prepare_and_post_send_desc(src, dst, proc, bytes,
                local_reg->lkey, remote_reg->rkey, 
                IBV_WR_ATOMIC_CMP_AND_SWP, OPENIB_UNLOCK);
        openib_waitall();
    } while (*(long *)(src) != l_state.rank + 1);
}

void openib_create_locks()
{
    // Create the locks and initialize them
    l_state.local_lock_buf = comex_malloc_local(sizeof(long));
    assert(l_state.local_lock_buf);

    l_state.atomic_lock_buf = (void **)malloc(l_state.size * sizeof(void *));
    assert(l_state.atomic_lock_buf);

    comex_malloc((l_state.atomic_lock_buf), sizeof(long), COMEX_GROUP_WORLD);

    *(long *)(l_state.atomic_lock_buf[l_state.rank]) = 0;
    *(long *)(l_state.local_lock_buf) = 0;
    
    MPI_Barrier(l_state.world_comm);
}


void openib_alloc_buf()
{
    l_state.acc_buf_len = 1048576;
    l_state.acc_buf = (void *)malloc(sizeof(char) * l_state.acc_buf_len);

    l_state.put_buf_len = 1048576;
    l_state.put_buf = (void *)malloc(sizeof(char) * l_state.put_buf_len);
   
    l_state.get_buf_len = 1048576;
    l_state.get_buf = (void *)malloc(sizeof(char) * l_state.get_buf_len);
    
    assert(l_state.acc_buf && l_state.put_buf && l_state.get_buf);

    void *info;
    
    info = openib_register_memory(l_state.acc_buf, l_state.acc_buf_len);
    assert(info);
    
    info = openib_register_memory(l_state.put_buf, l_state.put_buf_len);
    assert(info);
    
    info = openib_register_memory(l_state.get_buf, l_state.get_buf_len);
    assert(info);
}

int openib_waitproc(int proc)
{
    openib_waitall();
    return 0;
}

int openib_finalize()
{
    release_resources();
    return 0;
}

int openib_initialize()
{

    // Use the previously cached info
    me = l_state.rank;
    nprocs = l_state.size;
    assert(l_state.world_comm);

    // initialize the envs
    openib_init_envs();

    //Initialize the registration cache

    reg_cache_init(nprocs, 0);

    init_params();

    if(open_hca()) {
        release_resources();
        exit(1);
    }

    if(create_cq()) {
        release_resources();
        exit(1);
    }

    if(get_lid()) {
        release_resources();
        exit(1);
    }

    if(create_qp()) {
        release_resources();
        exit(1);
    }

    if(exch_addr()) {
        release_resources();
        exit(1);
    }

    if(connect_qp()) {
        release_resources();
        exit(1);
    }

    // Create network locks
    openib_create_locks();

    // Allocate buffers for one sided operations
    openib_alloc_buf();

    MPI_Barrier(l_state.world_comm);

    return 0;
}

