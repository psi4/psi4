#ifndef OPENIB_H_
#define OPENIB_H_

#include <stdint.h> /* for uint16_t, uint32_t */

#define OPENIB_LOCK 0
#define OPENIB_UNLOCK 1

/* Default values */
#define MEMCPY_ITERS            (10000)
#define MAX_PORT_NUM            (2)
#define DEFAULT_PORT            (1)
#define MEMCPY_SIZE             (16 * 1024 * 1024)
#define MAX_PORTS               (1)
#define MAX_QP_PER_PORT         (1)
#define MAX_SUBCHANNELS         (MAX_PORTS * MAX_QP_PER_PORT)
#define DEFAULT_MSG_SIZE        (8 * 1048576)
#define DEFAULT_BW_ITERS        (200)
#define DEFAULT_ITERATIONS      (10000)
#define DEFAULT_NUM_CQE         (40000)
#define DEFAULT_NUM_OUST_RECV   (0)
#define DEFAULT_NUM_OUST_SEND   (128)
#define DEFAULT_NUM_SGE         (1)
#define DEFAULT_WINDOW_SIZE     (32)
#define ALIGNMENT               (64)
#define STRIPING_THRESHOLD (16 * 1024)
#define MHZ 238
#define SKIP 10 

#define rdtsc(x) asm volatile("rdtsc" : "=A" (x))

/* Read and Write Barriers for PPC, required */
#ifdef _PPC64_
#   define STBAR()    asm volatile ("sync": : :"memory") /* ": : :" for C++ */
#   define READBAR()  asm volatile ("sync": : :"memory")
#   define WRITEBAR() asm volatile ("eieio": : :"memory")
#else
#   define STBAR()
#   define READBAR()
#   define WRITEBAR()
#endif

struct User_Opts
{
    int     msg_size;
    int     iterations;
    int     rdma_flag;
    int     sendrecv_flag;
    char    *ib_devname;
    int     all_msgs;
    int     align;
    int     num_cqe;
    int     send_wr;
    int     recv_wr;
    int     num_sge;
    int     latency;
    int     bandwidth;
    int     bibandwidth;
    int     window;
    int     num_ports;
    int     num_qp_per_port;
    int     striping_threshold;
    int     subchannels;
    int     regcost;
    int     memperf;
    int     default_port;
    int     use_apm;
    int     apm_test;
    int     test_nft;
    int     use_srq;
};

struct HCA
{
    struct ibv_device   *ib_dev;
    struct ibv_context  *context;
    struct ibv_pd       *pd;
    struct ibv_cq       *cq;
    struct ibv_srq      *srq_hndl;
    struct ibv_comp_channel *scomp_ch;
    struct ibv_comp_channel *rcomp_ch;
    void    *cntx;
};

struct Local_Buf
{
    char            *buf_original;
    char            *buf;
    char            *tmp;
    struct ibv_mr   *mr;
};

struct Remote_Buf
{
    char            **buf;
    uint32_t        *qp_num;
    uint16_t        *lid;
    uint32_t        *rkey;
};

struct RC_Conn
{
    struct ibv_qp       **qp;
    uint16_t            *lid;
    uint32_t        *qp_num;
};

//static struct ibv_srq * create_srq();

static inline int min(int a, int b)
{
    return (a < b ? a : b);
}

static inline int max(int a, int b)
{
    return (a > b ? a : b);
}

extern void post_recv_desc(int dest);
extern void init_params(void);
extern double my_wtime();
extern void release_resources(void);
extern int open_hca(void);
extern int reg_buf(void);
extern int create_cq(void);
extern int exch_addr(void);
extern int create_qp(void);
extern int get_lid(void);
extern int connect_qp(void);
extern void prepare_and_post_send_desc(void *src, void *dst,
        int dest, int len, int lkey, int rkey, int type, int lock_or_unlock);
extern int poll_cq(int n_sends, int n_recvs);
extern void openib_waitall();
extern void increment_outstanding();
extern void * openib_register_memory(void *buf, int len);
extern int openib_deregister_memory(void *buf);
extern int openib_put_nbi(void *src, void *dst, int bytes, int proc);
extern int openib_get_nbi(void *src, void *dst, int bytes, int proc);
extern void openib_network_lock(int proc);
extern void openib_network_unlock(int proc);
extern void openib_create_locks();
extern void openib_alloc_buf();
extern int openib_waitproc(int proc);
extern int openib_finalize();
extern int openib_initialize();

#endif /* OPENIB_H_ */
