#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* C and/or system headers */
#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>

/* 3rd party headers */
#include <mpi.h>
#include <portals4.h>

/* our headers */
#include "acc.h"
#include "comex.h"
#include "comex_impl.h"
#include "groups.h"

#define CHECK_PTL_RETVAL(retval) check_ptl_retval((retval), __FILE__, __LINE__)
#define CHECK_MPI_RETVAL(retval) check_mpi_retval((retval), __FILE__, __LINE__)
#if   SIZEOF_VOIDP == SIZEOF_MPI_AINT
#   define MPI_VOIDP MPI_AINT
#elif SIZEOF_VOIDP == SIZEOF_INT
#   define MPI_VOIDP MPI_INT
#elif SIZEOF_VOIDP == SIZEOF_LONG
#   define MPI_VOIDP MPI_LONG
#elif SIZEOF_VOIDP == SIZEOF_LONG_LONG
#   define MPI_VOIDP MPI_LONG_LONG_INT
#else
#   error cannot determine MPI_Datatype for void*
#endif

/* keep track of all nonblocking operations */
typedef struct {
    int in_use; /**< 1 == user-visible handle */
    int send;
    int reply;
    int ack;
} nb_t;
static nb_t *nb_state = NULL;
static int nb_index = 0;
static int nb_count_event = 0;
static int nb_count_event_processed = 0;
static int nb_count_send = 0;
static int nb_count_send_processed = 0;
static int nb_count_reply = 0;
static int nb_count_reply_processed = 0;
static int nb_count_ack = 0;
static int nb_count_ack_processed = 0;

/* exported state */
local_state l_state;

/* static state */
static int  initialized=0;  /* for comex_initialized(), 0=false */
static char skip_lock=0;    /* don't acquire or release lock */
static int  nb_max_outstanding = 0;
static void* acc_buf = NULL;
static void* rmw_buf = NULL;
static void* swp_buf = NULL;

/* static function declarations */
static void acquire_remote_lock(int proc);
static void release_remote_lock(int proc);

static inline void check_ptl_retval(int retval, const char *file, int line);
static inline const char * str_ptl_retval(int retval);
static inline void check_mpi_retval(int retval, const char *file, int line);
static inline const char * str_mpi_retval(int retval);

static inline void nb_process_event();
static inline void nb_wait_for_event_space();
static inline int  nb_get_handle_index();
static inline nb_t* nb_wait_for_handle();
static inline void nb_wait_for_send(nb_t *nb);
static inline void nb_wait_for_reply(nb_t *nb);
static inline void nb_wait_for_ack(nb_t *nb);
static inline void nb_wait_for_all(nb_t *nb);
static inline void nb_wait_all();
static inline void nb_put(void *src, void *dst, int bytes, int proc, nb_t *nb);
static inline void nb_get(void *src, void *dst, int bytes, int proc, nb_t *nb);
static inline void nb_acc(int datatype, void *scale, void *src, void *dst, int bytes, int proc, nb_t *nb);
static inline void nb_puts(
        void *src, int *src_stride, void *dst, int *dst_stride,
        int *count, int stride_levels, int proc, nb_t *nb);
static inline void nb_gets(
        void *src, int *src_stride, void *dst, int *dst_stride,
        int *count, int stride_levels, int proc, nb_t *nb);
static inline void nb_accs(
        int datatype, void *scale,
        void *src, int *src_stride, void *dst, int *dst_stride,
        int *count, int stride_levels, int proc, nb_t *nb);
static inline void nb_putv(comex_giov_t *iov, int iov_len, int proc, nb_t *nb);
static inline void nb_getv(comex_giov_t *iov, int iov_len, int proc, nb_t *nb);
static inline void nb_accv(int datatype, void *scale,
        comex_giov_t *iov, int iov_len, int proc, nb_t *nb);

int comex_init()
{
    int init_flag = 0;
    int status = 0;
    ptl_ni_limits_t ptl_ni_limits_requested;
    ptl_md_t ptl_md;
    ptl_le_t ptl_le;
    int nb_max_outstanding_initial = 0;
    
    if (initialized) {
        return 0;
    }
    initialized = 1;

    /* Assert MPI has been initialized */
    status = MPI_Initialized(&init_flag);
    CHECK_MPI_RETVAL(status);
    assert(init_flag);
    
    /* Duplicate the World Communicator */
    status = MPI_Comm_dup(MPI_COMM_WORLD, &(l_state.world_comm));
    CHECK_MPI_RETVAL(status);
    assert(l_state.world_comm); 

    /* My Rank */
    status = MPI_Comm_rank(l_state.world_comm, &(l_state.rank));
    CHECK_MPI_RETVAL(status);

    /* World Size */
    status = MPI_Comm_size(l_state.world_comm, &(l_state.size));
    CHECK_MPI_RETVAL(status);
    
    /* groups */
    comex_group_init();

    /* nonblocking state */
    nb_state = NULL; /* initialized later */
    nb_index = 0;
    nb_count_event = 0;
    nb_count_event_processed = 0;
    nb_count_send = 0;
    nb_count_send_processed = 0;
    nb_count_reply = 0;
    nb_count_reply_processed = 0;
    nb_count_ack = 0;
    nb_count_ack_processed = 0;
    nb_max_outstanding = COMEX_MAX_NB_OUTSTANDING; /* initial hint */

    /* env vars */
    {
        char *value = NULL;
        if ((value = getenv("COMEX_MAX_NB_OUTSTANDING")) != NULL) {
            nb_max_outstanding = atoi(value);
        }
    }

    /* portals4 */

    /* init portals */
    status = PtlInit();
    CHECK_PTL_RETVAL(status);

    /* init portals network */
    ptl_ni_limits_requested.max_entries = INT_MAX;
    ptl_ni_limits_requested.max_unexpected_headers = INT_MAX;
    ptl_ni_limits_requested.max_mds = INT_MAX;
    ptl_ni_limits_requested.max_eqs = INT_MAX;
    ptl_ni_limits_requested.max_cts = INT_MAX;
    ptl_ni_limits_requested.max_pt_index = INT_MAX;
    ptl_ni_limits_requested.max_iovecs = INT_MAX;
    ptl_ni_limits_requested.max_list_size = INT_MAX;
    ptl_ni_limits_requested.max_triggered_ops = INT_MAX;
    ptl_ni_limits_requested.max_msg_size = LONG_MAX;
    ptl_ni_limits_requested.max_atomic_size = LONG_MAX;
    ptl_ni_limits_requested.max_fetch_atomic_size = LONG_MAX;
    ptl_ni_limits_requested.max_waw_ordered_size = LONG_MAX;
    ptl_ni_limits_requested.max_war_ordered_size = LONG_MAX;
    ptl_ni_limits_requested.max_volatile_size = LONG_MAX;
    ptl_ni_limits_requested.features = PTL_TARGET_BIND_INACCESSIBLE;

    status = PtlNIInit(PTL_IFACE_DEFAULT,
            PTL_NI_NO_MATCHING | PTL_NI_LOGICAL,
            PTL_PID_ANY,
            &ptl_ni_limits_requested,
            &l_state.ptl_ni_limits,
            &l_state.ptl_ni_handle);
    CHECK_PTL_RETVAL(status);

    if (0 == l_state.rank) {
        printf("max_entries = %d (%d requested)\n", 
                l_state.ptl_ni_limits.max_entries,
                ptl_ni_limits_requested.max_entries);
        printf("max_unexpected_headers = %d (%d requested)\n", 
                l_state.ptl_ni_limits.max_unexpected_headers,
                ptl_ni_limits_requested.max_unexpected_headers);
        printf("max_mds = %d (%d requested)\n", 
                l_state.ptl_ni_limits.max_mds,
                ptl_ni_limits_requested.max_mds);
        printf("max_eqs = %d (%d requested)\n", 
                l_state.ptl_ni_limits.max_eqs,
                ptl_ni_limits_requested.max_eqs);
        printf("max_cts = %d (%d requested)\n", 
                l_state.ptl_ni_limits.max_cts,
                ptl_ni_limits_requested.max_cts);
        printf("max_pt_index = %d (%d requested)\n", 
                l_state.ptl_ni_limits.max_pt_index,
                ptl_ni_limits_requested.max_pt_index);
        printf("max_iovecs = %d (%d requested)\n", 
                l_state.ptl_ni_limits.max_iovecs,
                ptl_ni_limits_requested.max_iovecs);
        printf("max_list_size = %d (%d requested)\n", 
                l_state.ptl_ni_limits.max_list_size,
                ptl_ni_limits_requested.max_list_size);
        printf("max_triggered_ops = %d (%d requested)\n", 
                l_state.ptl_ni_limits.max_triggered_ops,
                ptl_ni_limits_requested.max_triggered_ops);
        printf("max_msg_size = %ld (%ld requested)\n",
                l_state.ptl_ni_limits.max_msg_size,
                ptl_ni_limits_requested.max_msg_size);
        printf("max_atomic_size = %ld (%ld requested)\n",
                l_state.ptl_ni_limits.max_atomic_size,
                ptl_ni_limits_requested.max_atomic_size);
        printf("max_fetch_atomic_size = %ld (%ld requested)\n",
                l_state.ptl_ni_limits.max_fetch_atomic_size,
                ptl_ni_limits_requested.max_fetch_atomic_size);
        printf("max_waw_ordered_size = %ld (%ld requested)\n",
                l_state.ptl_ni_limits.max_waw_ordered_size,
                ptl_ni_limits_requested.max_waw_ordered_size);
        printf("max_war_ordered_size = %ld (%ld requested)\n",
                l_state.ptl_ni_limits.max_war_ordered_size,
                ptl_ni_limits_requested.max_war_ordered_size);
        printf("max_volatile_size = %ld (%ld requested)\n",
                l_state.ptl_ni_limits.max_volatile_size,
                ptl_ni_limits_requested.max_volatile_size);
        printf("features = %u (%u requested)\n",
                l_state.ptl_ni_limits.features,
                ptl_ni_limits_requested.features);
    }

    /* establish physical to logical rank mapping */
    //status = PtlGetId(l_state.ptl_ni_handle, &l_state.ptl_process_id);
    status = PtlGetPhysId(l_state.ptl_ni_handle, &l_state.ptl_process_id);
    CHECK_PTL_RETVAL(status);

    /* keep track of all process ids */
    l_state.ptl_process_ids = (ptl_process_t*)malloc(
            sizeof(ptl_process_t)*l_state.size);
    assert(NULL != l_state.ptl_process_ids);
    l_state.ptl_process_ids[l_state.rank] = l_state.ptl_process_id;
    status = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
            l_state.ptl_process_ids, sizeof(ptl_process_t),
            MPI_CHAR, l_state.world_comm);
    CHECK_MPI_RETVAL(status);

    /* set the mapping (same on all procs) */
    status = PtlSetMap(l_state.ptl_ni_handle, l_state.size,
            l_state.ptl_process_ids);
    CHECK_PTL_RETVAL(status);

    /* get this proc's uid */
    status = PtlGetUid(l_state.ptl_ni_handle, &l_state.ptl_uid);
    CHECK_PTL_RETVAL(status);

    /* create event queue for initiator events */
    nb_max_outstanding_initial = nb_max_outstanding;
    status = PtlEQAlloc(l_state.ptl_ni_handle, nb_max_outstanding*2,
            &l_state.ptl_eq_handle);
    while (PTL_NO_SPACE == status && nb_max_outstanding > 0) {
        /* half the max and try again */
        nb_max_outstanding /= 2;
        status = PtlEQAlloc(l_state.ptl_ni_handle, nb_max_outstanding*2,
                &l_state.ptl_eq_handle);
    }
    assert(nb_max_outstanding > 0);
    CHECK_PTL_RETVAL(status);
    if (0 == l_state.rank) {
        printf("event queue size = %d (%d requested)\n",
                nb_max_outstanding,
                nb_max_outstanding_initial);
    }

    /* create nonblocking state registers */
    nb_state = (nb_t*)malloc(sizeof(nb_t) * nb_max_outstanding);
    assert(NULL != nb_state);
    (void)memset(nb_state, 0, sizeof(nb_t) * nb_max_outstanding);

    /* create memory descriptor */
    ptl_md.start = 0;
    ptl_md.length = PTL_SIZE_MAX;
    ptl_md.options = PTL_MD_UNORDERED;
    ptl_md.eq_handle = l_state.ptl_eq_handle;
    ptl_md.ct_handle = PTL_CT_NONE;
    status = PtlMDBind(l_state.ptl_ni_handle, &ptl_md, &l_state.ptl_md_handle);
    CHECK_PTL_RETVAL(status);

    /* allocate a single portal table entry as the target for all ops */
    status = PtlPTAlloc(l_state.ptl_ni_handle, 0, PTL_EQ_NONE, PTL_PT_ANY,
            &l_state.ptl_pt_index);
    assert(l_state.ptl_pt_index == 0);
    CHECK_PTL_RETVAL(status);
    l_state.ptl_pt_indexes = (ptl_pt_index_t*)malloc(
            sizeof(ptl_pt_index_t)*l_state.size);
    l_state.ptl_pt_indexes[l_state.rank] = l_state.ptl_pt_index;
    status = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
            l_state.ptl_pt_indexes, sizeof(ptl_pt_index_t),
            MPI_CHAR, l_state.world_comm);
    CHECK_MPI_RETVAL(status);

    /* allocate a single list entry for the single portal table index */
    ptl_le.start = NULL;
    ptl_le.length = PTL_SIZE_MAX;
    ptl_le.ct_handle = PTL_CT_NONE;
    ptl_le.uid = PTL_UID_ANY;
    ptl_le.options = PTL_LE_OP_PUT | PTL_LE_OP_GET;
    status = PtlLEAppend(l_state.ptl_ni_handle, l_state.ptl_pt_index, &ptl_le,
            PTL_PRIORITY_LIST, NULL, &l_state.ptl_le_handle);
    CHECK_PTL_RETVAL(status);

    /* mutexes */
    l_state.mutexes = NULL;
    l_state.local_mutex = NULL;
    l_state.num_mutexes = NULL;

    if (l_state.rank == l_state.size-1) {
        printf("using %d process(es)\n", l_state.size);
    }

    /* Synch - Sanity Check */
    MPI_Barrier(l_state.world_comm);

    return COMEX_SUCCESS;
}


int comex_init_args(int *argc, char ***argv)
{
    int init_flag;
    
    MPI_Initialized(&init_flag);
    
    if(!init_flag) {
        MPI_Init(argc, argv);
    }
    
    return comex_init();
}


int comex_initialized()
{
    return initialized;
}


void comex_error(char *msg, int code)
{
    fprintf(stderr,"[%d] Received an Error in Communication: (%d) %s\n",
            l_state.rank, code, msg);
    
    MPI_Abort(l_state.world_comm, code);
}


int comex_put(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group)
{
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;

    nb = nb_wait_for_handle();

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_put(src, dst, bytes, world_proc, nb);
    nb_wait_for_send(nb);

    return COMEX_SUCCESS;
}


int comex_get(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group)
{
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;

    nb = nb_wait_for_handle();

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_get(src, dst, bytes, world_proc, nb);
    nb_wait_for_reply(nb);

    return COMEX_SUCCESS;
}


int comex_acc(
        int datatype, void *scale,
        void *src, void *dst, int bytes,
        int proc, comex_group_t group)
{
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;

    nb = nb_wait_for_handle();

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_acc(datatype, scale, src, dst, bytes, world_proc, nb);
    nb_wait_for_send(nb);

    return COMEX_SUCCESS;
}


int comex_puts(
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels,
        int proc, comex_group_t group)
{
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;

    nb = nb_wait_for_handle();

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_puts(src, src_stride, dst, dst_stride, count, stride_levels, world_proc, nb);
    nb_wait_for_send(nb);

    return COMEX_SUCCESS;
}


int comex_gets(
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels,
        int proc, comex_group_t group)
{
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;

    nb = nb_wait_for_handle();

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_gets(src, src_stride, dst, dst_stride, count, stride_levels, world_proc, nb);
    nb_wait_for_reply(nb);

    return COMEX_SUCCESS;
}


int comex_accs(
        int datatype, void *scale,
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels,
        int proc, comex_group_t group)
{
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;

    nb = nb_wait_for_handle();

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_accs(datatype, scale,
            src, src_stride, dst, dst_stride, count, stride_levels, world_proc, nb);
    nb_wait_for_send(nb);

    return COMEX_SUCCESS;
}


int comex_putv(
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group)
{
    int i = 0;
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;

    nb = nb_wait_for_handle();

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_putv(iov, iov_len, world_proc, nb);
    nb_wait_for_send(nb);

    return COMEX_SUCCESS;
}


int comex_getv(
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group)
{
    int i = 0;
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;

    nb = nb_wait_for_handle();

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_getv(iov, iov_len, world_proc, nb);
    nb_wait_for_reply(nb);

    return COMEX_SUCCESS;
}


int comex_accv(
        int datatype, void *scale,
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group)
{
    int i = 0;
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;

    nb = nb_wait_for_handle();

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_accv(datatype, scale, iov, iov_len, world_proc, nb);
    nb_wait_for_send(nb);

    return COMEX_SUCCESS;
}


int comex_fence_all(comex_group_t group)
{
    return comex_wait_all(group);
}


int comex_fence_proc(int proc, comex_group_t group)
{
    return comex_wait_all(group);
}


/* comex_barrier is comex_fence_all + MPI_Barrier */
int comex_barrier(comex_group_t group)
{
    MPI_Comm comm;

    comex_fence_all(group);
    assert(COMEX_SUCCESS == comex_group_comm(group, &comm));
    MPI_Barrier(comm);

    return COMEX_SUCCESS;
}


void *comex_malloc_local(size_t size)
{
    int retval = 0;
    void *memptr = NULL;

    retval = posix_memalign(&memptr, sizeof(void*), size);
    if (0 != retval) {
        errno = retval;
        perror("comex_malloc_local: posix_memalign");
        MPI_Abort(l_state.world_comm, retval);
    }

    return memptr;
}


int comex_free_local(void *ptr)
{
    free(ptr);

    return COMEX_SUCCESS;
}


int comex_finalize()
{
    /* it's okay to call multiple times -- extra calls are no-ops */
    if (!initialized) {
        return;
    }

    initialized = 0;

    /* Make sure that all outstanding operations are done */
    comex_wait_all(COMEX_GROUP_WORLD);
    
    /* groups */
    comex_group_finalize();

    /* portals4 */
    PtlFini();

    MPI_Barrier(l_state.world_comm);

    // destroy the communicators
    MPI_Comm_free(&l_state.world_comm);

    return COMEX_SUCCESS;
}


int comex_wait_proc(int proc, comex_group_t group)
{
    return comex_wait_all(group);
}


int comex_wait(comex_request_t* hdl)
{
    int index = 0;
    nb_t *nb = NULL;

    assert(NULL != hdl);

    index = *(int*)hdl;
    assert(index >= 0);
    assert(index < nb_max_outstanding);
    nb = &nb_state[index];

    if (0 == nb->in_use) {
        fprintf(stderr, "{%d} comex_test Error: invalid handle\n",
                l_state.rank);
    }

    nb_wait_for_all(nb);

    nb->in_use = 0;

    return COMEX_SUCCESS;
}


int comex_test(comex_request_t* hdl, int *status)
{
    int index = 0;
    nb_t *nb = NULL;

    assert(NULL != hdl);

    index = *(int*)hdl;
    assert(index >= 0);
    assert(index < nb_max_outstanding);
    nb = &nb_state[index];

    if (0 == nb->in_use) {
        fprintf(stderr, "{%d} comex_test Error: invalid handle\n",
                l_state.rank);
    }

    if (nb->send == 0 && nb->reply == 0 && nb->ack == 0) {
        *status = 0;
        nb->in_use = 0;
    }
    else {
        *status = 1;
    }

    return COMEX_SUCCESS;
}


int comex_wait_all(comex_group_t group)
{
    nb_wait_all();

    return COMEX_SUCCESS;
}


int comex_nbput(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;
    comex_request_t _hdl = 0;

    nb = nb_wait_for_handle();
    _hdl = nb_get_handle_index();
    assert(&nb_state[_hdl] == nb);
    if (hdl) {
        *hdl = _hdl;
        nb->in_use = 1;
    }

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_put(src, dst, bytes, world_proc, nb);

    return COMEX_SUCCESS;
}


int comex_nbget(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;
    comex_request_t _hdl = 0;

    nb = nb_wait_for_handle();
    _hdl = nb_get_handle_index();
    assert(&nb_state[_hdl] == nb);
    if (hdl) {
        *hdl = _hdl;
        nb->in_use = 1;
    }

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_get(src, dst, bytes, world_proc, nb);

    return COMEX_SUCCESS;
}


int comex_nbacc(
        int datatype, void *scale,
        void *src, void *dst, int bytes,
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;
    comex_request_t _hdl = 0;

    nb = nb_wait_for_handle();
    _hdl = nb_get_handle_index();
    assert(&nb_state[_hdl] == nb);
    if (hdl) {
        *hdl = _hdl;
        nb->in_use = 1;
    }

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_acc(datatype, scale, src, dst, bytes, world_proc, nb);

    return COMEX_SUCCESS;
}


int comex_nbputs(
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels, 
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;
    comex_request_t _hdl = 0;

    nb = nb_wait_for_handle();
    _hdl = nb_get_handle_index();
    assert(&nb_state[_hdl] == nb);
    if (hdl) {
        *hdl = _hdl;
        nb->in_use = 1;
    }

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_puts(src, src_stride, dst, dst_stride, count, stride_levels, world_proc, nb);

    return COMEX_SUCCESS;
}


int comex_nbgets(
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels, 
        int proc, comex_group_t group,
        comex_request_t *hdl) 
{
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;
    comex_request_t _hdl = 0;

    nb = nb_wait_for_handle();
    _hdl = nb_get_handle_index();
    assert(&nb_state[_hdl] == nb);
    if (hdl) {
        *hdl = _hdl;
        nb->in_use = 1;
    }

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_gets(src, src_stride, dst, dst_stride, count, stride_levels, world_proc, nb);

    return COMEX_SUCCESS;
}


int comex_nbaccs(
        int datatype, void *scale,
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels,
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;
    comex_request_t _hdl = 0;

    nb = nb_wait_for_handle();
    _hdl = nb_get_handle_index();
    assert(&nb_state[_hdl] == nb);
    if (hdl) {
        *hdl = _hdl;
        nb->in_use = 1;
    }

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_accs(datatype, scale,
            src, src_stride, dst, dst_stride, count, stride_levels, world_proc, nb);

    return COMEX_SUCCESS;
}


int comex_nbputv(
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group,
        comex_request_t* hdl)
{
    int i = 0;
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;
    comex_request_t _hdl = 0;

    nb = nb_wait_for_handle();
    _hdl = nb_get_handle_index();
    assert(&nb_state[_hdl] == nb);
    if (hdl) {
        *hdl = _hdl;
        nb->in_use = 1;
    }

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_putv(iov, iov_len, world_proc, nb);

    return COMEX_SUCCESS;
}


int comex_nbgetv(
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group,
        comex_request_t* hdl)
{
    int i = 0;
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;
    comex_request_t _hdl = 0;

    nb = nb_wait_for_handle();
    _hdl = nb_get_handle_index();
    assert(&nb_state[_hdl] == nb);
    if (hdl) {
        *hdl = _hdl;
        nb->in_use = 1;
    }

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_getv(iov, iov_len, world_proc, nb);

    return COMEX_SUCCESS;
}


int comex_nbaccv(
        int datatype, void *scale,
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group,
        comex_request_t* hdl)
{
    int i = 0;
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;
    comex_request_t _hdl = 0;

    nb = nb_wait_for_handle();
    _hdl = nb_get_handle_index();
    assert(&nb_state[_hdl] == nb);
    if (hdl) {
        *hdl = _hdl;
        nb->in_use = 1;
    }

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    nb_accv(datatype, scale, iov, iov_len, world_proc, nb);

    return COMEX_SUCCESS;
}


int comex_rmw(
        int op, void *ploc, void *prem, int extra,
        int proc, comex_group_t group)
{
    ptl_process_t peer;
    nb_t *nb = NULL;
    int world_proc = -1;
    int status = 0;

    if (DEBUG) {
        printf("[%d] comex_rmw(op=%d, ploc(%p)=%d, prem(%p), extra=%d, proc=%d, group=%d)\n",
                l_state.rank, op, ploc, *((int*)ploc), prem, extra, proc, group);
    }

    nb = nb_wait_for_handle();
    nb_wait_for_event_space();
    nb_count_event += 2;
    nb_count_send += 1;
    nb_count_reply += 1;
    nb->send += 1;
    nb->reply += 1;

    status = comex_group_translate_world(group, proc, &world_proc);
    assert(COMEX_SUCCESS == status);

    /* allocate rmw buffer on first use */
    if (NULL == rmw_buf) {
        rmw_buf = comex_malloc_local(sizeof(long));
        assert(sizeof(long) <= l_state.ptl_ni_limits.max_fetch_atomic_size);
    }
    
    peer.rank = world_proc;

    if (op == COMEX_FETCH_AND_ADD) {
        *((int*)rmw_buf) = extra;
        if (DEBUG) {
            printf("[%d] comex_rmw ploc=%d rmw_buf=%d\n",
                    l_state.rank, *((int*)ploc), *((int*)rmw_buf));
        }
        status = PtlFetchAtomic(
                l_state.ptl_md_handle,
                (ptl_size_t)ploc,
                l_state.ptl_md_handle,
                (ptl_size_t)rmw_buf,
                sizeof(int),
                peer,
                l_state.ptl_pt_indexes[world_proc],
                (ptl_match_bits_t)0,
                (ptl_size_t)prem,
                nb,
                (ptl_hdr_data_t)0,
                PTL_SUM,
                PTL_INT32_T);
    }
    else if (op == COMEX_FETCH_AND_ADD_LONG) {
        *((long*)rmw_buf) = (long)extra;
        if (DEBUG) {
            printf("[%d] comex_rmw rmw_buf=%ld\n", l_state.rank, *((long*)rmw_buf));
        }
        status = PtlFetchAtomic(
                l_state.ptl_md_handle,
                (ptl_size_t)ploc,
                l_state.ptl_md_handle,
                (ptl_size_t)rmw_buf,
                sizeof(long),
                peer,
                l_state.ptl_pt_indexes[world_proc],
                (ptl_match_bits_t)0,
                (ptl_size_t)prem,
                nb,
                (ptl_hdr_data_t)0,
                PTL_SUM,
                PTL_INT64_T);
    }
    else if (op == COMEX_SWAP) {
        *((int*)rmw_buf) = *((int*)ploc);
        status = PtlSwap(
                l_state.ptl_md_handle,
                (ptl_size_t)ploc,
                l_state.ptl_md_handle,
                (ptl_size_t)rmw_buf,
                sizeof(int),
                peer,
                l_state.ptl_pt_indexes[world_proc],
                (ptl_match_bits_t)0,
                (ptl_size_t)prem,
                nb,
                (ptl_hdr_data_t)0,
                NULL,
                PTL_SWAP,
                PTL_INT32_T);
    }
    else if (op == COMEX_SWAP_LONG) {
        *((long*)rmw_buf) = *((long*)ploc);
        status = PtlSwap(
                l_state.ptl_md_handle,
                (ptl_size_t)ploc,
                l_state.ptl_md_handle,
                (ptl_size_t)rmw_buf,
                sizeof(long),
                peer,
                l_state.ptl_pt_indexes[world_proc],
                (ptl_match_bits_t)0,
                (ptl_size_t)prem,
                nb,
                (ptl_hdr_data_t)0,
                NULL,
                PTL_SWAP,
                PTL_INT64_T);
    }
    else  {
        assert(0);
    }

    nb_wait_for_all(nb);

    return COMEX_SUCCESS;
}


/* Mutex Operations */
int comex_create_mutexes(int num)
{
    int i=0;

    if (DEBUG) {
        printf("[%d] comex_create_mutexes(num=%d)\n",
                l_state.rank, num);
    }

    assert(NULL == l_state.mutexes);
    assert(NULL == l_state.local_mutex);
    assert(NULL == l_state.num_mutexes);

    /* every process knows how many mutexes created on every process */
    l_state.num_mutexes = (unsigned int*)malloc(l_state.size * sizeof(unsigned int));
    assert(l_state.num_mutexes);
    /* gather the counts */
    MPI_Allgather(&num, 1, MPI_INT,
            l_state.num_mutexes, 1, MPI_UNSIGNED, l_state.world_comm);

    /* create the 1 element buffer to hold a remote mutex */
    l_state.local_mutex = comex_malloc_local(sizeof(long));
    assert(l_state.local_mutex);
    /* init the local mutex holder to rank+1, indicating no mutex is held */
    *l_state.local_mutex = l_state.rank+1;
    MPI_Barrier(l_state.world_comm);

    /* create all of the mutexes */
    l_state.mutexes = (long**)malloc(l_state.size * sizeof(long*));
    assert(l_state.mutexes);
    comex_malloc((void**)l_state.mutexes, num*sizeof(long), COMEX_GROUP_WORLD);
    /* init all of my mutexes to 0 */
    for (i=0; i<num; ++i) {
        l_state.mutexes[l_state.rank][i] = 0;
    }

    MPI_Barrier(l_state.world_comm);

    return COMEX_SUCCESS;
}


int comex_destroy_mutexes()
{
    if (DEBUG) {
        printf("[%d] comex_destroy_mutexes()\n", l_state.rank);
    }

    MPI_Barrier(l_state.world_comm);

    /* you cannot free mutexes if one is in use */
    assert(*l_state.local_mutex == (long)(l_state.rank+1));
#ifndef NDEBUG
    {
        unsigned int i;
        for (i=0; i<l_state.num_mutexes[l_state.rank]; ++i) {
            long *mutexes = l_state.mutexes[l_state.rank];
            assert(mutexes[i] == 0);
        }
    }
#endif

    /* destroy mutex counts */
    free(l_state.num_mutexes);
    l_state.num_mutexes = NULL;

    /* destroy the 1 element buffer holding a remote mutex */
    comex_free_local(l_state.local_mutex);
    l_state.local_mutex = NULL;

    /* destroy the mutexes */
    comex_free(l_state.mutexes[l_state.rank], COMEX_GROUP_WORLD);
    free(l_state.mutexes);
    l_state.mutexes = NULL;

    MPI_Barrier(l_state.world_comm);

    return COMEX_SUCCESS;
}


int comex_lock(int mutex, int proc)
{
    ptl_process_t peer;
    nb_t *nb = NULL;
    int status = 0;

    if (DEBUG) {
        printf("[%d] comex_lock(mutex=%d, proc=%d)\n",
                l_state.rank, mutex, proc);
    }

    /* preconditions */
    assert(0 <= proc && proc < l_state.size);
    assert(0 <= mutex && ((unsigned int)mutex) < l_state.num_mutexes[proc]);
    /* you cannot lock mutexes if one is in use */
    assert(*l_state.local_mutex == (long)(l_state.rank+1));

    /* allocate rmw buffer on first use */
    if (NULL == rmw_buf) {
        rmw_buf = comex_malloc_local(sizeof(long));
        assert(sizeof(long) <= l_state.ptl_ni_limits.max_fetch_atomic_size);
    }
    if (NULL == swp_buf) {
        swp_buf = comex_malloc_local(sizeof(long));
        assert(sizeof(long) <= l_state.ptl_ni_limits.max_fetch_atomic_size);
    }

    peer.rank = proc;
    *((long*)rmw_buf) = l_state.rank + 1;
    *((long*)swp_buf) = 0; /* compare to 0; 0 means unlocked */

    do {
        nb = nb_wait_for_handle();
        nb_wait_for_event_space();
        nb_count_event += 2;
        nb_count_send += 1;
        nb_count_reply += 1;
        nb->send += 1;
        nb->reply += 1;

        status = PtlSwap(
                l_state.ptl_md_handle,
                (ptl_size_t)l_state.local_mutex,
                l_state.ptl_md_handle,
                (ptl_size_t)rmw_buf,
                sizeof(long),
                peer,
                l_state.ptl_pt_indexes[proc],
                (ptl_match_bits_t)0,
                (ptl_size_t)&l_state.mutexes[proc][mutex],
                nb,
                (ptl_hdr_data_t)0,
                swp_buf,
                PTL_CSWAP,
                PTL_INT64_T);

        nb_wait_for_all(nb);
    }
    while(*(l_state.local_mutex) != 0);

    return COMEX_SUCCESS;
}


int comex_unlock(int mutex, int proc)
{
    ptl_process_t peer;
    nb_t *nb = NULL;
    int status = 0;

    if (DEBUG) {
        printf("[%d] comex_unlock(mutex=%d, proc=%d)\n",
                l_state.rank, mutex, proc);
    }

    /* preconditions */
    assert(0 <= proc && proc < l_state.size);
    assert(0 <= mutex && ((unsigned int)mutex) < l_state.num_mutexes[proc]);
    /* you cannot unlock a mutex you haven't locked */
    assert(*l_state.local_mutex == 0);

    assert(NULL != rmw_buf);
    assert(NULL != swp_buf);

    peer.rank = proc;
    *((long*)rmw_buf) = 0; /* 0 == unlocked */
    *((long*)swp_buf) = l_state.rank + 1; /* value should be the one we put */

    do {
        nb = nb_wait_for_handle();
        nb_wait_for_event_space();
        nb_count_event += 2;
        nb_count_send += 1;
        nb_count_reply += 1;
        nb->send += 1;
        nb->reply += 1;

        status = PtlSwap(
                l_state.ptl_md_handle,
                (ptl_size_t)l_state.local_mutex,
                l_state.ptl_md_handle,
                (ptl_size_t)rmw_buf,
                sizeof(long),
                peer,
                l_state.ptl_pt_indexes[proc],
                (ptl_match_bits_t)0,
                (ptl_size_t)&l_state.mutexes[proc][mutex],
                nb,
                (ptl_hdr_data_t)0,
                swp_buf,
                PTL_CSWAP,
                PTL_INT64_T);

        nb_wait_for_all(nb);
    }
    while(*(l_state.local_mutex) == 0);
    assert(*(l_state.local_mutex) == (long)(l_state.rank + 1));

    return COMEX_SUCCESS;
}


int comex_malloc(void *ptrs[], size_t size, comex_group_t group)
{
    comex_igroup_t *igroup = NULL;
    MPI_Comm comm = MPI_COMM_NULL;
    int comm_rank = -1;
    int comm_size = -1;
    int status = 0;

    /* preconditions */
    assert(NULL != ptrs);
   
    igroup = comex_get_igroup_from_group(group);
    comm = igroup->comm;
    assert(MPI_COMM_NULL != comm);
    status = MPI_Comm_rank(comm, &comm_rank);
    CHECK_MPI_RETVAL(status);
    status = MPI_Comm_size(comm, &comm_size);
    CHECK_MPI_RETVAL(status);

    /* allocate and register segment */
    ptrs[comm_rank] = comex_malloc_local(sizeof(char)*size);
  
    /* exchange buffer address */
    status = MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, ptrs,
            1, MPI_VOIDP, comm);
    CHECK_MPI_RETVAL(status);

    status = MPI_Barrier(comm);
    CHECK_MPI_RETVAL(status);

    return COMEX_SUCCESS;
}


int comex_free(void *ptr, comex_group_t group)
{
    comex_igroup_t *igroup = NULL;
    MPI_Comm comm = MPI_COMM_NULL;
    int comm_rank = -1;
    int comm_size = -1;
    long **allgather_ptrs = NULL;
    int status = 0;

    /* preconditions */
    assert(NULL != ptr);

    igroup = comex_get_igroup_from_group(group);
    comm = igroup->comm;
    assert(MPI_COMM_NULL != comm);
    status = MPI_Comm_rank(comm, &comm_rank);
    CHECK_MPI_RETVAL(status);
    status = MPI_Comm_size(comm, &comm_size);
    CHECK_MPI_RETVAL(status);

    /* allocate receive buffer for exchange of pointers */
    allgather_ptrs = (long **)malloc(sizeof(void *) * comm_size);
    assert(allgather_ptrs);

    /* exchange of pointers */
    status = MPI_Allgather(&ptr, sizeof(void *), MPI_BYTE,
            allgather_ptrs, sizeof(void *), MPI_BYTE, comm);
    CHECK_MPI_RETVAL(status);

    /* TODO do something useful with pointers */

    /* remove my ptr from reg cache and free ptr */
    comex_free_local(ptr);
    free(allgather_ptrs);

    /* Is this needed? */
    status = MPI_Barrier(comm);
    CHECK_MPI_RETVAL(status);

    return COMEX_SUCCESS;
}


static void acquire_remote_lock(int proc)
{
    assert(0);
}


static void release_remote_lock(int proc)
{
    assert(0);
}


static inline void check_ptl_retval(int retval, const char *file, int line)
{
    if (PTL_OK != retval) {
        const char *msg = str_ptl_retval(retval);
        fprintf(stderr, "{%d} PTL Error: %s: line %d: %s\n",
                l_state.rank, file, line, msg);
        MPI_Abort(l_state.world_comm, retval);
    }
}


static inline const char * str_ptl_retval(int retval)
{
    const char *msg = NULL;

    switch(retval) {
        case PTL_OK              : msg = "PTL_OK"; break;
        case PTL_ARG_INVALID     : msg = "PTL_ARG_INVALID"; break;
        case PTL_CT_NONE_REACHED : msg = "PTL_CT_NONE_REACHED"; break;
        case PTL_EQ_DROPPED      : msg = "PTL_EQ_DROPPED"; break;
        case PTL_EQ_EMPTY        : msg = "PTL_EQ_EMPTY"; break;
        case PTL_FAIL            : msg = "PTL_FAIL"; break;
        case PTL_IGNORED         : msg = "PTL_IGNORED"; break;
        case PTL_IN_USE          : msg = "PTL_IN_USE"; break;
        case PTL_INTERRUPTED     : msg = "PTL_INTERRUPTED"; break;
        case PTL_LIST_TOO_LONG   : msg = "PTL_LIST_TOO_LONG"; break;
        case PTL_NO_INIT         : msg = "PTL_NO_INIT"; break;
        case PTL_NO_SPACE        : msg = "PTL_NO_SPACE"; break;
        case PTL_PID_IN_USE      : msg = "PTL_PID_IN_USE"; break;
        case PTL_PT_FULL         : msg = "PTL_PT_FULL"; break;
        case PTL_PT_EQ_NEEDED    : msg = "PTL_PT_EQ_NEEDED"; break;
        case PTL_PT_IN_USE       : msg = "PTL_PT_IN_USE"; break;
        default                  : msg = "DEFAULT"; break;
    }

    return msg;
}


static inline void check_mpi_retval(int retval, const char *file, int line)
{
    if (MPI_SUCCESS != retval) {
        const char *msg = str_mpi_retval(retval);
        fprintf(stderr, "{%d} MPI Error: %s: line %d: %s\n",
                l_state.rank, file, line, msg);
        MPI_Abort(l_state.world_comm, retval);
    }
}


static inline const char *str_mpi_retval(int retval)
{
    const char *msg = NULL;

    switch(retval) {
        case MPI_SUCCESS       : msg = "MPI_SUCCESS"; break;
        case MPI_ERR_BUFFER    : msg = "MPI_ERR_BUFFER"; break;
        case MPI_ERR_COUNT     : msg = "MPI_ERR_COUNT"; break;
        case MPI_ERR_TYPE      : msg = "MPI_ERR_TYPE"; break;
        case MPI_ERR_TAG       : msg = "MPI_ERR_TAG"; break;
        case MPI_ERR_COMM      : msg = "MPI_ERR_COMM"; break;
        case MPI_ERR_RANK      : msg = "MPI_ERR_RANK"; break;
        case MPI_ERR_ROOT      : msg = "MPI_ERR_ROOT"; break;
        case MPI_ERR_GROUP     : msg = "MPI_ERR_GROUP"; break;
        case MPI_ERR_OP        : msg = "MPI_ERR_OP"; break;
        case MPI_ERR_TOPOLOGY  : msg = "MPI_ERR_TOPOLOGY"; break;
        case MPI_ERR_DIMS      : msg = "MPI_ERR_DIMS"; break;
        case MPI_ERR_ARG       : msg = "MPI_ERR_ARG"; break;
        case MPI_ERR_UNKNOWN   : msg = "MPI_ERR_UNKNOWN"; break;
        case MPI_ERR_TRUNCATE  : msg = "MPI_ERR_TRUNCATE"; break;
        case MPI_ERR_OTHER     : msg = "MPI_ERR_OTHER"; break;
        case MPI_ERR_INTERN    : msg = "MPI_ERR_INTERN"; break;
        case MPI_ERR_IN_STATUS : msg = "MPI_ERR_IN_STATUS"; break;
        case MPI_ERR_PENDING   : msg = "MPI_ERR_PENDING"; break;
        case MPI_ERR_REQUEST   : msg = "MPI_ERR_REQUEST"; break;
        case MPI_ERR_LASTCODE  : msg = "MPI_ERR_LASTCODE"; break;
        default                : msg = "DEFAULT"; break;
    }

    return msg;
}


static inline void nb_process_event()
{
    ptl_event_t ptl_event;
    nb_t *nb = NULL;
    int status = 0;

    if (DEBUG) {
        printf("[%d] nb_process_event\n", l_state.rank);
    }

    assert(nb_count_event-nb_count_event_processed > 0);
#if 0
    status = PtlEQWait(l_state.ptl_eq_handle, &ptl_event);
    CHECK_PTL_RETVAL(status);
#else
    while (1) {
        unsigned int which = 1;
        status = PtlEQPoll(&l_state.ptl_eq_handle, 1, 500,
                &ptl_event, &which);
        if (PTL_OK == status) {
            assert(0 == which);
            break;
        }
        else if (PTL_EQ_EMPTY == status) {
            /* no event found, so try again */
            printf("[%d] nb_process_event timeout waiting for event; retry\n"
                    "send %d/%d reply %d/%d ack %d/%d all %d/%d\n",
                    l_state.rank,
                    nb_count_send_processed,  nb_count_send,  
                    nb_count_reply_processed, nb_count_reply, 
                    nb_count_ack_processed,   nb_count_ack,   
                    nb_count_event_processed, nb_count_event);
        }
        else {
            CHECK_PTL_RETVAL(status);
        }
    }
#endif

    assert(NULL != ptl_event.user_ptr);
    nb = (nb_t*)(ptl_event.user_ptr);
    assert(NULL != nb);

    /* check for error condition */
    if (PTL_NI_OK != ptl_event.ni_fail_type) {
        const char *msg = NULL;
        switch (ptl_event.ni_fail_type) {
            case PTL_NI_UNDELIVERABLE  : msg = "PTL_NI_UNDELIVERABLE"; break;
            case PTL_NI_PT_DISABLED    : msg = "PTL_NI_PT_DISABLED"; break;
            case PTL_NI_DROPPED        : msg = "PTL_NI_DROPPED"; break;
            case PTL_NI_PERM_VIOLATION : msg = "PTL_NI_PERM_VIOLATION"; break;
            case PTL_NI_OP_VIOLATION   : msg = "PTL_NI_OP_VIOLATION"; break;
            default                    : msg = "DEFAULT"; break;
        }
        fprintf(stderr, "{%d} PTL Error: %s: line %d: %s\n",
                l_state.rank, __FILE__, __LINE__, msg);
        MPI_Abort(l_state.world_comm, ptl_event.ni_fail_type);
    }

    /* decrement associated nonblocking handle */
    if (PTL_EVENT_REPLY == ptl_event.type) {
        if (DEBUG) {
            printf("[%d] nb_process_event PTL_EVENT_REPLY\n", l_state.rank);
        }
        assert(nb->reply > 0);
        nb_count_reply_processed += 1;
        nb->reply -= 1;
    }
    else if (PTL_EVENT_SEND == ptl_event.type) {
        if (DEBUG) {
            printf("[%d] nb_process_event PTL_EVENT_SEND\n", l_state.rank);
        }
        assert(nb->send > 0);
        nb_count_send_processed += 1;
        nb->send -= 1;
    }
    else if (PTL_EVENT_ACK == ptl_event.type) {
        if (DEBUG) {
            printf("[%d] nb_process_event PTL_EVENT_ACK\n", l_state.rank);
        }
        assert(nb->ack > 0);
        nb_count_ack_processed += 1;
        nb->ack -= 1;
    }
    else {
        fprintf(stderr, "{%d} nb_process_event Error: unrecognized event\n",
                l_state.rank);
        MPI_Abort(l_state.world_comm, ptl_event.type);
    }

    nb_count_event_processed += 1;
}


static inline void nb_wait_for_event_space()
{
    /* make sure we didn't overrun the event queue */
    assert((nb_count_event-nb_count_event_processed) <= nb_max_outstanding);

    /* we need at most space for two events */
    while ((nb_count_event-nb_count_event_processed) >= (nb_max_outstanding-2)) {
        nb_process_event();
    }
}


static inline int nb_get_handle_index()
{
    int value = 0;

    if (0 == nb_index) {
        value = nb_max_outstanding-1;
    }
    else {
        value = nb_index-1;
    }

    return value;
}


/** Poll event queue until resources become available. */
static inline nb_t* nb_wait_for_handle()
{
    nb_t *nb = NULL;
    int in_use_count = 0;

    /* find first handle that isn't associated with a user-level handle */
    /* make sure the handle we find has processed all events */
    /* the user can accidentally exhaust the available handles */
    do {
        ++in_use_count;
        if (in_use_count > nb_max_outstanding) {
            fprintf(stderr,
                    "{%d} nb_wait_for_handle Error: all user-level "
                    "nonblocking handles have been exhausted\n",
                    l_state.rank);
            MPI_Abort(l_state.world_comm, -1);
        }
        nb = &nb_state[nb_index++];
        nb_index %= nb_max_outstanding; /* wrap around if needed */
        nb_wait_for_all(nb);
    } while (nb->in_use);

    nb_wait_for_event_space();

    return nb;
}


static inline void nb_wait_for_send(nb_t *nb)
{
    assert(NULL != nb);

    while (nb->send > 0) {
        nb_process_event();
    }
}


static inline void nb_wait_for_reply(nb_t *nb)
{
    assert(NULL != nb);

    while (nb->reply > 0) {
        nb_process_event();
    }
}


static inline void nb_wait_for_ack(nb_t *nb)
{
    assert(NULL != nb);

    while (nb->ack > 0) {
        nb_process_event();
    }
}


static inline void nb_wait_for_all(nb_t *nb)
{
    if (DEBUG) {
        printf("[%d] nb_wait_for_all\n", l_state.rank);
    }
    assert(NULL != nb);

    while (nb->send > 0 || nb->reply > 0 || nb->ack > 0) {
        nb_process_event();
    }
}


static inline void nb_wait_all()
{
    assert(nb_count_event-nb_count_event_processed >= 0);

    while (nb_count_event-nb_count_event_processed > 0) {
        nb_process_event();
    }

#ifndef NDEBUG
    {
        size_t i=0;
        for (i=0; i<nb_max_outstanding; ++i) {
            assert(nb_state[i].send == 0);
            assert(nb_state[i].reply == 0);
            assert(nb_state[i].ack == 0);
        }
        assert(nb_count_event == nb_count_event_processed);
        assert(nb_count_send == nb_count_send_processed);
        assert(nb_count_reply == nb_count_reply_processed);
        assert(nb_count_ack == nb_count_ack_processed);
    }
#endif
}


static inline void nb_put(void *src, void *dst, int bytes, int proc, nb_t *nb)
{
    int status = 0;
    ptl_process_t peer;

    assert(NULL != src);
    assert(NULL != dst);
    assert(bytes > 0);
    assert(proc >= 0);
    assert(proc < l_state.size);
    assert(NULL != nb);

    if (l_state.rank == proc) {
        memcpy(dst, src, bytes);
        return;
    }

    peer.rank = proc;

    nb_wait_for_event_space();
    nb_count_event += 2;
    nb_count_send += 1;
    nb_count_ack += 1;
    nb->send += 1;
    nb->ack += 1;

    status = PtlPut(l_state.ptl_md_handle,
            (ptl_size_t)src,
            (ptl_size_t)bytes,
            PTL_ACK_REQ,
            peer,
            l_state.ptl_pt_indexes[proc],
            (ptl_match_bits_t)0,
            (ptl_size_t)dst,
            (void*)nb,
            (ptl_hdr_data_t)0);
    CHECK_PTL_RETVAL(status);
}


static inline void nb_get(void *src, void *dst, int bytes, int proc, nb_t *nb)
{
    int status = 0;
    ptl_process_t peer;

    assert(NULL != src);
    assert(NULL != dst);
    assert(bytes > 0);
    assert(proc >= 0);
    assert(proc < l_state.size);
    assert(NULL != nb);

    if (l_state.rank == proc) {
        memcpy(dst, src, bytes);
        return;
    }

    peer.rank = proc;

    nb_wait_for_event_space();
    nb_count_event += 1;
    nb_count_reply += 1;
    nb->reply += 1;

    status = PtlGet(l_state.ptl_md_handle,
            (ptl_size_t)dst,
            (ptl_size_t)bytes,
            peer,
            l_state.ptl_pt_indexes[proc],
            (ptl_match_bits_t)0,
            (ptl_size_t)src,
            (void*)nb);
    CHECK_PTL_RETVAL(status);
}


static inline void nb_acc(int datatype, void *scale,
        void *src, void *dst, int bytes, int proc, nb_t *nb)
{
    int status = 0;
    ptl_process_t peer;
    ptl_datatype_t ptl_datatype;
    char need_scale = 0;
    void *send_buf = NULL;

    assert(NULL != src);
    assert(NULL != dst);
    assert(bytes > 0);
    assert(proc >= 0);
    assert(proc < l_state.size);
    assert(NULL != nb);

    /* if max_atomic_size is surpassed, break up the transfer into strides */
    if (bytes > l_state.ptl_ni_limits.max_atomic_size) {
        int i = 0;
        int levels = 1;
        int xdim = l_state.ptl_ni_limits.max_atomic_size;
        int ydim = bytes / xdim;
        int remainder = bytes % xdim;
        for (i=0; i<ydim; ++i) {
            nb_acc(datatype, scale,
                    &((char*)src)[i*xdim],
                    &((char*)dst)[i*xdim],
                    xdim, proc, nb);
        }
        /* we may have bytes remaining */
        if (remainder != 0) {
            nb_acc(datatype, scale,
                   &((char*)src)[ydim*xdim],
                   &((char*)dst)[ydim*xdim],
                   remainder, proc, nb);
        }
        return;
    }

    assert(bytes <= l_state.ptl_ni_limits.max_atomic_size);

    /* allocate accumulate buffer on first use */
    if (NULL == acc_buf) {
        acc_buf = comex_malloc_local(l_state.ptl_ni_limits.max_atomic_size);
    }
    
    peer.rank = proc;

    nb_wait_for_event_space();
    nb_count_event += 2;
    nb_count_send += 1;
    nb_count_ack += 1;
    nb->send += 1;
    nb->ack += 1;

    switch (datatype) {
        case COMEX_ACC_INT: ptl_datatype = PTL_INT32_T; break;
        case COMEX_ACC_DBL: ptl_datatype = PTL_DOUBLE; break;
        case COMEX_ACC_FLT: ptl_datatype = PTL_FLOAT; break;
        case COMEX_ACC_CPL: ptl_datatype = PTL_FLOAT_COMPLEX; break;
        case COMEX_ACC_DCP: ptl_datatype = PTL_DOUBLE_COMPLEX; break;
        case COMEX_ACC_LNG: ptl_datatype = PTL_INT64_T; break;
        default: assert(0);
    }

    switch (datatype) {
        case COMEX_ACC_INT: need_scale = (*((int*)scale) != 1); break;
        case COMEX_ACC_DBL: need_scale = (*((double*)scale) != 1); break;
        case COMEX_ACC_FLT: need_scale = (*((float*)scale) != 1); break;
        case COMEX_ACC_CPL: need_scale = 1; break;
        case COMEX_ACC_DCP: need_scale = 1; break;
        case COMEX_ACC_LNG: need_scale = (*((long*)scale) != 1); break;
        default: assert(0);
    }

    if (need_scale) {
        /* local scaling of src */
        _scale(datatype, bytes, acc_buf, src, scale);
        send_buf = acc_buf;
    }
    else {
        send_buf = src;
    }

    status = PtlAtomic(l_state.ptl_md_handle,
            (ptl_size_t)send_buf,
            (ptl_size_t)bytes,
            PTL_ACK_REQ,
            peer,
            l_state.ptl_pt_indexes[proc],
            (ptl_match_bits_t)0,
            (ptl_size_t)dst,
            (void*)nb,
            (ptl_hdr_data_t)NULL,
            PTL_SUM,
            ptl_datatype);
    CHECK_PTL_RETVAL(status);

    if (need_scale) {
        nb_wait_for_send(nb); /* so we can reuse acc_buf */
    }
}


static inline void nb_puts(
        void *src, int *src_stride, void *dst, int *dst_stride,
        int *count, int stride_levels, int proc, nb_t *nb)
{
    int i, j;
    long src_idx, dst_idx;  /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];
    int status = 0;

    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++) {
        n1dim *= count[i];
    }

    /* calculate the destination indices */
    src_bvalue[0] = 0; src_bvalue[1] = 0; src_bunit[0] = 1; src_bunit[1] = 1;
    dst_bvalue[0] = 0; dst_bvalue[1] = 0; dst_bunit[0] = 1; dst_bunit[1] = 1;

    for(i=2; i<=stride_levels; i++) {
        src_bvalue[i] = 0;
        dst_bvalue[i] = 0;
        src_bunit[i] = src_bunit[i-1] * count[i-1];
        dst_bunit[i] = dst_bunit[i-1] * count[i-1];
    }

    status = PtlStartBundle(l_state.ptl_ni_handle);
    CHECK_PTL_RETVAL(status);

    /* index mangling */
    for(i=0; i<n1dim; i++) {
        src_idx = 0;
        dst_idx = 0;
        for(j=1; j<=stride_levels; j++) {
            src_idx += src_bvalue[j] * src_stride[j-1];
            if((i+1) % src_bunit[j] == 0) {
                src_bvalue[j]++;
            }
            if(src_bvalue[j] > (count[j]-1)) {
                src_bvalue[j] = 0;
            }
        }

        for(j=1; j<=stride_levels; j++) {
            dst_idx += dst_bvalue[j] * dst_stride[j-1];
            if((i+1) % dst_bunit[j] == 0) {
                dst_bvalue[j]++;
            }
            if(dst_bvalue[j] > (count[j]-1)) {
                dst_bvalue[j] = 0;
            }
        }
        
        nb_put((char *)src + src_idx, (char *)dst + dst_idx,
                count[0], proc, nb);
    }

    status = PtlEndBundle(l_state.ptl_ni_handle);
    CHECK_PTL_RETVAL(status);
}


static inline void nb_gets(
        void *src, int *src_stride, void *dst, int *dst_stride,
        int *count, int stride_levels, int proc, nb_t *nb)
{
    int i, j;
    long src_idx, dst_idx;  /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];
    int status = 0;

    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++) {
        n1dim *= count[i];
    }

    /* calculate the destination indices */
    src_bvalue[0] = 0; src_bvalue[1] = 0; src_bunit[0] = 1; src_bunit[1] = 1;
    dst_bvalue[0] = 0; dst_bvalue[1] = 0; dst_bunit[0] = 1; dst_bunit[1] = 1;

    for(i=2; i<=stride_levels; i++) {
        src_bvalue[i] = 0;
        dst_bvalue[i] = 0;
        src_bunit[i] = src_bunit[i-1] * count[i-1];
        dst_bunit[i] = dst_bunit[i-1] * count[i-1];
    }

    status = PtlStartBundle(l_state.ptl_ni_handle);
    CHECK_PTL_RETVAL(status);

    for(i=0; i<n1dim; i++) {
        src_idx = 0;
        for(j=1; j<=stride_levels; j++) {
            src_idx += src_bvalue[j] * src_stride[j-1];
            if((i+1) % src_bunit[j] == 0) {
                src_bvalue[j]++;
            }
            if(src_bvalue[j] > (count[j]-1)) {
                src_bvalue[j] = 0;
            }
        }

        dst_idx = 0;
        
        for(j=1; j<=stride_levels; j++) {
            dst_idx += dst_bvalue[j] * dst_stride[j-1];
            if((i+1) % dst_bunit[j] == 0) {
                dst_bvalue[j]++;
            }
            if(dst_bvalue[j] > (count[j]-1)) {
                dst_bvalue[j] = 0;
            }
        }
        
        nb_get((char *)src + src_idx, (char *)dst + dst_idx,
                count[0], proc, nb);
    }

    status = PtlEndBundle(l_state.ptl_ni_handle);
    CHECK_PTL_RETVAL(status);
}


static inline void nb_accs(
        int datatype, void *scale,
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels,
        int proc, nb_t *nb)
{
    int i, j;
    long src_idx, dst_idx;  /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];
    int status = 0;

    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++) {
        n1dim *= count[i];
    }

    /* calculate the destination indices */
    src_bvalue[0] = 0; src_bvalue[1] = 0; src_bunit[0] = 1; src_bunit[1] = 1;
    dst_bvalue[0] = 0; dst_bvalue[1] = 0; dst_bunit[0] = 1; dst_bunit[1] = 1;

    for(i=2; i<=stride_levels; i++) {
        src_bvalue[i] = 0;
        dst_bvalue[i] = 0;
        src_bunit[i] = src_bunit[i-1] * count[i-1];
        dst_bunit[i] = dst_bunit[i-1] * count[i-1];
    }

    status = PtlStartBundle(l_state.ptl_ni_handle);
    CHECK_PTL_RETVAL(status);

    /* index mangling */
    for(i=0; i<n1dim; i++) {
        src_idx = 0;
        dst_idx = 0;
        for(j=1; j<=stride_levels; j++) {
            src_idx += src_bvalue[j] * src_stride[j-1];
            if((i+1) % src_bunit[j] == 0) {
                src_bvalue[j]++;
            }
            if(src_bvalue[j] > (count[j]-1)) {
                src_bvalue[j] = 0;
            }
        }

        for(j=1; j<=stride_levels; j++) {
            dst_idx += dst_bvalue[j] * dst_stride[j-1];
            if((i+1) % dst_bunit[j] == 0) {
                dst_bvalue[j]++;
            }
            if(dst_bvalue[j] > (count[j]-1)) {
                dst_bvalue[j] = 0;
            }
        }
        
        nb_acc(datatype, scale, (char *)src + src_idx, (char *)dst + dst_idx,
                count[0], proc, nb);
    }

    status = PtlEndBundle(l_state.ptl_ni_handle);
    CHECK_PTL_RETVAL(status);
}


static inline void nb_putv(
        comex_giov_t *iov, int iov_len,
        int proc, nb_t *nb)
{
    int i = 0;
    int status = 0;

    status = PtlStartBundle(l_state.ptl_ni_handle);
    CHECK_PTL_RETVAL(status);

    for (i=0; i<iov_len; ++i) {
        int j;
        void **src = iov[i].src;
        void **dst = iov[i].dst;
        int bytes = iov[i].bytes;
        int limit = iov[i].count;
        for (j=0; j<limit; ++j) {
            nb_put(src[j], dst[j], bytes, proc, nb);
        }
    }

    status = PtlEndBundle(l_state.ptl_ni_handle);
    CHECK_PTL_RETVAL(status);
}


static inline void nb_getv(
        comex_giov_t *iov, int iov_len,
        int proc, nb_t *nb)
{
    int i = 0;
    int status = 0;

    status = PtlStartBundle(l_state.ptl_ni_handle);
    CHECK_PTL_RETVAL(status);

    for (i=0; i<iov_len; ++i) {
        int j;
        void **src = iov[i].src;
        void **dst = iov[i].dst;
        int bytes = iov[i].bytes;
        int limit = iov[i].count;
        for (j=0; j<limit; ++j) {
            nb_get(src[j], dst[j], bytes, proc, nb);
        }
    }

    status = PtlEndBundle(l_state.ptl_ni_handle);
    CHECK_PTL_RETVAL(status);
}


static inline void nb_accv(
        int datatype, void *scale,
        comex_giov_t *iov, int iov_len,
        int proc, nb_t *nb)
{
    int i = 0;
    int status = 0;

    status = PtlStartBundle(l_state.ptl_ni_handle);
    CHECK_PTL_RETVAL(status);

    for (i=0; i<iov_len; ++i) {
        int j;
        void **src = iov[i].src;
        void **dst = iov[i].dst;
        int bytes = iov[i].bytes;
        int limit = iov[i].count;
        for (j=0; j<limit; ++j) {
            nb_acc(datatype, scale, src[j], dst[j], bytes, proc, nb);
        }
    }

    status = PtlEndBundle(l_state.ptl_ni_handle);
    CHECK_PTL_RETVAL(status);
}


