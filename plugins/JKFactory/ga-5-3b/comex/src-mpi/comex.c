#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* C and/or system headers */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>

/* 3rd party headers */
#include <mpi.h>

/* our headers */
#include "comex.h"
#include "comex_impl.h"
#include "groups.h"
#include "acc.h"

#define MEMSET_AFTER_MALLOC 0
#define DEBUG 0
#define DEBUG_TO_FILE 0
#if DEBUG_TO_FILE
#   define printf(...) fprintf(my_file, __VA_ARGS__); fflush(my_file);
#endif


/* exported state */
local_state l_state;
int comex_me=-1;
int comex_nproc=-1;

/* static state */
static int initialized=0; /* for comex_initialized(), 0=false */
static char *fence_array=NULL;
#if DEBUG_TO_FILE
static FILE *my_file=NULL;
#endif

/* static function declarations */
static void  _mq_push(int dest, char *message, int count);
static int   _mq_test(void);
static void  _mq_pop(void);
static void  _gq_push(char *notify);
static void  _gq_push_item(get_t *item);
static void  _make_progress_if_needed(void);
static int   _my_isend(void *buf, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm, MPI_Request *request);
static void  _my_free(void *ptr);
static void* _my_malloc(size_t size);
static void* _my_memcpy(void *dest, const void *src, size_t n);
static int   _get_nbi(void *src, void *dst, int bytes, int proc);
static int   _put_nbi(void *src, void *dst, int bytes, int proc);
static void  _put_handler(header_t *header, char *payload);
static void  _get_request_handler(header_t *header, int proc);
static void  _get_response_handler(header_t *header, char *payload);
static void  _acc_handler(header_t *header, char *payload);
static void  _fence_request_handler(header_t *header, int proc);
static void  _fence_response_handler(header_t *header, int proc);
static void  _barrier_request_handler(header_t *header, int proc);
static void  _barrier_response_handler(header_t *header, int proc);
static void  _send_fence_message(int proc, int *notify);
static void  _fetch_and_add_request_handler(header_t *header, char *payload, int proc);
static void  _fetch_and_add_response_handler(header_t *header, char *payload);
static void  _swap_request_handler(header_t *header, char *payload, int proc);
static void  _swap_response_handler(header_t *header, char *payload);
static void  _lock_request_handler(header_t *header, int proc);
static void  _lock_response_handler(header_t *header);
static void  _unlock_request_handler(header_t *header, int proc);
static void  _lq_push(int rank, int id, char *notify);
static int   _lq_progress(void);


int _my_isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
        MPI_Comm comm, MPI_Request *request)
{
#if DEBUG
    {
        char *op;
        switch (((header_t*)buf)->operation) {
            case OP_PUT:                op="OP_PUT"; break;
            case OP_GET_REQUEST:        op="OP_GET_REQUEST"; break;
            case OP_GET_RESPONSE:       op="OP_GET_RESPONSE"; break;
            case OP_ACC_INT:            op="OP_ACC_INT"; break;
            case OP_ACC_DBL:            op="OP_ACC_DBL"; break;
            case OP_ACC_FLT:            op="OP_ACC_FLT"; break;
            case OP_ACC_CPL:            op="OP_ACC_CPL"; break;
            case OP_ACC_DCP:            op="OP_ACC_DCP"; break;
            case OP_ACC_LNG:            op="OP_ACC_LNG"; break;
            case OP_FENCE_REQUEST:      op="OP_FENCE_REQUEST"; break;
            case OP_FENCE_RESPONSE:     op="OP_FENCE_RESPONSE"; break;
            case OP_BARRIER_REQUEST:    op="OP_BARRIER_REQUEST"; break;
            case OP_BARRIER_RESPONSE:   op="OP_BARRIER_RESPONSE"; break;
            case OP_FETCH_AND_ADD_REQUEST:
                                        op="OP_FETCH_AND_ADD_REQUEST"; break;
            case OP_FETCH_AND_ADD_RESPONSE:
                                        op="OP_FETCH_AND_ADD_RESPONSE"; break;
            case OP_SWAP_REQUEST:       op="OP_SWAP_REQUEST"; break;
            case OP_SWAP_RESPONSE:      op="OP_SWAP_RESPONSE"; break;
            case OP_LOCK_REQUEST:       op="OP_LOCK_REQUEST"; break;
            case OP_LOCK_RESPONSE:      op="OP_LOCK_RESPONSE"; break;
            case OP_UNLOCK:             op="OP_UNLOCK"; break;
            default: assert(0);
        }
        printf("[%d] sending %s to %d\n", l_state.rank, op, dest);
    }
#endif
    assert(dest != MPI_ANY_SOURCE);
    return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
}


static void* _my_memcpy(void *dest, const void *src, size_t n)
{
#if DEBUG
    printf("[%d] _my_memcpy(dest=%p, src=%p, n=%zu)\n",
            l_state.rank, dest, src, n);
#endif
    return memcpy(dest, src, n);
}


static void* _my_malloc(size_t size)
{
    void *memptr=NULL;

    memptr = malloc(size);
    assert(memptr);
#if MEMSET_AFTER_MALLOC
    /* valgrind was reporting uninitialized values coming from _my_malloc
     * so use this memset to quiet valgrind for now. We may turn this
     * off after additional development/testing if we find the true
     * source of the uninitialized bytes. */
    memptr = memset(memptr, 0, size);
    assert(memptr);
#endif

#if DEBUG
    printf("[%d] _my_malloc(%zu) -> %p\n", l_state.rank, size, memptr);
#endif

    return memptr;
}


static void _my_free(void *ptr)
{
#if DEBUG
    printf("[%d] _my_free(%p)\n", l_state.rank, ptr);
#endif

    free(ptr);
}
#if DEBUG
#   define _my_free(ARG) _my_free(ARG); printf(#ARG)
#endif


static void _mq_push(int dest, char *message, int count)
{
    message_t *mq = NULL;

    mq = _my_malloc(sizeof(message_t));
    assert(mq);
    mq->next = NULL;
    mq->dest = dest;
    mq->message = message;
    _my_isend(message, count, MPI_CHAR, dest, COMEX_TAG,
            l_state.world_comm, &(mq->request));
    if (l_state.mq_tail) {
        l_state.mq_tail->next = mq;
        l_state.mq_tail = mq;
    }
    else {
        l_state.mq_head = mq;
        l_state.mq_tail = mq;
    }

    ++l_state.mq_size;
    assert(l_state.mq_size >= 0);
}


static int _mq_test()
{
    int flag;
    int rc;
    MPI_Status status;

    assert(l_state.mq_head);
    assert(l_state.mq_tail);
    assert(l_state.mq_head->request != MPI_REQUEST_NULL);

    rc = MPI_Test(&(l_state.mq_head->request), &flag, &status);
    assert(MPI_SUCCESS == rc);

    return flag;
}


static void _mq_pop()
{
    message_t *mq_to_delete = NULL;

    assert(l_state.mq_head);
    assert(l_state.mq_tail);
    assert(l_state.mq_head->request == MPI_REQUEST_NULL);

    mq_to_delete = l_state.mq_head;
    l_state.mq_head = l_state.mq_head->next;

    /* don't need message or request any longer */
    _my_free(mq_to_delete->message);
    _my_free(mq_to_delete);

    if (NULL == l_state.mq_head) {
        l_state.mq_tail = NULL;
    }

    --l_state.mq_size;
    assert(l_state.mq_size >= 0);
}


static void _gq_push(char *notify)
{
    get_t *item = _my_malloc(sizeof(get_t));
    item->next = NULL;
    item->notify_address = notify;
    _gq_push_item(item);
}


static void _gq_push_item(get_t *item)
{
    assert(item);
    if (l_state.gq_tail) {
        l_state.gq_tail->next = item;
        l_state.gq_tail = item;
    }
    else {
        l_state.gq_head = item;
        l_state.gq_tail = item;
    }
}


/* iterate over get queue and pop completed get requests
 * return true if there are get requests remaining */
static int _gq_test()
{
    get_t *old_gq = l_state.gq_head;
    l_state.gq_head = NULL;
    l_state.gq_tail = NULL;

    while (old_gq) {
        if (*((char*)old_gq->notify_address)) {
            _gq_push_item(old_gq);
            old_gq = old_gq->next;
        }
        else {
            get_t *need_to_free = old_gq;
            old_gq = old_gq->next;
            _my_free((char*)need_to_free->notify_address);
            _my_free(need_to_free);
        }
        if (l_state.gq_tail) {
            l_state.gq_tail->next = NULL;
        }
    }

    return NULL != l_state.gq_tail;
}


static void _make_progress_if_needed(void)
{
#if DEBUG
    printf("[%d] _make_progress_if_needed()\n", l_state.rank);
#endif

    if (l_state.mq_size >= MAX_NB_OUTSTANDING) {
        comex_make_progress();
    }
}


void comex_make_progress(void)
{
    int iprobe_flag=0;
    MPI_Status iprobe_status;
    int get_flag=0;
    int lock_flag=0;

#if DEBUG
    printf("[%d] comex_make_progress()\n", l_state.rank);
#endif

    do {
        iprobe_flag = 0;
        get_flag = 0;

        /* test for outgoing message completion only if we have outgoing
         * messages; we test from the head only */
        if (l_state.mq_head) {
            while (l_state.mq_head && _mq_test()) {
                _mq_pop();
            }
        }

        /* test for incoming get responses */
        if (_gq_test()) {
            get_flag = 1;
        }

        /* test for incoming messages */
        MPI_Iprobe(MPI_ANY_SOURCE, COMEX_TAG, l_state.world_comm,
                &iprobe_flag, &iprobe_status);
        if (iprobe_flag) {
            int length;
            char *message;
            char *payload;
            header_t *header;
            MPI_Status recv_status;

            /* allocate message buffer and get message */
            MPI_Get_count(&iprobe_status, MPI_CHAR, &length);
#if DEBUG
            printf("[%d] iprobe source=%d length=%d\n",
                    l_state.rank, iprobe_status.MPI_SOURCE, length);
#endif

            message = _my_malloc(length);
            MPI_Recv(message, length, MPI_CHAR,
                    iprobe_status.MPI_SOURCE, COMEX_TAG,
                    l_state.world_comm, &recv_status);
            header = (header_t*)message;
            payload = message + sizeof(header_t);
            /* dispatch message handler */
            switch (header->operation) {
                case OP_PUT:
                    _put_handler(header, payload);
                    break;
                case OP_GET_REQUEST:
                    _get_request_handler(header, iprobe_status.MPI_SOURCE);
                    break;
                case OP_GET_RESPONSE:
                    _get_response_handler(header, payload);
                    break;
                case OP_ACC_INT:
                case OP_ACC_DBL:
                case OP_ACC_FLT:
                case OP_ACC_CPL:
                case OP_ACC_DCP:
                case OP_ACC_LNG:
                    _acc_handler(header, payload);
                    break;
                case OP_FENCE_REQUEST:
                    _fence_request_handler(header, iprobe_status.MPI_SOURCE);
                    break;
                case OP_FENCE_RESPONSE:
                    _fence_response_handler(header, iprobe_status.MPI_SOURCE);
                    break;
                case OP_BARRIER_REQUEST:
                    _barrier_request_handler(header, iprobe_status.MPI_SOURCE);
                    break;
                case OP_BARRIER_RESPONSE:
                    _barrier_response_handler(header, iprobe_status.MPI_SOURCE);
                    break;
                case OP_FETCH_AND_ADD_REQUEST:
                    _fetch_and_add_request_handler(header, payload,
                            iprobe_status.MPI_SOURCE);
                    break;
                case OP_FETCH_AND_ADD_RESPONSE:
                    _fetch_and_add_response_handler(header, payload);
                    break;
                case OP_SWAP_REQUEST:
                    _swap_request_handler(header, payload,
                            iprobe_status.MPI_SOURCE);
                    break;
                case OP_SWAP_RESPONSE:
                    _swap_response_handler(header, payload);
                    break;
                case OP_LOCK_REQUEST:
                    _lock_request_handler(header, iprobe_status.MPI_SOURCE);
                    break;
                case OP_LOCK_RESPONSE:
                    _lock_response_handler(header);
                    break;
                case OP_UNLOCK:
                    _unlock_request_handler(header, iprobe_status.MPI_SOURCE);
                    break;
                default:
                    printf("[%d] header operation not recognized: %d\n",
                            l_state.rank, header->operation);
                    assert(0);
            }

            /* free the message buffer */
            _my_free(message);
        }

        lock_flag = _lq_progress();

        /* loop until we run out of incoming and outgoing messages */
#if DEBUG
        printf("[%d] iprobe_flag=%d || get_flag=%d || l_state.mq_size=%d lock_flag=%d\n",
                l_state.rank, iprobe_flag, get_flag, l_state.mq_size, lock_flag);
#endif
    } while (iprobe_flag || get_flag || l_state.mq_size);
}


static void _put_handler(header_t *header, char *payload)
{
#if DEBUG
    printf("[%d] _put_handler rem=%p loc=%p len=%d not=%p\n",
            l_state.rank,
            header->remote_address,
            header->local_address,
            header->length,
            header->notify_address);
#endif

    assert(OP_PUT == header->operation);
    _my_memcpy(header->remote_address, payload, header->length);
}


static void _get_request_handler(header_t *request_header, int proc)
{
    char *message;

#if DEBUG
    printf("[%d] _get_request_handler proc=%d\n", l_state.rank, proc);
#endif
#if DEBUG
    printf("[%d] request_header rem=%p loc=%p len=%d not=%p\n",
            l_state.rank,
            request_header->remote_address,
            request_header->local_address,
            request_header->length,
            request_header->notify_address);
#endif

    assert(OP_GET_REQUEST == request_header->operation);
    
    /* reuse the header, just change the OP */
    /* NOTE: it gets deleted anyway when we return from this handler */
    request_header->operation = OP_GET_RESPONSE;

    message = _my_malloc(sizeof(header_t) + request_header->length);
    _my_memcpy(message, request_header, sizeof(header_t));
    _my_memcpy(message+sizeof(header_t), request_header->remote_address,
            request_header->length);

    _mq_push(proc, message, sizeof(header_t) + request_header->length);
}


static void _get_response_handler(header_t *header, char *payload)
{
#if DEBUG
    printf("[%d] _get_response_handler rem=%p loc=%p len=%d not=%p\n",
            l_state.rank,
            header->remote_address,
            header->local_address,
            header->length,
            header->notify_address);
#endif

    assert(OP_GET_RESPONSE == header->operation);
    _my_memcpy(header->local_address, payload, header->length);
    *((char*)(header->notify_address)) = 0;
}


static void _acc_handler(header_t *header, char *payload)
{
    int op;
    int sizeof_scale;

#if DEBUG
    printf("[%d] _acc_handler\n", l_state.rank);
#endif

    switch (header->operation) {
        case OP_ACC_INT: sizeof_scale = sizeof(int);
                         op = COMEX_ACC_INT;
                         break;
        case OP_ACC_DBL: sizeof_scale = sizeof(double);
                         op = COMEX_ACC_DBL;
                         break;
        case OP_ACC_FLT: sizeof_scale = sizeof(float);
                         op = COMEX_ACC_FLT;
                         break;
        case OP_ACC_LNG: sizeof_scale = sizeof(long);
                         op = COMEX_ACC_LNG;
                         break;
        case OP_ACC_CPL: sizeof_scale = sizeof(SingleComplex);
                         op = COMEX_ACC_CPL;
                         break;
        case OP_ACC_DCP: sizeof_scale = sizeof(DoubleComplex);
                         op = COMEX_ACC_DCP;
                         break;
        default: assert(0);
    }
    _acc(op, header->length, header->remote_address,
            payload+sizeof_scale, payload);
}


static void _fence_request_handler(header_t *header, int proc)
{
    header_t *response_header = NULL;

#if DEBUG
    printf("[%d] _fence_request_handler proc=%d\n", l_state.rank, proc);
#endif

    /* preconditions */
    assert(header);

    /* allocate new header */
    response_header = _my_malloc(sizeof(header_t));
    assert(response_header);

    /* we make a copy since the original header will be free'd upon return */
    memcpy(response_header, header, sizeof(header_t));

    /* the only thing that changes is the operation */
    response_header->operation = OP_FENCE_RESPONSE;

    /* we send the header back to the originating proc */
    _mq_push(proc, (char*)response_header, sizeof(header_t));
}


static void _fence_response_handler(header_t *header, int proc)
{
#if DEBUG
    printf("[%d] _fence_response_handler proc=%d decrementing %d\n",
            l_state.rank, proc, (int)(*((int*)header->notify_address)));
#endif
    assert(header);
    --(*((int*)(header->notify_address)));
}


/* barrier requests are always queued for later */
static void _barrier_request_handler(header_t *header, int proc)
{
    barrier_t *new_barrier_request = NULL;

#if DEBUG
    printf("[%d] _barrier_request_handler proc=%d\n", l_state.rank, proc);
#endif

    /* create new queued barrier request */
    new_barrier_request = _my_malloc(sizeof(barrier_t));
    new_barrier_request->next = NULL;
    new_barrier_request->world_rank = proc;
    new_barrier_request->notify_address = header->notify_address;

    if (NULL == l_state.bq_head) {
        l_state.bq_head = new_barrier_request;
        l_state.bq_tail = new_barrier_request;
    }
    else {
        l_state.bq_tail->next = new_barrier_request;
        l_state.bq_tail = new_barrier_request;
    }
}


static void _barrier_response_handler(header_t *header, int proc)
{
#if DEBUG
    printf("[%d] _barrier_response_handler proc=%d\n", l_state.rank, proc);
#endif

    ++(*((char*)header->notify_address));
}


int comex_put(void *src, void *dst, int bytes, int proc, comex_group_t group)
{
    CHECK_GROUP(group,proc);
#if DEBUG
    printf("[%d] comex_put(src=%p, dst=%p, bytes=%d, proc=%d, group=%d)\n",
            l_state.rank, src, dst, bytes, proc, group);
#endif
    _put_nbi(src, dst, bytes, proc);
    comex_wait_proc(proc, group);
    return COMEX_SUCCESS;
}


int comex_get(void *src, void *dst, int bytes, int proc, comex_group_t group)
{
    CHECK_GROUP(group,proc);
#if DEBUG
    printf("[%d] comex_get(src=%p, dst=%p, bytes=%d, proc=%d)\n",
            l_state.rank, src, dst, bytes, proc);
#endif
    _get_nbi(src, dst, bytes, proc);
    comex_wait_proc(proc, group);
    return COMEX_SUCCESS;
}


static int _put_nbi(void *src, void *dst, int bytes, int proc)
{
    header_t header;
    char *message;

#if DEBUG
    printf("[%d] _put_nbi(src=%p, dst=%p, bytes=%d, proc=%d)\n",
            l_state.rank, src, dst, bytes, proc);
#endif

    /* Corner case */
    if (proc == l_state.rank) {
        _my_memcpy(dst, src, bytes);
        return COMEX_SUCCESS;
    }

    header.operation = OP_PUT;
    header.remote_address = dst;
    header.local_address = src;
    header.length = bytes;
    header.notify_address = NULL;

    message = _my_malloc(sizeof(header_t) + bytes);
    _my_memcpy(message, &header, sizeof(header_t));
    _my_memcpy(message+sizeof(header_t), src, bytes);

    fence_array[proc] = 1;

    _mq_push(proc, message, sizeof(header_t) + bytes);
    _make_progress_if_needed();

    return COMEX_SUCCESS;
}


/* TODO: need data type as a parameter before using this function! */
static int _acc_nbi(void *src, void *dst, int bytes, int proc,
        int datatype, void *scale)
{
    header_t header;
    char *message;
    int message_size = 0;
    int scale_size = 0;

#if DEBUG
    printf("[%d] _acc_nbi(scale=%p src=%p, dst=%p, bytes=%d, proc=%d)\n",
            l_state.rank, scale, src, dst, bytes, proc);
#endif

    switch (datatype) {
        case COMEX_ACC_INT:
            header.operation = OP_ACC_INT;
            scale_size = sizeof(int);
            break;
        case COMEX_ACC_DBL:
            header.operation = OP_ACC_DBL;
            scale_size = sizeof(double);
            break;
        case COMEX_ACC_FLT:
            header.operation = OP_ACC_FLT;
            scale_size = sizeof(float);
            break;
        case COMEX_ACC_CPL:
            header.operation = OP_ACC_CPL;
            scale_size = sizeof(SingleComplex);
            break;
        case COMEX_ACC_DCP:
            header.operation = OP_ACC_DCP;
            scale_size = sizeof(DoubleComplex);
            break;
        case COMEX_ACC_LNG:
            header.operation = OP_ACC_LNG;
            scale_size = sizeof(long);
            break;
        default: assert(0);
    }

    /* Corner case */
    if (proc == l_state.rank) {
        _acc(datatype, bytes, dst, src, scale);
        return COMEX_SUCCESS;
    }

    header.remote_address = dst;
    header.local_address = src;
    header.length = bytes;
    header.notify_address = NULL;

    message_size = sizeof(header_t) + scale_size + bytes;
    message = _my_malloc(message_size);
    _my_memcpy(message, &header, sizeof(header_t));
    _my_memcpy(message+sizeof(header_t), scale, scale_size);
    _my_memcpy(message+sizeof(header_t)+scale_size, src, bytes);

    fence_array[proc] = 1;

    _mq_push(proc, message, message_size);
    _make_progress_if_needed();

    return COMEX_SUCCESS;
}


static int _get_nbi(void *src, void *dst, int bytes, int proc)
{
    header_t *header = NULL;

#if DEBUG
    printf("[%d] _get_nbi(src=%p, dst=%p, bytes=%d, proc=%d)\n",
            l_state.rank, src, dst, bytes, proc);
#endif

    /* Corner case */
    if (proc == l_state.rank) {
        _my_memcpy(dst, src, bytes);
        return COMEX_SUCCESS;
    }

    header = _my_malloc(sizeof(header_t));
    assert(header);
    header->operation = OP_GET_REQUEST;
    header->remote_address = src;
    header->local_address = dst;
    header->length = bytes;
    header->notify_address = _my_malloc(sizeof(char));

    /* set the value to wait on for a get response */
    *((char*)(header->notify_address)) = 1;
    _gq_push(header->notify_address);

    /* get request is a header-only message */
    _mq_push(proc, (char*)header, sizeof(header_t));
    _make_progress_if_needed();

    return COMEX_SUCCESS;
}


int comex_puts(
        void *src_ptr, int *src_stride_ar,
        void *dst_ptr, int *dst_stride_ar,
        int *count, int stride_levels,
        int proc, comex_group_t group)
{
    int i, j;
    long src_idx, dst_idx;  /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];

    CHECK_GROUP(group,proc);
#if DEBUG
    if (stride_levels) {
        printf("[%d] comex_puts(src_ptr=%p, src_stride_ar[0]=%d, dst_ptr=%p, dst_stride_ar[0]=%d, count[0]=%d, stride_levels=%d, proc=%d)\n",
                l_state.rank, src_ptr, src_stride_ar[0], dst_ptr, dst_stride_ar[0], count[0], stride_levels, proc);
    }
    else {
        printf("[%d] comex_puts(src_ptr=%p, src_stride_ar=NULL, dst_ptr=%p, dst_stride_ar=NULL, count=NULL, stride_levels=%d, proc=%d)\n",
                l_state.rank, src_ptr, dst_ptr, stride_levels, proc);
    }
#endif

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

    /* index mangling */
    for(i=0; i<n1dim; i++) {
        src_idx = 0;
        dst_idx = 0;
        for(j=1; j<=stride_levels; j++) {
            src_idx += src_bvalue[j] * src_stride_ar[j-1];
            if((i+1) % src_bunit[j] == 0) {
                src_bvalue[j]++;
            }
            if(src_bvalue[j] > (count[j]-1)) {
                src_bvalue[j] = 0;
            }
        }

        for(j=1; j<=stride_levels; j++) {
            dst_idx += dst_bvalue[j] * dst_stride_ar[j-1];
            if((i+1) % dst_bunit[j] == 0) {
                dst_bvalue[j]++;
            }
            if(dst_bvalue[j] > (count[j]-1)) {
                dst_bvalue[j] = 0;
            }
        }
        
        _put_nbi((char *)src_ptr + src_idx, 
                (char *)dst_ptr + dst_idx, count[0], proc);
    }

    comex_wait_proc(proc, group);

    return COMEX_SUCCESS;
}


int comex_gets(
        void *src_ptr, int *src_stride_ar,
        void *dst_ptr, int *dst_stride_ar,
        int *count, int stride_levels,
        int proc, comex_group_t group)
{
    int i, j;
    long src_idx, dst_idx;  /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];

    CHECK_GROUP(group,proc);
#if DEBUG
    if (stride_levels) {
        printf("[%d] comex_gets(src_ptr=%p, src_stride_ar[0]=%d, dst_ptr=%p, dst_stride_ar[0]=%d, count[0]=%d, stride_levels=%d, proc=%d)\n",
                l_state.rank, src_ptr, src_stride_ar[0], dst_ptr, dst_stride_ar[0], count[0], stride_levels, proc);
    }
    else {
        printf("[%d] comex_gets(src_ptr=%p, src_stride_ar=NULL, dst_ptr=%p, dst_stride_ar=NULL, count=NULL, stride_levels=%d, proc=%d)\n",
                l_state.rank, src_ptr, dst_ptr, stride_levels, proc);
    }
#endif

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

    for(i=0; i<n1dim; i++) {
        src_idx = 0;
        for(j=1; j<=stride_levels; j++) {
            src_idx += src_bvalue[j] * src_stride_ar[j-1];
            if((i+1) % src_bunit[j] == 0) {
                src_bvalue[j]++;
            }
            if(src_bvalue[j] > (count[j]-1)) {
                src_bvalue[j] = 0;
            }
        }

        dst_idx = 0;
        
        for(j=1; j<=stride_levels; j++) {
            dst_idx += dst_bvalue[j] * dst_stride_ar[j-1];
            if((i+1) % dst_bunit[j] == 0) {
                dst_bvalue[j]++;
            }
            if(dst_bvalue[j] > (count[j]-1)) {
                dst_bvalue[j] = 0;
            }
        }
        
        _get_nbi((char *)src_ptr + src_idx, 
                (char *)dst_ptr + dst_idx, count[0], proc);
    }
    
    comex_wait_proc(proc, group);
    
    return COMEX_SUCCESS;
}


int comex_acc(
        int datatype, void *scale,
        void *src_ptr, void *dst_ptr, int bytes,
        int proc, comex_group_t group)
{

    comex_accs(datatype, scale, src_ptr, NULL, dst_ptr, 
            NULL, &bytes, 0, proc, group);
    return COMEX_SUCCESS;
}


int comex_accs(
        int datatype, void *scale,
        void *src_ptr, int *src_stride_ar,
        void *dst_ptr, int *dst_stride_ar,
        int *count, int stride_levels,
        int proc, comex_group_t group)
{
    int i, j;
    long src_idx, dst_idx;  /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];

    CHECK_GROUP(group,proc);
#if DEBUG
    if (stride_levels) {
        printf("[%d] comex_accs(src_ptr=%p, src_stride_ar[0]=%d, dst_ptr=%p, dst_stride_ar[0]=%d, count[0]=%d, stride_levels=%d, proc=%d)\n",
                l_state.rank, src_ptr, src_stride_ar[0], dst_ptr, dst_stride_ar[0], count[0], stride_levels, proc);
    }
    else {
        printf("[%d] comex_accs(src_ptr=%p, src_stride_ar=NULL, dst_ptr=%p, dst_stride_ar=NULL, count=NULL, stride_levels=%d, proc=%d)\n",
                l_state.rank, src_ptr, dst_ptr, stride_levels, proc);
    }
#endif

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

    /* index mangling */
    for(i=0; i<n1dim; i++) {
        src_idx = 0;
        dst_idx = 0;
        for(j=1; j<=stride_levels; j++) {
            src_idx += src_bvalue[j] * src_stride_ar[j-1];
            if((i+1) % src_bunit[j] == 0) {
                src_bvalue[j]++;
            }
            if(src_bvalue[j] > (count[j]-1)) {
                src_bvalue[j] = 0;
            }
        }

        for(j=1; j<=stride_levels; j++) {
            dst_idx += dst_bvalue[j] * dst_stride_ar[j-1];
            if((i+1) % dst_bunit[j] == 0) {
                dst_bvalue[j]++;
            }
            if(dst_bvalue[j] > (count[j]-1)) {
                dst_bvalue[j] = 0;
            }
        }
        
        _acc_nbi((char *)src_ptr + src_idx, 
                (char *)dst_ptr + dst_idx, count[0], proc,
                datatype, scale);
    }

    comex_wait_proc(proc, group);

    return COMEX_SUCCESS;
}


static void _send_fence_message(int proc, int *notify)
{
    header_t *header = NULL;

#if DEBUG
    printf("[%d] _send_fence_message(proc=%d notify=%p)\n",
            l_state.rank, proc, notify);
#endif

    /* create and prepare the header */
    header = _my_malloc(sizeof(header_t));
    assert(header);
    header->operation = OP_FENCE_REQUEST;
    header->remote_address = NULL;
    header->local_address = NULL;
    header->length = 0;
    header->notify_address = notify;

    _mq_push(proc, (char*)header, sizeof(header_t));
    _make_progress_if_needed();
}


int comex_fence_all(comex_group_t group)
{
    int p = 0;
    int count_before = 0;
    int count_after = 0;
    int *notify = NULL;

#if DEBUG
    printf("[%d] comex_fence_all()\n", l_state.rank);
#endif

    /* count how many fence messagse to send */
    for (p=0; p<l_state.size; ++p) {
        if (fence_array[p]) {
            ++count_before;
        }
    }

    /* check for no outstanding put/get requests */
    if (0 == count_before) {
        return COMEX_SUCCESS;
    }

    /* we allocate the notifier on the heap in order to avoid reentry
     * into this function causing the same stack address to be reused */
    notify = _my_malloc(sizeof(int));
    assert(notify);
    *notify = count_before;

    /* optimize by only sending to procs which we have outstanding messages */
    for (p=0; p<l_state.size; ++p) {
        if (p == l_state.rank) {
            continue;
        }
        if (fence_array[p]) {
            _send_fence_message(p, notify);
        }
    }

    while (*notify > 0) {
        comex_make_progress();
#if DEBUG
        printf("[%d] are we stuck? notify=%d\n", l_state.rank, *notify);
#endif
    }

    for (p=0; p<l_state.size; ++p) {
        if (fence_array[p]) {
            fence_array[p] = 0;
            ++count_after;
        }
    }
    assert(count_before == count_after);

    _my_free(notify);

    return COMEX_SUCCESS;
}


int comex_fence_proc(int proc, comex_group_t group)
{
    int *notify = NULL;

    CHECK_GROUP(group,proc);
#if DEBUG
    printf("[%d] comex_fence_proc(proc=%d,group=$d)\n",
            l_state.rank, proc, group);
#endif

    /* Corner case */
    if (proc == l_state.rank) {
        return COMEX_SUCCESS;
    }

    /* are there outstanding put/acc requests to proc? */
    if (!fence_array[proc]) {
        return COMEX_SUCCESS;
    }

    /* we allocate the notifier on the heap in order to avoid reentry
     * into this function causing the same stack address to be reused */
    notify = _my_malloc(sizeof(int));
    *notify = 1;

    _send_fence_message(proc, notify);

    while (*notify > 0) {
        comex_make_progress();
    }

    fence_array[proc] = 0;

    _my_free(notify);

    return COMEX_SUCCESS;
}


int comex_barrier(comex_group_t comex_group)
{
    MPI_Comm mpi_comm_default;
    MPI_Group mpi_group_world;
    MPI_Group mpi_group_default;
    comex_igroup_t *comex_igroup = NULL;
    int group_size = -1;
    int group_rank = -1;
    int group_proc = -1; /* used in for loop */
    int status = 0;
    char *group_mask; /* which procs are part of default group 0/1 */
    char *request_received_mask; /* which procs did we queue requests from */
    char *response_received_mask; /* which procs did we get ack from */
    int request_received_count = 0;
    int response_received_count = 0;

#if DEBUG
    printf("[%d] comex_barrier(%d)\n", l_state.rank, comex_group);
#endif

    /* get comm size and rank from group parameter */
    comex_igroup = comex_get_igroup_from_group(comex_group);
    mpi_group_default = comex_igroup->group;
    mpi_comm_default = comex_igroup->comm;
    status = MPI_Comm_size(mpi_comm_default, &group_size);
    assert(MPI_SUCCESS == status);
    status = MPI_Comm_rank(mpi_comm_default, &group_rank);
    assert(MPI_SUCCESS == status);

    /* get world group (for later rank translation) */
    status = MPI_Comm_group(l_state.world_comm, &mpi_group_world);
    assert(MPI_SUCCESS == status);

    /* create array to note which procs we should respond to */
    group_mask = _my_malloc(sizeof(char) * l_state.size);
    memset(group_mask, 0, sizeof(char) * l_state.size);

    /* create array to note which procs we've received requests from */
    request_received_mask = _my_malloc(sizeof(char) * l_state.size);
    memset(request_received_mask, 0, sizeof(char) * l_state.size);

    /* create array to note which procs we've received responses from */
    response_received_mask = _my_malloc(sizeof(char) * group_size);
    memset(response_received_mask, 0, sizeof(char) * group_size);

    /* send out barrier requests to all procs in default group */
#define USE_RING 1
#if USE_RING
    int ring_proc;
    for (ring_proc=0; ring_proc<group_size; ++ring_proc)
#else
    for (group_proc=0; group_proc<group_size; ++group_proc)
#endif
    {
        int world_rank = -1;
        header_t *header = NULL;
#if USE_RING
        group_proc = (group_rank+ring_proc)%group_size;
#endif

        status = MPI_Group_translate_ranks(
                mpi_group_default, 1, &group_proc,
                mpi_group_world, &world_rank);
        assert(MPI_SUCCESS == status);

        group_mask[world_rank] = 1;

        /* don't send barrier message to self */
        if (group_rank == group_proc) {
            request_received_mask[world_rank] = 1;
            continue;
        }

        header = _my_malloc(sizeof(header_t));
        header->operation = OP_BARRIER_REQUEST;
        header->remote_address = NULL;
        header->local_address = NULL;
        header->length = 0;
        header->notify_address = &response_received_mask[group_proc];

        /* barrier request is a header-only message */
#if DEBUG
        printf("[%d] sending barrier to %d\n", l_state.rank, world_rank);
#endif
        _mq_push(world_rank, (char*)header, sizeof(header_t));
        _make_progress_if_needed();
    }
    /* done with MPI_Group so free it */
    MPI_Group_free(&mpi_group_world);

    /* wait until we've queued all incoming requests */
    do {
        int p;
        barrier_t *qcur = NULL;

        comex_make_progress();

        qcur = l_state.bq_head;
        while (NULL != qcur) {
            if (group_mask[qcur->world_rank]) {
                request_received_mask[qcur->world_rank] = 1;
            }
            qcur = qcur->next;
        }

        request_received_count = 0;
        for (p=0; p<l_state.size; ++p) {
            if (1 == request_received_mask[p]) {
                ++request_received_count;
            }
        }
    } while (request_received_count != group_size);

    /* send response to self */
    request_received_mask[l_state.rank] = 0;
    response_received_mask[group_rank] = 1;

    /* now send out all (other) responses */
    do {
        int p;
        barrier_t *qcur = NULL;
        barrier_t *new_bq_head = NULL;
        barrier_t *new_bq_tail = NULL;

        comex_make_progress();

        qcur = l_state.bq_head;
        while (NULL != qcur) {
            if (group_mask[qcur->world_rank]
                    && 1 == request_received_mask[qcur->world_rank]) {
                barrier_t *last = NULL;
                header_t *response = NULL;

                request_received_mask[qcur->world_rank] = 0;
                response = _my_malloc(sizeof(header_t));
                response->operation = OP_BARRIER_RESPONSE;
                response->remote_address = NULL;
                response->local_address = NULL;
                response->length = 0;
                response->notify_address = qcur->notify_address;
                /* barrier response is a header-only message */
                _mq_push(qcur->world_rank, (char*)response, sizeof(header_t));
                _make_progress_if_needed();
                last = qcur;
                qcur = qcur->next;
                _my_free(last);
            }
            else {
                if (new_bq_tail) {
                    assert(NULL == new_bq_tail->next);
                    assert(new_bq_head);
                    new_bq_tail->next = qcur;
                    new_bq_tail = qcur;
                }
                else {
                    assert(NULL == new_bq_head);
                    new_bq_head = qcur;
                    new_bq_tail = qcur;
                }
                qcur = qcur->next;
                new_bq_tail->next = NULL;
            }
        }
        l_state.bq_head = new_bq_head;
        l_state.bq_tail = new_bq_tail;

        request_received_count = 0;
        for (p=0; p<l_state.size; ++p) {
            if (1 == request_received_mask[p]) {
                ++request_received_count;
            }
        }
    } while (0 != request_received_count);

    /* wait until we've received all responses */
    do {
        int p;

        comex_make_progress();

        response_received_count = 0;
        for (p=0; p<group_size; ++p) {
            if (1 == response_received_mask[p]) {
                ++response_received_count;
            }
        }
    } while (response_received_count != group_size);

    _my_free(group_mask);
    _my_free(request_received_mask);
    _my_free(response_received_mask);
#if DEBUG
    printf("[%d] comex_barrier(%d) DONE\n", l_state.rank, *comex_group);
#endif

    return COMEX_SUCCESS;
}


void comex_error(char *msg, int code)
{
    fprintf(stderr, msg);
    fprintf(stderr,"Received an Error in Communication\n");
    MPI_Abort(l_state.world_comm, code);
}


void *comex_malloc_local(size_t size)
{
    void *ptr = NULL;

    ptr = _my_malloc(size);

    return ptr;
}


int comex_free_local(void *ptr)
{
    /* preconditions */
    assert(NULL != ptr);

    /* free the memory */
    _my_free(ptr);

    return COMEX_SUCCESS;
}


int comex_init()
{
    int status;
    int init_flag;
    
    if (initialized) {
        return COMEX_SUCCESS;
    }
    initialized = 1;

    /* Assert MPI has been initialized */
    status = MPI_Initialized(&init_flag);
    assert(MPI_SUCCESS == status);
    assert(init_flag);
    
    /* Duplicate the World Communicator */
    status = MPI_Comm_dup(MPI_COMM_WORLD, &(l_state.world_comm));
    assert(MPI_SUCCESS == status);
    assert(l_state.world_comm); 

    /* My Rank */
    status = MPI_Comm_rank(l_state.world_comm, &(l_state.rank));
    assert(MPI_SUCCESS == status);
    comex_me = l_state.rank;

    /* World Size */
    status = MPI_Comm_size(l_state.world_comm, &(l_state.size));
    assert(MPI_SUCCESS == status);
    comex_nproc = l_state.size;
    
    /* groups */
    comex_group_init();

    /* static state */
    /* TODO: consider using calloc instead of malloc+memset? */
    fence_array = malloc(sizeof(char) * l_state.size);
    (void)memset(fence_array, 0, l_state.size);

#if DEBUG_TO_FILE
    {
        char pathname[80];
        sprintf(pathname, "trace.%d.log", l_state.rank);
        my_file = fopen(pathname, "w");
    }
#endif

    /* mutexes */
    l_state.mutexes = NULL;
    l_state.num_mutexes = 0;
    l_state.lq_head = NULL;
    l_state.lq_tail = NULL;

    /* Synch - Sanity Check */
    comex_barrier(COMEX_GROUP_WORLD);

    return COMEX_SUCCESS;
}


int comex_init_args(int *argc, char ***argv)
{
    int rc;
    int init_flag;
    
    MPI_Initialized(&init_flag);
    
    if(!init_flag)
        MPI_Init(argc, argv);
    
    rc = comex_init();
    return rc;
}


int comex_initialized()
{
    return initialized;
}


int comex_finalize()
{
    /* it's okay to call multiple times -- extra calls are no-ops */
    if (!initialized) {
        return COMEX_SUCCESS;
    }

    initialized = 0;

    /* Make sure that all outstanding operations are done */
    comex_barrier(COMEX_GROUP_WORLD);

    /* groups */
    comex_group_finalize();

    free(fence_array);

    /* destroy the communicators */
    MPI_Comm_free(&l_state.world_comm);

    return COMEX_SUCCESS;
}


int comex_nbput(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    int rc;
    rc = comex_put(src, dst, bytes, proc, group);
    return rc;
}


int comex_nbget(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    int rc;
    rc = comex_get(src, dst, bytes, proc, group);
    return COMEX_SUCCESS;
}


int comex_wait_proc(int proc, comex_group_t group)
{
    comex_wait_all(group);
    return COMEX_SUCCESS;
}


int comex_wait(comex_request_t* hdl)
{
    comex_wait_all(COMEX_GROUP_WORLD);
    return COMEX_SUCCESS;
}


int comex_test(comex_request_t* hdl, int *status)
{
    comex_wait_all(COMEX_GROUP_WORLD);
    *status = 0;
    return COMEX_SUCCESS;
}


int comex_wait_all(comex_group_t group)
{
    comex_make_progress();
    return COMEX_SUCCESS;
}


int comex_nbputs(
        void *src_ptr, int *src_stride_ar,
        void *dst_ptr, int *dst_stride_ar,
        int *count, int stride_levels, 
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    int rc;
    rc = comex_puts(src_ptr, src_stride_ar, dst_ptr, 
            dst_stride_ar, count, stride_levels, proc, group);
    return rc;

}


int comex_nbgets(
        void *src_ptr, int *src_stride_ar,
        void *dst_ptr, int *dst_stride_ar,
        int *count, int stride_levels, 
        int proc, comex_group_t group,
        comex_request_t *hdl) 
{
    int rc;
    rc = comex_gets(src_ptr, src_stride_ar, dst_ptr, 
            dst_stride_ar, count,stride_levels, proc, group);
    return rc;
}


int comex_nbaccs(
        int datatype, void *scale,
        void *src_ptr, int *src_stride_ar,
        void *dst_ptr, int *dst_stride_ar,
        int *count, int stride_levels,
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    int rc;
    rc = comex_accs(datatype, scale, src_ptr, src_stride_ar,
            dst_ptr, dst_stride_ar, count, stride_levels, proc, group);
    return rc;
}


/* Vector Calls */


int comex_putv(
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group)
{
    int i;
    for (i=0; i<iov_len; ++i) {
        int j;
        void **src = iov[i].src;
        void **dst = iov[i].dst;
        int bytes = iov[i].bytes;
        int count = iov[i].count;
        for (j=0; j<count; ++j) {
            comex_put(src[j], dst[j], bytes, proc, group);
        }
    }
    return COMEX_SUCCESS;
}


int comex_getv(
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group)
{
    int i;
    for (i=0; i<iov_len; ++i) {
        int j;
        void **src = iov[i].src;
        void **dst = iov[i].dst;
        int bytes = iov[i].bytes;
        int count = iov[i].count;
        for (j=0; j<count; ++j) {
            comex_get(src[j], dst[j], bytes, proc, group);
        }
    }
    return COMEX_SUCCESS;
}


int comex_accv(
        int datatype, void *scale,
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group)
{
    int i;
    for (i=0; i<iov_len; ++i) {
        int j;
        void **src = iov[i].src;
        void **dst = iov[i].dst;
        int bytes = iov[i].bytes;
        int count = iov[i].count;
        for (j=0; j<count; ++j) {
            comex_acc(datatype, scale, src[j], dst[j], bytes, proc, group);
        }
    }
    return COMEX_SUCCESS;
}


int comex_nbputv(
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group,
        comex_request_t* handle)
{
    int rc;
    rc = comex_putv(iov, iov_len, proc, group);
    return rc;
}


int comex_nbgetv(
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group,
        comex_request_t* handle)
{
    int rc;
    rc = comex_getv(iov, iov_len, proc, group);
    return rc;
}


int comex_nbaccv(
        int datatype, void *scale,
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group,
        comex_request_t* handle)
{
    int rc;
    rc = comex_accv(datatype, scale, iov, iov_len, proc, group);
    return rc;
}


int comex_rmw(
        int comex_op, void *ploc, void *prem, int extra,
        int proc, comex_group_t group)
{
    header_t *header = NULL;
    char *notify = NULL;
    char *payload = NULL;
    char *message = NULL;
    int length = 0;
    op_t op;
    long extra_long = (long)extra;

    CHECK_GROUP(group,proc);
#if DEBUG
    printf("[%d] comex_rmw(%d, %p, %p, %d, %d)\n",
            l_state.rank, comex_op, ploc, prem, extra, proc);
#endif

    switch (comex_op) {
        case COMEX_FETCH_AND_ADD:
            op = OP_FETCH_AND_ADD_REQUEST;
            length = sizeof(int);
            payload = _my_malloc(length);
            _my_memcpy(payload, &extra, length);
            break;
        case COMEX_FETCH_AND_ADD_LONG:
            op = OP_FETCH_AND_ADD_REQUEST;
            length = sizeof(long);
            payload = _my_malloc(length);
            _my_memcpy(payload, &extra_long, length);
            break;
        case COMEX_SWAP:
            op = OP_SWAP_REQUEST;
            length = sizeof(int);
            payload = _my_malloc(length);
            _my_memcpy(payload, ploc, length);
            break;
        case COMEX_SWAP_LONG:
            op = OP_SWAP_REQUEST;
            length = sizeof(long);
            payload = _my_malloc(length);
            _my_memcpy(payload, ploc, length);
            break;
        default: assert(0);
    }

    notify = _my_malloc(sizeof(char));

    /* create and prepare the header */
    header = _my_malloc(sizeof(header_t));
    assert(header);
    header->operation = op;
    header->remote_address = prem;
    header->local_address = ploc;
    header->length = length;
    header->notify_address = notify;

    /* set the value to wait on for a rmw response */
    *notify = 1;

    message = _my_malloc(sizeof(header_t) + length);
    _my_memcpy(message, header, sizeof(header_t));
    _my_memcpy(message+sizeof(header_t), payload, length);

    _mq_push(proc, (char*)message, sizeof(header_t)+length);
    
    while (*notify > 0) {
        comex_make_progress();
#if DEBUG
        printf("[%d] are we stuck? rmw notify=%d\n",
                l_state.rank, (int)*notify);
#endif
    }

    _my_free(notify);
    _my_free(payload);
    _my_free(header);

    return COMEX_SUCCESS;
}


static void _fetch_and_add_request_handler(header_t *request_header, char *payload, int proc)
{
    char *message;

#if DEBUG
    printf("[%d] _fetch_and_add_request_handler proc=%d\n", l_state.rank, proc);
#endif
#if DEBUG
    printf("[%d] request_header rem=%p loc=%p len=%d not=%p\n",
            l_state.rank,
            request_header->remote_address,
            request_header->local_address,
            request_header->length,
            request_header->notify_address);
#endif

    assert(OP_FETCH_AND_ADD_REQUEST == request_header->operation);
    
    /* reuse the header, just change the OP */
    /* NOTE: it gets deleted anyway when we return from this handler */
    request_header->operation = OP_FETCH_AND_ADD_RESPONSE;

    /* fetch */
    message = _my_malloc(sizeof(header_t) + request_header->length);
    _my_memcpy(message, request_header, sizeof(header_t));
    _my_memcpy(message+sizeof(header_t), request_header->remote_address,
            request_header->length);

    /* "add" */
    if (sizeof(int) == request_header->length) {
        *((int*)request_header->remote_address) += *((int*)payload);
    }
    else if (sizeof(long) == request_header->length) {
        *((long*)request_header->remote_address) += *((long*)payload);
    }
    else {
        assert(0);
    }

    _mq_push(proc, message, sizeof(header_t) + request_header->length);
}


static void  _fetch_and_add_response_handler(header_t *header, char *payload)
{
#if DEBUG
    printf("[%d] _fetch_and_add_response_handler rem=%p loc=%p len=%d not=%p\n",
            l_state.rank,
            header->remote_address,
            header->local_address,
            header->length,
            header->notify_address);
#endif

    assert(OP_FETCH_AND_ADD_RESPONSE == header->operation);
    _my_memcpy(header->local_address, payload, header->length);
    *((char*)(header->notify_address)) = 0;
}


static void _swap_request_handler(header_t *header, char *payload, int proc)
{
    char *message = NULL;
    void *tmp = NULL;

#if DEBUG
    printf("[%d] _swap_request rem=%p loc=%p len=%d not=%p\n",
            l_state.rank,
            header->remote_address,
            header->local_address,
            header->length,
            header->notify_address);
#endif

    assert(OP_SWAP_REQUEST == header->operation);
    assert(header->length > 0);

    /* need a temporary for the swap */
    tmp = _my_malloc(header->length);
    _my_memcpy(tmp, payload, header->length);
    _my_memcpy(payload, header->remote_address, header->length);
    _my_memcpy(header->remote_address, tmp, header->length);
    _my_free(tmp);

    /* reuse the header, just change the OP */
    /* NOTE: it gets deleted anyway when we return from this handler */
    header->operation = OP_SWAP_RESPONSE;

    message = _my_malloc(sizeof(header_t) + header->length);
    _my_memcpy(message, header, sizeof(header_t));
    _my_memcpy(message+sizeof(header_t), payload, header->length);

    _mq_push(proc, message, sizeof(header_t) + header->length);
}


static void  _swap_response_handler(header_t *header, char *payload)
{
#if DEBUG
    printf("[%d] _swap_response_handler rem=%p loc=%p len=%d not=%p\n",
            l_state.rank,
            header->remote_address,
            header->local_address,
            header->length,
            header->notify_address);
#endif

    assert(OP_SWAP_RESPONSE == header->operation);
    _my_memcpy(header->local_address, payload, header->length);
    *((char*)(header->notify_address)) = 0;
}


static void _lock_request_handler(header_t *header, int proc)
{
#if DEBUG
    printf("[%d] _lock_request_handler id=%d proc=%d\n",
            l_state.rank, header->length, proc);
#endif

    _lq_push(proc, header->length, header->notify_address);
}


static void _lock_response_handler(header_t *header)
{
#if DEBUG
    printf("[%d] _lock_response_handler id=%d\n",
            l_state.rank, header->length);
#endif

    *((char*)header->notify_address) = 0;
}


static void _unlock_request_handler(header_t *header, int proc)
{
#if DEBUG
    printf("[%d] _unlock_request_handler id=%d proc=%d\n",
            l_state.rank, header->length, proc);
#endif

    assert(header->length >= 0);
    assert(header->length < l_state.num_mutexes);
    /* make sure the unlock request came from the lock holder */
    assert(l_state.mutexes[header->length] == proc);
    l_state.mutexes[header->length] = -1;
}


static void _lq_push(int rank, int id, char *notify)
{
    lock_t *lock = NULL;

#if DEBUG
    printf("[%d] _lq_push rank=%d id=%d\n", l_state.rank, rank, id);
#endif

    lock = _my_malloc(sizeof(lock_t));
    lock->next = NULL;
    lock->rank = rank;
    lock->id = id;
    lock->notify_address = notify;

    if (l_state.lq_tail) {
#if DEBUG
    printf("[%d] _lq_push rank=%d id=%d to tail\n", l_state.rank, rank, id);
#endif
        assert(NULL == l_state.lq_tail->next);
        l_state.lq_tail->next = lock;
        l_state.lq_tail = lock;
    }
    else {
#if DEBUG
    printf("[%d] _lq_push rank=%d id=%d to head\n", l_state.rank, rank, id);
#endif
        assert(NULL == l_state.lq_head);
        l_state.lq_head = lock;
        l_state.lq_tail = lock;
    }
}


static int _lq_progress(void)
{
    int needs_progress = 0;
    lock_t *lock = NULL;
    lock_t *new_lock_head = NULL;
    lock_t *new_lock_tail = NULL;

#if DEBUG
    if (l_state.num_mutexes > 0) {
        printf("[%d] _lq_progress mutex[0]=%d\n",
                l_state.rank, l_state.mutexes[0]);
    }
    else {
        //printf("[%d] _lq_progress no mutexes\n", l_state.rank);
    }
#endif

    lock = l_state.lq_head;
    while (lock) {
        if (l_state.mutexes[lock->id] < 0) {
            lock_t *last = NULL;
            header_t *header = NULL;

            l_state.mutexes[lock->id] = lock->rank;
#if DEBUG
    printf("[%d] _lq_progress rank=%d now holds lock %d\n",
            l_state.rank, lock->rank, lock->id);
#endif
            header = _my_malloc(sizeof(header_t));
            header->operation = OP_LOCK_RESPONSE;
            header->remote_address = NULL;
            header->local_address = NULL;
            header->length = -1;
            header->notify_address = lock->notify_address;
            _mq_push(lock->rank, (char*)header, sizeof(header_t));
            last = lock;
            lock = lock->next;
            _my_free(last);
            needs_progress = 1;
        }
        else {
#if DEBUG
    printf("[%d] _lq_progress rank=%d could not acquire lock %d\n",
            l_state.rank, lock->rank, lock->id);
#endif
            if (new_lock_tail) {
                assert(new_lock_head);
                assert(NULL == new_lock_tail->next);
                new_lock_tail->next = lock;
                new_lock_tail = lock;
            }
            else {
                assert(NULL == new_lock_head);
                new_lock_head = lock;
                new_lock_tail = lock;
            }
            lock = lock->next;
            new_lock_tail->next = NULL;
        }
    }
    l_state.lq_head = new_lock_head;
    l_state.lq_tail = new_lock_tail;

    return needs_progress;
}


/* Mutex Operations */
int comex_create_mutexes(int num)
{
    int i=0;

    assert(0 <= num);
    assert(NULL == l_state.mutexes);
    assert(0 == l_state.num_mutexes);
    assert(NULL == l_state.lq_head);
    assert(NULL == l_state.lq_tail);

    l_state.num_mutexes = num;

    if (num > 0) {
        /* create all of the mutexes */
        l_state.mutexes = (int*)_my_malloc(num * sizeof(int));
        assert(l_state.mutexes);

        /* init all of my mutexes to unlocked */
        for (i=0; i<num; ++i) {
            l_state.mutexes[i] = -1;
        }
    }

    return COMEX_SUCCESS;
}


int comex_destroy_mutexes()
{
    int number_of_outstanding_locks = 0;

    /* make sure mutexes were previously created */
    assert(0 == l_state.num_mutexes || NULL != l_state.mutexes);

    /* fix race condition -- some procs reach destroy before servicing all
     * requests */
    comex_barrier(COMEX_GROUP_WORLD);

    /* you cannot free mutexes if one is in use or queued for use */
    do {
        int m;
        comex_make_progress();
        number_of_outstanding_locks = 0;
        for (m=0; m<l_state.num_mutexes; ++m) {
            if (l_state.mutexes[m] >= 0) {
                ++number_of_outstanding_locks;
            }
        }
    } while (l_state.lq_head || number_of_outstanding_locks > 0);

    assert(NULL == l_state.lq_head);
    assert(NULL == l_state.lq_tail);
#ifndef NDEBUG
    {
        int i;
        for (i=0; i<l_state.num_mutexes; ++i) {
            assert(l_state.mutexes[i] < 0);
        }
    }
#endif

    /* destroy mutex counts */
    l_state.num_mutexes = 0;

    /* destroy the mutexes */
    _my_free(l_state.mutexes);
    l_state.mutexes = NULL;

    return COMEX_SUCCESS;
}


int comex_lock(int mutex, int proc)
{
    header_t *header = NULL;
    char *notify = NULL;

#if DEBUG
    printf("[%d] comex_lock id=%d proc=%d\n", l_state.rank, mutex, proc);
#endif

    notify = _my_malloc(sizeof(char));
    *notify = 1;

    header = _my_malloc(sizeof(header_t));
    header->operation = OP_LOCK_REQUEST;
    header->remote_address = NULL;
    header->local_address = NULL;
    header->length = mutex;
    header->notify_address = notify;

    _mq_push(proc, (char*)header, sizeof(header_t));

    while (*notify) {
#if DEBUG
        printf("notify=%d\n", (int)*notify);
#endif
        comex_make_progress();
    }

    _my_free(notify);

    return COMEX_SUCCESS;
}


int comex_unlock(int mutex, int proc)
{
#if DEBUG
    printf("[%d] comex_unlock id=%d proc=%d\n", l_state.rank, mutex, proc);
#endif

    header_t *header = NULL;
    header = _my_malloc(sizeof(header_t));
    header->operation = OP_UNLOCK;
    header->remote_address = NULL;
    header->local_address = NULL;
    header->length = mutex;
    header->notify_address = NULL;

    _mq_push(proc, (char*)header, sizeof(header_t));
    _make_progress_if_needed();

    return COMEX_SUCCESS;
}


int comex_malloc(void **ptrs, size_t size, comex_group_t group)
{
    comex_igroup_t *igroup = NULL;
    MPI_Comm comm = MPI_COMM_NULL;
    int comm_rank = -1;
    int rc = MPI_SUCCESS; 

#if DEBUG
    printf("[%d] comex_malloc_group(ptrs=%p, size=%ld, ...)\n",
            l_state.rank, ptrs, size);
#endif

    /* preconditions */
    assert(ptrs);
    assert(group >= 0);
   
    igroup = comex_get_igroup_from_group(group);
    comm = igroup->comm;
    assert(comm != MPI_COMM_NULL);
    MPI_Comm_rank(comm, &comm_rank);

    /* allocate and register segment */
    ptrs[comm_rank] = comex_malloc_local(sizeof(char)*size);
  
    /* exchange buffer address */
    comex_barrier(group); /* end ARMCI epoch, enter MPI epoch */
    rc = MPI_Allgather(MPI_IN_PLACE, 0, MPI_LONG, ptrs, 1, MPI_LONG, comm);
    assert(MPI_SUCCESS == rc);
#if DEBUG
    {
        int i;
        int size;
        MPI_Comm_size(comm, &size);
        for (i=0; i<size; ++i) {
            printf("[%d] ptrs[%d]=%p\n", l_state.rank, i, ptrs[i]);
        }
    }
#endif

    /* TODO: This isn't needed! Right? The Allgather above is like a barrier */
    /*comex_barrier();*/

    return COMEX_SUCCESS;
}


int comex_free(void *ptr, comex_group_t group)
{
    comex_igroup_t *igroup = NULL;
    MPI_Comm comm = MPI_COMM_NULL;

    /* preconditions */
    assert(NULL != ptr);
    assert(group >= 0);

    igroup = comex_get_igroup_from_group(group);
    comm = igroup->comm;

    /* remove my ptr from reg cache and free ptr */
    comex_free_local(ptr);

    /* Synchronize: required by ARMCI semantics */
    comex_barrier(group);

    return COMEX_SUCCESS;
}

