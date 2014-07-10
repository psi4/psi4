#ifndef COMEX_IMPL_H_
#define COMEX_IMPL_H_

#include <mpi.h>

//#define MAX_NB_OUTSTANDING 1024
#define MAX_NB_OUTSTANDING 8
#define COMEX_TAG 27624

typedef enum {
    OP_PUT = 0,
    OP_GET_REQUEST,
    OP_GET_RESPONSE,
    OP_ACC_INT,
    OP_ACC_DBL,
    OP_ACC_FLT,
    OP_ACC_CPL,
    OP_ACC_DCP,
    OP_ACC_LNG,
    OP_FENCE_REQUEST,
    OP_FENCE_RESPONSE,
    OP_BARRIER_REQUEST,
    OP_BARRIER_RESPONSE,
    OP_FETCH_AND_ADD_REQUEST,
    OP_FETCH_AND_ADD_RESPONSE,
    OP_SWAP_REQUEST,
    OP_SWAP_RESPONSE,
    OP_LOCK_REQUEST,
    OP_LOCK_RESPONSE,
    OP_UNLOCK,
} op_t;

typedef struct {
    op_t operation;
    void *remote_address;
    void *local_address;
    int length; /**< length of message/payload not including header */
    void *notify_address;
} header_t;

typedef struct message_link {
    struct message_link *next;
    int dest;
    char *message;
    MPI_Request request;
} message_t;

typedef struct get_link {
    struct get_link *next;
    char *notify_address;
} get_t;

typedef struct barrier_link {
    struct barrier_link *next;
    int world_rank;
    void *notify_address;
} barrier_t;

typedef struct lock_link {
    struct lock_link *next;
    int rank;
    int id;
    void *notify_address;
} lock_t;

typedef struct {
    MPI_Comm world_comm;
    int rank;
    int size;

    /* buffers for locks */
    int *mutexes; /**< local mutexes */
    int num_mutexes; /**< how many mutexes on this process */

    /* a queue for outgoing messages */
    int mq_size;
    message_t *mq_head;
    message_t *mq_tail;

    /* a queue for get notifications */
    get_t *gq_head;
    get_t *gq_tail;

    /* a queue for barrier requests */
    barrier_t *bq_head;
    barrier_t *bq_tail;

    /* a queue for lock requests */
    lock_t *lq_head;
    lock_t *lq_tail;

} local_state;

extern local_state l_state;

extern void comex_make_progress(void);

#endif /* COMEX_IMPL_H_ */
