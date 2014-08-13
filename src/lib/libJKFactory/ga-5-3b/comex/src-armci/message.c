#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>

#include "armci.h"
#include "message.h"
#include "comex.h"

#if COMEX_NETWORK_MPI_TS
extern void comex_make_progress(void);
#define COMEX_TAG 27624
#endif

/* hacky alias for MPI_COMM_SELF */
static ARMCI_Group ARMCI_GROUP_SELF = -2;

extern int ARMCI_Default_Proc_Group;

/* for armci_msg_sel_scope */
static MPI_Datatype MPI_LONGLONG_INT;

static MPI_Comm wc()
{
    MPI_Comm comm;
    assert(COMEX_SUCCESS == comex_group_comm(COMEX_GROUP_WORLD, &comm));
    return comm;
}

/* undocumented, but used in GA to expose MPI_Comm */
MPI_Comm armci_group_comm(ARMCI_Group *group)
{
    MPI_Comm comm;
    assert(COMEX_SUCCESS == comex_group_comm(*group, &comm));
    return comm;
}


static MPI_Datatype armci_type_to_mpi_type(int type)
{
    MPI_Datatype mpi_dt;

    if (type == ARMCI_INT) {
        mpi_dt = MPI_INT;
    }
    else if (type == ARMCI_LONG) {
        mpi_dt = MPI_LONG;
    }
    else if (type == ARMCI_LONG_LONG) {
        mpi_dt = MPI_LONG_LONG;
    }
    else if (type == ARMCI_FLOAT) {
        mpi_dt = MPI_FLOAT;
    }
    else if (type == ARMCI_DOUBLE) {
        mpi_dt = MPI_DOUBLE;
    }
    else {
        assert(0);
    }

    return mpi_dt;
}


static MPI_Op armci_op_to_mpi_op(char *op)
{
    if (strncmp(op, "+", 1) == 0) {
        return MPI_SUM;
    }
    else if (strncmp(op, "max", 3) == 0) {
        return MPI_MAX;
    }
    else if (strncmp(op, "min", 3) == 0) {
        return MPI_MIN;
    }
    else if (strncmp(op, "*", 1) == 0) {
        return MPI_PROD;
    }
    else if (strncmp(op, "absmin", 6) == 0) {
        return MPI_MIN;
    }
    else if (strncmp(op, "absmax", 6) == 0) {
        return MPI_MAX;
    }
    else if (strncmp(op, "or", 2) == 0) {
        return MPI_BOR;
    }
    else {
        printf("Unsupported gop operation:%s\n",op);
        assert(0);
    }
}


static void do_abs(void *x, int n, int type)
{
#define ARMCI_ABS_INT(a)  (((a) >= 0)   ? (a) : (-(a)))
#define ARMCI_ABS_FLT(a)  (((a) >= 0.0) ? (a) : (-(a)))
#define DO_ABS(ARMCI_TYPE, C_TYPE, WHICH)       \
    if (type == ARMCI_TYPE) {                   \
        int i;                                  \
        C_TYPE *y = (C_TYPE *)x;                \
        for (i = 0; i < n; i++) {               \
            y[i] = ARMCI_ABS_##WHICH(y[i]);     \
        }                                       \
    }                                           \
    else
    DO_ABS(ARMCI_INT,       int,        INT)
    DO_ABS(ARMCI_LONG,      long,       INT)
    DO_ABS(ARMCI_LONG_LONG, long long,  INT)
    DO_ABS(ARMCI_FLOAT,     float,      FLT)
    DO_ABS(ARMCI_DOUBLE,    double,     FLT)
    {
        assert(0);
    }
#undef ARMCI_ABS_INT
#undef ARMCI_ABS_FLT
#undef DO_ABS
}


static MPI_Comm get_comm(ARMCI_Group *group)
{
    MPI_Comm comm;
    assert(COMEX_SUCCESS == comex_group_comm(*group, &comm));
    return comm;
}


static ARMCI_Group get_default_group()
{
    ARMCI_Group group;
    ARMCI_Group_get_default(&group);
    return group;
}


static MPI_Comm get_default_comm()
{
    ARMCI_Group group;
    ARMCI_Group_get_default(&group);
    return get_comm(&group);
}


static int get_default_rank()
{
    int rank;
    MPI_Comm_rank(get_default_comm(), &rank);
    return rank;
}


static void do_gop(void *x, int n, char* op, int type, ARMCI_Group group)
{
    MPI_Comm comm = MPI_COMM_NULL;
    int mpi_type_size = 0;
    MPI_Datatype mpi_type = MPI_DATATYPE_NULL;
    int rc = 0;
    void *result = NULL;
    MPI_Op mpi_op = MPI_OP_NULL;

    if (ARMCI_GROUP_SELF == group) {
        comm = MPI_COMM_SELF;
    }
    else {
        comm = get_comm(&group);
    }
    mpi_type = armci_type_to_mpi_type(type);
    MPI_Type_size(mpi_type, &mpi_type_size);
    mpi_op = armci_op_to_mpi_op(op);

    if (strncmp(op, "absmin", 6) == 0) {
        do_abs(x, n, type);
    }
    else if (strncmp(op, "absmax", 6) == 0) {
        do_abs(x, n, type);
    }
        
    result = malloc(n*mpi_type_size);
    assert(result);

    if (ARMCI_GROUP_SELF != group) {
        comex_barrier(group);
    }
    rc = MPI_Allreduce(x, result, n, mpi_type, mpi_op, comm); 
    assert(rc == MPI_SUCCESS);

    memcpy(x, result, mpi_type_size * n);
    free(result);
}


void armci_msg_bcast(void *buf, int len, int root)
{
    assert(buf != NULL);
    comex_barrier(ARMCI_Default_Proc_Group);
    MPI_Bcast(buf, len, MPI_BYTE, root, get_default_comm());
}


/* the payload is a struct with a union e.g.
 *
 * typedef struct {
 *    union val_t {double dval; int ival; long lval; long long llval; float fval;}v;
 *    Integer subscr[MAXDIM];
 *    DoubleComplex extra;
 *    SingleComplex extra2;
 * } elem_info_t;
 *
 * The key piece is the first sizeof(double) bytes. The rest of the struct
 * simply tags along for the communication and can be represented as a byte
 * stream.
 *
 * The 'n' parameter is the size of the entire payload i.e.
 * sizeof(struct elem_info_t).
 *
 * We really care which process has the min/max value and then bcast the
 * entire payload using the min/max answer as the root.
 */
void armci_msg_sel_scope(int scope, void *x, int n, char* op, int type, int contribute)
{
    static int initialized = 0;
    MPI_Op mpi_op = MPI_OP_NULL;
    MPI_Comm comm = get_default_comm();

    if (SCOPE_NODE == scope) {
        comm = MPI_COMM_SELF;
    }

    /* first time this function is called we establish the
     * long long w/ int type that MPI doesn't provide by default */
    if (!initialized) {
        int block[2];
        MPI_Aint disp[2];
        MPI_Datatype type[2];

        initialized = 1;
        type[0] = MPI_LONG_LONG;
        type[1] = MPI_INT;
        disp[0] = 0;
        disp[1] = sizeof(long long);
        block[0] = 1;
        block[1] = 1;
        MPI_Type_struct(2, block, disp, type, &MPI_LONGLONG_INT);
    }

    if (strncmp(op, "min", 3) == 0) {
        mpi_op = MPI_MINLOC;
    }
    else if (strncmp(op, "max", 3) == 0) {
        mpi_op = MPI_MAXLOC;
    }
    else {
        assert(0);
    }

#define SELECT(ARMCI_TYPE, C_TYPE, MPI_TYPE)                    \
    if (type == ARMCI_TYPE) {                                   \
        struct {                                                \
            C_TYPE val;                                         \
            int rank;                                           \
        } in, out;                                              \
        in.val = *((C_TYPE*)x);                                 \
        in.rank = get_default_rank();                           \
        if (SCOPE_NODE != scope) {                              \
            comex_barrier(ARMCI_Default_Proc_Group);            \
        }                                                       \
        MPI_Allreduce(&in, &out, 1, MPI_TYPE, mpi_op, comm);    \
        armci_msg_bcast(x, n, out.rank);                        \
    }                                                           \
    else
    SELECT(ARMCI_INT,       int,        MPI_2INT)
    SELECT(ARMCI_LONG,      long,       MPI_LONG_INT)
    SELECT(ARMCI_LONG_LONG, long long,  MPI_LONGLONG_INT)
    SELECT(ARMCI_FLOAT,     float,      MPI_FLOAT_INT)
    SELECT(ARMCI_DOUBLE,    double,     MPI_DOUBLE_INT)
    {
        assert(0);
    }
#undef SELECT
}


void armci_msg_bcast_scope(int scope, void* buffer, int len, int root)
{
    if (SCOPE_ALL == scope || SCOPE_MASTERS == scope) {
        armci_msg_bcast(buffer, len, root);
    }
    else if (SCOPE_NODE == scope) {
        assert(buffer != NULL);
        MPI_Bcast(buffer, len, MPI_BYTE, root, MPI_COMM_SELF);
    }
    else {
        assert(0);
    }
}


void armci_msg_brdcst(void* buffer, int len, int root)
{
    armci_msg_bcast(buffer, len, root);
}


/* there was a case in ghost update where a proc sent a message to itself */
static MPI_Request self_request = MPI_REQUEST_NULL;
static int self_request_flag = 0;

void armci_msg_snd(int tag, void* buffer, int len, int to)
{
    MPI_Comm comm = wc();
    MPI_Request request;
    MPI_Status status;
    int self;
    int flag = 0;
    int rc;

    rc = MPI_Comm_rank(comm, &self);
    assert(MPI_SUCCESS == rc);
#if COMEX_NETWORK_MPI_TS
    assert(COMEX_TAG != tag);
#endif
    rc = MPI_Isend(buffer, len, MPI_CHAR, to, tag, wc(), &request);
    assert(MPI_SUCCESS == rc);
    if (to == self) {
        /* make sure this proc hasn't already sent to itself -- we only allow
         * one outstanding self send */
        assert(!self_request_flag);
        self_request_flag = 1;
        self_request = request;
    }
    else {
        do {
            MPI_Test(&request, &flag, &status);
            assert(MPI_SUCCESS == rc);
#if COMEX_NETWORK_MPI_TS
            comex_make_progress();
#endif
        } while (!flag);
    }
}


void armci_msg_rcv(int tag, void* buffer, int len, int *msglen, int from)
{
    MPI_Comm comm = wc();
    MPI_Request request;
    MPI_Status status;
    int self;
    int flag = 0;
    int rc;

    rc = MPI_Comm_rank(comm, &self);
    assert(MPI_SUCCESS == rc);
#if COMEX_NETWORK_MPI_TS
    assert(COMEX_TAG != tag);
#endif
    rc = MPI_Irecv(buffer, len, MPI_CHAR, from, tag, wc(), &request);
    assert(MPI_SUCCESS == rc);
    if (from == self && self_request_flag) {
        do {
            MPI_Test(&self_request, &flag, &status);
            assert(MPI_SUCCESS == rc);
#if COMEX_NETWORK_MPI_TS
            comex_make_progress();
#endif
        } while (!flag);
        self_request = MPI_REQUEST_NULL;
        self_request_flag = 0;
    }
    do {
        MPI_Test(&request, &flag, &status);
        assert(MPI_SUCCESS == rc);
#if COMEX_NETWORK_MPI_TS
        comex_make_progress();
#endif
    } while (!flag);

    if(msglen) {
        rc = MPI_Get_count(&status, MPI_CHAR, msglen);
        assert(MPI_SUCCESS == rc);
    }
}


int armci_msg_rcvany(int tag, void* buffer, int len, int *msglen)
{
    MPI_Comm comm = wc();
    MPI_Request request;
    MPI_Status status;
    int rank;
    int flag = 0;
    int rc;

    rc = MPI_Comm_rank(comm, &rank);
    assert(MPI_SUCCESS == rc);
#if COMEX_NETWORK_MPI_TS
    assert(COMEX_TAG != tag);
#endif
    rc = MPI_Irecv(buffer, len, MPI_CHAR, MPI_ANY_SOURCE, tag, wc(), &request);
    assert(MPI_SUCCESS == rc);
    do {
        MPI_Test(&request, &flag, &status);
        assert(MPI_SUCCESS == rc);
#if COMEX_NETWORK_MPI_TS
        comex_make_progress();
#endif
    } while (!flag);

    if(msglen) {
        rc = MPI_Get_count(&status, MPI_CHAR, msglen);
        assert(MPI_SUCCESS == rc);
    }

    return (int)status.MPI_SOURCE;
}


void armci_msg_reduce(void *x, int n, char *op, int type)
{
    do_gop(x, n, op, type, get_default_group());
}


void armci_msg_reduce_scope(int scope, void *x, int n, char *op, int type)
{
    if (SCOPE_NODE == scope) {
        do_gop(x, n, op, type, ARMCI_GROUP_SELF);
    } else {
        do_gop(x, n, op, type, get_default_group());
    }
}


void armci_msg_gop_scope(int scope, void *x, int n, char* op, int type)
{
    if (SCOPE_NODE == scope) {
        do_gop(x, n, op, type, ARMCI_GROUP_SELF);
    } else {
        do_gop(x, n, op, type, get_default_group());
    }
}


void armci_msg_igop(int *x, int n, char* op)
{
    do_gop(x, n, op, ARMCI_INT, get_default_group());
}


void armci_msg_lgop(long *x, int n, char* op)
{
    do_gop(x, n, op, ARMCI_LONG, get_default_group());
}


void armci_msg_llgop(long long *x, int n, char* op)
{
    do_gop(x, n, op, ARMCI_LONG_LONG, get_default_group());
}


void armci_msg_fgop(float *x, int n, char* op)
{
    do_gop(x, n, op, ARMCI_FLOAT, get_default_group());
}


void armci_msg_dgop(double *x, int n, char* op)
{
    do_gop(x, n, op, ARMCI_DOUBLE, get_default_group());
}


void armci_exchange_address(void *ptr_ar[], int n)
{
    ARMCI_Group group;
    ARMCI_Group_get_default(&group);
    armci_exchange_address_grp(ptr_ar, n, &group);
}


void parmci_msg_barrier()
{
    comex_barrier(ARMCI_Default_Proc_Group);
    MPI_Barrier(get_default_comm());
}


void armci_msg_bintree(int scope, int* Root, int *Up, int *Left, int *Right)
{
    int root, up, left, right, index, nproc;

    assert(SCOPE_NODE != scope);
    assert(SCOPE_MASTERS != scope);

    root  = 0;
    nproc = armci_msg_nproc();
    index = armci_msg_me() - root;
    up    = (index-1)/2 + root; if( up < root) up = -1;
    left  = 2*index + 1 + root; if(left >= root+nproc) left = -1;
    right = 2*index + 2 + root; if(right >= root+nproc)right = -1;

    *Up = up;
    *Left = left;
    *Right = right;
    *Root = root;
}


int armci_msg_me()
{
    int rank;
    assert(comex_initialized());
    assert(wc() != MPI_COMM_NULL);
    MPI_Comm_rank(wc(), &rank);
    return rank;
}


int armci_msg_nproc()
{
    int size;
    assert(comex_initialized());
    assert(wc() != MPI_COMM_NULL);
    MPI_Comm_size(wc(), &size);
    return size;
}


void armci_msg_abort(int code)
{
    fprintf(stderr, "Exiting, Error in Communication\n");
    MPI_Abort(wc(), code);
}


void armci_msg_init(int *argc, char ***argv)
{
    int flag;
    MPI_Initialized(&flag);
    if(!flag) {
        MPI_Init(argc, argv);
    }
}


void armci_msg_finalize()
{
    int flag;
    MPI_Initialized(&flag);
    assert(flag);
    MPI_Finalize();
}


double armci_timer()
{
    return MPI_Wtime();
}


void armci_msg_clus_brdcst(void *buf, int len)
{
    assert(0);
}


void armci_msg_clus_igop(int *x, int n, char* op)
{
    assert(0);
}


void armci_msg_clus_fgop(float *x, int n, char* op)
{
    assert(0);
}


void armci_msg_clus_lgop(long *x, int n, char* op)
{
    assert(0);
}


void armci_msg_clus_llgop(long long *x, int n, char* op)
{
    assert(0);
}


void armci_msg_clus_dgop(double *x, int n, char* op)
{
    assert(0);
}


void armci_msg_group_gop_scope(int scope, void *x, int n, char* op, int type, ARMCI_Group *group)
{
    if (SCOPE_NODE == scope) {
        do_gop(x, n, op, type, ARMCI_GROUP_SELF);
    } else {
        do_gop(x, n, op, type, *group);
    }
}


void armci_msg_group_igop(int *x, int n, char* op, ARMCI_Group *group)
{
    do_gop(x, n, op, ARMCI_INT, *group);
}


void armci_msg_group_lgop(long *x, int n, char* op, ARMCI_Group *group)
{
    do_gop(x, n, op, ARMCI_LONG, *group);
}


void armci_msg_group_llgop(long long *x, int n, char* op, ARMCI_Group *group)
{
    do_gop(x, n, op, ARMCI_LONG_LONG, *group);
}


void armci_msg_group_fgop(float *x, int n, char* op, ARMCI_Group *group)
{
    do_gop(x, n, op, ARMCI_FLOAT, *group);
}


void armci_msg_group_dgop(double *x, int n,char* op, ARMCI_Group *group)
{
    do_gop(x, n, op, ARMCI_DOUBLE, *group);
}


void armci_exchange_address_grp(void *ptr_arr[], int n, ARMCI_Group *group)
{
    assert(0);
#if 0
    MPI_Datatype mpi_datatype;

    if (sizeof(void*) == sizeof(int)) {
        mpi_datatype = MPI_INT;
    }
    else if (sizeof(void*) == sizeof(long)) {
        mpi_datatype = MPI_LONG;
    }
    else if (sizeof(void*) == sizeof(long long)) {
        mpi_datatype = MPI_LONG_LONG;
    }
    else {
        assert(0);
    }

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
            ptr_ar, n, mpi_datatype, get_comm(group));
#endif
}


void parmci_msg_group_barrier(ARMCI_Group *group)
{
    comex_barrier(*group);
    MPI_Barrier(get_comm(group));
}


void armci_msg_group_bcast_scope(int scope, void *buf, int len, int root, ARMCI_Group *group)
{
    if (SCOPE_NODE == scope) {
        MPI_Bcast(buf, len, MPI_BYTE, root, MPI_COMM_SELF);
    }
    else {
        int root_sub;
        int err;
        /* NOTE: this function is passed a root which has been translated back
         * to the world group while the group passed in is the sub
         * communicator group... what a mess */
        MPI_Group group_world;
        MPI_Group group_sub;
        
        err = MPI_Comm_group(wc(), &group_world);
        assert(MPI_SUCCESS == err);
        err = MPI_Comm_group(get_comm(group), &group_sub);
        assert(MPI_SUCCESS == err);
        err = MPI_Group_translate_ranks(group_world, 1, &root,
                group_sub, &root_sub);
        comex_barrier(*group);
        err = MPI_Bcast(buf, len, MPI_BYTE, root_sub, get_comm(group));
        assert(MPI_SUCCESS == err);
    }
}


void armci_grp_clus_brdcst(void *buf, int len, int grp_master, int grp_clus_nproc,ARMCI_Group *mastergroup)
{
    assert(0);
}
