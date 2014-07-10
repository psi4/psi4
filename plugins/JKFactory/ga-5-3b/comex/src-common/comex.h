/* comex header file */
#ifndef _COMEX_H
#define _COMEX_H   

#include <mpi.h>

#include <stdlib.h>

#if defined(__cplusplus) || defined(c_plusplus)
extern "c" {
#endif

typedef struct {
    void **src; /**< array of source starting addresses */
    void **dst; /**< array of destination starting addresses */
    int count; /**< size of address arrays (src[count],dst[count]) */
    int bytes; /**< length in bytes for each src[i]/dst[i] pair */
} comex_giov_t;

typedef int comex_request_t;

typedef int comex_group_t;

#define COMEX_GROUP_WORLD 0
#define COMEX_GROUP_NULL -1
 
#define COMEX_SUCCESS 0
#define COMEX_FAILURE 1

#define COMEX_SWAP 10
#define COMEX_SWAP_LONG 11
#define COMEX_FETCH_AND_ADD 12
#define COMEX_FETCH_AND_ADD_LONG 13

#define COMEX_ACC_OFF 36
#define COMEX_ACC_INT (COMEX_ACC_OFF + 1)
#define COMEX_ACC_DBL (COMEX_ACC_OFF + 2)
#define COMEX_ACC_FLT (COMEX_ACC_OFF + 3)
#define COMEX_ACC_CPL (COMEX_ACC_OFF + 4)
#define COMEX_ACC_DCP (COMEX_ACC_OFF + 5)
#define COMEX_ACC_LNG (COMEX_ACC_OFF + 6)

#define COMEX_MAX_STRIDE_LEVEL 8

/**
 * Initialize comex.
 *
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_init();

/**
 * Initialize comex with command line arguments.
 *
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_init_args(int *argc, char ***argv);

/**
 * Test whether comex has been initialized.
 *
 * @return COMEX_SUCCESS if comex has been initialized
 *         COMEX_FAILURE if comex has not
 */
extern int comex_initialized();

/**
 * Terminate comex and clean up resources.
 *
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_finalize();

/**
 * Abort comex, printing the msg, and exiting with code.
 *
 * @param[in] msg the message to print
 * @param[in] code the code to exit with
 */
extern void comex_error(char *msg, int code);

/**
 * Create a new group from the given group and process ID list.
 *
 * The rank list selects the ranks from the given group to become members of
 * the new group. The ranks should be nonnegative and range from zero to the
 * size of the given group.
 *
 * This functions is collective only over the ranks within the rank list and
 * not over the entire original group.
 *
 * @param[in] n the number of ranks to select for the new group
 * @param[in] rank_list the list of ranks to select for the new group
 * @param[in] group the group to subset for the new group
 * @param[out] new_group the newly created group
 * @return COMEX_SUCCESS on success
 *         COMEX_FAILURE if a rank in the rank list is out of bounds
 */
extern int comex_group_create(
        int n, int *pid_list, comex_group_t group, comex_group_t *new_group);

/**
 * Marks the group for deallocation.
 *
 * @param[in] group group to be destroyed
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_group_free(comex_group_t group);

/**
 * Determines the rank of the calling process in the given group.
 *
 * @param[in] group group handle
 * @param[out] rank rank of the calling process in the group
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_group_rank(comex_group_t group, int *rank);

/**
 * Determines the size of the given group.
 *
 * @param[in] group group handle
 * @param[out] size number of processes in the group
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_group_size(comex_group_t group, int *size);

/**
 * Returns the MPI_Comm object backing the given group.
 *
 * The actual MPI_Comm object is returned, therefore do not call
 * MPI_Comm_free() on the returned communicator. This function is for
 * convenience to be able to MPI_Comm_dup() the returned MPI_Comm instance.
 *
 * @param[in] group group handle
 * @param[out] comm the communicator handle
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_group_comm(comex_group_t group, MPI_Comm *comm);

/**
 * Translates the ranks of processes in one group to those in another group.
 *
 * @param[in] n the number of ranks in the ranks_from and ranks_to arrays
 * @param[in] group_from the group to translate ranks from 
 * @param[in] ranks_from array of zer or more valid ranks in group_from
 * @param[in] group_to the group to translate ranks to 
 * @param[out] ranks_to array of corresponding ranks in group_to
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_group_translate_ranks(int n,
        comex_group_t group_from, int *ranks_from,
        comex_group_t group_to, int *ranks_to);

/**
 * Translate the given rank from its group to its corresponding rank in the
 * world group.
 *
 * Shorthand notation for common case.
 *
 * @param[in] group the group to translate from
 * @param[in] group_rank the rank to translate from
 * @param[out] world_rank the corresponding world rank
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_group_translate_world(
        comex_group_t group, int group_rank, int *world_rank);

/**
 * A collective communication and operations barrier.
 *
 * Ensures all comex communication has completed prior to performing the
 * operations barrier.
 *
 * @param[in] group the group to perform the collective barrier over
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_barrier(comex_group_t group);

/**
 * Contiguous Put.
 *
 * @param[in] src pointer to 1st segment at source
 * @param[in] dst pointer to 1st segment at destination
 * @param[in] bytes number of bytes to transfer
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_put(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group);

/**
 * Strided Put.
 *
 * @param[in] src pointer to 1st segment at source
 * @param[in] src_stride array of strides at source
 * @param[in] dst pointer to 1st segment at destination
 * @param[in] dst_stride array of strides at destination
 * @param[in] count number of units at each stride level count[0]=bytes
 * @param[in] stride_levels number of stride levels
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @return COMEX_SUCCESS on success
 */
extern int comex_puts(
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels,
        int proc, comex_group_t group);

/**
 * Vector Put.
 *
 * @param[in] darr descriptor array
 * @param[in] len length of descriptor array
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @return COMEX_SUCCESS on success
 */
extern int comex_putv(
        comex_giov_t *darr, int len,
        int proc, comex_group_t group);

/**
 * Nonblocking Contiguous Put.
 *
 * @param[in] src pointer to 1st segment at source
 * @param[in] dst pointer to 1st segment at destination
 * @param[in] bytes number of bytes to transfer
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @param[out] nb_handle nonblocking request object
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_nbput(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group,
        comex_request_t* nb_handle);

/**
 * Nonblocking Strided Put.
 *
 * @param[in] src pointer to 1st segment at source
 * @param[in] src_stride array of strides at source
 * @param[in] dst pointer to 1st segment at destination
 * @param[in] dst_stride array of strides at destination
 * @param[in] count number of units at each stride level count[0]=bytes
 * @param[in] stride_levels number of stride levels
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @param[out] nb_handle nonblocking request object
 * @return COMEX_SUCCESS on success
 */
extern int comex_nbputs(
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels,
        int proc, comex_group_t group,
        comex_request_t* nb_handle);

/**
 * Nonblocking Vector Put.
 *
 * @param[in] darr descriptor array
 * @param[in] len length of descriptor array
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @param[out] nb_handle nonblocking request object
 * @return COMEX_SUCCESS on success
 */
extern int comex_nbputv(
        comex_giov_t *darr, int len,
        int proc, comex_group_t group,
        comex_request_t* nb_handle);

/**
 * Contiguous Atomic Accumulate.
 *
 * @param[in] op operation
 * @param[in] scale factor x += scale*y
 * @param[in] src pointer to 1st segment at source
 * @param[in] dst pointer to 1st segment at destination
 * @param[in] bytes number of bytes to transfer
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @return COMEX_SUCCESS on success
 */
extern int comex_acc(
        int op, void *scale,
        void *src, void *dst, int bytes,
        int proc, comex_group_t group);

/**
 * Strided Atomic Accumulate.
 *
 * @param[in] op operation
 * @param[in] scale factor x += scale*y
 * @param[in] src pointer to 1st segment at source
 * @param[in] src_stride [stride_levels] array of strides at source
 * @param[in] dst pointer to 1st segment at destination
 * @param[in] dst_stride [stride_levels] array of strides at destination
 * @param[in] count [stride_levels+1] number of units at each stride level
 *            count[0]=bytes
 * @param[in] stride_levels number of stride levels
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @return COMEX_SUCCESS on success
 */
extern int comex_accs(
        int op, void *scale,
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels,
        int proc, comex_group_t group);

/**
 * Vector Atomic Accumulate.
 *
 * @param[in] op operation
 * @param[in] scale factor x += scale*y
 * @param[in] darr descriptor array
 * @param[in] len length of descriptor array
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @return COMEX_SUCCESS on success
 */
extern int comex_accv(
        int op, void *scale,
        comex_giov_t *darr, int len,
        int proc, comex_group_t group);

/**
 * Nonblocking Contiguous Atomic Accumulate.
 *
 * @param[in] op operation
 * @param[in] scale factor x += scale*y
 * @param[in] src pointer to 1st segment at source
 * @param[in] dst pointer to 1st segment at destination
 * @param[in] bytes number of bytes to transfer
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @param[out] nb_handle nonblocking request object
 * @return COMEX_SUCCESS on success
 */
extern int comex_nbacc(
        int op, void *scale,
        void *src, void *dst, int bytes,
        int proc, comex_group_t group,
        comex_request_t *nb_handle);

/**
 * Strided Atomic Accumulate.
 *
 * @param[in] op operation
 * @param[in] scale factor x += scale*y
 * @param[in] src pointer to 1st segment at source
 * @param[in] src_stride array of strides at source
 * @param[in] dst pointer to 1st segment at destination
 * @param[in] dst_stride array of strides at destination
 * @param[in] count number of units at each stride level count[0]=bytes
 * @param[in] stride_levels number of stride levels
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @param[out] nb_handle nonblocking request object
 * @return COMEX_SUCCESS on success
 */
extern int comex_nbaccs(
        int  op, void *scale,
        void *src, int *src_stride,
		void *dst, int *dst_stride,
		int *count, int stride_levels,
        int proc, comex_group_t group,
        comex_request_t* nb_handle);

/**
 * Vector Atomic Accumulate.
 *
 * @param[in] op operation
 * @param[in] scale factor x += scale*y
 * @param[in] darr descriptor array
 * @param[in] len length of descriptor array
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @param[out] nb_handle nonblocking request object
 * @return COMEX_SUCCESS on success
 */
extern int comex_nbaccv(
        int op, void *scale,
        comex_giov_t *darr, int len,
        int proc, comex_group_t group,
        comex_request_t* nb_handle);

/**
 * Contiguous Get.
 *
 * @param[in] src pointer to 1st segment at source
 * @param[in] dst pointer to 1st segment at destination
 * @param[in] bytes number of bytes to transfer
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_get(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group);

/**
 * Strided Get.
 *
 * @param[in] src pointer to 1st segment at source
 * @param[in] src_stride array of strides at source
 * @param[in] dst pointer to 1st segment at destination
 * @param[in] dst_stride array of strides at destination
 * @param[in] count number of units at each stride level count[0]=bytes
 * @param[in] stride_levels number of stride levels
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @return COMEX_SUCCESS on success
 */
extern int comex_gets(
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels,
        int proc, comex_group_t group);

/**
 * Vector Get.
 *
 * @param[in] darr descriptor array
 * @param[in] len length of descriptor array
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @return COMEX_SUCCESS on success
 */
extern int comex_getv(
        comex_giov_t *darr, int len,
        int proc, comex_group_t group);

/**
 * Nonblocking Contiguous Get.
 *
 * @param[in] src pointer to 1st segment at source
 * @param[in] dst pointer to 1st segment at destination
 * @param[in] bytes number of bytes to transfer
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @param[out] nb_handle nonblocking request object
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_nbget(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group,
        comex_request_t* nb_handle);

/**
 * Nonblocking Strided Get.
 *
 * @param[in] src pointer to 1st segment at source
 * @param[in] src_stride array of strides at source
 * @param[in] dst pointer to 1st segment at destination
 * @param[in] dst_stride array of strides at destination
 * @param[in] count number of units at each stride level count[0]=bytes
 * @param[in] stride_levels number of stride levels
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @param[out] nb_handle nonblocking request object
 * @return COMEX_SUCCESS on success
 */
extern int comex_nbgets(
        void *src, int *src_stride,
		void *dst, int *dst_stride,
		int *count, int stride_levels,
        int proc, comex_group_t group,
        comex_request_t *nb_handler);

/**
 * Nonblocking Vector Get.
 *
 * @param[in] darr descriptor array
 * @param[in] len length of descriptor array
 * @param[in] proc remote process(or) id
 * @param[in] group the calling process and remote process must belong to the
 *            same group
 * @param[out] nb_handle nonblocking request object
 * @return COMEX_SUCCESS on success
 */
extern int comex_nbgetv(
        comex_giov_t *darr, int len,
        int proc, comex_group_t group,
        comex_request_t* nb_handle);

/**
 * Collective allocation of registered memory and exchange of addresses.
 *
 * @param[out] ptr_arr array of memory addresses
 *             w.r.t. each process's address space
 * @param[in] bytes how many bytes to allocate locally
 * @param[in] group the group to which the calling process belongs
 * @return COMEX_SUCCESS on success
 */
extern int comex_malloc(
        void **ptr_arr, size_t bytes, comex_group_t group);

/**
 * Collective free of memory given the original local pointer.
 *
 * @param[in] ptr the original local memory allocated using comex_malloc
 * @param[in] group the group to which the calling process belongs
 * @return COMEX_SUCCESS on success
 */
extern int comex_free(void *ptr, comex_group_t group);

/**
 * Local (noncollective) allocation of registered memory.
 *
 * Using memory allocated here may have performance benefits when used as a
 * communication buffer.
 *
 * @param[in] bytes how many bytes to allocate locally
 * @return COMEX_SUCCESS on success
 */
extern void* comex_malloc_local(size_t bytes);

/**
 * Local (noncollective) free of memory allocated by comex_malloc_local.
 *
 * @param[in] the original local memory allocated using comex_malloc_local
 * @return COMEX_SUCCESS on success
 */
extern int comex_free_local(void *ptr);

/**
 * Flush all outgoing messages from me to the given proc.
 *
 * @param[in] proc the proc with which to flush outgoing messages
 * @return COMEX_SUCCESS on success
 */
extern int comex_fence_proc(int proc, comex_group_t group);

/**
 * Flush all outgoing messages to all procs.
 *
 * @return COMEX_SUCCESS on success
 */
extern int comex_fence_all(comex_group_t group);

/**
 * Collectively create num locks locally.
 *
 * Remote procs may create a different number of locks, including zero.
 *
 * This function is always collective on the world group.
 *
 * @param[in] num number of locks to create locally
 * @return COMEX_SUCCESS on success
 */
extern int comex_create_mutexes(int num);

/**
 * Collectively destroy all previously created locks.
 *
 * This function is always collective on the world group.
 *
 * @param[in] num number of locks to create locally
 * @return COMEX_SUCCESS on success
 */
extern int comex_destroy_mutexes();

/**
 * Lock the given mutex on the given proc.
 *
 * This function is always on the world group.
 *
 * @param[in] mutex the ID of the mutex to lock on proc
 * @param[in] the ID of the proc which owns the mutex
 *
 * @return COMEX_SUCCESS on success
 *         COMEX_FAILURE if given mutex or proc is out of range
 */
extern int comex_lock(int mutex, int proc);

/**
 * Unlock the given mutex on the given proc.
 *
 * This function is always on the world group.
 *
 * @param[in] mutex the ID of the mutex to unlock on proc
 * @param[in] the ID of the proc which owns the mutex
 *
 * @return COMEX_SUCCESS on success
 *         COMEX_FAILURE if given mutex or proc is out of range
 */
extern int comex_unlock(int mutex, int proc);

/**
 * Read-modify-write atomic operation.
 *
 * The operations may be one of
 *  - COMEX_SWAP
 *  - COMEX_SWAP_LONG
 *  - COMEX_FETCH_AND_ADD
 *  - COMEX_FETCH_AND_ADD_LONG
 *
 * For the swap operations, the extra parameter is not used. The values of the
 * ploc and prem locations are swapped.
 *
 * For the fetch and add operations, the extra parameter is also used to
 * indicate how much to increment the remote value. The original remove value
 * is returned in the ploc parameter.
 *
 * @param[in] op the operation to perform (see list above)
 * @param[in] ploc the value to update locally
 * @param[in] prem the value to update remotely
 * @param[in] extra for COMEX_FETCH_AND_ADD and COMEX_FETCH_AND_ADD_LONG, the
 *            amount to increment the remote value by
 * @param[in] proc remote process(or) id
 * @param[in] group group handle
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_rmw(
        int op, void *ploc, void *prem, int extra,
        int proc, comex_group_t group);

/**
 * Waits for completion of non-blocking comex operations with explicit handles.
 *
 * @param[in] nb_handle the handle
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_wait(comex_request_t *nb_handle);

/**
 * Checks completion status of non-blocking comex operations with explicit
 * handles.
 *
 * @param[in] nb_handle the handle
 * @param[out] status 0-completed, 1-in progress
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_test(comex_request_t *nb_handle, int *status);

/**
 * Wait for all outstanding implicit non-blocking operations to finish.
 *
 * @param[in] group group handle
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_wait_all(comex_group_t group);

/**
 * Wait for all outstanding implicit non-blocking operations to a particular
 * process to finish.
 *
 * @param[in] proc proc for which all the outstanding non-blocking operations
 * have to be completed
 * @param[in] group group handle
 * @return COMEX_SUCCESS on sucess
 */
extern int comex_wait_proc(int proc, comex_group_t group);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif /* _COMEX_H */
