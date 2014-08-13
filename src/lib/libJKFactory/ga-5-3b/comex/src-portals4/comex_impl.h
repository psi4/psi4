#ifndef COMEX_IMPL_H_
#define COMEX_IMPL_H_

#include <mpi.h>
#include <portals4.h>

#define COMEX_MAX_NB_OUTSTANDING 16

#define DEBUG 0
#define DEBUG_TO_FILE 0
#if DEBUG_TO_FILE
#   define printf(...) fprintf(l_state.my_file, __VA_ARGS__); fflush(l_state.my_file)
#else
#   define printf(...) fprintf(stderr, __VA_ARGS__); fflush(stderr)
#endif

typedef struct {
    MPI_Comm world_comm;
    int rank;
    int size;
    ptl_ni_limits_t ptl_ni_limits;
    ptl_handle_ni_t ptl_ni_handle;
    ptl_process_t ptl_process_id;
    ptl_process_t *ptl_process_ids;
    ptl_uid_t ptl_uid;
    ptl_pt_index_t ptl_pt_index;
    ptl_pt_index_t *ptl_pt_indexes;
    ptl_handle_eq_t ptl_eq_handle;
    ptl_handle_md_t ptl_md_handle;
    ptl_handle_le_t ptl_le_handle;
    long **mutexes; /**< all mutexes */
    long *local_mutex; /**< store the remote mutex value */
    unsigned int  *num_mutexes; /**< how many mutexes on each process */

#if DEBUG_TO_FILE
    FILE *my_file;
#endif
} local_state;

extern local_state l_state;

#endif /* COMEX_IMPL_H_ */
