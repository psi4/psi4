/**
 * Private header file for COMEX groups backed by MPI_comm.
 *
 * The rest of the COMEX_Group functions are defined in the public comex.h.
 *
 * @author Jeff Daily
 */
#ifndef _COMEX_GROUPS_H_
#define _COMEX_GROUPS_H_

#include <mpi.h>

/* dup of MPI_COMM_WORLD for internal MPI communication */
extern MPI_Comm COMEX_COMM_WORLD;

typedef struct group_link {
    struct group_link *next;
    comex_group_t id;
    MPI_Comm comm;
    MPI_Group group;
} comex_igroup_t;

extern void comex_group_init();
extern void comex_group_finalize();
extern comex_igroup_t* comex_get_igroup_from_group(comex_group_t group);

#endif /* _COMEX_GROUPS_H_ */
