/**
 * Private header file for comex groups backed by MPI_comm.
 *
 * The rest of the comex group functions are defined in the public comex.h.
 * These are the private group functions.
 *
 * @author Jeff Daily
 */
#ifndef _COMEX_GROUPS_H_
#define _COMEX_GROUPS_H_

#include <mpi.h>

#include "comex.h"

typedef struct group_link {
    struct group_link *next;
    comex_group_t id;
    MPI_Comm comm;
    MPI_Group group;
} comex_igroup_t;

extern void comex_group_init();
extern void comex_group_finalize();
extern comex_igroup_t* comex_get_igroup_from_group(comex_group_t group);

/* verify that proc is part of group
 * translate proc to world group
 * change group to world group */
#define CHECK_GROUP(GROUP,PROC) do {                            \
    int size;                                                   \
    assert(COMEX_SUCCESS == comex_group_size(GROUP,&size));     \
    assert(PROC >= 0);                                          \
    assert(PROC < size);                                        \
    if (COMEX_GROUP_WORLD != GROUP) {                           \
        int world_proc;                                         \
        comex_group_translate_world(GROUP, PROC, &world_proc);  \
        PROC = world_proc;                                      \
        GROUP = COMEX_GROUP_WORLD;                              \
    }                                                           \
} while(0)

#endif /* _COMEX_GROUPS_H_ */
