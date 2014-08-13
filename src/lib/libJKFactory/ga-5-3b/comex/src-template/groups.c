#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdlib.h>
#include <assert.h>

#include "comex.h"
#include "comex_impl.h"
#include "groups.h"


/* the HEAD of the group linked list */
comex_igroup_t *group_list = NULL;


/* static functions implemented in this file */
static void comex_create_group_and_igroup(
        comex_group_t *id, comex_igroup_t **igroup);
static void comex_igroup_finalize(comex_igroup_t *igroup);


/**
 * Return the comex igroup instance given the group id.
 *
 * The group linked list is searched sequentially until the given group
 * is found. It is an error if this function is called before
 * comex_group_init(). An error occurs if the given group is not found.
 */
comex_igroup_t* comex_get_igroup_from_group(comex_group_t id)
{
    comex_igroup_t *current_group_list_item = group_list;

    assert(group_list != NULL);
    while (current_group_list_item != NULL) {
        if (current_group_list_item->id == id) {
            return current_group_list_item;
        }
        current_group_list_item = current_group_list_item->next;
    }
    comex_error("comex group lookup failed", -1);

    return NULL;
}


/**
 * Creates and associates an comex group with an comex igroup.
 *
 * This does *not* initialize the members of the comex igroup.
 */
static void comex_create_group_and_igroup(
        comex_group_t *id, comex_igroup_t **igroup)
{
    comex_igroup_t *new_group_list_item = NULL;
    comex_igroup_t *last_group_list_item = NULL;

    /* find the last group in the group linked list */
    last_group_list_item = group_list;
    while (last_group_list_item->next != NULL) {
        last_group_list_item = last_group_list_item->next;
    }

    /* create, init, and insert the new node for the linked list */
    new_group_list_item = malloc(sizeof(comex_igroup_t));
    new_group_list_item->id = last_group_list_item->id + 1;
    new_group_list_item->comm = MPI_COMM_NULL;
    new_group_list_item->group = MPI_GROUP_NULL;
    new_group_list_item->next = NULL;
    last_group_list_item->next = new_group_list_item;

    /* return the group id and comex igroup */
    *igroup = new_group_list_item;
    *id = new_group_list_item->id;
}


int comex_group_rank(comex_group_t group, int *rank)
{
    int status;

    comex_igroup_t *igroup = comex_get_igroup_from_group(group);
    status = MPI_Group_rank(igroup->group, rank);
    if (status != MPI_SUCCESS) {
        comex_error("MPI_Group_rank: Failed ", status);
    }

    return COMEX_SUCCESS;
}


int comex_group_size(comex_group_t group, int *size)
{
    int status;

    comex_igroup_t *igroup = comex_get_igroup_from_group(group);
    status = MPI_Group_size(igroup->group, size);
    if (status != MPI_SUCCESS) {
        comex_error("MPI_Group_size: Failed ", status);
    }

    return COMEX_SUCCESS;
}


int comex_group_comm(comex_group_t group, MPI_Comm *comm)
{
    int status;

    comex_igroup_t *igroup = comex_get_igroup_from_group(group);
    *comm = igroup->comm;

    return COMEX_SUCCESS;
}


int comex_group_translate_world(comex_group_t group, int group_rank, int *world_rank)
{
    if (COMEX_GROUP_WORLD == group) {
        *world_rank = group_rank;
    }
    else {
        comex_igroup_t *igroup = comex_get_igroup_from_group(group);
        comex_igroup_t *world_igroup = comex_get_igroup_from_group(COMEX_GROUP_WORLD);
        int status = MPI_Group_translate_ranks(
                igroup->group, 1, &group_rank, world_igroup->group, world_rank);
        if (status != MPI_SUCCESS) {
            comex_error("MPI_Group_translate_ranks: Failed ", status);
        }
    }

    return COMEX_SUCCESS;
}


/**
 * Destroys the given comex igroup.
 */
static void comex_igroup_finalize(comex_igroup_t *igroup)
{
    int status;

    assert(igroup);

    if (igroup->group != MPI_GROUP_NULL) {
        status = MPI_Group_free(&igroup->group);
        if (status != MPI_SUCCESS) {
            comex_error("MPI_Group_free: Failed ", status);
        }
    }
    
    if (igroup->comm != MPI_COMM_NULL) {
        status = MPI_Comm_free(&igroup->comm);
        if (status != MPI_SUCCESS) {
            comex_error("MPI_Comm_free: Failed ", status);
        }
    }
}


int comex_group_free(comex_group_t id)
{
    comex_igroup_t *current_group_list_item = group_list;
    comex_igroup_t *previous_group_list_item = NULL;

    /* find the group to free */
    while (current_group_list_item != NULL) {
        if (current_group_list_item->id == id) {
            break;
        }
        previous_group_list_item = current_group_list_item;
        current_group_list_item = current_group_list_item->next;
    }
    /* make sure we found a group */
    assert(current_group_list_item != NULL);
    /* remove the group from the linked list */
    if (previous_group_list_item != NULL) {
        previous_group_list_item->next = current_group_list_item->next;
    }
    /* free the group */
    comex_igroup_finalize(current_group_list_item);
    free(current_group_list_item);

    return COMEX_SUCCESS;
}


int comex_group_create(
        int n, int *pid_list, comex_group_t id_parent, comex_group_t *id_child)
{
    int status;
    int grp_me;
    comex_igroup_t *igroup_child = NULL;
    MPI_Group      *group_child = NULL;
    MPI_Comm       *comm_child = NULL;
    comex_igroup_t *igroup_parent = NULL;
    MPI_Group      *group_parent = NULL;
    MPI_Comm       *comm_parent = NULL;

    /* create the node in the linked list of groups and */
    /* get the child's MPI_Group and MPI_Comm, to be populated shortly */
    comex_create_group_and_igroup(id_child, &igroup_child);
    group_child = &(igroup_child->group);
    comm_child  = &(igroup_child->comm);

    /* get the parent's MPI_Group and MPI_Comm */
    igroup_parent = comex_get_igroup_from_group(id_parent);
    group_parent = &(igroup_parent->group);
    comm_parent  = &(igroup_parent->comm);

    status = MPI_Group_incl(*group_parent, n, pid_list, group_child);
    if (status != MPI_SUCCESS) {
        comex_error("MPI_Group_incl: Failed ", status);
    }

    {
        MPI_Comm comm, comm1, comm2;
        int lvl=1, local_ldr_pos;
        MPI_Group_rank(*group_child, &grp_me);
        if (grp_me == MPI_UNDEFINED) {
            /* FIXME: keeping the group around for now */
            return COMEX_SUCCESS;
        }
        /* SK: sanity check for the following bitwise operations */
        assert(grp_me>=0);
        MPI_Comm_dup(MPI_COMM_SELF, &comm); /* FIXME: can be optimized away */
        local_ldr_pos = grp_me;
        while(n>lvl) {
            int tag=0;
            int remote_ldr_pos = local_ldr_pos^lvl;
            if (remote_ldr_pos < n) {
                int remote_leader = pid_list[remote_ldr_pos];
                MPI_Comm peer_comm = *comm_parent;
                int high = (local_ldr_pos<remote_ldr_pos)?0:1;
                MPI_Intercomm_create(
                        comm, 0, peer_comm, remote_leader, tag, &comm1);
                MPI_Comm_free(&comm);
                MPI_Intercomm_merge(comm1, high, &comm2);
                MPI_Comm_free(&comm1);
                comm = comm2;
            }
            local_ldr_pos &= ((~0)^lvl);
            lvl<<=1;
        }
        *comm_child = comm;
        /* cleanup temporary group (from MPI_Group_incl above) */
        MPI_Group_free(group_child);
        /* get the actual group associated with comm */
        MPI_Comm_group(*comm_child, group_child);
    }

    return COMEX_SUCCESS;
}


/**
 * Initialize group linked list. Prepopulate with world group.
 */
void comex_group_init() 
{
    /* create the head of the group linked list */
    assert(group_list == NULL);
    group_list = malloc(sizeof(comex_igroup_t));
    group_list->id = COMEX_GROUP_WORLD;
    group_list->next = NULL;

    /* save MPI world group and communicatior in COMEX_GROUP_WORLD */
    group_list->comm = l_state.world_comm;
    MPI_Comm_group(group_list->comm, &(group_list->group));
}


void comex_group_finalize()
{
    comex_igroup_t *current_group_list_item = group_list;
    comex_igroup_t *previous_group_list_item = NULL;

    /* don't free the world group (the list head) */
    current_group_list_item = current_group_list_item->next;

    while (current_group_list_item != NULL) {
        previous_group_list_item = current_group_list_item;
        current_group_list_item = current_group_list_item->next;
        comex_igroup_finalize(previous_group_list_item);
        free(previous_group_list_item);
    }

    /* ok, now free the world group, but not the world comm */
    MPI_Group_free(&(group_list->group));
    free(group_list);
    group_list = NULL;
}
