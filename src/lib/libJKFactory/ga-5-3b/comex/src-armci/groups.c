#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdlib.h>
#include <assert.h>

#include "armci.h"
#include "message.h"
#include "comex.h"

/* ARMCI has the notion of a default group and a world group. */
ARMCI_Group ARMCI_Default_Proc_Group = 0;


int ARMCI_Group_rank(ARMCI_Group *id, int *rank)
{
    return comex_group_rank(*id, rank);
}


void ARMCI_Group_size(ARMCI_Group *id, int *size)
{
    comex_group_size(*id, size);
}


int ARMCI_Absolute_id(ARMCI_Group *id, int group_rank)
{
    int world_rank;
    assert(COMEX_SUCCESS == 
            comex_group_translate_world(*id, group_rank, &world_rank));
    return world_rank;
}


void ARMCI_Group_set_default(ARMCI_Group *id) 
{
    ARMCI_Default_Proc_Group = *id;
}


void ARMCI_Group_get_default(ARMCI_Group *group_out)
{
    *group_out = ARMCI_Default_Proc_Group;
}


void ARMCI_Group_get_world(ARMCI_Group *group_out)
{
    *group_out = COMEX_GROUP_WORLD;
}


void ARMCI_Group_free(ARMCI_Group *id)
{
    comex_group_free(*id);
}


void ARMCI_Group_create_child(
        int n, int *pid_list, ARMCI_Group *id_child, ARMCI_Group *id_parent)
{
    comex_group_create(n, pid_list, *id_parent, id_child);
}


void ARMCI_Group_create(int n, int *pid_list, ARMCI_Group *group_out)
{
    comex_group_create(n, pid_list, ARMCI_Default_Proc_Group, group_out);
}
