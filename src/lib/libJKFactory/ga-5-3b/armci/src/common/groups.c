#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: groups.c,v 1.4.6.2 2007-08-15 08:37:16 manoj Exp $ */


#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif

#ifndef MPI
#  define MPI
#endif
#include "armcip.h"
#include "message.h"

#define DEBUG_ 0

MPI_Comm ARMCI_COMM_WORLD; /*dup of MPI_COMM_WORLD. Initialized first thing in ARMCI_Init*/

ARMCI_Group ARMCI_Default_Proc_Group = 0;
ARMCI_Group ARMCI_World_Proc_Group = 0;

typedef struct group_list_struct {
    ARMCI_Group group;
    ARMCI_iGroup igroup;
    struct group_list_struct *next;
} group_list_t;

group_list_t *group_list = NULL;

ARMCI_iGroup* armci_get_igroup_from_group(ARMCI_Group *group)
{
    group_list_t *current_group_list_item = group_list;

    assert(group_list != NULL);
    while (current_group_list_item != NULL) {
        if (current_group_list_item->group == *group) {
            return &current_group_list_item->igroup;
        }
        current_group_list_item = current_group_list_item->next;
    }
    armci_die("ARMCI_Group lookup failed", -1);
    return NULL;
}

static void armci_create_group_and_igroup(ARMCI_Group *group, ARMCI_iGroup **igroup)
{
    group_list_t *new_group_list_item = NULL;
    group_list_t *last_group_list_item = NULL;

    /* create the new group in the linked list */
    last_group_list_item = group_list;
    while (last_group_list_item->next != NULL) {
        last_group_list_item = last_group_list_item->next;
    }

    new_group_list_item = malloc(sizeof(group_list_t));
    new_group_list_item->group = last_group_list_item->group + 1;
    new_group_list_item->next = NULL;
    *igroup = &new_group_list_item->igroup;
    *group = new_group_list_item->group;
    last_group_list_item->next = new_group_list_item;
}

#ifdef ARMCI_GROUP
void ARMCI_Bcast_(void *buffer, int len, int root, ARMCI_Group *group)
{
  armci_msg_group_bcast_scope(SCOPE_ALL, buffer, len, 
			      ARMCI_Absolute_id(group, root),
			      group);
}
#else
void ARMCI_Bcast_(void *buffer, int len, int root, ARMCI_Comm comm)
{
    int result;
    MPI_Comm_compare(comm, ARMCI_COMM_WORLD, &result);
    if(result == MPI_IDENT)  armci_msg_brdcst(buffer, len, root); 
    else MPI_Bcast(buffer, len, MPI_BYTE, root, (MPI_Comm)comm);
}
#endif

int ARMCI_Group_rank(ARMCI_Group *group, int *rank)
{
    ARMCI_iGroup *igroup = armci_get_igroup_from_group(group);
#ifdef ARMCI_GROUP
    if(!igroup) return MPI_ERR_GROUP;
    *rank = igroup->grp_attr.grp_me;
    return MPI_SUCCESS;
#else
    return MPI_Group_rank((MPI_Group)(igroup->igroup), rank);
#endif
}

void ARMCI_Group_size(ARMCI_Group *group, int *size)
{
    ARMCI_iGroup *igroup = armci_get_igroup_from_group(group);
#ifdef ARMCI_GROUP
    *size = igroup->grp_attr.nproc;
#else
    MPI_Group_size((MPI_Group)(igroup->igroup), size);
#endif
}

int ARMCI_Absolute_id(ARMCI_Group *group,int group_rank)
{
    int abs_rank,status;
    ARMCI_iGroup *igroup = armci_get_igroup_from_group(group);
#ifdef ARMCI_GROUP
    assert(group_rank < igroup->grp_attr.nproc);
    return igroup->grp_attr.proc_list[group_rank];
#else
    MPI_Group grp;
    status = MPI_Comm_group(ARMCI_COMM_WORLD,&grp);
    MPI_Group_translate_ranks(igroup->igroup,1,&group_rank,grp,&abs_rank);
    return(abs_rank);
#endif
}

void ARMCI_Group_set_default(ARMCI_Group *group) 
{
    ARMCI_Default_Proc_Group = *group;
}

void ARMCI_Group_get_default(ARMCI_Group *group_out)
{
    *group_out = ARMCI_Default_Proc_Group;
}

void ARMCI_Group_get_world(ARMCI_Group *group_out)
{
    *group_out = ARMCI_World_Proc_Group;
}

static void get_group_clus_id(ARMCI_iGroup *igroup, int grp_nproc, 
                              int *grp_clus_id)
{
#ifdef ARMCI_GROUP
    int i;
    assert(grp_nproc<=igroup->grp_attr.nproc);
    for(i=0; i<grp_nproc; i++) {
      grp_clus_id[i] = armci_clus_id(igroup->grp_attr.proc_list[i]);
    }
#else
    int i, *ranks1, *ranks2;
    MPI_Group group2;
    
    /* Takes the list of processes from one group and attempts to determine
     * the corresponding ranks in a second group (here, ARMCI_COMM_WORLD) */

    ranks1 = (int *)malloc(2*grp_nproc*sizeof(int));
    ranks2 = ranks1 + grp_nproc;
    for(i=0; i<grp_nproc; i++) ranks1[i] = i;
    MPI_Comm_group(ARMCI_COMM_WORLD, &group2);
    MPI_Group_translate_ranks(igroup->igroup, grp_nproc, ranks1, group2, ranks2);
    
    /* get the clus_id of processes */
    for(i=0; i<grp_nproc; i++) grp_clus_id[i] = armci_clus_id(ranks2[i]);
    free(ranks1);
#endif
}

/**
 * Construct the armci_clus_t arrays of structs for a given group. 
 *
 * @param grp_nclus_nodes OUT #clus_info objects returned
 * @param igroup IN Group whose clus_info needs to be constructed
 * @return array of armci_clus_t objects
 */
static armci_clus_t* group_construct_clusinfo(int *grp_nclus_nodes, ARMCI_Group *group) {
  armci_clus_t *grp_clus_info=NULL;
  int i, *grp_clus_id, cluster, clus_id, grp_nproc, grp_nclus;
  ARMCI_iGroup *igroup = armci_get_igroup_from_group(group);

  ARMCI_Group_size(group, &grp_nproc);

  /* get the cluster_id of processes in the group */
  grp_clus_id = (int *)malloc(grp_nproc*sizeof(int));
  get_group_clus_id(igroup, grp_nproc, grp_clus_id);
       
  /* first find out how many cluster nodes we got for this group */
  grp_nclus=1;
  for(i=0; i<grp_nproc-1; i++) {
    if(grp_clus_id[i] != grp_clus_id[i+1]) ++grp_nclus;
  }
  *grp_nclus_nodes = grp_nclus;
  grp_clus_info = (armci_clus_t*)malloc(grp_nclus*sizeof(armci_clus_t));
  if(!grp_clus_info)armci_die("malloc failed for grp_clusinfo",grp_nclus);

  cluster = 1;
  clus_id = grp_clus_id[0];
  grp_clus_info[0].nslave = 1;
  grp_clus_info[0].master = 0;
  strcpy(grp_clus_info[0].hostname, armci_clus_info[clus_id].hostname);

  for(i=1; i<grp_nproc; i++) {
    if(grp_clus_id[i-1] == grp_clus_id[i]) 
      ++grp_clus_info[cluster-1].nslave;
    else {
      clus_id = grp_clus_id[i];
      grp_clus_info[cluster].nslave = 1;
      grp_clus_info[cluster].master = i;
      strcpy(grp_clus_info[cluster].hostname, 
	     armci_clus_info[clus_id].hostname);
      ++cluster;
    }
  }

  free(grp_clus_id);
  if(grp_nclus != cluster)
    armci_die("inconsistency processing group clusterinfo", grp_nclus);

#   if DEBUG_
  {
    int i,j;
    for(i=0; i<cluster;i++) {
      printf("%d: Cluster %d: Master=%d, #Slaves=%d, HostName=%s\n", 
	     grp_nclus, i, grp_clus_info[i].master, 
	     grp_clus_info[i].nslave, grp_clus_info[i].hostname);
      fflush(stdout);
    }
  }
#   endif  
  return grp_clus_info;
}

/**
 * Group cluster information "grp_clus_info" (similar to armci_clus_info)
 */
static void group_process_list(ARMCI_Group *group, 
                               armci_grp_attr_t *grp_attr) {
    ARMCI_iGroup *igroup = armci_get_igroup_from_group(group);
#ifndef ARMCI_GROUP
    ARMCI_Comm comm = igroup->icomm;
#endif

    int grp_me, grp_nproc, grp_nclus, grp_clus_me;
    armci_clus_t *grp_clus_info=NULL;
#ifdef CLUSTER
    int i, len, root=0;
#endif
    
#ifndef ARMCI_GROUP
    if(comm==MPI_COMM_NULL || igroup->igroup==MPI_GROUP_NULL) 
       armci_die("group_process_list: NULL COMMUNICATOR",0);
#endif
    
    ARMCI_Group_rank(group, &grp_me);
    ARMCI_Group_size(group, &grp_nproc);
    
#ifdef CLUSTER
#   ifdef ARMCI_GROUP
    /*all processes construct the clus_info structure in parallel*/
    grp_clus_info = group_construct_clusinfo(&grp_nclus, group);
#   else
    /* process 0 gets group cluster information: grp_nclus, grp_clus_info */
    if(grp_me == 0) {
      grp_clus_info = group_construct_clusinfo(&grp_nclus, group);
    }

    /* process 0 broadcasts group cluster information */
    len = sizeof(int);
    ARMCI_Bcast_(&grp_nclus, len, root, comm);
    
    if(grp_me){
       /* allocate memory */
       grp_clus_info = (armci_clus_t*)malloc(grp_nclus*sizeof(armci_clus_t));
       if(!armci_clus_info)armci_die("malloc failed for clusinfo",armci_nclus);
    }
    
    len = sizeof(armci_clus_t)*grp_nclus;
    ARMCI_Bcast_(grp_clus_info, len, root, comm);
#   endif
    /* determine current group cluster node id by comparing me to master */
    grp_clus_me =  grp_nclus-1;
    for(i =0; i< grp_nclus-1; i++) {
       if(grp_me < grp_clus_info[i+1].master){
          grp_clus_me=i;
          break;
       }
    }
#else /* !CLUSTER */
    grp_clus_me = 0;
    grp_nclus = 1;
    grp_clus_info = (armci_clus_t*)malloc(grp_nclus*sizeof(armci_clus_t));
    if(!grp_clus_info)armci_die("malloc failed for clusinfo",grp_nclus);
    strcpy(grp_clus_info[0].hostname, armci_clus_info[0].hostname);
    grp_clus_info[0].master=0;
    grp_clus_info[0].nslave=grp_nproc;
#endif /* CLUSTER */
#ifdef ARMCI_GROUP
    /*Set in ARMCI_Group_create. ARMCI_Group_rank is used before
      setting this field. So moving it there in the generic
      implementation.*/
#else
    grp_attr->grp_me        = grp_me;
#endif
    grp_attr->grp_clus_info = grp_clus_info;
    grp_attr->grp_nclus     = grp_nclus;
    grp_attr->grp_clus_me   = grp_clus_me;
}

/* attribute caching: group_cluster_information and memory_offset should 
   be cached in group data structure */
static void armci_cache_attr(ARMCI_Group *group) {
    armci_grp_attr_t *grp_attr;
    ARMCI_iGroup *igroup = armci_get_igroup_from_group(group);

    /* allocate storage for the attribute content. Note: Attribute contents 
       should be stored in persistent memory */
    grp_attr = &(igroup->grp_attr); 
    
    /* get group cluster information and  grp_attr */
    group_process_list(group, grp_attr);
}

armci_grp_attr_t *ARMCI_Group_getattr(ARMCI_Group *group)
{
    ARMCI_iGroup *igroup = armci_get_igroup_from_group(group);
    return(&(igroup->grp_attr));

}

static void armci_igroup_finalize(ARMCI_iGroup *igroup) {
#ifdef ARMCI_GROUP
    int world_me, i;

    world_me = armci_msg_me();
    for(i=0; i<igroup->grp_attr.nproc; i++) {
      if(igroup->grp_attr.proc_list[i] == world_me) {
	break;
      }
    }
    if(i==igroup->grp_attr.nproc) {
      return; /*not in group to be freed*/
    }

    assert(igroup);
    free(igroup->grp_attr.grp_clus_info);
    free(igroup->grp_attr.proc_list);
    igroup->grp_attr.nproc = 0;
#else
    int rv;

    assert(igroup);
    /*the following was causing seg fault*/
    /*free(igroup->grp_attr.grp_clus_info);*/
    
    rv=MPI_Group_free(&(igroup->igroup));
    if(rv != MPI_SUCCESS) armci_die("MPI_Group_free: Failed ",armci_me);
    
    if(igroup->icomm != MPI_COMM_NULL) {
      rv = MPI_Comm_free( (MPI_Comm*)&(igroup->icomm) );
      if(rv != MPI_SUCCESS) armci_die("MPI_Comm_free: Failed ",armci_me);
    }
#endif
}

void ARMCI_Group_free(ARMCI_Group *group) {
    group_list_t *current_group_list_item = group_list;
    group_list_t *previous_group_list_item = NULL;

    /* find the group to free */
    while (current_group_list_item != NULL) {
        if (current_group_list_item->group == *group) {
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
    armci_igroup_finalize(&current_group_list_item->igroup);
    free(current_group_list_item);
}

/*
  Create a child group for to the given group.
  @param n IN #procs in this group (<= that in group_parent)
  @param pid_list IN The list of proc ids (w.r.t. group_parent)
  @param group_out OUT Handle to store the created group
  @param group_parent IN Parent group 
 */
void ARMCI_Group_create_child(int n, int *pid_list, ARMCI_Group *group_out,
			      ARMCI_Group *grp_parent)
{
    int grp_me;
    ARMCI_iGroup *igroup = NULL;

#ifdef ARMCI_GROUP
    int i, world_me, parent_grp_me;
    armci_grp_attr_t *grp_attr = NULL;
#else
    int rv;
    ARMCI_iGroup *igroup_parent = NULL;
    MPI_Group *group_parent = NULL;
    MPI_Comm *comm_parent = NULL;
#endif

    armci_create_group_and_igroup(group_out, &igroup);

#ifdef ARMCI_GROUP
    grp_attr = &igroup->grp_attr;
    ARMCI_Group_rank(grp_parent, &parent_grp_me);
    for(i=0; i<n; i++) {
      if(pid_list[i] == parent_grp_me) {
        break;
      }
    }
    if(i==n) {
      /*this initialization is used in group free*/
      grp_attr->nproc=0;
      grp_attr->proc_list = NULL; 
      return; /*not in group to be created*/
    }
    for(i=0; i<n-1;i++) {
       if(pid_list[i] > pid_list[i+1]){
         armci_die("ARMCI_Group_create: Process ids are not sorted ",armci_me);
         break;
       }
    }
    grp_attr->grp_clus_info = NULL;
    grp_attr->nproc = n;
    grp_attr->proc_list = (int *)malloc(n*sizeof(int));
    assert(grp_attr->proc_list!=NULL);
    for(i=0; i<n; i++)  {
      grp_attr->proc_list[i] = ARMCI_Absolute_id(grp_parent,pid_list[i]); 
    }
    world_me = armci_msg_me();
    grp_attr->grp_me = grp_me = MPI_UNDEFINED;
    for(i=0; i<n; i++) {
      if(igroup->grp_attr.proc_list[i] == world_me) {
        grp_attr->grp_me = grp_me = i;
        break;
      }
    }
    if(grp_me != MPI_UNDEFINED) armci_cache_attr(group_out);
    armci_msg_group_barrier(group_out);
#else
    igroup_parent = armci_get_igroup_from_group(grp_parent);
    /* NOTE: default group is the parent group */
    group_parent = &(igroup_parent->igroup);
    comm_parent  = &(igroup_parent->icomm);

    rv=MPI_Group_incl(*group_parent, n, pid_list, &(igroup->igroup));
    if(rv != MPI_SUCCESS) armci_die("MPI_Group_incl: Failed ",armci_me);
    
    {
      MPI_Comm comm, comm1, comm2;
      int lvl=1, local_ldr_pos;
      MPI_Group_rank((MPI_Group)(igroup->igroup), &grp_me);
      if(grp_me == MPI_UNDEFINED) {
	igroup->icomm = MPI_COMM_NULL; /*FIXME: keeping the group around for now*/	
	return;
      }
      assert(grp_me>=0); /*SK: sanity check for the following bitwise operations*/
      MPI_Comm_dup(MPI_COMM_SELF, &comm); /*FIXME: can be optimized away*/
      local_ldr_pos = grp_me;
      while(n> lvl) {
	int tag=0;
	int remote_ldr_pos = local_ldr_pos^lvl;
	if(remote_ldr_pos < n) {
	  int remote_leader = pid_list[remote_ldr_pos];
	  MPI_Comm peer_comm = *comm_parent;
	  int high = (local_ldr_pos<remote_ldr_pos)?0:1;
	  MPI_Intercomm_create(comm, 0, peer_comm, remote_leader, tag, &comm1);
	  MPI_Comm_free(&comm);
	  MPI_Intercomm_merge(comm1, high, &comm2);
	  MPI_Comm_free(&comm1);
	  comm = comm2;
	}
	local_ldr_pos &= ((~0)^lvl);
	lvl<<=1;
      }
      igroup->icomm = comm;
      MPI_Group_free(&igroup->igroup); /*cleanup temporary group*/
      MPI_Comm_group(igroup->icomm, &igroup->igroup); /*the group associated with comm*/
      igroup->grp_attr.grp_clus_info=NULL;
      /* processes belong to this group should cache attributes */
      armci_cache_attr(group_out);
    }    

#endif
}

void ARMCI_Group_create(int n, int *pid_list, ARMCI_Group *group_out) {
  ARMCI_Group_create_child(n, pid_list, group_out, (ARMCI_Group *)&ARMCI_Default_Proc_Group);
}

void armci_group_init() 
{
#ifdef ARMCI_GROUP
    int i;
#else
    int grp_me;
#endif
    ARMCI_iGroup *igroup;

    /* Initially, World group is the default group */
    ARMCI_World_Proc_Group = 0;
    ARMCI_Default_Proc_Group = 0;

    /* create the head of the group linked list */
    assert(group_list == NULL);
    group_list = malloc(sizeof(group_list_t));
    group_list->group = ARMCI_World_Proc_Group;
    group_list->next = NULL;
    igroup = &group_list->igroup;

#ifdef ARMCI_GROUP
    /*setup the world proc group*/
    igroup->grp_attr.nproc = armci_msg_nproc();
    igroup->grp_attr.grp_me = armci_msg_me();
    igroup->grp_attr.proc_list = (int *)malloc(igroup->grp_attr.nproc*sizeof(int));
    assert(igroup->grp_attr.proc_list != NULL);
    for(i=0; i<igroup->grp_attr.nproc; i++) {
      igroup->grp_attr.proc_list[i] = i;
    } 
    igroup->grp_attr.grp_clus_info = NULL;
    armci_cache_attr(&ARMCI_World_Proc_Group);
#else
    /* save MPI world group and communicatior in ARMCI_World_Proc_Group */
    igroup->icomm = ARMCI_COMM_WORLD;
    MPI_Comm_group(ARMCI_COMM_WORLD, &(igroup->igroup));

    /* processes belong to this group should cache attributes */
    MPI_Group_rank((MPI_Group)(igroup->igroup), &grp_me);
    if(grp_me != MPI_UNDEFINED) {
       armci_cache_attr(&ARMCI_World_Proc_Group);
    }
#endif    
}

void armci_group_finalize() {
    group_list_t *current_group_list_item = group_list;
    group_list_t *previous_group_list_item = NULL;

    /* don't free the world group (the list head) */
    current_group_list_item = current_group_list_item->next;

    while (current_group_list_item != NULL) {
        previous_group_list_item = current_group_list_item;
        current_group_list_item = current_group_list_item->next;
        armci_igroup_finalize(&previous_group_list_item->igroup);
        free(previous_group_list_item);
    }
}

/*
  ISSUES:
  1. Make sure ARMCI_Group_free frees the attribute data structures 
  2. replace malloc with, kr_malloc using local_context.
*/
