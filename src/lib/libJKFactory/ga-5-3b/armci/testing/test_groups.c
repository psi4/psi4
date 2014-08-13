#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: test_groups.c,v 1.3 2004-11-22 20:29:53 manoj Exp $ */
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#elif HAVE_WINDOWS_H
#   include <windows.h>
#   define sleep(x) Sleep(1000*(x))
#endif

#include "armci.h"
#include "message.h"
/*#include "armcip.h"*/

#define MAXDIMS 7
#define MAXPROC 128
#define MINPROC 4

/***************************** macros ************************/
#define COPY(src, dst, bytes) memcpy((dst),(src),(bytes))
#define ABS(a)   (((a) >= 0) ? (a) : (-(a)))
#define MIN(a,b) (((a)<(b)) ? (a) : (b))

/***************************** global data *******************/
int me, nproc;
void *work[MAXPROC]; /* work array for propagating addresses */

#ifdef PVM
void pvm_init(int argc, char *argv[])
{
  int mytid, mygid, ctid[MAXPROC];
  int np, i;

  mytid = pvm_mytid();
  if ((argc != 2) && (argc != 1)) {
    goto usage;
  }
  if (argc == 1) {
    np = 1;
  }
  if (argc == 2)
    if ((np = atoi(argv[1])) < 1) {
      goto usage;
    }
  if (np > MAXPROC) {
    goto usage;
  }

  mygid = pvm_joingroup(MPGROUP);

  if (np > 1)
    if (mygid == 0) {
      i = pvm_spawn(argv[0], argv + 1, 0, "", np - 1, ctid);
    }

  while (pvm_gsize(MPGROUP) < np) {
    sleep(1);
  }

  /* sync */
  pvm_barrier(MPGROUP, np);

  printf("PVM initialization done!\n");

  return;

usage:
  fprintf(stderr, "usage: %s <nproc>\n", argv[0]);
  pvm_exit();
  exit(-1);
}
#endif

void create_array(void *a[], int elem_size, int ndim, int dims[])
{
  int bytes = elem_size, i, rc;

  assert(ndim <= MAXDIMS);
  for (i = 0; i < ndim; i++) {
    bytes *= dims[i];
  }

  rc = ARMCI_Malloc(a, bytes);
  assert(rc == 0);

  assert(a[me]);

}

void destroy_array(void *ptr[])
{
  ARMCI_Barrier();

  assert(!ARMCI_Free(ptr[me]));
}

#define GNUM_A 3
#define GNUM_B 2
#define ELEMS 10

/* to check if a process belongs to this group */
int chk_grp_membership(int rank, ARMCI_Group *grp, int *memberlist)
{
  int i, grp_size;
  ARMCI_Group_size(grp, &grp_size);
  for (i = 0; i < grp_size; i++) if (rank == memberlist[i]) {
      return 1;
    }
  return 0;
}

void test_one_group(ARMCI_Group *group, int *pid_list)
{
  int grp_me, grp_size;
  int i, j, src_proc, dst_proc;
  double *ddst_put[MAXPROC];
  double dsrc[ELEMS];
  int bytes, world_me;

  world_me = armci_msg_me();
  ARMCI_Group_rank(group, &grp_me);
  ARMCI_Group_size(group, &grp_size);
  if (grp_me == 0) {
    printf("GROUP SIZE = %d\n", grp_size);
  }
  printf("%d:group rank = %d\n", me, grp_me);

  src_proc = 0;
  dst_proc = grp_size - 1;

  bytes = ELEMS * sizeof(double);
  ARMCI_Malloc_group((void **)ddst_put, bytes, group);

  for (i = 0; i < ELEMS; i++) {
    dsrc[i] = i * 1.001 * (grp_me + 1);
  }
  for (i = 0; i < ELEMS; i++) {
    ddst_put[grp_me][i] = -1.0;
  }

  armci_msg_group_barrier(group);

  if (grp_me == src_proc) {
    /* NOTE: make sure to specify absolute ids in ARMCI calls */
    ARMCI_Put(dsrc, &ddst_put[dst_proc][0], bytes,
              ARMCI_Absolute_id(group, dst_proc));
  }

  armci_msg_group_barrier(group);
  /* NOTE: make sure to specify absolute ids in ARMCI calls */
  ARMCI_Fence(ARMCI_Absolute_id(group, dst_proc));
  sleep(1);


  /* Verify*/
  if (grp_me == dst_proc) {
    for (j = 0; j < ELEMS; j++) {
      if (ABS(ddst_put[grp_me][j] - j * 1.001 *(src_proc + 1)) > 0.1) {
        printf("\t%d: ddst_put[%d][%d] = %f and expected value is %f\n",
               me, grp_me, j, ddst_put[grp_me][j], j * 1.001 *(src_proc + 1));
        ARMCI_Error("groups: armci put failed...1", 0);
      }
    }
    printf("\n%d(%d): Test O.K. Verified\n", dst_proc, world_me);
  }
  armci_msg_group_barrier(group);
  ARMCI_Free_group(ddst_put[grp_me], group);
}


void test_groups()
{

  int pid_listA[MAXPROC]  = {0, 1, 2};
  int pid_listB[MAXPROC] = {1, 3};
  ARMCI_Group groupA, groupB;

  ARMCI_Barrier();

  ARMCI_Group_create(GNUM_A, pid_listA, &groupA); /* create group 1 */
  ARMCI_Group_create(GNUM_B, pid_listB, &groupB); /* create group 2 */


  /* ------------------------ GROUP A ------------------------- */
  if (chk_grp_membership(me, &groupA, pid_listA)) { /* group A */
    test_one_group(&groupA, pid_listA);
  }

  ARMCI_Barrier();

  /* ------------------------ GROUP B ------------------------- */
  if (chk_grp_membership(me, &groupB, pid_listB)) { /* group B */
    test_one_group(&groupB, pid_listB);
  }

  ARMCI_AllFence();
  ARMCI_Barrier();

  if (me == 0) {
    printf("O.K.\n");
    fflush(stdout);
  }
}

/**
 * Random permutation of 0..n-1 into an array.
 */
void random_permute(int *arr, int n)
{
  int i, j;
  int *vtmp = (int *)malloc(n * sizeof(int));
  assert(vtmp != NULL);
  for (i = 0; i < n; ++i) {
    vtmp[i] = -1;
  }
  for (i = 0; i < n; i++) {
    while (vtmp[j=(rand()%n)] != -1) /*no-op*/;
    assert(vtmp[j] == -1);
    vtmp[j] = 0;
    arr[i] = j;
  }
  free(vtmp);
}

int int_compare(const void *v1, const void *v2)
{
  int i1 = *(int *)v1, i2 = *(int *)v2;
  if (i1 < i2) {
    return -1;
  }
  if (i1 > i2) {
    return +1;
  }
  return 0;
}

/**
 * Test routine for non-collective process group management. This test
 * should not be used with MPI process group implementation.
 */
#define GROUP_SIZE 2
#define MAX_GROUPS (MAXPROC/GROUP_SIZE)
void test_groups_noncollective()
{
  int *pid_lists[MAX_GROUPS];
  int pids[MAXPROC];
  int i, nprocs, world_me;
  ARMCI_Group group;
  int *my_pid_list = NULL, my_grp_size = 0;
  int ngrps;

  ARMCI_Barrier();
  nprocs = armci_msg_nproc();
  world_me = armci_msg_me();

  random_permute(pids, nproc);

  ngrps = nprocs / GROUP_SIZE;

  for (i = 0; i < nprocs / GROUP_SIZE; i++) {
    pid_lists[i] = pids + (i * GROUP_SIZE);
  }

  for (i = 0; i < nprocs; i++) {
    if (pids[i] == world_me) {
      int grp_id = MIN(i / GROUP_SIZE, ngrps - 1);
      my_pid_list = pid_lists[grp_id];
      if (grp_id == ngrps - 1) {
        my_grp_size =  GROUP_SIZE + (nprocs % GROUP_SIZE);
      }
      else {
        my_grp_size = GROUP_SIZE;
      }
    }
  }

  qsort(my_pid_list, my_grp_size, sizeof(int), int_compare);

  ARMCI_Barrier();
  /*now create all these disjoint groups and test them in parallel*/

  ARMCI_Group_create(my_grp_size, my_pid_list, &group);

  test_one_group(&group, my_pid_list);

  ARMCI_Group_free(&group);

  ARMCI_AllFence();
  ARMCI_Barrier();

  if (world_me == 0) {
    printf("O.K.\n");
    fflush(stdout);
  }
}


int main(int argc, char *argv[])
{
  armci_msg_init(&argc, &argv);
  ARMCI_Init_args(&argc, &argv);
  nproc = armci_msg_nproc();
  me = armci_msg_me();

  /*    printf("nproc = %d, me = %d\n", nproc, me);*/

  if (nproc < MINPROC) {
    if (0 == me) {
      printf("Test needs at least %d processors (%d used)\n",
             MINPROC, nproc);
    }
    ARMCI_Barrier();
    armci_msg_finalize();
    exit(0);
  }
  if (nproc > MAXPROC) {
    if (0 == me) {
      printf("Test works for up to %d processors (%d used)\n",
             MAXPROC, nproc);
    }
    ARMCI_Barrier();
    armci_msg_finalize();
    exit(0);
  }

  if (me == 0) {
    printf("ARMCI test program (%d processes)\n", nproc);
    fflush(stdout);
    sleep(1);
  }

  if (me == 0) {
    printf("\n Testing ARMCI Groups!\n\n");
    fflush(stdout);
  }

  test_groups();

  ARMCI_AllFence();
  ARMCI_Barrier();
  if (me == 0) {
    printf("\n Collective groups: Success!!\n");
    fflush(stdout);
  }
  sleep(2);

#ifdef ARMCI_GROUP
  test_groups_noncollective();

  ARMCI_AllFence();
  ARMCI_Barrier();
  if (me == 0) {
    printf("\n Non-collective groups: Success!!\n");
    fflush(stdout);
  }
  sleep(2);
#endif

  ARMCI_Barrier();
  ARMCI_Finalize();
  armci_msg_finalize();
  return(0);
}
