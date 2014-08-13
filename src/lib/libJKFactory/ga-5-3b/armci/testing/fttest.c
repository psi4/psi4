#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif

#include "armci.h"
#include "message.h"

#define MAXPROCS 128
#define SIZE_    1024
#define FAILURE_SIZE_    512

int me, nproc;
double *ptr_arr[MAXPROCS];
long size;
void do_work(int sz)
{
  int i;
  static int d = 1;
  for (i = 0; i < sz; i++) {
    *((double *)(ptr_arr[me]) + i) = i + 1.12 * d++;
  }
}
static double time_array[100], time_array1[100], t1;
int main(int argc, char *argv[])
{
  int rc, i, j = 0, rid, ret;
  armci_ckpt_ds_t ckptds;
  ARMCI_Group grp;

  armci_msg_init(&argc, &argv);
  ARMCI_Init_args(&argc, &argv);
  nproc = armci_msg_nproc();
  me = armci_msg_me();

  if (me == 0) {
    if (nproc > MAXPROCS) {
      ARMCI_Error("nproc > MAXPROCS", nproc);
    }
    else {
      printf("ARMCI test program (%d processes)\n", nproc);
      fflush(stdout);
      sleep(1);
    }

  }
  armci_init_checkpoint2();
  ARMCI_Group_get_world(&grp);
  size = SIZE_;
  rc = ARMCI_Malloc((void **)ptr_arr, size * 8);
  printf("ARMCI test program (%d processes)\n", nproc);
  fflush(stdout);
  for (size = 1; size <= SIZE_; size *= 2) {
    t1 = MPI_Wtime();
    for (i = 0; i < 5; i++) {
      for (rc = 0; rc < 15; rc++) {
        do_work(size);
      }
    }
    time_array[j++] = MPI_Wtime() - t1;
    ARMCI_Barrier();
    printf("%d:done for size %ld\n", me, size);
    fflush(stdout);
  }

  (void)ARMCI_Ckpt_create_ds(&ckptds, 1);
  ckptds.ptr_arr[0] = ptr_arr[me];
  ckptds.sz[0] = SIZE_ * 8;
  rid = ARMCI_Ckpt_init(NULL, &grp, 1, 0, &ckptds);
  printf("%d: After ARMCI_Ckpt_init(): \n", me);

  j = 0;
  for (size = 128; size <= SIZE_; size *= 2) {

    int rc;
    int simulate_restart = 1;
    t1 = MPI_Wtime();

    ret = ARMCI_Ckpt(rid);
    if (ret == ARMCI_CKPT) {
      printf("%d: Performed CHECKPOINT @ size=%ld\n", me, size);
    }
    else if (ret == ARMCI_RESTART) {
      simulate_restart = 0;
      printf("%d: Performed RESTART @ size=%ld\n", me, size);
    }

    for (i = 0; i < 5; i++) {
      for (rc = 0; rc < 15; rc++)
        if (i == 3 && rc == 10) {
        }
      do_work(size);
    }

    time_array1[j++] = MPI_Wtime() - t1;
    sleep(1);

    if (simulate_restart && size == FAILURE_SIZE_) {
      printf("%d: Simulating FAILURE @ size = %d\n", me, size);
      ARMCI_Restart_simulate(rid, 1);
    }

    printf("%d: DONE for size=%ld regular=%f withckpt=%f\n\n",
           me, size, time_array[j-1], time_array1[j-1]);
    fflush(stdout);

  }

  ARMCI_Ckpt_finalize(rid);

  printf("Before Finalize()\n");
  ARMCI_Barrier();
  ARMCI_Finalize();
  armci_msg_finalize();
  return(0);
}
