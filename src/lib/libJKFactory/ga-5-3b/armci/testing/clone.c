#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdio.h>
#include <portals/portals3.h>
#include <portals/nal.h>
#include <mpi.h>

#define FORK_BEFORE_NI_INIT
#ifndef FORK_BEFORE_NI_INIT
#define FORK_AFTER_NI_INIT
#endif

char child_stack[256*1024];
char *child_stack_top = &child_stack[256*1024-1];
int iv;

int
server(void *arg)
{

  int ret;
  int num_interfaces;
  ptl_handle_ni_t nih;
  ptl_handle_eq_t eqh;
  ptl_ni_limits_t ptl_limits;
  ptl_event_t ev_t;
  ptl_event_t *ev = &ev_t;


  printf("IN SERVER\n");

  if ((ret = PtlInit(&num_interfaces)) != PTL_OK) {
    printf("%s: PtlInit failed: %d\n", FUNCTION_NAME, ret);
    exit(1);
  }
  printf("%s: PtlInit succeeds (%d:%d)\n", FUNCTION_NAME, ret, num_interfaces);

  if (((ret = PtlNIInit(
                IFACE_FROM_BRIDGE_AND_NALID(PTL_BRIDGE_UK, PTL_IFACE_SS),
                PTL_PID_ANY, NULL, &ptl_limits, &nih)) != PTL_OK) && (ret != PTL_IFACE_DUP)) {
    printf("%s: PtlNIInit failed: %d\n", FUNCTION_NAME, ret);
    exit(1);
  }
  printf("%s: PtlNIInit succeeds (%d)\n", FUNCTION_NAME, ret);

  if ((ret = PtlEQAlloc(nih, 4096, NULL, &eqh)) != PTL_OK) {
    printf("%s: PtlEQAlloc failed: %d\n",
           FUNCTION_NAME, ret);
    exit(1);
  }
  iv = 11;
  ret = PtlEQWait(nih, ev);
  printf("%s: PtlEQAlloc succeeds\n", FUNCTION_NAME);
  printf("%d\n", iv);
  iv = 13;

  while (1);
}

int
main(int argc, char **argv, char **envp)
{
  int ret;
  pid_t child;
  int status;
  iv = 12;

  child = clone(server, (void *)child_stack_top,
                CLONE_THREAD | CLONE_SIGHAND | CLONE_VM, NULL);

  if (child == -1) {
    perror("clone");
    exit(1);
  }
  printf("clone returns...(ret=%d)\n", child);
  while (iv != 11);
  printf("\nbetween after %d\n", iv);


  MPI_Init(&argc, &argv);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  printf("waiting...\n");
  waitpid(-1, &status, __WALL);
  printf("\nafter %d\n", iv);
  printf("done (%d)\n", status);
}

