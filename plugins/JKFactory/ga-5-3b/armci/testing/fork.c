#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

#include <portals/portals3.h>
#include <portals/nal.h>
#include <mpi.h>

int
main(int argc, char **argv, char **envp)
{
  int i, ret, npes;
  int num_interfaces;
  ptl_handle_ni_t nih;
  ptl_handle_eq_t eqh;
  ptl_ni_limits_t ptl_limits;
  pid_t child;
  ptl_process_id_t rnk;

  child = fork();
  if ((ret = PtlInit(&num_interfaces)) != PTL_OK) {
    printf("%s: PtlInit failed: %d\n", FUNCTION_NAME, ret);
    exit(1);
  }
  printf("%s: PtlInit succeeds (%d)\n", FUNCTION_NAME, ret);


  if ((ret = PtlNIInit(
               IFACE_FROM_BRIDGE_AND_NALID(PTL_BRIDGE_UK, PTL_IFACE_SS),
               PTL_PID_ANY, NULL, &ptl_limits, &nih)) != PTL_OK) {
    printf("%s: PtlNIInit 1 failed: %d\n", FUNCTION_NAME, ret);
  }

  if ((ret = PtlNIFini(nih)) != PTL_OK) {
    printf("%s: PtlNIFini failed: %d\n", FUNCTION_NAME, ret);
  }
  PtlFini();

  if ((ret = PtlInit(&num_interfaces)) != PTL_OK) {
    printf("%s: PtlInit failed: %d\n", FUNCTION_NAME, ret);
    exit(1);
  }
  if ((ret = PtlNIInit(
               IFACE_FROM_BRIDGE_AND_NALID(PTL_BRIDGE_UK, PTL_IFACE_SS),
               PTL_PID_ANY, NULL, &ptl_limits, &nih)) != PTL_OK) {
    printf("%s: PtlNIInit 2 failed: %d\n", FUNCTION_NAME, ret);
    exit(1);
  }
#if 0
  if ((ret = PtlNIInit(
               IFACE_FROM_BRIDGE_AND_NALID(PTL_BRIDGE_UK, PTL_IFACE_SS),
               PTL_PID_ANY, NULL, &ptl_limits, &nih)) != PTL_OK) {
    printf("%s: PtlNIInit failed: %d\n", FUNCTION_NAME, ret);
    exit(1);
  }
#endif
  printf("%s: PtlNIInit succeeds (%d)\n", FUNCTION_NAME, ret);

  if ((ret = PtlEQAlloc(nih, 4096, NULL, &eqh)) != PTL_OK) {
    printf("%s: PtlEQAlloc failed: %d(%d)\n",
           FUNCTION_NAME, ret, child);
    exit(1);
  }
  printf("%s: PtlEQAlloc succeeds (%d:%d)\n", FUNCTION_NAME, child, ret);


  if ((ret = PtlGetId(nih, &rnk)) != PTL_OK) {
    printf("%s: PtlGetId failed: %d(%d)\n",
           FUNCTION_NAME, ret, child);
    exit(1);
  }
  printf("%s: nid=%d pid=%d(%d)\n", FUNCTION_NAME, rnk.nid, rnk.pid, child);

  if (child) {
    MPI_Init(&argc, &argv);
    MPI_Finalize();
    printf("%s: mpi_init and finalize succeed(%d)\n", FUNCTION_NAME, child);
  }
}
