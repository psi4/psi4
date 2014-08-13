#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <portals/portals3.h>
#include <portals/nal.h>
#include <mpi.h>

#define FORK_BEFORE_NI_INIT
#ifndef FORK_BEFORE_NI_INIT
#define FORK_AFTER_NI_INIT
#endif


#ifndef PMI_SUCCESS
#define PMI_SUCCESS 0
#endif
extern int PMI_CNOS_Get_nidpid_map(void **);
int
main(int argc, char **argv, char **envp)
{
  int i, ret, *npes;
  int num_interfaces;
  ptl_handle_ni_t nih;
  ptl_handle_eq_t eqh;
  ptl_ni_limits_t ptl_limits;
  pid_t child;
  ptl_process_id_t rnk, *procid_map;
  int spv, *spawned = &spv;


  if ((ret = PtlInit(&num_interfaces)) != PTL_OK) {
    printf("%s: PtlInit failed: %d\n", FUNCTION_NAME, ret);
    exit(1);
  }
  printf("%s: PtlInit succeeds (%d)\n", FUNCTION_NAME, ret);

#ifdef FORK_BEFORE_NI_INIT
  child = fork();
#endif

  if ((ret = PtlNIInit(IFACE_FROM_BRIDGE_AND_NALID(PTL_BRIDGE_UK, PTL_IFACE_SS),
                       PTL_PID_ANY, NULL, &ptl_limits, &nih)) != PTL_OK) {
    printf("%s: PtlNIInit failed: %d\n", FUNCTION_NAME, ret);
    /*exit(1);*/
  }
  else {
    printf("%s: PtlNIInit succeeds (%d)\n", FUNCTION_NAME, ret);
  }

#ifdef FORK_AFTER_NI_INIT
  child = fork();
#endif

  if ((ret = PtlEQAlloc(nih, 4096, NULL, &eqh)) != PTL_OK) {
    printf("%s: PtlEQAlloc failed: %d(%d)\n",
           FUNCTION_NAME, ret, child);
    exit(1);
  }
  printf("%s: PtlEQAlloc succeeds (%d:%d)\n", FUNCTION_NAME, child, ret);

#if 1
  if (child) {
    MPI_Init(&argc, &argv);
  }

  if (child) {
    PMI_Init(spawned);
    printf("\n%d:spanwned=%d", child, *spawned);
    if ((ret = PMI_Get_size(npes)) != PMI_SUCCESS) {
      printf("%s: PMI_Get_size failed: %d\n", FUNCTION_NAME, ret);
      /*exit(1);*/
    }
    else {
      printf("%s: PMI_Get_size succeeds (%d)\n", FUNCTION_NAME, *npes);
    }
    /*procid_map = (ptl_process_id_t *)malloc(sizeof(ptl_process_id_t)*(*npes));
    if(procid_map==NULL)exit(1);*/
    if ((ret = PMI_CNOS_Get_nidpid_map(&procid_map)) != PMI_SUCCESS) {
      printf("Getting proc map failed (npes=%d)\n", *npes);
    }
    for (i = 0; i < *npes; i++) {
      printf("\npid=%d nid=%d npes=%d(%d)", procid_map[i].pid, procid_map[i].nid, *npes, child);
    }
  }
#endif

  if ((ret = PtlGetId(nih, &rnk)) != PTL_OK) {
    printf("%s: PtlGetId failed: %d(%d)\n",
           FUNCTION_NAME, ret, child);
    exit(1);
  }
  printf("%s: nid=%d pid=%d(%d)\n", FUNCTION_NAME, rnk.nid, rnk.pid, child);
  if (child) {
    MPI_Finalize();
    printf("%s: mpi_init and finalize succeed(%d)\n", FUNCTION_NAME, child);
  }

}
