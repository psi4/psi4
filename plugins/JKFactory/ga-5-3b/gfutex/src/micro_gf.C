#include <ios>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <mpi.h>
#include <ga.h>
#include <macdecls.h>

#include "gfutex.h"

using namespace globalFutures;

static void GFTest(int g_a, int ndims, int lo[], int hi[], void *arg);

static const int numTimes = 20000;

static int g_proc;
static GFHandle gfhndl;

struct gftestarg {
};

int main(int argc, char *argv[])
{
  int me, nproc;
  double t0, t1, t2, t3;

  MPI_Init(&argc, &argv);
  GA_Initialize();

  GFInitialize();

  me = GA_Nodeid();
  nproc = GA_Nnodes();

  if (me == 0) {
    std::cout << "Number of procs: " << nproc << std::endl;
    std::cout << "Number of nodes: " << GA_Cluster_nnodes() << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  int dims[] = { nproc };
  int lo[] = { nproc - 1 }, hi[] = { nproc - 1 };

  g_proc = NGA_Create(MT_INT, 1, dims, "procs", NULL);

  GA_Print_distribution(g_proc);

  gfhndl = GFRegister(GFTest, sizeof(gftestarg));

  gftestarg gfa;

  t0 = MPI_Wtime();
  if (me == 0)
    for (int i = 0; i < numTimes; i++)
      GFEnqueue(gfhndl, g_proc, 1, lo, hi, &gfa);
  t1 = MPI_Wtime();

  MPI_Barrier(MPI_COMM_WORLD);

  t2 = MPI_Wtime();
  GFAllQuiesce(gfhndl);
  t3 = MPI_Wtime();

  std::cout << "Proc: " << me << ", enqueue exec. time: " << t1 - t0 << ", quiesce time: " <<
    t3 - t2 << std::flush << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);

  t0 = MPI_Wtime();
  if (me == 0)
    for (int i = 0; i < numTimes; i++)
      GFExecute(gfhndl, g_proc, 1, lo, hi, &gfa);
  t1 = MPI_Wtime();

  MPI_Barrier(MPI_COMM_WORLD);

  t2 = MPI_Wtime();
  GFAllQuiesce(gfhndl);
  t3 = MPI_Wtime();

  std::cout << "Proc: " << me << ", execute exec. time: " << t1 - t0 << ", quiesce time: " <<
    t3 - t2 << std::flush << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);

  GA_Destroy(g_proc);

  GFFinalize();

  GA_Terminate();
  MPI_Finalize();

  return 0;
} // main

void GFTest(int g_a, int ndims, int lo[], int hi[], void *arg)
{
  gftestarg *gfa = static_cast<gftestarg *>(arg);
} // GFTest
