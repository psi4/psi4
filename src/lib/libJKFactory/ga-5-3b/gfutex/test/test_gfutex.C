#include <ios>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/types.h>

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <ga.h>
#include <macdecls.h>
#include <gpc.h>

#include <tbb/atomic.h>
#include <tbb/enumerable_thread_specific.h>

#include "gfutex.h"

using namespace tbb;
using namespace globalFutures;

typedef enumerable_thread_specific<int> TInfo;

static void GFTest(int g_a, int ndims, int lo[], int hi[], void *arg);

static const int numTimes = 2000;

static int g_a;
static GFHandle gfhndl;
static atomic<int> sum;
static int esum = 0;
static int *procs_contr = NULL;

static TInfo tinfo(0);

struct gftestarg {
  int val;
  int src;
};

int main(int argc, char *argv[])
{
  int me, nproc;
  int firstNLP;
  double t0, t1, t2, t3;

  MPI_Init(&argc, &argv);
  GA_Initialize();

  GFInitialize();

  me = GA_Nodeid();
  nproc = GA_Nnodes();

  procs_contr = new int[nproc];

  for (int i = 0; i < nproc; i++)
    procs_contr[i] = 0;

  firstNLP = nproc / GA_Cluster_nnodes();

  if (me == 0) {
    std::cout << "Number of procs: " << nproc << std::endl;
    std::cout << "Number of nodes: " << GA_Cluster_nnodes() << std::endl;
    std::cout << "Number of procs per node: " << firstNLP << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

#if 0
  for (int i = 0; i < nproc; i++) {
    if (me == i)
      std:: cout << "Proc: " << me << ", nodeid: " << GA_Cluster_nodeid() <<
	std::flush << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  sum = 0;

  for (int i = 0; i < nproc; i++)
    esum += (i + 1);

  int dims[] = { 1000, 1000 };
  int lo[] = { 10, 10 }, hi[] = { 10, 10 };

  g_a = NGA_Create(MT_DBL, 2, dims, "A", NULL);

  GA_Print_distribution(g_a);

  gfhndl = GFRegister(GFTest, sizeof(gftestarg));

  gftestarg gfa;

  gfa.val = me + 1;
  gfa.src = me;

#if 1
  t0 = MPI_Wtime();
  if (me >= 0)
    for (int i = 0; i < numTimes; i++)
      GFExecute(gfhndl, g_a, 2, lo, hi, &gfa);
  t1 = MPI_Wtime();

  MPI_Barrier(MPI_COMM_WORLD);

  t2 = MPI_Wtime();
  GFAllQuiesce(gfhndl);
  t3 = MPI_Wtime();

  if (sum != 0) {
    std::cout << "Proc: " << me << ", sum: " << std::fixed << sum << ", expected sum: " <<
      numTimes * esum << std::flush << std::endl;
    std::cout << "Proc: " << me << ", exec. time: " << t1 - t0 << ", quiesce time: " <<
      t3 - t2 << std::flush << std::endl;
  }
#endif

#if 0
  std::cout << "Proc: " << me << ", added: " << std::fixed << numTimes * gfa.val <<
    std::flush << std::endl;
#endif

  if (sum != 0)
    for (int i = 0; i < nproc; i++)
      std::cout << "Proc: " << me << ", procs_contr[" << std::dec << i << "]: " <<
	procs_contr[i] << std::flush << std::endl;

  std::cout << std::endl;

  for (TInfo::const_iterator iter = tinfo.begin(); iter != tinfo.end(); iter++)
    std::cout << *iter << " ";

  std::cout << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);

  GA_Destroy(g_a);

  GFFinalize();

  GA_Terminate();
  MPI_Finalize();

  delete [] procs_contr;

  return 0;
} // main

void GFTest(int g_a, int ndims, int lo[], int hi[], void *arg)
{
  gftestarg *gfa = static_cast<gftestarg *>(arg);

  sum += gfa->val;

  procs_contr[gfa->src] = 1;

  TInfo::reference ltinfo = tinfo.local();

  ltinfo++;
} // GFTest
