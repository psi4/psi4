#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstdlib>

#include <mpi.h>

#include "Comm.h"

using namespace std;
using namespace tascel;
using namespace tascel::comm;


void
comm::barrier() {
  MPI_Barrier(MPI_COMM_WORLD);
}

int
comm::nproc() {
  int rval;
  MPI_Comm_size(MPI_COMM_WORLD, &rval);
  return rval;
}

int
comm::me() {
  int rval;
  MPI_Comm_rank(MPI_COMM_WORLD, &rval);
  return rval;
}

