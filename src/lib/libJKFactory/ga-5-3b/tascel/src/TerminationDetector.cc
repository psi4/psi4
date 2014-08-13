#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <mpi.h>

#include "Comm.h"
#include "TerminationDetector.h"

using namespace std;
using namespace tascel;
using namespace tascel::comm;


void
TerminationDetector::resetUpPhase() {
  if (child[0] >= 0) {
    MPI_Irecv(&cvote[0], 1, MPI_INT, child[0], child[0], MPI_COMM_WORLD, &creq[0]);
  }
  if (child[1] >= 0) {
    MPI_Irecv(&cvote[1], 1, MPI_INT, child[1], child[1], MPI_COMM_WORLD, &creq[1]);
  }
  myvote = 1; //yes
  dirUp = false;
}

void
TerminationDetector::resetDownPhase() {
  if (parent >= 0) {
    MPI_Irecv(&pdecision, 1, MPI_INT, parent, parent, MPI_COMM_WORLD, &preq);
  }
  dirUp = true;
}

TerminationDetector::TerminationDetector() {
  child[0] = (2 * me() + 1 < nproc()) ? 2 * me() + 1 : -1;
  child[1] = (2 * me() + 2 < nproc()) ? 2 * me() + 2 : -1;
  parent = (me != 0) ? (me() - 1) / 2 : -1;
  pdecision = 0;
  resetUpPhase();
  resetDownPhase(); //this should after up. we start with the up phase
}

TerminationDetector::~TerminationDetector() {
  if (child[0] >= 0) {
    MPI_Cancel(&creq[0]);
  }
  if (child[1] >= 0) {
    MPI_Cancel(&creq[1]);
  }
  if (parent >= 0) {
    MPI_Cancel(&preq);
  }
}

void
TerminationDetector::vote(bool _myvote) {
  myvote = (_myvote == false) ? 0 : myvote;
}

void
TerminationDetector::propogateUp() {
  int decision = myvote;
  if (child[0] >= 0 && cvote[0] == 0) {
    decision = 0;
  }
  if (child[1] >= 0 && cvote[1] == 0) {
    decision = 0;
  }
  if (me() == 0) {
    pdecision = decision;
  }

  if (parent >= 0) {
    MPI_Send(&decision, 1, MPI_INT, parent, me(), MPI_COMM_WORLD);
  }
  resetUpPhase();
}

void
TerminationDetector::propogateDown() {
  if (child[0] >= 0) {
    MPI_Send(&pdecision, 1, MPI_INT, child[0], me(), MPI_COMM_WORLD);
  }
  if (child[1] >= 0) {
    MPI_Send(&pdecision, 1, MPI_INT, child[1], me(), MPI_COMM_WORLD);
  }
  resetDownPhase();
}

void
TerminationDetector::progress(bool _myvote) {
  vote(_myvote);
  if (pdecision == 1 || (me() == 0 && nproc() == 0 && myvote == 1)) {
    return; //termination detected. nothing more to do
  }
  if (dirUp) {
    int flag[2] = {1, 1};
    MPI_Status status[2];
    if (child[0] >= 0) {
      MPI_Test(&creq[0], &flag[0], &status[0]);
    }
    if (child[1] >= 1) {
      MPI_Test(&creq[1], &flag[1], &status[1]);
    }
    if (flag[0] && flag[1]) {
      propogateUp();
    }
  }
  /*Tiny optimization: propogate inverts dirUp. So both dirUp and
    !dirUp could be executed consecutively, especially for me()==0*/
  if (!dirUp) {
    int flag = 1;
    MPI_Status status;
    if (parent >= 0) {
      MPI_Test(&preq, &flag, &status);
    }
    if (flag) {
      propogateDown();
    }
  }
}

bool
TerminationDetector::hasTerminated() {
//   progress(myvote);
  return pdecision == 1;
}

