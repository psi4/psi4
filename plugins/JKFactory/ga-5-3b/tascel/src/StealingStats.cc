#include "StealingStats.h"
#include "Comm.h"
#include <mpi.h>
#include <cstdio>

using namespace tascel;
using namespace tascel::comm;

void
StealingStats::print() const {
  const int np = nproc();
  const int myid = me();

  if(myid==0) {
    printf(" pid #stealattempts #steals #tasks Steal-time Task-time td-time\n");
    printf("================================================================\n"); 
  }

  for(int i=0; i<np; i++) {
    if(myid == i) {
      printf("%4d %8d %8d %6d %10.2lf %10.2lf %10.2lf\n",
	     myid, 
	     numStealAttempts.count(),
	     numSteals.count(),
	     numTasks.count(),
	     stealTime.time(),
	     taskTime.time(),
	     tdTime.time());
      fflush(stdout);
    }
    barrier();
  }  
}

