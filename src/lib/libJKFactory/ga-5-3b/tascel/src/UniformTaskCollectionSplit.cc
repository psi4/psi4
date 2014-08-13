#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <cmath>
#include <cstdio>
#include <cstring>
#include <unistd.h>
#include <vector>

#include <armci.h>

#include "Comm.h"
#include "FuncReg.h"
#include "massert.h"
#include "Sleep.h"
#include "UniformTaskCollectionSplit.h"

using namespace std;
using namespace tascel;
using namespace tascel::comm;


UniformTaskCollectionSplit::UniformTaskCollectionSplit(
        const TaskCollProps &props)
  : UniformTaskCollection(props)
  , sq(tsk_size, max_ntsks) {
}

/*virtual */
UniformTaskCollectionSplit::~UniformTaskCollectionSplit() { }

void
UniformTaskCollectionSplit::process() {
  int p;
  bool got_work;
  barrier();
  TslFunc_t fn = frt.get(tfn);
  char buf[tsk_size];

  while (!sq.hasTerminated())  {
    while (sq.getTask(buf, tsk_size)) {
      stt.taskTime.startTimer();
      fn(this, buf, tsk_size, pldata, pldata_len, vector<void *>());
      stt.taskTime.stopTimer();
      stt.numTasks.inc();
    }
    stt.tdTime.startTimer();
    sq.td_progress();
    got_work = false;
    while(got_work == false && !sq.hasTerminated()) {
      do {
        p = rand() % nproc();
      }
      while (p == me());
      
      stt.stealTime.startTimer();
      got_work = sq.steal(p);
      stt.stealTime.stopTimer();
      stt.numStealAttempts.inc();
    }
    if(got_work) {
      stt.numSteals.inc();
    }
    
    if(sq.hasTerminated()) {
      stt.tdTime.stopTimer();
    }
  }
  barrier();
}

void
UniformTaskCollectionSplit::addTask(void *data, int dlen) {
  massert(dlen == tsk_size);
  sq.addTask(data, dlen);
  return;

error:
  throw TSL_ERR;
}

