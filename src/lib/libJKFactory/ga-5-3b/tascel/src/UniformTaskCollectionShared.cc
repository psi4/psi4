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
#include "UniformTaskCollectionShared.h"

using namespace std;
using namespace tascel;
using namespace tascel::comm;


UniformTaskCollectionShared::UniformTaskCollectionShared(
        const TaskCollProps &props)
  : UniformTaskCollection(props)
  , sq(tsk_size, max_ntsks) {
}

/* virtual */
UniformTaskCollectionShared::~UniformTaskCollectionShared() { }

void
UniformTaskCollectionShared::process() {
  int p;
  bool got_work;
  barrier();
  TslFunc_t fn = frt.get(tfn);
  char buf[tsk_size];
  while (!sq.hasTerminated())  {
    while (sq.getTask(buf, tsk_size)) {
      fn(this, buf, tsk_size, NULL, 0, vector<void *>());
    }
    sq.td_progress();
    got_work = false;
    while(got_work == false && !sq.hasTerminated()) {
      do {
        p = rand() % nproc();
      }
      while (p == me());
      got_work = sq.steal(p);
    }
  }
  barrier();
}

void
UniformTaskCollectionShared::addTask(void *data, int dlen) {
  massert(dlen == tsk_size);
  sq.addTask(data, dlen);
  return;

error:
  throw TSL_ERR;
}

