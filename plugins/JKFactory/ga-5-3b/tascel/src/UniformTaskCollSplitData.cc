#include "Comm.h"
#include "DataColl.h"
#include "FuncReg.h"
#include "UniformTaskCollSplitData.h"
#include "massert.h"

#include <armci.h>

#include <unistd.h>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

using namespace std;
using namespace tascel;
using namespace comm;

UniformTaskCollSplitData::UniformTaskCollSplitData(const TaskCollProps &props,
    const std::vector<DataColl*>& _colls, 
    const std::vector<AccessMode>& _modes, 
    const std::vector<int>& _idxlens,
    compute_index_t _ci_fn)
  : UniformTaskCollectionSplit(props)
  , ci_fn(_ci_fn)
  , colls(_colls.begin(), _colls.end())
  , modes(_modes.begin(), _modes.end())
  , idxlens(_idxlens.begin(), _idxlens.end())
{
  assert(colls.size() == modes.size());
  assert(colls.size() == idxlens.size());
}

/*virtual */
UniformTaskCollSplitData::~UniformTaskCollSplitData() { }

void
UniformTaskCollSplitData::setupDataBufs(vector<void *> &data_bufs, vector<int> &data_lens, const vector<void *> &idxs) {
  for(int i=0; i<idxlens.size(); i++) {
    int sz = colls[i]->getSize(idxs[i],idxlens[i]);
    char *dbuf = new char[sz];
    data_bufs.push_back(dbuf);
    data_lens.push_back(sz);

    switch(modes[i]) {
    case MODE_RONLY:
      colls[i]->get(idxs[i], idxlens[i], dbuf, sz);
      break;
    case MODE_RDWR:  
      colls[i]->get(idxs[i], idxlens[i], dbuf, sz);
      break;
    case MODE_ACC:
      memset(dbuf, 0, sz);
      break;
    }
  }
}

void
UniformTaskCollSplitData::cleanupDataBufs(vector<void *> &data_bufs, 
            vector<int> &data_lens, 
            const vector<void *> &idxs) {
  for(int i=0; i<idxlens.size(); i++) {
    int sz = data_lens[i];
    char *dbuf = (char *)data_bufs[i];

    switch(modes[i]) {
    case MODE_RONLY:
      break;
    case MODE_RDWR:  
      colls[i]->put(idxs[i], idxlens[i], dbuf, sz);
      break;
    case MODE_ACC:
      colls[i]->add(idxs[i], idxlens[i], dbuf, sz);
      break;
    }
    delete [] dbuf;
  }
  data_bufs.clear();
  data_lens.clear();
}

void
UniformTaskCollSplitData::process() {
  int p;
  bool got_work;
  barrier();
  TslFunc_t fn = frt.get(tfn);
  char buf[tsk_size];
  int tasksDone=0, stealAttempts=0, steals=0;
  vector<void *> idxs;
  for(int i=0; i<idxlens.size(); i++) {
    idxs.push_back(new char[idxlens[i]]);
  }

  while (!sq.hasTerminated())  {
    while (sq.getTask(buf, tsk_size)) {
      tasksDone += 1;
      vector<int> data_lens;
      vector<void *> data_bufs;
      for(int i=0; i<idxlens.size(); i++) {
        ci_fn(buf, tsk_size, pldata, pldata_len, 
        i, idxs[i], idxlens[i]);
      }
      setupDataBufs(data_bufs, data_lens, idxs);
      fn(this, buf, tsk_size, pldata, pldata_len, data_bufs);
      cleanupDataBufs(data_bufs, data_lens, idxs);
    }
    sq.td_progress();
    got_work = false;
    while(got_work == false && !sq.hasTerminated()) {
      do {
        p = rand() % nproc();
      }
      while (p == me());
      got_work = sq.steal(p);
      stealAttempts += 1;
    }
    if(got_work) {
      steals += 1;
    }
  }

//   printf("%d: processed %d tasks attempts=%d steals=%d\n", me(), tasksDone, stealAttempts, steals);
  barrier();
}



