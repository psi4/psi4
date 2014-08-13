#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstdio>
#include <cstring>
#include <vector>

#include <armci.h>

#include "Comm.h"
#include "FuncReg.h"
#include "massert.h"
#include "TaskCollProps.h"
#include "UniformTaskCollection.h"

using namespace std;
using namespace tascel;
using namespace comm;


UniformTaskCollection::UniformTaskCollection(const TaskCollProps &props)
  : max_ntsks(props.max_tasks)
  , tsk_size(props.task_size)
  , frt(props.frt)
  , tfn(props.tfn)
  , pldata(new char[props.pldata_len])
  , pldata_len(props.pldata_len) {
  memcpy(pldata, props.pldata, pldata_len);
}

/* virtual */
UniformTaskCollection::~UniformTaskCollection() {
  delete [] pldata;
}

/* virtual */
void UniformTaskCollection::printStats() const {
  stt.print();
}

