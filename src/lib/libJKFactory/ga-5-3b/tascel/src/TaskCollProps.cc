#if HAVE_CONFIG_H
# include "config.h"
#endif

#include "TaskCollProps.h"

using namespace tascel;

TaskCollProps::TaskCollProps()
  : max_tasks(0)
  , task_size(0)
  , frt()
  , tfn()
  , pldata(0)
  , pldata_len(0)
{
}

TaskCollProps& 
TaskCollProps::functions(
    const TslFunc &fn, const TslFuncRegTbl &frt) {
  this->frt = frt;
  this->tfn = fn;
  return *this;
}


TaskCollProps& 
TaskCollProps::taskSize(long tsk_size) {
  this->task_size = tsk_size;
  return *this;
}


TaskCollProps& 
TaskCollProps::maxTasks(long ntsks) {
  this->max_tasks = ntsks;
  return *this;
}


TaskCollProps& 
TaskCollProps::localData(void *pldata, long pldata_len) {
  this->pldata = static_cast<char*>(pldata);
  this->pldata_len = pldata_len;
  return *this;
}

