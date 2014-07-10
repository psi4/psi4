#ifndef __tascel_TaskCollProps_h__
#define __tascel_TaskCollProps_h__

#include "FuncReg.h"
#include "TaskCollProps.h"

namespace tascel {

  /**
   * TODO
   */
  class TaskCollProps {
    public:
      TaskCollProps();

      TaskCollProps& functions(const TslFunc &fn, const TslFuncRegTbl &frt);
      TaskCollProps& taskSize(long tsk_size);
      TaskCollProps& maxTasks(long ntsks);
      TaskCollProps& localData(void *pldata, long pldata_len);

    private:
      friend class UniformTaskCollection;
      friend class UniformTaskCollectionShared;
      friend class UniformTaskCollectionSplit;
      friend class UniformTaskCollSplitData;

      long max_tasks;
      long task_size;
      TslFuncRegTbl frt;
      TslFunc tfn;
      char *pldata;
      long pldata_len;
  };

}; /* tascel */

#endif /* __tascel_TaskCollProps_h__ */
