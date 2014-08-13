#ifndef __tascel_UniformTaskCollection_h__
#define __tascel_UniformTaskCollection_h__

#include "FuncReg.h"
#include "TaskCollProps.h"
#include "StealingStats.h"

namespace tascel {

  /**
   * Abstract base class representing a collection of uniform tasks.
   *
   * Here, uniform means that each task descriptor is the same size.  The
   * function table is allowed to have as many functions as needed, however
   * the entry function must be indicated.
   */
  class UniformTaskCollection {
    protected:
      const int            max_ntsks;  /**< max number of tasks in collection */
      const int            tsk_size;   /**< size of any one task */
      const TslFuncRegTbl &frt;        /**< function registration table */
      TslFunc              tfn;        /**< task function */
      char                *pldata;     /**< Process local data shared across tasks*/
      const int            pldata_len; /**< length of pldata buffer */
      StealingStats        stt;        /**< Stealing statistics*/
    public:
      /**
       * Constructs the task collection.
       *
       * @param[in] _fn entry function representing the task
       * @param[in] frt function registration table of which _fn is a member
       * @param[in] tsk_size size of any one task
       * @param[in] ntsks maximum number of tasks allowed in collection
       * @param[in] pldata Pointer to process-local data
       * @param[in] pldata_len Length of buffer pointed by pldata
       *
       * @pre _fn is a valid handle to a function in frt
       */
      UniformTaskCollection(const TaskCollProps &props);

      /**
       * Destroys the task collection.
       */
      virtual ~UniformTaskCollection();

      /**
       * Begins processing the tasks in this collection.
       *
       * This is a collective operation.
       */
      virtual void process() = 0;

      /**
       * Adds a task to this collection.
       *
       * @param[in] data task descriptor
       * @param[in] dlen size of the task descriptor
       *
       * @pre dlen == tsk_size
       */
      virtual void addTask(void *data, int dlen) = 0;

      /**
       * Print stealing statistics
       */
      virtual void printStats() const;

  }; /*UniformTaskCollection*/

}; /*tascel*/

#endif /*__tascel_UniformTaskCollection_h__*/

