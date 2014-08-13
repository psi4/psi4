#ifndef __tascel_UniformTaskCollectionSplit_h__
#define __tascel_UniformTaskCollectionSplit_h__

#include "Comm.h"
#include "FuncReg.h"
#include "SplitQueueOpt.h"
#include "TaskCollProps.h"
#include "UniformTaskCollection.h"

namespace tascel {

  /**
   * Implementation of the UniformTaskCollection using a SplitQueue.
   *
   * This is a thin wrapper around the SplitQueue.
   */
  class UniformTaskCollectionSplit : public UniformTaskCollection {
    protected:
      SplitQueueOpt sq; /**< the SplitQueue instance */

    public:
      /**
       * Constructs the UniformTaskCollectionSplit.
       *
       * @copydetails UniformTaskCollection::UniformTaskCollection(const TaskCollProps&)
       */
      UniformTaskCollectionSplit(const TaskCollProps& props);

      /**
       * Destroys the UniformTaskCollectionSplit.
       */
      virtual ~UniformTaskCollectionSplit();

      /**
       * @copybrief   UniformTaskCollection::process()
       * @copydetails UniformTaskCollection::process()
       */
      virtual void process();

      /**
       * @copybrief   UniformTaskCollection::addTask()
       * @copydetails UniformTaskCollection::addTask()
       */
      virtual void addTask(void *data, int dlen);

  }; /*UniformTaskCollectionSplit*/

}; /*tascel*/

#endif /*__tascel_UniformTaskCollectionSplit_h__*/

