#ifndef __tascel_UniformTaskCollectionShared_h__
#define __tascel_UniformTaskCollectionShared_h__

#include "Comm.h"
#include "FuncReg.h"
#include "SharedQueue.h"
#include "UniformTaskCollection.h"

namespace tascel {

  /**
   * Implementation of the UniformTaskCollection using a SharedQueue.
   *
   * This is a thin wrapper around the SharedQueue.
   */
  class UniformTaskCollectionShared : public UniformTaskCollection {
    private:
      SharedQueue sq; /**< the SharedQueue instance */

    public:
      /**
       * Constructs the UniformTaskCollectionShared.
       *
       * @copydetails UniformTaskCollection::UniformTaskCollection(const TaskCollProps&)
       */
      UniformTaskCollectionShared(const TaskCollProps &props);

      /**
       * Destroys the UniformTaskCollectionShared.
       */
      virtual ~UniformTaskCollectionShared();

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

  }; /*UniformTaskCollectionShared*/

}; /*tascel*/

#endif /*__tascel_UniformTaskCollectionShared_h__*/

