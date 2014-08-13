#ifndef __tascel_UniformTaskCollSplitData_h__
#define __tascel_UniformTaskCollSplitData_h__

#include "AccessMode.h"
#include "Comm.h"
#include "FuncReg.h"
#include "SplitQueueOpt.h"
#include "TaskCollProps.h"
#include "UniformTaskCollectionSplit.h"

namespace tascel {

  class DataColl;
  
  /**
   * TODO
   */
  typedef void (*compute_index_t)(void *dscr, int dscr_len, 
				 void *pldata, int pldata_len, 
				 int collid, void *idx, int idxlen);

  /**
   * Implementation of the UniformTaskCollection using a SplitQueue and data
   * movement.
   */
  class UniformTaskCollSplitData : public UniformTaskCollectionSplit {
    private:
      std::vector<DataColl*> colls;         /**< TODO */
      const std::vector<AccessMode> modes;  /**< TODO */
      const std::vector<int> idxlens;       /**< TODO */
      compute_index_t ci_fn;                /**< TODO */

      /**
       * TODO
       *
       * @param[in,out] data_bufs TODO
       * @param[in,out] data_lens TODO
       * @param[in] idxs TODO
       */
      void setupDataBufs(std::vector<void *> &data_bufs,
                         std::vector<int> &data_lens,
                         const std::vector<void *> &idxs);

      /**
       * TODO
       *
       * @param[in,out] data_bufs TODO
       * @param[in,out] data_lens TODO
       * @param[in] idxs TODO
       */
      void cleanupDataBufs(std::vector<void *> &data_bufs, 
			               std::vector<int> &data_lens, 
			               const std::vector<void *> &idxs);

    public:
      /**
       * Constructs the UniformTaskCollSplitData.
       *
       * @copydetails UniformTaskCollection::UniformTaskCollection(TslFunc,const TslFuncRegTbl&,int,int,void*,int)
       *
       * @param[in] colls TODO
       * @param[in] modes TODO
       * @param[in] idxlens TODO
       * @param[in] ci_fn TODO
       */
      UniformTaskCollSplitData(const TaskCollProps &props, 
			                   const std::vector<DataColl*>& colls, 
			                   const std::vector<AccessMode>& modes, 
			                   const std::vector<int>& idxlens,
			                   compute_index_t ci_fn);

      /**
       * Destroys the UniformTaskCollSplitData.
       */
      virtual ~UniformTaskCollSplitData();

      /**
       * @copybrief   UniformTaskCollection::process()
       * @copydetails UniformTaskCollection::process()
       */
      void process();

  }; /*UniformTaskCollSplitData*/
}; /*tascel*/

#endif /*__tascel_UniformTaskCollSplitData_h__*/

