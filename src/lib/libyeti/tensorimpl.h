#ifndef TENSORIMPL_H
#define TENSORIMPL_H

#include "tensor.h"

namespace yeti {

template <
    class ValidatePolicy,
    class MempoolRetrievePolicy,
    class MetaDataBranchRetrievePolicy,
    class TensorDataControllerRetrievePolicy,
    class DataBranchRetrievePolicy,
    class BranchFlushPolicy,
    class BranchReleasePolicy,
    class MetaDataNodeRetrievePolicy,
    class DataNodeRetrievePolicy,
    class DataStorageFlushPolicy,
    class ObsoletePolicy,
    class SyncPolicy,
    class DataClearPolicy
>
class TensorControllerTemplate :
    public TensorController,
    /** Policies */
    public ValidatePolicy,
    public MempoolRetrievePolicy,
    public MetaDataBranchRetrievePolicy,
    public TensorDataControllerRetrievePolicy,
    public DataBranchRetrievePolicy,
    public BranchFlushPolicy,
    public BranchReleasePolicy,
    public MetaDataNodeRetrievePolicy,
    public DataNodeRetrievePolicy,
    public DataStorageFlushPolicy,
    public ObsoletePolicy,
    public SyncPolicy,
    public DataClearPolicy
{

    public:
        TensorControllerTemplate(TensorBlock* block)
            : TensorController(block)
        {
        }

        void clear()
        {
            DataClearPolicy::clear(parent_block_);
        }

        void validate()
        {
            ValidatePolicy::validate(parent_block_);
        }

        void retrieve(
            TensorBranch* branch
        )
        {
            branch->set_parent(parent_block_);
            MempoolRetrievePolicy
                ::retrieve(branch);
            MetaDataBranchRetrievePolicy
                ::retrieve(branch);
            TensorDataControllerRetrievePolicy
                ::retrieve(branch);
            DataBranchRetrievePolicy
                ::retrieve(branch);
        }

        void flush(
            TensorBranch* branch
        )
        {
            BranchFlushPolicy
                ::flush(branch);
        }
        
        void release(
            TensorBranch* branch
        )
        {
            BranchReleasePolicy
                ::release(branch);
        }

        void retrieve(MetaDataNode* mdnode)
        {
            MetaDataNodeRetrievePolicy
                ::retrieve(mdnode);
        }

        bool need_recompute_data()
        {
            return DataNodeRetrievePolicy::need_recompute();
        }

        void retrieve(
            DataNode* dnode,
            const MetaDataNode* parent,
            const uli* indices
        )
        {
            DataNodeRetrievePolicy
                ::retrieve(dnode, parent, indices);
        }

        void flush_data()
        {
            DataStorageFlushPolicy
                ::flush(parent_block_);
        }

        void obsolete()
        {
            ObsoletePolicy
                ::obsolete(parent_block_);
        }

        void sync()
        {
            SyncPolicy::sync(parent_block_);
        }


};


}

#endif // TENSORIMPL_H
