#ifndef DATAPOLICIES_H
#define DATAPOLICIES_H

#include "tensor.hpp"

namespace yeti {

class DoNothingDataRetrieve {
    public:
        void retrieve(
            DataNode* node,
            const MetaDataNode* parent,
            const uli* indices
        );

        bool need_recompute() const;
};

class RecomputeDataRetrieve {
    private:
        template <typename data_t>
        void compute(
            TensorElementComputer* filler,
            DataNode* node,
            const uli* indices
        );
        
    public:
        void retrieve(
            DataNode* node,
            const MetaDataNode* parent,
            const uli* indices
        );

        bool need_recompute() const;
};

class DoNothingDataRelease {
    public:
        void release(
            DataNode* node,
            const MetaDataNode* parent
        );
};

class DoNothingDataStorageFlush {
    public:
        void flush(TensorBlock* block);
};

class DeleteAllDataStorageFlush {
    public:
        void flush(TensorBlock* block);
};

class CommitAllDataStorageFlush {
    public:
        void flush(TensorBlock* block);
};

class ClearAllDataStorageFlush {
    public:
        void flush(TensorBlock* block);
};

class DeleteAllDataClear {
    public:
        void clear(
            TensorBlock* block
        );
};

class ClearAllData {
    public:
        void clear(
            TensorBlock* block
        );
};

}

#endif // DATAPOLICIES_H
