#ifndef METADATAPOLICIES_H
#define METADATAPOLICIES_H

#include "tensor.hpp"

namespace yeti {

class DoNothingMetaDataRetrieve {
    public:
        void retrieve(
            MetaDataNode* mdnode
        );
};

class RealignMemoryPoolMetaDataRetrieve {
    public:
        void retrieve(
            MetaDataNode* mdnode
        );
};

class RecomputeMetaDataRetrieve {
    public:
        void retrieve(
            MetaDataNode* mdnode
        );
};

class SortMetaDataRetrieve {
    public:
        void retrieve(
            MetaDataNode* mdnode
        );
};

class DoNothingMetaDataRelease {
    public:
        void release(
            MetaDataNode* mdnode
        );
};

}

#endif // METADATAPOLICIES_H
