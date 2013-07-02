#ifndef yeti_data_hpp
#define yeti_data_hpp


namespace yeti {

class StorageBlock;
class InCoreBlock;
class DiskBuffer;
class CachedStorageBlock;
class LocalDiskBlock;

typedef boost::intrusive_ptr<StorageBlock> StorageBlockPtr;
typedef boost::intrusive_ptr<DiskBuffer> DiskBufferPtr;
typedef boost::intrusive_ptr<CachedStorageBlock> CachedDataBlockPtr;


}

#endif

