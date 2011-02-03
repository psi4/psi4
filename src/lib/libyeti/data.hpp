#ifndef yeti_data_hpp
#define yeti_data_hpp


namespace yeti {

class Data;
class DataBlock;
class DataMode;
class Buffer;
class DiskBuffer;
class CachedDataBlock;
class DataBlockFactory;
class SortedAccumulateBlock;

typedef boost::intrusive_ptr<DataBlock> DataBlockPtr;
typedef boost::intrusive_ptr<Buffer> BufferPtr;
typedef boost::intrusive_ptr<DiskBuffer> DiskBufferPtr;
typedef boost::intrusive_ptr<CachedDataBlock> CachedDataBlockPtr;
typedef boost::intrusive_ptr<DataBlockFactory> DataBlockFactoryPtr;

}

#endif

