#ifndef yeti_tensorblock_hpp
#define yeti_tensorblock_hpp

namespace yeti {


class TensorBlock;
class TensorController;

typedef boost::intrusive_ptr<TensorBlock> TensorBlockPtr;

typedef BlockMap<TensorBlock> TensorBlockMap;

/**
    This specifies a location, relative to the beginning of a memory
    pool, for metadata to beginning reading.
*/
typedef size_t tensor_block_vir_addr_unsigned_offset_t;

}

#endif

