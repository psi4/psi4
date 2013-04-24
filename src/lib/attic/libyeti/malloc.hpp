#ifndef yeti_malloc_hpp
#define yeti_malloc_hpp



namespace yeti {

class MemoryPool;

/// @ingroup MemoryManagement
typedef boost::intrusive_ptr<MemoryPool> MemoryPoolPtr;

}

#endif

