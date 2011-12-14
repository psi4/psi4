#ifndef yeti_cache_h
#define yeti_cache_h

#include "class.h"
#include "yetiobject.h"

#include "cache.hpp"
#include "data.hpp"
#include "tensor.hpp"
#include "thread.hpp"
#include "mallocimpl.h"
#include "tuple.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define MALLOC_CACHE_BLOCK 0

#define DEBUG_CACHE_USAGE 0

#define CACHE_BLOCK_PULL_ERROR 101

#define DEBUG_CACHE_HISTORY 0

#define NHISTORIES_CACHE 100

namespace yeti {


class Cachable :
    public YetiRuntimeCountable
{
    protected:
        bool in_destructor_;

    public:
        Cachable();

        Cachable(YetiRuntimeObject::thread_safety_flag_t flag);

        virtual void flush_from_cache() = 0;

        bool in_destructor() const;


};

struct CacheHistory {
    uli thread_number;
    Cachable* owner;
};

/** @ingroup MemoryManagement
    @class DataCacheEntry
    A single entry in a data cache which
    owns a particular block of memory of given size.
*/
struct DataCacheEntry :
    public Malloc<DataCacheEntry>
{

    public:
#if DEBUG_CACHE_HISTORY
        uli nentries_history_;
        
        CacheHistory histories_[NHISTORIES_CACHE];
#endif

        char* data;

        /** This holds a pointer to the data location */
        Cachable* owner;

        bool pulled;

        DataCacheEntry* next;

        DataCacheEntry* prev;

        uli threadnum;

    public:
        /**
            Create a data cache entry pointing at the given memory location
            @param d The memory location
        */
        DataCacheEntry(char* d);

        ~DataCacheEntry();

        void print(std::ostream& os) const;

        void* operator new(size_t size, MemoryPool* mempool);
        
        void operator delete(void* ptr, MemoryPool* mempool);

        void* operator new(size_t size);

        void operator delete(void* ptr);


};

#define NENTRIES_CACHE_MAX 10000

#define NHISTORIES_ALL 100000

/** @ingroup MemoryManagement
    @class DataCache
    A cache holding a set of memory blocks of uniform size.  This should not allocate
    its own memory, but receive the allocation from a LayeredDataCache.
*/
class DataCache :
    public smartptr::Countable,
    public MempoolVirtualAddressMalloc,
    public Malloc<DataCache>
{

    protected:
#if DEBUG_CACHE_HISTORY
        CacheHistory all_histories_[NHISTORIES_ALL];

        uli nhistories_;
#endif

        DataCacheEntry* all_entries_[NENTRIES_CACHE_MAX];

        DataCacheEntry* free_entries_start_;

        DataCacheEntry* used_entries_start_;

        DataCacheEntry* used_entries_end_;

        size_t blocksize_;

        uli nentries_;
        
        /**
            Start looking for an open cache entry at the offset
        */
        DataCacheEntry* find_free_entry();

        DataCacheEntry* find_used_entry();

        DataCacheEntry* allocate_new_entry();

        DataCache* cache_extension_;

        char* extension_block_;

        DataCache* parent_;

        DataCache(
            DataCache* parent,
            MemoryPool* mempool,
            uli nelements,
            size_t blocksize
        );

        DataCacheEntry* pull_no_lock(Cachable* owner);

        void append_used_entries(DataCacheEntry* head_entry);

        void append_free_entries(DataCacheEntry* head_entry);

        MultiThreadLock* lock_;

        std::string name_;

    public:
        /**
            @param Total size in bytes of available cache storage
            @param size_counts Key data size, value number of tiles with that size.
                               The data sizes will be used to determine which data
                               cache sizes are needed
        */
        DataCache(
            MultiThreadLock* lock,
            size_t blocksize,
            const std::string& name
        );

        ~DataCache();

        const std::string& name() const;

        DataCacheEntry* get_first_free_entry() const;

        DataCacheEntry* get_first_used_entry() const;

        DataCacheEntry* get_last_used_entry() const;

        void build_extension(uli nelements);

#if 0
        void delete_extension();
#endif

        void flush();

        /**
            @param entry
        */
        void free(DataCacheEntry* entry, Cachable* owner);

        void insert(DataCacheEntry* entry);

        /**
            This method is thread-safe and returns a locked cache entry
            @param entry
        */
        DataCacheEntry* pull(DataCacheEntry* entry, Cachable* owner);

        /**
            Pull a new data cache entry and assign the given data block
            as an owner. The cache block cannot be cleared while it is
            pulled and there assumes the block has been locked before being called.
            @param block The block to assign as owner to the new cache entry
            @return A new data cache entry
        */
        DataCacheEntry* pull(Cachable* owner);

        /**
            @return The number of blocks with no linked data
        */
        uli nfree() const;

        void lock();

        void unlock();

        /**
            @return The total number of blocks in the cache
        */
        uli ntotal() const;

        /**
            All cache blocks are of uniform size.
            @return The size of each block in the cache
        */
        uli blocksize() const;

#if DEBUG_CACHE_HISTORY
        void print_history();
#endif

        void print(std::ostream& os = std::cout) const;

        void* operator new(size_t size, MemoryPool* mempool);
        
        void operator delete(void* ptr, MemoryPool* mempool);

        void* operator new(size_t size);

        void operator delete(void* ptr);

};


}

#ifdef redefine_size_t
#undef size_t
#endif

#endif
