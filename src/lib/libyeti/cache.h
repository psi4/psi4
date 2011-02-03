#ifndef yeti_cache_h
#define yeti_cache_h

#include "class.h"
#include "yetiobject.h"

#include "cache.hpp"
#include "data.hpp"
#include "tensor.hpp"
#include "thread.hpp"
#include "mallocimpl.h"

namespace yeti {

/**
    @class DataCacheEntry
    A single entry in a data cache which owns a particular block of memory of given size.
*/
struct DataCacheEntry :
    public Malloc<DataCacheEntry>,
    public YetiRuntimeObject
{

    public:
        DataCacheEntry* next;

        DataCacheEntry* prev;

        char* data;

        /** This holds a pointer to the data location */
        CachedDataBlock* owner;

    public:
        /**
            Create a data cache entry pointing at the given memory location
            @param d The memory location
        */
        DataCacheEntry(char* d);

        ~DataCacheEntry();

        /**
            Clear the data block that might currently own the memory location
        */
        void clear();

        /**
            Assign this data cache entry to a new data block owner
            @param owner
        */
        void reassign(CachedDataBlock* owner);

        void print(std::ostream& os) const;

        void finalize_accumulate();

        void finalize_write();

};

/**
    @class DataCache
    A cache holding a set of memory blocks of uniform size.  This should not allocate
    its own memory, but receive the allocation from a LayeredDataCache.
*/
class DataCache :
    public YetiThreadedRuntimeObject
{

    private:
        friend class LayeredDataCache;

        /**
            The front of the linked list of cache blocks.  This will be the last
            block to popped off the queue for use.
        */
        DataCacheEntry* first_;

        /**
            The back of the linked list of cache blocks.  This will the first block
            to be popped off the queue for use.
        */
        DataCacheEntry* last_;

        uli blocksize_;

        uli n_;

        template <class T>
        void
        cache_action();

    public:
        /**
            Each data cache must have at least one entry for use.
            @param blocksize The uniform size of all blocks in the data cache
            @param datastart A block of data for use in initializing the first cache entry
        */
        DataCache(
            size_t blocksize,
            char* datastart
        );

        ~DataCache();

        /**
            Put the entry on the front of the list.  This is
            then the last block scheduled for reuse.
            This method is thread-safe.
            @param entry
        */
        void push_front(DataCacheEntry* entry);

        /**
            Put the entry on the back of the list.  This is
            then the first block scheduled for reuse.
            This method is thread-safe.
            @param entry
        */
        void push_back(DataCacheEntry* entry);

        /**
            Move this block to the front of the queue
            @param entry
        */
        void reorder(DataCacheEntry* entry);

        /**
            Pull entry from the cache. The cache block cannot be
            cleared while it is pulled out
            @param entry
        */
        void pull(DataCacheEntry* entry);

        /**
            Pull a new data cache entry and assign the given data block
            as an owner. The cache block cannot be cleared while it is
            pulled
            @param block The block to assign as owner to the new cache entry
            @return A new data cache entry
        */
        DataCacheEntry* pull(CachedDataBlock* block);

        /**
            @return The first cache entry in the queue.  This would be the last
                    entry to be allocated upon a pull.
        */
        DataCacheEntry* first() const;

        /**
            @return The first cache entry in the queue.  This would be the first
                    entry to be allocated upon a pull.
        */
        DataCacheEntry* last() const;

        /**
            @return The number of available blocks in cache
        */
        uli ncache() const;

        /**
            @return The number of blocks with no linked data
        */
        uli nfree() const;

        /**
            @return The total number of blocks in the cache
        */
        uli ntotal() const;

        /**
            All cache blocks are of uniform size.
            @return The size of each block in the cache
        */
        uli blocksize() const;

        void print(std::ostream& os = std::cout) const;

        /**
            Deletes all entries in the cache.  DataCache does
            not own the actual memory allocation.  This only deletes
            the DataCacheEntry objects.  LayeredDataCache owns the actual
            memory allocation.
        */
        void free();

        /**
            Clear all cache entries which have not been pulled.
            This releases any ownership of data blocks
        */
        void clear();

        /**
            Initialize the data cache with a single entry
            @param ptr The data pointer beginning the cache entry
        */
        void init(char* ptr);

        /**
            Add a new block to the cache beginning
            @param data The data pointer beginning the new cache entry
        */
        void allocate_block(char* data);

};

/**
    @class LayeredDataCache
    Top-level class for managing cached memory blocks.  A tensor will have a non-uniform block
    size, but each data cache - for simplicity - is designed to hold only a single block size
    to avoid memory fragmentation.  A layered data cache controls a common memory pool and
    distributes memory from that pool to various individual data caches of a given block size.
    Cache blocks are meant to be allocated "dynamically" as a contraction proceeds so that
    the contraction optimizes the distribution of cache block sizes itself.
*/
class LayeredDataCache :
    public YetiRuntimeCountable
{

    private:
        /**
            The total amount of memory remaining to be allocated to data blocks
        */
        uli totalsize_;

        /**
            The total amount of memory remaining to be allocated to data cache entries
        */
        uli remaining_;

        /**
            The list of data caches
        */
        std::list<DataCache*> caches_;

        /**
            The list of block sizes for the cache list
        */
        std::list<size_t> sizes_;

        /**
            Pointer to the block of data this actually malloc'd
        */
        char* data_;

        /**
            Pointer that gets incremented to prepare for the next
            call
        */
        char* next_;

        /**
            Given a list #sizes_, create data cache objects for each
            of the possible sizes
        */
        void init();

        template <class T>
        void
        cache_action();

    public:
        /**
            @param Total size in bytes of available cache storage
            @param size_counts Key data size, value number of tiles with that size.
                               The data sizes will be used to determine which data
                               cache sizes are needed
        */
        LayeredDataCache(
            size_t storage,
            const std::map<size_t, size_t>& size_counts
        );

        /**
            Empty constructor for debugging purposes
        */
        LayeredDataCache();

        ~LayeredDataCache();

        /**
            @param The desired cache block size
            @return The cache with the smallest allowable
                    cache size greather than n.  Returns null
                    if blocksize is too large for cache.
        */
        DataCache* get_cache(size_t n) const;

        /**
            Find a suitable data cache for the requested block size
            and attempt to allocate memory from the common pool to
            to this data cache.  If not enough memory left in the common
            pool, the appropriate cache is just returned without allocating
            any new slots.  Every data cache has at least one spot.
            @param The desired cache block size
        */
        DataCache* allocate_cache(size_t n);

        /**
            @return The total size of the memory pool available to the data cache
        */
        uli totalsize() const;

        void print(std::ostream& os = std::cout) const;

        /**
            Allocate a set of cache blocks. This is mostly used for debugging purposes
            @param n
            @param size
        */
        void allocate_blocks(uli n, uli size);

        /**
            Clear all data cache blocks of ownership.  This leaves all allocation from the
            memory pool in place, but forces all links to owner data blocks to be cleared
            leaving all data cache entries free for use.
        */
        void clear();

        /**
            Performs all the operations of #clear, but additionally deallocates all the
            memory given out from the memory pool.  All data cache objects
            will no longer have any available memory except an initial first entry.
        */
        void free();


        void finalize_accumulate();

        void finalize_write();

};

}

#endif
