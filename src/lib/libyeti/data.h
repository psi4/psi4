#ifndef yeti_dataabstract_h
#define yeti_dataabstract_h

#include "class.h"

#include "data.hpp"
#include "cache.hpp"
#include "sort.hpp"
#include "index.hpp"
#include "permutation.hpp"
#include "tensor.hpp"
#include "thread.hpp"
#include "filler.hpp"

#include "mallocimpl.h"
#include "yetiobject.h"
#include "cache.h"

#define MALLOC_CACHE_BLOCKS 0

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {


class StorageBlock :
    public smartptr::Countable,
    public Malloc<StorageBlock>
{

    protected:
        char* data_;

        size_t size_;


    public:

        StorageBlock(size_t size);

        virtual ~StorageBlock();

        char* data() const;

        size_t size() const;

        virtual void clear() = 0;

        virtual bool retrieve() = 0;

        virtual void release() = 0;

        virtual void commit() = 0;

        virtual bool is_cached() const = 0;

        virtual bool is_retrieved() const = 0;

        virtual void obsolete(Cachable* item) = 0;

        void memset();

};

/**
    @param MemoryBlock
    Block of data that is persistent in memory
*/
class InCoreBlock :
    public StorageBlock
{

    public:
        InCoreBlock(char* data, size_t size);

        void commit();

        void clear();

        void obsolete(Cachable* item);

        bool retrieve();

        bool is_cached() const;

        bool is_retrieved() const;

        void release();


};

/**
    @class CachedDataBlock
    Block of data that is allocated from a central cache system.
    Caching is used to minimize fetches from disk/remote locations.
*/
class CachedStorageBlock :
    public StorageBlock
{

    protected:
        friend class DataCacheEntry;

        /**
            A data cache holding blocks of the appropriate size
            for this data block
        */
        DataCache* cache_;

        /**
            A cache entry containing a block of data that
            this cache entry is linked to
        */
        DataCacheEntry* cache_entry_;

        Cachable* cache_item_;

        bool retrieved_;


    public:
        CachedStorageBlock(
            Cachable* parent,
            DataCache* cache
        );

        ~CachedStorageBlock();

        void assign_location(char* data);

        Cachable* get_cache_item() const;

        bool is_cached() const;

        bool retrieve();

        void release();

        void clear();

        void obsolete(Cachable* item);

        DataCacheEntry* get_entry() const;

        virtual void commit();

        bool is_retrieved() const;

};


class DiskBuffer :
    public YetiRuntimeCountable
{

    private:
        std::string filename_;

        int fileno_;

        size_t size_;

    public:
        DiskBuffer(const std::string& filename);

        ~DiskBuffer();

        size_t size() const;

        /**
            @return The buffer offset the region begins at
        */
        size_t allocate_region(size_t size);

        void read(size_t offset, size_t size, char* buffer);

        void write(size_t offset, size_t size, const char* buffer);

        const std::string& filename() const;

};

class LocalDiskBlock :
    public CachedStorageBlock
{
    private:
        DiskBufferPtr buffer_;

        size_t offset_;

    public:
        LocalDiskBlock(
            Cachable* parent,
            DataCache* cache,
            const DiskBufferPtr& buffer
        );

        void commit();

        void fetch_data();

};


}

#ifdef redefine_size_t
#undef size_t
#endif

#endif
