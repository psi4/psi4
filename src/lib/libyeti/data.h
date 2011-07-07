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
        typedef enum {
            in_core,
            in_cache,
            on_disk
        } storage_type_t;

        StorageBlock(size_t size);

        virtual ~StorageBlock();

        char* data() const;

        size_t size() const;

        virtual void clear(Cachable* owner) = 0;

        virtual bool retrieve() = 0;

        virtual void release() = 0;

        virtual void commit() = 0;

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

        bool retrieve();

        void release();

        void commit();

        void clear(Cachable* owner);

        void obsolete(Cachable* item);

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

        virtual void fetch_data();

    public:
        CachedStorageBlock(
            Cachable* parent,
            DataCache* cache
        );

        ~CachedStorageBlock();

        void assign_location(char* data);

        Cachable* get_cache_item() const;

        bool retrieve();

        void release();

        void clear(Cachable* owner);

        void free_cache_entry();

        void obsolete(Cachable* item);

        DataCacheEntry* get_entry() const;

        virtual void commit();

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

    protected:
        void fetch_data();

    public:
        LocalDiskBlock(
            Cachable* parent,
            DataCache* cache,
            const DiskBufferPtr& buffer
        );

        void commit();

};


}

#ifdef redefine_size_t
#undef size_t
#endif

#endif
