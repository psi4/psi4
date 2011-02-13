#include "cache.h"
#include "exception.h"
#include "env.h"
#include "data.h"
#include "tensor.h"
#include "tile.h"
#include "tuple.h"
#include "index.h"
#include "malloc.h"
#include "runtime.h"
#include "thread.h"

#include <algorithm>

using namespace yeti;
using namespace std;

#define MALLOC_CACHE_ENTRIES 0

DECLARE_MALLOC(DataCacheEntry);

DataCacheEntry::DataCacheEntry(
    char* d,
    uli posoffset
)
    :
    data(d),
    owner(0),
    offset(posoffset)
{
}

DataCacheEntry::~DataCacheEntry()
{
    clear();
}

void
DataCacheEntry::clear()
{
    if (owner)
    {
        if (owner->is_retrieved())
        {
            cerr << "cannot clear retrieved owner!" << endl;
            abort();
        }

        owner->clear();
        owner = 0;
    }
}

void
DataCacheEntry::finalize_accumulate()
{
    if (owner)
        owner->finalize_accumulate(NOT_THREADED);
}

void
DataCacheEntry::finalize_write()
{
    if (owner)
        owner->finalize_write(NOT_THREADED);
}

void
DataCacheEntry::reassign(CachedDataBlock* block)
{
    /** JJW 11/06/10
        On the off chance owner = d, owner must be zeroed first!
        Otherwise you assign the new pointer and then immediately zero it!
        ORDER MATTERS!
    */
    CachedDataBlock* old_owner = owner;
    owner = block;
    if (old_owner) //let the owner know that data is gone
    {

        old_owner->lock();
        if (old_owner->cache_entry_ == this) //I am responsible for clearing
        {
            #if YETI_SANITY_CHECK
            if (owner->is_retrieved())
            {
                void* entry = owner->cache_entry_;
                cerr << "unlocked? " << trylock() << endl;
                cerr << "old owner: " << (void*) owner << endl;
                cerr << "new owner: " << (void*) block << endl;
                cerr << "me: " << (void*) this << " owner:" << entry << endl;
                cerr << "cannot clear retrieved owner!" << endl;
                abort();
            }
            #endif
            old_owner->clear();
        }
        old_owner->unlock();
    }
    block->data()->reference(data);
}

void
DataCacheEntry::print(std::ostream& os) const
{
    os << stream_printf("This:%12p Offset:%ld Owner:%12p Data:%12p",
    this, offset, owner, data);
}

#define NCACHE_ENTRIES_MAX 100000

DataCache::DataCache(
    uli blocksize,
    char *datastart
) :
    nentries_(0),
    blocksize_(blocksize),
    entries_(new DataCacheEntry*[NCACHE_ENTRIES_MAX]),
    next_offset_(0)
{
    init(datastart);
}

DataCache::~DataCache()
{
    free();
    delete[] entries_;
}

void
DataCache::free()
{
    for (uli i=0; i < nentries_; ++i)
        delete entries_[i];
    nentries_ = 0;
}

void
DataCache::clear()
{
    for (uli i=0; i < nentries_; ++i)
    {
        DataCacheEntry* entry = entries_[i];
        if (!entry)
        {
            cerr << "cache should not be cleared until all entries have been returned" << endl;
            abort();
        }
        if (YetiRuntime::is_threaded_runtime() && entry->trylock())
        {
            cerr << "entry should not be unlocked!" << endl;
            abort();
        }
        entry->clear();
    }
}

void
DataCache::init(char* ptr)
{
    for (uli i=0; i < NCACHE_ENTRIES_MAX; ++i)
        entries_[i] = 0;

    allocate_block(ptr);
}

void
DataCache::allocate_block(char *ptr)
{
    if (nentries_ > NCACHE_ENTRIES_MAX)
    {
        cerr << "Cannot allocate any further cache entries" << endl;
        abort();
    }
    uli offset = nentries_;
    DataCacheEntry* entry = new DataCacheEntry(ptr, offset);
    entries_[nentries_] = entry;
    ++nentries_;
}

DataCacheEntry**
DataCache::entries() const
{
    return entries_;
}

DataCacheEntry*
DataCache::find_entry(uli offset)
{
    for (uli i=offset; i < nentries_; ++i)
    {
        DataCacheEntry* entry = entries_[i];
        if (entry && entry->trylock())
        {
            return entry;
        }
    }
    return 0;
}

DataCacheEntry*
DataCache::pull(CachedDataBlock* block, uli offset)
{
    DataCacheEntry* entry = find_entry(offset);
    if (!entry)
        entry = find_entry(0); //look again starting from zero

    if (!entry)
    {
        cerr << "too many active cache blocks.  no more free blocks" << endl;
        abort();
    }

    if (YetiRuntime::is_threaded_runtime() && entry->trylock())
    {
        cerr << "The cache entry should have been locked!" << endl;
        abort();
    }

    entry->reassign(block);
    entries_[entry->offset] = 0;
    return entry;
}

uli
DataCache::ncache() const
{
    uli ntot = 0;
    for (uli i=0; i < nentries_; ++i)
    {
        DataCacheEntry* entry = entries_[i];
        if (entry)
            ++ntot;

    }
    return ntot;
}

uli
DataCache::blocksize() const
{
    return blocksize_;
}

uli
DataCache::ntotal() const
{
    return nentries_;
}

uli
DataCache::nfree() const
{
    uli ntot = 0;
    for (uli i=0; i < nentries_; ++i)
    {
        DataCacheEntry* entry = entries_[i];
        if (entry && entry->owner == 0)
            ++ntot;
    }
    return ntot;
}

void
DataCache::insert(DataCacheEntry* entry)
{
    entries_[entry->offset] = entry;
}


void
DataCache::pull(DataCacheEntry* entry)
{
    entries_[entry->offset] = 0;
}

void
DataCache::print(std::ostream &os) const
{
    os << Env::indent << "Data Cache " << nentries_ << " blocks of size " << blocksize_;
    for (uli i=0; i < nentries_; ++i)
    {
        DataCacheEntry* entry = entries_[i];
        os << endl << Env::indent;
        entry->print(os);
    }
}



LayeredDataCache::LayeredDataCache(
    size_t storage,
    const std::map<size_t, size_t>& size_counts
) :
    data_(new char[storage]),
    totalsize_(storage)
{
    std::map<size_t, size_t>::const_iterator it = size_counts.begin();
    std::map<size_t, size_t>::const_iterator stop = size_counts.end();
    for ( ; it != stop; ++it)
    {
        size_t size = it->first;
        sizes_.push_back(size);
    }
    init();
}

#define DEBUG_CACHE_BLOCK_SIZE 500000
#define DEBUG_CACHE_TOTAL_MEM 10000000 //10 mb
LayeredDataCache::LayeredDataCache()
  :
    data_(new char[DEBUG_CACHE_TOTAL_MEM]),
    totalsize_(DEBUG_CACHE_TOTAL_MEM)
{
    sizes_.push_back(DEBUG_CACHE_BLOCK_SIZE);
    init();
}

LayeredDataCache::~LayeredDataCache()
{
    //delete all the vector entries
    list<DataCache*>::const_iterator it(caches_.begin());
    list<DataCache*>::const_iterator stop(caches_.end());
    for ( ; it != stop; ++it)
    {
        DataCache* cache(*it);
        delete cache;
    }
    delete[] data_;
}


void
LayeredDataCache::init()
{
    sizes_.sort();

    std::list<size_t>::const_iterator it(sizes_.begin());
    std::list<size_t>::const_iterator stop(sizes_.end());
    next_ = data_;
    remaining_ = totalsize_;
    for ( ; it != stop; ++it)
    {
        size_t size(*it);
        if (remaining_ < size)
        {
            raise(SanityCheckError, "insufficient memory to allocate cache bs");
        }
        DataCache* cache(new DataCache(size, next_));
        next_ += size;
        remaining_ -= size;
        caches_.push_back(cache);
    }
}

void
LayeredDataCache::allocate_blocks(size_t n, size_t size)
{
    DataCache* cache = get_cache(size);

    if (!cache)
    {
        raise(SanityCheckError, "no data cache available for requested size");
    }

    for (size_t i=0; i < n; ++i)
    {
        if (size > remaining_)
            raise(SanityCheckError, "insufficient size remaining for cache allocation");

        cache->allocate_block(next_);
        next_ += size;
        remaining_ -= size;
    }
}

size_t
LayeredDataCache::totalsize() const
{
    return totalsize_;
}

void
LayeredDataCache::free()
{
    list<DataCache*>::iterator it(caches_.begin());
    list<DataCache*>::iterator stop(caches_.end());
    list<size_t>::iterator it_size(sizes_.begin());
    next_ = data_;
    remaining_ = totalsize_;
    for ( ; it != stop; ++it, ++it_size)
    {
        size_t size(*it_size);
        DataCache* cache(*it);
        cache->free();
        cache->init(next_);
        next_ += size;
        remaining_ -= size;
    }
}

void
LayeredDataCache::clear()
{
    list<DataCache*>::iterator it(caches_.begin());
    list<DataCache*>::iterator stop(caches_.end());
    list<size_t>::iterator it_size(sizes_.begin());
    for ( ; it != stop; ++it, ++it_size)
    {
        size_t size(*it_size);
        DataCache* cache(*it);
        cache->clear();
    }
}

DataCache*
LayeredDataCache::allocate_cache(size_t size, uli& offset)
{
    DataCache* cache = get_cache(size);
    //this is the correct cache size

    uli blocksize = cache->blocksize();
    if (cache && remaining_ > blocksize) //we have enough space to allocate a new block
    {
        cache->allocate_block(next_);
        next_ += blocksize;
        remaining_ -= blocksize;
    }



    if (cache)
    {
        offset = cache->next_offset_;
        cache->next_offset_ = (offset + 1) % cache->nentries_;
    }

    return cache;
}

DataCache*
LayeredDataCache::get_cache(
    size_t n
) const
{
    list<DataCache*>::const_iterator it(caches_.begin());
    list<DataCache*>::const_iterator stop(caches_.end());

    list<size_t>::const_iterator it_size(sizes_.begin());

    for ( ; it != stop; ++it, ++it_size)
    {
        size_t size(*it_size);
        if (size >= n) //finally big enough
            return *it;
    }
    return 0;
}

void
LayeredDataCache::print(std::ostream &os) const
{
    os << Env::indent << "Layered Data Cache";
    //delete all the vector entries
    list<DataCache*>::const_iterator it(caches_.begin());
    list<DataCache*>::const_iterator stop(caches_.end());
    ++Env::indent;
    for ( ; it != stop; ++it)
    {
        DataCache* cache(*it);
        os << endl;
        cache->print(os);
    }
    --Env::indent;
}


template <class T>
void
LayeredDataCache::cache_action()
{
    list<DataCache*>::const_iterator it(caches_.begin());
    list<DataCache*>::const_iterator stop(caches_.end());
    for ( ; it != stop; ++it)
    {
        DataCache* cache(*it);
        cache->cache_action<T>();
    }
}

class FinalizeAccumulate {
    public:
        static void action(DataCacheEntry* entry)
        {
            entry->finalize_accumulate();
        }

};

class FinalizeWrite {
    public:
        static void action(DataCacheEntry* entry)
        {
            entry->finalize_write();
        }
};


template <class T>
void
DataCache::cache_action()
{
    for (uli i=0; i < nentries_; ++i)
    {
        DataCacheEntry* entry = entries_[i];
        if (!entry)
        {
            cerr << "Cannot perform action on cache in which all blocks have not been returned" << endl;
            abort();
        }
        T::action(entry);
    }
}


void
LayeredDataCache::finalize_accumulate()
{
    cache_action<FinalizeAccumulate>();
}

void
LayeredDataCache::finalize_write()
{
    cache_action<FinalizeWrite>();
}
