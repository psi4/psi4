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

DataCacheEntry::DataCacheEntry(char* d)
    :
    next(0),
    prev(0),
    data(d),
    owner(0)
{
    yeti_register_new(lock_.get());
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
    if (owner) //let the owner know that data is gone
    {
        owner->clear();
    }
    block->data()->reference(data);
    owner = block;
}

void
DataCacheEntry::print(std::ostream& os) const
{
    os << stream_printf("Prev:%12p This:%12p Next:%12p Owner:%12p Data:%12p", prev, this, next, owner, data);
}

DataCache::DataCache(
    size_t blocksize,
    char *datastart
) :
    first_(0),
    last_(0),
    n_(0),
    blocksize_(blocksize)
{
    yeti_register_new(lock_.get());
    init(datastart);
}

DataCache::~DataCache()
{
    free();
}

void
DataCache::free()
{
    DataCacheEntry* next = first_;
    while (next) //loop until hitting zero
    {
        DataCacheEntry* entry = next;
        next = next->next; //grab the next element
        delete entry;
    }
}

void
DataCache::clear()
{
    DataCacheEntry* next = first_;
    while (next) //loop until hitting zero
    {
        DataCacheEntry* entry = next;
        next = next->next; //grab the next element
        entry->clear();
    }
}

void
DataCache::init(char* ptr)
{
    first_ = new DataCacheEntry(ptr);
    last_ = first_;
    first_->prev = 0; //points nowhere
    first_->next = 0;
    n_ = 1;
}

void
DataCache::allocate_block(char *ptr)
{
    DataCacheEntry* entry = new DataCacheEntry(ptr);
    entry->retrieve(); //"retrieve" the newly created entry for pushing back
    //this will lock on push back
    push_back(entry);
    ++n_;
}

DataCacheEntry*
DataCache::pull(CachedDataBlock* block)
{
    lock();

    if (last_ == 0)
    {
        raise(SanityCheckError, "too many active cache blocks.  no more free blocks");
    }

    //grab temporary link to last value
    DataCacheEntry* tmp = last_;

    //move the last pointer back one and point it to the end
    last_ = last_->prev;

    if (last_) //only set to 0 if it still exists
        last_->next = 0;
    else
        first_ = 0; //nothing left

    tmp->reassign(block);
    tmp->retrieve();

    unlock();

    return tmp;
}

size_t
DataCache::ncache() const
{
    DataCacheEntry* next = first_;
    size_t count = 0;
    while (next) //loop until hitting zero
    {
        next = next->next; //grab the next element
        ++count;
    }
    return count;
}

size_t
DataCache::blocksize() const
{
    return blocksize_;
}

size_t
DataCache::ntotal() const
{
    return n_;
}

size_t
DataCache::nfree() const
{
    DataCacheEntry* next = first_;
    size_t count = 0;
    while (next) //loop until hitting zero
    {
        if (next->owner == 0)
            ++count;
        next = next->next; //grab the next element
    }
    return count;
}

DataCacheEntry*
DataCache::first() const
{
    return first_;
}

DataCacheEntry*
DataCache::last() const
{
    return last_;
}

void
DataCache::push_front(DataCacheEntry* entry)
{
    lock();
    if (!entry->is_retrieved())
    //another thread came in and did this while we were waiting
    {
        unlock();
        return;
    }

    if (first_) //set the pointer
    {
        first_->prev = entry;
        entry->next = first_;
    }
    else //this is now the only entry, first and last
    {
        last_ = entry;
        entry->next = 0;
    }

    first_ = entry;
    entry->prev = 0; //nothing ahead

    entry->release();

    unlock();
}

void
DataCache::push_back(DataCacheEntry* entry)
{
    lock();
    if (!entry->is_retrieved())
    //another thread came in and did this while we were waiting
    {
        unlock();
        return;
    }

    if (last_) //set the pointer
    {
        last_->next = entry;
        entry->prev = last_;
    }
    else //this is now the only entry, first and last
    {
        first_ = entry;
        entry->prev = 0;
    }

    last_ = entry;
    entry->next = 0; //nothing ahead
    entry->release();

    unlock();
}


void
DataCache::pull(DataCacheEntry* entry)
{
    lock();
    if (entry->is_retrieved())
    {
        unlock();
        return;
    }

    //swap entries between around this entry so as to pull it from the list
    if (entry == first_)
    {
        first_ = entry->next;
    }
    else
    {
        entry->prev->next = entry->next;
    }

    if (entry == last_)
    {
        last_ = entry->prev;
    }
    else
    {
        entry->next->prev = entry->prev;
    }

    entry->retrieve();

    unlock();
}

void
DataCache::print(std::ostream &os) const
{
    os << Env::indent << "Data Cache " << n_ << " blocks of size " << blocksize_;
    DataCacheEntry* next = first_;
    while (next) //loop until hitting zero
    {
        DataCacheEntry* entry = next;
        next = next->next; //grab the next element
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
            raise(SanityCheckError, "insufficient memory to allocate cache blocks");
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
    yeti_register_new(lock_.get());
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
LayeredDataCache::allocate_cache(size_t size)
{
    lock();
    DataCache* cache = get_cache(size);
    //this is the correct cache size


    if (cache && remaining_ > size) //we have enough space to allocate a new block
    {
        cache->allocate_block(next_);
        next_ += size;
        remaining_ -= size;
    }

    unlock();
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
    DataCacheEntry* next = first_;
    while (next) //loop until hitting zero
    {
        DataCacheEntry* entry = next;
        next = next->next; //grab the next element
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
