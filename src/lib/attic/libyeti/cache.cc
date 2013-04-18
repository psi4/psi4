#include "cache.h"
#include "exception.h"
#include "env.h"
#include "data.h"
#include "tensor.h"
#include "tensorblock.h"
#include "tuple.h"
#include "index.h"
#include "malloc.h"
#include "runtime.h"
#include "thread.h"

#include <algorithm>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define MALLOC_CACHE_ENTRIES 0

#define DEBUG_CACHE_PRINT 0

//#define name(x) static_cast<TensorBlock*>(x)->get_block_name()

DECLARE_MALLOC(DataCacheEntry);
DECLARE_MALLOC(DataCache);

DataCacheEntry::DataCacheEntry(
    char* d
)
    :
    data(d),
    owner(0),
    pulled(false),
    next(0),
    prev(0),
#if DEBUG_CACHE_HISTORY
    nentries_history_(0),
#endif
    threadnum(0)
{
}

DataCacheEntry::~DataCacheEntry()
{
}

void
DataCacheEntry::print(std::ostream& os) const
{
    os << stream_printf("This:%12p Owner:%12p Data:%12p",
                            this, owner, data);
}

void*
DataCacheEntry::operator new(size_t size)
{
    return Malloc<DataCacheEntry>::operator new(size);
}

void
DataCacheEntry::operator delete(void* ptr)
{
    Malloc<DataCacheEntry>::operator delete(ptr);
}

void*
DataCacheEntry::operator new(size_t size, MemoryPool* mempool)
{
    return mempool->get(size);
}

void
DataCacheEntry::operator delete(void* ptr, MemoryPool* mempool)
{
    //do nothing
}

// For making "extensions" to caches temporarily.  Currently, not used (though it has been debugged).
DataCache::DataCache(
    DataCache* parent,
    MemoryPool* mempool,
    uli nelements,
    size_t blocksize
) : 
    nentries_(nelements),
    blocksize_(blocksize),
#if DEBUG_CACHE_HISTORY
    nhistories_(0),
#endif
    free_entries_start_(0),
    used_entries_start_(0),
    used_entries_end_(0),
    cache_extension_(0),
    extension_block_(0),
    parent_(parent),
    lock_(parent->lock_)
{
    if (nelements == 0)
        yeti_throw(SanityCheckError, "cannot build a cache extension with zero entries");

    char* data = mempool->get(blocksize);
    DataCacheEntry* entry  = new (mempool) DataCacheEntry(data);
    free_entries_start_ = entry;
    for (uli i=1; i < nelements; ++i)
    {
        data = mempool->get(blocksize);
        DataCacheEntry* new_entry  = new (mempool) DataCacheEntry(data);
        entry->next = new_entry;
        new_entry->prev = entry;
        entry = new_entry;
    }
}

DataCache::DataCache(
    MultiThreadLock* lock,
    size_t blocksize,
    const std::string& name
) :
    nentries_(0),
    blocksize_(blocksize),
    free_entries_start_(0),
    used_entries_start_(0),
    used_entries_end_(0),
    cache_extension_(0),
    extension_block_(0),
    parent_(0),
    name_(name),
#if DEBUG_CACHE_HISTORY
    nhistories_(0),
#endif
    lock_(lock)
{
    YetiRuntime::register_cache(this);
}

const std::string&
DataCache::name() const
{
    return name_;
}

DataCacheEntry*
DataCache::allocate_new_entry()
{
    DataCacheEntry* entry = new DataCacheEntry(YetiRuntime::malloc(blocksize_));
    entry->next = 0;
    entry->prev = 0;
    all_entries_[nentries_] = entry;
    ++nentries_;
    YetiRuntime::register_allocation(this, blocksize_);
    return entry;
}

DataCache::~DataCache()
{
    YetiRuntime::unregister_cache(this);
#if DEBUG_CACHE_PRINT
    cout << stream_printf("Deleting data cache %p\n", this);
    cout.flush();
#endif
    if (parent_)
        yeti_throw(SanityCheckError, "Cannot delete cache extension");
    
#if YETI_SANITY_CHECK
    for (uli i=0; i < nentries_; ++i)
    {
        if (all_entries_[i]->pulled)
        {
            TensorBlock* block = static_cast<TensorBlock*>(all_entries_[i]->owner);        
            if (block)
                cerr << block->get_block_name().c_str() << endl;
            yeti_throw(SanityCheckError, "Not all cache entries returned to cache!");
        }
        //find the entry amongst everything
        DataCacheEntry* entry = free_entries_start_;
        while (entry && entry != all_entries_[i])
            entry = entry->next;

        if (entry != all_entries_[i])
        {
            entry = used_entries_start_;
            while (entry && entry != all_entries_[i])
                entry = entry->next;
        }

        if (entry != all_entries_[i])
        {
#if DEBUG_CACHE_HISTORY
            print_history();
#endif
            TensorBlock* block = static_cast<TensorBlock*>(all_entries_[i]->owner);        
            if (block)
                cerr << block->get_block_name().c_str() << endl;
            yeti_throw(SanityCheckError, "Not all cache entries returned to cache!");
        }
    }
#endif

    DataCacheEntry* entry = free_entries_start_;
    while(entry)
    {
        DataCacheEntry* old = entry;
        entry = entry->next;
        if (entry == free_entries_start_)
            yeti_throw(SanityCheckError, "deleting start again!");
        YetiRuntime::free(old->data, blocksize_);
        delete old;
    }

    entry = used_entries_start_;
    while (entry)
    {
        if (entry->owner)
        {
            entry->owner->lock();
            entry->owner->flush_from_cache();
            entry->owner->lock();
        }
        DataCacheEntry* old = entry;
        entry = entry->next;
        YetiRuntime::free(old->data, blocksize_);
        delete old;
    }
}

DataCacheEntry*
DataCache::find_free_entry()
{
    DataCacheEntry* entry = free_entries_start_;
    if (!entry) //no free entries
    {
        return 0;
    }

    free_entries_start_ = entry->next;
    if (free_entries_start_)
        free_entries_start_->prev = 0;

    return entry;
}

void
DataCache::flush()
{
#if DEBUG_CACHE_PRINT
    cout << stream_printf("Flushing cache %p\n", this);
    cout.flush();
#endif
    if (used_entries_start_ == 0) //nothing to do
        return;

    while (used_entries_start_) /** this is constantly getting added to free entries */
    {
        DataCacheEntry* entry = used_entries_start_;
        Cachable* owner = entry->owner;
        owner->lock();
        if (owner->is_retrieved())
            yeti_throw(SanityCheckError, "cache owner is being flush in part of code with retrieved blocks");
        owner->flush_from_cache();  //I think this is setting the cache entry to zero
        owner->unlock();

        entry->owner = 0;
    }

    //all used entries have become free entries in the owner->flush_from_cache calls
}

DataCacheEntry*
DataCache::find_used_entry()
{
    DataCacheEntry* entry = used_entries_start_;
    if (!entry)
        return 0;

    while (entry)
    {
        if (entry->owner->trylock())
        {
            Cachable* owner = entry->owner;
#if YETI_SANITY_CHECK
            if (owner->is_retrieved())
            {
                cerr << stream_printf("Failure in cache for block %s at %p\n", 
                    static_cast<TensorBlock*>(owner)->get_block_name().c_str(),
                    owner);
                cerr.flush();
                yeti_throw(SanityCheckError, "Owner is retrieved but still resides in cache");
            }
#endif
            try{
                owner->flush_from_cache(); //entry->owner = 0 after this call
                owner->unlock();
            } catch (int e)
            {
                abort();
            }

            //there is now a free entry
            return find_free_entry();
        }
        entry = entry->next;
    }
    return 0;
}


DataCacheEntry*
DataCache::pull(DataCacheEntry* entry, Cachable* owner)
{
    if (cache_extension_)
        return cache_extension_->pull(entry, owner);

    lock();
    if (entry->owner != owner) //the owner is not still valid
    {
        cout << stream_printf("****** Pulled new entry for %s *******\n",
            static_cast<TensorBlock*>(owner)->get_block_name().c_str());
        DataCacheEntry* new_entry = pull_no_lock(owner);
        unlock();
        return new_entry;
    }

    if (entry->prev)
    {
        entry->prev->next = entry->next;
    }
    else //beginning
    {
        used_entries_start_ = entry->next;
    }

    if (entry->next)
    {
        entry->next->prev = entry->prev;
    }
    else //end
    {
        used_entries_end_ = entry->prev;
    }

#if YETI_SANITY_CHECK
    if (used_entries_end_ && used_entries_end_->next)
    {
        yeti_throw(SanityCheckError, "cache alignment fail");
    }
    if (used_entries_start_ && used_entries_start_->prev)
    {
        yeti_throw(SanityCheckError, "cache alignment fail");
    }
#endif

    entry->pulled = true;

#if DEBUG_CACHE_PRINT
    cout << stream_printf("Repulled entry %p for owner %p on cache %p\n",
                entry, owner, this);
    cout.flush();
#endif

#if DEBUG_CACHE_HISTORY
    CacheHistory& his = entry->histories_[entry->nentries_history_];
    his.thread_number = YetiRuntime::get_thread_number();
    his.owner = owner;
    ++entry->nentries_history_;
    if (entry->nentries_history_ == NHISTORIES_CACHE)
        entry->nentries_history_ = 0;

    CacheHistory& main = all_histories_[nhistories_];
    main.thread_number = his.thread_number;
    main.owner = owner;
    ++nhistories_;
#endif

    unlock();


    return entry;
}

DataCacheEntry*
DataCache::pull_no_lock(Cachable* owner)
{
    DataCacheEntry* entry = find_free_entry();
    if (!entry) 
    {
       try {
        entry = find_used_entry();
       } catch (int e) {
          if (e == TENSOR_BLOCK_POLICY_EXCEPTION)
          {
            TensorBlock* block = static_cast<TensorBlock*>(owner);
            cerr << "Flush error in cache pull for " << block->get_block_name() << endl;
            throw e;
          }
        }
#if DEBUG_CACHE_PRINT
        cout << stream_printf("Pulled used entry %p for owner %p on cache %p\n",
                    entry, owner, this);
        cout.flush();
#endif
    }
    else
    {
#if DEBUG_CACHE_PRINT
        cout << stream_printf("Pulled free entry %p for owner %p on cache %p\n",
                    entry, owner, this);
        cout.flush();
#endif
    }

    if (!entry)
    { 
        entry = allocate_new_entry();
    }


    entry->owner = owner;
    entry->pulled = true;

#if DEBUG_CACHE_HISTORY
    CacheHistory& his = entry->histories_[entry->nentries_history_];
    his.thread_number = YetiRuntime::get_thread_number();
    his.owner = owner;
    ++entry->nentries_history_;
    if (entry->nentries_history_ == NHISTORIES_CACHE)
        entry->nentries_history_ = 0;

    CacheHistory& main = all_histories_[nhistories_];
    main.thread_number = his.thread_number;
    main.owner = owner;
    ++nhistories_;
#endif
        
    return entry;
}

DataCacheEntry*
DataCache::pull(Cachable* owner)
{
    if (cache_extension_)
        return cache_extension_->pull(owner);

    lock();
    DataCacheEntry* entry = pull_no_lock(owner);
    unlock();

    return entry;
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

void
DataCache::build_extension(uli nelements)
{

    if (parent_)
        yeti_throw(SanityCheckError, "cannot build a cache extension for non-parent cache");

    if (cache_extension_)
        return;

    if (nentries_ >= nelements)
        return; //nothing to do

#if DEBUG_CACHE_PRINT
    cout << stream_printf("Building extension on %p\n", this);
    cout.flush();
#endif

    uli nentries_new = nelements - nentries_;

    size_t extension_size =     
        nentries_new * blocksize_ + sizeof(DataCache) + nentries_new * sizeof(DataCacheEntry);
    extension_block_ = YetiRuntime::malloc(extension_size);
    MemoryPool mempool(extension_size, extension_block_);
    cache_extension_ = new (&mempool) DataCache(this, &mempool, nentries_new, blocksize_);

    if (used_entries_start_)
        cache_extension_->append_used_entries(used_entries_start_);
    if (free_entries_start_)
        cache_extension_->append_free_entries(free_entries_start_);
    cache_extension_->nentries_ += nentries_;

}

void
DataCache::append_used_entries(DataCacheEntry* head_entry)
{
    if (used_entries_end_)
    {
        used_entries_end_->next = head_entry;
        head_entry->prev = used_entries_end_;
    }
    else
    {
        //no entries
        used_entries_start_ = head_entry;
    }

    DataCacheEntry* entry = head_entry;
    while(entry->next)
    {
        entry = entry->next;
    }
    used_entries_end_ = entry;
}

void
DataCache::append_free_entries(DataCacheEntry* head_entry)
{
    if (free_entries_start_)
    {
        DataCacheEntry* entry = free_entries_start_;
        while (entry->next)
        {
            entry = entry->next;
        }
        entry->next = head_entry;
        head_entry->prev = entry;
    }
    else
    {
        free_entries_start_ = head_entry;
    }
}


#if 0
// Since this is not currently used, we haven't updated the new DataCache object
// to use extensions
void
DataCache::delete_extension()
{
    abort();
    if (cache_extension_)
    {
#if DEBUG_CACHE_PRINT
    cout << stream_printf("Deleting extension on %p\n", this);
    cout.flush();
#endif

        cache_extension_->flush();

        DataCacheEntry* entry = cache_extension_->free_entries_start_;
        DataCacheEntry* prev = 0;
        char* data_start = data_;
        char* data_stop = data_ + nentries_ * blocksize_;

        while (entry)
        {
        }

        delete cache_extension_->lock_;
        cache_extension_ = 0;
        extension_block_ = 0;
        delete[] extension_block_;
    }
    
}
#endif

uli
DataCache::nfree() const
{
    cerr << "not implemented" << endl;
    abort();
}

void
DataCache::insert(DataCacheEntry* entry)
{
    if (cache_extension_)
    {
        cache_extension_->insert(entry);
        return;
    }

    lock();

    entry->pulled  = false;
    DataCacheEntry* last = used_entries_end_;
    if (last == 0)
    {
        used_entries_start_ = entry;
        used_entries_end_ = entry;
        entry->prev = 0;
        entry->next = 0;
    }
    else
    {
        last->next = entry;
        entry->prev = last;
        entry->next = 0;
        used_entries_end_ = entry;
    }

#if YETI_SANITY_CHECK
    if (used_entries_end_->next)
    {
        yeti_throw(SanityCheckError, "Cache alignment fail");
    }
    if (used_entries_start_->prev)
    {
        yeti_throw(SanityCheckError, "Cache alignment fail");
    }
#endif

#if DEBUG_CACHE_PRINT
    cout << stream_printf("Inserted entry %p for owner %p on cache %p\n",
                entry, entry->owner, this);
    cout.flush();
#endif

    unlock();

}

void
DataCache::free(DataCacheEntry* entry, Cachable* owner)
{
#if YETI_SANITY_CHECK
    if (!entry->owner->is_locked())
        yeti_throw(SanityCheckError, "cannot free cache entry of unlocked owner");
    if (entry->pulled)
        yeti_throw(SanityCheckError, "cannot free a pulled entry");
#endif

    if (cache_extension_)
    {
        cache_extension_->free(entry, owner);
        return;
    }

    lock();

    if (entry->owner != owner)
    {
        unlock();
        return; //nothing left to free, something else is using it
    }

    

    /** yank out of the used entries */
    if (entry->prev)
        entry->prev->next = entry->next;
    else //beginning
        used_entries_start_ = entry->next;

    if (entry->next)
        entry->next->prev = entry->prev;
    else //end
        used_entries_end_ = entry->prev;

#if YETI_SANITY_CHECK
    if (used_entries_end_ && used_entries_end_->next)
    {
        yeti_throw(SanityCheckError, "Cache alignment fail");
    }
    if (used_entries_start_ && used_entries_start_->prev)
    {
        yeti_throw(SanityCheckError, "Cache alignment fail");
    }
#endif

    /** put into the free entries */
    DataCacheEntry* second = free_entries_start_;
    if (second == 0)
    {
        entry->next = 0;
        entry->prev = 0;
    }
    else
    {
        entry->next = second;
        second->prev = entry;
        entry->prev = 0;
    }
    entry->pulled  = false;
    entry->owner = 0;
    free_entries_start_ = entry;

#if YETI_SANITY_CHECK
    if (used_entries_end_ && used_entries_end_->next)
    {
        yeti_throw(SanityCheckError, "cache alignment fail");
    }
    if (used_entries_start_ && used_entries_start_->prev)
    {
        yeti_throw(SanityCheckError, "cache alignment fail");
    }
    if (free_entries_start_->prev)
    {
        yeti_throw(SanityCheckError, "cache alignment fail");
    }
#endif

#if DEBUG_CACHE_PRINT
    cout << stream_printf("Freed entry %p for owner %p on cache %p\n",
                entry, owner, this);
    cout.flush();
#endif

    unlock();
}

void
DataCache::lock()
{
    lock_->lock();
}

void
DataCache::unlock()
{
    lock_->unlock();
}

void
DataCache::print(std::ostream &os) const
{
    os << Env::indent << "Data Cache " << endl;
}

#if DEBUG_CACHE_HISTORY
void
DataCache::print_history()
{
    cout << stream_printf("Cache history for cache with %ld blocks of size %ld\n", nentries_, blocksize_);
    for (uli i=0; i < nhistories_; ++i)
    {
        CacheHistory& h = all_histories_[i];
        TensorBlock* block = static_cast<TensorBlock*>(h.owner);
        if (block)
            cout << stream_printf("Entry pulled by thread %ld by block %s\n", h.thread_number, block->get_block_name().c_str());
    }
    cout << "***************************" << endl;
}
#endif

void*
DataCache::operator new(size_t size)
{
    return Malloc<DataCache>::operator new(size);
}

void
DataCache::operator delete(void* ptr)
{
    Malloc<DataCache>::operator delete(ptr);
}

void*
DataCache::operator new(size_t size, MemoryPool* mempool)
{
    return mempool->get(size);
}

void
DataCache::operator delete(void* ptr, MemoryPool* mempool)
{
    //do nothing
}

Cachable::Cachable()
    : in_destructor_(false)
{
}

Cachable::Cachable(YetiRuntimeObject::thread_safety_flag_t flag)
 : YetiRuntimeCountable(flag), in_destructor_(false)
{
}

bool
Cachable::in_destructor() const
{
    return in_destructor_;
}


