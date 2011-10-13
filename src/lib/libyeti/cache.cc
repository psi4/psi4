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
    threadnum(0)
{
    if (!lock_)
        abort();
}

DataCacheEntry::~DataCacheEntry()
{
}

void
DataCacheEntry::assign(Cachable* new_owner)
{
#if YETI_SANITY_CHECK
    if (!is_locked())
        raise(SanityCheckError, "cannot assign to an unlocked cache entry");
#endif
    owner = new_owner;
}

void
DataCacheEntry::flush()
{
    if (owner)
        owner->flush_from_cache();
    owner = 0;
}

Cachable*
DataCacheEntry::get_owner() const
{
    return owner;
}

void
DataCacheEntry::print(std::ostream& os) const
{
    os << stream_printf("This:%12p Owner:%12p Data:%12p",
                            this, owner, data);
}

bool
DataCacheEntry::is_pulled() const
{
    return pulled;
}

DataCache::DataCache(
    size_t storage,
    size_t blocksize
) :
    nentries_(0),
    blocksize_(blocksize),
    data_(0),
    storage_(storage),
    free_entries_start_(0),
    used_entries_start_(0),
    used_entries_end_(0)
{
    data_ = YetiRuntime::malloc(storage);

    char* ptr = data_;

    nentries_ = storage / blocksize_;

    free_entries_start_ = new DataCacheEntry(ptr);
    DataCacheEntry* entry = free_entries_start_;
    ptr += blocksize_;
    for (uli i=1; i < nentries_; ++i, ptr += blocksize_)
    {
        DataCacheEntry* new_entry = new DataCacheEntry(ptr);
        entry->next = new_entry;
        new_entry->prev = entry;
        entry = new_entry;
    }
    used_entries_start_ = 0;
    used_entries_end_ = 0;
}

DataCache::~DataCache()
{
    DataCacheEntry* entry = free_entries_start_;
    while(entry)
    {
        DataCacheEntry* old = entry;
        entry = entry->next;
        if (entry == free_entries_start_)
            raise(SanityCheckError, "deleting start again!");
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
        delete old;
    }

    YetiRuntime::free(data_, storage_);
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
    entry->lock();
    return entry;
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
#if YETI_SANITY_CHECK
            if (entry->owner->is_retrieved())
                raise(SanityCheckError, "owner is retrieved but still resides in cache");
#endif
            entry->lock();
            entry->owner->flush_from_cache();
            entry->owner->unlock();
            if (entry->next)
            {
                entry->next->prev = entry->prev;
            }
            else
            {
                used_entries_end_ = entry->prev;
            }

            if (entry->prev)
            {
                entry->prev->next = entry->next;
            }
            else
            {
                used_entries_start_ = entry->next;
            }

#if YETI_SANITY_CHECK
            if (used_entries_end_ && used_entries_end_->next)
            {
                raise(SanityCheckError, "cache alignment fail");
            }
            if (used_entries_start_ && used_entries_start_->prev)
            {
                raise(SanityCheckError, "cache alignment fail");
            }
#endif

            return entry;
        }
        entry = entry->next;
    }
    return 0;
}

bool
DataCache::pull(DataCacheEntry* entry)
{
    bool test = entry->trylock();
    if (!test) //no dice... cache entry is no good
        return false;

    lock();
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
        raise(SanityCheckError, "cache alignment fail");
    }
    if (used_entries_start_ && used_entries_start_->prev)
    {
        raise(SanityCheckError, "cache alignment fail");
    }
#endif

    unlock();

    entry->pulled = true;

    return true;
}

DataCacheEntry*
DataCache::pull(Cachable* owner)
{

    lock();
    DataCacheEntry* entry = find_free_entry();
    if (!entry) entry = find_used_entry();

    unlock();

    if (!entry)
        raise(SanityCheckError, "too many cache blocks allocated");

#if YETI_SANITY_CHECK
    if (!entry->is_locked())
        raise(SanityCheckError, "how exactly did you manage to return an unlocked cache entry");
#endif

    entry->assign(owner);
    entry->pulled = true;

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

uli
DataCache::nfree() const
{
    cerr << "not implemented" << endl;
    abort();
}

void
DataCache::insert(DataCacheEntry* entry)
{
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
    entry->unlock();

#if YETI_SANITY_CHECK
    if (used_entries_end_->next)
    {
        raise(SanityCheckError, "cache alignment fail");
    }
    if (used_entries_start_->prev)
    {
        raise(SanityCheckError, "cache alignment fail");
    }
#endif

    unlock();
}

void
DataCache::free(DataCacheEntry* entry)
{
#if YETI_SANITY_CHECK
    if (!entry->owner->is_locked())
        raise(SanityCheckError, "cannot free cache entry of unlocked owner");
#endif

    lock();

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
        raise(SanityCheckError, "cache alignment fail");
    }
    if (used_entries_start_ && used_entries_start_->prev)
    {
        raise(SanityCheckError, "cache alignment fail");
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
        raise(SanityCheckError, "cache alignment fail");
    }
    if (used_entries_start_ && used_entries_start_->prev)
    {
        raise(SanityCheckError, "cache alignment fail");
    }
#endif

    unlock();
}

void
DataCache::print(std::ostream &os) const
{
    os << Env::indent << "Data Cache " << endl;
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