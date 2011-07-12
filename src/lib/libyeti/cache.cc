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

namespace yeti {
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define MALLOC_CACHE_ENTRIES 0

DECLARE_MALLOC(DataCacheEntry);
DECLARE_MALLOC(DataCache);

DataCacheEntry::DataCacheEntry(
    char* d,
    uli posoffset
)
    :
    data(d),
    owner(0),
    offset(posoffset),
    pulled(false)
{
}

DataCacheEntry::~DataCacheEntry()
{
}

void
DataCacheEntry::assign(Cachable* new_owner)
{
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
    os << stream_printf("This:%12p Offset:%ld Owner:%12p Data:%12p", this, offset, owner, data);
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
    storage_(0),
    offsets_(new uli[YetiRuntime::nthread()]),
    starts_(new uli[YetiRuntime::nthread()]),
    stops_(new uli[YetiRuntime::nthread()])
{
    data_ = YetiRuntime::malloc(storage);

    append(data_, storage);
}

DataCache::DataCache(size_t blocksize)
    :
    nentries_(0),
    blocksize_(blocksize),
    data_(0),
    storage_(0),
  offsets_(new uli[YetiRuntime::nthread()]),
  starts_(new uli[YetiRuntime::nthread()]),
  stops_(new uli[YetiRuntime::nthread()])
{
}

DataCache::~DataCache()
{
    for (uli i=0; i < nentries_; ++i)
    {
        delete entries_[i];
    }
    delete[] offsets_;
    delete[] starts_;
    delete[] stops_;

    YetiRuntime::free(data_, storage_);
}

void
DataCache::append(
    char *data,
    size_t storage
)
{
    if (data != data_)
    {
        raise(SanityCheckError, "cache not configured for multiple appends");
    }

    char* ptr = data;

    uli nentries_new = storage / blocksize_;

    uli index_start = nentries_;
    uli index_stop = nentries_ + nentries_new;

    if (index_stop > NCACHE_ENTRIES_MAX)
    {
        cerr << "Insufficient storage space in cache for all the cache blocks." << endl;
        cerr << "Need space for " << index_stop << " Have " << NCACHE_ENTRIES_MAX << endl;
        abort();
    }

    for (uli i=index_start; i < index_stop; ++i, ptr += blocksize_)
    {
        uli offset = i;
        DataCacheEntry* entry = new DataCacheEntry(ptr, offset);
        entries_[i] = entry;
        all_entries_[i] = entry;
        ++nentries_;
    }

    uli nthread = YetiRuntime::nthread();
    for (uli i=0; i < nthread; ++i)
    {
        uli start = (nentries_ * i) / nthread;
        uli stop = (nentries_ * (i + 1)) / nthread;
        offsets_[i] = start;
        starts_[i] = start;
        stops_[i] = stop;
    }

    storage_ += storage;
}


void
DataCache::flush()
{
    for (uli i=0; i < nentries_; ++i)
    {
        entries_[i]->flush();
    }
}

DataCacheEntry*
DataCache::find_entry(uli offset, uli stop, bool check_occupied)
{
    stop = stop > nentries_ ? nentries_ : stop;
    for (uli i=offset; i < stop; ++i)
    {
        DataCacheEntry* entry = entries_[i];
        if (entry)
        {
#if DEBUG_CACHE_USAGE
            dout << "checking cache entry " << entry->offset << endl;
#endif
            Cachable* owner = entry->get_owner();
            if (!owner)
            {
                if (entry->trylock())
                    return entry;
            }
            else
            {
                if (check_occupied && owner->trylock())
                {
                    if (owner->is_retrieved())
                    {
                        owner->unlock();
                        continue;
                    }
                    else if (owner->in_destructor())
                    {
                        cerr << "owner in destructor!" << endl;
                        abort();
                    }

                    bool got_lock = entry->trylock();
                    if (!got_lock)
                    {
                        raise(SanityCheckError, "cache locked owner but not entry");
                    }
                    return entry;
                }
            }
        }
    }
    return 0;
}

DataCacheEntry*
DataCache::pull_unused(uli threadnum)
{    
    DataCacheEntry* entry = find_entry(offsets_[threadnum], stops_[threadnum], false);
    if (!entry)
        entry = find_entry(starts_[threadnum], offsets_[threadnum], false);

    return entry;
}

DataCacheEntry*
DataCache::pull_any(uli threadnum)
{
    DataCacheEntry* entry = find_entry(offsets_[threadnum], stops_[threadnum], true);
    if (!entry)
        entry = find_entry(starts_[threadnum], offsets_[threadnum], true);
    return entry;
}

void
DataCache::print_details(std::ostream& os)
{
    for (uli i=0; i < nentries_; ++i)
    {
        DataCacheEntry* entry = all_entries_[i];
        os << "Cache entry " << i
             << " at pointer " << (void*) entry;
        if (!entry->is_pulled())
        {
            os << " is not pulled" << endl;
        }
        else
        {
            TensorBlock* block = dynamic_cast<TensorBlock*>(entry->get_owner());
            if (block)
            {
                 os << " pulled by block "
                    << block->get_index_string() << " on tensor "
                    << block->get_parent_tensor()->get_name() << endl;
            }
            else
               os << " pulled by something" << endl;
        }
    }
}

void
DataCache::process_entry(
    DataCacheEntry* entry,
    CachedStorageBlock* block,
    uli threadnum
)
{
    if (!entry)
    {
        uli nretrieved = 0;
        uli npulled = 0;
        uli nfree = 0;
        uli threadnum = YetiRuntime::get_thread_number();
        uli start = starts_[threadnum];
        uli stop = stops_[threadnum];
        for (uli i=start; i < stop; ++i)
        {
            DataCacheEntry* entry = entries_[i];
            if (!entry)
            {
                ++npulled;
            }
            else if (entry->trylock())
            {
                ++nfree;
                entry->unlock();
            }
            else
            {
                Cachable* owner = entry->owner;
                if (owner && owner->is_retrieved())
                    ++nretrieved;
            }
        }

        //loop the entire cache and see how many blocks are retrieve
        print_details(cerr);

        uli nblocks = stop - start;

        std::string errmsg = stream_printf(
            "Too many active cache blocks for thread %d!\n"
            "Total cache size is %d bytes.\n"
            "Cache block size is %d bytes.\n"
            "Total number of cache entries for thread is %d.\n"
            "Number of pulled cache entries is %d.\n"
            "Number of free cache entries is %d.\n"
            "Number of retrieved cache entries is %d.",
            threadnum,
            storage_, blocksize_,
            nblocks,
            npulled,
            nfree,
            nretrieved
        );
        raise(SanityCheckError, errmsg);
    }

    //if the entry has an existing cache entry, it will be locked
    Cachable* old_owner = entry->get_owner();

#if DEBUG_CACHE_USAGE
    TensorBlock* tblock = dynamic_cast<TensorBlock*>(old_owner);
    if (tblock)
    {
        dout << "flushing cache entry " << entry->offset
           << " on old owner " << tblock->get_index_string()
           << " at location " << (void*) tblock
           << endl;
    }
#endif

    if (old_owner)
    {
        old_owner->flush_from_cache();
        old_owner->unlock();
    }

    Cachable* new_owner = block->get_cache_item();
    entry->assign(new_owner);
    entries_[entry->offset] = 0;
    block->assign_location(entry->data);

    offsets_[threadnum] = entry->offset + 1;
}

DataCacheEntry*
DataCache::pull(CachedStorageBlock* block)
{
    if (block->size() > blocksize_)
    {
        raise(SanityCheckError,"Requested block size is too large!");
    }

    uli threadnum = YetiRuntime::get_thread_number();

    //look for unoccupied entries
    DataCacheEntry* entry = pull_unused(threadnum);
    if (!entry)
        entry = pull_any(threadnum);

    process_entry(entry, block, threadnum);
    entry->pulled = true;
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
    if (!entry->owner)
    {
        raise(SanityCheckError, "cache entry is getting inserted without an owner!");
    }
    entries_[entry->offset] = entry;
    entry->pulled = false;
}


void
DataCache::pull(DataCacheEntry* entry)
{
    entries_[entry->offset] = 0;
    entry->pulled = true;
}

void
DataCache::print(std::ostream &os) const
{
    os << Env::indent << "Data Cache "
        << nentries_ << " blocks of size " << blocksize_;
    for (uli i=0; i < nentries_; ++i)
    {
        DataCacheEntry* entry = entries_[i];
        os << endl << Env::indent;
        entry->print(os);
    }
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

} // end namespace yeti

