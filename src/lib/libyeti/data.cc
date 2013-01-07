#include "data.h"
#include "tensorblock.h"
#include "tensorbranch.h"
#include "node.h"
#include "index.h"
#include "sort.h"
#include "exception.h"
#include "env.h"
#include "dataimpl.h"
#include "malloc.h"
#include "runtime.h"
#include "thread.h"
#include "cache.h"
#include "tensor.h"

#include <sys/fcntl.h>
#include <sys/stat.h>
#include <cerrno>
#include <errno.h>
#include <cstdio>

using namespace yeti;
using namespace std;

#define name(x) static_cast<TensorBlock*>(x)->get_block_name()

DECLARE_PARENT_MALLOC(StorageBlock);
DECLARE_SUB_MALLOC(StorageBlock,InCoreBlock);
DECLARE_SUB_MALLOC(StorageBlock,LocalDiskBlock);
DECLARE_SUB_MALLOC(StorageBlock,CachedStorageBlock);

StorageBlock::StorageBlock(size_t size)
    :
    retrieved_(false),
    size_(size),
    data_(0)
{
}

size_t
StorageBlock::size() const
{
    return size_;
}

char*
StorageBlock::data() const
{
    return data_;
}

void
StorageBlock::memset()
{
    ::memset(data_, 0, size_);
}

StorageBlock::~StorageBlock()
{
}


InCoreBlock::InCoreBlock(char* data, size_t size)
    : StorageBlock(size)
{
    data_ = data;
    retrieved_ = true;
}

InCoreBlock::~InCoreBlock()
{
}

void
InCoreBlock::release()
{
}

bool
InCoreBlock::is_cached() const
{
    return false;
}

bool
InCoreBlock::retrieve()
{
#if MALLOC_CACHE_BLOCK
    if (data_ == 0)
    {
        data_= YetiRuntime::malloc(size_);
        return false;
    }
#else
    return true;
#endif
}

void
InCoreBlock::commit()
{
    yeti_throw(SanityCheckError, "in core block should never be committed");
}

void
InCoreBlock::clear()
{
#if MALLOC_CACHE_BLOCK
    if (data_)
    {
        YetiRuntime::free(data_, size_);
        data_ = 0;
    }
#else
    yeti_throw(SanityCheckError, "in core block should never be cleared");
#endif
}

bool
InCoreBlock::is_retrieved() const
{
    return true;
}

CachedStorageBlock::CachedStorageBlock(
    Cachable* cache_item,
    DataCache* cache
)
    :
    StorageBlock(cache->blocksize()),
    cache_item_(cache_item),
    cache_(cache),
    cache_entry_(0)
{
    data_ = 0;
    retrieved_ = false;
    if (size_ > cache_->blocksize())
        yeti_throw(SanityCheckError, "cache is not large enough to hold block");
}

CachedStorageBlock::~CachedStorageBlock()
{
    clear();
}

DataCacheEntry*
CachedStorageBlock::get_entry() const
{
    return cache_entry_;
}

bool
CachedStorageBlock::retrieve()
{
    if (retrieved_)
        return true;

#if YETI_SANITY_CHECK
    if (!cache_item_->is_locked() && !cache_item_->is_initialized())
        yeti_throw(SanityCheckError, "cache item is not locked in retrieve");
#endif

    if (data_) //we might be in cache...
    {
        DataCacheEntry* new_entry = cache_->pull(cache_entry_, cache_item_);
        retrieved_ = true;
        if (new_entry == cache_entry_)
        {
            return true;
        }
        else
        {
            cache_entry_ = new_entry;
            return false;
        }
    }

    //this returns a locked cache entry
    cache_entry_ = cache_->pull(cache_item_);
    data_ = cache_entry_->data;

    retrieved_ = true;
    return false;
}

void
CachedStorageBlock::release()
{
    if (!retrieved_)
        return;

    if (!cache_entry_->pulled)
    {
#if DEBUG_CACHE_HISTORY
        cache_->print_history();
        TensorBlock* block = static_cast<TensorBlock*>(cache_item_);
        cerr << stream_printf("Failure releasing %s\n", block->get_block_name().c_str());
        for (uli i=0; i < NHISTORIES_CACHE; ++i)
        {
            if (i == cache_entry_->nentries_history_)
                cout << "-------------------" << endl;
            CacheHistory& h = cache_entry_->histories_[i];
            TensorBlock* block = static_cast<TensorBlock*>(h.owner);
            if (block)
                cout << stream_printf("Entry pulled by thread %ld by block %s\n", h.thread_number, block->get_block_name().c_str());
        }
        cout.flush();
#endif
        yeti_throw(SanityCheckError, "Cannot insert unpulled cache entry");
    }

    cache_->insert(cache_entry_);
    retrieved_ = false;
}

bool
CachedStorageBlock::is_cached() const
{
    return data_;
}

bool
CachedStorageBlock::is_retrieved() const
{
    return retrieved_;
}

void
CachedStorageBlock::clear()
{
#if YETI_SANITY_CHECK
    if (!cache_item_->is_locked())
    {
        TensorBlock* block = static_cast<TensorBlock*>(cache_item_);
        block->controller_fail();
        yeti_throw(SanityCheckError, "Storage block is neither locked nor retrieved in clear");
    }
    if (cache_item_->is_retrieved())
        yeti_throw(SanityCheckError, "Cannot clear retrieved item");
#endif
    if (retrieved_)
    {
        cerr << stream_printf("Cannot clear retrieved cache storage block on node %d\n", YetiRuntime::me());
        TensorBlock* block = static_cast<TensorBlock*>(cache_item_);
        block->controller_fail();
        throw TENSOR_BLOCK_POLICY_EXCEPTION;
    }
    retrieved_ = false;

    if (cache_entry_) //I am responsible for clearing this
        cache_->free(cache_entry_, cache_item_);

    cache_entry_ = 0;
    data_ = 0;
}

void
CachedStorageBlock::assign_location(char* ptr)
{
    data_ = ptr;
}

void
CachedStorageBlock::commit()
{
    yeti_throw(SanityCheckError, "generic cached block should never be committed");
}

Cachable*
CachedStorageBlock::get_cache_item() const
{
    return cache_item_;
}

LocalDiskBlock::LocalDiskBlock(
    Cachable* parent,
    DataCache* cache,
    const DiskBufferPtr& buffer
) :
    CachedStorageBlock(parent, cache),
  buffer_(buffer),
  offset_(0)
{
    buffer->lock();
    offset_ = buffer->allocate_region(size_);
    buffer->unlock();
}

void
LocalDiskBlock::fetch_data()
{
    //read into data buffer
    buffer_->read(offset_, size_, data_);
}

void
LocalDiskBlock::commit()
{
    buffer_->write(offset_, size_, data_);
}


DiskBuffer::DiskBuffer(const std::string &filename)
    :
    filename_(filename),
    size_(0)
{

    fileno_ = open(filename_.c_str(),O_RDWR,0644);
    if (fileno_ == -1)
    {
        fileno_ = creat(filename_.c_str(), 0644);
    }
    else
    {
        struct stat buf;
        fstat(fileno_, &buf);
        size_ = buf.st_size;
    }

    if (fileno_ == -1)
    {
        cerr << "Could not create file " << filename_ << endl;
        abort();
    }

    //close and reopen to keep posix from complaining
    //about the existence of the file
    close(fileno_);
    fileno_ = open(filename_.c_str(),O_RDWR,0644);
    if (fileno_ == -1)
    {
        cerr << "Could not create file " << filename_ << endl;
        abort();
    }

}

DiskBuffer::~DiskBuffer()
{
    close(fileno_);
    ::remove(filename_.c_str());
}

uli
DiskBuffer::size() const
{
    return size_;
}

void
DiskBuffer::read(
    size_t offset,
    size_t size,
    char* data
)
{
    off_t pos = ::lseek(fileno_, offset, SEEK_SET);

    if (pos != offset)
    {
        cerr << "failed to seek file position" << endl;
        abort();
    }

    uli amt = ::read(fileno_, data, size);
    if (amt != size)
    {
        cerr << "failed to read correct amount in buffer" << endl;
        abort();
    }
}

void
DiskBuffer::write(
    size_t offset,
    size_t size,
    const char* data
)
{
    lseek(fileno_, offset, SEEK_SET);
    uli amt = ::write(fileno_, data, size);
    if (amt != size)
    {
        cerr << "failed to write correct amount in buffer" << endl;
        abort();
    }
}

size_t
DiskBuffer::allocate_region(size_t size)
{
    size_t offset = size_;
    if (size == 0)
        yeti_throw(SanityCheckError, "cannot allocate region of size 0");

    ::lseek(fileno_, size - 1, SEEK_END);
    char dummy[] = { 0 };
    ::write(fileno_, dummy, 1);
    size_ += size;

    return offset;
}


#ifdef redefine_size_t
#define size_t custom_size_t
#endif

