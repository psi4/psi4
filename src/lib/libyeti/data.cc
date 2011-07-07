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

using namespace yeti;
using namespace std;

DECLARE_PARENT_MALLOC(StorageBlock);
DECLARE_SUB_MALLOC(StorageBlock,InCoreBlock);
DECLARE_SUB_MALLOC(StorageBlock,LocalDiskBlock);
DECLARE_SUB_MALLOC(StorageBlock,CachedStorageBlock);

StorageBlock::StorageBlock(size_t size)
    :
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
}

void
InCoreBlock::release()
{
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
    raise(SanityCheckError, "in core block should never be committed");
}

void
InCoreBlock::clear(Cachable* owner)
{
#if MALLOC_CACHE_BLOCK
    if (data_)
    {
        YetiRuntime::free(data_, size_);
        data_ = 0;
    }
#else
    raise(SanityCheckError, "in core block should never be cleared");
#endif
}

void
InCoreBlock::obsolete(Cachable* item)
{
#if MALLOC_CACHE_BLOCK
    clear(item);
#else
    raise(SanityCheckError, "in core block should never be obsolete");
#endif
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
}

CachedStorageBlock::~CachedStorageBlock()
{
    clear(cache_item_);
}

DataCacheEntry*
CachedStorageBlock::get_entry() const
{
    return cache_entry_;
}

bool
CachedStorageBlock::retrieve()
{
    if (data_) //yay! we are in cache
    {
        bool test = cache_entry_->trylock();
        if (test && data_)
        {
            cache_->pull(cache_entry_);
            return true;
        }
    }

    //this returns a locked cache entry
    cache_entry_ = cache_->pull(this);
    fetch_data();
    return false;

}

void
CachedStorageBlock::release()
{
    cache_->insert(cache_entry_);
    cache_entry_->unlock();
}

void
CachedStorageBlock::clear(Cachable* owner)
{
    if (cache_entry_ && cache_entry_->trylock())
    {
#if DEBUG_CACHE_USAGE
        TensorBlock* tblock = dynamic_cast<TensorBlock*>(cache_entry_->owner);
        if (tblock)
        {
            dout << "clearing cache entry " << cache_entry_->offset
                    << " on block " << tblock->get_index_string()
                    << " on tensor " << tblock->get_parent_tensor()->get_name()
                    << " at location " << (void*) tblock
                    << endl;
        }
#endif

        if (cache_entry_->owner == owner)
            cache_entry_->owner = 0;

        cache_entry_->unlock();
    }

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
    raise(SanityCheckError, "generic cached block should never be committed");
}

void
CachedStorageBlock::fetch_data()
{
}

void
CachedStorageBlock::free_cache_entry()
{
    cache_entry_->owner = 0;
    release();
}

Cachable*
CachedStorageBlock::get_cache_item() const
{
    return cache_item_;
}

void
CachedStorageBlock::obsolete(Cachable* item)
{
    clear(item);
    data_ = 0;
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
        raise(SanityCheckError, "cannot allocate region of size 0");

    ::lseek(fileno_, size - 1, SEEK_END);
    char dummy[] = { 0 };
    ::write(fileno_, dummy, 1);
    size_ += size;

    return offset;
}


#ifdef redefine_size_t
#define size_t custom_size_t
#endif

