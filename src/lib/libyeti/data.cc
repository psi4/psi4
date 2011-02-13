#include "data.h"
#include "cache.h"
#include "tile.h"
#include "index.h"
#include "sort.h"
#include "exception.h"
#include "env.h"
#include "dataimpl.h"
#include "malloc.h"
#include "runtime.h"
#include "thread.h"

#include <sys/fcntl.h>
#include <sys/stat.h>
#include <cerrno>
#include <errno.h>

using namespace yeti;
using namespace std;


DECLARE_SUBMALLOC(Data, DataTemplate<double>);
DECLARE_SUBMALLOC(DataBlock,SortedBlock);

#define DEBUG_CACHE 0

DataMode::DataMode()
    : flag(DataMode::read)
{
}

Data::Data(size_t n)
    : n_(n)
{
}

Data::~Data()
{
}

size_t
Data::n() const
{
    return n_;
}

void
Data::accumulate(const double *data, double scale, const SortPtr &sort)
{
    raise(SanityCheckError, "data block not configured to accumulate double");
}

void
Data::accumulate(const int *data, int scale, const SortPtr &sort)
{
    raise(SanityCheckError, "data block not configured to accumulate int");
}

void
Data::accumulate(const float *data, float scale, const SortPtr &sort)
{
    raise(SanityCheckError, "data block not configured to accumulate float");
}

void
Data::accumulate(const quad *data, quad scale, const SortPtr &sort)
{
    raise(SanityCheckError, "data block not configured to accumulate float");
}

void
Data::assign(const double *data, const SortPtr &sort)
{
    raise(SanityCheckError, "data block not configured to assign double");
}

void
Data::assign(const quad *data, const SortPtr &sort)
{
    raise(SanityCheckError, "data block not configured to assign double");
}

void
Data::assign(const float *data, const SortPtr &sort)
{
    raise(SanityCheckError, "data block not configured to assign float");
}

void
Data::assign(const int *data, const SortPtr &sort)
{
    raise(SanityCheckError, "data block not configured to assign int");
}

void
Data::assign(Data* data, const SortPtr &sort)
{
    raise(SanityCheckError, "data block not configured to assign int");
}

void
Data::reference(void* data)
{
    raise(SanityCheckError, "data block not configured to reference char");
}


DataBlock::DataBlock(
    uli n
)
    :
  n_(n),
  data_(0)
{
    yeti_register_new(lock_.get());
}

DataBlock::~DataBlock()
{
    if (data_)
        delete data_;
}

void
DataBlock::print(std::ostream &os) const
{
    ++Env::indent;
    if (data_)
        data_->print(os);
    --Env::indent;
}

void
DataBlock::sort(const SortPtr &, void *buffer)
{
    raise(SanityCheckError, "only memory and disk blocks can be sorted");
}

size_t
DataBlock::n() const
{
    return n_;
}

void
DataBlock::allocate(TemplateInfo::type_t type)
{
    if (data_ && data_->type() == type)
        return; //already allocated

    data_type_switch(type, allocate);
}

Data*
DataBlock::data() const
{
    return data_;
}

void
DataBlock::set_buffer(uli offset, const BufferPtr& buffer)
{
}

MemoryBlock::MemoryBlock(
    uli n
) : DataBlock(n)
{
}

MemoryBlock::~MemoryBlock()
{
}

void
MemoryBlock::_retrieve(uli threadnum)
{
}

void
MemoryBlock::sort(const SortPtr& sort, void *buffer)
{
    data_->sort(sort, buffer);
}

void
MemoryBlock::_release(uli threadnum)
{
}

void
MemoryBlock::print(std::ostream &os) const
{
    os << Env::indent << "Memory Block at ";
    if (data_)
        os << data_->pointer() << endl;
    else
        os << 0 << endl;

    DataBlock::print(os);
}

CachedDataBlock::CachedDataBlock(
    const DataMode* mode,
    uli n,
    const LayeredDataCachePtr& cache
) : DataBlock(n),
    main_cache_(cache),
    cache_entry_(0),
    cache_(0),
    mode_(mode) //no cache block until we allocate
{
}

void
CachedDataBlock::finalize()
{
    if (is_retrieved())
    {
        if (YetiRuntime::is_threaded_runtime() && cache_entry_->trylock())
        {
            cerr << "cache entry should not be unlocked" << endl;
            abort();
        }

        cache_entry_->clear(); //clear the data pointer
        //put the cache entry at the end
        cache_->insert(cache_entry_);
        cache_entry_->unlock();
    }
    else if (cache_entry_ && cache_entry_->trylock())
    {
        if (data_ && data_->nonnull())
            cache_entry_->clear();
        cache_entry_->unlock();
    }
}

CachedDataBlock::~CachedDataBlock()
{
}

void
CachedDataBlock::init_block(uli threadnum)
{
    switch(mode_->flag)
    {
        case DataMode::read:        init_read(threadnum); break;
        case DataMode::write:       init_write(threadnum); break;
        case DataMode::accumulate:  init_accumulate(threadnum); break;
    }
}

bool
CachedDataBlock::in_cache() const
{
    return (cache_entry_ && cache_entry_->owner == this);
}

void
CachedDataBlock::_retrieve(uli threadnum)
{
    uli offset = 0;
    if (cache_ == 0) //fetch a cache object
    {
        //first time this block has been retrieved
        main_cache_->lock();
        cache_ = main_cache_->allocate_cache(data_->size(), offset);
        if (!cache_)
        {
            cerr << "No appropriate cache block found for size " << data_->size() << endl;
            abort();
        }
        main_cache_->unlock();
        cache_entry_ = cache_->pull(this, offset);
        init_block(threadnum);
        return;
    }

    if (data_->null())
    {
            //do nothing... I have been cleared
    }
    else if (cache_entry_->owner == this
             && cache_entry_->trylock()) //attempt to grab ownership!
    {
        //at this point we can guarantee the data is non-null
        //and that we own the cache entry
        //because the data block is locked
        cache_->pull(cache_entry_);
        return;
    }
    else
    {
        //I am not yet null, but I no longer own my cache entry
        //the data block that now owns my cache entry
        //is waiting to unlock me to clear my data
        //this means that it is my responsibility
        //to clear the data and move on
        offset = cache_entry_->offset + 1;
        clear();
        cache_entry_ = 0;
    }

    //this returns a locked cache entry
    cache_entry_ = cache_->pull(this, offset);
    init_block(threadnum);
}

void
CachedDataBlock::clear()
{
    if (is_retrieved())
    {
        cerr << "Cannot clear retrieved data block" << endl;
        abort();
    }

    switch(mode_->flag)
    {
        case DataMode::read:        finalize_read(NOT_THREADED); break;
        case DataMode::write:       finalize_write(NOT_THREADED); break;
        case DataMode::accumulate:  finalize_accumulate(NOT_THREADED); break;
    }
    data_->clear();
}

void
CachedDataBlock::_release(uli threadnum)
{
    if (YetiRuntime::is_threaded_runtime() && cache_entry_->trylock())
    {
        cerr << "The cache entry should have been locked!" << endl;
        abort();
    }
    //put the block back in the cache and unlock it
    cache_->insert(cache_entry_);
    cache_entry_->unlock();
}

void
CachedDataBlock::print(std::ostream& os) const
{
    os << Env::indent << "Cache:" << cache_entry_ << " Data:" << data_->pointer() << std::endl;
    DataBlock::print(os);
}

void
CachedDataBlock::init_read(uli threadnum)
{
}

void
CachedDataBlock::init_write(uli threadnum)
{
}

void
CachedDataBlock::init_accumulate(uli threadnum)
{
}

void
CachedDataBlock::finalize_read(uli threadnum)
{
}

void
CachedDataBlock::finalize_write(uli threadnum)
{
}

void
CachedDataBlock::finalize_accumulate(uli threadnum)
{
}

RecomputedBlock::RecomputedBlock(
    const ThreadedTileElementComputerPtr& filler,
    Tile* tile
) : CachedDataBlock(tile->config()->get_data_mode(), tile->ndata(), tile->config()->cache),
    filler_(filler),
    tile_(tile)
{
    if (YetiRuntime::is_threaded_runtime())
    {
        cerr << "this should not be threaded" << endl;
        abort();
    }
}

RecomputedBlock::RecomputedBlock(
    uli ndata,
    const ThreadedTileElementComputerPtr& filler,
    const LayeredDataCachePtr& cache,
    Tile* tile
) : CachedDataBlock(tile->config()->get_data_mode(), ndata, cache),
    filler_(filler),
    tile_(tile)
{
    if (YetiRuntime::is_threaded_runtime())
    {
        cerr << "this should not be threaded" << endl;
        abort();
    }
}

RecomputedBlock::~RecomputedBlock()
{
    CachedDataBlock::finalize();
}


void
RecomputedBlock::init_read(uli threadnum)
{
    filler_->compute(tile_, data_, threadnum);
}

void
RecomputedBlock::init_write(uli threadnum)
{
    raise(SanityCheckError, "init_write should never be called on a recomputed block");
}

void
RecomputedBlock::init_accumulate(uli threadnum)
{
    raise(SanityCheckError, "init_accumulate should never be called on a recomputed block");
}

void
RecomputedBlock::print(std::ostream &os) const
{
    os << Env::indent << "Recomputed Block" << endl;
    CachedDataBlock::print(os);
}

SortedBlock::SortedBlock(
    const DataMode* mode,
    const DataBlockPtr& parent,
    const LayeredDataCachePtr& cache,
    const IndexRangeTuplePtr& tuple,
    Permutation** perms,
    uli nperms
) :
    CachedDataBlock(mode, parent->n(), cache),
    parent_(parent),
    sort_(new Sort(*perms, tuple)),
    perms_(perms),
    nperms_(nperms),
    tuple_(tuple)
{
    allocate(parent_->data()->type());
}

SortedBlock::~SortedBlock()
{
    CachedDataBlock::finalize();
}

uli
SortedBlock::nperms() const
{
    return nperms_;
}

Permutation**
SortedBlock::perms() const
{
    return perms_;
}

void
SortedBlock::print(std::ostream &os) const
{
    os << Env::indent << "Sorted Block" << endl;
    CachedDataBlock::print(os);
}

void
SortedBlock::init_read(uli threadnum)
{
    sort_->configure(perms_[0], tuple_);
    parent_->retrieve(threadnum);
    data_->assign(parent_->data(), sort_);
    parent_->release(threadnum);
}

void
SortedBlock::init_accumulate(uli threadnum)
{
    //set the values to zero so we only accumulate new contributions
    data_->memset();
}

void
SortedBlock::finalize_accumulate(uli threadnum)
{
    if (data_->null())
        return;

    parent_->retrieve(threadnum);
    //loop through all the permutations in the list
    for (uli i=0; i < nperms_; ++i)
    {
        sort_->configure(perms_[i], tuple_);
        parent_->data()->accumulate(data_, sort_->get_permutation()->sign(), sort_);
    }
    parent_->release(threadnum);
}

SubsetDataBlock::SubsetDataBlock(uli n)
    : DataBlock(n),
    parent_(0),
    offset_(0)
{
}

void
SubsetDataBlock::print(std::ostream &os) const
{
    os << Env::indent << "Subset Block at offset " << offset_
        << " for parent " << (void*) parent_->data()->buffer()
        << endl;
    DataBlock::print(os);
}

void
SubsetDataBlock::_release(uli threadnum)
{
    parent_->release(threadnum);
}

void
SubsetDataBlock::_retrieve(uli threadnum)
{
    parent_->retrieve(threadnum);

    char* dataptr = parent_->data()->buffer() + offset_;
    data_->reference(dataptr);
}

void
SubsetDataBlock::configure(
    const DataBlockPtr &parent,
    uli offset
)
{
    if (parent_)
    {
        cerr << "Subset data block already configured" << endl;
        abort();
    }

    parent_ = parent;
    offset_ = offset;
    allocate(parent_->data()->type());
}

void
SubsetDataBlock::sort(const SortPtr &sort, void *buffer)
{
    data_->sort(sort, buffer);
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
    uli offset,
    uli size,
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
    uli offset,
    uli size,
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

void
DiskBuffer::allocate_region(uli offset, uli size)
{

    uli end = size + offset;
    if (size_ >= end)
        return; //this space already exists in the buffer

    if (size_ != offset)
        raise(SanityCheckError, "misaligned buffer allocation");

    if (size == 0)
        raise(SanityCheckError, "cannot allocate region of size 0");

    ::lseek(fileno_, size - 1, SEEK_END);
    char dummy[] = { 0 };
    ::write(fileno_, dummy, 1);
    size_ += size;
}

void
DiskBuffer::_retrieve(uli threadnum)
{
}

void
DiskBuffer::_release(uli threadnum)
{
}


LocalDiskBlock::LocalDiskBlock(
    const DataMode* mode,
    uli n,
    const LayeredDataCachePtr& cache
)    : CachedDataBlock(mode, n, cache),
  buffer_(0),
  offset_(0)
{
}

LocalDiskBlock::~LocalDiskBlock()
{
    CachedDataBlock::finalize();
}

void
LocalDiskBlock::init_read(uli threadnum)
{
    buffer_->read(offset_, data_->size(), data_->buffer());
}

void
LocalDiskBlock::init_write(uli threadnum)
{
    init_read(threadnum);
}

void
LocalDiskBlock::init_accumulate(uli threadnum)
{
    init_read(threadnum);
}

void
LocalDiskBlock::finalize_read(uli threadnum)
{
    //do nothing
}

void
LocalDiskBlock::finalize_write(uli threadnum)
{
    if (data_->null())
        return;

    buffer_->write(offset_, data_->size(), data_->buffer());
}

void
LocalDiskBlock::finalize_accumulate(uli threadnum)
{
    if (data_->null())
        return;

    buffer_->write(offset_, data_->size(), data_->buffer());
}

void
LocalDiskBlock::set_buffer(uli offset, const BufferPtr& buffer)
{
    offset_ = offset;
    buffer_ = buffer;
    buffer_->allocate_region(offset_, data_->size());
}


#define MAX_NBLOCKS_ALLOC_QUEUE 10000

DataBlockFactory::DataBlockFactory()
    : filler_(0),
    blocks_(new DataBlock*[MAX_NBLOCKS_ALLOC_QUEUE]),
    nblocks_(0)
{
    filler_ = new ThreadedTileElementComputer(new MemsetElementComputer);
}

DataBlockFactory::DataBlockFactory(const TileElementComputerPtr &filler)
    :
    filler_(new ThreadedTileElementComputer(filler)),
    blocks_(new DataBlock*[MAX_NBLOCKS_ALLOC_QUEUE]),
    nblocks_(0)
{
}

void
DataBlockFactory::allocate_blocks()
{
    init_allocation();
    for (uli i=0; i < nblocks_; ++i)
    {
        allocate(blocks_[i]);
    }
    //all blocks have now been allocated!
    nblocks_ = 0;
}

void
DataBlockFactory::configure(Tensor *tensor)
{
    //do nothing by default
}

void
DataBlockFactory::register_allocation(DataBlock *dblock)
{
    blocks_[nblocks_] = dblock;
    ++nblocks_;

    if (nblocks_ == MAX_NBLOCKS_ALLOC_QUEUE)
        allocate_blocks();
}

DataBlockFactory::~DataBlockFactory()
{
    delete[] blocks_;
}

ThreadedTileElementComputer*
DataBlockFactory::get_element_computer() const
{
    return filler_.get();
}

SubsetBlockFactory::SubsetBlockFactory(
    const DataBlockFactoryPtr& parent
) :
    parent_(parent)
{
}

SubsetBlockFactory::~SubsetBlockFactory()
{
}

void
SubsetBlockFactory::allocate(DataBlock* block)
{
    //do nothing
}

DataBlock*
SubsetBlockFactory::get_block(Tile *tile)
{
    return new SubsetDataBlock(tile->ndata());
}

DataBlockFactory*
SubsetBlockFactory::copy() const
{
    SubsetBlockFactory* factory
        = new SubsetBlockFactory(parent_);
    return factory;
}

void
SubsetBlockFactory::init_allocation()
{
    //do nothing
}

Data::storage_t
SubsetBlockFactory::storage_type() const
{
    return parent_->storage_type();
}



