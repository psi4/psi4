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
    const DataMode* mode,
    uli n
)
    :
  n_(n),
  data_(0),
  mode_(mode)
{
    yeti_register_new(lock_.get());
}

DataBlock::~DataBlock()
{
    if (data_)
        delete data_;
}

void
DataBlock::_retrieve(uli threadnum)
{
    switch(mode_->flag)
    {
        case DataMode::read:        init_read(threadnum); break;
        case DataMode::write:       init_write(threadnum); break;
        case DataMode::accumulate:  init_accumulate(threadnum); break;
    }
}

void
DataBlock::_release(uli threadnum)
{
    switch(mode_->flag)
    {
        case DataMode::read:        finalize_read(threadnum); break;
        case DataMode::write:       finalize_write(threadnum); break;
        case DataMode::accumulate:  finalize_accumulate(threadnum); break;
    }
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
DataBlock::init_read(uli threadnum)
{
}

void
DataBlock::init_write(uli threadnum)
{
}

void
DataBlock::init_accumulate(uli threadnum)
{
}

void
DataBlock::finalize_read(uli threadnum)
{
}

void
DataBlock::finalize_write(uli threadnum)
{
}

void
DataBlock::finalize_accumulate(uli threadnum)
{
}

void
DataBlock::set_buffer(uli offset, const BufferPtr& buffer)
{
}

MemoryBlock::MemoryBlock(
    const DataMode* mode,
    uli n
) : DataBlock(mode, n)
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
) : DataBlock(mode, n),
    main_cache_(cache),
    cache_entry_(0),
    cache_(0) //no cache block until we allocate
{
}

void
CachedDataBlock::finalize()
{
    if (is_retrieved())
    {
        //put the cache entry at the end
        cache_->push_back(cache_entry_);
        cache_entry_->clear(); //clear the data pointer
    }
    else if (data_ && data_->nonnull())
    {
        //clear the ownership and make sure block is on front of queue
        cache_->pull(cache_entry_);
        cache_->push_back(cache_entry_);
        cache_entry_->clear(); //make sure data is cleared
    }
}

CachedDataBlock::~CachedDataBlock()
{
}

void
CachedDataBlock::_retrieve(uli threadnum)
{
    if (cache_ == 0) //fetch a cache object and cache entry
    {
        cache_ = main_cache_->allocate_cache(data_->size());
        if (!cache_)
        {
            cerr << "No appropriate cache block found for size " << data_->size() << endl;
            abort();
        }
        cache_entry_ = cache_->pull(this);
        DataBlock::_retrieve(threadnum);
    }
    else if (data_->nonnull()) //I have a cache entry
    {
        //pull it out so we don't lose the data
        cache_->pull(cache_entry_);
    }
    else //fetch a new cache block from wherever it need be
    {
        cache_entry_ = cache_->pull(this);
        if (!cache_entry_)
        {
            std::cerr << "Cache size not large enough for computation. Need more memory." << std::endl;
            abort();
        }
        DataBlock::_retrieve(threadnum);
    }
}

void
CachedDataBlock::clear()
{
    DataBlock::_release(NOT_THREADED);
    data_->clear();
}

void
CachedDataBlock::_release(uli threadnum)
{
    //do not yet release here! only when cache blocks get cleared should
    //release be called

    //put on the front so its the last to be reused
    cache_->push_front(cache_entry_);
}

void
CachedDataBlock::print(std::ostream& os) const
{
    os << Env::indent << "Cache:" << cache_entry_ << " Data:" << data_->pointer() << std::endl;
    DataBlock::print(os);
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
    lseek(fileno_, offset, SEEK_SET);
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

void
DataBlockFactory::allocate(DataBlock* block)
{
    //do nothing
}

DataBlockFactory::~DataBlockFactory()
{
}

