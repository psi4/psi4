#include "runtime.h"
#include "thread.h"
#include "index.h"
#include "malloc.h"
#include "contraction.h"
#include "cache.h"

#include <libsmartptr/strop.h>
#include "data.h"
#include "tile.h"
#include "matrix.h"
#include "sort.h"


using namespace yeti;
using namespace std;

std::map<std::string, IndexRangePtr> YetiRuntime::ranges_;
size_t* YetiRuntime::max_range_sizes_ = 0;
size_t* YetiRuntime::max_block_sizes_ = 0;
usi YetiRuntime::max_depth_ = 0;
usi YetiRuntime::max_nindex_ = 0;
uli YetiRuntime::nthread_ = 0;
uli YetiRuntime::nproc_ = 0;
uli YetiRuntime::me_ = 0;
uli YetiRuntime::memory_ = 0;
uli YetiRuntime::cache_mem_ = 0;
uli YetiRuntime::remaining_mem_ = 0;
float YetiRuntime::matrix_multiply_cutoff = -12;
ThreadGroup* YetiRuntime::threadgrp_ = 0;
ThreadLock* YetiRuntime::printlock_ = 0;
ThreadLock* YetiRuntime::malloc_lock_ = 0;
bool YetiRuntime::threaded_runtime_ = false;

static FastMalloc matrix_generator_set_malloc("matrix generator set");
static FastMalloc matrix_generator_malloc("matrix generator");
static FastMalloc indexset_malloc("index set");
static FastMalloc tile_location_malloc("tile location");
static MallocOverride use_heap;

IndexRangeTuple*
YetiRuntime::get_index_tuple(const std::string& str)
{
    if (str.size() == 0)
    {
        return new IndexRangeTuple(0);
    }

    vector<std::string> str_tuple(0);
    smartptr::split(str_tuple, str, ", ");

    IndexRangeTuple* tuple
        = new IndexRangeTuple(str_tuple.size());

    usi i=0;
    vector<std::string>::const_iterator it(str_tuple.begin());
    vector<std::string>::const_iterator stop(str_tuple.end());
    for ( ; it != stop; ++it, ++i)
    {
        IndexRange* range = get_range( (*it).c_str() );
        tuple->set(i, range);
    }
    return tuple;
}

IndexRange*
YetiRuntime::get_range(const std::string& id)
{
    map<std::string, IndexRangePtr>::iterator it
        = ranges_.find(id);

    if (it == ranges_.end())
    {
        cerr << "Cannot find index range " << id << endl;
        abort();
    }

    return it->second.get();
}

void
YetiRuntime::init(
    uli me,
    uli nproc,
    usi nindex,
    uli nthread,
    uli memory,
    uli cache
)
{
    //the maximum number of indices to be allocated
    max_nindex_ = nindex;
    nproc_ = nproc;
    me_ = me;
    memory_ = memory;
    remaining_mem_ = memory;
    cache_mem_ = cache;

    init_sizes();
    init_malloc();
    init_threads(nthread);

    GlobalQueue::init();
}

void
YetiRuntime::init_sizes()
{
    if (ranges_.size() == 0)
    {
        cerr << "YetiRuntime has no index ranges" << endl;
        abort();
    }

    //figure out the maximum depth
    max_depth_ = 0;
    std::map<std::string, IndexRangePtr>::iterator start(ranges_.begin());
    std::map<std::string, IndexRangePtr>::iterator stop(ranges_.end());
    std::map<std::string, IndexRangePtr>::iterator it;
    for (it=start; it != stop; ++it)
    {
        const IndexRangePtr& range(it->second);
        usi depth = range->depth();
        if (depth > max_depth_)
            max_depth_ = depth;
    }

    //now loop the index ranges and determine the maximum range sizes at each level
    max_range_sizes_ = new uli[max_depth_ + 1];
    for (usi depth=0; depth <= max_depth_; ++depth)
    {
        uli maxsize = 0;
        for (it=start; it != stop; ++it)
        {
            const IndexRangePtr& range = it->second;
            uli size = range->nmax(depth);
            if (size > maxsize)
                maxsize = size;
        }
        max_range_sizes_[depth] = maxsize;
    }

    //now create the max block sizes
    max_block_sizes_ = new uli[max_depth_ + 1];
    for (usi depth=0; depth <= max_depth_; ++depth)
    {
        uli max_block_size = 1;
        uli max_range_size = max_range_sizes_[depth];
        for (usi idx=0; idx < max_nindex_; ++idx)
        {
            max_block_size *= max_range_size;
        }
        max_block_sizes_[depth] = max_block_size;
    }
}

void
YetiRuntime::init_malloc()
{
    indexset_malloc.init(indexset_size, nindexset);
    IndexRangeTuple::init_malloc(NMALLOC_BLOCKS_DEBUG);

    #define NINDEXERS 100000
    Indexer::init_malloc(NINDEXERS);

    Tile::init_malloc(NMALLOC_BLOCKS_DEBUG);
    TileMap::init_malloc(NMALLOC_BLOCKS_DEBUG / 5);
    Data::init_malloc(NMALLOC_BLOCKS_DEBUG);
    DataBlock::init_malloc(NMALLOC_BLOCKS_DEBUG);

    #define NCXN_TASKS_MALLOC 10000000
    Task::init_malloc(NCXN_TASKS_MALLOC);

    Matrix::init_malloc(NMALLOC_BLOCKS_DEBUG);
    MatrixMap::init_malloc(NMALLOC_BLOCKS_DEBUG / 5);

    Sort::init_malloc(NMALLOC_BLOCKS_DEBUG / 2);

    #define NTHREAD_BLOCKS_MALLOC 10000000
    ThreadLock::init_malloc(NTHREAD_BLOCKS_MALLOC);

    #define NCACHE_BLOCKS_MALLOC 100000
    DataCacheEntry::init_malloc(NCACHE_BLOCKS_MALLOC);


    #define PERMSET_SIZE 40
    #define NPERMSETS 100000
    matrix_generator_set_malloc.init(PERMSET_SIZE, NPERMSETS);

    #define NMATRIX_GENERATORS 10000000
    matrix_generator_malloc.init(sizeof(void*), NMATRIX_GENERATORS);

    #define NTILE_LOCATIONS 1000
    tile_location_malloc.init(sizeof(uli) * MAX_DEPTH, NTILE_LOCATIONS);

    //make sure all mallocs have been configured
#if USE_MALLOC_CLASS
    FastMalloc::check_queue();
#endif
}

bool
YetiRuntime::is_threaded_runtime()
{
    return threaded_runtime_;
}

void
YetiRuntime::set_threaded_runtime(bool flag)
{
    threaded_runtime_ = flag;
}

void
YetiRuntime::init_threads(uli nthread)
{
    if (threadgrp_) //overwriting previous thread initialization
    {
        ThreadEnvironment::deallocate();
        /**
        These were not allocated from the custom allocator classes
        so calling delete in the middle of runtime is not defined.
        Just leak the memory for now.
        delete printlock_;
        delete malloc_lock_;
        delete threadgrp_;
        */
    }
    else
    {
        printlock_ = new (use_heap) pThreadLock;
        malloc_lock_ = new (use_heap) pThreadLock;
    }

    threadgrp_ = new pThreadGroup(nthread);
    nthread_ = nthread;

    ThreadEnvironment::allocate();
}

void*
YetiRuntime::allocate(size_t n)
{
    if (n > remaining_mem_)
    {
        cerr << "insufficient memory for allocation" << endl;
        abort();
    }

    return ::malloc(n);
    remaining_mem_ -= n;
}

void
YetiRuntime::free(void* ptr, size_t n)
{
    ::free(ptr);
    remaining_mem_ += n;
}

usi
YetiRuntime::max_depth()
{
    return max_depth_;
}

usi
YetiRuntime::max_nindex()
{
    return max_nindex_;
}

const size_t*
YetiRuntime::max_block_sizes()
{
    return max_block_sizes_;
}

const size_t*
YetiRuntime::max_range_sizes()
{
    return max_range_sizes_;
}

uli
YetiRuntime::me()
{
    return me_;
}

uli
YetiRuntime::nproc()
{
    return nproc_;
}

uli
YetiRuntime::nthread()
{
    return nthread_;
}

void
YetiRuntime::reconfigure(uli me, uli nproc)
{
    me_ = me;
    nproc_ = nproc;
}

void
YetiRuntime::register_index_range(
    const std::string& id,
    IndexRange *range
)
{
    ranges_[id] = range;
    range->incref(); //make permanent for computation
}

void
YetiRuntime::register_index_range(
    IndexRange *range,
    const char *id1,
    const char *id2,
    const char *id3,
    const char *id4,
    const char *id5,
    const char *id6,
    const char *id7
)
{
    register_index_range(id1, range);

    if (id2) register_index_range(id2, range);
    if (id3) register_index_range(id3, range);
    if (id4) register_index_range(id4, range);
    if (id5) register_index_range(id5, range);
    if (id6) register_index_range(id6, range);
    if (id7) register_index_range(id7, range);

}

ThreadGroup*
YetiRuntime::get_thread_grp()
{
    return threadgrp_;
}

void
YetiRuntime::lock_print()
{
    printlock_->lock();
}

void
YetiRuntime::unlock_print()
{
    printlock_->unlock();
}

void
YetiRuntime::lock_malloc()
{
    if (malloc_lock_)
        malloc_lock_->lock();
}

void
YetiRuntime::unlock_malloc()
{
    if (malloc_lock_)
        malloc_lock_->unlock();
}


template <class data_t>
void
array_allocate(data_t* dst, data_t i, data_t j)
{
    dst[0] = i;
    dst[1] = j;
}

usi*
yeti::usi_allocate(usi i, usi j)
{
    usi* dst = yeti_malloc_perm();
    array_allocate<usi>(dst, i,j);
    return dst;
}

uli*
yeti::yeti_malloc_indexset()
{
    return reinterpret_cast<uli*>(indexset_malloc.malloc());
}

void
yeti::yeti_free_indexset(const uli* ptr)
{
    indexset_malloc.free(ptr);
}

void*
yeti::yeti_malloc_indexptr()
{
    return indexset_malloc.malloc();
}

void
yeti::yeti_free_indexptr(void *ptr)
{
    indexset_malloc.free(ptr);
}

usi*
yeti::yeti_malloc_perm()
{
    return reinterpret_cast<usi*>(indexset_malloc.malloc());
}

void
yeti::yeti_free_perm(usi* ptr)
{
    indexset_malloc.free(ptr);
}

void*
yeti::yeti_malloc_matrix_generator()
{
    return matrix_generator_malloc.malloc();
}

void
yeti::yeti_free_matrix_generator(void *ptr)
{
    matrix_generator_malloc.free(ptr);
}

void*
yeti::yeti_malloc_matrix_generator_set()
{
    return matrix_generator_set_malloc.malloc();
}

void
yeti::yeti_free_matrix_generator_set(void *ptr)
{
    matrix_generator_set_malloc.free(ptr);
}

uli*
yeti::yeti_malloc_tile_location()
{
    return reinterpret_cast<uli*>(tile_location_malloc.malloc());
}

void
yeti::yeti_free_tile_location(void *ptr)
{
    tile_location_malloc.free(ptr);
}

void
YetiRuntime::register_subrange(
    IndexRange *parent,
    IndexRange *child
)
{
}

void
YetiRuntime::finalize()
{
    GlobalQueue::finalize();
    ThreadEnvironment::deallocate();
    if (printlock_) //these came off standard malloc heap
        ::free(printlock_);
    if (malloc_lock_)
        ::free(malloc_lock_);
    if (threadgrp_) //this came from operator new
        delete threadgrp_;
}



