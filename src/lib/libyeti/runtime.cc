#include "runtime.h"
#include "thread.h"
#include "index.h"
#include "malloc.h"
#include "tensor.h"
#include "node.h"
#include "tensorblock.h"
#include "taskqueue.h"

#include <libsmartptr/strop.h>
#include "data.h"
#include "sort.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

using namespace yeti;
using namespace std;

static std::map<void*,size_t>* ptrs = 0;

#define RUNTIME_BLOCK_MALLOC 0

#if RUNTIME_BLOCK_MALLOC
#define NBLOCKS_RUNTIME_DATA 20000
typedef char runtime_data_block_t[200000];
DECLARE_MALLOC(runtime_data_block_t);
#endif

YetiOStream YetiRuntime::yetiout;
std::map<std::string, IndexDescrPtr> YetiRuntime::descrs_;
std::map<IndexRange*, std::map<IndexRange*, int> > YetiRuntime::valid_subranges_;
uli* YetiRuntime::max_range_sizes_ = 0;
uli* YetiRuntime::max_block_sizes_ = 0;
usi YetiRuntime::max_depth_ = 0;
usi YetiRuntime::max_nindex_ = 0;
uli YetiRuntime::nthread_ = 0;
uli YetiRuntime::nproc_ = 0;
uli YetiRuntime::me_ = 0;
size_t YetiRuntime::memory_ = 0;
size_t YetiRuntime::amount_mem_allocated_ = 0;
usi YetiRuntime::print_depth_ = 0;
float YetiRuntime::nonnull_cutoff = -20;
float YetiRuntime::matrix_multiply_cutoff = -20;
float YetiRuntime::print_cutoff = -6;
ThreadGroup* YetiRuntime::threadgrp_ = 0;
ThreadLock* YetiRuntime::printlock_ = 0;
ThreadLock* YetiRuntime::malloc_lock_ = 0;
bool YetiRuntime::threaded_runtime_ = false;
char* YetiRuntime::thread_stack_ = 0;
size_t YetiRuntime::thread_stack_size_ = 0;
bool YetiRuntime::print_cxn = false;
uli YetiRuntime::descr_id_count_ = DEFAULT_INDEX_DESCR_ID;

static FastMalloc matrix_generator_set_malloc("matrix generator set");
static FastMalloc matrix_generator_malloc("matrix generator");
static FastMalloc indexset_malloc("index set");
static FastMalloc tile_location_malloc("tile location");
static MallocOverride use_heap;


IndexDescr*
YetiRuntime::get_descr(const std::string& id)
{
    map<std::string, IndexDescrPtr>::iterator it
        = descrs_.find(id);

    if (it == descrs_.end())
    {
        cerr << "Cannot find index descr " << id << endl;
        abort();
    }

    return it->second.get();
}

IndexRange*
YetiRuntime::get_range(const std::string& id)
{
    map<std::string, IndexDescrPtr>::iterator it
        = descrs_.find(id);

    if (it == descrs_.end())
    {
        cerr << "Cannot find index range " << id << endl;
        abort();
    }

    return it->second->get_top_range();
}

void
YetiRuntime::init(
    uli me,
    uli nproc,
    usi nindex,
    uli nthread,
    size_t memory,
    size_t cache
)
{
    ptrs = new std::map<void*,size_t>;

    setenv("OMP_NUM_THREADS","1",1);

    //the maximum number of indices to be allocated
    max_nindex_ = nindex;
    nproc_ = nproc;
    me_ = me;
    memory_ = memory;
    amount_mem_allocated_ = 0;
    nthread_ = nthread;

    //validate that we have index ranges to work with
    if (descrs_.size() == 0)
    {
        cerr << "There are no index descrs during YetiRuntime::init()" << endl;
        abort();
    }

    {
        //validate that all index ranges have the same depth!
        //this is a new requirement for simplicity
        map<std::string, IndexDescrPtr>::const_iterator it(descrs_.begin());
        map<std::string, IndexDescrPtr>::const_iterator stop(descrs_.end());
        IndexDescrPtr topdescr = it->second;
        usi alignment_depth = topdescr->depth();
        ++it; //advance to next
        for ( ; it != stop; ++it)
        {
            IndexDescrPtr descr = it->second;
            usi depth = descr->depth();
            if (depth != alignment_depth)
            {
                cerr << "Not all index descriptors in Yeti aligned at same depth." << endl;
                cerr << "These two index descriptors differ: " << endl;
                topdescr->print(cerr); cerr << endl;
                descr->print(cerr); cerr << endl;
                abort();
            }
        }

        uli start = 0;
        uli nindex = 1;
        IndexRange* subrange = new IndexRange(start, nindex);
        IndexRangePtr dotrange = new IndexRange(subrange, alignment_depth, start);
        IndexDescrPtr dotprod_descr = new IndexDescr("dot product range", dotrange);
        YetiRuntime::register_index_descr("dotproduct", dotprod_descr);
    }

    init_sizes();
    init_malloc();
    init_threads(nthread);
    GlobalQueue::init();
}

void
YetiRuntime::init_sizes()
{
    if (descrs_.size() == 0)
    {
        cerr << "YetiRuntime has no index ranges" << endl;
        abort();
    }

    //figure out the maximum depth
    max_depth_ = 0;
    std::map<std::string, IndexDescrPtr>::iterator start(descrs_.begin());
    std::map<std::string, IndexDescrPtr>::iterator stop(descrs_.end());
    std::map<std::string, IndexDescrPtr>::iterator it;
    for (it=start; it != stop; ++it)
    {
        const IndexDescrPtr& range(it->second);
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
            const IndexDescrPtr& range = it->second;
            uli size = range->get_top_range()->nmax(depth);
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

    StorageBlock::init_malloc(NMALLOC_BLOCKS_DEBUG);
    DataStorageNode::init_malloc(NMALLOC_BLOCKS_DEBUG);
    TensorBlock::init_malloc(NMALLOC_BLOCKS_DEBUG / 5);
    TensorController::init_malloc(NMALLOC_BLOCKS_DEBUG * 4 / 5);
    MemoryPool::init_malloc(NMALLOC_BLOCKS_DEBUG / 5);

#if RUNTIME_BLOCK_MALLOC
    FastMallocTemplate<runtime_data_block_t>::init(NBLOCKS_RUNTIME_DATA);
#endif

    Tensor::init_malloc(1000);
    DataCache::init_malloc(1000);

    #define NCXN_TASKS_MALLOC 10000000
    Task::init_malloc(NCXN_TASKS_MALLOC);

    #define NSORT_OBJECTS 12
    Sort::init_malloc(NSORT_OBJECTS);

    #define NTHREAD_BLOCKS_MALLOC 10000000
    ThreadLock::init_malloc(NTHREAD_BLOCKS_MALLOC);

    #define NCACHE_BLOCKS_MALLOC 1000000
    DataCacheEntry::init_malloc(NCACHE_BLOCKS_MALLOC);

    #define PERMSET_SIZE 40
    #define NPERMSETS 100000
    matrix_generator_set_malloc.init(PERMSET_SIZE, NPERMSETS);

    #define NMATRIX_GENERATORS 10000000
    matrix_generator_malloc.init(sizeof(void*), NMATRIX_GENERATORS);

    #define NTILE_LOCATIONS 1000
    tile_location_malloc.init(sizeof(uli) * MAX_DEPTH * max_nindex_, NTILE_LOCATIONS);

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

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    int status = pthread_attr_getstacksize(&attr, &thread_stack_size_);
    if (status != 0)
    {
        cerr << "unable to get thread stack size" << endl;
        abort();
    }
    pthread_attr_destroy(&attr);

    thread_stack_ = YetiRuntime::malloc(thread_stack_size_ * nthread);

    threadgrp_ = new pThreadGroup(nthread);
    nthread_ = nthread;

    ThreadEnvironment::allocate();
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

const uli*
YetiRuntime::max_block_sizes()
{
    return max_block_sizes_;
}

const uli*
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
YetiRuntime::register_index_descr(
    const std::string& id,
    const IndexDescrPtr& descr
)
{
    descrs_[id] = descr;
}


void
YetiRuntime::register_index_range(
    IndexRange *range,
    const std::string& descr,
    const char *id1,
    const char *id2,
    const char *id3,
    const char *id4,
    const char *id5,
    const char *id6,
    const char *id7,
    const char *id8,
    const char *id9,
    const char *id10,
    const char *id11,
    const char *id12
)
{

    IndexDescrPtr indexdescr = new IndexDescr(descr, range);
    register_index_descr(
        indexdescr,
        id1, id2, id3, id4, id5, id6,
        id7, id8, id9, id10, id11, id12
    );
}

void
YetiRuntime::register_index_descr(
    const IndexDescrPtr& range,
    const char *id1,
    const char *id2,
    const char *id3,
    const char *id4,
    const char *id5,
    const char *id6,
    const char *id7,
    const char *id8,
    const char *id9,
    const char *id10,
    const char *id11,
    const char *id12
)
{
    register_index_descr(id1, range);

    if (id2) register_index_descr(id2, range);
    if (id3) register_index_descr(id3, range);
    if (id4) register_index_descr(id4, range);
    if (id5) register_index_descr(id5, range);
    if (id6) register_index_descr(id6, range);
    if (id7) register_index_descr(id7, range);
    if (id8) register_index_descr(id8, range);
    if (id9) register_index_descr(id9, range);
    if (id10) register_index_descr(id10, range);
    if (id11) register_index_descr(id11, range);
    if (id12) register_index_descr(id12, range);
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
yeti::yeti_malloc_indexset(uli threadnum)
{
    return reinterpret_cast<uli*>(indexset_malloc.malloc(threadnum));
}

void
yeti::yeti_free_indexset(const uli* ptr)
{
    indexset_malloc.free(ptr);
}

void*
yeti::yeti_malloc_indexptr(uli threadnum)
{
    return indexset_malloc.malloc(threadnum);
}

void
yeti::yeti_free_indexptr(void *ptr)
{
    indexset_malloc.free(ptr);
}

usi*
yeti::yeti_malloc_perm(uli threadnum)
{
    return reinterpret_cast<usi*>(indexset_malloc.malloc(threadnum));
}

void
yeti::yeti_free_perm(usi* ptr)
{
    indexset_malloc.free(ptr);
}

void*
yeti::yeti_malloc_matrix_generator(uli threadnum)
{
    return matrix_generator_malloc.malloc(threadnum);
}

void
yeti::yeti_free_matrix_generator(void *ptr)
{
    matrix_generator_malloc.free(ptr);
}

void*
yeti::yeti_malloc_matrix_generator_set(uli threadnum)
{
    return matrix_generator_set_malloc.malloc(threadnum);
}

void
yeti::yeti_free_matrix_generator_set(void* ptr)
{
    matrix_generator_set_malloc.free(ptr);
}

uli*
yeti::yeti_malloc_tile_location(uli threadnum)
{
    return reinterpret_cast<uli*>(tile_location_malloc.malloc(threadnum));
}

void
yeti::yeti_free_tile_location(void* ptr)
{
    tile_location_malloc.free(ptr);
}

size_t
YetiRuntime::get_thread_stack_size()
{
    return thread_stack_size_;
}

uli
YetiRuntime::get_thread_number()
{
    if (nthread_ == 1)
        return 0;

    int x; void* ptr = &x;

    if (ptr < thread_stack_)
        return 0;

    size_t distance = (size_t) ptr - (size_t) thread_stack_;
    uli threadnum = distance / thread_stack_size_;

    if (threadnum >= nthread_)
        return 0;

    return threadnum;
}

void*
YetiRuntime::get_thread_stack(uli threadnum)
{
    return thread_stack_ + threadnum * thread_stack_size_;
}


void
YetiRuntime::finalize()
{
    GlobalQueue::finalize();

    if (thread_stack_)
    {
        YetiRuntime::free(thread_stack_, nthread_ * thread_stack_size_);
    }

    ThreadEnvironment::deallocate();
    if (printlock_) //these came off standard malloc heap
    {
        ::free(printlock_);
        printlock_ = 0;
    }
    if (malloc_lock_)
    {
        ::free(malloc_lock_);
        malloc_lock_ = 0;
    }
    if (threadgrp_) //this came from operator new
        delete threadgrp_;
}

usi
YetiRuntime::print_depth()
{
    return print_depth_;
}

void
YetiRuntime::set_print_depth(usi depth)
{
    print_depth_ = depth;
}

IndexRange*
yeti::index_range(const std::string& str)
{
    return YetiRuntime::get_range(str);
}

void
YetiRuntime::register_subranges(
    IndexRange* parent, 
    IndexRange* subrange1, 
    IndexRange* subrange2, 
    IndexRange* subrange3, 
    IndexRange* subrange4, 
    IndexRange* subrange5, 
    IndexRange* subrange6
)
{
    std::map<IndexRange*, int>& range_map = valid_subranges_[parent];
    range_map[subrange1] = 1;
    
    if (subrange2) range_map[subrange2] = 1;
    if (subrange3) range_map[subrange3] = 1;
    if (subrange4) range_map[subrange4] = 1;
    if (subrange5) range_map[subrange5] = 1;
    if (subrange6) range_map[subrange6] = 1;
}

void
YetiRuntime::register_subranges(
    const char* parent,
    const char* subrange1,
    const char* subrange2,
    const char* subrange3,
    const char* subrange4,
    const char* subrange5,
    const char* subrange6
)
{

    IndexRange *_parent = 0,
            *_subrange1 = 0,
            *_subrange2 = 0,
            *_subrange3 = 0,
            *_subrange4 = 0,
            *_subrange5 = 0,
            *_subrange6 = 0;

    _parent = YetiRuntime::get_range(parent);
    _subrange1 = YetiRuntime::get_range(subrange1);

    if (subrange2) _subrange2 = get_range(subrange2);
    if (subrange3) _subrange3 = get_range(subrange3);
    if (subrange4) _subrange4 = get_range(subrange4);
    if (subrange5) _subrange5 = get_range(subrange5);
    if (subrange6) _subrange6 = get_range(subrange6);

    register_subranges(
        _parent,
        _subrange1,
        _subrange2,
        _subrange3,
        _subrange4,
        _subrange5,
        _subrange6
    );
}

void
YetiRuntime::register_equivalent_descrs(
    const char* s1,
    const char* s2,
    const char* s3,
    const char* s4,
    const char* s5,
    const char* s6,
    const char* s7,
    const char* s8
)
{
    ++descr_id_count_;

    YetiRuntime::get_descr(s1)->set_descr_id(descr_id_count_);
    YetiRuntime::get_descr(s2)->set_descr_id(descr_id_count_);
    if (s3) YetiRuntime::get_descr(s3)->set_descr_id(descr_id_count_);
    if (s4) YetiRuntime::get_descr(s4)->set_descr_id(descr_id_count_);
    if (s5) YetiRuntime::get_descr(s5)->set_descr_id(descr_id_count_);
    if (s6) YetiRuntime::get_descr(s6)->set_descr_id(descr_id_count_);
    if (s7) YetiRuntime::get_descr(s7)->set_descr_id(descr_id_count_);
    if (s8) YetiRuntime::get_descr(s8)->set_descr_id(descr_id_count_);
}

bool
YetiRuntime::is_valid_subrange(
    IndexRange* parent, 
    IndexRange* subrange
)
{
   std::map<IndexRange*, std::map<IndexRange*, int> >::const_iterator it
        = valid_subranges_.find(parent);

    if (it == valid_subranges_.end())
        return false;

    const std::map<IndexRange*, int>& range_map = it->second;
    std::map<IndexRange*, int>::const_iterator it2 = range_map.find(subrange);
    return it2 != range_map.end();
}


void get_memory(size_t mem, size_t& ngb, size_t& nmb)
{
    ngb = mem / 1e9;
    size_t remainder = mem - ngb * 1e9;
    nmb = remainder / 1e6;
}

char*
YetiRuntime::malloc(size_t size)
{
#if RUNTIME_BLOCK_MALLOC
    if (size > sizeof(runtime_data_block_t))
    {
        cerr << "size error!" << endl;
        abort();
    }
    void* ptr = FastMallocTemplate<runtime_data_block_t>::malloc->malloc();
    return reinterpret_cast<char*>(ptr);
#else
    YetiRuntime::lock_malloc();
    amount_mem_allocated_ += size;
    if (amount_mem_allocated_ > memory_)
    {
        cerr << "Allocation request exceeds total memory" << endl;
        size_t ngb, nmb;
        get_memory(memory_, ngb, nmb);
        cerr << stream_printf("Total memory: %ld GB %ld MB", ngb, nmb) << endl;
        get_memory(amount_mem_allocated_, ngb, nmb);
        cerr << stream_printf("Total allocated: %ld GB %ld MB", ngb, nmb) << endl;
        abort();
    }
    char* data = (char*) ::malloc(size);
    if (data == 0)
    {
        cerr << "Allocated failed for " << size << " bytes" << endl;
        size_t ngb, nmb;
        get_memory(memory_, ngb, nmb);
        cerr << stream_printf("Total memory: %ld GB %ld MB", ngb, nmb) << endl;
        abort();
    }
    (*ptrs)[data] = size;
    YetiRuntime::unlock_malloc();
    return data;
#endif
}

void
YetiRuntime::print_memory_allocation()
{
    size_t ngb, nmb;
    get_memory(amount_mem_allocated_, ngb, nmb);
    cout << stream_printf("Total allocated: %ld GB %ld MB", ngb, nmb) << endl;
}

void
YetiRuntime::free(void *ptr, size_t size)
{
 #if RUNTIME_BLOCK_MALLOC
    FastMallocTemplate<runtime_data_block_t>::malloc->free(ptr);
    return;
#else
    if(ptrs){
        std::map<void*,size_t>::iterator it = ptrs->find(ptr);
        if (it == ptrs->end())
        {
            cerr << "freeing non-malloc'd pointer " << ptr << endl;
            abort();
        }
        if (it->second != size)
        {
            cerr << "freed sizes don't match" << endl;
            abort();
        }
        YetiRuntime::lock_malloc();
        ::free(ptr);
        ptrs->erase(it);
        amount_mem_allocated_ -= size;
        YetiRuntime::unlock_malloc();
    }
#endif
}

YetiOStream&
yeti::operator<<(YetiOStream& os, const std::string& str)
{
    os.print(str);
    return os;
}

YetiOStream&
yeti::operator<<(YetiOStream& os, std::ostream& (*pf)(std::ostream&))
{
    os.print(pf);
    return os;
}

YetiOStream&
yeti::operator<<(YetiOStream& os, void* ptr)
{
    os.print(ptr);
    return os;
}

