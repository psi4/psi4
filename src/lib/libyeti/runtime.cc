#include "runtime.h"
#include "thread.h"
#include "index.h"
#include "malloc.h"
#include "tensor.h"
#include "node.h"
#include "tensorblock.h"
#include "taskqueue.h"
#include "aobasis.h"
#include "env.h"
#include "messenger.h"
#include "exception.h"

#include <libsmartptr/strop.h>
#include "data.h"
#include "sort.h"

#include <psiconfig.h>
#if HAVE_MPI
#include <mpi.h>
#endif

#if HAVE_MPQC
#include <util/group/message.h>
#include <chemistry/qc/zaptr12/matrixkits.h>
#endif

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

using namespace yeti;
using namespace std;

#define NTHREAD_COMM 2


#if USE_DEFAULT_THREAD_STACK
static void* thread_stacks[10];
#endif

YetiOStream YetiRuntime::yetiout;
std::map<std::string, IndexDescrPtr> YetiRuntime::descrs_;
std::map<IndexRange*, std::map<IndexRange*, int> > YetiRuntime::valid_subranges_;
uli* YetiRuntime::max_range_sizes_ = 0;
uli* YetiRuntime::max_block_sizes_ = 0;
usi YetiRuntime::max_depth_ = 0;
uli YetiRuntime::nthread_compute_ = 0;
uli YetiRuntime::nthread_comm_ = 0;
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
bool YetiRuntime::threaded_compute_ = false;
char* YetiRuntime::thread_stack_ = 0;
size_t YetiRuntime::thread_stack_size_ = 0;
bool YetiRuntime::print_cxn = false;
uli YetiRuntime::descr_id_count_ = DEFAULT_INDEX_DESCR_ID;
uli YetiRuntime::nblocks_to_allocate_ = 0;
std::map<std::string, AOBasisPtr> YetiRuntime::basis_sets_;
YetiRuntime::runtime_parallel_type_t YetiRuntime::parallel_type_;
YetiMPIMessenger* YetiRuntime::mpi_messenger_ = 0;
YetiMessenger* YetiRuntime::messenger_ = 0;

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
YetiRuntime::init_system_environment()
{
    const char* env = getenv("MPICH_MAX_SHORT_MSG_SIZE");
    if (env)
        cout << env << endl;
    setenv("MPICH_MAX_SHORT_MSG_SIZE","0",1);
    setenv("OMP_NUM_THREAD","1",1);
    const char* nthread_str = getenv("YETI_NUM_THREAD");
    if (nthread_str)
    {
        nthread_compute_ = atol(nthread_str);
        if (nthread_compute_ == 0)
        {
            cerr << "invalid thread specification " << nthread_str << endl;
            abort();
        }
    }
    else
    {
        nthread_compute_ = 1;
    }


    const char* memory_str = getenv("YETI_MEM");
    memory_ = 0;
    if (memory_str)
    {
        double memdbl = atof(memory_str);
        memory_ = memdbl;
        if (memory_ == 0) //formatted as an integer
        {
            memory_ = atol(memory_str);
        }
        if (memory_ == 0)
        {
            cerr << "invalid memory specification " << memory_str << endl;
            abort();
        }
    }
    else
    {
        memory_ = 1e8;
    }

    amount_mem_allocated_ = 0;

#include <psiconfig.h>
#if HAVE_MPI
    parallel_type_ = YetiRuntime::MPI;
    init_mpi();
#else
    parallel_type_ = YetiRuntime::Local;
#endif

    nthread_ = nthread_compute_ + nthread_comm_;

    init_malloc();

    init_threads();

    if (parallel_type_ == MPI)
        init_mpi_messenger();
    else
        init_local_messenger();

    GlobalQueue::init();

    TensorBlock::init_statics();
}

void
YetiRuntime::init_local_messenger()
{
    me_ = 0;
    nproc_ = 1;
    messenger_ = new YetiLocalMessenger;
    nthread_comm_ = 0;
}

void
YetiRuntime::init_molecule_environment()
{
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
        for (usi idx=0; idx < NINDEX; ++idx)
        {
            max_block_size *= max_range_size;
        }
        max_block_sizes_[depth] = max_block_size;
    }
}

void
YetiRuntime::init_malloc()
{
    /** Give up 1% of memory to allocate tensor blocks */
    uli block_mem_allocation = memory_ / 100;
    nblocks_to_allocate_ = block_mem_allocation / sizeof(TensorBlock);

    StorageBlock::init_malloc(nblocks_to_allocate_ * 5);
    DataStorageNode::init_malloc(nblocks_to_allocate_ * 5);
    TensorBlock::init_malloc(nblocks_to_allocate_);
    TensorController::init_malloc(nblocks_to_allocate_ * 4);
    MemoryPool::init_malloc(nblocks_to_allocate_);


    uli ntensors = 1000;
    Tensor::init_malloc(ntensors);
    DataCache::init_malloc(ntensors);

    uli ntasks = 1e5;
    Task::init_malloc(ntasks);

    uli nsort_objects = (nthread_compute_ + 1) * 3;
    Sort::init_malloc(nsort_objects);

    uli nthread_locks =  nblocks_to_allocate_ * 2;
    ThreadLock::init_malloc(nthread_locks);

    uli ncache_blocks = 1e6;
    DataCacheEntry::init_malloc(ncache_blocks);

    //make sure all mallocs have been configured
#if USE_MALLOC_CLASS
    FastMalloc::check_queue();
#endif
}

void
YetiRuntime::init_threads()
{
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    int status = pthread_attr_getstacksize(&attr, &thread_stack_size_);
    if (status != 0)
    {
        cerr << "unable to get thread stack size" << endl;
        abort();
    }

    if (thread_stack_size_ == 0)
    {
        cerr << "Default thread stack size is zero.  Go find a posix programmer and strangle him." << endl;
        abort();
    }

    pthread_attr_destroy(&attr);

#if USE_DEFAULT_THREAD_STACK
    thread_stack_ = 0;
    for (uli t=0; t < nthread_; ++t)
        thread_stacks[t] = 0;
#else
    thread_stack_ = YetiRuntime::malloc(thread_stack_size_ * nthread_);
#endif
    threadgrp_ = new pThreadGroup(nthread_compute_);
    printlock_ = new (use_heap) pThreadLock;
    malloc_lock_ = new (use_heap) pThreadLock;


    ThreadEnvironment::allocate();
}

usi
YetiRuntime::max_depth()
{
    return max_depth_;
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

size_t
YetiRuntime::memory()
{
    return memory_;
}

uli
YetiRuntime::nthread_compute()
{
    return nthread_compute_;
}

uli
YetiRuntime::nthread_comm()
{
    return nthread_comm_;
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
YetiRuntime::register_index_range(
    const IndexRangePtr& range,
    const std::string& name,
    const std::string& id

)
{
    IndexDescrPtr descr = new IndexDescr(name, range);
    register_index_descr(id, descr);
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

#if USE_DEFAULT_THREAD_STACK
    int x; 
    void* ptr = &x;
    size_t offset = (size_t) ptr;
    for (uli t=0; t < nthread_; ++t)
    {
        size_t stack = (size_t) thread_stacks[t];
        if (stack && offset < stack && (stack - offset) <= thread_stack_size_)
        {
            return t;
        }
    }
    return 0;
#else
    int x; void* ptr = &x;

    if (ptr < thread_stack_)
        return 0;

    size_t distance = (size_t) ptr - (size_t) thread_stack_;
    uli threadnum = distance / thread_stack_size_;

    if (threadnum >= nthread_ + NTHREAD_COMM)
        return 0;

    return threadnum;
#endif
}

void*
YetiRuntime::get_thread_stack(uli threadnum)
{
    return thread_stack_ + threadnum * thread_stack_size_;
}

#if USE_DEFAULT_THREAD_STACK
void
YetiRuntime::set_thread_stack(uli threadnum, void* stack)
{
    thread_stacks[threadnum] = stack;
}
#endif


void
YetiRuntime::finalize()
{
    TensorBlock::delete_statics();

    if (parallel_type_ == MPI)
    {
#include <psiconfig.h>
#if HAVE_MPI
        mpi_messenger_->wait_barrier();
        mpi_messenger_->stop();
        sc::MatrixKits::defaultkit = 0;
        sc::MatrixKits::blockkit = 0;
        sc::MessageGrp::set_default_messagegrp(0);

        int finalized;
        int err = MPI_Finalized(&finalized);
        if (err != MPI_SUCCESS)
        {
            cerr << "MPI Finalized call failed on node " << YetiRuntime::me() << endl;
            abort();
        }
        if (!finalized)
            MPI_Finalize();

        delete mpi_messenger_;
        mpi_messenger_ = 0;
#endif
    }
    else
    {
        delete messenger_;
        messenger_ = 0;
    }

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
    YetiRuntime::lock_malloc();
    ::free(ptr);
    amount_mem_allocated_ -= size;
    YetiRuntime::unlock_malloc();
#endif
}

bool
YetiRuntime::is_threaded_compute()
{
    return threaded_compute_;
}

void
YetiRuntime::set_threaded_compute(bool flag)
{
    threaded_compute_ = flag;
}

void
YetiRuntime::register_basis_set(const AOBasisPtr& basis)
{
    basis_sets_[basis->name()] = basis;
}

const AOBasisPtr&
YetiRuntime::get_basis(const std::string& id)
{
    return basis_sets_[id];
}

YetiOStream&
yeti::operator<<(YetiOStream& os, const std::string& str)
{
    os.print<std::string>(str);
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
    os.print<void*>(ptr);
    return os;
}

YetiMessenger*
YetiRuntime::get_messenger()
{
    return messenger_;
}

void
YetiRuntime::init_mpi()
{
#include <psiconfig.h>
#if HAVE_MPI
    int argc = 0;
    char** argv = new char*[argc + 1];
    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS)
    {
        Env::err0() << "MPI failed to init properly with error " << err << endl;
        abort();
    }

    int me, nproc;

    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    me_ = me;
    nproc_ = nproc;

    if (nproc_ > 1)
        nthread_comm_ = 2;
    else
        nthread_comm_ = 0;

    delete[] argv;
#endif
}

void
YetiRuntime::init_mpi_messenger()
{
#include <psiconfig.h>
#if HAVE_MPI
    if (nproc_ == 1) //mpi started with only 1 process
    {
        parallel_type_ = YetiRuntime::Local;
        MPI_Finalize();
        init_local_messenger();
    }
    else
    {
        mpi_messenger_ = new YetiMPIMessenger;
        messenger_ = mpi_messenger_;
        mpi_messenger_->start();
    }
#else
    raise(SanityCheckError, "calling init mpi but YETI not compiled for mpi");
#endif
}
