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

#undef raise
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>

#include <libsmartptr/strop.h>
#include "data.h"
#include "sort.h"

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


#if USE_DEFAULT_THREAD_STACK
static void* thread_stacks[10];
#endif

YetiOStream YetiRuntime::yetiout;
std::map<std::string, IndexDescrPtr> YetiRuntime::descrs_;
std::map<IndexRange*, std::map<IndexRange*, int> > YetiRuntime::valid_subranges_;
std::map<IndexDescr*, std::map<usi, std::string> > YetiRuntime::descr_letters_;
std::map<Tensor*, size_t> YetiRuntime::tensor_allocations_;
std::map<DataCache*, size_t> YetiRuntime::cache_allocations_;
uli* YetiRuntime::max_range_sizes_ = 0;
uli* YetiRuntime::max_block_sizes_ = 0;
uli YetiRuntime::cpu_mask_num_ = 0;
uli YetiRuntime::num_nodes_task_group_ = 0;
std::map<uli,uli> YetiRuntime::group_to_global_node_map_;
std::map<uli,uli> YetiRuntime::global_to_group_node_map_;
usi YetiRuntime::max_depth_ = 0;
uli YetiRuntime::nthread_ = 0;
uli YetiRuntime::nproc_ = 0;
uli YetiRuntime::me_ = 0;
size_t YetiRuntime::memory_ = 0;
size_t YetiRuntime::amount_mem_allocated_ = 0;
usi YetiRuntime::print_depth_ = 0;
float YetiRuntime::nonnull_cutoff = -8;
float YetiRuntime::matrix_multiply_cutoff = -12;
float YetiRuntime::print_cutoff = -6;
ThreadGroup* YetiRuntime::threadgrp_ = 0;
ThreadLock* YetiRuntime::printlock_ = 0;
ThreadLock* YetiRuntime::malloc_lock_ = 0;
ThreadLock* YetiRuntime::timer_lock_ = 0;
bool YetiRuntime::threaded_compute_ = false;
char* YetiRuntime::thread_stack_ = 0;
size_t YetiRuntime::thread_stack_size_ = 0;
bool YetiRuntime::print_cxn = false;
bool YetiRuntime::dynamic_load_balance_ = false;
uli YetiRuntime::descr_id_count_ = DEFAULT_INDEX_DESCR_ID;
uli YetiRuntime::nblocks_to_allocate_ = 0;
std::map<std::string, AOBasisPtr> YetiRuntime::basis_sets_;
YetiRuntime::runtime_parallel_type_t YetiRuntime::parallel_type_;
YetiMPIMessenger* YetiRuntime::mpi_messenger_ = 0;
YetiMessenger* YetiRuntime::messenger_ = 0;
#if HAVE_MPI
MPI_Comm YetiRuntime::mpicomm_;
#endif
#if COUNT_SCREENING_SKIPS
uli* YetiRuntime::screening_skips_ = 0;
#endif
usi YetiRuntime::nindex_min_distr_ = 0;
bool YetiRuntime::use_auto_distr_ = false;

/** Program exit variables */
pThreadLock exit_lock;
stringstream exit_sstr;

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
YetiRuntime::stack_print()
{
    void* array[20];
    int size = backtrace(array, 20); 
    char** symbols = backtrace_symbols(array, size);
    stringstream sstr;
    for (int i=0; i < size; ++i)
    {
        sstr << symbols[i] << "\n";
    }
    cerr << sstr.str();
}

void
YetiRuntime::exit()
{
    exit_lock.lock();
    cerr << stream_printf("YETI exiting main thread on node %ld\n", me_);
    stack_print();
    exit_lock.unlock();
    abort();
}

void
YetiRuntime::exit_with_signal(int sig)
{
    if (me_ != 0)
        abort();

    if (YetiRuntime::get_thread_number() != 0)
    {
        cerr << "YETI not exiting on main thread!" << endl;
        YetiRuntime::exit();
    }


    exit_lock.lock();
    if (sig == SIGSEGV)
    {
        cerr << stream_printf("YETI Segmentation Fault on Node %ld", YetiRuntime::me()) << endl;
        YetiRuntime::stack_print();
    }
    else if (sig == SIGINT)
    {
        cerr << stream_printf("YETI Interrupt on Node %ld", YetiRuntime::me()) << endl;
    }
    else if (sig == SIGKILL || sig == SIGTERM)
    {
        cerr << stream_printf("YETI Terminate on Node %ld", YetiRuntime::me()) << endl;
    }
    else if (sig == SIGFPE)
    {
        cerr << stream_printf("YETI Floating Point Exception on Node %ld", YetiRuntime::me()) << endl;
    }

    get_messenger()->send_exit_signal();

    get_thread_grp()->thread_crash();

    get_messenger()->thread_crash();

    exit_lock.unlock();

    sleep(1);

    YetiRuntime::exit();
}

void
YetiRuntime::thread_exit(int sig)
{
    if (sig != THREAD_KILL_SIGNAL)
    {
        cerr << "Thread received invalid signal " << sig << endl;
        YetiRuntime::exit();
    }
    exit_lock.lock();
    cerr << stream_printf("\nExit on thread %ld on node %ld\n", YetiRuntime::get_thread_number(), me_);
    YetiRuntime::stack_print();
    cerr << "\n\n";
    cerr.flush();
    exit_lock.unlock();
}

void
YetiRuntime::main_process_backtrace(int sig)
{
    if (sig != MAIN_PROCESS_BACKTRACE_SIGNAL)
    {
        cerr << "Main process backtrace received invalid signal " << sig << endl;
    }
    exit_lock.lock();
    cerr << stream_printf("Backtrace from node %ld\n", me_);
    YetiRuntime::stack_print();
    exit_lock.unlock();
    while(1);
}

void
YetiRuntime::init_system_environment()
{
    const char* nthread_str = getenv("YETI_NUM_THREAD");
    if (nthread_str)
    {
        nthread_ = atol(nthread_str);
        if (nthread_ == 0)
        {
            cerr << "invalid thread specification YETI_NUM_THREAD=" << nthread_str << endl;
            abort();
        }
    }
    else
    {
        nthread_ = 1;
    }

    signal(SIGSEGV, YetiRuntime::exit_with_signal);
    //signal(SIGINT, YetiRuntime::exit_with_signal);
    signal(SIGTERM, YetiRuntime::exit_with_signal);
    signal(SIGKILL, YetiRuntime::exit_with_signal);
    signal(MAIN_PROCESS_BACKTRACE_SIGNAL, YetiRuntime::main_process_backtrace);
    signal(THREAD_KILL_SIGNAL, YetiRuntime::thread_exit);
    signal(SIGFPE, YetiRuntime::exit_with_signal);

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
        memory_ = 1e9;
    }

    amount_mem_allocated_ = 0;

#if HAVE_MPI
    parallel_type_ = YetiRuntime::MPI;
    init_mpi();
#else
    parallel_type_ = YetiRuntime::Local;
#endif

    init_malloc();

    init_threads();

    if (parallel_type_ == MPI)
        init_mpi_messenger();
    else
        init_local_messenger();


    #define NUM_LINKED_NODES 20
    num_nodes_task_group_ = NUM_LINKED_NODES;
    uli num_nodes_max = YetiRuntime::nproc() - 1;
    if (num_nodes_task_group_ > num_nodes_max)
        num_nodes_task_group_ = num_nodes_max;

    for (uli i=0, p=YetiRuntime::me() + 1; i < num_nodes_task_group_; ++i, ++p)
    {
        uli node  = p % YetiRuntime::nproc();
        group_to_global_node_map_[i] = node;
        global_to_group_node_map_[node] = i;
        //cout << stream_printf("Node %d has group member %d -> %d\n", YetiRuntime::me(), i, node);
        //cout.flush();
    }

    const char* load_str = getenv("F_DYNAMIC_LOAD_BALANCE");
    if (load_str)
    {
        dynamic_load_balance_ = atoi(load_str);
    }

    GlobalQueue::init();

    TensorBlock::init_statics();

	uli mem_mb = amount_mem_allocated_ / 1e6;
	if (me_ == 0)
		cout << stream_printf("Core system malloc used %ld MB of memory", mem_mb) << endl;

    if (me_ == 0)
    {
        cout << stream_printf("Size of single tensor block is %ld bytes\n",
                                     sizeof(TensorBlock));
    }
}

void
YetiRuntime::init_local_messenger()
{
    me_ = 0;
    nproc_ = 1;
    messenger_ = new YetiLocalMessenger;
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

#if COUNT_SCREENING_SKIPS
    // Now that we have the max_depth_, initialize screening_skips_
    screening_skips_ = new uli[max_depth_ + 1];
    for(int i = 0; i <= max_depth_; ++i) screening_skips_[i] = 0;
#endif


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
    // TODO (Low priority) Make this expandable, which could be done with MallocSegment subsections of the Malloc'd symmetric heap which keep track of offsets on a given node making pointer serialization just as easy as it is now.
    /** Give up 1% of memory to allocate tensor blocks */
    uli block_mem_allocation = 2 * memory_ / 100;
    nblocks_to_allocate_ = block_mem_allocation / sizeof(TensorBlock);


    StorageBlock::init_malloc(nblocks_to_allocate_ * 5);
    DataStorageNode::init_malloc(nblocks_to_allocate_ * 5);
    TensorBlock::init_malloc(nblocks_to_allocate_);
    TensorController::init_malloc(50);
    MemoryPool::init_malloc(nblocks_to_allocate_);


    uli ncxns = 100;
    Contraction::init_malloc(ncxns);

    uli ntensors = 500;
    Tensor::init_malloc(ntensors);
    DataCache::init_malloc(ntensors);

    uli ntasks_cxn = 1e6;
    ContractionTask::init_malloc(ntasks_cxn);

    uli nsort_objects = nthread_ * 3;
    Sort::init_malloc(nsort_objects);

    uli nthread_locks =  nblocks_to_allocate_ * 2;
    ThreadLock::init_malloc(nthread_locks);

    uli ncache_blocks = ntensors * 20;
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
    uli nstacks = nthread_;
    if (YetiRuntime::nproc() > 1)
        ++nstacks; //another thread for recv
    uli nblocks = nstacks + 1; //stupid thing
    thread_stack_ = YetiRuntime::malloc(thread_stack_size_ * nblocks);
#endif
    threadgrp_ = new pThreadGroup(nthread_);
    printlock_ = new (use_heap) pThreadLock;
    malloc_lock_ = new (use_heap) pThreadLock;
    timer_lock_ = new (use_heap) pThreadLock;


    ThreadEnvironment::allocate();
}

void
YetiRuntime::start_timer(const char* str)
{
    timer_lock_->lock();
    timer::Timer::start(str);
    timer_lock_->unlock();
}

void
YetiRuntime::stop_timer(const char* str)
{
    timer_lock_->lock();
    timer::Timer::stop(str);
    timer_lock_->unlock();
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
YetiRuntime::num_nodes_task_group()
{
    return num_nodes_task_group_;
}

uli
YetiRuntime::group_node_number(uli global_node_number)
{
    return global_to_group_node_map_[global_node_number];
}

uli
YetiRuntime::global_node_number(uli group_node_number)
{
    return group_to_global_node_map_[group_node_number];
}

uli
YetiRuntime::nproc()
{
    return nproc_;
}

uli
YetiRuntime::cpu_mask_num()
{
    return cpu_mask_num_;
}

size_t
YetiRuntime::memory()
{
    return memory_;
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

    IndexDescr* descr = range.get();
    descr_letters_[descr][0] = id1;
    if (id2) descr_letters_[descr][1] = id2;
    if (id3) descr_letters_[descr][2] = id3;
    if (id4) descr_letters_[descr][3] = id4;
    if (id5) descr_letters_[descr][4] = id5;
    if (id6) descr_letters_[descr][5] = id6;
    if (id7) descr_letters_[descr][6] = id7;
    if (id8) descr_letters_[descr][7] = id8;
    if (id9) descr_letters_[descr][8] = id9;
    if (id10) descr_letters_[descr][9] = id10;
    if (id11) descr_letters_[descr][10] = id11;
    if (id12) descr_letters_[descr][11] = id12;
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
    if (nthread_ == 1 && nproc_ == 1)
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

    if (nproc_ > 1)
    {
        if (threadnum >= (nthread_ + 1))
            return 0;
    }
    else
    {
        if (threadnum >= nthread_)
            return 0;
    }

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
#if HAVE_MPI
        mpi_messenger_->wait_barrier();
        mpi_messenger_->stop();
#if HAVE_MPQC
        sc::MatrixKits::defaultkit = 0;
        sc::MatrixKits::blockkit = 0;
        sc::MessageGrp::set_default_messagegrp(0);
#endif

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

#if COUNT_SCREENING_SKIPS
    YetiRuntime::free(screening_skips_, max_depth_ * sizeof(uli));
#endif


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
#if TRACK_MALLOC_REQUEST
    size_t ngb, nmb;
    get_memory(size, ngb, nmb);
    if(nmb > 50) {
        cout << stream_printf("Large Malloc request: %ld GB %ld MB", ngb, nmb) << endl;
        get_memory(amount_mem_allocated_, ngb, nmb);
        cout << stream_printf("   Allocation total so far: %ld GB %ld MB\n-----------------------------", ngb, nmb) << endl;
    }
#endif
    amount_mem_allocated_ += size;
    if (amount_mem_allocated_ > memory_)
    {
        cerr << "Allocation request exceeds total memory" << endl;
        size_t ngb, nmb;
        get_memory(memory_, ngb, nmb);
        cerr << stream_printf("Total memory: %ld GB %ld MB", ngb, nmb) << endl;
        get_memory(amount_mem_allocated_, ngb, nmb);
        cerr << stream_printf("Total allocated: %ld GB %ld MB", ngb, nmb) << endl;

        {
            std::map<Tensor*, size_t>::const_iterator it = tensor_allocations_.begin();
            std::map<Tensor*, size_t>::const_iterator end = tensor_allocations_.end();
            for ( ; it != end; ++it)
            {
                size_t mem_mb = it->second / 1e6;
                cerr << stream_printf("Tensor %s allocated %ld MB\n", it->first->get_name().c_str(), mem_mb);
            }
        }
        {
            std::map<DataCache*, size_t>::const_iterator it = cache_allocations_.begin();
            std::map<DataCache*, size_t>::const_iterator end = cache_allocations_.end();
            for ( ; it != end; ++it)
            {
                size_t mem_mb = it->second / 1e6;
                cerr << stream_printf("Cache %s allocated %ld MB\n", it->first->name().c_str(), mem_mb);
            }
        }
        cerr.flush();

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
    Env::out0() << stream_printf("Total allocated: %ld GB %ld MB", ngb, nmb) << endl;
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
YetiRuntime::use_dynamic_load_balancing()
{
    return dynamic_load_balance_;
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

#if HAVE_MPI
MPI_Comm&
YetiRuntime::get_mpi_comm()
{
    return mpicomm_;
}
#endif

void
YetiRuntime::init_mpi()
{
#if HAVE_MPI
    int argc = 0;
    char** argv = new char*[argc + 1];

    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized)
    {
        int provided;
        int desired = MPI_THREAD_MULTIPLE;
        int err = MPI_Init_thread(&argc, &argv, desired, &provided);
        if (err != MPI_SUCCESS)
        {
            Env::err0() << "MPI failed to init properly with error " << err << endl;
            abort();
        }
    }

    MPI_Comm_dup(MPI_COMM_WORLD,&mpicomm_);

    int me, nproc;

    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    me_ = me;
    nproc_ = nproc;

    delete[] argv;
#endif

#if USE_THREAD_AFFINITY
    cpu_set_t cpuset;
    const char* ncores_str = getenv("YETI_NUM_CORES");
    uli ncores_per_socket = 0;
    if (ncores_str)
    {
        ncores_per_socket = atol(ncores_str);
        cpu_mask_num_ = ncores_per_socket * YetiRuntime::me();
        CPU_SET(cpu_mask_num_, &cpuset);
        sched_setaffinity(0, sizeof(cpu_set_t), &cpuset);
    }
    else
    {
        sched_getaffinity(0,sizeof(cpu_set_t),&cpuset);
        for (int i=0; i < 1000; ++i)
        {
            if (CPU_ISSET(i,&cpuset))
            {
                cpu_mask_num_ = i;
                break;
            }
        }
    }
    //cout << stream_printf("CPU %d has affinity on node %d\n", cpu_mask_num_, YetiRuntime::me());
    //cout.flush();
#endif
}

void
YetiRuntime::init_mpi_messenger()
{
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
    yeti_throw(SanityCheckError, "calling init mpi but YETI not compiled for mpi");
#endif
}

std::string
YetiRuntime::get_name(TensorIndexDescr* descr)
{
    if (descr->nindex() == 0)
        return "";

    std::string name = descr_letter(descr, 0);
    for (usi i=1; i < descr->nindex(); ++i)
        name += "," + descr_letter(descr, i);

    return name;
}

const std::string&
YetiRuntime::descr_letter(TensorIndexDescr* descr, usi index)
{
    IndexDescr* idx_descr = descr->get(index);
    return descr_letters_[idx_descr][index];
}

void
YetiRuntime::register_allocation(DataCache* cache, size_t size)
{
    cache_allocations_[cache] += size;
}

void
YetiRuntime::register_cache(DataCache* cache)
{
    cache_allocations_[cache] = 0;
}

void
YetiRuntime::unregister_cache(DataCache* cache)
{
    std::map<DataCache*,size_t>::iterator it = cache_allocations_.find(cache);
    cache_allocations_.erase(it);
}

void
YetiRuntime::register_allocation(Tensor* tensor, size_t size)
{
    tensor_allocations_[tensor] += size;
}

void
YetiRuntime::register_tensor(Tensor* tensor)
{
    tensor_allocations_[tensor] = 0;
}

void
YetiRuntime::unregister_tensor(Tensor* tensor)
{
    std::map<Tensor*,size_t>::iterator it = tensor_allocations_.find(tensor);
    tensor_allocations_.erase(it);
}

void
YetiRuntime::set_nindex_min_distr(usi nindex)
{
    nindex_min_distr_ = nindex;
    use_auto_distr_ = true;
}

bool
YetiRuntime::use_auto_distr()
{
    return use_auto_distr_;
}
