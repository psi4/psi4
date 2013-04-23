#include "taskqueue.h"
#include "sortimpl.h"
#include "exception.h"
#include "thread.h"
#include "runtime.h"
#include "contraction.h"
#include "tensorblock.h"
#include "tensor.h"

#include <algorithm>

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

using namespace yeti;
using namespace std;


#define DYNAMIC_LOAD_BALANCE 0


#define wtf //if (Malloc<TensorBlock>::get_object(43)->is_locked()) yeti_throw(SanityCheckError,"")

TaskQueue* GlobalQueue::queue_ = 0;
bool* GlobalQueue::node_acks_ = 0;
bool GlobalQueue::local_run_complete_ = true;
uli* GlobalQueue::remote_task_buffer_ = 0;
LoadBalanceQueue* GlobalQueue::head_queue_ = 0;
LoadBalanceQueue* GlobalQueue::tail_queue_ = 0;
ThreadLock* GlobalQueue::lock_ = 0;
double TaskQueue::thread_times_[1000];
static uli run_count = 0;

#define NQUEUES_LOAD_BALANCE 100
static LoadBalanceQueue load_balance_queues[NQUEUES_LOAD_BALANCE];
static uli load_balance_queue_num = 0;

static GlobalQueue global_queue;

#define MAX_NTASKS 4000000
#define MAX_NOWNERS 1000000

Task::Task()
    : 
    next(0)
{
}

Task::~Task()
{
}


TaskQueue::TaskQueue()
    :
    parent_(0),
    unexpected_head_(0),
    unexpected_tail_(0),
    active_(false)
{
}

TaskQueue::~TaskQueue()
{
}

// random generator function:
ptrdiff_t myrandom (ptrdiff_t i) { return rand()%i;}

// pointer object to it:
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;

#include <time.h>

void
TaskQueue::configure()
{
    active_ = true;
}

void
TaskQueue::clear()
{
#if YETI_SANITY_CHECK
    if (unexpected_tail_ || unexpected_head_)
        yeti_throw(SanityCheckError, "Unexpected tasks were not all run");
#endif
}

Task*
TaskQueue::get_next_task(uli threadnum) 
{
    lock();
    if (unexpected_head_)
    {
        Task* task = unexpected_head_;
        unexpected_head_ = 0;
        unexpected_tail_ = 0;
        unlock();
        return task;
    }

    Task* task = parent_->get_next_task();
    if (!task)
        active_ = false;
    unlock();

    return task;
}

bool
TaskQueue::worker_threads_active() const
{
    return active_;
}

void
TaskQueue::add_unexpected_task(Task* task)
{
    lock();
    if (active_)
    {
        if (unexpected_tail_)
        {
            unexpected_tail_->next = task;
            unexpected_tail_ = task;
        }
        else
        {
            unexpected_head_ = task;
            unexpected_tail_ = task;
        }
        unlock();
    }
    else //no more threads looking for tasks - I need to run it
    {
        unlock();
        uli threadnum = YetiRuntime::get_thread_number();
        task->run(threadnum);
    }
}

void
TaskQueue::print(std::ostream& os) const
{
    os << "Task Queue" << endl;
}

void
TaskQueue::reset_thread_times()
{
    for (uli i=0; i < YetiRuntime::nthread(); ++i)
    {
        thread_times_[i] = 0;
    }
}

void
TaskQueue::increment_thread_time(uli thr, double time)
{
    thread_times_[thr] += time;
}

double
TaskQueue::get_max_thread_time()
{
    double time = thread_times_[0];
    for (uli i=0; i < YetiRuntime::nthread(); ++i)
    {
        if (time < thread_times_[i])
            time = thread_times_[i];
    }
    return time;
}

GlobalQueue::GlobalQueue()
{
}

void
GlobalQueue::configure()
{
    local_run_complete_ = false;
    queue_->configure();
}

#define REMOTE_TASK_BUFFER_SIZE 100000
void
GlobalQueue::init()
{
    if (queue_)
        delete queue_;
    queue_ = new TaskQueue;
    queue_->incref();

    if (node_acks_)
        delete node_acks_;
    node_acks_ = new bool[YetiRuntime::num_nodes_task_group()];
    ::memset(node_acks_, 0, YetiRuntime::num_nodes_task_group() * sizeof(bool));

    if (remote_task_buffer_ == 0)
        remote_task_buffer_ = new uli[REMOTE_TASK_BUFFER_SIZE];

    if (lock_ == 0)
        lock_= new pThreadLock;
}

bool
GlobalQueue::is_node_complete(uli group_node_number)
{
    return node_acks_[group_node_number];
}

void
GlobalQueue::set_node_ack(uli global_node_number)
{
    uli group_number = YetiRuntime::group_node_number(global_node_number);
    node_acks_[group_number] = 1;
}

bool
GlobalQueue::local_run_complete()
{
    return local_run_complete_;
}

void
GlobalQueue::finalize()
{
    if (queue_)
        delete queue_;

    if (node_acks_)
        delete node_acks_;

    if (remote_task_buffer_)
        delete remote_task_buffer_;

    if (lock_)
        delete lock_;
}

void
GlobalQueue::run()
{


#if PRINT_TASK_TIMINGS
    cout << stream_printf("Running %ld tasks on node %ld\n",
                    queue_->ntasks(), YetiRuntime::me());
    cout.flush();
#endif

    ThreadGroup* thrgrp = YetiRuntime::get_thread_grp();

    for (uli i=0; i < YetiRuntime::nthread(); ++i)
    {
        Thread* thr = new TaskThread(i, queue_);
        thrgrp->add(thr);
    }

    double start = timer::Timer::getTime();

    TaskQueue::reset_thread_times();

    thrgrp->run();
    thrgrp->wait();
    queue_->clear();

    double max = TaskQueue::get_max_thread_time();
    cout << stream_printf("Total computation on node %ld is %8.4f\n", 
                    YetiRuntime::me(), max);

    /** At this point, no remote nodes should steal tasks from me */
    local_run_complete_ = true;
    
    double stop = timer::Timer::getTime();
    cout << stream_printf("Node task run %ld took %8.4f seconds \n",
                YetiRuntime::me(), stop - start);
    start = stop;


    if (YetiRuntime::use_dynamic_load_balancing())
    {
        for (uli plocal=0; plocal < YetiRuntime::num_nodes_task_group(); ++plocal)
        {
            uli pglobal = YetiRuntime::global_node_number(plocal);
            YetiRuntime::get_messenger()->send_task_notification(pglobal);
            while (!GlobalQueue::is_node_complete(plocal))
            {
                if (receive_dynamic_load_balance_queue())
                {
                    queue_->configure();
                    thrgrp->run();
                    thrgrp->wait();
                    queue_->clear();
                    YetiRuntime::get_messenger()->send_task_notification(pglobal);
                }
            }
        }

    }
    thrgrp->clear();

#if PRINT_TASK_TIMINGS
    cout << stream_printf("Node %ld run took %8.4f seconds idle %8.4f seconds and ran %ld tasks\n",    
                        YetiRuntime::me(), task_time, idle_time, queue_->task_run_count());
    cout.flush();
#endif


    /** this memset is safe since barrier follows */
    ::memset(node_acks_, 0, YetiRuntime::num_nodes_task_group() * sizeof(bool));

    queue_->parent_->finalize();
    queue_->parent_ = 0;

    stop = timer::Timer::getTime();
    cout << stream_printf("Node finalize task queue %ld took %8.4f seconds \n",
                YetiRuntime::me(), stop - start);
    start = stop;
}

TaskQueue*
GlobalQueue::get_task_queue()
{
    return queue_;
}

void
GlobalQueue::add(TaskParent* parent)
{
    if (queue_->parent_)
        yeti_throw(SanityCheckError, "Already have task parent");

    queue_->parent_ = parent;
}

void
GlobalQueue::print(std::ostream& os)
{
    queue_->print(os);
}


Task*
GlobalQueue::get_dynamic_load_balance_segment(uli& ntasks, uli& nentries, uli* data)
{
    if (local_run_complete_)
    {
        return 0;
    }
    return queue_->get_dynamic_load_balance_segment(ntasks, nentries, data);
}

#define NTASKS_MIN_LOAD_BALANCE 50
#define MIN_REMAINING_PERCENT_LOAD_BALANCE 10
#define LOAD_BALANCE_TASK_PERCENT 25

Task*
TaskQueue::get_dynamic_load_balance_segment(uli& ntasks, uli& nentries, uli* data)
{
    yeti_throw(SanityCheckError, "No dynamic load balancing yet");
    return 0;
}

LoadBalanceQueue*
GlobalQueue::get_dynamic_load_balance_queue()
{
    if (load_balance_queue_num == NQUEUES_LOAD_BALANCE)
        load_balance_queue_num = 0;

    LoadBalanceQueue* queue = &load_balance_queues[load_balance_queue_num];
    return queue;
}

void
GlobalQueue::add_dynamic_load_balance_queue(
    LoadBalanceQueue* queue
)
{
    lock_->lock();
    queue->next = 0;
    if (head_queue_ == 0)
    {
        head_queue_ = queue;
        tail_queue_ = queue;
    }
    else
    {
        tail_queue_->next = queue;
        tail_queue_ = queue;
    }

    ++load_balance_queue_num;
    lock_->unlock();
}

bool
GlobalQueue::receive_dynamic_load_balance_queue()
{
    if (head_queue_ == 0)
    {
        return false;
    }

    lock_->lock();
    uli entry = 0;
    const uli* dataptr = head_queue_->queue.data;
    uli nentries = head_queue_->nentries;
    head_queue_ = head_queue_->next;
    if (head_queue_ == 0)
        tail_queue_ = 0;
    lock_->unlock();
    

    YetiRuntime::start_timer("receive load balance");
    while(entry < nentries)
    {
        yeti_throw(SanityCheckError, "Dynamic load balancing not yet allowed");
    }
    YetiRuntime::stop_timer("receive load balance");
    return true;
}

TaskThread::TaskThread(
    uli threadnum,
    const TaskQueuePtr &queue
)
    : Thread(threadnum),
    queue_(queue)
{

}

void
TaskThread::run()
{
    Task* current_task = queue_->get_next_task(threadnum_);
    if (current_task)
        current_task->prefetch(threadnum_);
    Task* next = queue_->get_next_task(threadnum_);
    while (current_task)
    {
        Task* old = current_task;
        if (current_task->next)
        {
            current_task->next->prefetch(threadnum_);
            current_task->run(threadnum_);
            current_task = current_task->next;
        }
        else if (next)
        {
            next->prefetch(threadnum_);
            current_task->run(threadnum_);
            current_task = next;
            next = queue_->get_next_task(threadnum_);
        }
        else
        {
            current_task->run(threadnum_);
            current_task = 0; //nothing to do!
        }
        delete old;
    }
}

