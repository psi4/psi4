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

DECLARE_PARENT_MALLOC(Task);
DECLARE_SUB_MALLOC(Task,ContractionTask);

#define wtf //if (Malloc<TensorBlock>::get_object(43)->is_locked()) raise(SanityCheckError,"")

TaskQueue* GlobalQueue::queue_ = 0;

static GlobalQueue global_queue;

#define MAX_NTASKS 4000000
#define MAX_NOWNERS 1000000

#define DYNAMIC_LOAD_BALANCE_THREADS 0
#define RANDOMIZE_QUEUE 0
#define SORT_QUEUE 1

Task::~Task()
{
}

TaskQueue::TaskQueue()
    :
  workspace_(new TaskWorkspace)
{
}

TaskQueue::~TaskQueue()
{
    delete workspace_;
}

// random generator function:
ptrdiff_t myrandom (ptrdiff_t i) { return rand()%i;}

// pointer object to it:
ptrdiff_t (*p_myrandom)(ptrdiff_t) = myrandom;

#include <time.h>

void
TaskQueue::configure()
{
    if (workspace_->configured)
        return;

    uli ntasks = workspace_->ntasks;
    int* sorted_indices = workspace_->sorted_indices;
    TaskOwner** task_owners = workspace_->task_owners;
    uli* ntasks_per_owner = workspace_->ntasks_per_owner;
    uli* task_offsets = workspace_->task_offsets;
    Task** tasks = workspace_->tasks;
    void* sorted = workspace_->sorted;


    yeti::quicksort<uli>(
        workspace_->owner_numbers,
        sorted_indices,
        ntasks
    );


    TaskOwner** sorted_owners = reinterpret_cast<TaskOwner**>(sorted);
    for (uli i=0; i < ntasks; ++i)
        sorted_owners[i] = task_owners[sorted_indices[i]];
    for (uli i=0; i < ntasks; ++i)
        task_owners[i] = sorted_owners[i];


    Task** sorted_tasks = reinterpret_cast<Task**>(sorted);
    for (uli i=0; i < ntasks; ++i)
        sorted_tasks[i] = tasks[sorted_indices[i]];
    for (uli i=0; i < ntasks; ++i)
        tasks[i] = sorted_tasks[i];

    int last = -1;
    int owner = -1; //I actually mean to do this
    for (uli i=0; i < ntasks; ++i)
    {
        uli next = task_owners[i]->get_task_owner_number();
        if (next != last)
        {
            ++owner;
            ntasks_per_owner[owner] = 1;
            task_offsets[owner] = i;
        }
        else
        {
            ++ntasks_per_owner[owner];
        }
        last = next;
    }
    uli nowners = owner + 1; //off by one error

    //do a sanity check
    uli ntask_check = 0;
    for (uli i=0; i < nowners; ++i)
        ntask_check += ntasks_per_owner[i];

    if (ntask_check != ntasks)
        raise(SanityCheckError, "contraction misconfigured");

    workspace_->nowners = nowners;
    workspace_->configured = true;

    //set the run number for the calculation
    run_number_ = 0;

}

void
TaskQueue::run()
{
    configure();

    ThreadGroup* thrgrp = YetiRuntime::get_thread_grp();

    for (uli i=0; i < YetiRuntime::nthread_compute(); ++i)
    {
        Thread* thr = new TaskThread(i, this);
        thrgrp->add(thr);
    }

    thrgrp->run();
    thrgrp->wait();
    thrgrp->clear();

    clear();
}

void
TaskQueue::clear()
{
    workspace_->clear();
    workspace_->configured = false;

    std::list<TaskParentPtr>::iterator it(parents_.begin());
    std::list<TaskParentPtr>::iterator stop(parents_.end());
    for ( ; it != stop; ++it)
        (*it)->finalize();
    parents_.clear();
}

const uli*
TaskQueue::ntasks_per_owner() const
{
    return workspace_->ntasks_per_owner;
}

const uli*
TaskQueue::task_offsets() const
{
    return workspace_->task_offsets;
}

Task**
TaskQueue::get_tasks() const
{
    return workspace_->tasks;
}

TaskOwner**
TaskQueue::get_owners() const
{
    return workspace_->task_owners;
}

uli
TaskQueue::nowners()
{
    return workspace_->nowners;
}

uli
TaskQueue::ntasks()
{
    return workspace_->ntasks;
}

void
TaskQueue::add(TaskOwner* owner, Task* task)
{
    if (workspace_->ntasks >= MAX_NTASKS)
    {
        raise(SanityCheckError, "not enough tasks allocated for contraction");
    }

    if (workspace_->nowners >= MAX_NOWNERS)
    {
        raise(SanityCheckError, "not enough owner slots allocated for contraction");
    }

    workspace_->task_owners[workspace_->ntasks] = owner;
    workspace_->owner_numbers[workspace_->ntasks] = owner->get_task_owner_number();
    workspace_->tasks[workspace_->ntasks] = task;
    workspace_->configured = false;
    ++workspace_->ntasks;
}

void
TaskQueue::add(const TaskParentPtr& parent)
{
    parents_.push_back(parent);
}

void
TaskQueue::print(std::ostream& os)
{
    for (uli i=0; i < workspace_->ntasks; ++i)
    {
        workspace_->tasks[i]->print(os);
        os << std::endl;
    }
}

Task*
TaskQueue::pop()
{
    if (run_number_ == workspace_->ntasks)
        return 0;

    Task* task = workspace_->tasks[run_number_];
    ++run_number_;
    return task;
}

void
TaskQueue::get_next(
    Task **&taskstart,
    uli &ntasks
)
{
    lock_->lock();
    if (run_number_ == workspace_->nowners)
    {
        taskstart = 0;
        ntasks = 0;
        lock_->unlock();
        return;
    }

    taskstart = workspace_->tasks + workspace_->task_offsets[run_number_];
    ntasks = workspace_->ntasks_per_owner[run_number_];

    ++run_number_;

    lock_->unlock();
}

Task**
TaskQueue::get_task_list_for_owner(uli owner) const
{
    return workspace_->tasks + workspace_->task_offsets[owner];
}

uli
TaskQueue::get_ntasks_for_owner(uli owner) const
{
    return workspace_->ntasks_per_owner[owner];
}

TaskWorkspace::TaskWorkspace()
    :
    sorted_indices(0),
    ntasks(0),
    nowners(0),
    sorted(0),
    task_offsets(0),
    tasks(0),
    ntasks_per_owner(0),
    configured(false)
{
    task_owners = new TaskOwner*[MAX_NTASKS];
    sorted_indices = new int[MAX_NTASKS];
    tasks = new Task*[MAX_NTASKS];
    ntasks_per_owner = new uli[MAX_NOWNERS];
    task_offsets = new uli[MAX_NOWNERS];
    owner_numbers = new uli[MAX_NTASKS];
    sorted = YetiRuntime::malloc(MAX_NTASKS * sizeof(void*));
}

TaskWorkspace::~TaskWorkspace()
{
    delete[] task_owners;
    delete[] owner_numbers;
    delete[] sorted_indices;
    delete[] tasks;
    delete[] task_offsets;
    delete[] ntasks_per_owner;
    YetiRuntime::free(sorted, MAX_NTASKS * sizeof(void*));
}

void
TaskWorkspace::clear()
{
    for (uli i=0; i < ntasks; ++i)
        delete tasks[i];

    ::memset(task_owners, 0, sizeof(TaskOwner*) * ntasks);
    ::memset(tasks, 0, sizeof(TaskOwner*) * ntasks);
    ntasks = 0;
    nowners = 0;
}

GlobalQueue::GlobalQueue()
{
}

void
GlobalQueue::configure()
{
    queue_->configure();
}

void
GlobalQueue::init()
{
    if (queue_)
        delete queue_;

    queue_ = new TaskQueue;
    queue_->incref();
}

void
GlobalQueue::finalize()
{
    if (queue_)
        delete queue_;
}

void
GlobalQueue::run()
{
    queue_->run();
}

void
GlobalQueue::clear()
{
    queue_->clear();
}

uli
GlobalQueue::nowners()
{
    return queue_->nowners();
}

uli
GlobalQueue::ntasks()
{
    return queue_->ntasks();
}

void
GlobalQueue::add(TaskOwner* owner, Task* task)
{
    queue_->add(owner, task);
}

void
GlobalQueue::add(const TaskParentPtr& parent)
{
    queue_->add(parent);
}

void
GlobalQueue::print(std::ostream& os)
{
    queue_->print(os);
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
#if DYNAMIC_LOAD_BALANCE_THREADS
    Task** tasklist = 0;
    uli ntasks = 0;
    queue_->get_next(tasklist, ntasks);
    while (tasklist)
    {
        for (uli i=0; i < ntasks; ++i)
        {
            tasklist[i]->run(threadnum_);
        }
        queue_->get_next(tasklist, ntasks);
    }
#else

    Task** tasks = queue_->get_tasks();
    TaskOwner** owners = queue_->get_owners();
    uli ntasks_tot = queue_->ntasks();
    uli nowners = queue_->nowners();
    const uli* task_offsets = queue_->task_offsets();
    const uli* ntasks_per_owner = queue_->ntasks_per_owner();
    uli nthread = YetiRuntime::nthread_compute();
    uli taskstop = nowners < nthread ? 0 : nowners - nthread;

    if (nowners <= threadnum_ || ntasks_tot == 0) //no work for this thread to do
        return;

    /** Fetch the very first task */
    Task** first_task_subset = tasks + task_offsets[threadnum_];
    Task* first_initial_task = first_task_subset[0];
    first_initial_task->out_of_core_prefetch();

    uli taskidx = threadnum_;
    for ( ; taskidx < taskstop; taskidx += nthread)
    {
        Task** current_task_subset = tasks + task_offsets[taskidx];
        Task* current_initial_task = current_task_subset[0];

        Task** next_task_subset = tasks + task_offsets[taskidx + nthread];
        Task* next_initial_task = next_task_subset[0];
        next_initial_task->out_of_core_prefetch();

        run_task_subset(current_task_subset, ntasks_per_owner[taskidx], threadnum_);
    }

    Task** last_task_subset = tasks + task_offsets[taskidx];

    run_task_subset(last_task_subset, ntasks_per_owner[taskidx], threadnum_);

#endif
}

void
TaskThread::run_task_subset(Task** tasks, uli ntasks, uli threadnum)
{
    Task* prev_task = 0;
    for (uli taskidx=0; taskidx < ntasks; ++taskidx)
    {
        Task* current_task = tasks[taskidx];
        current_task->in_core_prefetch(prev_task);
        current_task->run(threadnum);
        prev_task = current_task;
    }
    tasks[0]->finalize_task_subset();
}


TaskOwner::TaskOwner()
    :
    task_owner_number_(0)
{
}

uli
TaskOwner::get_task_owner_number() const
{
    return task_owner_number_;
}