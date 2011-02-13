#include "tile.h"
#include "taskqueue.h"
#include "sortimpl.h"
#include "exception.h"
#include "thread.h"
#include "runtime.h"
#include "contraction.h"

#include <algorithm>

using namespace yeti;
using namespace std;

DECLARE_SUBMALLOC(Task,AccumulateTask);

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
    void* sorted = workspace_->sorted;
    uli* ntasks_per_owner = workspace_->ntasks_per_owner;
    uli* task_offsets = workspace_->task_offsets;
    Task** tasks = workspace_->tasks;

#if RANDOMIZE_QUEUE
    /* initialize random seed: */
    srand ( time(NULL) );
    for (uli i=0; i < ntasks; ++i)
        sorted_indices[i] = i;
    random_shuffle(sorted_indices, sorted_indices + ntasks, myrandom);
#elif SORT_QUEUE
    yeti::quicksort<TaskOwner*>(
        task_owners,
        sorted_indices,
        ntasks
    );
#else
    //just run through tasks without any sorting
    for (uli i=0; i < ntasks; ++i)
        sorted_indices[i] = i;
#endif

    TaskOwner** sorted_owners = reinterpret_cast<TaskOwner**>(sorted);

    for (uli i=0; i < ntasks; ++i)
        sorted_owners[i] = task_owners[sorted_indices[i]];

    TaskOwner* last = 0;
    TaskOwner** next = sorted_owners;
    uli owner = -1; //I actually mean to do this
    for (uli i=0; i < ntasks; ++i, ++next)
    {
        if (*next != last)
        {
            ++owner;
            ntasks_per_owner[owner] = 1;
            task_offsets[owner] = i;
        }
        else
        {
            ++ntasks_per_owner[owner];
        }

        last = *next;
    }
    uli nowners = owner + 1; //off by one error

    //do a sanity check
    uli ntask_check = 0;
    for (uli i=0; i < nowners; ++i)
        ntask_check += ntasks_per_owner[i];

    if (ntask_check != ntasks)
        raise(SanityCheckError, "contraction misconfigured");

    //now sort the tasks
    Task** sorted_tasks = reinterpret_cast<Task**>(sorted);
    for (uli i=0; i < ntasks; ++i)
        sorted_tasks[i] = tasks[sorted_indices[i]];

    for (uli i=0; i < ntasks; ++i)
        tasks[i] = sorted_tasks[i];

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

    for (uli i=0; i < YetiRuntime::nthread(); ++i)
    {
        Thread* thr = new TaskThread(i, this);
        thrgrp->add(thr);
    }

    thrgrp->run();
    thrgrp->wait();
    thrgrp->clear();


    //Task** tasks = reinterpret_cast<Task**>(workspace_->sorted);
    //for (uli i=0; i < workspace_->ntasks; ++i)
    //    tasks[i]->run();

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
TaskQueue::add(void* owner, Task* task)
{
    if (workspace_->ntasks >= MAX_NTASKS)
    {
        raise(SanityCheckError, "not enough tasks allocated for contraction");
    }

    if (workspace_->nowners >= MAX_NOWNERS)
    {
        raise(SanityCheckError, "not enough owner slots allocated for contraction");
    }

    workspace_->task_owners[workspace_->ntasks] = reinterpret_cast<TaskOwner*>(owner);
    workspace_->tasks[workspace_->ntasks] = task;

    ++workspace_->ntasks;

}

void
TaskQueue::add(const TaskParentPtr &parent)
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
    sorted = new TaskOwner*[MAX_NTASKS];
}

TaskWorkspace::~TaskWorkspace()
{
    delete[] task_owners;
    delete[] sorted_indices;
    delete[] tasks;
    delete[] task_offsets;
    delete[] ntasks_per_owner;
    delete[] sorted;
}

void
TaskWorkspace::clear()
{
    for (uli i=0; i < ntasks; ++i)
        delete tasks[i];

    ::memset(task_owners, 0, sizeof(TaskOwner*) * MAX_NOWNERS);
    ::memset(tasks, 0, sizeof(TaskOwner*) * MAX_NTASKS);
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
GlobalQueue::add(void* owner, Task* task)
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
    uli ntasks_tot = queue_->ntasks();
    uli nowners = queue_->nowners();
    const uli* task_offsets = queue_->task_offsets();
    const uli* ntasks_per_owner = queue_->ntasks_per_owner();
    uli nthread = YetiRuntime::nthread();


    for (uli i=0; i < nowners; ++i)
    {
        if (i % nthread != threadnum_)
            continue;

        uli ntasks = ntasks_per_owner[i];
        Task** tasklist = tasks + task_offsets[i];
        for (uli i=0; i < ntasks; ++i)
        {
            tasklist[i]->run(threadnum_);
        }
    }
#endif
}
