#ifndef yeti_TASKQUEUE_H
#define yeti_TASKQUEUE_H

#include "class.h"
#include "thread.h"
#include "yetiobject.h"

#include "taskqueue.hpp"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

/**
    @class Task
*/
class Task :
    public smartptr::Countable,
    public Malloc<Task>
{

    public:
        virtual ~Task();

        virtual void run(uli threadnum) = 0;

        virtual void print(std::ostream& os = std::cout) const = 0;

        virtual void out_of_core_prefetch() = 0;

        virtual void in_core_prefetch(Task* prev_task) = 0;

        virtual void finalize_task_subset() = 0;

};

class TaskParent :
    public YetiRuntimeCountable
{
    public:
        virtual void finalize() = 0;
};

/**
    @class TaskOwner
    Dummy type used in task workspace for defining an array of pointers.
    This is only used for clarity in reading code.
*/
class TaskOwner {

    protected:
        uli task_owner_number_;

    public:
        TaskOwner();

        uli get_task_owner_number() const;

};

/**
    @class TaskWorkspace
    Set of arrays and values for defining a task queue. A task owner is
    usually the alpha tensor (see Contraction::alpha_tensor_).
*/
class TaskWorkspace {

    public:
        TaskWorkspace();

        ~TaskWorkspace();

        void clear();

        /**
            The number of unique task owners
        */
        uli nowners;

        /**
            The total number of tasks
        */
        uli ntasks;

        /**
            Array giving the number of tasks for the ith task owner
        */
        uli* ntasks_per_owner;

        /**
            Array giving the task offset number for the first task belonging
            to a given owner, assuming tasks have been sorted
        */
        uli* task_offsets;

        /**
            The list of task owners for each task
        */
        TaskOwner** task_owners;

        /**
            The complete list of tasks
        */
        Task** tasks;

        /**
            Workspace array for performing sort operations
        */
        uli* owner_numbers;

        /**
            Workspace array for performing sort operations
        */
        int* sorted_indices;

        void* sorted;

        /**
            Whether the values here have been configured
        */
        bool configured;


};


/**
    @class TaskQueue
*/
class TaskQueue :
    public YetiRuntimeCountable
{

    private:
        TaskWorkspace* workspace_;

        std::list<TaskParentPtr> parents_;

        /**
            The entry number of the task in the current run
        */
        uli run_number_;

    public:
        TaskQueue();

        ~TaskQueue();

        /**
            Sort all tasks based on task owners to optimize
            the task ordering
        */
        void configure();

        /**
            Configure and run all tasks. This spawns threads
            if YetiRuntime has been configured for more than one thread.
        */
        void run();

        /**
            Clear all tasks in the queue.
        */
        void clear();

        /**
            @return The complete task list
        */
        Task** get_tasks() const;

        TaskOwner** get_owners() const;

        /**
            @return
        */
        const uli* ntasks_per_owner() const;

        /**
            @return
        */
        const uli* task_offsets() const;

        /**
            NOT THREAD SAFE! Pop a task off the queue.
        */
        Task* pop();

        /**
            Add a task corresponding to a given owner.
        */
        void add(TaskOwner* owner, Task* task);

        void add(const TaskParentPtr& parent);

        void print(std::ostream& os = std::cout);

        uli get_ntasks_for_owner(uli owner) const;

        Task** get_task_list_for_owner(uli owner) const;

        /**
            @return
        */
        uli nowners();

        /**
            @return
        */
        uli ntasks();

        /**
            Thread safe. Lock the queue and return the beginning of set of tasks
            belong to a particular owner.
            @param taskstart Reference return of tasklist beginning
            @param ntasks The number of tasks belonging to a task owner
        */
        void get_next(
            Task**& taskstart,
            uli& ntasks
        );

};


/**
    @class GlobalQueue

    A task queue with static implementation for allowing task queues
    to be built that are persistent in memory
*/
class GlobalQueue :
        public smartptr::Countable
{

    private:
        static TaskQueue* queue_;

    public:
        GlobalQueue();

        static void configure();

        static void run();

        static void init();

        static void finalize();

        static void clear();

        static void add(TaskOwner* owner, Task* task);

        static void add(const TaskParentPtr& parent);

        static void print(std::ostream& os = std::cout);

        static uli nowners();

        static uli ntasks();

};

/**
    @class TaskThread
    Implementation of Thread for running tasks from a queue.
*/
class TaskThread :
        public Thread
{

    protected:
        TaskQueuePtr queue_;

    private:
        void run_task_subset(Task** tasks, uli ntasks, uli threadnum);

    public:
        /**
            @param threadnum
            @param queue
        */
        TaskThread(
            uli threadnum,
            const TaskQueuePtr& queue
        );

        void run();

};

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif // TASKQUEUE_H
