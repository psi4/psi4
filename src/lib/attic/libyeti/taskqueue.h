#ifndef yeti_TASKQUEUE_H
#define yeti_TASKQUEUE_H

#include "class.h"
#include "thread.h"
#include "yetiobject.h"
#include "messenger.hpp"

#include "taskqueue.hpp"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

#define NELEMENTS_TASK_QUEUE 50000
struct TaskQueueList { uli data[NELEMENTS_TASK_QUEUE]; };

/**
    @class Task

    Abstract base class for tasks
*/
class Task :
    public smartptr::Countable
{

    public:
        Task* next; 

        typedef enum { ContractionType=1 } task_type_t;

        Task();

        virtual ~Task();

        virtual void run(uli threadnum) = 0;

        virtual void print(std::ostream& os = std::cout) const = 0;

        virtual void prefetch(uli threadnum) = 0;

        /**
            Append info to a data array for sending this
            task to another node
        */
        virtual uli append_info(uli* data) const = 0;

};

class TaskParent :
    public YetiRuntimeCountable
{
    public:
        virtual void finalize() = 0;

        virtual Task* get_next_task() = 0;

};

/**
    @class TaskQueue
*/
class TaskQueue :
    public YetiRuntimeCountable
{

    private:
        friend class GlobalQueue;

        Task* unexpected_head_;

        Task* unexpected_tail_;

        TaskParent* parent_;

        static double thread_times_[1000];

        bool active_;

    public:

        TaskQueue();

        ~TaskQueue();

        void print(std::ostream& os = std::cout) const;

        void configure();

        void clear();

        static void reset_thread_times();

        static void increment_thread_time(uli thr, double time);

        static double get_max_thread_time();

        /**
            Thread safe. Lock the queue and return the beginning of set of tasks
            belong to a particular owner.
            @param taskstart Reference return of tasklist beginning
            @param ntasks The number of tasks belonging to a task owner
        */
        Task* get_next_task(uli threadnum);

        Task* get_dynamic_load_balance_segment(uli& ntasks, uli& nentries, uli* data);
    
        /**
            If in the processo of running tasks, append it to the queue.
            If already finished, immediately run the task.
        */
        void add_unexpected_task(Task* task);

        bool worker_threads_active() const;
};


/**
    @class GlobalQueue

    A task queue with static implementation for allowing task queues
    to be built that are persistent in memory
*/

struct LoadBalanceQueue {
   LoadBalanceQueue* next; 

   TaskQueueList queue; 

   uli nentries;
};

class GlobalQueue :
        public smartptr::Countable
{

    private:
        static TaskQueue* queue_;

        static bool* node_acks_;

        static bool local_run_complete_;

        static uli* remote_task_buffer_;


        static LoadBalanceQueue* head_queue_;
        
        static LoadBalanceQueue* tail_queue_;

        static ThreadLock* lock_;

    public:
        GlobalQueue();

        static void configure();

        static void run();

        static void init();

        static void finalize();

        static void add_dynamic_load_balance_queue(LoadBalanceQueue* queue);

        static LoadBalanceQueue* get_dynamic_load_balance_queue();

        static bool receive_dynamic_load_balance_queue();

        static bool local_run_complete();

        static bool is_node_complete(uli group_node_number);
        
        static void set_node_ack(uli global_node_number);

        static void add(TaskParent* parent);

        static void print(std::ostream& os = std::cout);

        static TaskQueue* get_task_queue();

        static Task* get_dynamic_load_balance_segment(uli& ntasks, uli& nentries, uli* data);

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
