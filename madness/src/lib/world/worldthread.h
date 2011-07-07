/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/
#ifndef MADNESS_WORLD_WORLDTHREAD_H__INCLUDED
#define MADNESS_WORLD_WORLDTHREAD_H__INCLUDED

/// \file worldthread.h
/// \brief Implements Dqueue, Thread, ThreadBase and ThreadPool

#include <world/dqueue.h>
#include <vector>
#include <cstddef>
#include <pthread.h>

#ifndef _SC_NPROCESSORS_CONF
// Old macs don't have necessary support thru sysconf to determine the
// no. of processors so must use sysctl
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

namespace madness {

    // Forward decl.
    class Barrier;
    class ThreadPool;
    class WorldTaskQueue;
    class AtomicInt;
    void error(const char *msg);

    /// Simplified thread wrapper to hide pthread complexity

    /// If the thread is using any of the object state you cannot
    /// delete the object until the thread has terminated.
    ///
    /// The cleanest solution is to put the object on the heap and
    /// have the run method "delete this" at its end.
    class ThreadBase {
        friend class ThreadPool;
        static bool bind[3];
        static int cpulo[3];
        static int cpuhi[3];

        static void* main(void* self);

        int pool_num; ///< Stores index of thread in pool or -1
        pthread_t id;

        void set_pool_thread_index(int i) { pool_num = i; }

    public:

        /// Default constructor ... must invoke \c start() to actually begin the thread.
        ThreadBase() : pool_num(-1) { }

        virtual ~ThreadBase() { }

        /// You implement this to do useful work
        virtual void run() = 0;

        /// Start the thread running
        void start();

        /// A thread can call this to terminate its execution
        static void exit() { pthread_exit(0); }

        /// Get the pthread id of this thread (if running)
        const pthread_t& get_id() const { return id; }

        /// Get index of thread in ThreadPool (0,...,nthread-1) or -1 if not in ThreadPool
        int get_pool_thread_index() const { return pool_num; }

        /// Cancel this thread
        int cancel() const { return pthread_cancel(get_id()); }


        /// Get no. of actual hardware processors
        static int num_hw_processors();

        /// Specify the affinity pattern or how to bind threads to cpus
        static void set_affinity_pattern(const bool bind[3], const int cpu[3]);

        static void set_affinity(int logical_id, int ind=-1);
    }; // class ThreadBase

    /// Simplified thread wrapper to hide pthread complexity
    class Thread : public ThreadBase {
        void* (*f)(void *);
        void* args;

        void run() { f(args); }

    public:
        /// Default constructor ... must invoke \c start() to actually begin the thread.
        Thread() : f(0), args(0) {};

        /// Create a thread and start it running f(args)
        Thread(void* (*f)(void *), void* args=0)
                : f(f), args(args) {
            ThreadBase::start();
        }

        void start(void* (*f)(void *), void* args=0) {
            this->f = f;
            this->args = args;
            ThreadBase::start();
        }

        virtual ~Thread() {}
    };


    /// Contains attributes of a task

    /// \c generator : Setting this hints that a task will produce
    /// additional tasks and is used by the scheduler to
    /// increase/throttle parallelism. The default is false.
    ///
    /// \c stealable : Setting this indicates that a task may be
    /// migrated to another process for dynamic load balancing.  The
    /// default value is false.
    ///
    /// \c highpriority : indicates a high priority task. The default value is false.
    ///
    /// \c nthread : indicates number of threads. 0 threads is interpreted as 1 thread
    /// for backward compatibility and ease of specifying defaults. The default value
    /// is 0 (==1).
    class TaskAttributes {
        unsigned long flags;
    public:
    	static const unsigned long NTHREAD      = 0xff;          // Mask for nthread byte
        static const unsigned long GENERATOR    = 1ul<<8;        // Mask for generator bit
        static const unsigned long STEALABLE    = GENERATOR<<1;  // Mask for stealable bit
        static const unsigned long HIGHPRIORITY = GENERATOR<<2;  // Mask for priority bit

        TaskAttributes(unsigned long flags = 0) : flags(flags) {}

        TaskAttributes(const TaskAttributes& attr) : flags(attr.flags) {}

        virtual ~TaskAttributes() {}

        bool is_generator() const { return flags&GENERATOR; }

        bool is_stealable() const { return flags&STEALABLE; }

        bool is_high_priority() const { return flags&HIGHPRIORITY; }

        void set_generator(bool generator_hint) {
            if (generator_hint)
                flags |= GENERATOR;
            else
                flags &= ~GENERATOR;
        }

        void set_stealable(bool stealable) {
            if (stealable) flags |= STEALABLE;
            else flags &= ~STEALABLE;
        }

        void set_highpriority(bool hipri) {
            if (hipri)
                flags |= HIGHPRIORITY;
            else
                flags &= ~HIGHPRIORITY;
        }

        /// Are you sure this is what you want to call?

        /// Only call this for a \c TaskAttributes that is \em not a base class
        /// of a task object.
        ///
        /// If you are trying to set the number of threads in an \em existing
        /// task you should call \c TaskInterface::set_nthread() instead.
        /// No doubt there is some virtual/protected/something voodoo to prevent
        /// you from doing harm.
        void set_nthread(int nthread) {
            MADNESS_ASSERT(nthread>=0 && nthread<256);
            flags = (flags & (~NTHREAD)) | (nthread & NTHREAD);
        }

        int get_nthread() const {
        	int n = flags & NTHREAD;
        	if (n == 0)
        	    n = 1;
        	return n;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar & flags;
        }

        static TaskAttributes generator() {
            return TaskAttributes(GENERATOR);
        }

        static TaskAttributes hipri() {
            return TaskAttributes(HIGHPRIORITY);
        }

        static TaskAttributes multi_threaded(int nthread) {
            TaskAttributes t;
            t.set_nthread(nthread);
            return t;
        }
    };

    /// Used to pass info about thread environment into users task
    class TaskThreadEnv {
        const int _nthread; //< No. of threads collaborating on task
        const int _id;      //< Id of this thread (0,...,nthread-1)
        Barrier* _barrier;  //< Pointer to shared barrier, null if single thread

    public:
        TaskThreadEnv(int nthread, int id, Barrier* barrier)
            : _nthread(nthread), _id(id), _barrier(barrier)
        {}

        int nthread() const {return _nthread;}

        int id() const {return _id;}

        bool barrier() const;
    };


    /// Lowest level task interface

    /// The pool invokes run_multi_threaded() that does any necessary
    /// setup for multiple threads and then invokes the users \c run method.
    class PoolTaskInterface : public TaskAttributes {
    	friend class ThreadPool;

    private:
        Barrier* barrier;     //< Barrier, only allocated for multithreaded tasks
    	AtomicInt count;  //< Used to count threads as they start

    	/// Returns true for the one thread that should invoke the destructor
    	bool run_multi_threaded();

    public:
        PoolTaskInterface()
            : TaskAttributes()
            , barrier(0)
        {
            count = 0;
        }

        explicit PoolTaskInterface(const TaskAttributes& attr);

        /// Call this to reset the number of threads before the task is submitted

        /// Once a task has been constructed /c TaskAttributes::set_nthread()
        /// is insufficient because a multithreaded task includes a
        /// barrier that needs to know the number of threads.
        void set_nthread(int nthread);

        /// Override this method to implement a multi-threaded task

        /// \c info.nthread() will be the number of threads collaborating on this task
        ///
        /// \c info.id() will be the index of the current thread \c id=0,...,nthread-1
        ///
        /// \c info.barrier() will be a barrier for all of the threads, and returns
        /// \c true for the last thread to enter the barrier (other threads get false)
        virtual void run(const TaskThreadEnv& info) = 0;

        virtual ~PoolTaskInterface();
    };

    /// A no-op task used for various purposes
    class PoolTaskNull : public PoolTaskInterface {
    public:
        void run(const TaskThreadEnv& /*info*/) {}
        virtual ~PoolTaskNull() {}
    };


    /// A singleton pool of threads for dynamic execution of tasks.

    /// YOU MUST INSTANTIATE THE POOL WHILE RUNNING WITH JUST ONE THREAD
    class ThreadPool {
    private:
        friend class WorldTaskQueue;
        Thread *threads;              ///< Array of threads
        DQueue<PoolTaskInterface*> queue; ///< Queue of tasks
        int nthreads;		  ///< No. of threads
        volatile bool finish;              ///< Set to true when time to stop
        AtomicInt nfinished;

        static ThreadPool* instance_ptr;

        /// The constructor is private to enforce the singleton model
        ThreadPool(int nthread=-1);

        ThreadPool(const ThreadPool&);           // Verboten
        void operator=(const ThreadPool&);       // Verboten

        /// Get number of threads from the environment
        int default_nthread();

        /// Run next task ... returns true if one was run ... blocks if wait is true
        bool run_task(bool wait);

        bool run_tasks(bool wait);

        void thread_main(Thread* thread);

        /// Forwards thread to bound member function
        static void* pool_thread_main(void *v);


        /// Return a pointer to the only instance constructing as necessary
        static ThreadPool* instance(int nthread=-1);


    public:
        /// Please invoke while in single threaded environment
        static void begin(int nthread=-1);

        static void end();

        /// Add a new task to the pool
        static void add(PoolTaskInterface* task);

        template <typename opT>
        void scan(opT& op) {
            queue.scan(op);
        }


        /// Add a vector of tasks to the pool
        static void add(const std::vector<PoolTaskInterface*>& tasks);

        /// An otherwise idle thread can all this to run a task

        /// Returns true if one was run
        static bool run_task();

        /// Returns number of threads in the pool
        static std::size_t size();

        /// Returns number of tasks in the queue
        static std::size_t queue_size();

        /// Returns queue statistics
        static const DQStats& get_stats();

        ~ThreadPool() {}
    };

}

#endif // MADNESS_WORLD_WORLDTHREAD_H__INCLUDED
