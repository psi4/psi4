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

  $Id: $
*/

/// \file worldthread.h
/// \brief Implements Dqueue, Thread, ThreadBase and ThreadPool

#include <world/worldthread.h>
#include <world/worldprofile.h>
#include <world/worldexc.h>
#include <world/print.h>
#include <world/worldpapi.h>
#include <world/safempi.h>
#include <world/atomicint.h>
#include <cstring>

namespace madness {

    int ThreadBase::cpulo[3];
    int ThreadBase::cpuhi[3];
    bool ThreadBase::bind[3];

    ThreadPool* ThreadPool::instance_ptr = 0;

    void* ThreadBase::main(void* self) {
#ifdef HAVE_PAPI
        begin_papi_measurement();
#endif

        try {
            ((ThreadBase*)(self))->run();
        }
        catch (const MPI::Exception& e) {
            //        print(e);
            error("caught an MPI exception");
        }
        catch (const madness::MadnessException& e) {
            print(e);
            error("caught a MADNESS exception");
        }
        catch (const char* s) {
            print(s);
            error("caught a string exception");
        }
        catch (const std::string& s) {
            print(s);
            error("caught a string (class) exception");
        }
        catch (const std::exception& e) {
            print(e.what());
            error("caught an STL exception");
        }
        catch (...) {
            error("caught unhandled exception");
        }

#ifdef HAVE_PAPI
        end_papi_measurement();
#endif
        return 0;
    }

    /// Start the thread running
    void ThreadBase::start() {
        pthread_attr_t attr;
        // Want detached thread with kernel scheduling so can use multiple cpus
        // RJH ... why not checking success/failure????
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_DETACHED);
#ifndef HAVE_IBMBGP
        pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
#endif
        int result = pthread_create(&id, &attr, &ThreadBase::main, (void *) this);
        if (result) MADNESS_EXCEPTION("failed creating thread", result);

        pthread_attr_destroy(&attr);
    }


    /// Get no. of actual hardware processors
    int ThreadBase::num_hw_processors() {
#ifdef _SC_NPROCESSORS_CONF
        int ncpu = sysconf(_SC_NPROCESSORS_CONF);
        if (ncpu <= 0) MADNESS_EXCEPTION("ThreadBase: set_affinity_pattern: sysconf(_SC_NPROCESSORS_CONF)", ncpu);
#elif defined(HC_NCPU)
        int mib[2]={CTL_HW,HW_NCPU}, ncpu;
        size_t len = sizeof(ncpu);
        if (sysctl(mib, 2, &ncpu, &len, NULL, 0) != 0)
            MADNESS_EXCEPTION("ThreadBase: sysctl(CTL_HW,HW_NCPU) failed", 0);
        std::cout << "NCPU " << ncpu << std::endl;
#else
        int ncpu=1;
#endif
        return ncpu;
    }

    /// Specify the affinity pattern or how to bind threads to cpus
    void ThreadBase::set_affinity_pattern(const bool bind[3], const int cpu[3]) {
        memcpy(ThreadBase::bind, bind, 3*sizeof(bool));
        memcpy(ThreadBase::cpulo, cpu, 3*sizeof(int));

        int ncpu = num_hw_processors();

        // impose sanity and compute cpuhi
        for (int i=0; i<3; ++i) {
            if (cpulo[i] < 0) cpulo[i] = 0;
            if (cpulo[i] >= ncpu) cpulo[i] = ncpu-1;

            if (i<2 && bind[i]) cpuhi[i] = cpulo[i];
            else cpuhi[i] = ncpu-1;

            //std::cout << "PATTERN " << i << " " << bind[i] << " " << cpulo[i] << " " << cpuhi[i] << std::endl;
        }
    }

    void ThreadBase::set_affinity(int logical_id, int ind) {
        if (logical_id < 0 || logical_id > 2) {
            std::cout << "ThreadBase: set_affinity: logical_id bad?" << std::endl;
            return;
        }

        if (!bind[logical_id]) return;

        // If binding the main or rmi threads the cpu id is a specific cpu.
        //
        // If binding a pool thread, the cpuid is the lowest cpu
        // to be used.
        //
        // If floating any thread it is floated from the cpuid to ncpu-1

        int lo=cpulo[logical_id], hi=cpuhi[logical_id];

        if (logical_id == 2) {
            if (ind < 0) {
                std::cout << "ThreadBase: set_affinity: pool thread index bad?" << std::endl;
                return;
            }
            if (bind[2]) {
                int nnn = hi-lo+1;
                lo += (ind % nnn);
                hi = lo;
            }
        }

#ifndef ON_A_MAC
        cpu_set_t mask;
        CPU_ZERO(&mask);
        for (int i=lo; i<=hi; ++i) CPU_SET(i,&mask);
        if (sched_setaffinity(0, sizeof(mask), &mask) == -1) {
            perror("system error message");
            std::cout << "ThreadBase: set_affinity: Could not set cpu Affinity" << std::endl;
        }
        //else {
        //    printf("managed to set affinity\n");
        //}
#endif
    }

    bool TaskThreadEnv::barrier() const {
        if (_nthread == 1)
            return true;
        else {
            MADNESS_ASSERT(_barrier);
            return _barrier->enter(_id);
        }
    }

    /// Returns true for the one thread that should invoke the destructor
    bool PoolTaskInterface::run_multi_threaded() {
        // As a thread enters this routine it increments the shared counter
        // to generate a unique id without needing any thread-local storage.
        // A downside is this does not preserve any relationships between thread
        // numbering and the architecture ... more work ahead.
        int nthread = get_nthread();
        if (nthread == 1) {
            run(TaskThreadEnv(1,0,0));
            return true;
        }
        else {
            int id = count++;
            volatile bool barrier_flag;
            barrier->register_thread(id, &barrier_flag);

            run(TaskThreadEnv(nthread, id, barrier));

            return barrier->enter(id);
        }
    }

    PoolTaskInterface::PoolTaskInterface(const TaskAttributes& attr)
        : TaskAttributes(attr)
        , barrier(attr.get_nthread()>1 ? new Barrier(attr.get_nthread()) : 0)

    {
        count = 0;
    }

    void PoolTaskInterface::set_nthread(int nthread) {
        if (nthread != get_nthread()) {
            TaskAttributes::set_nthread(nthread);
            delete barrier;
            if (nthread > 1)
                barrier = new Barrier(nthread);
            else
                barrier = 0;

        }
    }

    PoolTaskInterface::~PoolTaskInterface() {
        delete barrier;
    }

    /// The constructor is private to enforce the singleton model
    ThreadPool::ThreadPool(int nthread) : nthreads(nthread), finish(false) {
        nfinished = 0;
        instance_ptr = this;
        if (nthreads < 0) nthreads = default_nthread();
        //std::cout << "POOL " << nthreads << std::endl;

        try {
            if (nthreads > 0)
                threads = new Thread[nthreads];
            else
                threads = 0;
        }
        catch (...) {
            MADNESS_EXCEPTION("memory allocation failed", 0);
        }

        for (int i=0; i<nthreads; ++i) {
            threads[i].set_pool_thread_index(i);
            threads[i].start(pool_thread_main, (void *)(threads+i));
        }
    }

    /// Get number of threads from the environment
    int ThreadPool::default_nthread() {
        int nthread;
        int shift = 0;
        char *cnthread = getenv("MAD_NUM_THREADS");
        // MAD_NUM_THREADS is total no. of application threads whereas
        // POOL_NTHREAD is just the number in the pool (one less)
        if (cnthread) shift = 1;
        if (cnthread == 0) cnthread = getenv("POOL_NTHREAD");

        if (cnthread) {
            int result = sscanf(cnthread, "%d", &nthread);
            if (result != 1)
                MADNESS_EXCEPTION("POOL_NTHREAD is not an integer", result);
            nthread -= shift;
        }
        else {
            nthread = ThreadBase::num_hw_processors();
            if (nthread < 2)
                nthread = 2;
            nthread = nthread - 1; // One less than # physical processors
        }
        return nthread;
    }

    /// Run next task ... returns true if one was run ... blocks if wait is true
    bool ThreadPool::run_task(bool wait) {
        if (!wait && queue.empty()) return false;
        std::pair<PoolTaskInterface*,bool> t = queue.pop_front(wait);
        // Task pointer might be zero due to stealing
        if (t.second && t.first) {
            PROFILE_BLOCK(working);
            if (t.first->run_multi_threaded())         // What we are here to do
                delete t.first;
        }
        return t.second;
    }

    bool ThreadPool::run_tasks(bool wait) {
        static const int nmax=128; // WAS 100 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DEBUG
        PoolTaskInterface* taskbuf[nmax];
        int ntask = queue.pop_front(nmax, taskbuf, wait);
        for (int i=0; i<ntask; ++i) {
            PROFILE_BLOCK(working);
            if (taskbuf[i]) { // Task pointer might be zero due to stealing
                if (taskbuf[i]->run_multi_threaded()) {
                delete taskbuf[i];
                }
            }
        }
        return (ntask>0);
    }

    void ThreadPool::thread_main(Thread* thread) {
        PROFILE_MEMBER_FUNC(ThreadPool);
        thread->set_affinity(2, thread->get_pool_thread_index());

#define MULTITASK
#ifdef  MULTITASK
        while (!finish) {
            run_tasks(true);
        }
#else
        while (!finish) {
            run_task(true);
        }
#endif
        nfinished++;
    }

    /// Forwards thread to bound member function
    void* ThreadPool::pool_thread_main(void *v) {
        instance()->thread_main((Thread*)(v));
        return 0;
    }


    /// Return a pointer to the only instance constructing as necessary
    ThreadPool* ThreadPool::instance(int nthread) {
        if (!instance_ptr) instance_ptr = new ThreadPool(nthread);
        return instance_ptr;
    }

    void ThreadPool::begin(int nthread) {
        instance(nthread);
    }

    void ThreadPool::end() {
        if (!instance_ptr) return;
        instance()->finish = true;
        for (int i=0; i<instance()->nthreads; ++i) {
            add(new PoolTaskNull);
        }
        while (instance_ptr->nfinished != instance_ptr->nthreads);
    }

    /// Add a new task to the pool
    void ThreadPool::add(PoolTaskInterface* task) {
        if (!task) MADNESS_EXCEPTION("ThreadPool: inserting a NULL task pointer", 1);
        int nthread = task->get_nthread();
        // Currently multithreaded tasks must be shoved on the end of the q
        // to avoid a race condition as multithreaded task is starting up
        if (task->is_high_priority() && nthread==1) {
            instance()->queue.push_front(task);
        }
        else {
            instance()->queue.push_back(task, nthread);
        }
    }

    void ThreadPool::add(const std::vector<PoolTaskInterface*>& tasks) {
        typedef std::vector<PoolTaskInterface*>::const_iterator iteratorT;
        for (iteratorT it=tasks.begin(); it!=tasks.end(); ++it) {
            add(*it);
        }
    }

    /// An otherwise idle thread can all this to run a task

    /// Returns true if one was run
    bool ThreadPool::run_task() {
        return instance()->run_tasks(false);
    }

    /// Returns number of threads in the pool
    size_t ThreadPool::size() {
        return instance()->nthreads;
    }

    /// Returns number of tasks in the queue
    size_t ThreadPool::queue_size() {
        return instance()->queue.size();
    }

    /// Returns a map of the pthreads ids to the unique integer thread ids
    pthread_t* ThreadPool::thread_id() {
        int nthreads = ThreadPool::size();
        pthread_t* id = new pthread_t[nthreads+1];

        for (int i=0; i<nthreads+1; i++) {
            id[i] = instance()->threads[i].get_id();
        }

        return id;
    }

    /// Returns queue statistics
    const DQStats& ThreadPool::get_stats() {
        return instance()->queue.get_stats();
    }

} // namespace madness
