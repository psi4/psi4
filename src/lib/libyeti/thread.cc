#include "threadimpl.h"
#include "exception.h"
#include "tensorblock.h"

#include <errno.h>

using namespace yeti;
using namespace std;

#ifdef redefine_size_t
#define size_t custom_size_t
#endif



DECLARE_PARENT_MALLOC(ThreadLock);
DECLARE_SUB_MALLOC(ThreadLock,NullThreadLock);
DECLARE_SUB_MALLOC(ThreadLock,pThreadLock);

list<ThreadWorkspaceAllocator*>* ThreadEnvironment::allocators_ = 0;


void*
pthread_run(void* thrptr)
{
    Thread* thr = static_cast<Thread*>(thrptr);
#if USE_DEFAULT_THREAD_STACK
    int x; void* ptr = &x;
    YetiRuntime::set_thread_stack(thr->get_thread_number(), ptr);
#endif
    thr->run();
    return 0;
}

void
ThreadEnvironment::allocate()
{
    if (allocators_ == 0)
        return; //no allocators were ever added

    list<ThreadWorkspaceAllocator*>::const_iterator it(allocators_->begin());
    list<ThreadWorkspaceAllocator*>::const_iterator stop(allocators_->end());
    for ( ; it != stop; ++it)
    {
        (*it)->allocate();
    }
}

void
ThreadEnvironment::deallocate()
{
    if (allocators_ == 0)
        return; //no allocators were ever added

    list<ThreadWorkspaceAllocator*>::const_iterator it(allocators_->begin());
    list<ThreadWorkspaceAllocator*>::const_iterator stop(allocators_->end());
    for ( ; it != stop; ++it)
    {
        (*it)->deallocate();
    }
}

void
ThreadEnvironment::add_allocator(ThreadWorkspaceAllocator *allocator)
{
    if (allocators_ == 0) //allocate
        allocators_ = new list<ThreadWorkspaceAllocator*>;
    allocators_->push_back(allocator);
}

ThreadLock::~ThreadLock()
{
}

pThreadLock::pThreadLock()
{
    int signal = pthread_mutexattr_init(&attr_);
    if (signal != 0)
    {
        cerr << "mutex attr init error " << signal << endl;
        abort();
    }

    signal = pthread_mutex_init(&mutex_, &attr_);
    if (signal != 0)
    {
        cerr << "mutex init error " << signal << endl;
        abort();
    }
}

pThreadLock::~pThreadLock()
{
    /** Ignore the signal for now since whatever fucktard wrote
      some of the pthread implementations doesn't know how turn off
      all of the locks. This often erroneously returns signal 16 EBUSY */
    int signal = pthread_mutex_destroy(&mutex_);
    //if (signal != 0)
    //{
    //    throw 20;
    //}

    signal = pthread_mutexattr_destroy(&attr_);
    if (signal != 0)
    {
        cerr << "mutex attr destroy error " << signal << endl;
        abort();
    }
}

void
pThreadLock::lock()
{
    int signal = pthread_mutex_lock(&mutex_);
    if (signal != 0)
    {
        cerr << "mutex lock error " << signal << " on " << &mutex_ << endl;
        abort();
    }
}

bool
pThreadLock::trylock()
{
    int signal = pthread_mutex_trylock(&mutex_);
    return signal == 0; //signal of 0 means lock successful
}

void
pThreadLock::unlock()
{
    int signal = pthread_mutex_unlock(&mutex_);
    if (signal != 0)
    {
        raise(SanityCheckError, "unlocking mutex that I don't own");
    }
}

pthread_mutex_t*
pThreadLock::mutex()
{
    return &mutex_;
}

void
pThreadLock::print(std::ostream &os) const
{
    os << "pThread Lock" << endl;
}

ThreadGroup::ThreadGroup(uli nthread)
  :
    nthread_max_(nthread),
    threads_(new Thread*[nthread]),
    nthread_(0)
{
}

ThreadGroup::~ThreadGroup()
{
    clear();
}

void
ThreadGroup::clear()
{
    for (uli i=0; i < nthread_; ++i)
    {
        delete threads_[i];
        threads_[i] = 0;
    }

    nthread_ = 0;
}

pThreadGroup::pThreadGroup(uli nthread)
    : ThreadGroup(nthread),
    pthreads_(new pthread_t[nthread]),
    attrs_(new pthread_attr_t[nthread])
{
    ::memset(threads_, 0, nthread * sizeof(void*));
    ::memset(pthreads_, 0, nthread * sizeof(void*));
    for (uli i=1; i < nthread; ++i)
    {
        int status = pthread_attr_init(&attrs_[i]);
        if (status != 0)
        {
            cerr << "unable to init thread attributes " << status << endl;
            abort();
        }

        //cpu_set_t cpuset;
        //CPU_ZERO(&cpuset);
        //CPU_SET(i, &cpuset);

        //pthread_attr_setaffinity_np(&attrs_[i], sizeof(cpu_set_t), &cpuset);


#if USE_DEFAULT_THREAD_STACK
#else
        size_t stack_size = YetiRuntime::get_thread_stack_size();
        void* stack_start = YetiRuntime::get_thread_stack(i);
        status = pthread_attr_setstack(&attrs_[i], stack_start, stack_size);
        if (status != 0)
        {
            cerr << "unable to set stack address " << stack_start 
                << " stack size " << stack_size
                << " on thread " << i << " with error "
                << status << endl;
            abort();
        }
#endif
    }
}

pThreadGroup::~pThreadGroup()
{
    for (uli i=0; i < nthread_; ++i)
    {
        pthread_attr_destroy(&attrs_[i]);
    }
    delete[] attrs_;
    delete[] threads_;
}

void
pThreadGroup::clear()
{
    ThreadGroup::clear();
}

void
pThreadGroup::add(Thread* thr)
{
    if (nthread_ >= nthread_max_)
    {
        raise(SanityCheckError, "too many threads added to thread group");
    }

    threads_[nthread_] = thr;
    ++nthread_;
}

void
pThreadGroup::run()
{
    YetiRuntime::set_threaded_compute(true);
    for (uli t=1; t < nthread_; ++t)
    {
        int status = pthread_create(&pthreads_[t], &attrs_[t],
                pthread_run, reinterpret_cast<void*>(threads_[t]));
        if (status != 0)
        {
            cerr << "thread " << t << " not created properly" << endl;
            abort();
        }
        //void* retval;
        //status = pthread_join(pthreads_[t], &retval);
    }
    threads_[0]->run();
}

void
pThreadGroup::wait()
{
    //return;
    void* retval;
    for (uli t=1; t < nthread_; ++t)
    {
        int status = pthread_join(pthreads_[t], &retval);
        if (status != 0)
        {
            cerr << "thread " << t << " returned error status" << endl;
            abort();
        }

        if (retval != 0)
        {
            cerr << "thread " << t << " returned error status" << endl;
            abort();
        }
    }
    YetiRuntime::set_threaded_compute(false);
}


ThreadLock*
pThreadGroup::get_lock()
{
    ThreadLock* lock;
    /** If we have multiple compute threads or we have a communication thread */
    if (nthread_max_ > 1 || YetiRuntime::nproc() > 1)
    {
        lock = new pThreadLock;
    }
    else
    {
        lock = new pThreadLock;
        //lock = new NullThreadLock;
    }

    return lock;
}

void
NullThreadLock::lock()
{
    //do nothing
}

void
NullThreadLock::unlock()
{
}

bool
NullThreadLock::trylock()
{
    heisenfxn(NullThreadLock::trylock);
    return true; //always free
}

void
NullThreadLock::print(std::ostream &os) const
{
    os << "NullThreadLock" << endl;
}

Thread::Thread(uli threadnum)
    : threadnum_(threadnum)
{
}

uli
Thread::get_thread_number() const
{
    return threadnum_;
}

Thread::~Thread()
{
}


