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
DECLARE_SUB_MALLOC(ThreadLock,MultiThreadLock);

list<ThreadWorkspaceAllocator*>* ThreadEnvironment::allocators_ = 0;

/** Program termination variables */
static sigset_t pthread_sigset;

void*
pthread_run(void* thrptr)
{
    int s = pthread_sigmask(SIG_BLOCK, &pthread_sigset, NULL);
    if (s != 0)
    {
        cerr << "Unable to set signal mask on thread" << endl;
    }

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
        yeti_throw(SanityCheckError, "Unlocking mutex that I don't own");
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

void
ThreadGroup::thread_crash()
{
}

pThreadGroup::pThreadGroup(uli nthread)
    : ThreadGroup(nthread),
    pthreads_(new pthread_t[nthread]),
    attrs_(new pthread_attr_t[nthread]),
    running_(false)
{
    /* Block SIGINT; other threads created by main() will inherit
    a copy of the signal mask. */
    sigemptyset(&pthread_sigset);
    sigaddset(&pthread_sigset, SIGQUIT);
    sigaddset(&pthread_sigset, SIGTERM);
    sigaddset(&pthread_sigset, SIGINT);
    sigaddset(&pthread_sigset, SIGKILL);
    sigaddset(&pthread_sigset, MAIN_PROCESS_BACKTRACE_SIGNAL);

    ::memset(threads_, 0, nthread * sizeof(void*));
    ::memset(pthreads_, 0, nthread * sizeof(void*));
    for (uli i=1; i < nthread; ++i)
    {
        int status = pthread_attr_init(&attrs_[i]);
        if (status != 0)
        {
            cerr << "Unable to init thread attributes " << status << endl;
            abort();
        }

#if USE_THREAD_AFFINITY
        uli offset = YetiRuntime::cpu_mask_num();
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(i + offset, &cpuset);

        pthread_attr_setaffinity_np(&attrs_[i], sizeof(cpu_set_t), &cpuset);
        //cout << stream_printf("Setting thread affinity of thread %ld on node %ld to %d\n", i, YetiRuntime::me(), offset + i);
        //cout.flush();

#endif


#if USE_DEFAULT_THREAD_STACK
#else
        size_t stack_size = YetiRuntime::get_thread_stack_size();
        void* stack_start = YetiRuntime::get_thread_stack(i);
        status = pthread_attr_setstack(&attrs_[i], stack_start, stack_size);
        if (status != 0)
        {
            cerr << "Unable to set stack address " << stack_start 
                << " stack size " << stack_size
                << " on thread " << i << " with error "
                << status << endl;
            abort();
        }
#endif
    }

#if USE_THREAD_AFFINITY
    cpu_set_t cpuset;
    uli offset = YetiRuntime::cpu_mask_num();
    CPU_ZERO(&cpuset);
    CPU_SET(offset, &cpuset);
    sched_setaffinity(0,sizeof(cpu_set_t), &cpuset);
    cout << stream_printf("Setting thread affinity of main process on node %ld to %d\n", YetiRuntime::me(), offset);
    cout.flush();
#endif
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
        yeti_throw(SanityCheckError, "too many threads added to thread group");
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
    running_ = true;

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
    running_ = false;
    YetiRuntime::set_threaded_compute(false);
}

void
pThreadGroup::thread_crash()
{
    if (!running_)
        return;

    for (uli t=1; t < nthread_; ++t)
    {
        int status = pthread_kill(pthreads_[t], THREAD_KILL_SIGNAL);
        if (status == EINVAL)
        {
            cerr << "Signal " << THREAD_KILL_SIGNAL << " is not a valid pThread kill signal" << endl;
        }
    }
    usleep(100);
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

MultiThreadLock::MultiThreadLock()
    : 
    thread_owner_(0),
    lock_count_(0)
{
}

void
MultiThreadLock::lock()
{
    uli threadnum = YetiRuntime::get_thread_number();
    if (lock_count_ > 0 && thread_owner_ == threadnum)
    {
        ++lock_count_;
    }
    else
    {
        lock_.lock();
        thread_owner_ = threadnum; //import to set things in this order
        lock_count_ = 1;
    }
}

bool
MultiThreadLock::trylock()
{
    uli threadnum = YetiRuntime::get_thread_number();
    if (lock_count_ > 0 && thread_owner_ == threadnum)
    {
        ++lock_count_;
        return true;
    }
    else
    {
        bool test = lock_.trylock();
        return test;
    }
}

void
MultiThreadLock::unlock()
{
    if (lock_count_ > 1)
    {
        --lock_count_;
    }
    else
    {
        lock_count_ = 0;
        lock_.unlock();
    }
}

void
MultiThreadLock::print(std::ostream &os) const
{
    os << "MultiThreadLock" << endl;
}

