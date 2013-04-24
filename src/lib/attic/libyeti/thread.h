#ifndef yeti_thread_h
#define yeti_thread_h
#include <list>

#include <pthread.h>

#include "class.h"
#include "mallocimpl.h"

#include "thread.hpp"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define USE_DEFAULT_THREAD_STACK 0

namespace yeti {

class ThreadWorkspaceAllocator {

    public:
        virtual void allocate() = 0;

        virtual void deallocate() = 0;

};

class ThreadLock :
    public Malloc<ThreadLock>
{

    public:
        virtual ~ThreadLock();

        virtual void lock() = 0;

        virtual void unlock() = 0;

        virtual bool trylock() = 0;

        virtual void print(std::ostream& os = std::cout) const = 0;

};

class NullThreadLock : 
    public ThreadLock 
{

    public:
        void lock();

        void unlock();

        bool trylock();

        void print(std::ostream& os = std::cout) const;

};

class pThreadLock : 
    public ThreadLock
{

    private:
        pthread_mutex_t mutex_;

        pthread_mutexattr_t attr_;

    public:
        pThreadLock();

        ~pThreadLock();

        void lock();

        void unlock();

        void incref();

        bool trylock();

        pthread_mutex_t* mutex();

        void print(std::ostream& os = std::cout) const;
};

class MultiThreadLock :
    public ThreadLock
{
    private:
        pThreadLock lock_;

        uli lock_count_;

        uli thread_owner_;

    public:
        MultiThreadLock();

        void lock();

        void unlock();

        bool trylock();

        void print(std::ostream& os = std::cout) const;
};


class ThreadEnvironment {

    private:
        static std::list<ThreadWorkspaceAllocator*>* allocators_;

        static uli nthread_;

        uli threadnum_;

    public:
        static void add_allocator(ThreadWorkspaceAllocator* allocator);

        static void allocate();

        static void deallocate();

};


class Thread {

    protected:
        uli threadnum_;

    public:
        Thread(uli threadnum);

        virtual ~Thread();

        virtual void run() = 0;

        uli get_thread_number() const;


};

class ThreadGroup {

    protected:
        Thread** threads_;

        uli nthread_max_;

        uli nthread_;

    public:
        ThreadGroup(uli nthread);

        virtual ~ThreadGroup();

        virtual void add(Thread* thr) = 0;

        virtual void run() = 0;

        virtual void wait() = 0;

        virtual ThreadLock* get_lock() = 0;

        virtual void clear();

        virtual void thread_crash();

};

class pThreadGroup :
    public ThreadGroup
{


    private:
        pthread_t* pthreads_;

        pthread_attr_t* attrs_;

        size_t stack_size_;

        bool running_;

    public:
        pThreadGroup(uli nthread);

        ~pThreadGroup();

        void add(Thread* thr);

        void run();

        void wait();

        void clear();

        void thread_crash();

        size_t get_stack_size();

        ThreadLock* get_lock();

};

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif

