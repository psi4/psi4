#ifndef yeti_thread_h
#define yeti_thread_h
#include <list>

#include <pthread.h>

#include "class.h"
#include "mallocimpl.h"

#include "thread.hpp"

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

};

class NullThreadLock : public ThreadLock {

    public:
        void lock();

        void unlock();

        bool trylock();

};

class pThreadLock : public ThreadLock {

    private:
        pthread_mutex_t mutex_;

        pthread_mutexattr_t attr_;

    public:
        pThreadLock();

        ~pThreadLock();

        void lock();

        void unlock();

        bool trylock();

        pthread_mutex_t* mutex();

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


};

class ThreadGroup {

    protected:
        Thread** threads_;

        uli nthread_max_;

        uli nthread_;

    public:
        ThreadGroup(uli nthread);

        virtual void add(Thread* thr) = 0;

        virtual void run() = 0;

        virtual void wait() = 0;

        virtual ThreadLock* get_lock() = 0;

        virtual void clear();

};

class pThreadGroup : public ThreadGroup {


    private:
        pthread_t* pthreads_;

        pthread_attr_t* attrs_;


    public:
        pThreadGroup(uli nthread);

        void add(Thread* thr);

        void run();

        void wait();

        void clear();

        ThreadLock* get_lock();

};

}

#endif

