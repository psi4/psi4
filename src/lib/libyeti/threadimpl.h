#ifndef thread_impl_h
#define thread_impl_h

#include "thread.h"
#include "runtime.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

#define DECLARE_STATIC_THREAD_WORKSPACE(x) \
static ThreadWorkspaceAllocatorTemplate<x> thr_workspace_##x; \
template<> x** ThreadWorkspaceAllocatorTemplate<x>::workspaces_ = 0


namespace yeti {

template <class T>
class ThreadWorkspaceAllocatorTemplate :
    public ThreadWorkspaceAllocator
{

    private:
        static T** workspaces_;

    public:
        ThreadWorkspaceAllocatorTemplate();

        static T* get_workspace(uli nthread);

        void allocate();

        void deallocate();

};

template <class T>
ThreadWorkspaceAllocatorTemplate<T>::ThreadWorkspaceAllocatorTemplate()
{
    ThreadEnvironment::add_allocator(this);
}

template <class T>
void
ThreadWorkspaceAllocatorTemplate<T>::allocate()
{
    uli nthread = YetiRuntime::nthread();
    workspaces_ = new T*[nthread];
    for (uli i=0; i < nthread; ++i)
    {
        workspaces_[i] = new T;
    }
}

template <class T>
void
ThreadWorkspaceAllocatorTemplate<T>::deallocate()
{
    uli nthread = YetiRuntime::nthread();
    for (uli i=0; i < nthread; ++i)
        delete workspaces_[i];
    delete[] workspaces_;
}

template <class T>
T*
ThreadWorkspaceAllocatorTemplate<T>::get_workspace(uli threadnum)
{
    return workspaces_[threadnum];
}

template <class T>
T*
get_workspace(uli threadnum)
{
    return ThreadWorkspaceAllocatorTemplate<T>::get_workspace(threadnum);
}

}

#ifdef redefine_size_t
#undef size_t
#endif

#endif

